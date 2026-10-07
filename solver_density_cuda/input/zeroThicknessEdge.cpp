#include "input/zeroThicknessEdge.hpp"
#include "input/speciesDB.hpp"   // speciesDB_sha256Hex (SHA-256 の共通実装)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <highfive/H5File.hpp>

namespace {

const char* kTag = "[zeroThicknessEdgeVelocity]";

[[noreturn]] void zteFail(const std::string& why)
{
    std::cerr << kTag << " ERROR: " << why << std::endl;
    std::exit(EXIT_FAILURE);
}

// 属性 (文字列) を読む。無ければ空。h5py の可変長文字列も読める (species_hash と同じ経路)。
std::string attrString(const HighFive::DataSet& ds, const std::string& name)
{
    std::string v;
    if (ds.hasAttribute(name)) {
        try { ds.getAttribute(name).read(v); } catch (const std::exception&) { v = "(読めない型)"; }
    }
    return v;
}

// 属性 (整数) を文字列で返す。無ければ空。
std::string attrIntString(const HighFive::DataSet& ds, const std::string& name)
{
    if (!ds.hasAttribute(name)) return std::string();
    try { long long v = 0; ds.getAttribute(name).read(v); return std::to_string(v); }
    catch (const std::exception&) { return attrString(ds, name); }
}

// ---- 格子署名 (版 zte-mesh-sig-v1)。tools/mark_zero_thickness_edges.py の mesh_signature() と同じ定義 ----
// 項目 (この順):
//   discretization               str  "node"
//   dim                          i8   (1)          vizBface の最大節点数が 2 → 2、3 以上 → 3、どの bcond にも無ければ 0
//   coord                        f8   (nNodes,3)   節点順の座標 (/MESH/COORD を倍精度にしたもの)
//   internal_face_cells          i8   (nNormalPlanes,2) 内部面 ip < nNormalPlanes の (iCells[0], iCells[1])
//   bcond/<physID>/vizBfaceSizes i8   (n)          physID の昇順。元の境界面接続 (無ければ長さ 0)
//   bcond/<physID>/vizBfaceNodes i8   (m)
// 各項目の行 = "<名前> <型> <形 (カンマ区切り)> <little-endian の値のバイト列の SHA-256>\n"。
// 署名 = "zte-mesh-sig-v1\n" と全行を連結した文字列の SHA-256 (16 進 64 桁)。/VALUE と /AUX は含めない
// (restart_field/interp_field で /VALUE が変わっても署名は変わらない)。HDF5 ファイルのバイト列は使わない (圧縮・配置で変わる)。
const char* kSigVersion = "zte-mesh-sig-v1";

template <typename T>
void appendLE(std::string& buf, T v)
{
    // x86 / CUDA ホストは little-endian。そのままのバイト列を積む
    buf.append(reinterpret_cast<const char*>(&v), sizeof(T));
}

std::string sigLine(const std::string& name, const char* dtype, const std::string& shape, const std::string& bytes)
{
    return name + " " + dtype + " " + shape + " " + speciesDB_sha256Hex(bytes) + "\n";
}

std::string meshSignature(const mesh& msh, const std::string& discretization)
{
    std::string text = std::string(kSigVersion) + "\n";
    text += sigLine("discretization", "str", std::to_string(discretization.size()), discretization);

    // 境界は physID の昇順 (readMesh の並びは HDF5 の名前順なので並べ直す)
    std::vector<size_t> order(msh.bconds.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return msh.bconds[a].physID < msh.bconds[b].physID; });
    long long maxFace = 0;
    for (const auto& bc : msh.bconds) for (const auto& f : bc.vizBfaceNodes) maxFace = std::max(maxFace, (long long)f.size());
    const long long dim = (maxFace == 0) ? 0 : (maxFace == 2 ? 2 : 3);
    { std::string b; appendLE<long long>(b, dim); text += sigLine("dim", "i8", "1", b); }

    {
        std::string b;
        b.reserve((size_t)msh.nNodes * 3 * sizeof(double));
        for (geom_int i = 0; i < msh.nNodes; ++i) {
            const auto& c = msh.nodes[i].coords;
            if (c.size() < 3) zteFail("節点 " + std::to_string(i) + " の座標が 3 成分でない");
            for (int d = 0; d < 3; ++d) appendLE<double>(b, (double)c[d]);
        }
        text += sigLine("coord", "f8", std::to_string((long long)msh.nNodes) + ",3", b);
    }
    {
        std::string b;
        b.reserve((size_t)msh.nNormalPlanes * 2 * sizeof(long long));
        for (geom_int ip = 0; ip < msh.nNormalPlanes; ++ip) {
            const auto& ic = msh.planes[ip].iCells;
            if (ic.size() != 2) zteFail("内部面 " + std::to_string(ip) + " の隣接 CV が 2 個でない");
            appendLE<long long>(b, (long long)ic[0]);
            appendLE<long long>(b, (long long)ic[1]);
        }
        text += sigLine("internal_face_cells", "i8", std::to_string((long long)msh.nNormalPlanes) + ",2", b);
    }
    for (const size_t k : order) {
        const auto& bc = msh.bconds[k];
        std::string bs, bn;
        long long m = 0;
        for (const auto& f : bc.vizBfaceNodes) {
            appendLE<long long>(bs, (long long)f.size());
            for (const geom_int n : f) { appendLE<long long>(bn, (long long)n); ++m; }
        }
        const std::string pre = "bcond/" + std::to_string((long long)bc.physID) + "/";
        text += sigLine(pre + "vizBfaceSizes", "i8", std::to_string((long long)bc.vizBfaceNodes.size()), bs);
        text += sigLine(pre + "vizBfaceNodes", "i8", std::to_string(m), bn);
    }
    return speciesDB_sha256Hex(text);
}

} // namespace

std::vector<flow_float> zteLoadVelocityWeight(const solverConfig& cfg, const mesh& msh, ZteLaunchInfo* info)
{
    std::vector<flow_float> w;
    const char* probe = std::getenv("FORGE_EDGE_MASK_PROBE");

    // ---- (1) 無効: 起動エコーだけ (旧経路 = g_zteVelW nullptr でビット同一) ----
    if (cfg.zeroThicknessEdgeVelocity == 0) {
        std::cout << "'zeroThicknessEdgeVelocity' effective: 0 (off; bit-identical to the flux without it)" << std::endl;
        if (probe != nullptr && *probe) {
            std::cout << "[FORGE_EDGE_MASK_PROBE] space.zeroThicknessEdgeVelocity が無効なので w は無い (全節点で従来の再構成)" << std::endl;
        }
        return w;
    }

    // ---- (2) 契約 (初版: plan §4「初版の契約」)。違反は全部まとめて出す ----
    {
        std::vector<std::string> bad;
        if (cfg.discretization != "node") {
            bad.push_back("mesh.discretization = '" + cfg.discretization + "' (node のみ。cell は双対 CV の節点値という前提が無い)");
        }
        if (cfg.solver != "SLAU" && cfg.solver != "SLAU2") {
            bad.push_back("solver = '" + cfg.solver + "' (SLAU / SLAU2 のみ。処置は SLAU_d の内部面だけに入っている)");
        }
        if (cfg.gpu != 1) {
            bad.push_back("gpu = " + std::to_string(cfg.gpu) + " (gpu: 1 のみ。CPU 経路は SLAU_d を通らない)");
        }
        if (cfg.isAxisymmetric != 0) {
            bad.push_back("isAxisymmetric = " + std::to_string(cfg.isAxisymmetric)
                          + " (未検証。2D の端の定義と軸ピンとの交差を検査するまで受け付けない)");
        }
        for (const auto& bc : msh.bconds) {
            if (bc.bcondKind == "periodic") {
                bad.push_back("周期境界 '" + bc.physName + "' (physID " + std::to_string(bc.physID)
                              + ")。周期節点は別 ID で同じ DOF として更新されるので、節点ごとの w が対応点で一致する保証が無い"
                              + " (periodicNode_d.cuh)");
            }
        }
        if (!bad.empty()) {
            std::ostringstream os;
            os << "space.zeroThicknessEdgeVelocity は node + SLAU/SLAU2 + gpu 1 + 非周期 + 非軸対称だけで有効にできる"
                  " (plan convection-zero-thickness-edge-reconstruction §4 初版の契約)。違反:";
            for (const auto& b : bad) os << "\n    - " << b;
            zteFail(os.str());
        }
    }

    // ---- (3) /AUX/<field> の読込と検査 ----
    const std::string& field = cfg.zeroThicknessEdgeVelocityField;
    const std::string path = "/AUX/" + field;
    const geom_int nC = msh.nCells;
    std::vector<float> w32;
    std::string shaAttr, sigAttr, sigVerAttr, gen, genVer, tags, tagIds, rings, nEdge, nMarked, nNodesAttr, created, mode;
    if (msh.nNodes != nC) {
        zteFail("節点数 " + std::to_string((long long)msh.nNodes) + " と CV 数 " + std::to_string((long long)nC)
                + " が違う (node 変換の h5 でない)");
    }
    try {
        HighFive::File file(cfg.meshFileName, HighFive::File::ReadOnly);
        if (!file.exist("AUX") || !file.getGroup("AUX").exist(field)) {
            zteFail("メッシュ h5 '" + cfg.meshFileName + "' に " + path + " が無い。tools/mark_zero_thickness_edges.py で"
                    " この格子 (変換後の h5) に作ること (格子を変えたら作り直す。/AUX は restart_field/interp_field が写さない)");
        }
        HighFive::DataSet ds = file.getDataSet(path);
        const HighFive::DataType t = ds.getDataType();
        if (t.getClass() != HighFive::DataTypeClass::Float || t.getSize() != 4) {
            zteFail(path + " は float32 でなければならない (道具が float32 で書き、SHA-256 も float32 のバイト列で取る)");
        }
        const std::vector<size_t> dims = ds.getDimensions();
        if (dims.size() != 1 || (geom_int)dims[0] != nC) {
            std::ostringstream os;
            os << path << " の形が (" ;
            for (size_t i = 0; i < dims.size(); ++i) os << (i ? "," : "") << dims[i];
            os << ") で、節点数 (node の CV 数) " << nC << " の 1 次元配列でない";
            zteFail(os.str());
        }
        ds.read(w32);
        shaAttr    = attrString(ds, "field_sha256");
        sigAttr    = attrString(ds, "mesh_signature");
        sigVerAttr = attrString(ds, "mesh_signature_version");
        gen        = attrString(ds, "generator");
        genVer     = attrString(ds, "generator_version");
        mode       = attrString(ds, "mode");
        tags       = attrString(ds, "tags");
        tagIds     = attrString(ds, "tag_physids");
        rings      = attrIntString(ds, "rings");
        nEdge      = attrIntString(ds, "n_edge_nodes");
        nMarked    = attrIntString(ds, "n_marked_nodes");
        nNodesAttr = attrIntString(ds, "n_nodes");
        created    = attrString(ds, "created");
    } catch (const std::exception& e) {
        zteFail("メッシュ h5 '" + cfg.meshFileName + "' の " + path + " を読めない: " + e.what());
    }

    // 値域・有限性 (最初の違反の位置と件数を出す)
    {
        geom_int nBad = 0, first = -1;
        for (geom_int i = 0; i < nC; ++i) {
            const float v = w32[i];
            if (!std::isfinite(v) || v < 0.0f || v > 1.0f) { if (first < 0) first = i; ++nBad; }
        }
        if (nBad > 0) {
            std::ostringstream os;
            os << path << " に非有限または [0,1] の外の値が " << nBad << " 個ある (最初は節点 " << first
               << " の " << (double)w32[first] << ")";
            zteFail(os.str());
        }
    }
    // 来歴: 値と格子の SHA-256 が属性と一致すること (手で値を書き換えた・別の格子の /AUX を写した、を止める)
    const std::string fieldSha = speciesDB_sha256Hex(std::string(reinterpret_cast<const char*>(w32.data()), w32.size() * sizeof(float)));
    if (shaAttr.empty()) zteFail(path + " に属性 field_sha256 が無い (道具 mark_zero_thickness_edges.py の出力でない)");
    if (shaAttr != fieldSha) {
        zteFail(path + " の値の SHA-256 " + fieldSha.substr(0, 16) + "… が属性 field_sha256 " + shaAttr.substr(0, 16)
                + "… と一致しない (値が書き換えられている)");
    }
    // 格子署名: ソルバが読み込んだ配列 (座標・内部面の接続・境界面の接続) から再計算して属性と照合する。
    // 節点番号の付け替え・別格子の /AUX の流用は、ここで時間更新の前に止まる
    if (sigAttr.empty() || sigVerAttr.empty()) {
        zteFail(path + " に属性 mesh_signature / mesh_signature_version が無い (どの格子の重みか確かめられない)");
    }
    if (sigVerAttr != kSigVersion) {
        zteFail(path + " の mesh_signature_version '" + sigVerAttr + "' はこのソルバの '" + kSigVersion + "' と違う");
    }
    const std::string meshSig = meshSignature(msh, cfg.discretization);
    if (sigAttr != meshSig) {
        zteFail(path + " の属性 mesh_signature " + sigAttr.substr(0, 16) + "… が、読み込んだ格子から再計算した署名 "
                + meshSig.substr(0, 16) + "… と一致しない (別の格子・節点番号の違う格子の重み)。この h5 で道具を回し直すこと");
    }

    // ---- (4) 統計と起動ログ ----
    w.resize(nC);
    geom_int nLt1 = 0, nZero = 0;
    double wmin = 1.0, wmax = 0.0, vol = 0.0;
    double lo[3] = { std::numeric_limits<double>::infinity(),  std::numeric_limits<double>::infinity(),  std::numeric_limits<double>::infinity() };
    double hi[3] = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
    for (geom_int i = 0; i < nC; ++i) {
        const float v = w32[i];
        w[i] = (flow_float)v;
        wmin = std::min(wmin, (double)v); wmax = std::max(wmax, (double)v);
        if (v < 1.0f) {
            ++nLt1;
            if (v == 0.0f) ++nZero;
            vol += (double)msh.cells[i].volume;
            for (int d = 0; d < 3; ++d) {
                const double x = (double)msh.cells[i].centCoords[d];
                lo[d] = std::min(lo[d], x); hi[d] = std::max(hi[d], x);
            }
        }
    }
    // 対象の内部面 (片側でも w<1 の節点を持つ面) と面の側の数。node の内部面は ip < nNormalPlanes
    geom_int nFaces = 0, nSides = 0;
    for (geom_int ip = 0; ip < msh.nNormalPlanes && ip < (geom_int)msh.planes.size(); ++ip) {
        const auto& ic = msh.planes[ip].iCells;
        bool hit = false;
        for (const geom_int c : ic) {
            if (c >= 0 && c < nC && w32[c] < 1.0f) { hit = true; ++nSides; }
        }
        if (hit) ++nFaces;
    }
    if (nLt1 == 0) {
        std::cout << kTag << " 注意: w が全節点で 1 (処置は恒等。対照用の重みか確認すること; mode=" << (mode.empty() ? "?" : mode) << ")" << std::endl;
    }
    std::cout << "'zeroThicknessEdgeVelocity' effective: 1 (field " << path << ", w<1 nodes " << nLt1
              << ", field_sha256 " << fieldSha.substr(0, 16) << ", mesh_signature " << meshSig.substr(0, 16) << ")" << std::endl;
    std::cout << kTag << " 生成元: " << (gen.empty() ? "?" : gen) << " " << (genVer.empty() ? "" : genVer)
              << " mode=" << (mode.empty() ? "?" : mode) << " tags=" << (tags.empty() ? "?" : tags)
              << " tag_physids=" << (tagIds.empty() ? "?" : tagIds)
              << " rings=" << (rings.empty() ? "?" : rings) << " n_edge_nodes(E)=" << (nEdge.empty() ? "?" : nEdge)
              << " n_marked_nodes(S)=" << (nMarked.empty() ? "?" : nMarked) << " n_nodes=" << (nNodesAttr.empty() ? "?" : nNodesAttr)
              << (created.empty() ? "" : " created=" + created) << std::endl;
    std::cout << kTag << " 検証: 値 [0,1]・有限、field_sha256 一致 (" << fieldSha << ")、mesh_signature 一致 ("
              << meshSig << ", " << kSigVersion << ")" << std::endl;
    std::cout << kTag << " w<1 の節点 " << nLt1 << " (うち w=0 が " << nZero << ") / " << nC
              << "、w の min/max = " << wmin << " / " << wmax
              << "、対象の内部面 " << nFaces << " (面の側 " << nSides << ")"
              << "、対象体積の合計 " << vol << std::endl;
    if (nLt1 > 0) {
        std::cout << kTag << " w<1 の節点の座標範囲: x [" << lo[0] << ", " << hi[0] << "] y [" << lo[1] << ", " << hi[1]
                  << "] z [" << lo[2] << ", " << hi[2] << "]" << std::endl;
    }
    if (!nMarked.empty() && nMarked != std::to_string((long long)nLt1)) {
        std::cout << kTag << " 注意: 属性 n_marked_nodes=" << nMarked << " と w<1 の節点数 " << nLt1
                  << " が違う (w を 0/1 以外で書いた道具か、別の数え方)" << std::endl;
    }
    if (info != nullptr) {
        info->enabled = 1;
        info->field = field;
        info->meshSignature = meshSig;
        info->meshSignatureVersion = kSigVersion;
        info->fieldSha256 = fieldSha;
        info->nZero = (long long)nZero; info->nLt1 = (long long)nLt1; info->nNodes = (long long)nC;
    }
    if (cfg.convMethod == 0) {
        std::cout << kTag << " 注意: convMethod 0 (1 次) では面の速度が既に節点値なので、この処置は作用しない" << std::endl;
    }

    // ---- 診断: 指定節点の w (例 FORGE_EDGE_MASK_PROBE="517160,517199") ----
    if (probe != nullptr && *probe) {
        std::string s(probe);
        size_t pos = 0;
        while (pos < s.size()) {
            size_t q = s.find(',', pos);
            if (q == std::string::npos) q = s.size();
            const std::string tok = s.substr(pos, q - pos);
            pos = q + 1;
            if (tok.empty()) continue;
            char* end = nullptr;
            const long long id = std::strtoll(tok.c_str(), &end, 10);
            if (end == tok.c_str() || *end != '\0') {
                std::cout << "[FORGE_EDGE_MASK_PROBE] '" << tok << "' は節点 ID でない" << std::endl;
                continue;
            }
            if (id < 0 || id >= (long long)nC) {
                std::cout << "[FORGE_EDGE_MASK_PROBE] 節点 " << id << ": 範囲外 (0.." << (nC - 1) << ")" << std::endl;
                continue;
            }
            const float v = w32[id];
            std::cout << "[FORGE_EDGE_MASK_PROBE] 節点 " << id << ": w=" << (double)v
                      << (v < 1.0f ? " (処置の対象)" : " (対象外 = 従来の再構成)")
                      << " 座標 (" << msh.cells[id].centCoords[0] << ", " << msh.cells[id].centCoords[1] << ", "
                      << msh.cells[id].centCoords[2] << ")" << std::endl;
        }
    }
    return w;
}
