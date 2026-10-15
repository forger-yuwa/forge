// 本段の監視用の出力 (協調停止・節点の状態・境界の質量流束)。説明は monitorOutputs.hpp と
// procedures/solver-settings.md「output」、設計は plans/active/tooling-nozzle-card-local-mesh.md §4.4 1′。
#include "output/monitorOutputs.hpp"

#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <highfive/H5File.hpp>

#include "cuda_forge/cudaWrapper.cuh"

using HighFive::File;

namespace monitorOutputs {

namespace {

const char* const kStopRequest = "STOP_REQUEST";
const char* const kStoppedAt = "STOPPED_AT";
const char* const kSnapshotRequest = "SNAPSHOT_REQUEST";
const char* const kSnapshotDone = "SNAPSHOT_DONE";
const char* const kSnapshotFailed = "SNAPSHOT_FAILED";
const char* const kMassfluxCsv = "boundary_massflux.csv";
// 節点の出力の変数 (全場の出力と同じ名前・同じ型 flow_float の値をそのまま写す)
const char* const kNodeProbeNames[] = {"ro", "Ux", "Uy", "Uz", "P", "T", "sonic"};

struct StopState {
    bool requested = false;
    int step = -1;
    std::string request;   // STOP_REQUEST の中身 (末尾の空白・改行は落とす)
    std::vector<std::string> writeErrors;   // 停止の step の出力の書き込みで出た例外 (空なら成功)
};
StopState& stopState() { static StopState s; return s; }

// 随時の書き出し (SNAPSHOT_REQUEST)。1 回の要求ごとに作り直す。
struct SnapshotState {
    bool active = false;      // この step で要求を処理中
    bool disabled = false;    // 要求を消せなかった (毎 step 同じ要求を繰り返さないよう以後は見ない)
    int step = -1;
    int handledStep = -1;     // 最後に締めた (SNAPSHOT_DONE / SNAPSHOT_FAILED を書いた) step
    std::string request;
    std::vector<std::string> writeErrors;
};
SnapshotState& snapState() { static SnapshotState s; return s; }

struct NodeProbeState {
    bool enabled = false;
    std::string file;
    std::vector<geom_int> ids;
};
NodeProbeState& nodeProbe() { static NodeProbeState s; return s; }

struct MassfluxBoundary {
    int physID;
    std::string name;
    std::vector<geom_int> planes;   // bc.iPlanes (node: 境界の半割の双対面の global plane index)
};
struct MassfluxState {
    bool enabled = false;
    double factor = 1.0;   // 軸対称 (axisymMethod 0) は 2π、それ以外は 1
    std::vector<MassfluxBoundary> bcs;
    std::ofstream ofs;
    int lastStep = -1;   // 最後に書けた行の step (協調停止の締めで突き合わせる)
};
MassfluxState& massflux() { static MassfluxState s; return s; }

// 協調停止・随時の書き出しの step に、4 種の出力のほかに書いたファイル (noteStepFiles; ソルバ内 CHT の固体の状態など)
struct ExtraFiles { int step = -1; std::vector<std::string> files; };
ExtraFiles& extraFiles() { static ExtraFiles s; return s; }

bool fileExists(const char* path)
{
    struct stat st;
    return ::stat(path, &st) == 0;
}

// 通常のファイルの大きさ [byte]。無い・通常のファイルでないときは -1。
long long regularFileSize(const std::string& path)
{
    struct stat st;
    if (::stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return -1;
    return (long long)st.st_size;
}

// テキストファイルの最後の空でない行 (# で始まる行も含む)。読めなければ空。
std::string lastNonEmptyLine(const std::string& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) return std::string();
    ifs.seekg(0, std::ios::end);
    const std::streamoff size = ifs.tellg();
    const std::streamoff n = std::min<std::streamoff>(size, 1 << 16);
    std::string tail((size_t)n, '\0');
    ifs.seekg(size - n);
    ifs.read(&tail[0], n);
    tail.resize((size_t)ifs.gcount());
    size_t end = tail.size();
    while (end > 0 && (tail[end - 1] == '\n' || tail[end - 1] == '\r')) --end;
    const size_t beg = tail.rfind('\n', end == 0 ? 0 : end - 1);
    return tail.substr(beg == std::string::npos ? 0 : beg + 1, end - (beg == std::string::npos ? 0 : beg + 1));
}

// 要求ファイルの中身 (先頭 4 KiB、末尾の空白・改行は落とす)
std::string readRequestText(const char* path)
{
    std::ifstream ifs(path, std::ios::binary);
    std::string txt;
    if (ifs) {
        char buf[4096];
        ifs.read(buf, sizeof(buf));
        txt.assign(buf, (size_t)ifs.gcount());
    }
    size_t n = txt.size();
    while (n > 0 && (txt[n - 1] == ' ' || txt[n - 1] == '\t' || txt[n - 1] == '\r' || txt[n - 1] == '\n')) --n;
    return txt.substr(0, n);
}

// CSV の行の先頭の欄 (step) を整数で。数でなければ false。
bool firstFieldInt(const std::string& line, long long& v)
{
    const size_t comma = line.find(',');
    const std::string f = line.substr(0, comma);
    if (f.empty()) return false;
    size_t pos = 0;
    try { v = std::stoll(f, &pos, 10); } catch (const std::exception&) { return false; }
    return pos == f.size();
}

// outStepInterval ごと・最後の step・force。境界の h5 (outputBconds_H5_XDMF) と同じ規則。
bool isOutputStep(const solverConfig& cfg, int iStep, bool force)
{
    if (force) return true;
    const bool lastStep = (cfg.nStepOuter > 0) && (iStep == cfg.nStepOuter);
    if (lastStep) return true;
    return (cfg.outStepInterval > 0) && (iStep % cfg.outStepInterval == 0) && (iStep >= cfg.outStepStart);
}

std::string trimRight(const std::string& s)
{
    size_t n = s.size();
    while (n > 0 && (s[n - 1] == ' ' || s[n - 1] == '\t' || s[n - 1] == '\r' || s[n - 1] == '\n')) --n;
    return s.substr(0, n);
}

std::string trimBoth(const std::string& s)
{
    size_t b = 0;
    while (b < s.size() && (s[b] == ' ' || s[b] == '\t' || s[b] == '\r')) ++b;
    return trimRight(s.substr(b));
}

std::string jsonEscape(const std::string& s)
{
    std::string o;
    o.reserve(s.size() + 8);
    for (unsigned char c : s) {
        switch (c) {
            case '"':  o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n";  break;
            case '\r': o += "\\r";  break;
            case '\t': o += "\\t";  break;
            default:
                if (c < 0x20) { char buf[8]; std::snprintf(buf, sizeof(buf), "\\u%04x", c); o += buf; }
                else o += static_cast<char>(c);   // UTF-8 はそのまま
        }
    }
    return o;
}

// 節点の番号の表を読む。1 行 1 つ、# 以降と空行は無視。数でない行・範囲外・重複は例外。
std::vector<geom_int> readNodeList(const std::string& path, geom_int nValid)
{
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("output.nodeProbeFile: '" + path + "' を開けない");
    std::vector<geom_int> ids;
    std::set<geom_int> seen;
    std::string line;
    int lineNo = 0;
    while (std::getline(ifs, line)) {
        ++lineNo;
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        const std::string t = trimBoth(line);
        if (t.empty()) continue;
        long long v = 0;
        size_t pos = 0;
        try {
            v = std::stoll(t, &pos, 10);
        } catch (const std::exception&) {
            pos = 0;
        }
        if (pos == 0 || pos != t.size())
            throw std::runtime_error("output.nodeProbeFile: '" + path + "' の " + std::to_string(lineNo)
                                     + " 行目 '" + t + "' は節点の番号 (整数 1 つ) でない");
        if (v < 0 || v >= (long long)nValid)
            throw std::runtime_error("output.nodeProbeFile: '" + path + "' の " + std::to_string(lineNo)
                                     + " 行目の節点 " + std::to_string(v) + " は範囲外 (0 〜 "
                                     + std::to_string((long long)nValid - 1) + ")");
        if (!seen.insert((geom_int)v).second)
            throw std::runtime_error("output.nodeProbeFile: '" + path + "' の " + std::to_string(lineNo)
                                     + " 行目の節点 " + std::to_string(v) + " が重複している");
        ids.push_back((geom_int)v);
    }
    if (ids.empty()) throw std::runtime_error("output.nodeProbeFile: '" + path + "' に節点が 1 つも無い");
    return ids;
}

// device (gpu 1) または host (gpu 0) の CV 値 [nCells] を buf に写す。var.c (host の出力用の写し) は書き換えない。
void fetchCellValues(const solverConfig& cfg, const mesh& msh, variables& var, const std::string& name,
                     std::vector<flow_float>& buf)
{
    buf.resize(msh.nCells);
    if (cfg.gpu == 1) {
        cudaWrapper::cudaMemcpy_D2H_wrapper(var.c_d.at(name), buf.data(), msh.nCells);
    } else {
        const std::vector<flow_float>& h = var.c.at(name);
        std::copy(h.begin(), h.begin() + msh.nCells, buf.begin());
    }
}

// step k の出力を確かめる (協調停止・随時の書き出しで共通)。流量の CSV は最後に書けた行と、ディスク上の最後の行の step = k。
// 書くべきファイル (全場・outputHDFflg 1 の境界の h5/xmf・節点の h5・CSV) がディスクにあって空でないこと。
// 足りないものは errs に足す。files・sizes は STOPPED_AT / SNAPSHOT_DONE に載せる一覧。
void checkStepOutputs(const mesh& msh, int k, std::vector<std::string>& errs,
                      std::vector<std::string>& files, std::vector<long long>& sizes)
{
    if (massflux().enabled) {
        long long last = -1;
        if (massflux().lastStep != k || !firstFieldInt(lastNonEmptyLine(kMassfluxCsv), last) || last != (long long)k)
            errs.push_back(std::string(kMassfluxCsv) + ": step " + std::to_string(k) + " の行が無い (最後の行の step "
                           + std::to_string(last) + ")");
    }
    const std::string ks = std::to_string(k);
    files.push_back("res_" + ks + ".h5");
    files.push_back("res_" + ks + ".xmf");
    for (const bcond& bc : msh.bconds) {
        if (bc.outputHDFflg != 1) continue;
        const std::string base = "res_" + bc.physName + "_" + std::to_string(bc.physID) + "_" + ks;
        files.push_back(base + ".h5");
        files.push_back(base + ".xmf");
    }
    if (nodeProbe().enabled) files.push_back("res_nodes_" + ks + ".h5");
    if (massflux().enabled) files.push_back(kMassfluxCsv);
    if (extraFiles().step == k) for (const auto& f : extraFiles().files) files.push_back(f);
}

void fileSizes(const std::vector<std::string>& files, std::vector<std::string>& errs, std::vector<long long>& sizes)
{
    sizes.clear();
    for (const auto& f : files) {
        const long long sz = regularFileSize(f);
        if (sz <= 0) errs.push_back(f + ": " + (sz < 0 ? "無い" : "空"));
        sizes.push_back(sz);
    }
}

std::string filesJson(const std::vector<std::string>& files, const std::vector<long long>& sizes, const char* sep)
{
    std::ostringstream js;
    for (size_t i = 0; i < files.size(); ++i)
        js << "{\"path\": \"" << jsonEscape(files[i]) << "\", \"bytes\": " << sizes[i] << "}"
           << (i + 1 < files.size() ? sep : "");
    return js.str();
}

// JSON lines のファイルに 1 行を書き足す。書けたら true。
bool appendLine(const char* path, const std::string& line)
{
    std::ofstream ofs(path, std::ios::out | std::ios::app);
    ofs << line << "\n";
    ofs.flush();
    return static_cast<bool>(ofs);
}

}  // namespace

void init(const solverConfig& cfg, const mesh& msh)
{
    // 前の run の残り (同じ run dir での再起動など) で 1 step 目に止まる・見張りが完了と取り違えるのを防ぐ。
    if (fileExists(kStopRequest))
        throw std::runtime_error(std::string("作業ディレクトリに ") + kStopRequest
                                 + " が既にある (前の run の残り?)。消してから起動する");
    if (fileExists(kStoppedAt))
        throw std::runtime_error(std::string("作業ディレクトリに ") + kStoppedAt
                                 + " が既にある (前の run の残り?)。消してから起動する");
    if (fileExists(kSnapshotRequest))
        throw std::runtime_error(std::string("作業ディレクトリに ") + kSnapshotRequest
                                 + " が既にある (前の run の残り?)。消してから起動する");

    if (!cfg.nodeProbeFile.empty()) {
        if (cfg.discretization != "node")
            throw std::runtime_error("output.nodeProbeFile は node の離散化 (mesh.discretization: node) だけで使える");
        // node では CV の番号 = primal の節点の番号 (全場の出力も VALUE を節点に Center='Node' で書く)
        const geom_int nValid = std::min(msh.nNodes, msh.nCells);
        NodeProbeState& np = nodeProbe();
        np.ids = readNodeList(cfg.nodeProbeFile, nValid);
        np.file = cfg.nodeProbeFile;
        np.enabled = true;
        std::cout << "[nodeProbe] " << cfg.nodeProbeFile << ": " << np.ids.size()
                  << " 節点 -> res_nodes_<step>.h5 (outStepInterval " << cfg.outStepInterval
                  << " ごと + 最後の step + 協調停止の step)" << std::endl;
    }

    if (cfg.boundaryMassflux == 1) {
        if (cfg.discretization != "node")
            throw std::runtime_error("output.boundaryMassflux は node の離散化だけで使える (cell の境界の massflux は未検証)");
        if (cfg.gpu != 1)
            throw std::runtime_error("output.boundaryMassflux は gpu: 1 だけで使える (CPU の経路の境界の massflux は未検証)");
        if (cfg.isAxisymmetric == 1 && cfg.axisymMethod != 0)
            throw std::runtime_error("output.boundaryMassflux は軸対称の axisymMethod 0 (面ベクトルの r 重み) だけに対応 "
                                     "(axisymMethod 1 は面が r で重みづけされず、実寸の流量を面の和から作れない)");
        MassfluxState& mf = massflux();
        // 軸対称 (axisymMethod 0) の面ベクトルは r で重みづけした 1 ラジアンあたりの量 (variables.cpp setStructuralVariables_d;
        // 2π は出力の側で掛ける規約、methods/axisymmetric/theory.md「2π の取扱い」)。
        mf.factor = (cfg.isAxisymmetric == 1) ? 6.283185307179586476925286766559 : 1.0;   // 2π
        for (const bcond& bc : msh.bconds) {
            if (bc.bcondKind == "periodic") continue;   // 周期面の流束は領域の内部の面
            mf.bcs.push_back({(int)bc.physID, bc.physName, bc.iPlanes});
        }
        mf.ofs.open(kMassfluxCsv, std::ios::out | std::ios::trunc);
        if (!mf.ofs) throw std::runtime_error(std::string("output.boundaryMassflux: ") + kMassfluxCsv + " を開けない");
        std::ostream& o = mf.ofs;
        o << "# " << kMassfluxCsv << " — forge output.boundaryMassflux: 1 (plans/active/tooling-nozzle-card-local-mesh.md §4.4 1'(d))\n";
        o << "# raw_<name>_<physID>: その境界の全面 (node: 境界の半割の双対面) の質量流束 massflux = |S_f|*rho_b*(U_b.n_f) の和\n";
        o << "#   (convectiveFlux_boundary_d が境界の流束に使った値そのもの; rho_b・U_b は境界の状態 bvar)。\n";
        o << "#   その step の最後の流束評価の値 (定常の陰解法では、その step の更新に使った step の頭の状態の流束)。\n";
        o << "# 符号: 面の向きは計算領域から外向き -> 流出が正、流入が負。\n";
        if (cfg.isAxisymmetric == 1) {
            o << "# 軸対称 (mesh.isAxisymmetric 1, axisymMethod 0): S_f は r で重みづけした 1 ラジアンあたりの面 -> raw は [kg/s/rad]、"
                 "mdot_<name>_<physID> = 2*pi*raw が実寸の流量 [kg/s]。\n";
        } else {
            o << "# 非軸対称: S_f は格子の面そのもの -> mdot_<name>_<physID> = raw [kg/s] (2 次元の格子なら奥行き 1 m あたり)。\n";
        }
        o << "# step: 全場 res_<step>.h5・境界の h5 と同じ番号 (residual_history.csv の step 列は step-1)。\n";
        o << "step";
        for (const auto& b : mf.bcs) o << ",mdot_" << b.name << "_" << b.physID << ",raw_" << b.name << "_" << b.physID;
        o << "\n";
        o.flush();
        if (!o) throw std::runtime_error(std::string("output.boundaryMassflux: ") + kMassfluxCsv + " の見出しを書けない");
        mf.enabled = true;
        std::cout << "[boundaryMassflux] " << kMassfluxCsv << ": " << mf.bcs.size() << " 境界 (factor "
                  << (cfg.isAxisymmetric == 1 ? "2*pi" : "1") << ", outStepInterval " << cfg.outStepInterval
                  << " ごと + 最後の step + 協調停止の step)" << std::endl;
    }
}

bool pollStopRequest(int iStep)
{
    StopState& s = stopState();
    if (s.requested) return true;
    if (!fileExists(kStopRequest)) return false;   // 毎 step の負担はこの stat 1 回だけ
    s.requested = true;
    s.step = iStep;
    s.request = readRequestText(kStopRequest);
    // 理由は改行を含みうるので、log を 1 行に保つため JSON と同じ規則でエスケープして出す
    std::cout << "[coopStop] step " << iStep << ": STOP_REQUEST を検出 (中身 \"" << jsonEscape(s.request)
              << "\") -> この step の全場・境界・節点・流量を書いて正常終了する (STOPPED_AT)" << std::endl;
    return true;
}

bool stopRequested() { return stopState().requested; }
int stoppedStep() { return stopState().step; }

void noteStopWriteError(const std::string& what)
{
    stopState().writeErrors.push_back(what);
    std::cerr << "[coopStop] ERROR: 停止の step の出力の書き込みに失敗: " << what << std::endl;
}

int finishStop(const solverConfig& cfg, const mesh& msh, bool residualStreamOk)
{
    const StopState& s = stopState();
    if (!s.requested) return 0;
    const int k = s.step;
    std::vector<std::string> errs = s.writeErrors;

    // 残差履歴: ストリームが健全で、ディスク上の最後の行の step 列 = k − 1 (全場の番号 k の step を記録した行)
    {
        const std::string rh = "residual_history.csv";
        long long last = -1;
        if (!residualStreamOk) errs.push_back(rh + ": ストリームが書き込みの失敗の状態");
        else if (!firstFieldInt(lastNonEmptyLine(rh), last) || last != (long long)k - 1)
            errs.push_back(rh + ": 最後の行の step 列が " + std::to_string(last) + " (期待 " + std::to_string(k - 1) + ")");
    }
    // 流量の CSV・全場・境界・節点のファイル (随時の書き出しと共通) と残差履歴
    std::vector<std::string> files;
    std::vector<long long> sizes;
    checkStepOutputs(msh, k, errs, files, sizes);
    files.push_back("residual_history.csv");
    fileSizes(files, errs, sizes);

    if (!errs.empty()) {
        for (const auto& e : errs) std::cerr << "[coopStop] ERROR: " << e << std::endl;
        std::cerr << "[coopStop] ERROR: step " << k << " の出力を確かめられないので " << kStoppedAt
                  << " を置かずに終了コード 3 で終える" << std::endl;
        return 3;
    }

    std::ostringstream js;
    js << "{\n"
       << "  \"status\": \"ok\",\n"
       << "  \"step\": " << k << ",\n"
       << "  \"residual_step\": " << (k - 1) << ",\n"
       << "  \"request\": \"" << jsonEscape(s.request) << "\",\n"
       << "  \"nStepOuter\": " << cfg.nStepOuter << ",\n"
       << "  \"files\": [\n    " << filesJson(files, sizes, ",\n    ") << "\n  ]\n}\n";
    const std::string tmp = std::string(kStoppedAt) + ".tmp";
    {
        std::ofstream ofs(tmp, std::ios::out | std::ios::trunc);
        ofs << js.str();
        ofs.flush();
        if (!ofs) {
            std::cerr << "[coopStop] ERROR: " << tmp << " を書けない (完了の通知に失敗) -> 終了コード 3" << std::endl;
            return 3;
        }
    }
    if (std::rename(tmp.c_str(), kStoppedAt) != 0) {
        std::cerr << "[coopStop] ERROR: " << tmp << " -> " << kStoppedAt << " の rename に失敗 (完了の通知に失敗) -> 終了コード 3" << std::endl;
        return 3;
    }
    return 0;
}

bool pollSnapshotRequest(int iStep)
{
    SnapshotState& s = snapState();
    s.active = false;
    if (s.disabled || !fileExists(kSnapshotRequest)) return false;   // 毎 step の負担はこの stat 1 回だけ
    s.active = true;
    s.step = iStep;
    s.request = readRequestText(kSnapshotRequest);
    s.writeErrors.clear();
    return true;
}

void noteSnapshotWriteError(const std::string& what)
{
    snapState().writeErrors.push_back(what);
}

void finishSnapshot(const mesh& msh)
{
    SnapshotState& s = snapState();
    if (!s.active) return;
    s.active = false;
    const int k = s.step;
    s.handledStep = k;
    std::vector<std::string> errs = s.writeErrors;
    std::vector<std::string> files;
    std::vector<long long> sizes;
    checkStepOutputs(msh, k, errs, files, sizes);
    fileSizes(files, errs, sizes);

    const std::string head = "{\"step\": " + std::to_string(k) + ", \"pid\": " + std::to_string((long long)::getpid())
                           + ", \"request\": \"" + jsonEscape(s.request) + "\", ";
    // 完了 (または失敗) の記録を先に書き、それから要求を消す: 見張りは「要求が消えた」を見たら必ず記録の行を読める。
    bool recorded = false;
    if (errs.empty()) {
        recorded = appendLine(kSnapshotDone, head + "\"files\": [" + filesJson(files, sizes, ", ") + "]}");
        if (recorded)
            std::cout << "[snapshot] step " << k << ": 全場・境界・節点・流量を書いた (中身 \"" << jsonEscape(s.request)
                      << "\", " << kSnapshotDone << " に記録、計算は続ける)" << std::endl;
        else
            errs.push_back(std::string(kSnapshotDone) + " に書き足せない");
    }
    if (!errs.empty()) {
        std::string ej;
        for (size_t i = 0; i < errs.size(); ++i) ej += (i ? ", \"" : "\"") + jsonEscape(errs[i]) + "\"";
        recorded = appendLine(kSnapshotFailed, head + "\"errors\": [" + ej + "]}");
        std::cerr << "[snapshot] ERROR: step " << k << " の書き出しに失敗 (計算は続ける): " << errs.front()
                  << (errs.size() > 1 ? " ほか" : "") << (recorded ? "" : " (SNAPSHOT_FAILED にも書けない)") << std::endl;
    }
    // 要求を消す (失敗でも消す: 同じ失敗を毎 step 繰り返さない)。消せなければ以後は要求を見ない。
    if (std::remove(kSnapshotRequest) != 0) {
        s.disabled = true;
        std::cerr << "[snapshot] ERROR: " << kSnapshotRequest << " を消せない -> この run では以後の要求を見ない" << std::endl;
        appendLine(kSnapshotFailed, head + "\"errors\": [\"" + std::string(kSnapshotRequest)
                                    + " を消せない; この run では以後の要求を見ない\"]}");
    }
}

bool snapshotHandledAt(int k) { return snapState().handledStep == k; }

void noteStepFiles(int k, const std::vector<std::string>& files)
{
    ExtraFiles& x = extraFiles();
    if (x.step != k) { x.step = k; x.files.clear(); }
    for (const auto& f : files)
        if (std::find(x.files.begin(), x.files.end(), f) == x.files.end()) x.files.push_back(f);
}

void writeNodeProbes(const solverConfig& cfg, const mesh& msh, variables& var, int iStep, bool force)
{
    const NodeProbeState& np = nodeProbe();
    if (!np.enabled || !isOutputStep(cfg, iStep, force)) return;

    const std::string fname = "res_nodes_" + std::to_string(iStep) + ".h5";
    File file(fname, File::ReadWrite | File::Truncate);
    file.createAttribute<int>("step", HighFive::DataSpace::From(iStep)).write(iStep);
    file.createAttribute<std::string>("nodeProbeFile", HighFive::DataSpace::From(np.file)).write(np.file);

    // 座標は全場の /MESH/COORD と同じ並び (x0,y0,z0,x1,...) の 1 次元
    std::vector<geom_float> coord;
    coord.reserve(3 * np.ids.size());
    for (const geom_int id : np.ids) {
        coord.push_back(msh.nodes[id].coords[0]);
        coord.push_back(msh.nodes[id].coords[1]);
        coord.push_back(msh.nodes[id].coords[2]);
    }
    file.createDataSet("/MESH/COORD", coord);
    file.createDataSet("/MESH/NODE_ID", np.ids);

    std::vector<flow_float> buf, vals(np.ids.size());
    for (const char* name : kNodeProbeNames) {
        fetchCellValues(cfg, msh, var, name, buf);
        for (size_t k = 0; k < np.ids.size(); ++k) vals[k] = buf[np.ids[k]];
        file.createDataSet(std::string("/VALUE/") + name, vals);
    }
}

void writeBoundaryMassflux(const solverConfig& cfg, const mesh& msh, variables& var, int iStep, bool force)
{
    MassfluxState& mf = massflux();
    if (!mf.enabled || !isOutputStep(cfg, iStep, force)) return;

    std::vector<flow_float> flux(msh.nPlanes);
    cudaWrapper::cudaMemcpy_D2H_wrapper(var.p_d.at("massflux"), flux.data(), msh.nPlanes);

    char buf[64];
    std::ostream& o = mf.ofs;
    o << iStep;
    for (const auto& b : mf.bcs) {
        double raw = 0.0;   // 面の順に double で足す (決定的)
        for (const geom_int ip : b.planes) raw += (double)flux[ip];
        std::snprintf(buf, sizeof(buf), ",%.17g", mf.factor * raw);
        o << buf;
        std::snprintf(buf, sizeof(buf), ",%.17g", raw);
        o << buf;
    }
    o << "\n";
    o.flush();
    // 書き込みの失敗を黙らせない (h5 の書き込みの失敗は HighFive が例外を投げるのと揃える)
    if (!o) throw std::runtime_error(std::string(kMassfluxCsv) + " への step " + std::to_string(iStep) + " の行の書き込みに失敗");
    mf.lastStep = iStep;
}

// 診断 FORGE_DIAG_OUTPUT_INVARIANCE 専用 (§5.1 #1c)。止めている間は (c)・(d) を書かず、停止・随時の書き出しの確認の一覧にも載せない
// (enabled を下ろすだけ。節点の表・CSV のストリームは持ったまま)。
namespace {
struct DiagSuspend { bool active = false; bool np = false; bool mf = false; };
DiagSuspend& diagSuspend() { static DiagSuspend s; return s; }
}  // namespace

void diagSuspendMonitorKeys(bool suspend)
{
    DiagSuspend& d = diagSuspend();
    if (suspend && !d.active) {
        d.np = nodeProbe().enabled;
        d.mf = massflux().enabled;
        nodeProbe().enabled = false;
        massflux().enabled = false;
        d.active = true;
    } else if (!suspend && d.active) {
        nodeProbe().enabled = d.np;
        massflux().enabled = d.mf;
        d.active = false;
    }
}

bool monitorKeysEnabled()
{
    const DiagSuspend& d = diagSuspend();
    return d.active ? (d.np || d.mf) : (nodeProbe().enabled || massflux().enabled);
}

}  // namespace monitorOutputs
