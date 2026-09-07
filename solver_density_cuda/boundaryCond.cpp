#include "boundaryCond.hpp"

#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/boundaryCond_d.cuh"
#include "cuda_forge/axisymmetricSource_d.cuh"
#include "cuda_forge/ransBoundary_d.cuh"
#include "cuda_forge/ransWallFunction_d.cuh"
#include "cuda_forge/fluct_variables_d.cuh"

#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <array>
#include <cuda_runtime.h>

using namespace std;

bcondConfFormat::bcondConfFormat(){};

void readBcondConfig(solverConfig& cfg , vector<bcond>& bconds)
{
    map<int,bcondConfFormat> bcondConfMap;
    const auto isInletKind = [](const std::string& kind) {
        return kind.rfind("inlet_", 0) == 0;
    };

    std::string bcondConfigFileName  = "bcondConfig.yaml";

    try {
        YAML::Node config = YAML::LoadFile(bcondConfigFileName);

        for(auto bcInYaml : config) {
            bcondConfFormat bcf;

            std::string bname = bcInYaml.first.as<std::string>();
            cout << "bname = " << bname << endl;

            // akirameta
            //if (bcInYaml.find("physID") == bcInYaml.end()) {
            //    // not found
            //    cerr << "Error: couldn't find physID of " << bname << " in bcondConfig.yaml\n";
            //    exit();
            //}

            int physID = bcInYaml.second["physID"].as<int>();

            std::map<std::string, int> inputInts_temp;
            for (auto bcInYaml_int : bcInYaml.second["ints"])
            {
                std::string key = bcInYaml_int.first.as<std::string>();
                int val = bcInYaml_int.second.as<int>();
                inputInts_temp[key] = val;
            }

            std::map<std::string, flow_float> inputFloats_temp;
            for (auto bcInYaml_float : bcInYaml.second["floats"])
            {
                std::string key = bcInYaml_float.first.as<std::string>();
                flow_float val = bcInYaml_float.second.as<flow_float>();
                inputFloats_temp[key] = val;
            }

            std::string kind = bcInYaml.second["kind"].as<std::string>();

            int outputHDFflg_temp = bcInYaml.second["outputHDFflg"].as<int>();

            bcf.physID = physID;
            bcf.physName = bname;
            bcf.kind = kind;
            bcf.inputInts = inputInts_temp;
            bcf.inputFloats = inputFloats_temp;
            bcf.outputHDFflg = outputHDFflg_temp;

            if (cfg.LESorRANS == 2 && isInletKind(kind)) {
                const bool hasK = inputFloats_temp.find("k") != inputFloats_temp.end();
                const bool hasOmega = inputFloats_temp.find("omega") != inputFloats_temp.end();
                if (!hasK || !hasOmega) {
                    cerr << "Error: inlet boundary '" << bname << "' of kind '" << kind
                         << "' must define floats 'k' and 'omega' when LESorRANS=2." << endl;
                    exit(EXIT_FAILURE);
                }
            }

            bcondConfMap[bcf.physID] = bcf;

            cout << "physID=" <<  bcf.physID << endl;
            cout << "kind=" << bcf.kind<< endl;
            cout << "outputHDF=" << bcf.outputHDFflg<< endl;
        }

    } catch(const YAML::BadFile& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);

    } catch(const YAML::ParserException& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);
    }

    // set bconds from cofig file (YAML file)
    for (bcond& bc : bconds)
    {
        cout << "bc.physID=" << bc.physID << endl;

        const auto confIt = bcondConfMap.find(bc.physID);
        if (confIt == bcondConfMap.end()) {
            // not found
            cerr << "Error: couldn't find physID " << bc.physID << " in bcondConfig.yaml\n";
            exit(EXIT_FAILURE);
        }
        const bcondConfFormat& bcf = confIt->second;

        cout << "bcf.kind=" << bcf.kind << endl;

        bc.bcondKind    = bcf.kind;
        bc.physName     = bcf.physName;
        bc.inputInts    = bcf.inputInts;
        bc.inputFloats  = bcf.inputFloats;
        bc.outputHDFflg = bcf.outputHDFflg;

        // copy value types to bcond
        const auto valueTypesIt = bcf.valueTypesOfBC.find(bcf.kind);
        if (valueTypesIt == bcf.valueTypesOfBC.end()) {
            cerr << "Error: unsupported boundary condition kind " << bcf.kind << " for physID " << bc.physID << "\n";
            exit(EXIT_FAILURE);
        }

        const map<string,int>& valueTypes = valueTypesIt->second;
        for (auto& vt : valueTypes)
        {
            bc.valueTypes[vt.first] = vt.second;
        }

        // 多成分 TP (M5): 超音速入口は全状態固定なので入口組成 Y_s も config から
        // 与える。inlet_* 種別に対して Y0..Y{n-1} を type-1 (uniform float read) として
        // 動的登録する。既存の単成分入口 config (Y 未指定) を壊さないよう、未指定なら
        // Y0=1, それ以外=0 を既定値とする (= 第 1 化学種のみの単成分入口)。
        if (cfg.nSpecies >= 2 && isInletKind(bcf.kind)) {
            for (int s = 0; s < cfg.nSpecies; s++) {
                const std::string yname = "Y" + std::to_string(s);
                bc.valueTypes[yname] = 1;             // uniform float read
                bc.bplaneValNames.push_back(yname);   // bcondInitVariables が確保・初期化する対象に
                if (bc.inputFloats.find(yname) == bc.inputFloats.end()) {
                    bc.inputFloats[yname] = (s == 0) ? 1.0 : 0.0;  // 既定: 単成分 (sp[0])
                }
            }
        }

        bc.bcondInitVariables(cfg.gpu); // allocate and set boundary variables
    }
};

// ---------------------------------------------------------------------------
// 入口分布プロファイル: inlet bcond の per-face 境界値 (bvar) を CSV テーブルから
// face 重心座標で補間してセットする。inlet_* の kernel は per-face の bvar (Uxb 等) を
// そのまま読むため、kernel 無改修で「一様でない入口分布 (壁法則 BL 等)」を実現できる。
//
// 有効化: bcondConfig.yaml の対象 inlet に `ints: {inletProfile: 1}` を付ける。
// テーブル: run dir の `inlet_profile_<physID>.csv`。1 行目ヘッダで補間方向と量を指定:
//   - 先頭の連続する x/y/z 列 = 補間座標。1 列 (例 `y`) なら 1D 線形補間、3 列 (`x y z`) なら 3D 最近傍。
//   - 残り列 = bvar 量名 (Ux Uy Uz Tt Pt Ts Ps ro k omega ...; その inlet が持つ bvar のみ反映)。
//   例 (1D-y):  `y Ux Uy Uz`\n`1.0 0 0 0`\n`1.1 12 0 0` ...
//   例 (3D):    `x y z Ux Uy Uz` ...
// methods/boundary.md, plans/active/boundary-inlet-profile.md 参照。
// ---------------------------------------------------------------------------
static std::vector<std::string> splitWS(const std::string& s)
{
    std::vector<std::string> out; std::istringstream iss(s); std::string t;
    while (iss >> t) out.push_back(t);
    return out;
}

void applyInletProfiles(solverConfig& cfg , mesh& msh)
{
    for (bcond& bc : msh.bconds)
    {
        const auto it = bc.inputInts.find("inletProfile");
        if (it == bc.inputInts.end() || it->second != 1) continue;

        const std::string fname = "inlet_profile_" + std::to_string(bc.physID) + ".csv";
        std::ifstream fin(fname);
        if (!fin) {
            std::cerr << "[applyInletProfiles] inletProfile=1 but file '" << fname
                      << "' not found (physID=" << bc.physID << ").\n";
            exit(EXIT_FAILURE);
        }

        // ---- ヘッダ: 先頭の連続 x/y/z = 補間座標、残り = 量名 ----
        std::string line;
        while (std::getline(fin, line)) { // 先頭の空行/コメントを飛ばす
            std::string t = line; size_t p = t.find_first_not_of(" \t\r\n");
            if (p == std::string::npos || t[p] == '#') continue;
            break;
        }
        std::vector<std::string> hdr = splitWS(line);
        std::vector<int> axisIdx;     // 補間軸 (0=x,1=y,2=z) の列
        size_t ncoord = 0;
        for (const auto& h : hdr) {
            if (h == "x") { axisIdx.push_back(0); ncoord++; }
            else if (h == "y") { axisIdx.push_back(1); ncoord++; }
            else if (h == "z") { axisIdx.push_back(2); ncoord++; }
            else break;
        }
        if (ncoord == 0) {
            std::cerr << "[applyInletProfiles] " << fname << ": header must start with x/y/z coordinate column(s).\n";
            exit(EXIT_FAILURE);
        }
        std::vector<std::string> qnames(hdr.begin() + ncoord, hdr.end());

        // ---- データ ----
        std::vector<std::array<double,3>> rowC;     // 座標 (使う軸のみ有効)
        std::vector<std::vector<double>>  rowQ;     // 量
        while (std::getline(fin, line)) {
            std::vector<std::string> tk = splitWS(line);
            if (tk.size() < hdr.size()) continue;
            std::array<double,3> c{0,0,0};
            for (size_t k = 0; k < ncoord; ++k) c[axisIdx[k]] = std::stod(tk[k]);
            std::vector<double> q(qnames.size());
            for (size_t k = 0; k < qnames.size(); ++k) q[k] = std::stod(tk[ncoord + k]);
            rowC.push_back(c); rowQ.push_back(q);
        }
        const int nrow = (int)rowC.size();
        if (nrow < 1) { std::cerr << "[applyInletProfiles] " << fname << ": no data rows.\n"; exit(EXIT_FAILURE); }

        // 1D の場合は補間軸で昇順ソート (線形補間用)
        const bool oneD = (ncoord == 1);
        std::vector<int> order(nrow);
        for (int i = 0; i < nrow; ++i) order[i] = i;
        if (oneD) {
            const int ax = axisIdx[0];
            std::sort(order.begin(), order.end(), [&](int a, int b){ return rowC[a][ax] < rowC[b][ax]; });
        }

        auto isInputQuantity = [&](const std::string& qn) {
            const auto vt = bc.valueTypes.find(qn);
            return vt != bc.valueTypes.end() && vt->second == 1 && bc.bvar.find(qn) != bc.bvar.end();
        };

        // ---- 各 inlet face で補間し bvar をセット ----
        for (size_t i = 0; i < bc.iPlanes.size(); ++i)
        {
            const geom_int ip = bc.iPlanes[i];
            const double fc[3] = { (double)msh.planes[ip].centCoords[0],
                                   (double)msh.planes[ip].centCoords[1],
                                   (double)msh.planes[ip].centCoords[2] };
            std::vector<double> qv(qnames.size());
            if (oneD) {
                const int ax = axisIdx[0];
                const double cf = fc[ax];
                // 端はクランプ、内部は線形補間 (order は昇順)
                if (cf <= rowC[order[0]][ax]) { qv = rowQ[order[0]]; }
                else if (cf >= rowC[order[nrow-1]][ax]) { qv = rowQ[order[nrow-1]]; }
                else {
                    int lo = 0;
                    while (lo < nrow-1 && rowC[order[lo+1]][ax] < cf) ++lo;
                    const int a = order[lo], b = order[lo+1];
                    const double ca = rowC[a][ax], cb = rowC[b][ax];
                    const double w = (cb > ca) ? (cf - ca)/(cb - ca) : 0.0;
                    for (size_t k = 0; k < qnames.size(); ++k) qv[k] = (1.0-w)*rowQ[a][k] + w*rowQ[b][k];
                }
            } else {
                // 3D 最近傍
                int best = 0; double bd = 1e300;
                for (int r = 0; r < nrow; ++r) {
                    double d = 0;
                    for (size_t k = 0; k < ncoord; ++k) { const int ax = axisIdx[k]; const double dd = fc[ax]-rowC[r][ax]; d += dd*dd; }
                    if (d < bd) { bd = d; best = r; }
                }
                qv = rowQ[best];
            }
            // bvar へ反映: この種別の「入力量」(valueTypes==1, config の一様値を読む量) だけ。
            // type 0 の量 (例 inlet_Pressure の Ux/Ps, inlet_uniformVelocity の Tt) は kernel が毎 step 書き換える
            // 出力/作業用で、CSV を書いても効かないので反映しない (ログの IGNORED に出す)。
            for (size_t k = 0; k < qnames.size(); ++k) {
                if (!isInputQuantity(qnames[k])) continue;
                auto bit = bc.bvar.find(qnames[k]);
                if (bit != bc.bvar.end() && i < bit->second.size()) bit->second[i] = (flow_float)qv[k];
            }
        }

        // ---- device へ再アップロード ----
        if (cfg.gpu == 1) {
            for (const auto& qn : qnames) {
                if (!isInputQuantity(qn)) continue;
                auto bit = bc.bvar.find(qn);
                auto dit = bc.bvar_d.find(qn);
                if (bit != bc.bvar.end() && dit != bc.bvar_d.end() && dit->second != nullptr) {
                    cudaMemcpy(dit->second, bit->second.data(),
                               bc.iPlanes.size()*sizeof(flow_float), cudaMemcpyHostToDevice);
                }
            }
        }
        // 適用できた列 (この種別の bvar に存在) と無視した列を分けて報告する。
        // 例: inlet_uniformVelocity に Tt 列を書いても kernel は読まない (bvar 無し) → ignored に出す。
        std::string applied, ignored;
        for (const auto& qn : qnames) {
            if (isInputQuantity(qn)) applied += " " + qn; else ignored += " " + qn;
        }
        std::cout << "[applyInletProfiles] physID=" << bc.physID << " kind=" << bc.bcondKind
                  << ": set " << qnames.size() << " quantities from " << fname
                  << " (" << (oneD ? "1D interp" : (std::to_string(ncoord) + "D nearest")) << ", " << nrow << " rows, "
                  << bc.iPlanes.size() << " faces). applied:" << (applied.empty() ? " (none)" : applied)
                  << (ignored.empty() ? "" : "  IGNORED (not an input quantity of this kind):" + ignored) << "\n";
        if (applied.empty()) {
            std::cerr << "[applyInletProfiles] " << fname << ": no column matches a boundary value of kind "
                      << bc.bcondKind << " (see procedures/inlet-profile.md).\n";
            exit(EXIT_FAILURE);
        }
        // 多成分入口: 組成列があれば各 face で 0<=Y_s<=1・ΣY_s≈1 を検査する (入口カーネルは負値クリップ+正規化するが、
        // node ピンは値をそのまま使うので、ここで弾く)。
        if (cfg.nSpecies >= 2) {
            bool anyY = false;
            for (const auto& qn : qnames) if (qn.size() >= 2 && qn[0] == 'Y' && isInputQuantity(qn)) anyY = true;
            if (anyY) {
                for (size_t i = 0; i < bc.iPlanes.size(); ++i) {
                    double ysum = 0.0;
                    for (int sidx = 0; sidx < cfg.nSpecies; ++sidx) {
                        const auto it2 = bc.bvar.find("Y" + std::to_string(sidx));
                        const double y = (it2 != bc.bvar.end()) ? (double)it2->second[i] : 0.0;
                        if (!(y >= -1.0e-6 && y <= 1.0 + 1.0e-6)) {
                            std::cerr << "[applyInletProfiles] " << fname << ": Y" << sidx << "=" << y
                                      << " at face " << i << " is outside [0,1].\n";
                            exit(EXIT_FAILURE);
                        }
                        ysum += y;
                    }
                    if (std::fabs(ysum - 1.0) > 1.0e-3) {
                        std::cerr << "[applyInletProfiles] " << fname << ": sum of Y_s = " << ysum << " at face " << i
                                  << " (must be 1 within 1e-3; write all species columns or use gen_inlet_profile.py).\n";
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }
}

void applyBconds(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p , fluct_variables& fluct)
{
    if (cfg.gpu == 0) {
        cerr << "Error: gpu=0 is not supported in applyBconds" << endl;
        exit(EXIT_FAILURE);
    } else if (cfg.gpu ==1) { // gpu
        // node-centered でも壁 ghost は残す (mirror で u_face=0 を与え、勾配/粘性の境界閉性を保つ)。
        // 壁ノードの u=0 厳密化は別途 Dirichlet (wall_flag + enforceWallNoSlip/zeroWallMomentumResidual) で行う。
        for (auto& bc : msh.bconds)
        {
            //if      (bc.bcondKind == "wall_isothermal") { wall_isothermal_d(cfg , bc , msh , var , mat_p); }
            if      (bc.bcondKind == "slip") slip_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p);
            else if (bc.bcondKind == "axis") slip_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p);
            else if (bc.bcondKind == "wall") wall_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p);
            else if (bc.bcondKind == "wall_isothermal") wall_isothermal_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p);
            else if (bc.bcondKind == "inlet_uniformVelocity") { inlet_uniformVelocity_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_fluctVelocity") { inlet_fluctVelocity_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p , fluct); }
            else if (bc.bcondKind == "outlet_statPress") { outlet_statPress_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_Pressure") { inlet_Pressure_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_Pressure_dir") { inlet_Pressure_dir_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "outflow") { outflow_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "periodic") { periodic_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
    } else {
        cerr << "Error: unsupported gpu flag " << cfg.gpu << endl;
        exit(EXIT_FAILURE);
    }
}

void applyRansScalarBoundaries(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!(cfg.gpu == 1 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) {
        return;
    }

    // automatic wall treatment: wall-function 生産 wf_pk を bc ループ前に全セル -1 (inactive) 化。
    // 各 wall bc の computeWallFrictionSST が自分の wall-adjacent セルだけ >=0 に埋める。
    if (cfg.wallTreatmentSST == 1) {
        initWallFunctionPk_d_wrapper(cfg , cuda_cfg , msh , var);
    }

    for (auto& bc : msh.bconds)
    {
        ransBoundary_d_wrapper(cfg , cuda_cfg , bc , msh , var);
    }

    // node k Dirichlet の状態ピン (全 wall bc の roK_wf が揃った後)。explicit RK3 でも実効化する
    // (従来は implicit 更新カーネルのみで、explicit では k が初期値凍結するバグ。ransWallFunction_d.cu)。
    applyNodeKwfStatePin_d_wrapper(cfg , cuda_cfg , msh , var);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

//void copyBcondsGradient(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p)
//{
//    if (cfg.gpu == 0) { // cpu
//        cerr << "Error: can't use with cpu" << endl;
//        exit(EXIT_FAILURE);
//
//    } else if (cfg.gpu ==1) { // gpu
//        for (auto& bc : msh.bconds)
//        {
//            copyBcondsGradient_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
//        }
//        gpuErrchk( cudaPeekAtLastError() );
//        gpuErrchkKernelSync();
//
//    }
//}
//