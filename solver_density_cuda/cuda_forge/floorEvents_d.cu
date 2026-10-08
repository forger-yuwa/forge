// 毎更新の EOS 床事象のカウンタ (出力専用)。定義は floorEvents_d.cuh 冒頭と
// plans/active/tooling-sern-te-wake-grid.md §4「床の判定」(2026-10-08 のやり直し)。
//
// 記録: run ディレクトリの floor_events.csv (1 ファイル、追記)。列は kHeader で固定。行の種類 (kind):
//   session_begin  forge の起動ごとに 1 行 (session は同じファイル内の通し番号。段階起動の各段・restart は別 session)。
//                  note に計測の条件 (経路・下限・節点数・cfg_fnv …) を書く。step の起点はこの起動の 0。
//   init           初期化の EOS (読み込んだ場が対象。定常でもここでの床は基準状態へ入る。更新事象とは別枠)
//   eos            時間ステップ内の EOS。step = その EOS を呼んだ外側 step (1 起点)、inner = その step 内の EOS の通し番号。
//                  定常陰解法 (path steady_implicit) では 1 step に 1 回で、q_index = step − 1 (= 直前の更新の結果を判定)。
//                  nStepInner の sweep は保存量を更新しないので更新として数えない (EOS も呼ばれない)。
//   aux            時間ステップの外で呼ばれた EOS (二相の残差監査など)。更新事象としては数えない
//   audit          終了時の監査。最後の更新の結果 Q_N をコピー上で同じ前処理・同じ判定にかける (q_index = N)
//   session_end    正常終了の印。note に監査の状態 (audit=done / unsupported:<理由>)
// 判定区間 (a, b] は、入口 Q_a = 行 q_index a (別枠)、Q_{a+1}..Q_{b−1} = 行 q_index a+1..b−1、Q_b = 監査行
// (b = N のとき) で覆う。行が無い・欠けた step は 0 件ではなく判定不能 (tools/check_floor_events.py)。
// 件数は 64 bit。整合しない値 (件数 > 節点数、ID の記録数 ≠ 件数) は overflow = 1。ID 列は種類ごとに
// FLOOR_EVENT_IDS 個まで (超えたら ids_truncated に種類名)。

#include "floorEvents_d.cuh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "cudaWrapper.cuh"
#include "thermo_d.cuh"
#include "speciesTransport_d.cuh"     // species_roY_device_ptr()
#include "nodeWallDirichlet_d.cuh"    // enforceWallNoSlip_d_wrapper / applyNodeIsothermalWallPin
#include "axisymmetricSource_d.cuh"   // enforceAxisSymmetry_d_wrapper
#include "dependentVariables_d.cuh"   // dependentVariables_d_wrapper / dependentVariablesTminTP

namespace {

const char* kFileName = "floor_events.csv";
const char* kHeader =
    "kind,session,step,inner,q_index,path,n_real,n_ghost,n_nonfinite_real,n_nonfinite_ghost,"
    "nT_real,nRho_real,nP_real,nT_ghost,nRho_ghost,nP_ghost,nT_uneval_real,nT_mismatch_real,n_near_real,"
    "dRhoE_T_sum,dRhoE_T_max,dRho_sum,dRho_max,dP_sum,dP_max,"
    "ids_T,ids_Rho,ids_P,ids_near,ids_truncated,overflow,note";
const double kNearBand = 1.0;   // 床近傍の幅 [K] (plan §4: T ≤ T_min + 1 K)

struct State {
    bool on = false;
    FILE* fp = nullptr;
    int session = 0;
    std::string path;               // steady_implicit / dual_time / explicit
    FloorEventPhase phase = FloorEventPhase::Init;
    int step = 0;
    int inner = 0;                  // phase 内の EOS の通し番号
    long long nReal = 0, nGhost = 0;
    FloorEventAcc* acc_d = nullptr;
    double* hT_d = nullptr;
    double tFloor = 0.0;
    std::string auditStatus = "not_run";
    // 試験用の注入 (FORGE_FLOOR_TEST_INJECT=step:node:mode:value; mode = de | rhoscale)
    int injStep = -1; long long injNode = -1; int injMode = 0; double injValue = 0.0;
};
State g;

std::string fmtG(double v)
{
    char b[64]; std::snprintf(b, sizeof(b), "%.17g", v); return b;
}

double bitsToDouble(unsigned long long u)
{
    double d; std::memcpy(&d, &u, sizeof(d)); return d;
}

void writeLine(const std::string& s)
{
    if (g.fp == nullptr) return;
    std::fputs(s.c_str(), g.fp); std::fputc('\n', g.fp);
    std::fflush(g.fp);   // 異常終了でもそこまでの行が残るように (行の欠けは判定不能として扱われる)
}

// 化学種ごとの h_s(T_min)/MW_s (thermo_cph_mix の項と同じ演算)
__global__ void floorEventsHTmin_d(const SpeciesThermo* sp, int n, double T, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    for (int i = 0; i < n; ++i) {
        double cpi, hi;
        thermo_cph_molar(sp[i], T, &cpi, &hi);
        out[i] = hi / sp[i].MW;
    }
}

// 試験用: 1 節点の保存量だけを書き換え、診断値を返す。mode 1: roe を ρ (e_mix(T_min) + ek + de) に置き換える
// (密度・運動量・組成は不変)。mode 2: その節点の保存量 (ρ, ρU, ρE, ρk, ρω, ρY_s) を同じ倍率 value で縮める
// (速度・比内部エネルギー・組成は不変で、密度と圧力だけが下がる)。
__global__ void floorEventsInjectScale_d(geom_int node, double value, int nSpecies, flow_float** roY,
                                         flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz,
                                         flow_float* roe, flow_float* roK, flow_float* roOmega, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    out[0] = ro[node]; out[1] = roe[node];
    const flow_float f = (flow_float)value;
    ro[node] *= f; roUx[node] *= f; roUy[node] *= f; roUz[node] *= f; roe[node] *= f;
    if (roK != nullptr) roK[node] *= f;
    if (roOmega != nullptr) roOmega[node] *= f;
    if (roY != nullptr) for (int s = 0; s < nSpecies; ++s) roY[s][node] *= f;
    out[2] = ro[node]; out[3] = roe[node]; out[4] = 0.0; out[5] = 0.0;
}

__global__ void floorEventsInject_d(geom_int node, int thermalMethod, flow_float gamma, flow_float cp, flow_float roMin,
                                    const SpeciesThermo* sp, int nSpecies, flow_float** roY, const double* hT, double T_min,
                                    double de, flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz,
                                    flow_float* roe, flow_float* roeN, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    const flow_float ro_temp = max(ro[node], roMin);
    const flow_float ux = roUx[node]/ro_temp, uy = roUy[node]/ro_temp, uz = roUz[node]/ro_temp;
    const flow_float ek = 0.5f*(ux*ux + uy*uy + uz*uz);
    double e_floor, cv;
    if (thermalMethod == 2) {
        double Y[THERMO_MAX_SPECIES];
        if (nSpecies <= 1 || roY == nullptr) { Y[0] = 1.0; }
        else {
            double ys = 0.0;
            for (int s = 0; s < nSpecies; ++s) { double y = (double)roY[s][node]/(double)ro_temp; if (y < 0.0) y = 0.0; Y[s] = y; ys += y; }
            const double inv = 1.0/(ys > 1.0e-30 ? ys : 1.0e-30);
            for (int s = 0; s < nSpecies; ++s) Y[s] *= inv;
        }
        const int n = (nSpecies <= 1 || roY == nullptr) ? 1 : nSpecies;
        const double R = thermo_R_mix(sp, n, Y);
        double hs = 0.0; for (int s = 0; s < n; ++s) hs += Y[s]*hT[s];
        e_floor = hs - R*T_min;
        double cpm, hm; thermo_cph_mix(sp, n, Y, T_min, &cpm, &hm);
        cv = cpm - R;
    } else {
        cv = (double)cp/(double)gamma;
        e_floor = cv*T_min;
    }
    const flow_float roe_old = roe[node];
    const flow_float roe_new = (flow_float)((double)ro_temp*(e_floor + (double)ek + de));
    roe[node] = roe_new;
    if (roeN != nullptr) roeN[node] = roe_new;
    out[0] = roe_old; out[1] = roe_new; out[2] = e_floor; out[3] = cv; out[4] = ro_temp; out[5] = ek;
}

std::string idsField(const FloorEventAcc& a, int kind, bool* truncated)
{
    const unsigned int n = a.nIds[kind];
    const unsigned int m = std::min<unsigned int>(n, (unsigned int)FLOOR_EVENT_IDS);
    std::vector<long long> v(a.ids[kind], a.ids[kind] + m);
    std::sort(v.begin(), v.end());
    std::string s;
    for (size_t i = 0; i < v.size(); ++i) { if (i) s += ";"; s += std::to_string(v[i]); }
    *truncated = (n > (unsigned int)FLOOR_EVENT_IDS);
    return s;
}

long long qIndexOf()
{
    if (g.phase == FloorEventPhase::Audit) return g.step;
    if (g.phase == FloorEventPhase::Step && g.path == "steady_implicit" && g.inner == 0) return (long long)g.step - 1;
    return -1;
}

const char* kindOf()
{
    switch (g.phase) {
        case FloorEventPhase::Init:  return "init";
        case FloorEventPhase::Step:  return "eos";
        case FloorEventPhase::Aux:   return "aux";
        case FloorEventPhase::Audit: return "audit";
    }
    return "aux";
}

}  // namespace

bool floorEventsEnabled() { return g.on; }

void floorEventsOpen(const solverConfig& cfg, const mesh& msh, const std::string& cfgFnv)
{
    if (cfg.floorEvents == 0) return;
    if (cfg.gpu != 1) {
        std::fprintf(stderr, "[floor-events] 拒否: output.floorEvents は GPU 経路 (gpu=1) のみ対応 (CPU の EOS には計数が無い)\n");
        std::exit(1);
    }
    g.path = (cfg.isImplicit == 1) ? (cfg.unsteady == 1 ? "dual_time" : "steady_implicit") : "explicit";
    g.nReal = (long long)msh.nCells; g.nGhost = (long long)msh.nCells_all - (long long)msh.nCells;
    g.tFloor = (cfg.thermalMethod == 2) ? dependentVariablesTminTP() : (double)cfg.tMin;

    // 既存の記録には追記する (段階起動の各段・restart の痕跡を消さない)。書式が違うファイルには混ぜない。
    int prior = 0;
    bool hasHeader = false;
    {
        std::ifstream in(kFileName);
        std::string line;
        while (in && std::getline(in, line)) {
            if (!hasHeader) {
                if (line != kHeader) {
                    std::fprintf(stderr, "[floor-events] 拒否: 既存の %s の見出しが現行の書式と違う (別の版の記録に追記しない)\n", kFileName);
                    std::exit(1);
                }
                hasHeader = true;
                continue;
            }
            if (line.rfind("session_begin,", 0) == 0) ++prior;
        }
    }
    g.fp = std::fopen(kFileName, "a");
    if (g.fp == nullptr) {
        std::fprintf(stderr, "[floor-events] 拒否: %s を開けない (記録の無い run を作らない)\n", kFileName);
        std::exit(1);
    }
    if (!hasHeader) writeLine(kHeader);   // 新規 (または空) のファイル
    else {
        // 前の起動が行の途中で落ちていたら、その行を閉じてから書く (欠けた行は判定側で読めない行 = 判定不能になる)
        std::ifstream in(kFileName, std::ios::binary | std::ios::ate);
        if (in && in.tellg() > 0) {
            in.seekg(-1, std::ios::end); char c = 0; in.get(c);
            if (c != '\n') { std::fputc('\n', g.fp); std::fflush(g.fp); }
        }
    }
    g.session = prior + 1;

    gpuErrchk( cudaMalloc(&g.acc_d, sizeof(FloorEventAcc)) );
    gpuErrchk( cudaMemset(g.acc_d, 0, sizeof(FloorEventAcc)) );
    if (cfg.thermalMethod == 2) {
        const int n = std::max(cfg.nSpecies, 1);
        gpuErrchk( cudaMalloc(&g.hT_d, sizeof(double)*THERMO_MAX_SPECIES) );
        gpuErrchk( cudaMemset(g.hT_d, 0, sizeof(double)*THERMO_MAX_SPECIES) );
        floorEventsHTmin_d<<<1, 1>>>(thermo_species_device_ptr(), n, g.tFloor, g.hT_d);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
    }
    if (const char* e = std::getenv("FORGE_FLOOR_TEST_INJECT"); e != nullptr && *e != '\0') {
        char mode[16] = {0};
        const bool ok = std::sscanf(e, "%d:%lld:%15[a-z]:%lf", &g.injStep, &g.injNode, mode, &g.injValue) == 4
                        && g.injStep >= 1 && g.injNode >= 0 && g.injNode < g.nReal
                        && (std::strcmp(mode, "de") == 0 || (std::strcmp(mode, "rhoscale") == 0 && g.injValue > 0.0));
        if (!ok) {
            std::fprintf(stderr, "[floor-events] 拒否: FORGE_FLOOR_TEST_INJECT は step:node:de:<J/kg> か step:node:rhoscale:<倍率>"
                                 " (step>=1, 0<=node<nCells): '%s'\n", e);
            std::exit(1);
        }
        g.injMode = (std::strcmp(mode, "de") == 0) ? 1 : 2;
    }
    g.on = true;

    bool periodic = (msh.nPeriodicMembers > 0);
    for (const auto& bc : msh.bconds) if (bc.bcondKind == "periodic") periodic = true;
    std::ostringstream note;
    note << "format=1;path=" << g.path << ";thermo=" << (cfg.thermalMethod == 2 ? "tp" : "cpg")
         << ";condensation=" << cfg.condensation << ";condEquilibrium=" << cfg.condEquilibrium
         << ";tfloor=" << fmtG(g.tFloor) << ";romin=" << fmtG((double)cfg.roMin) << ";pmin=" << fmtG((double)cfg.pMin)
         << ";near_band=" << fmtG(kNearBand) << ";ids_max=" << FLOOR_EVENT_IDS
         << ";discretization=" << cfg.discretization << ";periodic=" << (periodic ? 1 : 0) << ";axisymmetric=" << cfg.isAxisymmetric
         << ";nStepOuter=" << cfg.mainLoopCount() << ";prior_sessions=" << prior << ";cfg_fnv=" << cfgFnv
         << ";inject=" << (g.injStep > 0 ? (g.injMode == 1 ? "de" : "rhoscale") : "0");
    std::ostringstream row;
    row << "session_begin," << g.session << ",0,-1,-1," << g.path << "," << g.nReal << "," << g.nGhost
        << ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,,,,0," << note.str();
    writeLine(row.str());
    std::printf("[floor-events] 有効: %s に追記 (session %d, path %s, T_floor %.6g K, 実節点 %lld)\n",
                kFileName, g.session, g.path.c_str(), g.tFloor, g.nReal);
}

void floorEventsSetPhase(FloorEventPhase phase, int step)
{
    g.phase = phase; g.step = step; g.inner = 0;
}

FloorEventDev floorEventsKernelArgs()
{
    FloorEventDev d{nullptr, nullptr, 0.0, 0.0};
    if (!g.on) return d;
    gpuErrchk( cudaMemset(g.acc_d, 0, sizeof(FloorEventAcc)) );
    d.acc = g.acc_d; d.hTminMass = g.hT_d; d.tFloorTP = g.tFloor; d.nearBand = kNearBand;
    return d;
}

void floorEventsAfterEos()
{
    if (!g.on) return;
    FloorEventAcc a;
    gpuErrchk( cudaMemcpy(&a, g.acc_d, sizeof(FloorEventAcc), cudaMemcpyDeviceToHost) );
    bool tr[FE_ID_KINDS];
    const std::string idT = idsField(a, FE_ID_T, &tr[0]), idR = idsField(a, FE_ID_RHO, &tr[1]),
                      idP = idsField(a, FE_ID_P, &tr[2]), idN = idsField(a, FE_ID_NEAR, &tr[3]);
    std::string trunc;
    const char* nm[FE_ID_KINDS] = {"T", "Rho", "P", "near"};
    for (int k = 0; k < FE_ID_KINDS; ++k) if (tr[k]) { if (!trunc.empty()) trunc += ";"; trunc += nm[k]; }
    // 整合の検査: 実節点の件数 ≤ 実節点数、ghost の件数 ≤ ghost 数、ID の記録数 = 件数
    const unsigned long long R = (unsigned long long)g.nReal, G = (unsigned long long)g.nGhost;
    const bool overflow =
        a.nNonfinite[0] > R || a.nT[0] > R || a.nRho[0] > R || a.nP[0] > R || a.nTUneval > R || a.nTMismatch > R || a.nNear > R ||
        a.nNonfinite[1] > G || a.nT[1] > G || a.nRho[1] > G || a.nP[1] > G ||
        (unsigned long long)a.nIds[FE_ID_T] != a.nT[0] || (unsigned long long)a.nIds[FE_ID_RHO] != a.nRho[0] ||
        (unsigned long long)a.nIds[FE_ID_P] != a.nP[0] || (unsigned long long)a.nIds[FE_ID_NEAR] != a.nNear;
    std::ostringstream row;
    row << kindOf() << "," << g.session << "," << g.step << "," << g.inner << "," << qIndexOf() << "," << g.path << ","
        << g.nReal << "," << g.nGhost << "," << a.nNonfinite[0] << "," << a.nNonfinite[1] << ","
        << a.nT[0] << "," << a.nRho[0] << "," << a.nP[0] << "," << a.nT[1] << "," << a.nRho[1] << "," << a.nP[1] << ","
        << a.nTUneval << "," << a.nTMismatch << "," << a.nNear << ","
        << fmtG(a.sum_dRhoE_T) << "," << fmtG(bitsToDouble(a.max_dRhoE_T)) << ","
        << fmtG(a.sum_dRho) << "," << fmtG(bitsToDouble(a.max_dRho)) << ","
        << fmtG(a.sum_dP) << "," << fmtG(bitsToDouble(a.max_dP)) << ","
        << idT << "," << idR << "," << idP << "," << idN << "," << trunc << "," << (overflow ? 1 : 0) << ",";
    writeLine(row.str());
    ++g.inner;
}

void floorEventsAudit(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int lastStep)
{
    if (!g.on) return;
    // EOS 拘束形平衡 (condEquilibrium 2) の EOS は rog[0] (通常の配列) を書き換えるので、コピー上で閉じられない
    if (cfg.condensation == 1 && cfg.condEquilibrium == 2) {
        g.auditStatus = "unsupported:condEquilibrium2";
        std::printf("[floor-events] 終了時の監査は未対応 (condEquilibrium 2 の EOS が通常の配列 rog を書く) — 記録は判定不能になる\n");
        return;
    }
    const bool verify = (std::getenv("FORGE_FLOOR_AUDIT_VERIFY") != nullptr);
    const size_t nca = (size_t)msh.nCells_all;
    // 試験用: 監査の前後で全 cell 配列 (ghost 込み) をバイト比較する
    std::map<std::string, std::vector<flow_float>> before;
    if (verify) {
        for (auto& kv : var.c_d) if (kv.second) {
            auto& v = before[kv.first]; v.resize(nca);
            gpuErrchk( cudaMemcpy(v.data(), kv.second, nca*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        }
    }
    // EOS と境界ピンが読み書きする配列だけをコピーに差し替える (他の配列はこの区間で読まれも書かれもしない)
    static const char* kSwap[] = {"ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega",
                                  "P", "Ht", "sonic", "k", "omega", "T", "Ux", "Uy", "Uz", "gamma", "cp", "Rmix",
                                  "roN", "roUyN", "roeN"};   // roN/roUyN/roeN: 軸ピンが commit の基準も射影する
    std::vector<std::pair<std::string, flow_float*>> saved;
    auto restore = [&]() {
        for (auto& s : saved) {
            flow_float* scratch = var.c_d[s.first];
            var.c_d[s.first] = s.second;
            cudaFree(scratch);
        }
        saved.clear();
    };
    flow_float* injT = nullptr;
    try {
        for (const char* name : kSwap) {
            auto it = var.c_d.find(name);
            if (it == var.c_d.end() || it->second == nullptr) continue;
            flow_float* scratch = nullptr;
            gpuErrchk( cudaMalloc(&scratch, nca*sizeof(flow_float)) );
            gpuErrchk( cudaMemcpy(scratch, it->second, nca*sizeof(flow_float), cudaMemcpyDeviceToDevice) );
            saved.emplace_back(name, it->second);
            it->second = scratch;
        }
        // 時間ステップ冒頭 (assembleResidualPre) と同じ前処理: no-slip → 軸 → 等温壁、の後に同じ EOS (同じ判定)
        enforceWallNoSlip_d_wrapper(cfg , cuda_cfg , msh , var);
        enforceAxisSymmetry_d_wrapper(cfg , cuda_cfg , msh , var);
        applyNodeIsothermalWallPin(cfg , cuda_cfg , msh , var);
        floorEventsSetPhase(FloorEventPhase::Audit, lastStep);
        dependentVariables_d_wrapper(cfg , cuda_cfg , msh , var);
        if (g.injStep > 0) injT = var.c_d["T"];
        if (injT != nullptr) {
            flow_float t = 0.0f;
            gpuErrchk( cudaMemcpy(&t, injT + g.injNode, sizeof(flow_float), cudaMemcpyDeviceToHost) );
            std::printf("[floor-events] test audit copy: T[%lld] = %.9g K (T_floor %.9g K)\n", g.injNode, (double)t, g.tFloor);
        }
    } catch (...) {
        restore();
        throw;
    }
    restore();
    floorEventsSetPhase(FloorEventPhase::Aux, lastStep);
    g.auditStatus = "done";
    if (verify) {
        long nArr = 0, nDiff = 0; std::string first;
        for (auto& kv : before) {
            std::vector<flow_float> a(nca);
            gpuErrchk( cudaMemcpy(a.data(), var.c_d[kv.first], nca*sizeof(flow_float), cudaMemcpyDeviceToHost) );
            ++nArr;
            if (std::memcmp(a.data(), kv.second.data(), nca*sizeof(flow_float)) != 0) { ++nDiff; if (first.empty()) first = kv.first; }
        }
        std::printf("[floor-events] audit verify: %ld cell arrays compared byte-wise before/after the audit, %ld changed%s%s\n",
                    nArr, nDiff, nDiff ? " (first: " : "", nDiff ? (first + ")").c_str() : "");
    }
}

void floorEventsClose(int lastStep)
{
    if (!g.on) return;
    std::ostringstream row;
    row << "session_end," << g.session << "," << lastStep << ",-1,-1," << g.path << "," << g.nReal << "," << g.nGhost
        << ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,,,,0,audit=" << g.auditStatus;
    writeLine(row.str());
    std::fclose(g.fp); g.fp = nullptr;
    cudaFree(g.acc_d); g.acc_d = nullptr;
    if (g.hT_d) { cudaFree(g.hT_d); g.hT_d = nullptr; }
    g.on = false;
}

void floorEventsTestInject(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step)
{
    (void)cuda_cfg; (void)msh;
    if (!g.on || g.injStep != step) return;
    double* out_d = nullptr;
    gpuErrchk( cudaMalloc(&out_d, 6*sizeof(double)) );
    auto ptr = [&](const char* k) -> flow_float* { auto it = var.c_d.find(k); return it == var.c_d.end() ? nullptr : it->second; };
    if (g.injMode == 1) {
        floorEventsInject_d<<<1, 1>>>((geom_int)g.injNode, cfg.thermalMethod, cfg.gamma, cfg.cp, cfg.roMin,
                                     thermo_species_device_ptr(), cfg.nSpecies, species_roY_device_ptr(), g.hT_d, g.tFloor, g.injValue,
                                     ptr("ro"), ptr("roUx"), ptr("roUy"), ptr("roUz"), ptr("roe"), ptr("roeN"), out_d);
    } else {
        floorEventsInjectScale_d<<<1, 1>>>((geom_int)g.injNode, g.injValue, cfg.nSpecies, species_roY_device_ptr(),
                                          ptr("ro"), ptr("roUx"), ptr("roUy"), ptr("roUz"), ptr("roe"), ptr("roK"), ptr("roOmega"), out_d);
    }
    gpuErrchk( cudaPeekAtLastError() );
    double o[6];
    gpuErrchk( cudaMemcpy(o, out_d, sizeof(o), cudaMemcpyDeviceToHost) );
    cudaFree(out_d);
    if (g.injMode == 1)
        std::printf("[floor-events] test inject: step %d node %lld mode de roe %.9g -> %.9g (e_floor %.17g J/kg, cv(T_floor) %.9g J/kg/K, "
                    "rho %.9g, ek %.9g, de %.9g)\n", step, g.injNode, o[0], o[1], o[2], o[3], o[4], o[5], g.injValue);
    else
        std::printf("[floor-events] test inject: step %d node %lld mode rhoscale %.9g: ro %.9g -> %.9g, roe %.9g -> %.9g\n",
                    step, g.injNode, g.injValue, o[0], o[2], o[1], o[3]);
}
