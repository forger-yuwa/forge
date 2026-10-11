// 診断 FORGE_DIAG_COMMIT_LOSS=<N> (plans/active/architecture-float-state-double-geometry.md §4.5・§6 V4 の記録、既定 off)。
// 何を比べるかは commitLossDiag.hpp の説明。ここは集計のカーネルとホスト側の段取り・出力。

#include "cuda_forge/commitLossDiag.hpp"
#include "cuda_forge/cudaWrapper.cuh"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace commitLossDiag {

namespace detail {
bool g_enabled = false;
bool g_active  = false;
}

// ビット単位の一致 (+0 と −0 も区別する)。float のビルドでは double 版を使わないので、未使用の警告が出ない無名名前空間の外に置く。
__device__ inline bool sameBits(float a, float b)   { return __float_as_uint(a) == __float_as_uint(b); }
__device__ inline bool sameBits(double a, double b) { return __double_as_longlong(a) == __double_as_longlong(b); }

namespace {

constexpr int kNAcc = 6;   // n_req, n_lost, n_clip, S_req, S_act, S_lost
constexpr int kNReg = 2;   // 0 = interior (bnode_flag 0)、1 = boundary
constexpr int kStride = kNReg * kNAcc;
constexpr int kQK = 5, kQOmega = 6, kQY0 = 7;   // 量の番号: 0..4 = ρ..ρE、5 = ρk、6 = ρω、7.. = ρY_s

int g_N = 0;
int g_step = 0;            // 計っている step (iStep + 1)
geom_int g_nCells = 0;
const geom_int* g_bnode = nullptr;   // msh.bnode_flag_d (領域の分け方。借りるだけ)
std::vector<std::string> g_names;
std::vector<char> g_measured;        // この step で計った量
double* g_acc_d = nullptr;           // [nQ][kNReg][kNAcc]
flow_float* g_sstBefore_d[2] = {nullptr, nullptr};
flow_float* g_sstDq_d[2]     = {nullptr, nullptr};
double* g_spDq_d = nullptr;
std::ofstream g_csv;
bool g_noteFlowSkip = false;         // ρ..ρE を計らない理由を 1 回だけ書く

// 1 量ぶんの集計。warp ごとに足してから 1 回だけ atomicAdd する (毎 step ではなく N の倍数の step だけ起動する)。
// 床 (hasFloor) に当たる節点は Q_before + dq_req < 床 で判定する (カーネルの max(Q_before + dq, 床) が床を選ぶ条件)。
template <typename TD>
__global__ void commit_loss_reduce_d(geom_int nCells, const geom_int* bnode,
                                     const flow_float* before, const flow_float* after, const TD* dq,
                                     int hasFloor, double floorv, double* acc)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    double v[kStride];
    #pragma unroll
    for (int k = 0; k < kStride; ++k) v[k] = 0.0;
    if (ic < nCells) {
        const double d = static_cast<double>(dq[ic]);
        if (d != 0.0) {
            const int r = (bnode != nullptr && bnode[ic] != 0) ? 1 : 0;
            const flow_float qb = before[ic];
            const flow_float qa = after[ic];
            if (hasFloor != 0 && static_cast<double>(qb) + d < floorv) {
                v[r*kNAcc + 2] = 1.0;                                   // n_clip
            } else {
                const double act = static_cast<double>(qa) - static_cast<double>(qb);
                v[r*kNAcc + 0] = 1.0;                                   // n_req
                v[r*kNAcc + 1] = sameBits(qa, qb) ? 1.0 : 0.0;          // n_lost
                v[r*kNAcc + 3] = fabs(d);                               // S_req
                v[r*kNAcc + 4] = fabs(act);                             // S_act
                v[r*kNAcc + 5] = fabs(d - act);                         // 丸めで失われた分
            }
        }
    }
    #pragma unroll
    for (int k = 0; k < kStride; ++k) {
        double x = v[k];
        for (int off = 16; off > 0; off >>= 1) x += __shfl_down_sync(0xffffffffu, x, off);
        if ((threadIdx.x & 31) == 0 && x != 0.0) atomicAdd(&acc[k], x);
    }
}

template <typename TD>
void reduceQty(int iq, const flow_float* before, const flow_float* after, const TD* dq, int hasFloor, double floorv)
{
    if (iq < 0 || iq >= (int)g_names.size() || g_nCells <= 0) return;
    const int bs = 256;   // warp の和を取るので 32 の倍数
    const int nb = (int)((g_nCells + bs - 1) / bs);
    commit_loss_reduce_d<TD><<<nb, bs>>>(g_nCells, g_bnode, before, after, dq, hasFloor, floorv, g_acc_d + (size_t)iq * kStride);
    gpuErrchk( cudaPeekAtLastError() );
    g_measured[iq] = 1;
}

flow_float* devAlloc(geom_int n)
{
    flow_float* p = nullptr;
    gpuErrchk( cudaMalloc((void**)&p, (size_t)n * sizeof(flow_float)) );
    return p;
}

}  // namespace

void init(solverConfig& cfg, mesh& msh, variables& var)
{
    const char* e = std::getenv("FORGE_DIAG_COMMIT_LOSS");
    if (e == nullptr || *e == '\0') return;
    errno = 0;
    char* end = nullptr;
    const long n = std::strtol(e, &end, 10);
    if (errno != 0 || end == e || *end != '\0' || n < 0 || n > 1000000000L) {
        fprintf(stderr, "[commitLoss] FORGE_DIAG_COMMIT_LOSS='%s' は 0 以上の整数でない (N > 0 で有効)。止める\n", e);
        std::exit(EXIT_FAILURE);
    }
    if (n == 0) {
        printf("[commitLoss] FORGE_DIAG_COMMIT_LOSS=0: 無効\n");
        return;
    }
    // 計るのは定常の陰解法の commit だけ (main の advanceImplicitSteady が step の頭と終わりを知らせる)。
    // dual-time (applyBlockImplicitCorrectionInPlace) と陽解法・CPU の経路にはフックが無いので、黙って何も計らない run にしない。
    if (cfg.gpu != 1 || cfg.isImplicit != 1 || cfg.unsteady != 0) {
        fprintf(stderr, "[commitLoss] 対象外の経路: 定常の陰解法 (gpu 1・timeIntegration 11・unsteady 0) の commit だけを計る "
                        "(gpu %d, isImplicit %d, unsteady %d)。止める\n", cfg.gpu, cfg.isImplicit, cfg.unsteady);
        std::exit(EXIT_FAILURE);
    }

    g_N = (int)n;
    g_nCells = msh.nCells;
    g_bnode = msh.bnode_flag_d;
    g_names = {"ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"};
    const int nSp = (var.nSpeciesRegistered >= 2) ? var.nSpeciesRegistered : 0;
    for (int s = 0; s < nSp; ++s) g_names.push_back("roY" + std::to_string(s));
    g_measured.assign(g_names.size(), 0);

    gpuErrchk( cudaMalloc((void**)&g_acc_d, g_names.size() * kStride * sizeof(double)) );
    gpuErrchk( cudaMemset(g_acc_d, 0, g_names.size() * kStride * sizeof(double)) );
    for (int k = 0; k < 2; ++k) { g_sstBefore_d[k] = devAlloc(g_nCells); g_sstDq_d[k] = devAlloc(g_nCells); }
    if (nSp > 0) gpuErrchk( cudaMalloc((void**)&g_spDq_d, (size_t)g_nCells * sizeof(double)) );

    g_csv.open("commit_loss.csv", std::ios::out | std::ios::trunc);
    if (!g_csv) {
        fprintf(stderr, "[commitLoss] commit_loss.csv を開けない。止める\n");
        std::exit(EXIT_FAILURE);
    }
    g_csv << "step,qty,n_req,n_lost,S_req,S_act,S_lostfrac,n_clip,region\n";
    g_csv.flush();
    detail::g_enabled = true;

    printf("[commitLoss] 有効: step (= iStep + 1、res_<step>.h5 と同じ) が %d の倍数の step で、定常の陰解法の commit の直後に"
           "要求した更新 dq_req と Q_after − Q_before を比べる → commit_loss.csv\n", g_N);
    printf("[commitLoss]   ρ..ρE: %s commit (Q_before = roN、dq_req = implicitRelax %.6g を掛けた最終の補正)\n",
           (cfg.blockDPLUR == 1) ? "block DPLUR" : "スカラー DPLUR", (double)cfg.implicitRelax);
    if (cfg.qAccumulatorFP64 == 1)
        printf("[commitLoss]   ρ..ρE は計らない: qAccumulatorFP64 (commit は FP64 の正本に積む別の経路)\n");
    else if (cfg.blockDPLUR == 1 && cfg.updateGuardAlpha > (flow_float)0.0)
        printf("[commitLoss]   ρ..ρE は計らない: updateGuardAlpha %.6g > 0 (縮小率 s を写していないので dq_req が決まらない)\n",
               (double)cfg.updateGuardAlpha);
    printf("[commitLoss]   ρk・ρω: SST の point-implicit (Q_before = カーネル直前の値、dq_req = dk・dw、床 0 / 1e-20 は n_clip)\n");
    if (nSp > 0) {
        if (cfg.speciesImplicitCoupling == 2)
            printf("[commitLoss]   ρY_0..%d は計らない: speciesImplicitCoupling 2 (EOS 結合の最終 commit)\n", nSp - 1);
        else if (cfg.speciesImplicitCoupling == 1)
            printf("[commitLoss]   ρY_0..%d: scalar DPLUR の commit (Q_before = roY{s}N、dq_req = dq_roY{s}_old、床 0 は n_clip)\n", nSp - 1);
        else
            printf("[commitLoss]   ρY_0..%d: point-implicit の commit (Q_before = roY{s}N、dq_req = speciesImplicitRelax %.6g × δ、"
                   "床 0 は n_clip)\n", nSp - 1, (double)cfg.speciesImplicitRelax);
    }
    printf("[commitLoss]   比べる時点: commit のカーネルの直後 (等温壁のピン・SST の E_t 補正・化学種の再正規化・EOS の床の前)\n");
    printf("[commitLoss]   ログの書式: 量 n_req/n_lost S_act/S_req (interior = どの境界にも属さない節点)。boundary の節点と n_clip は CSV\n");
}

void beginStep(int iStep)
{
    if (!detail::g_enabled) return;
    g_step = iStep + 1;
    if (g_step % g_N != 0) return;
    gpuErrchk( cudaMemset(g_acc_d, 0, g_names.size() * kStride * sizeof(double)) );
    g_measured.assign(g_names.size(), 0);
    detail::g_active = true;
}

void endStep()
{
    if (!detail::g_active) return;
    detail::g_active = false;
    std::vector<double> h(g_names.size() * kStride, 0.0);
    gpuErrchk( cudaMemcpy(h.data(), g_acc_d, h.size() * sizeof(double), cudaMemcpyDeviceToHost) );

    static const char* kRegion[kNReg] = {"interior", "boundary"};
    std::ostringstream line;
    line << "[commitLoss] step " << g_step;
    for (size_t iq = 0; iq < g_names.size(); ++iq) {
        if (!g_measured[iq]) continue;
        for (int r = 0; r < kNReg; ++r) {
            const double* a = h.data() + iq * kStride + r * kNAcc;
            const long long nReq = (long long)a[0], nLost = (long long)a[1], nClip = (long long)a[2];
            const double sReq = a[3], sAct = a[4], sLost = a[5];
            g_csv << g_step << "," << g_names[iq] << "," << nReq << "," << nLost << ","
                  << std::setprecision(17) << sReq << "," << sAct << ",";
            if (sReq > 0.0) g_csv << (sLost / sReq); else g_csv << "nan";
            g_csv << "," << nClip << "," << kRegion[r] << "\n";
            if (r == 0) {
                line << " " << g_names[iq] << " " << nReq << "/" << nLost << " ";
                if (sReq > 0.0) line << std::setprecision(4) << (sAct / sReq); else line << "-";
            }
        }
    }
    g_csv.flush();
    printf("%s\n", line.str().c_str());
    fflush(stdout);
}

void finalize()
{
    if (!detail::g_enabled) return;
    detail::g_enabled = false;
    detail::g_active = false;
    if (g_csv.is_open()) g_csv.close();
    cudaFree(g_acc_d); g_acc_d = nullptr;
    for (int k = 0; k < 2; ++k) {
        cudaFree(g_sstBefore_d[k]); g_sstBefore_d[k] = nullptr;
        cudaFree(g_sstDq_d[k]);     g_sstDq_d[k] = nullptr;
    }
    if (g_spDq_d != nullptr) { cudaFree(g_spDq_d); g_spDq_d = nullptr; }
}

void flowCommit(solverConfig& cfg, mesh& msh, variables& var, flow_float* const* dq, flow_float guardAlpha)
{
    (void)msh;
    if (!detail::g_active) return;
    // qAccumulatorFP64: commit は FP64 の正本に積み、FP32 はその丸めの写しになる (基準が roN でない)。
    // updateGuardAlpha > 0: commit は s·dq を足すが s を写していない。どちらも dq_req が決まらないので計らない (起動時に書いた)。
    if (var.qacc_d[0] != nullptr || guardAlpha > (flow_float)0.0) {
        if (!g_noteFlowSkip) {
            printf("[commitLoss] ρ..ρE は計らない (%s)\n", (var.qacc_d[0] != nullptr) ? "qAccumulatorFP64" : "updateGuardAlpha > 0");
            g_noteFlowSkip = true;
        }
        return;
    }
    (void)cfg;
    static const char* kBefore[5] = {"roN", "roUxN", "roUyN", "roUzN", "roeN"};
    static const char* kAfter[5]  = {"ro", "roUx", "roUy", "roUz", "roe"};
    for (int k = 0; k < 5; ++k)
        reduceQty<flow_float>(k, var.c_d.at(kBefore[k]), var.c_d.at(kAfter[k]), dq[k], 0, 0.0);
}

void sstBefore(mesh& msh, variables& var, flow_float** dk, flow_float** dw)
{
    (void)msh;
    *dk = nullptr; *dw = nullptr;
    if (!detail::g_active) return;
    const size_t bytes = (size_t)g_nCells * sizeof(flow_float);
    gpuErrchk( cudaMemcpy(g_sstBefore_d[0], var.c_d.at("roK"),     bytes, cudaMemcpyDeviceToDevice) );
    gpuErrchk( cudaMemcpy(g_sstBefore_d[1], var.c_d.at("roOmega"), bytes, cudaMemcpyDeviceToDevice) );
    *dk = g_sstDq_d[0];
    *dw = g_sstDq_d[1];
}

void sstAfter(mesh& msh, variables& var)
{
    (void)msh;
    if (!detail::g_active) return;
    // 床は applySSTPointImplicit_d の max(·, 0) と max(·, 1e-20) (flow_float の値)
    reduceQty<flow_float>(kQK,     g_sstBefore_d[0], var.c_d.at("roK"),     g_sstDq_d[0], 1,
                          static_cast<double>(static_cast<flow_float>(0.0)));
    reduceQty<flow_float>(kQOmega, g_sstBefore_d[1], var.c_d.at("roOmega"), g_sstDq_d[1], 1,
                          static_cast<double>(static_cast<flow_float>(1.0e-20)));
}

double* speciesDqBuffer()
{
    return detail::g_active ? g_spDq_d : nullptr;
}

void speciesCommitPointImplicit(mesh& msh, int s, const flow_float* roYN, const flow_float* roY, const double* dq, double floor)
{
    (void)msh;
    if (!detail::g_active || dq == nullptr) return;
    reduceQty<double>(kQY0 + s, roYN, roY, dq, 1, floor);
}

void speciesCommitDPLUR(mesh& msh, int s, const flow_float* roYN, const flow_float* roY, const flow_float* dq, double floor)
{
    (void)msh;
    if (!detail::g_active || dq == nullptr) return;
    reduceQty<flow_float>(kQY0 + s, roYN, roY, dq, 1, floor);
}

}  // namespace commitLossDiag
