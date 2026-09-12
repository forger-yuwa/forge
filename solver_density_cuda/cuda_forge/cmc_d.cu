// =============================================================================
// 混合分率インフラ (Phase A of methods/chemistry_cmc.md): Bilger ξ̃ 診断、分散 ξ''² 輸送、χ̃ 診断。
// =============================================================================
#include "cmc_d.cuh"
#include "scalarTransport_d.cuh"
#include "chemistry_f32_d.cuh"
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>

// thrust::identity は CCCL 3 (CUDA 13) で削除されたので自前の述語 (CUDA 12/13 両対応)
struct CmcFlagNonZero { __host__ __device__ bool operator()(int f) const { return f != 0; } };
#include "chemistrySource_d.cuh"
#include "chemistry_d.cuh"
#include <vector>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <map>
#include <string>
#include <stdexcept>

namespace {

// Bilger 係数: β = Σ_s c_s Y_s、ξ = (β - β_O)/(β_F - β_O)
struct BilgerCoef {
    int nSpecies;
    double c[THERMO_MAX_SPECIES];
    double betaF, betaO;
};
__constant__ BilgerCoef g_bilger;
static BilgerCoef       g_bilger_host;
static flow_float**     g_Y_dev = nullptr;   // 原始質量分率 Y{s} の device ポインタ配列
static bool             g_ready = false;

__global__ void cmc_mixfrac_d(geom_int nCells, flow_float** Y, flow_float* xi)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    double beta = 0.0;
    for (int s = 0; s < g_bilger.nSpecies; ++s) beta += g_bilger.c[s] * (double)Y[s][ic];
    double z = (beta - g_bilger.betaO) / (g_bilger.betaF - g_bilger.betaO);
    z = fmin(fmax(z, 0.0), 1.0);
    xi[ic] = (flow_float)z;
}

// 1 スカラの Green-Gauss 勾配 (ransTransport の k/ω 版と同じ ghostless 規約)
__global__ void cmc_gradient_face_d(geom_int nPlanes, geom_int nCells, int isNode, geom_int* plane_cells,
                                    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz,
                                    flow_float* phi, flow_float* dpdx, flow_float* dpdy, flow_float* dpdz)
{
    const geom_int ip = blockDim.x * blockIdx.x + threadIdx.x;
    if (ip >= nPlanes) return;
    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];
    const geom_float f = fx[ip];
    flow_float pf;
    if (isNode != 0 && ic1 >= nCells) pf = phi[ic0];
    else                              pf = f * phi[ic0] + (1.0 - f) * phi[ic1];
    const geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip];
    atomicAdd(&dpdx[ic0],  sxx * pf); atomicAdd(&dpdy[ic0],  syy * pf); atomicAdd(&dpdz[ic0],  szz * pf);
    atomicAdd(&dpdx[ic1], -sxx * pf); atomicAdd(&dpdy[ic1], -syy * pf); atomicAdd(&dpdz[ic1], -szz * pf);
}
__global__ void cmc_gradient_div_vol_d(geom_int nCells, geom_float* vol, flow_float* dpdx, flow_float* dpdy, flow_float* dpdz)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const geom_float v = vol[ic];
    dpdx[ic] /= v; dpdy[ic] /= v; dpdz[ic] /= v;
}

// 分散のソース: 生成 2 μ_t/Sc_t |∇ξ̃|²、散逸 ρχ̃ = cChi β* ω (ρξ''²) (src_jac に対角)。χ̃ を診断出力。
__global__ void cmc_variance_source_d(geom_int nCells, geom_float* vol, flow_float* ro, flow_float* vis_turb, flow_float* omega,
                                      flow_float* dXidx, flow_float* dXidy, flow_float* dXidz,
                                      flow_float* roXiVar, flow_float* res_roXiVar, flow_float* src_jac,
                                      flow_float* chi, double Sc_t, double cChi, double betaStar)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double g2   = (double)dXidx[ic]*dXidx[ic] + (double)dXidy[ic]*dXidy[ic] + (double)dXidz[ic]*dXidz[ic];
    const double prod = 2.0 * (double)vis_turb[ic] / Sc_t * g2;                 // [kg/m3/s]
    const double rate = cChi * betaStar * fmax((double)omega[ic], 0.0);          // [1/s]
    const double rxv  = fmax((double)roXiVar[ic], 0.0);
    const double diss = rate * rxv;                                              // [kg/m3/s]
    res_roXiVar[ic] += (flow_float)((prod - diss) * (double)vol[ic]);
    src_jac[ic]     += (flow_float)rate;
    chi[ic] = (flow_float)(rate * rxv / fmax((double)ro[ic], 1.0e-30));         // χ̃ [1/s]
}

__global__ void cmc_primitive_d(geom_int nCells_all, flow_float* ro, flow_float* xi, flow_float* roXiVar, flow_float* xiVar)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells_all) return;
    const double r  = fmax((double)ro[ic], 1.0e-30);
    const double z  = (double)xi[ic];
    const double vmax = r * fmax(z * (1.0 - z), 0.0);                            // 実現可能性: ξ''² ≤ ξ̃(1-ξ̃)
    double rv = fmin(fmax((double)roXiVar[ic], 0.0), vmax);
    roXiVar[ic] = (flow_float)rv;
    xiVar[ic]   = (flow_float)(rv / r);
}

__global__ void cmc_dirichlet_zero_d(geom_int nb, int isNode, geom_int* bplane_cell, geom_int* bplane_cell_ghst,
                                     flow_float* rophi, flow_float* phi)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    const geom_int ig = bplane_cell_ghst[ib];
    rophi[ig] = 0.0; phi[ig] = 0.0;
    if (isNode != 0) { const geom_int ic = bplane_cell[ib]; rophi[ic] = 0.0; phi[ic] = 0.0; }   // 境界ノードへピン
}
__global__ void cmc_neumann_d(geom_int nb, geom_int* bplane_cell, geom_int* bplane_cell_ghst, flow_float* rophi, flow_float* phi)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    const geom_int ic = bplane_cell[ib], ig = bplane_cell_ghst[ib];
    rophi[ig] = rophi[ic]; phi[ig] = phi[ic];
}

ScalarTransportDesc buildVarianceDesc(variables& var, const solverConfig& cfg)
{
    ScalarTransportDesc d{};
    d.phi = var.c_d.at("xiVar");
    d.dphidx = var.c_d.at("dXiVardx"); d.dphidy = var.c_d.at("dXiVardy"); d.dphidz = var.c_d.at("dXiVardz");
    d.rho_phi = var.c_d.at("roXiVar"); d.rho_phi_N = var.c_d.at("roXiVarN"); d.rho_phi_M = var.c_d.at("roXiVarM");
    d.res_rho_phi = var.c_d.at("res_roXiVar"); d.res_rho_phi_m = var.c_d.at("res_roXiVar_m");
    d.src_jac = var.c_d.at("src_jac_xiVar"); d.transport_diag = var.c_d.at("transport_diag_xiVar");
    d.floor = static_cast<flow_float>(0.0);
    d.sigma = static_cast<flow_float>(1.0 / cfg.Sc_t);   // 有効拡散 = μ + μ_t/Sc_t
    d.diffusion = 1;
    return d;
}

void gradientOf(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                const char* phi, const char* gx, const char* gy, const char* gz)
{
    flow_float* grad_vol = (cfg.isAxisymmetric == 1) ? var.c_d.at("A_planar") : var.c_d.at("volume");
    flow_float* grad_sx  = (cfg.isAxisymmetric == 1) ? var.p_d.at("sx_planar") : var.p_d.at("sx");
    flow_float* grad_sy  = (cfg.isAxisymmetric == 1) ? var.p_d.at("sy_planar") : var.p_d.at("sy");
    flow_float* grad_sz  = (cfg.isAxisymmetric == 1) ? var.p_d.at("sz_planar") : var.p_d.at("sz");
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at(gx), 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at(gy), 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at(gz), 0, msh.nCells_all * sizeof(flow_float)));
    cmc_gradient_face_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nPlanes, msh.nCells, (cfg.discretization == "node") ? 1 : 0, msh.map_plane_cells_d,
        var.p_d.at("fx"), grad_sx, grad_sy, grad_sz, var.c_d.at(phi), var.c_d.at(gx), var.c_d.at(gy), var.c_d.at(gz));
    cmc_gradient_div_vol_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells, grad_vol, var.c_d.at(gx), var.c_d.at(gy), var.c_d.at(gz));
}

} // namespace

void cmcInit_d(solverConfig& cfg, variables& var)
{
    g_ready = false;
    if (!cmcMixfracEnabled(cfg, var)) return;
    const ReactionTable* rt = chemistry_table_host();
    if (rt == nullptr || rt->nElem == 0)
        throw std::runtime_error("[cmc] mixfrac needs element composition: mechanismFile must have species[].composition");
    if (rt->elemH < 0 && rt->elemC < 0)
        throw std::runtime_error("[cmc] Bilger mixture fraction needs H or C in the mechanism elements");

    // 種の分子量 (元素組成から) と Bilger 係数 c_s = Σ_e coef_e · atoms[s][e] W_e / W_s
    const int ns = rt->nSpecies;
    std::vector<double> W(ns, 0.0);
    for (int s = 0; s < ns; ++s) for (int e = 0; e < rt->nElem; ++e) W[s] += rt->atoms[s][e] * rt->elemW[e];
    auto coefOf = [&](int e) -> double {
        if (e == rt->elemC) return 2.0 / rt->elemW[e];
        if (e == rt->elemH) return 0.5 / rt->elemW[e];
        if (e == rt->elemO) return -1.0 / rt->elemW[e];
        return 0.0;
    };
    BilgerCoef& b = g_bilger_host;
    b.nSpecies = ns;
    for (int s = 0; s < ns; ++s) {
        double c = 0.0;
        if (W[s] > 0.0) for (int e = 0; e < rt->nElem; ++e) c += coefOf(e) * rt->atoms[s][e] * rt->elemW[e] / W[s];
        b.c[s] = c;
        if (W[s] <= 0.0) std::printf("[cmc] WARNING: species %s has no composition in the mechanism (treated as inert for xi)\n", cfg.speciesNames[s].c_str());
    }
    // 燃料/酸化剤流の β (モル分率 → 質量分率)
    auto betaOfStream = [&](const std::map<std::string, double>& X, const char* label) -> double {
        double mmix = 0.0;
        for (const auto& kv : X) {
            int s = -1; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i;
            if (s < 0) throw std::runtime_error(std::string("[cmc] mixfrac.") + label + ": species '" + kv.first + "' not in physProp.species");
            mmix += kv.second * W[s];
        }
        double beta = 0.0;
        for (const auto& kv : X) {
            int s = 0; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i;
            beta += b.c[s] * (kv.second * W[s] / mmix);
        }
        return beta;
    };
    b.betaF = betaOfStream(cfg.chemMixfracFuelX, "fuelX");
    b.betaO = betaOfStream(cfg.chemMixfracOxidX, "oxidizerX");
    if (std::fabs(b.betaF - b.betaO) < 1.0e-12) throw std::runtime_error("[cmc] fuelX and oxidizerX give the same Bilger beta");
    gpuErrchk( cudaMemcpyToSymbol(g_bilger, &b, sizeof(BilgerCoef)) );

    std::vector<flow_float*> hY(ns);
    for (int s = 0; s < ns; ++s) hY[s] = var.c_d.at("Y" + std::to_string(s));
    if (g_Y_dev) { cudaFree(g_Y_dev); g_Y_dev = nullptr; }
    gpuErrchk( cudaMalloc((void**)&g_Y_dev, ns * sizeof(flow_float*)) );
    gpuErrchk( cudaMemcpy(g_Y_dev, hY.data(), ns * sizeof(flow_float*), cudaMemcpyHostToDevice) );
    g_ready = true;
    std::printf("cmcInit_d: Bilger beta_F=%.6g beta_O=%.6g (nElem=%d, cChi=%.3g)\n", b.betaF, b.betaO, rt->nElem, cfg.chemCchi);
}

void cmcMixfrac_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    cmc_mixfrac_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, g_Y_dev, var.c_d.at("xi"));
    gradientOf(cfg, cuda_cfg, msh, var, "xi", "dXidx", "dXidy", "dXidz");
    gradientOf(cfg, cuda_cfg, msh, var, "xiVar", "dXiVardx", "dXiVardy", "dXiVardz");
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcVarianceTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at("res_roXiVar"), 0, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at("transport_diag_xiVar"), 0, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at("src_jac_xiVar"), 0, msh.nCells * sizeof(flow_float)));
    const ScalarTransportDesc desc = buildVarianceDesc(var, cfg);
    scalarTransportResidual_d(cfg, cuda_cfg, msh, var, desc);
    const double betaStar = 0.09;
    cmc_variance_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        (msh.volumePartial_d != nullptr) ? msh.volumePartial_d : var.c_d.at("volume"),
        var.c_d.at("ro"), var.c_d.at("vis_turb"), var.c_d.at("omega"),
        var.c_d.at("dXidx"), var.c_d.at("dXidy"), var.c_d.at("dXidz"),
        var.c_d.at("roXiVar"), var.c_d.at("res_roXiVar"), var.c_d.at("src_jac_xiVar"), var.c_d.at("chi"),
        (double)cfg.Sc_t, cfg.chemCchi, betaStar);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcVarianceTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    const ScalarTransportDesc desc = buildVarianceDesc(var, cfg);
    scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    // xi は最新の Y から (species primitive の後に呼ぶこと)
    cmc_mixfrac_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, g_Y_dev, var.c_d.at("xi"));
    cmc_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d.at("ro"), var.c_d.at("xi"), var.c_d.at("roXiVar"), var.c_d.at("xiVar"));
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void applyCmcBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    const int isNode = (cfg.discretization == "node") ? 1 : 0;
    for (auto& bc : msh.bconds) {
        if (bc.iPlanes.empty()) continue;
        const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
        const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
        if (isInlet) {
            cmc_dirichlet_zero_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(nb, isNode, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                                                                                 var.c_d.at("roXiVar"), var.c_d.at("xiVar"));
        } else {
            cmc_neumann_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(nb, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                                                                          var.c_d.at("roXiVar"), var.c_d.at("xiVar"));
        }
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!g_ready) return;
    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    gpuErrchk( cudaMemcpy(var.c_d.at("roXiVarN"), var.c_d.at("roXiVar"), bytes, cudaMemcpyDeviceToDevice) );
    gpuErrchk( cudaMemcpy(var.c_d.at("roXiVarM"), var.c_d.at("roXiVar"), bytes, cudaMemcpyDeviceToDevice) );
}

void cmcUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!g_ready) return;
    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    gpuErrchk( cudaMemcpy(var.c_d.at("roXiVarM"), var.c_d.at("roXiVar"), bytes, cudaMemcpyDeviceToDevice) );
}

// =============================================================================
// Phase B/C: 条件付きスカラー Q(η) — 混合線初期化・スライス輸送・η 拡散 (AMC) + 化学 (点陰解)・β-PDF 平均ソース
// =============================================================================
#include "thermo_d.cuh"
#include "speciesTransport_d.cuh"

namespace {

#define CMC_MAX_ETA 129

struct CmcParams {
    int nEta, nSpecies, chemOn, couple, kSt;
    double pdfFloor, Tmax, Tfreeze, dtScale, relax, alpha, dTmax, dTgate;
    double eta[CMC_MAX_ETA];       // η 格子
    double w[CMC_MAX_ETA];         // 台形重み (Σ=1)
    double G[CMC_MAX_ETA];         // AMC 形状 exp(-2 erfinv(2η-1)^2) (端は 0)
    double YF[THERMO_MAX_SPECIES], YO[THERMO_MAX_SPECIES];   // 燃料/酸化剤流の質量分率
    double hF, hO;                 // 燃料/酸化剤流の顕エンタルピー [J/kg]
};
__constant__ CmcParams g_cmc;
static CmcParams g_cmc_host;
static bool g_q_ready = false;
static int  g_nSlice = 0;
static geom_int g_nCellsAll = 0, g_nCells = 0;
static flow_float *g_Q = nullptr, *g_roQ = nullptr, *g_resQ = nullptr, *g_diagQ = nullptr, *g_sjQ = nullptr, *g_resmQ = nullptr;
static flow_float *g_gdum = nullptr;                       // 勾配ダミー (拡散カーネルは 1 次で勾配を使わない)
static flow_float *g_omega = nullptr, *g_qdot = nullptr, *g_jac = nullptr;   // PDF 平均ソース
static flow_float *g_ypdf = nullptr, *g_hpdf = nullptr, *g_tau = nullptr;   // couple 2: PDF 積分 Ỹ_s, h̃, 緩和時間 τ
static flow_float *g_dpdf = nullptr, *g_dpdf_prev = nullptr;
static flow_float *g_qrel = nullptr;
static flow_float *g_qdebt = nullptr;
static flow_float *g_gate = nullptr;    // couple 6/7: 燃焼領域ゲート g = clip((h_pdf − h_line(ξ_Ω))/(1500·dTgate), 0, 1)
static flow_float *g_qcap = nullptr;    // couple 5: リミッタ後の発熱 [J/m3] を陰的化学ソース経路 (res_roe) へ渡す   // couple 4: 1 step の ΔT 上限で出し切れなかった発熱の持ち越し [J/m3]   // couple 4: 条件付き点陰解で実際に進んだ反応の PDF 平均発熱 [J/kg] (この step)   // couple 4: PDF 積分の反応離脱 D_pdf,s = Ỹ_pdf,s − line_s(ξ_Ω) と前ステップ値
static int g_cmcStep = 0;
// 分割版 cmc_step の作業配列 (double): [η][node] の Ω, <χ|η>, T, Ω·q_rel, Ω·Q̇ と [s][node] の Ỹ_pdf, h̃_pdf, ω̄/J̄ 集計
static double *g_OmD = nullptr, *g_chiEtaD = nullptr, *g_TEtaD = nullptr, *g_qrelEtaD = nullptr, *g_QdEtaD = nullptr;
static double *g_ypdfD = nullptr, *g_hpdfD = nullptr, *g_omegaAccD = nullptr, *g_jacAccD = nullptr;
__constant__ SpeciesThermoF g_spf[THERMO_MAX_SPECIES];   // fp32 化学 (cmc.fp32) の係数表
__constant__ ReactionTableF g_rtf;
static int g_cmcFp32 = 1;

__device__ __forceinline__ size_t sl(int var, int k) { return (size_t)(var * g_cmc.nEta + k); }

// 混合線: Y(η) = η Y_F + (1-η) Y_O, h(η) = η h_F + (1-η) h_O
__global__ void cmc_q_init_d(geom_int nCells_all, const flow_float* ro, flow_float* Q, flow_float* roQ)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells_all) return;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    const double r = (double)ro[ic];
    for (int k = 0; k < ne; ++k) {
        const double e = g_cmc.eta[k];
        for (int s = 0; s < ns; ++s) {
            const double y = e * g_cmc.YF[s] + (1.0 - e) * g_cmc.YO[s];
            Q[sl(s, k) * nCells_all + ic] = (flow_float)y; roQ[sl(s, k) * nCells_all + ic] = (flow_float)(r * y);
        }
        const double h = e * g_cmc.hF + (1.0 - e) * g_cmc.hO;
        Q[sl(ns, k) * nCells_all + ic] = (flow_float)h; roQ[sl(ns, k) * nCells_all + ic] = (flow_float)(r * h);
    }
}

// 燃料/酸化剤流の顕エンタルピーを device thermo で評価 (1 スレッド)
__global__ void cmc_stream_enthalpy_d(const SpeciesThermo* sp, int ns, double TF, double TO, double* out)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    out[0] = thermo_h_mix(sp, ns, g_cmc.YF, TF);
    out[1] = thermo_h_mix(sp, ns, g_cmc.YO, TO);
}

// 境界: 入口は混合線 (ghost + node ピン)、他は ghost = owner。grid: (bplane × slice)
__global__ void cmc_q_boundary_d(geom_int nb, int nSlice, geom_int nCells_all, int isNode, int isInlet,
                                 const geom_int* bplane_cell, const geom_int* bplane_cell_ghst,
                                 const flow_float* ro, flow_float* Q, flow_float* roQ)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nb * nSlice) return;
    const int s = (int)(t / nb); const geom_int ib = (geom_int)(t % nb);
    const geom_int ic = bplane_cell[ib], ig = bplane_cell_ghst[ib];
    const size_t off = (size_t)s * nCells_all;
    if (isInlet) {
        const int var = s / g_cmc.nEta, k = s % g_cmc.nEta; const double e = g_cmc.eta[k];
        const double q = (var < g_cmc.nSpecies) ? (e * g_cmc.YF[var] + (1.0 - e) * g_cmc.YO[var]) : (e * g_cmc.hF + (1.0 - e) * g_cmc.hO);
        Q[off + ig] = (flow_float)q; roQ[off + ig] = (flow_float)((double)ro[ig] * q);
        if (isNode) { Q[off + ic] = (flow_float)q; roQ[off + ic] = (flow_float)((double)ro[ic] * q); }
    } else {
        Q[off + ig] = Q[off + ic]; roQ[off + ig] = roQ[off + ic];
    }
}

// 保存量 roQ = ρ̄ Q を現在の ρ̄ で作り直す (輸送 point-implicit の直前)。前ステップの ρ̄ で作った roQ を更新後の ρ̄ で割ると
// 平均密度の変化分だけ Q が偽に動く (凍結検証で条件付き T が 1045 K を超えた指紋) ので、更新の前後で同じ ρ̄ を使う。
__global__ void cmc_q_resync_d(geom_int nCells_all, int nSlice, const flow_float* ro, const flow_float* Q, flow_float* roQ)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nCells_all * nSlice) return;
    const int s = (int)(t / nCells_all); const geom_int ic = (geom_int)(t % nCells_all);
    const size_t i = (size_t)s * nCells_all + ic;
    roQ[i] = ro[ic] * Q[i];
}

// 非保存補正: 収束していない平均流では Σ_f ṁ_f ≠ 0 なので、保存形の移流残差は一様な Q でも −Q Σṁ_f = Q·res_ro が残る。
// 流れ側の密度更新 (block-DPLUR) と点陰解のスカラー更新は同じ Δρ を与えないため、これが毎ステップ Q を (ρ+Δρ_s)/(ρ+Δρ_f) 倍
// 偽に動かす (凍結検証で条件付き T が ±40 K ずれた指紋)。res_Q -= Q·res_ro で「Q の輸送方程式」(非保存形) に揃える。
__global__ void cmc_q_nonconservative_fix_d(geom_int nCells, int nSlice, geom_int nCells_all, const flow_float* res_ro, const flow_float* Q, flow_float* resQ)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nCells * nSlice) return;
    const int s = (int)(t / nCells); const geom_int ic = (geom_int)(t % nCells);
    const size_t i = (size_t)s * nCells_all + ic;
    resQ[i] -= Q[i] * res_ro[ic];
}

// 原始量 Q = roQ/ρ (輸送後)
__global__ void cmc_q_primitive_d(geom_int nCells, int nSlice, geom_int nCells_all, const flow_float* ro, const flow_float* roQ, flow_float* Q)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nCells * nSlice) return;
    const int s = (int)(t / nCells); const geom_int ic = (geom_int)(t % nCells);
    const size_t i = (size_t)s * nCells_all + ic;
    Q[i] = roQ[i] / fmaxf(ro[ic], 1.0e-30f);
}

// β-PDF の離散重み Ω_k (Σ=1)。分散が微小ならデルタ (隣接 2 点へ線形配分)。
__device__ void cmc_pdf_weights(double xi, double var, double* Om)
{
    const int ne = g_cmc.nEta;
    xi = fmin(fmax(xi, 0.0), 1.0);
    const double vmax = xi * (1.0 - xi);
    for (int k = 0; k < ne; ++k) Om[k] = 0.0;
    if (var < 1.0e-4 * vmax + 1.0e-10 || vmax < 1.0e-12) {
        int k = 0; while (k < ne - 2 && g_cmc.eta[k + 1] <= xi) ++k;
        const double d = g_cmc.eta[k + 1] - g_cmc.eta[k];
        const double f = fmin(fmax((xi - g_cmc.eta[k]) / d, 0.0), 1.0);
        Om[k] = 1.0 - f; Om[k + 1] = f;
        return;
    }
    double gam = vmax / fmin(var, 0.999 * vmax) - 1.0;
    gam = fmax(gam, 0.05);
    const double a = xi * gam, b = (1.0 - xi) * gam;
    double lnOm[CMC_MAX_ETA]; double lmax = -1.0e300;
    for (int k = 1; k < ne - 1; ++k) {
        const double e = g_cmc.eta[k];
        lnOm[k] = log(g_cmc.w[k]) + (a - 1.0) * log(e) + (b - 1.0) * log(1.0 - e);
        if (lnOm[k] > lmax) lmax = lnOm[k];
    }
    // 端ビン: ∫_0^{η1/2} η^{a-1} dη = (η1/2)^a / a (特異性の解析処理), 同様に η=1 側
    lnOm[0]      = a * log(0.5 * g_cmc.eta[1]) - log(a);
    lnOm[ne - 1] = b * log(0.5 * (1.0 - g_cmc.eta[ne - 2])) - log(b);
    if (lnOm[0] > lmax) lmax = lnOm[0];
    if (lnOm[ne - 1] > lmax) lmax = lnOm[ne - 1];
    double sum = 0.0;
    for (int k = 0; k < ne; ++k) { Om[k] = exp(lnOm[k] - lmax); sum += Om[k]; }
    for (int k = 0; k < ne; ++k) Om[k] /= sum;
}

// 小行列の Gauss 消去 (部分ピボット): A x = b を解く (A は n×n 行優先, 破壊)
__device__ bool cmc_solve_dense(int n, double* A, double* b)
{
    for (int c = 0; c < n; ++c) {
        int p = c; double pm = fabs(A[c*n + c]);
        for (int r = c + 1; r < n; ++r) if (fabs(A[r*n + c]) > pm) { pm = fabs(A[r*n + c]); p = r; }
        if (pm < 1.0e-300) return false;
        if (p != c) { for (int j = 0; j < n; ++j) { double t = A[c*n+j]; A[c*n+j] = A[p*n+j]; A[p*n+j] = t; } double t = b[c]; b[c] = b[p]; b[p] = t; }
        for (int r = c + 1; r < n; ++r) {
            const double f = A[r*n + c] / A[c*n + c];
            if (f == 0.0) continue;
            for (int j = c; j < n; ++j) A[r*n + j] -= f * A[c*n + j];
            b[r] -= f * b[c];
        }
    }
    for (int c = n - 1; c >= 0; --c) {
        double s = b[c];
        for (int j = c + 1; j < n; ++j) s -= A[c*n + j] * b[j];
        b[c] = s / A[c*n + c];
    }
    return true;
}

// node 毎: η 拡散 (陰的, AMC) → 各 η 点で化学 (点陰解) → PDF 平均 (ω̄, Q̇̄, J̄, Ỹ_pdf)
__global__ void cmc_step_d(geom_int nCells, geom_int nCells_all, const SpeciesThermo* sp, const ReactionTable* rt,
                           const flow_float* ro, const flow_float* P, const flow_float* xi, const flow_float* xiVar, const flow_float* chi,
                           const flow_float* dt_local, flow_float* const* roY_dev,
                           flow_float* Q, flow_float* roQ,
                           flow_float* omegaBar, flow_float* qdotBar, flow_float* jacBar, flow_float* cmc_dY, flow_float* cmc_TQmax,
                           flow_float* ypdfOut, flow_float* hpdfOut, flow_float* tauOut, flow_float* cmc_TQst, int doChem,
                           flow_float* cmc_xiOm, flow_float* cmc_xipdf, flow_float* dpdfOut, flow_float* qrelOut)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    const double rho = (double)ro[ic], p = (double)P[ic];
    const double dtau = (double)dt_local[ic] * g_cmc.dtScale;

    double Om[CMC_MAX_ETA], chiEta[CMC_MAX_ETA];
    cmc_pdf_weights((double)xi[ic], (double)xiVar[ic], Om);
    // AMC: <χ|η> = χ̃ G(η) / Σ Ω_k G_k
    double gnorm = 0.0;
    for (int k = 0; k < ne; ++k) gnorm += Om[k] * g_cmc.G[k];
    const double chim = fmax((double)chi[ic], 0.0);
    for (int k = 0; k < ne; ++k) chiEta[k] = (gnorm > 1.0e-30) ? chim * g_cmc.G[k] / gnorm : chim * g_cmc.G[k];

    // ---- η 拡散 (陰的三重対角, 端 Dirichlet) を各変数に ----
    // doChem==0 (初期化時の PDF 積分埋め) では拡散もスキップ: 初期化時は dt_local 未設定で Δτ がゴミ値になり、
    // 三重対角解が内点 Q を壊す (run_0092 プローブ: 端点だけ混合線・内点 0/NaN)。
    double a[CMC_MAX_ETA], bb[CMC_MAX_ETA], cc[CMC_MAX_ETA], d[CMC_MAX_ETA];
    for (int var = 0; var <= ns && doChem; ++var) {
        for (int k = 1; k < ne - 1; ++k) {
            const double hm = g_cmc.eta[k] - g_cmc.eta[k-1], hp = g_cmc.eta[k+1] - g_cmc.eta[k], hc = 0.5 * (hm + hp);
            const double D = 0.5 * chiEta[k] * dtau;
            a[k] = -D / (hm * hc); cc[k] = -D / (hp * hc); bb[k] = 1.0 - a[k] - cc[k];
            d[k] = (double)Q[sl(var, k) * nCells_all + ic];
        }
        const double q0 = (double)Q[sl(var, 0) * nCells_all + ic], qN = (double)Q[sl(var, ne-1) * nCells_all + ic];
        d[1] -= a[1] * q0; d[ne-2] -= cc[ne-2] * qN;
        // Thomas
        for (int k = 2; k < ne - 1; ++k) { const double m = a[k] / bb[k-1]; bb[k] -= m * cc[k-1]; d[k] -= m * d[k-1]; }
        d[ne-2] /= bb[ne-2];
        for (int k = ne - 3; k >= 1; --k) d[k] = (d[k] - cc[k] * d[k+1]) / bb[k];
        for (int k = 1; k < ne - 1; ++k) Q[sl(var, k) * nCells_all + ic] = (flow_float)d[k];
    }

    // ---- 化学 (各 η 点, 点陰解) + PDF 平均 ----
    double obar[THERMO_MAX_SPECIES], ypdf[THERMO_MAX_SPECIES];
    double jbar[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES];
    for (int s = 0; s < ns; ++s) { obar[s] = 0.0; ypdf[s] = 0.0; }
    for (int i = 0; i < ns * ns; ++i) jbar[i] = 0.0;
    double qbar = 0.0, TQmax = 0.0, hpdf = 0.0, TQst = 0.0, qrel = 0.0;   // qrel: この step の条件付き反応熱の PDF 平均 [J/kg]
    double Y[THERMO_MAX_SPECIES], omega[THERMO_MAX_SPECIES], dOdT[THERMO_MAX_SPECIES], J[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES];
    double A[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES], rhs[THERMO_MAX_SPECIES];
    for (int k = 0; k < ne; ++k) {
        double ysum = 0.0;
        for (int s = 0; s < ns; ++s) { double y = (double)Q[sl(s, k) * nCells_all + ic]; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
        if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
        double h = (double)Q[sl(ns, k) * nCells_all + ic];
        double T = thermo_T_from_h(sp, ns, Y, h, fmin(fmax(300.0 + h / 1200.0, 200.0), 4000.0), 200.0, g_cmc.Tmax);
        if (!(T == T)) { TQmax = nan(""); }          // 非有限は伝播させて診断で見える化 (fmax は NaN を無視する)
        if (T > TQmax) TQmax = T;
        if (k == g_cmc.kSt) TQst = T;
        for (int s = 0; s < ns; ++s) ypdf[s] += Om[k] * Y[s];
        hpdf += Om[k] * h;
        const bool interior = (k > 0 && k < ne - 1);
        if (!doChem || !g_cmc.chemOn || !interior || Om[k] < g_cmc.pdfFloor || T < g_cmc.Tfreeze) continue;
        const double rhoEta = p / (thermo_R_mix(sp, ns, Y) * T);
        double Qd = 0.0;
        chem_source(sp, rt, rhoEta, fmin(T, g_cmc.Tmax), Y, omega, &Qd, 2, J, dOdT);
        // PDF 平均 (更新前の状態で評価; 1 step の遅れは擬似時間収束で消える)
        for (int s = 0; s < ns; ++s) obar[s] += Om[k] * omega[s];
        for (int i = 0; i < ns * ns; ++i) jbar[i] += Om[k] * J[i];
        qbar += Om[k] * Qd;
        // 点陰解: (I - Δτ J) δY = Δτ ω/ρ_η
        for (int s = 0; s < ns; ++s) { for (int q = 0; q < ns; ++q) A[s*ns + q] = -dtau * J[s*ns + q]; A[s*ns + s] += 1.0; rhs[s] = dtau * omega[s] / rhoEta; }
        if (cmc_solve_dense(ns, A, rhs)) {
            for (int s = 0; s < ns; ++s) qrel -= Om[k] * (sp[s].h_datum / sp[s].MW) * rhs[s];   // −Σ c_s δY_s [J/kg]
            ysum = 0.0;
            for (int s = 0; s < ns; ++s) { Y[s] = fmax(Y[s] + rhs[s], 0.0); ysum += Y[s]; }
            if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
            h += dtau * Qd / rhoEta;
            for (int s = 0; s < ns; ++s) Q[sl(s, k) * nCells_all + ic] = (flow_float)Y[s];
            Q[sl(ns, k) * nCells_all + ic] = (flow_float)h;
        }
    }
    // 保存量 roQ = ρ̄ Q (輸送の重み), 診断, PDF 平均ソース
    for (int s = 0; s <= ns; ++s) for (int k = 0; k < ne; ++k) { const size_t i = sl(s, k) * nCells_all + ic; roQ[i] = (flow_float)(rho * (double)Q[i]); }
    double dy = 0.0;
    for (int s = 0; s < ns; ++s) { const double ym = (double)roY_dev[s][ic] / fmax(rho, 1.0e-30); dy = fmax(dy, fabs(ypdf[s] - ym)); }
    cmc_dY[ic] = (flow_float)dy; cmc_TQmax[ic] = (flow_float)TQmax; cmc_TQst[ic] = (flow_float)TQst;
    { double xo = 0.0, bpdf = 0.0; for (int k = 0; k < ne; ++k) xo += Om[k] * g_cmc.eta[k];
      for (int s = 0; s < ns; ++s) bpdf += g_bilger.c[s] * ypdf[s];
      cmc_xiOm[ic] = (flow_float)xo; cmc_xipdf[ic] = (flow_float)fmin(fmax((bpdf - g_bilger.betaO) / (g_bilger.betaF - g_bilger.betaO), 0.0), 1.0);
      for (int s = 0; s < ns; ++s) dpdfOut[(size_t)s * nCells + ic] = (flow_float)(ypdf[s] - (xo * g_cmc.YF[s] + (1.0 - xo) * g_cmc.YO[s])); }
    for (int s = 0; s < ns; ++s) ypdfOut[(size_t)s * nCells + ic] = (flow_float)ypdf[s];
    hpdfOut[ic] = (flow_float)hpdf;
    tauOut[ic]  = (flow_float)fmax((double)dt_local[ic] * g_cmc.relax, 1.0e-9);   // 下限: 初期化直後 (dt_local 未設定) でも有限
    for (int s = 0; s < ns; ++s) omegaBar[(size_t)s * nCells + ic] = (flow_float)obar[s];
    qdotBar[ic] = (flow_float)qbar; qrelOut[ic] = (flow_float)qrel;
    for (int i = 0; i < ns * ns; ++i) jacBar[(size_t)ic * ns * ns + i] = (flow_float)jbar[i];
}

// =============================================================================
// 分割版 cmc_step: 旧 cmc_step_d は node 1 thread が 41 η 点の化学 (9×9 点陰解) を直列に回すため占有率が低く
// 620 ms/step を占めていた (CMC_TIMING 計測)。同じ演算を node×η の thread に展開する。評価順序・意味は旧 kernel と同じ
// (拡散 → PDF 平均は更新前の Q → 各 η の化学 → 診断)。couple 1/2 の ω̄・J̄ だけ double atomicAdd で集計する (順序非決定;
// 生産の couple 5 は使わない)。
// =============================================================================
// (A) node 毎: β-PDF 重み Ω_k と AMC <χ|η>_k
__global__ void cmc_stepA_weights_d(geom_int nCells, const flow_float* xi, const flow_float* xiVar, const flow_float* chi,
                                    double* OmOut, double* chiOut)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int ne = g_cmc.nEta;
    double Om[CMC_MAX_ETA];
    cmc_pdf_weights((double)xi[ic], (double)xiVar[ic], Om);
    double gnorm = 0.0;
    for (int k = 0; k < ne; ++k) gnorm += Om[k] * g_cmc.G[k];
    const double chim = fmax((double)chi[ic], 0.0);
    for (int k = 0; k < ne; ++k) {
        OmOut[(size_t)k * nCells + ic] = Om[k];
        chiOut[(size_t)k * nCells + ic] = (gnorm > 1.0e-30) ? chim * g_cmc.G[k] / gnorm : chim * g_cmc.G[k];
    }
}

// (B) (node, 変数) 毎: η 拡散 (陰的三重対角, 端 Dirichlet)。対角優位なので float の Thomas で十分 (fp64 版と 1e-6 で一致)。
__global__ void cmc_stepB_etadiff_d(geom_int nCells, geom_int nCells_all, const flow_float* dt_local, const double* chiEta, flow_float* Q)
{
    const size_t idx = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    if (idx >= (size_t)(ns + 1) * nCells) return;
    const int var = (int)(idx / nCells); const geom_int ic = (geom_int)(idx % nCells);
    const float dtau = (float)((double)dt_local[ic] * g_cmc.dtScale);
    float a[CMC_MAX_ETA], bb[CMC_MAX_ETA], cc[CMC_MAX_ETA], d[CMC_MAX_ETA];
    for (int k = 1; k < ne - 1; ++k) {
        const float hm = (float)(g_cmc.eta[k] - g_cmc.eta[k-1]), hp = (float)(g_cmc.eta[k+1] - g_cmc.eta[k]), hc = 0.5f * (hm + hp);
        const float D = 0.5f * (float)chiEta[(size_t)k * nCells + ic] * dtau;
        a[k] = -D / (hm * hc); cc[k] = -D / (hp * hc); bb[k] = 1.0f - a[k] - cc[k];
        d[k] = Q[sl(var, k) * nCells_all + ic];
    }
    const float q0 = Q[sl(var, 0) * nCells_all + ic], qN = Q[sl(var, ne-1) * nCells_all + ic];
    d[1] -= a[1] * q0; d[ne-2] -= cc[ne-2] * qN;
    for (int k = 2; k < ne - 1; ++k) { const float m = a[k] / bb[k-1]; bb[k] -= m * cc[k-1]; d[k] -= m * d[k-1]; }
    d[ne-2] /= bb[ne-2];
    for (int k = ne - 3; k >= 1; --k) d[k] = (d[k] - cc[k] * d[k+1]) / bb[k];
    for (int k = 1; k < ne - 1; ++k) Q[sl(var, k) * nCells_all + ic] = d[k];
}

// (C) node 毎: PDF 平均 Ỹ_pdf, h̃_pdf (化学更新前の Q で評価 = 旧 kernel と同じ)
__global__ void cmc_stepC_average_d(geom_int nCells, geom_int nCells_all, const double* Om, const flow_float* Q, double* ypdfD, double* hpdfD)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    double ypdf[THERMO_MAX_SPECIES]; for (int s = 0; s < ns; ++s) ypdf[s] = 0.0;
    double hpdf = 0.0, Y[THERMO_MAX_SPECIES];
    for (int k = 0; k < ne; ++k) {
        const double om = Om[(size_t)k * nCells + ic];
        double ysum = 0.0;
        for (int s = 0; s < ns; ++s) { double y = (double)Q[sl(s, k) * nCells_all + ic]; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
        if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
        for (int s = 0; s < ns; ++s) ypdf[s] += om * Y[s];
        hpdf += om * (double)Q[sl(ns, k) * nCells_all + ic];
    }
    for (int s = 0; s < ns; ++s) ypdfD[(size_t)s * nCells + ic] = ypdf[s];
    hpdfD[ic] = hpdf;
}

// CMC 専用の T(h,Y) Newton 反転: thermo_T_from_h と同じ式だが lnT を種ループの外で 1 回だけ取り (fp64 log が 9 分の 1)、
// 前 step の T(η) をウォームスタートにする (通常 1–2 反復で収束)。範囲外は thermo_cph_molar にフォールバック。
__device__ __forceinline__ double cmc_T_from_h(const SpeciesThermo* sp, int n, const double* Y, double h, double T_guess, double T_min, double T_max)
{
    double T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it = 0; it < 20; ++it) {
        const double Ti = 1.0 / T, Ti2 = Ti * Ti, lnT = log(T), T2 = T * T, T3 = T2 * T, T4 = T3 * T;
        double cp = 0.0, h_T = 0.0;
        for (int i = 0; i < n; ++i) {
            double cpi, hi;
            if (T >= sp[i].Tlo && T <= sp[i].Thi) {
                const double* a = thermo_pick_coeffs(sp[i], T);
                cpi = THERMO_RU * (a[0]*Ti2 + a[1]*Ti + a[2] + a[3]*T + a[4]*T2 + a[5]*T3 + a[6]*T4);
                hi  = THERMO_RU * T * (-a[0]*Ti2 + a[1]*lnT*Ti + a[2] + a[3]*T/2.0 + a[4]*T2/3.0 + a[5]*T3/4.0 + a[6]*T4/5.0 + a[7]*Ti);
            } else { thermo_cph_molar(sp[i], T, &cpi, &hi); }
            cp += Y[i] * (cpi / sp[i].MW); h_T += Y[i] * (hi / sp[i].MW);
        }
        const double cpf = (cp > 1.0e-3 ? cp : 1.0e-3);
        double dT = (h_T - h) / cpf;
        if (dT >  0.5*T) dT =  0.5*T;
        if (dT < -0.5*T) dT = -0.5*T;
        T -= dT;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*T) break;
    }
    return T;
}

// (D0) (node, η) 毎: 条件付き温度 T(η) と化学を評価するかのフラグ
__global__ void cmc_stepD0_temp_d(geom_int nCells, geom_int nCells_all, const SpeciesThermo* sp, const double* Om, const flow_float* Q,
                                  double* TEta, double* qrelEta, double* QdEta, int* flag, int doChem)
{
    const size_t idx = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    if (idx >= (size_t)ne * nCells) return;
    const int k = (int)(idx / nCells); const geom_int ic = (geom_int)(idx % nCells);
    double Y[THERMO_MAX_SPECIES]; double ysum = 0.0;
    for (int s = 0; s < ns; ++s) { double y = (double)Q[sl(s, k) * nCells_all + ic]; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
    if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
    const double h = (double)Q[sl(ns, k) * nCells_all + ic];
    const double Tprev = TEta[idx];   // 前 step の T(η) (初回は 0 → 粗い推定にフォールバック)
    const double Tg = (Tprev > 200.0 && Tprev < g_cmc.Tmax) ? Tprev : fmin(fmax(300.0 + h / 1200.0, 200.0), 4000.0);
    const double T = cmc_T_from_h(sp, ns, Y, h, Tg, 200.0, g_cmc.Tmax);
    TEta[idx] = T; qrelEta[idx] = 0.0; QdEta[idx] = 0.0;
    const bool interior = (k > 0 && k < ne - 1);
    flag[idx] = (doChem && g_cmc.chemOn && interior && Om[idx] >= g_cmc.pdfFloor && T >= g_cmc.Tfreeze) ? 1 : 0;
}

// (D1) 活性対リスト毎: 化学 (点陰解)。qrelEta / QdEta は Ω 重み付き。
__global__ void cmc_stepD1_chem_d(geom_int nCells, geom_int nCells_all, int nActive, const int* list, const SpeciesThermo* sp, const ReactionTable* rt,
                                  const flow_float* P, const flow_float* dt_local, const double* Om, flow_float* Q, const double* TEta,
                                  double* qrelEta, double* QdEta, double* omegaAcc, double* jacAcc, int accumSrc, unsigned long long* omHist)
{
    const int ia = blockDim.x * blockIdx.x + threadIdx.x;
    if (ia >= nActive) return;
    const size_t idx = (size_t)list[ia];
    const int ns = g_cmc.nSpecies;
    const int k = (int)(idx / nCells); const geom_int ic = (geom_int)(idx % nCells);
    const double om = Om[idx], T = TEta[idx];
    if (omHist) atomicAdd(&omHist[(om < 1.0e-5) ? 0 : (om < 1.0e-4) ? 1 : (om < 1.0e-3) ? 2 : 3], 1ULL);
    double Y[THERMO_MAX_SPECIES]; double ysum = 0.0;
    for (int s = 0; s < ns; ++s) { double y = (double)Q[sl(s, k) * nCells_all + ic]; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
    if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
    double h = (double)Q[sl(ns, k) * nCells_all + ic];
    const double p = (double)P[ic], dtau = (double)dt_local[ic] * g_cmc.dtScale;
    const double rhoEta = p / (thermo_R_mix(sp, ns, Y) * T);
    double omega[THERMO_MAX_SPECIES], dOdT[THERMO_MAX_SPECIES], J[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES];
    double A[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES], rhs[THERMO_MAX_SPECIES];
    double Qd = 0.0;
    chem_source(sp, rt, rhoEta, fmin(T, g_cmc.Tmax), Y, omega, &Qd, 2, J, dOdT);
    if (accumSrc) {
        for (int s = 0; s < ns; ++s) atomicAdd(&omegaAcc[(size_t)s * nCells + ic], om * omega[s]);
        for (int i = 0; i < ns * ns; ++i) atomicAdd(&jacAcc[(size_t)ic * ns * ns + i], om * J[i]);
    }
    QdEta[idx] = om * Qd;
    for (int s = 0; s < ns; ++s) { for (int q = 0; q < ns; ++q) A[s*ns + q] = -dtau * J[s*ns + q]; A[s*ns + s] += 1.0; rhs[s] = dtau * omega[s] / rhoEta; }
    if (cmc_solve_dense(ns, A, rhs)) {
        double qr = 0.0;
        for (int s = 0; s < ns; ++s) qr -= om * (sp[s].h_datum / sp[s].MW) * rhs[s];
        qrelEta[idx] = qr;
        ysum = 0.0;
        for (int s = 0; s < ns; ++s) { Y[s] = fmax(Y[s] + rhs[s], 0.0); ysum += Y[s]; }
        if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
        h += dtau * Qd / rhoEta;
        for (int s = 0; s < ns; ++s) Q[sl(s, k) * nCells_all + ic] = (flow_float)Y[s];
        Q[sl(ns, k) * nCells_all + ic] = (flow_float)h;
    }
}

// (D0f/D1f) fp32 版 (cmc.fp32=1, 既定): 速度定数・Jacobian・T 反転を float、点陰解だけ double。式は D0/D1 と同じ。
__global__ void cmc_stepD0f_temp_d(geom_int nCells, geom_int nCells_all, const double* Om, const flow_float* Q,
                                   double* TEta, double* qrelEta, double* QdEta, int* flag, int doChem)
{
    const size_t idx = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    if (idx >= (size_t)ne * nCells) return;
    const int k = (int)(idx / nCells); const geom_int ic = (geom_int)(idx % nCells);
    float Y[THERMO_MAX_SPECIES]; float ysum = 0.0f;
    for (int s = 0; s < ns; ++s) { float y = Q[sl(s, k) * nCells_all + ic]; if (y < 0.0f) y = 0.0f; Y[s] = y; ysum += y; }
    if (ysum > 0.0f) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
    const float h = Q[sl(ns, k) * nCells_all + ic];
    const float Tmax = (float)g_cmc.Tmax, Tprev = (float)TEta[idx];
    const float Tg = (Tprev > 200.0f && Tprev < Tmax) ? Tprev : fminf(fmaxf(300.0f + h / 1200.0f, 200.0f), 4000.0f);
    const float T = chemf_T_from_h(g_spf, ns, Y, h, Tg, 200.0f, Tmax);
    TEta[idx] = (double)T; qrelEta[idx] = 0.0; QdEta[idx] = 0.0;
    const bool interior = (k > 0 && k < ne - 1);
    flag[idx] = (doChem && g_cmc.chemOn && interior && Om[idx] >= g_cmc.pdfFloor && (double)T >= g_cmc.Tfreeze) ? 1 : 0;
}

__global__ void cmc_stepD1f_chem_d(geom_int nCells, geom_int nCells_all, int nActive, const int* list,
                                   const flow_float* P, const flow_float* dt_local, const double* Om, flow_float* Q, const double* TEta,
                                   double* qrelEta, double* QdEta, double* omegaAcc, double* jacAcc, int accumSrc)
{
    const int ia = blockDim.x * blockIdx.x + threadIdx.x;
    if (ia >= nActive) return;
    const size_t idx = (size_t)list[ia];
    const int ns = g_cmc.nSpecies;
    const int k = (int)(idx / nCells); const geom_int ic = (geom_int)(idx % nCells);
    const double om = Om[idx];
    const float T = (float)TEta[idx];
    float Y[THERMO_MAX_SPECIES]; float ysum = 0.0f;
    for (int s = 0; s < ns; ++s) { float y = Q[sl(s, k) * nCells_all + ic]; if (y < 0.0f) y = 0.0f; Y[s] = y; ysum += y; }
    if (ysum > 0.0f) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
    float h = Q[sl(ns, k) * nCells_all + ic];
    const double dtau = (double)dt_local[ic] * g_cmc.dtScale;
    const float rhoEta = P[ic] / (chemf_R_mix(g_spf, ns, Y) * T);
    float omega[THERMO_MAX_SPECIES], Jf[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES], Qd = 0.0f;
    chemf_source(g_spf, &g_rtf, rhoEta, fminf(T, (float)g_cmc.Tmax), Y, omega, &Qd, Jf);
    if (accumSrc) {
        for (int s = 0; s < ns; ++s) atomicAdd(&omegaAcc[(size_t)s * nCells + ic], om * (double)omega[s]);
        for (int i = 0; i < ns * ns; ++i) atomicAdd(&jacAcc[(size_t)ic * ns * ns + i], om * (double)Jf[i]);
    }
    QdEta[idx] = om * (double)Qd;
    double A[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES], rhs[THERMO_MAX_SPECIES];
    for (int s = 0; s < ns; ++s) { for (int q = 0; q < ns; ++q) A[s*ns + q] = -dtau * (double)Jf[s*ns + q]; A[s*ns + s] += 1.0; rhs[s] = dtau * (double)omega[s] / (double)rhoEta; }
    if (cmc_solve_dense(ns, A, rhs)) {
        double qr = 0.0;
        for (int s = 0; s < ns; ++s) qr -= om * ((double)g_spf[s].h_datum / (double)g_spf[s].MW) * rhs[s];
        qrelEta[idx] = qr;
        ysum = 0.0f;
        for (int s = 0; s < ns; ++s) { Y[s] = fmaxf(Y[s] + (float)rhs[s], 0.0f); ysum += Y[s]; }
        if (ysum > 0.0f) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
        h += (float)(dtau * (double)Qd / (double)rhoEta);
        for (int s = 0; s < ns; ++s) Q[sl(s, k) * nCells_all + ic] = Y[s];
        Q[sl(ns, k) * nCells_all + ic] = h;
    }
}

// 全スライス一括の point-implicit 更新 (scalarTimeIntegration_d の timeIntegration 11 と同式: coef_N=1, coef_M=0, coef_Res=1,
//   ρφ = ρφ_N + (res Δτ/V) / (1 + Δτ(src_jac + diag/V)), floor なし)。410 回の小起動を 1 回にする。
__global__ void cmc_q_timeint_d(geom_int nCells, int nSlice, geom_int nCells_all, const flow_float* dt_local, const geom_float* vol,
                                flow_float* roQ, const flow_float* resQ, const flow_float* sjQ, const flow_float* diagQ)
{
    const size_t idx = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= (size_t)nSlice * nCells) return;
    const int s = (int)(idx / nCells); const geom_int ic = (geom_int)(idx % nCells);
    const size_t i = (size_t)s * nCells_all + ic;
    const flow_float dt_l = dt_local[ic]; const geom_float v = vol[ic];
    const flow_float fac = static_cast<flow_float>(1.0) + dt_l * (sjQ[i] + diagQ[i] / v);
    roQ[i] = roQ[i] + (resQ[i] * dt_l / v) / fac;
}

// (E) node 毎: roQ 再同期, 診断, PDF 平均ソースの書き出し (旧 kernel の末尾と同じ)
__global__ void cmc_stepE_finalize_d(geom_int nCells, geom_int nCells_all, const flow_float* ro, const flow_float* xi, const flow_float* dt_local, flow_float* const* roY_dev,
                                     const flow_float* Q, flow_float* roQ, const double* Om, const double* TEta, const double* qrelEta, const double* QdEta,
                                     const double* ypdfD, const double* hpdfD, const double* omegaAcc, const double* jacAcc,
                                     flow_float* omegaBar, flow_float* qdotBar, flow_float* jacBar, flow_float* cmc_dY, flow_float* cmc_TQmax,
                                     flow_float* ypdfOut, flow_float* hpdfOut, flow_float* tauOut, flow_float* cmc_TQst,
                                     flow_float* cmc_xiOm, flow_float* cmc_xipdf, flow_float* dpdfOut, flow_float* qrelOut, flow_float* gateOut)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    const double rho = (double)ro[ic];
    double TQmax = 0.0, TQst = 0.0, qrel = 0.0, qbar = 0.0, xo = 0.0;
    for (int k = 0; k < ne; ++k) {
        const size_t i = (size_t)k * nCells + ic;
        const double T = TEta[i];
        if (!(T == T)) TQmax = nan("");
        if (T > TQmax) TQmax = T;
        if (k == g_cmc.kSt) TQst = T;
        qrel += qrelEta[i]; qbar += QdEta[i]; xo += Om[i] * g_cmc.eta[k];
    }
    for (int s = 0; s <= ns; ++s) for (int k = 0; k < ne; ++k) { const size_t i = sl(s, k) * nCells_all + ic; roQ[i] = (flow_float)(rho * (double)Q[i]); }
    // 混合分率の整合補正 (2026-09-12): 離散 β-PDF 重みの 1 次モーメント ξ_Ω = Σ Ω_k η_k は端ビン (質量を η=0/1 の節点に載せる) と台形則の
    //   誤差で輸送 ξ̃ から希薄側 −0.5〜−2 %・濃厚側 +0.1 % ずれる。緩和ソース ω_s = ρ(Ỹ_pdf−Y)/τ_c はこの偏りを 1/τ_c の速さで平均場に
    //   刷り込み、希薄側で燃料元素の恒常的なシンクになって ξ が崩落する (run_0110: 50 ms で軸 ξ 0.5→0.008; 定常反復の run_0099/0101 も同じ)。
    //   目標組成・エンタルピーを混合線に沿って (ξ̃ − ξ_Ω) だけ平行移動し、Bilger ξ(Ỹ_pdf) = ξ̃ を厳密に保つ (反応離脱 D_pdf は不変)。
    const double xit = fmin(fmax((double)xi[ic], 0.0), 1.0), dxi = xit - xo;
    double dy = 0.0, bpdf = 0.0;
    for (int s = 0; s < ns; ++s) {
        const double yp = ypdfD[(size_t)s * nCells + ic] + dxi * (g_cmc.YF[s] - g_cmc.YO[s]);
        const double ym = (double)roY_dev[s][ic] / fmax(rho, 1.0e-30); dy = fmax(dy, fabs(yp - ym));
        bpdf += g_bilger.c[s] * yp;
        dpdfOut[(size_t)s * nCells + ic] = (flow_float)(yp - (xit * g_cmc.YF[s] + (1.0 - xit) * g_cmc.YO[s]));
        ypdfOut[(size_t)s * nCells + ic] = (flow_float)yp;
        omegaBar[(size_t)s * nCells + ic] = (flow_float)omegaAcc[(size_t)s * nCells + ic];
    }
    cmc_dY[ic] = (flow_float)dy; cmc_TQmax[ic] = (flow_float)TQmax; cmc_TQst[ic] = (flow_float)TQst;
    cmc_xiOm[ic] = (flow_float)xo; cmc_xipdf[ic] = (flow_float)fmin(fmax((bpdf - g_bilger.betaO) / (g_bilger.betaF - g_bilger.betaO), 0.0), 1.0);
    const double hp = hpdfD[ic] + dxi * (g_cmc.hF - g_cmc.hO);
    hpdfOut[ic] = (flow_float)hp;
    { const double hline = xit * g_cmc.hF + (1.0 - xit) * g_cmc.hO; gateOut[ic] = (flow_float)fmin(fmax((hp - hline) / (1500.0 * g_cmc.dTgate), 0.0), 1.0); }
    tauOut[ic]  = (flow_float)fmax((double)dt_local[ic] * g_cmc.relax, 1.0e-9);
    qdotBar[ic] = (flow_float)qbar; qrelOut[ic] = (flow_float)qrel;
    for (int i = 0; i < ns * ns; ++i) jacBar[(size_t)ic * ns * ns + i] = (flow_float)jacAcc[(size_t)ic * ns * ns + i];
}

// couple 3: 平均場の化学種・温度を PDF 積分値で上書き (文献 RANS-CMC の標準: 平均スカラーは CMC から診断)。
//   roY_s = ρ̄ Ỹ_pdf,s (正規化), T̃ = T(h̃_pdf, Ỹ_pdf), roe = ρ̄ (e_int(T̃,Ỹ) + |u|²/2), e_int = h − R_mix T (sensible datum)。
__global__ void cmc_overwrite_mean_d(geom_int nCells, const SpeciesThermo* sp, const flow_float* ypdf, const flow_float* hpdf,
                                     const flow_float* ro, const flow_float* roUx, const flow_float* roUy, const flow_float* roUz,
                                     flow_float* const* roY_dev, flow_float* roe, double Tmax)
{
    // 注意: 流れ更新直後に呼ばれるので原始量 (T, Ux, ...) は更新前の値のまま。保存量 (ρ, ρu, ρe, ρY) だけから現在の状態を
    // 再構成する (原始量を混ぜると α=0 の往復でも roe が壊れて発散した: run_0084_cmc_ab_a0)。
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int ns = g_cmc.nSpecies;
    const double rho = (double)ro[ic], a = g_cmc.alpha;
    if (!(rho > 0.0)) return;
    const double ux = (double)roUx[ic] / rho, uy = (double)roUy[ic] / rho, uz = (double)roUz[ic] / rho;
    const double ke = 0.5 * (ux * ux + uy * uy + uz * uz);
    const double eint_cur = (double)roe[ic] / rho - ke;
    double Yc[THERMO_MAX_SPECIES], Y[THERMO_MAX_SPECIES]; double ysum = 0.0, ycs = 0.0;
    for (int s = 0; s < ns; ++s) { double y = (double)roY_dev[s][ic] / rho; if (y < 0.0) y = 0.0; Yc[s] = y; ycs += y; }
    if (!(ycs > 0.0)) return;
    for (int s = 0; s < ns; ++s) Yc[s] /= ycs;
    const double Tc = thermo_T_from_e(sp, ns, Yc, eint_cur, 1000.0, 200.0, Tmax);
    if (!(Tc > 200.0) || !(Tc < Tmax)) return;   // 反転失敗なら触らない
    for (int s = 0; s < ns; ++s) { double y = (double)ypdf[(size_t)s * nCells + ic]; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
    if (!(ysum > 0.0)) return;
    for (int s = 0; s < ns; ++s) Y[s] = Yc[s] + a * (Y[s] / ysum - Yc[s]);
    const double hc = thermo_h_mix(sp, ns, Yc, Tc);
    const double h  = hc + a * ((double)hpdf[ic] - hc);
    const double T  = thermo_T_from_h(sp, ns, Y, h, Tc, 200.0, Tmax);
    if (!(T > 200.0) || !(T < Tmax)) return;
    const double eint = thermo_h_mix(sp, ns, Y, T) - thermo_R_mix(sp, ns, Y) * T;
    for (int s = 0; s < ns; ++s) roY_dev[s][ic] = (flow_float)(rho * Y[s]);
    roe[ic] = (flow_float)(rho * (eint + ke));
}

// couple 4: 組成だけを PDF 積分状態へ α ブレンドし、その組成変化に伴う反応熱 −Σ c_s Δ(ρY_s) をエネルギーに加える。
//   エンタルピー自体はブレンドしない (管壁の熱伝導・Le≠1 で平均 T は混合線 T(ξ) と一致しないため、h のブレンドは
//   リップで数百 K の非物理な冷却を作って発散した: run_0087)。平均エネルギー式が伝導・境界を担い、CMC は反応進行だけを渡す。
__global__ void cmc_blend_species_d(geom_int nCells, const SpeciesThermo* sp, const flow_float* ypdf, const flow_float* ro,
                                    flow_float* const* roY_dev, flow_float* roe, flow_float* chemQdot, const flow_float* dt_local,
                                    const flow_float* dpdf, flow_float* dpdf_prev, const flow_float* qrel, flow_float* qdebt, flow_float* qcap, int viaResidual,
                                    int mode, const flow_float* hpdf, const flow_float* xiOm, const flow_float* roUx, const flow_float* roUy, const flow_float* roUz,
                                    const flow_float* P)
{
    // couple 4: 組成を PDF 積分状態へ α ブレンド (熱を伴わない組成の同期)。反応熱は PDF 積分した反応離脱量
    //   D_pdf,s = Ỹ_pdf,s − line_s(ξ_Ω) の**時間変化**にだけ対応させる: ΔQ = −Σ c_s ρ (D_pdf − D_pdf,prev)。
    //   平均場の数値的な混合線からのずれ (種ごとのリミッタ差など) を補正する分に熱を計上すると、輸送が毎ステップずれを再生するので
    //   持続的な偽発熱になる (run_0093: T_Q 1392 K なのに平均 T 3477 K)。
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int ns = g_cmc.nSpecies;
    const double rho = (double)ro[ic], a = g_cmc.alpha;
    if (!(rho > 0.0)) return;
    double Yc[THERMO_MAX_SPECIES], Yp[THERMO_MAX_SPECIES]; double ycs = 0.0, yps = 0.0;
    for (int s = 0; s < ns; ++s) { double y = (double)roY_dev[s][ic] / rho; if (y < 0.0) y = 0.0; Yc[s] = y; ycs += y; }
    for (int s = 0; s < ns; ++s) { double y = (double)ypdf[(size_t)s * nCells + ic]; if (y < 0.0) y = 0.0; Yp[s] = y; yps += y; }
    if (!(ycs > 0.0) || !(yps > 0.0)) return;
    // 反応熱: 条件付き点陰解で実際に進んだ反応の PDF 平均 (qrel [J/kg], cmc_step_d) を平均密度で体積当たりに。
    // (D_pdf の時間差分は PDF 重みの移動でも ± に動き偽の吸発熱を作った: run_0091 v3 で −2e9 W/m3。)
    for (int s = 0; s < ns; ++s) {
        const double yn = Yc[s] / ycs + a * (Yp[s] / yps - Yc[s] / ycs);
        roY_dev[s][ic] = (flow_float)(rho * yn);
        const size_t i = (size_t)s * nCells + ic; dpdf_prev[i] = dpdf[i];
    }
    // 1 step の発熱を ΔT ≤ dTmax に制限し、残りは持ち越す (総発熱は保存)。擬似時間では条件付き化学が数百 step で燃え切り、
    // 流れが音響的に追従できない速さで熱が入って圧力が暴れる (run_0094: P 186 kPa → NaN) のを防ぐ。
    const double cvEst = 1000.0;   // [J/kg/K] 目安 (上限の判定にだけ使う)
    double qsrc = (double)qrel[ic];   // couple 4/5: 条件付き点陰解の反応熱の PDF 平均 [J/kg]
    if (mode == 6) {
        // couple 6: 平均 h̃ を PDF 積分 h̃_pdf へ α ブレンド (文献 RANS-CMC の「平均スカラーは Q から診断」に相当)。ただし燃焼領域だけに
        //   ゲート g = clip((h_pdf − h_line(ξ_Ω)) / (1500·dTgate), 0, 1) を掛け、未燃領域 (リップの壁伝熱・Le≠1 で h̃ が混合線から
        //   ずれる場所) は触らない。couple 5 の欠陥 = 熱は「条件付き化学が進んだ node の PDF 重み」でしか渡らず、上流の η_st で燃えて
        //   下流へ運ばれた Q の熱が平均場に届かない (run_0099: PDF 診断 T は実験と一致、平均 T は −120〜−220 K)。
        const double hp = (double)hpdf[ic], xo = fmin(fmax((double)xiOm[ic], 0.0), 1.0);
        const double hline = xo * g_cmc.hF + (1.0 - xo) * g_cmc.hO;
        const double g = fmin(fmax((hp - hline) / (1500.0 * g_cmc.dTgate), 0.0), 1.0);
        const double ux = (double)roUx[ic] / rho, uy = (double)roUy[ic] / rho, uz = (double)roUz[ic] / rho;
        const double hmean = (double)roe[ic] / rho - 0.5 * (ux * ux + uy * uy + uz * uz) + (double)P[ic] / rho;   // 顕エンタルピー [J/kg]
        qsrc = a * g * (hp - hmean);
    }
    double dQtot = rho * qsrc + (double)qdebt[ic];   // [J/m3]
    const double cap = rho * cvEst * g_cmc.dTmax;
    double dQ = dQtot;
    if (dQ >  cap) dQ =  cap;
    if (dQ < -cap) dQ = -cap;
    qdebt[ic] = (mode == 6) ? (flow_float)0.0 : (flow_float)(dQtot - dQ);   // couple 6 は状態緩和 (毎 step 目標から再計算) なので持ち越し不要。持ち越すと目標到達後も過冷却/過熱が続く (run_0102: ξ=0 の coflow が 770 K)
    if (viaResidual) { qcap[ic] = (flow_float)dQ; }   // couple 5: 次の残差組み立てで Q̇ = qcap/Δτ として res_roe へ (陰解法が圧力応答を処理)
    else            { roe[ic] += (flow_float)dQ; }   // couple 4: 直接加算
    const double dt = fmax((double)dt_local[ic], 1.0e-12);
    chemQdot[ic] = (flow_float)(dQ / dt);   // 診断: 等価な反応熱 [W/m3]
}

ScalarTransportDesc sliceDesc(int s, const solverConfig& cfg)
{
    const size_t off = (size_t)s * g_nCellsAll;
    ScalarTransportDesc d{};
    d.phi = g_Q + off; d.dphidx = g_gdum; d.dphidy = g_gdum; d.dphidz = g_gdum;
    d.rho_phi = g_roQ + off; d.rho_phi_N = g_roQ + off; d.rho_phi_M = g_roQ + off;
    d.res_rho_phi = g_resQ + off; d.res_rho_phi_m = g_resmQ + off;
    d.src_jac = g_sjQ + off; d.transport_diag = g_diagQ + off;
    d.floor = static_cast<flow_float>(-1.0e30);   // Q は負のエンタルピーもあり得る (フロアなし)
    d.sigma = static_cast<flow_float>(1.0 / cfg.Sc_t);
    d.diffusion = 1;
    return d;
}

} // namespace

// 区間計測: 環境変数 CMC_TIMING=<steps> で、その step 数毎に各区間の GPU 時間 (ms/step 平均) を stdout に出す
struct CmcTimer {
    int every = 0; int count = 0; cudaEvent_t ev[8]; double acc[7] = {0,0,0,0,0,0,0}; bool init = false;
    void begin() {
        if (!init) { const char* e = std::getenv("CMC_TIMING"); every = e ? std::atoi(e) : 0; if (every > 0) for (auto& x : ev) cudaEventCreate(&x); init = true; }
        if (every > 0) cudaEventRecord(ev[0]);
    }
    void mark(int i) { if (every > 0) cudaEventRecord(ev[i]); }
    void end(int n, const char* const* names) {
        if (every <= 0) return;
        cudaEventSynchronize(ev[n]);
        for (int i = 0; i < n; ++i) { float ms = 0; cudaEventElapsedTime(&ms, ev[i], ev[i + 1]); acc[i] += ms; }
        if (++count % every == 0) {
            std::printf("[cmc timing ms/step]");
            for (int i = 0; i < n; ++i) { std::printf(" %s=%.1f", names[i], acc[i] / every); acc[i] = 0; }
            std::printf("\n");
        }
    }
};
static CmcTimer g_timerR, g_timerU, g_timerS;
static int* g_activeFlag = nullptr; static int* g_activeList = nullptr; static int g_nActive = 0;   // 化学を評価する (node,η) 対の圧縮リスト (ワープ発散回避)
static unsigned long long* g_omHist = nullptr;   // 計測用: 活性対の Ω 分布 [<1e-5, <1e-4, <1e-3, ≥1e-3]

// デバッグプローブ: 環境変数 CMC_PROBE_NODE=<node index> [CMC_PROBE_EVERY=<steps>] で、その node の Q(η) を stdout に出す
static void cmcProbe(const char* tag)
{
    static int node = -2, every = 10, count = 0;
    if (node == -2) { const char* e = std::getenv("CMC_PROBE_NODE"); node = e ? std::atoi(e) : -1; const char* v = std::getenv("CMC_PROBE_EVERY"); if (v) every = std::atoi(v); }
    if (node < 0) return;
    if ((count++ % every) != 0) return;
    const int ns = g_cmc_host.nSpecies, ne = g_cmc_host.nEta;
    std::vector<flow_float> q((size_t)g_nSlice);
    for (int sl = 0; sl < g_nSlice; ++sl) gpuErrchk( cudaMemcpy(&q[sl], g_Q + (size_t)sl * g_nCellsAll + node, sizeof(flow_float), cudaMemcpyDeviceToHost) );
    std::printf("[cmc probe %s node %d] k eta Y_H2 Y_O2 Y_H2O Y_N2 h\n", tag, node);
    for (int k = 0; k < ne; k += 5) {
        std::printf("   %2d %.4f %.5f %.5f %.5f %.5f %.1f\n", k, g_cmc_host.eta[k], q[(size_t)(0*ne+k)], q[(size_t)(1*ne+k)], q[(size_t)(5*ne+k)], q[(size_t)(8*ne+k)], q[(size_t)(ns*ne+k)]);
    }
}

#include <fstream>
// Q(η) の永続化: 生バイナリ (flow_float × nSlice × nCells_all) + ヘッダ。restart で条件付き場を引き継ぐ (無ければ混合線から)。
void cmcQWrite(const std::string& path)
{
    if (!g_q_ready) return;
    std::vector<flow_float> h((size_t)g_nSlice * g_nCellsAll);
    gpuErrchk( cudaMemcpy(h.data(), g_Q, h.size() * sizeof(flow_float), cudaMemcpyDeviceToHost) );
    std::ofstream f(path, std::ios::binary);
    const int hdr[3] = { g_cmc_host.nEta, g_cmc_host.nSpecies, (int)g_nCellsAll };
    f.write((const char*)hdr, sizeof(hdr)); f.write((const char*)h.data(), h.size() * sizeof(flow_float));
}
static bool cmcQRead(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) { std::printf("cmcQRead: %s not found -> mixing-line init\n", path.c_str()); return false; }
    int hdr[3]; f.read((char*)hdr, sizeof(hdr));
    if (hdr[0] != g_cmc_host.nEta || hdr[1] != g_cmc_host.nSpecies || hdr[2] != (int)g_nCellsAll) {
        std::printf("cmcQRead: header mismatch (nEta %d/%d, ns %d/%d, n %d/%d) -> mixing-line init\n", hdr[0], g_cmc_host.nEta, hdr[1], g_cmc_host.nSpecies, hdr[2], (int)g_nCellsAll); return false;
    }
    std::vector<flow_float> h((size_t)g_nSlice * g_nCellsAll);
    f.read((char*)h.data(), h.size() * sizeof(flow_float));
    gpuErrchk( cudaMemcpy(g_Q, h.data(), h.size() * sizeof(flow_float), cudaMemcpyHostToDevice) );
    std::printf("cmcQRead: restored Q(eta) from %s\n", path.c_str());
    return true;
}

bool cmc_coupling_active() { return g_q_ready && g_cmc_host.couple != 0; }
int  cmc_coupling_mode()   { return g_q_ready ? g_cmc_host.couple : 0; }
const flow_float* cmc_ypdf_device_ptr() { return g_ypdf; }
const flow_float* cmc_hpdf_device_ptr() { return g_hpdf; }
const flow_float* cmc_tau_device_ptr()  { return g_tau; }
const flow_float* cmc_qcap_device_ptr() { return g_qcap; }
const flow_float* cmc_gate_device_ptr() { return g_gate; }
const flow_float* cmc_omega_device_ptr() { return g_omega; }
const flow_float* cmc_qdot_device_ptr()  { return g_qdot; }
const flow_float* cmc_jac_device_ptr()   { return g_jac; }

// cmc_step の起動: 既定は分割版 (A–E)。環境変数 CMC_STEP_LEGACY=1 で旧 1-kernel 版 (A/B 比較用)。
static void cmcStepLaunch(cudaConfig& cuda_cfg, variables& var, int doChem)
{
    static int legacy = -1;
    if (legacy < 0) { const char* e = std::getenv("CMC_STEP_LEGACY"); legacy = (e && std::atoi(e) != 0) ? 1 : 0; }
    if (legacy) {
        cmc_step_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            g_nCells, g_nCellsAll, thermo_species_device_ptr(), chemistry_table_device(),
            var.c_d.at("ro"), var.c_d.at("P"), var.c_d.at("xi"), var.c_d.at("xiVar"), var.c_d.at("chi"),
            var.c_d.at("dt_local"), species_roY_device_ptr(),
            g_Q, g_roQ, g_omega, g_qdot, g_jac, var.c_d.at("cmc_dY"), var.c_d.at("cmc_TQmax"),
            g_ypdf, g_hpdf, g_tau, var.c_d.at("cmc_TQst"), doChem, var.c_d.at("cmc_xiOm"), var.c_d.at("cmc_xipdf"), g_dpdf, g_qrel);
        return;
    }
    const int ns = g_cmc_host.nSpecies, ne = g_cmc_host.nEta, blk = 256;
    const int accumSrc = (g_cmc_host.couple == 1 || g_cmc_host.couple == 2) ? 1 : 0;
    auto grid = [&](size_t n) { return (unsigned)((n + blk - 1) / blk); };
    g_timerS.begin();
    cmc_stepA_weights_d<<<grid(g_nCells), blk>>>(g_nCells, var.c_d.at("xi"), var.c_d.at("xiVar"), var.c_d.at("chi"), g_OmD, g_chiEtaD);
    g_timerS.mark(1);
    if (doChem) cmc_stepB_etadiff_d<<<grid((size_t)(ns + 1) * g_nCells), blk>>>(g_nCells, g_nCellsAll, var.c_d.at("dt_local"), g_chiEtaD, g_Q);
    g_timerS.mark(2);
    cmc_stepC_average_d<<<grid(g_nCells), blk>>>(g_nCells, g_nCellsAll, g_OmD, g_Q, g_ypdfD, g_hpdfD);
    g_timerS.mark(3);
    if (accumSrc) { gpuErrchk( cudaMemsetAsync(g_omegaAccD, 0, (size_t)ns * g_nCells * sizeof(double)) ); gpuErrchk( cudaMemsetAsync(g_jacAccD, 0, (size_t)ns * ns * g_nCells * sizeof(double)) ); }
    const size_t nPair = (size_t)ne * g_nCells;
    if (g_cmcFp32) cmc_stepD0f_temp_d<<<grid(nPair), blk>>>(g_nCells, g_nCellsAll, g_OmD, g_Q, g_TEtaD, g_qrelEtaD, g_QdEtaD, g_activeFlag, doChem);
    else cmc_stepD0_temp_d<<<grid(nPair), blk>>>(g_nCells, g_nCellsAll, thermo_species_device_ptr(), g_OmD, g_Q, g_TEtaD, g_qrelEtaD, g_QdEtaD, g_activeFlag, doChem);
    {   // 活性対の圧縮リスト (thrust copy_if; 順序は idx 昇順で決定的)
        thrust::device_ptr<int> flag(g_activeFlag), list(g_activeList);
        auto end = thrust::copy_if(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>((int)nPair), flag, list, CmcFlagNonZero());
        g_nActive = (int)(end - list);
    }
    if (g_timerS.every > 0 && !g_omHist) { gpuErrchk( cudaMalloc((void**)&g_omHist, 4 * sizeof(unsigned long long)) ); gpuErrchk( cudaMemset(g_omHist, 0, 4 * sizeof(unsigned long long)) ); }
    if (g_nActive > 0 && g_cmcFp32) cmc_stepD1f_chem_d<<<grid((size_t)g_nActive), blk>>>(g_nCells, g_nCellsAll, g_nActive, g_activeList,
        var.c_d.at("P"), var.c_d.at("dt_local"), g_OmD, g_Q, g_TEtaD, g_qrelEtaD, g_QdEtaD, g_omegaAccD, g_jacAccD, accumSrc);
    else if (g_nActive > 0) cmc_stepD1_chem_d<<<grid((size_t)g_nActive), blk>>>(g_nCells, g_nCellsAll, g_nActive, g_activeList, thermo_species_device_ptr(), chemistry_table_device(),
        var.c_d.at("P"), var.c_d.at("dt_local"), g_OmD, g_Q, g_TEtaD, g_qrelEtaD, g_QdEtaD, g_omegaAccD, g_jacAccD, accumSrc, g_omHist);
    g_timerS.mark(4);
    cmc_stepE_finalize_d<<<grid(g_nCells), blk>>>(g_nCells, g_nCellsAll, var.c_d.at("ro"), var.c_d.at("xi"), var.c_d.at("dt_local"), species_roY_device_ptr(),
        g_Q, g_roQ, g_OmD, g_TEtaD, g_qrelEtaD, g_QdEtaD, g_ypdfD, g_hpdfD, g_omegaAccD, g_jacAccD,
        g_omega, g_qdot, g_jac, var.c_d.at("cmc_dY"), var.c_d.at("cmc_TQmax"),
        g_ypdf, g_hpdf, g_tau, var.c_d.at("cmc_TQst"), var.c_d.at("cmc_xiOm"), var.c_d.at("cmc_xipdf"), g_dpdf, g_qrel, g_gate);
    g_timerS.mark(5);
    { static const char* nm[] = {"A_weights", "B_etadiff", "C_average", "D_chem", "E_final"}; g_timerS.end(5, nm);
      if (g_timerS.every > 0 && (g_timerS.count % g_timerS.every) == 0) { unsigned long long hh[4] = {0,0,0,0}; gpuErrchk( cudaMemcpy(hh, g_omHist, sizeof(hh), cudaMemcpyDeviceToHost) ); gpuErrchk( cudaMemset(g_omHist, 0, sizeof(hh)) );
        const double tot = (double)(hh[0] + hh[1] + hh[2] + hh[3]) / g_timerS.every;
        std::printf("[cmc active (node,eta) pairs] %d of %zu (%.1f%%); Omega hist per step: <1e-5 %.0f, <1e-4 %.0f, <1e-3 %.0f, >=1e-3 %.0f\n", g_nActive, (size_t)ne * g_nCells,
                    100.0 * g_nActive / (double)((size_t)ne * g_nCells), hh[0] / (double)g_timerS.every, hh[1] / (double)g_timerS.every, hh[2] / (double)g_timerS.every, hh[3] / (double)g_timerS.every); (void)tot; } }
}

void cmcQInit_d(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    g_q_ready = false;
    if (cfg.chemCmc == 0 || !g_ready) return;
    const int ns = var.nSpeciesRegistered, ne = cfg.chemCmcNEta;
    if (ne > CMC_MAX_ETA) throw std::runtime_error("[cmc] nEta > CMC_MAX_ETA");
    CmcParams& c = g_cmc_host;
    std::memset(&c, 0, sizeof(CmcParams));
    c.nEta = ne; c.nSpecies = ns; c.chemOn = cfg.chemCmcChem; c.couple = cfg.chemCmcCouple;
    c.pdfFloor = cfg.chemCmcPdfFloor; c.Tmax = cfg.chemTmaxReaction; c.Tfreeze = cfg.chemFreezeBelowT; c.dtScale = cfg.chemCmcDtScale;
    c.relax = cfg.chemCmcRelax; c.alpha = cfg.chemCmcAlpha; c.dTmax = cfg.chemCmcDTmax; c.dTgate = cfg.chemCmcDTgate;
    for (int k = 0; k < ne; ++k) { const double s = (double)k / (double)(ne - 1); c.eta[k] = std::pow(s, cfg.chemCmcEtaPow); }
    c.eta[0] = 0.0; c.eta[ne - 1] = 1.0;
    { int kb = 0; double db = 1.0; for (int k = 0; k < ne; ++k) { const double dd = std::fabs(c.eta[k] - cfg.chemCmcXiSt); if (dd < db) { db = dd; kb = k; } } c.kSt = kb; }
    for (int k = 0; k < ne; ++k) {
        const double lo = (k == 0) ? 0.0 : 0.5 * (c.eta[k] + c.eta[k-1]);
        const double hi = (k == ne - 1) ? 1.0 : 0.5 * (c.eta[k] + c.eta[k+1]);
        c.w[k] = hi - lo;
        c.G[k] = 0.0;   // 端は 0、内点は下で AMC 形状
    }
    // AMC 形状 G(η) = exp(-2 [erfinv(2η-1)]^2): erfinv は host に無いので Newton で求める
    for (int k = 1; k < ne - 1; ++k) {
        const double y = 2.0 * c.eta[k] - 1.0; double x = 0.0;
        for (int it = 0; it < 60; ++it) { const double f = std::erf(x) - y; const double df = 2.0 / std::sqrt(M_PI) * std::exp(-x * x); x -= f / df; }
        c.G[k] = std::exp(-2.0 * x * x);
    }
    // 燃料/酸化剤流の質量分率 (fuelX/oxidizerX → Y)
    const ReactionTable* rt = chemistry_table_host();
    std::vector<double> W(ns, 0.0);
    for (int s = 0; s < ns; ++s) for (int e = 0; e < rt->nElem; ++e) W[s] += rt->atoms[s][e] * rt->elemW[e];
    auto toY = [&](const std::map<std::string, double>& X, double* Yout) {
        double mm = 0.0;
        for (const auto& kv : X) { int s = -1; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i; mm += kv.second * W[s]; }
        for (int s = 0; s < ns; ++s) Yout[s] = 0.0;
        for (const auto& kv : X) { int s = 0; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i; Yout[s] = kv.second * W[s] / mm; }
    };
    toY(cfg.chemMixfracFuelX, c.YF); toY(cfg.chemMixfracOxidX, c.YO);
    gpuErrchk( cudaMemcpyToSymbol(g_cmc, &c, sizeof(CmcParams)) );
    // 流れのエンタルピー (device thermo)
    double* hdev = nullptr; double hh[2];
    gpuErrchk( cudaMalloc((void**)&hdev, 2 * sizeof(double)) );
    cmc_stream_enthalpy_d<<<1, 32>>>(thermo_species_device_ptr(), ns, cfg.chemCmcTfuel, cfg.chemCmcTox, hdev);
    gpuErrchk( cudaMemcpy(hh, hdev, 2 * sizeof(double), cudaMemcpyDeviceToHost) ); cudaFree(hdev);
    c.hF = hh[0]; c.hO = hh[1];
    gpuErrchk( cudaMemcpyToSymbol(g_cmc, &c, sizeof(CmcParams)) );
    g_cmcFp32 = cfg.chemCmcFp32;
    if (g_cmcFp32) {   // fp32 係数表を __constant__ へ
        static SpeciesThermoF spf[THERMO_MAX_SPECIES]; static ReactionTableF rtf;
        chemf_build_tables(thermo_species_host(), ns, rt, spf, &rtf);
        gpuErrchk( cudaMemcpyToSymbol(g_spf, spf, sizeof(spf)) );
        gpuErrchk( cudaMemcpyToSymbol(g_rtf, &rtf, sizeof(ReactionTableF)) );
    }

    g_nSlice = (ns + 1) * ne; g_nCellsAll = msh.nCells_all; g_nCells = msh.nCells;
    const size_t bytes = (size_t)g_nSlice * g_nCellsAll * sizeof(flow_float);
    for (flow_float** pp : {&g_Q, &g_roQ, &g_resQ, &g_diagQ, &g_sjQ, &g_resmQ}) { if (*pp) cudaFree(*pp); gpuErrchk( cudaMalloc((void**)pp, bytes) ); gpuErrchk( cudaMemset(*pp, 0, bytes) ); }
    if (g_gdum) cudaFree(g_gdum); gpuErrchk( cudaMalloc((void**)&g_gdum, (size_t)g_nCellsAll * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_gdum, 0, (size_t)g_nCellsAll * sizeof(flow_float)) );
    if (g_omega) cudaFree(g_omega); gpuErrchk( cudaMalloc((void**)&g_omega, (size_t)ns * msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_omega, 0, (size_t)ns * msh.nCells * sizeof(flow_float)) );
    if (g_qdot) cudaFree(g_qdot); gpuErrchk( cudaMalloc((void**)&g_qdot, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_qdot, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_ypdf) cudaFree(g_ypdf); gpuErrchk( cudaMalloc((void**)&g_ypdf, (size_t)ns * msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_ypdf, 0, (size_t)ns * msh.nCells * sizeof(flow_float)) );
    if (g_hpdf) cudaFree(g_hpdf); gpuErrchk( cudaMalloc((void**)&g_hpdf, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_hpdf, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_dpdf) cudaFree(g_dpdf); gpuErrchk( cudaMalloc((void**)&g_dpdf, (size_t)ns * msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_dpdf, 0, (size_t)ns * msh.nCells * sizeof(flow_float)) );
    if (g_dpdf_prev) cudaFree(g_dpdf_prev); gpuErrchk( cudaMalloc((void**)&g_dpdf_prev, (size_t)ns * msh.nCells * sizeof(flow_float)) );
    if (g_qrel) cudaFree(g_qrel); gpuErrchk( cudaMalloc((void**)&g_qrel, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_qrel, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_qdebt) cudaFree(g_qdebt); gpuErrchk( cudaMalloc((void**)&g_qdebt, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_qdebt, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_qcap) cudaFree(g_qcap); gpuErrchk( cudaMalloc((void**)&g_qcap, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_qcap, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_gate) cudaFree(g_gate); gpuErrchk( cudaMalloc((void**)&g_gate, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_gate, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_tau)  cudaFree(g_tau);  gpuErrchk( cudaMalloc((void**)&g_tau,  (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_tau,  0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_jac) cudaFree(g_jac); gpuErrchk( cudaMalloc((void**)&g_jac, (size_t)msh.nCells * ns * ns * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_jac, 0, (size_t)msh.nCells * ns * ns * sizeof(flow_float)) );
    { auto allocD = [&](double** pp, size_t n) { if (*pp) cudaFree(*pp); gpuErrchk( cudaMalloc((void**)pp, n * sizeof(double)) ); gpuErrchk( cudaMemset(*pp, 0, n * sizeof(double)) ); };
      const size_t nk = (size_t)ne * msh.nCells, nsN = (size_t)ns * msh.nCells;
      allocD(&g_OmD, nk); allocD(&g_chiEtaD, nk); allocD(&g_TEtaD, nk); allocD(&g_qrelEtaD, nk); allocD(&g_QdEtaD, nk);
      allocD(&g_ypdfD, nsN); allocD(&g_hpdfD, (size_t)msh.nCells); allocD(&g_omegaAccD, nsN); allocD(&g_jacAccD, (size_t)msh.nCells * ns * ns); }
    for (int** pp : {&g_activeFlag, &g_activeList}) { if (*pp) cudaFree(*pp); gpuErrchk( cudaMalloc((void**)pp, (size_t)ne * msh.nCells * sizeof(int)) ); }

    cmc_q_init_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d.at("ro"), g_Q, g_roQ);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    cmcProbe("after-init-kernel");
    if (!cfg.chemCmcRestartQ.empty()) {
        if (cmcQRead(cfg.chemCmcRestartQ)) {   // roQ を復元した Q から作り直す
            const size_t n = (size_t)g_nCellsAll * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
            cmc_q_resync_d<<<grid, blk>>>(g_nCellsAll, g_nSlice, var.c_d.at("ro"), g_Q, g_roQ);
            gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        }
    }
    // PDF 積分値 (Ỹ_pdf, h̃_pdf, τ) を初期場で一度埋める (化学は走らせない)。これを怠ると初回の残差組み立てで
    // couple 2 の緩和ソースが Ỹ_pdf=0 に向かって全種を剥ぎ取り step 2 で NaN になる (run_0083 初版)。
    cmcStepLaunch(cuda_cfg, var, 0);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    gpuErrchk( cudaMemcpy(g_dpdf_prev, g_dpdf, (size_t)ns * msh.nCells * sizeof(flow_float), cudaMemcpyDeviceToDevice) );   // 初期の反応離脱 (=0 or 復元場の値) を基準に
    g_q_ready = true; g_cmcStep = 0;
    cmcProbe("init");
    std::printf("cmcQInit_d: nEta=%d nSlice=%d (%.1f MB x6), hF=%.4g hO=%.4g J/kg, couple=%d chem=%d\n",
                ne, g_nSlice, bytes / 1.0e6, c.hF, c.hO, c.couple, c.chemOn);
}

void cmcQTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_q_ready) return;
    const size_t bytes = (size_t)g_nSlice * g_nCellsAll * sizeof(flow_float);
    g_timerR.begin();
    CHECK_CUDA_ERROR(cudaMemset(g_resQ, 0, bytes)); CHECK_CUDA_ERROR(cudaMemset(g_diagQ, 0, bytes)); CHECK_CUDA_ERROR(cudaMemset(g_sjQ, 0, bytes));
    g_timerR.mark(1);
    for (int s = 0; s < g_nSlice; ++s) scalarTransportResidual_d(cfg, cuda_cfg, msh, var, sliceDesc(s, cfg));
    g_timerR.mark(2);
    cmcProbe("after-residual");
    { const size_t n = (size_t)g_nCells * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
      cmc_q_nonconservative_fix_d<<<grid, blk>>>(g_nCells, g_nSlice, g_nCellsAll, var.c_d.at("res_ro"), g_Q, g_resQ); }
    g_timerR.mark(3);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    { static const char* nm[] = {"memset", "residual", "nc_fix"}; g_timerR.end(3, nm); }
}

void applyCmcQBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_q_ready) return;
    const int isNode = (cfg.discretization == "node") ? 1 : 0;
    for (auto& bc : msh.bconds) {
        if (bc.iPlanes.empty()) continue;
        const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
        const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
        const size_t n = (size_t)nb * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
        cmc_q_boundary_d<<<grid, blk>>>(nb, g_nSlice, g_nCellsAll, isNode, isInlet ? 1 : 0, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                                         var.c_d.at("ro"), g_Q, g_roQ);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    cmcProbe("after-boundary");
}

void cmcQUpdate_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_q_ready) return;
    cmcProbe("before-transport");
    g_timerU.begin();
    // (0) roQ を現在の ρ̄ で作り直す (ρ̄ は直前の流れ更新で変わっている)
    { const size_t n = (size_t)g_nCellsAll * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
      cmc_q_resync_d<<<grid, blk>>>(g_nCellsAll, g_nSlice, var.c_d.at("ro"), g_Q, g_roQ); }
    g_timerU.mark(1);
    // (1) 物理空間の輸送 point-implicit (各スライス; N=M=現在値のエイリアス)
    if (cfg.timeIntegration == 11) {
        const size_t n = (size_t)g_nCells * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
        cmc_q_timeint_d<<<grid, blk>>>(g_nCells, g_nSlice, g_nCellsAll, var.c_d.at("dt_local"), var.c_d.at("volume"), g_roQ, g_resQ, g_sjQ, g_diagQ);
    } else {
        for (int s = 0; s < g_nSlice; ++s) scalarTimeIntegration_d(0, cfg, cuda_cfg, msh, var, sliceDesc(s, cfg));
    }
    g_timerU.mark(2);
    { const size_t n = (size_t)g_nCells * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
      cmc_q_primitive_d<<<grid, blk>>>(g_nCells, g_nSlice, g_nCellsAll, var.c_d.at("ro"), g_roQ, g_Q); }
    g_timerU.mark(3);
    cmcProbe("after-transport");
    if (cfg.outStepInterval > 0 && ((g_cmcStep + 1) % cfg.outStepInterval) == 0) cmcQWrite("cmc_Q.bin");   // 最新の Q を上書き保存 (restart 用)
    // (2) η 拡散 + 化学 + PDF 平均 (interval 毎)
    ++g_cmcStep;
    if (cfg.chemCmcInterval <= 1 || (g_cmcStep % cfg.chemCmcInterval) == 0) {
        cmcStepLaunch(cuda_cfg, var, 1);
        g_timerU.mark(4);
        if (g_cmc_host.couple >= 4 && g_cmc_host.couple != 7) {
            cmc_blend_species_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                g_nCells, thermo_species_device_ptr(), g_ypdf, var.c_d.at("ro"), species_roY_device_ptr(), var.c_d.at("roe"),
                var.c_d.at("chemQdot"), var.c_d.at("dt_local"), g_dpdf, g_dpdf_prev, g_qrel, g_qdebt, g_qcap,
                (g_cmc_host.couple >= 5) ? 1 : 0, g_cmc_host.couple, g_hpdf, var.c_d.at("cmc_xiOm"),
                var.c_d.at("roUx"), var.c_d.at("roUy"), var.c_d.at("roUz"), var.c_d.at("P"));
        }
        if (g_cmc_host.couple == 3) {
            cmc_overwrite_mean_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                g_nCells, thermo_species_device_ptr(), g_ypdf, g_hpdf,
                var.c_d.at("ro"), var.c_d.at("roUx"), var.c_d.at("roUy"), var.c_d.at("roUz"),
                species_roY_device_ptr(), var.c_d.at("roe"), g_cmc_host.Tmax);
        }
    } else { g_timerU.mark(4); }
    g_timerU.mark(5);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    { static const char* nm[] = {"resync", "timeint", "primitive", "cmc_step", "blend"}; g_timerU.end(5, nm); }
}
