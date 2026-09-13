#include "ransTransport_d.cuh"

#include "scalarTransport_d.cuh"

#include <array>

namespace {

inline bool ransTransportEnabled(const solverConfig& cfg)
{
    return cfg.LESorRANS == 2 && cfg.RANSmodel == 1;
}

// SST k/ω の Green-Gauss 勾配 (源項の交差拡散・F1 ブレンドで使用)。
__global__ void calc_scalar_gradient_face_d(
    geom_int nPlanes,
    geom_int nCells,
    int isNode,
    geom_int* plane_cells,
    geom_float* fx,
    geom_float* sx,
    geom_float* sy,
    geom_float* sz,
    flow_float* k,
    flow_float* omega,
    flow_float* dKdx,
    flow_float* dKdy,
    flow_float* dKdz,
    flow_float* dOmegadx,
    flow_float* dOmegady,
    flow_float* dOmegadz)
{
    geom_int ip = blockDim.x * blockIdx.x + threadIdx.x;
    if (ip >= nPlanes) return;

    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];

    geom_float f   = fx[ip];
    flow_float kf, wf;
    // node 境界面 (ic1=ghost): fx=dc2p/dcc が ghost mirror dcc≈0 で退化する。ghost を一切参照せず、
    // 境界ノード自身の値 (=境界上の値) を面値に使う (ghostless)。壁では omega[ic0]=omegab (ピン留め)・
    // k[ic0]=kb (zero-grad) なので実質 kb/omegab。Dirichlet/Neumann いずれも境界ノード値が境界値。
    if (isNode != 0 && ic1 >= nCells) {
        kf = k[ic0];
        wf = omega[ic0];
    } else {
        kf = f * k[ic0]     + (1.0f - f) * k[ic1];
        wf = f * omega[ic0] + (1.0f - f) * omega[ic1];
    }
    const geom_float sxx = sx[ip];
    const geom_float syy = sy[ip];
    const geom_float szz = sz[ip];

    atomicAdd(&dKdx[ic0],     sxx * kf);
    atomicAdd(&dKdy[ic0],     syy * kf);
    atomicAdd(&dKdz[ic0],     szz * kf);
    atomicAdd(&dOmegadx[ic0], sxx * wf);
    atomicAdd(&dOmegady[ic0], syy * wf);
    atomicAdd(&dOmegadz[ic0], szz * wf);

    atomicAdd(&dKdx[ic1],     -sxx * kf);
    atomicAdd(&dKdy[ic1],     -syy * kf);
    atomicAdd(&dKdz[ic1],     -szz * kf);
    atomicAdd(&dOmegadx[ic1], -sxx * wf);
    atomicAdd(&dOmegady[ic1], -syy * wf);
    atomicAdd(&dOmegadz[ic1], -szz * wf);
}

__global__ void calc_scalar_gradient_div_vol_d(
    geom_int nCells,
    geom_float* vol,
    flow_float* dKdx,
    flow_float* dKdy,
    flow_float* dKdz,
    flow_float* dOmegadx,
    flow_float* dOmegady,
    flow_float* dOmegadz)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const geom_float v = vol[ic];
    dKdx[ic]     /= v;
    dKdy[ic]     /= v;
    dKdz[ic]     /= v;
    dOmegadx[ic] /= v;
    dOmegady[ic] /= v;
    dOmegadz[ic] /= v;
}

// SST k/ω 2 変数ぶんの輸送記述子を構築する。
// floor: realizability 下限 (ρk≥0, ρω>1e-20)。sigma: 拡散係数スケール (σ_k=0.85, σ_ω=0.5)。
__global__ void fill_const_d(geom_int n, flow_float* a, flow_float v)
{
    const geom_int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) a[i] = v;
}

std::array<ScalarTransportDesc, 2> buildScalarDescs(variables& var, const solverConfig& cfg, cudaConfig& cuda_cfg, geom_int nCells)
{
    // sstSigmaBlend=1: σ_k = F1·0.85 + (1−F1)·1.0, σ_ω = F1·0.5 + (1−F1)·0.856 (Menter SST の正式ブレンド)。
    // 0 (既定): k-ω 側定数 0.85 / 0.5 (現行)。F1 は ransSource が書く sstF1 (初回は 1 で埋める)。
    flow_float* F1 = nullptr;
    if (cfg.sstSigmaBlend != 0 && var.c_d.count("sstF1")) {
        static bool inited = false;
        if (!inited) {
            fill_const_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(nCells, var.c_d["sstF1"], static_cast<flow_float>(1.0));
            inited = true;
        }
        F1 = var.c_d["sstF1"];
    }
    std::array<ScalarTransportDesc, 2> d = {{
        {var.c_d["k"],     var.c_d["dKdx"],     var.c_d["dKdy"],     var.c_d["dKdz"],     var.c_d["roK"], var.c_d["roKN"], var.c_d["roKM"], var.c_d["res_roK"], var.c_d["res_roK_m"], var.c_d["src_jac_k"], var.c_d["transport_diag_k"], static_cast<flow_float>(0.0), static_cast<flow_float>(0.85), 1},
        {var.c_d["omega"], var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"], var.c_d["roOmega"], var.c_d["roOmegaN"], var.c_d["roOmegaM"], var.c_d["res_roOmega"], var.c_d["res_roOmega_m"], var.c_d["src_jac_omega"], var.c_d["transport_diag_omega"], static_cast<flow_float>(1.0e-20), static_cast<flow_float>(0.5), 1}
    }};
    d[0].sigma2 = static_cast<flow_float>(1.0);   d[0].F1 = F1;
    d[1].sigma2 = static_cast<flow_float>(0.856); d[1].F1 = F1;
    return d;
}

}

void ransTransport_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roK"], 0.0f, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roOmega"], 0.0f, msh.nCells * sizeof(flow_float)));
    // 輸送ヤコビアン対角を毎 assembleResidual でゼロ初期化（advection/diffusion kernel で面ごと加算）。
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_k"], 0.0f, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_omega"], 0.0f, msh.nCells * sizeof(flow_float)));

    if (!ransTransportEnabled(cfg)) {
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        return;
    }

    const auto scalar_descs = buildScalarDescs(var, cfg, cuda_cfg, msh.nCells);

    // k/ω を 1 面ループで融合 (面幾何・massflux・μ の読みを共有)。面ごとの式は単一版と同一。
    scalarTransportResidualMulti_d(cfg, cuda_cfg, msh, var, scalar_descs.data(), 2);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void ransTimeIntegration_d_wrapper(int loop , solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!ransTransportEnabled(cfg)) {
        return;
    }

    const auto scalar_descs = buildScalarDescs(var, cfg, cuda_cfg, msh.nCells);

    for (const auto& desc : scalar_descs) {
        scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// sstEnergyIncludesK (plan turbulence-sst-energy-includes-k §4): 保存量 roe は平均流エネルギー E_m のまま保持し、
// エネルギー残差は E_t = E_m + ρk の流束で組む。流れの更新は E_t の更新と見なせるので、k 更新後に
// roe -= Δ(ρk) とすると E_t^{new} = E_t^{old} + Δ·R(E_t) が各 CV で厳密に成立する (k 式の点陰化・床置きの如何によらず)。
// Δ(ρk) の基準は流れの更新形に合わせる:
//   fromN=1 (explicit RK): 各 stage が N 状態から組み直すので Δ = roK − roKN。
//   fromN=0 (point-implicit / dual-time subiter): roe += dq の増分更新なので Δ = roK − (k 更新直前の roK)。
//   → begin で roK を退避し、end で差を引く。
namespace {
flow_float* g_sstEkPrev = nullptr; geom_int g_sstEkPrevN = 0;
__global__ void sst_energy_k_correction_d(geom_int nCells, flow_float* roe, flow_float* roK, flow_float* roKref)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) roe[ic] -= (roK[ic] - roKref[ic]);
}
}

void sstEnergyKCorrection_begin_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (cfg.sstEnergyIncludesK == 0 || !ransTransportEnabled(cfg)) return;
    if (g_sstEkPrev == nullptr || g_sstEkPrevN != msh.nCells) {
        if (g_sstEkPrev != nullptr) cudaFree(g_sstEkPrev);
        CHECK_CUDA_ERROR(cudaMalloc(&g_sstEkPrev, msh.nCells * sizeof(flow_float)));
        g_sstEkPrevN = msh.nCells;
    }
    CHECK_CUDA_ERROR(cudaMemcpy(g_sstEkPrev, var.c_d["roK"], msh.nCells * sizeof(flow_float), cudaMemcpyDeviceToDevice));
}

void sstEnergyKCorrection_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int fromN)
{
    if (cfg.sstEnergyIncludesK == 0 || !ransTransportEnabled(cfg)) return;
    flow_float* ref = (fromN != 0) ? var.c_d["roKN"] : g_sstEkPrev;
    if (ref == nullptr) return;
    sst_energy_k_correction_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(msh.nCells, var.c_d["roe"], var.c_d["roK"], ref);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void ransGradient_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!ransTransportEnabled(cfg)) {
        return;
    }

    flow_float* grad_vol = (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"];
    flow_float* grad_sx  = (cfg.isAxisymmetric == 1) ? var.p_d["sx_planar"] : var.p_d["sx"];
    flow_float* grad_sy  = (cfg.isAxisymmetric == 1) ? var.p_d["sy_planar"] : var.p_d["sy"];
    flow_float* grad_sz  = (cfg.isAxisymmetric == 1) ? var.p_d["sz_planar"] : var.p_d["sz"];

    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dKdx"],     0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dKdy"],     0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dKdz"],     0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dOmegadx"], 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dOmegady"], 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dOmegadz"], 0, msh.nCells_all * sizeof(flow_float)));

    calc_scalar_gradient_face_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nPlanes,
        msh.nCells,
        (cfg.discretization == "node") ? 1 : 0,
        msh.map_plane_cells_d,
        var.p_d["fx"],
        grad_sx, grad_sy, grad_sz,
        var.c_d["k"],
        var.c_d["omega"],
        var.c_d["dKdx"], var.c_d["dKdy"], var.c_d["dKdz"],
        var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"]);

    calc_scalar_gradient_div_vol_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        grad_vol,
        var.c_d["dKdx"], var.c_d["dKdy"], var.c_d["dKdz"],
        var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"]);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
