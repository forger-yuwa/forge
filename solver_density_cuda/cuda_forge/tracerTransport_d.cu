#include "tracerTransport_d.cuh"

#include "scalarTransport_d.cuh"
#include "passiveTransport_d.cuh"   // passiveScalarScheme 1: 化学種経路の受動種として処理

#include <string>

namespace {

constexpr flow_float kSmall = static_cast<flow_float>(1.0e-30);

inline bool tracerEnabled(const variables& var)
{
    return var.tracerRegistered != 0;
}

// roXi 用のスカラ輸送記述子 (condensation buildCondMomentDesc と同形)。floor=0, 拡散なし, ソースなし。
ScalarTransportDesc buildTracerDesc(variables& var)
{
    return ScalarTransportDesc{
        var.c_d["Xi"], nullptr, nullptr, nullptr,   // dphidx/y/z: 汎用拡散 (diffusion=0) 未使用
        var.c_d["roXi"], var.c_d["roXiN"], var.c_d["roXiM"],
        var.c_d["res_roXi"], var.c_d["res_roXi_m"],
        var.c_d["src_jac_Xi"], var.c_d["transport_diag_Xi"],
        static_cast<flow_float>(0.0), static_cast<flow_float>(0.0),
        0
    };
}

// 原始量 ξ = ρξ/ρ (全セル, ghost 含む)。実現可能性 0 <= ρξ <= ρ を同時にクランプする。
__global__ void tracer_primitive_d(
    geom_int nCells_all,
    flow_float* ro,
    flow_float* roXi,
    flow_float* Xi)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        const flow_float r = max(ro[ic], kSmall);
        flow_float v = roXi[ic];
        if (v < static_cast<flow_float>(0.0)) v = static_cast<flow_float>(0.0);
        if (v > r) v = r;
        roXi[ic] = v;
        Xi[ic]   = v / r;
    }
}

// Neumann (zero-gradient) ghost 充填。
__global__ void tracer_neumann_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* roXi,
    flow_float* Xi)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];
        roXi[ig] = roXi[ic];
        Xi[ig]   = Xi[ic];
    }
}

// Dirichlet ghost 充填 (入口): Xi[ig]=ξ_in, roXi[ig]=ρ[ig]·ξ_in。ρ[ig] は applyBconds が設定済みの ghost 密度。
// node-centered では境界半割面が ghost を読まず境界ノード自身の値を使う (scalarTransport の nodeBnd 分岐) ので、
// 化学種 (species_dirichlet_boundary_d) と同じく境界ノードを Dirichlet 値にピンし scalarDirichletPin で残差を除外する。
__global__ void tracer_dirichlet_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* ro,
    flow_float* Xib,
    flow_float* roXi,
    flow_float* Xi,
    flow_float* scalarDirichletPin,
    int isNode)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ig = bplane_cell_ghst[ib];
        flow_float xin = Xib[ib];
        if (xin < static_cast<flow_float>(0.0)) xin = static_cast<flow_float>(0.0);
        if (xin > static_cast<flow_float>(1.0)) xin = static_cast<flow_float>(1.0);
        Xi[ig]   = xin;
        roXi[ig] = ro[ig] * xin;
        if (isNode != 0) {
            const geom_int ic = bplane_cell[ib];
            Xi[ic]   = xin;
            roXi[ic] = ro[ic] * xin;
            scalarDirichletPin[ic] = static_cast<flow_float>(1.0);
        }
    }
}

// node 入口ピンノードの残差・ソース Jacobian を 0 化 (species_pin_residual_d と同形)。
__global__ void tracer_pin_residual_d(
    geom_int nCells,
    flow_float* res_roXi,
    flow_float* src_jac,
    flow_float* scalarDirichletPin)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells && scalarDirichletPin[ic] == static_cast<flow_float>(1.0)) {
        res_roXi[ic] = static_cast<flow_float>(0.0);
        src_jac[ic]  = static_cast<flow_float>(0.0);
    }
}

}  // namespace

void tracerPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!tracerEnabled(var)) return;
    if (passiveSchemeEnabled(cfg)) {
        // 受動種経路: primitive 段で保存量を書き換えるクランプは撤廃 (上下限は更新確定時の passive_bounds_d)。
        passivePrimitive_d_wrapper(cfg, cuda_cfg, msh, var, passive_tracer_index(), 1);
        return;
    }
    tracer_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells_all, var.c_d["ro"], var.c_d["roXi"], var.c_d["Xi"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void tracerBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var)
{
    (void)msh;
    if (!tracerEnabled(var)) return;
    if (bc.iPlanes.empty()) return;

    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
    const auto xbIt = bc.bvar_d.find("Xi");   // readBcondConfig が inlet_* に登録 (既定 0)
    if (passiveSchemeEnabled(cfg)) {
        // 受動種経路: 化学種の Dirichlet/Neumann カーネルを受動種ポインタで (入口は bvar Xi; node は境界ノードをピン)。
        passiveBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var, passive_tracer_index(),
                                  (isInlet && xbIt != bc.bvar_d.end()) ? xbIt->second : nullptr);
        return;
    }
    if (isInlet && xbIt != bc.bvar_d.end()) {
        tracer_dirichlet_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
            nb, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
            var.c_d["ro"], xbIt->second, var.c_d["roXi"], var.c_d["Xi"],
            var.c_d["scalarDirichletPin"], (cfg.discretization == "node") ? 1 : 0);
    } else {
        tracer_neumann_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
            nb, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
            var.c_d["roXi"], var.c_d["Xi"]);
    }
}

void applyTracerBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!tracerEnabled(var)) return;
    for (auto& bc : msh.bconds) tracerBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void tracerPinResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!tracerEnabled(var) || cfg.discretization != "node") return;
    tracer_pin_residual_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.c_d["res_roXi"], var.c_d["src_jac_Xi"], var.c_d["scalarDirichletPin"]);
    gpuErrchk( cudaPeekAtLastError() );
}

void tracerTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!tracerEnabled(var)) return;
    if (passiveSchemeEnabled(cfg)) {
        // 受動種経路: S3 面値 (SLAU) または 1 次風上 + トレーサ Fick 拡散 (粘性 run)。入口ピンは passivePinResidual (main) で全受動種一括。
        const int q = passive_tracer_index();
        passiveAdvection_d_wrapper(cfg, cuda_cfg, msh, var, q, 1);
        passiveDiffusion_d_wrapper(cfg, cuda_cfg, msh, var, q);
        return;
    }

    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roXi"], 0, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_Xi"], 0, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["src_jac_Xi"], 0, msh.nCells * sizeof(flow_float)));

    const ScalarTransportDesc desc = buildTracerDesc(var);
    scalarTransportResidual_d(cfg, cuda_cfg, msh, var, desc);
    tracerPinResidual_d_wrapper(cfg, cuda_cfg, msh, var);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void tracerTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!tracerEnabled(var)) return;
    if (passiveSchemeEnabled(cfg)) {
        // 受動種経路: point-implicit (passiveImplicitRelax) / scalar-DPLUR 増分 → 上下限 0<=ρξ<=ρ (更新済み ρ) と補正収支。
        passiveTracerUpdate_d_wrapper(loop, cfg, cuda_cfg, msh, var);
        passiveMirrorPeriodic_d_wrapper(cfg, cuda_cfg, msh, var);
        return;
    }
    const ScalarTransportDesc desc = buildTracerDesc(var);
    scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void tracerUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!tracerEnabled(var)) return;
    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    gpuErrchk( cudaMemcpy(var.c_d["roXiN"], var.c_d["roXi"], bytes, cudaMemcpyDeviceToDevice) );
    gpuErrchk( cudaMemcpy(var.c_d["roXiM"], var.c_d["roXi"], bytes, cudaMemcpyDeviceToDevice) );
}

void tracerUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!tracerEnabled(var)) return;
    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    gpuErrchk( cudaMemcpy(var.c_d["roXiM"], var.c_d["roXi"], bytes, cudaMemcpyDeviceToDevice) );
}
