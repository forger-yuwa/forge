#include "condensationTransport_d.cuh"
#include "condensationUpdateLimiter_d.cuh"   // cond_moment_update_limited_d (単体試験と共用)
#include "condensationSource_d.cuh"   // COND_PI, 物性 (消滅クランプ)
#include "condensationEOS_d.cuh"
#include "condensationSourceF_d.cuh"   // float 実体 (clamp の表評価)      // cond_clamp_vapor_pressure (蒸発塵判定の蒸気分圧)

#include "scalarTransport_d.cuh"

#include <string>
#include <vector>

// device rog (液相質量分率の保存量) ポインタ配列。二相 EOS が読む。condensationInit_d で構築。
static flow_float** g_rog_dev = nullptr;
static CondTablesF g_condTables;   // float 経路の物性表 (condensationInit_d で構築)
static int          g_nCond   = 0;

namespace {

constexpr flow_float kSmall = static_cast<flow_float>(1.0e-30);

inline bool condensationEnabled(const variables& var)
{
    return var.nCondSpeciesRegistered >= 1;
}

// 保存量名 consName ("rog_0" 等) から派生名を組み立てて ScalarTransportDesc を構築する
// (RANS buildScalarDescs / species buildSpeciesDesc と同形)。
// floor=0 (ρφ>=0)、sigma=0、diffusion=0 (Phase 1 は移流のみ)。src_jac=0 (ソースなし)。
ScalarTransportDesc buildCondMomentDesc(variables& var, const std::string& consName)
{
    const std::string prim = consName.substr(2); // 先頭 "ro" を除去 ("g_0","Q0_0")
    return ScalarTransportDesc{
        var.c_d[prim], nullptr, nullptr, nullptr,  // dphidx/y/z: 凝縮モーメントは汎用拡散(diffusion=0)未使用
        var.c_d[consName], var.c_d[consName+"N"], var.c_d[consName+"M"],
        var.c_d["res_"+consName], var.c_d["res_"+consName+"_m"],
        var.c_d["src_jac_"+prim], var.c_d["transport_diag_"+prim],
        static_cast<flow_float>(0.0), static_cast<flow_float>(0.0),
        0
    };
}

// 原始量 φ = ρφ/ρ を最大 4 モーメント同時に (起動数削減; plan condensation-float-speedup §4.2-7)。
struct CondPrimPtrs { flow_float* rophi[4]; flow_float* phi[4]; int n; };
__global__ void cond_primitive_multi_d(geom_int nCells_all, flow_float* ro, CondPrimPtrs P)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        const flow_float inv = static_cast<flow_float>(1.0) / max(ro[ic], kSmall);
        for (int k = 0; k < P.n; ++k) P.phi[k][ic] = P.rophi[k][ic] * inv;
    }
}
// 原始量 φ = ρφ/ρ (全セル, ghost 含む)。
__global__ void cond_primitive_d(
    geom_int nCells_all,
    flow_float* ro,
    flow_float* rophi,
    flow_float* phi)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        phi[ic] = rophi[ic] / max(ro[ic], kSmall);
    }
}

// 実現可能性クランプ: 液相保存量 rog と モーメント roQ0/1/2 を物理範囲に戻す。
//   0 ≤ rog ≤ roY_w (carrier: 利用可能な総水量) または ≤ 0.99 ρ (pure)。roQn ≥ 0。
// θ 律速は瞬間 Sg を抑えるが、陰解法の実効 Δt と dt_local の差で僅かに過凝縮しうるため
// 毎ステップここで硬クランプして二相 EOS・次ステップ source に渡す値を保証する。
// 蒸発 (condEvaporation=1) では液滴消滅の硬クランプも行う (plans/accepted/condensation-evaporation.md §4.4):
//   S=p_v/p_sat<=1 かつ液相あり かつ (Q0=0 [不整合] または r30<2*rmin) かつ g<=g_rm (潜熱飛びが小さい)
//   → rog,roQ0,roQ1,roQ2 を 0。source kernel の λ=0 は陰的緩和で厳密 0 に届かないためここで確定する。
//   T,P は前ステップ値 (dependentVariables 前) で判定する。核生成域 (S>1) には触れない。
__global__ void cond_realizability_clamp_d(
    geom_int nCells,
    flow_float* ro, flow_float* roY_w,   // roY_w: carrier の総水保存量 (pure では nullptr)
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    int evap, int condModel, double Rw, double rmin, double g_rm,
    flow_float* T, flow_float* P, CondPropOpts opts)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    // 実現可能性 g<=Y_w: TP carrier は roY_w (輸送), CPG carrier (空気) は定数 opts.Yw, pure は 0.99
    const flow_float gmax = (roY_w != nullptr) ? roY_w[ic] : ((opts.Yw > 0.0) ? (flow_float)opts.Yw*ro[ic] : (flow_float)0.99*ro[ic]);
    flow_float r = rog[ic];
    if (r < (flow_float)0.0) r = (flow_float)0.0;
    if (r > gmax)            r = gmax;
    rog[ic] = r;
    if (roQ0[ic] < (flow_float)0.0) roQ0[ic] = (flow_float)0.0;
    if (roQ1[ic] < (flow_float)0.0) roQ1[ic] = (flow_float)0.0;
    if (roQ2[ic] < (flow_float)0.0) roQ2[ic] = (flow_float)0.0;

    if (!evap) return;
    const double rod = (double)ro[ic];
    if (rod <= 1.0e-20) return;
    const double g  = (double)r/rod;
    if (g > g_rm) return;
    // g==0 で Q0/Q1/Q2 だけ残る「モーメント塵」(float アンダーフロー) も S<=1 では掃除する
    // (核生成域 S>1 は触らない: 新核の g がアンダーフローしていても Q0 は生かす)。
    const bool dust = (r <= (flow_float)0.0) &&
        (roQ0[ic] > (flow_float)0.0 || roQ1[ic] > (flow_float)0.0 || roQ2[ic] > (flow_float)0.0);
    if (r <= (flow_float)0.0 && !dust) return;
    const CondSpeciesProps cprops = condProps_make(condModel, opts);
    const double Td = (double)T[ic];
    // 蒸気分圧は source kernel (cond_vapor_state) と同じ定義: TP carrier=ρ(Y_w−g)R_wT, CPG carrier=ρ(Y_w,const−g)R_wT, pure=全圧 (codex 2026-09-13 M2)
    const double pv = cond_clamp_vapor_pressure(rod, g, (roY_w != nullptr) ? (double)roY_w[ic]/rod : -1.0, opts.Yw, Rw, Td, (double)P[ic]);
    if (pv > cond_psat(cprops, Td)) return;              // 過飽和: 消滅させない
    const double q0 = (double)roQ0[ic];
    bool remove = dust || (q0 <= 1.0e-30);
    if (!remove) {
        const double rho_l = cond_rho_cond(cprops, Td);
        const double r30 = cbrt(g/((4.0/3.0)*COND_PI*rho_l*q0/rod));
        remove = (r30 < 2.0*rmin);
    }
    if (remove) {
        rog[ic] = (flow_float)0.0; roQ0[ic] = (flow_float)0.0;
        roQ1[ic] = (flow_float)0.0; roQ2[ic] = (flow_float)0.0;
    }
}


// float 実体 (condFloat=1): 判定は閾値比較 (p_v<=p_sat, r30<2 r_min) なので ULP 差で消滅 step が 1 つずれ得る (plan §4.2-4)。
__global__ void cond_realizability_clamp_f_d(
    geom_int nCells,
    flow_float* ro, flow_float* roY_w,
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    int evap, float Rw, float rmin, float g_rm, float Yw_const,
    flow_float* T, flow_float* P, CondTablesF tb, CondSpeciesProps cpd)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float gmax = (roY_w != nullptr) ? roY_w[ic] : ((Yw_const > 0.0f) ? Yw_const*ro[ic] : 0.99f*ro[ic]);
    flow_float r = rog[ic];
    if (r < 0.0f) r = 0.0f;
    if (r > gmax) r = gmax;
    rog[ic] = r;
    if (roQ0[ic] < 0.0f) roQ0[ic] = 0.0f;
    if (roQ1[ic] < 0.0f) roQ1[ic] = 0.0f;
    if (roQ2[ic] < 0.0f) roQ2[ic] = 0.0f;
    if (!evap) return;
    const float rod = ro[ic];
    if (rod <= 1.0e-20f) return;
    const float g = r/rod;
    if (g > g_rm) return;
    const bool dust = (r <= 0.0f) && (roQ0[ic] > 0.0f || roQ1[ic] > 0.0f || roQ2[ic] > 0.0f);
    if (r <= 0.0f && !dust) return;
    const float Td = T[ic];
    // 蒸気分圧 (source kernel と同じ定義): TP carrier=ρ(Y_w−g)R_wT, CPG carrier=ρ(Y_w,const−g)R_wT, pure=全圧
    float pv;
    if (roY_w != nullptr)      { float yv = roY_w[ic]/rod - g; if (yv < 0.0f) yv = 0.0f; pv = rod*yv*Rw*Td; }
    else if (Yw_const > 0.0f)  { float yv = Yw_const - g;     if (yv < 0.0f) yv = 0.0f; pv = rod*yv*Rw*Td; }
    else                         pv = P[ic];
    const bool inTab = cond_tab_wet_ok(tb, Td);   // 表範囲外は旧 double 関数 (plan §5.1 #8)
    const float lnps = inTab ? cond_tab_lnpsat_f(tb, Td) : (float)log(cond_psat(cpd, (double)Td) > 1.0e-300 ? cond_psat(cpd, (double)Td) : 1.0e-300);
    if (pv > 0.0f && logf(pv) > lnps) return;   // 過飽和: 消滅させない
    const float q0 = roQ0[ic];
    bool remove = dust || (q0 <= 1.0e-30f);
    if (!remove) {
        const float rho_l = inTab ? cond_tab_rhol_f(tb, Td) : (float)cond_rho_cond(cpd, (double)Td);
        const float r30 = cbrtf(g/((4.0f/3.0f)*COND_PI_F*rho_l*q0/rod));
        remove = (r30 < 2.0f*rmin);
    }
    if (remove) { rog[ic] = 0.0f; roQ0[ic] = 0.0f; roQ1[ic] = 0.0f; roQ2[ic] = 0.0f; }
}

// Neumann (zero-gradient) ghost 充填: rophi[ig]=rophi[ic], phi[ig]=phi[ic]。
__global__ void cond_neumann_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* rophi,
    flow_float* phi)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];
        rophi[ig] = rophi[ic];
        phi[ig]   = phi[ic];
    }
}

// Dirichlet (dry 入口) ghost 充填: 液相モーメントは入口で 0 (ρφ=0, φ=0)。
__global__ void cond_dirichlet_zero_boundary_d(
    geom_int nb,
    geom_int* bplane_cell_ghst,
    flow_float* rophi,
    flow_float* phi)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ig = bplane_cell_ghst[ib];
        rophi[ig] = static_cast<flow_float>(0.0);
        phi[ig]   = static_cast<flow_float>(0.0);
    }
}

}  // namespace

void condensationInit_d(solverConfig& cfg, variables& var)
{
    (void)cfg;
    if (!condensationEnabled(var)) {
        g_rog_dev = nullptr;
        g_nCond = 0;
        return;
    }
    g_nCond = var.nCondSpeciesRegistered;

    // 各凝縮種の保存量 rog_{s} の device ポインタを集める (registerCondensation の命名規約)。
    std::vector<flow_float*> hrog(g_nCond);
    for (int s = 0; s < g_nCond; ++s) {
        hrog[s] = var.c_d["rog_" + std::to_string(s)];
    }
    const size_t pbytes = g_nCond * sizeof(flow_float*);
    gpuErrchk( cudaMalloc((void**)&g_rog_dev, pbytes) );
    gpuErrchk( cudaMemcpy(g_rog_dev, hrog.data(), pbytes, cudaMemcpyHostToDevice) );

    std::cout << "condensationInit_d: built device rog[] for nCondSpecies=" << g_nCond << "\n";
    // float 経路の物性表 (plans/active/condensation-float-speedup.md §4.2-1): 現行 double 関数から区分 3 次表を作り device へ。
    if (cfg.condFloat != 0) {
        CondTablesHost ht;
        cond_tables_build_host(condProps_make(cfg.condModel, cond_prop_opts(cfg)), ht);
        g_condTables = cond_tables_upload(ht);
        std::cout << "condensationInit_d: property tables for condFloat (model " << cfg.condModel << ", T0=" << ht.T0
                  << " h=" << ht.h << " n=" << ht.n << ", " << (6*ht.n*sizeof(float4))/1024 << " KB)\n";
    } else {
        g_condTables = CondTablesF();
        std::cout << "condensationInit_d: condFloat=0 (double condensation path)\n";
    }
}
const CondTablesF& cond_tables_device() { return g_condTables; }

flow_float** cond_rog_device_ptr() { return g_rog_dev; }
int          cond_num_species()    { return g_nCond; }

void condensationPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!condensationEnabled(var)) return;

    // 実現可能性クランプ (種ごと): 0 ≤ rog ≤ roY_w (carrier) / 0.99ρ (pure)、roQn ≥ 0。
    const bool carrier = (cfg.condGasSpecies >= 0);
    for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        flow_float* roY_w = carrier ? var.c_d["roY" + std::to_string(cfg.condGasSpecies)] : nullptr;
        const CondPropOpts opts = cond_prop_opts(cfg);
        const CondSpeciesProps cprops = condProps_make(cfg.condModel, opts);
        const double g_rm = 5.0e-7;   // 消滅硬クランプを許す g 上限 (潜熱飛び ΔT=gL/cv ≲ 1.5 K)
        if (cfg.condFloat != 0 && g_condTables.valid) {
            cond_realizability_clamp_f_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                msh.nCells, var.c_d["ro"], roY_w,
                var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
                cfg.condEvaporation, (float)cprops.R, (float)cfg.condEvapRmin, (float)g_rm, (float)opts.Yw,
                var.c_d["T"], var.c_d["P"], g_condTables, cprops);
        } else
        cond_realizability_clamp_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells, var.c_d["ro"], roY_w,
            var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
            cfg.condEvaporation, cfg.condModel, cprops.R, cfg.condEvapRmin, g_rm,
            var.c_d["T"], var.c_d["P"], opts);
    }

    {
        // 4 モーメントずつ 1 起動 (φ=ρφ/ρ は除算 1 回を逆数乗算に; 値は 1 ulp 以内)
        CondPrimPtrs P{}; P.n = 0;
        for (const auto& consName : var.condMomentConsNames) {
            const std::string prim = consName.substr(2);
            P.rophi[P.n] = var.c_d[consName]; P.phi[P.n] = var.c_d[prim]; ++P.n;
            if (P.n == 4) { cond_primitive_multi_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d["ro"], P); P.n = 0; }
        }
        if (P.n > 0) cond_primitive_multi_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d["ro"], P);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void condensationBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var)
{
    (void)cfg; (void)msh;
    if (!condensationEnabled(var)) return;
    if (bc.iPlanes.empty()) return;

    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    // 入口は dry (液相モーメント=0) の Dirichlet。他種別 (outlet/slip/wall/axis) は zero-gradient。
    const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;

    for (const auto& consName : var.condMomentConsNames) {
        const std::string prim = consName.substr(2);
        if (isInlet) {
            cond_dirichlet_zero_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nb,
                bc.map_bplane_cell_ghst_d,
                var.c_d[consName],
                var.c_d[prim]);
        } else {
            cond_neumann_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nb,
                bc.map_bplane_cell_d,
                bc.map_bplane_cell_ghst_d,
                var.c_d[consName],
                var.c_d[prim]);
        }
    }
}

void applyCondensationBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condensationEnabled(var)) return;

    for (auto& bc : msh.bconds) {
        condensationBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void condensationTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condensationEnabled(var)) return;

    for (const auto& consName : var.condMomentConsNames) {
        const std::string prim = consName.substr(2);
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_"+consName], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_"+prim], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["src_jac_"+prim], 0, msh.nCells * sizeof(flow_float)));
    }

    // 4 モーメント (ρg, ρQ0, ρQ1, ρQ2) の 1 次風上移流を 1 面ループに融合 (k/ω と同じ MultiScalarPtrs; 拡散なし)。
    // 面ごとの流束の算術は単独版と同一 (起動数 4→1; plan condensation-float-speedup §4.2-7)。
    std::vector<ScalarTransportDesc> descs;
    for (const auto& consName : var.condMomentConsNames) descs.push_back(buildCondMomentDesc(var, consName));
    scalarTransportResidualMulti_d(cfg, cuda_cfg, msh, var, descs.data(), (int)descs.size());

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void condensationTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condensationEnabled(var)) return;

    if (cfg.timeIntegration == 11 && cfg.condLimiterMode == 1) {
        // 更新クランプ経路 (plans/active/condensation-source-limiter-steady.md): 種ごとに 4 モーメントをまとめて更新。
        // condMomentConsNames の順序は registerCondensation の bases = {g, Q2, Q1, Q0} (種ごとに 4 本連続)。
        const int carrier = (cfg.condGasSpecies >= 0 || cfg.condVaporMassFraction > 0.0) ? 1 : 0;
        const CondPropOpts opts = cond_prop_opts(cfg);
        flow_float* roY_w = (carrier && cfg.condGasSpecies >= 0) ? var.c_d["roY" + std::to_string(cfg.condGasSpecies)] : nullptr;
        flow_float* cp_cell   = (cfg.thermalMethod == 2) ? var.c_d["cp"]   : nullptr;
        flow_float* Rmix_cell = (cfg.thermalMethod == 2) ? var.c_d["Rmix"] : nullptr;
        const double lam_min = 0.5;
        for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
            const std::string i = std::to_string(s);
            const std::string g = "rog_"+i, Q2 = "roQ2_"+i, Q1 = "roQ1_"+i, Q0 = "roQ0_"+i;
            cond_moment_update_limited_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                msh.nCells, var.c_d["dt_local"], var.c_d["volume"], var.c_d["ro"],
                roY_w, cfg.condVaporMassFraction, var.c_d["T"], cp_cell, Rmix_cell, cfg.cp, cfg.gamma,
                cfg.condModel, opts, cfg.condDgMaxStep, cfg.condDTmaxStep, lam_min,
                var.c_d[g+"N"], var.c_d[Q2+"N"], var.c_d[Q1+"N"], var.c_d[Q0+"N"],
                var.c_d["res_"+g], var.c_d["res_"+Q2], var.c_d["res_"+Q1], var.c_d["res_"+Q0],
                var.c_d["src_jac_g_"+i], var.c_d["src_jac_Q2_"+i], var.c_d["src_jac_Q1_"+i], var.c_d["src_jac_Q0_"+i],
                var.c_d["transport_diag_g_"+i], var.c_d["transport_diag_Q2_"+i], var.c_d["transport_diag_Q1_"+i], var.c_d["transport_diag_Q0_"+i],
                var.c_d[g], var.c_d[Q2], var.c_d[Q1], var.c_d[Q0],
                var.c_d["condLim_"+i], var.c_d["condClampCorr_"+i]);
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        return;
    }

    for (const auto& consName : var.condMomentConsNames) {
        const ScalarTransportDesc desc = buildCondMomentDesc(var, consName);
        scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void condensationUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!condensationEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (const auto& consName : var.condMomentConsNames) {
        gpuErrchk( cudaMemcpy(var.c_d[consName+"N"], var.c_d[consName], bytes, cudaMemcpyDeviceToDevice) );
        gpuErrchk( cudaMemcpy(var.c_d[consName+"M"], var.c_d[consName], bytes, cudaMemcpyDeviceToDevice) );
    }
}

void condensationUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!condensationEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (const auto& consName : var.condMomentConsNames) {
        gpuErrchk( cudaMemcpy(var.c_d[consName+"M"], var.c_d[consName], bytes, cudaMemcpyDeviceToDevice) );
    }
}
