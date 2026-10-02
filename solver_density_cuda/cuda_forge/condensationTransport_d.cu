#include "condensationTransport_d.cuh"
#include "gasPhaseComposition_d.cuh"   // gasPhaseLiquid (輸送物性の気相組成)
#include "condensationUpdateLimiter_d.cuh"   // cond_moment_update_limited_d (単体試験と共用)
#include "condensationRealizability_d.cuh"
#include "periodicNode_d.cuh"   // periodicNodeActive (実現可能性収支の root)
#include <algorithm>   // cond_realizability_clamp_d (double; 単体試験と共用)
#include <iterator>    // std::back_inserter (二相更新の診断 CSV)
#include "condensationSource_d.cuh"   // COND_PI, 物性 (消滅クランプ)
#include "condensationEOS_d.cuh"
#include "condensationSourceF_d.cuh"   // float 実体 (clamp の表評価)      // cond_clamp_vapor_pressure (蒸発塵判定の蒸気分圧)

#include "scalarTransport_d.cuh"
#include "passiveTransport_d.cuh"   // passiveScalarScheme 1: 化学種経路の受動種として移流・更新
#include "speciesTransport_d.cuh"   // 二相拡散 (twoPhaseDiffusion_d_wrapper / speciesRenormalizeTwoPhase_d_wrapper; #4e)
#include "twoPhaseDiffusion_d.cuh"   // tp_vl_update (非分割更新の 1 セル本体)

#include "input/speciesDB.hpp"   // 潜熱の気液ペア (speciesDB_current()->condensed; plan #10)

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// device rog (液相質量分率の保存量) ポインタ配列。二相 EOS が読む。condensationInit_d で構築。
static flow_float** g_rog_dev = nullptr;
static CondTablesF g_condTables;   // float 経路の物性表 (condensationInit_d で構築)
static int          g_nCond   = 0;

CondLatentRef cond_latent_pair_for(const solverConfig& cfg)
{
    if (!(cfg.condensation == 1 && cfg.condModel == COND_MODEL_H2O)) return CondLatentRef{};
    static bool           built = false, checked = false;
    static CondLatentPair cached;
    static CondLatentRef  ref;
    if (!built) {
        const ResolvedSpeciesDB* db = speciesDB_current();
        if (!db) db = &speciesDB_init(cfg);
        const ResolvedCondensed& c = db->condensed;
        if (!c.enabled) {
            std::fprintf(stderr, "[cond] ERROR: condensation of H2O needs the liquid phase H2O(L) paired with the gas in the species DB "
                                 "(speciesDB_resolve did not attach it; plan thermophysics-solver-owned-species-db #10)\n");
            std::exit(EXIT_FAILURE);
        }
        // datum は thermo_init_db と同じ条件 (TP かつ thermoHrefTemp>0) で、気相・液相に同じ Δa7 を掛ける
        const double Tref = (cfg.thermalMethod == 2 && cfg.thermoHrefTemp > 0.0) ? cfg.thermoHrefTemp : 0.0;
        cached = cond_latent_pair_make(c.gas, c.coeffs, c.Tlo, c.Thi, Tref);
        CondLatentPair* dev = nullptr;   // device 側の複製 (プロセス終了まで保持)
        gpuErrchk( cudaMalloc((void**)&dev, sizeof(CondLatentPair)) );
        gpuErrchk( cudaMemcpy(dev, &cached, sizeof(CondLatentPair), cudaMemcpyHostToDevice) );
        ref.h = &cached; ref.d = dev;
        built = true;
        const std::string where = (c.gasIndex >= 0) ? "species " + std::to_string(c.gasIndex) : std::string("built-in, not in the species list");
        std::printf("[cond] H2O latent heat L = h_v - h_l (plan #10): gas '%s' (%s), liquid %s %.2f-%.2f K, datum Tref=%g K;"
                    " L(150/200/250/300 K) = %.1f / %.1f / %.1f / %.1f J/kg\n",
                    c.gasName.c_str(), where.c_str(), c.name.c_str(), c.Tlo, c.Thi, Tref,
                    h2o_latent_pair(cached, 150.0), h2o_latent_pair(cached, 200.0), h2o_latent_pair(cached, 250.0), h2o_latent_pair(cached, 300.0));
    }
    // 気相が種 DB の device 係数 (thermo_init_db が datum を焼き込んだもの) とビット一致することを 1 回確かめる (TP のみ)
    const int gi = speciesDB_current() ? speciesDB_current()->condensed.gasIndex : -1;
    if (!checked && cfg.thermalMethod == 2 && gi >= 0 && gi < thermo_num_species() && thermo_species_host() != nullptr) {
        const SpeciesThermo& g = thermo_species_host()[gi];
        const bool same = thermo_same_coeffs(g, cached.gas);
        if (!same) {
            std::fprintf(stderr, "[cond] ERROR: the gas of the H2O latent-heat pair differs from species %d of the species DB after the datum "
                                 "offset (the latent heat must use the same evaluation as the species DB; plan #10)\n", gi);
            std::exit(EXIT_FAILURE);
        }
        checked = true;
    }
    return ref;
}

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
    std::cout << "[cond-corr] per-reason correction budgets (plan condensation-two-phase-transport §4.3) start at 0 in this process;"
                 " cumulative values do not carry over a restart\n";
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

// 輸送物性の気相組成に使う液 (plan condensation-two-phase-transport §4.1)。TP carrier (凝縮種が輸送種) の凝縮 run だけ
// {rog_0, condGasSpecies}。CPG carrier (condVaporMassFraction; 凝縮種は輸送種でない)・pure 凝縮・凝縮 OFF は {nullptr, -1}。
GasPhaseLiquid gasPhaseLiquid(const solverConfig& cfg, variables& var)
{
    if (!condensationEnabled(var) || cfg.thermalMethod != 2 || cfg.condGasSpecies < 0 || cfg.nSpecies < 2)
        return GasPhaseLiquid{nullptr, -1};
    return GasPhaseLiquid{var.c_d["rog_0"], cfg.condGasSpecies};
}
int          cond_num_species()    { return g_nCond; }

// モーメント不等式の違反数 (実現可能性射影の作動数; device カウンタ)。monitor ログが読んで 0 に戻す。
static int* g_realizViol_dev = nullptr;
int* condRealizViolCounter()
{
    if (g_realizViol_dev == nullptr) { gpuErrchk( cudaMalloc((void**)&g_realizViol_dev, 2*sizeof(int)) ); gpuErrchk( cudaMemset(g_realizViol_dev, 0, 2*sizeof(int)) ); }
    return g_realizViol_dev;
}
// 実現可能性クランプ (射影 + 非負化 + 液相上限 + 消滅) の成分別収支 [種 s][g,Q0,Q1,Q2]×(符号付き, 絶対)·V (root のみ; 全期間積算)。
static double* g_clampBudget_dev = nullptr; static int g_clampBudget_n = 0;
double* condClampBudget(int s)
{
    if (g_clampBudget_dev == nullptr) { g_clampBudget_n = 8; gpuErrchk( cudaMalloc((void**)&g_clampBudget_dev, (size_t)g_clampBudget_n*8*sizeof(double)) ); gpuErrchk( cudaMemset(g_clampBudget_dev, 0, (size_t)g_clampBudget_n*8*sizeof(double)) ); }
    return (s < g_clampBudget_n) ? g_clampBudget_dev + (size_t)s*8 : nullptr;
}
std::vector<double> condClampBudgetTotals(int nSpecies)
{
    std::vector<double> h((size_t)std::max(0, std::min(nSpecies, g_clampBudget_n))*8, 0.0);
    if (g_clampBudget_dev != nullptr && !h.empty()) gpuErrchk( cudaMemcpy(h.data(), g_clampBudget_dev, h.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    return h;
}
// 理由別の補正量 (§4.3, #2)。[s*COND_REASON_N + k]。和・数は全期間 (このプロセスの開始から; restart で 0)、VMIN/RN_MAX は区間値。
static double* g_condReasons_dev = nullptr; static int g_condReasons_n = 0;
static std::vector<double> g_condReasons_last;   // 前回ログ時点の値 (区間値 = 現在 − 前回)
static int g_condReasons_last_step = 0;
static void condReasonsResetInterval(int n)
{
    for (int s = 0; s < n; ++s) {
        const double init[3] = {1.0e300, 0.0, 1.0};
        gpuErrchk( cudaMemcpy(g_condReasons_dev + (size_t)s*COND_REASON_N + COND_REASON_VMIN,   &init[0], sizeof(double), cudaMemcpyHostToDevice) );
        gpuErrchk( cudaMemcpy(g_condReasons_dev + (size_t)s*COND_REASON_N + COND_REASON_RN_MAX, &init[1], sizeof(double), cudaMemcpyHostToDevice) );
        gpuErrchk( cudaMemcpy(g_condReasons_dev + (size_t)s*COND_REASON_N + COND_REASON_TU_MIN, &init[2], sizeof(double), cudaMemcpyHostToDevice) );
        gpuErrchk( cudaMemcpy(g_condReasons_dev + (size_t)s*COND_REASON_N + COND_REASON_TS_MIN, &init[2], sizeof(double), cudaMemcpyHostToDevice) );
    }
}
double* condCorrReasons(int s)
{
    if (g_nCond <= 0 || s < 0 || s >= g_nCond) return nullptr;
    if (g_condReasons_dev == nullptr) {
        g_condReasons_n = g_nCond;
        gpuErrchk( cudaMalloc((void**)&g_condReasons_dev, (size_t)g_condReasons_n*COND_REASON_N*sizeof(double)) );
        gpuErrchk( cudaMemset(g_condReasons_dev, 0, (size_t)g_condReasons_n*COND_REASON_N*sizeof(double)) );
        condReasonsResetInterval(g_condReasons_n);
        g_condReasons_last.assign((size_t)g_condReasons_n*COND_REASON_N, 0.0);
    }
    return g_condReasons_dev + (size_t)s*COND_REASON_N;
}
double* condCorrReasonsForRenormalize(const solverConfig& cfg, variables& var, int* iw)
{
    const GasPhaseLiquid L = gasPhaseLiquid(cfg, var);
    *iw = L.iw;
    return (L.iw >= 0) ? condCorrReasons(0) : nullptr;
}

int condRealizViolReadReset(int* degenerate)
{
    if (g_realizViol_dev == nullptr) { if (degenerate) *degenerate = 0; return 0; }
    int h[2] = {0, 0}; gpuErrchk( cudaMemcpy(h, g_realizViol_dev, 2*sizeof(int), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemset(g_realizViol_dev, 0, 2*sizeof(int)) );
    if (degenerate) *degenerate = h[1];
    return h[0];
}

// Q1/Q2 だけの実現可能性射影 (EOS 更新後の T; g は不変) + primitive の再同期 (plan §4.7 v5, plan-7 M4)。
void condensationRealizabilityProject_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    // 新しい実現可能性射影は**受動種経路 (passiveScalarScheme 1) 限定** (codex result-3 M3): 旧経路 (scheme 0) の
    // 状態を変えない契約 (§6-7 無影響) を守る。射影は Q1/Q2 を実際に動かす (例: Q0>0, g>0, Q1=Q2=0 → 単分散再初期化)
    // ので、丸め差ではなく state の変更になる。
    if (!condensationEnabled(var) || cfg.condRealizProject == 0 || cfg.passiveScalarScheme != 1) return;
    const CondPropOpts opts = cond_prop_opts(cfg);
    const int useTab = 0;   // 射影は step 末尾の 1 回なので物性は double の式で (float 表との差 [~1e-4] がゲートの許容 1e-6 を超えないように; plan-10 M4)
    for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        cond_realizability_project_only_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells, var.c_d["ro"], var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
            cfg.condModel, var.c_d["T"], opts, useTab, g_condTables,
            var.c_d["condClampCorrQ_"+i], condRealizViolCounter(), condClampBudget(s), periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr, var.c_d["volume"],
            condCorrReasons(s));
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    // primitive を再同期 (Q1, Q2)
    {
        CondPrimPtrs P{}; P.n = 0;
        for (const auto& consName : var.condMomentConsNames) {
            const std::string prim = consName.substr(2);
            P.rophi[P.n] = var.c_d[consName]; P.phi[P.n] = var.c_d[prim]; ++P.n;
            if (P.n == 4) { cond_primitive_multi_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d["ro"], P); P.n = 0; }
        }
        if (P.n > 0) cond_primitive_multi_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d["ro"], P);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void condensationPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!condensationEnabled(var)) return;

    // 実現可能性クランプ (種ごと): 0 ≤ rog ≤ roY_w (carrier) / 0.99ρ (pure)、roQn ≥ 0。
    // 実現可能性射影 (Q1/Q2) は定常・RK では各 step ここで、dual-time では sub-iter 内で作動させず物理 step 末尾 (EOS 更新後) の
    // condensationRealizabilityProject_d_wrapper で 1 回だけ (sub-iter 内で毎回射影すると反復値を揺らして収束床を作る; 2026-09-17 run_0263–0269)。
    // 射影は受動種経路 (scheme 1) 限定 (codex result-3 M3; 上の wrapper と同じ理由)。
    const int doProject = (cfg.condRealizProject != 0 && cfg.passiveScalarScheme == 1 && !(cfg.unsteady == 1 && cfg.dualTime == 1)) ? 1 : 0;
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
                var.c_d["T"], var.c_d["P"], g_condTables, cprops,
                var.c_d["condClampCorr_"+i], var.c_d["condClampCorrQ_"+i], condRealizViolCounter(),
                condClampBudget(s), periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr, var.c_d["volume"], doProject,
                condCorrReasons(s));
        } else
        cond_realizability_clamp_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells, var.c_d["ro"], roY_w,
            var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
            cfg.condEvaporation, cfg.condModel, cprops.R, cfg.condEvapRmin, g_rm,
            var.c_d["T"], var.c_d["P"], opts,
            var.c_d["condClampCorr_"+i], var.c_d["condClampCorrQ_"+i], condRealizViolCounter(),
            condClampBudget(s), periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr, var.c_d["volume"], doProject,
            condCorrReasons(s));
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

    if (passiveSchemeEnabled(cfg)) {
        // 受動種経路: 化学種の Dirichlet (ゼロ; node は境界ノードをピン) / Neumann カーネルを受動種ポインタで。
        const int q0 = passive_moment_index0();
        for (size_t k = 0; k < var.condMomentConsNames.size(); ++k) {
            passiveBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var, q0 + (int)k, nullptr);
        }
        return;
    }

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
    if (passiveSchemeEnabled(cfg)) {
        // 受動種経路 (§4.3): S3 面値 (SLAU) または 1 次風上。res_/transport_diag/src_jac のゼロ初期化と順序 (移流→ソース→更新) は同じ。
        passiveAdvection_d_wrapper(cfg, cuda_cfg, msh, var, passive_moment_index0(), (int)var.condMomentConsNames.size());
        // 二相拡散 (#4e): モーメントの残差ゼロ化の後に、化学種・液・Q・エネルギーへ同じ面流束を足す (無効構成は no-op)。
        twoPhaseDiffusion_d_wrapper(cfg, cuda_cfg, msh, var);
        return;
    }

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

    if (passiveSchemeEnabled(cfg)) {
        // 受動種経路 (§4.3): 更新クランプ (θ_u) は不変。増分は passiveImplicitCoupling 1 なら scalar-DPLUR sweep (dq_*_old)、
        // 0 なら point-implicit × passiveImplicitRelax。floor は更新カーネルで掛けず passive_bounds_d が担い補正収支を記録する。
        // 後段の実現可能性クランプ (condensationPrimitive_d_wrapper) は現行のまま。
        const int q0 = passive_moment_index0();
        const int nq = (int)var.condMomentConsNames.size();
        const bool useLimited = (cfg.timeIntegration == 11 && cfg.condLimiterMode == 1 && cfg.condEquilibrium == 0);
        const bool haveDq = passiveDPLURIncrement_d_wrapper(cfg, cuda_cfg, msh, var, q0, nq);
        const flow_float relax = (cfg.timeIntegration == 11) ? static_cast<flow_float>(cfg.passiveImplicitRelax) : static_cast<flow_float>(1.0);
        const flow_float dts   = scalarDtScale(cfg);
        if (useLimited) {
            const int carrier = (cfg.condGasSpecies >= 0 || cfg.condVaporMassFraction > 0.0) ? 1 : 0;
            const CondPropOpts opts = cond_prop_opts(cfg);
            flow_float* roY_w = (carrier && cfg.condGasSpecies >= 0) ? var.c_d["roY" + std::to_string(cfg.condGasSpecies)] : nullptr;
            flow_float* cp_cell   = (cfg.thermalMethod == 2) ? var.c_d["cp"]   : nullptr;
            flow_float* Rmix_cell = (cfg.thermalMethod == 2) ? var.c_d["Rmix"] : nullptr;
            const double lam_min = 0.5;
            for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
                const std::string i = std::to_string(s);
                const std::string g = "rog_"+i, Q2 = "roQ2_"+i, Q1 = "roQ1_"+i, Q0 = "roQ0_"+i;
                cond_moment_update_limited_passive_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                    msh.nCells, var.c_d["dt_local"], var.c_d["volume"], var.c_d["ro"],
                    roY_w, cfg.condVaporMassFraction, var.c_d["T"], cp_cell, Rmix_cell, cfg.cp, cfg.gamma,
                    cfg.condModel, opts, cfg.condDgMaxStep, cfg.condDTmaxStep, lam_min,
                    var.c_d[g+"N"], var.c_d[Q2+"N"], var.c_d[Q1+"N"], var.c_d[Q0+"N"],
                    var.c_d["res_"+g], var.c_d["res_"+Q2], var.c_d["res_"+Q1], var.c_d["res_"+Q0],
                    var.c_d["src_jac_g_"+i], var.c_d["src_jac_Q2_"+i], var.c_d["src_jac_Q1_"+i], var.c_d["src_jac_Q0_"+i],
                    var.c_d["transport_diag_g_"+i], var.c_d["transport_diag_Q2_"+i], var.c_d["transport_diag_Q1_"+i], var.c_d["transport_diag_Q0_"+i],
                    var.c_d[g], var.c_d[Q2], var.c_d[Q1], var.c_d[Q0],
                    var.c_d["condLim_"+i], var.c_d["condClampCorr_"+i], var.c_d["condClampCorrQ_"+i],
                    (double)relax, dts, /*applyFloor=*/0,
                    haveDq ? var.c_d["dq_"+g+"_old"]  : nullptr, haveDq ? var.c_d["dq_"+Q2+"_old"] : nullptr,
                    haveDq ? var.c_d["dq_"+Q1+"_old"] : nullptr, haveDq ? var.c_d["dq_"+Q0+"_old"] : nullptr,
                    /*boundByTheta=*/1,
                    passive_limCorr_cell_ptr(q0+4*s+0), passive_limCorr_cell_ptr(q0+4*s+1), passive_limCorr_cell_ptr(q0+4*s+2), passive_limCorr_cell_ptr(q0+4*s+3),
                    passive_lim_stats_ptr(q0+4*s+0) /* 連続 4 スロットではないので kernel 側は stride 8 で書く */, passive_periodic_root(cfg, msh),
                    // dual-time: dg_max/dT_max の閾値クランプは各物理 step の初回 sub-iter だけ (plan §5.1 #18); 定常は常に (不変)
                    (cfg.unsteady == 1 && cfg.dualTime == 1) ? ((cfg.dualTimeSubIter == 0) ? 1 : 0) : 1);
            }
        } else {
            if (loop == 0) {
                for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
                    const std::string i = std::to_string(s);
                    CHECK_CUDA_ERROR(cudaMemset(var.c_d["condClampCorr_"+i], 0, msh.nCells * sizeof(flow_float)));
                    CHECK_CUDA_ERROR(cudaMemset(var.c_d["condClampCorrQ_"+i], 0, msh.nCells * sizeof(flow_float)));
                }
            }
            for (size_t k = 0; k < var.condMomentConsNames.size(); ++k) {
                const std::string& consName = var.condMomentConsNames[k];
                if (haveDq) {
                    // DPLUR 増分の commit ρφ = ρφ_N + δ (floor なし)。
                    passiveCommitIncrement_d_wrapper(cfg, cuda_cfg, msh, var, q0 + (int)k);
                } else {
                    ScalarTransportDesc desc = buildCondMomentDesc(var, consName);
                    desc.floor = static_cast<flow_float>(-1.0e30);   // floor は passive_bounds_d
                    scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc, relax, dts);
                }
            }
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        const bool rec = passiveRecordStage(cfg, loop);
        passiveAddRhoTerm_d_wrapper(cfg, cuda_cfg, msh, var, q0, nq);              // + φ_N δρ (#19)
        passiveLimitIncrement_d_wrapper(cfg, cuda_cfg, msh, var, q0, nq, rec);   // θ_b 増分スケーリング (M5; 輸送増分 z に対して)
        passiveBounds_d_wrapper(cfg, cuda_cfg, msh, var, q0, nq, rec);           // ρφ >= 0 と補正収支 (最後の砦)
        passiveMirrorPeriodic_d_wrapper(cfg, cuda_cfg, msh, var);
        return;
    }

    if (cfg.timeIntegration == 11 && cfg.condLimiterMode == 1 && cfg.condEquilibrium == 0) {   // 平衡形 (1/2) は従来更新のまま (codex result M5)
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
                var.c_d["condLim_"+i], var.c_d["condClampCorr_"+i], var.c_d["condClampCorrQ_"+i]);
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        return;
    }

    // 旧経路 (mode 0 / 平衡形 / RK): 補正診断は実現可能性クランプ分だけを「このステップの量」として記録する
    // (更新 kernel の floor は未記録 = 旧経路では診断は部分的; codex 2026-09-16 result m1)。ステージ 0 でリセット。
    if (loop == 0) {
        for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
            const std::string i = std::to_string(s);
            CHECK_CUDA_ERROR(cudaMemset(var.c_d["condClampCorr_"+i], 0, msh.nCells * sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemset(var.c_d["condClampCorrQ_"+i], 0, msh.nCells * sizeof(flow_float)));
        }
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

// -----------------------------------------------------------------------------
// 理由別の補正量監視のログ (plans/active/condensation-two-phase-transport.md §4.3, #2)。monitorInterval ごと (main)。
//   1 行/凝縮種: 区間値 (前回ログからの差) と累積 (このプロセスの開始から; restart で 0) を、総液量 Σρg V (root のみ) との比で出す。
//   総液量 0 のときは絶対量を出し、比は「補正が非 0 なら 1」(既存の clampBudget と同じ約束) にする。
//   射影は Q1・Q2 を動かすので、比の分母は ΣρQ1 V・ΣρQ2 V。
//   理由: capViol = 蒸気上限違反 g > Y_w (実現可能性クランプ), negFloor = ρg < 0 → 0 (実現可能性クランプ),
//         incrLimit = 増分制限 (受動種経路の θ_u/θ_b が液 ρg から切った量 = [passive] limCorr の g), passFloor = 受動種の硬い floor (g),
//         proj = モーメント射影, renorm = 化学種再正規化が凝縮種 ρY_w に掛けた補正 Σ|Δ(ρY_w)|V と max|係数−1|,
//         removal = 液滴消滅 (物理)。WARN: removal 以外の区間値が比 1e-6 を超えたとき (閾値は #2 の基準取得で見直す)。
// -----------------------------------------------------------------------------
namespace {
__global__ void cond_liquid_totals_d(geom_int nCells, const flow_float* rog, const flow_float* q1, const flow_float* q2,
                                     const geom_float* vol, const geom_int* root, double* out, const flow_float* q0, const flow_float* rYw)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (root != nullptr && root[ic] != ic) return;
    const double V = (double)vol[ic];
    if (rog[ic] != (flow_float)0.0) atomicAdd(&out[0], (double)rog[ic]*V);
    if (q1[ic]  != (flow_float)0.0) atomicAdd(&out[1], (double)q1[ic]*V);
    if (q2[ic]  != (flow_float)0.0) atomicAdd(&out[2], (double)q2[ic]*V);
    if (q0[ic]  != (flow_float)0.0) atomicAdd(&out[3], (double)q0[ic]*V);   // 再正規化の Q0 補正の分母 (#1b-pre)
    if (rYw != nullptr && rYw[ic] != (flow_float)0.0) atomicAdd(&out[4], (double)rYw[ic]*V);   // 同 ρY_w (自分の総量で割る)
}

// θ の全更新を覆う集計 (#1b-pre (3)): セル配列 theta (condLim_s) の θ<1 のセル数・最小を区間・累積の reason スロットへ (root のみ; 計上だけ)。
__global__ void cond_theta_scan_d(geom_int nCells, const flow_float* theta, const geom_int* root, double* rs, int slotN, int slotMin, int slotLast, int slotCalls)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic == 0) atomicAdd(&rs[slotCalls], 1.0);
    if (ic >= nCells) return;
    if (root != nullptr && root[ic] != ic) return;
    const double t = (double)theta[ic];
    if (t < 1.0 || !(t == t)) {
        atomicAdd(&rs[slotN], 1.0); atomicAdd(&rs[slotLast], 1.0);
        cond_atomic_min_double(&rs[slotMin], (t == t) ? t : -1.0);   // NaN は −1 として残す
    }
}
}  // namespace

void condCorrectionLog_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int iStep)
{
    if (!condensationEnabled(var)) return;
    const int nsp = var.nCondSpeciesRegistered;
    if (condCorrReasons(0) == nullptr) return;
    std::vector<double> cur((size_t)g_condReasons_n*COND_REASON_N, 0.0);
    gpuErrchk( cudaMemcpy(cur.data(), g_condReasons_dev, cur.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    // 受動種経路の g の floor・増分制限 (stats: [2] floor abs 累積, [4] limCorr abs 累積)
    const bool passive = passiveSchemeEnabled(cfg) && passive_moment_index0() >= 0;
    const std::vector<double> pst = passive ? passiveFloorCorrTotals() : std::vector<double>();
    static std::vector<double> s_pst_last;
    if (s_pst_last.size() != pst.size()) s_pst_last.assign(pst.size(), 0.0);
    static double* tot_d = nullptr;
    if (tot_d == nullptr) gpuErrchk( cudaMalloc((void**)&tot_d, 5*sizeof(double)) );
    const geom_int* root = periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr;
    const int nstep = std::max(1, (iStep + 1) - g_condReasons_last_step);
    for (int s = 0; s < nsp && s < g_condReasons_n; ++s) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemset(tot_d, 0, 5*sizeof(double)) );
        flow_float* rYwTot = (cfg.condGasSpecies >= 0 && cfg.thermalMethod == 2 && var.c_d.count("roY" + std::to_string(cfg.condGasSpecies)))
                             ? var.c_d["roY" + std::to_string(cfg.condGasSpecies)] : nullptr;
        cond_liquid_totals_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells, var.c_d["rog_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i], var.c_d["volume"], root, tot_d, var.c_d["roQ0_"+i], rYwTot);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        double tot[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        gpuErrchk( cudaMemcpy(tot, tot_d, 5*sizeof(double), cudaMemcpyDeviceToHost) );
        const double* c = cur.data() + (size_t)s*COND_REASON_N;
        const double* l = g_condReasons_last.data() + (size_t)s*COND_REASON_N;
        auto relz = [](double v, double t) { return (t > 0.0) ? v/t : (v > 0.0 ? 1.0 : 0.0); };
        double pfI = 0.0, pfC = 0.0, plI = 0.0, plC = 0.0;
        if (passive) {
            const size_t q = (size_t)(passive_moment_index0() + 4*s);   // 順序 [g, Q2, Q1, Q0]
            if ((q+1)*8 <= pst.size()) {
                pfC = pst[q*8 + 2]; pfI = pfC - s_pst_last[q*8 + 2];
                plC = pst[q*8 + 4]; plI = plC - s_pst_last[q*8 + 4];
            }
        }
        const double capI = c[COND_REASON_CAP_SUM] - l[COND_REASON_CAP_SUM], negI = c[COND_REASON_NEG_SUM] - l[COND_REASON_NEG_SUM];
        const double rmI  = c[COND_REASON_RM_SUM]  - l[COND_REASON_RM_SUM];
        const double p1I  = c[COND_REASON_PQ1_SUM] - l[COND_REASON_PQ1_SUM], p2I = c[COND_REASON_PQ2_SUM] - l[COND_REASON_PQ2_SUM];
        const double rnI  = c[COND_REASON_RN_SUM]  - l[COND_REASON_RN_SUM];
        const double vmin = c[COND_REASON_VMIN];
        char vbuf[32], lbuf[64], fbuf[64];
        if (passive) {
            std::snprintf(lbuf, sizeof(lbuf), "%.3e (rel %.3e)", plI, relz(plI, tot[0]));
            std::snprintf(fbuf, sizeof(fbuf), "%.3e (rel %.3e)", pfI, relz(pfI, tot[0]));
        } else {
            std::snprintf(lbuf, sizeof(lbuf), "n/a (passiveScalarScheme 0)");
            std::snprintf(fbuf, sizeof(fbuf), "n/a (passiveScalarScheme 0)");
        }
        if (vmin < 1.0e299) std::snprintf(vbuf, sizeof(vbuf), "%.6e", vmin); else std::snprintf(vbuf, sizeof(vbuf), "n/a");
        printf("[cond-corr] step %d species %d interval(%d steps) | total liquid %.9e (Q1 %.6e Q2 %.6e)%s"
               " | capViol %.3e (rel %.3e) nodes %.0f minVaporFrac(pre-limit) %s | negFloor %.3e (rel %.3e) nodes %.0f"
               " | incrLimit %s | passFloor %s | proj Q1 %.3e (rel %.3e) Q2 %.3e (rel %.3e) nodes %.0f"
               " | renorm |dRhoYw| %.3e (rel %.3e) max|f-1| %.3e nodes %.0f | removal(physical) %.3e (rel %.3e) nodes %.0f\n",
               iStep + 1, s, nstep, tot[0], tot[1], tot[2], (tot[0] > 0.0) ? "" : " [liquid 0: values are absolute, rel = 1 if nonzero]",
               capI, relz(capI, tot[0]), c[COND_REASON_CAP_N] - l[COND_REASON_CAP_N], vbuf,
               negI, relz(negI, tot[0]), c[COND_REASON_NEG_N] - l[COND_REASON_NEG_N],
               lbuf, fbuf,
               p1I, relz(p1I, tot[1]), p2I, relz(p2I, tot[2]), c[COND_REASON_P_N] - l[COND_REASON_P_N],
               rnI, relz(rnI, tot[0]), c[COND_REASON_RN_MAX], c[COND_REASON_RN_N] - l[COND_REASON_RN_N],
               rmI, relz(rmI, tot[0]), c[COND_REASON_RM_N] - l[COND_REASON_RM_N]);
        printf("[cond-corr]   species %d cumulative since process start (restart resets to 0): capViol %.9e (rel %.3e) nodes %.0f"
               " | negFloor %.9e (rel %.3e) | incrLimit %.9e | passFloor %.9e | proj Q1 %.9e Q2 %.9e nodes %.0f"
               " | renorm |dRhoYw| %.9e nodes %.0f | removal(physical) %.9e (rel %.3e) nodes %.0f\n",
               s, c[COND_REASON_CAP_SUM], relz(c[COND_REASON_CAP_SUM], tot[0]), c[COND_REASON_CAP_N],
               c[COND_REASON_NEG_SUM], relz(c[COND_REASON_NEG_SUM], tot[0]), plC, pfC,
               c[COND_REASON_PQ1_SUM], c[COND_REASON_PQ2_SUM], c[COND_REASON_P_N],
               c[COND_REASON_RN_SUM], c[COND_REASON_RN_N], c[COND_REASON_RM_SUM], relz(c[COND_REASON_RM_SUM], tot[0]), c[COND_REASON_RM_N]);
        // #1b-pre (2)(3): 再正規化の成分別補正・係数偏差 (区間・累積) と θ・θ_src の全更新を覆う集計
        {
            auto I = [&](int k) { return c[k] - l[k]; };
            // rel は各成分の自分の総量 (ΣρY_w V, Σρg V, ΣρQ_n V) で割る (上の renorm の rel は総液量で割った従来の値)
            printf("[cond-corr]   species %d renorm components [record only: sums since last log / current own total; the gate is [renorm-gate]] | interval: max|f-1| %.3e, |dRhoYw| %.3e (rel %.3e) |dRhog| %.3e (rel %.3e)"
                   " |dQ2| %.3e (rel %.3e) |dQ1| %.3e (rel %.3e) |dQ0| %.3e (rel %.3e)"
                   " | cumulative: max|f-1| %.3e, |dRhoYw| %.9e |dRhog| %.9e |dQ2| %.9e |dQ1| %.9e |dQ0| %.9e\n",
                   s, c[COND_REASON_RN_MAX], rnI, relz(rnI, tot[4]), I(COND_REASON_RNG_SUM), relz(I(COND_REASON_RNG_SUM), tot[0]),
                   I(COND_REASON_RNQ2_SUM), relz(I(COND_REASON_RNQ2_SUM), tot[2]), I(COND_REASON_RNQ1_SUM), relz(I(COND_REASON_RNQ1_SUM), tot[1]),
                   I(COND_REASON_RNQ0_SUM), relz(I(COND_REASON_RNQ0_SUM), tot[3]),
                   c[COND_REASON_RN_MAXC], c[COND_REASON_RN_SUM], c[COND_REASON_RNG_SUM], c[COND_REASON_RNQ2_SUM], c[COND_REASON_RNQ1_SUM], c[COND_REASON_RNQ0_SUM]);
            printf("[cond-corr]   species %d theta over all updates | interval: update theta<1 %.0f cell-updates in %.0f updates (last update %.0f cells) min %.6g"
                   " | source theta_src<1 %.0f cell-evals in %.0f evals (last eval %.0f cells) min %.6g"
                   " | cumulative: update theta<1 %.0f in %.0f updates, theta_src<1 %.0f in %.0f evals\n",
                   s, I(COND_REASON_TU_N), I(COND_REASON_TU_CALLS), c[COND_REASON_TU_LAST], c[COND_REASON_TU_MIN],
                   I(COND_REASON_TS_N), I(COND_REASON_TS_CALLS), c[COND_REASON_TS_LAST], c[COND_REASON_TS_MIN],
                   c[COND_REASON_TU_N], c[COND_REASON_TU_CALLS], c[COND_REASON_TS_N], c[COND_REASON_TS_CALLS]);
        }
        // WARN: 液滴消滅以外の数値補正の区間値が総液量比 1e-6 超
        const double thr = 1.0e-6;
        struct { const char* name; double rel; } chk[] = {
            {"capViol", relz(capI, tot[0])}, {"negFloor", relz(negI, tot[0])},
            {"incrLimit", passive ? relz(plI, tot[0]) : 0.0}, {"passFloor", passive ? relz(pfI, tot[0]) : 0.0},
            {"projQ1", relz(p1I, tot[1])}, {"projQ2", relz(p2I, tot[2])}, {"renorm", relz(rnI, tot[0])}};
        for (const auto& k : chk)
            if (k.rel > thr)
                printf("[cond-corr] WARN step %d species %d: %s interval correction %.3e of total > %.0e\n", iStep + 1, s, k.name, k.rel, thr);
    }
    twoPhaseUpdateLog(cfg, iStep);   // 二相拡散の非分割更新の監視 (#4e; 無効構成は no-op)
    renormGateLog(cfg, iStep, false);   // 再正規化の受入ゲート: 前回ログからの全更新の max (#1b-pre)
    g_condReasons_last = cur;
    if (passive) s_pst_last = pst;
    g_condReasons_last_step = iStep + 1;
    condReasonsResetInterval(g_condReasons_n);
    fflush(stdout);
}

// -----------------------------------------------------------------------------
// 二相拡散の起動時検査・蒸気残差・非分割更新 (plans/active/condensation-two-phase-transport.md §4.2, §5.1 #4e)
// -----------------------------------------------------------------------------
void condTwoPhaseDiffusionValidate(const solverConfig& cfg)
{
    if (cfg.condAuditResidual == 1 && cfg.condTwoPhaseDiffusion == 0) {
        if (condResidualAuditActive(cfg))
            std::printf("[twophase-audit] condAuditResidual 1: the existing operator (species Fick diffusion; liquid and moments advected only) is re-evaluated in double at start and end\n");
        else
            std::printf("[twophase-audit] condAuditResidual 1 is inactive: needs TP carrier condensation (condGasSpecies >= 0, thermalMethod 2, nCondSpecies 1)\n");
    }
    if (cfg.condTwoPhaseDiffusion == 0) return;
    auto fail = [](const std::string& why) {
        std::fprintf(stderr, "Configuration Error: condensation.condTwoPhaseDiffusion 1 %s (plan condensation-two-phase-transport #4e: steady-only first version)\n", why.c_str());
        std::exit(EXIT_FAILURE);
    };
    if (cfg.condensation != 1) { std::printf("[twophase] condTwoPhaseDiffusion 1 is inactive: condensation is off\n"); return; }
    if (cfg.condGasSpecies < 0 || cfg.thermalMethod != 2 || cfg.nSpecies < 2) {
        // CPG carrier (condVaporMassFraction) と pure 凝縮は対象外 (§4.2; 液は拡散しない現行のまま)
        std::printf("[twophase] condTwoPhaseDiffusion 1 is inactive: %s (two-phase diffusion is not supported; liquid is not diffused)\n",
                    (cfg.condVaporMassFraction > 0.0) ? "CPG carrier" : "not a TP carrier (pure condensible)");
        return;
    }
    if (cfg.viscMethod == 0) { std::printf("[twophase] condTwoPhaseDiffusion 1 is inactive: inviscid (viscMethod 0)\n"); return; }
    if (cfg.unsteady == 1 && cfg.dualTime == 1) fail("cannot be combined with dual-time (unsteady 1, dualTime 1): the coupled FCT of vapour and liquid is not designed yet");
    if (cfg.unsteady != 0 || cfg.timeIntegration != 11) fail("requires steady implicit pseudo-time (unsteady 0, timeIntegration 11)");
    if (cfg.speciesImplicitCoupling == 2) fail("cannot be combined with speciesImplicitCoupling 2 (EOS cross coupling commits water before the liquid)");
    if (cfg.passiveScalarScheme != 1) fail("requires passiveScalarScheme 1 (moments on the passive-scalar path)");
    if (cfg.condEquilibrium != 0) fail("requires non-equilibrium condensation (condEquilibrium 0)");
    if (cfg.condLimiterMode != 1) fail("requires condLimiterMode 1 (the source enters the residual without theta)");
    if (cfg.nCondSpecies != 1) fail("supports one condensing species (nCondSpecies 1)");
    std::printf("[twophase] condTwoPhaseDiffusion 1: gas-phase molecular diffusion (z basis, upwind correction) + common turbulent mixing "
                "(species, liquid g, Q2/Q1/Q0, energy), unsplit vapour/liquid update with point-diagonal preconditioner, relax %.3g, "
                "condDgMaxStep %.3g, condDTmaxStep %.3g; residual column rms_roYv = res_roY%d - res_rog_0\n",
                cfg.condTwoPhaseRelax, cfg.condDgMaxStep, cfg.condDTmaxStep, cfg.condGasSpecies);
    if (cfg.discretization != "node")
        std::printf("[twophase] WARNING: cell discretization is unverified for two-phase diffusion (only the code path was reviewed)\n");
    if (cfg.passiveImplicitCoupling == 1 || cfg.speciesImplicitCoupling == 1)
        std::printf("[twophase] note: vapour, liquid and moments use the point-diagonal preconditioner (passiveImplicitCoupling / speciesImplicitCoupling "
                    "DPLUR is not used for them); the other species keep their update\n");
}

namespace {
__global__ void twophase_vapor_residual_d(geom_int nCells, const flow_float* res_w, const flow_float* res_g, flow_float* res_v)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) res_v[ic] = res_w[ic] - res_g[ic];
}

// stats (double[8], root のみ・ログ間隔の区間値): 0 Σ(1−θ)|δρv|V, 1 Σ(1−θ)|δρg|V, 2 θ<1 のセル数, 3 θ=0 のセル数,
//   4 Σ Q_cut V, 5 Σ v_round V, 6 更新したセル数。thetaMin: θ の最小 ×1e9 (atomicMin)。
__global__ void twophase_vl_update_d(
    geom_int nCells, flow_float* dt_local, flow_float dtScale, geom_float* vol, flow_float* ro, flow_float* T,
    flow_float* cp_cell, flow_float* Rmix_cell, int condModel, CondPropOpts opts, double dg_max, double dT_max, flow_float omega,
    flow_float* roYw, flow_float* res_w, flow_float* td_w,
    flow_float* g, flow_float* Q2, flow_float* Q1, flow_float* Q0,
    flow_float* N_g, flow_float* N_Q2, flow_float* N_Q1, flow_float* N_Q0,
    flow_float* res_g, flow_float* res_Q2, flow_float* res_Q1, flow_float* res_Q0,
    flow_float* sj_g, flow_float* sj_Q2, flow_float* sj_Q1, flow_float* sj_Q0,
    flow_float* td_g, flow_float* td_Q2, flow_float* td_Q1, flow_float* td_Q0,
    flow_float* diagLim, flow_float* diagCorrG, flow_float* diagCorrQ,
    double* stats, int* thetaMin, double* limStats_g, const geom_int* root,
    double* diag)   // 診断 (condTwoPhaseDiag; nullptr で書かない): TPD_* 列 × nCells
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double v = (double)vol[ic];
    const double dt = (double)(dt_local[ic]*dtScale);
    TpCellIn c;
    c.M = (float)(v/fmax(dt, 1.0e-30)); c.V = (float)v;
    c.Rw = res_w[ic]; c.Rg = res_g[ic]; c.RQ[0] = res_Q2[ic]; c.RQ[1] = res_Q1[ic]; c.RQ[2] = res_Q0[ic];
    c.Dv = td_w[ic]; c.Dg = td_g[ic]; c.DQ[0] = td_Q2[ic]; c.DQ[1] = td_Q1[ic]; c.DQ[2] = td_Q0[ic];
    c.sjg = sj_g[ic]; c.sjQ[0] = sj_Q2[ic]; c.sjQ[1] = sj_Q1[ic]; c.sjQ[2] = sj_Q0[ic];
    c.rYw = roYw[ic]; c.rg = N_g[ic]; c.rQ[0] = N_Q2[ic]; c.rQ[1] = N_Q1[ic]; c.rQ[2] = N_Q0[ic];
    c.rho = ro[ic]; c.omega = omega; c.dg_max = dg_max; c.dT_max = dT_max;
    // θ の ΔT 換算: 既存の更新クランプ (condensationUpdateLimiter_d.cuh) と同じ有効比熱
    {
        const CondSpeciesProps cprops = condProps_make(condModel, opts);
        const double rod = (double)ro[ic];
        const double g_old = (rod > 1.0e-20) ? fmax((double)N_g[ic]/rod, 0.0) : 0.0;
        const double Td = (double)T[ic];
        const double cvg = fmax((double)cp_cell[ic] - (double)Rmix_cell[ic], 1.0e-3);
        const double L  = cond_latent(cprops, Td);
        const double dL = (cond_latent(cprops, Td + 0.1) - cond_latent(cprops, Td - 0.1))/0.2;
        c.L = L; c.cveff = fmax(cvg + g_old*(cprops.R - dL), 1.0e-2*cvg);
    }
    TpCellOut o;
    tp_vl_update(c, o);
    if (diag != nullptr) {   // 更新前の格納値・制限前増分・θ を決めた制限・Q の残差 (読むだけ)
        const size_t n = (size_t)nCells;
        diag[0*n + ic] = (double)c.rYw - (double)c.rg; diag[1*n + ic] = (double)c.rg;
        diag[2*n + ic] = o.dv; diag[3*n + ic] = o.dg;
        for (int m = 0; m < TP_NQ; ++m) { diag[(4 + m)*n + ic] = o.dq[m]; diag[(8 + m)*n + ic] = c.RQ[m]; }
        diag[7*n + ic] = (double)o.reason; diag[11*n + ic] = o.theta;
        // #1b-r2: 更新に入った float の全残差 (組立 A) と前処理の分母
        diag[15*n + ic] = c.Rw; diag[16*n + ic] = c.Rg;
        diag[17*n + ic] = (double)(c.M + c.Dv); diag[18*n + ic] = (double)(c.M + c.Dg + c.V*c.sjg);
        for (int m = 0; m < TP_NQ; ++m) diag[(19 + m)*n + ic] = (double)(c.M + c.DQ[m] + c.V*c.sjQ[m]);
    }
    roYw[ic] = o.rYw; g[ic] = o.rg; Q2[ic] = o.rQ[0]; Q1[ic] = o.rQ[1]; Q0[ic] = o.rQ[2];
    diagLim[ic] = (flow_float)o.theta;
    const double rod = (double)ro[ic];
    diagCorrG[ic] = (flow_float)((rod > 1.0e-20) ? o.vround/rod : 0.0);   // 状態補正 (蒸気の丸め) [質量分率]; 後段の実現可能性クランプが同じ配列へ累積
    diagCorrQ[ic] = (flow_float)((o.qcut > 0.0) ? 1.0 : 0.0);
    if (root == nullptr || root[ic] == ic) {
        if (o.theta < 1.0) {
            atomicAdd(&stats[0], o.withheld_v*v); atomicAdd(&stats[1], o.withheld_g*v); atomicAdd(&stats[2], 1.0);
            if (o.theta == 0.0) atomicAdd(&stats[3], 1.0);
            atomicMin(thetaMin, (int)(o.theta*1.0e9));
            if (limStats_g != nullptr) { atomicAdd(&limStats_g[0], o.withheld_g*v); atomicAdd(&limStats_g[1], 1.0); }
        }
        if (o.qcut > 0.0) atomicAdd(&stats[4], o.qcut*v);
        if (o.vround > 0.0) atomicAdd(&stats[5], o.vround*v);
        atomicAdd(&stats[6], 1.0);
    }
}
}  // namespace

static double* g_tp_stats_dev = nullptr;
static int*    g_tp_thetaMin_dev = nullptr;

// ---- 診断 (condTwoPhaseDiag, #1b-r1; 読むだけ) ----
// セル列 (double × nCells): 0 ρv (更新前), 1 ρg (更新前), 2 δρv, 3 δρg, 4–6 δρQ2/Q1/Q0 (制限前・緩和後), 7 制限 (0 なし, 1 蒸気非負, 2 液非負, 3 dg_max, 4 dT_max),
//   8–10 Q2/Q1/Q0 の残差 (更新に入った全残差), 11 θ, 12 f−1 (再正規化), 13 Δ(ρY_w) (再正規化), 14 Δ(ρg) (再正規化),
//   15 R_w (float 組立 A), 16 R_g (A), 17 前処理の分母 蒸気 V/Δτ + D_v, 18 液 V/Δτ + D_g + V sj_g, 19–21 Q2/Q1/Q0 の分母 (#1b-r2)
#define TPD_NCOL 22
#define TPD_REC  (2 + TPD_NCOL)   // 記録: 更新番号, セル, 上の 15 列
#define TPD_WINDOW 200
static double* g_tpd_cell = nullptr;    // TPD_NCOL × nCells
static double* g_tpd_rec = nullptr;     // TPD_REC × cap
static unsigned long long* g_tpd_n = nullptr;
static int*    g_tpd_hit = nullptr;     // [0,n) θ<1 / [n,2n) θ=0 が区間内にあったセル (固有セル数)
static int*    g_tpd_cnt = nullptr;     // 区間の固有セル数 (2)
static size_t  g_tpd_cap = 0, g_tpd_nupd = 0, g_tpd_nupd_log = 0, g_tpd_ncells = 0;
static bool tpdOn(const solverConfig& cfg) { return cfg.condTwoPhaseDiag != 0 && condTwoPhaseDiffusionActive(cfg); }
// 診断の窓 (末尾 TPD_WINDOW 更新) に入っているか: 次に行う更新の番号 g_tpd_nupd で判定 (組立 B は更新の前に呼ぶので同じ番号)
bool twoPhaseDiagInWindow(const solverConfig& cfg)
{
    const size_t N = (size_t)std::max(0, cfg.mainLoopCount());
    return tpdOn(cfg) && (g_tpd_nupd + TPD_WINDOW >= N);
}
double* twoPhaseDiagRenormPtr() { return (g_tpd_cell != nullptr) ? g_tpd_cell + 12*g_tpd_ncells : nullptr; }
static void tpdAlloc(geom_int nCells)
{
    if (g_tpd_cell != nullptr) return;
    g_tpd_ncells = (size_t)nCells;
    g_tpd_cap = 400000;   // 記録の上限 (超えた分は数えるだけ; 終了時に報告)
    gpuErrchk( cudaMalloc((void**)&g_tpd_cell, TPD_NCOL*g_tpd_ncells*sizeof(double)) );
    gpuErrchk( cudaMemset(g_tpd_cell, 0, TPD_NCOL*g_tpd_ncells*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&g_tpd_rec, TPD_REC*g_tpd_cap*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&g_tpd_n, sizeof(unsigned long long)) ); gpuErrchk( cudaMemset(g_tpd_n, 0, sizeof(unsigned long long)) );
    gpuErrchk( cudaMalloc((void**)&g_tpd_hit, 2*g_tpd_ncells*sizeof(int)) ); gpuErrchk( cudaMemset(g_tpd_hit, 0, 2*g_tpd_ncells*sizeof(int)) );
    gpuErrchk( cudaMalloc((void**)&g_tpd_cnt, 2*sizeof(int)) );
}
namespace {
__global__ void tpd_collect_d(geom_int nCells, const double* cell, const geom_int* root, int mode, int inWindow, double upd,
                              double* rec, unsigned long long* nrec, unsigned long long cap, int* hit)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (root != nullptr && root[ic] != ic) return;
    const size_t n = (size_t)nCells;
    const double th = cell[11*n + ic];
    if (th < 1.0) hit[ic] = 1;
    if (th == 0.0) hit[n + ic] = 1;
    const bool sel = (mode == 2) ? (th < 1.0) : (th == 0.0);   // 1, 3: θ = 0、2: θ < 1
    if (!inWindow || !sel) return;
    const unsigned long long k = atomicAdd(nrec, 1ULL);
    if (k >= cap) return;
    double* r = rec + (size_t)k*TPD_REC;
    r[0] = upd; r[1] = (double)ic;
    for (int j = 0; j < TPD_NCOL; ++j) r[2 + j] = cell[(size_t)j*n + ic];
}
__global__ void tpd_count_d(geom_int nCells, int* hit, int* cnt)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (hit[ic]) { atomicAdd(&cnt[0], 1); hit[ic] = 0; }
    if (hit[nCells + ic]) { atomicAdd(&cnt[1], 1); hit[nCells + ic] = 0; }
}
}  // namespace

// ---- #1b-r2 (condTwoPhaseDiag 3): 組立 A (float, 更新に使用) と B (double, 診断) の比較 (host; 窓内の各更新) ----
struct TpbRow {   // 記録したセル 1 つ (更新ごと)
    long upd, cell; double x, y; int group; int inStop; double th; int reason;
    double rv, rg, RgA, RgB, RgAdv, RgDif, RgSrc, RgScaleLocal, RgErr, denG, dRg, RvA, RvB, denV, RwA, RwB;
    double RQA[3], RQB[3], denQ[3], tol[4]; int cls;
};
struct TpbUpd { long upd; long nStop; double maxIn[3], maxOut[3], tol[3]; long argmax[3]; int argInStop[3]; int reject; };
static std::vector<TpbRow> g_tpb_rows;
static std::vector<TpbUpd> g_tpb_upd;
static bool g_tpb_skip_noted = false;
static const char* TPB_GROUP[] = {"stop_dry_theta0", "max_w", "max_v", "max_g", "max_Q2", "max_Q1", "max_Q0"};
static const char* TPB_CLASS[] = {"n/a", "float_rounding_trigger", "discrete_term_negative", "indeterminate", "not_liquid_negative"};
static void tpdBProcess(solverConfig& cfg, mesh& msh, variables& var, size_t upd)
{
    const double *rBd, *ABd, *brkd; int NQ = 0;
    if (!twoPhaseDiagBGet(&rBd, &ABd, &brkd, &NQ)) {
        if (!g_tpb_skip_noted) { printf("[twophase-diag] condTwoPhaseDiag 3: assembly B unavailable (needs S3 face values with SLAU, non-periodic) — A/B rows skipped\n"); g_tpb_skip_noted = true; }
        return;
    }
    const size_t nc = (size_t)msh.nCells; const int ns = var.nSpeciesRegistered, iw = cfg.condGasSpecies;
    static std::vector<double> rB, AB, brk, cell; static std::vector<flow_float> pin, cx, cy;
    rB.resize(nc*NQ); AB.resize(nc*NQ); brk.resize(nc*3); cell.resize(TPD_NCOL*nc);
    gpuErrchk( cudaMemcpy(rB.data(), rBd, nc*NQ*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(AB.data(), ABd, nc*NQ*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(brk.data(), brkd, nc*3*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(cell.data(), g_tpd_cell, TPD_NCOL*nc*sizeof(double), cudaMemcpyDeviceToHost) );
    if (pin.size() != nc) {
        cx.resize(nc); cy.resize(nc);
        gpuErrchk( cudaMemcpy(cx.data(), var.c_d["ccx"], nc*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(cy.data(), var.c_d["ccy"], nc*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        pin.assign(nc, 0.0f);
        auto itp = var.c_d.find("scalarDirichletPin");
        if (cfg.discretization == "node" && itp != var.c_d.end() && itp->second != nullptr)
            gpuErrchk( cudaMemcpy(pin.data(), itp->second, nc*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    }
    auto C = [&](int col, size_t c) { return cell[(size_t)col*nc + c]; };
    const int qc[6] = {iw, ns, ns + 1, ns + 2, ns + 3, ns + 4};   // w, v, g, Q2, Q1, Q0
    const double eps = 1.1920928955078125e-7;
    const std::vector<double>& r0 = twoPhaseAuditR0();
    double tol[6]; long amax[6]; double vmax[6];
    for (int k = 0; k < 6; ++k) {
        double mA = 0.0; amax[k] = -1; vmax[k] = -1.0;
        for (size_t c = 0; c < nc; ++c) {
            if (pin[c] == 1.0f) continue;
            const double a = AB[c*NQ + qc[k]], v = std::fabs(rB[c*NQ + qc[k]]);
            if (a > mA) mA = a;
            if (v > vmax[k] || !(v == v)) { vmax[k] = v; amax[k] = (long)c; }
        }
        const double rr = (r0.size() == (size_t)NQ && r0[qc[k]] >= 0.0) ? 1.0e-7*r0[qc[k]] : 0.0;
        tol[k] = std::max(rr, 6.0*eps*mA);   // 監査と同じ許容 (変更しない)
    }
    std::vector<char> inS(nc, 0); long nS = 0;
    // 乾燥停止: θ = 0 かつ更新前 ρg = 0。試験用 FORGE_TPD3_STOPSET=lt1 は θ < 1 の全セルを停止集合にする (分類の経路を短い run で通すためだけ)
    static const bool lt1 = (std::getenv("FORGE_TPD3_STOPSET") != nullptr && std::string(std::getenv("FORGE_TPD3_STOPSET")) == "lt1");
    for (size_t c = 0; c < nc; ++c) if (pin[c] != 1.0f && (lt1 ? (C(11, c) < 1.0) : (C(11, c) == 0.0 && C(1, c) == 0.0))) { inS[c] = 1; ++nS; }
    TpbUpd u{}; u.upd = (long)upd; u.nStop = nS; u.reject = 1;
    for (int m = 0; m < 3; ++m) {
        const int q = qc[3 + m];
        double mi = 0.0, mo = 0.0;
        for (size_t c = 0; c < nc; ++c) { if (pin[c] == 1.0f) continue; const double v = std::fabs(rB[c*NQ + q]); if (inS[c]) mi = std::max(mi, v); else mo = std::max(mo, v); }
        u.maxIn[m] = mi; u.maxOut[m] = mo; u.tol[m] = tol[3 + m]; u.argmax[m] = amax[3 + m];
        u.argInStop[m] = (amax[3 + m] >= 0 && inS[(size_t)amax[3 + m]]) ? 1 : 0;
        if (u.argInStop[m] || !(mi <= tol[3 + m])) u.reject = 0;   // 直接支配説の棄却: 最大セルが集合外 かつ 集合内の Q 残差が許容以下 (全 Q)
    }
    g_tpb_upd.push_back(u);
    auto addRow = [&](size_t c, int group) {
        TpbRow r{}; r.upd = (long)upd; r.cell = (long)c; r.x = cx[c]; r.y = cy[c]; r.group = group; r.inStop = inS[c];
        r.th = C(11, c); r.reason = (int)C(7, c); r.rv = C(0, c); r.rg = C(1, c);
        r.RgA = C(16, c); r.RgB = rB[c*NQ + ns + 1]; r.RgAdv = brk[c*3 + 0]; r.RgDif = brk[c*3 + 1]; r.RgSrc = brk[c*3 + 2];
        r.RgScaleLocal = AB[c*NQ + ns + 1]; r.RgErr = 6.0*eps*r.RgScaleLocal;
        r.denG = C(18, c); r.dRg = C(3, c); r.RvA = (double)(float)((float)C(15, c) - (float)C(16, c)); r.RvB = rB[c*NQ + ns]; r.denV = C(17, c);
        r.RwA = C(15, c); r.RwB = rB[c*NQ + iw];
        for (int m = 0; m < 3; ++m) { r.RQA[m] = C(8 + m, c); r.RQB[m] = rB[c*NQ + ns + 2 + m]; r.denQ[m] = C(19 + m, c); }
        r.tol[0] = tol[2]; r.tol[1] = tol[3]; r.tol[2] = tol[4]; r.tol[3] = tol[5];
        if (group == 0) {   // 乾燥停止セルの判別 (plan #1b-r2 の規則を機械化; 局所演算誤差の尺度は監査と同じ 6ε × セルの Σ|項|)
            if (!(r.RgA < 0.0)) r.cls = 4;
            else if (r.RgB >= 0.0) r.cls = 1;                                   // A で負・B で非負
            else if (std::fabs(r.RgB) > r.RgErr) r.cls = 2;                      // B でも負で局所演算誤差を超える
            else r.cls = 3;                                                      // 誤差範囲で符号不定
        }
        g_tpb_rows.push_back(r);
    };
    for (size_t c = 0; c < nc; ++c) if (inS[c]) addRow(c, 0);
    for (int k = 0; k < 6; ++k) if (amax[k] >= 0) addRow((size_t)amax[k], 1 + k);
}

static void tpdBWrite(solverConfig& cfg)
{
    (void)cfg;
    if (g_tpb_upd.empty()) return;
    if (FILE* f = std::fopen("twophase_diag3_cells.csv", "w")) {
        std::fprintf(f, "update,cell,x,y,group,in_stop_set,theta,limit_reason,rv_pre,rg_pre,Rg_A_float,Rg_B_double,Rg_B_advection,Rg_B_diffusion,Rg_source,"
                        "Rg_local_scale,Rg_local_error_bound,precond_denom_g,d_rg_prelimit,Rv_A_float,Rv_B_double,precond_denom_v,Rw_A_float,Rw_B_double,"
                        "RQ2_A_float,RQ2_B_double,RQ1_A_float,RQ1_B_double,RQ0_A_float,RQ0_B_double,precond_denom_Q2,precond_denom_Q1,precond_denom_Q0,"
                        "tol_g,tol_Q2,tol_Q1,tol_Q0,stop_class,Rg_AB_gap,Rg_AB_gap_over_error_bound\n");
        for (const TpbRow& r : g_tpb_rows)
            std::fprintf(f, "%ld,%ld,%.9g,%.9g,%s,%d,%.9e,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
                            "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%s,%.9e,%.6e\n",
                         r.upd, r.cell, r.x, r.y, TPB_GROUP[r.group], r.inStop, r.th, r.reason, r.rv, r.rg, r.RgA, r.RgB, r.RgAdv, r.RgDif, r.RgSrc,
                         r.RgScaleLocal, r.RgErr, r.denG, r.dRg, r.RvA, r.RvB, r.denV, r.RwA, r.RwB,
                         r.RQA[0], r.RQB[0], r.RQA[1], r.RQB[1], r.RQA[2], r.RQB[2], r.denQ[0], r.denQ[1], r.denQ[2],
                         r.tol[0], r.tol[1], r.tol[2], r.tol[3], (r.group == 0) ? TPB_CLASS[r.cls] : "n/a",
                         std::fabs(r.RgA - r.RgB), (r.RgErr > 0.0) ? std::fabs(r.RgA - r.RgB)/r.RgErr : (r.RgA == r.RgB ? 0.0 : INFINITY));
        std::fclose(f);
    }
    if (FILE* f = std::fopen("twophase_diag3_qcompare.csv", "w")) {
        std::fprintf(f, "update,n_stop,maxin_Q2,maxout_Q2,tol_Q2,maxin_over_tol_Q2,maxout_over_tol_Q2,argmax_cell_Q2,argmax_in_stop_Q2,"
                        "maxin_Q1,maxout_Q1,tol_Q1,maxin_over_tol_Q1,maxout_over_tol_Q1,argmax_cell_Q1,argmax_in_stop_Q1,"
                        "maxin_Q0,maxout_Q0,tol_Q0,maxin_over_tol_Q0,maxout_over_tol_Q0,argmax_cell_Q0,argmax_in_stop_Q0,direct_dominance_rejected\n");
        for (const TpbUpd& u : g_tpb_upd) {
            std::fprintf(f, "%ld,%ld", u.upd, u.nStop);
            for (int m = 0; m < 3; ++m)
                std::fprintf(f, ",%.9e,%.9e,%.9e,%.6e,%.6e,%ld,%d", u.maxIn[m], u.maxOut[m], u.tol[m],
                             (u.tol[m] > 0.0) ? u.maxIn[m]/u.tol[m] : NAN, (u.tol[m] > 0.0) ? u.maxOut[m]/u.tol[m] : NAN, u.argmax[m], u.argInStop[m]);
            std::fprintf(f, ",%d\n", u.reject);
        }
        std::fclose(f);
    }
    // 判定の要約 (plan #1b-r2 の規則を機械的に適用; 規則の操作化は設計メモ §19)
    long cls[5] = {0, 0, 0, 0, 0}, nStopRows = 0, dom[3] = {0, 0, 0}, rej = 0, gapExceed = 0;   // gapExceed: |A−B| が局所の誤差尺度を超えた停止行 (float の非正規数域の項など)
    double minRgB = 0.0;
    for (const TpbRow& r : g_tpb_rows) {
        if (r.group != 0) continue;
        ++nStopRows; ++cls[r.cls];
        if (std::fabs(r.RgA - r.RgB) > r.RgErr) ++gapExceed;
        if (r.cls == 2) {   // B でも負: 最も負の項
            const double t[3] = {r.RgAdv, r.RgDif, r.RgSrc}; int k = 0; for (int j = 1; j < 3; ++j) if (t[j] < t[k]) k = j; ++dom[k];
            minRgB = std::min(minRgB, r.RgB);
        }
    }
    for (const TpbUpd& u : g_tpb_upd) rej += u.reject;
    const char* vs;
    if (nStopRows == 0) vs = "no_dry_stop_cells";
    else if (cls[2] > 0) vs = "discrete_terms_create_the_negative";
    else if (cls[1] > 0 && cls[3] == 0) vs = "float_assembly_rounding_triggers_the_stop";
    else if (cls[1] > 0) vs = "mixed_rounding_and_indeterminate";
    else vs = "indeterminate";
    const long W = (long)g_tpb_upd.size();
    const char* vd = (rej == W) ? "direct_dominance_rejected_in_all_updates" : (rej == 0 ? "direct_dominance_not_rejected_in_any_update" : "direct_dominance_rejected_in_some_updates");
    if (FILE* f = std::fopen("twophase_diag3_summary.csv", "w")) {
        std::fprintf(f, "key,value\nupdates,%ld\nfirst_update,%ld\nlast_update,%ld\nstop_cell_updates,%ld\n", W, g_tpb_upd.front().upd, g_tpb_upd.back().upd, nStopRows);
        std::fprintf(f, "class_float_rounding_trigger,%ld\nclass_discrete_term_negative,%ld\nclass_indeterminate,%ld\nclass_not_liquid_negative,%ld\n", cls[1], cls[2], cls[3], cls[4]);
        std::fprintf(f, "discrete_negative_most_negative_term_advection,%ld\ndiscrete_negative_most_negative_term_diffusion,%ld\ndiscrete_negative_most_negative_term_source,%ld\n", dom[0], dom[1], dom[2]);
        std::fprintf(f, "min_Rg_B_among_discrete_negative,%.9e\nverdict_sign_rule,%s\n", minRgB, vs);
        std::fprintf(f, "stop_rows_AB_gap_exceeds_local_error_bound,%ld\n", gapExceed);
        std::fprintf(f, "updates_direct_dominance_rejected,%ld\nverdict_direct_dominance,%s\n", rej, vd);
        std::fclose(f);
    }
    printf("[twophase-diag] wrote twophase_diag3_cells.csv, twophase_diag3_qcompare.csv, twophase_diag3_summary.csv (%ld updates, %ld stop cell-updates;"
           " classes rounding %ld / discrete %ld / indeterminate %ld; %s; direct dominance rejected in %ld of %ld updates)\n",
           W, nStopRows, cls[1], cls[2], cls[3], vs, rej, W);
    fflush(stdout);
}
static void twoPhaseStatsAlloc()
{
    if (g_tp_stats_dev != nullptr) return;
    gpuErrchk( cudaMalloc((void**)&g_tp_stats_dev, 8*sizeof(double)) );
    gpuErrchk( cudaMemset(g_tp_stats_dev, 0, 8*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&g_tp_thetaMin_dev, sizeof(int)) );
    const int one = 1000000000;
    gpuErrchk( cudaMemcpy(g_tp_thetaMin_dev, &one, sizeof(int), cudaMemcpyHostToDevice) );
}

void twoPhaseVaporResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condTwoPhaseDiffusionActive(cfg) || !condensationEnabled(var)) return;
    twophase_vapor_residual_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.c_d["res_roY" + std::to_string(cfg.condGasSpecies)], var.c_d["res_rog_0"], var.c_d["res_roYv"]);
    gpuErrchk( cudaPeekAtLastError() );
}

void twoPhaseHoldWater_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cuda_cfg;
    if (!condTwoPhaseDiffusionActive(cfg) || !condensationEnabled(var)) return;
    const std::string w = "roY" + std::to_string(cfg.condGasSpecies);
    gpuErrchk( cudaMemcpy(var.c_d[w], var.c_d[w+"N"], (size_t)msh.nCells_all*sizeof(flow_float), cudaMemcpyDeviceToDevice) );
}

void twoPhaseUpdate_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condTwoPhaseDiffusionActive(cfg) || !condensationEnabled(var)) return;
    twoPhaseStatsAlloc();
    if (tpdOn(cfg)) tpdAlloc(msh.nCells);
    const std::string w = "roY" + std::to_string(cfg.condGasSpecies);
    const std::string g = "rog_0", Q2 = "roQ2_0", Q1 = "roQ1_0", Q0 = "roQ0_0";
    const int q0 = passive_moment_index0();
    twophase_vl_update_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.c_d["dt_local"], scalarDtScale(cfg), var.c_d["volume"], var.c_d["ro"], var.c_d["T"],
        var.c_d["cp"], var.c_d["Rmix"], cfg.condModel, cond_prop_opts(cfg), cfg.condDgMaxStep, cfg.condDTmaxStep, (flow_float)cfg.condTwoPhaseRelax,
        var.c_d[w], var.c_d["res_"+w], var.c_d["transport_diag_Y" + std::to_string(cfg.condGasSpecies)],
        var.c_d[g], var.c_d[Q2], var.c_d[Q1], var.c_d[Q0],
        var.c_d[g+"N"], var.c_d[Q2+"N"], var.c_d[Q1+"N"], var.c_d[Q0+"N"],
        var.c_d["res_"+g], var.c_d["res_"+Q2], var.c_d["res_"+Q1], var.c_d["res_"+Q0],
        var.c_d["src_jac_g_0"], var.c_d["src_jac_Q2_0"], var.c_d["src_jac_Q1_0"], var.c_d["src_jac_Q0_0"],
        var.c_d["transport_diag_g_0"], var.c_d["transport_diag_Q2_0"], var.c_d["transport_diag_Q1_0"], var.c_d["transport_diag_Q0_0"],
        var.c_d["condLim_0"], var.c_d["condClampCorr_0"], var.c_d["condClampCorrQ_0"],
        g_tp_stats_dev, g_tp_thetaMin_dev, passive_lim_stats_ptr(q0), passive_periodic_root(cfg, msh),
        g_tpd_cell);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    // 再正規化 (係数を液・Q にも) → 受動種の最後の砦 (floor と収支の記録; 通常は無作用) → 周期ミラー
    speciesRenormalizeTwoPhase_d_wrapper(cfg, cuda_cfg, msh, var);
    if (g_tpd_cell != nullptr) {   // 診断: 末尾 200 更新の θ = 0 (θ < 1) セルを記録、区間の固有セルの印 (読むだけ)
        const int inWindow = twoPhaseDiagInWindow(cfg) ? 1 : 0;
        if (inWindow && cfg.condTwoPhaseDiag == 3) tpdBProcess(cfg, msh, var, g_tpd_nupd);   // #1b-r2: 組立 A/B の比較 (host; 読むだけ)
        tpd_collect_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(msh.nCells, g_tpd_cell, passive_periodic_root(cfg, msh), cfg.condTwoPhaseDiag,
            inWindow, (double)g_tpd_nupd, g_tpd_rec, g_tpd_n, (unsigned long long)g_tpd_cap, g_tpd_hit);
        gpuErrchk( cudaPeekAtLastError() );
        ++g_tpd_nupd;
    }
    passiveBounds_d_wrapper(cfg, cuda_cfg, msh, var, q0, (int)var.condMomentConsNames.size(), true);
    passiveMirrorPeriodic_d_wrapper(cfg, cuda_cfg, msh, var);
}

void twoPhaseUpdateLog(solverConfig& cfg, int iStep)
{
    if (!condTwoPhaseDiffusionActive(cfg) || g_tp_stats_dev == nullptr) return;
    double st[8]; int tmin = 0;
    gpuErrchk( cudaMemcpy(st, g_tp_stats_dev, 8*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(&tmin, g_tp_thetaMin_dev, sizeof(int), cudaMemcpyDeviceToHost) );
    printf("[twophase] step %d interval | updates %.0f | theta<1 cells %.0f (theta=0 %.0f) min theta %.6f | withheld vapour %.3e liquid %.3e"
           " | state corrections Qcut %.3e vround %.3e\n",
           iStep + 1, st[6], st[2], st[3], (st[2] > 0.0) ? tmin*1.0e-9 : 1.0, st[0], st[1], st[4], st[5]);
    if (g_tpd_cell != nullptr) {   // 診断: θ 制限の頻度を更新数で正規化し、区間の固有セル数を出す (#1b-r1)
        const size_t nu = g_tpd_nupd - g_tpd_nupd_log;
        int cnt[2] = {0, 0};
        gpuErrchk( cudaMemset(g_tpd_cnt, 0, 2*sizeof(int)) );
        tpd_count_d<<<(unsigned)((g_tpd_ncells + 255)/256), 256>>>((geom_int)g_tpd_ncells, g_tpd_hit, g_tpd_cnt);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        gpuErrchk( cudaMemcpy(cnt, g_tpd_cnt, 2*sizeof(int), cudaMemcpyDeviceToHost) );
        const double d = (nu > 0) ? (double)nu : 1.0;
        printf("[twophase-diag] step %d interval %zu updates | theta<1: %.0f cell-updates = %.3f cells per update, %d unique cells"
               " | theta=0: %.0f cell-updates = %.3f cells per update, %d unique cells\n",
               iStep + 1, nu, st[2], st[2]/d, cnt[0], st[3], st[3]/d, cnt[1]);
        g_tpd_nupd_log = g_tpd_nupd;
    }
    gpuErrchk( cudaMemset(g_tp_stats_dev, 0, 8*sizeof(double)) );
    const int one = 1000000000;
    gpuErrchk( cudaMemcpy(g_tp_thetaMin_dev, &one, sizeof(int), cudaMemcpyHostToDevice) );
}

// θ の全更新を覆う集計 (#1b-pre (3))。kind 0 = ソースの θ_src (condensationSource の直後; condLim_ はソースが全セルで書く),
// kind 1 = 更新の θ (θ_u の受動種更新クランプか二相の非分割更新の直後; condLim_ を更新カーネルが全セルで書く)。計上だけ。
void condThetaScan_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int kind)
{
    if (!condensationEnabled(var)) return;
    for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
        double* rs = condCorrReasons(s);
        if (rs == nullptr) continue;
        const int sN = kind ? COND_REASON_TU_N : COND_REASON_TS_N, sMin = kind ? COND_REASON_TU_MIN : COND_REASON_TS_MIN;
        const int sLast = kind ? COND_REASON_TU_LAST : COND_REASON_TS_LAST, sCalls = kind ? COND_REASON_TU_CALLS : COND_REASON_TS_CALLS;
        gpuErrchk( cudaMemset(rs + sLast, 0, sizeof(double)) );
        cond_theta_scan_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(msh.nCells, var.c_d["condLim_" + std::to_string(s)],
            periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr, rs, sN, sMin, sLast, sCalls);
    }
    gpuErrchk( cudaPeekAtLastError() );
}

// 終了時: 末尾 200 更新の θ = 0 (mode 2 は θ < 1) のセル記録を集計して CSV に書く (#1b-r1; 読むだけ)。
//   twophase_diag_theta0_cells.csv (mode 2 は theta_lt1): 頻度順の上位 500 セル。twophase_diag_theta0_summary.csv: 窓全体の要約 (key,value)。
void twoPhaseDiagWrite(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cuda_cfg;
    if (g_tpd_cell == nullptr) return;
    unsigned long long nrec = 0;
    gpuErrchk( cudaMemcpy(&nrec, g_tpd_n, sizeof(nrec), cudaMemcpyDeviceToHost) );
    const size_t nk = (size_t)std::min<unsigned long long>(nrec, g_tpd_cap);
    std::vector<double> rec(nk*TPD_REC);
    if (nk) gpuErrchk( cudaMemcpy(rec.data(), g_tpd_rec, nk*TPD_REC*sizeof(double), cudaMemcpyDeviceToHost) );
    const size_t nC = g_tpd_ncells;
    std::vector<flow_float> cx(nC), cy(nC), cz(nC);
    gpuErrchk( cudaMemcpy(cx.data(), var.c_d["ccx"], nC*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(cy.data(), var.c_d["ccy"], nC*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(cz.data(), var.c_d["ccz"], nC*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    const size_t N = g_tpd_nupd, w0 = (N > TPD_WINDOW) ? N - TPD_WINDOW : 0, W = N - w0;
    const char* tag = (cfg.condTwoPhaseDiag == 2) ? "theta_lt1" : "theta0";
    if (cfg.condTwoPhaseDiag == 3) tpdBWrite(cfg);   // #1b-r2 の CSV と判定の要約
    // 更新ごとのセル集合 (持続性) とセルごとの集計
    std::vector<std::vector<long>> setByUpd(W);
    // 窓内の最大 |値| と、その更新番号 (列ごと) を持つ。値の列は「最後の記録」ではない (codex 2026-10-02 relax-result: 最後の記録を代表値と読んだ誤り)。
    static const int KCOL[] = {2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 16};   // d_rv d_rg d_rQ2 d_rQ1 d_rQ0 res_Q2 res_Q1 res_Q0 f-1 renorm_dYw renorm_dg Rw_A Rg_A
    static const int NK = (int)(sizeof(KCOL)/sizeof(KCOL[0]));
    struct Agg { long n = 0, first = -1, last = -1, cur = 0, best = 0, persist = 0; long reason[5] = {0,0,0,0,0}; const double* lastRec = nullptr;
                 double mx[NK] = {}; long mxUpd[NK] = {}; bool any[NK] = {}; };
    std::vector<Agg> agg(nC);
    // 記録は更新番号の順に並ぶとは限らない (atomicAdd) ので更新番号で並べ直す
    std::vector<size_t> ord(nk); for (size_t k = 0; k < nk; ++k) ord[k] = k;
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { const double* ra = &rec[a*TPD_REC]; const double* rb = &rec[b*TPD_REC];
        return (ra[0] != rb[0]) ? ra[0] < rb[0] : ra[1] < rb[1]; });
    long reasonTot[5] = {0,0,0,0,0};
    for (size_t k : ord) {
        const double* r = &rec[k*TPD_REC];
        const long u = (long)r[0], c = (long)r[1];
        if (u < (long)w0 || c < 0 || (size_t)c >= nC) continue;
        setByUpd[(size_t)u - w0].push_back(c);
        Agg& a = agg[(size_t)c];
        if (a.last == u - 1) { ++a.cur; ++a.persist; } else a.cur = 1;
        a.best = std::max(a.best, a.cur);
        if (a.first < 0) a.first = u;
        a.last = u; ++a.n; a.lastRec = r;
        for (int j = 0; j < NK; ++j) {   // 符号付きの値を、|値| が窓内最大になった更新で持つ (同値は早い更新)
            const double v = r[2 + KCOL[j]];
            if (!a.any[j] || std::fabs(v) > std::fabs(a.mx[j]) || (!(v == v) && a.mx[j] == a.mx[j])) { a.mx[j] = v; a.mxUpd[j] = u; }   // 記録は更新順なので同値は早い更新
            a.any[j] = true;
        }
        const int rs = (int)r[2 + 7]; if (rs >= 0 && rs < 5) { ++a.reason[rs]; ++reasonTot[rs]; }
    }
    // 連続する更新でセル集合がどれだけ重なるか: mean_n |S_n ∩ S_{n−1}| / |S_n| (S_n 非空の n だけ)
    double pers = 0.0; long npers = 0; size_t maxPer = 0; long tot = 0;
    for (size_t i = 0; i < W; ++i) {
        std::sort(setByUpd[i].begin(), setByUpd[i].end());
        maxPer = std::max(maxPer, setByUpd[i].size()); tot += (long)setByUpd[i].size();
        if (i > 0 && !setByUpd[i].empty()) {
            std::vector<long> inter;
            std::set_intersection(setByUpd[i].begin(), setByUpd[i].end(), setByUpd[i-1].begin(), setByUpd[i-1].end(), std::back_inserter(inter));
            pers += (double)inter.size()/(double)setByUpd[i].size(); ++npers;
        }
    }
    std::vector<size_t> cells;
    for (size_t c = 0; c < nC; ++c) if (agg[c].n > 0) cells.push_back(c);
    std::sort(cells.begin(), cells.end(), [&](size_t a, size_t b) { return agg[a].n != agg[b].n ? agg[a].n > agg[b].n : a < b; });
    long nge90 = 0; for (size_t c : cells) if (W > 0 && (double)agg[c].n >= 0.9*(double)W) ++nge90;
    const std::string fc = std::string("twophase_diag_") + tag + "_cells.csv", fs = std::string("twophase_diag_") + tag + "_summary.csv";
    if (FILE* f = std::fopen(fc.c_str(), "w")) {
        std::fprintf(f, "rank,cell,x,y,z,n_updates,freq,first_update,last_update,max_consecutive,n_consecutive_with_previous,"
                        "reason_vapour_nonneg,reason_liquid_nonneg,reason_dg_max,reason_dT_max,reason_none");
        static const char* KN[] = {"d_rv","d_rg","d_rQ2","d_rQ1","d_rQ0","res_Q2","res_Q1","res_Q0","renorm_f_minus_1","renorm_d_rYw","renorm_d_rg","Rw_A","Rg_A"};
        for (int j = 0; j < NK; ++j) std::fprintf(f, ",winmaxabs_%s,winmaxabs_%s_update", KN[j], KN[j]);
        std::fprintf(f, ",last_update_rv_pre,last_update_rg_pre,last_update_theta\n");
        const size_t top = std::min<size_t>(cells.size(), 500);
        for (size_t i = 0; i < top; ++i) {
            const size_t c = cells[i]; const Agg& a = agg[c]; const double* r = a.lastRec + 2;
            std::fprintf(f, "%zu,%zu,%.9g,%.9g,%.9g,%ld,%.6g,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld",
                         i + 1, c, cx[c], cy[c], cz[c], a.n, (W ? (double)a.n/(double)W : 0.0), a.first, a.last, a.best, a.persist,
                         a.reason[1], a.reason[2], a.reason[3], a.reason[4], a.reason[0]);
            for (int j = 0; j < NK; ++j) std::fprintf(f, ",%.9e,%ld", a.mx[j], a.mxUpd[j]);
            std::fprintf(f, ",%.9e,%.9e,%.9e\n", r[0], r[1], r[11]);
        }
        std::fclose(f);
    }
    if (FILE* f = std::fopen(fs.c_str(), "w")) {
        std::fprintf(f, "key,value\n");
        std::fprintf(f, "selection,%s\nwindow_first_update,%zu\nwindow_updates,%zu\ntotal_updates,%zu\nrecords,%llu\nrecords_dropped,%llu\n",
                     tag, w0, W, N, nrec, (nrec > g_tpd_cap) ? nrec - g_tpd_cap : 0ULL);
        std::fprintf(f, "cell_updates,%ld\nunique_cells,%zu\nmean_cells_per_update,%.6g\nmax_cells_per_update,%zu\n",
                     tot, cells.size(), W ? (double)tot/(double)W : 0.0, maxPer);
        std::fprintf(f, "mean_fraction_also_in_previous_update,%.6g\ncells_in_ge_90pct_of_updates,%ld\n", npers ? pers/(double)npers : 0.0, nge90);
        std::fprintf(f, "reason_vapour_nonneg,%ld\nreason_liquid_nonneg,%ld\nreason_dg_max,%ld\nreason_dT_max,%ld\nreason_none,%ld\n",
                     reasonTot[1], reasonTot[2], reasonTot[3], reasonTot[4], reasonTot[0]);
        std::fprintf(f, "rows_in_cells_csv,%zu\n", std::min<size_t>(cells.size(), 500));
        std::fclose(f);
    }
    printf("[twophase-diag] wrote %s and %s (window %zu of %zu updates, %ld cell-updates, %zu unique cells, mean %.3f per update, persistence %.3f, records dropped %llu)\n",
           fc.c_str(), fs.c_str(), W, N, tot, cells.size(), W ? (double)tot/(double)W : 0.0, npers ? pers/(double)npers : 0.0,
           (nrec > g_tpd_cap) ? nrec - g_tpd_cap : 0ULL);
    fflush(stdout);
}
