#include "condensationTransport_d.cuh"
#include "gasPhaseComposition_d.cuh"   // gasPhaseLiquid (輸送物性の気相組成)
#include "condensationUpdateLimiter_d.cuh"   // cond_moment_update_limited_d (単体試験と共用)
#include "condensationRealizability_d.cuh"
#include "periodicNode_d.cuh"   // periodicNodeActive (実現可能性収支の root)
#include <algorithm>   // cond_realizability_clamp_d (double; 単体試験と共用)
#include "condensationSource_d.cuh"   // COND_PI, 物性 (消滅クランプ)
#include "condensationEOS_d.cuh"
#include "condensationSourceF_d.cuh"   // float 実体 (clamp の表評価)      // cond_clamp_vapor_pressure (蒸発塵判定の蒸気分圧)

#include "scalarTransport_d.cuh"
#include "passiveTransport_d.cuh"   // passiveScalarScheme 1: 化学種経路の受動種として移流・更新

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
        bool same = (g.MW == cached.gas.MW && g.Tlo == cached.gas.Tlo && g.Tmid == cached.gas.Tmid && g.Thi == cached.gas.Thi);
        for (int k = 0; k < 9; ++k) same = same && (g.low[k] == cached.gas.low[k]) && (g.high[k] == cached.gas.high[k]);
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
        const double init[2] = {1.0e300, 0.0};
        gpuErrchk( cudaMemcpy(g_condReasons_dev + (size_t)s*COND_REASON_N + COND_REASON_VMIN,   &init[0], sizeof(double), cudaMemcpyHostToDevice) );
        gpuErrchk( cudaMemcpy(g_condReasons_dev + (size_t)s*COND_REASON_N + COND_REASON_RN_MAX, &init[1], sizeof(double), cudaMemcpyHostToDevice) );
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
                                     const geom_float* vol, const geom_int* root, double* out)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (root != nullptr && root[ic] != ic) return;
    const double V = (double)vol[ic];
    if (rog[ic] != (flow_float)0.0) atomicAdd(&out[0], (double)rog[ic]*V);
    if (q1[ic]  != (flow_float)0.0) atomicAdd(&out[1], (double)q1[ic]*V);
    if (q2[ic]  != (flow_float)0.0) atomicAdd(&out[2], (double)q2[ic]*V);
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
    if (tot_d == nullptr) gpuErrchk( cudaMalloc((void**)&tot_d, 3*sizeof(double)) );
    const geom_int* root = periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr;
    const int nstep = std::max(1, (iStep + 1) - g_condReasons_last_step);
    for (int s = 0; s < nsp && s < g_condReasons_n; ++s) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemset(tot_d, 0, 3*sizeof(double)) );
        cond_liquid_totals_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells, var.c_d["rog_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i], var.c_d["volume"], root, tot_d);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        double tot[3] = {0.0, 0.0, 0.0};
        gpuErrchk( cudaMemcpy(tot, tot_d, 3*sizeof(double), cudaMemcpyDeviceToHost) );
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
    g_condReasons_last = cur;
    if (passive) s_pst_last = pst;
    g_condReasons_last_step = iStep + 1;
    condReasonsResetInterval(g_condReasons_n);
    fflush(stdout);
}
