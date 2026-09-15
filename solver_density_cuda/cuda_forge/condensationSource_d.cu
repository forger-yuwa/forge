#include "condensationTransport_d.cuh"  // wrapper 宣言 + 必要なクラス
#include "condensationSource_d.cuh"
#include "thermo_d.cuh"            // SpeciesThermo / thermo_cp_mass (Feder carrier 形)
#include "speciesTransport_d.cuh"  // species_roY_device_ptr()
#include "condensationEOS_d.cuh"       // cond_equilibrium_delta (緩和形平衡)
#include "condensationSourceF_d.cuh"   // float 実体 (物性表・対数 CNT; plans/active/condensation-float-speedup.md)

#include <string>

namespace {

inline bool condensationEnabled(const solverConfig& cfg)
{
    return cfg.condensation == 1 && cfg.nCondSpecies >= 1;
}

// cond_vapor_state は condensationSource_d.cuh へ移動 (単体試験 test_cond_kantrowitz_carrier (b3) が kernel と同じ蒸気状態再評価を呼ぶため)。

// cond_equilibrium_delta (緩和形平衡の Δ) は condensationEOS_d.cuh へ移動 (EOS 拘束形と単体テストで共用)。

#include "condensationSourceKernels_d.cuh"   // kernel 本体 (double / float 実体)
}  // namespace

void condensationSource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condensationEnabled(cfg)) return;

    const int carrier = (cfg.condGasSpecies >= 0 || cfg.condVaporMassFraction > 0.0) ? 1 : 0;   // TP carrier / CPG carrier (空気)
    const CondPropOpts opts = cond_prop_opts(cfg);
    const CondSpeciesProps cprops = condProps_make(cfg.condModel, opts);
    const double M  = cprops.M;
    const double Rw = cprops.R;
    const double Jmax   = 1.0e35;
    const double dg_max = cfg.condDgMaxStep;   // 1 更新あたりの上限 (mode 1 では更新クランプ、mode 0 では残差 θ)
    const double dT_max = cfg.condDTmaxStep;
    const int limiterMode = cfg.condLimiterMode;
    const double evapLamMin = 0.5;   // 蒸発: 1 step の半径縮小比の下限 (半減)。§5.1-3

    for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        flow_float* roY_w = nullptr;
        if (carrier) roY_w = var.c_d["roY" + std::to_string(cfg.condGasSpecies)];
        // CPG (thermalMethod!=2) では per-cell cp/Rmix 配列は未充填。nullptr で cfg.cp/γ フォールバック
        // (= carrier N2 の cp/R)。TP のみ per-cell 配列を渡す。
        flow_float* cp_cell   = (cfg.thermalMethod == 2) ? var.c_d["cp"]   : nullptr;
        flow_float* Rmix_cell = (cfg.thermalMethod == 2) ? var.c_d["Rmix"] : nullptr;
        // float 実体 (plans/active/condensation-float-speedup.md §4.2-6 の分岐表): condFloat=1 かつ平衡形/二温度でない (それらは double のまま)。
        const bool useFloat = (cfg.condFloat != 0) && cfg.condEquilibrium == 0 && cfg.condTwoTemp == 0 && cond_tables_device().valid;
        if (useFloat) {
            CondDoubleArgs dbl;
            dbl.opts = opts; dbl.sp = (cfg.thermalMethod == 2) ? thermo_species_device_ptr() : nullptr; dbl.condModel = cfg.condModel;
            dbl.Rw = Rw; dbl.M = M; dbl.twoTemp = cfg.condTwoTemp; dbl.gyarC = cfg.condGyarmathyC; dbl.evapRmin = cfg.condEvapRmin;
            dbl.evapLamMin = evapLamMin; dbl.Jmax = Jmax; dbl.dg_max = dg_max; dbl.dT_max = dT_max; dbl.cprops = cprops;
            condensation_source_f_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                msh.nCells,
                carrier, (float)Rw,
                cfg.condKantrowitz, cfg.condKantrowitzGammaMode, condProps_to_f(cprops), cond_tables_device(), (float)opts.Yw, dbl,
                (cfg.thermalMethod == 2) ? thermo_species_device_ptr_f() : nullptr, cfg.nSpecies,
                (cfg.thermalMethod == 2) ? species_roY_device_ptr() : nullptr, cfg.condGasSpecies,
                cfg.condGrowthModel, (float)cfg.condGyarmathyC,
                cfg.condEvaporation, (float)cfg.condEvapRmin, cfg.condEvapKelvin, (float)evapLamMin,
                cfg.cp, cfg.gamma,
                (float)dg_max, (float)dT_max, limiterMode,
                var.c_d["volume"], var.c_d["dt_local"],
                var.c_d["T"], var.c_d["P"], var.c_d["ro"], cp_cell, Rmix_cell,
                roY_w,
                var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
                var.c_d["res_rog_"+i], var.c_d["res_roQ0_"+i], var.c_d["res_roQ1_"+i], var.c_d["res_roQ2_"+i],
                var.c_d["src_jac_g_"+i], var.c_d["src_jac_Q0_"+i], var.c_d["src_jac_Q1_"+i], var.c_d["src_jac_Q2_"+i],
                var.c_d["condS_"+i], var.c_d["condDrdt_"+i], var.c_d["condR30_"+i], var.c_d["condTsat_"+i],
                var.c_d["condTheta_"+i], var.c_d["condLim_"+i]);
            continue;
        }
        condensation_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells,
            cfg.condModel, carrier, Rw, M,
            cfg.condKantrowitz, cfg.condKantrowitzGammaMode, opts,
            (cfg.thermalMethod == 2) ? thermo_species_device_ptr() : nullptr, cfg.nSpecies,
            (cfg.thermalMethod == 2) ? species_roY_device_ptr() : nullptr, cfg.condGasSpecies,
            cfg.condGrowthModel, cfg.condGyarmathyC, cfg.condTwoTemp,
            cfg.condEvaporation, cfg.condEvapRmin, cfg.condEvapKelvin, evapLamMin,
            cfg.condEquilibrium, cfg.condEqRelax, cfg.condEqDgMax, cfg.condEqDTmax,
            cfg.cp, cfg.gamma,
            Jmax, dg_max, dT_max, limiterMode,
            var.c_d["volume"], var.c_d["dt_local"],
            var.c_d["T"], var.c_d["P"], var.c_d["ro"], cp_cell, Rmix_cell,
            roY_w,
            var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
            var.c_d["res_rog_"+i], var.c_d["res_roQ0_"+i], var.c_d["res_roQ1_"+i], var.c_d["res_roQ2_"+i],
            var.c_d["src_jac_g_"+i], var.c_d["src_jac_Q0_"+i], var.c_d["src_jac_Q1_"+i], var.c_d["src_jac_Q2_"+i],
            var.c_d["condS_"+i], var.c_d["condDrdt_"+i], var.c_d["condR30_"+i], var.c_d["condTsat_"+i],
            var.c_d["condTheta_"+i], var.c_d["condLim_"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
