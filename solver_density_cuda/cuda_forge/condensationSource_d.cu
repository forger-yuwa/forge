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

// 相変化ソース kernel。pure (N2/CPG) と carrier (H2O/TP) の両対応。一温度 T_v=T_d=T。
// J,r*,dr/dt は現在セル状態から freeze。安定化: J 上限、dr/dt<0→0、r̄≤r* 成長停止、g≤g_max、
//   1 step の Δg・潜熱 ΔT・蒸気枯渇を θ で律速 (全モーメント同 θ)。src_jac=潜熱自己抑制 (g)。
__global__ void condensation_source_d(
    geom_int nCells,
    int condModel, int carrier, double Rw, double M,
    int kantrowitz, int kwGammaMode, CondPropOpts opts,   // opts: σ 倍率・N2 低温物性切替・kgas モデル・CPG carrier の Y_w (plans/accepted/condensation-air.md)
    const SpeciesThermo* sp, int nSpecies, flow_float* const* roYall, int condGasSpecies,   // Feder carrier 形の衝突項用 (TP のみ; CPG は nullptr)
    int growthModel, double gyarC, int twoTemp,
    int evap, double evapRmin, int evapKelvin, double evapLamMin,
    int eq, double eqRelax, double eqDgMax, double eqDTmax,
    flow_float cp_cpg, flow_float gamma_cpg,
    double Jmax, double dg_max, double dT_max,
    geom_float* vol, flow_float* dt_local,
    flow_float* T, flow_float* P, flow_float* ro, flow_float* cp_cell, flow_float* Rmix_cell,
    flow_float* roY_w,   // carrier: 凝縮気相種の保存量 ρY_w (pure では nullptr)
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    flow_float* res_rog, flow_float* res_roQ0, flow_float* res_roQ1, flow_float* res_roQ2,
    flow_float* sj_g, flow_float* sj_Q0, flow_float* sj_Q1, flow_float* sj_Q2,
    flow_float* diagS, flow_float* diagDrdt, flow_float* diagR30, flow_float* diagTsat,
    flow_float* diagTheta, flow_float* diagLim)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    sj_Q0[ic] = 0.0; sj_Q1[ic] = 0.0; sj_Q2[ic] = 0.0; sj_g[ic] = 0.0;
    diagS[ic] = 0.0; diagDrdt[ic] = 0.0; diagR30[ic] = 0.0; diagTsat[ic] = 0.0;
    diagTheta[ic] = 0.0; diagLim[ic] = 1.0;   // θ (非等温補正) / ソース律速係数 (1=律速なし)

    const double rod = (double)ro[ic];
    if (rod <= 1.0e-20) return;

    const CondSpeciesProps cprops = condProps_make(condModel, opts);   // σ 倍率・N2 低温物性・kgas モデルを反映
    const double Td = (double)T[ic];
    const double Pd = (double)P[ic];

    // 気相熱容量 (ΔT/src_jac 用)。carrier(TP)=cp_cell、pure(CPG)=cp_cpg。
    const double cpg = (carrier && cp_cell) ? (double)cp_cell[ic] : (double)cp_cpg;
    const double Rg  = (carrier && Rmix_cell) ? (double)Rmix_cell[ic] : ((double)gamma_cpg-1.0)*(double)cp_cpg/(double)gamma_cpg;
    const double cvg = (cpg - Rg > 1.0e-3) ? (cpg - Rg) : 1.0e-3;

    double g = (double)rog[ic]/rod; if (g < 0.0) g = 0.0;
    // carrier の凝縮種質量分率: TP は roY_w (輸送), CPG carrier (空気の N2 選択凝縮) は config 定数 opts.Yw
    const double Yw = (carrier && roY_w) ? (double)roY_w[ic]/rod : ((carrier && opts.Yw > 0.0) ? opts.Yw : 1.0);
    const double gmax = carrier ? Yw : 0.99;
    if (g > gmax) g = (gmax > 0.0 ? gmax : 0.0);
    double q0 = (double)roQ0[ic]; if (q0 < 0.0) q0 = 0.0;
    double q1 = (double)roQ1[ic]; if (q1 < 0.0) q1 = 0.0;
    double q2 = (double)roQ2[ic]; if (q2 < 0.0) q2 = 0.0;

    // 蒸気状態
    double pv, rho_v;
    cond_vapor_state(carrier, rod, Pd, Td, g, Yw, Rw, &pv, &rho_v);

    // Kantrowitz 用の比熱比: 既定は凝縮種 (蒸気) 自身の γ_v (純蒸気形の係数 2(γ_v−1)/(γ_v+1)=R_v/(c_v,v+R_v/2))。
    // kwGammaMode=1 は旧挙動 (セル気相混合 cp/cv。carrier では N2 支配 ≈1.40 で θ が過大) を A/B 用に残す。
    // Gyarmathy Kn 用の全圧 (carrier=Pd, pure=pv)
    const double gamma_gas = (kwGammaMode == 1) ? (cpg/cvg) : (cprops.cp/cprops.cv);
    const double p_gas = carrier ? Pd : pv;

    // 診断: 過飽和 S=p_v/p_sat(T)
    const double psat_T = cond_psat(cprops, Td);
    diagS[ic] = (flow_float)(pv/(psat_T > 1.0e-300 ? psat_T : 1.0e-300));
    diagTsat[ic] = (flow_float)cond_Tsat(cprops, pv, Td);

    // Feder carrier 形 (condKantrowitz 2/3) の衝突項: 種 DB (NASA-9 c_v(T), M_i) と種質量分率から種別に集計
    // (plans/accepted/condensation-kantrowitz-carrier.md §4.1)。TP carrier 以外 (CPG / pure) は純蒸気形 (carrierSum=0)。
    // 温度摂動 (src_jac) でも同じ値を使う (T 凍結)。
    CondNucCarrier car;
    car.a_v = 1.0; car.carrierSum = 0.0; car.cvv_tilde = cprops.cv*cprops.M/COND_RU;
    if (carrier && sp != nullptr && roYall != nullptr && condGasSpecies >= 0 && condGasSpecies < nSpecies) {
        const SpeciesThermo& sv = sp[condGasSpecies];
        const double Mv = sv.MW;
        car.cvv_tilde = (thermo_cp_mass(sv, Td) - thermo_R_species(sv))*Mv/COND_RU;
        car.a_v = ((Yw - g) > 0.0 ? (Yw - g) : 0.0)/Mv;
        double sum = 0.0;
        for (int i = 0; i < nSpecies; ++i) {
            if (i == condGasSpecies) continue;
            const double Yi = (double)roYall[i][ic]/rod;
            if (!(Yi > 0.0)) continue;
            const double Mi  = sp[i].MW;
            const double cvi = (thermo_cp_mass(sp[i], Td) - thermo_R_species(sp[i]))*Mi/COND_RU;
            sum += (Yi/Mi)*sqrt(Mv/Mi)*(cvi + 0.5);
        }
        car.carrierSum = sum;
    }
    if (kantrowitz && pv > psat_T) {
        const double gamma_kw = (kwGammaMode == 1) ? (cpg/cvg) : (cprops.cp/cprops.cv);
        diagTheta[ic] = (flow_float)cond_kantrowitz_theta(cprops, Td, log(pv/psat_T), kantrowitz, gamma_kw, &car);
    }

    // ---- 平衡凝縮・EOS 拘束形 (condEquilibrium=2): g は dependentVariables が (T,g) 同時反転で決めて rog に
    //      射影済み (plans/accepted/condensation-equilibrium-eos.md)。ここでは診断 (上の S, Tsat) だけ書き、rog の輸送残差を
    //      0 に凍結して point-implicit 更新 δ=0 にする (rms_rog 恒等 0)。Q0-Q2 は緩和形と同じく輸送のみ。
    if (eq == 2) {
        res_rog[ic] = (flow_float)0.0;
        sj_g[ic]    = (flow_float)0.0;
        return;
    }

    // ---- 平衡凝縮分岐 (condEquilibrium=1): g_eq への緩和。モーメント Q0-Q2 はソース 0 (輸送のみ) ----
    if (eq == 1) {
        const double yv0 = carrier ? (Yw - g) : (1.0 - g);
        const double D = cond_equilibrium_delta(cprops, rod, Td, g, (yv0 > 0.0 ? yv0 : 0.0), Rw, cvg);
        const double dt = (double)dt_local[ic];
        if (dt <= 0.0 || D == 0.0) return;
        double Sg = eqRelax*rod*D/dt;
        // θ 律速 (非平衡と共通の規約): |Δg|<=dg_max, 潜熱 |ΔT|<=dT_max, 蒸気/液の枯渇
        const double L  = cond_latent(cprops, Td);
        double theta = 1.0;
        {
            const double dg  = fabs(Sg)*dt/rod;
            const double dTl = dg*L/cvg;
            if (dg  > eqDgMax) theta = fmin(theta, eqDgMax/dg);
            if (dTl > eqDTmax) theta = fmin(theta, eqDTmax/dTl);
        }
        Sg *= theta;
        sj_g[ic] = (flow_float)(eqRelax*theta/dt);       // ∂S_g/∂(ρg) = -αθ/dt (負帰還)
        res_rog[ic] += (flow_float)(Sg*(double)vol[ic]);
        return;
    }

    // ---- 蒸発分岐 (condEvaporation=1, S<=1, 液相あり): 負成長 λ スケール + 液滴消滅 ----
    // plans/accepted/condensation-evaporation.md §4-5。核生成は S<=1 で 0 なので本分岐では評価しない。
    if (evap && g > 0.0 && pv <= psat_T) {
        const double dt = (double)dt_local[ic];
        double SQ0, SQ1, SQ2, Sg, r30, drdt;
        cond_evap_source(cprops, Td, pv, rod, g, q0, q1, q2, dt,
                         evapRmin, evapLamMin, dg_max, dT_max, cvg,
                         growthModel, p_gas, gyarC, evapKelvin,
                         &SQ0, &SQ1, &SQ2, &Sg, &r30, &drdt);
        diagDrdt[ic] = (flow_float)drdt; diagR30[ic] = (flow_float)r30;
        // src_jac (g 行のみ): 蒸気復帰 ∂S_g/∂(ρg) (carrier: p_v=ρ(Y_w-g)R_wT) と潜熱冷却
        // ∂S_g/∂T·∂T/∂(ρg) を数値微分で。どちらも負帰還 (S_g<0 が g↓/T↓で 0 に近づく) → sj_g>=0。
        // Q1 行の ∂S/∂Q1 は Kelvin 正帰還側なので陰的に入れない (sj_Q1=0)。
        if (Sg < 0.0) {
            const double L  = cond_latent(cprops, Td);
            double a0,a1,a2,ag,rr,dd;
            // (a) g 摂動 (蒸気状態も更新)
            const double dg = 1.0e-3*g;
            double pvg, rvg; cond_vapor_state(carrier, rod, Pd, Td, g - dg, Yw, Rw, &pvg, &rvg);
            cond_evap_source(cprops, Td, pvg, rod, g - dg, q0, q1, q2, dt,
                             evapRmin, evapLamMin, dg_max, dT_max, cvg,
                             growthModel, p_gas, gyarC, evapKelvin, &a0,&a1,&a2,&ag,&rr,&dd);
            const double dSgdrog = (Sg - ag)/(rod*dg);        // ∂S_g/∂(ρg)
            // (b) T 摂動
            const double dTp = 0.1;
            double pvT, rvT; cond_vapor_state(carrier, rod, Pd, Td+dTp, g, Yw, Rw, &pvT, &rvT);
            cond_evap_source(cprops, Td+dTp, pvT, rod, g, q0, q1, q2, dt,
                             evapRmin, evapLamMin, dg_max, dT_max, cvg,
                             growthModel, p_gas, gyarC, evapKelvin, &a0,&a1,&a2,&ag,&rr,&dd);
            const double dSgdT  = (ag - Sg)/dTp;
            const double dTdrog = (L - (carrier?Rw:Rg)*Td)/(rod*cvg);
            double sjg = -(dSgdrog + dSgdT*dTdrog);
            if (sjg < 0.0) sjg = 0.0;
            sj_g[ic] = (flow_float)sjg;
        }
        const double v = (double)vol[ic];
        res_roQ0[ic] += (flow_float)(SQ0*v);
        res_roQ1[ic] += (flow_float)(SQ1*v);
        res_roQ2[ic] += (flow_float)(SQ2*v);
        res_rog[ic]  += (flow_float)(Sg *v);
        return;
    }

    // 核生成・成長 (freeze, 亜臨界停止・蒸発 clamp)
    double J, rstar;
    cond_nucleation(cprops, Td, pv, rho_v, &J, &rstar, kantrowitz, gamma_gas, &car);
    if (J > Jmax) J = Jmax;
    if (J < 0.0)  J = 0.0;
    double r_bar = (q0 > 1.0e-30) ? (q1/q0) : rstar;
    double drdt = 0.0;
    if (q0 > 1.0e-30 && rstar > 0.0 && r_bar > rstar) {
        drdt = cond_growth(cprops, Td, pv, r_bar, rstar, growthModel, p_gas, gyarC, twoTemp);
        if (drdt < 0.0) drdt = 0.0;
    }
    diagDrdt[ic] = (flow_float)drdt;
    const double rho_l = cond_rho_cond(cprops, Td);
    const double r_nuc = COND_RNUC_FAC*rstar;   // わずかに超臨界で生み成長を起動 (r*ちょうどだと(1-r*/r)=0で停止)
    double SQ0 = J;
    double SQ1 = J*r_nuc + q0*drdt;
    double SQ2 = J*r_nuc*r_nuc + 2.0*q1*drdt;
    double Sg  = (4.0/3.0)*COND_PI*rho_l*(J*r_nuc*r_nuc*r_nuc + 3.0*q2*drdt);
    if (Sg < 0.0) Sg = 0.0;

    // θ 律速: Δg, 潜熱 ΔT, 蒸気枯渇 (carrier は利用可能蒸気=Yw-g)
    const double dt = (double)dt_local[ic];
    const double L  = cond_latent(cprops, Td);
    double theta = 1.0;
    if (Sg > 0.0 && dt > 0.0) {
        const double dg  = Sg*dt/rod;
        const double dTl = dg*L/cvg;                       // 潜熱による ΔT (気相 cv で割る)
        const double avail = carrier ? (Yw - g) : (1.0 - g);
        if (dg  > dg_max)    theta = fmin(theta, dg_max/dg);
        if (dTl > dT_max)    theta = fmin(theta, dT_max/dTl);
        if (avail > 0.0 && dg > 0.9*avail) theta = fmin(theta, 0.9*avail/fmax(dg,1.0e-300));
        else if (avail <= 0.0) theta = 0.0;
    }

    // src_jac: 潜熱自己抑制 (g↑→T↑→psat↑→S↓→S_g↓)。∂T/∂(ρg)=(L-R_w T)/(ρ c_vg)。∂S_g/∂T 数値。
    {
        const double dTp = 0.1;
        double pvp, rvp;
        cond_vapor_state(carrier, rod, Pd, Td+dTp, g, Yw, Rw, &pvp, &rvp);
        double a0,a1,a2,ag;
        cond_source_vector(cprops, Td+dTp, pvp, rvp, q0, q1, q2, &a0,&a1,&a2,&ag,
                           kantrowitz, growthModel, gamma_gas, p_gas, gyarC, twoTemp, &car);   // 本体と同じ carrier 項 (T 凍結)
        if (ag < 0.0) ag = 0.0;
        const double dSgdT  = (ag - Sg)/dTp;
        const double dTdrog = (L - (carrier?Rw:Rg)*Td)/(rod*cvg);
        double sjg = -theta*dSgdT*dTdrog;
        if (sjg < 0.0) sjg = 0.0;
        sj_g[ic] = (flow_float)sjg;
        if (q0 > 1.0e-30) {
            const double dq1 = (q1 > 0.0 ? 0.01*q1 : 1.0e-3);
            double b0,b1,b2,bg;
            cond_source_vector(cprops, Td, pv, rho_v, q0, q1+dq1, q2, &b0,&b1,&b2,&bg,
                               kantrowitz, growthModel, gamma_gas, p_gas, gyarC, twoTemp, &car);
            double sjq1 = -theta*(b1 - SQ1)/dq1;
            if (sjq1 < 0.0) sjq1 = 0.0;
            sj_Q1[ic] = (flow_float)sjq1;
        }
    }

    diagLim[ic] = (flow_float)theta;   // ソース律速係数 (核生成域で 1 に近いことを確認する診断)
    SQ0 *= theta; SQ1 *= theta; SQ2 *= theta; Sg *= theta;
    const double v = (double)vol[ic];
    res_roQ0[ic] += (flow_float)(SQ0*v);
    res_roQ1[ic] += (flow_float)(SQ1*v);
    res_roQ2[ic] += (flow_float)(SQ2*v);
    res_rog[ic]  += (flow_float)(Sg *v);
}


// ---------------------------------------------------------------------------------------------------------------
// float 実体 (condFloat=1, condEquilibrium=0, condTwoTemp=0 のとき wrapper が選ぶ)。上の double kernel の写し。
//   差: 物性は表、核生成は対数空間 (上限を対数で本体・摂動の両方に)、dry セル (S<=1, g=0, Q0=0) は診断だけ書いて早期退出、
//       T_sat は前 step 値から warm start。sj_*/diag* の初期化と res_* を触らない規約は double と同じ。
// ---------------------------------------------------------------------------------------------------------------
__global__ void condensation_source_f_d(
    geom_int nCells,
    int carrier, float Rw,
    int kantrowitz, int kwGammaMode, CondSpeciesPropsF cpf, CondTablesF tb, float Yw_const,
    const SpeciesThermoF* spf, int nSpecies, flow_float* const* roYall, int condGasSpecies,
    int growthModel, float gyarC,
    int evap, float evapRmin, int evapKelvin, float evapLamMin,
    float cp_cpg, float gamma_cpg,
    float dg_max, float dT_max,
    geom_float* vol, flow_float* dt_local,
    flow_float* T, flow_float* P, flow_float* ro, flow_float* cp_cell, flow_float* Rmix_cell,
    flow_float* roY_w,
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    flow_float* res_rog, flow_float* res_roQ0, flow_float* res_roQ1, flow_float* res_roQ2,
    flow_float* sj_g, flow_float* sj_Q0, flow_float* sj_Q1, flow_float* sj_Q2,
    flow_float* diagS, flow_float* diagDrdt, flow_float* diagR30, flow_float* diagTsat,
    flow_float* diagTheta, flow_float* diagLim)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const float Tsat_prev = diagTsat[ic];   // warm start (初期化前に読む)
    sj_Q0[ic] = 0.0f; sj_Q1[ic] = 0.0f; sj_Q2[ic] = 0.0f; sj_g[ic] = 0.0f;
    diagS[ic] = 0.0f; diagDrdt[ic] = 0.0f; diagR30[ic] = 0.0f; diagTsat[ic] = 0.0f;
    diagTheta[ic] = 0.0f; diagLim[ic] = 1.0f;
    const float rod = ro[ic];
    if (rod <= 1.0e-20f) return;
    const float Td = T[ic];
    const float Pd = P[ic];
    const float cpg = (carrier && cp_cell) ? cp_cell[ic] : cp_cpg;
    const float Rg  = (carrier && Rmix_cell) ? Rmix_cell[ic] : (gamma_cpg-1.0f)*cp_cpg/gamma_cpg;
    const float cvg = (cpg - Rg > 1.0e-3f) ? (cpg - Rg) : 1.0e-3f;
    float g = rog[ic]/rod; if (g < 0.0f) g = 0.0f;
    const float Yw = (carrier && roY_w) ? roY_w[ic]/rod : ((carrier && Yw_const > 0.0f) ? Yw_const : 1.0f);
    const float gmax = carrier ? Yw : 0.99f;
    if (g > gmax) g = (gmax > 0.0f ? gmax : 0.0f);
    float q0 = roQ0[ic]; if (q0 < 0.0f) q0 = 0.0f;
    float q1 = roQ1[ic]; if (q1 < 0.0f) q1 = 0.0f;
    float q2 = roQ2[ic]; if (q2 < 0.0f) q2 = 0.0f;
    float pv, rho_v;
    cond_vapor_state_f(carrier, rod, Pd, Td, g, Yw, Rw, &pv, &rho_v);
    const float gamma_gas = (kwGammaMode == 1) ? (cpg/cvg) : (cpf.cp/cpf.cv);
    const float p_gas = carrier ? Pd : pv;
    // 診断: S = p_v/p_sat, T_sat
    const float lnpsat = cond_tab_lnpsat_f(tb, Td);
    const float lnS = (pv > 0.0f) ? (logf(pv) - lnpsat) : -1.0e30f;
    diagS[ic]    = (pv > 0.0f) ? expf(lnS) : 0.0f;
    diagTsat[ic] = cond_Tsat_f(tb, pv, (Tsat_prev > 0.0f) ? Tsat_prev : Td);
    // Feder carrier 形の衝突項 (TP carrier のみ)
    CondNucCarrierF car;
    car.a_v = 1.0f; car.carrierSum = 0.0f; car.cvv_tilde = cpf.cv*cpf.M/COND_RU_F;
    if (carrier && spf != nullptr && roYall != nullptr && condGasSpecies >= 0 && condGasSpecies < nSpecies) {
        const SpeciesThermoF& sv = spf[condGasSpecies];
        const float Mv = sv.MW;
        car.cvv_tilde = thermo_cp_molar_f(sv, Td)/COND_RU_F - 1.0f;       // (c_p,mass − R_v) M_v/R_u = c_p,molar/R_u − 1
        car.a_v = ((Yw - g) > 0.0f ? (Yw - g) : 0.0f)/Mv;
        float sum = 0.0f;
        for (int i = 0; i < nSpecies; ++i) {
            if (i == condGasSpecies) continue;
            const float Yi = roYall[i][ic]/rod;
            if (!(Yi > 0.0f)) continue;
            const float Mi  = spf[i].MW;
            const float cvi = thermo_cp_molar_f(spf[i], Td)/COND_RU_F - 1.0f;
            sum += (Yi/Mi)*sqrtf(Mv/Mi)*(cvi + 0.5f);
        }
        car.carrierSum = sum;
    }
    if (kantrowitz && lnS > 0.0f) {
        const float gamma_kw = (kwGammaMode == 1) ? (cpg/cvg) : (cpf.cp/cpf.cv);
        diagTheta[ic] = cond_kantrowitz_theta_f(cpf, tb, Td, lnS, kantrowitz, gamma_kw, &car);
    }
    // dry セルの早期退出 (plan §4.2-3): S<=1 (J=0), g=0 (蒸発なし), Q0=0 (成長なし) → ソース・src_jac とも恒等 0。res_* は触らない。
    if (!(lnS > 0.0f) && g <= 0.0f && q0 <= 1.0e-30f) return;
    // ---- 蒸発分岐 (S<=1, 液相あり) ----
    if (evap && g > 0.0f && !(lnS > 0.0f)) {
        const float dt = dt_local[ic];
        float SQ0, SQ1, SQ2, Sg, r30, drdt;
        cond_evap_source_f(cpf, tb, Td, pv, rod, g, q0, q1, q2, dt,
                           evapRmin, evapLamMin, dg_max, dT_max, cvg,
                           growthModel, p_gas, gyarC, evapKelvin,
                           &SQ0, &SQ1, &SQ2, &Sg, &r30, &drdt);
        diagDrdt[ic] = drdt; diagR30[ic] = r30;
        if (Sg < 0.0f) {
            const float L = cond_tab_latent_f(tb, Td);
            float a0,a1,a2,ag,rr,dd;
            const float dg = 1.0e-3f*g;
            float pvg, rvg; cond_vapor_state_f(carrier, rod, Pd, Td, g - dg, Yw, Rw, &pvg, &rvg);
            cond_evap_source_f(cpf, tb, Td, pvg, rod, g - dg, q0, q1, q2, dt,
                               evapRmin, evapLamMin, dg_max, dT_max, cvg,
                               growthModel, p_gas, gyarC, evapKelvin, &a0,&a1,&a2,&ag,&rr,&dd);
            const float dSgdrog = (Sg - ag)/(rod*dg);
            const float dTp = 0.1f;
            float pvT, rvT; cond_vapor_state_f(carrier, rod, Pd, Td+dTp, g, Yw, Rw, &pvT, &rvT);
            cond_evap_source_f(cpf, tb, Td+dTp, pvT, rod, g, q0, q1, q2, dt,
                               evapRmin, evapLamMin, dg_max, dT_max, cvg,
                               growthModel, p_gas, gyarC, evapKelvin, &a0,&a1,&a2,&ag,&rr,&dd);
            const float dSgdT  = (ag - Sg)/dTp;
            const float dTdrog = (L - (carrier?Rw:Rg)*Td)/(rod*cvg);
            float sjg = -(dSgdrog + dSgdT*dTdrog);
            if (sjg < 0.0f) sjg = 0.0f;
            sj_g[ic] = sjg;
        }
        const float v = vol[ic];
        res_roQ0[ic] += SQ0*v;
        res_roQ1[ic] += SQ1*v;
        res_roQ2[ic] += SQ2*v;
        res_rog[ic]  += Sg *v;
        return;
    }
    // ---- 核生成・成長 ----
    float J, rstar;
    cond_nucleation_f(cpf, tb, Td, pv, rho_v, &J, &rstar, kantrowitz, gamma_gas, &car);   // 上限は関数内 (対数)
    if (J < 0.0f) J = 0.0f;
    float r_bar = (q0 > 1.0e-30f) ? (q1/q0) : rstar;
    float drdt = 0.0f;
    if (q0 > 1.0e-30f && rstar > 0.0f && r_bar > rstar) {
        drdt = cond_growth_f(cpf, tb, Td, pv, r_bar, rstar, growthModel, p_gas, gyarC);
        if (drdt < 0.0f) drdt = 0.0f;
    }
    diagDrdt[ic] = drdt;
    const float rho_l = cond_tab_rhol_f(tb, Td);
    const float r_nuc = COND_RNUC_FAC_F*rstar;
    float SQ0 = J;
    float SQ1 = J*r_nuc + q0*drdt;
    float SQ2 = J*r_nuc*r_nuc + 2.0f*q1*drdt;
    float Sg  = (4.0f/3.0f)*COND_PI_F*rho_l*(J*r_nuc*r_nuc*r_nuc + 3.0f*q2*drdt);
    if (Sg < 0.0f) Sg = 0.0f;
    const float dt = dt_local[ic];
    const float L  = cond_tab_latent_f(tb, Td);
    float theta = 1.0f;
    if (Sg > 0.0f && dt > 0.0f) {
        const float dg  = Sg*dt/rod;
        const float dTl = dg*L/cvg;
        const float avail = carrier ? (Yw - g) : (1.0f - g);
        if (dg  > dg_max)    theta = fminf(theta, dg_max/dg);
        if (dTl > dT_max)    theta = fminf(theta, dT_max/dTl);
        if (avail > 0.0f && dg > 0.9f*avail) theta = fminf(theta, 0.9f*avail/fmaxf(dg, 1.0e-30f));
        else if (avail <= 0.0f) theta = 0.0f;
    }
    {
        const float dTp = 0.1f;
        float pvp, rvp;
        cond_vapor_state_f(carrier, rod, Pd, Td+dTp, g, Yw, Rw, &pvp, &rvp);
        float a0,a1,a2,ag;
        cond_source_vector_f(cpf, tb, Td+dTp, pvp, rvp, q0, q1, q2, &a0,&a1,&a2,&ag,
                             kantrowitz, growthModel, gamma_gas, p_gas, gyarC, &car);
        if (ag < 0.0f) ag = 0.0f;
        const float dSgdT  = (ag - Sg)/dTp;
        const float dTdrog = (L - (carrier?Rw:Rg)*Td)/(rod*cvg);
        float sjg = -theta*dSgdT*dTdrog;
        if (sjg < 0.0f) sjg = 0.0f;
        sj_g[ic] = sjg;
        if (q0 > 1.0e-30f) {
            const float dq1 = (q1 > 0.0f ? 0.01f*q1 : 1.0e-3f);
            float b0,b1,b2,bg;
            cond_source_vector_f(cpf, tb, Td, pv, rho_v, q0, q1+dq1, q2, &b0,&b1,&b2,&bg,
                                 kantrowitz, growthModel, gamma_gas, p_gas, gyarC, &car);
            float sjq1 = -theta*(b1 - SQ1)/dq1;
            if (sjq1 < 0.0f) sjq1 = 0.0f;
            sj_Q1[ic] = sjq1;
        }
    }
    diagLim[ic] = theta;
    SQ0 *= theta; SQ1 *= theta; SQ2 *= theta; Sg *= theta;
    const float v = vol[ic];
    res_roQ0[ic] += SQ0*v;
    res_roQ1[ic] += SQ1*v;
    res_roQ2[ic] += SQ2*v;
    res_rog[ic]  += Sg *v;
}

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
    const double dg_max = 5.0e-3;
    const double dT_max = 1.0;
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
            condensation_source_f_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                msh.nCells,
                carrier, (float)Rw,
                cfg.condKantrowitz, cfg.condKantrowitzGammaMode, condProps_to_f(cprops), cond_tables_device(), (float)opts.Yw,
                (cfg.thermalMethod == 2) ? thermo_species_device_ptr_f() : nullptr, cfg.nSpecies,
                (cfg.thermalMethod == 2) ? species_roY_device_ptr() : nullptr, cfg.condGasSpecies,
                cfg.condGrowthModel, (float)cfg.condGyarmathyC,
                cfg.condEvaporation, (float)cfg.condEvapRmin, cfg.condEvapKelvin, (float)evapLamMin,
                cfg.cp, cfg.gamma,
                (float)dg_max, (float)dT_max,
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
            Jmax, dg_max, dT_max,
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
