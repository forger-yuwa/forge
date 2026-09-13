#pragma once
// 非平衡凝縮ソースの **float 実体** (plans/active/condensation-float-speedup.md §4.2-2, methods/condensation.md 実装 §9)。
// condensationSource_d.cuh の double 関数 (現行、ビット不変で残す) を、物性表 (condensationTables_d.cuh) と float 演算で写したもの。
// 設計上の差 (意図的):
//   * 核生成率は対数空間で組み、指数化前に ln J_max=ln(1e35) で上限を掛ける (本体・摂動評価とも)。m³=2.7e-77 が float 範囲外のため。
//   * 二温度 (twoTemp) は float 実体に無い (wrapper が double kernel を選ぶ)。平衡形 (condEquilibrium) も同様。
//   * 接尾辞なしリテラルを書かない (double 昇格を防ぐ; cuobjdump -sass で FP64 命令 0 を監査)。
#include "condensationSource_d.cuh"
#include "condensationTables_d.cuh"

#define COND_KB_F   1.380649e-23f
#define COND_RU_F   8.314462618f
#define COND_PI_F   3.14159265f
#define COND_LN_JMAX_F 80.5904783f   // ln(1e35)  (< ln FLT_MAX = 88.72)
#define COND_LN_JMIN_F (-80.0f)      // これ未満は J=0 (double 版の expo<-700 と同じ「物理的ゼロ」)
#define COND_RNUC_FAC_F 1.01f

struct CondNucCarrierF { float a_v; float carrierSum; float cvv_tilde; };

__host__ __device__ inline float cond_kantrowitz_theta_f(
    const CondSpeciesPropsF& cp, const CondTablesF& tb, float T, float lnS, int mode, float gamma_gas, const CondNucCarrierF* car)
{
    const float b = cond_tab_latent_f(tb, T) / (cp.R*T);
    if (mode == 1) {
        const float theta = (2.0f*(gamma_gas-1.0f)/(gamma_gas+1.0f)) * b * (b - 0.5f);
        return (theta > 0.0f) ? theta : 0.0f;
    }
    if (mode >= 2) {
        const float cvv = car ? car->cvv_tilde : cp.cv*cp.M/COND_RU_F;
        const float av  = car ? car->a_v : 1.0f;
        const float cs  = car ? car->carrierSum : 0.0f;
        const float qhat = b - 0.5f - ((mode == 3) ? lnS : 0.0f);
        const float den  = av*(cvv + 0.5f) + cs;
        if (!(den > 0.0f) || !(av > 0.0f)) return 0.0f;
        return av*qhat*qhat/den;
    }
    return 0.0f;
}

// 核生成 (対数空間 CNT)。J [1/(m³ s)] は ln J_max で上限、ln J < -80 で 0。
__host__ __device__ inline void cond_nucleation_f(
    const CondSpeciesPropsF& cp, const CondTablesF& tb, float T, float p_v, float rho_v,
    float* J_out, float* rstar_out,
    int kantrowitz = 0, float gamma_gas = 1.4f, const CondNucCarrierF* car = nullptr)
{
    *J_out = 0.0f; *rstar_out = 0.0f;
    if (!(p_v > 0.0f) || !(rho_v > 0.0f)) return;
    const float lnS = logf(p_v) - cond_tab_lnpsat_f(tb, T);
    if (!(lnS > 0.0f)) return;
    const float R     = cp.R;
    const float sigma = cond_tab_sigma_f(tb, T);
    const float rho_l = cond_tab_rhol_f(tb, T);
    if (!(sigma > 0.0f) || !(rho_l > 0.0f)) return;          // T_c 直下 σ=0: double 版も J=0, r*=0
    const float rstar = 2.0f*sigma / (rho_l*R*T*lnS);
    const float dG    = (4.0f/3.0f)*COND_PI_F*rstar*rstar*sigma;
    float lnJ = 0.5f*logf(2.0f*sigma/COND_PI_F) - 1.5f*cp.lnM + 2.0f*logf(rho_v) - logf(rho_l) - dG/(COND_KB_F*T);
    if (cp.model == COND_MODEL_N2) lnJ += -55.0f + 4270.0f/T;   // Iland 経験補正
    if (kantrowitz) lnJ -= log1pf(cond_kantrowitz_theta_f(cp, tb, T, lnS, kantrowitz, gamma_gas, car));
    if (lnJ > COND_LN_JMAX_F) lnJ = COND_LN_JMAX_F;
    *J_out = (lnJ > COND_LN_JMIN_F) ? expf(lnJ) : 0.0f;
    *rstar_out = rstar;
}

__host__ __device__ inline void cond_vapor_state_f(
    int carrier, float rod, float Pd, float Td, float g, float Yw, float Rw, float* pv, float* rho_v)
{
    if (carrier) {
        float yv = Yw - g; if (yv < 0.0f) yv = 0.0f;
        *rho_v = rod*yv;
        *pv    = rod*yv*Rw*Td;
    } else {
        float omg = 1.0f - g; if (omg < 0.0f) omg = 0.0f;
        *rho_v = rod*omg;
        *pv    = Pd;
    }
}

__host__ __device__ inline float cond_mean_free_path_f(const CondTablesF& tb, float T, float p, float R)
{
    const float pf = (p > 1.0f) ? p : 1.0f;
    return (cond_tab_mugas_f(tb, T)/pf) * sqrtf(0.5f*COND_PI_F*R*T);
}

// 成長率 dr/dt [m/s] (一温度のみ; twoTemp は double kernel)。
__host__ __device__ inline float cond_growth_f(
    const CondSpeciesPropsF& cp, const CondTablesF& tb, float T, float p_v, float r_bar, float rstar,
    int growthModel = 0, float p_gas = -1.0f, float gyarC = 3.18f)
{
    if (r_bar <= 0.0f) return 0.0f;
    const float R     = cp.R;
    const float lnps  = cond_tab_lnpsat_f(tb, T);
    const float rho_l = cond_tab_rhol_f(tb, T);
    const float pK    = (p_gas > 0.0f) ? p_gas : p_v;
    const float lnpv  = (p_v > 0.0f) ? logf(p_v) : -1.0e30f;
    if (growthModel == 1) {
        if (!(lnpv > lnps)) return 0.0f;
        const float driving = lnpv - lnps;
        const float L   = cond_tab_latent_f(tb, T);
        const float k   = cond_tab_kgas_f(tb, T);
        const float lam = cond_mean_free_path_f(tb, T, pK, R);
        const float Kn  = lam/(2.0f*r_bar);
        const float fac = (1.0f - rstar/r_bar) / (r_bar*(1.0f + gyarC*Kn));
        return (1.0f/rho_l)*fac*(k*R*T*T/(L*L))*driving;
    }
    if (cp.model == COND_MODEL_H2O) {
        const float sigma = cond_tab_sigma_f(tb, T);
        const float alpha = 1.0f;
        const float pd = expf(lnps + 2.0f*sigma/(rho_l*R*T*r_bar));   // Kelvin: p_sat·exp(2σ/(ρ_l R T r))
        return (alpha/rho_l)*(p_v - pd)/sqrtf(2.0f*COND_PI_F*R*T);
    } else {
        const float driving = lnpv - lnps;
        const float L   = cond_tab_latent_f(tb, T);
        const float k   = cond_tab_kgas_f(tb, T);
        const float lam = cond_mean_free_path_f(tb, T, pK, R);
        const float Kn  = lam/(2.0f*r_bar);
        const float fFS = (1.0f+2.0f*Kn)/(r_bar*(1.0f+3.42f*Kn+5.32f*Kn*Kn)) * (1.0f - rstar/r_bar);
        return (1.0f/rho_l)*fFS*(k*R*T*T/(L*L))*driving;
    }
}

__host__ __device__ inline void cond_source_vector_f(
    const CondSpeciesPropsF& cp, const CondTablesF& tb, float T, float p_v, float rho_v,
    float roQ0, float roQ1, float roQ2,
    float* SQ0, float* SQ1, float* SQ2, float* Sg,
    int kantrowitz = 0, int growthModel = 0, float gamma_gas = 1.4f, float p_gas = -1.0f,
    float gyarC = 3.18f, const CondNucCarrierF* car = nullptr)
{
    if (rho_v < 0.0f) rho_v = 0.0f;
    float J, rstar;
    cond_nucleation_f(cp, tb, T, p_v, rho_v, &J, &rstar, kantrowitz, gamma_gas, car);
    const float r_bar = (roQ0 > 1.0e-30f) ? (roQ1/roQ0) : rstar;
    float drdt = 0.0f;
    if (roQ0 > 1.0e-30f && rstar > 0.0f && r_bar > rstar) {
        drdt = cond_growth_f(cp, tb, T, p_v, r_bar, rstar, growthModel, p_gas, gyarC);
        if (drdt < 0.0f) drdt = 0.0f;
    }
    const float rho_l = cond_tab_rhol_f(tb, T);
    const float r_nuc = COND_RNUC_FAC_F*rstar;
    *SQ0 = J;
    *SQ1 = J*r_nuc + roQ0*drdt;
    *SQ2 = J*r_nuc*r_nuc + 2.0f*roQ1*drdt;
    *Sg  = (4.0f/3.0f)*COND_PI_F*rho_l*(J*r_nuc*r_nuc*r_nuc + 3.0f*roQ2*drdt);
}

// 蒸発率 dr/dt <= 0 (r で評価)。
__host__ __device__ inline float cond_evap_rate_f(
    const CondSpeciesPropsF& cp, const CondTablesF& tb, float T, float p_v, float r,
    int growthModel, float p_gas, float gyarC, int kelvin)
{
    if (r <= 0.0f) return 0.0f;
    const float R     = cp.R;
    const float lnps  = cond_tab_lnpsat_f(tb, T);
    const float rho_l = cond_tab_rhol_f(tb, T);
    float Ke_r = 0.0f;
    if (kelvin) {
        const float sigma = cond_tab_sigma_f(tb, T);
        Ke_r = 2.0f*sigma/(rho_l*R*T*r);
        if (Ke_r > 5.0f) Ke_r = 5.0f;
    }
    const float lnpd = lnps + Ke_r;
    const float pd   = expf(lnpd);
    if (p_v >= pd) return 0.0f;
    if (growthModel == 0 && cp.model == COND_MODEL_H2O) {
        const float alpha = 1.0f;
        return (alpha/rho_l)*(p_v - pd)/sqrtf(2.0f*COND_PI_F*R*T);
    }
    const float driving = ((p_v > 0.0f) ? logf(p_v) : -690.0f) - lnpd;   // < 0 (p_v=0: double は ln(1e-300 相当) → 同じ「大きな負」)
    const float L   = cond_tab_latent_f(tb, T);
    const float k   = cond_tab_kgas_f(tb, T);
    const float pK  = (p_gas > 0.0f) ? p_gas : p_v;
    const float lam = cond_mean_free_path_f(tb, T, pK, R);
    const float Kn  = lam/(2.0f*r);
    float fac;
    if (growthModel == 1) fac = 1.0f/(r*(1.0f + gyarC*Kn));
    else                  fac = (1.0f+2.0f*Kn)/(r*(1.0f+3.42f*Kn+5.32f*Kn*Kn));
    return (1.0f/rho_l)*fac*(k*R*T*T/(L*L))*driving;
}

__host__ __device__ inline float cond_evap_source_f(
    const CondSpeciesPropsF& cp, const CondTablesF& tb, float T, float p_v, float rod, float g,
    float q0, float q1, float q2, float dt,
    float rmin, float lam_min, float dg_max, float dT_max, float cvg,
    int growthModel, float p_gas, float gyarC, int kelvin,
    float* SQ0, float* SQ1, float* SQ2, float* Sg, float* r30_out, float* drdt_out)
{
    // λ=r_new/r30 は 1 に近い (|λ−1| ~ dr/dt·dt/r30 ~ 1e-8) ので float では λ を作らず δ=λ−1 で組む
    // (λ を作ると 1−1e-8 が 1.0f に丸まり S≡0 になる)。λ³−1 = δ(3+3δ+δ²), λ²−1 = δ(2+δ)。
    *SQ0 = 0.0f; *SQ1 = 0.0f; *SQ2 = 0.0f; *Sg = 0.0f; *r30_out = 0.0f; *drdt_out = 0.0f;
    if (g <= 0.0f || dt <= 0.0f) return 1.0f;
    const float lnps = cond_tab_lnpsat_f(tb, T);
    if (p_v > 0.0f && logf(p_v) > lnps) return 1.0f;      // 過飽和: 蒸発分岐ではない
    const float rho_l = cond_tab_rhol_f(tb, T);
    float delta;                                          // δ = λ − 1 (≤ 0)
    bool vanish = false;
    if (q0 <= 1.0e-30f) {
        vanish = true;
    } else {
        const float r30 = cbrtf(g/((4.0f/3.0f)*COND_PI_F*rho_l*q0/rod));
        *r30_out = r30;
        if (r30 < 2.0f*rmin) {
            vanish = true;
        } else {
            const float drdt = cond_evap_rate_f(cp, tb, T, p_v, r30, growthModel, p_gas, gyarC, kelvin);
            *drdt_out = drdt;
            delta = drdt*dt/r30;
            if (delta < lam_min - 1.0f) delta = lam_min - 1.0f;
        }
    }
    if (vanish) delta = -1.0f;
    {
        // Δg 律速と潜熱冷却律速: λ³ ≥ 1 − x, x = min(dg_max/g, dT_max c_v/(g L)) → δ ≥ (1−x)^{1/3} − 1 = expm1(log1p(−x)/3)
        float x = dg_max/g;
        const float xT = dT_max*cvg/(g*cond_tab_latent_f(tb, T));
        if (xT < x) x = xT;
        if (x < 1.0f) {
            const float delta_T = expm1f(log1pf(-x)/3.0f);
            if (delta < delta_T) delta = delta_T;
        }
        if (delta > 0.0f) delta = 0.0f;
    }
    const float lam = 1.0f + delta;
    *SQ0 = (lam <= 0.0f) ? (-q0/dt) : 0.0f;
    *SQ1 = delta*q1/dt;
    *SQ2 = delta*(2.0f + delta)*q2/dt;
    *Sg  = delta*(3.0f + delta*(3.0f + delta))*rod*g/dt;
    return lam;
}

// 飽和温度 T_sat(p_v): 表の ln p_sat とその解析微分で Newton (前 step の値を warm start に)。
__host__ __device__ inline float cond_Tsat_f(const CondTablesF& tb, float pv, float T_guess)
{
    if (!(pv > 1.0e-6f)) return 0.0f;
    const float lnpv = logf(pv);
    float T = (T_guess > 50.0f && T_guess < 1000.0f) ? T_guess : 250.0f;
    #pragma unroll 1
    for (int it = 0; it < 25; ++it) {
        float d;
        const float f = cond_tab_lnpsat_f(tb, T, &d) - lnpv;
        if (d < 1.0e-6f) d = 1.0e-6f;
        float dT = f/d;
        if (dT >  0.3f*T) dT =  0.3f*T;
        if (dT < -0.3f*T) dT = -0.3f*T;
        T -= dT;
        if (T < 50.0f) T = 50.0f;
        if (T > 1000.0f) T = 1000.0f;
        if (fabsf(dT) < 1.0e-3f) break;
    }
    return T;
}

// CPG 二相の面全エンタルピー (float 表版; double は condensationEOS_d.cuh cond_face_h_cpg と同じ規約: R_eff 非正/非有限は乾き面へ退避)。
__host__ __device__ inline float cond_face_h_cpg_f(const CondTablesF& tb, float cp_gas, float R_gas, float R_w,
                                                   float g_f, float p_f, float rho_f, float ek)
{
    float Reff = R_gas - g_f*R_w;
    if (!(Reff > 0.0f) || !isfinite(Reff)) { Reff = R_gas; g_f = 0.0f; }
    const float Tf = p_f/(rho_f*Reff);
    return cp_gas*Tf - g_f*cond_tab_latent_f(tb, Tf) + ek;
}
