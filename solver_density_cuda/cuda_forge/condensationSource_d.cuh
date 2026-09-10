#pragma once

#include "condensationProperties_d.cuh"

// 非平衡凝縮: 核生成 (CNT+Iland) と成長 (Goodheart) のソース項 device 関数 (Phase 2)。
// methods/condensation.md 2-3 節。内部 double。現状 N2 のみ (n2_* 物性)。
//
// モーメント保存ベクトル後半の相変化ソース (論文 Eq.1):
//   d(ρQ0)/dt = J
//   d(ρQ1)/dt = J r*   + ρQ0 dr/dt
//   d(ρQ2)/dt = J r*^2 + 2 ρQ1 dr/dt
//   d(ρg)/dt  = (4/3)π ρ_l (J r*^3 + 3 ρQ2 dr/dt)
// J: 核生成率 [1/(m^3 s)], r*: 臨界半径 [m], dr/dt: 平均成長率 [m/s] (平均半径 r̄=Q1/Q0 で評価)。
// 一温度: T_v=T_d=T (核生成・成長とも気相 T)。

#define COND_KB 1.380649e-23      // Boltzmann [J/K]
#define COND_NA 6.02214076e23     // Avogadro [1/mol]
#define COND_PI 3.141592653589793

// 核生成は臨界半径 r* ちょうどで生むと成長則の (1-r*/r) が 0 (不安定平衡) になり、
// 平均半径 r̄ が r* に張り付いて成長が起動しない。わずかに超臨界 r_nuc = COND_RNUC_FAC*r*
// で生むことで r̄>r* となり成長が立ち上がる (Q1/Q2/g の核生成項に適用)。
#define COND_RNUC_FAC 1.01

// 核生成 (CNT × 種ごと補正)。J [1/(m^3 s)], r* [m]。過飽和でなければ 0。cprops で N2/H2O を切替。
//   J_CNT = √(2σ/(π m^3))·(ρ_v²/ρ_l)·exp(-ΔG*/(k_B T)),  ΔG*=(4/3)π r*²σ,  r*=2σ/(ρ_l R T ln S)。
//   N2: ×exp(A+B/T) [Iland]。H2O: 等温 CNT (キャリア N2 がクラスタを熱平衡化し Kantrowitz 非等温抑制は
//        ~1 に量子化されるため、希薄水-N2 では等温近似。過抑制を避ける; 必要なら carrier-Kantrowitz を後段)。
// 核生成 (CNT × 種ごと補正 × 任意 Kantrowitz 非等温補正)。
//   kantrowitz!=0 のとき J を 1/(1+θ) 倍する (Feder/Kantrowitz 非等温補正, 純蒸気形):
//     θ = 2(γ_v-1)/(γ_v+1) · b(b-1/2),  b = L/(R_v T)  (γ_v = 凝縮種 (蒸気) 自身の比熱比。gamma_gas 引数に渡す。
//     2(γ_v−1)/(γ_v+1)=R_v/(c_v,v+R_v/2) は蒸気分子との衝突で持ち去れるエネルギー揺らぎの尺度なので蒸気の熱容量を使う)。
//   純蒸気では θ が大きく J を桁で抑える。キャリア気体 (N2) 中ではキャリア衝突が潜熱を奪い θ はさらに小さい (Feder 拡張, 未実装)。
//   2026-09-10 まではセル気相混合の γ (carrier では ≈1.40) を渡していた (旧挙動は condKantrowitzGammaMode=1)。
__host__ __device__ inline void cond_nucleation(
    const CondSpeciesProps& cp, double T, double p_v, double rho_v,
    double* J_out, double* rstar_out,
    int kantrowitz = 0, double gamma_gas = 1.4)
{
    const double psat = cond_psat(cp, T);
    const double S = p_v / (psat > 1.0e-300 ? psat : 1.0e-300);
    if (S <= 1.0 || rho_v <= 0.0) { *J_out = 0.0; *rstar_out = 0.0; return; }
    const double R     = cp.R;
    const double lnS   = log(S);
    const double sigma = cond_sigma(cp, T);
    const double rho_l = cond_rho_cond(cp, T);
    const double rstar = 2.0*sigma / (rho_l*R*T*lnS);
    const double dG    = (4.0/3.0)*COND_PI*rstar*rstar*sigma;
    const double m     = cp.M / COND_NA;
    const double Kcnt  = sqrt(2.0*sigma/(COND_PI*m*m*m)) * (rho_v*rho_v/rho_l);
    double corr = 1.0;
    if (cp.model == COND_MODEL_N2) corr = exp(-55.0 + 4270.0/T); // Iland 経験補正
    // H2O: corr=1 (carrier-thermalized 等温 CNT)
    if (kantrowitz) {
        const double b = cond_latent(cp, T) / (R*T);
        const double theta = (2.0*(gamma_gas-1.0)/(gamma_gas+1.0)) * b * (b - 0.5);
        corr /= (1.0 + (theta > 0.0 ? theta : 0.0));   // 非等温抑制 (θ<0 は掛けない)
    }
    const double expo = -dG/(COND_KB*T);
    const double J = (expo > -700.0) ? Kcnt*exp(expo)*corr : 0.0;
    *J_out = J;
    *rstar_out = rstar;
}

// 成長率 dr/dt [m/s]。一温度 T_d=T。
//   既定 (growthModel=0): N2=Goodheart, H2O=Hertz-Knudsen(自由分子, 質量輸送律速)。
//   growthModel=1 (Gyarmathy): 熱伝導律速の Gyarmathy(1982) 式 (極超音速ノズル凝縮で標準)。
//     dr/dt = (k R T²/L²) ln(S) /ρ_l · (1-r*/r)/(r(1+3.18 Kn)),  Kn=l/(2r) (キャリア気体の l)。
//     N2 Goodheart と前因子 (k R T² ln S /ρ_l L²) を共有し、Knudsen 内挿のみ 1/(1+3.18Kn) に差し替えた形。
//   p_gas: Kn 用の気相全圧 (<0 で p_v にフォールバック=純蒸気)。k_gas は n2_kgas(T) (キャリア N2 熱伝導率)。
// twoTemp!=0 のとき、Hertz-Knudsen 経路で液滴温度 T_d を準定常エネルギーバランス (Hill) で解き、
// 平衡蒸気圧 p_d を T_d で評価する (自己加熱で実効過飽和が下がり成長が遅くなる)。
//   潜熱解放 L·j(T_d) = 気相への熱伝導 h·(T_d−T_g),  j=α(p_v−p_d(T_d))/√(2πR T_g)
//   p_d(T_d)=p_sat(T_d)·exp(2σ/(ρ_l R T_d r)),  h=λ_g/(r(1+3.18Kn))。
// Gyarmathy(growthModel=1, 熱伝導律速)は元来 T_d≈T_s を内包するため twoTemp は適用しない。
__host__ __device__ inline double cond_growth(
    const CondSpeciesProps& cp, double T, double p_v, double r_bar, double rstar,
    int growthModel = 0, double p_gas = -1.0, double gyarC = 3.18, int twoTemp = 0)
{
    if (r_bar <= 0.0) return 0.0;
    const double R    = cp.R;
    const double psat = cond_psat(cp, T);
    const double rho_l= cond_rho_cond(cp, T);
    const double pK   = (p_gas > 0.0) ? p_gas : p_v;   // Kn 用は全圧 (carrier)、純蒸気では p_v
    if (growthModel == 1) {
        // Gyarmathy 熱伝導律速 (種共通)。過飽和 → 過冷却 ΔT_sub=(R T²/L)ln S を駆動力に。
        // Knudsen 補正 1/(1+gyarC·Kn): gyarC は Gyarmathy 標準 3.18 (config condGyarmathyC で可変)。
        if (p_v <= psat) return 0.0;
        const double driving = log(p_v / (psat > 1.0e-300 ? psat : 1.0e-300));
        const double L   = cond_latent(cp, T);
        const double k   = n2_kgas(T);                 // キャリア (N2) 熱伝導率
        const double lam = cond_mean_free_path(T, pK, R);
        const double Kn  = lam/(2.0*r_bar);
        const double fac = (1.0 - rstar/r_bar) / (r_bar*(1.0 + gyarC*Kn));
        return (1.0/rho_l)*fac*(k*R*T*T/(L*L))*driving;
    }
    if (cp.model == COND_MODEL_H2O) {
        // Hertz-Knudsen 質量律速。一温度 (twoTemp=0) は T_d=T_g、二温度は T_d を Newton で解く。
        const double sigma = cond_sigma(cp, T);
        const double alpha = 1.0; // 質量適応係数
        double Td = T;
        if (twoTemp) {
            const double L  = cond_latent(cp, T);
            const double kg = n2_kgas(T);                       // キャリア (N2) 熱伝導率
            const double lam= cond_mean_free_path(T, pK, R);
            const double Kn = lam/(2.0*r_bar);
            const double h  = kg/(r_bar*(1.0 + 3.18*Kn));       // 熱伝達係数 [W/m^2/K]
            const double c  = alpha/sqrt(2.0*COND_PI*R*T);      // 質量フラックス係数
            for (int it = 0; it < 25; ++it) {
                const double ps  = cond_psat(cp, Td);
                const double kel = exp(2.0*sigma/(rho_l*R*Td*r_bar));
                const double pd  = ps*kel;
                const double F   = L*c*(p_v - pd) - h*(Td - T); // =0 を解く
                const double dps = (cond_psat(cp, Td+0.1) - cond_psat(cp, Td-0.1))/0.2;
                const double dkel= kel*(-2.0*sigma/(rho_l*R*r_bar))/(Td*Td);
                const double dpd = dps*kel + ps*dkel;
                const double dF  = -L*c*dpd - h;
                double dTd = F/dF;
                if (dTd >  0.5*Td) dTd =  0.5*Td;
                if (dTd < -0.5*Td) dTd = -0.5*Td;
                Td -= dTd;
                if (Td < T) Td = T;                              // 凝縮中は T_d ≥ T_g
                if (fabs(dTd) < 1.0e-4) break;
            }
        }
        const double ps = cond_psat(cp, Td);
        const double pd = ps*exp(2.0*sigma/(rho_l*R*Td*r_bar));  // 液滴温度での平衡蒸気圧
        return (alpha/rho_l)*(p_v - pd)/sqrt(2.0*COND_PI*R*T);
    } else {
        const double driving = log(p_v / (psat > 1.0e-300 ? psat : 1.0e-300));
        const double L   = cond_latent(cp, T);
        const double k   = n2_kgas(T);
        const double lam = cond_mean_free_path(T, pK, R);
        const double Kn  = lam/(2.0*r_bar);
        const double fFS = (1.0+2.0*Kn)/(r_bar*(1.0+3.42*Kn+5.32*Kn*Kn)) * (1.0 - rstar/r_bar);
        return (1.0/rho_l)*fFS*(k*R*T*T/(L*L))*driving;
    }
}

// 相変化ソースベクトル S=(S_Q0,S_Q1,S_Q2,S_g)。モーメントは保存量 (ρQn)。
//   p_v: 凝縮種の蒸気分圧 (pure=P, carrier=ρ(Y-g)R_w T)。rho_v: 蒸気密度。cprops で種を切替。
__host__ __device__ inline void cond_source_vector(
    const CondSpeciesProps& cp, double T, double p_v, double rho_v,
    double roQ0, double roQ1, double roQ2,
    double* SQ0, double* SQ1, double* SQ2, double* Sg,
    int kantrowitz = 0, int growthModel = 0, double gamma_gas = 1.4, double p_gas = -1.0,
    double gyarC = 3.18, int twoTemp = 0)
{
    if (rho_v < 0.0) rho_v = 0.0;
    double J, rstar;
    cond_nucleation(cp, T, p_v, rho_v, &J, &rstar, kantrowitz, gamma_gas);
    const double r_bar = (roQ0 > 1.0e-30) ? (roQ1/roQ0) : rstar;
    // 成長は r̄>r* のときのみ・蒸発(<0)はしない (本体 kernel と整合; 亜臨界での発散を防ぐ)。
    double drdt = 0.0;
    if (roQ0 > 1.0e-30 && rstar > 0.0 && r_bar > rstar) {
        drdt = cond_growth(cp, T, p_v, r_bar, rstar, growthModel, p_gas, gyarC, twoTemp);
        if (drdt < 0.0) drdt = 0.0;
    }
    const double rho_l = cond_rho_cond(cp, T);
    const double r_nuc = COND_RNUC_FAC*rstar;   // わずかに超臨界で生み成長を起動
    *SQ0 = J;
    *SQ1 = J*r_nuc + roQ0*drdt;
    *SQ2 = J*r_nuc*r_nuc + 2.0*roQ1*drdt;
    *Sg  = (4.0/3.0)*COND_PI*rho_l*(J*r_nuc*r_nuc*r_nuc + 3.0*roQ2*drdt);
}

// =====================================================================================
// 蒸発 (plans/accepted/condensation-evaporation.md)。S=p_v/p_sat<=1 のセルで負成長・液滴消滅。
//   統一駆動力: 成長則の (1-r*/r) ln S ≡ ln S - K_e/r = ln(p_v/p_d(r)),  p_d=p_sat exp(K_e/r),
//   K_e=2σ/(ρ_l R T)。S<=1 では r* が定義されないので r* を経由せず p_d(r) で駆動力を評価する。
//   既定 (kelvin=0) は Kelvin 項を落とし平面 p_sat で駆動 (小半径ほど速く蒸発する正帰還を断つ。
//   質量収支には効かない: r0 から r まで縮んだ残量は (r/r0)^3)。
// =====================================================================================

// 蒸発率 dr/dt <= 0 [m/s] (半径 r で評価; 呼び出し側は体積平均半径 r30 を渡す)。p_v >= p_d なら 0。
__host__ __device__ inline double cond_evap_rate(
    const CondSpeciesProps& cp, double T, double p_v, double r,
    int growthModel, double p_gas, double gyarC, int kelvin)
{
    if (r <= 0.0) return 0.0;
    const double R     = cp.R;
    const double psat  = cond_psat(cp, T);
    const double rho_l = cond_rho_cond(cp, T);
    double Ke_r = 0.0;                                   // Kelvin 指数 K_e/r (kelvin=0 で 0)
    if (kelvin) {
        const double sigma = cond_sigma(cp, T);
        Ke_r = 2.0*sigma/(rho_l*R*T*r);
        if (Ke_r > 5.0) Ke_r = 5.0;                      // exp 発散ガード (r が極小のとき)
    }
    const double pd = psat*exp(Ke_r);
    if (p_v >= pd) return 0.0;
    if (growthModel == 0 && cp.model == COND_MODEL_H2O) {
        // Hertz-Knudsen 質量律速 (一温度)。
        const double alpha = 1.0;
        return (alpha/rho_l)*(p_v - pd)/sqrt(2.0*COND_PI*R*T);
    }
    // Goodheart (N2 既定) / Gyarmathy: 前因子 kRT²/(ρ_l L²) × 駆動力 ln(p_v/p_d) × Kn 補正。
    const double driving = log(p_v/(pd > 1.0e-300 ? pd : 1.0e-300));   // < 0
    const double L   = cond_latent(cp, T);
    const double k   = n2_kgas(T);
    const double pK  = (p_gas > 0.0) ? p_gas : p_v;
    const double lam = cond_mean_free_path(T, pK, R);
    const double Kn  = lam/(2.0*r);
    double fac;
    if (growthModel == 1) fac = 1.0/(r*(1.0 + gyarC*Kn));
    else                  fac = (1.0+2.0*Kn)/(r*(1.0+3.42*Kn+5.32*Kn*Kn));
    return (1.0/rho_l)*fac*(k*R*T*T/(L*L))*driving;
}

// 蒸発の 1 step 更新を「縮小比 λ=r_new/r30」で表し、monodisperse 整合ソース S_φ=(φ_new-φ)/dt を返す。
//   Q1→λQ1, Q2→λ²Q2, g→λ³g (Q0 は不変)。λ は次で下から律速:
//     λ >= lam_min (半径半減/step), g(1-λ³) <= dg_max, 潜熱冷却 g(1-λ³)L/c_v <= dT_max
//   r30 < 2 rmin (= λ_min r30 < rmin) なら中間状態を作らず λ=0 (全モーメント消滅、Q0 も 0)。
//   戻り値: λ。*r30_out に体積平均半径。S<=1 でないセル・液相なしセルは λ=1, S=0。
//   注意: q0,q1,q2 は保存量 ρQn、g は質量分率。S_g は ρg 単位。
__host__ __device__ inline double cond_evap_source(
    const CondSpeciesProps& cp, double T, double p_v, double rod, double g,
    double q0, double q1, double q2, double dt,
    double rmin, double lam_min, double dg_max, double dT_max, double cvg,
    int growthModel, double p_gas, double gyarC, int kelvin,
    double* SQ0, double* SQ1, double* SQ2, double* Sg, double* r30_out, double* drdt_out)
{
    *SQ0 = 0.0; *SQ1 = 0.0; *SQ2 = 0.0; *Sg = 0.0; *r30_out = 0.0; *drdt_out = 0.0;
    if (g <= 0.0 || dt <= 0.0) return 1.0;
    const double psat = cond_psat(cp, T);
    if (p_v > psat) return 1.0;                          // 過飽和: 蒸発分岐ではない
    const double rho_l = cond_rho_cond(cp, T);
    double lam;
    if (q0 <= 1.0e-30) {
        lam = 0.0;                                       // 液相はあるが液滴数 0 (不整合) → 消滅
    } else {
        const double r30 = cbrt(g/((4.0/3.0)*COND_PI*rho_l*q0/rod));  // q0/rod = Q0 [1/kg]
        *r30_out = r30;
        if (r30 < 2.0*rmin) {
            lam = 0.0;                                   // 次の半減で rmin を割る → 一括消滅
        } else {
            const double drdt = cond_evap_rate(cp, T, p_v, r30, growthModel, p_gas, gyarC, kelvin);
            *drdt_out = drdt;
            lam = 1.0 + drdt*dt/r30;
            if (lam < lam_min) lam = lam_min;
        }
    }
    // Δg 律速 (質量分率) と潜熱冷却律速: g(1-λ³) <= dg_max, g(1-λ³)L/c_v <= dT_max (λ=0 消滅にも適用)
    {
        double lam3_min = 1.0 - dg_max/g;
        const double lam3_T = 1.0 - dT_max*cvg/(g*cond_latent(cp, T));
        if (lam3_T > lam3_min) lam3_min = lam3_T;
        if (lam3_min > 0.0) {
            const double lam_T = cbrt(lam3_min);
            if (lam < lam_T) lam = lam_T;
        }
        if (lam > 1.0) lam = 1.0;
    }
    const double lam2 = lam*lam, lam3 = lam2*lam;
    *SQ0 = (lam <= 0.0) ? (-q0/dt) : 0.0;
    *SQ1 = (lam - 1.0)*q1/dt;
    *SQ2 = (lam2 - 1.0)*q2/dt;
    *Sg  = (lam3 - 1.0)*rod*g/dt;
    return lam;
}
