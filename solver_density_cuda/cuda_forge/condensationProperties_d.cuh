#pragma once

#include <cuda_runtime.h>
#include <math.h>   // host コンパイル時の pow/log/exp/sqrt (device は組込み)

// 非平衡凝縮: 凝縮種ごとの物性相関 (飽和蒸気圧・凝縮相密度・潜熱・表面張力)。
// 物性評価は exp/log で桁が飛ぶため内部は **double** で計算する (flow_float が float でも安全)。
//
// 【相の方針 (重要)】Lin 2014 のモデルは **液相 (過冷却=supercooled liquid)**。非平衡凝縮は
// 三重点以下でも準安定な過冷却液滴を作るため、本実装も**全域で液相フィットを使う** (固相昇華へは
// 切替えない)。r* = 2σ/(ρ_l RT ln(p_v/p_sat)) も液密度 ρ_l を使う。詳細は methods/condensation.md 8 節。
//
// 【低温外挿ガード (重要)】液相フィット (Jacobsen p_sat, 潜熱多項式) は ~40K 以下に外挿すると破綻する
// (p_sat が ~33K で非単調、潜熱が崩落)。有効域は概ね 45–126K。凝縮 ON では潜熱放出で活性域が ~40K 以上に
// 留まるが、過渡/dry セルが低温に落ちても破綻しないよう、物性評価温度を [T_PROP_FLOOR, Tc) にクランプする。
// クランプは「より低温=より低 p_sat」側に倒れる安全側 (過剰核生成を起こさない)。既知の近似として docs に明記。
//
// 【気相 thermo の制約】凝縮ケースの気相は **thermalMethod 0 (熱量的完全気体, cp/γ 一定)** を使うこと。
// NASA-9/CEA は ~200K 未満で無効 (出口 ~27K)。論文も calorically perfect gas。
//
// 種の切替は enum + 係数構造体 CondSpeciesProps で行う。N2 を最初に実装、H2O は係数差し替えで足す。

// 液相フィットの有効下限 (これ未満は外挿破綻するためクランプ)。
#define COND_T_PROP_FLOOR 45.0

enum CondPropModel {
    COND_MODEL_N2  = 0,
    COND_MODEL_H2O = 1,   // Phase 3
};

// 凝縮種の物性パラメータ (device へ value 渡しできる POD)。
struct CondSpeciesProps {
    int    model;   // CondPropModel
    double R;       // 蒸気の気体定数 [J/(kg·K)]
    double cv;      // 蒸気の定積比熱 [J/(kg·K)] (calorically perfect)
    double cp;      // 蒸気の定圧比熱 [J/(kg·K)]
    double Tt;      // 三重点温度 [K] (液/固 切替)
    double Tc;      // 臨界温度 [K]
    double M;       // 分子量 [kg/mol]
};

// N2 既定パラメータ。R=296.8, γ=1.4 → cv=R/(γ-1)=742, cp=γcv=1038.8。
__host__ __device__ inline CondSpeciesProps condProps_N2()
{
    CondSpeciesProps s;
    s.model = COND_MODEL_N2;
    s.R  = 296.8;
    s.cv = 742.0;
    s.cp = 1038.8;
    s.Tt = 63.15;
    s.Tc = 126.192;
    s.M  = 0.0280134;
    return s;
}

// 物性評価温度を液相フィット有効域 [T_PROP_FLOOR, Tc) にクランプ。
__host__ __device__ inline double cond_clamp_Tprop(double T, double Tc)
{
    double Tcl = T;
    if (Tcl < COND_T_PROP_FLOOR) Tcl = COND_T_PROP_FLOOR;
    const double Tmax = Tc - 0.5; // 臨界点直下で発散させない
    if (Tcl > Tmax) Tcl = Tmax;
    return Tcl;
}

// --- N2 潜熱 (蒸発) L(T) [J/kg] --- 式26 (4次多項式 MJ/kg)。psat の C-C 外挿でも使うため前方に置く。
__host__ __device__ inline double n2_latent(double T)
{
    const double Tcl = cond_clamp_Tprop(T, 126.192);
    const double p1 = -2.137e-8, p2 = 7.18e-6, p3 = -9.142e-4, p4 = 0.05069, p5 = -0.809;
    const double L = p1*Tcl*Tcl*Tcl*Tcl + p2*Tcl*Tcl*Tcl + p3*Tcl*Tcl + p4*Tcl + p5; // MJ/kg
    const double Lj = L * 1.0e6;
    return (Lj > 0.0) ? Lj : 0.0;
}

// Jacobsen 液飽和圧 (式22, atm→Pa)。有効域内の生評価。
__host__ __device__ inline double n2_psat_jacobsen(double Tcl)
{
    const double Tc = 126.192;
    const double n1 = 8394.409444, n2 = -1890.045259, n3 = -7.282229165;
    const double n4 = 0.01022850966, n5 = 5.556063825e-4, n6 = -5.944544662e-6;
    const double n7 = 2.715433932e-8, n8 = -4.879535904e-11, n9 = 509.5360824;
    const double dTc = Tc - Tcl;
    const double lnP = n1/Tcl + n2 + n3*Tcl + n4*pow(dTc, 1.95)
                     + n5*Tcl*Tcl*Tcl + n6*Tcl*Tcl*Tcl*Tcl + n7*Tcl*Tcl*Tcl*Tcl*Tcl
                     + n8*Tcl*Tcl*Tcl*Tcl*Tcl*Tcl + n9*log(Tcl);
    return exp(lnP) * 101325.0;
}

// --- N2 飽和蒸気圧 p_sat(T) [Pa] (過冷却液) ---
// 有効域 [T_switch, Tc) は Jacobsen (式22)。**T_switch 未満は Clausius–Clapeyron で物理外挿**:
//   p_sat(T) = p_sat(T_sw)·exp(-(L/R)(1/T - 1/T_sw))   (低温で単調に減少、過飽和を正しく与える)。
// 45K クランプ (cond_clamp_Tprop) では psat を凍結し過飽和 S=p/psat を潰してしまい核生成が起きないため、
// psat だけは C-C 外挿を使う (σ,ρ_l,L は緩変化なのでクランプのまま)。L_ref=L(T_sw)。
#define COND_PSAT_TSWITCH 50.0
__host__ __device__ inline double n2_psat(double T)
{
    const double Tc = 126.192;
    if (T >= COND_PSAT_TSWITCH) {
        double Tcl = (T < Tc - 0.5) ? T : (Tc - 0.5);  // 臨界直下のみ上クランプ
        return n2_psat_jacobsen(Tcl);
    }
    // C-C 外挿 (T < T_switch)
    const double Tsw  = COND_PSAT_TSWITCH;
    const double psw  = n2_psat_jacobsen(Tsw);
    const double Lref = n2_latent(Tsw);   // ~2.04e5 J/kg
    const double Rv   = 296.8;
    double Tlo = (T > 5.0) ? T : 5.0;     // 0 割回避
    return psw * exp(-(Lref/Rv)*(1.0/Tlo - 1.0/Tsw));
}

// --- N2 凝縮相 (液) 密度 ρ_l(T) [kg/m³] --- Nowak (式24)。
__host__ __device__ inline double n2_rho_cond(double T)
{
    const double Tc = 126.192, rhoc = 313.3;
    const double Tcl = cond_clamp_Tprop(T, Tc);
    double tau = 1.0 - Tcl/Tc;
    if (tau < 0.0) tau = 0.0;
    const double n1 = 1.48654237, n2 = -0.280476066, n3 = 0.0894143085, n4 = -0.119879866;
    const double lnr = n1*pow(tau, 0.3294) + n2*pow(tau, 4.0/6.0)
                     + n3*pow(tau, 16.0/6.0) + n4*pow(tau, 35.0/6.0);
    return rhoc * exp(lnr);
}

// --- N2 表面張力 (液) σ(T) [N/m] --- Stansfield (式28, dyn/cm)。
__host__ __device__ inline double n2_sigma(double T)
{
    const double Tc = 126.0, sig0 = 29.06; // dyn/cm
    const double Tcl = cond_clamp_Tprop(T, 126.192);
    double x = 1.0 - Tcl/Tc;
    if (x < 0.0) x = 0.0;
    return sig0 * pow(x, 1.247) * 1.0e-3; // dyn/cm -> N/m
}

// --- N2 気相輸送物性 (成長則の熱伝導 k と 平均自由行程 λ 用) ---
// 成長則 (Goodheart) は気相熱伝導 k(T) と Knudsen 数 Kn=λ/2r を要する。case/34 は非粘性 (viscMethod 0)
// で流れ場の μ/k を持たないため、凝縮ソース内で N2 の簡易モデルから別途評価する。低温外挿のため近似。

// N2 気相粘性 μ(T) [Pa·s] — Sutherland (μ0=1.663e-5 @273K, C=107K)。低温は外挿 (近似)。
__host__ __device__ inline double n2_mu_gas(double T)
{
    const double mu0 = 1.663e-5, T0 = 273.0, C = 107.0;
    double Tc = (T > 5.0) ? T : 5.0;
    return mu0 * (T0 + C)/(Tc + C) * pow(Tc/T0, 1.5);
}

// N2 気相熱伝導 k(T) [W/(m·K)] — k = μ cp/Pr (Pr=0.72, cp=1038.8)。
__host__ __device__ inline double n2_kgas(double T)
{
    const double cp = 1038.8, Pr = 0.72;
    return n2_mu_gas(T) * cp / Pr;
}

// 気相 平均自由行程 λ(T,p) [m] — kinetic theory λ = (μ/p)√(π R T/2)。R_N2=296.8。
__host__ __device__ inline double cond_mean_free_path(double T, double p, double R)
{
    const double pf = (p > 1.0) ? p : 1.0;
    return (n2_mu_gas(T)/pf) * sqrt(3.14159265358979*R*T/2.0);
}

// =====================================================================================
// H2O (水) 物性 — 過冷却液 (Wyslouzil ノズルは ~207K まで膨張、過冷却水滴)。
// psat は Murphy-Koop (2005) の過冷却液式 (123–332K で妥当)、他は標準フィット。
// =====================================================================================

// 水 飽和蒸気圧 p_sat(T) [Pa] — Murphy & Koop (2005) liquid (過冷却水, 123<T<332K)。
__host__ __device__ inline double h2o_psat(double T)
{
    double Tc = (T > 120.0) ? T : 120.0;
    const double lnp = 54.842763 - 6763.22/Tc - 4.210*log(Tc) + 0.000367*Tc
        + tanh(0.0415*(Tc - 218.8)) * (53.878 - 1331.22/Tc - 9.44523*log(Tc) + 0.014025*Tc);
    return exp(lnp); // Pa
}

// 水 液密度 ρ_l(T) [kg/m³] — 過冷却水の簡易フィット (200–300K で ~%)。
__host__ __device__ inline double h2o_rho_cond(double T)
{
    // 過冷却水: 277K で最大 ~1000、低温で僅かに低下。簡易: 線形近似 + フロア。
    double r = 1000.0 - 0.12*(277.0 - T); // ゆるい近似
    if (r < 920.0) r = 920.0;
    return r;
}

// 水 蒸発潜熱 L(T) [J/kg] = h_v(T) − h_l(T) を **CEA (NASA-9) の気相 H2O と液相 H2O(L) の全エンタルピー差**で作る
// (2026-08-18, ユーザ指示「CEA 式で L を逆算」)。
//   h_v: forge 種 DB と同じ CEA McBride–Gordon 2002 の H2O 200–1000 K 係数 (thermo_d.cu)。
//   h_l: CEA thermo.inp の H2O(L) 273.15–373.15 K 係数 (Cox 1989 / Haar 1984)。**273.15 K 未満は CEA に液相フィットが無く
//        (CEA は氷 H2O(cr))、多項式外挿は 250 K 以下で発散 (cp_l 230 K で 10 kJ/kgK, 200 K で 50 kJ/kgK) するので、
//        h_l を 273.15 K の値と勾配 cp_l(273.15)=4228 J/kgK で線形外挿**する (過冷却水の標準的な扱い、cp_l 一定)。
//   旧線形フィット 3.1485e6−2370 T は本構成と 250 K で 0.3 %、273 K で 0.05 % 一致していた (等価だったことの確認)。
//   L(250 K)=2.563 MJ/kg, L(273.15)=2.501 MJ/kg。氷 (昇華熱 2.835 MJ/kg) は使わない (飽和線も過冷却液 Murphy–Koop で統一)。
__host__ __device__ inline double h2o_nasa9_h_mass(const double* a, double T)
{
    // NASA-9: h/(RT) = -a1/T² + a2 lnT/T + a3 + a4 T/2 + a5 T²/3 + a6 T³/4 + a7 T⁴/5 + b1/T ; R_w = 8.314462618/0.0180153
    const double Rw = 8.314462618/0.0180153;
    const double hRT = -a[0]/(T*T) + a[1]*log(T)/T + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0
                     + a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[7]/T;
    return hRT*Rw*T;
}
__host__ __device__ inline double h2o_latent(double T)
{
    // 気相 H2O (CEA, 200–1000 K 区間; 1000 K 超は本用途で不要だが単調に外挿される)
    const double ag[8] = {-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03,
                          -7.342557370e-06, 4.955043490e-09,-1.336933246e-12,-3.303974310e+04};
    // 液相 H2O(L) (CEA, 273.15–373.15 K 区間)
    const double al[8] = { 1.326371304e+09,-2.448295388e+07, 1.879428776e+05,-7.678995050e+02,
                           1.761556813e+00,-2.151167128e-03, 1.092570813e-06, 1.101760476e+08};
    const double Tf = 273.15;
    double Tg = (T > COND_T_PROP_FLOOR) ? T : COND_T_PROP_FLOOR;
    if (Tg > 1000.0) Tg = 1000.0;
    const double hv = h2o_nasa9_h_mass(ag, Tg);
    double hl;
    if (Tg >= Tf) {
        const double Tl = (Tg < 373.15) ? Tg : 373.15;
        hl = h2o_nasa9_h_mass(al, Tl);
    } else {
        // 273.15 K 未満: cp_l(273.15) 一定で線形外挿 (過冷却水)
        const double h0  = h2o_nasa9_h_mass(al, Tf);
        const double cpl = (h2o_nasa9_h_mass(al, Tf + 0.5) - h2o_nasa9_h_mass(al, Tf - 0.5 + 1.0e-9))/1.0;  // ≈4228 J/kgK
        hl = h0 - cpl*(Tf - Tg);
    }
    double L = hv - hl;
    if (L < 1.5e6) L = 1.5e6;
    if (L > 3.5e6) L = 3.5e6;
    return L;
}

// 水 表面張力 σ(T) [N/m] — IAPWS 形を過冷却へ外挿 (Tc=647.096K)。
__host__ __device__ inline double h2o_sigma(double T)
{
    const double Tc = 647.096;
    double tau = (Tc - T)/Tc;
    if (tau < 0.0) tau = 0.0;
    double s = 0.2358 * pow(tau, 1.256) * (1.0 - 0.625*tau); // N/m
    if (s < 0.0) s = 0.0;
    return s;
}

// --- 種ディスパッチ (model で N2 / H2O を切替) ---
__host__ __device__ inline double cond_psat(const CondSpeciesProps& s, double T)
{
    return (s.model == COND_MODEL_H2O) ? h2o_psat(T) : n2_psat(T);
}
__host__ __device__ inline double cond_rho_cond(const CondSpeciesProps& s, double T)
{
    return (s.model == COND_MODEL_H2O) ? h2o_rho_cond(T) : n2_rho_cond(T);
}
__host__ __device__ inline double cond_latent(const CondSpeciesProps& s, double T)
{
    return (s.model == COND_MODEL_H2O) ? h2o_latent(T) : n2_latent(T);
}
__host__ __device__ inline double cond_sigma(const CondSpeciesProps& s, double T)
{
    return (s.model == COND_MODEL_H2O) ? h2o_sigma(T) : n2_sigma(T);
}

// 飽和温度 T_sat(p_v): ln p_sat(T) = ln p_v を Clausius-Clapeyron 勾配 d ln p/dT = L/(R T^2) の Newton で反転。
// 過冷却度 T_sat - T の診断用 (plans/accepted/condensation-equilibrium.md)。p_v <= 1e-6 Pa は 0 を返す。
__host__ __device__ inline double cond_Tsat(const CondSpeciesProps& s, double pv, double T_guess)
{
    if (!(pv > 1.0e-6)) return 0.0;
    const double lnpv = log(pv);
    double T = (T_guess > 50.0 && T_guess < 1000.0) ? T_guess : 250.0;
    #pragma unroll 1
    for (int it = 0; it < 25; ++it) {
        const double ps = cond_psat(s, T);
        if (!(ps > 1.0e-300)) { T *= 1.2; continue; }
        const double f  = log(ps) - lnpv;
        double dfdT = cond_latent(s, T)/(s.R*T*T);
        if (dfdT < 1.0e-6) dfdT = 1.0e-6;
        double dT = f/dfdT;
        if (dT >  0.3*T) dT =  0.3*T;
        if (dT < -0.3*T) dT = -0.3*T;
        T -= dT;
        if (T < 50.0) T = 50.0;
        if (T > 1000.0) T = 1000.0;
        if (fabs(dT) < 1.0e-4) break;
    }
    return T;
}

// H2O 既定パラメータ (キャリア+凝縮種)。R_H2O=461.5, M=0.0180153。
__host__ __device__ inline CondSpeciesProps condProps_H2O()
{
    CondSpeciesProps s;
    s.model = COND_MODEL_H2O;
    s.R  = 461.5;
    // 蒸気 cp/cv: NASA-9 (species_db.yaml H2O) の 200–300 K 平均 (cp 1851–1865, γ_v 1.329–1.332)。
    // Kantrowitz 非等温補正の γ_v=cp/cv=1.331 に使う (純蒸気形の係数 2(γ_v−1)/(γ_v+1)=R_v/(c_v,v+R_v/2))。
    // 気相 EOS は TP では NASA-9 を使うのでここは核生成補正専用 (plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)。
    s.cv = 1393.5;
    s.cp = 1855.0;
    s.Tt = 273.16;
    s.Tc = 647.096;
    s.M  = 0.0180153;
    return s;
}
