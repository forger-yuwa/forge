#pragma once

// =============================================================================
// thermo_d.cuh
//   多成分 thermally-perfect gas の熱力学コア。
//   - 比熱 cp(T) は NASA-9 係数多項式 (CEA 準拠, McBride-Gordon 2002, 2 温度域)。
//   - 混合則は ideal-gas mixing (質量重み cp/h、モル重み M)。
//   - エネルギーは NASA 絶対エンタルピー基準 (生成エンタルピー込み)。
//
//   設計メモ:
//   - flow_float は float (FP32)。NASA-9 多項式和・Newton 反転・混合和は
//     FP32 では破綻するため、ここでの内部計算は全て double で行う。
//   - 同じ関数を __host__ __device__ で提供し、初期条件 (host) と
//     カーネル (device) で同一の熱力学を保証する。
//   - 化学種データ (SpeciesThermo 配列) は device global memory に置き、
//     ポインタをカーネル引数で渡す (translation unit をまたぐ __constant__ は
//     -rdc を要するため避ける)。アクセスは全スレッド同係数読みで L2/読取専用
//     キャッシュに乗る。
// =============================================================================

#include "flowFormat.hpp"

#ifndef __CUDACC__
  #include <cmath>
  #define THERMO_HD inline
#else
  #define THERMO_HD __host__ __device__ inline
#endif

// 普遍気体定数 [J/(mol·K)]
#define THERMO_RU 8.314462618
// NASA-9 係数の標準状態圧力 p° [Pa] (CEA 規約 1 bar)。s = s0(T) − R ln(p/p°) の datum。
#define THERMO_P_STD 1.0e5
// コンパイル時の最大化学種数 (DB 容量の上限)
#define THERMO_MAX_SPECIES 16

// 1 化学種の NASA-9 係数 (2 温度域) と Lennard-Jones パラメータ。
// 係数配列 a[0..8] の規約:
//   cp/R = a0/T^2 + a1/T + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
//   h/RT = -a0/T^2 + a1 ln(T)/T + a2 + a3 T/2 + a4 T^2/3 + a5 T^3/4 + a6 T^4/5 + a7/T
//   s/R  = -a0/(2T^2) - a1/T + a2 ln(T) + a3 T + a4 T^2/2 + a5 T^3/3 + a6 T^4/4 + a8
struct SpeciesThermo {
    double MW;        // 分子量 [kg/mol]
    double sigma_LJ;  // Lennard-Jones 衝突直径 [Angstrom] (kinetic theory 輸送用)
    double eps_kB;    // Lennard-Jones ポテンシャル深さ ε/kB [K]
    double Tlo;       // 低温端 [K]
    double Tmid;      // 温度域境界 [K]
    double Thi;       // 高温端 [K]
    double low[9];    // Tlo  <= T <  Tmid の係数
    double high[9];   // Tmid <= T <= Thi  の係数
    double h_datum;   // sensible datum で係数から除いた絶対エンタルピー h_abs(Tref) [J/mol] (既定 0)。
                      // 反応流の反応熱 Q̇=−Σ(h_datum/W)ω と平衡定数 (H_abs=h+h_datum) が使う。
};

// -----------------------------------------------------------------------------
// 単一化学種の NASA-9 評価 (全て double)。範囲外は端でクランプし、
// エンタルピーは線形外挿 (h(T) ~ h(Tc) + cp(Tc)(T-Tc)) して衝撃波での暴走を防ぐ。
// -----------------------------------------------------------------------------

// 係数選択 (クランプ済み温度 Tc に対応する係数配列を返す)
THERMO_HD const double* thermo_pick_coeffs(const SpeciesThermo& sp, double Tc)
{
    return (Tc < sp.Tmid) ? sp.low : sp.high;
}

// cp_molar [J/(mol·K)] (クランプ温度で評価)
THERMO_HD double thermo_cp_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const double* a = thermo_pick_coeffs(sp, Tc);
    const double Ti  = 1.0/Tc;
    const double Ti2 = Ti*Ti;
    return THERMO_RU * ( a[0]*Ti2 + a[1]*Ti + a[2]
                       + a[3]*Tc + a[4]*Tc*Tc + a[5]*Tc*Tc*Tc + a[6]*Tc*Tc*Tc*Tc );
}

// h_molar [J/mol] (クランプ温度で評価, 絶対エンタルピー)
THERMO_HD double thermo_h_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const double* a = thermo_pick_coeffs(sp, Tc);
    const double Ti  = 1.0/Tc;
    const double Ti2 = Ti*Ti;
    const double lnT = log(Tc);
    // h/(R T) = -a0/T^2 + a1 ln(T)/T + a2 + a3 T/2 + a4 T^2/3 + a5 T^3/4 + a6 T^4/5 + a7/T
    const double hRT = -a[0]*Ti2 + a[1]*lnT*Ti + a[2]
                     + a[3]*Tc/2.0 + a[4]*Tc*Tc/3.0 + a[5]*Tc*Tc*Tc/4.0
                     + a[6]*Tc*Tc*Tc*Tc/5.0 + a[7]*Ti;
    return THERMO_RU * Tc * hRT;
}

// 範囲クランプ + 外挿込みの cp_molar [J/(mol·K)]
THERMO_HD double thermo_cp_molar(const SpeciesThermo& sp, double T)
{
    double Tc = T;
    if (Tc < sp.Tlo) Tc = sp.Tlo;
    if (Tc > sp.Thi) Tc = sp.Thi;
    return thermo_cp_molar_clamped(sp, Tc);
}

// 範囲クランプ + 線形外挿込みの h_molar [J/mol]
THERMO_HD double thermo_h_molar(const SpeciesThermo& sp, double T)
{
    if (T < sp.Tlo) {
        const double h0  = thermo_h_molar_clamped(sp, sp.Tlo);
        const double cp0 = thermo_cp_molar_clamped(sp, sp.Tlo);
        return h0 + cp0*(T - sp.Tlo);
    }
    if (T > sp.Thi) {
        const double h1  = thermo_h_molar_clamped(sp, sp.Thi);
        const double cp1 = thermo_cp_molar_clamped(sp, sp.Thi);
        return h1 + cp1*(T - sp.Thi);
    }
    return thermo_h_molar_clamped(sp, T);
}

// 単一化学種の比気体定数 R_s = Ru/MW [J/(kg·K)]
THERMO_HD double thermo_R_species(const SpeciesThermo& sp)
{
    return THERMO_RU / sp.MW;
}

// 単一化学種 質量基準: cp_s [J/(kg·K)], h_s [J/kg]
THERMO_HD double thermo_cp_mass(const SpeciesThermo& sp, double T)
{ return thermo_cp_molar(sp, T) / sp.MW; }

THERMO_HD double thermo_h_mass(const SpeciesThermo& sp, double T)
{ return thermo_h_molar(sp, T) / sp.MW; }

// -----------------------------------------------------------------------------
// 混合則 (質量分率 Y[0..n-1] を与える, ideal-gas mixing)
// -----------------------------------------------------------------------------

// 混合分子量の逆数 1/M_mix = Σ Y_s/MW_s [mol/kg]
THERMO_HD double thermo_inv_Mmix(const SpeciesThermo* sp, int n, const double* Y)
{
    double s = 0.0;
    for (int i=0;i<n;i++) s += Y[i]/sp[i].MW;
    return s;
}

// 混合比気体定数 R_mix = Ru·Σ Y_s/MW_s [J/(kg·K)]
THERMO_HD double thermo_R_mix(const SpeciesThermo* sp, int n, const double* Y)
{
    return THERMO_RU * thermo_inv_Mmix(sp, n, Y);
}

// 混合 cp [J/(kg·K)] (質量重み)
THERMO_HD double thermo_cp_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    double cp = 0.0;
    for (int i=0;i<n;i++) cp += Y[i]*thermo_cp_mass(sp[i], T);
    return cp;
}

// 混合絶対エンタルピー [J/kg] (質量重み)
THERMO_HD double thermo_h_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    double h = 0.0;
    for (int i=0;i<n;i++) h += Y[i]*thermo_h_mass(sp[i], T);
    return h;
}

// 混合内部エネルギー e_mix(T) = h_mix(T) - R_mix·T [J/kg]
THERMO_HD double thermo_e_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    return thermo_h_mix(sp, n, Y, T) - thermo_R_mix(sp, n, Y)*T;
}

// -----------------------------------------------------------------------------
// 性能最適化 (M6): cp と h を 1 スイープで評価する融合関数。
//   Newton 温度反転 (thermo_T_from_e/h) は毎反復 cp_mix と h_mix を別々の種ループ
//   で呼び、係数選択・T 冪・lnT を二重計算していた。両者を共有して半減させる。
//   値は thermo_cp_molar / thermo_h_molar (および _mass/_mix) と各々ビット同等
//   (同一演算順序: 範囲内は係数・T 冪・lnT を共有、範囲外は元関数へフォールバック)。
// -----------------------------------------------------------------------------
THERMO_HD void thermo_cph_molar(const SpeciesThermo& sp, double T, double* cp_out, double* h_out)
{
    if (T >= sp.Tlo && T <= sp.Thi) {
        const double* a = thermo_pick_coeffs(sp, T);
        const double Ti  = 1.0/T;
        const double Ti2 = Ti*Ti;
        const double lnT = log(T);
        const double T2 = T*T, T3 = T2*T, T4 = T3*T;
        *cp_out = THERMO_RU * ( a[0]*Ti2 + a[1]*Ti + a[2]
                              + a[3]*T + a[4]*T2 + a[5]*T3 + a[6]*T4 );
        const double hRT = -a[0]*Ti2 + a[1]*lnT*Ti + a[2]
                         + a[3]*T/2.0 + a[4]*T2/3.0 + a[5]*T3/4.0
                         + a[6]*T4/5.0 + a[7]*Ti;
        *h_out = THERMO_RU * T * hRT;
    } else {
        // 範囲外 (クランプ/線形外挿) は元の経路と完全一致させる。
        *cp_out = thermo_cp_molar(sp, T);
        *h_out  = thermo_h_molar(sp, T);
    }
}

// 混合 cp [J/(kg·K)] と h [J/kg] を 1 スイープで (質量重み)。
THERMO_HD void thermo_cph_mix(const SpeciesThermo* sp, int n, const double* Y, double T,
                              double* cp_out, double* h_out)
{
    double cp = 0.0, h = 0.0;
    for (int i=0;i<n;i++) {
        double cpi, hi;
        thermo_cph_molar(sp[i], T, &cpi, &hi);
        // 除算は thermo_cp_mass/thermo_h_mass と同一 (乗算化はビットがずれるため避ける)。
        cp += Y[i]*(cpi/sp[i].MW);
        h  += Y[i]*(hi/sp[i].MW);
    }
    *cp_out = cp;
    *h_out  = h;
}

// 混合比熱比 gamma_mix = cp/(cp-R)
THERMO_HD double thermo_gamma_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    const double cp = thermo_cp_mix(sp, n, Y, T);
    const double R  = thermo_R_mix(sp, n, Y);
    const double cv = cp - R;
    return cp / (cv > 1.0e-6 ? cv : 1.0e-6);
}

// -----------------------------------------------------------------------------
// 一般EOS陰解法ヤコビアン用の熱力学微分 (plan time_integration-general-eos-jacobian)。
//   1 セル状態 (組成 Y, 温度 T) から、固有系構築に必要な量を**同一 thermo 評価**でまとめて返す。
//   sonic[ic]/Ht[ic]/gamma[ic] を別経路で混ぜると今回直す不整合を再導入するため、ここで一括評価する。
//   κ = ∂p/∂(ρe)|_ρ = R/c_v = γ-1,  c² = γRT,  χ = c² - κ h  (CPG: h=c_pT で χ=0)。
//   h は NASA 絶対(または datum オフセット済)エンタルピー。e=h-RT, p=ρRT は呼び出し側で ρ を掛ける。
// -----------------------------------------------------------------------------
struct ThermoDerivatives {
    double R;      // 混合気体定数 [J/(kg K)]
    double cp;     // [J/(kg K)]
    double cv;     // = cp - R
    double gamma;  // cp/cv
    double h;      // 絶対エンタルピー h(T) [J/kg]
    double e;      // 内部エネルギー e(T) = h - R T [J/kg]
    double kappa;  // ∂p/∂(ρe)|_ρ = R/cv = γ-1
    double a2;     // 音速² c² = γ R T [m²/s²]
    double chi;    // ∂p/∂ρ|_{ρe} = c² - κ h  (CPG で 0)
};

THERMO_HD void thermo_derivatives_mix(const SpeciesThermo* sp, int n, const double* Y, double T,
                                      ThermoDerivatives* D)
{
    double cp, h;
    thermo_cph_mix(sp, n, Y, T, &cp, &h);   // cp,h を 1 スイープ (融合評価)
    const double R  = thermo_R_mix(sp, n, Y);
    const double cv = cp - R;
    const double cvs = (cv > 1.0e-6 ? cv : 1.0e-6);
    D->R     = R;
    D->cp    = cp;
    D->cv    = cv;
    D->gamma = cp / cvs;
    D->h     = h;
    D->e     = h - R*T;
    D->kappa = R / cvs;            // = γ-1
    D->a2    = D->gamma * R * T;   // c²
    D->chi   = D->a2 - D->kappa*h; // CPG: h=c_pT → χ=0
}

// -----------------------------------------------------------------------------
// 輸送係数 (kinetic theory, M3)。Chapman-Enskog 単成分 + Wilke/Mason-Saxena 混合。
//   - 衝突積分 Ω(2,2)*, Ω(1,1)* は Neufeld et al. (1972) の閉形式近似 (0.3<=T*<=100)。
//   - 単位は CGS 換算定数で Pa·s / W·m⁻¹·K⁻¹ / m²·s⁻¹ に整える。
//   - **内部 FP32 (M6)**: 輸送係数は桁落ちが無く (絶対エンタルピーのような大数+小増分が
//     無い) 単精度で十分。consumer GPU (FP64=FP32/64) では pow/exp/sqrt を毎セル
//     O(n_sp²) 回まわす本パスが律速になるため `expf/powf/sqrtf` で FP32 化する。
//     化学種データ (MW/σ/ε) は double のまま読み、評価のみ float。
//     viscMethod==0 (一定値) は本パスを通らないため非粘性回帰は不変。
// -----------------------------------------------------------------------------

// 換算衝突積分 Ω*(2,2) (粘性・熱伝導)。T* = T/(ε/kB)。
THERMO_HD float thermo_omega22(float Tstar)
{
    if (Tstar < 0.3f)   Tstar = 0.3f;
    if (Tstar > 100.0f) Tstar = 100.0f;
    return 1.16145f*powf(Tstar, -0.14874f)
         + 0.52487f*expf(-0.77320f*Tstar)
         + 2.16178f*expf(-2.43787f*Tstar);
}

// 換算衝突積分 Ω*(1,1) (拡散)。M4 の二元拡散係数で使用。
THERMO_HD float thermo_omega11(float Tstar)
{
    if (Tstar < 0.3f)   Tstar = 0.3f;
    if (Tstar > 100.0f) Tstar = 100.0f;
    return 1.06036f*powf(Tstar, -0.15610f)
         + 0.19300f*expf(-0.47635f*Tstar)
         + 1.03587f*expf(-1.52996f*Tstar)
         + 1.76474f*expf(-3.89411f*Tstar);
}

// 単成分 Chapman-Enskog 粘性 [Pa·s]。
//   μ = 2.6693e-6 √(M[g/mol]·T) / (σ[Å]²·Ω(2,2)*)
THERMO_HD float thermo_mu_species(const SpeciesThermo& sp, double T)
{
    const float Tf    = (float)T;
    const float epsf  = (float)(sp.eps_kB > 1.0e-30 ? sp.eps_kB : 1.0e-30);
    const float Tstar = Tf / epsf;
    const float Mg    = (float)(sp.MW * 1000.0);          // kg/mol -> g/mol
    const float sig2  = (float)(sp.sigma_LJ * sp.sigma_LJ);// Å²
    const float om    = thermo_omega22(Tstar);
    return 2.6693e-6f * sqrtf(Mg * Tf) / (sig2 * om);
}

// 単成分 熱伝導 [W/(m·K)] (modified Eucken: λ = μ(cp_mass + 1.25 R_s))。
THERMO_HD float thermo_lambda_species(const SpeciesThermo& sp, double T)
{
    const float mu = thermo_mu_species(sp, T);
    const float cp = (float)thermo_cp_mass(sp, T);
    const float R  = (float)thermo_R_species(sp);
    return mu * (cp + 1.25f*R);
}

// 質量分率 Y -> モル分率 X (X[i]=(Y_i/MW_i)/Σ_j Y_j/MW_j)。組成は double のまま (桁落ち無し・安価)。
THERMO_HD void thermo_X_from_Y(const SpeciesThermo* sp, int n, const double* Y, double* X)
{
    double s = 0.0;
    for (int i=0;i<n;i++) { X[i] = Y[i]/sp[i].MW; s += X[i]; }
    const double inv = 1.0/(s > 1.0e-300 ? s : 1.0e-300);
    for (int i=0;i<n;i++) X[i] *= inv;
}

// Wilke の φ_ij = [1+√(μ_i/μ_j)(M_j/M_i)^¼]² / √(8(1+M_i/M_j))。
THERMO_HD float thermo_wilke_phi(float mu_i, float mu_j, float M_i, float M_j)
{
    const float r   = sqrtf(mu_i/mu_j) * powf(M_j/M_i, 0.25f);
    const float num = (1.0f + r)*(1.0f + r);
    const float den = sqrtf(8.0f*(1.0f + M_i/M_j));
    return num/den;
}

// Wilke 混合粘性 [Pa·s]。X はモル分率 (double 配列を float で読む)。
THERMO_HD float thermo_mu_mix(const SpeciesThermo* sp, int n, const double* X, double T)
{
    if (n == 1) return thermo_mu_species(sp[0], T);
    float mu_s[THERMO_MAX_SPECIES];
    for (int i=0;i<n;i++) mu_s[i] = thermo_mu_species(sp[i], T);
    float mu = 0.0f;
    for (int i=0;i<n;i++) {
        const float Xi = (float)X[i];
        if (Xi <= 0.0f) continue;
        float denom = 0.0f;
        for (int j=0;j<n;j++)
            denom += (float)X[j]*thermo_wilke_phi(mu_s[i], mu_s[j], (float)sp[i].MW, (float)sp[j].MW);
        if (denom > 1.0e-30f) mu += Xi*mu_s[i]/denom;
    }
    return mu;
}

// Mason-Saxena (Wassiljewa, Wilke の φ 流用) 混合熱伝導 [W/(m·K)]。
THERMO_HD float thermo_lambda_mix(const SpeciesThermo* sp, int n, const double* X, double T)
{
    if (n == 1) return thermo_lambda_species(sp[0], T);
    float mu_s[THERMO_MAX_SPECIES];
    float la_s[THERMO_MAX_SPECIES];
    for (int i=0;i<n;i++) { mu_s[i] = thermo_mu_species(sp[i], T); la_s[i] = thermo_lambda_species(sp[i], T); }
    float la = 0.0f;
    for (int i=0;i<n;i++) {
        const float Xi = (float)X[i];
        if (Xi <= 0.0f) continue;
        float denom = 0.0f;
        for (int j=0;j<n;j++)
            denom += (float)X[j]*thermo_wilke_phi(mu_s[i], mu_s[j], (float)sp[i].MW, (float)sp[j].MW);
        if (denom > 1.0e-30f) la += Xi*la_s[i]/denom;
    }
    return la;
}

// 二元拡散係数 D_ij [m²/s] (Chapman-Enskog)。P [Pa]。M4 の混合平均拡散で使用。
//   D_ij = 1.8583e-7 √(T³(1/M_i+1/M_j)[1/(g/mol)]) / (P[atm]·σ_ij²[Å²]·Ω(1,1)*)  [cm²/s] -> m²/s
THERMO_HD float thermo_Dbinary(const SpeciesThermo& a, const SpeciesThermo& b, double T, double P)
{
    const float Tf   = (float)T;
    const float Mi   = (float)(a.MW*1000.0), Mj = (float)(b.MW*1000.0);  // g/mol
    const float sig  = (float)(0.5*(a.sigma_LJ + b.sigma_LJ));           // Å
    const float epsp = (float)(a.eps_kB * b.eps_kB);
    const float eps  = sqrtf(epsp > 1.0e-30f ? epsp : 1.0e-30f);        // K
    const float Tstar = Tf / eps;
    const float om   = thermo_omega11(Tstar);
    const float Patm = (float)(P / 101325.0);
    const float Dcm2 = 1.8583e-3f * sqrtf(Tf*Tf*Tf*(1.0f/Mi + 1.0f/Mj))
                       / (Patm * sig*sig * om);                          // cm²/s
    return Dcm2 * 1.0e-4f;                                               // -> m²/s
}

// 混合平均拡散係数 D_i [m²/s] (M4): D_i = (1-X_i)/Σ_{j≠i} X_j/D_ij。
THERMO_HD float thermo_Dmix_species(const SpeciesThermo* sp, int n, const double* X,
                                    int i, double T, double P)
{
    if (n == 1) return 0.0f;
    float denom = 0.0f;
    for (int j=0;j<n;j++) {
        if (j==i) continue;
        const float Dij = thermo_Dbinary(sp[i], sp[j], T, P);
        denom += (float)X[j]/(Dij > 1.0e-30f ? Dij : 1.0e-30f);
    }
    if (denom < 1.0e-30f) return thermo_Dbinary(sp[i], sp[i], T, P);
    return (1.0f - (float)X[i])/denom;
}

// -----------------------------------------------------------------------------
// 温度反転: 比内部エネルギー e [J/kg] と組成 Y から T を Newton 反転で求める。
//   f(T) = e_mix(T) - e = 0, f'(T) = cv_mix(T) = cp_mix(T) - R_mix (厳密微分)
//   初期値 T_guess (前ステップ T) でウォームスタート, 毎反復クランプ。
// -----------------------------------------------------------------------------
THERMO_HD double thermo_T_from_e(const SpeciesThermo* sp, int n, const double* Y,
                                 double e, double T_guess,
                                 double T_min, double T_max)
{
    const double R = thermo_R_mix(sp, n, Y);
    double T = T_guess;
    if (!(T > T_min)) T = T_min;      // NaN/未初期化ガード
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it=0; it<20; ++it) {
        double h_T, cp_T;
        thermo_cph_mix(sp, n, Y, T, &cp_T, &h_T);             // cp,h を 1 スイープ
        const double e_T = h_T - R*T;                         // e_mix(T)
        const double cv  = cp_T - R;                          // de/dT
        const double cvf = (cv > 1.0e-2*R ? cv : 1.0e-2*R);   // 0 割回避フロア
        double dT = (e_T - e)/cvf;
        // 1 反復あたりのステップを制限 (発散防止)
        if (dT >  0.5*T) dT =  0.5*T;
        if (dT < -0.5*T) dT = -0.5*T;
        T -= dT;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*T) break;
    }
    return T;
}

// 温度反転 (エンタルピー版): 比エンタルピー h [J/kg] と組成 Y から T を求める。
//   f(T) = h_mix(T) - h = 0, f'(T) = cp_mix(T) (厳密微分)。Roe 平均状態の有効 γ 用。
THERMO_HD double thermo_T_from_h(const SpeciesThermo* sp, int n, const double* Y,
                                 double h, double T_guess,
                                 double T_min, double T_max)
{
    double T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it=0; it<20; ++it) {
        double h_T, cp;
        thermo_cph_mix(sp, n, Y, T, &cp, &h_T);              // cp,h を 1 スイープ
        const double cpf = (cp > 1.0e-3 ? cp : 1.0e-3);
        double dT = (h_T - h)/cpf;
        if (dT >  0.5*T) dT =  0.5*T;
        if (dT < -0.5*T) dT = -0.5*T;
        T -= dT;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*T) break;
    }
    return T;
}

// 標準状態エントロピー S^o/R (圧力項を含まない). Roe/等エントロピー BC 用。
THERMO_HD double thermo_s_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const double* a = thermo_pick_coeffs(sp, Tc);
    const double Ti = 1.0/Tc; const double Ti2 = Ti*Ti; const double lnT = log(Tc);
    const double sR = -a[0]*Ti2/2.0 - a[1]*Ti + a[2]*lnT + a[3]*Tc
                    + a[4]*Tc*Tc/2.0 + a[5]*Tc*Tc*Tc/3.0 + a[6]*Tc*Tc*Tc*Tc/4.0 + a[8];
    return THERMO_RU * sR;
}
THERMO_HD double thermo_s0_mass(const SpeciesThermo& sp, double T)
{
    // 域外は cp 一定外挿 (h の線形外挿と整合): ds0/dT = cp/T → s0(T)=s0(Tb)+cp(Tb)·ln(T/Tb)。
    // 旧実装は s0 をクランプ (=const) しており h/cp の定 cp 外挿と熱力学的に非整合で、
    // sub-200K (低温ノズル) の等エントロピー関係を破綻させていた (TP 不安定の一因)。
    if (T < sp.Tlo) {
        const double s0  = thermo_s_molar_clamped(sp, sp.Tlo);
        const double cp0 = thermo_cp_molar_clamped(sp, sp.Tlo);
        return (s0 + cp0*log(T/sp.Tlo)) / sp.MW;
    }
    if (T > sp.Thi) {
        const double s1  = thermo_s_molar_clamped(sp, sp.Thi);
        const double cp1 = thermo_cp_molar_clamped(sp, sp.Thi);
        return (s1 + cp1*log(T/sp.Thi)) / sp.MW;
    }
    return thermo_s_molar_clamped(sp, T) / sp.MW;   // [J/(kg·K)]
}

// 単成分 (M1) の等エントロピー total→static 反転。
//   全温 Tt・全圧 Pt と局所マッハ M から静的 (Ts, Ps, ρ, |u|) を求める。
//     h(Tt) = h(Ts) + ½u²,  u = M·a(Ts),  a=√(γ(Ts)R Ts)
//     Ps = Pt·exp(-(s0(Tt)-s0(Ts))/R)
THERMO_HD void thermo_isentropic_from_total_single(
    const SpeciesThermo* sp, double Pt, double Tt, double mach,
    double* Ts_out, double* Ps_out, double* ro_out, double* umag_out)
{
    const double R  = thermo_R_species(sp[0]);
    const double h0 = thermo_h_mass(sp[0], Tt);
    double Ts = Tt/(1.0 + 0.2*mach*mach);   // CPG 近似の初期値
    if (Ts < 10.0) Ts = 10.0;
    #pragma unroll 1
    for (int it=0; it<25; ++it) {
        const double cp = thermo_cp_mass(sp[0], Ts);
        const double g  = cp/((cp-R) > 1.0e-6 ? (cp-R) : 1.0e-6);
        const double a2 = g*R*Ts;
        const double f  = thermo_h_mass(sp[0], Ts) + 0.5*mach*mach*a2 - h0;
        const double df = cp + 0.5*mach*mach*g*R;   // 近似微分 (dγ/dT 無視)
        double dT = f/df;
        if (dT >  0.4*Ts) dT =  0.4*Ts;
        if (dT < -0.4*Ts) dT = -0.4*Ts;
        Ts -= dT;
        if (Ts < 10.0) Ts = 10.0;
        if (Ts > Tt)   Ts = Tt;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*Ts) break;
    }
    const double cp = thermo_cp_mass(sp[0], Ts);
    const double g  = cp/((cp-R) > 1.0e-6 ? (cp-R) : 1.0e-6);
    const double a  = sqrt(g*R*Ts);
    *umag_out = mach*a;
    *Ps_out   = Pt*exp(-(thermo_s0_mass(sp[0],Tt) - thermo_s0_mass(sp[0],Ts))/R);
    *Ts_out   = Ts;
    *ro_out   = (*Ps_out)/(R*Ts);
}

// 混合 (質量分率 Y, 組成一定) の標準状態エントロピー s0_mix(T) = Σ Y_i s0_i(T) [J/(kg·K)]。
// 等エントロピー反転では組成が一定なので混合エントロピー項は差分でキャンセルし不要。
THERMO_HD double thermo_s0_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += Y[i]*thermo_s0_mass(sp[i], T);
    return s;
}

// 多成分版の等エントロピー total→static 反転 (Mach 指定版)。組成 Y は一定。
// thermo_isentropic_from_total_single の sp[0] を混合則 (R_mix/h_mix/cp_mix/s0_mix) に置換しただけ。
THERMO_HD void thermo_isentropic_from_total_mix(
    const SpeciesThermo* sp, int n, const double* Y, double Pt, double Tt, double mach,
    double* Ts_out, double* Ps_out, double* ro_out, double* umag_out)
{
    const double R  = thermo_R_mix(sp, n, Y);
    const double h0 = thermo_h_mix(sp, n, Y, Tt);
    double Ts = Tt/(1.0 + 0.2*mach*mach);   // CPG 近似の初期値
    if (Ts < 10.0) Ts = 10.0;
    #pragma unroll 1
    for (int it=0; it<25; ++it) {
        const double cp = thermo_cp_mix(sp, n, Y, Ts);
        const double g  = cp/((cp-R) > 1.0e-6 ? (cp-R) : 1.0e-6);
        const double a2 = g*R*Ts;
        const double f  = thermo_h_mix(sp, n, Y, Ts) + 0.5*mach*mach*a2 - h0;
        const double df = cp + 0.5*mach*mach*g*R;   // 近似微分 (dγ/dT 無視)
        double dT = f/df;
        if (dT >  0.4*Ts) dT =  0.4*Ts;
        if (dT < -0.4*Ts) dT = -0.4*Ts;
        Ts -= dT;
        if (Ts < 10.0) Ts = 10.0;
        if (Ts > Tt)   Ts = Tt;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*Ts) break;
    }
    const double cp = thermo_cp_mix(sp, n, Y, Ts);
    const double g  = cp/((cp-R) > 1.0e-6 ? (cp-R) : 1.0e-6);
    const double a  = sqrt(g*R*Ts);
    *umag_out = mach*a;
    *Ps_out   = Pt*exp(-(thermo_s0_mix(sp,n,Y,Tt) - thermo_s0_mix(sp,n,Y,Ts))/R);
    *Ts_out   = Ts;
    *ro_out   = (*Ps_out)/(R*Ts);
}

// 単成分の等エントロピー total→static 反転 (静圧 Ps 指定版)。
//   全温 Tt・全圧 Pt と静圧 Ps から静温 Ts・密度 ρ・速度 |u| を求める。
//     s0(Ts) = s0(Tt) + R·ln(Ps/Pt)   (等エントロピー: s0(T)-R ln P 一定)
//     |u| = √(2 max(h(Tt)-h(Ts), 0)),  ρ = Ps/(R Ts)
THERMO_HD void thermo_isentropic_from_total_Ps_single(
    const SpeciesThermo* sp, double Pt, double Tt, double Ps,
    double* Ts_out, double* ro_out, double* umag_out)
{
    const double R  = thermo_R_species(sp[0]);
    const double s_target = thermo_s0_mass(sp[0], Tt) + R*log((Ps > 1.0e-30 ? Ps : 1.0e-30)
                                                              /(Pt > 1.0e-30 ? Pt : 1.0e-30));
    double Ts = Tt;   // Ps<=Pt なら Ts<=Tt。Tt から下降側へ Newton。
    if (Ts < 10.0) Ts = 10.0;
    #pragma unroll 1
    for (int it=0; it<30; ++it) {
        const double f  = thermo_s0_mass(sp[0], Ts) - s_target;   // s0(Ts)-target
        const double df = thermo_cp_mass(sp[0], Ts)/Ts;            // ds0/dT = cp/T
        double dT = f/(df > 1.0e-12 ? df : 1.0e-12);
        if (dT >  0.4*Ts) dT =  0.4*Ts;
        if (dT < -0.4*Ts) dT = -0.4*Ts;
        Ts -= dT;
        if (Ts < 10.0) Ts = 10.0;
        if (Ts > Tt)   Ts = Tt;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*Ts) break;
    }
    const double dh = thermo_h_mass(sp[0], Tt) - thermo_h_mass(sp[0], Ts);
    *umag_out = sqrt(dh > 0.0 ? 2.0*dh : 0.0);
    *Ts_out   = Ts;
    *ro_out   = Ps/(R*Ts);
}

// thermo_isentropic_from_total_Ps_single の sp[0] を混合則 (R_mix/s0_mix/cp_mix/h_mix) に置換。
//   全温 Tt・全圧 Pt と静圧 Ps から静温 Ts・密度 ρ・速度 |u| を求める (多成分 TP)。
THERMO_HD void thermo_isentropic_from_total_Ps_mix(
    const SpeciesThermo* sp, int n, const double* Y, double Pt, double Tt, double Ps,
    double* Ts_out, double* ro_out, double* umag_out)
{
    const double R  = thermo_R_mix(sp, n, Y);
    const double s_target = thermo_s0_mix(sp, n, Y, Tt) + R*log((Ps > 1.0e-30 ? Ps : 1.0e-30)
                                                               /(Pt > 1.0e-30 ? Pt : 1.0e-30));
    double Ts = Tt;
    if (Ts < 10.0) Ts = 10.0;
    #pragma unroll 1
    for (int it=0; it<30; ++it) {
        const double f  = thermo_s0_mix(sp, n, Y, Ts) - s_target;
        const double df = thermo_cp_mix(sp, n, Y, Ts)/Ts;
        double dT = f/(df > 1.0e-12 ? df : 1.0e-12);
        if (dT >  0.4*Ts) dT =  0.4*Ts;
        if (dT < -0.4*Ts) dT = -0.4*Ts;
        Ts -= dT;
        if (Ts < 10.0) Ts = 10.0;
        if (Ts > Tt)   Ts = Tt;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*Ts) break;
    }
    const double dh = thermo_h_mix(sp, n, Y, Tt) - thermo_h_mix(sp, n, Y, Ts);
    *umag_out = sqrt(dh > 0.0 ? 2.0*dh : 0.0);
    *Ts_out   = Ts;
    *ro_out   = Ps/(R*Ts);
}

// =============================================================================
// device 側化学種データへのアクセス (thermo_d.cu が所有)
//   - thermo_init_db: solverConfig から DB を構築し device へアップロード
//   - thermo_species_device_ptr / thermo_num_species: カーネル引数用に取得
//   - thermo_species_host: host 計算 (初期条件) 用の CPU 側配列
// =============================================================================
class solverConfig;

void                 thermo_init_db(solverConfig& cfg);
const SpeciesThermo* thermo_species_device_ptr();   // device global memory
int                  thermo_num_species();
const SpeciesThermo* thermo_species_host();          // host array (length = thermo_num_species)

// =============================================================================
// 共通 helper: セル組成の取り出しと「所与温度 T での気体状態」評価
//   BC ゴースト構成 (wall_isothermal 等)・node 壁ノード温度ピンなど、
//   「セル ic の組成 (TP) / CPG 定数で、指定温度 T の R, e(T), cp(T), γ(T) が要る」
//   頻出パターンの重複排除 (2026-07-20)。boundaryCond_d.cu の bc_cell_Y はこの
//   thermo_cell_Y への転送で維持される。
// =============================================================================

// セル ic の質量分率 Y[] を roY から取り出して正規化する。多成分 (n>=2) で true。
// 単成分/roY 無しは false (Y は未定義のままなので、呼び手は sp[0] 単成分式を使う)。
#if defined(__CUDACC__)
__device__ inline bool thermo_cell_Y(flow_float* const* roY, int n, geom_int ic, double* Y)
{
    if (roY == nullptr || n < 2) return false;
    double s = 0.0;
    for (int k = 0; k < n; ++k) { Y[k] = (double)roY[k][ic]; s += Y[k]; }
    const double inv = 1.0 / (s > 1.0e-300 ? s : 1.0e-300);
    for (int k = 0; k < n; ++k) Y[k] *= inv;
    return true;
}

// 温度 T における気体状態係数 (R は温度非依存・組成依存)。
struct GasStateAtT {
    flow_float R;      // 気体定数 [J/(kg K)]
    flow_float e;      // 内部エネルギー e(T) [J/kg] (運動エネルギー含まず)
    flow_float cp;     // 定圧比熱 cp(T) [J/(kg K)]
    flow_float gamma;  // 比熱比 γ(T)
};

// セル ic の組成 (thermalMethod==2, NASA) または CPG 定数 (それ以外) で、温度 T の
// R / e / cp / γ を返す。TP 側は内部 double で評価 (既存 BC カーネルと同じ演算列)。
__device__ inline GasStateAtT thermo_state_at_T(
    int thermalMethod, flow_float ga, flow_float cp_const,
    const SpeciesThermo* sp, flow_float* const* roY, int nSpecies,
    geom_int ic, flow_float T)
{
    GasStateAtT g;
    if (thermalMethod == 2) {
        double Y[THERMO_MAX_SPECIES];
        const bool mix = thermo_cell_Y(roY, nSpecies, ic, Y);
        const int n = (nSpecies >= 1) ? nSpecies : 1;
        const double Td  = (double)T;
        const double Rd  = mix ? thermo_R_mix(sp, n, Y) : thermo_R_species(sp[0]);
        const double ed  = mix ? thermo_e_mix(sp, n, Y, Td)
                               : (thermo_h_mass(sp[0], Td) - Rd*Td);
        const double cpd = mix ? thermo_cp_mix(sp, n, Y, Td) : thermo_cp_mass(sp[0], Td);
        g.R     = (flow_float)Rd;
        g.e     = (flow_float)ed;
        g.cp    = (flow_float)cpd;
        g.gamma = (flow_float)(cpd/((cpd - Rd) > 1.0e-6 ? (cpd - Rd) : 1.0e-6));
    } else {
        g.R     = cp_const - cp_const/ga;
        g.e     = (cp_const/ga) * T;   // e = cv T, cv = cp/γ
        g.cp    = cp_const;
        g.gamma = ga;
    }
    return g;
}
#endif // __CUDACC__
