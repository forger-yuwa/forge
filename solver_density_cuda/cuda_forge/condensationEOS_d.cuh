#pragma once

#include "thermo_d.cuh"
#include "condensationProperties_d.cuh"

// 一温度 二相 温度反転 (carrier + condensible, thermally perfect; Phase 3 H2O)。
// 気相 = 全化学種 (蒸気とみなした NASA-9 混合)。凝縮種 (H2O) の液相分率 g。
//   混合内部エネルギー e_mix(T) = e_gas^全蒸気(Y,T) + g(R_w T - L_w(T))
//     (液 e_l = e_v + R_w T - L_w、凝縮ぶん g の補正)。
// e_mix(T)=e_in を Newton で解く。Y は総質量分率 (N2 + 総水)。g≤Y_H2O は呼び出し側で保証。
//   Rw = R_u/M_H2O (凝縮種の比気体定数)、cprops = 凝縮種物性 (model=H2O)。
__host__ __device__ inline double cond_T_from_e_carrier(
    const SpeciesThermo* sp, int nSp, const double* Y,
    double e_in, double g, double Rw, const CondSpeciesProps& cprops,
    double T_guess, double T_min, double T_max)
{
    const double Rmix = thermo_R_mix(sp, nSp, Y);
    double T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it = 0; it < 30; ++it) {
        double cp_T, h_T;
        thermo_cph_mix(sp, nSp, Y, T, &cp_T, &h_T);          // 全蒸気混合 cp,h
        const double e_allvap = h_T - Rmix*T;                 // 全蒸気の内部エネルギー
        const double cv = cp_T - Rmix;
        const double L  = cond_latent(cprops, T);
        const double dL = (cond_latent(cprops, T+0.1) - cond_latent(cprops, T-0.1)) / 0.2;
        const double e_mix = e_allvap + g*(Rw*T - L);
        const double demix = cv + g*(Rw - dL);
        const double demf  = (demix > 1.0e-2*cv) ? demix : 1.0e-2*cv;
        double dT = (e_mix - e_in)/demf;
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

// 非平衡凝縮: 一温度 二相 EOS の温度反転 (Phase 2)。methods/condensation.md 5 節。
//
// 【一温度近似】気相温度 T_v と液滴温度 T_d を分けず T_v=T_d=T とする (初期実装)。
// 将来は e=(1-g)e_v(T_v)+g e_l(T_d) と Hill T_d 式の 2 変数局所 Newton へ拡張する設計
// (本関数は T 1 変数 Newton; 拡張時は (T_v,T_d) の 2x2 Newton に置換)。
//
// 【混合内部エネルギー】論文形 e = (1-g)e_v(T) + g e_l(T)、液相 e_l(T)=e_v(T)-L(T) より
//   e = e_v(T) - g L(T)
// となる。気相 e_v(T) は forge の thermally-perfect thermo (NASA-9, thermo_cph_mix) で評価し
// **定比熱近似は使わない**。g は液相質量分率 (総和)。L は凝縮種の潜熱 (condensationProperties)。
//
// 【温度反転】保存量から e_in = roe/ρ - ½|u|² を作り、
//   G(T) = e_v(T) - g L(T) - e_in = 0,   G'(T) = c_v(T) - g L'(T)
// を Newton で解く。g=0 のとき G=e_v(T)-e_in で thermo_T_from_e と一致 → 単相 TP に厳密縮約
// (呼び出し側は g≈0 で従来 thermo_T_from_e を使い厳密同一を保証する)。
//
// 【圧力】実在気体補正なし。液滴は圧力を持たないとして p = (1-g)ρ R_v T。
//
// 現状 N2 のみ (潜熱 n2_latent)。多成分凝縮では Σ_sp g_sp L_sp(T) に一般化する (将来)。

// 一温度 二相 温度反転 (calorically perfect gas; 初期実装)。論文 Eq.16-18 (g0=0):
//   液相 e_l = e_v + R_v T - L(T) = c_p T - L  (∵ h_v-h_l=L, h_v=e_v+R_vT, h_l≈e_l)
//   e = (1-g)e_v + g e_l = (c_v + g R_v) T - g L(T)   (e_v=c_v T)
//   ⇒ T = (e_in + g L(T))/(c_v + g R_v)  (= Eq.18, Cv0=Cvv=c_v)
// **L の温度依存は入れる** (n2_latent(T))。g=0 で T=e_in/cv (従来 CPG と一致)。cv=cp/γ, R=(γ-1)cv。
__host__ __device__ inline double cond_T_from_e_cpg(
    double e_in, double g_tot, double cv, double R, double T_guess,
    const CondSpeciesProps& cprops)
{
    const double a = cv + g_tot*R;   // 実効熱容量 (T に対し一定)
    double T = T_guess;
    if (!(T > 1.0)) T = 1.0;
    #pragma unroll 1
    for (int it = 0; it < 30; ++it) {
        const double L  = cond_latent(cprops, T);
        const double dL = (cond_latent(cprops, T + 0.1) - cond_latent(cprops, T - 0.1)) / 0.2;
        const double G  = a*T - g_tot*L - e_in;
        const double Gp = a - g_tot*dL;
        const double Gpf = (Gp > 1.0e-2*a) ? Gp : 1.0e-2*a;
        double dT = G / Gpf;
        if (dT >  0.5*T) dT =  0.5*T;
        if (dT < -0.5*T) dT = -0.5*T;
        T -= dT;
        if (T < 1.0)    T = 1.0;
        if (T > 6000.0) T = 6000.0;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*T) break;
    }
    return T;
}

// 一温度 二相 温度反転 (thermally perfect; 後続拡張)。e_in [J/kg], g_tot=総液相質量分率, T_guess。
__host__ __device__ inline double cond_T_from_e_onetemp(
    const SpeciesThermo* sp, int nSp, const double* Y,
    double e_in, double g_tot,
    double T_guess, double T_min, double T_max)
{
    const double R = thermo_R_mix(sp, nSp, Y);
    double T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it = 0; it < 30; ++it) {
        double cp_T, h_T;
        thermo_cph_mix(sp, nSp, Y, T, &cp_T, &h_T);   // 気相 cp,h (TP)
        const double e_v = h_T - R*T;                  // 気相内部エネルギー
        const double cv  = cp_T - R;                   // de_v/dT
        const double L   = n2_latent(T);               // 潜熱 (N2)
        const double dL  = (n2_latent(T + 0.1) - n2_latent(T - 0.1)) / 0.2; // dL/dT (数値)
        // e = (1-g)e_v + g e_l, e_l = e_v + R T - L ⇒ e = e_v + g R T - g L
        const double G   = e_v + g_tot*R*T - g_tot*L - e_in;
        const double Gp  = cv + g_tot*R - g_tot*dL;
        const double Gpf = (Gp > 1.0e-2*R) ? Gp : 1.0e-2*R; // 0 割回避フロア
        double dT = G / Gpf;
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

// 二相 frozen 音速 (固定 g,Y; plans/active/condensation-kantrowitz-gamma-twophase-sonic.md §4.2, methods/condensation.md §5)。
//   一温度二相 EOS p=ρ R_eff T, e=e^全蒸気(T)+g(R_w T−L(T)) の g,Y 固定の等エントロピー微分:
//     c² = (∂p/∂ρ)_e + (p/ρ²)(∂p/∂e)_ρ = γ_2φ R_eff T,
//     c_p,2φ = c_p^全蒸気 − g L'(T)  (L'=c_p,v−c_l なので「g ぶんの蒸気 c_p を液 c_l に置換」),  c_v,2φ = c_p,2φ − R_eff.
//   carrier: R_eff=R_mix−g R_w、pure: R_eff=(1−g)R。c_v,2φ は温度反転 Newton の de/dT と同一。
//   戻り値 false = 熱力学的に不整合 (c_v,2φ<=0.05 c_p,2φ、N2 潜熱フィットの低温 L'>0 など) → 呼び出し側は旧式にフォールバック。
__host__ __device__ inline bool cond_twophase_sonic(
    double cp_allvap, double R_eff, double g, double dLdT, double T, double* gamma2, double* c2)
{
    const double cp2 = cp_allvap - g*dLdT;
    const double cv2 = cp2 - R_eff;
    if (!(cv2 > 0.05*cp2) || !(R_eff > 0.0) || !(T > 0.0)) return false;
    *gamma2 = cp2/cv2;
    *c2 = (*gamma2)*R_eff*T;
    return (*c2 > 0.0);
}

// 平衡凝縮 (condEquilibrium=1): 現セル状態から e 一定で液相を Δ だけ変えたときの
//   T(Δ) = T + Δ (L - R_w T)/(c_vg + g R_w)  (二相 EOS の線形化), p_v(Δ) = ρ (y_v0 - Δ) R_w T(Δ)
// が p_sat(T(Δ)) と釣り合う Δ ∈ [-g, y_v0] を二分法で解く (F は Δ について単調減少)。
// 戻り値 Δ (>0 凝縮, <0 蒸発)。plans/accepted/condensation-equilibrium.md。
__host__ __device__ inline double cond_equilibrium_delta(
    const CondSpeciesProps& cprops, double rod, double Td, double g, double yv0, double Rw, double cvg)
{
    const double L = cond_latent(cprops, Td);
    const double slope = (L - Rw*Td)/(cvg + g*Rw);          // dT/dΔ
    auto F = [&](double D) {
        double Tn = Td + D*slope; if (Tn < 50.0) Tn = 50.0;
        double yv = yv0 - D; if (yv < 0.0) yv = 0.0;
        return rod*yv*Rw*Tn - cond_psat(cprops, Tn);
    };
    double lo = -g, hi = yv0;
    if (hi <= 0.0) return (g > 0.0 ? -g : 0.0);
    const double F0 = F(0.0);
    if (F0 <= 0.0 && g <= 0.0) return 0.0;                 // 未飽和・液なし
    if (F(lo) < 0.0) return lo;                             // 全蒸発しても未飽和
    if (F(hi) > 0.0) return hi;                             // (理論上起きない) 全凝縮
    if (F0 > 0.0) lo = 0.0; else hi = 0.0;
    #pragma unroll 1
    for (int it = 0; it < 48; ++it) {
        const double mid = 0.5*(lo + hi);
        if (F(mid) > 0.0) lo = mid; else hi = mid;
        if (hi - lo < 1.0e-9*(fabs(yv0) + 1.0e-12)) break;
    }
    return 0.5*(lo + hi);
}

// ---------------------------------------------------------------------------------------------
// 平衡凝縮・EOS 拘束形 (condEquilibrium=2, plans/accepted/condensation-equilibrium-eos.md)
//   液相分率 g を輸送量ではなく状態量として、保存量 (ρ, e_in, Y_w) から (T, g) を同時反転する。
//     未飽和: g=0, e_gas(T)=e_in, p_v=ρ Y_w R_w T <= p_sat(T)
//     飽和  : e_gas(T)+g(R_w T-L)=e_in かつ ρ (Y_w-g) R_w T = p_sat(T),  0<g<=Y_w
//   G(g) = ρ(Y_w-g)R_w T(g) - p_sat(T(g)) は g について単調減少 (g↑: 蒸気↓, 潜熱で T↑→p_sat↑) で
//   G(0)>0, G(Y_w)=-p_sat<0 なので解は唯一。括弧付き Newton (解析 dG/dg) + 二分法退避。
//   T(g) の内側反転は呼び出し側が渡す関数 (carrier: cond_T_from_e_carrier / pure: onetemp, cpg)。
// ---------------------------------------------------------------------------------------------

// 共通反復器。tfun(g, T_guess) が T(g) を返す。pv_coef*(Yw-g)*T = 蒸気分圧 (pv_coef=ρR_w、pure は Yw=1)。
// dTdg(g,T) は dT/dg (Newton 用、粗くてよい: 括弧で保護)。戻り値 g、*T_out に T。
template <class TFun, class DTFun>
__host__ __device__ inline double cond_eq_solve_g(
    TFun tfun, DTFun dTdg, const CondSpeciesProps& cprops,
    double pv_coef, double Yw, double g_guess, double T_guess, double* T_out)
{
    // 1) 乾き判定
    double T0 = tfun(0.0, T_guess);
    const double G0 = pv_coef*Yw*T0 - cond_psat(cprops, T0);
    if (!(Yw > 0.0) || G0 <= 0.0) { *T_out = T0; return 0.0; }
    // 2) 括弧 [lo,hi] = [0, Yw] (G(lo)>0, G(hi)<0)。初期 g は前ステップ値をクランプ。
    double lo = 0.0, hi = Yw;
    double g = g_guess;
    if (!(g > 0.0)) g = 0.5*Yw;
    if (g >= Yw) g = 0.999*Yw;
    double T = T0;
    #pragma unroll 1
    for (int it = 0; it < 40; ++it) {
        T = tfun(g, T);
        const double ps = cond_psat(cprops, T);
        const double G  = pv_coef*(Yw - g)*T - ps;
        if (fabs(G) <= 1.0e-8*ps) break;
        if (G > 0.0) lo = g; else hi = g;
        if (hi - lo <= 1.0e-10*Yw) { g = 0.5*(lo + hi); break; }
        // Newton: dG/dg = -pv_coef*T + [pv_coef*(Yw-g) - dpsat/dT] * dT/dg
        const double dps = (cond_psat(cprops, T + 0.05) - cond_psat(cprops, T - 0.05))/0.1;
        const double dGdg = -pv_coef*T + (pv_coef*(Yw - g) - dps)*dTdg(g, T);
        double gn = (dGdg < 0.0) ? g - G/dGdg : -1.0;
        if (!(gn > lo && gn < hi)) gn = 0.5*(lo + hi);   // 括弧外 → 二分法
        g = gn;
    }
    T = tfun(g, T);
    *T_out = T;
    return g;
}

// carrier + condensible (TP): Y は総質量分率 (蒸気+液の総水を含む)、Yw = Y[condGasSpecies]。
__host__ __device__ inline double cond_equilibrium_Tg_carrier(
    const SpeciesThermo* sp, int nSp, const double* Y,
    double e_in, double rho, double Yw, double Rw, const CondSpeciesProps& cprops,
    double T_guess, double g_guess, double T_min, double T_max, double* T_out)
{
    const double Rmix = thermo_R_mix(sp, nSp, Y);
    auto tfun = [&](double g, double Tg) -> double {
        if (g > 1.0e-14) return cond_T_from_e_carrier(sp, nSp, Y, e_in, g, Rw, cprops, Tg, T_min, T_max);
        return thermo_T_from_e(sp, nSp, Y, e_in, Tg, T_min, T_max);
    };
    auto dTdg = [&](double g, double T) -> double {
        double cp_T, h_T; thermo_cph_mix(sp, nSp, Y, T, &cp_T, &h_T);
        const double cv = cp_T - Rmix;
        const double L  = cond_latent(cprops, T);
        const double dL = (cond_latent(cprops, T + 0.1) - cond_latent(cprops, T - 0.1))/0.2;
        double den = cv + g*(Rw - dL); if (den < 1.0e-2*cv) den = 1.0e-2*cv;
        return (L - Rw*T)/den;
    };
    return cond_eq_solve_g(tfun, dTdg, cprops, rho*Rw, Yw, g_guess, T_guess, T_out);
}

// pure-condensible (TP; 気相=凝縮種、p_v=ρ(1-g)R T)。
__host__ __device__ inline double cond_equilibrium_Tg_pure_tp(
    const SpeciesThermo* sp, int nSp, const double* Y,
    double e_in, double rho, const CondSpeciesProps& cprops,
    double T_guess, double g_guess, double T_min, double T_max, double* T_out)
{
    const double R = thermo_R_mix(sp, nSp, Y);
    auto tfun = [&](double g, double Tg) -> double {
        if (g > 1.0e-14) return cond_T_from_e_onetemp(sp, nSp, Y, e_in, g, Tg, T_min, T_max);
        return thermo_T_from_e(sp, nSp, Y, e_in, Tg, T_min, T_max);
    };
    auto dTdg = [&](double g, double T) -> double {
        double cp_T, h_T; thermo_cph_mix(sp, nSp, Y, T, &cp_T, &h_T);
        const double cv = cp_T - R;
        const double L  = cond_latent(cprops, T);
        const double dL = (cond_latent(cprops, T + 0.1) - cond_latent(cprops, T - 0.1))/0.2;
        double den = cv + g*(R - dL); if (den < 1.0e-2*cv) den = 1.0e-2*cv;
        return (L - R*T)/den;
    };
    return cond_eq_solve_g(tfun, dTdg, cprops, rho*R, 1.0, g_guess, T_guess, T_out);
}

// pure-condensible (CPG)。cv, R は定数。
__host__ __device__ inline double cond_equilibrium_Tg_pure_cpg(
    double e_in, double rho, double cv, double R, const CondSpeciesProps& cprops,
    double T_guess, double g_guess, double* T_out)
{
    auto tfun = [&](double g, double Tg) -> double {
        if (g > 1.0e-14) return cond_T_from_e_cpg(e_in, g, cv, R, Tg, cprops);
        double T = e_in/cv; return (T > 1.0) ? T : 1.0;
    };
    auto dTdg = [&](double g, double T) -> double {
        const double L  = cond_latent(cprops, T);
        const double dL = (cond_latent(cprops, T + 0.1) - cond_latent(cprops, T - 0.1))/0.2;
        double den = cv + g*(R - dL); if (den < 1.0e-2*cv) den = 1.0e-2*cv;
        return (L - R*T)/den;
    };
    return cond_eq_solve_g(tfun, dTdg, cprops, rho*R, 1.0, g_guess, T_guess, T_out);
}
