#pragma once
// =============================================================================
// twoPhaseFaceDiag_d.cuh — 診断 D1 (G1: 0 step の面作用素 A/B) の面ごとの中間量 (診断専用; 本番経路からは呼ばない)
//   plans/active/condensation-two-phase-default.md §5.1 #4、仕様 methods/condensation.md §7c、
//   設計の判断 notes/reviews/2026-10-04-twophase-g1-g3-diag-design-diagnose.md。
//
// 本番の面代数 tp_face_flux (twoPhaseDiffusion_d.cuh) は**変更しない**。ここでは同じ入力 TpFaceInT<R> から、tp_face_flux が
// 外へ出さない中間量 (気相組成 z、分子流束 j⁰、補正の風上 z_up) を同じ演算順で組み直し、
//   分子蒸気流束 j_v = j_v⁰ − z_up,v·Σj⁰ を補正前 (j_v⁰) と補正項 (z_up,v·Σj⁰) に分けて返す。
// 出力の J・Jv・Jl・JQ・q・Sm・up0 は呼び出し側が本番の tp_face_flux<R> を呼んで得る (本関数は Sm・up0 の一致を自己検査に使う)。
//
// 誤差尺度 (G1 (ii) の事前登録; plan §5.1 #4。結果を見て広げない):
//   A_l = |ct·geo|·(|ρg₀/ρ₀| + |ρg₁/ρ₁|)
//   A_v = |ρ_g,f D_w geo|·(|z_v0| + |z_v1|) + |z_up,v|·Σ_k |ρ_g,f D_k geo|·(|z_k0| + |z_k1|)
// 尺度は double 参照の値 (R = double で呼んだ結果) から作る。
// =============================================================================
#include "cuda_forge/twoPhaseDiffusion_d.cuh"   // TpFaceInT / TpFaceOutT / tp_max / TP_HD

// 面ごとの中間量 (R = float は本番と同じ float 演算、R = double は double 参照)
template <typename R> struct TpFaceDiagT {
    R rgf;                       // 面の気相密度 ρ_g,f
    R Sm;                        // Σ j⁰ (セル 0 へ入る向き; tp_face_flux の o.Sm と同じ演算)
    int up0;                     // 補正の z をセル 0 から取ったか (tp_face_flux の o.up0 と同じ判定)
    R jv0;                       // 補正前の分子蒸気流束 j_v⁰ = ρ_g,f D_w (z_v1 − z_v0) geo
    R jv_corr;                   // 補正項 z_up,v·Σj⁰
    R jv_mol;                    // 分子蒸気流束 j_v = j_v⁰ − z_up,v·Σj⁰ (tp_face_flux の jc と同じ式)
    R jv_turb;                   // 蒸気の乱流分 ct·geo·(v₁ − v₀) (v = (ρY_w − ρg)/ρ)
    R zv0, zv1, zupv;            // z_v (セル 0/1) と補正に使った z_up,v
};

// tp_face_flux の前半 (気相組成 z・分子流束 j⁰・Σj⁰・風上) を同じ演算順で組み直す。
// 演算は tp_face_flux (twoPhaseDiffusion_d.cuh) の該当行の写し。本番側を変えたら本関数も合わせること。
template <typename R>
TP_HD inline void tp_face_flux_diag(const TpFaceInT<R>& in, TpFaceDiagT<R>& d)
{
    const int n = in.n, iw = in.iw;
    const R rgas0 = tp_max<R>(in.rho0 - in.rg0, R(1.0e-30)), rgas1 = tp_max<R>(in.rho1 - in.rg1, R(1.0e-30));
    const R ig0 = R(1.0)/rgas0, ig1 = R(1.0)/rgas1;
    R z0[THERMO_MAX_SPECIES], z1[THERMO_MAX_SPECIES];
    R s0 = R(0), s1 = R(0);
    for (int k = 0; k < n; ++k) {
        const R a0 = (k == iw) ? (in.rY0[k] - in.rg0) : in.rY0[k];
        const R a1 = (k == iw) ? (in.rY1[k] - in.rg1) : in.rY1[k];
        z0[k] = a0*ig0; z1[k] = a1*ig1;
        s0 += z0[k]; s1 += z1[k];
    }
    const R n0 = R(1.0)/((s0 > R(1.0e-30)) ? s0 : R(1.0e-30)), n1 = R(1.0)/((s1 > R(1.0e-30)) ? s1 : R(1.0e-30));
    for (int k = 0; k < n; ++k) { z0[k] *= n0; z1[k] *= n1; }
    const R g = R(1.0) - in.f;
    const R rgf = in.f*rgas0 + g*rgas1;
    R jm[THERMO_MAX_SPECIES];
    R Sm = R(0.0);
    for (int k = 0; k < n; ++k) {
        jm[k] = rgf*in.D[k]*(z1[k] - z0[k])*in.geo;
        Sm += jm[k];
    }
    const int up0 = (Sm >= R(0.0)) ? 1 : 0;
    const R* zc = up0 ? z0 : z1;
    d.rgf = rgf; d.Sm = Sm; d.up0 = up0;
    if (iw >= 0) {
        const R ir0 = R(1.0)/tp_max<R>(in.rho0, R(1.0e-30)), ir1 = R(1.0)/tp_max<R>(in.rho1, R(1.0e-30));
        const R ctg = in.ct*in.geo;
        d.jv0 = jm[iw];
        d.jv_corr = zc[iw]*Sm;
        d.jv_mol = jm[iw] - zc[iw]*Sm;   // tp_face_flux の jc と同じ式 (FMA の縮約もコンパイラ任せで同じ形)
        const R v0 = (in.rY0[iw] - in.rg0)*ir0, v1 = (in.rY1[iw] - in.rg1)*ir1;
        d.jv_turb = ctg*(v1 - v0);
        d.zv0 = z0[iw]; d.zv1 = z1[iw]; d.zupv = zc[iw];
    } else {
        d.jv0 = d.jv_corr = d.jv_mol = d.jv_turb = R(0); d.zv0 = d.zv1 = d.zupv = R(0);
    }
}

// G1 (ii) の誤差尺度 A_l・A_v (double 参照の入力から)。z は tp_face_flux_diag と同じ正規化を double でやり直す。
TP_HD inline void tp_face_err_scales(const TpFaceInT<double>& in, double& A_l, double& A_v)
{
    const int n = in.n, iw = in.iw;
    const double r0 = (in.rho0 > 1.0e-30) ? in.rho0 : 1.0e-30, r1 = (in.rho1 > 1.0e-30) ? in.rho1 : 1.0e-30;
    const double ag = (in.ct*in.geo < 0.0) ? -(in.ct*in.geo) : in.ct*in.geo;
    const double g0 = in.rg0/r0, g1 = in.rg1/r1;
    A_l = ag*(((g0 < 0.0) ? -g0 : g0) + ((g1 < 0.0) ? -g1 : g1));
    A_v = 0.0;
    if (iw < 0) return;
    // z (tp_face_flux と同じ定義; double)
    const double rgas0 = tp_max<double>(in.rho0 - in.rg0, 1.0e-30), rgas1 = tp_max<double>(in.rho1 - in.rg1, 1.0e-30);
    double z0[THERMO_MAX_SPECIES], z1[THERMO_MAX_SPECIES], s0 = 0.0, s1 = 0.0;
    for (int k = 0; k < n; ++k) {
        z0[k] = ((k == iw) ? (in.rY0[k] - in.rg0) : in.rY0[k])/rgas0;
        z1[k] = ((k == iw) ? (in.rY1[k] - in.rg1) : in.rY1[k])/rgas1;
        s0 += z0[k]; s1 += z1[k];
    }
    const double n0 = 1.0/((s0 > 1.0e-30) ? s0 : 1.0e-30), n1 = 1.0/((s1 > 1.0e-30) ? s1 : 1.0e-30);
    double Sm = 0.0, sumA = 0.0;
    const double rgf = in.f*rgas0 + (1.0 - in.f)*rgas1;
    for (int k = 0; k < n; ++k) {
        z0[k] *= n0; z1[k] *= n1;
        const double c = tp_abs<double>(rgf*in.D[k]*in.geo);
        Sm += rgf*in.D[k]*(z1[k] - z0[k])*in.geo;
        sumA += c*(tp_abs<double>(z0[k]) + tp_abs<double>(z1[k]));
    }
    const double zupv = (Sm >= 0.0) ? z0[iw] : z1[iw];
    A_v = tp_abs<double>(rgf*in.D[iw]*in.geo)*(tp_abs<double>(z0[iw]) + tp_abs<double>(z1[iw])) + tp_abs<double>(zupv)*sumA;
}

// 監査 #4f (speciesTransport_d.cu twophase_audit_face_d) と同じ昇格: 本番の float 入力 (格納値と float 係数) をそのまま double に上げる。
// 差・除算・正規化は tp_face_flux<double> の中で double でやり直される。
TP_HD inline void tp_face_in_to_double(const TpFaceIn& in, TpFaceInT<double>& d)
{
    d.n = in.n; d.iw = in.iw;
    d.rho0 = in.rho0; d.rho1 = in.rho1; d.rg0 = in.rg0; d.rg1 = in.rg1; d.f = in.f; d.geo = in.geo; d.geo_abs = in.geo_abs;
    d.ct = in.ct; d.L = in.L;
    for (int s = 0; s < in.n; ++s) { d.rY0[s] = in.rY0[s]; d.rY1[s] = in.rY1[s]; d.D[s] = in.D[s]; d.h[s] = in.h[s]; }
    for (int m = 0; m < TP_NQ; ++m) { d.rQ0[m] = in.rQ0[m]; d.rQ1[m] = in.rQ1[m]; }
}
