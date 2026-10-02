#pragma once
// =============================================================================
// twoPhaseDiffusion_d.cuh — 凝縮 TP carrier の二相拡散 (面流束の共有) と蒸気/液の非分割更新 (定常専用初版)
//   plans/active/condensation-two-phase-transport.md §4.2・§5.1 #4e。仕様 methods/condensation.md §7c。
//   設計の記録 notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md §14 (§6.1 が更新 B の定義)。
//
// 面ごとの代数 (tp_face_flux) と 1 セルの更新 (tp_vl_update) を __host__ __device__ で置き、本番カーネル
// (speciesTransport_d.cu twophase_diffusion_d / condensationTransport_d.cu の更新) と GPU 単体試験
// (tests/unit/test_twophase_kernel.cu) が同じ関数を呼ぶ。演算はすべて float32 (格納と同じ精度)。
//
// 記号 (§4.2): 気相密度 ρ_g = ρ − ρg、気相内の質量分率 z_k = ρY_k/ρ_g (k ≠ 水)、z_v = (ρY_w − ρg)/ρ_g (Σz で正規化)。
//   分子 (気相のみ): j_k⁰ = −ρ_g,f D_k ∇z_k、補正 j_k = j_k⁰ − z_k,up Σj⁰ (z_up は補正流束 −Σj⁰ の流出側のセル)。
//   乱流 (全輸送量共通): J_t(φ) = −(μ_t/Sc_t) ∇φ (φ は質量当たり: Y_k, Y_w − g, g, Q_n/ρ)。
//   J_k = j_k + J_t(Y_k) (k ≠ 水)、J_v = j_v + J_t(Y_w − g)、J_l = J_t(g)、J_w = J_v + J_l、J_Qn = J_t(Q_n/ρ)。
//   エネルギー Σ_{k≠水} h_k J_k + h_v J_w − L(T_f) J_l (h_l = h_v − L; EOS e = e_gas + g(R_w T − L) と同じ基準)。
// 向き: 既存の species_diffusion_d と同じく「セル 0 へ入る向き」が正 (res[ic0] += J, res[ic1] −= J)。
// =============================================================================
#include "cuda_forge/thermo_d.cuh"   // THERMO_MAX_SPECIES

#ifndef TP_HD
#ifdef __CUDACC__
#define TP_HD __host__ __device__
#else
#define TP_HD
#endif
#endif

#define TP_NQ 3   // 液のモーメント Q2, Q1, Q0 (ρg とは別に乱流拡散する)

// 面の入力 (セル 0/1 の格納値と面の係数)。D[k] は気相内の分子拡散係数 D_k (乱流分を含まない)。
// R = float が本番 (TpFaceIn/TpFaceOut)。R = double は収束受入の独立残差監査 (#4f) が格納値を double に上げて同じ式を評価する。
template <typename R> struct TpFaceInT {
    int   n, iw;                                 // 化学種数、凝縮種 (水) の輸送種 index
    R rho0, rho1;                            // セルの ρ
    R rY0[THERMO_MAX_SPECIES], rY1[THERMO_MAX_SPECIES];   // ρY_k (水は総水分 ρY_w)
    R rg0, rg1;                              // 液 ρg
    R rQ0[TP_NQ], rQ1[TP_NQ];                // ρQ_n (順序 Q2, Q1, Q0)
    R f;                                     // 面の内挿重み (セル 0 側)
    R geo, geo_abs;                          // δ/dcc (流束の係数) と |δ|/dcc (対角)
    R D[THERMO_MAX_SPECIES];                 // D_k [m²/s]
    R ct;                                    // 乱流の係数 μ_t,f/Sc_t
    R h[THERMO_MAX_SPECIES];                 // h_k(T_f) [J/kg] (水は蒸気 h_v)
    R L;                                     // 潜熱 L(T_f) [J/kg]
};
using TpFaceIn = TpFaceInT<float>;

// 面の出力。J[k] は化学種の残差に足す流束 (水は総水分 J_w)。diag*[k] は輸送の点対角 (ρY_k 単位; 水は**蒸気**の対角)。
template <typename R> struct TpFaceOutT {
    R J[THERMO_MAX_SPECIES];
    R Jv, Jl, JQ[TP_NQ];
    R q;                                     // エネルギー流束
    R diag0[THERMO_MAX_SPECIES], diag1[THERMO_MAX_SPECIES];
    R diagt0, diagt1;                        // 液・Q の対角 (乱流のみ)
    R Sm;                                    // 補正前の分子流束の和 Σ j⁰ (セル 0 へ入る向き)
    int   up0;                                   // 補正の z をセル 0 から取ったか
};
using TpFaceOut = TpFaceOutT<float>;

template <typename R> TP_HD inline R tp_max(R a, R b) { return (a > b) ? a : b; }   // fmaxf と同値 (NaN を除く)
template <typename R> TP_HD inline R tp_abs(R a) { return (a < R(0)) ? -a : a; }

template <typename R>
TP_HD inline void tp_face_flux(const TpFaceInT<R>& in, TpFaceOutT<R>& o)
{
    const int n = in.n, iw = in.iw;
    // 気相組成 z (セルごと; Σz で正規化)
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

    // 分子流束 j⁰ と和
    R jm[THERMO_MAX_SPECIES];
    R Sm = R(0.0);
    for (int k = 0; k < n; ++k) {
        jm[k] = rgf*in.D[k]*(z1[k] - z0[k])*in.geo;
        Sm += jm[k];
    }
    // 補正の質量流束 −Σj⁰ (セル 0 へ入る向き) が負 = セル 0 から出る → セル 0 の z (風上)
    const int up0 = (Sm >= R(0.0)) ? 1 : 0;
    const R* zc = up0 ? z0 : z1;
    const R ir0 = R(1.0)/tp_max<R>(in.rho0, R(1.0e-30)), ir1 = R(1.0)/tp_max<R>(in.rho1, R(1.0e-30));
    const R ctg = in.ct*in.geo;

    R q = R(0.0);
    for (int k = 0; k < n; ++k) {
        const R jc = jm[k] - zc[k]*Sm;
        if (k == iw) {
            const R v0 = (in.rY0[k] - in.rg0)*ir0, v1 = (in.rY1[k] - in.rg1)*ir1;
            o.Jv = jc + ctg*(v1 - v0);
        } else {
            o.J[k] = jc + ctg*(in.rY1[k]*ir1 - in.rY0[k]*ir0);
            q += in.h[k]*o.J[k];
        }
        // 対角: 分子 (気相) + 乱流 + 補正の流出側
        const R dm = rgf*in.D[k]*in.geo_abs;
        o.diag0[k] = dm*ig0 + in.ct*in.geo_abs*ir0;
        o.diag1[k] = dm*ig1 + in.ct*in.geo_abs*ir1;
    }
    const R aS = tp_abs<R>(Sm);
    for (int k = 0; k < n; ++k) {
        if (up0) o.diag0[k] += aS*ig0; else o.diag1[k] += aS*ig1;
    }
    o.Jl = ctg*(in.rg1*ir1 - in.rg0*ir0);
    if (iw >= 0) {
        o.J[iw] = o.Jv + o.Jl;
        q += in.h[iw]*o.J[iw] - in.L*o.Jl;
    }
    for (int m = 0; m < TP_NQ; ++m) o.JQ[m] = ctg*(in.rQ1[m]*ir1 - in.rQ0[m]*ir0);
    o.q = q;
    o.diagt0 = in.ct*in.geo_abs*ir0;
    o.diagt1 = in.ct*in.geo_abs*ir1;
    o.Sm = Sm;
    o.up0 = up0;
}

// 1 セルの非分割更新 (設計メモ §6.1 の vl_limit_commit; #4c の C++ ハーネス test_twophase_real_source.cpp と同じ演算順)。
//   増分 (前処理の後): δρv = (R_w − R_g)/D_v、δρg = R_g/(D_g + V sj_g)、δρQ_n = R_Qn/(D_Qn + V sj_Qn)  (D_* = V/Δτ + 輸送の点対角)
//   緩和 ω を全増分に掛けてから制限 θ = min(θ_thr, θ_vg) を共通に掛け、Q は成分ごとに非負化、総水分は ρY_w + fl(θδρv + θδρg)。
struct TpCellIn {
    float M, V;                  // V/Δτ, V
    float Rw, Rg, RQ[TP_NQ];     // 全残差 (総水分にソースは無い; 液・Q はソース込み)
    float Dv, Dg, DQ[TP_NQ];     // 輸送の点対角 (V/Δτ を含まない)
    float sjg, sjQ[TP_NQ];       // ソースの点 Jacobian (体積当たり)
    float rYw, rg, rQ[TP_NQ];    // 更新前の格納値
    float rho;                   // 質量分率の換算に使う ρ (流れの更新後)
    float omega;                 // 緩和 (condTwoPhaseRelax)
    double dg_max, dT_max;       // condDgMaxStep / condDTmaxStep
    double L, cveff;             // 潜熱と有効定積比熱 (θ の ΔT 換算; 既存の更新クランプと同じ式で呼び出し側が作る)
};
struct TpCellOut {
    float rYw, rg, rQ[TP_NQ];
    double theta;
    double withheld_v, withheld_g;   // (1−θ)|δ|
    double qcut, vround;             // 状態を書き換えた補正 (Q の非負化、丸めによる蒸気の負)
};

// thetaRound: θ (double) を float にするときの丸め。1 = 切り上がったら 0 側の隣の float へ (安全側; 既定)、0 = 最近接 (#4e; 判別用)。
//   #4f (1) の 1 セル反例 (rYw 0.01f, rg 3e-6f, Rg −1.45e-4f): 最近接は θ_vg を 1 ulp 越えて液 ρg = −2.27e-13、安全側は +2.27e-13・補正 0 → 安全側を採用 (2026-10-02)。
//   安全側なら |fl(Th·δρg)| ≤ ρg (丸めは単調) なので ρg + fl(Th·δρg) ≥ 0 が成り立つ。
#ifndef TP_THETA_ROUND_DEFAULT
#define TP_THETA_ROUND_DEFAULT 1
#endif
TP_HD inline void tp_vl_update(const TpCellIn& c, TpCellOut& o, int thetaRound = TP_THETA_ROUND_DEFAULT)
{
    const float rv = c.Rw - c.Rg;                      // 全残差変換 (float の減算)
    float dv = rv/(c.M + c.Dv);
    float dg = c.Rg/(c.M + c.Dg + c.V*c.sjg);
    float dq[TP_NQ];
    for (int m = 0; m < TP_NQ; ++m) dq[m] = c.RQ[m]/(c.M + c.DQ[m] + c.V*c.sjQ[m]);
    if (c.omega != 1.0f) { dv *= c.omega; dg *= c.omega; for (int m = 0; m < TP_NQ; ++m) dq[m] *= c.omega; }
    // θ_thr (全増分に対する閾値) と θ_vg (蒸気・液の非負; 蒸気は更新前の ρv)
    double th = 1.0;
    const double rho = (double)c.rho;
    const double adg = (rho > 0.0) ? fabs((double)dg)/rho : 0.0;
    if (adg > 0.0) {
        th = fmin(th, c.dg_max/adg);
        const double adT = adg*c.L/c.cveff;
        if (adT > 0.0) th = fmin(th, c.dT_max/adT);
    }
    const double rvs = (double)c.rYw - (double)c.rg;
    if ((double)dv < 0.0) th = fmin(th, rvs/(-(double)dv));
    if ((double)dg < 0.0) th = fmin(th, (double)c.rg/(-(double)dg));
    if (!(th > 0.0)) th = 0.0;   // NaN も 0 (動かさない; 残差が下がらないことで監視に出る)
    float Th = (float)th;
    if (thetaRound != 0 && (double)Th > th) Th = nextafterf(Th, 0.0f);   // 共通の実効 θ を安全側へ (θ_vg の境界を越えない)
    o.theta = th;
    o.withheld_v = (1.0 - th)*fabs((double)dv);
    o.withheld_g = (1.0 - th)*fabs((double)dg);
    o.qcut = 0.0; o.vround = 0.0;
    for (int m = 0; m < TP_NQ; ++m) {
        float d = Th*dq[m];
        const float nq = c.rQ[m] + d;
        if (nq < 0.0f) { o.qcut += -(double)nq; d = -c.rQ[m]; }
        o.rQ[m] = c.rQ[m] + d;
    }
    const float gnew = c.rg + Th*dg;
    float wnew = c.rYw + (Th*dv + Th*dg);
    if (wnew - gnew < 0.0f) { o.vround += (double)(gnew - wnew); wnew = gnew; }
    o.rg = gnew; o.rYw = wnew;
}
