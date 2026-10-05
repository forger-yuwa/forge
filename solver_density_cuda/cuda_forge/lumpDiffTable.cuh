#pragma once
// =============================================================================
// lumpDiffTable.cuh
//   lump を含む化学種拡散の縮約 (thermo_Dmix_lumped_f) の二元係数の表引き。
//   plan thermophysics-solver-owned-species-db §4.4 確定版・§5.1 #7c、判定 §6 V4e、設計の根拠 notes/reviews/2026-10-05-lump-blanc-diffusion-diagnose.md。
//
//   表: 実種の組 (r < q) ごとに f(T) = D_rq·P [m²/s·Pa] (圧力によらない) の ln f を ln T の区分 3 次 Hermite で持つ。
//   分割: Neufeld Ω(1,1) のクランプ点 T* = 0.3・100 (T = 0.3ε_rq・100ε_rq) で最大 3 区間に分け、各区間の中を ln T 等間隔
//         (Δln T ≤ 1/LUMPDIFF_TAB_PER_LNT) に刻む。値と ln T 微分は double の式 (区間の式) で作る。
//   範囲: [LUMPDIFF_TAB_TMIN, LUMPDIFF_TAB_TMAX]。外は式 (thermo_Dbinary_raw_f) へ戻す (端でクランプしない)。
//   評価: 区間は元の T (float) とクランプ点の T で選び、位置は面ごとに 1 回の double の ln T から作る。f = f0·expf(p)、
//         p = u(c1 + u(c2 + u c3)) は |p| ≲ 1.5 Δln T の小さい量 (ln f をそのまま float で持つと |ln f| の丸めが乗る)。
// =============================================================================
#include "cuda_forge/thermo_d.cuh"

#define LUMPDIFF_TAB_TMIN     50.0
#define LUMPDIFF_TAB_TMAX     30000.0
#define LUMPDIFF_TAB_PER_LNT  128
#define LUMPDIFF_MAX_PAIRS    (THERMO_MAX_DIFF_REAL*(THERMO_MAX_DIFF_REAL-1)/2)

struct LumpDiffC4 { float x, y, z, w; };   // x = f0 (小区間左端の f)、y..w = ln(f/f0) の 3 次式の係数 (float4 は host 単体試験で使えないことがあるので独自型)

struct LumpDiffPairTab {
    int   nseg;
    float Tup[3];       // 区間 k の上端 T (最後は TMAX)
    double lnTa[3];     // 区間 k の下端 ln T (double: 位置の丸めが ln f の誤差に効くため)
    double invH[3];     // 1/Δln T
    int   m[3];         // 小区間数
    int   off[3];       // 係数配列の先頭
};

struct LumpDiffTabD {
    int                     on;
    const LumpDiffPairTab*  pair;   // [nPairs] (r < q の通し番号 r*nr − r(r+1)/2 + q − r − 1)
    const LumpDiffC4*       coef;   // ln f = c.x + u(c.y + u(c.z + u c.w))
    float                   Tmin, Tmax;
};

THERMO_HD int lumpdiff_pair_index(int r, int q, int nr) { return r*nr - r*(r + 1)/2 + (q - r - 1); }

// 表で D_rq を返す (範囲外は式)。lnT は double の log(T)。
THERMO_HD float lumpdiff_tab_D(const LumpDiffTabD& tb, const LumpDiffD& ld, int r, int q, float T, double lnT, float P)
{
    if (!(T >= tb.Tmin && T <= tb.Tmax))
        return thermo_Dbinary_raw_f(ld.MW[r], ld.sig[r], ld.eps[r], ld.MW[q], ld.sig[q], ld.eps[q], T, P);
    const LumpDiffPairTab& pt = tb.pair[lumpdiff_pair_index(r, q, ld.nReal)];
    int k = 0;
    while (k < pt.nseg - 1 && T > pt.Tup[k]) ++k;
    const double x = (lnT - pt.lnTa[k])*pt.invH[k];
    int i = (int)x;
    if (x < 0.0) i = 0;
    if (i > pt.m[k] - 1) i = pt.m[k] - 1;
    const float u = (float)(x - (double)i);
    const LumpDiffC4 c = tb.coef[pt.off[k] + i];
    return c.x*expf(u*(c.y + u*(c.z + u*c.w))) / P;
}

// thermo_Dmix_lumped_f の表引き版 (和の順序・純成分の扱いは同じ)。tb.on == 0 なら式の版へ。
THERMO_HD void thermo_Dmix_lumped_tab_f(const LumpDiffD& ld, const LumpDiffTabD& tb, int nLabel, const float* XL, float T, float P, float* Dl)
{
    if (!tb.on) { thermo_Dmix_lumped_f(ld, nLabel, XL, T, P, Dl); return; }
    if (nLabel <= 1) { for (int s = 0; s < nLabel; ++s) Dl[s] = 0.0f; return; }
    const int nr = ld.nReal;
    const double lnT = log((double)T);   // 面ごとに 1 回
    float Xr[THERMO_MAX_DIFF_REAL], num[THERMO_MAX_DIFF_REAL], den[THERMO_MAX_DIFF_REAL];
    for (int r = 0; r < nr; ++r) {
        float x = 0.0f;
        for (int s = 0; s < nLabel; ++s) x += XL[s]*ld.E[s][r];
        Xr[r] = x; num[r] = 0.0f; den[r] = 0.0f;
    }
    for (int r = 0; r < nr; ++r) {
        for (int q = r + 1; q < nr; ++q) {
            float d = lumpdiff_tab_D(tb, ld, r, q, T, lnT, P);
            d = (d > 1.0e-30f ? d : 1.0e-30f);
            num[r] += Xr[q]; den[r] += Xr[q]/d;
            num[q] += Xr[r]; den[q] += Xr[r]/d;
        }
    }
    for (int r = 0; r < nr; ++r) {
        num[r] = (den[r] < 1.0e-30f)
            ? thermo_Dbinary_raw_f(ld.MW[r], ld.sig[r], ld.eps[r], ld.MW[r], ld.sig[r], ld.eps[r], T, P)
            : num[r]/den[r];
    }
    for (int s = 0; s < nLabel; ++s) {
        float D = 0.0f;
        for (int r = 0; r < nr; ++r) if (ld.A[s][r] > 0.0f) D += ld.A[s][r]*num[r];
        Dl[s] = D;
    }
}

#include <cmath>
#include <vector>
// ---- host: 表の構築 (double の式) ----
namespace lumpdiff_detail {
inline double omega(double Ts) {
    return 1.06036*std::pow(Ts, -0.15610) + 0.19300*std::exp(-0.47635*Ts) + 1.03587*std::exp(-1.52996*Ts) + 1.76474*std::exp(-3.89411*Ts); }
inline double dlnOmega_dlnTs(double Ts) {
    const double d = -1.06036*0.15610*std::pow(Ts, -0.15610) - 0.19300*0.47635*Ts*std::exp(-0.47635*Ts)
                     - 1.03587*1.52996*Ts*std::exp(-1.52996*Ts) - 1.76474*3.89411*Ts*std::exp(-3.89411*Ts);
    return d/omega(Ts); }
// ln f と d ln f/d ln T (f = D·P [m²/s·Pa])。clampSide: -1 = T* ≤ 0.3 の区間、0 = 中間、+1 = T* ≥ 100 の区間 (区間の式で評価)
inline void lnf(double Ma, double sa, double ea, double Mb, double sb, double eb, double T, int clampSide, double* v, double* dv) {
    const double sig = 0.5*(sa + sb), eps = std::sqrt(ea*eb), Ts = T/eps;
    const double Tc = (clampSide < 0) ? 0.3 : (clampSide > 0 ? 100.0 : Ts);
    const double C = 1.8583e-3*1.0e-4*101325.0*std::sqrt(1.0/(Ma*1e3) + 1.0/(Mb*1e3))/(sig*sig);
    *v  = std::log(C) + 1.5*std::log(T) - std::log(omega(Tc));
    *dv = 1.5 - ((clampSide == 0) ? dlnOmega_dlnTs(Ts) : 0.0);
}
}  // namespace lumpdiff_detail

// ld (実種の MW・σ・ε) から表を作る。pairs・coef は呼び出し側が device へ上げる。
inline void lumpdiff_build_table(const LumpDiffD& ld, std::vector<LumpDiffPairTab>& pairs, std::vector<LumpDiffC4>& coef)
{
    using namespace lumpdiff_detail;
    const int nr = ld.nReal;
    pairs.assign(nr > 1 ? nr*(nr - 1)/2 : 0, LumpDiffPairTab{});
    coef.clear();
    for (int r = 0; r < nr; ++r) for (int q = r + 1; q < nr; ++q) {
        LumpDiffPairTab& pt = pairs[lumpdiff_pair_index(r, q, nr)];
        const double eps = std::sqrt((double)ld.eps[r]*(double)ld.eps[q]);
        double cuts[4]; int side[3]; int n = 0;
        double lo = LUMPDIFF_TAB_TMIN;
        const double k1 = 0.3*eps, k2 = 100.0*eps;
        const double bnd[2] = {k1, k2}; const int sd[3] = {-1, 0, +1};
        cuts[0] = lo;
        int sidx = (lo < k1) ? 0 : (lo < k2 ? 1 : 2);
        for (int b = sidx; b < 2; ++b) {
            if (bnd[b] >= LUMPDIFF_TAB_TMAX) break;
            side[n] = sd[b]; cuts[++n] = bnd[b]; sidx = b + 1;
        }
        side[n] = sd[sidx]; cuts[++n] = LUMPDIFF_TAB_TMAX;
        pt.nseg = n;
        for (int k = 0; k < n; ++k) {
            const double a = std::log(cuts[k]), b = std::log(cuts[k + 1]);
            int m = (int)std::ceil((b - a)*LUMPDIFF_TAB_PER_LNT); if (m < 1) m = 1;
            const double h = (b - a)/m;
            pt.Tup[k] = (float)cuts[k + 1]; pt.lnTa[k] = a; pt.invH[k] = 1.0/h; pt.m[k] = m; pt.off[k] = (int)coef.size();
            for (int i = 0; i < m; ++i) {
                double f0, d0, f1, d1;
                lnf(ld.MW[r], ld.sig[r], ld.eps[r], ld.MW[q], ld.sig[q], ld.eps[q], std::exp(a + i*h), side[k], &f0, &d0);
                lnf(ld.MW[r], ld.sig[r], ld.eps[r], ld.MW[q], ld.sig[q], ld.eps[q], std::exp(a + (i + 1)*h), side[k], &f1, &d1);
                LumpDiffC4 c;
                c.x = (float)std::exp(f0); c.y = (float)(d0*h);
                c.z = (float)(3.0*(f1 - f0) - (2.0*d0 + d1)*h);
                c.w = (float)(2.0*(f0 - f1) + (d0 + d1)*h);
                coef.push_back(c);
            }
        }
    }
}
