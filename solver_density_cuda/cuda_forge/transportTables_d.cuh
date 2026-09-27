#pragma once

// =============================================================================
// transportTables_d.cuh
//   種ごとの輸送物性 (physProp.transport) の表引き (float)。plans/active/thermophysics-solver-owned-species-db.md §5.1 #5t2-3、
//   設計の根拠 notes/reviews/2026-09-27-transport-tables-diagnose.md、仕様 methods/thermophysics.md。
//   前例は凝縮の表 cuda_forge/condensationTables_d.cuh (plans/accepted/condensation-float-speedup.md)。
//
//   表: 実種の ln μᵢ・ln λᵢ と、組 (両 kinetic の二元 Chapman–Enskog・CEA 相互作用) の ln ηᵢⱼ。剛体球の組は表を持たず、
//       実行時に種別表の μ から作る (transportMix_d.cuh の transport_eta_pair と同じ式)。
//   分割: 式の区間境界 (CEA・fit の区間、H2O の接続点 253.15/500/700 K、kinetic の T* クランプ点 0.3ε・100ε、
//         修正 Eucken の c_p の NASA Tlo/Tmid/Thi) で表を分割し、各分割区間の中を ln T 等間隔 (Δln T ≤ 1/256) に刻む。
//         各小区間は両端の値と ln T 微分から作る 3 次 Hermite。値・微分は**その分割区間の式**で評価する
//         (右の区間の左端でも右側の式を使う; 式は区間の中点で選び、端点ではその式を評価するだけ)。
//   選択: 区間の所属は**元の T と元の境界値**で決める (現行の T ≤ Thi / T < Tmid などの規約そのまま)。
//         float の ln T では選ばない: 999.99994/1000/1000.00006 K の ln T は float で同じ値になり、段差のある fit で左右を取り違える
//         (codex 2026-09-27 の算術反例)。境界 b と規約 (T ≤ b が左 / T < b が左) から「左の区間に属する最大の float」を
//         host で作り (TransportSegF::Tupper)、float の T と比べる (元の double 比較と同じ所属になる)。
//   範囲: 150–15000 K。外は段 2 の double 評価 (transport_mix_Y) へ委譲する (端の値でクランプしない)。
//   評価: 小区間の番号と位置 u は double の ln T から作り (1 セル 1 回の log と、表ごとに減算・乗算 1 回)、
//         係数・Hermite・exp・混合則は float。exp(p) は |p| ≤ 0.05 を保証して 4 次までの級数で評価する
//         (打ち切り ≤ 3e-9; 超える区間は構築時に細分する)。
//   組成: 輸送種の X_s ∝ max(ρY_s, 0)/M_s を float で作ってから展開行列を掛ける (段 2 と同じモル基底; 正規化の順は違うが値は同じ式)。
// =============================================================================

#include "cuda_forge/transportMix_d.cuh"

#define TRANSPORT_TAB_TMIN      150.0      // 表の下端 [K]
#define TRANSPORT_TAB_TMAX      15000.0    // 表の上端 [K]
#define TRANSPORT_TAB_PER_LNT   256        // 1 (ln T) あたりの小区間数の初期値 (Δln T ≤ 1/256)
#define TRANSPORT_TAB_PMAX      0.05       // 小区間内の |ln f − ln f0| の上限 (exp の級数の打ち切りを抑える)

// -----------------------------------------------------------------------------
// 評価 (host / device)
// -----------------------------------------------------------------------------

// 表の分割区間を元の T の所属規約で選ぶ (区間の数は多くても十数個なので線形探索)
THERMO_HD int transport_tab_seg(const TransportTablesF& tb, const TransportTabRefF& r, float T)
{
    int k = r.seg0;
    const int e = r.seg0 + r.nseg - 1;
    while (k < e && T > tb.seg[k].Tupper) ++k;
    return k;
}

// 分割区間 s の中の小区間番号 (係数配列の通し番号) と位置 u
THERMO_HD int transport_tab_locate(const TransportSegF& s, double lnT, float* u)
{
    const double x = (lnT - s.lnTa)*s.invH;
    int i = (int)x;
    if (x < 0.0) i = 0;             // 所属は T で決めてあるので、ここに来るのは ln の丸め (1 ulp 級) だけ
    if (i > s.m - 1) i = s.m - 1;
    *u = (float)(x - (double)i);
    return s.off + i;
}

THERMO_HD float transport_tab_hermite(const TransportHermiteF& c, float u)
{
    const float p = u*(c.c1 + u*(c.c2 + u*c.c3));
    return c.f0*(1.0f + p*(1.0f + p*(0.5f + p*(1.0f/6.0f + p*(1.0f/24.0f)))));
}

// 実種 r の μ [Pa s]・λ [W/(m K)] (範囲内の T を渡す)
THERMO_HD void transport_tab_species(const TransportTablesF& tb, int r, float T, double lnT, float* mu, float* lam)
{
    const TransportSegF& s = tb.seg[transport_tab_seg(tb, tb.spTab[r], T)];
    float u;
    const int j = transport_tab_locate(s, lnT, &u);
    *mu  = transport_tab_hermite(tb.spc[2*j],     u);
    *lam = transport_tab_hermite(tb.spc[2*j + 1], u);
}

// 組 (a < b, 通し番号 ip) の η_ab [Pa s]。eta_a, eta_b は種別表の μ (剛体球で使う)。
THERMO_HD float transport_tab_pair(const TransportTablesF& tb, int ip, float T, double lnT, float eta_a, float eta_b)
{
    const TransportTabRefF& r = tb.pairTab[ip];
    if (r.nseg == 0) {   // 剛体球 (cea2.f 5565–5570; transport_eta_pair と同じ式)
        const TransportPairMixF& pm = tb.pm[ip];
        const float d = 1.0f + sqrtf(pm.rsRatio*eta_a/eta_b);
        return pm.rsC*eta_a/(d*d);
    }
    const TransportSegF& s = tb.seg[transport_tab_seg(tb, r, T)];
    float u;
    const int j = transport_tab_locate(s, lnT, &u);
    return transport_tab_hermite(tb.pairc[j], u);
}

// CEA frozen 混合則 (float; 範囲内の T)。X: 実種のモル分率 (float, n 個)。
//   規則は transport_mix と同じ (X_i ≤ 0 の種は分子にも分母にも入れない、相手の X が 0 の組は η を評価しない)。
THERMO_HD void transport_mix_tab(const TransportTableD& t, const float* X, float T, float* mu, float* lam)
{
    const TransportTablesF& tb = t.tab;
    const int n = t.nReal;
    const double lnT = log((double)T);
    float eta[TRANSPORT_MAX_REAL_SPECIES], con[TRANSPORT_MAX_REAL_SPECIES];
    float sv[TRANSPORT_MAX_REAL_SPECIES], sc[TRANSPORT_MAX_REAL_SPECIES];
    for (int i = 0; i < n; ++i) {
        if (X[i] != 0.0f) transport_tab_species(tb, i, T, lnT, &eta[i], &con[i]);
        else { eta[i] = 0.0f; con[i] = 0.0f; }
        sv[i] = 0.0f;
        sc[i] = 0.0f;
    }
    for (int a = 0; a < n; ++a) {
        const bool pa = (X[a] > 0.0f);
        if (pa) { sv[a] += X[a]; sc[a] += X[a]; }   // φ_aa = ψ_aa = 1
        for (int b = a + 1; b < n; ++b) {
            const bool pb = (X[b] > 0.0f);
            const bool ua = pa && X[b] != 0.0f;
            const bool ub = pb && X[a] != 0.0f;
            if (!ua && !ub) continue;
            const int ip = transport_pair_index(a, b, n);
            const TransportPairMixF& pm = tb.pm[ip];
            const float ie = 1.0f/transport_tab_pair(tb, ip, T, lnT, eta[a], eta[b]);
            if (ua) {   // i = a, j = b
                const float phi = pm.kab*eta[a]*ie;
                sv[a] += phi*X[b];
                sc[a] += phi*pm.gab*X[b];
            }
            if (ub) {   // i = b, j = a
                const float phi = pm.kba*eta[b]*ie;
                sv[b] += phi*X[a];
                sc[b] += phi*pm.gba*X[a];
            }
        }
    }
    float m = 0.0f, l = 0.0f;
    for (int i = 0; i < n; ++i) {
        if (!(X[i] > 0.0f)) continue;
        m += eta[i]*X[i]/sv[i];
        l += con[i]*X[i]/sc[i];
    }
    *mu = m;
    *lam = l;
}

// セル・壁の共通入口 (表引き)。rY: 輸送種の ρY_s または Y_s (正規化不要; 負値は 0 に切る)。nY < 2 は単成分 (X = {1})。
//   th は範囲外で段 2 の double 評価へ委譲するときだけ使う (そのとき Y は Σ で正規化して double で渡す)。
THERMO_HD void transport_mix_Y_tab(const SpeciesThermo* th, const TransportTableD& t, int nY, const float* rY, float T,
                                   float* mu, float* lam)
{
    const TransportTablesF& tb = t.tab;
    if (!(T >= tb.Tmin && T <= tb.Tmax)) {   // 範囲外 (と NaN) は段 2 の double 評価
        double Y[THERMO_MAX_SPECIES];
        int nd = 1;
        if (nY >= 2) {
            double s = 0.0;
            for (int k = 0; k < nY; ++k) { Y[k] = (rY[k] > 0.0f) ? (double)rY[k] : 0.0; s += Y[k]; }
            const double inv = 1.0/(s > 1.0e-30 ? s : 1.0e-30);
            for (int k = 0; k < nY; ++k) Y[k] *= inv;
            nd = nY;
        }
        double md, ld;
        transport_mix_Y(th, t, nd, Y, (double)T, &md, &ld);
        *mu = (float)md;
        *lam = (float)ld;
        return;
    }
    float Xs[THERMO_MAX_SPECIES];
    if (nY < 2) {
        Xs[0] = 1.0f;
        for (int s = 1; s < t.nTransported; ++s) Xs[s] = 0.0f;
    } else {
        float sum = 0.0f;
        for (int s = 0; s < nY; ++s) { const float x = (rY[s] > 0.0f ? rY[s] : 0.0f)*tb.invMWs[s]; Xs[s] = x; sum += x; }
        const float inv = 1.0f/(sum > 1.0e-30f ? sum : 1.0e-30f);
        for (int s = 0; s < nY; ++s) Xs[s] *= inv;
    }
    float Xr[TRANSPORT_MAX_REAL_SPECIES];
    for (int r = 0; r < t.nReal; ++r) Xr[r] = 0.0f;
    for (int s = 0; s < t.nTransported; ++s) {
        const float xs = Xs[s];
        if (xs == 0.0f) continue;
        const float* e = tb.expandF + s*t.nReal;
        for (int r = 0; r < t.nReal; ++r) Xr[r] += xs*e[r];
    }
    transport_mix_tab(t, Xr, T, mu, lam);
}

// 単体値 (試験用 FORGE_TRANSPORT_TABLE_PROBE)。kind 0: 実種 idx の (μ, λ)、kind 1: 組 idx の (η, 0)。
//   範囲外は段 2 の double の単成分式 (委譲先と同じ) を float に丸めて返す。
THERMO_HD void transport_tab_single(const TransportTableD& t, int kind, int idx, float T, float* v0, float* v1)
{
    const TransportTablesF& tb = t.tab;
    const bool in = (T >= tb.Tmin && T <= tb.Tmax);
    const double lnT = log((double)T);
    if (kind == 0) {
        if (in) { transport_tab_species(tb, idx, T, lnT, v0, v1); return; }
        double m, l;
        transport_species(t.sp[idx], (double)T, &m, &l);
        *v0 = (float)m; *v1 = (float)l;
        return;
    }
    const int a = t.pairs[idx].a, b = t.pairs[idx].b;
    *v1 = 0.0f;
    if (in) {
        float ma, la, mb, lb;
        transport_tab_species(tb, a, T, lnT, &ma, &la);
        transport_tab_species(tb, b, T, lnT, &mb, &lb);
        *v0 = transport_tab_pair(tb, idx, T, lnT, ma, mb);
        return;
    }
    double eta[TRANSPORT_MAX_REAL_SPECIES], dummy;
    transport_species(t.sp[a], (double)T, &eta[a], &dummy);
    transport_species(t.sp[b], (double)T, &eta[b], &dummy);
    *v0 = (float)transport_eta_pair(t.pairs[idx], t.sp, eta, (double)T);
}

// -----------------------------------------------------------------------------
// host 側の構築 (double)。thermo_d.cu と単体試験 (host) から使う。
// -----------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

struct TransportTablesHost {
    std::vector<TransportSegF>     seg;
    std::vector<TransportHermiteF> spc, pairc;
    std::vector<TransportTabRefF>  spTab, pairTab;
    std::vector<TransportPairMixF> pm;
    std::vector<float>             invMWs, expandF;
    std::vector<std::vector<double>> segLo, segHi;   // 表ごとの分割区間の端 (double; 記録・試験用)。index: 実種 r、組は nReal + ip
    float Tmin = 0.0f, Tmax = 0.0f;
    std::string error;                               // 空でなければ構築失敗 (段 2 の double 評価のまま)
    size_t bytes() const
    {
        return seg.size()*sizeof(TransportSegF) + (spc.size() + pairc.size())*sizeof(TransportHermiteF)
             + (spTab.size() + pairTab.size())*sizeof(TransportTabRefF) + pm.size()*sizeof(TransportPairMixF)
             + (invMWs.size() + expandF.size())*sizeof(float);
    }
};

// 境界 1 つ: T = b、leftIncl = 1 なら T ≤ b が左の式 (0 なら T < b が左)、disc = 1 なら値が不連続になりうる境界
struct TransportTabBound {
    double b;
    int leftIncl, disc;
};

// ---- 区間を固定した (Tsel が選ぶ式の) double 評価。値は transportMix_d.cuh の各式と同じ形で書く ----
inline int transport_tab_fit_k(const TransportFitD& f, double Tsel)
{
    int k = f.n - 1;
    for (int i = 0; i < f.n; ++i) { if (Tsel <= f.Thi[i]) { k = i; break; } }
    return k;
}

inline double transport_tab_fit_forced(const TransportFitD& f, double T, double Tsel)
{
    const int k = transport_tab_fit_k(f, Tsel);
    const double lnT = log(T), iT = 1.0/T;
    return exp(f.A[k]*lnT + (f.B[k] + f.C[k]*iT)*iT + f.D[k]);
}

// T* のクランプは Tsel で決め、クランプしない区間では T/ε を使う
inline double transport_tab_Tstar_forced(double T, double Tsel, double eps)
{
    const double s = Tsel/eps;
    if (s < 0.3)   return 0.3;
    if (s > 100.0) return 100.0;
    return T/eps;
}

inline double transport_tab_cp_mass_forced(const SpeciesThermo& sp, double T, double Tsel)
{
    double Tc = T, Tcs = Tsel;
    if (Tsel < sp.Tlo) { Tc = sp.Tlo; Tcs = sp.Tlo; }
    if (Tsel > sp.Thi) { Tc = sp.Thi; Tcs = sp.Thi; }
    const double* a = (Tcs < sp.Tmid) ? sp.low : sp.high;
    const double Ti = 1.0/Tc, Ti2 = Ti*Ti;
    const double cpm = THERMO_RU*(a[0]*Ti2 + a[1]*Ti + a[2] + a[3]*Tc + a[4]*Tc*Tc + a[5]*Tc*Tc*Tc + a[6]*Tc*Tc*Tc*Tc);
    return cpm/sp.MW;
}

inline void transport_tab_species_forced(const SpeciesTransportD& s, double T, double Tsel, double* mu, double* lam)
{
    switch (s.model) {
    case TRANSPORT_MODEL_KINETIC: {
        const double Ts = transport_tab_Tstar_forced(T, Tsel, s.eps_kB);
        const double om = transport_omega22(Ts) + 0.2*s.deltaStar*s.deltaStar/Ts;
        const double m = TRANSPORT_CE_MU_CONST*sqrt(s.MW*1000.0*T)/(s.sigma_LJ*s.sigma_LJ*om);
        *mu = m;
        *lam = m*(transport_tab_cp_mass_forced(s.thermo, T, Tsel) + 1.25*THERMO_RU/s.MW);
        return;
    }
    case TRANSPORT_MODEL_H2O_IAPWS_CEA_V1: {
        double mI, lI, nm, nl;
        if (Tsel >= TRANSPORT_H2O_BLEND_HI) {
            *mu  = transport_tab_fit_forced(s.V, T, Tsel)*TRANSPORT_MICROPOISE_TO_PAS;
            *lam = transport_tab_fit_forced(s.C, T, Tsel)*TRANSPORT_CEA_COND_TO_SI;
            return;
        }
        if (Tsel < TRANSPORT_H2O_POWER_BELOW) {
            const double T0 = TRANSPORT_H2O_POWER_BELOW;
            transport_iapws_h2o(T0, &mI, &lI, &nm, &nl);
            const double r = log(T/T0);
            *mu  = mI*exp(nm*r);
            *lam = lI*exp(nl*r);
            return;
        }
        transport_iapws_h2o(T, &mI, &lI, &nm, &nl);
        if (Tsel <= TRANSPORT_H2O_BLEND_LO) { *mu = mI; *lam = lI; return; }
        const double s01 = (T - TRANSPORT_H2O_BLEND_LO)/(TRANSPORT_H2O_BLEND_HI - TRANSPORT_H2O_BLEND_LO);
        const double w   = s01*s01*(3.0 - 2.0*s01);
        const double mC  = transport_tab_fit_forced(s.V, T, Tsel)*TRANSPORT_MICROPOISE_TO_PAS;
        const double lC  = transport_tab_fit_forced(s.C, T, Tsel)*TRANSPORT_CEA_COND_TO_SI;
        *mu  = exp((1.0 - w)*log(mI) + w*log(mC));
        *lam = exp((1.0 - w)*log(lI) + w*log(lC));
        return;
    }
    default:   // CEA / FIT
        *mu  = transport_tab_fit_forced(s.V, T, Tsel)*TRANSPORT_MICROPOISE_TO_PAS;
        *lam = transport_tab_fit_forced(s.C, T, Tsel)*TRANSPORT_CEA_COND_TO_SI;
        return;
    }
}

inline double transport_tab_pair_forced(const TransportPairD& p, const SpeciesTransportD* sp, double T, double Tsel)
{
    if (p.kind == TRANSPORT_PAIR_CE) {
        const double Ts = transport_tab_Tstar_forced(T, Tsel, p.eps_ab);
        const double Mr = 2.0*sp[p.a].MW*sp[p.b].MW/(sp[p.a].MW + sp[p.b].MW)*1000.0;
        return TRANSPORT_CE_MU_CONST*sqrt(Mr*T)/(p.sigma_ab*p.sigma_ab*transport_omega22(Ts));
    }
    return transport_tab_fit_forced(p.V, T, Tsel)*TRANSPORT_MICROPOISE_TO_PAS;   // CEA 相互作用
}

// ---- 境界の列挙 ----
inline void transport_tab_fit_bounds(const TransportFitD& f, std::vector<TransportTabBound>& B)
{
    for (int k = 0; k + 1 < f.n; ++k) B.push_back({f.Thi[k], 1, 1});   // T ≤ Thi[k] が区間 k (最後の区間の Thi は外挿なので境界でない)
}

inline void transport_tab_Tstar_bounds(double eps, std::vector<TransportTabBound>& B)
{
    B.push_back({0.3*eps, 0, 0});     // T* < 0.3 がクランプ (左)
    B.push_back({100.0*eps, 1, 0});   // T* > 100 がクランプ (右)
}

inline std::vector<TransportTabBound> transport_tab_species_bounds(const SpeciesTransportD& s)
{
    std::vector<TransportTabBound> B;
    switch (s.model) {
    case TRANSPORT_MODEL_KINETIC:
        transport_tab_Tstar_bounds(s.eps_kB, B);
        B.push_back({s.thermo.Tlo, 0, 0});    // T < Tlo は c_p を Tlo でクランプ
        B.push_back({s.thermo.Tmid, 0, 1});   // T < Tmid が low の係数
        B.push_back({s.thermo.Thi, 1, 0});    // T > Thi は c_p を Thi でクランプ
        break;
    case TRANSPORT_MODEL_H2O_IAPWS_CEA_V1:
        B.push_back({TRANSPORT_H2O_POWER_BELOW, 0, 0});   // T < 253.15 が冪外挿
        B.push_back({TRANSPORT_H2O_BLEND_LO, 1, 0});      // T ≤ 500 が IAPWS
        B.push_back({TRANSPORT_H2O_BLEND_HI, 0, 0});      // T ≥ 700 が CEA
        transport_tab_fit_bounds(s.V, B);
        transport_tab_fit_bounds(s.C, B);
        break;
    default:
        transport_tab_fit_bounds(s.V, B);
        transport_tab_fit_bounds(s.C, B);
        break;
    }
    return B;
}

inline std::vector<TransportTabBound> transport_tab_pair_bounds(const TransportPairD& p)
{
    std::vector<TransportTabBound> B;
    if (p.kind == TRANSPORT_PAIR_CE) transport_tab_Tstar_bounds(p.eps_ab, B);
    else transport_tab_fit_bounds(p.V, B);
    return B;
}

// 左の区間に属する最大の float (leftIncl: T ≤ b、そうでなければ T < b)
inline float transport_tab_float_upper(double b, int leftIncl)
{
    float f = (float)b;
    if (leftIncl ? ((double)f > b) : ((double)f >= b)) f = std::nextafter(f, -INFINITY);
    return f;
}

// 表 1 つを作る。F(T, Tsel, out[nv]) は Tsel が選ぶ式を T で評価して nv 個の正値を返す。
//   T 範囲 [150, 15000] を境界で分割し、各区間を ln T 等間隔で刻む。
//   失敗 (同じ位置の不連続境界で所属規約が食い違う・非正/非有限の値) は err に書いて false。
template <class F>
inline bool transport_tab_build_one(std::vector<TransportTabBound> B, int nv, F fn, TransportTablesHost& H,
                                    std::vector<TransportHermiteF>& coef, TransportTabRefF& ref,
                                    std::vector<double>& segLo, std::vector<double>& segHi, float& TminReq, float& TmaxReq,
                                    const std::string& what)
{
    const double Tmin = TRANSPORT_TAB_TMIN, Tmax = TRANSPORT_TAB_TMAX;
    std::sort(B.begin(), B.end(), [](const TransportTabBound& x, const TransportTabBound& y) { return x.b < y.b; });
    // 同じ位置の境界をまとめる: 不連続な境界の規約を優先 (連続な境界では左右どちらの式でも値は同じ)
    std::vector<TransportTabBound> U;
    for (const auto& x : B) {
        if (!U.empty() && U.back().b == x.b) {
            TransportTabBound& u = U.back();
            if (x.disc && u.disc && x.leftIncl != u.leftIncl) {
                H.error = what + ": discontinuous boundaries at the same T with different membership";
                return false;
            }
            if (x.disc && !u.disc) { u.leftIncl = x.leftIncl; u.disc = 1; }
            continue;
        }
        U.push_back(x);
    }
    // 範囲の端にある境界: T = 150 が左 (範囲外) の式、T = 15000 が右の式なら、その点は表から外す
    std::vector<double> ends{Tmin};
    std::vector<int> incl;
    for (const auto& x : U) {
        if (x.b == Tmin && x.leftIncl) TminReq = std::max(TminReq, std::nextafter((float)Tmin, INFINITY));
        if (x.b == Tmax && !x.leftIncl) TmaxReq = std::min(TmaxReq, std::nextafter((float)Tmax, -INFINITY));
        if (x.b > Tmin && x.b < Tmax) { ends.push_back(x.b); incl.push_back(x.leftIncl); }
    }
    ends.push_back(Tmax);
    ref.seg0 = static_cast<int>(H.seg.size());
    ref.nseg = static_cast<int>(ends.size()) - 1;
    for (int k = 0; k + 1 < static_cast<int>(ends.size()); ++k) {
        const double Ta = ends[k], Tb = ends[k + 1];
        const double la = log(Ta), lb = log(Tb);
        const double Tsel = exp(0.5*(la + lb));   // この区間の式を選ぶ温度 (区間の内側)
        auto g = [&](double x, int v) {
            double out[2];
            fn(exp(x), Tsel, out);
            return log(out[v]);
        };
        int m = std::max(1, static_cast<int>(std::ceil((lb - la)*TRANSPORT_TAB_PER_LNT - 1.0e-9)));
        std::vector<TransportHermiteF> c;
        for (int refine = 0; ; ++refine) {
            c.assign(static_cast<size_t>(m)*nv, TransportHermiteF{0, 0, 0, 0});
            const double h = (lb - la)/m;
            const double d = 1.0e-3;   // 4 次の中心差分の刻み (ln T): 打切り ~d⁴ g⁽⁵⁾/30、丸め ~1e-16/d
            double pmax = 0.0, emax = 0.0;
            bool ok = true;
            for (int v = 0; v < nv && ok; ++v) {
                for (int j = 0; j < m; ++j) {
                    const double x0 = la + h*j, x1 = (j + 1 == m) ? lb : la + h*(j + 1);
                    const double T0 = (j == 0) ? Ta : exp(x0), T1 = (j + 1 == m) ? Tb : exp(x1);
                    double o0[2], o1[2];
                    fn(T0, Tsel, o0);
                    fn(T1, Tsel, o1);
                    if (!(o0[v] > 0.0) || !(o1[v] > 0.0) || !std::isfinite(o0[v]) || !std::isfinite(o1[v])) { ok = false; break; }
                    const double g0 = log(o0[v]), g1 = log(o1[v]);
                    const double d0 = (-g(x0 + 2*d, v) + 8*g(x0 + d, v) - 8*g(x0 - d, v) + g(x0 - 2*d, v))/(12*d);
                    const double d1 = (-g(x1 + 2*d, v) + 8*g(x1 + d, v) - 8*g(x1 - d, v) + g(x1 - 2*d, v))/(12*d);
                    const double c1 = h*d0;
                    const double c2 = 3.0*(g1 - g0) - h*(2.0*d0 + d1);
                    const double c3 = 2.0*(g0 - g1) + h*(d0 + d1);
                    pmax = std::max(pmax, std::fabs(c1) + std::fabs(c2) + std::fabs(c3));
                    // 小区間の中点で Hermite (double) と式を比べる (刻みの精度の自己点検)
                    const double gm = log([&]{ double o[2]; fn(exp(0.5*(x0 + x1)), Tsel, o); return o[v]; }());
                    const double hm = g0 + 0.5*(c1 + 0.5*(c2 + 0.5*c3));
                    emax = std::max(emax, std::fabs(hm - gm));
                    c[static_cast<size_t>(j)*nv + v] = TransportHermiteF{(float)o0[v], (float)c1, (float)c2, (float)c3};
                }
            }
            if (!ok) {
                H.error = what + ": non-positive or non-finite value in [" + std::to_string(Ta) + ", " + std::to_string(Tb) + "] K";
                return false;
            }
            if ((pmax <= TRANSPORT_TAB_PMAX && emax <= 1.0e-8) || refine >= 6) break;
            m *= 2;   // 精度不足の区間だけ細分する
        }
        TransportSegF s;
        s.lnTa = la;
        s.invH = m/(lb - la);
        s.Tupper = (k + 2 == static_cast<int>(ends.size())) ? (float)Tmax : transport_tab_float_upper(Tb, incl[k]);
        s.off = static_cast<int>(coef.size()/nv);
        s.m = m;
        s.pad = 0;
        H.seg.push_back(s);
        coef.insert(coef.end(), c.begin(), c.end());
        segLo.push_back(Ta);
        segHi.push_back(Tb);
    }
    return true;
}

// 全表を作る。sp/pairs: 解決済みの実種・組、MWs: 輸送種の分子量 [kg/mol] (熱物性と同じ値)、
//   expand: 輸送種 × 実種の展開行列 (double, 行優先)。
inline bool transport_tables_build_host(const std::vector<SpeciesTransportD>& sp, const std::vector<TransportPairD>& pairs,
                                        const std::vector<double>& MWs, const std::vector<double>& expand,
                                        TransportTablesHost& H)
{
    H = TransportTablesHost{};
    const int n = static_cast<int>(sp.size());
    float TminReq = (float)TRANSPORT_TAB_TMIN, TmaxReq = (float)TRANSPORT_TAB_TMAX;
    H.spTab.resize(n);
    H.segLo.resize(n + pairs.size());
    H.segHi.resize(n + pairs.size());
    for (int r = 0; r < n; ++r) {
        const SpeciesTransportD& s = sp[r];
        auto f = [&](double T, double Tsel, double* o) { transport_tab_species_forced(s, T, Tsel, &o[0], &o[1]); };
        if (!transport_tab_build_one(transport_tab_species_bounds(s), 2, f, H, H.spc, H.spTab[r], H.segLo[r], H.segHi[r],
                                     TminReq, TmaxReq, "species " + std::to_string(r)))
            return false;
    }
    H.pairTab.resize(pairs.size());
    H.pm.resize(pairs.size());
    for (size_t ip = 0; ip < pairs.size(); ++ip) {
        const TransportPairD& p = pairs[ip];
        const double Ma = sp[p.a].MW, Mb = sp[p.b].MW;
        TransportPairMixF& q = H.pm[ip];
        q.kab = (float)(2.0*Mb/(Ma + Mb));
        q.kba = (float)(2.0*Ma/(Ma + Mb));
        q.gab = (float)(1.0 + 2.41*(Ma - Mb)*(Ma - 0.142*Mb)/((Ma + Mb)*(Ma + Mb)));
        q.gba = (float)(1.0 + 2.41*(Mb - Ma)*(Mb - 0.142*Ma)/((Ma + Mb)*(Ma + Mb)));
        q.rsC = (float)(5.656854*sqrt(Mb/(Ma + Mb)));
        q.rsRatio = (float)sqrt(Mb/Ma);
        if (p.kind == TRANSPORT_PAIR_RIGID_SPHERE) { H.pairTab[ip] = TransportTabRefF{static_cast<int>(H.seg.size()), 0}; continue; }
        auto f = [&](double T, double Tsel, double* o) { o[0] = transport_tab_pair_forced(p, sp.data(), T, Tsel); o[1] = 1.0; };
        if (!transport_tab_build_one(transport_tab_pair_bounds(p), 1, f, H, H.pairc, H.pairTab[ip], H.segLo[n + ip], H.segHi[n + ip],
                                     TminReq, TmaxReq, "pair " + std::to_string(ip)))
            return false;
    }
    H.invMWs.resize(MWs.size());
    for (size_t s = 0; s < MWs.size(); ++s) H.invMWs[s] = (float)(1.0/MWs[s]);
    H.expandF.resize(expand.size());
    for (size_t k = 0; k < expand.size(); ++k) H.expandF[k] = (float)expand[k];
    H.Tmin = TminReq;
    H.Tmax = TmaxReq;
    return true;
}

// host ポインタのまま見え方を作る (単体試験用)
inline TransportTablesF transport_tables_view_host(const TransportTablesHost& H)
{
    TransportTablesF t;
    t.valid = 1;
    t.Tmin = H.Tmin; t.Tmax = H.Tmax;
    t.seg = H.seg.data(); t.spc = H.spc.data(); t.pairc = H.pairc.data();
    t.spTab = H.spTab.data(); t.pairTab = H.pairTab.data(); t.pm = H.pm.data();
    t.invMWs = H.invMWs.data(); t.expandF = H.expandF.data();
    return t;
}

// host 側の表 (thermo_d.cu が所有; 表を使っていなければ nullptr)。試験ハーネス (FORGE_TRANSPORT_TABLE_PROBE) が配置を書き出す。
const TransportTablesHost* thermo_transport_tables_host();
