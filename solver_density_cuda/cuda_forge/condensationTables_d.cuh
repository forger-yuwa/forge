#pragma once
// 凝縮物性の区分 3 次 Hermite 表 (float 評価)。plans/active/condensation-float-speedup.md §4.2-1、methods/condensation.md 実装 §9。
//
// 物性式 (Jacobsen 飽和圧・NASA-9 差の潜熱など) は float で直接評価すると相殺で 1e-4〜1e-2 ずれる (codex 2026-09-13 M1) ため、
// 現行の double 関数 (クランプ・低温外挿込み) からホストで double で表を作り、区間ごとの 4 係数 (float4) を kernel が引く。
//   f(T) ≈ c0 + c1 u + c2 u² + c3 u³,  u = (T − T0 − i h)/h ∈ [0,1],  i = 区間番号
// 係数は両端の値と**片側**微分 (接続点 N2 45/50/70 K, H2O 273.15 K を格子点に置く) の 3 次 Hermite。範囲外は端でクランプ。
// 対象: ln p_sat, L, σ, ρ_l, k_gas, μ_gas (平均自由行程用; n2_mu_gas 固定、現行と同じ)。
#include <cuda_runtime.h>
#include <vector_types.h>
#include "condensationProperties_d.cuh"

struct CondTableF {
    const float4* c = nullptr;   // n 区間 × (c0,c1,c2,c3)
    float T0 = 0.0f;
    float invH = 0.0f;
    int   n = 0;
};

struct CondTablesF {
    CondTableF lnpsat, latent, sigma, rhol, kgas, mugas;
    float Tmin = 0.0f, Tmax = 0.0f;   // 表の範囲 [T0, T0 + n h]
    int   valid = 0;                  // 0: 未構築 (float 経路は使えない)
};

// 凝縮種の定数 (float 写し)。kernel 先頭で CondSpeciesProps から 1 回作る。
struct CondSpeciesPropsF {
    int   model;
    float R, cv, cp, M, Tc, sigmaScale;
    float lnM;   // ln(M/N_A) [kg] を double で用意 (対数空間 CNT 用; float で m³ を作らない)
};
__host__ __device__ inline CondSpeciesPropsF condProps_to_f(const CondSpeciesProps& s)
{
    CondSpeciesPropsF f;
    f.model = s.model; f.R = (float)s.R; f.cv = (float)s.cv; f.cp = (float)s.cp; f.M = (float)s.M; f.Tc = (float)s.Tc;
    f.sigmaScale = (float)s.sigmaScale;
    f.lnM = (float)log(s.M / 6.02214076e23);
    return f;
}

// 表の評価 (値と、任意で dT 微分)。T は表範囲にクランプ。
__host__ __device__ inline float cond_tab_eval(const CondTableF& t, float T, float* dfdT = nullptr)
{
    float x = (T - t.T0) * t.invH;
    if (!(x > 0.0f)) x = 0.0f;                     // NaN も 0 側へ
    int i = (int)x;
    if (i >= t.n) { i = t.n - 1; x = (float)t.n; }
    float u = x - (float)i;
    if (u > 1.0f) u = 1.0f;
    const float4 c = t.c[i];
    if (dfdT) *dfdT = (c.y + u*(2.0f*c.z + 3.0f*u*c.w)) * t.invH;
    return c.x + u*(c.y + u*(c.z + u*c.w));
}
__host__ __device__ inline float cond_tab_lnpsat_f(const CondTablesF& tb, float T, float* d = nullptr) { return cond_tab_eval(tb.lnpsat, T, d); }
__host__ __device__ inline float cond_tab_psat_f  (const CondTablesF& tb, float T) { return expf(cond_tab_eval(tb.lnpsat, T)); }
__host__ __device__ inline float cond_tab_latent_f(const CondTablesF& tb, float T, float* d = nullptr) { return cond_tab_eval(tb.latent, T, d); }
__host__ __device__ inline float cond_tab_sigma_f (const CondTablesF& tb, float T) { return cond_tab_eval(tb.sigma, T); }
__host__ __device__ inline float cond_tab_rhol_f  (const CondTablesF& tb, float T) { return cond_tab_eval(tb.rhol, T); }
__host__ __device__ inline float cond_tab_kgas_f  (const CondTablesF& tb, float T) { return cond_tab_eval(tb.kgas, T); }
__host__ __device__ inline float cond_tab_mugas_f (const CondTablesF& tb, float T) { return cond_tab_eval(tb.mugas, T); }

// ---------------------------------------------------------------------------------------------
// ホスト側の構築 (double)。kernel を含む .cu と単体試験 (host) の両方から使う。
// ---------------------------------------------------------------------------------------------
// (host 専用の構築関数; device パスでも宣言は見える必要があるので #ifndef __CUDA_ARCH__ では囲まない)
#include <vector>
#include <cmath>

struct CondTablesHost {
    std::vector<float4> lnpsat, latent, sigma, rhol, kgas, mugas;
    double T0 = 0.0, h = 0.0; int n = 0;
};

// 格子の規約 (plan §4.2-1): N2 は T0=20, h=0.1, 上端 125.6 (接続点 45/50/70 が格子点)。H2O は T0=120.15, h=0.25, 上端 1200.15 (273.15 が格子点)。
inline void cond_tables_grid(int model, double* T0, double* h, int* n)
{
    if (model == COND_MODEL_H2O) { *T0 = 120.15; *h = 0.25; *n = 4320; }
    else                         { *T0 = 20.0;   *h = 0.1;  *n = 1056; }
}

template <class F>
inline void cond_tables_fill(std::vector<float4>& out, double T0, double h, int n, F f)
{
    out.resize(n);
    const double d = 1.0e-4;   // 片側差分の刻み [K] (double: 打切り O(d f'') ~1e-8 相対, 丸め ~1e-12)
    for (int i = 0; i < n; ++i) {
        const double Ta = T0 + h*i, Tb = Ta + h;
        const double fa = f(Ta), fb = f(Tb);
        const double da = (f(Ta + d) - fa)/d;          // 右側微分 (接続点で左区間へ漏れない)
        const double db = (fb - f(Tb - d))/d;          // 左側微分
        const double c0 = fa;
        const double c1 = h*da;
        const double c2 = 3.0*(fb - fa) - h*(2.0*da + db);
        const double c3 = 2.0*(fa - fb) + h*(da + db);
        out[i] = make_float4((float)c0, (float)c1, (float)c2, (float)c3);
    }
}

inline void cond_tables_build_host(const CondSpeciesProps& s, CondTablesHost& ht)
{
    cond_tables_grid(s.model, &ht.T0, &ht.h, &ht.n);
    cond_tables_fill(ht.lnpsat, ht.T0, ht.h, ht.n, [&](double T){ const double p = cond_psat(s, T); return log(p > 1.0e-300 ? p : 1.0e-300); });
    cond_tables_fill(ht.latent, ht.T0, ht.h, ht.n, [&](double T){ return cond_latent(s, T); });
    cond_tables_fill(ht.sigma,  ht.T0, ht.h, ht.n, [&](double T){ return cond_sigma(s, T); });
    cond_tables_fill(ht.rhol,   ht.T0, ht.h, ht.n, [&](double T){ return cond_rho_cond(s, T); });
    cond_tables_fill(ht.kgas,   ht.T0, ht.h, ht.n, [&](double T){ return cond_kgas(s, T); });
    cond_tables_fill(ht.mugas,  ht.T0, ht.h, ht.n, [&](double T){ return n2_mu_gas(T); });
}

// host ポインタのまま CondTablesF を組む (単体試験用)。
inline CondTablesF cond_tables_view_host(const CondTablesHost& ht)
{
    CondTablesF tb;
    auto mk = [&](const std::vector<float4>& v){ CondTableF t; t.c = v.data(); t.T0 = (float)ht.T0; t.invH = (float)(1.0/ht.h); t.n = ht.n; return t; };
    tb.lnpsat = mk(ht.lnpsat); tb.latent = mk(ht.latent); tb.sigma = mk(ht.sigma); tb.rhol = mk(ht.rhol); tb.kgas = mk(ht.kgas); tb.mugas = mk(ht.mugas);
    tb.Tmin = (float)ht.T0; tb.Tmax = (float)(ht.T0 + ht.h*ht.n); tb.valid = 1;
    return tb;
}

// device へ複製した CondTablesF を返す (呼び出し側が寿命管理; 通常は起動時に 1 回)。
inline CondTablesF cond_tables_upload(const CondTablesHost& ht)
{
    CondTablesF tb = cond_tables_view_host(ht);
    auto up = [&](CondTableF& t, const std::vector<float4>& v){
        float4* d = nullptr;
        cudaMalloc((void**)&d, v.size()*sizeof(float4));
        cudaMemcpy(d, v.data(), v.size()*sizeof(float4), cudaMemcpyHostToDevice);
        t.c = d;
    };
    up(tb.lnpsat, ht.lnpsat); up(tb.latent, ht.latent); up(tb.sigma, ht.sigma); up(tb.rhol, ht.rhol); up(tb.kgas, ht.kgas); up(tb.mugas, ht.mugas);
    return tb;
}

