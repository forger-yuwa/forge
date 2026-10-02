#pragma once
// =============================================================================
// renormGate_d.cuh — 再正規化の受入ゲートの計測 (plans/active/condensation-two-phase-transport.md §5.1 #1b-pre;
//   codex diagnose notes/reviews/2026-10-02-twophase-1b-final-diagnose.md)。計上だけで、再正規化の算術・書き込み値は変えない。
//
// 更新 n (再正規化の 1 回の呼び出し) ごと・成分 q ∈ {ρY_w, ρg, ρQ2, ρQ1, ρQ0} ごとに
//     C_q,n = Σ_i |q⁺_i − q⁻_i| V_i / Σ_i q⁻_i V_i      (q⁻・q⁺ は再正規化の直前・直後の格納値を double に上げたもの; root のみ)
// と、更新ごとの局所係数偏差 F_n = max_i |f_i − 1| (f = ρ/ΣρY) を作り、更新の履歴に 1 行 (RNG_ENTRY 個の double) ずつ残す。
// 分母 0 の成分は「絶対補正 0」を要求する: 分子も 0 なら C = 0、分子 > 0 なら C = +inf (不合格)。非有限はそのまま (不合格)。
// 判定 (host): 実更新数 N の最後の ceil(0.1 N) 更新の max_n C_q,n と max_n F_n がすべて有限かつ ≤ κ = 2 n_s ε₃₂。
// 窓は履歴から取るのでログの区切り (monitorInterval) に依らない。
// =============================================================================
#include <cmath>
#include <cstddef>
#include <vector>

#define RNG_NC    5                 // ρY_w, ρg, ρQ2, ρQ1, ρQ0
#define RNG_ACC   (2*RNG_NC + 1)    // 1 更新の device 集計: [0..4] 分子 Σ|Δq|V, [5..9] 分母 Σq⁻V, [10] max|f−1|
#define RNG_ENTRY (RNG_NC + 1)      // 履歴 1 行: C_q (5) と F

#if defined(__CUDACC__)
__device__ inline void rng_atomic_max_double(double* a, double v)
{
    unsigned long long* p = reinterpret_cast<unsigned long long*>(a);
    unsigned long long old = *p;
    // NaN は必ず残す (比較が偽になるので専用に書く)
    if (!(v == v)) { atomicExch(p, (unsigned long long)__double_as_longlong(v)); return; }
    while (v > __longlong_as_double((long long)old)) {
        const unsigned long long assumed = old;
        old = atomicCAS(p, assumed, (unsigned long long)__double_as_longlong(v));
        if (old == assumed) break;
    }
}
// 1 セルぶんを 1 更新の集計 acc に足す (root のセルだけで呼ぶ)。qm・qp は直前・直後の格納値 (float を double に上げた値)。
__device__ inline void rng_accumulate(double* acc, const double* qm, const double* qp, double V, double fdev)
{
    for (int k = 0; k < RNG_NC; ++k) {
        const double d = fabs(qp[k] - qm[k]);
        if (d != 0.0 || !(d == d)) atomicAdd(&acc[k], d*V);
        if (qm[k] != 0.0 || !(qm[k] == qm[k])) atomicAdd(&acc[RNG_NC + k], qm[k]*V);
    }
    rng_atomic_max_double(&acc[2*RNG_NC], fdev);
}
// 1 更新の集計を履歴の 1 行にする (1 スレッド)。
__global__ void rng_finalize_d(const double* acc, double* entry)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    for (int k = 0; k < RNG_NC; ++k) {
        const double num = acc[k], den = acc[RNG_NC + k];
        double C;
        if (den == 0.0) C = (num == 0.0) ? 0.0 : INFINITY;   // 総量 0 の成分は絶対補正 0 を要求
        else C = num/den;                                      // 非有限は非有限のまま (判定で不合格)
        entry[k] = C;
    }
    entry[RNG_NC] = acc[2*RNG_NC];
}
#endif

// ---- host: 窓の集計と判定 (ソルバのログと単体試験が同じ関数を使う) ----
struct RngWindow {
    double C[RNG_NC];   // max_n C_q,n (非有限が 1 つでもあれば NaN/inf のまま)
    double F;           // max_n F_n
    long   nonfinite;   // 窓内の非有限 (分母 0 で分子 > 0 の inf を含む) の数
    size_t begin, end;  // 履歴の [begin, end)
};
inline RngWindow rng_window_max(const std::vector<double>& hist, size_t begin, size_t end)
{
    RngWindow w{};
    for (int k = 0; k < RNG_NC; ++k) w.C[k] = 0.0;
    w.F = 0.0; w.nonfinite = 0; w.begin = begin; w.end = end;
    for (size_t n = begin; n < end; ++n) {
        const double* e = &hist[n*RNG_ENTRY];
        for (int k = 0; k <= RNG_NC; ++k) {
            const double v = e[k];
            double& m = (k < RNG_NC) ? w.C[k] : w.F;
            if (!std::isfinite(v)) { ++w.nonfinite; m = v; continue; }   // 非有限は最大値に残す
            if (std::isfinite(m) && v > m) m = v;
        }
    }
    return w;
}
// 実更新数 N の最後の ceil(0.1 N) 更新 (N ≥ 1 なら最低 1)。
inline size_t rng_final_window_begin(size_t N) { const size_t W = (N == 0) ? 0 : (size_t)std::ceil(0.1*(double)N); return N - W; }
inline bool rng_judge(const RngWindow& w, double kappa)
{
    if (w.end <= w.begin || w.nonfinite != 0) return false;
    for (int k = 0; k < RNG_NC; ++k) if (!(w.C[k] <= kappa)) return false;
    return (w.F <= kappa);
}
inline double rng_kappa(int nSpecies) { return 2.0*(double)nSpecies*1.1920928955078125e-7; }   // κ = 2 n_s ε₃₂
