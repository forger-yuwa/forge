#pragma once
#include "flowFormat.hpp"
#include <cmath>

// 理由別の補正量監視 (plans/active/condensation-two-phase-transport.md §4.3, #2)。凝縮種ごとに COND_REASON_N 個の double
// (device; condCorrReasons(s) が所有)。和・数は全期間の累積 (区間値は host が前回ログとの差で出す)、最小・最大は区間値 (ログごとに戻す)。
// 計上だけで、クランプの算術と書き込み値は変えない。root (node 周期) のみ・体積 V を掛けた量 [kg (軸対称は rad 当たり)]。
#define COND_REASON_CAP_SUM   0   // 蒸気上限違反 g > Y_w (0.99ρ / 定数 Y_w): Σ (ρg − 上限) V
#define COND_REASON_CAP_N     1   //   同・作動ノード数
#define COND_REASON_NEG_SUM   2   // 負値 floor ρg < 0 → 0: Σ |ρg| V
#define COND_REASON_NEG_N     3
#define COND_REASON_RM_SUM    4   // 液滴消滅 (物理; S ≤ 1 の小液滴・塵): Σ ρg V
#define COND_REASON_RM_N      5
#define COND_REASON_PQ1_SUM   6   // モーメント射影 (Q1, Q2 だけ動かす; 塵の Q1=Q2=0 化を含む): Σ |ΔρQ1| V
#define COND_REASON_PQ2_SUM   7   //   Σ |ΔρQ2| V
#define COND_REASON_P_N       8
#define COND_REASON_VMIN      9   // 制限前の最小蒸気分率 (ρY_w − ρg)/ρ (carrier のみ; 区間値, 初期 1e300)
#define COND_REASON_RN_SUM   10   // 化学種再正規化が凝縮種 ρY_w に掛けた補正: Σ |Δ(ρY_w)| V (species_renormalize_d)
#define COND_REASON_RN_MAX   11   //   max |係数 − 1| (区間値)
#define COND_REASON_RN_N     12   //   ρY_w が変わったノード数
#define COND_REASON_N        16

#if defined(__CUDACC__)
__device__ inline void cond_atomic_min_double(double* a, double v)
{
    unsigned long long* p = reinterpret_cast<unsigned long long*>(a);
    unsigned long long old = *p;
    while (v < __longlong_as_double((long long)old)) {
        const unsigned long long assumed = old;
        old = atomicCAS(p, assumed, (unsigned long long)__double_as_longlong(v));
        if (old == assumed) break;
    }
}
__device__ inline void cond_atomic_max_double(double* a, double v)
{
    unsigned long long* p = reinterpret_cast<unsigned long long*>(a);
    unsigned long long old = *p;
    while (v > __longlong_as_double((long long)old)) {
        const unsigned long long assumed = old;
        old = atomicCAS(p, assumed, (unsigned long long)__double_as_longlong(v));
        if (old == assumed) break;
    }
}
// クランプ前の上限違反・負値・制限前の蒸気分率を計上 (cond_realizability_clamp_{,f_}d の入口; 値は変えない)。
//   yv_pre: 制限前の蒸気分率 (carrier 以外は計上しない: 呼び出し側で have_yv = false)。
__device__ inline void cond_reason_pre(double* rs, double V, double r_in, double gmax, bool have_yv, double yv_pre)
{
    if (rs == nullptr) return;
    if (have_yv) cond_atomic_min_double(&rs[COND_REASON_VMIN], yv_pre);
    if (r_in < 0.0) { atomicAdd(&rs[COND_REASON_NEG_SUM], -r_in*V); atomicAdd(&rs[COND_REASON_NEG_N], 1.0); }
    const double rc = (r_in > 0.0) ? r_in : 0.0;
    if (rc > gmax) { atomicAdd(&rs[COND_REASON_CAP_SUM], (rc - gmax)*V); atomicAdd(&rs[COND_REASON_CAP_N], 1.0); }
}
__device__ inline void cond_reason_proj(double* rs, double V, double q1_in, double q2_in, double q1_out, double q2_out)
{
    if (rs == nullptr || (q1_in == q1_out && q2_in == q2_out)) return;
    atomicAdd(&rs[COND_REASON_PQ1_SUM], fabs(q1_out - q1_in)*V);
    atomicAdd(&rs[COND_REASON_PQ2_SUM], fabs(q2_out - q2_in)*V);
    atomicAdd(&rs[COND_REASON_P_N], 1.0);
}
#endif
