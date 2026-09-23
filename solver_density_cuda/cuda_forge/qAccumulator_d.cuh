// qAccumulator_d.cuh — 保存量の FP64 影アキュムレータ (plans/active/time_integration-fp64-accumulator.md §4.3)
//
// **問題**: 定常陰解法の commit は `Q = Q_N + dq` で、Q は float32。定常解へ近づくほど dq は小さくなり、
// |dq| < 1/2 ULP(Q) になった時点で加算が丸めで消え、残差が残っているのに場が動かなくなる。
// 実測 (case/56.gap_tp1187 の深いすきま): <|dq|> = 0.159 ULP、1 step で値が動く CV は 0.03 %、
// 実効更新は意図の 0.3 %。この状態では |ṁ| が step のべき乗則でしか減らない (step 倍加で 15 %)。
//
// **方針**: Q / Q_N は FP32 のまま**演算に使い**、内点 CV だけの FP64 配列 Qacc を**正本**として加える。
//   commit    : Qacc += dq (FP64・in-place) → Q = (float)Qacc
//   reconcile : Q != (float)Qacc のセルだけ Qacc = (double)Q
// reconcile は「**FP32 の writer が書き換えたセルはその値を採用し、それ以外は残余を保持する**」規則。
// 壁ピン・周期ミラー・化学種の roe += droe など commit 以外の writer と自然に両立する。
//
// **ビット同一性**: f32(a+b) == f32(f64(a)+f64(b)) が最近接丸めで成り立つ (20 万サンプルで実測 100 %)。
// したがって**残余ゼロなら ON は OFF とビット同一**で、ON 1 step と OFF 1 step の比較が安価な補助ゲートになる。
//
// Q の型は変えない (c_d["ro"] は 112 箇所から読まれる)。OFF 経路はカーネル 1 本も変わらない。
#ifndef Q_ACCUMULATOR_D_CUH
#define Q_ACCUMULATOR_D_CUH

#include "flowFormat.hpp"

// commit: dq を FP64 の正本へ積み、FP32 ミラーへ書き戻す。
// scale は updateGuardScale の縮小率 (呼び出し側で評価済み。既定 1)。
__device__ inline void qaccCommitCell(double& qacc, flow_float& q, flow_float dq, flow_float scale)
{
    qacc += (double)dq * (double)scale;
    q = (flow_float)qacc;
}

// reconcile: FP32 側が別の writer に書き換えられていたらそれを採用する。
// 戻り値は「採用した (= 他の writer が触った)」かどうか。呼び出し側でカウンタに積む。
__device__ inline bool qaccReconcileCell(double& qacc, flow_float q)
{
    if (q != (flow_float)qacc) { qacc = (double)q; return true; }
    return false;
}

#endif  // Q_ACCUMULATOR_D_CUH
