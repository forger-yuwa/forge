// axisSymmetryPin_d.cuh — node 軸対称の軸ノードの状態ピン (1 節点分の処理)
//   (enforceAxisSymmetry_d の本体。単体試験 tests/unit/test_qacc_axis.cu から直接駆動するため __device__ inline に切り出した。
//    plans/active/architecture-float-state-double-geometry.md §5.1 #18)
//
// 軸上 CV (R=0) で半径方向運動量を 0 に (軸対称の対称条件)。
// エネルギー整合のため、除去する半径方向運動エネルギー 0.5*roUy^2/ro を roe からも引く
// (これを怠ると内部エネルギー=圧力が過大になり発散する)。
//
// **作用する 3 つの状態**:
//   作用用の Q (ro, roUy, roe, Uy) : 残差と Jacobian を作る状態。commit の基準ではない (commit は Q_N か Qacc から作り直す)。
//   commit の基準 Q_N (roUyN, roeN) : OFF の commit (ro = roN + dq) が読む。射影しないと commit が stale な roUyN を巻き戻す。
//   FP64 の正本 Qacc (qacc_roUy, qacc_roe) : qAccumulatorFP64 の commit (Qacc += dq) が読む基準。nullptr = OFF。
// **Qacc には Q_N と同じ射影を当てる** (Q_N を射影する前の roUyN から作った運動エネルギー ke_N を Qacc_roe から引き、
// Qacc_roUy を 0 にする)。作用用の Q から作らないこと: 壁と軸に重なる節点では壁 no-slip が先に作用用の roUy を 0 にして
// いるので、作用用の Q から作ると Qacc の半径方向運動エネルギーが残る (codex diagnose 2026-10-11 Major 1 の反例)。
// 作用用の roe を Qacc から上書きしないこと (壁 no-slip が除いた運動エネルギーを戻してしまう)。ミラーへの書き戻しは commit が担う。
//
// **OFF の演算列は切り出し前と同一** (式の行は切り出し前の文字のまま。軸の flag の判定を関数の頭に移しただけ。
// float のビルドで OFF の浮動小数点の命令列が HEAD e24650a8 の enforceAxisSymmetry_d と同じことを SASS で確かめた)。
#ifndef AXIS_SYMMETRY_PIN_D_CUH
#define AXIS_SYMMETRY_PIN_D_CUH

#include "flowFormat.hpp"

__device__ inline void axisSymmetryPinCell
(
    geom_int ic, const geom_int* axis_flag,
    flow_float* ro, flow_float* roUy, flow_float* roe, flow_float* Uy,
    // 擬似時間ステップ開始値も同時にピン (implicit commit ro=roN+dq が stale roUyN を巻き戻さないため)。nullptr 可。
    flow_float* roUyN, flow_float* roeN, flow_float* roN,
    // FP64 の正本 (qAccumulatorFP64)。nullptr = OFF。確保は内点 (ic < nCells) だけ。
    double* qacc_roUy, double* qacc_roe
)
{
    if (axis_flag[ic] != 1) return;
    const flow_float r = ro[ic];
    if (r > (flow_float)0.0) {
        roe[ic] -= (flow_float)0.5 * roUy[ic] * roUy[ic] / r; // 半径方向 KE を除去
    }
    roUy[ic] = (flow_float)0.0;
    Uy[ic]   = (flow_float)0.0;
    if (roUyN != nullptr) {
        const flow_float rN = roN[ic];
        if (rN > (flow_float)0.0) {
            roeN[ic] -= (flow_float)0.5 * roUyN[ic] * roUyN[ic] / rN;
            // ON: 上の行と同じ式・同じガードの ke_N を、roUyN を 0 にする前に作って Qacc に当てる (値は上の行が引いた量と同じ)。
            // 上の行と変数を共有しないのは、共有すると OFF の機械語で負号の置き場所が変わるため (値は同じだが、OFF の
            // SASS を HEAD と同じ FMUL −0.5 / FADD の列に保つ)。
            if (qacc_roe != nullptr) qacc_roe[ic] -= (double)((flow_float)0.5 * roUyN[ic] * roUyN[ic] / rN);
        }
        roUyN[ic] = (flow_float)0.0;
        if (qacc_roUy != nullptr) qacc_roUy[ic] = 0.0;
    }
}

#endif  // AXIS_SYMMETRY_PIN_D_CUH
