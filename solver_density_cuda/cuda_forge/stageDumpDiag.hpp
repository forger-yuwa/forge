#pragma once

// 診断 FORGE_DIAG_STAGE_DUMP=<out.h5> (plans/active/architecture-float-state-double-geometry.md §5.1 #17・§6.30、既定 off)。
//
// 定常の陰解法の外側の step (擬似時間の 1 回の非線形更新) 1 つについて、平均流の保存量 5 つ (ro, roUx, roUy, roUz, roe) の
// 節点ごと (自前の CV、ic < nCells) の値を、更新の段ごとに double で <out.h5> に書く。書いた後も計算は続ける。
//   FORGE_DIAG_STAGE_STEP=<n> (任意、既定 1): 対象の step。step = iStep + 1 (res_<step>.h5・commit_loss.csv と同じ数え方)。
//
// 記録する量 (h5 の /<量>/{ro,roUx,roUy,roUz,roe}、どれも [nCells] の double。flow_float の値を広げただけ):
//   Q0      : step の頭 (advanceImplicitSteady の最初、組立の前) の ro..roe。前の step の末尾 (または初期化) の状態。
//   Q_asm   : 組立の後・線形 solve の直前の ro..roe (残差と Jacobian を作った状態)。組立の前処理が作業配列として上書きした値
//             (壁 no-slip・軸のピン・等温壁のピン・EOS の床と TP の roe の組み直し・境界条件) を含む。commit の基準ではない。
//   R       : 線形 solve の直前の res_ro..res_roe (組立と、あれば speciesImplicitCoupling 2 の res_roe への移項の後)。
//             block DPLUR / line Thomas が右辺として読む配列そのもの (軸・壁・等温壁の行の 0 化はカーネル内の写しに当たり、
//             この配列は変えない)。
//   b       : commit のカーネルが足し込む基準 roN..roeN (update_d.cu の ro = roN + d0 の roN)。commit の直後に読む
//             (カーネルは roN を書かない)。roN は前の step の updateVariablesOuter で ro を写したものだが、組立の前処理の
//             軸のピン (enforceAxisSymmetry_d) が軸の節点の roUyN・roeN を射影するので、Q0 と一致するとは限らない。
//   d       : commit のカーネルが b に足す更新 (block DPLUR は dq_block_old_0..4、スカラー DPLUR は dq_*_old。最後の sweep・
//             swap・周期のミラーの後で、カーネルが読む配列そのもの)。implicitRelax は solve のカーネル (点の 5×5) と
//             line Thomas の後退代入で掛け済み。EOS の床は commit に入らない (定常の陰解法の commit は ro = roN + d で、床は
//             組立の dependentVariables が作業配列 ro・roe に当てるだけ。その分は Q_asm − b に現れる)。
//             updateGuardAlpha > 0 (縮小率 s を写していない) と qAccumulatorFP64 (基準が FP64 の正本) は未対応 (下)。
//   q       : commit のカーネルの直後の ro..roe (SST・化学種・commit 後の等温壁のピン・updateVariablesOuter の前)。
//   q_final : step の末尾の ro..roe (commit 後の等温壁のピン、凍結しなければ SST の sstEnergyIncludesK の roe 補正、
//             updateVariablesOuter と出力の後) = 次の step の Q0。次の step の組立の前処理の上書き (EOS の床・境界条件) は
//             含まない (それは次の step の Q_asm に入り、commit には入らない)。
// 照合 (属性と標準出力): commit の式 q == fl(b + d) が全節点・5 量でビット一致するか (b と d が commit の実際の入力で
// あることの検査。float の和は double で足して float に丸めても同じ値になる)、b ≠ Q0・q_final ≠ q の値の数。
// メタ: /cells/cc64 (値の位置の double、あれば)、/cells/{bnode,axis,wall,iso_wall}_flag・line_prev・line_next (あれば)。
// 属性: step・iStep・flow_float_bytes・geom_float_bytes・path (lineImplicit = block DPLUR + line Thomas、blockDPLUR = ライン
// なしの block DPLUR (コード中の「point 経路」)、scalarDPLUR)・commit_kernel・d_source・設定 (implicitRelax ほか)・
// FORGE_FREEZE_TURB / FORGE_FREEZE_SPECIES・各量の時点の説明。
//
// commitLossDiag (FORGE_DIAG_COMMIT_LOSS) との関係: b・d・q はそれの Q_before・dq_req・Q_after と同じ配列を同じ時点
// (commit のラッパの中、カーネルの直後) で読み、対象外の経路 (qAccumulatorFP64・updateGuardAlpha > 0) も同じにした。
// 違いは、集計でなく節点ごとの値を書くこと、Q0・Q_asm・R・q_final を足すこと、指定の 1 step だけであること。
//
// 対象外 (「未対応」と標準出力に出して何も書かない。計算は止めない): gpu 0、陽解法、dual-time (unsteady 1)、
// qAccumulatorFP64 1、blockDPLUR 1 で updateGuardAlpha > 0。起動時に判定できない経路 (commit のフックを通らない、
// 1 step に commit が 2 回ある等) も step の末尾で「未対応」とする。
//
// 数値を変えないこと: off (環境変数なし) では各フックは分岐 1 つで、確保もカーネルの起動もしない。on でも読むだけ
// (既定のストリームの同期の cudaMemcpy D2H。デバイスの配列は書かず、カーネルの起動の順も変えない)。

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

namespace stageDumpDiag {

namespace detail {
extern bool g_enabled;   // 環境変数で有効にし、対象の経路だった run
extern bool g_active;    // 対象の step の中 (beginStep 〜 endStep)
}

// 有効か (main の step の頭で見る分岐)。
inline bool enabled() { return detail::g_enabled; }
// 対象の step の中か (solve の直前・commit のラッパが見る分岐)。
inline bool active() { return detail::g_active; }

// 起動時に 1 回 (main、時間更新の前)。環境変数を読み、対象の経路なら有効にする。対象外なら「未対応」と出して無効のまま。
void init(solverConfig& cfg, mesh& msh);

// 定常の外側の step の頭 (advanceImplicitSteady の最初)。iStep は 0 始まり。対象の step なら Q0 を読む。
void beginStep(solverConfig& cfg, mesh& msh, variables& var, int iStep);

// 線形 solve (blockDPLURSolve) の直前 (implicitNonlinearUpdate)。Q_asm と R を読む。
void beforeSolve(solverConfig& cfg, mesh& msh, variables& var);

// commit のカーネルの直後 (update_d.cu の commit のラッパ)。b = roN..roeN、d = dq[0..4]、q = ro..roe を読む。
// kernel・dqNames は属性に書く名前。guardAlpha は block の updateGuardAlpha (スカラー DPLUR は 0)。
void afterCommit(solverConfig& cfg, mesh& msh, variables& var, const char* kernel,
                 const char* const dqNames[5], flow_float guardAlpha);

// 定常の外側の step の末尾 (advanceImplicitSteady の最後)。q_final を読み、そろっていれば h5 を書いて 1 行出す。
void endStep(solverConfig& cfg, mesh& msh, variables& var);

// 終了時。対象の step に達しなかったらその旨を出す。
void finalize();

}  // namespace stageDumpDiag
