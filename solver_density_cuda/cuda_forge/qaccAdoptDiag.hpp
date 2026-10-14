#pragma once

// 診断 FORGE_DIAG_QACC_ADOPT=<接頭辞> (plans/active/architecture-float-state-double-geometry.md §5.1 #18・§6.34、既定 off)。
//
// qAccumulatorFP64 の reconcile (updateVariablesOuter_d の冒頭、update_d.cu) が FP32 側の値を「採用」する節点を、指定の step
// について節点ごとに h5 に書く。採用の条件は reconcile と同じ Q != (flow_float)Qacc (qaccReconcileCell)。reconcile の計器
// (monitorInterval ごとの採用数) は位置も量も持たないので、どの writer が触ったかを切り分けるために位置と量を残す。
//   FORGE_DIAG_QACC_ADOPT_STEPS=<step,...> (任意、既定 1): 対象の step (カンマ区切り)。step = iStep + 1 (res_<step>.h5 と同じ
//   数え方)。出力は step ごとに <接頭辞>_<step>.h5。
//
// 時点: advanceImplicitSteady の updateVariablesOuter (reconcile) の直前。commit・SST・化学種・commit 後の等温壁のピンの後。
// したがって mask はその step の reconcile が採用する節点そのもの (採用の計器を数える時点と同じ)。
//
// 記録する量 (h5):
//   /<量>/mask (uint8 [nCells])  : 1 = Q != (flow_float)Qacc (reconcile が採用する)。量は ro, roUx, roUy, roUz, roe。
//   /<量>/Q    (double [nCells]) : そのときの FP32 の Q (flow_float の値を広げただけ)。
//   /<量>/Qacc (double [nCells]) : そのときの FP64 の正本 (reconcile の前の値)。
//   /cells/cc64 (値の位置の double、あれば)、/cells/{axis,wall,iso_wall,bnode}_flag (あれば)。
// 属性: step・iStep・量ごとの採用数と、そのうち軸・壁・等温壁・内点 (bnode_flag 0)・軸と壁の両方の節点の数、
// 格子の軸・壁・等温壁・境界の節点の数と軸と壁の flag が重なる節点の数。
//
// 対象外: qAccumulatorFP64 0 (「何もしない」と出す)。フックは定常の陰解法 (advanceImplicitSteady) にだけ置く
// (qAccumulatorFP64 自体がその経路しか通さない)。
//
// 数値を変えないこと: off (環境変数なし・qAccumulatorFP64 0) ではフックは分岐 1 つ。on でも読むだけ (既定のストリームの同期の
// cudaMemcpy D2H。デバイスの配列は書かず、カーネルの起動の順も変えない。採用の計器 qaccAdopt_d も触らない)。

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

namespace qaccAdoptDiag {

namespace detail {
extern bool g_enabled;   // 環境変数で有効にし、qAccumulatorFP64 1 だった run
}

// 有効か (advanceImplicitSteady の updateVariablesOuter の直前で見る分岐)。
inline bool enabled() { return detail::g_enabled; }

// 起動時に 1 回 (main、qAccumulatorFP64 の確保と拒否の検査の後、時間更新の前)。
void init(solverConfig& cfg, mesh& msh, variables& var);

// reconcile (updateVariablesOuter) の直前。iStep は 0 始まり。対象の step なら Q と Qacc を読み、mask を作って h5 を書く。
void beforeReconcile(solverConfig& cfg, mesh& msh, variables& var, int iStep);

// 終了時。書かなかった対象の step があればその旨を出す。
void finalize();

}  // namespace qaccAdoptDiag
