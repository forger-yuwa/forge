#pragma once

// 診断 FORGE_DIAG_QUANTIZE_STATE=1 (plans/active/architecture-float-state-double-geometry.md §5.1 #19・#20・§6.36・§6.40、
// 既定 off)。
//
// FP64 のビルドで、定常の陰解法の各 step の末尾に保存量を float32 に丸めて double に戻す (十分性の試験: step の末尾の
// 量子化だけで float のビルドと同じ偏りが出るか)。step の中の演算・途中の値は double のまま。
//   対象: ro, roUx, roUy, roUz, roe、SST (main の scalarResidualEnabled) なら roK, roOmega、化学種 (nSpecies ≥ 2) なら
//   roY0..roY{n-1} すべて。全 CV (msh.nCells_all、ghost を含む)。ほかの配列 (履歴・原始量・補正) には触らない。
//   丸めない保存量 (遷移モデル・受動トレーサ・凝縮モーメント) があれば起動時に「注意」と出す (拒否はしない)。
//   時点: advanceImplicitSteady の commit 後の等温壁のピン (applyNodeIsothermalWallPin) の後・updateVariablesOuter の前。
//   したがって次の step の基準 (roN..roeN・roKN・roOmegaN と M) は丸めた値から写される。化学種の基準 roY{s}N・M は
//   次の step の化学種の更新の頭 (speciesUpdateOuter_d_wrapper) で roY から写される。
//
// 丸める量を選ぶ FORGE_DIAG_QUANTIZE_STATE_ONLY=<量,...> (任意、§5.1 #20。例 roe): 設定されていれば上の対象のうち名前を
// 挙げた量だけを丸める (並びは上の対象の順)。対象に無い名前 (この構成に無い量を含む)・空の値・空の要素・重複は起動時に
// stderr に理由を出して exit(1) (名前の前後の空白は取り除かない)。未設定なら #19 のまま (全対象を丸め、ログも同じ)。
// 起動時の行の「丸める量」の欄は「選んだ量だけ (roe)、FORGE_DIAG_QUANTIZE_STATE_ONLY。丸めない対象: ro, roUx, ...」になる。
//
// 拒否 (起動時に stderr に理由を出して exit(1)): float のビルド (sizeof(flow_float) != 8)、定常の陰解法以外
// (timeIntegration != 11 または unsteady != 0)、qAccumulatorFP64 1 との併用、cell (node のみ)、gpu 0 (カーネルで丸める)。
// 環境変数の値は "1" (有効) か "0" (無効) だけを受け付け、ほかの値は止める。
//
// 自己点検 FORGE_DIAG_QUANTIZE_STATE_CHECK=<step,...> (任意、既定 1。step = iStep + 1、res_<step>.h5 と同じ数え方):
//   (a) 対象の step の updateVariablesOuter の直後に、Q (対象の保存量)・Q_N (roN..roeN、SST なら roKN・roOmegaN)・
//       Q_M (roM..roeM、SST なら roKM・roOmegaM) を host に写し、float32 で表せない値 ((double)(float)x != x、NaN を含む)
//       の個数を量ごとに数えて 1 行出す:
//         [quantizeState] step N: 表せない値 ro 0 roUx 0 ... (Q) / roN 0 ... (Q_N) / roM 0 ... (Q_M) / 合計 0 (...)
//   (b) 化学種があれば、対象の step の次の step の化学種の更新の頭 (speciesUpdateOuter_d_wrapper の直後) に roY{s}N・
//       roY{s}M を数えて 1 行出す (step 番号は対象の step N のまま):
//         [quantizeState] step N: 表せない値 roY0N 0 ... (roY_N) / roY0M 0 ... (roY_M) / 合計 0 (...)
//   どれも 0 を期待する。0 でなければ、内点 (ic < nCells) に限った個数と、そこから node の入口のピン
//   (scalarDirichletPin 1) の節点を除いた個数を「注意」の行で続けて出す (計算は止めない)。(b) の roY{s}N は次の step の
//   組立の後に取られるので、組立の中で ρY を書く writer の値が入る: node の入口の化学種 Dirichlet
//   (speciesTransport_d.cu の species_dirichlet_boundary_d、内点の入口ノードと ghost に roY = ro·Y_in) と ghost の境界値。
//   これらは丸めた後の step の中の演算なので、(b) はそれらの節点で 0 にならないことがある。
//   FORGE_DIAG_QUANTIZE_STATE_ONLY で量を選んだとき:
//   - (a)・(b) は選んだ量とその履歴だけを数える (roe なら Q = roe、Q_N = roeN、Q_M = roeM。化学種 roY{s} を選べば
//     (a) の Q と (b) の roY{s}N・roY{s}M)。選ばなかった量の表せない値は数えない (0 を要求しない)。化学種を選ばなければ
//     (b) は出さない。空の群は行に出さない。
//   - (c) 対象の step の丸めのカーネルの直前と直後に、選ばなかった対象 (上の対象から選んだ量を除いたもの、全 CV) を
//     host に写し、生のビットが変わった値の個数を量ごとに 1 行出す (0 が期待。カーネルは選んだ配列にしか書かない):
//         [quantizeState] step N: 選ばなかった量の丸め前後の変化 ro 0 roUx 0 ... / 合計 0 (0 が期待。全 CV ...、...)
//
// 数値を変えないこと: off (環境変数なし) では各フックは分岐 1 つで、確保もカーネルの起動もしない。

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

struct cudaConfig;

namespace quantizeStateDiag {

namespace detail {
extern bool g_enabled;          // 環境変数で有効にし、拒否に当たらなかった run
extern bool g_speciesPending;   // (b) の点検を待っている step がある
}

// 有効か (advanceImplicitSteady の step の末尾・updateVariablesOuter の直後で見る分岐)。
inline bool enabled() { return detail::g_enabled; }
// (b) の点検を待っているか (化学種の更新の頭で見る分岐)。
inline bool speciesPending() { return detail::g_speciesPending; }

// 起動時に 1 回 (main、qAccumulatorFP64 の検査の後・時間更新の前)。sst は main の scalarResidualEnabled(cfg)
// (roK・roOmega を更新する構成か)。環境変数を読み、拒否に当たれば止める。
void init(solverConfig& cfg, mesh& msh, variables& var, bool sst);

// step の末尾 (advanceImplicitSteady の applyNodeIsothermalWallPin の後・updateVariablesOuter の前)。対象の保存量を
// 全 CV で float32 に丸める。iStep は 0 始まり。量を選んだときの自己点検の step なら (c) を数えて出す。
void quantize(cudaConfig& cuda_cfg, mesh& msh, variables& var, int iStep);

// updateVariablesOuter の直後。iStep は 0 始まり。自己点検の step なら (a) を数えて出す。
void afterUpdateOuter(mesh& msh, variables& var, int iStep);

// 化学種の更新の頭 (implicitNonlinearUpdate の speciesUpdateOuter_d_wrapper の直後)。iStep は 0 始まり。
// 前の step が自己点検の step なら (b) を数えて出す。
void afterSpeciesUpdateOuter(mesh& msh, variables& var, int iStep);

// 終了時。点検しなかった対象の step があればその旨を出す。
void finalize();

}  // namespace quantizeStateDiag
