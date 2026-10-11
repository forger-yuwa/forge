#pragma once

// 閉性の照合の診断 (plans/active/axisymmetric-freestream-hoop-gauge.md §4.6 の 3・5)。既定 off。
//
//   FORGE_DIAG_HOOP_CLOSURE=<out.h5>
//     起動の段取り (構造量・周期・ライン) が済んだ時点で 1 回だけ、デバイスに渡した最終の面ベクトル sx..sz・面積 ss・
//     A_planar・体積を読み、各 CV (実 CV、ゴーストを除く) の閉性を double で計算して <out.h5> に書く。
//     書いた後も計算は**続ける** (止めない)。壁・軸の射影 (nodeWallDirichlet・軸の Dirichlet) とは独立の、幾何だけの量。
//     面の向きは plane_cells (デバイスの map_plane_cells_d) の ic0 に +、ic1 に − (ic1 がゴーストの境界面は ic0 だけ)。
//
// 書く量 (h5 のパス):
//   /faces/{sx,sy,sz,ss} [nPlanes]   デバイスの値そのもの (flow_float)。新旧のバイナリ・キー 0/1 のビット比較用
//   /faces/plane_cells [2*nPlanes]   デバイスの接続
//   /faces/normlen_err [nPlanes]     ‖(sx,sy,sz)‖/ss − 1 (double)
//   /cells/{A_planar,volume} [nCells_all]  デバイスの値そのもの (A_planar は軸対称のときだけ)
//   /cells/{sum64_x,sum64_y} [nCells]      Σ±S_x、Σ±S_y (デバイスの値を double で足す)
//   /cells/{abs64_x,abs64_y} [nCells]      Σ|S_x|、Σ|S_y|
//   /cells/{E_x,E_y} [nCells]              E_x = |Σ±S_x|/A、E_y = |Σ±S_y − T_y|/A。A = A_planar (軸対称) か体積、
//                                          T_y = A_planar (軸対称の r 重み axisymMethod 0) か 0
//   /cells/cc64 [3*nCells]                 値の位置 (double の写し。領域の切り分け用、あれば)
// 属性: 設定 (isAxisymmetric・axisymMethod・axisRFloor・hoopAreaFromClosure・axisSegmentRWeight)、使った経路
// (segment_rweight_active)、ss の非有限・非正の数、法線の長さの検査の最大、E_x・E_y の最大、
// 100·ε64·(A + Σ|S の成分|) を超える CV の数 (FP64 のビルドでだけ意味がある丸めの規模)。

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

namespace hoopClosureDiag {

// FORGE_DIAG_HOOP_CLOSURE が設定されていれば書いて戻る (計算は続ける)。未設定なら何もしない。
void runIfRequested(const solverConfig& cfg, const mesh& msh, variables& var);

}  // namespace hoopClosureDiag
