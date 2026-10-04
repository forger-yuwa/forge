#pragma once

// =============================================================================
// node 等温壁のエネルギー境界を弱形式 (SU2 型) で課すための幾何 (mesh.nodeIsothermalEnergyBC=1)
//
// 仕様は methods/boundary.md「等温壁のエネルギー境界: 強制 と 弱形式」、
// 設計判断は plans/active/boundary-weak-isothermal-wall.md。
//
// 壁半割面ごとに、**第一内部点の DOF index** と**法線投影距離 d_1** を device に置く。
// 規約は conjugateWall::firstInterior と**同一** (数字が食い違うと切り分け不能になるため、
// 診断 q_compact と BC が同じ幾何を見ることを保証する)。
// **点間距離ではない** — 法線から 30° 傾いた辺で 13.4 % ずれる (codex plan レビュー M2)。
//
// 幾何は初期化時に 1 度だけ確定し、残差カーネルと陰解法が同じ値を参照する。
// 評価不能 (内部点なし / 整列度不足 / d_1 退化 / 選ばれた点が別の壁ノード) は**起動時に拒否**する。
// =============================================================================

#include "../flowFormat.hpp"
#include "../input/solverConfig.hpp"
#include "../mesh/mesh.hpp"

namespace weakIsoWall {

struct Geom {
    geom_int*   j_d  = nullptr;   // 第一内部点の DOF index (bplane 順)
    flow_float* d1_d = nullptr;   // 法線投影距離 [m]
    geom_int    n    = 0;
};

// 有効な構成か (node × nodeIsothermalEnergyBC=1 × wall_isothermal が存在)。
bool active(const solverConfig& cfg, const mesh& msh);

// 起動時の構成検査 (併用不可の組合せを拒否)。active でなければ何もしない。
void validate(const solverConfig& cfg, const mesh& msh);

// bcond ごとの幾何を初回に構築してキャッシュ。active でなければ空の Geom を返す。
const Geom& geom(const solverConfig& cfg, const mesh& msh, const bcond& bc);

// 陰解法の近似対角項の素材 g = Σ_{壁半割面} k_eff A_half / d_1 [W/K] (DOF ごと)。
// **残差カーネルと同じ k_eff・同じ幾何**で毎ステップ積む (codex plan レビュー M1/M2)。
// 対角への加算は timeIntegration 側で +g/(rho c_v) (forge は残差微分の符号を反転して組むので正号)。
flow_float* diagBuf(const mesh& msh);     // active でなければ nullptr
void        diagReset(const mesh& msh);   // viscousFlux の壁ループ前に 0 クリア

}  // namespace weakIsoWall
