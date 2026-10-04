#pragma once
// 凝縮 carrier (TP, 凝縮種は輸送種の 1 つ) の気相組成 (plans/active/condensation-two-phase-transport.md §4.1)。
//   輸送物性 (μ・λ・化学種拡散係数) は気相組成で評価する: 液 g = ρg/ρ を凝縮種の総水分から除き、
//     Y_s^gas = Y_s/(1−g) (s ≠ 凝縮種),  Y_w^gas = (Y_w − g)/(1−g)  (負は 0 に切る)。
//   (1−g) での再正規化は各入口の既存の正規化 (Y の和・X の和) が行うので、ここでは凝縮種から液を引いて 0 で切るだけにする。
//   全経路 (セルの表引き gas_transport_cell_rY・直接評価 gas_transport_cell_Y・壁 wmlesWallModel_d・FORGE_TRANSPORT_PROBE・
//   化学種拡散係数 species_diffusion_d) がこの関数を通る。液が無い (iw < 0 / liq ≤ 0) ときは何もしない = 現行とビット一致。
//   CPG carrier (空気凝縮) と pure 凝縮はここを通らない (iw = -1)。
#include "flowFormat.hpp"

// Y: 種ごとの量 (正規化前でも可; 単位は liq と揃える: ρY なら liq = ρg、Y なら liq = g)。
// iw: 凝縮種の輸送種 index (-1 = 液の除去なし)。
template <typename R>
__host__ __device__ inline void gas_phase_composition(R* Y, int iw, R liq)
{
    if (iw < 0 || !(liq > R(0))) return;
    const R v = Y[iw] - liq;
    Y[iw] = (v > R(0)) ? v : R(0);
}

// kernel に渡す液の所在 (host が gasPhaseLiquid で作る)。rog == nullptr / iw < 0 は液の除去なし。
struct GasPhaseLiquid {
    const flow_float* rog;   // 凝縮種の液保存量 ρg (nCells_all; ghost は凝縮の境界処理が埋める)
    int               iw;    // 凝縮種の輸送種 index
};

// 液の単位を揃えた取り出し (ρ 単位)。液が無い構成は 0。
__host__ __device__ inline flow_float gas_phase_liquid_rho(const GasPhaseLiquid& L, geom_int ic)
{
    return (L.rog != nullptr && L.iw >= 0) ? L.rog[ic] : (flow_float)0.0;
}

class solverConfig;
class variables;
// TP carrier の凝縮 run (thermalMethod 2・condGasSpecies ≥ 0・凝縮種登録あり) なら {rog_0, condGasSpecies}、それ以外は {nullptr, -1}。
GasPhaseLiquid gasPhaseLiquid(const solverConfig& cfg, variables& var);
