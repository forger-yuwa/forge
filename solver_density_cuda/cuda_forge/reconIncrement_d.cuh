#pragma once
// 再構成増分 Δ の**唯一の定義**。`phi_face = phiC + psi * Δ` の Δ を返す。
//
// 流束 (convectiveFlux_common_d.cuh の interp_*) と**リミッタ** (limiter_d.cu / limiterPeriodic_d.cuh) が
// 同じ関数を同じ引数で呼ぶためにここへ切り出した。両者で式や目標点が食い違うと、
// リミッタが保証する有界性は**実際に流束が使う値には掛からない**
// (plan convection-node-wall-reconstruction §4.8/§4.14)。
//
// `cpd` は目標点までのオフセット。**node は常にエッジ中点** (`±0.5*dcc`。convectiveFlux_d.cu の
// `g_reconEdgeMid`)、cell は双対面重心 `pc[ip] - cc[ic]`。
//
// 対象は `convMethod` 0 / 1 / 2 のみ。MINMOD (その他) は別式なので呼び出し側で除外すること。
#include "flowFormat.hpp"

__device__ __forceinline__ flow_float recon_increment(
    int scheme,
    flow_float phiC, flow_float phiD,
    flow_float dphidx, flow_float dphidy, flow_float dphidz,
    flow_float cpdx, flow_float cpdy, flow_float cpdz)
{
    if (scheme == 0 || scheme == -1) return (flow_float)0.0;               // 1 次風上: 増分なし
    const flow_float proj = dphidx*cpdx + dphidy*cpdy + dphidz*cpdz;
    if (scheme == 2) {                                                      // 3 次 MUSCL (k = 1/3)
        const flow_float k = (flow_float)(1.0/3.0);
        return (flow_float)0.5*k*(phiD-phiC) + ((flow_float)1.0-k)*proj;
    }
    return proj;                                                            // 2 次 MUSCL
}
