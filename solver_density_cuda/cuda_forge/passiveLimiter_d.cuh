#pragma once
// 受動種 (トレーサ・凝縮モーメント) 用のセル局所スケール Venkatakrishnan リミッタ
// (plans/active/species-passive-scalar-unification.md §4.1, codex plan-2 M1)。
//
// limiter_r1_d (limiter_d.cu) と同じ MLP 型の走査だが、差分をセル局所スケール
//   φ_ref = max(|φ_c|, max_nb |φ_nb|, φ_floor)
// で無次元化 (Δ̃ = Δ/φ_ref) してから同じ Venkat / Barth-Jespersen 関数 (limiterFunctions_d.cuh) に渡す。
// 現行 limiter_r1_d は Δ_+²Δ_- 等の 3 次積を float32 で組むため、モーメント Q0 ~ 1e15 では分子・分母が Inf
// (商は NaN) になる。無次元化後は |Δ̃| ≤ 2 なので 3 次積は常に有限。Venkat 関数はスケール不変
// (ε² を除く) なので ψ は変わらない。ε² は現行と同じ引数 (K³·V) を渡す (無次元差分に対しては実質的に小さい正則化)。
// limiter_d.cu (本番) と tests/unit/test_passive_scalar.cu (有限性の単体試験) が include する。
#include "flowFormat.hpp"

__global__ void limiter_r1_scaled_d
(
 int limiter_scheme,
 geom_int nCells,
 geom_int nNormalPlanes, geom_int* plane_cells,
 geom_int* cell_planes_index, geom_int* cell_planes,
 geom_float* vol, geom_float* ccx, geom_float* ccy, geom_float* ccz,
 geom_float* pcx, geom_float* pcy, geom_float* pcz,
 flow_float phi_floor,
 flow_float* Q, flow_float* limiter_Q,
 flow_float* dQdx, flow_float* dQdy, flow_float* dQdz
)
{
    const geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic0 >= nCells) return;
    if (limiter_scheme == 0) { limiter_Q[ic0] = 1.0f; return; }

    const geom_int index_st = cell_planes_index[ic0];
    const geom_int index_en = cell_planes_index[ic0+1];
    const flow_float Qc = Q[ic0];
    flow_float Q_max = Qc, Q_min = Qc;
    for (geom_int ilp = index_st; ilp < index_en; ++ilp) {
        const geom_int ip = cell_planes[ilp];
        if (ip >= nNormalPlanes) continue;
        const geom_int ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0;
        Q_max = max(Q_max, Q[ic1]);
        Q_min = min(Q_min, Q[ic1]);
    }
    // セル局所スケール (常に > 0)。
    const flow_float phi_ref = max(max(fabsf(Qc), max(fabsf(Q_max), fabsf(Q_min))), phi_floor);
    const flow_float inv_ref = 1.0f / phi_ref;
    const flow_float dp_max = (Q_max - Qc) * inv_ref;
    const flow_float dp_min = (Q_min - Qc) * inv_ref;
    const flow_float volume = vol[ic0];
    const flow_float gx = dQdx[ic0], gy = dQdy[ic0], gz = dQdz[ic0];

    flow_float lim = 1.0f;
    for (geom_int ilp = index_st; ilp < index_en; ++ilp) {
        const geom_int ip = cell_planes[ilp];
        if (ip >= nNormalPlanes) continue;
        const flow_float dcp_x = pcx[ip] - ccx[ic0];
        const flow_float dcp_y = pcy[ip] - ccy[ic0];
        const flow_float dcp_z = pcz[ip] - ccz[ic0];
        const flow_float delta_m = (gx*dcp_x + gy*dcp_y + gz*dcp_z) * inv_ref;   // (Qt − Qc)/φ_ref
        flow_float l;
        if (limiter_scheme == 1) l = barth_Jespersen_limiter(dp_max, dp_min, delta_m, volume);
        else                     l = venkata_limiter(dp_max, dp_min, delta_m, volume);
        lim = min(lim, l);
    }
    limiter_Q[ic0] = min(max(lim, 0.0f), 1.0f);
}
