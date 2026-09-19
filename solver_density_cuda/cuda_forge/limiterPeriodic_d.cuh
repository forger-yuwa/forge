#pragma once
#include "cuda_forge/reconIncrement_d.cuh"   // 再構成増分の唯一の定義 (流束と共有)
// 周期 node (median-dual の合併 CV) 用の 2 段リミッタ (plans/active/species-passive-scalar-unification.md §4.8, codex result-2 M3)。
//
// 周期 group (root + member) は 1 つの合併 CV だが、1 段の limiter_r1_d / limiter_r1_scaled_d は各 member の**部分 CV** の
// 内部面だけで極値 Q_max/Q_min・スケール φ_ref・ψ を作るので、状態を mirror しても周期対で ψ が一致しない
// (case/09 run_0089 の 4532/2982: limiter_Xi 0.985 vs 1.0)。ここでは
//   段 1 (limiter_extrema_d)      : 各 CV の極値 (自身 + 内部面隣接) を配列へ → 呼び出し側で周期 group の max/min gather + broadcast
//   段 2 (limiter_psi_merged_d)   : 合併極値 (受動種は φ_ref も合併極値から) で各 member の**自分の面**の ψ を評価 → 呼び出し側で group の min gather + broadcast
// に分ける。合併 CV の面集合 = 各 member の面の和なので、(max/min の合併, min の合併) は「合併 CV に 1 段 kernel を掛けた結果」と
// 演算順序まで一致する (tests/unit/test_periodic_limiter.cu で巻き付き鎖と等価性を検査)。
// 非周期・cell では limiter_d.cu の 1 段 kernel をそのまま使う (ビット不変)。
// SCALED=false は limiter_r1_d と同じ式 (Qt = Q + ∇Q·dcp, Δ = Qt − Q)、SCALED=true は limiter_r1_scaled_d と同じ式
// (Δ̃ = (∇Q·dcp)/φ_ref, φ_ref = max(|Qc|,|Q_max|,|Q_min|,φ_floor))。
// __global__ カーネルを持つので、ライブラリ側は limiter_d.cu からだけ include する (単体試験 TU は別実行体)。
#include "flowFormat.hpp"

// 段 1: 極値 (自身 + 内部面 ip < nNormalPlanes の隣接)。ghost・周期半割面は入らない。
__global__ void limiter_extrema_d
(
 geom_int nCells,
 geom_int nNormalPlanes, geom_int* plane_cells,
 geom_int* cell_planes_index, geom_int* cell_planes,
 flow_float* Q, flow_float* Q_max_out, flow_float* Q_min_out
)
{
    const geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic0 >= nCells) return;
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
    Q_max_out[ic0] = Q_max;
    Q_min_out[ic0] = Q_min;
}

// 段 2: 合併極値 (Q_max/Q_min: gather 済み) で自分の面の ψ を評価。limiter_scheme 1 = Barth-Jespersen, それ以外 = Venkatakrishnan。
template<bool SCALED>
__global__ void limiter_psi_merged_d
(
 int limiter_scheme,
 geom_int nCells,
 geom_int nNormalPlanes,
 geom_int* plane_cells,
 geom_int* cell_planes_index, geom_int* cell_planes,
 geom_float* vol, geom_float* ccx, geom_float* ccy, geom_float* ccz,
 geom_float* pcx, geom_float* pcy, geom_float* pcz,
 flow_float phi_floor,
 flow_float* Q, flow_float* Q_max_in, flow_float* Q_min_in,
 flow_float* limiter_Q,
 flow_float* dQdx, flow_float* dQdy, flow_float* dQdz,
 // 通常経路と同じ「流束一致」オプション (plan convection-node-wall-reconstruction §4.8)。既定 0 で式は変更前と同一。
 int matchRecon, int edgeMid, int convM,
 // 無次元化 Venkatakrishnan (plan §4.13 / codex plan-3 Critical 1)。既定 limScaled=0 で式は変更前と同一。
 int limScaled, flow_float qRef, flow_float eps2Coef, int lenArea, geom_float* A_planar
)
{
    const geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic0 >= nCells) return;
    if (limiter_scheme == 0) { limiter_Q[ic0] = 1.0f; return; }

    const geom_int index_st = cell_planes_index[ic0];
    const geom_int index_en = cell_planes_index[ic0+1];
    const flow_float Qc = Q[ic0];
    const flow_float Q_max = Q_max_in[ic0];
    const flow_float Q_min = Q_min_in[ic0];
    const flow_float volume = vol[ic0];
    const flow_float gx = dQdx[ic0], gy = dQdy[ic0], gz = dQdz[ic0];

    flow_float inv_ref = 1.0f;
    flow_float dp_max, dp_min;
    if (SCALED) {
        const flow_float phi_ref = max(max(fabsf(Qc), max(fabsf(Q_max), fabsf(Q_min))), phi_floor);
        inv_ref = 1.0f / phi_ref;
        dp_max = (Q_max - Qc) * inv_ref;
        dp_min = (Q_min - Qc) * inv_ref;
    } else {
        dp_max = Q_max - Qc;
        dp_min = Q_min - Qc;
    }

    flow_float lim = 1.0f;
    for (geom_int ilp = index_st; ilp < index_en; ++ilp) {
        const geom_int ip = cell_planes[ilp];
        if (ip >= nNormalPlanes) continue;
        flow_float dcp_x, dcp_y, dcp_z;
        geom_int ic1p = -1;
        if (matchRecon != 0 && edgeMid != 0) {        // node: 目標点 = エッジ中点 (流束と同じ)
            ic1p = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0;
            dcp_x = (flow_float)0.5*(ccx[ic1p] - ccx[ic0]);
            dcp_y = (flow_float)0.5*(ccy[ic1p] - ccy[ic0]);
            dcp_z = (flow_float)0.5*(ccz[ic1p] - ccz[ic0]);
        } else {
            dcp_x = pcx[ip] - ccx[ic0];
            dcp_y = pcy[ip] - ccy[ic0];
            dcp_z = pcz[ip] - ccz[ic0];
        }
        flow_float delta_m;
        if (matchRecon != 0) {
            if (ic1p < 0) ic1p = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0;
            // 流束と**同じ関数**で増分を作る (reconIncrement_d.cuh)。Qt を経由しないので桁落ちも無い
            delta_m = recon_increment(convM, Qc, Q[ic1p], gx, gy, gz, dcp_x, dcp_y, dcp_z);
            if (SCALED) delta_m *= inv_ref;
        } else if (SCALED) {
            delta_m = (gx*dcp_x + gy*dcp_y + gz*dcp_z) * inv_ref;              // limiter_r1_scaled_d と同式
        } else {
            const flow_float Qt = Qc + gx*dcp_x + gy*dcp_y + gz*dcp_z;          // limiter_r1_d と同式
            delta_m = Qt - Qc;
        }
        flow_float l;
        if (limScaled == 2 && limiter_scheme != 1) {
            // 比の形 (基準値も長さも不要)
            l = (fabsf(delta_m) > (flow_float)1.0e-20)
              ? venkata_limiter_ratio(dp_max, dp_min, delta_m, eps2Coef) : (flow_float)1.0;
        } else if (limScaled == 1 && limiter_scheme != 1) {
            // 変数ごとの固定参照で無次元化してから Venkatakrishnan (通常経路 limiter_r1_fused5_d と同式)
            const flow_float inv = (flow_float)1.0/qRef;
            const flow_float hi  = (lenArea != 0) ? sqrtf(A_planar[ic0]) : cbrtf(volume);
            const flow_float e2  = eps2Coef * hi*hi*hi;
            l = venkata_limiter_scaled(dp_max*inv, dp_min*inv, delta_m*inv, e2);
        } else if (limiter_scheme == 1) {
            l = barth_Jespersen_limiter(dp_max, dp_min, delta_m, volume);
        } else {
            l = venkata_limiter(dp_max, dp_min, delta_m, volume);
        }
        lim = min(lim, l);
    }
    limiter_Q[ic0] = min(max(lim, 0.0f), 1.0f);
}
