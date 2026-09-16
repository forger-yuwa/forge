#pragma once
// 化学種経路を受動種 (排気トレーサ・凝縮モーメント) と共用するカーネル群
// (plans/active/species-passive-scalar-unification.md §4.1)。speciesTransport_d.cu (本番) と
// tests/unit/test_passive_scalar.cu (単体試験) が include する。__global__ を含むので他のライブラリ TU からは include しない。
#include "flowFormat.hpp"

// S3: species 移流流束を **convectiveFlux が書いた face 組成** で組む (energy 流束と同一面組成)。
// 対角 transport_diag は 1 次風上のまま (defect-correction)。ΣY_face=1 なので Σ res_roY = res_ro。
// 受動種は同じカーネルを受動種ポインタ (nSpecies=nPassive, Yface=Pface) で呼ぶ。
__global__ void species_advection_faceY_d(
    geom_int nCells, geom_int nNormalHaloPlanes, geom_int* normal_halo_planes, geom_int* plane_cells,
    flow_float* ro, flow_float* massflux, int nSpecies, flow_float* Yface,
    flow_float** res_roY, flow_float** transport_diag,
    int isNode, flow_float** roY,
    int stride)   // Yface の面ストライド (化学種 = nSpecies; 受動種は部分範囲を Yface+q0, stride=nPassive で呼ぶ)
{
    geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih < nNormalHaloPlanes) {
        const geom_int ip  = normal_halo_planes[ih];
        const geom_int ic0 = plane_cells[2*ip+0];
        const geom_int ic1 = plane_cells[2*ip+1];
        const flow_float mdot = massflux[ip];
        const flow_float d0 = max(mdot, (flow_float)0.0) / max(ro[ic0], (flow_float)1.0e-30);
        const flow_float d1 = max(-mdot,(flow_float)0.0) / max(ro[ic1], (flow_float)1.0e-30);
        // node 境界半割面 (ic1=ghost): 主ループ (SLAU/ROE) は境界半割面を除外するため Yface[ip] が
        // 未書込 (stale)。node は ghost を読まない設計なので、境界ノード ic0 自身の組成を面組成に使う。
        const bool nodeBnd = (isNode != 0 && ic1 >= nCells);
        for (int s = 0; s < nSpecies; ++s) {
            const flow_float Yf = nodeBnd
                ? (roY[s][ic0] / max(ro[ic0], (flow_float)1.0e-30))
                : Yface[(size_t)ip*stride + s];     // 内部面 upwind は convectiveFlux 側で確定済み
            const flow_float flux = mdot * Yf;
            if (ic0 < nCells) { atomicAdd(&res_roY[s][ic0], -flux); atomicAdd(&transport_diag[s][ic0], d0); }
            if (ic1 < nCells) { atomicAdd(&res_roY[s][ic1],  flux); atomicAdd(&transport_diag[s][ic1], d1); }
        }
    }
}

// 受動種の更新確定時の上下限 (codex plan-2 M2): 更新済み密度 ro で 0 <= ρφ (<= ρ: upperIsRho) を適用し、
// 補正収支を記録する。
//   corrCell[ic] += |Δ(ρφ)|                      (セル配列 passiveFloorCorr_<prim>: 全期間の累積)
//   stats[0] += Σ Δ(ρφ)·V (下限, ≥0), stats[1] += Σ Δ(ρφ)·V (上限, ≤0), stats[2] += Σ|Δ(ρφ)|·V (累積),
//   stats[3] += Σ (ρφ)_after·V (現在の総量; 呼び出し側が毎回 0 にしてから呼ぶ)。
// 面クリップ (S3) は保存的なので保存量の補正はここだけ (更新 floor)。周期 node では合併体積 V が両側で
// 二重計上されるが、相対量 (補正/総量) は不変。
// root (node 周期の periodicRoot; 非周期は nullptr): 収支・総量は **周期 root だけを合併体積で** 数える (codex result M2;
// member を数えると合併 CV が二重計上される)。クランプ自体は全ノードに掛ける (member は後段のミラーで root に揃う)。
__global__ void passive_bounds_d(
    geom_int nCells, flow_float* rophi, int upperIsRho, flow_float* ro, geom_float* vol,
    flow_float* corrCell, double* stats, const geom_int* root)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    double lo = 0.0, hi = 0.0, ab = 0.0, tot = 0.0;
    if (ic < nCells) {
        const bool count = (root == nullptr) || (root[ic] == ic);
        const flow_float v0 = rophi[ic];
        flow_float v = v0;
        if (v < (flow_float)0.0) v = (flow_float)0.0;
        if (upperIsRho != 0) { const flow_float r = ro[ic]; if (v > r) v = r; }
        const double d = (double)v - (double)v0;
        const double V = (double)vol[ic];
        if (d != 0.0) {
            rophi[ic] = v;
            if (count && corrCell != nullptr) {   // corrCell==nullptr: 収支を記録しない (RK の中間ステージ)
                corrCell[ic] += (flow_float)fabs(d);
                if (d > 0.0) lo = d*V; else hi = d*V;
                ab = fabs(d)*V;
            }
        }
        if (count) tot = (double)v*V;
    }
    // block 縮約 → 1 block 1 回の atomicAdd (double)。
    __shared__ double sh[4][32];
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    for (int off = 16; off > 0; off >>= 1) {
        lo  += __shfl_down_sync(0xffffffffu, lo,  off);
        hi  += __shfl_down_sync(0xffffffffu, hi,  off);
        ab  += __shfl_down_sync(0xffffffffu, ab,  off);
        tot += __shfl_down_sync(0xffffffffu, tot, off);
    }
    if (lane == 0) { sh[0][wid] = lo; sh[1][wid] = hi; sh[2][wid] = ab; sh[3][wid] = tot; }
    __syncthreads();
    if (threadIdx.x == 0) {
        const int nw = (blockDim.x + 31) >> 5;
        double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        for (int w = 0; w < nw; ++w) { a0 += sh[0][w]; a1 += sh[1][w]; a2 += sh[2][w]; a3 += sh[3][w]; }
        if (a0 != 0.0) atomicAdd(&stats[0], a0);
        if (a1 != 0.0) atomicAdd(&stats[1], a1);
        if (a2 != 0.0) atomicAdd(&stats[2], a2);
        if (a3 != 0.0) atomicAdd(&stats[3], a3);
    }
}

// 受動種の増分スケーリング (codex result M5; 凝縮の θ_u と同じ概念): 更新後の候補 ρφ = ρφ_N + δ に対し
//   θ_b = min(1, allowed/|δ|),  allowed = (δ>0: 上限 [ρ, upperIsRho] − ρφ_N ; δ<0: ρφ_N − 0)
// で増分全体を縮め、確定状態を [0, ρ] (トレーサ) / ≥0 (モーメント) に保つ。ρφ_N 自体が範囲外 (流れ更新で ρ が減った等) なら
// θ_b=0 で N に留め、その後の硬い floor (passive_bounds_d, 最後の砦) が処理する。
//   limCell[ic] += (1−θ_b)|δ| (セル累積), stats[0] += Σ(1−θ_b)|δ|·V, stats[1] += 作動セル数, stats[2] += Σ(θ_b−1)δ·V (符号付き) (root のみ), thetaMinInt = min(θ_b·1e9)。
// 流れの密度更新と整合した受動種の増分 (plan §5.1 #19, 案C の ρY_N + z + Y_N δρ と同形): 候補 ρφ = ρφ_N + z に
//   φ_N·δρ = (ρφ_N/ρ_pre)·(ρ_new − ρ_pre)
// を加える (ρ_pre = 同じ (サブ) 反復の残差組み立て時の ρ, ρ_new = 流れ block 更新後の ρ)。z=0 なら ρφ = φ_N ρ_new で φ は不変。
__global__ void passive_add_rho_term_d(geom_int nCells, flow_float* rophi, const flow_float* rophiN, const flow_float* roPre, const flow_float* ro)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double rp = (double)roPre[ic];
    if (rp <= 0.0) return;
    rophi[ic] = (flow_float)((double)rophi[ic] + ((double)rophiN[ic]/rp)*((double)ro[ic] - rp));
}

// roPre != nullptr のとき増分の基点は φ_N ρ_new (= ρφ_N + φ_N δρ) で、制限するのは輸送増分 z だけ (基点自体は [0,ρ_new] 内)。
__global__ void passive_limit_increment_d(
    geom_int nCells, flow_float* rophi, const flow_float* rophiN, int upperIsRho, const flow_float* ro, const geom_float* vol,
    flow_float* limCell, double* stats, int* thetaMinInt, const geom_int* root, const flow_float* roPre)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    double N = (double)rophiN[ic];
    if (roPre != nullptr && roPre[ic] > (flow_float)0.0) N = N / (double)roPre[ic] * (double)ro[ic];
    const double d = (double)rophi[ic] - N;
    if (d == 0.0) return;
    double allowed;
    if (d > 0.0) allowed = (upperIsRho != 0) ? ((double)ro[ic] - N) : 1.0e300;
    else         allowed = N;
    if (allowed < 0.0) allowed = 0.0;
    double th = allowed / fabs(d);
    if (th >= 1.0) return;
    rophi[ic] = (flow_float)(N + th*d);
    const bool count = (root == nullptr) || (root[ic] == ic);
    if (count && limCell != nullptr) {   // limCell==nullptr: 収支を記録しない (RK の中間ステージ)
        const double amt = (1.0 - th)*fabs(d);
        limCell[ic] += (flow_float)amt;
        atomicAdd(&stats[0], amt*(double)vol[ic]);          // 絶対量 Σ(1−θ)|δ|V
        atomicAdd(&stats[1], 1.0);                          // 作動セル数
        atomicAdd(&stats[2], (th - 1.0)*d*(double)vol[ic]); // 符号付き (確定 − 候補) Σ(θ−1)δV: 保存収支の説明に使う
        atomicMin(thetaMinInt, (int)(th*1.0e9));
    }
}

// 受動種 (トレーサ) の Fick 拡散: J = ρ_f D (φ1−φ0)/dcc·δ (over-relaxed 法線), D = μ_f/(ρ_f Sc) + μt_f/(ρ_f Sc_t)。
// species_diffusion_d (定数 Schmidt 分岐) と同じ面幾何・同じ float32 評価だが ΣJ=0 補正・エンタルピー結合は無い
// (受動種は熱力学に入らない)。node 境界半割面は skip (species と同方針: Dirichlet はピン、Neumann は流束 0)。
__global__ void passive_diffusion_d(
    geom_int nCells, geom_int nNormalHaloPlanes, geom_int* normal_halo_planes, geom_int* plane_cells,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    flow_float* phi, flow_float* res, flow_float* transport_diag,
    flow_float* ro, flow_float* vis_lam, flow_float* vis_turb,
    flow_float Sc, flow_float Sc_t, int isNode)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0];
    const geom_int ic1 = plane_cells[2*ip+1];
    if (isNode != 0 && (ic0 >= nCells || ic1 >= nCells)) return;

    const flow_float f = fx[ip], g = 1.0f - f;
    const flow_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];
    const flow_float dccx = ccx[ic1] - ccx[ic0];
    const flow_float dccy = ccy[ic1] - ccy[ic0];
    const flow_float dccz = ccz[ic1] - ccz[ic0];
    const flow_float dcc  = sqrtf(dccx*dccx + dccy*dccy + dccz*dccz);
    const flow_float denom = dccx*sxx + dccy*syy + dccz*szz;
    const flow_float Dsafe = (fabsf(denom) < 1.0e-30f) ? ((denom>=0.0f)?1.0e-30f:-1.0e-30f) : denom;
    const flow_float delta = dcc * sss * sss / Dsafe;

    const flow_float ro0 = max(ro[ic0], (flow_float)1.0e-30f);
    const flow_float ro1 = max(ro[ic1], (flow_float)1.0e-30f);
    const flow_float ro_face  = f*ro0 + g*ro1;
    const flow_float mu_face  = f*vis_lam[ic0]  + g*vis_lam[ic1];
    const flow_float mut_face = f*vis_turb[ic0] + g*vis_turb[ic1];
    const flow_float Dt = (mut_face > 0.0f) ? mut_face/(ro_face*Sc_t) : 0.0f;
    const flow_float D  = mu_face/(ro_face*Sc) + Dt;
    const flow_float roD = ro_face * D;
    const flow_float inv_dcc = 1.0f/max(dcc, (flow_float)1.0e-30f);
    const flow_float J = roD * ((phi[ic1] - phi[ic0])*inv_dcc) * delta;
    const flow_float diag = roD * fabsf(delta) * inv_dcc;
    if (ic0 < nCells) { atomicAdd(&res[ic0],  J); atomicAdd(&transport_diag[ic0], diag/ro0); }
    if (ic1 < nCells) { atomicAdd(&res[ic1], -J); atomicAdd(&transport_diag[ic1], diag/ro1); }
}

// a[ic] += b[ic] (周期 gather 済みのクロス項を流れ RHS へ加える等)。
__global__ void passive_axpy1_d(geom_int n, flow_float* a, const flow_float* b)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < n) a[ic] += b[ic];
}

// 受動種の総量 Σ ρφ V (周期 root のみ; double) — 収支の独立照合 (開始前 / 全後処理後の確定状態) 用。
__global__ void passive_total_d(geom_int nCells, const flow_float* rophi, const geom_float* vol, const geom_int* root, double* out)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    double v = 0.0;
    if (ic < nCells && (root == nullptr || root[ic] == ic)) v = (double)rophi[ic]*(double)vol[ic];
    __shared__ double sh[32];
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    for (int off = 16; off > 0; off >>= 1) v += __shfl_down_sync(0xffffffffu, v, off);
    if (lane == 0) sh[wid] = v;
    __syncthreads();
    if (threadIdx.x == 0) { const int nw = (blockDim.x + 31) >> 5; double a = 0.0; for (int w = 0; w < nw; ++w) a += sh[w]; if (a != 0.0) atomicAdd(out, a); }
}
