#pragma once
// 受動種 (トレーサ・凝縮モーメント) の dual-time 物理 step 末尾の保存的 FCT 補正 (Zalesak 型;
// plans/active/species-passive-scalar-unification.md §4.7 v4, codex result-2 M2 / plan-3〜5)。
// speciesTransport_d.cu (本番 wrapper) と tests/unit/test_passive_fct.cu (単体試験) が include する。__global__ を含むので他の TU からは include しない。
//
// 記法: 面 ip=(ic0,ic1), ṁ>0 は ic0→ic1。面の ic0→ic1 向き全流束は F = ṁ φ_f − J_f (J は既存カーネルが ic0 へ加算する拡散流束)。
// 残差規約: 面流束 F の CV への寄与は s_f F (s = −1 [ic0], +1 [ic1])。境界半割面は ic0 だけ (外部との交換)。M = V/Δt。
//   BDF2 を流束形の BE にする (b = a + c):  (V/Δt)(q^{n+1} − q^n) = Σ_f s_f F^eff_f + S^eff V,
//     F^eff = (1/a) F(q^{n+1}) + (c/a) G^n,   S^eff V = (1/a) S V + (c/a) H^n,
//     G^n, H^n = 前 step の増分率 I^n = (V/Δt)(q^n − q^{n−1}) の面流束/局所への厳密分解 (確定状態から作る)。BDF1 は a=1, c=0。
//   q_H : sub-iter 終了状態 (終了状態で残差を再評価して ṁ, P_face, ソースを固定)
//   q_L : 低次 BE 陰解  (V/Δt)(q_L − q^n) = Σ s F_L(q_L) + S^eff V  (Jacobi; 受入は線形残差)
//   A^raw_f = F^eff_H,f − F_L,f(q_L)  (全輸送面)、基点 q_B = q_H − (Δt/V) Σ s A^raw = q_L + M^{-1}(r_L − r_H/a)
//   Zalesak (全輸送面; 境界は外部側 R=1): 前制限 A^pre → P± → 局所極値 {φ^n, φ_B} → Q±, R± → α_f
//   補正 q_C = q_H − (Δt/V) Σ s (A^raw − α A^pre);  実現流束 G^{n+1}_f = F_L,f(q_L) + α_f A^pre_f、H^{n+1} = (V/Δt)(q_final − q^n) − Σ s G^{n+1}
// 周期 node の gather (diag/nb/P±/base/corr/ΣsG は和, 極値は max/min) は呼び出し側 (wrapper) が行う。
#include "flowFormat.hpp"

// トレーサ Fick 拡散の面係数 c_f (passive_diffusion_d と同じ float32 評価): J = c_f (φ1 − φ0)。
__device__ __forceinline__ flow_float passive_fct_diff_coef(
    geom_int ip, geom_int ic0, geom_int ic1,
    const geom_float* ccx, const geom_float* ccy, const geom_float* ccz,
    const geom_float* fx, const geom_float* sx, const geom_float* sy, const geom_float* sz, const geom_float* ss,
    const flow_float* ro, const flow_float* vis_lam, const flow_float* vis_turb, flow_float Sc, flow_float Sc_t)
{
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
    const flow_float inv_dcc = 1.0f/max(dcc, (flow_float)1.0e-30f);
    return ro_face * D * inv_dcc * delta;
}

// 低次面流束 F_L (低次作用素は流れの BE 形連続式と整合する有効質量流束 ṁ^eff = (1/a)ṁ^{n+1} + (c/a)ṁ^eff,n で組む; plan-6 M1):
//   内部面: ṁ^eff φ_L,up − J_L。node 境界半割面: 流出 (ṁ^eff ≥ 0) は ṁ^eff φ_L,own、流入 (<0) は許容な外部状態 φ^n_own = q^n/ρ^n を定数 RHS に (plan-6 M2:
//   自セル値を対角に入れると M 行列でなくなる)。cell の ghost 面: 流入は ghost 値 (定数), 流出は自セル。
__device__ __forceinline__ flow_float passive_fct_lowflux(
    geom_int nCells, geom_int ip, geom_int ic0, geom_int ic1, flow_float meff, const flow_float* ro, const flow_float* roN,
    int q, int qDiff, const flow_float* cdiff, geom_int nNormalPlanes, int isNode, flow_float** qL, flow_float** qP)
{
    if (ic1 < nCells) {
        const geom_int up = (meff >= 0.0f) ? ic0 : ic1;
        flow_float F = meff * qL[q][up] / max(ro[up], (flow_float)1.0e-30);
        if (q == qDiff && cdiff != nullptr && ip < nNormalPlanes)
            F -= cdiff[ip] * (qL[q][ic1]/max(ro[ic1], (flow_float)1.0e-30) - qL[q][ic0]/max(ro[ic0], (flow_float)1.0e-30));
        return F;
    }
    if (isNode != 0) {
        if (meff >= 0.0f) return meff * qL[q][ic0] / max(ro[ic0], (flow_float)1.0e-30);
        return meff * qP[q][ic0] / max(roN[ic0], (flow_float)1.0e-30);
    }
    const geom_int up = (meff >= 0.0f) ? ic0 : ic1;
    return meff * qL[q][up] / max(ro[up], (flow_float)1.0e-30);
}

// (0) 診断: HO の BDF 残差 r_H = res^sp − (V/Δt)(a q − b qP + c qPP) の Σ r² と Σ (V a/Δt q)² を**成分ごと** (out2[2q], out2[2q+1]; root のみ)。
__global__ void passive_fct_rh_norm_d(
    geom_int nCells, const geom_float* vol, flow_float dt, flow_float a, flow_float b, flow_float c, int nq,
    flow_float** res, flow_float** q, flow_float** qP, flow_float** qPP, const geom_int* root, double* out2)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (root != nullptr && root[ic] != ic) return;
    const flow_float V = (flow_float)vol[ic];
    for (int qq = 0; qq < nq; ++qq) {
        const double r = (double)res[qq][ic] - (double)(V/dt)*((double)a*q[qq][ic] - (double)b*qP[qq][ic] + (double)c*qPP[qq][ic]);
        const double sc = (double)(V*a/dt)*(double)q[qq][ic];
        if (r != 0.0) atomicAdd(&out2[2*qq], r*r);
        if (sc != 0.0) atomicAdd(&out2[2*qq+1], sc*sc);
    }
}

// (1) 低次作用素の対角 (移流: 全輸送面; 拡散: トレーサのみ内部面) と拡散面係数。
__global__ void passive_fct_lo_diag_d(
    geom_int nCells, geom_int nNormalPlanes, geom_int nNormalHaloPlanes, const geom_int* normal_halo_planes, const geom_int* plane_cells,
    const flow_float* ro, const flow_float* meffFace, int isNode,
    int haveDiff, const geom_float* ccx, const geom_float* ccy, const geom_float* ccz,
    const geom_float* fx, const geom_float* sx, const geom_float* sy, const geom_float* sz, const geom_float* ss,
    const flow_float* vis_lam, const flow_float* vis_turb, flow_float Sc, flow_float Sc_t,
    flow_float* diagAdv, flow_float* diagDiff, flow_float* cdiff)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0];
    const geom_int ic1 = plane_cells[2*ip+1];
    const flow_float mdot = meffFace[ip];
    if (ic1 < nCells && ic0 < nCells) {
        if (mdot >= 0.0f) atomicAdd(&diagAdv[ic0],  mdot / max(ro[ic0], (flow_float)1.0e-30));
        else              atomicAdd(&diagAdv[ic1], -mdot / max(ro[ic1], (flow_float)1.0e-30));
        if (haveDiff != 0 && ip < nNormalPlanes) {
            const flow_float c = passive_fct_diff_coef(ip, ic0, ic1, ccx, ccy, ccz, fx, sx, sy, sz, ss, ro, vis_lam, vis_turb, Sc, Sc_t);
            cdiff[ip] = c;
            atomicAdd(&diagDiff[ic0], c / max(ro[ic0], (flow_float)1.0e-30));
            atomicAdd(&diagDiff[ic1], c / max(ro[ic1], (flow_float)1.0e-30));
        }
    } else if (ic0 < nCells) {
        if (mdot > 0.0f) atomicAdd(&diagAdv[ic0], mdot / max(ro[ic0], (flow_float)1.0e-30));   // 流出のみ対角 (流入は許容な外部状態の RHS: nb)
    }
}

// (2) Jacobi sweep の近傍和。
__global__ void passive_fct_lo_nb_d(
    geom_int nCells, geom_int nNormalPlanes, geom_int nNormalHaloPlanes, const geom_int* normal_halo_planes, const geom_int* plane_cells,
    const flow_float* ro, const flow_float* roN, const flow_float* meffFace, int isNode,
    int nq, int qDiff, const flow_float* cdiff, flow_float** qL, flow_float** rophi, flow_float** qP, flow_float** nb,
    int bndMode)   // 0: 内部近傍 + 境界の定数 RHS (通常), 1: 内部近傍のみ (作用素 L の非対角), 2: 境界の定数 RHS のみ (plan-8 M3 の上限診断用)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0];
    const geom_int ic1 = plane_cells[2*ip+1];
    const flow_float mdot = meffFace[ip];
    if (ic1 < nCells && ic0 < nCells) {
        if (bndMode == 2) return;
        const geom_int up = (mdot >= 0.0f) ? ic0 : ic1, dn = (mdot >= 0.0f) ? ic1 : ic0;
        const flow_float w = fabsf(mdot) / max(ro[up], (flow_float)1.0e-30);
        for (int q = 0; q < nq; ++q) atomicAdd(&nb[q][dn], w * qL[q][up]);
        if (qDiff >= 0 && ip < nNormalPlanes && cdiff != nullptr) {
            const flow_float c = cdiff[ip];
            if (c != 0.0f) {
                atomicAdd(&nb[qDiff][ic0], (c / max(ro[ic1], (flow_float)1.0e-30)) * qL[qDiff][ic1]);
                atomicAdd(&nb[qDiff][ic1], (c / max(ro[ic0], (flow_float)1.0e-30)) * qL[qDiff][ic0]);
            }
        }
    } else if (ic0 < nCells && mdot < 0.0f) {
        if (bndMode == 1) return;
        if (isNode != 0) { const flow_float w = -mdot / max(roN[ic0], (flow_float)1.0e-30); for (int q = 0; q < nq; ++q) atomicAdd(&nb[q][ic0], w * qP[q][ic0]); }   // 許容な外部状態 φ^n_own
        else             { const flow_float w = -mdot / max(ro[ic1], (flow_float)1.0e-30);  for (int q = 0; q < nq; ++q) atomicAdd(&nb[q][ic0], w * rophi[q][ic1]); }
    }
}

// (3) rhs = (V/Δt) q^n + S^eff V,  S^eff V = invA·resS + cOverA·H^n (H が nullptr なら H^n = (V/Δt)(qP − qPP) [restart の代替: 全部局所])。
//   rhs2 += Σ rhs² (root のみ; 受入判定の分母)。
__global__ void passive_fct_lo_rhs_d(
    geom_int nCells, const geom_float* vol, flow_float dt, flow_float invA, flow_float cOverA,
    int nq, flow_float** qP, flow_float** qPP, flow_float** resS, flow_float** H, const geom_int* root, flow_float** rhs, double* rhs2)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float V = (flow_float)vol[ic];
    const bool count = (root == nullptr) || (root[ic] == ic);
    for (int q = 0; q < nq; ++q) {
        const flow_float Hn = (H != nullptr) ? H[q][ic] : (V/dt)*(qP[q][ic] - qPP[q][ic]);
        const flow_float r = (V/dt)*qP[q][ic] + invA*(resS != nullptr ? resS[q][ic] : 0.0f) + cOverA*Hn;
        rhs[q][ic] = r;
        if (rhs2 != nullptr && count && r != 0.0f) atomicAdd(&rhs2[q], (double)r*(double)r);   // 成分ごと (plan-7 M5)
    }
}

// (4) 初期値: q_L = q_H を物理限界 ([0, ρ] / ≥0) にクリップ。
__global__ void passive_fct_lo_init_d(geom_int nCells, const flow_float* ro, int nq, int nUnit, flow_float** qH, flow_float** qL)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    for (int q = 0; q < nq; ++q) {
        flow_float v = max(qH[q][ic], (flow_float)0.0);
        if (q < nUnit) v = min(v, ro[ic]);
        qL[q][ic] = v;
    }
}

// (5) Jacobi 更新: q_L = (rhs + nb)/(V/Δt + diagAdv [+ diagDiff]); ピン行は q_H。res2 += Σ r_L² (更新前の反復値の線形残差; root のみ)。
__global__ void passive_fct_lo_solve_d(
    geom_int nCells, const geom_float* vol, flow_float oneOverDt,
    int nq, int qDiff, const flow_float* diagAdv, const flow_float* diagDiff,
    flow_float** rhs, flow_float** nb, flow_float** qH, const flow_float* pin, const geom_int* root,
    flow_float** qL, double* res2)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const bool pinned = (pin != nullptr) && (pin[ic] == (flow_float)1.0);
    const bool count = (root == nullptr) || (root[ic] == ic);
    const flow_float base = oneOverDt*(flow_float)vol[ic] + diagAdv[ic];
    for (int q = 0; q < nq; ++q) {
        if (pinned) { qL[q][ic] = qH[q][ic]; continue; }
        const flow_float d = max(base + ((q == qDiff && diagDiff != nullptr) ? diagDiff[ic] : 0.0f), (flow_float)1.0e-30);
        const flow_float rl = rhs[q][ic] + nb[q][ic] - d*qL[q][ic];
        if (res2 != nullptr && count && rl != 0.0f) atomicAdd(&res2[q], (double)rl*(double)rl);   // 成分ごと (plan-7 M5)
        qL[q][ic] = (rhs[q][ic] + nb[q][ic]) / d;
    }
}

// (6) 生の反拡散流束 A^raw = F^eff_H − F_L(q_L) (全輸送面; F^eff_H = invA·F_H + cOverA·G^n, F_H = ṁ P_face − J_H [内部] / ṁ φ_H,own [node 境界]) と
//     基点の集計 base[ic] += s_f A^raw。Aface[ip*stride+q] = A^raw (全面; G が nullptr なら G^n = 0)。cell の ghost 面は両者 1 次で同じなので 0。
__global__ void passive_fct_raw_d(
    geom_int nCells, geom_int nNormalPlanes, geom_int nNormalHaloPlanes, const geom_int* normal_halo_planes, const geom_int* plane_cells,
    const flow_float* ro, const flow_float* roN, const flow_float* massflux, const flow_float* meffFace, int isNode, flow_float invA, flow_float cOverA,
    int nq, int stride, const flow_float* Pface, int qDiff, const flow_float* cdiff, const flow_float* G,
    flow_float** qH, flow_float** qL, flow_float** qP, flow_float* Aface, flow_float** base)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0];
    const geom_int ic1 = plane_cells[2*ip+1];
    const flow_float mdot = massflux[ip];
    const flow_float meff = meffFace[ip];
    const bool interior = (ic0 < nCells && ic1 < nCells);
    if (!interior && !(ic0 < nCells && isNode != 0)) { for (int q = 0; q < nq; ++q) Aface[(size_t)ip*stride + q] = 0.0f; return; }
    for (int q = 0; q < nq; ++q) {
        flow_float FH;
        if (interior) {
            FH = mdot * Pface[(size_t)ip*stride + q];
            if (q == qDiff && cdiff != nullptr && ip < nNormalPlanes)
                FH -= cdiff[ip] * (qH[q][ic1]/max(ro[ic1], (flow_float)1.0e-30) - qH[q][ic0]/max(ro[ic0], (flow_float)1.0e-30));
        } else {
            FH = mdot * qH[q][ic0] / max(ro[ic0], (flow_float)1.0e-30);
        }
        const flow_float FL = passive_fct_lowflux(nCells, ip, ic0, ic1, meff, ro, roN, q, qDiff, cdiff, nNormalPlanes, isNode, qL, qP);
        const flow_float A = invA*FH + (G != nullptr ? cOverA*G[(size_t)ip*stride + q] : 0.0f) - FL;
        Aface[(size_t)ip*stride + q] = A;
        if (A != 0.0f) { atomicAdd(&base[q][ic0], -A); if (interior) atomicAdd(&base[q][ic1], A); }
    }
}

// (7) 基点 q_B = q_H − (Δt/V) base (gather 済み)。物理限界の逸脱量を stats[q*8+4] (絶対量·V, root のみ) に記録し、限界用にはクリップ値を書く。
__global__ void passive_fct_base_d(
    geom_int nCells, const geom_float* vol, flow_float dt, const flow_float* ro, int nq, int nUnit,
    flow_float** qH, flow_float** base, const geom_int* root, flow_float** qB, double* stats)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const bool count = (root == nullptr) || (root[ic] == ic);
    const flow_float V = max((flow_float)vol[ic], (flow_float)1.0e-30);
    for (int q = 0; q < nq; ++q) {
        const flow_float v0 = qH[q][ic] - (dt/V)*base[q][ic];
        flow_float v = max(v0, (flow_float)0.0);
        if (q < nUnit) v = min(v, ro[ic]);
        if (v != v0 && stats != nullptr && count) atomicAdd(&stats[(size_t)q*8 + 4], (double)fabsf(v - v0)*(double)V);
        qB[q][ic] = v;
    }
}

// (8) 前制限と P± (全輸送面; 境界は ic0 側だけ)。A^pre = 0 if A^raw (φ_B(ic1) − φ_B(ic0)) < 0 (内部面のみ判定)。stats[q*8+2] += Σ|前制限で落とした A|。
__global__ void passive_fct_prelimit_d(
    geom_int nCells, geom_int nNormalHaloPlanes, const geom_int* normal_halo_planes, const geom_int* plane_cells, const flow_float* ro,
    int nq, int stride, const flow_float* Aface, int prelimit, flow_float** qB,
    flow_float* Apre, flow_float** Pp, flow_float** Pm, double* stats)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0];
    const geom_int ic1 = plane_cells[2*ip+1];
    if (ic0 >= nCells) return;
    const bool interior = (ic1 < nCells);
    for (int q = 0; q < nq; ++q) {
        flow_float A = Aface[(size_t)ip*stride + q];
        if (prelimit != 0 && interior && A != 0.0f) {
            const flow_float d = qB[q][ic1]/max(ro[ic1], (flow_float)1.0e-30) - qB[q][ic0]/max(ro[ic0], (flow_float)1.0e-30);
            if (A*d < 0.0f) { if (stats != nullptr) atomicAdd(&stats[(size_t)q*8 + 2], (double)fabsf(A)); A = 0.0f; }
        }
        Apre[(size_t)ip*stride + q] = A;
        if (A > 0.0f)      { atomicAdd(&Pm[q][ic0], A);  if (interior) atomicAdd(&Pp[q][ic1], A); }   // ic0 が出し
        else if (A < 0.0f) { atomicAdd(&Pp[q][ic0], -A); if (interior) atomicAdd(&Pm[q][ic1], -A); }  // ic0 が受け
    }
}

// (9) 局所極値 (自身 + 内部面隣接) of {φ^n = q^n/ρ^n, φ_B = q_B/ρ}。
__global__ void passive_fct_extrema_d(
    geom_int nCells, geom_int nNormalPlanes, const geom_int* plane_cells, const geom_int* cell_planes_index, const geom_int* cell_planes,
    const flow_float* ro, const flow_float* roN, int nq, flow_float** qP, flow_float** qB, flow_float** phimax, flow_float** phimin)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const geom_int st = cell_planes_index[ic], en = cell_planes_index[ic+1];
    for (int q = 0; q < nq; ++q) {
        const flow_float a = qP[q][ic]/max(roN[ic], (flow_float)1.0e-30), b = qB[q][ic]/max(ro[ic], (flow_float)1.0e-30);
        flow_float mx = max(a, b), mn = min(a, b);
        for (geom_int ilp = st; ilp < en; ++ilp) {
            const geom_int ip = cell_planes[ilp];
            if (ip >= nNormalPlanes) continue;
            const geom_int jc = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic;
            if (jc >= nCells) continue;
            const flow_float ja = qP[q][jc]/max(roN[jc], (flow_float)1.0e-30), jb = qB[q][jc]/max(ro[jc], (flow_float)1.0e-30);
            mx = max(mx, max(ja, jb)); mn = min(mn, min(ja, jb));
        }
        phimax[q][ic] = mx; phimin[q][ic] = mn;
    }
}

// (10) R± = min(1, Q±/P±), Q± = (ρ φ_max − q_B, q_B − ρ φ_min)·V/Δt。nUnit 本は [0,1] と交差、残りは下限 0。ピン行は R=1。
__global__ void passive_fct_ratio_d(
    geom_int nCells, const geom_float* vol, flow_float oneOverDt, const flow_float* ro, int nq, int nUnit,
    flow_float** qB, flow_float** phimax, flow_float** phimin, flow_float** Pp, flow_float** Pm,
    const flow_float* pin, flow_float** Rp, flow_float** Rm)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const bool pinned = (pin != nullptr) && (pin[ic] == (flow_float)1.0);
    const flow_float r = max(ro[ic], (flow_float)1.0e-30);
    const flow_float VoverDt = oneOverDt*(flow_float)vol[ic];
    for (int q = 0; q < nq; ++q) {
        if (pinned) { Rp[q][ic] = 1.0f; Rm[q][ic] = 1.0f; continue; }
        flow_float mx = phimax[q][ic], mn = max(phimin[q][ic], (flow_float)0.0);
        if (q < nUnit) mx = min(mx, (flow_float)1.0);
        const flow_float L = qB[q][ic];
        const flow_float Qp = (r*mx - L) * VoverDt;
        const flow_float Qm = (L - r*mn) * VoverDt;
        const flow_float pp = Pp[q][ic], pm = Pm[q][ic];
        flow_float rp = 1.0f, rm = 1.0f;
        if (pp > 0.0f) rp = (Qp <= 0.0f) ? 0.0f : min((flow_float)1.0, Qp/pp);
        if (pm > 0.0f) rm = (Qm <= 0.0f) ? 0.0f : min((flow_float)1.0, Qm/pm);
        Rp[q][ic] = rp; Rm[q][ic] = rm;
    }
}

// (11) α_f, 補正の集計, 実現流束 G^{n+1} (全輸送面): corr[ic0] −= (A^raw − α A^pre), corr[ic1] += (…)。境界面は外部側 R=1 (ノード側だけ制約)。
//   Gout[ip*stride+q] = F_L(q_L) + α A^pre。stats[q*8+0] += Σ|A^raw − α A^pre|, [1] += 作動面数, [5] += 境界面の実現流束 G (符号付き, 流出正), [6] += 境界面で落とした |A^raw − α A^pre|。forceAlpha0 は単体試験用。
__global__ void passive_fct_apply_d(
    geom_int nCells, geom_int nNormalPlanes, geom_int nNormalHaloPlanes, const geom_int* normal_halo_planes, const geom_int* plane_cells,
    const flow_float* ro, const flow_float* roN, const flow_float* meffFace, int isNode,
    int nq, int stride, const flow_float* Aface, const flow_float* Apre, flow_float** Rp, flow_float** Rm, int forceAlpha0,
    int qDiff, const flow_float* cdiff, flow_float** qL, flow_float** qP,
    flow_float** corr, flow_float* Gout, double* stats)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0];
    const geom_int ic1 = plane_cells[2*ip+1];
    if (ic0 >= nCells) return;
    const bool interior = (ic1 < nCells);
    const flow_float meff = meffFace[ip];
    for (int q = 0; q < nq; ++q) {
        const flow_float Ar = Aface[(size_t)ip*stride + q];
        const flow_float Ap = Apre[(size_t)ip*stride + q];
        flow_float alpha;
        if (Ap > 0.0f)      alpha = interior ? min(Rp[q][ic1], Rm[q][ic0]) : Rm[q][ic0];   // ic0 が出し
        else if (Ap < 0.0f) alpha = interior ? min(Rp[q][ic0], Rm[q][ic1]) : Rp[q][ic0];   // ic0 が受け
        else                alpha = 0.0f;
        if (forceAlpha0 != 0) alpha = 0.0f;
        const flow_float aA = alpha*Ap;
        const flow_float Gf = passive_fct_lowflux(nCells, ip, ic0, ic1, meff, ro, roN, q, qDiff, cdiff, nNormalPlanes, isNode, qL, qP) + aA;
        if (Gout != nullptr) Gout[(size_t)ip*stride + q] = Gf;
        if (!interior && stats != nullptr) { atomicAdd(&stats[(size_t)q*8 + 5], (double)Gf); atomicAdd(&stats[(size_t)q*8 + 6], (double)fabsf(Ar - aA)); }   // 境界: 実現流束 (符号付き, 流出正) と落とした量
        const flow_float rA = Ar - aA;
        if (rA != 0.0f) {
            atomicAdd(&corr[q][ic0], -rA);
            if (interior) atomicAdd(&corr[q][ic1], rA);
            if (stats != nullptr) { atomicAdd(&stats[(size_t)q*8 + 0], (double)fabsf(rA)); atomicAdd(&stats[(size_t)q*8 + 1], 1.0); }
        }
    }
}

// (12) 確定: q = q_H − (Δt/V) corr。ピン行は据え置き、落とした交換量 |corr|·Δt を stats[q*8+3] に記録 (境界収支)。
__global__ void passive_fct_commit_d(
    geom_int nCells, const geom_float* vol, flow_float dt, int nq,
    flow_float** corr, const flow_float* pin, const geom_int* root, flow_float** rophi, double* stats)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const bool pinned = (pin != nullptr) && (pin[ic] == (flow_float)1.0);
    const bool count = (root == nullptr) || (root[ic] == ic);
    for (int q = 0; q < nq; ++q) {
        const flow_float c = corr[q][ic];
        if (c == 0.0f) continue;
        if (pinned) { if (stats != nullptr && count) atomicAdd(&stats[(size_t)q*8 + 3], (double)fabsf(c)*(double)dt); continue; }
        rophi[q][ic] -= (dt / max((flow_float)vol[ic], (flow_float)1.0e-30)) * c;
    }
}

// (13) 履歴の局所残り: 面ループで divG[ic] += s_f G_f (全輸送面; 周期は呼び出し側で和 gather) → H = (V/Δt)(q_final − q^n) − divG。
__global__ void passive_fct_divG_d(
    geom_int nCells, geom_int nNormalHaloPlanes, const geom_int* normal_halo_planes, const geom_int* plane_cells,
    int nq, int stride, const flow_float* G, flow_float** divG)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0];
    const geom_int ic1 = plane_cells[2*ip+1];
    for (int q = 0; q < nq; ++q) {
        const flow_float g = G[(size_t)ip*stride + q];
        if (g == 0.0f) continue;
        if (ic0 < nCells) atomicAdd(&divG[q][ic0], -g);
        if (ic1 < nCells) atomicAdd(&divG[q][ic1],  g);
    }
}
//   H = H_src + H_rem: H_src = invA·S V + cOverA·H_src^n (物理ソースの流束形履歴), H_rem = H − H_src (残差・クランプ・ピン由来の局所残り; 収支の「非物理」項)。
//   budget[q*4+0] += Σ H_src·Δt·(root), [1] += Σ H_rem·Δt (符号付き), [2] += Σ|H_rem|·Δt, [3] += Σ (V/Δt)(q_final − q^n)·Δt (総増分)。
__global__ void passive_fct_hist_local_d(geom_int nCells, const geom_float* vol, flow_float dt, flow_float invA, flow_float cOverA, int nq,
                                         flow_float** qfinal, flow_float** qP, flow_float** divG, flow_float** resS, flow_float** Hsrc, flow_float** H,
                                         const geom_int* root, double* budget)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float V = (flow_float)vol[ic];
    const bool count = (root == nullptr) || (root[ic] == ic);
    for (int q = 0; q < nq; ++q) {
        const flow_float inc = (V/dt)*(qfinal[q][ic] - qP[q][ic]);
        const flow_float Hn = inc - divG[q][ic];
        const flow_float Hs = invA*(resS != nullptr ? resS[q][ic] : 0.0f) + cOverA*Hsrc[q][ic];
        H[q][ic] = Hn; Hsrc[q][ic] = Hs;
        if (budget != nullptr && count) {
            const double rem = (double)(Hn - Hs)*(double)dt;
            if (Hs != 0.0f) atomicAdd(&budget[(size_t)q*4 + 0], (double)Hs*(double)dt);
            if (rem != 0.0) { atomicAdd(&budget[(size_t)q*4 + 1], rem); atomicAdd(&budget[(size_t)q*4 + 2], fabs(rem)); }
            if (inc != 0.0f) atomicAdd(&budget[(size_t)q*4 + 3], (double)inc*(double)dt);
        }
    }
}

// (15) 密度整合の診断 (plan-7 M6): E_ρ = (V/Δt)(ρ^{n+1} − ρ^n) − Σ s ṁ^eff (root のみ; Σ E_ρ² と Σ ((V/Δt)ρ)²)、
//      トレーサ上限条件 Lρ − f ≥ 0 の逸脱 (負の量の Σ·Δt と個数)。Lρ = (V/Δt + diagAdv + diagDiff)ρ − nb(ρ) は呼び出し側が nb を ρ で組んで渡す。
__global__ void passive_fct_density_closure_d(geom_int nCells, const geom_float* vol, flow_float dt, const flow_float* ro, const flow_float* roN,
                                              const flow_float* divMeff, const geom_int* root, double* out2)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (root != nullptr && root[ic] != ic) return;
    const double V = (double)vol[ic];
    const double e = (V/dt)*((double)ro[ic] - (double)roN[ic]) - (double)divMeff[ic];
    const double sc = (V/dt)*(double)ro[ic];
    if (e != 0.0) atomicAdd(&out2[0], e*e);
    if (sc != 0.0) atomicAdd(&out2[1], sc*sc);
}
//   Lρ = D ρ − nbInt(ρ) (内部近傍のみ), f_full = rhs + nbBnd (境界の定数 RHS |ṁ^eff| φ^n_own; トレーサ q の実値)。逸脱 (Lρ − f_full < 0) の量 (·Δt は呼び出し側) と個数。
__global__ void passive_fct_upper_margin_d(geom_int nCells, const geom_float* vol, flow_float oneOverDt, const flow_float* ro,
                                           const flow_float* diagAdv, const flow_float* diagDiff, int q, const flow_float* nbIntRho, const flow_float* nbBndQ, flow_float** rhs,
                                           const flow_float* pin, const geom_int* root, double* out)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (root != nullptr && root[ic] != ic) return;
    if (pin != nullptr && pin[ic] == (flow_float)1.0) return;
    const double Lrho = ((double)oneOverDt*(double)vol[ic] + (double)diagAdv[ic] + (diagDiff != nullptr ? (double)diagDiff[ic] : 0.0))*(double)ro[ic] - (double)nbIntRho[ic];
    const double m = Lrho - ((double)rhs[q][ic] + (nbBndQ != nullptr ? (double)nbBndQ[ic] : 0.0));
    if (m < 0.0) { atomicAdd(&out[0], -m); atomicAdd(&out[1], 1.0); }
}
__global__ void passive_fct_div_meff_d(geom_int nCells, geom_int nNormalHaloPlanes, const geom_int* normal_halo_planes, const geom_int* plane_cells,
                                       const flow_float* meff, flow_float* divMeff)
{
    const geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2*ip+0], ic1 = plane_cells[2*ip+1];
    const flow_float m = meff[ip];
    if (ic0 < nCells) atomicAdd(&divMeff[ic0], -m);
    if (ic1 < nCells) atomicAdd(&divMeff[ic1],  m);
}

// (14) 有効質量流束 ṁ^eff,{n+1} = invA·ṁ^{n+1} + cOverA·ṁ^eff,n (流れの BDF2 連続式の流束形; BDF1 は ṁ そのもの)。
__global__ void passive_fct_meff_d(geom_int nPlanes, flow_float invA, flow_float cOverA, const flow_float* massflux, flow_float* meff)
{
    const geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;
    if (ip < nPlanes) meff[ip] = invA*massflux[ip] + cOverA*meff[ip];
}
