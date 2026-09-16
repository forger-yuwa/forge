#pragma once
// 凝縮モーメントの更新クランプ kernel (condLimiterMode 1)。condensationTransport_d.cu と単体試験 tests/unit/test_cond_limiter_steady.cu が include する。
#include "flowFormat.hpp"
#include "cuda_forge/condensationProperties_d.cuh"   // CondPropOpts / condProps_make / cond_latent

// 更新クランプ (condLimiterMode 1; plans/active/condensation-source-limiter-steady.md §4.2-4)。
//   定常 point-implicit 更新 (timeIntegration 11) の 4 モーメントについて、floor 前の候補増分
//     δ_k = (res_k Δτ/V) / (1 + Δτ (src_jac_k + transport_diag_k/V))
//   を取り出し、Δg = δ_g/ρ (更新済みの流れ ρ を固定したモーメント修正量) から
//     θ_u = min(1, dg_max/|Δg|, dT_max/|ΔT|, avail/Δg [Δg>0], (1−λ_min³) g_old/|Δg| [Δg<0])
//   を作って **4 本の増分を同率で縮めてから** floor (≥0) を掛けて確定する。θ_u > 0 (更新を止める穴を作らない:
//   avail≤0 は残差側で S=0 なので候補増分は輸送分のみ)。収束時は δ→0 で θ_u→1・無作用 = 固定点は残差だけで決まる。
//   潜熱 ΔT は二相 EOS と同じ有効比熱 c_v,eff = c_v + g (R_w − dL/dT) で評価する。
//   診断: diagLim = θ_u、diagCorrG = floor による ρg の補正量 [質量分率]、diagCorrQ = Q0..Q2 の相対補正 (収束時 0 を確認する)。
// 実体 (__device__)。既存カーネル cond_moment_update_limited_d は (relax=1, dtScale=1, applyFloor=1, dq=nullptr) で呼び、
// 式・評価順は不変 (1.0 の乗算は厳密なのでビット不変)。受動種経路 (passiveScalarScheme 1) は
//   relax   : 候補増分の緩和 (passiveImplicitRelax; θ_u の評価前に掛ける = 緩和後の増分に対してクランプ)
//   dtScale : dt_local の倍率 (scalarCflMax)
//   applyFloor 0: floor (≥0) を掛けず候補値をそのまま書く (後段 passive_bounds_d が floor と補正収支を担う。診断 diagCorrG/Q は不変)
//   dq_*    : 非 nullptr なら候補増分を point-implicit で組まず、scalar-DPLUR sweep が作った増分 (緩和済み) をそのまま使う
//             (passiveImplicitCoupling 1)。θ_u クランプ・floor・診断は同じ。
__device__ __forceinline__ void cond_moment_update_limited_body(
    geom_int ic, flow_float* dt_local, geom_float* vol, flow_float* ro,
    flow_float* roY_w, double Yw_const, flow_float* T, flow_float* cp_cell, flow_float* Rmix_cell, flow_float cp_cpg, flow_float gamma_cpg,
    int condModel, const CondPropOpts& opts, double dg_max, double dT_max, double lam_min,
    flow_float* N_g, flow_float* N_Q2, flow_float* N_Q1, flow_float* N_Q0,
    flow_float* res_g, flow_float* res_Q2, flow_float* res_Q1, flow_float* res_Q0,
    flow_float* sj_g, flow_float* sj_Q2, flow_float* sj_Q1, flow_float* sj_Q0,
    flow_float* td_g, flow_float* td_Q2, flow_float* td_Q1, flow_float* td_Q0,
    flow_float* out_g, flow_float* out_Q2, flow_float* out_Q1, flow_float* out_Q0,
    flow_float* diagLim, flow_float* diagCorrG, flow_float* diagCorrQ,
    double relax, flow_float dtScale, int applyFloor,
    flow_float* dq_g, flow_float* dq_Q2, flow_float* dq_Q1, flow_float* dq_Q0,
    // codex result M5 (受動種経路のみ; boundByTheta 0 で従来と同一): 4 モーメントの非負を **共通 θ** で保証する増分縮小
    // θ_neg = min_k N_k/|d_k| (d_k<0 かつ N_k+d_k<0)。縮小量 (θ_u − θ)|d_k| を limCorr_k (セル累積) と limStats[8k] (Σ·V, root のみ; 受動種収支の stride 8) に記録。
    int boundByTheta, flow_float* limCorr_g, flow_float* limCorr_Q2, flow_float* limCorr_Q1, flow_float* limCorr_Q0,
    double* limStats, const geom_int* root)
{
    const double dt = (double)(dt_local[ic] * dtScale);
    const double v  = (double)vol[ic];
    const double rod = (double)ro[ic];
    // 候補増分 (floor 前) — runge_kutta_exp_scalar_d (coef 1/0/1) と同じ式
    auto cand = [&](flow_float* res, flow_float* sj, flow_float* td) -> double {
        const double fac = 1.0 + dt*((double)sj[ic] + (double)td[ic]/v);
        return ((double)res[ic]*dt/v)/fac;
    };
    double d_g, d_Q2, d_Q1, d_Q0;
    if (dq_g != nullptr) {
        d_g = (double)dq_g[ic]; d_Q2 = (double)dq_Q2[ic]; d_Q1 = (double)dq_Q1[ic]; d_Q0 = (double)dq_Q0[ic];
    } else {
        d_g  = cand(res_g,  sj_g,  td_g)  * relax;
        d_Q2 = cand(res_Q2, sj_Q2, td_Q2) * relax;
        d_Q1 = cand(res_Q1, sj_Q1, td_Q1) * relax;
        d_Q0 = cand(res_Q0, sj_Q0, td_Q0) * relax;
    }

    double theta = 1.0;
    if (rod > 1.0e-20 && dt > 0.0) {
        const double g_old = fmax((double)N_g[ic]/rod, 0.0);
        const double dg = d_g/rod;
        if (dg != 0.0) {
            const CondSpeciesProps cprops = condProps_make(condModel, opts);
            const double Td = (double)T[ic];
            const double cpg = (cp_cell != nullptr) ? (double)cp_cell[ic] : (double)cp_cpg;
            const double Rg  = (Rmix_cell != nullptr) ? (double)Rmix_cell[ic] : ((double)gamma_cpg-1.0)*(double)cp_cpg/(double)gamma_cpg;
            const double cvg = fmax(cpg - Rg, 1.0e-3);
            const double L   = cond_latent(cprops, Td);
            const double dL  = (cond_latent(cprops, Td + 0.1) - cond_latent(cprops, Td - 0.1))/0.2;
            const double cveff = fmax(cvg + g_old*(cprops.R - dL), 1.0e-2*cvg);
            const double adg = fabs(dg);
            const double adT = adg*L/cveff;
            if (adg > dg_max) theta = fmin(theta, dg_max/adg);
            if (adT > dT_max) theta = fmin(theta, dT_max/adT);
            if (dg > 0.0) {
                const double Yw = (roY_w != nullptr) ? (double)roY_w[ic]/rod : ((Yw_const > 0.0) ? Yw_const : 1.0);
                const double avail = (roY_w != nullptr || Yw_const > 0.0) ? (Yw - g_old) : (0.99 - g_old);
                if (avail > 0.0 && dg > avail) theta = fmin(theta, avail/dg);
            } else if (g_old > 0.0) {
                const double gmaxdrop = (1.0 - lam_min*lam_min*lam_min)*g_old;   // 蒸発: 半径半減/step 相当
                if (adg > gmaxdrop) theta = fmin(theta, gmaxdrop/adg);
            }
            if (!(theta > 1.0e-12)) theta = 1.0e-12;   // 常に正 (停止穴なし)
        }
    }
    if (boundByTheta != 0) {
        const double theta_u = theta;
        const double Nk[4] = {(double)N_g[ic], (double)N_Q2[ic], (double)N_Q1[ic], (double)N_Q0[ic]};
        const double dk[4] = {d_g, d_Q2, d_Q1, d_Q0};
        for (int k = 0; k < 4; ++k) {
            if (dk[k] < 0.0 && Nk[k] + theta*dk[k] < 0.0) {
                const double th = (Nk[k] > 0.0) ? Nk[k] / (-dk[k]) : 0.0;
                if (th < theta) theta = th;
            }
        }
        if (theta < theta_u) {
            const bool count = (root == nullptr) || (root[ic] == ic);
            flow_float* lc[4] = {limCorr_g, limCorr_Q2, limCorr_Q1, limCorr_Q0};
            for (int k = 0; k < 4; ++k) {
                const double amt = (theta_u - theta)*fabs(dk[k]);
                if (lc[k] != nullptr && count) lc[k][ic] += (flow_float)amt;
                if (limStats != nullptr && count) {
                    atomicAdd(&limStats[(size_t)k*8], amt*v);                              // 絶対量 (受動種 q=g..Q0 の収支スロット, stride 8)
                    atomicAdd(&limStats[(size_t)k*8 + 1], 1.0);                            // 作動セル数
                    atomicAdd(&limStats[(size_t)k*8 + 2], (theta - theta_u)*dk[k]*v);       // 符号付き (確定 − 候補)
                }
            }
        }
    }
    const double ng  = (double)N_g[ic]  + theta*d_g;
    const double nQ2 = (double)N_Q2[ic] + theta*d_Q2;
    const double nQ1 = (double)N_Q1[ic] + theta*d_Q1;
    const double nQ0 = (double)N_Q0[ic] + theta*d_Q0;
    if (applyFloor != 0) {
        out_g[ic]  = (flow_float)fmax(ng,  0.0);
        out_Q2[ic] = (flow_float)fmax(nQ2, 0.0);
        out_Q1[ic] = (flow_float)fmax(nQ1, 0.0);
        out_Q0[ic] = (flow_float)fmax(nQ0, 0.0);
    } else {
        out_g[ic]  = (flow_float)ng;
        out_Q2[ic] = (flow_float)nQ2;
        out_Q1[ic] = (flow_float)nQ1;
        out_Q0[ic] = (flow_float)nQ0;
    }
    diagLim[ic]  = (flow_float)theta;
    // 補正量の記録 (このステップの全補正の起点なのでリセット): G = floor による |Δρg|/ρ [質量分率], Q = Q0..Q2 の最大相対補正。
    // 後段の実現可能性クランプ (g≤Y_w / 0.99ρ, 負値, 液滴消滅) は cond_realizability_clamp_{,f_}d が同じ配列へ累積する (codex result M2)。
    diagCorrG[ic] = (flow_float)((ng < 0.0 && rod > 1.0e-20) ? (-ng/rod) : 0.0);
    double rq = 0.0;
    if (nQ2 < 0.0) rq = fmax(rq, 1.0); if (nQ1 < 0.0) rq = fmax(rq, 1.0); if (nQ0 < 0.0) rq = fmax(rq, 1.0);   // 負値 floor は 100 % 補正
    diagCorrQ[ic] = (flow_float)rq;
}

__global__ void cond_moment_update_limited_d(
    geom_int nCells, flow_float* dt_local, geom_float* vol, flow_float* ro,
    flow_float* roY_w, double Yw_const, flow_float* T, flow_float* cp_cell, flow_float* Rmix_cell, flow_float cp_cpg, flow_float gamma_cpg,
    int condModel, CondPropOpts opts, double dg_max, double dT_max, double lam_min,
    flow_float* N_g, flow_float* N_Q2, flow_float* N_Q1, flow_float* N_Q0,
    flow_float* res_g, flow_float* res_Q2, flow_float* res_Q1, flow_float* res_Q0,
    flow_float* sj_g, flow_float* sj_Q2, flow_float* sj_Q1, flow_float* sj_Q0,
    flow_float* td_g, flow_float* td_Q2, flow_float* td_Q1, flow_float* td_Q0,
    flow_float* out_g, flow_float* out_Q2, flow_float* out_Q1, flow_float* out_Q0,
    flow_float* diagLim, flow_float* diagCorrG, flow_float* diagCorrQ)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    cond_moment_update_limited_body(ic, dt_local, vol, ro, roY_w, Yw_const, T, cp_cell, Rmix_cell, cp_cpg, gamma_cpg,
        condModel, opts, dg_max, dT_max, lam_min, N_g, N_Q2, N_Q1, N_Q0, res_g, res_Q2, res_Q1, res_Q0,
        sj_g, sj_Q2, sj_Q1, sj_Q0, td_g, td_Q2, td_Q1, td_Q0, out_g, out_Q2, out_Q1, out_Q0, diagLim, diagCorrG, diagCorrQ,
        1.0, (flow_float)1.0, 1, nullptr, nullptr, nullptr, nullptr,
        0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
}

// 受動種経路 (passiveScalarScheme 1) 用: 緩和・dt 倍率・floor 省略・DPLUR 増分入力を持つ変種 (本文は同じ)。
__global__ void cond_moment_update_limited_passive_d(
    geom_int nCells, flow_float* dt_local, geom_float* vol, flow_float* ro,
    flow_float* roY_w, double Yw_const, flow_float* T, flow_float* cp_cell, flow_float* Rmix_cell, flow_float cp_cpg, flow_float gamma_cpg,
    int condModel, CondPropOpts opts, double dg_max, double dT_max, double lam_min,
    flow_float* N_g, flow_float* N_Q2, flow_float* N_Q1, flow_float* N_Q0,
    flow_float* res_g, flow_float* res_Q2, flow_float* res_Q1, flow_float* res_Q0,
    flow_float* sj_g, flow_float* sj_Q2, flow_float* sj_Q1, flow_float* sj_Q0,
    flow_float* td_g, flow_float* td_Q2, flow_float* td_Q1, flow_float* td_Q0,
    flow_float* out_g, flow_float* out_Q2, flow_float* out_Q1, flow_float* out_Q0,
    flow_float* diagLim, flow_float* diagCorrG, flow_float* diagCorrQ,
    double relax, flow_float dtScale, int applyFloor,
    flow_float* dq_g, flow_float* dq_Q2, flow_float* dq_Q1, flow_float* dq_Q0,
    int boundByTheta, flow_float* limCorr_g, flow_float* limCorr_Q2, flow_float* limCorr_Q1, flow_float* limCorr_Q0,
    double* limStats, const geom_int* root)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    cond_moment_update_limited_body(ic, dt_local, vol, ro, roY_w, Yw_const, T, cp_cell, Rmix_cell, cp_cpg, gamma_cpg,
        condModel, opts, dg_max, dT_max, lam_min, N_g, N_Q2, N_Q1, N_Q0, res_g, res_Q2, res_Q1, res_Q0,
        sj_g, sj_Q2, sj_Q1, sj_Q0, td_g, td_Q2, td_Q1, td_Q0, out_g, out_Q2, out_Q1, out_Q0, diagLim, diagCorrG, diagCorrQ,
        relax, dtScale, applyFloor, dq_g, dq_Q2, dq_Q1, dq_Q0,
        boundByTheta, limCorr_g, limCorr_Q2, limCorr_Q1, limCorr_Q0, limStats, root);
}
