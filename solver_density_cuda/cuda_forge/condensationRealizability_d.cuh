#pragma once
// 凝縮モーメントの実現可能性クランプ (double 実体)。condensationTransport_d.cu と単体試験 tests/unit/test_cond_limiter_steady.cu が include する
// (float 実体 cond_realizability_clamp_f_d は物性表を使うので condensationTransport_d.cu 内)。
#include "flowFormat.hpp"
#include "cuda_forge/condensationProperties_d.cuh"
#include "cuda_forge/condensationSource_d.cuh"   // COND_PI
#include "cuda_forge/condensationEOS_d.cuh"      // cond_clamp_vapor_pressure
#include "cuda_forge/condensationSourceF_d.cuh"  // float 実体の表評価 (cond_tab_*)

// 実現可能性クランプ: 液相保存量 rog と モーメント roQ0/1/2 を物理範囲に戻す。
//   0 ≤ rog ≤ roY_w (carrier: 利用可能な総水量) または ≤ 0.99 ρ (pure)。roQn ≥ 0。
// θ 律速は瞬間 Sg を抑えるが、陰解法の実効 Δt と dt_local の差で僅かに過凝縮しうるため
// 毎ステップここで硬クランプして二相 EOS・次ステップ source に渡す値を保証する。
// 蒸発 (condEvaporation=1) では液滴消滅の硬クランプも行う (plans/accepted/condensation-evaporation.md §4.4):
//   S=p_v/p_sat<=1 かつ液相あり かつ (Q0=0 [不整合] または r30<2*rmin) かつ g<=g_rm (潜熱飛びが小さい)
//   → rog,roQ0,roQ1,roQ2 を 0。source kernel の λ=0 は陰的緩和で厳密 0 に届かないためここで確定する。
//   T,P は前ステップ値 (dependentVariables 前) で判定する。核生成域 (S>1) には触れない。
__global__ void cond_realizability_clamp_d(
    geom_int nCells,
    flow_float* ro, flow_float* roY_w,   // roY_w: carrier の総水保存量 (pure では nullptr)
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    int evap, int condModel, double Rw, double rmin, double g_rm,
    flow_float* T, flow_float* P, CondPropOpts opts,
    flow_float* diagCorrG, flow_float* diagCorrQ)   // 補正量の記録 (nullptr 可): |Δρg|/ρ [質量分率] の累積, Q の最大相対補正 (codex result M2)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    // 実現可能性 g<=Y_w: TP carrier は roY_w (輸送), CPG carrier (空気) は定数 opts.Yw, pure は 0.99
    const flow_float gmax = (roY_w != nullptr) ? roY_w[ic] : ((opts.Yw > 0.0) ? (flow_float)opts.Yw*ro[ic] : (flow_float)0.99*ro[ic]);
    const flow_float r_in = rog[ic], q0_in = roQ0[ic], q1_in = roQ1[ic], q2_in = roQ2[ic];
    flow_float r = r_in;
    if (r < (flow_float)0.0) r = (flow_float)0.0;
    if (r > gmax)            r = gmax;
    rog[ic] = r;
    if (roQ0[ic] < (flow_float)0.0) roQ0[ic] = (flow_float)0.0;
    if (roQ1[ic] < (flow_float)0.0) roQ1[ic] = (flow_float)0.0;
    if (roQ2[ic] < (flow_float)0.0) roQ2[ic] = (flow_float)0.0;
    auto record = [&]() {
        if (diagCorrG == nullptr) return;
        const double rod0 = (double)ro[ic] > 1.0e-20 ? (double)ro[ic] : 1.0e-20;
        diagCorrG[ic] += (flow_float)(fabs((double)rog[ic] - (double)r_in)/rod0);
        double rq = 0.0;
        const flow_float qin[3] = {q0_in, q1_in, q2_in}; const flow_float qout[3] = {roQ0[ic], roQ1[ic], roQ2[ic]};
        for (int k = 0; k < 3; ++k) { const double den = fabs((double)qin[k]) > 1.0e-30 ? fabs((double)qin[k]) : 1.0e-30;
            const double rel = fabs((double)qout[k] - (double)qin[k])/den; if (qout[k] != qin[k] && rel > rq) rq = rel; }
        if ((flow_float)rq > diagCorrQ[ic]) diagCorrQ[ic] = (flow_float)rq;
    };

    if (!evap) { record(); return; }
    const double rod = (double)ro[ic];
    if (rod <= 1.0e-20) { record(); return; }
    const double g  = (double)r/rod;
    if (g > g_rm) { record(); return; }
    // g==0 で Q0/Q1/Q2 だけ残る「モーメント塵」(float アンダーフロー) も S<=1 では掃除する
    // (核生成域 S>1 は触らない: 新核の g がアンダーフローしていても Q0 は生かす)。
    const bool dust = (r <= (flow_float)0.0) &&
        (roQ0[ic] > (flow_float)0.0 || roQ1[ic] > (flow_float)0.0 || roQ2[ic] > (flow_float)0.0);
    if (r <= (flow_float)0.0 && !dust) { record(); return; }
    const CondSpeciesProps cprops = condProps_make(condModel, opts);
    const double Td = (double)T[ic];
    // 蒸気分圧は source kernel (cond_vapor_state) と同じ定義: TP carrier=ρ(Y_w−g)R_wT, CPG carrier=ρ(Y_w,const−g)R_wT, pure=全圧 (codex 2026-09-13 M2)
    const double pv = cond_clamp_vapor_pressure(rod, g, (roY_w != nullptr) ? (double)roY_w[ic]/rod : -1.0, opts.Yw, Rw, Td, (double)P[ic]);
    if (pv > cond_psat(cprops, Td)) { record(); return; }   // 過飽和: 消滅させない
    const double q0 = (double)roQ0[ic];
    bool remove = dust || (q0 <= 1.0e-30);
    if (!remove) {
        const double rho_l = cond_rho_cond(cprops, Td);
        const double r30 = cbrt(g/((4.0/3.0)*COND_PI*rho_l*q0/rod));
        remove = (r30 < 2.0*rmin);
    }
    if (remove) {
        rog[ic] = (flow_float)0.0; roQ0[ic] = (flow_float)0.0;
        roQ1[ic] = (flow_float)0.0; roQ2[ic] = (flow_float)0.0;
    }
    record();
}

// float 実体 (condFloat=1): 判定は閾値比較 (p_v<=p_sat, r30<2 r_min) なので ULP 差で消滅 step が 1 つずれ得る (plan §4.2-4)。
__global__ void cond_realizability_clamp_f_d(
    geom_int nCells,
    flow_float* ro, flow_float* roY_w,
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    int evap, float Rw, float rmin, float g_rm, float Yw_const,
    flow_float* T, flow_float* P, CondTablesF tb, CondSpeciesProps cpd,
    flow_float* diagCorrG, flow_float* diagCorrQ)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float gmax = (roY_w != nullptr) ? roY_w[ic] : ((Yw_const > 0.0f) ? Yw_const*ro[ic] : 0.99f*ro[ic]);
    const flow_float r_in = rog[ic], q0_in = roQ0[ic], q1_in = roQ1[ic], q2_in = roQ2[ic];
    flow_float r = r_in;
    if (r < 0.0f) r = 0.0f;
    if (r > gmax) r = gmax;
    rog[ic] = r;
    if (roQ0[ic] < 0.0f) roQ0[ic] = 0.0f;
    if (roQ1[ic] < 0.0f) roQ1[ic] = 0.0f;
    if (roQ2[ic] < 0.0f) roQ2[ic] = 0.0f;
    auto record = [&]() {
        if (diagCorrG == nullptr) return;
        const float rod0 = ro[ic] > 1.0e-20f ? ro[ic] : 1.0e-20f;
        diagCorrG[ic] += fabsf(rog[ic] - r_in)/rod0;
        float rq = 0.0f;
        const flow_float qin[3] = {q0_in, q1_in, q2_in}; const flow_float qout[3] = {roQ0[ic], roQ1[ic], roQ2[ic]};
        for (int k = 0; k < 3; ++k) { const float den = fabsf(qin[k]) > 1.0e-30f ? fabsf(qin[k]) : 1.0e-30f;
            const float rel = fabsf(qout[k] - qin[k])/den; if (qout[k] != qin[k] && rel > rq) rq = rel; }
        if (rq > diagCorrQ[ic]) diagCorrQ[ic] = rq;
    };
    if (!evap) { record(); return; }
    const float rod = ro[ic];
    if (rod <= 1.0e-20f) { record(); return; }
    const float g = r/rod;
    if (g > g_rm) { record(); return; }
    const bool dust = (r <= 0.0f) && (roQ0[ic] > 0.0f || roQ1[ic] > 0.0f || roQ2[ic] > 0.0f);
    if (r <= 0.0f && !dust) { record(); return; }
    const float Td = T[ic];
    // 蒸気分圧 (source kernel と同じ定義): TP carrier=ρ(Y_w−g)R_wT, CPG carrier=ρ(Y_w,const−g)R_wT, pure=全圧
    float pv;
    if (roY_w != nullptr)      { float yv = roY_w[ic]/rod - g; if (yv < 0.0f) yv = 0.0f; pv = rod*yv*Rw*Td; }
    else if (Yw_const > 0.0f)  { float yv = Yw_const - g;     if (yv < 0.0f) yv = 0.0f; pv = rod*yv*Rw*Td; }
    else                         pv = P[ic];
    const bool inTab = cond_tab_wet_ok(tb, Td);   // 表範囲外は旧 double 関数 (plan §5.1 #8)
    const float lnps = inTab ? cond_tab_lnpsat_f(tb, Td) : (float)log(cond_psat(cpd, (double)Td) > 1.0e-300 ? cond_psat(cpd, (double)Td) : 1.0e-300);
    if (pv > 0.0f && logf(pv) > lnps) { record(); return; }   // 過飽和: 消滅させない
    const float q0 = roQ0[ic];
    bool remove = dust || (q0 <= 1.0e-30f);
    if (!remove) {
        const float rho_l = inTab ? cond_tab_rhol_f(tb, Td) : (float)cond_rho_cond(cpd, (double)Td);
        const float r30 = cbrtf(g/((4.0f/3.0f)*COND_PI_F*rho_l*q0/rod));
        remove = (r30 < 2.0f*rmin);
    }
    if (remove) { rog[ic] = 0.0f; roQ0[ic] = 0.0f; roQ1[ic] = 0.0f; roQ2[ic] = 0.0f; }
    record();
}

