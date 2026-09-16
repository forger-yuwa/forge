#pragma once
// 凝縮モーメントの実現可能性クランプ (double 実体)。condensationTransport_d.cu と単体試験 tests/unit/test_cond_limiter_steady.cu が include する
// (float 実体 cond_realizability_clamp_f_d は物性表を使うので condensationTransport_d.cu 内)。
#include "flowFormat.hpp"
#include "cuda_forge/condensationProperties_d.cuh"
#include "cuda_forge/condensationSource_d.cuh"   // COND_PI
#include "cuda_forge/condensationEOS_d.cuh"      // cond_clamp_vapor_pressure
#include "cuda_forge/condensationSourceF_d.cuh"  // float 実体の表評価 (cond_tab_*)


// モーメント実現可能性の射影 (plan §4.7 v5, codex plan-5 M3/M4, plan-6 M3): 無次元 x = Q1/(Q0 r), y = Q2/(Q0 r²), r = (Q3/Q0)^{1/3}。

// 半径 r=(Q3/Q0)^{1/3} と無次元モーメント x=Q1/(Q0 r), y=Q2/(Q0 r²) を**指数を分離して**作る (codex result-4 M1)。
//   ρg が極小だと q3 = ρg/(4/3 π ρ_l) が underflow して r=0 になり、素直に割ると x,y が ∞/NaN になって
//   「修復できないので触らない」= 実現不能状態が残る。log 空間で組めば r が表現できない領域でも x,y は有限に出る。
// 戻り値 false: Q1, Q2 がともに 0 以下 (= 液滴なし; 呼び出し側で退化として扱う)。
__host__ __device__ inline bool cond_moment_xy(double q3, double q0, double q1, double q2, double& x, double& y)
{
    if (!(q3 > 0.0) || !(q0 > 0.0)) return false;
    const double lr = (log(q3) - log(q0))/3.0;                 // ln r
    x = (q1 > 0.0) ? exp(log(q1) - log(q0) - lr)       : 0.0;
    y = (q2 > 0.0) ? exp(log(q2) - log(q0) - 2.0*lr)   : 0.0;
    if (!(x == x) || !(y == y)) return false;
    if (x > 1.0e300) x = 1.0e300;
    if (y > 1.0e300) y = 1.0e300;
    return true;
}


// 射影後の (x, y) を保存量 Q1, Q2 へ書き戻す (codex result-5 M2)。q3/q0 が非正規化数に丸められると
// rr = cbrt(q3/q0) は桁が落ちて (x, y) を復元できない (射影直後なのに x が 14 % も 1 を超える例がある) ので、
// 直接の掛け算は q3/q0 が正規化数のときだけ使い、それ以外は対数で戻す。
__host__ __device__ inline void cond_moment_writeback(double q3, double q0, double x, double y, double& q1_out, double& q2_out)
{
    const double ratio = q3/q0;
    const double rr = cbrt(ratio);
    if (ratio >= 2.2250738585072014e-308 && rr > 0.0 && rr < 1.0e300) {
        q1_out = q0*rr*x; q2_out = q0*rr*rr*y;
    } else {
        const double lq0 = log(q0), lr = (log(q3) - lq0)/3.0;
        q1_out = (x > 0.0) ? exp(log(x) + lq0 + lr)       : 0.0;
        q2_out = (y > 0.0) ? exp(log(y) + lq0 + 2.0*lr)   : 0.0;
    }
}

// 許容領域 A = {0 ≤ x ≤ 1, x² ≤ y ≤ √x} (Hankel H1, H2 ⪰ 0)。判定は相対許容 eps (境界上の整合状態 [単分散 (1,1)] は触らない)。
//   特異不整合 (x または y が厳密 0 / アンダーフロー なのに Q3 > 0) は (Q0, g) 保存の単分散 (1,1) へ再初期化 (閾値なし: 小さい正の x は領域内部でもあり得る)。
//   それ以外の違反は A への最近点 (ユークリッド距離; x, y 両方を動かす) — 境界 y=√x / y=x² / 角 (1,1) の候補から最小距離を選ぶ (連続)。
// 戻り値: 0 = 変更なし, 1 = 最近点射影, 2 = 退化の単分散再初期化。
__host__ __device__ inline int cond_realizability_project(double& x, double& y, double eps)
{
    if (!(x == x) || !(y == y) || fabs(x) > 1.0e300 || fabs(y) > 1.0e300) return 0;   // 非有限 (r=0 由来の ∞ など) は触らない (2026-09-17 の NaN 退行の再発防止)
    // 特異不整合 (plan-7 M3): Q3>0 なのに Q1 または Q2 が厳密 0 / アンダーフロー (x, y ≤ 1e-30) → 非負半径分布では Q1=0 ⇒ Q3=0 なので実現不能。
    // 不等式検査より先に扱い、(Q0, g) 保存の単分散 (1,1) へ再初期化 (明示した修復方針; 数を別に log)。
    if (x <= 1.0e-30 || y <= 1.0e-30) { x = 1.0; y = 1.0; return 2; }
    const bool viol = (x > 1.0 + eps) || (y < x*x*(1.0 - eps)) || (y*y > x*(1.0 + eps));
    if (!viol) return 0;
    // 最近点の候補 (plan-8 M4: 小さい正の状態でも連続・最小): 単純クランプ候補 (x, clamp(y)), (clamp(x), y) を上界として持ち、
    // 境界 y=√t / y=t² 上の最近点を局所スケールの区間 [0, hi] で相対停止の黄金分割で探し、最小距離の候補を採る。
    const double x0 = x, y0 = y;
    // 初期候補は**必ず許容領域内の点** (単純クランプ) にする (codex result-3 の追跡で見つけた欠陥: x, y が極端に大きいと
    // すべての候補で距離が inf に overflow し、どれも採用されずに違反状態のまま「射影した」と返っていた。
    // 例: x=1.7e155, y=1.7e116 の塵セル [rog 1e-193] が case/44 run_0415 の確定場に残った)。
    double bx = x0, by = y0;
    bx = fmin(fmax(bx, 0.0), 1.0); by = fmin(fmax(by, bx*bx), sqrt(bx));
    double bd = (bx-x0)*(bx-x0) + (by-y0)*(by-y0);
    if (!(bd == bd) || bd > 1.0e300) bd = 1.0e300;   // 距離が overflow しても後続候補と比較できるようにする
    auto consider = [&](double cx, double cy) {
        cx = fmin(fmax(cx, 0.0), 1.0); cy = fmin(fmax(cy, cx*cx), sqrt(cx));
        const double d = (cx-x0)*(cx-x0) + (cy-y0)*(cy-y0);
        if (d < bd) { bd = d; bx = cx; by = cy; }
    };
    consider(x0, y0);                       // y をクランプ (x は保持)
    consider(fmin(fmax(y0*y0, 0.0), 1.0), y0);   // y=√x の枝で x を合わせる
    consider(sqrt(fmax(y0, 0.0)), y0);           // y=x² の枝で x を合わせる
    consider(1.0, 1.0);
    auto nearest = [&](int branch) {
        const double scale = fmax(fmax(x0, y0*y0), fmax(sqrt(fmax(y0, 0.0)), 1.0e-300));
        double lo = 0.0, hi = fmin(1.0, 8.0*scale + 1.0e-300); const double gr = 0.6180339887498949;
        auto dist2 = [&](double t) { const double yy = (branch == 0) ? sqrt(t) : t*t; return (t-x0)*(t-x0) + (yy-y0)*(yy-y0); };
        double a = hi - gr*(hi-lo), b = lo + gr*(hi-lo), fa = dist2(a), fb = dist2(b);
        for (int it = 0; it < 200 && (hi - lo) > 1.0e-14*hi; ++it) { if (fa < fb) { hi = b; b = a; fb = fa; a = hi - gr*(hi-lo); fa = dist2(a); } else { lo = a; a = b; fa = fb; b = lo + gr*(hi-lo); fb = dist2(b); } }
        const double t = 0.5*(lo+hi);
        consider(t, (branch == 0) ? sqrt(t) : t*t);
    };
    nearest(0); nearest(1);
    x = bx; y = by;
    return 1;
}

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
    flow_float* diagCorrG, flow_float* diagCorrQ,   // 補正量の記録 (nullptr 可): |Δρg|/ρ [質量分率] の累積, Q の最大相対補正 (codex result M2)
    int* realizViol,                                // [0] 最近点射影の作動数, [1] 退化の単分散再初期化数 (nullptr 可; plan species-passive-scalar-unification §4.7 v5)
    double* budget, const geom_int* root, const geom_float* vol,   // 成分別収支 [g,Q0,Q1,Q2]×(符号付き, 絶対)·V (root のみ; nullptr 可)
    int doProject)                                  // 1: 実現可能性射影も行う (定常/RK の各 step; dual-time は sub-iter 内 0 で step 末尾に project_only を呼ぶ)
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
    // 実現可能性の射影 (cond_realizability_project; plan §4.7 v5)。Q0>0, g>0 のみ。(Q0, g) は保存し Q1, Q2 だけ動かす。
    // Q0=0 または g=0 で他が正の「塵」は下の消滅規則 (S≤1) に任せる (核生成域は触らない)。
    {
        const double q0 = (double)roQ0[ic], rg = (double)r;
        if (q0 > 0.0 && rg > 0.0) {
            const CondSpeciesProps cp0 = condProps_make(condModel, opts);
            const double rho_l = cond_rho_cond(cp0, (double)T[ic]);
            const double q3 = rg / ((4.0/3.0)*COND_PI*rho_l);
            const double rr = cbrt(q3/q0);
            double x = 0.0, y = 0.0;
            const bool haveXY = cond_moment_xy(q3, q0, (double)roQ1[ic], (double)roQ2[ic], x, y);   // 指数分離 (codex result-4 M1)
            if (doProject != 0 && !haveXY) {
                // Q3 が表現できない塵 (r=0): 実現可能な唯一の状態は Q1=Q2=0
                if (!(roQ1[ic] == (flow_float)0.0 && roQ2[ic] == (flow_float)0.0)) {
                    roQ1[ic] = (flow_float)0.0; roQ2[ic] = (flow_float)0.0;
                    if (realizViol != nullptr) atomicAdd(realizViol + 1, 1);
                }
            }
            const int kind = (doProject != 0 && haveXY) ? cond_realizability_project(x, y, 1.0e-6) : 0;
            if (kind != 0) {
                double q1n = 0.0, q2n = 0.0; cond_moment_writeback(q3, q0, x, y, q1n, q2n);
                roQ1[ic] = (flow_float)q1n; roQ2[ic] = (flow_float)q2n;
                if (realizViol != nullptr) atomicAdd(realizViol + (kind == 2 ? 1 : 0), 1);
            }
        }
    }
    auto record = [&]() {
        if (diagCorrG == nullptr) return;
        const double rod0 = (double)ro[ic] > 1.0e-20 ? (double)ro[ic] : 1.0e-20;
        diagCorrG[ic] += (flow_float)(fabs((double)rog[ic] - (double)r_in)/rod0);
        double rq = 0.0;
        const flow_float qin[3] = {q0_in, q1_in, q2_in}; const flow_float qout[3] = {roQ0[ic], roQ1[ic], roQ2[ic]};
        for (int k = 0; k < 3; ++k) { const double den = fabs((double)qin[k]) > 1.0e-30 ? fabs((double)qin[k]) : 1.0e-30;
            const double rel = fabs((double)qout[k] - (double)qin[k])/den; if (qout[k] != qin[k] && rel > rq) rq = rel; }
        if ((flow_float)rq > diagCorrQ[ic]) diagCorrQ[ic] = (flow_float)rq;
        if (budget != nullptr && (root == nullptr || root[ic] == ic)) {
            const double V = (vol != nullptr) ? (double)vol[ic] : 1.0;
            const double d[4] = {(double)rog[ic] - (double)r_in, (double)roQ0[ic] - (double)q0_in, (double)roQ1[ic] - (double)q1_in, (double)roQ2[ic] - (double)q2_in};
            for (int k = 0; k < 4; ++k) if (d[k] != 0.0) { atomicAdd(&budget[2*k], d[k]*V); atomicAdd(&budget[2*k+1], fabs(d[k])*V); }
        }
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
    flow_float* diagCorrG, flow_float* diagCorrQ, int* realizViol,
    double* budget, const geom_int* root, const geom_float* vol, int doProject)
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
    {   // 実現可能性の射影 (double 実体と同じ規則; ρ_l は表 [範囲外は double 関数]; 射影は double で評価)
        const float q0 = roQ0[ic];
        if (q0 > 0.0f && r > 0.0f) {
            // 射影の物性は double の式 (表との差 ~1e-4 がゲートの許容 1e-6 を超えないように; 判定は同じ式のゲートと一致させる)
            const double rho_l = cond_rho_cond(cpd, (double)T[ic]);
            const double q3 = (double)r / ((4.0/3.0)*COND_PI*rho_l);
            const double rr = cbrt(q3/(double)q0);
            double x = 0.0, y = 0.0;
            const bool haveXY = cond_moment_xy(q3, (double)q0, (double)roQ1[ic], (double)roQ2[ic], x, y);   // 指数分離 (codex result-4 M1)
            if (doProject != 0 && !haveXY) {
                if (!(roQ1[ic] == 0.0f && roQ2[ic] == 0.0f)) {
                    roQ1[ic] = 0.0f; roQ2[ic] = 0.0f;
                    if (realizViol != nullptr) atomicAdd(realizViol + 1, 1);
                }
            }
            const int kind = (doProject != 0 && haveXY) ? cond_realizability_project(x, y, 1.0e-6) : 0;
            if (kind != 0) {
                double q1n = 0.0, q2n = 0.0; cond_moment_writeback(q3, (double)q0, x, y, q1n, q2n);
                roQ1[ic] = (float)q1n; roQ2[ic] = (float)q2n;
                if (realizViol != nullptr) atomicAdd(realizViol + (kind == 2 ? 1 : 0), 1);
            }
        }
    }
    auto record = [&]() {
        if (diagCorrG == nullptr) return;
        const float rod0 = ro[ic] > 1.0e-20f ? ro[ic] : 1.0e-20f;
        diagCorrG[ic] += fabsf(rog[ic] - r_in)/rod0;
        float rq = 0.0f;
        const flow_float qin[3] = {q0_in, q1_in, q2_in}; const flow_float qout[3] = {roQ0[ic], roQ1[ic], roQ2[ic]};
        for (int k = 0; k < 3; ++k) { const float den = fabsf(qin[k]) > 1.0e-30f ? fabsf(qin[k]) : 1.0e-30f;
            const float rel = fabsf(qout[k] - qin[k])/den; if (qout[k] != qin[k] && rel > rq) rq = rel; }
        if (rq > diagCorrQ[ic]) diagCorrQ[ic] = rq;
        if (budget != nullptr && (root == nullptr || root[ic] == ic)) {
            const double V = (vol != nullptr) ? (double)vol[ic] : 1.0;
            const double d[4] = {(double)rog[ic] - (double)r_in, (double)roQ0[ic] - (double)q0_in, (double)roQ1[ic] - (double)q1_in, (double)roQ2[ic] - (double)q2_in};
            for (int k = 0; k < 4; ++k) if (d[k] != 0.0) { atomicAdd(&budget[2*k], d[k]*V); atomicAdd(&budget[2*k+1], fabs(d[k])*V); }
        }
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


// Q1/Q2 だけの実現可能性射影 (g・Q0 は不変; dual-time FCT の後処理で EOS 更新後の T を使って呼ぶ; plan-7 M4)。budget は [g,Q0,Q1,Q2]×(符号付き,絶対)·V。
__global__ void cond_realizability_project_only_d(
    geom_int nCells, flow_float* ro, flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    int condModel, flow_float* T, CondPropOpts opts, int useTab, CondTablesF tb,
    flow_float* diagCorrQ, int* realizViol, double* budget, const geom_int* root, const geom_float* vol)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double q0 = (double)roQ0[ic], rg = (double)rog[ic];
    if (!(q0 > 0.0 && rg > 0.0)) return;
    const CondSpeciesProps cp0 = condProps_make(condModel, opts);
    const flow_float Tq = T[ic];
    const double rho_l = (useTab != 0 && cond_tab_wet_ok(tb, Tq)) ? (double)cond_tab_rhol_f(tb, Tq) : cond_rho_cond(cp0, (double)Tq);
    const double q3 = rg / ((4.0/3.0)*COND_PI*rho_l);
    const double rr = cbrt(q3/q0);
    double x = 0.0, y = 0.0;
    const flow_float q1_in = roQ1[ic], q2_in = roQ2[ic];
    int kind = 0;
    if (!cond_moment_xy(q3, q0, (double)q1_in, (double)q2_in, x, y)) {
        // Q3 が表現できない (underflow) / x,y が作れない: 液滴半径 0 の塵 → 実現可能な唯一の状態は Q1=Q2=0
        if (q1_in == (flow_float)0.0 && q2_in == (flow_float)0.0) return;
        roQ1[ic] = (flow_float)0.0; roQ2[ic] = (flow_float)0.0; kind = 2;   // 収支・診断は下の共通経路で記録する
    } else {
        kind = cond_realizability_project(x, y, 1.0e-6);
        if (kind == 0) return;
        double q1n = 0.0, q2n = 0.0; cond_moment_writeback(q3, q0, x, y, q1n, q2n);
        roQ1[ic] = (flow_float)q1n; roQ2[ic] = (flow_float)q2n;
    }
    if (realizViol != nullptr) atomicAdd(realizViol + (kind == 2 ? 1 : 0), 1);
    if (diagCorrQ != nullptr) {
        double rq = 0.0;
        const double d1 = fabs((double)roQ1[ic] - (double)q1_in)/fmax(fabs((double)q1_in), 1.0e-30), d2 = fabs((double)roQ2[ic] - (double)q2_in)/fmax(fabs((double)q2_in), 1.0e-30);
        rq = fmax(d1, d2); if ((flow_float)rq > diagCorrQ[ic]) diagCorrQ[ic] = (flow_float)rq;
    }
    if (budget != nullptr && (root == nullptr || root[ic] == ic)) {
        const double V = (vol != nullptr) ? (double)vol[ic] : 1.0;
        const double d[2] = {(double)roQ1[ic] - (double)q1_in, (double)roQ2[ic] - (double)q2_in};
        for (int k = 0; k < 2; ++k) if (d[k] != 0.0) { atomicAdd(&budget[2*(2+k)], d[k]*V); atomicAdd(&budget[2*(2+k)+1], fabs(d[k])*V); }
    }
}
