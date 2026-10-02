// test_twophase_real_source.cpp — 実際の凝縮ソース (condFloat) と採用予定の二相更新の接続検証 (定常擬似時間, 1D)
//   plans/active/condensation-two-phase-transport.md §5.1 #4c。合格条件は notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md §10
//   (実行前に固定)。CFD ではない (ρ・U 固定の 1D ダクト)。
//
// build (host のみ; nvcc で __host__ __device__ 関数を host 側に実体化する。GPU は使わない):
//   cmake -DIN=solver_density_cuda/data/species/forge_species_v1.yaml -DOUT=<GEN>/forge_species_data.hpp -P solver_density_cuda/cmake/embed_species_data.cmake
//   cd solver_density_cuda
//   nice -n 19 g++ -O2 -std=c++17 -c -I . -I <GEN> input/speciesDB.cpp -o <OUT>/speciesDB.o   (speciesTransportDB.cpp も同様)
//   cp tests/unit/test_twophase_real_source.cpp <OUT>/test_twophase_real_source.cu   (nvcc の -x cu は .o にも掛かるので拡張子で .cu 扱いにする)
//   nice -n 19 nvcc -O2 --expt-relaxed-constexpr -std=c++17 -I . -I <GEN> -o <OUT>/test_twophase_real_source \
//       <OUT>/test_twophase_real_source.cu <OUT>/speciesDB.o <OUT>/speciesTransportDB.o -lyaml-cpp
//   <OUT>/test_twophase_real_source [--quick]
//
// 呼ぶ実関数 (host/device; condensationSourceF_d.cuh): cond_vapor_state_f, cond_nucleation_f, cond_growth_f, cond_source_vector_f,
//   cond_evap_source_rate_f, cond_tab_*_f (物性表 condensationTables_d.cuh)。EOS は cond_twophase_resid (condensationEOS_d.cuh)。
// ソースの組み立て (θ_src・src_jac) は float カーネル condensationSourceKernels_d.cuh:296-560 (__global__ で host から呼べない) の本体を
//   同じ順序・同じ関数でここに写す (cond_source_cell_f)。違い: 蒸発の src_jac は double 実体ではなく float の差分 (蒸発セルの Jacobian; 前処理にしか効かない)。
//
// 問題 (§10.1): 1D 一定断面ダクト、ρ・U 固定、N=40・dx=0.5 mm、入口 Dirichlet・出口ゼロ勾配、定常擬似時間 (Δτ = 5 dx/U)。
//   輸送量 ρY_N2, ρY_w (総水分), ρg, ρQ0..2, ρE (内部エネルギー)。1 次風上の移流 + 乱流拡散 Γ = μ_t/Sc_t (全量共通) + 気相の分子拡散
//   (z 基準, 風上の補正; D 一定) + エネルギーは移流 + 種のエンタルピー流束 (h_N2 J_N2 + h_v J_w − L J_l)。熱伝導なし。
//   入口状態: case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/res_48000.h5 の軸近傍 (|y| < 1.5e-4) で g が初めて 1e-6 を超える点
//   (節点 18335, x = 10.65 mm): g 1.190e-6, S 78.78, T 221.350 K, P 23837.6 Pa, ρ 0.360648, Y_w 0.01095, Q0 8.6487e15, Q1 1.8865e7,
//   Q2 6.3338e-2, U 369.99 m/s, μ_t 1.825e-5 (抽出は本ファイル末尾のコメントの python)。気相は内蔵 N2 (run の MIXDRY の代わり) + H2O。
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cfloat>
#include <vector>
#include <string>
#include <algorithm>
#include "cuda_forge/condensationTables_d.cuh"
#include "cuda_forge/condensationSourceF_d.cuh"
#include "cuda_forge/condensationEOS_d.cuh"
#include "tests/unit/cond_latent_test_helper.cuh"

static int g_fail = 0;
static void verdict(bool ok, const char* msg) { printf("[%s] %s\n", ok ? "PASS" : "FAIL", msg); if (!ok) ++g_fail; }

// ------------------------------------------------------------------ 定数・入口状態
static const int    NC = 40;
static const double DX = 0.5e-3;
static double CFL = 5.0;   // 判定は 5 (§10.1)。感度 (判定外) で 50・500 も回す
static const double RHO = 0.360648, U = 369.99, MUT = 1.825e-5, SCT = 0.9, DMOL = 9.0e-5;   // DMOL: N2–H2O, 23.8 kPa, 221 K の目安 (一定)
static const double T_IN = 221.350, YW_IN = 0.01095, G_IN = 1.190e-6, Q0_IN = 8.6487e15, Q1_IN = 1.8865e7, Q2_IN = 6.3338e-2;
static const double DG_MAX = 5.0e-3;   // B_LIMITS (forge 既定 condDgMaxStep)
static double DT_MAX = 1.0;            // forge 既定 condDTmaxStep。#4d の制限作動 A/B だけが 0.01 K に替える
static const double EPS32 = 1.1920928955078125e-7;
enum { K_N2 = 0, K_W, K_L, K_Q0, K_Q1, K_Q2, K_E, NVAR };   // 格納変数 (ρY_N2, ρY_w, ρg, ρQ0, ρQ1, ρQ2, ρE)
static const char* VNAME[] = {"N2", "w", "l", "Q0", "Q1", "Q2", "E"};

struct Ctx {
    SpeciesThermo sp[2];          // [N2, H2O]
    CondSpeciesProps cp;          // H2O (double, 潜熱の気液ペア)
    CondSpeciesPropsF cpf;
    CondTablesHost ht; CondTablesF tb;
    double Rw, RN2;
};

// EOS: (ρY, ρg, ρE) → T (double Newton; cond_twophase_resid の残差 ≤ 1e-9|e| + 0.05 J/kg)
static double eos_T(const Ctx& c, double rN2, double rYw, double rg, double rE, double Tg, bool* ok)
{
    const double Y[2] = {rN2/RHO, rYw/RHO};
    const double g = rg/RHO, e = rE/RHO;
    double T = Tg;
    for (int it = 0; it < 60; ++it) {
        double dG;
        const double G = cond_twophase_resid(c.sp, 2, Y, T, g, c.Rw, 1, c.cp, e, &dG);
        if (fabs(G) <= 1e-9*fabs(e) + 0.05) { *ok = true; return T; }
        T -= G/dG;
    }
    *ok = false; return T;
}
static double e_of_T(const Ctx& c, double rN2, double rYw, double rg, double T)
{
    const double Y[2] = {rN2/RHO, rYw/RHO};
    return cond_twophase_resid(c.sp, 2, Y, T, rg/RHO, c.Rw, 1, c.cp, 0.0) * RHO;   // e(T) − 0
}

// ------------------------------------------------------------------ ソース (float カーネル本体の写し)
struct Src { float Sg, SQ0, SQ1, SQ2, sjg, sjq1, theta; float dSg_drv, dSg_drg; bool evap; bool oot; };   // oot: 物性表の範囲外 (未対応 → FAIL)
static Src cond_source_cell_f(const Ctx& c, float rod, float Td, float rYw, float rog, float q0, float q1, float q2, bool jac2)
{
    Src o{}; o.theta = 1.0f;
    const float Rw = (float)c.Rw;
    float Yw = rYw/rod;
    float g = rog/rod; if (g < 0.0f) g = 0.0f; if (g > Yw) g = (Yw > 0.0f ? Yw : 0.0f);
    if (q0 < 0.0f) q0 = 0.0f; if (q1 < 0.0f) q1 = 0.0f; if (q2 < 0.0f) q2 = 0.0f;
    // 圧力と気相の比熱 (TP: P = ρ T (R_mix − g R_w), cp_cell = 全蒸気混合の c_p)
    const double Yd[2] = {1.0 - (double)Yw, (double)Yw};
    double cpm, hm; thermo_cph_mix(c.sp, 2, Yd, (double)Td, &cpm, &hm);
    const double Rmix = thermo_R_mix(c.sp, 2, Yd);
    const float Pd = (float)((double)rod*(double)Td*(Rmix - (double)g*c.Rw));
    const float cvg = (float)std::max(cpm - Rmix, 1e-3);
    float pv, rho_v; cond_vapor_state_f(1, rod, Pd, Td, g, Yw, Rw, &pv, &rho_v);
    const float gamma_gas = c.cpf.cp/c.cpf.cv;   // condKantrowitzGammaMode 0 (既定)
    const float p_gas = Pd;
    const float lnS = (pv > 0.0f) ? (logf(pv) - cond_tab_lnpsat_f(c.tb, Td)) : -1.0e30f;
    CondNucCarrierF car; car.a_v = 1.0f; car.carrierSum = 0.0f; car.cvv_tilde = c.cpf.cv*c.cpf.M/COND_RU_F;
    const int kw = 1, growth = 0; const float gyarC = 3.18f;
    if (!(Td >= c.tb.Tmin && Td + 0.1f <= c.tb.TwetMax)) { o.oot = true; return o; }   // ゼロソースで続けず、呼び出し側が FAIL にする (#4d)
    if (!(lnS > 0.0f) && g <= 0.0f && q0 <= 1.0e-30f) return o;   // dry
    if (g > 0.0f && !(lnS > 0.0f)) {
        // 蒸発 (limiterMode 1: 率形)
        float r30, drdt;
        cond_evap_source_rate_f(c.cpf, c.tb, Td, pv, rod, g, q0, q1, q2, 1.0e-9f, growth, p_gas, gyarC, 0, &o.SQ0, &o.SQ1, &o.SQ2, &o.Sg, &r30, &drdt);
        o.evap = true;
        // src_jac (蒸発): float 差分 (g と T) — カーネルは double 実体で差分する (前処理にしか効かない)
        const float dgp = 1.0e-3f*fmaxf(g, 1.0e-9f);
        float pvg, rvg; cond_vapor_state_f(1, rod, Pd, Td, g - dgp, Yw, Rw, &pvg, &rvg);
        float a0, a1, a2, agg; cond_evap_source_rate_f(c.cpf, c.tb, Td, pvg, rod, g - dgp, q0, q1, q2, 1.0e-9f, growth, p_gas, gyarC, 0, &a0, &a1, &a2, &agg, &r30, &drdt);
        float pvT, rvT; cond_vapor_state_f(1, rod, Pd, Td + 0.1f, g, Yw, Rw, &pvT, &rvT);
        float agT; cond_evap_source_rate_f(c.cpf, c.tb, Td + 0.1f, pvT, rod, g, q0, q1, q2, 1.0e-9f, growth, p_gas, gyarC, 0, &a0, &a1, &a2, &agT, &r30, &drdt);
        const float L = cond_tab_latent_f(c.tb, Td);
        const float dSgdrog = (o.Sg - agg)/(rod*dgp);
        const float dTdrog = (L - Rw*Td)/(rod*cvg);
        const float sjg = -(dSgdrog + (agT - o.Sg)/0.1f*dTdrog);
        o.sjg = sjg > 0.0f ? sjg : 0.0f;
        o.dSg_drg = dSgdrog + (agT - o.Sg)/0.1f*dTdrog;
        o.dSg_drv = 0.0f;
        if (jac2) {   // 蒸気 (Y_w を ρv ぶん動かす; g 固定、T 固定)
            const float dv = 1.0e-3f*fmaxf(Yw - g, 1.0e-9f);
            float pvv, rvv; cond_vapor_state_f(1, rod, Pd, Td, g, Yw + dv, Rw, &pvv, &rvv);
            float bg; cond_evap_source_rate_f(c.cpf, c.tb, Td, pvv, rod, g, q0, q1, q2, 1.0e-9f, growth, p_gas, gyarC, 0, &a0, &a1, &a2, &bg, &r30, &drdt);
            o.dSg_drv = (bg - o.Sg)/(rod*dv);
        }
        return o;
    }
    // 核生成・成長
    float J, rstar; cond_nucleation_f(c.cpf, c.tb, Td, pv, rho_v, &J, &rstar, kw, gamma_gas, &car);
    if (J < 0.0f) J = 0.0f;
    float r_bar = (q0 > 1.0e-30f) ? (q1/q0) : rstar;
    float drdt = 0.0f;
    if (q0 > 1.0e-30f && rstar > 0.0f && r_bar > rstar) { drdt = cond_growth_f(c.cpf, c.tb, Td, pv, r_bar, rstar, growth, p_gas, gyarC); if (drdt < 0.0f) drdt = 0.0f; }
    const float rho_l = cond_tab_rhol_f(c.tb, Td);
    const float r_nuc = COND_RNUC_FAC_F*rstar;
    float SQ0 = J, SQ1 = J*r_nuc + q0*drdt, SQ2 = J*r_nuc*r_nuc + 2.0f*q1*drdt;
    float Sg = (4.0f/3.0f)*COND_PI_F*rho_l*(J*r_nuc*r_nuc*r_nuc + 3.0f*q2*drdt);
    if (Sg < 0.0f) Sg = 0.0f;
    const float L = cond_tab_latent_f(c.tb, Td);
    float theta = 1.0f;
    const float avail = Yw - g;
    if (Sg > 0.0f && avail <= 0.0f) theta = 0.0f;   // limiterMode 1: 蒸気枯渇のみ
    {
        const float dTp = 0.1f;
        float pvp, rvp; cond_vapor_state_f(1, rod, Pd, Td + dTp, g, Yw, Rw, &pvp, &rvp);
        float a0, a1, a2, ag; cond_source_vector_f(c.cpf, c.tb, Td + dTp, pvp, rvp, q0, q1, q2, &a0, &a1, &a2, &ag, kw, growth, gamma_gas, p_gas, gyarC, &car);
        if (ag < 0.0f) ag = 0.0f;
        const float dSgdT = (ag - Sg)/dTp;
        const float dTdrog = (L - Rw*Td)/(rod*cvg);
        const float sjg = -theta*dSgdT*dTdrog;
        o.sjg = fmaxf(sjg, 0.0f);
        o.dSg_drg = theta*dSgdT*dTdrog;
        if (q0 > 1.0e-30f) {
            const float dq1 = (q1 > 0.0f ? 0.01f*q1 : 1.0e-3f);
            float b0, b1, b2, bg; cond_source_vector_f(c.cpf, c.tb, Td, pv, rho_v, q0, q1 + dq1, q2, &b0, &b1, &b2, &bg, kw, growth, gamma_gas, p_gas, gyarC, &car);
            o.sjq1 = fmaxf(-theta*(b1 - SQ1)/dq1, 0.0f);
        }
        if (jac2) {   // ∂S_g/∂ρv (Y_w を動かす; g・T 固定)
            const float dv = 1.0e-3f*fmaxf(Yw - g, 1.0e-9f);
            float pvv, rvv; cond_vapor_state_f(1, rod, Pd, Td, g, Yw + dv, Rw, &pvv, &rvv);
            float b0, b1, b2, bg; cond_source_vector_f(c.cpf, c.tb, Td, pvv, rvv, q0, q1, q2, &b0, &b1, &b2, &bg, kw, growth, gamma_gas, p_gas, gyarC, &car);
            if (bg < 0.0f) bg = 0.0f;
            o.dSg_drv = theta*(bg - Sg)/(rod*dv);
        }
    }
    o.theta = theta;
    o.SQ0 = SQ0*theta; o.SQ1 = SQ1*theta; o.SQ2 = SQ2*theta; o.Sg = Sg*theta;
    return o;
}

// ------------------------------------------------------------------ 残差の組み立て (精度 F)
template <typename F> struct State { std::vector<F> v[NVAR]; State() { for (auto& a : v) a.assign(NC, F(0)); } };
struct Assembled {
    std::vector<double> r[NVAR + 1];   // + 蒸気 (index NVAR) = r_w − r_l
    std::vector<double> A[NVAR + 1];   // 項の大きさ (床の尺度)
    std::vector<double> diag[NVAR + 1];// 輸送の点対角 (M を除く)
    std::vector<Src> src;
};
static const int K_V = NVAR;

template <typename F>
static void assemble(const Ctx& c, const State<F>& st, const std::vector<double>& T, const State<F>& inlet, Assembled& out, bool jac2)
{
    for (int q = 0; q <= NVAR; ++q) { out.r[q].assign(NC, 0.0); out.A[q].assign(NC, 0.0); out.diag[q].assign(NC, 0.0); }
    out.src.assign(NC, Src{});
    const F rho = F(RHO), mdot = F(RHO*U), Gam = F(MUT/SCT/DX), Dm = F(DMOL/DX);
    // セル値 (ghost = 入口 index -1)
    auto cell = [&](int i, int q) -> F { return (i < 0) ? inlet.v[q][0] : st.v[q][i]; };
    auto Tc = [&](int i) -> double { return (i < 0) ? T_IN : T[i]; };
    std::vector<F> R[NVAR + 1]; std::vector<double> Aa[NVAR + 1];
    for (int q = 0; q <= NVAR; ++q) { R[q].assign(NC, F(0)); Aa[q].assign(NC, 0.0); }
    // 面 f = i−1/2 (i = 0..NC-1; i=0 は入口面) と出口面
    for (int i = 0; i <= NC; ++i) {
        const int L = i - 1, Rr = i;           // 左右セル (Rr = NC は出口: 流出のみ)
        F flux[NVAR + 1] = {};
        // 移流 (ṁ > 0, 風上 = 左)
        const F yN2 = cell(L, K_N2)/rho, yw = cell(L, K_W)/rho, gl = cell(L, K_L)/rho, e = cell(L, K_E)/rho;
        flux[K_N2] = mdot*yN2; flux[K_L] = mdot*gl;
        const F yv = (cell(L, K_W) - cell(L, K_L))/rho;
        flux[K_V] = mdot*yv; flux[K_W] = flux[K_V] + flux[K_L];
        flux[K_Q0] = mdot*(cell(L, K_Q0)/rho); flux[K_Q1] = mdot*(cell(L, K_Q1)/rho); flux[K_Q2] = mdot*(cell(L, K_Q2)/rho);
        flux[K_E] = mdot*e;
        double fabsAdv[NVAR + 1]; for (int q = 0; q <= NVAR; ++q) fabsAdv[q] = fabs((double)flux[q]);
        F dif[NVAR + 1] = {};
        if (Rr < NC) {
            // 拡散 (乱流: 全量共通; 分子: 気相 z、補正は風上)
            auto d = [&](int q) { return -Gam*((cell(Rr, q) - cell(L, q))/rho); };
            const F gL = cell(L, K_L), gR = cell(Rr, K_L);
            const F rgL = rho - gL, rgR = rho - gR;
            const F zvL = (cell(L, K_W) - gL)/rgL, zvR = (cell(Rr, K_W) - gR)/rgR;
            const F zNL = cell(L, K_N2)/rgL, zNR = cell(Rr, K_N2)/rgR;
            const F rgf = F(0.5)*(rgL + rgR);
            const F j0v = -rgf*Dm*(zvR - zvL), j0N = -rgf*Dm*(zNR - zNL);
            const F S = j0v + j0N;
            const bool upL = (S <= F(0));   // 補正流束 −S の流出側
            const F jv = j0v - (upL ? zvL : zvR)*S, jN = j0N - (upL ? zNL : zNR)*S;
            const F JtN = d(K_N2), Jtl = d(K_L);
            const F Jtv = -Gam*(((cell(Rr, K_W) - gR) - (cell(L, K_W) - gL))/rho);
            dif[K_N2] = jN + JtN; dif[K_L] = Jtl; dif[K_V] = jv + Jtv; dif[K_W] = dif[K_V] + dif[K_L];
            dif[K_Q0] = d(K_Q0); dif[K_Q1] = d(K_Q1); dif[K_Q2] = d(K_Q2);
            const double Tf = 0.5*(Tc(L) + Tc(Rr));
            const double hN = thermo_h_mass(c.sp[0], Tf), hv = thermo_h_mass(c.sp[1], Tf), Lf = cond_latent(c.cp, Tf);
            dif[K_E] = F(hN*(double)dif[K_N2] + hv*(double)dif[K_W] - Lf*(double)dif[K_L]);
            const double dE_abs = fabs(hN*(double)dif[K_N2]) + fabs(hv*(double)dif[K_W]) + fabs(Lf*(double)dif[K_L]);
            for (int q = 0; q <= NVAR; ++q) { flux[q] += dif[q]; fabsAdv[q] += (q == K_E) ? dE_abs : fabs((double)dif[q]); }
            // 点対角 (輸送): 左セルは流出 ṁ/ρ と拡散、右セルは拡散
            const double gd = (double)Gam/RHO, md = (double)Dm*(double)rgf;
            for (int q = 0; q <= NVAR; ++q) {
                double dl = gd, dr = gd;
                if (q == K_V || q == K_N2) { dl += md/(double)(q == K_V ? rgL : rgL); dr += md/(double)rgR; }
                if (q == K_E) { dl = 0.0; dr = 0.0; }
                if (L >= 0) out.diag[q][L] += dl;
                out.diag[q][Rr] += dr;
            }
            const double sabs = fabs((double)S);
            if (L >= 0 && upL) { out.diag[K_V][L] += sabs/(double)rgL; out.diag[K_N2][L] += sabs/(double)rgL; }
            if (!upL) { out.diag[K_V][Rr] += sabs/(double)rgR; out.diag[K_N2][Rr] += sabs/(double)rgR; }
        }
        if (L >= 0) for (int q = 0; q <= NVAR; ++q) { out.diag[q][L] += RHO*U/RHO; }
        for (int q = 0; q <= NVAR; ++q) {
            if (L >= 0) { R[q][L] -= flux[q]; Aa[q][L] += fabsAdv[q]; }
            if (Rr < NC) { R[q][Rr] += flux[q]; Aa[q][Rr] += fabsAdv[q]; }
        }
    }
    // ソース (体積 V = DX; 液・Q の残差にだけ入る。蒸気は r_v = r_w − r_l で −S を受ける)
    for (int i = 0; i < NC; ++i) {
        const Src s = cond_source_cell_f(c, (float)RHO, (float)T[i], (float)st.v[K_W][i], (float)st.v[K_L][i],
                                         (float)st.v[K_Q0][i], (float)st.v[K_Q1][i], (float)st.v[K_Q2][i], jac2);
        out.src[i] = s;
        const F V = F(DX);
        R[K_L][i] += F(s.Sg)*V; R[K_Q0][i] += F(s.SQ0)*V; R[K_Q1][i] += F(s.SQ1)*V; R[K_Q2][i] += F(s.SQ2)*V;
        Aa[K_L][i] += fabs((double)s.Sg)*DX; Aa[K_Q0][i] += fabs((double)s.SQ0)*DX; Aa[K_Q1][i] += fabs((double)s.SQ1)*DX; Aa[K_Q2][i] += fabs((double)s.SQ2)*DX;
        Aa[K_V][i] += fabs((double)s.Sg)*DX;
    }
    for (int q = 0; q < NVAR; ++q) for (int i = 0; i < NC; ++i) { out.r[q][i] = (double)R[q][i]; out.A[q][i] = Aa[q][i]; }
    for (int i = 0; i < NC; ++i) { out.r[K_V][i] = (double)(R[K_W][i] - R[K_L][i]); out.A[K_V][i] = Aa[K_V][i]; }   // 全残差変換 (精度 F の減算)
}

// ------------------------------------------------------------------ 反復
template <typename F> static F ulp_of(F x) { return std::nextafter(x, (F)INFINITY) - x; }

struct RunOpt { double omega; bool jac2; int cap; const char* tag; };
struct RunOut {
    bool converged = false; int iters = 0; double ratio_max = 0.0; double comp_ratio[NVAR + 1] = {};
    double theta_min = 1.0; long theta_lt1 = 0; double theta_src_min = 1.0; int theta_final_lt1 = 0, theta_src_final_lt1 = 0;
    double withheld_v = 0, withheld_g = 0, qcut = 0, vround = 0, qcut_last10 = 0, vround_last10 = 0;
    double min_rv = INFINITY, min_rg = INFINITY, min_Q = INFINITY;
    double ulp_max[NVAR + 1] = {}; double frac_sub_half = 0.0; bool stagnate = false;
    int realiz_viol = 0; double g_out = 0, S_out = 0, T_out = 0, x_onset = -1; double sumY_err = 0;
    std::vector<std::string> hist;
    // #4d: 補正は更新ごとに記録し、実際に行った更新数の末尾 10 % (切り上げ、最低 1) で集計する
    std::vector<double> qcut_it, vround_it; int tail_n = 0;
    int theta_lt1_updates = 0, theta_last_it = -1;   // θ<1 のセルがあった更新の数・最後にあった更新
    // #4d: 非有限値・表範囲外は集計で落とさず、検出した時点で打ち切って FAIL
    long oot = 0; std::string bad;
};
// 合否 (§10.2 + #4d): 収束 (独立残差比 ≤1、全成分有限)・非負・末尾 10 % の補正 0・最後の θ と θ_src が 1・非有限/表範囲外なし
static bool accept(const RunOut& x)
{
    if (!x.converged || !x.bad.empty() || x.oot != 0) return false;
    for (int q = 0; q <= NVAR; ++q) if (!std::isfinite(x.comp_ratio[q]) || x.comp_ratio[q] > 1.0) return false;
    return x.min_rv >= 0 && x.min_rg >= 0 && x.min_Q >= 0 && x.qcut_last10 == 0 && x.vround_last10 == 0
        && x.theta_src_final_lt1 == 0 && x.theta_final_lt1 == 0;
}

// 組み立て結果の検査 (#4d): 表範囲外のセル・非有限の残差/尺度/ソースを見つけたら ro.bad に記録して false
static bool scan_assembled(const Assembled& A_, const char* where, int it, RunOut& ro)
{
    char b[192];
    for (int i = 0; i < NC; ++i) {
        const Src& s = A_.src[i];
        if (s.oot) { ++ro.oot; snprintf(b, sizeof b, "%s: 物性表の範囲外 (セル %d, 反復 %d)", where, i, it); ro.bad = b; return false; }
        const float sv[] = {s.Sg, s.SQ0, s.SQ1, s.SQ2, s.sjg, s.sjq1, s.theta, s.dSg_drv, s.dSg_drg};
        for (float v_ : sv) if (!std::isfinite(v_)) { snprintf(b, sizeof b, "%s: 非有限のソース (セル %d, 反復 %d)", where, i, it); ro.bad = b; return false; }
        for (int q = 0; q <= NVAR; ++q)
            if (!std::isfinite(A_.r[q][i]) || !std::isfinite(A_.A[q][i]) || !std::isfinite(A_.diag[q][i])) {
                snprintf(b, sizeof b, "%s: 非有限の残差/尺度 %s (セル %d, 反復 %d)", where, q == K_V ? "v" : VNAME[q], i, it); ro.bad = b; return false; }
    }
    return true;
}

template <typename F>
static RunOut run(const Ctx& c, const RunOpt& o)
{
    RunOut ro;
    State<F> st, inlet;
    inlet = State<F>();
    for (auto& a : inlet.v) a.assign(1, F(0));
    const double rYw = RHO*YW_IN, rg = RHO*G_IN, rN2 = RHO - rYw;
    const double rE = e_of_T(c, rN2, rYw, rg, T_IN);
    const double init[NVAR] = {rN2, rYw, rg, RHO*Q0_IN, RHO*Q1_IN, RHO*Q2_IN, rE};
    for (int q = 0; q < NVAR; ++q) { inlet.v[q][0] = F(init[q]); for (int i = 0; i < NC; ++i) st.v[q][i] = F(init[q]); }
    std::vector<double> T(NC, T_IN);
    const double M = DX/(CFL*DX/U);   // V/Δτ
    Assembled a;
    // 初期残差 r0 と停止の床 (床は各反復の項の大きさで評価)
    double r0[NVAR + 1] = {};
    const int qlist[] = {K_N2, K_W, K_V, K_L, K_Q0, K_Q1, K_Q2, K_E};
    std::vector<int> lastSub(0);
    int subHalfStreak = 0;
    for (int it = 0; it <= o.cap; ++it) {
        // EOS (格納値から double Newton) — 反復の T と独立残差の T は同じ定義
        for (int i = 0; i < NC; ++i) {
            bool ok; T[i] = eos_T(c, (double)st.v[K_N2][i], (double)st.v[K_W][i], (double)st.v[K_L][i], (double)st.v[K_E][i], T[i], &ok);
            if (!ok || !std::isfinite(T[i])) { char b[128]; snprintf(b, sizeof b, "EOS 失敗/非有限 T (セル %d, 反復 %d)", i, it); ro.bad = b; ro.iters = it; return ro; }
        }
        // 独立残差 (float64 で組み直し; 停止判定)
        State<double> s64; for (int q = 0; q < NVAR; ++q) for (int i = 0; i < NC; ++i) s64.v[q][i] = (double)st.v[q][i];
        State<double> in64; for (auto& v_ : in64.v) v_.assign(1, 0.0); for (int q = 0; q < NVAR; ++q) in64.v[q][0] = (double)inlet.v[q][0];
        Assembled a64; assemble<double>(c, s64, T, in64, a64, false);
        if (!scan_assembled(a64, "独立残差", it, ro)) { ro.iters = it; return ro; }
        bool pass = true; double rmax = 0.0;
        for (int q : qlist) {
            double m = 0.0, Am = 0.0; for (int i = 0; i < NC; ++i) { m = std::max(m, fabs(a64.r[q][i])); Am = std::max(Am, a64.A[q][i]); }
            if (it == 0) r0[q] = m;
            const double tol = std::max(1e-7*r0[q], 6.0*EPS32*Am);
            const double ratio = (tol > 0) ? m/tol : (m == 0 ? 0.0 : INFINITY);
            ro.comp_ratio[q] = ratio; rmax = std::max(rmax, ratio);
            if (!(ratio <= 1.0)) pass = false;   // NaN も不合格側へ
        }
        ro.ratio_max = rmax;
        if (it % 100 == 0 || pass || it == o.cap) {
            char b[512]; int n = snprintf(b, sizeof b, "it %5d:", it);
            for (int q : qlist) n += snprintf(b + n, sizeof b - n, " %s %.2e", (q == K_V ? "v" : VNAME[q]), ro.comp_ratio[q]);
            ro.hist.push_back(b);
        }
        if (pass && it > 0) { ro.converged = true; ro.iters = it; break; }
        if (it == o.cap) { ro.iters = it; break; }
        // 反復の残差 (精度 F) と前処理
        assemble<F>(c, st, T, inlet, a, o.jac2);
        if (!scan_assembled(a, "反復残差", it, ro)) { ro.iters = it; return ro; }
        bool allSub = true; int nsub = 0, ntot = 0; int thetaLastLt1 = 0;
        double ul[NVAR + 1] = {};
        double qcut_this = 0.0, vround_this = 0.0;
        for (int i = 0; i < NC; ++i) {
            const Src& s = a.src[i];
            const F Mf = F(M);
            // 増分 (前処理の後、精度 F)
            F d[NVAR + 1];
            const F rv = F(a.r[K_V][i]), rl = F(a.r[K_L][i]);
            const F Dv = Mf + F(a.diag[K_V][i]), Dg = Mf + F(a.diag[K_L][i]);
            F dv, dg;
            if (o.jac2) {
                const F avv = F(DX)*F(s.dSg_drv), avg = F(DX)*F(s.dSg_drg);   // V ∂S/∂ρv, V ∂S/∂ρg
                const F P11 = Dv + avv, P12 = avg, P21 = -avv, P22 = Dg - avg;
                const F det = P11*P22 - P12*P21;
                dv = (P22*rv - P12*rl)/det; dg = (P11*rl - P21*rv)/det;
            } else {
                dv = rv/Dv; dg = rl/(Dg + F(DX)*F(s.sjg));
            }
            d[K_V] = dv; d[K_L] = dg;
            d[K_N2] = F(a.r[K_N2][i])/(Mf + F(a.diag[K_N2][i]));
            d[K_Q0] = F(a.r[K_Q0][i])/(Mf + F(a.diag[K_Q0][i]));
            d[K_Q1] = F(a.r[K_Q1][i])/(Mf + F(a.diag[K_Q1][i]) + F(DX)*F(s.sjq1));
            d[K_Q2] = F(a.r[K_Q2][i])/(Mf + F(a.diag[K_Q2][i]));
            d[K_E] = F(a.r[K_E][i])/(Mf + F(U));
            // 緩和 (前処理の後・limiter の前; 全増分)
            if (o.omega != 1.0) { const F w = F(o.omega); for (int q = 0; q <= NVAR; ++q) if (q != K_W) d[q] *= w; }
            // vl_limit_commit
            const double rvs = (double)st.v[K_W][i] - (double)st.v[K_L][i];
            double th = 1.0;
            const double adg = fabs((double)d[K_L])/RHO;
            if (adg > 0) {
                th = std::min(th, DG_MAX/adg);
                double cpm, hm; const double Yd[2] = {(double)st.v[K_N2][i]/RHO, (double)st.v[K_W][i]/RHO};
                thermo_cph_mix(c.sp, 2, Yd, T[i], &cpm, &hm);
                const double gg = (double)st.v[K_L][i]/RHO;
                const double cveff = cpm - thermo_R_mix(c.sp, 2, Yd) + gg*(c.Rw - (cond_latent(c.cp, T[i] + 0.1) - cond_latent(c.cp, T[i] - 0.1))/0.2);
                th = std::min(th, DT_MAX/(adg*cond_latent(c.cp, T[i])/cveff));
            }
            if ((double)d[K_V] < 0) th = std::min(th, rvs/(-(double)d[K_V]));
            if ((double)d[K_L] < 0) th = std::min(th, (double)st.v[K_L][i]/(-(double)d[K_L]));
            th = std::max(th, 0.0);
            const F Th = F(th);
            {   // 増分・θ の非有限 (#4d)
                bool fin = std::isfinite(th);
                for (int q = 0; q <= NVAR; ++q) if (q != K_W && !std::isfinite((double)d[q])) fin = false;
                if (!fin) { char b[128]; snprintf(b, sizeof b, "非有限の増分/θ (セル %d, 反復 %d)", i, it); ro.bad = b; ro.iters = it; return ro; }
            }
            ro.theta_min = std::min(ro.theta_min, th); if (th < 1.0) ++ro.theta_lt1;
            ro.theta_src_min = std::min(ro.theta_src_min, (double)s.theta);
            ro.withheld_v += (1.0 - th)*fabs((double)d[K_V]); ro.withheld_g += (1.0 - th)*fabs((double)d[K_L]);
            // Q: 共通 θ と成分ごとの非負化
            for (int q : {K_Q0, K_Q1, K_Q2}) {
                F dq = Th*d[q]; F nq = st.v[q][i] + dq;
                if (nq < F(0)) { ro.qcut += -(double)nq; qcut_this += -(double)nq; dq = -st.v[q][i]; }
                d[q] = dq;
            }
            const F gnew = st.v[K_L][i] + Th*d[K_L];
            F wnew = st.v[K_W][i] + (Th*d[K_V] + Th*d[K_L]);
            if (wnew - gnew < F(0)) { ro.vround += (double)(gnew - wnew); vround_this += (double)(gnew - wnew); wnew = gnew; }
            // ULP の記録 (commit 前の値基準)
            const F dw = wnew - st.v[K_W][i], dgc = gnew - st.v[K_L][i];
            auto rec = [&](int q, F incr, F x) { const double u = (double)ulp_of<F>(x); const double r_ = (u > 0) ? fabs((double)incr)/u : 0.0;
                                                 ul[q] = std::max(ul[q], r_); ++ntot; if (r_ < 0.5) ++nsub; if (r_ >= 0.5) allSub = false; };
            rec(K_W, Th*d[K_V] + Th*d[K_L], st.v[K_W][i]); rec(K_L, Th*d[K_L], st.v[K_L][i]);
            rec(K_N2, d[K_N2], st.v[K_N2][i]); rec(K_Q0, d[K_Q0], st.v[K_Q0][i]); rec(K_Q1, d[K_Q1], st.v[K_Q1][i]); rec(K_Q2, d[K_Q2], st.v[K_Q2][i]);
            rec(K_E, d[K_E], st.v[K_E][i]);
            (void)dw; (void)dgc;
            st.v[K_L][i] = gnew; st.v[K_W][i] = wnew;
            st.v[K_N2][i] += d[K_N2]; st.v[K_E][i] += d[K_E];
            for (int q : {K_Q0, K_Q1, K_Q2}) st.v[q][i] += d[q];
            for (int q = 0; q < NVAR; ++q) if (!std::isfinite((double)st.v[q][i])) {   // 状態の非有限 (#4d)
                char b[128]; snprintf(b, sizeof b, "非有限の状態 %s (セル %d, 反復 %d)", VNAME[q], i, it); ro.bad = b; ro.iters = it; return ro; }
            ro.min_rv = std::min(ro.min_rv, (double)(st.v[K_W][i] - st.v[K_L][i]));
            ro.min_rg = std::min(ro.min_rg, (double)st.v[K_L][i]);
            for (int q : {K_Q0, K_Q1, K_Q2}) ro.min_Q = std::min(ro.min_Q, (double)st.v[q][i]);
            if (th < 1.0) ++thetaLastLt1;
        }
        for (int q = 0; q <= NVAR; ++q) ro.ulp_max[q] = ul[q];
        ro.frac_sub_half = ntot ? (double)nsub/ntot : 0.0;
        subHalfStreak = allSub ? subHalfStreak + 1 : 0;
        ro.stagnate = (subHalfStreak >= 100);
        ro.theta_final_lt1 = thetaLastLt1;   // 最後に行った更新の θ<1 セル数 (次の反復で収束したら、これが最終の更新)
        if (thetaLastLt1 > 0) { ++ro.theta_lt1_updates; ro.theta_last_it = it; }
        ro.qcut_it.push_back(qcut_this); ro.vround_it.push_back(vround_this);
    }
    // 補正の末尾 10 % 集計 (実際に行った更新数が基準; 上限 o.cap ではない)
    {
        const int n = (int)ro.qcut_it.size();
        ro.tail_n = std::max(1, (int)std::ceil(0.1*n));
        ro.qcut_last10 = ro.vround_last10 = 0.0;
        for (int k = std::max(0, n - ro.tail_n); k < n; ++k) { ro.qcut_last10 += ro.qcut_it[k]; ro.vround_last10 += ro.vround_it[k]; }
    }
    // 最終状態の監視
    {
        Assembled af; State<double> s64; for (int q = 0; q < NVAR; ++q) for (int i = 0; i < NC; ++i) s64.v[q][i] = (double)st.v[q][i];
        State<double> in64; for (auto& v_ : in64.v) v_.assign(1, 0.0); for (int q = 0; q < NVAR; ++q) in64.v[q][0] = (double)inlet.v[q][0];
        assemble<double>(c, s64, T, in64, af, false);
        if (!scan_assembled(af, "最終状態", ro.iters, ro)) return ro;
        for (int i = 0; i < NC; ++i) {
            if (af.src[i].theta < 1.0f) ++ro.theta_src_final_lt1;
            const double rhol = cond_rho_cond(c.cp, T[i]);
            const double Q3 = s64.v[K_L][i]/((4.0/3.0)*M_PI*rhol);
            if (s64.v[K_Q2][i]*s64.v[K_Q2][i] > s64.v[K_Q1][i]*Q3*(1.0 + 1e-6)) ++ro.realiz_viol;
            if (ro.x_onset < 0 && s64.v[K_L][i]/RHO > 1e-4) ro.x_onset = (i + 0.5)*DX;
            ro.sumY_err = std::max(ro.sumY_err, fabs(s64.v[K_N2][i] + s64.v[K_W][i] - RHO)/RHO);
        }
        const int o_ = NC - 1;
        ro.g_out = s64.v[K_L][o_]/RHO; ro.T_out = T[o_];
        const double Yw = s64.v[K_W][o_]/RHO, g = ro.g_out;
        const double Yd[2] = {1.0 - Yw, Yw};
        const double Pd = RHO*T[o_]*(thermo_R_mix(c.sp, 2, Yd) - g*c.Rw);
        const double pv = RHO*(Yw - g)*c.Rw*T[o_]; (void)Pd;
        ro.S_out = pv/cond_psat(c.cp, T[o_]);
    }
    return ro;
}

template <typename F>
static RunOut report(const Ctx& c, const RunOpt& o, bool judged)
{
    const RunOut r = run<F>(c, o);
    const bool ok = accept(r);
    char b[2400];
    snprintf(b, sizeof b,
        "%s: %s%s%s (反復 %d / 上限 %d)、独立残差の最大比 %.3f [N2 %.2f w %.2f v %.2f l %.2f Q0 %.2f Q1 %.2f Q2 %.2f E %.2f]、"
        "min ρv %.3e・ρg %.3e・ρQ %.3e、θ 最小 %.3f (θ<1 のセル·反復 %ld、最後の更新で θ<1 のセル %d)、θ_src 最小 %.3f (最終で <1 のセル %d)、"
        "θ<1 のあった更新 %d (最後 反復 %d)、保留 Σ(1−θ)|δ| 蒸気 %.2e / 液 %.2e、状態補正 Q_cut %.2e (末尾 %d 更新 %.2e)・v_round %.2e (同 %.2e)、"
        "最終反復の増分/ULP 最大 [w %.2g l %.2g N2 %.2g Q0 %.2g Q1 %.2g Q2 %.2g E %.2g]・0.5 ULP 未満の割合 %.2f%s、"
        "実現可能性違反 %d、出口 g %.4e・S %.3g・T %.2f K、g > 1e-4 の位置 %.2f mm、|ΣρY−ρ|/ρ %.1e",
        o.tag, r.converged ? "収束" : (r.bad.empty() ? "上限到達" : "打ち切り"), r.bad.empty() ? "" : " — ", r.bad.c_str(), r.iters, o.cap, r.ratio_max,
        r.comp_ratio[K_N2], r.comp_ratio[K_W], r.comp_ratio[K_V], r.comp_ratio[K_L], r.comp_ratio[K_Q0], r.comp_ratio[K_Q1], r.comp_ratio[K_Q2], r.comp_ratio[K_E],
        r.min_rv, r.min_rg, r.min_Q, r.theta_min, r.theta_lt1, r.theta_final_lt1, r.theta_src_min, r.theta_src_final_lt1,
        r.theta_lt1_updates, r.theta_last_it, r.withheld_v, r.withheld_g, r.qcut, r.tail_n, r.qcut_last10, r.vround, r.vround_last10,
        r.ulp_max[K_W], r.ulp_max[K_L], r.ulp_max[K_N2], r.ulp_max[K_Q0], r.ulp_max[K_Q1], r.ulp_max[K_Q2], r.ulp_max[K_E],
        r.frac_sub_half, r.stagnate ? " (最後 100 反復すべて 0.5 ULP 未満 = 格納丸めによる停滞)" : "",
        r.realiz_viol, r.g_out, r.S_out, r.T_out, r.x_onset*1e3, r.sumY_err);
    if (judged) verdict(ok, b); else printf("[INFO] %s\n", b);
    for (size_t k = 0; k < r.hist.size(); k += std::max<size_t>(1, r.hist.size()/8)) printf("    %s\n", r.hist[k].c_str());
    if (!r.hist.empty()) printf("    %s\n", r.hist.back().c_str());
    return r;
}

int main(int argc, char** argv)
{
    const bool quick = (argc > 1 && !strcmp(argv[1], "--quick"));
    Ctx c;
    {
        ResolvedSpeciesDB db = speciesDB_resolve(std::vector<std::string>{"N2", "H2O"}, "");
        c.sp[0] = db.species[0]; c.sp[1] = db.species[1];
        for (int k = 0; k < 2; ++k) if (!(c.sp[k].invMW > 0.0)) c.sp[k].invMW = 1.0/c.sp[k].MW;
        c.cp = cond_test_props_H2O(false);
        cond_tables_build_host(c.cp, c.ht); c.tb = cond_tables_view_host(c.ht);
        c.cpf = condProps_to_f(c.cp);
        c.Rw = c.cp.R; c.RN2 = THERMO_RU/c.sp[0].MW;
    }
    const int cap = quick ? 3000 : 20000;
    printf("=== #4c 実ソース (condFloat) × 採用更新: 1D ダクト定常擬似時間 (N %d, dx %.1f mm, CFL %.0f, 上限 %d) ===\n", NC, DX*1e3, CFL, cap);
    printf("[INFO] 入口: T %.3f K, ρ %.6f, Y_w %.5f, g %.3e, Q0 %.4e, Q1 %.4e, Q2 %.4e, U %.2f, μ_t %.3e (run_0482 res_48000 節点 18335)\n",
           T_IN, RHO, YW_IN, G_IN, Q0_IN, Q1_IN, Q2_IN, U, MUT);
    // 判定 (float32): 緩和 × 結合
    RunOut r[2][2];
    const double om[2] = {1.0, 0.5};
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) {
        char tag[128]; snprintf(tag, sizeof tag, "float32 ω=%.1f %s", om[a], b ? "2×2 連成" : "対角 (現行 src_jac 形)");
        RunOpt o{om[a], b == 1, cap, tag};
        r[a][b] = report<float>(c, o, true);
    }
    // 対照 (float64 格納; 判定外)
    for (int b = 0; b < 2; ++b) { char tag[128]; snprintf(tag, sizeof tag, "float64 対照 ω=1.0 %s", b ? "2×2 連成" : "対角"); RunOpt o{1.0, b == 1, cap, tag}; report<double>(c, o, false); }
    // 感度 (判定外; §10 の固定条件の外): 擬似時間刻みを大きくしてソースを硬くする
    for (double cfl : {50.0, 500.0}) {
        CFL = cfl;
        for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) {
            char tag[160]; snprintf(tag, sizeof tag, "感度 CFL %.0f float32 ω=%.1f %s", cfl, om[a], b ? "2×2 連成" : "対角");
            RunOpt o{om[a], b == 1, cap, tag};
            report<float>(c, o, false);
        }
    }
    CFL = 5.0;
    // 比較条件 (§10.4)
    auto pass = [](const RunOut& x) { return accept(x); };
    int ia;
    const bool p10 = pass(r[0][0]), p05 = pass(r[1][0]);
    if (p10 && !(r[0][0].qcut + r[0][0].vround > 0 && p05)) ia = 0; else if (p05) ia = 1; else ia = -1;
    if (ia < 0) printf("[判定] 緩和: ω=1・ω=0.5 とも不合格 (対角) → 採用せず上位へ\n");
    else printf("[判定] 緩和: ω=%.1f を採る (対角で ω=1 %s・ω=0.5 %s)\n", om[ia], p10 ? "合格" : "不合格", p05 ? "合格" : "不合格");
    const int ib = (ia < 0) ? 0 : ia;
    const bool pd = pass(r[ib][0]), p2 = pass(r[ib][1]);
    if (!pd && p2) printf("[判定] 2×2 連成 (ω=%.1f): 対角 不合格・2×2 合格 → 2×2 を必須\n", om[ib]);
    else if (pd && p2) printf("[判定] 2×2 連成 (ω=%.1f): 両方合格 → 必須にしない (反復 対角 %d / 2×2 %d)\n", om[ib], r[ib][0].iters, r[ib][1].iters);
    else if (pd && !p2) printf("[判定] 2×2 連成 (ω=%.1f): 対角 合格・2×2 不合格 → 2×2 は採らない\n", om[ib]);
    else printf("[判定] 2×2 連成 (ω=%.1f): 両方不合格 → 上位へ\n", om[ib]);
    // #4d 制限作動 A/B (変更は DT_MAX だけ; 同じ初期値・CFL 5・ω 1・対角・上限 cap)
    {
        printf("\n=== #4d 制限作動 A/B (DT_MAX のみ変更; float32, CFL %.0f, ω 1, 対角, 上限 %d) ===\n", CFL, cap);
        RunOut ab[2]; const double dtm[2] = {1.0, 0.01};
        for (int k = 0; k < 2; ++k) {
            DT_MAX = dtm[k];
            char tag[128]; snprintf(tag, sizeof tag, "A/B %s: DT_MAX %.2g K", k ? "B" : "A", dtm[k]);
            RunOpt o{1.0, false, cap, tag};
            ab[k] = report<float>(c, o, true);
        }
        DT_MAX = 1.0;
        const bool pA = accept(ab[0]), pB = accept(ab[1]), actB = ab[1].theta_lt1 > 0;
        if (!ab[1].bad.empty() || !ab[0].bad.empty()) printf("[判定] A/B: 打ち切り (A %s / B %s) → 判定不能\n", ab[0].bad.empty() ? "-" : ab[0].bad.c_str(), ab[1].bad.empty() ? "-" : ab[1].bad.c_str());
        else if (!actB) printf("[判定] A/B: B で θ が常に 1 → 判別不成立 (θ_src の作動検証にもならない)\n");
        else if (pA && pB) printf("[判定] A/B: A・B とも達成、B で θ<1 %ld セル·反復 → 当該条件で更新制限が反復を妨げる懸念を除外\n", ab[1].theta_lt1);
        else if (pA && !pB) printf("[判定] A/B: A だけ達成 → 制限・commit を調べる\n");
        else printf("[判定] A/B: A 不達成 (A %s / B %s) → 規則の想定外、上位へ\n", pA ? "合格" : "不合格", pB ? "合格" : "不合格");
        if (ab[1].bad.empty() && !actB) { printf("[FAIL] B で θ<1 が発生しない (判別不成立)\n"); ++g_fail; }
        printf("[INFO] B の θ_src 最小 %.3f (θ_src<1 は蒸気枯渇でしか起きない; DT_MAX は θ_src に作用しない)\n", ab[1].theta_src_min);
    }
    printf("\n%s\n", g_fail ? "FAIL あり" : "ALL PASS");
    return g_fail ? 1 : 0;
}

// 入口状態の抽出 (本体ワークツリー, 読み取りのみ):
//   import h5py, numpy as np
//   f = h5py.File('case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/res_48000.h5'); V = f['VALUE']; X = f['MESH/COORD'][:].reshape(-1, 3)
//   n = V['T'].shape[0]; x, y = X[:n, 0], X[:n, 1]; idx = np.where(np.abs(y) < 1.5e-4)[0]; idx = idx[np.argsort(x[idx])]
//   i = [k for k in idx if V['g_0'][k] > 1e-6][0]   # → 18335 (x = 10.65 mm)
