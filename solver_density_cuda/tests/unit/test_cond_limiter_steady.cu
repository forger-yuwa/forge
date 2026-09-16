// test_cond_limiter_steady.cu — 凝縮ソースの Δτ 不変性と更新クランプの単体試験 (plans/active/condensation-source-limiter-steady.md §5-6)。
//   build: nvcc --expt-relaxed-constexpr -I. -o test_cond_limiter_steady tests/unit/test_cond_limiter_steady.cu
//   (a) condLimiterMode 1: 状態を固定して dt_local を 1e-7 / 1e-3 に振っても、double / float 実体とも res_* と sj_* が
//       ビット一致する (凝縮・蒸発・枯渇を含む H2O carrier と N2 pure の格子)。mode 0 では少なくとも 1 セルで差が出る (旧挙動の確認)。
//   (b) 更新クランプ kernel: 上限内 θ_u=1 で増分不変、(c) 潜熱 ΔT 超過で 4 本同率縮小、(d) 残差 0 で無作用、
//       (f) 蒸気枯渇 θ_u=avail/Δg、負増分の半径半減上限、g=0 への負増分は floor で 0 (補正量 >0)、float 極小増分で有限。
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/thermo_d.cuh"
namespace {
#include "cuda_forge/condensationSourceKernels_d.cuh"
}
#include "cuda_forge/condensationUpdateLimiter_d.cuh"
#include "cuda_forge/condensationRealizability_d.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)
template <class T> static T* up(const std::vector<T>& v) { T* d = nullptr; cudaMalloc((void**)&d, v.size()*sizeof(T)); cudaMemcpy(d, v.data(), v.size()*sizeof(T), cudaMemcpyHostToDevice); return d; }
template <class T> static std::vector<T> down(const T* d, size_t n) { std::vector<T> v(n); cudaMemcpy(v.data(), d, n*sizeof(T), cudaMemcpyDeviceToHost); return v; }
static SpeciesThermo mk(double MW,double sig,double eps,const double lo[9],const double hi[9]){
    SpeciesThermo s; s.MW=MW; s.sigma_LJ=sig; s.eps_kB=eps; s.Tlo=200.0; s.Tmid=1000.0; s.Thi=6000.0; s.h_datum=0.0; s.invMW=1.0/MW;
    for(int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; } return s; }
static SpeciesThermoF toF(const SpeciesThermo& s){
    SpeciesThermoF f; f.MW=(float)s.MW; f.invMW=(float)(1.0/s.MW); f.R=(float)(THERMO_RU/s.MW);
    f.sigma_LJ=(float)s.sigma_LJ; f.eps_kB=(float)s.eps_kB; f.Tlo=(float)s.Tlo; f.Tmid=(float)s.Tmid; f.Thi=(float)s.Thi;
    for(int k=0;k<9;k++){ f.low[k]=(float)s.low[k]; f.high[k]=(float)s.high[k]; } return f; }
struct State { double T, P, ro, Yw, g, q0, q1, q2; };

// ---- (a) Δτ 不変性 ----
struct Outs { std::vector<flow_float> rr, r0, r1, r2, sg, s1, dL; };
static Outs run_source(int model, int carrier, const std::vector<State>& st, int useFloat, int limiterMode, double dtval, int evap)
{
    const int n = (int)st.size();
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(model, o);
    CondTablesHost ht; cond_tables_build_host(cp, ht); const CondTablesF tb = cond_tables_upload(ht);
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double H2Olo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double H2Ohi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    std::vector<SpeciesThermo> sp = { mk(0.0280134,3.621,97.53,N2lo,N2hi), mk(0.0180153,2.605,572.4,H2Olo,H2Ohi) };
    for (auto& s : sp) { const double hr = thermo_h_molar(s, 298.15); s.low[7] += -hr/THERMO_RU; s.high[7] += -hr/THERMO_RU; }
    std::vector<SpeciesThermoF> spf = { toF(sp[0]), toF(sp[1]) };
    const int nSp = 2, cgs = 1;
    std::vector<flow_float> T(n),P(n),ro(n),cpc(n),Rm(n),roY0(n),roY1(n),rog(n),q0(n),q1(n),q2(n),vol(n,1.0e-9f),dt(n,(flow_float)dtval);
    for (int i=0;i<n;i++){ const State& s=st[i]; T[i]=(flow_float)s.T; P[i]=(flow_float)s.P; ro[i]=(flow_float)s.ro; rog[i]=(flow_float)(s.ro*s.g);
        q0[i]=(flow_float)s.q0; q1[i]=(flow_float)s.q1; q2[i]=(flow_float)s.q2; roY1[i]=(flow_float)(s.ro*s.Yw); roY0[i]=(flow_float)(s.ro*(1.0-s.Yw));
        double Y[2]={1.0-s.Yw, s.Yw}; double c,h; thermo_cph_mix(sp.data(),2,Y,s.T,&c,&h); cpc[i]=(flow_float)c; Rm[i]=(flow_float)thermo_R_mix(sp.data(),2,Y); }
    flow_float *dT=up(T),*dP=up(P),*dro=up(ro),*dcp=up(cpc),*dRm=up(Rm),*dY0=up(roY0),*dY1=up(roY1),*dvol=up(vol),*ddt=up(dt);
    std::vector<flow_float*> hY = {dY0, dY1}; flow_float** dYall = up(hY);
    SpeciesThermo* dsp = up(sp); SpeciesThermoF* dspf = up(spf);
    auto zeros = [&](){ return up(std::vector<flow_float>(n, 0.0f)); };
    flow_float *A_rog=up(rog),*A_q0=up(q0),*A_q1=up(q1),*A_q2=up(q2),*rr=zeros(),*r0=zeros(),*r1=zeros(),*r2=zeros(),*sg=zeros(),*s0=zeros(),*s1=zeros(),*s2=zeros();
    flow_float *dS=zeros(),*dD=zeros(),*dR=zeros(),*dTs=zeros(),*dTh=zeros(),*dLm=zeros();
    const int blk = 128, grd = (n + blk - 1)/blk;
    const double Rw = cp.R, M = cp.M; const float cp_cpg = 1038.8f, gamma_cpg = 1.4f;
    if (!useFloat) {
        condensation_source_d<<<grd,blk>>>(n, model, carrier, Rw, M, 1, 0, o, carrier ? dsp : nullptr, nSp, carrier ? dYall : nullptr, carrier ? cgs : -1,
            0, 3.18, 0, evap, 1.0e-9, 0, 0.5, 0, 1.0, 5.0e-3, 10.0, cp_cpg, gamma_cpg, 1.0e35, 5.0e-3, 1.0, limiterMode,
            dvol, ddt, dT, dP, dro, carrier ? dcp : nullptr, carrier ? dRm : nullptr, carrier ? dY1 : nullptr,
            A_rog, A_q0, A_q1, A_q2, rr, r0, r1, r2, sg, s0, s1, s2, dS, dD, dR, dTs, dTh, dLm);
    } else {
        CondDoubleArgs dbl; dbl.opts = o; dbl.sp = carrier ? dsp : nullptr; dbl.condModel = model; dbl.Rw = Rw; dbl.M = M; dbl.twoTemp = 0;
        dbl.gyarC = 3.18; dbl.evapRmin = 1.0e-9; dbl.evapLamMin = 0.5; dbl.Jmax = 1.0e35; dbl.dg_max = 5.0e-3; dbl.dT_max = 1.0; dbl.cprops = cp;
        condensation_source_f_d<<<grd,blk>>>(n, carrier, (float)Rw, 1, 0, condProps_to_f(cp), tb, 0.0f, dbl, carrier ? dspf : nullptr, nSp, carrier ? dYall : nullptr, carrier ? cgs : -1,
            0, 3.18f, evap, 1.0e-9f, 0, 0.5f, cp_cpg, gamma_cpg, 5.0e-3f, 1.0f, limiterMode,
            dvol, ddt, dT, dP, dro, carrier ? dcp : nullptr, carrier ? dRm : nullptr, carrier ? dY1 : nullptr,
            A_rog, A_q0, A_q1, A_q2, rr, r0, r1, r2, sg, s0, s1, s2, dS, dD, dR, dTs, dTh, dLm);
    }
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; }
    Outs out; out.rr=down(rr,n); out.r0=down(r0,n); out.r1=down(r1,n); out.r2=down(r2,n); out.sg=down(sg,n); out.s1=down(s1,n); out.dL=down(dLm,n);
    return out;
}
static void test_dt_invariance(const char* name, int model, int carrier, const std::vector<State>& st)
{
    for (int useFloat : {0, 1}) {
        Outs a = run_source(model, carrier, st, useFloat, 1, 1.0e-7, 1);
        Outs b = run_source(model, carrier, st, useFloat, 1, 1.0e-3, 1);
        int nd = 0, nAct = 0, nNF = 0;
        for (size_t i = 0; i < st.size(); ++i) {
            const bool act = (a.rr[i] != 0.0f || a.r0[i] != 0.0f || a.r1[i] != 0.0f || a.r2[i] != 0.0f);
            if (act) ++nAct;
            if (a.rr[i]!=b.rr[i] || a.r0[i]!=b.r0[i] || a.r1[i]!=b.r1[i] || a.r2[i]!=b.r2[i] || a.sg[i]!=b.sg[i] || a.s1[i]!=b.s1[i]) ++nd;
            if (!std::isfinite((double)a.rr[i]) || !std::isfinite((double)a.sg[i])) ++nNF;
        }
        printf("  [%s %s] mode 1: dt 1e-7 vs 1e-3 -> %d/%zu cells differ (active sources %d), non-finite %d\n", name, useFloat ? "float" : "double", nd, st.size(), nAct, nNF);
        CHECK(nd == 0, "%s %s: mode 1 residual depends on dt_local", name, useFloat ? "float" : "double");
        CHECK(nAct > 0, "%s %s: no active source cells (test grid too weak)", name, useFloat ? "float" : "double");
        CHECK(nNF == 0, "%s %s: non-finite outputs", name, useFloat ? "float" : "double");
        // 旧 mode 0 は dt に依存する (旧挙動の確認; 大 dt で θ<1 のセルが存在)
        Outs c = run_source(model, carrier, st, useFloat, 0, 1.0e-7, 1);
        Outs d = run_source(model, carrier, st, useFloat, 0, 1.0e-3, 1);
        int nd0 = 0; for (size_t i = 0; i < st.size(); ++i) if (c.rr[i] != d.rr[i]) ++nd0;
        printf("  [%s %s] mode 0 (legacy): %d cells differ between dt 1e-7 and 1e-3 (expected > 0)\n", name, useFloat ? "float" : "double", nd0);
        CHECK(nd0 > 0, "%s %s: legacy mode 0 shows no dt dependence (test not discriminating)", name, useFloat ? "float" : "double");
    }
}

// ---- (b)-(f) 更新クランプ kernel ----
struct Cell { double dt, vol, ro, roYw, T, Ng, NQ2, NQ1, NQ0, rg, rQ2, rQ1, rQ0; };
static void run_limiter(const std::vector<Cell>& cs, int carrier, std::vector<double>& theta, std::vector<double>& corr, std::vector<double>& corrQ_out, std::vector<std::array<double,4>>& outv, std::vector<std::array<double,4>>& delta)
{
    const int n = (int)cs.size();
    std::vector<flow_float> dt(n),vol(n),ro(n),roYw(n),T(n),Ng(n),NQ2(n),NQ1(n),NQ0(n),rg(n),rQ2(n),rQ1(n),rQ0(n),z(n,0.0f);
    for (int i=0;i<n;i++){ const Cell& c=cs[i]; dt[i]=c.dt; vol[i]=c.vol; ro[i]=c.ro; roYw[i]=c.roYw; T[i]=c.T; Ng[i]=c.Ng; NQ2[i]=c.NQ2; NQ1[i]=c.NQ1; NQ0[i]=c.NQ0; rg[i]=c.rg; rQ2[i]=c.rQ2; rQ1[i]=c.rQ1; rQ0[i]=c.rQ0; }
    flow_float *ddt=up(dt),*dvol=up(vol),*dro=up(ro),*dYw=up(roYw),*dT=up(T),*dNg=up(Ng),*dNQ2=up(NQ2),*dNQ1=up(NQ1),*dNQ0=up(NQ0);
    flow_float *drg=up(rg),*drQ2=up(rQ2),*drQ1=up(rQ1),*drQ0=up(rQ0);
    auto zeros=[&](){ return up(z); };
    flow_float *sjg=zeros(),*sj2=zeros(),*sj1=zeros(),*sj0=zeros(),*tdg=zeros(),*td2=zeros(),*td1=zeros(),*td0=zeros();
    flow_float *og=zeros(),*o2=zeros(),*o1=zeros(),*o0=zeros(),*lim=zeros(),*cor=zeros(),*corQ=zeros();
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    cond_moment_update_limited_d<<<(n+127)/128,128>>>(n, ddt, dvol, dro, carrier ? dYw : nullptr, 0.0, dT, nullptr, nullptr, 1220.7f, 1.315f,
        COND_MODEL_H2O, o, 5.0e-3, 1.0, 0.5, dNg, dNQ2, dNQ1, dNQ0, drg, drQ2, drQ1, drQ0, sjg, sj2, sj1, sj0, tdg, td2, td1, td0, og, o2, o1, o0, lim, cor, corQ);
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; }
    auto L=down(lim,n), C=down(cor,n), CQ=down(corQ,n), G=down(og,n), Q2=down(o2,n), Q1=down(o1,n), Q0=down(o0,n);
    theta.resize(n); corr.resize(n); corrQ_out.resize(n); outv.resize(n); delta.resize(n);
    for (int i=0;i<n;i++){ theta[i]=L[i]; corr[i]=C[i]; corrQ_out[i]=CQ[i]; outv[i]={G[i],Q2[i],Q1[i],Q0[i]};
        const Cell& c=cs[i]; const double f=c.dt/c.vol; delta[i]={c.rg*f, c.rQ2*f, c.rQ1*f, c.rQ0*f}; }   // sj=td=0 → δ = res dt/V
}
static void test_update_limiter()
{
    // 空気+水蒸気 (ρ 0.1, Y_w 0.0377, T 240 K)。res は [kg/m3/s]×V 単位: δ = res dt/V。
    const double ro=0.1, Yw=0.0377, T=240.0, dt=1.0e-5, vol=1.0e-6;
    std::vector<Cell> cs;
    // (b) 上限内: Δg = 1e-4 (< dg_max 5e-3, ΔT = 1e-4·2.5e6/~700 ≈ 0.36 K < 1 K)
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*1e-3, 1e-2, 1e2, 1e12,  ro*1e-4*vol/dt, 1e-3*vol/dt, 1e1*vol/dt, 1e11*vol/dt});
    // (c) 潜熱超過: Δg = 2e-3 → ΔT ≈ 7 K > 1 K → θ ≈ 1/7 (dg は 5e-3 未満)
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*1e-3, 1e-2, 1e2, 1e12,  ro*2e-3*vol/dt, 1e-3*vol/dt, 1e1*vol/dt, 1e11*vol/dt});
    // (d) 残差 0
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*1e-3, 1e-2, 1e2, 1e12,  0,0,0,0});
    // (f1) 蒸気枯渇: g_old = Yw − 5e-5, Δg 候補 = 1e-4 > avail 5e-5 → θ = 0.5 (ΔT 0.36 K は非作動)
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*(Yw-5e-5), 1e-2, 1e2, 1e12,  ro*1e-4*vol/dt, 1e-3*vol/dt, 1e1*vol/dt, 1e11*vol/dt});
    // (f2) 負増分の半径半減上限: g_old 1e-3, Δg = −9.5e-4 (> 0.875 g_old = 8.75e-4) → θ = 8.75e-4/9.5e-4 (dg/dT は非作動: ΔT 3.4 K! → dT_max/ΔT = 0.29 が先に効く)
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*1e-3, 1e-2, 1e2, 1e12,  -ro*9.5e-4*vol/dt, -1e-3*vol/dt, -1e1*vol/dt, 0});
    // (f3) 負増分 (小): g_old 1e-4, Δg = −9.5e-5 (> 0.875e-4) → θ = 0.875e-4/9.5e-5 (ΔT 0.34 K, dg 非作動)
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*1e-4, 1e-2, 1e2, 1e12,  -ro*9.5e-5*vol/dt, -1e-3*vol/dt, -1e1*vol/dt, 0});
    // (f4) g_old = 0 への負増分 (輸送の負流入): θ=1 (dg 1e-5 は上限内)、floor で 0、補正量 = 1e-5
    cs.push_back({dt,vol,ro,ro*Yw,T, 0.0, 0.0, 0.0, 0.0,  -ro*1e-5*vol/dt, 0,0,0});
    // (i) 密度が半分になった後の更新 (ρ_new=0.05, N は ρ_old=0.1 の保存量): δ_g=0 → Δg=0 → θ=1 (Δg は更新済み密度で測る相変化分のみ; plan §4.2-4)
    cs.push_back({dt,vol,0.05,0.05*Yw,T, ro*1e-3, 1e-2, 1e2, 1e12,  0,0,0,0});
    // (f5) float 極小増分
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*1e-3, 1e-2, 1e2, 1e12,  ro*1e-30*vol/dt, 1e-30, 1e-30, 1e-30});
    auto near = [](double a, double b, double rtol){ return std::fabs(a-b) <= rtol*std::max(std::fabs(a), std::fabs(b)) + 1e-30; };
    // (j) Q の負値 floor: Q1 に −大 の増分 → out_Q1=0, corrQ=1 (100 % 補正)。g は上限内で θ=1
    cs.push_back({dt,vol,ro,ro*Yw,T, ro*1e-3, 1e-2, 1e2, 1e12,  0, 0, -2e2*vol/dt, 0});
    std::vector<double> th, corr, corrQ; std::vector<std::array<double,4>> out, del;
    run_limiter(cs, 1, th, corr, corrQ, out, del);
    for (int i = 0; i < 9; ++i) CHECK(corrQ[i] == 0.0, "case %d: corrQ should be 0 (no Q floor), got %g", i, corrQ[i]);
    printf("  (j) negative Q1 increment: theta=%.4f out_Q1=%.3g corrQ=%.3g (expected 1)\n", th[9], out[9][2], corrQ[9]);
    CHECK(near(th[9],1.0,1e-6) && out[9][2] == 0.0 && corrQ[9] == 1.0, "(j) Q floor not recorded as 100 %% correction");
    // (b)
    printf("  (b) theta=%.6f corr=%.3g\n", th[0], corr[0]);
    CHECK(near(th[0],1.0,1e-6), "(b) theta != 1");
    for (int k=0;k<4;k++) CHECK(near(out[0][k]-(double)(flow_float)( (k==0? cs[0].Ng : k==1? cs[0].NQ2 : k==2? cs[0].NQ1 : cs[0].NQ0)), del[0][k], 1e-4), "(b) increment k=%d changed", k);
    // (c)
    { const double Ldc = 2.5e6/(1220.7-285.0); const double dg = 2e-3; const double exp_th = 1.0/(dg*Ldc); printf("  (c) theta=%.4f (expected ~%.4f), corr=%.3g\n", th[1], exp_th, corr[1]);
      CHECK(th[1] < 0.3 && th[1] > 0.05, "(c) theta not in latent-limited range");
      const double Nk[4]={cs[1].Ng,cs[1].NQ2,cs[1].NQ1,cs[1].NQ0}; double ratio[4];
      for (int k=0;k<4;k++) ratio[k]=(out[1][k]-(double)(flow_float)Nk[k])/del[1][k];
      for (int k=1;k<4;k++) CHECK(near(ratio[k], ratio[0], 1e-3), "(c) moment %d not scaled by the same theta (%.5f vs %.5f)", k, ratio[k], ratio[0]);
      CHECK(near(ratio[0], th[1], 1e-3), "(c) g increment ratio %.5f != theta %.5f", ratio[0], th[1]); }
    // (d)
    printf("  (d) theta=%.6f out_g=%.6g N_g=%.6g\n", th[2], out[2][0], cs[2].Ng);
    CHECK(near(th[2],1.0,1e-6) && near(out[2][0], (double)(flow_float)cs[2].Ng, 1e-7), "(d) zero residual is not a no-op");
    // (f1)
    printf("  (f1) theta=%.6f (expected 0.5)\n", th[3]); CHECK(near(th[3],0.5,2e-2), "(f1) vapor exhaustion theta != avail/dg");
    CHECK(out[3][0] <= (double)(flow_float)(ro*Yw)*(1.0+1e-6), "(f1) g exceeds available vapour");
    // (f2): ΔT = 9.5e-4·L/c_v,eff ≈ 3.4 K → dT 律速 (0.29) が半径半減 (0.92) より強い
    printf("  (f2) theta=%.4f\n", th[4]); CHECK(th[4] < 0.5 && th[4] > 0.1, "(f2) negative-increment theta not latent-limited as expected");
    // (f3)
    printf("  (f3) theta=%.4f (expected 0.875e-4/9.5e-5 = %.4f)\n", th[5], 0.875e-4/9.5e-5); CHECK(near(th[5], 0.875e-4/9.5e-5, 2e-2), "(f3) evaporation radius-halving bound wrong");
    // (f4)
    printf("  (f4) theta=%.4f out_g=%.3g corr=%.3g (expected 1e-5)\n", th[6], out[6][0], corr[6]);
    CHECK(near(th[6],1.0,1e-6) && out[6][0] == 0.0 && near(corr[6], 1e-5, 1e-2), "(f4) floor/correction wrong");
    // (i)
    printf("  (i) density halved, zero increment: theta=%.4f out_g=%.6g (N_g %.6g)\n", th[7], out[7][0], cs[7].Ng);
    CHECK(near(th[7],1.0,1e-6) && near(out[7][0], (double)(flow_float)cs[7].Ng, 1e-7), "(i) density change with zero increment must be a no-op (theta 1)");
    // (f5)
    printf("  (f5) theta=%.4f out finite=%d\n", th[8], std::isfinite(out[8][0]) && std::isfinite(out[8][3]));
    CHECK(near(th[8],1.0,1e-6) && std::isfinite(out[8][0]) && std::isfinite(out[8][3]), "(f5) tiny increment not finite / theta != 1");
}


__global__ void evap_rate_f_kernel(CondSpeciesPropsF cp, CondTablesF tb, float T, float pv, float rod, float g, float q0, float q1, float q2, float* out)
{
    float S0,S1,S2,Sg,r30,drdt;
    cond_evap_source_rate_f(cp, tb, T, pv, rod, g, q0, q1, q2, 1.0e-9f, 0, pv, 3.18f, 0, &S0,&S1,&S2,&Sg,&r30,&drdt);
    out[0]=S0; out[1]=S1; out[2]=S2; out[3]=Sg; out[4]=r30; out[5]=drdt;
}
// ---- (g) 蒸発ソース (一様 ṙ 形) の値そのものを多分散モーメントで検証 (codex result M4) ----
static void test_evap_source_values()
{
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o);
    const double T = 250.0, rod = 0.1, Yw = 0.0377, S = 0.5;
    const double pv = S*cond_psat(cp, T), rho_l = cond_rho_cond(cp, T);
    // 半径 r と 2r を同数 (N/2 ずつ): Q0=N, Q1=N·1.5r, Q2=N·2.5r², g = (4/3)πρ_l N (r³+8r³)/2 /ρ
    const double N = 1.0e15, r = 3.0e-8;
    const double q0 = N, q1 = N*1.5*r, q2 = N*2.5*r*r;
    const double g  = (4.0/3.0)*COND_PI*rho_l*N*4.5*r*r*r/rod;
    double S0,S1,S2,Sg,r30,drdt;
    cond_evap_source_rate(cp, T, pv, rod, g, q0, q1, q2, 1.0e-9, 0, pv, 3.18, 0, &S0,&S1,&S2,&Sg,&r30,&drdt);
    const double r30_exp = cbrt(4.5)*r;
    const double drdt_exp = cond_evap_rate(cp, T, pv, r30_exp, 0, pv, 3.18, 0);
    auto near = [](double a, double b, double rtol){ return std::fabs(a-b) <= rtol*std::max(std::fabs(a), std::fabs(b)); };
    printf("  (g) r30 %.4e (exp %.4e) drdt %.4e  S_Q1/(q0 drdt) %.6f  S_Q2/(2 q1 drdt) %.6f  S_g/(4πρl q2 drdt) %.6f  S_Q0 %g\n",
           r30, r30_exp, drdt, S1/(q0*drdt), S2/(2.0*q1*drdt), Sg/(4.0*COND_PI*rho_l*q2*drdt), S0);
    CHECK(near(r30, r30_exp, 1e-12) && drdt < 0.0 && near(drdt, drdt_exp, 1e-12), "(g) r30/drdt mismatch");
    CHECK(near(S1, q0*drdt, 1e-12) && near(S2, 2.0*q1*drdt, 1e-12) && near(Sg, 4.0*COND_PI*rho_l*q2*drdt, 1e-12) && S0 == 0.0, "(g) evaporation source is not the uniform-rdot moment form");
    // 旧 λ スケール極限 (a=ṙ/r30: 3aρg) との差が多分散では 0 でないことを記録 (monodisperse では一致)
    const double a = drdt/r30; printf("  (g) polydisperse: S_g(uniform rdot)/S_g(self-similar) = %.4f (monodisperse would be 1)\n", Sg/(3.0*a*rod*g));
    // float 版も同式 (物性表は device メモリなので 1 スレッド kernel で評価する)
    CondTablesHost ht; cond_tables_build_host(cp, ht); const CondTablesF tb = cond_tables_upload(ht);
    float* dout = up(std::vector<float>(6, 0.0f));
    evap_rate_f_kernel<<<1,1>>>(condProps_to_f(cp), tb, (float)T, (float)pv, (float)rod, (float)g, (float)q0, (float)q1, (float)q2, dout);
    cudaDeviceSynchronize();
    auto fo = down(dout, 6); const float f1 = fo[1], f2 = fo[2], fg = fo[3];
    printf("  (g) float: S_Q1 rel %.2e S_Q2 rel %.2e S_g rel %.2e\n", std::fabs(f1-S1)/std::fabs(S1), std::fabs(f2-S2)/std::fabs(S2), std::fabs(fg-Sg)/std::fabs(Sg));
    CHECK(std::fabs(f1-S1)/std::fabs(S1) < 2e-3 && std::fabs(fg-Sg)/std::fabs(Sg) < 2e-3, "(g) float evaporation source differs from double");
}

// ---- (h) 輸送とソースが非ゼロで釣り合う 1 セル固定点: Δτ を変えても固定点が同じ (codex result M7) ----
// 1 セル: 残差 = 流入 β(ρφ_in − ρφ)·V + 実ソース kernel (condensation_source_d, mode 1) → 更新クランプで反復。
// 上流は dry (φ_in=0) なので凝縮域では流入希釈と成長が釣り合う非自明な固定点になる。
static void test_one_cell_fixed_point()
{
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o);
    CondTablesHost ht; cond_tables_build_host(cp, ht); const CondTablesF tb = cond_tables_upload(ht); (void)tb;
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double H2Olo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double H2Ohi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    std::vector<SpeciesThermo> sp = { mk(0.0280134,3.621,97.53,N2lo,N2hi), mk(0.0180153,2.605,572.4,H2Olo,H2Ohi) };
    for (auto& s : sp) { const double hr = thermo_h_molar(s, 298.15); s.low[7] += -hr/THERMO_RU; s.high[7] += -hr/THERMO_RU; }
    SpeciesThermo* dsp = up(sp);
    printf("  (h) setup ok\n"); fflush(stdout);
    const double T = 232.0, Yw = 0.0377, S = 20.0, Rw = cp.R, Rmix = 285.0;
    const double pv = S*cond_psat(cp, T), rod = pv/(Yw*Rw*T);
    double Y[2]={1.0-Yw, Yw}; double cpc,h; thermo_cph_mix(sp.data(),2,Y,T,&cpc,&h); const double Rm = thermo_R_mix(sp.data(),2,Y);
    const double V = 1.0e-6, beta = 2000.0;   // 流入率 [1/s] (滞留時間 0.5 ms)
    // dt は 1e-6 / 1e-5 / 1e-4 (100 倍): flow_float=float の反復は増分が値の ~1e-7 倍を割ると停滞するので、極端に小さい dt (1e-7) は
    // 残差が残ったまま止まる (相対 ~3e-4)。固定点の Δτ 非依存はこの停滞床より大きい範囲で見る。
    const double dts[3] = {1.0e-6, 1.0e-5, 1.0e-4};
    std::vector<double> fixed[3]; double thetas[3];
    for (int k = 0; k < 3; ++k) {
        const double dt = dts[k];
        std::vector<flow_float> one(1);
        auto arr = [&](double v){ return up(std::vector<flow_float>(1,(flow_float)v)); };
        flow_float *dT=arr(T),*dP=arr(rod*Rm*T),*dro=arr(rod),*dcp=arr(cpc),*dRm=arr(Rm),*dY0=arr(rod*(1-Yw)),*dY1=arr(rod*Yw),*dvol=arr(V),*ddt=arr(dt);
        std::vector<flow_float*> hY={dY0,dY1}; flow_float** dYall=up(hY);
        flow_float *rog=arr(0),*q0=arr(0),*q1=arr(0),*q2=arr(0),*Ng=arr(0),*NQ0=arr(0),*NQ1=arr(0),*NQ2=arr(0);
        flow_float *rr=arr(0),*r0=arr(0),*r1=arr(0),*r2=arr(0),*sg=arr(0),*s0=arr(0),*s1=arr(0),*s2=arr(0),*td=arr(0);
        flow_float *dS=arr(0),*dD=arr(0),*dR=arr(0),*dTs=arr(0),*dTh=arr(0),*dLm=arr(0),*cG=arr(0),*cQ=arr(0);
        const double rawres_tol = 1e-4;   // 固定点判定: |res_k|/(β V |ρφ_k|) = 輸送速度で正規化した残差 (Δτ 非依存; float の停滞床 ~1e-5)
        double rawres_last = 0.0;
        printf("  (h) dt %.0e arrays allocated\n", dt); fflush(stdout);
        double prev[4]={0,0,0,0}; int it=0; double chg=1.0;
        for (it = 0; it < 400000 && chg > 1e-13; ++it) {
            // 残差: 流入 (dry 上流) β(0 − ρφ)V + ソース
            std::vector<flow_float> cur(4); cudaMemcpy(&cur[0], rog, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&cur[1], q0, sizeof(flow_float), cudaMemcpyDeviceToHost);
            cudaMemcpy(&cur[2], q1, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&cur[3], q2, sizeof(flow_float), cudaMemcpyDeviceToHost);
            flow_float rt[4]; for (int m=0;m<4;m++) rt[m] = (flow_float)(-beta*(double)cur[m]*V);
            cudaMemcpy(rr, &rt[0], sizeof(flow_float), cudaMemcpyHostToDevice); cudaMemcpy(r0, &rt[1], sizeof(flow_float), cudaMemcpyHostToDevice);
            cudaMemcpy(r1, &rt[2], sizeof(flow_float), cudaMemcpyHostToDevice); cudaMemcpy(r2, &rt[3], sizeof(flow_float), cudaMemcpyHostToDevice);
            flow_float tdv = (flow_float)(beta*V); cudaMemcpy(td, &tdv, sizeof(flow_float), cudaMemcpyHostToDevice);   // transport_diag = βV (陰的流出)
            cudaMemcpy(Ng, rog, sizeof(flow_float), cudaMemcpyDeviceToDevice); cudaMemcpy(NQ0, q0, sizeof(flow_float), cudaMemcpyDeviceToDevice);
            cudaMemcpy(NQ1, q1, sizeof(flow_float), cudaMemcpyDeviceToDevice); cudaMemcpy(NQ2, q2, sizeof(flow_float), cudaMemcpyDeviceToDevice);
            condensation_source_d<<<1,1>>>(1, COND_MODEL_H2O, 1, Rw, cp.M, 1, 0, o, dsp, 2, dYall, 1,
                0, 3.18, 0, 1, 1.0e-9, 0, 0.5, 0, 1.0, 5.0e-3, 10.0, (float)cpc, 1.315f, 1.0e35, 5.0e-3, 1.0, 1,
                dvol, ddt, dT, dP, dro, dcp, dRm, dY1, rog, q0, q1, q2, rr, r0, r1, r2, sg, s0, s1, s2, dS, dD, dR, dTs, dTh, dLm);
            cond_moment_update_limited_d<<<1,1>>>(1, ddt, dvol, dro, dY1, 0.0, dT, dcp, dRm, (float)cpc, 1.315f, COND_MODEL_H2O, o, 5.0e-3, 1.0, 0.5,
                Ng, NQ2, NQ1, NQ0, rr, r2, r1, r0, sg, s2, s1, s0, td, td, td, td, rog, q2, q1, q0, dLm, cG, cQ);
            // 実際の更新全経路: 更新クランプ → 実現可能性クランプ (g<=Y_w, 負値, 消滅) → 次の残差
            cond_realizability_clamp_d<<<1,1>>>(1, dro, dY1, rog, q0, q1, q2, 1, COND_MODEL_H2O, Rw, 1.0e-9, 5.0e-7, dT, dP, o, cG, cQ, nullptr, nullptr, nullptr, nullptr, 1);
            { cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("  (h) CUDA error at it %d: %s\n", it, cudaGetErrorString(e)); fflush(stdout); ++g_fail; break; } }
            { std::vector<flow_float> rs(4); cudaMemcpy(&rs[0], rr, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&rs[1], r0, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&rs[2], r1, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&rs[3], r2, sizeof(flow_float), cudaMemcpyDeviceToHost);
              std::vector<flow_float> st(4); cudaMemcpy(&st[0], rog, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&st[1], q0, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&st[2], q1, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&st[3], q2, sizeof(flow_float), cudaMemcpyDeviceToHost);
              rawres_last = 0.0; for (int m=0;m<4;m++) rawres_last = std::max(rawres_last, std::fabs((double)rs[m])/(beta*V*std::max(std::fabs((double)st[m]),1e-30))); }   // Δτ を含まない正規化: |res|/(β V |ρφ|) (codex result-3 m4)
            std::vector<flow_float> nw(4); cudaMemcpy(&nw[0], rog, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&nw[1], q0, sizeof(flow_float), cudaMemcpyDeviceToHost);
            cudaMemcpy(&nw[2], q1, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&nw[3], q2, sizeof(flow_float), cudaMemcpyDeviceToHost);
            chg = 0.0; for (int m=0;m<4;m++){ const double d = std::fabs((double)nw[m]-prev[m])/std::max(std::fabs((double)nw[m]),1e-30); chg = std::max(chg,d); prev[m]=nw[m]; }
            if (it % 50000 == 0 || it < 3) { float th; cudaMemcpy(&th, dLm, sizeof(float), cudaMemcpyDeviceToHost); thetas[k]=th; }
        }
        fixed[k].assign(prev, prev+4);
        float th, cg, cq; cudaMemcpy(&th, dLm, sizeof(float), cudaMemcpyDeviceToHost); cudaMemcpy(&cg, cG, sizeof(float), cudaMemcpyDeviceToHost); cudaMemcpy(&cq, cQ, sizeof(float), cudaMemcpyDeviceToHost);
        printf("  (h) dt %.0e: %d iterations, fixed point g=%.6e Q0=%.6e Q1=%.6e Q2=%.6e, theta %.4f, raw residual (rel update/step) %.2e, clamp corr G %.2e Q %.2e\n", dt, it, prev[0]/rod, prev[1], prev[2], prev[3], th, rawres_last, cg, cq);
        CHECK(rawres_last < rawres_tol, "(h) raw residual not at the fixed point (%.2e)", rawres_last);
        CHECK(cg == 0.0f && cq == 0.0f, "(h) hard clamps active at the fixed point (G %g Q %g)", cg, cq);
        CHECK(prev[0] > 0.0 && prev[1] > 0.0, "(h) trivial fixed point (no condensation) — test state not condensing");
        CHECK(th > 0.999f, "(h) theta_u not 1 at the fixed point");
    }
    double worst = 0.0; for (int a=0;a<3;a++) for (int b=a+1;b<3;b++) for (int m=0;m<4;m++) worst = std::max(worst, std::fabs(fixed[a][m]-fixed[b][m])/std::max(std::fabs(fixed[b][m]),1e-30));
    printf("  (h) fixed-point max pairwise relative difference over dt 1e-6/1e-5/1e-4: %.3e\n", worst);
    CHECK(worst < 2e-4, "(h) fixed point depends on dt (%.3e)", worst);
}


// ---- (k) 小液滴の反例 (codex result-3 M1): 整合した単分散 g=1e-5, r30=1.5e-9 (< 2 r_min=2e-9, g>g_rm=5e-7) が、ソース→更新→消滅を
//      通して蒸発し続け、g≤g_rm で消滅クランプが確定すること (旧 r30<2 r_min→S=0 では止まっていた)。
static void test_small_droplet_evaporates()
{
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o);
    CondTablesHost ht; cond_tables_build_host(cp, ht); const CondTablesF tb = cond_tables_upload(ht);
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double H2Olo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double H2Ohi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    std::vector<SpeciesThermo> sp = { mk(0.0280134,3.621,97.53,N2lo,N2hi), mk(0.0180153,2.605,572.4,H2Olo,H2Ohi) };
    for (auto& s : sp) { const double hr = thermo_h_molar(s, 298.15); s.low[7] += -hr/THERMO_RU; s.high[7] += -hr/THERMO_RU; }
    std::vector<SpeciesThermoF> spf = { toF(sp[0]), toF(sp[1]) };
    SpeciesThermo* dsp = up(sp); SpeciesThermoF* dspf = up(spf);
    const double T = 250.0, Yw = 0.0377, S = 0.5, Rw = cp.R, g0 = 1.0e-5;
    const double pv = S*cond_psat(cp, T), rho_l = cond_rho_cond(cp, T);
    const double rod = pv/((Yw - g0)*Rw*T);   // kernel は p_v=ρ(Y_w−g)R_wT で蒸気状態を作るので、S=0.5 になる ρ を選ぶ
    double Y[2]={1.0-Yw, Yw}; double cpc,h; thermo_cph_mix(sp.data(),2,Y,T,&cpc,&h); const double Rm = thermo_R_mix(sp.data(),2,Y);
    const double V = 1.0e-6, dt = 1.0e-6;
    // 反例 2 種: (A) 整合した単分散 r30=1.5e-9 (< 2 r_min, g>g_rm)、(B) Q0=Q1=Q2=0 で g=1e-5 (不整合; codex result-4 M1)
    struct Case { const char* name; double r30; };
    for (const Case cs : { Case{"r30=1.5e-9 (consistent)", 1.5e-9}, Case{"Q0=Q1=Q2=0 (inconsistent)", 0.0} })
    for (int useFloat : {0, 1}) {
        const double q0 = (cs.r30 > 0.0) ? g0*rod/((4.0/3.0)*COND_PI*rho_l*cs.r30*cs.r30*cs.r30) : 0.0, q1 = q0*cs.r30, q2 = q0*cs.r30*cs.r30;
        auto arr = [&](double v){ return up(std::vector<flow_float>(1,(flow_float)v)); };
        flow_float *dT=arr(T),*dP=arr(rod*Rm*T),*dro=arr(rod),*dcp=arr(cpc),*dRm=arr(Rm),*dY0=arr(rod*(1-Yw)),*dY1=arr(rod*Yw),*dvol=arr(V),*ddt=arr(dt);
        std::vector<flow_float*> hY={dY0,dY1}; flow_float** dYall=up(hY);
        flow_float *rog=arr(rod*g0),*Q0=arr(q0),*Q1=arr(q1),*Q2=arr(q2),*Ng=arr(0),*NQ0=arr(0),*NQ1=arr(0),*NQ2=arr(0);
        flow_float *rr=arr(0),*r0=arr(0),*r1=arr(0),*r2=arr(0),*sg=arr(0),*s0=arr(0),*s1=arr(0),*s2=arr(0),*td=arr(0);
        flow_float *dS=arr(0),*dD=arr(0),*dR=arr(0),*dTs=arr(0),*dTh=arr(0),*dLm=arr(0),*cG=arr(0),*cQ=arr(0);
        CondDoubleArgs dbl; dbl.opts = o; dbl.sp = dsp; dbl.condModel = COND_MODEL_H2O; dbl.Rw = Rw; dbl.M = cp.M; dbl.twoTemp = 0;
        dbl.gyarC = 3.18; dbl.evapRmin = 1.0e-9; dbl.evapLamMin = 0.5; dbl.Jmax = 1.0e35; dbl.dg_max = 5.0e-3; dbl.dT_max = 1.0; dbl.cprops = cp;
        int it = 0; float g = (float)g0, sgv = 0.0f, drdt0 = 0.0f, qs[3] = {1,1,1};
        for (it = 0; it < 200000; ++it) {
            flow_float z = 0.0f; for (flow_float* p : {rr,r0,r1,r2}) cudaMemcpy(p, &z, sizeof(flow_float), cudaMemcpyHostToDevice);
            cudaMemcpy(Ng, rog, sizeof(flow_float), cudaMemcpyDeviceToDevice); cudaMemcpy(NQ0, Q0, sizeof(flow_float), cudaMemcpyDeviceToDevice);
            cudaMemcpy(NQ1, Q1, sizeof(flow_float), cudaMemcpyDeviceToDevice); cudaMemcpy(NQ2, Q2, sizeof(flow_float), cudaMemcpyDeviceToDevice);
            if (useFloat)
                condensation_source_f_d<<<1,1>>>(1, 1, (float)Rw, 1, 0, condProps_to_f(cp), tb, 0.0f, dbl, dspf, 2, dYall, 1,
                    0, 3.18f, 1, 1.0e-9f, 0, 0.5f, (float)cpc, 1.315f, 5.0e-3f, 1.0f, 1,
                    dvol, ddt, dT, dP, dro, dcp, dRm, dY1, rog, Q0, Q1, Q2, rr, r0, r1, r2, sg, s0, s1, s2, dS, dD, dR, dTs, dTh, dLm);
            else
                condensation_source_d<<<1,1>>>(1, COND_MODEL_H2O, 1, Rw, cp.M, 1, 0, o, dsp, 2, dYall, 1,
                    0, 3.18, 0, 1, 1.0e-9, 0, 0.5, 0, 1.0, 5.0e-3, 10.0, (float)cpc, 1.315f, 1.0e35, 5.0e-3, 1.0, 1,
                    dvol, ddt, dT, dP, dro, dcp, dRm, dY1, rog, Q0, Q1, Q2, rr, r0, r1, r2, sg, s0, s1, s2, dS, dD, dR, dTs, dTh, dLm);
            cond_moment_update_limited_d<<<1,1>>>(1, ddt, dvol, dro, dY1, 0.0, dT, dcp, dRm, (float)cpc, 1.315f, COND_MODEL_H2O, o, 5.0e-3, 1.0, 0.5,
                Ng, NQ2, NQ1, NQ0, rr, r2, r1, r0, sg, s2, s1, s0, td, td, td, td, rog, Q2, Q1, Q0, dLm, cG, cQ);
            if (useFloat)
                cond_realizability_clamp_f_d<<<1,1>>>(1, dro, dY1, rog, Q0, Q1, Q2, 1, (float)Rw, 1.0e-9f, 5.0e-7f, 0.0f, dT, dP, tb, cp, cG, cQ, nullptr, nullptr, nullptr, nullptr, 1);
            else
                cond_realizability_clamp_d<<<1,1>>>(1, dro, dY1, rog, Q0, Q1, Q2, 1, COND_MODEL_H2O, Rw, 1.0e-9, 5.0e-7, dT, dP, o, cG, cQ, nullptr, nullptr, nullptr, nullptr, 1);
            cudaDeviceSynchronize();
            cudaMemcpy(&g, rog, sizeof(flow_float), cudaMemcpyDeviceToHost); g /= (float)rod;
            if (it == 0) { cudaMemcpy(&sgv, rr, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&drdt0, dD, sizeof(flow_float), cudaMemcpyDeviceToHost); }
            if (g <= 0.0f) { cudaMemcpy(&qs[0], Q0, sizeof(float), cudaMemcpyDeviceToHost); cudaMemcpy(&qs[1], Q1, sizeof(float), cudaMemcpyDeviceToHost); cudaMemcpy(&qs[2], Q2, sizeof(float), cudaMemcpyDeviceToHost); break; }
        }
        printf("  (k) %s %s: step-0 S_g/V=%.3e (drdt %.3e), liquid gone after %d steps (g=%.2e, Q0/Q1/Q2 after removal %g/%g/%g)\n", cs.name, useFloat ? "float " : "double", sgv/V, drdt0, it+1, g, qs[0], qs[1], qs[2]);
        CHECK(sgv < 0.0f && drdt0 < 0.0f, "(k) %s %s: liquid does not evaporate (source %g, drdt %g)", cs.name, useFloat ? "float" : "double", sgv, drdt0);
        CHECK(g <= 0.0f && it < 200000 && qs[0] == 0.0f && qs[1] == 0.0f && qs[2] == 0.0f, "(k) %s %s: liquid/moments never removed (g=%g after %d steps)", cs.name, useFloat ? "float" : "double", g, it);
    }
}
// ---- (l) 実現可能性の最小補正 (plan species-passive-scalar-unification §4.7 v4, codex plan-5 M3/M4): 許容領域 x≤1, x²≤y≤√x への境界クランプ。
static void test_realizability_projection()
{
    printf("[l] moment realizability minimal correction (x=Q1/(Q0 r), y=Q2/(Q0 r^2), r=(Q3/Q0)^(1/3))\n");
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o);
    const double T = 250.0, rho_l = cond_rho_cond(cp, T), ro = 1.0;
    const double Q0 = 1.0e14, r = 5.0e-8;
    const double g = (4.0/3.0)*COND_PI*rho_l*Q0*r*r*r;   // Q3 = Q0 r^3 → 単分散半径 r
    struct S { const char* name; double x, y; bool expectChange; double expectRelCorr; };
    const S cases[] = {
        {"monodisperse (1,1)", 1.0, 1.0, false, 0.0},
        {"polydisperse interior (0.8,0.7)", 0.8, 0.7, false, 0.0},
        {"violation y=x^2(1-1e-3)", 0.8, 0.64*(1.0-1.0e-3), true, 1.0e-3},
        {"violation y=x^2(1-1e-7) (within tol)", 0.8, 0.64*(1.0-1.0e-7), false, 0.0},
        {"violation y^2=x(1+1e-3)", 0.5, std::sqrt(0.5*(1.0+1.0e-3)), true, 5.0e-4},
        {"x>1 (1.01)", 1.01, 1.0, true, 1.0e-2},
        {"near-monodisperse (1-1e-4, 1+1e-4)", 1.0-1.0e-4, 1.0+1.0e-4, true, 3.0e-4},
        {"singular boundary (1,1,1,8)-type x=0.5,y=x^2", 0.5, 0.25, false, 0.0},
        {"degenerate x=0 (Q1=0, Q2>0, Q3>0)", 0.0, 0.3, true, 1.0e30},
        {"degenerate (0,0) (Q1=Q2=0, Q3>0)", 0.0, 0.0, true, 1.0e30},
        {"small positive x=1e-3, y=x^2(1-2e-6) (tiny violation -> nearest point, continuous)", 1.0e-3, 1.0e-6*(1.0-2.0e-6), true, 1.0e-5},
        {"small positive interior x=5e-4, y=1e-4", 5.0e-4, 1.0e-4, false, 0.0},
        {"tiny x=1e-14, y=x^2(1-2e-6) (nearest point must stay tiny)", 1.0e-14, 1.0e-28*(1.0-2.0e-6), true, 1.0e-5},
    };
    for (const S& c : cases) {
        std::vector<flow_float> h_ro{(flow_float)ro}, h_g{(flow_float)(ro*g)}, h_q0{(flow_float)(ro*Q0)}, h_q1{(flow_float)(ro*Q0*r*c.x)}, h_q2{(flow_float)(ro*Q0*r*r*c.y)};
        std::vector<flow_float> h_T{(flow_float)T}, h_P{101325.0f}, z{0.0f};
        flow_float *dro=up(h_ro),*dg=up(h_g),*d0=up(h_q0),*d1=up(h_q1),*d2=up(h_q2),*dT=up(h_T),*dP=up(h_P),*cG=up(z),*cQ=up(z);
        int hv[2] = {0, 0}; int* dv = nullptr; cudaMalloc((void**)&dv, 2*sizeof(int)); cudaMemcpy(dv, hv, 2*sizeof(int), cudaMemcpyHostToDevice);
        cond_realizability_clamp_d<<<1,1>>>(1, dro, nullptr, dg, d0, d1, d2, 0, COND_MODEL_H2O, cp.R, 1.0e-9, 5.0e-7, dT, dP, o, cG, cQ, dv, nullptr, nullptr, nullptr, 1);
        cudaDeviceSynchronize();
        flow_float q1n, q2n; cudaMemcpy(&q1n, d1, sizeof(flow_float), cudaMemcpyDeviceToHost); cudaMemcpy(&q2n, d2, sizeof(flow_float), cudaMemcpyDeviceToHost);
        cudaMemcpy(hv, dv, 2*sizeof(int), cudaMemcpyDeviceToHost); const int changed = hv[0] + hv[1];
        const double xn = q1n/(ro*Q0*r), yn = q2n/(ro*Q0*r*r);
        const double rel = std::max(std::fabs(xn - c.x)/std::max(c.x, 1e-30), std::fabs(yn - c.y)/std::max(c.y, 1e-30));
        printf("   %-45s (x,y) %.7f,%.7f -> %.7f,%.7f  corrected=%d (degenerate %d) rel change %.2e\n", c.name, c.x, c.y, xn, yn, hv[0], hv[1], rel);
        CHECK((changed != 0) == c.expectChange, "(l) %s: corrected flag %d (expected %d)", c.name, changed, (int)c.expectChange);
        if (!c.expectChange) CHECK(q1n == h_q1[0] && q2n == h_q2[0], "(l) %s: state changed although admissible", c.name);
        else if (c.expectRelCorr < 1.0) CHECK(rel <= 3.0*c.expectRelCorr, "(l) %s: correction %.2e larger than the violation scale %.2e", c.name, rel, c.expectRelCorr);
        if (c.x == 0.0) CHECK(hv[1] == 1 && xn == 1.0 && yn == 1.0, "(l) %s: degenerate state must reinit to monodisperse", c.name);
        if (c.x > 0.0 && c.expectChange && c.expectRelCorr < 1.0) CHECK(hv[1] == 0, "(l) %s: small positive state must not be reinitialised", c.name);
        CHECK(xn <= 1.0 + 1e-6 && yn >= xn*xn*(1 - 1e-6) && yn*yn <= xn*(1 + 1e-6), "(l) %s: result outside admissible region", c.name);
    }
    // 極端な (x, y): 候補距離が overflow しても必ず許容領域内へ戻ること (codex result-3 の追跡で見つけた欠陥;
    // case/44 run_0415 の確定場に x=1.7e155 の塵セルが残っていた)
    {
        const double ex[][2] = {{1.68e155, 1.68e116}, {5.2e152, 4.2e86}, {1.0e100, 1.0e-5}, {1.0e200, 1.0e200}};
        for (const auto& e : ex) {
            double x = e[0], y = e[1];
            const int rc = cond_realizability_project(x, y, 1.0e-6);
            printf("   extreme (%.2e,%.2e) -> (%.6f,%.6f) code %d\n", e[0], e[1], x, y, rc);
            CHECK(rc != 0, "(l) extreme (%g,%g): must be corrected", e[0], e[1]);
            CHECK(x >= 0.0 && x <= 1.0 + 1e-12 && y >= x*x*(1.0 - 1e-12) && y*y <= x*(1.0 + 1e-12),
                  "(l) extreme (%g,%g) -> (%g,%g) outside admissible region", e[0], e[1], x, y);
        }
    }
    // 半径が表現できない塵 (Q3 underflow) でも x,y を指数分離で作れること (codex result-4 M1)
    {
        const double cases[][4] = {   // q3, q0, q1, q2
            {1.0e-320, 1.0, 1.0e-100, 1.0e-200},
            {1.0e-300, 1.0e-10, 1.0e-120, 1.0e-240},
            {1.0e-200, 1.0, 2.0e-70, 1.0e-140},
        };
        for (const auto& c : cases) {
            double x = -1.0, y = -1.0;
            const bool okxy = cond_moment_xy(c[0], c[1], c[2], c[3], x, y);
            printf("   underflowing radius: q3 %.1e q0 %.1e -> x %.6e y %.6e (ok %d)\n", c[0], c[1], x, y, (int)okxy);
            CHECK(okxy && x == x && y == y && x >= 0.0 && y >= 0.0, "(l) q3=%g q0=%g: x,y must be formed without underflow", c[0], c[1]);
            const int rc = cond_realizability_project(x, y, 1.0e-6);
            CHECK(x >= 0.0 && x <= 1.0 + 1e-12 && y >= x*x*(1.0 - 1e-12) && y*y <= x*(1.0 + 1e-12),
                  "(l) q3=%g: projected state outside admissible region (code %d, x %g y %g)", c[0], rc, x, y);
        }
    }
    // 射影 → 保存量への書き戻し → 再判定 の往復 (codex result-5 M2): q3/q0 が非正規化数に丸められるケースでも
    // 書き戻した保存量から作り直した (x, y) が許容領域に入ること
    {
        const double cases[][4] = {   // q3, q0, q1, q2
            {1.0e-310, 3.0e13, 0.0, 0.0},
            {1.0e-320, 1.0, 1.0e-100, 1.0e-200},
            {1.0e-200, 1.0e5, 2.0e-70, 1.0e-140},
            {1.0e-30,  1.0e14, 1.0e-32, 1.0e-50},
        };
        for (const auto& c : cases) {
            double x = 0.0, y = 0.0;
            const bool okxy = cond_moment_xy(c[0], c[1], c[2], c[3], x, y);
            double q1n = 0.0, q2n = 0.0;
            if (okxy) { cond_realizability_project(x, y, 1.0e-6); cond_moment_writeback(c[0], c[1], x, y, q1n, q2n); }
            else { q1n = 0.0; q2n = 0.0; }
            double x2 = 0.0, y2 = 0.0;
            const bool ok2 = cond_moment_xy(c[0], c[1], q1n, q2n, x2, y2);
            printf("   round trip: q3 %.1e q0 %.1e -> (x,y) %.6f,%.6f -> Q1 %.3e Q2 %.3e -> (x,y) %.6f,%.6f\n", c[0], c[1], x, y, q1n, q2n, x2, y2);
            if (ok2 && (q1n > 0.0 || q2n > 0.0))
                CHECK(x2 <= 1.0 + 1e-6 && y2 >= x2*x2*(1.0 - 1e-6) && y2*y2 <= x2*(1.0 + 1e-6),
                      "(l) round trip q3=%g q0=%g: rebuilt (x,y)=(%g,%g) outside admissible region", c[0], c[1], x2, y2);
        }
    }
    // Q0=0, g>0 (核生成域の g アンダーフロー相当): 触らない
    { std::vector<flow_float> h_ro{1.0f}, h_g{1.0e-12f}, h_q0{0.0f}, h_q1{0.0f}, h_q2{0.0f}, h_T{250.0f}, h_P{101325.0f}, z{0.0f};
      flow_float *dro=up(h_ro),*dg=up(h_g),*d0=up(h_q0),*d1=up(h_q1),*d2=up(h_q2),*dT=up(h_T),*dP=up(h_P),*cG=up(z),*cQ=up(z);
      int hv[2] = {0, 0}; int* dv = nullptr; cudaMalloc((void**)&dv, 2*sizeof(int)); cudaMemcpy(dv, hv, 2*sizeof(int), cudaMemcpyHostToDevice);
      cond_realizability_clamp_d<<<1,1>>>(1, dro, nullptr, dg, d0, d1, d2, 0, COND_MODEL_H2O, cp.R, 1.0e-9, 5.0e-7, dT, dP, o, cG, cQ, dv, nullptr, nullptr, nullptr, 1);
      cudaDeviceSynchronize(); cudaMemcpy(hv, dv, 2*sizeof(int), cudaMemcpyDeviceToHost);
      flow_float gn; cudaMemcpy(&gn, dg, sizeof(flow_float), cudaMemcpyDeviceToHost);
      CHECK(hv[0] == 0 && hv[1] == 0 && gn == 1.0e-12f, "(l) Q0=0,g>0 dust must not be projected (count %d/%d, g %g)", hv[0], hv[1], gn); printf("   Q0=0,g>0: untouched (count %d)\n", hv[0]+hv[1]); }
}

int main()
{
    test_realizability_projection();
    // (a) Δτ 不変性: H2O carrier (T × S × Q0 × g) と N2 pure
    { std::vector<State> st; CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
      const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o); const double Yw = 0.0377, Rw = cp.R, Rmix = 285.0;
      for (double T = 220.0; T <= 290.0; T += 10.0) for (double S : {0.5, 0.99, 1.01, 2.0, 30.0})
        for (double q0 : {0.0, 1.0e13, 1.0e16}) for (double g : {0.0, 1.0e-5, 1.0e-2, 0.03}) for (double rb : {1.0e-8, 1.0e-7}) {
            if (g >= Yw) continue;
            const double pv = S*cond_psat(cp, T); const double ro = pv/((Yw - g)*Rw*T);
            st.push_back({T, ro*Rmix*T, ro, Yw, g, q0, q0*rb, q0*rb*rb});
        }
      // 蒸気枯渇 (g = Yw − 1e-9): 残差側で S=0 になるべき
      for (double T : {230.0, 250.0}) { const double g = Yw - 1.0e-9; const double pv = 30.0*cond_psat(cp, T); const double ro = pv/((Yw - g)*Rw*T);
          st.push_back({T, ro*Rmix*T, ro, Yw, g, 1.0e16, 1.0e9, 1.0e2}); }
      printf("== (a) H2O TP carrier: %zu states ==\n", st.size());
      test_dt_invariance("H2O", COND_MODEL_H2O, 1, st); }
    { std::vector<State> st; CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
      const CondSpeciesProps cp = condProps_make(COND_MODEL_N2, o);
      for (double T = 40.0; T <= 100.0; T += 10.0) for (double S : {0.5, 0.99, 1.1, 10.0})
        for (double q0 : {0.0, 1.0e15}) for (double g : {0.0, 1.0e-3, 5.0e-2}) for (double rb : {3.0e-8, 5.0e-7}) {
            const double P = S*cond_psat(cp, T); if (P > 1.0e6) continue; const double ro = P/((1.0 - g)*cp.R*T);
            st.push_back({T, P, ro, 1.0, g, q0, q0*rb, q0*rb*rb}); }
      printf("== (a) N2 pure CPG: %zu states ==\n", st.size());
      test_dt_invariance("N2", COND_MODEL_N2, 0, st); }
    printf("== (b)-(f) update limiter kernel ==\n");
    test_update_limiter();
    printf("== (g) evaporation source values (polydisperse) ==\n");
    test_evap_source_values();
    printf("== (h) one-cell transport+source fixed point vs dt ==\n");
    test_one_cell_fixed_point();
    printf("== (k) small droplet evaporation through source -> update clamp -> realizability clamp ==\n");
    test_small_droplet_evaporates();
    printf("%s (%d failures)\n", g_fail ? "FAILED" : "ALL PASS", g_fail);
    return g_fail ? 1 : 0;
}
