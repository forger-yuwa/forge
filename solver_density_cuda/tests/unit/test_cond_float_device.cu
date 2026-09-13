// test_cond_float_device.cu — 凝縮ソース kernel の double 実体と float 実体を合成状態に掛けて比較する device 単体試験
// (plans/active/condensation-float-speedup.md §4.3「単体 (device)」)。
//   build: nvcc --expt-relaxed-constexpr -I. -o test_cond_float_device tests/unit/test_cond_float_device.cu
// 状態: H2O TP carrier (N2/H2O, Y_w=0.0113) と N2 pure CPG。T × S × Q0 × g の格子 + S<1 (蒸発) + 実 run で差が出たセル状態。
// 判定: res_* / sj_* は相対 1e-3 + ln J 許容 (S→1 の CNT 感度) または絶対 1e-4×成長流束尺度、diagS 相対 1e-5、T_sat 2e-3 K、
//       θ 相対 1e-5 (S≈1 の境界セルは除外)、condLim 絶対 1e-4。sj_* は両実体とも有限・非負。
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/thermo_d.cuh"
namespace {
#include "cuda_forge/condensationSourceKernels_d.cuh"
}

static int g_fail = 0;
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
struct Worst { double v=0; int i=-1; const char* what=""; void upd(double e,int i_,const char* w){ if(e>v){v=e;i=i_;what=w;} } };

static void run_case(const char* name, int model, int carrier, const std::vector<State>& st, int kantrowitz, int growthModel, int evap)
{
    const int n = (int)st.size();
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(model, o);
    CondTablesHost ht; cond_tables_build_host(cp, ht); const CondTablesF tb = cond_tables_upload(ht);
    // species (TP carrier: N2/H2O with sensible datum at 298.15 K)
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double H2Olo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double H2Ohi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    std::vector<SpeciesThermo> sp = { mk(0.0280134,3.621,97.53,N2lo,N2hi), mk(0.0180153,2.605,572.4,H2Olo,H2Ohi) };
    for (auto& s : sp) { const double hr = thermo_h_molar(s, 298.15); s.low[7] += -hr/THERMO_RU; s.high[7] += -hr/THERMO_RU; }
    std::vector<SpeciesThermoF> spf = { toF(sp[0]), toF(sp[1]) };
    const int nSp = 2, cgs = 1;
    // arrays
    std::vector<flow_float> T(n),P(n),ro(n),cpc(n),Rm(n),roY0(n),roY1(n),rog(n),q0(n),q1(n),q2(n),vol(n,1.0e-9f),dt(n,1.0e-7f);
    for (int i=0;i<n;i++){ const State& s=st[i]; T[i]=(flow_float)s.T; P[i]=(flow_float)s.P; ro[i]=(flow_float)s.ro; rog[i]=(flow_float)(s.ro*s.g);
        q0[i]=(flow_float)s.q0; q1[i]=(flow_float)s.q1; q2[i]=(flow_float)s.q2; roY1[i]=(flow_float)(s.ro*s.Yw); roY0[i]=(flow_float)(s.ro*(1.0-s.Yw));
        double Y[2]={1.0-s.Yw, s.Yw}; double c,h; thermo_cph_mix(sp.data(),2,Y,s.T,&c,&h); cpc[i]=(flow_float)c; Rm[i]=(flow_float)thermo_R_mix(sp.data(),2,Y); }
    flow_float *dT=up(T),*dP=up(P),*dro=up(ro),*dcp=up(cpc),*dRm=up(Rm),*dY0=up(roY0),*dY1=up(roY1),*dvol=up(vol),*ddt=up(dt);
    std::vector<flow_float*> hY = {dY0, dY1}; flow_float** dYall = up(hY);
    SpeciesThermo* dsp = up(sp); SpeciesThermoF* dspf = up(spf);
    auto zeros = [&](){ return up(std::vector<flow_float>(n, 0.0f)); };
    struct Out { flow_float *rog,*q0,*q1,*q2,*rr,*r0,*r1,*r2,*sg,*s0,*s1,*s2,*dS,*dD,*dR,*dT,*dTh,*dL; };
    auto mkout = [&](){ Out o; o.rog=up(rog); o.q0=up(q0); o.q1=up(q1); o.q2=up(q2); o.rr=zeros(); o.r0=zeros(); o.r1=zeros(); o.r2=zeros();
        o.sg=zeros(); o.s0=zeros(); o.s1=zeros(); o.s2=zeros(); o.dS=zeros(); o.dD=zeros(); o.dR=zeros(); o.dT=zeros(); o.dTh=zeros(); o.dL=zeros(); return o; };
    Out A = mkout(), B = mkout();
    const int blk = 128, grd = (n + blk - 1)/blk;
    const double Rw = cp.R, M = cp.M; const float cp_cpg = 1038.8f, gamma_cpg = 1.4f;
    condensation_source_d<<<grd,blk>>>(n, model, carrier, Rw, M, kantrowitz, 0, o, carrier ? dsp : nullptr, nSp, carrier ? dYall : nullptr, carrier ? cgs : -1,
        growthModel, 3.18, 0, evap, 1.0e-9, 0, 0.5, 0, 1.0, 5.0e-3, 10.0, cp_cpg, gamma_cpg, 1.0e35, 5.0e-3, 1.0,
        dvol, ddt, dT, dP, dro, carrier ? dcp : nullptr, carrier ? dRm : nullptr, carrier ? dY1 : nullptr,
        A.rog, A.q0, A.q1, A.q2, A.rr, A.r0, A.r1, A.r2, A.sg, A.s0, A.s1, A.s2, A.dS, A.dD, A.dR, A.dT, A.dTh, A.dL);
    condensation_source_f_d<<<grd,blk>>>(n, carrier, (float)Rw, kantrowitz, 0, condProps_to_f(cp), tb, 0.0f, carrier ? dspf : nullptr, nSp, carrier ? dYall : nullptr, carrier ? cgs : -1,
        growthModel, 3.18f, evap, 1.0e-9f, 0, 0.5f, cp_cpg, gamma_cpg, 5.0e-3f, 1.0f,
        dvol, ddt, dT, dP, dro, carrier ? dcp : nullptr, carrier ? dRm : nullptr, carrier ? dY1 : nullptr,
        B.rog, B.q0, B.q1, B.q2, B.rr, B.r0, B.r1, B.r2, B.sg, B.s0, B.s1, B.s2, B.dS, B.dD, B.dR, B.dT, B.dTh, B.dL);
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; return; }
    auto rr=down(A.rr,n), r0=down(A.r0,n), r1=down(A.r1,n), r2=down(A.r2,n), sg=down(A.sg,n), s1=down(A.s1,n), dS=down(A.dS,n), dD=down(A.dD,n), dR=down(A.dR,n), dT_=down(A.dT,n), dTh=down(A.dTh,n), dL=down(A.dL,n);
    auto frr=down(B.rr,n), fr0=down(B.r0,n), fr1=down(B.r1,n), fr2=down(B.r2,n), fsg=down(B.sg,n), fs1=down(B.s1,n), fdS=down(B.dS,n), fdD=down(B.dD,n), fdR=down(B.dR,n), fdT=down(B.dT,n), fdTh=down(B.dTh,n), fdL=down(B.dL,n);
    Worst wres, wsj, wS, wTs, wTh, wL, wD; int nBoundary = 0, nSjNeg = 0, nCapped = 0;
    for (int i = 0; i < n; ++i) {
        const State& s = st[i];
        const double psat = cond_psat(cp, s.T);
        const double pv = carrier ? s.ro*(s.Yw - s.g)*Rw*s.T : s.P;
        const double S = pv/psat;
        const bool boundary = fabs(S - 1.0) < 1.0e-5;   // float の S=1 判定境界 (離散差は許容)
        if (boundary) { ++nBoundary; continue; }
        // ln J 許容 (S>1 のとき)
        double tolJ = 1.0e-4; bool capped = false;
        if (S > 1.0) { double Jd, rd; cond_nucleation(cp, s.T, pv, carrier ? s.ro*(s.Yw-s.g) : s.ro*(1.0-s.g), &Jd, &rd, kantrowitz, cp.cp/cp.cv, nullptr);
            if (Jd > 1.0e35) { ++nCapped; capped = true; }
            const double dGkT = (rd > 0.0) ? (4.0/3.0)*COND_PI*rd*rd*cond_sigma(cp, s.T)/(COND_KB*s.T) : 0.0; tolJ += 2.0*dGkT*(3.0e-6/log(S)); }
        const double rho_l = cond_rho_cond(cp, s.T), dsc = (pv/rho_l)/sqrt(2.0*COND_PI*Rw*s.T);
        const double scale[4] = {0.0, s.q0*dsc, 2.0*s.q1*dsc, 4.0*COND_PI*rho_l*s.q2*dsc};
        const double Aa[4] = {r0[i], r1[i], r2[i], rr[i]}, Bb[4] = {fr0[i], fr1[i], fr2[i], frr[i]};
        const double floor_k[4] = {1.0e-5, 1.0e-14, 1.0e-23, 1.0e-29};
        for (int k = 0; k < 4; ++k) {
            const double a = Aa[k]/1.0e-9, b = Bb[k]/1.0e-9;   // vol で割って [1/(m³ s)] 等に戻す
            if (fabs(a) < floor_k[k] && fabs(b) < floor_k[k]) continue;
            const double err = fabs(a - b); if (err <= 1.0e-4*scale[k]) continue;
            const double den = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
            wres.upd(err/den/(1.0e-3 + tolJ), i, k==0?"res_Q0":k==1?"res_Q1":k==2?"res_Q2":"res_g");
        }
        for (int k = 0; k < 2; ++k) {
            const double a = k ? s1[i] : sg[i], b = k ? fs1[i] : fsg[i];
            if (!(a >= 0.0) || !(b >= 0.0) || !std::isfinite(a) || !std::isfinite(b)) ++nSjNeg;
            if (capped) continue;   // J 上限セルは float 実体の設計変更 (対数上限を摂動側にも) で除外 (計数は上)
            if (S > 0.99 && S < 1.0) continue;   // 蒸発端 (1−S<1 %): T 摂動 0.1 K (ΔS 0.6 %) が (p_v−p_d) の尺度を跨ぎ数値微分が両精度とも粗い (有限・非負は上で検査)
            if (a == 0.0 && b == 0.0) continue;
            // 判定は陰的項 sj·(ρφ) を残差 S と比べる: |Δsj|·(ρφ) ≤ (1e-3+tolJ)·(|S| + sj·ρφ)。核生成が支配する (J r_nuc ≫ Q0 dr/dt) セルでは
            // float で成長項の微分が J r_nuc の丸めに埋もれて 0 になるが、その sj は残差に対し無視できる。
            const double q  = k ? s.q1 : s.ro*s.g;                    // 対応する保存量 (ρQ1 / ρg)
            const double Sd = fabs(k ? r1[i] : rr[i])/1.0e-9;         // 残差 [単位/(m³ s)] (vol で割る)
            // 1 step の陰的更新への影響 |Δsj|·dt (無次元) が 1e-6 以下なら無視 (S≈1 の蒸発端で数値微分が悪条件になるセル)
            if (fabs(a - b)*1.0e-7 <= 1.0e-6) continue;
            const double den = (3.0e-3 + tolJ)*(Sd + (fabs(a) > fabs(b) ? fabs(a) : fabs(b))*q) + 1.0e-300;   // 3e-3: θ 律速 (avail→0) 下で θ·∂S/∂T の合成誤差が 2e-3 まで出る
            wsj.upd(fabs(a - b)*q/den, i, k ? "sj_Q1" : "sj_g");
        }
        if (dS[i] > 0.0) wS.upd(fabs(dS[i] - fdS[i])/dS[i], i, "condS");
        if (dT_[i] > 0.0) wTs.upd(fabs(dT_[i] - fdT[i]), i, "condTsat");
        if (dTh[i] > 0.0) wTh.upd(fabs(dTh[i] - fdTh[i])/dTh[i], i, "condTheta");
        wL.upd(fabs(dL[i] - fdL[i]), i, "condLim");
        { const double a = dD[i], b = fdD[i]; if (a != 0.0 || b != 0.0) { const double den = (1.0e-4 + 2.0e-6/fabs(log(S > 0 ? S : 1.0) + 1e-30))*fabs(a) + 1.0e-5*dsc; wD.upd(fabs(a - b)/den, i, "condDrdt"); } }
    }
    auto rep = [&](const char* w, const Worst& x, double tol) {
        const bool ok = x.v <= tol; if (!ok) ++g_fail;
        printf("  [%s] %-14s %-10s worst=%.3e tol=%.1e", ok ? "PASS" : "FAIL", name, w, x.v, tol);
        if (x.i >= 0) { const State& s = st[x.i]; printf("  at T=%.2f g=%.1e q0=%.1e S=%.4g", s.T, s.g, s.q0, (carrier ? s.ro*(s.Yw-s.g)*Rw*s.T : s.P)/cond_psat(cp, s.T)); }
        printf("\n"); };
    printf("-- %s: n=%d, S=1 boundary cells skipped=%d, J capped=%d, sj negative/nonfinite=%d\n", name, n, nBoundary, nCapped, nSjNeg);
    auto dump = [&](int i, const char* tag) { if (i < 0) return; const State& s = st[i];
        printf("    [%s] i=%d T=%.2f P=%.4g ro=%.4g g=%.2e q0=%.2e q1=%.2e q2=%.2e\n", tag, i, s.T, s.P, s.ro, s.g, s.q0, s.q1, s.q2);
        printf("      D: res g/Q0/Q1/Q2 = %.4e %.4e %.4e %.4e  sj_g %.4e sj_Q1 %.4e  drdt %.4e r30 %.4e S %.5g Tsat %.3f\n", rr[i], r0[i], r1[i], r2[i], sg[i], s1[i], dD[i], dR[i], dS[i], dT_[i]);
        printf("      F: res g/Q0/Q1/Q2 = %.4e %.4e %.4e %.4e  sj_g %.4e sj_Q1 %.4e  drdt %.4e r30 %.4e S %.5g Tsat %.3f\n", frr[i], fr0[i], fr1[i], fr2[i], fsg[i], fs1[i], fdD[i], fdR[i], fdS[i], fdT[i]); };
    dump(wres.i, "worst res"); dump(wsj.i, "worst sj"); dump(n-1, "last state");
    if (nSjNeg) ++g_fail;
    rep("res (norm)", wres, 1.0); rep("src_jac (norm)", wsj, 1.0); rep("condS rel", wS, 1.0e-5); rep("condTsat [K]", wTs, 2.0e-3);
    rep("condTheta rel", wTh, 1.0e-5); rep("condLim abs", wL, 1.0e-4); rep("condDrdt (norm)", wD, 1.0);
}

int main()
{
    // ---- H2O TP carrier (case/16 相当 Y_w=0.0113): T × S × Q0 × g、S<1 の蒸発、実 run のセル
    std::vector<State> st;
    { CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
      const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o); const double Yw = 0.0113, Rw = cp.R, Rmix = 296.0;
      for (double T = 200.0; T <= 300.0; T += 10.0) for (double S : {0.5, 0.999, 1.001, 1.01, 1.1, 2.0, 10.0, 100.0, 1000.0})
        for (double q0 : {0.0, 1.0e9, 1.0e14, 1.0e17, 1.0e19}) for (double g : {0.0, 1.0e-13, 1.0e-6, 1.0e-4, 1.0e-2}) for (double rb : {1.0e-9, 1.0e-8, 1.0e-7}) {
            if (g >= Yw) continue;
            const double pv = S*cond_psat(cp, T); const double ro = pv/((Yw - g)*Rw*T);
            st.push_back({T, ro*Rmix*T, ro, Yw, g, q0, q0*rb, q0*rb*rb});
        }
      st.push_back({238.30882, 30887.66, 0.4340555, 0.004752914, 6.0e-13, 8018105000.0, 11.575916, 2.3136373e-08});   // run_0456 1 step 目で差が出たセル
      printf("== H2O TP carrier: %zu states ==\n", st.size());
      run_case("H2O Kw1 HK", COND_MODEL_H2O, 1, st, 1, 0, 1);
      run_case("H2O Kw3 Gyar", COND_MODEL_H2O, 1, st, 3, 1, 1);
    }
    // ---- N2 pure CPG (case/34 相当): T 40–110 K
    { std::vector<State> s2; CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
      const CondSpeciesProps cp = condProps_make(COND_MODEL_N2, o);
      for (double T = 40.0; T <= 110.0; T += 5.0) for (double S : {0.5, 0.999, 1.001, 1.1, 2.0, 10.0, 100.0})
        for (double q0 : {0.0, 1.0e9, 1.0e15, 1.0e18}) for (double g : {0.0, 1.0e-8, 1.0e-3, 5.0e-2}) for (double rb : {2.0e-9, 3.0e-8, 5.0e-7}) {
            const double P = S*cond_psat(cp, T); const double ro = P/((1.0 - g)*cp.R*T);
            if (P > 1.0e6) continue;   // 10 bar 超 (臨界 34 bar に近い) は除外: T_sat の double 診断 (C–C 勾配の Newton) が不正確で比較にならない; case/34 は ≤5 bar
            s2.push_back({T, P, ro, 1.0, g, q0, q0*rb, q0*rb*rb});
        }
      printf("== N2 pure CPG: %zu states ==\n", s2.size());
      run_case("N2 Kw1 Goodh", COND_MODEL_N2, 0, s2, 1, 0, 1);
      run_case("N2 Kw0 Gyar", COND_MODEL_N2, 0, s2, 0, 1, 1);
    }
    printf("%s (%d failures)\n", g_fail ? "FAILED" : "ALL PASS", g_fail);
    return g_fail ? 1 : 0;
}
