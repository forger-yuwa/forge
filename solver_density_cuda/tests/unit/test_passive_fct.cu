// test_passive_fct.cu — dual-time 物理 step 末尾の保存的 FCT 補正 (passiveFct_d.cuh, plan §4.7 v5) の単体試験
//   (plans/active/species-passive-scalar-unification.md §4.7 / §6-2; codex result-2 M2, plan-3〜6)。
//   build: nvcc --expt-relaxed-constexpr -arch=sm_86 -I. -o test_passive_fct tests/unit/test_passive_fct.cu
//   (i)  滑らかなガウス (CFL 0.5, BDF2, 前制限なし): 変更 (>1e-6) は極値近傍の ≤6 セルだけ、変更量 ≤5e-3 (FCT の極値クリップ; 裾の 1e-7 級の変更は数えない)
//   (ii) ステップ (CFL 2, BDF1 / BDF2 [流束形履歴 G/H]): 補正前後の総量 (float 1e-6, host double 参照 1e-12), 0≤φ≤1, 局所限界内, device/host 一致 ≤1e-5,
//        低次 BE 解 q_L の有界性
//   (iii) α≡0 → q_C = q_B, q_B − q_L = M^{-1}(r_L − r_H/a) の範囲
//   (iv) plan-4 M1 の 3 CV 反例 (BE, CFL 1): 生流束の引き戻しで負値なし
//   (v)  plan-5 M2 の非周期 3 CV 反例 (node 逆流境界, M=0.5): 低次作用素が M 行列 (流入は φ^n の RHS) で q_L 有界、q_C ≤ 局所上限
//   (vi) 無流束 1 CV の BDF2 履歴 q^n=0.9, q^{n-1}=0.5: 流束形なら G^n が履歴を運び基点は有界; 局所形 (G/H 欠落相当) は逸脱 1/30 が記録される
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/limiterFunctions_d.cuh"
#include "cuda_forge/passiveLimiter_d.cuh"
#include "cuda_forge/passiveKernels_d.cuh"
#include "cuda_forge/passiveFct_d.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)
template <class T> static T* up(const std::vector<T>& v) { T* d = nullptr; cudaMalloc((void**)&d, v.size()*sizeof(T)); cudaMemcpy(d, v.data(), v.size()*sizeof(T), cudaMemcpyHostToDevice); return d; }
template <class T> static std::vector<T> down(const T* d, size_t n) { std::vector<T> v(n); cudaMemcpy(v.data(), d, n*sizeof(T), cudaMemcpyDeviceToHost); return v; }

// 1 次元周期セル鎖: 面 ip=(ip, (ip+1)%n)、全面内部。
struct Ring {
    int n; float h;
    std::vector<geom_int> plane_cells, cpi, cp, nhp;
    std::vector<geom_float> vol, ccx, ccy, ccz, pcx, pcy, pcz, fx;
    Ring(int n_, float h_) : n(n_), h(h_) {
        plane_cells.assign(2*n, 0);
        for (int ip = 0; ip < n; ++ip) { plane_cells[2*ip] = ip; plane_cells[2*ip+1] = (ip+1)%n; }
        cpi.assign(n+1, 0); cp.clear();
        for (int ic = 0; ic < n; ++ic) { cpi[ic] = (geom_int)cp.size(); cp.push_back((ic+n-1)%n); cp.push_back(ic); }
        cpi[n] = (geom_int)cp.size();
        for (int ip = 0; ip < n; ++ip) nhp.push_back(ip);
        vol.assign(n, h); ccx.resize(n); ccy.assign(n, 0.f); ccz.assign(n, 0.f);
        for (int ic = 0; ic < n; ++ic) ccx[ic] = (ic+0.5f)*h;
        pcx.resize(n); pcy.assign(n, 0.f); pcz.assign(n, 0.f); fx.assign(n, 0.5f);
        for (int ip = 0; ip < n; ++ip) pcx[ip] = (ip+1)*h;
    }
};

__global__ void face_recon_d(int nPlanes, int nCells, const geom_int* plane_cells, const geom_float* ccx, const geom_float* pcx,
                             const float* phi, const float* gx, const float* lim, const float* massflux, float* Pface)
{
    const int ip = blockDim.x*blockIdx.x + threadIdx.x;
    if (ip >= nPlanes) return;
    const int ic0 = plane_cells[2*ip], ic1 = plane_cells[2*ip+1];
    float pl = phi[ic0] + lim[ic0]*gx[ic0]*(pcx[ip]-ccx[ic0]);
    float pr = phi[ic1] + lim[ic1]*gx[ic1]*(pcx[ip]-ccx[ic1]);
    if (ip == nPlanes-1) { pl = phi[ic0]; pr = phi[ic1]; }
    pl = fminf(fmaxf(pl, 0.f), 1.f); pr = fminf(fmaxf(pr, 0.f), 1.f);
    Pface[ip] = (massflux[ip] >= 0.f) ? pl : pr;
}
__global__ void dual_update_d(int n, float V, float dt, float dtau, float a, float b, float c, const float* qP, const float* qPP,
                              const float* res, const float* td, float* q, double* resnorm)
{
    const int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= n) return;
    const float r = res[ic] - (V/dt)*(a*q[ic] - b*qP[ic] + c*qPP[ic]);
    q[ic] += r / (V/dtau + td[ic] + V*a/dt);
    atomicAdd(resnorm, (double)r*(double)r);
}

// FCT 1 step を実行する汎用 device ハーネス (任意の面リスト; 周期 gather なし)。入力: q_H (dq), q^n (dqP), G^n/H^n (nullptr 可), Pface, massflux, meff。
struct FctDev {
    int n, nPlanes, nNormal, nHalo; geom_int *dpc,*dcpi,*dcp,*dnhp; geom_float *dvol,*dccx,*dfx; float *dro,*droN;
    float *dqL,*drhs,*dnb,*dPp,*dPm,*dmx,*dmn,*dRp,*dRm,*dcorr,*dbase,*dqB,*dH,*dqP,*dqPP,*dq,*ddA,*ddD,*dcd,*dAf,*dAp,*dG,*dGn,*dmeff;
    float **pq,**pqP,**pqPP,**pqL,**prhs,**pnb,**pPp,**pPm,**pmx,**pmn,**pRp,**pRm,**pcorr,**pbase,**pqB,**pH;
    double *dst,*dch,*drhs2; float* pin;
    FctDev(int n_, int nPlanes_, int nNormal_, const std::vector<geom_int>& pc, const std::vector<geom_int>& cpi, const std::vector<geom_int>& cp, const std::vector<geom_int>& nhp,
           const std::vector<geom_float>& vol, const std::vector<geom_float>& ccx, const std::vector<float>& ro, int nGhost)
        : n(n_), nPlanes(nPlanes_), nNormal(nNormal_), nHalo((int)nhp.size()) {
        dpc=up(pc); dcpi=up(cpi); dcp=up(cp); dnhp=up(nhp); dvol=up(vol); dccx=up(ccx); dfx=up(std::vector<geom_float>(nPlanes_,0.5f));
        dro=up(ro); droN=up(ro);
        std::vector<float> z((size_t)n_+nGhost, 0.f), zp((size_t)nPlanes_, 0.f);
        dqL=up(z);drhs=up(z);dnb=up(z);dPp=up(z);dPm=up(z);dmx=up(z);dmn=up(z);dRp=up(z);dRm=up(z);dcorr=up(z);dbase=up(z);dqB=up(z);dH=up(z);dqP=up(z);dqPP=up(z);dq=up(z);ddA=up(z);ddD=up(z);
        dcd=up(zp);dAf=up(zp);dAp=up(zp);dG=up(zp);dGn=up(zp);dmeff=up(zp);
        auto P=[&](float* a){ std::vector<float*> v={a}; return up(v); };
        pq=P(dq);pqP=P(dqP);pqPP=P(dqPP);pqL=P(dqL);prhs=P(drhs);pnb=P(dnb);pPp=P(dPp);pPm=P(dPm);pmx=P(dmx);pmn=P(dmn);pRp=P(dRp);pRm=P(dRm);pcorr=P(dcorr);pbase=P(dbase);pqB=P(dqB);pH=P(dH);
        dst=up(std::vector<double>(8,0.0)); dch=up(std::vector<double>(1,0.0)); drhs2=up(std::vector<double>(1,0.0)); pin=nullptr;
    }
    // a,c: BDF 係数; useHist: G/H (device 配列 dG/dH) を使う; massflux/meff は dmeff/mf に入れて渡す。戻り: q_C を dq に書く。
    void run(float dt, float a, float c, bool useHist, const float* dmf, const float* dPface, int isNode, int prelimit, int forceAlpha0, double tol, int maxSweep) {
        const float invA = useHist ? 1.f/a : 1.f, cOverA = useHist ? c/a : 0.f;
        std::vector<float> mfh = down(dmf, (size_t)nPlanes), meffh = down(dmeff, (size_t)nPlanes);
        for (int ip = 0; ip < nPlanes; ++ip) meffh[ip] = invA*mfh[ip] + cOverA*meffh[ip];
        cudaMemcpy(dmeff, meffh.data(), nPlanes*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemset(ddA, 0, n*sizeof(float)); cudaMemset(ddD, 0, n*sizeof(float)); cudaMemset(dcd, 0, nPlanes*sizeof(float)); cudaMemset(drhs2, 0, sizeof(double));
        passive_fct_lo_diag_d<<<(nHalo+127)/128,128>>>(n, nNormal, nHalo, dnhp, dpc, dro, dmeff, isNode, 0, dccx, dccx, dccx, dfx, dccx, dccx, dccx, dccx, dro, dro, 1.f, 1.f, ddA, ddD, dcd);
        passive_fct_lo_rhs_d<<<(n+127)/128,128>>>(n, dvol, dt, invA, cOverA, 1, pqP, pqPP, nullptr, useHist ? pH : nullptr, nullptr, prhs, drhs2);
        passive_fct_lo_init_d<<<(n+127)/128,128>>>(n, dro, 1, 1, pq, pqL);
        double rhs2; cudaMemcpy(&rhs2, drhs2, sizeof(double), cudaMemcpyDeviceToHost);
        for (int sw = 0; sw < maxSweep; ++sw) {
            cudaMemset(dnb, 0, n*sizeof(float)); cudaMemset(dch, 0, sizeof(double));
            passive_fct_lo_nb_d<<<(nHalo+127)/128,128>>>(n, nNormal, nHalo, dnhp, dpc, dro, droN, dmeff, isNode, 1, -1, dcd, pqL, pq, pqP, pnb);
            passive_fct_lo_solve_d<<<(n+127)/128,128>>>(n, dvol, 1.f/dt, 1, -1, ddA, ddD, prhs, pnb, pq, pin, nullptr, pqL, dch);
            double r2; cudaMemcpy(&r2, dch, sizeof(double), cudaMemcpyDeviceToHost);
            if (std::sqrt(r2) <= tol*(std::sqrt(rhs2) + 1e-30)) break;
        }
        cudaMemset(dbase, 0, n*sizeof(float)); cudaMemset(dPp, 0, n*sizeof(float)); cudaMemset(dPm, 0, n*sizeof(float)); cudaMemset(dcorr, 0, n*sizeof(float)); cudaMemset(dst, 0, 8*sizeof(double));
        passive_fct_raw_d<<<(nHalo+127)/128,128>>>(n, nNormal, nHalo, dnhp, dpc, dro, droN, dmf, dmeff, isNode, invA, cOverA, 1, 1, dPface, -1, nullptr, useHist ? dG : nullptr, pq, pqL, pqP, dAf, pbase);
        passive_fct_base_d<<<(n+127)/128,128>>>(n, dvol, dt, dro, 1, 1, pq, pbase, nullptr, pqB, dst);
        passive_fct_prelimit_d<<<(nHalo+127)/128,128>>>(n, nHalo, dnhp, dpc, dro, 1, 1, dAf, prelimit, pqB, dAp, pPp, pPm, dst);
        passive_fct_extrema_d<<<(n+127)/128,128>>>(n, nNormal, dpc, dcpi, dcp, dro, droN, 1, pqP, pqB, pmx, pmn);
        passive_fct_ratio_d<<<(n+127)/128,128>>>(n, dvol, 1.f/dt, dro, 1, 1, pqB, pmx, pmn, pPp, pPm, pin, pRp, pRm);
        passive_fct_apply_d<<<(nHalo+127)/128,128>>>(n, nNormal, nHalo, dnhp, dpc, dro, droN, dmeff, isNode, 1, 1, dAf, dAp, pRp, pRm, forceAlpha0, -1, nullptr, pqL, pqP, pcorr, dGn, dst);
        passive_fct_commit_d<<<(n+127)/128,128>>>(n, dvol, dt, 1, pcorr, pin, nullptr, pq, dst);
        cudaDeviceSynchronize();
    }
    // 確定状態から G^{n+1}, H^{n+1}
    void finish(float dt) {
        std::swap(dG, dGn);
        cudaMemset(dbase, 0, n*sizeof(float));
        passive_fct_divG_d<<<(nHalo+127)/128,128>>>(n, nHalo, dnhp, dpc, 1, 1, dG, pbase);
        passive_fct_hist_local_d<<<(n+127)/128,128>>>(n, dvol, dt, 1, pq, pqP, pbase, pH);
        cudaDeviceSynchronize();
    }
};

struct Case { std::vector<float> qH, qNew, qL, qP, qAlpha0, qB; double massH, massNew, massHost; float minL, maxL; std::vector<double> host; int nModified; float maxMod; };

// host double 参照 (周期鎖, 流束形履歴 G/H, 前制限なし)。
static std::vector<double> host_fct(int n, double V, double dt, double a, double c, bool useHist, const std::vector<float>& qH, const std::vector<float>& qP,
                                    const std::vector<float>& Pface, const std::vector<float>& G, const std::vector<float>& H, double mdot, double* massOut, double* minL, double* maxL)
{
    const double invA = useHist ? 1.0/a : 1.0, cA = useHist ? c/a : 0.0;
    const double meff = invA*mdot + cA*mdot;   // 一様流: ṁ^eff = ṁ
    std::vector<double> qL(n), rhs(n), nb(n), diag(n, V/dt + meff);
    for (int i = 0; i < n; ++i) { qL[i] = std::min(1.0, std::max(0.0, (double)qH[i])); rhs[i] = (V/dt)*qP[i] + (useHist ? cA*H[i] : 0.0); }
    for (int it = 0; it < 200000; ++it) {
        std::fill(nb.begin(), nb.end(), 0.0);
        for (int ip = 0; ip < n; ++ip) nb[(ip+1)%n] += meff*qL[ip];
        double r2 = 0.0, f2 = 0.0;
        for (int i = 0; i < n; ++i) { const double rl = rhs[i]+nb[i]-diag[i]*qL[i]; r2 += rl*rl; f2 += rhs[i]*rhs[i]; qL[i] = (rhs[i]+nb[i])/diag[i]; }
        if (std::sqrt(r2) <= 1e-15*std::sqrt(f2)) break;
    }
    *minL = *std::min_element(qL.begin(), qL.end()); *maxL = *std::max_element(qL.begin(), qL.end());
    std::vector<double> A(n), Pp(n, 0.0), Pm(n, 0.0), mx(n), mn(n), Rp(n), Rm(n), corr(n, 0.0), base(n, 0.0), qB(n);
    for (int ip = 0; ip < n; ++ip) {
        const int i0 = ip, i1 = (ip+1)%n;
        const double FH = invA*mdot*(double)Pface[ip] + (useHist ? cA*G[ip] : 0.0), FL = meff*qL[i0];
        A[ip] = FH - FL; base[i0] -= A[ip]; base[i1] += A[ip];
    }
    for (int i = 0; i < n; ++i) qB[i] = std::min(1.0, std::max(0.0, (double)qH[i] - (dt/V)*base[i]));
    for (int ip = 0; ip < n; ++ip) { const int i0 = ip, i1 = (ip+1)%n; const double Af = A[ip];
        if (Af > 0) { Pm[i0] += Af; Pp[i1] += Af; } else if (Af < 0) { Pp[i0] -= Af; Pm[i1] -= Af; } }
    for (int i = 0; i < n; ++i) {
        const int l = (i+n-1)%n, r = (i+1)%n;
        mx[i] = std::min(1.0, std::max({(double)qP[i], qB[i], (double)qP[l], qB[l], (double)qP[r], qB[r]}));
        mn[i] = std::max(0.0, std::min({(double)qP[i], qB[i], (double)qP[l], qB[l], (double)qP[r], qB[r]}));
        const double Qp = (mx[i]-qB[i])*V/dt, Qm = (qB[i]-mn[i])*V/dt;
        Rp[i] = (Pp[i] > 0) ? ((Qp <= 0) ? 0.0 : std::min(1.0, Qp/Pp[i])) : 1.0;
        Rm[i] = (Pm[i] > 0) ? ((Qm <= 0) ? 0.0 : std::min(1.0, Qm/Pm[i])) : 1.0;
    }
    for (int ip = 0; ip < n; ++ip) { const int i0 = ip, i1 = (ip+1)%n; const double Ar = A[ip]; if (Ar == 0.0) continue;
        const double al = (Ar > 0) ? std::min(Rp[i1], Rm[i0]) : std::min(Rp[i0], Rm[i1]);
        corr[i0] -= (Ar - al*Ar); corr[i1] += (Ar - al*Ar); }
    std::vector<double> q(n); double m = 0.0;
    for (int i = 0; i < n; ++i) { q[i] = (double)qH[i] - (dt/V)*corr[i]; m += q[i]*V; }
    *massOut = m; return q;
}

static Case run_case(const Ring& m, int profile, float cfl, int bdf, int nStep, int prelimit)
{
    const int n = m.n; const float h = m.h, u = 1.f, dt = cfl*h/u, V = h;
    std::vector<float> ro(n, 1.f), q(n, 0.f), phi(n), gx(n), lim(n, 1.f), res(n), td(n), mflux(n, u);
    for (int i = 0; i < n; ++i) { const float x = (i+0.5f)*h; q[i] = (profile == 0) ? expf(-powf((x-40.f*h)/(8.f*h), 2.f)) : ((i >= 20 && i < 40) ? 1.f : 0.f); }
    FctDev D(n, n, n, m.plane_cells, m.cpi, m.cp, m.nhp, m.vol, m.ccx, ro, 0);
    cudaMemcpy(D.dq, q.data(), n*sizeof(float), cudaMemcpyHostToDevice);
    geom_int *dpc=D.dpc,*dcpi=D.dcpi,*dcp=D.dcp,*dnhp=D.dnhp; geom_float *dvol=D.dvol,*dccx=D.dccx,*dccy=up(m.ccy),*dccz=up(m.ccz),*dpcx=up(m.pcx),*dpcy=up(m.pcy),*dpcz=up(m.pcz);
    float *dphi=up(phi),*dgx=up(gx),*dgy=up(std::vector<float>(n,0.f)),*dgz=up(std::vector<float>(n,0.f)),*dlim=up(lim),*dres=up(res),*dtd=up(td),*dmf=up(mflux),*dPface=up(std::vector<float>(n,0.f));
    std::vector<float*> hres={dres}, htd={dtd}; float **dres2=up(hres), **dtd2=up(htd);
    double* drn=up(std::vector<double>(1,0.0));
    Case C; std::vector<float> Pface_h, G_h, H_h;
    for (int step = 0; step < nStep; ++step) {
        cudaMemcpy(D.dqPP, D.dqP, n*sizeof(float), cudaMemcpyDeviceToDevice);
        cudaMemcpy(D.dqP, D.dq, n*sizeof(float), cudaMemcpyDeviceToDevice);
        const bool b2 = (bdf == 2 && step >= 1);
        const float a = b2 ? 1.5f : 1.f, b = b2 ? 2.f : 1.f, c = b2 ? 0.5f : 0.f;
        for (int k = 0; k < 600; ++k) {   // sub-iter (擬似時間 dτ = dt) を収束させる
            auto qh = down(D.dq, (size_t)n);
            for (int i = 0; i < n; ++i) phi[i] = qh[i]/ro[i];
            for (int i = 0; i < n; ++i) gx[i] = (phi[(i+1)%n]-phi[(i+n-1)%n])/(2.f*h);
            cudaMemcpy(dphi, phi.data(), n*sizeof(float), cudaMemcpyHostToDevice); cudaMemcpy(dgx, gx.data(), n*sizeof(float), cudaMemcpyHostToDevice);
            limiter_r1_scaled_d<<<(n+127)/128,128>>>(2, n, n, dpc, dcpi, dcp, dvol, dccx, dccy, dccz, dpcx, dpcy, dpcz, 1.0e-30f, dphi, dlim, dgx, dgy, dgz);
            face_recon_d<<<(n+127)/128,128>>>(n, n, dpc, dccx, dpcx, dphi, dgx, dlim, dmf, dPface);
            cudaMemset(dres, 0, n*sizeof(float)); cudaMemset(dtd, 0, n*sizeof(float)); cudaMemset(drn, 0, sizeof(double));
            species_advection_faceY_d<<<(n+127)/128,128>>>(n, n, dnhp, dpc, D.dro, dmf, 1, dPface, dres2, dtd2, 0, D.pq, 1);
            dual_update_d<<<(n+127)/128,128>>>(n, V, dt, dt, a, b, c, D.dqP, D.dqPP, dres, dtd, D.dq, drn);
            cudaDeviceSynchronize();
            double rn; cudaMemcpy(&rn, drn, sizeof(double), cudaMemcpyDeviceToHost);
            if (std::sqrt(rn) < 1e-10) break;
        }
        // 終了状態で P_face を再評価 (本番の assembleResidual 相当)
        { auto qh = down(D.dq, (size_t)n); for (int i = 0; i < n; ++i) phi[i] = qh[i]/ro[i]; for (int i = 0; i < n; ++i) gx[i] = (phi[(i+1)%n]-phi[(i+n-1)%n])/(2.f*h);
          cudaMemcpy(dphi, phi.data(), n*sizeof(float), cudaMemcpyHostToDevice); cudaMemcpy(dgx, gx.data(), n*sizeof(float), cudaMemcpyHostToDevice);
          limiter_r1_scaled_d<<<(n+127)/128,128>>>(2, n, n, dpc, dcpi, dcp, dvol, dccx, dccy, dccz, dpcx, dpcy, dpcz, 1.0e-30f, dphi, dlim, dgx, dgy, dgz);
          face_recon_d<<<(n+127)/128,128>>>(n, n, dpc, dccx, dpcx, dphi, dgx, dlim, dmf, dPface); cudaDeviceSynchronize(); }
        C.qH = down(D.dq, (size_t)n); Pface_h = down(dPface, (size_t)n); C.qP = down(D.dqP, (size_t)n);
        G_h = down(D.dG, (size_t)n); H_h = down(D.dH, (size_t)n);
        const bool useHist = b2;
        if (step == nStep-1) {   // (iii) α≡0
            D.run(dt, a, c, useHist, dmf, dPface, 0, prelimit, 1, 1e-8, 20000);
            C.qAlpha0 = down(D.dq, (size_t)n); C.qB = down(D.dqB, (size_t)n); C.qL = down(D.dqL, (size_t)n);
            cudaMemcpy(D.dq, C.qH.data(), n*sizeof(float), cudaMemcpyHostToDevice);
        }
        D.run(dt, a, c, useHist, dmf, dPface, 0, prelimit, 0, 1e-8, 20000);
        if (step == nStep-1) {
            double mh; double mnL, mxL;
            C.host = host_fct(n, V, dt, a, c, useHist, C.qH, C.qP, Pface_h, G_h, H_h, u, &mh, &mnL, &mxL);
            C.massHost = mh; C.minL = *std::min_element(C.qL.begin(), C.qL.end()); C.maxL = *std::max_element(C.qL.begin(), C.qL.end());
            C.massH = 0.0; for (int i = 0; i < n; ++i) C.massH += (double)C.qH[i]*V;
            C.qNew = down(D.dq, (size_t)n); C.massNew = 0.0; C.nModified = 0; C.maxMod = 0.f;
            for (int i = 0; i < n; ++i) { C.massNew += (double)C.qNew[i]*V; const float d = fabsf(C.qNew[i]-C.qH[i]); if (d > 1.0e-6f) { ++C.nModified; C.maxMod = fmaxf(C.maxMod, d); } }
        }
        D.finish(dt);
    }
    return C;
}

int main()
{
    Ring m(200, 1.0e-3f);
    const float V = m.h;
    printf("[a] smooth gaussian, CFL 0.5, BDF2 (flux-form history), 8 steps, no prelimit: unmodified cells are bit-identical, modified only near the extremum\n");
    {
        Case C = run_case(m, 0, 0.5f, 2, 8, 0);
        printf("   modified cells %d, max |q_C - q_H| %.3e (peak %.3f); mass q_H %.9e -> q_C %.9e (rel %.2e)\n", C.nModified, C.maxMod, *std::max_element(C.qH.begin(), C.qH.end()), C.massH, C.massNew, (C.massNew-C.massH)/C.massH);
        CHECK(C.nModified <= 6, "smooth case: FCT modified %d cells by >1e-6 (expected only the extremum neighbourhood)", C.nModified);
        CHECK(C.maxMod <= 5.0e-3f, "smooth case: modification %.3e too large (extremum clipping should stay <0.5%% of the peak)", C.maxMod);
        CHECK(std::fabs(C.massNew - C.massH) <= 1e-6*C.massH, "smooth case: mass not conserved");
    }
    for (int bdf = 1; bdf <= 2; ++bdf) {
        printf("[b] step profile, CFL 2, BDF%d, 6 steps: conservation / bounds / host double reference / alpha==0 identity\n", bdf);
        Case C = run_case(m, 1, 2.0f, bdf, 6, 0);
        float mnH = 1e30f, mxH = -1e30f, mnN = 1e30f, mxN = -1e30f;
        for (int i = 0; i < m.n; ++i) { mnH = fminf(mnH, C.qH[i]); mxH = fmaxf(mxH, C.qH[i]); mnN = fminf(mnN, C.qNew[i]); mxN = fmaxf(mxN, C.qNew[i]); }
        double dmax = 0.0, mhost = 0.0; for (int i = 0; i < m.n; ++i) { dmax = std::max(dmax, std::fabs(C.host[i] - (double)C.qNew[i])); mhost += C.host[i]*V; }
        int viol = 0;
        for (int i = 0; i < m.n; ++i) {
            const int l = (i+m.n-1)%m.n, r = (i+1)%m.n;
            const float mx = std::min(1.f, std::max({C.qP[i], C.qB[i], C.qP[l], C.qB[l], C.qP[r], C.qB[r]}));
            const float mn = std::max(0.f, std::min({C.qP[i], C.qB[i], C.qP[l], C.qB[l], C.qP[r], C.qB[r]}));
            if (C.qNew[i] > mx + 1e-6f || C.qNew[i] < mn - 1e-6f) ++viol;
        }
        double dA0 = 0.0, dBL = 0.0; for (int i = 0; i < m.n; ++i) { dA0 = std::max(dA0, std::fabs((double)C.qAlpha0[i] - (double)C.qB[i])); dBL = std::max(dBL, std::fabs((double)C.qB[i] - (double)C.qL[i])); }
        printf("   mass: q_H %.9e | q_C device %.9e (rel %.2e) | host double %.15e (rel to q_H %.2e)\n", C.massH, C.massNew, (C.massNew-C.massH)/C.massH, mhost, (mhost-C.massH)/C.massH);
        printf("   q_H in [%.4e, %.4e] -> q_C in [%.4e, %.4e]; q_L in [%.4e, %.4e]; local-bound violations %d; max|device-host| %.3e; alpha==0: |q_C-q_B| %.2e, |q_B-q_L| %.2e\n",
               mnH, mxH, mnN, mxN, C.minL, C.maxL, viol, dmax, dA0, dBL);
        CHECK(std::fabs(C.massNew - C.massH) <= 1e-6*C.massH, "BDF%d: device mass not conserved (rel %.3e)", bdf, (C.massNew-C.massH)/C.massH);
        CHECK(std::fabs(mhost - C.massH) <= 1e-12*C.massH, "BDF%d: host double reference not conserved (rel %.3e)", bdf, (mhost-C.massH)/C.massH);
        CHECK(mnN >= -1e-6f && mxN <= 1.f + 1e-6f, "BDF%d: q_C out of [0,1]: [%g, %g]", bdf, mnN, mxN);
        CHECK(viol == 0, "BDF%d: %d cells outside local extrema", bdf, viol);
        CHECK(dmax <= 1e-5, "BDF%d: device vs host double %.3e", bdf, dmax);
        CHECK(C.minL >= -1e-6f && C.maxL <= 1.f + 1e-6f, "BDF%d: low-order BE solution not bounded [%g, %g]", bdf, C.minL, C.maxL);
        CHECK(dA0 <= 1e-6, "BDF%d: alpha==0 does not give q_B (%.3e)", bdf, dA0);
        CHECK(dBL <= 1e-3, "BDF%d: base deviates from low-order solution by %.3e", bdf, dBL);
    }
    printf("[c] plan-4 M1 counterexample: 3-CV periodic, BE, CFL 1, hist [1,0,0], P_face [0,0.5,0]\n");
    {
        Ring r3(3, 1.0f); std::vector<float> ro(3, 1.f);
        FctDev D(3, 3, 3, r3.plane_cells, r3.cpi, r3.cp, r3.nhp, r3.vol, r3.ccx, ro, 0);
        std::vector<float> qH{1.f, -0.5f, 0.5f}, qP{1.f, 0.f, 0.f}, Pf{0.f, 0.5f, 0.f}, mf(3, 1.f);
        cudaMemcpy(D.dq, qH.data(), 3*sizeof(float), cudaMemcpyHostToDevice); cudaMemcpy(D.dqP, qP.data(), 3*sizeof(float), cudaMemcpyHostToDevice);
        float *dmf=up(mf), *dPf=up(Pf);
        D.run(1.f, 1.f, 0.f, false, dmf, dPf, 0, 0, 0, 1e-12, 2000);
        auto qL = down(D.dqL, 3), qB = down(D.dqB, 3), qC = down(D.dq, 3);
        printf("   q_L = [%.6f %.6f %.6f] (exact 4/7 2/7 1/7), q_B = [%.6f %.6f %.6f], q_C = [%.6f %.6f %.6f], mass %.6f\n", qL[0], qL[1], qL[2], qB[0], qB[1], qB[2], qC[0], qC[1], qC[2], qC[0]+qC[1]+qC[2]);
        CHECK(std::fabs(qL[0]-4.f/7) < 1e-5 && std::fabs(qL[1]-2.f/7) < 1e-5 && std::fabs(qL[2]-1.f/7) < 1e-5, "3-CV low-order solution wrong");
        CHECK(qC[0] >= -1e-6f && qC[1] >= -1e-6f && qC[2] >= -1e-6f && qC[0] <= 1.f+1e-6f, "3-CV: bounds violated after correction");
        CHECK(std::fabs(qC[0]+qC[1]+qC[2]-1.f) < 1e-6f, "3-CV: mass not conserved");
    }
    printf("[d] plan-5 M2 counterexample: non-periodic 3 CV (node), reverse-flow inlet mdot=-1 (own value in HO), internal +1, outlet +1, V/dt=0.5\n");
    {
        // planes: 0:(0,1) 1:(1,2) internal; 2:(0,3=ghostL) mdot -1; 3:(2,4=ghostR) mdot +1
        std::vector<geom_int> pc{0,1, 1,2, 0,3, 2,4}, cpi{0,2,4,6}, cp{0,2, 0,1, 1,3}, nhp{0,1,2,3};
        std::vector<geom_float> vol(3, 1.f), ccx{0.5f,1.5f,2.5f,-0.5f,3.5f}; std::vector<float> ro(5, 1.f);
        FctDev D(3, 4, 2, pc, cpi, cp, nhp, vol, ccx, ro, 2);
        std::vector<float> qH{0.9f, 0.7f, 0.5f}, qP{0.5f, 0.3f, 0.5f}, Pf{0.7f, 0.5f, 0.f, 0.f}, mf{1.f, 1.f, -1.f, 1.f};
        cudaMemcpy(D.dq, qH.data(), 3*sizeof(float), cudaMemcpyHostToDevice); cudaMemcpy(D.dqP, qP.data(), 3*sizeof(float), cudaMemcpyHostToDevice);
        float *dmf=up(mf), *dPf=up(Pf);
        D.run(2.f, 1.f, 0.f, false, dmf, dPf, 1, 0, 0, 1e-12, 2000);
        auto qL = down(D.dqL, 3), qB = down(D.dqB, 3), qC = down(D.dq, 3); auto st = down(D.dst, 8);
        // 低次 (M 行列, 流入は φ^n=0.5 の RHS): 0.5(qL0-0.5) = 0.5 - qL0 → qL0 = 0.5; 0.5(qL1-0.3) = qL0 - qL1 → 13/30; 0.5(qL2-0.5) = qL1 - qL2 → 41/90
        printf("   q_L = [%.6f %.6f %.6f] (exact 0.5, 13/30, 41/90), q_B = [%.6f %.6f %.6f], q_C = [%.6f %.6f %.6f], boundary exchange %.4e\n",
               qL[0], qL[1], qL[2], qB[0], qB[1], qB[2], qC[0], qC[1], qC[2], st[5]*2.0);
        CHECK(std::fabs(qL[0]-0.5f) < 1e-5 && std::fabs(qL[1]-13.f/30) < 1e-5 && std::fabs(qL[2]-41.f/90) < 1e-5, "boundary 3-CV low-order solution wrong");
        CHECK(qC[0] <= 0.5f + 1e-6f, "boundary 3-CV: CV0 exceeds its local max 0.5 (%g)", qC[0]);
        CHECK(qC[0] >= -1e-6f && qC[1] >= -1e-6f && qC[2] >= -1e-6f && qC[1] <= 1.f && qC[2] <= 1.f, "boundary 3-CV: bounds violated");
    }
    printf("[e] no-flux single CV, BDF2 history q^n=0.9, q^(n-1)=0.5\n");
    {
        std::vector<geom_int> pc, cpi{0,0}, cp, nhp; std::vector<geom_float> vol(1, 1.f), ccx{0.5f}; std::vector<float> ro(1, 1.f);
        FctDev D(1, 0, 0, pc, cpi, cp, nhp, vol, ccx, ro, 0);
        D.nPlanes = 0; D.nHalo = 0;
        std::vector<float> qH{31.f/30}, qP{0.9f}, qPP{0.5f};
        cudaMemcpy(D.dq, qH.data(), sizeof(float), cudaMemcpyHostToDevice); cudaMemcpy(D.dqP, qP.data(), sizeof(float), cudaMemcpyHostToDevice); cudaMemcpy(D.dqPP, qPP.data(), sizeof(float), cudaMemcpyHostToDevice);
        // 局所形の履歴 (G=0, H=(V/dt)(q^n-q^{n-1})=0.4) → 基点 31/30 → 物理限界逸脱 1/30 が記録される (制限不能)
        std::vector<float> H{0.4f}; cudaMemcpy(D.dH, H.data(), sizeof(float), cudaMemcpyHostToDevice);
        float* dmf=up(std::vector<float>(1,0.f)); float* dPf=up(std::vector<float>(1,0.f));
        D.run(1.f, 1.5f, 0.5f, true, dmf, dPf, 0, 0, 0, 1e-12, 100);
        auto qC = down(D.dq, 1), qB = down(D.dqB, 1); auto st = down(D.dst, 8);
        printf("   local-form history: q_C %.6f (unchanged), clipped base %.6f, recorded base violation %.6f (expected 1/30=0.0333)\n", qC[0], qB[0], st[4]);
        CHECK(std::fabs(st[4] - 1.0/30) < 1e-4, "no-flux BDF2: base violation not recorded (%g)", st[4]);
        CHECK(std::fabs(qC[0] - 31.f/30) < 1e-6f, "no-flux BDF2: FCT must not change a fluxless cell (floor handles it)");
    }
    printf(g_fail ? "FAILED (%d)\n" : "ALL PASS\n", g_fail);
    return g_fail ? 1 : 0;
}
