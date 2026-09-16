// test_passive_scalar.cu — 受動種 (トレーサ・凝縮モーメント) の化学種経路の単体試験
//   (plans/active/species-passive-scalar-unification.md §4.1 / §6-2; codex plan-2 M1/M2)。
//   build: nvcc --expt-relaxed-constexpr -arch=sm_86 -I. -o test_passive_scalar tests/unit/test_passive_scalar.cu
//   (a) 無次元化 Venkat リミッタ limiter_r1_scaled_d: g/Q2/Q1/Q0 の実測スケール (1e-2 / 2e1 / 1e8 / 9e14)・ゼロ近傍・極値で
//       ψ が有限かつ [0,1]、スケール不変 (ψ(κφ)=ψ(φ))、単調域で ψ≈1・極値で ψ<1。対照: 旧式 (3 次積) は Q0 スケールで Inf/NaN。
//   (b) 1 次元ステップ状分布の移流 (一様 ρ=1, u=1; S3 = MUSCL2 + ψ_P, 1 次 = 風上) を species_advection_faceY_d + point-implicit
//       更新 + passive_bounds_d で回し: 保存 (∫ρφ の変化 = 記録した補正量; 補正前は 1e-6 で保存)、更新確定後 0≤ρξ≤ρ、補正収支の
//       記録 (lo/hi/abs/total) の整合、S3 の遷移幅 (L1 誤差) が 1 次より小さいこと。
//   (c) passive_bounds_d の収支の厳密性 (既知の負値・超過値で lo/hi/abs/total が解析値と一致)。
#include <cstdio>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/limiterFunctions_d.cuh"
#include "cuda_forge/passiveLimiter_d.cuh"
#include "cuda_forge/passiveKernels_d.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)
template <class T> static T* up(const std::vector<T>& v) { T* d = nullptr; cudaMalloc((void**)&d, v.size()*sizeof(T)); cudaMemcpy(d, v.data(), v.size()*sizeof(T), cudaMemcpyHostToDevice); return d; }
template <class T> static std::vector<T> down(const T* d, size_t n) { std::vector<T> v(n); cudaMemcpy(v.data(), d, n*sizeof(T), cudaMemcpyDeviceToHost); return v; }

// 1 次元セル鎖 (h 間隔, nCells 実セル + 両端 ghost)。plane ip: 内部 ip<nCells-1 は (ip, ip+1)、境界 2 面は (0, ghostL), (nCells-1, ghostR)。
struct Chain {
    int n; float h;
    std::vector<geom_int> plane_cells, cpi, cp, nhp;
    std::vector<geom_float> vol, ccx, ccy, ccz, pcx, pcy, pcz, fx;
    int nNormal, nPlanes;
    Chain(int n_, float h_) : n(n_), h(h_) {
        nNormal = n-1; nPlanes = n+1;
        plane_cells.assign(2*nPlanes, 0);
        for (int ip = 0; ip < nNormal; ++ip) { plane_cells[2*ip] = ip; plane_cells[2*ip+1] = ip+1; }
        plane_cells[2*nNormal] = 0;   plane_cells[2*nNormal+1] = n;     // 左 ghost = n
        plane_cells[2*nNormal+2] = n-1; plane_cells[2*nNormal+3] = n+1; // 右 ghost = n+1
        cpi.assign(n+1, 0); cp.clear();
        for (int ic = 0; ic < n; ++ic) {
            cpi[ic] = (geom_int)cp.size();
            if (ic > 0) cp.push_back(ic-1); else cp.push_back(nNormal);
            if (ic < n-1) cp.push_back(ic); else cp.push_back(nNormal+1);
        }
        cpi[n] = (geom_int)cp.size();
        for (int ip = 0; ip < nPlanes; ++ip) nhp.push_back(ip);
        vol.assign(n+2, h); ccx.resize(n+2); ccy.assign(n+2, 0.f); ccz.assign(n+2, 0.f);
        for (int ic = 0; ic < n; ++ic) ccx[ic] = (ic+0.5f)*h;
        ccx[n] = -0.5f*h; ccx[n+1] = (n+0.5f)*h;
        pcx.resize(nPlanes); pcy.assign(nPlanes, 0.f); pcz.assign(nPlanes, 0.f); fx.assign(nPlanes, 0.5f);
        for (int ip = 0; ip < nNormal; ++ip) pcx[ip] = (ip+1)*h;
        pcx[nNormal] = 0.f; pcx[nNormal+1] = n*h;
    }
};

// (a) リミッタ: 指定分布 Q (nCells 実セル) と勾配 (中心差分) で ψ を計算。
static std::vector<float> run_limiter(const Chain& m, const std::vector<float>& Qh, int scheme, float volume_eps)
{
    std::vector<float> Q(m.n+2, 0.f), gx(m.n+2, 0.f), gy(m.n+2, 0.f), gz(m.n+2, 0.f), lim(m.n+2, 1.f);
    for (int i = 0; i < m.n; ++i) Q[i] = Qh[i];
    for (int i = 0; i < m.n; ++i) { const float ql = (i>0)?Q[i-1]:Q[i], qr = (i<m.n-1)?Q[i+1]:Q[i]; gx[i] = (qr-ql)/((i>0&&i<m.n-1)?2.f*m.h:m.h); }
    std::vector<geom_float> vol(m.n+2, volume_eps);
    geom_int *dpc=up(m.plane_cells),*dcpi=up(m.cpi),*dcp=up(m.cp);
    geom_float *dvol=up(vol),*dccx=up(m.ccx),*dccy=up(m.ccy),*dccz=up(m.ccz),*dpcx=up(m.pcx),*dpcy=up(m.pcy),*dpcz=up(m.pcz);
    float *dQ=up(Q),*dgx=up(gx),*dgy=up(gy),*dgz=up(gz),*dlim=up(lim);
    limiter_r1_scaled_d<<<(m.n+127)/128,128>>>(scheme, m.n, m.nNormal, dpc, dcpi, dcp, dvol, dccx, dccy, dccz, dpcx, dpcy, dpcz,
        1.0e-30f, dQ, dlim, dgx, dgy, dgz);
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; }
    auto out = down(dlim, (size_t)m.n);
    cudaFree(dpc); cudaFree(dcpi); cudaFree(dcp); cudaFree(dvol); cudaFree(dccx); cudaFree(dccy); cudaFree(dccz); cudaFree(dpcx); cudaFree(dpcy); cudaFree(dpcz);
    cudaFree(dQ); cudaFree(dgx); cudaFree(dgy); cudaFree(dgz); cudaFree(dlim);
    return out;
}

// 旧式 (limiter_r1_d と同じ差分をそのまま Venkat に渡す) を device で評価 → Q0 スケールで Inf/NaN になる対照。
__global__ void venkat_raw_probe_d(float delta_p, float delta_m, float volume, float* out)
{
    out[0] = venkata_limiter(delta_p, -delta_p, delta_m, volume);
}

static void test_limiter()
{
    printf("[a] scaled Venkat limiter finiteness / scale invariance\n");
    const int n = 32; const float h = 1.0e-3f;   // 双対 CV 体積 ~1e-9 m3 相当 (case/44 の最小 4.5e-8 より厳しい側)
    Chain m(n, h);
    // 実測スケール (case/44 run_0170 res_24000: g 1.6e-2, Q2 2.2e1, Q1 1.3e8, Q0 8.7e14) + 極小 (1e-25; φ_floor 1e-30 より上) + 極大 (1e30;
    // 勾配 φ/h ~ 1e33 が float に収まる範囲)。φ_floor 以下 (1e-30) は下の別試験で有限性だけを見る。
    const double scales[] = {1.0, 1.6e-2, 2.2e1, 1.3e8, 8.7e14, 1.0e-25, 1.0e30};
    // 分布: 滑らかな山 + ステップ + ゼロ帯 (極値・単調・定数域を含む)
    std::vector<float> base(n);
    for (int i = 0; i < n; ++i) {
        if (i < 8) base[i] = 0.f;
        else if (i < 16) base[i] = (float)(0.5*(1.0 - cos(3.14159265*(i-8)/8.0)));
        else if (i < 24) base[i] = 1.f;
        else base[i] = (float)(0.2 + 0.1*sin(i));
    }
    std::vector<float> ref;
    for (double sc : scales) {
        std::vector<float> Q(n); for (int i = 0; i < n; ++i) Q[i] = (float)(base[i]*sc);
        for (int scheme : {2, 1}) {
            auto psi = run_limiter(m, Q, scheme, (float)(h*h*h));
            for (int i = 0; i < n; ++i) CHECK(std::isfinite(psi[i]) && psi[i] >= 0.f && psi[i] <= 1.f, "scale %.1e scheme %d cell %d psi=%g not finite/in [0,1]", sc, scheme, i, psi[i]);
            if (scheme == 2) {
                if (sc == 1.0) ref = psi;
                else {
                    float dmax = 0.f; for (int i = 0; i < n; ++i) dmax = fmaxf(dmax, fabsf(psi[i]-ref[i]));
                    CHECK(dmax < 1e-5f, "scale %.1e: psi differs from scale 1 by %g", sc, dmax);   // ε²=V は無次元差分に対し無視できる
                    printf("   scale %.1e: max|psi - psi(scale 1)| = %.2e\n", sc, dmax);
                }
            }
        }
    }
    // 定性: 定数域 (i=2, 20) は ψ=1、極値 (i=15,16 付近の山頂・ステップ肩) は ψ<1
    CHECK(ref[2] == 1.f && ref[20] == 1.f, "constant region psi != 1 (%g %g)", ref[2], ref[20]);
    CHECK(ref[16] < 0.999f || ref[15] < 0.999f, "extremum psi not limited (%g %g)", ref[15], ref[16]);
    // ゼロ場: ψ=1 (φ_ref=floor, 差分 0)
    { std::vector<float> Z(n, 0.f); auto psi = run_limiter(m, Z, 2, (float)(h*h*h)); for (int i = 0; i < n; ++i) CHECK(psi[i] == 1.f, "zero field psi[%d]=%g", i, psi[i]); }
    // φ_floor 以下のスケール (1e-30·base): φ_ref が floor に張り付き ψ はスケール 1 と一致しないが有限・[0,1] であること
    { std::vector<float> Q(n); for (int i = 0; i < n; ++i) Q[i] = (float)(base[i]*1.0e-30); auto psi = run_limiter(m, Q, 2, (float)(h*h*h));
      for (int i = 0; i < n; ++i) CHECK(std::isfinite(psi[i]) && psi[i] >= 0.f && psi[i] <= 1.f, "floor-regime psi[%d]=%g", i, psi[i]); }
    // 対照: 旧式は Q0 スケールの差分 (1e14, 1e13) で分子・分母が Inf → 商 NaN
    { float* d = nullptr; cudaMalloc((void**)&d, sizeof(float));
      venkat_raw_probe_d<<<1,1>>>(1.0e14f, 1.0e13f, 1.0e-9f, d); cudaDeviceSynchronize();
      float v; cudaMemcpy(&v, d, sizeof(float), cudaMemcpyDeviceToHost);
      printf("   reference: raw venkata_limiter(1e14, 1e13) = %g (expected non-finite: the unscaled form overflows at Q0 scale)\n", v);
      CHECK(!std::isfinite(v), "raw Venkat at Q0 scale is finite (%g): the scaled variant would not be needed", v);
      // 同じ差分を無次元化 (÷1e14) すると有限
      venkat_raw_probe_d<<<1,1>>>(1.0f, 0.1f, 1.0e-9f, d); cudaDeviceSynchronize(); cudaMemcpy(&v, d, sizeof(float), cudaMemcpyDeviceToHost);
      CHECK(std::isfinite(v) && v > 0.f, "scaled Venkat(1, 0.1) = %g", v); cudaFree(d); }   // 関数自体は >1 を返し得る (カーネル側で [0,1] にクリップ)
}

// (b) 移流試験用: point-implicit forward-Euler 更新 (runge_kutta_exp_scalar_d の coef 1/0/1・floor なし版) と 1 次面値・S3 面値。
__global__ void pi_update_d(int n, const float* rophiN, const float* res, const float* td, float dt, float vol, float* rophi)
{
    const int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < n) { const float fac = 1.f + dt*(td[ic]/vol); rophi[ic] = rophiN[ic] + (res[ic]*dt/vol)/fac; }
}
// MUSCL 2 次 (interp_MUSCL_2nd と同式: φ_f = φ_C + ψ ∇φ·(x_f − x_C)) で L/R を作り、下限 0・上限 1 でクリップし upwind を書く。
__global__ void face_recon_d(int nPlanes, int nCells, const geom_int* plane_cells, const geom_float* ccx, const geom_float* pcx,
                             const float* phi, const float* gx, const float* lim, const float* massflux, int order, float* Pface)
{
    const int ip = blockDim.x*blockIdx.x + threadIdx.x;
    if (ip >= nPlanes) return;
    const int ic0 = plane_cells[2*ip], ic1 = plane_cells[2*ip+1];
    float pl = phi[ic0], pr = phi[ic1];
    if (order == 2 && ic1 < nCells) {
        pl = phi[ic0] + lim[ic0]*gx[ic0]*(pcx[ip]-ccx[ic0]);
        pr = phi[ic1] + lim[ic1]*gx[ic1]*(pcx[ip]-ccx[ic1]);
    }
    pl = fminf(fmaxf(pl, 0.f), 1.f); pr = fminf(fmaxf(pr, 0.f), 1.f);
    Pface[ip] = (massflux[ip] >= 0.f) ? pl : pr;
}

struct AdvResult { std::vector<float> phi; double l1; int width; double mass0, massRaw, massFinal; double lo, hi, ab, tot; float minphi, maxphi; };

static AdvResult run_advection(const Chain& m, int order, int nStep, float cfl)
{
    const int n = m.n; const float h = m.h, u = 1.f, dt = cfl*h/u;
    std::vector<float> ro(n+2, 1.f), rophi(n+2, 0.f), rophiN(n+2, 0.f), phi(n+2, 0.f), gx(n+2, 0.f), lim(n+2, 1.f), res(n, 0.f), td(n, 0.f), corr(n, 0.f);
    std::vector<float> mflux(m.nPlanes, u*1.f);   // ρuA, A=1 (面積 1, V = h·1)
    mflux[m.nNormal] = -u;   // 左境界面 (0, ghostL) は ic0→ic1 が外向きなので流入は負
    for (int i = 8; i < 16; ++i) rophi[i] = 1.f;   // ステップ (ξ=1) を左寄りに置く (右端に達しない範囲で回す)
    double mass0 = 0.0; for (int i = 0; i < n; ++i) mass0 += rophi[i]*m.vol[i];
    geom_int *dpc=up(m.plane_cells),*dcpi=up(m.cpi),*dcp=up(m.cp),*dnhp=up(m.nhp);
    geom_float *dvol=up(m.vol),*dccx=up(m.ccx),*dccy=up(m.ccy),*dccz=up(m.ccz),*dpcx=up(m.pcx),*dpcy=up(m.pcy),*dpcz=up(m.pcz);
    float *dro=up(ro),*drophi=up(rophi),*drophiN=up(rophiN),*dphi=up(phi),*dgx=up(gx),*dgy=up(std::vector<float>(n+2,0.f)),*dgz=up(std::vector<float>(n+2,0.f)),*dlim=up(lim),*dres=up(res),*dtd=up(td),*dcorr=up(corr),*dmf=up(mflux);
    float *dPface=up(std::vector<float>(m.nPlanes,0.f));
    std::vector<float*> hres={dres}, htd={dtd}, hro={drophi}; float **dres2=up(hres), **dtd2=up(htd), **drophi2=up(hro);
    std::vector<double> stats(4,0.0); double* dstats=up(stats);
    double massRaw = mass0;
    for (int it = 0; it < nStep; ++it) {
        cudaMemcpy(drophiN, drophi, (n+2)*sizeof(float), cudaMemcpyDeviceToDevice);
        // primitive (ghost = 0: 入口 φ=0 Dirichlet, 出口 Neumann はステップが届かないので 0 のまま)
        auto rh = down(drophi, (size_t)n+2); for (int i = 0; i < n+2; ++i) phi[i] = rh[i]/ro[i];
        for (int i = 0; i < n; ++i) { const float ql=(i>0)?phi[i-1]:phi[n], qr=(i<n-1)?phi[i+1]:phi[n+1]; gx[i]=(qr-ql)/(2.f*h); }
        cudaMemcpy(dphi, phi.data(), (n+2)*sizeof(float), cudaMemcpyHostToDevice); cudaMemcpy(dgx, gx.data(), (n+2)*sizeof(float), cudaMemcpyHostToDevice);
        if (order == 2) limiter_r1_scaled_d<<<(n+127)/128,128>>>(2, n, m.nNormal, dpc, dcpi, dcp, dvol, dccx, dccy, dccz, dpcx, dpcy, dpcz, 1.0e-30f, dphi, dlim, dgx, dgy, dgz);
        face_recon_d<<<(m.nPlanes+127)/128,128>>>(m.nPlanes, n, dpc, dccx, dpcx, dphi, dgx, dlim, dmf, order, dPface);
        cudaMemset(dres, 0, n*sizeof(float)); cudaMemset(dtd, 0, n*sizeof(float));
        species_advection_faceY_d<<<(m.nPlanes+127)/128,128>>>(n, m.nPlanes, dnhp, dpc, dro, dmf, 1, dPface, dres2, dtd2, 0, drophi2, 1);
        pi_update_d<<<(n+127)/128,128>>>(n, drophiN, dres, dtd, dt, h, drophi);
        cudaDeviceSynchronize();
        { auto r = down(drophi, (size_t)n); double s = 0.0; for (int i = 0; i < n; ++i) s += r[i]*m.vol[i]; massRaw = s; }   // 補正前の総量 (最終 step のみ使う)
        cudaMemset(dstats+3, 0, sizeof(double));
        passive_bounds_d<<<(n+127)/128,128>>>(n, drophi, 1, dro, dvol, dcorr, dstats, nullptr);
        cudaDeviceSynchronize();
    }
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; }
    AdvResult R; R.phi = down(drophi, (size_t)n); auto st = down(dstats, (size_t)4);
    R.lo = st[0]; R.hi = st[1]; R.ab = st[2]; R.tot = st[3]; R.mass0 = mass0; R.massRaw = massRaw;
    R.massFinal = 0.0; R.minphi = 1e30f; R.maxphi = -1e30f;
    for (int i = 0; i < n; ++i) { R.massFinal += R.phi[i]*m.vol[i]; R.minphi = fminf(R.minphi, R.phi[i]); R.maxphi = fmaxf(R.maxphi, R.phi[i]); }
    // 解析解 (平行移動したステップ) との L1 誤差。point-implicit forward-Euler は増分を 1/(1+dt·td/V) = 1/(1+CFL) に減衰する
    // 擬似時間スキームなので、実効の移動量は nStep·CFL·h/(1+CFL) (1 次・S3 とも対角は同じ 1 次風上)。
    const double shift = nStep*dt*u/(1.0 + cfl); R.l1 = 0.0; R.width = 0;
    for (int i = 0; i < n; ++i) { const double x = (i+0.5)*h; const double ex = (x >= 8*h+shift && x < 16*h+shift) ? 1.0 : 0.0; R.l1 += fabs(R.phi[i]-ex)*h;
        if (R.phi[i] > 0.1f && R.phi[i] < 0.9f) ++R.width; }
    cudaFree(dpc); cudaFree(dcpi); cudaFree(dcp); cudaFree(dnhp); cudaFree(dvol); cudaFree(dccx); cudaFree(dccy); cudaFree(dccz); cudaFree(dpcx); cudaFree(dpcy); cudaFree(dpcz);
    cudaFree(dro); cudaFree(drophi); cudaFree(drophiN); cudaFree(dphi); cudaFree(dgx); cudaFree(dgy); cudaFree(dgz); cudaFree(dlim); cudaFree(dres); cudaFree(dtd); cudaFree(dcorr); cudaFree(dmf); cudaFree(dPface);
    cudaFree(dres2); cudaFree(dtd2); cudaFree(drophi2); cudaFree(dstats);
    return R;
}

static void test_advection()
{
    printf("[b] 1-D step advection: conservation / bounds / correction budget / S3 sharper than 1st order\n");
    const int n = 200; const float h = 1.0e-3f; Chain m(n, h);
    const int nStep = 100; const float cfl = 0.5f;   // 移動 50 セル (右端 200 には届かない)
    AdvResult r1 = run_advection(m, 1, nStep, cfl);
    AdvResult r2 = run_advection(m, 2, nStep, cfl);
    for (int order = 1; order <= 2; ++order) {
        const AdvResult& R = (order == 1) ? r1 : r2;
        const double budget = R.lo + R.hi;   // 符号付き補正の和 = 補正が総量に与えた変化 (最終 step 分は massRaw→massFinal)
        printf("   order %d: mass0 %.9e final %.9e (raw before last floor %.9e) | corr lo %.3e hi %.3e abs %.3e total %.9e | phi in [%.3e, %.3e] | L1 %.4e width(0.1<phi<0.9) %d cells\n",
               order, R.mass0, R.massFinal, R.massRaw, R.lo, R.hi, R.ab, R.tot, R.minphi, R.maxphi, R.l1, R.width);
        // 保存: 補正の符号付き総和で総量の変化を説明できる (面クリップは保存的、更新は保存的)。
        CHECK(fabs((R.massFinal - R.mass0) - budget) <= 1e-6*R.mass0, "order %d: mass change %.3e != recorded corrections %.3e", order, R.massFinal-R.mass0, budget);
        CHECK(fabs(R.tot - R.massFinal) <= 1e-6*R.mass0, "order %d: stats total %.9e != final mass %.9e", order, R.tot, R.massFinal);
        CHECK(R.minphi >= 0.f && R.maxphi <= 1.f, "order %d: phi out of [0,1] after commit (%g, %g)", order, R.minphi, R.maxphi);
        CHECK(R.ab >= fabs(budget) - 1e-12, "order %d: abs corr %.3e < |signed| %.3e", order, R.ab, budget);
        if (order == 1) CHECK(R.ab == 0.0, "1st-order upwind at CFL 0.5 should need no floor correction (abs %.3e)", R.ab);
    }
    CHECK(r2.l1 < 0.75*r1.l1, "S3 not sharper: L1 S3 %.4e vs 1st %.4e", r2.l1, r1.l1);
    CHECK(r2.width < r1.width, "S3 transition not narrower: %d vs %d cells", r2.width, r1.width);
    printf("   L1 error: 1st %.4e, S3 %.4e (ratio %.3f); transition width 1st %d, S3 %d cells\n", r1.l1, r2.l1, r2.l1/r1.l1, r1.width, r2.width);
}

// (c) passive_bounds_d の収支: 既知の負値・超過値
static void test_bounds_budget()
{
    printf("[c] passive_bounds_d budget on known values\n");
    const int n = 300; std::vector<float> ro(n, 2.f), v(n), corr(n, 0.f); std::vector<geom_float> vol(n, 0.5f);
    for (int i = 0; i < n; ++i) v[i] = (i % 3 == 0) ? -0.25f : ((i % 3 == 1) ? 2.5f : 1.0f);   // 100 個ずつ: 負 / ρ 超過 / 正常
    float *dro=up(ro),*dv=up(v),*dcorr=up(corr); geom_float* dvol=up(vol); std::vector<double> st(4,0.0); double* dst=up(st);
    passive_bounds_d<<<(n+127)/128,128>>>(n, dv, 1, dro, dvol, dcorr, dst, nullptr); cudaDeviceSynchronize();
    auto s = down(dst, (size_t)4); auto out = down(dv, (size_t)n); auto c = down(dcorr, (size_t)n);
    const double lo = 100*0.25*0.5, hi = -100*0.5*0.5, ab = lo - hi, tot = (100*0.0 + 100*2.0 + 100*1.0)*0.5;
    CHECK(fabs(s[0]-lo) < 1e-9 && fabs(s[1]-hi) < 1e-9 && fabs(s[2]-ab) < 1e-9 && fabs(s[3]-tot) < 1e-9, "stats lo %g hi %g abs %g tot %g (expected %g %g %g %g)", s[0], s[1], s[2], s[3], lo, hi, ab, tot);
    for (int i = 0; i < n; ++i) {
        const float ex = (i%3==0) ? 0.f : ((i%3==1) ? 2.f : 1.f), ec = (i%3==0) ? 0.25f : ((i%3==1) ? 0.5f : 0.f);
        CHECK(out[i] == ex && c[i] == ec, "cell %d: out %g corr %g (expected %g %g)", i, out[i], c[i], ex, ec);
    }
    // モーメント (上限なし): 超過値はそのまま
    cudaMemcpy(dv, v.data(), n*sizeof(float), cudaMemcpyHostToDevice); cudaMemset(dst, 0, 4*sizeof(double)); cudaMemset(dcorr, 0, n*sizeof(float));
    passive_bounds_d<<<(n+127)/128,128>>>(n, dv, 0, dro, dvol, dcorr, dst, nullptr); cudaDeviceSynchronize();
    s = down(dst, (size_t)4); out = down(dv, (size_t)n);
    CHECK(fabs(s[0]-lo) < 1e-9 && s[1] == 0.0 && fabs(s[2]-lo) < 1e-9, "moment stats lo %g hi %g abs %g", s[0], s[1], s[2]);
    CHECK(out[1] == 2.5f, "moment upper bound must not apply (out %g)", out[1]);
    printf("   stats: lo %.4f hi %.4f abs %.4f total %.4f (tracer) — OK\n", lo, hi, ab, tot);
    cudaFree(dro); cudaFree(dv); cudaFree(dcorr); cudaFree(dvol); cudaFree(dst);
}

// (d) 流れの密度更新と整合した増分 (plan §5.1 #19): z=0 のとき ρφ = φ_N ρ_new で φ は不変、増分制限は基点 φ_N ρ_new に対して掛かる。
static void test_rho_term()
{
    printf("[d] passive_add_rho_term_d: z=0 keeps phi, limiter base is phi_N*rho_new\n");
    const int n = 256; std::vector<float> roPre(n), ro(n), N(n), cand(n), vol(n, 1.f), lim(n, 0.f);
    for (int i = 0; i < n; ++i) { roPre[i] = 1.0f + 0.3f*sinf(0.1f*i); ro[i] = roPre[i]*(1.0f + 0.2f*cosf(0.05f*i)); N[i] = roPre[i]*(0.1f + 0.8f*(i%7)/6.0f); cand[i] = N[i]; }
    float *dpre=up(roPre),*dro=up(ro),*dN=up(N),*dc=up(cand),*dlim=up(lim); geom_float* dvol=up(vol); std::vector<double> st(8,0.0); double* dst=up(st); int one=1000000000; int* dth=nullptr; cudaMalloc((void**)&dth,sizeof(int)); cudaMemcpy(dth,&one,sizeof(int),cudaMemcpyHostToDevice);
    passive_add_rho_term_d<<<(n+127)/128,128>>>(n, dc, dN, dpre, dro); cudaDeviceSynchronize();
    auto out = down(dc, (size_t)n); float worst = 0.f;
    for (int i = 0; i < n; ++i) worst = fmaxf(worst, fabsf(out[i]/ro[i] - N[i]/roPre[i]));
    CHECK(worst < 2e-7f, "z=0: phi changed by %g", worst);
    printf("   z=0: max|phi_new - phi_N| = %.2e\n", worst);
    // z=0 の候補に対する増分制限 (基点 φ_N ρ_new) は無作用
    passive_limit_increment_d<<<(n+127)/128,128>>>(n, dc, dN, 1, dro, dvol, dlim, dst, dth, nullptr, dpre); cudaDeviceSynchronize();
    auto s2 = down(dst, (size_t)8); auto out2 = down(dc, (size_t)n); float dd = 0.f; for (int i = 0; i < n; ++i) dd = fmaxf(dd, fabsf(out2[i]-out[i]));
    CHECK(s2[0] == 0.0 && s2[1] == 0.0 && dd == 0.f, "limiter acted on z=0 candidate (amount %g cells %g maxΔ %g)", s2[0], s2[1], dd);
    // z>0 で上限超過: 基点 φ_N ρ_new から ρ_new までが allowed
    for (int i = 0; i < n; ++i) cand[i] = N[i]/roPre[i]*ro[i] + 0.5f*ro[i];   // z = 0.5 ρ_new (多くのセルで超過)
    cudaMemcpy(dc, cand.data(), n*sizeof(float), cudaMemcpyHostToDevice); cudaMemset(dst, 0, 8*sizeof(double));
    passive_limit_increment_d<<<(n+127)/128,128>>>(n, dc, dN, 1, dro, dvol, dlim, dst, dth, nullptr, dpre); cudaDeviceSynchronize();
    out2 = down(dc, (size_t)n); int bad = 0; for (int i = 0; i < n; ++i) if (out2[i] > ro[i]*(1.f+1e-6f) || out2[i] < 0.f) ++bad;
    CHECK(bad == 0, "%d cells out of [0,rho] after limiting with rho term", bad);
    printf("   z=0.5rho: all cells within [0,rho] after theta_b (limited cells %.0f)\n", down(dst,(size_t)8)[1]);
    cudaFree(dpre); cudaFree(dro); cudaFree(dN); cudaFree(dc); cudaFree(dlim); cudaFree(dvol); cudaFree(dst); cudaFree(dth);
}

int main()
{
    test_limiter();
    test_advection();
    test_bounds_budget();
    test_rho_term();
    printf(g_fail ? "FAILED (%d)\n" : "ALL PASS\n", g_fail);
    return g_fail ? 1 : 0;
}
