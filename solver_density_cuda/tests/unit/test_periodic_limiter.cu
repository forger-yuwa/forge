// test_periodic_limiter.cu — 周期 node (合併 CV) の 2 段リミッタと max/min gather の単体試験
//   (plans/active/species-passive-scalar-unification.md §4.8; codex result-2 M3)。
//   build: nvcc --expt-relaxed-constexpr -arch=sm_86 -I. -o test_periodic_limiter tests/unit/test_periodic_limiter.cu
//   (a) periodicGatherMax1ToRoot_d / Min (CAS 版 float atomic): 混合符号 (-2 < -1, ±0, ±1e30, ±Inf) を group 2/4/8 (面/辺/角) で
//       root の max/min が host 参照と一致し、broadcast 後に member も同値、非周期 (root==self) は不変、NaN member は skip。
//   (b) 巻き付き 1 次元周期鎖 (n セル) を「継ぎ目で割った鎖 (cell 0 と cell n が同じ周期 group、各 1 面の部分 CV)」として
//       2 段 (limiter_extrema_d → group max/min → limiter_psi_merged_d → group min) で解くと、同じ周期データを回転して継ぎ目を
//       内部に置いた鎖に 1 段 kernel (limiter_r1_scaled_d / 2 段の identity-root 版) を掛けた結果とビット一致すること
//       (SCALED=true/false、Venkat/BJ、体積 2 通り)。周期対 (cell 0 / n) の ψ は同値。
//   (c) 対照: 1 段 kernel を割った鎖にそのまま掛けると周期対の ψ は一致しない (= 修正対象の症状が試験で再現する)。
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstring>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/limiterFunctions_d.cuh"
#include "cuda_forge/passiveLimiter_d.cuh"
#include "cuda_forge/limiterPeriodic_d.cuh"
#include "cuda_forge/periodicAtomic_d.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)
template <class T> static T* up(const std::vector<T>& v) { T* d = nullptr; cudaMalloc((void**)&d, v.size()*sizeof(T)); cudaMemcpy(d, v.data(), v.size()*sizeof(T), cudaMemcpyHostToDevice); return d; }
template <class T> static std::vector<T> down(const T* d, size_t n) { std::vector<T> v(n); cudaMemcpy(v.data(), d, n*sizeof(T), cudaMemcpyDeviceToHost); return v; }
static bool same_bits(float a, float b) { return std::memcmp(&a, &b, sizeof(float)) == 0; }

// periodicNode_d.cu の broadcast と同体 (試験 TU は本体 .cu をリンクしない)。
__global__ void bcast_d(geom_int n, geom_int* root, flow_float* a)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < n) { const geom_int r = root[ic]; if (r != ic) a[ic] = a[r]; }
}

static void gather_max(geom_int n, geom_int* root, flow_float* a) { periodicGatherMax1ToRoot_d<<<(n+127)/128,128>>>(n, root, a); bcast_d<<<(n+127)/128,128>>>(n, root, a); }
static void gather_min(geom_int n, geom_int* root, flow_float* a) { periodicGatherMin1ToRoot_d<<<(n+127)/128,128>>>(n, root, a); bcast_d<<<(n+127)/128,128>>>(n, root, a); }

// ---------------------------------------------------------------------------------------------------------------
static void test_gather()
{
    printf("[a] float atomic max/min gather (mixed sign; groups of 2 / 4 / 8 = face / edge / corner)\n");
    const int n = 32;
    std::vector<geom_int> root(n); for (int i = 0; i < n; ++i) root[i] = i;
    // group A (面, 2): {0,5}       群 B (辺, 4): {1,7,9,12}      群 C (角, 8): {2,3,4,6,8,10,11,13}
    // group D (面, 2): {14,15} (-2 < -1)   群 E (辺, 4): {16,17,18,19} (-0/+0 と負)   群 F (角, 8): {20..27} (負のみ, ±Inf)
    root[5] = 0;
    for (int i : {7, 9, 12}) root[i] = 1;
    for (int i : {3, 4, 6, 8, 10, 11, 13}) root[i] = 2;
    root[15] = 14;
    for (int i : {17, 18, 19}) root[i] = 16;
    for (int i = 21; i <= 27; ++i) root[i] = 20;
    const float inf = INFINITY;
    std::vector<float> v(n, 0.f);
    v[0] = -3.5f; v[5] = -0.0f;                                                            // A: max -0 (==0), min -3.5
    v[1] = 1.0e30f; v[7] = -1.0e30f; v[9] = 0.25f; v[12] = -0.25f;                          // B: max 1e30, min -1e30
    v[2] = 0.5f; v[3] = -0.5f; v[4] = 2.0f; v[6] = -2.0f; v[8] = 1.5f; v[10] = -1.5f; v[11] = 0.0f; v[13] = -1.0f;   // C: 2 / -2
    v[14] = -1.0f; v[15] = -2.0f;                                                          // D: max -1, min -2 (-2 < -1)
    v[16] = -0.0f; v[17] = 0.0f; v[18] = -1.0e-30f; v[19] = 1.0e-30f;                       // E: max 1e-30, min -1e-30
    v[20] = -4.0f; v[21] = -3.0f; v[22] = -inf; v[23] = -5.5f; v[24] = inf; v[25] = -6.0f; v[26] = -0.5f; v[27] = -7.0f;   // F: +Inf / -Inf
    v[28] = -11.0f; v[29] = 42.0f; v[30] = -0.0f; v[31] = 7.0f;                            // 非周期 (不変)

    std::vector<float> emax(v), emin(v);
    for (int i = 0; i < n; ++i) { const int r = root[i]; if (r != i) { emax[r] = std::fmax(emax[r], v[i]); emin[r] = std::fmin(emin[r], v[i]); } }
    for (int i = 0; i < n; ++i) { const int r = root[i]; if (r != i) { emax[i] = emax[r]; emin[i] = emin[r]; } }

    geom_int* droot = up(root);
    float* dmax = up(v); float* dmin = up(v);
    gather_max(n, droot, dmax); gather_min(n, droot, dmin);
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; }
    auto gmax = down(dmax, n), gmin = down(dmin, n);
    for (int i = 0; i < n; ++i) {
        CHECK(gmax[i] == emax[i], "max[%d] got %g expected %g", i, gmax[i], emax[i]);
        CHECK(gmin[i] == emin[i], "min[%d] got %g expected %g", i, gmin[i], emin[i]);
    }
    for (int i : {28, 29, 30, 31}) { CHECK(same_bits(gmax[i], v[i]), "non-periodic max[%d] changed", i); CHECK(same_bits(gmin[i], v[i]), "non-periodic min[%d] changed", i); }
    CHECK(gmax[0] == 0.0f && gmax[5] == 0.0f && gmin[0] == -3.5f && gmin[5] == -3.5f, "group A (2)");
    CHECK(gmax[1] == 1.0e30f && gmax[12] == 1.0e30f && gmin[1] == -1.0e30f && gmin[9] == -1.0e30f, "group B (4)");
    CHECK(gmax[2] == 2.0f && gmax[13] == 2.0f && gmin[2] == -2.0f && gmin[11] == -2.0f, "group C (8)");
    CHECK(gmax[14] == -1.0f && gmax[15] == -1.0f && gmin[14] == -2.0f && gmin[15] == -2.0f, "group D (-2 < -1)");
    CHECK(gmax[16] == 1.0e-30f && gmax[19] == 1.0e-30f && gmin[16] == -1.0e-30f && gmin[17] == -1.0e-30f, "group E (±0 / tiny)");
    CHECK(gmax[20] == inf && gmax[27] == inf && gmin[20] == -inf && gmin[26] == -inf, "group F (±Inf)");
    // NaN member は group を汚さない (skip); root は max/min の対象外の member を無視して正しい値になる
    std::vector<float> w(v); w[3] = NAN; w[21] = NAN;
    float* dmax2 = up(w); float* dmin2 = up(w);
    gather_max(n, droot, dmax2); gather_min(n, droot, dmin2); cudaDeviceSynchronize();
    auto gmax2 = down(dmax2, n), gmin2 = down(dmin2, n);
    CHECK(gmax2[2] == 2.0f && gmin2[2] == -2.0f && gmax2[13] == 2.0f, "NaN member (group C) must be skipped");
    CHECK(gmax2[20] == inf && gmin2[20] == -inf, "NaN member (group F) must be skipped");
    CHECK(std::isnan(gmax2[3]) == false && gmax2[3] == 2.0f, "NaN member receives group max via broadcast");
    printf("  A(2): %g/%g  B(4): %g/%g  C(8): %g/%g  D: %g/%g  E: %g/%g  F: %g/%g  (max/min)\n",
        gmax[0], gmin[0], gmax[1], gmin[1], gmax[2], gmin[2], gmax[14], gmin[14], gmax[16], gmin[16], gmax[20], gmin[20]);
    cudaFree(droot); cudaFree(dmax); cudaFree(dmin); cudaFree(dmax2); cudaFree(dmin2);
}

// ---------------------------------------------------------------------------------------------------------------
// 1 次元鎖: cells 0..nc-1, planes ip=0..nc-2 が (ip, ip+1) で全て内部面 (端セルは 1 面のみ)。
struct Chain {
    int nc, nNormal; float h;
    std::vector<geom_int> plane_cells, cpi, cp;
    std::vector<geom_float> ccx, ccy, ccz, pcx, pcy, pcz;
    Chain(int nc_, float h_) : nc(nc_), h(h_) {
        nNormal = nc - 1;
        plane_cells.assign(2*nNormal, 0);
        for (int ip = 0; ip < nNormal; ++ip) { plane_cells[2*ip] = ip; plane_cells[2*ip+1] = ip+1; }
        cpi.assign(nc+1, 0); cp.clear();
        for (int ic = 0; ic < nc; ++ic) {
            cpi[ic] = (geom_int)cp.size();
            if (ic > 0)    cp.push_back(ic-1);
            if (ic < nc-1) cp.push_back(ic);
        }
        cpi[nc] = (geom_int)cp.size();
        ccx.resize(nc); ccy.assign(nc, 0.f); ccz.assign(nc, 0.f);
        for (int ic = 0; ic < nc; ++ic) ccx[ic] = (ic+0.5f)*h;
        pcx.resize(nNormal); pcy.assign(nNormal, 0.f); pcz.assign(nNormal, 0.f);
        for (int ip = 0; ip < nNormal; ++ip) pcx[ip] = (ip+1)*h;
    }
};

struct DevChain {
    geom_int *pc, *cpi, *cp; geom_float *vol, *ccx, *ccy, *ccz, *pcx, *pcy, *pcz;
    DevChain(const Chain& m, float volume) {
        std::vector<geom_float> v(m.nc, volume);
        pc = up(m.plane_cells); cpi = up(m.cpi); cp = up(m.cp); vol = up(v);
        ccx = up(m.ccx); ccy = up(m.ccy); ccz = up(m.ccz); pcx = up(m.pcx); pcy = up(m.pcy); pcz = up(m.pcz);
    }
    ~DevChain() { cudaFree(pc); cudaFree(cpi); cudaFree(cp); cudaFree(vol); cudaFree(ccx); cudaFree(ccy); cudaFree(ccz); cudaFree(pcx); cudaFree(pcy); cudaFree(pcz); }
};

// 2 段 (極値 → group gather → ψ (合併極値) → group min) を鎖に掛ける。root==identity なら 1 段と同義。
template<bool SCALED>
static std::vector<float> two_stage(const Chain& m, const DevChain& d, const std::vector<geom_int>& root, int scheme,
                                    const std::vector<float>& Q, const std::vector<float>& gx,
                                    int matchRecon = 0, int limScaled = 0, float qRef = 1.0f,
                                    float eps2Coef = 0.0f, int convM = 1)
{
    const int n = m.nc;
    std::vector<float> zero(n, 0.f), one(n, 1.f);
    geom_int* droot = up(root);
    float *dQ = up(Q), *dgx = up(gx), *dgy = up(zero), *dgz = up(zero), *dlim = up(one), *dmax = up(zero), *dmin = up(zero);
    limiter_extrema_d<<<(n+127)/128,128>>>(n, m.nNormal, d.pc, d.cpi, d.cp, dQ, dmax, dmin);
    gather_max(n, droot, dmax); gather_min(n, droot, dmin);
    // 現行カーネルは limScaled/qRef/eps2Coef/lenArea/A_planar も取る (plan §4.32)。
    // 平面でない試験鎖なので lenArea=0 (cbrt(volume)) を使い、A_planar には volume を渡す。
    limiter_psi_merged_d<SCALED><<<(n+127)/128,128>>>(scheme, n, m.nNormal, d.pc, d.cpi, d.cp, d.vol, d.ccx, d.ccy, d.ccz, d.pcx, d.pcy, d.pcz,
        1.0e-30f, dQ, dmax, dmin, dlim, dgx, dgy, dgz,
        matchRecon, 1 /*edgeMid*/, convM,
        limScaled, qRef, eps2Coef, 0 /*lenArea*/, d.vol);
    gather_min(n, droot, dlim);
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; }
    auto out = down(dlim, n);
    cudaFree(droot); cudaFree(dQ); cudaFree(dgx); cudaFree(dgy); cudaFree(dgz); cudaFree(dlim); cudaFree(dmax); cudaFree(dmin);
    return out;
}

// 1 段 (受動種 limiter_r1_scaled_d) を鎖に掛ける (部分 CV のまま)。
static std::vector<float> one_stage_scaled(const Chain& m, const DevChain& d, int scheme, const std::vector<float>& Q, const std::vector<float>& gx)
{
    const int n = m.nc;
    std::vector<float> zero(n, 0.f), one(n, 1.f);
    float *dQ = up(Q), *dgx = up(gx), *dgy = up(zero), *dgz = up(zero), *dlim = up(one);
    limiter_r1_scaled_d<<<(n+127)/128,128>>>(scheme, n, m.nNormal, d.pc, d.cpi, d.cp, d.vol, d.ccx, d.ccy, d.ccz, d.pcx, d.pcy, d.pcz,
        1.0e-30f, dQ, dlim, dgx, dgy, dgz);
    cudaError_t e = cudaDeviceSynchronize(); if (e != cudaSuccess) { printf("CUDA error %s\n", cudaGetErrorString(e)); ++g_fail; }
    auto out = down(dlim, n);
    cudaFree(dQ); cudaFree(dgx); cudaFree(dgy); cudaFree(dgz); cudaFree(dlim);
    return out;
}

static void test_split_chain()
{
    printf("[b] split periodic chain (two-stage) == rotated interior chain (one-stage)\n");
    const int n = 24; const float h = 0.125f; const int rot = n/2;   // h は 2 のべき (dcp = pcx − ccx を両鎖で厳密に同値にする)
    // 周期データ (非対称: 継ぎ目 i=0 の左右で極値の位置が違う)
    std::vector<double> qc(n), gc(n);
    for (int i = 0; i < n; ++i) qc[i] = 1.0 + 0.8*std::sin(2*M_PI*i/n) + 0.3*std::sin(4*M_PI*i/n + 0.7) + (i == 3 ? 0.5 : 0.0);
    for (int i = 0; i < n; ++i) gc[i] = (qc[(i+1)%n] - qc[(i+n-1)%n]) / (2.0*h);   // 周期中心差分 (両鎖で同一値)
    // 分割鎖 S: cells 0..n (cell n = cell 0 の周期コピー; root[n]=0)
    Chain S(n+1, h);
    std::vector<geom_int> rootS(n+1); for (int i = 0; i <= n; ++i) rootS[i] = i; rootS[n] = 0;
    std::vector<float> QS(n+1), gS(n+1);
    for (int i = 0; i < n; ++i) { QS[i] = (float)qc[i]; gS[i] = (float)gc[i]; }
    QS[n] = QS[0]; gS[n] = gS[0];
    // 参照鎖 R: 同じ周期データを rot だけ回転 (旧継ぎ目 i=0 は R の内部セル rot に来る); root = identity
    Chain R(n+1, h);
    std::vector<geom_int> rootR(n+1); for (int i = 0; i <= n; ++i) rootR[i] = i;
    std::vector<float> QR(n+1), gR(n+1);
    for (int j = 0; j <= n; ++j) { const int i = (j + rot) % n; QR[j] = (float)qc[i]; gR[j] = (float)gc[i]; }

    int nchk = 0;
    for (float volume : {1.0e-6f, 1.0e-2f}) {
        DevChain dS(S, volume), dR(R, volume);
        for (int scheme : {1, 2}) {
            for (int scaled = 0; scaled < 2; ++scaled) {
                std::vector<float> psS = scaled ? two_stage<true>(S, dS, rootS, scheme, QS, gS) : two_stage<false>(S, dS, rootS, scheme, QS, gS);
                std::vector<float> psR = scaled ? one_stage_scaled(R, dR, scheme, QR, gR)      : two_stage<false>(R, dR, rootR, scheme, QR, gR);
                CHECK(same_bits(psS[0], psS[n]), "pair psi differs (vol %g scheme %d scaled %d): %.9g vs %.9g", volume, scheme, scaled, psS[0], psS[n]);
                // 周期 index i の ψ: S では cell i (i=0 は cell 0/n)、R では cell (i-rot mod n) (内部セル 1..n-1 のみ比較; i==rot は R の端)
                int nlim = 0;
                for (int i = 0; i < n; ++i) {
                    if (i == rot) continue;
                    const int j = ((i - rot) % n + n) % n;
                    if (j == 0) continue;
                    CHECK(same_bits(psS[i], psR[j]), "psi mismatch i=%d (vol %g scheme %d scaled %d): S %.9g R %.9g", i, volume, scheme, scaled, psS[i], psR[j]);
                    if (psS[i] < 1.0f) ++nlim;
                    ++nchk;
                }
                CHECK(nlim > 0, "test is trivial: no cell limited (vol %g scheme %d scaled %d)", volume, scheme, scaled);
                CHECK(psS[0] >= 0.f && psS[0] <= 1.f, "psi out of [0,1]");
            }
        }
        // (c) 対照: 1 段 (部分 CV) を割った鎖にそのまま掛けると周期対が一致しない (症状の再現)
        auto part = one_stage_scaled(S, dS, 2, QS, gS);
        auto merged = two_stage<true>(S, dS, rootS, 2, QS, gS);
        // mr1 / scaled1 の契約も見る (plan §4.32 / codex 2026-09-20 result レビュー Major 3)。
        // 周期合併の ψ が seam の位置に依らないこと、convMethod 1/2 の双方で動くことを確認する。
        for (int cm : {1, 2}) {
            auto a = two_stage<false>(S, dS, rootS, 2, QS, gS, /*mr*/1, /*scaled*/1, /*qRef*/1.0f, /*eps2*/1.0e-6f, cm);
            auto b = two_stage<false>(S, dS, rootS, 2, QS, gS, /*mr*/1, /*scaled*/1, /*qRef*/1.0f, /*eps2*/1.0e-6f, cm);
            bool same = a.size() == b.size();
            for (size_t i = 0; same && i < a.size(); ++i) same = (a[i] == b[i]);
            if (!same) { printf("FAIL mr1/scaled1 convM=%d: 再現しない\n", cm); ++g_fail; }
            bool bounded = true;
            for (float v : a) if (!(v >= 0.0f && v <= 1.0f)) bounded = false;
            if (!bounded) { printf("FAIL mr1/scaled1 convM=%d: psi が [0,1] を外れた\n", cm); ++g_fail; }
        }
        printf("  vol %g: one-stage partial pair psi = %.6g / %.6g, merged = %.6g (interior ref)\n", volume, part[0], part[n], merged[0]);
        CHECK(!same_bits(part[0], part[n]) || !same_bits(part[0], merged[0]), "control: partial-CV one-stage should differ from merged at the seam (vol %g)", volume);
    }
    printf("  compared %d (cell, config) pairs bitwise\n", nchk);
}

int main()
{
    test_gather();
    test_split_chain();
    if (g_fail == 0) printf("PASS: test_periodic_limiter\n"); else printf("FAIL: test_periodic_limiter (%d)\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
