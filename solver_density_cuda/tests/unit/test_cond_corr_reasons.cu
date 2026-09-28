// test_cond_corr_reasons.cu — 凝縮の理由別の補正量監視 (plans/active/condensation-two-phase-transport.md §4.3, §5.1 #2) の単体試験。
//   build: nvcc --expt-relaxed-constexpr -I. -o test_cond_corr_reasons tests/unit/test_cond_corr_reasons.cu
//   実現可能性クランプ (double / float 実体) と Q1/Q2 射影を、計上先あり (reasons) となし (nullptr) で同じ入力に掛け、
//   (1) 書き込み値 (rog, Q0..Q2, 診断) がビット一致 (計上は振る舞いを変えない)、
//   (2) 理由別の計上値が手計算と一致 (蒸気上限違反・負値 floor・射影・液滴消滅・制限前の最小蒸気分率・体積の重み) を確かめる。
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/thermo_d.cuh"
namespace {
#include "cuda_forge/condensationSourceKernels_d.cuh"
}
#include "cuda_forge/condensationRealizability_d.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } else { printf("  ok: " __VA_ARGS__); printf("\n"); } } while (0)
template <class T> static T* up(const std::vector<T>& v) { T* d = nullptr; cudaMalloc((void**)&d, v.size()*sizeof(T)); cudaMemcpy(d, v.data(), v.size()*sizeof(T), cudaMemcpyHostToDevice); return d; }
template <class T> static std::vector<T> down(const T* d, size_t n) { std::vector<T> v(n); cudaMemcpy(v.data(), d, n*sizeof(T), cudaMemcpyDeviceToHost); return v; }

struct Cell { double Yw, g, q0, x, y; };   // 単分散半径 r0 を基準に Q1 = Q0 r0 x, Q2 = Q0 r0² y

static bool near(double a, double b, double rt) { return std::fabs(a - b) <= rt*std::max(std::fabs(a), std::fabs(b)) + 1e-300; }

int main()
{
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o);
    CondTablesHost ht; cond_tables_build_host(cp, ht); const CondTablesF tb = cond_tables_upload(ht);
    const double T = 250.0, ro = 1.0, rho_l = cond_rho_cond(cp, T);
    const double r0 = 5.0e-8;
    auto q0_of = [&](double g, double r) { return g/((4.0/3.0)*COND_PI*rho_l*r*r*r); };
    // 0: 正常 / 1: 蒸気上限違反 g > Y_w / 2: 負値 g < 0 / 3: 射影 (x = 1.01) / 4: 液滴消滅 (S ≤ 1, r30 = 1 nm < 2 r_min)
    const double gRm = 1.0e-7, rRm = 1.0e-9;
    std::vector<Cell> cells = {
        {0.01, 0.005, q0_of(0.005, r0), 1.0, 1.0},
        {0.01, 0.012, q0_of(0.012, r0), 1.0, 1.0},
        {0.01, -1.0e-4, 0.0, 0.0, 0.0},
        {0.01, 0.005, q0_of(0.005, r0), 1.01, 1.0},
        {gRm + 1.0e-5, gRm, q0_of(gRm, rRm), 1.0, 1.0},
    };
    const int n = (int)cells.size();
    std::vector<flow_float> hro(n, (flow_float)ro), hYw(n), hg(n), h0(n), h1(n), h2(n), hT(n, (flow_float)T), hP(n, 101325.0f), hV(n), z(n, 0.0f);
    for (int i = 0; i < n; ++i) {
        const Cell& c = cells[i];
        const double r = (i == 4) ? rRm : r0;
        hYw[i] = (flow_float)(ro*c.Yw); hg[i] = (flow_float)(ro*c.g); h0[i] = (flow_float)(ro*c.q0);
        h1[i] = (flow_float)(ro*c.q0*r*c.x); h2[i] = (flow_float)(ro*c.q0*r*r*c.y); hV[i] = (flow_float)(1.0 + i);
    }
    const flow_float* dro = up(hro); flow_float* dYw = up(hYw); flow_float* dT = up(hT); flow_float* dP = up(hP); flow_float* dV = up(hV);
    const double g_rm = 5.0e-7, rmin = 1.0e-9;

    for (int useF = 0; useF < 2; ++useF) {
        printf("== realizability clamp (%s) ==\n", useF ? "float" : "double");
        std::vector<std::vector<flow_float>> outs[2];
        std::vector<double> rs;
        for (int withR = 0; withR < 2; ++withR) {
            flow_float *dg = up(hg), *d0 = up(h0), *d1 = up(h1), *d2 = up(h2), *cG = up(z), *cQ = up(z);
            std::vector<double> hb(8, 0.0); double* db = up(hb);
            std::vector<double> hr(COND_REASON_N, 0.0); hr[COND_REASON_VMIN] = 1.0e300; double* dr = up(hr);
            int hv[2] = {0, 0}; int* dv = nullptr; cudaMalloc((void**)&dv, 2*sizeof(int)); cudaMemcpy(dv, hv, 2*sizeof(int), cudaMemcpyHostToDevice);
            if (useF)
                cond_realizability_clamp_f_d<<<1, 32>>>(n, (flow_float*)dro, dYw, dg, d0, d1, d2, 1, (float)cp.R, (float)rmin, (float)g_rm, 0.0f,
                                                       dT, dP, tb, cp, cG, cQ, dv, db, nullptr, dV, 1, withR ? dr : nullptr);
            else
                cond_realizability_clamp_d<<<1, 32>>>(n, (flow_float*)dro, dYw, dg, d0, d1, d2, 1, COND_MODEL_H2O, cp.R, rmin, g_rm,
                                                     dT, dP, o, cG, cQ, dv, db, nullptr, dV, 1, withR ? dr : nullptr);
            cudaDeviceSynchronize();
            outs[withR] = {down(dg, n), down(d0, n), down(d1, n), down(d2, n), down(cG, n), down(cQ, n)};
            std::vector<double> bb = down(db, 8);
            if (withR) rs = down(dr, COND_REASON_N);
            static std::vector<double> b0; if (!withR) b0 = bb; else CHECK(std::memcmp(b0.data(), bb.data(), 8*sizeof(double)) == 0, "clampBudget identical with/without reasons");
        }
        bool same = true;
        for (size_t k = 0; k < 6; ++k) same = same && std::memcmp(outs[0][k].data(), outs[1][k].data(), n*sizeof(flow_float)) == 0;
        CHECK(same, "written state (rog, Q0, Q1, Q2, condClampCorr, condClampCorrQ) bit-identical with/without reasons");
        const std::vector<flow_float>& g1 = outs[1][0];
        CHECK(g1[1] == hYw[1] && g1[2] == 0.0f && g1[4] == 0.0f, "clamp results: capped %.6g (Y_w %.6g), negative -> %.3g, removed -> %.3g",
              (double)g1[1], (double)hYw[1], (double)g1[2], (double)g1[4]);
        const double capExp = ((double)hg[1] - (double)hYw[1])*hV[1], negExp = -(double)hg[2]*hV[2], rmExp = (double)hg[4]*hV[4];
        CHECK(near(rs[COND_REASON_CAP_SUM], capExp, 1e-12) && rs[COND_REASON_CAP_N] == 1.0,
              "capViol %.9e (expected %.9e) nodes %.0f", rs[COND_REASON_CAP_SUM], capExp, rs[COND_REASON_CAP_N]);
        CHECK(near(rs[COND_REASON_NEG_SUM], negExp, 1e-12) && rs[COND_REASON_NEG_N] == 1.0,
              "negFloor %.9e (expected %.9e) nodes %.0f", rs[COND_REASON_NEG_SUM], negExp, rs[COND_REASON_NEG_N]);
        CHECK(near(rs[COND_REASON_RM_SUM], rmExp, 1e-12) && rs[COND_REASON_RM_N] == 1.0,
              "removal %.9e (expected %.9e) nodes %.0f", rs[COND_REASON_RM_SUM], rmExp, rs[COND_REASON_RM_N]);
        // 射影は cell 3 (x = 1.01) と cell 1 (上限で g を削ると r が縮み x > 1 になる) の 2 点。期待値は出力から (Q1, Q2 は負値 floor 済みの入力との差)
        double p1 = 0.0, p2 = 0.0; int pn = 0;
        for (int i = 0; i < n; ++i) {
            const double a1 = std::max((double)h1[i], 0.0), a2 = std::max((double)h2[i], 0.0);
            const double d1 = std::fabs((double)outs[1][2][i] - a1), d2 = std::fabs((double)outs[1][3][i] - a2);
            if (i == 4) continue;   // 液滴消滅は射影ではない (Q を 0 にする; 消滅として計上)
            if (d1 > 0.0 || d2 > 0.0) { p1 += d1*hV[i]; p2 += d2*hV[i]; ++pn; }
        }
        CHECK(rs[COND_REASON_P_N] == (double)pn && pn == 2 && near(rs[COND_REASON_PQ1_SUM], p1, 1e-12) && near(rs[COND_REASON_PQ2_SUM], p2, 1e-12),
              "projection Q1 %.6e (expected %.6e) Q2 %.6e (expected %.6e) nodes %.0f (expected %d: x=1.01 and the capped cell)",
              rs[COND_REASON_PQ1_SUM], p1, rs[COND_REASON_PQ2_SUM], p2, rs[COND_REASON_P_N], pn);
        const double vmExp = ((double)hYw[1] - (double)hg[1])/ro;
        CHECK(rs[COND_REASON_VMIN] == vmExp, "min vapor fraction before limiting %.9e (expected %.9e)", rs[COND_REASON_VMIN], vmExp);
    }
    // Q1/Q2 だけの射影 (dual-time の step 末尾)
    {
        printf("== project only ==\n");
        std::vector<flow_float> r1[2], r2[2]; std::vector<double> rs;
        for (int withR = 0; withR < 2; ++withR) {
            flow_float *dg = up(hg), *d0 = up(h0), *d1 = up(h1), *d2 = up(h2), *cQ = up(z);
            std::vector<double> hr(COND_REASON_N, 0.0); double* dr = up(hr);
            cond_realizability_project_only_d<<<1, 32>>>(n, (flow_float*)dro, dg, d0, d1, d2, COND_MODEL_H2O, dT, o, 0, tb, cQ, nullptr, nullptr, nullptr, dV,
                                                         withR ? dr : nullptr);
            cudaDeviceSynchronize();
            r1[withR] = down(d1, n); r2[withR] = down(d2, n);
            if (withR) rs = down(dr, COND_REASON_N);
        }
        CHECK(std::memcmp(r1[0].data(), r1[1].data(), n*sizeof(flow_float)) == 0 && std::memcmp(r2[0].data(), r2[1].data(), n*sizeof(flow_float)) == 0,
              "project-only: Q1/Q2 bit-identical with/without reasons");
        const double p1 = std::fabs((double)r1[1][3] - (double)h1[3])*hV[3];
        CHECK(rs[COND_REASON_P_N] >= 1.0 && rs[COND_REASON_PQ1_SUM] >= p1 && p1 > 0.0, "project-only: projection nodes %.0f Q1 %.6e (cell 3 alone %.6e)",
              rs[COND_REASON_P_N], rs[COND_REASON_PQ1_SUM], p1);
    }
    printf(g_fail ? "FAILED (%d)\n" : "ALL PASS (%d failures)\n", g_fail);
    return g_fail ? 1 : 0;
}
