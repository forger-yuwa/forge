// =============================================================================
// test_cond_latent_ab.cu — V7(f) の 0 step A/B: 潜熱の供給関数だけを替えて float 表と利用先の差を測る
//   plans/active/thermophysics-solver-owned-species-db.md §5.1 #10・§6 V7(f)、判断 notes/reviews/2026-10-06-v7f-latent-regression-diagnose.md
//   A = #10 以前の h2o_latent (test_cond_latent_pair.cu の移植と同じ式)、B = 気液ペア方式 (h2o_latent_pair)。
//   B の気相・液相の MW は A と同じ 0.0180153 に固定 (#13-3 の MW 変更を混ぜない; 表を作った後の MW 比掛けでは代用しない)。
//   同じ cond_tables_fill・同じ格子 (cond_tables_grid H2O) で潜熱表を作り、
//   (1) 係数 c0..c3 の成分別ビット不一致数 (区間全体・200 K 以上・200 K をまたぐ 199.90–200.15 K)、
//   (2) device の cond_tab_eval で L・L' (200–400 K を 1e-3 K 刻み 200001 点 + 200 K と float 両隣) の不一致点数・最大差、
//   (3) 利用先: 二相 EOS の熱容量 cp2 = cp_allvap − g L' と、Kantrowitz 補正の L/(R T) の差 (固定状態の組)、
//   (4) 正の対照: 150 K で double の ΔL ≈ −465.920 J/kg、
//   を記録する (判定ではなく記録; 第 1 仮説「float 表の L' に差が残る」の採否を出す)。
// ビルド (AWS; test_cond_latent_pair.cu と同じ):
//   cmake -DIN=data/species/forge_species_v1.yaml -DOUT=<GEN>/forge_species_data.hpp -P cmake/embed_species_data.cmake
//   nvcc -std=c++17 --expt-relaxed-constexpr -arch=sm_86 -I . -I <GEN> -o <OUT>/test_cond_latent_ab tests/unit/test_cond_latent_ab.cu
//        -x none input/speciesDB.cpp input/speciesTransportDB.cpp -lyaml-cpp
// =============================================================================
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include "tests/unit/cond_latent_test_helper.cuh"
#include "cuda_forge/condensationTables_d.cuh"

static const double MW_OLD = 0.0180153;

// #10 以前のソルバの h2o_latent (test_cond_latent_pair.cu の legacy_h2o_latent_v0 と同じ)
static double legacy_h(const double* a, double T)
{
    const double Rw = 8.314462618/0.0180153;
    const double hRT = -a[0]/(T*T) + a[1]*log(T)/T + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0 + a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[7]/T;
    return hRT*Rw*T;
}
static double legacy_L(double T)
{
    const double ag[8] = {-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03,
                          -7.342557370e-06, 4.955043490e-09,-1.336933246e-12,-3.303974310e+04};
    const double al[8] = { 1.326371304e+09,-2.448295388e+07, 1.879428776e+05,-7.678995050e+02,
                           1.761556813e+00,-2.151167128e-03, 1.092570813e-06, 1.101760476e+08};
    const double Tf = 273.15;
    double Tg = (T > 45.0) ? T : 45.0;
    if (Tg > 1000.0) Tg = 1000.0;
    const double hv = legacy_h(ag, Tg);
    double hl;
    if (Tg >= Tf) hl = legacy_h(al, (Tg < 373.15) ? Tg : 373.15);
    else hl = legacy_h(al, Tf) - (legacy_h(al, Tf + 0.5) - legacy_h(al, Tf - 0.5 + 1.0e-9))*(Tf - Tg);
    double L = hv - hl;
    if (L < 1.5e6) L = 1.5e6;
    if (L > 3.5e6) L = 3.5e6;
    return L;
}

__global__ void eval_dev(CondTableF ta, CondTableF tb, const float* T, float* La, float* Lb, float* da, float* db, int n)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    La[i] = cond_tab_eval(ta, T[i], &da[i]);
    Lb[i] = cond_tab_eval(tb, T[i], &db[i]);
}

static unsigned fbits(float x) { unsigned u; memcpy(&u, &x, 4); return u; }

int main()
{
    // B: 気液ペア (Tref 298.15 = 生産の thermoHrefTemp) を MW_OLD で組み直す
    ResolvedSpeciesDB db = speciesDB_resolve(std::vector<std::string>{"H2O"}, "");
    speciesDB_attachCondensed(db, "H2O(L)", "H2O", true);
    ResolvedCondensed c = db.condensed;
    printf("gas MW (DB) %.10g -> fixed %.10g\n", c.gas.MW, MW_OLD);
    c.gas.MW = MW_OLD; c.gas.invMW = 1.0/MW_OLD;
    const CondLatentPair pB = cond_latent_pair_make(c.gas, c.coeffs, c.Tlo, c.Thi, 298.15);
    auto fA = [&](double T){ return legacy_L(T); };
    auto fB = [&](double T){ return h2o_latent_pair(pB, T); };

    printf("== (4) positive control: double ΔL = B − A ==\n");
    for (double T : {120.0, 150.0, 199.9, 200.0, 200.1, 207.5, 228.0, 250.0, 300.0, 373.15})
        printf("  T %7.2f K: L_A %.9e  L_B %.9e  ΔL %+.6e J/kg (rel %+.2e)\n", T, fA(T), fB(T), fB(T) - fA(T), (fB(T) - fA(T))/fA(T));

    double T0, h; int n;
    cond_tables_grid(COND_MODEL_H2O, &T0, &h, &n);
    std::vector<float4> cA, cB;
    cond_tables_fill(cA, T0, h, n, fA);
    cond_tables_fill(cB, T0, h, n, fB);
    printf("== (1) table coefficients (grid T0 %.2f h %.2f n %d): bitwise mismatches c0 c1 c2 c3 ==\n", T0, h, n);
    auto cnt = [&](double Tlo, double Thi, const char* tag) {
        int m[4] = {0, 0, 0, 0}, tot = 0;
        for (int i = 0; i < n; ++i) {
            const double Ta = T0 + h*i, Tb = Ta + h;
            if (Tb <= Tlo || Ta >= Thi) continue;
            ++tot;
            m[0] += fbits(cA[i].x) != fbits(cB[i].x); m[1] += fbits(cA[i].y) != fbits(cB[i].y);
            m[2] += fbits(cA[i].z) != fbits(cB[i].z); m[3] += fbits(cA[i].w) != fbits(cB[i].w);
        }
        printf("  %-28s intervals %5d : %d %d %d %d\n", tag, tot, m[0], m[1], m[2], m[3]);
    };
    cnt(120.0, 200.0, "120-200 K");
    cnt(199.90, 200.15, "199.90-200.15 K (crossing)");
    cnt(200.0, 400.0, "200-400 K");
    cnt(400.0, 1200.0, "400-1200 K");

    printf("== (2) device cond_tab_eval: L and L' (200-400 K every 1e-3 K + 200 K neighbours) ==\n");
    std::vector<float> T;
    for (int k = 0; k <= 200000; ++k) T.push_back((float)(200.0 + 1e-3*k));
    T.push_back(nextafterf(200.0f, 0.0f)); T.push_back(200.0f); T.push_back(nextafterf(200.0f, 1000.0f));
    const int N = (int)T.size();
    float4 *dA, *dB; float *dT, *La, *Lb, *Da, *Db;
    cudaMalloc(&dA, n*sizeof(float4)); cudaMalloc(&dB, n*sizeof(float4));
    cudaMemcpy(dA, cA.data(), n*sizeof(float4), cudaMemcpyHostToDevice); cudaMemcpy(dB, cB.data(), n*sizeof(float4), cudaMemcpyHostToDevice);
    cudaMalloc(&dT, N*sizeof(float)); cudaMemcpy(dT, T.data(), N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMalloc(&La, N*sizeof(float)); cudaMalloc(&Lb, N*sizeof(float)); cudaMalloc(&Da, N*sizeof(float)); cudaMalloc(&Db, N*sizeof(float));
    CondTableF ta, tb;
    ta.c = dA; ta.T0 = (float)T0; ta.invH = (float)(1.0/h); ta.n = n;
    tb = ta; tb.c = dB;
    eval_dev<<<(N + 255)/256, 256>>>(ta, tb, dT, La, Lb, Da, Db, N);
    cudaDeviceSynchronize();
    std::vector<float> hLa(N), hLb(N), hDa(N), hDb(N);
    cudaMemcpy(hLa.data(), La, N*4, cudaMemcpyDeviceToHost); cudaMemcpy(hLb.data(), Lb, N*4, cudaMemcpyDeviceToHost);
    cudaMemcpy(hDa.data(), Da, N*4, cudaMemcpyDeviceToHost); cudaMemcpy(hDb.data(), Db, N*4, cudaMemcpyDeviceToHost);
    int nL = 0, nD = 0; double mL = 0, mD = 0, mDrel = 0;
    for (int i = 0; i < N; ++i) {
        if (fbits(hLa[i]) != fbits(hLb[i])) { ++nL; mL = fmax(mL, fabs((double)hLa[i] - hLb[i])); }
        if (fbits(hDa[i]) != fbits(hDb[i])) { ++nD; mD = fmax(mD, fabs((double)hDa[i] - hDb[i])); mDrel = fmax(mDrel, fabs((double)hDa[i] - hDb[i])/fabs((double)hDa[i])); }
    }
    printf("  points %d: L mismatches %d (max |ΔL| %.3e J/kg), L' mismatches %d (max |ΔL'| %.3e J/(kg K), max rel %.2e)\n", N, nL, mL, nD, mD, mDrel);

    printf("== (3) downstream with the table values (fixed states) ==\n");
    // 二相熱容量の L' の項 g·L' と Kantrowitz の L/(R_v T) (R_v = Ru/MW_OLD) を、L・L' の旧新で比べる
    const double Rv = 8.314462618/MW_OLD;
    double mcp = 0, mK = 0;
    CondTableF hA, hB;
    hA.c = cA.data(); hA.T0 = (float)T0; hA.invH = (float)(1.0/h); hA.n = n;
    hB = hA; hB.c = cB.data();
    for (double Tq : {207.5, 210.0, 228.0, 250.0, 280.0, 300.0}) for (double g : {1e-4, 1e-3, 1.1e-2}) {
        float dA_, dB_;
        const float lA = cond_tab_eval(hA, (float)Tq, &dA_);
        const float lB = cond_tab_eval(hB, (float)Tq, &dB_);
        mcp = fmax(mcp, fabs(g*((double)dB_ - dA_)));
        mK  = fmax(mK, fabs(((double)lB - lA)/(Rv*Tq)));
    }
    printf("  max |Δ(g L')| over T {207.5..300} x g {1e-4,1e-3,1.1e-2}: %.3e J/(kg K)   max |Δ(L/(R_v T))|: %.3e\n", mcp, mK);
    printf("DONE\n");
    return 0;
}
