// =============================================================================
// test_thermo_intervals_bitexact.cu — 区間可変 (plan thermophysics-solver-owned-species-db §5.1 #13-1 G1-b) の前後比較用ダンプ。
//   内蔵 7 種 (N2 O2 CO2 H2O Ar He AIR) と 7 種の等分混合について、T = {200, 999.99994, 1000, 1000.00006, 6000, 7000} K と
//   datum (thermoHrefTemp 0 / 298.15) の全組で、device と host の両方で
//     double: cp, h, s°, cph (融合), T_from_e (warm 0.9T / 1.1T), T_from_h, 気相潜熱の h_v (cond_gas_h_mass), 研磨 cph
//     float : cp, h, cph (融合), T_from_e_f, ハイブリッド反転 (T・cp・h・T_f)
//     輸送 : 旧 kinetic 経路の μ・λ (thermo_mu_species / thermo_lambda_species; λ が c_p を読む)
//   を 16 進 (%a) で 1 行ずつ出す。**同じソースを本変更前 (2 区間固定) のヘッダと本変更後のヘッダでビルド**し、出力を diff して 0 ulp を確かめる
//   (THERMO_MAX_INTERVALS の有無で構造体の組み方だけを切り替える; 評価関数は両版で同名)。
//
// ビルド (共通データの埋め込みヘッダを先に生成; SRC は本変更前 / 後のどちらかの solver_density_cuda):
//   cmake -DIN=$SRC/data/species/forge_species_v1.yaml -DOUT=$GEN/forge_species_data.hpp -P $SRC/cmake/embed_species_data.cmake
//   nvcc -O3 -std=c++17 -arch=sm_86 --expt-relaxed-constexpr -I $SRC -I $GEN -o thermo_bitexact \
//       solver_density_cuda/tests/unit/test_thermo_intervals_bitexact.cu $SRC/input/speciesDB.cpp $SRC/input/speciesTransportDB.cpp -lyaml-cpp
//   ./thermo_bitexact > out.txt   (新旧を diff; 差が無ければ 0 ulp)
// =============================================================================
#include "input/speciesDB.hpp"
#include "cuda_forge/thermo_d.cuh"
#include "cuda_forge/condensationProperties_d.cuh"

#include <cstdio>
#include <string>
#include <vector>
#include <cuda_runtime.h>

static const char* kNames[] = {"N2", "O2", "CO2", "H2O", "Ar", "He", "AIR"};
static const int kNS = 7;
static const double kT[] = {200.0, 999.99994, 1000.0, 1000.00006, 6000.0, 7000.0};
static const int kNT = 6;
#define NQ 26   // 1 点あたりの出力数

static void applyDatum(SpeciesThermo& s, double Tref)
{
    const double h_ref = thermo_h_molar(s, Tref);
    const double da7 = -h_ref/THERMO_RU;
#ifdef THERMO_MAX_INTERVALS
    thermo_add_a7(s, da7);
#else
    s.low[7] += da7; s.high[7] += da7;
#endif
    s.h_datum = h_ref;
}

static SpeciesThermoF toF(const SpeciesThermo& s)
{
#ifdef THERMO_MAX_INTERVALS
    return thermo_to_float(s);
#else
    SpeciesThermoF f;
    f.MW = (float)s.MW; f.invMW = (float)(1.0/s.MW); f.R = (float)(THERMO_RU/s.MW);
    f.sigma_LJ = (float)s.sigma_LJ; f.eps_kB = (float)s.eps_kB;
    f.Tlo = (float)s.Tlo; f.Tmid = (float)s.Tmid; f.Thi = (float)s.Thi;
    for (int k = 0; k < 9; ++k) { f.low[k] = (float)s.low[k]; f.high[k] = (float)s.high[k]; }
    return f;
#endif
}

// case: 0..kNS-1 は単成分 (sp+case, n=1)、kNS は 7 種の等分混合 (n=kNS)。
__host__ __device__ inline void evalPoint(const SpeciesThermo* sp, const SpeciesThermoF* spf, int n, const double* Y, const float* Yf,
                                          double T, double* q)
{
    int j = 0;
    const float Tf = (float)T;
    // double
    q[j++] = thermo_cp_mix(sp, n, Y, T);
    q[j++] = thermo_h_mix(sp, n, Y, T);
    q[j++] = thermo_s0_mix(sp, n, Y, T);
    { double cp, h; thermo_cph_mix(sp, n, Y, T, &cp, &h); q[j++] = cp; q[j++] = h; }
    { double cp, h; thermo_cph_mix_polish(sp, n, Y, T, &cp, &h); q[j++] = cp; q[j++] = h; }
    const double R = thermo_R_mix(sp, n, Y);
    const double e = thermo_h_mix(sp, n, Y, T) - R*T;
    q[j++] = thermo_T_from_e(sp, n, Y, e, 0.9*T, 50.0, 20000.0);
    q[j++] = thermo_T_from_e(sp, n, Y, e, 1.1*T, 50.0, 20000.0);
    q[j++] = thermo_T_from_h(sp, n, Y, thermo_h_mix(sp, n, Y, T), 0.9*T, 50.0, 20000.0);
    q[j++] = cond_gas_h_mass(sp[0], T);   // 潜熱の気相 h_v (種 DB と同じ評価; 単成分 sp[0])
    // float
    { float s = 0.0f; for (int i = 0; i < n; ++i) s += Yf[i]*thermo_cp_mass_f(spf[i], Tf); q[j++] = s; }
    q[j++] = thermo_h_mix_f(spf, n, Yf, Tf);
    { float cp, h; thermo_cph_mix_f(spf, n, Yf, Tf, &cp, &h); q[j++] = cp; q[j++] = h; }
    const float Rf = thermo_R_mix_f(spf, n, Yf);
    { float cp, h; thermo_cph_mix_f(spf, n, Yf, Tf, &cp, &h); const float ef = h - Rf*Tf;
      q[j++] = thermo_T_from_e_f(spf, n, Yf, ef, 0.9f*Tf, 50.0f, 20000.0f, nullptr, 12); }
    { double cpH, hH, TfH;
      const double Th = thermo_T_from_e_hybrid(sp, spf, n, Y, Yf, e, 0.9*T, 50.0, 20000.0, &cpH, &hH, &TfH, 8);
      q[j++] = Th; q[j++] = cpH; q[j++] = hH; q[j++] = TfH; }
    // 輸送 (旧 kinetic 経路; λ = μ(c_p + 1.25 R) が c_p を読む)
    q[j++] = thermo_mu_species(sp[0], T);
    q[j++] = thermo_lambda_species(sp[0], T);
    { double cp, h; thermo_cph_molar(sp[0], T, &cp, &h); q[j++] = cp; q[j++] = h; }
    q[j++] = thermo_s0_mass(sp[0], T);
    q[j++] = thermo_cp_molar(sp[0], T);
}

__global__ void kEval(const SpeciesThermo* sp, const SpeciesThermoF* spf, const double* Yall, const float* Yfall, double* out)
{
    const int c = blockIdx.x;          // case 0..kNS
    const int t = threadIdx.x;         // T index
    if (t >= kNT) return;
    const double Ts[] = {200.0, 999.99994, 1000.0, 1000.00006, 6000.0, 7000.0};
    const int n = (c < kNS) ? 1 : kNS;
    const int off = (c < kNS) ? c : 0;
    evalPoint(sp + off, spf + off, n, Yall + c*kNS, Yfall + c*kNS, Ts[t], out + (c*kNT + t)*NQ);
}

int main()
{
    const auto b = speciesDB_builtin();
    for (int datum = 0; datum < 2; ++datum) {
        std::vector<SpeciesThermo> sp(kNS);
        for (int i = 0; i < kNS; ++i) { sp[i] = b.at(kNames[i]); sp[i].invMW = 1.0/sp[i].MW; if (datum) applyDatum(sp[i], 298.15); }
        std::vector<SpeciesThermoF> spf(kNS);
        for (int i = 0; i < kNS; ++i) spf[i] = toF(sp[i]);
        const int nC = kNS + 1;
        std::vector<double> Y(nC*kNS, 0.0);
        std::vector<float> Yf(nC*kNS, 0.0f);
        for (int c = 0; c < kNS; ++c) { Y[c*kNS] = 1.0; Yf[c*kNS] = 1.0f; }
        for (int i = 0; i < kNS; ++i) { Y[kNS*kNS + i] = 1.0/kNS; Yf[kNS*kNS + i] = (float)(1.0/kNS); }
        // host
        std::vector<double> hq(nC*kNT*NQ);
        for (int c = 0; c < nC; ++c)
            for (int t = 0; t < kNT; ++t) {
                const int n = (c < kNS) ? 1 : kNS, off = (c < kNS) ? c : 0;
                evalPoint(sp.data() + off, spf.data() + off, n, Y.data() + c*kNS, Yf.data() + c*kNS, kT[t], hq.data() + (c*kNT + t)*NQ);
            }
        // device
        SpeciesThermo* dsp; SpeciesThermoF* dspf; double* dY; float* dYf; double* dq;
        cudaMalloc(&dsp, kNS*sizeof(SpeciesThermo)); cudaMalloc(&dspf, kNS*sizeof(SpeciesThermoF));
        cudaMalloc(&dY, Y.size()*sizeof(double)); cudaMalloc(&dYf, Yf.size()*sizeof(float)); cudaMalloc(&dq, hq.size()*sizeof(double));
        cudaMemcpy(dsp, sp.data(), kNS*sizeof(SpeciesThermo), cudaMemcpyHostToDevice);
        cudaMemcpy(dspf, spf.data(), kNS*sizeof(SpeciesThermoF), cudaMemcpyHostToDevice);
        cudaMemcpy(dY, Y.data(), Y.size()*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(dYf, Yf.data(), Yf.size()*sizeof(float), cudaMemcpyHostToDevice);
        kEval<<<nC, 32>>>(dsp, dspf, dY, dYf, dq);
        std::vector<double> dqh(hq.size());
        if (cudaMemcpy(dqh.data(), dq, hq.size()*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess || cudaGetLastError() != cudaSuccess) {
            std::fprintf(stderr, "CUDA error\n"); return 1;
        }
        cudaFree(dsp); cudaFree(dspf); cudaFree(dY); cudaFree(dYf); cudaFree(dq);
        for (int c = 0; c < nC; ++c)
            for (int t = 0; t < kNT; ++t) {
                const char* nm = (c < kNS) ? kNames[c] : "MIX7";
                std::printf("datum=%d %s T=%a host", datum, nm, kT[t]);
                for (int k = 0; k < NQ; ++k) std::printf(" %a", hq[(c*kNT + t)*NQ + k]);
                std::printf("\ndatum=%d %s T=%a dev ", datum, nm, kT[t]);
                for (int k = 0; k < NQ; ++k) std::printf(" %a", dqh[(c*kNT + t)*NQ + k]);
                std::printf("\n");
            }
    }
    return 0;
}
