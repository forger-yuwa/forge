// =============================================================================
// test_dmix_complement.cu — 混合平均拡散係数の分子の補数形 (plans/active/condensation-two-phase-transport.md §5.1 #3b) の
//   係数ハーネス (CFD 0 step)。判定は driver の test_dmix_complement.py が行い、本体は評価値を書き出すだけ。
//   A = 現行分子 1−X_i (#3b 以前の thermo_Dmix_species_f を下の dmix_A に移植)
//   B = 補数形 Σ_{j≠i} X_j (ソルバの thermo_Dmix_species_f そのもの)
//   分母・二元係数 (thermo_Dbinary_f)・ガード (n==1 → 0, denom<1e-30 → 自己拡散)・コンパイル条件は A/B で共通。
//   組成は species_Dmix_probe_d (speciesTransport_d.cu) と同じ: Y = max(ρY/ρ, 0) → 和で正規化 → 気相組成 (液 gl = ρg/ρ/ΣY を
//   凝縮種から引く) → thermo_X_from_Y_f。
//
// 入力 (テキスト):
//   1 行目: nSp nStates
//   続く nSp 行: MW[kg/mol] sigma_LJ[Å] eps_kB[K]   (double; ソルバの float ミラー thermo_d.cu と同じく float へ丸める)
//   続く nStates 行: n idx_0..idx_{n-1} iw T P ro rog roY_0..roY_{n-1}   (iw < 0 は液なし; T,P,ro,rog,roY は float へ丸める;
//                    n = 1 の行は n==1 ガード (A/B とも 0) の確認)
// 出力 (テキスト, 1 状態 1 行, %.9g = float の往復一致): X_0..X_{n-1} DA_0..DA_{n-1} DB_0..DB_{n-1}
//
// ビルド/実行 (driver が行う):
//   nvcc -std=c++17 -O2 -arch=sm_86 --expt-relaxed-constexpr -I solver_density_cuda -o test_dmix_complement
//       solver_density_cuda/tests/unit/test_dmix_complement.cu && ./test_dmix_complement in.txt out.txt
// =============================================================================
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "cuda_forge/thermo_d.cuh"
#include "cuda_forge/gasPhaseComposition_d.cuh"

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA %s: %s\n", #x, cudaGetErrorString(e_)); exit(2); } } while (0)

// A: #3b 以前の thermo_Dmix_species_f (分子 1−X_i)。比較の基準としてだけ移植。
__device__ float dmix_A(const SpeciesThermoF* sp, int n, const float* X, int i, float T, float P)
{
    if (n == 1) return 0.0f;
    float denom = 0.0f;
    for (int j=0;j<n;j++) {
        if (j==i) continue;
        const float Dij = thermo_Dbinary_f(sp[i], sp[j], T, P);
        denom += X[j]/(Dij > 1.0e-30f ? Dij : 1.0e-30f);
    }
    if (denom < 1.0e-30f) return thermo_Dbinary_f(sp[i], sp[i], T, P);
    return (1.0f - X[i])/denom;
}

struct State {
    int n, idx[THERMO_MAX_SPECIES], iw;
    float T, P, ro, rog, roY[THERMO_MAX_SPECIES];
};

__global__ void eval(int nStates, const SpeciesThermoF* all, const State* st, float* X_out, float* DA, float* DB)
{
    const int k = blockDim.x*blockIdx.x + threadIdx.x;
    if (k >= nStates) return;
    const State& s = st[k];
    const int n = s.n;
    SpeciesThermoF sp[THERMO_MAX_SPECIES];
    for (int q = 0; q < n; ++q) sp[q] = all[s.idx[q]];
    // species_Dmix_probe_d と同じ組成
    const flow_float inv_ro = 1.0f/max(s.ro, (flow_float)1.0e-30f);
    flow_float Yf[THERMO_MAX_SPECIES], X[THERMO_MAX_SPECIES];
    flow_float ysum = 0.0f;
    for (int q = 0; q < n; q++) { flow_float y = s.roY[q]*inv_ro; if (y < 0.0f) y = 0.0f; Yf[q] = y; ysum += y; }
    const flow_float yinv = 1.0f/(ysum>1.0e-30f?ysum:1.0e-30f);
    for (int q = 0; q < n; q++) Yf[q] *= yinv;
    const flow_float gl = (s.iw >= 0) ? s.rog*inv_ro*yinv : 0.0f;
    // species_transport_X_f (speciesTransport_d.cu) と同じ分岐
    if (s.iw >= 0 && gl > 0.0f) {
        flow_float Yg[THERMO_MAX_SPECIES];
        for (int q = 0; q < n; ++q) Yg[q] = Yf[q];
        gas_phase_composition(Yg, s.iw, gl);
        thermo_X_from_Y_f(sp, n, Yg, X);
    } else {
        thermo_X_from_Y_f(sp, n, Yf, X);
    }
    for (int q = 0; q < n; q++) {
        X_out[(size_t)k*THERMO_MAX_SPECIES + q] = X[q];
        DA[(size_t)k*THERMO_MAX_SPECIES + q] = dmix_A(sp, n, X, q, s.T, s.P);
        DB[(size_t)k*THERMO_MAX_SPECIES + q] = thermo_Dmix_species_f(sp, n, X, q, s.T, s.P);
    }
}

int main(int argc, char** argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s in.txt out.txt\n", argv[0]); return 2; }
    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
    int nSp = 0, nStates = 0;
    if (fscanf(f, "%d %d", &nSp, &nStates) != 2) return 2;
    std::vector<SpeciesThermoF> all(nSp);
    for (int i = 0; i < nSp; ++i) {
        double MW, sig, eps;
        if (fscanf(f, "%lf %lf %lf", &MW, &sig, &eps) != 3) return 2;
        SpeciesThermoF s{};   // 係数は使わない (0)
        s.MW = (float)MW; s.invMW = (float)(1.0/MW); s.R = (float)(THERMO_RU/MW);   // thermo_d.cu の float ミラーと同じ丸め
        s.sigma_LJ = (float)sig; s.eps_kB = (float)eps;
        all[i] = s;
    }
    std::vector<State> st(nStates);
    for (int k = 0; k < nStates; ++k) {
        State& s = st[k];
        if (fscanf(f, "%d", &s.n) != 1 || s.n < 1 || s.n > THERMO_MAX_SPECIES) return 2;
        for (int q = 0; q < s.n; ++q) if (fscanf(f, "%d", &s.idx[q]) != 1) return 2;
        double T, P, ro, rog;
        if (fscanf(f, "%d %lf %lf %lf %lf", &s.iw, &T, &P, &ro, &rog) != 5) return 2;
        s.T = (float)T; s.P = (float)P; s.ro = (float)ro; s.rog = (float)rog;
        for (int q = 0; q < s.n; ++q) { double v; if (fscanf(f, "%lf", &v) != 1) return 2; s.roY[q] = (float)v; }
    }
    fclose(f);

    SpeciesThermoF* all_d; State* st_d; float *X_d, *DA_d, *DB_d;
    const size_t m = (size_t)nStates*THERMO_MAX_SPECIES;
    CK(cudaMalloc(&all_d, nSp*sizeof(SpeciesThermoF)));
    CK(cudaMalloc(&st_d, nStates*sizeof(State)));
    CK(cudaMalloc(&X_d, m*sizeof(float))); CK(cudaMalloc(&DA_d, m*sizeof(float))); CK(cudaMalloc(&DB_d, m*sizeof(float)));
    CK(cudaMemcpy(all_d, all.data(), nSp*sizeof(SpeciesThermoF), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(st_d, st.data(), nStates*sizeof(State), cudaMemcpyHostToDevice));
    eval<<<(nStates + 127)/128, 128>>>(nStates, all_d, st_d, X_d, DA_d, DB_d);
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    std::vector<float> X(m), DA(m), DB(m);
    CK(cudaMemcpy(X.data(), X_d, m*sizeof(float), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(DA.data(), DA_d, m*sizeof(float), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(DB.data(), DB_d, m*sizeof(float), cudaMemcpyDeviceToHost));

    FILE* o = fopen(argv[2], "w");
    if (!o) return 2;
    for (int k = 0; k < nStates; ++k) {
        const int n = st[k].n;
        const size_t b = (size_t)k*THERMO_MAX_SPECIES;
        for (int q = 0; q < n; ++q) fprintf(o, "%.9g ", X[b+q]);
        for (int q = 0; q < n; ++q) fprintf(o, "%.9g ", DA[b+q]);
        for (int q = 0; q < n; ++q) fprintf(o, "%.9g%s", DB[b+q], q + 1 < n ? " " : "\n");
    }
    fclose(o);
    CK(cudaFree(all_d)); CK(cudaFree(st_d)); CK(cudaFree(X_d)); CK(cudaFree(DA_d)); CK(cudaFree(DB_d));
    return 0;
}
