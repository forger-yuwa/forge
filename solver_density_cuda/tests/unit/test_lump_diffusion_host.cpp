// test_lump_diffusion_host.cpp — lump を含む化学種拡散の縮約 thermo_Dmix_lumped_f (cuda_forge/thermo_d.cuh) の host 単体試験
//   plan thermophysics-solver-owned-species-db §4.4 確定版, #7b, §6 V4a/V4b/V4e。参照式は tests/unit/test_lump_diffusion_reduction.py と同じ。
//   (1) 非 lump ラベルの係数が「同じ実種を全部輸送種にした config」の現行 thermo_Dmix_species_f とビット一致 (va3: MIXDRY + H2O の H2O)
//   (2) 全ラベルの係数が独立 double 参照と相対 ≤1e-5 (va3・SERN の EXH/AMB 重複・純成分・微量種・T* 両端)
//   (3) 非重複の恒等式 (double 参照で lump の補正後流束 = 構成実種の和、外部種 = full、≤1e-12)
// build: g++ -O2 -std=c++17 -I. -o test_lump_diffusion_host tests/unit/test_lump_diffusion_host.cpp
//        (solver_density_cuda で。CUDA は使わない)
#include "cuda_forge/thermo_d.cuh"
#include <cstdio>
#include <cmath>
#include <map>
#include <string>
#include <vector>

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } else { printf("  ok  : " __VA_ARGS__); printf("\n"); } } while (0)

struct Real { double MW, sig, eps; };
// 設計側共通データの既定 ljSource [gri30, svehla1962] の値 (test_lump_diffusion_reduction.py が読むものと同じ)
static const std::map<std::string, Real> RDB = {
    {"N2", {0.0280134, 3.621, 97.53}}, {"H2O", {0.01801528, 2.605, 572.4}}, {"H2", {0.00201588, 2.92, 38.0}},
    {"AR", {0.039948, 3.33, 136.5}},   {"OH", {0.01700734, 2.75, 80.0}},    {"O2", {0.0319988, 3.458, 107.4}},
    {"NO", {0.0300061, 3.621, 97.53}}, {"H", {0.00100794, 2.05, 145.0}},    {"O", {0.0159994, 2.75, 80.0}},
    {"CO2", {0.0440095, 3.763, 244.0}}, {"CO", {0.0280101, 3.65, 98.1}}};

using Comp = std::vector<std::pair<std::string, double>>;   // lump の構成 (モル)。非 lump は {名前, 1}

struct Setup {
    LumpDiffD ld{};
    std::vector<std::string> reals;
    std::vector<double> ML;                 // ラベルのモル質量
    std::vector<std::vector<double>> E, A;  // double 版
};

static Setup build(const std::vector<Comp>& labels)
{
    Setup S;
    auto idx = [&](const std::string& k) {
        for (size_t r = 0; r < S.reals.size(); ++r) if (S.reals[r] == k) return (int)r;
        const int r = (int)S.reals.size(); S.reals.push_back(k);
        const Real& R = RDB.at(k);
        S.ld.MW[r] = (float)R.MW; S.ld.sig[r] = (float)R.sig; S.ld.eps[r] = (float)R.eps;
        return r;
    };
    S.E.assign(labels.size(), std::vector<double>(THERMO_MAX_DIFF_REAL, 0.0));
    S.A = S.E;
    for (size_t s = 0; s < labels.size(); ++s) {
        double t = 0.0, M = 0.0;
        for (auto& c : labels[s]) t += c.second;
        for (auto& c : labels[s]) M += c.second / t * RDB.at(c.first).MW;
        S.ML.push_back(M);
        for (auto& c : labels[s]) {
            const int r = idx(c.first);
            const double x = c.second / t;
            S.E[s][r] += x; S.A[s][r] += x * RDB.at(c.first).MW / M;
            S.ld.E[s][r] += (float)x; S.ld.A[s][r] += (float)(x * RDB.at(c.first).MW / M);
        }
    }
    S.ld.nReal = (int)S.reals.size();
    return S;
}

// 独立 double 参照 (実種ごとに q を回す素直な形)
static double omega11d(double Ts) { Ts = std::fmin(std::fmax(Ts, 0.3), 100.0);
    return 1.06036 * std::pow(Ts, -0.15610) + 0.19300 * std::exp(-0.47635 * Ts) + 1.03587 * std::exp(-1.52996 * Ts) + 1.76474 * std::exp(-3.89411 * Ts); }
static double dbin(const Real& a, const Real& b, double T, double P) {
    const double sig = 0.5 * (a.sig + b.sig), eps = std::sqrt(a.eps * b.eps);
    return 1.8583e-3 * std::sqrt(T * T * T * (1.0 / (a.MW * 1e3) + 1.0 / (b.MW * 1e3))) / ((P / 101325.0) * sig * sig * omega11d(T / eps)) * 1e-4; }
static std::vector<double> refD(const Setup& S, const std::vector<double>& XL, double T, double P, std::vector<double>* Dr_out = nullptr)
{
    const int nr = (int)S.reals.size(), nl = (int)XL.size();
    std::vector<double> Xr(nr, 0.0), Dr(nr);
    for (int r = 0; r < nr; ++r) for (int s = 0; s < nl; ++s) Xr[r] += XL[s] * S.E[s][r];
    for (int r = 0; r < nr; ++r) {
        double num = 0, den = 0;
        for (int q = 0; q < nr; ++q) if (q != r) { const double d = dbin(RDB.at(S.reals[r]), RDB.at(S.reals[q]), T, P); num += Xr[q]; den += Xr[q] / d; }
        Dr[r] = (den < 1e-300) ? dbin(RDB.at(S.reals[r]), RDB.at(S.reals[r]), T, P) : num / den;
    }
    std::vector<double> D(nl, 0.0);
    for (int s = 0; s < nl; ++s) for (int r = 0; r < nr; ++r) D[s] += S.A[s][r] * Dr[r];
    if (Dr_out) *Dr_out = Dr;
    return D;
}
static std::vector<double> X_from_Y(const Setup& S, const std::vector<double>& Y) {
    std::vector<double> w(Y.size()); double t = 0;
    for (size_t s = 0; s < Y.size(); ++s) { w[s] = Y[s] / S.ML[s]; t += w[s]; }
    for (auto& v : w) v /= t;
    return w; }

static double max_rel_float_vs_ref(const Setup& S, const std::vector<double>& Y, double T, double P, const char* tag)
{
    const std::vector<double> XL = X_from_Y(S, Y);
    float Xf[THERMO_MAX_SPECIES], Dl[THERMO_MAX_SPECIES];
    for (size_t s = 0; s < XL.size(); ++s) Xf[s] = (float)XL[s];
    thermo_Dmix_lumped_f(S.ld, (int)XL.size(), Xf, (float)T, (float)P, Dl);
    const std::vector<double> D = refD(S, XL, T, P);
    double m = 0;
    for (size_t s = 0; s < D.size(); ++s) m = std::fmax(m, std::fabs(Dl[s] - D[s]) / D[s]);
    (void)tag;
    return m;
}

int main()
{
    printf("[V4b] 非 lump ラベルの係数 = 全実種 config の現行 thermo_Dmix_species_f (ビット一致)\n");
    const Comp mixdry = {{"N2", 6.64860e-1}, {"O2", 2.16072e-1}, {"AR", 7.97588e-3}, {"CO2", 4.90034e-2}};
    Setup va3 = build({mixdry, {{"H2O", 1.0}}});
    {
        // full config: [N2, O2, AR, CO2, H2O] (実種の順序は lump の展開順と同じ)
        SpeciesThermoF spf[5] = {};
        const char* nm[5] = {"N2", "O2", "AR", "CO2", "H2O"};
        for (int i = 0; i < 5; ++i) { const Real& R = RDB.at(nm[i]); spf[i].MW = (float)R.MW; spf[i].sigma_LJ = (float)R.sig; spf[i].eps_kB = (float)R.eps; }
        int nbad = 0, ntot = 0;
        for (double T : {200.0, 300.0, 1000.0, 2500.0}) for (double Yw : {1e-8, 0.04, 0.2, 0.9}) {
            const std::vector<double> XL = X_from_Y(va3, {1 - Yw, Yw});
            float Xf[2] = {(float)XL[0], (float)XL[1]}, Dl[2];
            thermo_Dmix_lumped_f(va3.ld, 2, Xf, (float)T, 101325.0f, Dl);
            float Xr[5];
            for (int r = 0; r < 5; ++r) { float x = 0.0f; for (int s = 0; s < 2; ++s) x += Xf[s] * va3.ld.E[s][r]; Xr[r] = x; }
            const float Dfull = thermo_Dmix_species_f(spf, 5, Xr, 4, (float)T, 101325.0f);
            ++ntot; if (Dfull != Dl[1]) { ++nbad; printf("    T %g Y_H2O %g: lumped %.9g vs full %.9g\n", T, Yw, Dl[1], Dfull); }
        }
        CHECK(nbad == 0, "va3 H2O: %d/%d 状態でビット一致", ntot - nbad, ntot);
    }

    printf("[V4e] float 実装 vs 独立 double 参照 (相対 ≤1e-5)\n");
    {
        double m = 0;
        for (double T : {200.0, 300.0, 1000.0, 2500.0, 6000.0}) for (double Yw : {1e-8, 0.04, 0.5, 1.0 - 1e-7}) m = std::fmax(m, max_rel_float_vs_ref(va3, {1 - Yw, Yw}, T, 101325.0, "va3"));
        CHECK(m <= 1e-5, "va3 MIXDRY + H2O: 最大相対差 %.2e", m);
        const Comp exh = {{"N2", 0.6389}, {"H2O", 0.32695}, {"H2", 0.01251}, {"AR", 0.00767}, {"OH", 0.00604}, {"O2", 0.00398},
                          {"NO", 0.00206}, {"H", 0.00125}, {"O", 0.00037}, {"CO2", 0.00021}, {"CO", 0.00005}};
        const Comp amb = {{"N2", 0.78084}, {"O2", 0.20946}, {"AR", 0.00934}};
        Setup sern = build({exh, amb});
        CHECK(sern.ld.nReal == 11, "SERN: EXH ∪ AMB の実種 %d (重複 N2/O2/AR は 1 つにまとまる)", sern.ld.nReal);
        m = 0;
        for (double T : {50.0, 300.0, 1000.0, 2000.0, 30000.0}) for (double Ye : {0.0, 1e-8, 0.1, 0.5, 0.9, 1.0}) m = std::fmax(m, max_rel_float_vs_ref(sern, {Ye, 1 - Ye}, T, 1.0e5, "sern"));
        CHECK(m <= 1e-5, "SERN EXH/AMB (T* 両端・純成分・微量含む): 最大相対差 %.2e", m);
        // 純成分 1 実種 + 複数ラベル: 自己拡散
        Setup pure = build({{{"N2", 1.0}}, {{"N2", 0.5}, {"N2", 0.5}}});
        CHECK(pure.ld.nReal == 1, "同じ実種だけのラベル 2 つは実種 1 (自己拡散)");
        m = max_rel_float_vs_ref(pure, {0.3, 0.7}, 300.0, 101325.0, "pure");
        CHECK(m <= 1e-5, "純成分の自己拡散 D_rr: 相対差 %.2e", m);
    }

    printf("[V4a] 非重複の恒等式 (double 参照; 外部種 = full、lump = 構成和、≤1e-12)\n");
    {
        double worst = 0;
        for (double T : {300.0, 1000.0}) for (double Yw : {0.04, 0.2}) {
            const std::vector<double> Y = {1 - Yw, Yw}, gY = {-1.0, 1.0};
            const std::vector<double> XL = X_from_Y(va3, Y);
            std::vector<double> Dr; const std::vector<double> D = refD(va3, XL, T, 101325.0, &Dr);
            double J[2], S = 0; for (int s = 0; s < 2; ++s) { J[s] = -D[s] * gY[s]; S += J[s]; }
            double Js[2]; for (int s = 0; s < 2; ++s) Js[s] = J[s] - Y[s] * S;
            const int nr = va3.ld.nReal; std::vector<double> Yr(nr), gYr(nr), Jr(nr); double Sr = 0;
            for (int r = 0; r < nr; ++r) { for (int s = 0; s < 2; ++s) { Yr[r] += va3.A[s][r] * Y[s]; gYr[r] += va3.A[s][r] * gY[s]; } Jr[r] = -Dr[r] * gYr[r]; Sr += Jr[r]; }
            double lumpsum = 0, scale = 0; std::vector<double> Jrs(nr);
            for (int r = 0; r < nr; ++r) { Jrs[r] = Jr[r] - Yr[r] * Sr; scale += Jrs[r] * Jrs[r]; if (r < 4) lumpsum += Jrs[r]; }
            scale = std::sqrt(scale);
            worst = std::fmax(worst, std::fmax(std::fabs(Js[1] - Jrs[4]), std::fabs(Js[0] - lumpsum)) / scale);
        }
        CHECK(worst <= 1e-12, "va3: 最大規格化誤差 %.2e", worst);
    }
    printf(g_fail ? "\nFAIL (%d)\n" : "\nALL PASS\n", g_fail);
    return g_fail ? 1 : 0;
}
