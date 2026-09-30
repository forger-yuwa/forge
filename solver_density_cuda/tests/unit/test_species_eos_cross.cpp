// =============================================================================
// test_species_eos_cross.cpp  (host-only 単体検証ハーネス)
//   案C EOS クロス応答 species_eos_cross_response() の解析式を、圧力の
//   中央差分 (有限差分は検証専用) と突き合わせる。
//
//   検証点: p(ρ, e固定, Y) = ρ R_mix(Y) T(e,Y) を組成接空間方向 d (Σ d=0) に
//   ε 摂動し、δp_FD = [p(Y+εd) - p(Y-εd)]/(2ε) を解析 δp_Y (z=ρd) と比較。
//   保存内部エネルギー e は固定 (T は e から毎回反転)。
//
//   build: nvcc/g++ host コンパイル。tests/unit/run_unit.sh 参照。
// =============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "../../cuda_forge/species_eos_coupling_d.cuh"

// NASA-9 (低 200-1000K / 高 1000-6000K)。O2/N2 は gen_ternary_ic.py と同一。
static SpeciesThermo makeSp(double MW, const double lo[9], const double hi[9]) {
    SpeciesThermo s{}; s.MW = MW; s.sigma_LJ = 3.6; s.eps_kB = 97.0;
    thermo_set_nasa9_2(s, 200.0, 1000.0, 6000.0, lo, hi);
    return s;
}
static const double O2_lo[9] = {-3.425563420e4,4.847000970e2,1.119010961,4.293889240e-3,-6.836300520e-7,-2.023372700e-9,1.039040018e-12,-3.391454870e3,1.849699470e1};
static const double O2_hi[9] = {-1.037939022e6,2.344830282e3,1.819732036,1.267847582e-3,-2.188067988e-7,2.053719572e-11,-8.193467050e-16,-1.689010929e4,1.738716506e1};
static const double N2_lo[9] = {2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1};
static const double N2_hi[9] = {5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1};
// He: 単原子, cp/R=2.5 一定 (a2=2.5, a7 で h datum)。CEA He。
static const double He_lo[9] = {0,0,2.5,0,0,0,0,-745.375,0.928723974};
static const double He_hi[9] = {0,0,2.5,0,0,0,0,-745.375,0.928723974};
// H2O: CEA NASA-9。
static const double H2O_lo[9] = {-3.947960830e4,5.755731020e2,9.317826530e-1,7.222712860e-3,-7.342557370e-6,4.955043490e-9,-1.336933246e-12,-3.303974310e4,1.724205775e1};
static const double H2O_hi[9] = {1.034972096e6,-2.412698562e3,4.646110780,2.291998307e-3,-6.836830480e-7,9.426468930e-11,-4.821580530e-15,-1.384286509e4,-7.978148510e0};

static int g_fail = 0;
static void check(const std::string& name, double analytic, double fd, double tol) {
    const double denom = std::max(1.0, std::fabs(analytic) + std::fabs(fd));
    const double rel = std::fabs(analytic - fd) / denom;
    const bool ok = rel < tol;
    if (!ok) ++g_fail;
    printf("  [%s] %-28s analytic=% .6e  FD=% .6e  rel=%.2e  %s\n",
           ok ? "PASS" : "FAIL", name.c_str(), analytic, fd, rel, ok ? "" : "<-- FAIL");
}

// p(ρ, e固定, Y): T を e から反転して p=ρ R_mix T。
static double pressure(const std::vector<SpeciesThermo>& sp, const std::vector<double>& Y,
                       double ro, double e, double Tguess) {
    const int n = (int)sp.size();
    const double T = thermo_T_from_e(sp.data(), n, Y.data(), e, Tguess, 50.0, 8000.0);
    return ro * thermo_R_mix(sp.data(), n, Y.data()) * T;
}

// 1 ケース: 状態(sp,Y0,ro,T0) と組成方向 dir(Σ=0) で解析 δp と 中央差分 δp を比較。
static void runCase(const std::string& name, const std::vector<SpeciesThermo>& sp,
                    std::vector<double> Y0, double ro, double T0, std::vector<double> dir) {
    const int n = (int)sp.size();
    // 方向を Σ dir = 0 に射影 (現組成で分配, 案C の接空間射影と同じ規約)
    double s = 0.0; for (int i = 0; i < n; ++i) s += dir[i];
    for (int i = 0; i < n; ++i) dir[i] -= Y0[i] * s;   // dir is δY 方向 (Σ=0)
    // 保存内部エネルギー e (= e_mix(T0,Y0)) を固定
    const double e = thermo_e_mix(sp.data(), n, Y0.data(), T0);
    // 解析 δp_Y: z = ρ δY = ρ dir
    std::vector<double> z(n); for (int i = 0; i < n; ++i) z[i] = ro * dir[i];
    double dpY = 0.0, dTY = 0.0;
    species_eos_cross_response(sp.data(), n, Y0.data(), T0, ro, z.data(), &dpY, &dTY);
    // 中央差分: Y ± ε dir, e/ρ 固定で p を評価
    const double eps = 1.0e-6;
    std::vector<double> Yp(n), Ym(n);
    for (int i = 0; i < n; ++i) { Yp[i] = Y0[i] + eps * dir[i]; Ym[i] = Y0[i] - eps * dir[i]; }
    const double pp = pressure(sp, Yp, ro, e, T0);
    const double pm = pressure(sp, Ym, ro, e, T0);
    const double dp_fd = (pp - pm) / (2.0 * eps);
    // δT も検証 (T(e,Y) の中央差分)
    const double Tp = thermo_T_from_e(sp.data(), n, Yp.data(), e, T0, 50.0, 8000.0);
    const double Tm = thermo_T_from_e(sp.data(), n, Ym.data(), e, T0, 50.0, 8000.0);
    const double dT_fd = (Tp - Tm) / (2.0 * eps);
    printf("Case %s (T0=%.0fK ro=%.3f):\n", name.c_str(), T0, ro);
    check(name + ":dp", dpY, dp_fd, 1.0e-5);
    check(name + ":dT", dTY, dT_fd, 1.0e-5);
}

int main() {
    SpeciesThermo He = makeSp(0.0040026, He_lo, He_hi);
    SpeciesThermo O2 = makeSp(0.0319988, O2_lo, O2_hi);
    SpeciesThermo N2 = makeSp(0.0280134, N2_lo, N2_hi);
    SpeciesThermo H2O = makeSp(0.0180153, H2O_lo, H2O_hi);

    // He/O2/N2 (Cutler): 冷流 300K と中温 800K, He コア寄りの組成
    runCase("HeO2N2_300", {He,O2,N2}, {0.5,0.1,0.4}, 0.5, 300.0, {1.0,-0.5,-0.5});
    runCase("HeO2N2_800", {He,O2,N2}, {0.5,0.1,0.4}, 0.5, 800.0, {1.0,-0.5,-0.5});
    runCase("HeO2N2_air", {He,O2,N2}, {0.05,0.22,0.73}, 1.1, 300.0, {0.3,-0.1,-0.2});
    // hot N2/H2O (wys 系: 生成エンタルピーの大きい H2O)
    runCase("N2H2O_1500", {N2,H2O}, {0.95,0.05}, 0.3, 1500.0, {-1.0,1.0});
    runCase("N2H2O_500",  {N2,H2O}, {0.99,0.01}, 0.6, 500.0,  {-1.0,1.0});
    // binary N2/O2 (空気近傍)
    runCase("N2O2_300",   {N2,O2}, {0.767,0.233}, 1.18, 300.0, {1.0,-1.0});

    printf("\n%s (%d failures)\n", g_fail == 0 ? "ALL PASS" : "SOME FAILED", g_fail);
    return g_fail == 0 ? 0 : 1;
}
