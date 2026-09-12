// =============================================================================
// test_cond_kantrowitz_carrier.cu  (host + device 単体検証; nvcc)
//   plans/active/condensation-kantrowitz-carrier.md §5(5)/§6:
//   (a) 純蒸気極限の解析値: mode 2 = (b−½)²/(c̃+½), mode 3 = 同 + 表面項 (1e-12)
//   (b) Wysłouzil 条件 (230 K, ln S 3.4, Y_w 0.01095, N2 carrier): θ1=167.4 / θ2=4.05 / θ3=2.98 (±1 %)
//   (c) 掃引 T 200–260 K × ln S 2–6 × Y_w 0.005–0.05 × g/Y_w 0–0.99: θ 有限・非負、Yv→0 で θ→0、Yc→0 で純蒸気形、
//       局所 J 序列 J0 >= J3 >= J2 >= J1、q̂>0 (適用域)
//   (d) mode 0/1 と sigmaScale 1.0 で旧関数値とビット同一 (mode 1 の旧式を再実装して比較)
//   (e) 同じ入力を device で評価し host double と一致 (1e-12)
//   build: nvcc --expt-relaxed-constexpr -I. -o test_kwc tests/unit/test_cond_kantrowitz_carrier.cu
// =============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "../../cuda_forge/condensationSource_d.cuh"
#include "../../cuda_forge/thermo_d.cuh"

static int nfail = 0;
static void check(const char* name, double a, double b, double tol) {
    const double rel = std::fabs(a - b)/std::max(1e-300, std::fabs(a) + std::fabs(b))*2.0;
    const bool ok = rel < tol; if (!ok) ++nfail;
    printf("  [%s] %-46s a=% .8e b=% .8e rel=%.2e\n", ok ? "PASS" : "FAIL", name, a, b, rel);
}
static SpeciesThermo makeSp(double MW, const double lo[9], const double hi[9]) {
    SpeciesThermo s{}; s.MW = MW; s.sigma_LJ = 3.6; s.eps_kB = 97.0; s.Tlo = 200.0; s.Tmid = 1000.0; s.Thi = 6000.0;
    for (int i = 0; i < 9; ++i) { s.low[i] = lo[i]; s.high[i] = hi[i]; } return s;
}
static const double N2_lo[9] = {2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1};
static const double N2_hi[9] = {5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1};
static const double H2O_lo[9] = {-3.947960830e4,5.755731020e2,9.317826530e-1,7.222712860e-3,-7.342557370e-6,4.955043490e-9,-1.336933246e-12,-3.303974310e4,1.724205775e1};
static const double H2O_hi[9] = {1.034972096e6,-2.412698562e3,4.646110780,2.291998307e-3,-6.836830480e-7,9.426468930e-11,-4.821580530e-15,-1.384286509e4,-7.978148510e0};

// kernel 側と同じ集計 (condensationSource_d.cu の car 構築を host で再現)
static CondNucCarrier build_car(const SpeciesThermo* sp, int n, int iv, const double* Y, double g, double T) {
    CondNucCarrier c; const double Mv = sp[iv].MW;
    c.cvv_tilde = (thermo_cp_mass(sp[iv], T) - thermo_R_species(sp[iv]))*Mv/COND_RU;
    c.a_v = std::max(Y[iv] - g, 0.0)/Mv; double sum = 0.0;
    for (int i = 0; i < n; ++i) { if (i == iv || !(Y[i] > 0.0)) continue;
        const double Mi = sp[i].MW, cvi = (thermo_cp_mass(sp[i], T) - thermo_R_species(sp[i]))*Mi/COND_RU;
        sum += (Y[i]/Mi)*std::sqrt(Mv/Mi)*(cvi + 0.5); }
    c.carrierSum = sum; return c;
}
// 旧 mode 1 (2026-09-10 版) の再実装
static double theta1_old(const CondSpeciesProps& cp, double T, double gamma_gas) {
    const double b = cond_latent(cp, T)/(cp.R*T); const double th = (2.0*(gamma_gas-1.0)/(gamma_gas+1.0))*b*(b-0.5); return th > 0.0 ? th : 0.0;
}
__global__ void theta_dev(CondSpeciesProps cp, double T, double lnS, int mode, double gam, CondNucCarrier car, double p_v, double rho_v, double* out) {
    out[0] = cond_kantrowitz_theta(cp, T, lnS, mode, gam, &car);
    double J, r; cond_nucleation(cp, T, p_v, rho_v, &J, &r, mode, gam, &car); out[1] = J; out[2] = r;
}

int main() {
    const CondSpeciesProps h2o = condProps_H2O();
    SpeciesThermo sp[2] = { makeSp(0.0280134, N2_lo, N2_hi), makeSp(0.0180153, H2O_lo, H2O_hi) };
    const double gv = h2o.cp/h2o.cv;
    printf("== (a) pure-vapor limits ==\n");
    { const double T = 230.0, lnS = 3.4; const double b = cond_latent(h2o, T)/(h2o.R*T);
      CondNucCarrier pure; pure.a_v = 1.0; pure.carrierSum = 0.0; pure.cvv_tilde = h2o.cv*h2o.M/COND_RU;
      check("mode2 pure = (b-1/2)^2/(c+1/2)", cond_kantrowitz_theta(h2o, T, lnS, 2, gv, &pure), (b-0.5)*(b-0.5)/(pure.cvv_tilde+0.5), 1e-12);
      check("mode3 pure = (b-1/2-lnS)^2/(c+1/2)", cond_kantrowitz_theta(h2o, T, lnS, 3, gv, &pure), (b-0.5-lnS)*(b-0.5-lnS)/(pure.cvv_tilde+0.5), 1e-12);
      check("mode2 nullptr == pure struct", cond_kantrowitz_theta(h2o, T, lnS, 2, gv, nullptr), cond_kantrowitz_theta(h2o, T, lnS, 2, gv, &pure), 1e-15);
      printf("      b=%.3f c~v,v=%.3f  theta1=%.3f theta2(pure)=%.3f theta3(pure)=%.3f\n", b, pure.cvv_tilde,
             cond_kantrowitz_theta(h2o, T, lnS, 1, gv, nullptr), cond_kantrowitz_theta(h2o, T, lnS, 2, gv, &pure), cond_kantrowitz_theta(h2o, T, lnS, 3, gv, &pure)); }
    printf("== (b) Wyslouzil condition ==\n");
    { const double T = 230.0, lnS = 3.4, Y[2] = {0.98905, 0.01095}; const CondNucCarrier car = build_car(sp, 2, 1, Y, 0.0, T);
      const double t1 = cond_kantrowitz_theta(h2o, T, lnS, 1, gv, &car), t2 = cond_kantrowitz_theta(h2o, T, lnS, 2, gv, &car), t3 = cond_kantrowitz_theta(h2o, T, lnS, 3, gv, &car);
      printf("      a_v=%.4e carrierSum=%.4e cvv~=%.3f  theta1=%.3f theta2=%.4f theta3=%.4f  J/Jiso: %.4f %.4f %.4f\n", car.a_v, car.carrierSum, car.cvv_tilde, t1, t2, t3, 1/(1+t1), 1/(1+t2), 1/(1+t3));
      check("theta1 = 167.4 (codex recomputation)", t1, 167.404, 1e-2); check("theta2 = 4.05", t2, 4.046, 1e-2); check("theta3 = 2.98", t3, 2.982, 1e-2);
      // (e) device
      double *d; cudaMalloc(&d, 3*sizeof(double)); double h[3];
      const double pv = std::exp(lnS)*cond_psat(h2o, T), rv = pv/(h2o.R*T);
      for (int m = 1; m <= 3; ++m) { theta_dev<<<1,1>>>(h2o, T, lnS, m, gv, car, pv, rv, d); cudaMemcpy(h, d, 3*sizeof(double), cudaMemcpyDeviceToHost);
        double J, r; cond_nucleation(h2o, T, pv, rv, &J, &r, m, gv, &car);
        char nm[64]; snprintf(nm, 64, "device theta mode %d == host", m); check(nm, h[0], cond_kantrowitz_theta(h2o, T, lnS, m, gv, &car), 1e-12);
        snprintf(nm, 64, "device J mode %d == host", m); check(nm, h[1], J, 1e-12); }
      cudaFree(d); }
    printf("== (c) sweep ==\n");
    { int bad = 0, ntot = 0; double qmin = 1e300, tmax = 0;
      for (double T : {200.0, 215.0, 230.0, 245.0, 260.0}) for (double lnS : {2.0, 3.0, 4.0, 5.0, 6.0}) for (double Yw : {0.005, 0.011, 0.02, 0.05}) for (double f : {0.0, 0.5, 0.9, 0.99}) {
        const double Y[2] = {1.0 - Yw, Yw}; const double g = f*Yw; const CondNucCarrier car = build_car(sp, 2, 1, Y, g, T);
        const double pv = std::exp(lnS)*cond_psat(h2o, T), rv = pv/(h2o.R*T);
        double th[4], J[4], r;
        for (int m = 0; m < 4; ++m) { th[m] = cond_kantrowitz_theta(h2o, T, lnS, m, gv, &car); cond_nucleation(h2o, T, pv, rv, &J[m], &r, m, gv, &car); }
        const double b = cond_latent(h2o, T)/(h2o.R*T); qmin = std::min(qmin, b - 0.5 - lnS); tmax = std::max(tmax, th[1]);
        bool ok = std::isfinite(th[2]) && std::isfinite(th[3]) && th[2] >= 0 && th[3] >= 0 && th[3] <= th[2] && th[2] <= th[1] && J[0] >= J[3] && J[3] >= J[2] && J[2] >= J[1];
        // Yv → 0 で θ → 0 (単調減少)
        const CondNucCarrier car99 = build_car(sp, 2, 1, Y, 0.999*Yw, T); if (!(cond_kantrowitz_theta(h2o, T, lnS, 3, gv, &car99) <= th[3])) ok = false;
        ++ntot; if (!ok) { ++bad; printf("      FAIL T=%g lnS=%g Yw=%g g/Yw=%g th=%g %g %g %g\n", T, lnS, Yw, f, th[0], th[1], th[2], th[3]); }
      }
      // Yc → 0 で純蒸気形に一致
      { const double T = 230.0, lnS = 3.4, Y[2] = {0.0, 1.0}; const CondNucCarrier car = build_car(sp, 2, 1, Y, 0.0, T); const double b = cond_latent(h2o, T)/(h2o.R*T);
        check("Yc->0: mode2 = pure Feder", cond_kantrowitz_theta(h2o, T, lnS, 2, gv, &car), (b-0.5)*(b-0.5)/(car.cvv_tilde+0.5), 1e-12); }
      { const double T = 230.0, lnS = 3.4, Y[2] = {0.98905, 0.01095}; const CondNucCarrier car = build_car(sp, 2, 1, Y, 0.01095, T);
        check("Yv->0: theta3 -> 0", cond_kantrowitz_theta(h2o, T, lnS, 3, gv, &car) + 1.0, 1.0, 1e-15); }
      printf("      %d/%d states ok; min q^ = %.2f (>0: Feder 形の適用域内), max theta1 = %.1f\n", ntot - bad, ntot, qmin, tmax);
      check("sweep all ok", bad == 0 ? 1.0 : 0.0, 1.0, 1e-12); check("q^ > 0 over sweep", qmin > 0 ? 1.0 : 0.0, 1.0, 1e-12); }
    printf("== (d) legacy paths bit-identical ==\n");
    { int nb = 0; for (double T : {200.0, 230.0, 260.0}) for (double lnS : {2.0, 4.0, 6.0}) {
        const double pv = std::exp(lnS)*cond_psat(h2o, T), rv = pv/(h2o.R*T); double J1, r1, J0, r0, Jn, rn;
        cond_nucleation(h2o, T, pv, rv, &J1, &r1, 1, 1.40, nullptr); const double th_old = theta1_old(h2o, T, 1.40);
        cond_nucleation(h2o, T, pv, rv, &J0, &r0, 0, 1.40, nullptr); cond_nucleation(h2o, T, pv, rv, &Jn, &rn, 0, 1.40, nullptr);
        if (J1 != J0*(1.0/(1.0+th_old))) ++nb;   // 旧式: corr=1; corr/=(1+θ); J=K e^x corr → J0*(1/(1+θ)) と厳密一致
        if (cond_sigma(h2o, T) != h2o.sigmaScale*h2o_sigma(T)) ++nb; }
      check("mode1 == old formula (bitwise) and sigmaScale 1.0 identity", nb == 0 ? 1.0 : 0.0, 1.0, 1e-12);
      CondSpeciesProps hs = h2o; hs.sigmaScale = 1.03; check("sigmaScale 1.03 scales sigma", cond_sigma(hs, 230.0), 1.03*h2o_sigma(230.0), 1e-15); }
    printf("%s (%d failures)\n", nfail ? "FAILED" : "ALL PASS", nfail); return nfail ? 1 : 0;
}
