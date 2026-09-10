// =============================================================================
// test_cond_sonic.cpp  (host-only 単体検証)
//   plans/active/condensation-kantrowitz-gamma-twophase-sonic.md §6:
//   (a) CPG pure (N2 / H2O 物性): 一温度二相 EOS p(ρ,e,g) (cond_T_from_e_cpg) の有限差分
//       (∂p/∂ρ)_e + (p/ρ²)(∂p/∂e)_ρ と cond_twophase_sonic の c² が一致 (rel<1e-5)
//   (b) TP carrier (NASA-9 N2 + H2O, cond_T_from_e_carrier) で同様
//   (c) g=0 で従来式 γ_mix R_mix T と一致 (diff 0)
//   (b') 掃引 T×Y_w×g で c_v,2φ>0, c²>0, g>0 では c_2φ<c^全蒸気 (厳密不等号), g=0 は等号; L' 刻み依存; float32 γ−1 精度
//   (d) Kantrowitz θ: γ_v=1.331 (H2O) は旧 1.40 比で 0.85 倍
//   (f) resolveCondSonicModel の選択条件
//   build: nvcc -x cu --expt-relaxed-constexpr -o test_cond_sonic tests/unit/test_cond_sonic.cpp
// =============================================================================
#include <cstdio>
#include <cmath>
#include "../../cuda_forge/condensationEOS_d.cuh"
#include "../../cuda_forge/condensationSource_d.cuh"
#include "../../input/condSonicResolve.hpp"

static int nfail = 0;
static void check(const char* name, double a, double b, double tol) {
    const double rel = std::fabs(a - b)/std::max(1e-300, std::fabs(a) + std::fabs(b))*2.0;
    const bool ok = rel < tol; if (!ok) ++nfail;
    printf("  [%s] %-40s a=% .8e b=% .8e rel=%.2e\n", ok ? "PASS" : "FAIL", name, a, b, rel);
}

static SpeciesThermo makeSp(double MW, const double lo[9], const double hi[9]) {
    SpeciesThermo s{}; s.MW = MW; s.sigma_LJ = 3.6; s.eps_kB = 97.0;
    s.Tlo = 200.0; s.Tmid = 1000.0; s.Thi = 6000.0;
    for (int i = 0; i < 9; ++i) { s.low[i] = lo[i]; s.high[i] = hi[i]; }
    return s;
}
static const double N2_lo[9] = {2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1};
static const double N2_hi[9] = {5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1};
static const double H2O_lo[9] = {-3.947960830e4,5.755731020e2,9.317826530e-1,7.222712860e-3,-7.342557370e-6,4.955043490e-9,-1.336933246e-12,-3.303974310e4,1.724205775e1};
static const double H2O_hi[9] = {1.034972096e6,-2.412698562e3,4.646110780,2.291998307e-3,-6.836830480e-7,9.426468930e-11,-4.821580530e-15,-1.384286509e4,-7.978148510e0};

// ---- (a) CPG pure: p(ρ,e,g) = (1-g) ρ R T(e,g) ----
static double p_cpg(double rho, double e, double g, double cv, double R, const CondSpeciesProps& cp, double Tg) {
    const double T = cond_T_from_e_cpg(e, g, cv, R, Tg, cp);
    return (1.0 - g)*rho*R*T;
}
static void test_cpg(const CondSpeciesProps& cp, double gamma, double T, double rho, double g, const char* tag) {
    const double cv = cp.cp/gamma, R = (gamma - 1.0)*cv;
    // e from T (EOS の順方向): e = (cv + g R) T - g L(T)
    const double e = (cv + g*R)*T - g*cond_latent(cp, T);
    const double p = p_cpg(rho, e, g, cv, R, cp, T);
    const double dr = 1e-4*rho, de = 1e-4*std::fabs(e) + 1.0;
    const double dpdr = (p_cpg(rho + dr, e, g, cv, R, cp, T) - p_cpg(rho - dr, e, g, cv, R, cp, T))/(2*dr);
    const double dpde = (p_cpg(rho, e + de, g, cv, R, cp, T) - p_cpg(rho, e - de, g, cv, R, cp, T))/(2*de);
    const double c2_fd = dpdr + p/(rho*rho)*dpde;
    const double dL = (cond_latent(cp, T + 0.1) - cond_latent(cp, T - 0.1))/0.2;
    double gam2, c2;
    cond_twophase_sonic(cp.cp, (1.0 - g)*R, g, dL, T, &gam2, &c2);
    char nm[96]; snprintf(nm, sizeof nm, "CPG %s T=%.0f g=%.3f c2 FD", tag, T, g);
    check(nm, c2, c2_fd, 1e-5);
    if (g == 0.0) { snprintf(nm, sizeof nm, "CPG %s g=0 == gamma R T", tag); check(nm, c2, gamma*R*T, 1e-15); }
    printf("      gamma_2ph=%.5f (gas %.3f)  c=%.3f (gas-only %.3f)\n", gam2, gamma, std::sqrt(c2), std::sqrt(gamma*R*T));
}

// ---- (b) TP carrier: p = ρ T (R_mix - g R_w), e = e_allvap(T) + g (R_w T - L) ----
struct TPctx { const SpeciesThermo* sp; int n; const double* Y; double Rw; CondSpeciesProps cp; };
static double p_tp(const TPctx& c, double rho, double e, double g, double Tg) {
    const double T = cond_T_from_e_carrier(c.sp, c.n, c.Y, e, g, c.Rw, c.cp, Tg, 50.0, 3000.0);
    return rho*T*(thermo_R_mix(c.sp, c.n, c.Y) - g*c.Rw);
}
static void test_tp(const TPctx& c, double T, double rho, double g) {
    double cpm, hm; thermo_cph_mix(c.sp, c.n, c.Y, T, &cpm, &hm);
    const double Rmix = thermo_R_mix(c.sp, c.n, c.Y);
    const double e = (hm - Rmix*T) + g*(c.Rw*T - cond_latent(c.cp, T));
    const double p = p_tp(c, rho, e, g, T);
    const double dr = 1e-4*rho, de = 1e-4*std::fabs(e) + 1.0;
    const double dpdr = (p_tp(c, rho + dr, e, g, T) - p_tp(c, rho - dr, e, g, T))/(2*dr);
    const double dpde = (p_tp(c, rho, e + de, g, T) - p_tp(c, rho, e - de, g, T))/(2*de);
    const double c2_fd = dpdr + p/(rho*rho)*dpde;
    const double dL = (cond_latent(c.cp, T + 0.1) - cond_latent(c.cp, T - 0.1))/0.2;
    double gam2, c2;
    cond_twophase_sonic(cpm, Rmix - g*c.Rw, g, dL, T, &gam2, &c2);
    char nm[96]; snprintf(nm, sizeof nm, "TP carrier T=%.0f g=%.4f c2 FD", T, g);
    check(nm, c2, c2_fd, 2e-5);   // Newton 収束 1e-3 K のぶん FD が粗い
    const double gmix = cpm/(cpm - Rmix);
    if (g == 0.0) check("TP carrier g=0 == gmix Rmix T", c2, gmix*Rmix*T, 1e-15);
    printf("      gamma_2ph=%.5f (allvap %.5f)  c=%.3f (allvap %.3f, %+.2f %%)\n", gam2, gmix, std::sqrt(c2), std::sqrt(gmix*Rmix*T),
           (std::sqrt(c2)/std::sqrt(gmix*Rmix*T) - 1.0)*100.0);
}

int main() {
    printf("== (a) CPG pure ==\n");
    const CondSpeciesProps n2 = condProps_N2(), h2o = condProps_H2O();
    for (double g : {0.0, 0.01, 0.1, 0.2}) test_cpg(n2, 1.4, 80.0, 0.05, g, "N2");
    for (double g : {0.0, 0.01, 0.1}) test_cpg(h2o, h2o.cp/h2o.cv, 230.0, 0.02, g, "H2O");
    printf("== (b) TP carrier (N2 + H2O) ==\n");
    SpeciesThermo sp[2] = { makeSp(0.0280134, N2_lo, N2_hi), makeSp(0.0180153, H2O_lo, H2O_hi) };
    const double Y[2] = {0.98905, 0.01095};
    TPctx c{sp, 2, Y, h2o.R, h2o};
    for (double g : {0.0, 0.005, 0.0109}) test_tp(c, 210.0, 0.15, g);
    test_tp(c, 260.0, 0.3, 0.005);
    printf("== (d) Kantrowitz theta ==\n");
    {
        const double T = 230.0, b = cond_latent(h2o, T)/(h2o.R*T);
        auto th = [&](double gm){ return 2.0*(gm - 1.0)/(gm + 1.0)*b*(b - 0.5); };
        const double gv = h2o.cp/h2o.cv;
        printf("      b=%.3f gamma_v=%.4f theta_v=%.2f theta_1.40=%.2f  J ratio (new/old)=%.4f\n",
               b, gv, th(gv), th(1.40), (1 + th(1.40))/(1 + th(gv)));
        check("gamma_v H2O in [1.325,1.335]", std::fabs(gv - 1.330) < 0.005 ? 1.0 : 0.0, 1.0, 1e-12);
        double J1, r1, J2, r2;
        const double pv = 200.0, rhov = pv/(h2o.R*T);
        cond_nucleation(h2o, T, pv, rhov, &J1, &r1, 1, gv);
        cond_nucleation(h2o, T, pv, rhov, &J2, &r2, 1, 1.40);
        check("J(gamma_v)/J(1.40) == (1+th_old)/(1+th_v)", J1/J2, (1 + th(1.40))/(1 + th(gv)), 1e-12);
    }
    printf("== (b') sweep T x Y_w x g, L' step dependence, float32 gamma-1 ==\n");
    {
        const double Ts[] = {200.0, 230.0, 260.0, 272.15, 274.15, 300.0, 350.0};
        const double Yws[] = {0.005, 0.011, 0.05};
        const double gf[] = {0.0, 0.25, 0.5, 1.0};
        int nbad = 0, ntot = 0; double worst_step = 0.0, worst_f32 = 0.0;
        for (double T : Ts) for (double Yw : Yws) for (double f : gf) {
            const double Yc[2] = {1.0 - Yw, Yw}; const double g = f*Yw;
            double cpm, hm; thermo_cph_mix(sp, 2, Yc, T, &cpm, &hm);
            const double Rmix = thermo_R_mix(sp, 2, Yc), Reff = Rmix - g*h2o.R;
            const double dL = (cond_latent(h2o, T + 0.1) - cond_latent(h2o, T - 0.1))/0.2;
            double g2, c2; const bool okf = cond_twophase_sonic(cpm, Reff, g, dL, T, &g2, &c2);
            const double gmix = cpm/(cpm - Rmix), c2v = gmix*Rmix*T;
            const double cv2 = cpm - g*dL - Reff;
            bool ok = okf && cv2 > 0.0 && c2 > 0.0 && ((g == 0.0) ? (c2 == c2v) : (c2 < c2v));
            // L' 刻み依存 (0.05 / 0.2)
            for (double h : {0.05, 0.2}) {
                const double dLh = (cond_latent(h2o, T + h) - cond_latent(h2o, T - h))/(2*h);
                double g2h, c2h; cond_twophase_sonic(cpm, Reff, g, dLh, T, &g2h, &c2h);
                const double rel = std::fabs(c2h - c2)/c2; worst_step = std::max(worst_step, rel);
                const double tol = (std::fabs(T - 273.15) < 1.5) ? 1e-4 : 1e-6;   // 液相多項式/外挿の接続点近傍は緩める
                if (rel > tol) { ok = false; printf("      step-dep fail T=%.2f g=%.4f h=%.2f rel=%.2e\n", T, g, h, rel); }
            }
            const float gf32 = (float)g2; const double relf = std::fabs(((double)gf32 - 1.0) - (g2 - 1.0))/(g2 - 1.0);
            worst_f32 = std::max(worst_f32, relf); if (relf > 1e-6) ok = false;
            ++ntot; if (!ok) { ++nbad; printf("      FAIL T=%.2f Yw=%.3f g=%.5f cv2=%.3f c2=%.3e c2v=%.3e\n", T, Yw, g, cv2, c2, c2v); }
        }
        printf("      %d/%d states ok, worst L' step dependence %.2e, worst float32 (gamma-1) rel %.2e\n", ntot - nbad, ntot, worst_step, worst_f32);
        check("sweep all states ok", nbad == 0 ? 1.0 : 0.0, 1.0, 1e-12);
    }
    printf("== (f) resolveCondSonicModel ==\n");
    {
        std::string r, w; using V = std::vector<std::string>;
        const V ok_bc = {"inlet_Pressure", "outflow", "wall", "slip"};
        check("auto verified -> 1", resolveCondSonicModel(-1, 1, 2, 1, 1, 0, ok_bc, r, w), 1.0, 1e-12);
        check("auto condEquilibrium 1 -> 0", resolveCondSonicModel(-1, 1, 2, 1, 1, 1, ok_bc, r, w), 0.0, 1e-12);
        check("auto condEquilibrium 2 -> 0", resolveCondSonicModel(-1, 1, 2, 1, 1, 2, ok_bc, r, w), 0.0, 1e-12);
        check("auto outlet_statPress -> 0", resolveCondSonicModel(-1, 1, 2, 1, 1, 0, V{"inlet_Pressure", "outlet_statPress", "wall"}, r, w), 0.0, 1e-12);
        check("auto wall_isothermal -> 0", resolveCondSonicModel(-1, 1, 2, 1, 1, 0, V{"inlet_Pressure", "outflow", "wall_isothermal"}, r, w), 0.0, 1e-12);
        check("auto CPG -> 0", resolveCondSonicModel(-1, 1, 0, 1, 1, 0, ok_bc, r, w), 0.0, 1e-12);
        check("auto pure -> 0", resolveCondSonicModel(-1, 1, 2, -1, 0, 0, ok_bc, r, w), 0.0, 1e-12);
        check("auto condensation off -> 0", resolveCondSonicModel(-1, 0, 2, 1, 1, 0, ok_bc, r, w), 0.0, 1e-12);
        check("explicit 0 kept", resolveCondSonicModel(0, 1, 2, 1, 1, 0, ok_bc, r, w), 0.0, 1e-12);
        check("explicit 1 unverified kept (+warn)", resolveCondSonicModel(1, 1, 2, 1, 1, 2, ok_bc, r, w), 1.0, 1e-12);
        check("explicit 1 unverified has warn", w.empty() ? 0.0 : 1.0, 1.0, 1e-12);
        resolveCondSonicModel(1, 1, 2, 1, 1, 0, ok_bc, r, w);
        check("explicit 1 verified no warn", w.empty() ? 1.0 : 0.0, 1.0, 1e-12);
    }
    printf("%s (%d failures)\n", nfail ? "FAILED" : "ALL PASS", nfail);
    return nfail ? 1 : 0;
}
