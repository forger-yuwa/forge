// =============================================================================
// test_cond_air.cpp  (host 単体検証; nvcc -x cu --expt-relaxed-constexpr -I. -o test_cond_air tests/unit/test_cond_air.cpp)
//   plans/active/condensation-air.md §5 (7):
//   (a) N2 潜熱/飽和圧の低温整合: 70 K / 50 K の C0 接続 (値連続)、片側微分 (c_l>0 ⇔ L'<0)、p_sat の単調増、旧一式 (lowT=0) が旧関数と bitwise、
//       38 K の p_sat 比 (新/旧) = 0.518 (閉形式)
//   (b) CPG carrier EOS の往復: e(T,g) → T (括弧付き Newton) → e、g∈[0,0.99 Y_w], T∈[25,125] K、codex の反例 (g=0.75, T=122 K) で収束
//   (c) pure (Y_w=1, R_w=R) が pure 式と一致、carrier の p, e が定義式どおり
//   (d) SLAU CPG 二相面エンタルピー h_f=c_p T_f − g L(T_f) (T_f=p_f/(ρ_f R_eff)) が同じ面状態の EOS h=e+p/ρ と一致
//   (e) slip bvar エネルギー ρE_b = ρE_i − ½ρU_n² (法線速度除去) と ghost ρE_g=ρE_i (反射)
//   (f) 空気定数 (二成分 0.79/0.21): M 28.850, R 288.19, Y_N2 0.7671, R_air − Y_w R_N2 > 0; cond_kgas の切替
// =============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "../../cuda_forge/condensationEOS_d.cuh"
#include "../../cuda_forge/condensationSource_d.cuh"

static int nfail = 0;
static void check(const char* name, double a, double b, double tol) {
    const double rel = std::fabs(a - b)/std::max(1e-300, std::fabs(a) + std::fabs(b))*2.0;
    const bool ok = rel < tol; if (!ok) ++nfail;
    printf("  [%s] %-52s a=% .8e b=% .8e rel=%.2e\n", ok ? "PASS" : "FAIL", name, a, b, rel);
}
static void checkb(const char* name, bool ok) { if (!ok) ++nfail; printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name); }

int main() {
    const double cpv = 1038.8, RN2 = 296.8;
    CondPropOpts newo{1, 1, 2000.0, 0, 1.0, -1.0}, oldo{0, 0, 2000.0, 0, 1.0, -1.0};
    const CondSpeciesProps n2n = condProps_make(COND_MODEL_N2, newo), n2o = condProps_make(COND_MODEL_N2, oldo);
    printf("== (a) N2 latent / psat low-T ==\n");
    check("L C0 at 70 K (new = poly)", cond_latent(n2n, 70.0), n2_latent_poly(70.0), 1e-14);
    check("L C0 at 70 K (new 69.999 ~ 70.001)", cond_latent(n2n, 69.999), cond_latent(n2n, 70.001), 1e-4);
    check("L'(70-) = cpv - cl = -961.2", (cond_latent(n2n, 70.0) - cond_latent(n2n, 69.9))/0.1, cpv - 2000.0, 1e-9);
    printf("      L'(70+) poly = %.1f J/kg/K (C0 only; derivative discontinuous by design)\n", (cond_latent(n2n, 70.1) - cond_latent(n2n, 70.0))/0.1);
    bool clpos = true; for (double T = 25.0; T <= 125.0; T += 0.5) { const double dL = (cond_latent(n2n, T + 0.05) - cond_latent(n2n, T - 0.05))/0.1; if (!(cpv - dL > 0.0)) { clpos = false; printf("      c_l<=0 at T=%.1f (L'=%.1f)\n", T, dL); } }
    checkb("c_l = cpv - L' > 0 over 25-125 K (new)", clpos);
    bool clold = true; for (double T = 25.0; T <= 60.0; T += 0.5) { const double dL = (cond_latent(n2o, T + 0.05) - cond_latent(n2o, T - 0.05))/0.1; if (!(cpv - dL > 0.0)) clold = false; }
    checkb("old polynomial has c_l<=0 somewhere below 60 K (documented defect)", !clold);
    check("psat C0 at 50 K (new 50-1e-9 ~ 50+1e-9)", cond_psat(n2n, 50.0 - 1e-9), cond_psat(n2n, 50.0 + 1e-9), 1e-8);
    printf("      (ln psat)'(50-) new = %.4f, (50+) Jacobsen = %.4f  1/K (C0 only)\n", (std::log(cond_psat(n2n, 50.0)) - std::log(cond_psat(n2n, 49.9)))/0.1, (std::log(cond_psat(n2n, 50.1)) - std::log(cond_psat(n2n, 50.0)))/0.1);
    bool mono = true; double prev = cond_psat(n2n, 25.0); for (double T = 25.5; T <= 125.0; T += 0.5) { const double p = cond_psat(n2n, T); if (!(p > prev)) { mono = false; printf("      psat not monotone at %.1f\n", T); } prev = p; }
    checkb("psat monotone increasing 25-125 K (new)", mono);
    check("psat(38 K) new/old = 0.5176 (closed form)", cond_psat(n2n, 38.0)/cond_psat(n2o, 38.0), 0.5176, 2e-3);
    check("old set == legacy n2_latent (bitwise)", cond_latent(n2o, 45.2), n2_latent(45.2), 1e-15);
    check("old set == legacy n2_psat (bitwise)", cond_psat(n2o, 38.0), n2_psat(38.0), 1e-15);
    { CondPropOpts o15 = newo; o15.liquidCp = 1500.0; const CondSpeciesProps n15 = condProps_make(COND_MODEL_N2, o15);
      check("cl=1500: L'(60) = cpv-1500", (cond_latent(n15, 60.05) - cond_latent(n15, 59.95))/0.1, cpv - 1500.0, 1e-8); }

    printf("== (b) CPG carrier EOS round trip (bracketed Newton) ==\n");
    const double gam = 1.4, cp_air = 1008.7, R_air = (gam - 1.0)*cp_air/gam, cv_air = cp_air/gam, Yw = 0.7671;
    int nbad = 0, ntot = 0, nnotok = 0;
    for (double T = 25.0; T <= 125.0; T += 5.0) for (double f : {0.0, 0.1, 0.5, 0.9, 0.99}) {
        const double g = f*Yw; const double e = (cv_air + g*RN2)*T - g*cond_latent(n2n, T);
        for (double Tg : {T, 40.0, 300.0, 1.5}) {   // 悪い初期推定も含む
            bool ok = false; const double Tr = cond_T_from_e_cpg(e, g, cv_air, RN2, Tg, n2n, &ok);
            const double eb = (cv_air + g*RN2)*Tr - g*cond_latent(n2n, Tr);
            ++ntot; if (!ok) ++nnotok;
            const double tolE = 1e-9*std::fabs(e) + 0.05, tolT = tolE/(cv_air + g*RN2);   // 解法の停止条件と同じ (e 残差 → T 誤差)
            if (!ok || std::fabs(Tr - T) > 2.0*tolT || std::fabs(eb - e) > tolE) { ++nbad; if (nbad <= 5) printf("      FAIL T=%.1f g=%.4f Tg=%.1f -> Tr=%.6f ok=%d de=%.3e\n", T, g, Tg, Tr, (int)ok, eb - e); }
        }
    }
    printf("      %d/%d inversions exact (ok=false: %d)\n", ntot - nbad, ntot, nnotok); checkb("EOS round trip all exact", nbad == 0);
    { // codex counterexample (pure-like R=RN2): g=0.75, T=122 K
      const double g = 0.75, T = 122.0, cv = 742.0; const double e = (cv + g*RN2)*T - g*cond_latent(n2n, T);
      bool ok = false; const double Tr = cond_T_from_e_cpg(e, g, cv, RN2, 300.0, n2n, &ok); check("codex case g=0.75 T=122 K converges", Tr, T, 1e-6); checkb("  ok flag true", ok); }
    printf("== (c) pure limit / carrier definitions ==\n");
    { const double T = 40.0, g = 0.1, rho = 0.5; const double e_pure = (cv_air + g*R_air)*T - g*cond_latent(n2n, T);
      bool ok; const double Tr = cond_T_from_e_cpg(e_pure, g, cv_air, R_air, 45.0, n2n, &ok); check("pure (Rw=R) inversion", Tr, T, 1e-10);
      const double p_c = rho*T*(R_air - g*RN2), pv = rho*(Yw - g)*RN2*T; printf("      carrier p=%.3f Pa, p_v=%.3f Pa (p_v/p=%.4f)\n", p_c, pv, pv/p_c);
      checkb("p_v < p and R_eff > 0", pv < p_c && (R_air - g*RN2) > 0.0); }
    printf("== (d) SLAU face enthalpy consistency ==\n");
    { const double g = 0.1, Tcell = 40.0, rho_f = 0.6, Tf = 45.0; const double Reff = R_air - g*RN2; const double pf = rho_f*Reff*Tf;
      // 新: T_f = p_f/(ρ_f R_eff), h_f = cp T_f − g L(T_f)
      const double Tf_rec = pf/(rho_f*Reff); const double hf = cp_air*Tf_rec - g*cond_latent(n2n, Tf_rec);
      const double e_eos = (cv_air + g*RN2)*Tf - g*cond_latent(n2n, Tf); const double h_eos = e_eos + pf/rho_f;
      check("h_f (face state) == e+p/rho at same state", hf, h_eos, 1e-12);
      const double h_old = gam*pf/((gam - 1.0)*rho_f) + g*(cp_air*Tcell - cond_latent(n2n, Tcell));
      printf("      old mixed-temperature form: h_old=%.2f vs EOS %.2f (diff %.1f J/kg; the inconsistency codex found)\n", h_old, h_eos, h_old - h_eos); }
    printf("== (e) slip energies ==\n");
    { const double rho = 0.5, u = 300.0, v = 50.0, w = 0.0, e = 30000.0; const double nx = 0.0, ny = 1.0, nz = 0.0; const double Un = u*nx + v*ny + w*nz;
      const double roe_i = rho*(e + 0.5*(u*u + v*v + w*w)); const double ug = u - 2*Un*nx, vg = v - 2*Un*ny, wg = w - 2*Un*nz; const double ub = u - Un*nx, vb = v - Un*ny, wb = w - Un*nz;
      check("ghost: reflected |u| gives roe_g == roe_i", rho*(e + 0.5*(ug*ug + vg*vg + wg*wg)), roe_i, 1e-14);
      check("bvar: roe_b == roe_i - 0.5 rho Un^2", rho*(e + 0.5*(ub*ub + vb*vb + wb*wb)), roe_i - 0.5*rho*Un*Un, 1e-14); }
    printf("== (f) air constants / kgas dispatch ==\n");
    { const double yN2 = 0.79, yO2 = 0.21, MN2 = 28.0134, MO2 = 31.9988; const double M = yN2*MN2 + yO2*MO2; const double Y = yN2*MN2/M;
      check("M_air (binary) = 28.850", M, 28.850, 1e-3); check("R_air = 288.19", 8314.462618/M, 288.19, 1e-4); check("Y_N2 = 0.7671", Y, 0.7671, 1e-3);
      check("cp_air = 3.5 R_air = 1008.7", 3.5*8314.462618/M, 1008.7, 1e-3); checkb("R_air - Yw R_N2 > 0", 8314.462618/M - Y*RN2 > 0);
      CondPropOpts oa = newo; oa.gasKgasModel = 1; const CondSpeciesProps na = condProps_make(COND_MODEL_N2, oa);
      check("cond_kgas model 0 == n2_kgas", cond_kgas(n2n, 60.0), n2_kgas(60.0), 1e-15); check("cond_kgas model 1 == air_kgas", cond_kgas(na, 60.0), air_kgas(60.0), 1e-15);
      printf("      k_gas(60 K): N2 %.4e, air %.4e W/m/K\n", n2_kgas(60.0), air_kgas(60.0)); }
    printf("%s (%d failures)\n", nfail ? "FAILED" : "ALL PASS", nfail); return nfail ? 1 : 0;
}
