// =============================================================================
// test_cond_air.cpp  (host 単体検証; nvcc -x cu --expt-relaxed-constexpr -I. -o test_cond_air tests/unit/test_cond_air.cpp)
//   plans/accepted/condensation-air.md §5 (7):
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
#include <limits>
#include <algorithm>
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
    printf("== (g) SLAU CPG two-phase face enthalpy in float32 (kernel inputs are flow_float) ==\n");
    {
        // 実カーネル (convectiveFlux_slau_d.inc.cuh) と同じ経路: P_L, ro_L, g, velocity2 は float に丸めてから double で cond_face_h_cpg。
        // 判定: float 入力による誤差は入力丸めの伝播 (~1e-7 相対) に収まり、旧 (セル温度混在) の ~1 kJ/kg ずれは出ない。
        int nb = 0, nt = 0; double worst = 0.0;
        for (double T = 30.0; T <= 120.0; T += 10.0) for (double f : {0.0, 0.3, 0.9}) for (double P : {200.0, 2.0e3, 5.0e5}) {
            const double g = f*Yw, Reff = R_air - g*RN2, rho = P/(Reff*T), u2 = 2.0*250.0*250.0;
            const double h_ref = cp_air*T - g*cond_latent(n2n, T) + 0.5*u2;                          // 面状態から解析的に
            const float Pf = (float)P, rf = (float)rho, gf = (float)g, u2f = (float)u2;
            const double h_f32 = cond_face_h_cpg(n2n, cp_air, R_air, RN2, (double)gf, (double)Pf, (double)rf, 0.5*(double)u2f);
            const double tolh = 4.0e-6*std::fabs(h_ref) + 1.0e-6*(cp_air*T + g*cond_latent(n2n, T));   // float 入力 (T_f=p/(ρR) の ~2 ulp) の伝播
            ++nt; const double d = std::fabs(h_f32 - h_ref); if (d > worst) worst = d;
            if (d > tolh) { ++nb; if (nb <= 5) printf("      FAIL T=%.0f g=%.3f P=%.0f: h32=%.4f href=%.4f (d=%.3e J/kg)\n", T, g, P, h_f32, h_ref, d); }
        }
        printf("      %d/%d states within float-input tolerance (worst |dh| = %.3e J/kg)\n", nt - nb, nt, worst);
        checkb("face enthalpy float32 inputs == analytic face state (no ~1 kJ/kg cell-T mixing error)", nb == 0);
        // pure 極限 (R_w=R_gas, g=0): 単相 γp/((γ−1)ρ)+ek と一致
        const double P = 800.0, rho = P/(R_air*40.0);
        check("g=0 pure limit == gamma p/((gamma-1) rho) + ek", cond_face_h_cpg(n2n, cp_air, R_air, R_air, 0.0, P, rho, 100.0), gam*P/((gam - 1.0)*rho) + 100.0, 1e-12);
        // Y_w 定数キャリアで g=Y_w (全 N2 凝縮) でも R_eff = R_O2 分が残り有限
        const double hfull = cond_face_h_cpg(n2n, cp_air, R_air, RN2, Yw, P, P/((R_air - Yw*RN2)*40.0), 0.0);
        checkb("g=Y_w: finite and equals cp T - Yw L(T)", std::isfinite(hfull) && std::fabs(hfull - (cp_air*40.0 - Yw*cond_latent(n2n, 40.0))) < 1e-6*std::fabs(hfull));
    }
    printf("== (h) realizability clamp: dust check uses vapor partial pressure (codex 2026-09-13 M2) ==\n");
    {
        // 状態: 空気 CPG carrier, T=40 K, ρ で p=1.2 p_sat(N2) だが p_v = (Y_w−g) ρ R_N2 T < p_sat → 「未飽和」で塵を掃除できるべき。
        const double T = 40.0, ps = cond_psat(n2n, T), g = 0.5*Yw, Reff = R_air - g*RN2, rod = 1.2*ps/(Reff*T), P = rod*Reff*T;
        const double pv_cpg  = cond_clamp_vapor_pressure(rod, g, -1.0, Yw, RN2, T, P);
        const double pv_pure = cond_clamp_vapor_pressure(rod, g, -1.0, -1.0, RN2, T, P);
        const double pv_tp   = cond_clamp_vapor_pressure(rod, g, Yw, -1.0, RN2, T, P);
        printf("      P/psat=%.3f  pv_cpg/psat=%.3f  pv_tp/psat=%.3f\n", P/ps, pv_cpg/ps, pv_tp/ps);
        check("CPG carrier pv = rho (Yw-g) Rw T", pv_cpg, rod*(Yw - g)*RN2*T, 1e-14);
        check("TP carrier (roY_w) gives the same pv", pv_tp, pv_cpg, 1e-14);
        check("pure: pv = total P", pv_pure, P, 1e-14);
        checkb("CPG carrier: P > psat (old rule keeps dust) but pv < psat (new rule removes it)", P > ps && pv_cpg < ps);
        check("g -> Yw: pv -> 0 (clamped, not negative)", cond_clamp_vapor_pressure(rod, Yw*1.01, -1.0, Yw, RN2, T, P), 0.0, 1e-300);
    }
    printf("== (i) inversion rejects abnormal inputs (codex 2026-09-13 result-2 M1) ==\n");
    {
        const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
        struct Case { const char* name; double e, g, Tg; } cases[] = {
            {"e=+Inf", +inf, 0.1, 40.0}, {"e=-Inf", -inf, 0.1, 40.0}, {"e=NaN", nan, 0.1, 40.0}, {"g=NaN", 1.0e4, nan, 40.0},
            {"T_guess=NaN", 1.0e4, 0.1, nan}, {"e beyond T=6000 K (unreachable)", (cv_air + 0.1*RN2)*7000.0, 0.1, 40.0},
            {"e below T=1 K (unreachable)", (cv_air + 0.1*RN2)*0.2 - 0.1*cond_latent(n2n, 1.0), 0.1, 40.0} };
        int nb = 0;
        for (const Case& c : cases) {
            bool ok = true; const double Tr = cond_T_from_e_cpg(c.e, c.g, cv_air, RN2, c.Tg, n2n, &ok);
            const bool good = (!ok) && std::isfinite(Tr);
            if (!good) ++nb;
            printf("      %-36s ok=%d T=%g -> %s\n", c.name, (int)ok, Tr, good ? "rejected" : "ACCEPTED (BUG)");
        }
        checkb("all abnormal inputs return ok=false with finite T", nb == 0);
        bool ok = false; const double T0 = 40.0, e0 = (cv_air + 0.1*RN2)*T0 - 0.1*cond_latent(n2n, T0);
        check("normal input still converges (ok=true)", cond_T_from_e_cpg(e0, 0.1, cv_air, RN2, 300.0, n2n, &ok), T0, 1e-9); checkb("  ok flag true", ok);
    }
    printf("== (j) face enthalpy near depletion (R_eff < 1) is the EOS state, no floor (codex 2026-09-13 result-2 M3) ==\n");
    {
        // 受付条件 (R_air − Y_w R_w > 0) だけを満たす Y_w=g=0.999 (R_eff=0.26): 面温度は EOS 温度に一致し h_f = e + p/ρ + ek。
        // 気相定数 2 組: 空気 (R 288.19; Y_w≤0.97 まで受付) と codex の反例 (γ 1.4, cp 1038.67 → R_gas=R_N2: Y_w=g=0.999 で R_eff=0.297)
        int nb = 0, nlow = 0; double worst = 0.0, reff_min = 1e300;
        for (int k = 0; k < 2; ++k) {
            const double Rg = (k == 0) ? R_air : RN2, cpg = 3.5*Rg, cvg = 2.5*Rg;
            for (double Ywx : {0.7671, 0.95, 0.99, 0.999}) for (double f : {0.9, 0.999, 1.0}) for (double T : {30.0, 40.0, 60.0}) {
                const double g = f*Ywx, Reff = Rg - g*RN2; if (!(Reff > 0.0)) continue;   // 受付条件 R_gas − Y_w R_w > 0 の範囲だけ
                reff_min = std::min(reff_min, Reff); if (Reff < 1.0) ++nlow;
                const double rho = 1.0, P = rho*Reff*T, e = (cvg + g*RN2)*T - g*cond_latent(n2n, T);
                const double h = cond_face_h_cpg(n2n, cpg, Rg, RN2, g, P, rho, 0.0), href = e + P/rho;
                const double d = std::fabs(h - href)/std::fabs(href); worst = std::max(worst, d);
                if (d > 1e-12) { ++nb; printf("      FAIL Rgas=%.2f Yw=%.4f g=%.5f T=%.0f Reff=%.4f h=%.4f href=%.4f\n", Rg, Ywx, g, T, Reff, h, href); }
            }
        }
        printf("      worst rel |h_f - (e+p/rho)| = %.2e over accepted states; min R_eff hit = %.4f (%d states with R_eff < 1)\n", worst, reff_min, nlow);
        checkb("sweep reaches R_eff < 1 (old floor would have altered these)", nlow > 0);
        checkb("h_f == e + p/rho for all accepted states incl. R_eff < 1", nb == 0);
        // 異常 (g_f > Y_w → R_eff <= 0, または非有限): 乾き面へ退避して有限
        const double hbad = cond_face_h_cpg(n2n, cp_air, R_air, RN2, 1.2, 800.0, 0.1, 0.0);
        checkb("R_eff <= 0 falls back to the dry face (finite, = cp T_dry)", std::isfinite(hbad) && std::fabs(hbad - cp_air*(800.0/(0.1*R_air))) < 1e-9*hbad);
        checkb("non-finite g_f falls back to the dry face", std::isfinite(cond_face_h_cpg(n2n, cp_air, R_air, RN2, std::numeric_limits<double>::quiet_NaN(), 800.0, 0.1, 0.0)));
    }
    // (g) 蒸気 c_p,v の出どころ (2026-09-14): Kirchhoff の傾き L' = c_p,v − c_l は **凝縮する蒸気** の c_p を使う。
    //     pure-condensible CPG だけ config の physProp.cp を渡す (CondPropOpts::gasCp)。carrier では渡さない。
    printf("== (g) vapour c_p source for the Kirchhoff slope (below 70 K) ==\n");
    {
        CondPropOpts o{1, 1, 2000.0, 0, 1.0, -1.0};          // gasCp 未指定 = 内蔵 1038.8
        CondPropOpts og = o; og.gasCp = 1050.0;               // config 由来で上書き
        const CondSpeciesProps a = condProps_make(COND_MODEL_N2, o), b = condProps_make(COND_MODEL_N2, og);
        check("gasCp unset keeps the built-in N2 vapour cp", a.cp, 1038.8, 1e-12);
        check("gasCp set overrides the vapour cp", b.cp, 1050.0, 1e-12);
        for (double T : {80.0, 100.0}) check("T>=70 K is the polynomial (cp_v irrelevant)", cond_latent(a, T), cond_latent(b, T), 1e-14);
        const double T1 = 45.0, dcp = 1050.0 - 1038.8;
        check("below 70 K the slope moves by (cp_v_new - cp_v_old)", cond_latent(b, T1) - cond_latent(a, T1), dcp*(T1 - 70.0), 1e-10);
        // 液比熱は全域で正 (c_l = c_p,v - L' > 0)
        int nneg = 0;
        for (double T = 25.0; T <= 125.0; T += 0.25) {
            const double dL = (cond_latent(a, T + 0.05) - cond_latent(a, T - 0.05))/0.1;
            if (a.cp - dL <= 0.0) ++nneg;
        }
        checkb("implied liquid cp = cp_v - dL/dT stays positive over 25-125 K", nneg == 0);
    }

    // (h) H2O: 気相 h の有効域外は種 DB と同じ「端点 cp 一定の線形外挿」(2026-09-14)。
    //     生の多項式評価だと dL/dT が 200 K 未満で温度依存になり、含意される液比熱が 4228 からずれる。
    printf("== (h) H2O latent heat below the gas fit range (200 K) ==\n");
    {
        const CondSpeciesProps w = condProps_H2O();
        const double cpl_ref = 4228.268;   // CEA H2O(L) の 273.15 K 解析 cp
        double worst = 0.0;
        for (double T = 130.0; T <= 265.0; T += 5.0) {
            const double dL = (cond_latent(w, T + 0.05) - cond_latent(w, T - 0.05))/0.1;
            const double cpv = 1851.2;     // NASA-9 H2O(g) の 200 K 値 (200 K 未満は定 cp)
            if (T < 199.0) worst = std::max(worst, std::fabs((cpv - dL) - cpl_ref)/cpl_ref);
        }
        printf("      max |c_l,implied - 4228.27|/4228.27 over 130-199 K = %.2e\n", worst);
        checkb("implied liquid cp is the constant 4228 below 200 K (gas h is linear there)", worst < 2e-4);
        // 200 K の接続: 片側微分が一致する (定 cp 外挿は c_p(200) で接続するので C1)
        const double dm = (cond_latent(w, 199.9) - cond_latent(w, 199.8))/0.1;
        const double dp = (cond_latent(w, 200.2) - cond_latent(w, 200.1))/0.1;
        check("dL/dT is continuous across the 200 K join", dm, dp, 5e-5);
    }

    // (i) H2O 潜熱の 373 K 超 Watson 外挿 (2026-09-14, plans/active/condensation-h2o-latent-supercritical.md §6)。
    //     旧実装は h_l を 373.15 K でクランプしており L が温度とともに増加していた。
    printf("== (i) H2O latent heat above 373.15 K (Watson extrapolation) ==\n");
    {
        const CondSpeciesProps w = condProps_H2O();
        const double Tc = 647.096, T0 = 373.15, Tfr = 646.15;
        // 1. 373.15 K で値が連続 / 373.15-647 K で単調非増加
        // 接続点の両側を近づけて比べる (±1e-6 K だと傾き ~3000 J/(kg K) 分の 2.4e-9 が残り、連続性の判定にならない)
        check("value is continuous at the 373.15 K anchor",
              cond_latent(w, T0 - 1e-9), cond_latent(w, T0 + 1e-9), 1e-11);
        {
            int nup = 0; double prev = cond_latent(w, T0);
            for (double T = T0 + 0.05; T <= Tc; T += 0.05) {
                const double v = cond_latent(w, T); if (v > prev + 1e-6) ++nup; prev = v;
            }
            checkb("L is monotonically non-increasing over 373.15-647.1 K", nup == 0);
        }
        // 2. IAPWS 飽和線 (h_g - h_f) の慣用表値との比較。参照は 100/150/200/250/300/350 degC。
        {
            const double Tref[6] = {373.15, 423.15, 473.15, 523.15, 573.15, 623.15};
            const double Lref[6] = {2257e3, 2114e3, 1940e3, 1716e3, 1404e3, 893e3};
            double worst = 0.0, ss = 0.0;
            for (int i = 0; i < 6; ++i) {
                const double e = (cond_latent(w, Tref[i]) - Lref[i])/Lref[i];
                worst = std::max(worst, std::fabs(e)); ss += e*e;
            }
            const double rms = std::sqrt(ss/6.0);
            printf("      vs IAPWS 6 points: max |err| = %.2f %%, rms = %.2f %%\n", 100*worst, 100*rms);
            checkb("within 3 % of the IAPWS saturation line (max)", worst < 0.03);
            checkb("within 2 % of the IAPWS saturation line (rms)", rms < 0.02);
            // 指数 0.38 が 0.28-0.40 の中で rms 最小であること (計画 §4.2 の探索を回帰させる)
            const double L0 = cond_latent(w, T0);
            double best = 1e30; double bestn = 0.0;
            for (double n = 0.28; n <= 0.401; n += 0.01) {
                double s2 = 0.0;
                for (int i = 0; i < 6; ++i) {
                    const double Lw = L0*std::pow((Tc - Tref[i])/(Tc - T0), n);
                    const double e = (Lw - Lref[i])/Lref[i]; s2 += e*e;
                }
                if (s2 < best) { best = s2; bestn = n; }
            }
            printf("      exponent search over 0.28-0.40: best n = %.2f\n", bestn);
            checkb("Watson exponent 0.38 minimises the rms error", std::fabs(bestn - 0.38) < 0.005);
        }
        // 4. 凍結帯: T_fr 以上は一定、全域で L>0
        check("L is frozen above 646.15 K (647 K)",  cond_latent(w, 647.0),  cond_latent(w, Tfr), 1e-14);
        check("L is frozen above 646.15 K (1200 K)", cond_latent(w, 1200.0), cond_latent(w, Tfr), 1e-14);
        {
            int nz = 0;
            for (double T = 130.0; T <= 1200.0; T += 0.5) if (!(cond_latent(w, T) > 0.0)) ++nz;
            checkb("L stays strictly positive over 130-1200 K", nz == 0);
        }
        // 5. 成長則の有限性 (codex plan M1 の再現ケースを含む)
        {
            int nbad = 0; double worst_abs = 0.0;
            for (int gm = 0; gm <= 1; ++gm)
              for (double T : {400.0, 600.0, 646.0, 650.0, 700.0})
                for (double r : {1.0e-8, 1.0e-7}) {
                    const double pv = 1.0e5, pg = 1.0e5;   // p_v < p_sat (蒸発側) になる高温
                    const double d = cond_growth(w, T, pv, r, 0.0, gm, pg);
                    if (!std::isfinite(d)) { ++nbad; printf("      FAIL gm=%d T=%.0f r=%.0e drdt=%g\n", gm, T, r, d); }
                    else worst_abs = std::max(worst_abs, std::fabs(d));
                }
            printf("      max |dr/dt| over the hot sweep = %.4g m/s\n", worst_abs);
            checkb("growth/evaporation rate stays finite above 373 K (no 1/L^2 blow-up)", nbad == 0);
        }
        // 6. T_sat の往復 (codex plan M3: 根から離れた初期推定でも収束すること)
        {
            int nbad = 0; double worst = 0.0;
            for (double T : {300.0, 450.0, 550.0, 640.0})
                for (double guess : {250.0, 500.0, 900.0}) {
                    const double ps = cond_psat(w, T), Ts = cond_Tsat(w, ps, guess);
                    const double e = std::fabs(Ts - T);
                    worst = std::max(worst, e);
                    if (!(e < 1.0e-3)) { ++nbad; printf("      FAIL T=%.1f guess=%.0f -> Tsat=%.6f\n", T, guess, Ts); }
                }
            printf("      worst |Tsat(psat(T)) - T| = %.2e K\n", worst);
            checkb("Tsat inverts psat from far-away guesses (300-640 K)", nbad == 0);
        }
    }

    printf("%s (%d failures)\n", nfail ? "FAILED" : "ALL PASS", nfail); return nfail ? 1 : 0;
}
