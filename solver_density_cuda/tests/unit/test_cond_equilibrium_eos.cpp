// =============================================================================
// test_cond_equilibrium_eos.cpp  (host-only 単体検証)
//   plans/accepted/condensation-equilibrium-eos.md §6: EOS 拘束形平衡 (T,g) 同時反転の 0-D 検証。
//   case/44 旧条件 (H2O 2.87 wt% in 擬似種 MIXDRY) の膨張状態を模した (ρ, e_in) で
//     (1) 未飽和状態 → g=0、T が単相反転と一致
//     (2) 過飽和状態 → g>0、p_v(T,g)=p_sat(T) (S=1) と e_mix(T,g)=e_in を同時に満たす
//     (3) 緩和形 cond_equilibrium_delta を反復した固定点と g が一致 (同じ平衡)
//     (4) g_guess の初期値によらず同じ解、単調性 (e_in を下げると g 増)
//   build: nvcc -x cu --expt-relaxed-constexpr -I solver_density_cuda -o test_cond_eq_eos tests/unit/test_cond_equilibrium_eos.cpp
//   (plan #10 以降: H2O の潜熱に共通データの気液ペアが要るので、-I <solver_density_cuda> -I <埋め込みヘッダの生成先> と
//    solver_density_cuda/input/speciesDB.cpp solver_density_cuda/input/speciesTransportDB.cpp -lyaml-cpp を足す;
//    生成は tests/unit/cond_latent_test_helper.cuh 冒頭)
// =============================================================================
#include <cstdio>
#include <cmath>
#include "../../cuda_forge/condensationEOS_d.cuh"
#include "tests/unit/cond_latent_test_helper.cuh"   // H2O 潜熱の気液ペア (plan thermophysics-solver-owned-species-db #10)

static int nfail = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::printf("  FAIL: %s\n", msg); ++nfail; } else { std::printf("  ok  : %s\n", msg); } } while(0)

// NASA-9 (CEA) 係数: H2O と、MIXDRY 代わりに N2 (組成の詳細は本テストの目的に無関係)。
static SpeciesThermo makeN2() {
    SpeciesThermo s{}; s.MW = 0.0280134; s.sigma_LJ = 3.621; s.eps_kB = 97.53; s.Tlo = 200; s.Tmid = 1000; s.Thi = 6000;
    const double lo[9] = {2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double hi[9] = {5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    for (int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; } return s;
}
static SpeciesThermo makeH2O() {
    SpeciesThermo s{}; s.MW = 0.0180153; s.sigma_LJ = 2.605; s.eps_kB = 572.4; s.Tlo = 200; s.Tmid = 1000; s.Thi = 6000;
    const double lo[9] = {-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double hi[9] = {1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    for (int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; } return s;
}

int main()
{
    SpeciesThermo sp[2] = { makeN2(), makeH2O() };
    const int nSp = 2;
    const double Y[2] = { 0.9713353, 0.0286647 };   // 旧条件の H2O 質量分率
    const double Yw = Y[1];
    const CondSpeciesProps cprops = cond_test_props_H2O(false);
    const double Rw = cprops.R;
    const double Rmix = thermo_R_mix(sp, nSp, Y);
    const double Tmin = 50.0, Tmax = 6000.0;

    auto e_gas = [&](double T){ double cp,h; thermo_cph_mix(sp,nSp,Y,T,&cp,&h); return h - Rmix*T; };
    auto e_mix = [&](double T,double g){ return e_gas(T) + g*(Rw*T - cond_latent(cprops,T)); };

    // --- (1) 未飽和: T=282 K, P=5040 Pa (新条件出口相当) ---
    {
        const double T = 282.0, P = 5040.0; const double rho = P/(Rmix*T);
        const double e_in = e_gas(T);
        double Tn; const double g = cond_equilibrium_Tg_carrier(sp,nSp,Y,e_in,rho,Yw,Rw,cprops,300.0,0.0,Tmin,Tmax,&Tn);
        std::printf("(1) T=%g rho=%g -> g=%.3e T=%.4f  S=%.3f\n", T, rho, g, Tn, rho*Yw*Rw*Tn/cond_psat(cprops,Tn));
        CHECK(g == 0.0, "(1) unsaturated -> g=0");
        CHECK(std::fabs(Tn-T) < 1.0e-2, "(1) T equals single-phase inversion");
    }
    // --- (2) 過飽和: T=244 K, P=4640 Pa (旧条件 dry の出口軸 = S 3.9) を e 一定で平衡化 ---
    double g_ref = 0.0, T_ref = 0.0, rho_ref = 0.0, e_ref = 0.0;
    {
        const double T = 244.0, P = 4640.0; const double rho = P/(T*(Rmix));   // 乾き状態 (g=0) の密度
        const double e_in = e_gas(T);
        double Tn; const double g = cond_equilibrium_Tg_carrier(sp,nSp,Y,e_in,rho,Yw,Rw,cprops,T,0.0,Tmin,Tmax,&Tn);
        const double pv = rho*(Yw-g)*Rw*Tn, ps = cond_psat(cprops,Tn);
        std::printf("(2) dry T=%g S0=%.2f -> g=%.5f (%.1f %% of Yw) T=%.3f  S=%.6f  e_mix-e_in=%.3e J/kg\n",
                    T, rho*Yw*Rw*T/cond_psat(cprops,T), g, 100*g/Yw, Tn, pv/ps, e_mix(Tn,g)-e_in);
        CHECK(g > 0.0 && g < Yw, "(2) saturated -> 0<g<Yw");
        CHECK(std::fabs(pv/ps - 1.0) < 1.0e-6, "(2) S=1 (p_v=p_sat)");
        CHECK(std::fabs(e_mix(Tn,g)-e_in) < 1.0e-3*std::fabs(e_in) + 50.0, "(2) energy consistent");
        CHECK(Tn > T && Tn < T + 40.0, "(2) latent heat warms the gas by O(10 K)");
        g_ref = g; T_ref = Tn; rho_ref = rho; e_ref = e_in;
    }
    // --- (3) 緩和形の固定点と一致: cond_equilibrium_delta を g に足し込みながら T を再反転する反復 ---
    {
        double g = 0.0, T = 244.0;
        for (int it = 0; it < 200; ++it) {
            double cp,h; thermo_cph_mix(sp,nSp,Y,T,&cp,&h); const double cvg = cp - Rmix;
            const double D = cond_equilibrium_delta(cprops, rho_ref, T, g, Yw - g, Rw, cvg);
            g += D; if (g < 0) g = 0; if (g > Yw) g = Yw;
            T = cond_T_from_e_carrier(sp,nSp,Y,e_ref,g,Rw,cprops,T,Tmin,Tmax);
            if (std::fabs(D) < 1.0e-12) break;
        }
        std::printf("(3) relaxation fixed point: g=%.5f T=%.3f  | eos: g=%.5f T=%.3f\n", g, T, g_ref, T_ref);
        CHECK(std::fabs(g - g_ref) < 1.0e-5*Yw, "(3) g matches relaxation-form fixed point");
        CHECK(std::fabs(T - T_ref) < 1.0e-2, "(3) T matches relaxation-form fixed point");
    }
    // --- (4) 初期値非依存 + 単調性 ---
    {
        double T1,T2,T3;
        const double ga = cond_equilibrium_Tg_carrier(sp,nSp,Y,e_ref,rho_ref,Yw,Rw,cprops,244.0,0.5*Yw,Tmin,Tmax,&T1);
        const double gb = cond_equilibrium_Tg_carrier(sp,nSp,Y,e_ref,rho_ref,Yw,Rw,cprops,300.0,0.999*Yw,Tmin,Tmax,&T2);
        const double gc = cond_equilibrium_Tg_carrier(sp,nSp,Y,e_ref-5.0e3,rho_ref,Yw,Rw,cprops,244.0,0.0,Tmin,Tmax,&T3);
        std::printf("(4) g(guess .5Yw)=%.6f g(guess .999Yw)=%.6f | e_in-5kJ: g=%.5f T=%.3f\n", ga, gb, gc, T3);
        CHECK(std::fabs(ga-g_ref) < 1.0e-7 && std::fabs(gb-g_ref) < 1.0e-7, "(4) independent of g_guess");
        CHECK(gc > g_ref, "(4) lower e_in (colder) -> more condensate");
    }
    // --- (5) pure CPG (N2 相当の簡易チェック): 飽和線上 ---
    {
        const CondSpeciesProps cn = condProps_N2();
        const double cv = 743.0, R = 296.8; const double T = 60.0, rho = 1.0;    // p=17.8 kPa > p_sat(60 K)≈6.7 kPa (過飽和)
        double Tn; const double g = cond_equilibrium_Tg_pure_cpg(cv*T, rho, cv, R, cn, T, 0.0, &Tn);
        const double S = rho*(1-g)*R*Tn/cond_psat(cn,Tn);
        std::printf("(5) pure CPG N2: g=%.4f T=%.3f S=%.6f\n", g, Tn, S);
        CHECK(g > 0.0 && std::fabs(S-1.0) < 1.0e-6, "(5) pure CPG saturated -> S=1");
    }
    std::printf(nfail ? "FAILED (%d)\n" : "ALL PASSED\n", nfail);
    return nfail ? 1 : 0;
}
