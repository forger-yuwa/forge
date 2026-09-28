// =============================================================================
// test_cond_evaporation.cpp  (host-only 単体検証)
//   plans/accepted/condensation-evaporation.md §6-1: 0-D 蒸発の検証。
//   S<1 の一様場に monodisperse 液滴雲 (r0, g0) を置き、cond_evap_source の λ 更新を
//   時間積分して (1) HK 解析解 r(t)=r0+ṙt, g(t)=g0(r/r0)^3 との一致、(2) 実現可能性
//   (Q1^2<=Q0Q2, 全モーメント同時消滅)、(3) 消滅 r30<2rmin、(4) λ_min/dT_max 律速、
//   (5) Kelvin on/off の質量収支差、を確認する。
//   build: nvcc -x cu --expt-relaxed-constexpr -o test_cond_evap tests/unit/test_cond_evaporation.cpp
//   (plan #10 以降: H2O の潜熱に共通データの気液ペアが要るので、-I <solver_density_cuda> -I <埋め込みヘッダの生成先> と
//    solver_density_cuda/input/speciesDB.cpp solver_density_cuda/input/speciesTransportDB.cpp -lyaml-cpp を足す;
//    生成は tests/unit/cond_latent_test_helper.cuh 冒頭)
//          (または g++ -x c++ -DCUDA_HOST_ONLY ... は __host__ __device__ 修飾のため nvcc 推奨)
// =============================================================================
#include <cstdio>
#include <cmath>
#include "../../cuda_forge/condensationSource_d.cuh"
#include "tests/unit/cond_latent_test_helper.cuh"   // H2O 潜熱の気液ペア (plan thermophysics-solver-owned-species-db #10)

static int nfail = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::printf("  FAIL: %s\n", msg); ++nfail; } else { std::printf("  ok  : %s\n", msg); } } while(0)

int main()
{
    const CondSpeciesProps cp = cond_test_props_H2O(false);
    const double Rw = cp.R;
    // 場: T=290K, S=0.7 (run_0065 圧縮帯相当), キャリア ρ=0.28 (P≈23 kPa)
    const double T = 290.0, S = 0.7, rod = 0.28;
    const double psat = cond_psat(cp, T);
    const double pv0  = S*psat;
    const double rho_l = cond_rho_cond(cp, T);
    const double cvg = 900.0;
    // 液滴雲: r0=300nm, g0=1.2% → Q0 = g/(4/3 π ρ_l r0^3) [1/kg]
    const double r0 = 3.0e-7, g0 = 0.012;
    const double Q0 = g0/((4.0/3.0)*COND_PI*rho_l*r0*r0*r0);
    std::printf("T=%g psat=%g pv=%g r0=%g g0=%g Q0=%.3e\n", T, psat, pv0, r0, g0, Q0);

    // --- (1) HK 解析解との一致 (小 dt, 律速なし; pv 固定=無限蒸気槽) ---
    {
        const double drdt_ref = (1.0/rho_l)*(pv0 - psat)/std::sqrt(2.0*COND_PI*Rw*T);
        double q0 = rod*Q0, q1 = q0*r0, q2 = q0*r0*r0, g = g0;
        const double dt = 1.0e-6;   // |ṙ|dt/r0 ≈ 2e-3
        double t = 0.0;
        for (int n = 0; n < 200; ++n) {
            double S0,S1,S2,Sg,r30,dr;
            cond_evap_source(cp, T, pv0, rod, g, q0,q1,q2, dt, 1e-9, 0.5, 1e9, 1e9, cvg, 0, -1.0, 3.18, 0,
                             &S0,&S1,&S2,&Sg,&r30,&dr);
            q0 += S0*dt; q1 += S1*dt; q2 += S2*dt; g += Sg*dt/rod; t += dt;
        }
        const double r_an = r0 + drdt_ref*t;
        const double g_an = g0*std::pow(r_an/r0, 3);
        const double r30 = std::cbrt(g/((4.0/3.0)*COND_PI*rho_l*q0/rod));
        const double r10 = q1/q0;
        std::printf("(1) t=%.2e drdt_ref=%.3e  r30=%.4e r_an=%.4e (%.2e)  g=%.5e g_an=%.5e (%.2e)  r10/r30=%.5f\n",
                    t, drdt_ref, r30, r_an, r30/r_an-1, g, g_an, g/g_an-1, r10/r30);
        CHECK(std::fabs(r30/r_an-1) < 5e-3, "r30 follows HK analytic r0+drdt*t (<0.5%)");
        CHECK(std::fabs(g/g_an-1) < 2e-2,   "g follows g0 (r/r0)^3 (<2%; explicit Euler)");
        CHECK(std::fabs(r10/r30-1) < 1e-3,  "monodisperse consistency r10=r30 preserved");
        CHECK(q1*q1 <= q0*q2*(1.0+1e-9),    "realizability Q1^2<=Q0Q2");
    }

    // --- (2)(3) 大 dt: λ_min 半減律速 → 一括消滅 (r30<2rmin) ---
    {
        double q0 = rod*Q0, q1 = q0*r0, q2 = q0*r0*r0, g = g0;
        const double dt = 1.0e-2;   // |ṙ|dt >> r0 → λ_min 律速
        int n = 0, nremove = -1;
        for (; n < 40; ++n) {
            double S0,S1,S2,Sg,r30,dr;
            const double lam = cond_evap_source(cp, T, pv0, rod, g, q0,q1,q2, dt, 1e-9, 0.5, 1e9, 1e9, cvg, 0, -1.0, 3.18, 0,
                             &S0,&S1,&S2,&Sg,&r30,&dr);
            if (n == 0) { std::printf("(2) step0 lam=%.3f (expect 0.5)\n", lam); CHECK(std::fabs(lam-0.5)<1e-12, "lam clamped to lam_min=0.5"); }
            q0 += S0*dt; q1 += S1*dt; q2 += S2*dt; g += Sg*dt/rod;
            if (lam <= 0.0) { nremove = n; break; }
        }
        std::printf("(3) removal at step %d (expect ~log2(300nm/2nm)=7..8), after: q0=%.2e q1=%.2e q2=%.2e g=%.2e\n", nremove, q0,q1,q2,g);
        CHECK(nremove >= 7 && nremove <= 9, "complete removal after ~8 halvings");
        CHECK(std::fabs(q0) < 1e-12*rod*Q0 && std::fabs(g) < 1e-12*g0, "all moments zeroed together on removal");
    }

    // --- (4) dT_max 律速: 1 step の潜熱冷却 <= dT_max ---
    {
        double q0 = rod*Q0, q1 = q0*r0, q2 = q0*r0*r0, g = g0;
        const double dt = 1.0e-2, dT_max = 1.0;
        double S0,S1,S2,Sg,r30,dr;
        cond_evap_source(cp, T, pv0, rod, g, q0,q1,q2, dt, 1e-9, 0.5, 5e-3, dT_max, cvg, 0, -1.0, 3.18, 0,
                         &S0,&S1,&S2,&Sg,&r30,&dr);
        const double dg = -Sg*dt/rod;
        const double dT = dg*cond_latent(cp,T)/cvg;
        std::printf("(4) dg=%.3e dT=%.3f K (limit %.1f)\n", dg, dT, dT_max);
        CHECK(dT <= dT_max*(1+1e-9), "latent cooling per step limited to dT_max");
    }

    // --- (5) Kelvin on/off: 蒸発総量は同じ (消滅まで積分)、時間だけ変わる ---
    {
        for (int kelvin = 0; kelvin <= 1; ++kelvin) {
            double q0 = rod*Q0, q1 = q0*r0, q2 = q0*r0*r0, g = g0;
            const double dt = 2.0e-6; double t = 0.0; double g_at_10nm = -1;
            for (int n = 0; n < 2000000; ++n) {
                double S0,S1,S2,Sg,r30,dr;
                const double lam = cond_evap_source(cp, T, pv0, rod, g, q0,q1,q2, dt, 1e-9, 0.5, 1e9, 1e9, cvg, 0, -1.0, 3.18, kelvin,
                                 &S0,&S1,&S2,&Sg,&r30,&dr);
                if (g_at_10nm < 0 && r30 > 0 && r30 < 1.0e-8) g_at_10nm = g;
                q0 += S0*dt; q1 += S1*dt; q2 += S2*dt; g += Sg*dt/rod; t += dt;
                if (lam <= 0.0) break;
            }
            std::printf("(5) kelvin=%d: t_evap=%.4e s, g remaining when r30<10nm = %.2e (=%.1e of g0)\n", kelvin, t, g_at_10nm, g_at_10nm/g0);
            CHECK(g_at_10nm/g0 < 1e-4, "mass below 10 nm is <1e-4 of g0 (Kelvin irrelevant to mass budget)");
        }
    }

    // --- (6) S>1 では蒸発分岐に入らない ---
    {
        double S0,S1,S2,Sg,r30,dr;
        const double lam = cond_evap_source(cp, T, 1.5*psat, rod, g0, rod*Q0, rod*Q0*r0, rod*Q0*r0*r0, 1e-6, 1e-9, 0.5, 1e9, 1e9, cvg, 0, -1.0, 3.18, 0,
                                            &S0,&S1,&S2,&Sg,&r30,&dr);
        CHECK(lam == 1.0 && Sg == 0.0, "S>1: no evaporation source (lam=1)");
    }
    // --- (7) N2 Goodheart 経路も負・有限 ---
    {
        const CondSpeciesProps n2 = condProps_N2();
        const double Tn = 60.0; const double ps = cond_psat(n2, Tn);
        const double dr0 = cond_evap_rate(n2, Tn, 0.5*ps, 1e-7, 0, 0.5*ps, 3.18, 0);
        const double dr1 = cond_evap_rate(n2, Tn, 0.5*ps, 1e-7, 1, 0.5*ps, 3.18, 0);
        std::printf("(7) N2 S=0.5 T=60K r=100nm: Goodheart drdt=%.3e Gyarmathy drdt=%.3e\n", dr0, dr1);
        CHECK(dr0 < 0 && std::isfinite(dr0) && dr1 < 0 && std::isfinite(dr1), "N2 Goodheart/Gyarmathy evaporation rate negative & finite");
    }

    std::printf(nfail ? "RESULT: %d FAIL\n" : "RESULT: ALL PASS\n", nfail);
    return nfail ? 1 : 0;
}
