// =============================================================================
// test_cond_latent_pair.cu — H2O 潜熱を種 DB の気液ペアで作る (plans/active/thermophysics-solver-owned-species-db.md §4.8, #10; §6 V7)
//   (a)  datum 不変性: thermoHrefTemp を 0 / 298.15 / 150 / 1234.5 K に変えても L(T) が 120–400 K の全点 (0.5 K 刻み) で相対 1e-12 以内
//        (反例 M6: 液相に同じ Δa7 を掛けないと崩れることも確かめる)
//   (a') 気相 h_v がソルバの種 DB 評価 (thermo_init_db と同じ datum を焼き込んだ thermo_h_mass) とビット一致 (50–6500 K, 区切り ±1e-9 K)
//   (b)  既知の気液差: L(298.15 K) = CEA thermo.inp の H2O (200–1000 K) と H2O(L) (273.15–373.15 K) の絶対エンタルピー差
//        (thermo.inp の係数を独立に書き写した参照; 相対 1e-12)。ヘッダの Hf 差 (−241826 + 285830 J/mol) との差は情報として出す
//   (c') float 表: 200 K (気相の区間端) とその float 両隣、200 K を含む表の小区間 [199.65, 200.40] K で L 相対 ≤ 2e-6 と
//        |ΔL'| ≤ 2e-4|L'| + 0.1 J/(kg K) を直接判定 (test_cond_float.cpp の接続点除外に頼らない)
//   (d)  g>0 の場で保存エネルギーを datum 変換 (0 → 298.15 K) した後、二相 EOS の反転 T・P・二相音速が保たれる
//   (新旧) 新 L と #10 以前の h2o_latent (試験内に移植) の差を 150–373.15 K で表にし、200 K 以上は丸め程度 (相対 1e-12)、
//        150 K で −465.920 J/kg (codex diagnose の検算値) を確かめる
//   (未設定) 気液ペアの無い H2O の cond_latent は NaN、device でも host と同じ値
//
// ビルド/実行 (共通データの埋め込みヘッダを先に生成; tests/unit/cond_latent_test_helper.cuh):
//   cmake -DIN=solver_density_cuda/data/species/forge_species_v1.yaml -DOUT=/tmp/forge_species_gen/forge_species_data.hpp
//       -P solver_density_cuda/cmake/embed_species_data.cmake
//   nvcc -std=c++17 --expt-relaxed-constexpr -I solver_density_cuda -I /tmp/forge_species_gen -o /tmp/test_cond_latent_pair
//       solver_density_cuda/tests/unit/test_cond_latent_pair.cu -x none solver_density_cuda/input/speciesDB.cpp
//       solver_density_cuda/input/speciesTransportDB.cpp -lyaml-cpp && /tmp/test_cond_latent_pair
// 規約: [PASS]/[FAIL] を出し、失敗があれば非ゼロ終了。
// =============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "tests/unit/cond_latent_test_helper.cuh"
#include "cuda_forge/condensationTables_d.cuh"
#include "cuda_forge/condensationEOS_d.cuh"

static int g_fail = 0;
static void check(bool ok, const char* what, double val, double tol)
{
    printf("  [%s] %-88s worst=%.3e tol=%.1e\n", ok ? "PASS" : "FAIL", what, val, tol);
    if (!ok) ++g_fail;
}

// #10 以前のソルバの h2o_latent (比較の基準としてだけ移植; ソルバからは撤去済み)
static double legacy_h(const double* a, double T)
{
    const double Rw = 8.314462618/0.0180153;
    const double hRT = -a[0]/(T*T) + a[1]*log(T)/T + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0 + a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[7]/T;
    return hRT*Rw*T;
}
static double legacy_h2o_latent_v0(double T)
{
    const double ag[8] = {-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03,
                          -7.342557370e-06, 4.955043490e-09,-1.336933246e-12,-3.303974310e+04};
    const double al[8] = { 1.326371304e+09,-2.448295388e+07, 1.879428776e+05,-7.678995050e+02,
                           1.761556813e+00,-2.151167128e-03, 1.092570813e-06, 1.101760476e+08};
    const double Tf = 273.15;
    double Tg = (T > 45.0) ? T : 45.0;
    if (Tg > 1000.0) Tg = 1000.0;
    const double hv = legacy_h(ag, Tg);
    double hl;
    if (Tg >= Tf) hl = legacy_h(al, (Tg < 373.15) ? Tg : 373.15);
    else hl = legacy_h(al, Tf) - (legacy_h(al, Tf + 0.5) - legacy_h(al, Tf - 0.5 + 1.0e-9))*(Tf - Tg);
    double L = hv - hl;
    if (L < 1.5e6) L = 1.5e6;
    if (L > 3.5e6) L = 3.5e6;
    return L;
}

// 種 DB の気相 (datum 前) に thermo_init_db と同じ演算で datum を焼き込む
static SpeciesThermo with_datum(SpeciesThermo s, double Tref)
{
    s.invMW = 1.0/s.MW;
    if (Tref > 0.0) { const double h_ref = thermo_h_molar(s, Tref); const double da7 = -h_ref/THERMO_RU; s.low[7] += da7; s.high[7] += da7; s.h_datum = h_ref; }
    return s;
}

__global__ void latent_dev(CondLatentRef r, const double* T, double* L, int n)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) L[i] = h2o_latent(r, T[i]);
}

int main()
{
    // ---------------------------------------------------------------- (a) datum 不変性
    printf("== (a) datum invariance of L = h_v - h_l ==\n");
    const CondLatentPair p0 = cond_test_latent_pair(0.0);
    {
        double worst = 0.0, worstU = 0.0;
        for (double Tref : {298.15, 150.0, 1234.5}) {
            const CondLatentPair p1 = cond_test_latent_pair(Tref);
            for (double T = 120.0; T <= 400.0 + 1e-9; T += 0.5) {
                worst = fmax(worst, fabs(h2o_latent_pair(p1, T) - h2o_latent_pair(p0, T))/h2o_latent_pair(p0, T));
                // クランプ前の差 (h_v − h_l) でも見る (クランプで隠れないように)
                const double u0 = cond_gas_h_mass(p0.gas, T) - cond_liquid_h_mass(p0, T);
                const double u1 = cond_gas_h_mass(p1.gas, T) - cond_liquid_h_mass(p1, T);
                worstU = fmax(worstU, fabs(u1 - u0)/fabs(u0));
            }
        }
        check(worst <= 1e-12, "L(T) Tref in {298.15,150,1234.5} vs 0, 120-400 K every 0.5 K (rel)", worst, 1e-12);
        check(worstU <= 1e-12, "unclamped h_v - h_l, same set (rel)", worstU, 1e-12);
        // 反例 (codex plan M6): 気相だけ datum を掛け、液相に同じ Δa7 を掛けない誤った置換 → L(300 K) が桁違いになる
        CondLatentPair bad = cond_test_latent_pair(298.15);
        bad.hShift = 0.0;
        const double ub = cond_gas_h_mass(bad.gas, 300.0) - cond_liquid_h_mass(bad, 300.0);
        printf("      counter-example (liquid without the shared shift): h_v - h_l at 300 K = %.4e J/kg (correct %.4e)\n",
               ub, cond_gas_h_mass(p0.gas, 300.0) - cond_liquid_h_mass(p0, 300.0));
        check(fabs(ub - (cond_gas_h_mass(p0.gas, 300.0) - cond_liquid_h_mass(p0, 300.0))) > 1e6,
              "counter-example is detected (difference > 1 MJ/kg)", fabs(ub), 1e6);
    }

    // ---------------------------------------------------------------- (a') 気相 = 種 DB の評価 (ビット一致)
    printf("== (a') gas h_v == species-DB thermo_h_mass (bitwise) ==\n");
    {
        ResolvedSpeciesDB db = speciesDB_resolve(std::vector<std::string>{"N2", "H2O"}, "");
        speciesDB_attachCondensed(db, "H2O(L)", "H2O", true);
        long nbad = 0, n = 0; double worst = 0.0;
        for (double Tref : {0.0, 298.15}) {
            const SpeciesThermo g = with_datum(db.species[1], Tref);
            const CondLatentPair p = cond_test_latent_pair(Tref);
            bool same = (g.MW == p.gas.MW && g.Tlo == p.gas.Tlo && g.Tmid == p.gas.Tmid && g.Thi == p.gas.Thi);
            for (int k = 0; k < 9; ++k) same = same && g.low[k] == p.gas.low[k] && g.high[k] == p.gas.high[k];
            if (!same) ++nbad;
            std::vector<double> Ts;
            for (double T = 50.0; T <= 6500.0; T += 0.37) Ts.push_back(T);
            for (double b : {200.0, 1000.0, 6000.0}) for (double e : {-1e-9, 0.0, 1e-9}) Ts.push_back(b + e);
            for (double T : Ts) {
                const double a = cond_gas_h_mass(p.gas, T), r = thermo_h_mass(g, T);
                ++n; if (a != r) { ++nbad; worst = fmax(worst, fabs(a - r)/fabs(r)); }
            }
        }
        char b[160]; snprintf(b, sizeof b, "coefficients after datum and h(T) bit-identical (%ld points, Tref 0/298.15)", n);
        check(nbad == 0, b, worst, 0.0);
        // 外部 DB の種は起動時の検査で通るのは内蔵とビット一致のときだけ (speciesDB_attachCondensed; test_species_record_host で検査)
    }

    // ---------------------------------------------------------------- (b) 既知の気液差 (298.15 K)
    printf("== (b) L(298.15 K) = CEA H2O - H2O(L) absolute enthalpy difference ==\n");
    {
        // thermo.inp (McBride 2002; header 9/09/04) の係数を独立に書き写した参照: H2O 200–1000 K、H2O(L) 273.15–373.15 K
        const double gas[8] = {-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03, -7.342557370e-06,
                               4.955043490e-09, -1.336933246e-12, -3.303974310e+04};
        const double liq[8] = {1.326371304e+09, -2.448295388e+07, 1.879428776e+05, -7.678995050e+02, 1.761556813e+00,
                               -2.151167128e-03, 1.092570813e-06, 1.101760476e+08};
        const double T = 298.15, MW = 0.0180153;
        auto hmol = [](const double* a, double T) {
            return 8.314462618*T*(-a[0]/(T*T) + a[1]*log(T)/T + a[2] + a[3]*T/2 + a[4]*T*T/3 + a[5]*T*T*T/4 + a[6]*T*T*T*T/5 + a[7]/T);
        };
        const double Lref = (hmol(gas, T) - hmol(liq, T))/MW;
        double worst = 0.0;
        for (double Tref : {0.0, 298.15}) worst = fmax(worst, fabs(h2o_latent_pair(cond_test_latent_pair(Tref), T) - Lref)/Lref);
        check(worst <= 1e-12, "L(298.15) vs thermo.inp H2O - H2O(L) (independent literals; Tref 0 and 298.15)", worst, 1e-12);
        const double LHf = (-241826.0 + 285830.0)/MW;
        printf("      info: L(298.15) = %.6f J/kg; header Hf difference (-241826 + 285830 J/mol)/MW = %.6f J/kg (diff %.3f J/kg = NASA fit reproduction)\n",
               h2o_latent_pair(p0, T), LHf, h2o_latent_pair(p0, T) - LHf);
    }

    // ---------------------------------------------------------------- (新旧) 現行 h2o_latent との差
    printf("== new L vs pre-#10 h2o_latent (150-373.15 K) ==\n");
    printf("      %8s %16s %16s %14s\n", "T [K]", "L_new [J/kg]", "L_old [J/kg]", "dL [J/kg]");
    {
        double worstHi = 0.0;
        for (double T : {120.0, 150.0, 160.0, 170.0, 180.0, 190.0, 195.0, 199.0, 199.9, 200.0, 200.1, 210.0, 230.0, 250.0, 273.15, 298.15, 300.0, 350.0, 373.15}) {
            const double Ln = h2o_latent_pair(p0, T), Lo = legacy_h2o_latent_v0(T);
            printf("      %8.2f %16.6f %16.6f %14.6f\n", T, Ln, Lo, Ln - Lo);
            if (T >= 200.0) worstHi = fmax(worstHi, fabs(Ln - Lo)/Lo);
        }
        for (double T = 200.0; T <= 373.15 + 1e-9; T += 0.05) worstHi = fmax(worstHi, fabs(h2o_latent_pair(p0, T) - legacy_h2o_latent_v0(T))/legacy_h2o_latent_v0(T));
        check(worstHi <= 1e-12, "200-373.15 K: same coefficients -> rounding only (rel, every 0.05 K)", worstHi, 1e-12);
        const double d150 = h2o_latent_pair(p0, 150.0) - legacy_h2o_latent_v0(150.0);
        check(fabs(d150 - (-465.9201)) <= 1e-3, "dL(150 K) = -465.9201 J/kg (codex diagnose check value)", fabs(d150 - (-465.9201)), 1e-3);
    }

    // ---------------------------------------------------------------- (c') float 表 (200 K 近傍)
    printf("== (c') float table around the gas interval boundary 200 K ==\n");
    {
        CondSpeciesProps s = condProps_H2O();
        s.lat.h = &p0;
        CondTablesHost ht; cond_tables_build_host(s, ht);
        const CondTablesF tb = cond_tables_view_host(ht);
        const int i200 = (int)floor((200.0 - ht.T0)/ht.h);
        printf("      table cell containing 200 K: [%.2f, %.2f] K (T0 %.2f, h %.2f)\n", ht.T0 + ht.h*i200, ht.T0 + ht.h*(i200 + 1), ht.T0, ht.h);
        std::vector<float> Ts = {nextafterf(200.0f, 0.0f), 200.0f, nextafterf(200.0f, 1000.0f)};
        for (double T = 199.65; T <= 200.40 + 1e-9; T += 0.001) Ts.push_back((float)T);
        double wL = 0.0, wD = 0.0;
        for (float Tf : Ts) {
            if (!cond_tab_wet_ok(tb, Tf)) { printf("      unexpected: table not used at %.4f K\n", Tf); ++g_fail; continue; }
            float dLf = 0.0f;
            const double L = cond_tab_latent_f(tb, Tf, &dLf), Ld = cond_latent(s, (double)Tf);
            wL = fmax(wL, fabs(L - Ld)/Ld);
            const double d = 1e-3;
            const double dLd = (cond_latent(s, (double)Tf + d) - cond_latent(s, (double)Tf - d))/(2*d);
            wD = fmax(wD, fabs((double)dLf - dLd)/(2e-4*fabs(dLd) + 0.1));   // ≤ 1 で合格
        }
        check(wL <= 2e-6, "L rel (200 K, float neighbours, [199.65, 200.40] K every 1 mK)", wL, 2e-6);
        check(wD <= 1.0, "|dL'| / (2e-4|L'| + 0.1 J/kg/K) same points", wD, 1.0);
    }

    // ---------------------------------------------------------------- (d) datum 変換後の二相状態
    printf("== (d) wet state: datum conversion of the conserved energy keeps T, P and two-phase sound speed ==\n");
    {
        ResolvedSpeciesDB db = speciesDB_resolve(std::vector<std::string>{"N2", "H2O"}, "");
        double worstT = 0.0, worstP = 0.0, worstC = 0.0;
        for (double T : {150.0, 199.95, 230.0, 280.0}) for (double g : {1e-4, 0.01, 0.04}) {
            const double Y[2] = {0.95, 0.05};
            const double rho = 0.5;
            SpeciesThermo sp0[2] = {with_datum(db.species[0], 0.0), with_datum(db.species[1], 0.0)};
            SpeciesThermo sp1[2] = {with_datum(db.species[0], 298.15), with_datum(db.species[1], 298.15)};
            CondSpeciesProps c0 = condProps_H2O(), c1 = condProps_H2O();
            const CondLatentPair q0 = cond_test_latent_pair(0.0), q1 = cond_test_latent_pair(298.15);
            c0.lat.h = &q0; c1.lat.h = &q1;
            const double Rw = c0.R, Rmix = thermo_R_mix(sp0, 2, Y);
            const double e0 = thermo_e_mix(sp0, 2, Y, T) + g*(Rw*T - cond_latent(c0, T));
            // 保存エネルギーの datum 変換: e1 = e0 + Σ Y_s (h_s,1 − h_s,0) (液相分は総水の Y に含まれ、液相 h も同じ定数だけ動く)
            double de = 0.0;
            for (int s = 0; s < 2; ++s) de += Y[s]*(thermo_h_mass(sp1[s], T) - thermo_h_mass(sp0[s], T));
            const double e1 = e0 + de;
            const double T0 = cond_T_from_e_carrier(sp0, 2, Y, e0, g, Rw, c0, 0.9*T, 50.0, 6000.0);
            const double T1 = cond_T_from_e_carrier(sp1, 2, Y, e1, g, Rw, c1, 0.9*T, 50.0, 6000.0);
            worstT = fmax(worstT, fabs(T1 - T0));
            const double P0 = rho*T0*(Rmix - g*Rw), P1 = rho*T1*(Rmix - g*Rw);
            worstP = fmax(worstP, fabs(P1 - P0)/P0);
            double cp0, h0, cp1, h1, g20, c20, g21, c21;
            thermo_cph_mix(sp0, 2, Y, T0, &cp0, &h0); thermo_cph_mix(sp1, 2, Y, T1, &cp1, &h1);
            const double dL0 = (cond_latent(c0, T0 + 0.1) - cond_latent(c0, T0 - 0.1))/0.2;
            const double dL1 = (cond_latent(c1, T1 + 0.1) - cond_latent(c1, T1 - 0.1))/0.2;
            const bool ok0 = cond_twophase_sonic(cp0, Rmix - g*Rw, g, dL0, T0, &g20, &c20);
            const bool ok1 = cond_twophase_sonic(cp1, Rmix - g*Rw, g, dL1, T1, &g21, &c21);
            if (!ok0 || !ok1) { ++g_fail; printf("      sonic failed at T=%.2f g=%.1e\n", T, g); continue; }
            worstC = fmax(worstC, fabs(sqrt(c21) - sqrt(c20))/sqrt(c20));
        }
        check(worstT <= 1e-6, "inverted T after datum conversion vs before [K] (T 150-280 K, g 1e-4..0.04)", worstT, 1e-6);
        check(worstP <= 1e-9, "P = rho T (R_mix - g R_w) (rel)", worstP, 1e-9);
        check(worstC <= 1e-9, "two-phase frozen sound speed (rel)", worstC, 1e-9);
    }

    // ---------------------------------------------------------------- 未設定と device
    printf("== missing pair -> NaN; device == host ==\n");
    {
        const CondSpeciesProps s = condProps_H2O();
        const double Lnan = cond_latent(s, 250.0);
        check(std::isnan(Lnan), "H2O cond_latent without the gas-liquid pair is NaN", std::isnan(Lnan) ? 0.0 : 1.0, 0.0);
        const CondLatentRef r = cond_test_latent_ref(298.15, true);
        std::vector<double> Th;
        for (double T = 100.0; T <= 1300.0; T += 0.73) Th.push_back(T);
        const int n = (int)Th.size();
        double *dT = nullptr, *dL = nullptr;
        cudaMalloc((void**)&dT, n*sizeof(double)); cudaMalloc((void**)&dL, n*sizeof(double));
        cudaMemcpy(dT, Th.data(), n*sizeof(double), cudaMemcpyHostToDevice);
        latent_dev<<<(n + 127)/128, 128>>>(r, dT, dL, n);
        std::vector<double> Ld(n);
        const cudaError_t e = cudaMemcpy(Ld.data(), dL, n*sizeof(double), cudaMemcpyDeviceToHost);
        double worst = (e == cudaSuccess) ? 0.0 : 1.0;
        for (int i = 0; i < n && e == cudaSuccess; ++i) worst = fmax(worst, fabs(Ld[i] - h2o_latent(r, Th[i]))/h2o_latent(r, Th[i]));
        check(worst <= 1e-13, "device h2o_latent (ref.d) vs host (ref.h), 100-1300 K (rel)", worst, 1e-13);
        cudaFree(dT); cudaFree(dL);
    }

    printf("%s (%d failure%s)\n", g_fail ? "FAILED" : "ALL PASS", g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
