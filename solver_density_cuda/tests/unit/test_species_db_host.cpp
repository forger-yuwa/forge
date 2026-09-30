// =============================================================================
// test_species_db_host.cpp — 化学種 DB の host 側解決 (input/speciesDB.cpp) の単体試験 (GPU 不要)
//   (1) 内蔵 DB の MW (N2 / H2O) と名前解決 (完全一致・大文字小文字無視・未知種名の拒否・yaml 上書きの source)
//   (2) X→Y 換算: va3 組成 (X: H2O 6.09135e-2, N2 6.64860e-1, O2 2.16072e-1, AR 7.97588e-3, CO2 4.90034e-2, Σ≠1)
//       → Y_H2O = 0.03769539643 (rtol 1e-9; 内蔵 MW で再計算した値) と mole↔mass 往復
//   (3) bcondConfig floats の拒否条件 (X/Y 混在, X の種欠落, 負値, 非有限, 総和 0, 未知 index, Y の負値, |ΣY−1|>1e-3)
//   (4) 凝縮種の名前→index (ResolvedSpeciesDB::index; solverConfig::read の condensationSpecies と同じ大文字小文字無視の照合)
//   (5) 共通データの別名 (AR/HE/h2o/WATER/Air/air) が canonical と同値、CO と Co が別種、外部 DB は同じキーだけ上書き (plan #4)
//
// ビルド/実行 (単一 TU + speciesDB.cpp; 共通データの埋め込みヘッダを先に生成する):
//   cmake -DIN=solver_density_cuda/data/species/forge_species_v1.yaml -DOUT=/tmp/forge_species_gen/forge_species_data.hpp \
//       -P solver_density_cuda/cmake/embed_species_data.cmake
//   g++ -O1 -std=c++17 -I solver_density_cuda -I /tmp/forge_species_gen solver_density_cuda/tests/unit/test_species_db_host.cpp \
//       solver_density_cuda/input/speciesDB.cpp solver_density_cuda/input/speciesTransportDB.cpp -lyaml-cpp -o /tmp/test_species_db_host && /tmp/test_species_db_host
// 規約: [PASS]/[FAIL] を出し、失敗があれば非ゼロ終了。
// =============================================================================
#include "input/speciesDB.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

static int g_fail = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s %s\n", ok ? "[PASS]" : "[FAIL]", what.c_str());
    if (!ok) ++g_fail;
}

// floats ノードを YAML 文字列から作り、例外が出ることを判定する (メッセージも表示)
static void checkThrows(const std::string& floatsYaml, const ResolvedSpeciesDB& db, const std::string& what)
{
    try {
        YAML::Node n = YAML::Load(floatsYaml);
        (void)bcondSpeciesMassFractions(n, db, "inlet_test");
        check(false, what + " (no exception)");
    } catch (const std::exception& e) {
        check(true, what + ": " + e.what());
    }
}

int main()
{
    // ---- (1) 内蔵 DB と名前解決 ----
    {
        ResolvedSpeciesDB db = speciesDB_resolve({"N2", "H2O"}, "");
        check(db.size() == 2, "resolve builtin [N2,H2O] -> 2 species");
        check(std::fabs(db.MW(0) - 0.0280134) < 1e-15, "builtin N2 MW = 0.0280134");
        check(std::fabs(db.MW(1) - 0.0180153) < 1e-15, "builtin H2O MW = 0.0180153");
        check(db.source[0] == "builtin" && db.source[1] == "builtin", "source = builtin");
        check(db.index("H2O") == 1 && db.index("h2o") == 1 && db.index("N2") == 0 && db.index("XX") == -1,
              "index(): exact + case-insensitive lookup, -1 for unknown");

        ResolvedSpeciesDB db2 = speciesDB_resolve({"Ar", "co2", "he"}, "");
        check(db2.size() == 3 && std::fabs(db2.MW(0) - 0.039948) < 1e-15 && std::fabs(db2.MW(1) - 0.0440095) < 1e-15
              && std::fabs(db2.MW(2) - 0.0040026) < 1e-15, "aliases / case-insensitive builtin names (Ar, co2, he)");

        bool thrown = false;
        try { (void)speciesDB_resolve({"N2", "XENON"}, ""); } catch (const std::exception&) { thrown = true; }
        check(thrown, "unknown species name is rejected");

        thrown = false;
        try { (void)speciesDB_resolve({"N2", "n2"}, ""); } catch (const std::exception&) { thrown = true; }
        check(thrown, "duplicate species name (case-insensitive) is rejected");

        // yaml 上書き: MIXDRY (新規) と H2O (内蔵と同値) を持つ DB ファイル
        const std::string dbfile = "/tmp/test_species_db_host_db.yaml";
        {
            std::ofstream f(dbfile);
            f << "MIXDRY:\n  MW: 0.0298687837\n  LJ_sigma: 3.59\n  LJ_eps_kB: 111.7\n  Tlo: 200.0\n  Tmid: 1000.0\n  Thi: 6000.0\n"
                 "  nasa9_low: [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]\n  nasa9_high: [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]\n"
                 "H2O:\n  MW: 0.0180153\n  Tlo: 200.0\n  Tmid: 1000.0\n  Thi: 6000.0\n"
                 "  nasa9_low: [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]\n  nasa9_high: [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]\n";
        }
        ResolvedSpeciesDB db3 = speciesDB_resolve({"MIXDRY", "H2O", "N2"}, dbfile);
        check(db3.size() == 3 && db3.source[0] == "file" && db3.source[1] == "file" && db3.source[2] == "builtin",
              "speciesDBFile overlay: MIXDRY/H2O from file, N2 builtin");
        check(std::fabs(db3.MW(0) - 0.0298687837) < 1e-15 && db3.species[1].coef[0][0] == 1.0, "file entries override builtin coefficients");
    }

    // ---- (2) X→Y 換算 (va3) ----
    {
        ResolvedSpeciesDB db = speciesDB_resolve({"H2O", "N2", "O2", "AR", "CO2"}, "");
        const std::vector<double> X = {6.09135e-2, 6.64860e-1, 2.16072e-1, 7.97588e-3, 4.90034e-2};
        std::vector<double> MW(5);
        for (int s = 0; s < 5; ++s) MW[s] = db.MW(s);
        const std::vector<double> Y = speciesMoleToMass(X, MW);
        // 期待値は内蔵 MW (H2O 0.0180153, N2 0.0280134, O2 0.0319988, AR 0.039948, CO2 0.0440095) で
        // Y_k = X_k M_k / Σ X_j M_j を倍精度で再計算した値 (plan §6 (b) の 0.03769539643 と一致)。
        const double Yh2o_expected = 0.03769539643469918;
        check(std::fabs(Y[0] - Yh2o_expected) / Yh2o_expected < 1e-9, "va3: Y_H2O = 0.03769539643 (rtol 1e-9)");
        double sum = 0.0; for (double v : Y) sum += v;
        check(std::fabs(sum - 1.0) < 1e-14, "va3: sum(Y) = 1 (input sum X = 0.998825 is normalized)");
        const std::vector<double> Xb = speciesMassToMole(Y, MW);
        double sx = 0.0; for (double v : X) sx += v;
        double err = 0.0;
        for (int s = 0; s < 5; ++s) err = std::fmax(err, std::fabs(Xb[s] - X[s] / sx) / (X[s] / sx));
        check(err < 1e-12, "mole->mass->mole round trip (rtol 1e-12)");

        // bcond 経路 (floats ノード) で同じ換算
        YAML::Node n = YAML::Load("{X0: 6.09135e-2, X1: 6.64860e-1, X2: 2.16072e-1, X3: 7.97588e-3, X4: 4.90034e-2, Pt: 1.0}");
        const std::vector<double> Yb = bcondSpeciesMassFractions(n, db, "inlet");
        check(Yb.size() == 5 && std::fabs(Yb[0] - Yh2o_expected) / Yh2o_expected < 1e-9, "bcond floats X0..X4 -> Y (same value)");
    }

    // ---- (3) bcond 拒否条件 ----
    {
        ResolvedSpeciesDB db = speciesDB_resolve({"N2", "H2O"}, "");
        checkThrows("{X0: 0.9, X1: 0.1, Y0: 0.5}", db, "reject X/Y mixing: ");
        checkThrows("{X0: 0.9}", db, "reject X with missing species: ");
        checkThrows("{X0: 0.9, X1: -0.1}", db, "reject negative X: ");
        checkThrows("{X0: 0.9, X1: .nan}", db, "reject non-finite X: ");
        checkThrows("{X0: 0.0, X1: 0.0}", db, "reject sum(X)=0: ");
        checkThrows("{X0: 0.9, X1: 0.1, X5: 0.0}", db, "reject unknown index X5: ");
        checkThrows("{Y0: 0.9, Y1: -0.1}", db, "reject negative Y: ");
        checkThrows("{Y0: 0.9, Y1: 0.2}", db, "reject |sum(Y)-1| > 1e-3: ");
        checkThrows("{Y1: 0.1}", db, "reject Y1 only (default Y0=1 makes sum 1.1): ");
        checkThrows("{Y0: 0.9, Y7: 0.1}", db, "reject unknown index Y7: ");
        // 受理: Y 明示 (和 1)、Y0 のみ (= 1)、何も無し (空 → 呼び手が既定補完)
        std::vector<double> Y = bcondSpeciesMassFractions(YAML::Load("{Y0: 0.96230460, Y1: 0.03769540, Pt: 1.0}"), db, "in");
        check(Y.size() == 2 && std::fabs(Y[1] - 0.03769540) < 1e-12, "accept explicit Y0/Y1 (sum 1)");
        Y = bcondSpeciesMassFractions(YAML::Load("{Y0: 1.0}"), db, "in");
        check(Y.size() == 2 && Y[0] == 1.0 && Y[1] == 0.0, "accept Y0=1 only (Y1 defaults to 0)");
        Y = bcondSpeciesMassFractions(YAML::Load("{Pt: 1.0, Tt: 300.0}"), db, "in");
        check(Y.empty(), "neither X nor Y -> empty (caller applies default Y0=1)");
        Y = bcondSpeciesMassFractions(YAML::Load("{X0: 0.5, X1: 0.5}"), db, "in");
        check(Y.size() == 2 && std::fabs(Y[0] - 0.0280134 / (0.0280134 + 0.0180153)) < 1e-14, "X 50/50 -> Y by MW ratio");
    }

    // ---- (4) 凝縮種の名前→index ----
    {
        // 順序入替 (先頭 / 中間 / 末尾) で H2O の index が名前から正しく引ける
        for (const auto& names : std::vector<std::vector<std::string>>{{"H2O", "N2", "O2"}, {"N2", "H2O", "O2"}, {"N2", "O2", "H2O"}}) {
            ResolvedSpeciesDB db = speciesDB_resolve(names, "");
            int expect = 0; for (int s = 0; s < 3; ++s) if (names[s] == "H2O") expect = s;
            check(db.index("H2O") == expect && db.index("h2o") == expect,
                  "condensing species index by name (H2O at position " + std::to_string(expect) + ")");
        }
    }

    // ---- (5) 共通データ (data/species/forge_species_v1.yaml) の別名と大小文字の区別 (plan #4) ----
    {
        auto same = [](const SpeciesThermo& a, const SpeciesThermo& b) {
            return thermo_same_coeffs(a, b) && a.sigma_LJ == b.sigma_LJ && a.eps_kB == b.eps_kB;
        };
        const auto b = speciesDB_builtin();
        check(b.size() == 13, "builtin keys: 7 canonical IDs + 6 aliases (" + std::to_string(b.size()) + ")");
        const std::vector<std::pair<std::string, std::string>> aliases = {
            {"AR", "Ar"}, {"HE", "He"}, {"h2o", "H2O"}, {"WATER", "H2O"}, {"Air", "AIR"}, {"air", "AIR"}};
        for (const auto& a : aliases) {
            const ResolvedSpeciesDB r = speciesDB_resolve(std::vector<std::string>{a.first}, "");
            check(b.count(a.first) == 1 && same(b.at(a.first), b.at(a.second)) && same(r.species[0], b.at(a.second))
                  && r.names[0] == a.first, "alias " + a.first + " -> " + a.second + " (same values; name kept as written)");
        }
        check(!speciesDB_builtinDataName().empty() && speciesDB_builtinDataSha256().size() == 64,
              "embedded data provenance: " + speciesDB_builtinDataName() + " sha256 " + speciesDB_builtinDataSha256().substr(0, 16));

        // 大小文字だけ違う別種 (CEA の CO と Co) を外部 DB で与え、完全一致で別々に解決されること
        const std::string dbfile = "/tmp/test_species_db_host_co.yaml";
        {
            std::ofstream f(dbfile);
            f << "CO:\n  MW: 0.0280101\n  nasa9_low: [1,2,3,4,5,6,7,8,9]\n  nasa9_high: [1,2,3,4,5,6,7,8,9]\n"
                 "Co:\n  MW: 0.0589332\n  nasa9_low: [9,8,7,6,5,4,3,2,1]\n  nasa9_high: [9,8,7,6,5,4,3,2,1]\n"
                 "AR:\n  MW: 0.05\n  nasa9_low: [1,2,3,4,5,6,7,8,9]\n  nasa9_high: [1,2,3,4,5,6,7,8,9]\n";
        }
        const ResolvedSpeciesDB rco = speciesDB_resolve({"N2", "CO"}, dbfile);
        const ResolvedSpeciesDB rCo = speciesDB_resolve({"N2", "Co"}, dbfile);
        check(rco.MW(1) == 0.0280101 && rco.species[1].coef[0][0] == 1.0 && rCo.MW(1) == 0.0589332 && rCo.species[1].coef[0][0] == 9.0
              && rco.names[1] == "CO" && rCo.names[1] == "Co", "CO and Co resolve to different species (exact match)");
        // 外部 DB は同じキーだけを上書きする (従来どおり): AR は file、Ar は内蔵
        const ResolvedSpeciesDB rAR = speciesDB_resolve({"AR"}, dbfile);
        const ResolvedSpeciesDB rAr = speciesDB_resolve({"Ar"}, dbfile);
        check(rAR.source[0] == "file" && rAR.MW(0) == 0.05 && rAr.source[0] == "builtin" && rAr.MW(0) == 0.039948,
              "speciesDBFile key AR overrides only AR (Ar stays builtin)");
        // 従来の大小文字無視 (互換) は残る
        const ResolvedSpeciesDB rci = speciesDB_resolve({"co2", "n2"}, "");
        check(same(rci.species[0], b.at("CO2")) && same(rci.species[1], b.at("N2")), "legacy case-insensitive fallback (co2, n2)");
    }

    std::printf("%s (%d failures)\n", g_fail == 0 ? "ALL PASS" : "FAILED", g_fail);
    return g_fail == 0 ? 0 : 1;
}
