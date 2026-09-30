// =============================================================================
// dump_species_builtin.cpp — ソルバ内蔵 species DB の値を 16 進浮動小数 (%a) で JSON に書き出す (GPU 不要)
//   plans/active/thermophysics-solver-owned-species-db.md §5.1 #4 (共通データ化, 値は変えない) の前後比較用。
//   (a) speciesDB_builtin() の全キー (別名込み) の MW・温度区切り・LJ・全係数
//   (b) 名前解決 speciesDB_resolve({name}, "") の結果 (config が使う綴りと大小文字違い)
//   (c) 内蔵種だけの互換性ハッシュ (thermoHrefTemp 0 / 298.15)
//   (d) LJ: null の内蔵種 (#13-2) を kinetic 輸送・LJ の混合平均拡散に使ったときの起動時の可否
//   比較と基準の保存は tests/unit/test_species_data_bitexact.py が行う (ここは出力だけ)。
//
// ビルド: test_species_data_bitexact.py がビルドする (共通データの埋め込みヘッダの生成を含む)。
// =============================================================================
#include "input/speciesDB.hpp"
#include "input/speciesTransportDB.hpp"
#include "input/solverConfig.hpp"

#include <cstdio>
#include <string>
#include <vector>

static void hexArr(const double* a, int n)
{
    std::printf("[");
    for (int i = 0; i < n; ++i) std::printf("%s\"%a\"", i ? ", " : "", a[i]);
    std::printf("]");
}

static void entry(const SpeciesThermo& s)
{
    if (s.nInt == 2) {
        // 2 区間は #13-1 前と同じ JSON (基準 data/species_builtin_baseline.json と比較する)
        std::printf("{\"MW\": \"%a\", \"Tlo\": \"%a\", \"Tmid\": \"%a\", \"Thi\": \"%a\", \"LJ_sigma\": \"%a\", \"LJ_eps_kB\": \"%a\", \"low\": ",
                    s.MW, s.Tlo, s.Tbrk[0], s.Thi, s.sigma_LJ, s.eps_kB);
        hexArr(s.coef[0], 9);
        std::printf(", \"high\": ");
        hexArr(s.coef[1], 9);
        std::printf("}");
        return;
    }
    // 区間可変 (plan #13-1): 全境界と区間ごとの係数
    double Tb[THERMO_MAX_INTERVALS + 1];
    for (int k = 0; k <= s.nInt; ++k) Tb[k] = thermo_bound(s, k);
    std::printf("{\"MW\": \"%a\", \"Tbounds\": ", s.MW);
    hexArr(Tb, s.nInt + 1);
    std::printf(", \"LJ_sigma\": \"%a\", \"LJ_eps_kB\": \"%a\"", s.sigma_LJ, s.eps_kB);
    for (int k = 0; k < s.nInt; ++k) { std::printf(", \"coef%d\": ", k); hexArr(s.coef[k], 9); }
    std::printf("}");
}

int main()
{
    const auto db = speciesDB_builtin();
    std::printf("{\n\"builtin\": {\n");
    bool first = true;
    for (const auto& kv : db) {
        std::printf("%s  \"%s\": ", first ? "" : ",\n", kv.first.c_str());
        entry(kv.second);
        first = false;
    }
    std::printf("\n},\n\"resolve\": {\n");
    // config が使う綴り + 大小文字違い (従来の大小文字無視で解決される名前)
    const std::vector<std::string> names = {
        "N2", "O2", "AR", "Ar", "ar", "CO2", "co2", "HE", "He", "he", "H2O", "h2o", "WATER", "water", "Water",
        "AIR", "Air", "air", "n2", "o2", "H2", "CO", "Co", "XENON"};
    first = true;
    for (const auto& n : names) {
        std::printf("%s  \"%s\": ", first ? "" : ",\n", n.c_str());
        try {
            const ResolvedSpeciesDB r = speciesDB_resolve(std::vector<std::string>{n}, "");
            std::printf("{\"source\": \"%s\", \"entry\": ", r.source[0].c_str());
            entry(r.species[0]);
            std::printf("}");
        } catch (const std::exception&) {
            std::printf("\"not found\"");
        }
        first = false;
    }
    std::printf("\n},\n\"compat_hash\": {\n");
    const std::vector<std::vector<std::string>> sets = {
        {"N2", "H2O"}, {"N2", "O2", "AR", "CO2", "H2O"}, {"Ar", "He", "AIR", "WATER"}, {"air"}};
    first = true;
    for (const auto& set : sets) {
        std::string key;
        for (const auto& n : set) key += (key.empty() ? "" : ",") + n;
        const ResolvedSpeciesDB r = speciesDB_resolve(set, "");
        std::printf("%s  \"%s|0\": \"%s\",\n  \"%s|298.15\": \"%s\"", first ? "" : ",\n",
                    key.c_str(), speciesDB_compatHash(r, 0.0).c_str(), key.c_str(), speciesDB_compatHash(r, 298.15).c_str());
        first = false;
    }
    std::printf("\n},\n\"lj_use\": {\n");
    // (d) LJ: null の内蔵種 (#13-2) を LJ を読む使い方に回したときの起動時の可否 ("ok" / "rejected: <メッセージ>")。
    //     kinetic 輸送 (speciesTransportDB_resolve) と LJ の混合平均拡散 (speciesDB_resolve(cfg); 化学種 2 以上・viscMethod != 0・
    //     speciesDiffusionMethod 1)。test_species_data_bitexact.py が期待と照合する。
    auto verdict = [](auto&& fn) -> std::string {
        try { fn(); return "ok"; } catch (const std::exception& e) {
            std::string m = e.what(), o;
            for (char c : m) { if (c == '"' || c == '\\') o += '\\'; o += c; }
            return "rejected: " + o;
        }
    };
    auto kinetic = [&](const std::string& nm) {
        return verdict([&] {
            ResolvedSpeciesDB r = speciesDB_resolve(std::vector<std::string>{"N2", nm}, "");
            speciesTransportDB_resolve(r, {{"N2", "cea"}, {nm, "kinetic"}}, "");
        });
    };
    auto diffusion = [&](const std::vector<std::string>& names, const std::vector<SpeciesLumpSpec>& lumps, int visc, int diff) {
        return verdict([&] {
            solverConfig cfg;
            cfg.thermalMethod = 2; cfg.speciesNames = names; cfg.speciesLumps = lumps;
            cfg.nSpecies = static_cast<int>(names.size()); cfg.viscMethod = visc; cfg.speciesDiffusionMethod = diff;
            (void)speciesDB_resolve(cfg);
        });
    };
    SpeciesLumpSpec mixKr;
    mixKr.name = "MIXKR"; mixKr.basis = "mole"; mixKr.members = {"N2", "Kr"}; mixKr.fractions = {0.5, 0.5};
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"kinetic Kr (LJ null)", kinetic("Kr")},
        {"kinetic N (LJ from the Cantera table)", kinetic("N")},
        {"kinetic O2", kinetic("O2")},
        {"diffusion N2+Kr visc1 diff1", diffusion({"N2", "Kr"}, {}, 1, 1)},
        {"diffusion N2+Kr visc1 diff0", diffusion({"N2", "Kr"}, {}, 1, 0)},
        {"diffusion N2+Kr visc0 diff1", diffusion({"N2", "Kr"}, {}, 0, 1)},
        {"diffusion MIXKR(N2,Kr)+O2 visc1 diff1", diffusion({"MIXKR", "O2"}, {mixKr}, 1, 1)},
        {"diffusion N2+N visc1 diff1", diffusion({"N2", "N"}, {}, 1, 1)},
    };
    first = true;
    for (const auto& c : cases) {
        std::printf("%s  \"%s\": \"%s\"", first ? "" : ",\n", c.first.c_str(), c.second.c_str());
        first = false;
    }
    std::printf("\n}\n}\n");
    return 0;
}
