// =============================================================================
// dump_species_builtin.cpp — ソルバ内蔵 species DB の値を 16 進浮動小数 (%a) で JSON に書き出す (GPU 不要)
//   plans/active/thermophysics-solver-owned-species-db.md §5.1 #4 (共通データ化, 値は変えない) の前後比較用。
//   (a) speciesDB_builtin() の全キー (別名込み) の MW・温度区切り・LJ・全係数
//   (b) 名前解決 speciesDB_resolve({name}, "") の結果 (config が使う綴りと大小文字違い)
//   (c) 内蔵種だけの互換性ハッシュ (thermoHrefTemp 0 / 298.15)
//   比較と基準の保存は tests/unit/test_species_data_bitexact.py が行う (ここは出力だけ)。
//
// ビルド: test_species_data_bitexact.py がビルドする (共通データの埋め込みヘッダの生成を含む)。
// =============================================================================
#include "input/speciesDB.hpp"

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
    std::printf("{\"MW\": \"%a\", \"Tlo\": \"%a\", \"Tmid\": \"%a\", \"Thi\": \"%a\", \"LJ_sigma\": \"%a\", \"LJ_eps_kB\": \"%a\", \"low\": ",
                s.MW, s.Tlo, s.Tmid, s.Thi, s.sigma_LJ, s.eps_kB);
    hexArr(s.low, 9);
    std::printf(", \"high\": ");
    hexArr(s.high, 9);
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
    std::printf("\n}\n}\n");
    return 0;
}
