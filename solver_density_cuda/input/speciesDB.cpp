#include "input/speciesDB.hpp"
#include "input/solverConfig.hpp"

// 共通 species データの埋め込み (ビルド時生成; cmake/embed_species_data.cmake)。
// CMake では cuda_forge の生成ディレクトリ、単独ビルドでは -I <生成先> で見つける。
#include "forge_species_data.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {

// 内蔵 DB の 1 種を構築する。
// NASA-9 係数 (a0..a8): cp/R = a0/T^2 + a1/T + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
// MW [kg/mol], sigma_LJ [Angstrom], eps_kB [K]。Tb[0..nInt] は区間の境界、a[k] は区間 k の係数 (区間可変, plan #13-1)。
// 値は共通データ data/species/forge_species_v1.yaml から (各 species・各定数の出典は同ファイルの source/LJ.source と
// methods/thermophysics.md 「内蔵 species DB の一覧と出典」)。
SpeciesThermo makeSpecies(double MW, double sigma, double eps_kB, int nInt, const double* Tb, const double (*a)[9])
{
    SpeciesThermo s{};
    s.MW = MW; s.sigma_LJ = sigma; s.eps_kB = eps_kB;
    thermo_set_intervals(s, nInt, Tb, a);
    s.h_datum = 0.0;
    s.invMW = 1.0/MW;
    return s;
}

// 区間の境界 (Tlo, 区切り..., Thi) を double の列で返す。
std::vector<double> boundsOf(const SpeciesThermo& s)
{
    std::vector<double> b;
    for (int k = 0; k <= s.nInt; ++k) b.push_back(thermo_bound(s, k));
    return b;
}

// 区間数・境界の検査 (空 = 正常)。1 <= nInt <= THERMO_MAX_INTERVALS、境界は有限・正・狭義単調増加。
std::string intervalProblem(int nInt, const std::vector<double>& Tb)
{
    if (nInt < 1 || nInt > THERMO_MAX_INTERVALS) {
        return "has " + std::to_string(nInt) + " temperature intervals; the solver supports 1 to "
               + std::to_string(THERMO_MAX_INTERVALS) + " (THERMO_MAX_INTERVALS)";
    }
    if (static_cast<int>(Tb.size()) != nInt + 1) return "needs " + std::to_string(nInt + 1) + " interval bounds";
    for (int k = 0; k <= nInt; ++k) {
        if (!std::isfinite(Tb[k]) || !(Tb[k] > 0.0)) return "has a non-positive or non-finite interval bound";
        if (k > 0 && !(Tb[k-1] < Tb[k])) return "intervals must be contiguous and increasing";
    }
    return "";
}

// yaml ノードから 9 係数を読む
bool read9(const YAML::Node& n, double out[9])
{
    if (!n || !n.IsSequence() || n.size() != 9) return false;
    for (int i=0;i<9;i++) out[i] = n[i].as<double>();
    return true;
}

std::string toUpper(std::string s)
{
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// 大文字小文字無視で map を引く (完全一致を優先)。
std::map<std::string, SpeciesThermo>::const_iterator findSpecies(const std::map<std::string, SpeciesThermo>& db,
                                                                  const std::string& name)
{
    auto it = db.find(name);
    if (it != db.end()) return it;
    const std::string up = toUpper(name);
    for (auto jt = db.begin(); jt != db.end(); ++jt) {
        if (toUpper(jt->first) == up) return jt;
    }
    return db.end();
}

// bcond floats のキー "X<digits>" / "Y<digits>" を index に分解する (該当しなければ -1)。
int fractionIndex(const std::string& key, char prefix)
{
    if (key.size() < 2 || key[0] != prefix) return -1;
    for (size_t i = 1; i < key.size(); ++i) if (!std::isdigit(static_cast<unsigned char>(key[i]))) return -1;
    return std::atoi(key.c_str() + 1);
}

std::unique_ptr<ResolvedSpeciesDB> g_current;   // speciesDB_init の保持先

// ---- 共通 species データ (埋め込み YAML) の解析 ----
//   スキーマは data/species/forge_species_v1.yaml の冒頭。#4 (値は変えない) では legacy_builtin に solver を含む
//   エントリだけを内蔵種とする (移行前のハードコード 7 種と同じ集合)。
#define FORGE_SPECIES_DATA_SCHEMA "forge_species_data_v1"

struct BuiltinEntry {
    std::string              id;        // canonical ID (大小文字を区別)
    std::vector<std::string> aliases;   // 完全一致の別名
    SpeciesThermo            sp;
};

std::vector<BuiltinEntry> parseBuiltinData()
{
    const std::string where = std::string("[speciesDB] built-in species data ") + kForgeSpeciesDataName;
    YAML::Node root;
    try {
        root = YAML::Load(kForgeSpeciesDataYaml);
    } catch (const std::exception& e) {
        throw std::runtime_error(where + " is not valid YAML: " + e.what());
    }
    if (!root["schema"] || root["schema"].as<std::string>() != FORGE_SPECIES_DATA_SCHEMA) {
        throw std::runtime_error(where + ": schema must be '" FORGE_SPECIES_DATA_SCHEMA "'");
    }
    if (!root["species"] || !root["species"].IsSequence()) throw std::runtime_error(where + ": 'species' list is missing");
    std::vector<BuiltinEntry> out;
    std::map<std::string, std::string> seen;   // 名前 (ID と別名) → 所有 ID。完全一致で重複を拒否
    for (const auto& n : root["species"]) {
        const std::string id = n["id"].as<std::string>();
        auto claim = [&](const std::string& nm) {
            auto it = seen.find(nm);
            if (it != seen.end()) throw std::runtime_error(where + ": name '" + nm + "' is used by both '" + it->second + "' and '" + id + "'");
            seen[nm] = id;
        };
        claim(id);
        std::vector<std::string> aliases;
        if (n["aliases"]) for (const auto& a : n["aliases"]) { aliases.push_back(a.as<std::string>()); claim(aliases.back()); }
        bool solver = false;
        if (n["legacy_builtin"]) for (const auto& t : n["legacy_builtin"]) solver = solver || (t.as<std::string>() == "solver");
        if (!solver) continue;
        const std::string phase = n["phase"] ? n["phase"].as<std::string>() : "";
        if (phase != "gas") throw std::runtime_error(where + ": '" + id + "' phase '" + phase + "' (only gas is supported as a built-in species)");
        const YAML::Node iv = n["intervals"];
        // 区間可変 (1..THERMO_MAX_INTERVALS; plan #13-1)。上限を超える種は起動時に拒否する。
        if (!iv || !iv.IsSequence() || iv.size() < 1) {
            throw std::runtime_error(where + ": '" + id + "' needs a list of temperature intervals");
        }
        const int nInt = static_cast<int>(iv.size());
        if (nInt > THERMO_MAX_INTERVALS) throw std::runtime_error(where + ": '" + id + "' " + intervalProblem(nInt, {}));
        double a[THERMO_MAX_INTERVALS][9];
        std::vector<double> Tb;
        for (int k = 0; k < nInt; ++k) {
            if (!read9(iv[k]["coeffs"], a[k])) {
                throw std::runtime_error(where + ": '" + id + "' needs 9 NASA-9 coefficients per interval");
            }
            const double lo = iv[k]["Tlo"].as<double>(), hi = iv[k]["Thi"].as<double>();
            if (k == 0) Tb.push_back(lo);
            else if (lo != Tb.back()) throw std::runtime_error(where + ": '" + id + "' intervals must be contiguous and increasing");
            Tb.push_back(hi);
        }
        {
            const std::string bad = intervalProblem(nInt, Tb);
            if (!bad.empty()) throw std::runtime_error(where + ": '" + id + "' " + bad);
        }
        const double MW = n["MW"].as<double>();
        if (!(MW > 0.0) || !std::isfinite(MW)) throw std::runtime_error(where + ": '" + id + "' has invalid MW");
        // LJ: null は「輸送データなし」。内蔵種での扱い (輸送に使うと拒否) は plan #6/#7 なので、現時点では内蔵種に置かない。
        const YAML::Node lj = n["LJ"];
        if (!lj || lj.IsNull()) {
            throw std::runtime_error(where + ": '" + id + "' has no LJ data (LJ: null); built-in species without transport data are not supported yet (plan #6/#7)");
        }
        out.push_back({id, aliases, makeSpecies(MW, lj["sigma"].as<double>(), lj["eps_kB"].as<double>(), nInt, Tb.data(), a)});
    }
    return out;
}

// 起動時に 1 回だけ解析する (失敗は std::runtime_error; speciesDB_init はメッセージを出して exit)。
const std::vector<BuiltinEntry>& builtinEntries()
{
    static const std::vector<BuiltinEntry> es = parseBuiltinData();
    return es;
}

// ---- 共通データの凝縮相エントリ (phase: condensed; plan #10 §4.8) ----
//   内蔵の気相種 (legacy_builtin: solver) とは別に読む。気液ペアの気相 (pair_of) は内蔵の気相種でなければならない。
#define SPECIES_CONDENSED_BELOW "linear_cp_fd_at_Tlo"   // 実装が受け付ける延長規約 (condensationProperties_d.cuh cond_latent_pair_make)
#define SPECIES_CONDENSED_ABOVE "hold_at_Thi"

struct CondensedEntry {
    std::string id;       // canonical ID ("H2O(L)")
    std::string pairOf;   // 気相の canonical ID ("H2O")
    double      MW = 0.0, Tlo = 0.0, Thi = 0.0;
    double      coeffs[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::string below, above;
};

std::vector<CondensedEntry> parseCondensedData()
{
    const std::string where = std::string("[speciesDB] built-in species data ") + kForgeSpeciesDataName;
    const YAML::Node root = YAML::Load(kForgeSpeciesDataYaml);   // 構文・スキーマは parseBuiltinData が先に検査済み
    std::vector<CondensedEntry> out;
    for (const auto& n : root["species"]) {
        if (!n["phase"] || n["phase"].as<std::string>() != "condensed") continue;
        CondensedEntry c;
        c.id = n["id"].as<std::string>();
        const std::string w = where + ": condensed '" + c.id + "'";
        if (!n["pair_of"]) throw std::runtime_error(w + " has no pair_of (the gas species it pairs with)");
        c.pairOf = n["pair_of"].as<std::string>();
        c.MW = n["MW"].as<double>();
        const YAML::Node iv = n["intervals"];
        if (!iv || !iv.IsSequence() || iv.size() != 1 || !read9(iv[0]["coeffs"], c.coeffs)) {
            throw std::runtime_error(w + " must have exactly 1 temperature interval with 9 NASA-9 coefficients");
        }
        c.Tlo = iv[0]["Tlo"].as<double>(); c.Thi = iv[0]["Thi"].as<double>();
        if (!(c.Tlo > 0.0 && c.Tlo < c.Thi)) throw std::runtime_error(w + ": invalid interval");
        const YAML::Node ex = n["extension"];
        c.below = (ex && ex["below"]) ? ex["below"].as<std::string>() : "";
        c.above = (ex && ex["above"]) ? ex["above"].as<std::string>() : "";
        if (c.below != SPECIES_CONDENSED_BELOW || c.above != SPECIES_CONDENSED_ABOVE) {
            throw std::runtime_error(w + ": extension {below: '" + c.below + "', above: '" + c.above + "'} is not implemented (only below: "
                                     SPECIES_CONDENSED_BELOW ", above: " SPECIES_CONDENSED_ABOVE ")");
        }
        out.push_back(c);
    }
    return out;
}

const std::vector<CondensedEntry>& condensedEntries()
{
    builtinEntries();   // 共通データの構文・スキーマ検査を先に通す
    static const std::vector<CondensedEntry> es = parseCondensedData();
    return es;
}

} // anonymous namespace

int ResolvedSpeciesDB::index(const std::string& name) const
{
    for (int s = 0; s < size(); ++s) if (names[s] == name) return s;
    const std::string up = toUpper(name);
    for (int s = 0; s < size(); ++s) if (toUpper(names[s]) == up) return s;
    return -1;
}

std::map<std::string, SpeciesThermo> speciesDB_builtin()
{
    // 共通データ (data/species/forge_species_v1.yaml, ビルド時に埋め込み) の solver 集合を、canonical ID と
    // 別名の両方をキーにして返す (移行前のハードコード表と同じキー集合・同じ値; tests/unit/test_species_data_bitexact.py)。
    std::map<std::string, SpeciesThermo> db;
    for (const auto& e : builtinEntries()) {
        db[e.id] = e.sp;
        for (const auto& a : e.aliases) db[a] = e.sp;
    }
    return db;
}

const std::string& speciesDB_builtinDataName()
{
    static const std::string s = kForgeSpeciesDataName;
    return s;
}

const std::string& speciesDB_builtinDataSha256()
{
    static const std::string s = kForgeSpeciesDataSha256;
    return s;
}

ResolvedSpeciesDB speciesDB_resolve(const std::vector<std::string>& namesIn, const std::string& dbFile)
{
    return speciesDB_resolve(namesIn, dbFile, std::vector<SpeciesLumpSpec>{});
}

namespace {

std::string g17(double v);   // 下の「解決済み記録」節で定義 (%.17g)

// 構成種の重複検査用のキー: 内蔵種は canonical ID (別名・大小文字違いを同一視)、それ以外は大文字化した名前。
std::string memberKey(const std::string& resolvedKey, bool fromFile)
{
    if (!fromFile) {
        for (const auto& e : builtinEntries()) {
            if (e.id == resolvedKey) return e.id;
            for (const auto& a : e.aliases) if (a == resolvedKey) return e.id;
        }
    }
    return toUpper(resolvedKey);
}

std::string breakpoints(const SpeciesThermo& s)
{
    std::ostringstream o;
    o << std::setprecision(17);
    for (int k = 0; k <= s.nInt; ++k) o << (k ? "/" : "") << thermo_bound(s, k);
    return o.str();
}

// 区間の境界と数が完全に同じか
bool sameIntervals(const SpeciesThermo& a, const SpeciesThermo& b)
{
    if (a.nInt != b.nInt) return false;
    for (int k = 0; k <= a.nInt; ++k) if (thermo_bound(a, k) != thermo_bound(b, k)) return false;
    return true;
}

// 構成種 m の [Tlo,Thi] の外 (below: T<Tlo、それ以外 T>Thi) の延長 (cp 一定・h 線形・s° 対数; thermo_d.cuh の外挿) を
// NASA-9 の 1 区間として表す: cp/R = a2、h/RT = a2 + a7/T、s/R = a2 ln T + a8。端の値は thermo_*_clamped と同じ評価。
void extensionCoeffs(const SpeciesThermo& m, bool below, double a[9])
{
    const double Tb = below ? m.Tlo : m.Thi;
    const double cpR = thermo_cp_molar_clamped(m, Tb)/THERMO_RU;
    const double hR  = thermo_h_molar_clamped(m, Tb)/THERMO_RU;
    const double sR  = thermo_s_molar_clamped(m, Tb)/THERMO_RU;
    for (int i = 0; i < 9; ++i) a[i] = 0.0;
    a[2] = cpR;
    a[7] = hR - cpR*Tb;
    a[8] = sR - cpR*std::log(Tb);
}

// lump を合成する (plan thermophysics-solver-owned-species-db §4.2 #6a; 式は設計側 composition.lump_entry と同じ)。
//   x_k: lump 内モル分率 (mole は正規化、mass は Y_k/M_k を正規化)
//   MW  = Σ x_k M_k、NASA-9 係数 = Σ x_k a_k (区間ごと。係数に線形なので cp/R, h/RT, s°/R は構成種のモル加重和と厳密に一致)
//   LJ  = Σ Y_k σ_k, Σ Y_k ε_k (Y_k = x_k M_k / MW)。**暫定**: 設計側の現行規約に合わせただけで根拠は無い。
//         plan #7 で粘性・熱伝導は実種展開、拡散は Blanc の法則に置き換える。
//   区切り: 構成種の区間がすべて同じなら同じ区間で畳む (SPECIES_LUMP_SYNTHESIS; #6a と同じ演算)。違えば区切りの和集合で畳み
//   (SPECIES_LUMP_SYNTHESIS_UNION; plan #6b → #13-1)、構成種が自分の [Tlo,Thi] の外にある和集合区間では外挿 (extensionCoeffs) を
//   その区間の係数として足す (lump の外挿も構成種の外挿の和と一致する)。和集合の区間数が THERMO_MAX_INTERVALS を超えたら拒否。
//   synthesis に使った規約文字列を返す (記録・互換性ハッシュに入る)。
SpeciesThermo synthesizeLump(const std::string& name, const std::vector<double>& x, const std::vector<SpeciesThermo>& m,
                             std::string& synthesis)
{
    SpeciesThermo s{};
    double MW = 0.0;
    for (size_t k = 0; k < m.size(); ++k) MW += x[k]*m[k].MW;
    if (!(MW > 0.0) || !std::isfinite(MW)) throw std::runtime_error("[speciesDB] lump '" + name + "': synthesized MW is not positive");
    s.MW = MW;
    bool same = true;
    for (size_t k = 1; k < m.size(); ++k) same = same && sameIntervals(m[0], m[k]);
    double a[THERMO_MAX_INTERVALS][9] = {};
    std::vector<double> Tb;
    if (same) {
        synthesis = SPECIES_LUMP_SYNTHESIS;
        Tb = boundsOf(m[0]);
        for (size_t k = 0; k < m.size(); ++k)
            for (int j = 0; j < m[0].nInt; ++j)
                for (int i = 0; i < 9; ++i) a[j][i] += x[k]*m[k].coef[j][i];
    } else {
        synthesis = SPECIES_LUMP_SYNTHESIS_UNION;
        for (const auto& mk : m) for (double b : boundsOf(mk)) Tb.push_back(b);
        std::sort(Tb.begin(), Tb.end());
        Tb.erase(std::unique(Tb.begin(), Tb.end()), Tb.end());
        const int nU = static_cast<int>(Tb.size()) - 1;
        if (nU > THERMO_MAX_INTERVALS) {
            std::string lst;
            for (size_t k = 0; k < Tb.size(); ++k) lst += (k ? "/" : "") + g17(Tb[k]);
            throw std::runtime_error("[speciesDB] lump '" + name + "': the union of the constituents' temperature breakpoints (" + lst
                                     + " K) makes " + std::to_string(nU) + " intervals, more than THERMO_MAX_INTERVALS="
                                     + std::to_string(THERMO_MAX_INTERVALS) + "; keep such species as separate species.");
        }
        for (int j = 0; j < nU; ++j) {
            for (size_t k = 0; k < m.size(); ++k) {
                double ak[9];
                if (Tb[j+1] <= m[k].Tlo)     extensionCoeffs(m[k], true, ak);
                else if (Tb[j] >= m[k].Thi)  extensionCoeffs(m[k], false, ak);
                else { const double* c = m[k].coef[thermo_interval(m[k], Tb[j])]; for (int i = 0; i < 9; ++i) ak[i] = c[i]; }
                for (int i = 0; i < 9; ++i) a[j][i] += x[k]*ak[i];
            }
        }
    }
    thermo_set_intervals(s, static_cast<int>(Tb.size()) - 1, Tb.data(), a);
    double sig = 0.0, eps = 0.0;
    for (size_t k = 0; k < m.size(); ++k) {
        const double Yk = x[k]*m[k].MW/MW;
        sig += Yk*m[k].sigma_LJ;
        eps += Yk*m[k].eps_kB;
    }
    s.sigma_LJ = sig; s.eps_kB = eps;
    s.h_datum = 0.0;
    s.invMW = 1.0/MW;
    return s;
}

} // anonymous namespace

std::string speciesDB_identityKey(const std::string& dbKey, bool fromFile)
{
    return memberKey(dbKey, fromFile);
}

ResolvedSpeciesDB speciesDB_resolve(const std::vector<std::string>& namesIn, const std::string& dbFile,
                                    const std::vector<SpeciesLumpSpec>& lumps)
{
    std::map<std::string, SpeciesThermo> db = speciesDB_builtin();
    std::map<std::string, bool> fromFile;   // 上書き/追加された種名

    // 外部 DB ファイルがあれば内蔵 DB を上書き/追加
    if (!dbFile.empty()) {
        YAML::Node root;
        try {
            root = YAML::LoadFile(dbFile);
        } catch (const std::exception& e) {
            throw std::runtime_error("[speciesDB] failed to read speciesDBFile '" + dbFile + "': " + e.what());
        }
        for (auto it = root.begin(); it != root.end(); ++it) {
            std::string name = it->first.as<std::string>();
            YAML::Node s = it->second;
            SpeciesThermo sp{};
            sp.MW       = s["MW"]       ? s["MW"].as<double>()       : 0.0;
            sp.sigma_LJ = s["LJ_sigma"] ? s["LJ_sigma"].as<double>() : 3.6;
            sp.eps_kB   = s["LJ_eps_kB"]? s["LJ_eps_kB"].as<double>(): 97.0;
            const std::string w = "[speciesDB] species '" + name + "' in DB file '" + dbFile + "'";
            bool ok = false;
            if (s["Tbounds"] || s["nasa9_intervals"]) {
                // 区間可変の書式 (plan #13-1; 解決済み記録の nInt != 2 の種と同じキー): Tbounds [Tlo, 区切り..., Thi] と
                // nasa9_intervals [[a0..a8] x nInt]。2 区間の書式 (Tlo/Tmid/Thi, nasa9_low/nasa9_high) と混ぜられない。
                if (s["nasa9_low"] || s["nasa9_high"] || s["Tmid"] || s["Tlo"] || s["Thi"]) {
                    throw std::runtime_error(w + ": Tbounds/nasa9_intervals cannot be mixed with Tlo/Tmid/Thi/nasa9_low/nasa9_high");
                }
                const YAML::Node tb = s["Tbounds"], iv = s["nasa9_intervals"];
                if (!tb || !tb.IsSequence() || !iv || !iv.IsSequence()) throw std::runtime_error(w + ": needs both Tbounds and nasa9_intervals lists");
                const int nInt = static_cast<int>(iv.size());
                std::vector<double> Tb;
                for (const auto& v : tb) Tb.push_back(v.as<double>());
                const std::string bad = intervalProblem(nInt, Tb);
                if (!bad.empty()) throw std::runtime_error(w + " " + bad);
                double a[THERMO_MAX_INTERVALS][9];
                for (int k = 0; k < nInt; ++k)
                    if (!read9(iv[k], a[k])) throw std::runtime_error(w + ": nasa9_intervals[" + std::to_string(k) + "] needs 9 coefficients");
                thermo_set_intervals(sp, nInt, Tb.data(), a);
                ok = true;
            } else {
                double lo[9], hi[9];
                if (read9(s["nasa9_low"], lo) && read9(s["nasa9_high"], hi)) {
                    thermo_set_nasa9_2(sp, s["Tlo"]  ? s["Tlo"].as<double>()  : 200.0,
                                           s["Tmid"] ? s["Tmid"].as<double>() : 1000.0,
                                           s["Thi"]  ? s["Thi"].as<double>()  : 6000.0, lo, hi);
                    ok = true;
                }
            }
            if (ok) {
                if (!(sp.MW > 0.0) || !std::isfinite(sp.MW)) {
                    throw std::runtime_error("[speciesDB] species '" + name + "' has invalid MW in DB file '" + dbFile + "'");
                }
                sp.h_datum = 0.0;
                sp.invMW = 1.0/sp.MW;
                db[name] = sp;
                fromFile[name] = true;
            } else {
                std::cerr << "[speciesDB] species '" << name
                          << "' in DB file lacks valid nasa9_low/nasa9_high (need 9 coeffs) or Tbounds/nasa9_intervals; ignored" << std::endl;
            }
        }
    }

    std::vector<std::string> names = namesIn;
    if (names.empty()) names.push_back("N2");  // ダミー (calorically-perfect 用)

    if (static_cast<int>(names.size()) > THERMO_MAX_SPECIES) {
        throw std::runtime_error("[speciesDB] nSpecies=" + std::to_string(names.size())
                                 + " exceeds THERMO_MAX_SPECIES=" + std::to_string(THERMO_MAX_SPECIES));
    }

    auto notFound = [&](const std::string& nm, const std::string& what) {
        std::string avail;
        for (const auto& kv : db) avail += (avail.empty() ? "" : " ") + kv.first;
        return std::runtime_error("[speciesDB] " + what + " '" + nm + "' not found in built-in DB nor speciesDBFile"
                                  + (dbFile.empty() ? std::string("") : " '" + dbFile + "'")
                                  + " (available: " + avail + ")");
    };

    // lump 指定の対応付け (lump 名は names にちょうど 1 回現れる)
    std::map<std::string, const SpeciesLumpSpec*> lumpOf;
    for (const auto& lp : lumps) {
        if (lumpOf.count(lp.name)) throw std::runtime_error("[speciesDB] lump '" + lp.name + "' is defined twice in physProp.species");
        if (std::find(names.begin(), names.end(), lp.name) == names.end()) {
            throw std::runtime_error("[speciesDB] lump '" + lp.name + "' is not in the species list");
        }
        lumpOf[lp.name] = &lp;
    }

    ResolvedSpeciesDB out;
    for (const auto& nm : names) {
        if (out.index(nm) >= 0) {
            throw std::runtime_error("[speciesDB] species '" + nm + "' is listed twice in physProp.species");
        }
        auto lit = lumpOf.find(nm);
        if (lit == lumpOf.end()) {
            auto it = findSpecies(db, nm);
            if (it == db.end()) throw notFound(nm, "species");
            out.names.push_back(nm);
            out.species.push_back(it->second);
            out.source.push_back(fromFile.count(it->first) ? "file" : "builtin");
            out.lumps.push_back(ResolvedLump{});
            out.dbKey.push_back(it->first);
            continue;
        }
        // ---- lump: 検査 → lump 内モル分率 → 合成 ----
        const SpeciesLumpSpec& lp = *lit->second;
        const std::string where = "[speciesDB] lump '" + nm + "'";
        {
            auto ct = findSpecies(db, nm);
            if (ct != db.end()) {
                throw std::runtime_error(where + ": the name collides with " + std::string(fromFile.count(ct->first) ? "speciesDBFile" : "built-in")
                                         + " species '" + ct->first + "' (names are matched case-insensitively); choose another lump name");
            }
        }
        if (lp.basis != "mole" && lp.basis != "mass") {
            throw std::runtime_error(where + ": basis must be 'mole' or 'mass' (got '" + lp.basis + "')");
        }
        if (lp.members.empty() || lp.members.size() != lp.fractions.size()) {
            throw std::runtime_error(where + ": needs at least one constituent with a fraction");
        }
        ResolvedLump rl;
        rl.basis = lp.basis;
        rl.members = lp.members;
        rl.input = lp.fractions;
        std::map<std::string, std::string> seenKey;   // 重複検査 (内蔵の別名・大小文字違いを同一視)
        double sum = 0.0;
        for (size_t k = 0; k < lp.members.size(); ++k) {
            const std::string& mn = lp.members[k];
            const double v = lp.fractions[k];
            if (!std::isfinite(v) || !(v > 0.0)) {
                std::ostringstream o;
                o << where << ": fraction of '" << mn << "' must be positive and finite (got " << std::setprecision(17) << v << ")";
                throw std::runtime_error(o.str());
            }
            auto it = findSpecies(db, mn);
            if (it == db.end()) throw notFound(mn, "lump '" + nm + "' constituent");
            const bool file = fromFile.count(it->first) > 0;
            const std::string key = memberKey(it->first, file);
            auto sk = seenKey.find(key);
            if (sk != seenKey.end()) {
                throw std::runtime_error(where + ": constituent '" + mn + "' is the same species as '" + sk->second + "' (listed twice)");
            }
            seenKey[key] = mn;
            rl.memberSpecies.push_back(it->second);
            rl.memberSource.push_back(file ? "file" : "builtin");
            rl.memberDbKey.push_back(it->first);
            sum += v;
        }
        if (std::fabs(sum - 1.0) >= 1.0e-3) {
            std::cerr << where << " WARNING: " << lp.basis << " fractions sum to " << std::setprecision(10) << sum
                      << " (differs from 1 by >= 1e-3); normalized within the lump" << std::endl;
        }
        // lump 内モル分率
        const size_t nm_ = lp.members.size();
        rl.x.assign(nm_, 0.0);
        if (lp.basis == "mole") {
            for (size_t k = 0; k < nm_; ++k) rl.x[k] = lp.fractions[k]/sum;
        } else {
            double den = 0.0;
            for (size_t k = 0; k < nm_; ++k) { rl.x[k] = (lp.fractions[k]/sum)/rl.memberSpecies[k].MW; den += rl.x[k]; }
            for (size_t k = 0; k < nm_; ++k) rl.x[k] /= den;
        }
        out.names.push_back(nm);
        out.species.push_back(synthesizeLump(nm, rl.x, rl.memberSpecies, rl.synthesis));
        out.source.push_back("lump");
        out.lumps.push_back(rl);
        out.dbKey.push_back("");
    }
    for (auto& s : out.species) s.invMW = 1.0/s.MW;
    return out;
}

ResolvedSpeciesDB speciesDB_resolve(const solverConfig& cfg)
{
    ResolvedSpeciesDB db = speciesDB_resolve(cfg.speciesNames, cfg.speciesDBFile, cfg.speciesLumps);
    // 凝縮種は lump に入れられない (独立種として置く; methods/thermophysics.md §1b.2)。
    // 凝縮種名は solverConfig::read が condensationSpecies / condGasSpecies / 単成分 TP から解決済み (condGasSpeciesName)。
    if (cfg.condensation == 1) {
        std::vector<std::string> cond;
        if (!cfg.condGasSpeciesName.empty()) cond.push_back(cfg.condGasSpeciesName);
        if (!cfg.condensationSpecies.empty()) cond.push_back(cfg.condensationSpecies);
        auto keyOf = [&](const std::string& n) {
            // 内蔵種は canonical ID (H2O/WATER, Ar/AR を同一視)、それ以外は大文字化
            const std::string up = toUpper(n);
            for (const auto& e : builtinEntries()) {
                if (toUpper(e.id) == up) return e.id;
                for (const auto& a : e.aliases) if (toUpper(a) == up) return e.id;
            }
            return up;
        };
        for (int s = 0; s < db.size(); ++s) {
            if (!db.isLump(s)) continue;
            for (const auto& mn : db.lumps[s].members) {
                for (const auto& c : cond) {
                    if (keyOf(mn) == keyOf(c)) {
                        throw std::runtime_error("[speciesDB] lump '" + db.names[s] + "': constituent '" + mn
                                                 + "' is the condensing species (condensation: 1, '" + c
                                                 + "'); a condensing species cannot be lumped — list it as a separate species");
                    }
                }
            }
        }
    }
    // 凝縮種 H2O の液相 (気液ペア; plan #10 §4.8)。TP はペアの気相種が種リストに要る (潜熱の h_v = 種 DB の気相そのもの)、
    // CPG (thermalMethod != 2) は種 DB に H2O が無いので共通データの内蔵気相を使う。凝縮 OFF・N2 では何もしない (記録はバイト不変)。
    if (cfg.condensation == 1 && cfg.condModel == 1) {
        speciesDB_attachCondensed(db, "H2O(L)", cfg.thermalMethod == 2 ? cfg.condGasSpeciesName : std::string(), cfg.thermalMethod == 2);
    }
    // 種ごとの輸送物性の出所 (physProp.transport; plan #5t2 段 1)。書かれていなければ何もしない (記録・ハッシュは従来のまま)。
    if (!cfg.speciesTransport.empty()) speciesTransportDB_resolve(db, cfg.speciesTransport, cfg.speciesDBFile);
    return db;
}

namespace {

// 熱物性 (MW・区間・係数; withLJ なら LJ も) の差を従来の順序 (MW, 区間, LJ, 係数) で emit(key, a の値, b の値) に渡す。
// 2 区間同士は #13-1 前と同じキー名 (Tlo/Tmid/Thi, nasa9_low[k]/nasa9_high[k])、それ以外は nInt / Tbounds[k] / nasa9_intervals[j][k]。
template <class Emit>
void thermoDiff(const SpeciesThermo& a, const SpeciesThermo& b, bool withLJ, Emit emit)
{
    auto d = [&](const std::string& key, double x, double y) { if (x != y) emit(key, g17(x), g17(y)); };
    const bool two = (a.nInt == 2 && b.nInt == 2);
    d("MW", a.MW, b.MW);
    if (two) {
        d("Tlo", a.Tlo, b.Tlo); d("Tmid", a.Tbrk[0], b.Tbrk[0]); d("Thi", a.Thi, b.Thi);
    } else if (a.nInt != b.nInt) {
        emit("nInt", std::to_string(a.nInt) + " (" + breakpoints(a) + " K)", std::to_string(b.nInt) + " (" + breakpoints(b) + " K)");
    } else {
        for (int k = 0; k <= a.nInt; ++k) d("Tbounds[" + std::to_string(k) + "]", thermo_bound(a, k), thermo_bound(b, k));
    }
    if (withLJ) { d("LJ_sigma", a.sigma_LJ, b.sigma_LJ); d("LJ_eps_kB", a.eps_kB, b.eps_kB); }
    if (two) {
        for (int k = 0; k < 9; ++k) d("nasa9_low[" + std::to_string(k) + "]", a.coef[0][k], b.coef[0][k]);
        for (int k = 0; k < 9; ++k) d("nasa9_high[" + std::to_string(k) + "]", a.coef[1][k], b.coef[1][k]);
    } else if (a.nInt == b.nInt) {
        for (int j = 0; j < a.nInt; ++j)
            for (int k = 0; k < 9; ++k)
                d("nasa9_intervals[" + std::to_string(j) + "][" + std::to_string(k) + "]", a.coef[j][k], b.coef[j][k]);
    }
}

// 内蔵種の canonical ID (別名・大小文字違いを同一視)。内蔵に無ければ大文字化した名前。
std::string canonicalKeyCI(const std::string& n)
{
    const std::string up = toUpper(n);
    for (const auto& e : builtinEntries()) {
        if (toUpper(e.id) == up) return e.id;
        for (const auto& a : e.aliases) if (toUpper(a) == up) return e.id;
    }
    return up;
}

// 気液ペアの基準契約: 熱力学に効く値 (MW・区間・全区間の 9 係数) の最初の差 (空 = 一致)。LJ・datum・invMW は見ない。
std::string firstThermoDiff(const SpeciesThermo& a, const SpeciesThermo& b)
{
    std::string r;
    thermoDiff(a, b, false, [&](const std::string& key, const std::string& x, const std::string& y) {
        if (r.empty()) r = key + " " + x + " vs built-in " + y;
    });
    return r;
}

} // anonymous namespace

void speciesDB_attachCondensed(ResolvedSpeciesDB& db, const std::string& id, const std::string& gasName, bool requireInList)
{
    const std::string where = "[speciesDB] condensed phase '" + id + "'";
    const CondensedEntry* ce = nullptr;
    for (const auto& c : condensedEntries()) if (c.id == id) ce = &c;
    if (!ce) throw std::runtime_error(where + " is not in the built-in species data " + std::string(kForgeSpeciesDataName));
    const BuiltinEntry* ge = nullptr;
    for (const auto& e : builtinEntries()) if (e.id == ce->pairOf) ge = &e;
    if (!ge) throw std::runtime_error(where + ": pair_of '" + ce->pairOf + "' is not a built-in gas species");
    if (ce->MW != ge->sp.MW) {
        throw std::runtime_error(where + ": MW " + g17(ce->MW) + " differs from its gas pair '" + ce->pairOf + "' MW " + g17(ge->sp.MW)
                                 + " (the pair is converted to mass with the same MW)");
    }
    ResolvedCondensed c;
    c.enabled = true;
    c.name = ce->id; c.pairOf = ce->pairOf;
    c.MW = ce->MW; c.Tlo = ce->Tlo; c.Thi = ce->Thi;
    for (int k = 0; k < 9; ++k) c.coeffs[k] = ce->coeffs[k];
    c.below = ce->below; c.above = ce->above;
    c.gasIndex = -1; c.gasName = ce->pairOf; c.gas = ge->sp;
    if (requireInList) {
        int gi = -1;
        if (!gasName.empty()) {
            gi = db.index(gasName);
            if (gi < 0) throw std::runtime_error(where + ": condensing gas species '" + gasName + "' is not in physProp.species");
        } else {
            for (int s = 0; s < db.size() && gi < 0; ++s) {
                if (!db.isLump(s) && canonicalKeyCI(db.dbKey[s]) == ce->pairOf) gi = s;
            }
            if (gi < 0) {
                throw std::runtime_error(where + ": its gas pair '" + ce->pairOf + "' is not in physProp.species (thermalMethod 2 evaluates the "
                                         "latent heat from the species DB entry of the condensing gas; list it as a species)");
            }
        }
        if (db.isLump(gi)) throw std::runtime_error(where + ": the condensing gas '" + db.names[gi] + "' is a lump");
        const bool file = (db.source[gi] == "file");
        if (canonicalKeyCI(db.dbKey[gi]) != ce->pairOf) {
            throw std::runtime_error(where + ": condensing gas species '" + db.names[gi] + "' (DB key '" + db.dbKey[gi]
                                     + "') is not its gas pair '" + ce->pairOf + "'");
        }
        // 気液ペアの基準契約 (§4.8, 2 回目 M2): 液相 H2O(L) は内蔵の気相 H2O と同じ絶対基準 (CEA) のペア。外部 DB が気相を
        // 上書きしていても値が内蔵と完全に同じなら整合を確かめられるので通し、違えば拒否する (L がその差だけずれるため)。
        const std::string diff = firstThermoDiff(db.species[gi], ge->sp);
        if (!diff.empty()) {
            throw std::runtime_error(where + ": the condensing gas '" + db.names[gi] + "' comes from " + (file ? "speciesDBFile" : db.source[gi])
                                     + " and differs from the built-in gas '" + ce->pairOf + "' that the liquid phase pairs with (" + diff
                                     + "). The latent heat L = h_v - h_l needs the gas and liquid on the same absolute (CEA) basis; "
                                     "a gas entry that cannot be checked against the pair is refused with condensation ON "
                                     "(remove the gas override from speciesDBFile, or give it exactly the built-in coefficients and MW).");
        }
        c.gasIndex = gi; c.gasName = db.names[gi]; c.gas = db.species[gi];
    }
    db.condensed = c;
}

const ResolvedSpeciesDB& speciesDB_init(const solverConfig& cfg)
{
    try {
        g_current = std::make_unique<ResolvedSpeciesDB>(speciesDB_resolve(cfg));
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return *g_current;
}

const ResolvedSpeciesDB* speciesDB_current()
{
    return g_current.get();
}

std::vector<double> speciesMoleToMass(const std::vector<double>& X, const std::vector<double>& MW)
{
    if (X.size() != MW.size()) throw std::runtime_error("[speciesDB] mole->mass: X and MW length mismatch");
    double denom = 0.0;
    for (size_t k = 0; k < X.size(); ++k) denom += X[k]*MW[k];
    if (!(denom > 0.0) || !std::isfinite(denom)) throw std::runtime_error("[speciesDB] mole->mass: sum(X_j M_j) must be positive");
    std::vector<double> Y(X.size());
    for (size_t k = 0; k < X.size(); ++k) Y[k] = X[k]*MW[k]/denom;
    return Y;
}

std::vector<double> speciesMassToMole(const std::vector<double>& Y, const std::vector<double>& MW)
{
    if (Y.size() != MW.size()) throw std::runtime_error("[speciesDB] mass->mole: Y and MW length mismatch");
    double denom = 0.0;
    for (size_t k = 0; k < Y.size(); ++k) denom += Y[k]/MW[k];
    if (!(denom > 0.0) || !std::isfinite(denom)) throw std::runtime_error("[speciesDB] mass->mole: sum(Y_j/M_j) must be positive");
    std::vector<double> X(Y.size());
    for (size_t k = 0; k < Y.size(); ++k) X[k] = (Y[k]/MW[k])/denom;
    return X;
}

std::vector<double> bcondSpeciesMassFractions(const YAML::Node& floats, const ResolvedSpeciesDB& db,
                                              const std::string& bname)
{
    const int n = db.size();
    std::map<int, double> X, Y;
    if (floats && floats.IsMap()) {
        for (auto it = floats.begin(); it != floats.end(); ++it) {
            const std::string key = it->first.as<std::string>();
            const int ix = fractionIndex(key, 'X');
            const int iy = fractionIndex(key, 'Y');
            if (ix < 0 && iy < 0) continue;
            double v;
            try {
                v = it->second.as<double>();   // double で読む (flow_float 化は換算後)
            } catch (const std::exception&) {
                throw std::runtime_error("boundary '" + bname + "': floats." + key + " is not a number");
            }
            if (!std::isfinite(v)) throw std::runtime_error("boundary '" + bname + "': floats." + key + " is not finite");
            if (v < 0.0) throw std::runtime_error("boundary '" + bname + "': floats." + key + " is negative (" + std::to_string(v) + ")");
            const int idx = (ix >= 0) ? ix : iy;
            if (idx >= n) {
                throw std::runtime_error("boundary '" + bname + "': floats." + key + " refers to species index " + std::to_string(idx)
                                         + " but physProp.species has only " + std::to_string(n) + " entries");
            }
            if (ix >= 0) X[ix] = v; else Y[iy] = v;
        }
    }
    if (!X.empty() && !Y.empty()) {
        throw std::runtime_error("boundary '" + bname + "': X{s} (mole fractions) and Y{s} (mass fractions) must not be mixed on one boundary");
    }
    std::vector<double> out;
    if (!X.empty()) {
        // X 指定時は全種必須 (既定補完しない)。
        std::vector<double> Xv(n), MW(n);
        for (int s = 0; s < n; ++s) {
            auto it = X.find(s);
            if (it == X.end()) {
                throw std::runtime_error("boundary '" + bname + "': mole fractions require every species; X" + std::to_string(s)
                                         + " (" + db.names[s] + ") is missing");
            }
            Xv[s] = it->second; MW[s] = db.MW(s);
        }
        double sum = 0.0;
        for (double v : Xv) sum += v;
        if (!(sum > 0.0)) throw std::runtime_error("boundary '" + bname + "': sum of X{s} must be positive");
        out = speciesMoleToMass(Xv, MW);   // 内部で Σ X M > 0 も検査
    } else if (!Y.empty()) {
        out.assign(n, 0.0);
        for (int s = 0; s < n; ++s) {
            auto it = Y.find(s);
            out[s] = (it != Y.end()) ? it->second : ((s == 0) ? 1.0 : 0.0);   // 従来の既定補完
        }
        double sum = 0.0;
        for (double v : out) sum += v;
        if (std::fabs(sum - 1.0) > 1.0e-3) {
            std::ostringstream oss;
            oss << "boundary '" << bname << "': sum of Y{s} = " << std::setprecision(10) << sum
                << " differs from 1 by more than 1e-3 (unspecified species default to Y0=1, others 0; give every Y{s} or use X{s})";
            throw std::runtime_error(oss.str());
        }
    }
    return out;
}

void speciesDB_printTable(const solverConfig& cfg, const ResolvedSpeciesDB& db)
{
    std::cout << "[species] table (physProp.species order; MW [kg/mol]; source builtin|file"
              << (cfg.speciesDBFile.empty() ? "" : " '" + cfg.speciesDBFile + "'") << ")\n";
    std::cout << "[species]   built-in data: " << speciesDB_builtinDataName() << " (embedded, sha256 "
              << speciesDB_builtinDataSha256().substr(0, 16) << ")\n";
    for (int s = 0; s < db.size(); ++s) {
        std::cout << "[species]   " << std::setw(2) << s << "  " << std::setw(8) << std::left << db.names[s] << std::right
                  << "  MW=" << std::setprecision(10) << db.MW(s) << "  " << db.source[s] << "\n";
    }
    // lump (擬似種) の中身と合成結果 (h は datum 適用前の絶対基準; thermoHrefTemp>0 なら thermo_init_db が全種同じ規約で移す)
    for (int s = 0; s < db.size(); ++s) {
        if (!db.isLump(s)) continue;
        const ResolvedLump& l = db.lumps[s];
        const SpeciesThermo& sp = db.species[s];
        std::cout << "[species]   lump " << db.names[s] << " (basis " << l.basis << "; synthesized at startup: " << l.synthesis << ")\n";
        for (size_t k = 0; k < l.members.size(); ++k) {
            std::cout << "[species]     " << std::setw(8) << std::left << l.members[k] << std::right
                      << "  input=" << std::setprecision(17) << l.input[k]
                      << "  x=" << std::setprecision(17) << l.x[k]
                      << "  MW=" << std::setprecision(10) << l.memberSpecies[k].MW << "  " << l.memberSource[k] << "\n";
        }
        std::cout << "[species]     MW=" << std::setprecision(17) << sp.MW << " kg/mol, T breakpoints "
                  << breakpoints(sp) << " K\n";
        std::cout << "[species]     LJ sigma=" << std::setprecision(17) << sp.sigma_LJ << " A eps/kB=" << sp.eps_kB
                  << " K (PROVISIONAL: mass-fraction mean of constituents, as the design-side generator; to be replaced by"
                  << " constituent expansion / Blanc diffusion in plan #7)\n";
        for (double T : {298.15, 1000.0, 2000.0}) {
            std::cout << "[species]     T=" << std::setprecision(6) << T << " K: cp=" << std::setprecision(10) << thermo_cp_mass(sp, T)
                      << " J/(kg K), h_abs=" << thermo_h_mass(sp, T) << " J/kg\n";
        }
    }
    if (cfg.condensation == 1) {
        if (cfg.condGasSpecies >= 0 && cfg.condGasSpecies < db.size()) {
            std::cout << "[species]   condensing species: " << db.names[cfg.condGasSpecies]
                      << " (condGasSpecies=" << cfg.condGasSpecies << ")\n";
        } else {
            std::cout << "[species]   condensing species: pure condensible (condGasSpecies=-1)\n";
        }
        if (db.condensed.enabled) {
            const ResolvedCondensed& c = db.condensed;
            std::cout << "[species]   liquid phase: " << c.name << " (" << speciesDB_builtinDataName() << ", " << g17(c.Tlo) << "-" << g17(c.Thi)
                      << " K; below " << c.below << ", above " << c.above << ") paired with gas '" << c.gasName << "'"
                      << (c.gasIndex >= 0 ? " (species " + std::to_string(c.gasIndex) + ")" : std::string(" (built-in; not in the species list)"))
                      << ": latent heat L = h_v - h_l with the same datum and MW (plan #10)\n";
        }
    }
    if (db.transport.enabled) {
        // 種ごとの輸送物性の出所 (physProp.transport; plan #5t2 段 1)。GPU はまだ現行の輸送経路を使う (段 2 で接続)。
        const ResolvedTransport& tr = db.transport;
        std::cout << "[species]   transport (physProp.transport; " << TRANSPORT_RECORD_SCHEMA << ", data "
                  << speciesTransportDB_dataName() << " sha256 " << speciesTransportDB_dataSha256().substr(0, 16)
                  << ", trans.inp sha256 " << speciesTransportDB_transInpSha256().substr(0, 16) << ")\n";
        std::cout << "[species]     mixing: " << TRANSPORT_MIXING_RULE << "\n";
        for (int r = 0; r < tr.nReal(); ++r) {
            double m3, l3, m10, l10;
            transport_species(tr.sp[r], 300.0, &m3, &l3);
            transport_species(tr.sp[r], 1000.0, &m10, &l10);
            std::cout << "[species]     " << std::setw(2) << r << "  " << std::setw(8) << std::left << tr.realName[r] << std::right
                      << "  " << tr.modelName[r] << "  (" << tr.dataSource[r] << ")  mu(300K)=" << std::setprecision(6) << m3
                      << " Pa s, lambda(300K)=" << l3 << " W/(m K), mu(1000K)=" << m10 << ", lambda(1000K)=" << l10 << "\n";
        }
        for (int a = 0; a < tr.nReal(); ++a)
            for (int b = a + 1; b < tr.nReal(); ++b)
                std::cout << "[species]     eta_ij " << tr.realName[a] << "-" << tr.realName[b] << ": "
                          << tr.pairSource[transport_pair_index(a, b, tr.nReal())] << "\n";
        for (const auto& n : tr.notes) std::cout << "[species]     NOTE: " << n << "\n";
    }
    std::cout << "[species]   tracer: " << (cfg.tracerEnabled() ? cfg.tracer + " (roXi transported; inlet floats Xi)" : "none") << "\n";
}

// =============================================================================
// 解決済み記録と内容照合 (plans/active/thermophysics-solver-owned-species-db.md §4.3, #3a)
// =============================================================================
namespace {

// ---- SHA-256 (FIPS 180-4)。外部ライブラリに依存しない (変換器・単体試験も同じ TU を単独でビルドする) ----
struct Sha256 {
    uint32_t h[8];
    unsigned char buf[64];
    uint64_t len = 0;
    size_t   nbuf = 0;
    static uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
    Sha256() {
        const uint32_t h0[8] = {0x6a09e667u,0xbb67ae85u,0x3c6ef372u,0xa54ff53au,0x510e527fu,0x9b05688cu,0x1f83d9abu,0x5be0cd19u};
        for (int i = 0; i < 8; ++i) h[i] = h0[i];
    }
    void block(const unsigned char* p) {
        static const uint32_t k[64] = {
            0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
            0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
            0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
            0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
            0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
            0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
            0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
            0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u};
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = (uint32_t(p[4*i]) << 24) | (uint32_t(p[4*i+1]) << 16) | (uint32_t(p[4*i+2]) << 8) | uint32_t(p[4*i+3]);
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^ (w[i-15] >> 3);
            const uint32_t s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^ (w[i-2] >> 10);
            w[i] = w[i-16] + s0 + w[i-7] + s1;
        }
        uint32_t a=h[0], b=h[1], c=h[2], d=h[3], e=h[4], f=h[5], g=h[6], hh=h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const uint32_t ch = (e & f) ^ (~e & g);
            const uint32_t t1 = hh + S1 + ch + k[i] + w[i];
            const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t t2 = S0 + mj;
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=hh;
    }
    void update(const unsigned char* p, size_t n) {
        len += n;
        while (n > 0) {
            const size_t t = std::min(n, size_t(64) - nbuf);
            std::copy(p, p + t, buf + nbuf);
            nbuf += t; p += t; n -= t;
            if (nbuf == 64) { block(buf); nbuf = 0; }
        }
    }
    std::string hex() {
        const uint64_t bits = len * 8;
        const unsigned char pad = 0x80;
        update(&pad, 1);
        const unsigned char zero = 0;
        while (nbuf != 56) update(&zero, 1);
        unsigned char lb[8];
        for (int i = 0; i < 8; ++i) lb[i] = static_cast<unsigned char>(bits >> (56 - 8*i));
        update(lb, 8);
        static const char* hx = "0123456789abcdef";
        std::string out;
        for (int i = 0; i < 8; ++i)
            for (int j = 28; j >= 0; j -= 4) out.push_back(hx[(h[i] >> j) & 0xf]);
        return out;
    }
};

// %.17g (C の printf と Python の '%.17g' % x は同じ文字列を返す; 往復で double がビット一致する)
std::string g17(double v)
{
    char b[64];
    std::snprintf(b, sizeof(b), "%.17g", v);
    return b;
}

// YAML の二重引用符文字列 (\ と " をエスケープ)
std::string yq(const std::string& s)
{
    std::string o = "\"";
    for (char c : s) {
        if (c == '\\' || c == '"') o.push_back('\\');
        o.push_back(c);
    }
    return o + "\"";
}

// lump の構成種 1 つ (互換性テキストの材料)
struct RecordLumpMember {
    std::string   name;
    double        x;      // lump 内モル分率 (正規化済み)
    SpeciesThermo sp;     // 構成種の絶対基準係数
};

// 記録 1 種分 (互換性テキストの材料)
struct RecordEntry {
    std::string   name;
    std::string   phase;
    SpeciesThermo sp;
    bool          isLump = false;
    std::string   synthesis;                  // lump の合成規約 (SPECIES_LUMP_SYNTHESIS)
    std::vector<RecordLumpMember> members;    // lump の構成種
};

// nInt == 2 の種は #13-1 前と一字一句同じ行 (T: Tlo Tmid Thi / low / high)。それ以外は T: に全境界、区間ごとに coef[k] の行
// (その記録は schema・外挿規約の文字列も別値になる: recordSchemaOf / recordExtrapolationOf)。
void compatLineSpecies(std::ostringstream& o, const std::string& tag, const SpeciesThermo& s)
{
    o << tag << ".MW: " << g17(s.MW) << "\n";
    o << tag << ".T:";
    for (int k = 0; k <= s.nInt; ++k) o << " " << g17(thermo_bound(s, k));
    o << "\n";
    o << tag << ".LJ: " << g17(s.sigma_LJ) << " " << g17(s.eps_kB) << "\n";
    if (s.nInt == 2) {
        o << tag << ".low:";
        for (int k = 0; k < 9; ++k) o << " " << g17(s.coef[0][k]);
        o << "\n";
        o << tag << ".high:";
        for (int k = 0; k < 9; ++k) o << " " << g17(s.coef[1][k]);
        o << "\n";
        return;
    }
    for (int j = 0; j < s.nInt; ++j) {
        o << tag << ".coef[" << j << "]:";
        for (int k = 0; k < 9; ++k) o << " " << g17(s.coef[j][k]);
        o << "\n";
    }
}

// 互換性テキスト本体。書式は tools/forge_species.py compat_text と一字一句同じにすること。
// transportLines: 輸送ブロック (physProp.transport を書いた run だけ; plan #5t2)。空なら何も足さない (従来とバイト一致)。
// condensedLines: 凝縮種の液相 (凝縮 ON・H2O の run だけ; plan #10)。空なら何も足さない (従来とバイト一致)。種の後・輸送の前に置く。
std::string compatTextRaw(const std::string& schema, const std::string& datum, const std::string& extrap,
                          double Tref, const std::vector<RecordEntry>& es,
                          const std::vector<std::string>& transportLines = std::vector<std::string>{},
                          const std::vector<std::string>& condensedLines = std::vector<std::string>{})
{
    std::ostringstream o;
    o << "schema: " << schema << "\n";
    o << "datum: " << datum << "\n";
    o << "thermoHrefTemp: " << g17(Tref) << "\n";
    o << "extrapolation: " << extrap << "\n";
    o << "nSpecies: " << es.size() << "\n";
    for (size_t i = 0; i < es.size(); ++i) {
        const std::string tag = "species[" + std::to_string(i) + "]";
        o << tag << ": name=" << es[i].name << " phase=" << es[i].phase << "\n";
        compatLineSpecies(o, tag, es[i].sp);
        // lump だけ追記する (lump の無い DB のテキストは #6a 前とバイト一致)。basis・入力の分率・source は入れない。
        if (es[i].isLump) {
            o << tag << ".lump: n=" << es[i].members.size() << " synthesis=" << es[i].synthesis << "\n";
            for (size_t k = 0; k < es[i].members.size(); ++k) {
                const std::string mt = tag + ".lump[" + std::to_string(k) + "]";
                o << mt << ": name=" << es[i].members[k].name << " x=" << g17(es[i].members[k].x) << "\n";
                compatLineSpecies(o, mt, es[i].members[k].sp);
            }
        }
    }
    for (const auto& l : condensedLines) o << l << "\n";
    for (const auto& l : transportLines) o << l << "\n";
    return o.str();
}

// 凝縮種の液相の互換性行 (plan #10)。書式は tools/forge_species.py condensed_compat_lines と一字一句同じにすること。
//   ペアの気相種の係数は species[gas_index] の行に既にあるので、ここは液相の値と規約 (延長・潜熱・datum) だけ。
//   rule/latent/datum は規約文字列 (現在のコードでは SPECIES_CONDENSED_*; 記録から再計算するときは記録に書いた文字列)。
std::vector<std::string> condensedCompatLines(const ResolvedCondensed& c,
                                              const std::string& rule = SPECIES_CONDENSED_EXTENSION,
                                              const std::string& latent = SPECIES_CONDENSED_LATENT,
                                              const std::string& datum = SPECIES_CONDENSED_DATUM)
{
    std::vector<std::string> v;
    if (!c.enabled) return v;
    const std::string t = "condensed[0]";
    v.push_back(t + ": name=" + c.name + " phase=condensed pair_of=" + c.pairOf + " gas_index=" + std::to_string(c.gasIndex));
    v.push_back(t + ".MW: " + g17(c.MW));
    v.push_back(t + ".T: " + g17(c.Tlo) + " " + g17(c.Thi));
    std::string co = t + ".coeffs:";
    for (int k = 0; k < 9; ++k) co += " " + g17(c.coeffs[k]);
    v.push_back(co);
    v.push_back(t + ".extension: below=" + c.below + " above=" + c.above);
    v.push_back(t + ".rule: " + rule);
    v.push_back(t + ".latent: " + latent);
    v.push_back(t + ".datum: " + datum);
    return v;
}

// 記録に nInt != 2 の種 (lump の構成種を含む) が 1 つでもあるか (plan #13-1: そのときだけ schema・外挿規約が別値)
bool hasNonTwoInterval(const ResolvedSpeciesDB& db)
{
    for (int s = 0; s < db.size(); ++s) {
        if (db.species[s].nInt != 2) return true;
        if (db.isLump(s)) for (const auto& m : db.lumps[s].memberSpecies) if (m.nInt != 2) return true;
    }
    return false;
}

// 記録のスキーマ名 (輸送ブロックがあれば v2; nInt != 2 の種があれば区間可変の別名)
std::string recordSchemaOf(const ResolvedSpeciesDB& db)
{
    const bool nint = hasNonTwoInterval(db);
    if (db.transport.enabled) return nint ? SPECIES_RECORD_SCHEMA_TRANSPORT_NINT : SPECIES_RECORD_SCHEMA_TRANSPORT;
    return nint ? SPECIES_RECORD_SCHEMA_NINT : SPECIES_RECORD_SCHEMA;
}

// 外挿・区間選択の規約文字列 (nInt != 2 の種があれば区間可変の規約)
const char* recordExtrapolationOf(const ResolvedSpeciesDB& db)
{
    return hasNonTwoInterval(db) ? SPECIES_RECORD_EXTRAPOLATION_NINT : SPECIES_RECORD_EXTRAPOLATION;
}

const std::vector<std::string>& transportLinesOf(const ResolvedSpeciesDB& db)
{
    static const std::vector<std::string> none;
    return db.transport.enabled ? db.transport.compatLines : none;
}

std::vector<RecordEntry> entriesOf(const ResolvedSpeciesDB& db)
{
    std::vector<RecordEntry> es;
    for (int s = 0; s < db.size(); ++s) {
        RecordEntry e;
        e.name = db.names[s]; e.phase = "gas"; e.sp = db.species[s];   // 液相は #10 まで含めない
        if (db.isLump(s)) {
            const ResolvedLump& l = db.lumps[s];
            e.isLump = true;
            e.synthesis = l.synthesis;
            for (size_t k = 0; k < l.members.size(); ++k) e.members.push_back({l.members[k], l.x[k], l.memberSpecies[k]});
        }
        es.push_back(e);
    }
    return es;
}

// 記録の係数ブロック (種・構成種で共通)。indent は行頭の空白。
// nInt == 2 は #13-1 前と同じキー (Tlo/Tmid/Thi, nasa9_low/nasa9_high)。それ以外は Tbounds と nasa9_intervals (外部 DB と同じキー)。
void recordSpeciesBlock(std::ostringstream& o, const std::string& indent, const SpeciesThermo& sp)
{
    o << indent << "MW: " << g17(sp.MW) << "\n";
    if (sp.nInt == 2) {
        o << indent << "Tlo: " << g17(sp.Tlo) << "\n";
        o << indent << "Tmid: " << g17(sp.Tbrk[0]) << "\n";
        o << indent << "Thi: " << g17(sp.Thi) << "\n";
    } else {
        o << indent << "Tbounds: [";
        for (int k = 0; k <= sp.nInt; ++k) o << (k ? ", " : "") << g17(thermo_bound(sp, k));
        o << "]\n";
    }
    o << indent << "LJ_sigma: " << g17(sp.sigma_LJ) << "\n";
    o << indent << "LJ_eps_kB: " << g17(sp.eps_kB) << "\n";
    if (sp.nInt == 2) {
        o << indent << "nasa9_low: [";
        for (int k = 0; k < 9; ++k) o << (k ? ", " : "") << g17(sp.coef[0][k]);
        o << "]\n";
        o << indent << "nasa9_high: [";
        for (int k = 0; k < 9; ++k) o << (k ? ", " : "") << g17(sp.coef[1][k]);
        o << "]\n";
    } else {
        o << indent << "nasa9_intervals:\n";
        for (int j = 0; j < sp.nInt; ++j) {
            o << indent << "  - [";
            for (int k = 0; k < 9; ++k) o << (k ? ", " : "") << g17(sp.coef[j][k]);
            o << "]\n";
        }
    }
}

// 記録の 1 種 (YAML ノード) から係数を読む (2 区間の書式と区間可変の書式の両方)
void readRecordSpecies(const YAML::Node& n, const std::string& what, SpeciesThermo& s)
{
    s = SpeciesThermo{};
    s.MW = n["MW"].as<double>();
    s.sigma_LJ = n["LJ_sigma"].as<double>(); s.eps_kB = n["LJ_eps_kB"].as<double>();
    if (n["Tbounds"]) {
        std::vector<double> Tb;
        for (const auto& v : n["Tbounds"]) Tb.push_back(v.as<double>());
        const YAML::Node iv = n["nasa9_intervals"];
        const int nInt = (iv && iv.IsSequence()) ? static_cast<int>(iv.size()) : 0;
        const std::string bad = intervalProblem(nInt, Tb);
        if (!bad.empty()) throw std::runtime_error("species '" + what + "' " + bad);
        double a[THERMO_MAX_INTERVALS][9];
        for (int k = 0; k < nInt; ++k)
            if (!read9(iv[k], a[k])) throw std::runtime_error("species '" + what + "' lacks 9 nasa9 coefficients in interval " + std::to_string(k));
        thermo_set_intervals(s, nInt, Tb.data(), a);
        return;
    }
    double lo[9], hi[9];
    if (!read9(n["nasa9_low"], lo) || !read9(n["nasa9_high"], hi))
        throw std::runtime_error("species '" + what + "' lacks 9 nasa9 coefficients");
    thermo_set_nasa9_2(s, n["Tlo"].as<double>(), n["Tmid"].as<double>(), n["Thi"].as<double>(), lo, hi);
}

bool readFileBytes(const std::string& path, std::string& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
}

std::unique_ptr<SpeciesRecordInfo> g_record;   // ソルバ起動時の記録 (出力属性用)

} // anonymous namespace

std::string speciesDB_sha256Hex(const std::string& bytes)
{
    Sha256 h;
    h.update(reinterpret_cast<const unsigned char*>(bytes.data()), bytes.size());
    return h.hex();
}

std::string speciesDB_compatText(const ResolvedSpeciesDB& db, double Tref)
{
    return compatTextRaw(recordSchemaOf(db), SPECIES_RECORD_DATUM, recordExtrapolationOf(db), Tref, entriesOf(db),
                         transportLinesOf(db), condensedCompatLines(db.condensed));
}

std::string speciesDB_compatHash(const ResolvedSpeciesDB& db, double Tref)
{
    return speciesDB_sha256Hex(speciesDB_compatText(db, Tref));
}

std::string speciesDB_recordText(const ResolvedSpeciesDB& db, double Tref, const std::string& dbFile,
                                 const std::string& inputField, const std::string& inputStatus)
{
    std::ostringstream o;
    o << "# forge resolved species record: OUTPUT (provenance), not an input.\n";
    o << "# plans/active/thermophysics-solver-owned-species-db.md 4.3 / methods/thermophysics.md 1b.4\n";
    o << "# compat_hash = SHA-256 of the canonical text rebuilt from schema/datum/thermoHrefTemp/extrapolation/species\n";
    o << "#   (source and provenance excluded; see tools/forge_species.py compat_text).\n";
    o << "# integrity   = SHA-256 of this whole file (res_*.h5 attribute species_record_sha256).\n";
    o << "schema: " << yq(recordSchemaOf(db)) << "\n";
    o << "compat_hash: " << yq(speciesDB_compatHash(db, Tref)) << "\n";
    o << "datum: " << yq(SPECIES_RECORD_DATUM) << "\n";
    o << "thermoHrefTemp: " << g17(Tref) << "\n";
    o << "extrapolation: " << yq(recordExtrapolationOf(db)) << "\n";
    o << "species:\n";
    for (int s = 0; s < db.size(); ++s) {
        const SpeciesThermo& sp = db.species[s];
        o << "  - index: " << s << "\n";
        o << "    name: " << yq(db.names[s]) << "\n";
        o << "    phase: \"gas\"\n";
        o << "    source: " << yq(db.source[s]) << "\n";
        recordSpeciesBlock(o, "    ", sp);
        if (db.isLump(s)) {
            // lump: 上の係数は起動時に合成した値。構成 (config の basis と分率、正規化したモル分率) と構成種の係数を書く。
            const ResolvedLump& l = db.lumps[s];
            o << "    lump:\n";
            o << "      basis: " << yq(l.basis) << "\n";
            o << "      synthesis: " << yq(l.synthesis) << "\n";
            o << "      members:\n";
            for (size_t k = 0; k < l.members.size(); ++k) {
                o << "        - name: " << yq(l.members[k]) << "\n";
                o << "          source: " << yq(l.memberSource[k]) << "\n";
                o << "          fraction_input: " << g17(l.input[k]) << "\n";
                o << "          x: " << g17(l.x[k]) << "\n";
                recordSpeciesBlock(o, "          ", l.memberSpecies[k]);
            }
        }
    }
    if (db.condensed.enabled) {
        // 凝縮種の液相 (plan #10): 係数は datum 前 (実行時はペアの気相と同じ Δa7 を a7 に足す)。規約文字列も書く
        // (記録から互換性テキストを作り直すときはここの文字列を使う)。source は来歴 (共通データの版は provenance に無いので名前だけ)。
        const ResolvedCondensed& c = db.condensed;
        o << "condensed:\n";
        o << "  - name: " << yq(c.name) << "\n";
        o << "    phase: \"condensed\"\n";
        o << "    source: " << yq("builtin " + speciesDB_builtinDataName()) << "\n";
        o << "    pair_of: " << yq(c.pairOf) << "\n";
        o << "    gas_index: " << c.gasIndex << "\n";
        o << "    gas_name: " << yq(c.gasName) << "\n";
        o << "    MW: " << g17(c.MW) << "\n";
        o << "    Tlo: " << g17(c.Tlo) << "\n";
        o << "    Thi: " << g17(c.Thi) << "\n";
        o << "    nasa9: [";
        for (int k = 0; k < 9; ++k) o << (k ? ", " : "") << g17(c.coeffs[k]);
        o << "]\n";
        o << "    extension: {below: " << yq(c.below) << ", above: " << yq(c.above) << "}\n";
        o << "    rule: " << yq(SPECIES_CONDENSED_EXTENSION) << "\n";
        o << "    latent: " << yq(SPECIES_CONDENSED_LATENT) << "\n";
        o << "    datum: " << yq(SPECIES_CONDENSED_DATUM) << "\n";
    }
    if (db.transport.enabled) {
        // 輸送ブロック (plan #5t2): 互換性ハッシュに入る正規化行をそのまま並べる (記録から同じテキストを再構成できる)。
        // 種ごとのモデルと係数・η_ij の出所と係数・混合則の版・H2O の接続規約・lump の展開行列を含む。
        o << "transport_compat:\n";
        for (const auto& l : db.transport.compatLines) o << "  - " << yq(l) << "\n";
    }
    o << "provenance:\n";
    o << "  speciesDBFile: " << yq(dbFile) << "\n";
    o << "  input_field: " << yq(inputField) << "\n";
    o << "  input_status: " << yq(inputStatus) << "\n";
    if (db.transport.enabled) {
        // 輸送の来歴 (互換性ハッシュには入れない): 埋め込みデータの版と、実種・組ごとのデータの出典
        const ResolvedTransport& tr = db.transport;
        o << "  transport_data: " << yq(speciesTransportDB_dataName() + " sha256 " + speciesTransportDB_dataSha256()
                                        + "; CEA trans.inp sha256 " + speciesTransportDB_transInpSha256()) << "\n";
        o << "  transport_sources:\n";
        for (int r = 0; r < tr.nReal(); ++r)
            o << "    - " << yq("real[" + std::to_string(r) + "] " + tr.realName[r] + " (physProp.transport key '" + tr.realConfigKey[r]
                                + "'): " + tr.modelName[r] + " -- " + tr.dataSource[r]) << "\n";
        for (int a = 0; a < tr.nReal(); ++a)
            for (int b = a + 1; b < tr.nReal(); ++b)
                o << "    - " << yq("pair[" + std::to_string(a) + "," + std::to_string(b) + "] " + tr.realName[a] + "-" + tr.realName[b]
                                    + ": " + tr.pairSource[transport_pair_index(a, b, tr.nReal())]) << "\n";
    }
    return o.str();
}

SpeciesRecordInfo speciesDB_writeRecord(const ResolvedSpeciesDB& db, double Tref, const std::string& dbFile,
                                        const std::string& inputField, const std::string& inputStatus,
                                        int inputUnverified, const std::string& dir)
{
    namespace fs = std::filesystem;
    SpeciesRecordInfo r;
    const std::string text = speciesDB_recordText(db, Tref, dbFile, inputField, inputStatus);
    r.compatHash = speciesDB_compatHash(db, Tref);
    r.recordSha256 = speciesDB_sha256Hex(text);
    r.inputUnverified = inputUnverified;
    const std::string base  = "resolved_species_" + r.compatHash.substr(0, 16);
    // 同名があって全文一致ならそのまま、違えば完全性ハッシュ付きの別名 (既存は上書きしない)
    for (const std::string& name : {base + ".yaml", base + "_" + r.recordSha256.substr(0, 16) + ".yaml"}) {
        const fs::path p = fs::path(dir) / name;
        std::string old;
        if (fs::exists(p)) {
            if (readFileBytes(p.string(), old) && old == text) { r.recordFile = name; return r; }
            if (name == base + ".yaml") {
                // 名前は同じ互換ハッシュでも来歴 (input_status 等) が違う正常な場合と、取り違えの場合がある。後者は警告。
                const std::string oldSha = speciesDB_sha256Hex(old);
                std::vector<std::string> d = speciesDB_diffRecord(p.string(), db, Tref);
                if (!d.empty()) {
                    std::cerr << "[species] WARNING: existing record " << p.string() << " (sha256 " << oldSha.substr(0, 16)
                              << ") does not match its name / the current species (" << d.front() << "); writing a separate record\n";
                }
                continue;
            }
            throw std::runtime_error("[species] record " + p.string() + " exists with different content (integrity collision?)");
        }
        std::ofstream f(p, std::ios::binary);
        if (!f) throw std::runtime_error("[species] cannot write record " + p.string());
        f << text;
        f.close();
        if (!f) throw std::runtime_error("[species] failed to write record " + p.string());
        r.recordFile = name;
        return r;
    }
    throw std::runtime_error("[species] could not place record " + base);
}

std::vector<std::string> speciesDB_diffRecord(const std::string& recordPath, const ResolvedSpeciesDB& db, double Tref)
{
    std::vector<std::string> d;
    YAML::Node root;
    try {
        root = YAML::LoadFile(recordPath);
    } catch (const std::exception& e) {
        d.push_back("record " + recordPath + " is not readable YAML: " + e.what());
        return d;
    }
    std::vector<RecordEntry> es;
    std::string schema, datum, extrap, hashInFile;
    double TrefRec = 0.0;
    std::vector<std::string> trLines;   // 輸送ブロック (v2 の記録だけ)
    ResolvedCondensed cRec;             // 凝縮種の液相 (記録にあるときだけ enabled; plan #10)
    std::string cRule, cLatent, cDatum;
    try {
        schema = root["schema"].as<std::string>();
        datum  = root["datum"].as<std::string>();
        extrap = root["extrapolation"].as<std::string>();
        hashInFile = root["compat_hash"].as<std::string>();
        TrefRec = root["thermoHrefTemp"].as<double>();
        for (const auto& n : root["species"]) {
            RecordEntry e;
            e.name  = n["name"].as<std::string>();
            e.phase = n["phase"].as<std::string>();
            readRecordSpecies(n, e.name, e.sp);
            if (n["lump"]) {
                const YAML::Node l = n["lump"];
                e.isLump = true;
                e.synthesis = l["synthesis"].as<std::string>();
                for (const auto& m : l["members"]) {
                    RecordLumpMember mm;
                    mm.name = m["name"].as<std::string>();
                    mm.x = m["x"].as<double>();
                    readRecordSpecies(m, e.name + "." + mm.name, mm.sp);
                    e.members.push_back(mm);
                }
            }
            es.push_back(e);
        }
        if (root["transport_compat"]) for (const auto& l : root["transport_compat"]) trLines.push_back(l.as<std::string>());
        if (root["condensed"]) {
            // 凝縮種の液相 (plan #10)。1 種だけ (H2O(L))。
            const YAML::Node cn = root["condensed"];
            if (!cn.IsSequence() || cn.size() != 1) throw std::runtime_error("condensed must be a list of one entry");
            const YAML::Node n = cn[0];
            cRec.enabled = true;
            cRec.name = n["name"].as<std::string>();
            cRec.pairOf = n["pair_of"].as<std::string>();
            cRec.gasIndex = n["gas_index"].as<int>();
            cRec.gasName = n["gas_name"].as<std::string>();
            cRec.MW = n["MW"].as<double>();
            cRec.Tlo = n["Tlo"].as<double>(); cRec.Thi = n["Thi"].as<double>();
            if (!read9(n["nasa9"], cRec.coeffs)) throw std::runtime_error("condensed '" + cRec.name + "' lacks 9 nasa9 coefficients");
            cRec.below = n["extension"]["below"].as<std::string>();
            cRec.above = n["extension"]["above"].as<std::string>();
            cRule = n["rule"].as<std::string>(); cLatent = n["latent"].as<std::string>(); cDatum = n["datum"].as<std::string>();
        }
    } catch (const std::exception& e) {
        d.push_back("record " + recordPath + " is malformed: " + e.what());
        return d;
    }
    // 記録の自己整合: 中身から互換ハッシュを作り直して記録内の値と比べる (改竄・取り違えの検出)
    const std::string rehash = speciesDB_sha256Hex(compatTextRaw(schema, datum, extrap, TrefRec, es, trLines,
                                                                 condensedCompatLines(cRec, cRule, cLatent, cDatum)));
    if (rehash != hashInFile) {
        d.push_back("record " + recordPath + ": compat_hash in file " + hashInFile.substr(0, 16)
                    + " != recomputed from its content " + rehash.substr(0, 16) + " (edited or corrupted record)");
    }
    if (schema != recordSchemaOf(db)) d.push_back("schema: field '" + schema + "' vs current '" + recordSchemaOf(db) + "'");
    {
        // 輸送ブロック (plan #5t2): 行単位で比べ、最初の差を示す
        const std::vector<std::string>& cur = transportLinesOf(db);
        if (trLines != cur) {
            size_t k = 0;
            while (k < trLines.size() && k < cur.size() && trLines[k] == cur[k]) ++k;
            d.push_back("transport: field '" + (k < trLines.size() ? trLines[k] : std::string("(end)")) + "' vs current '"
                        + (k < cur.size() ? cur[k] : std::string("(end)")) + "'");
        }
    }
    if (datum != SPECIES_RECORD_DATUM) d.push_back("datum convention: field '" + datum + "' vs current");
    if (extrap != recordExtrapolationOf(db)) d.push_back("extrapolation: field '" + extrap + "' vs current '" + recordExtrapolationOf(db) + "'");
    if (TrefRec != Tref) d.push_back("thermoHrefTemp: field " + g17(TrefRec) + " vs current " + g17(Tref));
    {
        // 凝縮種の液相 (plan #10): 有無・名前・ペア・MW・区間・係数・延長規約・規約文字列
        const ResolvedCondensed& cc = db.condensed;
        if (cRec.enabled != cc.enabled) {
            d.push_back(std::string("condensed (liquid phase for the latent heat): field ")
                        + (cRec.enabled ? "has '" + cRec.name + "'" : "has none (record written before plan #10, or condensation OFF: the latent-heat model of the field is not recorded)")
                        + " vs current " + (cc.enabled ? "'" + cc.name + "'" : "none"));
        } else if (cc.enabled) {
            const std::string t = "condensed " + cc.name;
            auto cs = [&](const std::string& key, const std::string& a, const std::string& b) {
                if (a != b) d.push_back(t + "." + key + ": field '" + a + "' vs current '" + b + "'");
            };
            auto cn = [&](const std::string& key, double a, double b) {
                if (a != b) d.push_back(t + "." + key + ": field " + g17(a) + " vs current " + g17(b));
            };
            cs("name", cRec.name, cc.name); cs("pair_of", cRec.pairOf, cc.pairOf);
            cs("gas_index", std::to_string(cRec.gasIndex), std::to_string(cc.gasIndex));
            cn("MW", cRec.MW, cc.MW); cn("Tlo", cRec.Tlo, cc.Tlo); cn("Thi", cRec.Thi, cc.Thi);
            for (int k = 0; k < 9; ++k) cn("nasa9[" + std::to_string(k) + "]", cRec.coeffs[k], cc.coeffs[k]);
            cs("extension.below", cRec.below, cc.below); cs("extension.above", cRec.above, cc.above);
            cs("rule", cRule, SPECIES_CONDENSED_EXTENSION); cs("latent", cLatent, SPECIES_CONDENSED_LATENT); cs("datum", cDatum, SPECIES_CONDENSED_DATUM);
        }
    }
    std::vector<std::string> nf, nc;
    for (const auto& e : es) nf.push_back(e.name);
    nc = db.names;
    if (nf != nc) {
        std::string a, b;
        for (const auto& x : nf) a += (a.empty() ? "" : ",") + x;
        for (const auto& x : nc) b += (b.empty() ? "" : ",") + x;
        d.push_back("species list/order: field [" + a + "] vs current [" + b + "]");
        return d;
    }
    const std::vector<RecordEntry> cur = entriesOf(db);
    auto cmpSpecies = [&](const std::string& n, const SpeciesThermo& f, const SpeciesThermo& c) {
        thermoDiff(f, c, true, [&](const std::string& key, const std::string& x, const std::string& y) {
            d.push_back(n + "." + key + ": field " + x + " vs current " + y);
        });
    };
    for (int s = 0; s < db.size(); ++s) {
        const std::string& n = db.names[s];
        if (es[s].phase != "gas") d.push_back(n + ".phase: field " + es[s].phase + " vs current gas");
        cmpSpecies(n, es[s].sp, db.species[s]);
        // lump の構成 (有無・合成規約・構成種の名前と順序・モル分率・構成種の係数)
        if (es[s].isLump != cur[s].isLump) {
            d.push_back(n + ": field is " + (es[s].isLump ? "a lump" : "not a lump") + " vs current " + (cur[s].isLump ? "a lump" : "not a lump"));
            continue;
        }
        if (!cur[s].isLump) continue;
        if (es[s].synthesis != cur[s].synthesis) d.push_back(n + ".lump.synthesis: field '" + es[s].synthesis + "' vs current '" + cur[s].synthesis + "'");
        std::string a, b;
        for (const auto& m : es[s].members) a += (a.empty() ? "" : ",") + m.name;
        for (const auto& m : cur[s].members) b += (b.empty() ? "" : ",") + m.name;
        if (a != b) { d.push_back(n + ".lump members: field [" + a + "] vs current [" + b + "]"); continue; }
        for (size_t k = 0; k < cur[s].members.size(); ++k) {
            const std::string mn = n + ".lump." + cur[s].members[k].name;
            if (es[s].members[k].x != cur[s].members[k].x)
                d.push_back(mn + ".x: field " + g17(es[s].members[k].x) + " vs current " + g17(cur[s].members[k].x));
            cmpSpecies(mn, es[s].members[k].sp, cur[s].members[k].sp);
        }
    }
    return d;
}

bool speciesDB_checkInputField(const ResolvedSpeciesDB& db, double Tref,
                               const std::string& fieldHash, const std::string& fieldRecordSha, int fieldUnverified,
                               const std::string& fieldPath, const std::vector<std::string>& searchDirs,
                               bool allowUnverified, std::string& inputStatus, int& inputUnverified, std::string& msg)
{
    namespace fs = std::filesystem;
    const std::string own = speciesDB_compatHash(db, Tref);
    msg.clear();
    if (fieldHash.empty()) {
        if (allowUnverified) {
            inputStatus = "unverified_env";
            inputUnverified = 1;
            msg = "input field '" + fieldPath + "' has no species_hash attribute (UNVERIFIED: cannot check which species "
                  "properties produced it). Allowed for this invocation by FORGE_ALLOW_UNVERIFIED_SPECIES=1; "
                  "outputs carry species_input_unverified=1.";
            return true;
        }
        msg = "input field '" + fieldPath + "' has no species_hash attribute: UNVERIFIABLE (the species properties that "
              "produced this field are unknown; fields written before species records existed, or initial fields not yet "
              "stamped by their generator). Refusing to start. Either:\n"
              "  (1) regenerate the initial field with a generator that stamps species attributes (it resolves this run with "
              "forge --resolve-species and writes species_hash; copies/restarts made with tools/restart_field.py or "
              "tools/interp_field.py from a stamped field inherit them), or\n"
              "  (2) if you have checked that the field was produced with the same species/datum as this run, allow it for THIS "
              "invocation only:  FORGE_ALLOW_UNVERIFIED_SPECIES=1 forge   (not a config key; outputs are marked "
              "species_input_unverified=1 and restarts from them inherit the mark).\n"
              "  If the species set/DB changed, convert the field with tools/convert_species_field.py instead.";
        return false;
    }
    if (fieldHash == own) {
        inputUnverified = (fieldUnverified == 1) ? 1 : 0;
        inputStatus = inputUnverified ? "unverified_inherited" : "verified";
        msg = "input field species_hash matches (" + own.substr(0, 16) + ")"
              + std::string(inputUnverified ? "; the input itself descends from an unverified start (mark inherited)" : "");
        return true;
    }
    std::ostringstream o;
    o << "input field '" << fieldPath << "' was produced with different species properties:\n"
      << "  field species_hash   " << fieldHash << "\n"
      << "  current species_hash " << own << "\n";
    // 入力側の記録を探して差を出す (完全性ハッシュが属性と一致するものだけ信用する)
    std::string found;
    std::vector<std::string> badIntegrity;
    const std::string prefix = "resolved_species_" + fieldHash.substr(0, std::min<size_t>(16, fieldHash.size()));
    for (const auto& dir : searchDirs) {
        std::error_code ec;
        if (!fs::is_directory(dir, ec)) continue;
        for (const auto& ent : fs::directory_iterator(dir, ec)) {
            const std::string fn = ent.path().filename().string();
            if (fn.rfind(prefix, 0) != 0 || ent.path().extension() != ".yaml") continue;
            std::string bytes;
            if (!readFileBytes(ent.path().string(), bytes)) continue;
            if (!fieldRecordSha.empty() && speciesDB_sha256Hex(bytes) == fieldRecordSha) { found = ent.path().string(); break; }
            badIntegrity.push_back(ent.path().string());
        }
        if (!found.empty()) break;
    }
    if (!found.empty()) {
        o << "  differences (field record " << found << " vs current):\n";
        for (const auto& x : speciesDB_diffRecord(found, db, Tref)) o << "    " << x << "\n";
    } else if (!badIntegrity.empty()) {
        o << "  record file(s) for the field hash exist but their integrity SHA-256 does not match the field attribute "
          << "species_record_sha256 (" << fieldRecordSha.substr(0, 16) << "): mixed-up or edited record; the coefficient differences cannot be identified:\n";
        for (const auto& x : badIntegrity) o << "    " << x << "\n";
    } else {
        o << "  no record file " << prefix << "*.yaml found next to the input field or in the run directory; "
          << "the coefficient differences cannot be identified.\n";
    }
    o << "  Refusing to start (a mismatch is never allowed by FORGE_ALLOW_UNVERIFIED_SPECIES). "
      << "Use the matching species DB / thermoHrefTemp, or convert the field with tools/convert_species_field.py.";
    msg = o.str();
    return false;
}

void speciesDB_setCurrentRecord(const SpeciesRecordInfo& rec)
{
    g_record = std::make_unique<SpeciesRecordInfo>(rec);
}

const SpeciesRecordInfo* speciesDB_currentRecord()
{
    return g_record.get();
}
