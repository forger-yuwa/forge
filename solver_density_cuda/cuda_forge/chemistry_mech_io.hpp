#pragma once
// =============================================================================
// chemistry_mech_io.hpp — 反応機構ファイル (Cantera YAML サブセット) → ReactionTable (host 専用)
//
//   対応: units {length, quantity, activation-energy}, phases[0].species (省略可),
//         reactions[]: equation, rate-constant {A,b,Ea}, type: three-body, efficiencies,
//         default-efficiency, duplicate。falloff (Lindemann/Troe) は Phase 2 で対応 (現状はエラー)。
//   化学種 index は流れ側 (cfg.speciesNames) の順序に合わせる。機構に現れる種が流れ側に無ければエラー、
//   流れ側にだけある種 (不活性) は許容する。
// =============================================================================

#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include "chemistry_d.cuh"

namespace chem_io {

inline std::string trim(const std::string& s)
{
    const auto a = s.find_first_not_of(" \t"); if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t");  return s.substr(a, b - a + 1);
}

// 化学種名の照合 (AR/Ar, HE/He の別名を吸収)
inline int findSpecies(const std::vector<std::string>& names, const std::string& q)
{
    auto up = [](std::string s){ for (auto& c : s) c = (char)toupper((unsigned char)c); return s; };
    for (size_t i = 0; i < names.size(); ++i) if (names[i] == q) return (int)i;
    for (size_t i = 0; i < names.size(); ++i) if (up(names[i]) == up(q)) return (int)i;
    return -1;
}

struct Side { std::vector<std::pair<int,int>> terms; bool M = false; bool falloffM = false; };

inline Side parseSide(const std::string& side, const std::vector<std::string>& names, const std::string& eq)
{
    Side out;
    std::string s = side;
    // "(+M)" は falloff 記法
    const auto fm = s.find("(+M)");
    if (fm != std::string::npos) { out.falloffM = true; s.erase(fm, 4); }
    std::stringstream ss(s); std::string tok;
    std::vector<std::string> parts; std::string cur;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '+' && (i + 1 >= s.size() || s[i+1] != '+')) { parts.push_back(cur); cur.clear(); }
        else cur += s[i];
    }
    parts.push_back(cur);
    for (auto& p : parts) {
        std::string t = trim(p);
        if (t.empty()) continue;
        if (t == "M") { out.M = true; continue; }
        int nu = 1;
        // 先頭の整数係数 ("2 OH")
        size_t k = 0; while (k < t.size() && isdigit((unsigned char)t[k])) ++k;
        if (k > 0 && k < t.size() && t[k] == ' ') { nu = std::stoi(t.substr(0, k)); t = trim(t.substr(k)); }
        const int idx = findSpecies(names, t);
        if (idx < 0) throw std::runtime_error("[chemistry] species '" + t + "' in reaction '" + eq + "' is not in physProp.species");
        bool merged = false;
        for (auto& e : out.terms) if (e.first == idx) { e.second += nu; merged = true; }
        if (!merged) out.terms.push_back({idx, nu});
    }
    return out;
}

// 単位換算係数: A [(length^3/quantity)^(order-1)/s] → [(m^3/mol)^(order-1)/s], Ea → [J/mol]
struct Units { double vol = 1.0e-6; /* cm^3 → m^3 */ double qty = 1.0; /* mol */ double Ea = 4.184; /* cal/mol → J/mol */ };

inline Units parseUnits(const YAML::Node& root)
{
    Units u;
    if (!root["units"]) return u;
    const YAML::Node un = root["units"];
    if (un["length"]) {
        const std::string L = un["length"].as<std::string>();
        if (L == "cm") u.vol = 1.0e-6; else if (L == "m") u.vol = 1.0;
        else throw std::runtime_error("[chemistry] unsupported units.length: " + L);
    }
    if (un["quantity"]) {
        const std::string Q = un["quantity"].as<std::string>();
        if (Q == "mol") u.qty = 1.0; else if (Q == "kmol") u.qty = 1.0e3;
        else throw std::runtime_error("[chemistry] unsupported units.quantity: " + Q);
    }
    if (un["activation-energy"]) {
        const std::string E = un["activation-energy"].as<std::string>();
        if (E == "cal/mol") u.Ea = 4.184; else if (E == "kcal/mol") u.Ea = 4184.0;
        else if (E == "J/mol") u.Ea = 1.0; else if (E == "kJ/mol") u.Ea = 1000.0;
        else if (E == "K") u.Ea = THERMO_RU;
        else throw std::runtime_error("[chemistry] unsupported units.activation-energy: " + E);
    }
    return u;
}

// 機構ファイルを読み ReactionTable を組む。names = 流れ側化学種名 (index 順)。
inline void loadMechanism(const std::string& path, const std::vector<std::string>& names, ReactionTable& rt,
                          std::vector<std::string>* equations = nullptr)
{
    YAML::Node root = YAML::LoadFile(path);
    const Units u = parseUnits(root);
    std::memset(&rt, 0, sizeof(ReactionTable));
    rt.nSpecies = (int)names.size();
    if (rt.nSpecies > THERMO_MAX_SPECIES) throw std::runtime_error("[chemistry] too many species");
    if (!root["reactions"]) throw std::runtime_error("[chemistry] no 'reactions' in " + path);

    int r = 0;
    for (const auto& rx : root["reactions"]) {
        if (r >= CHEM_MAX_REACTIONS) throw std::runtime_error("[chemistry] too many reactions (CHEM_MAX_REACTIONS)");
        const std::string eq = rx["equation"].as<std::string>();
        std::string type = rx["type"] ? rx["type"].as<std::string>() : "elementary";
        if (type == "falloff" || type == "chemically-activated" || type == "pressure-dependent-Arrhenius" || type == "Chebyshev")
            throw std::runtime_error("[chemistry] reaction type '" + type + "' not supported yet (Phase 2): " + eq);

        bool reversible = true; size_t pos; std::string lhs, rhs;
        if ((pos = eq.find("<=>")) != std::string::npos) { lhs = eq.substr(0, pos); rhs = eq.substr(pos + 3); }
        else if ((pos = eq.find("=>")) != std::string::npos) { lhs = eq.substr(0, pos); rhs = eq.substr(pos + 2); reversible = false; }
        else if ((pos = eq.find(" = ")) != std::string::npos) { lhs = eq.substr(0, pos); rhs = eq.substr(pos + 3); }
        else throw std::runtime_error("[chemistry] cannot parse equation: " + eq);

        const Side L = parseSide(lhs, names, eq), R = parseSide(rhs, names, eq);
        if (L.falloffM || R.falloffM) throw std::runtime_error("[chemistry] falloff '(+M)' not supported yet (Phase 2): " + eq);
        if ((int)L.terms.size() > CHEM_MAX_SIDE || (int)R.terms.size() > CHEM_MAX_SIDE)
            throw std::runtime_error("[chemistry] too many species on one side (CHEM_MAX_SIDE): " + eq);
        const bool thirdBody = (L.M || R.M || type == "three-body");

        rt.nR[r] = (int)L.terms.size(); rt.nP[r] = (int)R.terms.size();
        int order = thirdBody ? 1 : 0; double dnu = 0.0;
        for (int j = 0; j < rt.nR[r]; ++j) { rt.rIdx[r][j] = L.terms[j].first; rt.rNu[r][j] = L.terms[j].second; order += L.terms[j].second; dnu -= L.terms[j].second; }
        for (int j = 0; j < rt.nP[r]; ++j) { rt.pIdx[r][j] = R.terms[j].first; rt.pNu[r][j] = R.terms[j].second; dnu += R.terms[j].second; }
        rt.dnu[r] = dnu;
        rt.thirdBody[r] = thirdBody ? 1 : 0;
        rt.reversible[r] = reversible ? 1 : 0;

        const YAML::Node rc = rx["rate-constant"];
        if (!rc) throw std::runtime_error("[chemistry] 'rate-constant' missing: " + eq);
        const double A = rc["A"].as<double>();
        const double b = rc["b"] ? rc["b"].as<double>() : 0.0;
        const double Ea = rc["Ea"].as<double>();
        // A: (vol/qty)^(order-1)/s → (m^3/mol)^(order-1)/s
        rt.A[r]  = A * std::pow(u.vol / u.qty, order - 1);
        rt.b[r]  = b;
        rt.Ea[r] = Ea * u.Ea;

        const double defEff = rx["default-efficiency"] ? rx["default-efficiency"].as<double>() : 1.0;
        for (int s = 0; s < rt.nSpecies; ++s) rt.eff[r][s] = thirdBody ? defEff : 0.0;
        if (rx["efficiencies"]) {
            for (const auto& kv : rx["efficiencies"]) {
                const std::string nm = kv.first.as<std::string>();
                const int idx = findSpecies(names, nm);
                if (idx < 0) continue;    // 流れ側に無い種の効率は無視 (不在なので寄与ゼロ)
                rt.eff[r][idx] = kv.second.as<double>();
            }
        }
        if (equations) equations->push_back(eq);
        ++r;
    }
    rt.nReac = r;
}

} // namespace chem_io
