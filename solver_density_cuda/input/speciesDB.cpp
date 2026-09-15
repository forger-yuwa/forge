#include "input/speciesDB.hpp"
#include "input/solverConfig.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {

// 内蔵 DB の 1 種を構築する。
// NASA-9 係数 (a0..a8): cp/R = a0/T^2 + a1/T + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
// MW [kg/mol], sigma_LJ [Angstrom], eps_kB [K]。
// 各 species・各定数の出典 (CEA / Svehla / forge 擬似種 AIR の区別) は
// methods/thermophysics.md 「内蔵 species DB の一覧と出典」を参照。
SpeciesThermo makeSpecies(double MW, double sigma, double eps_kB,
                          double Tlo, double Tmid, double Thi,
                          const double low[9], const double high[9])
{
    SpeciesThermo s;
    s.MW = MW; s.sigma_LJ = sigma; s.eps_kB = eps_kB;
    s.Tlo = Tlo; s.Tmid = Tmid; s.Thi = Thi;
    for (int i=0;i<9;i++){ s.low[i]=low[i]; s.high[i]=high[i]; }
    s.h_datum = 0.0;
    s.invMW = 1.0/MW;
    return s;
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
    std::map<std::string, SpeciesThermo> db;

    // ---- N2 (CEA) ----
    {
        const double low[9] = { 2.210371497e+04,-3.818461820e+02, 6.082738360e+00,
                               -8.530914410e-03, 1.384646189e-05,-9.625793620e-09,
                                2.519705809e-12, 7.108460860e+02,-1.076003744e+01 };
        const double high[9]= { 5.877124060e+05,-2.239249073e+03, 6.066949220e+00,
                               -6.139685500e-04, 1.491806679e-07,-1.923105485e-11,
                                1.061954386e-15, 1.283210415e+04,-1.586640027e+01 };
        db["N2"] = makeSpecies(0.0280134, 3.621, 97.53, 200.0, 1000.0, 6000.0, low, high);
    }
    // ---- O2 (CEA) ----
    {
        const double low[9] = {-3.425563420e+04, 4.847000970e+02, 1.119010961e+00,
                                4.293889240e-03,-6.836300520e-07,-2.023372700e-09,
                                1.039040018e-12,-3.391454870e+03, 1.849699470e+01 };
        const double high[9]= {-1.037939022e+06, 2.344830282e+03, 1.819732036e+00,
                                1.267847582e-03,-2.188067988e-07, 2.053719572e-11,
                               -8.193467050e-16,-1.689010929e+04, 1.738716506e+01 };
        db["O2"] = makeSpecies(0.0319988, 3.458, 107.4, 200.0, 1000.0, 6000.0, low, high);
    }
    // ---- Ar (単原子, cp/R = 2.5) ----
    {
        const double low[9] = { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 4.379674910e+00 };
        const double high[9]= { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 4.379674910e+00 };
        db["AR"] = makeSpecies(0.039948, 3.330, 136.5, 200.0, 1000.0, 6000.0, low, high);
        db["Ar"] = db["AR"];
    }
    // ---- CO2 (CEA) ----
    {
        const double low[9] = { 4.943650540e+04,-6.264116010e+02, 5.301725240e+00,
                                2.503813816e-03,-2.127308728e-07,-7.689988780e-10,
                                2.849677801e-13,-4.528198460e+04,-7.048279440e+00 };
        const double high[9]= { 1.176962419e+05,-1.788791477e+03, 8.291523190e+00,
                               -9.223156780e-05, 4.863676880e-09,-1.891053312e-12,
                                6.330036590e-16,-3.908350590e+04,-2.652669281e+01 };
        db["CO2"] = makeSpecies(0.0440095, 3.763, 244.0, 200.0, 1000.0, 6000.0, low, high);
    }
    // ---- He (単原子) ----
    {
        const double low[9] = { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 9.287239740e-01 };
        const double high[9]= { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 9.287239740e-01 };
        db["HE"] = makeSpecies(0.0040026, 2.551, 10.22, 200.0, 1000.0, 6000.0, low, high);
        db["He"] = db["HE"];
    }
    // ---- H2O (CEA, McBride-Gordon 2002)。非平衡凝縮 (Wyslouzil) のキャリア+凝縮種で使用 ----
    {
        const double low[9] = {-3.947960830e+04, 5.755731020e+02, 9.317826530e-01,
                                7.222712860e-03,-7.342557370e-06, 4.955043490e-09,
                               -1.336933246e-12,-3.303974310e+04, 1.724205775e+01 };
        const double high[9]= { 1.034972096e+06,-2.412698562e+03, 4.646110780e+00,
                                2.291998307e-03,-6.836830480e-07, 9.426468930e-11,
                               -4.822380530e-15,-1.384286509e+04,-7.978148510e+00 };
        db["H2O"] = makeSpecies(0.0180153, 2.605, 572.4, 200.0, 1000.0, 6000.0, low, high);
        db["h2o"] = db["H2O"];
        db["WATER"] = db["H2O"];
    }
    // ---- AIR (単成分擬似空気: 二原子 cp/R≈3.5 + 空気の MW) ----
    // CEA の「Air」混合に近い熱力学を単成分で近似したい場合の簡便種。
    // 低温で cp≈1004.5 J/kgK となるよう定数 cp/R=3.5 (γ=1.4) の擬似多項式。
    {
        const double c[9] = { 0.0,0.0,3.5,0.0,0.0,0.0,0.0,-1.0431373e+03, 3.0 };
        db["AIR"] = makeSpecies(0.0289647, 3.711, 78.6, 200.0, 1000.0, 6000.0, c, c);
        db["Air"] = db["AIR"];
        db["air"] = db["AIR"];
    }

    return db;
}

ResolvedSpeciesDB speciesDB_resolve(const std::vector<std::string>& namesIn, const std::string& dbFile)
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
            sp.Tlo      = s["Tlo"]      ? s["Tlo"].as<double>()      : 200.0;
            sp.Tmid     = s["Tmid"]     ? s["Tmid"].as<double>()     : 1000.0;
            sp.Thi      = s["Thi"]      ? s["Thi"].as<double>()      : 6000.0;
            double lo[9], hi[9];
            if (read9(s["nasa9_low"], lo) && read9(s["nasa9_high"], hi)) {
                for (int i=0;i<9;i++){ sp.low[i]=lo[i]; sp.high[i]=hi[i]; }
                if (!(sp.MW > 0.0) || !std::isfinite(sp.MW)) {
                    throw std::runtime_error("[speciesDB] species '" + name + "' has invalid MW in DB file '" + dbFile + "'");
                }
                sp.h_datum = 0.0;
                sp.invMW = 1.0/sp.MW;
                db[name] = sp;
                fromFile[name] = true;
            } else {
                std::cerr << "[speciesDB] species '" << name
                          << "' in DB file lacks valid nasa9_low/nasa9_high (need 9 coeffs); ignored" << std::endl;
            }
        }
    }

    std::vector<std::string> names = namesIn;
    if (names.empty()) names.push_back("N2");  // ダミー (calorically-perfect 用)

    if (static_cast<int>(names.size()) > THERMO_MAX_SPECIES) {
        throw std::runtime_error("[speciesDB] nSpecies=" + std::to_string(names.size())
                                 + " exceeds THERMO_MAX_SPECIES=" + std::to_string(THERMO_MAX_SPECIES));
    }

    ResolvedSpeciesDB out;
    for (const auto& nm : names) {
        if (out.index(nm) >= 0) {
            throw std::runtime_error("[speciesDB] species '" + nm + "' is listed twice in physProp.species");
        }
        auto it = findSpecies(db, nm);
        if (it == db.end()) {
            std::string avail;
            for (const auto& kv : db) avail += (avail.empty() ? "" : " ") + kv.first;
            throw std::runtime_error("[speciesDB] species '" + nm + "' not found in built-in DB nor speciesDBFile"
                                     + (dbFile.empty() ? std::string("") : " '" + dbFile + "'")
                                     + " (available: " + avail + ")");
        }
        out.names.push_back(nm);
        out.species.push_back(it->second);
        out.source.push_back(fromFile.count(it->first) ? "file" : "builtin");
    }
    for (auto& s : out.species) s.invMW = 1.0/s.MW;
    return out;
}

ResolvedSpeciesDB speciesDB_resolve(const solverConfig& cfg)
{
    return speciesDB_resolve(cfg.speciesNames, cfg.speciesDBFile);
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
    for (int s = 0; s < db.size(); ++s) {
        std::cout << "[species]   " << std::setw(2) << s << "  " << std::setw(8) << std::left << db.names[s] << std::right
                  << "  MW=" << std::setprecision(10) << db.MW(s) << "  " << db.source[s] << "\n";
    }
    if (cfg.condensation == 1) {
        if (cfg.condGasSpecies >= 0 && cfg.condGasSpecies < db.size()) {
            std::cout << "[species]   condensing species: " << db.names[cfg.condGasSpecies]
                      << " (condGasSpecies=" << cfg.condGasSpecies << ")\n";
        } else {
            std::cout << "[species]   condensing species: pure condensible (condGasSpecies=-1)\n";
        }
    }
    std::cout << "[species]   tracer: " << (cfg.tracerEnabled() ? cfg.tracer + " (roXi transported; inlet floats Xi)" : "none") << "\n";
}
