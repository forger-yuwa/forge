#include "input/speciesDB.hpp"
#include "input/solverConfig.hpp"

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

// 記録 1 種分 (互換性テキストの材料)
struct RecordEntry {
    std::string   name;
    std::string   phase;
    SpeciesThermo sp;
};

// 互換性テキスト本体。書式は tools/forge_species.py compat_text と一字一句同じにすること。
std::string compatTextRaw(const std::string& schema, const std::string& datum, const std::string& extrap,
                          double Tref, const std::vector<RecordEntry>& es)
{
    std::ostringstream o;
    o << "schema: " << schema << "\n";
    o << "datum: " << datum << "\n";
    o << "thermoHrefTemp: " << g17(Tref) << "\n";
    o << "extrapolation: " << extrap << "\n";
    o << "nSpecies: " << es.size() << "\n";
    for (size_t i = 0; i < es.size(); ++i) {
        const SpeciesThermo& s = es[i].sp;
        o << "species[" << i << "]: name=" << es[i].name << " phase=" << es[i].phase << "\n";
        o << "species[" << i << "].MW: " << g17(s.MW) << "\n";
        o << "species[" << i << "].T: " << g17(s.Tlo) << " " << g17(s.Tmid) << " " << g17(s.Thi) << "\n";
        o << "species[" << i << "].LJ: " << g17(s.sigma_LJ) << " " << g17(s.eps_kB) << "\n";
        o << "species[" << i << "].low:";
        for (int k = 0; k < 9; ++k) o << " " << g17(s.low[k]);
        o << "\n";
        o << "species[" << i << "].high:";
        for (int k = 0; k < 9; ++k) o << " " << g17(s.high[k]);
        o << "\n";
    }
    return o.str();
}

std::vector<RecordEntry> entriesOf(const ResolvedSpeciesDB& db)
{
    std::vector<RecordEntry> es;
    for (int s = 0; s < db.size(); ++s) es.push_back({db.names[s], "gas", db.species[s]});   // 液相は #10 まで含めない
    return es;
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
    return compatTextRaw(SPECIES_RECORD_SCHEMA, SPECIES_RECORD_DATUM, SPECIES_RECORD_EXTRAPOLATION, Tref, entriesOf(db));
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
    o << "schema: " << yq(SPECIES_RECORD_SCHEMA) << "\n";
    o << "compat_hash: " << yq(speciesDB_compatHash(db, Tref)) << "\n";
    o << "datum: " << yq(SPECIES_RECORD_DATUM) << "\n";
    o << "thermoHrefTemp: " << g17(Tref) << "\n";
    o << "extrapolation: " << yq(SPECIES_RECORD_EXTRAPOLATION) << "\n";
    o << "species:\n";
    for (int s = 0; s < db.size(); ++s) {
        const SpeciesThermo& sp = db.species[s];
        o << "  - index: " << s << "\n";
        o << "    name: " << yq(db.names[s]) << "\n";
        o << "    phase: \"gas\"\n";
        o << "    source: " << yq(db.source[s]) << "\n";
        o << "    MW: " << g17(sp.MW) << "\n";
        o << "    Tlo: " << g17(sp.Tlo) << "\n";
        o << "    Tmid: " << g17(sp.Tmid) << "\n";
        o << "    Thi: " << g17(sp.Thi) << "\n";
        o << "    LJ_sigma: " << g17(sp.sigma_LJ) << "\n";
        o << "    LJ_eps_kB: " << g17(sp.eps_kB) << "\n";
        o << "    nasa9_low: [";
        for (int k = 0; k < 9; ++k) o << (k ? ", " : "") << g17(sp.low[k]);
        o << "]\n";
        o << "    nasa9_high: [";
        for (int k = 0; k < 9; ++k) o << (k ? ", " : "") << g17(sp.high[k]);
        o << "]\n";
    }
    o << "provenance:\n";
    o << "  speciesDBFile: " << yq(dbFile) << "\n";
    o << "  input_field: " << yq(inputField) << "\n";
    o << "  input_status: " << yq(inputStatus) << "\n";
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
            SpeciesThermo& s = e.sp;
            s = SpeciesThermo{};
            s.MW = n["MW"].as<double>();
            s.Tlo = n["Tlo"].as<double>(); s.Tmid = n["Tmid"].as<double>(); s.Thi = n["Thi"].as<double>();
            s.sigma_LJ = n["LJ_sigma"].as<double>(); s.eps_kB = n["LJ_eps_kB"].as<double>();
            if (!read9(n["nasa9_low"], s.low) || !read9(n["nasa9_high"], s.high))
                throw std::runtime_error("species '" + e.name + "' lacks 9 nasa9 coefficients");
            es.push_back(e);
        }
    } catch (const std::exception& e) {
        d.push_back("record " + recordPath + " is malformed: " + e.what());
        return d;
    }
    // 記録の自己整合: 中身から互換ハッシュを作り直して記録内の値と比べる (改竄・取り違えの検出)
    const std::string rehash = speciesDB_sha256Hex(compatTextRaw(schema, datum, extrap, TrefRec, es));
    if (rehash != hashInFile) {
        d.push_back("record " + recordPath + ": compat_hash in file " + hashInFile.substr(0, 16)
                    + " != recomputed from its content " + rehash.substr(0, 16) + " (edited or corrupted record)");
    }
    if (schema != SPECIES_RECORD_SCHEMA) d.push_back("schema: field '" + schema + "' vs current '" SPECIES_RECORD_SCHEMA "'");
    if (datum != SPECIES_RECORD_DATUM) d.push_back("datum convention: field '" + datum + "' vs current");
    if (extrap != SPECIES_RECORD_EXTRAPOLATION) d.push_back("extrapolation: field '" + extrap + "' vs current");
    if (TrefRec != Tref) d.push_back("thermoHrefTemp: field " + g17(TrefRec) + " vs current " + g17(Tref));
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
    for (int s = 0; s < db.size(); ++s) {
        const SpeciesThermo& f = es[s].sp;
        const SpeciesThermo& c = db.species[s];
        const std::string& n = db.names[s];
        auto cmp = [&](const std::string& key, double a, double b) {
            if (a != b) d.push_back(n + "." + key + ": field " + g17(a) + " vs current " + g17(b));
        };
        if (es[s].phase != "gas") d.push_back(n + ".phase: field " + es[s].phase + " vs current gas");
        cmp("MW", f.MW, c.MW);
        cmp("Tlo", f.Tlo, c.Tlo); cmp("Tmid", f.Tmid, c.Tmid); cmp("Thi", f.Thi, c.Thi);
        cmp("LJ_sigma", f.sigma_LJ, c.sigma_LJ); cmp("LJ_eps_kB", f.eps_kB, c.eps_kB);
        for (int k = 0; k < 9; ++k) cmp("nasa9_low[" + std::to_string(k) + "]", f.low[k], c.low[k]);
        for (int k = 0; k < 9; ++k) cmp("nasa9_high[" + std::to_string(k) + "]", f.high[k], c.high[k]);
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
              "stamped by their generator).\n"
              "  - If you have checked that the field was produced with the same species/datum as this run, allow it for THIS "
              "invocation only:  FORGE_ALLOW_UNVERIFIED_SPECIES=1 forge   (not a config key; outputs are marked species_input_unverified=1)\n"
              "  - If the species set/DB changed, convert the field with tools/convert_species_field.py.";
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
