#include "input/speciesTransportDB.hpp"
#include "input/speciesDB.hpp"

// 共通 species データと輸送データ (CEA trans.inp 由来) の埋め込み (cmake/embed_species_data.cmake)。
#include "forge_species_data.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>
#include <stdexcept>

// =============================================================================
// speciesTransportDB.cpp — physProp.transport の解決 (plan thermophysics-solver-owned-species-db #5t2 段 1)
//   式は cuda_forge/transportMix_d.cuh、規約の文字列は input/speciesTransportDB.hpp。
// =============================================================================

namespace {

#define FORGE_TRANSPORT_DATA_SCHEMA "forge_transport_data_v1"

using Row = std::array<double, 6>;   // Tlo, Thi, A, B, C, D

struct TransEntry {
    std::string      id;          // trans.inp の名前 (単成分) / "a/b" (相互作用)
    std::string      reference;
    std::vector<Row> V, C;
};

struct TransData {
    std::string                        transSha;   // 元の trans.inp の SHA-256
    std::vector<TransEntry>            species;
    std::vector<std::pair<std::pair<std::string, std::string>, TransEntry>> pairs;
};

std::string upper(std::string s)
{
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::string g17(double v)
{
    char b[64];
    std::snprintf(b, sizeof(b), "%.17g", v);
    return b;
}

std::vector<Row> readRows(const YAML::Node& n, const std::string& what)
{
    std::vector<Row> out;
    if (!n) return out;
    if (!n.IsSequence()) throw std::runtime_error(what + ": must be a list of [Tlo, Thi, A, B, C, D]");
    for (const auto& r : n) {
        if (!r.IsSequence() || r.size() != 6) throw std::runtime_error(what + ": each interval must be [Tlo, Thi, A, B, C, D]");
        Row x;
        for (int i = 0; i < 6; ++i) {
            try {
                x[i] = r[i].as<double>();
            } catch (const std::exception&) {
                throw std::runtime_error(what + ": interval value is not a number");
            }
            if (!std::isfinite(x[i])) throw std::runtime_error(what + ": interval value is not finite");
        }
        out.push_back(x);
    }
    return out;
}

TransData parseTransData()
{
    const std::string where = std::string("[transport] built-in transport data ") + kForgeTransportDataName;
    YAML::Node root;
    try {
        root = YAML::Load(kForgeTransportDataYaml);
    } catch (const std::exception& e) {
        throw std::runtime_error(where + " is not valid YAML: " + e.what());
    }
    if (!root["schema"] || root["schema"].as<std::string>() != FORGE_TRANSPORT_DATA_SCHEMA) {
        throw std::runtime_error(where + ": schema must be '" FORGE_TRANSPORT_DATA_SCHEMA "'");
    }
    TransData d;
    d.transSha = root["provenance"]["cea_trans_inp"]["sha256"].as<std::string>();
    for (const auto& n : root["species"]) {
        TransEntry e;
        e.id = n["id"].as<std::string>();
        e.reference = n["reference"] ? n["reference"].as<std::string>() : "";
        e.V = readRows(n["V"], where + " " + e.id + ".V");
        e.C = readRows(n["C"], where + " " + e.id + ".C");
        d.species.push_back(e);
    }
    for (const auto& n : root["interactions"]) {
        TransEntry e;
        const std::string a = n["pair"][0].as<std::string>(), b = n["pair"][1].as<std::string>();
        e.id = a + "/" + b;
        e.reference = n["reference"] ? n["reference"].as<std::string>() : "";
        e.V = readRows(n["V"], where + " " + e.id + ".V");
        e.C = readRows(n["C"], where + " " + e.id + ".C");
        d.pairs.push_back({{a, b}, e});
    }
    return d;
}

const TransData& transData()
{
    static const TransData d = parseTransData();
    return d;
}

// 内蔵種の双極子モーメント [D] (共通データの LJ.dipole; 無ければ 0)
const std::map<std::string, double>& builtinDipoles()
{
    static const std::map<std::string, double> m = [] {
        std::map<std::string, double> out;
        const YAML::Node root = YAML::Load(kForgeSpeciesDataYaml);
        for (const auto& n : root["species"]) {
            const YAML::Node lj = n["LJ"];
            if (lj && lj.IsMap() && lj["dipole"]) out[n["id"].as<std::string>()] = lj["dipole"].as<double>();
        }
        return out;
    }();
    return m;
}

// trans.inp の単成分エントリ (完全一致 → 大小文字無視; 生成器が大小文字だけ違う種名が無いことを保証)
const TransEntry* findTransSpecies(const std::string& name)
{
    for (const auto& e : transData().species) if (e.id == name) return &e;
    const std::string up = upper(name);
    for (const auto& e : transData().species) if (upper(e.id) == up) return &e;
    return nullptr;
}

const TransEntry* findTransPair(const std::string& a, const std::string& b)
{
    for (const auto& p : transData().pairs) {
        if ((p.first.first == a && p.first.second == b) || (p.first.first == b && p.first.second == a)) return &p.second;
    }
    return nullptr;
}

TransportFitD toFit(const std::vector<Row>& rows, const std::string& what)
{
    if (rows.empty()) throw std::runtime_error(what + ": no intervals");
    if (static_cast<int>(rows.size()) > TRANSPORT_MAX_FIT_INTERVALS) {
        throw std::runtime_error(what + ": " + std::to_string(rows.size()) + " intervals exceed TRANSPORT_MAX_FIT_INTERVALS="
                                 + std::to_string(TRANSPORT_MAX_FIT_INTERVALS));
    }
    TransportFitD f{};
    f.n = static_cast<int>(rows.size());
    for (int k = 0; k < f.n; ++k) {
        const Row& r = rows[k];
        if (!(r[0] > 0.0) || !(r[0] < r[1])) throw std::runtime_error(what + ": interval " + std::to_string(k) + " needs 0 < Tlo < Thi");
        if (k > 0 && !(rows[k - 1][1] <= r[0])) throw std::runtime_error(what + ": intervals must be increasing and not overlap");
        f.Tlo[k] = r[0]; f.Thi[k] = r[1]; f.A[k] = r[2]; f.B[k] = r[3]; f.C[k] = r[4]; f.D[k] = r[5];
    }
    return f;
}

std::string fitLine(const TransportFitD& f)
{
    std::ostringstream o;
    o << "n=" << f.n;
    for (int k = 0; k < f.n; ++k) {
        o << " [" << g17(f.Tlo[k]) << " " << g17(f.Thi[k]) << " " << g17(f.A[k]) << " " << g17(f.B[k]) << " "
          << g17(f.C[k]) << " " << g17(f.D[k]) << "]";
    }
    return o.str();
}

// 実種 (lump 展開後) の候補
struct RealSp {
    std::string              identity;   // speciesDB_identityKey
    std::vector<std::string> spellings;  // config に現れた綴り (種名・構成種名)
    std::string              dbKey;      // DB のキー
    bool                     fromFile = false;
    SpeciesThermo            thermo{};
    std::string              thermoFrom; // "species[s]" / "species[s].lump[k]"
};

const char* kModels = "cea, kinetic, fit, custom:h2o_iapws_cea_v1";

} // anonymous namespace

const std::string& speciesTransportDB_dataName()
{
    static const std::string s = kForgeTransportDataName;
    return s;
}

const std::string& speciesTransportDB_dataSha256()
{
    static const std::string s = kForgeTransportDataSha256;
    return s;
}

const std::string& speciesTransportDB_transInpSha256()
{
    return transData().transSha;
}

void speciesTransportDB_resolve(ResolvedSpeciesDB& db, const std::vector<std::pair<std::string, std::string>>& spec,
                                const std::string& dbFile)
{
    db.transport = ResolvedTransport{};
    if (spec.empty()) return;
    const std::string W = "[transport] physProp.transport";

    // ---- 1. 実種の列 (lump は構成種へ展開、同じ実種は 1 つにまとめて分率を合算) ----
    std::vector<RealSp> reals;
    ResolvedTransport tr;
    tr.enabled = true;
    tr.expand.resize(db.size());
    auto addReal = [&](int s, const std::string& spelling, const std::string& key, bool file, const SpeciesThermo& th,
                       const std::string& from, double x) {
        const std::string id = speciesDB_identityKey(key, file);
        int r = -1;
        for (int q = 0; q < static_cast<int>(reals.size()); ++q) {
            if (reals[q].identity == id && reals[q].fromFile == file) { r = q; break; }
        }
        if (r < 0) {
            RealSp rs;
            rs.identity = id; rs.dbKey = key; rs.fromFile = file; rs.thermo = th; rs.thermoFrom = from;
            reals.push_back(rs);
            r = static_cast<int>(reals.size()) - 1;
        }
        if (std::find(reals[r].spellings.begin(), reals[r].spellings.end(), spelling) == reals[r].spellings.end())
            reals[r].spellings.push_back(spelling);
        tr.expand[s].push_back({r, x});
    };
    for (int s = 0; s < db.size(); ++s) {
        if (db.isLump(s)) {
            const ResolvedLump& l = db.lumps[s];
            for (size_t k = 0; k < l.members.size(); ++k) {
                addReal(s, l.members[k], l.memberDbKey[k], l.memberSource[k] == "file", l.memberSpecies[k],
                        "species[" + std::to_string(s) + "].lump[" + std::to_string(k) + "]", l.x[k]);
            }
        } else {
            addReal(s, db.names[s], db.dbKey[s], db.source[s] == "file", db.species[s], "species[" + std::to_string(s) + "]", 1.0);
        }
    }
    const int n = static_cast<int>(reals.size());
    if (n > TRANSPORT_MAX_REAL_SPECIES) {
        throw std::runtime_error(W + ": " + std::to_string(n) + " real species after lump expansion exceed TRANSPORT_MAX_REAL_SPECIES="
                                 + std::to_string(TRANSPORT_MAX_REAL_SPECIES));
    }
    auto realList = [&]() {
        std::string o;
        for (const auto& r : reals) o += (o.empty() ? "" : ", ") + r.spellings.front();
        return o;
    };

    // ---- 2. 指定キー → 実種 ----
    std::vector<int> specOf(n, -1);
    for (int k = 0; k < static_cast<int>(spec.size()); ++k) {
        const std::string& key = spec[k].first;
        const std::string idk = speciesDB_identityKey(key, false);
        std::vector<int> hit;
        for (int r = 0; r < n; ++r) {
            bool m = (reals[r].identity == idk) || (upper(reals[r].identity) == upper(key));
            for (const auto& sp : reals[r].spellings) m = m || (upper(sp) == upper(key));
            if (m) hit.push_back(r);
        }
        if (hit.empty()) {
            for (int s = 0; s < db.size(); ++s) {
                if (db.isLump(s) && upper(db.names[s]) == upper(key)) {
                    std::string mem;
                    for (const auto& m : db.lumps[s].members) mem += (mem.empty() ? "" : ", ") + m;
                    throw std::runtime_error(W + ": '" + key + "' is a lump; give the model of each constituent instead (" + mem + ")");
                }
            }
            throw std::runtime_error(W + ": '" + key + "' is not a species of this run (real species after lump expansion: " + realList() + ")");
        }
        if (hit.size() > 1) {
            throw std::runtime_error(W + ": '" + key + "' matches more than one real species (built-in and speciesDBFile entries of the same name); spell it exactly");
        }
        const int r = hit.front();
        if (specOf[r] >= 0) {
            throw std::runtime_error(W + ": species '" + reals[r].spellings.front() + "' is given twice ('" + spec[specOf[r]].first
                                     + "' and '" + key + "')");
        }
        specOf[r] = k;
    }
    {
        std::string miss;
        for (int r = 0; r < n; ++r) if (specOf[r] < 0) miss += (miss.empty() ? "" : ", ") + reals[r].spellings.front();
        if (!miss.empty()) {
            throw std::runtime_error(W + ": no transport model for species " + miss
                                     + ". Every species (lump constituents individually) needs one of: " + kModels);
        }
    }

    // ---- 3. 実種ごとのモデル ----
    YAML::Node dbRoot;
    if (!dbFile.empty()) {
        try {
            dbRoot = YAML::LoadFile(dbFile);
        } catch (const std::exception& e) {
            throw std::runtime_error(W + ": failed to read speciesDBFile '" + dbFile + "': " + e.what());
        }
    }
    auto dbEntry = [&](const RealSp& r) -> YAML::Node {
        if (!r.fromFile || !dbRoot) return YAML::Node();
        return dbRoot[r.dbKey];
    };
    std::vector<std::string> transName(n);   // trans.inp の名前 (無ければ空; 相互作用の検索に使う)
    bool anyKinetic = false, anyH2O = false;
    for (int r = 0; r < n; ++r) {
        const RealSp& rs = reals[r];
        const std::string& model = spec[specOf[r]].second;
        const std::string nm = rs.spellings.front();
        const std::string where = W + " '" + spec[specOf[r]].first + "': " + model;
        SpeciesTransportD d{};
        d.MW = rs.thermo.MW;
        d.thermo = rs.thermo;
        const std::string tname = rs.fromFile ? rs.dbKey : rs.identity;
        const TransEntry* te = findTransSpecies(tname);
        transName[r] = te ? te->id : "";
        std::string src;
        double dip = 0.0;
        if (model == "cea" || model == "custom:h2o_iapws_cea_v1") {
            if (model == "custom:h2o_iapws_cea_v1") {
                if (rs.identity != "H2O") {
                    throw std::runtime_error(where + ": this custom model is defined only for H2O (species '" + nm + "')");
                }
                d.model = TRANSPORT_MODEL_H2O_IAPWS_CEA_V1;
                anyH2O = true;
            } else {
                d.model = TRANSPORT_MODEL_CEA;
            }
            if (!te) {
                throw std::runtime_error(where + ": species '" + nm + "' has no entry in CEA trans.inp (" + speciesTransportDB_dataName()
                                         + "); choose kinetic or fit for it");
            }
            if (te->V.empty() || te->C.empty()) {
                throw std::runtime_error(where + ": CEA trans.inp entry '" + te->id + "' lacks viscosity or conductivity data; choose kinetic or fit");
            }
            d.V = toFit(te->V, where + " trans.inp " + te->id + " V");
            d.C = toFit(te->C, where + " trans.inp " + te->id + " C");
            src = "trans.inp " + te->id + " (" + te->reference + ")";
            if (d.model == TRANSPORT_MODEL_H2O_IAPWS_CEA_V1) {
                src = "IAPWS 2008/2011 dilute gas below 500 K + " + src + " above 700 K";
                tr.notes.push_back(nm + " (custom:h2o_iapws_cea_v1): IAPWS dilute-gas formulas are used as is down to 150 K; "
                                   "below 253.15 K this is outside their formal range (IAPWS 2008/2011: 253.15-1173.15 K); "
                                   "below 150 K a power law matched to the 150 K log-slope (C1)");
            } else {
                const double lo = std::min(d.V.Tlo[0], d.C.Tlo[0]), hi = std::max(d.V.Thi[d.V.n - 1], d.C.Thi[d.C.n - 1]);
                char rng[96];
                std::snprintf(rng, sizeof(rng), "%.6g-%.6g K", lo, hi);
                tr.notes.push_back(nm + " (cea): data " + rng + "; outside, the nearest interval is extrapolated (as CEA, cea2.f TRANIN)");
            }
        } else if (model == "kinetic") {
            d.model = TRANSPORT_MODEL_KINETIC;
            anyKinetic = true;
            if (rs.fromFile) {
                const YAML::Node e = dbEntry(rs);
                if (!e || !e["LJ_sigma"] || !e["LJ_eps_kB"]) {
                    throw std::runtime_error(where + ": species '" + nm + "' from speciesDBFile has no explicit LJ_sigma/LJ_eps_kB "
                                             "(the default LJ placeholder is not used for transport)");
                }
                if (e["LJ_dipole"]) dip = e["LJ_dipole"].as<double>();
            } else {
                auto it = builtinDipoles().find(rs.identity);
                if (it != builtinDipoles().end()) dip = it->second;
            }
            d.sigma_LJ = rs.thermo.sigma_LJ;
            d.eps_kB = rs.thermo.eps_kB;
            if (!(d.sigma_LJ > 0.0) || !(d.eps_kB > 0.0) || !std::isfinite(d.sigma_LJ) || !std::isfinite(d.eps_kB) || !(dip >= 0.0)) {
                throw std::runtime_error(where + ": invalid LJ parameters / dipole for '" + nm + "'");
            }
            // δ* = μ_d² / (2 ε σ³) (CGS: μ_d [D] = 1e-18 statC cm, ε = ε/k_B · k_B [erg], σ [Å] = 1e-8 cm)
            const double mu_d = dip*1.0e-18, sig = d.sigma_LJ*1.0e-8;
            d.deltaStar = mu_d*mu_d/(2.0*d.eps_kB*1.380649e-16*sig*sig*sig);
            src = std::string(rs.fromFile ? "speciesDBFile" : "built-in") + " LJ sigma=" + g17(d.sigma_LJ) + " eps/kB=" + g17(d.eps_kB)
                  + (dip > 0.0 ? " dipole=" + g17(dip) + " D" : "");
        } else if (model == "fit") {
            d.model = TRANSPORT_MODEL_FIT;
            const YAML::Node e = dbEntry(rs);
            if (!rs.fromFile || !e || !e["transport_fit"] || !e["transport_fit"].IsMap()) {
                throw std::runtime_error(where + ": species '" + nm + "' needs 'transport_fit: {V: [[Tlo,Thi,A,B,C,D],...], C: [...]}' in its "
                                         "speciesDBFile entry (CEA form: ln f = A lnT + B/T + C/T^2 + D; V in microPoise, C in microW/(cm K))");
            }
            const YAML::Node f = e["transport_fit"];
            for (auto it = f.begin(); it != f.end(); ++it) {
                const std::string k = it->first.as<std::string>();
                if (k != "V" && k != "C" && k != "reference") {
                    throw std::runtime_error(where + ": unknown key '" + k + "' in transport_fit (allowed: V, C, reference)");
                }
            }
            d.V = toFit(readRows(f["V"], where + " transport_fit.V"), where + " transport_fit.V");
            d.C = toFit(readRows(f["C"], where + " transport_fit.C"), where + " transport_fit.C");
            src = "speciesDBFile transport_fit" + (f["reference"] ? " (" + f["reference"].as<std::string>() + ")" : std::string(""));
        } else if (model.rfind("custom:", 0) == 0) {
            throw std::runtime_error(where + ": unknown custom model (known: custom:h2o_iapws_cea_v1)");
        } else {
            throw std::runtime_error(where + ": unknown transport model (allowed: " + kModels + ")");
        }
        tr.realName.push_back(rs.identity);
        tr.realConfigKey.push_back(spec[specOf[r]].first);
        tr.modelName.push_back(model);
        tr.dataSource.push_back(src);
        tr.thermoFrom.push_back(rs.thermoFrom);
        tr.dipoleDebye.push_back(dip);
        tr.sp.push_back(d);
    }

    // ---- 4. 相互作用粘性 η_ij の出所 (上三角 a<b) ----
    tr.pairs.assign(n*(n - 1)/2, TransportPairD{});
    tr.pairSource.assign(n*(n - 1)/2, "");
    bool anyRigid = false;
    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            const int k = transport_pair_index(a, b, n);
            TransportPairD p{};
            p.a = a; p.b = b;
            const SpeciesTransportD& A = tr.sp[a];
            const SpeciesTransportD& B = tr.sp[b];
            if (A.model == TRANSPORT_MODEL_KINETIC && B.model == TRANSPORT_MODEL_KINETIC) {
                p.kind = TRANSPORT_PAIR_CE;
                p.sigma_ab = 0.5*(A.sigma_LJ + B.sigma_LJ);
                p.eps_ab = std::sqrt(A.eps_kB*B.eps_kB);
                tr.pairSource[k] = "chapman_enskog";
            } else {
                const TransEntry* pe = (!transName[a].empty() && !transName[b].empty()) ? findTransPair(transName[a], transName[b]) : nullptr;
                if (pe && !pe->V.empty()) {
                    p.kind = TRANSPORT_PAIR_CEA;
                    p.V = toFit(pe->V, W + " trans.inp " + pe->id + " V");
                    tr.pairSource[k] = "trans.inp " + pe->id + " (" + pe->reference + ")";
                } else {
                    p.kind = TRANSPORT_PAIR_RIGID_SPHERE;
                    tr.pairSource[k] = "rigid_sphere";
                    anyRigid = true;
                }
            }
            tr.pairs[k] = p;
        }
    }
    if (anyRigid) tr.notes.push_back("pairs without CEA interaction data use the CEA rigid-sphere estimate (accuracy of mixed-model pairs not verified)");
    if (anyKinetic) tr.notes.push_back("kinetic: lambda is the modified Eucken form (Warnatz form undecided, plan section 10)");

    // ---- 5. 互換性ハッシュ・記録の正規化行 ----
    auto& L = tr.compatLines;
    L.push_back(std::string("transport: schema=") + TRANSPORT_RECORD_SCHEMA);
    L.push_back(std::string("transport.mixing: ") + TRANSPORT_MIXING_RULE);
    L.push_back(std::string("transport.eta_ij: ") + TRANSPORT_ETA_IJ_RULE);
    L.push_back(std::string("transport.fit: ") + TRANSPORT_FIT_RULE);
    if (anyKinetic) L.push_back(std::string("transport.kinetic: ") + TRANSPORT_KINETIC_RULE);
    if (anyH2O)     L.push_back(std::string("transport.h2o_iapws_cea_v1: ") + TRANSPORT_H2O_V1_RULE);
    L.push_back("transport.nReal: " + std::to_string(n));
    for (int r = 0; r < n; ++r) {
        const std::string t = "transport.real[" + std::to_string(r) + "]";
        const SpeciesTransportD& d = tr.sp[r];
        L.push_back(t + ": name=" + tr.realName[r] + " model=" + tr.modelName[r]);
        L.push_back(t + ".MW: " + g17(d.MW));
        if (d.model == TRANSPORT_MODEL_KINETIC) {
            L.push_back(t + ".kinetic: sigma=" + g17(d.sigma_LJ) + " eps_kB=" + g17(d.eps_kB) + " dipole_debye=" + g17(tr.dipoleDebye[r])
                        + " delta_star=" + g17(d.deltaStar) + " cp_from=" + tr.thermoFrom[r]);
        } else {
            L.push_back(t + ".V: " + fitLine(d.V));
            L.push_back(t + ".C: " + fitLine(d.C));
        }
    }
    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            const TransportPairD& p = tr.pairs[transport_pair_index(a, b, n)];
            std::string l = "transport.pair[" + std::to_string(a) + "," + std::to_string(b) + "]: ";
            if (p.kind == TRANSPORT_PAIR_CE)        l += "kind=chapman_enskog sigma=" + g17(p.sigma_ab) + " eps_kB=" + g17(p.eps_ab);
            else if (p.kind == TRANSPORT_PAIR_CEA)  l += "kind=cea_interaction V: " + fitLine(p.V);
            else                                    l += "kind=rigid_sphere";
            L.push_back(l);
        }
    }
    for (int s = 0; s < db.size(); ++s) {
        std::string l = "transport.expand[" + std::to_string(s) + "]: name=" + db.names[s] + " ->";
        for (const auto& e : tr.expand[s]) l += " " + std::to_string(e.first) + ":" + g17(e.second);
        L.push_back(l);
    }
    db.transport = tr;
}

void speciesTransportDB_mixture(const ResolvedTransport& tr, const std::vector<double>& Xs, double T,
                                double& mu, double& lam, std::vector<double>* Xreal,
                                std::vector<double>* mu_i, std::vector<double>* lam_i)
{
    if (!tr.enabled) throw std::runtime_error("[transport] speciesTransportDB_mixture: transport is not resolved (physProp.transport absent)");
    if (Xs.size() != tr.expand.size()) throw std::runtime_error("[transport] speciesTransportDB_mixture: mole fraction count mismatch");
    const int n = tr.nReal();
    std::vector<double> X(n, 0.0);
    for (size_t s = 0; s < Xs.size(); ++s)
        for (const auto& e : tr.expand[s]) X[e.first] += Xs[s]*e.second;
    std::vector<double> m(n), l(n);
    transport_mix(tr.sp.data(), tr.pairs.data(), n, X.data(), T, &mu, &lam, m.data(), l.data());
    if (Xreal) *Xreal = X;
    if (mu_i) *mu_i = m;
    if (lam_i) *lam_i = l;
}
