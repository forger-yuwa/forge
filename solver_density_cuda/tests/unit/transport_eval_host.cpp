// =============================================================================
// transport_eval_host.cpp — physProp.transport の host 評価ドライバ (GPU 不要; plan thermophysics-solver-owned-species-db #5t2 段 1)
//   input/speciesDB.cpp + input/speciesTransportDB.cpp の解決結果と cuda_forge/transportMix_d.cuh の式で、
//   与えた T・組成の μ・λ・実種ごとの単成分値・η_ij を JSON で出す。試験 test_species_transport.py が独立参照と比べる。
//
//   transport_eval_host SPEC.yaml
//     SPEC.yaml:
//       species: [N2, H2O]                         # physProp.species と同じ書式 (lump は {name, lump, basis})
//       transport: {N2: cea, H2O: cea}             # physProp.transport と同じ
//       speciesDBFile: path                        # 任意
//       thermoHrefTemp: 0                          # 任意 (互換性ハッシュ用)
//       ljSource: [gri30, svehla1962]              # 任意 (physProp.ljSource と同じ; 無指定は既定; plan #14)
//       states: [{T: 400, X: [0.5, 0.5]}, ...]     # X は physProp.species の順のモル分率
//   出力 (stdout, 1 行の JSON): {"ok": true, "real": [...], "model": [...], "pair_kind": [...], "compat_hash": ..., "compat_lines": [...],
//                                "record": "<記録全文>", "notes": [...], "data_source": [...], "dipole_debye": [...],
//                                "states": [{"T", "mu", "lam", "Xreal", "mu_i", "lam_i", "eta_ij", "D_ij" (輸送種の上三角, 101325 Pa)}]}
//   失敗: {"ok": false, "error": "..."} を出して終了コード 2。
//
// ビルド (共通データの埋め込みヘッダを先に生成; test_species_transport.py が行う):
//   cmake -DIN=solver_density_cuda/data/species/forge_species_v1.yaml -DOUT=<gen>/forge_species_data.hpp
//       -P solver_density_cuda/cmake/embed_species_data.cmake
//   g++ -O1 -std=c++17 -I solver_density_cuda -I <gen> solver_density_cuda/tests/unit/transport_eval_host.cpp
//       solver_density_cuda/input/speciesDB.cpp solver_density_cuda/input/speciesTransportDB.cpp -lyaml-cpp -o transport_eval_host
// =============================================================================
#include "input/speciesDB.hpp"

#include <yaml-cpp/yaml.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

static std::string js(const std::string& s)
{
    std::string o = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') { o.push_back('\\'); o.push_back(c); }
        else if (c == '\n') o += "\\n";
        else if (static_cast<unsigned char>(c) < 0x20) { char b[8]; std::snprintf(b, sizeof(b), "\\u%04x", c); o += b; }
        else o.push_back(c);
    }
    return o + "\"";
}

static std::string jnum(double v)
{
    char b[64];
    std::snprintf(b, sizeof(b), "%.17g", v);
    return b;
}

template <class T, class F>
static std::string jarr(const std::vector<T>& v, F f)
{
    std::string o = "[";
    for (size_t i = 0; i < v.size(); ++i) o += (i ? ", " : "") + f(v[i]);
    return o + "]";
}

int main(int argc, char** argv)
{
    if (argc != 2) { std::fprintf(stderr, "usage: transport_eval_host SPEC.yaml\n"); return 1; }
    try {
        const YAML::Node spec = YAML::LoadFile(argv[1]);
        std::vector<std::string> names;
        std::vector<SpeciesLumpSpec> lumps;
        for (const auto& sn : spec["species"]) {
            if (sn.IsScalar()) { names.push_back(sn.as<std::string>()); continue; }
            SpeciesLumpSpec lp;
            lp.name = sn["name"].as<std::string>();
            lp.basis = sn["basis"].as<std::string>();
            for (auto it = sn["lump"].begin(); it != sn["lump"].end(); ++it) {
                lp.members.push_back(it->first.as<std::string>());
                lp.fractions.push_back(it->second.as<double>());
            }
            names.push_back(lp.name);
            lumps.push_back(lp);
        }
        std::vector<std::pair<std::string, std::string>> tspec;
        if (spec["transport"]) {
            for (auto it = spec["transport"].begin(); it != spec["transport"].end(); ++it)
                tspec.emplace_back(it->first.as<std::string>(), it->second.as<std::string>());
        }
        const std::string dbFile = spec["speciesDBFile"] ? spec["speciesDBFile"].as<std::string>() : "";
        const double Tref = spec["thermoHrefTemp"] ? spec["thermoHrefTemp"].as<double>() : 0.0;

        std::vector<std::string> ljSource;
        if (spec["ljSource"]) for (const auto& s : spec["ljSource"]) ljSource.push_back(s.as<std::string>());
        ResolvedSpeciesDB db = speciesDB_resolve(names, dbFile, lumps, ljSource);
        speciesTransportDB_resolve(db, tspec, dbFile);
        const ResolvedTransport& tr = db.transport;

        std::string o = "{\"ok\": true";
        o += ", \"compat_hash\": " + js(speciesDB_compatHash(db, Tref));
        o += ", \"record\": " + js(speciesDB_recordText(db, Tref, dbFile, "", "not_checked_resolve_only"));
        o += ", \"transport_enabled\": " + std::string(tr.enabled ? "true" : "false");
        if (tr.enabled) {
            o += ", \"real\": " + jarr(tr.realName, js);
            o += ", \"model\": " + jarr(tr.modelName, js);
            o += ", \"pair_source\": " + jarr(tr.pairSource, js);
            std::vector<int> kinds;
            for (const auto& p : tr.pairs) kinds.push_back(p.kind);
            o += ", \"pair_kind\": " + jarr(kinds, [](int k) { return std::to_string(k); });
            o += ", \"MW_real\": " + jarr(tr.sp, [](const SpeciesTransportD& s) { return jnum(s.MW); });
            o += ", \"delta_star\": " + jarr(tr.sp, [](const SpeciesTransportD& s) { return jnum(s.deltaStar); });
            o += ", \"compat_lines\": " + jarr(tr.compatLines, js);
            o += ", \"notes\": " + jarr(tr.notes, js);
            o += ", \"data_source\": " + jarr(tr.dataSource, js);
            o += ", \"dipole_debye\": " + jarr(tr.dipoleDebye, jnum);
            o += ", \"states\": [";
            bool first = true;
            for (const auto& st : spec["states"]) {
                const double T = st["T"].as<double>();
                std::vector<double> X;
                for (const auto& x : st["X"]) X.push_back(x.as<double>());
                double mu, lam;
                std::vector<double> Xr, mi, li;
                speciesTransportDB_mixture(tr, X, T, mu, lam, &Xr, &mi, &li);
                const int n = tr.nReal();
                std::vector<double> eij;
                for (int a = 0; a < n; ++a)
                    for (int b = a + 1; b < n; ++b)
                        eij.push_back(transport_eta_pair(tr.pairs[transport_pair_index(a, b, n)], tr.sp.data(), mi.data(), T));
                // 輸送種どうしの二元拡散係数 (M4 の thermo_Dbinary; 101325 Pa)。双極子を読まないことの確認用 (#14-L1b)
                std::vector<double> dij;
                for (int a = 0; a < db.size(); ++a)
                    for (int b = a + 1; b < db.size(); ++b) dij.push_back(thermo_Dbinary(db.species[a], db.species[b], T, 101325.0));
                o += std::string(first ? "" : ", ") + "{\"T\": " + jnum(T) + ", \"mu\": " + jnum(mu) + ", \"lam\": " + jnum(lam)
                     + ", \"Xreal\": " + jarr(Xr, jnum) + ", \"mu_i\": " + jarr(mi, jnum) + ", \"lam_i\": " + jarr(li, jnum)
                     + ", \"eta_ij\": " + jarr(eij, jnum) + ", \"D_ij\": " + jarr(dij, jnum) + "}";
                first = false;
            }
            o += "]";
        }
        o += "}";
        std::cout << o << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cout << "{\"ok\": false, \"error\": " << js(e.what()) << "}" << std::endl;
        return 2;
    }
}
