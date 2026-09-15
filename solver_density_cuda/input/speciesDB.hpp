#pragma once

// =============================================================================
// speciesDB.hpp
//   化学種熱物性 DB の host 側解決 (GPU 非依存)。
//   - 内蔵 DB (NASA-9 + Lennard-Jones; 出典は methods/thermophysics.md) に
//     cfg.speciesDBFile (yaml) を上書き/追加し、cfg.speciesNames の順で解決する。
//   - thermo_init_db (device アップロード) と convertGmshToForge (GPU 無し) の両方が
//     同じ関数を使うので、名前→MW の解決結果が solver と変換器で食い違わない。
//   - bcondConfig の X{s} (モル分率) → Y{s} (質量分率) 換算と検証もここに置く
//     (plans/active/thermophysics-cea-mole-fraction-species.md §4.3)。
// =============================================================================

#include <string>
#include <vector>
#include <map>

#include "cuda_forge/thermo_d.cuh"   // SpeciesThermo (host では inline 関数のみ)

namespace YAML { class Node; }
class solverConfig;

struct ResolvedSpeciesDB {
    std::vector<std::string>   names;    // cfg.speciesNames の順 (index s を定義)
    std::vector<SpeciesThermo> species;  // 同順。datum オフセット前の絶対基準係数
    std::vector<std::string>   source;   // 同順。"builtin" | "file"

    int size() const { return static_cast<int>(names.size()); }
    // 種名 → index。大文字小文字を無視 (無ければ -1)。
    int index(const std::string& name) const;
    double MW(int s) const { return species.at(s).MW; }
};

// 内蔵 DB を返す (キーは内蔵の別名込み: AR/Ar, HE/He, H2O/h2o/WATER, AIR/Air/air)。
std::map<std::string, SpeciesThermo> speciesDB_builtin();

// names を内蔵 DB + dbFile (空なら内蔵のみ) で解決する。未知種名・不正 DB は std::runtime_error。
// 名前照合は完全一致を優先し、無ければ大文字小文字無視で照合する。
ResolvedSpeciesDB speciesDB_resolve(const std::vector<std::string>& names, const std::string& dbFile);

// cfg.speciesNames / cfg.speciesDBFile で解決する。calorically-perfect (species 未指定) では N2 ダミー 1 種。
ResolvedSpeciesDB speciesDB_resolve(const solverConfig& cfg);

// cfg.read() 直後に呼び、解決結果をプロセス内に保持する (thermo_init_db / readBcondConfig が参照)。
// 失敗はメッセージを出して exit する。
const ResolvedSpeciesDB& speciesDB_init(const solverConfig& cfg);

// speciesDB_init で保持した DB。未解決なら nullptr。
const ResolvedSpeciesDB* speciesDB_current();

// モル分率 X → 質量分率 Y: Y_k = X_k M_k / Σ_j X_j M_j。X, MW は同じ長さ。Σ X_j M_j <= 0 は std::runtime_error。
std::vector<double> speciesMoleToMass(const std::vector<double>& X, const std::vector<double>& MW);
// 質量分率 Y → モル分率 X: X_k = (Y_k/M_k) / Σ_j (Y_j/M_j)。ログ用。
std::vector<double> speciesMassToMole(const std::vector<double>& Y, const std::vector<double>& MW);

// bcondConfig の floats ノードから入口組成 Y{s} (db.size() 個) を作る。
//   - X{s} が 1 つでもあれば全種必須。double で読み、負値・非有限・総和 0・未知 index はエラー。
//   - X{s} と Y{s} の混在はエラー。
//   - Y{s} だけのとき: 未指定は既定 (s==0: 1, 他: 0) で補完し、負値・非有限・|ΣY−1|>1e-3 はエラー。
//   - どちらも無ければ空 vector (呼び手が既定補完する)。
// エラーは std::runtime_error (境界名 bname を含む)。
std::vector<double> bcondSpeciesMassFractions(const YAML::Node& floats, const ResolvedSpeciesDB& db,
                                              const std::string& bname);

// 起動ログ: 種表 (name, MW, source) と凝縮種・トレーサの状態。
void speciesDB_printTable(const solverConfig& cfg, const ResolvedSpeciesDB& db);
