#pragma once

// =============================================================================
// speciesDB.hpp
//   化学種熱物性 DB の host 側解決 (GPU 非依存)。
//   - 内蔵 DB (NASA-9 + Lennard-Jones; 共通データ data/species/forge_species_v1.yaml をビルド時に埋め込む。
//     出典は同ファイルと methods/thermophysics.md) に
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
#include "input/speciesLump.hpp"     // SpeciesLumpSpec (config の lump 指定)

namespace YAML { class Node; }
class solverConfig;

// lump (擬似種) の合成規約 (記録と互換性ハッシュに入る; plan thermophysics-solver-owned-species-db §4.2 #6a)。
//   NASA-9 係数と MW は lump 内モル分率 x_k の線形結合 (固定組成なら cp/h/s° は構成種の和と厳密に一致)。
//   LJ は**暫定**で質量分率の単純平均 (設計側 composition.lump_entry と同じ; plan #7 で実種展開に置き換える)。
//   構成種の温度区切り (Tlo/Tmid/Thi) がすべて同じ場合だけ合成する (区切りの違う種は plan #6b)。
#define SPECIES_LUMP_SYNTHESIS "nasa9 and MW mole-fraction weighted (equal breakpoints only); LJ mass-fraction mean (provisional, plan #7)"

// 解決済み lump の中身 (lump でない種は members が空)。
struct ResolvedLump {
    std::string                basis;          // config の basis ("mole" | "mass")
    std::vector<std::string>   members;        // 構成種名 (config に書いた綴り・順序)
    std::vector<double>        input;          // config に書いた分率そのまま
    std::vector<double>        x;              // lump 内モル分率 (正規化済み; 合成の重み)
    std::vector<SpeciesThermo> memberSpecies;  // 構成種の絶対基準係数 (datum 前)
    std::vector<std::string>   memberSource;   // "builtin" | "file"

    bool empty() const { return members.empty(); }
};

struct ResolvedSpeciesDB {
    std::vector<std::string>   names;    // cfg.speciesNames の順 (index s を定義)
    std::vector<SpeciesThermo> species;  // 同順。datum オフセット前の絶対基準係数 (lump は合成後)
    std::vector<std::string>   source;   // 同順。"builtin" | "file" | "lump"
    std::vector<ResolvedLump>  lumps;    // 同順。lump でない種は空

    int size() const { return static_cast<int>(names.size()); }
    bool isLump(int s) const { return s >= 0 && s < static_cast<int>(lumps.size()) && !lumps[s].empty(); }
    // 種名 → index。大文字小文字を無視 (無ければ -1)。
    int index(const std::string& name) const;
    double MW(int s) const { return species.at(s).MW; }
};

// 内蔵 DB を返す。値は共通データ data/species/forge_species_v1.yaml (ビルド時に埋め込み、起動時に解析) の
// legacy_builtin: solver の種で、キーは canonical ID と別名の両方 (Ar/AR, He/HE, H2O/h2o/WATER, AIR/Air/air)。
// 共通データが壊れていれば std::runtime_error。
std::map<std::string, SpeciesThermo> speciesDB_builtin();

// 埋め込んだ共通データのファイル名と全文の SHA-256 (来歴・ログ用。互換性ハッシュには入れない)。
const std::string& speciesDB_builtinDataName();
const std::string& speciesDB_builtinDataSha256();

// names を内蔵 DB + dbFile (空なら内蔵のみ) で解決する。未知種名・不正 DB は std::runtime_error。
// 名前照合: 完全一致 (外部 DB のキー、内蔵の canonical ID と別名; 外部 DB は同じキーだけを上書き) → 無ければ
// 従来の大小文字無視 (互換; canonical ID への移行と完全一致化は plan #8)。
// 解決結果の names は config に書いた名前のまま (互換性ハッシュに入る)。
ResolvedSpeciesDB speciesDB_resolve(const std::vector<std::string>& names, const std::string& dbFile);

// lump 指定付きの解決 (names に lump の名前も含む; lumps[i].name が names のどれかに一致する)。
//   lump の構成種は内蔵 DB + dbFile から上と同じ規則で解決し、起動時に 1 種へ合成する (SPECIES_LUMP_SYNTHESIS)。
//   拒否 (std::runtime_error): 分率が非正・非有限、構成種が空・重複・未知、lump 名が内蔵種/外部 DB の種名と衝突 (大小文字無視)、
//   basis が mole|mass 以外、構成種の温度区切りが揃わない (plan #6b が必要)。分率の総和が 1 から 1e-3 以上外れたら警告して正規化。
//   凝縮種を構成種に入れる検査は cfg を受ける版 (speciesDB_resolve(cfg)) が行う。
ResolvedSpeciesDB speciesDB_resolve(const std::vector<std::string>& names, const std::string& dbFile,
                                    const std::vector<SpeciesLumpSpec>& lumps);

// cfg.speciesNames / cfg.speciesLumps / cfg.speciesDBFile で解決する。calorically-perfect (species 未指定) では N2 ダミー 1 種。
// 凝縮 ON (condensation: 1) で凝縮種を lump の構成種に入れていたら拒否する。
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

// 起動ログ: 種表 (name, MW, source) と凝縮種・トレーサの状態。lump は中身 (分率)・MW・LJ (暫定)・
// 参照温度 (298.15/1000/2000 K) での cp・h (datum 前の絶対基準) も出す。
void speciesDB_printTable(const solverConfig& cfg, const ResolvedSpeciesDB& db);

// =============================================================================
// 解決済み記録と内容照合 (plans/active/thermophysics-solver-owned-species-db.md §4.3, #3a;
// 仕様 methods/thermophysics.md §1b.4)。
//   - 互換性ハッシュ: 種の順序・名前・相・MW・datum 適用前 (絶対基準) の全係数と温度区間・外挿規約・LJ・
//     thermoHrefTemp と datum 規約・スキーマ版を speciesDB_compatText で正規化 (浮動小数は %.17g) した文字列の SHA-256。
//     lump は合成後の値に加えて、合成規約・構成種の名前・lump 内モル分率 (正規化済み)・構成種の MW・区間・LJ・係数を入れる
//     (#6a)。basis と config に書いた分率そのもの・source は記録にだけ書く。lump の無い DB のテキストは #6a 前とバイト一致。
//     source (builtin/file)・ファイルパスは入れない (来歴として記録にだけ書く)。
//   - 完全性ハッシュ: 記録ファイル resolved_species_<互換16桁>[_<完全性16桁>].yaml 全文の SHA-256。
//   - 各 res_*.h5 (境界出力を含む) のルート属性: species_hash (互換性, 全長) / species_record_sha256 (完全性) /
//     species_record_file / species_input_unverified (0|1)。
//   - Python 側の再計算は tools/forge_species.py (compat_text / load_record)。書式を変えるときは両方を同時に変え、
//     スキーマ版 (SPECIES_RECORD_SCHEMA) を上げる。
//   CPG (thermalMethod != 2) は記録・照合の対象外。液相 (凝縮種の液) は #10 まで含めない。
// =============================================================================
#define SPECIES_RECORD_SCHEMA "forge_resolved_species_v1"
// 外挿規約 (thermo_d.cuh: 区間外は cp を端でクランプ、h は端の cp で線形外挿、T < Tmid で low 係数)。
#define SPECIES_RECORD_EXTRAPOLATION "nasa9_2interval; low if T<Tmid; cp clamped at Tlo/Thi; h linear with end cp outside [Tlo,Thi]"
// datum 規約 (thermo_d.cu: thermoHrefTemp>0 のとき両区間の a7 に -h_abs(Tref)/Ru を加算)。記録の係数は加算前。
#define SPECIES_RECORD_DATUM "coefficients are absolute (before datum); runtime adds -h_abs(Tref)/Ru to a7 of every interval when thermoHrefTemp>0"

struct SpeciesRecordInfo {
    std::string compatHash;     // 互換性ハッシュ (64 桁 hex)
    std::string recordSha256;   // 記録ファイル全文の SHA-256 (64 桁 hex)
    std::string recordFile;     // 記録ファイル名 (run ディレクトリ相対)
    int         inputUnverified = 0;   // 1: 未検証の入力場から開始した (env 許可または入力の印を継承)
};

// SHA-256 (hex 小文字 64 桁)。
std::string speciesDB_sha256Hex(const std::string& bytes);

// 互換性ハッシュの元になる正規化テキスト。Tref は有効な datum 温度 (thermoHrefTemp>0 ならその値、それ以外は 0)。
std::string speciesDB_compatText(const ResolvedSpeciesDB& db, double Tref);
std::string speciesDB_compatHash(const ResolvedSpeciesDB& db, double Tref);

// 記録ファイルの全文。inputStatus は来歴 (verified / unverified_env / unverified_inherited / not_checked_resolve_only)。
std::string speciesDB_recordText(const ResolvedSpeciesDB& db, double Tref, const std::string& dbFile,
                                 const std::string& inputField, const std::string& inputStatus);

// 記録を dir に書く。同名があり全文が一致すればそのまま使い、違えば resolved_species_<互換16>_<完全性16>.yaml に書く
// (既存ファイルは上書きしない)。書き込み失敗は std::runtime_error。
SpeciesRecordInfo speciesDB_writeRecord(const ResolvedSpeciesDB& db, double Tref, const std::string& dbFile,
                                        const std::string& inputField, const std::string& inputStatus,
                                        int inputUnverified, const std::string& dir);

// 入力場の照合。thermalMethod==2 のときだけ呼ぶ。属性の値 (無ければ空文字 / -1) を受け取り、
//   一致 → true (inputStatus を設定)、不一致・照合不能 → メッセージを msg に入れて false。
//   allowUnverified (env FORGE_ALLOW_UNVERIFIED_SPECIES=1) は属性なしの場だけを通す (不一致は通さない)。
//   searchDirs から入力側の記録 (resolved_species_<互換16>*.yaml; 完全性ハッシュが属性と一致するもの) を探し、見つかれば種・キー単位の差を msg に入れる。
bool speciesDB_checkInputField(const ResolvedSpeciesDB& db, double Tref,
                               const std::string& fieldHash, const std::string& fieldRecordSha, int fieldUnverified,
                               const std::string& fieldPath, const std::vector<std::string>& searchDirs,
                               bool allowUnverified, std::string& inputStatus, int& inputUnverified, std::string& msg);

// 記録ファイル (path) と db の種・キー単位の差 (空なら差なし)。記録の完全性・互換性ハッシュの自己整合も検査して差に含める。
std::vector<std::string> speciesDB_diffRecord(const std::string& recordPath, const ResolvedSpeciesDB& db, double Tref);

// ソルバが起動時に書いた記録 (出力 h5 の属性用)。未設定 (CPG・resolve 前) は nullptr。
void speciesDB_setCurrentRecord(const SpeciesRecordInfo& rec);
const SpeciesRecordInfo* speciesDB_currentRecord();
