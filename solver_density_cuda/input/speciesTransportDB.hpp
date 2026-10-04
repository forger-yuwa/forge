#pragma once

// =============================================================================
// speciesTransportDB.hpp
//   種ごとの輸送物性の出所 (config physProp.transport) の解決・検査・記録用テキスト (host, GPU 非依存)。
//   plans/active/thermophysics-solver-owned-species-db.md §4.3b・§4.3c・§5.1 #5t2 (段 1)、仕様 methods/thermophysics.md。
//
//   physProp:
//     species: [{name: MIXDRY, lump: {N2: 0.7, O2: 0.3}, basis: mole}, H2O]
//     transport: {N2: cea, O2: cea, H2O: "custom:h2o_iapws_cea_v1"}    # 実種 (lump は構成種名) ごとに必須
//
//   モデル: cea (CEA trans.inp の種別フィット; 埋め込みデータ data/species/forge_transport_v1.yaml)
//           kinetic (LJ + Chapman–Enskog、双極子があれば Brokaw 補正、λ は修正 Eucken)
//           fit (speciesDBFile のエントリの transport_fit: {V: [[Tlo,Thi,A,B,C,D],...], C: [...]}; CEA と同じ形・単位)
//           custom:h2o_iapws_cea_v1 (H2O 専用; IAPWS 希薄気体 ↔ CEA を 500–700 K で接続)
//   式と規約は cuda_forge/transportMix_d.cuh。
//
//   段 1 の約束: physProp.transport が書かれているときだけ解決し、ResolvedSpeciesDB::transport に入れる
//   (記録と互換性ハッシュに輸送ブロックを追記; 書かれていない config の記録はバイト不変)。
//   段 2 (#5t2-2): thermo_init_db (cuda_forge/thermo_d.cu) が実種表・組の表・展開行列を device へ上げ、viscMethod 2 の
//   セル・壁が transport_mix_Y で使う (viscMethod 0/1 と physProp.transport の併用は main.cpp が起動を止める)。
// =============================================================================

#include <string>
#include <utility>
#include <vector>

#include "cuda_forge/transportMix_d.cuh"

// 記録・互換性ハッシュに入れる規約の版 (変えたら版を上げる)
#define TRANSPORT_RECORD_SCHEMA  "forge_transport_v1"
#define TRANSPORT_MIXING_RULE    "cea_frozen_v1 (cea2.f TRANP: phi_ij=2 M_j eta_i/((M_i+M_j) eta_ij), psi_ij=phi_ij(1+2.41(M_i-M_j)(M_i-0.142M_j)/(M_i+M_j)^2), phi_ii=psi_ii=1)"
#define TRANSPORT_ETA_IJ_RULE    "v1: both kinetic -> binary Chapman-Enskog (sigma=(si+sj)/2, eps=sqrt(ei ej), no polar term); else CEA trans.inp interaction; else CEA rigid sphere (cea2.f 5565-5570, i=a j=b)"
#define TRANSPORT_FIT_RULE       "ln f = A lnT + B/T + C/T^2 + D; V: f=eta[microPoise], C: f=lambda[microW/(cm K)]; interval = first with T<=Thi else last (cea2.f TRANIN kt)"
#define TRANSPORT_KINETIC_RULE   "mu=2.6693e-6 sqrt(M[g/mol] T)/(sigma^2 Omega22), Neufeld 1972 Omega22 with T* clamped to [0.3,100], Brokaw polar +0.2 delta*^2/T* (delta*=mu_d^2/(2 eps sigma^3), CGS); lambda=mu(cp+1.25 R) (modified Eucken)"
#define TRANSPORT_H2O_V1_RULE    "custom:h2o_iapws_cea_v1: T<=500 IAPWS dilute (mu0 IAPWS 2008 H=1.67752,2.20462,0.6366564,-0.241605; lambda0 IAPWS 2011 L=2.443221e-3,1.323095e-2,6.770357e-3,-3.454586e-3,4.096266e-4; Tc=647.096); T>=700 CEA; 500<T<700 ln f=(1-w)ln f_IAPWS+w ln f_CEA, w=3s^2-2s^3, s=(T-500)/200; T<253.15 f(253.15)(T/253.15)^n, n=dln f/dln T of IAPWS at 253.15 (IAPWS formal range T>=253.15)"

struct ResolvedTransport {
    bool enabled = false;                          // physProp.transport が書かれていたとき true

    // 実種 (lump 展開後、同じ実種は 1 つ; index r)。順序は「physProp.species の順に、lump は構成種の順」で初出順。
    std::vector<std::string>        realName;      // 実種の同一性キー (内蔵は canonical ID、外部 DB は大文字化したキー)
    std::vector<std::string>        realConfigKey; // physProp.transport に書かれたキー (綴りそのまま)
    std::vector<std::string>        modelName;     // cea | kinetic | fit | custom:h2o_iapws_cea_v1
    std::vector<std::string>        dataSource;    // 係数の出所 (例 "trans.inp H2O (SENGERS & WATSON ...)", "speciesDBFile transport_fit")
    std::vector<std::string>        thermoFrom;    // KINETIC の c_p に使う熱物性の位置 (例 "species[1]", "species[0].lump[2]")
    std::vector<double>             dipoleDebye;   // KINETIC: 双極子モーメント [D] (0 = 非極性)
    std::vector<SpeciesTransportD>  sp;            // 評価用 (transportMix_d.cuh)

    // 異種の組 (上三角 a<b, transport_pair_index の順)
    std::vector<TransportPairD>     pairs;
    std::vector<std::string>        pairSource;    // "chapman_enskog" | "trans.inp A/B (ref)" | "rigid_sphere"

    // 展開行列: 輸送する種 s (physProp.species の順) → [(実種 r, lump 内モル分率 x)] (lump でない種は [(r, 1)])
    std::vector<std::vector<std::pair<int, double>>> expand;

    // 互換性ハッシュ・記録に入れる正規化行 (resolver が作る; speciesDB.cpp が本文に追記する)
    std::vector<std::string>        compatLines;
    // 起動ログ用の注意 (適用域外など)
    std::vector<std::string>        notes;

    int nReal() const { return static_cast<int>(sp.size()); }
};

struct ResolvedSpeciesDB;   // input/speciesDB.hpp

// db の輸送指定を解決して db.transport に入れる (spec が空なら何もしない = enabled false)。
//   spec: physProp.transport の (キー, モデル名) を config に書いた順で。dbFile: physProp.speciesDBFile。
//   拒否 (std::runtime_error): 実種に指定が無い / 指定が実種に対応しない (lump 名を書いた・未知・同じ実種に 2 回) /
//     未知のモデル名 / custom: を対象外の種に使った / cea で trans.inp に V・C の両方が無い / kinetic で LJ が無い /
//     fit で speciesDBFile に transport_fit が無い・不正 / 実種数が TRANSPORT_MAX_REAL_SPECIES を超える。
void speciesTransportDB_resolve(ResolvedSpeciesDB& db, const std::vector<std::pair<std::string, std::string>>& spec,
                                const std::string& dbFile);

// 埋め込んだ輸送データのファイル名・SHA-256・元の trans.inp の SHA-256 (来歴・ログ用)
const std::string& speciesTransportDB_dataName();
const std::string& speciesTransportDB_dataSha256();
const std::string& speciesTransportDB_transInpSha256();

// 混合物の μ [Pa s]・λ [W/(m K)] (host)。Xs: 輸送する種のモル分率 (db.size() 個; lump を含む)。
//   Xreal (任意): 展開後の実種モル分率の出力 (nReal 個)。mu_i/lam_i (任意): 実種ごとの単成分値。
void speciesTransportDB_mixture(const ResolvedTransport& tr, const std::vector<double>& Xs, double T,
                                double& mu, double& lam, std::vector<double>* Xreal = nullptr,
                                std::vector<double>* mu_i = nullptr, std::vector<double>* lam_i = nullptr);
