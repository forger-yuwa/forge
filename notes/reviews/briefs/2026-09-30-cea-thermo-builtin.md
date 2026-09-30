# 諮問ブリーフ: 熱物性 (NASA-9) を CEA thermo.inp から生成してソルバ内蔵にする — 設計の確認 (2026-09-30)

作業ツリー `/home/sano/work/forge-species` (`feature/species-transport`, HEAD 8d175ecb 以降)。plan: `plans/active/thermophysics-solver-owned-species-db.md` (§4.1・§4.2・§5.1 #4/#5/#6b・§6 V2/V3/V3f)。
AGENTS.md エスカレーション 1 (§4 の設計方針を新規に書く) と 6 (`cuda_forge/thermo_d` の `SpeciesThermo`・float 表を変える)。
ユーザは diagnostician (Fable) の意見を指名した (codex スイッチは立っているが、ユーザ指定を優先)。

## 問い (1 つ)

ユーザ決定 (2026-09-30): 「熱物性も CEA thermo.inp から生成し、輸送物性の 66 種と揃える (区間数可変 #6b を先に)」。
下の素案で plan §4 を書き実装に進んでよいか。**穴があれば最も重い 1 つに絞って直し方を示し、段の分け方と各段の事前に固定すべき合格条件を示してほしい。**

## 観測事実

- 共通データ `solver_density_cuda/data/species/forge_species_v1.yaml` は 13 種 (気相 12 + `H2O(L)`)。ソルバは `legacy_builtin: [solver]` の 7 種 (N2/O2/CO2/H2O/Ar/He/AIR) だけを内蔵として読む (`input/speciesDB.cpp:117-119`)。
  CO/H2/OH/H/NO/O は `legacy_builtin: [design]` のみで、ソルバでは外部 DB (`species_db_external.yaml`) で渡す必要がある (SERN が 2026-09-30 に指摘、当方の依頼文の誤り)。`legacy_builtin` は plan #4 で「過渡欄」とした。
- 輸送物性は `forge_transport_v1.yaml` に CEA trans.inp から 66 種・相互作用 41 組 (生成器 `tools/cea_trans_to_forge_transport.py`)。
- `SpeciesThermo` (`cuda_forge/thermo_d.cuh:42-54`) は 2 区間固定 (`Tlo/Tmid/Thi`, `low[9]/high[9]`)。共通データの解析も 2 区間を強制 (`speciesDB.cpp:124-127`)。float 表 `SpeciesThermoF` あり。CEA は多くの種で 3 区間 (200–1000–6000–20000 K)、区切りが種で違うものもある。
- 既存の CEA→species_db 生成器 `tools/cea_thermo_to_species_db.py` があり、LJ の表が Python `LJ_PARAMS` と 6 種で食い違う (plan #5 行)。
- 現行共通データと thermo.inp 直読みの既知差 (`deviations`): H2O MW 0.0180153 (thermo.inp 18.01528, 1.1e-6)、He MW (5e-7)、Ar 高温区間が単原子理想の低温係数の流用 (thermo.inp と別物)、内蔵 `AIR` (cp/R 3.5) と CEA `Air` の別名衝突。
- 互換性ハッシュは内容のみ (出所・パスを含まない)。既存 run のハッシュ・記録のバイト不変を試験で守ってきた (`test_species_data_bitexact.py`, `test_species_record_*`)。
- 気液ペア (H2O/H2O(L)) の契約: 凝縮 ON で外部 DB の気相 H2O は内蔵とビット一致のときだけ許可 (#10)。
- lump は起動時に構成種から合成 (#6a, 区切りが揃う種のみ; 区切りの違う種を畳むのが #6b)。

## 素案 (当方)

1. **生成器**: `tools/cea_thermo_to_forge_species.py` で thermo.inp から気相全種 (少なくとも trans.inp の 66 種を含む) を `forge_species_v2.yaml` (または v1 を置き換え) に生成。区間は CEA のまま (1〜3 区間、区切り温度も CEA のまま)。LJ は輸送側 (trans.inp は LJ を持たない) の出典を 1 つに決めて付ける (#5)。既存 13 種の手保守エントリ (H2O(L), AIR, deviations) は生成物に上書きマージ。
2. **区間可変 (#6b)**: `SpeciesThermo` を `nInt` + `Tbound[MAXI+1]` + `coef[MAXI][9]` (MAXI=3 程度) に。デバイスの h/cp/s 評価・float 表・datum・温度反転・lump 合成 (区切りの和集合) を追随。上限超過は起動時拒否。
3. **内蔵の扱い**: `legacy_builtin` の絞り込みを撤去し、共通データの全気相種をソルバ内蔵に。外部 DB による上書きは残す (凝縮の気液ペア契約は維持)。
4. **既存 run との互換**: 既存 7 種の値は今回変えない (deviations の寄せは #5 で別判断)。新たに内蔵になる 6 種は、SERN が外部 DB で渡している係数と同一ならハッシュ不変 (内容ハッシュのため)。3 区間目を足した種は内容が変わるのでハッシュが変わる → 既存 run の種は 2 区間のまま据え置くか、全種を CEA 3 区間に揃えて一斉にハッシュを変えるか、が論点。
5. **試験**: V2 (合成と係数)、V3/V3f (区間可変の double/float)、既存 7 種のビット一致、生成器の thermo.inp 往復 (全種・全係数ビット一致)、Python `semiperfect.py`/`composition.py` との一致。

## 当方の懸念 (棄却してよい)

- 4 の論点: 既存 run の再現 (ハッシュ不変) と「CEA そのもの」を両立できない。ユーザの目的は「CEA が正解、ソルバが持つ」。
- `AIR` 別名と CEA `Air` の衝突 (#5)。CEA 全種を入れると名前空間が広がり、既存 config の種名 (大小文字違い・別名) の解決が変わる危険。
- GPU 構造体の拡大 (`SpeciesThermo` を値渡し・定数メモリに置いている箇所) によるレジスタ・定数メモリ上限 (前例: 潜熱の値渡しでレジスタ 126→166 で起動不能)。
- 20000 K 区間を有効にするか (#5 の未決)。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`、`thermo.inp`・`trans.inp` 全体を読まない (grep で該当行だけは可)。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの (パスは `/home/sano/work/forge-species/` 基準)

- `sed -n '40,130p;225,245p;280,300p' plans/active/thermophysics-solver-owned-species-db.md`
- `sed -n '1,60p' solver_density_cuda/data/species/forge_species_v1.yaml` と grep で各種の `legacy_builtin`/`deviations`
- `sed -n '30,120p' solver_density_cuda/cuda_forge/thermo_d.cuh`、grep で `SpeciesThermo` の利用箇所 (値渡し・`__constant__`)
- `sed -n '90,200p' solver_density_cuda/input/speciesDB.cpp`
- `sed -n '1,80p' solver_density_cuda/tools/cea_thermo_to_species_db.py`、`sed -n '1,60p' solver_density_cuda/tools/cea_trans_to_forge_transport.py`
