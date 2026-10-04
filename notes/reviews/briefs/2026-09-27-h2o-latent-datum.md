# 諮問ブリーフ: H2O 潜熱を種 DB の気相と液相 H2O(L) の差で作る (#10) — 実装前の設計確認 (2026-09-27)

作業ツリー `/home/sano/work/forge-species` (`feature/species-transport`)。plan: `plans/active/thermophysics-solver-owned-species-db.md` §4.8・§5.1 #10・§6 V7。
AGENTS.md エスカレーション 6 (凝縮カーネルの数値の振る舞いを変える編集の前)。§4.8 は codex plan 段レビュー 2 回を通過済み (§6.1)。

## 問い (1 つ)

下の実装素案で #10 に着手してよいか。**穴があれば最も重い 1 つに絞って直し方を示し、V7 (§6) の合否で足りないものを挙げてほしい。**

## 観測事実

- 現行 `cuda_forge/condensationProperties_d.cuh:239-270` `h2o_latent(T)`: 気相 H2O の CEA 200–1000 K 係数を**再ハードコード**し、T を [COND_T_PROP_FLOOR, 1000] にクランプして多項式を 200 K 未満へそのまま外挿。
  液相 H2O(L) の CEA 273.15–373.15 K 係数を持ち、273.15 K 未満は cp_l(273.15)≈4228 一定の線形外挿、373.15 K 超は 373.15 K の値で頭打ち。L を [1.5e6, 3.5e6] にクランプ。
- 種 DB (`solver_density_cuda/data/species/forge_species_v1.yaml:95-111`) の気相 H2O は同じ係数だが、200 K 未満は種 DB の外挿規約 (thermo_d 側) で評価され、`h2o_latent` と 120 K で 2.39 kJ/kg、150 K で 0.47 kJ/kg ずれる (codex 検算, §4.8)。
  種 DB は `thermoHrefTemp` による datum オフセット (`thermo_d.cuh:51` `h_datum`) を持つ。
- 共通データにはまだ液相 H2O(L) エントリが無い (`phase: condensed` はヘッダの説明だけ)。
- L を使う経路 (§4.8): `cond_latent` (`condensationProperties_d.cuh:299`) → `condensationSourceKernels_d.cuh`・`condensationUpdateLimiter_d.cuh:68-69` (dL/dT を ±0.1 K 差分)・表生成 `condensationTables_d.cuh:106`・float 退避 `condensationSourceF_d.cuh:283`・二相 EOS/音速 `condensationEOS_d.cuh:184`、Python `convert_species_field.py:21,59-68` (湿り場のエネルギー再構築)。
- 気相 MW 0.0180153 と CEA thermo.inp 18.01528 g/mol の相対差 1.1e-6 (yaml deviations)。絶対エンタルピー (−13.4 MJ/kg 級) に掛かると L で ~15 J/kg。

## 素案 (当方)

1. 共通データに `H2O(L)` (`phase: condensed`, `pair_of: H2O`, CEA 273.15–373.15 K 係数, 延長規約 `cp_const_below: 273.15`, `clamp_above: 373.15`) を生成器から追加。外部 DB で凝縮種の気相 H2O を上書きした場合、凝縮 ON なら起動時に拒否 (§4.8 気液ペア契約)。
2. 起動時 (host) に、解決済みの気相 H2O (種 DB の区間・外挿規約・datum オフセット込み) と H2O(L) (**同じ datum オフセット**・気相と同じ MW で質量換算) から、
   `CondSpeciesProps` に液相係数と気相の参照 (区間係数のコピー) を詰める。デバイスの `h2o_latent` は `h_v = 種 DB と同じ評価関数` − `h_l` で作る (再ハードコード撤去)。L のクランプ [1.5e6, 3.5e6] は維持。
3. 表 (`condFloat: 1`) は新 `cond_latent` から再生成 (生成器は関数を呼ぶだけなので変更不要のはず)。dL/dT は現行どおり差分。
4. Python `convert_species_field.py` の `h2o_latent` 移植を、`forge --resolve-species` の記録 (液相係数を含める) から作る実装に変更。
5. 解決済み記録・互換性ハッシュに液相係数・延長規約を入れる (凝縮 OFF の config はバイト不変)。
6. 試験は V7 (a)–(f)。湿潤回帰の基準 run は着手時に固定 (候補 case/44 va3 入口 Tt 分布 noneq `run_0127` 系、AWS)。

## 当方の懸念 (棄却してよい)

- 気相の評価を種 DB と共有すると、200 K 未満の外挿規約が L に入り、onset 付近 (150–200 K) の L が 0.5 kJ/kg 級で変わる。これは「統一分」として記録だけでよいか、それとも onset に効くので合否に入れるべきか。
- 液相を「気相と同じ MW」で質量換算する (CEA の 18.01528 でなく) のは、ペアの一貫性として正しいか。
- float 表 (区間境界 200 K が表の内部に入る) で、種 DB の区間切替 (値は連続・傾きは不連続かも) を表の節点に置く必要があるか (前例: 輸送表は境界ごとに分割)。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `trans.inp`・`thermo.inp` 全体を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの (パスは `/home/sano/work/forge-species/` 基準)

- `sed -n '180,215p;280,310p' plans/active/thermophysics-solver-owned-species-db.md` と §5.1 #10 行 (grep '^| 10 ')
- `sed -n '20,60p;200,320p' solver_density_cuda/cuda_forge/condensationProperties_d.cuh`
- `sed -n '1,140p' solver_density_cuda/cuda_forge/condensationTables_d.cuh`
- `sed -n '270,300p' solver_density_cuda/cuda_forge/condensationSourceF_d.cuh`
- `sed -n '160,200p' solver_density_cuda/cuda_forge/condensationEOS_d.cuh`
- `sed -n '40,80p;300,340p' solver_density_cuda/cuda_forge/thermo_d.cuh`
- `sed -n '1,130p' solver_density_cuda/tools/convert_species_field.py`
- `sed -n '1,30p;90,112p' solver_density_cuda/data/species/forge_species_v1.yaml`
- `sed -n '50,80p;180,210p' solver_density_cuda/tests/unit/test_cond_float.cpp`
