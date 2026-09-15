# 燃焼ガス組成を CEA (NASA-9) の各種で、モル分率で指定できるようにする (擬似種 MIXDRY はユーザ指定の選択肢に)

## メタ

- **area**: `thermophysics / tooling (design chain)`
- **status**: `draft`
- **related_docs**:
  - `methods/thermophysics.md` (多成分 TP 熱物性: 実装 §5 設定)
  - `methods/design/overview.md` (設計チェーンの species 分割 `evaluate.tp_species`)
  - `methods/condensation.md` 実装 §7b (carrier を擬似種に畳む運用)
- **related_plans**: [tooling-nozzle-tp-split-h2o-condensation.md](../accepted/tooling-nozzle-tp-split-h2o-condensation.md) (MIXDRY の起源),
  [chemistry-finite-rate-h2.md](chemistry-finite-rate-h2.md) Phase 0 (CEA `thermo.inp` → `species_db.yaml` 変換ツール), [condensation-kantrowitz-carrier.md](../accepted/condensation-kantrowitz-carrier.md) (多成分 carrier の衝突項)
- **created**: `2026-09-15`
- **owner**: Claude (session 011spCYH)

## 1. 目的

現状、設計チェーン (`forge_design`) は燃焼ガス組成 (問題 YAML `gas.species`, **質量分率**) を CFD 用に
(a) 全部を 1 つの擬似種 `MIX` に畳む (`tp_species: pseudo`) か (b) H₂O 以外を擬似種 `MIXDRY` に畳み H₂O だけ独立種にする
(`tp_species: split_h2o`) かの 2 択で、`species_db.yaml` に NASA-9 係数の質量分率線形混合を書き出している。
擬似種の中身が run の外から見えず (「MIXDRY が何者かわからない」)、モル分率で条件を渡された場合は手で質量分率に換算している。
完了時: (1) 問題 YAML と forge 入力で**組成をモル分率でも指定**でき、(2) **各ガス種を CEA (NASA-9) の係数でそのまま独立種**として
CFD に渡せ (`tp_species: full`)、(3) 擬似種に畳む経路は**名前と中身をユーザが指定**する明示的オプション (`tp_species: lumped`) として残る。
どの経路でも `species_db.yaml` の各エントリに**由来 (CEA 種名 / 擬似種なら構成種とモル分率)** がコメントで残る。

## 2. スコープ

- **やる** (設計チェーン `design/forge_design`):
  - 問題 YAML: `gas.composition_basis: mole | mass` (既定 `mass` = 後方互換) と `gas.species` の解釈。モル→質量換算は `gas/semiperfect.py` の 1 関数に集約し、正規化 (Σ≠1 の入力) をログに出す。
  - `evaluate.tp_species`: `pseudo` (現行) / `split_h2o` (現行, `lumped` の別名に) / **`lumped`** (`tp_lump: {name: MIXDRY, keep: [H2O]}` で畳む種と名前を指定) / **`full`** (全種独立)。
  - `species_db.yaml` の出力: 各種の係数は既存 `SPECIES_NASA9` (CEA `thermo.inp` 転記) から。**`gas.species_db: <path>`** を指定した場合はそのファイル
    (`solver_density_cuda/tools/cea_thermo_to_species_db.py thermo.inp --species ...` で生成した CEA 直読み DB) を優先し、run dir にコピーする。
    擬似種エントリには構成種とモル/質量分率をコメントで書く。
  - 凝縮種 index (`condGasSpecies`) は species 順序から自動決定 (`tp_keep_species` → `gas.condensing_species`, 既定 H2O)。
  - IC (`paste_isentropic_ic`)・BC (`bcondConfig` の `Y{s}`)・`gen_inlet_profile.py` (`--X NAME=EXPR` を追加) を N 種対応に。
- **やる** (forge 本体, 小):
  - `bcondConfig` / `initial` の組成を `X{s}` (モル分率) でも受け付ける (`solverConfig`/`boundaryCond` の読み込みで `Y{s}` に換算。`Y` と `X` の混在はエラー)。
  - 起動ログに species 表 (名前, MW, 入口 Y と X) を出す。
- **やる** (後処理): `axis_csv_va.py` / `total_quantities.py` / ParaView `Forge Saturation` の H₂O index を `physProp.species` から引く (現状 `Y1` 決め打ち)。
- **やらない**: 化学反応 (frozen 組成のまま)、輸送係数の kinetic 混合則の変更、3 温度域 NASA-9、`thermo_d.cu` 内蔵 DB の変更。

## 3. 関連 docs と前提

- forge は既に N 種 TP (`physProp.species: [...]`, `speciesDBFile`) を解き、化学ブランチでは 13 種で運用実績がある。凝縮 carrier 形
  (`condGasSpecies`, Feder 衝突項) も `nSpecies` ループで書かれており 2 種に限定していない (`condensationSourceKernels_d.cuh` L79–90)。
  したがって forge 側の変更は入力の便宜 (`X{s}`) とログだけで、物理・カーネルは触らない。
- `SPECIES_NASA9` (`design/forge_design/gas/semiperfect.py`) は CEA2 `thermo.inp` (McBride–Gordon 2002) からの転記: N2 O2 CO2 H2O AR H2 OH H NO O CO。
  他の種は `cea_thermo_to_species_db.py` で `thermo.inp` (ローカル `.venv-cea/nasa_cea/` にビルド済 CEA2 の同梱ファイル) から生成する。
- 熱力学は質量分率線形混合が厳密なので、`full` と `lumped`/`pseudo` は **dry では同一解** (split_h2o 導入時に軸 M 差 ≤1e-4 を確認済)。
  差は輸送方程式の本数 (2D で +3 本, step 時間 +10–20 % 見込み) と、凝縮 carrier 衝突項が構成種ごとに評価される点 (擬似種では
  MW_mix で代表していた: [condensation-kantrowitz-carrier.md](../accepted/condensation-kantrowitz-carrier.md) §4 の「多成分擬似種の厳密扱い」が `full` で自動的に満たされる)。

## 4. 設計方針

### 4.1 組成の入力と換算 (単一ソース)

問題 YAML:
```yaml
gas:
  model: semiperfect
  composition_basis: mole          # mole | mass (既定 mass)
  species: {H2O: 6.09135e-2, N2: 6.64860e-1, O2: 2.16072e-1, AR: 7.97588e-3, CO2: 4.90034e-2}   # Σ 0.998825 → 正規化
  condensing_species: H2O
  species_db: null                 # 省略時は内蔵 SPECIES_NASA9; パスを書けば CEA 直読み DB を使う
evaluate:
  tp_species: full                 # pseudo | lumped | full (split_h2o は lumped {name: MIXDRY, keep: [H2O]} の別名)
  tp_lump: {name: MIXDRY, keep: [H2O]}
```
換算 $Y_k = X_k M_k / \sum_j X_j M_j$ は `semiperfect.mole_to_mass(X, MW)` に集約し、`GasSemiPerfect` は従来どおり質量分率で動く
(設計 MOC 側の熱力学は不変)。`prepare_info.json` に `gas.X` と `gas.Y` の両方を残す。

### 4.2 species_db.yaml の生成

- `full`: `species:` = YAML の順序 (凝縮種は `condensing_species` の index を `condGasSpecies` に)。各エントリは `SPECIES_NASA9[k]` か `gas.species_db` の同名エントリ。
  `# source: CEA thermo.inp (McBride-Gordon 2002), transcribed in forge_design.gas.semiperfect` のコメント。
- `lumped`: `tp_lump.keep` 以外を `tp_lump.name` に畳む (現行 `mixture_pseudo_species_split` を一般化)。エントリに
  `# lumped: {N2: 0.7089, O2: 0.2304, AR: 0.0085, CO2: 0.0522} (mole fractions within the lump)` を書く。
- `pseudo`: 現行のまま (`MIX`) + 同様のコメント。

### 4.3 forge 側のモル分率入力

`bcondConfig.yaml` の `floats: {X0: .., X1: ..}` と `solverConfig.initial` / IC 生成側で `X{s}` を受けたら、`speciesDBFile` の MW で
`Y{s}` に換算して既存経路に流す (`solverConfig.cpp` の species 読み込み直後)。`Y` と `X` が同じ境界に混在したらエラー。
`inletProfile` CSV の列名も `X_H2O` を許す (`gen_inlet_profile.py` が Y に換算して書く; カーネルは Y のまま)。

### 4.4 後処理の species index

`Y1` 決め打ちのスクリプト (`axis_csv_va.py`, `total_quantities.py` の凝縮警告, ParaView `Forge Saturation` の既定配列名) は
`solverConfig.yaml` の `physProp.species` から凝縮種 index を引く共通関数 (`solver_density_cuda/tools/forge_species.py`) を使う。
ParaView フィルタは配列名をユーザが選べるので既定値のヘルプ文だけ更新する。

## 5. 実装ステップ

1. `methods/thermophysics.md` 実装 §5 に `X{s}` 入力と species ログ、`methods/design/overview.md` に `tp_species` 3 モードと `composition_basis` を追記。
2. `design/forge_design/gas/semiperfect.py`: `mole_to_mass` / `mass_to_mole`、`mixture_pseudo_species_split` の一般化 (`name`, `keep` 任意)、DB コメント。
3. `design/forge_design/probdef.py`: `composition_basis`, `condensing_species`, `species_db` の読み込みと検証。
4. `design/forge_design/evaluate/runner_axismach.py` (+ `runner.py`, `runner_sern.py` の共通部): `_tp_species_list` / `_tp_species_Y` / `_apply_gas_to_config` を 3 モード対応に。
5. `solver_density_cuda/input/solverConfig.cpp`, `boundaryCond.cpp`: `X{s}` → `Y{s}` 換算、species 表ログ。
6. `solver_density_cuda/tools/gen_inlet_profile.py` `--X`、`tools/forge_species.py`、`axis_csv_va.py` / `total_quantities.py` の index 参照。
7. 回帰 (§6)、docs 同期、`design/tests/run_gas_tests.py` に換算と `full`≡`lumped` の等価性テストを追加。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | codex plan レビュー | §4 の入力仕様 (YAML キー名, X/Y の扱い) と §6 の等価性ゲート。**実装はレビュー採否を反映してから** |
| 2 | docs 先行更新 | ステップ 1 |
| 3 | 設計チェーン実装 | ステップ 2–4 |
| 4 | forge 入力 `X{s}` | ステップ 5 |
| 5 | 後処理・ツール | ステップ 6 |
| 6 | 回帰 run と codex result レビュー | §6 |

## 6. 検証

- **単体**: `design/tests/run_gas_tests.py` に (a) mole↔mass 往復 1e-12、(b) va3 の指定モル分率 → 質量分率が README 記載値 (H2O 0.0376954 …) と 1e-6 で一致、
  (c) `full` と `lumped` の混合 cp(T)/h(T) が 200–3000 K で 1e-10 以内 (線形混合の厳密性)。
- **検証ケース** (case/44 va3 M4.19, node Euler TP):
  1. **dry 等価性**: `run_0126` (MIXDRY+H2O, 一様 Tt) を `full` (5 種) で再実行 → 軸 M 差 ≤ 1e-4、ṁ 差 ≤ 1e-4 (split_h2o 導入時と同じゲート)。
  2. **凝縮 (carrier 5 種)**: `run_0131` プロトコル (入口 Tt 分布, cfl 0.5) を `full` で再実行 → onset x 差 ≤ 0.1 r_t、出口 g 平均差 ≤ 2 %
     (Feder 衝突項が構成種ごとになる分の差は `condKantrowitz 1` では出ないはず; `condKantrowitz 2` でも別途記録)。
  3. **モル分率入力**: 同じ run を `composition_basis: mole` の指定で作り、`species_db.yaml`/`bcondConfig.yaml` がビット同一。
  4. **forge `X{s}`**: bcond を `X0..X4` で書いた run が `Y` 版と `res_*.h5` で同一 (ノイズ床以内)。
  5. **CEA 直読み DB**: `cea_thermo_to_species_db.py` で生成した DB を `gas.species_db` に渡し、内蔵転記と係数が一致 (`--check`) → run 同一。
- **判定基準**: 上のゲート。step 時間の増分 (+3 輸送式) を記録 (合否には含めない)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-15` | (実行中: `notes/reviews/2026-09-15-thermophysics-cea-mole-fraction-species-plan.md`) | — | — |

## 7. 影響範囲

- `design/forge_design/{gas/semiperfect.py, probdef.py, evaluate/runner*.py}`, `design/tests/run_gas_tests.py`
- `solver_density_cuda/input/solverConfig.cpp`, `boundaryCond.cpp` (入力換算のみ), `tools/{gen_inlet_profile.py, forge_species.py, total_quantities.py, paraview/forge_filters.py}`
- 既存 problem YAML は `composition_basis` 省略で従来どおり (質量分率)。既存 run の `species_db.yaml` は不変。
- docs: `methods/thermophysics.md` §5, `methods/design/overview.md` (species 分割節), `methods/condensation.md` §7b (擬似種運用は選択肢の 1 つと明記), `procedures/inlet-profile.md` (`--X`)。

## 8. 完了条件

- [ ] 関連 methods を更新済み
- [ ] 実装・§6 の検証 1–5 を満たす
- [ ] codex レビュー 2 回を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status: done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-15` — 初稿 (ユーザ要望: MIXDRY の中身を明示・ユーザ指定可能に、CEA ベースでモル分率指定)。

## 10. 未確定事項

- forge 側 `X{s}` を bcond だけでなく `inletProfile` CSV でも受けるか (§4.3 は受ける前提)。
- `pseudo` (MIX 1 種) を残すか。設計スイープの速度用途があるので残す前提。
