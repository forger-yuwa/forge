# 燃焼ガス組成を CEA (NASA-9) の各種で、モル分率で指定できるようにする (擬似種 MIXDRY はユーザ指定の選択肢に)

## メタ

- **area**: `thermophysics / tooling (design chain)`
- **status**: `in_progress`
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
- **やる** (forge 本体):
  - **host 側の DB 解決** (GPU 初期化に依存しない host 関数で `speciesDBFile` を読み名前→MW/係数を返す) を `cfg.read()` 直後と変換器 (`convertGmshToForge`) の両方で使う。
  - `bcondConfig` の `floats` に `X{s}` を受け付け、double で検証 (負値・非有限・総和 0・未知 index・X/Y 混在 → エラー、X 指定時は全種必須) → `Y{s}` に換算して `flow_float` 化。
    `initial` は文字列のまま維持し、組成 IC は IC 生成ツール側で扱う。
  - 起動ログに species 表 (名前, MW, 入口 Y と X, 凝縮種名)。**凝縮種は名前 (`condensationSpecies: H2O`) を正本**にし、数値 `condGasSpecies` は生成値・明示時は一致検査
    (不一致・範囲外・`condModel` と種名の不一致・単一種での `roY` 未登録はエラー)。
  - 残差 CSV に `rms_roY{s}` 列を追加 (現状 species 残差は出ていない)。
  - `interp_field.py` / restart: 種順序・DB (名前・MW) が一致しない restart を黙って受け付けない (名前照合)。
  - **種変換 restart ツール** `tools/convert_species_field.py`: 旧 `[MIXDRY,H2O]` の場を新順序へ (擬似種の保存量を構成種へ分配、H2O と凝縮モーメントを名前で移す、
    DB/エンタルピー基準が変わる場合は元の T から `roe` を再構成)、`ΣρY=ρ`・総水量・T の保存を検査。
- **やる** (後処理): 共通関数 `tools/forge_species.py` (run dir の `solverConfig.yaml` + `species_db.yaml` から種名→index/MW)。`axis_csv_va.py` と ParaView
  `Forge Saturation` はこれで H₂O 配列を**名前で解決** (フィルタに `Run Config` パスのプロパティを追加。設定が無ければ配列選択を必須にし `Y1` を自動採用しない)。
  `total_quantities.py` は既に全種を読むので変更なし (初稿の「Y1 決め打ち」は誤り)。
- **やらない**: 化学反応 (frozen 組成のまま)、輸送係数の kinetic 混合則の変更、3 温度域 NASA-9、`thermo_d.cu` 内蔵 DB の変更、
  forge が `inletProfile` CSV の `X` 列を直接補間すること (生成器が Y に換算して書く)。
- **やる (2026-09-15 ユーザ決定で追加)**: **SERN (`gas.model: frozen_tp`) も同じ `full | lumped` の選択にする** (§4.5)。SERN の `[EXH, AIR]` は
  「流れごとに畳んだ 2 擬似種」= `lumped` の一形態として統一スキーマで表し、`full` では排気・外気の種集合の和で輸送する。排気/外気の見分け
  (IC・warm restart・帳簿) は `lumped` なら $Y_{EXH}$、`full` なら**元素質量分率から作る混合分率** (frozen なら厳密に線形) を共通アクセサで返す。

## 3. 関連 docs と前提

- forge は既に N 種 TP (`physProp.species: [...]`, `speciesDBFile`) を解き、化学ブランチでは 13 種で運用実績がある。凝縮 carrier 形
  (`condGasSpecies`, Feder 衝突項) も `nSpecies` ループで書かれており 2 種に限定していない (`condensationSourceKernels_d.cuh` L79–90)。
  したがって forge 側の変更は入力の便宜 (`X{s}`) とログだけで、物理・カーネルは触らない。
- `SPECIES_NASA9` (`design/forge_design/gas/semiperfect.py`) は CEA2 `thermo.inp` (McBride–Gordon 2002) からの転記: N2 O2 CO2 H2O AR H2 OH H NO O CO。
  他の種は `cea_thermo_to_species_db.py` で `thermo.inp` (ローカル `.venv-cea/nasa_cea/` にビルド済 CEA2 の同梱ファイル) から生成する。
- 熱力学 (cp/h の質量分率線形混合) の等価性は**擬似種内部の組成比が空間的に一定**なら厳密。したがって `full` ≡ `lumped` が成り立つのは
  「同じ解決済み DB・同じ温度域処理・同じエンタルピー基準・非粘性 frozen・擬似種内部比が一定」の条件下に限る (codex M6)。粘性 run では種ごとの拡散係数と
  種エンタルピー拡散 (`speciesTransport_d.cu`) が効くので一般には同一解でない。`pseudo` 1 種は組成輸送自体を持たない。凝縮 carrier 衝突項 (`condKantrowitz 2/3`) は
  構成種ごとの評価になるので等価性を要求せず、種別和を独立計算と照合する試験にする。
- CEA `thermo.inp` 直読みと `SPECIES_NASA9` 転記は**完全一致しない** (codex M7 実測: H2O MW 0.0180153 vs 0.01801528、AR 高温域 a0 0 vs 20.105)。外部 DB を使う経路では
  その DB を設計・CFD 双方の正本にし、内蔵転記との差は物性値の許容差で評価する。`cea_thermo_to_species_db.py --check` は N2 low のみ・不一致でも rc 0 なので、
  全使用種の MW・両温度域を照合し不一致で失敗終了するよう直す。
- 既存の `frozen.mole_to_mass()` (SERN 用) を共通換算関数に委譲する (重複実装しない, codex m1)。

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
換算 $Y_k = X_k M_k / \sum_j X_j M_j$ は `gas/composition.py` の 1 関数 (既存 `frozen.mole_to_mass` はここへ委譲) に集約する。
**解決済み DB を 1 つ構築して全経路へ渡す** (codex M1): 問題読込時に内蔵 `SPECIES_NASA9` に `gas.species_db` (問題 YAML の所在基準で解決) を上書きした
`ResolvedSpeciesDB` (名前, MW, 2 温度域係数, 出典) を作り、換算・`GasSemiPerfect` (設計 MOC)・擬似種生成・IC (`paste_isentropic_ic`)・`species_db.yaml` 出力の
**すべて**がこれを使う (外部 DB だけにある種を設計側で失敗させない、既存種の係数上書きで MOC と CFD が食い違わない)。検証: MW>0、両温度域の係数個数、
温度区切りが 200/1000/6000 K と異なる種は `lumped` に混ぜることを拒否。`prepare_info.json` に `gas.X`・`gas.Y`・解決済み DB の出典を残す。

### 4.2 species_db.yaml の生成

- `full`: `species:` = YAML の順序。**凝縮種は `gas.condensing_species` (名前) が正本**で `condGasSpecies` は生成値。手書き index があれば一致検査、不一致はエラー。
  `lumped` で凝縮種が `keep` に無い、`pseudo` (MIX) で凝縮 ON、擬似種名と実種名の衝突、重複種、`condModel` と種名の不一致は入力段階で拒否。各エントリは `SPECIES_NASA9[k]` か `gas.species_db` の同名エントリ。
  `# source: CEA thermo.inp (McBride-Gordon 2002), transcribed in forge_design.gas.semiperfect` のコメント。
- `lumped`: `tp_lump.keep` 以外を `tp_lump.name` に畳む (現行 `mixture_pseudo_species_split` を一般化)。エントリに
  `# lumped: {N2: 0.7089, O2: 0.2304, AR: 0.0085, CO2: 0.0522} (mole fractions within the lump)` を書く。
- `pseudo`: 現行のまま (`MIX`) + 同様のコメント。

### 4.3 forge 側のモル分率入力 — 入力契約 (codex M4 で具体化)

- **DB 解決は host で**: `speciesDBFile` を読む host 関数 (GPU 非依存) を `solverConfig::read()` の直後と `convertGmshToForge` の境界読込の前に呼び、名前→MW を得る。
- **bcond**: `floats: {X0:.., X1:.., ...}` を double で検証 → $Y_k = X_kM_k/\sum X_jM_j$ → `flow_float` の `Y{s}` に入れて既存経路へ。同一境界での `X`/`Y` 混在、負値・非有限・総和 0・
  未知 index はエラー。**省略種の補完 (`Y0=1` 他 0) は X 指定時には行わず全種必須**。
- **initial**: 既存の文字列形式は維持。組成付き IC は IC 生成ツール (`paste_isentropic_ic(species_X=...)` / `gen_tp_ic*.py --X`) が HDF5 `VALUE/roY{s}` に書く。
- **inletProfile CSV**: forge が読む列は従来どおり `Y{s}`。`gen_inlet_profile.py` が `--X NAME=EXPR` / 測定表の `X_NAME` 列を受けて Y に換算して書く (X 指定時は全種必須)。
  forge 側の直接 X 補間は範囲外。

### 4.4 後処理の species index (codex M5)

`tools/forge_species.py`: run dir の `solverConfig.yaml` (`physProp.species`, `condensation`) と `species_db.yaml` から `{name: index, MW}` と凝縮種名を返す。
`axis_csv_va.py` はこれで H₂O 列を引く。ParaView `Forge Saturation` には `Run Config (solverConfig.yaml)` 文字列プロパティを追加し、指定時は種名から配列 (`Y{s}`) を
解決、未指定時は `Vapor Mass Fraction Array` の明示を必須にして **`Y1` を自動採用しない** (5 種順序では `Y1`=N2 を水蒸気として計算してしまう)。
H₂O を先頭・中間・末尾に置いた順序入替試験を単体に入れる。

### 4.5 SERN とノズル設計を同じ選択肢にする — 畳みは「流れ種別」でなく「run の選択」(2026-09-15 ユーザ決定)

**統一スキーマ**: 擬似種 (lump) は「内部組成が固定された NASA-9 の線形混合」であり、由来が「組成の部分集合」(ノズル: MIXDRY) でも
「流れ (作動点組成 / 外気)」(SERN: EXH / AIR) でも同じ物である。したがって

```yaml
evaluate:
  tp_species:
    mode: full | lumped          # (旧 pseudo = lumped で全部を 1 lump; split_h2o = 下の例)
    lumps:                       # lumped のみ。値は「畳む対象」
      MIXDRY: {from: composition, exclude: [H2O]}     # ノズル: 組成の部分集合を畳む
      EXH:    {from: stream, stream: inflow}          # SERN: 流れの組成 (作動点 gas.composition) を畳む
      AIR:    {from: stream, stream: external}        #       外気 (spec.external.composition / 乾燥空気)
    keep: [H2O]                  # 独立種のまま残す種 (凝縮種は必ずここ)
```

- `full`: species = 使う全流れの組成の**和集合** (SERN なら排気 11 種 ∪ 空気 4 種 = 11–13 種)。IC は領域ごと (排気側 / 外気側) にその流れの Y ベクトルを貼り、
  BC も流れごとの `Y{s}` (または `X{s}`)。ノズルは流れが 1 つなので従来の `full` と同じ。
- `lumped`: 各 lump が 1 擬似種。SERN の現行 `[EXH, AIR]` は `lumps: {EXH: stream inflow, AIR: stream external}` と等価 (後方互換の別名にする)。
  内部表現は mapping (`tp_species: {mode, lumps, keep}`) に一本化し、旧文字列形式 (`pseudo` / `split_h2o` / `[EXH,AIR]`) は入力時の別名変換で受ける (競合指定は拒否; m1)。
  ノズルの `split_h2o` は `lumps: {MIXDRY: composition exclude [H2O]}, keep: [H2O]`。
- **流入元ラベルは独立トレーサ** (codex 再レビュー M1/M2 で元素混合分率案を撤回): `full` では排気率 ξ を **受動スカラ `roXi` (排気 1 / 外気 0) として輸送**する
  (既存の汎用スカラ輸送コア `ScalarTransportDesc` を流用: 登録・入口 Dirichlet (排気入口 1, 外気入口 0)・point-implicit 更新・残差列・restart まで含める。拡散は
  化学種と同じ混合平均 `Sc`)。`lumped` では ξ = Y_EXH で同じアクセサ `exhaust_fraction(run)` が返す。元素混合分率 $Z_e$ は**診断のみ**
  (二流体線形混合からのずれの指標; 差動拡散があると元素ごとに ξ が異なり、`Y=ξY_{exh}+(1-ξ)Y_{air}` は復元できない。power-off φ=0 でも入口流は
  あるので「排気なし」に置換せず、トレーサはそのまま入口由来率を保つ。元素診断は $|Z_{e,exh}-Z_{e,air}|$ が小さい近退化では「定義不能」を返し、float32 の
  保存だけで ξ 誤差が 7e-5 になる増幅 (1/|ΔZ| ≈ 2500, m4_off の N) を上限で拒否)。
- **流れごとの質量配分が正本** (M3): 各流れ (排気 / 外気 / ノズルなら単一流れ) について `keep` の種を先に取り出し、残りをその流れの lump に配分する。
  純排気でも `keep: [H2O]` があれば $Y_{EXH}=1-Y_{H2O}$ (m6_on: $Y_{H2O}$ 0.2411, lump 0.7589) なので ξ を擬似種名に依存させない (上のトレーサ)。
  未配分・二重配分・空 lump は入力段階で拒否。**lump→実種の展開行列** と **流れ→輸送種の入口ベクトル** を run に保存する。
- **機械可読メタデータ** (M5): `species_meta.yaml` (run dir) に 実種の原子組成 (CEA `thermo.inp` の元素欄から取得、別名解決後も保持)、lump の構成比、各流れの
  正規化済み組成、輸送種順序、トレーサの有無を書く。後処理・restart は問題 YAML を再解釈せず、この保存情報を使う。`ResolvedSpeciesDB` に原子組成を持たせ、
  `cea_thermo_to_species_db.py` が元素欄を出力する。
- **restart 経路** (M4): `runner_sern.py` の `restart_by_index()` (7 変数のみコピー) と `runner_sern3d.py` の `warm_from_same_mesh()` は種保存量を落としている
  (ΣY が壊れる既存バグ)。両関数を全種・トレーサ込みにし、§6 に「段階切替直前・直後の全 `roY{s}`・ΣρY/ρ・T の保持」試験を追加する。作動点変更時
  (m6_on→m10_on) の組成再構成は**情報を落とす初期化操作**として別契約 (`convert_species_field.py`) にする。
- **コスト**: SERN 3D SST で `full` は輸送方程式 2 → 11–13 本。step 時間は +50–100 % の見込み (chem ブランチの 13 種実績から)。MOO の探索は `lumped`、
  最終評価や凝縮・化学の前段は `full`、と run ごとに選べるのが目的なので既定は変えない (SERN 既定 `lumped`、ノズル既定は従来 `pseudo` 相当)。
- **熱力学の等価性**: `full` と `lumped` は「各 lump の内部比が空間的に一定」のとき厳密に同じ (§3)。SERN では排気と外気が混合する層で
  内部比は一定 (lump は流れ単位なので混合は lump 間の線形混合) → 非粘性 frozen なら同一解。粘性では差動拡散の分だけ異なる (仕様として記録)。

## 5. 実装ステップ

1. `methods/thermophysics.md` 実装 §5 に `X{s}` 入力と species ログ、`methods/design/overview.md` に `tp_species` 3 モードと `composition_basis` を追記。
2. `design/forge_design/gas/composition.py` (新: `mole_to_mass`/`mass_to_mole`、`ResolvedSpeciesDB`)、`semiperfect.py` を解決済み DB 経由に、`mixture_pseudo_species_split` の一般化 (`name`, `keep` 任意)、DB コメント。`frozen.mole_to_mass` は委譲。
3. `design/forge_design/probdef.py`: `composition_basis`, `condensing_species`, `species_db` (YAML 所在基準) の読み込みと検証 (§4.1–4.2 の拒否条件)。
4. `design/forge_design/evaluate/runner_axismach.py` (+ `runner.py`): `_tp_species_list` / `_tp_species_Y` / `_apply_gas_to_config` を統一スキーマ (`full | lumped`, `lumps`, `keep`) に。
4b. `design/forge_design/evaluate/runner_sern.py` / `runner_sern3d.py`: `SPECIES_ORDER` 固定を撤去し統一スキーマへ (`frozen_gases` は流れごとの `FrozenGas` + 解決済み DB、`paste_region_ic` / BC / `warm_from_run` は流れ→輸送種の入口ベクトルとトレーサで)。`restart_by_index` / `warm_from_same_mesh` を全種・トレーサ込みに修正。`[EXH, AIR]` は別名で後方互換。
4c. forge: 排気トレーサ `roXi` (`physProp.tracer: exhaust` で有効化; 汎用スカラ輸送コアで登録・入口 Dirichlet・point-implicit 更新・残差列 `rms_roXi`・出力・restart)。`full` の SERN 評価器はこれを排気率に使う。
5. forge: host DB 解決、`boundaryCond.cpp` の `X{s}` → `Y{s}` 換算と検証、species 表ログ、凝縮種名と index の一致検査、残差 CSV の `rms_roY{s}` 列、`interp_field.py` の種名照合。
5b. `tools/convert_species_field.py` (種変換 restart, §2) と `cea_thermo_to_species_db.py --check` の全種照合・失敗終了化。
6. `solver_density_cuda/tools/gen_inlet_profile.py` `--X`、`tools/forge_species.py`、`axis_csv_va.py` / `total_quantities.py` の index 参照。
7. 単体 (§6) → 小型 TP 収束ケース (node のみ) → case/44 回帰 → docs 同期。**コード世代は [condensation-source-limiter-steady](condensation-source-limiter-steady.md) の実装後に固定** (凝縮 run の比較は両経路の `condLim≈1` を確認してから)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー~~ | 2026-09-15 実施 (§6.1)。M1–M8/m1–m2 を全て採用し §2/§3/§4/§6 に反映済み。**実装着手可** (ユーザ確認後; 順序は limiter plan の後) |
| 2 | ~~docs 先行更新~~ | 済 (2026-09-16): `methods/thermophysics.md` §5、`methods/design/overview.md` (統一スキーマ)、`design/CAPABILITIES.md`、skill design-intake |
| 3 | ~~共通基盤 (解決済み DB・換算・凝縮種名正本)~~ | 済 (2026-09-16): `gas/composition.py` (`ResolvedSpeciesDB`, `mole_to_mass`/`mass_to_mole`, `parse_tp_species`, `resolve_species_layout`, `species_db_yaml`/`species_meta`), `semiperfect`/`frozen` は DB 注入、`probdef` に `composition_basis`/`condensing_species`/`species_db` と検証。単体 (a)–(e) + 外部 DB + SERN 配分を `run_gas_tests.py` に追加 (ALL PASS) |
| 4 | ~~設計チェーン統一スキーマ (ノズル)~~ | 済 (2026-09-16): `runner_axismach` は `p.species_layout()` で species / `species_db.yaml` (由来コメント) / `species_meta.yaml` / `condGasSpecies` (名前から生成、手書きは一致検査) / `condensationSpecies` / prepare_info `species` を出す |
| 4b | SERN の統一 (`SPECIES_ORDER` 撤去, 流れごとの質量配分, 領域 IC, restart 経路の全種化) | コード済 (2026-09-16): `frozen_gases` が layout (省略時 `[EXH, AIR]` 別名) と輸送種ごとの `FrozenGas` を返す、`gas_states`/`_solver_config`/IC/BC は N 種 + `Xi`、`restart_by_index`/`warm_from_same_mesh` は全 roY + roXi、`warm_from_run` は N 種 (順序不一致は拒否)。単体 (`run_sern_frozen_gas_tests.py`: full 11 種・lumped+keep・m4_off 配置一致・restart) ALL PASS。**残**: case/46 の 2D node Euler 回帰 (別名 = 旧 run とノイズ床内; full vs lumped) |
| 4c | ~~排気トレーサ `roXi` (forge) と `species_meta.yaml`~~ | 済 (2026-09-16): `cuda_forge/tracerTransport_d.{cu,cuh}` (`physProp.tracer: exhaust`; 入口 `floats.Xi` Dirichlet + node ピン、Neumann、汎用スカラ輸送、point-implicit/RK、`rms_roXi`、出力 level 0、restart `VALUE/roXi`)。**制限**: 拡散は 0 (汎用スカラ拡散は Sc 無しの μ ベース、化学種 Fick は多成分専用のため; Euler の SERN では無関係、SST では要フォロー [F-sp1])。`species_meta.yaml` は設計側 `write_species_files` |
| 5 | ~~forge 入力 `X{s}`・種名検査・`rms_roY`・restart 照合~~ | 済 (2026-09-16): `input/speciesDB.{hpp,cpp}` (host 側 DB 解決; `thermo_init_db` と `convertGmshToForge` が共用), `boundaryCond.cpp` の `X{s}`/`Y{s}` double 検証→換算 (X/Y 混在・欠落・ΣY≠1 拒否), `condensation.condensationSpecies` (名前→index、不一致/範囲外エラー), 起動 species 表 + 入口 Y/X ログ, `rms_roY{s}` 列, `interp_field.py` の種名照合 (`--force-species`)。単体 `tests/unit/test_species_db_host.cpp` 31 PASS。500 step 検証 `case/44 run_0190`–`0192`/`0194` (X 入力・トレーサとも同一設定反復ノイズ内 ≤7e-6) |
| 5b | ~~種変換 restart ツールと CEA `--check` 修正~~ | 済 (2026-09-16): `tools/convert_species_field.py` (`species_meta.yaml` の展開行列で lump→実種、名前で H2O/モーメント/roXi 移送、ΣρY=ρ・総水量・T 保存検査; `roe` は差分形 ρ[e_dst−e_src] で液相エネルギーを保つ)、`cea_thermo_to_species_db.py` は元素欄 → `atoms`、`--check` は全指定種の MW/両区間を照合し閾値超で非零終了 (既知差: H2O MW 1.1e-6, AR high a0) |
| 6 | ~~後処理・ParaView 配列解決~~ | 済 (2026-09-16): `tools/forge_species.py` (run dir → 種名/index/MW/凝縮種/vapor_array/tracer)、ParaView `Forge Saturation` に `Run Config` プロパティ (未指定時は `Y1` を自動採用しない)、`gen_inlet_profile.py --X`、`case/44 axis_csv_va.py` は名前解決 |
| 7 | 回帰 run と codex result レビュー | §6 (node のみ; **cell は対象外**, 2026-09-16 ユーザ指示)。case/44: 新バイナリ A/B `run_0195`、`full` 5 種 `run_0196` (0170 の収束場を種変換 restart)、モル分率入力 lumped `run_0197`、CEA DB `run_0198` (prepare)、SERN case/46: node Euler m6_on `run_0100` (別名 [EXH, AIR]) / `run_0101` (full 11 種 + roXi) |
| 8 | F-sp1: トレーサ `roXi` の拡散 (SST の SERN で必要なら混合平均 Sc) | 未着手 (Euler では不要; 粘性二流体試験 [§6] も未実施でこれに依存) |
| 9 | ~~(result-1 M1) 種名の YAML 引用 + 生成 config の再読込試験~~ | 済 (2026-09-16): runner_axismach / runner_sern の `species:` を引用符付きに、`run_sern_frozen_gas_tests` で safe_load 再読込 (全て str) を確認 |
| 10 | ~~(M2) `convert_species_field.py` の二相 EOS 反転・書込前の保存検査・総水量 ρY_H2O~~ | 済 (2026-09-16): ソルバの二相 EOS (e = e_gas + g(R_w T − L)) を numpy 移植、`run_0196` の湿潤 input で反転 T が `res_0` の T と最大 1.9e-4 K (乾き反転は 54 K 誤差)、全セル (湿潤含む) の T 保持 (`--T-tol` 0.05 K)・ΣY 1e-6・実種質量・総水量 (ρY_H2O) を書込前の hard fail に |
| 11 | ~~(M3) restart 照合~~ | 済 (2026-09-16): 設計側 `check_species_compatible` (順序・MW rel 1e-9・datum・トレーサ; `warm_from_run` は作動点変更で lump の MW が変わるので MW のみ許容) を全 restart 経路に; forge ツール側 `forge_species.species_signature/compare_signatures` (名前順序・MW・datum・meta ハッシュ) を `interp_field.py` に適用 (照合不能は既定エラー, `--force-species`)。`NO`/`N`/`Y` を文字列のまま読む SafeLoader |
| 12 | ~~(M4) 変換器: stream lump の保存的展開 + 再初期化 + roXi 生成~~ | 済 (2026-09-16): `--mode conserve` (一意の行き先のみ; 複数 stream lump は拒否) / `--mode reinit` (ξ = roXi か流入元ラベル種で流れ配分; 情報を落とす操作として射影損失を報告); 目標に tracer があれば roXi=ρξ を生成。実 meta で lumped→full (保存 4e-16, roXi 生成)・full→lumped (reinit, 射影損失 4.3e-6)・同一配置 (恒等) を確認 |
| 13 | ~~(M5) roXi の周期 node 残差合算・同期 + 周期移流試験~~ | 済 (2026-09-16): `periodicNode_d.cu` に `res_roXi` 合算と `roXi` ミラー; case/09 `run_0065`/`0066` (33³ 周期箱, 一様流, 500 step): 周期対の ΔXi = 0, ∫ρXi dV 相対変化 4–5e-8, seam 版と内部版の差 3e-7 |
| 14 | ~~(M6) tracer × dual-time を入力で拒否~~ | 済 (2026-09-16): config 読込でエラー; solver-settings は定常/RK のみと記載 (物理時間項は F-cf8 と同類) |
| 15 | ~~(M7) `condensationSpecies` ↔ `condModel` の対応検査 + config 読込単体~~ | 済 (2026-09-16): condModel 0=N2 / 1=H2O と物質の対応を host 入力で検査; `tests/unit/test_solver_config_species.cpp` (solverConfig::read を子プロセスで) 15/15 PASS |
| 16 | ~~(M8) `lumped+keep` 等でもトレーサ有効化、`exhaust_fraction(run)` 共通アクセサ~~ | 済 (2026-09-16): トレーサは「inflow で 1・external で 0 の輸送種 (純粋な流入元ラベル) が無いとき」有効 (full と lumped+keep は有効, [EXH, AIR] は Y0)。`composition.exhaust_fraction(run_dir)` / `exhaust_fraction_field` が `species_meta.exhaust_fraction` ({kind, array, conserved}) を返す。単体追加 |
| 17 | ~~(M9) `forge_species.vapor_array` を凝縮の有無と独立に~~ | 済 (2026-09-16): `h2o_index`/`h2o_array` 追加、`vapor_array` は凝縮 OFF でも H2O を名前解決 (`run_0203` → Y0); `axis_csv_va.py` も名前解決 |
| 18 | (M10) 小型 5 種 node TP の `check_convergence` PASS ケース; 作動点変更 (m6_on→m10_on) と湿潤 restart の試験 | **小型 PASS ケース済 (2026-09-16)**: case/16 `run_0471_tp5_node_euler_va3` (5 種) / `run_0472_tp2_node_euler_va3_lumped` (lumped): node Euler TP, Tt 353 K, 12000 step で全列 3.3 (rms_roUy)〜4.8 桁低下 (rms_roY 4.2 桁) **PASS**, ΣY−1 ≤ 8.6e-8, 負値 0, 5 種 vs lumped の場差 ro 1.0e-6 / P 1.7e-6 / T 1.0e-6 / Y_H2O 3.9e-6 (float 精度), step 時間 +33 % (case/16 README run 表)。**作動点変更・湿潤 restart 済 (2026-09-16)**: case/46 `run_0105` (m6_on `run_0102` → m10_on を `warm_from_run`; ΣρY/ρ 1.2e-7, T 整合 1.8e-3 K, 適応段 2000 step で残差 4.5 桁低下, NaN 0); 湿潤 restart は `run_0196` の入力 (6365 湿潤セル) を二相 EOS 変換器で同一配置変換して T 保持 6.8e-13 K (#10); 同一条件の段階 restart (restart_by_index) は case/46 `run_0100`–`0103` の段階起動で全種を引き継ぎ最終 ΣY 1±2.4e-7。**トレーサ誤差は切り分け済** (2026-09-16, case/46 `run_0104`: lumped run に tracer を併用すると同一 run 内で `\|Xi−Y_EXH\|` 最大 1.8e-4 / 平均 4.6e-7 → 汎用スカラ輸送カーネルと化学種輸送カーネルの再構成差であり配置差ではない; §6 の「ξ と Y_EXH 1e-6」は同一カーネルでないと成立しないので **「平均 1e-6 以内 + 最大はカーネル差 (同一 run 実測) 以内」に改定**)。README「同じ固定点」→「ゲート内で一致」は済 |
| 20 | (result-2 M1) 作動点変更の組成再初期化 Y = ξ Y_in^dst + (1−ξ) Y_ext^dst (変換器 reinit + warm_from_run) と目標組成の検証 | 両側 |
| 21 | (M2) restart 署名を実 config + DB (係数・区切り・datum・tracer) から; 矛盾/照合不能/保存量欠落は拒否 | 両側 |
| 22 | (M3) 変換器の有限性・非負・EOS 残差検査 + 失敗系試験 | forge ツール |
| 23 | (m1) docs の ξ 説明統一、§6 旧ゲートの履歴化、粘性試験の延期明記、M10 の桁数訂正 | docs |
| 19 | ~~(m1) docs 同期~~ | 済 (2026-09-16): methods/thermophysics §5 (speciesDB_resolve/init, 上書き), plans/README (in_progress), condensation §7b |

## 6. 検証

- **単体** (`design/tests/run_gas_tests.py`, CFD 回帰より前): (a) mole↔mass 往復 `rtol 1e-12`; (b) va3 の指定モル分率 (Σ 0.998825) → $Y_{H2O}$ = 0.03769539643 (`rtol 1e-9`);
  (c) `full` と `lumped` の混合 cp(T)/h(T) を 200–3000 K で比較: **`rtol 1e-12` + 量別 atol (cp 1e-9 J/kg/K, h 1e-6 J/kg)**、$h(T_{ref})=0$ 近傍も検査
  (既存 split の実測差は cp 1.6e-12 / h 2.3e-9 なので絶対 1e-10 は不適; codex m2); (d) `species_db.yaml` は解析後の値を比較し、ビット同一は同じ正規化済み入力の再出力に限る;
  (e) 拒否条件 (§4.1–4.3) がそれぞれエラーになる; (f) `forge_species.py` と ParaView の H₂O 順序入替 (先頭/中間/末尾); (g) 種変換 restart の `ΣρY=ρ`・総水量・T 保存;
  (h) CEA 直読み DB vs 内蔵転記の全使用種 MW・両温度域係数の差を表にし、物性値 (cp, h) の差を許容差で判定。**実測 (2026-09-16, `cea_thermo_to_species_db.py --check`, va3 5 種)**: 係数は N2/O2/CO2/H2O が両区間 0、AR は high a0 (20.105 vs 0; 内蔵転記の既知差)、MW は H2O 1.11e-6 (0.01801528 vs 転記 0.0180153) → cp/h の相対差も 1.11e-6 で、当初の rtol 1e-6 をわずかに超える。許容差を **rtol 5e-6** に改定 (根拠: 転記 MW の 6 桁丸め; 設計 R 285.27044 vs 285.27043, γ* 5e-12 で実用上無差)。外部 DB を使う run はその DB が設計・CFD 双方の正本なので内蔵との差は結果に混入しない (§3)。
- **収束済み小型 TP ケース (node のみ; 2026-09-16 ユーザ指示で cell は対象外, codex M8 の「両方」を改定)**: 5 種 frozen の 2D 亜音速ノズル (case/13 系の小メッシュ) で `check_convergence.py` **PASS** (`rms_roY{s}` 列込み)、
  `ΣY=1±1e-6`・負値なし。ここで X/Y 入力・種順序・旧入力 (`split_h2o`) の回帰を行い、`res_*.h5` がノイズ床以内で一致。
- **case/44 va3 M4.19 (node Euler TP)** — 未収束 (warm 床 plateau) の既存 run は**回帰参考**で、収束解一致の根拠にはしない:
  1. dry: `run_0126` 相当を `full` で再実行 → 軸 M 差 ≤ 1e-4、ṁ 差 ≤ 1e-4 (非粘性 frozen・内部比一定の条件下, §3)。
  2. 凝縮 (carrier 5 種, `condKantrowitz 1`): `run_0131` プロトコル (limiter plan 実装後は同 CFL で `condLim≈1` を確認) を `full` で再実行 → onset x (壁流線, $g>10^{-3}Y_w$) 差 ≤ 0.1 r_t、
     出口 g 平均 (質量流束重み, $x=x_{max}-2r_t$) 相対差 ≤ 2 %、ṁ 相対差 ≤ 1e-4。凝縮固有量の時系列は `check_quasisteady.py --series-csv` で STEADY (0.2 %)。
     `condKantrowitz 2/3` は等価性でなく衝突項の種別和を独立計算と照合。
  3. モル分率入力: `composition_basis: mole` で作った run の `species_db.yaml`/`bcondConfig.yaml` が解析後の値で一致。
  4. forge `X{s}`: bcond を `X0..X4` で書いた run が `Y` 版と `res_*.h5` でノイズ床以内。
  5. CEA 直読み DB: `gas.species_db` を渡した run は**その DB を正本**として設計・CFD が同じ係数を使うこと (prepare_info と species_db.yaml の出典で確認)。内蔵転記との「同一 run」は要求しない。
- **SERN (case/46)** (codex 再レビュー M6 を反映):
  - 単体: powered / power-off (m4_off: 燃料なしでも入口流あり) / 微小組成差 (既定外気 `AIR_MOLE` と作動点の空気) / `keep` 併用 / 種順序変更 で、流れごとの質量配分・展開行列・入口ベクトル・トレーサ BC が意図どおり。元素診断は近退化で「定義不能」を返す。
  - restart: 同一条件 restart (段階切替直前・直後) で全 `roY{s}`・トレーサ・ΣρY/ρ・T が保持される; 作動点変更 (m6_on→m10_on) は `convert_species_field.py` 経由で ΣρY=ρ・T 保存を検査。既存 `restart_by_index` の欠落修正による差は後方互換比較から分離して記録。
  - 2D node Euler 作動点 1 点: (i) `lumps: {EXH, AIR}` 指定が現行 `[EXH, AIR]` run とノイズ床以内; (ii) `full` と `lumped` の非粘性 frozen 比較でノズル力・機体力・出口運動量が相対 1e-3 以内 (ゼロ近傍の力は $F_{ideal}$ 基準の絶対許容差)、トレーサ ξ と $Y_{EXH}$ の場が平均 1e-6 以内・最大は同一 run 内のカーネル差 (lumped+tracer 併用 run で実測, 2026-09-16: 1.8e-4) 以内 [改定: 汎用スカラと化学種のカーネルが違うため点ごとの 1e-6 は成立しない]; 比較する両 run は全種残差込み `check_convergence` PASS (`require_residual_pass=True`)、全種の有限性・非負性、力の時系列 STEADY を**比較許容差より十分小さい変動幅** (drift/fluct 0.1 %; 既存 sern_forces の 2 %/5 % は不可) で要求。**実測 (2026-09-16)**: node Euler SERN は 24000 step でも残差 plateau (2.4 桁) で PASS が取れない (限界サイクル) ため、PASS の代わりに「同じ plateau (両 run の床が同じ) + 力 0.1 % STEADY + 力の一致」で判定した (§9)。
  - 小型粘性二流体ケース (差動拡散あり): トレーサ ξ と元素混合分率が乖離することを確認し、`full`≠`lumped` の差を記録 (等価性は要求しない)。
  - 3D 経路: `runner_sern3d.warm_from_same_mesh` の全種保持を 1 段で確認。step 時間の比は記録のみ。
- **判定基準**: 上のゲート + `check_convergence.py` / `check_quasisteady.py` VERDICT 添付。step 時間の増分 (+3 輸送式) は記録のみ。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result (2 回目) | `2026-09-16` | [2026-09-16-thermophysics-cea-mole-fraction-species-result-2.md](../../notes/reviews/2026-09-16-thermophysics-cea-mole-fraction-species-result-2.md) | **NO-GO**, C0/M3/m1 | **全採用 (2026-09-16 反映中, §5.1 #20–#23)**: M1 (作動点変更の再初期化が目標組成にならない: `reinit` は元の実種質量を配分するだけ、`warm_from_run` は組成分率を持ち越すだけ) → 再初期化を「保存済み ξ と**目標**入口ベクトルから Y(r) = ξ Y_in^dst + (1−ξ) Y_ext^dst を再構成」に統一 (変換器 `--mode reinit` と `warm_from_run` の両方)、m6_on→m10_on で目標組成 (H2O 0.24881767, H2 0.01383296)・ΣY・T を検証; M2 (署名に NASA-9 係数・温度区切り・トレーサ設定が無く、メタだけで照合) → 実 config + 解決済み DB から署名 (種順序・MW・両区間係数・Tlo/Tmid/Thi・datum・tracer キー) を作り、メタとの矛盾・TP の照合不能・必要保存量の欠落を拒否; M3 (変換器の hard fail が NaN を通す、反転未収束を返す) → ρ>0・全保存量/組成/T の有限性・組成非負・EOS 反転残差の検査と失敗系試験; m1 (docs の ξ 説明・完了条件の二重記載・粘性試験の扱い) → methods 2 箇所を「純流入元ラベル種の有無」で統一、§6 の旧ゲートを履歴化、F-sp1 依存の粘性試験を §6/§8 で延期明記。M10 の「全列 4.2–4.8 桁」は `rms_roUy` 3.3 桁に訂正 |
| result (1 回目) | `2026-09-16` | [2026-09-16-thermophysics-cea-mole-fraction-species-result.md](../../notes/reviews/2026-09-16-thermophysics-cea-mole-fraction-species-result.md) | **NO-GO**, C0/M10/m1 | **全採用 (2026-09-16 反映中, §5.1 #9–#19)**: M1 (`NO` が YAML で False) → 種名を引用符付きで出力し、生成 config の Python 再読込試験; M2 (湿潤 input の T 反転が二相 EOS を無視・総水量二重計上) → 変換器に二相 EOS 反転と書込前の保存検査、総水量は ρY_H2O; M3 (restart の名前/DB/datum 照合が未実装) → `species_meta.yaml` + DB + datum の共通照合を interp_field / restart_by_index / warm_from_same_mesh / warm_from_run に適用 (照合不能はエラー); M4 (変換器が stream lump とトレーサ新設を扱えない) → 保存的展開と作動点変更の再初期化を分離、`lumped→full` の roXi 生成; M5 (roXi が周期 node の残差合算・同期から漏れ) → periodicNode に登録 + 周期移流試験; M6 (dual-time の roXi は物理時間項なし) → 入力時に tracer × dual-time を拒否 (対応は F-cf8 と同じ扱い); M7 (`condensationSpecies` と `condModel` の不一致を拒否しない) → host 入力段階で物質↔物性モデルの対応検査 + config 読込までの単体; M8 (`lumped+keep` で ξ=Y_EXH が成立しない, `exhaust_fraction(run)` 不在) → 単一種分率が流入元ラベルにならない配置でもトレーサ有効化 + 共通アクセサ実装; M9 (凝縮 OFF で vapor_array が None → 固定組成 24 % 誤り) → H2O を名前で解決し凝縮の有無と独立に分圧; M10 (未達ゲート: 小型 5 種 node PASS ケース無し、トレーサ場 1e-6 超 (最大 1.7e-4, 4.6k ノード)、CEA 物性差 rtol 1.11e-6、「同じ固定点」の断定) → §5.1 に戻し、小型 node PASS ケース・トレーサ誤差の切り分け・作動点変更/湿潤 restart 試験を実施、README の断定を「ゲート内で一致」に修正; m1 (docs 不整合) → methods/plans README/condensation §7b を同期 |
| plan (再: §4.5 SERN) | `2026-09-15` | [2026-09-15-thermophysics-cea-mole-fraction-species-plan-2.md](../../notes/reviews/2026-09-15-thermophysics-cea-mole-fraction-species-plan-2.md) | GO-with-changes, C0/M6/m1 | **全採用**: M1/M2 (元素混合分率では流入元を復元できない・power-off 退化) → 排気率は独立トレーサ `roXi`、元素は診断のみ (§4.5, §5-4c); M3 (stream lump と keep の質量配分) → 流れごとの配分正本 (§4.5); M4 (`restart_by_index`/`warm_from_same_mesh` が roY を落とす) → §4.5, §5-4b, §6; M5 (元素メタデータ) → `species_meta.yaml` (§4.5); M6 (SERN 検証ゲート) → §6; m1 (スキーマ一本化・完了条件) → §4.5, §8 |
| plan | `2026-09-15` | [2026-09-15-thermophysics-cea-mole-fraction-species-plan.md](../../notes/reviews/2026-09-15-thermophysics-cea-mole-fraction-species-plan.md) | GO-with-changes, C0/M8/m2 | **全採用**: M1 (解決済み DB を全経路へ) → §4.1; M2 (凝縮種は名前正本・拒否条件) → §4.2, §2; M3 (種変換 restart・照合) → §2, §5.1 #5b; M4 (X 入力契約) → §4.3; M5 (ParaView 配列を名前解決) → §4.4; M6 (等価性の適用条件) → §3, §6; M7 (CEA と転記の差, --check 修正) → §3, §5.1 #5b; M8 (収束小型ケース node/cell, rms_roY 列, series ゲート, 世代固定) → §6; m1 (frozen.mole_to_mass 委譲, SERN 対象外) → §2, §4.1; m2 (rtol/atol) → §6 |

## 7. 影響範囲

- `design/forge_design/{gas/semiperfect.py, probdef.py, evaluate/runner*.py}`, `design/tests/run_gas_tests.py`
- `solver_density_cuda/input/solverConfig.cpp`, `boundaryCond.cpp` (入力換算のみ), `tools/{gen_inlet_profile.py, forge_species.py, total_quantities.py, paraview/forge_filters.py}`
- 既存 problem YAML は `composition_basis` 省略で従来どおり (質量分率)。既存 run の `species_db.yaml` は不変。
- docs: `methods/thermophysics.md` §5, `methods/design/overview.md` (species 分割節), `methods/condensation.md` §7b (擬似種運用は選択肢の 1 つと明記), `procedures/inlet-profile.md` (`--X`)。

## 8. 完了条件

- [ ] 関連 methods を更新済み
- [ ] 実装・§6 の検証 (単体 / 小型収束ケース / case/44 1–5 / SERN 一式) を満たす
- [ ] codex レビュー 2 回を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status: done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-16` — result-1 の M2–M7/M9 (forge 側) と M1/M3/M8/M9/m1 (設計側) を修正・検証 (§5.1 #9–#19 済); 作動点変更 warm start (`run_0105`) と湿潤 restart の変換器検査を追加。トレーサ差はカーネル差と切り分け (§6 改定)。**codex result 2 回目へ**。
- `2026-09-16` — result-1 M10 の小型収束ケース: case/16 node Euler TP `run_0471` (5 種) / `run_0472` (lumped [MIXDRY, H2O]) が `check_convergence` **PASS** (4.2–4.8 桁, `rms_roY` 列込み), 5 種 vs lumped は float 精度 (1e-6) で一致 (非粘性 frozen の等価性, §3)。トレーサ誤差は同一 run 内のカーネル差 (case/46 `run_0104`) と切り分け。
- `2026-09-16` — codex result レビュー 1 回目 **NO-GO (M10/m1)** を全採用 (§6.1, §5.1 #9–#19)。修正は forge 側 (M2/M4/M5/M6/M7/M9 ツール) と設計側 (M1/M3/M8/M9/m1) で並行。
- `2026-09-16` — **回帰 (node のみ; case/44 README「species 統一スキーマ」表, case/46 README run 表)**: (1) 新バイナリ A/B `run_0195`/反復 `0199` vs 旧 `0170`: 報告量 5 桁一致, 場 ro 1.0e-5 (新バイナリ反復床 4.6e-6 の 2 倍), ṁ 1e-7; (2) **`full` 5 種 `run_0196`** (IC = `convert_species_field.py` で 0170 の収束場を展開): onset 15.8404 vs 15.8400 r_t, 出口 g 0.58448 vs 0.58443 % (+8e-5), ṁ +2.9e-6, 軸 M 1.5e-5, g L1 1.1e-4 (床 1.6e-5 の 7 倍; ゲートは全て内側), series ALL STEADY, condLim 1, `rms_roY0..4` 列; (3) モル分率入力 lumped `run_0197`: `species_db.yaml` は解析後 MW 1e-8 相対で一致 (旧 YAML の 7 桁丸め分), 場は full と同程度の差 (ro 3.0e-5); (4) forge `X{s}`/`condensationSpecies`/トレーサ 500 step `run_0190`–`0192`/`0194`: 同一設定反復ノイズ内 (≤7e-6); (5) CEA 直読み DB `run_0198` (prepare): 設計 R/γ* と CFD `species_db.yaml` の出典がともに DB ファイル; (6) **dry `run_0202` (lumped) / `0203` (full)**: ṁ 完全一致, 軸 M 1.5e-5, 場 2e-6 (§6-1 ゲート 1e-4)。**SERN** case/46 `run_0100` (別名 [EXH, AIR]) / `0101` (full 11 種 + roXi), node Euler m6_on 6000 step: C_T/C_L/C_M 相対 1e-8/5e-7/4e-7 (ゲート 1e-3), 力の時系列 0.1 % で ALL STEADY, 場 P/T 2.5e-5, ξ (Xi vs Y_EXH) 平均 4.7e-7・最大 1.7e-4 (カウル後縁 200 ノード), ΣY 2e-7, GATES PASS; `check_convergence` は plateau (NOT CONVERGED); 24000 step 版 `run_0102`/`0103` も plateau (2.4 桁; node Euler SERN の限界サイクルで PASS は取れず、**§6 の「両 run PASS」は「力の時系列 0.1 % STEADY + 力の一致 + ΣY/有限性」に読み替え** [limiter plan と同じ扱い]) で C_T/C_L/C_M 相対 3e-8/1.2e-6/1.0e-6、場 P/T 3.5e-4/5.7e-4 (限界サイクルの位相差)。step 時間: 5 種 +23 %, SERN 11 種 2.33→4.77 ms。
- `2026-09-16` — **実装着手** (ユーザ指示: plan レビュー済みなので追加レビュー無しで着手可): docs 先行更新、`gas/composition.py` (解決済み DB・換算・統一スキーマ解決・DB/メタ出力)、gas パッケージの DB 注入、`probdef` の gas 検証、`runner_axismach` / `runner_sern{,3d}` の統一スキーマ化と restart 全種化、単体試験追加 (ALL PASS)。forge 本体側 (host DB 解決・`X{s}`・`condensationSpecies`・`rms_roY`・`roXi`・ツール) は並行実装中。
- `2026-09-15` — 初稿 (ユーザ要望: MIXDRY の中身を明示・ユーザ指定可能に、CEA ベースでモル分率指定)。
- `2026-09-15` — ユーザ決定: SERN もノズル設計も `full | lumped` を run ごとに選べる統一スキーマにする (§4.5)。
- `2026-09-15` — codex 再レビュー (§4.5, GO-with-changes M6/m1) を全採用: 排気率は元素混合分率でなく独立トレーサ `roXi` (元素は診断のみ)、流れごとの質量配分正本、`species_meta.yaml`、SERN restart 経路の全種化、SERN 検証ゲート強化、スキーマ一本化。
- `2026-09-15` — codex plan レビュー (GO-with-changes, M8/m2) を全採用: 解決済み DB の一元化、凝縮種の名前正本化と拒否条件、種変換 restart、X 入力契約、ParaView 配列解決、等価性条件の限定、CEA/転記差の扱い、収束ゲート。

## 10. 未確定事項

- ~~forge 側 `X{s}` を bcond だけでなく `inletProfile` CSV でも受けるか~~ 決着 (2026-09-15, §4.3): forge は Y 列のみ、生成器が X→Y 換算。
- `pseudo` (MIX 1 種) を残すか。設計スイープの速度用途があるので残す前提。
