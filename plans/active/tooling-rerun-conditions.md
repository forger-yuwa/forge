# rerun_conditions: 形状を固定したまま入口条件・背圧・組成だけ変えて回す run を作る

## メタ

- **area**: `tooling`
- **status**: `draft`
- **related_docs**:
  - [`procedures/nozzle-design-workflow.md`](../../procedures/nozzle-design-workflow.md) §3a (本ツールが置き換える手作業の手順)
  - [`.claude/skills/nozzle-design/SKILL.md`](../../.claude/skills/nozzle-design/SKILL.md)
- **related_plans**:
  - [`tooling-nozzle-cfd-pinned-initial-line.md`](tooling-nozzle-cfd-pinned-initial-line.md) (case/45 の生産 run — 本ツールの最初の利用先)
  - [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (将来の campaign 層。本ツールは run 単位の小道具で、それとは独立)
- **created**: `2026-10-06`
- **owner**: `Claude (主セッション) / ユーザ`

## 1. 目的

作成済みのノズル (既存 run の `nozzle.h5` に物理壁の格子が入っている) で、**形状を変えずに**入口条件 (Pt・Tt・H2O 分率・乱流量)・
背圧・乾き成分の組成だけを変えた計算を回したい (ユーザ 2026-10-06「単純に入口条件変えて計算回したい」)。
problem YAML を書き換えて runner/deltastar_loop に通すと `design_chain` が壁を作り直してしまうので、run ディレクトリを直接複製・
編集する手順 (§3a) が要る。手作業だと複製漏れ (`resolved_species_*.yaml` を忘れて forge が起動しない前例)・化学種の互換判断の誤り・
変更記録の欠落が起きやすいので、run の準備 (複製・条件の書き換え・初期場の投入・記録) をツールにする。

## 2. スコープ

- **やる**: run の準備まで (複製 → bcond/solverConfig の書き換え → 初期場 → 記録)。forge の起動・判定はしない (既存の `run_staged_ns`/`run_staged` と判定ツールに任せる)。
- **やらない**: 形状の変更・再設計 (§1・§2 の設計作業)、格子の変更 (cross-mesh)、凝縮 ON/OFF の切り替え (種構成が変わる — `convert_species_field` の既存手順)、複数条件のスイープ管理 (campaign 層)。

## 3. 関連 docs と前提 (観測事実)

- run ディレクトリの入力: `nozzle.h5` (格子 + 初期場、`meshFileName`/`valueFileName`)、`solverConfig.yaml`、`bcondConfig.yaml` (1 行 1 境界の flow 形式、
  例 `inlet: {physID: 1, kind: inlet_Pressure, …, floats: {Y0: 0.9142, Y1: 0.0858, Pt: 5500000.0, Tt: 1600.0, k: 1.0, omega: 18000.0}}`、
  `outlet: {… kind: outlet_statPress, floats: {Ps: 2237.0, Pt: 2237.0, Tt: 300.0}}`)、`species_meta.yaml`・`resolved_species_*.yaml` (化学種の解決記録)、
  `probe.yaml`、`prepare_info.json` (`nozzle_report` が `scale_m`・`mesh.ni`・`Md`・`x_E`・`ic_from`・`dstar_source` を読む)、`wall_*.csv`、`MESH_QUALITY.txt`。
- `solverConfig.yaml` の `physProp.species` は `[{name: MIXDRY, lump: {CO2: …, O2: …, N2: …}, basis: mole}, H2O]`。Pt・Tt に依存する基準値は config に無い (case/45 run_0094 で確認)。
- 同一格子の初期場: `solver_density_cuda/tools/restart_field.py SRC_res.h5 DST_input.h5` は保存量を index コピーし、`forge --resolve-species` で宛先の化学種を解決して
  互換ハッシュが一致しなければ**書き込まずに停止**する (lump の成分比を変えると不一致)。種構成が違うときは `convert_species_field.py --mode conserve|reinit`。
- 一様初期場から始めない (AGENTS.md)。細分格子 NS は本段 CFL 5 で発散する (plan cfd-pinned §5.1 #11a) → 段階起動。

## 4. 設計方針 (diagnostician 2026-10-06 を採用)

1. **置き場所と CLI**: `solver_density_cuda/tools/rerun_conditions.py`。
   `REF_RUN NEW_RUN [--res res_N.h5] [--Pt] [--Tt] [--Y NAME=v | --Y1 v] [--balance NAME] [--k] [--omega] [--Ps] [--Tw | --keep-Tw] [--steps] [--out-interval] [--cfl] [--scale-ic pt|none] [--forge BIN] [--dry-run]`。
   forge は起動しない。`--lump` は受理するが下記 4 の理由を示して**停止** (NEW_RUN を作らない)。
   `--Y` の NAME→index は `solverConfig.yaml` の `physProp.species` から引く。`--steps`/`--out-interval` は `nStepOuter % outStepInterval == 0` を強制、
   `--cfl` は `cfl: X, cfl_pseudo: X` がちょうど 1 回当たることを検査。`--res` 既定は `res_[0-9]*.h5` の最大番号 (参照 config の nStepOuter と違えば警告)。
   SRC が float64 なら `restart_field.py --keep-src-dtype`。
2. **複製は許可リスト**: `nozzle.h5`・`nozzle.xmf`・`bcondConfig.yaml`・`solverConfig.yaml`・`species_meta.yaml`・`resolved_species_*.yaml`・`species_db_external.yaml`・
   `probe.yaml`・`prepare_info.json`・`wall_*.csv`・`target_axis_M.csv`・`delta_r_initial.*`・`MESH_QUALITY.txt`。`EXTENDS.json`・provenance・ログ・VERDICT・series・report は持ち込まない。
   NEW_RUN が既にあれば失敗、途中失敗 (restart_field の失敗を含む) は NEW_RUN ごと削除。
3. **書き換え**: bcond は対象行の floats だけを正規表現で置換 (他の行はバイト一致)、書き換え後に YAML で読み直して値を検証。inlet が `inlet_Pressure` 以外は停止。
   Y は**全種を書き** Σ=1 を 1e−12 で検査 (forge の起動検査は 1e−3 で、入口カーネルが黙って正規化するため頼らない; `speciesDB.cpp:992-1004`・`boundaryCond_d.cu:1064-1069`)。
   残差を吸収する種は 2 種なら自動、3 種以上は `--balance` 必須。`--Ps` は出口の Ps と Pt を同時に書く (Tt 300 は据え置き)。
   等温壁 (`wall_isothermal`) + `--Tt` は `--Tw`/`--keep-Tw` の明示が無ければ停止。`species_meta.yaml` の `streams.inflow.Y_transport`・`Y` を同期 (`X` は実種 MW が記録から取れるときだけ、取れなければキーを落として注記)。
4. **初期場**: `restart_field.py REF_res NEW/nozzle.h5 --dst-run NEW [--forge]` (VERDICT OK 行を要求)。
   **Pt を変えたら既定で 9 保存量を × f (= Pt_new/Pt_ref)** — 保存量は ρ に線形なので T・U・Y・k・ω 不変・P だけ f 倍 (Euler では厳密な相似、NS は Re だけ変わる)。
   検査: 各量 d_new/d_ref が f と相対 1e−6、比 (roe/ro 等) が参照と相対 1e−6。`--scale-ic none` で無効化。f < 0.1 または > 10 は警告。
   **凝縮 block がある参照 run ではスケーリング禁止** (過飽和度が変わる)。**Tt・組成はスケールしない** (semi-perfect で相似不成立・等温壁・ω 床の T 依存で検証できる不変量が無い)。
   **lump 変更は v1 対象外で停止**: restart_field (互換ハッシュ不一致)・convert_species_field conserve (実種ごとの ρY 保存検査で拒否)・reinit (ξ の出所が無い)・forge 自身 (旧ハッシュの属性) の全経路が拒否する (`forge_species.py:1098-1107`、`convert_species_field.py:424-437, 711-719, 395-403, 659-661`)。
5. **記録**: `NEW_RUN/RERUN_CONDITIONS.json` (参照 run/res、変更前後の全値、scale_ic と f、restart_field の VERDICT 行、forge の sha256、ツールの commit、`recommended_stages`)。
   `prepare_info.json` は幾何 (`scale_m`・`mesh`・`Md`・`x_E`) を据え置き、`ic_from` を参照 res で上書き、`rerun_of` を追加 (`initializer.Pt/Tt` は壁を設計した条件として残す)。
6. **入口 k・ω は変えない** (指定時のみ上書き; 設計チェーン自身も Pt・Tt によらず k 1.0・ω 18000 固定、`runner_wt.py:280-285`)。
7. **`recommended_stages`**: Pt のみ・スケール ON・|ln f| ≤ ln 1.25・凝縮なし → `none` (参照 cfl); それ以外 → `full` (細分格子は本段 cfl 1)。閾値 1.25 は未検証 (§6 (ii) で f=0.8 を確かめる)。
8. **Euler 参照**: Tt・組成を変える rerun は、同条件の Euler 参照 rerun (例 run_0114) を対で作る (`nozzle_report` の流量比と δ_E の E 法は Euler 参照を使う)。

## 5. 実装ステップ

1. 諮問 (diagnostician、ユーザ指定) で §4 の要判断 4 点と §6 を確定。
2. codex plan 段レビュー。
3. 実装 + 単体試験 (`solver_density_cuda/tools/test_rerun_conditions.py`)。
4. 検証 run (§6) と文書同期 (procedures §3a・skill)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4 の要判断と §6 の確定 — 判断: 2026-10-06 diagnostician「lump と Tt スケールは v1 から外し、Pt だけの厳密スケーリングを既定 ON、検証は粗格子 run_0094 系列で無変更再実行と Pt 0.8 のスケールあり/なし A/B」 | 完了 | F |
| 2 | codex plan 段レビュー | `codex_review.py --stage plan` | O |
| 3 | 実装と単体試験 | 合格: 試験 FAIL 0、書き換え後の bcond を YAML で読み直して値一致、他の行がバイト一致 | O |
| 4 | 検証 run | §6 (i)(ii): run_0094 res_6000 起点、cfl 5・12000 step | O |
| 6 | (将来、本 plan 外) lump 変更対応 | `convert_species_field` に「lump 名集合が同一で輸送 Y を保ち roe を DB 差分で再構成、ΣY と T だけ検査」モードを足す | F |
| 5 | 文書同期 | procedures §3a・skill nozzle-design・AGENTS の文書表 | O |

## 6. 検証 (diagnostician 2026-10-06、事前登録)

- **単体** (`solver_density_cuda/tools/test_rerun_conditions.py`、fixture は run_0094 の入力): (a) 無変更 → bcond・solverConfig バイト一致 (`--steps/--out-interval` 指定時はその 2 トークンだけ差)、
  restart_field の `VERDICT: OK (9 量を移した、SRC とビット一致)`、9 量が参照 res と `np.array_equal`、属性 `species_hash` 一致; (b) `--Pt 4.4e6` → inlet 行の差分は Pt トークン 1 個・YAML 再読込で一致・他 3 行バイト一致・スケール検査合格;
  (c) `--lump` → 非ゼロ終了・NEW_RUN なし; (d) inlet が `inlet_uniformVelocity` → 非ゼロ終了・NEW_RUN なし; (e) `--Y H2O=0.09` → Y0 = 0.91 (Σ=1 を 1e−12)、3 種 fixture で `--balance` なし → 停止、
  `species_meta.yaml` の Y_transport [0.91, 0.09]; (f) 等温壁 fixture で `--Tt` のみ → 停止; (g) 凝縮 block ありで `--Pt` (スケール既定) → 停止; (h) restart_field を失敗させる → NEW_RUN が残らない。
- **検証 run** (AWS、参照 `case/45.isobutane_m6_d155/run_0094_ns_c2pin_pass2_ext6k/res_6000.h5`、cfl 5・12000 step・1000 ごと = run_0104 と同条件、`run_staged_ns(stages="none")`):
  - (i) 無変更: NaN 0; δ_E(x_F=95.208) が `check_quasisteady --series-csv` STEADY かつ末尾 5 枚平均が run_0104 の 0.725280 に |Δ| ≤ 0.01 %;
    出口コア M 末尾 5 枚平均が 5.999265 に |Δ| ≤ 2e−5 で合格、(2e−5, 1e−4] は保留 (2 本目で再実行揺れを測る)、> 1e−4 は FAIL (ツール欠陥を疑う); check_convergence は plateau 許容、DIVERGED は FAIL。
  - (ii) Pt 4.4e6 (f = 0.8): 腕 A `--scale-ic pt`、腕 B `--scale-ic none`、他は (i) と同一。両腕 NaN 0・δ_E と出口 M が STEADY; A−B の末尾 5 枚平均差 δ_E ≤ 0.02 %・M ≤ 5e−5;
    δ_E(A)/0.725280 ∈ [1.02, 1.10] (乱流 δ ∝ Re^−1/5 → 1.046 に ±2 倍の幅、推測)、出口コア M < 5.99926、`nozzle_report --euler run_0086` の `mdot_ratio_vs_euler` ∈ [0.7985, 0.7995]、壁解像 PASS。
    解釈: B が発散/TRANSIENT-UNSETTLED で A が STEADY → スケール既定 ON 確定; A が発散 → スケール却下 (`none` + `full` を既定に); 両腕 STEADY で一致 → 既定 ON のまま記録; 両腕 STEADY で不一致 > 帯 → 両腕 6000 step を 1 回延長、なお不一致なら保留。
- **方針が誤りと言える観測**: (ii) 腕 A の発散、または A−B が帯を超えたまま延長後も不一致; (i) の M 差 > 1e−4。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-06` | diagnostician (Fable、ユーザ指定) — ブリーフ [`notes/reviews/briefs/2026-10-06-rerun-conditions-plan.md`](../../notes/reviews/briefs/2026-10-06-rerun-conditions-plan.md) | C2/M6/m6 | 全件採用 → §4・§6 確定。C: lump の自動切替却下 (v1 停止)・Tt スケール却下。M: Pt スケール既定 ON (凝縮 run 禁止)・Y 全種書き Σ=1 1e−12・検証参照を run_0094 系列へ・複製許可リスト・記録の拡充・k/ω 自動調整しない |

## 7. 影響範囲

新規ツール 1 本と試験。既存ツール・ソルバは変えない。

## 8. 完了条件

- [ ] 諮問・codex レビュー 2 回を記録
- [ ] 単体試験と検証 run が §6 を満たす
- [ ] procedures §3a・skill を同期
- [ ] `status: done`、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-10-06` — 起票 (ユーザ依頼「rerun_conditions.py つくってほしい。diagnostician に諮ってね」)。§4 は草案、要判断 4 点を diagnostician に諮る。
- `2026-10-06` — diagnostician (ユーザ指定) の判断を全件採用し §4・§6 を確定 (§6.1)。lump 変更と Tt スケーリングは v1 対象外、Pt スケーリングは既定 ON (凝縮 run は禁止)、検証は粗格子 run_0094 系列。
