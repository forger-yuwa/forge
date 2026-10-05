# rerun_conditions: 形状を固定したまま入口条件・背圧・組成だけ変えて回す run を作る

## メタ

- **area**: `tooling`
- **status**: `in_progress`
- **related_docs**:
  - [`procedures/nozzle-design-workflow.md`](../../procedures/nozzle-design-workflow.md) §3a (本ツールが置き換える手作業の手順)
  - [`.claude/skills/nozzle-design/SKILL.md`](../../.claude/skills/nozzle-design/SKILL.md)
  - [`methods/design/overview.md`](../../methods/design/overview.md) (対応入力・初期場変換・拒否条件の節を追加する)
- **related_plans**:
  - [`tooling-nozzle-cfd-pinned-initial-line.md`](tooling-nozzle-cfd-pinned-initial-line.md) (case/45 の生産 run — 本ツールの最初の利用先)
  - [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (将来の campaign 層。本ツールは run 単位の小道具で、それとは独立)
- **created**: `2026-10-06`
- **owner**: `Claude (主セッション) / ユーザ`

## 1. 目的

作成済みのノズル (既存 run の `nozzle.h5` に物理壁の格子が入っている) で、**形状を変えずに**入口条件 (Pt・Tt・H2O 分率・乱流量)・
背圧・輸送種の定義を保った入口分率 (H2O) だけを変えた計算を回したい (乾き成分 lump の組成変更は v1 で拒否) (ユーザ 2026-10-06「単純に入口条件変えて計算回したい」)。
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

## 4. 設計方針 (diagnostician 2026-10-06 → codex plan 段 M1〜m8 を diagnostician が採否して確定)

1. **置き場所と CLI**: `solver_density_cuda/tools/rerun_conditions.py`。
   `REF_RUN NEW_RUN [--res res_N.h5] [--Pt] [--Tt] [--Y NAME=v | --Y1 v] [--balance NAME] [--k] [--omega] [--Ps P | --keep-Ps] [--Tw T | --keep-Tw] [--steps] [--out-interval] [--cfl] [--scale-ic none|pt] [--forge BIN] [--dry-run]`。
   forge は起動しない。`--lump` は受理するが停止 (4 参照)。`--Y` の NAME→index は `solverConfig.yaml` の `physProp.species` から。`--steps`/`--out-interval` は `nStepOuter % outStepInterval == 0` を強制、
   `--cfl` は `cfl: X, cfl_pseudo: X` がちょうど 1 回当たることを検査。`--res` 既定は `res_[0-9]*.h5` の最大番号 (参照 config の nStepOuter と違えば警告)。SRC が float64 なら `restart_field.py --keep-src-dtype`。
   **Pt を変えたら `--Ps` か `--keep-Ps` が必須** (無ければ停止し、参照の出口圧 P_exit_ref [参照 res の出口断面 (最終 x の節点列) の内部節点の静圧 P の中央値。壁・軸の BC 節点は除く。節点列が取れなければ出口 BC の節点と同じ x の節点で代替し、それも取れなければ停止せず `P_exit_ref: null` と警告を記録。`res_outlet_*` の Ps は課した値そのものなので使わない (2026-10-06 改訂)] に対する `Ps/(f·P_exit_ref)` を表示)。
   node の出口は壁列が常に亜音速で Ps を見るため、Pt だけ下げて Ps 据え置きにすると出口列の不安定要因になる (`boundaryCond_d.cu:610-625`、run_0094 の出口壁側 3 節点 M 0.06/0.03/0)。
2. **入力契約 (v1 の対応入力、作成前に検査)**: 単一の `inlet_Pressure` (floats に `Y{s}` 形式、`inletProfile` 無し)、`outlet_statPress`、`wall`/`wall_isothermal`、`axis`、`meshFileName == valueFileName == "nozzle.h5"` (run 内相対)。
   `X{s}` 形式・`inletProfile: 1`・複数 inlet・外部参照・未知の入力依存は NEW_RUN を作る前に拒否。
   **複製は許可リスト**: `nozzle.h5`・`nozzle.xmf`・`bcondConfig.yaml`・`solverConfig.yaml`・`species_meta.yaml`・`resolved_species_*.yaml`・`species_db_external.yaml`・`probe.yaml`・`prepare_info.json`・`wall_*.csv`・`target_axis_M.csv`・`delta_r_initial.*`・`MESH_QUALITY.txt`。
   系譜・ログ・VERDICT・series・report は持ち込まない。NEW_RUN が既にあれば失敗、途中失敗は NEW_RUN ごと削除。
3. **書き換え**: bcond は対象行の floats だけを正規表現で置換 (他の行はバイト一致)、書き換え後に YAML で読み直して検証。Y は**全種を書き** Σ=1 を 1e−12 で検査
   (forge の起動検査は 1e−3 で入口カーネルが黙って正規化する; `speciesDB.cpp:992-1004`・`boundaryCond_d.cu:1064-1069`)。吸収種は 2 種なら自動、3 種以上は `--balance` 必須。
   `--Ps` は出口の Ps と Pt を同時に書く (Tt 300 は据え置き)。等温壁 + `--Tt` は `--Tw`/`--keep-Tw` が無ければ停止。`species_meta.yaml` の `streams.inflow.Y_transport`・`Y` を同期 (`X` は実種 MW が取れるときだけ)。
4. **初期場**: (1) config から**必要保存量集合**を決める = `forge_species.required_conserved` (ro/roU/roe/roY*/roXi) + SST → roK・roOmega + 遷移 → roGamma・roReth + 凝縮 → rog_*・roQ*_* + トレーサ → roXi。
   SRC/DST の存在・shape・有限・ρ>0・0 ≤ roY/ro ≤ 1+1e−6・|ΣroY − ro| ≤ 1e−6·ro を restart_field の前に検査し、DST `/VALUE` にこの集合と `wall_dist` 以外があれば拒否 (codex M1)。
   (2) `restart_field.py REF_res NEW/nozzle.h5 --dst-run NEW [--forge]` (VERDICT OK 行必須)。
   (3) **`--scale-ic pt` のときだけ** (明示 opt-in) 必要保存量を全部 × f (= Pt_new/Pt_ref)。これは **T・U・Y・k・ω を保つ初期場変換**であって境界値問題の相似ではない (Ps・壁温は別の BC)。
   検査 `allclose(d_new, f·d_ref, rtol=1e−6, atol=0)` でゼロはゼロのまま (比で検査しない — ゼロ成分で NaN)。凝縮 block・Tt 変更・組成変更が同時にあればスケール不可 (指定されても停止)。
   (4) lump の組成変更は v1 で停止: restart_field (互換ハッシュ不一致)・convert_species_field conserve (実種ごとの ρY 保存検査)・reinit (ξ の出所なし)・forge (旧ハッシュ属性) の全経路が拒否する。
5. **記録**: `NEW_RUN/RERUN_CONDITIONS.json` (参照 run/res、変更前後の全値、必要保存量集合、scale_ic と f・スケール検査の結果、P_exit_ref と `Ps/(f·P_exit_ref)`、restart_field の VERDICT 行、forge の sha256、ツールの commit、`recommended_stages`)。
   `prepare_info.json` は幾何を据え置き、`ic_from` を参照 res で上書き、`rerun_of` を追加。
6. **入口 k・ω は変えない** (指定時のみ上書き; 設計チェーン自身も固定値、`runner_wt.py:280-285`)。
7. **`recommended_stages`**: 変更なし → `none` (参照 cfl); 条件を 1 つでも変えた → `full` (細分格子は本段 cfl 1)。「スケール IC + none」を Pt のみ変更の推奨に昇格するかは §6 (ii) の結果で別途判断 (§5.1 #8)。
8. **Euler 参照**: Pt のみの変更でも、δ_E 評価には同条件の Euler rerun を対で作る (旧 Pt の Euler 参照だと edge 帯の |傾き| 判定が変わり δ が −4 % 動く合成反例; `deltastar.py:243-250`)。旧条件参照の `mdot_ratio_vs_euler` は診断量として記録のみ。
9. **前提作業: `run_staged_ns` の改修** (codex M4、design 側コード、ソルバ数値は触らない): StageManifest 記録 (`stage_manifest.py` 既存 API)・段ごとの `residual_history_<tag>.csv` 保持・`convMethod: [12]` → 0 の正規表現・段終了ゲート (最終 res の必要保存量が有限・ρ>0、非有限なら次段へ進まず停止)。

## 5. 実装ステップ

1. 諮問 (diagnostician、ユーザ指定) で §4 の要判断 4 点と §6 を確定。
2. codex plan 段レビュー。
3. 実装 + 単体試験 (`solver_density_cuda/tools/test_rerun_conditions.py`)。
4. 検証 run (§6) と文書同期 (procedures §3a・skill)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4 の要判断と §6 の確定 — 判断: 2026-10-06 diagnostician「lump・Tt スケールは v1 から外す」→ 同日 codex plan 段 M1〜m8 を diagnostician が採否し全件採用 (Pt スケールは opt-in、条件変更は full 既定、検証は粗格子でツールのゲート・run_0117 系列で物理受入れ) | 完了 | F |
| 2 | `run_staged_ns` 改修 (§4.9) | `design/forge_design/evaluate/runner_axismach.py`: StageManifest・段ごとの残差履歴・`convMethod: [12]`→0・段終了ゲート。合格: 段ゲートの単体 (NaN 注入で次段へ進まない)、既存の design/tests が通る | O |
| 3 | `rerun_conditions.py` 実装と単体 | §4.1〜4.8、§6 単体 (a)〜(m)。合格: `test_rerun_conditions.py` FAIL 0 | O |
| 4 | `methods/design/overview.md` 追記 | 対応入力・初期場変換 (T・U・Y・k・ω を保つ)・拒否条件 | O |
| 5 | 検証 run (i)(ii)(iv) | §6、AWS、粗格子 run_0094 系列。case README 台帳に run を追記 (主セッション) | O |
| 6 | 文書同期 | `procedures/nozzle-design-workflow.md` §3a・skill `nozzle-design` | O |
| 7 | codex result 段レビュー | `codex_review.py --stage result` | O |
| 8 | 「スケール IC + none」を Pt のみ変更の推奨に昇格するか | §6 (ii) の結果で判断 | F |
| 9 | (本 plan 外) lump 変更対応・物理受入れ (iii) | `convert_species_field` に lump 同名・輸送 Y 保持のモード; (iii) は最初の生産利用で | F |

## 6. 検証 (事前登録; diagnostician 2026-10-06 の改訂版)

**参照の使い分け**: ツールのゲート (再現性・入力契約・full 経路) = 粗格子 `case/45.isobutane_m6_d155/run_0094_ns_c2pin_pass2_ext6k/res_6000.h5` 系列 (壁解像 FAIL 29.4 % は既知として記録)。
物理の受入れ = 壁解像 PASS の `run_0117_ns_recal_final_ext/res_60000.h5` 系列 ((iii)、最初の生産利用で実施し本 plan の完了条件には含めない)。

**準定常の判定規約 (全 run 共通)**: 系列は res_1000..res_N (res_0 は含めない)。`check_quasisteady.py --series-csv <run>/quantities_series.csv --series-cols <列> --tail 0.4 --min-snaps 10` (n=12 → 末尾 5 枚)。
- δ_E(x_F) [`delta_E`]: `--drift 5e-5 --osc 1e-4` (末尾幅 ≤ 0.01 %)、直前窓差 |mean(8000–12000) − mean(3000–7000)| ≤ 0.01 %。
- 出口コア M [`exitM_A`]: `--drift 3e-6 --osc 3e-6` (末尾幅 ≤ 1.8e−5)、直前窓差 ≤ 2e−5。
- 流量 [`mdot`、2π∫ρU_x r dr の中央値]: `--drift 1e-5 --osc 2e-5`。
- STEADY 以外は不合格・延長 1 回 (6000 step) まで、なお不達は保留。check_convergence は NOT CONVERGED (plateau) を記録のみ、DIVERGED は FAIL。`--from-floor` は使わない。

**単体** (`solver_density_cuda/tools/test_rerun_conditions.py`、fixture は run_0094 の入力):
(a) 無変更 → bcond・solverConfig バイト一致 (`--steps/--out-interval` 指定時はその 2 トークンだけ差)、restart_field の VERDICT OK、必要保存量が参照 res と `np.array_equal`、`species_hash` 一致;
(b) `--Pt 4.4e6 --Ps 1789.6` → inlet/outlet 行の差分は当該トークンだけ・YAML 再読込で一致・他の行バイト一致; (c) `--lump` → 非ゼロ終了・NEW_RUN なし; (d) inlet が `inlet_uniformVelocity` → 拒否;
(e) `--Y H2O=0.09` → Y0 = 0.91 (Σ=1 を 1e−12)、3 種 fixture で `--balance` なし → 停止、`species_meta.yaml` の Y_transport [0.91, 0.09]; (f) 等温壁 + `--Tt` のみ → 停止; (g) 凝縮 block + `--scale-ic pt` → 停止;
(h) restart_field を失敗させる → NEW_RUN が残らない; (i′) 必要保存量が SRC に欠ける → 拒否; (j) 3 種 + トレーサ fixture で `--scale-ic pt` → roY2・roXi も f 倍、ゼロ成分はゼロ;
(k) `inletProfile: 1`・X{s} 形式・inlet 2 本・`valueFileName` が別名 → 作成前に拒否; (l) `--Pt` のみ (Ps 指定なし) → 停止し `Ps/(f·P_exit_ref)` を表示; (m) 無変更 → `recommended_stages: none`、`--Tt` → `full`。
`run_staged_ns` 改修の単体: 段ゲートが非有限の res で次段へ進まない (soft 段の res に NaN を注入)。

**検証 run** (AWS):
- (i) 無変更の再実行 (run_0094 → `--steps 12000 --out-interval 1000`、`stages="none"`、cfl 5): NaN 0; δ_E STEADY (規約) かつ末尾 5 枚平均が run_0104 の 0.725280 に |Δ| ≤ 0.01 %;
  出口コア M 末尾平均が 5.999265 に |Δ| ≤ 2e−5 (〜1e−4 は保留、> 1e−4 は FAIL)。Euler 参照は run_0086 (同条件)。
- (ii) Pt 0.8 倍 (4.4e6 Pa、`--Ps 1789.6` = 0.8 × 2237、cfl 5・12000 step・`stages="none"` を実験として明示):
  - 腕 E (Euler 対参照): run_0086 → `--Pt 4.4e6 --Ps 1789.6 --scale-ic pt`、`run_staged(stages="none", cfl 2)` 6000 step。ゲート: 出口コア M 末尾平均が run_0086 に |Δ| ≤ 2e−5、流量が 0.8 × run_0086 に相対 1e−4 (状態変換が Euler 解を保つことの直接検証)。E が不合格なら (ii) 全体を保留し再諮問。
  - 腕 A (`--scale-ic pt`) / 腕 B (`--scale-ic none`)、δ_E は腕 E を参照に抽出。ゲート: 両腕 NaN 0・δ_E と M が STEADY; |A−B| が δ_E ≤ 0.02 %・M ≤ 5e−5。
    壁解像は記録のみ (既知 FAIL; y1+>1 割合・最大値が run_0094 の 29.4 %/13.4 より下がることを記録)。診断量 (合否外): `mdot_ratio_vs_euler` (run_0086 基準、予想 ≈ 0.7995)、δ_E(A)/0.725280 (予想 1.02〜1.10)、M の向き (予想: 低下)。
  - 解釈 (先に固定): B 不安定・A STEADY → 「スケール IC + none」を Pt のみ変更の opt-in 推奨として §4.7 に追記 (既定は full のまま); A 発散 → スケール機能を削除候補に; 両腕 STEADY で一致 → 記録のみ; 不一致 → 延長 1 回 → 保留。
- (iv) `full` 経路の結合試験: run_0094 → `--Tt 1500 --Y H2O=0.10 --keep-Ps`、`run_staged_ns(stages="full")` (soft 3000・mid 3000・本段 cfl 5・6000)。
  ゲート: `stage_manifest.json` に 3 段、各段の `residual_history_<tag>.csv` が有限、段ゲートで停止なし、本段 NaN 0、実効入力 (bcond の Tt/Y1/Y0、`species_meta` の Y_transport) が要求値、本段の δ_E・M が STEADY (値の合否なし)。
  (Tt 1500・H2O 0.10 は経路試験用の任意値で報告に使わない。)
- (iii) 物理受入れ (最初の生産利用で実施、本 plan の完了条件外; 担当 F で条件決定): 参照 run_0117 res_60000、同条件 Euler rerun (run_0114 から) + NS rerun (`stages="full"`、本段 cfl 1・60000 step、Pt 変更なら `--Ps` 必須)。
  ゲート: NaN 0; δ_E・出口コア M・流量が規約で STEADY; 壁解像 PASS; 値は報告 (参照値が無いので合否なし)。
- **方針が誤りと言える観測**: (i) の M 差 > 1e−4 または δ_E 差 > 0.01 %; (iv) で段ゲートが非有限を通した、または実効入力が要求値と違う; (ii) 腕 E の不合格。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-06` | diagnostician (Fable、ユーザ指定) — ブリーフ [`notes/reviews/briefs/2026-10-06-rerun-conditions-plan.md`](../../notes/reviews/briefs/2026-10-06-rerun-conditions-plan.md) | C2/M6/m6 | 全件採用 → §4・§6 確定。C: lump の自動切替却下 (v1 停止)・Tt スケール却下。M: Pt スケール既定 ON (凝縮 run 禁止)・Y 全種書き Σ=1 1e−12・検証参照を run_0094 系列へ・複製許可リスト・記録の拡充・k/ω 自動調整しない |
| plan | `2026-10-06` | [`notes/reviews/2026-10-06-tooling-rerun-conditions-plan.md`](../../notes/reviews/2026-10-06-tooling-rerun-conditions-plan.md) | GO-with-changes, C0/M7/m1 | 採否は diagnostician (同日 2 回目) が判断し全件採用: M1 必要保存量集合とゼロ対応の検査、M2 Pt スケールは opt-in・条件変更は full 既定・Pt 変更時は Ps 指定必須 (前回の既定 ON を撤回)、M3 入力契約、M4 run_staged_ns 改修を前提作業に、M5 同条件 Euler 参照、M6 準定常の閾値を量ごとに事前登録、M7 物理受入れは壁解像 PASS の run_0117 系列に分離、m8 methods と台帳 → §4・§5.1・§6 |

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
- `2026-10-06` — codex plan 段レビュー (GO-with-changes, C0/M7/m1) の採否を diagnostician (同日 2 回目、ユーザ指定) に諮り全件採用。diagnostician は前回の「Pt スケール既定 ON + stages none」を撤回 (出口壁列が Ps を見るので Pt だけ下げると出口不安定要因、合成反例で旧 Pt の Euler 参照は δ を −4 % 動かす)。§4・§5.1・§6 を確定版に置換、status → in_progress。
- `2026-10-06` — §5.1 #2・#3・#4 を実装 (implementer)。#2 `runner_axismach.run_staged_ns`: 段ごとに `stage_manifest.json` (`S1_soft`/`S2_mid`/`R<i>_cfl<c>`/`main`) と `residual_history_<tag>.csv` を残す・前段の 1 次化を `convMethod: [12]\b`→0・段終了ゲート `stage_gate` (最終 res の必要保存量が有限・ρ>0 でなければ次段へ進まず RuntimeError、restart_field へ渡さない)。段の CFL・step 数は不変。#3 `solver_density_cuda/tools/rerun_conditions.py` (§4.1〜4.8) と `test_rerun_conditions.py`: 単体 (a)〜(m)・段ゲート・restart を模擬した作成経路で FAIL 0。ただし (a)(j) の**作成経路** (実 restart_field) は、ローカルに `--resolve-species` を持つ forge が無いため SKIP (理由つき; `--force-species` では通していない) — AWS など新しい forge (FORGE_BIN) で再実行が要る。#4 `methods/design/overview.md`「既存 run の条件変更 (rerun_conditions)」節。
- `2026-10-06` — §4.1 の P_exit_ref の定義を改訂 (implementer の指摘を主セッションが採用: `res_outlet_*` の Ps は課した値 2237 そのもので比較にならない)。参照 res の出口断面の内部節点 (壁・軸の BC 節点を除く) の P の中央値に変更し、`rerun_conditions.exit_pressure_ref` と単体試験を更新 (run_0094 res_6000 で P_exit_ref = 2241.998 Pa、95/97 節点; FAIL 0・SKIP 2)。
- `2026-10-06` — **単体 FAIL 0・SKIP 0 (AWS)**: ローカルでは --resolve-species を持つ forge が無く (a)(j) の作成経路を SKIP したが、AWS で `FORGE_BIN=~/forge-wallfit-bin/.../forge` を渡して全項目合格。**検証 run 投入** (`case/45.isobutane_m6_d155/run_rerun_validation.sh`、量の時系列 `rerun_series.py`): (i) `run_0119_rerun_ctrl`、(ii) `run_0120_rerun_euler_pt08` (E)・`run_0121_rerun_pt08_scale` (A)・`run_0122_rerun_pt08_noscale` (B)、(iv) `run_0123_rerun_fullpath`。§6 からの変更 (投入前): E と (iv) は 6000 step だと準定常判定の 10 枚に届かないので出力間隔を 500 に (判定規約 `--tail 0.4 --min-snaps 10` を守るため; step 数は不変)。
- `2026-10-06` — 検証バッチが Euler 対参照 (run_0086、壁 `kind: slip`) で「対応外の境界種別」により作成前に停止 (forge 未起動; run_0119 は作成済み)。§4.2 の入力契約に `slip` (Euler の滑り壁、壁温なし) を追加し、P_exit_ref の除外にも含めた。単体 (f′) 追加、FAIL 0。バッチは既存の run_0119 を作り直さずに再開。
- `2026-10-06` — 再開した検証バッチが Euler 対参照の作成で再停止: run_0086 の `nozzle.h5` に乱流モデルなしでも変換器が作る roK/roOmega の入れ物があり、§4.4 (1) の「必要保存量と wall_dist 以外は拒否」に当たった。乱流モデルなし (Euler・層流) の run に限り roK/roOmega を「使われない余りの量」として許し警告に記録 (スケールしない; 他の未知の量は従来どおり拒否)。単体 (f″) 追加、FAIL 0。
