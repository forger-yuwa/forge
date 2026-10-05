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

## 4. 設計方針 (草案 — diagnostician 諮問前)

1. **置き場所と呼び方**: `solver_density_cuda/tools/rerun_conditions.py` (forge の run ディレクトリに対する汎用ツール; axismach 専用にしない)。
   `python3 solver_density_cuda/tools/rerun_conditions.py REF_RUN NEW_RUN [--res res_N.h5] [--Pt P] [--Tt T] [--Y H2O=0.09 | --Y1 0.09] [--k K] [--omega W] [--Ps P] [--lump CO2=..,O2=..,N2=..] [--steps N] [--out-interval N] [--cfl C] [--dry-run]`
2. **複製**: 入力ファイルの許可リスト (上記) だけを複製し、`res_*`・ログ・VERDICT・report は持ち込まない。`NEW_RUN` が既にあれば失敗。
3. **書き換え**: bcond は対象行の `floats` だけを書き換え、他の行・キー・順序・書式は保持する (YAML を読み直して書き出すと flow 形式と精度が崩れる — 正規表現での置換 + 書き換え後に YAML として読み直して値を検証)。
   入口の境界種別が `inlet_Pressure` でなければ停止 (対応外を明示)。`Y` は合計 1 を検査 (H2O を変えたら乾き成分 Y0 = 1 − ΣY_other に自動で合わせる、を既定とするか要判断)。
4. **初期場**: 参照 run の最終 `res_*.h5` (または `--res`) を `restart_field.py` で `NEW_RUN/nozzle.h5` へ。lump を変えたときは restart_field が拒否するので、
   `convert_species_field.py --mode reinit` (or conserve) に切り替えるか停止するかは要判断。**条件の大きな変更に対する初期場のスケーリング** (Pt 比で保存量を一様にスケール、Tt 比で速度・温度をスケール) を入れるかは要判断 (Euler では Pt スケールは厳密な相似、NS では Re が変わる)。
5. **記録**: `NEW_RUN/RERUN_CONDITIONS.json` に参照 run・参照 res・変更前後の値・初期場の方法・ツールの git commit を書き、`prepare_info.json` に `rerun_of`・`ic_from` を追記。
6. **乱流の入口量** (k, ω): Pt・Tt を変えたとき自動で合わせるか (例: 入口の乱流強度・粘性比を保つ) は要判断。既定は「変えない、指定があれば上書き」。

## 5. 実装ステップ

1. 諮問 (diagnostician、ユーザ指定) で §4 の要判断 4 点と §6 を確定。
2. codex plan 段レビュー。
3. 実装 + 単体試験 (`solver_density_cuda/tools/test_rerun_conditions.py`)。
4. 検証 run (§6) と文書同期 (procedures §3a・skill)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4 の要判断と §6 の確定 | diagnostician に諮る (ユーザ指定) | F |
| 2 | codex plan 段レビュー | `codex_review.py --stage plan` | O |
| 3 | 実装と単体試験 | 合格: 試験 FAIL 0、書き換え後の bcond を YAML で読み直して値一致、他の行がバイト一致 | O |
| 4 | 検証 run | §6 | O |
| 5 | 文書同期 | procedures §3a・skill nozzle-design・AGENTS の文書表 | O |

## 6. 検証 (草案)

- **単体**: (a) 条件を何も変えずに作った run が参照 run の入力とバイト一致 (bcond・solverConfig) し、初期場が参照の最終場とビット一致; (b) Pt だけ変えると bcond の当該 float だけが変わる; (c) lump を変えると化学種の扱いが §4.4 の決定どおり; (d) 入口種別が違う bcond で停止。
- **検証 run** (case/45、AWS): 参照 `run_0116_ns_recal_final` (or 延長 run_0117) から (i) 無変更の再実行 → 参照の量 (出口コア M・δ_E) を再現 (forge の再実行揺れの範囲 — 幅は参照の末尾 5 枚の幅で事前に決める); (ii) Pt を例えば 0.8 倍 → 完走 (NaN なし)・量が準定常。合否の数値は諮問で決める。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

新規ツール 1 本と試験。既存ツール・ソルバは変えない。

## 8. 完了条件

- [ ] 諮問・codex レビュー 2 回を記録
- [ ] 単体試験と検証 run が §6 を満たす
- [ ] procedures §3a・skill を同期
- [ ] `status: done`、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-10-06` — 起票 (ユーザ依頼「rerun_conditions.py つくってほしい。diagnostician に諮ってね」)。§4 は草案、要判断 4 点を diagnostician に諮る。
