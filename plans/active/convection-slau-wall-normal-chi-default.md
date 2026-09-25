# node の既定を `slauWallNormalChi: 1` にする (三値 auto)

## メタ

- **area**: `convection`
- **status**: `draft`
- **related_docs**:
  - [`methods/convection/theory.md`](../../methods/convection/theory.md) (SLAU の $\chi$、「既知の限界」「対策」「既定」節)
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) §1.0a
- **related_plans**:
  - [`convection-slau-wall-normal-chi.md`](../accepted/convection-slau-wall-normal-chi.md) (実装と受入)
  - [`convection-slau-wall-normal-chi-usage-rule.md`](../accepted/convection-slau-wall-normal-chi-usage-rule.md) (適用規則。本 plan の決定で §4.1 を置き換える)
  - [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (R5o-chi: 3D 格子収束は委譲先)
  - [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) (設計チェーンの runner)
- **created**: `2026-09-26`
- **owner**: `CFD Dev`

## 1. 目的

**ユーザ決定 (2026-09-26)**: node では `nodeWallDirichlet: 1` と `slauWallNormalChi: 1` を既定にする。理由 (ユーザ): 傾向はほぼ変わらず
(2D 生産 m6_on で 5 列とも生産許容内、case/16・case/48 はノイズ同程度)、壁 CV の排出による発散を防げる。
`nodeWallDirichlet` の既定は既に 1 (`solverConfig.hpp:443`) なので、本 plan で変えるのは `slauWallNormalChi` の既定のみ。
2026-09-25 の `diagnostician` 判断「既定化しない」(usage-rule plan §4.1) をユーザ決定で置き換える。

## 2. スコープ

- **やる**: 三値 auto の導入と解決・起動エコー (§4.1–4.2)、`stage_manifest` の実効値化 (§4.3)、ツールと設計チェーンの追随 (§4.4)、docs (§4.5)、
  ビット不変の確認 (B0) と標準 node ケースの A/B (B1)。
- **やらない**: `SLAU_d` カーネルの変更 (面流束の式は前 plan のまま)。3D の格子収束・固定点 (sern-3d R5o-chi)。カウル衝撃足への影響の評価
  (usage-rule plan で判定保留、衝撃衝突構成で flag 1 を使うときの注意として docs に残す)。SLAU のレジスタ上限 (既存の別問題)。

## 3. 関連 docs と前提

- 前 plan の受入 (試験した条件で回帰許容内) と usage-rule plan §6.2 (2D 生産の帯判定、診断可能性、衝撃足は判定保留)。
- 実装の制約: `solverConfig.cpp:680-693` は `slauWallNormalChi: 1` を cell・`nodeWallDirichlet ≠ 1`・solver が SLAU/SLAU2 以外で**起動エラー**にする。
  → 既定は構成に応じて解決する必要がある。
- `stage_manifest.py` は前 plan #12 で「`"0"` ならキーを落とす (省略 ≡ 0)」にしている。既定を変えると省略の意味が変わる。
- 未検証域: 周期・軸対称の生産規模 (V6 は小規模のみ)、凝縮/二相、衝撃衝突点。

## 4. 設計方針 (2026-09-26 `diagnostician`)

### 4.1 三値と解決

`space.slauWallNormalChi` は **省略 (= auto、内部 −1) / 0 / 1**。全キー読込後に `resolveSlauWallNormalChi()` で
**auto → 1 iff `discretization == "node"` ∧ `nodeWallDirichlet == 1` ∧ `solver ∈ {SLAU, SLAU2}`、それ以外 → 0**。
明示 1 の検証エラー (cpp:680-693) は不変。明示 0 は旧挙動。cell・非 SLAU で auto → 0 は**静かな解決 + エコー** (既定なので警告は出さない)。

### 4.2 起動エコー (常に 1 行)

`'slauWallNormalChi' effective: 1 (auto: node+nodeWallDirichlet+SLAU)` / `0 (auto: cell)` / `0 (auto: solver=ROE)` /
`0 (auto: nodeWallDirichlet=0)` / `1 (explicit)` / `0 (explicit)`。manifest・`RUN_PROVENANCE`・`diag_applicability.py` はこの行を正本にする。

### 4.3 `stage_manifest`

- `YAML_HARD_DEFAULTS` の「`"0"` ならキーを落とす」を**廃止**し、**`space.slauWallNormalChi.effective` を常に `"0"` / `"1"` で書く**。
- 解決順: `forge_run.log` のエコー > YAML 明示値 > YAML 省略なら manifest 側で同じ規則 (`mesh.discretization`・`mesh.nodeWallDirichlet`・`solver` から)
  で推定し `inferred: true` を付ける。
- **旧 manifest (キー無し) は `"0"` と解釈** (当時の既定)。
- **同一 YAML でも、バイナリ更新をまたいで継続した run はエコーが違えば別区間**。

### 4.4 ツールと設計チェーン

- `check_solver_config.py`: node + SLAU で省略なら INFO「既定が 1 (2026-09-26)。旧結果の再現は 0 を明記」。
- `migrate_solver_config.py`: **書き換えない** (黙って 0 を足すと旧/新の区別が消える)。
- `design/forge_design/evaluate/runner*.py` が生成する config には**明示 `slauWallNormalChi: 1`** を書き、問題 YAML `mesh.slau_wall_normal_chi`
  (既定 1) で 0 に落とせるようにする (設計 DB の旧エントリと区別するため)。

### 4.5 docs

- `recommended-settings.md` §1.0a: 「既定 1 (2026-09-26 ユーザ決定)。旧挙動の再現は 0 を明記。診断 3 条件は『0 に落とす/落とさない』の判断材料として残す」。
  §9 旧設定表に「省略 = 0 (〜2026-09-25)」を追加。
- `methods/convection/theory.md` の「対策」「既定」節 (本 plan 起票時に更新済み、実装完了で「実装完了までは既定 0」の注記を外す)。
- `solverConfig.hpp:515` のコメント。前 plan 2 本 (accepted) に「ユーザ決定で既定化 (本 plan)」の 1 行。

### 4.6 未検証域の扱い

周期・軸対称・凝縮/二相・化学種は auto に**含める** ($\chi_n$ は面局所で体積ソースを持たず、mask は `wall_flag`。V6 小規模 PASS)。
除外は cell・非 SLAU・`nodeWallDirichlet ≠ 1` のみ (物理的に対象外)。未検証域は §6 の A/B を通すまで plan を `in_progress` に留める。
**不合格時の処置は §6 に先に固定** (その域を auto の除外条件に足し、明示 1 に戻す)。

## 5. 実装ステップ

1. codex plan 段。
2. `solverConfig.{hpp,cpp}` の三値化・解決・エコー。
3. `stage_manifest.py` の実効値化と試験の書き換え。
4. `check_solver_config.py` の INFO、runner の明示キーと問題 YAML キー。
5. B0 (ローカル、1 step) → B1 (case/05 ローカル、他は AWS)。
6. docs、codex result、accepted。

### 5.1 残作業 (優先順)

**計算資源**: B0 と case/05 はローカル (1 step・小規模)。B1 の case/36・case/44・case/46 は AWS (ユーザ指示: ローカルで大きな計算をかけない。
AWS は他セッションと共有なので起動前に `aws_instance.sh status`/`busy`、手動 stop はしない)。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段 | `codex_review.py --stage plan`。Critical/Major は `diagnostician` に諮る | F |
| 2 | 三値化 + 解決 + エコー | `solver_density_cuda/input/solverConfig.{hpp,cpp}`。合格: 単体でエコー 6 パターン + 明示 1 の 3 エラー不変、ビルド成功 | O |
| 3 | manifest 実効値 + 推定 + 旧互換 | `solver_density_cuda/tools/stage_manifest.py`、`test_stage_manifest_wall_normal_chi.py`。合格: node+SLAU で「省略 ≡ 1 ≠ 0」、cell・非 SLAU で「省略 ≡ 0」、旧 manifest (キー無し) ≡ 0、エコーが違えば別区間 の各試験 PASS、`test_gate_bad_input.py` PASS | O |
| 4 | `check_solver_config` INFO、runner の明示キー + 問題 YAML キー | `solver_density_cuda/tools/check_solver_config.py`、`design/forge_design/evaluate/runner_sern*.py`。合格: 生成 config に 1 が出る / 問題 YAML で 0 に落ちる | O |
| 5 | **B0 ビット不変** | 面流束 `FORGE_DUMP_MASSFLUX` 1 step (前 plan V5 の枠)。§6 B0 (a)–(e) | O |
| 6 | **B1 標準 node ケースの A/B** | §6 B1 の表。各 case README の run 表に追記 | O (結論 F) |
| 7 | docs | §4.5 の 4 ファイル。合格: `check_plans.py` PASS | O |
| 8 | codex result → accepted | | F |

## 6. 検証

**判定基準はすべて測る前に固定する。許容は既存のものを再利用し、新設しない。**

### B0 ビット不変 (面流束、1 step)

| # | 比較 | 合格 |
| --- | --- | --- |
| (a) | node + SLAU で省略 (auto → 1) vs 明示 1 | `massflux` ビット同一 |
| (b) | 明示 0 vs 実装前 commit のバイナリ | ビット同一 |
| (c) | cell で省略 vs 実装前 | ビット同一 |
| (d) | 非 SLAU (case/09 KEEP、周期) で省略 vs 実装前 | ビット同一 |
| (e) | 明示 1 の検証エラー 3 種 (cell / `nodeWallDirichlet: 0` / 非 SLAU) | 従来どおりエラー |

### B1 標準 node ケースの A/B (新既定 = 省略 vs 明示 0、同一収束場から分岐、同 step 数、両側 `check_convergence` を併記)

| ケース | 量 | 許容 (出典) | 準定常 |
| --- | --- | --- | --- |
| case/05 sod (SLAU, `convMethod: 2`) | 厳密 Riemann 解との L1 誤差 | limiter plan の既存ゲート値そのまま (悪化なし) | — (非定常、同 step) |
| case/09 TGV (KEEP、周期) | B0 (d) でビット同一なら run 不要 | — | — |
| case/36 node SST 擬似衝撃波 | `check_quasisteady --quantity shock,asym`、壁 $p/p_0$ | 衝撃位置差 ≤ 局所格子 1 つ (前 plan V3 の位置規則)、壁 $p/p_0$ の差区間 ≤ 0.5 % (前 plan case/16 V3) | 各 run STEADY (`--drift 0.001 --osc 0.0025`) |
| case/44 軸対称 + TP + 凝縮 | 軸 M ゲート (case の生産ゲート)、onset 位置、壁 $p/p_0$ | case/44 の既存ゲート値、onset 差 ≤ 局所格子 1 つ、壁 $p/p_0$ 0.5 % | 同上 |
| case/46 2D 3 作動点 (設計チェーン) | runner の GATES、5 列 (`C_T`, `C_T_with_shear`, `C_L`, `C_L_with_shear`, `C_M`) | GATES PASS、差区間が R5n 帯内 (0.002 / 0.002 / 0.05)。m6_on は usage-rule plan §6.2 の結果を再利用 | usage-rule plan §4.3 の閾値 |
| case/48・case/16 | 前 plan V3 の結果を引用 (再 run なし) | — | — |
| 収束悪化 | 全 A/B で flag 1 側の `rms_*` 末尾平均 ≤ 2× flag 0 (usage-rule plan (c)) | — | — |

- **不合格時の処置 (先に固定)**: 量が帯外 → その域 (周期 / 軸対称 / 凝縮) を auto の除外条件に足し、明示 1 に戻す。
  判定不能 (DRIFTING) → 1 回だけ倍に延長、それでも不能なら「未検証のまま auto から除外」。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

### 6.2 結果

(未実施)

## 7. 影響範囲

- `solver_density_cuda/input/solverConfig.{hpp,cpp}` (設定の解決のみ。カーネルは不変)。
- **既存の node + SLAU の config で `slauWallNormalChi` を省略しているもの**は、新バイナリで挙動が変わる (flag 1 になる)。旧結果の再現には `0` を明記する。
- `stage_manifest.py`・`check_solver_config.py`・設計チェーンの runner。
- docs: `methods/convection/theory.md`、`procedures/recommended-settings.md`、前 plan 2 本。

## 8. 完了条件

- [ ] `methods/convection/theory.md` の既定節を実装完了に合わせて更新 (「実装完了までは既定 0」の注記を外す)
- [ ] 実装・検証完了 (§6 の B0・B1)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — 初稿。ユーザ決定 (node の既定を `slauWallNormalChi: 1` に) を受けて起票。§4/§6 は `diagnostician` の設計。
