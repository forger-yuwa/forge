# 段の区間判定を起動記録に結び付ける (実効設定・バイナリ id・起動順)

## メタ

- **area**: `その他` (tooling: `solver_density_cuda/tools/stage_manifest.py`、`main.cpp` の起動記録)
- **status**: `draft` (§4・§6 確定 2026-09-27、codex plan 段が次)
- **related_docs**:
  - [`procedures/solver-settings.md`](../../procedures/solver-settings.md) (`mesh.scalarGradient` 節、`slauWallNormalChi` の区間の扱い)
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) §1.0a
  - [`AGENTS.md`](../../AGENTS.md) 「収束確認」の段階起動 run の区間判定
- **related_plans**:
  - [`gradient-scalar-lsq-unification.md`](gradient-scalar-lsq-unification.md) §4.4・§5.1 #2g・#6 (本 plan の完了が `scalarGradient` 既定化の前提。S4 の試験を本 plan で行う)
  - `convection-slau-wall-normal-chi-default` (同じ欠陥が `slauWallNormalChi` の既定化にも効いている、codex M7)
- **created**: `2026-09-27`
- **owner**: `sano`

## 1. 目的

`stage_manifest.py --segment` は「数値設定 (方程式・BC・離散化) が同一の区間」だけを連結して収束判定する。
ところが既定が変わる設定 (`slauWallNormalChi`、次は `scalarGradient`) について、現状は
**同じ YAML を旧バイナリ → 新バイナリで続けて起動すると、両方の段に最後の起動の実効値が付き、別の作用素の段が 1 区間に連結される**
(codex 2026-09-26 gradient-scalar-lsq plan レビュー M7 が実証)。既定化の前にこれを直し、段ごとの実効設定とバイナリで区間を切る。

## 2. スコープ

- **やる**: 起動記録 (`forge_launches.jsonl`) へのバイナリ id の追加、段と起動の対応付け (設定ハッシュ + 起動順)、実効 `scalarGradient`・`slauWallNormalChi`・バイナリ id を hard キーに、対応の付かない段 (legacy / inferred) の扱い、回帰試験 (gradient plan の S4 を含む)、設計 DB (`runner_sern.py`) の行に実効 `scalarGradient` を持たせる。
- **やらない**: `scalarGradient` の既定値の変更そのもの (gradient plan #6)。既存 run の manifest の書き換え (読み取り時の正規化で扱う)。

## 3. 関連 docs と前提 (2026-09-27 調査)

- 起動記録: `main.cpp` の `appendLaunchRecord` が起動ごとに 1 行 (`time`、`cfg_fnv` = solverConfig.yaml の FNV-1a 64、`bcond_fnv`、`exe_size`、`exe_mtime`、`slauWallNormalChi` と由来、`scalarGradient` と由来)。**バイナリの内容を表す id (commit や sha) が無い**。`run_case.sh` は `RUN_PROVENANCE.txt` に `forge_sha256` を書くが、起動記録とは結び付いていない。
- `stage_manifest.py`: `load_launches` が **`cfg_fnv` → 最後の起動の chi 値** の辞書を作る (`:152-163`)。同じ `cfg_fnv` の段は全部その 1 値で上書きされる (`_normalize`、`:279-295`)。`scalarGradient` は読んでいない (2026-09-26 に YAML の hard キーとしてだけ暫定追加、省略と明示 gg は別キー)。
- 由来 (`chi_source`): `explicit` / `inferred` (YAML から推定) / `legacy` (旧 manifest でキー無し = 当時の既定) / `launch` (起動記録で確定)。推定の段と確定の段は値が同じでも連結しない (`segments`、`:297-310`)。
- 設計 DB: `design/forge_design/evaluate/runner_sern.py` の `FLAG_POLICY = "2026-09-26"` と `_last_launch_chi` (最後の起動の chi)。

## 4. 設計方針 (2026-09-27 codex (diagnose) `notes/reviews/2026-09-27-stage-manifest-launch-binding-design-diagnose.md` を全件採用して改訂)

**原則**: 設定が同じことを起動が同じことと扱わない (codex 第 1 仮説: 同じ `cfg_fnv` の起動実効値 [0, 1] が段では [1, 1]・1 区間になることを再現。設定ハッシュ + 起動順の単調対応も、中間の起動記録が欠けると s2→L3 にずれる)。

1. **起動 ID を発行し、残差履歴と直接結ぶ**: forge は起動ごとに一意な `launch_id` (起動時刻 + PID + 乱数) を作り、(a) 起動記録 `forge_launches.jsonl`、(b) `forge_run.log` の先頭、(c) 残差履歴の横に置く `residual_history.launch.jsonl` (その起動が履歴に書いた step 範囲つき。1 つの履歴に複数起動が追記されたら行を増やす) に書く。段を記録する側 (`StageManifest.add`) は、段の履歴ファイルの横の `.launch.jsonl` から `launch_id` の列を読んで段に保存する。**設定ハッシュ (`cfg_fnv`・`bcond_fnv`) は照合用** (manifest にも `bcond_fnv` を保存)、対応付けには使わない。
2. **バイナリ id**: 正本は実行ファイルの **SHA-256 全 64 桁** (`/proc/self/exe` を起動時に 1 回読む。読み取り時間も起動記録に書く)。git rev + dirty は補助 (CMake 埋め込み)。取得に失敗したら `exe_sha256: null` とし、**null 同士を同一とみなさない** (その段は区間判定不能)。
3. **hard キー**: 段に対応した起動の実効 `scalarGradient`・`slauWallNormalChi`・`exe_sha256` を key に入れる。**SHA が違えば実効値が同じでも別区間** (保守側)。
4. **2 種類の由来を分ける**: 記録の根拠 `record_source` (`launch` = 起動 ID で確定 / `inferred` = YAML から推定 / `legacy` / `unknown`) と、設定値の決まり方 `value_origin` (`default` / `explicit` / `auto`、起動記録の `*_source`)。`launch` で後者を潰さない。
5. **legacy と unknown**: 起動 ID の無い段で、対応する旧バイナリ・設定原本で裏付けられる場合だけ `legacy` (例: 旧 manifest のキー無し chi は当時の既定 0)。**`scalarGradient` は記録が無ければ `unknown`** (記録していないことは gg で走った証拠ではない)。**`unknown` の段をまたぐ自動連結は止める** (`--segment` は判定不能を返し、人が区間を明示する)。
6. **設計 DB**: `runner_sern.py` の行は、**評価した成果物 (res を書いた起動) の `launch_id`・`exe_sha256`・実効値**を持つ (最後に成功した起動ではない)。`FLAG_POLICY` の選別はこれで行う。

## 5. 実装ステップ

1. §4 の確定 (上位判断) と codex plan 段レビュー。
2. `CMakeLists.txt` と `main.cpp` (起動記録)、`stage_manifest.py` (`load_launches` を起動の列に、`_normalize` を単調対応に)、回帰試験。
3. `runner_sern.py` の列追加。
4. 検証 (§6) と codex result 段レビュー。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 (**完了 2026-09-27**、判断: codex (diagnose) `notes/reviews/2026-09-27-stage-manifest-launch-binding-design-diagnose.md`・全件採用: 起動 ID と履歴の直接対応、unknown で連結停止、DB は成果物の起動、完全一致検査) | §4・§6 の確定 | 上位判断 | F |
| 2 | codex plan 段 | `codex_review.py --stage plan` | F |
| 3 | 実装 | §5 の 2・3 | O |
| 4 | 検証 | §6 の T1–T7 (AWS で実バイナリ 2 種) | O (結論 F) |
| 5 | codex result 段 | | F |

## 6. 検証 (2026-09-27 codex (diagnose) 採用で改訂。**区間数だけでなく、各段の対応 `launch_id`・実効値・`record_source`・`value_origin` を完全一致で検査**)

- **T1 (gradient plan S4 の中核)**: 同じ YAML (`scalarGradient` 省略) を、既定 gg のバイナリ → 既定だけ lsq に変えた試験ビルド (生産用と分離、差分・SHA・起動記録を保存) で 2 段起動 → 2 区間、段 1 = gg/default、段 2 = lsq/default、各段の `launch_id` が実際の起動と一致。
- **T2**: 同じバイナリで明示 gg → 明示 lsq → 2 区間。明示 lsq → 明示 lsq → 1 区間。
- **T3 (実効値と由来)**: 旧バイナリ省略 / 新バイナリ省略 / 明示 gg / 明示 lsq / **cell で明示 lsq (実効 gg)** の 5 通りで実効値・`record_source`・`value_origin` が期待と完全一致。
- **T4 (欠落耐性、誤対応 0 件が合格)**: 同一設定の 3 段・3 起動で、(a) 中間の起動記録 L2 を欠落させると s1→L1・s2→unknown・s3→L3 で、s2 をまたぐ自動連結を拒否、(b) 余分な失敗起動 (履歴を書かずに落ちた) が混ざっても対応がずれない、(c) 同じ段の再試行 (1 つの履歴に 2 起動) は起動別の step 範囲で記録され、実効値が違えば区間判定不能。
- **T5 (chi)**: T1 と同じ手順で `slauWallNormalChi` の既定が違う 2 バイナリ → 2 区間、各段の chi 実効値が一致。
- **T6 (SHA)**: 実効値が同じで SHA だけ違う 2 段 → 2 区間。`exe_sha256: null` の段は連結しない。
- **T7 (設計 DB)**: DB 行の `launch_id`・`exe_sha256`・実効値が評価対象の res を書いた起動のもの。**成果物生成後に別設定の失敗起動があっても行の由来が変わらない**。`FLAG_POLICY` の旧行が選別から外れる。
- 既存の回帰試験 (`test_stage_manifest_wall_normal_chi.py`、`test_stage_manifest_scalar_gradient.py`、`test_gate_bad_input.py`) が PASS のまま (chi テストの「同一 YAML に異なる cfg_fnv」の fixture は実態に合わせて直す)。
- 起動時の SHA 読み取り時間を AWS と WSL で記録する (許容の判断材料。閾値は設けない)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/CMakeLists.txt`、`main.cpp` (起動記録に 2 項目追加。数値は変えない)、`tools/stage_manifest.py`、`design/forge_design/evaluate/runner_sern.py`。
- 既存 run の `--segment` 判定: 起動記録があり同じ設定の段が複数ある run では、区間の切れ方が変わりうる (より細かく分かれる側)。

## 8. 完了条件

- [ ] §4・§6 確定、codex plan
- [ ] 実装・T1–T7 PASS
- [ ] codex result
- [ ] `status: done`、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-27` — codex (diagnose) で案を却下・修正 (Major 4/Minor 1、全件採用): 設定ハッシュ + 起動順の対応を起動 ID と履歴の直接対応に、記録無しの scalarGradient は unknown、DB は成果物の起動、試験は対応 ID・実効値・由来の完全一致と欠落耐性 (T4)。

- `2026-09-27` — 起票 (gradient-scalar-lsq-unification #2g、Phase 1 合格を受けて)。§4・§6 は案。
