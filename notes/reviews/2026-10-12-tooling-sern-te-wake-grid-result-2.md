# codex レビュー: tooling-sern-te-wake-grid (result)

- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-12
- **commit**: `3c640956` (feature/sern-design)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M1/m0
- **focus**: 2026-10-12 の result 再レビュー (NO-GO、M1 化学種の残差列・M2 残差の時系列の完全性) への対応 R2i・R2j が塞いでいるか。生産の標準経路で、必須の証拠が欠けた評価が受理・延長・学習に入らない状態になったか。移管した項目 (共有ツールの plan #12〜#15・親 R15・R7b・2D plan) の追跡先と適用範囲。status done にして accepted へ移してよいか。
- **extra**: `case/46.sern_design/README.md`, `notes/investigations/2026-10-12-r2c/run_1100_r2c_lsw10_m6_judgment/R2C_SUMMARY.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
R2i・R2j は、新規評価での残差欠落の受理・自動延長を防いでいます。  
ただし、旧ゲートで PASS になった不完全な評価が、保存済み台帳から学習に入る経路が残っています。

1. **Major — 残差検査の改訂が、保存済み評価の学習可否に反映されない**

   **根拠:** [`driver_sern.py:260`](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:260) の `_learnable` は、保存された PASS・評価方式・床事象・判定区間の要約を確認しますが、R2i・R2j による残差検査の実施を要求しません。さらに [`driver_sern.py:476`](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:476) は残差の VERDICT 文字列だけを台帳へ渡し、`input`・`series` の検査結果を保存していません。`FLAG_POLICY` も `2026-09-27` のままです。

   **再現:** `run_1100` の生の残差をメモリ内で変更し、他の証拠を保存結果に固定しました。旧実装 `9e5f0ffd` の判定から `_op_summary` で台帳行を構成し、現行 `_learnable` に評価方式・床・判定区間の要件をすべて渡した結果です。

   | 残差の入力 | 現行ゲートで再評価 | 旧 PASS 行を現行の学習判定へ |
   |---|---|---|
   | 完全 | PASS | 採用 |
   | `rms_roY0` 欠落 | `JUDGMENT_UNVERIFIABLE` | **採用** |
   | 化学種２列とも欠落 | `JUDGMENT_UNVERIFIABLE` | **採用** |
   | 先頭／末尾の１行だけ | `JUDGMENT_UNVERIFIABLE` | **採用** |
   | 途中の step 欠測／重複 | `JUDGMENT_UNVERIFIABLE` | **採用** |

   これは元の `run_1100` が不完全だったという指摘ではありません。**旧判定を持つ不完全な評価を、再判定なしで採れる**ことの再現です。学習配列と Pareto の双方がこの関数を使います（[`driver_sern.py:560`](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:560)、同 `:598`）。

   **対案:** 評価方式 ID とは別に受理検査の版と残差入力の検査済み要約を保存し、学習時に要求してください。版・証拠のない旧行は再判定まで除外し、現行ゲートで再判定した行だけ復帰させます。上表を学習・Pareto まで通す回帰試験も必要です。

確認できた範囲では、次は支持できます。

- **R2i・R2j の新規受理・延長経路:** 上記の不完全入力は拒否され、完全な入力では１回目が延長可、２回目が PASS のままです。
- **`case/46.sern_design/run_1100_r2c_lsw10_m6_judgment/`:** 生の残差２本は試験登録の SHA-256 と一致。再実行した `check_convergence` は両回とも **`NOT CONVERGED (stalled/plateau)`**、２回目の力履歴の `check_quasisteady` は **`OVERALL: ALL STEADY`**。保存系列から１回目 UNSETTLED・２回目 PASS を再現しました。[run 索引](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:409)
- **回帰・文書:** `run_sern_mesh3d_tests.py` は **ALL PASS**、対象 plan の lint は **PASS**。共有ツール #12〜#15・親 R15/R7b・2D plan の移管先は存在します。ただし、上記の旧判定の失効は残作業に登録されていません。

**推奨:** この学習経路の修正と回帰確認を完了してから `accepted` へ移してください。今回の指摘は読み取り専用レビューのため **plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 0
