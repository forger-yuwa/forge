# codex レビュー: tooling-sern-te-wake-grid (result)

- **plan**: [`plans/accepted/tooling-sern-te-wake-grid.md`](../../plans/accepted/tooling-sern-te-wake-grid.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-12
- **commit**: `857adf6d` (feature/sern-design)
- **codex**: effort `high`, 4.7 min, rc=0
- **判定**: **GO**, 指摘 C0/M0/m0
- **focus**: 2026-10-12 の result 再レビュー 3 回目 (NO-GO、M1: 旧判定の PASS 行が残差検査の版を問われずに学習・Pareto に入る) への対応 R2k が塞いでいるか。生産の標準経路と学習経路で、必須の証拠が欠けた評価が受理・延長・学習に入らない状態になったか。移管した項目の追跡先。status done にして accepted へ移してよいか。
- **extra**: `case/46.sern_design/README.md`, `notes/investigations/2026-10-12-r2c/run_1100_r2c_lsw10_m6_judgment/R2C_SUMMARY.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO**  
R2k は前回の Major を解消しています。旧判定の不完全な評価は学習・Pareto から除外され、完全な入力を現行ゲートで再判定した行だけが採用されました。今回の対象範囲で追加指摘はありません。

**指摘一覧: Critical / Major / Minor ともになし。**

確認した根拠は次のとおりです。

- **実装と学習経路**  
  [`sern_gates.py:65`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:65) の受理検査の版と残差証拠を、[`driver_sern.py:303`](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:303) が一律に要求しています。学習配列 `_XF` と Pareto `summary` の両方が、この判定を通ります。旧判定を要約し直しただけでは現行版が付与されません。

  `run_1100` の生の残差をメモリ内で変更し、他の証拠を保存結果に固定して再現しました。

  | 残差入力 | 旧 PASS 行の現行学習判定 | 現行ゲートで再判定後 |
  |---|---|---|
  | 完全 | 除外 | PASS・採用 |
  | `rms_roY0` 欠落 | 除外 | 判定不能・除外 |
  | 化学種２列欠落 | 除外 | 判定不能・除外 |
  | 先頭／末尾１行だけ | 両方除外 | 両方判定不能・除外 |
  | 途中の step 欠測／重複 | 両方除外 | 両方判定不能・除外 |

  不完全な６種類は延長対象にも入りません。学習配列・Pareto に残ったのは、現行ゲートで再判定した完全な１行だけです。

- **標準経路の実測**  
  `case/46.sern_design/run_1100_r2c_lsw10_m6_judgment/` の保存系列から、１回目 **UNSETTLED**、２回目 **PASS** を再現しました。再実行した判定は、残差が両回とも **NOT CONVERGED (stalled/plateau)**、２回目の力履歴が **OVERALL: ALL STEADY**、床事象が **VERDICT: PASS**。判定区間は最後の session の `(10000, 20000]` です。収束したとは認定しません。[run 索引と成果物](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:409)

- **格子感度と回帰**  
  `case/46.sern_design/run_1084_tewake_g3_L10_m10/` と `case/46.sern_design/run_1086_tewake_g4_L10_cont40k/` の保存履歴から、§6 #2・#3 の PASS と、４力量の平均差が登録許容内であることを再確認しました。保存された品質判定は g3 **PASS**、g4 の入力格子 `run_1082_tewake_g4_L10_m10/` は **SOFT-PASS**、局所変形の判定は両方 **TE_WAKE_LOCAL PASS**。残差は両方 **NOT CONVERGED** の記録です。[比較対象の索引](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:401)  
  `run_sern_mesh3d_tests.py` は **ALL PASS**、対象 plan の lint も **PASS** でした。

- **文書と残作業**  
  R2k の仕様・影響は methods と plan に記録されています。共有ツール #12〜#15、親 plan R15・R7b、2D、ghost 初期化の移管先も確認しました。今回の GO は、移管済み課題の解決や 3D 最適化の再開を認定するものではありません。

通常の R2k 試験全体は、一時ディレクトリ作成が read-only 制約で拒否され、再実行できませんでした。上表は独立したメモリ内試験の結果です。GPU 計算は再実行していません。

**推奨:** 今回のレビュー記録を反映し、`status: done`・索引・参照リンクを同期して `accepted` へ移してください。本レビューではファイルを変更しておらず、**plan 未反映**です。

指摘数: Critical 0 / Major 0 / Minor 0
