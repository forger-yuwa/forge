# codex レビュー: boundary-cht-conjugate-flat-plate (result)

- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **stage**: `result` (diff base `feature/cht-axisym-graetz`)
- **date**: 2026-10-01
- **commit**: `e904cb94` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m3
- **extra**: `plans/accepted/boundary-cht-conjugate-benchmarks.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

収束・準定常・主判定の主要数値は裏付けられ、n32/n64 の限定結果は支持できます。  
ただし、登録された固体効果の判定と領域感度評価が未完了です。以下を補ってから `accepted` に移してください。

1. **Major — C2 の「固体効果」の判定量が登録と異なる**

   **根拠:** 発注元 [plan:99](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:99) は、評価窓内の「軸方向伝導あり／なしの界面温度差」を検出能力の根拠にしています。しかし [eval_conj.py:307](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:307) は `max|Q_ax|/Q_tot` を熱流束の比較許容 `0.03` と比べています。断面積分した熱量比と局所熱流束の誤差許容は、直接比較できる判定量ではありません。「効果/許容 1.9〜2.0」は登録した検出能力の証明になっていません。

   `case/65.conjugate_flat_plate/run_0027_c2_n32_lim0/` の最終流れ場を固定し、参照解の固体軸方向伝導だけを切り替えて再計算すると、参照 level 1→2 で次を得ました。

   - 窓内 `max|Δθ|`: **0.013205 → 0.013032**
   - 窓内 `max|Δq|/q_mean`: **0.055915 → 0.054517**

   効果自体は確認できますが、現行評価器はこの差と、その不確かさを判定していません。

   **対案:** C2 の採用格子 n32/n64 で、同じ凍結流れ場・評価窓による伝導あり／なしの差を評価し、差の不確かさ込みで登録条件を確認する。`Q_ax/Q_tot` は補助指標として残してください。

2. **Major — C の領域感度が、別掲も残作業登録もなく抜けている**

   **根拠:** 発注元 [plan:107](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:107) は「C の上境界を1.5倍」にする感度評価を要求しています。事後改訂は比較の `U` から外して**別掲する**もので、評価自体の削除ではありません。

   一方、[eval_conj.py:96](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:96) は C の延長を「後継 plan で扱う」として拒否し、[同:238](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:238) では A だけを実行します。今回の後継 plan の事後改訂一覧・残作業表にも、この未実施項目がありません。

   **対案:** 延長部分の凍結流れ場・源項の与え方を明記して上境界延長の参照計算を行い、温度・熱流束への感度を別掲する。同一有限領域に限定した主判定の `U` と混ぜる必要はありません。

3. **Minor — 「格子間差」は実際には各格子での参照解との差**

   **根拠:** [plan:29](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:29) の C1 `0.72→0.23→0.085 %`、C2 `2.23→0.73→0.31 %` は、`EVAL_CONJ.txt` の **forge−参照解**の差です。各参照解がそれぞれの forge 流れ場を使うため、forge 自身の格子間差とは区別が必要です。

   保存 CSV から、隣接する forge 格子の共通節点で直接比較すると、熱流束差は次のように減少しました。規格化は各比較の共通窓における細格子の平均熱流束です。

   - C1、`run_0024→0026→0022`: **0.523 % → 0.111 %**
   - C2、`run_0025→0027→0023`: **1.588 % → 0.350 %**

   **対案:** 現行数値を「参照解との差」に改称し、forge 同士の格子間比較を別に記録する。細分化による減少という結論は維持できます。

4. **Minor — 判定ファイルを上書きする欠陥が実行チェーンに残っている**

   **根拠:** [lim0_chain.sh:61](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/lim0_chain.sh:61) は、ツール自身も保存する `CHT_INTERFACE_VERDICT.txt` に標準出力をリダイレクトし、さらに旧閾値で再実行します。[check_cht_interface.py:59](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:59) の保存処理と競合します。今回の成果物は `lim0_regate.sh` で修復されていますが、チェーンを再利用すると再発します。

   **対案:** チェーン本体を修正し、旧閾値の結果を別名で保存してから登録閾値の判定を正本へ保存する。修復スクリプトと処理を共通化してください。

5. **Minor — 閉鎖時の文書同期と未完了項目の扱いが不足している**

   **根拠:** [plans/README.md:27](/home/sano/work/forge-cht/plans/README.md:27) は依然 `draft`・旧結果のままです。[case README:3](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:3) は「登録文そのまま」と記載し、まだ存在しない `accepted` 側を参照しています。また [plan:221](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:221) の非一様固体温度での丸め分離は未完了です。

   **対案:** 移動時に索引・リンク・事後改訂の記述を同期する。丸め分離と `limiter: 2` の停滞原因について、限定閉鎖の対象外とする理由、継続先または保留状態を残作業表へ明記してください。

照合できた結果も明記します。[run 一覧](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:9) の `run_0022〜0027` は、圧縮残差をメモリ上で展開して公式判定関数を再実行し、全件 **`PASS (converged)`** を再現しました。準定常も保存系列の全列で **`STEADY`**、G-if は登録閾値で全件 **`PASS`**、旧閾値では④だけ **`NOT CONVERGED`**。最終場の `VALUE/*` に NaN/Inf はありません。保存された G-cons は全件 **`VERDICT: PASS`**、3格子の品質記録も **`VERDICT: PASS`** です。

C2 n32 の主評価を書き込みなしで再実行し、**`VERDICT: PASS`** と記載数値を再現しました。C2 n16 の4水準結果 **`FAIL`（2.3536 %＋0.97918 %）** も保存 CSV と整合します。参照解の自己検査は **`VERDICT: PASS`**。指定 diff に CUDA/C++ ソルバ本体・`design`・既定設定の変更はなく、共有準定常ツールの既定動作は旧版との200系列比較で一致しました。FP32・cell・周期 seam の精度検証へは拡張できません。

**推奨は、既存 run を維持し、①登録された固体効果の評価 → ②領域感度の別掲 → ③格子比較・チェーン・文書の修正、の順に補完して、n32/n64 の限定結果として閉じることです。** 現時点でソルバ本体の変更や全面再走を要求する根拠はありません。

ファイル変更なし。指摘・提案は plan 未反映です。

指摘数: Critical 0 / Major 2 / Minor 3
