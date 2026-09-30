# 諮問: case A result レビュー 2 回目の採否と accepted 移動の可否

## 観測事実
- レビュー記録: `notes/reviews/2026-09-30-boundary-cht-conjugate-benchmarks-result-2.md` (GO-with-changes, C0/M1/m2)。commit `64312448`。
- `case/64.conjugate_pipe_wall/ab_levels3/EVAL_CONJ_a2_r16_L3.txt`: A2 r16 の加熱区間の熱流束は 差 1.8751e-02 + U 2.2475e-03 > 許容 0.02 (差単独は許容内)。
- 同 A1 r16 (`EVAL_CONJ_a1_r16_L3.txt`): 加熱区間の温度 差 1.2486e-02 (差単独で許容 0.01 超)。
- `run_0008_a2_r16/CONVERGENCE_CHECK.txt` は step 32300 で NOT CONVERGED (安全停止)。同一格子の Df_scale 5 収束対照は無い。
- 最新の再評価 (評価版 70b155cd) で A2 r64 は Q_up を含め全項目 PASS。plan §1 はそれ以前の「判定不能」を書いたまま。

## 期待値と出典
- 合否規則は plan §4.6/§6 (|差| + U ≤ 許容、U > 許容/3 は判定不能)。U の構成は §4.6 で事後改訂 (延長感度 (d) を U から外し別記)。

## 私の採否案
1. M1 採用: A2 r16 の原因帰属を「示唆する」に弱め、「差が許容を超える」「参照側でなく」を削除 (A1 r16 は差単独で許容超なので従来の表現を維持)。FAIL は維持。追加 run なし。
2. m2 採用: 「偽発火」を削除し、「Df_scale 5 で更新量増大を検知して安全停止、Df_scale 20 の別 run で収束、同一格子の Df 5 収束対照は無い」とする (methods/boundary.md・plan #6)。
3. m3 採用: plan §1 の結論を最新評価で書き直し (r32・r64 の A1/A2 全項目 PASS、事後改訂した比較条件である旨)、README の「登録文そのまま」を「事後改訂 (§4.6) を含む」に、#3c を済に、§6.1 に result 1・2 回目の行、plans/README の説明を更新。
4. 以上の反映後に status done・accepted 移動 (A の限定結果、C は後継 plan)。

## 質問
- 上記の採否に異論はあるか。特に M1 の弱めた文面で閉じてよいか、あるいは A2 r16 の追加 A/B (参照 4 水準等) を要するか。
- accepted 移動の前に他に必要なことはあるか。
