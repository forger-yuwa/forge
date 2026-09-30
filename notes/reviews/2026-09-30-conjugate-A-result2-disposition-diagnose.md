# codex 諮問 (diagnose): conjugate-A-result2-disposition

- **brief**: [`notes/reviews/briefs/2026-09-30-conjugate-A-result2-disposition.md`](../../notes/reviews/briefs/2026-09-30-conjugate-A-result2-disposition.md)
- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **date**: 2026-09-30
- **commit**: `64312448` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 2.8 min, rc=0
- **結論**: **追加計算を要求せず、上記を文書全体に反映してから、A の「事後改訂条件・FP64・node・流れ場固定」の限定結果として `accepted` へ移すことを推奨する。**
- **extra**: `case/64.conjugate_pipe_wall/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 指摘 | 採否 | 根拠と必要な処置 |
|---|---|---|
| **Major M1：A2 r16 の原因断定** | **採用。追加 run 不要** | [L3 評価:4](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/ab_levels3/EVAL_CONJ_a2_r16_L3.txt:4) は `0.018751 + 0.0022475 = 0.0209985 > 0.02`。**FAIL は正しいが、差単独は許容内**。「追加細分後も不確かさ込みで FAIL。forge の粗格子誤差の寄与を示唆するが、参照側の寄与は排除していない」とする。 |
| **Minor m2：安全停止の偽発火扱い** | **採用** | [停止条件:1044](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1044) は更新量の連続増大と累積増幅を検出する。`run_0008_a2_r16/CONVERGENCE_CHECK.txt` は step 32300 で `NOT CONVERGED`。提案文面でよい。「壁温が落ち着いていた」を偽発火の根拠にしない。 |
| **Minor m3：結論・条件・完了記録の不整合** | **採用、修正対象を補足** | [A2 r64 評価:8](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/run_0010_a2_r64/EVAL_CONJ_L2.txt:8) は `Q_up/Q_tot` も PASS。提案された同期に加え、**参照水準数の誤記、§8 完了条件、§9 変更ログ、移動で切れるリンク**も直す。`L2` は3水準、`L3` は4水準である（[評価器:235](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:235)）。 |

結論: **追加計算を要求せず、上記を文書全体に反映してから、A の「事後改訂条件・FP64・node・流れ場固定」の限定結果として `accepted` へ移すことを推奨する。**

第 1 仮説: A2 r16 の差には forge の粗格子誤差が大きく寄与するが、不合格の原因を参照側から完全には分離できていない。 確度: **中**  
　根拠: `case/64.conjugate_pipe_wall/` 配下の `run_0011_a2_r16_df20`、`run_0009_a2_r32`、`run_0010_a2_r64` の `EVAL_CONJ_L2.txt:4` では、加熱区間熱流束の差が `0.017807 → 0.0043997 → 0.0010746`。一方、r16 の参照追加細分後も差は `0.018751` で、許容 `0.02` 未満。**格子依存は原因の示唆、FAIL は不確かさ込みの判定**として分ける。  
　反証条件: forge の保存場を固定したまま参照側の誤差を独立に減らし、差そのものが大幅に消えるなら、forge 側が主因という仮説は撤回する。合否だけの反転では、差の主因までは決まらない。

第 2・第 3 仮説: 追加しない。

判別 A/B: **既実施の参照3水準→4水準を採用し、再実行しない。** 変更点は参照の最大細分段のみ、forge は step 600000 の同じ保存場を固定する。判別は「PASS に変われば追加1段で不合格解消／U が上限内でも FAIL なら『追加1段で解消する』を棄却」。実測は後者。**『参照誤差一般を排除した』とは読めない。**

やらない方がよいこと: 合格するまで参照を細分すること、安全ガードを緩めること、A2 r16 の FAIL を PASS に変更すること。今回の限定閉鎖に `Df_scale=5` の追加対照は不要。

呼び出し側の前提への異議:

- 「参照4水準等の追加 A/B」は既に実施済み。`*_L3.txt:1` にも「参照4水準」とある。
- A1 は差単独の許容超過を維持してよい。[A1 L3:3](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/ab_levels3/EVAL_CONJ_a1_r16_L3.txt:3) は `差 − U = 0.0114126 > 0.01` でもある。ただし、A1 でも「参照側の誤差が無い」という断定までは不要。
- 移動前には、A/C の相互リンクも更新する。[後継 plan:11](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:11) は現在、同じディレクトリの A plan を参照している。C の未完了条件を保持し、旧判定は履歴として残す。

不足情報: **この限定閉鎖の判断に必要な不足はない。** 本番6本について、圧縮 CSV の読み込みだけを適応し、公式判定関数で全保存量の `PASS (converged)`（step 0–599999）と全系列の `STEADY`（50時点、12010–600000、末尾半分）を再確認した。G-if・G-cons の保存記録も全 PASS。run 索引は [case README](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/README.md:9)。中間場からの系列再生成は今回行っていない。

ファイル変更なし。**plan 未反映**。
