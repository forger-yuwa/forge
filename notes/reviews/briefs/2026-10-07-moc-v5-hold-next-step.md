# 諮問: MOC の軸処理の V5 (Euler) が登録どおり「保留」になった — 次の手

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 3 (事前登録の比較が判定不能) と 4 (plan に無い再実行へ進む前)。
plan: `plans/active/discretization-moc-axis-limit-and-corrector.md` (§6 V5 と「追加登録」、§9 の 2026-10-07 の記録を全文読むこと)。作業ツリー `/home/sano/work/forge-integ-1005` (HEAD 05519291)。

## 観測事実

- 評価器 `case/45.isobutane_m6_d155/moc_v5_euler_eval.py` (sha256 f7efa2ca…、登録 9a09c0a5 どおり、commit 3e20accf) を AWS で回した。出力 `case/45.isobutane_m6_d155/_band_ab/moc_v5_euler_eval.json`・`.log` (手元に複製済み)。各 run の時系列 `_band_ab/v5_series/<run>/wallfit_series_v5.csv` と `QUASISTEADY_wallfit_v5.txt`・`CONVERGENCE_VERDICT_segment.txt`。
- **総合: 保留 (前提不成立)**。7 本すべて (腕 B run_0143〜0145、腕 M run_0150〜0152、ISEN run_0153) が準定常の前提を満たさない。
  - 腕 B・腕 M: 6 量の多くが窓 (本段 6000〜18000、13 枚) の末尾 5 枚または全 13 枚で DRIFTING / OSCILLATING (一部 TRANSIENT-UNSETTLED)。exit_M_dev は E4 の絶対許容内の例外に当たらない (末尾 10 枚の最大絶対値 腕 B 2.1〜2.6e-5、腕 M 5.9〜6.8e-5、全 13 枚の幅 1.5〜1.6e-4、例外の上限 1.8e-5)。
  - 腕 B・腕 M の本段の残差は NOT CONVERGED stalled/plateau (rms_ro 0.9 dec、本段の約 3300 step で下げ止まり、ほぼ水平)。
  - ISEN: 残差は still converging (2.6 dec、全列が下がり続けている)。6 量とも DRIFTING (P 傾きは末尾 OSCILLATING)、exit_M_dev の全 13 枚の最大絶対値 0.0032。
- **参考値 (前提不成立なので判定ではない)**: 6 量の (D + 2SE)/Δq は M 波 0.29・P 波 0.18・オーバーシュート −0.72・出口規格化オーバーシュート −0.97・\|P 傾き\| −1.33・exit_M_dev 0.17 (負は腕 M が小さい向き)。出口コア M は腕 B 5.99999、腕 M 6.00003 (D +4.0e-5、\|D\| + 2SE 5.7e-5)。ISEN − 腕 M の差は参考値。
- 腕 B は 2026-10-06 の E′ 判定 (`plans/accepted/tooling-nozzle-throat-monotone-r2.md` §6 E′、`throat_mono_practical_eval.py`) で、準定常を前提にせず「許容幅未満」とされ、単調壁の採用に使われた。今回の V5 の前提を当てると、腕 B 自身が保留になる。
- **計算コスト**: Euler G1 (194,000 節点) は AWS で 4 本並列のとき約 85 step/s。本段 18000 step が約 4 分。延長は安い。
- 登録 (追加登録の最後の箇条): 「いずれかの run が固定窓の条件を満たさなければ、今回の V5 は保留とする。結果を見て窓を動かしたり、合格するまで延長したりしない。追加の検証は別に登録する。」
- 関連の経験 (メモリ): 残差のプラトーが停滞か有界な振動かは、保存量の同期・傾きの有意性・振幅で測る (`residual-plateau-bounded-oscillation`)。延長 run は再開直後に残差が跳ねるので、A/B の両腕を同じ回数だけ再開して揃える (skill forge-aws-run §3)。

## 期待値と出典

- V5 の目的 (plan §6 V5): 生産 Euler 格子で、新しい壁が今の単調壁より悪くないこと (E′ の許容幅) と、出口較正を変える必要があるか、IC 依存が無いか。
- 腕 B・腕 M の Euler は、収束済みの run_0114 から 18000 step 回したもの。残差は約 3300 step 以降水平。

## 問い

1. この保留をどう扱うか。候補:
   (a) 新しく登録して 7 本を延長する (再開の回数をそろえ、本段の長さ・新しい固定窓・停止条件を結果を見る前に決める)。延長で準定常の前提を満たす見込みはあるか。腕 B・M は残差が 3300 step 以降水平なので、量の DRIFTING / OSCILLATING が有界な振動 (リミットサイクル) なら延長しても STEADY にならない。
   (b) 量の時系列が有界な振動かゆっくりした過渡かを先に測る (CFD 0 step、手元の時系列で)。有界な振動なら、準定常の前提の決め方 (振動を平均±振幅で扱う、など) を新しく登録し直す。
   (c) ISEN だけ延長する。
   (d) ユーザ判断に回す (E′ と同じ「実務判定」で進めるか)。
2. 手元の時系列 (`_band_ab/v5_series/`) から、腕 B・腕 M の量の変化はゆっくりした過渡か、有界な振動か。判別の方法と、その結果で上の候補がどう分かれるか。
3. どの候補でも、新しい登録に書くべき条件 (延長の長さ・窓・停止・再開の揃え方・準定常の判定の決め方・参考値の扱い)。

## 読んでよいもの

- 上記 plan、`plans/accepted/tooling-nozzle-throat-monotone-r2.md` (§6 E1〜E4・E′、§9)
- `case/45.isobutane_m6_d155/_band_ab/{moc_v5_euler_eval.json, moc_v5_euler_eval.log, moc_v5_ic_inspection.json, wallfit_series_v5.json}`、`_band_ab/v5_series/` 以下
- `case/45.isobutane_m6_d155/{moc_v5_euler_eval.py, moc_v5_euler.py, run_moc_v5_euler.sh, eval_wallfit_euler.py, throat_mono_judge.py}`、`solver_density_cuda/tools/check_quasisteady.py`
- `.claude/skills/forge-aws-run/SKILL.md` (§3 restart)

編集は禁止。
