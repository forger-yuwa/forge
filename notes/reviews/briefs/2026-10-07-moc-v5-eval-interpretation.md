# 諮問: MOC の軸処理の V5 (Euler) の判定器の解釈 6 点 — 結果を見る前の確定

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 3 (事前登録の比較の合否条件) と 1 (§6 の解釈の確定)。**run は投入済みだが結果はまだ誰も見ていない**。ここで決めた解釈を plan §6 V5 に書いてから結果を読む。
plan: `plans/active/discretization-moc-axis-limit-and-corrector.md` (§6 V5 を全文読むこと)。作業ツリー `/home/sano/work/forge-integ-1005` (HEAD ba4b57ef)。

## 観測事実

- V5 の run (AWS、投入 2026-10-07): 腕 M = run_0150〜0152 (単調壁 + `moc_axis_limit: analytic`・`moc_corrector: converge`、IC = run_0114 res_6000 の検証付き番号写像、上限 6.21 µm)、IC 依存の確認 = run_0153 (同じ壁・格子、等エントロピー IC から段階起動)。比較相手は腕 B = run_0143〜0145 (単調壁、2026-10-06 完走、IC = run_0114 から番号写像、上限 1 µm)。段はどれも soft 3000 → 本段 18000、窓は本段 6000〜18000 の 13 枚。
- 判定器 `case/45.isobutane_m6_d155/moc_v5_euler_eval.py` と負例テスト `test_moc_v5_euler_eval.py` (69 件 FAIL 0)。起動 `run_moc_v5_euler.sh`、準備・IC の検査 `moc_v5_euler.py`。
- 乾式の IC の検査 (腕 M): 移動の最大 6.199 µm (x = 0.406 r_t、壁から 4 層目)、中央値 0.83 µm、99 % 点 5.64 µm。反転 0。C1〜C5 すべて成立。最近傍なら 5,948 節点で層を取り違える。
- 腕 B の E′ 判定 (`throat_mono_practical_eval.py`、monotone plan §6 E′) は、前提検査に準定常 (check_quasisteady) を**入れていない** (残差・壁・IC だけ)。13 枚の窓平均の時間変動を SE に入れる方法で、振動は SE で吸収した。monotone plan の E4 (各 run の単独判定) は「末尾 5 枚 STEADY でない量は判定に入れず保留」と近零量の「絶対許容内」を持つ。
- plan §6 V5 の前提の文言: 「全 run の判定区間の check_convergence が pass か既知の plateau、評価量の check_quasisteady (末尾 5 枚)。欠損・DIVERGED は保留」。
- 代用の run (run_0114 の 6 枚) で判定器を通すと、overshoot η0.1 が OSCILLATING、|出口コア M − 6| が DRIFTING だった。

## 実装担当が入れた解釈 (plan の文言に無い)

1. 準定常: 6 量と出口コア M のすべてで末尾 5 枚 STEADY を必須にした。OSCILLATING・DRIFTING は保留。近零量の「絶対許容内」は適用せず記録だけ。腕 B にも同じ前提を課すかは未確認 (腕 B の 3 本の各量の状態はまだ見ていない)。
2. IC 依存 (run_0153 と腕 M の差) を両側で判定 (\|D_ic\| + 2SE ≤ Δq)。片側だと ISEN の方が良い向きの大差を見逃すため。
3. 出口較正を 3 区分: \|D\| + 2SE ≤ 1e-4 なら据え置き、\|D\| − 2SE > 1e-4 ならやり直し、それ以外は判別不能で保留 (較正は変えない)。
4. 評価の X_E・X_F は 7 本とも腕 B の r1 の値を使い、各 run 自身の値との差 ≤ 1e-6 r_t を検査 (新しい閾値。腕 M と腕 B の X_F の差は 9.9e-9 r_t で、既存の e3_geometry の 1e-9 では止まる)。
5. 番号写像の上限: C1 (ic_index_map の 7 検査)・C2 (列の x 不変)・C3 (η = r/r_w を保つ縮尺)・C4 (壁の移動が設計壁の差で説明できる)・C5 (メッシュ品質) と、構造から決まる上界 max √(tol_x² + (\|Δr_w\| + 2 ulp)²) で決めた (= 6.2063 µm)。
6. ISEN (run_0153) は soft 3000 → 本段 18000 で、窓 6000〜18000 までに定常化するか分からない (腕 B・腕 M は収束済みの run_0114 から始める)。

## 問い

1. 準定常の前提をどう定めるか。候補: (a) 実装どおり全量 STEADY 必須 (腕 B にも課す)、(b) E′ と同じく前提に入れず、窓平均と SE で扱う (DRIFTING だけ保留にする等の中間案を含む)、(c) E4 の絶対許容内を使う。V5 は E′ の方法で比べると事前登録したこと、腕 B の E′ では準定常を前提にしていないことを踏まえて。
2. 解釈 2〜5 は妥当か。
3. ISEN が窓までに定常化しない場合の扱い (延長するか、IC 依存の判定を保留にするか、窓をずらすか) を、結果を見る前にどう決めるか。
4. 上の決定を plan §6 V5 に書くときの文言 (結果を見る前の登録として十分か)。

## 読んでよいもの

- 上記 plan、`plans/accepted/tooling-nozzle-throat-monotone-r2.md` (§6 E1〜E4・E′、§9)
- `case/45.isobutane_m6_d155/{moc_v5_euler_eval.py, test_moc_v5_euler_eval.py, moc_v5_euler.py, run_moc_v5_euler.sh, throat_mono_practical_eval.py, throat_mono_judge.py, eval_wallfit_euler.py}`
- `solver_density_cuda/tools/check_quasisteady.py`

run の結果ファイルは読まないこと (まだ出ていないし、結果を見ずに決めるための諮問)。編集は禁止。
