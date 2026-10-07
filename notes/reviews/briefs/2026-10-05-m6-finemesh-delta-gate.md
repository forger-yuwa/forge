# 諮問: #11 の「細分前後で δ_E(x_F) の差 ≤ 1 %」が 1.077 % で FAIL — 細分メッシュを生産格子として ② へ進めてよいか

日付 2026-10-05。諮問先: codex (diagnose) (`~/.config/forge/diagnose-backend` = codex)。エスカレーション条件 3 (plan §5.1 #11 に事前登録した許容差との比較が FAIL)。
plan: `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` §5.1 #11・#11a、§6 生産 pass、§9 末尾 3 行 (前回諮問 `notes/reviews/2026-10-05-m6-finemesh-divergence-diagnose.md` の採用と #11a の結果)。
commit: ブランチ `feature/nozzle-wall-fit-and-pipeline` 最新 (AWS も同 commit)。バイナリは run_0089〜0103 で同一 (main 統合ビルド)。

## 観測事実

- #11a (事前登録どおり): 起点 run_0098 の本段開始場、2 次のまま cfl 5 は step 35 で発散 (run_0101)、cfl 1 は 1000 step 健全 (run_0102)。登録解釈「CFL 依存、メッシュ修正必須説を退ける」。
- #11 ① の続行 `case/45.isobutane_m6_d155/run_0103_ns_finemesh_pass_cfl1`: 同じ起点、cfl 1、60000 step (累積 CFL を run_0092 の 5 × 12000 にそろえた)、5000 step ごと保存。
  - `check_convergence.py`: NOT CONVERGED (stalled/plateau)、低下 rms_ro 2.1・roUx 2.3・roUy 2.2・roe 2.1・roK 2.8 (falling)・roOmega 2.6 桁 — run_0092 (現格子) も plateau だった。
  - `check_quasisteady.py` (既定量): machmax・pmax STEADY、shock は衝撃波なしで非有限 (適用外)。
  - δ_E(x_F) (E 法 `extract_and_merge(band_select="edge")`、x_F 95.208、Euler 参照 run_0086): res_40000〜60000 = 0.733123 / 0.733128 / 0.733109 / 0.733116 / 0.733096 (r_t 単位)、5000 step ごとの変化 +0.0006/−0.0025/+0.0010/−0.0027 % (登録 ≤ 0.1 %/累積 CFL 20000 → 合格)。
  - 壁解像 `check_wall_resolution.py --over-frac 5`: **PASS** (y1+ > 1 の面積 3.6 %、最大 1.50 は index 0)。run_0094 は FAIL 29.4 %・最大 13.39。
  - **細分前後**: 同じ物理壁 (run_0092 の壁、r_t 76.7531 mm、k_f は c2pin_solve_pass2.json) で δ_E(x_F) run_0103 res_60000 = 0.733096、run_0094 res_6000 (現格子 ni 1250 × nj 97、wall_first_frac 4.5e-5、cfl 5、run_0092 から通算 18000 step) = 0.725282 → **差 +1.077 %** (登録 ≤ 1 % → FAIL)。
- 格子差: ni 1250 → 2000 (`throat_refine` 3 → 4、`throat_width` 1.5 → 3)、第 1 セル 4.5e-5 (一様) → 上流 1.3e-5・スロート 4.5e-6・下流 1.3e-5、nj 97 同じ (`_radial_fracs` が断面全体の等比を解き直すので軸側間隔 0.0752 → 0.0888 r_w)。cfl も 5 → 1 で違う (定常解は CFL に依らない前提)。
- 換算: δ_E 1 % は x_F で 0.0073 r_t = 0.56 mm (出口半径の登録 0.775 ± 0.1 mm より大きい)。② は細分格子の δ_E で k_f・r_t を解き直すので、生産設計は細分格子の δ に合わせて作られる。

## 当方の仮説・提案

- 1.08 % は「現格子の誤差」(壁解像 FAIL の格子が δ を 1 % 過小評価) を表し、細分格子が格子収束していることの証明ではない。ゲートの意図は「細分で δ_E が大きく動かない = 格子依存が小さい」の確認だったが、2 水準では細分格子側の誤差は分からない。
- 案 P: 細分格子を生産格子として ② へ進む。1.08 % を「格子依存の実測 (粗 → 細)」として記録し、ゲートは FAIL のまま残す。
- 案 Q: 3 水準目 (さらに細かい格子、例えば第 1 セル半分・ni 3000) で δ_E を測り、細分 → 超細分の差が 1 % 未満 (望ましくは粗 → 細の差の 1/2 未満) を確かめてから ②。
- 案 R: 差の成分を切り分ける (ni だけ・第 1 セルだけの格子) — コストは安い (1 本 ~4 分、cfl 1 × 60000 step)。

## 問い

Q1. このゲート FAIL の解釈と、② へ進む条件をどう決めるか (P/Q/R)。結果を見てから合格条件を作らない形で。
Q2. CFL 5 (粗) と CFL 1 (細) の違いが δ_E の差に混ざる可能性 (定常解が CFL に依存する経路 — 例えば plateau の振動の平均位置) を排除する必要があるか。必要なら最小の確認 (例: 粗格子を cfl 1・60000 step で回し直して δ_E を比べる)。
Q3. ③ (最終 NS) の合格条件のうち「細分前後で δ_E(x_F) の差 ≤ 1 %」をどう扱うか。

## 読んでよいもの

上記 plan、`case/45.isobutane_m6_d155/README.md`、`case/45.isobutane_m6_d155/{run_finemesh_pin.sh,problem_d155_ns_finemesh_pin.yaml,problem_d155_ns_c2pin_final.yaml,c2pin_solve.py}`、
`design/forge_design/meshing/mesh2d.py`、`design/forge_design/feedback/deltastar_loop.py` (`extract_and_merge`)、`notes/reviews/2026-10-05-m6-finemesh-divergence-diagnose.md`。run はローカルに無い (上の数値が全て)。編集は禁止。
