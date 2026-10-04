# 諮問ブリーフ: ON で液滴モーメントの更新を実現可能性の射影がほぼ打ち消す現象の原因と扱い

エスカレーション条件 2・4・7。plan: `plans/active/condensation-two-phase-default.md` §5.1 #4g3b (観測)、#5 (G2: ON の正式 corr-gate が射影で FAIL)。
成果物: `notes/investigations/2026-10-04-twophase-g3/evidence/G3B_RESULT.txt`、診断 h5 は AWS `case/16.nozzle_wys/run_0580_g3b_update_on/tp_update.h5` (ON)・`run_0581_g3b_update_off/tp_update.h5` (OFF) (1 更新、G2 の最終場 run_0561/run_0567 から、build 9be0d73c)。

## 観測事実
- 1 更新の全域収支 (ΣV·): ON の Q1 は制限後の増分 −0.771 に対し射影 +0.733、Q2 は +7.64e-9 に対し −7.72e-9。OFF の射影は Q1 −2.2e-32 (実質 0)。
- ON で Q1 の射影が働いた節点 4953、位置: x 7.7〜91.6 mm (中央 54.5)、壁距離 0.001〜0.053 mm (中央 0.013 mm) = 壁直近の層。相対補正 abs(ΔQ1)/abs(Q1) 中央 0.0018、95 % 1.49。OFF は 136 節点 (壁距離中央 0.079 mm)。
- 段階別の実現可能性 (射影が働いた節点の最大値; Q3 = g/(4/3·π·1000) と仮定した比):
  ON: upd_start まで Q1²/(Q0·Q2) ≤ 1.000000、Q2²/(Q1·Q3) ≤ 1.0023 (Q3 の換算の仮定のずれと思われる) → **tp_commit (二相の非分割更新の commit) の直後に Q1²/(Q0·Q2) 最大 1.045、Q2²/(Q1·Q3) = inf (g = 0 で Q1・Q2 > 0)** → tp_renorm・passive_floor でもそのまま → 射影の後 1.000001。
  OFF: cm_commit (モーメントの更新) の直後に inf、add_rho で Q1²/(Q0·Q2) 7.0e6 → 射影の後 1.00001 (ただし値が微小で収支への寄与は ~1e-32)。
- G1 (iii): 壁直近帯では ON の液は微小な壁向きの乱流流束を持つ。#4j: ON のみ液が閾値を超える節点は壁直近 (中央 0.054 mm)、そこは未飽和 (S 中央 0.53) で液滴は蒸発中。
- 報告量 (7 量) は G2 で STEADY・レシピ非依存。補正ゲートの commit/floor/clamp ≤ κ、射影は 1.0〜1.3e-6 > κ (記録のみの扱い)。

## 前提 (コード)
二相の更新 `tp_vl_update` (`solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh`)、DPLUR (`condensationTransport_d.cu`)、再正規化を液・Q に掛ける (`speciesTransport_d.cu` species_renormalize_liquid_d)、射影 `condensationRealizability_d.cuh` cond_realizability_clamp_f_d、凝縮ソース (蒸発の一様 ṙ 形、`condensationSourceKernels_d.cuh`)。

## 仮説 (私の見立て、未検証)
H1: 二相の更新は各モーメント (g, Q2, Q1, Q0) をそれぞれ別の対角 (乱流の対角は共通だが、ソースのヤコビアンがモーメントごとに違う) と別の制限で解くので、更新後の組は「あり得る組」の凸結合にならず、蒸発が進む壁直近で破れる。
H2: 蒸発で g を 0 にする制限 (θ や床) が g だけに掛かり、Q1・Q2 に同じ比率で掛からない (g = 0 で Q1・Q2 > 0 の観測)。
H3: 再正規化の係数の液・Q への適用。

## 問い
1. 原因の最有力はどれか、判別する最小の A/B (保存場から 1 更新、または短い run) は。
2. 既定化の判断への影響: 報告量は STEADY・レシピ非依存だが、毎更新で射影がモーメントの更新を打ち消す状態は「固定点」と言えるか。既定化を止めるべき欠陥か、受け入れて記録すべき近似か。止めるなら最小の修正方針 (例: モーメントに共通の制限係数、蒸発で g を消すときに同じ比率で Q を消す、更新を実現可能性を保つ形に) は。
3. G3 の合否 (局所の数値補正ゲート) にこれをどう反映するか。
読んでよいファイル: 上記、plan 2 本 (`condensation-two-phase-default.md`、`condensation-two-phase-transport.md`)、`methods/condensation.md` (§7c と蒸発・射影の節)、`solver_density_cuda/cuda_forge/` の上記ファイル、`notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`。
