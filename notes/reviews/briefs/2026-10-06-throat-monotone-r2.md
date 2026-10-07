# 諮問: joint 壁のスロート直後 r″ の山を単調拘束で消す方針と、その検証計画 (§4・§6) の妥当性

日付 2026-10-06。諮問先 codex (diagnose) — `~/.config/forge/diagnose-backend` = codex。
エスカレーション条件 1 (plan §4 設計方針・§6 検証計画を新規に書く)。
plan: `plans/active/tooling-nozzle-throat-monotone-r2.md` (全文を読むこと)。作業ツリー `/home/sano/work/forge-integ-1005`、ブランチ `feature/nozzle-wall-fit-and-pipeline`、commit c821b71e + 本 plan (未 commit)。

## 観測事実 (形状のみ、CFD なし)

出典はすべて `case/45.isobutane_m6_d155/throat_r2_explainer.py` (commit c821b71e) の出力 `throat_r2_explainer.json`。対象は最終問題 `problem_d155_ns_finemesh_recal_final.yaml` の design_chain。

- 現行の joint 壁 (λ 1e-9、始点 r=1・r′=0・r″=1/R=0.5 をハード拘束) の r″:
  - [−0.5, 0) の上流 Hermite は 0.411 → 0.5。
  - x = 0.01435 で最大 0.51016 (1/R を超えるのは x 0〜0.0336)。
  - r‴ は [0, 0.05) で +1.02 → −0.92。以降は単調に減少し、x ≈ 2.15 に浅い極小がある。[0, 3] の r″ の極値は 2 個。
- 物理壁 (run_0117 の δ_r、`PhysicalNozzleWall` 解析経路): r″ の最大 0.51135 (x 0.01435)、スロート x_t −0.00155。
- MOC 壁点の区間平均曲率 (Δtan θ / Δx): 0.5151, 0.4840, 0.4748, 0.4641 …。
  - 第 1 点 x₁ = 0.02525 の θ₁ = 0.74519°、円弧角 atan(x₁/R) = 0.72330°。
  - 第 1 区間の自己整合 c₀ (割線角 − 平均 θ) = 0.0068°。
- λ 依存 (現行の拘束): λ 1e-10 / 1e-9 / 1e-8 で次のとおり。
  - r″ の最大: 0.5212 / 0.5102 / 0.5028
  - 点上の最大ずれ (角度): 0.0073 / 0.0124 / 0.0195°
  - 点上の最大ずれ (半径): 0.32 / 0.60 / 0.99 µm (r_t 76.654 mm)
- 試し C (同じ目的関数・拘束に、台が [0, 1.5] にかかる r‴ の B-spline 係数 ≤ 0 を追加。有効制約法):
  - 全 λ で r″ の最大 0.5000。
  - 点上の最大ずれ: 0.0219 / 0.0219 / 0.0220°、1.04 / 1.05 / 1.11 µm。いずれも第 1 点 x = 0.025 で生じる。
  - 有効制約は 4 本、解き直しは 5 回 (λ 1e-9)。
  - 現行との差 (λ 1e-9): 半径最大 0.50 µm (x ≈ 0.06)、壁角最大 0.0105°。0.05 µm を超えるのは x 0.017〜0.12、x > 0.3 では 1e-3 µm。
  - 平らな区間 (x < 0.03) から下降に移る所で、r‴ が 0 → −0.9 (x 0.04〜0.06)。
- 計算格子 (生産 NS/Euler G1) の壁節点間隔は、スロートで 0.0155 r_t (1.19 mm)。山の区間には 2〜3 節点。

## 期待値と出典

- 生産判定の枠組み: verification-m6 §5.1 #15・#16、CFD ピン plan の V3 (非劣化、評価器 `eval_wallfit_euler.py --fixed-coef`、Δq 表)。
- 過去の A/B: V0 → V4 (始点の r′・r″ 自由、角度 0.09° 級) で、試験部の P 傾き η0 が +0.154 → −0.231 %pt に動いた (交絡あり)。

## 実施済み

形状の試算のみ。実装・CFD は未実施。

## 仮説・方針

plan §4・§6 のとおり。不等式は凸包性による十分条件。区間は [0, 1.5]。λ は 1e-9 のまま。形状ゲート S1〜S7 → Euler A/B (G1 格子、各腕 3 回、18000 step、IC は run_0114 最終場から: A は restart_field、B は interp_field) → 非劣化で生産化。NS は再計算しない。

## 問い

1. 方針: 始点 r″ = 1/R を保ったまま単調拘束で山を消すことは、MOC 点群との関係で妥当か。
   - 第 1 点に 0.022° のずれを集中させることの是非。
   - 代案 (例: 第 1 区間の点の重みを下げる、λ をスロート近傍だけ強める、[0, b] の b の選び方) と比べてどうか。
   - 形状上の要求はユーザ決定で、流れの改善は主張しない。
2. 試し C の x 0.04〜0.06 の r‴ の急変 (0 → −0.9) は、山と比べて「滑らかさ」の観点で後退ではないか。ゲート S に r‴ (または r⁗) の条件を足すべきか。
3. S1〜S7 の閾値は、試し C の値を見てから決めている (S4 の 0.02°・2e-5 r_t、S5 の 0.001°)。結果を見てから合格条件を作っていることになるか。どう直すべきか。
4. E1 の IC の非対称 (A は restart_field、B は interp_field) は判定に効くか。18000 step 後の末尾で比べるので無視できると見ているが、別の手当て (例: 両腕とも interp_field、または座標差 ≤ 1e-6 r_t を許して index コピーする restart オプション) が要るか。AGENTS.md は同一メッシュへの interp_field を禁じている。
5. E3 の非劣化判定 (V3 と同じ Δq・U) は、今回の形状差 (0.5 µm、0.01°) に対して検出力が足りるか。各腕 3 回は妥当か。
6. 生産化で NS を再計算しない判断は妥当か。どの条件なら NS が要るか。

## 読んでよいもの

- 上記 plan
- `plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md` (§4.0、§6 V3、§9)
- `plans/active/verification-m6-axis-wave-mesh-su2.md` (§5.1 #15・#16、§9 2026-10-05)
- `design/forge_design/geometry/wall_axismach.py` (`joint_fit_wall`、`JointFitCFDWall`、`PhysicalNozzleWall`)
- `design/forge_design/evaluate/runner_axismach.py` (`design_chain`)
- `case/45.isobutane_m6_d155/{throat_r2_explainer.py, eval_wallfit_euler.py, euler_grid_ab.py, problem_d155_euler_pin_G1_recal.yaml, README.md}`
- `solver_density_cuda/tools/restart_field.py`

run データは手元に無い (上の数値がすべて)。編集は禁止。
