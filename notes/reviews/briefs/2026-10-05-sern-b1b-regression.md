# 諮問: B1b (3D 双対の境界半割面重心を double に) の回帰で C_M だけ事前閾値を超過 — 扱いと次 (2026-10-05)

関連 plan: `plans/active/tooling-sern-mesh-blocking.md` §5.1 B1b・R1 (事前登録と結果)、前回諮問 `notes/reviews/2026-10-05-sern-cowl-blunt-te-route-diagnose.md`。
原本: `notes/investigations/2026-10-05-sern-b1b/B1B_VERDICT.txt`。コード: `solver_density_cuda/mesh/gmshReader.hpp::buildMedianDual3D` (commit ef74e342)。

## 事実
- 原因の特定: 境界半割面の面重心だけが保存済みの float の primal 面重心を使い、内部双対面は節点から double で計算していた。修正は境界側も節点から double で計算 (1 か所)。
- 接続模型: 閉性 FAIL 866 CV → 0 (max 7.41e-5 → 5.55e-8)。生産 g3 (exact A): FAIL 408 CV (max 2.83e-4) → 0 (5.66e-8)。座標・接続・体積同一、境界半割面 154298 面が変化 (面ベクトル最大相対 0.22 %)。
- CFD 回帰 (事前登録): run_1050 最終場から旧格子 run_1067 / 新格子 run_1068 各 20000 step、同じ forge。両者 GATES PASS・STEADY・窓条件 OK。D = |差| + U (U = Σ max(a, d, 1e-6)): C_T 5.6e-6、C_T_with_shear 5.6e-6、C_L 2.4e-5 (閾値 5e-5 内)、C_M 差 −4.2e-5・U 4.81e-4・D 5.23e-4 > 5e-4。事前規則「超えたら記録して諮る」。
- 同型の前例: R7a の写像対照で C_M が U の振幅だけで閾値に達し「判定不能」とした (codex 2026-10-05)。

## 問い
1. この回帰の結論 (B1b の修正を生産の変換器として採用してよいか) と、C_M の扱い。
2. 2D の双対構築 (`buildMedianDual`) に同じ欠陥 (境界半割面だけ保存済みの面重心を使う) があるか確認すべきか、その最小の確認方法。
3. 次の工程 (B4: 有限厚カウルの SERN 全体ブロッキング) に進む前に B1b でほかに閉じるべきこと。
