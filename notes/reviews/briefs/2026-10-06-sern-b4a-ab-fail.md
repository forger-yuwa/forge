# 諮問: B4a 判別 A/B が A・B とも FAIL — 仕様 (隣接間隔比「全方向・継ぎ目を含む」・cowl_side の h1e) と 19 ブロック構成の衝突 (2026-10-06)

関連 plan: `plans/active/tooling-sern-mesh-blocking.md` §4.14、§5.1 B4-0/B4-1/B4a、§6 (共通ゲート、特に 4 群)、§6.2 (事前登録の A/B・数値仕様・形状ゲート)。
前回: `notes/reviews/2026-10-06-sern-b4-full-blocking-diagnose.md` (B4a A/B を提案)、`notes/reviews/2026-10-06-tooling-sern-mesh-blocking-plan.md` (M2: 「全水準で実測比 ≤ 1.2」)。
エスカレーション条件 3 (§6 に事前に書いた合格条件が FAIL)。

## 観測事実 (作業ツリーの未 commit 実装、report は `notes/investigations/2026-10-06-sern-b4a/jm_{A,B}_report.json`)
- 実装: `case/46.sern_design/cad/hex_junction_model.py` (MOC 輪郭 B-spline、カウル外壁 = 内壁法線オフセット t・後端 x=L_cowl 平面、側壁厚 0.005 H、断面フィレット 0、全 x 位置を station 化、後流 [L, L+0.02H] Δx ≤ 1e-3 H、G 1.2 固定、検査は純関数 `first_layer_check`/`adjacent_spacing_check`/`wake_dx_check`/`shape_check`)。試験 `design/tests/run_sern_junction_check_tests.py` ALL PASS (30 項目)。
- A (t 0.02 H) 302 万節点・RSS 6.4 GB、B (t 0.005 H) 286 万節点・RSS 6.1 GB。**両者 VERDICT FAIL、落ちたゲートは同じ 2 つで板厚に依存しない**。
- PASS: 位相・タグ (タグ無し/二重 0)、負 Jacobian 0 (f64/f32)、skew max 0.567、AR max 679、後流 Δx 1.000000000001e-3 H、形状ゲート全部 (cowl_in/ramp 距離 0、法線板厚誤差 ≤ 3e-17 H、端面 |Δx| 4e-16 H、cowl_base 外側交点 y 誤差 0)、8 壁タグの第一層 (比 0.974–1.018、実層数 ≥ 期待、**壁法線方向の隣接比 max 1.196**)。
- FAIL 1 (隣接間隔比、§6.2 で「全方向・継ぎ目を含む」と書いた): 面を共有する全ヘキサ対で面直交辺長比の最大 **208**、>1.2 は A 76,388 / B 74,939 面対。継ぎ目別: U1b|U2b・U1a|U2a・CW1|CW2 208 (z=W/2 の継ぎ目で一様 NZ 0.033 H と、板厚を両端 h1 で切る NSW 分布が接する)、N_TOP|N_SIDE・N_BOT|N_SIDE 204、N_CORE|N_TOP/BOT 97 (リング外縁 0.011 H とコア e7/e11 の両端細分 h1、§4.11-3 で許容した伝播)、Sb_a|U2a 9.6、Sc_b|Sc_c A 4.5/B 9.55、その他 ≤ 3.7。**x 方向**は 1.366 (x≈0.784 H のコア内、`e7/e11 = max(mc, h5·…)` のクランプ切替で断面分布の x 微分が不連続)。
- FAIL 2: `cowl_side` を端面扱い (h1e 1e-3 H) にしたが生成は h1 → 比 0.16 = h1/h1e。実装者の指摘: A 区間で z=z_o 面上に `cowl_side` (y∈[y_cl,y_c]) と `sidewall_out` (y∈[y_c,y_r]) が同一平面で隣接し S の b 列の z 分布 (hz32) を共有するので、片方だけ h1e にすると他方か SW|Sb_de の継ぎ目が壊れる。
- 旧 B1 模型 (§4.11 PASS) は壁法線方向の層しか比を見ておらず、断面内の継ぎ目の比はもともと同程度に大きかったはず (未測定)。

## 期待値と出典
- plan §6-4 (2026-09-22): 「壁法線・端面法線・両端細分の各方向で隣接間隔比 ≤ 1.2 (区間境界を含む)」。
- plan §6.2 (2026-10-06、plan レビュー M2 を受けて私が書いた): 「実測の隣接間隔比 ≤ 1.2 (全方向・継ぎ目を含む、許容の上乗せなし)」、端面の第一内部点距離 h1e は cowl_base・sidewall_end・cowl_side。

## 問い
1. 隣接比ゲートの正しい範囲: §6-4 の方向限定 (壁法線・端面法線・両端細分) に戻すのが妥当か、それとも断面内の継ぎ目 208 は node 双対 FV の精度/安定に効くので構成 (分布の伝播・ブロック配置) を変えるべきか。変えるなら最小の変更は何か。継ぎ目の比を許容するなら上限値と根拠 (既存 PASS 模型・生産 run の実績で測れるもの) をどう決めるか。
2. `cowl_side` (カウル板の側端面、z = 側壁外面と同一平面) を h1 の壁として扱ってよいか (位相上 h1e は不可という主張の妥当性も)。
3. x 方向 1.366 のクランプ平滑化は中位で直してよいか。
4. 事前登録の分岐表 (A PASS/B FAIL、A も FAIL) はこの結果に当てはまらない (板厚非依存の仕様衝突)。どう記録し、何をもって B4b へ進むか。
