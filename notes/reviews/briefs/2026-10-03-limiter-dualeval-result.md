# 諮問ブリーフ: ψ 差し替えの残差二重評価の結果の解釈

エスカレーション条件 3 (事前登録の前提が満たされず判定不能になりうる) と 7 (解釈の確定前)。
plan: `plans/active/limiter-inlet-column-oscillation.md` §5.1 #5d (事前登録)・#5dr (結果、数値はすべてここ)。
実装: commit 05c5ce36 / aef6c1e1 の `solver_density_cuda/main.cpp` (pdeBeforeA / pdeAfterA)、`solver_density_cuda/cuda_forge/limiter_d.cu` (psi_override_d)。
前回諮問: `notes/reviews/2026-10-03-limiter-psi-freeze-diag-diagnose.md` (二重評価の提案と判定 D≥1/≤0.1)。

## 観測事実
#5dr のとおり。D(res_ro) ≈ 1.06 で前半・後半とも再現、雑音 A'−A は 1e-5。前提「組立を重ねても保存量がバイト一致」は崩れた:
内部の壁近傍の roe が 1〜2 ulp、ghost の roUx/ro/roY が変わる。

## 問い
1. 前提の逸脱 (壁近傍 roe の 1〜2 ulp と ghost) は、B−A を ψ の差し替えの効果として読むことを妨げるか。A'−A の雑音 (B−A の 1e-5〜1e-3) がこの逸脱の影響も含んだ上限になっている、という私の読みは正しいか。どの組立段 (壁射影 `enforceWallNoSlip`・EOS・BC の ghost 更新) が原因と考えられるか、コードで確認してほしい。
2. 判定: res_ro・roUx・roe・roY・roOmega で D ≥ 1 (両窓)、roK 0.94、roUy 0.33。事前登録の基準で「第 1 仮説 (入口 2 列の ψ の変動は、同じ状態の入口残差へ直接大きく寄与) を支持」と書いてよいか。
3. D ≈ 1 は「差し替えで残差が残差と同程度変わる」だが、方向 (減るか増えるか) は言っていない。B と A のノルム比 (S_B/S_A) を出すべきか (CSV から計算可能: S_B = ‖A + (B−A)‖ は直接無いので要追加?)。
4. 次の一手: 長期凍結 (F-in) へ進むべきか、それとも原因側 (入口境界列の極値判定が数 ulp で決まる #4r) を直接変える A/B (例: 入口境界節点を近傍集合から外す・入口列の ψ を 1 次化) に進むべきか。どちらもカーネル変更なので別 plan の要否も。
読んでよいファイル: plan、上記ソース、`notes/investigations/2026-10-03-limiter-inlet/*.py`、`case/16.nozzle_wys/README.md`。
