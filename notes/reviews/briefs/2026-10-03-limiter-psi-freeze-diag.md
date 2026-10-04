# 諮問ブリーフ: ψ 凍結の診断オプション (因果の向きの試験) の設計点検

エスカレーション条件 6 (cuda_forge の数値経路に触れる編集の前) と 1。plan: `plans/active/limiter-inlet-column-oscillation.md` §5.1 #5b (設計案と事前登録の判定を記載済み)。ユーザ了承済み。

## 背景
前回諮問 `notes/reviews/2026-10-03-limiter-inlet-next-ab-diagnose.md` の「状態・BC・勾配から切り離して ψ だけを変えたときの残差応答」を実現する。
#5a で ρ の再現誤差は時点ずれと確定 (`face_probe.py --ro-lag`、4 窓とも ≤1.93e-7)。

## 設計案
`limiter_d.cu` の `limiter_d_wrapper` 末尾 (流れ 5 変数の ψ 計算後、reconT・diag の前) に、環境変数で有効化する診断:
- `FORGE_DIAG_PSI_FREEZE_CALL=N`: N 回目の呼び出しで `limiter_ro/Ux/Uy/Uz/P` (nCells_all) を device バッファへ保存。
- 以後の呼び出し: ψ を通常どおり計算した後、`ccx < FORGE_DIAG_PSI_FREEZE_XMAX` (省略時は全節点) の節点だけ保存値で上書き。
- 未設定なら何もしない (ビット一致)。ログを起動時・凍結時に 1 行。

## 問い
1. この介入は「状態・BC・勾配を変えずに ψ だけを固定」になっているか。`limiter_d_wrapper` が 1 外反復に何回呼ばれるか (定常陰解法 timeIntegration 11、block-DPLUR、nStepInner 5) と、ψ を読む他の経路 (陰解法の Jacobian、境界、periodic、passive) を確認して、呼び出し回数 N と凍結範囲の決め方が正しいか。
2. 凍結すると作用素が別の固定点を持つ。判定 (§5.1 #5b) の比の基準と、F-all の「L0 床の 3 倍以内」の基準は妥当か。10000 step で足りるか (L0 は L1 の床から 1.5 dec 下がるのに 48000 step)。
3. F-in の x 閾値 (−0.059 m で入口 2 列 x=−0.060/−0.05937 を含み 3 列目 −0.05875 を含まない) は妥当か。
4. 実装上の罠 (nCells と nCells_all、周期 node 経路、ghost、restart 時の呼び出し回数) は。
読んでよいファイル: plan、`solver_density_cuda/cuda_forge/limiter_d.cu`、`solver_density_cuda/main.cpp`、`timeIntegration_d.cu`、`convection/` 配下。
