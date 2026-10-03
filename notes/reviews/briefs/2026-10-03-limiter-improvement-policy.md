# 諮問ブリーフ: リミッタの定常収束停滞 — 改善方針 (§4.5) の点検

エスカレーション条件 1 (§4 の方針を新規に書く)。plan: `plans/active/limiter-inlet-column-oscillation.md` §4.5・§5.1 #7–#9。
文献調査: `notes/investigations/limiter-unstructured-convergence-survey.md` (サブエージェントの調査。主張ごとに確認度合いの印あり)。
これまでの実測: 同 plan §3.1、§5.1 #2r・#4r・#5er。関連決定: `plans/active/limiter-config-simplify.md` §4.3 (ψ 凍結は入れない)・§4.7 (既定 limiterScaled 1 / K 0.05)。
現在仕様: `methods/limiter.md` (SU2 に関する誤記を本日訂正)。実装: `solver_density_cuda/cuda_forge/limiter_d.cu:344-384`、`limiterFunctions_d.cuh:43`。

## 問い
1. 方針 1 (ε の局所 h_i 依存をやめる: (a) 領域一定 ε̂ / (b) Wang 型) は、観測 (入口列で ε̂ ≈ 1 ulp、極値判定で ψ が 0/1) への対処として筋が良いか。forge が局所 h_i を入れた経緯 (limiter-config-simplify / convection-node-wall-reconstruction の §4.11–§4.13) と、それを外したときに失うもの (メッシュ寸法への不変性、細かい格子での平滑化の消失を狙った設計意図?) を確認してほしい。
2. (a)(b) の具体形と既定値の決め方 (ℓ や K' をどう選ぶか、SU2 の無次元化時の相対 ε 0.011 を目安にしてよいか)。精度ゲート (Sod 厳密解・近傍逸脱・標準ケース) と case/16 の床で、事前登録すべき合格条件。
3. 方針 2 (ψ 凍結) を、limiter-config-simplify §4.3 の決定を覆して入れる根拠として、この調査で足りるか。入れるなら FUN3D 型 (破綻面だけ再計算) と SU2 型 (保持) のどちらが forge (陰解法・restart の運用) に合うか。
4. 優先順 (ε → 凍結) は妥当か。逆にすべき理由はあるか。
読んでよいファイル: 上記すべて、`solver_density_cuda/cuda_forge/` 配下、`procedures/verification/README.md`。
