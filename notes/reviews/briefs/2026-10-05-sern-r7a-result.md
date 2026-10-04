# 諮問: SERN R7a (3D `L_sw` 2 水準 × 3 反復) の結果解釈と R7b への進め方 (2026-10-05)

関連 plan: `plans/active/tooling-nozzle-sern-chain.md` §5.1 R7a (事前登録と結果)・R7b、§8。前回諮問: `notes/reviews/2026-10-05-sern-next-priority-diagnose.md` (R7a の設計)。
原本: `notes/investigations/2026-10-05-sern-r7a/R7A_VERDICT.txt`。run: AWS `~/forge-r8/case/46.sern_design/run_1034`–`1042` (case/46 README)。

## 観測事実
- 条件: 3D g3・m6_on・生産 YAML (side_far farfield、種ごとの輸送)、リミッタ基準値固定、バイナリ eb0ea331 (main 77318d0e 系)、`FORGE_CUDA_BLOCKSIZE=128`。A = `L_sw` 0.8 (側壁 sidewall_in の x 最大 0.080910 m)、B = 1.0 (0.100433 m)。他の形状・格子は同一 (B は 803 節点多い = 延長した側壁の双子節点)。
- 初期場: A は run_1017 (同じ A 格子の生産 run) の最終場から restart_field。B は格子を作り直し同じ最終場から interp_field (最近傍距離 0、双子節点は両側とも同じ元の値)。
- A 3 本: 20000 step で窓条件 OK。B 3 本: 20000 step で窓条件未達 (C_T_with_shear 前窓差 2.6e-4、C_M 3.5e-3) → 事前規則で +20000 → 全量 OK。全 9 本 GATES PASS。
- 末尾 10000 step 平均 (B は延長 run): C_T_with_shear A 0.8705560 / B 0.8726419 (B − A +0.002086)、N_A 3.29e-6・N_B 3.59e-6 (N = max(3r, a, d, 1e-6))、E 6.89e-6 → |Δ| − E = 0.002079 > 0.002 → 事前規則の分岐 A。
- その他: C_T +0.003434、C_L +0.008975、C_M −0.213185 (A −1.6247 / B −1.8378)。各群の再現性は C_M で r 3–7e-5、a 2.6–2.8e-4。
- 参考 (既存): g3→g4 の格子差 G (同じ farfield・新輸送、#4f) は C_T_with_shear −0.00019、C_M +0.0167。§8 の総許容 C_T・C_L 2e-3、C_M 5e-2。

## 問い
1. 分岐 A の結論 (固定格子 g3・m6_on・二水準で識別できる) の言える範囲と言い過ぎ。余裕が検出目標の 4 % しかないこと、B の起点の補間 (双子節点) と 20000 step 過渡、を考えたときに追加の確認が要るか (例: B を A 側からでなく B 自身の段階起動で作った別起点の 1 本、または 4 本目)。
2. C_M の −0.213 (§8 許容の 4 倍、作動点別 C_M 窓は未定) をどう扱うか。L_sw が推力 +0.24 % と同時にピッチを大きく動かすこと自体は設計情報だが、格子差 (g3/g4) や領域感度の裏付けなしに B の C_M を信用してよいか。
3. R7b (3D 最適化への接続) に進む前の最小の確認は何か (格子差 G を B 形状でも取る、他作動点、L_sw の物理 station 化)。
