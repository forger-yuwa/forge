# 諮問: 方向別の擬似 dt (`lineDtDirectional`) が冷却壁の格子で発散する理由

日付 2026-10-08。諮問先 codex (diagnose)。エスカレーション条件 2 (同じ対処で 2 回発散) と 4 (原因を書く前)。
ユーザ: 「方向別の刻みというのがなぜうまくいかないのか、すまんが追求したい」。
plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #27 (全文)、`plans/accepted/time_integration-line-implicit-viscous-v2.md` (2026-09-03 の case/45 の記録)。
作業ツリー `/home/sano/work/forge-integ-1005` (commit 3b8966f1)。forge は FP64 のビルド (`~/forge-wallfit-bin-fp64`、AWS)。

## 観測事実

共通の設定 (300 K 等温壁の冷却ノズル case/45、run_0183 の res_100000 からビット一致で restart、リミッタの基準値は run_0183 の値に固定):
node・軸対称・`nodeWallDirichlet 1`・SLAU (slauWallNormalChi 1 自動)・2 次 (convMethod 1、limiter 2)・SST (dilatationCorrection 2、katoLaunder 1)・
陰解法 block-DPLUR (timeIntegration 11)・`nStepInner 5`・`implicitRelax 0.7`・`lowMachPrecond 0`・unsteady 0・燃焼ガス TP。
格子は ni 4719 × nj 121 の構造 (node = i·121 + j、j = 120 が壁)、近壁は壁法線、第一層 y1/局所半径 ≈ 3.4e-7、近壁の縦横比 (流れ方向の幅 / 壁法線の幅) は縮流部〜スロートで 3500〜4400。

| run (`case/45.isobutane_m6_d155/` AWS) | 設定 | 結果 |
| --- | --- | --- |
| run_0191 | point (lineImplicit 0)、cfl 4 | 40000 step 安定 (NaN なし)。rms_roOmega は上下するが有界 |
| run_0200 | `lineImplicit 1` + `lineDtDirectional 1`、cfl 4 | 65 step で非有限 (detectNaN が 66 で停止) |
| run_0202 | 同、**cfl 1** | 129 step で非有限 |
| run_0203 | run_0200 と同じ設定の再実行 + `FORGE_DUMP_LEDGER` (列 36〜52 × 壁から 0〜30 層、527 節点、毎 step の状態) | 66 step で非有限 (run_0200 と同じ = 決定的) |
| run_0201 | `lineImplicit 1` だけ (方向別 dt なし)、cfl 4 | 走行中・安定 (1 step のコストは point の約 2 倍) |

- ライン: `[lineImplicit] lines=4719 covered CVs=570999/570999 (100.0%) maxLen=121` (1 列 = 1 ライン、壁から軸まで)。
- 壊れた場所 (run_0200・0202 とも): x/r_t −10.8〜−10.2 (列 39〜49) のライン 11 本で、ρ が壁から軸まで全層で非有限 (ライン解が 1 本丸ごと壊れる)。
- 帳簿ダンプ (run_0203) の振れ始め:
  - step 2 で最大の相対変化は ρ 4e-4 (列 36、壁から 6 層目)。
  - step 3〜8 で最大は**壁から 1 層目** (壁節点のすぐ内側) の列 37〜44: Uy (半径方向速度、|Δ| / |Ux|) 1.6e-2 → 8e-2、ρ・T・ρω が 1e-3 → 1e-2。
  - step 8 の ρ の変化の上位 12 点は全て列 37〜44 の 1〜2 層目。
  - その後は**周期およそ 16 step の、振幅が育つ低周波の振動** (1 step ごとの符号反転ではない)。列 44・5 層目の ρ の増分の符号は `++++++++-------++++++++--------+++++++++------++++++++++------++`、
    ρ は (4 step ごと) 31.07, 31.19, 31.35, 31.13, 30.97, 31.78, 32.37, 30.30, 29.82, 36.9, 47.4, 27.7, 25.2, 50.0, 59.9, 21.3, 24.0。
  - step 24 以降は Ux・Uy の 1 層目の変化が O(1) を超え、step 56〜65 で 12〜14 層目まで広がる。
- 形状 (壊れた列に特有なもの):

| 列 | x/r_t (壁) | 壁の傾き | 流れ方向の格子間隔 (r_t) | 1 層目の縦横比 | 壁から 10 層目の M | 中心の M |
| --- | --- | --- | --- | --- | --- | --- |
| 5〜20 | −12.3〜−11.8 | 0° | 0.037 | 3500〜4400 | 0.007 | 0.023 |
| 35〜49 | −11.0〜−10.2 | −1.7〜−5.5° (曲がり始め) | **0.061 (領域で最大)** | 3900〜4000 | 0.005〜0.006 | 0.03 |
| 120〜200 | −7.3〜−5.4 | −26〜−36° | 0.023〜0.034 | 4000 | 0.005〜0.009 | 0.03 |
| 1800 | 0.03 (スロート) | 1° | 0.0016 | 4300 | 0.35 | 0.92 |

- 実装 (`solver_density_cuda/cuda_forge/setDT_d.cu` 133〜180 行): 節点の cfl = 面の cfl_pln の**最大**。`lineDtDirectional` では `line_prev/next` に一致する内部面を最大から外す。
  壁節点の境界半割面は外さない (壁節点の Δτ は壁法線の音響で縛られたまま)。1 層目以降の Δτ は流れ方向の面だけで決まる (おおよそ Δx/(|u|+c)、縦横比の分だけ伸びる)。
- 過去 (2026-09-03、`plans/accepted/time_integration-line-implicit-viscous-v2.md`): case/45 の断熱・y+≈2 の格子 (別格子、FP32) で、directional は cfl 6・8 とも step 25 で発散。
  種は x/r_t 71〜89 の下流域・全断面 (壁至近 3 %・近軸 1 %) で、「既知の streamwise 内部モードが Δτ の上がった下流域で先に点火した」と記録。
  同じ格子で line だけ (directional なし) は 1 step 1.81 倍で末尾 roe −11 %、ni2 (nStepInner 2) は発散 (off-line lag に sweep ≥ 3 が要る)。
  case/39 (DDES、dual-time) では directional は cp4〜cp8 で安定 — dual-time の BDF の対角が保護していた、と推測されている (未確認)。

## 期待値

- 定常解は Δτ の取り方によらない。方向別 dt は、ラインで厳密に解く方向の λ を Δτ の制約から外すので、壁近くの遅い過渡 (縮流部の壁から 1〜60 層目に残る連続の残差、低 M × 高縦横比) を速めるはずだった。

## 実施済みの操作

- 上表の 4 本。CFL を 4 → 1 に下げても同じ場所で発散 (step 65 → 129)。
- point の cfl 8 (run_0199) は発散しないが、スロートの壁から 7 層目で残差が 60〜120 倍に張り付いた (別の現象)。

## 仮説 (確かめていない)

1. **off-line (流れ方向) の lag の不安定**: Δτ が縦横比倍に伸びると対角の V/Δτ が小さくなり、DPLUR の off-line の Jacobi 的な緩和 (nStepInner 5、relax 0.7) が低周波のモードを増幅する。流れ方向の格子間隔が最大の列で Δτ が最大になるので、そこで先に点火する。
2. **壁節点と 1 層目の Δτ の不整合**: ラインの解で、Δτ の小さい壁節点 (Dirichlet) と Δτ が数千倍の 1 層目が結合している。1 層目の更新が壁節点の拘束と食い違う。
3. **LHS と RHS の不整合 (defect correction)**: LHS は 1 次の Jacobian、RHS は 2 次の SLAU (低 M では SLAU の圧力の散逸が LHS と違う)。V/Δτ が小さいと不整合が増幅される。壁近くは M ≈ 0.005。
4. 2026-09-03 の下流の発散と同じ型 (定常では Δτ を伸ばすと保護が無い) が、この格子では縮流部の近壁で先に出た。

## 問い

1. 観測 (1 層目から、Uy 主体、周期 16 step の成長、流れ方向の格子間隔が最大の列、CFL によらない) と最も整合する機構は何か。上の仮説のどれか、または別のものか。
2. 機構を切り分ける最小の A/B は何か。既存のスイッチ (nStepInner、implicitRelax、1 次の RHS (convMethod 0)、limiter、slauWallNormalChi、FORGE_DUMP_LEDGER の faces) だけでできるものを優先し、コードの変更 (例: 方向別の Δτ の伸びに上限 R を付ける、壁の 1 層目だけ除外しない) が要るものは分けて示す。
3. 収束を速める目的にとって、この方向 (directional dt の安定化) を追う価値はあるか。line だけ (1 step 約 2 倍) で得をするには、1 step あたり 2 倍以上の収束の速さが要る。

## 読んでよいもの

- 上記 2 つの plan、`solver_density_cuda/cuda_forge/setDT_d.cu`、`timeIntegration_d.cu` (line の Thomas と DPLUR の対角・sweep)、`methods/time_integration/implementation.md` の「line-implicit」節、`procedures/solver-settings.md` の「lineImplicit」
- `case/45.isobutane_m6_d155/README.md` (run_0183・0190〜0203 の行)
