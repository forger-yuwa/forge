# 諮問 (2 回目): 方向別の擬似 dt (`lineDtDirectional`) の発散 — 単因子の A/B 3 本でも止まらない

日付 2026-10-09。諮問先: diagnostician (Fable) と codex (diagnose) の両方 (ユーザ指示「diagnostician と codex 両方に聴いて」)。
エスカレーション条件 2 (同じ対処で発散 4 回) と 4 (原因を書く前)。
前回のブリーフ `notes/reviews/briefs/2026-10-08-linedir-divergence.md` と前回の codex の回答 `notes/reviews/2026-10-08-linedir-divergence-diagnose.md` を**先に全文読むこと**。
plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #27 (「方向別の dt の発散の追究」以降)。作業ツリー `/home/sano/work/forge-integ-1005` (commit bec196f6)。
ユーザ: 「方向別の刻みというのがなぜうまくいかないのか、すまんが追求したい」「続けてよい」。

## 観測事実 (前回のブリーフ以降)

全 run 共通: case/45 の 300 K 等温壁の冷却ノズル、run_0183 の res_100000 からビット一致の restart、リミッタの基準値は run_0183 の値に固定、
node・軸対称・`nodeWallDirichlet 1`・SLAU (`slauWallNormalChi` 自動 = 1)・limiter 2・SST (dilat 2、KL 1)・block-DPLUR (timeIntegration 11)・
`lineImplicit 1`・`lineDtDirectional 1`・cfl_pseudo 4・implicitRelax 0.7・lowMachPrecond 0・燃焼ガス TP・FP64 のビルド。

| run | 前回 (run_0203) からの単一の変更 | 結果 |
| --- | --- | --- |
| run_0203 | — (convMethod 1、nStepInner 5、implicitSolvePrecision 0) | 65 step で非有限 (step 66 で detectNaN) |
| run_0204 | `implicitSolvePrecision: 1` (行列の組立て・解法を double) | **65 step** で非有限。振れ始めの場所・大きさ・周期が run_0203 とほぼ一致 (下表) |
| run_0206 | `nStepInner: 15` | **65 step** で非有限。rms の推移も run_0203 とほぼ同じ (step 20: rms_ro 1.30e-4 vs 1.34e-4、step 40: 1.03e-3 vs 1.09e-3、step 60: 5.94e-3 vs 5.95e-3) |
| run_0207 | `convMethod: 0` (RHS を 1 次) | 72 step で非有限。開始時の残差は 1 次の作用素なので大きい (rms_ro 4.9e-5) が、その後の成長は同じ型 (step 40: 1.25e-3、step 60: 5.3e-3) |
| (前回) run_0202 | cfl_pseudo 1 | 129 step で非有限、同じ場所 |
| run_0201 | `lineDtDirectional` なし (lineImplicit だけ) | 10000 step 安定。質量の欠損 (Σ res_ro×2π) は 5000/10000 step で 1.9098 / 1.8202 で、point (run_0191) の 1.9168 / 1.8197 と同じ |

run_0203 と run_0204 の帳簿 (列 36〜52 × 壁から 0〜30 層、毎 step の entry 状態):

| step | ρ の最大の相対変化 (run_0203, float) | (run_0204, double) |
| --- | --- | --- |
| 4 | 2.86e-3 @ 列 43・壁から 1 層目 | 2.86e-3 @ 列 43・1 層目 |
| 8 | 1.09e-2 @ 列 43・1 層目 | 1.09e-2 @ 列 43・1 層目 |
| 16 | 9.46e-3 @ 列 43・1 層目 | 9.36e-3 @ 列 43・1 層目 |
| 32 | 7.60e-2 @ 列 44・2 層目 | 7.36e-2 @ 列 44・2 層目 |
| 48 | 4.55e-1 @ 列 48・6 層目 | 4.42e-1 @ 列 47・6 層目 |

列 44・5 層目の ρ (4 step ごと) — float: 31.07 31.19 31.35 31.13 30.97 31.78 32.36 30.30 29.82 36.88 47.40 27.70 25.23 50.01 59.86 21.30 23.99、
double: 31.07 31.19 31.35 31.13 30.98 31.78 32.34 30.30 29.88 36.82 46.34 27.52 25.68 50.48 60.13 21.21 22.16。

前回の帳簿の解析 (run_0203): step 3〜8 で最大の変化は壁から 1 層目の列 37〜44。Uy (|ΔUy|/|Ux|) が 1.6e-2 (step 3) → 8e-2 (step 6) と先に動き、ρ・T・ρω が 1e-3 → 1e-2 で続く。
step 24 以降は 1 層目の Ux・Uy の変化が O(1) を超え、step 56〜65 で 12〜14 層目まで広がる。振動は周期約 16 step で振幅が育つ (1 step ごとの符号反転ではない)。

帳簿の生データ (AWS、現在インスタンス停止中): `~/forge-wallfit/case/45.isobutane_m6_d155/run_020{3,4,6,7}_*/ledger.csv` (call,tag,node,field,value。tag は entry・after_eos_bc・res_before_conv・res_after_conv・res_after_rans_transport・res_after_species・res_after_sources・res_after_viscous・res_final)、
`ledger.csv.faces` (SLAU の面の記録: ro_L/R・P_L/R・Pf_L/R・U_L/R・c_hat・M_hat・chi・chi_mass・Vn_p・Vn_m・p_tilde_r・mdot・F_*・fx・limiter 等、列名は 1 行目)、節点 ID は `ledger_nodes.txt` (node = i·121 + j、j = 120 が壁)。
`res_nan_*.h5` と 20 step ごとの `res_*.h5` (全保存量の残差 res_* を含む) もある。

## 実施済みの判定 (事前登録どおり)

- 「float の組立てが主因」→ 外す (run_0204 が同じ位置・周期・成長)。
- 「off-line の lag (sweep 不足)」→ 外す (run_0206、nStepInner を 3 倍にしても同じ step で発散、残差の推移も同じ)。
- 「LHS (1 次 FVS) と RHS (2 次) の次数の不整合」→ 外す (run_0207、RHS を 1 次にしても 72 step で発散)。
- 「CFL によらない」とは書かない (試した 2 値 4・1 では安定化しなかった)。

## 仮説 (呼び出し側、確かめていない)

A. 陰的部分 (近似の Jacobian) が、壁と 1 層目の間の何かを表していない。Δτ が小さい間は V/Δτ が覆い、方向別の dt で 1 層目の Δτ が縦横比 (約 4000) 倍になると露出する。候補:
   - SLAU の壁隣接面の圧力差による質量流束 (`slauWallNormalChi` 1 で χ を面法線成分で組む; 圧力束は変えない) が LHS に無い。
   - 壁節点の行の扱い (運動量・エネルギーは単位行化、密度の行は残る — 前回 codex 指摘、`timeIntegration_d.cu` 1029・1053 行) と、1 層目の大きな Δτ の組み合わせ。
   - 等温壁: 壁節点の T を 300 K に固定する BC 適用 (after_eos_bc) が、1 層目の陰的更新と整合しない。
   - 軸対称の hoop の源項 (p/r) や SST の源項が LHS に無い (segregated)。
B. 壁節点の Δτ は境界半割面の音響で縛られたまま小さく、1 層目は数千倍。ライン解の中で隣接する 2 節点の時間対角が数千倍違うことそのものが、ライン上の低周波モードを不安定にする。
C. 2026-09-03 の case/45 (断熱・y+≈2・FP32) の発散 (下流 x 71〜89 の全断面で step 25) と同じ型 — 定常では Δτ を伸ばすと擬似時間の対角という保護が無くなり、近似の LHS で表していないモードが育つ (dual-time の DDES では BDF の対角が保護するので安定)。

## 問い

1. 4 本の結果 (精度・sweep 数・RHS の次数のどれを変えても同じ場所・ほぼ同じ step で発散、CFL 4 → 1 で 65 → 129 step に遅れるだけ) と最も整合する機構は何か。上の A〜C のどれか、または別か。
2. 次の単因子の A/B (既存のスイッチだけで、200 step 以内) を優先順に 2〜3 本挙げてほしい。候補の例: `slauWallNormalChi: 0` (明示)、`implicitRelax` (0.7 → 0.3)、等温壁 → 断熱壁 (同じ開始場)、乱流を凍結する設定があればそれ、別の対流スキーム (SLAU 以外)、`lineViscCoupling: 1`。それぞれの予測 (止まる / 止まらない) を事前に書けるように。
3. 帳簿の面の記録 (`ledger.csv.faces`) や各段の残差 (res_after_conv / res_after_viscous / res_final) から、CFD を回さずに切り分けられることは何か (例: 1 層目で成長を駆動している残差の項を特定する方法)。
4. コードの変更が要る対策 (例: 方向別の Δτ の伸びに上限 R を付ける、壁の 1 つ内側だけ除外しない) を試す価値はあるか、試すならどれを先に。高速化の目的 (欠損の時定数は point の cfl 4 で約 10 万 step) に照らして判断してほしい。

## 読んでよいもの

- 前回のブリーフと回答 (上記)、plan §5.1 #27、`plans/accepted/time_integration-line-implicit-viscous-v2.md`
- `solver_density_cuda/cuda_forge/setDT_d.cu`、`timeIntegration_d.cu` (DPLUR の対角・壁の行・line の Thomas・sweep)、`convection/convectiveFlux_slau_d.inc.cuh` (slauWallNormalChi)、`nodeWallDirichlet_d.cu`、`methods/time_integration/implementation.md` の「line-implicit」節、`methods/boundary.md` (等温壁)
- `case/45.isobutane_m6_d155/README.md` (run_0183・0190〜0207 の行)、`case/45.isobutane_m6_d155/cold_cfl.py`
