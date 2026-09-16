# 09. Taylor-Green Vortex (3D)

3 次元 Taylor-Green 渦 (TGV)。全方向周期境界の立方体 $[0,2\pi]^3$、構造化 hex $32^3$ (節点 $33^3$)。
低散逸スキームの **運動エネルギー (KE)・エントロピー総和の保存性** を見るための検証ケース。

- メッシュ: [`mesh/Taylor-Green.geo`](mesh/Taylor-Green.geo) (`nx=ny=33`, z 方向 32 層 extrude)。
- 初期条件: `initial: "Taylor-Green"` (無次元, $M_0=0.4$, $\rho_0=1$, $P_0=1/\gamma\approx0.714$)。
  $u=M_0\sin x\cos y\cos z$, $v=-M_0\cos x\sin y\cos z$, $w=0$。
- 数値: pure KEEP (`solver: KEEP`, `keepDissipation: 0` = Roe 散逸無し・非散逸中心流束のみ), 陽解法 RK4 (`timeIntegration: 4`), 乱流モデル無し。

## 重要な設定上の注意

- **無次元・低圧ケースなので EOS フロアを下げること**: 既定の圧力フロア `pMin=1.0` Pa は本ケースの圧力場
  ($P\in[0.65,0.78]$) を全域クランプし、初期場を破壊する (step0 で $P\equiv1.0$)。`physProp.pMin/roMin/tMin`
  を小さく (例 `1e-6`) 指定する。floor の config 化は [`plans/accepted/thermophysics-eos-positivity-floor-config.md`](../../plans/accepted/thermophysics-eos-positivity-floor-config.md) 参照。
- **node モードのメッシュは node 設定で変換すること**: `discretization: "node"` のとき `convertGmshToForge` は
  median-dual を構築して primal を置換する。cell 設定で変換した h5 を node run に流用すると primal hex 上で
  解いてしまい不正 (双対 nCV=35937 ≠ primal hex 32768)。各 node run dir で自身の node config を使って変換する。

## 計算 run 一覧

全 run 共通: pure KEEP は `keepDissipation: 0`、陽解法 RK4、乱流なし、`pMin: 1e-6`、出力間隔 50 step。

全 run で全運動量 $\sum\rho\mathbf{u}\,V$ (理論=0) と全質量を保存量チェックに用いる。**node の res は周期マージ双対
セルを重複格納する** (VALUE 長 35937 = ノード数、一意 CV は 32³=32768) ため、後処理は各境界 CV を多重度
$m=2^{(境界方向数)}$ で割って一意 1 個分に補正する (`plot_ke_entropy_history.py` の `cv_weight`)。cell は補正不要。

以下は **cell 周期 partnerCellID device 転送バグ修正後** ([`plans/accepted/boundary-cell-periodic-conservation.md`](../../plans/accepted/boundary-cell-periodic-conservation.md)) の結果。
保存量チェックは全運動量 $\sum\rho\mathbf{u}\,V$ (理論=0) と全質量。**node の res は周期マージ双対セルを重複格納する**
(VALUE 長 35937 = ノード数、一意 CV は 32³=32768) ため、後処理は各境界 CV を多重度 $m=2^{(境界方向数)}$ で割る
(`plot_ke_entropy_history.py` の `cv_weight`)。cell は補正不要。

| run_* | 離散化 | 設定差分 | 主要結果 (多重度補正後) | 状態 |
| `run_0060_periodic_seam_cond_source` / `run_0061_periodic_seam_cond_source_oldbin` | **node 周期 seam の凝縮ソース積分体積の試験** (plan condensation-source-limiter-steady, codex result M6): 一様静止・過飽和 N2 (60 K, S≈5, CPG) を 1 step (timeIntegration 11, cfl 0.01) 回し、`Q0·V/dt` (=積分ソース/合併体積) を seam ノードと内部で比較。0060 = 新バイナリ (部分体積 `volumePartial_d`)、0061 = 旧バイナリ (forge-cond 2026-09-13, 合併体積) | **旧: seam/内部 = 2.000 (面) 〜 8.0 (角) の二重計上、新: 1.000 (分位 0.9999998〜1.0000002)**。最終バイナリ `run_0062_periodic_seam_cond_source_final` / 最終確定 `run_0063_…_final2` も 1.000000。seam の dt_local は部分体積由来で半分 (別件, 固定点には無関係) | active (seam 試験の記録) |
| `run_0065_periodic_tracer_seam` / `run_0066_periodic_tracer_interior` | node 周期 (33³ 箱) | **受動トレーサ `roXi` の周期 seam 試験** (plan thermophysics-cea-mole-fraction-species, codex result M5): 一様流 u=10 m/s (KEEP + advective gauge), RK4 dt 1.7e-4, 500 step, Xi₀=½(1+cos(x−x₀)) を seam 中心 (0065) / 内部 x=π (0066) に置く (判定 `run_0065/analyze_seam.py`, run dir 内) | 周期対 (1089 対 × 3 方向) の max\|ΔXi\| = 0、∫ρXi dV の相対変化 3.7e-8 / 5.0e-8、seam 版と内部版 (π シフト) の差 3.0e-7 (float 床)、流れは一様のまま 1e-7 → seam の残差合算・状態同期が正しい (解析解との差 4e-2 は 1 次風上の数値拡散で両者同一) | active |
| `run_0067_passiveA_tracer_seam_1st` / `run_0068_passiveA_tracer_seam_s3` | node 周期 (33³ 箱) | **受動種経路 (`passiveScalarScheme 1`) の周期 seam 試験** (plan species-passive-scalar-unification §6-2, Phase A 2026-09-16): `run_0065` 入力 (KEEP RK4, 一様流 u=10) を scheme 1 で (`0067` = 1 次風上; `0068` = SLAU `convMethod 1` + `limiter 2` + `speciesFaceReconstruction 2` = 受動種 S3 [ψ_P 無次元化 Venkat]). 判定 `run_0067/analyze_seam2.py` (周期対 ΔXi・∫ρXi dV・解析平行移動・run_0065 との差) | 周期対 (1089 対 × 3 方向) の max\|ΔXi\| = **0** (両方)、∫ρXi dV 相対変化 3.7e-8 (1 次) / 3.1e-7 (S3)、0≤roXi/ρ≤1 (S3 max 0.99806)、floor 補正 0 (両方)、流れは一様のまま (Ux 10.000000, P 101324.99)。1 次は legacy `run_0065` と max\|ΔXi\| 2.4e-7 (float 床)。S3 は解析解との差 8.9e-3 (1 次 4.0e-2 の 1/4.5; seam 近傍 1.4e-3 = 内部より小) → seam の勾配除外+gather・面値の周期整合が正しい | active |
| `run_0069_passiveA_cond_seam_1st` / `run_0070_passiveA_cond_seam_s3` | node 周期 | **凝縮ソース seam 試験 (`run_0064` プロトコル) を受動種経路で** (scheme 1; `0069` = 1 次, `0070` = S3 [`convMethod 1`, `limiter 2`, SFR 2]); 1 step, `Q0·V/dt` の seam/内部比 | 比 = **1.000000** (面 0.9999931–1.0000056 / 辺 0.9999959–1.0000037 / 角 0.9999976–1.0000015; `run_0064` の再判定は 0.9999974–1.0000031) → 受動種経路でも seam の積分ソース・輸送対角 gather は内部と一致。roQ0/rog ≥ 0、floor 補正 0 | active |
| `run_0071_passiveA_diff_gauss33_D0` / `run_0072_passiveA_diff_gauss33_D` / `run_0073_passiveA_diff_gauss64_D0` / `run_0074_passiveA_diff_gauss64_D` (2026-09-17, `forge-bin-passive-A`) | node 周期 (33³ / 64³ [`run_0040` の `Taylor-Green_64.h5`]) | **トレーサ拡散の解析解試験 (plan species-passive-scalar-unification §6-3, F-sp1)**: `run_0068` 入力 (SLAU convMethod 1 + limiter 2 + SFR 2 = 受動種 S3, RK4 unsteady, 一様流 u=10, CPG `nSpecies==1` + `tracer`), ξ₀ = exp(−(x−x₀)²/2s₀²), s₀=0.6, x₀=π−u·t_end (最終中心 x=π), t_end 0.085 s (33³: dt 1.7e-4 × 500; 64³: dt 8.5e-5 × 1000)。D0 = viscMethod 0 (拡散なし基準); D = viscMethod 1 (Sutherland μ(294.2 K)=1.8193e-5) + `Sc 1.0337e-5` → D = μ/(ρSc) = 1.4667 m²/s (2Dt = 0.2493)。判定 `scratchpad/diffB.py` (分散成長・解析ガウス σ²=s₀²+2Dt との誤差) | **解析解と一致・2 次収束**: 分散成長 Δvar(D)−Δvar(D0) = 0.2470 (33³, −0.93 %) / 0.2488 (64³, −0.23 %) vs 2Dt 0.2493 (誤差比 4.1 = 2 次); 解析ガウスとの差 L∞ 6.4e-3 / 1.5e-3 (比 4.2), L2 2.95e-3 / 7.4e-4 (峰 0.769 基準); 中心 3.1417 (π), ∫ξdx 1.5038 (解析 1.5040)。∫ρξ dV 保存 1.1e-7 / 1.5e-7 (D 有: floor 補正 0), D0 は 4.8e-5 / 1.5e-5 (S3 の負の undershoot を下限 0 に floor する補正 = 累積 2.1e-2 / 6.5e-3 [総量 75.8 / 588] が全て)。S3 の数値拡散 (D0 の分散成長) 0.0106 / 0.0023 (比 4.6, D_num 相当 0.062 / 0.014)。流れは一様のまま (Ux 10.00000, P 101324.99)。**注意**: `passiveDiffusion_d_wrapper` は `viscMethod == 0` で return する (化学種拡散 `speciesTransport_d.cu:742` と同じ規約) ので、定数粘性 viscMethod 0 + visc>0 ではトレーサ拡散が入らない (最初の試行 visc 1.76 は D0 とビット同一だった; 記録) | active (F-sp1 の解析解試験) |
| `run_0079_passiveB_dt_tracer_cpg_smoke` | node 周期 (33³ 箱) | **nSpecies==1 (CPG) + トレーサの dual-time smoke** (Phase B §6-6-iv): `run_0067` 入力 (KEEP + 移流ゲージ, cos ξ) を `timeIntegration 11, unsteady 1, dualTime 1, blockDPLUR 1, control 0`, dt 1e-3 (物理 CFL_c≈1.7), `nSubIterDualTime 10`, cfl_pseudo 5, `passiveScalarScheme 1` で 100 step | NaN 0、流れは一様のまま (Ux 10.00000, P 101324.99; 流れ残差は float 床 1e-12)、ξ は移流して 0≤ξ≤1 ([0.0467, 0.9533] = 1 次風上の拡散)、**トレーサの BDF 込み sub-iter 残差は毎 step 3.8–4.8 桁低下** (rms_roXi 4.4e-2 → 9.3e-7) | active (Phase B 検証) |
| `run_0080`–`run_0086_passiveB_order_*` | node 周期 (33³ 箱) | **dual-time の時間精度 (Phase B §6-6-ii)**: CPG + 化学種 [N2, O2] (`speciesDBFile` は case/16 の DB; 熱力学は CPG のまま = 受動的組成) + トレーサ, 一様流 u=10 (KEEP + 移流ゲージ), IC = ガウス ξ = exp(−(x−π)²/2·0.5²) と Y0 = 0.5+0.3·G (`TG_gauss.h5`), 最終時刻 0.2 s を dt 4e-3 / 2e-3 / 1e-3 (50/100/200 step, 物理 CFL_c 6.9/3.5/1.7), `nSubIterDualTime 20`, cfl_pseudo 5, 1 次移流 (SFR 0): `0080/0081/0082` = bdfOrder 2, `0083` = dt 2e-3 で nSub 40, `0084/0085/0086` = bdfOrder 1 (判定 `analyze_order.py` 相当: L2 差 3 水準の観測次数) | **BDF2: ξ 2.03, Y0 2.03 (\|Δt−Δt/2\| 6.8e-4 / 2.0e-4, \|Δt/2−Δt/4\| 1.7e-4 / 5.0e-5); BDF1: ξ 0.91, Y0 0.91**。sub-iter 倍増 (40 vs 20) の差 L2 1.1e-7 (ξ) / 2.6e-8 (Y0) ≤ 最小水準差の 1/1000。各物理 step の BDF 込み sub-iter 残差: rms_roXi 4.7–6.0 桁、rms_roY0 4.7–5.2 桁低下 (流れは一様のまま float 床 1e-13、Ux 10.00000 / P 101324.99)。NaN 0 | active (Phase B 検証, 時間次数) |
| `run_0087`/`run_0098_passiveC_m2_step_seam_s3_final` / `run_0088_passiveC_m2_rerun0067` / `run_0089_passiveC_m2_rerun0068` | node 周期 (33³ 箱) | **codex result M2 (周期 root のみで収支・総量を数える) + M5 (増分スケーリング θ_b) の検証**: `0087`/`0098` = `run_0068` プロトコル (SLAU S3, RK4 500 step) で **ξ のステップ (\|x\|<0.8 の seam 跨ぎ)** を移流 (`0087` = 中間ステージも収支に載せた版, `0098` = 確定ステージのみ = 最終仕様); `0088`/`0089` = `run_0067`/`0068` 入力を新バイナリで再取得 | **総量 (ログ) = root のみ積分**: `0088` 148.8301 (旧ログ 166.246 = 重複計上) vs 後処理 root-only 148.8301 (差 2.1e-7 = float 加算順)。**収支の恒等式 (`0098`)**: root-only ∫ρξ の変化 −4.5376e-1 (S3 のステップ肩で ξ>1 の候補を θ_b で縮めた分) = floor lo+hi (0) + 符号付き limCorr −4.5380e-1 → 不一致 4.1e-5 (総量比 5e-7); floor は 0 (θ_b が全て担う; θ_min 0 = plateau 上で候補が上限外), 0≤ξ≤1。`0087` (中間ステージ込み) は恒等式が 0.9 ずれる → 確定ステージのみ記録に修正。再取得: `0088` ∫ρξ 相対変化 −3.7e-8, `0089` +3.1e-7 (旧 `0067`/`0068` の 3.7e-8 / 3.1e-7 と同値; 収支 0 のまま) | active (codex M2/M5) |
| `run_0090_passiveC_m1_dt_c1_interior` / `run_0091_passiveC_m1_dt_c1_seam` / `run_0092`/`run_0093_…_c2_*` | node 周期 | **codex result M1 (周期 scalar-DPLUR の非対角 gather) の陰解法試験**: `run_0081` プロトコル (dual-time tI 11, dt 2e-3, nSub 20, ガウス ξ + ガウス組成 [N2,O2] CPG) を `speciesImplicitCoupling 1` + `passiveImplicitCoupling 1` で、ガウスを内部 (x=π) / seam (x=0) に置いて π シフト等価を比較; `0092`/`0093` = coupling 2 | **coupling 1**: 周期対 max\|Δ\| = 0 (ξ, Y0 とも両 run)、root-only ∫ρξ 相対変化 1.2e-6 / 1.4e-6・∫ρY0 −1.7e-7 / −5.2e-8、**π シフト等価 (内部 vs seam 中心) max\|Δξ\| 6.5e-7 / max\|ΔY0\| 4.5e-7 (L2 1.1e-7)** = float 床 → 合併 CV の更新作用素が内部と一致; sub-iter 低下 rms_roY0 4.6–5.2 桁, rms_roXi 4.7–6.0 桁 (両 run 同値)。**coupling 2 は step 4–5 で NaN**: 案C の EOS クロス応答 (`species_eos_cross_response`, NASA 熱力学) は thermalMethod 2 前提で、この CPG 受動組成の箱では成立しない (試験設定の限界; 周期 gather 経路自体は coupling 1 と共通の `dplurSweepOnce`) | active (codex M1) |
| `run_0094_passiveC_m3_cont200` / `run_0095_…_restart100` / `run_0096_…_restart_dtchanged` / `run_0097_…_restart_histdeleted` | node 周期 | **codex result M3 (checkpoint 有効性) の非定常試験**: `run_0081` プロトコル (dual-time, ガウス ξ + 組成, BDF2, dt 2e-3) を 200 step 連続 vs (a) `res_100` からの通常 restart, (b) dt を 1e-3 に変更した restart, (d) `/CHECKPOINT/roXiP` を削除した restart | (a) `history restored ... nHistoryValid=2`: 連続 res_200 との差 roXi 5.1e-7 / roY0 4.6e-7 (ρ, ρu は 0 = 一様流); (b) ログ `history NOT restored (physical dt differs from the checkpoint (file 0.002000 vs run 0.001000)) ... BDF1` → 全系 BDF1 再開 (連続との差 3.2e-4 / 7.2e-5 = dt 変更 + BDF1 初手の分); (d) `history NOT restored (missing dataset /CHECKPOINT/roXiP)` → BDF1 再開 (差 2.9e-4 / 6.6e-5)。layout に dt/bdfOrder/passiveScalarScheme/speciesImplicitCoupling を含め、scheme 0 では受動種履歴を書かない | active (codex M3) |
| `run_0099`–`run_0104_passiveD_tp_*` (SLAU nStepInner 1: step 3–13 で NaN, 不採用) / **`run_0105_passiveD_tp_c2_u0_interior` / `run_0106_…_c2_u0_seam` / `run_0107`/`0108_…_c1_u0_*` / `run_0109_…_c2_u10_interior` / `run_0110_…_c2_u10_seam`** / `run_0111`–`0114_…_u0_suth_*` | node 周期 (33³ 箱) | **plan §5.1 #17: TP 周期箱の coupling 2 試験** — `thermalMethod 2`, species [N2, H2O] (case/16 の DB), `thermoHrefTemp 298.15`, 一様 P=101325 Pa / T=298.15 K (roe = −P + KE), 組成ガウス Y_H2O = 0.1+0.1·G (幅 0.5) を内部 (x=π) / seam (x=0) に置く (ΣY=1, ρ は組成で変わる), トレーサ ξ=G, SLAU 1 次, dual-time (tI 11, dt 1e-3, nSub 20, cfl_pseudo 2, **nStepInner 4**: 1 sweep では周期箱で発散), `passiveImplicitCoupling 1`, `speciesImplicitCoupling 2` (対照 1), 100 step。`0105`–`0110`: `viscMethod 0` + `visc 1.0` (= 化学種拡散は viscMethod≠0 でのみ有効なので拡散なし); `0111`–`0114`: `viscMethod 1` (Sutherland) dt 5e-3 | NaN 0 (`0105`–`0114`)。**u=10, coupling 2 (`0109`/`0110`; EOS クロス項は mdot≠0 で作動)**: 周期対 max\|Δ\| = 0 (ρ, ρY0, ρY1, ρξ), root-only ∫ρY 相対変化 ≤6e-7・∫ρ ≤1.6e-8, **π シフト等価 Y1/ρ 6.1e-7 (L2 8.5e-8), ρ 8.3e-7, ξ 9.7e-6** = float 床、sub-iter 低下 rms_ro 3.0 / rms_roY1 3.9–4.4 / rms_roXi 3.7–4.8 桁、流れは一様のまま (Ux 9.9997–10.0002, T 298.149–298.150)。u=0 (`0105`–`0108`): 移流も拡散も無く組成は不変 (Y1 0.10000–0.20000)、周期対 0、等価 Y1 6.0e-7/7.0e-7・ρ 6e-7/3.6e-7 (ξ 2e-5/4e-5 = φ_N δρ 項の float 累積)、残差は床で桁低下なし = **EOS クロス項は u=0 では作動しない** (mdot=0)。Sutherland u=0 (`0111`–`0114`, 0.5 s): 拡散長 1.5e-3 ≪ 幅 0.5 で組成変化は 1e-5 級、代わりに float 音響ノイズ (Ux ±0.007, ρ 3.8e-5) が両 run で異なり等価 Y1 7.0e-5 / 3.4e-5 = 流れノイズ支配 → 物理が見えない設定 (記録のみ)。結論: coupling 2 の周期整合は u=10 の系列で float 床まで確認 | active (#17) |
| --- | --- | --- | --- | --- |
| `run_0009_node_keep_pure_eul` | node (median-dual) | pure KEEP, **非粘性** | **質量厳密・運動量 ~1e-7・KE 0.3%・エントロピー ~1e-5 保存**。教科書通りの KEEP 保存性 | active ✅ |
| `run_0010_cell_pure_eul` | cell (primal hex 32768) | pure KEEP, **非粘性** | **修正後: 質量厳密・運動量 ~1e-7・KE 0.4%・エントロピー ~1e-5 保存**。node と一致 (修正前は step358 発散) | active ✅ |
| `run_0007_node_keep_pure_visc` | node (median-dual) | pure KEEP, 粘性 (Re≈160) | 質量厳密・運動量 ~1e-7 保存。**KE 物理減衰 K/K0→0.644、エントロピー +1.3e-2 (第二法則)** | active ✅ |
| `run_0008_cell_pure_visc` | cell (primal hex) | pure KEEP, 粘性 | **修正後: 運動量 ~1e-7 保存・KE 物理減衰 K/K0→0.661** (node とほぼ一致; 残差は primal/dual メッシュ差)。修正前は KE×8 スプリアス増殖 | active ✅ |
| `run_keep`, `run_keep_M0.1` | cell | 旧参照入力 (旧スキーマ config) | — | ref |
| `run_0053_cond_periodic_conservation` (ローカル, 2026-09-13) | 凝縮モーメント移流の融合 (plan condensation-float-speedup §5.1 #12b): run_0052 config + N2 CPG 凝縮 `condEquilibrium: 2` (Q0–Q2 はソース 0 の移流のみ) + 非一様 ρQ IC、ラベル merged/merged1/merged1b (融合前) と condf/condf1 (融合後), 1 step と 20 step | 1 step の Σ ρQ_n V 変化は両者 5.60e-6 で一致 (差 ≤5e-10 = 反復ノイズ); 20 step は両者 +1.1e-3 (既存挙動) | active (回帰) |
| `run_0115_passiveE_m3_lim_s3_interior` / `run_0116_passiveE_m3_lim_s3_seam` / `run_0117_…_interior_preM3` / `run_0118_…_seam_preM3` | node 周期 | **codex result-2 M3 (周期 node のリミッタを group で統合, plan §4.8) の検証**: `run_0090`/`0091` プロトコル (dual-time tI 11, dt 2e-3, nSub 20, ガウス ξ + ガウス組成 [N2,O2] CPG, coupling 1) を SLAU `convMethod 1` + `limiter 2` + `speciesFaceReconstruction 2` (受動種 S3 + 2 次流れ) で、`0115`/`0116` = 2 段リミッタ (極値 group max/min → ψ group min) の新バイナリ、`0117`/`0118` = 修正前スナップショット (`forge-bin-passive-preM3`) の対照。評価 `run_0115/analyze_m3_periodic.py` | **周期対 max\|Δ\| (`limiter_Xi`)**: 修正前 2.5e-2 (内部) / **1.07e-1** (seam; 群 [285,316,347,378] = 1.0/0.893/0.893/1.0) → 修正後 **0 (厳密)**; ξ/Y0/ρ/P/U と `limiter_ro`/`Ux`/`P` の周期対差は前後とも 0 (状態 mirror 済; ρ 一様なので流れ側 ψ は自明に 1)。**π シフト等価 (内部 vs seam)**: `limiter_Xi` max\|Δ\| 1.04e-1 (L2 1.8e-2) → **7.4e-6 (L2 9.8e-7)**; ξ 6.3e-7 → 6.0e-7・Y0 4.8e-7 (float 床, 変化なし)。root-only ∫ρξ dV 相対変化 +1.16e-6 / +1.43e-6 (前後同値)、∫ρY0 −1.7e-7 / −5.2e-8、NaN なし。`check_convergence` は dual-time で従来どおり `rms_ro` 1e-13 桁 (ρ 一様=機械零) の RISING / `rms_roXi` 1.1 桁 STALLED (前後同値; 定常判定は非適用)。壁時間 57 ms/step (修正前 56–60; 2 段化で遅くならず)。`0115`/`0116` は最終バイナリ (NaN skip 付き atomic) で再取得済み。`limiter_Y0` は `output_cellValNames` に無く出力不可 (化学種 ψ の合併は単体試験 `tests/unit/test_periodic_limiter.cu` の非スケール版で担保) | active (codex result-2 M3) |
| `run_0119_passiveE_m3_lim_s3_um10_interior` / `run_0120_…_um10_seam` / `run_0121`–`run_0126_…_robump*` | node 周期 | **§4.8 の追加検証 (codex plan レビュー 3)**: `0119`/`0120` = `0115`/`0116` と同一で **u = −10** (IC `TG_gauss_um10*.h5`: roUx=−10ρ, roe は P 不変で補正, `uRef [−10,0,0]`; seam を逆向きに横切る)。`0121`/`0122` (+20 % ρ バンプ, ξ と同中心, P 一様) / `0124`/`0125` (+2 %) = 流れ側 ψ_ρ を非自明にする試み; `0123`/`0126` = 同入力の修正前バイナリ対照 | **u = −10**: 周期対 max\|Δ\| = **0 (厳密)** for ξ, Y0, ρ, P, T, U と `limiter_Xi`/`limiter_ro`/`limiter_Ux`/`limiter_Uy`/`limiter_Uz`/`limiter_P` (両 run); π シフト等価 ξ 5.7e-7 / Y0 4.8e-7 (float 床), `limiter_Xi` 8.5e-6 (L2 8.9e-7); root-only ∫ρξ 相対変化 +1.12e-6 / +1.33e-6, ∫ρY0 −2.0e-7 / −4.6e-8; NaN なし (ρ 一様のため流れ側 ψ は 1 で自明)。**ρ バンプは +20 % も +2 % も step 5–6 で `ro` NaN (箱の角ノードから; step 0 の擬似反復で rms_roe が 0.12→0.18 と増大)**、修正前バイナリ (`0123`/`0126`) も同 step で同一に発散 → この dual-time coupling-1 プロトコルは非一様 ρ の IC を受け付けない (§4.8 とは無関係; ψ_ρ の非自明な周期 CFD 検証は別プロトコルが要る = 残件)。`0121`–`0126` は `res_nan_*.h5` のみ | `0119`/`0120` active, `0121`–`0126` 破棄予定 (発散記録) |
| `run_0127_passiveF_fct_step_seam` / `run_0128_passiveF_nofct_step_seam` / `run_0129_passiveF_fct_step_interior` (旧ビルド 11:55) → **`run_0136_passiveF_fct_step_seam_v2` / `run_0137_passiveF_nofct_step_seam_v2` / `run_0138_passiveF_fct_step_interior_v2`** (現行バイナリ 208be829) | node 周期 (33³ 箱) | **plan species-passive-scalar-unification §4.7 / §6-2: dual-time seam ステップ + 保存的 FCT** (2026-09-16): `run_0091` プロトコル (dual-time tI 11, dt 2e-3, nSub 20, coupling 1/1, ガウス組成 [N2,O2] CPG) に `run_0098` のステップ ξ (\|x\|<0.8, seam 跨ぎ; IC 合成 `run_0127/gen_step_fct_ic.py`, interior 版は π シフト) を載せ、**`solver: SLAU`** + `convMethod 1` + `speciesFaceReconstruction 2` (受動種 S3) + **`limiter 1` (Barth)**、160 step (= π 移動)。`0128`/`0137` = `passiveFct 0`。評価 `run_0127/analyze_fct_step.py`。**注**: (a) `run_0115`–`0120` から複製すると `solver: KEEP` のままで **受動種 S3 も FCT も無効** (`passiveFctActive` は `cfg.solver` SLAU を要求; log に S3 有効/無効の表示なし) — 最初の試行はこれで FCT 無作動だった; (b) 指定の `limiter 2` (Venkat) は **SLAU `convMethod 1/2` + node 周期 + dual-time BDF2 で一様流が丸めから指数成長し t≈0.22 s (step 115–145) で NaN** (流れのみでも再現、cfl_pseudo/nStepInner/lowMachPrecond/SLAU2/dt に依らず; `limiter 0` も同じ、`convMethod 0` と Barth は安定; 切り分けは `_probe_passiveF/`) → Barth に変更 | NaN 0 (全 run)。**FCT (0136)**: root-only ∫ρξ dV 相対変化 **+1.37e-6** (旧ビルド 0127 +1.36e-6; `run_0091` 同プロトコルの float 床 1.2–1.4e-6 と同水準だが厳密ゲート 1e-6 は僅かに超過)、no-FCT (0137) −1.16e-6; 0≤ξ≤1: FCT [5.6e-10, 0.99324], no-FCT [0, 0.99839]; 周期対 max\|Δ\| = **0** (ξ, ρ, P, U, `limiter_Xi`); π シフト等価 (0138 vs 0136) ξ 2.9e-6 (L2 3.8e-7), Y0 1.2e-7, ρ/P/U 0, `limiter_Xi` 5.2e-4; floor/limCorr 0 (Barth ψ_P はステップで無作動)。fctCorr: dropped antidiffusion 0.283 (rel 3.4e-3, 4.6M 面), baseViol 0, pinCorr/bnd 0, **remSigned 1.13e-4 (=∫ 変化) remAbs 5.7e-4 (rel 6.8e-6)**, qL 5 sweep rel-res 5e-7, HO 残差 4e-8, 「tracer upper-bound margin violation」累積 9.5e-6 (30642 cells)。`check_passive_budget.py` は **`FAIL(NO_FCT_RECORD)`** = ツールの regex が `fctCorr %-8s` の空白詰めに合わず (要 `\s+`) → 手計算では remAbs rel 6.8e-6 で SUM>tol。壁時間 FCT 61.8–65.1 / no-FCT 63.3–67.5 ms/step (反復 2 本ずつ; 差はノイズ内 <5 %)。同一設定の反復差 roXi 3.7e-6 (`0130` vs `0132` res_100; node FCT 経路の非決定性床) | `0136`–`0138` active; `0127`–`0129` 破棄予定 (旧ビルド, 数値は 0136 系と一致) |
| `run_0130_passiveF_fct_ckpt100` / `run_0131_passiveF_fct_restart100` / `run_0132_passiveF_fct_cont200` / `run_0133_passiveF_fct_restart_from_nofct` (旧 0128 の res) → `run_0139_passiveF_fct_restart_from_nofct_v2` | node 周期 | **§4.4/§4.7 流束形履歴 (G/H/mEff/Hsrc) の restart**: `run_0136` 設定で 100 step checkpoint → `res_100.h5` (`/CHECKPOINT/roXi_fctG,_fctH,_fctHsrc,passive_fctMeff` あり) から 100 step 再開 vs 連続 200 step; `0133`/`0139` = `passiveFct 0` の res (G/H 無し) から `passiveFct 1` で再開 | **バグ: FCT 有効 run の checkpoint は復元されない** — `0131` log `history NOT restored (layout mismatch (file '...;passiveFct=1' vs run '...;passiveFct=0'))`: 読み側 (main.cpp の layout 文字列) が `passiveFctActive(cfg)` を `g_Pface_dev` 確保前に評価して 0 になる (書き側は 1)。結果 BDF1 再開で連続との差 roXi **8.4e-4** (期待 ~1e-6; 「flux-form history restored」行は出ない)。`0139` (G/H 無し checkpoint) は layout が偶然一致 (`passiveFct=0` 同士) して `history restored (10 levels)` + `[passiveFct] WARNING: flux-form history (G/H) missing for a BDF2 step: using the local form ...` → 完走 (NaN 0, ξ∈[1.4e-8, 0.986])、ただし「history missing → BDF1 再開」ではなく BDF2 + 局所形 H で続行 (remAbs rel 1.1e-2 に跳ねる = 制限不能分)。`0130` vs `0132` の res_100 (同一設定) 差 roXi 3.7e-6 | `0130`–`0132`, `0139` active (バグ記録); `0133` 破棄予定 (旧ビルド res から) |
| `run_0134_passiveF_fct_step_seam_venkat100` / `run_0135_passiveF_nofct_step_seam_venkat100` + `_probe_passiveF/*` | node 周期 | **補足: 指定どおり `limiter 2` (Venkat) で 100 step (0.2 s, 流れ発散の手前)** FCT 1 / 0; `_probe_passiveF/` = 上記 SLAU 発散の切り分け (p*/q*/r*: cfl_pseudo, SFR, coupling, convMethod, KEEP, dt, lowMachPrecond, SLAU2, limiter 0/1, conv 2, 流れのみ) と壁時間反復 (t_*) | Venkat ψ_P ではステップ肩で **θ_b (limCorr) が毎 sub-iter 作動: cumulative abs 6.5 / signed −5.0 (rel 7.8e-2, thetaMin 0)** → root-only ∫ρξ **−0.69 %** (83.717→83.140; no-FCT −0.70 %) = FCT を掛けても sub-iter 内の非保存縮小が残る (§4.7 の「limCorr ≤1e-6」ゲート FAIL; Barth では 0)。流れは 100 step まで一様 (Ux 偏差 4e-6) | `0134`/`0135` active (記録); `_probe_passiveF/` 破棄予定 |
| **`run_0156_passiveG_fct_step_seam` / `run_0157_passiveG_nofct_step_seam` / `run_0158_passiveG_fct_step_interior` / `run_0159_passiveG_fct_step_seam_venkat`** (同設定の先行 `run_0140`–`0143` [20:34 ビルド] / `run_0148`–`0151` [20:46 ビルド] は場が反復ノイズ内 [roXi max\|Δ\| ≤3e-6] で同一だが、ゲート更新 [initial-total 記録・nominal time] の要件を満たさないので破棄予定) | node 周期 (33³ 箱) | **検証 2 巡目 (§6-2, 2026-09-16; バイナリ 20:58:22 = commit b437a45a; ゲートは 22b7f0ba)**: `run_0136`–`0138` プロトコルの再取得 (dual-time tI 11, dt 2e-3, nSub 20, SLAU convMethod 1, Barth `limiter 1`, SFR 2, coupling 1/1, ステップ ξ `gen_step_fct_ic.py`, 160 step)。`0156` FCT (既定 1) / `0157` `passiveFct 0` / `0158` interior (π シフト IC) / `0159` = `0156` の **`limiter 2` (Venkat) 100 step** (1 巡目 `0134` は θ_b で −0.69 %) | NaN 0 (4 run, `nonfinite 0`), WARNING 0。root-only ∫ρξ dV 相対変化: **FCT `0156` +1.368e-6** (1 巡目 `0136` +1.37e-6 と同値 = Barth では θ_b は元々無作動 [limCorr 0]), no-FCT `0157` −1.156e-6, interior `0158` +1.134e-6, **Venkat `0159` +1.888e-6 (θ_b off で limCorr 0; 1 巡目の −0.69 % は解消)**。0≤ξ≤1 (`check_passive_field.py` 4 run **PASS**): `0156` [5.6e-10, 0.99324], `0157` [0, 0.99839], `0159` [0, 0.9999949]。周期対 max\|Δ\| = **0** (ξ, roXi, Y0, ρ, P, Ux, `limiter_Xi`; 全 run)。π シフト等価 (`0158` vs `0156`) ξ 3.0e-6 (L2 3.7e-7), Y0 1.2e-7, ρ/P/U 0, `limiter_Xi` 5.2e-4。`check_passive_budget.py` (22b7f0ba): **全 4 run FAIL** — `0156` SUM 7.30e-6 = dropped rel 3.56e-3 / base 0 / pin 0 / upper 1.14e-7 / **remainder (\|H_rem\| 累積) 7.19e-6** / 境界 0 / closure 7.9e-10 / total-vs-increment 3.0e-9 / qL 1.0e-6 / HO 4.5e-7; `0158` 7.36e-6 (remainder 7.25e-6); `0159` 5.07e-6 (floor 3.0e-10, base 3.0e-10, upper 3.2e-7, remainder 4.75e-6); `0157` (conservative モード) NOT_CONSERVED(−1.16e-6) = 総量ドリフト自体が 1e-6 を僅かに超える。remAbs (6.0e-4 絶対) は remSigned (1.15e-4 = 総量変化) の 5 倍で no-FCT の総量ドリフト (−1.2e-6) と同桁 → float32 の H=V/Δt(q^{n+1}−q^n)−ΣG 恒等式の丸め床 (1 巡目 `run_0091` の床 1.2–1.4e-6 と同水準) で、**ゲート 1e-6 は float の累積丸めでは届かない** (要判断: 許容 1e-5 か double 集計)。壁時間 FCT 61.5 / 62.9 / 57.6 / 62.7 ms/step (`0156`/`0158`/`0160`/`0162`) vs no-FCT 57.0 (`0157`) → FCT オーバーヘッド ≈ +8 % (反復ノイズ ±5 %) | `0156`–`0159` active; `0140`–`0143`, `0148`–`0151` 破棄予定 (旧ビルド写し) |
| `run_0160_passiveG_fct_ckpt100` / `run_0161_passiveG_fct_restart100` / `run_0162_passiveG_fct_cont200` / `run_0163_passiveG_fct_restart_from_nofct` (先行 `run_0144`–`0147`, `0152`–`0155` = 旧ビルドの同設定, 結果同一, 破棄予定) | node 周期 | **検証 2 巡目 (§4.7 流束形履歴の restart)**: `0156` 設定で 100 step checkpoint (`res_100.h5` に `roXi_fctG/_fctH/_fctHsrc`, `passive_fctMeff` あり, layout `passiveFct=1`) → `0161` はそこから 100 step 再開 vs `0162` 連続 200 step; `0163` = `0157` (`passiveFct 0`) の `res_160.h5` から `passiveFct 1` で再開 | NaN 0 (4 run)。**`0161` は依然として復元されない**: log `[dual-time] passive FCT flux-form history missing or inconsistent (size mismatch): all systems start from P=PP=current, first physical step is BDF1` (「(G/H/mEff) restored」行は出ない) → 連続 `0162` res_200 との差 **max\|ΔroXi\| 8.35e-4** (L2 3.6e-4; 期待 ~1e-6; 同設定反復ノイズ `0160` vs `0162` res_100 3.1e-6)、ρ/roY0 は 0。**真因 (2 巡目で残った同型バグ)**: `solver_density_cuda/cuda_forge/speciesTransport_d.cu` `passiveFctHistoryFromHost()` 先頭の `if (!passiveFctActive(cfg)) return false;` — `passiveFctActive` は `g_Pface_dev != nullptr` を要求するが、`initDualTimeHistory` (main.cpp:1230) の時点では `passive_Pface_alloc` 未実行で nullptr → false を返し、読み側 (main.cpp:1050) は `miss` が空なので「size mismatch」と表示する。1 巡目の修正は layout 文字列を `passiveFctConfigured` に変えただけで、復元関数の判定は `passiveFctActive` のまま (要 `passiveFctConfigured` 化 or 復元前に `fctAlloc`/Pface 確保)。`0163`: log `history NOT restored (layout mismatch (file '...passiveFct=0' vs run '...passiveFct=1'))` → 全系 BDF1 で再開、完走 (ξ∈[1.5e-8, 0.98543]; 「履歴欠落 → BDF1」の契約は満たすが文言は layout mismatch で「FCT history missing」ではない)。budget: `0160` FAIL (remainder 4.6e-6), `0161` FAIL (remainder 4.3e-6 + **field ゲートが restart run の totalTime 0.4 ≠ nominal nStepOuter×dt 0.2 で落ちる = ツールが restart 継続を想定していない**), `0163` FAIL (remainder 4.3e-6), `0162` FAIL (remainder 8.9e-6) — いずれも上記 float 床; field ゲート (bounds) は `0160`/`0162`/`0163` PASS | `0160`–`0163` active (バグ記録); `0144`–`0147`, `0152`–`0155` 破棄予定 |
| `run_0164_passiveG_fct_restart100_fixed` | node 周期 | **restart 試験 B (修正後)**: `run_0160` `res_100.h5` (FCT 流束形履歴付き) から `passiveFctHistoryFromHost` を `passiveFctConfigured` で gate した binary (commit 8d2661b1) で 100 step 再開、`run_0162` 連続 200 step と比較 | log `[dual-time] passive FCT flux-form history (G/H/mEff) restored` + `history restored … nHistoryValid=2`。連続 run との差 **max\|ΔroXi\| 3.2e-6** (L2 5.8e-5; 修正前 `0161` 8.35e-4)、ρ/roU/roY/roe は **0** (ビット一致)。`check_passive_budget.py` **PASS** (restart 開始時刻 0.2 を IC の checkpoint から取得, sum 0, remainder 4.3e-6 ≤ 6e-8×100…) | active (§6-6 (iv) restart) |

### KEEP 陰解法 (block-DPLUR dual-time) 物理CFL掃引 (2026-06-29)

KEEP を陰解法化するときの flux Jacobian 方針 ([`plans/active/convection-keep-revive-node.md`](../../plans/active/convection-keep-revive-node.md) §7)
の検証。**新規コードゼロ** (`solver: KEEP` + `dualTime: 1` + `blockDPLUR: 1`, `timeIntegration: 11`) で、LHS は
既存の Roe 分割ヤコビアン (`accumulate_split_jacobian_cf`) を流用。全 cell・非粘性・共通 T≈3.5。explicit は
物理CFL≈0.05 (dt=0.007)、implicit は dt を上げて物理CFLを掃引。掃引解析は
[`analyze_implicit_sweep.py`](analyze_implicit_sweep.py) → [`implicit_sweep_ke_entropy.png`](implicit_sweep_ke_entropy.png)。

| run_* | 物理CFL | 主要設定差分 | K/K0 (T末, 理想=1) | dS/\|S0\| | 状態 |
| --- | --- | --- | --- | --- | --- |
| `run_0011_cell_keep_expl_ref` | 0.05 | explicit RK4 (基準) | 1.0033 | -1.1e-6 | active ✅ |
| `run_0012_cell_keep_impl_cfl005` | 0.05 | implicit, dt=explicitと一致 | **1.0033 (explicitと一致)** | -5.7e-7 | active ✅ |
| `run_0013_cell_keep_impl_cfl02` | 0.2 | implicit dt=0.028 | 1.0032 | -8.0e-7 | active ✅ |
| `run_0014_cell_keep_impl_cfl05` | 0.5 | implicit dt=0.07 | 1.0032 | 1.8e-6 | active ✅ |
| `run_0015_cell_keep_impl_cfl1` | 1.0 | implicit dt=0.14 | 1.0031 | 1.1e-5 | active ✅ |
| `run_0016_cell_keep_impl_cfl2` | 2.0 | implicit dt=0.28 | 1.0020 | 4.7e-5 | active ✅ |
| `run_0017_cell_keep_impl_cfl4` | 4.0 | implicit dt=0.56 | 0.9983 | 1.8e-4 | active ✅ |
| `run_0018_cell_keep_impl_cfl8` | 8.0 | implicit dt=1.12 | 0.9865 | 6.1e-4 | active ✅ |
| `run_0019_cell_keep_impl_cfl16` | 16 | implicit dt=2.24 | 0.9542 | 1.8e-3 | active ✅ |
| `run_0020_cell_keep_impl_cfl2_sub50` | 2.0 | CFL2 で nSubIterDualTime 20→50 | 1.0020 (sub20と完全一致) | 4.7e-5 | active ✅ |
| `run_0021_cell_keep_impl_cfl16_sub80` | 16 | CFL16 で nSubIterDualTime 20→80 | 0.9536 (sub20とほぼ同一) | 1.8e-3 | active ✅ |

### keepDissType (ES 散逸レイヤ) L2 較正 run ([plan](../../plans/accepted/convection-keep-es-dissipation.md))

| `run_*` | 目的・設定差分 | KE drop (500step, M0.4, RK4) | 状態 |
| --- | --- | --- | --- |
| `run_0022_cell_keep_dissoff_l2` | keepDissType=0 (新バイナリで無散逸経路の不変確認) | **−0.33%** (run_0011 参照と同挙動, 差は atomicAdd ノイズ級) | ref (不変確認) |
| `run_0023_cell_keep_diss015_l2` | keepDissType=1, σ=0.15, lowMachPrecond=1 | **8.37%** (市松6桁減衰の対価) | ref (σ 較正上限) |
| `run_0024_cell_keep_diss005_l2` | 同 σ=0.05 | **2.71%** (市松~4桁減衰で十分 → **既定採用**) | ref (σ 較正・既定根拠) |
| `run_0025_cell_keep_dissmat005_l2` | **keepDissType=2 (matrix ES)**, σ=0.05 | **1.36%** (scalar の半分。市松減衰は同等以上=case/35 run_0018) — せん断/エントロピー波に \|Un\| のみの選択性が効く | ref (Step 2 matrix L2) |
| `run_0026_cell_keep_dissmat005fullc_l2` | matrix σ=0.05 **フル c** (c' 無し) | **4.35%** (同一 σ 比較では c' 比 3.2 倍) — ただし下記 σ 掃引で「σ 弱化すれば同等 Pareto 点に到達」と判明 (同一 σ 比較は設計比較としては誤導) | ref (c' 切り分け・同一σ) |
| `run_0027_cell_keep_dissmat_fullc_s0015_l2` | フル c **σ=0.015** (σ 弱化で c' 代替) | **1.10%** (市松 3.9e-8 = case/35 run_0021)。**単一マッハ領域では c' 不要・σ 弱化で十分**。c' の残存価値はマッハ混在流の面ごと自動調整のみ | ref (σ掃引・c'代替検証) |
| `run_0028_cell_keep_dissmat_fullc_s004_l2` | フル c σ=0.04 | 3.43% (σ にほぼ線形) | ref (σ掃引) |

**所見**: (1) KEEP+block-DPLUR は **CFL≈16 (explicit の320倍) まで NaN なしで完走** — Roe-Jacobian LHS の流用で
陰解法が安定動作。(2) 物理CFL≤2 では explicit と KE/エントロピー保存が一致 (matched CFL.05 は K/K0=1.0033 で
explicit と4桁一致) → **LHS の Roe 散逸は収束解を汚染しない**。(3) CFL を上げると KE が単調減衰 (CFL16 で
-4.6%)・エントロピー増。**ただしこの劣化は内部反復不足 (LHS散逸漏れ) ではなく BDF2 dual-time の時間離散誤差**が
支配: CFL2 で内部残差は 6.8e-5→3e-9 と完全収束し sub20≡sub50、CFL16 でも sub20≈sub80 で KE が回復しない
(サブ反復を増やしても改善しない=温度離散誤差は irreducible)。→ **LHS の Jacobian 選択は解の質にほぼ無関係**で、
大 dt での精度限界は外側の時間積分 (BDF2) が決める。LES では物理 dt を時間精度が許す範囲に保つこと。

成果物: KE・エントロピー時間履歴 [`ke_entropy_history.png`](ke_entropy_history.png)、ポストスクリプト
[`plot_ke_entropy_history.py`](plot_ke_entropy_history.py) (旧 `plot_taylorGreen.py` の現行版=res ベース・多重度補正付き)。
各 run の VERDICT は `CONVERGENCE_VERDICT.txt` (非定常 TGV のため定常収束ツールは `NOT CONVERGED`=正常)。

### 結論

- **pure KEEP は cell・node とも KE/エントロピー/運動量を保存する** (修正後)。非粘性で KE 0.3–0.4%・エントロピー
  ~1e-5・運動量 ~1e-7、粘性で KE 物理減衰・エントロピー物理増加。cell と node が一致 (差は primal hex vs median-dual)。
- **根本原因 (修正済)**: cell 全周期で運動量が線形注入され seam が非周期化していたのは、**`bint_d["partnerCellID"]`
  (device) が host から未転送**だったため (`setPeriodicPartner` は host のみ充填、device コピーは yaml uniform 経路だけ)。
  periodic_d が未初期化 device partnerCellID を読み ghost が誤値に。`setPeriodicPartner` 直後に H2D コピーで根治。
  node は host の partnerCellID を直接使う (`buildPeriodicNodeGroups`) ため無傷だった。float/double・ghost 更新
  タイミングは無関係 (切り分け済)。
- 切り分けの決め手: ghost セル値 vs 解析 TGV(ghost 重心) が修正前 6112/6144 不一致 (誤差 0.43) → 修正後 0 (1.5e-7)。
- 補足: 過去 cell で「KEEP」と思っていた流束は legacy で `if(false)` 無効化されており実体は **MUSCL+Roe**。本物の
  KEEP 中心流束は `Revive KEEP` (79b4e67) で初めて有効化された (この periodic バグが顕在化したのもそのため)。
- 旧 `run_keep` が現バイナリで回らないのは solver 退行ではなく①config スキーマ進化②圧力フロア (`pMin` 既定 1.0 Pa
  が無次元低圧場をクランプ) であり、modern config + `pMin` 引き下げで解消する。

### node L2 / L3 (WALE 込み粘性 TGV Re=1600) — keepDiss 検証続き

| `run_*` | 目的・設定差分 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_0029_node_keep_dissmat005_l2` | **node L2**: node TGV M0.4 + matrix σ=0.05 (run_0009 σ=0 と比較) | KE cost **+1.371%**/500step = cell (+1.356%) と一致。**散逸レイヤは node 内部で cell 同等** | ref (node L2 PASS) |
| `run_0030_cell_keep_wale_re1600_s0` | **L3 基準**: 粘性 Re=1600 (visc 2.5e-4, Pr0.71) + WALE, keepDissType=0, 3600step (t*≈10) | ε* ピーク 0.066 @ t*=10.1、K/K0(10)=0.71。層流期はほぼ無散逸 (正)。32³ ゆえピーク過小・遅め | ref (L3 基準) |
| `run_0031_cell_keep_wale_re1600_diss` | L3 + matrix **σ=0.05** | ピーク 0.083 @ t*=8.1、K/K0(10)=0.49。**層流期 (t*<4) に 5-6% 食う = 滑らか場への連続ドレインが強すぎ**。LES 用には過大 | ref (L3 σ過大の記録) |
| `run_0032_cell_keep_wale_re1600_diss0015` | L3 + matrix **σ=0.015** | ピーク 0.080 @ t*=8.4、K/K0(10)=0.60、層流期損失 2.6% | ref (L3 σ掛引) |
| `run_0033_cell_keep_wale_re1600_diss002` | L3 + matrix **σ=0.02** | ピーク 0.082 @ t*=8.4、K/K0(10)=**0.594**、層流期 3.1% | ref (L3 σ掛引) |
| `run_0034_cell_keep_wale_re1600_diss003` | L3 + matrix **σ=0.03** | ピーク 0.083 @ t*=8.4、K/K0(10)=**0.557 (DNS帯 0.50-0.57 内)**、層流期 4.1% | ref (L3 σ掛引) |

L3 の判定は定性ゲート (ピークが t*≈9±1・カーブ滑らか・NaN なし): 全 run 通過。32³ では層流期と終値を同時に満たす σ は無い (解像度律速)。~~DNS 参考値 (K/K0(10)≈0.50-0.57...) はアンカー2点+積分近似で ±0.05 幅~~ **← この近似帯は誤り** (下の 64³ 節の実データ参照: 真値 K/K0(10)=0.596。32³ σ=0.03 を「帯内」とした判定は近似帯の誤りによるもので、実データでは σ=0.02 の 0.594 がほぼ厳密一致だった)。

### 64³ L3 定量化 (2026-07-19, 実 DNS データで σ 確定)

DNS 参照を実データ化: [`ref_dns/TGV_Re1600.dat`](ref_dns/README.md) (**Dairay+2017 JCP, Incompact3d 512³**, t=0..20 の Ek/ε_t 表; 真のアンカー **K/K0(10)=0.596, ε* ピーク=0.1029 @ t*=8.98**)。`plot_dissipation_rate.py --dns` で重ね描き可。

| run_* | σ | K/K0(4) 誤差 | K/K0(10) 誤差 | ε* ピーク (DNS 0.103) | ピーク t* (DNS 8.98) | 状態 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0035_cell_keep_wale64_re1600_s0` | 0 | **+0.8%** | +1.7% | 0.085 | **8.96** | ref (64³ L3 基準) |
| `run_0036_cell_keep_wale64_re1600_diss002` | 0.02 | −0.6% | **−1.5%** | 0.081 | 8.68 | ref (64³ σ掛引・**推奨確定**) |
| `run_0037_cell_keep_wale64_re1600_diss003` | 0.03 | −1.3% | −2.8% | 0.079 | 8.68 | ref (64³ σ掛引) |
| `run_0038_node_keep_wale64_re1600_diss002` | 0.02 (**node**) | −0.7% | −2.1% | 0.082 | **8.96** | ref (64³ node L3) |
| `run_0039_node_keep_wale64_re1600_s0` | 0 (**node**) | +0.7% | −0.4% | 0.092 | 9.24 | ref (64³ node σ=0 基準) |
| `run_0040_node_keep_wale64_re1600_diss002_jump` | 0.02 (**node, keepDissJump=1**) | **+0.7%** | **−0.8%** | **0.088** | **8.96** | ref (**recon-jump L3**) |
| `run_0041_node_keep_wale64_re1600_diss002_jump2` | 0.02 (**node, keepDissJump=2**) | +0.7% | −1.4% | 0.089 | **8.96** | ref (**sign-property L3**) |
| `run_0042_node_keep_wale64fix_s0` | 0 (**修正後の本物 WALE**) | −0.3% | −3.4% | 0.087 | **7.84 (早すぎ)** | ref (real-WALE 基準) |
| `run_0043_node_keep_wale64fix_diss002_jump2` | 0.02 jump2 (**本物 WALE 併用**) | −0.4% | −3.6% | 0.085 | 7.84 | ref (real-WALE+ES) |
| `run_0044_node_keep_sigma64_s0` | 0 (**σ-model** LESmodel:2) | −0.5% | −5.5% | 0.082 | 7.84 | ref (σ-model 基準) |
| `run_0045_node_keep_sigma64_diss002_jump2` | 0.02 jump2 (**σ-model 併用**) | −0.5% | −5.6% | 0.082 | 7.84 | ref (σ-model+ES) |
| `run_0046_sst_hdt_uniform_ek1` / `run_0047_sst_hdt_uniform_ek0` / `run_0048_sst_hdt_sine_ek1` | **一様乱流減衰 (SST の E_t 保存試験)**: run_0007 の node 周期メッシュ流用、SLAU + SST、u=0・ρ=1.4・P=1・ω=1、k=0.1 一様 (uniform) / k=0.1(1+0.5 sin x) (sine)、RK4 explicit dt 0.005 × 400。ek1=`sstEnergyIncludesK: 1`, ek0=0 (plan [turbulence-sst-energy-includes-k](../../plans/active/turbulence-sst-energy-includes-k.md) §6.1) | **ek1: ΣV(E_m+ρk) 変化 ≤1.3e-7、max\|u\| 1.9e-7 (一様 p* で偽の力なし)、k 0.100→0.0847 の減少が c_v ΔT に 5 桁一致 (Δe=+0.01534=−Δk)。ek0: E_t −0.81 % (k の散逸が消える)、T 不変。sine: p* 勾配で u≈0.02 が立つが E_t 変化 ≤1.1e-7・質量 4e-8** | active (保存試験) |
| `run_0049_sst_hdt_uniform_ek1_dual` / `run_0050_sst_hdt_sine_ek1_dual` | 同上を **dual-time (timeIntegration 11 + dualTime 1, block-DPLUR)** で。エネルギー行の BDF に ρk を含める修正後 | E_t 変化 uniform +3.5e-7 / sine −3.1e-6 (subiter 収束で決まる)、c_p ΔT = +0.02148 (explicit と同一)。修正前は E_t −0.81 % で T 不変だった | active (保存試験) |
| `run_0051_monitor_dualtime` | run_0050 の複製 (100 step)。**console モニタ行** ([plans/accepted/architecture-runtime-monitor-line.md](../../plans/accepted/architecture-runtime-monitor-line.md)) の dual-time 陰解法 (node) 検証: `t`/`dt`/`maxCFL` (物理 CFL) + ms/step + 残差要約。`_baseline/` は旧バイナリ比較用 | 完走・NaN 0。`t/dt/maxCFL 0.13` 表示、log 4330→230 行、CSV 旧×2/新×2 列別比較 (`csv_compare.txt`) 合格 (1.2e-4 vs 床 1.3e-4) | 破棄予定 (ログ形式検証) |

**結論 (σ 確定)**: 64³ では層流期・終値とも全 σ≤0.03 で ±3% 内 = 32³ の「両立不可」は解像度律速と確定。σ=0 (+1.7%) と σ=0.02 (−1.5%) が DNS をほぼ対称に挟み、**解像 LES の推奨は σ=0.02** (市松頑健性 L1 を持ちつつ KE 追従 ≤1.6%)。σ=0.03 は 64³ では一律に劣る (旧「0.03=帯内」判定は撤回)。TGV 追従だけなら最適 σ≈0.01 だが、σ は物理較正ノブではなく市松ロバスト性の床 — 解像度を上げるほど最適値は下がる (32³: ~0.02 → 64³: ~0.01) ので、「L1 を満たす最小 σ」で選ぶ。ε* ピーク高さは 64³+WALE の解像限界で全 run 8 割どまり (LES として正常)。成果物 `dissipation_rate_L3_64.png`。

**node L3 (run_0038)**: 同条件 node (median-dual 65³=274625 CV, 周期 DOF 合併 12481=厳密) は
cell とカーブがほぼ重なり (K/K0(10) で node −2.1% vs cell −1.5%、ピーク時刻は node 8.96 が DNS
8.98 により近い)、**node×σ=0.02 も解像 LES として cell 同等に検証済**。KE 集計は周期 slave CV の
ミラー重複を除くため `plot_dissipation_rate.py --node-mesh` (半開区間 dedupe) を使うこと
(素朴な総和は周期面を 2-8 倍過重み)。成果物 `dissipation_rate_L3_64_node.png`。

**再構成ジャンプ散逸 (keepDissJump=1, run_0039/0040)**: matrix 散逸の Δw を再構成後ジャンプで
組むと ([recon-jump plan](../../plans/accepted/convection-keep-diss-recon-jump.md))、σ=0.02 の
KE コストが**終値 −2.1%→−0.8%・層流期 −0.7%→+0.7% (=σ=0 と同一・実質ゼロ)** に縮み、遷移期の
ε\* の早すぎる立ち上がりも解消、ピーク時刻 8.96 (σ=0 単独は 9.24 と遅い) は維持 — つまり
**市松頑健性と σ=0 並みの KE 追従を両立**。市松減衰は無傷 (case/35 run_0039/0040)。
**node LES の新推奨: `keepDissType: 2, keepDissCoeff: 0.02, keepDissJump: 2`**。
成果物 `dissipation_rate_L3_64_node_jump.png`。

**sign-property クリップ (keepDissJump=2, run_0041)**: 再構成ジャンプの特性射影を生ジャンプ射影と
minmod し、各波で積 ≥0 を構造保証 = **エントロピー散逸性が証明付きで復活** (TeCNO 型)。KE コストは
jump=1 (−0.8%) と生 (−2.7%) の間の **−1.4%**、ピーク時刻 8.96・層流期 +0.7% は jump=1 と同一。
市松減衰も無傷 (case/35 run_0042)。**証明と性能の両立点として jump=2 を推奨既定とする**
(jump=1 は証明なしの最軽量オプションとして残置)。成果物 `dissipation_rate_L3_64_node_jump2.png`。

**重大訂正 (2026-07-19, WALE 不活性バグ発覚)**: 本節までの全「WALE 併用」run (32³ run_0030-0034・
64³ run_0035-0041, cell/node とも) は、**WALE が不活性 (vis_turb≡0) の ILES だった**。真因は
① 壁なしメッシュで wall_dist≡0 → Ls=min(κd,CwΔ)=0、② WALE の Sd テンソルが成分 2 乗の誤式
(行列 2 乗が正; 純せん断で Sd=0 になるべき性質を破っていた)。両方修正済
([turbulence-wale-fix](../../plans/accepted/turbulence-wale-fix.md))。**数値結果は「ILES
(KEEP+分子粘性+ES散逸)」として全て有効** — ラベルのみ訂正。修正後の本物 WALE で再検証した結果
(run_0042/0043): **WALE は遷移期に先回りして散逸し (ピーク t*=7.84 ≪ DNS 8.98)、終値も −3.4% と
ILES より悪化**。64³ TGV 級の解像/遷移流では **WALE off の ILES + matrix ES (σ=0.02, jump=2) が
最良 (−1.4%, ピーク 8.96)** = **解像 LES の推奨は LESorRANS: 0 + ES 散逸**。WALE は「本当に
未解像の高 Re 乱流」用オプションとし、その適用検証は今後の課題。成果物
`dissipation_rate_L3_64_walefix.png`。

**σ-model (run_0044/0045)**: 静的モデル改良版 (LESmodel:2) も検証したが、層流期 ν_t は設計どおり
WALE 比 −25% になるものの圧縮性 TG の早期 3 次元化でゼロ性質の効く時間が短く、発達後は C_σ=1.35 の
散逸が勝って **−5.5% (WALE より悪い)・ピーク 7.84**。**64³ TGV 級では静的 SGS はどちらも逆効果、
ILES+ES (σ=0.02, jump=2) が最良のまま**。σ-model の本領 (壁乱流・回転流・未解像高 Re) の検証は今後
([plan](../../plans/accepted/turbulence-sigma-model.md))。成果物 `dissipation_rate_L3_64_sigma.png`。
| `run_0052_perf_regress_node_periodic_dualtime` (ローカル, 2026-09-12) | 高速化ブランチ `feature/perf-3d-speedup` の回帰: run_0049 (node 周期, dual-time 陰解法, SST, sstEnergyIncludesK=1) の res_400 から 100 step 継続、基準バイナリ ×2 vs 最終 (`tools/perf_regress.py`) | P/T/k/ω/ρK/ρΩ/μt の差 ≤4.3e-7 で base×base と同等。速度 Ux/Uy は base×base 自体が 0.32 (減衰乱流の run 間ばらつき) で比較不能 | active (回帰) |
