# 3D node SST 生産計算の速度向上 (FP64 律速の除去・陰解法 sweep の軽量化)

## メタ

- **area**: `architecture / performance`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/architecture/performance.md`](../../methods/architecture/performance.md) (性能プロファイルと数値精度方針の現在仕様)
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md)
  - [`procedures/development-environment.md`](../../procedures/development-environment.md) (速度評価は native)
- **related_plans**:
  - [`thermophysics-multicomponent-tpgas.md`](../accepted/thermophysics-multicomponent-tpgas.md) (M6: 輸送係数 FP32 化・cp+h 融合・Rmix キャッシュは実装済。本 plan はその後に残る double 演算が対象)
  - [`gpu-implicit-plan.md`](../accepted/gpu-implicit-plan.md) (block-DPLUR の構造)
  - [`tooling-cloud-gpu-env.md`](tooling-cloud-gpu-env.md) (AWS g5 環境)
- **created**: `2026-09-12`
- **owner**: `Claude (branch feature/perf-3d-speedup)`

## 1. 目的

case/16 の 3D node SST 多成分 TP 生産計算 (run_0234 相当, 2.37 M 節点 / 7.04 M 双対面) は A10G で
**82.85 ms/step** (12000 step ≈ 17 分) かかる。律速を計測で同定し、**収束解を変えずに** 1 step の
所要時間を段階的に下げる。目標は同一 case で **≤ 40 ms/step (2 倍以上)**。さらに陰解法側
(sweep 数・CFL) の工夫で「収束までの壁時計」も短縮する。

## 2. スコープ

- **やる**: 計測基盤 (`tools/bench_steps.sh`, nsys/ncu 手順の記録)、FP64 混入の除去 (リテラル昇格・熱力学の
  面評価)、化学種拡散カーネルの再構成、リミッタの間接呼び出し除去、block-DPLUR sweep の対角キャッシュと占有率改善、
  小カーネル融合、節点の RCM 再番号付け (変換時 opt-in, 2026-09-12 に追加)。各項目は A/B 計測 + 場の一致確認を伴う。
- **やらない**: MPI 多 GPU 化、
  スキーム自体の変更 (SLAU/SST の式は不変)。倍精度が意図的に必要な箇所 (幾何前処理の閉性、双対体積・重心 float32 桁落ち対策
  [`node-yp1-dual-geometry-float32-fix`]、周期・軸対称の閉性) は触らない。

## 3. 関連 docs と前提

- 速度評価は native で行う ([`procedures/development-environment.md`](../../procedures/development-environment.md))。
  基準環境は AWS g5.xlarge (A10G, CC 8.6, FP64 の加算/乗算/FMA スループットは FP32 の 2/128 = **1/64**)。ローカル RTX 3060 も同じ CC 8.6 で FP64 比は同じ。
- 計測 run: `case/16.nozzle_wys/run_0400_perf_baseline` (AWS `~/forge-perf`)。run_0234 の入力 (メッシュ h5・config) を複製し、
  `res_12000.h5` を index コピーで IC にした収束場からの継続 (100 step / 30 step)。`res_*.h5` は書かない。
- 既存の性能上の工夫 (残差の非同期 reduction、`FORGE_KERNEL_SYNC` 既定 off、LSQ 係数事前計算、リミッタ 5 変数融合) は
  維持する。

## 4. 設計方針

### 4.1 計測結果 (2026-09-12, A10G, run_0400_perf_baseline, 100 step)

`FORGE_PROFILE=1` セクション別 (ms/step): turbulence_model (k/ω + 化学種輸送・拡散・ソース) 19.7 / convective_flux 16.5 /
time_integr (block-DPLUR 5 sweep) 14.7 / limiter 9.7 / viscous_flux 4.9 / calc_gradient 4.5 / set_dt 1.5 / gas_properties 1.0。
host 側オーバーヘッドは無視できる (合計 = 82.6/82.85 ms)。

nsys カーネル別 (30 step 平均, ms/step) と ncu (sudo 必須) の律速指標:

| カーネル | ms/step | SM% | MEM% | FP64 pipe | 占有率 | 判定 |
| --- | --- | --- | --- | --- | --- | --- |
| `SLAU_d` | 16.2 | 80 | 14 | **81 %** | 33 % (110 regs) | FP64 律速 |
| `implicit_defect_correction_block_d<float>` ×5 | 2.9 ×5 = 14.9 | 15 | 83 | 0 | 28 % (128 regs) | メモリ律速 (long-scoreboard stall 81 %) |
| `species_diffusion_d` | 13.6 | 87 | 10 | **87 %** | 33 % | FP64 律速 (面状態 Y/T/P・J_s・h_s(T_f) の演算が double。`thermo_Dbinary/Dmix` 自体は M6 で float 化済) |
| `limiter_r1_fused5_d` | 9.5 | — | — | — | — | 関数ポインタ経由の間接呼び出し (要確認) |
| `dependentVariables_d` | 6.1 | 85 | 9 | **86 %** | 32 % | FP64 律速 (TP Newton 全量 double) |
| `viscousFlux_d` | 4.7 | 78 | 31 | **78 %** | 33 % | FP64 律速 (意図的 double 無し → リテラル昇格) |
| `lsqPreGrad_internal_d` | 2.7 | — | — | — | — | gather |
| `setCFL_pln_d` | 1.2 | 86 | 53 | **86 %** | 90 % | FP64 律速 (リテラル昇格) |
| `scalar_diffusion_first_order_d` ×2 | 0.83 ×2 | 19 | 94 | 19 | 88 % | メモリ律速 (健全) |

**結論**: 1 step の約半分 (~42 ms) が FP64 パイプで消費されている。原因は 2 系統:
(a) **接尾辞なし浮動小数リテラル** (`0.5*x`, `2.0/3.0` 等) が `flow_float`=float の式を double に昇格させる
(`viscousFlux_d` に意図的 double は無いのに FP64 78 %)。(b) **熱力学 (NASA-9) の面ごと/セルごと double 評価**
(`thermo_h_mix` を SLAU の面ごと L/R で、`species_diffusion_d` は面状態の組立と h_s(T_f) 結合が double、`dependentVariables_d` の Newton)。
block-DPLUR は逆にメモリ律速で、sweep ごとに対角 5×5 と近傍幾何 (ccx/ccy/ccz) を再構築している。

### 4.2 方針

1. **数値精度の原則** (methods/architecture/performance.md に明文化): 状態・残差は float32。カーネル内の一時演算は
   float32 を既定とし、`double` は「桁落ちが実測で問題になる箇所」に限定して明示する。リテラルは `f` 接尾辞
   (または `static_cast<flow_float>`) を付ける。既存の意図的 double (閉性・幾何・Newton の初期実装) は本 plan の A/B で
   float 化の可否を検証してから置換する。
2. **熱力学の面評価を float 化 (離散式は不変)**: `thermo_*` に float 版 (`thermo_h_mass_f` 等、係数は double 保持のまま
   float へキャスト、クランプ・外挿の分岐は同一) を追加し、SLAU の `h_mix(Y, T_face)`・化学種拡散の面状態組立と
   エンタルピー結合 `Σ h_s(T_f) J_s` で使う。**評価点は現行どおり面状態 (T_f, P_f, Y_f)** — セル値を評価してから補間する案は
   `h(fT0+(1−f)T1) ≠ f h(T0)+(1−f)h(T1)` で離散式が変わるため採らない (codex plan レビュー M1)。
3. **`dependentVariables_d` の温度反転はハイブリッド** (`physProp.thermoFloat: 1`, opt-in): float Newton (warm start, 最大 12 反復)
   で ~1e-6·T まで寄せ、double の Newton を 1 段だけ当てて研磨する (`thermo_T_from_e_hybrid`)。double 評価は cph_mix 1 回
   (従来は反復数+1 回)。cp/h は研磨点の double 値から Taylor で組む。単体検証 `tools/test_thermo_float.cpp` (DB 4 種・組成 8 通り・
   50–6000 K・区間境界・datum 有無、厳密 double 参照比): **誤差 ≤1.0e-8·T** (float の T 格納 ulp 6e-8 未満) で全 PASS。
   純 float 反転は ~1e-6·T で、abs datum の H2O は 20 反復に張り付くため採らず、`thermoFloat` は `thermoHrefTemp>0` を必須にした
   (codex M2 対応)。凝縮 (二相 EOS) セルは従来経路。
4. **block-DPLUR**: sweep 0 で組んだ対角 5×5 を保存し sweep ≥1 は対角組立と近傍幾何読みを省略する。**適用は float・
   point-DPLUR 経路 (implicitSolvePrecision 0, lineImplicit 0) に限定**、キャッシュは各 `blockDPLURSolve` 呼び出しで更新、
   軸/壁/等温壁の拘束行と RHS 処理・周期ミラー・ピボット失敗処理は各 sweep で従来どおり (codex M6)。現行は
   `__launch_bounds__(128)` (min blocks 指定なし, 128 regs, 占有率 28 %) で、占有率変更はレジスタ・spill・実測で評価する。
5. **リミッタ**: `limiter_scheme` を template 引数にして関数ポインタ呼び出しを除去 (ビット同一)。
6. **小物**: `fill_limiter_d` ×5 → 1 カーネル、`gasProperties_d` の `pow(x,1.5)` を float `powf` に (`x*sqrt(x)` はビット同一でない, 別評価)、
   `convectiveFlux_d_wrapper` の毎ステップ `cudaMemcpyToSymbol` ×5 を初回だけに (ビット同一)。
7. **陰解法側** (収束までの壁時計): 上記で 1 step が軽くなった後、`nStepInner` (5→3) と `cfl_pseudo`・`implicitRelax` の
   組合せで「同じ残差到達までの壁時計」を比較する。数値解は同じ不動点なので収束解は不変、経路のみ変わる。

### 4.3 判定基準 (数値の同一性・収束・速度)

- **非決定性の前提**: 面流束は float `atomicAdd` で蓄積されるため同一バイナリでも全場のビット一致はしない。
  基準バイナリ同士 (base×base) の 100 step 継続で **run-to-run ノイズ床**を測る (2026-09-12 実測, `run_0401_perf_verify`
  `cmp_base2.txt`: 場のスケール max|a| で正規化した最大差は vis_turb 1.85e-3 / Uz 4.8e-4 / h0 1.9e-4 / 他 ≤6e-5)。
- **メッシュ品質の悪いセルでの確認 (ユーザ指摘 2026-09-12)**: float 化で高 AR・歪みセルが悪化しないかを、差 |Δ| を**壁距離ビン**
  (壁第 1 層 y₁ 0.6 µm・AR ~1200 が最小ビン) ごとに max / 99.9 パーセンタイルで集計し base×base ノイズと比較する
  (`cmp_by_walldist.py`)。実測 (`run_0401_perf_verify`, lit2 / thermof): 全ビン・全場で p99.9 がノイズの ≤1.5 倍、壁第 1 層への
  集中なし (max は単一節点の非決定性で 2〜9 倍ぶれる)。今回触っていないもの: 双対体積/重心・LSQ 擬似逆・閉性 (double のまま)、
  状態量・幾何量の格納精度 (元から float32)。
- **場の一致 (2026-09-12 改訂, codex result M3 採用)**: 各項目後に同 IC から 100 step 継続し、`VALUE/*` 全量の正規化最大差
  (`tools/perf_regress.py cmp`) を見る。float 化は熱力学の丸めを 1e-6 相対で意図的に変えるので、**基準は絶対値**: A = ρ/P/T/ρY/音速/μ ≤ 1e-5、
  B = 速度成分 (**|U| ベクトル尺度で正規化**; 成分の max で割ると 2D の Uz や小さい成分が過大評価される [codex result-2 M2])・ρe/h0/k/ω/ρK/ρΩ/凝縮量 ≤ 1e-4、
  C = vis_turb ≤ 1e-2。診断量 (res_*, limiter_*, dt_local, wall_dist 等) は除外。run-to-run ノイズ床 (基準×基準 ≥3 本) は参考値として併記し、ノイズの 2 倍を超えた
  変数は項目別に原因を書く (2D node TP 300 step: ρ 2.3 倍・Ux 2.2 倍 = 絶対 1.2e-6 / 2.4e-6、float 熱力学の系統差)。
  ビット一致は「固定入力に対する単一カーネル出力」にだけ要求する (リミッタ template 化など)。
- **収束と準定常** (codex M5): 最終バイナリで run_0234 config を同 IC から 12000 step 走らせ、基準 run と**両方**が
  `check_convergence.py` で同じ verdict 区分 (全保存量の到達残差が同桁以上) であること、`check_quasisteady.py`
  (pmax/machmax) が両方 STEADY であること、壁 p/p0 の差 ≤ 0.1 % (壁圧の時系列も末尾 2000 step で頭打ち) を必須にする。
  未収束なら比較を確定しない。
- **12000 step 本 run の判定 (2026-09-12 実測, 表現は codex result M4 採用で改訂)**: run_0410 (thermoFloat, nStepInner 5) / run_0411 (3) /
  run_0413 (4) / run_0414 (4, cfl 8) は run_0234 と同 IC・同 config。**いずれも (基準 run_0234 も) `check_convergence` は NOT CONVERGED
  (plateau: roUy/roK 頭打ち・他 falling) で、「収束解」の比較ではない**。比較は「同一プラトー状態での有限時間継続の差」として報告する:
  到達残差は同桁 (rms_ro 1.2e-11 vs 6.4e-12, roe 6.0e-6 vs 5.0e-6, roOmega 1.8e-2 同値)、`check_quasisteady` pmax/machmax ALL STEADY、
  **壁 p/p0 の時系列** (3000 step 間隔スナップショット, x 10–95 mm): 末尾 2 枚の差 3.7e-4 (run_0413) / 2.2e-4 (run_0234) で両者とも同じ
  ゆっくりした漂いが残り、バイナリ間の差 (最終スナップショット比 3.4e-6〜8.1e-5) はその漂いより 1 桁小さい。
- **速度** (codex M4): 合否は `FORGE_PROFILE=0` (セクション同期なし)、warm-up 後、専用 run ディレクトリで基準/変更版を
  交互に 2 回以上回して ms/step (時間ループ内の壁時計 `Time = ...`; 初期化・I/O は含まない) の中央値で判定する。
  `FORGE_PROFILE=1` / nsys / ncu は原因分析用。FP64 稼働率 = 削減可能時間ではない。
- **CFL/sweep 比較**: 同じ到達残差 (全保存量) までの壁時計で評価する。

## 5. 実装ステップ

1. 計測基盤: `solver_density_cuda/tools/bench_steps.sh` (済)、ncu は `sudo` で `--log-file` に CSV、`cuobjdump -sass` で
   DFMA/DMUL 数を数える手順を methods/architecture/performance.md に記録。
2. §4.2-5 リミッタ template 化 (`limiter_d.cu`) → A/B (ビット同一)。
3. §4.2-1 リテラル昇格除去: `viscousFlux_d.cu`, `setDT_d.cu`, `convectiveFlux_slau_d.inc.cuh`, `convectiveFlux_common_d.cuh`,
   `convectiveFlux_boundary_d.inc.cuh`, `scalarTransport_d.cu`, `calcGradient_d.cu`, `ransSource_d.cu`, `turbulent_viscosity_d.cu`,
   `limiter_d.cu`, `boundaryCond_d.cu` (順に、各ファイルごとに SASS の DFMA 数と ms/step を確認)。
4. §4.2-2 化学種拡散: per-cell 前計算カーネル (`species_cell_props_d`: h_s(T_c), D_s^{lam}(T_c,P_c,X_c)) + 面カーネル float 化。
5. §4.2-2 SLAU TP 面エンタルピー: `thermo_h_mix_f` (float) と `thermo_cph_*_f`。
6. §4.2-3 `dependentVariables_d` float Newton (config キー付き)。
7. §4.2-4 block-DPLUR 対角キャッシュ / 占有率 / LU 保存。
8. §4.2-6 小物。
9. §4.2-7 nStepInner / CFL の壁時計比較 (run_0234 config の 12000 step を新バイナリで再実行)。
10. codex result レビュー → accepted へ。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 0 | 基準の固定 (codex M4/M7) | 基準バイナリ `~/forge-bin/forge_base`・入力ハッシュ・生ログ・ノイズ床 (`run_0401_perf_verify/cmp_base2.txt`) を保存 (済)。run_0234 の VERDICT/メッシュ品質/旧キー監査/IC 確認は README 表を参照。`bench_steps.sh` は専用 run・`FORGE_PROFILE=0` 既定・ワイルドカード削除撤去・失敗即終了に修正 |
| 1 | リミッタ template 化 (§4.2-5) | 済 (d331a29b)。81.9 ms/step (−1 ms)。関数ポインタは主因でなく、venkata の `2.0*` 昇格 (double 除算) が主因 → #2 で解消 |
| 2 | リテラル昇格除去 (§4.2-1) | 済: batch1 (limiter/viscous/setDT, e2fe1ee) 70.55 → batch2 (SLAU/common/boundary/scalar/gradient/ransSource/turb_visc/gasProperties/depVar CPG 経路, 5c1d7455) **67.4 ms/step** (交互 2 回 67.42/67.35 vs base 82.37/82.36)、ノイズ床内 (`cmp_lit2.txt`) |
| 3 | 化学種拡散の float 化 (§4.2-2, 面状態評価のまま) | 済 (fad80e54): `SpeciesThermoF` ミラー + `thermo_*_f`、species_diffusion 13.6→1.64 ms |
| 4 | SLAU TP 面エンタルピー float (§4.2-2) | 済 (fad80e54): SLAU 16.2→2.92 ms。#3+#4 で **44.0 ms/step** (43.97/44.02 vs base 82.37/82.38)、場はノイズ床内 (`cmp_thermof.txt`: vis_turb 2.17e-3 = 床の 1.17 倍, 他 ≤ 床) |
| 5 | block-DPLUR 対角キャッシュ・占有率 (§4.2-4) | **占有率実験も却下**: `BLOCK_DPLUR_MINBLOCKS` 4 (regs≤128) は不変 44.4/43.9、6 (regs≤85, spill) は 62.9 ms/step。|
| 5' | (旧 #5 記録) | 対角キャッシュは実装したが **不採用** (0fb8ee73, `blockDPLURDiagCache` 既定 0): 44.0→46.5 ms/step と遅化 (25 floats/cell の保存+4 読込 ≈1.2 GB/step > 省ける gather)。sweep はレイテンシ律速 (占有率 28 %) → `BLOCK_DPLUR_MINBLOCKS` で占有率実験中。次候補: dq の AoS 化 (5 配列→stride 5), nStepInner 5→3 の壁時計比較 (#8) |
| 6 | 面ループ融合・小物 (§4.2-6) | **済** (d2663103/f8de4ae0/600e03c3): k/ω 移流+拡散と化学種移流を多スカラー面カーネルに融合、`speciesFaceReconstruction=0` で読者の無い ∇Y を省略 (−1.2 ms)、DPLUR/LSQ の読み取り専用引数を `const __restrict__`。A/B `fuse3`: **35.5/36.0 ms/step** (base 82.34/82.35), 場はノイズ床内 (`cmp_fuse3.txt`)。小物 (d5637b31: DPLUR の未読 rhs 書込撤去, `fill_limiter5_d`, 定数転送の変化時のみ化) → A/B `small4`: **34.7 ms/step** (34.71/34.69 vs base 82.38/82.37), ノイズ床内 (`cmp_small4.txt`) |
| 7 | dependentVariables ハイブリッド反転 (§4.2-3) | 済 (045f4b5e): 単体検証 PASS (≤1e-8·T)。3D A/B (`thermoFloat: 1`): **38.4 ms/step** (38.44/38.41 vs base 82.41/82.37), 場はノイズ床内 (`cmp_thermofl.txt` vis_turb 1.90e-3, 他 ≤ 床)。**ユーザ決定 (2026-09-12): 既定 1** (datum 無しは自動 0 + 警告、明示 1 で datum 無しはエラー) |
| 8 | 陰解法 sweep 数の壁時計比較 (§4.2-7) | **済**: run_0410/0411/0412 (nStepInner 5/3/2, run_0234 と同 IC・config + thermoFloat, 12000 step)。**3 は 5 と全残差列が一致** (末尾平均 rms_ro 1.22e-11/1.23e-11, roe 6.11e-6 同値, 5 の到達残差に着く step も同じ ±60) で 38.9→34.6 ms/step (−11 %)。**2 は発散** (rms_roe 8.7e-6→4.5e-2 上昇, roOmega 79 で停滞)。→ **ユーザ決定 (2026-09-12): 推奨レシピは `nStepInner: 4`**。確認 run_0413 (sweep 4, `thermoFloat` 既定 1): 末尾残差が 5 と同値 (rms_ro 1.23e-11, roe 6.11e-6, roOmega 1.90e-2)、pmax/machmax STEADY、run_0234 比の壁 p/p0 差 ≤3.4e-6。(壁時計 68.5 ms/step は別セッションの chem run と GPU 共有のため無効)。**cfl_pseudo 6→8** (run_0414, sweep 4): 33.0 ms/step (クリーン)、cfl 6 の到達残差に着く step 11751→11222 (−4.5 %)、末尾残差は同等〜わずかに低い (roe 5.96e-6, roOmega 1.68e-2)、ALL STEADY、run_0234 比の壁 p/p0 差 8.1e-5 (cfl 6 の 3.4e-6 より大きいが ≪0.1 %)。効果が小さいので**推奨レシピは cfl 6 のまま** (8 は任意、+5 %)。|
| 9 | AoS gather (**不採用で確定**, opt-in 残置) | (a) block-DPLUR の近傍 dq を stride-8 AoS (`blockDPLURDqPack`, 874e54d1) + loop 0 の gather 省略 (dq_old≡0)、(b) 原始量 (ro,Ux,Uy,Uz,P,T) の AoS パック (`mesh.primPack`) を applyBconds 後に組み LSQ 勾配とリミッタが 1 セクタで gather。いずれも同じ値を別レイアウトで読むだけ (ビット同一)。**不採用 (既定 0, opt-in 記録)**: ローカル RTX 3060 の 3D 257k 節点 (`run_0452_perf_bench3d_coarse`, 200 step ×3) で dq パック 22.84/23.15 vs off 22.66/22.72、原始量パック on/off 22.9–23.1 vs 22.8–22.9 ms/step と差なし〜微増。gather は L2 ヒットで、セクタ数削減より追加書込が勝つ。loop 0 の gather 省略のみ残す |
| 10 | 節点の RCM 再番号付け (**A10G −3.3 %, opt-in で確定**) | `convertGmshToForge` に `mesh.renumber: rcm` (makeMesh 前に nodes と要素 iNodes を並べ替え、`/MESH/RENUMBER_PERM` に new→old を保存; `tools/permute_res_h5.py` で旧番号の res を移植)。ローカル coarse 3D 257k: 帯域 256481→3024、同一変換器の非 RCM と場の差はノイズ級 (k 5.5e-4, vis_turb 8e-4, 他 ≤1.4e-4)、速度 −1〜−4 % (ローカル GPU は別セッションと共有で不安定)。**注意**: 旧 h5 (2026-09-08 以前の変換) とは wall_dist 定義が 6.75 % 違うので比較は同一変換器で。2.37M ローカル RTX 3060 (`run_0455_perf_big_{norcm,rcm}`, run_0234/res_12000 を perm 移植, 100 step ×2): base 177/173 → head 63.9/66.0、RCM 65.0/62.0、RCM+dq 66.1/68.2、RCM+prim 65.6/64.1 ms/step = **RCM・パックとも差なし** (帯域 2369972→12000 でも gather は L2 で吸収済、sweep はレイテンシ律速)。**A10G (`run_0415_perf_rcm_baseline`, 100 step ×2, 同値再現)**: norcm 33.87 / **rcm 32.74 (−3.3 %)** / rcm+dq 35.34 / rcm+prim 33.41 / dq 36.58 ms/step → RCM は小さいが再現する利得、パック 2 種は逆効果 (既定 0 確定)。`mesh.renumber: rcm` は opt-in (既定 none): 節点順が変わるので旧 res の restart は `tools/permute_res_h5.py` (`/MESH/RENUMBER_PERM`) 経由 |
| 11 | 最終バイナリでの長時間比較と定常性の VERDICT (codex result-2 M1/M3) | **済**: `forge_final` (ソース 9333df5c = ソルバは ac9262e8、sha256 49f5db32…, 入力 sha256 は `run_0416_final_sweep4/INPUT_SHA256.txt`) で run_0416 (12000 step, nStepInner 4, 既定 thermoFloat): **32.62 ms/step** (A10G, クリーン)、check_convergence NOT CONVERGED plateau (基準 run_0234 と同区分: roUy/roK 頭打ち)、check_quasisteady pmax/machmax ALL STEADY、壁 p/p0 時系列 (2000 step 間隔 6 枚) の隣接差 2.2e-3→1.0e-3→5.4e-4→3.3e-4→2.2e-4 (単調減衰、末尾 <1e-3/2000 step = STEADY 判定, `WALL_PP0_SERIES_VERDICT.txt`)、run_0234 最終比 3.1e-6。証拠は `case/16.nozzle_wys/_aws_perf_evidence/run_0416_final_sweep4/`。「収束解不変」は主張しない (§4.3) |
| 12 | line-implicit 経路の回帰 (codex result-2 M4) | **済**: case/39 `run_0100_perf_regress_node_lineimplicit_dualtime` (node 周期 1.58 M 節点・dual-time 20 サブ反復・lineImplicit 1・lineKFreeze 1・SST-DDES, run_diag_lineimp2_frz2/res_100 から 20 step): 基準 ×2 vs 最終 **PASS 18/18** (ノイズ比判定, `perf_regress.py cmp --noise`)、1533→1365 ms/step |

## 6. 検証

- **単体 / ビルド**: AWS `~/forge-perf/solver_density_cuda/build` (Release, sm_86, nvcc 13.2) とローカル native (`.build-native`)。
- **速度**: `tools/bench_steps.sh run_0400_perf_baseline 100 <label>` の ms/step (起動込み) を各項目の前後で記録 (§9)。
- **同一性**: 各項目後に `run_04xx` で 100 step 継続 → `h5diff` (ビット同一項目) / 最大相対差 (float 化項目)。
- **収束**: 最終バイナリで run_0234 config を 12000 step 再実行 (`run_0410_perf_final_12k`)、`check_convergence.py` VERDICT と
  `wall_pp0.csv` を run_0234 と比較。
- **標準検証ケース**: case/20.naca_ml node/cell、case/08.bump verify (`verify/run_verification.sh`)、case/16 2D SST 継続。
- **ローカル回帰 (2026-09-12 実施, RTX 3060, 基準 = 0512823d の native ビルド, 各 300 step 継続, base×2 でノイズ床)**:
  case/20 `run_perf_regress_cell_slau_explicit` (cell, CPG 層流, 陽解法 RK3, IC run_slau/res_4000) / `run_perf_regress_cell_tpair_explicit`
  (cell, TP 空気, 陽解法, IC run_tp_air/res_2000) / case/16 `run_0450_perf_regress_node2d_tp` (node 2D TP SST 陰解法, IC run_0213/res_24000) /
  `run_0451_perf_regress_cell2d_tp` (cell 2D TP SST 陰解法, IC run_0196/res_24000)。new と new+`thermoFloat:1` の全場差はいずれも
  base×base ノイズと同程度 (P/T/ρ/U/k/ω/ρY ≤ 1.5 倍)。注記: 2D 平面 (naca) で面外速度 Uz が基準 3.7e-10 → 新 2.8e-6 m/s
  (主流 ~300 m/s の 1e-8, 旧 double 中間演算で打ち消していた z 成分の float 丸め)。物理的影響なし。
  `run_cmp_node_sst`/`run_cmp_cell_sst` は旧キー `LESorRANS` で現行 config 非互換のため使わず。
- **最終版 (ac9262e8, WALE double 復帰・研磨反復化後) の追加回帰 (2026-09-12, codex result M5 採用, `tools/perf_regress.py`)**:
  case/09 `run_0052_perf_regress_node_periodic_dualtime` (node 周期・dual-time・SST・E_t): P/T/k/ω/μt ≤4.3e-7 = base×base 同等
  (速度は base×base 自体 0.32 の減衰乱流ばらつきで比較不能); case/23 `run_0100_perf_regress_cell_axisym_tp_implicit` (cell 軸対称 TP
  陰解法 WALE): 全場 ≤3.4e-6 = base×base 同等; case/16 `run_0456_perf_regress_node2d_cond` (node 2D SST 凝縮): 全場 ≤1.6e-5、
  g_0/rog_0 3〜4e-6 = base×base 同等; `run_0450` は基準 5 本で再計測 (ρ 2.3 倍・Ux 2.2 倍 = 絶対 1.2e-6/2.4e-6, §4.3 絶対基準内)。
  line-implicit は §5.1 #12 (PASS)。判定ツールは `tools/perf_regress.py cmp --noise <基準ラベル>` (絶対基準 or 2×ノイズ床で場ごとに合否、勾配・診断量は除外):
  2026-09-12 再集計 — node 2D TP (基準 5 本) / cell 2D TP / node 2D 凝縮 / naca cell CPG 陽 / naca cell TP 陽 (基準 4 本; ρ がノイズの 2.0 倍 = 2.2e-5) /
  cell 軸対称 TP WALE / node 周期 dual-time (速度は base×base 0.28 の減衰乱流ばらつきでノイズ比判定) / node line-implicit: **全て PASS**。
  温度反転の単体試験 (`tools/test_thermo_float.cpp`): 誤差 ≤1.0e-8·T、float 格納 roe の 10 往復ドリフトは従来 double 反転と同値 (float 格納固有, 7e-8〜4.4e-7·T)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result | `2026-09-12` | [`notes/reviews/2026-09-12-performance-3d-node-sst-speedup-result.md`](../../notes/reviews/2026-09-12-performance-3d-node-sst-speedup-result.md) | **NO-GO**, C0/M5/m2 | M1 採用 (WALE/SIGMA の float pow を double に戻す, ac9262e8) / M2 採用 (研磨を収束まで最大 3 段, 単体試験を本番 12 反復・冷間開始・float 組成・再格納ドリフト込みに: errHyb ≤3.4e-10·T PASS) / M3 採用 (基準 5 本でノイズ再計測、§4.3 を絶対値基準に改訂、「全てノイズ床内」撤回) / M4 採用 (「収束解不変」撤回、同一プラトーの有限時間差と壁圧時系列で記述) / M5 採用 (AWS 証拠を `case/16.nozzle_wys/_aws_perf_evidence/` に回収: bench ログ・cmp・VERDICT・壁圧・残差 CSV・バイナリ/メッシュ sha256; 周期/軸対称/凝縮の回帰は §6 に追記) / m6 採用 (bench の上書き拒否、手順書を実装に合わせ更新) / m7 採用 (thermophysics.md / performance.md / §2 / §5.1 / §8 / README 索引を整理)。修正後に result 段を再実行する |
| result (2 回目) | `2026-09-12` | [`notes/reviews/2026-09-12-performance-3d-node-sst-speedup-result-2.md`](../../notes/reviews/2026-09-12-performance-3d-node-sst-speedup-result-2.md) | **NO-GO**, C0/M4/m2 | M1 採用 (recommended-settings の「収束解不変」撤回, 未収束区分の明記) / M2 採用 (速度の正規化を \|U\| 尺度に、絶対基準 or 2×ノイズ床の判定に確定、全回帰を再集計し EXCEED 0) / M3 採用 (最終バイナリ run_0416 + sha256 + 壁圧時系列 VERDICT) / M4 採用 (line-implicit 回帰 PASS) / m5 採用 (float 格納往復ドリフト試験, double 反転と同値) / m6 採用 (performance.md 研磨段数・判定基準, plans/README, §5.1 見出し整理)。3 回目を実行 |
| plan | `2026-09-12` | [`notes/reviews/2026-09-12-performance-3d-node-sst-speedup-plan.md`](../../notes/reviews/2026-09-12-performance-3d-node-sst-speedup-plan.md) | GO-with-changes, C0/M7/m2 | M1 採用 (セル前計算→面状態 float 評価に変更, §4.2-2) / M2 採用 (Newton float は後段+単体検証, §5.1 #7) / M3 採用 (ノイズ床基準, §4.3) / M4 採用 (PROFILE=0・交互実行, §4.3, bench_steps.sh) / M5 採用 (両 run PASS+STEADY+壁圧, §4.3) / M6 採用 (float point 経路限定+回帰, §4.2-4) / M7 採用 (専用 run・削除撤去, §5.1 #0) / m8 採用 (M6 関係を §4.1 に明記) / m9 採用 (1/64, launch_bounds 記述訂正) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/`: `limiter_d.cu`, `viscousFlux_d.cu`, `setDT_d.cu`, `convection/*.inc.cuh`, `convectiveFlux_common_d.cuh`,
  `scalarTransport_d.cu`, `speciesTransport_d.cu`, `thermo_d.cuh`, `dependentVariables_d.cu`, `timeIntegration_d.cu`, `gasProperties_d.cu`,
  `calcGradient_d.cu`, `ransSource_d.cu`, `turbulent_viscosity_d.cu`, `boundaryCond_d.cu`
- `solver_density_cuda/main.cpp` (blockDPLURSolve), `input/solverConfig.*` (`thermoFloat`)
- `solver_density_cuda/tools/bench_steps.sh` (新規)
- docs: `methods/architecture/performance.md` (新規), `methods/index.md`, `procedures/development-environment.md` (ncu sudo 注記)

## 8. 完了条件

- [ ] `methods/architecture/performance.md` を作成し index に登録
- [ ] §5.1 #1–#10 の実装と A/B (§4.3 の判定) 完了 (#5/#9 は不採用で確定、#10 RCM は opt-in)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status: done`、§9 変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-12` — 起票。ベースライン計測 (run_0400_perf_baseline @A10G 82.85 ms/step)、nsys/ncu で FP64 律速を同定 (§4.1)。
- `2026-09-12` — codex plan レビュー (GO-with-changes, M7/m2) を全件採用し §4.2/§4.3/§5.1 を改訂。リミッタ template 化 (−1 ms) と batch1 リテラル修正 (82.85→70.55 ms/step, ノイズ床内) を実測。
- `2026-09-12` — DPLUR 対角キャッシュ (§4.2-4) は 46.5 ms/step と遅化したため既定 0 (opt-in 記録)。sweep は gather 数でなくレイテンシ律速。`__launch_bounds__` minBlocks 4/6 も不変/悪化。
- `2026-09-12` — ユーザ決定: `thermoFloat` 既定 1、推奨 `nStepInner` 4 (recommended-settings.md 反映)。run_0413 (sweep 4) / run_0414 (cfl 8) で確認 (§5.1 #8)。RCM 再番号付けは A10G −3.3 % (opt-in)、AoS パック 2 種は逆効果で既定 0 (§5.1 #9, #10)。
- `2026-09-12` — codex result レビュー NO-GO (M5/m2) を全件採用: WALE double 復帰・研磨反復化・判定基準の絶対値化・「収束解不変」の撤回・証拠回収 (§6.1)。最終バイナリ (ac9262e8) は A10G 33.9 ms/step、場は基準内。
- `2026-09-12` — codex result レビュー 2 回目 NO-GO (M4/m2) を全件採用: 最終バイナリで run_0416 (32.6 ms/step, 壁圧時系列 STEADY, 証拠+sha256)、判定を「絶対基準 or 2×ノイズ床」に確定し全回帰を再集計 (全 PASS)、line-implicit 回帰追加 (PASS)、単体試験に float 格納ドリフト (double 反転と同値)、recommended-settings/methods/plans README の撤回済み記述を修正。
- `2026-09-12` — 小物 3 件で **34.7 ms/step** (small4)。累積 82.4→34.7 (2.37 倍)、nStepInner 3 併用で ≈30.5 ms/step (2.7 倍)。
- `2026-09-12` — 12000 step 本 run (run_0410/0411) は run_0234 と壁 p/p0 差 ≤1e-5・STEADY・同 verdict 区分 (§4.3)。面ループ融合+∇Y 省略+restrict で **35.7 ms/step** (fuse3)。sweep 数比較: nStepInner 3 は 5 と残差経路一致で 34.6 ms/step、2 は発散 (§5.1 #8)。累積: 82.4 → 35.7 (2.3 倍)、nStepInner 3 併用で ≈31.5 ms/step (2.6 倍)。
- `2026-09-12` — ハイブリッド温度反転 (`thermoFloat: 1`) で **38.4 ms/step** (累積 82.4→38.4, 2.15 倍)。ローカル回帰 4 ケース (node/cell × CPG/TP × 陽/陰) は全てノイズ床内 (§6)。k/ω・化学種の面ループ融合、未使用 ∇Y の省略、DPLUR/LSQ の `__restrict__` を実装 (A/B 待ち)。sweep 数比較 run_0410–0412 投入。
- `2026-09-12` — batch2 リテラル修正 67.4 ms/step、thermo float ミラー (SLAU h_mix + 化学種拡散) で **44.0 ms/step** (−47 %)。再プロファイル: block-DPLUR 5 sweep 14.9 (34 %) / dependentVariables 6.6 / SLAU 2.9 / lsqPreGrad 2.8 / k-ω 輸送 3.3 / species_diffusion 1.6 / viscous 1.5 / limiter 1.25。次は DPLUR 対角キャッシュ・占有率、dependentVariables float Newton (単体検証付き)、k/ω 面ループ融合。
