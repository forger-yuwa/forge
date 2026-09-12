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
  小カーネル融合。各項目は A/B 計測 + 場の一致確認を伴う。
- **やらない**: MPI 多 GPU 化、メッシュ再番号付け (locality) は候補として記録するのみ (§5.1 末尾)、
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
3. **`dependentVariables_d` の Newton (後段・検証付き)**: float 版 `thermo_T_from_e_f` は **本 plan の後段**に置き、
   double 参照との単体検証 (使用 DB・組成端点・微量成分・50–6000 K・区間境界・外挿・datum 有無、温度誤差 + エネルギー残差 +
   反復上限 + 反転→再構成ドリフト) を通してから opt-in キー `physProp.thermoFloat` で入れる (codex M2)。未検証条件は double のまま。
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
- **場の一致**: 各項目後に同 IC から 100 step 継続し、`VALUE/*` 全量 (ro,P,T,Ux,Uy,Uz,roe,roY*,Y*,k,ω,vis_turb,h0 …) の
  正規化最大差がノイズ床の **2 倍以内**なら「収束解不変」の候補、超えたら原因を切り分ける (`cmp_h5.py`)。
  ビット一致は「固定入力に対する単一カーネル出力」にだけ要求する (リミッタ template 化など)。
- **収束と準定常** (codex M5): 最終バイナリで run_0234 config を同 IC から 12000 step 走らせ、基準 run と**両方**が
  `check_convergence.py` で同じ verdict 区分 (全保存量の到達残差が同桁以上) であること、`check_quasisteady.py`
  (pmax/machmax) が両方 STEADY であること、壁 p/p0 の差 ≤ 0.1 % (壁圧の時系列も末尾 2000 step で頭打ち) を必須にする。
  未収束なら比較を確定しない。
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
| 2 | リテラル昇格除去 (§4.2-1) | batch1 (limiter/viscous/setDT, e2fe1ee) 済 → 70.55 ms/step、ノイズ床内。batch2 (SLAU/common/boundary/scalar/gradient/ransSource/turb_visc/gasProperties/depVar CPG 経路) 適用中 |
| 3 | 化学種拡散の float 化 (§4.2-2, 面状態評価のまま) | `speciesTransport_d.cu`: 面状態 Y_f/T_f/P_f・J_s・補正・`thermo_h_mass_f(T_f)` を float に |
| 4 | SLAU TP 面エンタルピー float (§4.2-2) | `thermo_d.cuh` float 版 + `convectiveFlux_slau_d.inc.cuh` |
| 5 | block-DPLUR 対角キャッシュ・占有率 (§4.2-4, float point 経路限定) | `timeIntegration_d.cu` / `main.cpp blockDPLURSolve`; Taylor-Green (周期)・軸対称・等温壁の回帰を追加 |
| 6 | 小物 (§4.2-6) | `limiter_d.cu` fill 融合, `gasProperties_d.cu` powf, `convectiveFlux_d.cu` の毎ステップ `cudaMemcpyToSymbol` |
| 7 | dependentVariables float Newton (§4.2-3, 後段) | 単体検証ツール (double 参照) → opt-in `physProp.thermoFloat` |
| 8 | 陰解法 sweep/CFL の壁時計比較 (§4.2-7) | run_0234 config で 12000 step、同一到達残差までの壁時計 |
| 9 | 候補 (未着手): メッシュ再番号付け (RCM/Hilbert) で gather の L2 ヒット率改善 | L2 hit 62–75 % の改善余地。変換器側 |

## 6. 検証

- **単体 / ビルド**: AWS `~/forge-perf/solver_density_cuda/build` (Release, sm_86, nvcc 13.2) とローカル native (`.build-native`)。
- **速度**: `tools/bench_steps.sh run_0400_perf_baseline 100 <label>` の ms/step (起動込み) を各項目の前後で記録 (§9)。
- **同一性**: 各項目後に `run_04xx` で 100 step 継続 → `h5diff` (ビット同一項目) / 最大相対差 (float 化項目)。
- **収束**: 最終バイナリで run_0234 config を 12000 step 再実行 (`run_0410_perf_final_12k`)、`check_convergence.py` VERDICT と
  `wall_pp0.csv` を run_0234 と比較。
- **標準検証ケース**: case/20.naca_ml node/cell、case/08.bump verify (`verify/run_verification.sh`)、case/16 2D SST 継続。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
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
- [ ] §5.1 #1–#8 の実装と A/B (§4.3 の判定) 完了
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status: done`、§9 変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-12` — 起票。ベースライン計測 (run_0400_perf_baseline @A10G 82.85 ms/step)、nsys/ncu で FP64 律速を同定 (§4.1)。
- `2026-09-12` — codex plan レビュー (GO-with-changes, M7/m2) を全件採用し §4.2/§4.3/§5.1 を改訂。リミッタ template 化 (−1 ms) と batch1 リテラル修正 (82.85→70.55 ms/step, ノイズ床内) を実測。
