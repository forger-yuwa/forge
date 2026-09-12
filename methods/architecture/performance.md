# 性能プロファイルと数値精度の方針

## 1. 目的

`solver_density_cuda` の 1 ステップがどこで時間を使っているか (現状のプロファイル)、計測の手順、
そしてカーネル内の数値精度 (float32 / double) をどう使い分けるかの現在仕様をまとめる。
設計判断の経緯は [`plans/active/performance-3d-node-sst-speedup.md`](../../plans/active/performance-3d-node-sst-speedup.md)。

## 2. 計測手順 (native)

速度評価は native で行う ([`procedures/development-environment.md`](../../procedures/development-environment.md))。

1. **ms/step**: `solver_density_cuda/tools/bench_steps.sh <template_run> <nsteps> <label>`。template の入力 (config・bcond・species・
   probe・メッシュ h5 のハードリンク) から専用ディレクトリ `<template_run>_bench/<label>_n<N>/` を作り、`nStepOuter=<N>`・出力なし
   (`outStepInterval` 巨大) で回して `Time = ... (wall, N steps, X ms/step)` を `bench_<日時>.log` に残す (バイナリ・入力の sha256 も記録)。
   既存 label には上書きしない (別 label か `FORGE_BENCH_OVERWRITE=1`)。既定 `FORGE_PROFILE=0`; `FORGE_PROFILE=1` はセクション別計時用で
   同期が入るので合否判定には使わない。`FORGE_CFG_SUB="old=new"` で config を差し替えて A/B できる。計時は時間ループ全体 (`main.cpp` の
   壁時計; 初期化・I/O は含まない) で warm-up 区間の除外はしていないので、A/B は基準/変更版を**交互に 2 回以上**回し中央値で比べる
   (100 step で 1 step 目の影響は 1 % 未満)。
2. **カーネル別時間**: `nsys profile --trace=cuda --sample=none --cpuctxsw=none -o <out> forge` →
   `nsys stats --report cuda_gpu_kern_sum --format csv <out>.nsys-rep`。
3. **律速の種別** (演算/メモリ/FP64): `ncu` は GPU 性能カウンタの権限が要るので **`sudo -E` で起動**し、
   `--csv --log-file <csv>` に書く (stdout に混ぜると forge の出力と混ざる)。見る指標:
   `sm__throughput` / `gpu__compute_memory_throughput` (どちらが天井か)、`sm__pipe_fp64_cycles_active` (FP64 パイプ稼働率;
   CC 8.6 では FP64 の加算/乗算/FMA は FP32 の 2/128 = 1/64 スループットなので数 % でも支配的になる)、`sm__warps_active` (占有率)、
   `launch__registers_per_thread`、`lts__t_sector_hit_rate` (L2 ヒット率)、`smsp__warp_issue_stalled_long_scoreboard_*` (メモリ待ち)。
4. **FP64 命令の混入確認**: `cuobjdump -sass forge` でカーネルごとに `DFMA/DMUL/DADD/F2F` を数える。
   `float` 式に接尾辞なしリテラル (`0.5*x`) があると double に昇格し、これらが現れる。

## 3. 現状プロファイル (2026-09-12, A10G, case/16 3D node SST TP 2.37 M 節点)

ベースライン 82.85 ms/step → 高速化後 **34.7 ms/step** (`thermoFloat: 1`, nStepInner 5)。`run_0400_perf_baseline_bench/prof_*`。

| カーネル | 基準 ms/step | 最終 ms/step | 処置 / 律速 |
| --- | --- | --- | --- |
| `SLAU_d` (対流流束, 7.04 M 面) | 16.2 | 2.9 | リテラル f、面エンタルピー h_mix(Y_f,T_f) を `SpeciesThermoF` (float 係数) で評価 |
| `implicit_defect_correction_block_d` ×5 sweep | 14.9 | 10.9 | `__restrict__`・未読 rhs 撤去。メモリ (レイテンシ) 律速、占有率 28 %。対角キャッシュ / launch_bounds は逆効果 |
| `species_diffusion_d` | 13.6 | 1.6 | 面状態組立と h_s(T_f) を float (離散式不変) |
| `limiter_r1_fused5_d` | 9.5 | 1.2 | venkata の `2.0*` (double 除算) → f、関数ポインタ除去 |
| `dependentVariables_d` | 6.4 | 3.9 | ハイブリッド温度反転 (float Newton + double 1 段研磨) |
| `viscousFlux_d` | 4.7 | 1.5 | リテラル f のみ |
| `lsqPreGrad_internal_d` | 2.7 | 2.8 | gather 律速 (原始量 AoS パックで改善見込み) |
| k/ω 移流+拡散・化学種移流 | 3.3 | 2.3 | 多スカラー面カーネルに融合 (6 起動→3) |
| `species_gradient_d` (+normalize) | 1.2 | 0 | 読者が無い (speciesFaceReconstruction=0) ので省略 |
| その他 (setDT・ソース・境界・reduction・更新) | ~10 | ~8 | — |

**陰解法 sweep 数**: 同 IC・同 config の 12000 step で `nStepInner` 3 と 5 は全残差列が一致、2 は発散 (run_0410–0412)。推奨レシピは 4
(procedures/recommended-settings.md)。

## 4. 数値精度の方針

- **状態・残差・勾配は float32** (`flow_float`)。GPU メモリ帯域と CC 8.6 の FP64 スループット (1/64) のため、
  カーネル内の一時演算も **float32 を既定**とする。
- **リテラルは必ず float にする**: `0.5f`、`static_cast<flow_float>(2.0/3.0)` など。`flow_float x = 0.5*y` は
  double 演算 + 2 回の変換になり、面ループでは 1 命令が FP32 の 64 倍のコストになる (実測: limiter/viscousFlux/setDT の接尾辞付与だけで 82.85→70.55 ms/step)。
  `block-DPLUR` カーネル (`static_cast<ST>` で徹底) が手本。
- **double を使う箇所は明示し、根拠を残す**: 幾何前処理 (双対体積・重心・閉性 Σr_f S_f の桁落ち対策)、
  周期・軸対称の閉性、凝縮 EOS の (T,g) 同時反転など、桁落ちが実測で問題になった箇所に限る。
  熱力学 (NASA-9) の**面ごと**評価は float 版 (`SpeciesThermoF`) を使う。評価点は従来どおり面状態 (T_f, P_f, Y_f) で、セル値の補間に
  置き換えない (離散式が変わる: codex plan レビュー M1)。WALE/SIGMA の高次べき乗 (`pow(x, 5/2)` 等) は float では高勾配で overflow するため double のまま。
- **TP の温度反転** (`physProp.thermoFloat`, 既定 1): float Newton (warm start, 最大 12 反復) → double Newton 1 段研磨。double 評価は cph_mix 1 回。
  誤差 ≤1e-8·T (float の T 格納 ulp 6e-8 未満; `tools/test_thermo_float.cpp`)。`thermoHrefTemp>0` が前提 (絶対 datum の H2O は float 段が収束しない)。
- **gather 系の試み (結果)**: 近傍 dq / 原始量の stride-8 AoS パック (`blockDPLURDqPack`, `mesh.primPack`) は A10G で +0.7〜+2.7 ms/step の逆効果、
  節点 RCM 再番号付け (`mesh.renumber: rcm`) は −3.3 %。gather は既に L2 で吸収されておりレイテンシ律速 (占有率 28 %) が残る。パックは既定 0、RCM は opt-in。
- 面流束の `atomicAdd` 蓄積のため同一バイナリでも全場はビット一致しない。精度変更は、基準バイナリ同士の run-to-run
  ノイズ床 (場のスケールで正規化した最大差) を先に測り、その 2 倍以内に収まることを A/B で示してから採用する (判定基準は plan §4.3)。
