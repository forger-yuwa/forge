# 性能プロファイルと数値精度の方針

## 1. 目的

`solver_density_cuda` の 1 ステップがどこで時間を使っているか (現状のプロファイル)、計測の手順、
そしてカーネル内の数値精度 (float32 / double) をどう使い分けるかの現在仕様をまとめる。
設計判断の経緯は [`plans/active/performance-3d-node-sst-speedup.md`](../../plans/active/performance-3d-node-sst-speedup.md)。

## 2. 計測手順 (native)

速度評価は native で行う ([`procedures/development-environment.md`](../../procedures/development-environment.md))。

1. **ms/step**: `solver_density_cuda/tools/bench_steps.sh <run_dir> <nsteps> <label>`。run の `solverConfig.yaml` を
   `nStepOuter=<nsteps>`・出力なしに書き換え (`solverConfig.yaml.orig` を保存)、`FORGE_PROFILE=1` で回して
   `Time = ... (wall, N steps, X ms/step)` と "Runtime Profile Summary" (セクション別) を `bench_<label>_n<N>.log` に残す。
   起動 (メッシュ読込・初期化) を相殺するには 2 種類の step 数で回して差分を取る。
2. **カーネル別時間**: `nsys profile --trace=cuda --sample=none --cpuctxsw=none -o <out> forge` →
   `nsys stats --report cuda_gpu_kern_sum --format csv <out>.nsys-rep`。
3. **律速の種別** (演算/メモリ/FP64): `ncu` は GPU 性能カウンタの権限が要るので **`sudo -E` で起動**し、
   `--csv --log-file <csv>` に書く (stdout に混ぜると forge の出力と混ざる)。見る指標:
   `sm__throughput` / `gpu__compute_memory_throughput` (どちらが天井か)、`sm__pipe_fp64_cycles_active` (FP64 パイプ稼働率;
   CC 8.6 では FP64 は FP32 の 1/32 スループットなので数 % でも支配的になる)、`sm__warps_active` (占有率)、
   `launch__registers_per_thread`、`lts__t_sector_hit_rate` (L2 ヒット率)、`smsp__warp_issue_stalled_long_scoreboard_*` (メモリ待ち)。
4. **FP64 命令の混入確認**: `cuobjdump -sass forge` でカーネルごとに `DFMA/DMUL/DADD/F2F` を数える。
   `float` 式に接尾辞なしリテラル (`0.5*x`) があると double に昇格し、これらが現れる。

## 3. 現状プロファイル (2026-09-12, A10G, case/16 3D node SST TP 2.37 M 節点)

ベースライン 82.85 ms/step (`run_0400_perf_baseline`, 収束場から 100 step 継続)。

| カーネル | ms/step | 律速 |
| --- | --- | --- |
| `SLAU_d` (対流流束, 7.04 M 面) | 16.2 | FP64 (面ごとの `thermo_h_mix` double + リテラル昇格) |
| `implicit_defect_correction_block_d` ×5 sweep | 14.9 | メモリ (gather, 占有率 28 %) |
| `species_diffusion_d` | 13.6 | FP64 (全演算 double) |
| `limiter_r1_fused5_d` | 9.5 | 関数ポインタ経由の呼び出し |
| `dependentVariables_d` | 6.1 | FP64 (TP Newton) |
| `viscousFlux_d` | 4.7 | FP64 (リテラル昇格のみ) |
| `lsqPreGrad_internal_d` | 2.7 | gather |
| その他 (k/ω 輸送・ソース、setDT、gasProperties、境界、残差 reduction) | ~15 | — |

host 側 (残差ログ・モニタ) のオーバーヘッドは無視できる。

## 4. 数値精度の方針

- **状態・残差・勾配は float32** (`flow_float`)。GPU メモリ帯域と CC 8.6 の FP64 スループット (1/32) のため、
  カーネル内の一時演算も **float32 を既定**とする。
- **リテラルは必ず float にする**: `0.5f`、`static_cast<flow_float>(2.0/3.0)` など。`flow_float x = 0.5*y` は
  double 演算 + 2 回の変換になり、面ループでは 1 命令が FP32 の 64 倍のコストになる。
  `block-DPLUR` カーネル (`static_cast<ST>` で徹底) が手本。
- **double を使う箇所は明示し、根拠を残す**: 幾何前処理 (双対体積・重心・閉性 Σr_f S_f の桁落ち対策)、
  周期・軸対称の閉性、凝縮 EOS の (T,g) 同時反転など、桁落ちが実測で問題になった箇所に限る。
  熱力学 (NASA-9) の**面ごと**評価は float 版を使い、セルごとに前計算できる量 (h_s(T_c), D_s) はセル配列に置く。
- 精度変更は「ビット同一」か「相対差 ≤ 1e-4 で収束解不変」のどちらかを A/B で示してから採用する
  (判定基準は plan §4.3)。
