# 諮問ブリーフ: 状態 float・差を取る幾何の量を double の座標から作る設計 (2026-10-10)

- 依頼者: ユーザ「float化よろ」 (前回の諮問 `notes/reviews/2026-10-10-m64-and-float-diagnose.md` の「案 (b) を先に」を受けて)
- plan: `plans/active/architecture-float-state-double-geometry.md` (draft、§4・§6 が案)
- 調査: `notes/investigations/2026-10-10-geometry-precision-inventory.md` (座標の使われ方の全表)
- コード: `fd2d4c3c` 以降 (ブランチ `feature/nozzle-wall-fit-and-pipeline`)

## 0. 問い

1. **§4.2 の量の選び方は足りているか、過剰か。** 面ごとに e = cc1 − cc0、|e|、δ/|e| = |S|²/|e·S|、r0・r1 = pc − cc を double で作って float で渡す。読み込み時の量 (LSQ・closure・壁関数の代表距離・d1・delta_les) も double の座標から作る。
   - 抜けている使い道はないか。node 離散化の境界の半割面、周期の継ぎ目、軸の節点、3D の扱いなど。
   - 「辺中点の再構成は 0.5·e で足りる」という整理は正しいか。
2. **§4.3 の既定の扱い。** 案 A (切り替えなし、常に新しい経路) と案 B (opt-in) のどちらにすべきか。案 A なら、FP64 のビルドで「丸めの順序の範囲で同一」(§6 V2) をどう確かめれば十分か。
3. **§6 の確かめ方と基準を、事前登録できる形にしてほしい。** 特に次の 3 つ。
   - V3 (同じ状態の残差) の物差し
   - V4 (到達) の閾値 (R3 と同じ θ_r 0.05 %・Q_w 0.1 % でよいか)
   - V6 (標準ケースの回帰) の対象ケースと物差し
4. **状態の commit の精度 (`qAccumulatorFP64` の軸対称対応) を、この plan に含めるべきか。** 呼び出し側は「V4 で停滞が見えたら別に立てる」としている。
5. 実装の順序とリスク。カーネル 20 余りを一度に変えるか、段階に分けて各段で V2 を回すか。

## 1. 観測事実

- **速さ** (line-implicit-speed §6.21、`case/45.isobutane_m6_d155/run_0383_t5_*`): B0 の構成 (値 0 のライン・LAYOUT2) で、FP64 のビルド (flow・geom とも double) 32.0 ms/step、float のビルド 20.5 ms/step。
  - float は 1000 step 後の残差が FP64 の約 10 倍 (`rms_ro` 9.5e-5)。
- **冷却壁の第一層** (メモ [cooled-wall-mesh-precision]、2026-10-08):
  - float の座標では、第一層厚 / 半径 = f のとき第一層が約 f / 1.2e-7 ulp。f = 3.4e-7 で第一層厚が最大 28 % ずれた。
  - case/45 の第一層は 0.06〜12 µm、半径は 0.1〜0.8 m (半径方向の格子は比 1.07〜1.09 の等比)。
- **座標の使われ方** (調査の §0〜§5):
  - 幾何の配列は `flow_float` の表に入っていて、`geom_float ≠ flow_float` ではビルドが通らない。
  - ソルバは HDF5 の座標を `geom_float` に丸めて読む。
  - 毎 step、カーネルが絶対座標の差 (dcc・pc − cc・壁関数の y) を float で作る。
- **B0 の水準での 1 step の変化と float の刻み** (`run_0387_dqulp_B0`、FP64、B0 の最終状態から 10 step、9 → 10):
  - |Δq| < ½ ULP_f32(q) の節点の割合: ρ 17.7 %・ρu 8.9 %・ρv 3.5 %・ρE 11.9 %・ρk 3.4 %・ρω 5.9 %。
  - 比の中央値は 2.2・4.9・28・4.1・31・20。
- **`qAccumulatorFP64`**: 軸対称は `main.cpp:3371-3374` で拒否している。理由は `enforceAxisSymmetry` が commit の基準を射影するのに Qacc に効かないこと (幾何の精度とは関係ない)。

## 2. 期待値と出典

- FP64 のビルドでは、カーネルの差も double なので、新しい経路は今と丸めの順序の範囲で一致するはず (呼び出し側の推論。未確認)。
- 到達の閾値の候補は line-implicit-speed §6.18 の R3 (θ_r 0.05 %、Q_w 0.1 %)。
- 回帰の物差しは「基準コミットの double のビルドを真値として、新旧の float の距離を比べる」(メモ [float-regress-double-truth])。

## 3. 再現条件

- case/45 の B0: `cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext <run> --cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --cap 50`。FP64 のバイナリは `~/forge-linespeed-fp64` (lineM_fp64)、float は `~/forge-linespeed-f32`。

## 4. 実施済み

- 調査 (上記)。float と FP64 の 1 step の時間。|Δq|/ULP の測定。plan の draft。

## 5. 呼び出し側の仮説・案 (検証していない)

- **H1**: 精度を壊しているのは、毎 step の絶対座標の差 (特に壁節点と第一内部節点の辺の δ/dcc) である。読み込み時の量 (面積・体積) は、FP64 の変換器が double で作って float に丸めるだけなので、相対精度で足りる。
- **H2**: 状態 float の commit の丸めは、B0 の水準の揺れ (刻みの 2〜30 倍) に埋もれ、場全体の停滞は起きない。ただし θ_r のゆっくりした漂いが進むかは不明。
- **案**: 案 A (常に新しい経路)。実装は 2 段 (1: 読み込みと粘性・拡散・対角、2: 再構成・リミタ・読み込み時の量) に分け、各段で V2 を回す。

## 6. 読んでよいファイル

- `plans/active/architecture-float-state-double-geometry.md`
- `notes/investigations/2026-10-10-geometry-precision-inventory.md`
- `notes/reviews/2026-10-10-m64-and-float-diagnose.md`
- `methods/architecture/overview.md` (§6.2a)、`methods/time_integration/implementation.md` (commit の丸め)
- コード (grep で該当箇所だけ):
  - `solver_density_cuda/mesh/mesh.cpp`、`solver_density_cuda/variables.cpp`・`variables.hpp`
  - `solver_density_cuda/cuda_forge/viscousFlux_d.cu`・`scalarTransport_d.cu`・`calcGradient_d.cu`・`calcStructualVariables_d.cu`
  - `solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh`・`limiter_d.cu`・`ransWallFunction_d.cu`
  - `solver_density_cuda/cuda_forge/timeIntegration_d.cu` (8000 行超)
- 巨大なファイル (`*.h5`・ログ・`residual_history.csv`・`notes/reviews/*.log`) は読まないこと。
