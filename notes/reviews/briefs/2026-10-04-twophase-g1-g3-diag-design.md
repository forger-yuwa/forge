# 諮問ブリーフ: 二相拡散 既定化 plan の G1 (0 step 作用素 A/B) と G3 (閉じた収支) の診断出力の設計

エスカレーション条件 6 (cuda_forge に診断を足す編集の前) と 1。plan: `plans/active/condensation-two-phase-default.md` §5.1 #4、§6 G1・G3 (事前登録済み)。
親 plan: `plans/active/condensation-two-phase-transport.md` §5.1 #4k。ユーザ了承済み (「ほい」2026-10-04)。

## 前提 (実装の現状)
- ON の面流束: `speciesTransport_d.cu:2045` `tp_build_face_in` → `twoPhaseDiffusion_d.cuh` `tp_face_flux<R>` (R = float 本番、double は収束受入の独立残差監査 #4f が格納値を double に上げて同じ式を評価: `speciesTransport_d.cu:2330` 付近)。出力 `TpFaceOutT`: J[k] (水は総水分)、Jv、Jl、JQ[3]、q (エネルギー)、diag、Sm、up0。
- OFF の面流束: `species_diffusion_d` (`speciesTransport_d.cu:224`) が面ごとに J_s = −ρ_f (D_s + D_t) ∇Y_s·S と Σh_s J_s を計算し、残差へ atomicAdd するだけ (面の値は外に出ない)。液・Q は拡散しない。
- 組立は前半 (射影・EOS・BC・勾配) / 後半 (リミッタ以降) に分割済み (`main.cpp` `assembleResidualPre/Post`、limiter plan の診断で導入)。残差の段ごとの取得の前例: `captureNodeIsothermalEnergyResidual` (main.cpp の組立途中で res を写す既定 no-op の診断)。
- 補正の記録: `[twophase-corr-gate]` (commit・floor・clamp・projection・removal を成分別に)、`[renorm-gate]`、`condCorrReasons` (理由別の補正量)。
- 共通入力: `case/16.nozzle_wys/run_0520_twophase_off_16k` res_16000 (OFF の最終場、AWS)。ON の最終場は `run_0521_twophase_dplur_nn0_cont` res_16000。

## 設計案
**D1 (G1 用、面ごとの流束ダンプ)**: 環境変数 `FORGE_DIAG_TP_FACES=<出力パス>` (既定 off) で、指定 step の組立の後に 1 回だけ、全通常面について次を書く (HDF5):
- 面の id・両端節点・面重心・幾何 (geo, geo_abs, 面積)、ct (μt,f/Sc_t)、D_k、両端の ρ・ρY・ρg・ρQ・T・P。
- **OFF の作用素 A**: species_diffusion_d の面の流束と同じ値 (J_s、エネルギー)。species_diffusion_d の面計算を `__device__` 関数に切り出し (本番の数値は不変)、診断カーネルから同じ関数を呼ぶ。
- **ON の作用素 B (float)**: tp_face_flux<float> の出力 (J, Jv, Jl, JQ, q) と、分子・乱流の成分を分けた値 (Jl は乱流のみ、蒸気の分子流束 = ∇z_v 駆動)。
- **double 参照**: 同じ格納入力・同じ係数で tp_face_flux<double> を評価した値 (監査 #4f と同じ経路)。差し引き前の項の大きさ (ρ_g,f D abs(z) など) も書く (G1 の絶対誤差尺度用)。
- 実行は ON の設定で起動し、A は診断だけで評価する (状態は ON の更新に使わない、0 step = nStepOuter 1 で出力のみ)。
G1 の判定 (事前登録済み) は Python で: (i) 乱流域で abs(J_w^A) ≤ 0.05·max abs(J_l^B)、(ii) B の float と double が絶対誤差尺度内。

**D2 (G3 用、閉じた収支)**: 環境変数 `FORGE_DIAG_TP_BUDGET=<出力パス>` で、組立後半の各段 (対流流束・二相拡散/化学種拡散・凝縮ソース・その他) の直後に成分別 (ρY_w, ρv = ρY_w − ρg, ρg, ρQn) の res を写して差を取り、段ごとの節点寄与を出す。あわせて更新の補正 (再正規化・射影・床・clamp) の成分別・節点別の量を同じ step で記録。Python で制御体積 (全域、x の帯、壁近傍帯) ごとに 境界流束 (移流+拡散)・相変化ソース (正負別)・残差・補正 を積算し、不整合 < 観測された液流束差の 10 % を判定。ON・OFF の両方の最終場で。

## 問い
1. D1 は G1 の事前登録 (codex plan M4: 同一面・同一格納入力・同一係数で独立な double 参照、絶対誤差尺度) を満たすか。tp_face_flux<double> は同じ式の倍精度評価で、式そのものの独立検証ではない (式は host ハーネス `tests/unit/test_twophase_kernel.cu` と Python 参照 `test_twophase_diffusion_harness.py` で照合済み)。これで足りるか、Python で別実装すべきか。
2. species_diffusion_d の面計算を関数に切り出すのは本番 (OFF) の数値を変えないか (float の評価順・FMA)。ビット不変をどう確かめるか (G0 と同じプロトコルで足りるか)。
3. D2 の段の切り方 (どの kernel の後で写すか) と、ρv の扱い (派生量なので res_roY_w − res_rog で作る)、補正の取得方法。閉じた収支の式 (定常の残差ゼロでない場で何が閉じるべきか) を明示してほしい。
4. 0 step (nStepOuter 1) で G1・G3 の両方を出してよいか。ON の最終場 (run_0521) に対しても同じ診断を回すべきか。
5. 実装の罠 (周期・境界半割面の skip、atomicAdd、ghost、node の境界節点)。

読んでよいファイル: 上記 plan 2 本、`methods/condensation.md` §7c、`solver_density_cuda/cuda_forge/speciesTransport_d.cu`・`twoPhaseDiffusion_d.cuh`・`condensationTransport_d.cu`・`condensationSourceKernels_d.cuh`、`solver_density_cuda/main.cpp`、`notes/investigations/2026-10-03-twophase-onoff/*.py`、`case/16.nozzle_wys/analyze_liquid_diffusion_error.py`。
