# 諮問ブリーフ: 入口列のリミッタの揺れ — cfl A/B と面単位検査の解釈、次の A/B

エスカレーション条件 4 (承認済み手順に無い次の run へ進む前) と 7 (解釈の確定前)。plan: `plans/active/limiter-inlet-column-oscillation.md` §5.1 #2r・#4r・#5。

## 観測事実
- #2r: cfl_pseudo 2→1 (参照値固定・同一 restart・48000 step) で床 0.81 倍、入口残差ノルム 0.83 倍、ψ の揺れ同程度。§4.3 の判定 B。報告量は不変。
- #4r: `notes/investigations/2026-10-03-limiter-inlet/face_probe.py` でソルバ式 (`limiter_d.cu:344-384`, scaled=1, node, convMethod 1) を再計算し、Ux/Uy/P は出力と 2e-7 以内で一致。ρ は一致しない (原因未確認)。
  入口境界列の壁近傍で、制限面の δ± がしばしば厳密に 0、P の δ± の中央値は float32 の 2.5〜22 ulp、ε̂ ≈ 5e-8 ≈ P の 1 ulp/q_ref。制限面の入れ替わり率 0.85〜1.0/step。cfl 2・1、開始・終了の 4 窓で同じ。
- 前の諮問: `notes/reviews/2026-10-03-limiter-inlet-column-plan-diagnose.md` (第 2 仮説 = H-R)。

## 提案する次の A/B (事前登録前)
(a) `venkatK` 0.05 → 0.2 (ε̂ = (K h/L)^1.5 が 8 倍)。config のみ。
(b) 全域 FP64 ビルド (同 commit 8ce8da3a、`flowFormat.hpp` の typedef だけ変更、AWS で native build) で K 0.05・cfl 2 のまま。
共通: IC run_0524 res_48000 (FP64 は --keep-src-dtype ではなく float の場から)、参照 4 値固定、48000 step、開始・終了の毎 step 窓、§4.3 と同じ判定 (R_入口・全域ノルム・ψ の揺れ)。

## 問い
1. #4r の観測から「H-R と整合」と書いてよいか。ρ の不一致は解釈に影響するか (`limiter_d.cu` の prim 経路 :325-334 で ρ スロットに何が入るか確認してほしい)。
2. (a)(b) は H-R の因果を判別できるか。H-R が真なら (b) で揺れと床が消えるはずだが、FP64 化は流束・勾配・EOS も変えるので交絡しないか。より絞った試験 (例: リミッタ関数だけ double) が config/ビルドで可能か。
3. 判定基準の数値 (主因 = 1/5 以下、否定 = 2 倍以内) を (a)(b) にもそのまま使ってよいか。
4. 揺れが残差の床を作っている (因果の向き) をどう確かめるか。
読んでよいファイル: 上記 plan と記録、`solver_density_cuda/cuda_forge/limiter_d.cu`、`limiterFunctions_d.cuh`、`flowFormat.hpp`、`notes/investigations/2026-10-03-limiter-inlet/*.py`。
