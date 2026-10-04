# 諮問ブリーフ: #1b A/B 4000 step の結果 — B の残差跳ね上がりと θ=0 の大量発生 (2026-10-02)

作業ツリー `/home/sano/work/forge-species` (HEAD 6a7a865d でビルド、AWS)。plan `plans/active/condensation-two-phase-transport.md` §5.1 #1b (事前登録 確定済み)・#4e/#4f/#1b-pre。設計メモ `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` §14–§17。AGENTS.md エスカレーション 3 (事前登録の受入未達)・4 (手順にない修正の前)。

## 観測事実 (run は AWS `~/forge-species-tp/case/16.nozzle_wys/`; 事前登録どおり作成・投入、共通初期場 run_0482 res_48000、A・B の初期 VALUE は全一致)

- `run_0505_twophase_ab_A` (旧作用素 + 監査) と `run_0506_twophase_ab_B` (二相拡散 ON + 監査)、各 4000 step、rc 0、NaN 0、各 ~5 分。
- A: check_convergence 全列 STALLED (rms_ro init 1.80e-7 → fin 9.51e-8, 0.3 dec; rms_roe 0.5 dec; rms_roK 1.1 dec falling)。renorm-gate FINAL FAIL (max|f−1| 1.27e-4 > κ 4.77e-7; 成分別 rhoYw 9.6e-8、液・Q 0)。監査 NOT CONVERGED (終了時 Q1/Q0 比 ~131)。
- B: 最初の 200 step で全残差が約 2 桁跳ね上がり以後 4000 step まで平坦 (rms_ro 1.80e-7 → 2.43e-5 @200 → 2.56e-5 @3999; rms_roUx 2.92e-5 → 1.48e-2; rms_roe 2.38 → 4.10; rms_roOmega 5.4e2 → 4.5e3 RISING; rms_roQ1 4.2e5 → 4.7e6; rms_roYv 9.7e-6 → 3.3e-6)。renorm-gate FINAL FAIL (max|f−1| 8.53e-3、成分別 rhoYw 9.9e-5・rhog 8.6e-4・Q 7〜8e-4)。監査 NOT CONVERGED (終了時 Q2/Q1/Q0 比 ~2.6e5〜3.0e5、開始時 1.2e4)。
- B の `[twophase]` 区間ログ: step 1 は θ<1 セル 0、step 3801/3901/4000 の区間で θ<1 セル 7.3〜7.7e4 (うち **θ=0 が 2.8〜2.9e4**)、min θ 0、保留量 蒸気 7e-8・液 6e-8、Qcut 4〜6e-6、vround ~3e-16。
- 監査の自己検査 (r_double vs r_float) は全成分 ≤0.92 ε·max A (監査は正しく同じ離散式を組んでいる)。

## ログ抜粋

```
[twophase] condTwoPhaseDiffusion 1: gas-phase molecular diffusion (z basis, upwind correction) + common turbulent mixing (species, liquid g, Q2/Q1/Q0, energy), unsplit vapour/liquid update with point-diagonal preconditioner, relax 1, condDgMaxStep 0.005, condD
[twophase] note: vapour, liquid and moments use the point-diagonal preconditioner (passiveImplicitCoupling / speciesImplicitCoupling DPLUR is not used for them); the other species keep their update
[twophase-audit] step 0 (initial state: r0): two-phase (key on) operator, double re-evaluation from the stored state (EOS via assembleResidual, face fluxes in double, source per cell; energy not covered); pinned nodes excluded 118
[twophase-audit]   roY0     max|r| 7.251e-06  max A 4.304e-02  r0 7.251e-06  ratio 2.356e+02 (<= 1)  max|r_double - r_float|/(eps max A) 0.659
[twophase-audit]   roY1(w)  max|r| 3.464e-07  max A 4.765e-04  r0 3.464e-07  ratio 1.016e+03 (<= 1)  max|r_double - r_float|/(eps max A) 0.653
[twophase-audit]   v        max|r| 7.578e-06  max A 4.759e-04  r0 7.578e-06  ratio 2.226e+04 (<= 1)  max|r_double - r_float|/(eps max A) 0.667
[twophase-audit]   g        max|r| 7.238e-06  max A 4.725e-04  r0 7.238e-06  ratio 2.142e+04 (<= 1)  max|r_double - r_float|/(eps max A) 0.566
[twophase-audit]   Q2       max|r| 5.628e-02  max A 6.555e+00  r0 5.628e-02  ratio 1.200e+04 (<= 1)  max|r_double - r_float|/(eps max A) 0.546
[twophase-audit]   Q1       max|r| 2.762e+06  max A 4.020e+08  r0 2.762e+06  ratio 9.603e+03 (<= 1)  max|r_double - r_float|/(eps max A) 0.719
[twophase-audit]   Q0       max|r| 2.251e+14  max A 2.624e+16  r0 2.251e+14  ratio 1.199e+04 (<= 1)  max|r_double - r_float|/(eps max A) 0.796
[cond-corr]   species 0 theta over all updates | interval: update theta<1 0 cell-updates in 1 updates (last update 0 cells) min 1 | source theta_src<1 0 cell-evals in 2 evals (last eval 0 cells) min 1 | cumulative: update theta<1 0 in 1 updates, theta_src<1 0 
[twophase] step 1 interval | updates 41640 | theta<1 cells 0 (theta=0 0) min theta 1.000000 | withheld vapour 0.000e+00 liquid 0.000e+00 | state corrections Qcut 2.305e-44 vround 0.000e+00
[twophase-audit]   Q2       max|r| 1.573e+00  max A 8.455e+00  r0 5.628e-02  ratio 2.602e+05 (<= 1)  max|r_double - r_float|/(eps max A) 0.733
[twophase-audit]   Q1       max|r| 1.135e+08  max A 5.351e+08  r0 2.762e+06  ratio 2.965e+05 (<= 1)  max|r_double - r_float|/(eps max A) 0.463
[twophase-audit]   Q0       max|r| 7.695e+15  max A 3.535e+16  r0 2.251e+14  ratio 3.044e+05 (<= 1)  max|r_double - r_float|/(eps max A) 0.429
[twophase-audit] VERDICT: NOT CONVERGED (component ratios max|r| / max(1e-7 r0, 6 eps max A); non-finite no)
[twophase-audit]   Q1       max|r| 3.779e+04  max A 4.021e+08  r0 1.861e+04  ratio 1.314e+02 (<= 1)  max|r_double - r_float|/(eps max A) 0.692
[twophase-audit]   Q0       max|r| 2.453e+12  max A 2.624e+16  r0 1.212e+12  ratio 1.307e+02 (<= 1)  max|r_double - r_float|/(eps max A) 0.924
[twophase-audit] VERDICT: NOT CONVERGED (component ratios max|r| / max(1e-7 r0, 6 eps max A); non-finite no)
```

## 問い (1 つ)

B のこの状態の真因の第 1 仮説と、それを判別する 1 つの A/B (事前に合否を書ける形) を示してほしい。事前登録どおり 8000 step の再実行 (`run_0507_twophase_ab_A8k`/`run_0508_twophase_ab_B8k`) を投入済みだが、平坦な 4000 step の後で延長が意味を持つかも判断してほしい。候補 (棄却してよい): (i) θ=0 で止まったセルが恒久的に更新を止め (#4b の S10 型)、原方程式残差が残ったまま擬似時間が回っている、(ii) 二相拡散の再正規化 (液・Q に係数を掛ける) と DPLUR 側の化学種更新の不整合で連続式と化学種が噛み合わない (max|f−1| 8.5e-3)、(iii) 点対角前処理が本番の CFL (cfl_pseudo) に対して不安定、(iv) 乱流 (ω RISING) との結合。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない (run は AWS 上で手元に無い)。推奨は 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- plan、設計メモ §14–§17、`solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh`、`speciesTransport_d.cu` (grep `twophase`・`renormalize`)、`condensationTransport_d.cu` (grep `twophase_vl_update`・`twoPhaseHoldWater`)、`main.cpp` (grep `twoPhase`)、`case/16.nozzle_wys/prepare_twophase_ab.sh`
