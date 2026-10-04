# 諮問ブリーフ: 凝縮二相拡散 #4c (実 condFloat ソース接続) の結果とカーネル着手の可否 (2026-10-02)

作業ツリー `/home/sano/work/forge-species`。plan: `plans/active/condensation-two-phase-transport.md` §5.1 #4〜#4c。前回諮問: `notes/reviews/2026-10-02-twophase-diffusion-4b-diagnose.md`。設計メモ: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` (§10 事前固定の合格条件、§11 結果)。ハーネス: `solver_density_cuda/tests/unit/test_twophase_real_source.cpp`。

**問い**:
- (Q1) §10 の事前固定条件で 4 variant とも合格した。比較条件どおり「ω = 1、2×2 は必須にしない」として、定常専用の初版カーネル (#4) に着手してよいか。
- (Q2) 固定条件では limiter も θ_src も作動しなかった。制限が効く領域の検証 (CFL を大きくする等) を、着手の前に要るとするか。
  感度試験では CFL ≥ 50 で ω = 1 が Q0 のリミットサイクルで上限到達し、ω = 0.5 だけが収束した。
- (Q3) 核生成の自己 Jacobian (sj_Q0 = 0、現行カーネルも同じ) を前処理に入れるべきか、緩和で扱うか。
- (Q4) 蒸気残差の停止判定を「独立に評価した残差」(float64 または同等の精度) で行うことを、カーネルの契約に書くか (S8 の f32 停止で v の比が 1.09〜1.12)。

**観測事実** (ホストのみ。CFD は無い):
- 実ソースの関数はすべて host/device で、ホストから呼べた。セルの組み立て (カーネル) は `__global__` なので、同じ順序の写しを使った。
- 固定条件 (float32、CFL 5、上限 20000):
  - 69 / 70 / 161 / 161 反復で合格した (ω=1 対角 / ω=1 2×2 / ω=0.5 対角 / ω=0.5 2×2)。
  - 独立残差の比 (v を含む) は ≤0.90。最終の増分は 0.29〜40 ULP (ω=0.5 では 0.5 ULP 未満が 93〜95 %)。
  - θ も θ_src も 1、補正 0、非負、不動点は同じ (出口 g 2.6033e-3)。
- 感度 (判定外):
  - CFL 50 / 500 で ω=1 は上限到達 (Q0 の比 7.8〜11.6、周期的、増分は Q0 で 110〜160 ULP)。ω=0.5 は 131〜134 反復で収束した。
  - 2×2 連成はどの条件でも差を作らなかった。
- S8: 監査に v を足すと、float32 停止は v の比 1.09〜1.12 で不合格。float64 再評価の停止なら合格。

**期待値と出典**: §10 の事前固定の合格条件と比較条件。diagnose 2026-10-02 (4b) の採否表。

**実施済みの操作**:
- C++ ハーネス `tests/unit/test_twophase_real_source.cpp` を新設した。
  - 入口状態は `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/res_48000.h5` の軸近傍の節点 18335 から取った。
  - 問題は ρ・U を固定した 1D ダクトで、気相は内蔵 N2 (run の MIXDRY の代わり) と H2O。
- Python ハーネスの Minor 修正 3 件 (§11.4)。
- カーネル・plan・methods は変更していない。

**仮説** (検証していない):
- H1: 固定条件で制限が作動しないのは、CFL 5 の擬似刻みでは 1 反復の Δg が dg_max (5e-3) と dT_max (1 K) を下回るからである。
- H2: CFL ≥ 50 のリミットサイクルは、核生成ソース (J の強い T・S 依存) の自己 Jacobian が前処理に無いことで起きる。
  CFD の局所 CFL (cfl_pseudo 4〜8 前後) がこの 1D の CFL のどこに当たるかは、未換算。
- H3: 蒸気・液の 2×2 連成は、この実ソースでは蒸気側の結合が弱いので効かない。
  #4b の S9 で効いたのは、人工の成長則の硬さ (Δt/τ = 2、ρg_ref 1e-6) のためである。

**呼び出し側の前提で疑わしいもの**:
- 1D の CFL を CFD の擬似時間刻みにどう対応させるか (ここでは Δτ = CFL·dx/U、音速を含めていない)。
- 制限が効かない条件での合格を、初版の着手条件として十分とみなすか。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は問いごとに 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` の §10〜§12、plan 全体、前回の諮問記録、`solver_density_cuda/tests/unit/test_twophase_real_source.cpp`、`cuda_forge/condensationSourceKernels_d.cuh:296-560`・`condensationSourceF_d.cuh`・`condensationUpdateLimiter_d.cuh`
