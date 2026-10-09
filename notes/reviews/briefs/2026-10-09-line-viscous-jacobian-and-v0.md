# 諮問: (1) ライン面の薄層の粘性・熱伝導 Jacobian の設計、(2) V0 の判定不能の扱い、(3) run_0224 の解釈

日付 2026-10-09。諮問先 codex (diagnose、`~/.config/forge/diagnose-backend` = codex)。作業ツリー `/home/sano/work/forge-integ-1005`。
エスカレーション条件: 1 (§4・§6 を新規に書く)、3 (事前登録の比較が判定不能)、6 (cuda_forge の数値の変更の前)、7 (result の解釈の前)。
plan: `plans/active/time_integration-line-viscous-jacobian.md` (新規、draft)、`plans/active/time_integration-implicit-thermal-jacobian.md` (§6・§6.0)。
ユーザ: 「K_ij に粘性の寄与は入れたほうがいいんじゃない?」「プラン書いて諮っていいよ」「まずは粘性込みのヤコビアンを試してほしい」。

## (1) 粘性・熱伝導 Jacobian (`lineViscCoupling: 2`)

### 観測事実

- ライン陰解法の K は `block_dplur_jacobian_d.cuh` の `accumulate_split_jacobian_cf` の −A⁻ (固有値 λ₁ = V + c、λ₂ = V (重複 3)、λ₅ = V − c の固有分解、V = 面の法線速度) を単位ベクトルで列ごとに取り出したもの
  (`timeIntegration_d.cu` の `isLineFace && storeLU`)。粘性はライン面でもライン外の面でも対角のスカラー 2ν_eff δ/dcc (全 5 行) だけで、隣との結合は入っていない (`lineViscCoupling 0`)。
- case/45 (300 K の冷却ノズル、y1+ < 1、近壁の縦横比 約 4000、ライン 4719 本 × 121 節点) で、ライン + 方向別の擬似 dt (`lineDtDirectional`) は 65 step で発散した (run_0203)。
  キー 5 (エネルギー行に熱伝導の T の Jacobian を対角だけ) + 上限 R = 50 (Δτ ≤ 50 Δτ_point) で安定 (run_0223/0224)。
  1 層目の拡散数 νΔτ/Δn² は point で約 0.006、方向別で約 23。FREEZE_TURB の A/B (run_0222) で SST と平均流の両方が不安定に寄与。
- `lineViscCoupling 1` (スカラー α I を 5 行すべて、対角 2α → α) + 方向別は 20 step で発散 (run_0211、キー 0)。
- 残差 (`viscousFlux_d.cu`、node、`heatCorrSU2 0`): τ = μ_f δ/dcc (u₁ − u₀) + 面平均の勾配の項、q = k_f δ/dcc (T₁ − T₀) + 勾配の項、
  `res_ρu[ic0] += τ`、`res_ρE[ic0] += τ·u_f + q` (u_f = f u₀ + (1−f) u₁)、ic1 は逆符号 (atomicAdd)。

### 呼び出し側の設計 (plan §4.1)

ライン面 f (自節点 i、ライン上の隣 j) で β = μ_f δ/dcc、κ = k_f δ/dcc、P = I + ⅓ n̂ n̂ᵀ、f_i、Δu = u_j − u_i、ū = f_i u_i + (1−f_i) u_j:
連続の行 0、運動量の行 D_i += β P ∂u_i/∂Q_i・K_ij += β P ∂u_j/∂Q_j、エネルギーの行 D_i += κ ∂T_i/∂Q_i + β (P ū − f_i P Δu)ᵀ ∂u_i/∂Q_i・
K_ij += κ ∂T_j/∂Q_j + β (P ū + (1−f_i) P Δu)ᵀ ∂u_j/∂Q_j。符号系は D ΔQ_k − K_prev ΔQ_{k−1} − K_next ΔQ_{k+1} = rhs、D = V/Δτ − ∂R_i/∂Q_i、K = ∂R_i/∂Q_j。
隣 j が壁 (速度 Dirichlet) なら運動量・仕事の K を 0、等温壁なら熱伝導の K を 0。ライン面のスカラー 2ν δ/dcc は入れない。ライン外の面は従来どおり。

### 問い

1. 式・符号・壁の扱いに誤りはないか。P (法線の 4/3) と粘性の仕事を入れるのは妥当か (残差の勾配の項は入れない)。
2. 連続の行がライン面の人工の減衰 (従来 2α) を失うことの危険。キー 1 で起きた「行ごとの減衰の基準の食い違い」(温度を変えない補正で行 4 だけ減衰が抜け、高速域に偽の ΔT) と同じ種類の問題は起きないか。
3. run_0211 (`lineViscCoupling 1`) の発散の見立て (連続の行への拡散・ρE での結合・対角の半減) は妥当か。値 2 でも同じ壊れ方をする可能性は。
4. 検証 (plan §6: U0・U2 純伝導・V-n1 上限なし 2000 step・V-n2 run_0223 と同条件の 15000 step) の判定に穴はないか。SST が分離したままで上限なしを試す意味はあるか。

## (2) V0 (thermal-jacobian plan §6、事前登録 = 既定で旧・新バイナリの res_20 がビット一致)

- 結果: 4 組 (point / directional × `implicitSolvePrecision` 0/1) ともビット不一致 (run_0225〜0232、`_band_ab/cold_pair/V0_bitident.json`)。
- 事後の物差し (run_0233〜0248、`v0_repeat.sh`・`v0_compare.py` → `V0_repeat.json`): 同じバイナリの再実行でも 20 step で 0/3 組・1 step で 0/1 組がビット不一致
  (粘性・対流の残差を atomicAdd で足す)。ρ の場の差の相対 RMS: point ST0 20 step で旧×旧 3.3〜3.5e-9・新×新 3.0〜3.6e-9・旧×新 3.1〜3.6e-9、
  directional ST0 20 step で 1.0〜1.2e-7・0.85〜1.0e-7・0.90〜1.2e-7、point 1 step で 8.6e-13・8.6e-13・2.2〜8.4e-13、directional 1 step で 1.4e-14・1.6e-14・1.6〜1.9e-14。
- 問い: 判定不能の扱い。代わりの基準 (例: 旧×新の差が同じバイナリの再実行の差の分布に入る、または atomic を決定的にするデバッグ経路でビット比較) を事後に置くことの是非と、置くならどの形か
  (「旧の範囲に入る」規則は同分布でも高率で落ちるという過去の教訓がある)。V2・V3 の run_0191 との比較を進めてよい条件。

## (3) run_0224 の解釈 (thermal-jacobian plan §6.0)

- run_0223 (directional + キー 5 + 上限 50、15000 step) → run_0224 (同設定、60000 step)。通算 75000 step で、θ_r (x = 40/70/94) 0.07761 / 0.11453 / 0.12727
  は run_0217 (point cfl 4、通算 64 万 step、事前登録の水準に入っている) の −0.16 / −0.21 / −0.24 %、Q_w は +0.03 %、欠損 0.0206 kg/s (point 0.0089)。
  θ_r(70) の 5000 step ごとの増分は比 0.66〜0.73 で減衰 (等比外挿 0.1149)、末尾 2 万 step のドリフト +1.4 % (水準 0.05 % は未達)。
- check_convergence (0〜59999): NOT CONVERGED、rms_ro (3.66e-6 → 1.00e-5、最大 1.47e-5) と rms_roUy (3.85e-4 → 1.29e-3) が RISING、k・ω は下降。場所は未確認。
- 事前登録どおり延長中 (run_0252、同設定 60000 step)。水準に達したら切り戻し (同じ新バイナリでキー 0・上限 0 の point cfl 4、2〜4 万 step)。
- 問い: 残差の RISING をこの段階でどう扱うべきか (延長を続けてよいか、場所の特定を先にするか)。point の水準との一致を「同じ不動点」と言うための条件に追加すべきものは。

## 読んでよいもの

- 上記 2 つの plan、`plans/active/time_integration-line-implicit-speed.md` (同じカーネルの速度、参考)
- `solver_density_cuda/cuda_forge/timeIntegration_d.cu` (`implicit_defect_correction_block_d` の粘性の対角とライン K、`lineThomas*`)、`block_dplur_jacobian_d.cuh`、`viscousFlux_d.cu`
- `methods/time_integration/implementation.md` (line-implicit・v2・熱伝導 Jacobian)、`case/45.isobutane_m6_d155/README.md` (run_0203〜0253 の行)
- 過去の諮問: `notes/reviews/2026-10-09-thermjac-key5-diagnose.md`、`notes/reviews/2026-10-09-linedir-divergence-2-diagnose.md`
