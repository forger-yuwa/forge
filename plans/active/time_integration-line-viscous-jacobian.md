# ライン陰解法の K と対角に薄層の粘性・熱伝導の Jacobian を入れる (`lineViscCoupling: 2`)

## メタ

- **area**: `time_integration`
- **status**: `draft`
- **related_docs**:
  - `methods/time_integration/implementation.md` の「line-implicit」「v2」(`lineViscCoupling`) と「エネルギー行の熱伝導 Jacobian」
- **related_plans**:
  - [`time_integration-implicit-thermal-jacobian.md`](time_integration-implicit-thermal-jacobian.md) (キー 5 と上限 R、§6.0 の run_0223/0224)
  - [`time_integration-line-implicit-speed.md`](time_integration-line-implicit-speed.md) (同じカーネルの速度の改善)
  - [`../accepted/time_integration-line-implicit-viscous-v2.md`](../accepted/time_integration-line-implicit-viscous-v2.md) (`lineViscCoupling: 1` = スカラーの結合)
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) §5.1 #27 (冷却壁の遅い過渡)
- **created**: `2026-10-09`
- **owner**: Claude (ユーザ指示 2026-10-09「K_ij に粘性の寄与は入れたほうがいいんじゃない」「プラン書いて諮っていいよ」「まずは粘性込みのヤコビアンを試してほしい」)

## 1. 目的

ライン陰解法の近傍行列 K は対流流束の固有値分割の負の側 (−A⁻) だけで、固有値は音波の V ± c とせん断・エントロピーの V (V = 面の法線速度) である。
壁法線のラインでは壁際で V ≈ 0 なので、K が結合するのは音波 (圧力と法線速度) だけで、**せん断 (接線の速度) とエントロピー (等圧での T・ρ) のモードはライン内で結合していない**。
これらを壁法線につなぐのは物理的には粘性と熱伝導だけで、LHS では対角のスカラー 2ν_f δ_f/dcc_f しか入っていない。
冷却ノズル (case/45) の遅い過渡は近壁の T・ρ と境界層の速度なので、薄層の粘性・熱伝導の Jacobian をライン面の K と対角に入れて、
(1) 方向別の擬似 dt の上限 R を外しても安定か、(2) 遅い過渡の step 数がさらに縮むかを確かめる。

## 2. スコープ

- **やる**: ライン面 (Thomas が厳密に解く面) の粘性・熱伝導の薄層 Jacobian (`lineViscCoupling: 2`)、純伝導の試験 (U2)、case/45 の A/B。
- **やらない**: ライン外の面 (従来のスカラーのまま)、SST の k・ω のライン化 (別議題)、Thomas の速度 ([time_integration-line-implicit-speed](time_integration-line-implicit-speed.md))。

## 3. 関連 docs と前提

- 残差 (`viscousFlux_d.cu` の内部面、node・既定 `heatCorrSU2 0`): 面 f (ic0 → ic1) で
  τ = μ_f δ/dcc (u₁ − u₀) + (面平均の勾配による転置・発散・非直交の項)、q = k_f δ/dcc (T₁ − T₀) + (勾配の項)、
  `res_ρu[ic0] += τ`、`res_ρE[ic0] += τ·u_f + q` (u_f = f u₀ + (1 − f) u₁)、ic1 は符号を反転。μ_f = f μ₀ + (1−f) μ₁ (層流 + μ_t)、k_f = `tc_face` (層流 k + c_p μ_t/Pr_t の面の値)。
- LHS の形 (`implicit_defect_correction_block_d` と `lineThomas*`): 節点 i の行は D_i ΔQ_i − Σ_{line j} K_ij ΔQ_j = rhs_i、D = V/Δτ + Σ A⁺ + (粘性の対角)、K = −A⁻ (+ `lineViscCoupling 1` の α I)。
  拡散の残差 R_i = α (Q_j − Q_i) なら D += α、K += α の符号系。
- `lineViscCoupling 1` (スカラー α I を 5 行すべて、対角 2α → α) は case/45 の directional で 20 step で発散した (run_0211)。連続の行にも ρ の拡散を入れ、エネルギーの行を ρE で結合する点が物理の Jacobian と違う (仮説、確かめていない)。

## 4. 設計方針 (案。codex の plan 段と諮問の前)

### 4.1 式 (ライン面 f、自節点 i・ライン上の隣 j、`isLineFace` のとき)

β = μ_f δ/dcc、κ = k_f δ/dcc (残差と同じ面の値)、n̂ = 面の単位法線、P = I + ⅓ n̂ n̂ᵀ (薄層の応力の法線成分 4/3)、f_i = 自節点の補間の重み (ic0 なら f、ic1 なら 1 − f)、
Δu = u_j − u_i、ū = f_i u_i + (1 − f_i) u_j。保存量 Q = (ρ, ρu, ρv, ρw, ρE) について

  ∂u/∂Q = (1/ρ) [−u, I₃, 0]、 ∂T/∂Q = (γ/c_p)(1/ρ) [−(e − ½|u|²), −uᵀ, 1] (e = ρE/ρ − ½|u|²、キー 1 と同じ)。

- 連続の行: 0 (粘性の流束がない)。
- 運動量の行: D_i += β P ∂u_i/∂Q_i、 K_ij += β P ∂u_j/∂Q_j。
- エネルギーの行: D_i += κ ∂T_i/∂Q_i + β (P ū − f_i P Δu)ᵀ ∂u_i/∂Q_i、 K_ij += κ ∂T_j/∂Q_j + β (P ū + (1 − f_i) P Δu)ᵀ ∂u_j/∂Q_j
  (粘性の仕事 τ·ū、τ = β P Δu の微分)。
- 隣 j が壁の節点 (速度の Dirichlet) なら運動量の項と仕事の項の K を 0、等温壁の節点なら熱伝導の項の K を 0 (境界条件で値が固定され Δu_j = 0・ΔT_j = 0)。
  自節点が壁の節点のときは従来どおり `rowDec` の行が単位行になる。
- ライン面では従来のスカラー 2ν δ/dcc (全行) を入れない (上の式で置き換える)。ライン外の面は従来どおり (キー `implicitThermalJacobian` のビットもライン外の面にだけ効く)。
- 物性 (μ、k、c_p、γ) は凍結する。勾配の項 (転置・発散・非直交) は入れない (薄層近似)。

### 4.2 予想される危険

- 連続の行はライン面の粘性の減衰を失う (物理として正しいが、従来は 2α の人工の減衰があった)。V/Δτ が小さい方向別では ρ の行は対流 (音波) の結合と A⁺ だけになる。
- k・ω は従来どおり分離した点陰解法で同じ Δτ を使うので、上限 R を外すと SST 側の不安定 (run_0222 で寄与を確認) は残る可能性がある。
- 勾配の項を入れないので、非直交が強い面や軸対称の τ_θθ では LHS と残差の差が残る。

### 4.3 実装

`timeIntegration_d.cu` の粘性の対角のブロック (`isLineFace && lineViscCoupling != 0` の分岐) に値 2 の経路を足す。K は `storeLU` の sweep で `Kprev`/`Knext` に書き、
対角は D に足す。設定の検査 (`solverConfig.cpp`): 値 0〜2、2 は node・`lineImplicit 1`・`blockDPLUR 1`・`timeIntegration 11`・`lowMachPrecond < 2` を要求、
`implicitThermalJacobian` との併用は許す (値 1 の併用拒否は残す)。

## 5. 実装ステップ

1. `solverConfig.cpp` の検査と `methods/time_integration/implementation.md` の追記。
2. カーネル (§4.3)。FP64 ビルド (AWS) と FP32 ビルドのコンパイル。
3. §6 の検証。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段と諮問 | §4.1 の式・符号・壁の扱い、§6 の判定 | F |
| 2 | 実装 (§4.3) | 合格: FP64・FP32 のビルドが通る、`lineViscCoupling 0/1` の経路は不変 (U0) | O |
| 3 | U2・V-n1・V-n2 (§6) | | O |
| 4 | result 段のレビューと採否 | | F |

## 6. 検証 (事前登録の案)

- **U0 (既定で不変)**: `lineViscCoupling 0` で旧・新バイナリの 1 step・20 step の場の差が、同じバイナリの再実行の差 (V0_repeat.json) と同じ幅。
- **U2 (純伝導、Jacobian が効くか)**: case/52 の流体の板 (U1 と同じ、下 300 K・上 350 K、静止) をライン陰解法 + 方向別で回す。
  比較: (a) `lineViscCoupling 0`、(b) 2、ともに cfl 5 と 50。合格: (b) が L∞ < 0.05 K に入る step 数が point キー 0 cfl 20 の最良 (2500 step) の 1/10 以下、NaN なし。
  (a) より少ない step で入ること。
- **V-n1 (元の発散が止まるか)**: run_0183 の res_100000 から directional cfl 4・上限なし・キー 5 + `lineViscCoupling 2`、2000 step・200 ごと。
  合格: NaN なし、全残差列の最大が開始の 10 倍以内 (run_0203 は 66 step で発散、run_0221 は 580 倍)。
- **V-n2 (速さ)**: run_0223 と同じ (上限 50・キー 5、15000 step・500 ごと) + `lineViscCoupling 2`。合格: 欠損が run_0223 の同じ step 以下 (7500・15000)、
  θ_r(70) が run_0223 より point の水準 (run_0217、0.11477) に近い、ms/step の増分 10 % 以内。
- V-n1 が合格したら上限なしで 15000 step (V-n3) を回し、V-n2 と同じ量で比べる。最終の判断 (本線に使うか) は thermal-jacobian plan §6.0 と同じ条件 (状態の水準・切り戻し) で行う。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/timeIntegration_d.cu`、`solver_density_cuda/input/solverConfig.cpp`。`lineViscCoupling 0/1` の経路は変えない。

## 変更ログ

- 2026-10-09: 起票 (draft)。
