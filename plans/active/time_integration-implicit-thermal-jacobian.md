# block DPLUR のエネルギー行に熱伝導の Jacobian を入れる (`implicitThermalJacobian`)

## メタ

- **area**: `time_integration`
- **status**: `draft`
- **related_docs**:
  - `methods/time_integration/implementation.md` の「エネルギー行の熱伝導 Jacobian (`implicitThermalJacobian`)」(本 plan で追加)
  - `procedures/solver-settings.md` (キーの説明、実装時に追加)
- **related_plans**:
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) §5.1 #27 (発見元: 方向別の擬似 dt の発散の追究)
  - [`../accepted/time_integration-line-implicit-viscous-v2.md`](../accepted/time_integration-line-implicit-viscous-v2.md) (方向別の dt と 2026-09-03 の case/45 の発散)
  - [`../accepted/time_integration-general-eos-jacobian.md`](../accepted/time_integration-general-eos-jacobian.md) (TP の block DPLUR の Jacobian)
- **created**: `2026-10-09`
- **owner**: Claude (ユーザ指示 2026-10-09「カーネル修正やって」)

## 1. 目的

冷却ノズル (case/45、300 K 等温壁、y1+ < 1、近壁の縦横比 約 4000) で、ライン陰解法に方向別の擬似 dt (`lineDtDirectional`) を足すと
縮流部の壁から 1 層目で等圧のエントロピーのモードが育って発散する。単因子の A/B (§3) から、これは「固定温度の壁への熱伝導が駆動するモードを、
LHS の粘性の対角 (保存量 ρE にかかるスカラー) が抑えない」ためと見ている。エネルギー行を温度で組んだ熱伝導の Jacobian にして、
(1) この発散が止まるか、(2) 止まるなら方向別の dt で冷却壁の遅い過渡 (縮流部の質量の欠損、時定数は point の cfl 4 で約 13 万 step) が速くなるかを確かめる。

## 2. スコープ

- **やる**: `implicit_defect_correction_block_d` (ST = float / double) の内部の node 間面の粘性の対角のうち、エネルギー行だけを
  熱伝導の Jacobian (§4) に置き換える opt-in のキー `time.deltaT.implicitThermalJacobian` (既定 0 = ビット同一)。設定の検査、methods・procedures の記載、検証 (§6)。
- **やらない**: 近傍との熱伝導の結合 (非対角、ライン内の K)、局所 μ(T) の使用 (従来どおり `physProp.visc` 定数)、
  scalar 版 (`implicit_defect_correction_d`)・前処理版 (`lowMachPrecond>=2`)・cell 離散化、等温壁の拘束に合わせた Jacobian
  (Δ(ρE)_w = e_w·Δρ_w の畳み込み、diagnostician の第 1 仮説) — 本 plan の結果を見て別項目にする。既定化はしない (検証後にユーザ判断)。

## 3. 関連 docs と前提

発見の経緯と証拠は tooling-nozzle-isothermal-wall-chain §5.1 #27 (「方向別の dt の発散の追究」以降)。要点:

| run (case/45、run_0183 の res_100000 から、directional・cfl 4) | 変更 | 結果 |
| --- | --- | --- |
| run_0203 | — | 65 step で発散 (決定的) |
| run_0204 | `implicitSolvePrecision 1` | 65 step、壊れ方も同じ |
| run_0206 | `nStepInner 15` | 65 step |
| run_0207 | `convMethod 0` | 72 step |
| run_0210 | `FORGE_FREEZE_TURB=1` | 64 step |
| run_0211 | `lineViscCoupling 1` | 20 step (悪化) |
| run_0209 | 冷却壁を断熱壁に | 200 step まで有限、近壁の振動は減衰 (壁温 300 → 440 K の過渡を含む) |

- 帳簿 (run_0203、after_eos_bc): 1 層目は Δρ/ρ ≈ −ΔT/T、Δp/p はその 1/25 (等圧のエントロピーのモード)。1 層目のエネルギーの行の段別の残差は、
  熱伝導 (粘性の段) が育ち (step 8/16/24/32: +159/−197/+723/−1066)、対流の段は横ばい (±76〜148)。
- 諮問: codex 2 回 (`notes/reviews/2026-10-08-linedir-divergence-diagnose.md`、`notes/reviews/2026-10-09-linedir-divergence-2-diagnose.md`)、
  diagnostician 1 回 (第 2 仮説 = 粘性・熱伝導の陰的部分の不足、`timeIntegration_d.cu` 937〜951 行)。
- 現行の粘性の対角: `timeIntegration_d.cu` の block カーネル (2ν_eff·δ/dcc を `add_identity_scaled` で 5 行に同じ量)、ν_eff = (`physProp.visc` + μ_t)/ρ。

## 4. 設計方針

内部の node 間面 (`!(isNode && !has_nbr)`、従来の粘性の対角を課す面と同じ) で、`implicitThermalJacobian == 1` かつその面が
`lineViscCoupling` の対象でないとき:

- 運動量の行 (1〜3) と連続の行 (0) は従来どおり `viscous_diag = 2ν_eff·δ/dcc` を対角に足す。**行 0 にも従来どおり足す** (従来と同じ減衰を保つ。連続の式に粘性は無いが、従来の LHS の性質を変えない)。
- エネルギー行 (4) は `viscous_diag` を足さず、代わりに `diag_block[4][j] += Λ^T·∂e/∂Q_j` (j = 0..4):
  - Λ^T = γ_i·(μ_lam/Pr + μ_t/Pr_t)·δ/dcc。μ_lam = `laminar_visc` (= `physProp.visc`)、μ_t = `vis_turb[ic]`、Pr = `physProp.prandtlLam`、Pr_t = `turbulence.turbulentPrandtl` (カーネル引数を 3 つ足す: フラグ・Pr・Pr_t)。
  - e = ρE/ρ − ½|u|²、∂e/∂ρ = −(e − ½|u|²)/ρ、∂e/∂(ρu_k) = −u_k/ρ、∂e/∂(ρE) = 1/ρ (`roe[ic]`・`ro[ic]`・`Ux/Uy/Uz[ic]` から)。
  - 導出: 熱伝導の残差 k(T_j − T_i)δ/dcc の Q_i による微分の符号を反転したもの。k∂T/∂Q = (k/c_v)∂e/∂Q、k/c_v = γ(μ/Pr + μ_t/Pr_t) (凍結 γ、c_p = γc_v)。T・c_v・R を持たずに組め、TP の e の基準点によらない。
- 等温壁の節点 (行 4 は後段で単位行化される) と壁の節点の運動量の行は従来どおり上書きされる (変更なし)。
- 設定の検査 (`solverConfig.cpp`): 値は 0/1。1 のときは `timeIntegration 11`・`blockDPLUR 1`・`lowMachPrecond < 2` を要求し、`lineViscCoupling 1` との併用は拒否。
- 既定 0 ではコードの経路がビット同一になるよう、分岐は `if (thermalJac)` の中だけに置く。

**期待 (仮説)**: 1 層目の等圧のモード (Δρ ≠ 0、ΔT ≠ 0、Δ(ρE) は小) に対して、行 4 に ∂e/∂ρ の交差項が入るので ΔT に比例する減衰が働く。
1 層目の熱拡散数 (Δτ·κ/Δy²) は point の Δτ で約 0.01、方向別の Δτ で 0.5 を大きく超えるので、陽的に扱われていた熱伝導が陰的になる。

## 5. 実装ステップ

1. `solver_density_cuda/input/solverConfig.{hpp,cpp}`: キー `time.deltaT.implicitThermalJacobian` (int、既定 0) と検査。
2. `solver_density_cuda/cuda_forge/timeIntegration_d.cu`: block カーネルに引数 3 つ (フラグ、Pr、Pr_t) を足し、粘性の対角の分岐 (§4)。`FORGE_BDPLUR_ARGS` に追加。
3. `procedures/solver-settings.md`: キーの説明 (opt-in、効く範囲、併用不可)。
4. FP64 のビルド (AWS、`~/forge-wallfit-bin-fp64` と同じ作り方で新しい作業ツリー) と FP32 のビルド (コンパイルの確認)。
5. 検証 (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段のレビューと上位の諮問 | `codex_review.py --stage plan` と diagnostician に §4・§6 を諮る。Critical/Major の採否を §6.1 に | F |
| 2 | 実装 (§5 の 1〜3) | 合格: FP32・FP64 ともビルドが通る、既定 0 で V0 がビット同一 | O |
| 3 | 検証 V0〜V4 (§6) | 合格条件は §6 | O |
| 4 | result 段のレビューと採否、既定化の判断 (ユーザ) | | F |

## 6. 検証 (事前登録、2026-10-09)

ビルド: 本 plan の commit を FP64 (typedef double + 座標の stod、`cold_pair.py` の FORGE_SHA を新しいバイナリで登録し直す) でビルドする。

- **V0 (既定でビット同一)**: キー無し (= 0) で run_0203 の設定から directional を外した point の cfl 4 を 20 step、旧バイナリ (65be5e28…) と新バイナリで回し、
  res_20.h5 の全 `VALUE/*` がビット一致すること。
- **V1 (発散が止まるか、目的 1)**: run_0203 と同じ (run_0183 の res_100000、directional、cfl 4、nStepInner 5、relax 0.7、基準値固定、同じ 527 節点の帳簿) で `implicitThermalJacobian: 1` だけを足して 200 step。
  合格: 200 step まで有限、帯 (列 37〜49・1〜15 層) の 1 step ごとの ρ の増分の振動成分 (16 step の移動平均を引いた RMS) が最後の 32 step の窓 3 本で増えない。
  同じ位置・周期・成長なら仮説を外す。発散が遅れるだけなら保留。
- **V2 (速くなるか、目的 2; V1 が合格のときだけ)**: V1 の設定で 20000 step・5000 ごと (`extraFields: [res_ro, volume]`)。
  合格: 質量の欠損 (Σ res_ro×2π) が 5000/10000/20000 step で point の cfl 4 (run_0191: 1.917/1.820/1.660) より小さく、かつ壁時計あたりでも速い (ms/step を同じ負荷で記録)。
  NaN・DIVERGED なし。θ_r (x = 40/94) と δ_loc を記録 (定常解の比較には使わない)。
- **V3 (point での後退がないか)**: point の cfl 4 (run_0191 と同じ、run_0183 の res_100000 から) に `implicitThermalJacobian: 1` を足して 20000 step。
  合格: 欠損が run_0191 の同じ step の値より 2 % 以上悪くならない、NaN・DIVERGED なし、rms_roOmega が開始時の 10 倍を超えない。
- **V4 (断熱壁で後退がないか)**: run_0181 (断熱) の res_100000 から point・cfl 1 (run_0181 と同じ設定) に `implicitThermalJacobian: 1` を足して 5000 step。
  合格: NaN なし、全残差列の末尾 1000 step の中央値が開始時 (run_0181 の末尾) の 2 倍を超えない。
- 定常解: RHS を変えないので不動点は同じ (論証)。短い run では確かめられないので、V2・V3 の θ_r・δ_loc の差は記録のみ。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/timeIntegration_d.cu` (block カーネル)、`solver_density_cuda/input/solverConfig.{hpp,cpp}`
- 既定 0 のため既存のケース・実行手順への影響なし
- `methods/time_integration/implementation.md` (追加済み)、`procedures/solver-settings.md`

## 8. 完了条件

- [x] 関連 `methods/time_integration/` の現在仕様を更新済み
- [ ] 実装・検証完了 (本計画の §6 を満たす)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [ ] ファイルを `plans/active/` → `plans/accepted/` (superseded なら `archived/`) へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期 (移動元・移動先)

## 9. 変更ログ

- `2026-10-09` — 初稿 (ユーザ指示「カーネル修正やって」)。
