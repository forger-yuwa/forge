# block DPLUR のエネルギー行に熱伝導の Jacobian を入れる (`implicitThermalJacobian`)

## メタ

- **area**: `time_integration`
- **status**: `in_progress`
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
  熱伝導の Jacobian (§4) に置き換える opt-in のキー `time.deltaT.implicitThermalJacobian` (ビットマスク、既定 0 = ビット同一)。
  ビット 1 = 熱伝導の Jacobian、ビット 2 = 等温壁の節点のエネルギー行を拘束の行 [−e_w, 0, 0, 0, 1] にする (第 1 仮説の検査用、§4.2)。設定の検査、methods・procedures の記載、検証 (§6)。
- **やらない**: 近傍との熱伝導の結合 (非対角、ライン内の K)、局所 μ(T) の使用 (従来どおり `physProp.visc` 定数)、
  scalar 版 (`implicit_defect_correction_d`)・前処理版 (`lowMachPrecond>=2`)・cell 離散化、k の温度微分。既定化はしない (検証後にユーザ判断)。

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

### 4.1 熱伝導の Jacobian (ビット 1)

内部の node 間面 (`isNode && has_nbr`。従来の粘性の対角を課す面のうち node 間のもの。cell は設定で拒否) で、
`implicitThermalJacobian & 1` かつその面が `lineViscCoupling` の対象でないとき:

- 連続の行 (0) と運動量の行 (1〜3) は従来どおり `viscous_diag = 2ν_eff·δ/dcc` を対角に足す (従来の LHS の性質を変えない)。
- エネルギー行 (4) は `viscous_diag` を足さず、代わりに `diag_block[4][j] += Λ^T·(γ_i/c_p,i)·∂e/∂Q_j` (j = 0..4):
  - **Λ^T = k_face·δ/dcc**、k_face は残差 (`viscousFlux_d.cu` 264〜266 行) と同じ式: `f·thermCond[ic0] + (1−f)·thermCond[ic1] + (f·cp[ic0] + (1−f)·cp[ic1])·(f·μt[ic0] + (1−f)·μt[ic1])/Pr_t`、f = `fx[ip]`。
    物性は凍結 (k の温度微分は入れない)。熱伝導の非直交補正は従来どおり lag (残差だけ)。
  - c_p,i = `cp[ic]`、γ_i = `gamma_arr[ic]` (c_v = c_p/γ、凍結 γ)。
  - e = ρE/ρ − ½|u|²、∂e/∂ρ = −(e − ½|u|²)/ρ、∂e/∂(ρu_k) = −u_k/ρ、∂e/∂(ρE) = 1/ρ (`roe`・`ro`・`Ux/Uy/Uz`)。
  - 導出: 熱伝導の残差 k_face(T_j − T_i)δ/dcc の Q_i による微分の符号を反転したもの。∂T/∂Q = ∂e/∂Q / c_v は TP でも厳密
    (`dependentVariables_d.cu` 104〜112 行が Y を正規化してから T を反転するので、ρY 凍結で ρ が動いても Y は不変 — diagnostician 2026-10-09)。
    e の基準点によらず Λ^T·Δe/c_v = Λ^T·ΔT になる。
- 壁の節点 (node の壁 DOF) と内部節点の間の面も `has_nbr` なので、壁への熱伝導は 1 層目の対角に入る (追加不要)。
  等温壁の節点自身の行 4 は後段 (`timeIntegration_d.cu` 1055 行〜) で単位行に上書きされる (変更なし)。
- 設定の検査 (`solverConfig.cpp`): 値は 0〜3。≠0 のときは `discretization node`・`timeIntegration 11`・`blockDPLUR 1`・`lowMachPrecond < 2` を要求し、
  `lineViscCoupling 1` と `nodeIsothermalEnergyBC 1` との併用は拒否 (後者は weakIsoDiag が (4,4) に CPG の c_v のスカラーを足すので行の形と混ざる)。
- 既定 0 ではコードの経路をビット同一にする: 新しい配列の読み・一時変数は `if (thermalJac & 1)` の中だけに置き、従来の `add_identity_scaled` の経路は残す。

**期待 (仮説)**: 1 層目の等圧のエントロピーのモード (ΔT ≠ 0) に対して、行 4 が Λ^T·ΔT に比例する対角の減衰を持つ。
対角だけの陰的拡散は減衰 Jacobi と同型で、V/Δτ → 0 の極限でも r = (RHS の係数)/(LHS の係数) < 1/ω (ω = implicitRelax 0.7、1/ω ≈ 1.43) なら安定 (diagnostician)。
LHS の係数を残差と同じ k_face にするのはこのため。

### 4.2 等温壁の拘束の行 (ビット 2、第 1 仮説の検査用)

`implicitThermalJacobian & 2` のとき、等温壁の節点の行 4 を単位行 `[0,0,0,0,1]`・rhs 0 から拘束の行 `[−e_w, 0, 0, 0, 1]`・rhs 0 に替える
(e_w = roe[ic]/ro[ic]、壁は u = 0)。線形系の解が dq_roe_w = e_w·dq_ro_w となり、後段の壁温のピン (ρe = ρ·e(T_w)) の再構成と一致する。
`rowDec[4]` (K の行 4 のゼロ化) はそのまま。

## 5. 実装ステップ

1. `solver_density_cuda/input/solverConfig.{hpp,cpp}`: キー `time.deltaT.implicitThermalJacobian` (int、既定 0) と検査。
2. `solver_density_cuda/cuda_forge/timeIntegration_d.cu`: block カーネルに引数 3 つ (フラグ、Pr、Pr_t) を足し、粘性の対角の分岐 (§4)。`FORGE_BDPLUR_ARGS` に追加。
3. `procedures/solver-settings.md`: キーの説明 (opt-in、効く範囲、併用不可)。
4. FP64 のビルド (AWS、`~/forge-wallfit-bin-fp64` と同じ作り方で新しい作業ツリー) と FP32 のビルド (コンパイルの確認)。
5. 検証 (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段のレビューと上位の諮問 | 判断: 2026-10-09 codex GO-with-changes (C0/M4/m1) と diagnostician (Major 2・Minor 3) — 採否は §6.1、§4・§6 に反映済み | F |
| 2 | 実装 (§5 の 1〜3、§4.1・§4.2) | 合格: FP32・FP64 ともビルドが通る、`-Xptxas -v` のレジスタ・スピルを前後で記録、既定 0 で V0 がビット同一 | O |
| 3 | 検証 U1・V0〜V4 (§6) | 合格条件は §6 | O |
| 4 | result 段のレビューと採否、既定化の判断 (ユーザ) | | F |

## 6. 検証 (事前登録、2026-10-09、レビュー反映後)

ビルド: 本 plan の commit を FP64 (typedef double。座標の stod は既に HEAD にある) でビルドする。旧バイナリ `~/forge-wallfit-bin-fp64` (sha256 65be5e28…) は
e2696d8f0 + 同じ 2 つの手の変更で、HEAD との solver の差は `gmshReader.hpp` の stod だけ (同じ内容) — 確認済み。FP32 もビルドする。
`cold_pair.py` のバイナリの照合は新しいバイナリを登録し直す。ST (= `implicitSolvePrecision`) と FP64/FP32 のビルドは別の軸として表に分けて書く。

- **U1 (Jacobian そのものの検査、純伝導)**: case/52 のテンプレートの流体の板 (一様格子・低圧) を、上下とも等温壁 (300 K / 350 K、CPG、静止、node、陰解法) にして、
  初期 T 300 K 一様から回す。キー 0 と 1 を、cfl_pseudo 5 (memory の安定値) で比べ、さらにキー 1 は cfl_pseudo 20・50 も回す。
  量: 壁の熱流束 q の解析解 kΔT/H からの差が 0.5 % 以内に入る step 数と壁時計、NaN、線形解の失敗の件数。
  合格: キー 1 が同じ cfl で発散せず、0.5 % に入る step 数がキー 0 以下。キー 1 で cfl を上げて同じ q に少ない step で入れば Jacobian が効いている。
  キー 1 だけ発散すれば符号・大きさの誤り。FP32 のビルド (ST float) でもキー 1 を 1 本回し、NaN と線形解の失敗がないこと。
- **V0 (既定でビット同一)**: キー無しで (a) run_0203 の設定そのもの (ライン・方向別・cfl 4) と (b) point の cfl 4 を 20 step、
  `implicitSolvePrecision` 0 と 1 の両方で、旧バイナリと新バイナリの res_20.h5 の全 `VALUE/*` がビット一致すること。**V2・V3 を run_0191 と比べるのは V0 が合格のときだけ**。
- **V1 (発散が止まるか、目的 1)**: run_0203 と同じ (run_0183 の res_100000、directional、cfl 4、nStepInner 5、relax 0.7、基準値固定、同じ 527 節点の帳簿) で `implicitThermalJacobian: 1` だけを足して 200 step。
  合格: 200 step まで有限、線形解の失敗 0 件、帯 (列 37〜49・1〜15 層) の 1 step ごとの ρ の増分の振動成分の RMS が最後の 32 step の窓 3 本で増えない、
  かつ全残差列と物理量 (ρ・T・P > 0、T ≤ 全温) が悪化していない。帳簿の 1 層目のエネルギー行の粘性の段の成長が消えていることも確認する (消えずに止まったなら別の減衰)。
  同じ位置・周期・成長なら **V1b** (キー 2 = 拘束の行だけ、同じ条件) を回す: V1b で止まれば第 1 仮説 (壁の拘束)、V1b でも同じなら両仮説を外し Δτ の伸びの上限 (方向別 dt 側) に戻る。
  発散が遅れるだけなら保留。なお、今回の対角だけの近似で止まらなくても「熱伝導の陰的扱い全般」は棄却しない (codex m5)。
- **V1c (内訳、修正と独立)**: run_0203 と同じ設定 (キー 0) を `FORGE_WI_FORCE_DIAG=1` で 70 step 回し、1 層目のエネルギーの粘性の段を熱伝導 (`wi_eheat`) と粘性仕事 (`wi_ework`) に分けて記録する (codex m5)。
- **V2 (速くなるか、目的 2; V1 が合格のときだけ)**: V1 の設定で 20000 step・5000 ごと (`extraFields: [res_ro, volume]`)。
  量: 欠損の大きさ |2πΣres_ro| / ṁ_in (ṁ_in は数値流束の入口流量) と局所の残差のノルム (RMS)。
  **主判定は V3 (同じ LHS の point) との比較**: 同じ欠損の水準 (V3 の 20000 step の値) に達する壁時計が V3 より短いこと (ライン陰解法は 1 step 約 2 倍なので、step 数でなく時間で比べる)。
  双方が 20000 step 内にその水準へ達しないときは、共通の経過時間での値で比べ、結論は「過渡の改善」に限る。run_0191 (旧 LHS の point) は参考。
- **V3 (point での後退がないか)**: point の cfl 4 (run_0191 と同じ、run_0183 の res_100000 から) に `implicitThermalJacobian: 1` を足して 20000 step。
  合格: |欠損|/ṁ_in が run_0191 の同じ step の値より 2 % 以上悪くならない、NaN・DIVERGED なし、rms_roOmega が開始時の 10 倍を超えない。
- **V4 (断熱壁で後退がないか)**: run_0181 (断熱) の res_100000 から point・cfl 1 (run_0181 と同じ設定) に `implicitThermalJacobian: 1` を足して 5000 step。
  合格: NaN なし、全残差列の末尾 1000 step の中央値が開始時 (run_0181 の末尾) の 2 倍を超えない。
- **判定の道具**: V1〜V4 は `check_convergence.py --segment` の VERDICT と、V2・V3 は 5 枚以上の出力があれば `check_quasisteady.py` (欠損・θ_r) も保存する。
  短い試験の `NOT CONVERGED`・`DRIFTING` は許すが、「安定化・過渡の加速」と「収束・定常解の一致」は別の判定として書く。
- 定常解: RHS を変えないので不動点は同じ (根の維持の論証であって、根への到達の保証ではない — codex)。V2・V3 の θ_r・δ_loc の差は記録のみ。

### 6.0 結果と追加の事前登録 (2026-10-09)

- **U1 (純伝導)**: キー 0 / 1 の 0.05 K に入る step = cfl 5: 8500 / 9500、cfl 20: 2500 / 3500。cfl 50 はキー 1・キー 0 とも負の P・ρ で壊れる (NaN にならない)。
  → 事前の基準 (同じ cfl で step 数がキー 0 以下) は**不合格**、符号・大きさの誤りの基準 (キー 1 だけ壊れる) には当たらない。
  **ただし U1 は目的に合っていなかった** (ユーザ指摘 2026-10-09「このシンプルな試験がノズル CFD の状況をきれいに表しているとは限らない」):
  1 step の熱拡散数 D = Δτ·κ/Δy² は U1 で 0.05〜0.5 (V/Δτ が支配し、熱伝導の LHS の形はほとんど効かない)、ノズル 1 層目は point cfl 4 で約 0.006、方向別 cfl 4 で約 23。
  U1 の遅れは D ≪ 1 で対角の大きさが少し変わった影響で、D ≫ 1 の効きは試していない。「point で利点がない」は U1 からは言えない (point のノズルは D ≈ 0.006 で効く余地が無い、が正しい)。
- **V1 (run_0212、キー 1、directional cfl 4)**: 200 step まで有限 (キー 0 は 65 step で NaN)、最初の 40 step は残差が下がるが 60 step から育ち、200 step で rms_ro 8e-3 (開始の 1000 倍) → **不合格**。
  帳簿: 1 層目の Δρ/ρ (step 8) 1.09e-2 → 4.2e-3、1 層目のエネルギー行の熱伝導・粘性の段は ±40〜150 で育たない (キー 0 は ±1000 まで)。60 step 以降は**対流の段**が育つ (5 層目 −228 → −795、1 層目 −334)。
  → 熱伝導が駆動する 1 層目のモードは抑えた (第 2 仮説を支持)、残る成長は対流の結合 (仮説、未確定)。
- **V1b (run_0213、キー 2 = 拘束の行だけ)**: 72 step で NaN → 第 1 仮説の単独説を外す。
- **ノズルでの追加 (事前登録、ユーザ指示「まずノズルで試したら?」)**: キー 1 + directional を cfl 2 (`run_0214_ns_coldmesh_tw300_linedir_tj1_cfl2`) と cfl 1 (`run_0215_ns_coldmesh_tw300_linedir_tj1_cfl1`) で 2000 step・500 ごと (run_0183 の res_100000 から、同じ帳簿)。
  合格 (安定): 2000 step まで有限、全残差列の最後の 500 step の中央値が開始時の 3 倍以下で、最後の 500 step で上がり続けない。合格したものは 20000 step に延ばして V2 と同じ量 (|欠損|/ṁ_in、壁時計) で V3 と比べる。
  **V3** (`run_0216_ns_coldmesh_tw300_cfl4_tj1`、point cfl 4 + キー 1、20000 step・5000 ごと) は §6 のとおり。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-09 | [2026-10-09-time_integration-implicit-thermal-jacobian-plan.md](../../notes/reviews/2026-10-09-time_integration-implicit-thermal-jacobian-plan.md) | GO-with-changes, C0/M4/m1 | 全件採用: M1 係数を残差と同じ面の k_face と c_p/γ に (§4.1)、M2 node 限定 (`isNode && has_nbr` と設定の検査)、M3 U1 (純伝導)・V0 の範囲の拡大・線形解の失敗 0 件を §6 に、M4 欠損の大きさ/ṁ_in・局所ノルム・同じ水準までの壁時計で判定 (§6 V2)、m5 原因の書き方を仮説に・V1c で熱伝導と粘性仕事を分ける |
| plan (上位の諮問) | 2026-10-09 | diagnostician (Fable、本節に要約) | 式は残差と符号・幾何で整合、Major 2 (k の定義 = codex M1 と同じ、Jacobian 単体の検査が無い = U1)・Minor 3 | 全件採用: ∂T/∂ρ は厳密と書き直し (§4.1)、V1b (拘束の行、ビット 2) を事前登録、V2 の主比較を V3 に、V0 に run_0203 の設定と ST 1 を追加、レジスタ・スピルを記録、`nodeIsothermalEnergyBC` を併用拒否。係数は diagnostician の自節点の値でなく codex の面の値 (残差の式そのもの) を採る |

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
- `2026-10-09` — codex plan 段のレビューと diagnostician の諮問を反映 (§4・§6 改訂、ビットマスク化と V1b・U1・V1c の追加)。status in_progress。
