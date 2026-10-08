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

### 4.3 行 4 に従来のスカラーも残す (ビット 4、キー 5 = 1 + 4; 2026-10-09 ユーザ判断「じゃあ b」)

キー 1 の不安定 (§6.0) の仮説: 行 0〜3 は従来のスカラー A = Σ 2ν_eff δ/dcc (スペクトル半径の近似、物理の微分ではない人工的な減衰) のまま、行 4 だけを
温度の微分 B·ρc_v·ΔT に替えたので、行ごとに歩幅の縮め方の基準が食い違う。温度を変えない補正 (例: R = R₀(1, u, v, w, E)) を解くと
ΔT = (E/c_v)(Δρ/ρ)·A/(V/Δτ + B) の偽の温度の補正が出る (V/Δτ で規格化)。x = 70 の μ_t 最大の層 (A/(V/Δτ) 0.71、B/(V/Δτ) 0.54、E/c_v ≈ 1740 K) では Δρ/ρ 1 % ごとに約 8 K。
壁際 (u ≈ 0) では E/c_v が e(T_w)/c_v(T_w) 程度で小さい。キー 1 の V3 の NaN の場所 (μ_t/μ 最大 474 の層) と一致する。

ビット 4 では行 4 にも A を残し、温度の項 B·∂(ρc_vT)/∂Q をその上に足す: 行 4 = A·e₅ᵀ + B·[−(e−½|u|²), −u, −v, −w, 1]。
全行に同じ人工的な減衰 A がかかるので、温度を変えない補正は従来どおり向きを保って縮み (ΔT = 0)、温度が変わる補正にだけ熱伝導の見込みが加わる。
期待: point で従来 (キー 0) と同等に安定、directional で 1 層目の熱伝導のモードは引き続き抑える。対流の段の成長 (V1 で 60 step 以降) が残るかは別問題。

### 4.4 方向別 dt の伸びの上限 R (`time.deltaT.lineDtDirectionalCap`、2026-10-09)

SST 凍結の判別 (run_0222、§6.0) で SST と平均流の両方が寄与 (中間) と出たので、diagnostician の事前の分岐どおり上限 R を入れる。
`setCFL_cell_d` (setDT_d.cu) で、line 上の節点について line 面を除いた max (方向別) に加えて全面の max (point と同じ) も取り、
cfl = max(cfl_方向別, cfl_全面/R) とする (Δτ ∝ 1/cfl なので Δτ = min(Δτ_方向別, R·Δτ_point))。平均流と SST は同じ `dt_local` を使うので両方に効き、
壁の節点 (境界半割面で律速) と 1 層目の Δτ の比も R 倍までに収まる。R = 0 (既定) で従来どおりビット同一。lineImplicit 1・lineDtDirectional 1 が要る (設定で検査)。
1 層目の熱拡散数 D (方向別で約 23) を 0.5 以下にするには R ≲ 80 (diagnostician)。

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
| 3 | 検証 U1・V0〜V4 (§6) | U1・V1・V1b・ノズルの cfl 2/1・V3 は不合格 (§6.0)。V0 は判定不能 (forge 自体が再実行でビット一致しない、§6.0)。代わりの基準 (同じバイナリの再実行の差の分布と比べる) の採否を codex 諮問で決める | F |
| 4 | result 段のレビューと採否 (不採用で閉じる案)、コードを残すか戻すかの判断 (ユーザ) | | F |

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
- **(b) キー 5 の事前登録 (2026-10-09、§4.3。codex diagnose [記録](../../notes/reviews/2026-10-09-thermjac-key5-diagnose.md) を反映して改訂)**:
  新しいバイナリ (06d1b149 + typedef double、sha256 985aca0f…、`cold_cfl.py` の `thermjac5_fp64`) で、**同じバイナリ・同じ開始場 (run_0183 の res_100000)・point cfl 4・nStepInner 5・relax 0.7・基準値固定で、キー 1 と 5 だけを変える A/B** を 2000 step:
  `run_0218_ns_coldmesh_tw300_cfl4_tj1b` (キー 1、新バイナリでの再現) と `run_0219_ns_coldmesh_tw300_cfl4_tj5` (キー 5)。場は 200 step ごと (全保存量の残差を含む)、
  帳簿は下流の壊れた領域 (列 4295〜4320 × 壁から 76〜92 層目) を毎 step (最初の 400 回)、after_eos_bc の位相で δρ/ρ・δT/T・δP/P・速度の補正を初期の固定尺度で測る。
  判定: A (キー 1) が run_0216 と同じく育ち、B (キー 5) が 2000 step まで育たなければ「行 4 の A を戻せば救える」を支持 (第 1 仮説を支持、ただし粘性仕事の第 2 仮説は残る)。
  B でも同じ場所で育てば「A の欠落が主因」の単独説を棄却。A が再現しなければ判別不能 (バイナリ・設定・restart の同一性に戻る)。
  「育たない」の定義: NaN・非物理値 (ρ・P・T ≤ 0) なし、全残差列の末尾 500 step の中央値が開始の 3 倍以下かつ末尾 500 step の傾きが有意に正でない、
  帳簿の領域の振動成分の振幅が末尾 3 窓 (各 100 step) で増えない。check_convergence の区間と VERDICT も保存する。合格は「2000 step で成長を検出しなかった」に限る (収束・高速化とは呼ばない)。
  run_0191 (旧バイナリのキー 0) との後退判定 (V3 と同じ基準) は V0 合格の後。Vb-dir (directional + キー 5) はこの A/B の後に決める。
- **(b) の結果 (2026-10-09)**: `run_0218` (キー 1、新バイナリ) は 711 step で NaN (run_0216 の再現)。**`run_0219` (キー 5) は 2000 step で成長を検出しなかった**
  (NaN・非物理値なし、末尾 500 step の中央値/開始 = ρ 0.62・ρu 0.89・ρv 0.37・ρE 0.64・k 1.23・ω 1.38、log10 の傾き ≤ +0.011 (ω だけ +0.084)、下流の領域の 200 step ごとの ρ の変化 1.7〜2.9e-5 で横ばい;
  check_convergence (全区間、stage_manifest が無いので --segment 不可) = NOT CONVERGED (stalled/plateau)、RISING・DIVERGED なし)。→ 「行 4 の A を戻せば救える」を支持 (第 2 仮説は残る)。
  **Vb-dir `run_0220` (directional cfl 4 + キー 5、200 step)**: 有限、最初の 100 step は残差が下がり (rms_ro 8.0e-6 → 4.2e-6)、110 step から上がる (199 step で 1.8e-5)。事前の基準 (残差が悪化しない) では不合格。
  **ただし帳簿と場の比較で解釈が変わる**: 1〜5 層目の Δρ/ρ・ΔT/T は振動でなく単調 (列 40・1 層目で 20 step +0.22 % → 200 step +2.2 %、Δp/p ≈ 0)、
  縮流部の近壁 (1〜60 層) の ΔT/T は point cfl 4 の 44 万 step の変化 (run_0183 → run_0208) と符号が 100 % 一致・相関 0.62、
  欠損 (Σ res_ro×2π) は 200 step で 2.026 → 1.423 (point cfl 4 は 4 万 step で 1.398) — **遅い過渡を約 200 倍/step で早送りしている** (仮説: 残差の増加は過渡の速さによる)。
  なお point cfl 4 の 44 万 step で縮流部の近壁 (列 40、10〜20 層) は ρ +42〜45 %・T −30 % 変わった (run_0183 は近壁で定常から遠かった)。
  **延長の事前登録**: `run_0221_ns_coldmesh_tw300_linedir_tj5_long` = run_0220 と同じ設定で 5000 step・500 ごと (`extraFields: [res_ro, volume]`)。
  合格 (安定に過渡を進める): NaN・非物理値なし、欠損 (Σ res_ro×2π、大きさ) が 500 step ごとに減り続ける (増加 2 回連続で不合格)、全残差列が開始の 100 倍を超えない、
  縮流部の近壁の点 (列 40、1・5・20 層) の ρ・T が point cfl 4 の 44 万 step の値の方向へ進み、振動 (500 step 間の符号反転の繰り返し) が育たない。
  合格したら、point cfl 4 の延長 (run_0217) と同じ欠損の水準・θ_r で比べ、方向別 + キー 5 の場を本線に切り替えるかをユーザと決める。
  **結果 (`run_0221`、69 ms/step)**: NaN・非物理値なし。残差は 500 step で rms_ro 4.7e-3 (開始の 580 倍)・rms_roOmega 1.9e5 まで振れ、3000 step で開始以下に戻り、5000 step で rms_ro 3.9e-6・rms_roUy 4.9e-4・rms_roOmega 427 (いずれも開始の約半分)。
  欠損は 2.03 → 4.88 (500) → 2.31 (1000) → 1.38 (2000) → 1.11 (3000) → 1.00 (5000)。近壁 (列 40) は行き過ぎて戻る (20 層の ρ: +25 % (500) → +53 % (1000) → +30 % (5000))。
  → **事前の基準では不合格** (残差が開始の 100 倍を超えた)。ただし過渡の後は回復し、欠損 1.0 kg/s までの壁時計は point cfl 4 (約 9.5 万 step) の約 1/8。最後の 1000 step の減り方は −4.7 %/1000 step (point は約 −0.6 %)。
  解釈 (振れが不動点の同じ過渡か不安定の兆候か、振れの抑え方、本線に使う条件) は diagnostician に諮る (2026-10-09、ユーザ指示「あきらめたくない、diagnostician にきいて」)。
  **diagnostician の回答 (2026-10-09、要約)**: 「過渡を早送りしているだけ」(呼び出し側の上の仮説) は**支持しない**: 欠損 (状態の関数) が point の軌道の最大値 2.03 を超えて 4.88 になったので point の軌道を離れている。
  100 step の減少の後に一定率 (e 倍 約 55〜60 step) で 3 桁育つのは、その状態での作用素の線形不安定の形。回復は状態依存か非線形の飽和。「200 倍/step」「壁時計 1/8」は GPU 共有中の ms/step と軌道外の状態の比較なので**速さの根拠に使わない (撤回)**。
  見落としていた危険 (第 1 仮説、確度 中): SST の点陰解法は平均流と同じ `dt_local` を使い (`update_d.cu` 379〜392)、生産は lagged (`ransSource_d.cu` 277〜284) — 方向別 Δτ で Newton 風の大股になり μ_t を通じて平均流を揺らす。キー 5 はこの経路を直していない。
  run_0210 (キー 0 の FREEZE_TURB) はキー 5 の遅い成長の除外証拠にならない。第 2 仮説: 平均流側 (壁節点の Δτ の遅れによる P_w の不整合、粘性仕事の欠落)。第 3 仮説 (前線の通過、確度 低)。
  **判別 A/B (事前登録)**: `run_0222_ns_coldmesh_tw300_linedir_tj5_freeze` = run_0221 と同じ + `FORGE_FREEZE_TURB=1` だけ、1000 step・100 ごと (`extraFields: [res_ro, volume]`)、対照 run_0221。
  A (110〜500 step の log10(rms_ro) の傾き ≤ 0.003/step、500 step で rms_ro < 8e-5、欠損 < 2.03 で単調) → 第 1 仮説 (対策: SST の更新だけ point の Δτ)。
  B (傾き ≥ 0.01/step、500 step で rms_ro ≥ 1e-4) → 第 1 仮説を外し第 2 仮説 (対策: Δτ = min(Δτ_dir, R·Δτ_point)、R = 50 から)。中間 (0.003〜0.01) → 両方が寄与 (R 上限から)。凍結の実効性は壁・入口ピン以外の roK/roOmega が res_0 と一致することで確認。
  並行して CFD 0 step で run_0221 の残差の履歴 (110〜500 step の傾き・R²・周期、どの列が先に立ち上がるか) と res_500〜2000 の res_ro の地図 (場所が固定か移動か) を出す。
  **判別の結果 (2026-10-09)**: run_0222 (FREEZE_TURB) の 110〜500 step の log10(rms_ro) の傾き +0.0052/step (run_0221 は +0.0086、R² 0.983)、500 step の rms_ro 2.4e-4 (run_0221 は 4.7e-3)、
  欠損 1.67 (100) → 1.28 (300) → 1.25 (500) → 1.22 (1000) で行き過ぎなし → **中間 (SST と平均流の両方が寄与)**。凍結の実効性: roK の res_0 からの最大相対変化は 1000 step で 2.9e-2 (壁・入口ピンの除外は未確認)。
  run_0221 の残差の履歴: 先に立ち上がったのは rms_ro・rms_roUy (3 倍を超えた step 184) で、rms_roOmega は 231、rms_roK は 455 (第 1 仮説の反証条件の側)。自己相関に周期約 7.4 step の成分。
  |res_ro|/V の最大は 500〜2000 step で縮流部の入口寄りの壁の節点 (列 35〜42、0 層) に固定、3000 step 以降はスロートの壁の節点 (列 1650) に移る。
  → 事前の分岐どおり上限 R (§4.4) に進む。**事前登録**: `run_0223_ns_coldmesh_tw300_linedir_tj5_cap50` = run_0221 と同じ (directional cfl 4・キー 5) + `lineDtDirectionalCap: 50`、run_0183 の res_100000 から 15000 step・500 ごと。
  合格 (diagnostician 案): (i) 全残差列が各出力で開始の 10 倍以内、(ii) 欠損が 2.03 を一度も超えず 500 step ごとに減る、(iii) 列 40・1/5/20 層の ρ・T が反転しない (出力間の増分の符号反転 0 回)、
  (iv) 欠損 1.0 kg/s に 15000 step 以内に達する。(iv) だけ満たさなければ「安定化したが速さが無い」で不採用。ms/step はその時の GPU の共有状況と一緒に記録し、速さの比較は専有 GPU で測り直す。
  **結果 (`run_0223`、新バイナリ b4771052 + typedef double、sha256 35e498b1…、36.45 ms/step・GPU はこの 1 本)**: 15000 step、非有限 0、非物理値 0。
  (i) 各列の最大/開始 = ρ 1.00・ρu 1.01・ρv 1.62・ρE 1.00・k 2.73・ω 3.43 → 合格。(ii) 欠損の最大 1.413 (2.03 を超えない) だが 500 → 1000 step で 1.406 → 1.413 と 1 回増えた → **不合格 (1 回)**。
  (iii) 列 40・1/5/20 層の ρ・T の増分の符号反転 0 回 → 合格。(iv) 欠損 1.0 に 7500 step で到達 → 合格。欠損は 0.98 (7500) → 0.78 (15000)。
  近壁 (列 40) の最終は ρ +11.1 %・+32.0 %・+43.2 % (1/5/20 層)、point cfl 4 の 44 万 step 後は +11.4 %・+33.3 %・+45.3 % (同じ向き・行き過ぎなし)。
  速さ (専有に近い GPU の ms/step): 欠損 1.0 まで 7500 × 36 ms ≈ 4.5 分 vs point cfl 4 の約 9.5 万 × 26 ms ≈ 42 分 (約 9 倍)。欠損 0.8 付近の減り方は −3.1 %/1000 step vs 約 −0.7 % (壁時計で約 3 倍、進むほど差が縮む)。
  **次 (diagnostician の条件)**: `run_0224_ns_coldmesh_tw300_linedir_tj5_cap50_ext` = run_0223 の res_15000 から同じ設定で 60000 step・5000 ごと (定常の水準 = 欠損 ≤ 0.1 kg/s かつ θ_r のドリフト ≤ 0.05 %/2 万 step まで)。
  水準に達したら切り戻し (同じ新バイナリでキー 0・cap 0 の point cfl 4、2〜4 万 step) で θ_r・欠損・Q_w がプラトーの振幅を超えて動かないことを見る。V0 (旧/新バイナリ・キー 0 で 20 step のビット一致) も並行して回す。
  本線に使う条件 (diagnostician): 同じ状態水準 (欠損 ≤ 入口の 0.1 %、θ_r のドリフト ≤ 0.05 %/2 万 step) に達した後の窓の時間平均 ± 振幅で point と比べ、directional の最終場から point に戻して 2〜4 万 step で動かないこと (切り戻し)、V0、専有 GPU の ms/step で同じ状態水準までの壁時計。
- **V0 の結果 (2026-10-09、`run_0225`〜`run_0232`、`v0_bitident.sh`)**: 4 組 (point / directional × ST 0 / 1) とも旧・新バイナリの res_20 が**ビット不一致** (`_band_ab/cold_pair/V0_bitident.json`)。
  forge は面の流束などを `atomicAdd` で足すので、同じバイナリの再実行でもビット一致しない可能性がある。そこで**事後に**同じバイナリの再実行を測った
  (`run_0233`〜`run_0248`、`v0_repeat.sh`・`v0_compare.py` → `_band_ab/cold_pair/V0_repeat.json`):
  同じバイナリでも 20 step で 0/3 組、1 step で 0/1 組がビット不一致 (1 step の T だけは全組一致)。場の差の相対 RMS (ρ) は
  point ST 0・20 step で旧 × 旧 3.3〜3.5e-9・新 × 新 3.0〜3.6e-9・旧 × 新 3.1〜3.6e-9、directional ST 0・20 step で 1.0〜1.2e-7・0.85〜1.0e-7・0.90〜1.2e-7、
  point・1 step で 8.6e-13・8.6e-13・2.2〜8.4e-13、directional・1 step で 1.4e-14・1.6e-14・1.6〜1.9e-14 (ρE・ρω・T も同じ傾向)。
  → **事前登録の基準 (ビット一致) は forge 自体が満たさないので判定不能**。旧 × 新の差は 4 群とも同じバイナリの再実行の差と同じ幅に入っている (観測)。
  代わりの基準 (再実行の差の分布と比べる) は結果を見た後の変更になるので、採否は codex 諮問の後に決める (plan 未決、§5.1 #3)。
- **point の水準 (run_0217、`cold_series.py`)**: run_0208 の res_200000 から 20 万 step で欠損 0.0538 → 0.0089 kg/s、θ_r (x = 40/70/94) の終値 0.07774 / 0.11477 / 0.12757、
  ドリフト +0.024 / +0.021 / +0.019 %/2 万 step、Q_w 9.3025 MW (−0.028 %/2 万 step) → 事前登録の状態の水準に入っている (残差 CSV は 113604 step で途切れており、収束判定はその区間だけ)。
  run_0183 の res_100000 (= run_0223 の開始) の θ_r(70) は 0.10077、point の 44 万 step 後 (run_0217 の開始) は 0.11412、64 万 step 後は 0.11477 (+13.9 %)。
  directional の run_0223 は 15000 step で 0.10213 (+1.3 %)、2500 step ごとの増分は 1.5〜3.4e-4 で増えている。Q_w は 13.53 → 10.55 MW (point の終値 9.30 MW)。
  近壁の ρ (列 40) は point の 44 万 step 後と同程度まで進んだが、下流の θ_r・Q_w はまだ point の水準から離れている (観測。原因は未特定)。
- **codex の指摘の採否 (2026-10-09 diagnose)**: 見立ては局所モデルで確認 (キー 1 で Δρ/ρ 1 % あたり ΔT 8.0 K、キー 5 で 0 K) だが、局所行列の固有値はすべて正で発散の証明ではない → 実 run の A/B で確かめる (採用)。
  粘性仕事 (τ·u) も新しい行 4 に入っていない (第 2 仮説、採用して記録)。U1 は高速 TP の交差項を検査できない (採用、Jacobian の単体照合は未消化として残す)。
  「熱伝導が駆動するモードを抑えた」は解釈 (観測は「粘性の段の成長が弱まった」) と書き分ける (採用)。「TP の壁際で E/c_v ≈ T」は一般に成り立たない (採用、e(T_w)/c_v(T_w) で書く)。

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
- `2026-10-09` — 実装 (1ad0b9bb)、FP64 ビルド (sha256 9ffc4d1e…)。U1・V1・V1b・ノズルの cfl 2/1・V3 はいずれも不合格 (§6.0)。キー 1 はノズルを不安定化させる → 不採用の方向。
