forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 計画立案時レビュー (stage = plan)

対象の plan は下に全文を貼る (`plans/active/time_integration-implicit-thermal-jacobian.md`)。これから実装に入る前の段階なので、次を順に評価せよ。

1. **目的とスコープ** (§1, §2): 解こうとしている課題は正しく同定されているか。既に解決済み/別 plan と重複していないか
   (`plans/README.md` と `plans/accepted/` を確認)。
2. **設計方針** (§4): 数学的・数値的に健全か。forge の既存構造 (`methods/architecture/overview.md`、該当 `methods/<area>/`)
   と整合するか。node/cell 両離散化、float32、陰解法 (block-DPLUR)、周期・軸対称などの既知の落とし穴に抵触しないか。
   代替案と比べて費用対効果は妥当か。
3. **実装ステップと残作業表** (§5, §5.1): 順序・粒度は妥当か。抜けている前提 (メッシュ品質、IC、段階起動) はないか。
4. **検証計画** (§6): 判定基準は定量的か。検証ケースの選択は `procedures/verification/README.md` と整合するか。
   「収束」「一致」を何で判定するかが書かれているか。
5. **見落としているリスク**: 我々が気づいていない構造的問題があれば挙げよ。

最後に「この計画で実装に進んでよいか」を **GO / GO-with-changes / NO-GO** の 1 語で判定し、
GO-with-changes なら実装前に直すべき点を優先順で列挙すること。

## 重点

§6.2 (2026-10-10 に追加した事前登録: point 仕上げでの面エンタルピーの精度の A/B) だけを点検する。§4・§6.0 は既に実装・検証済みで対象外。点検してほしいこと: (1) 主判定の共通の評価器を d (面エンタルピーを double で評価) にした理由が妥当か (諮問の記録と引き継ぎメモ notes/sessions/2026-10-10-handoff-face-enthalpy-audit.md は「A の設定の 1 step だけの run」を示唆していた)。(2) 判定の閾値 (10 %、ノイズ = max(再評価の差, 3 本の軌道の幅) の 3 倍、範囲が分かれる条件、床の判定の傾き −0.01 桁/500 step) と 3 本ずつの反復で足りるか。(3) ゲート (step 0 の表示と同じ状態の評価の一致 1e-6、設定の一致、状態の sha256) の抜け。(4) fh_floor.sh と fh_floor_judge.py が §6.2 の文言どおりか (台本は AWS の ~/forge-wallfit/case/45.isobutane_m6_d155 で動く。cold_pair.py・nozzle 系のモジュールはリポジトリの同じディレクトリにある)。出発の run_0354 の末尾の rms_roe は step ごとに 6.51〜6.54 を行き来する。point の 1 step は約 19 ms、2000 step は約 40 s。

## plan 全文 (`plans/active/time_integration-implicit-thermal-jacobian.md`)

```markdown
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
| 5 | point 仕上げにおける残差評価の精度の監査 (移管、2026-10-09) | [time_integration-line-viscous-jacobian](time_integration-line-viscous-jacobian.md) §6.16 で、TP の SLAU の面エンタルピーを double にすると同じ状態のエネルギーの残差の評価が ‖R_B − R_A‖/‖R_A‖ = 5.5e-3 変わった。停滞への寄与を測るなら (codex 2026-10-09): 同じ point 仕上げの保存状態・同じバイナリ (lineH 以降)・同じ設定から A = 既定、B = `FORGE_DIAG_FACE_H_DOUBLE=1` だけを変えて各 2000 step (全残差は毎 step、開始・1000・1500・2000 step の状態を保存)。末尾 500 step の全残差の水準・傾きに加え、保存した状態を共通の評価の精度で再評価する。B だけで共通の評価のエネルギーの残差が 10 % 以上下がり再評価のノイズを十分上回る → この期間の停滞への寄与を支持、自方式の表示だけが変わり共通の評価で 10 % 未満 → 主要因説を支持しない、まだ減衰中 → 床の判定は不能 (自動で延長しない)。回すかは本線の評価の後に決める (未着手)。**別のセッションへ引き継ぎ (2026-10-10)**: [`notes/sessions/2026-10-10-handoff-face-enthalpy-audit.md`](../../notes/sessions/2026-10-10-handoff-face-enthalpy-audit.md) (float 化の判断材料)。**2026-10-10 着手**: §6.2 に事前登録 (出発 run_0354 の res_40000、各腕 3 本、共通の評価は評価器 d)、codex plan 段の後に回す | F |

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
- **run_0224 の結果 (2026-10-09、run_0223 の res_15000 から 60000 step、35.9 ms/step)**: NaN なし。check_convergence (0〜59999) は NOT CONVERGED で、
  rms_ro (3.66e-6 → 1.00e-5、最大 1.47e-5) と rms_roUy (3.85e-4 → 1.29e-3) が RISING、rms_roK・rms_roOmega は下降。
  欠損 0.657 (5000) → 0.0206 kg/s (60000)。θ_r (x = 40/70/94) 0.07761 / 0.11453 / 0.12727 は run_0217 (point、通算 64 万 step) の終値の −0.16 / −0.21 / −0.24 %、
  Q_w 9.3055 MW は +0.03 %。θ_r(70) の 5000 step ごとの増分は 1.7e-3 (最大) → 2.1e-4 で、比 0.66〜0.73 で減っている (等比の外挿で 0.1149、point の終値 0.1148)。
  ただし末尾 2 万 step のドリフトは θ_r で +1.4 %、Q_w で −0.71 % で、**事前登録の水準 (0.05 %/2 万 step) には未達**。
  directional の通算は run_0223 + run_0224 = 75000 step (約 45 分)、point は 64 万 step (19〜26 ms/step で 3.4〜4.6 時間、GPU の共有の有無が混ざる)。
  → 事前登録どおり水準まで延長する (`run_0252`、同じ設定で 60000 step)。残差の RISING (場所は未確認) は延長と切り戻しの判定で併せて見る。
- **codex 諮問 (2026-10-09、[記録](../../notes/reviews/2026-10-09-line-viscous-jacobian-and-v0-diagnose.md)) の採否**: 全件採用。
  (a) V0: 事前登録の基準 (ビット一致) は不合格のまま、「回帰の帰属は判定不能」と記録する。代わりの基準は別の版として、全 `VALUE/*`・ST 0/1・固定の尺度・
  独立 run を単位とした許容差を先に登録し、独立のデータで判定する。それまで run_0191 (旧バイナリ) との正式の後退判定はしない (同じ新バイナリの中の比較だけ)。
  (b) run_0224 は「point に近づいた未収束の過渡」で、同じ不動点の証拠ではない (θ_r・Q_w は末尾 2 万 step で DRIFTING、run_0217 は判定窓の標本不足)。
  (c) 延長の判断より先に残差の場所を調べる → 下記。(d) Q_w は `cold_series.py` が CPG の定数で再構成した共通の後処理の指標で、MW の値や熱収支の根拠に使わない。
  (e) FREEZE_TURB (run_0222) の「SST と平均流の両方が寄与」は凍結の実効性の確認が未完了なので仮の結論に下げる。
- **run_0252 の結果 (run_0224 の res_60000 から同じ設定で 60000 step、通算 135000 step)**: NaN なし。check_convergence (0〜59999): NOT CONVERGED (stalled/plateau、RISING なし)。
  `cold_series.py`: 欠損 |·| ≤ 0.004 kg/s、θ_r (x = 40/70/94) 0.07787 / 0.11493 / 0.12773、ドリフト −0.010 / −0.002 / −0.003 %/2 万 step、Q_w 9.2931 (−0.005 %) → **事前登録の状態の水準に入った**。
  run_0217 (point、通算 64 万 step) の終値との差は θ_r(70) +0.14 %、Q_w −0.10 %。
  **残差の場所 (`cold_resloc.py`、res_ro の絶対値の和 Σ|res_ro|×2π)**: point の run_0217 は 0.60 kg/s で、ほぼ全部が x_w −5〜0 の壁の節点 (層 0) で最大は列 1640 付近に固定。
  directional の run_0224 は 1.1 → 3.9 → 2.95 kg/s、run_0252 は 2.0〜3.7 kg/s (point の 4〜6 倍) で、最大の場所は x_w −0.3〜+0.8 (列 1600〜2260) を動き、
  x_w 0〜5・5〜20・20 以降の区間でも point の約 10 倍。符号付きの和 (欠損) は打ち消し合いで 0 に近い。→ 積分量は point と 0.1 % 台で近いが、局所の残差は directional のほうが大きく揺れている (観測、原因は未特定)。
- **事前登録 (2026-10-09、run の前)**:
  (i) 切り戻し `run_0262_ns_coldmesh_tw300_cutback_point` = run_0252 の res_60000 から、同じ新バイナリ (35e498b1…) でキー 0・上限なし・ライン 0 の point cfl 4 (run_0217 と同じ設定)、40000 step・5000 ごと。
  合格: θ_r (x = 40/70/94) と Q_w (後処理の指標) の開始からの変化が、run_0252 の末尾 2 万 step の振幅 (最大 − 最小) 以内、または `check_quasisteady` (末尾 2 万 step、5 点、drift 0.0005) で STEADY かつ
  開始の値との差が 0.05 % 以内。欠損 |·| ≤ 0.1 kg/s。Σ|res_ro| の推移 (point の 0.60 kg/s に下がるか) は記録する。
  (ii) point 側の標本の追加 `run_0263_ns_coldmesh_tw300_cfl4_ext4` = run_0217 の res_200000 から同じ設定で 40000 step・5000 ごと。末尾 2 万 step (5 点) で `check_quasisteady` を出し、(i) と同じ量を並べる。
  (iii) 判断: (i) が合格なら「directional の最終場は point の作用素でも動かない」、不合格で point の終値側へ動けば「別の不動点か、point 側の遅い緩和」で、(ii) と併せて読む。
- **切り戻しの結果 (2026-10-09、新バイナリ 35e498b1、40000 step・5000 ごと)**:
  (i) `run_0262_ns_coldmesh_tw300_cutback_point`: θ_r (x = 40/70/94) 0.077869 / 0.114934 / 0.127731 → 0.077795 / 0.114840 / 0.127643 (−0.096 / −0.081 / −0.069 %)、Q_w +0.0093 % (2026-10-09 に系列から再計算して訂正。当初 +0.001 % と転記、codex 指摘)、
  欠損 |·| ≤ 0.004 kg/s。θ_r(70) の 5000 step ごとの増分は −2.5e-5 → 0.0 (最後の 5000 step)。`check_quasisteady` (20000〜40000、5 点、drift 0.0005): θ_r(40)・Q_w STEADY、
  θ_r(70)・θ_r(94) TRANSIENT-UNSETTLED (末尾でまだ動いている)。Σ|res_ro| は 2.4 (run_0252) → 0.48〜0.53 kg/s (point の run_0263 は 0.60)。
  → **事前登録の基準 (開始からの変化 0.05 % 以内) は不合格** (θ_r が 0.07〜0.10 % 動いた)。
  (ii) `run_0263_ns_coldmesh_tw300_cfl4_ext4`: θ_r 0.077738 / 0.114768 / 0.127568 → 0.077767 / 0.114804 / 0.127605 (+0.04 / +0.03 / +0.03 %)、Q_w 9.3025 → 9.2987 MW、
  `check_quasisteady` (20000〜40000、5 点): ALL STEADY (ただし単調に +0.014 %/2 万 step)。
  両者の終値の差: θ_r +0.036 / +0.031 / +0.030 %、Q_w −0.05 % (切り戻し − point)。両側から互いに近づいている (観測)。解釈 (同じ不動点か) は上位の諮問の後に書く (plan 未確定)。
  **codex 諮問 (2026-10-09、[記録](../../notes/reviews/2026-10-09-line-viscous-jacobian-divergence-diagnose.md)) の採否 (全件採用)**: 切り戻しは「同じ不動点」の証明になっていない
  (切り戻し側の θ_r(70・94) が TRANSIENT-UNSETTLED、変化は登録の許容を超える)。Σ|res_ro| は 2.417 → 0.493 kg/s に下がるが RMS(res_ro/V) は 5.84e4 → 4.57e5 に増えるので、
  局所の残差全般の改善とは言えない。run_0263 の外挿の漸近値は θ_r 0.0777773 / 0.114823 / 0.127632 (推定であって収束の証明ではない)。
  → directional (値 0、キー 5、上限 50) は「point で仕上げる前提の加速の候補」に留める。採用には、仕上げ込みで同じ終了条件 (全残差・局所の残差のノルムと分布・目的量の共通の条件) までの壁時計の短縮を示す。
- **codex の指摘の採否 (2026-10-09 diagnose)**: 見立ては局所モデルで確認 (キー 1 で Δρ/ρ 1 % あたり ΔT 8.0 K、キー 5 で 0 K) だが、局所行列の固有値はすべて正で発散の証明ではない → 実 run の A/B で確かめる (採用)。
  粘性仕事 (τ·u) も新しい行 4 に入っていない (第 2 仮説、採用して記録)。U1 は高速 TP の交差項を検査できない (採用、Jacobian の単体照合は未消化として残す)。
  「熱伝導が駆動するモードを抑えた」は解釈 (観測は「粘性の段の成長が弱まった」) と書き分ける (採用)。「TP の壁際で E/c_v ≈ T」は一般に成り立たない (採用、e(T_w)/c_v(T_w) で書く)。

### 6.2 事前登録: point 仕上げでの面エンタルピーの精度の A/B (§5.1 #5、2026-10-10、run の前)

**目的**: point 仕上げの残差の停滞に、TP の SLAU の面エンタルピーの float の評価 (`convectiveFlux_slau_d.inc.cuh` の `thermo_h_mix_f`) が寄与しているかを、
表示の残差の定義が変わっただけの効果と分けて測る。float 化 (plan architecture-float-state-double-geometry) で「どこを double に残すか」の判断材料にする。
設計は codex 諮問 ([記録](../../notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md) の「判別 A/B」) のとおりで、本節は条件と閾値を具体化する。

- **出発の状態 S**: `case/45.isobutane_m6_d155/run_0354_m9_L5cut/res_40000.h5`。速度 plan §6.14 の L5 の point 仕上げ (粘性入りのライン 11.5 万 step の後に point 4 万 step) の最終の場で、
  `check_convergence` は全列 STALLED (plateau、低下 0.0〜0.5 桁)。その末尾の `rms_roe` は step ごとに 6.51〜6.54 を行き来している (1 step で ±0.4 % 程度)。
- **設定**: run_0354 と同じ (point・cfl 4・`implicitThermalJacobian` 0・ライン 0・緩和 0.7・リミッタの基準値は run_0183 で固定・`extraFields: [res_ro]`)。
  `cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext <run> --cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --extra res_ro --field-from <状態>` で作り、
  `solverConfig.yaml` が run_0354 と `nStepOuter`・`outStepInterval` の他は一致することを判定のゲートにする。
- **バイナリ**: 全 run で `lineM_fp64` (`~/forge-linespeed-fp64`、sha256 05ad8bdf…、FP64)。台本は sha256 を確かめてから起動する。
  run_0354 自身は lineL (e10a195d…) で回したが、S は状態として使うだけで、比較はすべて同じバイナリの中で行う。
- **腕**: A = 既定、B = `FORGE_DIAG_FACE_H_DOUBLE=1` だけを変える (起動ログの `[DIAG] FORGE_DIAG_FACE_H_DOUBLE` の行が B にだけあることを確かめる)。
  軌道の再現のばらつきを測るため、各腕 3 本を A1・B1・A2・B2・A3・B3 の順に逐次に回す。各 2000 step、全残差は毎 step (`residual_history.csv` の `outer_begin` の行)、場は 500 step ごと。
- **共通の評価**: 保存した状態から同じ設定で 1 step だけの run を回し、step 0 の `outer_begin` の行 (更新前の残差 R(Q)) を読む。評価器は f (切替なし) と d (切替あり) の 2 つ。
  状態は S、A1・B1 の 1000・1500・2000 step、A2・A3・B2・B3 の 2000 step の 11 個 × 評価器 2 つ = 22 本。
  再評価のノイズとして、S・A1 の 2000・B1 の 2000 を各評価器でもう 1 回ずつ評価する (6 本)。
  状態は `_fh_states/<名前>/` に元の res へのシンボリックリンクだけを置いて `--field-from` に渡す (prep は restart_field でビット一致を確かめて移す)。
  判定では、各評価 run の `COLD_PAIR.json` の `parent_res_sha256` が状態の sha256 と一致することを確かめる。
- **主判定の量**: 評価器 d の `rms_roe` (面エンタルピーを double で評価した残差)。
  評価器 f の床には float の評価そのものの揺れが入るので、B の反復が残差の小さい状態に進んでも f では見えない可能性がある。そのため d を共通の物差しにする。f の結果も同じ表に出す (副)。
- **判定** (`fh_floor_judge.py`、run の前に commit する。E_X(状態) = 評価器 X の `rms_roe`):
  - ノイズの尺度: σ_eval = 同じ状態の 2 回の評価の相対差の最大、w_A・w_B = 2000 step の状態の 3 本の (最大 − 最小)/平均。noise = max(σ_eval, w_A, w_B)。r = 平均 E_d(B, 2000) / 平均 E_d(A, 2000)。
  - **支持** (この期間の停滞への寄与を支持): B の 3 本の最大 ≤ 0.9 × A の 3 本の最小、かつ 1 − r ≥ 3 × noise。
  - **支持しない** (この期間の主要因説を支持しない): r > 0.9。このとき B の自方式の表示 (末尾 500 step の中央値の 3 本の平均) が A の表示の 0.9 倍以下なら、「表示の定義の違いだけ」と書く。
  - **判別不能**: 上のどちらにも当たらない (平均は 10 % 以上違うが範囲が分かれない、またはノイズに埋もれる)。
  - **床**: B の 3 本の自方式の表示 (= 評価器 d) の、末尾 500 step (1500〜1999) の log10(`rms_roe`) の最小二乗の傾き × 500 が、3 本とも負で平均が −0.01 桁以下なら「減衰中 = 床の判定は不能」。
    そうでなければ「この期間で頭打ち」。どちらでも**自動で延長しない**。A1・B1 の 1000・1500・2000 step の E_d は記録だけにする。
  - **ゲート** (外れたら判定しない): 全 run の起動・最後の step・切替の表示・設定の一致・状態の sha256。2000 step の場の ρ・P・T が有限で正。
    A1 の step 0 の表示と E_f(S)、B1 の step 0 の表示と E_d(S) の相対差がそれぞれ 1e-6 以下 (同じ状態・同じ評価器なので、外れたら評価の手順が表示と違う残差を読んでいる)。
- **記録だけにするもの**: 同じ状態 S での E_d(S)/E_f(S) (評価の精度だけで表示がどれだけ変わるか、全列)、他の列 (ρ・運動量・k・ω・化学種) の比、
  各腕の `check_convergence --segment` の VERDICT (区間 0〜1999)。θ_r・Q_w などの目的量は判定しない (2000 step では遅いモードが動かない)。
- **run** (case/45 の 05xx): 軌道 `run_0500_fh_a1`・`run_0501_fh_b1`・`run_0502_fh_a2`・`run_0503_fh_b2`・`run_0504_fh_a3`・`run_0505_fh_b3`、
  評価 `run_0510`〜`run_0531_fhe_<状態>_{f,d}`・`run_0532`〜`run_0537_fhe_<状態>_{fr,dr}`。台本 `fh_floor.sh` (判定はしない)。
  評価 run は残差の CSV だけを使い、場は台本の中で消す。軌道 run の途中の場と `nozzle.h5` は判定の後に消し、res_2000 は A1・B1 だけ残す。
- **やらないこと**: 自動の延長、面エンタルピーの double を既定にすること (速度・回帰への影響を測っていない)、この A/B だけで float 化の範囲を決めること。

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
```

## 参考: `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md`

```
# codex 諮問 (diagnose): line-viscous-faceh

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-faceh.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-faceh.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `00cc0831` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.6 min, rc=0
- **結論**: **§6.15を限定付きの支持として記録し、値2・3の探索を終了して、既定の面エンタルピー精度を維持した本線の「方向別dt＋上限＋point仕上げ」の総壁時計評価へ戻る。**
- **extra**: `case/45.isobutane_m6_d155/_jprobe/fh_judge.json`, `case/45.isobutane_m6_d155/_jprobe/s0p7hb_compare.txt`, `plans/active/time_integration-line-implicit-speed.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0、Major 3）

| 論点 | 採否 | 根拠と対案 |
|---|---|---|
| ①「精度依存を支持、H-c の説明としても支持」 | **採用** | 保存された集計では、両側とも差分再現・再評価ノイズのゲートを通過。J_t の変化は **0.96576**、近似作用素との相対差は **0.97476 → 0.38926**、比は **0.39934**。§6.15 の分岐を満たす。[判定値](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_jprobe/fh_judge.json:5)。対象は **S0・p7・ライン2183・自由なエネルギー行・ε = 1e-6** に限定する。 |
| 「h の凍結で ṁ∂h が抜けた。近似作用素の欠陥ではなかった」 | **Major・要再検証** | 切替は温度入力だけでなく、組成・係数・演算・戻り値の精度も変える。[float係数生成](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/thermo_d.cuh:133)、[入力変換](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/thermo_d.cuh:1079)。面ごとの Δh と流束微分の内訳は未測定で、**凍結と相殺欠落という機構の確定には足りない**。Bでも相対差 **38.9%** が残る。「約36倍の大きさの隔たりには面エンタルピー評価精度が強く寄与した」と記録する。 |
| ②「収束した S0 で残差が 0.55% 変化し、残差床の可能性」 | **Major・前提を訂正して課題登録を採用** | S0 の出所 `case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/` は、保存された判定で **NOT CONVERGED (stalled/plateau)**。区間は `main`、連結300000行で、`rms_roOmega` は RISING。[判定記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/gates_aws.json:32)。また **0.005456 は ‖R_B−R_A‖/‖R_A‖** であり、残差ノルムが0.55%増減したという意味ではない。[集計式](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jprobe.py:253)。現時点では「固定状態の残差評価が変わる」まで。 |
| B の有限振幅誤差 0.232、実補正の予測差 0.406 | **Major・却下** | 両側の `pp` に旧 float 経路の `run_0317_jp_s0p7_pp` を指定している。[指定箇所](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jprobe_fh.sh:42)。Bでは **R_A(Q+0.7p)−R_B(Q)** を比較してしまう。[計算箇所](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jprobe.py:175)。Bの `finite_vs_linear` と `op*_secant_rel` は無効として掲載を外す。**§6.15 の正式判定はこれらを参照しないため維持できる**。今回は追加評価せず区切る。 |
| ③ 値2・3を保留し、本線へ戻る | **採用** | [§6.15](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:381)の停止条件どおり。面エンタルピーのdouble化による発散回避・長期収束・速度改善は測っていない。既定値と本線の残差評価精度を維持する。 |
| ④ plan を閉じ、値2・3を診断として残す | **条件付き採用** | コード残置は[既存の処置](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:241)と整合する。ただし[§5.1](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:94)の **U0・FP32未実施**を完了扱いにしない。resultレビューで各項目の処置を確定した後、現役の「本線不採用・診断として残置」という判断文書として **accepted** に移すことを推奨する。 |

結論: **§6.15を限定付きの支持として記録し、値2・3の探索を終了して、既定の面エンタルピー精度を維持した本線の「方向別dt＋上限＋point仕上げ」の総壁時計評価へ戻る。**

第 1 仮説: **面エンタルピーのfloat評価に含まれる量子化が、今回のエネルギー行の大きな方向微分と近似作用素との差に強く寄与した。** 確度: **高。ただし詳細機構は中。**

根拠: [切替コード](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:399)は面エンタルピー評価を変更し、固定した方向・作用素との比較で ‖J_t p‖ が **148.862 → 5.27835（約1/28.2）**、相対差が **0.97476 → 0.38926** になった。差分幅の再現性は、量子化された評価経路の局所応答を測れていることを示しても、連続な熱力学モデルの微分を保証しない。

反証条件: 元データの監査で、入力状態・方向・実効設定が共通でない、または保存集計が再現せず、有効な差分で J_t の変化が10%未満だった場合。**「ṁ∂h が消えた」という狭い説明は、面ごとの応答を測っていないため未確定のまま残す。**

第 2 仮説: **残る38.9%の差には、近似作用素と再構成を含む実残差とのモデル差や、他の混合精度経路が寄与する。** 確度: 中。今回スカラー対角を引いても **0.38926 → 0.38874** なので、この対象ではそのスカラーだけを主要因とする説明は支持されない。

第 3 仮説: **面エンタルピー精度がpoint仕上げの残差停滞にも寄与する。** 確度: 低。固定状態の残差差だけでは、反復の到達限界との因果を示せない。

判別 A/B: **今回追加では回さず、別課題を再開するときの一組として登録する。**

- 登録先は `plans/active/time_integration-implicit-thermal-jacobian.md` **§5.1** の「point仕上げにおける残差評価精度の監査」。速度planから参照し、今回の粘性Jacobian planには移管先を残す。
- 同じpoint仕上げの保存状態・同じバイナリ・同じ設定から、A＝現行float、B＝`FORGE_DIAG_FACE_H_DOUBLE=1` のみ変更。**各2000 step、全残差を毎step**、開始・1000・1500・2000 stepの状態を保存する。
- 末尾500 stepの全残差の水準・傾きに加え、保存状態を**共通の残差評価精度でも再評価**する。表示する残差の定義が違うだけの改善を分離する。
- **Bだけで共通評価のエネルギー残差が10%以上低下し、再評価ノイズを十分上回るなら**、この期間の停滞への寄与を支持する。**自方式の表示残差だけが変わり、共通評価では10%未満なら**、この期間の主要因説を支持しない。まだ減衰中なら残差床の判定は不能とし、自動延長しない。収束・定常性は別途各判定ツールで扱う。

やらない方がよいこと:

- 今回の支持を、§6.13の正式判定「S0・S1とも(d)判別不能」の上書きに使う。
- 面エンタルピーを直ちにdoubleへ既定化する、値2・3を本線へ戻す、熱伝導Kの列削除を続ける。
- 「datum込みだからhが大きい」を実測なしに原因へ採用する。係数にはsensible datumが反映されており、対象面でのhと微分項の大きさは未測定。
- 不合格・未実施項目を消してplanを成功完了にする。

呼び出し側の前提への異議: **S0は収束済みではない。精度依存の支持は、近似作用素の健全性や発散原因の確定ではない。Bの有限振幅指標には評価経路の混在がある。**

planを閉じる際は、§6.15の結果、上記の留保、§6.13の(d)、各試験の不合格・判定不能・未実施の処置、診断コードの残置範囲、再開条件を残す。U0は既定経路の回帰確認として扱い、FP32など延期する検証は適用範囲の制限と移管先を明記する。§5.1に残る「実装経路は正しい」「1 ulp」等の古い要約も、本文の訂正に合わせて更新してからresultレビューへ渡す。

不足情報: ローカルには `fh_judge.json` と比較TXTがある一方、元のNPZ、ライン行列、対象run、実効設定、今回の比較JSONはない。**集計の閾値判定とコードは確認したが、元配列からの独立再計算・AWSバイナリの同一性確認はできていない。** `forge`は起動せず、ファイルも変更していない。**plan未反映。**
```

## 参考: `case/45.isobutane_m6_d155/fh_floor.sh`

```
#!/bin/bash
# plan time_integration-implicit-thermal-jacobian §6.2 (2026-10-10): point 仕上げでの面エンタルピーの精度の A/B (§5.1 #5)。
# 出発 S = run_0354_m9_L5cut/res_40000.h5 (point 仕上げの最終の場)。設定は run_0354 と同じ (point・cfl 4・キー 0・リミッタの基準値は run_0183)。
# A = 既定 (面エンタルピーは float)、B = FORGE_DIAG_FACE_H_DOUBLE=1。
#   軌道 6 本 (A1 B1 A2 B2 A3 B3 の順、各 2000 step・500 ごと)。
#   共通の評価 28 本: 保存した状態から同じ設定で 1 step。step 0 の outer_begin の行が更新前の残差 R(Q)。
#     評価器 f = 切替なし、d = 切替あり。11 状態 × 2 + 再評価のノイズ 3 状態 × 2。
# バイナリは lineM_fp64 (sha256 05ad8bdf…) に固定。判定は fh_floor_judge.py (この台本は判定しない)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge COLD_ALT_BINARY=lineM_fp64
SRC=run_0183_ns_coldmesh_tw300_ext
S0=run_0354_m9_L5cut
ST=_fh_states
LOG=$PWD/fh_floor.log
B="--cfl 4 --limiter-ref-from $SRC --extra res_ro"

echo "== 開始 $(date -Is)" >> $LOG
if pgrep -x forge > /dev/null; then echo "他の forge が走っている — 止める ($(pgrep -x forge | tr '\n' ' '))" >> $LOG; exit 1; fi
[ "$(sha256sum $FORGE_BIN | cut -c1-16)" = 05ad8bdf6100874f ] || { echo "FORGE_BIN の sha256 が lineM_fp64 でない — 止める" >> $LOG; exit 1; }

# 切替の表示を確かめる (B と評価器 d は表示がある、A と評価器 f は無い)
switch_ok() {  # switch_ok <run> <fh>
  if [ "$2" = 1 ]; then grep -q "FORGE_DIAG_FACE_H_DOUBLE" $1/forge_run.log
  else ! grep -q "FORGE_DIAG_FACE_H_DOUBLE" $1/forge_run.log; fi
}

state() {  # state <名前> <run> <res>: 状態だけを置いたディレクトリ (prep の --field-from は最後の res を取る)
  mkdir -p $ST/$1
  [ -f $2/$3 ] || { echo "状態 $2/$3 が無い" >> $LOG; return 1; }
  ln -sfn ../../$2/$3 $ST/$1/$3
  sha256sum $2/$3 | cut -d' ' -f1 > $ST/$1/SHA256
}

traj() {  # traj <run> <fh>
  local r=$1 fh=$2
  python3 cold_cfl.py prep $SRC $r --steps 2000 --out 500 $B --field-from $S0 >> $LOG 2>&1 || { echo "prep $r 失敗" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  ( [ "$fh" = 1 ] && export FORGE_DIAG_FACE_H_DOUBLE=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r fh=$fh rc=$rc $(date -Is) last=$(awk -F, '$3=="outer_begin"{s=$1} END{print s}' $r/residual_history.csv)" >> $LOG
  switch_ok $r $fh || { echo "$r: 切替の表示が腕と合わない — 止める" >> $LOG; return 1; }
  [ -f $r/res_2000.h5 ] || { echo "$r: res_2000.h5 が無い — 止める" >> $LOG; return 1; }
}

ev() {  # ev <run> <状態の名前> <fh>
  local r=$1 sd=$2 fh=$3
  python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 $B --field-from $ST/$sd >> $LOG 2>&1 || { echo "prep $r 失敗" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  ( [ "$fh" = 1 ] && export FORGE_DIAG_FACE_H_DOUBLE=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r state=$sd fh=$fh rc=$rc $(date -Is) step0=$(awk -F, '$1==0&&$3=="outer_begin"{print $8}' $r/residual_history.csv)" >> $LOG
  switch_ok $r $fh || { echo "$r: 切替の表示が評価器と合わない — 止める" >> $LOG; return 1; }
  grep -q "^0,-1,outer_begin" $r/residual_history.csv || { echo "$r: step 0 の行が無い — 止める" >> $LOG; return 1; }
  rm -f $r/*.h5    # 評価の run は残差の CSV だけを使う (場は消す)
}

state s0 $S0 res_40000.h5 || exit 1
i=0
for arm in a1:0 b1:1 a2:0 b2:1 a3:0 b3:1; do
  n=${arm%%:*}; fh=${arm##*:}
  traj run_$(printf %04d $((500 + i)))_fh_$n $fh || exit 1
  i=$((i + 1))
done
state a1k1000 run_0500_fh_a1 res_1000.h5 && state a1k1500 run_0500_fh_a1 res_1500.h5 && state a1k2000 run_0500_fh_a1 res_2000.h5 &&
state b1k1000 run_0501_fh_b1 res_1000.h5 && state b1k1500 run_0501_fh_b1 res_1500.h5 && state b1k2000 run_0501_fh_b1 res_2000.h5 &&
state a2k2000 run_0502_fh_a2 res_2000.h5 && state b2k2000 run_0503_fh_b2 res_2000.h5 &&
state a3k2000 run_0504_fh_a3 res_2000.h5 && state b3k2000 run_0505_fh_b3 res_2000.h5 || exit 1

k=510
for sd in s0 a1k1000 a1k1500 a1k2000 b1k1000 b1k1500 b1k2000 a2k2000 a3k2000 b2k2000 b3k2000; do
  ev run_$(printf %04d $k)_fhe_${sd}_f $sd 0 || exit 1; k=$((k + 1))
  ev run_$(printf %04d $k)_fhe_${sd}_d $sd 1 || exit 1; k=$((k + 1))
done
for sd in s0 a1k2000 b1k2000; do   # 再評価のノイズ (同じ状態・同じ評価器をもう 1 回)
  ev run_$(printf %04d $k)_fhe_${sd}_fr $sd 0 || exit 1; k=$((k + 1))
  ev run_$(printf %04d $k)_fhe_${sd}_dr $sd 1 || exit 1; k=$((k + 1))
done
echo "== 終了 $(date -Is)" >> $LOG
touch fh_floor.done
```

## 参考: `case/45.isobutane_m6_d155/fh_floor_judge.py`

```
#!/usr/bin/env python3
"""point 仕上げでの面エンタルピーの精度の A/B の判定 (plan time_integration-implicit-thermal-jacobian §6.2、2026-10-10 事前登録)。

fh_floor.sh の run (軌道 run_0500〜0505、評価 run_0510〜0537) を読み、§6.2 の規則どおりに判定して
_band_ab/cold_pair/fh_floor_judge.json に書く。規則を変えるときは plan §6.2 を先に改訂する (結果を見てから変えない)。

  主判定の量: 評価器 d (面エンタルピーを double で評価) の rms_roe を、A・B の 2000 step の状態で比べる。
  E_X(状態) = 評価 run の step 0 の outer_begin の行 (更新前の残差)。X = f (既定) / d (切替あり)。
"""
import csv
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
COLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
Q = "rms_roe"
S0, S0_RES = "run_0354_m9_L5cut", "res_40000.h5"
TRAJ = {"a1": ("run_0500_fh_a1", 0), "b1": ("run_0501_fh_b1", 1), "a2": ("run_0502_fh_a2", 0),
        "b2": ("run_0503_fh_b2", 1), "a3": ("run_0504_fh_a3", 0), "b3": ("run_0505_fh_b3", 1)}
STATES = ["s0", "a1k1000", "a1k1500", "a1k2000", "b1k1000", "b1k1500", "b1k2000", "a2k2000", "a3k2000", "b2k2000", "b3k2000"]
REPEAT = ["s0", "a1k2000", "b1k2000"]
# 事前登録の閾値 (§6.2)
DROP = 0.10          # 10 % の低下
NOISE_MULT = 3.0     # 低下はノイズの 3 倍以上
DECAY_SLOPE = -0.01  # 末尾 500 step の log10(rms_roe) の傾き × 500 [桁]
GATE_SELF = 1e-6     # 軌道の step 0 の表示と同じ状態・同じ評価器の評価の相対差の上限
TAIL = (1500, 1999)
EXPECT_DIFF = {"output", "space/limiterARef", "space/limiterPRef", "space/limiterRoRef",
               "time/deltaT/cfl", "time/deltaT/cfl_pseudo", "time/last/nStepOuter"}


def eval_runs():
    """(状態, 評価器, 回) → run 名。fh_floor.sh と同じ並び。"""
    out, k = {}, 510
    for sd in STATES:
        for x in ("f", "d"):
            out[(sd, x, 1)] = f"run_{k:04d}_fhe_{sd}_{x}"
            k += 1
    for sd in REPEAT:
        for x in ("f", "d"):
            out[(sd, x, 2)] = f"run_{k:04d}_fhe_{sd}_{x}r"
            k += 1
    return out


def outer_rows(run):
    rows = {}
    with open(HERE / run / "residual_history.csv") as f:
        for r in csv.DictReader(f):
            if r["phase"] == "outer_begin":
                rows[int(r["step"])] = {c: float(r[c]) for c in COLS}
    return rows


def flat(d, p=""):
    out = {}
    for k, v in d.items():
        kk = f"{p}/{k}" if p else str(k)
        if isinstance(v, dict):
            out.update(flat(v, kk))
        else:
            out[kk] = v
    return out


def cfg_diff(run, ref):
    a = flat(yaml.safe_load((HERE / run / "solverConfig.yaml").read_text()))
    b = flat(yaml.safe_load((HERE / ref / "solverConfig.yaml").read_text()))
    skip = ("time/last/nStepOuter", "time/outStepInterval")
    return sorted(k for k in set(a) | set(b) if k not in skip and a.get(k) != b.get(k))


def switch_on(run):
    return "FORGE_DIAG_FACE_H_DOUBLE" in (HERE / run / "forge_run.log").read_text(errors="replace")


def rel(a, b):
    return abs(a - b) / (0.5 * (abs(a) + abs(b)))


def slope500(rows, q=Q):
    s = np.array([k for k in range(TAIL[0], TAIL[1] + 1)], dtype=float)
    y = np.log10(np.array([rows[int(k)][q] for k in s]))
    return float(np.polyfit(s, y, 1)[0] * 500.0)


def main():
    gates, rec = [], {"plan": "time_integration-implicit-thermal-jacobian §6.2", "quantity": Q}
    ev = eval_runs()
    # --- ゲート: 起動・切替・設定・状態の同一性 ---
    traj_rows = {}
    for name, (run, fh) in TRAJ.items():
        rows = outer_rows(run)
        traj_rows[name] = rows
        last = max(rows)
        gates.append((f"{run}: 最後の step 1999", last == 1999, last))
        gates.append((f"{run}: 切替の表示 = {fh}", switch_on(run) == bool(fh), switch_on(run)))
        d = cfg_diff(run, S0)
        gates.append((f"{run}: 設定が run_0354 と同じ (出力の間隔・step 数を除く)", d == [], d))
        cp = json.loads((HERE / run / "COLD_PAIR.json").read_text())
        gates.append((f"{run}: config_diff", set(cp["config_diff"]) <= EXPECT_DIFF | {"time/outStepInterval"}, cp["config_diff"]))
        with h5py.File(HERE / run / "res_2000.h5", "r") as h:
            bad = {v: int(np.sum(~np.isfinite(h["VALUE"][v][...]) | (h["VALUE"][v][...] <= 0))) for v in ("ro", "P", "T")}
        gates.append((f"{run}: res_2000 の ρ・P・T が有限・正", sum(bad.values()) == 0, bad))
    E = {}
    for (sd, x, n), run in ev.items():
        rows = outer_rows(run)
        gates.append((f"{run}: step 0 の行", 0 in rows, sorted(rows)[:3]))
        gates.append((f"{run}: 切替の表示 = {x}", switch_on(run) == (x == "d"), switch_on(run)))
        d = cfg_diff(run, S0)
        gates.append((f"{run}: 設定が run_0354 と同じ", d == [], d))
        cp = json.loads((HERE / run / "COLD_PAIR.json").read_text())
        want = (HERE / "_fh_states" / sd / "SHA256").read_text().strip()
        gates.append((f"{run}: 状態 {sd} の sha256", cp["parent_res_sha256"] == want, cp["parent_res_sha256"][:16]))
        E[(sd, x, n)] = rows.get(0)
    ok = all(g[1] for g in gates)
    # 自己一致: 軌道の step 0 の表示 (A1 = f、B1 = d) と、S の評価 (同じ状態・同じ評価器)
    self_f = rel(traj_rows["a1"][0][Q], E[("s0", "f", 1)][Q])
    self_d = rel(traj_rows["b1"][0][Q], E[("s0", "d", 1)][Q])
    gates.append(("A1 の step 0 の表示と E_f(S) の相対差 ≤ 1e-6", self_f <= GATE_SELF, self_f))
    gates.append(("B1 の step 0 の表示と E_d(S) の相対差 ≤ 1e-6", self_d <= GATE_SELF, self_d))
    ok = ok and self_f <= GATE_SELF and self_d <= GATE_SELF
    rec["gates"] = [{"check": g[0], "ok": bool(g[1]), "value": g[2]} for g in gates]
    rec["gates_ok"] = bool(ok)
    # --- 量 ---
    rec["E"] = {f"{sd}|{x}|{n}": E[(sd, x, n)] for (sd, x, n) in E}
    sig = {x: max(rel(E[(sd, x, 1)][Q], E[(sd, x, 2)][Q]) for sd in REPEAT) for x in ("f", "d")}
    rec["sigma_eval"] = sig
    rec["eval_effect_at_S"] = {c: E[("s0", "d", 1)][c] / E[("s0", "f", 1)][c] for c in COLS}
    per = {}
    for x in ("f", "d"):
        A = [E[(f"a{i}k2000", x, 1)][Q] for i in (1, 2, 3)]
        Bv = [E[(f"b{i}k2000", x, 1)][Q] for i in (1, 2, 3)]
        mA, mB = float(np.mean(A)), float(np.mean(Bv))
        wA, wB = (max(A) - min(A)) / mA, (max(Bv) - min(Bv)) / mB
        noise = max(sig[x], wA, wB)
        r = mB / mA
        if max(Bv) <= (1 - DROP) * min(A) and (1 - r) >= NOISE_MULT * noise:
            v = "SUPPORT"
        elif r > 1 - DROP:
            v = "NOT_SUPPORT"
        else:
            v = "INDETERMINATE"
        per[x] = {"A": A, "B": Bv, "mean_A": mA, "mean_B": mB, "ratio": r, "w_A": wA, "w_B": wB, "noise": noise, "verdict": v,
                  "all_cols_ratio": {c: float(np.mean([E[(f"b{i}k2000", x, 1)][c] for i in (1, 2, 3)]) /
                                         np.mean([E[(f"a{i}k2000", x, 1)][c] for i in (1, 2, 3)])) for c in COLS}}
    rec["by_evaluator"] = per
    # 自方式の表示 (各腕の CSV の末尾 500 step)
    disp = {}
    for name, rows in traj_rows.items():
        disp[name] = {"tail_median": {c: float(np.median([rows[k][c] for k in range(TAIL[0], TAIL[1] + 1)])) for c in COLS},
                      "tail_slope500": {c: slope500(rows, c) for c in COLS if all(rows[k][c] > 0 for k in range(TAIL[0], TAIL[1] + 1))},
                      "step0": rows[0]}
    rec["display"] = disp
    dA = np.mean([disp[f"a{i}"]["tail_median"][Q] for i in (1, 2, 3)])
    dB = np.mean([disp[f"b{i}"]["tail_median"][Q] for i in (1, 2, 3)])
    rec["display_ratio_B_over_A"] = float(dB / dA)
    # 減衰中か (B の表示 = 評価器 d の末尾 500 step の傾き)
    sl = [disp[f"b{i}"]["tail_slope500"][Q] for i in (1, 2, 3)]
    decaying = float(np.mean(sl)) <= DECAY_SLOPE and all(s < 0 for s in sl)
    rec["B_tail_slope500"] = sl
    rec["floor"] = "DECAYING (床の判定は不能、自動で延長しない)" if decaying else "FLAT (この期間で頭打ち)"
    rec["B1_Ed_series"] = {k: E[(f"b1k{k}", "d", 1)][Q] for k in (1000, 1500, 2000)}
    rec["A1_Ed_series"] = {k: E[(f"a1k{k}", "d", 1)][Q] for k in (1000, 1500, 2000)}
    # 主判定 (評価器 d)
    if not ok:
        rec["VERDICT"] = "INVALID (ゲート不合格: 判定しない)"
    else:
        v = per["d"]["verdict"]
        txt = {"SUPPORT": "この期間の停滞への寄与を支持",
               "NOT_SUPPORT": "この期間の主要因説を支持しない",
               "INDETERMINATE": "判別不能 (平均は 10 % 以上違うが範囲が分かれない、またはノイズに埋もれる)"}[v]
        if v == "NOT_SUPPORT" and rec["display_ratio_B_over_A"] <= 1 - DROP:
            txt += " (自方式の表示だけが 10 % 以上下がった = 表示の定義の違い)"
        rec["VERDICT"] = f"{v}: {txt}; 床: {rec['floor']}"
    out = HERE / "_band_ab" / "cold_pair" / "fh_floor_judge.json"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=float))
    print(json.dumps({k: rec[k] for k in ("gates_ok", "VERDICT", "sigma_eval", "display_ratio_B_over_A", "B_tail_slope500")},
                     ensure_ascii=False, indent=1, default=float))
    for x in ("f", "d"):
        p = per[x]
        print(f"評価器 {x}: A {['%.5g' % a for a in p['A']]} B {['%.5g' % b for b in p['B']]} 比 {p['ratio']:.4f} ノイズ {p['noise']:.2e} → {p['verdict']}")
    for g in gates:
        if not g[1]:
            print("ゲート不合格:", g[0], g[2])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
```

## 参考: `case/45.isobutane_m6_d155/cold_cfl.py`

```
"""冷却壁の腕の擬似 CFL の試行と、質量の収支の監査の準備 (plan tooling-nozzle-isothermal-wall-chain §5.1 #27)。
cold_pair.py の prep_ext (延長) と同じ手順 (restart_field でビット一致・FP64 の型のまま、段なし) で、変えてよい設定を
nStepOuter・cfl・cfl_pseudo・outStepInterval と output.extraFields に限る (implicitRelax は変えない — codex 2026-10-08。例外は --relax を明示した §6.13 の run)。
走行中の run が cold_pair.py を使っているので、cold_pair.py は書き換えずにここで包む。

usage (AWS の case dir、別バイナリは COLD_ALT_BINARY=<キー> を前に付ける):
  python3 cold_cfl.py prep <src_run> <run> --steps N --cfl C [--out 5000] [--extra res_ro,volume] [--limiter-ref-from <run>]
    --limiter-ref-from: リミッタの基準値 (limiterRoRef・limiterPRef・limiterARef) を指定した run の forge_run.log の値に固定する。
      既定 (自動) では開始場から決まるので、restart した run は親と別の作用素になる (forge の警告; run_0182 → run_0183 で a_ref が 1.6 % 違った)。
    --line dir|only: 壁法線のライン陰解法 (lineImplicit 1)。dir は方向別の擬似 dt (lineDtDirectional 1) も足す (§5.1 #27 の試行)。
    --isp 0|1: time.deltaT.implicitSolvePrecision (陰解法の行列の組立て・解法の精度。既定 0 = float は FP64 のビルドでも float)。
  python3 cold_pair.py run <run>       (投入は既存の run_one; FORGE_DUMP_MASSFLUX は投入側の環境変数で渡す)
  python3 cold_cfl.py run <run>        (同じ run_one を、COLD_ALT_BINARY の登録を効かせて呼ぶ)
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_pair as CP  # noqa: E402
NS = CP.NS

# 別バイナリの登録 (plan time_integration-implicit-thermal-jacobian §6 の検証用)。環境変数 COLD_ALT_BINARY にキーを書くと、
# cold_pair.binary_record の照合 (FORGE_SHA・FP64_TREE) をこの登録値に差し替える。変換器は従来の FP64 のもの (CONV_SHA) のまま。
ALT_BINARIES = {
    # キー: (forge の sha256, ソースの作業ツリー)
    # implicitThermalJacobian の検証用: commit 1ad0b9bb + typedef double (座標の stod は HEAD に入っている)、2026-10-09 AWS でビルド
    "thermjac_fp64": ("9ffc4d1efcec6ca4799f4026ef418931cf99bd27ad61f29d264cbfc4639018b7", "~/forge-thermjac-fp64"),
    # 06d1b149 (ビット 4 を追加) + typedef double、2026-10-09 に同じ作業ツリーで再ビルド (上のバイナリは上書きされて残っていない)
    "thermjac5_fp64": ("985aca0f2ebba912851042ec9abc8b7deab9707b20f19f757641dc9b3521215e", "~/forge-thermjac-fp64"),
    # b4771052 (lineDtDirectionalCap を追加) + typedef double、2026-10-09 に同じ作業ツリーで再ビルド
    "thermjac_cap_fp64": ("35e498b14b5f3cdaa09b754f7bbaa6455f54631fd2edf1ad4929964e8828d4bc", "~/forge-thermjac-fp64"),
    # 4d394a71 (ライン上の節点の対角を storeLU の sweep 以外で組まない、plan time_integration-line-implicit-speed 案 A) + typedef double、2026-10-09 AWS で新しい作業ツリーにビルド
    "linespeed_fp64": ("d8b06ebcb91cfc3151cfcedf441b3e79f20f4c3702aa89e287801cb1c4f1d10b", "~/forge-linespeed-fp64"),
    # 5ab83056 (lineViscCoupling 2 = 薄層の粘性・熱伝導の Jacobian、案 A を含む) + typedef double、2026-10-09 AWS でビルド
    # (上のバイナリは 2026-10-09 に下の lineB で上書きされて残っていない)
    "linevisc_fp64": ("6631a87eb1279b2fa4eff22080378ac12db4937890ac4470cf0d9185daf55653", "~/forge-linevisc-fp64"),
    # 64ed2cd6 (ライン内で並列にした Thomas・新旧の比較モード・ライン行列の書き出し、値 2 と案 A を含む) + typedef double、~/forge-linevisc-fp64 を上書きしてビルド
    # (上の lineB は下の lineB2 で上書きされて残っていない)
    "lineB_fp64": ("7a7e9eb9150104babc5b1c9877e15eb826eaa46fde8f2fe96dce74b65ad61f6d", "~/forge-linevisc-fp64"),
    # 874ae90a (並列 Thomas の代入を行の分担に、行のポインタを 1 回だけ選ぶ、書き出しのフラグの修正) + typedef double、~/forge-linevisc-fp64 を上書き
    # (上の lineB2 は下の lineB3 で上書きされて残っていない)
    "lineB2_fp64": ("f3765961a601cff3f63e5ed86e0f0cf84d62138f504ba4a03dc317330ea45fad", "~/forge-linevisc-fp64"),
    # 3764852d (lu5 の行の入れ替えを静的な添字にしてレジスタに置く、既定は 1 ライン 1 スレッド、並列は FORGE_LINE_PAR=1) + typedef double
    # (上の lineB3 は下の lineD で上書きされて残っていない)
    "lineB3_fp64": ("9a094160caa9b3a68cbad20cd7b035a529f8f31a3de34fcf82590039040f6d3f", "~/forge-linevisc-fp64"),
    # d5538001 (lu5 を元に戻す。既定 = 案 A + 値 2 の経路 + 並列 Thomas は FORGE_LINE_PAR=1 + 比較・書き出しのデバッグ) + typedef double
    # (上の lineD は下の lineE で上書きされて残っていない)
    "lineD_fp64": ("9ef80d5f148f73b1ee330a141a6f5fde31737a996ea0b5611e749e79b9442e8d", "~/forge-linevisc-fp64"),
    # 6d49738e (診断用の lineViscCoupling 3 を追加) + typedef double、~/forge-linevisc-fp64 を上書き
    # (上の lineE は下の lineF で上書きされて残っていない)
    "lineE_fp64": ("8e989dc688b7ccef70c444aecf06fb9ec406fa0a7f4d2ad0dc8158cef7851178", "~/forge-linevisc-fp64"),
    # 4fdc2c0e (診断のマスク FORGE_LVC_TERMS) + typedef double、~/forge-linevisc-fp64 を上書き
    # (上の lineF は下の lineG で上書きされて残っていない)
    "lineF_fp64": ("89ea94385435c16babde906465b58a54748d0e4b2d50b5edc661b46903ae19b0", "~/forge-linevisc-fp64"),
    # 846727be (マスクのビット 8 = 熱伝導の K の密度の列を外す) + typedef double、~/forge-linevisc-fp64 を上書き
    "lineG_fp64": ("c0b649052fdcf6092c7c73317563d77b570df60aca31e1bf709cf47586000c8c", "~/forge-linevisc-fp64"),
    # 5197e00e (診断の FORGE_DIAG_FACE_H_DOUBLE = TP の面エンタルピーを double で、plan time_integration-line-viscous-jacobian §6.15) + typedef double、
    # ~/forge-linevisc-fp64 を上書き (上の lineG は残っていない。切替なしの残差が lineG と同じことは run_0323_jph_a_q0 で確かめる)
    "lineH_fp64": ("561813b3564473420824ccca302b795bf9c344f69776f6381a43cb23e0ee6456", "~/forge-linevisc-fp64"),
    # 270f1d75 (opt-in の逆行列の保存 FORGE_LINE_INV と比較の経路の非有限・全ラインの η、plan time_integration-line-implicit-speed §6.2) + typedef double、~/forge-linevisc-fp64 を上書き (lineH は残っていない)
    "lineI_fp64": ("b467c3d7a1b32b595ac0ed8cc5db33f434329ef39865d727abfd11ef835b7661", "~/forge-linevisc-fp64"),
    # d1d0de7c (opt-in の float の Thomas FORGE_LINE_F32=1/2 と比較の経路の非有限の検査、plan time_integration-line-implicit-speed §6.4) + typedef double、~/forge-linevisc-fp64 を上書き (lineI は残っていない)
    "lineJ_fp64": ("3d045221b8677ca108af365a4012010c3e8f0e65e5d9e3ef53827d9cef963c0f", "~/forge-linevisc-fp64"),
    # 00da938b (opt-in の Thomas の配列の並べ替え FORGE_LINE_LAYOUT=1 とビット列の比較・因子の比較、plan time_integration-line-implicit-speed §6.7) + typedef double、~/forge-linevisc-fp64 を上書き (lineJ は残っていない)
    "lineK_fp64": ("c62eaf5a3910f9b1573ebbb1e740727d3c81b0a356a9352ae167ed67e93eec27", "~/forge-linevisc-fp64"),
    # f9be0c4f (opt-in の並べ替え + 連鎖の短縮 FORGE_LINE_LAYOUT=2、plan time_integration-line-implicit-speed §6.10) + typedef double、~/forge-linevisc-fp64 を上書き (lineK は残っていない)
    "lineL_fp64": ("e10a195dbbec9bd4b67650348985d68d6a5354e1f6884ab5e73f11974f085e3d", "~/forge-linevisc-fp64"),
    # LAYOUT2 を既定にしたコミット (plan time_integration-line-implicit-speed §6.21) + typedef double、~/forge-linespeed-fp64 (f9be0c4f に timeIntegration_d.cu を入れてビルド)
    "lineM_fp64": ("05ad8bdf6100874f39ae39a5d9ea08613223e87e1e4fdc5270baaeaad5206b11", "~/forge-linespeed-fp64"),
}


def _use_alt_binary():
    import os
    k = os.environ.get("COLD_ALT_BINARY", "")
    if not k:
        return None
    if k not in ALT_BINARIES:
        raise SystemExit(f"COLD_ALT_BINARY={k} は登録されていない ({sorted(ALT_BINARIES)}) — 止める")
    sha, tree = ALT_BINARIES[k]
    CP.FORGE_SHA = sha
    CP.FP64_TREE = Path(tree).expanduser()
    return k


def limiter_refs(run: Path) -> dict:
    """forge_run.log の最後の limiterRoRef・limiterPRef・limiterARef (forge が貼るよう促す行) を読む。"""
    import re
    ref = {}
    for line in (run / "forge_run.log").read_text(errors="replace").splitlines():
        m = re.match(r"^\[limiter\]\s+(limiterRoRef|limiterPRef|limiterARef):\s*([0-9.eE+-]+)\s*$", line)
        if m:
            ref[m.group(1)] = m.group(2)
    if set(ref) != {"limiterRoRef", "limiterPRef", "limiterARef"}:
        raise SystemExit(f"{run.name}/forge_run.log にリミッタの基準値の 3 行がそろっていない: {ref}")
    return ref


def prep(src: Path, run: Path, steps: int, cfl: float, out_int: int, extra: list[str], ref_from: Path | None = None,
         line: str = "", isp: int | None = None, inner: int | None = None, conv: int | None = None, itj: int | None = None,
         cap: float | None = None, lvc: int | None = None, field_from: Path | None = None, relax: float | None = None) -> dict:
    NS.check_dry_env(False)
    binrec = CP.binary_record()
    if not NS.RUN_RE.match(run.name) or run.exists():
        raise SystemExit(f"{run} の名前が不正か既にある — 止める")
    srec = NS.jload(src / CP.RECORD)
    rs = NS.res_files(src)
    if not rs:
        raise SystemExit(f"{src} に res が無い")
    src_h5 = rs[-1]
    if field_from is not None:              # 設定は src、場は別の run の最終の res (同じ格子。切り戻し試験用)
        fr = NS.res_files(field_from)
        if not fr:
            raise SystemExit(f"{field_from} に res が無い")
        src_h5 = fr[-1]
    ys = NS.yaml_strict()
    ptext = (src / "solverConfig.yaml").read_text()
    ctext = ys.replace_scalars(ptext, {NS.NSTEP: str(int(steps)), NS.CFL: repr(float(cfl)), NS.CFLP: repr(float(cfl)),
                                       NS.OUTINT: str(int(out_int))})
    pcfg = ys.load(ptext)
    if extra:
        if "output" in pcfg:
            raise SystemExit("親の設定に output がある — extraFields の足し方を決めていないので止める")
        ctext = ctext.rstrip("\n") + "\noutput: {level: 1, extraFields: [" + ", ".join(extra) + "]}\n"
    refs = limiter_refs(ref_from) if ref_from is not None else {}
    if refs:
        if any(k in pcfg.get("space", {}) for k in refs):
            raise SystemExit("親の設定に既にリミッタの基準値がある — 止める")
        if ctext.count("space: {") != 1:
            raise SystemExit("space が 1 行のフロー形式でない — 基準値の足し方を決めていないので止める")
        ctext = ctext.replace("space: {", "space: {" + ", ".join(f"{k}: {v}" for k, v in refs.items()) + ", ")
    # line: "" = なし、"dir" = lineImplicit + lineDtDirectional、"only" = lineImplicit だけ (方向別の擬似 dt なし)
    line_keys = {"dir": {"lineImplicit": 1, "lineDtDirectional": 1}, "only": {"lineImplicit": 1},
                 "dirvisc": {"lineImplicit": 1, "lineDtDirectional": 1, "lineViscCoupling": 1}, "": {}}[line]
    if line_keys:                           # 壁法線のライン陰解法 + 方向別の擬似 dt (procedures/solver-settings.md「lineImplicit」)
        dt = pcfg["time"]["deltaT"]
        if int(dt.get("blockDPLUR", 0)) != 1 or int(pcfg["time"].get("timeIntegration", 0)) != 11 or int(dt.get("lowMachPrecond", 0)) >= 2:
            raise SystemExit("lineImplicit は timeIntegration 11 + blockDPLUR 1 + lowMachPrecond < 2 専用 — 止める")
        if any(k in dt for k in line_keys) or ctext.count("deltaT: {") != 1:
            raise SystemExit("deltaT に既にライン陰解法のキーがあるか、deltaT が 1 つのフロー形式でない — 止める")
        ctext = ctext.replace("deltaT: {", "deltaT: {" + ", ".join(f"{k}: {v}" for k, v in line_keys.items()) + ", ")
    if isp is not None:                     # 陰解法の行列の組立て・解法の精度 (0 = float、1 = double; 既定 0 は FP64 のビルドでも float)
        if "implicitSolvePrecision" in pcfg["time"]["deltaT"] or ctext.count("deltaT: {") != 1:
            raise SystemExit("deltaT に既に implicitSolvePrecision があるか、deltaT が 1 つのフロー形式でない — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{implicitSolvePrecision: {int(isp)}, ")
    if itj is not None:                     # implicitThermalJacobian (plan time_integration-implicit-thermal-jacobian、別バイナリ)
        if "implicitThermalJacobian" in pcfg["time"]["deltaT"] or ctext.count("deltaT: {") != 1:
            raise SystemExit("deltaT に既に implicitThermalJacobian があるか、deltaT が 1 つのフロー形式でない — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{implicitThermalJacobian: {int(itj)}, ")
    if cap is not None:                     # 方向別 dt の伸びの上限 R (lineDtDirectionalCap、line dir が要る)
        if not line_keys.get("lineDtDirectional"):
            raise SystemExit("--cap は --line dir と一緒に使う — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{lineDtDirectionalCap: {float(cap)!r}, ")
    if lvc is not None:                     # lineViscCoupling (2 = 薄層の粘性 Jacobian、plan time_integration-line-viscous-jacobian、line が要る)
        if not line_keys.get("lineImplicit") or "lineViscCoupling" in pcfg["time"]["deltaT"]:
            raise SystemExit("--lvc は --line と一緒に使う (親に lineViscCoupling があっても止める) — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{lineViscCoupling: {int(lvc)}, ")
    if relax is not None:                   # implicitRelax を明示して変える (plan time_integration-line-implicit-speed §5.1 #19・§6.13 の事前登録の run だけ)
        if ctext.count("deltaT: {") != 1 or "implicitRelax" not in pcfg["time"]["deltaT"]:
            raise SystemExit("deltaT が 1 つのフロー形式でないか親に implicitRelax が無い — 止める")
        import re as _re
        ctext, nsub = _re.subn(r"implicitRelax:\s*[0-9.eE+-]+", f"implicitRelax: {float(relax)!r}", ctext)
        if nsub != 1: raise SystemExit("implicitRelax の書き換えが 1 か所でない — 止める")
    one = {}
    if inner is not None:
        one[("time", "nStepInner")] = str(int(inner))
    if conv is not None:
        one[NS.CONVP] = str(int(conv))
    if one:                                  # 単因子の切り分け用 (nStepInner・convMethod)
        ctext = ys.replace_scalars(ctext, one)
    allowed = ({NS.NSTEP, NS.CFL, NS.CFLP, NS.OUTINT} | set(one) | ({("output",)} if extra else set())
               | {("space", k) for k in refs} | {("time", "deltaT", k) for k in line_keys}
               | ({("time", "deltaT", "implicitSolvePrecision")} if isp is not None else set())
               | ({("time", "deltaT", "implicitThermalJacobian")} if itj is not None else set())
               | ({("time", "deltaT", "lineDtDirectionalCap")} if cap is not None else set())
               | ({("time", "deltaT", "lineViscCoupling")} if lvc is not None else set())
               | ({("time", "deltaT", "implicitRelax")} if relax is not None else set()))
    diff = set(NS.MK.diff_paths(pcfg, ys.load(ctext)))
    if not diff <= allowed:
        raise SystemExit(f"許していない設定の差がある: {sorted(diff - allowed)} — 止める")
    if relax is None and ys.load(ctext)["time"]["deltaT"]["implicitRelax"] != pcfg["time"]["deltaT"]["implicitRelax"]:
        raise SystemExit("implicitRelax が変わった — 止める")
    if relax is not None and float(ys.load(ctext)["time"]["deltaT"]["implicitRelax"]) != float(relax):
        raise SystemExit("implicitRelax が指定どおりでない — 止める")
    run.mkdir(parents=True)
    for fn in NS.EXT_COPY + ("wall_repr.json", "bcondConfig.yaml", "species_meta.yaml"):
        if (src / fn).is_file():
            shutil.copy2(src / fn, run / fn)
    for p in sorted(src.glob("resolved_species_*.yaml")):
        shutil.copy2(p, run / p.name)
    (run / "solverConfig.yaml").write_text(ctext)
    cmd = [sys.executable, str(NS.TOOLS / "restart_field.py"), str(src_h5), str(run / "nozzle.h5"), "--dst-run", str(run), "--keep-src-dtype"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=NS.runner()._ENV)
    (run / "restart_field.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0 or "ビット一致" not in (r.stdout + r.stderr):
        print((r.stdout + r.stderr)[-3000:])
        raise SystemExit(f"restart_field がビット一致を確認していない (rc {r.returncode}) — 止める")
    info = NS.jload(run / "prepare_info.json")
    info.update(stages={"stages": "none", "ramp": None, "ramp_steps": 1000}, extends=src.name, restart_from=f"{src_h5.parent.name}/{src_h5.name}")
    NS.jdump(run / "prepare_info.json", info)
    keep = ("plan", "kind", "problem", "problem_sha256", "delta_r_csv", "euler_ref", "implicit_relax", "mesh_checks",
            "geometry_vs_production", "wall_thermal", "bcond_wall")
    rec = {**{k: srec[k] for k in keep if k in srec},
           "tool": "cold_cfl.py prep", "plan_item": "§5.1 #27", "created": NS.now(), "git_head": NS.git_head(), "binary": binrec,
           "stages": "none", "parent": src.name, "field_from": field_from.name if field_from is not None else None, "parent_res": src_h5.name, "parent_res_sha256": NS.sha256_file(src_h5),
           "ext_steps": int(steps), "cfl_main": float(cfl), "cfl_parent": srec.get("cfl_main"), "out_interval": int(out_int),
           "extra_fields": extra, "limiter_ref_from": ref_from.name if ref_from is not None else None, "limiter_refs": refs, "line_keys": line_keys, "implicit_solve_precision": isp, "n_step_inner": inner, "conv_method": conv, "implicit_thermal_jacobian": itj, "line_dt_directional_cap": cap, "line_visc_coupling": lvc, "implicit_relax_override": relax,
           "config_diff": sorted("/".join(p) for p in diff),
           "restart_field_tail": (r.stdout + r.stderr).strip().splitlines()[-1:], "nozzle_sha256_after_prep": NS.sha256_file(run / "nozzle.h5")}
    NS.jdump(run / CP.RECORD, rec)
    print(f"[cold_cfl prep] {run.name} ← {src.name}/{src_h5.name}: cfl {cfl}・{steps} step・出力 {out_int} ごと・extra {extra}; "
          f"{rec['restart_field_tail']}")
    return rec


if __name__ == "__main__":
    alt = _use_alt_binary()
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("run"); p.add_argument("run")
    p = sp.add_parser("prep"); p.add_argument("src"); p.add_argument("run")
    p.add_argument("--steps", type=int, required=True); p.add_argument("--cfl", type=float, required=True)
    p.add_argument("--out", type=int, default=5000); p.add_argument("--extra", default="")
    p.add_argument("--limiter-ref-from", default=None)
    p.add_argument("--line", choices=("dir", "only", "dirvisc"), default="",
                   help="dir = lineImplicit 1 + lineDtDirectional 1、only = lineImplicit 1 だけ、dirvisc = dir + lineViscCoupling 1")
    p.add_argument("--isp", type=int, choices=(0, 1), default=None, help="time.deltaT.implicitSolvePrecision を書く")
    p.add_argument("--inner", type=int, default=None, help="time.nStepInner を変える")
    p.add_argument("--conv", type=int, default=None, help="space.convMethod を変える (0 = 1 次)")
    p.add_argument("--itj", type=int, default=None, help="time.deltaT.implicitThermalJacobian を書く (別バイナリ COLD_ALT_BINARY が要る)")
    p.add_argument("--cap", type=float, default=None, help="time.deltaT.lineDtDirectionalCap を書く (別バイナリが要る)")
    p.add_argument("--field-from", default=None, help="場だけをこの run の最終の res から取る (設定は src、切り戻し試験用)")
    p.add_argument("--lvc", type=int, default=None, help="time.deltaT.lineViscCoupling を書く (2 = 薄層の粘性 Jacobian、別バイナリが要る)")
    p.add_argument("--relax", type=float, default=None, help="time.deltaT.implicitRelax を書き換える (§6.13 の事前登録の run だけ。既定は親のまま)")
    a = ap.parse_args()
    if a.cmd == "run":                      # cold_pair.run_one を (別バイナリの登録を効かせて) 呼ぶ
        sys.exit(CP.run_one(HERE / a.run))
    prep(HERE / a.src, HERE / a.run, a.steps, a.cfl, a.out, [s for s in a.extra.split(",") if s],
         HERE / a.limiter_ref_from if a.limiter_ref_from else None, a.line, a.isp, a.inner, a.conv, a.itj, a.cap, a.lvc, HERE / a.field_from if a.field_from else None, a.relax)
```

## 出力形式

1. 冒頭に **判定 (GO / GO-with-changes / NO-GO)** と 3 行以内の要約。
2. 指摘一覧 (Critical → Major → Minor の順、番号付き。各項目に根拠と対案)。
3. 推奨 (1 つに絞る)。
4. 末尾に `指摘数: Critical N / Major N / Minor N` の 1 行。
