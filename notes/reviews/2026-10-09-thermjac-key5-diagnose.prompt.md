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

## 依頼: 診断・設計判断の諮問 (stage = diagnose)

あなたは forge の**診断・設計判断係**である。呼び出し側は実装と run を進めている別のモデル (Claude) で、
**もっともらしい真因に飛びつく前に**あなたに諮っている。仕事は手を動かすことではなく、**次の一手を 1 つに絞ること**。

### 前提
- あなたは呼び出し側の会話を見ていない。下のブリーフと、自分で読んだファイルだけが根拠になる。
  足りなければ推測で埋めずに「何が足りないか」を返す。
- ブリーフは「観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説」に分かれて渡される約束である。
  **観測事実と呼び出し側の解釈が混ざっていたら、まず分け直す**。呼び出し側の要約より、run の数値・コード・
  設定ファイルを自分で確かめた内容を優先する。
- forge を起動しない。`python3` による `residual_history.csv` / `res_*.h5` の読み取りは**統計量だけ**を出す
  (全量ダンプ・長いログ全文をコンテキストに流さない。`*.log`・`*.vtu`・`plans/README.md` は読まない)。

### 診断の作法
1. **「除外済み」というラベルを信用せず、潰した証拠を確認する** (run パス・設定差分・判定区間・VERDICT)。
   証拠が足りない・判定期間が短い・変えた設定が実際には効いていない (YAML の階層違い等) なら**候補へ戻す**。
   証拠が十分な候補は出し直さない。
2. **症状と原因を分ける**。`detectNaN` が指す変数は結果であって原因ではない (EOS 床 → 負密度 → 圧力暴走 → ω の実績)。
   後処理のアーチファクト (2 列混在の抽出、`centCoords` の置換、ソルバ `ypls` の退化) を先に疑う。
3. **このリポジトリで繰り返された真因**を照合する: 投入設定の不整合 (IC と BC、亜音速に超音速 BC)、
   押し出し 2 ノード spanwise、float32 桁落ち (双対幾何・r 重み)、stale build、cross-mesh IC の基底不一致、
   絶対値のゼロ割ガード、境界ノードの凍結、YAML キーの階層違いで黙って無視される設定。
4. 仮説は**確度順に最大 3 つ**。第 1 仮説には根拠を `ファイル:行` か run の数値で付ける。示せないものは「未確認」と明記。
5. **判別する A/B を 1 つだけ**提案する。安く短く回せて、結果がどちらに出ても仮説が 1 つ消えるもの。
   「A なら仮説 1、B なら仮説 2」を先に書く (結果を見てから解釈を作らない)。
6. 少数点の一致・短い窓の値・未収束のトランジェント同士の比較を根拠にしない。

### 設計判断 (plan §4・§6、codex 指摘の採否、result 段の解釈) を諮られたとき
- 採否は指摘ごとに「採用 / 却下 / 要再検証」と理由。根拠が示されていない指摘は自分で該当箇所を読んでから判定する。
- 検証計画は「何が出たら方針が誤りと言えるか」が定量的に書かれているかを見る。
- 既定値の変更・opt-in 機能の削除は、plan の処置欄とユーザ決定の履歴を確認してから判断する
  (「opt-in 残置」は削除対象でない)。
- result 段の解釈は、主張ごとに根拠 run・判定ツールの VERDICT・判定区間が揃っているかを確かめる
  (過渡ピークを定常値と、抽出アーチファクトを物理と誤認した実績は「予想どおり」に見える場面で起きた)。

あなたの結論は**仮説**であって確定ではない。呼び出し側はこの A/B を回して確かめ、plan への反映も呼び出し側が行う。

## ブリーフ (`notes/reviews/briefs/2026-10-09-thermjac-key5.md`)

# 諮問: 熱伝導 Jacobian のキー 1 がノズルを不安定化させた理由の見立てと、対策 (b) キー 5 の妥当性

日付 2026-10-09。諮問先 codex (diagnose)。エスカレーション条件 4 (原因を書く前) と 6 (cuda_forge の数値の振る舞いの変更)。
plan: `plans/active/time_integration-implicit-thermal-jacobian.md` (§4.1・§4.3・§6・§6.0 を全文読むこと)。作業ツリー `/home/sano/work/forge-integ-1005` (commit 06d1b149)。
コード: `solver_density_cuda/cuda_forge/timeIntegration_d.cu` の `implicit_defect_correction_block_d` の粘性の対角 (ビット 1・4 の分岐、`thermalJac`)。
ユーザ: 「じゃあ b」(行 4 に従来のスカラーも残す案)。

## 観測事実 (§6.0 にも記載)

- キー 1 (行 4 を温度の熱伝導の Jacobian に置き換え、行 0〜3 は従来のスカラー A = Σ 2ν_eff δ/dcc):
  - point cfl 4 (run_0216、case/45 の 300 K 冷却ノズル、run_0183 の res_100000 から): 最初の 200 step は従来 (run_0191) より速く残差が下がる (rms_roe 10.4 → 4.5) が、
    300 step 前後から運動量の残差が育ち 922 step で NaN。非有限は x_w 69.9〜70.5 の壁から 79〜89 層目。キー 0 (run_0191) は同じ条件で 40000 step 安定。
  - その場所 (run_0208 の res_200000、キー 0 の収束に近い場、列 4307・壁から 88 層目): μ_t/μ 474 (列内の最大)、ρ 0.031 kg/m³、T 261 K、|u| 1710 m/s、c 324 m/s、
    e = ρE/ρ − ½|u|² = −1.16e5 J/kg (TP の顕熱の基準)、½|u|² = 1.46e6 J/kg、γ 1.377、c_v 775。
    概算で A/(V/Δτ) = 0.71、B/(V/Δτ) = 0.54 (B = Σ k_f δ/dcc /(ρ c_v)、k_f は乱流の分を含む面の値)。
  - directional (cfl 4/2/1) + キー 1 は 1 層目の熱伝導の段の成長を止めたが、数百 step で対流の段が育って発散 (run_0212・0214・0215)。
- 純伝導 (U1、静止、一様格子、D = Δτκ/Δy² ≤ 0.5) ではキー 0 と 1 は同じ解に収束し、キー 1 は同じ cfl で 12〜40 % 遅い。cfl 50 は両方とも負の P・ρ で壊れる。

## 呼び出し側の見立て (仮説、確かめていない)

行 0〜3 は物理の微分でない人工的な減衰 A (スペクトル半径の近似) のまま、行 4 だけを温度の微分 B·ρc_v·ΔT に替えたので、行ごとに歩幅の縮め方の基準が食い違う。
温度を変えない残差 R = R₀(1, u, v, w, E) を対角ブロックだけで解くと (V/Δτ = 1 に規格化、対流 A⁺ は無視):
Δρ = R₀/(1 + A)、ΔρE は行 4 から ΔρE(1 + B) = E R₀ + B(e + ½|u|²)Δρ、よって ΔT = (E/c_v)(Δρ/ρ)·A/(1 + B)。
x = 70 の値で E/c_v ≈ 1740 K、A/(1 + B) ≈ 0.46 → Δρ/ρ 1 % ごとに約 8 K の偽の ΔT。壁際 (u ≈ 0) では E/c_v ≈ T で小さい。
NaN の場所 (μ_t 最大・高速) と一致する。

## 対策 (b) = キー 5 (ビット 1 + 4、実装済み 06d1b149)

行 4 = A·e₅ᵀ (従来のスカラー) + B·[−(e−½|u|²), −u, −v, −w, 1] (温度の項)。行 0〜3 は従来どおり。
全行に同じ A がかかるので、温度を変えない補正は従来どおり向きを保って縮み、温度が変わる補正にだけ熱伝導の見込みが加わる。
事前登録 (plan §6.0): Vb-point (point cfl 4 + キー 5、20000 step、早い判定 2000 step)、Vb-dir (directional cfl 4 + キー 5、200 step、帳簿)。

## 問い

1. 上の見立て (行ごとの減衰の基準の食い違いが、高速・高 μ_t の所で偽の ΔT を作る) は、観測 (point の NaN の場所・時期、directional の対流の段の成長) と整合するか。対流 A⁺ を入れると結論は変わるか。別の機構の候補は。
2. (b) は見立てに対する対策として妥当か。A と B の両方が行 4 にかかることで、温度のモードの減衰が二重になる影響 (過剰な減衰で遅くなるだけか、別の不整合を生むか)。
3. 事前登録の Vb-point・Vb-dir の判定に穴はないか。(b) が point で安定でも directional で対流の段が育つ場合、次に何を見るべきか。

## 読んでよいもの

- 上記 plan、`methods/time_integration/implementation.md` の「エネルギー行の熱伝導 Jacobian」節、`solver_density_cuda/cuda_forge/timeIntegration_d.cu`・`viscousFlux_d.cu`
- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #27 (経緯)、`case/45.isobutane_m6_d155/README.md` (run_0191・0203〜0216 の行)
- 過去の諮問: `notes/reviews/2026-10-09-linedir-divergence-2-diagnose.md`、`notes/reviews/2026-10-09-time_integration-implicit-thermal-jacobian-plan.md`

## 関連 plan 全文 (`plans/active/time_integration-implicit-thermal-jacobian.md`)

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
壁際 (u ≈ 0) では E/c_v ≈ T で小さい。キー 1 の V3 の NaN の場所 (μ_t/μ 最大 474 の層) と一致する。

ビット 4 では行 4 にも A を残し、温度の項 B·∂(ρc_vT)/∂Q をその上に足す: 行 4 = A·e₅ᵀ + B·[−(e−½|u|²), −u, −v, −w, 1]。
全行に同じ人工的な減衰 A がかかるので、温度を変えない補正は従来どおり向きを保って縮み (ΔT = 0)、温度が変わる補正にだけ熱伝導の見込みが加わる。
期待: point で従来 (キー 0) と同等に安定、directional で 1 層目の熱伝導のモードは引き続き抑える。対流の段の成長 (V1 で 60 step 以降) が残るかは別問題。

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
| 3 | 検証 U1・V0〜V4 (§6) | U1・V1・V1b・ノズルの cfl 2/1・V3 は不合格 (§6.0)。残り: V0 (既定でビット同一) — コードを残すなら必須 | O |
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
- **(b) キー 5 の事前登録 (2026-10-09、§4.3)**: 新しいバイナリ (ビット 4 を足した commit) で
  **Vb-point** `run_0218_ns_coldmesh_tw300_cfl4_tj5` = point cfl 4 + キー 5、run_0183 の res_100000 から 20000 step・5000 ごと (`extraFields: [res_ro, volume]`)。
  早い判定: 2000 step まで NaN なし、かつ全残差列が開始時の 10 倍を超えない (キー 1 は 300 step から育ち 922 で NaN)。
  合格: 20000 step で V3 と同じ条件 (|欠損|/ṁ_in が run_0191 の同じ step より 2 % 以上悪くならない、NaN・DIVERGED なし、rms_roOmega が開始の 10 倍以下)。
  **Vb-dir** `run_0219_ns_coldmesh_tw300_linedir_tj5` = directional cfl 4 + キー 5、200 step・20 ごと・同じ 527 節点の帳簿。
  判定は V1 と同じ (200 step まで有限・帯の振動成分が最後の 32 step 窓 3 本で増えない・全残差と物理量が悪化しない)。帳簿で熱伝導の段と対流の段を V1 と並べる。
- **ノズルの結果 (2026-10-09)**: run_0214 (directional・キー 1・cfl 2) は 234 step で rms_ro が開始の 1343 倍、run_0215 (cfl 1) は 157 step で 1296 倍 → 発散の途中で止めた (**不合格**)。
  **V3 (run_0216、point cfl 4 + キー 1) は 922 step で NaN** (対照の run_0191 = キー 0 は 40000 step 安定) → **不合格**。最初の 200 step は残差が速く下がる (rms_roe 10.4 → 4.5) が、
  300 step 前後から運動量の残差が育ち、非有限は試験部 x_w 69.9〜70.5 の壁から 79〜89 層目 (超音速の中心部) で出た。
  仮説 (未確定): 新しいエネルギー行の ∂T/∂ρ = −(e − ½|u|²)/(ρc_v)・∂T/∂(ρu) = −u/(ρc_v) の交差項が高速域で大きく、近傍の結合を入れない対角だけの近似が高速域で行列の性質を崩す。
  **まとめ**: キー 1 は近壁の熱伝導のモード (1 層目) を抑えるが、ノズル全体では方向別でも point でも不安定化させる → **不採用**。V0 (既定でビット同一) は未実施 — コードを残すか戻すかの前に確認する。

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

## 参考: `solver_density_cuda/cuda_forge/timeIntegration_d.cu`

```
#include "timeIntegration_d.cuh"
#include "weakIsothermalWall_d.cuh"
#include "lowMachPrecond_d.cuh"   // Phase 4: β (lowMachBeta)・c' (lowMachCprime) device ヘルパ
#include "cuda_forge/eos_jacobian_d.cuh"  // 一般EOS固有系 (eos_split_jacobian_general_closed)。precond 経路で使用
#include "cuda_forge/block_dplur_jacobian_d.cuh"  // block_dplur::accumulate_split_jacobian_cf (共有ヘッダ)
#include <cstdlib>

// 診断 (env FORGE_AXIS_DIAG_ALPHA, 既定 0=不変): 近軸 (r→0) の半径方向音響モード安定化。
// roUy 対角に α·A_planar·c を加える。revolved 軸面積 (r_f·S→0) が落とす半径音響スペクトル半径を
// planar 面積 A_pl で補うもの。near-axis radial-momentum 不安定 (case/28 TP) の最小修正診断。
__device__ float g_axisDiagAlpha = 0.0f;

namespace block_dplur {

template<typename T>
__device__ __forceinline__ void zero5(T* vec)
{
    #pragma unroll
    for (int i = 0; i < 5; ++i) {
        vec[i] = 0.0;
    }
}

template<typename T>
__device__ __forceinline__ void add_identity_scaled(T mat[5][5], T scale)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        mat[row][row] += scale;
    }
}

// accumulate_split_jacobian_cf は cuda_forge/block_dplur_jacobian_d.cuh へ移設 (host/device 共有・Level3
// 単体テスト tools/test_eos_jacobian.cpp から呼ぶため)。block_dplur:: で参照する。

__device__ __forceinline__ void multiply_add_5x5_vec(
    const flow_float mat[5][5],
    const flow_float vec[5],
    flow_float out[5]
)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        flow_float sum = 0.0;
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            sum += mat[row][col] * vec[col];
        }
        out[row] += sum;
    }
}

__device__ __forceinline__ void copy5(const flow_float* src, flow_float* dst)
{
    #pragma unroll
    for (int i = 0; i < 5; ++i) {
        dst[i] = src[i];
    }
}

template<typename T>
__device__ __forceinline__ void zero5x5(T mat[5][5])
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            mat[row][col] = 0.0;
        }
    }
}

__device__ __forceinline__ void add_scaled_5x5(flow_float dst[5][5], const flow_float src[5][5], flow_float scale)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            dst[row][col] += scale * src[row][col];
        }
    }
}

__device__ __forceinline__ void build_abs_jacobian(
    flow_float gamma,
    flow_float nx,
    flow_float ny,
    flow_float nz,
    flow_float u,
    flow_float v,
    flow_float w,
    flow_float H,
    flow_float c,
    flow_float abs_jac[5][5]
)
{
    zero5x5(abs_jac);

    const flow_float sonic = max(c, static_cast<flow_float>(1.0e-8));
    const flow_float ek = 0.5 * (u * u + v * v + w * w);
    const flow_float U = u * nx + v * ny + w * nz;
    const flow_float chi = (gamma - 1.0) / sonic;
    const flow_float inv_sqrt2 = static_cast<flow_float>(0.7071067811865475244);

    flow_float lambda[5] = {
        fabs(U + sonic),
        fabs(U),
        fabs(U),
        fabs(U),
        fabs(U - sonic)
    };

    flow_float R[5][5];
    flow_float L[5][5];

    R[0][0] = inv_sqrt2 / sonic;             R[0][1] = ny / sonic;               R[0][2] = nz / sonic;               R[0][3] = nx / sonic;               R[0][4] = inv_sqrt2 / sonic;
    R[1][0] = (u / sonic + nx) * inv_sqrt2; R[1][1] = u * ny / sonic + nz;      R[1][2] = u * nz / sonic - ny;      R[1][3] = u * nx / sonic;           R[1][4] = (u / sonic - nx) * inv_sqrt2;
    R[2][0] = (v / sonic + ny) * inv_sqrt2; R[2][1] = v * ny / sonic;           R[2][2] = v * nz / sonic + nx;      R[2][3] = v * nx / sonic - nz;      R[2][4] = (v / sonic - ny) * inv_sqrt2;
    R[3][0] = (w / sonic + nz) * inv_sqrt2; R[3][1] = w * ny / sonic - nx;      R[3][2] = w * nz / sonic;           R[3][3] = w * nx / sonic + ny;      R[3][4] = (w / sonic - nz) * inv_sqrt2;
    R[4][0] = (ek / sonic + 1.0 / chi + U) * inv_sqrt2;
    R[4][1] = ek * ny / sonic + nz * u - nx * w;
    R[4][2] = ek * nz / sonic + nx * v - ny * u;
    R[4][3] = ek * nx / sonic + ny * w - nz * v;
    R[4][4] = (ek / sonic + 1.0 / chi - U) * inv_sqrt2;

    L[0][0] = (chi * ek - U) * inv_sqrt2;   L[0][1] = (-chi * u + nx) * inv_sqrt2; L[0][2] = (-chi * v + ny) * inv_sqrt2; L[0][3] = (-chi * w + nz) * inv_sqrt2; L[0][4] = chi * inv_sqrt2;
    L[1][0] = ny * (-chi * ek + sonic) - nz * u + nx * w;
    L[1][1] = ny * chi * u + nz;            L[1][2] = ny * chi * v;              L[1][3] = ny * chi * w - nx;         L[1][4] = -ny * chi;
    L[2][0] = nz * (-chi * ek + sonic) - nx * v + ny * u;
    L[2][1] = nz * chi * u - ny;            L[2][2] = nz * chi * v + nx;         L[2][3] = nz * chi * w;              L[2][4] = -nz * chi;
    L[3][0] = nx * (-chi * ek + sonic) - ny * w + nz * v;
    L[3][1] = nx * chi * u;                 L[3][2] = nx * chi * v - nz;         L[3][3] = nx * chi * w + ny;         L[3][4] = -nx * chi;
    L[4][0] = (chi * ek + U) * inv_sqrt2;   L[4][1] = (-chi * u - nx) * inv_sqrt2; L[4][2] = (-chi * v - ny) * inv_sqrt2; L[4][3] = (-chi * w - nz) * inv_sqrt2; L[4][4] = chi * inv_sqrt2;

    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            flow_float sum = 0.0;
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                sum += R[row][k] * lambda[k] * L[k][col];
            }
            abs_jac[row][col] = sum;
        }
    }
}

// LU-SGS の通量分割 Jacobian を同時構築する。
//   a_plus = A^+ = R Λ^+ L,  k_off = -A^- = R(-Λ^-)L = ½(|A|-A)。
// 検証済みの一般EOS閉形式 eos_split_jacobian_general_closed (eos_jacobian_d.cuh, Level1〜3 検証済) を流用し
// CPG/TP を統一 (旧来の CPG 専用べた打ち R/L=RL≠I の近似を廃止)。H は実全エンタルピー Ht[ic]、
// κ=γ−1、χ_eos=c²−κh。CPG では χ_eos≈0 で従来の CPG 固有系に簡約 (収束先は不変・厳密 Jacobian で僅かに向上)。
// double で組み立て float へ格納 (precond カーネルは元々 double 相当)。標準経路 accumulate_split_jacobian_cf は
// 別実装 (CPG ビット不変) のまま。
__device__ __forceinline__ void build_jacobian_split(
    flow_float gamma,
    flow_float nx, flow_float ny, flow_float nz,
    flow_float u, flow_float v, flow_float w,
    flow_float H, flow_float c,
    flow_float a_plus[5][5], flow_float k_off[5][5]
)
{
    const double sonic = (double)max(c, static_cast<flow_float>(1.0e-8));
    const double ek    = 0.5*((double)u*u + (double)v*v + (double)w*w);
    const double kappa = (double)gamma - 1.0;
    const double chi   = sonic*sonic - kappa*((double)H - ek);   // χ_eos = c²−κh (CPG で ≈0)
    double Ap[5][5], Am[5][5];
    eos_split_jacobian_general_closed((double)u,(double)v,(double)w, (double)nx,(double)ny,(double)nz,
                                      sonic, (double)H, kappa, chi, Ap, Am);   // Ap=A⁺, Am=A⁻
    #pragma unroll
    for (int i=0;i<5;++i)
        #pragma unroll
        for (int j=0;j<5;++j){ a_plus[i][j]=(flow_float)Ap[i][j]; k_off[i][j]=(flow_float)(-Am[i][j]); }
}

template<typename T>
__device__ __forceinline__ bool solve_5x5(T mat[5][5], T rhs[5], T sol[5])
{
    #pragma unroll
    for (int col = 0; col < 5; ++col) {
        int pivot = col;
        T pivot_abs = fabs(mat[col][col]);
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const T candidate = fabs(mat[row][col]);
            if (candidate > pivot_abs) {
                pivot = row;
                pivot_abs = candidate;
            }
        }

        if (pivot_abs < static_cast<T>(1.0e-20)) {
            zero5(sol);
            return false;
        }

        if (pivot != col) {
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                const T tmp = mat[col][k];
                mat[col][k] = mat[pivot][k];
                mat[pivot][k] = tmp;
            }
            const T rhs_tmp = rhs[col];
            rhs[col] = rhs[pivot];
            rhs[pivot] = rhs_tmp;
        }

        const T inv_pivot = static_cast<T>(1.0) / mat[col][col];
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const T factor = mat[row][col] * inv_pivot;
            mat[row][col] = 0.0;
            #pragma unroll
            for (int k = col + 1; k < 5; ++k) {
                mat[row][k] -= factor * mat[col][k];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    for (int row = 4; row >= 0; --row) {
        T sum = rhs[row];
        #pragma unroll
        for (int col = row + 1; col < 5; ++col) {
            sum -= mat[row][col] * sol[col];
        }
        sol[row] = sum / mat[row][row];
    }

    return true;
}

// 5×5 を 2 つの RHS について同時に解く (部分ピボット Gauss 消去を 1 回・float)。
// Phase 4 (lowMachPrecond=2) の Sherman-Morrison 解法で D0⁻¹b と D0⁻¹g を同じ分解で得るのに使う。
// D0 は物理ブロックで良条件 (既存 0/1 カーネルが float で解いているのと同形) ゆえ float で十分。
__device__ __forceinline__ bool solve_5x5_2rhs(flow_float mat[5][5],
                                               flow_float r1[5], flow_float r2[5],
                                               flow_float s1[5], flow_float s2[5])
{
    #pragma unroll
    for (int col = 0; col < 5; ++col) {
        int pivot = col;
        flow_float pivot_abs = fabs(mat[col][col]);
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const flow_float candidate = fabs(mat[row][col]);
            if (candidate > pivot_abs) { pivot = row; pivot_abs = candidate; }
        }

        if (pivot_abs < static_cast<flow_float>(1.0e-20)) {
            zero5(s1); zero5(s2);
            return false;
        }

        if (pivot != col) {
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                const flow_float tmp = mat[col][k]; mat[col][k] = mat[pivot][k]; mat[pivot][k] = tmp;
            }
            flow_float t = r1[col]; r1[col] = r1[pivot]; r1[pivot] = t;
            t = r2[col]; r2[col] = r2[pivot]; r2[pivot] = t;
        }

        const flow_float inv_pivot = static_cast<flow_float>(1.0) / mat[col][col];
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const flow_float factor = mat[row][col] * inv_pivot;
            mat[row][col] = 0.0;
            #pragma unroll
            for (int k = col + 1; k < 5; ++k) mat[row][k] -= factor * mat[col][k];
            r1[row] -= factor * r1[col];
            r2[row] -= factor * r2[col];
        }
    }

    for (int row = 4; row >= 0; --row) {
        flow_float sum1 = r1[row], sum2 = r2[row];
        #pragma unroll
        for (int col = row + 1; col < 5; ++col) { sum1 -= mat[row][col] * s1[col]; sum2 -= mat[row][col] * s2[col]; }
        const flow_float inv = static_cast<flow_float>(1.0) / mat[row][row];
        s1[row] = sum1 * inv;
        s2[row] = sum2 * inv;
    }

    return true;
}

__device__ __forceinline__ void load_block_vec(
    geom_int ic,
    flow_float* v0,
    flow_float* v1,
    flow_float* v2,
    flow_float* v3,
    flow_float* v4,
    flow_float out[5]
)
{
    out[0] = v0[ic];
    out[1] = v1[ic];
    out[2] = v2[ic];
    out[3] = v3[ic];
    out[4] = v4[ic];
}

__device__ __forceinline__ void store_block_vec(
    geom_int ic,
    const flow_float in[5],
    flow_float* v0,
    flow_float* v1,
    flow_float* v2,
    flow_float* v3,
    flow_float* v4
)
{
    v0[ic] = in[0];
    v1[ic] = in[1];
    v2[ic] = in[2];
    v3[ic] = in[3];
    v4[ic] = in[4];
}

}

__global__ void runge_kutta_exp_4th_d
// see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
// N: previous outer step , M: previous inner loop
( 
 int loop, 
 flow_float coef_DT,
 flow_float coef_Res,

 flow_float dt ,
 flow_float* dt_local ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol ,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* res_ro_m,
 flow_float* res_roUx_m,
 flow_float* res_roUy_m,
 flow_float* res_roUz_m,
 flow_float* res_roe_m

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    geom_float v = vol[ic];

    if (ic < nCells) {
        flow_float dt_l = dt_local[ic];

        if (loop == 0) {
            res_ro_m[ic]   = 0.0;
            res_roUx_m[ic] = 0.0;
            res_roUy_m[ic] = 0.0;
            res_roUz_m[ic] = 0.0;
            res_roe_m[ic]  = 0.0;
        }
        // N: previous outer step , M: previous inner loop
        res_ro_m[ic]   += coef_Res*res_ro[ic]*dt_l/v;
        res_roUx_m[ic] += coef_Res*res_roUx[ic]*dt_l/v;
        res_roUy_m[ic] += coef_Res*res_roUy[ic]*dt_l/v;
        res_roUz_m[ic] += coef_Res*res_roUz[ic]*dt_l/v;
        res_roe_m[ic]  += coef_Res*res_roe[ic]*dt_l/v;

        if (loop < 3) {
            ro[ic]   = roN[ic]   +coef_DT*res_ro[ic]*dt_l/v;
            roUx[ic] = roUxN[ic] +coef_DT*res_roUx[ic]*dt_l/v;
            roUy[ic] = roUyN[ic] +coef_DT*res_roUy[ic]*dt_l/v;
            roUz[ic] = roUzN[ic] +coef_DT*res_roUz[ic]*dt_l/v;
            roe[ic]  = roeN[ic]  +coef_DT*res_roe[ic]*dt_l/v;
        } else {
            res_ro[ic]   = res_ro_m[ic]*v/dt_l ;
            res_roUx[ic] = res_roUx_m[ic]*v/dt_l ;
            res_roUy[ic] = res_roUy_m[ic]*v/dt_l ;
            res_roUz[ic] = res_roUz_m[ic]*v/dt_l ;
            res_roe[ic]  = res_roe_m[ic]*v/dt_l ;

            ro[ic]   = roN[ic]  + res_ro_m[ic] ;
            roUx[ic] = roUxN[ic]+ res_roUx_m[ic] ;
            roUy[ic] = roUyN[ic]+ res_roUy_m[ic] ;
            roUz[ic] = roUzN[ic]+ res_roUz_m[ic] ;
            roe[ic]  = roeN[ic] + res_roe_m[ic] ;
        }
    }
}

__global__ void runge_kutta_exp_d
// see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
// N: previous outer step , M: previous inner loop
( 
 int loop,
 flow_float coef_N,
 flow_float coef_M,
 flow_float coef_Res,

 flow_float dt ,
 flow_float* dt_local ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol ,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        geom_float v = vol[ic];
        flow_float dt_l = dt_local[ic];
        // N: previous outer step , M: previous inner loop
        ro[ic]   = coef_N*roN[ic]   + coef_M*roM[ic]   + coef_Res*res_ro[ic]*dt_l/v;
        roUx[ic] = coef_N*roUxN[ic] + coef_M*roUxM[ic] + coef_Res*res_roUx[ic]*dt_l/v;
        roUy[ic] = coef_N*roUyN[ic] + coef_M*roUyM[ic] + coef_Res*res_roUy[ic]*dt_l/v;
        roUz[ic] = coef_N*roUzN[ic] + coef_M*roUzM[ic] + coef_Res*res_roUz[ic]*dt_l/v;
        roe[ic]  = coef_N*roeN[ic]  + coef_M*roeM[ic]  + coef_Res*res_roe[ic]*dt_l/v;
    }
}
__global__ void implicit_defect_correction_d
(
 int loop,
 flow_float dt,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)。軸対称ソース Jacobian 用

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol,
 geom_int* plane_cells,
 geom_int* cell_planes_index,
 geom_int* cell_planes,
 geom_float* ccx,
 geom_float* ccy,
 geom_float* ccz,
 geom_float* sx,
 geom_float* sy,
 geom_float* sz,
 geom_float* ss,

 // variables
 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float* roN,
 flow_float* roUxN,
 flow_float* roUyN,
 flow_float* roUzN,
 flow_float* roeN,

 flow_float laminar_visc,
 flow_float* vis_turb,
 flow_float* sonic,
 flow_float* Ux,
 flow_float* Uy,
 flow_float* Uz,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* corr_ro_old,
 flow_float* corr_roUx_old,
 flow_float* corr_roUy_old,
 flow_float* corr_roUz_old,
 flow_float* corr_roe_old,

 flow_float* corr_ro_new,
 flow_float* corr_roUx_new,
 flow_float* corr_roUy_new,
 flow_float* corr_roUz_new,
 flow_float* corr_roe_new,

 // 軸対称ソース項 Jacobian 用 (CPG/TP 共通)
 int isAxisymmetric,
 flow_float* A_planar,
 flow_float axisRFloor,
 // dual-time 物理時間項の対角係数 a/Δt（定常は 0）
 flow_float unsteady_diag
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        geom_float v = vol[ic];
        flow_float dt_l = dt_local[ic];
        const flow_float density = max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float velocity_x = Ux[ic];
        const flow_float velocity_y = Uy[ic];
        const flow_float velocity_z = Uz[ic];
        const flow_float local_sonic = max(sonic[ic], static_cast<flow_float>(0.0));
        const flow_float nu_eff = (laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;

        if (loop == 0) {
            corr_ro_old[ic] = 0.0;
            corr_roUx_old[ic] = 0.0;
            corr_roUy_old[ic] = 0.0;
            corr_roUz_old[ic] = 0.0;
            corr_roe_old[ic] = 0.0;
        }

        flow_float diag_face_sum = 0.0;
        flow_float neighbor_ro = 0.0;
        flow_float neighbor_roUx = 0.0;
        flow_float neighbor_roUy = 0.0;
        flow_float neighbor_roUz = 0.0;
        flow_float neighbor_roe = 0.0;
        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const flow_float face_area = max(ss[ip], static_cast<flow_float>(1.0e-30));
            const flow_float advective_radius = fabs(
                velocity_x * sx[ip] + velocity_y * sy[ip] + velocity_z * sz[ip]
            ) / face_area + local_sonic;

            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float dcc_x = ccx[other_ic] - ccx[ic];
            const flow_float dcc_y = ccy[other_ic] - ccy[ic];
            const flow_float dcc_z = ccz[other_ic] - ccz[ic];
            const flow_float dcc = max(
                sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float dcc_dot_s = max(
                fabs(dcc_x * sx[ip] + dcc_y * sy[ip] + dcc_z * sz[ip]),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<flow_float>(1.0e-30));
            // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる。
            // 旧 face_area·(2ν/delta) は面積 delta を長さ²扱いし ≈2ν に潰れ、軸対称近軸で r 重みを失い、
            // ゼロ面積(対称)面にもスプリアス項を載せていた (residual は ip<nNormalPlanes で軸面を除外)。
            const flow_float viscous_diag = static_cast<flow_float>(2.0) * nu_eff * delta / dcc;
            const flow_float face_coeff = face_area * advective_radius + viscous_diag;
            const flow_float offdiag_coeff = static_cast<flow_float>(0.5) * face_coeff;
            diag_face_sum += face_coeff;

            if (other_ic < nCells) {
                neighbor_ro += offdiag_coeff * corr_ro_old[other_ic];
                neighbor_roUx += offdiag_coeff * corr_roUx_old[other_ic];
                neighbor_roUy += offdiag_coeff * corr_roUy_old[other_ic];
                neighbor_roUz += offdiag_coeff * corr_roUz_old[other_ic];
                neighbor_roe += offdiag_coeff * corr_roe_old[other_ic];
            }
        }

        const flow_float diag = max(
            static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)))
                + static_cast<flow_float>(v) * unsteady_diag + diag_face_sum,
            static_cast<flow_float>(1.0e-30)
        );
        const flow_float inv_diag = 1.0 / diag;

        // 軸対称フープ源 (res_roUy += (P-τθθ)·A_planar, axisymmetricSource_d.cu) の Jacobian 対角成分を
        // roUy 方程式の対角に陰化する。これが無いと近軸の剛フープ源 (∝1/r) が陽 (lagged) 扱いになり、
        // block-DPLUR が安定な CFL でも scalar が発散する (block は diag_block[2][2] で陰化済。
        // 切り分けは case/29 README / plan time_integration-scalar-dplur-axisym-source.md)。
        // block 版 diag_block[2][2] と同形: A_pl·((γ-1)u_y + 2μ/(ρ r_eff))。
        // γ は per-cell gamma_arr[ic] (TP=γ_mix(T) / CPG=cfg.gamma) を使い thermally perfect でも整合。
        // scalar 対角の正値性 (対角優位) を保つため非負側のみ加える (defect-correction の不動点は不変)。
        flow_float diag_roUy = diag;
        if (isAxisymmetric == 1 &&
            !(axisRFloor > (flow_float)0.0 && ccy[ic] < axisRFloor)) {
            const flow_float A_pl = max(A_planar[ic], static_cast<flow_float>(1.0e-30));
            const flow_float r_eff = max(static_cast<flow_float>(v) / A_pl, static_cast<flow_float>(1.0e-30));
            const flow_float g1 = gamma_arr[ic] - static_cast<flow_float>(1.0);
            const flow_float mu_total = laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0));
            const flow_float hoop = static_cast<flow_float>(2.0) * mu_total / (density * r_eff);
            const flow_float src_diag = A_pl * (g1 * velocity_y + hoop);
            diag_roUy = diag + max(src_diag, static_cast<flow_float>(0.0));
        }
        const flow_float inv_diag_roUy = 1.0 / diag_roUy;

        const flow_float jacobi_ro = (res_ro[ic] + neighbor_ro) * inv_diag;
        const flow_float jacobi_roUx = (res_roUx[ic] + neighbor_roUx) * inv_diag;
        const flow_float jacobi_roUy = (res_roUy[ic] + neighbor_roUy) * inv_diag_roUy;
        const flow_float jacobi_roUz = (res_roUz[ic] + neighbor_roUz) * inv_diag;
        const flow_float jacobi_roe = (res_roe[ic] + neighbor_roe) * inv_diag;

        corr_ro_new[ic] = implicit_relax * jacobi_ro;
        corr_roUx_new[ic] = implicit_relax * jacobi_roUx;
        corr_roUy_new[ic] = implicit_relax * jacobi_roUy;
        corr_roUz_new[ic] = implicit_relax * jacobi_roUz;
        corr_roe_new[ic] = implicit_relax * jacobi_roe;
    }
}

// 5×5 行列を複数保持しレジスタ消費が大きいため、block 上限を超えないよう __launch_bounds__ で
// 1 block あたりスレッド数を 128 に制限する（起動時の "too many resources" を回避）。
#define BLOCK_DPLUR_THREADS 128
// 占有率実験用: __launch_bounds__ の最小常駐ブロック数 (既定 1 = 従来どおり制限なし, 128 regs → 占有率 ~28 %)。
// -DBLOCK_DPLUR_MINBLOCKS=n でビルドすると regs ≤ 65536/(128 n) に制限され (spill と引き換えに) 占有率が上がる。
#ifndef BLOCK_DPLUR_MINBLOCKS
#define BLOCK_DPLUR_MINBLOCKS 1
#endif
// block-DPLUR の閉形式 FVS 版。線形 solve の内部精度を ST (float 既定 / double で軸対称近軸を根治) で
// テンプレート化。残差/状態 (flow_float=float) を ST へキャストして取り込み、R/L を作らず閉形式で
// diag/nbr を畳み、ST で in-place 5×5 solve、補正を float dq_new へ書戻す (混合精度 iterative refinement)。
// 詳細: plans/archived/precision-mixed-axisym.md。
template<typename ST>
__global__ void __launch_bounds__(BLOCK_DPLUR_THREADS, BLOCK_DPLUR_MINBLOCKS) implicit_defect_correction_block_d
(
 int loop,
 flow_float dt,
 const flow_float* __restrict__ dt_local,
 flow_float implicit_relax,
 const flow_float* __restrict__ gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)。frozen-coefficient Jacobian 用
 int thermallyPerfect,    // 1: TP 固有系 (実 Ht・χ_eos=c²−κh, κ=γ−1), 0: CPG 閉形式 (従来・ビット不変)

 geom_int nCells_all , geom_int nCells,
 const geom_float* __restrict__ vol,
 const geom_int* __restrict__ plane_cells,
 const geom_int* __restrict__ cell_planes_index,
 const geom_int* __restrict__ cell_planes,
 const geom_float* __restrict__ ccx,
 const geom_float* __restrict__ ccy,
 const geom_float* __restrict__ ccz,
 const geom_float* __restrict__ sx,
 const geom_float* __restrict__ sy,
 const geom_float* __restrict__ sz,
 const geom_float* __restrict__ ss,

 const flow_float* __restrict__ ro,
 const flow_float* __restrict__ roUx,
 const flow_float* __restrict__ roUy,
 const flow_float* __restrict__ roUz,
 const flow_float* __restrict__ roe,

 flow_float laminar_visc,
 const flow_float* __restrict__ vis_turb,
 const flow_float* __restrict__ sonic,
 const flow_float* __restrict__ Ux,
 const flow_float* __restrict__ Uy,
 const flow_float* __restrict__ Uz,
 const flow_float* __restrict__ Ht,

 const flow_float* __restrict__ res_ro,
 const flow_float* __restrict__ res_roUx,
 const flow_float* __restrict__ res_roUy,
 const flow_float* __restrict__ res_roUz,
 const flow_float* __restrict__ res_roe,

 const flow_float* __restrict__ dq_old_0,
 const flow_float* __restrict__ dq_old_1,
 const flow_float* __restrict__ dq_old_2,
 const flow_float* __restrict__ dq_old_3,
 const flow_float* __restrict__ dq_old_4,

 flow_float* dq_new_0,
 flow_float* dq_new_1,
 flow_float* dq_new_2,
 flow_float* dq_new_3,
 flow_float* dq_new_4,

 flow_float* rhs_0,
 flow_float* rhs_1,
 flow_float* rhs_2,
 flow_float* rhs_3,
 flow_float* rhs_4,

 flow_float* diag_00, flow_float* diag_01, flow_float* diag_02, flow_float* diag_03, flow_float* diag_04,
 flow_float* diag_10, flow_float* diag_11, flow_float* diag_12, flow_float* diag_13, flow_float* diag_14,
 flow_float* diag_20, flow_float* diag_21, flow_float* diag_22, flow_float* diag_23, flow_float* diag_24,
 flow_float* diag_30, flow_float* diag_31, flow_float* diag_32, flow_float* diag_33, flow_float* diag_34,
 flow_float* diag_40, flow_float* diag_41, flow_float* diag_42, flow_float* diag_43, flow_float* diag_44,

 // 軸対称ソースヤコビアン用（isAxisymmetric==1 のときのみ使用）
 int isAxisymmetric,
 const flow_float* __restrict__ A_planar,

 // 軸対称 r 床 (axisymMethod==0): ccy < axisRFloor の帯は hoop ソース不課につき Jacobian も課さない。
 flow_float axisRFloor,

 // dual-time 物理時間項の対角係数 a/Δt（定常は 0）
 flow_float unsteady_diag,

 // node-centered 軸対称: 軸上 CV で半径方向運動量 (roUy, index2) 行を decouple する (nullptr 可)。
 // SU2 流の対称面を Jacobian 内で課す = solve の外で状態を手術せず一貫して dq_roUy=0 を得る。
 const geom_int* __restrict__ axis_flag,     // (未使用: 旧 nodeAxisDirichlet の全 5 行 decouple。常に nullptr)
 // node × 軸対称: 軸ノードで roUy 行 (index 2) のみ単位行化 (nullptr で無効)。
 const geom_int* __restrict__ axis_ur_flag,

 // axisymMethod==1 (isAxisymmetric enc==2) の軸ソース Jacobian ガード: 軸上ノード (==1) はソース 0 なので
 // Jacobian も加えない。decouple 用 axis_flag (nodeAxisDirichlet ゲート) とは独立に渡す (nullptr 可)。
 const geom_int* __restrict__ axis_flag_src,

 // node-centered 壁 no-slip: 壁ノードで運動量3行 (index1=roUx,2=roUy,3=roUz) を decouple する (nullptr 可)。
 // SU2 `DeleteValsRowi` 相当。残差射影だけでは block-DPLUR が壁運動量を連成したまま dq≠0 を返し速度 drift
 // するのを防ぐ。連続(行0)・エネルギー(行4)は保持。methods/discretization.md §7.2.1。
 const geom_int* __restrict__ wall_flag,

 // node-centered 等温壁: 壁ノードでエネルギー行 (index4=roe) を decouple する (nullptr 可)。
 // 壁ノード T ピン (applyNodeIsothermalWallPin / WMLES 等温 pin) と対。ピンで状態を上書きしながら
 // エネルギー行を連成したまま解くと Jacobian 不整合で発散する (2026-07-20 純伝導検証で実測)。
 const geom_int* __restrict__ iso_wall_flag,
 // 等温壁エネルギー境界の弱形式 (mesh.nodeIsothermalEnergyBC=1) の近似対角項の素材
 // g = Σ_{壁半割面} k_eff A_half / d_1 [W/K] (viscousFlux_wall_d が残差と同じ k_eff・幾何で積む)。
 // 厳密微分は対角 0 で温度微分は内部点の列にある。ここで足すのは SU2 型の**近似対角**である
 // (codex plan レビュー M1)。forge は残差微分の符号を反転して行列を組むので **+** で足す。
 // CPG: ΔA[4][4] = + g/(ρ c_v)。nullptr なら何もしない (既定はビット不変)。
 const flow_float* __restrict__ weakIsoDiag,
 flow_float weakIsoCv,

 // node-centered 弱形式 (Phase 2, 5e): node モードはゴーストセルを使わない。境界半割面 (has_nbr=false=ゴースト
 // 側) をこの node-to-node Jacobian ループから完全に除外する (continue)。境界ノードは物理境界上に乗るため
 // node→ghost が退化 (dcc≈0) し粘性対角 2ν·delta/dcc が爆発→対角巨大→dq≈0 で境界ノードが凍結する (出口 BL
 // 崩壊・残差プラトーの真因)。境界の対流/粘性は弱形式カーネルが残差側で担う。cell (isNode=0) は境界ゴースト
 // が法線方向に正しく置かれ非退化なので従来どおり境界面も処理する。
 int isNode,

 // --- line-implicit (plans/active/time_integration-line-implicit.md) ---
 // line_prev/next != nullptr で有効。ライン CV は点解せず、diag (storeLU 時)・rhs (毎 sweep)・
 // 近傍行列 K (storeLU 時) を保存して Thomas カーネルに委ねる。
 const geom_int* line_prev,
 const geom_int* line_next,
 flow_float* Kprev,
 flow_float* Knext,
 int storeLU,
 int lineViscCoupling,
 // 対角キャッシュ (plans/active/performance-3d-node-sst-speedup.md §4.2-4): 1 のとき loop==0 で組んだ対角 5×5
 // (拘束行の単位行化込み) を diag_** に保存し、loop>0 は対角組立・粘性対角・軸対称 Jacobian・近傍幾何読みを省略して
 // 保存値を読む。状態は sweep 中凍結なので結果はビット同一。線形 solve・rhs 拘束・近傍積は毎 sweep 従来どおり。
 // 呼び出し側で float・point 経路 (line_prev==nullptr) に限定する。
 int useDiagCache,
 // 近傍 dq の AoS 版 (stride 8 floats = 32 B セクタ整列, [0..4] を使用)。非 nullptr のとき近傍 gather は
 // dq_pack_old から float4+float の 2 ロード (SoA 5 配列の 5 ロード = 5 セクタから 1 セクタへ)。dq_pack_new には
 // dq_new と同じ値を書く。SoA の dq_new_* も従来どおり書く (commit・周期ミラー・診断が読む)。
 // line-implicit / node 周期 (SoA だけを書き換える経路) では呼び出し側が nullptr を渡す。
 const flow_float* __restrict__ dq_pack_old,
 flow_float* dq_pack_new,
 // エネルギー行の熱伝導 Jacobian (plans/active/time_integration-implicit-thermal-jacobian.md、ビットマスク、0 で従来どおり):
 //   ビット 1: 内部の node 間面でエネルギー行の粘性対角を k_face·δ/dcc·(γ/cp)·∂e/∂Q に置き換える (k_face は残差と同じ式)。
 //   ビット 2: 等温壁の節点のエネルギー行を拘束の行 [−e_w,0,0,0,1] にする。
 //   ビット 4 (ビット 1 と併用): 行 4 の従来のスカラー 2ν_eff·δ/dcc も残し、温度の項はその上に足す。
 // thermCondArr・cpArr・fxArr は thermalJac のビット 1 が立っているときだけ読む (それ以外は nullptr でよい)。
 int thermalJac,
 const flow_float* __restrict__ thermCondArr,
 const flow_float* __restrict__ cpArr,
 const flow_float* __restrict__ fxArr,
 flow_float Prt
)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        // 状態/残差 (flow_float=float) を solve 精度 ST へキャストして取り込む (混合精度)。
        const ST gamma = static_cast<ST>(gamma_arr[ic]);   // 局所 γ
        const ST v = static_cast<ST>(vol[ic]);
        const ST dt_l = static_cast<ST>(dt_local[ic]);
        const ST density = max(static_cast<ST>(ro[ic]), static_cast<ST>(1.0e-30));
        const ST velocity_x = static_cast<ST>(Ux[ic]);
        const ST velocity_y = static_cast<ST>(Uy[ic]);
        const ST velocity_z = static_cast<ST>(Uz[ic]);
        const ST local_sonic = max(static_cast<ST>(sonic[ic]), static_cast<ST>(1.0e-8));
        const ST local_Ht = static_cast<ST>(Ht[ic]);   // 一般EOS固有系のエネルギー成分 (TP)
        // line-implicit: 自 CV の役割 (ラインに載るか) と decouple 行マスク (保存 K の行ゼロ化用)
        const geom_int lp = (line_prev != nullptr) ? line_prev[ic] : (geom_int)(-1);
        const geom_int ln_ = (line_next != nullptr) ? line_next[ic] : (geom_int)(-1);
        const bool onLine = (lp >= 0) || (ln_ >= 0);
        // 対角キャッシュを読む sweep か (loop>0 かつ line に載らない CV)。
        const bool cached = (useDiagCache != 0) && (loop > 0) && !onLine;
        ST nu_eff = static_cast<ST>(0.0);
        if (!cached) nu_eff = (static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0))) / density;
        bool rowDec[5] = {false, false, false, false, false};
        if (onLine) {
            if (axis_ur_flag != nullptr && axis_ur_flag[ic] == 1) rowDec[2] = true;
            if (wall_flag != nullptr && wall_flag[ic] == 1) { rowDec[1] = true; rowDec[2] = true; rowDec[3] = true; }
            if (iso_wall_flag != nullptr && iso_wall_flag[ic] == 1) rowDec[4] = true;
        }

        // (loop==0 の dq_old ゼロ化は blockDPLURSolve の cudaMemset が担う。dq_old は本カーネルでは読み取り専用
        //  (const __restrict__) にして read-only キャッシュ経路を許す。dq_new とは別バッファ = 別名無し。)

        ST diag_block[5][5];
        block_dplur::zero5x5(diag_block);
        if (!cached) {
            block_dplur::add_identity_scaled(diag_block, static_cast<ST>(v / max(dt_l, static_cast<ST>(1.0e-30))));
            // dual-time: 物理時間項 a·V/Δt を対角へ（定常は unsteady_diag==0）。
            block_dplur::add_identity_scaled(diag_block, v * static_cast<ST>(unsteady_diag));
        }

        ST rhs[5] = {
            static_cast<ST>(res_ro[ic]),
            static_cast<ST>(res_roUx[ic]),
            static_cast<ST>(res_roUy[ic]),
            static_cast<ST>(res_roUz[ic]),
            static_cast<ST>(res_roe[ic])
        };
        ST neighbor_accum[5];
        block_dplur::zero5(neighbor_accum);

        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const ST face_area = max(static_cast<ST>(ss[ip]), static_cast<ST>(1.0e-30));
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;

            // 格納法線 (sx,sy,sz) は ic0→ic1。ic が neighbor 側 (ic1) のとき符号反転。
            const ST nsign = (ic0 == ic) ? static_cast<ST>(1.0) : static_cast<ST>(-1.0);
            const ST nx = nsign * static_cast<ST>(sx[ip]) / face_area;
            const ST ny = nsign * static_cast<ST>(sy[ip]) / face_area;
            const ST nz = nsign * static_cast<ST>(sz[ip]) / face_area;

            // 閉形式 FVS (R/L を作らず rank-2 外積で A⁺S を diag・k_off·sdq を nbr へ)。
            // 対流項: has_nbr=false (境界半割面) のとき自セル状態+面法線から A⁺S を対角に積むだけで、ゴースト
            // セルの状態/中心は一切読まない (ghostless)。よって node モードでもこの対流寄与は残す (境界ノードの
            // 流出 Jacobian=陰的安定化に必要。除くと rms_roUy 等が発散した)。
            const bool has_nbr = (other_ic < nCells);
            // ライン面: dq_old の lag 参照をスキップ (Thomas が厳密連成) — sdq=0 で対角 A⁺ だけ積む
            const bool isLineFace = onLine && has_nbr && (other_ic == lp || other_ic == ln_);
            ST sdq[5] = {static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0)};
            // loop==0 は dq_old≡0 (blockDPLURSolve の memset) なので gather を省く (寄与は厳密に 0 = ビット同一)。
            if (has_nbr && !isLineFace && loop > 0) {
                if (dq_pack_old != nullptr) {
                    const float4 q4 = *reinterpret_cast<const float4*>(dq_pack_old + (size_t)other_ic * 8);
                    const flow_float q5 = dq_pack_old[(size_t)other_ic * 8 + 4];
                    sdq[0] = face_area * static_cast<ST>(q4.x);
                    sdq[1] = face_area * static_cast<ST>(q4.y);
                    sdq[2] = face_area * static_cast<ST>(q4.z);
                    sdq[3] = face_area * static_cast<ST>(q4.w);
                    sdq[4] = face_area * static_cast<ST>(q5);
                } else {
                    sdq[0] = face_area * static_cast<ST>(dq_old_0[other_ic]);
                    sdq[1] = face_area * static_cast<ST>(dq_old_1[other_ic]);
                    sdq[2] = face_area * static_cast<ST>(dq_old_2[other_ic]);
                    sdq[3] = face_area * static_cast<ST>(dq_old_3[other_ic]);
                    sdq[4] = face_area * static_cast<ST>(dq_old_4[other_ic]);
                }
            }
            if (cached) {
                block_dplur::accumulate_split_jacobian_cf<ST, false>(
                    gamma, nx, ny, nz, velocity_x, velocity_y, velocity_z,
                    local_sonic, local_Ht, thermallyPerfect != 0,
                    face_area, has_nbr, sdq, diag_block, neighbor_accum
                );
            } else {
                block_dplur::accumulate_split_jacobian_cf<ST>(
                    gamma, nx, ny, nz, velocity_x, velocity_y, velocity_z,
                    local_sonic, local_Ht, thermallyPerfect != 0,
                    face_area, has_nbr, sdq, diag_block, neighbor_accum
                );
            }
            // K 行列の列抽出 (状態凍結ゆえ storeLU=loop0 のみ): nbr 寄与は sdq に線形なので
            // 単位ベクトル×face_area で列が得られる (対角へは dummy に捨てる)。decouple 行は 0。
            if (isLineFace && storeLU != 0) {
                flow_float* Kdst = (other_ic == lp) ? Kprev : Knext;
                ST ddum[5][5];
                ST kcol[5];
                ST evec[5];
                for (int j = 0; j < 5; ++j) {
                    block_dplur::zero5x5(ddum);
                    block_dplur::zero5(kcol);
                    #pragma unroll
                    for (int q = 0; q < 5; ++q) evec[q] = static_cast<ST>(0.0);
                    evec[j] = face_area;
                    block_dplur::accumulate_split_jacobian_cf<ST>(
                        gamma, nx, ny, nz, velocity_x, velocity_y, velocity_z,
                        local_sonic, local_Ht, thermallyPerfect != 0,
                        face_area, true, evec, ddum, kcol);
                    #pragma unroll
                    for (int i = 0; i < 5; ++i)
                        Kdst[(size_t)ic * 25 + i * 5 + j] =
                            rowDec[i] ? (flow_float)0.0 : static_cast<flow_float>(kcol[i]);
                }
            }

            // 粘性対角: node モードはゴーストセルを使わない。境界半割面 (has_nbr=false=ゴースト側) では
            // dcc 計算 (ccx[ghost] 読み) も viscous_diag も行わない。境界ノードは境界面上に乗るため node→ghost が
            // 退化 (dcc≈0) し 2ν·delta/dcc が爆発→対角巨大→dq≈0 で境界ノードが凍結する (出口 BL 崩壊・残差
            // プラトーの真因)。境界粘性は弱形式カーネルが残差側で担う。内部 node-to-node 面のみ粘性対角を課す。
            // cell モード (isNode=0) は境界ゴーストが法線方向に正しく置かれ非退化なので従来どおり境界面も課す。
            if (!cached && !(isNode != 0 && !has_nbr)) {
                const ST dcc_x = static_cast<ST>(ccx[other_ic]) - static_cast<ST>(ccx[ic]);
                const ST dcc_y = static_cast<ST>(ccy[other_ic]) - static_cast<ST>(ccy[ic]);
                const ST dcc_z = static_cast<ST>(ccz[other_ic]) - static_cast<ST>(ccz[ic]);
                const ST dcc = max(sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z), static_cast<ST>(1.0e-30));
                const ST dcc_dot_s = max(
                    fabs(dcc_x * static_cast<ST>(sx[ip]) + dcc_y * static_cast<ST>(sy[ip]) + dcc_z * static_cast<ST>(sz[ip])),
                    static_cast<ST>(1.0e-30)
                );
                const ST delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<ST>(1.0e-30));
                // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる
                // (旧 face_area·(2ν/delta) は ≈2ν に潰れ近軸で r 重み喪失・ゼロ面積面にスプリアス。詳細は site1 コメント)。
                const ST viscous_diag = static_cast<ST>(2.0) * nu_eff * delta / dcc;
                if (isLineFace && lineViscCoupling != 0) {
                    // v2 (plans/active/time_integration-line-implicit-viscous-v2.md): line 面は
                    // スカラー粘性結合 K += α·I (α=ν_eff·δ/dcc) と対にし、対角は 2α→α に置換して
                    // 真の 1D 拡散行 [−α, 2α, −α] を line 内で完成させる (off-line 面は従来 2α のまま)。
                    const ST alpha = static_cast<ST>(0.5) * viscous_diag;
                    block_dplur::add_identity_scaled(diag_block, alpha);
                    if (storeLU != 0) {
                        flow_float* Kdst = (other_ic == lp) ? Kprev : Knext;
                        #pragma unroll
                        for (int i = 0; i < 5; ++i)
                            if (!rowDec[i]) Kdst[(size_t)ic * 25 + i * 5 + i] += static_cast<flow_float>(alpha);
                    }
                } else if ((thermalJac & 1) != 0 && isNode != 0 && has_nbr) {
                    // 熱伝導の Jacobian (implicitThermalJacobian ビット 1): 連続・運動量の行は従来どおりスカラー、
                    // エネルギー行は熱伝導の残差 k_face·(T_j − T_i)·δ/dcc の Q_i による微分の符号反転 = Λ^T·(γ/cp)·∂e/∂Q。
                    // k_face は viscousFlux_d.cu の tc_face と同じ式 (f 補間の層流 k + 面 cp × 面 μ_t / Pr_t)、物性・γ は凍結。
                    // ∂e/∂ρ = −(e − ½|u|²)/ρ、∂e/∂(ρu_k) = −u_k/ρ、∂e/∂(ρE) = 1/ρ (e = ρE/ρ − ½|u|²; TP でも Y は正規化済みなので厳密)。
                    #pragma unroll
                    // ビット 4: 行 4 にも従来のスカラー (スペクトル半径の近似) を残し、その上に温度の項を足す。
                    // 行 0〜3 と同じ一律の減衰を行 4 にも保つので、行ごとに歩幅の縮め方が違う状態 (温度を変えない
                    // 補正でエネルギー行だけ減衰が抜け、高速域で偽の ΔT を作る) を避ける (plan §4.3)。
                    const int nScalarRows = ((thermalJac & 4) != 0) ? 5 : 4;
                    for (int r = 0; r < nScalarRows; ++r) diag_block[r][r] += viscous_diag;
                    const ST f = static_cast<ST>(fxArr[ip]);
                    const ST omf = static_cast<ST>(1.0) - f;
                    const ST cp_face = f * static_cast<ST>(cpArr[ic0]) + omf * static_cast<ST>(cpArr[ic1]);
                    const ST mut_face = f * static_cast<ST>(vis_turb[ic0]) + omf * static_cast<ST>(vis_turb[ic1]);
                    const ST k_face = f * static_cast<ST>(thermCondArr[ic0]) + omf * static_cast<ST>(thermCondArr[ic1])
                                    + cp_face * mut_face / static_cast<ST>(Prt);
                    const ST q2 = velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z;
                    const ST e_int = static_cast<ST>(roe[ic]) / density - static_cast<ST>(0.5) * q2;
                    const ST cfac = (k_face * delta / dcc) * (gamma / max(static_cast<ST>(cpArr[ic]), static_cast<ST>(1.0e-30))) / density;
                    diag_block[4][0] += -cfac * (e_int - static_cast<ST>(0.5) * q2);
                    diag_block[4][1] += -cfac * velocity_x;
                    diag_block[4][2] += -cfac * velocity_y;
                    diag_block[4][3] += -cfac * velocity_z;
                    diag_block[4][4] += cfac;
                } else {
                    block_dplur::add_identity_scaled(diag_block, viscous_diag);
                    // **診断専用 A/B** (codex 2026-09-21、既定はコンパイルから除外されビット不変)。
                    // 粘性対角 2ν_eff·delta/dcc は `add_identity_scaled` で**5 行すべてに同じ量**が入る。
                    // エネルギー行の真の拡散 Jacobian はこれと違う (完全気体で rho,rho u を固定すると
                    // ∂T/∂(ρE)=1/(ρc_v) なので k∂T/∂(ρE)=γα、さらに交差微分がある) が、
                    // **不足率は α/ν=1.11-1.39 からは証明できない** (codex 指摘)。
                    // ここで調べるのは「非収束状態の 2 節点指標がエネルギー対角に感度を持つか」だけで、
                    // 1.4 倍は物理係数の再現ではなく試験強度である。
                    // 効いた場合も「Pr 補正が正しい」ではなく「指標が陰解法作用素に依存する」と結論する。
                    // plan boundary-conjugate-heat-transfer §5.1 #43。
#if defined(FORGE_TEST_ENERGY_VISCOUS_DIAG)
                    diag_block[4][4] += static_cast<ST>(0.4) * viscous_diag;
#endif
                }
            }
        }

        #pragma unroll
        for (int i = 0; i < 5; ++i) {
            rhs[i] += neighbor_accum[i];
        }

        // 軸対称ソース項のヤコビアンを対角ブロックに加える（roUy 行 = index 2）。詳細は実装ドキュメント参照。
        // axisRFloor 帯 (r 床, ソース不課) は Jacobian も課さない。
        if (!cached && isAxisymmetric == 1 &&
            !(static_cast<ST>(axisRFloor) > static_cast<ST>(0.0) && static_cast<ST>(ccy[ic]) < static_cast<ST>(axisRFloor))) {
            const ST A_pl = static_cast<ST>(A_planar[ic]);
            const ST r_eff = max(v / max(A_pl, static_cast<ST>(1.0e-30)), static_cast<ST>(1.0e-30));
            const ST g1 = gamma - static_cast<ST>(1.0);
            const ST q2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
            const ST mu_total = static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0));
            const ST hoop = static_cast<ST>(2.0) * mu_total / (density * r_eff);
            diag_block[2][0] += -A_pl * (static_cast<ST>(0.5)*g1*q2 + hoop * velocity_y);
            diag_block[2][1] += A_pl * (g1 * velocity_x);
            diag_block[2][2] += A_pl * (g1 * velocity_y + hoop);
            diag_block[2][3] += A_pl * (g1 * velocity_z);
            diag_block[2][4] += -A_pl * g1;
            // 診断: 近軸半径音響スペクトル半径 α·A_pl·c を roUy 対角に補う (FORGE_AXIS_DIAG_ALPHA>0 のみ)。
            diag_block[2][2] += static_cast<ST>(g_axisDiagAlpha) * A_pl * local_sonic;
        } else if (!cached && isAxisymmetric == 2) {
            // SU2 流 (axisymMethod==1) 非粘性軸対称ソースの解析 Jacobian (CSourceAxisymmetric_Flow 移植,
            // 行/列 = [ro, roUx, roUy, roe] → forge [0,1,2,4])。forge 対角は -∂S/∂U = +SU2 jacobian。
            // 軸ノード (axis_flag_src==1) と y≤eps はソース 0 のためスキップ。γ は frozen (gamma_arr)。
            const ST y = static_cast<ST>(ccy[ic]);
            const bool onAxisSrc = (axis_flag_src != nullptr && axis_flag_src[ic] == 1);
            if (!onAxisSrc && y > static_cast<ST>(1.0e-12)) {
                const ST yv = v / y;
                const ST g1 = gamma - static_cast<ST>(1.0);
                const ST uu = velocity_x, ww = velocity_y;
                const ST q2d = uu*uu + ww*ww;
                const ST et = static_cast<ST>(roe[ic]) / density;   // 比全エネルギー
                diag_block[0][2] += yv;
                diag_block[1][0] += yv * (-uu * ww);
                diag_block[1][1] += yv * ww;
                diag_block[1][2] += yv * uu;
                diag_block[2][0] += yv * (-ww * ww);
                diag_block[2][2] += yv * static_cast<ST>(2.0) * ww;
                diag_block[4][0] += yv * (-gamma * ww * et + g1 * ww * q2d);
                diag_block[4][1] += yv * (-g1 * uu * ww);
                diag_block[4][2] += yv * (gamma * et - static_cast<ST>(0.5) * g1 * (q2d + static_cast<ST>(2.0) * ww * ww));
                diag_block[4][4] += yv * (gamma * ww);
                // 粘性軸対称ソースの stiff 主対角: S_roUy ∋ -V·2μ_tot·v/y² → -∂S/∂(ρv) = +V·2μ/(ρy²)。
                // 近軸第一列 (y~1e-4) で極めて stiff で、これを lag すると implicit が喉部近軸で
                // limit cycle 化し rms_ro ~1e-5 で頭打ちになる (explicit は 3e-7 到達 = 空間は健全)。
                const ST mu_tot_ax = static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0));
                diag_block[2][2] += yv * static_cast<ST>(2.0) * mu_tot_ax / (density * y);
            }
        }

        // node × 軸対称: 軸ノードの半径運動量行のみ decouple (dq_roUy=0)。状態は enforceAxisSymmetry がピン。
        if (axis_ur_flag != nullptr && axis_ur_flag[ic] == 1) {
            if (!cached) {
                for (int jj = 0; jj < 5; ++jj) diag_block[2][jj] = static_cast<ST>(0.0);
                diag_block[2][2] = static_cast<ST>(1.0);
            }
            rhs[2] = static_cast<ST>(0.0);
        }

        // SU2 `DeleteValsRowi` 相当の壁 no-slip Dirichlet: 壁ノードで運動量3行 (index 1,2,3) を単位行に
        // 置換し rhs=0 → solve が一貫して dq_roUx=dq_roUy=dq_roUz=0 を返す。連続(0)・エネルギー(4)行は
        // 保持され ρ,ρe は保存式で発展、圧力は EOS が復元 (CPG/TP 共通)。残差射影だけでは block-DPLUR が
        // 壁運動量を連成し dq≠0 を返して壁速度が drift する問題を Jacobian 整合で根治する。
        if (wall_flag != nullptr && wall_flag[ic] == 1) {
            for (int row = 1; row <= 3; ++row) {
                if (!cached) {
                    for (int jj = 0; jj < 5; ++jj) diag_block[row][jj] = static_cast<ST>(0.0);
                    diag_block[row][row] = static_cast<ST>(1.0);
                }
                rhs[row] = static_cast<ST>(0.0);
            }
        }

        // 弱形式の等温壁 (nodeIsothermalEnergyBC=1): エネルギー行は残したまま、壁寄与の近似対角を足す。
        // iso_wall_flag が nullptr になっているので下の単位行化とは排他。
        if (weakIsoDiag != nullptr && !cached) {
            const ST g = static_cast<ST>(weakIsoDiag[ic]);
            if (g > static_cast<ST>(0.0)) {
                const ST rho = static_cast<ST>(max(ro[ic], (flow_float)1.0e-30));
                diag_block[4][4] += g / (rho * static_cast<ST>(weakIsoCv));
            }
        }

        // 等温壁ノード: エネルギー行 (4) も単位行に置換し dq_roe=0 → 壁ノード T は pin (applyBconds 位相) が
        // 一意に決める。連続 (0) 行は保持 (ρ は保存式で発展し P=ρRTw が追従)。
        if (iso_wall_flag != nullptr && iso_wall_flag[ic] == 1) {
            if (!cached) {
                for (int jj = 0; jj < 5; ++jj) diag_block[4][jj] = static_cast<ST>(0.0);
                diag_block[4][4] = static_cast<ST>(1.0);
                // implicitThermalJacobian ビット 2: 拘束の行 Δ(ρE)_w − e_w·Δρ_w = 0 (壁温のピン ρE = ρ·e(T_w) と一致、壁は u = 0)。
                if ((thermalJac & 2) != 0) diag_block[4][0] = -static_cast<ST>(roe[ic]) / density;
            }
            rhs[4] = static_cast<ST>(0.0);
        }

        // 対角キャッシュ: loop>0 は保存値を読む / loop==0 (useDiagCache かつ line 外) は組んだ対角を保存する。
        // ST=float・point 経路に限定して呼ばれる (呼び出し側ゲート) ので、保存/読込で丸めは発生しない (ビット同一)。
        if (cached) {
            diag_block[0][0]=static_cast<ST>(diag_00[ic]); diag_block[0][1]=static_cast<ST>(diag_01[ic]); diag_block[0][2]=static_cast<ST>(diag_02[ic]); diag_block[0][3]=static_cast<ST>(diag_03[ic]); diag_block[0][4]=static_cast<ST>(diag_04[ic]);
            diag_block[1][0]=static_cast<ST>(diag_10[ic]); diag_block[1][1]=static_cast<ST>(diag_11[ic]); diag_block[1][2]=static_cast<ST>(diag_12[ic]); diag_block[1][3]=static_cast<ST>(diag_13[ic]); diag_block[1][4]=static_cast<ST>(diag_14[ic]);
            diag_block[2][0]=static_cast<ST>(diag_20[ic]); diag_block[2][1]=static_cast<ST>(diag_21[ic]); diag_block[2][2]=static_cast<ST>(diag_22[ic]); diag_block[2][3]=static_cast<ST>(diag_23[ic]); diag_block[2][4]=static_cast<ST>(diag_24[ic]);
            diag_block[3][0]=static_cast<ST>(diag_30[ic]); diag_block[3][1]=static_cast<ST>(diag_31[ic]); diag_block[3][2]=static_cast<ST>(diag_32[ic]); diag_block[3][3]=static_cast<ST>(diag_33[ic]); diag_block[3][4]=static_cast<ST>(diag_34[ic]);
            diag_block[4][0]=static_cast<ST>(diag_40[ic]); diag_block[4][1]=static_cast<ST>(diag_41[ic]); diag_block[4][2]=static_cast<ST>(diag_42[ic]); diag_block[4][3]=static_cast<ST>(diag_43[ic]); diag_block[4][4]=static_cast<ST>(diag_44[ic]);
        } else if (useDiagCache != 0 && !onLine) {
            diag_00[ic]=static_cast<flow_float>(diag_block[0][0]); diag_01[ic]=static_cast<flow_float>(diag_block[0][1]); diag_02[ic]=static_cast<flow_float>(diag_block[0][2]); diag_03[ic]=static_cast<flow_float>(diag_block[0][3]); diag_04[ic]=static_cast<flow_float>(diag_block[0][4]);
            diag_10[ic]=static_cast<flow_float>(diag_block[1][0]); diag_11[ic]=static_cast<flow_float>(diag_block[1][1]); diag_12[ic]=static_cast<flow_float>(diag_block[1][2]); diag_13[ic]=static_cast<flow_float>(diag_block[1][3]); diag_14[ic]=static_cast<flow_float>(diag_block[1][4]);
            diag_20[ic]=static_cast<flow_float>(diag_block[2][0]); diag_21[ic]=static_cast<flow_float>(diag_block[2][1]); diag_22[ic]=static_cast<flow_float>(diag_block[2][2]); diag_23[ic]=static_cast<flow_float>(diag_block[2][3]); diag_24[ic]=static_cast<flow_float>(diag_block[2][4]);
            diag_30[ic]=static_cast<flow_float>(diag_block[3][0]); diag_31[ic]=static_cast<flow_float>(diag_block[3][1]); diag_32[ic]=static_cast<flow_float>(diag_block[3][2]); diag_33[ic]=static_cast<flow_float>(diag_block[3][3]); diag_34[ic]=static_cast<flow_float>(diag_block[3][4]);
            diag_40[ic]=static_cast<flow_float>(diag_block[4][0]); diag_41[ic]=static_cast<flow_float>(diag_block[4][1]); diag_42[ic]=static_cast<flow_float>(diag_block[4][2]); diag_43[ic]=static_cast<flow_float>(diag_block[4][3]); diag_44[ic]=static_cast<flow_float>(diag_block[4][4]);
        }

        if (onLine) {
            // line-implicit: 点解せず Thomas カーネル用に保存する。
            //   diag: 状態凍結ゆえ storeLU (loop==0) のみ / rhs: ライン外 lag 込みなので毎 sweep。
            //   dq_new は Thomas が上書きする (保険で前回反復値を置く)。
            if (storeLU != 0) {
                diag_00[ic]=static_cast<flow_float>(diag_block[0][0]); diag_01[ic]=static_cast<flow_float>(diag_block[0][1]); diag_02[ic]=static_cast<flow_float>(diag_block[0][2]); diag_03[ic]=static_cast<flow_float>(diag_block[0][3]); diag_04[ic]=static_cast<flow_float>(diag_block[0][4]);
                diag_10[ic]=static_cast<flow_float>(diag_block[1][0]); diag_11[ic]=static_cast<flow_float>(diag_block[1][1]); diag_12[ic]=static_cast<flow_float>(diag_block[1][2]); diag_13[ic]=static_cast<flow_float>(diag_block[1][3]); diag_14[ic]=static_cast<flow_float>(diag_block[1][4]);
                diag_20[ic]=static_cast<flow_float>(diag_block[2][0]); diag_21[ic]=static_cast<flow_float>(diag_block[2][1]); diag_22[ic]=static_cast<flow_float>(diag_block[2][2]); diag_23[ic]=static_cast<flow_float>(diag_block[2][3]); diag_24[ic]=static_cast<flow_float>(diag_block[2][4]);
                diag_30[ic]=static_cast<flow_float>(diag_block[3][0]); diag_31[ic]=static_cast<flow_float>(diag_block[3][1]); diag_32[ic]=static_cast<flow_float>(diag_block[3][2]); diag_33[ic]=static_cast<flow_float>(diag_block[3][3]); diag_34[ic]=static_cast<flow_float>(diag_block[3][4]);
                diag_40[ic]=static_cast<flow_float>(diag_block[4][0]); diag_41[ic]=static_cast<flow_float>(diag_block[4][1]); diag_42[ic]=static_cast<flow_float>(diag_block[4][2]); diag_43[ic]=static_cast<flow_float>(diag_block[4][3]); diag_44[ic]=static_cast<flow_float>(diag_block[4][4]);
            }
            rhs_0[ic] = static_cast<flow_float>(rhs[0]);
            rhs_1[ic] = static_cast<flow_float>(rhs[1]);
            rhs_2[ic] = static_cast<flow_float>(rhs[2]);
            rhs_3[ic] = static_cast<flow_float>(rhs[3]);
            rhs_4[ic] = static_cast<flow_float>(rhs[4]);
            dq_new_0[ic] = dq_old_0[ic];
            dq_new_1[ic] = dq_old_1[ic];
            dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic];
            dq_new_4[ic] = dq_old_4[ic];
        } else {
        // diag_block を破壊して in-place で解く (solve_mat コピー排除)。
        ST correction[5] = {static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0)};
        const bool ok = block_dplur::solve_5x5(diag_block, rhs, correction);
        if (!ok) {
            block_dplur::zero5(correction);
        }

        const ST relax = static_cast<ST>(implicit_relax);
        #pragma unroll
        for (int i = 0; i < 5; ++i) {
            correction[i] *= relax;
        }

        // 古典 DPLUR: float dq_new へ書戻し。Q への commit は applyBlockImplicitCorrection。
        dq_new_0[ic] = static_cast<flow_float>(correction[0]);
        dq_new_1[ic] = static_cast<flow_float>(correction[1]);
        dq_new_2[ic] = static_cast<flow_float>(correction[2]);
        dq_new_3[ic] = static_cast<flow_float>(correction[3]);
        dq_new_4[ic] = static_cast<flow_float>(correction[4]);
        if (dq_pack_new != nullptr) {
            *reinterpret_cast<float4*>(dq_pack_new + (size_t)ic * 8) =
                make_float4(static_cast<flow_float>(correction[0]), static_cast<flow_float>(correction[1]),
                            static_cast<flow_float>(correction[2]), static_cast<flow_float>(correction[3]));
            dq_pack_new[(size_t)ic * 8 + 4] = static_cast<flow_float>(correction[4]);
        }
        // rhs_** の診断書き出しは撤去 (読者なし。5 配列×sweep の書込 ≈240 MB/step を節約, 2026-09-12)。
        // line 経路 (上の onLine 分岐) は Thomas カーネルが rhs を読むので従来どおり書く。
        }
    }
}

// =============================================================================
// Phase 4 (a): 完全 Γ⁻¹A 低マッハ前処理の block DPLUR (lowMachPrecond>=2; 2=RHS+LHS / 3=LHS-only)。
// 既存 implicit_defect_correction_block_d (lowMachPrecond 0/1) とは別カーネルにし、
// 0/1 経路のレジスタ・ビットを一切変えない。
// 保存形は (Γ_c·V/Δτ' + A_c)ΔQ = -R で、**前処理は擬似時間項 Γ_c のみ**。フラックス A_c は
// 既存と同じ物理厳密 FVS (a_plus=A_c⁺・k_off=-A_c⁻) をそのまま使う (非前処理が正しい)。収束は
// (Δτ'/V)Γ_c⁻¹A_c の固有値 λ' で一様に前処理され、スカラー Δτ'=cell/ρ' で効く。
//   - 擬似時間項: Γ_c·V/Δτ' (Δτ' は setDT 側で前処理スペクトル半径から拡大した dt_local)。
//   - フラックス: 物理 a_plus を対角・k_off を近傍 (既存 block と同一)。
//   - 物理 BDF 項 a·V/Δt·I は非前処理 (dual-time 所有)。
// **Sherman-Morrison 解法**: Γ_c=I+α g rᵀ がランク1なので D=D0+γ g rᵀ (D0=物理ブロック・良条件)。
//   D0 を float で 2 RHS 同時 (solve_5x5_2rhs) に解き、悪条件 ~1/β は分母スカラーのみ double に隔離。
//   FP64 を回避して 0/1 カーネルに近い速度。β=1 (超音速) で Γ_c=I・Δτ'=Δτ・フラックス同一ゆえ現行と解一致。
// 理論: methods/time_integration/theory.md「低マッハ前処理固有系」、計画 §5 Phase 4。
__global__ void __launch_bounds__(BLOCK_DPLUR_THREADS) implicit_defect_correction_block_precond_d
(
 int loop,
 flow_float dt,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)
 flow_float precondEps,

 geom_int nCells_all , geom_int nCells,
 geom_float* vol,
 geom_int* plane_cells,
 geom_int* cell_planes_index,
 geom_int* cell_planes,
 geom_float* ccx, geom_float* ccy, geom_float* ccz,
 geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,

 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,

 flow_float laminar_visc,
 flow_float* vis_turb,
 flow_float* sonic,
 flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* Ht,

 flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe,

 flow_float* dq_old_0, flow_float* dq_old_1, flow_float* dq_old_2, flow_float* dq_old_3, flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,

 int isAxisymmetric,
 flow_float* A_planar,
 flow_float axisRFloor,
 flow_float unsteady_diag
)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float gamma = gamma_arr[ic];   // 局所 γ (TP の γ_mix。CPG は cfg.gamma で不変)
        const geom_float v = vol[ic];
        const flow_float dt_l = dt_local[ic];
        const flow_float density = max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float vx = Ux[ic];
        const flow_float vy = Uy[ic];
        const flow_float vz = Uz[ic];
        const flow_float local_sonic = max(sonic[ic], static_cast<flow_float>(1.0e-8));
        const flow_float local_enthalpy = max(Ht[ic], static_cast<flow_float>(1.0e-8));
        const flow_float nu_eff = (laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;
        const flow_float velMag = sqrt(vx*vx + vy*vy + vz*vz);
        const flow_float beta = lowMachBeta(local_sonic, velMag, precondEps);

        if (loop == 0) {
            dq_old_0[ic] = 0.0; dq_old_1[ic] = 0.0; dq_old_2[ic] = 0.0; dq_old_3[ic] = 0.0; dq_old_4[ic] = 0.0;
        }

        // 対角ブロックは Γ_c=I+α g rᵀ がランク1ゆえ D = D0 + γ g rᵀ と書ける:
        //   D0 = V/Δτ'·I + a·V/Δt·I + Σ A_c⁺ S + 粘性 + 軸対称  (物理ブロック・良条件・float 可)
        //   γ g rᵀ = (V/Δτ')·α g rᵀ                             (Γ_c 前処理寄与、悪条件 ~1/β の源)
        // Sherman-Morrison: x = y - [γ(rᵀy)/(1+γ(rᵀz))] z,  y=D0⁻¹b, z=D0⁻¹g。
        // D0 を float で 2 RHS 同時に解き、悪条件は分母スカラー(double)に隔離 → FP64 を回避 (RTX 等で高速)。
        flow_float D0[5][5];
        block_dplur::zero5x5(D0);
        const flow_float v_over_dtau = static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)));
        block_dplur::add_identity_scaled(D0, v_over_dtau);
        block_dplur::add_identity_scaled(D0, static_cast<flow_float>(v) * unsteady_diag);  // dual-time BDF (非前処理)

        flow_float b[5] = { res_ro[ic], res_roUx[ic], res_roUy[ic], res_roUz[ic], res_roe[ic] };
        flow_float nbr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const flow_float face_area = max(ss[ip], static_cast<flow_float>(1.0e-30));
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float nsign = (ic0 == ic) ? static_cast<flow_float>(1.0) : static_cast<flow_float>(-1.0);
            const flow_float nx = nsign * sx[ip] / face_area;
            const flow_float ny = nsign * sy[ip] / face_area;
            const flow_float nz = nsign * sz[ip] / face_area;

            // フラックスは物理の厳密 FVS。前処理は時間項のみ (保存形 (Γ_c V/Δτ'+A_c)ΔQ=-R)。
            flow_float a_plus[5][5];
            flow_float k_off[5][5];
            block_dplur::build_jacobian_split(gamma, nx, ny, nz, vx, vy, vz,
                                              local_enthalpy, local_sonic, a_plus, k_off);
            block_dplur::add_scaled_5x5(D0, a_plus, face_area);

            const flow_float dcc_x = ccx[other_ic] - ccx[ic];
            const flow_float dcc_y = ccy[other_ic] - ccy[ic];
            const flow_float dcc_z = ccz[other_ic] - ccz[ic];
            const flow_float dcc = max(sqrt(dcc_x*dcc_x + dcc_y*dcc_y + dcc_z*dcc_z), static_cast<flow_float>(1.0e-30));
            const flow_float dcc_dot_s = max(fabs(dcc_x*sx[ip] + dcc_y*sy[ip] + dcc_z*sz[ip]), static_cast<flow_float>(1.0e-30));
            const flow_float delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<flow_float>(1.0e-30));
            // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる
            // (旧 face_area·(2ν/delta) は ≈2ν に潰れ近軸で r 重み喪失・ゼロ面積面にスプリアス。詳細は site1 コメント)。
            const flow_float viscous_diag = static_cast<flow_float>(2.0) * nu_eff * delta / dcc;
            block_dplur::add_identity_scaled(D0, viscous_diag);

            // 近傍 += k_off S ΔQ_nbr (= -A_c⁻ S ΔQ_nbr)。
            if (other_ic < nCells) {
                flow_float dqn[5];
                block_dplur::load_block_vec(other_ic, dq_old_0, dq_old_1, dq_old_2, dq_old_3, dq_old_4, dqn);
                #pragma unroll
                for (int i = 0; i < 5; ++i) dqn[i] *= face_area;
                block_dplur::multiply_add_5x5_vec(k_off, dqn, nbr);
            }
        }

        #pragma unroll
        for (int i = 0; i < 5; ++i) b[i] += nbr[i];

        // 軸対称ソースヤコビアン (物理ブロック D0 へ float で加算、既存 block と同式)。
        // axisRFloor 帯 (r 床, ソース不課) は Jacobian も課さない。
        if (isAxisymmetric == 1 &&
            !(axisRFloor > (flow_float)0.0 && ccy[ic] < axisRFloor)) {
            const flow_float A_pl = A_planar[ic];
            const flow_float r_eff = max(v / max(A_pl, static_cast<flow_float>(1.0e-30)), static_cast<flow_float>(1.0e-30));
            const flow_float g1 = gamma - static_cast<flow_float>(1.0);
            const flow_float q2 = vx*vx + vy*vy + vz*vz;
            const flow_float mu_total = laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0));
            const flow_float hoop = static_cast<flow_float>(2.0) * mu_total / (density * r_eff);
            // ∂P/∂Q の第 1 成分は一般 EOS で χ_eos + κ e_k (χ_eos = c² − κ h)。CPG では χ_eos=0 で
            // 従来式 ½κq² にビット一致。TP (thermalMethod 2) では χ_eos≠0 で、これを落とすと
            // 対流 Jacobian (rvec, 下記) と軸ソース Jacobian が不整合になり、ホップ項が支配する
            // 軸近傍で block-DPLUR が発散する (case/42 run_0020: 一定 cp 種・陽解法では完走、
            // 実 NASA-9 + 陰解法のみ発散 → 2026-08-17 に特定)。
            const flow_float chi_hoop = local_sonic*local_sonic
                                      - g1*(local_enthalpy - static_cast<flow_float>(0.5)*q2);
            D0[2][0] += -A_pl * (chi_hoop + static_cast<flow_float>(0.5)*g1*q2 + hoop*vy);
            D0[2][1] +=  A_pl * (g1*vx);
            D0[2][2] +=  A_pl * (g1*vy + hoop);
            D0[2][3] +=  A_pl * (g1*vz);
            D0[2][4] += -A_pl * g1;
        }

        // Γ_c のランク1寄与: g=(1,u,v,w,H_t), r=∂p/∂Q=(χ_eos+κek,-κu,-κv,-κw,κ), γ=(V/Δτ')·(1-β)/(βc²)。
        // CPG では H_t=c²/(γ-1)+ek・χ_eos=0 で従来式に簡約 (ビット不変)。TP は実 H_t(=local_enthalpy) と
        // χ_eos=c²−κh を使う (build_jacobian_split と同じ一般EOS整合)。κ=γ-1。
        const flow_float ek = static_cast<flow_float>(0.5) * velMag * velMag;
        const flow_float gm1 = gamma - static_cast<flow_float>(1.0);
        const flow_float Htot = local_enthalpy;   // 実 Ht[ic] (CPG/TP 統一。CPG も Ht[ic]=ek+c²/(γ-1))
        const flow_float chi_eos = local_sonic*local_sonic - gm1*(local_enthalpy - ek);  // c²−κh (CPG で ≈0)
        flow_float gvec[5] = { static_cast<flow_float>(1.0), vx, vy, vz, Htot };
        const flow_float rvec[5] = { chi_eos + gm1*ek, -gm1*vx, -gm1*vy, -gm1*vz, gm1 };
        const double dbeta = static_cast<double>(beta);
        const double alpha = (1.0 - dbeta) / (dbeta * static_cast<double>(local_sonic) * static_cast<double>(local_sonic));
        const double gam = static_cast<double>(v_over_dtau) * alpha;   // = (V/Δτ')·α

        // D0 を float で 2 RHS 同時に解く: y=D0⁻¹b, z=D0⁻¹g。
        flow_float y[5], z[5];
        const bool ok = block_dplur::solve_5x5_2rhs(D0, b, gvec, y, z);

        // Sherman-Morrison スカラー (悪条件 1/β はここだけ double): x = y - [γ(rᵀy)/(1+γ(rᵀz))] z。
        double ry = 0.0, rz = 0.0;
        #pragma unroll
        for (int i = 0; i < 5; ++i) { ry += static_cast<double>(rvec[i]) * y[i]; rz += static_cast<double>(rvec[i]) * z[i]; }
        const double denom = 1.0 + gam * rz;
        const double sfac = (fabs(denom) > 1.0e-300) ? gam * ry / denom : 0.0;

        flow_float correction[5];
        #pragma unroll
        for (int i = 0; i < 5; ++i)
            correction[i] = ok ? static_cast<flow_float>((static_cast<double>(y[i]) - sfac * static_cast<double>(z[i]))
                                                         * static_cast<double>(implicit_relax))
                               : static_cast<flow_float>(0.0);

        block_dplur::store_block_vec(ic, correction, dq_new_0, dq_new_1, dq_new_2, dq_new_3, dq_new_4);
    }
}

// block DPLUR の sweep 間バッファ入れ替え。ドライバ側から各 sweep 後に明示的に呼ぶ
// （旧実装は wrapper 内部で暗黙に swap していたが、古典 DPLUR では制御フローを明示化する）。
// 近傍 dq の AoS バッファ (stride 8)。wrapper で nCells_all に合わせて確保し、swap で old/new を入れ替える。
static flow_float* g_dqPackOld = nullptr;
static flow_float* g_dqPackNew = nullptr;
static geom_int    g_dqPackN   = 0;

void swapBlockImplicitCorrectionBuffers(variables& var)
{
    std::swap(g_dqPackOld, g_dqPackNew);
    std::swap(var.c_d["dq_block_old_0"], var.c_d["dq_block_new_0"]);
    std::swap(var.c_d["dq_block_old_1"], var.c_d["dq_block_new_1"]);
    std::swap(var.c_d["dq_block_old_2"], var.c_d["dq_block_new_2"]);
    std::swap(var.c_d["dq_block_old_3"], var.c_d["dq_block_new_3"]);
    std::swap(var.c_d["dq_block_old_4"], var.c_d["dq_block_new_4"]);
}

// scalar 対角版 (blockDPLUR==0) の sweep 間バッファ入れ替え。block 版と同様にドライバ側から呼ぶ。
void swapScalarImplicitCorrectionBuffers(variables& var)
{
    std::swap(var.c_d["dq_ro_old"],   var.c_d["dq_ro_new"]);
    std::swap(var.c_d["dq_roUx_old"], var.c_d["dq_roUx_new"]);
    std::swap(var.c_d["dq_roUy_old"], var.c_d["dq_roUy_new"]);
    std::swap(var.c_d["dq_roUz_old"], var.c_d["dq_roUz_new"]);
    std::swap(var.c_d["dq_roe_old"],  var.c_d["dq_roe_new"]);
}

//TODO __global__ void runge_kutta_dual_explicit_d
//TODO // see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
//TODO // N: previous outer step , M: previous inner loop
//TODO ( 
//TODO  geom_int dt ,
//TODO 
//TODO  // mesh structure
//TODO  geom_int nCells_all , geom_int nCells,
//TODO  geom_float* vol ,
//TODO 
//TODO  // variables
//TODO  flow_float* ro  ,
//TODO  flow_float* roUx  ,
//TODO  flow_float* roUy  ,
//TODO  flow_float* roUz  ,
//TODO  flow_float* roe  ,
//TODO 
//TODO  flow_float* roN ,
//TODO  flow_float* roUxN ,
//TODO  flow_float* roUyN ,
//TODO  flow_float* roUzN ,
//TODO  flow_float* roeN ,
//TODO 
//TODO  flow_float* roM ,
//TODO  flow_float* roUxM ,
//TODO  flow_float* roUyM ,
//TODO  flow_float* roUzM ,
//TODO  flow_float* roeM ,
//TODO 
//TODO  flow_float* res_ro,
//TODO  flow_float* res_roUx,
//TODO  flow_float* res_roUy,
//TODO  flow_float* res_roUz,
//TODO  flow_float* res_roe,
//TODO 
//TODO  flow_float* res_ro_dual ,
//TODO  flow_float* res_roUx_dual ,
//TODO  flow_float* res_roUy_dual ,
//TODO  flow_float* res_roUz_dual ,
//TODO  flow_float* res_roe_dual 
//TODO 
//TODO )
//TODO {
//TODO     geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
//TODO 
//TODO     geom_float v = vol[ic];
//TODO 
//TODO     if (ic < nCells) {
//TODO         // N: previous outer step , M: previous inner loop
//TODO         res_ro_dual[ic]   = -(ro[ic]-roN[ic])*v/dt     + res_ro[ic];
//TODO         res_roUx_dual[ic] = -(roUx[ic]-roUxN[ic])*v/dt + res_roUx[ic];
//TODO         res_roUy_dual[ic] = -(roUy[ic]-roUyN[ic])*v/dt + res_roUy[ic];
//TODO         res_roUz_dual[ic] = -(roUz[ic]-roUzN[ic])*v/dt + res_roUz[ic];
//TODO         res_roe_dual[ic]  = -(roe[ic]-roeN[ic])*v/dt   + res_roe[ic];
//TODO     }
//TODO     __syncthreads();
//TODO }

void timeIntegration_d_wrapper(int loop , solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , int lineStoreK)
{
    // 軸対称エンコード: 0=非軸対称 / 1=r 重み方式 (hoop Jacobian) / 2=SU2 流 planar+ソース (SU2 4x4 Jacobian)。
    // ==1 判定しかしない旧カーネル (scalar/lowmach) は 2 のとき軸対称 Jacobian を持たない (source lag, 定常解不変)。
    // 診断 (env): FORGE_DIAG_SU2JAC_OFF=1 で SU2 ソース Jacobian を落とす (ソースは残す = 完全 lag)。
    static const bool diagSu2JacOff = (getenv("FORGE_DIAG_SU2JAC_OFF") != nullptr);
    int axisymEnc = (cfg.isAxisymmetric == 1) ? ((cfg.axisymMethod == 1) ? 2 : 1) : 0;
    if (diagSu2JacOff && axisymEnc == 2) axisymEnc = 0;
    // 診断 near-axis 安定化係数を env から 1 度だけ device へ設定 (既定 0 = 不変)。
    static bool s_axisAlphaInit = false;
    if (!s_axisAlphaInit) {
        float a = 0.0f;
        if (const char* e = getenv("FORGE_AXIS_DIAG_ALPHA")) a = static_cast<float>(atof(e));
        cudaMemcpyToSymbol(g_axisDiagAlpha, &a, sizeof(float));
        s_axisAlphaInit = true;
    }
    if (cfg.timeIntegration == 4) { // 4th order runge kutta
        runge_kutta_exp_4th_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            loop, 
            cfg.coef_DT_4thRunge[loop],
            cfg.coef_Res_4thRunge[loop],
            cfg.dt ,
            var.c_d["dt_local"],

            // mesh structure
            msh.nCells_all , msh.nCells ,
            var.c_d["volume"],

            // basic variables
            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] ,
            var.c_d["res_ro_m"], var.c_d["res_roUx_m"], var.c_d["res_roUy_m"], var.c_d["res_roUz_m"] , var.c_d["res_roe_m"] 
        ) ;

    } else if (cfg.timeIntegration == 1 or cfg.timeIntegration == 3) { // explicit
        runge_kutta_exp_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            loop,
            cfg.coef_N[loop],
            cfg.coef_M[loop],
            cfg.coef_Res[loop],
            cfg.dt , 
            var.c_d["dt_local"],

            // mesh structure
            msh.nCells_all , msh.nCells ,
            var.c_d["volume"],

            // basic variables
            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] 
        ) ;
    } else if (cfg.timeIntegration == 11) { // implicit defect-correction with diagonal Jacobian approximation
        if (cfg.blockDPLUR == 1) {
            // レジスタ過多のため専用の小さい block サイズで起動（__launch_bounds__ と整合）。
            const int block_threads = BLOCK_DPLUR_THREADS;
            const int block_grid = (msh.nCells_all + block_threads - 1) / block_threads;
            if (cfg.lowMachPrecond >= 2) {
              // Phase 4: 完全 Γ⁻¹A 前処理の倍精度カーネル (dt_local は前処理 Δτ' に拡大済)。
              // lowMachPrecond==2: RHS 散逸 c' (Phase 1) + 本 LHS 前処理。
              // lowMachPrecond==3: LHS 前処理のみ (RHS 散逸は c_hat で =0 と不変)。本カーネルは
              //   擬似時間項 Γ_c のみ前処理しフラックス A_c は非前処理ゆえ、収束解は前処理なしと一致する。
              implicit_defect_correction_block_precond_d<<<block_grid , block_threads>>>(
                loop, cfg.dt, var.c_d["dt_local"], cfg.implicitRelax, var.c_d["gamma"], cfg.precondEps,
                msh.nCells_all, msh.nCells, var.c_d["volume"],
                msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
                var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
                cfg.visc, var.c_d["vis_turb"], var.c_d["sonic"],
                var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["Ht"],
                var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"],
                var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
                var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
                axisymEnc,
                (cfg.isAxisymmetric == 1) ? ((cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1) ? var.c_d["A_closure_y"] : var.c_d["A_planar"]) : var.c_d["volume"],
                cfg.axisRFloor,
                cfg.unsteadyDiagCoef
              );
            } else {
            // implicitSolvePrecision: 0=float (既定・高速), 1=double (軸対称近軸の根治, 遅い)。
            // 同じテンプレートカーネルを ST=float/double で起動。引数は共通 (FORGE_BDPLUR_ARGS)。
            #define FORGE_BDPLUR_ARGS \
                loop, cfg.dt, var.c_d["dt_local"], cfg.implicitRelax, var.c_d["gamma"], \
                (cfg.thermalMethod == 2 ? 1 : 0), \
                msh.nCells_all, msh.nCells, var.c_d["volume"], \
                msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d, \
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"], \
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"], \
                var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"], \
                cfg.visc, var.c_d["vis_turb"], var.c_d["sonic"], \
                var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["Ht"], \
                var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"], \
                var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"], \
                var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"], \
                var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"], \
                var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"], \
                var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"], \
                var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"], \
                var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"], \
                var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"], \
                axisymEnc, (cfg.isAxisymmetric == 1) ? ((cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1) ? var.c_d["A_closure_y"] : var.c_d["A_planar"]) : var.c_d["volume"], cfg.axisRFloor, cfg.unsteadyDiagCoef, \
                nullptr,  /* axis_flag: 旧 nodeAxisDirichlet の全 5 行 decouple (撤去) */ \
                ((cfg.discretization == "node" && cfg.isAxisymmetric == 1) ? msh.axis_flag_d : nullptr),  /* axis_ur_flag: 軸ノードの roUy 行 decouple (常時) */ \
                ((cfg.discretization == "node" && cfg.isAxisymmetric == 1) ? msh.axis_flag_d : nullptr),  /* axis_flag_src: SU2 流 (enc==2) の軸ソース Jacobian ガード (軸ノードはソース 0) */ \
                ((cfg.discretization == "node" && cfg.nodeWallDirichlet == 1) ? msh.wall_flag_d : nullptr),  /* wall_flag: 壁運動量3行 decouple */ \
                ((cfg.discretization == "node" && cfg.nodeWallDirichlet == 1 && cfg.nodeIsothermalEnergyBC != 1) ? msh.iso_wall_flag_d : nullptr),  /* iso_wall_flag: 等温壁 roe 行 decouple (T ピンと対)。弱形式 (nodeIsothermalEnergyBC=1) では単位行化も rowDec も外す (plan boundary-weak-isothermal-wall §4.3) */ \
                (weakIsoWall::active(cfg, msh) ? weakIsoWall::diagBuf(msh) : nullptr),  /* weakIsoDiag: 弱形式の近似対角 */ \
                (flow_float)(cfg.cp / max(cfg.gamma, 1.0e-30)),  /* weakIsoCv = cp/gamma = c_v (CPG) */ \
                ((cfg.discretization == "node") ? 1 : 0),  /* isNode: 5e 境界半割面の粘性対角スキップ */ \
                ((cfg.lineImplicit == 1) ? msh.line_prev_d : nullptr), \
                ((cfg.lineImplicit == 1) ? msh.line_next_d : nullptr), \
                msh.line_Kprev_d, msh.line_Knext_d, (((loop == 0) && (lineStoreK != 0)) ? 1 : 0), cfg.lineViscCoupling,  /* line-implicit */ \
                ((cfg.implicitSolvePrecision == 0 && cfg.lineImplicit == 0 && cfg.blockDPLURDiagCache != 0) ? 1 : 0),  /* useDiagCache: float・point 経路のみ */ \
                (usePack ? (const flow_float*)g_dqPackOld : nullptr), (usePack ? g_dqPackNew : nullptr),  /* 近傍 dq の AoS 版 */ \
                cfg.implicitThermalJacobian,  /* エネルギー行の熱伝導 Jacobian / 等温壁の拘束の行 (ビットマスク) */ \
                ((cfg.implicitThermalJacobian & 1) ? var.c_d["thermCond"] : nullptr), \
                ((cfg.implicitThermalJacobian & 1) ? var.c_d["cp"] : nullptr), \
                ((cfg.implicitThermalJacobian & 1) ? var.p_d["fx"] : nullptr), \
                cfg.turbulentPrandtl
            // 近傍 dq の AoS 経路: line-implicit と node 周期 (SoA だけを直接書き換える) では使わない。
            const bool usePack = (cfg.lineImplicit == 0) && (cfg.blockDPLURDqPack != 0) &&
                                 !(cfg.discretization == "node" && msh.periodicRoot_d != nullptr && msh.nPeriodicMembers > 0);
            if (usePack && (g_dqPackOld == nullptr || g_dqPackN != msh.nCells_all)) {
                if (g_dqPackOld) { cudaFree(g_dqPackOld); cudaFree(g_dqPackNew); }
                const size_t nb = (size_t)msh.nCells_all * 8 * sizeof(flow_float);
                gpuErrchk(cudaMalloc((void**)&g_dqPackOld, nb)); gpuErrchk(cudaMalloc((void**)&g_dqPackNew, nb));
                gpuErrchk(cudaMemset(g_dqPackOld, 0, nb)); gpuErrchk(cudaMemset(g_dqPackNew, 0, nb));
                g_dqPackN = msh.nCells_all;
            }
            if (cfg.implicitSolvePrecision == 1)
                implicit_defect_correction_block_d<double><<<block_grid , block_threads>>>(FORGE_BDPLUR_ARGS);
            else
                implicit_defect_correction_block_d<float><<<block_grid , block_threads>>>(FORGE_BDPLUR_ARGS);
            #undef FORGE_BDPLUR_ARGS
            }
        } else {
            implicit_defect_correction_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                loop,
                cfg.dt,
                var.c_d["dt_local"],
                cfg.implicitRelax,
                var.c_d["gamma"],
                msh.nCells_all,
                msh.nCells,
                var.c_d["volume"],
                msh.map_plane_cells_d,
                msh.map_cell_planes_index_d,
                msh.map_cell_planes_d,
                var.c_d["ccx"],
                var.c_d["ccy"],
                var.c_d["ccz"],
                var.p_d["sx"],
                var.p_d["sy"],
                var.p_d["sz"],
                var.p_d["ss"],
                var.c_d["ro"],
                var.c_d["roUx"],
                var.c_d["roUy"],
                var.c_d["roUz"],
                var.c_d["roe"],
                var.c_d["roN"],
                var.c_d["roUxN"],
                var.c_d["roUyN"],
                var.c_d["roUzN"],
                var.c_d["roeN"],
                cfg.visc,
                var.c_d["vis_turb"],
                var.c_d["sonic"],
                var.c_d["Ux"],
                var.c_d["Uy"],
                var.c_d["Uz"],
                var.c_d["res_ro"],
                var.c_d["res_roUx"],
                var.c_d["res_roUy"],
                var.c_d["res_roUz"],
                var.c_d["res_roe"],
                var.c_d["dq_ro_old"],
                var.c_d["dq_roUx_old"],
                var.c_d["dq_roUy_old"],
                var.c_d["dq_roUz_old"],
                var.c_d["dq_roe_old"],
                var.c_d["dq_ro_new"],
                var.c_d["dq_roUx_new"],
                var.c_d["dq_roUy_new"],
                var.c_d["dq_roUz_new"],
                var.c_d["dq_roe_new"],
                axisymEnc,
                (cfg.isAxisymmetric == 1) ? ((cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1) ? var.c_d["A_closure_y"] : var.c_d["A_planar"]) : var.c_d["volume"],
                cfg.axisRFloor,
                cfg.unsteadyDiagCoef
            );
        }
        // 古典 DPLUR: buffer swap と Q への commit はドライバ側 (main.cpp blockDPLURSolve /
        // applyBlockImplicitCorrection) で明示的に行う。ここでは sweep カーネルの起動のみ。
//TODO    } else if (cfg.timeIntegration == 10) { // implicit (m-time stepping & explicit scheme)
//TODO        runge_kutta_dual_explicit_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
//TODO            cfg.dt , 
//TODO
//TODO            // mesh structure
//TODO            msh.nCells_all , msh.nCells ,
//TODO            var.c_d["volume"],
//TODO
//TODO            // basic variables
//TODO            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
//TODO            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
//TODO            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
//TODO            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] ,
//TODO            var.c_d["res_ro_m"], var.c_d["res_roUx_m"], var.c_d["res_roUy_m"], var.c_d["res_roUz_m"] , var.c_d["res_roe_m"] 
//TODO        ) ;
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}
// =============================================================================
// line-implicit: ライン block-Thomas (plans/active/time_integration-line-implicit.md)。
// sweep カーネルが保存した diag (loop0 凍結)・K (loop0)・rhs (毎 sweep, ライン外 lag 込み) から、
// 各ラインの block 三重対角系
//   D_k ΔQ_k − Kprev_k ΔQ_{k-1} − Knext_k ΔQ_{k+1} = rhs_k
// を前進消去+後退代入で厳密に解き dq_new を上書きする。1 ライン = 1 スレッド (v1)、内部 double。
// scratch: W_k = D̃_k⁻¹ Knext_k (25/cell), y_k = D̃_k⁻¹ b̃_k (5/cell)。
namespace line_implicit {

__device__ __forceinline__ bool lu5_factor(double A[5][5], int piv[5])
{
    for (int col = 0; col < 5; ++col) {
        int pv = col; double pa = fabs(A[col][col]);
        for (int r = col + 1; r < 5; ++r) { const double c = fabs(A[r][col]); if (c > pa) { pv = r; pa = c; } }
        if (pa < 1.0e-30) return false;
        piv[col] = pv;
        if (pv != col) for (int k = 0; k < 5; ++k) { const double tmp = A[col][k]; A[col][k] = A[pv][k]; A[pv][k] = tmp; }
        const double inv = 1.0 / A[col][col];
        for (int r = col + 1; r < 5; ++r) {
            const double f = A[r][col] * inv;
            A[r][col] = f;                      // L を下三角に格納
            for (int k = col + 1; k < 5; ++k) A[r][k] -= f * A[col][k];
        }
    }
    return true;
}

__device__ __forceinline__ void lu5_solve(const double A[5][5], const int piv[5], double x[5])
{
    // LAPACK getrs 流: ① 行交換を全て先に適用 (LASWP) ② 単位下三角 L 前進代入 ③ U 後退代入。
    // 交換と代入をインタリーブする書き方は、後段ピボットが L 部分も行交換する getrf 形格納と
    // 非整合で誤解を返す (2026-09-02 に numpy 照合で確認済みの罠)。
    for (int col = 0; col < 5; ++col) {
        if (piv[col] != col) { const double tmp = x[col]; x[col] = x[piv[col]]; x[piv[col]] = tmp; }
    }
    for (int col = 0; col < 5; ++col) {
        for (int r = col + 1; r < 5; ++r) x[r] -= A[r][col] * x[col];
    }
    for (int r = 4; r >= 0; --r) {
        double s = x[r];
        for (int c = r + 1; c < 5; ++c) s -= A[r][c] * x[c];
        x[r] = s / A[r][r];
    }
}

} // namespace line_implicit

__global__ void lineThomas_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev, const flow_float* Knext,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 const flow_float* dq_old_0, const flow_float* dq_old_1, const flow_float* dq_old_2, const flow_float* dq_old_3, const flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax,
 double* Wd, double* yd
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    bool fail = false;

    // ---- 前進消去 ----
    for (geom_int p = b; p < e && !fail; ++p) {
        const geom_int ic = line_cells[p];
        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            // M -= Kprev·W_{k-1},  b̃_k = b_k + Kprev·y_{k-1}
            // (標準形 L=−Kprev, U=−Knext につき符号は加算側に出る)
            double Kp[5][5];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j)
                    Kp[i][j] = (double)Kprev[(size_t)ic * 25 + i * 5 + j];
            for (int i = 0; i < 5; ++i) {
                double bacc = 0.0;
                for (int m = 0; m < 5; ++m) bacc += Kp[i][m] * yd[(size_t)icm * 5 + m];
                bk[i] += bacc;
                for (int j = 0; j < 5; ++j) {
                    double macc = 0.0;
                    for (int m = 0; m < 5; ++m) macc += Kp[i][m] * Wd[(size_t)icm * 25 + m * 5 + j];
                    M[i][j] -= macc;
                }
            }
        }
        int piv[5];
        if (!line_implicit::lu5_factor(M, piv)) { fail = true; break; }
        line_implicit::lu5_solve(M, piv, bk);                  // y_k
        for (int i = 0; i < 5; ++i) yd[(size_t)ic * 5 + i] = bk[i];
        if (p + 1 < e) {                                        // W_k = M⁻¹·Knext_k
            for (int j = 0; j < 5; ++j) {
                double col[5];
                for (int i = 0; i < 5; ++i) col[i] = (double)Knext[(size_t)ic * 25 + i * 5 + j];
                line_implicit::lu5_solve(M, piv, col);
                for (int i = 0; i < 5; ++i) Wd[(size_t)ic * 25 + i * 5 + j] = col[i];
            }
        }
    }

    // ---- 後退代入 (relax を掛けて dq_new へ) / 失敗時は前回反復値を保持 ----
    if (fail) {
        for (geom_int p = b; p < e; ++p) {
            const geom_int ic = line_cells[p];
            dq_new_0[ic] = dq_old_0[ic]; dq_new_1[ic] = dq_old_1[ic]; dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic]; dq_new_4[ic] = dq_old_4[ic];
        }
        return;
    }
    double dq[5];
    {
        const geom_int ic = line_cells[e - 1];
        for (int i = 0; i < 5; ++i) dq[i] = yd[(size_t)ic * 5 + i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
    }
    for (geom_int p = e - 2; p >= b; --p) {
        const geom_int ic = line_cells[p];
        double nx[5];
        for (int i = 0; i < 5; ++i) {
            double acc = yd[(size_t)ic * 5 + i];
            for (int m = 0; m < 5; ++m) acc += Wd[(size_t)ic * 25 + i * 5 + m] * dq[m];
            nx[i] = acc;
        }
        for (int i = 0; i < 5; ++i) dq[i] = nx[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
        if (p == b) break;   // geom_int が unsigned の場合の負回りガード
    }
}

// v2 (plans/active/time_integration-line-implicit-viscous-v2.md): factor/solve 分離。
// 前進消去の M̃_k = D_k − Kprev·W_{k−1} の構築・LU 分解・W_k = M̃⁻¹Knext は rhs に依存しない
// (D, K は storeLU 時に凍結) ので、storeLU のタイミングで 1 回だけ行い LU/piv/W を保存する。
// sweep 毎の solve は保存済み因子での代入 (Kp·y 25 積 + LASWP 前進/後退) だけになる。
// モノリシック版 (lineThomas_d) は毎 sweep この 5 列 solve + Kp·W (625 積) を再計算しており、
// これが DDES A/B での step 単価 2.44 倍の主犯 — 分離は厳密 (近似ゼロ) の最適化。
__global__ void lineThomasFactor_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev, const flow_float* Knext,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 double* Wd, double* LUd, signed char* pivd, unsigned char* faild
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    faild[l] = 0;

    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            double Kp[5][5];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j)
                    Kp[i][j] = (double)Kprev[(size_t)ic * 25 + i * 5 + j];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j) {
                    double macc = 0.0;
                    for (int m = 0; m < 5; ++m) macc += Kp[i][m] * Wd[(size_t)icm * 25 + m * 5 + j];
                    M[i][j] -= macc;
                }
        }
        int piv[5];
        if (!line_implicit::lu5_factor(M, piv)) { faild[l] = 1; return; }
        for (int i = 0; i < 5; ++i) {
            pivd[(size_t)ic * 5 + i] = (signed char)piv[i];
            for (int j = 0; j < 5; ++j) LUd[(size_t)ic * 25 + i * 5 + j] = M[i][j];
        }
        if (p + 1 < e) {                                        // W_k = M̃⁻¹·Knext_k
            for (int j = 0; j < 5; ++j) {
                double col[5];
                for (int i = 0; i < 5; ++i) col[i] = (double)Knext[(size_t)ic * 25 + i * 5 + j];
                line_implicit::lu5_solve(M, piv, col);
                for (int i = 0; i < 5; ++i) Wd[(size_t)ic * 25 + i * 5 + j] = col[i];
            }
        }
    }
}

__global__ void lineThomasSolve_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev,
 const double* Wd, const double* LUd, const signed char* pivd, const unsigned char* faild,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 const flow_float* dq_old_0, const flow_float* dq_old_1, const flow_float* dq_old_2, const flow_float* dq_old_3, const flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax,
 double* yd
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    if (faild[l] != 0) {
        for (geom_int p = b; p < e; ++p) {
            const geom_int ic = line_cells[p];
            dq_new_0[ic] = dq_old_0[ic]; dq_new_1[ic] = dq_old_1[ic]; dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic]; dq_new_4[ic] = dq_old_4[ic];
        }
        return;
    }

    // ---- 前進 (保存因子で代入のみ) ----
    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            for (int i = 0; i < 5; ++i) {
                double bacc = 0.0;
                for (int m = 0; m < 5; ++m)
                    bacc += (double)Kprev[(size_t)ic * 25 + i * 5 + m] * yd[(size_t)icm * 5 + m];
                bk[i] += bacc;
            }
        }
        double M[5][5];
        int piv[5];
        for (int i = 0; i < 5; ++i) {
            piv[i] = (int)pivd[(size_t)ic * 5 + i];
            for (int j = 0; j < 5; ++j) M[i][j] = LUd[(size_t)ic * 25 + i * 5 + j];
        }
        line_implicit::lu5_solve(M, piv, bk);
        for (int i = 0; i < 5; ++i) yd[(size_t)ic * 5 + i] = bk[i];
    }

    // ---- 後退代入 (relax を掛けて dq_new へ) ----
    double dq[5];
    {
        const geom_int ic = line_cells[e - 1];
        for (int i = 0; i < 5; ++i) dq[i] = yd[(size_t)ic * 5 + i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
    }
    for (geom_int p = e - 2; p >= b; --p) {
        const geom_int ic = line_cells[p];
        double nx[5];
        for (int i = 0; i < 5; ++i) {
            double acc = yd[(size_t)ic * 5 + i];
            for (int m = 0; m < 5; ++m) acc += Wd[(size_t)ic * 25 + i * 5 + m] * dq[m];
            nx[i] = acc;
        }
        for (int i = 0; i < 5; ++i) dq[i] = nx[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
        if (p == b) break;   // geom_int が unsigned の場合の負回りガード
    }
}

// 診断: ライン CV を「保存済み diag/rhs の点解」だけで更新する (K/Thomas 不使用)。
// FORGE_LINE_DEBUG_POINT=1 で有効。格納 (diag/rhs) の正しさと Thomas 本体の切り分け用。
__global__ void lineDebugPoint_d
(
 geom_int nLines, const geom_int* line_offsets, const geom_int* line_cells,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    for (geom_int p = line_offsets[l]; p < line_offsets[l + 1]; ++p) {
        const geom_int ic = line_cells[p];

        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        int piv[5];
        if (line_implicit::lu5_factor(M, piv)) {
            line_implicit::lu5_solve(M, piv, bk);
        } else {
            for (int i = 0; i < 5; ++i) bk[i] = 0.0;
        }
        dq_new_0[ic] = (flow_float)(implicit_relax * bk[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * bk[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * bk[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * bk[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * bk[4]);
    }
}

// v2: モノリシック版へ戻す退避スイッチ (FORGE_LINE_MONO=1)。既定は factor/solve 分離。
static bool lineMonoEnabled() {
    static const bool v = [](){ const char* e = getenv("FORGE_LINE_MONO"); return e && atoi(e) != 0; }();
    return v;
}

// factor 位相: storeLU した sweep の直後に 1 回だけ呼ぶ (blockDPLURSolve が管理)。
void lineThomasFactor_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (msh.nImplicitLines <= 0) return;
    if (lineMonoEnabled()) return;   // モノリシック時は毎 sweep の lineThomas_d が全てやる
    static const bool dbgPoint = [](){ const char* e = getenv("FORGE_LINE_DEBUG_POINT"); return e && atoi(e) != 0; }();
    static const bool dbgNoop = [](){ const char* e = getenv("FORGE_LINE_NOOP"); return e && atoi(e) != 0; }();
    if (dbgNoop || dbgPoint) return;
    const int threads = 64;
    const int grid = (int)((msh.nImplicitLines + threads - 1) / threads);
    lineThomasFactor_d<<<grid, threads>>>(
        msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d,
        msh.line_Kprev_d, msh.line_Knext_d,
        var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"],
        var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"],
        var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"],
        var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"],
        var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"],
        msh.line_W_d, msh.line_LU_d, msh.line_piv_d, msh.line_fail_d);
    gpuErrchk( cudaPeekAtLastError() );
}

void lineThomas_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (msh.nImplicitLines <= 0) return;
    static const bool dbgPoint = [](){ const char* e = getenv("FORGE_LINE_DEBUG_POINT"); return e && atoi(e) != 0; }();
    static const bool dbgNoop = [](){ const char* e = getenv("FORGE_LINE_NOOP"); return e && atoi(e) != 0; }();
    if (dbgNoop) return;   // 切り分け: ライン CV は dq 据え置き (sweep 内 placeholder のまま)
    if (dbgPoint) {
        const int threads0 = 64;
        const int grid0 = (int)((msh.nImplicitLines + threads0 - 1) / threads0);
        lineDebugPoint_d<<<grid0, threads0>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d,
            var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"],
            var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"],
            var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"],
            var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"],
            var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"],
            var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
            var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
            cfg.implicitRelax);
        gpuErrchk( cudaPeekAtLastError() );
        return;
    }
    const int threads = 64;
    const int grid = (int)((msh.nImplicitLines + threads - 1) / threads);
    if (!lineMonoEnabled()) {
        // v2 既定: 保存済み LU/piv/W での代入のみ (factor は lineThomasFactor_d_wrapper が実施済み)。
        lineThomasSolve_d<<<grid, threads>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d,
            msh.line_Kprev_d,
            msh.line_W_d, msh.line_LU_d, msh.line_piv_d, msh.line_fail_d,
            var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
            var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
            var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
            cfg.implicitRelax,
            msh.line_y_d);
        gpuErrchk( cudaPeekAtLastError() );
        return;
    }
    lineThomas_d<<<grid, threads>>>(
        msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d,
        msh.line_Kprev_d, msh.line_Knext_d,
        var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"],
        var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"],
        var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"],
        var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"],
        var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"],
        var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
        var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
        var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
        cfg.implicitRelax,
        msh.line_W_d, msh.line_y_d);
    gpuErrchk( cudaPeekAtLastError() );
}
```

## 参考: `methods/time_integration/implementation.md`

```
# 時間積分 — 実装

forge の時間積分・更新カーネルの実装とソース対応をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/cuda_forge/timeIntegration_d.cuh`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cuh) | カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/timeIntegration_d.cu`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) | RK / 陰解法カーネル本体 (約 950 行) |
| [`solver_density_cuda/cuda_forge/update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) | 外側・内側ループ前後の Q コピー、陰解法補正反映 |
| [`solver_density_cuda/cuda_forge/implicitCorrection_d.cu`](../../solver_density_cuda/cuda_forge/implicitCorrection_d.cu) | dual-time 陽 (簡易) スキームの補正 |
| [`solver_density_cuda/cuda_forge/setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu) | 局所 $\Delta t_{\text{loc}}$ 計算 |
| [`solver_density_cuda/update.cpp`](../../solver_density_cuda/update.cpp) | 旧 CPU 経路 (現状未使用) |

## エントリポイント

[`timeIntegration_d_wrapper`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L780)
が共通入口。`cfg.timeIntegration` の値に応じて次を呼ぶ。

| 値 | 呼び出し |
| --- | --- |
| `4` | `runge_kutta_exp_4th_d` |
| `1`, `3` | `runge_kutta_exp_d` |
| `11`, `blockDPLUR == 1` | `implicit_defect_correction_block_d`（5×5 block DPLUR、LU-SGS $A^\pm$。**推奨既定**） |
| `11`, `blockDPLUR == 0` | `implicit_defect_correction_d`（scalar 対角版＝スペクトル半径。軽量だが擬似 CFL が低く収束が遅い） |

> `solverConfig::initTimeIntegrationScheme` の `case 11` は `blockDPLUR ∈ {0,1}` を受理（それ以外を throw）。両者とも古典 DPLUR 制御フロー（残差固定 + `nStepInner` sweep + 単一 commit）に対応。`unsteady == 1`（dual-time）は dispatcher 側で throw する（本体未実装）。

> **軸対称ソース項の Jacobian**: scalar 版 `implicit_defect_correction_d` も block と整合させ、軸対称フープ源 `res_roUy += (P − τ_θθ)·A_planar` の Jacobian 対角成分 `A_pl·((γ−1)u_y + 2μ/(ρ r_eff))`（非負側、per-cell γ で TP 整合）を roUy 方程式の対角に陰化する。ただし scalar は式間連成 (源 Jacobian の非対角) を表現できないため、軸近傍以外が律速の強膨張ケース（例 `case/29` 出口コーナーの 2 次 MUSCL オーバーシュート）は救えない。**2 次精度の陰解法は block DPLUR 推奨**、scalar DPLUR は 1 次（起動・ロバスト用）に限るのが実用指針（切り分けは [`time_integration-scalar-dplur-axisym-source.md`](../../plans/accepted/time_integration-scalar-dplur-axisym-source.md) / `case/29.bell_vs_conical/README.md`）。

## ループ全体 ([`main.cpp`](../../solver_density_cuda/main.cpp))

`advanceOneStep` は巨大 lambda を廃し、`StepContext`（cfg, cuda_cfg, msh, mat_ns, var, fluct,
pprobes, profiler, residual_logger, implicit_diag_logger, iStep を束ねた参照集約構造体）を
受け取る自由関数群に分解する。スキームは 3 階層構造で、dual-time が後付けできる形にする。

```
advanceOneStep(ctx):                       // dispatcher
  if (isImplicit):
     if (unsteady) advanceImplicitDualTime(ctx)   // 後続フェーズ: 今は throw
     else          advanceImplicitSteady(ctx)
  else            advanceExplicitRK(ctx)

assembleResidual(ctx, stage):              // 残差組み立ての単一情報源（旧 assembleCurrentState）
  updateInner → dependentVariables → gasProperties → applyBconds → applyRansScalarBoundaries
  → calcGradient → axisymmetricGeomTerms → limiter → ducrosSensor → turbulent_viscosity
  → convectiveFlux → ransTransport → ransGradient + ransSource
  → axisymmetricSource → viscousFlux
  → [addUnsteadyTimeTerm(ctx)]            // dual-time の BDF 物理時間項フック（定常は no-op）

advanceExplicitRK(ctx):                    // tI 1/3/4（挙動不変）
  for iloop in perStepIterationCount():
     updateVariablesInner; assembleResidual(ctx,iloop+1); logResidualSnapshot
     timeIntegration_d_wrapper(iloop); ransTimeIntegration_d_wrapper(iloop)
  updateVariablesOuter; writeStepOutputs; setDT; logOuterEnd

implicitNonlinearUpdate(ctx):              // 定常・dual-time 共有の核
  assembleResidual(ctx, 1)                 // 残差・フラックスは 1 回（ransSource が src_jac_k/ω も出力）
  setDT_d_wrapper                          // 局所擬似時間 dτ（diag の V/Δτ）
  blockDPLURSolve(ctx)                     // 下記（平均流 5 式）
  applyBlockImplicitCorrection(ctx)        // Q = Q_baseline + dq を 1 回 commit
  if scalarResidualEnabled:                // RANS(SST) のとき
     applySSTPointImplicit(ctx)            // k/ω を segregated point-implicit で更新（凍結解除）

blockDPLURSolve(ctx):                      // 古典 DPLUR 線形ソルバ（res・Q 固定）
  for iSweep in nStepInner:
     implicit_defect_correction_block_d(...)   // 固定 res + lagged dq_old → dq_new
     swapBlockImplicitCorrectionBuffers(var)

advanceImplicitSteady(ctx):                // 定常: 擬似時間=メインループ、1 更新/step の縮退形
  updateVariablesInner; logOuterBegin
  implicitNonlinearUpdate(ctx); logResidualSnapshot
  updateVariablesOuter; writeStepOutputs; setDT; logOuterEnd
```

`updateVariablesOuter_d` は外側ループ開始時に `Q_N`, `Q_M` の両方を現在値に揃え、
`updateVariablesInner_d` は内側ステージ後の `Q_M` のみを更新する。

## RK カーネル詳細

### `runge_kutta_exp_d` (Jameson 多段)

ステージ係数 `coef_N, coef_M, coef_Res` を内部ループ index `loop` で参照し、

```cpp
Q[ic] = coef_N * Q_N[ic] + coef_M * Q_M[ic]
      + coef_Res * res[ic] * dt_local[ic] / vol[ic];
```

を全成分に適用。`dt_local` は `setDT_d_wrapper` で事前に書き込まれている。

### `runge_kutta_exp_4th_d`

低 storage 4 段 4 次 RK。`loop == 0` で残差累積バッファ `res_*_m` をゼロクリアし、
各段で `res_*_m += coef_Res * res * dt_local / vol`。
`loop < 3` の中間段では `Q = Q_N + coef_DT * res * dt_local / vol`、
最終段 (`loop == 3`) で `Q = Q_N + res_*_m` として確定。

### `runge_kutta_exp_scalar_d`（スカラー k/ω の陽解法 RK ＋ point-implicit 源項）

[`scalarTransport_d.cu`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu) の
`ransTimeIntegration_d_wrapper`（[`ransTransport_d.cu`](../../solver_density_cuda/cuda_forge/ransTransport_d.cu)）が平均流 RK と同じ段で k/ω を別カーネルで積分する
（`timeIntegration==1/3` は `runge_kutta_exp_scalar_d`、`==4` は `runge_kutta_exp_scalar_4th_d`）。

RANS (SST) の消散項・輸送項は stiff なため、`timeIntegration==1/3` の更新は**残差増分のみ源項+輸送ヤコビアンで減衰**する:

```cpp
const flow_float fac = 1.0 + coef_Res * dt_l * (src_jac[ic] + transport_diag[ic] / v); // ≥ 1
rho_phi[ic] = coef_N * rho_phi_N[ic] + coef_M * rho_phi_M[ic]
            + (coef_Res * res_rho_phi[ic] * dt_l / v) / fac;
```

- `src_jac`（消散 $\beta^\*\omega,\,2\beta\omega$）は `ScalarTransportDesc` 経由で k=`src_jac_k`、ω=`src_jac_omega`
  （[`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) が毎 `assembleResidual` で出力）。
- `transport_diag`（移流+拡散の対角 $\Lambda^{T}_\phi$ [m³/s]）は `scalar_advection`/`scalar_diffusion` カーネルが
  面ループで集計（k=`transport_diag_k`、ω=`transport_diag_omega`）。$V$ で割って $[1/s]$ 化して `fac` に入る。
- 非 RANS (`model` が `sst*` 以外: 乱流なし/LES) では両者 0 で `fac=1`、従来の純陽的更新と一致（LES は無影響）。
- 減衰係数は `applySSTPointImplicit_d` の対角 $D_\phi=V/\Delta\tau+V\cdot\text{src\_jac}+\text{transport\_diag}$ に $\Delta\tau/V$ を掛けた形と整合。
- `runge_kutta_exp_scalar_4th_d`（4 次）は本減衰未適用＝RANS 非対応のまま（平均流 4 次自体が RANS 未検証）。
- 理論は [theory.md](theory.md) §"陽解法 RK での point-implicit 源項"。

## 陰解法カーネル詳細

### `implicit_defect_correction_block_d`（block DPLUR、LU-SGS $A^\pm$ 分割）

5×5 ブロック版。`block_dplur` 名前空間に補助 device 関数群を持つ。カーネルは線形 solve の内部精度
`ST` (float/double) で `template<typename ST>` 化され、状態/残差 (float) を `ST` へキャストして取り込む
(混合精度。下記「閉形式 FVS と混合精度」参照)。

- `accumulate_split_jacobian_cf<T>`（**既定の閉形式 FVS**）— $R,\Lambda,L$ の 5×5×5 三重積を作らず、
  **音響右/左固有ベクトルのみ**で $M(g)=g_2 I+(g_1-g_2)\,r_1\!\otimes l_1+(g_5-g_2)\,r_5\!\otimes l_5$
  ($=R\,\mathrm{diag}(g)\,L$) を構成し、$a^{+}=M(\Lambda^{+})$ を対角へ・$k_{\rm off}=M((-\Lambda)^{+})$ を近傍へ
  **直接畳み込む**（$a^{+}/k_{\rm off}$ を materialize しない）。`build_jacobian_split` (下記) と軸対称 (nz=0) で
  数値等価かつ ~10% 高速 (レジスタ・スピル削減)。
- `build_jacobian_split` — （旧経路・precond 版が使用）固有分解 $R,\Lambda,L$ から **対角用 $A^{+}=R\,\Lambda^{+}L$** と
  **RHS 近傍用 $K=-A^{-}=R\,(-\Lambda^{-})L=\tfrac12(|\widetilde A|-\widetilde A)$** を同時に返す
  （$\Lambda^{+}=\max(\Lambda,0)$、$-\Lambda^{-}=\max(-\Lambda,0)$、共に非負）。
- `add_identity_scaled`, `add_scaled_5x5` — $V/\Delta\tau\,I$・粘性対角・面寄与の加算。
- `solve_5x5` — 部分ピボット付き Gauss 消去（`diag` を破壊して in-place）。`|pivot| < 1e-20` でゼロ解にフォールバック。
- `multiply_add_5x5_vec` — 行列ベクトル積。

各セルで近傍寄与を集約し、対角 $D_i = V/\Delta\tau\,I + \sum_f A^{+}_f S_f + \sum_f \Lambda^{\nu}_f\,I$、
RHS = $-\mathbf R + \sum_f K_f S_f \cdot \Delta\mathbf Q_{\text{nbr}}^{\text{old}}$（$K_f=-A^{-}_f$）を構築し、
$\Delta\mathbf Q_{\text{new}} = D_i^{-1}\,\text{RHS}$ を解く。`cfg.implicitRelax` で $\Delta\mathbf Q$ を緩和。

ここで粘性対角は $\Lambda^{\nu}_f = 2\nu_f\,\dfrac{|S_f|^2}{\Delta\mathbf{cc}_f\cdot S_f}$（$\nu_f=(\mu_{\rm lam}+\mu_t)/\rho$）。
これは粘性流束 residual ([`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu)) の
法線拡散項 $\mu_f(\Delta U/|\Delta\mathbf{cc}|)\,\delta$（$\delta=|\Delta\mathbf{cc}|\,|S_f|^2/(\Delta\mathbf{cc}_f\cdot S_f)$）の
**Jacobian の大きさと整合**する（$\Lambda^{\nu}_f = 2\nu_f\,\delta/|\Delta\mathbf{cc}|$）。

> **2026-06-14 修正（粘性対角の幾何是正）**: 旧コードは粘性対角を $|S_f|\cdot(2\nu_f/\delta)$ と書いていたが、
> $\delta$ は**面積次元** ($\approx|S_f|$) なので $|S_f|$ が約分されて $\approx 2\nu_f$ に潰れ、(1) 軸対称近軸で
> 本来 $\propto r$ で消えるべき内側面の寄与を過大評価し、(2) residual に無いゼロ面積(対称/軸)面にも
> スプリアス項を載せていた（residual 側は `ip<nNormalPlanes` で境界面を除外）。これが **float block-DPLUR が
> 軸対称近軸第一セルの $U_r$ を収束させきれず固着する真因**だった。上記の residual 整合形
> $2\nu_f\,|S_f|^2/(\Delta\mathbf{cc}_f\cdot S_f)$ に是正すると $\propto r$ でゼロ面積面では消え、**float のまま固着が解消**
> （case 29 laminar conical 第一セル $U_r$: $+1.4\to+17.9$ で double solve と一致、1 次では未収束→収束）。
> 修正は `timeIntegration_d.cu` の scalar (`implicit_defect_correction_d`) / block (`implicit_defect_correction_block_d`) /
> precond (`implicit_defect_correction_block_precond_d`) の 3 箇所。**LHS のみの変更**で defect-correction の
> 定常解は不変（planar 回帰 bump で base/fix 場が $L2\sim10^{-5}$ 一致・RANS で残差レベル同一を確認）。

#### エネルギー行の熱伝導 Jacobian (`implicitThermalJacobian`、2026-10-09、既定 0 = 従来どおり)

上の粘性対角 $\Lambda^{\nu}_f I$ はエネルギー行にも同じスカラーを $\Delta(\rho E)$ に掛ける。熱伝導の残差
$k_f\,(T_j-T_i)\,\delta/|\Delta\mathbf{cc}|$ が反応するのは $T$ なので、$\rho E$ がほとんど動かずに $\rho$ と $T$ が入れ替わる
等圧のエントロピーのモードをこの対角は抑えない。普段は $V/\Delta\tau$ がそれを隠すが、方向別の擬似 dt
(`lineDtDirectional`) で $V/\Delta\tau$ が縦横比の分だけ小さくなると、300 K の等温壁の 1 層目でこのモードが育つ
(plan [tooling-nozzle-isothermal-wall-chain](../../plans/active/tooling-nozzle-isothermal-wall-chain.md) §5.1 #27)。

`time.deltaT.implicitThermalJacobian` はビットマスク (既定 0 = 従来どおり、ビット同一)。block DPLUR (`implicit_defect_correction_block_d`、ST = float/double とも)、
node 離散化だけで効く。

**ビット 1 (熱伝導の Jacobian)**: 内部の node 間面について、エネルギー行 (行 4) の $\Lambda^{\nu}_f$ を熱伝導の Jacobian に置き換える
(連続・運動量の行は従来どおり $\Lambda^{\nu}_f$):

$$
D_i[4,:] \mathrel{+}= \Lambda^{T}_f\,\frac{\gamma_i}{c_{p,i}}\,\frac{\partial e}{\partial \mathbf Q_i},\qquad
\Lambda^{T}_f=k_f\,\frac{\delta}{|\Delta\mathbf{cc}|},\qquad
k_f = f\,k_0+(1-f)\,k_1+\bigl(f\,c_{p,0}+(1-f)\,c_{p,1}\bigr)\frac{f\,\mu_{t,0}+(1-f)\,\mu_{t,1}}{Pr_t},
$$

$$
e=\frac{\rho E}{\rho}-\tfrac12|\mathbf u|^2,\quad
\frac{\partial e}{\partial\rho}=-\frac{e-\tfrac12|\mathbf u|^2}{\rho},\quad
\frac{\partial e}{\partial(\rho u_k)}=-\frac{u_k}{\rho},\quad
\frac{\partial e}{\partial(\rho E)}=\frac1\rho .
$$

$k_f$ は熱伝導の残差 ([`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) の `tc_face`) と同じ式 ($k$ = `thermCond`、$c_p$ = `cp`、$f$ = `fx`)。
物性と $\gamma$ は凍結し ($c_v=c_p/\gamma$)、熱伝導の非直交補正と近傍との結合 (非対角) は従来どおり入れない。
$\partial T/\partial\mathbf Q=(\gamma/c_p)\,\partial e/\partial\mathbf Q$ は TP でも厳密 (T を求める前に Y を正規化するので、$\rho Y$ を凍結して $\rho$ が動いても Y は変わらない)。
$e$ の基準点 (TP の `thermoHrefTemp`) によらず $\Lambda^{T}_f\,\Delta T$ の減衰になる。壁の節点と内部節点の間の面も node 間面なので、壁への熱伝導は 1 層目の対角に入る。

**ビット 4 (ビット 1 と併用、キー 5)**: 行 4 にも従来のスカラー $\Lambda^{\nu}_f$ を残し、温度の項はその上に足す。ビット 1 だけだと、行 0〜3 (スカラー) と行 4 (温度) で減衰の基準が食い違い、温度を変えない補正で行 4 だけ減衰が抜けて高速域に偽の $\Delta T$ を作る (μ_t の大きい境界層の外側で point の cfl 4 が 922 step で発散した)。

**ビット 2 (等温壁の拘束の行)**: 等温壁の節点の行 4 を単位行から $[-e_w,0,0,0,1]$ (rhs 0、$e_w=\rho E_w/\rho_w$) に替え、
$\Delta(\rho E)_w=e_w\,\Delta\rho_w$ を線形系の中で満たす (後段の壁温のピン $\rho E=\rho\,e(T_w)$ と一致させる)。原因の切り分け用。

**LHS だけの変更なので定常解 (不動点) は変わらない** (到達の保証ではない)。`lineViscCoupling: 1`・`nodeIsothermalEnergyBC: 1`・`lowMachPrecond>=2`・`blockDPLUR 0`・
`timeIntegration` ≠ 11・cell 離散化との併用は起動時に拒否する。検証と経緯は plan
[time_integration-implicit-thermal-jacobian](../../plans/active/time_integration-implicit-thermal-jacobian.md)。

> **2026-06 修正**: 旧コードは対角に $A^{+}$ ではなく $|\widetilde A|$ を、近傍に $-A^{-}$ ではなく $+|\widetilde A|$ を
> 使っていた（符号付き分割でなく絶対値の誤用）。対角が upwind 自己 Jacobian と不一致・近傍結合が逆符号となり、
> block DPLUR は収束せず発散していた。`build_jacobian_split` による $A^{+}/{-}A^{-}$ 分割でこれを修正。
> `res_*` の符号は $-\mathbf R$（陽解法 `runge_kutta_exp_d` の `Q=Q_N+\mathrm{res}\,\Delta t/V` と整合）なので、
> カーネルでは `rhs` を `res_*` で初期化し近傍寄与 $K_f S_f\,\Delta\mathbf Q_{\text{nbr}}$ を**加算**する。

**古典 DPLUR の構造（重要）**: カーネルは `dq_block_new` の生成のみを行い、
`Q`（ro..roe）を**インライン更新しない**。`blockDPLURSolve` が `nStepInner` 回 sweep を回し、
各 sweep 後にドライバ側で [`swapBlockImplicitCorrectionBuffers`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cuh)
（block 専用 swap、`.cuh` に公開）で `dq_block_old <-> dq_block_new` を入れ替える。
全 sweep 後に [`applyBlockImplicitCorrection`](../../solver_density_cuda/cuda_forge/update_d.cu)
（`applyScalarImplicitCorrection` と対称、`update_d.cu` に新設）で `Q = Q_baseline + dq_block` を
**1 度だけ** commit する。残差 `res_*` と $|\widetilde A_f|$ は sweep 中固定（matrix-free のため固定 Q から毎 sweep 再構築してよい）。

#### commit の丸め — 定常解の到達限界を決める (2026-09-23)

commit は `update_d.cu` で `ro[ic] = roN[ic] + d0`（`d0` = `dq_block_old_0`）である。
`Q` が `flow_float`（既定 float32）なので、**$|dq| < \tfrac12\,\mathrm{ULP}(Q)$ になった時点で加算は丸めで消え、
反復はそこで進まなくなる**。定常解へ近づくほど $dq$ は小さくなるので、これは**収束の到達限界**そのものである。

`updateGuardScale`（同ファイル）は $\rho$ か $e_i$ を $\alpha$ 倍未満に落とす更新だけを半減列で縮める
局所 under-relax であり、$dq/\rho \ll 1$ の領域では発動しない（$s$=1）。したがってこの丸めを緩和しない。

**実測例**（`case/56.gap_tp1187`、M7 の深いすきま、深さ $z/W>10$ の 16607 CV）:

| 量 | 値 |
| --- | --- |
| $\langle\rho\rangle$ | 0.017755 |
| 1 ULP (float32) | 1.86e-9 |
| $\langle\lvert dq\rvert\rangle$ | 2.97e-10 = **0.159 ULP** |
| 1 step で値が動く CV | **0.03 %** |
| 実効 $\langle d\rho\rangle$ / 意図した $dq$ | **0.3 %** |

この状態では、残差が系統的に残っているのに場が動かない。すきま断面を通る正味の質量流束
$\lvert\dot m\rvert$（定常解ならゼロ）は **step のべき乗則** $\propto \mathrm{step}^{-0.23}$ でしか減らず、
step を 2 倍にしても 15 % しか下がらない。同じ場を**倍精度ビルド**で継続すると**幾何級数**（25k step ごとに
2.81 分の 1）に変わり、減衰区間で 4.54e-7 から 4 桁以上落ちる。
~~「300k step で 4.56e-12」~~ **撤回** (2026-09-23): 4.56e-12 は 1.9M の谷で、2.0M では 1.48e-11 に戻る。
床は**平坦でなく**、到達最小レベルと振れ幅で書くこと (1.234e-11 ± 2.5e-12、振れ幅 1.50 倍)。

**使い方** (2026-09-24): `time.deltaT.qAccumulatorFP64: 1` (既定 0、**`time:` 直下ではない**)。
対応するのは **GPU (`gpu: 1`)・node 離散化**かつ `timeIntegration: 11` かつ `unsteady: 0`、軸対称でない、`sstEnergyIncludesK: 0`、
node 周期でない場合のみで、それ以外は起動時に拒否する。**`Qacc` は checkpoint されない** = **restart は残余を失う**。失う量は ½ ULP 分:
commit は `Qacc += dq` のあと必ず `Q = (flow_float)Qacc` とするため $\lvert Q_{acc}-Q\rvert\le\tfrac12\mathrm{ULP}(Q)$ が
構造上いつでも成り立ち、**残余は 1 ULP を超えて溜まらない**。実測の $\langle\lvert dq\rvert\rangle$ = 0.159 ULP/step から
restart の代償は case/56 の $dq$/ULP 比で**平均 3 step 分に相当する**。
**ただしこれは「$\tfrac12\mathrm{ULP}\div\langle\lvert dq\rvert\rangle$」という割り算であって、
符号相殺を含む進捗の損失そのものではない**。言えるのは —
言えるのは「このケース・この $dq$/ULP 比・100k step に 1 回の restart では、差がノイズ床の中」まで。
**ただし一般化しないこと**: $dq$ が小さいほど相当 step 数は増え、**頻回 restart は機能を丸ごと消す**
($Q=1$, $dq=0.125$ ULP を 100 回: 連続は 12 ULP 動くが毎 step restart では **0 ULP**)。case/56 で 100k から再開した軌道は、連続で回した軌道と
通算 125k–200k の 4 点すべてで **run 間ノイズ床 (絶対 1.3e-9) の中**にあり区別できない
(`run_0030_restart_cost`, 2026-09-24)。したがって `/QACC` の出力は行わない。

**切り分けの指標**: $\lvert dq\rvert/\mathrm{ULP}(Q)$ が O(1) を下回っていないか。下回っていれば、
sweep 数（`nStepInner`）を増やしても `lineImplicit` を入れても改善しない（どちらも $dq$ を精緻にするだけで、
その $dq$ が表現できない）。実測でも 4→16 sweep が ±20 % 以内、line-implicit は壁時計あたり 1 桁悪化した。

詳細と対処の設計は [`plans/active/time_integration-fp64-accumulator.md`](../../plans/active/time_integration-fp64-accumulator.md)。

#### 閉形式 FVS と混合精度 (`implicitSolvePrecision`)

`accumulate_split_jacobian_cf<T>` は固有ベクトル行列 $R,L$ を陽に作らず、$\mathrm{diag}(g)-g_2 I$ が
shear/entropy の 3 モードを消すことを使って **音響右/左固有ベクトル $r_1,r_5,l_1,l_5$ のみ**で
$M(g)=R\,\mathrm{diag}(g_1,g_2,g_2,g_2,g_5)\,L = g_2 I+(g_1-g_2)\,r_1 l_1^\top+(g_5-g_2)\,r_5 l_5^\top$
を構成し、$a^{+}=M(\Lambda^{+})$ を対角へ・$k_{\rm off}=M((-\Lambda)^{+})$ を近傍へ直接畳み込む。
$R\,\Lambda\,L$ の 5×5×5 三重積と $a^{+}/k_{\rm off}/\text{solve\_mat}$ 保持を排除し、float 陰解法で ~10% 高速
(レジスタ・スピル削減)。**軸対称・平面 (nz=0) では legacy `build_jacobian_split` と数値厳密一致**
(一般 3D は legacy の $R,L$ が厳密逆行列でない=$RL\neq I$ ため僅差だが、固有値を厳密に $\max(\lambda,0)$ とする
valid な FVS で defect-correction の定常解は不変)。

カーネルは線形 solve の内部精度 `ST` で `template<typename ST>` 化。`solverConfig` の
**`implicitSolvePrecision`** (`time.deltaT`、既定 `0`) で切替える:

- `0` (float): 状態/残差 (float) のまま float で組立・solve（既定・高速）。
- `1` (double): 状態/残差を **double へキャスト**して Jacobian 構築・5×5 solve・近傍 sweep を **double** で行い、
  補正 $\Delta\mathbf Q$ を float `dq_new` へ書戻す（混合精度 iterative refinement）。

**動機と位置づけ（重要・2026-06-14 更新）**: 当初 float32 の block-DPLUR が軸対称 近軸第一セルで平均速度 `Uy` を
収束させきれず偽固着する (laminar conical で `Uy` が物理値 $-15$ でなく $-0.6$) 問題に対し、`implicitSolvePrecision=1`
(線形 solve を double) を root-fix とした。**その後、真因は精度ではなく上記「粘性対角の幾何不整合」であることが判明**
(粘性が無い Euler は float でも固着せず、固着は粘性 LHS 由来。詳細は本節冒頭「粘性対角の幾何是正」)。
粘性対角を residual 整合形に直せば **float のまま固着が解消**し double solve は不要。
したがって `implicitSolvePrecision=1` は**根治ではなく、悪条件 LHS を倍精度で押し切る検証/保険用の手段**として残す
(幾何是正後は通常 `0` で良い)。RTX 3060 では FP64=FP32 の 1/32 ゆえ ~×2.8 遅い。
切り分け・速度の詳細は [`.github/plans/precision-mixed-axisym.md`](../../plans/archived/precision-mixed-axisym.md)
と [`.github/plans/architecture-axisym-axis-singularity.md`](../../plans/accepted/architecture-axisym-axis-singularity.md)。
現状 `blockDPLUR=1`・`lowMachPrecond` 0/1 経路のみ対応 (precond>=2 / scalar 版は float のまま)。

**低マッハ前処理 (LHS 固有値) は不採用** (2026-06 検証・実装後 revert)。`build_jacobian_split` の
固有値 `lambda[5]={U+c,U,U,U,U-c}` を前処理固有値 `U±c'` に差し替える案を実装・検証したが、
block DPLUR では**対角優位性の源である大きい音響固有値を縮めてしまい有害**（フラックス散逸前処理
単独で安定だった `eps=0.15` すら発散させ、安定 `eps` 範囲を狭めた。収束加速も根治もなし）。よって
LHS は従来の $A^\pm$（物理音速 `sonic`）のまま。低マッハ前処理は `SLAU_d` の散逸スケールにのみ適用する。
根拠・データは [`theory.md`](theory.md#低マッハ前処理固有値-weisssmith--試行したが不採用) と計画
[`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md) §9。

#### 一般EOS固有系 (TP gas, `thermalMethod==2`) — 実装・検証完了 (閉形式)

閉形式 `accumulate_split_jacobian_cf` は `inv_chi=sonic/(γ-1)` で全エンタルピーを $K+c^2/(\gamma-1)$ と再構成し、
接触波エネルギー成分に $K$ を使う。これは **CPG 専用**で TP では真の $\partial\mathbf F/\partial\mathbf Q$ と一致しない
(理論は [theory.md](theory.md) 「一般EOS固有系」、FD 検証で TP 誤差 469%)。修正方針 (plan
[`time_integration-general-eos-jacobian.md`](../../plans/accepted/time_integration-general-eos-jacobian.md)):

- **分岐保持**: `thermalMethod==0` は現行閉形式のまま (ビット不変・回帰基準)。`==2` のみ一般EOS固有系へ。
  数値 LU は閉形式と演算順序が違うため CPG ビット一致は望めない → CPG 経路を残すのが回帰構成。
- **thermo helper**: `thermo_d.cuh` に `ThermoDerivatives{p,h,cp,cv,R,kappa,chi,a2}` を 1 セル状態から返す
  `thermo_derivatives_mix` を追加。$c^2=\gamma RT,\ \kappa=\gamma-1,\ h=e+RT$ を**同一 $T$・同一組成・同一 NASA
  data・同一 datum**で評価 (`sonic[ic]`/`Ht[ic]`/`gamma[ic]` の別経路混在を避ける)。$\chi=c^2-\kappa h$。
- **固有系+LU**: 一般EOS右固有ベクトル ($\mathbf r_\mp,\mathbf r_c,\mathbf r_{sk}$, 実 $H_t$・接触 $H_t-c^2/\kappa$) を
  構築し、**部分ピボット付き double 5×5 LU** で $L=R^{-1}$ を解く ($R$ 構築・LU・$R\Lambda L$ 積算は試作段階で
  double; $H_t-c^2/\kappa=K-\chi/\kappa$ が大きな二数の差で float だと条件数誤差)。$A^\pm=R\Lambda^\pm L$。
- **セル/面の不混在**: DPLUR の対角・非対角ブロックは同一セル $i$ の $(\rho_i,\mathbf u_i,H_{t,i},c_i,\kappa_i,\chi_i)$ で
  構築 (RHS 数値流束が面/Roe 平均でも可、LHS は近似ヤコビアン)。$H_t$ だけ面・$c,\kappa$ セルの混成は不可。
- **検証 (デバイス単体 → ノズル)**: Level1 `‖LR−I‖`/`‖RΛL−A_FD‖` (CPG/TP 250・1000・高温/M=0/亜音速/音速近傍/
  超音速/一般方向)、Level2 split `A⁺+A⁻=A`・法線反転、Level3 DPLUR 行列作用、Level4 ノズル cfl 2/5/20/50/100
  (発散 step・収束率・残差床・出口M・massflux・全エンタルピー・Newton 反転回数)。**残差床は別トラック**で EOS
  反転誤差 ($\max_i|e(T_i)-e_{\text{target},i}|/\max(|e_{\text{target}}|,e_{\rm ref})$) と切り分け、機械ゼロは保証しない。
- **実装 (閉形式・確定)**: `accumulate_split_jacobian_cf` (共有ヘッダ `block_dplur_jacobian_d.cuh`) の音響右固有ベクトル
  エネルギーを実 `Ht`、左密度成分に `χ_eos=c²−κh` を加える 3 項改変 (`thermallyPerfect` 分岐、CPG ビット不変)。
  数値 L=R⁻¹ は検証参照のみ。**結果**: TP cfl 上限 2→≥100・残差 9.6e-8→4e-11、Level1/2/3 PASS (`tools/test_eos_jacobian.cpp`)。
- **precond=2 経路 (一般EOS統一・CPG 回帰確認済・TP end-to-end 未検証)**: `implicit_defect_correction_block_precond_d` は
  TP で 3 箇所 CPG 仮定 (`build_jacobian_split` の R[4]・Γ_c の g ベクトル `Htot`・∂p/∂Q の `rvec[0]` の χ 欠落) だった。
  **`build_jacobian_split` を検証済み `eos_split_jacobian_general_closed` を呼ぶ薄いラッパに統一** (CPG 専用べた打ち R/L=近似を廃止)、
  Γ_c の g/r も実 Ht・χ_eos=c²−κh に統一し `thermallyPerfect` 分岐を撤去 (CPG は χ_eos≈0 で簡約)。precond=2 が標準経路と同じ
  検証済み固有系を共有。**CPG precond=2 回帰確認** (`case/23` inviscid precond=2 eps0.15: step4000 rms_ro 1.60e-5 vs 旧 1.38e-5,
  NaN 無, 同一収束域)。ただし **TP precond=2 の収束を正で示す低マッハ TP ケースが無く** (超音速 wys は CPG でも precond=2 発散)、
  TP の end-to-end 検証は未。**TP は標準経路 (precond=0/1) が検証済・推奨**。

### `implicit_defect_correction_block_precond_d`（Phase 4: 完全 $\Gamma^{-1}A$ 前処理・`lowMachPrecond>=2`）

上の LHS 固有値前処理 (不採用) と異なり、**前処理を一貫させた別カーネル**。`blockDPLUR==1 && lowMachPrecond>=2`
のとき wrapper がこちらを起動する (既存 `implicit_defect_correction_block_d` は 0/1 専用・ビット/レジスタ不変)。

- **`lowMachPrecond=2`**: RHS のフラックス散逸も `c'` に是正 (`SLAU_d`) しつつ本 LHS 前処理を併用 → 低マッハ域の
  **収束解そのものを変える** (自励振動フロアの根治)。
- **`lowMachPrecond=3` (LHS-only)**: 本 LHS 前処理だけを行い RHS 散逸は `c_hat` のまま (`SLAU_d` 側で `==1||==2`
  のみ `c'` を使う)。本カーネルは保存形 $(\Gamma_c V/\Delta\tau' + A_c)\Delta Q=-R$ で**前処理は時間項のみ・$A_c$ と
  $R$ は非前処理**なので、**収束解は `lowMachPrecond=0` とビット一致 (保存性・解ともに不変)**。純粋に低マッハ
  剛性を除去して収束を加速する LHS 操作として使う。`setDT` の $\Delta\tau'$ 拡大も `>=2` で両モード共通に効く。

- **Sherman-Morrison 解法 (FP64 回避・高速)**: $\Gamma_c=I+\alpha g r^\top$ がランク 1 ゆえ対角ブロックは
  $D=D_0+\gamma g r^\top$ ($D_0$=V/Δτ'·I+物理 FVS+粘性+軸対称=既存 block と同形・良条件、$\gamma=(V/\Delta\tau')\alpha$)。
  $D_0$ を **float** で 2 RHS 同時 ([`solve_5x5_2rhs`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)) に解き
  $y=D_0^{-1}b,\ z=D_0^{-1}g$、$x=y-[\gamma(r^\top y)/(1+\gamma(r^\top z))]z$。悪条件 $\sim1/\beta$ は分母スカラーのみ double。
  consumer GPU (FP64=FP32/64) でも物理 block とほぼ同速 (per-step 16.9 vs 17.8ms)。$g,r,\alpha$ はカーネル内インライン。
- **時間項**: $\Gamma_c\,V/\Delta\tau'$ (上の $\gamma g r^\top$ 寄与)。$\Delta\tau'$ は `setDT_d` の
  `setDTlocal_precond_scale_d` が `dt_local *= (|u|+c)/ρ'` で拡大済 ($\rho'$=前処理スペクトル半径)。
- **フラックス分割**: **物理の厳密 FVS** をそのまま使う (対角に $a^{+}=A_c^{+}$、近傍に $k_{\rm off}=-A_c^{-}$。既存 block と同一)。
  保存形 $(\Gamma_c V/\Delta\tau' + A_c)\Delta Q=-R$ より前処理は**時間項のみ**で、$A_c$ は残差の真のヤコビアンゆえ非前処理が正しい
  ([`theory.md`](theory.md) 参照)。収束は $(\Delta\tau'/V)\Gamma_c^{-1}A_c$ の固有値 $\lambda'$ で一様に前処理される。
  ※当初フラックスも $\hat A_\Gamma=\Gamma_c^{-1}A_c$ のスペクトル半径分割で前処理したが、別系で過散逸ゆえ撤回 (上記が正)。
- **その他**: 粘性スペクトル半径・軸対称ソースヤコビアン・dual-time 物理 BDF 項 (非前処理) も倍精度で踏襲。
- **`SLAU_d`** は `lowMachPrecond==1||==2` で `c'` 散逸を使う (==2 は散逸是正を併用、==3 は RHS 散逸を触らず `c_hat`)。
  $\beta=1$ で $\Gamma_c=I$・$\Delta\tau'=\Delta\tau$・フラックス同一ゆえ現行カーネルと同一組み立て (倍精度の丸め差 ~1e-7、解一致)。

> **検証結果 (2026-06-09・採用)**: `case/23.axi_nozzle` で**低マッハ自励振動を根治**。Phase1 (物理 LHS) が発散した
> $\epsilon=0.05$ を前処理 LHS が安定化し、chamber 圧振幅 (M<0.08, 4k–20k) を 0.882%→**0.087% (定常収束・振動消滅)**。
> 安定 `cfl_pseudo` も m1~1→m2~5-7 と拡大 (収束加速は per-step 2.54× で等 wall-clock 互角。価値は根治)。データは計画 §9。

### `implicit_defect_correction_d`（scalar 対角版＝スペクトル半径）

平均流 5 式を、5×5 ブロックの代わりに**スカラー対角**で解く軽量版（`blockDPLUR == 0`）。
各セルで対角 $D = V/\Delta\tau + \sum_f(|U_n|+c+\rho^\nu)S_f$、off-diagonal $\tfrac12(|U_n|+c+\rho^\nu)S_f$ で
`dq_*_new = relax·(res_* + Σ offdiag·dq_*_nbr^{old}) / D`（5 変数とも同じスカラー $D$）を作り、
[`blockDPLURSolve`](../../solver_density_cuda/main.cpp) が `nStepInner` 回 sweep（各 sweep 後に
`swapScalarImplicitCorrectionBuffers`）、最後に [`applyScalarImplicitCorrection`](../../solver_density_cuda/update.cpp)
で `Q = Q_N + dq_*_old` を 1 回 commit する（block 版と同じ古典 DPLUR 制御フロー）。
スペクトル半径は符号不変なので block の $A^\pm$ 法線符号問題は持たない。

**block 版との比較（2026-06, `case/20.naca_ml`）**: 収束先は同一（収束場から再開すると同じ解を保持、
壁面静圧 平均 0.02% 一致）。ただし近似ヤコビアンが粗いため**安定 `cfl_pseudo` が大幅に低い**
（scalar ≲ 1〜2、block は 20〜50）。supercritical 始動では block が cfl_pseudo=20 で 4000 step / 25s で
roe→0.5 に収束する一方、scalar は cfl_pseudo=1 で 12000 step でも収束せず大きな過渡オーバーシュート
（roe ピーク ~186 vs block/explicit ~85）を示す。よって **block DPLUR が既定、scalar 対角版は
5×5 を避けたい軽量用途・低レジスタ用途向けのフォールバック**と位置づける。

### `applySSTPointImplicit_d`（SST k-ω の segregated point-implicit）

平均流 commit 後に呼ぶスカラー陰解法。源項+輸送項のヤコビアン対角を陰化して k/ω の凍結を解き、
壁近傍の陽的輸送 stiff 性を緩和して安定 `cfl_pseudo` を一桁以上引き上げる。

- 消散ヤコビアン `src_jac_k=β*ω`・`src_jac_omega=2βω` は
  [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) `rans_sst_source_d` が出力。
- 輸送ヤコビアン `transport_diag_k`/`transport_diag_omega` [m³/s] は
  [`scalarTransport_d.cu`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu) の `scalar_advection_first_order_d`
  （1次風上 $\sum_f\max(\pm\dot m,0)/\rho$）と `scalar_diffusion_first_order_d`（$\sum_f(\mu_{\text{face}}/\rho)|\delta|/dcc$）が
  面ループで atomicAdd 集計する（`ransTransport_d_wrapper` 冒頭で毎 `assembleResidual` ゼロ初期化）。
- [`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) の `applySSTPointImplicit_d` が各セルで
  $D_\phi = V/\Delta\tau + V\cdot\text{src\_jac}_\phi + \text{transport\_diag}_\phi$、
  $\delta(\rho\phi)=\text{relax}\cdot\text{res}_{\rho\phi}/D_\phi$、$\rho\phi=\max(\rho\phi^{N}+\delta,\ \text{floor})$ を適用
  （$\rho k\ge0$, $\rho\omega>0$）。`dt_local` は平均流と共用。
- [`main.cpp`](../../solver_density_cuda/main.cpp) `implicitNonlinearUpdate` で `scalarResidualEnabled` のとき
  `applyBlockImplicitCorrection` 直後に `applySSTPointImplicit` を呼ぶ。生産・近傍 ΔQ は `res` に含む lagged。
- 近傍 ΔQ 結合なしの純 point-implicit（消散+輸送の対角のみ陰化）。defect-correction のため定常解は不変。
  理論・検証は [theory.md](theory.md) §"輸送項 (移流+拡散) の point-implicit 対角"。

## 局所時間刻み

`setDT_d_wrapper` ([`setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu)) が
セル中心スペクトル半径を集計して `dt_local[ic] = CFL * V / λ_max` を書き込む。
`cfg.dt`, `cfg.cfl` で挙動を制御。

**setDT の低マッハ前処理は不採用** (2026-06 検証)。`setCFL_pln_d` のスペクトル半径の音速 `sonic` を
前処理音速 `c'` に置換する案は、低マッハ域で `dt_local` を増大させ陰解法対角 `V/Δτ` を縮め、block DPLUR の
対角優位性を崩して発散させた。よって `setCFL_pln_d` は従来の `sonic` のまま。低マッハ前処理は対流フラックス
([`convectiveFlux_d.cu`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) `SLAU_d`) の散逸スケールにのみ
適用する。詳細は計画 [`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md) §9。

## 入出力

入力: 残差 `res_ro, res_roUx, …`、$\Delta t_{\text{loc}}$ (`dt_local`)、
過去ステージの `Q_N, Q_M`、対角ヤコビアン構築用の `Ux, Uy, Uz, sonic, Ht, vis_turb`、
ステージ係数 (`cfg.coef_*`)。

出力: 更新後の `Q = (ro, roUx, roUy, roUz, roe)`。
陰解法では補助バッファ `dq_*_old/new`, `dq_block_old/new_k`, `diag_block_*`, `rhs_block_k`。

## 非定常 dual-time 陰解法（実装済み 2026-06）

`advanceImplicitDualTime` ([`main.cpp`](../../solver_density_cuda/main.cpp)) が 1 物理ステップを担当する。
使用条件 `unsteady=1, dualTime=1, timeIntegration=11, blockDPLUR=1, time.deltaT.control=0`（それ以外は throw）。

1 物理ステップの流れ:
1. **時間レベルシフト** `shiftDualTimeLevels_d_wrapper`: `roNN←roN`, `roN←ro`（$\mathbf Q^{n-1}\!\leftarrow\!\mathbf Q^n$, $\mathbf Q^n\!\leftarrow$現在）。
2. **BDF 係数**: 初回ステップ or `bdfOrder==1` は BDF1 $(a,b,c)=(1,1,0)$、以降 BDF2 $(\tfrac32,2,\tfrac12)$。
   `cfg.unsteadyDiagCoef = a/\Delta t` を設定（陰解法カーネルが対角に $V\cdot$この係数を加える）。
3. **擬似時間サブ反復** `nSubIterDualTime` 回:
   - `assembleResidual` → `addUnsteadyTimeTerm_d_wrapper(a,b,c)` で `res_* -= (V/\Delta t)(a\mathbf Q - b\mathbf Q^n + c\mathbf Q^{n-1})`
     （`include_scalar` で k/ω も）。
   - `setDT`（擬似 $\Delta\tau$）→ `blockDPLURSolve`（対角に $V\,a/\Delta t$ 込み）。
   - **in-place commit** `applyBlockImplicitCorrectionInPlace_d_wrapper`（$\mathbf Q\mathrel{+}=\delta\mathbf Q$。`roN`=$\mathbf Q^n$ 固定のため）。
     RANS のとき `applySSTPointImplicit`（こちらも in-place）。
4. `cfg.totalTime += \Delta t`、`unsteadyDiagCoef=0` リセット、出力。

実装した CUDA（[`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu)）: `addUnsteadyTimeTerm_d`,
`applyBlockImplicitCorrectionInPlace_d`, `shiftDualTimeLevels_d`。block/scalar/SST 各カーネルに対角の物理時間係数
`unsteady_diag` (`cfg.unsteadyDiagCoef`) を追加。第 2 時間レベル `roNN`/`roKNN` 系は [`variables.hpp`](../../solver_density_cuda/variables.hpp) に登録済。

## 一様体積力と質量流量一定制御（bodyForce / bodyForceCtrl）

周期境界系（チャネル・周期丘）を駆動する空間一様体積力。理論的には周期分解
$p_\mathrm{total} = -\beta(t)\,x + p'$ の非周期線形成分 $\beta$ を運動量ソース項に移項したものであり、
「平均圧力勾配駆動」と数学的に同一（局所の圧力勾配変動は解かれる $p'$ が担う）。
実装は [`bodyForce_d.cu`](../../solver_density_cuda/cuda_forge/bodyForce_d.cu):
運動量に $f_i V$、エネルギーに $(\mathbf f\cdot\mathbf u)V$ を residual へ加算（仕事項を落とすとエネルギー収支が破れる）。
閉じた周期系では体積力の仕事で系が加熱し続けるため、等温壁で排熱して温度を定常化させる。

### 固定値駆動（`bodyForce: [fx, fy, fz]`）

チャネルのように目標 $u_\tau$ から $f_x=\rho u_\tau^2/\delta$ を先験的に決められる場合はこれで足りる。

### 質量流量一定制御（`bodyForceCtrl: 1`）

周期丘のように断面が $x$ で変わる系では目標がバルク速度（=質量流量）で与えられ、抗力（壁摩擦+圧力抗力）が
先験的に分からないため、$\beta(t)$ を毎物理ステップ調整する（Benocci & Pinelli 1990 の圧縮性版）。
制御量は**体積平均 streamwise 運動量密度** $\langle\rho u_x\rangle_V = \frac1V\int \rho u_x\,dV$
（質量保存の下で統計定常では任意 $x$ 断面の質量流束と等価。目標値はケース側で
$\langle\rho u_x\rangle_V^{\,t} = \rho_0 U_b A_\mathrm{crest} L_x / V$ と換算して与える）。

更新則は 2 時刻の運動量収支から抗力を推定する deadbeat 型:

$$
f_x^{\,n} = f_x^{\,n-1} + \gamma\,\frac{M_t - 2M^n + M^{n-1}}{V\,\Delta t}
$$

（$M=\int\rho u_x\,dV$、$M_t$ は目標、$\gamma$=`bodyForceCtrlRelax` 既定 1.0。
初回は $M^{n-1}:=M^n$ とし P 項のみで立ち上がる。）
導出: $M^{n+1}=M^n+\Delta t(f_x V - D)$ で $M^{n+1}=M_t$ を課し、抗力 $D$ を前ステップの実収支
$D^{n-1}=f_x^{n-1}V-(M^n-M^{n-1})/\Delta t$ で推定して代入したもの。

実装 ([`bodyForce_d.cu`](../../solver_density_cuda/cuda_forge/bodyForce_d.cu) の `bodyForceCtrlUpdate`):

- 呼び出しは `advanceOneStep` 冒頭（物理ステップ境界）で 1 回。dual-time ではサブ反復を通じて $f_x$ 固定。
- $M^n$ は `thrust::inner_product(roUx, volume)`（once-per-step の D2H 同期 1 スカラー、コスト無視可）。
- 制御は $x$ 成分のみ。`cfg.bodyForceX` を書き換える（`bodyForceY/Z` は固定値のまま）。
- 全 CV 体積 $V$ は初回に `thrust::reduce` でキャッシュ。
- 履歴を `bodyforce_history.csv`（step, time, fx, M, M_target）へ 10 step ごとに追記。
- **`unsteady: 1` 必須**（物理 $\Delta t$ を使うため。定常局所 dt では throw）。定常 RANS の段階起動では
  固定 `bodyForce` を使い、非定常へ引き継ぐ際に `bodyForceCtrl: 1` へ切り替える運用。
- 陰解法での扱いは固定値版と同じ明示ソース（ステップ内定数のため Jacobian 不要。$\mathbf f\cdot\mathbf u$
  仕事項の対角寄与は微小につき省略）。

## 陰的更新の正値性ガード (`updateGuardAlpha`, 2026-09-02)

計画: [`plans/active/time_integration-update-positivity-guard.md`](../../plans/active/time_integration-update-positivity-guard.md)。

block-DPLUR の commit (`applyBlockImplicitCorrection_d` / 同 InPlace) で、セルごとに
縮小率 $s\in\{1,1/2,\dots,1/32\}$ (半減列) を選び $q\leftarrow q_N+s\Delta q$ とする:

$$\rho(q_N+s\Delta q)\ge\alpha\rho(q_N),\qquad e_i(q_N+s\Delta q)\ge\alpha e_i(q_N),\qquad e_i=\rho e-\tfrac12|\rho\mathbf u|^2/\rho$$

$\alpha$ = `time.deltaT.updateGuardAlpha` (既定 0.0 = OFF・ビット同一迂回、推奨 0.5)。
P でなく $(\rho, e_i)$ を使うのは TP で EOS 反転を避けつつ P 崩落と等価な検知をするため
(CPG では $P=(\gamma-1)e_i$ で厳密に等価)。5 回半減しても満たせないセルは $s=1/32$ で
commit し、EOS 床は最後の防波堤として残す。既に $\rho\le0$/$e_i\le0$ のセルは対象外。
目的は「P アンダーシュート → pMin 床洗浄 → NaN」(case/45 run_0015_cflsweep で特定した
cfl_pseudo 上限の律速) を、全域 `implicitRelax` なしにセル局所で抑えること。

## line-implicit (`lineImplicit`, 2026-09-02)

計画: [`plans/active/time_integration-line-implicit.md`](../../plans/active/time_integration-line-implicit.md)。

壁法線ライン (積層方向の CV 鎖、壁 CV 種の greedy 構築) 上の隣接結合を、point-DPLUR の
lag から **block 三重対角の直接解 (block-Thomas, 1 ライン 1 スレッド・内部 double)** に
昇格する。sweep カーネルはライン CV について点解せず diag/rhs/近傍行列 K (単位ベクトル列
抽出 = 点経路と同一 Jacobian) を保存し、各 sweep 直後の `lineThomas_d` が dq_new を上書き
する。反復に残る lag はライン外 (流れ方向) のみになる。**ただし case/45 M6 ノズルの実測では cfl 上限は不変** (律速は壁法線でなく streamwise lag と判明) で、効果は同一設定の収束 −14 %/step に留まる — 詳細と罠 (lu5 ピボット/printf 引数上限) は plan 参照。`lineImplicit: 1`
(既定 0 = 挙動不変)。blockDPLUR==1 専用、lowMachPrecond>=2 と併用不可。

### v2: factor/solve 分離・K 凍結・粘性結合・粘性 dt 割引 (2026-09-02)

計画: [`plans/accepted/time_integration-line-implicit-viscous-v2.md`](../../plans/accepted/time_integration-line-implicit-viscous-v2.md)。

- **factor/solve 分離 (常時, 厳密)**: D̃ の LU 分解・W=D̃⁻¹Knext・Kprev·W は rhs に依存しない
  ため `lineThomasFactor_d` (storeLU 時 1 回) と `lineThomasSolve_d` (毎 sweep, 保存因子で代入
  のみ) に分離。v1 モノリシックの毎 sweep 再分解が DDES で 2.44× だったコストの主犯で、分離
  だけで 1.59× に落ちる。`FORGE_LINE_MONO=1` で旧動作。因子は `line_LU_d`/`line_piv_d`。
- **`lineKFreeze: 1`**: dual-time サブ反復間で K 抽出・LU 分解を凍結 (subiter 0 のみ構築)。
  LHS 近似の強化で収束経路のみ変化 (実測で収束軌道は非凍結と一致)。コスト 1.32× まで低減。
- **`lineViscCoupling: 1`**: line 面にスカラー粘性結合 K += α·I (α=ν_eff·δ/dcc)、対角は
  2α→α で line 内に真の拡散行 [−α, 2α, −α] を完成。**圧縮性 pseudo-dt の壁法線律速は音響
  (λ_visc/λ_ac = 2ν/(Δn·c) ≪ 1) なので効果は僅差** — 意味を持つのは Δn < 2ν/c の超極薄セルのみ。
- **`lineViscousDtRelief: θ`**: on-line セルの擬似 dt 粘性スペクトル半径を (1−θ) 倍 (`setDT_d`
  で面ごとに割引、対流+音響分は残す)。θ=1 でも安定 (上と同じ理由で利得も僅差)。
- **`lineDtDirectional: 1`**: 方向別 dt — line 面 (Thomas が厳密に解く結合) の λ を音響込みで
  CFL の max から除外し、Δτ を off-line 面 (lag 側) の λ だけで決める。**注意: 除外は
  `line_prev/next` に一致する内部面のみで、壁ノードの境界半割面は残る** — 壁 CV 自身の Δτ は
  境界面の音響制約のまま、Δτ が streamwise 基準 (×AR) に伸びるのは壁の 1 つ内側以降。
  case/39 DDES で cp4 安定・ωバースト最大 9.5→6.15 (directional による低減, 帰属機構は未分離)。
- **実測の価値は pseudo-CFL 引き上げ** (case/39 ny160 DDES): point は cfl_pseudo 2 で発散、
  line は 8 まで安定。**cfl_pseudo 4 + nSub 13 で point (cp1+nSub20) の 0.88 倍時間・同品質**、
  同時間ならより深い収束・ωバースト低減。定常 (M6) と dual-time DDES で律速モードが違う点に注意。
  サブ反復収縮の残る律速候補は off-line lag / segregated SST / 2次 KEEP RHS×1次 FVS LHS の
  defect-correction 不整合の 3 者。

## 既知の TODO / 注意点

- 非定常 dual-time 陰解法（`tI==11 && unsteady==1 && dualTime==1`）は実装済（2026-06、`blockDPLUR==1` のみ、物理 $\Delta t$ 固定 `control=0`）。`implicitCorrection_d.cu` の `dualtime_explicit_d` は SLAU/Roe 用の別系統補助で本流とは独立（未使用）。
- scalar 対角陰解法（`tI==11 && blockDPLUR==0`）は有効（2026-06）。block より低 `cfl_pseudo` で収束も遅いため既定は block（上記比較参照）。
- block 陰解法でも SST(k/ω) は `applySSTPointImplicit` で segregated point-implicit 更新され、**凍結しない**（2026-06）。消散+輸送(移流+拡散)の対角を陰化、生産・近傍 ΔQ は lagged。
- `matrix mat_ns` は陰解法では未使用だが非陰解法のシグネチャに残るため `StepContext` に保持（除去は別途）。
- 旧 CPU の [`update.cpp`](../../solver_density_cuda/update.cpp) は使用されていない。
- 局所時間刻みの粘性スペクトル半径寄与は `setDT_d` 側で個別に実装されている。詳細は同ファイルを参照。
```

## 出力形式 (この形のまま)

```
結論: <次にやる一手を 1 文で>
第 1 仮説: <内容>  確度: <高/中/低>
  根拠: <ファイル:行 / run パスと数値>
  反証条件: <何が観測されたらこの仮説は誤りか>
第 2・第 3 仮説: <あれば 1 行ずつ>
判別 A/B: <変える設定 1 点、回す長さ、見る量>  → A なら … / B なら …
やらない方がよいこと: <呼び出し側が取りそうな誤った一手>
呼び出し側の前提への異議: <ブリーフの枠組み・除外判断・指標の定義で受け入れなかったものと理由。無ければ「無し」>
不足情報: <あれば>
```
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
