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

## ブリーフ (`notes/reviews/briefs/2026-10-10-faceh-floor-result.md`)

# 諮問ブリーフ: point 仕上げでの面エンタルピーの精度の A/B (plan time_integration-implicit-thermal-jacobian §6.2) の結果の解釈 (2026-10-10)

日付 2026-10-10。諮問先 codex (diagnose)。AGENTS.md の条件 7 (result の解釈を確定する前)。
plan: `plans/active/time_integration-implicit-thermal-jacobian.md` の §5.1 #5 と §6.2 (事前登録)・§6.1 (codex plan 段の採否)。
判定の出力: `case/45.isobutane_m6_d155/_band_ab/cold_pair/fh_floor_judge.json`。判定器 `case/45.isobutane_m6_d155/fh_floor_judge.py`、台本 `fh_floor.sh`。
元の諮問: `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md` (この A/B の設計の出典)。

## 観測事実

- 判定器のゲート 315 項目はすべて合格 (終了コード、step 0〜1999 の行の連続、全行の有限性、設定・入力・格子・出発の場の一致、切替の表示、最終場の有限・正)。
  自己一致: A1 の step 0 の表示と E_f(S) の相対差 2.5e-14、B1 の step 0 の表示と E_d(S) の相対差 3.3e-15。
- 同じ状態での評価の精度の効果 (E_d/E_f、`rms_roe`): S で 0.999977、A1 の 2000 step の状態で 0.999961。他の列はすべて 1.000000 (6 桁)。
  E_f(S) = 6.538800870347、E_d(S) = 6.538648285966。再評価のばらつき σ_eval は f 1.9e-14、d 2.0e-14。
- 2000 step の状態の共通の評価 (`rms_roe`):

  | | 1 本目 | 2 本目 | 3 本目 | 平均の比 B/A | 3 本の幅 (A / B) |
  | --- | --- | --- | --- | --- | --- |
  | 評価器 d: A | 6.5124 | 6.4988 | 6.4876 | 0.99013 | 0.38 % / 1.44 % |
  | 評価器 d: B | 6.4769 | 6.3845 | 6.4449 | | |
  | 評価器 f: A | 6.5127 | 6.4991 | 6.4879 | 0.99018 | 0.38 % / 1.45 % |
  | 評価器 f: B | 6.4778 | 6.3846 | 6.4456 | | |

  noise = max(σ_eval, w_A, w_B) = 1.44 %、u = 3 × noise = 4.3 %、低下 1 − r = 0.99 %。
  B の 3 本の最大 (6.477) は A の 3 本の最小 (6.488) より低い (範囲は分かれている) が、差は約 1 %。
- 全列の B/A (評価器 d、2000 step の 3 本の平均): ρ 0.983、ρu 0.996、ρv 0.983、ρE 0.990、k 1.010、ω 1.017、Y0・Y1 0.983。
- 自方式の表示 (末尾 500 step の `rms_roe` の中央値): A 6.511 / 6.554 / 6.545、B 6.459 / 6.437 / 6.456、比 B/A 0.987。範囲は分かれている。
- 出発 S に対する比 (評価器 d、`rms_roe`): A1 は 1000 step 1.0148・1500 step 0.9966・2000 step 0.9960、B1 は 0.9992・0.9739・0.9906。
- ドリフト: 6 本とも「許容内」(末尾 500 step の傾き × 500 は −0.0073〜+0.0075 桁、窓の水準の差 −0.0011〜+0.0032 桁)。傾きを除いた step ごとの揺れは 0.8〜1.0 %。
- `check_convergence --segment` (区間 main、0〜1999): 6 本とも NOT CONVERGED (stalled/plateau)。
- 判定器の出力: `NOT_SUPPORT: この期間の主要因説を支持しない; ドリフト: A この窓で有意なドリフトを検出せず・B この窓で有意なドリフトを検出せず; f と d の判定の一致: True`。

## 期待値と出典

- 事前登録 (§6.2): 支持 = 低下 − u ≥ 0.10 かつ範囲が分かれる、支持しない = 低下 + u < 0.10、それ以外は判別不能。主判定は評価器 d の `rms_roe`。
- 参考: 粘性ヤコビアン plan §6.16 (`plans/accepted/time_integration-line-viscous-jacobian.md`) では、収束していない状態 S0 (run_0183 の res_100000) で ‖R_B − R_A‖/‖R_A‖ = 5.5e-3 (エネルギーの行)。
  今回は残差の場を書き出していないので、差のベクトルのノルムは測っていない。RMS の比だけを測った。
  直交に近い差 e があるとき、RMS の比は 1 + O(|e|²/|R|²) なので、RMS の比 2e-5 は差のベクトル 0.5 % 程度と矛盾しない (推定であって測定ではない)。

## 再現条件

- 出発 S: `case/45.isobutane_m6_d155/run_0354_m9_L5cut/res_40000.h5` (L5 の point 仕上げの最終の場、`check_convergence` 全列 STALLED)。
- 設定: run_0354 と同じ (point・cfl 4・`implicitThermalJacobian` 0・ライン 0・緩和 0.7・リミッタの基準値は run_0183 で固定)。バイナリ lineM_fp64 (sha256 05ad8bdf…、FP64)。
- A = 既定、B = `FORGE_DIAG_FACE_H_DOUBLE=1`。軌道は各 2000 step・3 本 (`run_0500_fh_a1`〜`run_0505_fh_b3`)、共通の評価は 1 step の run (`run_0510`〜`run_0537_fhe_*`)。
- コード: ブランチ `feature/faceh-audit-viscjac-close`、commit c51f9fc4 (台本・判定器は run の前に commit した。forge のソースは変えていない)。

## 実施済みの操作と結果

- 1 回目の投入は prep で停止した (状態のディレクトリに化学種の記録が無く restart_field が拒否、forge は起動していない)。台本を直して (c51f9fc4) 再投入した。
- 判定の後に、評価 run の場と軌道の `nozzle.h5`・`res_0`・`res_500`、A2・A3・B2・B3 の `res_1000`・`res_1500` を消した。残っているのは A1・B1 の res_1000・1500・2000、全軌道の res_2000 と壁・出口の出力。

## 仮説 (呼び出し側の読み。確定していない)

- H1: 面エンタルピーの float の評価は、この point 仕上げの残差の床の主要因ではない (登録の判定どおり)。同じ状態での評価の差は `rms_roe` で 2e-5 で、床 (step ごとの揺れ 1 %) より 3 桁小さい。
- H2: B の約 1 % 低い水準 (3 本とも、表示でも共通の評価でも、範囲が分かれる) は、評価の差 (2e-5) では説明できない大きさなので、B の反復が少し違う状態に進んだことを示すかもしれない。
  ただし 2000 step・3 本で、事前登録の検出対象 (10 %) の 1/10 であり、登録した判定には使っていない。
- H3: float 化 (plan architecture-float-state-double-geometry) への含意は、「面エンタルピーを double に残す理由は、この point 仕上げの残差の床からは出てこない」まで。float の状態のビルドには他の float の経路 (`thermoFloat` の γ・音速、輸送物性の表引き) があるので、それらの寄与は別。

## 問い

1. 登録の判定 (NOT_SUPPORT、ドリフトは許容内) を、この記録の範囲で書いてよいか。書くときの限定 (状態・期間・評価器・FP64 ビルド) に抜けはあるか。
2. H2 の約 1 % の差をどう記録すべきか (登録外の観測として書く、書かない、追加の確認が要る)。追加の run を提案するなら、判断に要るものだけにしてほしい (長い run は要らない)。
3. H3 の float 化への含意の書き方は妥当か。元のセッション (float 化の担当) に申し送る文面として、主張と限定を示してほしい。
4. §5.1 #5 を「判定済み」として閉じてよいか。閉じる前に要るものはあるか。

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
| 5 | point 仕上げにおける残差評価の精度の監査 (移管、2026-10-09) | [time_integration-line-viscous-jacobian](../accepted/time_integration-line-viscous-jacobian.md) §6.16 で、TP の SLAU の面エンタルピーを double にすると同じ状態のエネルギーの残差の評価が ‖R_B − R_A‖/‖R_A‖ = 5.5e-3 変わった。停滞への寄与を測るなら (codex 2026-10-09): 同じ point 仕上げの保存状態・同じバイナリ (lineH 以降)・同じ設定から A = 既定、B = `FORGE_DIAG_FACE_H_DOUBLE=1` だけを変えて各 2000 step (全残差は毎 step、開始・1000・1500・2000 step の状態を保存)。末尾 500 step の全残差の水準・傾きに加え、保存した状態を共通の評価の精度で再評価する。B だけで共通の評価のエネルギーの残差が 10 % 以上下がり再評価のノイズを十分上回る → この期間の停滞への寄与を支持、自方式の表示だけが変わり共通の評価で 10 % 未満 → 主要因説を支持しない、まだ減衰中 → 床の判定は不能 (自動で延長しない)。回すかは本線の評価の後に決める (未着手)。**別のセッションへ引き継ぎ (2026-10-10)**: [`notes/sessions/2026-10-10-handoff-face-enthalpy-audit.md`](../../notes/sessions/2026-10-10-handoff-face-enthalpy-audit.md) (float 化の判断材料)。**2026-10-10 着手**: §6.2 に事前登録 (出発 run_0354 の res_40000、各腕 3 本、共通の評価は評価器 d)、codex plan 段の後に回す | F |

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

### 6.2 事前登録: point 仕上げでの面エンタルピーの精度の A/B (§5.1 #5、2026-10-10、run の前。codex plan 段の採否を反映済み)

**目的**: point 仕上げの残差の停滞に、TP の SLAU の面エンタルピーの float の評価 (`convectiveFlux_slau_d.inc.cuh` の `thermo_h_mix_f`) が寄与しているかを、
表示の残差の定義が変わっただけの効果と分けて測る。float 化 (plan architecture-float-state-double-geometry) で「どこを double に残すか」の判断材料にする。
設計は codex 諮問 ([記録](../../notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md) の「判別 A/B」) のとおりで、本節は条件と閾値を具体化する。
各腕 3 本・2000 step は、大きな効果 (10 % 以上) を探す**探索試験**であり、効果の大きさの精度は保証しない。

- **出発の状態 S**: `case/45.isobutane_m6_d155/run_0354_m9_L5cut/res_40000.h5`。速度 plan §6.14 の L5 の point 仕上げ (粘性入りのライン 11.5 万 step の後に point 4 万 step) の最終の場で、
  `check_convergence` は全列 STALLED (plateau、低下 0.0〜0.5 桁)。その末尾の `rms_roe` は step ごとに 6.51〜6.54 を行き来している (1 step で ±0.4 % 程度)。
- **設定**: run_0354 と同じ (point・cfl 4・`implicitThermalJacobian` 0・ライン 0・緩和 0.7・リミッタの基準値は run_0183 で固定・`extraFields: [res_ro]`)。
  `cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext <run> --cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --extra res_ro --field-from _fh_states/<状態>` で作る。
  軌道も評価も、状態のファイルへのシンボリックリンクだけを置いた `_fh_states/<状態>/` を渡す (prep は渡したディレクトリの最後の res を restart_field でビット一致の確認つきで移すので、入力を 1 つに固定する)。
- **バイナリ**: 全 run で `lineM_fp64` (`~/forge-linespeed-fp64`、sha256 05ad8bdf…、FP64)。台本は sha256 を確かめてから起動する。
  run_0354 自身は lineL (e10a195d…) で回したが、S は状態として使うだけで、比較はすべて同じバイナリの中で行う。
- **腕**: A = 既定、B = `FORGE_DIAG_FACE_H_DOUBLE=1` だけを変える。軌道の再現のばらつきを測るため、各腕 3 本を A1・B1・A2・B2・A3・B3 の順に逐次に回す。
  各 2000 step、全残差は毎 step (`residual_history.csv` の `outer_begin` の行)、場は 500 step ごと。
- **共通の評価**: 保存した状態から同じ設定で 1 step だけの run を回し、step 0 の `outer_begin` の行 (更新前の残差 R(Q)) を読む。評価器は f (切替なし) と d (切替あり) の 2 つ。
  状態は S、A1・B1 の 1000・1500・2000 step、A2・A3・B2・B3 の 2000 step の 11 個 × 評価器 2 つ = 22 本。
  固定した状態での実行のばらつきとして、S・A1 の 2000・B1 の 2000 を各評価器でもう 1 回ずつ評価する (6 本)。
- **主判定の量**: 評価器 d の `rms_roe` (面エンタルピーを double で評価した離散残差)。
  既定の経路は FP64 ビルドでも T・Y を float にして面エンタルピーを float で返す (`thermo_d.cuh` 1079 行付近) ので、その丸めを避ける d を両腕の状態の共通の物差しにする。
  引き継ぎメモの「A の設定の 1 step だけの run」は簡便な方法の提案で、諮問は共通の評価器を f に固定していない (codex plan 段 2026-10-10)。
  f の結果も同じ表に出し、f と d で判定が食い違えばそれも報告する。主張は「d で評価した離散残差の改善」に限り、物理の精度や収束の床の改善には広げない。
- **判定** (`fh_floor_judge.py`、run の前に commit する。E_X(状態) = 評価器 X の `rms_roe`):
  - ノイズの尺度: σ_eval = 同じ状態の 2 回の評価の相対差の最大、w_A・w_B = 2000 step の状態の 3 本の (最大 − 最小)/平均、noise = max(σ_eval, w_A, w_B)、不確かさの幅 u = 3 × noise
    (3 本の範囲の 3 倍は統計的な 3σ ではない。探索試験の目安)。r = 平均 E_d(B, 2000) / 平均 E_d(A, 2000)、低下 = 1 − r。
  - **支持** (この期間の停滞への寄与を支持): 低下 − u ≥ 0.10、かつ B の 3 本の最大 ≤ 0.9 × A の 3 本の最小。
  - **支持しない** (この期間の主要因説を支持しない): 低下 + u < 0.10。
  - **判別不能**: どちらにも当たらない (10 % の境界を不確かさの幅がまたぐ、または範囲が分かれない)。
  - 支持しない・判別不能のときに B の自方式の表示 (末尾 500 step の中央値の 3 本の平均) が A の表示の 0.9 倍以下なら、「自方式の表示では 10 % 以上下がったが、共通の評価では基準に届かない」と書く。
    「表示の定義の違いだけ」と書くのは、|1 − r| ≤ u (共通の評価で状態の差がノイズ内) かつ、表示の比が同じ状態 (A の 2000 step) の E_d/E_f と u 以内で一致するときに限る。
  - **ドリフト** (各腕の各反復): 末尾 500 step (1500〜1999) の log10(`rms_roe`) の最小二乗の傾き × 500 と、末尾の窓とその前の窓 (1000〜1499) の中央値の差 (桁) を出す。
    両方が −0.01 桁以下なら「減衰」、両方が +0.01 桁以上なら「増加」、両方の絶対値が 0.01 桁未満なら「許容内」、それ以外は「混在」。3 本で分類が揃わなければ「反復間不一致」。
    「許容内」でも結論は「この窓で有意なドリフトを検出せず」までとし、床に達したとは言わない。「減衰」は床の判定が不能 (自動で延長しない)。
    A1・B1 の 1000・1500・2000 step の E_f・E_d を S に対する比で併記する (判定には使わない)。傾きを除いた step ごとの揺れも記録する。
  - **ゲート** (1 つでも外れたら量を判定せず INVALID): 全 run の `RUN_RC` が 0。軌道の `outer_begin` の行が step 0〜1999 で一意・連続。全 run の全行の `rms_*` が有限・非負。
    評価 run の step 0 の行が 1 行で、判定に使う値が正。切替の表示 (`[DIAG] FORGE_DIAG_FACE_H_DOUBLE`) が B と評価器 d にだけある。
    `solverConfig.yaml` が run_0354 と `nStepOuter`・`outStepInterval` の他は一致。境界条件・化学種・壁の入力 (`bcondConfig.yaml`・`species_meta.yaml`・`resolved_species_*.yaml`・`wall_*.csv`・`wall_repr.json` など) が run_0354 と sha256 で一致。
    格子 (`nozzle.h5` の MESH の全データセットのハッシュ) が run_0354 と一致。全 run の `COLD_PAIR.json` の出発の場 (`field_from`・`parent_res_sha256`) が指定の状態と一致し、S の sha256 が run_0354 の res_40000 と一致。
    2000 step の場の保存量・P・T が有限で、ρ・P・T が正。最後に、A1 の step 0 の表示と E_f(S)、B1 の step 0 の表示と E_d(S) の相対差がそれぞれ 1e-6 以下
    (保守的な停止の基準。外れたら再評価のばらつきと切り分けるまで判定せず、「別の残差を読んだ」とは断定しない)。
- **記録だけにするもの**: 同じ状態での E_d/E_f (評価の精度だけで表示がどれだけ変わるか、全状態・全列)、他の列 (ρ・運動量・k・ω・化学種) の比、
  各腕の `check_convergence --segment` の VERDICT (区間 0〜1999、`CONVERGENCE_SEGMENT.txt`。短い試験の NOT CONVERGED は許し、DIVERGED・判定不能は分けて書く)。
  θ_r・Q_w などの目的量は判定しない (2000 step では遅いモードが動かない)。
- **run** (case/45 の 05xx): 軌道 `run_0500_fh_a1`・`run_0501_fh_b1`・`run_0502_fh_a2`・`run_0503_fh_b2`・`run_0504_fh_a3`・`run_0505_fh_b3`、
  評価 `run_0510`〜`run_0531_fhe_<状態>_{f,d}`・`run_0532`〜`run_0537_fhe_<状態>_{fr,dr}`。台本 `fh_floor.sh` (判定はしない。失敗したらその場で止まる)。
  場は判定の後に消す: 評価 run は h5 をすべて、軌道 run は途中の場と `nozzle.h5` を消し、res_2000 は A1・B1 だけ残す。
- **やらないこと**: 自動の延長、面エンタルピーの double を既定にすること (速度・回帰への影響を測っていない)、この A/B だけで float 化の範囲を決めること。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-09 | [2026-10-09-time_integration-implicit-thermal-jacobian-plan.md](../../notes/reviews/2026-10-09-time_integration-implicit-thermal-jacobian-plan.md) | GO-with-changes, C0/M4/m1 | 全件採用: M1 係数を残差と同じ面の k_face と c_p/γ に (§4.1)、M2 node 限定 (`isNode && has_nbr` と設定の検査)、M3 U1 (純伝導)・V0 の範囲の拡大・線形解の失敗 0 件を §6 に、M4 欠損の大きさ/ṁ_in・局所ノルム・同じ水準までの壁時計で判定 (§6 V2)、m5 原因の書き方を仮説に・V1c で熱伝導と粘性仕事を分ける |
| plan (上位の諮問) | 2026-10-09 | diagnostician (Fable、本節に要約) | 式は残差と符号・幾何で整合、Major 2 (k の定義 = codex M1 と同じ、Jacobian 単体の検査が無い = U1)・Minor 3 | 全件採用: ∂T/∂ρ は厳密と書き直し (§4.1)、V1b (拘束の行、ビット 2) を事前登録、V2 の主比較を V3 に、V0 に run_0203 の設定と ST 1 を追加、レジスタ・スピルを記録、`nodeIsothermalEnergyBC` を併用拒否。係数は diagnostician の自節点の値でなく codex の面の値 (残差の式そのもの) を採る |
| plan (§6.2 の追加) | 2026-10-10 | [2026-10-10-time_integration-implicit-thermal-jacobian-plan.md](../../notes/reviews/2026-10-10-time_integration-implicit-thermal-jacobian-plan.md) | GO-with-changes, C0/M5/m1 | 全件採用 (根拠の箇所を確かめ、判定器は人工入力で再試験): M1 終了コード・行の連続性・全行の有限性・最終場の全保存量のゲートを追加し、外れたら INVALID、評価 run の場は判定後に消す。M2 軌道も `_fh_states/s0` から始め、出発の場の sha256・境界条件などの入力・格子のハッシュを run_0354 と照合。M3 床の判定を各反復の「減衰・増加・許容内・混在」と窓の水準の差に分け、「許容内」でも「有意なドリフトを検出せず」までに限定。M4 不確かさの幅 u = 3 × noise を「支持しない」側にも課し、10 % の境界をまたげば判別不能。M5 「表示の定義の違いだけ」は、共通の評価の差がノイズ内かつ同じ状態の f/d の差で表示の差を説明できるときに限る。m6 各軌道の `check_convergence --segment` を台本で保存し判定器が記録。d を主の評価器にする方針は codex も妥当とした |

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
- `2026-10-10` — §5.1 #5 (point 仕上げでの面エンタルピーの精度の A/B) を引き継いで §6.2 に事前登録。codex plan 段 (GO-with-changes、C0/M5/m1) を全件採用してゲートと判定の分岐を直した。台本 `fh_floor.sh`・判定 `fh_floor_judge.py` は run の前に commit。粘性 Jacobian plan の accepted への移動に合わせてリンクを更新。
```

## 参考: `case/45.isobutane_m6_d155/_band_ab/cold_pair/fh_floor_judge.json`

```
{
 "plan": "time_integration-implicit-thermal-jacobian §6.2",
 "quantity": "rms_roe",
 "gates": [
  {
   "check": "S の sha256 が run_0354/res_40000.h5 と一致",
   "ok": true,
   "value": "4118dc99ce143303"
  },
  {
   "check": "run_0500_fh_a1: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0500_fh_a1: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0500_fh_a1: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0500_fh_a1: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0500_fh_a1: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0500_fh_a1: outer_begin の行が step 0〜1999 で一意・連続",
   "ok": true,
   "value": [
    0,
    1999,
    2000,
    []
   ]
  },
  {
   "check": "run_0500_fh_a1: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0500_fh_a1: 切替の表示 = 0",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0500_fh_a1: res_2000 の保存量・P・T が有限",
   "ok": true,
   "value": {
    "ro": 0,
    "roUx": 0,
    "roUy": 0,
    "roUz": 0,
    "roe": 0,
    "roK": 0,
    "roOmega": 0,
    "roY0": 0,
    "roY1": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0500_fh_a1: res_2000 の ρ・P・T が正",
   "ok": true,
   "value": {
    "ro": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0501_fh_b1: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0501_fh_b1: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0501_fh_b1: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0501_fh_b1: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0501_fh_b1: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0501_fh_b1: outer_begin の行が step 0〜1999 で一意・連続",
   "ok": true,
   "value": [
    0,
    1999,
    2000,
    []
   ]
  },
  {
   "check": "run_0501_fh_b1: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0501_fh_b1: 切替の表示 = 1",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0501_fh_b1: res_2000 の保存量・P・T が有限",
   "ok": true,
   "value": {
    "ro": 0,
    "roUx": 0,
    "roUy": 0,
    "roUz": 0,
    "roe": 0,
    "roK": 0,
    "roOmega": 0,
    "roY0": 0,
    "roY1": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0501_fh_b1: res_2000 の ρ・P・T が正",
   "ok": true,
   "value": {
    "ro": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0502_fh_a2: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0502_fh_a2: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0502_fh_a2: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0502_fh_a2: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0502_fh_a2: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0502_fh_a2: outer_begin の行が step 0〜1999 で一意・連続",
   "ok": true,
   "value": [
    0,
    1999,
    2000,
    []
   ]
  },
  {
   "check": "run_0502_fh_a2: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0502_fh_a2: 切替の表示 = 0",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0502_fh_a2: res_2000 の保存量・P・T が有限",
   "ok": true,
   "value": {
    "ro": 0,
    "roUx": 0,
    "roUy": 0,
    "roUz": 0,
    "roe": 0,
    "roK": 0,
    "roOmega": 0,
    "roY0": 0,
    "roY1": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0502_fh_a2: res_2000 の ρ・P・T が正",
   "ok": true,
   "value": {
    "ro": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0503_fh_b2: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0503_fh_b2: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0503_fh_b2: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0503_fh_b2: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0503_fh_b2: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0503_fh_b2: outer_begin の行が step 0〜1999 で一意・連続",
   "ok": true,
   "value": [
    0,
    1999,
    2000,
    []
   ]
  },
  {
   "check": "run_0503_fh_b2: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0503_fh_b2: 切替の表示 = 1",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0503_fh_b2: res_2000 の保存量・P・T が有限",
   "ok": true,
   "value": {
    "ro": 0,
    "roUx": 0,
    "roUy": 0,
    "roUz": 0,
    "roe": 0,
    "roK": 0,
    "roOmega": 0,
    "roY0": 0,
    "roY1": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0503_fh_b2: res_2000 の ρ・P・T が正",
   "ok": true,
   "value": {
    "ro": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0504_fh_a3: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0504_fh_a3: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0504_fh_a3: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0504_fh_a3: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0504_fh_a3: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0504_fh_a3: outer_begin の行が step 0〜1999 で一意・連続",
   "ok": true,
   "value": [
    0,
    1999,
    2000,
    []
   ]
  },
  {
   "check": "run_0504_fh_a3: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0504_fh_a3: 切替の表示 = 0",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0504_fh_a3: res_2000 の保存量・P・T が有限",
   "ok": true,
   "value": {
    "ro": 0,
    "roUx": 0,
    "roUy": 0,
    "roUz": 0,
    "roe": 0,
    "roK": 0,
    "roOmega": 0,
    "roY0": 0,
    "roY1": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0504_fh_a3: res_2000 の ρ・P・T が正",
   "ok": true,
   "value": {
    "ro": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0505_fh_b3: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0505_fh_b3: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0505_fh_b3: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0505_fh_b3: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0505_fh_b3: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0505_fh_b3: outer_begin の行が step 0〜1999 で一意・連続",
   "ok": true,
   "value": [
    0,
    1999,
    2000,
    []
   ]
  },
  {
   "check": "run_0505_fh_b3: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0505_fh_b3: 切替の表示 = 1",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0505_fh_b3: res_2000 の保存量・P・T が有限",
   "ok": true,
   "value": {
    "ro": 0,
    "roUx": 0,
    "roUy": 0,
    "roUz": 0,
    "roe": 0,
    "roK": 0,
    "roOmega": 0,
    "roY0": 0,
    "roY1": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0505_fh_b3: res_2000 の ρ・P・T が正",
   "ok": true,
   "value": {
    "ro": 0,
    "P": 0,
    "T": 0
   }
  },
  {
   "check": "run_0510_fhe_s0_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0510_fhe_s0_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0510_fhe_s0_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0510_fhe_s0_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0510_fhe_s0_f: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0510_fhe_s0_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0510_fhe_s0_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0510_fhe_s0_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0510_fhe_s0_f: 判定に使う値が正",
   "ok": true,
   "value": 6.538800870347329
  },
  {
   "check": "run_0511_fhe_s0_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0511_fhe_s0_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0511_fhe_s0_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0511_fhe_s0_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0511_fhe_s0_d: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0511_fhe_s0_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0511_fhe_s0_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0511_fhe_s0_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0511_fhe_s0_d: 判定に使う値が正",
   "ok": true,
   "value": 6.538648285965989
  },
  {
   "check": "run_0512_fhe_a1k1000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0512_fhe_a1k1000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0512_fhe_a1k1000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0512_fhe_a1k1000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0512_fhe_a1k1000_f: 出発の場が状態 a1k1000",
   "ok": true,
   "value": [
    "a1k1000",
    "res_1000.h5",
    "9b47db7a07bdb472"
   ]
  },
  {
   "check": "run_0512_fhe_a1k1000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0512_fhe_a1k1000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0512_fhe_a1k1000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0512_fhe_a1k1000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.635434953120642
  },
  {
   "check": "run_0513_fhe_a1k1000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0513_fhe_a1k1000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0513_fhe_a1k1000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0513_fhe_a1k1000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0513_fhe_a1k1000_d: 出発の場が状態 a1k1000",
   "ok": true,
   "value": [
    "a1k1000",
    "res_1000.h5",
    "9b47db7a07bdb472"
   ]
  },
  {
   "check": "run_0513_fhe_a1k1000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0513_fhe_a1k1000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0513_fhe_a1k1000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0513_fhe_a1k1000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.635502449588447
  },
  {
   "check": "run_0514_fhe_a1k1500_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0514_fhe_a1k1500_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0514_fhe_a1k1500_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0514_fhe_a1k1500_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0514_fhe_a1k1500_f: 出発の場が状態 a1k1500",
   "ok": true,
   "value": [
    "a1k1500",
    "res_1500.h5",
    "b6b495c034e96d5c"
   ]
  },
  {
   "check": "run_0514_fhe_a1k1500_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0514_fhe_a1k1500_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0514_fhe_a1k1500_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0514_fhe_a1k1500_f: 判定に使う値が正",
   "ok": true,
   "value": 6.516350494521116
  },
  {
   "check": "run_0515_fhe_a1k1500_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0515_fhe_a1k1500_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0515_fhe_a1k1500_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0515_fhe_a1k1500_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0515_fhe_a1k1500_d: 出発の場が状態 a1k1500",
   "ok": true,
   "value": [
    "a1k1500",
    "res_1500.h5",
    "b6b495c034e96d5c"
   ]
  },
  {
   "check": "run_0515_fhe_a1k1500_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0515_fhe_a1k1500_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0515_fhe_a1k1500_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0515_fhe_a1k1500_d: 判定に使う値が正",
   "ok": true,
   "value": 6.516568875096715
  },
  {
   "check": "run_0516_fhe_a1k2000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0516_fhe_a1k2000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0516_fhe_a1k2000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0516_fhe_a1k2000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0516_fhe_a1k2000_f: 出発の場が状態 a1k2000",
   "ok": true,
   "value": [
    "a1k2000",
    "res_2000.h5",
    "ebe7e7a270378444"
   ]
  },
  {
   "check": "run_0516_fhe_a1k2000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0516_fhe_a1k2000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0516_fhe_a1k2000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0516_fhe_a1k2000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.512654161893749
  },
  {
   "check": "run_0517_fhe_a1k2000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0517_fhe_a1k2000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0517_fhe_a1k2000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0517_fhe_a1k2000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0517_fhe_a1k2000_d: 出発の場が状態 a1k2000",
   "ok": true,
   "value": [
    "a1k2000",
    "res_2000.h5",
    "ebe7e7a270378444"
   ]
  },
  {
   "check": "run_0517_fhe_a1k2000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0517_fhe_a1k2000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0517_fhe_a1k2000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0517_fhe_a1k2000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.512400268351112
  },
  {
   "check": "run_0518_fhe_b1k1000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0518_fhe_b1k1000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0518_fhe_b1k1000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0518_fhe_b1k1000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0518_fhe_b1k1000_f: 出発の場が状態 b1k1000",
   "ok": true,
   "value": [
    "b1k1000",
    "res_1000.h5",
    "9e8c8bf45b17a82a"
   ]
  },
  {
   "check": "run_0518_fhe_b1k1000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0518_fhe_b1k1000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0518_fhe_b1k1000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0518_fhe_b1k1000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.533379514685974
  },
  {
   "check": "run_0519_fhe_b1k1000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0519_fhe_b1k1000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0519_fhe_b1k1000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0519_fhe_b1k1000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0519_fhe_b1k1000_d: 出発の場が状態 b1k1000",
   "ok": true,
   "value": [
    "b1k1000",
    "res_1000.h5",
    "9e8c8bf45b17a82a"
   ]
  },
  {
   "check": "run_0519_fhe_b1k1000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0519_fhe_b1k1000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0519_fhe_b1k1000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0519_fhe_b1k1000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.533289979154832
  },
  {
   "check": "run_0520_fhe_b1k1500_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0520_fhe_b1k1500_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0520_fhe_b1k1500_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0520_fhe_b1k1500_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0520_fhe_b1k1500_f: 出発の場が状態 b1k1500",
   "ok": true,
   "value": [
    "b1k1500",
    "res_1500.h5",
    "66ccf5edd30556e2"
   ]
  },
  {
   "check": "run_0520_fhe_b1k1500_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0520_fhe_b1k1500_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0520_fhe_b1k1500_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0520_fhe_b1k1500_f: 判定に使う値が正",
   "ok": true,
   "value": 6.368478885485925
  },
  {
   "check": "run_0521_fhe_b1k1500_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0521_fhe_b1k1500_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0521_fhe_b1k1500_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0521_fhe_b1k1500_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0521_fhe_b1k1500_d: 出発の場が状態 b1k1500",
   "ok": true,
   "value": [
    "b1k1500",
    "res_1500.h5",
    "66ccf5edd30556e2"
   ]
  },
  {
   "check": "run_0521_fhe_b1k1500_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0521_fhe_b1k1500_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0521_fhe_b1k1500_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0521_fhe_b1k1500_d: 判定に使う値が正",
   "ok": true,
   "value": 6.367791335071682
  },
  {
   "check": "run_0522_fhe_b1k2000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0522_fhe_b1k2000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0522_fhe_b1k2000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0522_fhe_b1k2000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0522_fhe_b1k2000_f: 出発の場が状態 b1k2000",
   "ok": true,
   "value": [
    "b1k2000",
    "res_2000.h5",
    "a3aa192d5c0a72d8"
   ]
  },
  {
   "check": "run_0522_fhe_b1k2000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0522_fhe_b1k2000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0522_fhe_b1k2000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0522_fhe_b1k2000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.477839300828546
  },
  {
   "check": "run_0523_fhe_b1k2000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0523_fhe_b1k2000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0523_fhe_b1k2000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0523_fhe_b1k2000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0523_fhe_b1k2000_d: 出発の場が状態 b1k2000",
   "ok": true,
   "value": [
    "b1k2000",
    "res_2000.h5",
    "a3aa192d5c0a72d8"
   ]
  },
  {
   "check": "run_0523_fhe_b1k2000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0523_fhe_b1k2000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0523_fhe_b1k2000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0523_fhe_b1k2000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.476892277117087
  },
  {
   "check": "run_0524_fhe_a2k2000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0524_fhe_a2k2000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0524_fhe_a2k2000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0524_fhe_a2k2000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0524_fhe_a2k2000_f: 出発の場が状態 a2k2000",
   "ok": true,
   "value": [
    "a2k2000",
    "res_2000.h5",
    "0ce8c2c64d5da252"
   ]
  },
  {
   "check": "run_0524_fhe_a2k2000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0524_fhe_a2k2000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0524_fhe_a2k2000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0524_fhe_a2k2000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.499090431533056
  },
  {
   "check": "run_0525_fhe_a2k2000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0525_fhe_a2k2000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0525_fhe_a2k2000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0525_fhe_a2k2000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0525_fhe_a2k2000_d: 出発の場が状態 a2k2000",
   "ok": true,
   "value": [
    "a2k2000",
    "res_2000.h5",
    "0ce8c2c64d5da252"
   ]
  },
  {
   "check": "run_0525_fhe_a2k2000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0525_fhe_a2k2000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0525_fhe_a2k2000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0525_fhe_a2k2000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.498760427167691
  },
  {
   "check": "run_0526_fhe_a3k2000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0526_fhe_a3k2000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0526_fhe_a3k2000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0526_fhe_a3k2000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0526_fhe_a3k2000_f: 出発の場が状態 a3k2000",
   "ok": true,
   "value": [
    "a3k2000",
    "res_2000.h5",
    "c6547a9058102336"
   ]
  },
  {
   "check": "run_0526_fhe_a3k2000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0526_fhe_a3k2000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0526_fhe_a3k2000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0526_fhe_a3k2000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.487879493511477
  },
  {
   "check": "run_0527_fhe_a3k2000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0527_fhe_a3k2000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0527_fhe_a3k2000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0527_fhe_a3k2000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0527_fhe_a3k2000_d: 出発の場が状態 a3k2000",
   "ok": true,
   "value": [
    "a3k2000",
    "res_2000.h5",
    "c6547a9058102336"
   ]
  },
  {
   "check": "run_0527_fhe_a3k2000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0527_fhe_a3k2000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0527_fhe_a3k2000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0527_fhe_a3k2000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.487647427214041
  },
  {
   "check": "run_0528_fhe_b2k2000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0528_fhe_b2k2000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0528_fhe_b2k2000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0528_fhe_b2k2000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0528_fhe_b2k2000_f: 出発の場が状態 b2k2000",
   "ok": true,
   "value": [
    "b2k2000",
    "res_2000.h5",
    "5600200b9d3bcbd0"
   ]
  },
  {
   "check": "run_0528_fhe_b2k2000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0528_fhe_b2k2000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0528_fhe_b2k2000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0528_fhe_b2k2000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.384635482129856
  },
  {
   "check": "run_0529_fhe_b2k2000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0529_fhe_b2k2000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0529_fhe_b2k2000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0529_fhe_b2k2000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0529_fhe_b2k2000_d: 出発の場が状態 b2k2000",
   "ok": true,
   "value": [
    "b2k2000",
    "res_2000.h5",
    "5600200b9d3bcbd0"
   ]
  },
  {
   "check": "run_0529_fhe_b2k2000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0529_fhe_b2k2000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0529_fhe_b2k2000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0529_fhe_b2k2000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.38448525790721
  },
  {
   "check": "run_0530_fhe_b3k2000_f: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0530_fhe_b3k2000_f: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0530_fhe_b3k2000_f: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0530_fhe_b3k2000_f: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0530_fhe_b3k2000_f: 出発の場が状態 b3k2000",
   "ok": true,
   "value": [
    "b3k2000",
    "res_2000.h5",
    "14e1a88b8076e75a"
   ]
  },
  {
   "check": "run_0530_fhe_b3k2000_f: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0530_fhe_b3k2000_f: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0530_fhe_b3k2000_f: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0530_fhe_b3k2000_f: 判定に使う値が正",
   "ok": true,
   "value": 6.445575039217898
  },
  {
   "check": "run_0531_fhe_b3k2000_d: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0531_fhe_b3k2000_d: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0531_fhe_b3k2000_d: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0531_fhe_b3k2000_d: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0531_fhe_b3k2000_d: 出発の場が状態 b3k2000",
   "ok": true,
   "value": [
    "b3k2000",
    "res_2000.h5",
    "14e1a88b8076e75a"
   ]
  },
  {
   "check": "run_0531_fhe_b3k2000_d: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0531_fhe_b3k2000_d: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0531_fhe_b3k2000_d: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0531_fhe_b3k2000_d: 判定に使う値が正",
   "ok": true,
   "value": 6.444935810214566
  },
  {
   "check": "run_0532_fhe_s0_fr: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0532_fhe_s0_fr: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0532_fhe_s0_fr: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0532_fhe_s0_fr: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0532_fhe_s0_fr: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0532_fhe_s0_fr: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0532_fhe_s0_fr: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0532_fhe_s0_fr: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0532_fhe_s0_fr: 判定に使う値が正",
   "ok": true,
   "value": 6.538800870347452
  },
  {
   "check": "run_0533_fhe_s0_dr: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0533_fhe_s0_dr: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0533_fhe_s0_dr: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0533_fhe_s0_dr: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0533_fhe_s0_dr: 出発の場が状態 s0",
   "ok": true,
   "value": [
    "s0",
    "res_40000.h5",
    "4118dc99ce143303"
   ]
  },
  {
   "check": "run_0533_fhe_s0_dr: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0533_fhe_s0_dr: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0533_fhe_s0_dr: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0533_fhe_s0_dr: 判定に使う値が正",
   "ok": true,
   "value": 6.538648285966021
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: 出発の場が状態 a1k2000",
   "ok": true,
   "value": [
    "a1k2000",
    "res_2000.h5",
    "ebe7e7a270378444"
   ]
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0534_fhe_a1k2000_fr: 判定に使う値が正",
   "ok": true,
   "value": 6.512654161893741
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: 出発の場が状態 a1k2000",
   "ok": true,
   "value": [
    "a1k2000",
    "res_2000.h5",
    "ebe7e7a270378444"
   ]
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0535_fhe_a1k2000_dr: 判定に使う値が正",
   "ok": true,
   "value": 6.512400268350982
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: 出発の場が状態 b1k2000",
   "ok": true,
   "value": [
    "b1k2000",
    "res_2000.h5",
    "a3aa192d5c0a72d8"
   ]
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: 切替の表示 = f",
   "ok": true,
   "value": false
  },
  {
   "check": "run_0536_fhe_b1k2000_fr: 判定に使う値が正",
   "ok": true,
   "value": 6.477839300828452
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: RUN_RC = 0",
   "ok": true,
   "value": "0"
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: 境界条件・化学種・壁などの入力が run_0354 と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: 格子 (MESH) が run_0354 と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: 出発の場が状態 b1k2000",
   "ok": true,
   "value": [
    "b1k2000",
    "res_2000.h5",
    "a3aa192d5c0a72d8"
   ]
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: step 0 の outer_begin が 1 行",
   "ok": true,
   "value": [
    [
     0
    ],
    []
   ]
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: 全行の rms_* が有限・非負",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: 切替の表示 = d",
   "ok": true,
   "value": true
  },
  {
   "check": "run_0537_fhe_b1k2000_dr: 判定に使う値が正",
   "ok": true,
   "value": 6.47689227711699
  },
  {
   "check": "A1 の step 0 の表示と E_f(S) の相対差 ≤ 1e-6",
   "ok": true,
   "value": 2.4585592550271594e-14
  },
  {
   "check": "B1 の step 0 の表示と E_d(S) の相対差 ≤ 1e-6",
   "ok": true,
   "value": 3.2600441468238118e-15
  }
 ],
 "gates_ok": true,
 "E": {
  "s0|f|1": {
   "rms_ro": 3.69970431394538e-06,
   "rms_roUx": 0.003121715636196732,
   "rms_roUy": 0.0003868401681075493,
   "rms_roUz": 0.0,
   "rms_roe": 6.538800870347329,
   "rms_roK": 0.0002899758303757727,
   "rms_roOmega": 3251.612630613857,
   "rms_roY0": 3.382192595459101e-06,
   "rms_roY1": 3.174273951971575e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "s0|d|1": {
   "rms_ro": 3.699704313945344e-06,
   "rms_roUx": 0.003121715636196659,
   "rms_roUy": 0.0003868401681075231,
   "rms_roUz": 0.0,
   "rms_roe": 6.538648285965989,
   "rms_roK": 0.0002899758303757726,
   "rms_roOmega": 3251.612630613712,
   "rms_roY0": 3.382192595459111e-06,
   "rms_roY1": 3.174273951971492e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k1000|f|1": {
   "rms_ro": 3.755429566987045e-06,
   "rms_roUx": 0.003160963918424999,
   "rms_roUy": 0.0003914697495266336,
   "rms_roUz": 0.0,
   "rms_roe": 6.635434953120642,
   "rms_roK": 0.0002881543423578108,
   "rms_roOmega": 3288.772626339693,
   "rms_roY0": 3.433136327965832e-06,
   "rms_roY1": 3.222085943315155e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k1000|d|1": {
   "rms_ro": 3.755429566986994e-06,
   "rms_roUx": 0.003160963918425157,
   "rms_roUy": 0.0003914697495266874,
   "rms_roUz": 0.0,
   "rms_roe": 6.635502449588447,
   "rms_roK": 0.0002881543423578109,
   "rms_roOmega": 3288.772626339244,
   "rms_roY0": 3.433136327965841e-06,
   "rms_roY1": 3.22208594331518e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k1500|f|1": {
   "rms_ro": 3.686598714338017e-06,
   "rms_roUx": 0.003122656457294422,
   "rms_roUy": 0.0003985235225601041,
   "rms_roUz": 0.0,
   "rms_roe": 6.516350494521116,
   "rms_roK": 0.0002882781684617978,
   "rms_roOmega": 3231.960456931614,
   "rms_roY0": 3.370212722178093e-06,
   "rms_roY1": 3.163030535563602e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k1500|d|1": {
   "rms_ro": 3.686598714338038e-06,
   "rms_roUx": 0.003122656457294405,
   "rms_roUy": 0.0003985235225601257,
   "rms_roUz": 0.0,
   "rms_roe": 6.516568875096715,
   "rms_roK": 0.0002882781684617991,
   "rms_roOmega": 3231.9604569316,
   "rms_roY0": 3.370212722178132e-06,
   "rms_roY1": 3.163030535563595e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k2000|f|1": {
   "rms_ro": 3.674254260713269e-06,
   "rms_roUx": 0.003113306827588208,
   "rms_roUy": 0.0003924010156009334,
   "rms_roUz": 0.0,
   "rms_roe": 6.512654161893749,
   "rms_roK": 0.0002871064075878131,
   "rms_roOmega": 3276.54771178772,
   "rms_roY0": 3.358925580652281e-06,
   "rms_roY1": 3.152437265569644e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k2000|d|1": {
   "rms_ro": 3.674254260713299e-06,
   "rms_roUx": 0.003113306827588387,
   "rms_roUy": 0.0003924010156009257,
   "rms_roUz": 0.0,
   "rms_roe": 6.512400268351112,
   "rms_roK": 0.0002871064075878134,
   "rms_roOmega": 3276.5477117876,
   "rms_roY0": 3.358925580652279e-06,
   "rms_roY1": 3.152437265569583e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k1000|f|1": {
   "rms_ro": 3.650048441134663e-06,
   "rms_roUx": 0.003144080062173102,
   "rms_roUy": 0.0003853327938060909,
   "rms_roUz": 0.0,
   "rms_roe": 6.533379514685974,
   "rms_roK": 0.0002892383130284789,
   "rms_roOmega": 3297.266485590319,
   "rms_roY0": 3.336796497166009e-06,
   "rms_roY1": 3.131668556718653e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k1000|d|1": {
   "rms_ro": 3.650048441134652e-06,
   "rms_roUx": 0.003144080062172927,
   "rms_roUy": 0.0003853327938060475,
   "rms_roUz": 0.0,
   "rms_roe": 6.533289979154832,
   "rms_roK": 0.0002892383130284795,
   "rms_roOmega": 3297.266485590329,
   "rms_roY0": 3.33679649716603e-06,
   "rms_roY1": 3.131668556718618e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k1500|f|1": {
   "rms_ro": 3.553368270255621e-06,
   "rms_roUx": 0.003089676429589154,
   "rms_roUy": 0.0003808002146701567,
   "rms_roUz": 0.0,
   "rms_roe": 6.368478885485925,
   "rms_roK": 0.0002902154549669925,
   "rms_roOmega": 3242.096624600204,
   "rms_roY0": 3.248406473211662e-06,
   "rms_roY1": 3.04871226645922e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k1500|d|1": {
   "rms_ro": 3.553368270255678e-06,
   "rms_roUx": 0.00308967642958931,
   "rms_roUy": 0.0003808002146702004,
   "rms_roUz": 0.0,
   "rms_roe": 6.367791335071682,
   "rms_roK": 0.0002902154549669916,
   "rms_roOmega": 3242.096624600061,
   "rms_roY0": 3.248406473211661e-06,
   "rms_roY1": 3.048712266459239e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k2000|f|1": {
   "rms_ro": 3.643390393577328e-06,
   "rms_roUx": 0.003102712452654722,
   "rms_roUy": 0.000393774213941897,
   "rms_roUz": 0.0,
   "rms_roe": 6.477839300828546,
   "rms_roK": 0.0002896141351855703,
   "rms_roOmega": 3320.033940488261,
   "rms_roY0": 3.330706552599267e-06,
   "rms_roY1": 3.125952988530733e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k2000|d|1": {
   "rms_ro": 3.643390393577253e-06,
   "rms_roUx": 0.003102712452654942,
   "rms_roUy": 0.0003937742139419189,
   "rms_roUz": 0.0,
   "rms_roe": 6.476892277117087,
   "rms_roK": 0.0002896141351855706,
   "rms_roOmega": 3320.033940488263,
   "rms_roY0": 3.330706552599231e-06,
   "rms_roY1": 3.125952988530696e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a2k2000|f|1": {
   "rms_ro": 3.671892928482832e-06,
   "rms_roUx": 0.003104443272456606,
   "rms_roUy": 0.0003908087683983788,
   "rms_roUz": 0.0,
   "rms_roe": 6.499090431533056,
   "rms_roK": 0.0002867848204002292,
   "rms_roOmega": 3256.191270929935,
   "rms_roY0": 3.356766150863537e-06,
   "rms_roY1": 3.150410585676603e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a2k2000|d|1": {
   "rms_ro": 3.671892928482853e-06,
   "rms_roUx": 0.003104443272456765,
   "rms_roUy": 0.0003908087683984255,
   "rms_roUz": 0.0,
   "rms_roe": 6.498760427167691,
   "rms_roK": 0.0002867848204002285,
   "rms_roOmega": 3256.191270930059,
   "rms_roY0": 3.356766150863552e-06,
   "rms_roY1": 3.150410585676629e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a3k2000|f|1": {
   "rms_ro": 3.655313429217482e-06,
   "rms_roUx": 0.003109620264185232,
   "rms_roUy": 0.0003865463649431292,
   "rms_roUz": 0.0,
   "rms_roe": 6.487879493511477,
   "rms_roK": 0.0002867768792274841,
   "rms_roOmega": 3257.612928234795,
   "rms_roY0": 3.341610769315766e-06,
   "rms_roY1": 3.136186873832958e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a3k2000|d|1": {
   "rms_ro": 3.655313429217531e-06,
   "rms_roUx": 0.003109620264185058,
   "rms_roUy": 0.0003865463649431234,
   "rms_roUz": 0.0,
   "rms_roe": 6.487647427214041,
   "rms_roK": 0.0002867768792274836,
   "rms_roOmega": 3257.612928234576,
   "rms_roY0": 3.341610769315745e-06,
   "rms_roY1": 3.136186873832961e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b2k2000|f|1": {
   "rms_ro": 3.561443180723182e-06,
   "rms_roUx": 0.00308066760545585,
   "rms_roUy": 0.000373848255826916,
   "rms_roUz": 0.0,
   "rms_roe": 6.384635482129856,
   "rms_roK": 0.0002897123690935014,
   "rms_roOmega": 3315.785139385565,
   "rms_roY0": 3.255789156902594e-06,
   "rms_roY1": 3.055641103266747e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b2k2000|d|1": {
   "rms_ro": 3.561443180723226e-06,
   "rms_roUx": 0.003080667605455897,
   "rms_roUy": 0.0003738482558269107,
   "rms_roUz": 0.0,
   "rms_roe": 6.38448525790721,
   "rms_roK": 0.0002897123690935012,
   "rms_roOmega": 3315.785139385651,
   "rms_roY0": 3.255789156902654e-06,
   "rms_roY1": 3.055641103266682e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b3k2000|f|1": {
   "rms_ro": 3.606417022007553e-06,
   "rms_roUx": 0.003102863486896804,
   "rms_roUy": 0.0003819059308616188,
   "rms_roUz": 0.0,
   "rms_roe": 6.445575039217898,
   "rms_roK": 0.0002896236169893312,
   "rms_roOmega": 3316.392468541205,
   "rms_roY0": 3.296907916003194e-06,
   "rms_roY1": 3.094232106667142e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b3k2000|d|1": {
   "rms_ro": 3.606417022007517e-06,
   "rms_roUx": 0.00310286348689687,
   "rms_roUy": 0.0003819059308615966,
   "rms_roUz": 0.0,
   "rms_roe": 6.444935810214566,
   "rms_roK": 0.0002896236169893315,
   "rms_roOmega": 3316.392468540863,
   "rms_roY0": 3.296907916003145e-06,
   "rms_roY1": 3.094232106667199e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "s0|f|2": {
   "rms_ro": 3.699704313945355e-06,
   "rms_roUx": 0.0031217156361968,
   "rms_roUy": 0.000386840168107559,
   "rms_roUz": 0.0,
   "rms_roe": 6.538800870347452,
   "rms_roK": 0.0002899758303757734,
   "rms_roOmega": 3251.612630613906,
   "rms_roY0": 3.382192595459116e-06,
   "rms_roY1": 3.17427395197152e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "s0|d|2": {
   "rms_ro": 3.699704313945307e-06,
   "rms_roUx": 0.003121715636196435,
   "rms_roUy": 0.0003868401681075276,
   "rms_roUz": 0.0,
   "rms_roe": 6.538648285966021,
   "rms_roK": 0.0002899758303757733,
   "rms_roOmega": 3251.612630613606,
   "rms_roY0": 3.382192595459142e-06,
   "rms_roY1": 3.174273951971615e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k2000|f|2": {
   "rms_ro": 3.674254260713308e-06,
   "rms_roUx": 0.003113306827588162,
   "rms_roUy": 0.000392401015600953,
   "rms_roUz": 0.0,
   "rms_roe": 6.512654161893741,
   "rms_roK": 0.0002871064075878133,
   "rms_roOmega": 3276.547711787721,
   "rms_roY0": 3.358925580652271e-06,
   "rms_roY1": 3.152437265569585e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "a1k2000|d|2": {
   "rms_ro": 3.674254260713288e-06,
   "rms_roUx": 0.003113306827588372,
   "rms_roUy": 0.0003924010156009615,
   "rms_roUz": 0.0,
   "rms_roe": 6.512400268350982,
   "rms_roK": 0.0002871064075878137,
   "rms_roOmega": 3276.547711787603,
   "rms_roY0": 3.3589255806523e-06,
   "rms_roY1": 3.152437265569699e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k2000|f|2": {
   "rms_ro": 3.643390393577298e-06,
   "rms_roUx": 0.003102712452654819,
   "rms_roUy": 0.0003937742139419386,
   "rms_roUz": 0.0,
   "rms_roe": 6.477839300828452,
   "rms_roK": 0.0002896141351855706,
   "rms_roOmega": 3320.033940488188,
   "rms_roY0": 3.330706552599263e-06,
   "rms_roY1": 3.125952988530732e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  },
  "b1k2000|d|2": {
   "rms_ro": 3.64339039357735e-06,
   "rms_roUx": 0.003102712452654885,
   "rms_roUy": 0.0003937742139418812,
   "rms_roUz": 0.0,
   "rms_roe": 6.47689227711699,
   "rms_roK": 0.0002896141351855707,
   "rms_roOmega": 3320.033940488143,
   "rms_roY0": 3.330706552599271e-06,
   "rms_roY1": 3.125952988530762e-07,
   "rms_dq_ro": 0.0,
   "rms_dq_roUx": 0.0,
   "rms_dq_roUy": 0.0,
   "rms_dq_roUz": 0.0,
   "rms_dq_roe": 0.0
  }
 },
 "sigma_eval": {
  "f": 1.8880648422584315e-14,
  "d": 2.0048249848895883e-14
 },
 "eval_effect_d_over_f": {
  "s0": {
   "rms_ro": 0.9999999999999902,
   "rms_roUx": 0.9999999999999766,
   "rms_roUy": 0.9999999999999322,
   "rms_roe": 0.999976664776254,
   "rms_roK": 0.9999999999999997,
   "rms_roOmega": 0.9999999999999554,
   "rms_roY0": 1.0000000000000029,
   "rms_roY1": 0.999999999999974
  },
  "a1k1000": {
   "rms_ro": 0.9999999999999865,
   "rms_roUx": 1.00000000000005,
   "rms_roUy": 1.0000000000001374,
   "rms_roe": 1.0000101721241006,
   "rms_roK": 1.0000000000000002,
   "rms_roOmega": 0.9999999999998636,
   "rms_roY0": 1.0000000000000027,
   "rms_roY1": 1.0000000000000078
  },
  "a1k1500": {
   "rms_ro": 1.0000000000000056,
   "rms_roUx": 0.9999999999999946,
   "rms_roUy": 1.0000000000000542,
   "rms_roe": 1.0000335127117215,
   "rms_roK": 1.0000000000000044,
   "rms_roOmega": 0.9999999999999958,
   "rms_roY0": 1.0000000000000115,
   "rms_roY1": 0.9999999999999978
  },
  "a1k2000": {
   "rms_ro": 1.0000000000000082,
   "rms_roUx": 1.0000000000000573,
   "rms_roUy": 0.9999999999999802,
   "rms_roe": 0.999961015350067,
   "rms_roK": 1.000000000000001,
   "rms_roOmega": 0.9999999999999634,
   "rms_roY0": 0.9999999999999993,
   "rms_roY1": 0.9999999999999806
  },
  "b1k1000": {
   "rms_ro": 0.999999999999997,
   "rms_roUx": 0.9999999999999444,
   "rms_roUy": 0.9999999999998873,
   "rms_roe": 0.999986295678838,
   "rms_roK": 1.000000000000002,
   "rms_roOmega": 1.000000000000003,
   "rms_roY0": 1.0000000000000064,
   "rms_roY1": 0.9999999999999889
  },
  "b1k1500": {
   "rms_ro": 1.000000000000016,
   "rms_roUx": 1.0000000000000504,
   "rms_roUy": 1.0000000000001148,
   "rms_roe": 0.9998920385186783,
   "rms_roK": 0.999999999999997,
   "rms_roOmega": 0.9999999999999559,
   "rms_roY0": 0.9999999999999998,
   "rms_roY1": 1.000000000000006
  },
  "b1k2000": {
   "rms_ro": 0.9999999999999795,
   "rms_roUx": 1.0000000000000708,
   "rms_roUy": 1.0000000000000555,
   "rms_roe": 0.9998538056182811,
   "rms_roK": 1.0000000000000009,
   "rms_roOmega": 1.0000000000000004,
   "rms_roY0": 0.9999999999999892,
   "rms_roY1": 0.9999999999999881
  },
  "a2k2000": {
   "rms_ro": 1.0000000000000056,
   "rms_roUx": 1.0000000000000513,
   "rms_roUy": 1.0000000000001195,
   "rms_roe": 0.9999492229922261,
   "rms_roK": 0.9999999999999976,
   "rms_roOmega": 1.0000000000000382,
   "rms_roY0": 1.0000000000000044,
   "rms_roY1": 1.0000000000000082
  },
  "a3k2000": {
   "rms_ro": 1.0000000000000135,
   "rms_roUx": 0.999999999999944,
   "rms_roUy": 0.999999999999985,
   "rms_roe": 0.9999642307941032,
   "rms_roK": 0.9999999999999983,
   "rms_roOmega": 0.9999999999999327,
   "rms_roY0": 0.9999999999999937,
   "rms_roY1": 1.000000000000001
  },
  "b2k2000": {
   "rms_ro": 1.0000000000000122,
   "rms_roUx": 1.000000000000015,
   "rms_roUy": 0.9999999999999858,
   "rms_roe": 0.9999764709789515,
   "rms_roK": 0.9999999999999994,
   "rms_roOmega": 1.000000000000026,
   "rms_roY0": 1.0000000000000184,
   "rms_roY1": 0.9999999999999787
  },
  "b3k2000": {
   "rms_ro": 0.99999999999999,
   "rms_roUx": 1.0000000000000213,
   "rms_roUy": 0.9999999999999419,
   "rms_roe": 0.999900826691266,
   "rms_roK": 1.000000000000001,
   "rms_roOmega": 0.9999999999998969,
   "rms_roY0": 0.9999999999999852,
   "rms_roY1": 1.0000000000000182
  }
 },
 "by_evaluator": {
  "f": {
   "A": [
    6.512654161893749,
    6.499090431533056,
    6.487879493511477
   ],
   "B": [
    6.477839300828546,
    6.384635482129856,
    6.445575039217898
   ],
   "mean_A": 6.499874695646095,
   "mean_B": 6.436016607392101,
   "ratio": 0.9901754893372376,
   "decrease": 0.009824510662762354,
   "w_A": 0.003811560921146116,
   "w_B": 0.01448160009277177,
   "noise": 0.01448160009277177,
   "u": 0.04344480027831531,
   "verdict": "NOT_SUPPORT",
   "all_cols_ratio": {
    "rms_ro": 0.9827104755720203,
    "rms_roUx": 0.9955907380519177,
    "rms_roUy": 0.9827077221775691,
    "rms_roe": 0.9901754893372376,
    "rms_roK": 1.009622773265842,
    "rms_roOmega": 1.0165325658295805,
    "rms_roY0": 0.9827091931149381,
    "rms_roY1": 0.9827091931147427
   }
  },
  "d": {
   "A": [
    6.512400268351112,
    6.498760427167691,
    6.487647427214041
   ],
   "B": [
    6.476892277117087,
    6.38448525790721,
    6.444935810214566
   ],
   "mean_A": 6.4996027075776155,
   "mean_B": 6.4354377817462876,
   "ratio": 0.9901278695455461,
   "decrease": 0.009872130454453898,
   "w_A": 0.0038083621800780424,
   "w_B": 0.014359088277099594,
   "noise": 0.014359088277099594,
   "u": 0.043077264831298784,
   "verdict": "NOT_SUPPORT",
   "all_cols_ratio": {
    "rms_ro": 0.9827104755720051,
    "rms_roUx": 0.9955907380519359,
    "rms_roUy": 0.9827077221775365,
    "rms_roe": 0.9901278695455461,
    "rms_roK": 1.0096227732658436,
    "rms_roOmega": 1.0165325658295767,
    "rms_roY0": 0.9827091931149364,
    "rms_roY1": 0.982709193114741
   }
  }
 },
 "f_d_concordant": true,
 "display": {
  "a1": {
   "tail_median": {
    "rms_ro": 3.6690747416727567e-06,
    "rms_roUx": 0.003124426454908325,
    "rms_roUy": 0.00039039660814514556,
    "rms_roe": 6.510580831571201,
    "rms_roK": 0.0002876444957581432,
    "rms_roOmega": 3242.9910822935162,
    "rms_roY0": 3.3541747409175627e-06,
    "rms_roY1": 3.147978481406638e-07
   },
   "step0": {
    "rms_ro": 3.699704313945369e-06,
    "rms_roUx": 0.003121715636196634,
    "rms_roUy": 0.0003868401681075896,
    "rms_roUz": 0.0,
    "rms_roe": 6.538800870347489,
    "rms_roK": 0.0002899758303757728,
    "rms_roOmega": 3251.612630613622,
    "rms_roY0": 3.382192595459102e-06,
    "rms_roY1": 3.174273951971518e-07,
    "rms_dq_ro": 0.0,
    "rms_dq_roUx": 0.0,
    "rms_dq_roUy": 0.0,
    "rms_dq_roUz": 0.0,
    "rms_dq_roe": 0.0
   },
   "drift": {
    "slope500": -0.0042474834957147625,
    "level_tail_vs_prev": -0.0011149984271937856,
    "osc_rel": 0.009355134016451229,
    "class": "許容内"
   },
   "segment_verdict": "NOT CONVERGED"
  },
  "b1": {
   "tail_median": {
    "rms_ro": 3.608519895553817e-06,
    "rms_roUx": 0.0031135299139064187,
    "rms_roUy": 0.0003820761302702006,
    "rms_roe": 6.4585530237904845,
    "rms_roK": 0.0002898589099890644,
    "rms_roOmega": 3262.1831909451967,
    "rms_roY0": 3.298815075452476e-06,
    "rms_roY1": 3.0960220244218096e-07
   },
   "step0": {
    "rms_ro": 3.699704313945344e-06,
    "rms_roUx": 0.003121715636196673,
    "rms_roUy": 0.0003868401681075867,
    "rms_roUz": 0.0,
    "rms_roe": 6.53864828596601,
    "rms_roK": 0.0002899758303757733,
    "rms_roOmega": 3251.612630613771,
    "rms_roY0": 3.382192595459139e-06,
    "rms_roY1": 3.174273951971508e-07,
    "rms_dq_ro": 0.0,
    "rms_dq_roUx": 0.0,
    "rms_dq_roUy": 0.0,
    "rms_dq_roUz": 0.0,
    "rms_dq_roe": 0.0
   },
   "drift": {
    "slope500": 0.007503492490443751,
    "level_tail_vs_prev": 0.001773268121701624,
    "osc_rel": 0.008373954047622206,
    "class": "許容内"
   },
   "segment_verdict": "NOT CONVERGED"
  },
  "a2": {
   "tail_median": {
    "rms_ro": 3.7011609836723416e-06,
    "rms_roUx": 0.0031312214004866725,
    "rms_roUy": 0.0003944853672363632,
    "rms_roe": 6.55352307290419,
    "rms_roK": 0.00028760775364448064,
    "rms_roOmega": 3243.551863139194,
    "rms_roY0": 3.383539491433139e-06,
    "rms_roY1": 3.175538048170959e-07
   },
   "step0": {
    "rms_ro": 3.699704313945341e-06,
    "rms_roUx": 0.003121715636196654,
    "rms_roUy": 0.0003868401681075663,
    "rms_roUz": 0.0,
    "rms_roe": 6.538800870347448,
    "rms_roK": 0.0002899758303757731,
    "rms_roOmega": 3251.612630613924,
    "rms_roY0": 3.38219259545906e-06,
    "rms_roY1": 3.174273951971506e-07,
    "rms_dq_ro": 0.0,
    "rms_dq_roUx": 0.0,
    "rms_dq_roUy": 0.0,
    "rms_dq_roUz": 0.0,
    "rms_dq_roe": 0.0
   },
   "drift": {
    "slope500": -0.0035374968553794914,
    "level_tail_vs_prev": 0.003221514883109167,
    "osc_rel": 0.010102050208697293,
    "class": "許容内"
   },
   "segment_verdict": "NOT CONVERGED"
  },
  "b2": {
   "tail_median": {
    "rms_ro": 3.5993377400110187e-06,
    "rms_roUx": 0.0031098127644034043,
    "rms_roUy": 0.00038477356511000587,
    "rms_roe": 6.4367300198328135,
    "rms_roK": 0.0002898735348027722,
    "rms_roOmega": 3262.8818747312694,
    "rms_roY0": 3.2904513347177282e-06,
    "rms_roY1": 3.088172440572918e-07
   },
   "step0": {
    "rms_ro": 3.699704313945342e-06,
    "rms_roUx": 0.003121715636196844,
    "rms_roUy": 0.0003868401681075219,
    "rms_roUz": 0.0,
    "rms_roe": 6.538648285966071,
    "rms_roK": 0.0002899758303757732,
    "rms_roOmega": 3251.612630613714,
    "rms_roY0": 3.382192595459121e-06,
    "rms_roY1": 3.174273951971554e-07,
    "rms_dq_ro": 0.0,
    "rms_dq_roUx": 0.0,
    "rms_dq_roUy": 0.0,
    "rms_dq_roUz": 0.0,
    "rms_dq_roe": 0.0
   },
   "drift": {
    "slope500": -0.007293884640044055,
    "level_tail_vs_prev": 0.0006239438567942884,
    "osc_rel": 0.008501549390896944,
    "class": "許容内"
   },
   "segment_verdict": "NOT CONVERGED"
  },
  "a3": {
   "tail_median": {
    "rms_ro": 3.6898812521424957e-06,
    "rms_roUx": 0.003132324730023873,
    "rms_roUy": 0.000396492658662321,
    "rms_roe": 6.544462431357358,
    "rms_roK": 0.0002874449031480455,
    "rms_roOmega": 3241.853540742747,
    "rms_roY0": 3.3732206900899795e-06,
    "rms_roY1": 3.1658535901201596e-07
   },
   "step0": {
    "rms_ro": 3.699704313945389e-06,
    "rms_roUx": 0.003121715636196416,
    "rms_roUy": 0.0003868401681075653,
    "rms_roUz": 0.0,
    "rms_roe": 6.538800870347385,
    "rms_roK": 0.0002899758303757726,
    "rms_roOmega": 3251.612630613772,
    "rms_roY0": 3.382192595459088e-06,
    "rms_roY1": 3.174273951971519e-07,
    "rms_dq_ro": 0.0,
    "rms_dq_roUx": 0.0,
    "rms_dq_roUy": 0.0,
    "rms_dq_roUz": 0.0,
    "rms_dq_roe": 0.0
   },
   "drift": {
    "slope500": -0.0033009628212249263,
    "level_tail_vs_prev": 0.002373503791957498,
    "osc_rel": 0.009428386397529302,
    "class": "許容内"
   },
   "segment_verdict": "NOT CONVERGED"
  },
  "b3": {
   "tail_median": {
    "rms_ro": 3.6153348402757703e-06,
    "rms_roUx": 0.0031138186827136677,
    "rms_roUy": 0.00038642467718558166,
    "rms_roe": 6.45617069506297,
    "rms_roK": 0.00028978987380676527,
    "rms_roOmega": 3262.2903321500253,
    "rms_roY0": 3.305068709265736e-06,
    "rms_roY1": 3.1018912191398485e-07
   },
   "step0": {
    "rms_ro": 3.699704313945329e-06,
    "rms_roUx": 0.003121715636196766,
    "rms_roUy": 0.0003868401681075565,
    "rms_roUz": 0.0,
    "rms_roe": 6.538648285966008,
    "rms_roK": 0.0002899758303757728,
    "rms_roOmega": 3251.612630613823,
    "rms_roY0": 3.38219259545911e-06,
    "rms_roY1": 3.17427395197153e-07,
    "rms_dq_ro": 0.0,
    "rms_dq_roUx": 0.0,
    "rms_dq_roUy": 0.0,
    "rms_dq_roUz": 0.0,
    "rms_dq_roe": 0.0
   },
   "drift": {
    "slope500": -0.002186187349652276,
    "level_tail_vs_prev": 0.0015769558813011177,
    "osc_rel": 0.009352723094224346,
    "class": "許容内"
   },
   "segment_verdict": "NOT CONVERGED"
  }
 },
 "display_ratio_B_over_A": 0.9868877411666435,
 "eval_effect_at_A2000": 0.9999581563787988,
 "drift": {
  "a": {
   "classes": [
    "許容内",
    "許容内",
    "許容内"
   ],
   "class": "許容内",
   "text": "この窓で有意なドリフトを検出せず"
  },
  "b": {
   "classes": [
    "許容内",
    "許容内",
    "許容内"
   ],
   "class": "許容内",
   "text": "この窓で有意なドリフトを検出せず"
  }
 },
 "series_vs_S": {
  "a1": {
   "f": {
    "1000": 1.0147785633313193,
    "1500": 0.9965665912953516,
    "2000": 0.9960012991721231
   },
   "d": {
    "1000": 1.0148125666631034,
    "1500": 0.996623245370658,
    "2000": 0.9959857119595784
   }
  },
  "b1": {
   "f": {
    "1000": 0.9991708945158523,
    "1500": 0.9739521070852925,
    "2000": 0.9906769496842708
   },
   "d": {
    "1000": 0.9991805176579603,
    "1500": 0.973869683240033,
    "2000": 0.9905552331081258
   }
  }
 },
 "VERDICT": "NOT_SUPPORT: この期間の主要因説を支持しない; ドリフト: A この窓で有意なドリフトを検出せず・B この窓で有意なドリフトを検出せず; f と d の判定の一致: True"
}
```

## 参考: `case/45.isobutane_m6_d155/fh_floor_judge.py`

```
#!/usr/bin/env python3
"""point 仕上げでの面エンタルピーの精度の A/B の判定 (plan time_integration-implicit-thermal-jacobian §6.2、2026-10-10 事前登録)。

fh_floor.sh の run (軌道 run_0500〜0505、評価 run_0510〜0537) を読み、§6.2 の規則どおりに判定して
_band_ab/cold_pair/fh_floor_judge.json に書く。規則を変えるときは plan §6.2 を先に改訂する (結果を見てから変えない)。

  主判定の量: 評価器 d (面エンタルピーを double で評価) の rms_roe を、A・B の 2000 step の状態で比べる。
  E_X(状態) = 評価 run の step 0 の outer_begin の行 (更新前の残差)。X = f (既定) / d (切替あり)。
  ゲートが 1 つでも外れたら量を判定せず INVALID を出す。
"""
import csv
import hashlib
import json
import math
import re
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
# 入力の同一性を run_0354 と比べるファイル (prep が run_0183 から複製するもの + 解決済みの化学種)
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv",
              "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
CONS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1", "P", "T"]
# 事前登録の閾値 (§6.2)
DROP = 0.10          # 検出の対象にする低下 (10 %)
NOISE_MULT = 3.0     # 不確かさの幅 u = 3 × noise (統計的な 3σ ではない。探索試験の目安)
DRIFT = 0.01         # 傾き × 500 と窓の水準差の許容 [桁] (0.01 桁 ≈ 2.3 %)
GATE_SELF = 1e-6     # 軌道の step 0 の表示と、同じ状態・同じ評価器の評価の相対差の上限
TAIL, PREV = (1500, 1999), (1000, 1499)


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


def read_csv(run):
    """全行 (phase ごと) と outer_begin の行 (step → 値)、行の重複・非有限・負を返す。"""
    allrows, outer, dup, bad = [], {}, [], []
    with open(HERE / run / "residual_history.csv") as f:
        rd = csv.DictReader(f)
        rcols = [c for c in rd.fieldnames if c.startswith("rms_")]
        for r in rd:
            vals = {c: float(r[c]) for c in rcols}
            if any((not math.isfinite(v)) or v < 0 for v in vals.values()):
                bad.append((r["step"], r["phase"]))
            allrows.append(r)
            if r["phase"] == "outer_begin":
                s = int(r["step"])
                if s in outer:
                    dup.append(s)
                outer[s] = vals
    return outer, dup, bad, len(allrows)


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


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def input_diff(run, ref):
    out = [fn for fn in SAME_FILES if (HERE / ref / fn).is_file() != (HERE / run / fn).is_file()
           or ((HERE / ref / fn).is_file() and sha(HERE / ref / fn) != sha(HERE / run / fn))]
    rs = lambda d: {p.name: sha(p) for p in sorted((HERE / d).glob("resolved_species_*.yaml"))}
    if rs(run) != rs(ref):
        out.append("resolved_species_*.yaml")
    return out


def switch_on(run):
    return "FORGE_DIAG_FACE_H_DOUBLE" in (HERE / run / "forge_run.log").read_text(errors="replace")


def rc_ok(run):
    p = HERE / run / "RUN_RC"
    return p.is_file() and p.read_text().strip() == "0"


def rel(a, b):
    return abs(a - b) / (0.5 * (abs(a) + abs(b)))


def seg_verdict(run):
    p = HERE / run / "CONVERGENCE_SEGMENT.txt"
    if not p.is_file():
        return "判定不能 (CONVERGENCE_SEGMENT.txt が無い)"
    m = re.search(r"->\s*([A-Z][A-Z ]+[A-Z])", p.read_text(errors="replace"))
    return m.group(1) if m else "判定不能 (VERDICT の行が無い)"


def drift(rows):
    """末尾 500 step の log10(rms_roe) の傾き × 500 と、末尾の窓とその前の窓の中央値の差 [桁]、振動の大きさ。"""
    s = np.arange(TAIL[0], TAIL[1] + 1, dtype=float)
    y = np.log10([rows[int(k)][Q] for k in s])
    c = np.polyfit(s, y, 1)
    slope = float(c[0] * 500.0)
    lvl = float(np.log10(np.median([rows[k][Q] for k in range(*TAIL)] + [rows[TAIL[1]][Q]]) /
                         np.median([rows[k][Q] for k in range(*PREV)] + [rows[PREV[1]][Q]])))
    osc = float(np.std(y - np.polyval(c, s)) * math.log(10))   # 傾きを除いた残りの相対の揺れ (step ごと)
    if slope <= -DRIFT and lvl <= -DRIFT:
        cls = "減衰"
    elif slope >= DRIFT and lvl >= DRIFT:
        cls = "増加"
    elif abs(slope) < DRIFT and abs(lvl) < DRIFT:
        cls = "許容内"
    else:
        cls = "混在"
    return {"slope500": slope, "level_tail_vs_prev": lvl, "osc_rel": osc, "class": cls}


def main():
    gates, rec = [], {"plan": "time_integration-implicit-thermal-jacobian §6.2", "quantity": Q}
    ev = eval_runs()
    st = HERE / "_fh_states"
    mesh_ref = (st / "MESH_SHA_REF.txt").read_text().strip()
    s0sha = (st / "s0" / "SHA256").read_text().strip()
    gates.append(("S の sha256 が run_0354/res_40000.h5 と一致", s0sha == sha(HERE / S0 / S0_RES), s0sha[:16]))

    def common(run, state_name):
        gates.append((f"{run}: RUN_RC = 0", rc_ok(run), (HERE / run / "RUN_RC").read_text().strip() if (HERE / run / "RUN_RC").is_file() else None))
        d = cfg_diff(run, S0)
        gates.append((f"{run}: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)", d == [], d))
        d = input_diff(run, S0)
        gates.append((f"{run}: 境界条件・化学種・壁などの入力が run_0354 と同じ", d == [], d))
        m = (HERE / run / "MESH_SHA.txt").read_text().strip() if (HERE / run / "MESH_SHA.txt").is_file() else ""
        gates.append((f"{run}: 格子 (MESH) が run_0354 と同じ", m == mesh_ref, m[:16]))
        cp = json.loads((HERE / run / "COLD_PAIR.json").read_text())
        want = (st / state_name / "SHA256").read_text().strip()
        gates.append((f"{run}: 出発の場が状態 {state_name}", cp["field_from"] == state_name and cp["parent_res_sha256"] == want,
                      [cp["field_from"], cp["parent_res"], cp["parent_res_sha256"][:16]]))

    # --- ゲート: 軌道 ---
    traj_rows = {}
    for name, (run, fh) in TRAJ.items():
        common(run, "s0")
        rows, dup, bad, n = read_csv(run)
        traj_rows[name] = rows
        gates.append((f"{run}: outer_begin の行が step 0〜1999 で一意・連続", sorted(rows) == list(range(2000)) and not dup,
                      [min(rows, default=None), max(rows, default=None), len(rows), dup[:3]]))
        gates.append((f"{run}: 全行の rms_* が有限・非負", not bad, bad[:3]))
        gates.append((f"{run}: 切替の表示 = {fh}", switch_on(run) == bool(fh), switch_on(run)))
        with h5py.File(HERE / run / "res_2000.h5", "r") as h:
            nonfin = {v: int(np.sum(~np.isfinite(h["VALUE"][v][...]))) for v in CONS}
            nonpos = {v: int(np.sum(h["VALUE"][v][...] <= 0)) for v in ("ro", "P", "T")}
        gates.append((f"{run}: res_2000 の保存量・P・T が有限", sum(nonfin.values()) == 0, nonfin))
        gates.append((f"{run}: res_2000 の ρ・P・T が正", sum(nonpos.values()) == 0, nonpos))
    # --- ゲート: 評価 ---
    E = {}
    for (sd, x, n), run in ev.items():
        common(run, sd)
        rows, dup, bad, nrow = read_csv(run)
        gates.append((f"{run}: step 0 の outer_begin が 1 行", 0 in rows and not dup, [sorted(rows)[:3], dup]))
        gates.append((f"{run}: 全行の rms_* が有限・非負", not bad, bad[:3]))
        gates.append((f"{run}: 切替の表示 = {x}", switch_on(run) == (x == "d"), switch_on(run)))
        e = rows.get(0)
        gates.append((f"{run}: 判定に使う値が正", e is not None and all(e[c] > 0 for c in COLS if c != "rms_roUz"), e and e[Q]))
        E[(sd, x, n)] = e
    ok = all(g[1] for g in gates)
    if ok:   # 自己一致 (保守的な停止の基準。外れても「別の残差を読んだ」とは断定しない)
        self_f = rel(traj_rows["a1"][0][Q], E[("s0", "f", 1)][Q])
        self_d = rel(traj_rows["b1"][0][Q], E[("s0", "d", 1)][Q])
        gates.append(("A1 の step 0 の表示と E_f(S) の相対差 ≤ 1e-6", self_f <= GATE_SELF, self_f))
        gates.append(("B1 の step 0 の表示と E_d(S) の相対差 ≤ 1e-6", self_d <= GATE_SELF, self_d))
        ok = all(g[1] for g in gates)
    rec["gates"] = [{"check": g[0], "ok": bool(g[1]), "value": g[2]} for g in gates]
    rec["gates_ok"] = bool(ok)
    if not ok:
        rec["VERDICT"] = "INVALID (ゲート不合格: 量を判定しない)"
        (HERE / "_band_ab" / "cold_pair" / "fh_floor_judge.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
        for g in gates:
            if not g[1]:
                print("ゲート不合格:", g[0], g[2])
        print(rec["VERDICT"])
        return 1
    # --- 量 ---
    rec["E"] = {f"{sd}|{x}|{n}": E[(sd, x, n)] for (sd, x, n) in E}
    sig = {x: max(rel(E[(sd, x, 1)][Q], E[(sd, x, 2)][Q]) for sd in REPEAT) for x in ("f", "d")}
    rec["sigma_eval"] = sig
    rec["eval_effect_d_over_f"] = {sd: {c: E[(sd, "d", 1)][c] / E[(sd, "f", 1)][c] for c in COLS if E[(sd, "f", 1)][c] > 0} for sd in STATES}
    per = {}
    for x in ("f", "d"):
        A = [E[(f"a{i}k2000", x, 1)][Q] for i in (1, 2, 3)]
        Bv = [E[(f"b{i}k2000", x, 1)][Q] for i in (1, 2, 3)]
        mA, mB = float(np.mean(A)), float(np.mean(Bv))
        wA, wB = (max(A) - min(A)) / mA, (max(Bv) - min(Bv)) / mB
        noise = max(sig[x], wA, wB)
        u = NOISE_MULT * noise
        r = mB / mA
        dec = 1.0 - r
        if dec - u >= DROP and max(Bv) <= (1 - DROP) * min(A):
            v = "SUPPORT"
        elif dec + u < DROP:
            v = "NOT_SUPPORT"
        else:
            v = "INDETERMINATE"
        per[x] = {"A": A, "B": Bv, "mean_A": mA, "mean_B": mB, "ratio": r, "decrease": dec, "w_A": wA, "w_B": wB,
                  "noise": noise, "u": u, "verdict": v,
                  "all_cols_ratio": {c: float(np.mean([E[(f"b{i}k2000", x, 1)][c] for i in (1, 2, 3)]) /
                                         np.mean([E[(f"a{i}k2000", x, 1)][c] for i in (1, 2, 3)]))
                                     for c in COLS if all(E[(f"a{i}k2000", x, 1)][c] > 0 for i in (1, 2, 3))}}
    rec["by_evaluator"] = per
    rec["f_d_concordant"] = per["f"]["verdict"] == per["d"]["verdict"]
    # 自方式の表示 (各腕の CSV) と、反復ごとのドリフト
    disp = {}
    for name, rows in traj_rows.items():
        disp[name] = {"tail_median": {c: float(np.median([rows[k][c] for k in range(TAIL[0], TAIL[1] + 1)])) for c in COLS},
                      "step0": rows[0], "drift": drift(rows), "segment_verdict": seg_verdict(TRAJ[name][0])}
    rec["display"] = disp
    dA = float(np.mean([disp[f"a{i}"]["tail_median"][Q] for i in (1, 2, 3)]))
    dB = float(np.mean([disp[f"b{i}"]["tail_median"][Q] for i in (1, 2, 3)]))
    rec["display_ratio_B_over_A"] = dB / dA
    ee_A = float(np.mean([E[(f"a{i}k2000", "d", 1)][Q] / E[(f"a{i}k2000", "f", 1)][Q] for i in (1, 2, 3)]))
    rec["eval_effect_at_A2000"] = ee_A
    floor = {}
    for arm in ("a", "b"):
        cl = [disp[f"{arm}{i}"]["drift"]["class"] for i in (1, 2, 3)]
        c = cl[0] if len(set(cl)) == 1 else "反復間不一致"
        floor[arm] = {"classes": cl, "class": c,
                      "text": {"許容内": "この窓で有意なドリフトを検出せず", "減衰": "減衰中 = 床の判定は不能 (自動で延長しない)",
                               "増加": "増加 (頭打ちではない)"}.get(c, "判定不能 (" + c + ")")}
    rec["drift"] = floor
    rec["series_vs_S"] = {f"{arm}1": {x: {k: E[(f"{arm}1k{k}", x, 1)][Q] / E[("s0", x, 1)][Q] for k in (1000, 1500, 2000)}
                                      for x in ("f", "d")} for arm in ("a", "b")}
    # 主判定 (評価器 d)
    v = per["d"]["verdict"]
    txt = {"SUPPORT": "この期間の停滞への寄与を支持 (d で評価した離散残差の改善に限る)",
           "NOT_SUPPORT": "この期間の主要因説を支持しない",
           "INDETERMINATE": "判別不能 (10 % の境界を不確かさの幅がまたぐ、または範囲が分かれない)"}[v]
    if v != "SUPPORT" and rec["display_ratio_B_over_A"] <= 1 - DROP:
        u = per["d"]["u"]
        only = abs(1 - per["d"]["ratio"]) <= u and abs(rec["display_ratio_B_over_A"] - ee_A) <= u
        txt += ("; 自方式の表示では 10 % 以上下がったが共通の評価 (d) では基準に届かない" +
                (" — 同じ状態の f/d の差で表示の差を説明できるので、表示の定義の違いだけ" if only else ""))
    rec["VERDICT"] = f"{v}: {txt}; ドリフト: A {floor['a']['text']}・B {floor['b']['text']}; f と d の判定の一致: {rec['f_d_concordant']}"
    out = HERE / "_band_ab" / "cold_pair" / "fh_floor_judge.json"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=float))
    print(rec["VERDICT"])
    print("sigma_eval", sig, "表示の比 B/A", round(rec["display_ratio_B_over_A"], 5), "A の 2000 での d/f", round(ee_A, 5))
    for x in ("f", "d"):
        p = per[x]
        print(f"評価器 {x}: A {['%.5g' % a for a in p['A']]} B {['%.5g' % b for b in p['B']]} 比 {p['ratio']:.4f} u {p['u']:.3e} → {p['verdict']}")
    for name in TRAJ:
        dd = disp[name]["drift"]
        print(f"{name}: 傾き×500 {dd['slope500']:+.4f} 窓の差 {dd['level_tail_vs_prev']:+.4f} 揺れ {dd['osc_rel']:.2e} {dd['class']} / {disp[name]['segment_verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
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
