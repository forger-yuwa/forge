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

## ブリーフ (`notes/reviews/briefs/2026-10-08-linedir-divergence.md`)

# 諮問: 方向別の擬似 dt (`lineDtDirectional`) が冷却壁の格子で発散する理由

日付 2026-10-08。諮問先 codex (diagnose)。エスカレーション条件 2 (同じ対処で 2 回発散) と 4 (原因を書く前)。
ユーザ: 「方向別の刻みというのがなぜうまくいかないのか、すまんが追求したい」。
plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #27 (全文)、`plans/accepted/time_integration-line-implicit-viscous-v2.md` (2026-09-03 の case/45 の記録)。
作業ツリー `/home/sano/work/forge-integ-1005` (commit 3b8966f1)。forge は FP64 のビルド (`~/forge-wallfit-bin-fp64`、AWS)。

## 観測事実

共通の設定 (300 K 等温壁の冷却ノズル case/45、run_0183 の res_100000 からビット一致で restart、リミッタの基準値は run_0183 の値に固定):
node・軸対称・`nodeWallDirichlet 1`・SLAU (slauWallNormalChi 1 自動)・2 次 (convMethod 1、limiter 2)・SST (dilatationCorrection 2、katoLaunder 1)・
陰解法 block-DPLUR (timeIntegration 11)・`nStepInner 5`・`implicitRelax 0.7`・`lowMachPrecond 0`・unsteady 0・燃焼ガス TP。
格子は ni 4719 × nj 121 の構造 (node = i·121 + j、j = 120 が壁)、近壁は壁法線、第一層 y1/局所半径 ≈ 3.4e-7、近壁の縦横比 (流れ方向の幅 / 壁法線の幅) は縮流部〜スロートで 3500〜4400。

| run (`case/45.isobutane_m6_d155/` AWS) | 設定 | 結果 |
| --- | --- | --- |
| run_0191 | point (lineImplicit 0)、cfl 4 | 40000 step 安定 (NaN なし)。rms_roOmega は上下するが有界 |
| run_0200 | `lineImplicit 1` + `lineDtDirectional 1`、cfl 4 | 65 step で非有限 (detectNaN が 66 で停止) |
| run_0202 | 同、**cfl 1** | 129 step で非有限 |
| run_0203 | run_0200 と同じ設定の再実行 + `FORGE_DUMP_LEDGER` (列 36〜52 × 壁から 0〜30 層、527 節点、毎 step の状態) | 66 step で非有限 (run_0200 と同じ = 決定的) |
| run_0201 | `lineImplicit 1` だけ (方向別 dt なし)、cfl 4 | 走行中・安定 (1 step のコストは point の約 2 倍) |

- ライン: `[lineImplicit] lines=4719 covered CVs=570999/570999 (100.0%) maxLen=121` (1 列 = 1 ライン、壁から軸まで)。
- 壊れた場所 (run_0200・0202 とも): x/r_t −10.8〜−10.2 (列 39〜49) のライン 11 本で、ρ が壁から軸まで全層で非有限 (ライン解が 1 本丸ごと壊れる)。
- 帳簿ダンプ (run_0203) の振れ始め:
  - step 2 で最大の相対変化は ρ 4e-4 (列 36、壁から 6 層目)。
  - step 3〜8 で最大は**壁から 1 層目** (壁節点のすぐ内側) の列 37〜44: Uy (半径方向速度、|Δ| / |Ux|) 1.6e-2 → 8e-2、ρ・T・ρω が 1e-3 → 1e-2。
  - step 8 の ρ の変化の上位 12 点は全て列 37〜44 の 1〜2 層目。
  - その後は**周期およそ 16 step の、振幅が育つ低周波の振動** (1 step ごとの符号反転ではない)。列 44・5 層目の ρ の増分の符号は `++++++++-------++++++++--------+++++++++------++++++++++------++`、
    ρ は (4 step ごと) 31.07, 31.19, 31.35, 31.13, 30.97, 31.78, 32.37, 30.30, 29.82, 36.9, 47.4, 27.7, 25.2, 50.0, 59.9, 21.3, 24.0。
  - step 24 以降は Ux・Uy の 1 層目の変化が O(1) を超え、step 56〜65 で 12〜14 層目まで広がる。
- 形状 (壊れた列に特有なもの):

| 列 | x/r_t (壁) | 壁の傾き | 流れ方向の格子間隔 (r_t) | 1 層目の縦横比 | 壁から 10 層目の M | 中心の M |
| --- | --- | --- | --- | --- | --- | --- |
| 5〜20 | −12.3〜−11.8 | 0° | 0.037 | 3500〜4400 | 0.007 | 0.023 |
| 35〜49 | −11.0〜−10.2 | −1.7〜−5.5° (曲がり始め) | **0.061 (領域で最大)** | 3900〜4000 | 0.005〜0.006 | 0.03 |
| 120〜200 | −7.3〜−5.4 | −26〜−36° | 0.023〜0.034 | 4000 | 0.005〜0.009 | 0.03 |
| 1800 | 0.03 (スロート) | 1° | 0.0016 | 4300 | 0.35 | 0.92 |

- 実装 (`solver_density_cuda/cuda_forge/setDT_d.cu` 133〜180 行): 節点の cfl = 面の cfl_pln の**最大**。`lineDtDirectional` では `line_prev/next` に一致する内部面を最大から外す。
  壁節点の境界半割面は外さない (壁節点の Δτ は壁法線の音響で縛られたまま)。1 層目以降の Δτ は流れ方向の面だけで決まる (おおよそ Δx/(|u|+c)、縦横比の分だけ伸びる)。
- 過去 (2026-09-03、`plans/accepted/time_integration-line-implicit-viscous-v2.md`): case/45 の断熱・y+≈2 の格子 (別格子、FP32) で、directional は cfl 6・8 とも step 25 で発散。
  種は x/r_t 71〜89 の下流域・全断面 (壁至近 3 %・近軸 1 %) で、「既知の streamwise 内部モードが Δτ の上がった下流域で先に点火した」と記録。
  同じ格子で line だけ (directional なし) は 1 step 1.81 倍で末尾 roe −11 %、ni2 (nStepInner 2) は発散 (off-line lag に sweep ≥ 3 が要る)。
  case/39 (DDES、dual-time) では directional は cp4〜cp8 で安定 — dual-time の BDF の対角が保護していた、と推測されている (未確認)。

## 期待値

- 定常解は Δτ の取り方によらない。方向別 dt は、ラインで厳密に解く方向の λ を Δτ の制約から外すので、壁近くの遅い過渡 (縮流部の壁から 1〜60 層目に残る連続の残差、低 M × 高縦横比) を速めるはずだった。

## 実施済みの操作

- 上表の 4 本。CFL を 4 → 1 に下げても同じ場所で発散 (step 65 → 129)。
- point の cfl 8 (run_0199) は発散しないが、スロートの壁から 7 層目で残差が 60〜120 倍に張り付いた (別の現象)。

## 仮説 (確かめていない)

1. **off-line (流れ方向) の lag の不安定**: Δτ が縦横比倍に伸びると対角の V/Δτ が小さくなり、DPLUR の off-line の Jacobi 的な緩和 (nStepInner 5、relax 0.7) が低周波のモードを増幅する。流れ方向の格子間隔が最大の列で Δτ が最大になるので、そこで先に点火する。
2. **壁節点と 1 層目の Δτ の不整合**: ラインの解で、Δτ の小さい壁節点 (Dirichlet) と Δτ が数千倍の 1 層目が結合している。1 層目の更新が壁節点の拘束と食い違う。
3. **LHS と RHS の不整合 (defect correction)**: LHS は 1 次の Jacobian、RHS は 2 次の SLAU (低 M では SLAU の圧力の散逸が LHS と違う)。V/Δτ が小さいと不整合が増幅される。壁近くは M ≈ 0.005。
4. 2026-09-03 の下流の発散と同じ型 (定常では Δτ を伸ばすと保護が無い) が、この格子では縮流部の近壁で先に出た。

## 問い

1. 観測 (1 層目から、Uy 主体、周期 16 step の成長、流れ方向の格子間隔が最大の列、CFL によらない) と最も整合する機構は何か。上の仮説のどれか、または別のものか。
2. 機構を切り分ける最小の A/B は何か。既存のスイッチ (nStepInner、implicitRelax、1 次の RHS (convMethod 0)、limiter、slauWallNormalChi、FORGE_DUMP_LEDGER の faces) だけでできるものを優先し、コードの変更 (例: 方向別の Δτ の伸びに上限 R を付ける、壁の 1 層目だけ除外しない) が要るものは分けて示す。
3. 収束を速める目的にとって、この方向 (directional dt の安定化) を追う価値はあるか。line だけ (1 step 約 2 倍) で得をするには、1 step あたり 2 倍以上の収束の速さが要る。

## 読んでよいもの

- 上記 2 つの plan、`solver_density_cuda/cuda_forge/setDT_d.cu`、`timeIntegration_d.cu` (line の Thomas と DPLUR の対角・sweep)、`methods/time_integration/implementation.md` の「line-implicit」節、`procedures/solver-settings.md` の「lineImplicit」
- `case/45.isobutane_m6_d155/README.md` (run_0183・0190〜0203 の行)

## 関連 plan 全文 (`plans/active/tooling-nozzle-isothermal-wall-chain.md`)

```markdown
# ノズル設計チェーンの等温壁化: 冷却壁 (T_w 指定) の低 Re SST 検証 (超音速平板 → ノズル × SU2) + 等温 NS の δ\* 反復 + 壁温影響の評価方針

## メタ

- **area**: `tooling / boundary layer / boundary`
- **status**: `in_progress`  <!-- 2026-09-12 起票。S0 (plan) → S1 配管 → S2 平板 → S3 ノズル CPG × SU2 → S4 生産 TP 等温 δ* 反復 → S5 方針 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「壁の熱境界条件 (断熱 / 等温) と壁温影響の評価」節 (本計画と同時に起草)
  - [`methods/boundary.md`](../../methods/boundary.md) `wall_isothermal` (ゴースト構成・node 壁ノード温度ピン)
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5 (壁処理。本計画は low-Re `wallTreatmentSST: 0` のみ)
- **related_plans**:
  - 親: [`../accepted/tooling-nozzle-deltastar-core-matched-euler.md`](../accepted/tooling-nozzle-deltastar-core-matched-euler.md) (δ\* 反復の生産形。§9 の残件「等温壁の NS 実行 (V3 の CFD 部分)」を本計画が引き受ける)
  - [`../active/turbulence-sst-thermal-flux-model.md`](turbulence-sst-thermal-flux-model.md) (壁関数 × 等温壁の Kader q_w。**本計画では使わない** — 圧縮性冷却壁で +87 % 過大の既知限界、§3)
  - [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) (⑤ SERN。等温化は共通 bcond 配管で自動追随、検証は本計画の対象外)
  - [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md) (①② 風洞チェーン = 本計画の主対象)
- **created**: `2026-09-12`
- **owner**: `sano`

## 1. 目的

ノズル設計チェーンの NS 評価 (風洞 ①②④ の δ\* 反復、⑤ SERN の力評価) は現状**断熱壁**で回している。実機の壁は冷却 (水冷銅 / 再生冷却) か
有限熱容量の壁で、設計上は「等温壁 ($T_w$ 指定)」か「CHT」のどちらかになる。本計画で得る状態:

1. **等温壁の低 Re SST (壁関数なし) が、冷却された超音速乱流境界層の厚さ ($\delta^*, \theta, H$)・摩擦・熱流束を理論/経験式と SU2 の両方に対して再現する**ことを、ノズル出口相当の条件の平板で確認済み ($T_w/T_{aw}$ 0.26 と 1.0、y⁺ 掃引つき)。
2. ノズル形状で等温 NS が回り、**同一メッシュ・同一 BC の SU2 と境界層厚さ・壁熱流束が一致**する (CPG チェーン)。
3. **等温 NS で δ\* 反復が回る** (積分法初期壁 → NS → 帯局所抽出 → 壁更新、全段が同じ $T_w$ を見る) — 生産 TP 風洞 (case/44 va3, $T_w$=300 K) でゲート達成。
4. **最適設計で壁温影響をどう評価するか**の方針 (§4.6) が plan に確定し、チェーンの台帳に「壁温感度」が標準出力として入る。

## 2. スコープ

- **やる**
  - 問題定義 YAML に壁熱境界条件 `spec.wall_thermal` を追加し、bcond (`wall_isothermal` + `Ts`)・積分法初期壁 (`thermal_bc`)・帳簿の 3 箇所を**単一のソースから**駆動する。
  - 超音速冷却壁平板の検証ケース (新規 `case/48.flat_plate_cooled_m4`): forge node 低 Re SST (断熱 / $T_w$ 300 K / 中間) + y₁⁺ 掃引 + SU2 同一メッシュ比較 + 理論/経験式との照合ツール。
  - ノズル CPG チェーン (case/45 `problem_d155_cpg_ns.yaml` 系) で等温 NS × SU2 等温の同一メッシュ比較。
  - 生産 TP 風洞 (case/44 va3 M4.19) で等温 δ\* 反復 (pass 0 積分法 + pass 1) と壁温感度台帳。
  - `methods/design/overview.md`・`procedures/recommended-settings.md`・`procedures/verification/` の同期。
- **やらない**
  - **壁関数 (`wallTreatmentSST: 1`) の等温壁** — ユーザ判断 (2026-09-12): 壁関数の準備状況が悪く、まず壁関数なしで確立する。Kader q_w の圧縮性補正は [`turbulence-sst-thermal-flux-model.md`](turbulence-sst-thermal-flux-model.md) §8 のまま。
  - **CHT (固体伝導との連成)** — 本計画は「等温 = $T_w$ 既知」まで。$T_w(x)$ が未知で強く連成する場合の「弱 CHT ループ」(1D 壁伝導/冷却モデルで $q_w \to T_w(x)$ を反復) は §4.6 で方針だけ決め、別 plan で実装する。
  - $T_w(x)$ 分布の NS 入力 (`Tw_table`) — 積分法初期壁は既対応だが forge の `wall_isothermal` は bcond 単位の定数 `Ts` しか持たない。初版は定数 $T_w$。分布は `inletProfile` と同型の `wallProfile` CSV として別途 (§5.1)。
  - SERN (⑤) の等温評価の検証 — bcond 配管の共通化で `runner_sern.py` も同じ `wall_thermal` を読めるようにするが、検証 run は SERN の R1–R7 の後。
  - 化学非平衡・凝縮との組合せ (等温壁 × 凝縮 ON は case/44 の凝縮 restart 手順で後続)。

## 3. 関連 docs と前提 (既存資産の照合)

| 部品 | 既存資産 | 状態・ギャップ |
| --- | --- | --- |
| 等温壁 BC | `wall_isothermal` (floats `Ts`, `Ux/Uy/Uz`)。cell: ゴースト $T_R = 2T_w - T_L$ (2026-07-20 修正、純伝導厳密解 +0.02 %)。node: 壁ノード温度ピン + `res_roe` 0 化 + DPLUR 行 decouple ([`methods/boundary.md`](../../methods/boundary.md) 「node 等温壁の壁ノード温度ピン」) | **再利用**。case/40 `run_0048` (node y⁺1, $T_s$ 1000 K) で運用実績。y₁⁺≲1 の step1 発散は cross-mesh IC の P 段差が真因で解決済 ([[node-isothermal-wall-thin-cell-mass-source]]) — **同一メッシュ index コピー warm start なら問題なし** |
| 低 Re SST | `wallTreatmentSST: 0` (`_config_sst_node`, ノズル NS の生産設定) | 再利用。平板の既往検証は M 0.2・$T_w$ 320 K (case/26 `run_0023`/`run_0027`: q_w が Colburn と 1.5 %) のみ。**超音速・強冷却は未検証** — 本計画 S2 |
| 壁関数 × 等温 | `sstEnergyWallFunction: 1` (Kader q_w)。case/40 ベル部 (M≈4, $T_w/T_{aw}$≈0.4) で **+87 % 過大** | **使わない** (§2) |
| 積分法初期壁 | `feedback/deltastar_integral.py` (`thermal_bc.mode: adiabatic | prescribed_temperature`, `Tw` / `Tw_table`)。V3 単体試験のみ | 再利用。**NS 側と同じ $T_w$ を渡す配管が無い** (現状 `deltastar_loop` の既定は断熱固定) — S1 |
| δ\* 抽出 | `metrics/deltastar.py::deltastar_from_run` (ρu 質量収支、帯局所参照) | 熱境界条件に依らない定義なので**そのまま使える**。冷却で $\rho u$ 欠損が減り δ\* が小さく (強冷却では負にも) なるのは物理 |
| bcond 生成 | `runner_wt._bcond(p, euler)`: `wall_kind = "slip" if euler else "wall"` 固定。`runner_sern.py` L169 も同形 | S1 で `p.spec.wall_thermal` を読み `wall_isothermal` + `Ts` を書く |
| メッシュ y⁺ | `wall_first_frac` (CPG 6.5e-5 = y⁺≈2 断熱、TP 4.5e-5 = y⁺≈1.4)。AR ≤ 1000 ゲート | **冷却で y⁺ は上がる** (§4.2)。第一セルを詰めると AR が上がるので、平板の y₁⁺ 掃引で許容 y₁⁺ を決めてからノズルの `wall_first_frac`/`ni` を選ぶ |
| SU2 | `.external/su2/bin/SU2_CFD` v8.5、平板 (case/26 `run_0049`) とノズル CPG (case/45 `run_0012`) の SST cfg・msh→su2 変換 ([`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md)) | 再利用。等温は `MARKER_ISOTHERMAL= ( wall, 300.0 )`。**素 SST (`dilatationCorrection: 0, katoLaunder: 0`) で比較する** (case/45 `run_0013`: 素 SST で δ99 ≤3 %・θ ≤0.4 %・δ\* ≤1.3 % 一致、`dilatationCorrection: 2` は BL −16 % のモデル形式差) |
| 後処理 | `tools/flatplate_bl.py` (θ/δ\*/K–S $C_f$/壁法則)、case/26 `tools/cf_node.py`・`cf_retheta_analysis.py`、case/45 `compare_bl_su2.py` | 平板は**圧縮性版** (van Driest II, 回復温度, Crocco–Busemann, van Driest 変換) を `case/48/tools/` に新設し `flatplate_bl.py` の抽出を流用 |

## 4. 設計方針

### 4.1 問題定義: `spec.wall_thermal` を単一ソースにする

```yaml
spec:
  wall_thermal: {mode: isothermal, Tw: 300.0}      # 既定 (省略時) = {mode: adiabatic}
```

- `runner_wt._bcond(p, euler)`: `mode == isothermal` かつ NS のとき `wall: {kind: wall_isothermal, floats: {Ux: 0, Uy: 0, Uz: 0, Ts: <Tw>}}`。Euler は従来どおり `slip`。`runner_sern.py` の壁行も同じヘルパを通す。
- `feedback/deltastar_loop.run_pass0_integral`: `thermal_bc` は**常に** `spec.wall_thermal` から作る (`isothermal` → `{"mode": "prescribed_temperature", "Tw": Tw}`)。`--init-thermal` の上書きは**廃止** (codex m1: NS と積分法が別の壁温を読む状態を作らない。不一致の初期化実験が要るなら明示の例外として run に記録する)。
- `prepare_ns` の `prepare_info.json` に `wall_thermal` を記録し、`collect` の metrics に **壁熱流束の積分 $Q_w=\int q_w\,2\pi r\,ds$ ($ds=\sqrt{1+(dr/dx)^2}\,dx$; 平面は単位幅 $\int q_w\,ds$) と $q_w(x)$ のピーク位置** を追加 (等温のときのみ。断熱は 0)。
- **低 Re の $q_w$ は解像勾配から後処理で取る** (codex M2 採用): 既存の bvar `qwall` は壁関数経路専用で低 Re では 0 のまま (case/40 `run_0048` の `res_wall_3_12000.h5` で全点 0 を確認)。node では壁ノード $T_w$ と壁法線方向の第 1・第 2 内点から 2 次片側差分で $\partial T/\partial n$ を作り $q_w = -\lambda_w \partial T/\partial n$ ($\lambda_w = \mu(T_w) c_p/Pr$)。符号は壁へ入る向きを正。**閉合検証**: 平板で $\int q_w\,ds$ と入口・出口の全エンタルピー流束差 (node の壁 Dirichlet が落とす壁エネルギー残差込み) を突き合わせ 5 % 以内。
- **符号付き δ\*** (codex M1 採用): 生産抽出器 (`metrics/deltastar.py` L379 `negative_deficit` hard 不合格、L443 平滑化 `positive=True` が負値を 0 に丸める) は冷却壁で破綻する。$T_w < T_e$ の区間 (ノズルのチャンバ〜スロート: $T_e$≈1000 K に対し $T_w$ 300 K) では $\rho_w/\rho_e$≈3 で**質量欠損は負** (δ\* < 0 = 壁が実効的に外へ動く) が物理。S1 で抽出ゲート・等価半径変換・平滑化 (`positive=False`)・壁更新を符号付きに拡張し、正/零/負を横断する回帰 (合成プロファイル) を追加する。
- 変換器 (`convertGmshToForge`) の wall_dist は `wall`/`wall_isothermal` 両方の bcond を壁とみなすので変更不要 (`runner_wt.py` L210 の注記どおり)。

### 4.2 冷却壁と y⁺ (メッシュ要件)

壁単位の $y^+ = y_1\sqrt{\rho_w\tau_w}/\mu_w$ は、同じ第一セル高さ $y_1$ に対して**冷却で上がる** ($\rho_w \propto 1/T_w$ で増え、$\mu_w \propto T_w^{0.7}$ で減る。$\tau_w$ も冷却で増える)。
ノズル出口相当 (M 4.19, $T_e$ 283 K, $T_{aw}$≈1170 K) で $T_w$ 300 K なら $\rho_w$ ×3.9、$\mu_w$ ×0.39、$\tau_w$ ×1.2〜1.4 の見積りで **$y^+$ は $\sqrt{3.9\times(1.2\text{〜}1.4)}/0.39$ = ×5.5〜6** (codex m2 で算術を訂正。実測は case/48 run A の y₁⁺ 0.10 [3 µm, 断熱] と run B で確定する) — 断熱で y⁺1 のメッシュは冷却壁で y⁺5〜6 になり low-Re SST の前提を外れる。
対策は第一セルを詰めることだが AR ゲート (≤1000) と競合するので、**S2 の平板で冷却壁の y₁⁺ 掃引 (0.5 / 1 / 2 / 4) を取り、δ\*・q_w が y₁⁺ に依らなくなる上限を実測してからノズルの `wall_first_frac` と `ni` を決める**。
run ごとに実測 y₁⁺ (壁 $\tau_w, \rho_w, \mu_w$ から) を台帳に出し、上限超えは `SUSPECT` にする。

### 4.3 検証 1: 超音速冷却壁平板 (`case/48.flat_plate_cooled_m4`)

**条件** = case/44 va3 風洞の試験部壁を模す: 空気 CPG (γ 1.4, R 287, Sutherland, Pr 0.72, Pr_t 0.9)、$M_e$ 4.19, $P_e$ 5037 Pa, $T_e$ 283 K ($T_{aw} = T_e(1 + r\frac{\gamma-1}{2}M_e^2)$, $r = Pr^{1/3}$ → ≈1170 K)、$Re/m ≈ 5\times10^6$、平板長 1 m ($Re_L$ 5e6, $Re_\theta$ 帯 ≈ 2000–6000)。
メッシュ = case/26 `flat_plate_planar.geo` を派生 (平面 2D、node。押し出し 2 ノードは MUSCL 散逸消滅で不可 [[node-2node-spanwise-muscl-zero-dissipation]])、上流 slip 助走 0.1 m、上面 slip (前縁波の反射は $x$≈1.6 m で板外)、入口 = 超音速一様 Dirichlet、出口 = `outlet_statPress` ($P_s = P_e$ 一致)。
壁法線は冷却壁で y₁⁺ 0.5 になる第一セル (≈1 µm) を基準に、y₁⁺ 1 / 2 / 4 の粗化メッシュを同じ `.geo` のパラメータで生成。

**run 行列** (すべて node・低 Re SST・SLAU・MUSCL 2 次・陰解法; 段階起動は [`procedures/divergence-and-startup.md`](../../procedures/divergence-and-startup.md)):

| run | 壁 | 目的 |
| --- | --- | --- |
| A | `wall` (断熱) | 基準。$T_w = T_{aw}$ の実測 (回復係数) |
| B | `wall_isothermal` $T_w$ 300 K ($T_w/T_{aw}$ 0.26) | 主対象 (強冷却) |
| C | `wall_isothermal` $T_w$ 700 K (≈0.6) | 中間点 — $T_w$ 依存の傾きを 3 点で取る |
| B-y⁺ | B を y₁⁺ 1 / 2 / 4 メッシュで | §4.2 の許容 y₁⁺ 決定 |
| A-plain / B-plain | A / B で素 SST (下の対応表) | SU2-A / SU2-B の**基準対** (codex M7: 素 SST は A/B 両方に要る)。生産設定 (dilatation 2) との差は別途 |
| SU2-A / SU2-B | 同一 `.su2`、SST-2003m、`MARKER_HEATFLUX 0` / `MARKER_ISOTHERMAL 300` | コード間比較 |

**SST 条件の対応表 (codex M7 採用: 2 キーだけでは「同じ SST」にならない)**: forge 素 SST = `dilatationCorrection: 0, katoLaunder: 0, sstOmegaProdFromPk: 0, sstSigmaBlend: 0, sstEnergyIncludesK: 0, sstNodeWallKPin: 1` (2026-09-08 の既定変更 2 件を明示的に旧値へ)、`wallTreatmentSST: 0`、乱流輸送 1 次風上 ↔ SU2 `KIND_TURB_MODEL= SST` + `SST_OPTIONS= V2003m`, `MUSCL_TURB= NO`。物性: Sutherland ($\mu_0$ 1.716e-5, $T_0$ 273.0, $S$ 111.0 — forge の定数に SU2 を合わせる)、$Pr$ 0.72 (constant-Pr 伝導)、$Pr_t$ 0.9。入口 k/ω は同じ $(k_\infty, \omega_\infty)$ = (75 m²/s², 26000 1/s) (SU2 は TI 0.5 % / $\mu_t/\mu$ 10.1 で同値)。**流れ方向の格子感度** (nx 1000 vs 1500) を 1〜3 % 比較の前提として 1 回取る。

**理論・経験式 (何を再現できれば合格か)**:

1. **摩擦 $C_f(Re_\theta; M_e, T_w/T_{aw})$ = van Driest II** (Hopkins & Inouye 1971 が冷却壁データで最良と評価した標準形、散布 ±10 %)。非圧縮基準は case/26 で採用済みの Kármán–Schoenherr、変換は
   $F_c = (T_{aw}/T_e - 1)/(\sin^{-1}\alpha + \sin^{-1}\beta)^2$, $F_\theta = \mu_e/\mu_w$, $C_f = C_{f,i}(F_\theta Re_\theta)/F_c$
   ($\alpha, \beta$ は $T_w/T_e, T_{aw}/T_e$ の標準式)。判定は 2 段: (a) 絶対値 ±10 % (VD-II 自身の散布)、(b) **比 $C_f^{B}/C_f^{A}$ が VD-II の比と ±5 %** (forge 固有の −6 % 級バイアス [[reichardt-5pct-gap-not-forge]] は比で相殺。壁温影響そのものの検証)。
2. **熱流束**: $St = q_w/[\rho_e u_e c_p (T_{aw} - T_w)]$、Reynolds アナロジー係数 $2St/C_f$ が 1.0〜1.2 (Chi–Spalding 1.16、Colburn $Pr^{-2/3}$ 1.24 を上限側) に入る。$T_{aw}$ は run A の実測壁温を使う (回復係数の検証を兼ねる: $r$ = 0.88–0.90)。
3. **積分厚さ $\delta^*, \theta, H$**: (a) **CONTUR 積分法** (`deltastar_integral.py` の平面極限 $r_w \to \infty$、同じ $T_w$) の $\delta^*(x), \theta(x)$ と ±15 % (積分法の精度)。**比 $\delta^{*B}/\delta^{*A}$ は ±10 %** — δ\* チェーンの初期壁がそのまま冷却壁でも使える証拠になる。(b) **温度–速度関係は診断** (codex M5 採用: 古典 Walz/CB 形は CONTUR の閉包と同じ式なので独立参照にならず、冷却壁では Duan–Martín 型 [線形項係数 $C_T$=0.8259, Chen–Gan–Fu JFM 2025 式 1.1a] と $u/u_e$=0.5 で約 8 % 違う)。forge の $T(u)$ を Walz 形と Duan–Martín 形の両方と重ね、**どちらに近いか**を記録する (合否にしない)。CONTUR 平面版は巨大 $r_w$ の流用ではなく、一定外縁条件で $d\theta/dx = C_f/2$ を積分する平板入口 (`flat_plate_integral`) を用意する。(c) $H$ の圧縮性関係 (Walz): $H = H_i T_w/T_e + (T_{aw}/T_e - 1)\cdot(\ldots)$ を CB 求積で作った値と比較 (診断)。
4. **速度分布**: van Driest 変換 $u^+_{VD}$ が対数則 (κ 0.41, B 5.0) に乗るか (診断。強冷却では Trettel–Larsson 変換のほうが良いことが知られており、両方を図示するが合否にしない)。
5. **SU2 同一メッシュ**: 素 SST 同士で $C_f(x)$・$q_w(x)$・$\delta^*(x)$・$\theta(x)$ が **3 % 以内** (x = 0.3–0.9 m 平均。case/45 断熱の実績 δ\* ≤1.3 %)。生産設定 (dilatation 2) との差は「モデル形式差」として台帳に別掲。
6. **ゲート** (codex M3 採用: 「全列 falling/flat」だけでは高い残差の停滞も通るので定量化する): (i) `check_convergence.py` PASS、または warm 床のときは**全列の最終残差が断熱基準 run A の床以下**かつ NaN 0 (元の `NOT CONVERGED` 表示は残す)。(ii) 報告量の時系列判定は既存の `check_quasisteady.py --quantity theta,cf_retheta` (平板専用) に加え、本 case の `tools/cooled_plate_eval.py --series` で $q_w$・$Q_w$・δ\*・θ・$C_f$ の全スナップショット時系列を出し、**末尾 50 % の drift が比較公差の 1/3 以下** (δ\*/θ: 1 %、$q_w$: 1.5 %、$C_f$: 1.5 %) を STEADY とする。(iii) `check_mesh_quality.py` PASS、(iv) 実測 y₁⁺ ≤ 掃引で決めた上限 (超過は生産ゲートで不合格)、(v) **起動ゲート**: step 0/1 の壁 P・ρ・T・ω が有限で、壁 P が自由流の 0.5〜2 倍 (codex M6)。

### 4.4 検証 2: ノズル形状 (CPG) × SU2 等温

case/45 CPG チェーン (`problem_d155_cpg_ns.yaml`, forge `run_0013_cpg_ns_plainsst` ↔ SU2 `run_0012_su2_sst` が断熱で一致済) に `wall_thermal: isothermal 300 K` を足し、
forge (素 SST, `run_0013` から index コピー warm start) と SU2 (`MARKER_ISOTHERMAL`, `run_0012` から restart) を同じメッシュで回す。
メッシュは §4.2 の結果で決める: 現行 y⁺≈2 (断熱) は冷却で y⁺≈10 になるので、**平板の掃引で許容 y₁⁺ を超えるなら `wall_first_frac` を詰めた新メッシュを両者に使う** (AR 超過時は `ni` を増やす)。
**引き継ぎ手順は 3 種を分ける** (codex M6 採用: 「同一 index なら安全」は温度変更と再メッシュに当てはまらない): (a) **同一メッシュで壁温変更** — 温度ピンは ρ を保って $P=\rho R T_w$ を再設定するので 1170→300 K の切替直後に壁圧が 0.26 倍に落ちる過渡が出る。切替後は必ず soft 段 (1 次, cfl 0.5, 2000 step) を挟み、起動ゲート (step 0/1 の壁 P/ρ/T/ω) を通す。必要なら中間温度 (700 K) を経由する。(b) **同一トポロジで形状変更** (δ\* 反復の pass 間) — index コピー。(c) **解像度変更** (`wall_first_frac`/`ni`) — `interp_field.py` (原始変数補間) + soft 段。SU2 側は `SU2_SOL`/内蔵補間でなく新メッシュで cold start (段階 CFL) し、両者の収束を独立に確認する。
比較量 = `compare_bl_su2.py` の 4 ステーション (x/r_t 40/60/80/94) の δ99/δ\*/θ に $q_w(x)$ と $Q_w$ を追加。判定: δ\* ≤3 %・θ ≤1 %・$q_w$ ≤5 % (等温は温度場が新たに効くので断熱より緩める)。

### 4.5 生産: TP 風洞 (case/44 va3) の等温 δ\* 反復

2 段に分ける (codex M4 採用: 等温で初期壁を作り直すと形状補正が壁温影響を打ち消し、「壁温だけの感度」にならない):

- **S4a 固定形状の壁温感度**: 断熱の生産 run `run_0107` (case/44 V4 pass 0) と**同じ物理壁・同じメッシュ**で `wall_isothermal` 300 K を回す (手順 §4.4 (a): 同一メッシュ + soft 段)。台帳 = δ\*_exit / ṁ 比 / 出口面コア M / 軸 M 波 / **出口の全温・全圧分布** (壁熱移動で全エンタルピー・全圧が変わる分。「風洞では δ\* 経由でしか効かない」は撤回) / $Q_w$ / $q_w$ ピーク位置 / 実測 y₁⁺。これが §4.6 の「感度ブラケット」の実測。
- **S4b 300 K での再設計**: `problem_va3_M4.19_Lc8_iso300.yaml` (`wall_thermal` 300 K) で生産レシピ [[deltastar-production-recipe]] を回す: pass 0 = 積分法初期壁 (`thermal_bc` = 300 K、符号付き δ\*) + NS (Euler 基準 `run_0091`) → pass 1 (ω=1.0)。ゲートは同じ (|ṁ_NS/ṁ_E − 1| ≤ 0.3 %、出口面コア M ±0.1 %)。台帳には **S4a の感度と S4b の回復量を別項目**で記す。
- 積分法 (`Tw` 指定) は候補選別の近似にとどめ、最終採否はエネルギー込みの NS 評価に置く。

### 4.6 最適設計で壁温影響をどう評価するか (方針の提案)

壁温は目的関数に **(i) δ\* (実効輪郭 → 出口 M・一様性・推力係数)、(ii) $C_f$ (摩擦損失 → SERN の $C_T$)、(iii) $q_w$ (熱負荷 → 冷却設計の制約)** の 3 経路で効く。
一方 $T_w$ は設計者が選ぶ量ではなく**冷却方式と運転で決まる環境量** (断熱 ↔ 300 K の間のどこか、しかも $x$ 分布)。したがって:

1. **$T_w$ は dv にしない。作動点と同じ「環境シナリオ」として扱う** (SERN の多作動点重み付けと同じ枠組み)。既定シナリオ = {断熱, 等温 $T_w^{\rm nom}$ (冷却設計の公称値。水冷銅なら 300–400 K)}。
2. **まず感度ブラケットを取る**: 公称形状 1 点で断熱と $T_w^{\rm nom}$ の 2 run を回し、目的量の差 $\Delta f = f(T_w^{\rm nom}) - f({\rm ad})$ を台帳に出す (§4.5 が風洞の実測)。
   - $|\Delta f|$ が設計公差より小さい (風洞: 出口コア M ±0.1 %・ṁ 比 0.3 %; SERN: $\Delta C_T$ < 0.002 [SERN plan R6(c)]) → **断熱で設計し、台帳に壁温ロバスト性を記す**だけでよい。
   - 大きい → **公称 $T_w$ で設計する** (本計画の等温 δ\* 反復)。オフノミナル ($T_w$ の不確かさ ±ΔT) は同じブラケットで再評価し、公差内なら終了。
   - それでも公差外 (壁温不確かさが目的量を支配する) → **ロバスト設計**: 既存 MOO の目的ベクトルに $T_w$ シナリオを作動点として束ね (期待値 or 最悪値)、パレートを取る。
3. **δ\* 経路の近似評価が使える場面**: 風洞の一様性目的は主に δ\* 経由で壁温を見る (熱負荷は制約側) が、壁熱移動は出口の全温・全圧分布も変える (S4a で実測)。δ\* の $T_w$ 依存は積分法 (CONTUR, `Tw` 指定) で NS なしに見積れるので、**MOO の内側では積分法で $T_w$ 感度を先に篩い、NS は採用点だけ**にする (S2 で積分法の比 δ\*_cold/δ\*_ad が ±10 % で当たることを確認するのはこのため)。最終採否は NS。
4. **CHT が要るのはいつか**: $T_w(x)$ が未知で、しかも δ\*・$q_w$ が $T_w(x)$ の分布形に敏感なとき (薄肉・再生冷却・局所ホットスポット)。その場合も**フル CHT (固体伝導ソルバの連成) の前に「弱 CHT ループ」** — NS の $q_w(x)$ → 1D 壁伝導 + 冷却剤熱伝達モデル → $T_w(x)$ → `wallProfile` で NS 再実行 — を推奨する。等温壁機構と積分法の `Tw_table` がそのまま使え、固体側は解析式なので実装コストが小さい。フル CHT へ進む判断は**反復の収束/不収束では決めない** (反復が収束しないことからモデル不足は判定できない — codex 2026-09-19)。判断基準は [boundary-conjugate-heat-transfer](boundary-conjugate-heat-transfer.md) §4.8 の**モデル感度と厚さ方向近似の評価** (`local1d` vs `shell2d` の差、$\mathrm{Bi}_t=ht/k_s>0.1$、角部の熱橋) に置く。
5. **台帳の標準項目** (全チェーン共通): `wall_thermal` / 実測 y₁⁺ (max) / $Q_w$ / $q_w$ ピーク位置 / δ\*_exit / 目的量の断熱比。

### 4.7 CONTUR の熱力学・輸送の整合 (第 1 層、2026-10-08 起案)

背景と諮問: §5.1 #10 (壁温が頻繁に変わる前提、§8-4)。codex diagnose 2026-10-08 ([記録](../../notes/reviews/2026-10-08-contur-property-temperature-diagnose.md)) の 2 層分離に従う。
第 1 層は「NS と同じ気体・同じ輸送物性で書く」整合で、選択の余地を持たせない。第 2 層 (圧縮性変換の参照温度、h(v) の 2 次分布・r・a・N の閉包) は今の van Driest II 型を基準に据え置く。
§5.1 #10a の A/B で、熱閉包の書き方だけで冷却壁の δ_r が最大 2.8 % 動くことを確かめた (第 1 仮説を支持)。

1. **熱閉包をエンタルピーで書く**: $h(v) = h_w + a(h_{aw}-h_w)v + [h_e - a(h_{aw}-h_w) - h_w]v^2$、$T = h^{-1}(h)$、$\rho/\rho_e = T_e/T$ (境界層内で圧力・組成が一定)。
   $h(T)$ は縁の状態と同じ気体 (semi-perfect なら NASA-9 の $c_p(T)$ の積分)。CPG では今の温度形と式として一致する。
2. **断熱壁の回復をエンタルピーで**: $h_{aw} = h_e + r(h_0 - h_e)$、$r = 0.72^{1/3}$。r は「NS の断熱壁温に合う回復の閉包」として記録し、分子 Pr とは呼ばない (#10d)。
3. **粘性を NS と同じにする** (#10c): `gas.transport` があれば NS の viscMethod 2 と同じ CEA の種ごとの粘性と混合則 (データは `solver_density_cuda/data/species/forge_transport_v1.yaml`)。
   TP で `gas.transport` が無ければ止める。CPG と明示した旧方式だけ Sutherland。効く場所は μ_e・μ_w (F_Rδ)・N の R_δ・入口 θ₀・R_θi の床。
4. **運動量式の加速の項を気体に整合させる** (#10b): Eq. 61 の $(2 - M^2 + H)/(M(1 + (\gamma-1)M^2/2))\,dM/dx$ を
   $(2 + H - M^2)\,d\ln u_e/dx$ に置き換える (等エントロピー・組成一定なら $d\ln\rho_e = -M^2 d\ln u_e$)。
   上流の亜音速枝は面積から $d\ln u_e/dx = -(2 r_w'/r_w)/(1 - M^2)$ (γ を使わない)、下流は縁の $u_e(M)$ を気体から作って微分する。CPG では今と一致する。
5. **版で切り替える**: `deltastar_initializer.closure_version` (`contur_v1` = 今 / `contur_v2` = 1〜4)。既定は当面 `contur_v1` のまま
   (今の生産の壁をビット同一で再現できる)。既定を `contur_v2` に替えるのは §6 V-c45 で熱閉包の形が NS で裏付けられた後に、ユーザの決定で行う
   (そのとき生産の YAML は `contur_v1` を明示して再現性を残す)。版は `delta_r_initial.json` の settings と prepare_info に残す。
6. **較正の流用の制限** (#11) は別項目。`contur_v2` で較正し直した k_f は `contur_v1` の k_f と混ぜない。

## 5. 実装ステップ

1. **S1 配管** — `design/forge_design/probdef.py` (`spec.wall_thermal` 既定・検証)、`evaluate/runner_wt.py::_bcond` (+ 共通ヘルパ `wall_bcond_line`)、`evaluate/runner_sern.py` L169/L181、`feedback/deltastar_loop.py` (initializer 既定を spec から)、`evaluate/runner_axismach.py::prepare_ns/collect` (帳簿: `wall_thermal`, y₁⁺, $Q_w$)。
2. **S2 平板** — `case/48.flat_plate_cooled_m4/`: `mesh/flat_plate_cooled.geo` (パラメータ化 y₁)、`gen_runs.py` (config/bcond/IC 生成 + 段階起動)、`tools/cooled_plate_eval.py` (VD-II / RAF / CONTUR 平面 / CB / VD 変換 / SU2 読み込み)、SU2 cfg。README に run 一覧。
3. **S3 ノズル CPG × SU2** — case/45 に `problem_d155_cpg_ns_iso300.yaml`、forge run + SU2 run、`compare_bl_su2.py` に $q_w$ 追加。
4. **S4 生産 TP 等温 δ\* 反復** — case/44 `problem_va3_M4.19_Lc8_iso300.yaml`、`deltastar_loop` pass 0/1、壁温感度台帳 (`case/44/README.md`)。
5. **S5 docs** — `methods/design/overview.md` 節の実装同期、`procedures/recommended-settings.md` の壁行に等温レシピ (y₁⁺ 要件)、`procedures/verification/48-flat-plate-cooled.md`、`design/CAPABILITIES.md`。

### 5.1 残作業 (優先順)

**残作業の正本はこの表**。`notes/sessions/` の引き継ぎ文書には写しとポインタだけを置く。

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~S1 配管~~ **実装済 (2026-09-12)**: `Problem.wall_thermal` / `wall_thermal_bc_integral` / `wall_bcond_line`、`runner_wt._bcond`・`runner_sern` の壁行、`prepare_ns` (thermal_bc を spec から強制・`prepare_info.wall_thermal`)、`deltastar_loop` (`--init-thermal` 廃止)。回帰: Euler `run_0037`/`run_0091`・NS `run_0107` の壁行が文字列一致 | §5-1。`wall_thermal` 未指定は断熱で**ビット同一** (回帰: Euler `run_0037`/`run_0091` に加え **NS** の `run_0107` の bcond 文字列・積分法初期壁が一致すること [codex m1]) |
| 1a | ~~**符号付き δ\*** (codex M1)~~ **実装済 (2026-09-12)**: `_negative_delta_r` (参照 ρu を壁値で外へ延長し $r_{eff}^2 = r_w^2 - D/(\pi q_w)$)、`negative_deficit` を soft 化、P-spline `positive` を `wall_thermal` で切替 (断熱は従来どおり ≥0)。`design/tests/run_deltastar_tests.py` に正/零/負回帰 (負 δ_r 復元 2〜3 %) ALL PASS | `metrics/deltastar.py` の `negative_deficit` ゲート・等価半径変換・平滑化 `positive`・`feedback/deltastar_loop` の壁更新を符号付きに。正/零/負の合成プロファイル回帰 |
| 1b | ~~**低 Re 熱流束の後処理 + 閉合** (codex M2)~~ **平板側は実装済** (`cooled_plate_eval.py`: 2 次片側差分 $q_w$、閉合 B 1.017 / C 1.018)。残 = ノズル軸対称 ($ds$ 重み) の `metrics/extract.py` | `case/48/tools/cooled_plate_eval.py` (平板) と `metrics/extract.py` (ノズル軸対称 $ds$ 重み) に $q_w$ 抽出 (2 次片側差分)。平板で全エンタルピー流束差との閉合 5 % |
| 1c | ~~**定量ゲートの実装** (codex M3)~~ **平板側は実装済** (`cooled_plate_eval.py --series`: 末尾 50 % drift, δ\*/θ 1 %・$q_w$/$C_f$ 1.5 %)。残 = ノズル用 (出口コア M・ṁ 比・$Q_w$ の時系列) と起動ゲート | `cooled_plate_eval.py --series` (末尾 50 % drift 公差)、warm 床の定量条件、起動ゲート (step 0/1 壁 P/ρ/T/ω) |
| 2 | S2 平板 — **A/B/C 完了 (2026-09-12, case/48 run_0004/0005/0006)**: B ($T_w$ 300 K): $C_f$/VD-II 0.97–0.99・$2St/C_f$ 1.16・δ\* CONTUR 比 +3.5〜6 %・θ ±0 %・閉合 1.017・series STEADY。C (700 K): 0.91–0.95・1.15・+8〜16 %・±3 %・1.018。A (断熱): 0.87–0.92 (低 $Re_\theta$ 1400–3500)。**比 B/A の $C_f$ は 1.50 (VD-II 1.35–1.39, +8 %: 目標 ±5 % 外。C/A は 1.20 vs 1.16 で内)** — forge の断熱側の不足 (既知の Reichardt 系ギャップ) が冷却で縮む。SU2 の比で切り分ける。**y₁⁺ 掃引完了 (run_0007/0008/0009)**: 冷却壁 y₁⁺ 0.46→0.9/1.8/3.6 で δ\*/θ は −1/−2/−2.2 %、$C_f$ −1/−2.6/−3.1 %、$q_w$ −1.3/−3.7/−7.0 % → **δ\* チェーンは y₁⁺ ≤ 3.6 で 2 % 内 (許容上限 3.5)、熱負荷は y₁⁺ ≤ 1 (1.5 %) / ≤ 2 (4 %)**。**SU2 平板 (中間, it≈3000)**: B の積分量 CD (∫Cf dx) forge/SU2 = 0.998・HF (∫q_w dx) 0.982 (forge 生産 SST vs SU2 V2003m)。**A-plain/B-plain 完了 (run_0010/0011)**: 生産 SST との差は $C_f$ +1〜1.6 %・δ\* +1 % (平板では dilatation 補正は 1〜2 %)。**B-plain vs SU2-B: CD +0.8 % / HF −0.8 %** (積分量; 3 % 目標内)。**SU2-B ステーション比較 (it 5000)**: $C_f$ 1.000、$q_w$ 1.008、θ 0.998、**δ\* +4 %** (3 % 目標をわずかに超過: H 4.66 vs 4.48 の壁近傍密度分布差。$C_f$/θ/$q_w$ が一致しているので離散化でなく温度–速度関係のモデル差)。**SU2-A ステーション比較**: $C_f$ 0.997–1.000, θ 0.98–1.02, δ\* 0.96–1.04 (SU2-A は $T_w$ 1145 K でまだ上昇中)。**nx 感度 (run_0012, nx 1500)**: $C_f$ +0.2 %, δ\*/θ ±0.1 % → 流れ方向は収束。**S2 の残は無し** | §4.3。合否は §6。**起動レシピ (2026-09-12 実測)**: 一様 IC + SST は前縁で step 114 NaN → 層流暖機 (model none, 1 次, cfl 0.2, 2000) → SST soft (1 次, cfl 0.3) → mid (1 次, cfl 1.0) → 2 次ランプ cfl 0.5/1/2 (各 2000) → 本段 cfl 2 (cfl 4 は前縁 x≈1〜3 mm で P 床→NaN)。IC は壁近傍 tanh ランプ (δ₀ 300 µm) |
| 3 | S3 ノズル CPG × SU2 等温 | §4.4。**メッシュ確定 (2026-09-12, `case/44/mesh_probe_iso.py`)**: 冷却壁の y₁⁺ 倍率は run_0107 の壁データから 5.0–5.4 (スロート y₁⁺ 7→37)。`mesh2d` に x 依存の第一セル (`wall_first_frac_throat` + 上下流ブレンド) を追加し、**ni 4001 / nj 113 / wff 1.8e-5 / wfft 2.5e-6 / throat_refine 30 で AR max 846 PASS、y₁⁺_cold 1.5–4.5 (452k 節点)**。AR ≤ 1000 を守ると ni ≥ 3600 が要る (y₁⁺ 1 にするには ni 20000 級で不可: §8-3)。対象は case/45 でなく **case/44 va (Tt 1060, M4.19)** に変更 (SU2 の規模; 断熱 SU2 も同メッシュで取り直す = codex M7)。**投入 (2026-09-12)**: `run_iso_chain_cpg.py` — forge `run_0112` (CPG 断熱・素 SST, 完了: 出口コア M 4.1813, y₁⁺ スロート 0.88) → SU2 `su2_va_cpg_ad` (forge 場から restart_flow_in.csv を座標対応で生成 [float32 座標ずれ 2.4e-7 m、対応一意]、`READ_BINARY_RESTART= NO`、6000 it ≈ 2 s/it) 実行中 → forge `run_0113` (300 K 同壁, 完了: 出口コア M 4.1982 = 断熱比 +0.40 %, TP 系と同じ傾き) → SU2 `su2_va_cpg_iso300` (実行中)。**SU2 6000 it 完了・暫定比較 (`compare_nozzle_su2.py`, 2026-09-12)**: 300 K 対 (run_0113 / su2_va_cpg_iso300): ṁ 0.999、出口コア M +0.28 %、q_w ピーク +4 %、局所 q_w +7〜10 %、δ\* +6〜8 %、θ 0〜9 %、τ_w +3〜6 %、δ99 一致 (両者 y₁⁺ 2〜4)。断熱対: ṁ +0.2 %、出口コア M +0.28 %、δ\* −2〜+9 %、τ_w −2 %、ただし**両コードとも壁温が未発達 (出口 T_w 739 / 784 K vs T_aw ≈ 970 K)** = 熱場の緩和が遅い。**SU2 は未収束 (rms[RhoE] 300 K +1.15 / 断熱 −0.89)** → 各 10000 it を継続中 (`_su2_cont.sh`, 順次)。§6 の判定 (δ\* ≤3 %・q_w ≤5 %) は継続後に出す。現時点の +7 % は 3 % 目標外で、y₁⁺ 2〜4 の影響 (平板掃引: y₁⁺ 1.8→3.6 で q_w −4〜−7 %) と SU2 未収束が交絡 |
| 4 | S4a 固定形状の壁温感度 (同一壁・同一メッシュ, 断熱 vs 300 K) → S4b 300 K 再設計 | §4.5 (codex M4)。台帳は感度と回復量を別項目。**投入 (2026-09-12)**: `case/44/run_iso_chain.py` → run_0108 (断熱, 細分メッシュ, run_0107 と同じ壁 `delta_r_initial.csv`) → run_0109 (300 K 同壁) → run_0110 (300 K 積分法初期壁 pass 0, 符号付き δ\*)。**run_0108 完了**: 品質 PASS (AR 846)、実測 y₁⁺ max 2.8 (上流収縮部) / スロート 0.61 / 出口 0.39 (断熱)、ṁ_NS/ṁ_E 0.9988 (run_0107 0.9992)、出口コア M 4.1869 (run_0107 4.1900, −0.07 % = メッシュ感度)、出口 T0 1060.0 K / P0 1.139 MPa (h0 逆算)、`wall_thermal_ledger.py --series` STEADY (M_core drift 0.000 %)。台帳ツール `case/44/wall_thermal_ledger.py` (y₁⁺(x)・q_w(x)・Q_w・出口 T0/P0・ṁ 比・series)。**run_0109 (S4a-B, 300 K 同壁) 完了 (24000 step, 残差 still converging, series: M_core drift 0.07 %・Q_w 9.6 % = 熱場が未静定 → run_0114 で +36000 継続)**。**固定形状の壁温感度 (暫定, run_0108 → run_0109)**: 出口コア M 4.1869 → **4.2071 (+0.48 %)**、ṁ_NS/ṁ_E 0.9988 → 1.0011 (+0.23 %)、出口 T0 コア 1060 K 不変・質量平均 1060 → 1048 K (−1.2 %, BL の熱損失)、P0 コア不変、**Q_w 7.8 MW、q_w ピーク 2.83 MW/m² @ x −0.41 r_t (スロート直前)**、実測 y₁⁺ スロート 3.0 / 出口 1.8 / 上流収縮部 max 10.8。抽出 δ_r (同壁) は冷却で 0.6〜0.75 倍 (出口 0.091 vs 0.127 r_t)、スロートで ≈0 (−2e-5)。**壁温は風洞の出口 M を設計公差 (±0.1 %) の 5 倍動かす → §4.6-2 の判定は「公称 T_w で設計」側**。run_0110 (S4b pass 0): 積分法初期壁 (Tw 300) は δ_r 出口 0.115 (断熱 0.122, −6 %) と冷却効果を過小評価 (CFD 抽出は −29 %) → pass 1 (run_0115, 符号付き抽出) で補正する設計どおり。**run_0114 (run_0109 +36000 step) 完了**: 出口コア M 4.2016 (24k 時 4.2071 → まだ −0.13 % 動く、series drift 0.03 % STEADY)、ṁ 比 1.0024、Q_w 6.40 MW (tail drift 6 %: 熱場はなお緩和中)、出口質量平均 T0 1043.8 K。**固定形状の壁温感度 = 出口 M +0.35〜0.48 %、ṁ +0.35 %、質量平均 T0 −1.5 %** (符号・桁は確定、Q_w は ±10 % 幅)。残差は rms_roOmega が末尾で上昇 (スロート極薄セル 0.5 µm の ω が要因の疑い、要確認) |**run_0110 完了**: ṁ_NS/ṁ_E 1.0017 (ゲート 0.3 % 内)、出口コア M 4.1998 (+0.23 %: ゲート ±0.1 % 外 → pass 1 へ)、Q_w 7.36 MW、q_w ピーク 2.74 MW/m²。帯局所抽出 (符号付き) はスロートで δ_r ≈ 0 (入力 0.00115)、下流で入力の 0.79 倍、上流収縮部で負値 251 点を保持 (hard 不合格 0)。**run_0115 (pass 1, ω 1.0, 36000 step) 完了 — 等温 δ\* 反復の生産ゲート達成**: ṁ_NS/ṁ_E **1.0002** (≤0.3 %)、出口面コア M **4.1909 (+0.02 %, ±0.1 % 内)**、series M_core drift 0.02 % STEADY・Q_w 6.8 MW (drift 8.6 %, 熱場は緩和中)、固定点比 (抽出/入力) 0.94〜0.99 (x ≥ 8)、スロート近傍はさらに縮む (x=1 で 0.38、スロート −3e-4)。→ **S4b 完了 (pass 1 で固定点近傍)**。残 = 熱負荷 Q_w の静定 (継続 run) と S3 (SU2 対) |
| 5 | S5 docs・CAPABILITIES | §5-5 |
| 9 | ~~**AR 緩和の A/B**~~ **完了 (2026-09-12)**: AR 4140 メッシュ (run_0116 断熱 / run_0117 300 K, 24000 step) vs AR 846 (run_0108 / run_0109, 同 step): 断熱 ṁ 比 −0.002 %・出口コア M −0.02 %・T0 同一; 300 K ṁ 比 −0.02 %・出口コア M +0.02 %・Q_w +10 % (7.8→8.6 MW)・q_w ピーク +8 % (y₁⁺ スロート 3.0→1.04, 出口 1.8→0.81 の壁解像向上分で、平板掃引の y₁⁺ 依存 [3.6→0.5 で q_w +7 %] と整合)。残差は全列 falling、発散なし。**→ 壁法線構造層の AR ≤ 5000 は平均流に無害 (0.02 %) で熱負荷はむしろ改善** (§8-3 確定) | `problem_va_R2_LU6_Lc8_ns_ar5k{,_iso300}.yaml` (ni 2401 / nj 121 / wff 8e-6 / wfft 8e-7, AR max 4140, スロート y₁⁺_cold ≈ 1.0, 290k 節点) で run_0116 (断熱, 同壁, IC run_0108) → run_0117 (300 K 同壁) → run_0118 (pass 0)。判定: run_0108/0114 との差 (出口 M ≤0.1 %, ṁ ≤0.1 %, δ_r ≤5 %, Q_w ≤10 %) と発散/残差床の有無 | `run_iso_chain_ar5k.py` (2026-09-12 投入) |
| 6 | ~~`wallProfile` CSV ($T_w(x)$ 分布の NS 入力)~~ **[boundary-conjugate-heat-transfer](boundary-conjugate-heat-transfer.md) へ移管 (2026-09-19)** | 同 plan §4.5 / §5.1 #4。`applyInletProfiles` を一般化して `ints: {wallProfile: 1}` で `wall_isothermal` の per-face `Ts` を埋める (Ts は `valueTypes==1` なのでカーネル無改修で効く) |
| 7 | ~~弱 CHT ループ (1D 壁伝導 + 冷却剤モデル)~~ **[boundary-conjugate-heat-transfer](boundary-conjugate-heat-transfer.md) へ移管 (2026-09-19)** | §4.6-4 の方針を同 plan が引き取り、Phase 1 (外部弱連成) → Phase 2 (ソルバ内薄肉シェル) の 2 段に具体化した。**種別は `wall_isothermal` のままで `ints: {conjugate: 1}` 属性**とする (新種別 `wall_conjugate` は codex レビューで撤回: `iso_wall_flag`・T ピン・粘性壁・壁距離・DPLUR エネルギー行切離しの 5 経路から漏れるため) |
| 8 | SERN (⑤) の等温評価 | 配管は #1 で共通化。検証は SERN plan R1–R7 の後 |
| 10 | **CONTUR の物性をどの温度・どのモデルで評価するか (2026-10-08 起票、担当 F)** — 壁温が頻繁に変わる前提 (§8-4) で、積分法の初期壁が壁温に正しく追随するようにする。**判断: 2026-10-08 codex diagnose ([記録](../../notes/reviews/2026-10-08-contur-property-temperature-diagnose.md)) — 「2 層分離は修正付き採用。第 1 層 = 熱力学・輸送の整合 (同じ組成・h(T)・EOS・輸送物性)、第 2 層 = 閉包 (h(u) の 2 次分布・回復係数 r・a・N) と適用範囲。まず熱閉包の表現だけを温度形→エンタルピー形に替える CFD 0 step の A/B (#10a)。#8d の再較正と第 2 層の方式選定は保留。第 2 層は当面 van Driest II 型を基準に維持」**。指摘 10 件 (Major 9・Minor 1) の採否は下の #10a〜#10f と #11・#12 に書いた。§4.7 は #10a の結果を見てから書く | 実測 (case/45、CFD 0 step、[verification-m6-axis-wave-mesh-su2](verification-m6-axis-wave-mesh-su2.md) §9 の 2026-10-08): (i) 断熱壁温 T_aw = T_e(1 + r(γ_e − 1)/2 M²) が局所 γ_e のため全温を超える (試験部 1654〜1663 K、Tt 1600 K、NS 1470〜1482 K)。エンタルピー形 h_aw = h_e + r(h_0 − h_e) (r = 0.72^{1/3}) なら NS と +1.5〜4 K。(ii) T_aw の式の違いで δ_r が 3〜5.5 % 動き (等温でも Eq. 69 の (T_aw − T_w) 項で残る)、現行の式では冷やすと出口 δ_r が厚くなる (+0.95 %、全温基準では −1.1 %)。(iii) μ は空気の Sutherland で、NS (種ごとの CEA) と 250 K で −4.4 %、1470 K で +8.9 % ずれる。案: 第 1 層 = NS と同じ気体での整合 (エンタルピー形の Walz・h_aw・NS と同じ μ(T))、第 2 層 = 圧縮性変換の参照温度を壁温を跨いだ交差検証で選ぶ |
| 10a | ~~**熱閉包の A/B (CFD 0 step、2026-10-08 事前登録、担当 O)** — codex diagnose の判別 A/B~~ **済み (2026-10-08): 第 1 仮説を支持** | A = 今の温度形、B = エンタルピー形: h_aw = h_e + r(h_0 − h_e)、h(v) = h_w + a(h_aw − h_w)v + [h_e − a(h_aw − h_w) − h_w]v²、T(v) = h⁻¹(h(v))。変える因子は熱閉包の表現だけ。r = 0.72^{1/3}・a = 1・k_f = k_N = 1・今の μ・入口 θ₀・縁条件・運動量式・求積は両腕で固定。断熱 / 1000 / 600 / 300 K を入口から出口まで積分。**前提**: CPG 極限で A と B の差 < 0.1 %、積分精度 (rtol 1e-6 と 1e-8) の差 < 0.1 %。**判定 (事前登録)**: 試験部 [40, 94] で、B による 300 K/断熱の δ_r 比の変化 (B の比 / A の比 − 1) の絶対値の最大が 1 % 以上なら第 1 仮説 (今の熱閉包は冷却による δ_r の変化に 1 % 以上の系統差を生む) を支持、1 % 未満なら棄却。NS (断熱) との差が減るかは別欄に記録するだけで、再較正はしない。スクリプト `case/45.isobutane_m6_d155/delta_contur_compare.py hform`。**結果 (2026-10-08、`_band_ab/delta_contur/hform_ab.json`)**: 前提は満たした (CPG 極限の A/B 差 最大 0.020 %、rtol の差 最大 0.0074 %)。**判定: 第 1 仮説を支持** — 300 K/断熱の δ_r 比の変化は最大 2.80 % (x = 40、基準 1 %)。600 K で 1.43 %、1000 K で 0.18 %。温度形 (A) では 300 K に冷やすと試験部で δ_r が +1.0〜+4.8 % 厚くなり、エンタルピー形 (B) では +1.9 % (x = 40) → −1.5 % (x = 94)。どちらが NS に近いかは冷却壁の NS が無いので未確定。記録のみ (NS 断熱、k_f = 1): 出口の δ_E/δ_C は A 1.047 → B 1.032、試験部の振れは 2.93 → 3.17 %。次 (第 1 層の実装) は新しい設計の壁を変えるので、ユーザの了承待ち |
| 10b | **運動量積分式の圧力勾配項の定比熱近似 (codex Major、採用、担当 O、#10a の後)** | Eq. 61 の `(2 − M² + H)/(M(1 + (γ−1)/2 M²)) dM/dx` (`deltastar_integral.py` の `rhs`) は一定 γ の d ln u_e/dM を局所 γ_e で使っている。case/45 の気体で、この式と熱力学的に整合な d ln u_e/dx の比は M 3 で 0.964、M 6 で 0.987 (2026-10-08 手元で再計算、codex の 0.964/0.985 と一致)。`d ln ρ_e/dx + (2 + H) d ln u_e/dx + d ln r_w/dx` で書く。上流の `dMdx` の解析式も同じ近似なので、片側だけ替えない |
| 10c | **μ を NS と同じ輸送モデルに (codex Major、採用、担当 O)** | `gas.transport` があれば混合気の μ(T) (NS の viscMethod 2 と同じ CEA 輸送物性)。TP で `gas.transport` が無ければ止める (`probdef.transport_for_ns` の契約と同じ)。Sutherland は CPG と明示した旧方式だけ。効く場所は μ_e・μ_w・N の Re・入口 θ₀・Re_θi の床。注意: F_Rδ·Re_θc = ρ_e u_e θ_c/μ_w なので、F_Rδ の 11 % 差がそのまま摩擦の誤差にはならない |
| 10d | **r は「NS と同じ Pr」ではない (codex Major、採用・記録)** | viscMethod 2 は `prandtlLam` を使わない (μ と λ を混合輸送から直接計算)。混合気の Pr は 0.749〜0.758。r = 0.72^{1/3} は「NS の断熱壁温に合う回復の閉包」として記録し、分子 Pr とは呼ばない。分子 Pr を NS に近いという理由で 0.72 に置き換えることもしない |
| 10e | **訂正: 「δ の不足と傾きは T_aw の誤差では説明できない」は言い過ぎ (codex Major、採用)** | `tw` は壁温だけを NS 値に替え、温度分布の中の T_aw は今の式のままだった。`taw` も全温基準の近似で、エンタルピー形の分布ではない。言えるのは「壁温だけを替えても傾きは消えない」まで。熱閉包全体は原因の候補に残る ([verification-m6-axis-wave-mesh-su2](verification-m6-axis-wave-mesh-su2.md) §9 の 2026-10-08 に訂正を記録) |
| 10f | **θ の抽出の単位の誤り (codex Major、採用、修正済み 2026-10-08)** | `delta_contur_compare.py extract` が CONTUR の θ (既に θ/r_t) をさらに r_t で割っていた (13.04 倍)。δ 比・C_f 比の集計には使っていない。コードは直した。既存の `extract.npz` の `cont_*_theta_rt` は読むときに r_t を掛けて戻す |
| 11 | **較正の壁温条件の記録とガード (2026-10-08 起票、担当 F)** | `deltastar_initializer` (k_f 等) は較正した壁温条件でしか意味を持たない (case/45 で NS の壁温を与えるだけで出口合わせの k_f が 1.0565 → 1.0682)。較正時の `wall_thermal` を記録し、prepare で `spec.wall_thermal` と照合する案。**判断: 2026-10-08 codex diagnose — 修正付き採用**: 壁温だけの照合は弱く、全 prepare を一律に止めるのは強すぎる。較正係数の**流用と生産採用**を、熱条件・組成と物性 DB・輸送モデル・熱閉包の版・入口 θ₀・検証済みの作動範囲に結びつけて制限する。較正していない初期壁での検証計算はできるようにする。壁温を変えた後の NS 1 回は再評価の開始であって、再較正の完了の保証ではない |
| 12 | **壁温分布の単一ソース (2026-10-08 起票、担当 O、#10 の後)** | `spec.wall_thermal` に分布 (表) を足し、NS の `wallProfile` CSV (ソルバは実装済み、`methods/boundary.md`「壁温分布の入力」) と CONTUR の `Tw_table` の両方をそこから作る。**判断: 2026-10-08 codex diagnose — 採用**: 物理長 [m] (スロート原点) を正本にし、CONTUR に渡すときだけ x/r_t に変換する。同じ条件かの判定は平均温度でなく、座標原点・単位・補間/外挿の規則を含む区分関数の一致で行う。表が範囲を覆わないときに `np.interp` の端値保持で黙って外挿しない。`Tw_table` は局所の T_w を代入するだけで熱境界層の発達の履歴を持たないので、急な冷却の開始や加熱への切り替えでの妥当性は別に確かめる |
| 13 | **冷却壁用の格子の表 (2026-10-08、担当 O)** | `meshing/mesh2d.py` に `wall_first_frac_table` (第一セル厚/局所半径を log 線形の表で) と `x_density_table` (x 方向の相対密度の表) を追加 (どちらも opt-in、表が範囲を覆わなければ例外)。既存の格子はビット同一 (4 構成で確認: 生産の throat/blend・一様・axis_cap・axis_gap)。`runner_axismach._mesh_params_from` が NS の mesh ブロックから渡す。冷却で y1+ は断熱の 8.5〜9.4 倍 (case/45、case/44 の 5〜5.4 倍より大きい) で、smoothstep 1 本の第一セルでは y1+ ≤ 1 と AR ≤ 5000 を同時に満たせなかった (AR 最大の見込み 19000) |
| 14 | **冷却壁の NS の対 (2026-10-08、担当 O、判定は F)** | §6 V-c45。`case/45.isobutane_m6_d155/cold_pair_mesh.py` (格子の設計: 断熱 NS run_0179 の y1+ 分布 × 局所の冷却倍率 → y1+_cold ≤ 0.8・AR ≤ 4500 の見込み、ni 4496 × nj 121)、`cold_pair.py`・`run_cold_pair.sh` (準備・実行)。run: `run_0181_ns_coldmesh_ad`・`run_0182_ns_coldmesh_tw300` |
| 15 | **冷却壁の格子の作り直し (2026-10-08 codex plan M1・M2、担当 O、方式の決定は F)** | 初版 (ni 4496) は品質 FAIL (AR 最大 6014、float32 座標)。倍精度でも傾斜壁の接線長を見ていなかった (縮流部 AR 5488)。接線長を入れた版 (ni 4719) は倍精度で AR 最大 4573 だが、**格子の線が半径方向なので傾斜壁ではスキューがある** (縮流部で最大 0.44)。AGENTS.md の AR ≤ 5000 の例外は「壁法線に沿いスキューの無い層」だけなので、スキューの有る層は AR ≤ 1000。**生産の格子も同じ状態** (スキュー > 0.2 かつ AR > 1000 が 2193 セル、最大 3571、縮流部 x −7.7〜−0.8; 記録)。厳密に守ると x 方向の細分で ni 14653 (177 万節点、スキュー ≤ 0.1 を「無い」とした場合)、縮流部だけなら 106 万節点。代案は近壁の層を壁法線に沿わせる格子 (mesh2d の改修)。どちらにするかは諮問とユーザ判断。y1+ の見積もりは法線距離 (半径方向間隔 × cos θ_w) で行い、投入前に変換後の格子 (FP64) で第一内部点の距離と y1+ の見積もりを検査する |
| 16 | **FP64 のビルド (2026-10-08、諮問中)** | 生産と同じソース (`~/forge-wallfit-bin`、HEAD e2696d8f0) の `flowFormat.hpp` の typedef 4 行だけを double にした `~/forge-wallfit-bin-fp64` (forge sha256 65be5e28…)。float32 の座標では第一セル厚 / 局所半径 3.4e-7 が 3.6 ulp で、第一セル厚が最大 25 % ずれる (倍精度の格子を float32 に丸めて実測)。速度は生産の格子 (19.4 万節点) で FP32 3.44 ms/step → FP64 6.94 ms/step (2.0 倍)。残差の床は rms_ro 3.4e-7 → 7.4e-8 (400 step の比較)。**codex diagnose 2026-10-08 ([記録](../../notes/reviews/2026-10-08-cold-pair-precision-diagnose.md)、全件採用)**: typedef だけでは不十分 — `mesh/gmshReader.hpp` が座標を `stof` で読むので double のビルドでも入力時に float に丸まる → double のビルドだけ `stod` に (float のビルドは従来どおり、c096d66c)。msh の座標は 17 桁 (`mesh.msh_digits: 17`、10 桁だと第一層厚が 0.22 % 動く)。品質は厳密な `VERDICT: PASS` (SOFT-PASS を通さない)。「同じ FP64 なら比の誤差は消える」「生産の格子は 40 ulp なので影響なし」は却下 (生産の第一層厚の量子化 −2.1〜+2.0 %、流れの量への影響は未評価として記録)。**A/B (CFD 0 step、事前登録)**: 同じ FP64 typedef・同じ 17 桁の msh・同じ変換設定で、座標の読み込みだけ A = `stof` (`convertGmshToForge_A_stof`)、B = `stod`。生成時の倍精度座標に対し、B の第一層厚の相対誤差 ≤ 1e-6・非正の層厚なし・品質 `VERDICT: PASS` なら採用。B が超えたら残る切り詰めを追い、冷却 NS は投入しない。スクリプト `case/45.isobutane_m6_d155/fp64_reader_ab.py` |
| 17 | **壁解像の判定を面積で (2026-10-08 codex plan M3、担当 O)** | `check_wall_resolution.py` の超過率・評価率は点数の割合。生産の断熱では点数 3.6 % に対し面積 (2πr ds) では約 0.2 %。境界面の実面積で数え、領域 (縮流部・スロート・試験部) ごとの値と超過位置を出す。V-c45 の壁解像のゲートは面積と領域別で登録し直す |
| 18 | **V-c45 の判定基準の書き直し (2026-10-08 codex plan M4・M5、担当 O、投入前)** | ゲートに両腕の δ_E (x = 40・70・94)・断熱壁温・冷却側の Q_w の定常を加え、`check_quasisteady` の条件 (窓の 5 枚すべて、drift・osc の許容) を固定する。plateau を許す場合も残差床の上限を数値で決める。R_NS = 各腕の 5 枚平均の δ_E の比 (共通の x 点・補間・平滑化を固定)。不確かさは相対で揃え、時間 (5 枚の幅)・抽出 (帯の係数の感度)・CONTUR の積分と表 (≤ 0.1 %)・格子 (生産の格子との断熱の δ_E の差を目安) を合成する。R_NS の区間から e_A・e_B の上下限を作り、**区間が分離したときだけ支持**、1 % 以内の主張も上限で判定 |
| 19 | **熱閉包の A/B と contur_v2 全体を分ける (2026-10-08 codex plan M6、担当 O)** | V-c45 が判定するのは熱閉包の形 (温度形 / エンタルピー形) だけ。k_f は診断用の固定値として両腕に同じ値 (1 と生産の 1.0541) を使い、新版の較正値とは扱わない。`delta_contur_compare.py hform` は `cf_scale` を渡しておらず常に k_f = 1 だった → 生産の k_f でも計算するよう直す。contur_v2 全体 (熱閉包 + 粘性 + 加速の項) の評価と再較正は別に行い、既定を替える条件にする |
| 20 | **contur_v2 の試験 (2026-10-08 codex plan m7・m8、担当 O)** | v1 のビット再現、同じ物性・同じ微分の条件での CPG 極限 (熱閉包の一致)、NS と同じ物性への一致 (μ を独立参照実装と照合) を別の試験にする。Sutherland の定数は今の CONTUR (273.15 / 110.4) と NS (273 / 111) で違うので「CPG + Sutherland で v1 と一致」とは書かない。上流の面積の式は M = 1 で 0/0 になるので音速点の近くの接続を決める。B 腕の積分・表の細分の感度を確かめ、温度・エンタルピーの表の範囲外は例外にする。**実装・試験済み (2026-10-08)**: `feedback/deltastar_integral.py` に `closure_version` (`contur_v1` 既定 = 従来どおり / `contur_v2`)、`gas/transport.py` (NS と同じ CEA の混合気の μ)、`runner_axismach._closure_kw`・`contur_mu_fn` (TP は gas.transport 必須、CPG は NS の Sutherland 1.716e-5/273/111)。`design/tests/run_contur_v2_tests.py` 全件合格: (a) v1 がビット同一、(b1) CPG 極限の式の一致 (同じ点で 1.6e-15)、(b2) CPG 極限の積分の v2 − v1 が v1 自身の縁の格子の誤差以下、(c) 生産の燃焼ガスで v2 の断熱壁温が NS と試験部で最大 0.281 %、(d) 表の範囲外は例外、(e) 混合気の μ が独立参照と 2.8e-10。**(b) の初版の基準「積分した δ_r が相対 1e-4」は不合格だった** (断熱 2.1e-4、300 K 5.2e-4): 式は (b1) で一致し、差は縁の状態の補間と数値微分の離散化から来る (v1 自身の縁の格子 4000 → 64000 の差が 1.0〜2.3e-3、積分の許容差 1e-6 → 1e-9 の差が 2〜3e-4)。1e-4 は積分器の誤差より小さく、試験の設計の誤りとして (b1)・(b2) に分けた。既存の design の試験 (deltastar・mesh_params・mesh_euler・wall_single_bspline・pw_upstream_poly・run_tests) も合格 (mesh_euler の (d) は新しい格子のキーを使う 2 本を「変更前と同一」の対象から外し、キーが入ることを別に確認) |
| 21 | **抽出の座標の A/B (2026-10-08 事前登録、codex diagnose、担当 O、CFD 0 step)** | 判断: 2026-10-08 codex diagnose ([記録](../../notes/reviews/2026-10-08-cold-pair-result-diagnose.md)) — 「判定不能を維持。閉包を変える前に抽出の座標の A/B を」。近壁を壁法線にした格子は列の中で x が変わる (x ≈ 40 で列内の幅 0.153 r_t) のに、抽出器 (`metrics/deltastar.py`) は列全体を `x[i, 0]` として積分していた。A = 今の抽出 (列をそのまま)、B = 各行を x 方向に補間して一定 x の断面にそろえる (`_load_structured(remap_constant_x=True)`、opt-in)。Euler 参照・帯 E の選び方・求積・符号・平滑化は同じ。対象は run_0181・run_0182 の判定窓の 5 枚 (CFD 0 step)。**判定 (事前登録)**: 5 枚・試験部 [40, 94] (0.25 刻み) の max \|R_B/R_A − 1\| が 1 % 以上なら「抽出の座標の不整合の寄与を支持」、1 % 未満なら棄却 (棄却しても摩擦の閉包が原因とは確定しない)。`cold_pair.py extract-b`・`ab-compare` → `_band_ab/cold_pair/extract_ab.json`。**結果 (2026-10-08、`_band_ab/cold_pair/extract_ab.json`)**: max \|R_B/R_A − 1\| = **0.060 %** (5 枚・試験部; 平均 0.019〜0.020 %) で**第 1 仮説を棄却**。各腕の δ_E は A/B で最大 0.64 % (断熱)・0.70 % (300 K)、平均 0.13〜0.15 % 違うが、比ではほぼ打ち消す (R は x = 40 で 0.7250 → 0.7254、x = 94 で 0.7869 → 0.7869)。冷却で δ_E が 21〜27 % 薄くなるのは抽出の座標の不整合ではない |
| 22 | **V-c45 の残りの指摘の採否 (2026-10-08 codex diagnose)** | 採用: judge のゲートの集約 (NaN 検査の記録を読む、収束判定の文言で判定不能・RISING・DIVERGED を不合格に、記録の欠落も不合格) と、U の移動平均の端を有効点数で割る (x = 40 の U 0.85 % → 1.24 %) — `cold_pair.py` を修正済み。ゲート 2 (残差の床): 「格子・精度の違う run の絶対値との比較」は不適切という問題提起は採用、私の代案 (同じ run の中で横ばい + 派生量の STEADY) は却下 (高い残差での停滞を通す) — **同じ格子・精度で、残差の正規化・保存の収支・対象量の精度と結びついた上限が要る** (未設計)。ゲート 4: 縮流部の超過は不成立のまま、緩めるには格子感度の裏付けが要る。「熱閉包の形を選ぶ」問いは置き換えない (相対の優劣と絶対の精度は別の問い; 今の `decision_by_kf` の「支持」は条件付きの計算結果で正式な結論ではない)。**設計の手順: 壁温を変えたら NS で再評価は必須、壁の作り直しは設計公差を外れたとき** (§4.6 の感度の判定と同じ)。contur_v2 は opt-in のまま。基準を改訂するときは旧 V-c45 と「判定不能」を残し、改訂日・観測済みのデータ・理由・新しい条件を別版で書き、同じデータの再解析は事後解析と明記する。延長 (run_0183) の結果だけでは旧基準を通せない (断熱側の残差の床の不成立が残る) |
| 23 | **摩擦の冷却応答 (2026-10-08、記録、CFD 0 step)** | NS の壁面せん断 (接線、`res_wall_3_100000.h5`) を CONTUR と同じ非粘性の縁の状態 ½ρ_e u_e² で割った C_f を両腕で比較 (`_band_ab/cold_pair/cf_compare.json`)。試験部の平均で **C_f の比 (300 K/断熱) は NS 1.487、CONTUR 温度形 1.513・エンタルピー形 1.466・v2 1.436 — 冷却による摩擦の増え方はどの形も ±4 % で合う**。絶対値は NS/CONTUR (生産の k_f) が断熱 0.89〜0.94・300 K 0.91〜0.92。NS の 300 K の C_f は判定窓の中で −0.4 %。CONTUR は冷却で θ が 1.9〜2.0 倍・H が 0.52〜0.53 倍になって δ* がほぼ変わらない。**ずれ (NS の δ_E の比 0.75 vs CONTUR 約 1.0) は摩擦でなく、θ の発達か δ*/θ の関係の側** (どちらかは NS の θ が要る; 原因は未確定) |
| 24 | **NS の θ の抽出 (2026-10-08、記録、担当 O)** | 一定 x の断面 (`remap_constant_x`) で、δ_E と同じ帯の外縁を縁として運動量の欠損から θ を作り、両腕の θ の比と δ_E/θ の比を CONTUR (θ 1.9〜2.0 倍、H 0.52〜0.53 倍) と比べる。δ_E/θ は CONTUR の H と定義が違う (codex diagnose) ので、比どうしだけを比べる。**結果 (2026-10-08、run_0181/0182 の res_100000、`_band_ab/cold_pair/theta_*.npz`)**: 300 K/断熱の比は x = 40 → 94 で、NS の θ_r が 2.66 → 2.40 (CONTUR は温度形 2.05 → 1.94、エンタルピー形 1.94 → 1.85、v2 1.93 → 1.83)、NS の δ_loc/θ_r が 0.25 → 0.33 (CONTUR の H の比 0.51〜0.53)、NS の δ_loc が 0.66 → 0.79 (抽出の δ_r 0.73 → 0.79 と整合)。NS の δ_loc/θ_r の値は断熱 11.6〜16.3 (CONTUR の H 12〜14.8 と同程度)、300 K 2.9〜5.4 (CONTUR の 300 K の H は約 7)。**冷却による形状係数の下がり方が NS (1/3〜1/4) と CONTUR (1/2) で違い、θ の増え方も NS が約 3 割大きい** (摩擦の比は ±4 % で合う; θ の差は運動量式の加速の項に H が入ることと整合)。冷却平板 case/48 (一定 M) では CONTUR の δ*・θ が NS と 6 % 以内だった。原因は未確定 (候補: 強い加速の下での形状の閉包 — 速度分布の N(R_δ) と温度分布の組み合わせ — と縮流部からの履歴)。300 K の run はまだ DRIFTING なので、延長 run_0183 で測り直す |
| 25 | **冷却壁で効く CONTUR の係数の試算 (2026-10-08、記録、CFD 0 step)** | ユーザ「NS の反復をしながら CONTUR の係数を調整する感じで合ってる?」への回答として測った。生産の CONTUR (温度形) の 300 K/断熱の δ_r の比を係数を 1 つずつ動かして NS の R (延長後、x = 40/60/80/94 で 0.726/0.762/0.785/0.793) と比べた: 今 (a 1・k_N 1) 1.05→1.01・ずれ最大 45 %、a 0.4 で 0.84→0.81・16 %、a 0.2 で 0.71→0.69・13 % (C_f の比 2.08)、k_N 1.6 で 1.04→1.00・43 %、k_f 0.75 で 0.78→0.76・8 % (C_f の比 1.13、NS 1.49)。**NS の比は下流へ 1 に近づくのに、CONTUR はどの係数でも平らか下流で下がる — 係数 1 つでは冷却壁の x 方向の分布を再現できない**。k_N はほぼ効かず、a は効くが摩擦も動かす。出口だけ k_f で合わせると x = 40 で 1 割近くずれ、摩擦も NS から離れる。どの形で冷却壁の壁を作るか (NS の δ_E を直接使う / CONTUR に x の傾きの補正を足す / 冷却の閉包を直す) は未決 (ユーザ判断・諮問) |
| 26 | **冷却ノズルの SU2 との照合 (2026-10-08 ユーザ指示「SU2 の計算も同時に進めて。SU2 と比較するために。dilatation 補正は本来は入れるのが正しいと思う」、事前登録、担当 O、判定は F)** | 目的: 「forge が間違っている説」の切り分け (冷却で δ_E が 21〜27 % 薄くなる、が forge の解き方のせいでないか)。SU2 は燃焼ガスの多成分を同じ形で解けないので **CPG (γ 1.27354、R 292.13、Sutherland 1.716e-5/273/111、Pr 0.72、Pr_t 0.9) でそろえる**。壁と格子は生産の TP の対 (run_0181/0182) と同じ (CPG の設計は CFD ピンの凍結源が TP なので拒否される — 2026-10-08 確認) で、**CFD だけを CPG にする** (準備済みの run の設定を書き換え: thermalMethod 0・viscMethod 1・cp 1360・γ 1.27354、化学種を外す)。初期値は収束した TP の解 (run_0181 / run_0183 の最終場) の ρ・U・P・k・ω を CPG の保存量に組み直し、段階起動 full。SU2 は同じ格子 (`nozzle.msh` を 17 桁のまま su2 形式に変換) で同じ節点を解く。**run**: forge (FP64、AWS) = CPG × {断熱, 300 K} × {素の SST (dilatationCorrection 0・katoLaunder 0、SU2 と同じ形), 生産の SST (dilatationCorrection 2・katoLaunder 1)} の 4 本、SU2 v8.5 (手元の PC、6 スレッド × 2) = CPG × {断熱, 300 K} の 2 本 (SST V2003m、ROE + MUSCL、`MARKER_ISOTHERMAL` 300 K、Sutherland 273/111)。**比べる量 (Euler の参照が要らないもの)**: SU2 の解を同じ節点で forge の形式に写し、(a) 帯の外縁 y_b(x) を全 run で共通 (TP の断熱 run_0181 の抽出の band_y_b) にした δ_loc と θ_r (`theta_diag` と同じ式、一定 x の断面)、(b) 壁の C_f と q_w、(c) x = 40・70・94 の断面の ρu・T の分布。**判定 (事前登録)**: forge 素 SST と SU2 で、冷却の比 R_loc = δ_loc(300 K)/δ_loc(断熱) の差が試験部 [40, 94] で max \|R_forge/R_SU2 − 1\| ≤ 2 %、各腕の δ_loc・θ_r の差 ≤ 3 %、C_f の差 ≤ 3 % (前回の断熱 CPG の照合の実績: δ* ≤ 1.3 %・θ ≤ 0.4 %・δ99 ≤ 3 %) なら「forge の解き方は冷却ノズルでも SU2 と一致」。超えたら forge 側の原因を調べる (諮問)。記録のみ: forge 生産 SST − 素 SST の R_loc の差 (dilatation 補正の効き)、forge CPG と TP の R_loc の差 (気体の効き)、q_w。**比べる前提**: 全 run で NaN なし・RISING/DIVERGED なし・比べる量 (各腕の δ_loc の x = 40/70/94 と Q_w) が判定窓 5 枚で STEADY (SU2 は最後の 5 つの出力)。満たさなければ延長してから比べる (同じ回数ずつ)。未収束どうしを「一致」と呼ばない。**codex diagnose 2026-10-08 ([記録](../../notes/reviews/2026-10-08-cold-su2-crosscheck-diagnose.md)) の採否 (全件採用) と改訂 (投入前)**: (1) 物性: TP の設定は `thermCondMethod` を消しているので、CPG では `thermCondMethod: 1`・`prandtlLam: 0.72` を明示し `transport`・`species` を外す (prep-cpg が実効の physProp を検査して記録)。(2) SST の形: SU2 V2003m は ω の生産を制限後の P_k から作り σ を F1 でブレンドする (SU2 v8.5 `turb_sources.hpp`・`turb_diffusion.hpp` で確認) — forge の `sstOmegaProdFromPk: 1`・`sstSigmaBlend: 1` (今の既定) に当たる。**§4.3 の対応表 (両方 0) は誤りとして訂正**し、素の SST = dilatationCorrection 0・katoLaunder 0・他は既定。入口の k/ω (forge は k 1・ω 18000 固定、SU2 は TI 0.077・粘性比 12.1)、壁の ω、エネルギー中の k (forge 0、SU2 は含む — 2026-10-08 にソースで確認) の違いは実効条件の表に記録。SU2 の出口は 2237 Pa (旧 cfg の 1389 Pa は流用しない)。(3) 結論の範囲は「この CPG 条件・格子・比較量で forge と SU2 の差が許容内か」まで (TP の h(T)・組成・混合気の物性と冷却の相互作用は検証の外)。不一致でも forge の誤りと即断しない。(4) 帯の外縁: 300 K の y_b は断熱の 1.21〜1.24 倍なので、**y_b = 各 x で TP の断熱 (run_0181) と 300 K (run_0183) の y_b の大きい方**、感度として ×1.25 も記録。δ_loc・θ_r は固定した観測量であって δ_E の検証の代わりではない。(5) 許容差は暫定の照合目標: R_loc ≤ 2 %、各腕の δ_loc・θ_r ≤ 3 %、C_f ≤ 3 %、**q_w ≤ 5 % (§4.4 の目標) を判定に入れる**。(6) 比べる前提に**保存の収支**を加える: 入口と出口の質量流量の差 ≤ 0.1 %、全エンタルピーの収支 (入口 − 出口 − 壁の入熱) ≤ 0.1 % of 入口 (両コード)。比べる量 (δ_loc・θ_r・C_f・断熱の壁温・Q_w) が判定窓 5 枚で STEADY。延長は各コードのゲートで決め、反復数はそろえない。(7) 生産の SST の対は dilatation と Kato-Launder の 2 因子を同時に変えるので「生産の SST 設定一式の感度」と呼ぶ **(8) 量の測り方 (2026-10-08、どちらのコードの結果も見る前に固定)**: 後処理は `case/45.isobutane_m6_d155/cold_xcheck.py` (forge は res、SU2 は restart の CSV を同じ節点の (ni, nj) 構造に並べ、同じ関数で処理)。δ_loc・θ_r は `theta_diag` と同じ式 (run_0181 res_100000 で TP の θ の npz と全桁一致を確認)、y_b は (4) のとおり。C_f = τ_w / (½ρ_e u_e²)、縁の値は y_b での断面の値を壁の x へ補間。τ_w = μ(T_w)\|∂u_t/∂n\|、q_w = μ(T_w) c_p/Pr ∂T/∂n (壁から流体へ向かう法線、流体から壁へ入る向きを正)、微分は壁節点と壁法線上の第 1・第 2 内部節点の 2 次の片側差分、μ は Sutherland 1.716e-5/273/111。Q_w = ∮ q_w 2πr ds。C_f・q_w・δ_loc・θ_r の差は試験部 [40, 94] の 0.25 刻みの点ごとの \|forge/SU2 − 1\| の最大で判定する (q_w は 300 K の腕だけ)。収支は入口と出口の列 (出口の列は壁法線の層で一定 x でない) を通る流束 2π∫(F_x dr − F_r dx) r。全エンタルピーは h0 = c_p T + ½\|u\|² に、**SU2 だけ k を足す** (SU2 の全エネルギーは k を含む: `CNSVariable::SetPrimVar` で静的エネルギー = E − ½u² − k を確認。forge は sstEnergyIncludesK 0)。準定常は判定窓の各スナップショットで δ_loc・θ_r・C_f の x = 40/70/94、断熱の T_w の x = 40/70/94、300 K の Q_w を `check_quasisteady.classify` (5 枚、drift・osc 0.1 %) にかける。**SU2 の初期値**は forge の run_0184/0185 の段 S1 直後の場 (README 参照)。SU2 の C_f・熱流束の表面出力は compact restart のため無いので、SU2 の履歴の `HF` と自前の Q_w の差を記録する (後処理の式の確認、判定には使わない) |
| 27 | **冷却壁の腕の質量流量の欠損と擬似 CFL の試行 (2026-10-08、ユーザ指示「cfl 挙げて計算してみてはいかが」、事前登録、担当 O、判定は F)** | **観測 (CFD 0 step)**: 300 K の腕は列ごとの質量流量 (節点値の台形積分、`cold_xcheck.py` と同じ式) が縮流部で単調に減る。TP run_0183 (通算 200000 step) で入口 100.41 → スロート直後 98.52 → 出口 98.47 kg/s、CPG run_0185 (本段 40000) で 100.67 → 98.63。欠損 (入口 − スロート直後) は run_0182 の 80000/90000/100000 で 3.55/3.28/3.05、run_0183 で 2.04/1.96/1.89 と縮むが、10000 step ごとの比は 0.92 → 0.965 と鈍る。断熱の腕は ±0.1 % で平ら。中心差分のセル収支では欠損は縮流部 (x < −5) の壁から 6〜30 層目に集まる (測り方の打ち切り誤差を含む)。原因は未確定。**codex diagnose 2026-10-08 ([記録](../../notes/reviews/2026-10-08-cold-arm-mass-imbalance-diagnose.md)) の採否**: (M1) 欠損を蓄積と決めない — 採用 (収支の監査を下記で行う)。(M2) 総質量の step 差は kg/s の収支でなく、Σ CFL もコード間の時間尺度でない — 採用 (各コードが独立にゲートを満たすことで比べる; memory の「Σ CFL をそろえる」はこの照合に使わない)。(M3) CFL 4〜8 の延長は旧記録だけでは正当化できない — **ユーザ判断で試行に替える** (下記; 短い A/B の場の一致を定常解の CFL 不変の証明に使わない、relax は同時に変えない、は採用)。(M4) R_NS・δ・θ・C_f を定常の比較に使わない — 採用 (V-c45 は判定不能のまま、値は未収束場の参考)。(M5) 断熱の全エンタルピーの収支 −1.03e-3 を求積誤差として通さない — 採用 (離散の収支と後処理の指標を分けて監査する; forge の h0 は `VALUE/h0` を使う — `cold_xcheck.py` を直す)。(M6) `cold_xcheck.py` のゲートは収束判定の実行失敗・判定行の欠落・判定不能を明示的に不合格にする — 採用。**run (AWS、FP64、`cold_cfl.py prep` = restart_field でビット一致・段なし、変えるのは nStepOuter・cfl・cfl_pseudo・outStepInterval・output.extraFields だけ、implicitRelax 0.7 のまま)**: `run_0190_ns_coldmesh_tw300_cfl2` (cfl 2) と `run_0191_ns_coldmesh_tw300_cfl4` (cfl 4) を run_0183 の res_100000 から 40000 step、5000 ごと、`extraFields: [res_ro, volume]`。収支の監査: `run_0192_ns_coldmesh_tw300_audit` (run_0183 から) と `run_0193_ns_coldmesh_ad_audit` (run_0181 から) を cfl 1 のまま 2 step・毎 step 出力・`extraFields: [res_ro, volume]`・`FORGE_DUMP_MASSFLUX` (最初の評価の面の質量流束と状態)。**判定 (事前登録)**: (i) 各 CFL: NaN・DIVERGED・RISING が出たらその CFL は使えない (発散した step と場所を記録)。(ii) 速さ: 欠損 (ṁ_0 − ṁ_2000、台形) の 10000 step ごとの比を cfl 1 の run_0183 (0.962・0.965) と並べて記録する。(iii) 冷却の腕の収束場として使える条件: 入口と出口の質量流量の差 ≤ 0.1 % (監査で台形の指標が数値流束と食い違うと分かったら、数値流束の側で判定し、旧指標の値も残す)、Q_w と δ_E (x = 40/70/94) が最後の 5 枚で STEADY、check_convergence に RISING・DIVERGED・判定不能なし。満たした CFL が複数あれば、互いの比較量の差を記録する (CFL 不変の証明には使わない)。(iv) 監査: 領域全体と縮流部 (列 i ≤ 2000) で、入口の流束 − 切断面の流束 と Σ res_ro を照合する。閉合の許容は入口流量の 1e-6。閉合し欠損が Σ res_ro で説明されれば「未収束の離散収支」を支持、台形だけに欠損があれば流量の評価の問題、閉合しなければ非保存を原因候補にする (codex の事前の解釈どおり) **監査の結果 (2026-10-08、`_band_ab/cold_pair/audit_run_019{2,3}_*.json`、`cold_audit.py`)**: 監査 run は親の最終場からビット一致、リミッタの基準値も親の値に固定 (下記)。節点ごとの恒等式 res_ro = −Σ(面の数値流束) は相対 3e-13 で成立 (符号・状態の対応を確認)。**断熱 (run_0193 ← run_0181)**: 数値流束で入口 98.0753 → 出口 98.0729 kg/s (×2π、差 0.0025 %)。台形の指標は入口の列で数値流束より 0.10 % 小さく、他の列で 0.00〜0.02 % — **台形の ±0.1 % は測り方の誤差**。**300 K (run_0192 ← run_0183)**: 数値流束で入口 100.499 → x_w −5.4: 99.171 → −1.3: 98.573 → スロート直後 98.524 → 出口 98.474。入口から各断面までの Σ res_ro (×2π) は 1.329・1.926・1.976・2.026 で、欠損と 1e-15 で閉合。壁・軸の境界を通る流束は 0。→ 事前の解釈により **「壁から質量が消える」は棄却、「未収束の離散収支」を支持** (減衰しきるかは未証明)。残差 (流入超過) は縮流部 x_w < −1.3 の壁から 1〜60 層目 (j 60〜119) に広がり、6〜30 層目が最大。そこのセル体積は中心部の 1/100〜1/1000。密度の更新位相ごとの上書き量 (codex の追加の検査) は未実施。**収支のゲートの指標の改訂 (事後、旧値も残す)**: forge は台形をやめて数値流束 (massflux のダンプ、または境界面の流束) で判定する。SU2 は数値流束を出さないので、入口の列の台形の誤差 (0.1 %) を避ける扱いを比べる前に決める (未決)。**リミッタの基準値** (limiterRoRef・PRef・ARef) は既定で開始場から自動で決まるので、restart した run は親と別の作用素になる (forge の警告)。run_0182 → run_0183 で a_ref が 354.11 → 359.77 (1.6 %)、つまり V-c45 の延長は作用素を変えていた。以後の延長・試行・監査は `cold_cfl.py --limiter-ref-from <親>` で親の値に固定する (run_0190〜0195)。影響の実測 (同じ場で基準値だけを run_0182 / run_0183 の値に替えた 1 step、`run_0196`/`run_0197`): 残差の変化はいまの残差の RMS の 2e-6 (ρ)〜3.5e-5 (ρe) で、V-c45 の冷却比と質量の欠損には効かない。定常解のずれは未測定 (記録は plan convection-node-wall-reconstruction §4.27 追記)。**AWS の停止 (2026-10-08 20:40〜21:07)**: CPG の 4 本に試行と監査の 4 本を足して forge が 8 本になり (FP64・57 万節点で 1 本約 2.5 GB、g5.xlarge 16 GB・スワップなし)、OS ごと応答しなくなった (journal は 20:40:07 で途切れ、OOM の記録なし)。ユーザの承認で強制停止 → 起動。CPG の 4 本は res_60000 まで残り、試行と監査は 1 step も進んでいなかったので作り直した。同時に回すのは 4 本まで。 **CFL 4 の延長と CFL 8 の試行 (2026-10-08 21:40 ユーザ指示「cfl4 で延ばして、cfl8 も試してみて」、事前登録)**: 途中経過 (15000 step) で欠損 (Σ res_ro×2π) は 2.026 → CFL 2: 1.851、CFL 4: 1.738 (5000 step ごと 0.95 前後、やや鈍化)、δ_loc はほぼ不変、θ_r は全 CFL で単調増 (CFL 4 で 5000 step ごと +0.25 %)。`run_0199_ns_coldmesh_tw300_cfl8` = run_0183 の res_100000 から cfl 8・40000 step・5000 ごと (CFL 2・4 と同じ開始場、基準値も同じ)。`run_0198_ns_coldmesh_tw300_cfl4_ext` = run_0191 の res_40000 から cfl 4 のまま 200000 step・10000 ごと (基準値は run_0183 の値のまま)、`extraFields: [res_ro, volume]`。判定は (i)〜(iii) と同じ: NaN・DIVERGED・RISING なし、欠損 (数値流束の Σ res_ro) ≤ 0.1 % of 入口、Q_w・δ_E (x = 40/70/94)・θ_r (同) が最後の 5 枚で STEADY。延長の判定区間は run_0198 だけ (作用素は同じだが restart の跳ねを含むので、その旨を書く)。CFL 8 は (i) で使えるかを判定し、使えるなら欠損の縮み方を CFL 1・2・4 と並べて記録する (延長への切り替えは結果を見てから決める)。 **CFL 8 の結果 (2026-10-08)**: NaN なし・場は物理的 (T 226〜1604 K) だが、開始直後に残差が 60〜120 倍に跳ねて横ばい (rms_ro 8e-6 → 4.7e-4、rms_roUx 0.0035 → 0.42)。5000 step の場で |res| はスロート (x_w −0.03) の壁から 7 層目に集中し、x_w −1.3〜3 の Σ|res| は 70 (CFL 4 は 0.34)。正味の欠損は 1.193 と速く減るが、試験部の δ_loc が +0.6〜0.7 %、θ_r (x = 94) が −1.1 % CFL 4 からずれる。→ **CFL 8 はスロートの壁近くで振動して収束しない見込み、使わない** (ユーザ判断 2026-10-08、7786 step で停止)。 **ライン陰解法の試行 (2026-10-08 22:00 ユーザ指示「ライン陰解法も試してみて」、事前登録)**: 遅さの見立て — 縮流部は M 0.015〜0.1 (入口 0.02) で、擬似 dt は u + c で決まるので CFL 1 では中心部の通り抜けだけで約 2 万 step。壁近くは縦横比 (最大 4500) の分さらに遅く (局所 dt が壁法線の幅で決まる)、監査の欠損もそこ (壁から 1〜60 層目) に溜まる。`lineImplicit: 1` + `lineDtDirectional: 1` は壁法線の結合を直接解き、壁の 1 つ内側以降の Δτ を壁法線でなく流れ方向の幅で決める (procedures/solver-settings.md)。手順書の「定常 M6 ノズルでは使わない」(2026-09-03) は断熱・y+≈2 の格子で CFL の上限を見たもので、この格子の収束速度は未測定。`run_0200_ns_coldmesh_tw300_cfl4_line` = run_0183 の res_100000 から cfl 4・relax 0.7・`lineImplicit 1`・`lineDtDirectional 1`・40000 step・5000 ごと・基準値は run_0183 に固定 (`cold_cfl.py --line`)。比べる相手は同じ開始場・同じ CFL の point の run_0191。**判定**: (i) NaN・DIVERGED・RISING なし。(ii) 欠損 (Σ res_ro×2π) と θ_r (x = 40/70/94) の 5000 step ごとの値を run_0191 と並べ、同じ step 数と同じ壁時計時間の両方で比べる (ms/step を同じ負荷で記録)。(iii) 試験部の δ_loc・θ_r の run_0191 との差を記録する (CFL・解法の不変性の証明には使わない)。効くと判定する条件: 20000 step で欠損が run_0191 の同じ step の値より小さく、かつ壁時計あたりでも速い。効けば 300 K の腕の延長をライン陰解法に切り替えるかを提案する (切り替えは結果を見てから決める)。 **結果 1 回目 (run_0200)**: 65 step で発散 (detectNaN が 66 step で `res_nan_66.h5`)。1 step 目から rms_ro が振動しながら膨らむ (8e-6 → 4e-5 (step 4) → 8.8e-3 (step 64))。壊れたのは縮流部の入口寄り x_w −10.8〜−10.2 (列 39〜49) のライン 11 本で、壁から軸まで全層の ρ が非有限。手前の列 38・44 の壁から 8〜13 層目で ρ が 2.5 倍・T が 6 倍に跳ねていた。原因は未特定 (方向別の dt が壁近くの Δτ を伸ばしすぎたのか、ライン解そのものか)。**次の 2 本 (ユーザ判断 2026-10-08「両方を短く試す」、事前登録)**: `run_0201_ns_coldmesh_tw300_cfl4_lineonly` = 方向別の dt を外し `lineImplicit 1` だけ、cfl 4。`run_0202_ns_coldmesh_tw300_cfl1_linedir` = `lineImplicit 1` + `lineDtDirectional 1`、cfl 1。どちらも run_0183 の res_100000 から 10000 step・5000 ごと・基準値は run_0183。判定: NaN・DIVERGED・RISING の有無、5000・10000 step の欠損 (Σ res_ro×2π) と θ_r を point の run_0191 (cfl 4) と CFL 1 の run_0183 の傾向 (0.965/10000 step) と並べる。ms/step も記録。どちらも発散したらライン陰解法の試行は止める (発散 2 回で手順外の対処はしない、AGENTS のエスカレーション条件 2)。 **結果 2 回目 (run_0202、方向別の dt あり・cfl 1)**: 129 step で発散 (`res_nan_130.h5`)、壊れた場所は run_0200 と同じ x_w −10.8〜−10.2 (列 39〜49) のライン 11 本。rms_ro は 1 step 目から膨らむ (8e-6 → 8.8e-5 (step 8) → 4.1e-3 (step 32))。→ **方向別の dt は CFL によらず同じ場所で壊れるので、この格子では使わない** (発散 2 回、原因は未特定、深追いしない)。ライン陰解法だけ (run_0201、cfl 4) は走行中。 **CFL 2・4 の試行の結果 (40000 step、2026-10-08 22:05)**: 欠損 (Σ res_ro×2π) は 2.026 → CFL 2: 1.626 (1 万 step ごと 0.946 倍)、CFL 4: 1.398 (0.911 倍、5000 step ごとの比 0.946 → 0.959 と鈍る)。NaN なし。check_convergence (本段区間): CFL 2 = NOT CONVERGED (stalled/plateau、k/ω は still converging、rms_roOmega は 2.8 万 step まで 1700〜1966 → 最後は 40)。**CFL 4 = NOT CONVERGED (stalled/plateau) で rms_roOmega が RISING** (2〜3 万 step で 40 まで下がった後、3.2 万 step から再び上がり最後の 4000 step は中央値 704・最大 1109)。→ 事前登録 (i) の文言どおりなら CFL 4 は不可。**ユーザ判断 (2026-10-08)「rising となったら NG なの? まだわからないんじゃないの」で、(i) を事後に改訂する** (旧判定「CFL 4 は RISING で不可」は残す): ω の残差は CFL 2 でも高い時期と低い時期を行き来しており、末尾の傾向 1 区間だけで不可にしない。新しい不可の条件 = NaN・DIVERGED、または rms_roOmega が開始時 (906) の 10 倍を超える、または 1 万 step の区間の中央値が 3 区間続けて上がる。使える条件は (iii) のまま (欠損 ≤ 0.1 % of 入口、Q_w・δ_E・θ_r が最後の 5 枚で STEADY) に「ω の残差が有界」を足す。これは結果を見た後の改訂 (事後) である。CFL 4 の延長 `run_0198_ns_coldmesh_tw300_cfl4_ext` (run_0191 の res_40000 から 200000 step・10000 ごと) は 22:03 に投入済みで、そのまま続ける。CPG 断熱の続き (run_0194・0195) は 22:03 にユーザ判断で止めたが、止めた時点で 39305 / 38950 step と完了の直前だった (確認せずに止めた; 通算 100000 の最終場は無く、通算 80000 の res_20000 が最後)。 | O | 監査済み・試行は実行中 |
| 28 | **ノズルの runner でリミッタの基準値を問題ごとに固定する (2026-10-08 ユーザ判断「追って対応でいい」、担当 O、後回し)** | `runner_axismach` に `runner_sern` と同じ `evaluate.limiter_ref: {length, ro, p, a}` を入れ、全段・延長・δ* の反復・断熱と冷却の両腕の solverConfig に `limiterRefLength`・`limiterRoRef`・`limiterPRef`・`limiterARef` を書く。値は生産の run が自動で決めた値を写す (作用素を今とほぼ同じに保つ)。背景と実測は #27 と plan convection-node-wall-reconstruction §4.27 追記。着手前に方針 (値の決め方・L_ref の扱い) を plan に書いて codex に諮る | O | 未着手 (後回し) |

## 6. 検証

- **単体 / 回帰**: `wall_thermal` 省略で `_bcond` 出力が現行と文字列一致 (Euler: case/45 `run_0037` / case/44 `run_0091`、NS: case/44 `run_0107` の `bcondConfig.yaml` と diff 0、積分法初期壁も一致)。`deltastar_integral` の `prescribed_temperature` 単体 (既存 V3) と新設の平板入口 (`flat_plate_integral`) の単体。符号付き δ\* の正/零/負回帰。
- **平板 (S2)** — 合否ライン (§4.3):
  - $C_f$: VD-II 絶対 ±10 %、比 B/A ±5 %。$2St/C_f$ ∈ [0.95, 1.25]。回復係数 (run A) 0.88–0.90。
  - δ\*/θ: CONTUR (平板入口) 絶対 ±15 %、比 B/A ±10 %。温度–速度関係 (Walz / Duan–Martín) は診断 (どちらに近いかを記録)。
  - $q_w$ 閉合: $\int q_w\,ds$ と入口・出口の全エンタルピー流束差が 5 % 以内 (抽出の検証)。
  - SU2 素 SST 同一メッシュ (A-plain↔SU2-A, B-plain↔SU2-B): $C_f, q_w, \delta^*, \theta$ ≤3 %。流れ方向格子感度 (nx 1000→1500) の差が同じ 3 % 未満であること。
  - y₁⁺ 掃引: δ\*・q_w の y₁⁺ 依存が ≤2 % になる上限 y₁⁺ を決める (期待 1〜2)。
  - 全 run: `check_convergence.py` / `check_quasisteady.py --quantity theta,cf_retheta` / `cooled_plate_eval.py --series` / `check_mesh_quality.py` の VERDICT を README に貼る (§4.3-6 の定量条件)。
- **ノズル CPG (S3)**: δ\* ≤3 %・θ ≤1 %・$q_w$ ≤5 % (4 ステーション)。壁温 = 300 K がノード値で再現 (ピン)。
- **生産 TP (S4)**: 生産ゲート (ṁ 比 ≤0.3 %、出口面コア M ±0.1 %) 達成、`check_quasisteady --series` STEADY、壁温感度台帳。
- **V-c45 冷却壁の NS の対 (2026-10-08 事前登録、§5.1 #14・#15〜#19 で改訂、投入前)**: case/45 の生産の壁 (`_band_ab/prod_confirm/prep/wall_repr.json`、run_0167 の入力をビット同一で再現した準備) を
  冷却壁用の格子 (`cold_pair_mesh.py`、ni 4719 × nj 121、近壁は壁法線 `wall_normal_layer [0.02, 0.8]`、msh 17 桁) に載せ、FP64 のビルドで断熱 (`run_0181_ns_coldmesh_ad`) と 300 K (`run_0182_ns_coldmesh_tw300`) を回す。
  **この判定は「この格子・FP64」の条件での比較**であり、生産の FP32 の格子に対する予測精度の主張ではない (codex diagnose 2026-10-08)。
  - **バイナリ**: ソース = 生産のバイナリの commit e2696d8f0 + `flowFormat.hpp` の typedef 4 行 (double) + `mesh/gmshReader.hpp` の座標読み込み (double のビルドは `stod`、c096d66c と同じ差分)。
    Release・CUDA arch 86・CXX flags は生産と同じ。forge sha256 65be5e28ca1aed9f…、変換器 sha256 ac88861fa04ebdf4…。両 run で同じものを使い、run の記録に全差分と sha256 を残す。
  - **投入の前提 (不成立なら回さない)**: 品質は厳密な `VERDICT: PASS` (`--ar-max 5000`、SOFT-PASS 不可)。スキュー > 0.1 かつ AR > 1000 のセルが 0。
    変換後の座標の第一層厚が生成時の倍精度座標と相対 1e-6 以内 (`fp64_reader_ab.py` で B 腕は誤差 0 を確認済み)。物理壁の差 |Δr| ≤ 1e-8 r_t。壁距離が変換し直しと相対 1e-6 以内。
    壁の bcond が断熱 `wall` / 300 K `wall_isothermal` Ts 300。熱境界条件以外 (格子・物性 DB・実効設定) が 2 本で同じ。
  - **ゲート (結果を使う前提、不成立なら判定不能として延長・見直し)**:
    1. NaN・Inf なし (全段の残差と全 res)。
    2. `check_convergence.py --segment` で RISING・DIVERGED なし。plateau は許すが VERDICT の文言をそのまま残し「収束」と書かない。
       残差の床の上限: 本段の最後の 5000 step の rms_ro・rms_roUx・rms_roUy・rms_roe・rms_roK・rms_roOmega の中央値が、生産の run_0179 の同じ量以下。
    3. 準定常: 判定窓 80000〜100000 の 5 枚すべてで、x = 40, 41, …, 94 (1 r_t ごと) の両腕の δ_E、R_NS、断熱の壁温、300 K の壁の熱流束の積分 Q_w を
       `check_quasisteady.classify` (5 枚全部、drift ≤ 0.1 %・osc ≤ 0.1 %) で STEADY。
    4. 壁解像 (`check_wall_resolution.py --weight area`、300 K の run、断熱も記録): y1+ > 1 の面積割合が全壁で ≤ 5 %、試験部 [40, 94] で ≤ 1 %、x ∈ [−1, 40) で ≤ 5 %、
       縮流部 x < −1 で ≤ 10 % (入口の角 0.05 r_t は別に記録)。超過位置を記録する。
  - **判定の量**: R_NS(x) = (5 枚平均の δ_E,300 K) / (5 枚平均の δ_E,断熱)、δ_E は生産の抽出の経路 (`extract_and_merge`、帯 E、300 K は符号付き) の平滑後の値、x は [40, 94] の 0.25 r_t 刻み。
    CONTUR の予測 R_A (温度形)・R_B (エンタルピー形) は熱閉包だけを替えた `delta_contur_compare.py hform` の経路で、同じ壁・同じ k_f を両腕に使う。k_f は診断用の固定値として 1 と 1.0541 の両方。
  - **不確かさ U(x) (相対、線形和)**: u_t = 5 枚それぞれの比 R_k/R_NS − 1 の最大の絶対値。u_ext = 抽出の感度 (`delta_r_sens`、帯の係数 1.125 / 2.25) による δ_E の相対の振れの最大を、
    2 本で足したもの (2 r_t の移動平均)。u_c = 0.001 (CONTUR の積分・表、A/B の前提で ≤ 0.1 %)。U = u_t + u_ext + u_c。格子・精度を替えた感度
    (この格子の断熱の δ_E と生産の FP32 の δ_E の比) は記録するが U には入れない (判定の範囲を「この格子・FP64」に限るため)。
  - **判定**: 各 x で R_NS が [R_NS(1 − U), R_NS(1 + U)] を動くときの |R_m/R − 1| の下限・上限を取り、e_m の下限 = 下限の最大、上限 = 上限の最大 (m = A, B)。
    **e_B の上限 < e_A の下限ならエンタルピー形を支持**、e_A の上限 < e_B の下限なら温度形を支持、それ以外は判定保留。k_f の 2 値で判定が違えば保留。
    勝った側の e の上限が 1 % 以下なら「この格子で冷却の効果を 1 % 以内で当てる」、超えるなら「冷却の効果は当てきれない (壁温を変えたら NS の δ* 反復が要る)」と記録する。
  - **記録のみ**: θ・H・C_f・q_w の 300 K/断熱の比 (NS と CONTUR A/B)、NS の断熱壁温と CONTUR の T_aw、断熱の run と生産 (run_0179、FP32・別格子) の δ_E の差、
    contur_v2 全体 (熱閉包 + 粘性 + 加速の項) の予測 (§5.1 #19、判定には使わない)。
  - **延長の決め方 (2026-10-08、結果の 1 回目の後・延長の前に追記)**: ゲート 3 が不成立の腕だけを、その run の最後の res から `cold_pair.py prep-ext`
    (restart_field、同一格子・ビット一致・FP64 の型のまま、nStepOuter だけ変える、段なし) で 100000 step 延長する。判定窓は**延長した run の最後の 5 枚**。
    STEADY だった腕 (断熱) は元の判定窓をそのまま使う。延長してもゲート 3 が不成立なら、その時点の値と等比外挿の漸近値を記録して判定は保留し、
    さらに延長するかは諮問で決める。ゲート 2 (残差の床) とゲート 4 (縮流部の壁解像) の不成立の扱いは諮問で決める (結果を見た後に基準を緩めない)。
  - **1 回目の結果 (2026-10-08、`_band_ab/cold_pair/V_c45.json`)**: **判定不能 (ゲート不成立)**。
    ゲート 1 (NaN) 合格。ゲート 2: 両 run とも `NOT CONVERGED (stalled/plateau)`、RISING なし、**残差の床が生産の run_0179 を上回る** (rms_ro 7.0e-6・9.1e-6 vs 3.4e-7 ほか全量、10〜40 倍)。
    ゲート 3: 断熱の腕は δ_E・壁温とも全点 STEADY、**300 K の腕は δ_E が全点 DRIFTING** (5000 step ごとの増分が 0.65〜0.85 倍ずつ減る単調増加)、
    **Q_w も DRIFTING** (−12.34 → −11.49 MW、増分の比 0.91〜0.92)。ゲート 4: 300 K の全壁 1.8 %・試験部 0 %・[−1, 40) 0 % は合格、**縮流部 (入口の角を除く) 24 % で上限 10 % を超える** (最大 y1+ 1.23)。
    記録 (判定には使わない): R_NS = 0.725 (x = 40) → 0.787 (x = 94)、**冷却で δ_E が 21〜27 % 薄くなる**。CONTUR の予測は温度形 1.05 → 1.01、エンタルピー形 1.02 → 0.99、
    contur_v2 1.01 → 0.98 で、どの形でも冷却の効果をほぼ再現しない (e は 0.40〜0.46)。300 K の δ_E の等比外挿の残りは x = 94 で約 +0.17 mm (R にして +0.004)。
    格子と精度を替えた感度 (この格子の断熱の δ_E / 生産の FP32 の δ_E) は試験部の平均 +0.32 %・最大 0.52 %。原因は未確定 (諮問の前に書かない)。
    **codex diagnose (2026-10-08、[記録](../../notes/reviews/2026-10-08-cold-pair-result-diagnose.md))**: 判定不能を維持。「冷却で 21〜27 % 薄くなる」「CONTUR が冷却の効果を外す」は保存した窓の抽出値とモデルの差であって、定常解に対する確定した予測誤差ではない。次は抽出の座標の A/B (§5.1 #21)。採否は §5.1 #22。
  - **2 回目 (延長、2026-10-08、300 K = `run_0183_ns_coldmesh_tw300_ext` の最後の 5 枚、断熱 = run_0181 の元の窓、`_band_ab/cold_pair/V_c45.json`; 1 回目は `V_c45_v1.json`)**: **判定不能 (ゲート不成立)**。
    ゲート 1 合格 (延長も 21 枚 CLEAN)。ゲート 3: STEADY でないのは x = 40 の δ_E・R (TRANSIENT-UNSETTLED) と Q_w (−9.79 → −9.55 MW、DRIFTING、増分の比 0.956、等比外挿の残り約 1.3 MW) の 3 つ (1 回目は 111)。
    ゲート 2: 延長の run で rms_roOmega が RISING (`check_convergence --segment`: init 1.17e3・peak 1.92e3・fin 9.1e2)、残差の床は不成立のまま。ゲート 4: 縮流部 23 % (上限 10 %) で不成立のまま。
    記録: R_NS = 0.726 (x = 40)・0.762 (60)・0.785 (80)・0.793 (94)、1 回目からの変化は +0.1〜+0.8 %。e は温度形 0.44〜0.46・エンタルピー形 0.40〜0.41 (k_f 1 と 1.0541)。U 0.6〜1.4 %。
    延長後の θ の比 (`theta_run_0183_*.npz`) は x = 40 → 94 で θ_r 2.85 → 2.53、δ_loc/θ_r 0.23 → 0.32、δ_loc 0.66 → 0.80 (1 回目より θ の差がさらに大きい)。**延長しても「冷却で δ_E が 21〜27 % 薄くなり、CONTUR はどの形でもそれを再現しない」は変わらない**。
    さらに延長するか、ゲート 2・4 をどう作り直すかは未決 (諮問・ユーザ判断)。

### 6.1 レビュー記録 (codex)

[`AGENTS.md`](../../AGENTS.md) 「codex レビュー」の記録表 (`solver_density_cuda/tools/codex_review.py`)。

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-12 | [2026-09-12-tooling-nozzle-isothermal-wall-chain-plan.md](../../notes/reviews/2026-09-12-tooling-nozzle-isothermal-wall-chain-plan.md) | GO-with-changes, C0/M7/m2 | **全件採用**: M1 符号付き δ\* → §4.1 + §5.1-1a / M2 低 Re $q_w$ 後処理 + $ds$ 重み + 閉合 → §4.1 + §5.1-1b / M3 定量ゲート (`--quantity cf`・`--series` の誤記も訂正) → §4.3-6 + §5.1-1c / M4 固定形状感度と再設計の分離 → §4.5 (S4a/S4b) / M5 温度–速度関係は診断 + Duan–Martín + 平板積分入口 → §4.3-3 / M6 引き継ぎ 3 種 + 起動ゲート → §4.4 / M7 A-plain 追加 + SST 条件表 + nx 感度 → §4.3 / m1 `--init-thermal` 廃止 → §4.1 / m2 y⁺ ×5.5〜6 → §4.2。スポット検証 (codex): `run_0048` の `qwall` 全点 0、`run_0013`/`run_0107`/`run_0048` は既定条件で NOT CONVERGED (stalled) — 本計画の warm 床条件はこれを踏まえて定量化 |
| plan | 2026-10-08 | [2026-10-08-tooling-nozzle-isothermal-wall-chain-plan.md](../../notes/reviews/2026-10-08-tooling-nozzle-isothermal-wall-chain-plan.md) (§4.7・§5.1 #13・§6 V-c45 に集中) | GO-with-changes, C0/M6/m2 | **全件採用** (§5.1 #15〜#20): M1 格子の AR (傾斜壁の接線長・半径方向の線はスキュー有りで AR 5000 の例外外) → #15 / M2 float32 の数 ulp・法線距離と半径方向間隔の区別 → #15・#16 (FP64 ビルドの諮問中) / M3 壁解像は点数でなく面積で・領域別 → #17 / M4 比だけでなく両腕の δ_E・断熱壁温・Q_w の定常・残差床の上限 → #18 / M5 e と u の尺度・不確かさの寄与・区間の分離で判定 → #18 / M6 熱閉包だけの A/B と contur_v2 全体を分ける・hform の k_f が常に 1 → #19 / m7 CPG の一致と v1 再現・Sutherland 定数・M=1 の 0/0 → #20 / m8 B 腕の精度確認・表の範囲外 → #20。V-c45 の判定基準は #18 を反映して投入前に書き直す |

## 7. 影響範囲

- `design/forge_design/{probdef.py, evaluate/runner_wt.py, evaluate/runner_sern.py, evaluate/runner_axismach.py, feedback/deltastar_loop.py}`
- 新規 `case/48.flat_plate_cooled_m4/`、case/45・case/44 の問題 YAML 変種と run
- docs: `methods/design/overview.md` (新節)、`procedures/recommended-settings.md` (壁行)、`procedures/verification/README.md` + 新ファイル、`design/CAPABILITIES.md`
- forge 本体は変更しない (等温壁・低 Re SST は既存)。変更が要ると分かった場合 (例: 冷却壁での EOS 床・ピンの不具合) は本 plan に §4.7 として追記してから触る

## 8. 未確定事項 (ユーザ確認)

1. **公称 $T_w$**: 初版は 300 K (ユーザ発言 2026-09-12「壁面温度 300 K くらいでまずは」)。実機の冷却方式が決まったら $T_w^{\rm nom}$ を差し替える。
2. **CHT の要否**: §4.6-4 の判断基準 (弱 CHT ループが収束するか) で決める。フル CHT の実装は本計画の外。
3. ~~**AR ゲートと y₁⁺**: 冷却壁で y₁⁺ ≤ 上限を守ると AR が 1000 を超える形状があり得る。その場合 `ni` を増やす (計算時間 ×数倍) か、AR 上限を「壁法線方向の構造格子は AR 2000 まで可」に緩めるかは S2/S3 の結果でユーザ判断。~~ **決着 (2026-09-12, ユーザ決定「AR 上限緩めようか」)**: 壁法線の構造格子層に限り **AR ≤ 5000** まで緩和 (AGENTS.md / calculation-workflow / recommended-settings に反映、`check_mesh_quality.py --ar-max`、問題 YAML `mesh.ar_max`)。裏付け = §5.1-9 の A/B (AR 846 メッシュ run_0108/0114 vs AR 4140 メッシュ run_0116/0117): **確定 — 出口 M・ṁ の差 0.02 %、熱負荷は y₁⁺ 改善分 +8〜10 %、発散なし。冷却壁ノズルの生産メッシュは `problem_va_R2_LU6_Lc8_ns_ar5k*.yaml` (AR ≤ 5000, スロート y₁⁺ ≈ 1) とする。**
4. **壁温は頻繁に変わる前提で設計チェーンを作る (2026-10-08 ユーザ発言)**: 「これから壁温設定はガンガン変わり得る」「等温壁にするかもしれないし、壁温分布を与えるかも」「CONTUR で物性をどの温度で評価するかは任意性がある話なので、どうするべきかは考えて」。対応は §5.1 #10〜#12。

## 9. 完了条件

- [ ] `methods/design/overview.md` の新節を実装と同期
- [ ] S1–S5 完了、§6 の判定を満たす
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映
- [ ] `status: done`、§10 変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 10. 変更ログ

- `2026-10-08` — **壁温が頻繁に変わる前提 (§8-4、ユーザ発言) を受けて §5.1 #10〜#12 を起票**。case/45 の実測で、CONTUR の断熱壁温の式が燃焼ガスで全温を超えること (+190 K)・その式の違いで δ_r が 3〜5.5 % 動くこと・μ が NS と最大 9 % ずれること・出口合わせの k_f が壁温条件で変わることを確認した (数値は §5.1 #10、詳細は [verification-m6-axis-wave-mesh-su2](verification-m6-axis-wave-mesh-su2.md) §9 の 2026-10-08)。物性の評価の方針は codex (diagnose) に諮問中で、§4 は未反映。
- `2026-10-08` — **codex (diagnose) に諮った**: [`notes/reviews/2026-10-08-contur-property-temperature-diagnose.md`](../../notes/reviews/2026-10-08-contur-property-temperature-diagnose.md) — 結論「係数と輸送物性を固定し、熱閉包だけを温度形からエンタルピー形へ替える CFD 0 step の A/B を先に」。指摘 10 件は全件採用 (§5.1 #10〜#12 に採否)。運動量式の項の近似 (#10b) と θ の単位 (#10f) は手元で再現してから採用した。第 2 層 (圧縮性変換の参照温度) の選定は保留: case/48 の断熱基準の壁温ドリフトと case/44 run_0114 の熱負荷の変動のため、今のデータでは決めない。後で足す NS の最小候補は case/45 の生産形状の 300 K (メッシュが冷却壁の壁解像を満たせば 1 本、満たさなければ同じメッシュで断熱と 300 K の 2 本)。#10a を事前登録した。
- `2026-10-08` — **#10a の A/B を実施: 第 1 仮説を支持** (300 K/断熱の δ_r 比の変化 最大 2.80 %、基準 1 %)。熱閉包の書き方だけで冷却壁の δ_r が数 % 動く。詳細は §5.1 #10a。説明ページ https://claude.ai/artifact/5dFnxMbYCD4DzbmXLk6QX7 (CONTUR の式と較正の係数)。
- `2026-09-12` — **AR 上限を壁法線構造層で ≤5000 に緩和 (ユーザ決定)**。AGENTS.md・手順書・runner (`mesh.ar_max`) に反映。A/B (§5.1-9) で平均流 0.02 %・熱負荷 +8〜10 % (y₁⁺ 改善分)・発散なしを確認し確定。冷却壁の生産メッシュを AR 5000 版 (スロート y₁⁺ ≈ 1) に切替。
- `2026-09-12` — **S4 完了 (case/44)**: S4a 固定形状の壁温感度 (run_0108 断熱 vs run_0109/0114 300 K, 同壁同メッシュ): 出口コア M +0.35〜0.48 %・ṁ +0.35 %・出口質量平均 T0 −1.5 %・Q_w 6.4〜7.8 MW (q_w ピーク 2.3〜2.8 MW/m² @ スロート直前)。S4b 300 K 再設計 (run_0110 pass 0 → run_0115 pass 1): ṁ 比 1.0002・出口コア M +0.02 % で生産ゲート達成、積分法初期壁は冷却効果を過小評価するが pass 1 で固定点近傍。S3 は run_0112 (CPG 断熱) 完了・SU2 対は実行中。
- `2026-09-12` — S1 配管・符号付き δ\*・平板評価ツール実装、case/48 run A/B/C 完了 (§5.1-2 に数値)。y₁⁺ の冷却倍率は実測 ×5.7 (3 µm: 断熱 0.08 → 300 K 0.46)。
- `2026-09-12` — codex plan 段レビュー (GO-with-changes, C0/M7/m2) を §6.1 に記録、**全件採用**して §4.1/§4.2/§4.3/§4.4/§4.5/§5.1/§6 を改訂 (符号付き δ\*、低 Re $q_w$ 後処理と閉合、定量ゲート、S4a/S4b 分離、SST 条件表、引き継ぎ 3 種)。case/48 の起動レシピ実測 (層流暖機 + 2 次 cfl ランプ、本段 cfl 2) を §5.1-2 に記録。
- `2026-09-12` — 起票 (ユーザ依頼: ノズル設計の壁を断熱から等温/CHT へ。等温 low-Re SST を平板で理論・経験式と SU2 に対して検証 → ノズルで SU2 と比較 → 等温 NS で δ\* 反復 → 壁温影響の評価方針)。ユーザ判断 2 件を反映: 壁関数は使わない (準備状況が悪い)、冷却で y⁺ が変わる点を平板の y₁⁺ 掃引で扱う。
```

## 参考: `plans/accepted/time_integration-line-implicit-viscous-v2.md`

```
# time_integration: line-implicit v2 試作 (粘性結合 + K 凍結) — 上限確認の限定実験

- status: done
- 起票: 2026-09-02
- related_docs: [methods/time_integration/implementation.md](../../methods/time_integration/implementation.md)
- 先行: [plans/accepted/time_integration-line-implicit.md](../accepted/time_integration-line-implicit.md) (v1: 対流のみ、DDES A/B で不採用)

## 目的 (製品化ではない)

v1 の DDES A/B (case/39 ny160, run_diag_lineimp_*) で確定したのは「**現行 v1 は割に合わない**」
(サブ反復収束 +8% / step 単価 2.44 倍) であって line-implicit の到達点ではない。本 plan は
**上限確認のための一段だけ**を区切って実施し、「現実装が遅かった」のか「この離散化・ケースでは
line-implicit の上限自体が低い」のかを決着させる。

損益分岐: 現収束性能 (nSub 20→17 相当) のままなら line のサブ反復単価を ctrl の **~1.18 倍以下**に
落とす必要がある。K 凍結だけで届かなければ、粘性結合による大幅な nSub 削減が必須。

## 実験項目

1. **K 凍結 (`lineKFreeze`)**: K ブロック抽出 (単位ベクトル×5 の Jacobian 列展開 = v1 の主コスト仮説)
   を物理 step の最初のサブ反復のみ実行し、以後のサブ反復で再利用。純粋な計算コスト下限を測る。
   近似の質: defect-correction の LHS はもともと近似であり、サブ反復間の状態変化は O(ΔQ)。
2. **壁法線 line へのスカラー粘性結合**: line 面に限り K_prev/K_next へ −α_f·I を加算
   (α_f = ν_eff·δ/dcc, 既存 `viscous_diag`=2α と同じメトリック)。対角は既存の 2α を維持
   (行和優位を保つ安全側。真の対称形 [−α, 2α, −α] が line 内で完成する)。
3. **粘性 CFL 割引 (`lineViscousDtRelief` θ∈[0,1])**: on-line セルに限り `setDT` の粘性
   スペクトル半径項 2ν_eff/(ρ·dx_min) を (1−θ) 倍。θ=0 (現状) → 0.5 → 0.8 → 1.0 で
   安定限界と必要 nSub を測る。dx_min は壁近傍で壁法線=line 方向なので方向整合は近似的に成立。
4. **同一プロトコル比較**: case/39 ny160・発達場 (run_0023 res_16000)・300 step 窓・
   総 wall time で ctrl (point) と比較。測るもの: (a) K 凍結後のサブ反復単価比、
   (b) θ を上げたときの安定限界、(c) 同等残差に必要な nSub、(d) 総時間が ctrl を下回るか。

## 判定

- **勝ち筋**: (K凍結後単価) × (必要 nSub / 20) < 1 × ctrl 単価 → line-implicit は本実装
  (完全粘性 5×5 Jacobian・adaptive nSub) を検討する価値がある。
- **負け筋**: θ を上げても nSub が削れない/不安定 → 「対流 defect-correction が律速で、
  粘性 line 陰化の上限は低い」と結論し、v2 は記録して閉じる。

## 変更ログ

- 2026-09-02: 起票。
- 2026-09-02: 実装・検証完了 (**勝ち筋で決着 — 「現実装が遅かった」が正**)。case/39 ny160 DDES
  発達場・300 step 窓 (`run_diag_lineimp2_*`, IC=run_0023 res_16000, ctrl=point 345 s):
  - **コスト分解**: v1 モノリシック 843 s (2.44×) の主犯は **Thomas の毎 sweep 再分解**
    (rhs 非依存の LU + W + Kprev·W 625 積を 5 回/subiter 重複)。factor/solve 分離 (厳密) で
    547 s (1.59×)、`lineKFreeze` (サブ反復間凍結) で **456 s (1.32×) = コスト下限**。
    split/kfreeze は v1 と収束軌道が一致 (分離の厳密性・凍結の無害性を同時検証)。
  - **粘性結合・dt 割引はこのケースでは僅差** (roe +0.02 桁)。理由: 壁第一セル (Δn≈3.5e-5 m)
    の λ_visc/λ_acoustic = 2ν/(Δn·c) ≈ **0.02** — 圧縮性 pseudo-dt の壁法線律速は
    **音響 (c/Δn) であり v1 の対流・音響 line が既に陰化している**。粘性 Δn² 律速が主役に
    なるのは Δn < 2ν/c (本ケース y+≈0.02 相当) の超極薄セルのみ。θ=1.0 でも安定。
  - **本命は pseudo-CFL 引き上げ**: point は cfl_pseudo 2 で発散 (step 99)、**line は 2/4/8 全て
    300 step 安定**。cp4 で ctrl@20 品質に subiter 12-13 到達 (ctrl は 19)、cp8+nSub20 は
    roUx 17 倍深い・ωバースト最大 10.8→8.5 に低減。
  - **実測の勝ち**: `cp4 + nSub13` = **302 s (ctrl の 0.88 倍) で ctrl 同水準の step 終端残差**
    (rms_roUx 5.6e-7 vs 8.0e-7)。同時間 (nSub15) なら ctrl@20 より深い。
  - **判定**: M6 定常ノズル (streamwise 対流律速 → line 無効) と DDES dual-time (壁法線音響律速
    → line 有効) で**律速モードが違う**。DDES 生産推奨 = `lineImplicit:1 + lineKFreeze:1 +
    cfl_pseudo 4 + nSub 13-15` (粘性結合/割引は任意)。
  - **未検証の注意**: 300 step 窓のみ — nSub 削減のバースト余裕 (run_0013 の教訓: バースト最悪値は
    乱流発達と共に成長) は長時間 run で要検証。cp8×nSub 削減の組合せ、adaptive early-exit は今後。
- 2026-09-03: **方向別 dt (`lineDtDirectional: 1`) 追試** — ユーザー指摘「陰化した方向の λ は
  Δτ を縛らなくてよいはず」の実証。line 面の λ (音響込み) を CFL の max から完全除外し、
  壁セルの Δτ を off-line (streamwise) 基準 (×AR≈74, BDF 物理項が対角支配する Newton 的
  レジーム) にする。結果 (`run_diag_lineimp2_dir_cp4*`, 300 step):
  - **cp4+nSub20 で安定** (478 s)。収束は line_cp4 と同等〜微改善 (roUx 3.13→3.37 桁)、
    **ωバースト最大 9.5→6.15 に低減** (壁 Δτ 拡大が k-ω 隔離更新の突発を均す)。
  - **nSub8 は品質不足** (0.77/1.35 桁 < ctrl@20 の 1.53/2.17) — 方向別 dt でもサブ反復の
    収縮率は上がらない。**収縮律速はもはや壁擬似時間でなく off-line lag / SST segregated 結合**。
    ctrl@20 品質に必要な nSub は 12-13 のまま。
  - 運用示唆: 同品質の最速点は cp4+nSub13 (~0.87×) で不変だが、**directional はバースト余裕を
    +50% 積み増す**ので nSub 削減時の安全マージンとして併用推奨。
- 2026-09-03: **記録の補正 (Codex レビュー 3 点)**:
  1. **「壁セル Δτ が AR≈74 倍」は現実装では不成立** — `setDT` の directional 除外は
     `line_prev/next` に一致する**内部 line 面だけ**で、壁ノードの**境界半割面は CFL の max に
     残る**。壁 CV は内部 line 面と境界面が同じ V/S (実測 1.408e-5 m) なので、境界面が同じ音響
     制約を残す。**AR 倍化は壁の 1 つ内側以降の line CV のみ**。壁境界面 λ の除外可否
     (Dirichlet 行 decouple 済みなら安全かもしれない) は未検証の将来項目。
  2. ωバースト半減の帰属は「directional dt により半減 (実測)」に留める — 「壁 Δτ 拡大が均した」
     は未分離の推論 (1. の通り壁セル自身の Δτ は伸びていない)。
  3. 収縮律速の候補は **3 者**: off-line lag / segregated SST / **2次 KEEP RHS×1次 FVS LHS の
     defect-correction 不整合** (M6 の cfl×relax 飽和と同根の可能性)。
- 2026-09-03: **推奨構成の直接 A/B (directional cp4 × nSub13/15)** — 結果は下記追記。
- 2026-09-03: **推奨構成の直接 A/B 完了** (`run_diag_lineimp2_dir_cp4_nsub13/15`, 300 step):
  - `directional cp4 + nSub13`: **317 s = ctrl の 0.92 倍**、step 終端品質は同等〜良
    (roUx 4.4e-7 vs ctrl 8.0e-7, roe 2.45e-4 vs 2.13e-4)、ωバースト最大 6.15 (ctrl 10.8)。
  - `directional cp4 + nSub15`: 363 s (1.05 倍) で全量 ctrl より良い (roUx 2.0e-7, roe 1.33e-4)。
  - 非 directional nSub13 (302 s, 0.88 倍) より directional は +5% 遅いが、バースト余裕 −35%
    (9.5→6.15) を買う取引。**生産候補 = directional cp4+nSub13 (速度優先) / nSub15 (品質・余裕優先)**。
  - 残作業 (本採用の条件): 生産再開時に nSub15 で長時間 (数万 step) のバースト余裕検証を 1 本。
- 2026-09-03: **壁境界半割面 Λ の除外実験 (`lineDtWallRelief`, opt-in 診断)** — Codex 提案の
  3 分岐実験。wall 種 bcond の境界面のみ明示フラグ (19642 面) で setDT の max から除外
  (inlet/outlet/periodic は残す。LHS の境界 A⁺ 対角は不変)。結果は**分岐③: 発散**
  (`run_diag_lineimp2_dirwall_cp4`, step ~80-100 で ro が非有限 → detectNaN 停止。
  対照 `dir_cp4` 再走は完走)。最初に落ちたのは **ro** (種の位置は dump 削除により未特定)。
  **解釈の補正 (Codex レビュー)**: node 強壁の RHS 壁対流流束は F_b=(0, pn, 0) で壁質量流束は
  恒等的にゼロ → **連続行に「欠けている物理的境界 Jacobian」はほぼ無い** (境界面の A⁺S 対角は
  もともと人工的な対角強化 [gpu-implicit-plan.md 参照])。今回証明されたのは
  **「壁面 Λ 由来の pseudo-time 質量項 (V/Δτ) が現行反復の壁 CV 減衰として必要」**まで。
  除外時の発散原因は境界 Jacobian 欠落とは確定できず、LHS/RHS 不整合・K 凍結・非線形更新過大を
  含む複数候補のまま。将来壁 Δτ を広げたい場合の手順: 発散直前のセル位置・δρ/ρ・V/Δτ を確認
  → 局所緩和 / 正値性ガード / K 凍結解除で原因を分ける (境界 Jacobian 追加はその後)。
- 2026-09-03: **dt_local 壁距離プロファイル (dir_cp4, `dt_profile.csv`)** が Codex 指摘を定量確認:
  壁ノード帯 (wd<5e-5) の dt_med **1.9e-6** vs 内側 **1.4-1.7e-5 (~8 倍差)** — directional の
  AR 倍恩恵は第一内点以降のみで、壁 CV 自身は境界面 Λ で絞られたまま。内側は ~7.5·dt で
  BDF 物理項支配レジーム到達済み。⇒ **「壁擬似時間律速の完全除去」は未達のまま**が正確で、
  収縮律速 3 候補 (off-line lag / segregated SST / defect-correction 不整合) の切り分けは、
  壁端点を除去できない以上、別経路 (例: SST 連成陰化 or LHS 2 次化) からになる。
- 2026-09-03: **収縮律速の消去法完了 (診断 4 腕, 各 100 step 窓 20-95, 対照 dir_cp4)**:
  ① `FORGE_FREEZE_TURB=1`: roe/roUx 収縮 **2.07/3.37 桁で対照と完全一致** → segregated SST は
  平均流収縮を妨げていない (注: μt は凍結 k,ω から毎回再計算されるので「ほぼ固定」)。
  ② `implicitRelaxSST` 0.7/1.0: 平均流不変。**ω 収縮はむしろ悪化** (0.39→0.26→0.05 桁) =
  ω は relax 0.5 で縁辺安定。SST 2×2 line 化は ω 品質・バースト向けで平均流 nSub は減らない。
  ③ `nStepInner` 5→10: **全チェックポイントで収縮一致** → 線形系は 5 sweep で解き切れている
  (off-line lag は近似 LHS の線形解の中で処理済み)。
  ⇒ **平均流のサブ反復収縮律速 = defect-correction 不整合 (1次 FVS LHS × 2次 KEEP+ES RHS) に確定**。
  根治候補: line K の FD 化 (RHS 整合をライン方向だけ入れる, v2 機構流用可) / JFNK
  (DPLUR+line を前処理に) / adaptive early-exit は律速に依らず ~20-30% 得。
- 2026-09-03: **訂正 (Codex レビュー・重大)**: 上記「消去法完了」の freeze A/B は**無効だった** —
  `FORGE_FREEZE_TURB` のゲートは定常経路 (`implicitNonlinearUpdate`) 専用で、dual-time の
  `advanceImplicitDualTime` は `applySSTPointImplicit` を無条件に呼んでいた (同一経路の 2 回走で
  一致は自明)。**「defect-correction 不整合が単独原因として確定」は撤回** (commit d3b55d6d の
  結論を本追補で差し替え)。dual-time にもゲート+起動ログを追加して再走 (`run_diag_lineimp2_frz2`):
  - 凍結の実証: res_100 の roK/roOmega は**内部 (wd>1e-4) で IC とビット一致**。変化は壁ピン帯
    (wd≤3.2e-5, roOmega の 2.5% ノード = `nodeOmegaWfDirichlet` 等の壁 Dirichlet — 両腕共通) のみ。
  - **結果: 平均流収縮は active と実質一致** (roe@20 2.07/2.07, roUx 3.37/3.36 — 0.01 差が
    「今回は本当に別経路」の傍証)。
  ⇒ 支持される結論は「**この窓では SST フィードバックは平均流収縮の主律速でない**」まで。
  nStepInner 一致の読みも訂正: 「線形系を解き切った」でなく「**同じ近似 LHS をさらに解いても
  非線形収縮は改善しない**」(sweep 増は無価値、の実用結論は不変)。
  残候補: defect-correction 不整合 / implicitRelax=0.5 の緩和上限 / lineKFreeze の古い LHS /
  off-line Jacobian の忠実度 / 壁 CV の pseudo-time 制約 → ir0.7・kfreeze-off の追試で切り分け続行。
- 2026-09-03: **残候補 probe 2 本 (各 100 step, 対照 dir_cp4)**:
  - `implicitRelax` 0.5→0.7 (`run_diag_lineimp2_ir07`): **発散** (NaN)。relax 0.5 は任意の減衰でなく
    **安定必須** — 近似 LHS の defect-correction が要求する緩和で、収縮率の上限 (毎 subiter ≤半歩)
    を課すが外せない。⇒ この候補は「defect-correction 不整合」に吸収される (LHS が忠実なら
    relax→1 で Newton 級収縮が許されるはず)。
  - `lineKFreeze` off = K/LU 毎 subiter 再構築 (`run_diag_lineimp2_nokfrz`): **収縮完全一致**
    (roe@20 2.07 / roUx@20 3.37)。**「K 凍結の古い LHS」候補は消去** — 凍結はコスト削減のみで
    収縮に無害と直接実証。
  ⇒ 総括: 安く動かせる要素 (SST/sweep 数/K 鮮度/relax) を全て振っても収縮は 2.07/3.37 に固着。
  **残るのは近似 LHS の忠実度ファミリー (1次 FVS×2次 KEEP の defect-correction 不整合 +
  off-line Jacobian 忠実度 + 壁 CV pseudo-time 制約)** で、単独犯の特定はこの窓では不能。
  打ち手は変わらず: adaptive early-exit (無条件) / line K の FD 化 / JFNK。
- 2026-09-03: **CFL-nSub マップ + 残差局在の測定 (Codex 提案の最安切り分け, run_diag_cflmap_*)**:
  - **Phase 1 (cp マップ, directional+nSub20 固定, 各 100 step)**: cp8 OK / cp12 OK / **cp14 発散
    @step3 / cp16 @step2 / cp24 @step1** → 上限は [12,14)。**収縮は cp8 で飽和** (roe@10:
    cp4 1.22 / cp8 1.27 / cp12 1.28) — BDF 項が対角支配に入った後は Δτ 増が効かない構図。
    振動兆候ゼロ (roe 増加 subiter 率 0.000)。min ro=1.16/min P=1.02e5 で EOS 床とは無縁
    (M6 の床洗浄とは終端症状が別物)。**運用最良点 = cp8** (上限へ 1.5 倍の余裕)。
  - **発散モードの指紋 (cp14 の inner 履歴)**: (ro, roUy, roe) が **subiter 間で単調増幅**
    (step0 subiter11 から成長・×1.6/3subiter)、roUz/roK は無傷 → 壁法線 (y) の音響/質量-
    エネルギー系の反復モード。NaN 拡散は 1 step で 86-100% に及ぶため種セル特定には
    per-subiter 検査が必要 (将来項目)。
  - **Phase 3 (残差局在, cp8 + `FORGE_OUT_RESIDUALS=1` [main.cpp に env 追加, res_*/dq_block_new_*
    を h5 出力へ])**: 相対補正要求 η=|dq|/max(|Q|,Q_ref) は**丘頂直前のせん断層 (x/h≈8.6,
    wd≈0.09-0.11h) に全 5 式が同座で局在**。壁ノード帯は p99 で BL 上部より小さく **壁 CV は
    支配的でない**。max η ~4e-5 と絶対値も健全。⇒ 判定木では **off-line (streamwise/spanwise)
    枝** — line K の忠実化だけでは動かず、次に効かせるなら streamwise 第 2 ライン族 (x 周期 →
    cyclic Thomas) か、物理活性域の要求として受容 (現状 cp8+nSub13 で実用充分) の二択。
- 2026-09-03: **診断修正 + 局所収縮率 g + cp8 nSub 直接 A/B (Codex 補正 3 点への対応)**:
  - 診断修正: `FORGE_OUT_RESIDUALS` の dq 出力を **dq_block_old_*** (swap 後の最終補正) に修正
    (旧 dq_block_new_* は 1 sweep 前 — 定性結論は不変だが正した)。`FORGE_RESID_SNAP=1` を追加
    (subiter0 の res/dq を res_*_m / dq_*_new スロットへ退避し出力 → g を場で測れる)。
  - **局所収縮率 g=|dq_final|/|dq_sub0| (cp8, step100, run_diag_cflmap_cp8_g)**: 中央値 0.000-0.004、
    p99 ≤0.039、**増幅 (g>1) ノードはゼロ**。帯別も壁ノード帯 0.010 vs コア 0.004 の緩い勾配のみ。
    ⇒ **局所的な数値律速は存在しない — 収縮率は全域一様**。判定木は「off-line 枝」でなく
    **「全域一様 → 1次 FVS LHS × 2次 KEEP RHS の演算子不整合が本命」** に確定 (先の off-line 判定は
    η 最大位置の誤読で撤回済み)。streamwise 第 2 ライン族は根拠を失い保留。
  - **cp8 × nSub 直接 A/B (300 step, 同一 IC)**: cp8+nSub13 = **315 s (ctrl の 0.91 倍) で
    全指標 ctrl 同等以上** (ro 6.6e-10 / roUx 3.3e-7 / roe 2.13e-4, ωバースト 6.02)。
    cp8+nSub12 = 291 s (0.84 倍) だが roe 3.0e-4 (+42%) の妥協。cp4+nSub13 とはほぼ等価。
  - **生産推奨 (最終)**: `lineImplicit1 + lineKFreeze1 + lineDtDirectional1 + cfl_pseudo 8 + nSub 13`
    (同品質 0.91×・バースト −44%)。速度優先なら nSub12 (0.84×, roe 微妥協)。
    さらなる短縮は演算子不整合ファミリー (JFNK / LHS 散逸整合) か adaptive early-exit。
- 2026-09-03: **訂正 (Codex レビュー): g 解析のマスク欠陥 — 「全域一様・FVS-KEEP 確定」を撤回**。
  先の「増幅ゼロ」は分母フィルタ (ini>p50) が「初期 dq が小さく後から成長する点」を定義から
  除外していた。マスクなし再計算 (`analyze_g.py`, マスク明記の再現スクリプトとして case に残置):
  **roe g>1 = 15,242 ノード / roUy 3,972** (最終側マスク |dq_fin|>1%max でも roe 7,613 / roUy 1,135)。
  成長点の空間分布が示唆的: **丘頂域 (x/h 0-1, 7-9) に集中し、roUy 成長点は壁至近
  (wd/h med 0.004)** — **cp14 発散モードと同じ (ro,roUy,roe) ファミリー・同じ場所**。
  ⇒ 支持される言明: 「subiter0 で活発な領域の大部分は良く収縮する。一方、遅れて補正が成長する
  点が丘頂・近壁に実在し、その方向・原因は未判定」。off-line 確定でも FVS-KEEP 確定でもない。
  `FORGE_RESID_SNAP=<m>` を subiter 指定可能に拡張し、m=0/5/10/15/19 の履歴で成長点の
  軌跡を追う (run_diag_cflmap_cp8_g{,5,10,15})。
- 2026-09-03: **長回し検証の条件を訂正**: 最終採用を nSub13 とするなら**長回しも nSub13 で行う**
  (nSub15 の長回しでは nSub13 のバースト耐性は証明できない — run_0013 の nSub12 が step3839 で
  破綻した前歴があるため重要)。
- 2026-09-03: **subiter 軌跡解析 (m=0/5/10/15/19, run_diag_cflmap_cp8_g{,5,10,15}) — 「成長点」の正体は床**:
  同一 run 内の頑健な比較 (m0 と final は同ファイル) で、roe の g>1 集合は **subiter0 で対照の
  1/60 (8.6e-3 vs 5.1e-1) = 最初からほぼ収束済みの点**であり、subiter19 には両集合が**同じ dq 床
  (~1.6e-2) に合流**する (roUy も同型: 床 ~2.5e-5)。⇒ **cp8 に局所増幅モードは存在しない**。
  サブ反復は ~15 回で全域共通の dq 床に到達して止まる — これが収縮の実像で、nSub13 前後が
  最適という運用結果とも整合。**床の起源 (LHS/RHS 不整合のリミットサイクル / DDES 乱流の
  物理揺らぎ / relax 上限) は未分離** — ここから先は凍結場 (物理揺らぎ排除) での床測定などが要る。
  注意: run 間の最終 dq 再現性は中央値 60% 差 (深収束 dq はジッタ支配) — 集合統計のみ有効、
  点単位のクロス run 比較は不可 (`analyze_g.py` に明記)。
- 2026-09-03: **定常 M6 ノズルへの v2 適用可否 (run_0015_cflsweep/tp_cfl6r07_line2 等, warm 2000 step)**:
  **不可 — v2 でも定常側は損のまま**。① line v2 (directional なし): 末尾残差 −10% だが
  **step 単価 1.92 倍** (15.9→30.6 s) — 2D 小メッシュはライン 1250 本=1250 スレッドで占有率が
  壊滅し、DDES (9821 本, 1.32×) と逆に v1 比でも重い。② **directional は cfl6/8 とも発散 @step25**:
  壁法線音響で絞られていた BL セルの小 Δτ が実は安定マージンだった (dual-time の BDF 対角の
  ような保護が定常には無い)。⇒ **v2 は dual-time DDES 専用の武器**と結論。定常 RANS の生産は
  従来どおり point-DPLUR + cfl6/8+relax0.7。
- 2026-09-03: **訂正 (Codex レビュー) + back-to-back 正式 A/B (同一バイナリ・同一 IC,
  run_0015_cflsweep/tp_cfl6r07_{point2,line2b,linemono,line_ni3,line_ni2,linedir2})**:
  前項の「v2 でも定常側は損」「DDES 専用」は言い過ぎだった。正式値:
  | 腕 (cfl6+r0.7, warm 2000 step) | Time | 末尾 roe |
  | point (ni5) | 15.99 s | 0.626 |
  | **split line (ni5)** | **28.97 s (1.81×)** | 0.556 (−11%) |
  | mono line (ni5, FORGE_LINE_MONO=1) | 43.90 s (2.75×) | 0.556 |
  | split line ni3 | 25.08 s (1.57×) | 0.577 |
  | split line ni2 | 発散 @785 | — |
  - **factor/solve 分離は定常でも有効** (2.75×→1.81×, −34%) — 「効かない」は誤り。ただし最良
    チューニング (ni3) でも 1.57× コスト vs roe −8% で**この case では採算不成立**。ni2 発散 =
    off-line lag に sweep ≥3 が必要。「1250 line の占有率壊滅」は未プロファイルの仮説に格下げ。
  - **directional 発散の帰属を訂正**: 種は **x/r_t 71-89 の下流域・全断面** (壁至近 3%・近軸 1%,
    res_nan_26.h5 保持) — 「BL セルの露出」でなく**既知の streamwise 内部モード**が、directional で
    Δτ の上がった下流域セルで先に点火したもの (cp6/cp8 同 step の理由も状態依存で説明がつく)。
    定常律速 = streamwise 対流、の既存結論を補強。
  - **結論の正確な形**: case/45 M6 定常の現設定では non-directional v2 は採算不成立・directional は
    使用不可。**他の定常高 AR 問題への一般化はしない** (3D 定常や別ケースは未検証)。
```

## 参考: `solver_density_cuda/cuda_forge/setDT_d.cu`

```
//#if defined(__CUDA_ARCH__)
//#pragma message "content of __CUDA_ARCH__: " __CUDA_ARCH__
//#endif

#include "setDT_d.cuh"
#include "lowMachPrecond_d.cuh"   // Phase 4: β・c' で Δτ' を前処理スペクトル半径に合わせて拡大


__global__ void setCFL_pln_d
( 
 flow_float dt,
 //flow_float dt_local,
 flow_float visc,
 flow_float* vis_turb  ,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 //flow_float* cfl_pseudo ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln
 //flow_float* cfl_pseudo_pln

)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    if (ip < nPlanes) {

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        //__syncthreads();

        geom_float f = fx[ip];

        geom_float vol0 = vol[ic0];
        geom_float vol1 = vol[ic1];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];


        geom_float dx0 = vol0/sss;
        geom_float dx1 = vol1/sss;
        geom_float dx_min = min(dx0,dx1);

        flow_float Ux0 = Ux[ic0];
        flow_float Uy0 = Uy[ic0];
        flow_float Uz0 = Uz[ic0];

        flow_float Ux1 = Ux[ic1];
        flow_float Uy1 = Uy[ic1];
        flow_float Uz1 = Uz[ic1];

        flow_float US  = (f*Ux0 + (1.0f-f)*Ux1)*sxx
                        +(f*Uy0 + (1.0f-f)*Uy1)*syy
                        +(f*Uz0 + (1.0f-f)*Uz1)*szz;

        flow_float rof = f*ro[ic0] + (1.0f-f)*ro[ic1];
        flow_float v_turb = f*vis_turb[ic0] + (1.0f-f)*vis_turb[ic1];
        flow_float lambda = abs(US)/sss + sonic[ic0] + 2.0f*(visc+v_turb)/(rof*dx_min);

        cfl_pln[ip] = dt*lambda/dx_min;

        //cfl_pseudo_pln[ip] = dt_local*lambda/dx_min;
    }
}

__global__ void setCFL_cell_d
( 
 int dtControl, flow_float cfl_target, flow_float cfl_pseudo_target,
 int dualTime, int unsteady, 
 flow_float dt,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* cell_planes_index, geom_int* cell_planes,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 flow_float* dt_local  ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln,
 //flow_float* cfl_pseudo_pln

 // 軸対称 near-axis 安定化 (isAxisymmetric==1 かつ axisBeta>0 のときのみ作用)
 int isAxisymmetric,
 flow_float* A_planar,
 flow_float axisBeta,

 // v2 line-implicit 粘性 CFL 割引 (plans/active/time_integration-line-implicit-viscous-v2.md):
 // on-line セルに限り、面 λ の粘性項 2ν_eff/(ρ·dx_min) 由来の CFL を (1−θ) 倍に割引く。
 // 壁近傍セルの dx_min は壁法線 (=line 方向) なので方向整合は近似的に成立。θ=0 で完全不変。
 geom_int* plane_cells,
 flow_float visc_lam,
 flow_float* vis_turb,
 flow_float lineReliefTheta,
 const geom_int* line_prev,
 const geom_int* line_next,
 // 方向別 dt (lineDtDirectional==1): line 面 (Thomas が厳密に解く結合) の λ を CFL の max から
 // 完全除外し、Δτ を off-line 面 (lag 側) の λ だけで決める。壁セルの Δτ は streamwise 基準 (×AR)。
 int lineDtDirectional,
 // 診断 (lineDtWallRelief==1): wall 種境界半割面の λ も on-line セルの max から除外 (壁端点律速の切り分け)。
 const unsigned char* plane_wall,
 int lineDtWallRelief
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {

        geom_int index_st = cell_planes_index[ic];
        geom_int index_en = cell_planes_index[ic+1];
        geom_int np = index_en - index_st;

        cfl[ic] = 0.0f;
        //cfl_pseudo[ic] = 0.0;

        const bool onLine = (line_prev != nullptr && (line_prev[ic] >= 0 || line_next[ic] >= 0));
        const bool lineRelief = (lineReliefTheta > (flow_float)0.0 && onLine);

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            geom_int ip = cell_planes[ilp];

            if (lineDtWallRelief != 0 && onLine && plane_wall != nullptr && plane_wall[ip] != 0) continue;
            if (lineDtDirectional != 0 && onLine) {
                const geom_int ic0 = plane_cells[2*ip+0];
                const geom_int ic1 = plane_cells[2*ip+1];
                const geom_int other = (ic0 == ic) ? ic1 : ic0;
                if (other == line_prev[ic] || other == line_next[ic]) continue;  // line 面は除外
            }
            flow_float cfl_p = cfl_pln[ip];
            if (lineRelief) {
                // setCFL_pln_d と同一式で面の粘性 CFL を再計算し θ 分を引く (対流+音響分は必ず残る)。
                const geom_int ic0 = plane_cells[2*ip+0];
                const geom_int ic1 = plane_cells[2*ip+1];
                const flow_float f = fx[ip];
                const flow_float sss = ss[ip];
                const flow_float dx_min = min(vol[ic0]/sss, vol[ic1]/sss);
                const flow_float rof = f*ro[ic0] + ((flow_float)1.0-f)*ro[ic1];
                const flow_float v_turb = f*vis_turb[ic0] + ((flow_float)1.0-f)*vis_turb[ic1];
                const flow_float cfl_visc = dt * (flow_float)2.0 * (visc_lam + v_turb)
                                            / (rof * dx_min * dx_min);
                cfl_p = max(cfl_p - lineReliefTheta * cfl_visc, (flow_float)0.0);
            }
            cfl[ic] = max(cfl[ic], cfl_p);
        }

        // 軸対称 near-axis 半径音響スペクトル半径を加える: λ_axis = β·(|u_r|+c)·A_planar。
        // revolved 軸面積 (r_f·S→0) が落とす半径音響モードを planar 面積で補う。face 項と同じ
        // 無次元化 (dt·λ/V) で cfl へ加算 → 近軸で Δτ=cfl_pseudo·dt/cfl ∝ cfl_pseudo·r/(|u_r|+c)。
        // V=r·A_planar (per-radian) なので A_planar/V=1/r。LHS 時間項 v/Δτ が全保存式で自動的に強まる。
        if (isAxisymmetric == 1 && axisBeta > (flow_float)0.0) {
            const flow_float u_r      = fabsf(Uy[ic]);                                 // 半径方向速度 (axisym: y=radial)
            const flow_float lam_axis = axisBeta * (u_r + sonic[ic]) * A_planar[ic];   // [V/time]
            cfl[ic] += dt * lam_axis / max(vol[ic], (flow_float)1.0e-30);
        }

        //if (dualTime == 1 or unsteady == 0) {
        //    dt_local[ic] = cfl_pseudo_target*dt/cfl[ic];
        //} else {
        //    dt_local[ic] = dt;
        //}
    }
}

__global__ void setDTlocal_pseudo_cell_d
( 
 flow_float cfl_pseudo_target,
 flow_float dt,

 // mesh structure
 geom_int nCells,
 flow_float* cfl ,
 flow_float* dt_local
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float local_cfl_rate = cfl[ic] / max(dt, static_cast<flow_float>(1.0e-30));
        dt_local[ic] = cfl_pseudo_target / max(local_cfl_rate, static_cast<flow_float>(1.0e-30));
    }
}

__global__ void setDTlocal_uniform_cell_d
( 
 flow_float dt,
 geom_int nCells,

 flow_float* dt_local

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        dt_local[ic] = dt;
    }
}

// Phase 4 (lowMachPrecond>=2): 擬似時間刻みを前処理スペクトル半径に合わせて拡大する。
// 既存 dt_local は物理スペクトル半径 λ_phys=|u|+c 基準。前処理後の律速は
// λ'=½(1+β)|u|+c' (低マッハで小) なので、dt_local ×= λ_phys/λ' で擬似 CFL を λ' 基準に合わせる。
// ε フロアにより λ'≥ε·c 程度に下限が付き、倍率は ~1/ε で有界 (発散しない)。
// LHS を前処理した本モードでのみ整合的 (Phase 2 で LHS 非前処理のまま setDT 前処理して破綻した轍を踏まない)。
__global__ void setDTlocal_precond_scale_d
(
 flow_float precondEps,
 geom_int nCells,
 flow_float* ro,
 flow_float* Ux, flow_float* Uy, flow_float* Uz,
 flow_float* sonic,
 flow_float* dt_local
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float c = max(sonic[ic], static_cast<flow_float>(1.0e-8));
        const flow_float vmag = sqrt(Ux[ic]*Ux[ic] + Uy[ic]*Uy[ic] + Uz[ic]*Uz[ic]);
        const flow_float beta = lowMachBeta(c, vmag, precondEps);
        const flow_float cprime = lowMachCprime(c, vmag, vmag, precondEps);
        const flow_float lam_phys = vmag + c;
        const flow_float lam_prec = static_cast<flow_float>(0.5)*(static_cast<flow_float>(1.0)+beta)*vmag + cprime;
        const flow_float ratio = lam_phys / max(lam_prec, static_cast<flow_float>(1.0e-30));
        dt_local[ic] *= max(ratio, static_cast<flow_float>(1.0));
    }
}


// 診断: 出口近傍の局所 dt キャップ (壁∩出口コーナー不安定の切り分け用, env ゲート)。
//   FORGE_DT_OUTLET_SCALE (0<s<1) と FORGE_DT_OUTLET_XMIN [m] の両方を与えると、
//   x > XMIN のセルの dt_local を s 倍に縮める。既定 (未設定) は完全に不変。
__global__ void scaleOutletDt_d(geom_int nCells, const flow_float* ccx,
                                flow_float xmin, flow_float scale, flow_float* dt_local)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells && ccx[ic] > xmin) dt_local[ic] *= scale;
}

void setDT_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , bool adaptDt , bool printCfl)
{
    setCFL_pln_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        cfg.dt,
        //cfg.dt_local,
        cfg.visc,
        var.c_d["vis_turb"] ,

        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["cfl"]  , 
        //var.c_d["cfl_pseudo"] , 
        var.c_d["sonic"]  , 
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  ,

        var.p_d["cfl_pln"]
        //var.p_d["cfl_pseudo_pln"]
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    // 注: 同一 default stream なので後続カーネルは順序保証される。per-step 同期削減のため
    // 中間の cudaDeviceSynchronize は撤去し peek のみとする (誤差チェックは末尾/次同期点で担保)。

    setCFL_cell_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
        cfg.dtControl, cfg.cfl, cfg.cfl_pseudo,
        cfg.dualTime, cfg.unsteady,
        cfg.dt,
        //cfg.dt_local,

        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes, msh.map_cell_planes_index_d , msh.map_cell_planes_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["cfl"]  , 
        var.c_d["dt_local"] ,
        var.c_d["sonic"]  ,
        var.c_d["Ux"]  ,
        var.c_d["Uy"]  ,
        var.c_d["Uz"]  ,

        var.p_d["cfl_pln"] ,
        //var.p_d["cfl_pseudo_pln"]

        // 軸対称 near-axis 安定化
        cfg.isAxisymmetric,
        (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"],
        cfg.axisTimestepBeta,

        // v2 line-implicit 粘性 CFL 割引 (lineViscousDtRelief=0 で完全不変)
        msh.map_plane_cells_d,
        cfg.visc,
        var.c_d["vis_turb"],
        (cfg.lineImplicit == 1) ? cfg.lineViscousDtRelief : (flow_float)0.0,
        (cfg.lineImplicit == 1) ? msh.line_prev_d : nullptr,
        (cfg.lineImplicit == 1) ? msh.line_next_d : nullptr,
        (cfg.lineImplicit == 1) ? cfg.lineDtDirectional : 0,
        (cfg.lineImplicit == 1) ? msh.plane_wall_flag_d : nullptr,
        0   /* 壁境界半割面の λ 除外は不採用 (発散した診断スイッチ) */
    ) ;
    gpuErrchk( cudaPeekAtLastError() );

    if (cfg.dualTime == 1 or cfg.unsteady == 0) {
        setDTlocal_pseudo_cell_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            cfg.cfl_pseudo,
            cfg.dt,
            msh.nCells,
            var.c_d["cfl"],
            var.c_d["dt_local"]
        );

    } else {
        setDTlocal_uniform_cell_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            cfg.dt,
            msh.nCells,
            var.c_d["dt_local"]
        );
    }

    // Phase 4: 完全前処理モード (LHS) では擬似時間刻みを前処理スペクトル半径基準に拡大する。
    // lowMachPrecond==2 (RHS+LHS) / ==3 (LHS のみ) の双方で擬似時間項は前処理されるため両方で適用する。
    if (cfg.lowMachPrecond >= 2) {
        setDTlocal_precond_scale_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            cfg.precondEps,
            msh.nCells,
            var.c_d["ro"],
            var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
            var.c_d["sonic"],
            var.c_d["dt_local"]
        );
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
    }

    // max cfl の host 読み出し (thrust::max_element は D2H + 同期)。
    //  - printCfl: モニタ行表示用に cfg.monitorCflMax へ格納するため host read が必要 (呼び出し側が
    //    unsteady==1 かつ monitor step のときだけ true にする)。
    //  - dt 適応 (adaptDt && dtControl==1): 原理的には device 上で完結できる (cfl_max も cfg.dt も device 化すれば
    //    host 同期不要) が、現状 cfg.dt が host スカラ (多数カーネルに値渡し・dual-time BDF 係数等で host 使用) の
    //    ため host で計算している。よってこの条件は「現実装の都合」であり本質的要請ではない。
    //    device-resident dt 化は explicit/dual-time の毎ステップ適応 (dtControl==1) で per-step 同期を消せるが、
    //    定常 implicit は dt_local=cfl_pseudo·dx/λ で cfg.dt が打ち消され不影響なので利得なし (別 plan 候補)。
    // 両条件とも false なら host 同期は一切発生しない (後続カーネルは同一 default stream で順序保証)。
    // 出口近傍 dt キャップ (診断, env ゲート): dt_local 確定後に縮める。既定 (未設定) は不変。
    {
        static const double dtOutletScale = [](){ const char* e = getenv("FORGE_DT_OUTLET_SCALE"); return e ? atof(e) : 0.0; }();
        static const double dtOutletXmin  = [](){ const char* e = getenv("FORGE_DT_OUTLET_XMIN");  return e ? atof(e) : 0.0; }();
        if (dtOutletScale > 0.0f && dtOutletScale < 1.0f) {
            scaleOutletDt_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                msh.nCells, var.c_d["ccx"], (flow_float)dtOutletXmin, (flow_float)dtOutletScale, var.c_d["dt_local"]);
            gpuErrchk( cudaPeekAtLastError() );
        }
    }

    const bool needHostRead = (adaptDt && cfg.dtControl == 1) || printCfl;
    if (needHostRead) {
        thrust::device_ptr<flow_float> d_ptr = thrust::device_pointer_cast(var.c_d["cfl"]);
        const flow_float cfl_max = *(thrust::max_element(d_ptr, d_ptr + msh.nCells));

        // 印字は console モニタ行 (main.cpp StepMonitor) が担う。cfl_max とそれを評価した dt を対で格納する
        // (下の適応で cfg.dt が変わっても表示の対応が崩れないように、適応前に取る)。
        if (printCfl) {
            cfg.monitorCflMax = cfl_max;
            cfg.monitorCflDt  = cfg.dt;
        }

        if (adaptDt && cfg.dtControl == 1) { // cfl based time control
            flow_float cfl_target = cfg.cfl;
            cfg.dt = cfg.dt*cfl_target/cfl_max;

            cfg.dt = max(cfg.dt, cfg.dt_min);
            cfg.dt = min(cfg.dt, cfg.dt_max);
        }

    }

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
