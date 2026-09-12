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
4. **CHT が要るのはいつか**: $T_w(x)$ が未知で、しかも δ\*・$q_w$ が $T_w(x)$ の分布形に敏感なとき (薄肉・再生冷却・局所ホットスポット)。その場合も**フル CHT (固体伝導ソルバの連成) の前に「弱 CHT ループ」** — NS の $q_w(x)$ → 1D 壁伝導 + 冷却剤熱伝達モデル → $T_w(x)$ → `wallProfile` で NS 再実行 — を推奨する。等温壁機構と積分法の `Tw_table` がそのまま使え、固体側は解析式なので実装コストが小さい。フル CHT は弱ループが収束しない (壁内の軸方向伝導が支配的) と分かってから別 plan。
5. **台帳の標準項目** (全チェーン共通): `wall_thermal` / 実測 y₁⁺ (max) / $Q_w$ / $q_w$ ピーク位置 / δ\*_exit / 目的量の断熱比。

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
| 3 | S3 ノズル CPG × SU2 等温 | §4.4。**メッシュ確定 (2026-09-12, `case/44/mesh_probe_iso.py`)**: 冷却壁の y₁⁺ 倍率は run_0107 の壁データから 5.0–5.4 (スロート y₁⁺ 7→37)。`mesh2d` に x 依存の第一セル (`wall_first_frac_throat` + 上下流ブレンド) を追加し、**ni 4001 / nj 113 / wff 1.8e-5 / wfft 2.5e-6 / throat_refine 30 で AR max 846 PASS、y₁⁺_cold 1.5–4.5 (452k 節点)**。AR ≤ 1000 を守ると ni ≥ 3600 が要る (y₁⁺ 1 にするには ni 20000 級で不可: §8-3)。対象は case/45 でなく **case/44 va (Tt 1060, M4.19)** に変更 (SU2 の規模; 断熱 SU2 も同メッシュで取り直す = codex M7) |
| 4 | S4a 固定形状の壁温感度 (同一壁・同一メッシュ, 断熱 vs 300 K) → S4b 300 K 再設計 | §4.5 (codex M4)。台帳は感度と回復量を別項目。**投入 (2026-09-12)**: `case/44/run_iso_chain.py` → run_0108 (断熱, 細分メッシュ, run_0107 と同じ壁 `delta_r_initial.csv`) → run_0109 (300 K 同壁) → run_0110 (300 K 積分法初期壁 pass 0, 符号付き δ\*)。**run_0108 完了**: 品質 PASS (AR 846)、実測 y₁⁺ max 2.8 (上流収縮部) / スロート 0.61 / 出口 0.39 (断熱)、ṁ_NS/ṁ_E 0.9988 (run_0107 0.9992)、出口コア M 4.1869 (run_0107 4.1900, −0.07 % = メッシュ感度)、出口 T0 1060.0 K / P0 1.139 MPa (h0 逆算)、`wall_thermal_ledger.py --series` STEADY (M_core drift 0.000 %)。台帳ツール `case/44/wall_thermal_ledger.py` (y₁⁺(x)・q_w(x)・Q_w・出口 T0/P0・ṁ 比・series)。**run_0109 (S4a-B, 300 K 同壁) 完了 (24000 step, 残差 still converging, series: M_core drift 0.07 %・Q_w 9.6 % = 熱場が未静定 → run_0114 で +36000 継続)**。**固定形状の壁温感度 (暫定, run_0108 → run_0109)**: 出口コア M 4.1869 → **4.2071 (+0.48 %)**、ṁ_NS/ṁ_E 0.9988 → 1.0011 (+0.23 %)、出口 T0 コア 1060 K 不変・質量平均 1060 → 1048 K (−1.2 %, BL の熱損失)、P0 コア不変、**Q_w 7.8 MW、q_w ピーク 2.83 MW/m² @ x −0.41 r_t (スロート直前)**、実測 y₁⁺ スロート 3.0 / 出口 1.8 / 上流収縮部 max 10.8。抽出 δ_r (同壁) は冷却で 0.6〜0.75 倍 (出口 0.091 vs 0.127 r_t)、スロートで ≈0 (−2e-5)。**壁温は風洞の出口 M を設計公差 (±0.1 %) の 5 倍動かす → §4.6-2 の判定は「公称 T_w で設計」側**。run_0110 (S4b pass 0): 積分法初期壁 (Tw 300) は δ_r 出口 0.115 (断熱 0.122, −6 %) と冷却効果を過小評価 (CFD 抽出は −29 %) → pass 1 (run_0115, 符号付き抽出) で補正する設計どおり |
| 5 | S5 docs・CAPABILITIES | §5-5 |
| 6 | `wallProfile` CSV ($T_w(x)$ 分布の NS 入力) | 弱 CHT ループ (§4.6-4) の前提。`inletProfile` と同型 (bcond floats `wallProfile: path`)。別 plan 候補 |
| 7 | 弱 CHT ループ (1D 壁伝導 + 冷却剤モデル) | §4.6-4。#6 の後、別 plan |
| 8 | SERN (⑤) の等温評価 | 配管は #1 で共通化。検証は SERN plan R1–R7 の後 |

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
- **判定に使う run パスは case README の run 一覧に同期** (AGENTS.md)。

### 6.1 レビュー記録 (codex)

[`AGENTS.md`](../../AGENTS.md) 「codex レビュー」の記録表 (`solver_density_cuda/tools/codex_review.py`)。

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-12 | [2026-09-12-tooling-nozzle-isothermal-wall-chain-plan.md](../../notes/reviews/2026-09-12-tooling-nozzle-isothermal-wall-chain-plan.md) | GO-with-changes, C0/M7/m2 | **全件採用**: M1 符号付き δ\* → §4.1 + §5.1-1a / M2 低 Re $q_w$ 後処理 + $ds$ 重み + 閉合 → §4.1 + §5.1-1b / M3 定量ゲート (`--quantity cf`・`--series` の誤記も訂正) → §4.3-6 + §5.1-1c / M4 固定形状感度と再設計の分離 → §4.5 (S4a/S4b) / M5 温度–速度関係は診断 + Duan–Martín + 平板積分入口 → §4.3-3 / M6 引き継ぎ 3 種 + 起動ゲート → §4.4 / M7 A-plain 追加 + SST 条件表 + nx 感度 → §4.3 / m1 `--init-thermal` 廃止 → §4.1 / m2 y⁺ ×5.5〜6 → §4.2。スポット検証 (codex): `run_0048` の `qwall` 全点 0、`run_0013`/`run_0107`/`run_0048` は既定条件で NOT CONVERGED (stalled) — 本計画の warm 床条件はこれを踏まえて定量化 |

## 7. 影響範囲

- `design/forge_design/{probdef.py, evaluate/runner_wt.py, evaluate/runner_sern.py, evaluate/runner_axismach.py, feedback/deltastar_loop.py}`
- 新規 `case/48.flat_plate_cooled_m4/`、case/45・case/44 の問題 YAML 変種と run
- docs: `methods/design/overview.md` (新節)、`procedures/recommended-settings.md` (壁行)、`procedures/verification/README.md` + 新ファイル、`design/CAPABILITIES.md`
- forge 本体は変更しない (等温壁・低 Re SST は既存)。変更が要ると分かった場合 (例: 冷却壁での EOS 床・ピンの不具合) は本 plan に §4.7 として追記してから触る

## 8. 未確定事項 (ユーザ確認)

1. **公称 $T_w$**: 初版は 300 K (ユーザ発言 2026-09-12「壁面温度 300 K くらいでまずは」)。実機の冷却方式が決まったら $T_w^{\rm nom}$ を差し替える。
2. **CHT の要否**: §4.6-4 の判断基準 (弱 CHT ループが収束するか) で決める。フル CHT の実装は本計画の外。
3. **AR ゲートと y₁⁺**: 冷却壁で y₁⁺ ≤ 上限を守ると AR が 1000 を超える形状があり得る。その場合 `ni` を増やす (計算時間 ×数倍) か、AR 上限を「壁法線方向の構造格子は AR 2000 まで可」に緩めるかは S2/S3 の結果でユーザ判断。

## 9. 完了条件

- [ ] `methods/design/overview.md` の新節を実装と同期
- [ ] S1–S5 完了、§6 の判定を満たす
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映
- [ ] `status: done`、§10 変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 10. 変更ログ

- `2026-09-12` — S1 配管・符号付き δ\*・平板評価ツール実装、case/48 run A/B/C 完了 (§5.1-2 に数値)。y₁⁺ の冷却倍率は実測 ×5.7 (3 µm: 断熱 0.08 → 300 K 0.46)。
- `2026-09-12` — codex plan 段レビュー (GO-with-changes, C0/M7/m2) を §6.1 に記録、**全件採用**して §4.1/§4.2/§4.3/§4.4/§4.5/§5.1/§6 を改訂 (符号付き δ\*、低 Re $q_w$ 後処理と閉合、定量ゲート、S4a/S4b 分離、SST 条件表、引き継ぎ 3 種)。case/48 の起動レシピ実測 (層流暖機 + 2 次 cfl ランプ、本段 cfl 2) を §5.1-2 に記録。
- `2026-09-12` — 起票 (ユーザ依頼: ノズル設計の壁を断熱から等温/CHT へ。等温 low-Re SST を平板で理論・経験式と SU2 に対して検証 → ノズルで SU2 と比較 → 等温 NS で δ\* 反復 → 壁温影響の評価方針)。ユーザ判断 2 件を反映: 壁関数は使わない (準備状況が悪い)、冷却で y⁺ が変わる点を平板の y₁⁺ 掃引で扱う。
