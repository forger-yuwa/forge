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

## ブリーフ (`notes/reviews/briefs/2026-10-12-bell-m-deltastar-extraction-acceptance.md`)

# 諮問ブリーフ: ③ベル (M) の δ* 抽出の受理 — 帯位置の感度の線が不受理、較正への進め方

- 日付: 2026-10-12
- 依頼元: 主セッション (Opus)。エスカレーション条件 3 (事前に書いた受理線との比較が FAIL)・1。諮問は codex (diagnose) — ユーザ指示
- 作業ツリー: `/home/sano/work/forge-thruster` (HEAD `f9e4d6a7`)
- 対象: plan `plans/active/tooling-nozzle-bell-two-mode.md` §5.1 #7 の「(7b) 較正の設計 = codex (diagnose) 2026-10-12」(受理線: 抽出不確かさ ≤ 5 %、L70 の出口の誤差 ≤ 10 %・補正区間の相対 RMS ≤ 10 %・局所 ≤ 20 %)。前回の記録 `notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md`

## 1. 実施したこと (`f9e4d6a7`、試験 `design/tests/run_bell_deltastar_tests.py` 31 件 PASS、①の既定の抽出・contur_v1/v2 の既定はビット同一)

- 抽出を一定 x の断面に (`metrics/bell_deltastar.py` の `extract_bell`、`remap_constant_x=True` 既定)。流量の帳簿も一定 x で取り直し (NS 268.336 / Euler 268.870 kg/s、比 0.998015、断面間の散らばり 0.03 %。列のままの帳簿は 0.96769・散らばり 1.8 % で③では使えない)。
- 帯位置の感度: 実際に選ばれた帯位置 y_b を倍率 (0.75・1.5・2.0) で固定して測り直す (`deltastar_from_core_matched_euler(yb_sens_factors=…)`)。抽出不確かさ u = 倍率ごとの主系列との相対差の最大。
- contur_v2 を③に接続 (`EdgeConditions(edge_M_full=…)` で縁の M を Euler run_0018 の壁の節点から全域、③の壁、RocketExhaust の MOC 気体、NS と同じ μ(T)、断熱、θ0 は仮想平板 x_v = L_pipe·r_t = 0.14 m)。
- δ_r の時系列の台本 (`case/68.bell_two_mode/rao_m/run_bell_dstar.py series`) — まだ回していない (drift・osc の閾値が未登録)。
- L70 は台本だけ作り、値は読まずに出力を消した (検証の盲検性のため)。

## 2. 観測事実 (出力ファイル: `case/68.bell_two_mode/rao_m/dstar_l100_uncert_bands.csv`・`contur_l100_ratio.csv`・`dstar_l100_ns_run0019_constx.json`)

**抽出不確かさ (L100、run_0019 res_60000 / run_0018 res_30000)**: 登録の判定 (x ≥ 0 で u ≤ 5 %) は**不受理** — x ≥ 0 の 789 断面すべてで超え、最大 62 % (x 10.89)。

| x の区間 | n | u の最大 | 0.75 倍の差 | 外向き (1.5・2.0) の最大 |
|---|---|---|---|---|
| [0, 0.5) | 75 | 0.117 | −11.7〜−5.5 % | 9.1 % |
| [0.5, 1) | 73 | 0.168 | −13.5〜−7.8 % | 16.8 % |
| [1, 2) | 123 | 0.158 | −15.8〜−7.4 % | 4.1 % |
| [2, 4) | 140 | 0.556 | −55.6〜−12.4 % | 3.7 % |
| [4, 8) | 97 | 0.557 | −55.7〜−12.2 % | 3.1 % |
| [8, 16) | 153 | 0.620 | −62.0〜−7.7 % | 3.2 % |
| [16, 出口] | 128 | 0.604 | −60.4〜−8.5 % | 2.5 % |

- 外向きだけなら 5 % 超は 89 断面で、すべて x < 0.891。x ≥ 0.89 では最大 4.1 %。
- 選ばれた帯位置 y_b/r_w: x 0.5〜2 で 0.02 (下限 `band_min_frac`)、4 で 0.025、8 で 0.031、16 で 0.049、出口で 0.061。
- 実装者の仮説 (未検証): adaptive の帯は「最初に条件を満たす位置」で止まるので 0.75 倍は構造上ほぼ境界層の中に入り、内向きの差は抽出の曖昧さでなく境界層そのものを測っている。x 0.25〜1 は帯が下限に張り付き、外へ動かすと +7〜17 % 増える — 主系列がこの区間で過小の可能性 (前回の第 2 仮説の再現)。帯を 2 倍に固定した系列で比をとると x ≥ 1 で中央値 1.0375・範囲 1.012〜1.083。

**contur_v2 (未較正) と NS の比 δ_r,NS/δ_r,CONTUR**: x 0.5 で 0.987、1 で 1.040、2 で 1.043、4 で 1.042、8 で 1.035、16 で 1.022、出口で 1.014。x ≥ 1 で中央値 1.0285、最小 0.9901・最大 1.0435、相対標準偏差 1.25 %、最良の定数 1.0279 で相対 RMS 1.26 %・局所最大 3.8 %。x 1〜4 の 1.04 から出口の 1.014 へ単調に下がる。CONTUR の感度: 縁の上流を 1D にすると x ≥ 1 で 1.6 %、θ0 × 4 で 0.7 %。

**その他**: δ_r の定常性 (M3) は未判定 (時系列の閾値が未登録)。回復係数の Pr は 0.72 のまま (燃焼室の凍結 Pr ≈ 0.709)。

## 3. 問い (推奨を 1 つずつ)

1. **抽出不確かさの定義**: 登録の倍率 (0.75・1.5・2.0) の内向き 0.75 は「境界層の中に入る」ので不確かさの尺度として不適切か。外向きの倍率 (例 1.5・2.0、帯が境界層の外にあることの確認) だけで測るのが正しい定義か。それとも帯の選び方 (adaptive → edge など) を変えるべきか。結果を見てから線を変えることにならないよう、どう扱うか (定義の欠陥としての登録の訂正か、不受理のまま扱いを変えるか)。
2. **x < 0.9 の区間** (帯が下限に張り付き外向きで +7〜17 %): 較正区間から外す (例 x ≥ 1) か、下限を変えるか。
3. **較正 kδ への進め方**: 未較正で比 1.03 ± 1.3 % (x ≥ 1)。定数 1 個の較正で受理線 (L70 で出口 ≤ 10 %・RMS ≤ 10 %・局所 ≤ 20 %) は十分に余裕があるが、抽出の主系列そのものが帯の位置で数 % 動く (2 倍固定で中央値 1.0375)。kδ の値と不確かさをどう決めるか。
4. **δ_r の時系列 (M3)**: drift・osc の閾値と判定の対象 (代表 x 数点と出口、全体の L2) の事前登録。
5. 次の一手 1 つとその判定 (L70 の盲検の検証を開ける条件を含めて)。

## 4. 読んでよいもの

- `plans/active/tooling-nozzle-bell-two-mode.md` (§4.4・§5.1 #7)
- `notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md` (前回)
- `design/forge_design/metrics/bell_deltastar.py`、`design/forge_design/metrics/deltastar.py` (帯の選択 232 行付近、`yb_sens_factors` 634・660 行付近)
- `case/68.bell_two_mode/rao_m/dstar_l100_uncert_bands.csv`・`contur_l100_ratio.csv`・`dstar_l100_ns_run0019_constx.json`
- python で抽出を呼び直してよい (手元に run_0019・run_0018・run_0020・run_0023 の最終場がある — **L70 [run_0020・run_0023] の δ_r の値は読まない**、検証の盲検性のため)。ファイルは書かない。禁止: `*.log`, `*.pdf`, `plans/README.md`

## 関連 plan 全文 (`plans/active/tooling-nozzle-bell-two-mode.md`)

```markdown
# ③ベル: 特性曲線で辿るモード (TP の Rao) と形状を直接変えるモード (C2 B-spline) の比較

## メタ

- **area**: `tooling / optimization`
- **status**: `draft`  <!-- 2026-10-11 起票。同日 codex diagnose NO-GO → 改訂 → plan 段 GO-with-changes (C0/M8/m1) → 全件採用で改訂。段階方式: 共通評価器の受理 → (D) → S1 ゲートを通した (M) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (実装着手前に「ベルの 2 モード」節を足す — §5 ステップ 0)
  - [`notes/investigations/nozzle-top-internal-shock-diagnosis.md`](../../notes/investigations/nozzle-top-internal-shock-diagnosis.md) (TOP の内部衝撃波。Rd を変えても強さはほぼ不変)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) §3 確定事項 5〜11、§5.1 #3・#4
  - 兄弟: [`tooling-nozzle-ambient-separation.md`](tooling-nozzle-ambient-separation.md) (**別セッションの担当**、2026-10-11 ユーザ決定。本計画は外部領域と剥離評価が使える前提で進め、海面作動の比較 §6 V4 はその受理後。推力の会計は同じ定義)
  - 資産元: [`../accepted/tooling-nozzle-moo-loop.md`](../accepted/tooling-nozzle-moo-loop.md) (Phase 2)、[`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) §4.3 (平面の縁条件 `rao_planar.py`)・§4.5 (凍結 TP `gas/frozen.py`)、[`../accepted/tooling-nozzle-cfd-pinned-initial-line.md`](../accepted/tooling-nozzle-cfd-pinned-initial-line.md) (CFD ピンの初期線)
- **created**: `2026-10-11`
- **owner**: `sano`

## 1. 目的

③ベルの壁を決める**設計手順**として、(M) 特性曲線法で Rao の最大推力輪郭を作る手順と、(D) 壁の形を少数の幾何変数で直接振る手順を、
同じ評価器 (TP のロケット排気・同じ物理寸法・同じ格子規約・同じゲート・同じ推力定義) で比べる。比べるのは「完成した設計手順」で、
設計点 (真空・付着流) の $\eta$–$L$、形状が作れない点の率、較正・探索を含む総費用を比べ、外部領域の受理後は海面作動の剥離と内部衝撃波も比べる。
結果で、ベル (と④デュアルベルの基部) の設計モードを決める。

## 2. スコープ

- **やる**
  - ③の評価器を TP のロケット排気にする (基準ケース: LOX/LH2・燃焼室圧 10 MPa・O/F 6・面積比 $\varepsilon$ 50、CEA の燃焼室平衡組成を凍結、断熱壁)
  - 両モード共通の C2 の壁 (スロートの曲率は 1 つの値 $1/R_t$) と、共通の寸法決定 (§4.2)
  - (M) 凍結 TP の熱力学 ($h(T)$, $s(T,p)$) で導いた軸対称 Rao 条件と生成器、曲率連続の kernel、端点を任意に指定できる当てはめ、δ\* の一発補正
  - (D) 拡大部の 5 次 B-spline (内部の形を 2 変数で振る)
  - 比較の run、TOP (Phase 2) の参照点の再評価
- **やらない**
  - 海面作動の比較の実装 (兄弟 plan の外部領域を待つ。§6 V4 は登録だけ)
  - 冷却壁 (当面断熱、親 §3-10)、平衡流 (凍結のみ)
  - 壁静圧分布を指定する逆 MOC、TIC、軸マッハ分布 (M 側の追加の族として後で検討)

## 3. 関連 docs と前提

- **Phase 2 の実績**: TOP ($\theta_n,\theta_e,L$) の MOO で 77 評価、Rao 80 % のチャート点が同じ長さの最適 $\eta$ の 0.133 % 下 (判定帯 ±0.15 %)。固定長で $\theta_n$ 24.5〜28.5°・$\theta_e$ 8〜11° が平坦な台地。評価器は CPG・`viscMethod: 1` (空気の Sutherland)・SST・断熱壁・出口面で領域を打ち切り (`evaluate/runner.py:90,103,124`)。
- **既存の部品**: 軸対称・semi-perfect の MOC 単位過程 `geometry/moc_kernel.py`、Hall の遷音速級数 `geometry/transonic.py`、CFD ピンの初期線 `feedback/cfd_initial_line.py` (凍結源の形の照合は U→T 5 次 Hermite が前提)、位置と壁角の同時当てはめ `geometry/wall_axismach.py:999` `joint_fit_wall` (**始点の拘束は $r'=x_0/R,\ r''=1/R$ に固定** — `wall_axismach.py:1061`)、δ\* の積分法 `feedback/deltastar_integral.py` (粘性係数は空気の Sutherland、温度は回復係数の近似)、RANS からの δ\* 抽出 `metrics/deltastar.py`、凍結 TP `gas/frozen.py`、平面の縁条件 `geometry/rao_planar.py` (自前導出扱い — 同ファイル 12 行)。
- **輸送物性**: `viscMethod: 1` は空気の Sutherland。燃焼ガスには `viscMethod: 2` (種ごとの輸送物性、`physProp.transport` 必須、`thermalMethod: 2` のみ — [`procedures/solver-settings.md`](../../procedures/solver-settings.md) の `physProp.viscMethod`) を使う。
- **Rao の軸対称の条件**: Cain 2010 式 4.1–4.3 (軸対称形)、Shyne & Keith 1990 式 11–15 (平面版)。`papers/nozzle_design/` (git 管理外)。CPG の導出は TP の最適性の保証にならない (codex 2026-10-11 B2)。
- **内部衝撃波**: case/40 で TOP の Rd を 0.4→1.5 にしても内部の ΔM は −0.43→−0.40 でほぼ不変。C2 にしても内部衝撃波が消える保証は無い。

## 4. 設計方針

### 4.1 評価器 (両モード共通)

- **ガス**: CEA で燃焼室平衡組成を求めて凍結し、NASA-9 の擬似種 1 つ (`EXH`) の TP (`thermalMethod: 2`、`thermoHrefTemp: 298.15`)。輸送物性は `viscMethod: 2` で CEA 形 (lump は構成実種ごと)。全温は CEA の燃焼室温度。
- **壁**: 断熱、SST 低 Re ($y_1^+\le1$、`wallTreatmentSST: 0`)。SST の項は推奨レシピ。
- **TP 化の範囲** (plan レビュー M1): ソルバ設定だけでなく、推力抽出 (`metrics/extract.py:38` は CPG の圧力復元)、1D 理想推力と 1D チョーク流量 (`opt/driver.py:84` は定比熱)、IC (`evaluate/runner.py:197`)、壁の平滑化 (`runner.py:243`)、最適化ドライバのゲートを TP にする。同じ解決済みの組成・物性 DB・エンタルピー基準を MOC 側の `GasSemiPerfect` と CFD に渡す (`FrozenGas` には MOC が要る `nu`・`mach_of_nu` が無い)。CEA との照合は出口温度だけでなく $\dot m_{1D}$・出口速度・出口圧力・$F_{1D}$ まで行い、反応物の状態と凍結位置 (燃焼室) を固定する。
- **作動点 (第 1 段)**: 真空の設計点。領域は出口面までで、出口は超音速の流出 (`outflow`。現行 runner の `outlet_statPress: Ps=p_ambient` は残さない — `runner.py:126`)。 作動点は親 §3-14 と同じ「高度 (または $p_a$) と重みのリスト」で書き、真空は `p_a: 0` の要素とする。標準大気 (US 1976) は④と共用の `gas/atmosphere.py` (新規、所有はこの系列)。
- **推力** (plan レビュー M2): ノズル内壁・入口面・出口面で閉じた検査体積 CV の全境界を、外向き法線で $\Pi_x=\int[\rho u_x(\mathbf u\cdot\mathbf n)+(p-p_a)n_x-(\boldsymbol\tau_{eff}\mathbf n)_x]\,dA$ ($dA=2\pi r\,ds$) に統一し、定常で総和 0 になることから壁会計 (入口面 + 内壁) と出口面会計 (出口面、粘性応力を含む) を導く。両者の差を CV 残差とし、$F_{1D}$ で正規化する。壁の摩擦の符号は node/cell の違い (`metrics/sern_forces.py:73`) を確かめてから使う。兄弟 plan §4.4 と同じ定義で、本計画はこの定義を自前で持つ。$\eta_{vac}=F_{vac}/F_{1D,vac}(\varepsilon)$ ($F_{1D}$ は同じ凍結 TP で $\varepsilon$ まで等エントロピー膨張した 1D 推力)。
- **ユーザ決定 (2026-10-11)**: (1) 推力を壁会計と出口面会計の両方で出し、**2 つが大きく変わらないこと (CV 残差 $|F_{wall}-F_{exit}|/F_{1D}$ の絶対値が小さいこと) を収束判定の一つにする**。閾値と、残差の VERDICT・推力系の STEADY との組み合わせ方は diagnostician の回答を見て決める (主セッション案 1e-4)。(2) メッシュの AR は **5000 以下を基準**にする。AGENTS.md の例外 (壁法線に沿った、スキューの無い構造格子の境界層は ≤ 5000) の条件を満たすよう、近壁を壁法線の層にする方向で検討する (今の格子は半径方向の線が $x$ 一定で、拡大部の近壁セルにスキューがある)。以下のゲートの記述は、この 2 点を反映して改める。
- **ゲートの組み替え (2026-10-11 diagnostician、ユーザ決定 (1)(2) を反映、[記録](../../notes/reviews/2026-10-11-bell-2b-gates-diagnostician.md))** — 下の旧記述 (残差 PASS 必須) を置き換える。§6 V-gate の継続 run が A なら driver の定数をこれに変える:
  - (i) 残差: `check_convergence --segment` の VERDICT は台帳に記録。FAIL は DIVERGED と有意な末尾上昇だけ。`NOT CONVERGED (stalled)` は併記 (⑤ SERN §4.39 と同型)。基準ケースのプラトーは ω 残差の高周波ジッタと、遅い単調トランジェント (出口側壁の熱的発達・収縮部の μt) で、推力には 1e-4 以下しか効かない。
  - (ii) 定常: 壁会計・出口面会計・$\eta$ (出口面)・$\eta$ (壁) を drift・osc 1.5e-4 で、ṁ を 1e-3 で `check_quasisteady --series-csv` が STEADY。単調な列は漸近値 (最終比) も 1.5e-4 以内 (`classify` の実効の単調許容は 7.5e-5)。cv_residual は 0 を横切るので列から外す。本段の出力間隔は本段長の 1/20 以下 (末尾窓に 8 点以上)。
  - (iv) **収束 (ユーザ決定)**: 最終スナップショットの $|F_{wall}-F_{exit}|/F_{1D}\le1\times10^{-4}$。(ii) と組で使う (単一スナップショットの |cv| は両会計の交差点で偽合格する — run_0003 の −5.2e-7 は交差点で、閉じ代の床は両会計の漸近値の差 ≈ 5e-5)。旧 1e-3 の精度ゲートはこれに吸収。
  - (D) $\eta$ は出口面会計を主、壁会計を照合 (台帳には両方と cv)。
  - (v) 格子: 近壁を**壁法線の層**にして AGENTS.md の例外 (AR ≤ 5000) を使う (ユーザ決定 (2))。条件: 入口面・出口面は平面のまま (出口の数ステーションで半径線へ戻す)、AR > 1000 の全セルでスキュー ≤ 0.1 を自前で検査して README に「AR 緩和 (≤5000)」と書く (`check_mesh_quality.py --ar-max 5000` は全域一律なので条件を保証しない)。y1+ ≤ 1 にはスロート第一セル約 0.6 µm・収縮部約 2 µm (ni 1201 のまま AR ≈ 1600 / 1000〜1200)。y1+ 5.7 → 1.6 で $\eta_{exit}$ の差 6e-6・摩擦の差は $F_{1D}$ の 1e-4 — 新格子との差を足して格子感度の裏付けにする。
  - `space.pRef` は**置かない** (燃焼室 1e7 Pa に合わせると出口 1.4e4 Pa で相対誤差が約 700 倍。広い圧力域のノズルに代表静圧は無い — 推奨レシピ §5 からの逸脱の理由)。
- **ゲート** (旧記述、plan レビュー M3): 受理条件は、(i) 同じ実効設定の区間 (`StageManifest`) で `check_convergence` が **PASS** (現行 `opt/driver.py:173` の「DIVERGED でなければ可」は採らない)、(ii) 両推力会計・$\eta$・$\dot m$・CV 残差の時系列が `check_quasisteady.py --series-csv` で **STEADY**、(iii) Phase 2 の物理ゲート ($\dot m/\dot m_{1D}$、$\eta$ の範囲)、(iv) CV 残差 ≤ 0.1 %、(v) メッシュ品質と局所 $y_1^+$。polish は `restart_field.py` による**同一メッシュの継続**に替える (現行 `opt/polish.py:41` は `runner.prepare` 経由で補間し直すので継続ではない)。~~最初の TP の seed は段階起動、形状が変わる warm start は `interp_field.py` の cross-mesh restart。~~ (下の「起動」で置き換え)
- **起動 (2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-d-moo-warmstart-diagnostician.md))**: TP の評価は**全点を冷間の段階起動** (準 1D 等エントロピー IC → soft → mid → 本段)。cross-mesh の warm start は MOO では使わない。`interp_field.py` は生座標の最近傍なので (`interp_field.py:212-216`)、形状が変わると新旧の壁の間の帯に旧い壁の値 (u≈0・回復温度) が入り、本段で出口の壁の角から逆流が育って発散した (`case/68.bell_two_mode/run_0010_moo_d/doe_001` → 同じ形の冷間 `run_0011_tp_d_doe001_cold` は全ゲート PASS。帯の由来は仕様からの推論)。soft・mid の 4000 step は帯を洗えない。段の長さは同じなので費用は変わらず、評価が評価順に依存しなくなる。同一メッシュの継続 (polish) は従来どおり `restart_field.py`。
- **内部衝撃波の指標**: 診断ノートの方法 (各 $x$ 断面で Mach の最大降下の位置を追う) で最大の $|\Delta M|$ を記録する。
- **Euler の設定** (V2 で使う): `viscMethod: 0`、`visc: 0.0`、`thermCond: 0.0`、`turbulence.model: none`、物理壁はすべて `slip`。

### 4.2 共通の壁の構成と寸法 (C2)

- 入口直管 → 収縮部 (U→T 5 次 Hermite、スロートで $r''=1/R_t$。①の axis-Mach チェーンと同じ部品) → スロート → 拡大部。拡大部の作り方だけがモードで違う。
- **寸法の spec** (両モード共通、キャンペーン内で固定): スロート半径 $r_t$ (物理単位)、$R_t$、収縮比、収縮部長さ $L_U$ と入口直管長 $L_{pipe}$、**物理壁の** $\varepsilon$ と $L$。値は §8-1。収縮部の単調条件 $L_U^2/[R_t(r_U-r_t)]\le20$ (`geometry/wall_walldriven.py:99`) を満たすこと (収縮比 3 で $L_U/r_t\le3.826$ [$R_t=1$]、$\le2.420$ [$R_t=0.4$])。
- **最終の物理壁の拘束** (両モード、plan レビュー M6): $r(0)=r_t$、$r'(0)=0$、$r''(0)=1/R_t$、$r(L)=r_t\sqrt\varepsilon$、全域で最小半径がスロート、収縮部との C2 接続。
- 曲率は $\kappa=r''/(1+r'^2)^{3/2}$ と $r''$ を区別して扱う (スロートでは $r'=0$ なので一致)。

### 4.3 (D) 直接モード

拡大部 $r(x)$, $x\in[0,L]$ を 5 次 B-spline にし、係数は「拘束を満たし $\int(r''')^2dx$ を最小にする」線形問題で一意に決める。拘束は

- スロート: $r=1,\ r'=0,\ r''=1/R_t$
- 最大壁角の点 $x_m$: $r'=\tan\theta_{max},\ r''=0$
- 出口: $r=\sqrt\varepsilon,\ r'=\tan\theta_e$

**設計変数 = $(L,\ \theta_e,\ \theta_{max})$ の 3 個** (TOP の $\theta_n,\theta_e,L$ と同じ次元)。**$x_m$ は設計変数にせず内部で決める** (diagnostician 2026-10-11): 成立する $x_m$ は $L$ に比例せず絶対値 1.1〜2.1 $r_t$ の帯 ($R_t$ に縛られる) にあり、$x_m/L$ で振ると長い $L$ で帯を外す (V0-a で成立 8/135 の原因)。$x_m$ を 0.5〜4.0 $r_t$ で走査して成立区間を求め、その中で $\int(r''')^2$ 最小の $x_m$ を黄金分割で決める (決定的)。成立区間が空なら INFEASIBLE と理由を返す。成立条件: $r'\ge0$、壁角の極大が 1 つ、曲率の上限。破れたら INFEASIBLE。この壁がそのまま物理壁 (δ\* 補正なし)。拘束 7 個は、拘束行列が満行ランクならスロートの 3 拘束で目的関数の零空間 (2 次多項式) が除かれるので一意に解ける。**ノットの規約・係数の数・設計変数の範囲・曲率の上限・区間内の極値の検査・無次元座標**は §5.1 #3 で固定する (plan レビュー M8)。

### 4.4 (M) 特性曲線モード (TP の Rao)

- **初期線**: スロートの C⁻ 特性線。場の出所は **CFD ピン** (同じ収縮部・$R_t$・ガスの Euler 場を 1 本回して凍結) を既定、Hall は検算。凍結源の形の照合を本計画の収縮部と $R_t$ に合わせて広げる。
- **MOC の経路**: ①風洞の生産で検証済みの配列版 `interior_vec` (`axis_limit="analytic"`・`corrector="converge"`、対ごとの分類 `stats["status"]` を不受理の理由に) に乗せ、スカラー版 `KernelMOC` は使わない (2026-10-11 ユーザ指示「風洞のスクリプトを参考に」)。マーチの組み方は風洞の逆 MOC `geometry/moc_inverse.py` を手本にし、軸対称の検証には風洞の放射源流の試験 (`design/tests/moc_axis_limit_radial.py`) を使う。風洞のコードの既定は変えない。
- **kernel**: スロートから終点 B までの膨張部。曲率がスロートの $1/R_t$ から連続に変わる曲線族で、B での曲率 $\kappa_B$ をパラメータに持つ。古典の円弧 kernel も同じ生成器で選べるようにする (回帰試験用)。
- **Rao 条件 (TP)**: 凍結 TP の $h(T)$・$s(T,p)$ を使い、質量流量と長さの拘束付きで推力の変分をとって制御面 DE の条件を導く (S1、§5.1 #4)。比熱比を陽に含まない形に書けるかは導出の結果で判断し、仮定しない。**固定 $\varepsilon$ (出口半径の固定) が端点の変分をどう制限するか**も導出で定める。与えた $(L,\varepsilon)$ に対して B と D を求め、領域 BDE を MOC で埋め、質量流量で壁流線 BE を取る。
- **S1 の導出 (2026-10-11、diagnostician [記録](../../notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md))**: 凍結 TP でも制御面 DE の条件は古典と同じ — φ = θ + μ (DE は C⁺)、$\Phi_2=V\cos(\theta-\mu)/\cos\mu$ と $\Phi_3=r\rho V^2\sin^2\theta\tan\mu$ が DE 上で一定 — で、比熱比は陽に現れない (使う熱力学は $dp=-\rho V\,dV$・$d\rho=-\rho M^2 dV/V$ だけ)。**固定 $(L,\varepsilon)$ では縁条件 (Cain 4.3) は課さない** (端点項が消える。停留輪郭は $p_a$ に依らない。$F(E)$ は感度として報告)。D に追加条件は要らない。TP 固有の $k=d\ln a/d\ln V$ は C2 の条件と Goursat 領域の出だしにだけ入る。**構成**: kernel $\kappa(s)=1/R_t+(\kappa_B-1/R_t)S(s/s_B)$ (古典円弧も選べる) → kernel 場 (`forward_moc`) → **BD は B から下流の C⁻ を net 場の補間で追跡** (θ_B > μ_B では net の前線が逆向き) → DE を **r (弧長) をパラメータに** $\Phi_2,\Phi_3$ 一定で引き (θ = μ で M は極値をとるが線は正則に続く — 下の訂正)、流量が $\dot m_{BD}(s_D)$ に達する点を E → 入れ子の 1 次元根 (内 $s_D$: $r_E=\sqrt\varepsilon$、中 $\theta_B$: $z_E=L$、外 $\kappa_B$: C2) → 領域 BDE を Goursat で埋め、壁流線を質量流量で取る。**C2 の根の関数は $\tilde g(\kappa_B)=(d\theta/d\sigma)_{DE,D}-(d\theta/d\sigma)_{kernel\,C^+,D}$** (B の曲率連続と等価、流線の 2 階微分より安定)。上の $g(\kappa_B)$ は根の後の検算に使う。複数根は MOC 推力最大。~~折り返し・$|c|>1$・根が括れないは INFEASIBLE と理由。~~ **訂正 (2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-s1-generator-holes-diagnostician.md))**: DE 上で $c=\Phi_2\cos\mu/V=\cos(\theta-\mu)\le1$ は恒等的で、$c=1$ (θ = μ) は M をパラメータにした線の折り返しにすぎない (θ = μ では C⁺ の適合式と $\Phi_3$ 一定の両方が $d\theta/dr=-\tan\theta/(2r)$ を与え、線は正則に続く)。DE の解なしは **θ ≤ 0・θ+μ ≥ 90°・M が表の外 (TP は M_max 16)・M_cap まで流量に達しない** の 4 条件 (「dr ≤ 0」は C⁺ では θ+μ ≥ 90° と同値)。候補の受理条件に **$\theta_E\le\theta_B$** を足す (B の後で壁が外へ開いて戻る偽枝 — L_k 60 % 付近の θ_B ≳ 37°・θ_E ≈ 40° — を切る。「θ_E ≤ 20°」は (D) の設計箱に揃える比較の箱としてだけ使う)。`solve_thetaB` の枝は添字でなく **η_D の連続性**で追う。根が括れない L_k は「(M) で作れない L_k」として理由つきで記録する。
- **(M) は円弧 kernel・B で C1 (2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-rao-mainbranch-c2-diagnostician.md)、ユーザ了承)**: 上の「外 $\kappa_B$: C2」と下の「C2 接続」を撤回する。kernel の終端曲率 $\kappa_B$ を振る族には C2 の根が存在しない — $\kappa_B$ の影響は BD (約 6.7 $r_t$) を伝わる間にほぼ失われ (Δκ_B 0.6 で kernel 側の dθ/dσ が 0.006)、D での DE の横断微分は D の状態だけで決まる。部品試験 (`case/68.bell_two_mode/rao_m/probe_c2_l100.py`・`probe_c2_l70.py`、主セッションで再現) で主枝の根の $\tilde g$ は L100 で −0.069 (κ_B 1.00)・−0.067 (0.80)・−0.067 (0.59)・−0.063 (0.39)、L70 で −0.17〜−0.56、符号変化なし (事前登録の判定 A)。**ユーザ決定 (2026-10-11)**: C2 は目標のまま、解が無ければ (M) は古典どおり円弧 (C1) でよい — 推力用ノズルは Rao・TOP とも C1 (曲率の跳びあり) が標準。したがって (M) = 円弧 kernel (曲率 $1/R_t$ 一定)、B で接線連続、曲率の跳び $\Delta\kappa_B=\kappa_W(B^+)-1/R_t$ と $\tilde g$ を $L_k$ ごとに記録する。スロート・収縮部との C2 (§4.2) は不変。(D) は B-spline のまま C2。
- **主枝の受理規則 (同上)**: 受理条件は **det_D > 0** (DE の陰関数の行列式。det_D > 0 ⟺ DE が D で流れを軸向きに曲げ始める — 古典の枝の延長の物理的な定義)。上の θ_E ≤ θ_B は θ_B ≥ 40.8° の偽枝 (det_D < 0) を通すので**記録に格下げ** (de = dθ/dσ|_{DE,D} < 0 も記録)。det_D > 0 の候補が θ_E > θ_B になったら方針の穴として諮る。枝の結合は η_D の連続性かつ同じ det_D の符号。
- **S1 を独立の受理ゲートにする** (plan レビュー M5): kernel の族、$\kappa_B$ の探索域、B・D の解の枝の追跡、複数根の選択規則を定義し、各 $L_k$ で根と分解能依存性を確かめる。停留性の試験は摂動幅の掃引・MOC の細分化・二次変分の符号まで見る。証明できる範囲が限られる間は「指定した kernel と拘束の下の停留輪郭」と呼び、「最大推力」と呼ばない。S1 を通るまで (M) を比較キャンペーンに入れない。 **(2026-10-12: S1 を通過 — §5.1 #5 の札つき)**
- **スロート流量の定義**: CFD ピンの $m^*$ は CFD の保存量の積分ではなく等エントロピー写像の流量 (`methods/design/overview.md` の CFD ピン節)。V1 の「スロート流量」は $m^*$ と定義し、凍結源の CFD 流量との差を別に記録する。
- ~~**C2 接続**~~ (撤回、上の「(M) は円弧 kernel・B で C1」): 符号付きの関数 $g(\kappa_B)=\kappa_{Rao}(B^+;\kappa_B)-\kappa_B$ の根を求める ($\kappa_{Rao}(B^+)$ は kernel を $\kappa_B$ で作ったときの Rao 輪郭の出だしの曲率)。**各評価で B・D を解き直し**、$g$ の符号を挟めない候補は「C2 の成立未確認」として台帳に残し、評価に回さない。BE の MOC 壁点は、**始点の $r,r',r''$ と終点の $r,r'$ を任意に指定できる**よう `joint_fit_wall` を拡張した当てはめで B-spline にする。B を跨いで当てはめ直す方式は採らない。
- **δ\* 補正の扱い (ユーザ決定 2026-10-12、[記録](../../notes/reviews/2026-10-12-bell-m-deltastar-route-diagnostician.md))**: δ\* 補正は (M) の**設計者の選択肢** (補正する / しない)。補正の方法は **δ\* の理論式 (積分法など) を NS の計算結果で較正する形** (下の段落の向き — 積分法 + $k_f$ の較正 — で、較正の元を NS の結果にする)。補正しない場合は非粘性壁をそのまま物理壁として (D) と同じ NS recipe で評価する (停留輪郭なので推力の損は δ\* の 2 次、≈ 1e-4 η の見込み — diagnostician)。補正で狙いに近づくのは壁圧・壁沿いのマッハ数・出口の分布 (1 次)、推力・有効出口面積・摩擦損失は補正で大きくは変わらない (④・海面作動では壁圧が効く)。較正の設計 (理論式・較正の形・合格の線・出口半径とスロートの扱い) は η_M0 の NS 1 本 (ピン版 L100 の非粘性壁) の結果を添えて上位に諮ってから。
- **δ\* と物理寸法** (plan レビュー M6・M7): 積分法を TP にする — 粘性係数を評価器と同じ $\mu(T)$ に、運動量厚さの式の速度勾配を TP の縁速度から得た $d\ln u_e/dx$ に (現行 `feedback/deltastar_integral.py:208` と 116 行は CPG の換算で、単一 H2O・$T_t$ 3500 K・$T$ 1000 K の合成例で −6 % ずれる — codex の検算)、回復温度をエンタルピーから。残す閉包の近似を明記する。$k_f$ は較正用の形状と検証用の形状を分けて受理する。非粘性壁を法線にずらしたあと、§4.2 の最終の物理壁の拘束を課した当てはめで確定し、拘束への射影で δ\* 補正をどれだけ戻したかを位置・壁角・曲率の残差として記録して受理値を決める。帰還ループは回さない。
- **設計変数 = $L$ の 1 個**。指定した kernel の下で長さごとに輪郭が 1 本 (S1 で一意性を確かめる) なので MOO はせず、$L_k$ で評価する。

### 4.5 比較の手順

- 両モードを §8-1 の同じ寸法 spec で回す。比較用の**共通の長さ $L_k$ を 5〜7 点**決め、(M) はその $L_k$ で、(D) は MOO (LHS 40 点 + EHVI infill 20〜30 点) の後に各 $L_k$ で $L$ を固定した $(\theta_e,\theta_{max})$ の 2 変数の再最適化 (Phase 2 の固定長最適化と同じ) を行って、同じ $L_k$ の値どうしを比べる (前線の補間で比べない)。
- TOP (Phase 2 の $\theta_n,\theta_e,L$) を 3 点だけ新しい評価器で回し、過去の結果との橋にする。
- **第一の判定量 = 前進 MOC ゲート** (2026-10-11 diagnostician の推奨を採用): §8-1 の寸法 spec・凍結 TP・初期線 (ゲートは CFD を使わないので Hall。C2 スロートの単一曲率なので Hall の前提に合う。CFD ピンとの差は後で記録) で、**壁を与えた前進 MOC** (①の `interior_vec`、analytic・converge) が出口まで衝撃波・同族特性線の交差・網の反転なしに場を構成できるか。実装と判定は④ §4.4 の基部受理ゲートと**同一モジュール** (§5.1 #3b、所有はこの系列、剥離 plan の TC1A の MOC もこれを使う)。交差の検出は `geometry/moc_diagnostics.py` の非隣接 C± 交差の検査を使う。内部衝撃波の指標はこのゲートの連続量として残す。
- **(M) 壁の前進 MOC ゲートの所見 (2026-10-12 diagnostician [記録](../../notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md))**: 前進 MOC ゲートは (M) の Rao 輪郭そのものを不受理にする。これは軸近傍の判定不能 (j ≤ 2) とは別の物理の所見: (M) の輪郭は DE の下 (L100: x ≥ 15.8 r_t = 0.70 L・r ≤ 0.34、L70: x ≥ 9.0〜9.4 = 0.57〜0.60 L・r ≤ 2.9) で壁 BE 直後の圧縮 C⁻ が軸近くで集束して交差する — Rao の構成が制御しない領域の内部衝撃波。推力会計 (制御面・壁) はこの領域を使わない (`dod_clean`)。TOP (x 3.9〜5.0・r 0.9) より下流・軸寄り。④の基部・海面作動の判断材料には位置 x_s/L と r の範囲を併記する。「(M) は内部衝撃波が無い」とは書かない。V2 の Euler run でこの位置に圧縮・衝撃が出るかを記録する。
- **2 モードとも提供し、設計者が選ぶ** (親 §3-16、2026-10-11 ユーザ判断): 問題定義 YAML の `geometry.bell_mode` (仮名: `moc_rao` / `direct_bspline`、旧 `top` は照合用) でモードを選ぶ。比較はどちらかに絞るためでなく、**選ぶための判断材料**を設計者に示すために行う。主な材料は「十分良い形」(§6 V3 の定義) が得られるまでの時間で、両モードと基準 (TOP チャート + CFD 1 評価) で比べる。時間は 1 回の設計キャンペーンでかかる分 (CFD ピンの Euler・$k_f$ の較正などキャンペーンごとの準備 + 評価回数 × 1 評価の時間、同じ GPU の実時間と評価回数の両方) で数え、TP の Rao の導出などの**一度きりの開発費は別に記録**して設計時間に入れない。直接モードは MOO の台帳から「その時点までの最良 $\eta$」の推移を取り、十分良い形に達した評価回数を読む。前進 MOC ゲートは④の基部と海面作動のための判定量として残す。
- **ゲートを MOO の制約に使うときは連続の余裕で** (diagnostician 2026-10-11): 二値の合否は出口面 x = L で崖になる (合体点が x = L をまたぐ) ので、C⁻ 族の間隔比の最小 ($x\le L$) ≥ 0.2 を制約にし、C⁺ 族は指標が直るまで「n 321・161 の両方で反転 0」。C⁺ 族の間隔比の指標は初期線側の前線を基準にするアーチファクト (`moc_forward.py:854-870`、形状によらず 0.84〜0.89) なので、MOO の前に正規化を直す (§5.1 #3e)。
- **η は拒否権に格下げ**: (D) だけが審判 (粘性 CFD) に対して最適化されるので、同じ $L$ で $\eta_D\ge\eta_M$ は探索誤差の範囲で構成的に成り立ち、Phase 2 の 0.133 % から V3 は高い確率で「同等」を返す。$\eta$ は「$|\Delta\eta|-U_k>1.5\times10^{-3}$ で負けた側を落とす」ためだけに使い、方式選択の主判定量にしない。
- 併せて比べる量: $\Delta\eta(L_k)=\eta_D(L_k)-\eta_M(L_k)$、形状が作れない点の率、総費用 (評価回数・時間と、(M) の CFD ピン・$k_f$ 較正の run を含める)、内部衝撃波の指標、(M) の $g$ の根の残差と当てはめ残差。

## 5. 実装ステップ

段階方式 (plan レビューと diagnostician の推奨): **共通評価器 (2a 済み) → (D) の壁生成 (MOO はまだ回さない) → 前進 MOC ゲート (#3b) を (D) の格子点と TOP 最適形に適用 (§6 V0、CFD なし) → 結果で S1 と (D) の MOO の順を決める**。2b (推力会計・ゲート) は CFD 評価の前提なので並行して進める。(M) が S1 を通る前に大規模 MOO を投入しない。

0. `methods/design/overview.md` に「ベルの 2 モード」節を書く (実装前)
1. 共通評価器: 推力定義 (§4.1)・寸法 spec・輸送物性を確定し、TP 化の全経路 (§4.1) とゲート (§4.1) を実装して受理 (§5.1 #2)
2. (D) の壁生成と単体試験、(D) の MOO
3. (M) の Rao 条件の導出と S1 ゲート → 生成器 → 単体試験 (§6 V1)
4. CFD ピンの凍結源と Hall 検算
5. 当てはめの拡張、δ\* の積分法の TP 化と $k_f$ の較正・検証、拘束付きの物理壁の確定
6. 比較の run (case/68、新設) と判定

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1g | ~~#2b のゲートの扱い~~ | **判断: 2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-2b-gates-diagnostician.md) — 残差 PASS 必須を取り下げ、最終 \|cv\| ≤ 1e-4 + 推力系 STEADY (1.5e-4) + 残差は DIVERGED/上昇だけ FAIL、η は出口面、壁法線の近壁層で AR ≤ 5000、pRef は置かない。前提を #2c の継続 run で検証** | F |
| 1v0 | ~~V0-a の解釈~~ | **判断: 2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-v0-gate-diagnostician.md) — V0-a は判定不能 (標本退化)、$x_m$ を内部決定にして回し直す、初期線は Hall n ≤ 321 + ガード、軸の θ_r 代用は合否に効いていない** (主張の $x_m$ の帯は主セッションで再確認) | F |
| 1 | ~~§4・§6 を上位に諮る~~ | **判断: 2026-10-11 codex (diagnose) [記録](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) — NO-GO。C1・B1〜B4 を採用して改訂** (採否は §6.1)。着手順 TP 化 → D → M は採用 | F |
| 1b | ~~寸法 spec (仮値) で plan 段レビュー~~ | **判断: 2026-10-11 codex `--stage plan` [記録](../../notes/reviews/2026-10-11-tooling-nozzle-bell-two-mode-plan.md) — GO-with-changes (C0/M8/m1)、全件採用** (採否は §6.1)。寸法の仮値に桁の破綻なし ($r_e$ 0.990 m、$L_k$ 1.903〜3.172 m、必要 $C_F$ 1.624) | F |
| 2a | ~~共通評価器: TP のガスと 1D 理想値・runner の配線~~ **完了 (2026-10-11、implementer、主セッションで再確認)** | `gas/rocket_exhaust.py` (新規: CEA2 → 排気ファイル `forge_rocket_exhaust_v1`、契約 `RocketExhaust` が同じ Y・種 DB・$T_c$ から MOC の `GasSemiPerfect`・lump `EXH`・`FrozenGas` を作る、1D 理想値 `ideal_1d`)、`probdef.py` の `gas.exhaust`、`evaluate/runner.py` の TP 経路 (`gas.model: semiperfect` + `gas.exhaust` のときだけ: config・`outflow`・IC・壁平滑化・段)。**結果**: 単体試験 `design/tests/run_bell_tp_tests.py` ALL PASS (主セッションで再実行) — $\varepsilon$ 50 で 1D と CEA 凍結の差は出口温度 −5.2e-6・$\dot m_{1D}$ +1.19e-5・出口速度 −3.2e-6・出口圧力 −2.1e-5・$F_{1D,vac}$ +8.7e-6 (合格線 1e-3)、基準ケース $\dot m_{1D}$ 270.87 kg/s・$F_{1D,vac}$ 1159.1 kN・$T_c$ 3523.96 K。CPG 経路は case/40 の 17 本で config がビット同一、design の試験 31/33 PASS (落ちる 2 本 `run_sern_gates_tests`・`run_sern_moc_tests` は変更前の HEAD でも同じく落ちる既知 — sern-chain R13)。起動確認 [`case/68.bell_two_mode/run_0001_tp_soft_smoke`](../../case/68.bell_two_mode/README.md): soft 1000 step、残差・全場に NaN/Inf 0 件 (主セッションで再確認)、`NOT CONVERGED (stalled)` は soft 段なので想定内、メッシュ品質 PASS (AR 836・skew 0.417)。**主セッションの判断**: CEA の生成物を forge の熱物性・輸送データにある H2・O2・H2O・OH・H・O に限る (HO2・H2O2・O3 を除いた影響はモル分率 4.5e-5・$T_c$ −0.17 K、排気ファイルの `omitted_species` に記録) | O |
| 2b | 共通評価器: 推力の CV 会計とゲート **実装済み (2026-10-11、`b630ca99`) — ゲート (i)(ii)(v) が基準ケースを通さず、diagnostician に諮問中** ([ブリーフ](../../notes/reviews/briefs/2026-10-11-bell-2b-gates-diagnostician.md)): CV 残差 −5.2e-7 (`case/68.bell_two_mode/run_0003_tp_cv_continue`、η 0.979220・ṁ/ṁ_1D 0.99255)、推力系は STEADY、配線回帰 |Δη| ≤ 1.5e-5、残差は 2.2〜2.7 桁 (ω 1.6〜1.7) で stalled、cv_residual は相対ドリフトで DRIFTING、スロート y1+ 最大 1.60 (1201×121) | §4.1 の後半。触るもの: `metrics/extract.py` (TP の圧力復元、全境界の $\Pi_x$ による壁会計と出口面会計、CV 残差)、`opt/driver.py` (TP の $\dot m_{1D}$・$F_{1D}$、ゲート (i)〜(v))、`opt/polish.py` (`restart_field.py` による同一メッシュの継続)、runner の `StageManifest`。合格: 全流れ解で CV 残差 ≤ 0.1 %、ゲートの負例試験 (STALLED・DRIFTING を落とす)、**配線の回帰** (同一バイナリ・同一入力場・同一実効設定で新旧の評価コードの $\eta$ 差 ≤ 3e-4、CPG 経路)。過去 run の再現は別扱い。**2a からの引き継ぎ**: (a) TP の同一メッシュでの段の継続 (`restart_field.py` + `StageManifest`。今の TP の `warm_from` は cross-mesh 用)、(b) driver が評価ごとに YAML を書き出すと相対パスの `gas.exhaust` が切れる、(c) スロートの $y_1^+$ が約 5 (第一セル 3.5 µm) — $\le1$ には約 0.7 µm が要り、runner が `wall_first_frac_throat` を読むようにする、(d) `space.pRef` の値 (推奨レシピ §5 は軸対称で動作静圧を求めるが、ノズルは 1e7→1.4e4 Pa と幅があり値が一意でない。既存の軸対称ノズル runner も書いていない — 決める前に上位に諮る)、(e) 起動直後に収縮部 ($x/r_t$ −4.0〜−0.6) の P・T が $P_c$・$T_c$ を超える (step 250 で 1.159e7 Pa・3622 K → step 1000 で 1.068e7 Pa・3567 K と減少中。収束解で残るかを確かめる、解釈はまだしない) | O |
| 3 | (D) の壁生成 **実装済み (2026-10-11、`c636e716`)**: `geometry/wall_bell_c2.py`、拘束の残差 ≤ 2.7e-12・C2・ビット同一・負例 PASS。**V0 の格子 135 点のうち成立は 8 点** (最小ジャーク解が $x_m$ の後で $\theta$ を非単調にする) — 族の見直しを diagnostician に諮問中。 | §4.3。**規約 (2026-10-11 主セッションで固定)**: 座標は $r_t$ で無次元 ($x$ はスロート基準)。拡大部 $[0,L]$ を等間隔 16 区間の clamped 5 次 B-spline (係数 21 個)。拘束 7 個の下で $\int_0^L(r''')^2dx$ を最小化 (KKT、拘束行列の階数を検査)。設計変数の範囲: $\theta_{max}$ 20〜45°、$\theta_e$ 2〜20°、$L$ は $L_k$。$x_m$ は内部決定 (§4.3: 走査 0.5〜4.0 $r_t$ を 141 点 → 成立区間で $J$ 最小を黄金分割、同一入力でビット同一、区間が空なら INFEASIBLE と理由)。成立条件: 2000 点の標本で $r'\ge0$、$\theta'$ の符号変化がちょうど 1 回、$|\kappa|\le1/R_t$。収縮部は既存の U→T 5 次 Hermite ($L_U$ 3.0・$L_{pipe}$ 1.0、$R_t$ 1.0 で単調条件 $\le3.826$ を満たす)。合格: 拘束の残差 ≤ 1e-10、C2、同一入力でビット同一、成立条件の負例試験、収縮部の単調条件 | O |
| 2c | ~~ゲートの組み替えの検証~~ **結果 A (2026-10-11、主セッションで判定)**: `case/68.bell_two_mode/run_0005_tp_cv_continue2` (run_0003 から +24000、本段計 60000、設定変更なし) で $\eta_{exit}$ 0.9791941 (予測漸近値 0.979183 との差 1.1e-5 ≤ 1.5e-4)、継続区間 24 スナップショットの $|cv|$ 最大 4.25e-5 ≤ 1e-4 (cv は 4.5e-6 → 4.25e-5 で床 ≈ 5e-5 に落ち着いた)。`check_quasisteady --series-csv` (drift・osc 1.5e-4) で壁会計・出口面会計・$\eta$・$\eta_{wall}$ が STEADY、ṁ (1e-3) STEADY。残差は `NOT CONVERGED (stalled)` (併記)。NaN/Inf 0。→ §4.1 の組み替えを driver に反映する (#2c' ) | O (実行) / F (解釈) |
| 2c' | ~~ゲートの組み替えを driver に反映~~ **完了 (2026-10-11)**: run_0005・0007・0003 が (i)〜(iv) 通過、run_0004 は η TRANSIENT、run_0002 は cv −6.6e-4 で不合格 (想定どおり)。末尾窓 8 点以上を (ii) の条件に、ṁ の漸近値は 1e-3 で判定 | `opt/driver.py` の (i) を「DIVERGED と有意な末尾上昇だけ FAIL、NOT CONVERGED は併記」、(ii) を推力系 4 列 1.5e-4・ṁ 1e-3 (cv_residual は列から外す、単調列の漸近値も 1.5e-4)、(iv) を最終 $|cv|\le1\times10^{-4}$、$\eta$ は出口面を主。本段の出力間隔を本段長の 1/20 以下に。合格: run_0003・0004・0005 の判定が §4.1 の想定どおり (0005 は PASS、0004 は η が末尾で動いて TRANSIENT、0002 は cv で落ちる)、負例試験の更新 | O | run_0003 を同一メッシュ・同一設定で +24000 継続 (`case/68.bell_two_mode/run_0005_tp_cv_continue2`、本段計 60000、出力 1000 step ごと、**変える設定は無し**)。**A**: $|\eta_{exit}(48000)-0.979183|\le1.5\times10^{-4}$ かつ継続区間の全スナップショットで $|cv|\le1\times10^{-4}$ → §4.1 の組み替えを driver の定数に反映。**B**: いずれか外れる → 壁会計の内訳 (摩擦か入口面か) を分けてから次の 1 変更 (cfl 4→6 か relax 0.7→1.0) を諮る | O (実行) / F (解釈) |
| 2d | ~~壁法線の近壁層~~ **完了 (2026-10-11)**: `case/68.bell_two_mode/run_0008_tp_wn_staged` (TOP、1201×121) で y1+ 最大 0.82・AR 最大 2700・AR > 1000 のスキュー最大 0.006、全ゲート PASS。run_0007 との差 η_exit −6.8e-5・摩擦 −8.0e-5 F_1D (y1+ 5.7 → 0.82 の格子感度)。`mesh2d` に `wall_normal_end_stations` (出口・入口面を平面に戻す) を追加、`check_wall_resolution.py` は float32 の壁ダンプ × float64 の格子 (merge で変換器が倍精度化) に対応 (`test_gate_bad_input`・`test_wall_resolution` PASS) | #2e の merge で入る `mesh2d` の `wall_normal_layer` ($d_n$・$d_b$、2026-10-08 の冷却壁の格子用) を使う (自前で書かない)。runner の格子キーの受け渡しに `wall_normal_layer` を足す。合格: スロート y1+ ≤ 1 (`check_wall_resolution.py`)、AR ≤ 5000 かつ AR > 1000 の全セルでスキュー ≤ 0.1、入口・出口面は平面、旧格子との $\eta_{exit}$ と摩擦の差を記録 | O |
| 2e | **近軸の hoop 閉性の修正を取り込む** — **ユーザ決定 (2026-10-11)「修正は必ず取り込む」** | 修正 0b4dff4e (`mesh.axisSegmentRWeight`、既定 1。変換器が `/PLANES/rSurfVect` を書く) は `feature/nozzle-wall-fit-and-pipeline` と `feature/nozzle-core-grid` にだけあり、main には未統合。**取り込み方**: push 済みの `feature/nozzle-wall-fit-and-pipeline` を作業ツリー `forge-thruster-merge` のブランチ `integrate/thruster-m6-hoopfix` で merge (衝突なし、253 commit。M6 の作業ツリーには触らない、調整板の依存関係に記載)。同じ merge で `mesh2d` の `wall_normal_layer` (壁法線の近壁層、#2d に使う) と倍精度の格子座標も入る。**回帰** (取り込み後、記録): design の全試験、風洞の MOC 試験、ベル TP の基準ケースを新しい変換器で変換し直して同じ段階起動・継続 (run_0003/0005 と同じ長さ) を回し、$\eta_{exit}$・cv・ṁ の差を記録 (合否でなく効きの記録)。その後 `feature/thruster-nozzle-design` に merge。**回帰の結果 (2026-10-11、作業ツリー `forge-thruster-merge`)**: ビルド成功、変換器が `/PLANES/rSurfVect` を書き forge が `axisSegmentRWeight: ON` を出す。基準ケースを新しい変換器で作り直して同じ長さ (`case/68.bell_two_mode/run_0006_tp_hoopfix_staged` → `run_0007_tp_hoopfix_cont48k`、本段 12000 + 48000) — 修正前の run_0005 との差: $\eta_{exit}$ +5.1e-6、$\eta_{wall}$ −3.4e-5、ṁ +5.0e-6、摩擦 +6e-5、**CV 残差 4.25e-5 → 3.7e-6** (壁会計が −39 N 動き、2 つの会計の食い違いが 1 桁縮んだ。後半 24 点の |cv| 最大 3.5e-5)。推力系 4 列 STEADY (1.5e-4)、残差は NOT CONVERGED (stalled、併記)、NaN/Inf 0 | O |
| 2f | CPG の CV 残差 3〜5e-3 (記録) | `outlet_statPress` の亜音速分岐と節点 P の不一致で説明できる見込み。確認は run 不要: 出口列で $M_n<1$ の節点 P を Ps に置き換えて再積分し、≤ 1e-3 に落ちれば確定。TP の受理には関係しない | O |
| 3f | ~~runner への (D) の組み込み~~ **完了 (2026-10-11)**: `case/68.bell_two_mode/run_0009_tp_d_l80_tm25_te10_staged` (L80・θ_max 25°・θ_e 10°) が全ゲート PASS、η_exit 0.971897、前進 MOC ゲート pass (C⁻ 比 0.961)。**観察**: ω 残差のプラトーが run_0008 (TOP) の約 200 倍、rms_ro・rms_roe も約 5 倍高い (推力系は STEADY)。TOP 基準ケースでの「プラトーは推力に 1e-4 以下」という性格づけが (D) にも当てはまるかは未確認 — (D) の MOO で推力の STEADY が崩れる点が出たら調べる。以下は元の記述: **runner への (D) の組み込み** (2026-10-11、V0 を閉じた後の本題 — ユーザ了承「やっていいよ」) | `evaluate/runner.py`・`geometry/wall.py` 側に `geometry.bell_mode: direct_bspline` を足し、`design_bell_c2` の壁 (C2 スロート・U→T Hermite の収縮部) で CFD 評価できるようにする (TOP・hermite は従来どおり)。INFEASIBLE は理由つきで driver に返す。前進 MOC ゲートの指標 (C⁻ 間隔比の最小・反転・交差) を評価ごとに台帳へ記録する。**③の真空の比較ではゲートを制約にせず記録だけ** (ゲートの用途は④の基部と海面作動。制約に使うのは④の MOO から、#3e の後)。合格: TP 基準ケースで (D) の 1 点が段階起動から §4.1 の組み替えたゲート (#2c') を通る | O |
| 3d | ~~(D) の成立条件の穴と区間の端~~ **完了 (2026-10-11)**: 成立条件 (ii) を x_m に結んだ符号条件に、区間の端を二分法 (1e-6 r_t) で詰めた。接触型 12 点が INFEASIBLE、残り 29 点 OK、x_m は連続 (θ_e 0.01° 刻みで \|Δx_m\| ≤ 7.95e-4)、試験 4 本 PASS。残る注意: J 最小で 29 点中 18 点が成立集合の境界 (詰めた端) で確定し、標本 2000 点の間で tol を超える場合がある (L90・25°・10° で +7.6e-7) — MOO で問題になったら端から余裕を取る。以下は元の記述 (diagnostician 2026-10-11) | `wall_bell_c2.py`: 成立条件 (ii) を「$r''\ge-tol$ on $[0,x_m]$ かつ $r''\le tol$ on $[x_m,L]$」($tol$ = `RESID_TOL`) に置き換える (x_m で r'' が 0 に接して符号を変えない壁が通っていた — `_sign_changes` が微小値を飛ばすため。V0-b の 2 点で最大壁角が θ_max を超え、細い成立区間 10 点も同型)。成立区間の端を走査点と隣の不成立点の間で二分法で詰めてから黄金分割 (x_m の 0.025 刻みの量子化を消す)。J 最小の規則は変えない。合格: 接触型の負例 (旧 V0-b の 2 点と sliver 10 点) が INFEASIBLE、既存 29 点が OK のまま、θ_e を 0.01° 振って x_m が連続 (0.025 の跳びなし)。続けて V0-c (§6) | O |
| 3e | C⁺ 族の間隔比の指標の正規化 | `moc_forward.py:854-870` の基準前線から初期線近傍を除く。合格: 形状によらない最小値 (n 321 で 0.84 @ x 0.93、n 161 で 0.885 @ x 0.20) が消え、曲げ戻しの負例で C⁺ 族の詰まりを映す。(D) の MOO でゲートを制約に使う前に | O |
| 3c | **V0 の回し直しの準備 (diagnostician 2026-10-11)** | (a) #3 の $x_m$ の内部決定、(b) #3b のガード 3 点: Hall 初期線 ($R_t\le1$) で n > 321 を拒否 (321 可・641 不可を docstring に)、前線 1 に反転がある場合は壁の不受理でなく「判定不能 (初期線不整合)」、不受理の根拠の点に軸からの添字 j を記録し j ≤ 2 の点だけなら「判定不能 (軸近傍)」。(c) V0-a の出力を `v0a_*` に退避してから回し直す。合格: 既存の試験 + 新しい試験 (内部決定の決定性・空区間の INFEASIBLE、ガード 3 点の負例) PASS | O |
| 3b | **前進 MOC ゲート (共用モジュール)** **実装済み (2026-10-11、`c636e716`)**: `geometry/moc_forward.py`、放射源流で壁 2 次・流量 2 次、負例 PASS、風洞の試験不変。未決: Hall 初期線 (R_t 1) の適合式との不整合 1.1e-2 rad で n ≥ 641 が破綻、前進の軸の θ_r は 1 次の代用 — 諮問中 ([ブリーフ](../../notes/reviews/briefs/2026-10-11-bell-v0-gate-diagnostician.md))。 | 壁を与えた前進 MOC (①の `interior_vec`、analytic・converge、対の分類を不受理の理由に、`moc_diagnostics.py` の交差検査、亜音速入力の検査) を新規モジュールに書く。③のゲート・④の基部受理・剥離 plan の TC1A の全流れ壁圧が共用し、**所有はこの系列** (別セッションには使い方だけ渡す)。合格: 風洞の放射源流の試験 (`design/tests/moc_axis_limit_radial.py`) を前進でも再現、衝撃波を要する壁と亜音速入力の負例で不受理の理由が正しく返る、`moc_kernel.py`・`moc_inverse.py` の既定はビット不変 | O |
| 4 | ~~(M) TP の Rao 条件の導出と S1 ゲート~~ **判断: 2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md) — 条件は γ-free で古典と同形、固定 (L, ε) で縁条件は落ち、κ_B・θ_B・s_D の入れ子の 1 次元根で閉じる。成立の見込み高。TIC への切り替えの目安: 5 本の L_k のうち 2 本以上で根が括れない / T2 が n で収束しない / T4 の c₁ が負対照の 5 % に落ちない ((ii)(iii) は 1 回修正してから)** | §4.4。凍結 TP の $h(T)$・$s(T,p)$ で拘束付き変分 (固定 $\varepsilon$ の端点条件を含む)。kernel の族・$\kappa_B$ の探索域・B と D の枝の追跡・複数根の選択規則を定義。平面縮約が `rao_planar.py` と一致するかも確かめる。S1 を通るまで (M) を比較に入れない | F |
| 5 | (M) 生成器と単体試験 **WIP (`56dc8013`) — 実装者が規則の穴で停止 → diagnostician 判断 (2026-10-11 [記録](../../notes/reviews/2026-10-11-bell-s1-generator-holes-diagnostician.md)): |c| 規則を撤回、DE を r パラメータに。T1b・T2・T1d・accessor は PASS (円弧 kernel・EXH・n 161 で、旧規則のまま L_k 100 % だけ成立: θ_B 29.76°・θ_E 5.26°、推力 2 会計の差 −1.39e-4 @ n 161 / −4.4e-5 @ n 321)** **→ (a)〜(d) 実装済み (`8f97ec76`)、事前登録の判定は A (2026-10-11)**: 円弧 kernel・EXH・R_t 1・ε 50・n 161 で L_k 70・80・90・100 % に主枝の根 (θ_B 36.09/33.03/31.12/29.76°、θ_E 10.91/8.63/6.76/5.26°、θ_max 36.10/33.35/31.77/30.71°、推力 2 会計の差 −1.20〜−1.26e-4)、60 % は根なし (主枝は θ_B 37.6° で det_D → 0⁺ に終わり x_E の最小 15.18 > 13.59 — 「(M) で作れない L_k」1 本、TIC の閾値 2 本に届かない)。跨ぐ D の T2 は 3.5e-12 以下、試験 25 件 PASS (主セッションで再実行)。結果 `case/68.bell_two_mode/rao_m/rao_lk_n161.csv`・`rao_l60_scan.csv`。**T3〜T6 の前に上位に諮問中** ([ブリーフ](../../notes/reviews/briefs/2026-10-11-bell-rao-mainbranch-c2.md)): θ_E ≤ θ_B は θ_B ≳ 40.8° の偽枝 (det_D < 0) を通す、C2 の g̃ が c にほぼ依らない (L 100 % で c·R_t² 0/−0.35 に −0.0690/−0.0677)、T3 は n 161 で 1.2e-4、1 設計 約 30 分 (C2 は約 5 時間/本)、`solve_kappaB` の外側の枝は添字のまま **→ 判断 (2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-rao-mainbranch-c2-diagnostician.md)、ユーザ了承)**: C2 は撤回 (linear 族に根なし、部品試験で判定 A を再現)、(M) = 円弧 kernel・B で C1、受理条件は det_D > 0。**次**: (e) 受理条件を det_D > 0 に (θ_E ≤ θ_B・de < 0 は記録)、`design_rao` の既定を円弧に・`solve_kappaB` は外側の根を取らない (c は明示パラメータ)、(f) T6: 円弧で L100・L70 × n 241/321 (161 は済み) — T3 を同じ run で判定、Δκ_B と g̃ を記録、(g) T4・T1a′ は T6 の後 **→ (e)(f) 完了 (2026-10-12)**: 受理条件 det_D > 0 (θ_E ≤ θ_B・de < 0 は記録、det_D > 0 で θ_E > θ_B は `policy_hole` に記録)、`design_rao` の既定は円弧、`solve_kappaB` は c を明示で受けるだけ。T6 = 円弧・L100・L70 × n 161/241/321 + n 321 の dr 半分・Goursat の m・q 倍 (`case/68.bell_two_mode/rao_m/rao_t6.csv`)。全 10 本で候補 1・det_D > 0・policy_hole 0、n 161 は `rao_lk_n161.json` とビット一致。**T3 (事前登録の基準) は満たす**: 推力 2 会計の差 L100 1.25e-4 → 6.15e-5 → 4.07e-5 (次数 1.75/1.43)、L70 1.20e-4 → 5.72e-5 → 2.87e-5 (1.82/2.40)、閉合 L100 2.70e-5・L70 9.33e-5 (n 321、次数 2.0)。**T6**: θ_B・θ_E・θ_max・F_cs は 3 水準で単調 (θ の観測次数 1.8〜2.5、外挿 L100 θ_B 29.7567°・θ_E 5.2555°、L70 36.0929°・10.9085°)、dr 半分は全量で 1e-14 以下 (DE の刻みは律速しない)、Goursat の m・q 倍で 2 会計の差 −4.07e-5 → −1.15e-5 (L100) — θ_max・F_wall・Δκ_B は m に律速される (κ_W(B⁺) を B から 6 点で当てはめるため)。**T5 の記録**: Δκ_B L100 −0.92712/−0.92519/−0.92428、L70 −0.99659/−0.99624/−0.99604、g̃ L100 −0.06901 (n で 3e-5 以内)・L70 −0.16876/−0.16867/−0.16865。**次 (g)**: 生産の Goursat を m = q = n にする (費用 2 s、2 会計の差が約 1/3)、T4 停留性・T1a′・T1c **→ (g) 完了 (2026-10-12、`431c7e2f`)**: Goursat の既定 m = q = n、前進 MOC で壁を解き直す推力の評価器 `forward_thrust` (放射源流で 2 次)。**T1a′ は満たす** (前進 MOC ゲートを通った (D) の最良が (M) より 2.33e-3 [L100]・1.86e-3 [L70] F_1D 低い、`rao_t1ap.csv`)。**T1c は 1e-13 で一致** (言い換え後の基準を満たす)。**T4 は判定不能** (上記、diagnostician 2026-10-12 [記録](../../notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md)) → **(h) A/B (事前登録)**: 壁 = `rao_t6_parts/L100_n321_mq2_wall.npz` の BE、初期線 = n 321 の Hall マーチの整合な前線 k₀ を r 等間隔に N ∈ {321, 641, 1281} へ 3 次補間 (新規関数、前線の status が軸まで OK であることと `_check_line` を先に確かめる)、族 sin_c1 k1、ε ±{1, 2, 4}e-3 の 7 点、スプラインと台形の両会計で c₁ε₀・c₂ε₀²・ε*・F0 の N 収束。**A**: N 1281 で \|c₁ε₀\| ≤ 1e-7 かつ両会計の c₁ε₀ の差 ≤ 1e-7 かつ 641 → 1281 で \|c₁ε₀\| が減る → k2・k3・L70・ctrl_plus を N 1281 で回して T4 PASS、S1 を閉じる。**B**: 両会計の差 ≤ 1e-7 だが \|c₁ε₀\| が 641 と 1281 で 30 % 以内に一致し両方 > 1e-7 → T4 FAIL (壁は ε* だけ非停留)、ε* と懸かる推力を記録し S1 の扱いはユーザ判断。**C**: N 1281 でも 2 会計の差 > 1e-7 → 評価器はこの線を解像できない、T4 は恒久に判定不能、B と同じ扱い。S1 は (h) の結果まで閉じない **→ (h) の判定 A (2026-10-12、`e6e01cb5`、diagnostician 追加諮問 2)**: 再標本化前線 N 321/641/1281 で c₁ε₀ スプライン 2.92/2.57/2.62e-8・台形 7.56/3.68/2.96e-8、2 会計の差 3.4e-9、線 9.9e-8 に対し比 0.013〜0.015、c₂ε₀² −1.979e-6、ε* 2.6〜3.0e-5 r_t、懸かる推力 ≤ 1.3e-11 F_1D。A3 (641 → 1281 で減る) は解像後に成り立たない条件だった (登録の欠陥、判定は 2 会計の一致幅で)。**次 (h′)**: k2・k3・L70 (`L70_n321_mq2`) と負対照 (既存 `rao_t4_parts/gen_*_n161` の `ctrl_plus` 壁、生成器は回し直さない) を N 1281 で、任意で k₀ を 1 本下流に変えた k1。全て線以内 → S1 を閉じる (札は記録のとおり、「最大推力」とは書かない)。1 つでも外れる → 2 会計差が線以下なら再諮問、線以上なら N 2561 を 1 本。**欠陥票**: Hall のマーチの最初の前線の軸点 θ_r の符号振動による 2Δ の鋸歯の帯 (r 0.62〜0.74、θ 1.45e-3 rad) が kernel 場・BD・Goursat・壁流線 (x 3.4〜5.3) に乗る — #3b・#6 で、CFD ピン初期線で再発するかと (M) の壁の x 3.4〜5.3 の差を記録。T1a′ は前任の評価器の床の 4 桁上なので再計算不要 **→ (h′) 完了・S1 を閉じた (2026-10-12)**: N 1281 の再標本化前線で sin_c1 の全組が線の内側 — 比 \|c₁ε₀\|/\|c₂ε₀²\| は L100 k1/k2/k3 で 0.013/0.010/0.005、L70 で 0.006/0.005/0.002 (スプライン、台形も ≤ 0.015)、c₂ < 0 は全行、2 会計の差 1.9e-9〜1.3e-8 (L70 k3 だけ 1.33e-8、線の 1/43)。負対照 (既存 n161 の ctrl_plus 壁) の c₁ε₀ は 4.35e-6 (L100)・6.25e-6 (L70) で床の 1000 倍以上、(M) との比は 0.60〜0.68 % (L100)・0.21〜0.26 % (L70) で線 5 % の内側 (`case/68.bell_two_mode/rao_m/rao_t4h_follow.csv`)。任意の k₀ 感度は回していない — 壁点 x ≈ 0.45 の前線は θ > μ (最大 +0.05°) で r が単調でなく初期線にならない (鋸歯の帯と重なる、欠陥票に添える)。**S1 の札**: 「(M) は指定した kernel (円弧、Hall n 321 初期線 — 鋸歯の帯あり) と拘束 (L, ε, ṁ) の下の停留輪郭。sin_c1 k1〜3 で \|c₁ε₀\|/\|c₂ε₀²\| ≤ 0.05 (最大 0.015)、ε* ≤ 3e-5 r_t、懸かる推力 ≤ 2e-11 F_1D、c₂ < 0、評価器の床 3e-9」。「最大推力」とは書かない。**次**: #6 (CFD ピン・T1a の記録・鋸歯の帯が CFD ピンで再発するか)、V2 (Euler で MOC 推力との整合 0.2 %)、#7 (当てはめ・δ* の TP 化・物理壁) | §4.4 の構成と §6 V1 の T1–T6。新規 `geometry/rao_axisym.py` (`kernel_wall`・`kernel_field`・`trace_bd`・`rao_line`・`solve_D`・`solve_thetaB`・`solve_kappaB`・`goursat_fill`・`wall_from_massflow`・`thrust_accounts`・`design_rao`)、等エントロピーの accessor (`isentrope(M)` → T, V, ρ, p, a, μ, k — `gas/semiperfect.py` は他ブランチも触るので新規モジュールの関数として足すのを優先)、試験 `design/tests/run_rao_axisym_tests.py`。Rao チャートの出典 (γ・R_t・長さの定義) は case/29 の資産で確かめる。**次 (2026-10-11 判断)**: (a) `rao_line` の DE を r (弧長) パラメータに、解なしは θ ≤ 0・θ+μ ≥ 90°・M が表の外・流量に達しない、(b) `solve_thetaB` の枝を η_D の連続性で追い、受理条件 $\theta_E\le\theta_B$、(c) T1b・T2 に折り返しを跨ぐ D の組 (`synth_bd` M_B 2.2・θ_B 32°、折り返し点で dθ/dr = −tan θ/(2r))、(d) 円弧 kernel・EXH・R_t 1・ε 50・n 161 で L_k 60〜100 % の 5 本の表 (根の有無・θ_B・θ_E・θ_max・κ_W(B⁺)−κ_B・推力 2 会計の差、60 % は θ_B 36〜46° を細かく)。**A**: 4 本以上で主枝の根 + 跨ぐ D の T2 が 1e-5 → T3〜T5 (C2)・T6 (n 161/241/321 + n 321 の dM/m/q 細分)・T1a′。**B**: 2 本以上で根なし → §8-4 の TIC (諮問なし)。T1a は記録に格下げして #6 へ。C2 の根が κ_B < 0 に落ちたら明記 | O |
| 6 | CFD ピンと Hall 検算 **— W1 完了 (`7e89630d`)、W2 完了 (2026-10-12、V2 の Hall 版は FAIL・原因は初期線の流量 — §6 V2)** **— W3 (L100) 完了 (2026-10-12、`6c755591`)**: (a) 照合 run_0015/0016/0017 とも 15 項目通過 (照合できないのは gamma_hall だけ — ピン線は Hall の M″ を使わない)。`cfd_initial_line` の③の照合は以前は gas を何も比べずに「照合済み」にしていた — 排気の記録・熱力学の参照値で照合するよう直した (①はビット同一)。(b) ṁ ドリフト 3.0e-7。(c) 格子 ×1.5 (run_0017) とのピン線の差 m*/ṁ_1D −7.2e-7・θ 4.7e-6 rad → 線 1e-4 の内側、run_0015 を凍結源に採用。**ピン線** m*/ṁ_1D 0.9926224 (Euler 0.9926297 との差 −7.4e-6、Hall 0.9897056 より +2.9e-3)、C⁻ 適合式の不整合 +9.4e-5 rad (Hall −1.09e-2)、**鋸歯の帯は再発しない** (帯の θ のずれ 2.5e-5 vs Hall 1.45e-3)。**ピン版 L100 vs Hall 版** (m = q = 320): θ_B 29.8009° (+0.0445°)、θ_E 5.2506° (−0.0048°)、θ_max 30.7264° (+0.012°)、壁 Δr 最大 9.46e-4 r_t (目安 1e-3 の 95 %、余裕は小さい) → 目安の内側で V3 の U_k に入れない。η_cs 0.990345 (Hall 版 0.987530)、推力 2 会計の差 −1.30e-5、閉合 6.8e-6、Δκ_B −0.9256 (`case/68.bell_two_mode/rao_m/rao_pin_compare.csv`)。L70 のピン版は未作成。ピン版 V2 = `run_0018_tp_m_euler_l100_pin` (B、実行中) | 合格: 凍結源の照合 PASS、Hall との初期線の差を記録 (合否なし)。**T1a の記録をここで行う** (diagnostician 2026-10-11): case/29 の R_t 0.382 の TOP 壁の Euler 1 本を凍結源に、CPG γ {1.2, 1.3, 1.4} で (M) を回してチャート値 (θn 38.5°・θe 7.5° @ ε 50) がその帯のどこにあるかを記録。n 641 の T6 もここで記録として足す。**段取り (2026-10-12、主セッション)**: (W1) runner に (M) の壁 (`geometry.bell_mode: moc_rao`、生成器の壁ファイル [`rao_t6_parts/*_wall.npz` の形式] を読む — 生成器は回さない) と TP の Euler 設定 (§4.1「Euler の設定」: `viscMethod: 0`・`visc 0`・`thermCond 0`・`turbulence none`・物理壁 `slip`) と Euler 用の格子 (壁際を細分しない — case/45 の教訓、NS 用の壁法線の近壁層は使わない) を足す。(W2) B で Euler 2 本: Hall 初期線の (M) 壁 L100・L70 (`*_n321_mq2_wall.npz`)。用途は V2 (Euler の推力と MOC の F_cs が相対 0.2 % 以内) と CFD ピンの凍結源 (Hall 版の壁の場 — ピン壁自身の Euler は循環定義になるので使わない)。(W3) `CFDPinnedThroat` の線を `rao_axisym` の初期線に通し (`design_rao` に初期線を渡せるように)、CFD ピン版の (M) を L100・L70 で生成 → Hall 版との差 (θ_B・θ_E・θ_max・F_cs・壁の x 3.4〜5.3 の差) と鋸歯の帯が再発するかを記録 (合否なし)。(W4) T1a の記録 (case/29 の R_t 0.382、CPG γ 3 水準) と n 641 の T6 は W3 の後。**W3 の前段 (diagnostician 2026-10-12)**: 凍結源 = run_0015 (L100)、凍結の前に (a) `check_source_matches` PASS、(b) ṁ の末尾 10 スナップショット (15000〜30000) の相対ドリフト ≤ 1e-5、(c) 格子確認 = run_0015 を ni 1801・nj 181 (×1.5) で 1 本、ピン線の m*/ṁ_1D の差 ≤ 1e-4 かつ線上の θ の最大差 ≤ 1e-4 rad なら現格子を凍結源に (外れたら細格子を凍結源に)。run_0016 から取った線との差を照合の記録に。**W3 の記録項目**: ピン線 m*/ṁ_1D (期待 0.9926 ± 2e-4) と Hall 0.989706 の差、Hall 版とピン版の (M) の差 (θ_B・θ_E・θ_max・壁 x 3.4〜5.3・F_cs — Δθ_B・Δθ_E ≤ 0.1°・壁 Δr ≤ 1e-3 r_t なら V3 の U_k に入れない)、鋸歯の帯が再発するか | O |
| 7 | 当てはめの拡張・δ\* の TP 化・物理壁 **— 2026-10-12 ユーザ決定で向きを確定**: 補正は設計者の選択肢、方法は理論式 + NS 較正 (§4.4)。(7a) η_M0 = ピン版 L100 の非粘性壁を (D) の recipe で NS (補正なしの評価値 + 較正の元データ)。(7b) 較正の設計を上位に諮る (理論式・較正の形・合格の線・ε_inv とスロート)。(7c) 実装と検証形状での確認 **— (7a) 完了 (2026-10-12)**: (M) の NS 4 本 (補正なし) `run_0019`〜`0022` は全ゲート PASS、η_exit L70 0.970936・L80 0.977196・L90 0.980506・L100 0.982155。**(7b) 較正の設計 = codex (diagnose) 2026-10-12 [記録](../../notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md) を全件採用**: (M1) 抽出は一定 x の断面で (`remap_constant_x=True` — NS 格子の列は x が最大 0.282 r_t ずれ、列のまま抽出すると途中の δ_r が最大 31 % 小さく出る、主セッションで再現: x 1.0 で 0.00317 → 0.00416)。(M2) 帯の感度は実際の帯位置 y_b を動かして見る (同じ壁では delta_in = 0 で k_band の感度系列が主系列とビット一致していた)。(M3) δ_r の定常性を時系列で判定 (推力の STEADY は δ* に及ばない)。(M4) 補正しない損を 2 次で小さいと決めつけない — 補正の効きを符号つきで測る。(M5) 抽出値は一定 x の半径方向等価厚さ δ_r — 補正は固定 x の半径補正に統一。(m1) 理論式は既存の `contur_v2` (μ(T)・エンタルピーの回復温度・TP の d ln u_e/dx 込み) を③に接続して検証 (新たな TP 化は不要)。**較正の形**: δ̂_r = kδ·δ_r,CONTUR の定数 1 個を L100 で較正、L70 で独立検証 (ピン版 L70 の壁の Euler 基準が要る — Hall 版の Euler で代用しない)。**事前登録の受理線**: 抽出不確かさ ≤ 5 %、L70 の出口の誤差 ≤ 10 %・補正区間の相対 RMS ≤ 10 %・局所 ≤ 20 % (分母に抽出不確かさの床)、外れたら定数 1 個の移植性を棄却 (L70 を追加較正して検証成功にしない)。**補正の幾何**: ε_inv = (√ε − δ_r,e/r_t)² (今の配列で約 48.12、生成器を回し直す約 25 分は設計費)、物理スロートは固定。**補正の効きの測定**: L100 で補正壁を同じ NS recipe で評価、主判定は補正後の非粘性目標に対する壁圧・境界層外縁 Mach・出口分布の誤差 (壁圧誤差を補正前から 30 % 以上減らし、他の指標を不確かさ以上に悪化させない)、η は改善必須にせず Δη と不確かさを報告。**次**: (7b-1) 抽出の修正と受理 (③の run 情報の `x_E` 欠け・一定 x・帯位置の感度・設定と感度の保存・δ_r の時系列の準定常判定)、(7b-2) ピン版 L70 の壁の Euler 1 本 (B)、(7b-3) `contur_v2` の③への接続・未較正の分布・L100 で kδ・L70 で検証、(7b-4) 補正壁 L100 の生成と NS | §4.2・§4.4。合格: 拡張した当てはめが旧 `joint_fit_wall` の拘束を与えたときビット一致、積分法の速度勾配が TP の数値微分と相対 1e-4 以内、$k_f$ の検証形状 (較正に使わない形状) での出口 δ\* の差を記録、最終壁が §4.2 の拘束を 1e-6 以内で満たす、射影の残差 (位置・壁角・曲率) が受理値以内 (値は #4 の後に事前登録) | O |
| 8a | **(D) の MOO (2026-10-11 着手、AWS B、ユーザ了承)** — doe_000 PASS (η_exit 0.96895、1192 s)。doe_001〜018 は warm start の `interp_field` が化学種の署名を解決できず 10 秒ずつ SUSPECT で落ちた (共用 `forge_species.py` の分子量表に lump の構成種 H・O・OH が無かった) → `6081993` で共通データ `data/species/forge_species_v1.yaml` から引くよう修正、driver は stderr を台帳に残し SUSPECT が 3 回続いたら止まる。失敗 18 行は `ledger_infra_fail_20261011.jsonl` に退避して台帳から除き、同じキャンペーンで再開 (B の forge-8a の単独 run の終了を待って) **→ 停止 (2026-10-11 21:55 JST)**: 再開後の warm start 2 点がどちらも問題 — doe_001 (θmax 31.96°・θe 12.39°・L 16.20) は本段で発散 (出口の壁の角から逆流、step 41558 で NaN)、doe_002 (37.10°・2.63°・20.09、η_exit 0.96572) は壁解像で FAIL (y1+ 最大 1.096、1 超え 7.3 %)。冷間の A/B `run_0011` が PASS → §4.1「起動」・#8a′ へ。run_0010 の台帳は記録として残し、新しい点を混ぜない | `case/68.bell_two_mode/problem_bell_tp_d_moo.yaml`、キャンペーン `case/68.bell_two_mode/run_0010_moo_d` (B の `~/forge-thruster`、ブランチ `feature/thruster-nozzle-design`、RelWithDebInfo・sm_86 のビルド)。LHS 40 点 (候補 120 点のうち成立 57 点から先頭 40) + EHVI infill 10 回 × 2 点、seed 0、基準点 (−0.90, 23.5)、目的 = −η_exit と L/r_t。1 評価 = 段階起動 (soft 2000 → mid 2000 → 本段 60000)、ゲートは §4.1 の組み替え後。前進 MOC ゲートは記録だけ。同時は 1 本 (共有 GPU)。評価ごとの「その時点の最良 η」の推移を設計時間の判断材料 (§6 V3) に使う | O (実行) / F (解釈) |
| 8a′ | **(D) の MOO のやり直し (事前登録、diagnostician 2026-10-11 [記録](../../notes/reviews/2026-10-11-bell-d-moo-warmstart-diagnostician.md))** **— (2) は判定 B → ×0.5 に確定 (2026-10-11)**: `run_0012_tp_d_corner_wf175` (×0.7) は収束系のゲートを全部通ったが y1+ 最大 1.207・1 超え 6.2 % (x 1.477 m の再圧縮域) で壁解像 FAIL。`run_0013_tp_d_corner_wf125` (×0.5) は全ゲート PASS (y1+ 最大 0.878・超過 0 %、η_exit 0.927032、×0.7 との差 −4.2e-5、AR 最大 3366)。範囲の角で冷間起動が通ることも確認。問題 YAML の `wall_first_frac` を 1.25e-5 にして本番 `run_0014_moo_d` を投入。**(4) 格子感度の記録** (`run_0014_moo_d`、1.25e-5、冷間): doe_000 η_exit 0.968926 (run_0010/doe_000 [2.5e-5] 0.968948 との差 −2.2e-5)、doe_001 0.970223 (run_0011 [2.5e-5] 0.970260 との差 −3.7e-5) — 期待どおり \|Δη\| ≲ 1e-4。どちらも全ゲート PASS・warm_from 空、1 評価 約 190 s (GPU 単独)。**(6) DOE の RETRYABLE (2026-10-12)**: doe_024 のみ (θ_max 27.17°・θ_e 2.02°・L 18.57、η_exit 0.97214 [「DOE 最良」は誤り — 下の訂正]、ゲート (i) で rms_roOmega が RISING、他ゲート PASS・推力系 STEADY)。判断 (diagnostician、記録 `notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md` の追加諮問): `check_convergence.py` の spike 条件 (末尾 2 % の最大 > 窓中央値の 10 倍) がプラトーのジッタ (p90/中央値 ≈ 15〜30) で常時真になり、rising の判定が窓平均の比に退化した誤検出。**確認 (主セッション、B 上)**: spike True・slow False・末尾 2 窓の平均比 1.189 に対し中央値比 0.739 (下降)、末尾最大 3.07e5 @ step 59898 (窓中央値 543) → 誤検出で確定。扱い (a): 設定を変えず同一メッシュで本段 +24000 (`runner --continue-from`) = `run_0014_moo_d/retry_doe_024`。**事前登録の判定**: **A** `--segment` (本段 84000) が RISING でなく他ゲート PASS・推力系 STEADY → `retry_doe_024` を PASS として台帳と Pareto に入れる。**B** 再び RISING だが 3000 step 窓の中央値が末尾 2 窓で増えていない (比 ≤ 1.05) かつ slow False → ツールの既知の弱点として記録し外したまま。**C** slow True または中央値が 2 窓続けて増える → 有意、外したまま ω の上昇の場所を見てから諮る。判定ツールはキャンペーン中は変えない。同型が 2 件目で欠陥票 (plan `tooling-convergence-and-wall-resolution-gates`、候補: spike を p99 か全窓最大の比で、窓比を中央値か log 平均で)、直したらキャンペーンの全点を同じツールで再判定。キャンペーン中に見るものに「RISING だけで落ちた件数と 3 行 (中央値比・spike/slow・末尾最大の step)」を足す。**→ 判定 A (2026-10-12)**: `retry_doe_024` の `--segment` (S3_main + S4_main_p、84000 step) は STALLED (rising なし、ω の末尾 2 窓の中央値比 0.849)、CV 3.4e-6・メッシュ・壁解像 PASS。推力系は継続分 8 点だけだと末尾窓の点数不足 (4 < 8) で機械的に FAIL — 親 20 点 + 継続 8 点の連結 (`retry_doe_024/thrust_history_concat.csv`) で ALL STEADY (η の漸近値 0.972136、最終比 −0.000 %)。台帳に PASS として追記し、ドライバを評価の切れ目で止めて再開した (infill iter 6 から、`add_retry_and_restart.sh`)。継続の run に `thrust_history` の連結が要ることは driver の `--gate-only` の欠け (継続 run の準定常ゲートが親の時系列を見ない) — 記録。**訂正 (2026-10-12、主セッション)**: doe_024 の η_exit 0.97214 を「DOE 最良」と書いたのは主セッションの読み違い — 実際は DOE 40 点中 13 位で Pareto 外 (doe_006・doe_037・inf_03_0・inf_09_1 に支配される)。doe_024 を見た時点で doe_008 (0.97935) が既に出ていた。この誤りは諮問のブリーフにも入っており、ドライバを止めて再開した理由 (最良点がサロゲートから抜ける) は誤った前提だった。誤検出の確定と継続 run の判定 A そのものは前提に依らない。。**完了 (2026-10-12)**: `run_0014_moo_d` 61 評価で PASS 60、Pareto 21 点 (HV 0.73097、L/r_t 13.59〜21.69 で η_exit 0.952〜0.981)。最良 η の推移は case/68 README (V3 の設計時間の材料: 9 評価目で 0.97935、27 評価目で 0.98001、59 評価目で 0.98083)。RISING だけで落ちた件数は 1 (doe_024、欠陥票は 2 件目で起票の規則なので未起票)。解釈 (V3 の比較) は (M) の CFD 評価の後 | (1) **driver**: TP は全点冷間 (`_evaluate_tp` は `warm_from=None`、`--seed-res` を無視)、停止規則 = 同じ fail_class (SUSPECT・RETRYABLE) が 3 評価続く、または DOE の SUSPECT + RETRYABLE が DOE 点数の 20 % を超える (実装済み、`run_bell_cv_gate_tests.py`)。(2) **第一セル厚の判別 A/B** (投入前): 範囲の角 θ_max 45°・θ_e 2°・L 13.595 r_t (成立、(θ_max − θ_e)/L = 3.16 で DOE の成立 57 点の最大 2.31 より厳しい) を `mesh.wall_first_frac` 1.75e-5 (×0.7、スロート側は据え置き) で冷間 1 本 = `case/68.bell_two_mode/run_0012_tp_d_corner_wf175`。**A**: y1+ 最大 ≤ 0.9 かつ超過 0 % → ×0.7 で本番。**B**: それ以外 → 同じ点を 1.25e-5 (×0.5) で `run_0013_tp_d_corner_wf125` (品質 VERDICT と AR 緩和の条件を再確認) → PASS なら ×0.5 で本番。本段で発散したら recipe を諮ってから。(3) **本番**: 新しいディレクトリ (A なら `run_0013_moo_d`、B なら `run_0014_moo_d`)、同じ seed・bounds・LHS 40 + EHVI 10 × 2。(4) **記録**: 新しい doe_000・doe_001 と run_0010/doe_000・run_0011 の η_exit の差を第一セル厚の格子感度として (期待 \|Δη\| ≲ 1e-4)。(5) **キャンペーン中に見るもの**: fail_class の dv 空間の分布 (角に集中するか)、RETRYABLE の失敗段 (soft/mid = IC、本段 = 2 次/SST)、infill が同じ失敗領域を繰り返し提案していないか、推力系の STEADY が崩れる点 (#3f の ω プラトー)。RETRYABLE は再投入せずサロゲートにも入れない — DOE の後にまとめて見て、集中していれば recipe 変更を諮ってから手動で再評価し同じ台帳に `retry_` タグで足す | O (実行) / F (解釈) |
| 8 | 比較の run と判定 | §6 V2・V3。case/68 の README に run 一覧。**V3 の前提 (diagnostician 2026-10-12)**: ピン版 (M) の V2 PASS と、#7 の扱い (δ* 補正済みの物理壁で NS 評価するか、非粘性壁のまま NS で評価して δ* の効きを U_k の δ* 項で持つか) を plan に書いてから (未決なら諮問)。(M) の NS 評価は (D) と同じ格子族・recipe (第一セル 1.25e-5・冷間・本段 60000)。**V3 の列 (codex 2026-10-12)**: (D)・(M 補正なし)・(M 理論式較正の補正あり) の 3 列、補正なしを基準列に残す。同じ L_k で比べ、(D) の Pareto の内挿値は正式判定に使わない (固定 L の再最適化)。較正の NS・Euler 基準・検証・補正形状の再生成はキャンペーン準備費として数える (ガス・Re・熱条件をまたぐ再利用は未検証なので一度きりの開発費にしない)。補正の有無による性能差と U_k は分け、二重計上しない | O (実行) / F (解釈) |
| 9 | 海面作動の比較 | 兄弟 plan の外部領域が受理されてから。§6 V4 に参照・判定量を書いてから回す | F |

## 6. 検証

**判定基準は run の前にここで固定する。** codex diagnose・plan 段と diagnostician (2026-10-11) の指摘を反映した。

- **V0-a (2026-10-11、`case/68.bell_two_mode/v0_moc_gate/` の `v0a_*`) = 判定不能 (標本が退化)** — (D) 成立 8/135、$L_k$ 90・100 % で成立 0、7/8 が $x_m/L$ 0.1 に集中 (diagnostician 2026-10-11 [記録](../../notes/reviews/2026-10-11-bell-v0-gate-diagnostician.md))。
- **V0-b (2026-10-11、`v0b_*`) = 中間 (解釈は確定、diagnostician [記録](../../notes/reviews/2026-10-11-bell-v0b-diagnostician.md))**: 成立 31/60 で有効。n 321 で通過 18・不受理 12・判定不能 1 → 60.0 %。θ_max 違反 2 点の除外・n 161 との不一致 2 点・出口直前の弱い不受理 2 点・marginal 1 点のどの再割当てでも 54〜69 % で分類は変わらない。健全性の確認は「40° 全不受理」が一致、他 2 パターンは不一致 — 原因は (D) の壁で圧縮波の合体点が出口面の前後にあること (弱い不受理群の開始 x/L 0.87〜0.97)。TOP は 5 本とも不受理。言ってよい範囲: 「θ_max 25°・θ_e 10° の (D) は L80・L90 で通過、同じ格子点の TOP は全 L で不受理」。言わない: 「(D) は内部衝撃波が無い」「(D) が TOP より良い」「Rao の台地を含む」(Phase 2 の台地は ε 9・L 5.97 r_t で ε 50 の $L_k$ に対応しない)。
- **V0-c 結果 = A (2026-10-11、`v0c_*`) → V0 を閉じた**: 成立 29/60 (有効)、n 321 で通過 15・不受理 13・判定不能 1 → 53.6 % で中間、n 161 との不一致 2 (V0-b と同じ点)、V0-b の 29 点で合否の反転 1 (L90・32.5°・5°、出口面 x/L 0.998 の交差)、TOP 0/5。読み方を変えても A (`judgement.reading_sensitive: false`)。事前登録どおり (D) の MOO の準備と S1 を並行、④の基部は (D) 先行。
- **V0-c (事前登録、#3d の後に記録を差し替える回帰)**: 同じ 60 点・n 321/161・同じ閾値。**A**: 分類が中間のまま、かつ n 161 との不一致 ≤ 2 → V0 を閉じ、(D) の MOO の準備と S1 を並行、④の基部は (D) 先行。**B**: 分類が変わる・有効条件が破れる・既存 29 点のうち合否が 4 点以上反転 → 再諮問。
- **V0 (前進 MOC ゲートの A/B、CFD なし、事前登録。V0-b として回し直す)**: 対象 = (D) の格子点 ($\theta_{max}$ {20, 25, 32.5, 40}° × $\theta_e$ {5, 10, 15}° × 各 $L_k$ = 60 点、$x_m$ は内部決定、成立条件を満たす点のみ) と、Phase 2 の TOP 最適形 ($\theta_n$ 25.8°・$\theta_e$ 10°) を $\varepsilon$ 50・本計画のスロートに置き直したもの。各点でゲートの合否と内部衝撃波の指標を記録する。**A**: (D) の成立点の 80 % 以上がゲートを通る → ④の基部は (D) から先に使い、(M) の S1 は優先度を下げて続ける (両モード提供の方針なので (M) は取りやめない)。**B**: 通る点が 20 % 以下 → ④の基部は (M) 必須、S1 (または §8-4 の TIC) を最優先に繰り上げ、(D) は③で使う場合に「海面作動・デュアルベル基部には不向き」と設計者に示す。**中間** → S1 と (D) の MOO を並行し、優先順は通過した領域が Rao の台地を含むかで決める (結果を添えて上位に諮る)。**有効条件**: 各 $L_k$ に成立点 ≥ 3 かつ合計 ≥ 15、満たさなければ「判定不能」(次は最小ジャークの選択子が不等式を見ない件を QP 化として諮る)。主設定は Hall n 321、n 161 で合否が一致すること。C⁻ 間隔比の最小 < 0.2 の通過点は「marginal」と注記して分類には入れる。**ゲートの健全性の確認** (既知のパターン): $\theta_{max}$ 40° は全 $L$ で不受理、L60 は 32.5° でも不受理、32.5°・L ≥ 70 % は通過 — 崩れたら「ゲートが帯内の $x_m$ に敏感」として報告する。**TOP の書き方**: 「前進 MOC ゲートは TOP を $L_k$ 5 本すべてで不受理、最初の C⁻ 交差は x 3.9〜5.0 $r_t$・r 0.86〜0.96 $r_t$、case/40 で forge・SU2 が確認した TOP の内部衝撃波と向きが一致」まで (本ノズルで衝撃波を CFD で確認したとは書かない、交差後の間隔比を強さとして書かない、真空 $\eta$ に結びつけない)。§4.5 の TOP 3 点の CFD で衝撃トレースの立ち上がり x とゲートの交差開始 x を併記する (±1.5 $r_t$ なら整合、合否にしない)。どちらでも V3 の $\eta$ 比較は「同等幅内」を既定の予想として併記し、方式選択の根拠にしない。

- **V1 (M の生成器)**:
  - 回帰: CPG・円弧 kernel で古典 Rao を再現 (Rao チャートの $\theta_n,\theta_e$ と ±1°、チャートの比熱比に合わせる)。TP 生成器を定比熱にすると CPG 生成器と 1e-8 で一致。
  - **(diagnostician 2026-10-11 の T1–T6 で置き換え)** ~~T1a CPG・円弧 kernel・**チャートの R_t (0.382) と γ** で θ_n・θ_e が ±1°。~~ **T1a は記録 (合否なし) に格下げ** (diagnostician 2026-10-11: チャートの γ は出典不明 — case/29 `mesh/make_nozzle.py:41-48`・README に記録なし。R_t 0.382 は Hall が精度低下域 [`transonic.py:126-129`] で前進 MOC が通らない)。#6 の中で case/29 の R_t 0.382 の TOP 壁の Euler 1 本を凍結源に CPG γ {1.2, 1.3, 1.4} で回し、チャート値がその帯のどこにあるかを記録する。**T1a′ (必要条件)**: 同じ $L_k$・ε で前進 MOC ゲートを通った (D) 壁の MOC 推力 ≤ (M) の MOC 推力 (破れたら停留点が最大でない証拠)。T1b 平面縮約が `rao_planar.theta_lip` と 1e-8。~~T1c 自由 ε で得た (L, ε) を固定 ε で解き直してビット一致。~~ **T1c (2026-10-12 改)**: 自由 ε の解を固定 ε = ε_free で解き直した差が根の許容 XTOL 1e-12 以内 (根の方程式が違うのでビット一致は期待しない。観測 ≤ 3.3e-13 — `case/68.bell_two_mode/rao_m/rao_t1c.csv`)。T1d 定比熱の擬似種で CPG 生成器と 1e-8。**T2 包絡線の恒等式** $dF/ds_D=2\pi[\Phi_3\,dL/ds_D+F(E)\,dr_E/ds_D]$ が相対 1e-5 で M 格子を倍にして収束。T3 保存: 閉合 $|W_m-E|$ と推力の 2 会計 (制御面 / 初期線 + kernel 壁 + BE 壁) の差 ≤ 1e-4 で n とともに減る (「DE の流量 = スロート流量」は構成で恒等なので試験にしない)。T4 停留性 (下の項目を具体化) **— 2026-10-12: Hall n ≤ 321 の評価器 (`forward_thrust`) では判定不能** (c₁ε₀ の床は壁圧の求積が支配: スプラインと台形の会計の c₁ε₀ の差が L100 で 1.22e-6 / 4.6e-7 / 3.1e-7 [n_eval 161/241/321]、線 0.05|c₂ε₀²| ≈ 1e-7 の 3〜12 倍。c₂ < 0 は 138/138)。線は変えず、#4 の「1 回修正」を評価器に使う = 整合な前線の再標本化で初期線を N 1281 まで (§5.1 #5 (h))。記録項目に ε* = |c₁/(2c₂)| [r_t]・懸かる推力 c₁²/(4|c₂|)/F_1D・スプラインと台形の差・c₂ の符号、`ctrl_minus` (別の (L′, ε′) の Rao 輪郭) は正対照として足す。負対照 `ctrl_plus` は ε が +1.6 % 変わる (σ の族で x = L・√ε を両方通るのは (M) だけ) — 注記して採用、比の線は負対照の c₁ が床の 10 倍以上ある k1 でだけ使う。登録の内容: 壁 BE に $\eta_k$ ($k$=1..3) を $\pm\{1,2,4\}\times10^{-3}r_t$ で足し、推力 $F(\epsilon)$ を 3 次で当てはめ $|c_1\epsilon_0|\le0.05|c_2\epsilon_0^2|$、負対照 (D を 2 % ずらした非 Rao 線) の $c_1$ の 5 % 以下、n で減る。~~T5 C2: $|\tilde g|\le10^{-8}$ rad/$r_t$、$\kappa_W(B^+)-\kappa_B$ は n 161/321 の変化の 3 倍以内。~~ **T5 (2026-10-11 改)**: C2 は課さない — $\tilde g$ と $\Delta\kappa_B=\kappa_W(B^+)-1/R_t$ を n 161/241/321 で記録し、n の変化の 3 倍以内に収まること。**T3 (2026-10-11 改)**: T6 の系列で判定 — 生産の n 321 で推力 2 会計の差 ≤ 1e-4 かつ n 161/241/321 で単調に減り次数 ≥ 1 (閉合 |W_m − E| も同じ)。n 161 の 1.2e-4 は記録で PASS と書かない。321 で > 1e-4 か次数 < 1 なら FAIL として諮る。~~T6 分解能: n 161/321/641 で θ_max・θ_e・F が単調に変わる。~~ **T6 分解能 (2026-10-11 改、対象は円弧 kernel の L100・L70 — 成立域の両端)**: n {161, 241, 321} で θ_max・θ_e・F が単調に変わり Richardson 外挿が取れる。加えて n 321 固定で dM を半分・m・q を倍にし、kernel 場と DE/Goursat の分解能を分けて見る (n 641 は Hall のガードで不可、CFD ピンの後に記録として足す)。**T1b・T2 は折り返し (θ = μ) を跨ぐ D の組を必須にする** (例: `synth_bd` で M_B 2.2・θ_B 32°、折り返し点で dθ/dr = −tan θ/(2r) を点検)。
  - (旧) 保存: 制御面 DE を通る流量がスロートの流量と相対 1e-4 以内。制御面の推力と「入口 + 壁」の推力が相対 1e-4 以内。
  - **TP の停留性**: $c_p$ が実際に変わる条件 (基準ケースの EXH) で、拘束 ($\dot m$、$L$、$\varepsilon$) を保つ壁の摂動族を作り、MOC で求めた推力の一次変分が 0 であること。摂動振幅 $\epsilon$ を掃引して推力変化が $O(\epsilon^2)$ になること、MOC を細分しても結論が変わらないこと、二次変分の符号 (負なら局所最大、そうでなければ「停留」とだけ呼ぶ) を記録する。任意の形状に勝つことを最適性の証明として使わない。
  - スロート流量は CFD ピンの $m^*$ (等エントロピー写像の流量) と定義し、凍結源の CFD 流量との差を別に記録する。
  - C2 接続: (T5 に置き換え) 旧「$|g(\kappa_B)|\le10^{-6}$」は流線の 2 階微分では達成できない。
- **V2 (評価器の整合)**: (M) の非粘性壁を §4.1 の Euler 設定で回し、推力が MOC の値と相対 0.2 % 以内 (粗い整合確認)。
  - **結果 (Hall 版の壁、2026-10-12、記録 [diagnostician](../../notes/reviews/2026-10-12-bell-m-v2-euler-diagnostician.md))**: `case/68.bell_two_mode/run_0015_tp_m_euler_l100`・`run_0016_tp_m_euler_l70` (AWS B、全ゲート PASS) で F_Euler/F_MOC − 1 = +0.287 % (L100)・+0.288 % (L70) → **線 0.2 % の外で FAIL**。内訳: ṁ_Euler/ṁ_Hall − 1 = +0.295 % (Euler ṁ/ṁ_1D 0.99263 vs Hall n 321 の線 0.989706)、(F/ṁ)_Euler/(F/ṁ)_MOC − 1 = −8.6e-5 / −7.3e-5。原因は初期線 (Hall、R_t 1) の流量 — §4.4 の「CFD ピン既定・Hall は検算」の検算として機能した。「MOC の膨張が 0.3 % ずれている」とは書かない。Hall 版の V2 は記録、**合否が立つのは CFD ピン版の V2**。
  - **ピン版 V2 の事前登録 (2026-10-12)**: (主) \|F_Euler/F_MOC − 1\| ≤ 0.2 %、(内訳 1) \|ṁ_MOC(ピン線)/ṁ_Euler − 1\| ≤ 2e-4、(内訳 2) \|(F/ṁ)_Euler/(F/ṁ)_MOC − 1\| ≤ 5e-4。**A** 主が通る → PASS (内訳も記録)。**B** 主が外れ内訳 1 は通る → 膨張部の不一致 (鋸歯の帯・Euler 格子・軸近傍の衝撃の扱い) として再諮問。**C** 内訳 1 が外れる → ピンの写像か凍結源の照合の問題 (#6 の照合に戻る)。
  - **結果 (ピン版、2026-10-12) = A → V2 PASS**: `case/68.bell_two_mode/run_0018_tp_m_euler_l100_pin` (AWS B、全ゲート PASS・CV −5.7e-7・残差 STALLED 併記・推力系 ALL STEADY) で (主) F_Euler/F_MOC − 1 = +1.70e-5 (線 0.2 %)、(内訳 1) ṁ_MOC/ṁ_Euler − 1 = −7.3e-6 (線 2e-4)、(内訳 2) (F/ṁ)_Euler/(F/ṁ)_MOC − 1 = +9.7e-6 (線 5e-4)。Euler η_exit 0.9903618・ṁ/ṁ_1D 0.9926296、MOC η_cs 0.9903449・0.9926224。Hall 版の +0.29 % は初期線の流量だけが原因だったことが確かめられた。記録: runner の Euler の metrics に MOC との差 (`thrust_exit_rel_diff_cs` 等) が空で出た — `--gate-only` の経路で `moc_reference` が拾われていない欠け (上の値は WALL_DESIGN.json と thrust_history.csv から手で計算)。
- **V3 (比較、事前登録の分類)**:
  - (V3 は 2 モードのどちらかを選ぶ判定ではない。両モードを提供し、結果は設計者向けの判断材料の表にする — 親 §3-16。)
  - **設計時間 (判断材料の主項目)**: 「十分良い形」= 制約 (ゲート・ゲート (i)〜(v)) を満たし、各 $L_k$ で $\eta\ge\eta_{best}(L_k)-\delta_{good}$ ($\eta_{best}$ は両モードを通じて得た最良、$\delta_{good}$ は仮 3e-3 = 0.3 % — ユーザ確認待ち)。両モードと基準 (TOP チャート + CFD 1 評価) について、十分良い形に達するまでの評価回数と実時間 (同じ GPU、キャンペーンごとの準備を含む) を表にする。合否ではなく比較の表として出し、一度きりの開発費は別欄。
  - 各 $L_k$ の $\Delta\eta(L_k)$ の不確かさ $U_k$ (plan レビュー M4、δ\* 項は diagnostician): **格子** = 各 $L_k$ で採用した 2 形状の組を系統格子 3 水準で評価し、$\Delta\eta$ **そのもの**の格子収束から推定 (各モードの格子差の最大値は使わない)。**反復** = 判定区間を延長した履歴と漸近値の推定から。**探索** = (D) の固定長再最適化を独立の初期点・seed で繰り返し、予算を追加したときの改善量から。探索の上限を保証できない場合、結論は「登録した予算で得た設計どうしの比較」に限定する。**δ\*** = (M) だけが背負う δ\* 補正・$k_f$ 較正・拘束への射影の誤差。$k_f$ の検証形状 (較正に使わない形状) での RANS 抽出との差から推力への寄与を見積もる。4 項は保守的に足し合わせる。
  - **$|\Delta\eta|-U_k>1.5\times10^{-3}$ → 「差あり」** (符号で優位なモードを記録)。**$|\Delta\eta|+U_k\le1.5\times10^{-3}$ → 「1.5e-3 の工学的同等幅の中で同等」**。他は判定不能。1.5e-3 は Phase 2 の Rao 照合の判定帯と同じ値を工学的な同等幅として使う。
  - 併記 (合否なし): 形状が作れない点の率、総費用、内部衝撃波の指標、$g$ と当てはめの残差。
- **V4 (海面作動)**: 兄弟 plan の外部領域の受理後に、参照・判定量 (剥離位置、海面の推力、内部衝撃波) をここに書いてから回す (未登録)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-11` | [`2026-10-11-thruster-nozzle-child-plans-diagnose.md`](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) (兄弟 plan と一括) | NO-GO, C1/M4 (本 plan 分: 共通 C1・B1〜B4) | **C1 採用** (兄弟 plan §4.4 の閉じた CV を共通定義に)。**B1 採用** ($g(\kappa_B)$ を明示、各評価で B・D を解き直す、当てはめの端点拘束を一般化 — `wall_axismach.py:1061` の始点固定を確認、$\kappa$ と $r''$ を区別)。**B2 採用** (凍結 TP の $h,s$ で導出、TP の一次変分試験、CPG は回帰)。**B3 採用** (相似倍率では $\varepsilon$ が不変なのは自明 — 拘束付き当てはめで物理壁を確定、「完成した設計手順の比較」、総費用を比較)。**B4 採用** (Euler を明示、共通 $L_k$ で比較、$U_k$ の構成、「差なし」→「同等幅内で同等」)。前提への異議 (`viscMethod: 1` は空気の Sutherland で定数ではない) を採用し §3 を訂正。0.1〜0.2 % の予測は根拠不足として撤回 |
| diagnose (codex) | `2026-10-12` | [`2026-10-12-bell-m-deltastar-calibration-diagnose.md`](../../notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md) | C0/M5/m1: 抽出を一定 x に・帯位置の感度・δ_r の定常性・効きを測る・δ_r は半径方向・contur_v2 を使う | **全件採用** (M1 の 31 % を主セッションで再現)。§5.1 #7・#8。以降の諮問は codex (ユーザ指示 2026-10-12) |
| diagnostician | `2026-10-12` | [`2026-10-12-bell-m-deltastar-route-diagnostician.md`](../../notes/reviews/2026-10-12-bell-m-deltastar-route-diagnostician.md) | δ* の扱い: (C) 生産・(A) 棚上げ・(B) 昇格先を推奨 | **一部採用** — ユーザ決定で補正は設計者の選択肢、方法は理論式 + NS 較正 ((A) の棚上げは不採用)。§4.4・§5.1 #7 |
| diagnostician | `2026-10-12` | [`2026-10-12-bell-m-v2-euler-diagnostician.md`](../../notes/reviews/2026-10-12-bell-m-v2-euler-diagnostician.md) | Hall 版 V2 は FAIL (原因は初期線の流量)、ピン版 V2 を事前登録、凍結源の前に 3 点、V3 の前に #7 の扱い | **全件採用**。§6 V2・§5.1 #6・#8 |
| diagnostician | `2026-10-12` | [`2026-10-12-bell-rao-t4-stationarity-diagnostician.md`](../../notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md) | T4 判定不能 → 線は変えず評価器を 1 回修正 (A/B 事前登録)、ゲートの所見の文言、T1c の言い換え | **全件採用** (スプラインと台形の c₁ の差を `rao_t4.json` から主セッションで再計算)。選択肢 (a) 格下げを並べたことへの異議を受け入れる。§4.5・§5.1 #5 (h)・§6 V1 |
| diagnostician | `2026-10-11` | [`2026-10-11-bell-rao-mainbranch-c2-diagnostician.md`](../../notes/reviews/2026-10-11-bell-rao-mainbranch-c2-diagnostician.md) | C2 の根なし → (M) は円弧・B で C1、受理は det_D > 0、T3 は T6 の系列で | **全件採用** (部品試験を主セッションで再現、C1 はユーザ了承)。§4.4・§5.1 #5・§6 V1 T3・T5・T6 |
| diagnostician | `2026-10-11` | [`2026-10-11-bell-d-moo-warmstart-diagnostician.md`](../../notes/reviews/2026-10-11-bell-d-moo-warmstart-diagnostician.md) | 全点冷間・新ディレクトリは採用、第一セル厚は角の 1 本で判別、停止規則 | **全件採用** (`interp_field.py:212-216` の生座標の最近傍と `_warm_for` を主セッションで確認)。§4.1「起動」・§5.1 #8a・#8a′ |
| diagnostician | `2026-10-11` | [`2026-10-11-bell-s1-generator-holes-diagnostician.md`](../../notes/reviews/2026-10-11-bell-s1-generator-holes-diagnostician.md) | 前回の |c| 規則を撤回、DE を r パラメータに、T1a 格下げ・T1a′ 追加、T6 改、θ_E ≤ θ_B | **全件採用** (c = cos(θ−μ) の恒等・θ = μ での dθ/dr の一致を主セッションで再導出、チャートの γ 不明を case/29 で確認)。§4.4・§5.1 #5・#6・§6 V1・§8-4 |
| diagnostician | `2026-10-11` | [`2026-10-11-bell-v0-gate-diagnostician.md`](../../notes/reviews/2026-10-11-bell-v0-gate-diagnostician.md) | V0-a 判定不能、x_m 内部決定 | **全件採用** (x_m の帯を主セッションで再確認)。§4.3・§5.1 #3c・§6 V0 |
| diagnostician | `2026-10-11` | [`2026-10-11-bell-2b-gates-diagnostician.md`](../../notes/reviews/2026-10-11-bell-2b-gates-diagnostician.md) | ゲート組み替え、継続 run 1 本 | **全件採用** (cv の交差・hoop 修正の欠落を主セッションで確認)。§4.1・§5.1 #2c〜#2f |
| diagnostician | `2026-10-11` | [`2026-10-11-thruster-strategy-diagnostician.md`](../../notes/reviews/2026-10-11-thruster-strategy-diagnostician.md) (3 plan の横断) | 推奨 1 (ゲート先行)、食い違い 6、危険 2 | **全件採用** (根拠を主セッションで確認: `moc_diagnostics.py:88-91,201`、`moo-loop.md:102-103`)。§4.5・§5・§5.1 #3b・§6 V0・V3・§8-4 に反映。ゲートの初期線は Hall (主セッションの判断) |
| plan | `2026-10-11` | [`2026-10-11-tooling-nozzle-bell-two-mode-plan.md`](../../notes/reviews/2026-10-11-tooling-nozzle-bell-two-mode-plan.md) | GO-with-changes, C0/M8/m1 | **全件採用** (根拠のコード行は主セッションで確認: `extract.py:38` CPG 圧力、`driver.py:84` 定比熱の流量、`driver.py:173` の DIVERGED 以外受理、`polish.py:41` の再 prepare、`deltastar_integral.py:208` の CPG 換算)。M1 → §4.1 TP 化の範囲・§5.1 #2、M2 → §4.1 推力 (全境界の $\Pi_x$、出口の粘性応力)、M3 → §4.1 ゲート・polish の同一メッシュ継続・`outflow`、M4 → §6 V3 の $U_k$、M5 → §4.4 S1 ゲート・§6 V1 の停留性・$m^*$ の定義、M6 → §4.2 最終の物理壁の拘束、M7 → §4.4 δ\* の TP 化・$k_f$ の較正/検証の分離、M8 → §4.2 $L_U$・$L_{pipe}$ と単調条件・§4.3 と §5.1 #3 の基底の固定、m1 → §5.1 #2 の配線回帰。段階方式 (共通評価器 → D → S1 を通した M) を §5 に採用 |

## 7. 影響範囲

- 新規: `design/forge_design/geometry/rao_axisym.py` (仮)、`geometry/wall_bell_c2.py` (仮)、`case/68.bell_two_mode/` (新設)
- 変更: `evaluate/runner.py` (TP・CEA 組成・`viscMethod: 2`・低 Re 格子・Euler 設定)、`geometry/wall_axismach.py` (当てはめの端点拘束の一般化。既定の挙動は不変)、`feedback/deltastar_integral.py` ($\mu(T)$)、`feedback/cfd_initial_line.py` (照合の形)、`probdef.py`、`design/CAPABILITIES.md`、`methods/design/overview.md`
- ソルバは変更しない予定

## 8. 未確定事項

1. **寸法 spec** — **仮値で plan 段レビューに出す (2026-10-11)**: スロート半径 $r_t=0.14$ m (推力 1000 kN 級、$A_t=F/(p_c C_F)$)、$R_t=1.0\,r_t$ (0.4 $r_t$ [TOP の Rd] は感度として 1 組)、収縮比 3、比較の長さ $L_k$ = 等面積比 15° 円錐の長さの 60・70・80・90・100 % (1.903〜3.172 m)、$L_U$・$L_{pipe}$ は単調条件の内側で未定。ユーザが別の値を決めたら差し替える
2. (M) の kernel の曲線族と、B・D の反復の組み方 (S1 で導出)
3. 輸送物性の出所 (CEA の trans.inp の種がそろうか、lump の扱い) — 2a で H2・O2・H2O・OH・H・O に決着 (H2O は custom、他は cea)
4. **S1 が閉じない場合の代替 (事前登録、diagnostician R2)**: §5.1 #4 の導出と S1 ゲートが上位への諮問 2 回で閉じなければ、(M) の族を **TIC** (①の semi-perfect 軸マッハ/理想ノズルのチェーンで作って切り詰める。導出が要らない) に切り替えて先に使う。TIC は同じ長さで TOC より理論的に推力が低いので $\eta$ の比較は (M) に不利に偏るが、方式選択は V0 のゲートで行うので影響しない。検証データ (TC1A・P6.2 の基部) はいずれも TIC
   - **S1 の諮問は 2 回を消化** (2026-10-11: 導出・|c| 規則の撤回)。§5.1 #5 の r パラメータの修正後に、円弧 kernel・EXH・R_t 1・ε 50・n 161 で L_k 60〜100 % の 5 本を回し、**2 本以上で主枝 ($\theta_E\le\theta_B$) の根が括れなければ、追加の諮問なしで TIC に切り替える**。旧規則 (|c| > 1 は解なし) での 4/5 根なしは数えない。修正後の結果に「もう 1 回直してから」の例外は足さない

## 9. 変更ログ

- `2026-10-12` — (M) の NS 4 本 (補正なし) が全ゲート PASS。δ* の較正の設計を codex に諮り全件採用: 抽出を一定 x に直してから、contur_v2 の定数 1 個を L100 で較正・L70 で検証、補正の効きは壁圧の誤差で測る、V3 は 3 列。以降の諮問は codex (ユーザ指示)。
- `2026-10-12` — ユーザ決定: (M) の δ* 補正は設計者の選択肢、方法は理論式を NS 結果で較正。次は η_M0 (非粘性壁の NS) → 較正の設計を諮問。
- `2026-10-12` — CFD ピン版の (M) L100 で V2 PASS (推力の差 +1.7e-5、事前登録の判定 A)。Hall 版の +0.29 % は初期線の流量だった。
- `2026-10-12` — (M) の W1 (runner の moc_rao・TP の Euler) と W2 (Euler 2 本)。V2 の Hall 版は +0.29 % で線の外 — 原因は Hall 初期線の流量 (流量あたりの推力は 1e-4 で一致)。ピン版の V2 を事前登録して W3 へ。
- `2026-10-12` — (D) の MOO `run_0014_moo_d` 完了 (61 評価、PASS 60、Pareto 21 点)。doe_024 を「DOE 最良」と書いた誤りを訂正 (13 位・Pareto 外)。
- `2026-10-12` — **S1 を閉じた**: T4 の follow (k2・k3・L70・負対照) が全て線の内側。(M) は「指定した kernel と拘束の下の停留輪郭」(ε* ≤ 3e-5 r_t) として比較に入る。次は #6 CFD ピン・V2 Euler 整合・#7 物理壁。
- `2026-10-12` — T4 の (h) A/B は判定 A (c₁ε₀ 2.6〜3.0e-8 < 線 9.9e-8、2 会計 3.4e-9 で一致)。follow (k2・k3・L70・負対照) の後に S1 を閉じる。Hall の初期線の鋸歯の帯を欠陥票に。
- `2026-10-12` — (M) の (g): T1a′・T1c は満たし、T4 は評価器の床で判定不能。線は変えず評価器を 1 回修正する A/B (h) を事前登録。(M) 壁は前進 MOC で DE の下の軸近くに内部衝撃波 (推力会計に影響なし) — §4.5 に所見。
- `2026-10-12` — (M) 生成器 (e)(f): 受理条件 det_D > 0・円弧が既定。T6 で n 161/241/321 が単調、T3 は事前登録の基準を満たす (n 321 で推力 2 会計の差 4.1e-5・2.9e-5)。次は T4 停留性・T1a′。
- `2026-10-11` — (M) の C2 を撤回 (kernel の終端曲率は D での判定量をほとんど動かさず根が無い、部品試験で確認)。ユーザ了承のうえで (M) は古典どおり円弧 kernel・B で C1 (曲率の跳びは記録)。主枝の受理規則を det_D > 0 に。L_k 5 本は事前登録の A (70〜100 % に根、60 % は作れない)。
- `2026-10-11` — (D) の MOO の warm start が本段で発散 (doe_001)、冷間の A/B は PASS → diagnostician の判断で TP は全点冷間に (driver を変更)、停止規則を追加、第一セル厚は範囲の角の 1 本 (run_0012) で決めてから新しいディレクトリで回し直す (#8a′)。
- `2026-10-11` — (M) の生成器で S1 の構成手順の誤り (DE の |c| > 1 を解なしとする規則) が見つかり、diagnostician が撤回 → DE を r パラメータに、T1a を記録に格下げ・T1a′ を追加、T6 を n 161/241/321 に、§8-4 は修正後の 5 本で判定 (諮問 2 回は消化)。(D) の MOO は warm start の化学種の署名の欠け (共用 `forge_species.py`) で 18 点を空費 → 修正して同じキャンペーンで再開。
- `2026-10-11` — #2c'・#2d・#3f 完了: ゲートの定数を driver に、壁法線の近壁層で y1+ ≤ 1 (AR ≤ 5000)、(D) を runner で評価。run_0008 (TOP)・run_0009 ((D) の 1 点) とも全ゲート PASS。次は (D) の MOO と (M) の生成器 (#5、実装中)。
- `2026-10-11` — S1 の導出を記録 (diagnostician): 凍結 TP でも Rao の制御面条件は古典と同形で γ-free、固定 (L, ε) で縁条件は落ちる、κ_B・θ_B・s_D の入れ子の根、C2 は g̃、§6 V1 を T1–T6 に置き換え。#5 (生成器) に進む。
- `2026-10-11` — #3d 完了、V0-c は事前登録の判定で A → **V0 を閉じた**。
- `2026-10-11` — 本題に進む (ユーザ了承): #3d の後に V0 を閉じ、#2c' (ゲートの定数)・#2d (壁法線の近壁層)・#3f (runner への (D) の組み込み) と S1 (TP の Rao の導出) を並行。③の真空の比較ではゲートは記録のみ。
- `2026-10-11` — V0-b を中間で確定 (diagnostician、再割当てで 54〜69 %)。#3d (成立条件の穴・区間の端)・#3e (C⁺ 指標)・V0-c の事前登録・MOO の制約は連続の余裕、を追加。
- `2026-10-11` — ユーザ決定「修正は必ず取り込む」: 近軸 hoop 修正 (0b4dff4e) を、修正を含む `feature/nozzle-wall-fit-and-pipeline` の merge で取り込む (#2e)。同じ merge の `wall_normal_layer` を #2d に使う。
- `2026-10-11` — #2c の事前登録の検証は **A** (run_0005: $\eta_{exit}$ は予測漸近値から 1.1e-5、継続区間の $|cv|\le4.25\times10^{-5}$)。組み替えを driver に反映する #2c' を追加。
- `2026-10-11` — #2b のゲートを組み替える方針を採用 (diagnostician): 残差は記録、収束は最終 |cv| ≤ 1e-4 + 推力系 STEADY 1.5e-4、η は出口面、壁法線の近壁層で AR ≤ 5000、pRef は置かない。#2c (継続 run で前提を検証、事前登録)・#2d (壁法線層)・#2e (hoop 修正の欠落)・#2f (CPG の cv) を追加。
- `2026-10-11` — V0-a を判定不能 (標本退化) として記録 (diagnostician)。(D) の $x_m$ を内部決定に (設計変数は 3 個)、V0 の有効条件・marginal・健全性の確認・TOP の書き方を事前登録、#3c (回し直しの準備) を追加。
- `2026-10-11` — ユーザ決定: 壁と出口面の推力の一致を収束判定の一つに、AR は 5000 以下を基準 (§4.1)。
- `2026-10-11` — #2b・#3・#3b を実装 (implementer、主セッションで試験と主要な数値を再確認)。V0 (CFD なし) は暫定で「中間」((D) 成立 8/135・通過 3/8、TOP 0/5)。ゲート (i)(ii)(v) が基準ケースを通さない件と、V0 の解釈・(D) の族・Hall 初期線の件を diagnostician に諮問中 (結果を見てから §4・§6 を直す)。
- `2026-10-11` — ユーザ判断 (親 §3-16): 2 モードとも提供し設計者が選ぶ (`geometry.bell_mode`)。比較は選ぶための判断材料で、主項目は「十分良い形に達するまでの評価回数・実時間」($\delta_{good}$ 仮 0.3 %)。V0 の A/B は優先順だけを変え、どちらのモードも取りやめない。
- `2026-10-11` — diagnostician (Fable) の横断レビュー [記録](../../notes/reviews/2026-10-11-thruster-strategy-diagnostician.md) を全件採用: 第一の判定量を前進 MOC ゲートに (η は拒否権)、順序を「(D) の壁 → ゲートの A/B (V0、CFD なし) → S1 と MOO の順を決める」に、共用モジュール #3b (所有はこの系列)、$U_k$ に δ\* 項、作動点を高度のリスト (真空は `p_a: 0`)、S1 が閉じない場合の TIC の代替をトリガ付きで §8-4 に。
- `2026-10-11` — §5.1 #2a 完了 (TP の排気・1D 理想値・runner の TP 経路。1D と CEA 凍結の差 ≤ 2.1e-5、CPG 経路ビット不変、起動確認 run_0001)。2b への引き継ぎ 5 件を #2b に追記。
- `2026-10-11` — codex plan 段 (GO-with-changes, C0/M8/m1) を全件採用して改訂: TP 化を評価経路全体 (推力抽出・理想推力・チョーク流量・IC・壁平滑化・ドライバ) に、推力を全境界の運動量の式に (出口の粘性応力)、ゲートを PASS + 推力系の STEADY に・polish を同一メッシュの継続に・真空は `outflow`、$U_k$ を $\Delta\eta$ 自体の 3 水準格子収束・反復の漸近・探索の再現性に、S1 を独立ゲートに、最終の物理壁の拘束と射影の残差、δ\* の積分法の TP 化、$L_U$・$L_{pipe}$ と単調条件。段階方式。
- `2026-10-11` — codex diagnose (NO-GO) を反映して改訂: 推力を兄弟 plan の閉じた CV に、Euler 設定の明示、`viscMethod: 2` (燃焼ガスの輸送物性)、TP の Rao を凍結 TP の $h,s$ で導出し一次変分で検査、C2 接続を符号付き $g(\kappa_B)$ の根探しに、当てはめの端点拘束の一般化、δ\* 後の物理壁を拘束付き当てはめで確定 (相似倍率は撤回)、比較を共通 $L_k$ で行い不確かさ $U_k$ を構成、「完成した設計手順の比較」と総費用。
- `2026-10-11` — 起票 (親計画 §5.1 #3・#4)。§4・§6 は上位に諮る前の下書き。
```

## 参考: `notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md`

```
# codex 諮問 (diagnose): bell-m-deltastar-calibration

- **brief**: [`notes/reviews/briefs/2026-10-12-bell-m-deltastar-calibration.md`](../../notes/reviews/briefs/2026-10-12-bell-m-deltastar-calibration.md)
- **plan**: [`plans/active/tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md)
- **date**: 2026-10-12
- **commit**: `939d1b64` (feature/thruster-nozzle-design)
- **codex**: effort `high`, 5.6 min, rc=0
- **結論**: **L100 の抽出器を一定 x の断面と実効的な帯位置感度に対応させ、較正に使える δ_r(x) を確定する。**
- **extra**: `case/68.bell_two_mode/README.md`, `notes/reviews/2026-10-12-bell-m-deltastar-route-diagnostician.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0件、Major 5件、Minor 1件）。**次は較正や補正壁の NS 評価より先に、δ* 抽出を修正・受理することを推奨する。**

| 項目 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 現在の δ*(x) を較正データに使う | **要再検証／Major M1** | 保存済み配列は `remap_constant_x=False` の平滑化結果と最大差 6.5e−16 で再現した。しかし NS 格子の列内 x 幅は最大 0.282 r_t。一定 x に再標本化すると拡大部の抽出値は最大 **31.35%**変わる。既存の `remap_constant_x=True` を使う。[抽出コード:74](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:74)、[呼出箇所:441](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:441) |
| 「全1201断面 ok」なので抽出は受理済み | **却下／Major M2** | 同じ物理壁なので `delta_in=0`。帯の開始位置は `max(k_band·delta_in, band_min_frac·rw)` となり、`k_band` を変えても変わらない。実測でも**全1201断面で感度2系列が主系列とビット一致**した。対案は、選択された実際の帯位置 `y_b` を変え、境界層外で抽出値が安定するか検査すること。[帯位置:232](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:232)、[感度計算:579](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:579) |
| 推力系の STEADY を δ* の定常性にも適用する | **却下／Major M3** | `run_0019` の準定常判定対象は推力・η・ṁだけ。δ*(x) の時系列はない。残差判定も `NOT CONVERGED (stalled/plateau…)` で、ω はまだ下降している。対案は、修正した抽出器による δ_r の時系列を別途 `check_quasisteady.py --series-csv` で判定すること。 |
| 「補正しない損は2次で約1e−4」 | **却下／Major M4** | 停留性試験は固定された ṁ・L・ε の下での非粘性壁摂動に限定される。粘性の排除効果、出口の有効半径変更、摩擦変化まで一次項が消える証拠ではない。[plan:167](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:167)。対案は補正による Δη を符号も含めて測定し、事前に小さいと決めないこと。 |
| 理論式＋NS較正を選択肢として提供する | **採用** | ユーザ決定を維持する。基礎式は既存の **CONTUR `contur_v2`**、較正はまず無次元の定数1個とする。平板相関や任意の x 依存補正には進まない。 |
| 積分法は CPG なので新たな TP 化が必要 | **却下／Minor m1** | 現コードには、同じ μ(T) の注入、エンタルピーによる温度分布・回復温度、TP の d ln u_e/dx を使う `contur_v2` がある。必要なのはベルへの接続と検証。[積分法:134](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:134)、[回復温度:186](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:186)、[運動量式:277](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:277) |
| 抽出 δ* をそのまま法線オフセットに使う | **却下／Major M5** | 抽出値は**一定 x における半径方向等価厚さ δ_r**であり、法線厚さではない。[抽出定義:410](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:410)。対案は固定 x の半径補正に統一すること。法線変位を使うなら一次近似で δ_n = δ_r cos θ とし、x の移動も扱う。 |

結論: **L100 の抽出器を一定 x の断面と実効的な帯位置感度に対応させ、較正に使える δ_r(x) を確定する。**

第 1 仮説: 現在の配列で理論式を較正すると、格子列の取り違えによる系統誤差を較正係数へ吸収する。 **確度: 高**

  根拠: 同じ最終場を読み、断面の扱いだけを変更した再計算で次を得た。値は δ_r/r_t。これは**最終スナップショットの診断値であり、定常値とは認定していない**。

| x/r_t | 従来の列による抽出 | 一定 x に補間した抽出 |
|---:|---:|---:|
| 0.5 | 0.0017765 | 0.0021987 |
| 1.0 | 0.0031796 | 0.0041637 |
| 2.0 | 0.0061992 | 0.0079916 |
| 4.0 | 0.0137714 | 0.0164927 |
| 8.0 | 0.0325620 | 0.0365670 |
| 16.0 | 0.0822393 | 0.0861515 |
| 出口 22.6575 | 0.1344135 | 0.1344135 |

  根拠 run は `case/68.bell_two_mode/run_0019_tp_m_ns_l100/` と Euler 基準 `case/68.bell_two_mode/run_0018_tp_m_euler_l100_pin/`。出口面は平面なので出口値は変わらず、**出口の約1.9%は今回の再標本化でも消えない。一方、較正に重要な途中の分布は変わる。**

  反証条件: 独立に一定 x 断面を構成した抽出が旧配列を再現し、既存の再標本化側に誤りが見つかること。出口だけの一致は反証にならない。

第 2 仮説: 適応帯が境界層の縁を十分外しておらず、局所的な抽出不足が残る。 **確度: 中**。一定 x に直した後も、x≈0.5 で実際の帯位置を1.953倍にすると δ_r は約9.3%増える。`band_rms` が小さいだけでは排除できない。

第 3 仮説: 推力が定常でも境界層分布には時間変化が残る。 **確度: 低・未確認**。手元には各 run の最終場1枚しかなく、δ_r の時系列による判定ができない。

判別 A/B: **変更は `remap_constant_x=False → True` の1点、追加計算は0 step。** 同じ NS・Euler 最終場から拡大部全域の δ_r を再抽出する。

- A：全域で差が既存の抽出感度目安5%以内なら、断面取り違えを支配的候補から外す。
- B：局所的に5%を超えるなら、旧配列を較正から外して抽出修正を先行する。

本レビューでこの比較を実施し、**最大31.35%で B に相当する結果**を得た。これはレビュー内の診断であり、事前登録済みの検証結果とは扱わない。

やらない方がよいこと: 旧配列に x 依存の係数を当てて誤差を消すこと、現在の配列を直接使った補正壁 NS を先に投入すること、推力の STEADY だけで較正データも定常と認定すること。

呼び出し側の前提への異議:

1. **効きの見積もり。** 仮に旧見積もりの二次則と係数をそのまま信じても、0.5%から1.901%への変更は **約14.46倍**で、1e−4 は約 **1.45e−3**になる。ただしこれは条件付きの算術換算で、推力改善の予測ではない。係数・摂動方向・拘束が対応していないため、現段階では改善の符号も未確定とする。Euler と NS の η 差全体を「補正で取り戻せる損失」とも扱わない。

2. **理論式と較正。** `contur_v2` の出力を半径方向等価厚さの定義に揃え、まず **δ̂_r = kδ·δ_r,CONTUR** の定数1個を L100 の分布で較正する。これは出力の較正係数であり、`cf_scale` と同一視しない。入口の運動量厚さ・発達履歴は固定して記録し、スロートを幾何的に補正しないことを「境界層厚さゼロ」と取り違えない。L70 は係数を固定した独立検証に使う。同じ壁の Euler 基準が必要で、Hall 版 L70 の Euler を黙って代用しない。

   受理線の提案は、較正前の抽出不確かさを5%以下に抑えたうえで、**L70 の出口誤差10%以下、登録した補正区間の相対 RMS 誤差10%以下、局所誤差20%以下**。相対誤差の分母には抽出不確かさによる床を置く。これらは新たな事前登録案であり、達成済みではない。失敗したら定数1個の移植性を棄却し、L70 を追加較正して「検証成功」に変えない。

3. **出口半径とスロート。** 固定 x の半径補正なら、ε_inv = (√ε_phys − δ_r,e/r_t)² を**凍結した一発補正の幾何条件**として採用する。今回の保存配列では約 **48.12**。生成器をこの ε_inv で回し直す約25分は必要な設計費とする。物理スロート半径・位置・曲率は固定し、補正の立上げでそれらを壊さないようにする。ただし、スロートの排除効果や流量差が両モードで完全に共通とは仮定せず測定する。物理 ε の固定だけでは、流れの有効 ε は固定されない。

4. **補正効果の測定。** 抽出と較正の受理後、L100 で同じ物理 L・ε・r_t、同じ NS recipe の補正壁を評価する。主判定は、事前に固定した補正後の非粘性目標に対する壁圧・境界層外縁 Mach・出口分布の誤差。壁上の NS Mach は no-slip により比較対象にならない。提案する線は、主指標の壁圧誤差を補正前から30%以上減らし、他の登録指標を不確かさ以上に悪化させないこと。η は改善必須にせず、Δη と不確かさを報告する。NS 抽出厚さを直接使う壁は診断用であり、理論式＋較正という提供方式の代替にはしない。

5. **V3 と費用。** `(D)`、`(M・補正なし)`、`(M・理論式較正補正あり)` の3列で比較し、既存の補正なしを基準列として残す。同じ L_k の比較を守り、D の Pareto 内挿値を正式判定に使わない。[plan:102](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:102)。較正 NS・Euler 基準・検証・補正形状の再生成は、現在は**キャンペーン準備費**として数える。ガス・Re・熱条件をまたぐ再利用を検証していないため、一度きりの開発費には移さない。また補正有無による実測性能差と、固定した設計の数値的不確かさ U_k は分け、性能差そのものを不確かさに加えて二重計上しない。

不足情報:

- δ_r、壁圧、出口分布の時系列と、それらの準定常 VERDICT。
- 抽出時の全設定・帯位置・感度系列。保存 NPZ は `x,d,rw` だけで、抽出の受理証跡が不足している。
- L70 の独立検証用データと、`contur_v2` のベル条件での未較正分布。

今回、推力履歴から `check_quasisteady.py` を再実行し、両 run とも **`OVERALL: ALL STEADY`** を確認した。`run_0019` の η は最終 0.9821548、単調減少の推定漸近値 **0.982141**。保存された残差 VERDICT は両 run とも **`NOT CONVERGED (stalled/plateau…)`**、判定区間は `S3_main`。残差 CSV が手元にないため残差判定の再実行はしていない。保存されたメッシュ品質 VERDICT は両者 **PASS**、NS の壁解像も **PASS**。最終場の主要8変数では NaN/Inf は0件だった。これらは δ_r の定常性を保証しない。

ファイル変更・forge 起動は行っていない。抽出再計算では診断プロセス内だけで `x_E` 欠落による流量帳簿処理を省略し、旧 NumPy の `trapz` を使用した。**plan 未反映**。反映先は `tooling-nozzle-bell-two-mode.md` §4.4・§5.1 #7・§6 V3。
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
