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
- **ゲート** (plan レビュー M3): 受理条件は、(i) 同じ実効設定の区間 (`StageManifest`) で `check_convergence` が **PASS** (現行 `opt/driver.py:173` の「DIVERGED でなければ可」は採らない)、(ii) 両推力会計・$\eta$・$\dot m$・CV 残差の時系列が `check_quasisteady.py --series-csv` で **STEADY**、(iii) Phase 2 の物理ゲート ($\dot m/\dot m_{1D}$、$\eta$ の範囲)、(iv) CV 残差 ≤ 0.1 %、(v) メッシュ品質と局所 $y_1^+$。polish は `restart_field.py` による**同一メッシュの継続**に替える (現行 `opt/polish.py:41` は `runner.prepare` 経由で補間し直すので継続ではない)。最初の TP の seed は段階起動、形状が変わる warm start は `interp_field.py` の cross-mesh restart。
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

**設計変数 = $(L,\ \theta_e,\ \theta_{max},\ x_m/L)$ の 4 個**。成立条件: $r'\ge0$、壁角の極大が 1 つ、曲率の上限。破れたら INFEASIBLE。この壁がそのまま物理壁 (δ\* 補正なし)。拘束 7 個は、拘束行列が満行ランクならスロートの 3 拘束で目的関数の零空間 (2 次多項式) が除かれるので一意に解ける。**ノットの規約・係数の数・設計変数の範囲・曲率の上限・区間内の極値の検査・無次元座標**は §5.1 #3 で固定する (plan レビュー M8)。

### 4.4 (M) 特性曲線モード (TP の Rao)

- **初期線**: スロートの C⁻ 特性線。場の出所は **CFD ピン** (同じ収縮部・$R_t$・ガスの Euler 場を 1 本回して凍結) を既定、Hall は検算。凍結源の形の照合を本計画の収縮部と $R_t$ に合わせて広げる。
- **MOC の経路**: ①風洞の生産で検証済みの配列版 `interior_vec` (`axis_limit="analytic"`・`corrector="converge"`、対ごとの分類 `stats["status"]` を不受理の理由に) に乗せ、スカラー版 `KernelMOC` は使わない (2026-10-11 ユーザ指示「風洞のスクリプトを参考に」)。マーチの組み方は風洞の逆 MOC `geometry/moc_inverse.py` を手本にし、軸対称の検証には風洞の放射源流の試験 (`design/tests/moc_axis_limit_radial.py`) を使う。風洞のコードの既定は変えない。
- **kernel**: スロートから終点 B までの膨張部。曲率がスロートの $1/R_t$ から連続に変わる曲線族で、B での曲率 $\kappa_B$ をパラメータに持つ。古典の円弧 kernel も同じ生成器で選べるようにする (回帰試験用)。
- **Rao 条件 (TP)**: 凍結 TP の $h(T)$・$s(T,p)$ を使い、質量流量と長さの拘束付きで推力の変分をとって制御面 DE の条件を導く (S1、§5.1 #4)。比熱比を陽に含まない形に書けるかは導出の結果で判断し、仮定しない。**固定 $\varepsilon$ (出口半径の固定) が端点の変分をどう制限するか**も導出で定める。与えた $(L,\varepsilon)$ に対して B と D を求め、領域 BDE を MOC で埋め、質量流量で壁流線 BE を取る。
- **S1 を独立の受理ゲートにする** (plan レビュー M5): kernel の族、$\kappa_B$ の探索域、B・D の解の枝の追跡、複数根の選択規則を定義し、各 $L_k$ で根と分解能依存性を確かめる。停留性の試験は摂動幅の掃引・MOC の細分化・二次変分の符号まで見る。証明できる範囲が限られる間は「指定した kernel と拘束の下の停留輪郭」と呼び、「最大推力」と呼ばない。S1 を通るまで (M) を比較キャンペーンに入れない。
- **スロート流量の定義**: CFD ピンの $m^*$ は CFD の保存量の積分ではなく等エントロピー写像の流量 (`methods/design/overview.md` の CFD ピン節)。V1 の「スロート流量」は $m^*$ と定義し、凍結源の CFD 流量との差を別に記録する。
- **C2 接続**: 符号付きの関数 $g(\kappa_B)=\kappa_{Rao}(B^+;\kappa_B)-\kappa_B$ の根を求める ($\kappa_{Rao}(B^+)$ は kernel を $\kappa_B$ で作ったときの Rao 輪郭の出だしの曲率)。**各評価で B・D を解き直し**、$g$ の符号を挟めない候補は「C2 の成立未確認」として台帳に残し、評価に回さない。BE の MOC 壁点は、**始点の $r,r',r''$ と終点の $r,r'$ を任意に指定できる**よう `joint_fit_wall` を拡張した当てはめで B-spline にする。B を跨いで当てはめ直す方式は採らない。
- **δ\* と物理寸法** (plan レビュー M6・M7): 積分法を TP にする — 粘性係数を評価器と同じ $\mu(T)$ に、運動量厚さの式の速度勾配を TP の縁速度から得た $d\ln u_e/dx$ に (現行 `feedback/deltastar_integral.py:208` と 116 行は CPG の換算で、単一 H2O・$T_t$ 3500 K・$T$ 1000 K の合成例で −6 % ずれる — codex の検算)、回復温度をエンタルピーから。残す閉包の近似を明記する。$k_f$ は較正用の形状と検証用の形状を分けて受理する。非粘性壁を法線にずらしたあと、§4.2 の最終の物理壁の拘束を課した当てはめで確定し、拘束への射影で δ\* 補正をどれだけ戻したかを位置・壁角・曲率の残差として記録して受理値を決める。帰還ループは回さない。
- **設計変数 = $L$ の 1 個**。指定した kernel の下で長さごとに輪郭が 1 本 (S1 で一意性を確かめる) なので MOO はせず、$L_k$ で評価する。

### 4.5 比較の手順

- 両モードを §8-1 の同じ寸法 spec で回す。比較用の**共通の長さ $L_k$ を 5〜7 点**決め、(M) はその $L_k$ で、(D) は MOO (LHS 40 点 + EHVI infill 20〜30 点) の後に各 $L_k$ で $L$ を固定した $(\theta_e,\theta_{max},x_m)$ の再最適化 (Phase 2 の固定長最適化と同じ) を行って、同じ $L_k$ の値どうしを比べる (前線の補間で比べない)。
- TOP (Phase 2 の $\theta_n,\theta_e,L$) を 3 点だけ新しい評価器で回し、過去の結果との橋にする。
- **第一の判定量 = 前進 MOC ゲート** (2026-10-11 diagnostician の推奨を採用): §8-1 の寸法 spec・凍結 TP・初期線 (ゲートは CFD を使わないので Hall。C2 スロートの単一曲率なので Hall の前提に合う。CFD ピンとの差は後で記録) で、**壁を与えた前進 MOC** (①の `interior_vec`、analytic・converge) が出口まで衝撃波・同族特性線の交差・網の反転なしに場を構成できるか。実装と判定は④ §4.4 の基部受理ゲートと**同一モジュール** (§5.1 #3b、所有はこの系列、剥離 plan の TC1A の MOC もこれを使う)。交差の検出は `geometry/moc_diagnostics.py` の非隣接 C± 交差の検査を使う。内部衝撃波の指標はこのゲートの連続量として残す。
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
| 1 | ~~§4・§6 を上位に諮る~~ | **判断: 2026-10-11 codex (diagnose) [記録](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) — NO-GO。C1・B1〜B4 を採用して改訂** (採否は §6.1)。着手順 TP 化 → D → M は採用 | F |
| 1b | ~~寸法 spec (仮値) で plan 段レビュー~~ | **判断: 2026-10-11 codex `--stage plan` [記録](../../notes/reviews/2026-10-11-tooling-nozzle-bell-two-mode-plan.md) — GO-with-changes (C0/M8/m1)、全件採用** (採否は §6.1)。寸法の仮値に桁の破綻なし ($r_e$ 0.990 m、$L_k$ 1.903〜3.172 m、必要 $C_F$ 1.624) | F |
| 2a | ~~共通評価器: TP のガスと 1D 理想値・runner の配線~~ **完了 (2026-10-11、implementer、主セッションで再確認)** | `gas/rocket_exhaust.py` (新規: CEA2 → 排気ファイル `forge_rocket_exhaust_v1`、契約 `RocketExhaust` が同じ Y・種 DB・$T_c$ から MOC の `GasSemiPerfect`・lump `EXH`・`FrozenGas` を作る、1D 理想値 `ideal_1d`)、`probdef.py` の `gas.exhaust`、`evaluate/runner.py` の TP 経路 (`gas.model: semiperfect` + `gas.exhaust` のときだけ: config・`outflow`・IC・壁平滑化・段)。**結果**: 単体試験 `design/tests/run_bell_tp_tests.py` ALL PASS (主セッションで再実行) — $\varepsilon$ 50 で 1D と CEA 凍結の差は出口温度 −5.2e-6・$\dot m_{1D}$ +1.19e-5・出口速度 −3.2e-6・出口圧力 −2.1e-5・$F_{1D,vac}$ +8.7e-6 (合格線 1e-3)、基準ケース $\dot m_{1D}$ 270.87 kg/s・$F_{1D,vac}$ 1159.1 kN・$T_c$ 3523.96 K。CPG 経路は case/40 の 17 本で config がビット同一、design の試験 31/33 PASS (落ちる 2 本 `run_sern_gates_tests`・`run_sern_moc_tests` は変更前の HEAD でも同じく落ちる既知 — sern-chain R13)。起動確認 [`case/68.bell_two_mode/run_0001_tp_soft_smoke`](../../case/68.bell_two_mode/README.md): soft 1000 step、残差・全場に NaN/Inf 0 件 (主セッションで再確認)、`NOT CONVERGED (stalled)` は soft 段なので想定内、メッシュ品質 PASS (AR 836・skew 0.417)。**主セッションの判断**: CEA の生成物を forge の熱物性・輸送データにある H2・O2・H2O・OH・H・O に限る (HO2・H2O2・O3 を除いた影響はモル分率 4.5e-5・$T_c$ −0.17 K、排気ファイルの `omitted_species` に記録) | O |
| 2b | 共通評価器: 推力の CV 会計とゲート | §4.1 の後半。触るもの: `metrics/extract.py` (TP の圧力復元、全境界の $\Pi_x$ による壁会計と出口面会計、CV 残差)、`opt/driver.py` (TP の $\dot m_{1D}$・$F_{1D}$、ゲート (i)〜(v))、`opt/polish.py` (`restart_field.py` による同一メッシュの継続)、runner の `StageManifest`。合格: 全流れ解で CV 残差 ≤ 0.1 %、ゲートの負例試験 (STALLED・DRIFTING を落とす)、**配線の回帰** (同一バイナリ・同一入力場・同一実効設定で新旧の評価コードの $\eta$ 差 ≤ 3e-4、CPG 経路)。過去 run の再現は別扱い。**2a からの引き継ぎ**: (a) TP の同一メッシュでの段の継続 (`restart_field.py` + `StageManifest`。今の TP の `warm_from` は cross-mesh 用)、(b) driver が評価ごとに YAML を書き出すと相対パスの `gas.exhaust` が切れる、(c) スロートの $y_1^+$ が約 5 (第一セル 3.5 µm) — $\le1$ には約 0.7 µm が要り、runner が `wall_first_frac_throat` を読むようにする、(d) `space.pRef` の値 (推奨レシピ §5 は軸対称で動作静圧を求めるが、ノズルは 1e7→1.4e4 Pa と幅があり値が一意でない。既存の軸対称ノズル runner も書いていない — 決める前に上位に諮る)、(e) 起動直後に収縮部 ($x/r_t$ −4.0〜−0.6) の P・T が $P_c$・$T_c$ を超える (step 250 で 1.159e7 Pa・3622 K → step 1000 で 1.068e7 Pa・3567 K と減少中。収束解で残るかを確かめる、解釈はまだしない) | O |
| 3 | (D) の壁生成 | §4.3。先にノットの規約・係数の数・設計変数の範囲・曲率の上限・極値の検査・無次元座標を本表に書く。合格: 拘束の残差 ≤ 1e-10、C2、同一入力でビット同一、成立条件の負例試験、収縮部の単調条件 | O |
| 3b | **前進 MOC ゲート (共用モジュール)** | 壁を与えた前進 MOC (①の `interior_vec`、analytic・converge、対の分類を不受理の理由に、`moc_diagnostics.py` の交差検査、亜音速入力の検査) を新規モジュールに書く。③のゲート・④の基部受理・剥離 plan の TC1A の全流れ壁圧が共用し、**所有はこの系列** (別セッションには使い方だけ渡す)。合格: 風洞の放射源流の試験 (`design/tests/moc_axis_limit_radial.py`) を前進でも再現、衝撃波を要する壁と亜音速入力の負例で不受理の理由が正しく返る、`moc_kernel.py`・`moc_inverse.py` の既定はビット不変 | O |
| 4 | (M) TP の Rao 条件の導出と S1 ゲート | §4.4。凍結 TP の $h(T)$・$s(T,p)$ で拘束付き変分 (固定 $\varepsilon$ の端点条件を含む)。kernel の族・$\kappa_B$ の探索域・B と D の枝の追跡・複数根の選択規則を定義。平面縮約が `rao_planar.py` と一致するかも確かめる。S1 を通るまで (M) を比較に入れない | F |
| 5 | (M) 生成器と単体試験 | §6 V1 | O |
| 6 | CFD ピンと Hall 検算 | 合格: 凍結源の照合 PASS、Hall との初期線の差を記録 (合否なし) | O |
| 7 | 当てはめの拡張・δ\* の TP 化・物理壁 | §4.2・§4.4。合格: 拡張した当てはめが旧 `joint_fit_wall` の拘束を与えたときビット一致、積分法の速度勾配が TP の数値微分と相対 1e-4 以内、$k_f$ の検証形状 (較正に使わない形状) での出口 δ\* の差を記録、最終壁が §4.2 の拘束を 1e-6 以内で満たす、射影の残差 (位置・壁角・曲率) が受理値以内 (値は #4 の後に事前登録) | O |
| 8 | 比較の run と判定 | §6 V2・V3。case/68 の README に run 一覧 | O (実行) / F (解釈) |
| 9 | 海面作動の比較 | 兄弟 plan の外部領域が受理されてから。§6 V4 に参照・判定量を書いてから回す | F |

## 6. 検証

**判定基準は run の前にここで固定する。** codex diagnose・plan 段と diagnostician (2026-10-11) の指摘を反映した。

- **V0 (前進 MOC ゲートの A/B、CFD なし、事前登録)**: 対象 = (D) の格子点 ($\theta_{max}$ × $\theta_e$ × $x_m/L$ の 3×3×3 を各 $L_k$、成立条件を満たす点のみ) と、Phase 2 の TOP 最適形 ($\theta_n$ 25.8°・$\theta_e$ 10°) を $\varepsilon$ 50・本計画のスロートに置き直したもの。各点でゲートの合否と内部衝撃波の指標を記録する。**A**: (D) の成立点の 80 % 以上がゲートを通る → ③④とも (D) を基部の主候補にして MOO へ進み、S1 は研究課題として後送。**B**: 通る点が 20 % 以下 → ④の基部は (M) 必須、S1 (または §8-4 の TIC) を最優先に繰り上げ、(D) は③の真空専用の候補に限る。**中間** → 結果を添えて上位に諮る。どちらでも V3 の $\eta$ 比較は「同等幅内」を既定の予想として併記し、方式選択の根拠にしない。

- **V1 (M の生成器)**:
  - 回帰: CPG・円弧 kernel で古典 Rao を再現 (Rao チャートの $\theta_n,\theta_e$ と ±1°、チャートの比熱比に合わせる)。TP 生成器を定比熱にすると CPG 生成器と 1e-8 で一致。
  - 保存: 制御面 DE を通る流量がスロートの流量と相対 1e-4 以内。制御面の推力と「入口 + 壁」の推力が相対 1e-4 以内。
  - **TP の停留性**: $c_p$ が実際に変わる条件 (基準ケースの EXH) で、拘束 ($\dot m$、$L$、$\varepsilon$) を保つ壁の摂動族を作り、MOC で求めた推力の一次変分が 0 であること。摂動振幅 $\epsilon$ を掃引して推力変化が $O(\epsilon^2)$ になること、MOC を細分しても結論が変わらないこと、二次変分の符号 (負なら局所最大、そうでなければ「停留」とだけ呼ぶ) を記録する。任意の形状に勝つことを最適性の証明として使わない。
  - スロート流量は CFD ピンの $m^*$ (等エントロピー写像の流量) と定義し、凍結源の CFD 流量との差を別に記録する。
  - C2 接続: 採用した候補で $|g(\kappa_B)|\le10^{-6}$ (曲率は $1/r_t$ 単位)。
- **V2 (評価器の整合)**: (M) の非粘性壁を §4.1 の Euler 設定で回し、推力が MOC の値と相対 0.2 % 以内 (粗い整合確認)。
- **V3 (比較、事前登録の分類)**:
  - (V3 の $\eta$ は拒否権。方式の選択は V0。)
  - 各 $L_k$ の $\Delta\eta(L_k)$ の不確かさ $U_k$ (plan レビュー M4、δ\* 項は diagnostician): **格子** = 各 $L_k$ で採用した 2 形状の組を系統格子 3 水準で評価し、$\Delta\eta$ **そのもの**の格子収束から推定 (各モードの格子差の最大値は使わない)。**反復** = 判定区間を延長した履歴と漸近値の推定から。**探索** = (D) の固定長再最適化を独立の初期点・seed で繰り返し、予算を追加したときの改善量から。探索の上限を保証できない場合、結論は「登録した予算で得た設計どうしの比較」に限定する。**δ\*** = (M) だけが背負う δ\* 補正・$k_f$ 較正・拘束への射影の誤差。$k_f$ の検証形状 (較正に使わない形状) での RANS 抽出との差から推力への寄与を見積もる。4 項は保守的に足し合わせる。
  - **$|\Delta\eta|-U_k>1.5\times10^{-3}$ → 「差あり」** (符号で優位なモードを記録)。**$|\Delta\eta|+U_k\le1.5\times10^{-3}$ → 「1.5e-3 の工学的同等幅の中で同等」**。他は判定不能。1.5e-3 は Phase 2 の Rao 照合の判定帯と同じ値を工学的な同等幅として使う。
  - 併記 (合否なし): 形状が作れない点の率、総費用、内部衝撃波の指標、$g$ と当てはめの残差。
- **V4 (海面作動)**: 兄弟 plan の外部領域の受理後に、参照・判定量 (剥離位置、海面の推力、内部衝撃波) をここに書いてから回す (未登録)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-11` | [`2026-10-11-thruster-nozzle-child-plans-diagnose.md`](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) (兄弟 plan と一括) | NO-GO, C1/M4 (本 plan 分: 共通 C1・B1〜B4) | **C1 採用** (兄弟 plan §4.4 の閉じた CV を共通定義に)。**B1 採用** ($g(\kappa_B)$ を明示、各評価で B・D を解き直す、当てはめの端点拘束を一般化 — `wall_axismach.py:1061` の始点固定を確認、$\kappa$ と $r''$ を区別)。**B2 採用** (凍結 TP の $h,s$ で導出、TP の一次変分試験、CPG は回帰)。**B3 採用** (相似倍率では $\varepsilon$ が不変なのは自明 — 拘束付き当てはめで物理壁を確定、「完成した設計手順の比較」、総費用を比較)。**B4 採用** (Euler を明示、共通 $L_k$ で比較、$U_k$ の構成、「差なし」→「同等幅内で同等」)。前提への異議 (`viscMethod: 1` は空気の Sutherland で定数ではない) を採用し §3 を訂正。0.1〜0.2 % の予測は根拠不足として撤回 |
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

## 9. 変更ログ

- `2026-10-11` — diagnostician (Fable) の横断レビュー [記録](../../notes/reviews/2026-10-11-thruster-strategy-diagnostician.md) を全件採用: 第一の判定量を前進 MOC ゲートに (η は拒否権)、順序を「(D) の壁 → ゲートの A/B (V0、CFD なし) → S1 と MOO の順を決める」に、共用モジュール #3b (所有はこの系列)、$U_k$ に δ\* 項、作動点を高度のリスト (真空は `p_a: 0`)、S1 が閉じない場合の TIC の代替をトリガ付きで §8-4 に。
- `2026-10-11` — §5.1 #2a 完了 (TP の排気・1D 理想値・runner の TP 経路。1D と CEA 凍結の差 ≤ 2.1e-5、CPG 経路ビット不変、起動確認 run_0001)。2b への引き継ぎ 5 件を #2b に追記。
- `2026-10-11` — codex plan 段 (GO-with-changes, C0/M8/m1) を全件採用して改訂: TP 化を評価経路全体 (推力抽出・理想推力・チョーク流量・IC・壁平滑化・ドライバ) に、推力を全境界の運動量の式に (出口の粘性応力)、ゲートを PASS + 推力系の STEADY に・polish を同一メッシュの継続に・真空は `outflow`、$U_k$ を $\Delta\eta$ 自体の 3 水準格子収束・反復の漸近・探索の再現性に、S1 を独立ゲートに、最終の物理壁の拘束と射影の残差、δ\* の積分法の TP 化、$L_U$・$L_{pipe}$ と単調条件。段階方式。
- `2026-10-11` — codex diagnose (NO-GO) を反映して改訂: 推力を兄弟 plan の閉じた CV に、Euler 設定の明示、`viscMethod: 2` (燃焼ガスの輸送物性)、TP の Rao を凍結 TP の $h,s$ で導出し一次変分で検査、C2 接続を符号付き $g(\kappa_B)$ の根探しに、当てはめの端点拘束の一般化、δ\* 後の物理壁を拘束付き当てはめで確定 (相似倍率は撤回)、比較を共通 $L_k$ で行い不確かさ $U_k$ を構成、「完成した設計手順の比較」と総費用。
- `2026-10-11` — 起票 (親計画 §5.1 #3・#4)。§4・§6 は上位に諮る前の下書き。
