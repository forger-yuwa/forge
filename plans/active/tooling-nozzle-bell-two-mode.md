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
- **ゲート** (旧記述、plan レビュー M3): 受理条件は、(i) 同じ実効設定の区間 (`StageManifest`) で `check_convergence` が **PASS** (現行 `opt/driver.py:173` の「DIVERGED でなければ可」は採らない)、(ii) 両推力会計・$\eta$・$\dot m$・CV 残差の時系列が `check_quasisteady.py --series-csv` で **STEADY**、(iii) Phase 2 の物理ゲート ($\dot m/\dot m_{1D}$、$\eta$ の範囲)、(iv) CV 残差 ≤ 0.1 %、(v) メッシュ品質と局所 $y_1^+$。polish は `restart_field.py` による**同一メッシュの継続**に替える (現行 `opt/polish.py:41` は `runner.prepare` 経由で補間し直すので継続ではない)。最初の TP の seed は段階起動、形状が変わる warm start は `interp_field.py` の cross-mesh restart。
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
- **S1 の導出 (2026-10-11、diagnostician [記録](../../notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md))**: 凍結 TP でも制御面 DE の条件は古典と同じ — φ = θ + μ (DE は C⁺)、$\Phi_2=V\cos(\theta-\mu)/\cos\mu$ と $\Phi_3=r\rho V^2\sin^2\theta\tan\mu$ が DE 上で一定 — で、比熱比は陽に現れない (使う熱力学は $dp=-\rho V\,dV$・$d\rho=-\rho M^2 dV/V$ だけ)。**固定 $(L,\varepsilon)$ では縁条件 (Cain 4.3) は課さない** (端点項が消える。停留輪郭は $p_a$ に依らない。$F(E)$ は感度として報告)。D に追加条件は要らない。TP 固有の $k=d\ln a/d\ln V$ は C2 の条件と Goursat 領域の出だしにだけ入る。**構成**: kernel $\kappa(s)=1/R_t+(\kappa_B-1/R_t)S(s/s_B)$ (古典円弧も選べる) → kernel 場 (`forward_moc`) → **BD は B から下流の C⁻ を net 場の補間で追跡** (θ_B > μ_B では net の前線が逆向き) → DE を M をパラメータに $\Phi_2,\Phi_3$ から引き流量が $\dot m_{BD}(s_D)$ に達する点を E → 入れ子の 1 次元根 (内 $s_D$: $r_E=\sqrt\varepsilon$、中 $\theta_B$: $z_E=L$、外 $\kappa_B$: C2) → 領域 BDE を Goursat で埋め、壁流線を質量流量で取る。**C2 の根の関数は $\tilde g(\kappa_B)=(d\theta/d\sigma)_{DE,D}-(d\theta/d\sigma)_{kernel\,C^+,D}$** (B の曲率連続と等価、流線の 2 階微分より安定)。上の $g(\kappa_B)$ は根の後の検算に使う。複数根は MOC 推力最大。折り返し・$|c|>1$・根が括れないは INFEASIBLE と理由。
- **S1 を独立の受理ゲートにする** (plan レビュー M5): kernel の族、$\kappa_B$ の探索域、B・D の解の枝の追跡、複数根の選択規則を定義し、各 $L_k$ で根と分解能依存性を確かめる。停留性の試験は摂動幅の掃引・MOC の細分化・二次変分の符号まで見る。証明できる範囲が限られる間は「指定した kernel と拘束の下の停留輪郭」と呼び、「最大推力」と呼ばない。S1 を通るまで (M) を比較キャンペーンに入れない。
- **スロート流量の定義**: CFD ピンの $m^*$ は CFD の保存量の積分ではなく等エントロピー写像の流量 (`methods/design/overview.md` の CFD ピン節)。V1 の「スロート流量」は $m^*$ と定義し、凍結源の CFD 流量との差を別に記録する。
- **C2 接続**: 符号付きの関数 $g(\kappa_B)=\kappa_{Rao}(B^+;\kappa_B)-\kappa_B$ の根を求める ($\kappa_{Rao}(B^+)$ は kernel を $\kappa_B$ で作ったときの Rao 輪郭の出だしの曲率)。**各評価で B・D を解き直し**、$g$ の符号を挟めない候補は「C2 の成立未確認」として台帳に残し、評価に回さない。BE の MOC 壁点は、**始点の $r,r',r''$ と終点の $r,r'$ を任意に指定できる**よう `joint_fit_wall` を拡張した当てはめで B-spline にする。B を跨いで当てはめ直す方式は採らない。
- **δ\* と物理寸法** (plan レビュー M6・M7): 積分法を TP にする — 粘性係数を評価器と同じ $\mu(T)$ に、運動量厚さの式の速度勾配を TP の縁速度から得た $d\ln u_e/dx$ に (現行 `feedback/deltastar_integral.py:208` と 116 行は CPG の換算で、単一 H2O・$T_t$ 3500 K・$T$ 1000 K の合成例で −6 % ずれる — codex の検算)、回復温度をエンタルピーから。残す閉包の近似を明記する。$k_f$ は較正用の形状と検証用の形状を分けて受理する。非粘性壁を法線にずらしたあと、§4.2 の最終の物理壁の拘束を課した当てはめで確定し、拘束への射影で δ\* 補正をどれだけ戻したかを位置・壁角・曲率の残差として記録して受理値を決める。帰還ループは回さない。
- **設計変数 = $L$ の 1 個**。指定した kernel の下で長さごとに輪郭が 1 本 (S1 で一意性を確かめる) なので MOO はせず、$L_k$ で評価する。

### 4.5 比較の手順

- 両モードを §8-1 の同じ寸法 spec で回す。比較用の**共通の長さ $L_k$ を 5〜7 点**決め、(M) はその $L_k$ で、(D) は MOO (LHS 40 点 + EHVI infill 20〜30 点) の後に各 $L_k$ で $L$ を固定した $(\theta_e,\theta_{max})$ の 2 変数の再最適化 (Phase 2 の固定長最適化と同じ) を行って、同じ $L_k$ の値どうしを比べる (前線の補間で比べない)。
- TOP (Phase 2 の $\theta_n,\theta_e,L$) を 3 点だけ新しい評価器で回し、過去の結果との橋にする。
- **第一の判定量 = 前進 MOC ゲート** (2026-10-11 diagnostician の推奨を採用): §8-1 の寸法 spec・凍結 TP・初期線 (ゲートは CFD を使わないので Hall。C2 スロートの単一曲率なので Hall の前提に合う。CFD ピンとの差は後で記録) で、**壁を与えた前進 MOC** (①の `interior_vec`、analytic・converge) が出口まで衝撃波・同族特性線の交差・網の反転なしに場を構成できるか。実装と判定は④ §4.4 の基部受理ゲートと**同一モジュール** (§5.1 #3b、所有はこの系列、剥離 plan の TC1A の MOC もこれを使う)。交差の検出は `geometry/moc_diagnostics.py` の非隣接 C± 交差の検査を使う。内部衝撃波の指標はこのゲートの連続量として残す。
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
| 2c' | ゲートの組み替えを driver に反映 | `opt/driver.py` の (i) を「DIVERGED と有意な末尾上昇だけ FAIL、NOT CONVERGED は併記」、(ii) を推力系 4 列 1.5e-4・ṁ 1e-3 (cv_residual は列から外す、単調列の漸近値も 1.5e-4)、(iv) を最終 $|cv|\le1\times10^{-4}$、$\eta$ は出口面を主。本段の出力間隔を本段長の 1/20 以下に。合格: run_0003・0004・0005 の判定が §4.1 の想定どおり (0005 は PASS、0004 は η が末尾で動いて TRANSIENT、0002 は cv で落ちる)、負例試験の更新 | O | run_0003 を同一メッシュ・同一設定で +24000 継続 (`case/68.bell_two_mode/run_0005_tp_cv_continue2`、本段計 60000、出力 1000 step ごと、**変える設定は無し**)。**A**: $|\eta_{exit}(48000)-0.979183|\le1.5\times10^{-4}$ かつ継続区間の全スナップショットで $|cv|\le1\times10^{-4}$ → §4.1 の組み替えを driver の定数に反映。**B**: いずれか外れる → 壁会計の内訳 (摩擦か入口面か) を分けてから次の 1 変更 (cfl 4→6 か relax 0.7→1.0) を諮る | O (実行) / F (解釈) |
| 2d | 壁法線の近壁層 (§4.1 (v)) | #2e の merge で入る `mesh2d` の `wall_normal_layer` ($d_n$・$d_b$、2026-10-08 の冷却壁の格子用) を使う (自前で書かない)。runner の格子キーの受け渡しに `wall_normal_layer` を足す。合格: スロート y1+ ≤ 1 (`check_wall_resolution.py`)、AR ≤ 5000 かつ AR > 1000 の全セルでスキュー ≤ 0.1、入口・出口面は平面、旧格子との $\eta_{exit}$ と摩擦の差を記録 | O |
| 2e | **近軸の hoop 閉性の修正を取り込む** — **ユーザ決定 (2026-10-11)「修正は必ず取り込む」** | 修正 0b4dff4e (`mesh.axisSegmentRWeight`、既定 1。変換器が `/PLANES/rSurfVect` を書く) は `feature/nozzle-wall-fit-and-pipeline` と `feature/nozzle-core-grid` にだけあり、main には未統合。**取り込み方**: push 済みの `feature/nozzle-wall-fit-and-pipeline` を作業ツリー `forge-thruster-merge` のブランチ `integrate/thruster-m6-hoopfix` で merge (衝突なし、253 commit。M6 の作業ツリーには触らない、調整板の依存関係に記載)。同じ merge で `mesh2d` の `wall_normal_layer` (壁法線の近壁層、#2d に使う) と倍精度の格子座標も入る。**回帰** (取り込み後、記録): design の全試験、風洞の MOC 試験、ベル TP の基準ケースを新しい変換器で変換し直して同じ段階起動・継続 (run_0003/0005 と同じ長さ) を回し、$\eta_{exit}$・cv・ṁ の差を記録 (合否でなく効きの記録)。その後 `feature/thruster-nozzle-design` に merge。**回帰の結果 (2026-10-11、作業ツリー `forge-thruster-merge`)**: ビルド成功、変換器が `/PLANES/rSurfVect` を書き forge が `axisSegmentRWeight: ON` を出す。基準ケースを新しい変換器で作り直して同じ長さ (`case/68.bell_two_mode/run_0006_tp_hoopfix_staged` → `run_0007_tp_hoopfix_cont48k`、本段 12000 + 48000) — 修正前の run_0005 との差: $\eta_{exit}$ +5.1e-6、$\eta_{wall}$ −3.4e-5、ṁ +5.0e-6、摩擦 +6e-5、**CV 残差 4.25e-5 → 3.7e-6** (壁会計が −39 N 動き、2 つの会計の食い違いが 1 桁縮んだ。後半 24 点の |cv| 最大 3.5e-5)。推力系 4 列 STEADY (1.5e-4)、残差は NOT CONVERGED (stalled、併記)、NaN/Inf 0 | O |
| 2f | CPG の CV 残差 3〜5e-3 (記録) | `outlet_statPress` の亜音速分岐と節点 P の不一致で説明できる見込み。確認は run 不要: 出口列で $M_n<1$ の節点 P を Ps に置き換えて再積分し、≤ 1e-3 に落ちれば確定。TP の受理には関係しない | O |
| 3f | **runner への (D) の組み込み** (2026-10-11、V0 を閉じた後の本題 — ユーザ了承「やっていいよ」) | `evaluate/runner.py`・`geometry/wall.py` 側に `geometry.bell_mode: direct_bspline` を足し、`design_bell_c2` の壁 (C2 スロート・U→T Hermite の収縮部) で CFD 評価できるようにする (TOP・hermite は従来どおり)。INFEASIBLE は理由つきで driver に返す。前進 MOC ゲートの指標 (C⁻ 間隔比の最小・反転・交差) を評価ごとに台帳へ記録する。**③の真空の比較ではゲートを制約にせず記録だけ** (ゲートの用途は④の基部と海面作動。制約に使うのは④の MOO から、#3e の後)。合格: TP 基準ケースで (D) の 1 点が段階起動から §4.1 の組み替えたゲート (#2c') を通る | O |
| 3d | ~~(D) の成立条件の穴と区間の端~~ **完了 (2026-10-11)**: 成立条件 (ii) を x_m に結んだ符号条件に、区間の端を二分法 (1e-6 r_t) で詰めた。接触型 12 点が INFEASIBLE、残り 29 点 OK、x_m は連続 (θ_e 0.01° 刻みで \|Δx_m\| ≤ 7.95e-4)、試験 4 本 PASS。残る注意: J 最小で 29 点中 18 点が成立集合の境界 (詰めた端) で確定し、標本 2000 点の間で tol を超える場合がある (L90・25°・10° で +7.6e-7) — MOO で問題になったら端から余裕を取る。以下は元の記述 (diagnostician 2026-10-11) | `wall_bell_c2.py`: 成立条件 (ii) を「$r''\ge-tol$ on $[0,x_m]$ かつ $r''\le tol$ on $[x_m,L]$」($tol$ = `RESID_TOL`) に置き換える (x_m で r'' が 0 に接して符号を変えない壁が通っていた — `_sign_changes` が微小値を飛ばすため。V0-b の 2 点で最大壁角が θ_max を超え、細い成立区間 10 点も同型)。成立区間の端を走査点と隣の不成立点の間で二分法で詰めてから黄金分割 (x_m の 0.025 刻みの量子化を消す)。J 最小の規則は変えない。合格: 接触型の負例 (旧 V0-b の 2 点と sliver 10 点) が INFEASIBLE、既存 29 点が OK のまま、θ_e を 0.01° 振って x_m が連続 (0.025 の跳びなし)。続けて V0-c (§6) | O |
| 3e | C⁺ 族の間隔比の指標の正規化 | `moc_forward.py:854-870` の基準前線から初期線近傍を除く。合格: 形状によらない最小値 (n 321 で 0.84 @ x 0.93、n 161 で 0.885 @ x 0.20) が消え、曲げ戻しの負例で C⁺ 族の詰まりを映す。(D) の MOO でゲートを制約に使う前に | O |
| 3c | **V0 の回し直しの準備 (diagnostician 2026-10-11)** | (a) #3 の $x_m$ の内部決定、(b) #3b のガード 3 点: Hall 初期線 ($R_t\le1$) で n > 321 を拒否 (321 可・641 不可を docstring に)、前線 1 に反転がある場合は壁の不受理でなく「判定不能 (初期線不整合)」、不受理の根拠の点に軸からの添字 j を記録し j ≤ 2 の点だけなら「判定不能 (軸近傍)」。(c) V0-a の出力を `v0a_*` に退避してから回し直す。合格: 既存の試験 + 新しい試験 (内部決定の決定性・空区間の INFEASIBLE、ガード 3 点の負例) PASS | O |
| 3b | **前進 MOC ゲート (共用モジュール)** **実装済み (2026-10-11、`c636e716`)**: `geometry/moc_forward.py`、放射源流で壁 2 次・流量 2 次、負例 PASS、風洞の試験不変。未決: Hall 初期線 (R_t 1) の適合式との不整合 1.1e-2 rad で n ≥ 641 が破綻、前進の軸の θ_r は 1 次の代用 — 諮問中 ([ブリーフ](../../notes/reviews/briefs/2026-10-11-bell-v0-gate-diagnostician.md))。 | 壁を与えた前進 MOC (①の `interior_vec`、analytic・converge、対の分類を不受理の理由に、`moc_diagnostics.py` の交差検査、亜音速入力の検査) を新規モジュールに書く。③のゲート・④の基部受理・剥離 plan の TC1A の全流れ壁圧が共用し、**所有はこの系列** (別セッションには使い方だけ渡す)。合格: 風洞の放射源流の試験 (`design/tests/moc_axis_limit_radial.py`) を前進でも再現、衝撃波を要する壁と亜音速入力の負例で不受理の理由が正しく返る、`moc_kernel.py`・`moc_inverse.py` の既定はビット不変 | O |
| 4 | ~~(M) TP の Rao 条件の導出と S1 ゲート~~ **判断: 2026-10-11 diagnostician [記録](../../notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md) — 条件は γ-free で古典と同形、固定 (L, ε) で縁条件は落ち、κ_B・θ_B・s_D の入れ子の 1 次元根で閉じる。成立の見込み高。TIC への切り替えの目安: 5 本の L_k のうち 2 本以上で根が括れない / T2 が n で収束しない / T4 の c₁ が負対照の 5 % に落ちない ((ii)(iii) は 1 回修正してから)** | §4.4。凍結 TP の $h(T)$・$s(T,p)$ で拘束付き変分 (固定 $\varepsilon$ の端点条件を含む)。kernel の族・$\kappa_B$ の探索域・B と D の枝の追跡・複数根の選択規則を定義。平面縮約が `rao_planar.py` と一致するかも確かめる。S1 を通るまで (M) を比較に入れない | F |
| 5 | (M) 生成器と単体試験 | §4.4 の構成と §6 V1 の T1–T6。新規 `geometry/rao_axisym.py` (`kernel_wall`・`kernel_field`・`trace_bd`・`rao_line`・`solve_D`・`solve_thetaB`・`solve_kappaB`・`goursat_fill`・`wall_from_massflow`・`thrust_accounts`・`design_rao`)、等エントロピーの accessor (`isentrope(M)` → T, V, ρ, p, a, μ, k — `gas/semiperfect.py` は他ブランチも触るので新規モジュールの関数として足すのを優先)、試験 `design/tests/run_rao_axisym_tests.py`。Rao チャートの出典 (γ・R_t・長さの定義) は case/29 の資産で確かめる | O |
| 6 | CFD ピンと Hall 検算 | 合格: 凍結源の照合 PASS、Hall との初期線の差を記録 (合否なし) | O |
| 7 | 当てはめの拡張・δ\* の TP 化・物理壁 | §4.2・§4.4。合格: 拡張した当てはめが旧 `joint_fit_wall` の拘束を与えたときビット一致、積分法の速度勾配が TP の数値微分と相対 1e-4 以内、$k_f$ の検証形状 (較正に使わない形状) での出口 δ\* の差を記録、最終壁が §4.2 の拘束を 1e-6 以内で満たす、射影の残差 (位置・壁角・曲率) が受理値以内 (値は #4 の後に事前登録) | O |
| 8 | 比較の run と判定 | §6 V2・V3。case/68 の README に run 一覧 | O (実行) / F (解釈) |
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
  - **(diagnostician 2026-10-11 の T1–T6 で置き換え)** T1a CPG・円弧 kernel・**チャートの R_t (0.382) と γ** で θ_n・θ_e が ±1°。T1b 平面縮約が `rao_planar.theta_lip` と 1e-8。T1c 自由 ε で得た (L, ε) を固定 ε で解き直してビット一致。T1d 定比熱の擬似種で CPG 生成器と 1e-8。**T2 包絡線の恒等式** $dF/ds_D=2\pi[\Phi_3\,dL/ds_D+F(E)\,dr_E/ds_D]$ が相対 1e-5 で M 格子を倍にして収束。T3 保存: 閉合 $|W_m-E|$ と推力の 2 会計 (制御面 / 初期線 + kernel 壁 + BE 壁) の差 ≤ 1e-4 で n とともに減る (「DE の流量 = スロート流量」は構成で恒等なので試験にしない)。T4 停留性 (下の項目を具体化): 壁 BE に $\eta_k$ ($k$=1..3) を $\pm\{1,2,4\}\times10^{-3}r_t$ で足し、推力 $F(\epsilon)$ を 3 次で当てはめ $|c_1\epsilon_0|\le0.05|c_2\epsilon_0^2|$、負対照 (D を 2 % ずらした非 Rao 線) の $c_1$ の 5 % 以下、n で減る。T5 C2: $|\tilde g|\le10^{-8}$ rad/$r_t$、$\kappa_W(B^+)-\kappa_B$ は n 161/321 の変化の 3 倍以内。T6 分解能: n 161/321/641 で θ_max・θ_e・F が単調に変わる。
  - (旧) 保存: 制御面 DE を通る流量がスロートの流量と相対 1e-4 以内。制御面の推力と「入口 + 壁」の推力が相対 1e-4 以内。
  - **TP の停留性**: $c_p$ が実際に変わる条件 (基準ケースの EXH) で、拘束 ($\dot m$、$L$、$\varepsilon$) を保つ壁の摂動族を作り、MOC で求めた推力の一次変分が 0 であること。摂動振幅 $\epsilon$ を掃引して推力変化が $O(\epsilon^2)$ になること、MOC を細分しても結論が変わらないこと、二次変分の符号 (負なら局所最大、そうでなければ「停留」とだけ呼ぶ) を記録する。任意の形状に勝つことを最適性の証明として使わない。
  - スロート流量は CFD ピンの $m^*$ (等エントロピー写像の流量) と定義し、凍結源の CFD 流量との差を別に記録する。
  - C2 接続: (T5 に置き換え) 旧「$|g(\kappa_B)|\le10^{-6}$」は流線の 2 階微分では達成できない。
- **V2 (評価器の整合)**: (M) の非粘性壁を §4.1 の Euler 設定で回し、推力が MOC の値と相対 0.2 % 以内 (粗い整合確認)。
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

## 9. 変更ログ

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
