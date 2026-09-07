# ノズル排除厚さ更新計画: 積分法初期壁 + 固定 Euler 基準 CFD 反復

## メタ

- **area**: `tooling / design / boundary layer`
- **status**: `done`  <!-- 2026-09-04 起票・同日 V0–V4 完了 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「粘性 δ\* 補正」節 (本計画で書き換え)
- **related_plans**:
  - 旧: [tooling-nozzle-axismach-viscous-deltastar.md](../accepted/tooling-nozzle-axismach-viscous-deltastar.md) (A12: 相関 δ\* + ρU_x 最大縁抽出 + x_lo 継ぎはぎ — 本計画で生産経路から置換)
  - 旧: [tooling-nozzle-axismach-physical-throat.md](../accepted/tooling-nozzle-axismach-physical-throat.md) (A13: 真のスロート探索・上流 Hermite 再生成 — 機構は維持)
  - 根拠: [notes/investigations/nozzle-deltastar-throat-review.md](../../notes/investigations/nozzle-deltastar-throat-review.md) (スロート δ\* 3〜12 倍過大の実測)
  - 影響: [design-isobutane-m6-d155.md](../accepted/design-isobutane-m6-d155.md) (Md トリム 6.0144 は本計画完了後に撤回)
- **created**: `2026-09-04`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

排除厚さ補正を次の 2 段構成に置き換える。

1. **初回 NS 用の壁**は、入口から前進積分する境界層積分法で推定する (初期値生成専用)。
2. **NS 計算後**は、反復中固定した設計 Euler 場と NS 場の差から、**スロートを含む全域**の半径方向等価排除厚さ
   $\delta_r(x)$ を単一の抽出法で取り、固定点反復する。

温度縁・質量流束最大縁・$x_{lo}$・上下流の方式切替・相関とのブレンド・Md トリム・law 側 Mach 帰還は生産経路から廃止する。

到達目標 (case/45): NS/Euler 質量流量比 $1.000\pm0.003$、出口コア M $6.00\pm0.006$ を**トリムなし**で満たし、
初期壁の与え方 (相関 / 積分法) に依らず同じ固定点に収束する。

## 2. スコープ

- **やる**
  - 積分法初期推定器 `feedback/deltastar_integral.py` (断熱壁 / 指定壁温、CPG / semi-perfect)
  - 固定 Euler 基準抽出器 `metrics/deltastar.py::deltastar_from_core_matched_euler`
  - `PhysicalNozzleWall` の半径方向オフセット化、`prepare_ns` / `collect` の統合 (固定 Euler 参照・緩和・帳簿)
  - 単体試験 `design/tests/run_deltastar_tests.py`
  - V0 (既存 run で抽出器成立性) → V1 (case/45 相関初期壁から反復) → V2 (積分法初期壁) → V3 (等温) → V4 (移植性)
- **やらない (非スコープ・生産経路から廃止)**
  - 温度縁 / 質量流束最大縁の生産抽出 (旧 `deltastar_from_run` は比較用に残置)
  - $x_{lo}$、上流相関補完、積分法と CFD のブレンド
  - Md トリム、law 側 Mach 帰還
  - Sasman–Cresci 閉包そのもの (§4.1 参照。閉包差し替え口は残す)

## 3. 関連 docs と前提

- 診断ノートの事実: 与えたスロート δ\* は NS 実効値 ($0.0013$〜$0.0021\,r_t$) の 3〜12 倍。NS 質量流量が Euler 比 +0.8〜3.7 %。
  1D 感度 $d\ln M/d\ln(A/A^*)\approx 1/5$〜$1/6$ で試験部 M −0.2〜−0.7 %。観測残差と ±0.1 % で一致。
- 固定点では **NS のコア = 「設計壁 $r_{inv}$ を壁とする非粘性流」そのもの**。したがって同じ物理 $r$ で設計 Euler と一致する。
  この事実が抽出器の定義 (§4.3) の根拠。
- 設計 Euler run は全 case に存在する (case/45 `run_0001`, case/42 `run_0096`/`run_0055`, case/44 `run_0005`)。
- 構造格子 (node, index $= i\,n_j + j$) の列 $i$ は $x$ 一定線なので、断面積分は列の台形則で厳密に取れる (補間なし)。

## 4. 設計方針

### 4.1 初期推定器: CONTUR 型運動量積分法 (Sasman–Cresci の代替)

ユーザ計画は Sasman–Cresci (SC) を指定していたが、SC 閉包 (モーメント運動量方程式のせん断積分・形状係数関係) は原論文
(Sasman & Cresci, AIAA J. 1966) が手元に無く再現できない。同じ役割 (入口から前進積分、加速による薄化を表現、断熱/等温両対応)
を持ち、かつ**手元の一次資料で式が確定している** Sivells CONTUR (AEDC-TR-78-63 §5, Eq. 61–90,
[papers/nozzle_design/](../../papers/nozzle_design/)) の方法を実装する。CONTUR は超音速・極超音速風洞ノズルの粘性補正の
標準実装であり、目的に対して SC より劣る理由はない。閉包は `closure:` で差し替え可能にし、SC は原論文入手後の追加とする。

軸対称 von Kármán 運動量積分 (Eq. 61):

$$\frac{d\theta}{dx} + \theta\left[\frac{2 - M^2 + H}{M\,(1+\frac{\gamma-1}{2}M^2)}\frac{dM}{dx} + \frac{1}{r_w}\frac{dr_w}{dx}\right] = \frac{C_f}{2}\sec\phi_w$$

- プロファイル: $q/q_e = (z/\delta)^{1/N}$ (Eq. 67)、$\rho/\rho_e = T_e/T$、
  $T = T_w + a(T_{aw}-T_w)\,q/q_e + [T_e - a(T_{aw}-T_w) - T_w](q/q_e)^2$ (Eq. 69; $a=1$ Crocco / $a=0$ 放物)。
  断熱壁は $T_w = T_{aw} = T_e(1 + r_f\frac{\gamma-1}{2}M^2)$, $r_f = Pr^{1/3}$。
- $\theta,\ \delta^*$ は横曲率重み $(1 - z\cos\phi_w/r_w)$ 込み (Eq. 62/63)、$H=\delta^*/\theta$。
- 摩擦: $C_f = C_{f_i}/F_c$, $R_{\theta_i} = F_{R_\delta} R_{\theta_c}$ (Spalding–Chi 形, Eq. 70–73)、
  $F_c = [\int_0^1 (\rho/\rho_e)^{1/2} d(q/q_e)]^{-2}$、$F_{R_\delta} = \mu_e/\mu_w$、
  $C_{f_i} = 0.0773/[(\log R_{\theta_i}+4.561)(\log R_{\theta_i}-0.546)]$ (Eq. 75)。$R_{\theta_c}$ は平板形 $\theta_c$ (Eq. 74)。
- $N = N(R_\delta)$ は CONTUR Fig. 10 の曲線をテーブル化 (Eq. 77–82 の連鎖の結果)。
- 各ステーションで $\theta$ (状態量) から $\delta, N$ を内部反復で決め、$\delta^*_n$ (壁法線) を出す。
- 半径方向補正へ変換: $\delta_r = \delta^*_n/\cos\theta_w$。初回壁 $r_{phys} = r_{inv} + \delta_r$。
- 入口初期条件: YAML `theta0`/`H0` 直接指定、未指定時は仮想発達長 $x_v$ (既定 = 入口直管長) の乱流平板 $\theta_0 = 0.036\,x_v Re_{x_v}^{-0.2}$。
  初期条件と採用理由は出力 JSON へ記録し、感度 (θ0 ×0.5/×2) を診断として出す。
- 熱境界条件: `thermal_bc.mode: adiabatic | prescribed_temperature` (`Tw` 定数または `Tw_table: [[x, Tw], ...]`)。
  同じ積分器へ壁温・物性 (Sutherland、TP は `gas.T_of_M`・$\gamma(T)$) を供給する。
- **合格条件は CFD 一致率ではない**: 全域で有限・非負、壁が滑らか、スロート探索と曲率ゲートを通る、初回 NS が立つ、
  最終固定点が初期推定器に依存しない。

### 4.2 座標対応 (Euler 場を NS 断面へ写す)

固定点では NS コアは設計壁 $r_{inv}$ の非粘性流に一致するので、対応は**同じ物理 $r$** が正解である。反復途中は NS の有効壁
$r_{ref}(x) = r_{w,NS}(x) - \delta_{r,in}(x)$ (この pass に与えた補正を差し引いた壁) が設計壁に対応する。よって

$$\eta_{NS} = \frac{r}{r_{w,NS} - \delta_{r,in}},\qquad \eta_E = \frac{r}{r_{inv}}$$

とし、同じ $\eta$ で対応させる (固定点で $\eta_{NS}=\eta_E=r/r_{inv}$)。
物理壁 $r_{w,NS}$ で正規化すると固定点でも $r_{w,NS}/r_{inv}$ (M6 で 1.5〜7 %) の半径引き伸ばしが残り、
コア勾配の強い膨張部で $\delta_r$ の 10〜20 % 級の偽欠損を生むため採らない (ユーザ計画からの修正点 1)。

- 走査は $x$ 一定断面の半径方向、流束は $q=\rho U_x$。
- NS 各列の不等間隔半径格子を積分格子にし、Euler 場を同じ $x$ (列間線形)・同じ $\eta$ (列内線形) で補間する。
  $\eta_E > 1$ (Euler 壁の外側) は Euler 壁値で定数外挿。
- Euler メッシュ・解は反復中固定、HDF5 は書き換えない。$x$ の無次元化は設計スロート半径 `scale_m` で固定。

### 4.3 参照プロファイルと等価半径 (2026-09-04 V1 pass 1 で修正: 帯局所参照)

**当初案 (内側 30 % コアの単一倍率 α)** は V1 pass 1 で不採用になった。同じ $x$ でコア全体を比べると、
上流の壁 δ 誤差が特性線に沿って約 $r_w/\tan\mu \approx 45\,r_t$ 下流の軸へ運ぶ波 (コア形状差) をその $x$ の欠損に
取り込む。pass 0 の $x$=20〜50 で旧抽出より 15 % 大きな値が出て、pass 1 の壁に $+0.03\,r_t$ のこぶを作り、軸 M に
$+0.8/-1.6\,\%$ の波を立てた。pass 1 の場から同じ方式で抽出すると δ が ±25 % 動き、反復は収縮しない
(case/45 `run_0019_ns_cm_pass1`)。

**採用 (帯局所参照)**: 境界層のすぐ外側の帯 $y\in[y_b, 2y_b]$、$y_b=\max(4\,\delta_{in}(x),\ 0.02\,r_w)$ で
比 $q_{NS}/q_E$ を $y$ の 1 次で LSQ フィットし、境界層域 $y<y_b$ へ延長して $q_{ref}$ を作る (縁判定は不要、
$\delta_{in}$ は「NS 壁 − Euler 壁」で常に既知)。

$$D = 2\pi\int_{y<y_b}(q_{ref}-q_{NS})\,r\,dr,\qquad
2\pi\int_{r_{eff}}^{r_w} q_{ref}\,r\,dr = D\ \Rightarrow\ \delta_r = r_w - r_{eff}$$

$D$ が帯下端までの容量を超えれば $y_b$ を 2 倍して再試行 (2 回まで)。内側 30 % の α とコア RMS は**診断**に降格。
実測: pass 0 / pass 1 / v1 の 3 つの異なる壁の場から同じ $\delta_r(x)$ が ±1 % (x ≥ 0) で出る = 壁に依らない境界層量。
内部名は `delta_r_equiv`、関数は `metrics/deltastar.py::band_local_deficit` (旧 `core_matched_deficit` は比較用に残置)。

### 4.4 抽出ゲート

断面ごとに CSV/JSON へ: `x, delta_r_raw, delta_r_smooth, delta_r_use, delta_in, alpha, core_rms_noaxis,
delta_r_sens1/2 (帯幅 k=3/6), band_y_b, band_slope, band_rms, band_deficit_share, ok, hard_ok` (+ reason)。

- **hard** (値を使わず前回値保持): $D<0$ / 求根解なし / NaN。
- **soft** (値は使い記録): 帯フィット RMS > 0.5 % / 帯幅感度 > max(5 %, 0.001 $r_t$) / 帯内残留欠損 > 10 % of $D$。
平滑化は hard-ok の生値だけで、範囲外は端値 (外挿しない)。

**δ_r の平滑化は 5 次 P-spline** (`metrics/deltastar.py::smooth_delta_quintic`, 2026-09-04 ユーザ指摘で変更): 壁は $r_{inv}+\delta_r$ を
補間 5 次 B-spline で通すので、δ_r の 2 階微分のノイズがそのまま壁曲率に乗る。3 次平滑化 (λ=1e-3) の壁は曲率の高周波成分が
V1 pass 3 で 8e-3 [1/r_t] (設計壁の曲率 2〜6e-3 より大きい)、旧 v3 でも 1.4e-3 あった。抽出 raw を 5 次 B-spline + 係数 3 階差分ペナルティ
(ノットは $\sqrt{x}$ で等間隔・物理 2 r_t 相当、λ=1、hard 不合格は重み 0) で平滑化すると曲率ノイズ 2.2e-4 (積分法初期壁 3.9e-4 より滑らか)、
フィット残差 4〜5 % (抽出ノイズ水準)。積分法の初期 δ_r も同じ平滑化を通す。

**壁更新の単調性ガード** (`deltastar_loop.extract_and_merge`, 2026-09-04): 新しい物理壁 $r_{inv}+\delta_{next}$ が $x\ge0.5$ で
非単調なら P-spline の λ を 10 倍ずつ (最大 1000 倍) 強めて再判定し、使った λ を CSV ヘッダと summary に記録する
(run_0022 → pass 1 で x=90–94 の隣接列が帯を 1.25 倍刻みでトグルし δ が ±1 % ジグザグ → `PhysicalNozzleWall.validate` 不合格。λ=0.1 で解消)。
帯位置を x 方向に平滑化して固定する 2 パス方式 (`band_smooth_x`) も実装したが、4 場のばらつきを 3.3 → 5.4 % に悪化させたため既定 off。

**pass 0 では診断扱い、固定点で本ゲート** (ユーザ計画からの修正点 2): pass 0 は有効輪郭が 2 % 太い状態で
コア形状 RMS 1 % を広域で破り得る。破った断面も値は使い `ok=false` で記録し、合否判定は最終固定点で行う。

### 4.5 平滑化と壁更新

全域を同じ抽出法で処理 ($x_{lo}$ なし、相関補完なし、積分法混入なし、下流切替なし、端点クリップ外挿なし)。
平滑化は壁の波打ちを防ぐ最低限 (平滑化スプライン、生値との差を必ず出力)。

$$\delta_{in}^{k+1} = (1-\omega)\,\delta_{in}^{k} + \omega\,\delta_{ext}^{k},\qquad r_{phys}^{k+1} = r_{inv} + \delta_{in}^{k+1}$$

$\omega=0.5$ 初期値、安定なら 1.0。壁更新は**半径方向** (法線オフセット + $x$ シフトを廃止)。真の物理スロート探索と
上流 Hermite 再生成は維持する。

### 4.6 runner 統合と帳簿

- `prepare_ns(problem, run_dir, ic_from, euler_ref=, delta_r_csv=, omega=, initializer=)`
  - 初期壁: `deltastar_initializer` (YAML) または CLI。`dstar_blend` は削除。
  - `delta_r_csv` は全域 `delta_r_equiv.csv` を直接使用。
  - `prepare_info.json`: 初期推定器設定・熱境界条件・固定 Euler 参照・コア範囲・緩和係数・**実際に使ったスロート補正量**。
- `collect`: Euler/NS 質量流量比、質量流量由来の等価スロート補正量 $r_{t,W} - \sqrt{\dot m_{NS}/\dot m_E}$、
  抽出ゲート集計、固定点差 (max/median of $|\delta^{k+1}-\delta^k|/\delta^k$)、初期推定値と最終値の差。

## 5. 実装ステップ

1. 文書: 本 plan、`methods/design/overview.md` δ\* 節、`plans/README.md`、調査ノート §6–8。
2. `metrics/deltastar.py`: `deltastar_from_core_matched_euler` + 断面 CSV/JSON 出力 + ゲート。
3. V0: 既存 run (case/42 M5/M6・case/44・case/45) で抽出器を実行し成立性を判定 (追加計算なし)。
4. `geometry/wall_axismach.py`: `PhysicalNozzleWall(offset="radial", delta_r_x=...)`。
   `evaluate/runner_axismach.py`: `prepare_ns` 統合・`collect` 帳簿・反復ドライバ `iterate_deltastar`。
5. V1: case/45 `run_0019_ns_cm_pass0` (相関初期壁 = run_0004 の壁と同一入力, Md 6.0) → 抽出 → `run_0020_ns_cm_pass1` → 固定点。
6. `feedback/deltastar_integral.py` (CONTUR 法) + 単体試験。
7. V2: case/45 積分法初期壁から pass0/pass1、同じ固定点への収束を確認。
8. V3: 等温 (単体 + 初回 NS 開始性)。V4: case/42 M5 / case/44 へ設定不変で適用。
9. 文書同期、Md トリムの撤回、commit/push。

## 6. 検証

- **単体** (`design/tests/run_deltastar_tests.py`): $q_{NS}=c\,q_E$ で $\delta_r=0$ / 既知の壁側欠損から円環厚さ復元 /
  非一様 Euler プロファイルでも復元 / 壁半径が異なっても $\eta$ 写像で復元 / 点数・伸長率収束 / $\alpha$ で振幅差除去 /
  コア形状差をゲートで検出 / 25–35 % 感度 / 積分法 CPG 断熱・等温 / 等温 → 断熱の連続性 / 初期条件感度 / 法線→半径変換。
- **V0**: 全域で有限・非負、$\alpha(x)$ 滑らか、コア形状差小、コア範囲感度小、累積欠損が壁側で増える、下流で旧抽出器と整合、
  スロートで質量流量由来の等価値と整合。
- **V1 主ゲート**: NS/Euler 質量流量比 $1.000\pm0.003$、出口コア M $6.00\pm0.006$。品質・NaN・収束・準定常は AGENTS どおり。
- **V2**: 相関初期壁と積分法初期壁が同じ固定点 ($\delta_r$・質量流量・出口 M)。
- **V4**: case/42 M5 と case/44 にコア 30 %・平滑化・緩和係数を変えず適用。

## 7. 影響範囲

- `design/forge_design/{metrics/deltastar.py, feedback/deltastar_integral.py, geometry/wall_axismach.py, evaluate/runner_axismach.py}`
- `design/tests/run_deltastar_tests.py`
- `methods/design/overview.md`、`plans/README.md`、case/45・42・44 README run 表
- 旧 `feedback/deltastar.py` (平板相関) は比較・回帰用に残置。case 配下の `build_dstar_v3_*.py` は非推奨 (履歴として残す)。

## 8. 完了条件

- [x] 全域で同じ CFD 抽出法を使用 / 積分法は初回壁のみ / 断熱・等温を同じ実装で扱える (等温は単体試験まで; V3 の NS 実行は対象 case が無く未実施)
- [x] Euler は NS 反復中固定 / NS メッシュ変更へ補間で対応 / 半径方向補正へ統一
- [x] case/45 で質量流量と出口 M のゲートを満たす (run_0020/0021/0022/0023) / 初期推定方法によらず同じ固定点 (V1 pass3 vs V2 pass1 3 % 以内)
- [x] case/42 M5 と case/44 で移植性を確認 (run_0108 / run_0107: 設定不変で pass 0 がゲート達成)
- [x] methods・plan・case README を同期、検証済みマイルストーンごとに commit/push

## 9. 変更ログ

- `2026-09-04` — 起票 (ユーザ計画を採択。修正 3 点: η 正規化を有効壁で / pass 0 ゲートは診断扱い / SC → CONTUR 運動量積分法)。
- `2026-09-04` — 抽出器・単体試験・loop ドライバ実装、V0 成立 (スロート δ_r 0.0015 = 質量流量実効値, 全 case)。
  V1 pass 1 (`run_0019`, ω=0.5, コア全体 α): 質量流量比 1.020→1.010、x≤30 の残差 −0.4→−0.1 % に改善したが試験部軸 M に
  +0.8/−1.6 % の波 → 抽出器を帯局所参照 (§4.3) に変更。pass 0/1/v1 の場で δ_r が ±1 % 一致。pass 2 (`run_0020`, ω=1.0) 投入。
  積分法初期推定器 (CONTUR) 実装: スロート δ*_n 0.0011 (抽出 0.0015, 旧相関 0.017)、x≤20 で抽出と 5〜20 % 一致、θ0 無記憶。
- `2026-09-04` — **V1 pass 2 (`run_0020`, 帯局所抽出, ω=1.0) で主ゲート達成**: ṁ_NS/ṁ_E = 0.9999 (スロート補正 0.0014 r_t = 質量流量実効値)、
  出口面コア M 6.0008〜6.0011 (+0.01 %)、**Md 6.0 トリムなし**。軸 M は x=20–75 で Euler +0.05〜0.17 %。残: x≥80 の軸上に ±0.5〜1 % の波
  (pass 1→2 の壁変化 x≈40 [−0.02 r_t] が特性線で末端軸へ)。pass 2 の抽出は入力と x=40–60 で +4 % → pass 3 (`run_0021`, ω=0.5) で収縮確認中。
  積分法初期壁の prepare 確認: スロート +0.0011 r_t、品質 SOFT-PASS (V2 は pass 3 の後)。
- `2026-09-04` — **帯を適応化** (k=1.5 の境界層内から始め、縁候補帯の比の変化 <1 % になるまで 1.25 倍ずつ広げる; 縁候補の比がコア比より
  5 % 低ければ内部と判定して広げる): 4 つの壁の異なる NS 場での δ_r ばらつき 6.5 % → 3.3 % (旧 ρu 最大縁抽出 ≤2.8 % と同水準、スロートでも有効)。
  既知バイアス: 縁より外の尾部の取りこぼしで −3 % (べき乗則合成; 出口 M で ≤0.05 %)。フィット帯を外へずらす (fit_from 1.5) と尾部は取れるが
  コア波を拾い実場ばらつき 8 % → 不採用 (既定 fit_from=1.0)。pass 3 (`run_0021`, ω=0.5): ṁ 1.0000、出口面コア M ±0.08 %、軸波 ±1 → ±0.5 %。
- `2026-09-04` — **V2 pass 0 (`run_0022`, 積分法初期壁のみ)**: ṁ 0.9991、出口面コア M +0.04〜0.06 %、軸 M 全域 +0.1 % で平坦 (波なし)、
  抽出との差 ±2 % (x≤60)。**初期壁だけで主ゲートを満たす**。pass 1 (`run_0023`, ω=1.0) で固定点を確認中。
- `2026-09-04` — **V2 pass 1 (`run_0023`, ω=1.0)**: ṁ 0.9999、出口面コア M +0.01〜0.02 %、固定点差 p10–p90 0.997〜1.002 (収束)、
  軸 M −0.14〜+0.45 %。V1 pass 3 と δ_r が x=20–60 で 3 % 以内 → **初期推定方法に依らず同じ固定点** (§8 完了条件)。
  生産推奨 = 積分法初期壁 + 1 pass (pass 0 単独でも主ゲートは満たす)。V4 (case/42 M5 `run_0108`, case/44 `run_0107`: 積分法 pass 0) 投入。
- `2026-09-04` — **V4 (移植性, 設定不変の積分法 pass 0)**: case/42 M5 `run_0108` ṁ 0.9990・出口面コア M +0.02〜0.05 %・軸 ≤0.19 % (v3: 1.0164 / −0.3〜−0.5 % / 波 −1 %);
  case/44 `run_0107` ṁ 0.9992・0.00〜+0.03 %・軸 ≤0.36 % (v3: 1.008 / −0.2〜−0.3 %)。**status: done、accepted へ移動**。
  残件 (別 plan 候補): 等温壁の NS 実行 (V3 の CFD 部分; 対象 case なし)、抽出の −3 % 尾部バイアスの是正、
  試験部末端の軸波 (±0.5 %; 壁 δ 変化の特性線応答) の抑制 = 反復は ω≤0.5 か初期壁のみで止める運用。
  旧 A12 記述と case/45 Md トリム (6.0144) は本計画で置換: 生産形は `run_0023` (Md 6.0 トリムなし)。
- `2026-09-04` — **δ_r 平滑化を 5 次 P-spline に変更して検証** (`run_0024` V1 系 ω=0.5 / `run_0025` V2 系 ω=1.0): 壁曲率の高周波ノイズ
  7.7e-3 → 3.8e-4 [1/r_t]。ゲートは維持 (run_0025: ṁ 1.0002・出口面コア M +0.01 %・固定点 0.999〜1.006)。**試験部の軸 M の波 (x≈70 で +0.4 %) は
  滑らかな壁でも同じ** → 壁の凸凹ではなく、積分法壁と抽出壁の ~1 % の δ 差 (x≈15–25) に対する軸の集束応答。コア平均 M は平坦なので許容。
  生産形 = `run_0025` (積分法初期壁 + 抽出 1 pass + P-spline)。
- `2026-09-04` — **NS の軸 M の「一貫した +0.1 %」は軸ノードだけの数値スパイク**: run_0022/0023/0025 で r=0 ノードの P が第 1 内点より 0.5〜0.7 % 低く
  (T −0.3 %, Ux +0.03 %) M が +0.1〜0.17 %。Euler run では同じ効果が 0.14 % (M +0.02 %) と小さい (軸側の第 1 節点間隔 0.31 vs NS 0.68 r_t)。
  r/r_w ≥ 0.05 のコアは x=50/90 とも Euler 設計と ±0.05 % で一致 (旧 v3 は全コアが −0.3〜−0.4 %)。→ 「NS の軸指標を第 1 内点で読む」案は**保留** (ユーザ判断 2026-09-04: 特性曲線法の設計基準 = 軸値と不整合になる)。軸ノード値を指標のまま維持し、+0.1 % は既知の軸スパイクとして併記する。
  軸ノードの P ディップ自体は node 軸対称の軸特異性 (`nodeAxisDirichlet` は TP 不可、`axisymMethod:1` は出口角で発散 — 既知) の残件。
- `2026-09-04` — **起動レシピの A/B (ユーザ提案: soft/mid 不要? cfl 5? CFL ランプ?)**: forge に CFL ランプ機能は無いので runner に
  `run_staged_ns(stages="full"|"none"|"ramp", ramp=(1,2,3.5), ramp_steps=1000)` と `prepare_ns(cfl_main=, implicit_relax=)`、loop CLI に
  `--stages/--cfl/--implicit-relax/--ramp/--ramp-steps` を追加。run_0025 の固定点から (a) `none` + cfl 5 + relax 0.7 (`run_0026`) と
  (b) `ramp` (`run_0027`): metrics は run_0025 と同一 (出口面コア M 6.0000、ṁ 1.0002)、NaN 0・STEADY。(a) は全残差列が cfl1/30000 step の
  最終水準に 5000〜8900 step で到達し出口 M は 8000 step で凍結 → **warm start の生産レシピ = stages none, cfl 5, implicitRelax 0.7, 12000 step
  (NS ≈ 95 s、従来 195 s)**。ランプは利得なし (cold start 用に残置)。cold start (中継 → 初回 NS) を cfl 5 で立てられるかは未検証 (`full` を既定のまま)。
- `2026-09-04` — **中継の要否と Euler の CFL (ユーザ提案)**: (1) NS pass 0 を Euler 場から直接 y+~1 メッシュへ (`--init-integral --ic-from <Euler>`, 既定 full 起動)
  → 中継ありと同等 (`run_0029`: ṁ 0.9999、出口面コア M +0.00〜0.04 %) = **y+~50 の中継 run は不要**。soft/mid を抜いて cfl 5 で一気に立てると step 1 で NaN
  (`run_0030`, 削除) → cold start には soft/mid が必要。soft/mid + 本段 cfl 5 + relax 0.7 12000 step (`run_0031`) は cfl1/24000 と同等。
  (2) Euler は cfl 6 + implicitRelax 0.7 (soft のみ + 本段 12000, `run_0028`) で run_0001 と同一 metrics・ALL PASS (残差到達 4700〜8900 step)。8000 step (`run_0032`) は
  metrics 同一だが残差未達判定。**新レシピ: 設計 → Euler (soft + 12000 @cfl6/r0.7) → NS pass 0 (中継なし, full 起動, 本段 cfl5/r0.7 12000) → 必要なら pass 1 (none, cfl5/r0.7 12000)**。
  Euler runner に `--cfl/--implicit-relax/--stages/--ic-from` を追加。注: 本日午後の壁時計は別セッションの GPU 共有で 2〜3 倍に膨らんでおり要再計測。
- `2026-09-04` — **上流配管の影響と r_t 最終補正**: 直管 L_pipe 0.5 → 10 r_t (`run_0033/0034`) でスロート以降の δ_r は 0.3 % 以内・出口 M 同一 →
  収縮部で境界層履歴が消えるため配管仕様に依らない (積分法の θ0 無記憶と一致)。`solve_rt` (積分法 Newton / 実測 δ の Re^-0.2 補正) を実装し、
  1 回目 (生抽出 δ, 77.05 mm, `run_0036`) は出口半径 0.7753 m と 0.3 mm 残り → 壁に載せる平滑化後 δ を使う修正で 77.02 mm (`run_0037/0038`): 出口半径 0.7750 m、
  ṁ 1.0002、出口面コア M 6.0002。**補正後の NS 解を成果物として保持** (run_0038)。
