# 特性型の遠方境界 `farfield` (node)

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) 「特性型の遠方境界 (`farfield`)」
- **related_plans**: [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §5.1 R4d (動機)
- **created**: `2026-09-27`
- **owner**: `Claude (feature/sern-design)`

## 1. 目的

外部流の計算領域を有限で打ち切る境界として、波を外へ通し外気の状態を入れる特性型の遠方境界を node に足す。
動機は SERN 3D の側方境界 (R4d): `side_far` が `slip` (質量を通さない壁) で、遠方面を z 2.50 → 3.42 H に動かすと
C_L が 7.1e-4、C_M が 0.021 動いた (許容 5e-4 / 5e-3 を超過。3.42 H より外では頭打ち)。遠方境界を入れれば、狭い領域でも
遠方面の位置に依存しない解が得られるかを確かめる。

## 2. スコープ

**2026-09-27 ユーザ決定: 「レビューどおり全部やる」** (汎用の遠方境界として受理できる水準まで検証する。SERN 用に絞る案は採らない)。

- **やる**: node・密度ベース・**SLAU/SLAU2** の境界半割面流束 (専用カーネル)、CPG と TP (単一・多成分 lump)、SST の $k,\omega$ (`sstEnergyIncludesK` 両設定)、化学種、
  定常 (陽解法・block-DPLUR) と非定常 (dual-time)。面流束の診断ダンプ。runner の `side_far` / `top_out` 切替。
- **やらない (起動時に拒否し、拒否試験を置く)**: cell 方式、ROE/HLLE/KEEP (初版は SLAU 系のみ)、凝縮・トレーサ・受動種・遷移モデル、軸対称、周期と同じ節点を共有する配置、
  動的格子。粘性の遠方境界流束は 0 (現行の非壁境界と同じ)。陰解法の境界 Jacobian は現行の ghostless 対角 A⁺ のまま (近似前処理)。

## 3. 関連 docs と前提

- 理論は [`methods/boundary.md`](../../methods/boundary.md) 「特性型の遠方境界」。参照実装は SU2 `CEulerSolver::BC_Far_Field`
  (`.external/su2-src/SU2_CFD/src/solvers/CEulerSolver.cpp:4822`): 特性で境界状態を構成し、**内部状態と構成状態を数値流束 (近似 Riemann) に渡す**。
- 既存の node 境界流束 `convectiveFlux_boundary_d` は、質量流束を R 状態から作り、流出時の運動量・エネルギーを内部状態で風上化する簡略形で、
  node 境界半割面のスカラー移流は常に内部値を使う (§4.0)。**本 BC はこれらを使わず専用の経路を持つ**。
- 化学種等の境界は `inlet_` 接頭辞、RANS は明示リスト (未知の種別は何もしない) で分岐する。

## 4. 設計方針

改訂履歴: 初稿 → plan 段 NO-GO (C1/M6/m1、全件採用) → plan-2 NO-GO (C0/M5/m2、全件採用: 本版)。

### 4.0 既存経路の事実 (コードで確認、2026-09-27)

- `convectiveFlux_boundary_d.inc.cuh:192` 以降: 質量流束は境界状態 (R)、流出時の運動量・エネルギーは内部状態で風上化。構成状態の物理流束でも近似 Riemann 流束でもない。
- node 境界半割面のスカラー移流は流入・流出によらず内部節点の値 (`scalarTransport_d.cu:166` `ext_is_self`、`passiveKernels_d.cuh` `nodeBnd`)。入口で組成が入るのは節点ピンによる。
- スカラー境界処理は対流流束より前 (`main.cpp:1458`)。SST 全エネルギーの境界 k は `inlet` 接頭辞のときだけ bvar (`convectiveFlux_d.cu:413`)。
- SST の k 式には体積ソース $(P_k-D_k)V$ がある (`ransSource_d.cu:210, 248`)。`sstEnergyIncludesK` でも保存配列 `roe` は平均流エネルギーで、全エネルギーは `roe + roK`。

### 4.1 入力

`bcondConfig.yaml`: `{kind: farfield, floats: {ro, Ux, Uy, Uz, Ps, k, omega, Y0, ...}}` = 自由流 (`inlet_uniformVelocity` と同じキー)。多成分は `Y{s}`/`X{s}` 必須、RANS は `k`, `omega` 必須。
`valueTypesOfBC["farfield"]` を新設し、自由流の値 (type 1) と構成した境界状態 (別名の bvar) を分けて持つ。

### 4.2 境界状態の構成 (SU2 と同じ判定、frozen-γ 契約)

1 境界半割面 = 1 スレッド、`ic` = 境界節点、$\hat{\mathbf n}$ = 外向き単位法線。

- **γ の契約**: 面ごとに単一の $\gamma^*$ = 内部の $\gamma_i$ (CPG は `cfg.gamma`、TP は `gamma_cell[ic]`)。音速・Riemann 不変量・エントロピー・$\rho,P$ の復元をすべて $\gamma^*$ で行う
  ($c_i=\sqrt{\gamma^*P_i/\rho_i}$、$c_\infty=\sqrt{\gamma^*P_\infty/\rho_\infty}$、$s=\rho^{\gamma^*}/P$ [SU2 の定義])。内部エネルギーだけ TP の実物性 ($T_b=P_b/(\rho_bR_{mix}(Y_b))$ → `thermo_state_at_T`)。
  誤差は V2d (独立参照解) で測る。
- **判定は自由流の法線速度 $Q_n=U_{n,\infty}$ と、自由流の実物性で計算した固定の音速 $a_\infty$ だけで行う** (面ごとに固定 = 時間で切り替わらない。plan-2 M1)。**判定用の $a_\infty$ は起動時に自由流の $T_\infty,Y_\infty$ から実物性の $\gamma_{mix}$ で一度だけ計算し、再構成用の $c_\infty$ ($\gamma^*$ で作る) とは別物** (plan-3 M1: 判定に $\gamma^*=\gamma_i$ の $c_\infty$ を使うと、内部の温度・組成の変化で判定が切り替わり、質量流束が 34 % 跳ぶ反例あり):
  - $R^+$: $Q_n > -a_\infty$ なら $U_{n,i}+2c_i/(\gamma^*-1)$、それ以外 (超音速流入) は $U_{n,\infty}+2c_\infty/(\gamma^*-1)$。
  - $R^-$: $Q_n > a_\infty$ (超音速流出) なら $U_{n,i}-2c_i/(\gamma^*-1)$、それ以外は $U_{n,\infty}-2c_\infty/(\gamma^*-1)$。
  - $U_{n,b}=(R^++R^-)/2$、$c_b=(\gamma^*-1)(R^+-R^-)/4$。
  - $Q_n>0$ (流出面): 接線速度・$s$・組成・$k,\omega$ は内部。$Q_n\le0$ (流入面・平行面): 自由流。
  - $\rho_b=(s\,c_b^2/\gamma^*)^{1/(\gamma^*-1)}$、$P_b=\rho_bc_b^2/\gamma^*$。
  - $c_b\le0$ または非有限: その面を内部状態 ($U_b=U_i$) に置換し、面 ID・理由を診断カウンタに記録 (**検証の評価区間では発動 0 を要求**)。
- **平行面で境界に達したプルーム** (SERN 側方): $Q_n=0$ なので構成状態の $s$・組成は自由流だが、面流束は下の近似 Riemann 流束で**質量流束の符号で風上化**されるので、
  流出する内部の組成・エンタルピーは内部値のまま出る。

### 4.3 面流束 (farfield 専用カーネル `farfield_flux_d`)

- 面流束 = **SLAU/SLAU2 の数値流束 $F(U_L=U_i, U_R=U_b)$** (1 次、再構成なし)。内部面の `SLAU_d` と同じ式 ($\dot m$ の圧力差項・$\chi$・$\tilde p$) を、
  L/R の状態を引数に取る `__device__` 関数として切り出して使う (内部面の `SLAU_d` は従来のまま = ビット不変。切り出し関数は同じ式の複製で、V0u で内部面の値と一致を確認)。
  TP の面エンタルピーは L/R それぞれの組成・温度で NASA (内部面と同じ `thermo_h_mix_f`)。
- 同じカーネルが `massflux[ip]` と、スカラーの面値 (質量流束の符号で風上: $\dot m>0$ なら内部、$\dot m<0$ なら $U_b$ の $Y,k,\omega$) を書く。
- スカラー移流: node 境界半割面の `ext_is_self` / `nodeBnd` 経路に、**farfield 面だけ面値配列を読む**分岐 (面フラグで判定、既存面はビット不変)。一次輸送と S3 (化学種) の両経路。
  呼び出し順は対流流束 → スカラー移流 = 同じ評価時点。
- `sstEnergyIncludesK`: $H^*=H+\tfrac53k$、$p^*=p+\tfrac23\rho k$ の $k$ を L = 内部 $k_i$、R = $U_b$ の $k_b$ として SLAU に渡す。スカラー $k$ の面値も同じ $k_b$/$k_i$ を同じ風上で使う。
- ピンしない。
- **面流束の診断ダンプ** (plan-2 M5): env `FORGE_DUMP_FARFIELD=<path>` で、farfield 面ごとに (面 ID、評価回、$\dot m$、運動量 3、エネルギー、化学種、$k,\omega$ の面流束、構成状態、置換フラグ) を書く。
  収支検証 (V2b) と SERN の境界帳簿はこのダンプを積分する (所有節点の値から流束を組み直さない)。ダンプには面積ベクトル $\mathbf S$、圧力基準 `pRef`、残差に投入した流束 (運動量は $\tilde p - p_{ref}$ 基準) をそのまま書く。SERN 帳簿 (外気圧 $p_a$ 基準) へ渡すときは運動量に $(p_{ref}-p_a)\mathbf S$ を加える (plan-3 m7)。

### 4.4 その他

- ディスパッチ・読込 (`boundaryCond.{hpp,cpp}`): `farfield` を追加、Y/X と k/omega の読込・検査を「入口または farfield」に広げる (`inletCornerWall` は入口専用のまま)。§2 の非対応構成は起動時に拒否。
- 陰解法: 境界半割面は現行の ghostless 対角 A⁺ のまま。
- runner: `evaluate.side_far_kind` / `top_out_kind` に `farfield`、`top_out_kind: outflow` が slip に落ちる不具合も直す。生成 YAML の実効 BC を照合する試験。

## 5. 実装ステップ

1. `boundaryCond.{hpp,cpp}`: `valueTypesOfBC["farfield"]`・読込・ディスパッチ・起動時の拒否。
2. 境界状態の構成を `__host__ __device__` 関数に (§4.2)。ホスト単体試験から呼ぶ。
3. SLAU の面流束を L/R 状態を取る `__device__` 関数として複製し、`farfield_flux_d` を書く (§4.3)。`convectiveFlux_d.cu` の境界ループで farfield 面を振り分ける。
4. スカラー: `scalarTransport_d.cu`・`passiveKernels_d.cuh`・`speciesTransport_d.cu` の farfield 面値分岐、`ransBoundary_d.cu` に farfield (何もしない分岐を明示)。
5. 診断ダンプ `FORGE_DUMP_FARFIELD`、収支の積分ツール `solver_density_cuda/tools/farfield_balance.py`。
6. 独立参照解: `solver_density_cuda/tools/ref1d_euler_tp.py` (1 次元 Euler 有限体積、HLLC、NASA-9 lump 物性、長い領域 = 境界の影響が評価窓に届かない)。
7. runner (`runner_sern3d.py`)・`sern_momentum.py` (farfield 面は診断ダンプを使う)。
8. docs: `methods/boundary.md` (計画中 → 実装済み、ディスパッチ表)、`procedures/recommended-settings.md`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~codex plan 段レビュー~~ **済 (2026-09-27 NO-GO C1/M6/m1、全件採用)** | | F |
| 1b | ~~plan-2~~ **済 (2026-09-27 NO-GO C0/M5/m2、全件採用で本版)** | M1 → §4.2 (判定を $Q_n$ に固定) + §4.3 (近似 Riemann 流束)、M2 → V2d、M3 → V2b、M4 → V2a、M5 → §4.3 ダンプ、m6 → V3 の位置付け、m7 → §2/§3 統一・拒否試験・置換 0 回 | F |
| 1c | ~~plan-3~~ **済 (2026-09-27 NO-GO C0/M6/m1、全件採用)** | | F |
| 1d | codex plan 段 plan-4 | 本版 | F |
| 2 | 実装 (§5 の 1–5)、V0・V0u | AWS でビルド | O |
| 3 | 独立参照解 (§5 の 6) と V1–V2 | AWS | O |
| 4 | SERN V3 | AWS。R4d へ反映 (restart は `r4d_common_restart.py` 型の index コピー、メッシュ品質、判定区間、case README の run 索引) | O |
| 5 | codex result 段 → accepted | | F |

## 6. 検証

すべて AWS。合否は実行前に固定。判定区間を明記し `check_convergence.py`・`check_quasisteady.py` の出力原本を run に残す。
**置換 (§4.2 の $c_b\le0$ 等) は各試験の評価区間で 0 回**を必要条件とする。
V1–V2 は**境界機能の受入れ**、V3 は**SERN での配置 (側方幅) の採否**。

- **V0 既存境界の不変**: farfield を含まない構成 (run_0971 設定) で、同一初期状態からの初回 `massflux` と状態ダンプが変更前バイナリとビット一致。更新後保存量は旧バイナリ 3 回反復の再現性幅以内。
- **V0u 単体試験** (ホスト): (i) 決定表の各分岐、(ii) **連続性**: $U_{n,i}$ を ±1 付近・0 付近で 1e-6 刻みに振り、面流束 5 成分の変化が入力変化に比例 (隣接点の差 ≤ 1e-4 × 流束規模) — plan-2 M1 の反例 ($U_{n,\infty}=2$、$U_{n,i}$ 0.999999/1.000001) を含む、
  (iii) 極限 (plan-3 M2: 構成状態と数値流束を分ける): 超音速流出で $U_b=U_i$、超音速流入で $U_b=U_\infty$ (相対 1e-6)。面流束はどちらも $F_{num}(U_i,U_b)$ に一致 (相対 1e-6)。物理流束との一致は $U_i=U_\infty$ (一様) のときだけ要求 ($F_{num}(U,U)=F(U)$、相対 1e-6)、(iii') TP の分類の固定: $T_i,Y_i$ を変えて $\gamma_i$ が $Q_n=a$ を横切る状況でも分岐が変わらないこと (plan-3 M1 の反例: $\gamma^*$ 1.389999/1.390001 で流束の差 ≤ 1e-4 × 流束規模)、(iv) 切り出した SLAU 関数が内部面 `SLAU_d` と同じ入力で同じ流束 (相対 1e-6)、
  (v) 非対応構成の起動拒否 (cell・ROE・凝縮・トレーサ・遷移・軸対称・周期共有) がそれぞれエラー終了。
- **V1 自由流保持** (定常、block-DPLUR): 一様流の 3D hex 直方体 (z 3 層以上)、全境界 farfield、2000 step。(a) CPG M 0.5、(b) TP 2 種 lump (SERN 外気) M 6、(c) (b) を面に 30° 傾ける、
  (d) (b) + SST。合格 (a)–(c): 全節点で $|P/P_\infty-1|,|\rho/\rho_\infty-1|,|\mathbf u-\mathbf u_\infty|/|\mathbf u_\infty|\le10^{-5}$、$|Y-Y_\infty|\le10^{-6}$。(d) は **一様状態での輸送残差の打消し**だけを見る (plan-3 M3: SST は一様流でも消滅項 $-\beta^*\rho k\omega$、$-\beta\rho\omega^2$ があり一様値は保たれない): 初回評価で帳簿ダンプの段別残差のうち対流 + $k,\omega$ 輸送の寄与が全節点で $|\cdot|\le10^{-6}\times$ (面の $\dot m k_\infty$ 規模)。ソース込みの時間発展は V2b で見る。
- **V2a 音響反射 (非定常 dual-time)**: 3D 薄板チャネル (z 3 層以上、長さ 1 m、Δx 5 mm)。
  - 本試験: 一様流 M 0.3 (CPG、空気)、ガウス圧力パルス (振幅 1e-3 P∞、半値幅 40 mm = 8 セル以上)、下流端 farfield。参照 = 同じ背景流で長さ 3 m (反射が評価窓内に戻らない)。
    評価点 = 下流端から 0.2 m 上流、評価窓 = 入射パルス通過後から 0.8 m / (c−u) までの反射到達時刻窓。**反射振幅 (参照との差の最大) ≤ 0.05 × 入射振幅**。
  - 対照 (試験が反射を検出できること): M 0 (静止) で下流端 slip → 反射 ≥ 0.9 × 入射。
  - 時間積分: `unsteady: 1`・dual-time、物理 CFL ≤ 1 (音響)。**時間精度の確認** (plan-3 M6、[`recommended-settings.md`](../../procedures/recommended-settings.md) の dual-time 節): `nSub` を倍にした run と Δt を半分にした run で、評価点の圧力時系列の差が合格許容の 1/5 以下 (= 入射振幅の 1 %)。非定常の精度判定は定常用 VERDICT と別ファイルで保存する。
  - SLAU と SLAU2 の両方。
- **V2b 保存収支** (定常、block-DPLUR): V1 の箱で内部を自由流と異なる状態 (Y 0.13、T 600 K、k・ω 10 倍) から始め、(i) 流出配置、(ii) 流入配置、`sstEnergyIncludesK` 0/1 の 4 本。
  - **離散恒等式 (毎評価)**: 全節点の対流残差の和 = −(farfield 面流束の和) (内部面は相殺)。診断ダンプ (§4.3) と帳簿ダンプ (`FORGE_DUMP_LEDGER`、全節点に印) で照合、相対 1e-6 (float32 の和の丸め幅以内)。
  - **定常後の全体収支** (plan-3 M4): 各保存量で残差の正本を $R=-\sum_{faces}F+\sum_{nodes}S$ (外向き流束は残差に $-F$、生成ソースは $+S$) とし、定常で $R\to0$ を見る。$S$ は $k,\omega$ の体積ソース (帳簿の段別残差 `res_after_sources` − `res_after_species`) で、他の保存量は 0。`sstEnergyIncludesK: 1` の全エネルギーは `res_roe` が既に全エネルギーの残差なので **`res_roK` を足さない** ([`turbulence-sst-energy-includes-k.md`](../accepted/turbulence-sst-energy-includes-k.md) §分割保持)。規格化は**事前登録した非零の代表量** (ゼロ成分を落とさない): 質量・化学種 $\rho_\infty|\mathbf u_\infty|A$、運動量 $(\rho_\infty|\mathbf u_\infty|^2+P_\infty)A$、エネルギー $\rho_\infty|\mathbf u_\infty|H_\infty A$、$k$ $\rho_\infty|\mathbf u_\infty|k_\infty A$、$\omega$ $\rho_\infty|\mathbf u_\infty|\omega_\infty A$ ($A$ = farfield 面の総面積)。**$|R|$/代表量 ≤ 1e-4** を全成分に。
  - (ii) で内部が自由流へ置き換わること ($|Y-Y_\infty|\le10^{-4}$)。
- **V2e 時間積分経路の網羅** (plan-3 M6): V1(b) と V2b(i) を **定常陽解法** でも回す (合格条件は同じ)。V2a の本試験に SST を足し、dual-time で `sstEnergyIncludesK` 0/1 の両方 (反射 ≤ 0.05、NaN・置換 0)。
- **V2c 超音速の斜め衝撃波** (定常): 3D 薄板、M 2.5、半角 10° ウェッジ、衝撃が上側境界に当たる配置。上側境界 (A) slip / (B) farfield / (C) 上方に、反射衝撃が評価線の下流端より後ろに着く高さまで広げた slip
  (共通領域の格子は同一、`z_append` 型)。評価線 = ウェッジ下流の壁面、全線。合格: $\max|p_B-p_C|\le0.02\,\Delta p_{shock}$ かつ A は 0.1 Δp 以上。
- **V2d TP 近似の独立参照** (非定常、plan-3 M5 で初期値と誤差予算を固定): 1 次元問題を独立コード `ref1d_euler_tp.py` (HLLC、NASA-9 lump) と forge の薄板 (境界 = farfield、対照に同格子で 3 倍長の forge 長領域) で解く。格子 Δx 5 mm、長さ 1 m、自由流 = SERN 外気 (T 220 K、P 2851 Pa、Y_EXH 0)。2 種類に分ける:
  - **V2d-1 接触波 (温度・組成の塊の流出/流入)**: 初期 = P・u 一様 (P∞、u = ±0.5 c∞)、中央にガウス型 (半値幅 0.1 m) の T 600 K・Y_EXH 0.13 の塊 (ρ は P・T・Y から)。連続系では圧力擾乱ゼロの移流。合格: 塊が境界を通過する時間窓で、境界から 5 セル内側の $|P-P_\infty|/P_\infty\le10^{-3}$、T・Y の時系列が参照と一致 (T は塊の振幅の 2 %、Y は 0.002 以内)。
  - **V2d-2 TP 音響**: 同じ高温・高 Y の一様状態 (T 600 K、Y 0.13) の中を進む圧力パルス (振幅 1e-3 P) が、外気 (T 220 K、Y 0) を自由流に持つ farfield に当たる (= 内外の γ が違う)。合格: 反射振幅 ≤ 0.05 × 入射 (参照との差で)。
  - 誤差予算: 独立参照は Δx を 1/4 にした解との差が合格許容の 1/5 以下であること (参照の収束確認)。forge 長領域との比較で内部の離散化誤差と境界の誤差を分ける。
  - CPG・単成分 TP・多成分 TP の 3 物性で行う。
- **V3 SERN での配置 (側方幅)**: 境界機能の受入れ (V1–V2 合格) の後。g3、同一新バイナリ。
  - V3a 同一格子の対照: 遠方面 2.50 H で `side_far` のみ slip ↔ farfield (他の BC は生成 YAML で照合)。run_0986 最終場から各 20000 step (同一格子 restart)。R4d と同じ ε・D・窓条件で差を**記録** (合否なし)。
  - V3b farfield の幅系列: 2.50 / 3.42 / 4.35 H (`z_append`、restart は `r4d_common_restart.py` 型、メッシュ品質 PASS)。
    **採否**: 隣り合う幅の D ≤ ε (C_T・C_T_with_shear・C_L 5e-4、C_M 5e-3) が成立する最小の幅を「farfield での必要幅」とする。結論は「試験した幅系列で許容内」に限定し、無限遠との一致とは言わない。
  - SERN の帳簿 (C_T 等の境界寄与を使う場合) は診断ダンプの面流束を使う。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan.md) | NO-GO, C1/M6/m1 | **全件採用** (C1・M2 はコードで再確認: `scalarTransport_d.cu:166` の `ext_is_self`、`convectiveFlux_boundary_d.inc.cuh` の流出時内部風上)。§4 を「面流束を自前で組む専用境界」に改訂、§6 に V0u/V2a–c/V3a–b。§5.1 #1 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-2.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-2.md) | NO-GO, C0/M5/m2 | 構造 (専用流束 + 同時刻スカラー面値) は妥当。**全件採用 (2026-09-27、ユーザ決定「レビューどおり全部やる」)**: M1 → §4.2 判定を自由流 $Q_n$ に固定 + §4.3 近似 Riemann (SLAU) 流束、M2 → V2d (独立 1D 参照)、M3 → V2b (離散恒等式 + SST ソース込み全体収支)、M4 → V2a (dual-time、長領域参照、M0 の slip 対照)、M5 → §4.3 面流束ダンプ、m6 → V3 を配置の採否に、m7 → §2/§3 統一・拒否試験・置換 0 回 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-3.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-3.md) | NO-GO, C0/M6/m1 | **全件採用**: M1 → §4.2 判定用 $a_\infty$ を自由流の実物性で固定・V0u(iii')、M2 → V0u(iii) 構成状態と数値流束の極限を分離、M3 → V1(d) を輸送残差の打消しに、M4 → V2b の符号・代表量・res_roK、M5 → V2d-1/2 と誤差予算、M6 → nSub 倍・Δt 半減、V2e、m7 → §4.3 ダンプの圧力基準 |

## 7. 影響範囲

- `solver_density_cuda/boundaryCond.{hpp,cpp}`、`cuda_forge/boundaryCond_d.cu`、`cuda_forge/speciesTransport_d.cu`、`cuda_forge/ransBoundary_d.cu`
- `design/forge_design/evaluate/runner_sern3d.py`、`design/forge_design/metrics/sern_momentum.py` (`OPEN_KINDS`)
- 既存ケース: 変更なし (新種別を書かなければビット不変、V0)
- docs: `methods/boundary.md`、`methods/index.md` (見出しのみ)、`procedures/recommended-settings.md`

## 8. 完了条件

- [ ] 関連 `methods/` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6 V0–V3)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`
- [ ] [`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-27` — plan-3 NO-GO (C0/M6/m1) を全件採用 (判定用音速の固定、合格条件の修正)。
- `2026-09-27` — plan-2 NO-GO (C0/M5/m2) を全件採用し全面改訂 (ユーザ決定: 検証はレビューどおり全部)。
- `2026-09-27` — codex plan 段 NO-GO (C1/M6/m1) を全件採用し §4/§6 を改訂。
- `2026-09-27` — 初稿 (ユーザ「遠方境界入れたらすっきりかもね。やってみますか」)。`methods/boundary.md` に理論節を追加し、`outflow` の説明 (実装は全量コピー) を訂正。
