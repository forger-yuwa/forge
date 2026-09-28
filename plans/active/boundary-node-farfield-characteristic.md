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

- **やる**: node・密度ベース・内部面 SLAU/SLAU2 の構成で、**境界半割面だけ HLLC** の専用カーネル、CPG と TP (単一・多成分 lump)、SST の $k,\omega$ (`sstEnergyIncludesK` 両設定)、化学種、
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

### 4.2 方式: 境界半割面の外側状態を特性で作り、境界面だけ HLLC で解く (2026-09-28、ホスト 1D 試作で決定)

**用語 (node 方式)**: node にはゴーストセル (値を持つ外側の計算点) は無い。境界節点そのものが自由度で、境界に接する半分の双対 CV を持つ。
本 BC は、境界半割面の近似 Riemann 問題に渡す**外側状態 $U_R$** を流束計算のその場で作るだけで (弱形式、SU2 の vertex 型遠方境界と同じ)、値を保持・更新する点は作らない。
スカラーもゴーストでなく面ごとの値の配列で渡す。既存の node 境界も同じく、境界節点の値と外側状態 (`bvar`) から半割面の流束を作っている。

**経緯**: 初稿〜plan-6 の机上レビューで (a) Riemann 不変量の混合は接触波で擾乱、(b) 局所線形化の特性組み立てで密度・組成の側を流向で切り替えると流束が跳ぶ、
(c) 外側に自由流をそのまま置いて SLAU で解くと超音速流出でも外側に影響され音響反射 14–18 %、と指摘された。ホスト 1D 試作
[`solver_density_cuda/tools/farfield_proto1d.py`](../../solver_density_cuda/tools/farfield_proto1d.py) (node の境界半 CV、内部 SLAU + MUSCL、RK3、CPG) で候補を比較し、
plan-7 の反例 (M1: 高温の内部を出ていく音響が外気そのままだと 24–33 % 反射 / M2: 音速の算術平均の波速で HLLC の波速順序が壊れる) も試作に入れた:

| 外側状態 + 面流束 | 超音速流出の極限 | 連続性 (1e-4 刻み) | 音響反射 一様 M 0 (Δx 5/2.5/1.25 mm) | 一様 M 0.3 | **高温内部 (T 600/220、u 0.5)** | 接触波 |
| --- | --- | --- | --- | --- | --- | --- |
| 外気そのまま + SLAU | **+25 %** | 滑らか | **61 / 70 / 78 %** | **12 / 15 / 17 %** | — | 擾乱 2e-14 |
| 外気そのまま + HLLC (Roe 平均の波速) | 0 | 滑らか | 0.50 / 0.24 / 0.08 % | 0.15 / 0.09 / 0.04 % | **24 / 30 / 33 %** | 2e-14 |
| 外気そのまま + HLLC (Davis 波速) | 0 | 滑らか | 0.51 / 0.26 / 0.10 % | 0.14 / 0.08 / 0.02 % | **31 / 38 / 43 %** | 2e-14 |
| 特性組み立て (側を流向で切替) + HLLC | 0 | **跳び 0.14** | 0.52 / 0.26 / 0.10 % | 0.14 / 0.08 / 0.02 % | — | 2e-14 |
| 特性外側状態 (線形 $w^\pm$) + HLLC (Davis 波速) [plan-7 版] | 0 (1 点のみ) | 滑らか | 0.52 / 0.26 / 0.10 % | 0.14 / 0.08 / 0.02 % | 0.15 / 0.08 / 0.03 % | 2e-14 |
| **特性外側状態 (TRRS + 超音速の滑らかな重み) + HLLC (Davis) [採用]** | **0 (掃引で差 0)** | **滑らか** | **0.50 / 0.24 / 0.08 %** | **0.15 / 0.09 / 0.04 %** | **0.15 / 0.10 / 0.04 %** | **4e-14** |

対照: M 0 で右端 slip の反射 95.7 % (試験が反射を検出できる)。極端な比 (内部 P 1e4 倍・T 10 倍) に共通速度 −3〜3 を加えた掃引でも、採用案は隣接差 ≤ 1e-6、HLL 退避 0、全て有限。
これはホスト試作の結果で、forge 本体の run ではない。

→ **採用 (plan-8 対策版、試作の `charghost2_hllcd`): 外側状態 $U_R$ を次で作り、面流束 = $\mathbf F_{\mathrm{HLLC}}(U_L=U_i,U_R)$**:

- **圧力・法線速度 (音響部)**: 内部のエントロピーのまま、内部 (外向き特性) と擬似外側 ($P_\infty,U_{n,\infty}$、内部エントロピー) の間の **2 膨張波近似 (TRRS、Toro 9.4)**:
  $z=(\gamma-1)/(2\gamma)$、$c_{po}=c_i(P_\infty/P_i)^{z}$、
  $P_R=\Big[\dfrac{c_i+c_{po}-\tfrac{\gamma-1}{2}(U_{n,\infty}-U_{n,i})}{c_iP_i^{-z}+c_{po}P_\infty^{-z}}\Big]^{1/z}$、$U_{n,R}=U_{n,i}+\dfrac{2c_i}{\gamma-1}\big(1-(P_R/P_i)^z\big)$。
  分子 $B=c_i+c_{po}-\tfrac{\gamma-1}{2}(U_{n,\infty}-U_{n,i})>0$ なら $P_R>0$ (plan-8 M1: 線形の $w^\pm$ 式は内外の速度が逆向きに大きいと負圧になった)。**$B\le10^{-3}(c_i+c_{po})$ (真空・近真空。外気と内部が法線方向に音速の数倍で離れる) はその面を $U_R=U_i$ (内部の物理流束) にして真空置換として数える** (plan-9 M1。下限へのクランプはしない。構成直後に $\rho_R,P_R$ が正・有限でない場合も同じ置換と計数。**評価区間では 0 回**)。小振幅では線形式 $w^\pm=P\pm ZU_n$ ($Z=\rho_ic_i$) に一致するので、
  外へ出る音波は内部の $Z$ のまま抜ける (高温内部の反射 0.04–0.15 %)。
- **密度 (エントロピー)・接線速度・組成・$k,\omega$**: 常に自由流側 ($\rho_R=\rho_\infty(P_R/P_\infty)^{1/\gamma}$)。流出/流入は HLLC の接触波速度 $S_*$ の符号が選ぶ (側の切替分岐なし)。
- **超音速の境目は滑らかな重みでつなぐ** (plan-8 M2/M3。硬い分岐は流束の跳びになる): $s(x)=$ smoothstep (0→1)、帯幅 $b=0.1$。**混ぜるのは原始変数 $(\rho,\ \mathbf u$ [3 成分]$,\ P,\ Y_s,\ k,\ \omega)$ を同じ重みで** (plan-9 M2: 保存量で混ぜると別の流束になる)。EOS・音速・エネルギーは混ぜた後の同じ状態から作り、**流入時に運ぶスカラーの面値も混ぜた後の $\phi_R$** (帯の外では $\phi_\infty$ に一致)。
  - 自由流の超音速流入側: $w_q=s\big(((-1+b)-Q_n/a_\infty)/b\big)$ で $U_R\leftarrow(1-w_q)U_R+w_qU_\infty$ ($Q_n/a_\infty\le-1$ で $U_R=U_\infty$)。
  - 内部の超音速流出側 (**最後に適用 = 優先**): $w_i=s\big((M_{n,i}-(1-b))/b\big)$ で $U_R\leftarrow(1-w_i)U_R+w_iU_i$ ($M_{n,i}\ge1$ で $U_R=U_i$ → $F=F(U_i)$ 厳密。内部の特性がすべて外向きなら外の情報は入らない)。
  - 帯の中 ($0.9\le M_{n,i}<1$ 等) は特性の数え方の中間で、厳密な特性境界ではない (近音速で流出する面だけに関わる)。
- **波速 (plan-7 M2)**: Davis: $S_L=\min(U_{n,L}-c_L,U_{n,R}-c_R)$、$S_R=\max(U_{n,L}+c_L,U_{n,R}+c_R)$。$S_L\le S_*\le S_R$ と星状態の密度 $>0$ を検査し、外れたら HLL に退避して数える (**評価区間では 0 回**)。
  外側状態の $\rho,P$ も構成直後に検査する (正・有限。外れたら同じく数え、その面は HLL$(U_i, U_i)$ = 内部の物理流束)。
- **TP**: $P,\rho,E,c$ は各側の実物性。TRRS と $\rho_R$ の等エントロピー関係は frozen の $\gamma$ (内部側 $\gamma_i$、自由流側 $\gamma_\infty$) で近似し、V0h/V2d で誤差を測る。

**試作での確認 (CPG、plan-9 のゲート)**: 真空の反例 ($U_{n,i}=-0.9$、$U_{n,\infty}=10$) は float64・float32 とも真空置換 1 で有限。独立掃引 14,787 点 (速度・密度比 0.1/1/10・圧力比 0.1/1/10) で float64・float32 とも非有限 0・外側状態の負値 0・HLL 退避 0、真空置換 1,017 (すべて法線方向に音速の数倍で離れる点)。混合帯の逆流例で外側状態 $(\rho,u,P,Y)=(5.5,-1.025,3.93,0.065)$、質量流束 −5.117、流入する種の面値 0.065。
**試作での確認 (CPG、plan-8 のゲート)**: 反例 ($U_{n,i}=-0.9$、$U_{n,\infty}=+0.9$) で有限・退避 0、正値性掃引 24,843 点 ($U_{n,i},U_{n,\infty}\in[-0.99,0.99]$、$P_i/P_\infty$ 0.1/1/10) で非有限 0・退避 0、
内部 $M_{n,i}=1.1$ で自由流の速度 −2〜2・圧力 0.2〜5 倍を掃引して $\max|F/F(U_i)-1|=0$、$Q_n/a_\infty=-1\pm h$ の流束差 7.5e-4 / 7.5e-6 / 7.5e-8 ($h$ 1e-4/1e-6/1e-8、$\propto h$ = 連続)、
内部 $M$ 0.8→1.2 の 1e-5 刻みの隣接差 4.6e-5、音響反射 (一様 M 0: 0.50/0.24/0.08 %、M 0.3: 0.15/0.09/0.04 %、高温内部: 0.15/0.10/0.04 %)、接触波の擾乱 4e-14。

- **既知の限界** (plan-6 M2): TP 多成分では温度・組成の違う気体の保存形の混合だけで圧力が 0.5–0.65 % ずれる (境界と無関係)。境界の試験とは分けて判定する (§6 V0p・V2d)。

### 4.3 面流束 (farfield 専用カーネル `farfield_flux_d`)

- 面流束 = $\mathbf F_{\mathrm{HLLC}}(U_L=U_i, U_R;\hat{\mathbf n})|S|$ ($U_R$ は §4.2、1 次、再構成なし)。3D では法線成分で HLLC を解き、接線速度は星状態で各側の値を保つ (Toro 10.4)。
  HLLC を L/R 状態を引数に取る **`__host__ __device__` の純粋関数**にする (GPU 配列・診断・残差加算は外側)。
- 化学種・$k,\omega$: HLLC の質量流束の符号 (= 接触波速度 $S_*$ の符号) で風上化: $\dot m>0$ なら内部値、$\dot m<0$ なら外側状態の $Y_R,k_R,\omega_R$ (§4.2 の混合後。帯の外では自由流の値)。
  (HLLC の星状態の質量流束は $\rho_*S_*$ なので、スカラー $\dot m\,\phi_{\text{風上}}$ は HLLC の多成分拡張 [星状態で $Y$ を保存] と同じ。)
- 同じカーネルが `massflux[ip]` と、スカラーの面値配列を書く。スカラー移流は node 境界半割面の `ext_is_self` / `nodeBnd` 経路に、**farfield 面だけ面値配列を読む**分岐
  (面フラグで判定、既存面はビット不変)。一次輸送と S3 (化学種) の両経路。呼び出し順は対流流束 → スカラー移流 = 同じ評価時点。
- `sstEnergyIncludesK`: エネルギー・圧力に $k$ を含める経路では、L = $k_i$、R = $k_\infty$ で $E^*=E+\rho k$、$p^*=p+\tfrac23\rho k$ として HLLC に渡す。スカラー $k$ の面値も同じ $\dot m$ の符号。
- ピンしない。非有限の流束はその面の流束を 0 にして面 ID・理由を診断カウンタへ (**検証の評価区間では 0 回**)。
- **面流束の診断ダンプ** (plan-2 M5): env `FORGE_DUMP_FARFIELD=<path>` で、farfield 面ごとに (面 ID、評価回、$\dot m$、運動量 3、エネルギー、化学種、$k,\omega$ の面流束、面積ベクトル $\mathbf S$、圧力基準 `pRef`)。
  運動量は残差に投入した $p-p_{ref}$ 基準のまま書き、SERN 帳簿 (外気圧 $p_a$ 基準) へ渡すときは $(p_{ref}-p_a)\mathbf S$ を加える (plan-3 m7)。
- **帳簿ダンプ (`FORGE_DUMP_LEDGER`) の拡張** (plan-5 m4): 化学種を実際の種数から列挙、全節点指定、一括転送、評価回・段の対応付け。

### 4.4 その他

- ディスパッチ・読込 (`boundaryCond.{hpp,cpp}`): `farfield` を追加、Y/X と k/omega の読込・検査を「入口または farfield」に広げる (`inletCornerWall` は入口専用のまま)。§2 の非対応構成は起動時に拒否。
- 陰解法: 境界半割面は現行の ghostless 対角 A⁺ のまま。
- runner: `evaluate.side_far_kind` / `top_out_kind` に `farfield`、`top_out_kind: outflow` が slip に落ちる不具合も直す。生成 YAML の実効 BC を照合する試験。

## 5. 実装ステップ

1. `boundaryCond.{hpp,cpp}`: `valueTypesOfBC["farfield"]`・読込 (Y/X、k/omega)・ディスパッチ (自由流状態の実物性計算)・起動時の拒否。
2. HLLC の `__host__ __device__` 純粋関数 (CPG/TP、3D 法線、$k$ 込みの $E^*,p^*$ 分岐)。ホスト単体試験から呼ぶ。
3. `farfield_flux_d` (§4.3) と、`convectiveFlux_d.cu` の境界ループでの振り分け (farfield 面は既存の `convectiveFlux_boundary_d` を通らない)。
4. スカラー: `scalarTransport_d.cu`・`passiveKernels_d.cuh`・`speciesTransport_d.cu` の farfield 面値分岐、`ransBoundary_d.cu` に farfield (何もしない分岐を明示)。
5. 診断ダンプ `FORGE_DUMP_FARFIELD`、帳簿ダンプの拡張、収支の積分ツール `solver_density_cuda/tools/farfield_balance.py`。
6. 独立参照解 `solver_density_cuda/tools/ref1d_euler_tp.py` (1D、HLLC、NASA-9 lump、長い領域)。1D 試作 `farfield_proto1d.py` を TP に拡張して流用してよい。
7. runner (`runner_sern3d.py`)・`sern_momentum.py` (farfield 面は診断ダンプを使う)。
8. docs: `methods/boundary.md` (計画中 → 実装済み、ディスパッチ表)、`procedures/recommended-settings.md`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~codex plan 段レビュー~~ **済 (2026-09-27 NO-GO C1/M6/m1、全件採用)** | | F |
| 1b | ~~plan-2~~ **済 (2026-09-27 NO-GO C0/M5/m2、全件採用で本版)** | M1 → §4.2 (判定を $Q_n$ に固定) + §4.3 (近似 Riemann 流束)、M2 → V2d、M3 → V2b、M4 → V2a、M5 → §4.3 ダンプ、m6 → V3 の位置付け、m7 → §2/§3 統一・拒否試験・置換 0 回 | F |
| 1c | ~~plan-3~~ **済 (2026-09-27 NO-GO C0/M6/m1、全件採用)** | | F |
| 1d | ~~plan-4~~ **済 (2026-09-27 NO-GO C0/M3/m2、全件採用)** | | F |
| 1e | ~~plan-5~~ **済 (2026-09-27 NO-GO C0/M2/m3、全件採用)** | | F |
| 1f | ~~plan-6~~ **済 (2026-09-27 NO-GO C0/M3/m1、全件採用)** | | F |
| 1g | ~~ホスト 1D 試作で境界流束の候補を比較~~ **済 (2026-09-28)**。経過: 自由流そのまま + HLLC (plan-7 で高温内部の反射 24–33 % と判明) → 線形の特性外側状態 + HLLC (plan-8 で負圧・超音速流出の独立性・超音速流入の跳び) → **TRRS + 滑らかな重み + HLLC (現採用、§4.2)**。未通過: TP 版 (V0h) (plan-6 M1) | `solver_density_cuda/tools/farfield_proto1d.py`: node 境界半 CV を持つ 1D Euler (CPG、内部 2 次 MUSCL・境界 1 次、RK3)。候補 = (a) 自由流ゴースト + SLAU、(b) 自由流ゴースト + HLLC、(c) 局所線形化の特性振幅で組み立て + HLLC、(d) 同 + SLAU。測る量 = (1) 音響反射率 (M 0 / 0.3、長領域との差、Δx 3 水準)、(2) 超音速流出の極限 (外側状態を変えても流束が内部の物理流束に一致)、(3) 接触波 (P・u 同じで T 違い) の擾乱、(4) 流向反転・音速通過の連続性。合格 = (1) ≤ 5 %、(2) 相対 1e-6、(3) 圧力擾乱 ≤ 1e-4 P∞、(4) 流束の跳びなし。結果で §4 を決め直し、plan-7 に回す | F/O |
| 1h | ~~plan-7~~ **済 (2026-09-28 NO-GO C0/M2/m1、全件採用、試作で対策を確認)** | | F |
| 1i | ~~plan-8~~ **済 (2026-09-28 NO-GO C0/M3/m1、全件採用、試作で確認)** | | F |
| 1k | ~~plan-9~~ **済 (2026-09-28 NO-GO C0/M2/m1、全件採用、試作で確認)** | | F |
| 1l | codex plan 段 plan-10 | 本版 | F |
| 1j | V0h ホストゲート (TP 拡張) | 本体実装の前 | O |
| 2 | 実装 (§5 の 1–5)、V0・V0u | AWS でビルド | O |
| 3 | 独立参照解 (§5 の 6) と V1–V2 | AWS | O |
| 4 | SERN V3 | AWS。R4d へ反映 (restart は `r4d_common_restart.py` 型の index コピー、メッシュ品質、判定区間、case README の run 索引) | O |
| 5 | codex result 段 → accepted | | F |

## 6. 検証

すべて AWS。合否は実行前に固定。判定区間を明記し `check_convergence.py`・`check_quasisteady.py` の出力原本を run に残す。
**非有限流束の置換 (§4.3) は各試験の評価区間で 0 回**を必要条件とする。
V1–V2 は**境界機能の受入れ**、V3 は**SERN での配置 (側方幅) の採否**。

- **V0 既存境界の不変**: farfield を含まない構成 (run_0971 設定) で、同一初期状態からの初回 `massflux` と状態ダンプが変更前バイナリとビット一致。更新後保存量は旧バイナリ 3 回反復の再現性幅以内。
- **V0h ホストゲート (本体実装の前、plan-7 推奨)**: `farfield_proto1d.py` を TP (単成分・多成分 lump、既存 `FrozenGas`) に拡張し、§4.2 の表の全項目 + V2d-2 (TP の高温内部を出ていく音響、反射 ≤ 5 %) + 波速順序・星状態の妥当性・退避 0 を確認してから CUDA に進む。
- **V0p TP の保存形混合 (前提、境界と無関係)** (plan-6 M2): 単一 CV に高温側 (T 600 K、Y 0.13) と外気 (T 220 K、Y 0) の保存量を 0.9:0.1 で混ぜ、EOS から復元した圧力の誤差を記録する (期待 +0.65 %)。これは境界の合否に入れず、V2d の「長領域自身の誤差」の説明に使う。
- **V0u 単体試験** (ホスト、HLLC の `__host__ __device__` 関数): (i) **連続性**: $U_{n,i}$・$\rho_i$・$P_i$・接線速度を、法線 Mach ±1・0 付近で 1e-6 刻みに掃引し、面流束 5 成分・スカラー面流束の隣接点の差 ≤ 1e-4 × 流束規模 (plan-2〜5 の反例入力をすべて含む)、(ii) **接触波**: CPG と TP で P・$U_n$ 同じ・T 600/220 K・Y 0.13/0 の内部/外気、流出と流入の両方で面流束 = 風上側の物理流束 (相対 1e-6)、(iii) **流向の全組合せ**: $Q_n$ の符号 × 実際の $\dot m$ の符号 × ゼロ通過で、$\dot m<0$ のスカラー面値が外側状態の値 (帯の外では外気値)・$\dot m>0$ で内部値、(iii'') **真空と float32** (plan-9 M1): 真空の反例と、速度 ($U_{n,i}$ −3〜3、$U_{n,\infty}$ −3〜10)・密度比・圧力比の独立掃引を float64 と float32 で行い、非有限 0・外側状態の負値 0 を確認、真空置換と HLL 退避の回数を別々に記録 (有限であることと置換が無いことを分けて判定)、(iii''') **混合帯の逆流**: 内部 $(1,0.95,1/\gamma,Y\,0.13)$・自由流 $(10,-3,10/\gamma,Y\,0)$ で、流入の種の面値が外側状態の $Y_R$ (0.065) と一致すること、TP・SST (k 込み) でも同じ、(iv) 一様 ($U_i=U_\infty$) で $F_{num}=F(U)$ (相対 1e-6)、(v) 非対応構成の起動拒否 (cell・ROE・凝縮・トレーサ・遷移・軸対称・周期共有) がそれぞれエラー終了。
- **V0k CUDA 一致**: 同じ入力を GPU カーネルに与え、ホストの HLLC 関数と面流束が相対 1e-6 で一致。
- **V1 自由流保持** (定常、block-DPLUR): 一様流の 3D hex 直方体 (z 3 層以上)、全境界 farfield、2000 step。(a) CPG M 0.5、(b) TP 2 種 lump (SERN 外気) M 6、(c) (b) を面に 30° 傾ける、
  (d) (b) + SST。合格 (a)–(c): 全節点で $|P/P_\infty-1|,|\rho/\rho_\infty-1|,|\mathbf u-\mathbf u_\infty|/|\mathbf u_\infty|\le10^{-5}$、$|Y-Y_\infty|\le10^{-6}$。(d) は **一様状態での輸送残差の打消し**だけを見る (plan-3 M3: SST は一様流でも消滅項 $-\beta^*\rho k\omega$、$-\beta\rho\omega^2$ があり一様値は保たれない): 初回評価で帳簿ダンプの段別残差のうち対流 + $k,\omega$ 輸送の寄与が全節点で、保存量ごとの規模 (その節点に接する面流束の絶対和) の 1e-5 以下 ($\rho,\rho\mathbf u,\rho E$ は各流束、$k$ は $\dot m k_\infty$、$\omega$ は $\dot m\omega_\infty$。plan-4 m5)。ソース込みの時間発展は V2b で見る。
- **V2a 音響反射 (非定常 dual-time)**: 3D 薄板チャネル (z 3 層以上、長さ 1 m、Δx 5 mm)。
  - 本試験: 一様流 M 0.3 (CPG、空気)、ガウス圧力パルス (振幅 1e-3 P∞、半値幅 40 mm = 8 セル以上)、下流端 farfield。参照 = 同じ背景流で長さ 3 m (反射が評価窓内に戻らない)。
    評価点 = 下流端から 0.2 m 上流、評価窓 = 入射パルス通過後から 0.8 m / (c−u) までの反射到達時刻窓。**反射振幅 (参照との差の最大) ≤ 0.05 × 入射振幅**。
  - 対照 (試験が反射を検出できること): M 0 (静止) で下流端 slip → 反射 ≥ 0.9 × 入射。
  - 時間積分: `unsteady: 1`・dual-time、物理 CFL ≤ 1 (音響)。**時間精度の確認** (plan-3 M6、[`recommended-settings.md`](../../procedures/recommended-settings.md) の dual-time 節): `nSub` を倍にした run と Δt を半分にした run で、評価点の圧力時系列の差が合格許容の 1/5 以下 (= 入射振幅の 1 %)。非定常の精度判定は定常用 VERDICT と別ファイルで保存する。
  - SLAU と SLAU2 の両方。
- **V2b 保存収支** (定常、block-DPLUR): V1 の箱で内部を自由流と異なる状態 (Y 0.13、T 600 K、k・ω 10 倍) から始め、(i) 流出配置、(ii) 流入配置、`sstEnergyIncludesK` 0/1 の 4 本。
  - **離散恒等式 (毎評価)**: 全節点の対流残差の和 = −(farfield 面流束の和) (内部面は相殺)。診断ダンプ (§4.3) と帳簿ダンプ (`FORGE_DUMP_LEDGER`、全節点に印) で照合。**判定 (plan-4 m5)**: |差| ≤ 1e-5 × (その保存量の全面流束の絶対和) (正味和で割らない。float32 の atomicAdd 積算の丸め幅を含む。ホスト集計は double)。
  - **定常後の全体収支** (plan-3 M4): 各保存量で残差の正本を $R=-\sum_{faces}F+\sum_{nodes}S$ (外向き流束は残差に $-F$、生成ソースは $+S$) とし、定常で $R\to0$ を見る。$S$ は $k,\omega$ の体積ソース (帳簿の段別残差 `res_after_sources` − `res_after_species`) で、他の保存量は 0。`sstEnergyIncludesK: 1` の全エネルギーは `res_roe` が既に全エネルギーの残差なので **`res_roK` を足さない** ([`turbulence-sst-energy-includes-k.md`](../accepted/turbulence-sst-energy-includes-k.md) §分割保持)。規格化は**事前登録した非零の代表量** (ゼロ成分を落とさない): 質量・化学種 $\rho_\infty|\mathbf u_\infty|A$、運動量 $(\rho_\infty|\mathbf u_\infty|^2+P_\infty)A$、エネルギー $\rho_\infty|\mathbf u_\infty|H_\infty A$、$k$ $\rho_\infty|\mathbf u_\infty|k_\infty A$、$\omega$ $\rho_\infty|\mathbf u_\infty|\omega_\infty A$ ($A$ = farfield 面の総面積)。**$|R|$/代表量 ≤ 1e-4** を全成分に。
  - (ii) で内部が自由流へ置き換わること ($|Y-Y_\infty|\le10^{-4}$)。
- **V2f 自由流分類と局所流況が違う面** (plan-5 M2): V1 の箱で、自由流は面から超音速流出 ($Q_n=2a_\infty$) だが、内部の初期状態を局所に逆流 ($U_{n,i}=-0.2a$) させた帯を置く。合格 (plan-6 M3): **実際の $\dot m<0$ の面**で化学種面値が外気値・$\dot m>0$ の面で内部値 (内部速度の符号ではなく流束の符号で判定)、NaN・置換 0、定常後は自由流へ戻る。逆流帯の強さは、HLLC で実際に $\dot m<0$ になる面が生じることをホスト関数で事前に確かめて決める。
- **V2e 時間積分経路の網羅** (plan-3 M6): V1(b) と V2b(i) を **定常陽解法** でも回す (合格条件は同じ)。V2a の本試験に SST を足し、dual-time で `sstEnergyIncludesK` 0/1 の両方 (反射 ≤ 0.05、NaN・置換 0)。
- **V2c 超音速の斜め衝撃波** (定常): 3D 薄板、M 2.5、半角 10° ウェッジ、衝撃が上側境界に当たる配置。上側境界 (A) slip / (B) farfield / (C) 上方に、反射衝撃が評価線の下流端より後ろに着く高さまで広げた slip
  (共通領域の格子は同一、`z_append` 型)。評価線 = ウェッジ下流の壁面、全線。合格: $\max|p_B-p_C|\le0.02\,\Delta p_{shock}$ かつ A は 0.1 Δp 以上。
- **V2d TP での境界の独立参照** (非定常、plan-3 M5 で初期値と誤差予算を固定): 1 次元問題を独立コード `ref1d_euler_tp.py` (HLLC、NASA-9 lump) と forge の薄板 (境界 = farfield、対照に同格子で 3 倍長の forge 長領域) で解く。格子 Δx 5 mm、長さ 1 m、自由流 = SERN 外気 (T 220 K、P 2851 Pa、Y_EXH 0)。2 種類に分ける:
  - **V2d-1 接触波 (温度・組成の塊の流出/流入)**: 初期 = P・u 一様 (P∞、u = ±0.5 c∞)、中央にガウス型 (半値幅 0.1 m) の T 600 K・Y_EXH 0.13 の塊 (ρ は P・T・Y から)。連続系では圧力擾乱ゼロの移流。合格: 塊が境界を通過する時間窓で、境界から 5 セル内側の $|P-P_\infty|/P_\infty\le10^{-3}$、T・Y の時系列が参照と一致 (T は塊の振幅の 2 %、Y は 0.002 以内)。
  - **V2d-2 TP 音響**: 同じ高温・高 Y の一様状態 (T 600 K、Y 0.13) の中を進む圧力パルス (振幅 1e-3 P) が、外気 (T 220 K、Y 0) を自由流に持つ farfield に当たる (= 内外の γ が違う)。**パルスなしの対照を必須** (plan-4 M1: 内外の状態差だけで境界が擾乱を出さないこと。評価点の $|P-P_0|$ ≤ 0.01 × パルス振幅) を先に通し、その上で反射振幅 ≤ 0.05 × 入射 (参照との差で)。
  - 誤差予算: 独立参照は Δx を 1/4 にした解との差が合格許容の 1/5 以下であること (参照の収束確認)。**判定を 2 つに分ける** (plan-6 M2): (a) forge 長領域 vs 独立参照 = 既存の内部離散化 (TP 保存形混合を含む) の誤差 — 記録のみ、合格条件に入れない。(b) forge 短領域 vs forge 長領域 = **境界が加える誤差** — 上の合格条件はこちらに適用する。
  - CPG・単成分 TP・多成分 TP の 3 物性で行う。
- **V3 SERN での配置 (側方幅)**: 境界機能の受入れ (V1–V2 合格) の後。g3、同一新バイナリ。
  - V3a 同一格子の対照: 遠方面 2.50 H で `side_far` のみ slip ↔ farfield (他の BC は生成 YAML で照合)。run_0986 最終場から各 20000 step (同一格子 restart)。R4d と同じ ε・D・窓条件で差を**記録** (合否なし)。
  - V3b farfield の幅系列: 2.50 / 3.42 / 4.35 H (`z_append`、restart は `r4d_common_restart.py` 型、メッシュ品質 PASS)。
    **採否** (plan-4 M3・plan-5 m5): 候補幅を**それより広いすべての試験済み幅** (少なくとも 1 つ必須) と比較し、全 4 量で D ≤ ε (C_T・C_T_with_shear・C_L 5e-4、C_M 5e-3) を満たす最小の幅を「farfield での必要幅」とする。最大幅しか残らない (比較相手が無い) ときは「必要幅未確定」とし、系列を広げる。結論は「試験した幅系列で許容内」に限定し、無限遠との一致とは言わない。
  - SERN の帳簿 (C_T 等の境界寄与を使う場合) は診断ダンプの面流束を使う。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan.md) | NO-GO, C1/M6/m1 | **全件採用** (C1・M2 はコードで再確認: `scalarTransport_d.cu:166` の `ext_is_self`、`convectiveFlux_boundary_d.inc.cuh` の流出時内部風上)。§4 を「面流束を自前で組む専用境界」に改訂、§6 に V0u/V2a–c/V3a–b。§5.1 #1 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-2.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-2.md) | NO-GO, C0/M5/m2 | 構造 (専用流束 + 同時刻スカラー面値) は妥当。**全件採用 (2026-09-27、ユーザ決定「レビューどおり全部やる」)**: M1 → §4.2 判定を自由流 $Q_n$ に固定 + §4.3 近似 Riemann (SLAU) 流束、M2 → V2d (独立 1D 参照)、M3 → V2b (離散恒等式 + SST ソース込み全体収支)、M4 → V2a (dual-time、長領域参照、M0 の slip 対照)、M5 → §4.3 面流束ダンプ、m6 → V3 を配置の採否に、m7 → §2/§3 統一・拒否試験・置換 0 回 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-3.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-3.md) | NO-GO, C0/M6/m1 | **全件採用**: M1 → §4.2 判定用 $a_\infty$ を自由流の実物性で固定・V0u(iii')、M2 → V0u(iii) 構成状態と数値流束の極限を分離、M3 → V1(d) を輸送残差の打消しに、M4 → V2b の符号・代表量・res_roK、M5 → V2d-1/2 と誤差予算、M6 → nSub 倍・Δt 半減、V2e、m7 → §4.3 ダンプの圧力基準 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-4.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-4.md) | NO-GO, C0/M3/m2 | **全件採用**: M1 → §4.2 を局所線形化の特性振幅 (Whitfield–Janus/Blazek) に変更・V0u(vi)・V2d-2 のパルスなし対照、M2 → エントロピー部を実際の流向 ($U_{n,b}$、食い違い時は $\dot m$) で選択・V0u(vii)、M3 → V3b の採否を全広幅と比較、m4 → 明示 $Y,T$ の熱物性関数・V0u(viii)、m5 → 規格化の保存量別・絶対和 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-5.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-5.md) | NO-GO, C0/M2/m3 | **全件採用**: M1 (密度側の選び直しが SLAU の $\dot m$ と自己整合しない)・M2 (固定の超音速分類が局所逆流の情報を失う) → §4.2 を「外側に自由流をそのまま置き SLAU で解く」に変更 (組み立て・分類を撤去)、V0u(i)(iii)・V2f。m3 → 流束関数を `__host__ __device__` に・V0k、m4 → 帳簿ダンプ拡張、m5 → V3b の最大幅規則 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-6.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-6.md) | NO-GO, C0/M3/m1 | **全件採用**。M1 (自由流ゴースト + SLAU は超音速流出でも外側速度が流束を 25 % 変え、1D 線形化で音響反射 14–18 % [格子細分で下がらない]) → 境界だけ HLLC を候補にし、**全体実装の前にホストの 1D 試作で候補を比較** (§5.1 #1g)。M2 (TP 多成分の保存形混合で圧力 +0.5〜0.65 %、境界と無関係の既存性質) → V0u に単一 CV 更新 + EOS 復元の前提試験、V2d で長領域自身の誤差と短領域の追加誤差を分離。M3 → V2f の判定対象を「実際の ṁ<0 の面」に。m4 → §5・§6・methods・README の同期 |
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-7.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-7.md) | NO-GO, C0/M2/m1 | **全件採用**。M1 (外気そのままの外側状態は、高温の内部を出ていく音響を 24–33 % 反射) → §4.2 外側状態の音響部を境界節点の $Z$ で線形化した特性量に (試作で 0.03–0.15 %)、V0h ホストゲートに TP 版を追加。M2 (音速の算術平均の波速で $S_*>S_R$・流束の跳び) → Davis 波速 + 順序・星密度の検査 + HLL 退避と計数 (試作の極端な比の掃引で退避 0・跳びなし)。m3 → §7・methods の同期 |
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-8.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-8.md) | NO-GO, C0/M3/m1 | **全件採用、試作で対策を確認**: M1 (線形 $w^\pm$ で負圧・NaN) → TRRS (内部エントロピー) で常に正、M2 (内部超音速流出でも外気で流束が 19 % 変わる) → 内部 $M_{n,i}\ge1$ で $U_R=U_i$ に滑らかに寄せ最後に適用、M3 (自由流の超音速流入判定で跳び 0.215) → 滑らかな重み、m4 → §5.1 #1g と README を同期 |
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-9.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-9.md) | NO-GO, C0/M2/m1 | **全件採用、試作で確認**: M1 (TRRS も真空生成で非正、試作のクランプは float32 で NaN・未計数) → 真空判定 $B\le10^{-3}(c_i+c_{po})$ で $U_R=U_i$ に置換して計数・V0u(iii'') に float32 と独立掃引、M2 (混合変数が未定義で EOS とスカラーの状態が分離) → 原始変数 $(\rho,\mathbf u,P,Y,k,\omega)$ を同じ重みで混合・流入スカラーは $\phi_R$・V0u(iii''')、m3 → methods 一覧表を同期 |

## 7. 影響範囲

- ソルバ: `solver_density_cuda/boundaryCond.{hpp,cpp}`、`cuda_forge/boundaryCond_d.cu` (自由流状態の起動時計算)、`cuda_forge/convection/convectiveFlux_d.cu` (境界ループの振り分け)、
  新規 `cuda_forge/convection/farfieldFlux_d.inc.cuh` (HLLC・外側状態・`farfield_flux_d`)、`cuda_forge/scalarTransport_d.cu`・`cuda_forge/passiveKernels_d.cuh`・`cuda_forge/speciesTransport_d.cu` (farfield 面値)、
  `cuda_forge/ransBoundary_d.cu`、`cuda_forge/convection/convectiveFlux_d.cu` の帳簿ダンプ拡張。
- ツール: `solver_density_cuda/tools/farfield_proto1d.py` (ホスト試作・ゲート)、`farfield_balance.py`、`ref1d_euler_tp.py`。
- 設計側: `design/forge_design/evaluate/runner_sern3d.py`、`design/forge_design/metrics/sern_momentum.py`。
- 既存ケース: 変更なし (新種別を書かなければ既存経路はビット不変、V0)。
- docs: `methods/boundary.md` (一覧表・理論節・ディスパッチ表)、`methods/index.md` (見出しのみ)、`procedures/recommended-settings.md`。

## 8. 完了条件

- [ ] 関連 `methods/` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6 V0–V3)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`
- [ ] [`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-28` — plan-9 NO-GO (C0/M2/m1) を全件採用。真空・近真空は $U_R=U_i$ に置換して計数、混合は原始変数を同じ重みで、流入スカラーは混合後の外側状態の値。試作に float32 と独立掃引のゲートを追加し合格。
- `2026-09-28` — plan-8 NO-GO (C0/M3/m1) を全件採用。外側状態の音響部を TRRS (内部エントロピー、常に正) に、超音速の境目を滑らかな重み (内部の超音速流出を優先) に。試作で plan-8 の全ゲート合格。
- `2026-09-28` — plan-7 NO-GO (C0/M2/m1) を全件採用。外側状態の音響部を境界節点の $Z$ で線形化した特性量、密度・組成は常に自由流側、波速は Davis + 検査 + HLL 退避。試作で全項目合格。「ゴースト」の表記を node の弱形式の外側状態に改めた (ユーザ指摘)。
- `2026-09-28` — ホスト 1D 試作で 5 候補を比較し、自由流ゴースト + 境界面だけ HLLC に決定 (§4.2 の表)。§5・§6 を同期 (V0p・V2d の判定分離・V2f の判定対象)。
- `2026-09-28` — plan-6 NO-GO (C0/M3/m1) を全件採用。机上の反復をやめ、ホスト 1D 試作で境界流束の候補を比較してから §4 を決める (§5.1 #1g)。
- `2026-09-27` — plan-5 NO-GO (C0/M2/m3) を全件採用。境界状態の組み立てをやめ、外側に自由流を置いて SLAU で解く方式に変更 (frozen-γ も不要に)。
- `2026-09-27` — plan-4 NO-GO (C0/M3/m2) を全件採用。境界閉包を Riemann 不変量の混合から局所線形化の特性振幅に変更 (接触波で擾乱を出さない)。
- `2026-09-27` — plan-3 NO-GO (C0/M6/m1) を全件採用 (判定用音速の固定、合格条件の修正)。
- `2026-09-27` — plan-2 NO-GO (C0/M5/m2) を全件採用し全面改訂 (ユーザ決定: 検証はレビューどおり全部)。
- `2026-09-27` — codex plan 段 NO-GO (C1/M6/m1) を全件採用し §4/§6 を改訂。
- `2026-09-27` — 初稿 (ユーザ「遠方境界入れたらすっきりかもね。やってみますか」)。`methods/boundary.md` に理論節を追加し、`outflow` の説明 (実装は全量コピー) を訂正。
