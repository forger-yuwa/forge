# 特性型の遠方境界 `farfield` (node)

## メタ

- **area**: `boundary`
- **status**: `in_progress`
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
- `sstEnergyIncludesK`: エネルギー・圧力に $k$ を含める経路では、**L = $k_i$、R = 混合・置換後の外側状態の $k_R$** (§4.2。帯の外では $k_\infty$) で $E^*=E+\rho k$、$p^*=p+\tfrac23\rho k$ として HLLC に渡す。移流する $k$ の面値 (流入時) も同じ $k_R$ (plan-10 M1: 右側のエネルギーを $k_\infty$ で作ると混合帯でエネルギー流束が 3.4 % ずれる)。
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
| 1l | ~~plan-10~~ **済 (2026-09-28 GO-with-changes C0/M2/m1、全件採用)** → 実装着手可 (V0h ホストゲートから) | | F |
| 1j | ~~V0h ホストゲート (TP 拡張)~~ **済 (2026-09-29、合格)**: `solver_density_cuda/tools/farfield_proto1d_tp.py` (SERN の EXH/AIR 2 擬似種、NASA-9、TP SLAU 内部 + TP HLLC 境界)。超音速流出の極限 (外気の速度 −2〜2 a・圧力 0.2〜5 倍を掃引) で $\max|F/F(U_i)-1|=0$ (外気そのまま + SLAU 16、+ HLLC 40)。連続性 (高温排気の内部、1e-4 a 刻み) 隣接差 ≤ 4e-4、退避 0。正値性 9,477 点で非有限 0・HLL 退避 0・真空置換 242 (極端な離反のみ)。音響反射: 一様な外気 M 0 0.50/0.24 %・M 0.3 0.15/0.09 % (Δx 5/2.5 mm)、**高温排気 (T 600 K・Y 0.13) を出ていく音響 0.16/0.10 %** (外気そのまま + HLLC 32/39 %)、パルスなし対照 3e-11。接触波で境界が加える誤差 4.4e-6/1.3e-6 P (長領域自身 8.0e-6/1.5e-6 P)。V0p (保存形の混合 0.9:0.1) の圧力誤差 +0.654 % (境界と無関係、記録) | O |
| 2 | ~~実装 (§5 の 1–5)、V0・V0u・V0k~~ **済 (2026-09-29)**。V0: farfield を含まない SERN 構成 (run_0971 設定) で初回 massflux・状態ダンプが旧バイナリとビット一致 (case/46 run_0991_ff_v0_old/new)。V0u/V0k: `tests/unit/test_farfield_flux.cu` 全件 PASS — (iv) 一様 1.3e-16、(ii) 接触波 9e-16、(i) 連続性 隣接差 1.1e-5 (≤1e-4)・置換/退避 0、(iii) 243 組合せで帯外の外側値 = 外気・有界、(iii'') 独立掃引 f64/f32 で非有限 0・負値 0 (真空置換 CPG 2144/16875・TP 54/10125、HLL 退避 0)、(iii''') $Y_R$ 0.065・$k_R$ 0.075、V0k 56 面×CPG/TP で相対 2.7e-8。**実装中の不具合**: 外側状態の構造体代入 `R = f` が GPU 上で組成 Y と k・ω を 0 にしていた (ホストでは正常)。明示コピー `ff_copy` と値初期化に替え、V0k で旧版 29/56 面不一致 → 修正版 0 を確認。**V0u (v) 済 (2026-09-29)**: `case/58.farfield_verification/v0u_reject.py` — cell・ROE・凝縮 (CPG N2 担体)・トレーサ・遷移・軸対称・周期共有の 7 構成が拒否メッセージで非 0 終了、対照は完走 (`V0U_REJECT_VERDICT.txt`) | AWS でビルド | O |
| 3 | 独立参照解 (§5 の 6) と V1–V2。**V1 済 (2026-09-29、全 4 本合格、置換・退避 0)**: (a) CPG M0.5 7.5e-7 (`case/58.farfield_verification/run_0001_v1a`)、(b) TP M6 6.7e-6 (`run_0005_v1b_fix`)、(c) 30° 傾き 6.3e-6 (`run_0006_v1c_fix`)、(d) SST 輸送残差の打消し 5.8e-7 (`run_0004_v1d`、帳簿 全節点)。(b)(c) のずれは step 0 で 1e-6 (float32 EOS 往復) → 500 step で 6e-6 に達し以後横ばい。旧 `run_0002_v1b`/`run_0003_v1c` は不具合版 (破棄予定)。**V2a 済 (2026-09-29、全判定合格)**: dual-time、Δx 5 mm・dt 1.25e-6 (音響 CFL 0.11)・nSub 20・cfl_pseudo 12。反射率 M0.3 SLAU **0.099 %** (`run_0040`/`0041`)・SLAU2 0.118 % (`0042`/`0043`)・M0 0.203 % (`0048`/`0047`)、対照 M0 slip 91.3 % (`0046`)。時間精度 nSub 40 との差 0.15 %・dt/2 との差 0.58 % (≤ 1 %)。置換・退避 0。**経緯**: (1) 初回 dt 5e-6 (`run_0010`–`0018`) は dt/2 との差 7.6 % で時間精度不足、(2) dt 1.25e-6 (`run_0020`–`0028`) で差の最大が入射パルス通過時刻 (反射でない) に 1.1 % 出た → 原因は**リミッタ基準値の自動決定** (`limiterRefLength` = 領域の対角長、ρ・P・a = 初期場平均) が短 (L 1.0)・長 (L 3.0) で違い離散化が別物になっていたこと。自由流の値で明示固定して解消。プローブ出力の桁 (既定 6 桁で P が 1 Pa 刻み) も 10 桁にした (`probe/point_probes.cu`、出力専用)。**V3 への含意**: 幅系列も同じ自動決定なので、V3 では `limiterRefLength`・`limiterRoRef`・`limiterPRef`・`limiterARef` を全幅で同一値に固定する (R4d の幅系列は未固定 — 4.35/6.19 H の差 C_L 2.5e-5 から主因ではないと見るが未検証)。<br>**V2b 済 (2026-09-29、6 本 × 3 判定 = 18 件合格)**: `run_0030`–`0035` (+ `_bal` = 最終場から restart した 1 評価)。離散恒等式 (初回評価) 全成分 |差|/規模 ≤ 1e-9、Σ種 = 質量 相対 ≤ 2e-8、定常後の全体収支 |R|/代表量 ≤ 6e-9、Y_EXH 最大 5e-36 (置換)。残差は float 床まで 4.4 桁低下 (`rms_roY0` は恒等的 0 で判定不能 = 完全置換)。**帳簿の読み方の訂正**: 各残差配列は担当段の中でゼロに戻されるので、2 回目以降の評価では「段の前後の差」は寄与にならない (段の前に前回の最終残差が残る)。`eval_v2b.py` は担当段直後の値を寄与とする。<br>**V2f 済 (2026-09-29、合格)**: `run_0036_v2f` — 評価 1 で xmax の 169 面すべてが ṁ<0、その面の外側組成 = 外気、化学種の恒等式で流入面の Y = 外気 (Σ 1.6e-10)・流出面 (評価 2・3) で内部値 (相対 1e-10)、置換・退避 0、定常後は自由流に戻る (Y 1e-36)。<br>**V2e 済 (2026-09-29、合格)**: 音響 + SST dual-time 反射 0.12 % (ek0、`run_0052`/`0053`)・0.33 % (ek1、`0054`/`0055`)、V1b 陽解法 (RK3、局所 dt) ずれ 8.7e-6 (`run_0050`)、V2b(i) 陽解法 恒等式・収支・置換 PASS、残差 5.7 桁低下 (`run_0051` + `_bal`)。<br>**V2c FAIL (2026-09-29、保持)**: `run_0063_v2c_A`/`0064_B`/`0065_C` (ランプは 0.7 m で終わる形。出口まで伸ばした `run_0060`–`0062` は A が流路閉塞で発散、破棄)。B vs C 最大 0.087 Δp (x 0.21、ランプ角直後)、角・出口端を除いても 0.048 Δp (基準 0.02)。A vs C 2.61 Δp (対照 OK)。角の値は A = B (6 桁一致)・C だけ違う、cfl 0.5 (`0066`/`0067`) と再変換 (`0068`) で不変。上面の境界節点で B が C より 0.10 Δp 低く、弱い反射波 (内部 0.06 Δp) が出口寄りの下面に届く。**原因未分離** (codex diagnose 採用)。<br>**V2d**: 接触波 V2d-1 は 3 物性 PASS (`run_0070`/`71`・`75`/`76`・`80`/`81`)。**V2d-2 パルスなし対照 FAIL (保持)**: cpg 2.4e-4 Pa PASS、tp1 1.64 Pa・tp2 2.49 Pa (基準 0.0285 Pa)。長領域でも同じ (`run_0085` 2.4888 Pa)、左端の自由流を内部と同じ高温にすると 0.0000 Pa (`run_0086`/`0087`、隔離対照として別記) → 左端流入で作られる接触面の TP 保存形混合が音源という仮説 (未確定)。V2d-2 の反射 (0.07–0.23 %) は左端擾乱が反射より先に評価点に届くため隔離試験になっていない → 再設計。 | AWS | O |
| 3a | **V2c 角の作用素 A/B** (codex diagnose 2026-09-29): 低領域 (A 格子、上面 slip) と高領域 (C 格子、上面 slip) に、C の最終場を座標対応でそのまま与え、1 step の初回評価だけを比較。対象 = 下面 x 0.195–0.230 と再構成に要る近傍。照合 = 双対体積・面積ベクトル・`after_eos_bc` 状態・面の再構成状態と流束 (帳簿 `.faces`)・`res_after_conv`。**判定 (事前固定)**: 各保存量の残差差 / その節点に接する面流束の絶対和 ≤ 1e-5 (ゼロ規模は自由流基準)。超えたら局所作用素の交絡を採用し最初に異なる幾何・再構成・流束へ絞る、全て許容内なら作用素差を棄却し反復過程へ移る。どちらでも収束・境界性能は判定しない。**結果 (2026-09-29)**: 作用素は同一 (`run_0092`/`0093`、幾何・状態・面の差 0、残差差 ≤ 7.7e-8、`V2C_OPAB_VERDICT.txt`) → 第 1 仮説棄却。C の最終場から低領域 B・A・C 自身を 6000 step 継続 (`run_0094`–`0096`) → 角は 3 本とも C の値を保持。同じ初期場での B vs C は全線 0.0210 Δp (出口端)、除くと 0.0080。codex diagnose 2 回目 ([記録](../../notes/reviews/2026-09-29-farfield-v2c-corner-and-v2d2-redesign-diagnose.md)、全件採用) の絶対残差 A/B (`run_0097`/`0098`、同じ B 格子に 2 つの最終場、事前閾値 \|R\|/Σ\|F\| ≤ 1e-5、slip 節点の運動量は接線成分、運動量は 3 成分共通の規模、規模 0 は自由流基準): **角 (x 0.195–0.23) は 2 状態とも全 5 成分 ≤ 5.2e-7**、しかし**場全体では約 2000 節点が超過し最大はランプ終端の凸角 (x 0.715、下面) で 4.6e-3 / 7.3e-3** → 規則どおり「2 つとも離散定常解」を棄却、非零残差が残る凸角を追う (`V2C_RESAB_VERDICT.txt`)。**次**: 凸角の停滞の扱いを諮る (#3c) | AWS | O |
| 3c | **V2c 凸角 (ランプ終端 x 0.7) の残差停滞**: 形状側 (slip の凸角) の問題として扱うか、境界と無関係に node slip 凸角の既知欠陥か、を切り分ける。案を codex に諮ってから §6 V2c の手順 (共通初期場の追加隔離試験を含む) を事前登録し直す | AWS | F |
| 3b | **V2d-2 の再設計** (codex diagnose 2026-09-29、2 回目で採用確定: `--left-hot` を右端流出の隔離試験として正式化、短・長の両方でパルスなし対照、3 物性、評価器は入射窓で振幅・反射到達窓で短長差を測る形に直す): 右端流出を隔離する対照 (左端 = 内部と同じ高温状態) を正式な隔離試験として別記、元の配置の FAIL は保持。反射の評価は左端擾乱の到達前に反射が評価点へ着く配置に直す (§6 の該当節を改訂してから回す)。左端流入の TP 混合仮説は境界を含まない試験で確かめる | AWS | F |
| 4 | SERN V3 | AWS。R4d へ反映 (restart は `r4d_common_restart.py` 型の index コピー、メッシュ品質、判定区間、case README の run 索引) | O |
| 5 | codex result 段 → accepted | | F |

## 6. 検証

すべて AWS。合否は実行前に固定。判定区間を明記し `check_convergence.py`・`check_quasisteady.py` の出力原本を run に残す。
**非有限流束の置換 (§4.3) は各試験の評価区間で 0 回**を必要条件とする。
V1–V2 は**境界機能の受入れ**、V3 は**SERN での配置 (側方幅) の採否**。

- **V0 既存境界の不変**: farfield を含まない構成 (run_0971 設定) で、同一初期状態からの初回 `massflux` と状態ダンプが変更前バイナリとビット一致。更新後保存量は旧バイナリ 3 回反復の再現性幅以内。
- **V0h ホストゲート (本体実装の前、plan-7 推奨)**: `farfield_proto1d.py` を TP (単成分・多成分 lump、既存 `FrozenGas`) に拡張し、§4.2 の表の全項目 + V2d-2 (TP の高温内部を出ていく音響、反射 ≤ 5 %) + 波速順序・星状態の妥当性・退避 0 を確認してから CUDA に進む。
- **V0p TP の保存形混合 (前提、境界と無関係)** (plan-6 M2): 単一 CV に高温側 (T 600 K、Y 0.13) と外気 (T 220 K、Y 0) の保存量を 0.9:0.1 で混ぜ、EOS から復元した圧力の誤差を記録する (期待 +0.65 %)。これは境界の合否に入れず、V2d の「長領域自身の誤差」の説明に使う。
- **V0u 単体試験** (ホスト、HLLC の `__host__ __device__` 関数): (i) **連続性**: $U_{n,i}$・$\rho_i$・$P_i$・接線速度を、法線 Mach ±1・0 付近で 1e-6 刻みに掃引し、面流束 5 成分・スカラー面流束の隣接点の差 ≤ 1e-4 × 流束規模 (plan-2〜5 の反例入力をすべて含む)、(ii) **接触波**: CPG と TP で P・$U_n$ 同じ・T 600/220 K・Y 0.13/0 の内部/外気、流出と流入の両方で面流束 = 風上側の物理流束 (相対 1e-6)、(iii) **流向の全組合せ**: $Q_n$ の符号 × 実際の $\dot m$ の符号 × ゼロ通過で、$\dot m<0$ のスカラー面値が外側状態の値 (帯の外では外気値)・$\dot m>0$ で内部値、(iii'') **真空と float32** (plan-9 M1): 真空の反例と、速度 ($U_{n,i}$ −3〜3、$U_{n,\infty}$ −3〜10)・密度比・圧力比の独立掃引を float64 と float32 で行い、非有限 0・外側状態の負値 0 を確認、真空置換と HLL 退避の回数を別々に記録 (有限であることと置換が無いことを分けて判定)、(iii''') **混合帯の逆流**: 内部 $(1,0.95,1/\gamma,Y\,0.13)$・自由流 $(10,-3,10/\gamma,Y\,0)$ で、流入の種の面値が外側状態の $Y_R$ (0.065) と一致すること、TP・SST (k 込み) でも同じ。`sstEnergyIncludesK: 1` で $k_i=0.14$、$k_\infty=0.01$ を与え、運動量・全エネルギー・$\rho k$ の流束が同じ $k_R$ (0.075) から作られていること (plan-10 M1)、(iv) 一様 ($U_i=U_\infty$) で $F_{num}=F(U)$ (相対 1e-6)、(v) 非対応構成の起動拒否 (cell・ROE・凝縮・トレーサ・遷移・軸対称・周期共有) がそれぞれエラー終了。
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
  - **配置の定義 (2026-09-29、実行前に確定)**: 代表量が $|\mathbf u_\infty|$ を含むので静止外気は使わない。自由流 = SERN 外気 (TP、220 K、2851 Pa、Y_EXH 0、SST k 479.653・ω 119844)、M 0.5 の通過流。(i) **流出配置** = +x 方向 (xmin だけ流入、xmax 流出、側面 4 面は $Q_n=0$ の接線流 = SERN 側面と同じ状況)、(ii) **流入配置** = 方向 $(\cos30°\cos20°,\sin30°\cos20°,\sin20°)$ (xmin・ymin・zmin の 3 面から斜めに流入、残り 3 面から流出)。初期値は内部全域を P∞・u∞ のまま T 600 K・Y_EXH 0.13・k・ω 10 倍。4 本 = (i)(ii) × `sstEnergyIncludesK` 0/1 (化学種は `speciesFaceReconstruction` 0)、加えて化学種の 2 経路用に (i)(ii) × `speciesFaceReconstruction` 2 (`sstEnergyIncludesK` 0) の 2 本。離散恒等式は初回評価 (非一様の初期場) の帳簿 + 面ダンプで、定常後の全体収支は最終場から同一格子 restart した 1 評価の帳簿で見る。
  - **離散恒等式 (毎評価)**: 全節点の対流残差の和 = −(farfield 面流束の和) (内部面は相殺)。診断ダンプ (§4.3) と帳簿ダンプ (`FORGE_DUMP_LEDGER`、全節点に印) で照合。**判定 (plan-4 m5)**: |差| ≤ 1e-5 × (その保存量の全面流束の絶対和) (正味和で割らない。float32 の atomicAdd 積算の丸め幅を含む。ホスト集計は double)。
  - **定常後の全体収支** (plan-3 M4): 各保存量で残差の正本を $R=-\sum_{faces}F+\sum_{nodes}S$ (外向き流束は残差に $-F$、生成ソースは $+S$) とし、定常で $R\to0$ を見る。$S$ は $k,\omega$ の体積ソース (帳簿の段別残差 `res_after_sources` − `res_after_species`) で、他の保存量は 0。`sstEnergyIncludesK: 1` の全エネルギーは `res_roe` が既に全エネルギーの残差なので **`res_roK` を足さない** ([`turbulence-sst-energy-includes-k.md`](../accepted/turbulence-sst-energy-includes-k.md) §分割保持)。規格化は**事前登録した非零の代表量** (ゼロ成分を落とさない): 質量・化学種 $\rho_\infty|\mathbf u_\infty|A$、運動量 $(\rho_\infty|\mathbf u_\infty|^2+P_\infty)A$、エネルギー $\rho_\infty|\mathbf u_\infty|H_\infty A$、$k$ $\rho_\infty|\mathbf u_\infty|k_\infty A$、$\omega$ $\rho_\infty|\mathbf u_\infty|\omega_\infty A$ ($A$ = farfield 面の総面積)。**$|R|$/代表量 ≤ 1e-4** を全成分に。
  - (ii) で内部が自由流へ置き換わること ($|Y-Y_\infty|\le10^{-4}$)。
  - **化学種の 2 経路** (plan-10 M2): 流入・流出の両配置を `speciesFaceReconstruction: 0` (一次輸送) と `2` (S3、`convMethod ≥ 1`・`speciesImplicitCoupling: 1`) の両方で回し、実効設定を run に保存する。各経路で、流入組成が外側状態の値であること、Σ種流束 = 質量流束 (相対 1e-6)、残差との離散恒等式を確認する (一次経路だけ直した実装を通さないため)。
- **V2f 自由流分類と局所流況が違う面** (plan-5 M2): V1 の箱で、自由流は面から超音速流出 ($Q_n=2a_\infty$) だが、内部の初期状態を局所に逆流 ($U_{n,i}=-0.2a$) させた帯を置く。 **配置 (2026-09-29 確定)**: TP 外気 +x M 2 (xmax で $Q_n=2a_\infty$)、初期値は x ≥ 0.8 m の帯だけ $U_x=-0.95a$・Y_EXH 0.13 (他は自由流)。**逆流の強さはホスト関数で決めた** (CPG、同密度・同圧): $Q_n=2a_\infty$ の自由流は TRRS の J⁻ で外向きに強く引くので、$U_{n,i}=-0.2a$ では $\dot m=+0.29\rho c$ (流出のまま)、−0.5a で +0.13、−0.8a で −0.0065、**−0.95a で −0.018** (置換・退避 0)。帳簿を最初の 20 評価記録し、実際に $\dot m<0$ となった xmax の面を面ダンプから数える。合格 (plan-6 M3): **実際の $\dot m<0$ の面**で化学種面値が外気値・$\dot m>0$ の面で内部値 (内部速度の符号ではなく流束の符号で判定)、NaN・置換 0、定常後は自由流へ戻る。逆流帯の強さは、HLLC で実際に $\dot m<0$ になる面が生じることをホスト関数で事前に確かめて決める。
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
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-10.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-10.md) | **GO-with-changes**, C0/M2/m1 | **全件採用**: M1 → §4.3 の SST 右状態を混合後の $k_R$ に統一・V0u(iii''') に k 込みの照合、M2 → V2b に `speciesFaceReconstruction` 0/2 の両経路、m3 → methods・説明ページの非反射の主張を「法線方向の小振幅音波」に限定し、斜入射の線形反射係数 $(1-\cos\theta)/(1+\cos\theta)$ (45° で約 17 %) を明記。codex はホストで反射率・独立掃引を再現 |
| diagnose | 2026-09-29 | [2026-09-29-farfield-v2c-v2d-fail-diagnose.md](../../notes/reviews/2026-09-29-farfield-v2c-v2d-fail-diagnose.md) (brief [2026-09-29-farfield-v2c-v2d-fail.md](../../notes/reviews/briefs/2026-09-29-farfield-v2c-v2d-fail.md)) | V2c・V2d-2 の FAIL を保持、Major 4 | **全件採用**: (1) V2c を「境界の実力 ≈6 % 反射」と確定しない (評価器は反射係数を分離していない)、(2) 角を評価から除外しない (原因を同定・修正して全線で判定し直す)、(3) Z = ρU/√(M²−1) への変更は現段階で却下 (上面 Mₙ≈0.36 で小振幅・接線流の仮定外)、(4) V2d-2 を「右端擾乱ゼロ」で合格にしない (left-hot は隔離対照として別記)。V2d-2 の反射 PASS も左端の擾乱が反射より先に評価点に届くため隔離試験でない → 再設計。次の一手 = V2c 低・高領域に同一状態を与えた 1 評価の作用素 A/B (§5.1 #3a) |

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

- `2026-09-29` — V2e 合格。V2c・V2d-2 が事前基準で FAIL → codex diagnose (全件採用、FAIL 保持、§5.1 #3a・#3b)。R4d 幅系列はリミッタ自動基準値 (p_ref 最大 10 %・L_ref) が幅ごとに違い、別バイナリ (`sglsq/forge_2fa3826c`) だったことを記録 (V3 は新バイナリ・基準値固定で回す)。
- `2026-09-29` — V0u (v)・V2a・V2b・V2f 合格 (§5.1 #2・#3)。V2a でリミッタ基準値の自動決定が領域長に依存し比較を汚すことを発見 → V3 で固定する。
- `2026-09-29` — GPU 実装完了。V0 (ビット一致)・V0u/V0k (単体試験 `tests/unit/test_farfield_flux.cu`)・V1 (a)–(d) 合格。GPU 上の構造体代入で組成・k・ω が落ちる不具合を V1b で検出し修正 (§5.1 #2)。
- `2026-09-29` — V0h (TP ホストゲート) 合格 (§5.1 #1j)。GPU 実装に着手。
- `2026-09-28` — plan-10 **GO-with-changes** (C0/M2/m1) を全件採用。SST の右状態を混合後の $k_R$ に、化学種の一次・S3 両経路を V2b に、非反射の主張を法線方向に限定 (斜入射 45° で線形理論でも約 17 %)。status を in_progress へ。
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
