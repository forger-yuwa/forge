# 共役熱伝達 (CHT) の導入方針 — 保存的な界面契約と段階連成

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) — 壁 / 等温壁 (cell ゴースト・node 壁ノード T ピン) の現在仕様
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md) — 残差組立てと時間積分 (連成フックの位置)
  - [`methods/diffusion.md`](../../methods/diffusion.md) — 粘性・熱流束の離散 (面物性・非直交補正)
  - [`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md) — SU2 比較手順 (**CHT は multizone なので手順の追加が要る**, §6 V2)
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) §1・§1.2 — 現行レシピと段階起動
- **related_plans**:
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) §4.6-4 / §5.1 #6・#7 — **本計画の親**。弱 CHT ループと `wallProfile` を移管済み
  - [`tooling-energy-balance-diagnostics.md`](tooling-energy-balance-diagnostics.md) — **本計画の合格条件** (codex C1 採用で「audit」から格上げ)。拘束反力を含む実効界面熱量はこれが無いと取れない
  - [`case-hypersonic-gap-heating-validation.md`](case-hypersonic-gap-heating-validation.md) — 発注元の 1 つ。壁温実測を `wallProfile` で与える**流体側検証**を担う
  - [`case-plate-annular-cavity-m5.md`](case-plate-annular-cavity-m5.md) — case/49。**別セッション所有につき参照のみ**
- **related_notes**:
  - [`notes/investigations/cht-validation-case-survey.md`](../../notes/investigations/cht-validation-case-survey.md) — 一次適用先の公知データ調査 (§4.9 の根拠)
- **created**: `2026-09-19`
- **owner**: `sano`

## 1. 目的

forge の壁は現在 **断熱**か**等温 ($T_w$ 既知・bcond 単位の定数)** の 2 択で、すきま / 深いキャビティ・冷却ノズル壁では
$T_w(x)$ は**与件ではなく解の一部**である。完了時に得られる状態:

- **保存的な界面契約が固定されている**: 連成に渡す熱量が「物理境界流束 + 拘束反力」として定義され、
  流体・固体・接合部の熱収支が**定量ゲート**として閉じている (codex C1)。
- **弱連成 (Phase 1) が回る**: 外部ループで共役解が得られ、1 次元純伝導の**解析解**で検証済み。
- **ソルバ内連成 (Phase 2) が回る**: 1 run 内で共役定常解が得られ、Phase 1 と一致する。
- **公知試験データで当たっている**: 超音速出口のタービン翼 (§4.9) で、測定された壁温分布を予測できる。

## 2. スコープ

- **やる**
  - **Phase 0**: 保存的な界面契約の確立 (実効界面熱量・符号・帰属・接合)、診断出力、`wallProfile` (壁温分布入力)。
  - **Phase 1**: 外部弱連成ループ。固体バックエンドは **`local1d` / `shell2d` / `fem2d` (一般 2D 領域)** の 3 つ。
  - **Phase 2**: ソルバ内連成 (`wall_isothermal` + `ints: {conjugate: 1}`、**node 限定**、host 側固体求解)。
  - 検証 V1–V6 と、界面専用の収束ゲート。
- **やらない** (別 plan / 条件付き)
  - **Phase 3 = 固体ゾーンをソルバ内で解く**: §4.8 の基準を満たしたら別 plan。
  - **cell 離散化での連成の実装保証** (codex M4): cell は `vizBfaceNodes` が空でシェルとの 1 対 1 対応を主張できない。
    `wallProfile` (面重心補間) までは cell も対象、連成は node 限定。
  - 表面間放射、非定常 CHT (thin-skin 過渡)、熱応力、接触熱抵抗の同定 (§4.10)。

## 3. 前提 — 既存資産の棚卸し (本計画 + codex レビューで実機確認)

| 事実 | 出典 | 影響 |
| --- | --- | --- |
| `Ts` は **per-face bvar**。起動時に YAML の一様値で 1 度埋め、以後カーネルは書かない | [`mesh/mesh.cpp`:64](../../solver_density_cuda/mesh/mesh.cpp) | 面ごとに違う $T_w$ がそのまま効く。連成の入口 |
| cell ゴースト ([`boundaryCond_d.cu`:308](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu)) も node 壁ノード T ピン ([`nodeWallDirichlet_d.cu`:75](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu)) も `Tsb[ib]` を読む | 同左 | CUDA カーネルのアルゴリズム改修なしで壁温分布を課せる |
| **`applyInletProfiles` は種別を限定しておらず、`inletProfile` フラグだけを見る** | [`boundaryCond.cpp`:244](../../solver_density_cuda/boundaryCond.cpp) (codex M5) | 一般化は**小さい**。ただし補間は**面重心** ([:312](../../solver_density_cuda/boundaryCond.cpp)) で、node が温度を課すのは**ノード位置**。case/48 `run_0011` で実測 **0.679 mm のずれ** |
| 壁ノードは T ピン後に `res_roe` を 0 化する | [`nodeWallDirichlet_d.cu`:90](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu) | **実効の熱授受は拘束反力を含む**。片側差分でも拘束反力単独でもない (codex C1) |
| 解像壁の熱流束は node = $k_{\rm eff}\nabla T\cdot\mathbf S$ / cell = ゴースト差分。内部 W–I 面は**面物性 + 非直交補正**の別離散 | [`viscousFlux_d.cu`:529](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu), [:223](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) | 「壁の熱流束」に複数の離散が混在する。§4.3 で正本を 1 つ決める |
| 壁種別は**名前で直書き判定**されている: `iso_wall_flag` [`mesh.cpp`:929](../../solver_density_cuda/mesh/mesh.cpp) / T ピン [`nodeWallDirichlet_d.cu`:130](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu) / 粘性壁 [`viscousFlux_d.cu`:954](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) / 壁距離 [`calcWallDistance_kdtree.cpp`:129](../../solver_density_cuda/input/calcWallDistance_kdtree.cpp) / DPLUR エネルギー行切離し [`timeIntegration_d.cu`:1025](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) | **新種別 `wall_conjugate` を足すとこれら全部から漏れる** (codex M4) → 種別は `wall_isothermal` のまま、連成は属性 |
| node は bcond ごとに境界**ノード**単位で 1 bplane。半割面ベクトルは**合算**され、面積は合成ベクトルの大きさ | [`gmshReader.hpp`:1618,1641](../../solver_density_cuda/mesh/gmshReader.hpp), [:2250](../../solver_density_cuda/mesh/gmshReader.hpp) | 角・曲面で $|\sum_f\mathbf S_f|\ne\sum_f|\mathbf S_f|$。**シェルの幾何を node の合成量から作ってはいけない** (codex M6) |
| primal 境界面の接続は `bc.vizBfaceNodes` (node のみ) | [`mesh/mesh.hpp`:93](../../solver_density_cuda/mesh/mesh.hpp) | facet 単位の幾何・荷重の正本にする (温度 DOF の共有とは分離) |
| **既存 run に監査材料が無い**: case/48 `run_0011` と case/44 `run_0114` の壁ダンプは `qwall` が全値 0、$T_1,d_1,k_{\rm eff}$・拘束反力は未保存。両 run とも既定条件で `NOT CONVERGED` (`run_0114` は `rms_roUy`/`rms_roOmega` が RISING) | codex M9 (保存済みダンプを全時刻確認) | §5.1 #2 は**新規 run で**やる。旧 run を合格根拠にしない |
| ただし codex の実測: `run_0011` $x\approx0.5$ m で **コンパクト差分 96.184 / 2 次片側差分 98.820 kW/m² (差 2.67 %)**、両者 `STEADY` | codex M9 | **差分形式の選択だけで V2 の 3 % 許容に迫る**。界面定義を先に固定する理由 |
| `check_convergence.py` は既知の保存量列しか見ない。`stage_manifest.py` は config テキストしか見ない | [`check_convergence.py`:45](../../solver_density_cuda/tools/check_convergence.py), [`stage_manifest.py`:75](../../solver_density_cuda/tools/stage_manifest.py) (codex M10) | 界面残差の列を足すだけでは検査されない。外部入力 (`wallProfile`/`solid.json`) の変更も区間に出ない |
| SU2 は CHT を持つが **multizone** (固体ゾーン + `MARKER_CHT_INTERFACE`)。既定は `DIRECT_TEMPERATURE_ROBIN_HEATFLUX`。界面転送は $(T_w,q,k,k/d)$、$dT/dn=(T_w-T_1)/d$ のコンパクト差分 | [`CConfig.cpp`:2786](../../.external/su2-src/Common/src/CConfig.cpp), [`CConjugateHeatInterface.cpp`:58](../../.external/su2-src/SU2_CFD/src/interfaces/cht/CConjugateHeatInterface.cpp) | 先例かつ外部基準。**ただし `CNSSolver.cpp`:593 の引数順 (`thermal_conductivity, dist_ij`) と :718 の呼出しが逆** (codex M8)。`AVERAGED_*` の係数に効くので、比較は `DIRECT_*` で行う |
| 閉じた系を定常局所 dt で回すと音響過渡で数 step 発散する | [[steady-localdt-acoustic-instability]] | **V1 (純伝導の閉じた箱) は dual-time で定常化させる** |
| float32: 1000 K 付近の 1 ULP は約 $6.1\times10^{-5}$ K | codex M10 | $\Delta T_w=0$ でも収支は閉じていないことがある。ゲートは温度量子化で必要熱流束を解像できるかを含める |

## 4. 設計方針

### 4.1 作戦の比較 — 何を採り、何を採らないか

| 作戦 | 中身 | 判定 |
| --- | --- | --- |
| **A. 外部弱連成** | forge run → 界面熱量 → 外部固体求解 → `wallProfile` → warm restart、を反復。固体は `local1d`/`shell2d`/`fem2d` | **採用 (Phase 1)**。ソルバ改修を最小にして物理・更新式・検証を確立する |
| **B. ソルバ内連成** | `wall_isothermal` + `conjugate` 属性。host で固体を double で解き $K$ step ごとに `bvar_d["Ts"]` を更新 | **採用 (Phase 2, node 限定)**。restart・起動・I/O の反復コストが消える。**流体の熱場緩和は消えない** (codex m11) |
| **C. 固体を別ゾーンで解く** | SU2 と同型 | **条件付き (Phase 3)**。§4.8 |
| **C'. 同一メッシュに固体セルを混ぜる** | `regionId` で分岐 | **不採用**。理由は**既存実装への侵襲性と検証費用** — 対流流束・限界子・SST・EOS 床・陰解法・出力の全経路に固体分岐が要り、全既存ケースに回帰リスクを撒く。median-dual の界面ノードでは双対 CV が両材料にまたがるので DOF の置き方に追加設計が要る (**ただしこれは原理的な不可能性ではない** — codex m11 を採用して表現を弱めた)。SU2 も頂点中心だが CHT は別ゾーン方式 |
| **D. 外部 FEM (CalculiX 等)** | A の固体側を差し替え | **A のバックエンド選択肢**。熱応力が要るときに使う |

### 4.2 界面の更新式 — **固定点を保存する形で書く** (codex M2 採用)

固体の離散作用素を $A_s$ (面内伝導 + 背面 Robin、SPD)、荷重を $b_s$、流体から固体へ入る**節点積分熱量**を
$Q_f(T)$ [W] とすると、共役定常解は

$$A_s T^{*} = b_s + Q_f(T^{*})$$

を満たす。反復は**両辺に同じ $D_f T$ を足した**次の形で行う ($D_f$ は節点対角の推定コンダクタンス [W/K]):

$$\left(A_s + D_f\right) T^{k+1} = b_s + Q_f(T^{k}) + D_f\,T^{k}$$

- **$D_f$ が何であっても固定点は上の共役解のまま**。$D_f$ は収束速度だけを決める。
  面内伝導を解いた後に局所平均を掛ける (旧稿) と固定点がずれるので採らない。
- **$D_f$ の取り方**: **初期推定値**を $D_f^{(0)} = k_{\rm eff}A/d_1$ とする。
  **これは上界ではない** (codex 2 巡目 #2 の反例): 実効界面応答 $H=-\partial Q_f/\partial T_w$ には
  壁 CV の接線輸送と内部場の応答が入り、非対角成分を持つ。
  直交 median-dual ($x=\{0,2\}$, $y=\{0,4,8\}$, $k$=1, 側面断熱, 上端固定, 下端 2 節点が界面, AR 2) では
  $H=\begin{pmatrix}1.1806&-1.0556\\-1.0556&1.1806\end{pmatrix}$ に対し $D_f=0.25I$, $A_s=0.1I$ で
  反復行列の固有値が $\{0.357,\,-5.675\}$ = **スペクトル半径 5.675 (発散)**。
  → **$D_f$ は「当たりをつける初期値」であり、安定性は受理判定で担保する** (次項)。
  **節点ごとのセカント推定 $\Delta Q_f/\Delta T$ を既定にしてはいけない**: 真の応答 $-\partial Q_f/\partial T$ は
  非対角結合を持つ行列で、成分ごとのセカントはこれを**過小評価する**。
  codex の反例 — 真の対角 $[2,2]$ に対し成分セカントは $[1,1]$ となり、
  反復行列 $(A_s+D_f)^{-1}(D_f-H)$ のスペクトル半径が **1.818 (発散)**。
  → セカントは **$D_f$ を増やす方向にだけ**使う。減らす方向は受理判定つきでのみ許す。
- **受理判定はメリット関数 + line search で行う (最大ノルムの単調減少を要求しない)** — codex 3 巡目 #1 採用。
  未緩和残差を $r^k = A_sT^k-b_s-Q_f(T^k)$ とし、**比較の間は重みを固定した**

  $$\Phi(r) = r^{\mathsf T}\left(A_s+D_f\right)^{-1} r$$

  の**降下**で受理する。$\Phi$ が減らなければ **line search** (更新量を $\beta\in\{1,\tfrac12,\tfrac14,\dots\}$ で縮める)
  → それでも降下しなければ $D_f$ を増やす → 再試行上限に達したら**停止して失敗を報告**する。
  - **最大ノルムの単調減少を受理条件にしてはいけない**: codex の反例 —
    $A_s=0.1I$, $H=\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$ (SPD), $D_f=s\,\mathrm{diag}(1,10)$, $r=(1,1)^{\mathsf T}$ で
    試行後の最大ノルムは $s$=1,2,4,1024 に対し 1.701 / 1.372 / 1.192 / 1.0008 と**常に 1 を超え**、
    $D_f$ を何度倍増しても受理できない。しかし $s$=1 の固定点反復のスペクトル半径は **0.960 < 1 で収束する**。
    **発散を止める問題ではなく、収束する反復を受理規則が止めていた**。同じ反例で $\Phi$ は $1.0081\to0.8831$ と減少する。
  - $\Phi$ の SPD 性は伝導支配の問題でのみ期待できる (一般の流体応答には保証が無い)。
    **降下しない場合の方向変更・再試行上限・量子化による停滞の扱い**を仕様に書き、
    **停滞 ($\Delta\Phi$ が丸め以下) を合格にしない**。
  - **V1b にこの反例と非一様 $D_f$ を入れる**。交番モードだけの試験では検出できない。
  - 局所最大ノルムは**最終合格ゲート (G-if)** に残す (受理規則としては使わない)。
- **旧稿の誤り (codex M2)**: 0 次元の $\omega^\star=1/(1+\mathrm{Bi})$ で「1 反復収束」と書いたが、
  $h$ と $k_f/d_1$ は同じ量ではなく、流体を再収束させれば $T_1$ も応答する。
  反例 (codex): $h=100$, $g_f=k_f/d_1=1000$, $g_s=1/R_s=100$ で増幅率 $(g_f-h)/(g_f+g_s)=0.818$ (1 反復収束ではない)。
  **「Robin だから無条件最適」とは書かない**。言えるのは「$D_f \ge$ 実効応答なら単調・非発散」まで。
- **未収束流体の応答は「小さい」とは限らない** (codex 2 巡目 #2): 第一内部点を凍結すれば応答は $kA/d_1$ に近く、
  再収束させると流体層全体の抵抗が効いて**小さくなる**場合がある。向きを決め打ちせず、$K$ 依存は実測する (§6 V6)。
- **スカラー $D_f$ だけでは非対角な流体応答を捌けない (2026-09-19 実測)**。
  `solid_shell.py` の T5 (流体応答 $H$ = 三重対角 SPD, 対角 100 / 非対角 −50 [W/K]) では、
  **素の固定点反復は 200 反復で収束せず** (最終誤差 30 K)、**Anderson 加速 (深さ 5) で 57 反復・誤差 5.6e-8 K**。
  → **加速器を既定に含める** (Anderson。Aitken はその特別な場合)。安全装置は必須:
  (i) 加速候補は**メリット関数が降下しなければ棄却**して素の反復に戻る、(ii) 係数が発散したら無効化、
  (iii) **$D_f$ を変えたら履歴を捨てる** (重みが変わると比較が濁る)。
- **$\Phi$ の比較は $D_f$ を固定した区間内で行う** (2026-09-19 実測): $D_f$ を動かしながら $\Phi$ を比べると、
  **重み変更による見かけの降下**で受理が通り、$D_f$ が発散的に増えて反復が停滞する (実測: $D_f$ が $10^8$ まで増え、
  $\Phi$ も $\Delta T$ も減るのに残差は 1 のまま)。**収束判定は物理量** (max$|\Delta T|$ [K] と
  max$|r|$/スケール) で行い、$\Phi$ は受理判定だけに使う。
- 参考 (**目安であって判定に使わない**): 0 次元の $\mathrm{Bi}=hR_s$ は薄肉金属 (0.4 mm, 15 W/mK, $h$=500) で 0.013、
  断熱材 (5 mm, 1 W/mK) で 2.5。素の Dirichlet–Neumann ($D_f=0$) は後者で発散する。

### 4.3 界面熱量の定義 — **保存的な実効熱量が正本** (codex C1 採用)

**符号規約を先に固定する** (依存先 plan [`tooling-energy-balance-diagnostics.md`](tooling-energy-balance-diagnostics.md) §4.1 と同じ):
面流束 $F^E$ は**流体 CV から外向きを正**、拘束反力 $C$ は**流体への供給を正**。このとき壁 1 節点で
**流体から固体へ入る熱量**は

$$Q_{f,i} \;=\; \underbrace{\textstyle\sum_{f\in\partial_w} F^{E}_{if}}_{\text{壁面の物理境界流束 (外向き正)}} \;-\; \underbrace{C_i}_{\text{拘束反力 (流体への供給が正)}}$$

である (外向き流束は流体から出て行く分、$C<0$ は拘束が奪った分。どちらも固体側の受取になる)。
定常の Dirichlet 行では $C_i=-R_i^{raw}$。

- **検算** (codex 2 巡目が検出): $\sum F^E=80$, $C=-20$ のとき $Q_f = 80-(-20) = 100$。
  **旧稿の $-(\sum F + C)$ は $-60$ になり、符号も大きさも誤り**だった。
- 符号規約は**散文で保証しない**。V1 (§6, 解析解で $q$ が既知) で**符号と絶対値を再現すること**を
  実装の解除条件にし、`conjugate_interface.hpp` の冒頭コメントを唯一の正本にする。

- **旧稿の誤り (codex C1)**: 「コンパクト差分形を界面の一次とする」は、壁 CV に実際に入った熱
  ([`viscousFlux_d.cu`:529](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) の再構成勾配 + 内部 W–I 面の別離散
  [:223](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu)) と一致しないので、**熱量が閉じない解を合格させてしまう**。撤回する。
- **コンパクト差分形 $k_{\rm eff}(T_1-T_w)/d_1$** (固体向き正。**旧稿は符号が逆だった**) は
  **$D_f$ の推定と精度診断にのみ**使う。
- **再構成勾配形 (`qwall_b`)** と **2 次片側差分**は診断として併記する。codex 実測で
  コンパクト 96.184 / 2 次片側 98.820 kW/m² (2.67 % 差) なので、**どれを見ているかを列名に書く**。
- **依存関係を格上げし、対応範囲のマイルストーンまで書く** (codex 2 巡目 #1):
  拘束反力の採取は [`tooling-energy-balance-diagnostics.md`](tooling-energy-balance-diagnostics.md) が提供するが、
  同 plan は初版で **dual-time・周期・軸対称を明示的に拒否**している ([:32](tooling-energy-balance-diagnostics.md))。
  一方で本計画は V1 に dual-time、V5 に周期翼列を要求し、全 run に G-cons を課している。
  **「依存先の完了」だけでは V1 も V5 も判定できない**ので、次を解除マイルストーンとして登録する:
  **(a) dual-time の診断解除を V1 の前、(b) 周期の診断解除を V5 の前**。軸対称と壁モデルは
  それぞれの解除試験ができるまで**本計画の合格対象外**と明記する
  (§6 V6 の「壁関数併用」も**解除後の追加試験**に分離した — codex 3 巡目 #2)。
  **提供側 ([`tooling-energy-balance-diagnostics.md`](tooling-energy-balance-diagnostics.md) §5.1) にも次を登録する**
  (登録しないと解除が誰の作業でもなくなる):
  1. $R^{raw}$ は**壁残差射影より前の空間残差**と固定する
     ([`main.cpp`:1814](../../solver_density_cuda/main.cpp) は空間残差を組んだ後に BDF 項を足すので、
     **後から `res_roe` を読むだけでは復元できない**)。
  2. 同一評価状態について、**実際の BDF 係数・履歴・体積**から $D_t(VE)$ を計算し $C=D_t(VE)-R^{raw}$ とする。
  3. **蓄積が非零の試験**で BDF1 / BDF2 の符号と単位を検証する。
  4. 周期は **root 単位の 1 回集計**を採用し、部分 CV と壁反力の対応を出力する。
- **過渡では拘束反力の式が違う**: $C_i = D_t(V_iE_i) - R_i^{raw}$ であり、定常の $C=-R^{raw}$ を
  瞬時入熱として使ってはいけない (codex 2 巡目 #1)。dual-time で採るときは**採取位相 (BDF 項を含む/含まない)**
  と単位を固定し、診断側の契約と一致させる。
  それまでに走らせる run は**予備**であり、合格判定には使わない。
- **収支ゲート (規格化まで定義する — codex 2 巡目 #5)**: 同一状態で
  $\varepsilon = \big|\sum_i Q_{f,i} - (\text{固体正味入熱}) \big|$ を、
  **正味量ではなく $\max\!\big(\sum_i|Q_{f,i}|,\; Q_{\rm floor}\big)$ で規格化**する
  (正負が相殺する構成で分母が消えるため)。$Q_{\rm floor}$ は**ケースごとに計算前に登録**する絶対床。
  合格は $\varepsilon/\text{分母} \le 0.5\,\%$ かつ $\varepsilon \le$ 事前登録の絶対許容。
  接合部 (異なる `conjugateGroup` の境目) の授受も明示的に計上し、**共有角の反力は 1 回だけ**数える。
- **面積重みと軸対称**: 界面積分は **host・double**。軸対称は $r$ 重み ([[axisym-rweight-closure-fp32]])。
  **面内伝導作用素にも軸対称の重みが要る** (最終積分に $r$ を掛けるだけでは足りない — codex M6)。

### 4.4 固体モデル

#### 4.4a 抵抗の定義 (codex M3 採用)

未知数は**ガス側表面温度** $T_w$ とする。したがって背面環境 $T_b$ までの**全抵抗**は

$$R_{\rm tot} \;=\; \frac{t}{k_s} \;+\; R_{\rm back},\qquad
R_{\rm back} = \begin{cases} 0 & \text{背面等温 (}T_b\text{ を与える)}\\ 1/h_c & \text{冷却剤}\\ \sum_i t_i/k_i & \text{多層}\\ \infty & \text{断熱}\end{cases}$$

であり、**背面等温でも $t/k_s$ を落としてはいけない** (旧稿は $R_{\rm back}=0$ を「背面等温」とし、
$T_w=T_b$ に退化する誤りだった)。シェル方程式は

$$\nabla_{\!s}\!\cdot\!\left(k_s t\,\nabla_{\!s} T_w\right) \;+\; q_{\rm gas} \;-\; \frac{T_w-T_b}{R_{\rm tot}} \;=\; 0 .$$

- **断熱・孤立系** ($R_{\rm tot}=\infty$ かつ端部断熱) では作用素が定数零空間を持ち、正味入熱が非零なら**定常解が無い**。
  → その構成は**起動時に拒否**するか、適合条件 ($\sum Q=0$) と零空間の固定を明示する。
  **「CG・対称正定値」は Robin 項がある構成に限る**と書く。
- $k_s(T)$・$t(x)$・多層は上の定義に沿って入れる。温度依存は Picard の内側で更新する。

#### 4.4b 幾何と DOF (codex M6 採用)

- **温度 DOF はノード共有、幾何と熱荷重は primal facet 単位**で保持する (`bc.vizBfaceNodes`)。
  node の合成半割面ベクトル ([`gmshReader.hpp`:2250](../../solver_density_cuda/mesh/gmshReader.hpp)) は
  $|\sum_f\mathbf S_f|$ なので**面積として使わない**。
- **四角形は 2 通りの対角線分割を 1/2 ずつ使う** (2026-09-19 実測)。片方の対角線だけで割ると
  集中面積が頂点間で非対称になり、**境界節点に O(1) の荷重不均衡**が残って**全体の収束次数が 2 次 → 1 次**に落ちる
  (帯フィン問題で rate 1.00 → 対称化で 2.00)。
- 定義を明文化して、外れる構成は**起動時に拒否**する: 第一内部点の選び方、面内伝導の非直交補正、
  端部 BC (対称 / 断熱 / 指定温度)、周期境界の同一視、軸対称の面内伝導 ($r$ 重み)、曲率の扱い。
- **`conjugateGroup`**: 連成する壁の集合。グループ内は 1 つの連立。
  **共有 CV に接する「温度を拘束する全ての壁」を走査して競合を起動時に拒否する** (codex 3 巡目 #4):
  - 旧稿は「異なる `conjugateGroup` の共有ノード」しか見ておらず、**相手が素の `wall_isothermal` の場合を捕捉できない**。
    温度ピンは bcond ごとに順に適用され ([`nodeWallDirichlet_d.cu`:154](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu))、
    角ノードは複数 bcond に重複する ([`gmshReader.hpp`:2292](../../solver_density_cuda/mesh/gmshReader.hpp)) ので、
    連成側 800 K と非連成側 300 K が共有すると**最後に適用した方が勝つ**。
  - 同一グループ内は **global CV ID ごとに固体 DOF を一意に対応**させ、
    $Q_j=\sum_{\text{当該 CV の壁面}}F^E-C_j$ を**一度だけ**転送する。
  - 試験は 2 本: **正常な共有角の保存試験**と、**連成壁 × 非連成等温壁の競合が拒否される試験**。
    **どちらも Phase 0 に置く**。

#### 4.4c バックエンド

| mode | 内容 | 用途 |
| --- | --- | --- |
| `local1d` | 点ごとの 1D 抵抗 ($k_s t\to0$) | 面内伝導の寄与を測る基準、薄肉の第一近似 |
| `shell2d` | §4.4a の面内 2D シェル | すきまライナ、冷却ノズル壁 (Phase 2 の生産経路) |
| `fem2d` | **一般 2D 固体領域**の伝導 (三角形 FEM、内部孔は Robin 境界) | **公知データ検証 (§4.9 のタービン翼) に必須**。Phase 1 (外部) のみ |

#### 4.4d `fem2d` の連成契約 (codex 2 巡目 #3 採用)

`fem2d` では未知数が表面温度だけでなく**固体の全節点温度 $u$** になるので、§4.2 の更新式を制限作用素つきで書く。
界面抽出を $E$ (固体節点 → 界面節点)、固体剛性を $K_s$、荷重を $b_s$ として

$$\left(K_s + E^{\mathsf T} D_f E\right)u^{k+1} \;=\; b_s + E^{\mathsf T}\!\left[\,Q_f(Eu^{k}) + D_f\,Eu^{k}\right]$$

- **初版は流体と固体の外周節点を一致させる** (補間を挟まない)。一致しない構成は起動時に拒否。
- $Q_f$ は**積分済みの節点荷重 [W] (平面 2D は W/m)**。$E^{\mathsf T}$ で載せるときに**面積を再乗算しない**。
- 共有角の拘束反力は依存診断側で**接合部の別勘定**になる ([`tooling-energy-balance-diagnostics.md`](tooling-energy-balance-diagnostics.md) §4.1)。
  **どの固体 DOF が受け取るか**を契約に書き、二重計上を禁じる。
- **単体検証** (直列抵抗・フィンだけでは足りない): (i) **孔に Robin 条件を持つ円環伝導の解析解**、
  (ii) **非一様な界面荷重**、(iii) **共有角の保存試験** (角を跨いで熱量が保存すること)。

### 4.5 Phase 0/1: 界面契約の実装と外部弱連成

- **診断出力 (opt-in)**: 壁ダンプに $T_1$, $d_1$, $k_{\rm eff}$, `q_compact`, `q_recon`, `q_eff` (拘束反力込み) を足す。
  既定 run のコストを上げない。**既存 run には無いので、監査は新規 run で行う** (codex M9)。
- **`wallProfile`**: 既存 `applyInletProfiles` は種別非依存なので**フラグ名の追加と座標評価の修正が主**。
  - **node はノード座標、cell は面重心**で補間する (codex M5。実測ずれ 0.679 mm)。
  - CHT 内部の転送は**安定なノード ID を正本**にし、座標補間は外部 CSV 入力のときだけ使う。
  - **verify は `bvar` の再出力では不十分**。非一様プロファイル試験では `VALUE/T` (場の温度) と EOS 整合まで見る。
- **ループ** `tools/cht_loop.py`: warm restart (同一メッシュ index コピー) → 界面量抽出 → 固体求解 → §4.2 の更新 → 次 run。
  収束判定は §6 の界面ゲート。**旧 run を根拠にしない**。
- **コストの正直な記述** (codex m11): 外部ループで消えるのは**起動・I/O・スケジューリング**であり、
  **流体の熱場緩和 (case/44 で 36000 step 後も $Q_w$ が 6 % 動く) は Phase 2 でも残る**。
  Phase 2 の速度向上は主張でなく**測定項目**にする。

### 4.6 Phase 2: ソルバ内連成 (node 限定)

- **種別は `wall_isothermal` のまま**、`ints: {conjugate: 1}` で連成を有効化する (codex M4)。
  これで `iso_wall_flag` ([`mesh.cpp`:929](../../solver_density_cuda/mesh/mesh.cpp))・T ピン・粘性壁・壁距離・
  DPLUR のエネルギー行切離し ([`timeIntegration_d.cu`:1025](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)) の
  既存経路にそのまま乗る。**新種別は作らない**。
- **フック位置は行番号でなく関数・評価位相で決める** (codex 2 巡目 #4)。
  旧稿が指した `main.cpp`:1684 付近は **`advanceExplicitRK` ([:1689](../../solver_density_cuda/main.cpp)) のステージループ**で、
  **主対象の定常陰解法は `advanceImplicitSteady` ([:1750](../../solver_density_cuda/main.cpp)) → `implicitNonlinearUpdate`
  ([:1596](../../solver_density_cuda/main.cpp)) を通る**ため、旧稿どおりに実装すると陰解法で連成が作動しない。
  - **定常陰解法**: 完了したステップ数で $K$ を数え、**同一状態で採った $(Q_f, T_w)$ の対**を固体へ渡し、
    **次の残差組立ての前**に新しい `Ts` を適用する。
  - **dual-time** ([`advanceImplicitDualTime`](../../solver_density_cuda/main.cpp):1769) は**別契約**とし、
    subiteration の途中で壁温を変えない (物理時間ステップの境界でのみ更新)。初版では dual-time 連成は対象外。
- **暖機**: $N_{\rm warm}$ step は等温固定。段階起動の**本段でのみ**連成を有効化する。
- **restart / 区間**: 固体状態を `conjugate_state.h5` に保存。`stage_manifest.py` の区間識別に
  **連成の有効/無効・$K$・$D_f$ 方式・`solid.json` と `wallProfile` の内容ハッシュ**を含める (codex M10)。
- **実装保証は node のみ**。cell で `conjugate: 1` を指定したら起動時に拒否する。

### 4.7 (欠番 — 旧 §4.7 は §4.8 に統合)

### 4.8 Phase 3 (固体ゾーン) に進む条件

Phase 2 の実測で次のいずれかが示されたとき、**別 plan** を起票する。

1. `local1d` と `shell2d` の差が目的量の許容を超え、かつ厚さ方向 1D 近似が破れる ($\mathrm{Bi}_t=ht/k_s>0.1$)。
2. 角部・リップで固体内 2D/3D 熱橋が支配的で、シェル近似誤差が V6 の格子・モデル感度より大きい。
3. 厚肉で背面条件を 1 点に縮約できない (§4.9 のタービン翼はここに該当するので、**Phase 1 の `fem2d` で扱う**)。

### 4.9 一次適用先 — 公知試験データ (ユーザ指示 2026-09-19)

**指示**: 「一次の適用先は**公知文献の試験結果があるもの**」「**超音速系だとなお良い**」。
調査は [`notes/investigations/cht-validation-case-survey.md`](../../notes/investigations/cht-validation-case-survey.md)。

- **決定**: **NASA CR-168015 (Hylton et al., 1983; NTRS 19830020105) のタービン翼**を一次適用先とする。
  - **主 = Mark II 翼 case 5411 = run 42、出口 $M_2$=1.04** (超音速出口。衝撃 × 境界層 × 熱伝達)。
  - **従 = C3X 翼 case 4411 = run 108、出口 $M_2$=0.90** (亜音速出口)。先に通す足場。
    **run 番号と $M_2$ は codex 2 巡目が報告の表 VIII・IX で確認**したもの (本計画でも一次資料で再確認する)。
  - **報告自身のデータ整理が 2D 伝導 FEM を使っている**ので、`fem2d` バックエンドの選択は原典と整合する。
  - 選定理由: 固体条件 (材料・$k_s(T)$・冷却孔) が公開され、**壁温が結果として測られている**。
    極超音速側には「$T_w$ を解として測った」公知試験がほとんど無い (数値研究が主) — 調査メモ §2。
- **注意 (検証設計に効く)**:
  1. 固体は薄肉シェルでなく**冷却孔を持つ厚い 2D 固体** → **`fem2d` バックエンドが必須**、Phase 1 で実施。
  2. **入力不確かさは HTC だけではない** (codex 2 巡目 #6): 報告 p.21 は内部 HTC を**入口助走補正つき相関**で与え、
     **計測断面の冷却剤温度は入口・出口の実測からの推定**である。さらに表 VII の**温度比不確かさは ±2 %** で
     HTC の不確かさとは別項目。→ 帯は (冷却孔ごとの HTC 相関) + (冷却剤温度推定) + (材料物性) + (±2 %) を
     **別々に積んで計算前に登録**し、**合うまで広げない**。
  3. 翼列なので**周期境界**が要る (node の周期は修正履歴あり)。超音速出口は衝撃があるので乱流モデル依存が大きい。
- **役割分担と、その限界** (codex 2 巡目 #6): CHT の誤差 = 流体側の $q_w$ 誤差 + 連成側の誤差。
  ただし**平板・すきまの流体側検証は、翼の圧力勾配・遷移・衝撃による流体側誤差を拘束しない**。
  → 同じ翼・同じ試験条件の中で **(a) 実測壁温を与えた流体計算 → (b) 公開条件での固体単独検証 → (c) 壁温未知の CHT**
  の順に分解する (§6 V5)。**壁温が合っても、流体 HTC と内部冷却条件の誤差相殺は排除できない**ので、
  (a)(b) を通さずに (c) の一致を成果として報告しない。
  case/49・case/51 への誤差幅は**その用途のモデル感度から別途**評価する (翼の結果を直接移植しない)。

### 4.10 既知の限界 (合否の外。報告には必ず書く)

- **表面間放射**: $\sigma T^4$ は **500 K で 3.5 / 1000 K で 57 / 1500 K で 287 kW/m²**
  (旧稿は 1500 K を 57 と誤記 — codex 2 巡目が検出)。深いキャビティでは無視できない。解かないが**桁を併記**する。
- **非定常 CHT**: 固体の熱時定数は流体の $10^3$〜$10^6$ 倍。準定常連成はコストが別次元なので扱わない。
- **float32**: 1000 K で 1 ULP $\approx6.1\times10^{-5}$ K。**host を double にしても失われた温度差は戻らない** (codex M10)。
  ゲートは「必要な熱流束を温度量子化で解像できるか」を含める。
- **接触熱抵抗**: 入力できるが同定はしない。

## 5. 実装ステップ

**docs を先に更新する** (codex m12、[`AGENTS.md`](../../AGENTS.md) 開発フロー)。

1. ~~**S0 docs 先行 + 索引の同期**~~ — **完了 (2026-09-19)**。
   [`methods/boundary.md`](../../methods/boundary.md) に「共役熱伝達 (CHT)」節を追加
   (対応範囲 / 界面量の定義と符号 / 固体モデルと `fem2d` 契約 / 反復と受理判定 / `wallProfile` /
   **起動時に拒否する構成 4 件** / 診断出力とゲート)。冒頭に「仕様確定・未実装」を明示。
   [`methods/index.md`](../../methods/index.md) の境界条件行に状態を追記。
   [`plans/README.md`](../README.md) と親 plan の旧方針も同期済み。以下は元の記述: [`methods/boundary.md`](../../methods/boundary.md) に「共役熱伝達」節 —
   界面契約 (§4.3)・固体式 (§4.4)・対応範囲 (node 限定・拒否条件) を**実装前に**書き、
   [`methods/index.md`](../../methods/index.md) を同期。
   **撤回済みの旧方針が索引・親 plan に残っているので同時に直す** (codex 2 巡目 #7):
   [`plans/README.md`](../README.md) の「コンパクト差分形を一次」「Bi 非依存の安定化」、
   [親 plan](tooling-nozzle-isothermal-wall-chain.md) §5.1 の `wall_conjugate`。
2. ~~**S1 界面診断**~~ — **完了 (2026-09-19)**。`conjugateWall.{hpp,cpp}` 新設 (符号規約の正本)、
   `output: {interfaceDiag: 1}` で壁ダンプに `iface_T1/d1/keff/q_compact/q_recon/q_2nd/ok/align`。
   host 側で作るので device `bvar` を汚さない。`q_eff` (拘束反力込み) は
   [`tooling-energy-balance-diagnostics`](tooling-energy-balance-diagnostics.md) の完成後。
3. ~~**S2 `wallProfile`**~~ — **完了 (2026-09-19)**。`applyInletProfiles` を `applyBoundaryProfiles` に
   一般化し `applyWallProfiles` を追加。**壁は値を課す位置 (node はノード座標) で補間**、入口は従来どおり
   face 重心 (既存 run のビット不変)。verify は場の `VALUE/T` で確認済み。
4. **S3 幾何と接合** — **部分完了 (2026-09-19)**: **共有 CV の壁温競合の起動時エラーを実装**
   (`checkWallTemperatureSharing`。壁温を陽に扱う run でのみ走り、既定 run の挙動は不変)。
   **残り**: facet 単位の幾何・荷重、`conjugateGroup`、軸対称・周期の規定と未対応構成の拒否。
5. ~~**S4 固体ソルバ**~~ — **完了 (2026-09-19/20)**: `tools/solid_shell.py` (`local1d`/`shell2d`) と
   `tools/solid_fem2d.py` (一般 2D + Robin 孔 + Schur 縮約)。単体は解析解 (`test_solid_shell.py` /
   `test_solid_fem2d.py` とも **PASS (all)**)。
6. ~~**S5 外部ループ**~~ — **完了 (2026-09-19)**: `tools/cht_loop.py` (1 反復 = forge 1 回。
   §4.2 の更新式 + Anderson + 受理/退避、`wallProfile` の書き出し、warm start、`cht_history.csv`)。
6.5. **S5b 界面ゲートの実装**: G-if / G-cons の判定 (局所ノルム・絶対許容・連続反復数) と区間ハッシュ。
   **V1 の判定より前に用意する** (codex 2 巡目 #5。旧稿は S6 に置いていた)。
7. **V1 / V2 / V3** (§6)。**ここを通るまで Phase 2 に進まない**。
   V1 は**依存診断の dual-time 解除**が前提 (§4.3)。
8. **S6 ソルバ内連成** — **Phase 2a 完了 (2026-09-19)**: `conjugate:` ブロック + `ints: {conjugate: 1}`、
   **`local1d` (点ごとの 1 次元抵抗、行列不要)** を C++ 実装。更新は抵抗加重平均で、`interval` step ごとに
   **ステップ完了後** (`advanceOneStep` 末尾 = 次の残差組立ての前) に適用。起動時拒否 (node 以外 / dual-time /
   `mode != local1d` / 背面断熱 / 第一内部点が定まらない壁)。再開は `conjugate_Tw_<physID>.csv` → `wallProfile`。
   **残り**: ソルバ内の面内伝導 (`shell2d`)、`stage_manifest` の区間 (外部入力ハッシュ)、界面ゲートの出力。
9. **V4 / V5 / V6**。
10. **S7 docs 整合確認**: 実装と S0 の差分、`procedures/su2-cross-check.md` に CHT (multizone) 手順、
    `procedures/recommended-settings.md`、`procedures/verification/`、`design/CAPABILITIES.md`。

### 5.1 残作業 (優先順)

**順序 (codex 推奨)**: ①実効壁熱量と符号 → ②固体抵抗と固定点保存の更新式 → ③等温壁経路・座標・接合の配管 → ④V1/V2 とゲート。

**S0 完了 (2026-09-19)**: #20 (受理判定)・#22 (ゲートの式)・#23 (共有角の拒否)・#24 (`fem2d` 契約)・
#26 (陰解法フックの対応範囲) の**仕様は [`methods/boundary.md`](../../methods/boundary.md) に反映済み**。
#21 は提供側 [`tooling-energy-balance-diagnostics`](tooling-energy-balance-diagnostics.md) §5.1 #7・#8 に登録済み。
#27 は完了。**以下の表に残るのは実装・検証の作業**である。

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | **codex 2 巡目レビュー** | 1 巡目 **NO-GO** (C1/M9/m2) を全件採用して改訂済み。**再レビューを通してから実装に入る** (§6.1) |
| 2 | 界面契約の確定 (C1) | §4.3。$Q_f$ の定義・符号・帰属・接合部の計上。`methods/boundary.md` に先に書く (S0) |
| 3 | 依存: 拘束反力の採取 | [`tooling-energy-balance-diagnostics`](tooling-energy-balance-diagnostics.md) の完了を **Phase 1 の合格条件**として登録 |
| 4 | 界面量の実測差 (M9) | **新規 run** で コンパクト / 2 次片側 / 再構成勾配 / 実効 を比較。旧 run (`run_0011`, `run_0114`) は予備資料。既知: コンパクト vs 2 次片側 = 2.67 % |
| 5 | 更新式 (M2) | §4.2。$(A_s+D_f)T^{k+1}=b_s+Q_f(T^k)+D_fT^k$、セカント $D_f$ の上下限、Aitken の 3 つの安全装置 |
| 6 | 固体抵抗 (M3) | §4.4a。$R_{\rm tot}=t/k_s+R_{\rm back}$、断熱孤立系の零空間処理、SPD 成立条件の明記 |
| 7 | 配管 (M4) | §4.6。`wall_isothermal` + `conjugate` 属性。`iso_wall_flag`・T ピン・粘性壁・壁距離・DPLUR 切離しの 5 経路を確認。cell は拒否 |
| 8 | `wallProfile` 座標 (M5) | §4.5。node はノード座標。verify は `VALUE/T` まで |
| 9 | 幾何・接合 (M6) | §4.4b。facet 単位の幾何、`conjugateGroup`、**角ノード試験を Phase 0 に前倒し**、軸対称面内伝導、周期 |
| 10 | 固体ソルバ | ~~`local1d` / `shell2d`~~ **実装・検証済み (2026-09-19)**: `tools/solid_shell.py` + `tools/test_solid_shell.py` (直列抵抗 機械精度 5.7e-14 K / フィン **格子収束 rate 2.00** / 零空間の拒否 / $k_s(T)$ / 連成反復)。**残り**: `fem2d` (一般 2D 領域) |
| 11 | V1 (M7) | §6。**1 次元純伝導 (静止流体層 + 固体層)** に固定。厚さから $R_f$ を厳密定義。dual-time で定常化 |
| 12 | V2 (M8) | §6。SU2 CHT は multizone。**まず 1D スラブで成立**させ、版・固体メッシュ・界面対応・物性・両ゾーン残差を固定。`procedures/su2-cross-check.md` に手順追加 |
| 13 | 界面ゲート (M10) | §6。既存ツールの拡張: 未緩和の局所界面不釣合い・固体方程式残差・**絶対**温度更新量・温度量子化の解像性。区間識別に外部入力のハッシュ |
| 14 | ~~外部ループ~~ **完了 (2026-09-19)**: `tools/cht_loop.py` + `case/52.conjugate_slab` (V1 PASS) | 残り: `--flux q_eff` (拘束反力込み) を依存診断の完成後に接続 |
| 15 | ソルバ内連成 | §4.6。~~`local1d`~~ **Phase 2a 完了 (2026-09-19)**。残り: `shell2d` の C++ 化、`stage_manifest` 区間、界面ゲートの出力、dual-time 契約 |
| 16 | **一次適用先 (公知データ)** | §4.9。**着手 (2026-09-20)**: CR-168015 入手 (`papers/cht/`, 追跡外) → **翼型座標と試験条件を抽出・検証済み** (`case/53.c3x_vane_cht/tools/extract_vane_data.py`: cm/in の 2.54 則で全点検証、OCR が壊した 11 点は画像目視で補修、欠測があれば落ちる)。C3X 78 点 / Mark II 60 点、run 108 (M2 0.90) と run 42 (M2 1.04)。**報告自身の不整合 2 件を記録** (表 VIII/IX の SI 圧力列が psia と 51.7 倍ずれ / C3X 点 29 の inch 誤植)。**残り**: $k_s(T)$・冷却孔配置・孔ごとの HTC と冷却剤温度・測定壁温分布・カスケード幾何 (表 IV) → メッシュ |
| 17 | Phase 3 判定 | §4.8 の 3 条件を V6 の数値で評価し、要否を結論づける |
| 18 | docs 同期 | S0 / S7 |
| 20 | **受理・退避仕様 (停止しないこと)** (3 巡目 #1) | §4.2。メリット関数 $\Phi=r^{\mathsf T}(A_s+D_f)^{-1}r$ + line search、再試行上限、停滞の不合格化。V1b に反例と非一様 $D_f$ |
| 21 | **依存診断の解除契約** (2 巡目 #1 + 3 巡目 #2) | dual-time 解除を **V1 前**、周期解除を **V5 前**。**提供側 plan の残作業表にも** $R^{raw}$ の定義・$D_t(VE)$ の算出・BDF 符号試験・周期 root 集計を登録する |
| 22 | **ゲートの数値仕様** (2 巡目 #5 + 3 巡目 #3) | **一部実装 (2026-09-19)**: `tools/cht_wall_series.py` が壁ダンプから界面量の時系列 (`Tw_*`, `q_total`, `imbalance_max/rel`) を書き、`check_quasisteady.py --series-csv` に渡せる。**残り**: 絶対/相対/連続回数の事前登録と自動判定、区間ハッシュ。**実測の注意**: `imbalance_rel` は 1e-4 級の微小量なので相対 drift 判定は意味を持たず `DRIFTING` になる → **絶対値 (W/m²) で報告する**。| §6 G-cons/G-if。規格化 ($\sum|Q|$ と絶対床)、局所面積尺度、絶対×相対×連続回数、準定常は **drift と osc の両方**を比較許容の 1/5。**Phase 1 判定前に実装** |
| 23 | **共有角の唯一の所有者** (3 巡目 #4) | §4.4b。~~温度を拘束する全壁を走査して拒否~~ **実装・検証済み (2026-09-19)**: `checkWallTemperatureSharing` が競合 CV・physID・各 $T_w$ を出して exit 1 (case/48 で `sym` を 500 K 等温壁にした拒否試験)。**残り**: global CV ID で固体 DOF を一意化、$Q_j$ の 1 回転送 (連成実装時) |
| 24 | **`fem2d` の連成契約と単体試験** (2 巡目 #3) | ~~§4.4d~~ **実装・検証済み (2026-09-20)**: `tools/solid_fem2d.py` (線形三角形 FE + Robin 辺 + **内部節点を消去した Schur 補元**で界面作用素をシェルと同じ顔にする) + `tools/test_solid_fem2d.py` **PASS (all)**: 孔 Robin の円環解析解 (**rate 1.97**)、Schur = 全系解 (5e-12)、非一様荷重、連成 **7 反復**で厳密固定点。~~`cht_loop` からの選択~~ **完了 (2026-09-20)**: `--solid-mode fem2d` + 界面を座標一致で 1 対 1 に強制 (内挿しない)。**残り**: 共有角の保存試験 (界面が複数 bcond に分かれる構成) |
| 25 | **V5 各段の入力と合格条件** (2 巡目 #6 + 3 巡目 #5) | §4.9・§6 V5。(a) 実測 $T_w$ → $h$/熱流束/壁圧、(b) ~~外周 Dirichlet + 孔 Robin の固体単独~~ → **設計変更 (2026-09-20)**: **冷却孔ごとの冷却剤温度・流量は報告に無い** (方法のみ記載、付録 A にも欄が無い)。(b) は**公開量から内部条件を逆算する**段に変更し、同定の不確かさを帯に別項目で積む。(c) CHT。帯の**合成規則**まで事前登録 |
| 26 | **陰解法フックの設計** (2 巡目 #4) | §4.6。`advanceImplicitSteady`/`implicitNonlinearUpdate` 経路で $K$ を数える。dual-time は別契約 (初版は対象外) |
| 27 | **索引・親 plan の残り** (3 巡目 #6) | 指定箇所は同期済み。**親 plan §4.6-4 の「弱ループが収束しない = 軸方向伝導が支配的」**を本計画 §4.8 のモデル感度基準に置き換える |
| 28 | ~~**C3X 翼列の SU2 対照**~~ **完了 (2026-09-20)**: 同一メッシュ・同一 BC・一様壁温 566 K で SU2 8.5 低 Re SST と比較。負圧面層流域は SU2 も **+45.8 %** (forge +42.7 %) で **forge vs SU2 は −3.0 %/6.0 %**、遷移後は **−0.6 %/2.1 %**。**層流域の過大は forge 固有でなく低 Re SST 共通**と確定。残る solver 間差 (正圧面 −6.0 %/7.1 %) は別件。| 同一メッシュ・同一の実測壁温分布で SU2 低 Re SST を回し、$h$ の分布を突き合わせる。**V5 段 (a) の差が forge 固有かモデル共通かを切り分ける**ため ([[reichardt-5pct-gap-not-forge]] と同じ性格の差が出る可能性)。手順は [`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md) |
| 29 | **遷移モデル** (2026-09-20, V5 の差の主因。**Mark II で必須と確定**) | 低 Re SST に $\gamma$–$Re_\theta$ 等が無いため負圧面前縁の 層流域で $h$ が +40 % になり、CHT 壁温に +21 K 効く。**CHT の合否とは分けて扱う** (帯の別項目) が、Mark II (超音速出口) では遷移位置がさらに効くので、少なくとも**遷移位置を与えた感度計算**を V6 前に行う。実測 (2026-09-20): C3X 負圧面層流域 $h$ +40.5 %/$T_w$ +21.4 K、**Mark II は +74.8 %/+42.8 K** |
| 30 | **C1 保存的界面熱量の接続** | **定常 node 等温壁について完了 (2026-09-20)**: `iface_q_eff` = $(R^{raw}-F_w)/A$ を壁ダンプに出し、`cht_loop --flux q_eff` が使えるようになった。素材は `ifaceFw` (壁半割面が `res_roe` に入れた寄与) と `ifaceRraw` (**壁残差射影の直前**の `res_roe`。後から読んでも 0 しか出ない)。**V1 で検算済み** (解析解 81.3090 W/m² に対し 81.1837 = −0.154 %、`q_compact` −0.155 % と同等 = 符号・絶対値とも再現)。**効く**: case/53 C3X 480 節点で `q_compact` との差が bias +2.26 % / rms 3.28 % / 局所最大 19.1 %、連成壁温で +11.5 → +13.5 K。仕様は [`methods/boundary.md`](../../methods/boundary.md)。**残り**: (a) dual-time ($C=D_t(VE)-R^{raw}$)、(b) 周期壁ノード (root 単位集計)、(c) 軸対称の $r$ 重み、(d) G-cons ゲート。いずれも該当構成では `NaN` を出す (case/53 は壁ノードが周期に属さないので 480/480 有限) |
| 30b | **G-cons (収支ゲート)** | **完了 (2026-09-20)**: `tools/check_cht_balance.py`。同一状態で $\varepsilon=\lvert\sum Q_f-Q_{\rm solid}\rvert$ を **$\max(\sum\lvert Q_f\rvert, Q_{\rm floor})$ で規格化** (正味量で割らない)、`--q-floor` は**ケースごとに事前登録**。`iface_q_eff` が `NaN` の節点が 1 つでもあれば不合格。実測: `run_0022_cht_qeff` (1 反復 4000 step) は **FAIL 2.77 %** で G-if の `res_rel` 3.5 % と整合。16000 step では **PASS 0.0272 %**。**同一状態で定義だけ変えると `q_compact` 2.66 % / `q_recon` 2.66 % が FAIL、`q_2nd` 0.29 % / `q_eff` 0.027 % が PASS** = C1 の実証。V1 (`local1d`, ソルバ内連成) は **0.000125 %** で PASS。`local1d` 対応済み (`--solid-mode local1d`) |
| 31 | **M2 受理・退避の状態一貫性** | 棄却時に $T$・$Q_f$・固体状態を**同じ評価点の組**で保存/復元する。$D_f$ を変えたら基準メリットを再計算し履歴を捨てる。初期壁温は**実際に課した分布**から取る (実 run では課 512–612 K / 仮定 566 K 一様だった)。line search と再試行上限も未実装 |
| 32 | **M3 $D_f$ の安定条件** | `hA` は初期推定に留め、「安全率 2 で十分」の主張は撤回。安定性は #31 の受理処理で担保する。試験は低周波だけでなく**交番温度摂動**と CFD 緩和長依存 |
| 33 | **M4 `fem2d` の局所 $k_s(T)$** | `FixedPointDriver.advance()` が `recover_interior()` を呼ばず、全域が $k_s(\overline{T_w})$ になっている。各評価温度で内部温度と物性を自己整合させ、内部残差も判定する。温度依存円環で `driver` と全系求解を照合 (codex 実測 494.82 vs 492.75 K) |
| 34 | **M5 収束ゲート (G-if) の実装** | **完了 (2026-09-20)**: 規格化を $\max\lvert Q_f\rvert$ のみに直し (反例は `test_solid_shell.py` **T7**: 旧 `res_rel` 3.33e-6 に対し物理的不釣合い 100 %)、**絶対残差 `res_abs` [W]・相対残差 `res_rel`・固体内部残差 `res_solid`・温度更新 `dTw`・退避していないこと**を独立に満たし `--n-consec` 回連続したときだけ収束とする。許容は `cht_loop --tol-abs-W / --tol-rel / --tol-solid / --tol-K` で**事前登録**。`res_solid` は `Fem2DOperator.interior_residual`。`cht_history.csv` に 3 列追加。仕様は [`methods/boundary.md`](../../methods/boundary.md) |
| 35 | **M6 V5 の格下げと独立検証** | V5 は「**同定条件下での整合性評価**」。$h_c$ は公開されていない (公開は $T_c$ と流量のみ) ので、孔位置の再構成と相関の選択が $h_c$ に入る。**「連成は正しい」「差の主因は遷移モデル」を結論として書かない** (遷移は有力仮説)。切り分けは #28 SU2 対照 + #29 遷移感度 + 格子感度 |
| 36 | **M7 局所量の準定常判定** | 報告する**局所量**にも事前登録した許容を当てる。積分 `q_total` の `STEADY` を局所分布の保証に使わない |
| 37 | **M8 共有角の所有者検査** | 初期温度の一致でなく**共有 CV の拘束所有者**を検査する。`conjugateGroup` と一意温度 DOF、未対応構成の明示拒否 |
| 38 | **M9 同一メッシュ再開** | `cht_loop` の warm start が `interp_field` (2D 最近傍)。同一性を確認して index コピーにする。3D spanwise 非一様場で再開不変を検証 |
| 39 | **M10 cell の `wallProfile` 座標** | `boundaryCond.cpp` が cell で内部セル重心を使っている。面重心 (`msh.planes[ip].centCoords`) に戻す |
| 40 | **m11 docs と残作業表の同期** | `methods/index.md` / `design/CAPABILITIES.md` が「未実装」のまま。#10/#24 の記述矛盾、`stage_manifest` の外部入力ハッシュ未実装 |
| 41 | **forge と SU2 の残差 (層流壁熱伝達)** | 遷移後 −1.6 %/2.6 % に対し **正圧面 −6.7 %/7.6 %・層流域 −4.0 %/6.4 %**。**消した可能性** (2026-09-20, case/53 README「残った forge–SU2 差の切り分け」): 熱流束の定義 (SU2 も同じコンパクト差分に揃えた)・幾何 ($d_1$ 同一)・層流物性 (Sutherland に対し 0.14 % 以内)・第一内点の $\mu_t$ (両者 0)・`sstEnergyIncludesK` (SU2 形にしても改善せず、**棄却**)。**残る事実**: 速度は ±1 % 一致なのに $T(y)$ だけ壁から 20 µm でずれる。層流単独対照は forge が非定常で**不成立**。次の候補: 対流スキーム (SLAU+MUSCL vs ROE+Venkat) の近壁温度への効き / node 壁半 CV のエネルギー流束 / 前縁からの BL 発達履歴 |
| 42 | ~~**`outlet_statPress` が擬似 CFL を律速する**~~ → **撤回 (2026-09-20)**。`cfl_pseudo` 4.0 の発散は事実だが **(a) 原因は出口 BC でなく 2 次再構成** (1 次なら安定、リミッタを Barth/`venkatK` 0.01/旧経路に替えても全滅、ゴースト緩和 w=0.3 も無効)、**(b) `implicitRelax` 0.7 や `nStepInner` 30 で上限は外せるが速くならない** — SST 本番で 20000 step の最終 `rms_ro` は現行 `cfl_pseudo` 0.5 が **4.6e−06**、cfl 2.0 が 4.5e−03、cfl 8+relax 0.7 が 1.6e−03、cfl 4+inner 30 が 7.0e−03 (所要時間はほぼ同じ)。**現行設定が最良**で、生産 run の律速ではない。詳細は case/53 README |
| 19 | codex result レビュー | `done` にする前 |

## 6. 検証

- **固体ソルバの単体 (完了, 2026-09-19)**: `python3 solver_density_cuda/tools/test_solid_shell.py` → **PASS (all)**。
  T1 直列抵抗 = 解析解に機械精度 (5.7e-14 K)、T2 フィン = **格子収束 rate 2.00** (4 段)、
  T3 断熱孤立系の拒否、T4 $k_s(T)$、T5 連成反復が非対角応答でも厳密固定点に到達 (57 反復, 5.6e-8 K)、
  **T5b = codex 3 巡目 #1 の反例を回帰試験化** (最大ノルムは 1.000→1.701 と増えるが $\Phi$ は 1.0081→0.8831 と減る)、
  T5c = 素の固定点反復では収束しないこと (加速器が飾りでない根拠)。
- **V4 (Phase1 ↔ Phase2) の実測 (2026-09-19)** — run: `case/52.conjugate_slab/run_0003_v1_insolver_cont/`。
  **ソルバ内連成 (`local1d`)** (`check_convergence.py`: `run_0002` は **`PASS (converged)`**、
  継続の `run_0003` は `NOT CONVERGED (stalled/plateau)` — 静止純伝導で残差が初期から 1e-10 のため
  低下桁数が取れないことによる。合否は解析解照合・界面静定・流体緩和の 3 点で見る: case/52 README):
  $T_w$=**316.2371 K** (解析解比 **−0.152 % of rise**)、
  **両側 $q$ の不一致 0.0002 %** (81.1856 / 81.1855 W/m²)、更新量 1.5e-5 K で静定。
  外部ループ (shell2d, 316.2659 K) との差は **0.029 K = 0.18 % of rise** で V4 の許容 (≤0.3 %) 内。
  差の出どころは**同じ壁温での流体側の離散解の差** (固定壁温の緩和試験で $q$ が ±0.2 % 動く) であり、
  連成の定式化の差ではない (両者とも自分の界面条件を 0.005 % 以下で満たしている)。
  再開経路 (`conjugate_Tw_<physID>.csv` → `wall_profile_<physID>.csv`) も同 run で通した。
- **外部ループを乱流平板に当てた結果 (2026-09-19, 未収束)** — run: `case/48.flat_plate_cooled_m4/run_0017_cht_shell/`。
  12 反復 (各 6000 step) で `res_rel` が 0.28–0.31 から下がらず、$T_w$ 最大値が 818–863 K で振動して**未収束**。
  前縁近傍で $h$ が桁で変わるため $D_f^{(0)}=k_{\rm eff}A/d_1$ が過大 (平均 68.8 W/K) で過減衰になったことと、
  各反復の CFD が 6000 step では収束しきっていないことが効いている。
  **「外部ループは検証用、生産はソルバ内」という §4.5 の位置づけを実測で裏づけた形**。
  **同じ問題をソルバ内連成は静定させられた** (`run_0018_cht_insolver` 30000 step → `run_0019_cht_insolver_cont` +60000 step、
  合計 6.7 分)。**`check_quasisteady.py --series-csv wall_series.csv` → `STEADY` (ALL STEADY、漸近値 = 最終値 +0.000 %)**:
  $T_w$ 平均 **578.30 K** / x=0.5 m **569.87 K** / 範囲 **551.63–986.15 K**、$Q_w$ **60.39 kW/m**、
  **両側 $q$ の不一致 局所最大 31.6 W/m² = 相対 2.1e-4**、更新量 1.6e-3 K。
  ただし **`check_convergence.py` は `NOT CONVERGED (stalled/plateau)`** で、**`--from-floor run_0011` は REFUSED**
  (参照 run_0011 自体が PASS しない = この case には収束した参照床が無い)。
  → **case/48 系列は残差ベースでは合格にできない**ので、派生量の定常性と界面不釣合いで判定する。
  界面量の時系列は `tools/cht_wall_series.py` が壁ダンプから作る (G-if の入力)。
  以下は 30000 step 時点の途中値の記録:
  (`check_convergence.py` VERDICT は **`NOT CONVERGED (stalled/plateau)`** — `rms_roe` が 2.1 dec で flat、
  `rms_ro`/`rms_roK`/`rms_roOmega` は falling。**親の run_0011 が NOT CONVERGED なのを継承**しており、
  以下は収束解ではなく 30000 step 時点の値):
  $T_w$ = 548 K (x=1.0) … 970 K (前縁)、更新量 0.043 → 0.028 K/更新で単調減衰、
  **両側 $q$ の不一致 中央値 0.13 % / 最大 0.60 %**。`q_2nd` と `q_compact` の差は 1.59 %
  (§4.3 の「どの定義か明示する」が実ケースでも効く桁)。
- **V1 合格 (2026-09-19)** — run: `case/52.conjugate_slab/run_0001_v1_slab/` (11 反復 × 10000 step)。
  `python3 case/52.conjugate_slab/verify_v1.py run_0001_v1_slab` → **VERDICT PASS**:
  $T_w$ = **316.2659 K** (解析解 316.2618、温度上昇 16.27 K に対し **0.025 %** ≤ 0.5 %)、
  **両側 $q$ の不一致 0.0053 %** (81.3252 / 81.3296 W/m² ≤ 0.1 %)、$q$ の解析解差 0.0199 %。
  外部ループは `dTw_max` 7.2e-5 K・`res_rel` 7.2e-6 が 2 反復連続で許容以下。
  **この case で分かった運用上の要件** (README に記載):
  (a) **壁クラスタメッシュでは擬似時間の熱緩和が終わらない** (d₁=8.5e-5 で 40000 step でも $q$ が定常値の 6 倍) →
  一様メッシュにする。(b) **低圧にして $k/(\rho c_p)$ を上げる**と緩和 step 数がほぼ比例して減る
  (101325 Pa → 1013.25 Pa で 20000 step 未収束 → 10000 step で解析解 −0.16 %)。解析解は圧力に依らない。
  (c) **case/24 の旧設定 `cfl_pseudo: 100` は現 HEAD では step 100 以内に NaN** (P/T 床に張り付く EOS 床洗浄。
  1 次・リミッタ無しでも同じなのでスキーム依存ではない)。本 case は 5 で安定。**これは本計画の変更とは無関係の既存挙動**。
- **S1–S3 のコード検証 (完了, 2026-09-19)** — run は `case/48.flat_plate_cooled_m4/`:
  - `run_0013_iface_base` (基準, interfaceDiag 0) / `run_0014_iface_diag` (1) /
    `run_0015_wallprofile` (+ `wallProfile`) / `run_0016_iface_base_rep` (反復=ノイズ床)。各 200 step warm start。
  - **診断は解を動かさない**: 場の相対差 base↔diag 1.9e-6 に対し、**同一設定の反復間 (ノイズ床) が 1.8e-6**
    (atomicAdd の非決定性)。
  - **第一内部点**: $d_1$=3.0001 µm = 第一層厚と一致、align 1.000、1001/1001 点評価可。
  - **3 形式の差**: `q_recon` は `q_compact` と **1.6e-7 相対で一致** (壁法線に整列した node メッシュでは
    再構成勾配がコンパクト差分に帰着する)。**`q_2nd` (2 次片側) は中央値 2.79 %・最大 3.1 % 違う**
    → codex が `run_0011` で測った 2.67 % は**後処理の差分形式の差**であり、ソルバ内部の不整合ではない。
    滑らかな分布では 2 次片側が正確なので、**カーネルの壁熱流束はこの解像度で ~3 % の 1 次打ち切り誤差**を持つ
    (§4.10 の誤差予算に入れる)。**斜交メッシュでの 3 形式の分離は未測定** (残作業)。
  - **`wallProfile`**: $T_w=300+100x$ が**場の `VALUE/T`** に入る (最大差 0.068 K)。
  - **共有角の拒否**: `sym` を 500 K 等温壁にすると前縁の 1 CV を検出して exit 1。
  - これらは**収束を主張する run ではない** (200 step のコード検証)。NaN/Inf 無しは確認済み。
- **単体 / ビルド**: 固体ソルバの単体 (1D 直列抵抗・面内フィン)。**格子収束で連続解析解に近づくこと**と
  **離散残差が機械精度で 0 になること**を**別の試験**にする (codex M7)。`wallProfile` の回帰 (ビット不変)。
  角ノード競合の起動時エラー試験 (L 字の 2 bcond)。
- **検証ケース**: 新規 `case/52.conjugate_slab/` (V1/V2)、`case/48.flat_plate_cooled_m4/` (V3)、
  新規 `case/53.c3x_vane_cht/`・`case/54.markii_vane_cht/` (V5)、`case/51` (V6、親 plan が作る)。

| # | 検証 | 内容 | 合格ライン |
| --- | --- | --- | --- |
| V1 ✅ | 1D 純伝導の共役解 | **静止流体層 (厚さ $d$, $k_f$) + 固体層 ($t$, $k_s$)** の 1 次元。$R_f=d/k_f$, $R_s=t/k_s$。<br>$T_w=\dfrac{T_\infty/R_f+T_b/R_s}{1/R_f+1/R_s}$ (**旧稿は重みが逆だった** — codex M7)。<br>検算: $R_f$=0.002, $R_s$=0.005, $T_\infty$=1000, $T_b$=300 K → **$T_w$=800 K, 両側 $q$=100 kW/m²**。<br>閉じた系なので **dual-time で定常化** ([[steady-localdt-acoustic-instability]]) → **依存診断の dual-time 解除が前提** (§4.3, codex 2 巡目 #1) | $T_w$ が解析解と **0.5 % 以内**、両側の $q$ 不一致が **0.1 % 以内**。**採取位相を明記** |
| V1b | 反復の性質 | (i) §4.2 の**限定モデル** (0 次元/1 次元) で $D_f$ と $\mathrm{Bi}$ を振る、(ii) **分布問題**: 面内伝導つき + **接線方向の交番温度モード** + **非一様 $D_f$** + **codex 3 巡目 #1 の反例** ($A_s=0.1I$, $H=[[2,-1],[-1,2]]$, $D_f=s\,\mathrm{diag}(1,10)$) | (i) 理論どおりの増幅率、(ii) **メリット関数 $\Phi$ + line search の受理**が交番モード・非一様 $D_f$・上記反例のすべてで**停止せず収束**する ($\Phi$ 単調減少、最大ノルムは途中で増えてよい)。**$D_f$ 単独の「上界だから安全」は主張しない**。(i) の結果を (ii) の保証に使わない |
| V2 | SU2 CHT | **まず V1 と同じ 1D スラブ**で SU2 multizone (固体ゾーン + `MARKER_CHT_INTERFACE`, `DIRECT_TEMPERATURE_*`) と一致させる。その後 case/48 へ拡張し、**固体格子収束とシェル近似誤差を別に測る** | 1D スラブ: $T_w$ **≤0.2 %**。case/48: $T_w$ ≤2 %・$\int q\,ds$ ≤3 %。**両ゾーンの残差と版・メッシュ・物性を記録** |
| V3 | 等温との整合 | case/48 で「共役解の $T_w(x)$ を `wallProfile` で与えた等温 run」と共役 run を比較 | $q_w,\delta^*,C_f$ が **≤0.5 %** |
| V4 | Phase 1 ↔ Phase 2 | 同一問題を外部ループとソルバ内で | $T_w$ **≤0.3 %**、$Q_w$ **≤1 %** |
| V5 | **公知データ (超音速)** | §4.9。**3 段に分けて誤差を分離する** (codex 2 巡目 #6): **(a) 実測壁温を与えた流体計算** (流体側の $h$ を当てる) → **(b) 公開条件による固体単独検証** (`fem2d` の検算) → **(c) 壁温を未知とした CHT**。`case/54` Mark II run 42 ($M_2$=1.04) を主、`case/53` C3X run 108 ($M_2$=0.90) を先行 | **各段の入力と合格を先に固定する** (codex 3 巡目 #5)。**(a)**: 実測 $T_w$ を課し、**原典と同じ温度基準で定義した $h$**・熱流束・壁圧を比較点ごとに比較。**(b)**: 実測 $T_w$ を**外周 Dirichlet**、公開冷却条件を孔 Robin として固体単独で解き、**外周の反力熱流束**を比較する (= **原典のデータ処理の再現検査**。入力した壁温に一致することを成果にしない)。**(c)**: 壁温を未知に戻して CHT。<br>(a)(b) が通ってから (c) を判定し、(c) は測定点 $T_w$ が**事前登録した帯**の中。帯は 冷却孔ごとの HTC 相関 / 計測断面の冷却剤温度推定 / 材料物性 / **表 VII の温度比 ±2 %** を項目別に立て、**合成規則 (単純和か二乗和か、相関の扱い) を計算前に決める**。**合うまで帯を広げる運用を禁止** |
| V6 | すきま適用と感度 | case/51 の薄肉ライナ。$K$ (50/200)、$D_f$ 方式、`local1d` vs `shell2d`。**壁関数併用は §4.3 の対象外宣言に合わせ、診断解除後の追加試験に分離** (codex 3 巡目 #2) | $K$・$D_f$ 依存が $T_w$ で **≤0.5 %**。面内項の寄与と放射の桁 (§4.10) を数値で併記 |
| G-cons | **熱収支** (C1) | 同一状態で 流体側 $\sum Q_{f,i}$ = 固体正味入熱 = 背面 + 端部流出 + 接合授受 | **熱量の 0.5 % 以内**。**これを満たさない run は合格させない** |
| G-if | **界面収束** (M10 + 2 巡目 #5) | 3 つを**別々に**判定する: ① 未緩和残差 $r=A_sT-b_s-Q_f$ の**局所最大ノルム**、② 固体内部残差、③ **絶対**温度更新量 $\max\|\Delta T_w\|$ [K]。既存 `check_convergence.py` は保存量列しか見ないので**拡張が要る** | **数式で書く** (codex 3 巡目 #3): ① は**局所面積で規格化した熱流束尺度** $\max_i |r_i|/A_i$ [W/m²] で見る (積分荷重 [W] の絶対閾値だけだと細分化で緩む)。合格は **(絶対条件 $\max_i|r_i|/A_i \le \epsilon_{\rm abs}$) かつ (相対条件 $\max_i|r_i|/\max_i|Q_{f,i}| \le \epsilon_{\rm rel}$) が $n_{\rm consec}$ 回連続**。欠損・非有限・**量子化による停滞は不合格**。$(\epsilon_{\rm abs},\epsilon_{\rm rel},n_{\rm consec})$ は **V1 の実測値を実装仕様の基準として固定**する。**実装は Phase 1 の合格判定より前** |
| G1–G4 | 収束 / 準定常 / 壁解像 / メッシュ | `check_convergence.py --segment` / `check_quasisteady.py` / `check_wall_resolution.py --target 1 --over-frac 0` / `check_mesh_quality.py` | VERDICT を README と応答に貼る。**旧 run の NOT CONVERGED を継承しない** |

**準定常許容は drift と振動幅の両方を締める** (codex 3 巡目 #3): `check_quasisteady.py` は
drift と振動幅を**別々に**判定する ([:293](../../solver_density_cuda/tools/check_quasisteady.py)、既定 5 % / 10 %)。
codex の実測: 末尾 `[99,101,101,99]` の系列は **drift を 0.04 % に締めても `STEADY` (fluct 2 %)**、
`osc` も 0.04 % にして初めて `OSCILLATING` になる。
→ **`--drift` と `--osc` の両方を比較許容の 1/5 以下**にし、$T_w$・実効 $Q_f$ だけでなく
**比較に使う局所温度・局所熱量**にも同じ判定を適用する。

**外したときの扱い**: 界面定義 (§4.3 の 4 形式)・更新式・固体抵抗・壁解像・固体格子のどれが効いたかを
切り分けてから結論する。**切り分け前に「一致」と書かない**。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result | `2026-09-20` | [`notes/reviews/2026-09-20-boundary-conjugate-heat-transfer-result.md`](../../notes/reviews/2026-09-20-boundary-conjugate-heat-transfer-result.md) | **NO-GO**, C1/M9/m1 | **全件採用** (反例つきで再現されており、うち 4 件は自分でも独立に気づいた): C1 保存的界面熱量 $Q_f=\sum F^E-C$ 未実装 → §5.1 #30 / M2 棄却時に別状態の $T$ と $Q_f$ を混ぜる + 初回の壁温不一致 (課したのは実測分布、ドライバは一様) → #31 / M3 `Df_safety=2` は上界でない (反例 固有値 −5.67) → #32 / M4 `fem2d` 外部連成が $k_s(T)$ を局所で解いていない (2.07 K 差で converged) → #33 / M5 収束ゲートの規格化が $\max(|Q_f|,|b|)$ で不釣合い 100 % でも合格 → #34 / **M6 V5 は同定データへの再適合で独立検証でない → §4.9 と case README を格下げ、「差の主因は遷移」を仮説へ撤回** → #35 / M7 積分量の `STEADY` が局所量を保証しない (C3X 局所 $q$ 402 点中 255 点 `DRIFTING`) → #36 / M8 共有角検査が連成中の競合を防げない → #37 / M9 同一メッシュ再開が 2D 最近傍 → #38 / M10 cell の `wallProfile` が壁面重心でなく内部セル重心 → #39 / m11 docs と残作業表の同期 → #40。**`accepted` への移動は取り消し、`active` のまま保存的連成の検証を先にやる** |
| 自由形式 (V5 実行中の診断) | `2026-09-20` | [`notes/reviews/2026-09-20-codex-c3x-cascade-diagnosis.md`](../../notes/reviews/2026-09-20-codex-c3x-cascade-diagnosis.md) + [`…-codex-node-isothermal-roe-frozen.md`](../../notes/reviews/2026-09-20-codex-node-isothermal-roe-frozen.md) | 判定なし (診断依頼)。実装欠陥 3 件 + 切り分け | **採用**: (1) node 等温壁の `roe` 凍結 → `main.cpp` の壁ピン順序修正 (EOS/物性/RANS 壁の前 + 陰解法最終更新後)。(2) `inlet_Pressure_dir` の根号・方向ゼロ割 → `boundaryCond_d.cu` ガード追加。(3) `convMethod: 1, limiter: 0` は 2 次 → 起動段を `convMethod: 0` に。(4) SU2 対照 → §5.1 #28 |
| plan (3 巡目) | `2026-09-19` | [`notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-4.md`](../../notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-4.md) | **GO-with-changes**, C0/M5/m1 (2 巡目 #3 一部・#4・#6・#7 は**解消**) | **全件採用**: #1 受理判定が収束反復を止める反例 → §4.2 をメリット関数 + line search へ / #2 提供側 plan への解除契約登録 + V6 壁関数の分離 → §4.3 + §5.1 #21 / #3 drift と osc の両方 + G-if の数式化 → §6 / #4 非連成等温壁との共有角競合 → §4.4b / #5 V5 (a)(b) の入力・比較・合格と不確かさ合成 → §6 V5 / m6 親 plan の旧判断 → §5.1 #27 |
| plan (2 巡目) | `2026-09-19` | [`notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-3.md`](../../notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-3.md) | **NO-GO**, C0/M6/m1 (1 巡目の C1/M4/M5/M9/m11/m12 は**解消**と再評価) | **全件採用**: #1 依存診断の対応範囲 (dual-time/周期の解除マイルストーン・過渡の $C$) → §4.3 + §5.1 #20 / #2 $D_f$ は上界でない (反例 5.675) → §4.2 + #21 / #3 `fem2d` の連成契約と単体試験 → §4.4d + #22 / #4 陰解法フック (`advanceImplicitSteady` 経路) → §4.6 + #23 / #5 ゲートの数値仕様と前倒し → §6 + #24 / #6 V5 の 3 段分解と不確かさ項目 → §4.9 + §6 V5 + #25 / #7 索引・親 plan の同期 → §5 S0 + #26 |
| plan (2 巡目 初回, **中断**) | `2026-09-19` | [`notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-2.md`](../../notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-2.md) (codex の利用上限で最終メッセージ無し。検算ログのみ) | 判定なし。検算で 3 件検出 | **全件採用**: 界面熱量の符号・組合せ ($\sum F-C$; 旧稿 $-60$ vs 正 $100$) → §4.3 / 成分セカント $D_f$ は過小評価でスペクトル半径 1.818 → §4.2 / 放射 1500 K は 287 kW/m² → §4.10。**上限解除後に 2 巡目をやり直す** |
| plan | `2026-09-19` | [`notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan.md`](../../notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan.md) | **NO-GO**, C1/M9/m2 | **全件採用**: C1 実効界面熱量 → §4.3 + G-cons / M2 固定点保存の更新式 → §4.2 / M3 $R_{\rm tot}$ と零空間 → §4.4a / M4 `wall_isothermal` + 属性・node 限定 → §4.6 + §2 / M5 ノード座標補間と verify → §4.5 / M6 facet 幾何・接合・軸対称 → §4.4b / M7 V1 の式と問題設定 → §6 V1 / M8 SU2 multizone・1D 先行 → §6 V2 / M9 旧 run は根拠にしない → §3 + §5.1 #4 / M10 界面ゲートと区間ハッシュ → §6 G-if / m11 不採用理由と性能主張の修正 → §4.1 + §4.5 / m12 docs 先行 → §5 S0 |

## 7. 影響範囲

- **forge 本体**: `boundaryCond.{hpp,cpp}` (`conjugate` 属性・プロファイル座標)、新規 `conjugateWall.{hpp,cpp}`、
  `input/solverConfig.{hpp,cpp}` (`conjugate:` 節)、`main.cpp` (フック)、`output/output.cpp` + `variables.hpp` (界面診断)、
  `cuda_forge/viscousFlux_d.cu` (診断書込みのみ)。
  **確認が要る既存経路** (codex M4): `mesh.cpp`:929 (`iso_wall_flag`)、`nodeWallDirichlet_d.cu`:130、
  `viscousFlux_d.cu`:954、`calcWallDistance_kdtree.cpp`:129、`timeIntegration_d.cu`:1025。
- **ツール**: `tools/solid_shell.py`, `tools/solid_fem2d.py`, `tools/cht_loop.py`,
  `tools/check_convergence.py` (界面列)、`tools/stage_manifest.py` (外部入力ハッシュ)。
- **ケース**: `case/52.conjugate_slab/`, `case/53.c3x_vane_cht/`, `case/54.markii_vane_cht/`, case/48 の派生。
- **docs**: `methods/boundary.md` (**実装前**) + `methods/index.md`、`procedures/su2-cross-check.md` (CHT multizone)、
  `procedures/recommended-settings.md`、`procedures/verification/`、`design/CAPABILITIES.md`。

## 8. 未確定事項 (ユーザ確認)

1. ~~一次の適用先をどれにするか~~ **決着 (2026-09-19, §4.9)**: ユーザ指示「公知文献の試験結果があるもの」
   「超音速系だとなお良い」により、**NASA CR-168015 の Mark II 翼 (超音速出口) を主・C3X を従**とする。
2. **固体物性の与え方**: 材料 DB を作るか、ケースごとに `solid.json` 直書きか (初版は直書きを想定)。
3. **放射**: §4.10 のとおり解かない。深いキャビティで桁が無視できないと分かったら別 plan にするか。

## 9. 変更ログ

- `2026-09-19` — 初稿。作戦 A–D の比較、界面契約、連成の安定性、薄肉シェル、検証 V1–V5 を定義。
- `2026-09-19` — codex plan レビュー 1 巡目 **NO-GO** (C1/M9/m2) を**全件採用して全面改訂**:
  界面熱量を「物理境界流束 + 拘束反力」に変更 (§4.3)、固定点を保存する更新式へ (§4.2)、
  $R_{\rm tot}=t/k_s+R_{\rm back}$ と零空間処理 (§4.4a)、facet 単位の幾何と接合 (§4.4b)、
  `wall_isothermal` + `conjugate` 属性・node 限定 (§4.6)、プロファイルのノード座標補間 (§4.5)、
  V1 の解析解の重み修正と問題設定の固定、V2 の 1D 先行と SU2 multizone、界面ゲート新設 (§6)、
  docs 先行 (§5)、性能・不採用理由の主張を緩和 (§4.1)。
- `2026-09-19` — ユーザ指示により**一次適用先を公知試験データ (超音速出口のタービン翼) に確定** (§4.9、§8-1 決着)。
  固体バックエンドに `fem2d` (一般 2D 領域) を追加 (§4.4c)。
- `2026-09-19` — codex 2 巡目 (利用上限で中断) の検算 3 件を採用: 界面熱量の符号・組合せ (§4.3)、
  成分セカント $D_f$ の過小評価と発散例 (§4.2)、放射の桁の誤記 (§4.10)。**2 巡目は上限解除後にやり直す**。
- `2026-09-19` — codex plan 2 巡目 **NO-GO** (C0/M6/m1) を**全件採用**。1 巡目の Critical は解消と再評価された。
  主な変更: $D_f=k_{\rm eff}A/d_1$ を「上界」から**初期推定値**へ格下げし受理判定 + 退避を既定に (§4.2)、
  依存診断の **dual-time / 周期の解除マイルストーン**と過渡の拘束反力式 (§4.3)、
  `fem2d` の連成契約 $(K_s+E^{\mathsf T}D_fE)u=\dots$ と単体試験 (§4.4d)、
  **陰解法フックの位置** (`advanceImplicitSteady` 経路。旧稿は explicit RK を指していた, §4.6)、
  G-cons/G-if の**数値仕様と Phase 1 前倒し** (§6)、V5 の **3 段分解**と run 42/108・$M_2$ 1.04/0.90 (§4.9)。
- `2026-09-19` — codex plan 3 巡目 **GO-with-changes** (C0/M5/m1) を**全件採用**。
  受理判定を最大ノルム単調減少から**メリット関数 $\Phi=r^{\mathsf T}(A_s+D_f)^{-1}r$ + line search** へ
  (旧規則は収束する反復を棄却する反例あり, §4.2)、提供側 plan への**解除契約の登録** (§4.3)、
  **非連成等温壁との共有角競合の拒否**と global CV ID による所有 (§4.4b)、
  準定常は **drift と osc の両方**・G-if の数式化 (§6)、**V5 (a)(b) の入力と合格条件・不確かさ合成規則** (§6 V5)。
  → **実装着手は §5.1 #20–#27 を確定させてから** (Phase 2 は V1–V3 合格後)。
- `2026-09-19` — **S0 完了**: [`methods/boundary.md`](../../methods/boundary.md) に「共役熱伝達 (CHT)」節を追加し、
  界面契約・固体モデル・受理判定・`wallProfile`・**起動時拒否条件**・ゲートを実装前の契約として固定。
  `methods/index.md` の状態欄と、提供側 plan への解除マイルストーン登録も完了。
- `2026-09-19` — **S1/S2 完了・S3 部分完了**。`conjugateWall.{hpp,cpp}` (符号規約の正本 + 第一/第二内部点 +
  界面診断 + 共有 CV の壁温競合検査)、`output.interfaceDiag`、`applyWallProfiles` を実装し、
  case/48 の 4 run (`run_0013`–`run_0016`) で検証。**主な実測**: 診断は解をノイズ床内でしか動かさない /
  $d_1$ が第一層厚と一致 / **`q_recon` ≡ `q_compact` (1.6e-7)** で、~2.7 % の食い違いは
  **2 次片側差分との差 = 後処理の選択**だった / `wallProfile` が場に入る / 共有角の拒否が効く。
- `2026-09-19` — **S4 (固体シェルソルバ) 実装・検証**。`tools/solid_shell.py` (`local1d`/`shell2d`、
  $R_{tot}=t/k_s+R_{back}$、零空間検出、$k_s(T)$、固定点保存の連成反復 + メリット関数受理 + Anderson 加速) と
  `tools/test_solid_shell.py` (解析解 5 種)。**設計に反映した実測 2 件**:
  (a) **四角形の片側対角線分割は境界で O(1) の荷重不均衡を作り収束次数を 1 次に落とす** → 両対角線の平均で 2.00 に回復。
  (b) **スカラー $D_f$ では非対角な流体応答を捌けず素の固定点反復は収束しない** → Anderson 加速を既定に、
  かつ **$\Phi$ の比較は $D_f$ 固定区間で**・収束判定は物理量で (重み変更による見かけの降下で $D_f$ が発散する)。
- `2026-09-19` — **S5 (外部弱連成ループ) 実装・V1 合格**。`tools/cht_loop.py` と `case/52.conjugate_slab`
  (`gen_mesh.py` / `verify_v1.py` / `solid.json`)。`run_0001_v1_slab`: $T_w$ 誤差 **0.025 %**・
  **両側 $q$ 不一致 0.0053 %** で **VERDICT PASS** (11 反復)。
  検証ケースの設計要件 (一様メッシュ・低圧化・`cfl_pseudo` 5) を README に記録。
- `2026-09-19` — **Phase 2a (ソルバ内連成, `local1d`) 実装・検証**。`conjugate:` ブロックと
  `ints: {conjugate: 1}`、抵抗加重平均の更新、起動時の拒否条件、`conjugate_Tw_*.csv` による再開。
  case/52 で **$T_w$ 誤差 −0.152 % of rise・両側 $q$ 不一致 0.0002 %**、外部ループとの差 0.18 % of rise (V4 許容内)。
  **外部ループは乱流平板 (case/48 run_0017) では 12 反復で未収束**で、生産はソルバ内という位置づけを実測で確認。
- `2026-09-20` — 平板のソルバ内連成を**静定まで回した** (`run_0019_cht_insolver_cont`, +60000 step)。
  `check_quasisteady` **ALL STEADY** ($T_w$ x=0.5 m 569.87 K、$Q_w$ 60.39 kW/m、両側不一致 31.6 W/m² = 2.1e-4)。
  **`check_convergence` は `NOT CONVERGED (stalled/plateau)`、`--from-floor` は参照未収束のため REFUSED**
  → case/48 系列は残差では合格にできないと確定。界面時系列ツール `tools/cht_wall_series.py` を追加。
- `2026-09-20` — **`fem2d` 実装・検証** (`tools/solid_fem2d.py` + `test_solid_fem2d.py`)。
  界面作用素は**内部節点を消去した Schur 補元**にしてシェルと同じ契約 ([W/K], [W]) に揃えた。
  円環 (Robin 孔) の解析解に対し **rate 1.97**、Schur = 全系解 5e-12、
  **物理的な (対角優位な) 流体応答なら連成は 7 反復**で厳密固定点 → §4.2 の反例は最悪ケースであることを確認。
  これで **V5 (Mark II / C3X) の固体側が揃った**。
- `2026-09-20` — **V5 の一次資料を抽出** (`case/53.c3x_vane_cht/`): NASA CR-168015 から C3X 78 点 / Mark II 60 点の
  翼型座標と試験条件 (run 108 = $M_2$ 0.90 / run 42 = $M_2$ 1.04) を**検証つきで**取り出した。
  スキャン OCR なので `cm == 2.54 in` を全点に課し、通らない 11 点は**ページ画像を読んで補修**、
  欠測があればスクリプトが落ちるようにした。**報告の不整合 2 件** (表 VIII/IX の SI 圧力列、C3X 点 29 の inch 値) を
  case README に記録。残りは $k_s(T)$・冷却孔条件・測定壁温・カスケード幾何。
- `2026-09-20` — V5 のデータ整備を継続: **カスケード幾何 (表 IV)・冷却孔 10 本 (図 6/7)・
  材料 $k_s(T)$ (報告に無いため外部出典)・run 108 の実測壁温と $h$ (75 点)** を抽出・検査。
  正規化は **$T_w/811\,\mathrm{K}$ と $h_0$=1135 W/m²K**、$h$ は**ガス全温基準**と確認。
  **設計変更**: **冷却孔の冷却剤条件は公開されていない**ので、V5 段 (b) は「公開条件で固体を解く」から
  **「公開量 (実測 $T_w$ と $h$) から内部条件を逆算し、その不確かさを帯に積む」**に変更した (§5.1 #25)。
- `2026-09-20` — **V5 の固体側が立った** (CFD 不要の段): `case/53` に固体メッシュ生成
  (`gen_solid_mesh.py`: 外周 240 節点を流体と共有できる形で配置、孔は物理曲線) と
  **内部条件の逆算** (`infer_internal_bc.py`) を実装。run 108 で **$h_c$ = 313–4207 W/m²K**
  ($T_c$=300 K)、**外表面温度の残差 RMS 10.7 K**、熱収支は機械精度で一致。
  **V5 の帯に入る新事実**: (a) $h_c$–$T_c$ の縮退 ($T_c$ 280–350 K で RMS ほぼ不変、400 K 以上は非物理)、
  (b) **残差床 10.7 K はマッピングに鈍感**でモデル (10 孔 Robin + 2D 伝導) の限界、
  (c) 後縁 35/240 節点はデータ範囲外で検証対象外。
  → **段 (c) で「実測と一致」を主張できる下限がこの 10 K 級**であることを、計算前に登録する。
- `2026-09-20` — **V5 段 (a) の翼列 run を投入するなかで forge 本体の欠陥 2 件を直した** (どちらも codex の
  自由形式診断 [`notes/reviews/2026-09-20-codex-c3x-cascade-diagnosis.md`](../../notes/reviews/2026-09-20-codex-c3x-cascade-diagnosis.md) と独立再現):
  (1) **node 等温壁の `roe` が凍結**し、EOS が密度床に落ちた温度 (9.19e6 K) で $\mu$ が 151 倍・壁 $\omega$ が
  8.1e15 になっていた → `main.cpp` で壁ピンを `dependentVariables`/`gasProperties`/RANS 壁の**前**へ移し、
  陰解法の最終更新後にも適用。(2) **`inlet_Pressure_dir` が $P_c>P_t$ で根号負・方向ベクトル長 0 で 0/0** →
  `boundaryCond_d.cu` にガード ([`methods/boundary.md`](../../methods/boundary.md) 「例: 全圧固定流入」)。
  あわせて **`inlet_uniformVelocity` は全温を固定しない**ことを実測で確認 (case/53 で入口 $T_t$ が
  786→738 K)。$T_g$ 基準の $h$ を比較する V5 では**全圧入口を使う**と決めた (§4.9)。
- `2026-09-20` — **V5 段 (a)/(c) を実施** (case/53 C3X, run 108)。
  - **段 (a) 実測 $T_w$ → $h$** (`run_0003_smoothgeom`): 負圧面**遷移後 +4.0 % / rms 11.0 %**、
    正圧面 +10.6 %、**負圧面層流域 +40.5 %** (低 Re SST に遷移モデルが無いため)、全体 +16.6 %。
    入口 $M_1$=0.167 (報告 0.17)、出口 $M_2$=0.901 (報告 0.90)。
  - **段 (c) 連成 CHT** (`run_0004_cht`, `cht_loop --solid-mode fem2d`): $T_w$ 523.7–617.2 K
    (実測 512.6–611.8 K)、**バイアス +11.9 K / RMS 18.8 K** ($s/S\le0.87$)。遷移後の負圧面は
    **+6.8 K / rms 12.1 K** で、段 (b) 同定の残差床 **9.9 K** とほぼ同桁。
  - **帯の照合**: 内部条件の推定由来は小さい ($T_c$ 280→350 K で 2.2 K、$k_s$±3 % で 2.0 K;
    `run_0005_band_*`)。**差の主因は CFD 側 (遷移モデル無し)** と結論。→ §5.1 #28 (SU2 対照) と
    #29 (遷移モデル) を残作業に置く。
  - 途中で直したもの: 翼型の折れ角 60° (`smooth_profile.py`、円弧中心を $R$ 固定 LS に)、
    `cht_loop` の $D_f$ 初期値 (第一セル伝導 → **熱伝達係数 × 面積**。31 倍過大で反復が止まっていた)、
    `cht_loop` の `fem2d` バックエンド接続 (§5.1 #24 の残り)。
- `2026-09-20` — **Mark II (超音速出口, run 42) も V5 段 (a)/(c) を通した** (case/54)。
  - 実測の転記 (68 点、PDF p.150)。**300 dpi 2 値スキャンで判読不能な 9 セルは空欄**にして
    推測で埋めない (検査は欠測率 15 % 超で不合格)。位置列は同じ翼の他 run (43) と一致で相互検査。
  - 段 (a): **正圧面 +0.8 %**、負圧面 層流域 **+74.8 %**、遷移後 +13.7 %。
  - 段 (c): **正圧面 −0.2 K / rms 12.0 K** (同定の残差床 7.78 K)、負圧面 層流域 +42.8 K、全体 +16.6 K。
  - **結論 (両翼共通)**: 連成そのものは正圧面と C3X 遷移後で確認できた。差の主因は
    **遷移モデルの欠如**で、超音速負圧面 (Mark II) は層流区間が長いぶん C3X の倍効く。
    → §5.1 #29 を「V6 前に遷移位置を与えた感度計算」から**必須**へ格上げ。
  - 途中で直したもの: 周期境界の対応付けに**距離と 1 対 1 性の検査**を追加 (`mesh.cpp`。
    ピッチの単位間違いが無検査で通り「完走したが M=1264」になった)、閉塞翼列の**背圧段階起動**、
    超音速出口の**下流ブロックを短くする** (長いと出口面がまだらになり閉包が毎 step 入れ替わる)。
- `2026-09-20` — **#28 SU2 対照を実施し、遷移仮説が確定した**。同一 `.geo` からの同一メッシュ・
  同一 BC・一様壁温 566 K で SU2 8.5 (RANS/SST 低 Re, $y^+\le1.85$) と比較。**負圧面層流域は
  SU2 も +45.8 % 過大** (forge +42.7 %)、**forge vs SU2 は −3.0 %/6.0 %**、遷移後は **−0.6 %/2.1 %**。
  → 2026-09-20 に「仮説」へ格下げした「差の主因は遷移モデルの欠如」は、**独立ソルバで裏付けられた**
  (case/53 README「SU2 対照」節)。壁総熱量も forge 48.05 kW/m vs SU2 49.3 kW/m (+2.6 %)。
  **残る solver 間差** (正圧面 −6.0 %/7.1 %、層流域 −3.0 %/6.0 %) を §5.1 #41 に起票。
- `2026-09-20` — **#41 (forge–SU2 の残差) を切り分け、2 つの仮説を棄却した**。熱流束の定義・幾何・
  層流物性・第一内点の $\mu_t$ を交絡から外したうえで、`sstEnergyIncludesK: 1` (SU2 形の
  エネルギー式) を試したが **−4.31 % → −4.50 %** で改善せず棄却 (`run_0010_energyK`)。
  層流単独対照は **forge が定常解に落ちない** (負圧面剥離で非定常、残差 0.0 dec) ため不成立
  (`run_0011_laminar`、破棄予定)。**残る事実は「速度は ±1 % 一致なのに $T(y)$ が壁から 20 µm で
  系統的にずれる」**こと。次の候補は対流スキーム・node 壁半 CV のエネルギー流束・BL 発達履歴。
- `2026-09-20` — **「forge が層流で収束しない」は欠陥ではなく物理だった** (ユーザ指摘
  「層流で収束しない? まずくない?」→ codex 診断
  [`notes/reviews/2026-09-20-codex-laminar-nonconvergence.md`](../../notes/reviews/2026-09-20-codex-laminar-nonconvergence.md) →
  dual-time で確認)。後縁 $Re_D$=4.2e4、dual-time ($\Delta t$=50 ns, 200 µs) で**後縁 2 mm 下流の
  2 点が逆位相・PSD ピーク 29.3 kHz ($St$=0.29)** のカルマン渦列を確認 (`run_0018_lam_dualtime`)。
  **$\Delta t$ を 25 ns に半減しても ピークは 29.9 kHz** (差 2 %、分解能 3.3 kHz 以内、振幅も 5 % 以内)
  で**時間解像できている** (`run_0019_lam_dt_half`、`maxCFL` も 14.5–14.9 で推奨域)。$St$=0.295。
  `cfl_pseudo` を 0.05 まで下げても `nStepInner` を 30 にしても残差が動かないことと整合する。
  **あわせて 2 つの記述を訂正**: (a) 変動は負圧面剥離でなく**後縁後流**、(b) **SU2 も残差基準では
  `NOT CONVERGED`** (`rms[RhoE]` +0.465)、定常なのは積分量。ただし場の静定度は SU2 44 Pa に対し
  forge 2934 Pa で **67 倍**違う。**SST の非収束は別物** (スナップショット間 $P$ rms 差が
  15.9 / 1.07 Pa と 200–2700 倍静か)。
  副産物: `cfl_pseudo` 4.0 は**出口 (`outlet_statPress`) で圧力床に落ちて発散**する (層流・SST とも)。
  → §5.1 #42。
- `2026-09-20` — **#42 を撤回した**。「出口 BC が擬似 CFL を律速し生産 run の時間を食っている」は誤り。
  発散の**場所**は出口だが**原因は 2 次再構成** (1 次なら `cfl_pseudo` 4.0 で安定、リミッタ 4 種すべて
  発散、`outlet_statPress` のゴースト緩和も無効)。`implicitRelax` 0.7 / `nStepInner` 30 で上限は
  外せるものの、**SST 本番 20000 step の最終残差は現行 `cfl_pseudo` 0.5 が 4.6e−06 に対し
  cfl 2.0 で 4.5e−03・cfl 8+relax 0.7 で 1.6e−03 と 300–1500 倍悪く、所要時間はほぼ同じ**。
  現行設定が最良。速度を詰めるなら擬似時間の CFL 以外を見る。
- `2026-09-20` — **C1 (保存的界面熱量) を定常 node 等温壁について実装した** (§5.1 #30)。
  `iface_q_eff` $=(R^{raw}-F_w)/A$ を壁ダンプに追加し、`cht_loop --flux q_eff` を通した。
  採取点は 2 つ: 壁半割面が `res_roe` に入れた寄与 (`viscousFlux_d.cu`) と、
  **壁残差射影の直前**の `res_roe` (`nodeWallDirichlet_d.cu`)。後者は 0 化されるので
  この位置で採るしかない。**V1 で符号・絶対値を再現** (解析解 81.3090 に対し 81.1837 W/m²、−0.154 %)。
  **コンパクト差分との差は無視できない**: case/53 で bias +2.26 % / rms 3.28 % / 局所 19.1 %、
  連成壁温が +11.5 → +13.5 K 動く。→ **V1 と V5 は「保存的 CHT」の枠で語れるようになった**
  (定常・node・等温壁・非周期壁ノードの範囲で)。dual-time / 周期壁ノード / 軸対称は `NaN` を出す。
  副産物: 新しい bvar は `mesh.hpp` の `bplaneValNames` にも足さないと確保されず、
  D2H コピーで segfault する (`valueTypes` 側の宣言だけでは足りない)。
- `2026-09-20` — **ゲートを 2 つ実装した** (§5.1 #34 / #30b)。**G-if**: 規格化を $\max|Q_f|$ のみに直し
  (codex M5 の反例を `test_solid_shell.py` T7 に回帰試験化)、絶対残差・相対残差・**固体内部残差**・
  温度更新・退避していないことを独立に満たすときだけ収束とする。**G-cons**: `check_cht_balance.py` が
  $\varepsilon/\max(\sum|Q_f|,Q_{\rm floor})$ で判定し、`iface_q_eff` に `NaN` があれば不合格。
  **効いている**: `run_0022_cht_qeff` (1 反復 4000 step) は G-cons **FAIL 2.77 %**・G-if `res_rel` 3.5 % と
  一致し、1 反復 16000 step に上げると `res_abs` 487 → 2.4 W・`res_solid` 1.02 → 1.2e-4 W まで落ちた
  (`run_0023_cht_tight`)。**律速は CFD 側のノイズ床**で、固体でも反復でもない。
- `2026-09-20` — **G-cons を V1 とソルバ内連成にも通し、C1 を実証した**。
  V1 (`case/52`, `local1d`) は不釣合い **0.000125 %** (2e-6 / 1.62 W)。
  case/53 の**同一状態で界面熱量の定義だけ変える**と、`q_eff` 0.027 % / `q_2nd` 0.29 % が PASS、
  **`q_compact` 2.66 % と `q_recon` 2.66 % が FAIL**。**コンパクト差分で連成すると、固体が持ち去る熱を
  流体が供給していない解を通してしまう**ことを数値で示した (codex C1 の主張どおり)。
