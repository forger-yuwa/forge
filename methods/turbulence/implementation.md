# SST 乱流モデルの実装方針

## 1. 目的

本章では、forge の `solver_density_cuda` に低 Re SST を組み込む際の
実装責務と変更箇所を整理する。理論は
[theory.md](theory.md) を参照し、本章ではコード上の配置と責務分担に絞る。

## 2. 実装の基本方針

- 初期導入は **explicit advection-only** に限定する
- 軸対称は SST 固有の diffusion / source を子 plan に分離する
- 既存 LES (WALE) と laminar 経路を壊さない形で dispatcher を拡張する
- 既存の `wall_dist`, `vis_turb`, 勾配再構成、出力系を最大限流用する
- 5 方程式 NS 本体と scalar transport を分離し、`k`, `omega` は後者の最初の 2 本として扱う
- scalar diffusion は face ベースで実装し、`vis_lam + sigma * vis_turb` を使う
- `scalarTransport_d.*` は共通輸送コアに限定し、SST 固有の boundary / source / closure は model-specific layer に分ける

## 3. 主要変更点

### 3.1 状態変数

`variables.hpp` のセル変数に次を追加する。

- 保存変数: `roK`, `roOmega`
- stage / backup 変数: `roKN`, `roOmegaN`, `roKM`, `roOmegaM` など
- 勾配: `dKdx`, `dKdy`, `dKdz`, `dOmegadx`, `dOmegady`, `dOmegadz`
- 残差: `res_roK`, `res_roOmega`

初期フェーズでは implicit 用の `dq_*`, `rhs_block_*`, `diag_block_*` は
拡張しない。

### 3.1.1 初期乱流場 (`kInit` / `omegaInit`)

IC (`valueFileName`) に `roK`/`roOmega` が無い場合 (非 SST 場・メッシュ変換直後など)、
`variables::readValueHDF5` は `roK = ρ·kInit`, `roOmega = ρ·omegaInit` で初期化する
(`solverConfig` の `turbulence.kInit` / `turbulence.omegaInit`、既定 0)。

**`ω=0` 初期化は `μ_t = ρ a₁ k / max(a₁ ω, S F₂)` が ill-posed で SST が cold start から
発散しやすい** (特に node `wallTreatmentSST=1`: ω₀=0 だと inlet コーナーで ω が暴走)。
非 SST IC から SST を始めるときは inlet と同じ freestream 値 (例 `kInit=1.0`, `omegaInit=1000.0`)
を設定して 0 初期化を避ける。既定 0 では `ρ·0 = 0` で従来動作 (ビット不変)。
k/ω を持つ収束場から restart する場合は不要 (IC の roK/roOmega をそのまま使う)。

### 3.2 設定

`solverConfig` の turbulence 設定は、現在の `LESorRANS` と `LESmodel` に加えて、
RANS 用モデル選択と自由流乱流量入力を持てるよう拡張する。

初期実装で想定する入力項目:

- `LESorRANS`: `0=no model`, `1=LES`, `2=RANS`
- `LESmodel`: 既存 WALE 用
- `RANSmodel`: `1=SST`
- `turbIntensity` または直接 `kInf`, `omegaInf`
- `turbLengthScale` または等価な自由流尺度

### 3.3 メインループ

メインループの基本順序は
[methods/architecture/overview.md](../architecture/overview.md) を維持する。

初回マイルストーンの概念的な順序は次である。

1. 保存変数から primitive 量を更新
2. 境界条件適用
3. NS 用の勾配計算
4. SST から `vis_turb` を更新
5. NS の対流流束
6. generic scalar transport による $k$, $\omega$ の対流流束
7. model-specific scalar boundary / source layer
8. NS の粘性流束
9. explicit time integration

既存コードでは `convectiveFlux_d_wrapper` と `viscousFlux_d_wrapper` が
NS のみを想定しているため、初期実装では `scalarTransport_d.*` を追加し,
まず 1 次 upwind の advection-only を別 wrapper として差し込んだ。
その後、同じ wrapper に face ベースの diffusion を追加した。
ただし source と壁面生成ロジックは `scalarTransport_d.*` に埋め込まず、
`ransBoundary_d.*` と `ransSource_d.*` のような model-specific layer に分離する。
さらに k/ω への適用 (descriptor 構築・有効化判定・勾配) も共通コアから切り出し、
`ransTransport_d.*` にまとめた。

### 3.4 渦粘性

`turbulent_viscosity_d_wrapper` は現在 WALE 専用である。
ここを次の dispatcher に再編する。

- `LESorRANS == 0`: `vis_turb = 0`
- `LESorRANS == 1`: 既存 WALE
- `LESorRANS == 2 && RANSmodel == 1`: SST

**node-centered での WALE (`LESmodel==1`)**: `WALE_d` は DOF (cell/node) ごとに速度勾配・体積・
`wall_dist` を読み `vis_turb[ic]` を書くだけで `cfg.discretization` に依存しない。したがって SST と
同じ `dimGrid_normalcell` グリッドで起動すれば node の双対 CV (median-dual) も自動で網羅する
(以前は `dimGrid_cell` 固定で node では起動されず `vis_turb=0` = LES OFF 相当だった)。LES 時の
`vis_turb` は粘性流束で $\mu_{\rm eff}=\mu+\mu_{\rm sgs}$ として効き、RANS (SST) と同じ `vis_turb`
配列を共有する。node + KEEP + WALE の LES 構成は
[`plans/active/convection-keep-revive-node.md`](../../plans/active/convection-keep-revive-node.md)。

**2026-07-19 の 2 バグ修正** ([turbulence-wale-fix](../../plans/accepted/turbulence-wale-fix.md)):
① **壁なしメッシュで不活性**: converter が `wall_dist≡0` のまま → $L_s=\min(\kappa d, C_w\Delta)=0$
→ `vis_turb≡0`。solver 読込後に「全 CV で ≤0」なら 1e30 充填で修正 (壁ありメッシュ不変)。
② **Sd テンソル誤式**: $S^d_{ij}$ は速度勾配の**行列 2 乗** $\bar g^2_{ij}=g_{ik}g_{kj}$ から作る
(Nicoud-Ducros 1999) が、成分ごとの 2 乗になっていた (純せん断で $S^d=0$ の設計特性を破る)。
**適用指針**: 64³ TGV Re=1600 (h≈8η) では本物の WALE は遷移を先食いし ILES より悪化 —
**解像/遷移流は `LESorRANS: 0` + KEEP matrix ES 散逸 (σ=0.02, keepDissJump=2) を推奨**。
WALE は未解像の高 Re 乱流用オプション (適用検証は未実施)。

**σ-model (`LESmodel: 2`, Nicoud+2011)**: $\nu_t=(C_\sigma\Delta)^2\,\sigma_3(\sigma_1-\sigma_2)(\sigma_2-\sigma_3)/\sigma_1^2$
($\sigma_i$=速度勾配の特異値, $C_\sigma=1.35$, $\Delta=V^{1/3}$, wall_dist 不要)。2成分流・純せん断・
剛体回転・等方/軸対称伸長で厳密に $\nu_t=0$ (検証: `tools/verify_sigma_model.py`)。ただし 64³ TGV
では発達後の散逸が WALE より強く (K/K0(10) −5.5%)、**解像/遷移流の推奨は依然 ILES+ES**。
本領は壁乱流・回転流・未解像高 Re (適用検証は今後; [turbulence-sigma-model](../../plans/accepted/turbulence-sigma-model.md))。

### 3.5 境界条件

境界条件は既存の `wall`, `wall_isothermal`, `inlet_*`, `outlet_*`, `slip` などに
対して turbulence scalar を追加定義する。

ただし責務は 2 層に分ける。

- 既存 NS BC 層: 速度・圧力・温度など主流れ変数を処理する
- turbulence scalar BC 層: `kb`, `omegab` の生成と `kg`, `omegag` の ghost 反射を処理する

初期方針:

- 壁: `wall` / `wall_isothermal` では turbulence scalar BC 層がまず `kb`, `omegab` を生成し,
	ghost セルへ反映する。等距離 ghost を前提にすると、Dirichlet wall 値 $\phi_w$ に対して
	$\phi_g = 2\phi_w - \phi_c$ を使うため、$k_w = 0$ なら $k_g = -k_c$ になる。
	`\omega_w` は壁距離と粘性から与え、ghost 側は $\omega_g = 2\omega_w - \omega_c$ で閉じる。
- 流入 / farfield: 既知の `kb`, `omegab` を与え、ghost 側は同じ反射式で閉じる
- 流出: 1 次外挿ベース
- `slip` / `axis` / `periodic`: no-gradient / copy を使う

現行コードの対応は次のとおり。

- NS BC 層が既存の `wall_d` / `wall_isothermal_d` / `inlet_*` / `outlet_*` / `slip_d` / `periodic_d` を持つ
- turbulence scalar BC はこの後段に差し込む `ransBoundary_d.*` へ移していく
- `bcondConfFormat` の `wall` / `wall_isothermal` に `kb` と `omegab`、`inlet_*` に `k` と `omega` を持てる

この説明は、壁面がセル中心と ghost セルの中点にある等距離配置を前提にしたもの。
現在の実装では、`wall_d` / `wall_isothermal_d` が `kb` / `omegab` を壁値として受け取り、
ghost セル側へ反射して入れる。
将来、非対称な距離配置にするなら、理論側の距離比の式へ切り替える。

**node 入口 (Dirichlet スカラー境界) の保存量整合と残差除外** (`rans_dirichlet_scalar_boundary_d`):
node-centered では境界ノードが実 DOF なので、`omega[ic]=omegab`/`k[ic]=kb` の **primitive ピンだけでは
保存量 `roOmega`/`roK` が解かれ続けて不整合**になる (実測 `roOmega/ro` vs `omega` 888727、`roK/ro` vs `k` 7729)。
さらにその残差 `res_roOmega`/`res_roK` が `rms` を汚染し、入口×壁コーナーで最大化して node の `rms_roOmega`
プラトーを作る。よって壁ピン (§3.7) と同形に、node 入口では **保存量 `roOmega[ic]=ρ·omegab`/`roK[ic]=ρ·kb` を
整合ピン + `res_roOmega`/`res_roK`/`src_jac_*` を 0 化** する (フラグ `scalarDirichletPin` で識別、ransSource で除外)。
cell モードは ghost 経由で正しく課されるため不変。設計詳細は
[`plans/active/turbulence-node-inlet-dirichlet-conserved.md`](../../plans/active/turbulence-node-inlet-dirichlet-conserved.md)。

### 3.6 出力と診断

`output_cellValNames` に次を追加する。

- `roK`, `roOmega`
- 必要なら primitive として `k`, `omega`
- 既存 `vis_turb`

`gatherResidualSnapshot` は現状 5 方程式固定なので、初回段階では既存の
主流れ residual を壊さないことを優先し、必要なら scalar residual は別ログで扱う。

### 3.7 automatic (enhanced / y⁺ 非依存) wall treatment

理論は [`theory.md`](theory.md) §6.5。`wallTreatmentSST`
(`solverConfig` の `turbulence.wallTreatmentSST`, 0:wall-resolved / **1:automatic [既定]**)。
**既定は 1** (2026-06-28 に 0→1 へ変更, user 指示)。`0` で §6.1 の wall-resolved 型 (`60ν/β₁y²`) に戻せる。
あわせて node Dirichlet 既定も ON: `mesh.nodeWallDirichlet=1`・`turbulence.nodeKwfDirichlet=1` (§3.7 適用先・theory §6.5(e))。
**注**: 既定変更で y⁺~1 wall-resolved 前提の cell 検証 (case/26 $C_f$, flat-plate 回帰 等) は automatic に切替わるため要再検証。

責務分担とソース対応:

- **`ransWallFunction_d.*` (新設)**: 壁面ごとに Reichardt 則の逆解き (Newton 数回) で
	摩擦速度 $u_\tau$ を解き、`bc.bvar_d` の `utau` / `ypls` / `twall_*` を埋める。
	$u_\tau$ は速度・`wall_dist`・$\nu$ のみで決まるため勾配計算前に評価でき、後段の
	$\omega$ 壁 BC と運動量壁せん断の両方で共有する。`applyRansScalarBoundaries`
	(=`ransBoundary_d_wrapper`) の壁ループ直前で `wallTreatmentSST==1` 時のみ呼ぶ。
	`ransWallFunction_d.*` はあわせて新セル変数 `wf_pk`
	(wall-function 生産 $P_k = \rho u_\tau^4/\nu\cdot g(1-g)$, $g=du^+/dy^+(y^+_1)$) を
	wall-adjacent セルに書く (非 wall セルは `-1`)。
- **`ransBoundary_d.*` の wall scalar カーネル**: `wallTreatmentSST`, `utau`, `roOmega` を受け、
	mode 1 で $\omega_w = \sqrt{\omega_{vis}^2 + \omega_{log}^2}$
	($\omega_{vis} = 6\nu/\beta_1 y^2$, $\omega_{log} = u_\tau/(\sqrt{\beta^*}\kappa y)$) を
	**wall-adjacent セル値 `omega[ic]`/`roOmega[ic]` にピン留め** (§6.5(b)) し、$k$ は
	zero-gradient ($k_g=k_c$, §6.5(b'))。mode 0 は現行 $60\nu/\beta_1 y^2$ と $k_w=0$ Dirichlet。
- **`ransSource_d.*`**: `wallTreatmentSST` と `wf_pk` を受け、mode 1 かつ `wf_pk>=0` のセルで
	標準 P_k を `wf_pk` に置換 (§6.5(d))。さらに ω ピン留めセルの `res_roOmega`/`src_jac_omega` を
	0 化 (Dirichlet セルの残差は 0; `rms_roOmega` 汚染回避)。
- **`viscousFlux_d.*` の `viscousFlux_wall_d`**: `wallTreatmentSST` と `utau` を受け、
	mode 1 で接線せん断を modeled $\boldsymbol{\tau}_w = \rho u_\tau^2 \hat{\mathbf e}_t$ に
	置換 (法線粘性項・熱流束は不変、no-slip なので壁せん断仕事 0)。`twall_*_b` / `ypls_b` は
	この modeled 値で上書き出力。mode 0 は現行の分子勾配式。

`utau` は `wall` / `wall_isothermal` の `bvar` 初期化リスト (`boundaryCond.hpp`) と
`mesh.hpp` の `bplaneValNames` マスターリストの**両方**に追加する (片方だと device 未確保で
illegal memory access)。`wf_pk` は `variables.hpp` の `cellValNames` に登録し、
`boundaryCond.cpp` の `applyRansScalarBoundaries` が bc ループ前に
`initWallFunctionPk_d_wrapper` で全セル `-1` に初期化する。`ypls`, `twall_x/y/z` は既存。

`wallTreatmentSST==1` の自律収束: ω ピン留め (μ_t 上限) + wall-function P_k + k zero-gradient の
3 点が揃うと、$k$ は生産・消滅平衡 $P_k=\beta^*\rho k\omega_w$ から平衡値 $u_\tau^2/\sqrt{\beta^*}$
に収束し暴走しない。いずれか欠けると k 暴走 (zero-grad 単独) または過小 (P_k 解像勾配のまま) になる。

**適用先 (cell / node, theory.md §6.5(e))**: **P_k 置換**と **ω ピン/decouple**で適用先が分かれる。

- **cell**: `bplane_cell[ib]` = 壁隣接第一セルが唯一の近壁 DOF。`wf_pk`(P_k 置換)/ω ピン/`res_roOmega`・
  `src_jac_omega` の 0 化/`applySSTPointImplicit` の ω decouple をすべて `ic`=第一セルに当てる (現行どおり・不変)。
- **node**: 壁ノード ($y=0$) と第一内層ノードの 2 DOF がある。
  - **P_k 置換 (`wf_pk`)**: `compute_wall_friction_sst_d` が `wf_pk[ic]`(壁ノード)に加えて **`wf_pk[irep]`(第一内層
    ノード = $u_\tau$ 用に選ぶ Normal_Neighbor)** にも同じ wall-function 生産を書く。`ransSource` の P_k 置換は
    `wf_pk>=0` 判定なので両ノードに効く。壁ノードだけだと第一内層に標準解像 P_k=μ_t S² が残り近壁 k 暴走
    (`case/18.backstep` 段差コーナーで実証)。
  - **ω ピン/残差ゼロ化/decouple**: **壁ノードのみ** (現行どおり)。`rans_wall_scalar_boundary_d` が `omega[ic]`
    (壁ノード, y=`wall_y_eff`) にピン、`applySSTPointImplicit_d` が `wall_flag` で `dω=0`。`ransSource` の
    `res_roOmega`/`src_jac_omega` 0 化は **node では `wall_flag==1` 限定** (cell は `wf_pk>=0`)。`wf_pk` が第一内層
    にも付くので、ω ゼロ化を `wf_pk>=0` のままにすると第一内層ノードの solved ω まで凍結して崩壊するため、
    `wall_flag` で壁ノードに限定する。
  - **重要 (やってはいけない)**: ω ピンを第一内層ノードへ**移す** (壁アンカを外す) と、凹コーナーで複数壁が同じ
    第一内層ノードを共有してピン値が race し、壁ノード ω アンカが外れて ω 崩壊 → k 暴走 (実測で悪化、不採用)。
  - **opt-in 拡張 `nodeOmegaWfDirichlet` (2026-07-20, 既定 0)**: 壁アンカを**維持したまま**第一内層ノードにも
    ρω_w (Menter ブレンド) を**追加**ピンする (SU2 `SetTurbVars_WF` は第一点の k と ω の両方を設定する)。
    有効時は同ノードで SST の Bradshaw/strain リミッタも迂回し νt=ρk/ω (=構成上 ρν_t,wall) を使う
    (ジャンプせん断 S がリミッタ経由で νt を頭打ちにする正帰還ロックの遮断。case/38 周期チャネルで実測・
    log 則張替でも再形成される構造不安定と確認)。凹コーナー race の懸念 (上記) は壁アンカ維持により失敗形が
    異なるが未検証のため、**角の無いケース (周期チャネル等) でのみ opt-in** し、コーナー系で使う場合は
    backstep 級の再検証を先に行うこと。explicit RK3 でも状態ピンが効くよう `applyNodeKwfStatePin`
    (applyBconds 位相) が k/ω を固定する (従来の k ピンは implicit 更新カーネル内のみで explicit では
    k が初期値凍結するバグがあった — 同日修正)。
  - **残課題**: node 壁ノードは cell 第一セルより壁から遠く ω ピンが ~1/4 低いため、再付着近傍で μ_t ピークが残る
    (cell 990 に対し ~5800、局所)。場平均・x_R は cell/SU2 整合。
  - 関連 plan: [`turbulence-node-wall-function-coverage.md`](../../plans/active/turbulence-node-wall-function-coverage.md)。

**熱的閉包 `sstThermalWallFunction` (theory §6.5(f), 既定 0=OFF; 2026-08-11: mode 3 defect-flux を
正式採用 — 断熱壁 node 生産 run は `sstThermalWallFunction: 3` を指定する。以下の mode 1
[output-only] は比較用 fallback)**:
automatic wall treatment は熱側の壁法則を持たず、粗壁メッシュの断熱壁温出力が回復温度より
~200 K 冷える (case/40)。Crocco 型 $T_{aw}$ (`Taw_diag`, `ransWallFunction_d.cu` で毎ステップ計算)
を**境界出力専用のモデル壁面温度**として使う:

- 実装は `applySstThermalWallFunction` (`nodeWallDirichlet_d.cu`) の `set_wall_taw_output_d` が
  断熱壁 bplane の bvar `Ts` に `Taw_diag[ic]` を書くだけ (node/cell 共通)。境界勾配は
  owner-state-only (§6.2/§7.2.2, [`architecture-node-boundary-gradient-dof-only.md`](../../plans/accepted/architecture-node-boundary-gradient-dof-only.md))
  のため bvar `Ts` は場に効く経路を持たず、**真に出力専用** (これが成立するのは owner-state 化後。
  旧弱閉包時代は境界勾配が bvar を読んでおり「出力専用」は不正確だった — theory 参照)。
- 出力の区別: `res_wall_*` の `Ts` = モデル壁面温度 $T_{aw}$ / volume field `VALUE/T` = 壁半 CV の
  計算状態温度 $T[W]$。報告時は必ず区別する。
- W–I 内部粘性拡散は常に DOF 状態 (`T[I]-T[W]` + DOF 勾配)。断熱壁の境界半割面熱流束は厳密 0
  (`viscousFlux_wall_d` の `adiabaticWall` 分岐)。壁ノードの `res_roe` はゼロ化せず通常どおり解く。

**却下した 2 方式 (実測発散, 復活させない — theory §6.5(f) 詳述)**: ①状態ピン + res_roe ゼロ化
(過熱暴走, run_0038/0039 旧版)。②W–I compact 端点への $T_{aw}$ 注入 (全辺/代表辺限定とも、
壁半 CV の復元項喪失による異常冷却で発散: step1000 で $T[W]$ min 1100.7→27.1/15.5 K,
run_0053/0057)。**W–I 内部辺にモデル温度を混ぜてはならない**。

計画: [`turbulence-sst-thermal-wall-function.md`](../../plans/archived/turbulence-sst-thermal-wall-function.md)
(初代・弱閉包, superseded) →
[`turbulence-sst-adiabatic-taw-fluxmodel.md`](../../plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md)
(流束注入の試行と棄却の記録 + output-only fallback §0, done) →
[`turbulence-sst-su2-taw-coupling.md`](../../plans/accepted/turbulence-sst-su2-taw-coupling.md)
(SU2 式熱結合の experimental mode 2, in_progress)。

**mode 3 (`sstThermalWallFunction: 3`, defect-flux 閉包, 2026-08-11 正式採用)**: 断熱壁 × node ×
`wallTreatmentSST==1`。`compute_wall_friction_sst_d` が壁ノードに $\vec H = H_T\,\hat m$
($H_T=\rho_{rep}c_p u_\tau/T^+$, $\hat m$=内向き壁法線, 淀みは $\lambda_{rep}/y$) を
`Taw_HTnx/y/z` へ格納 (非壁 0=inactive)。`viscousFlux_d` は片端のみ $|\vec H|>0$ の W–I 辺で、
エネルギー行の (伝導+仕事) を $F=\pm(\vec S\!\cdot\!\vec H)(T_{aw}-T_W)$ に置換 (± 保存、運動量
不変)。$T_W$ は状態 `Ts[W]` を直接参照 (overlay 不要)。係数 ($u_\tau$, $T^+$, Taw) はステージ頭の
壁関数評価値で lag。mode 0/1/2 は nullptr でビット不変。theory §6.5(f) の mode 3 節参照。

**experimental mode 2 (`sstThermalWallFunction: 2`, 2026-08-11, 未採用)**: SU2 式結合。
`Taw_Prim_Overlay` (per-node overlay, `variables.hpp`) を `compute_wall_friction_sst_d`
(tawCoupled ゲート: node × 断熱壁 × mode 2) が `Taw_diag` と同値で書き、`viscousFlux_d` が
overlay 端点辺のみ SU2 corrected-gradient + 算術平均 $k_{eff}$/面速度で熱流束・粘性仕事を評価。
同ゲートで壁ノードにも `roK_wf` Dirichlet を拡張 (SU2 `SetTurbVars_WF` 両点書込相当)。
mode 0/1 は全ゲート nullptr でビット不変 (run_0061 で回帰確認)。**動的試験 (case/40
run_0059/0062/0063) で壁ノード単調冷却が止まらず現形は未採用** — 詳細は plan §10/§10b。
実験 env `FORGE_TAW_WALL_MUT_SCALE` (overlay 辺の壁側 μt スケール, 既定 1=無効) は診断専用で
恒久実装禁止。

**エネルギー壁関数 `sstEnergyWallFunction` (theory §6.5(g), 既定 0=OFF)**: 等温壁
(`wall_isothermal`) × `wallTreatmentSST==1` で、壁隣接の伝導熱流束を Kader 型 $q_w$ に
モデル置換する (冷却壁の熱負荷予測)。実装は `compute_wall_friction_sst_d`
(`ransWallFunction_d.cu`) が SST の収束済み $u_\tau, y^+, u^+$ と代表点物性から
$q_w = \rho\,c_p\,u_\tau (T_w - T_{aw,\mathrm{rep}})/T^+$ (`wallLaw_kader_tplus`, WMLES と
共有) を計算し、node は壁ノードの `Qw_Wall` に格納 → `viscousFlux_d` の AddQWall (W–I 面
xor 置換, WMLES 等温壁と同一機構) が消費、cell は bvar `qwall` → `viscousFlux_wall_d` が
壁面熱流束を $q_w S$ に置換。淀み/低 $y^+$ ($U_t\!\to\!0$ or $y^+<0.1$) は純伝導式へ退避。
状態層 (node 壁 $T$ ピン・`res_roe` 0 化・DPLUR roe decouple) は既存共有 (変更なし)。計画:
[`turbulence-sst-thermal-flux-model.md`](../../plans/active/turbulence-sst-thermal-flux-model.md)。

### 3.8 SST-DES (DDES / IDDES) length scale 修正

理論は [`theory.md`](theory.md) §8 (IDDES は §8.6)、計画は
[`turbulence-iddes-sst.md`](../../plans/active/turbulence-iddes-sst.md)。`DESmode`
(`solverConfig` の `turbulence.DESmode`, 0:RANS[既定] / 1:DDES / 2:IDDES) で
opt-in する。**`DESmode==0` では DES カーネルを一切呼ばず、`ransSource` も従来 $D_k$ 分岐を
通るため既存 SST とビット不変** (検証で確認済み)。

設定・変数・ソース対応:

- **`solverConfig.hpp/.cpp`**: `DESmode` (0..2 検証)、`C_DES_kw=0.78`, `C_DES_ke=0.61` を追加。
	`wallTreatmentSST` のパターンに倣い `getOptionalValidatedValue` で読む。
- **`variables.hpp`**: `cellValNames` に静的量 `delta_les` (Δmax) と毎 step 診断 `l_des` /
	`fd_shield` / `rd_des` を追加 (自動確保)。`output_cellValNames` にも 4 つ追加し HDF5 出力する
	(診断必須・float32 で $f_d$ の 0/1 張り付きは NaN より危険)。
- **`variables.cpp` (`setStructuralVariables_d`)**: 幾何セットアップ時に host で
	$\Delta_\mathrm{max}=\max_{隣接面}|\mathbf{cc}_{ic}-\mathbf{cc}_{jc}|$ を面ループで 1 回計算し
	`delta_les` へ H2D 転送。実行時 `ccx`/`ccy`/`ccz` (= CV 重心) のみ使い primal mesh を直接
	参照しないので node-centered でも双対 CV の Δ が得られる。
- **`turbulent_viscosity_d.cu` (新カーネル `compute_des_length_d` + device 関数
	`compute_rd_ddes`/`compute_fd_ddes`)**: $r_d$, $f_d$, $l_\mathrm{DDES}$ を計算し
	`l_des`/`fd_shield`/`rd_des` に書く。**依存順が最重要**: $r_d$ は $\nu_t=$`vis_turb`/ρ を
	読むため、`turbulent_viscosity_d_wrapper` 内で `sst_eddy_viscosity_d` (vis_turb 計算) の
	**直後**に `DESmode>0` でのみ呼ぶ。$l_\mathrm{RANS}=\sqrt{k}/(\beta^{*}\omega)$ ($\beta^{*}=0.09$,
	§8.1 の整合)。float32 ガード (grad/分母 floor + $r_d$ clamp[0,10]) を内蔵。
	カーネルは `cfg.discretization` を参照せず CV 抽象 (volume/wall_dist/勾配/vis_turb) のみに
	依存させ、node 化を検証だけで済むようにする。
- **`ransSource_d.cu` (`rans_sst_source_d`)**: 引数に `int DESmode`, `flow_float* l_des` を追加。
	`DESmode>0` で $D_k=\rho k^{3/2}/l_\mathrm{des}$, 陰解法対角 $\partial D_k/\partial(\rho k)=1.5\sqrt{k}/l_\mathrm{des}$
	に分岐。`DESmode==0` は従来 $\beta^{*}\rho k\omega$ / $\beta^{*}\omega$ をそのまま (ビット不変)。
	$\omega$ 方程式・wall-function P_k 置換・ω ピン留めは DES と独立に不変。

- **IDDES (`DESmode==2`)**: 同カーネル `compute_des_length_d` 内の分岐で theory §8.6 の
	完全関数群 ($f_B/f_{e1}/f_{e2}/f_t/f_l/f_e/f_{dt}/\tilde f_d$、
	$\Delta_\mathrm{IDDES}=\min(C_w\max(d_w,h_\mathrm{max}),h_\mathrm{max})$) を計算し
	$l_\mathrm{IDDES}$ を `l_des` に書く。**新規の幾何量・事前計算は不要**
	($\Delta_\mathrm{IDDES}$ が `wall_dist` と `delta_les` だけで組める Gritskevich 簡略形のため)。
	$r_{dt}$/$r_{dl}$ は `compute_rd_ddes` を流用 (clamp[0,10] は tanh 飽和域なので影響なし)。
	定数 $C_t=1.87, C_l=5.0, C_{dt1}=20, C_w=0.15$ は compile-time。診断: `fd_shield`←$\tilde f_d$、
	`rd_des`←$r_{dt}$ (モードで意味替え)、$f_e$ は専用フィールド `fe_iddes` (HDF5 出力)。
	$f_B$/$\Delta_\mathrm{IDDES}$ は後処理で `wall_dist`/`delta_les` から厳密再計算する。
	コード化前の関数検証は `tools/verify_iddes_functions.py` (VERDICT: PASS)。

呼び出しフロー (main.cpp): `calcGradient` → `turbulent_viscosity` (vis_turb→l_des) →
`ransTransport` → `ransSource` (l_des を読む)。

> **build 注意**: `solverConfig.hpp` に `DESmode` を追加したため、差分ビルドだと CUDA obj を
> 取りこぼし step0 で dt=0/NaN 凍結することがある (既知の罠)。config 変更後は full rebuild する。

## 4. Generic Scalar Transport 基盤

`k`, `omega` の輸送は、将来の passive scalar や species でも再利用できる
generic scalar transport 層で扱う。

実装上の責務分担は次とする。

- NS 本体: `convectiveFlux_d.*`, `viscousFlux_d.*`
- scalar 共通輸送コア (物理非依存の移流・拡散・時間積分): `scalarTransport_d.*`
- RANS 適用層 (k/ω への共通コア適用・k/ω 勾配): `ransTransport_d.*`
- model-specific boundary: `ransBoundary_d.*` など
- model-specific source / closure: `ransSource_d.*` など

`scalarTransport_d.*` は `ScalarTransportDesc` を受ける汎用 launch ヘルパ
(`scalarTransportResidual_d` / `scalarTimeIntegration_d`) のみを公開し、
k/ω 固有の descriptor 構築・有効化判定・勾配計算は `ransTransport_d.*` 側に置く。
将来 species など別の scalar を足す場合も同じ共通コアを再利用する。

初回は advection-only とし、descriptor は内部実装用の metadata として扱う。
ユーザ入力に `enabled` を追加するのではなく、`LESorRANS == 2 && RANSmodel == 1`
のとき `k`, `omega` を自動有効化する。

一方で source は core に埋め込まず、main loop で `ransTransport_d_wrapper()` の後に
model-specific source wrapper を並べる。

diffusion では `k` と `omega` に別々の係数を与え、`vis_lam` と `vis_turb`
から face 係数を組み立てる。

## 5. 軸対称 SST の幾何項

軸対称 SST の幾何項は子 plan
[`architecture-axisym-sst.md`](../../plans/accepted/architecture-axisym-sst.md)
で扱う。理論は [theory.md §7](theory.md) を参照。実装上の要点は次のとおり。

### 5.1 対流・拡散

対流・拡散は B 流儀の $r$ 重み付き面積 (`A_planar`, `sx_planar`) に幾何効果が
畳み込まれており、追加 source は不要 (theory.md §7.1)。
`ransGradient_d_wrapper` と `scalar_diffusion_first_order_d` は
すでに `A_planar` / `sx_planar` を使うため、コード変更なしで軸対称に対応する。
拡散の $\theta\theta$ 寄与が本当に発散形を保つかは、子 plan で
半径方向 1D 解析解との単体比較で確認する。

### 5.2 生産項・渦粘性のフープひずみ

`rans_sst_source_d` ([`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu))
と `sst_eddy_viscosity_d`
([`turbulent_viscosity_d.cu`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu))
は現在 planar 速度勾配 9 成分のみから $S^2$ を組む。軸対称では
フープひずみ $S_{\theta\theta} = u_r/r$ が欠落するため (theory.md §7.2)、次を追加する。

- 両 kernel に `int isAxisymmetric` と `flow_float* axisym_uy_over_r` を引数追加。
- `isAxisymmetric == 1` のとき
  $S^2 \mathrel{+}= 2\,(\texttt{axisym\_uy\_over\_r}[ic])^2$ を加える。
- `axisym_uy_over_r` は NS の
  [`axisymmetricSource_d.cu`](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)
  が毎ステップ更新済みなので、SST 側で再計算しない。呼び出し順序は
  `axisymmetricDiagnostics_d` → SST source の順を `main.cpp` で保証する。

### 5.3 圧縮性 (dilatation) 補正

生産項の圧縮性補正は `solverConfig` の `dilatationCorrection` で段階制御する
(theory.md §7.3)。

| 値 | 内容 |
| --- | --- |
| `0` | off (非圧縮形 $P_k = \mu_t S^2$) |
| `1` | (A) deviatoric のみ: $S^2 \mathrel{-}= \tfrac23(\nabla\!\cdot\!\mathbf u)^2$ |
| `2` | (A)+(B): さらに $P_k \mathrel{-}= \tfrac23\rho k(\nabla\!\cdot\!\mathbf u)$ を加え $P_k\ge0$ クリップ (**既定値**) |

既定は `2`。SST kernel は `LESorRANS==2 && RANSmodel==1` でのみ走るため、
laminar / LES ケースには影響しない。従来の非圧縮挙動が必要な場合のみ明示的に `0` を指定する。

実装は `rans_sst_source_d`
([`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu)) 内で完結する。

- 発散は `divU = dUxdx + dUydy + dUzdz`、軸対称時は
  さらに `+= axisym_uy_over_r[ic]`（= 完全発散 `axisym_divU` と一致）。
  既存配列を流用し新規計算しない。
- (A): $P_k, P_\omega$ に使う $S^2$ から $\tfrac23(\nabla\!\cdot\!\mathbf u)^2$ を引く。
- (B): $k$ 方程式の生産にのみ $-\tfrac23\rho k(\nabla\!\cdot\!\mathbf u)$ を加え、
  生産リミッタ適用後に `Pk = max(Pk, 0)`。$\omega$ 生産には入れない。
- `dilatationCorrection == 0` では従来挙動とビット一致。

検証は run_0091(off) を基準に run_0092(A) → run_0093(A+B) と段階比較する。

### 5.3.1 Kato–Launder 生産補正 (`katoLaunder`)

strain-based 生産の stagnation/加速アノマリー (theory.md §7.5) を抑えるオプション。
`solverConfig` の `turbulence.katoLaunder` で制御する。

| 値 | 内容 |
| --- | --- |
| `0` | 標準 $P_k=\mu_t S^2$ (**既定**, 既存挙動とビット一致) |
| `1` | Kato–Launder $P_k=\mu_t S\,\Omega$ ($\Omega=\sqrt{2\Omega_{ij}\Omega_{ij}}=|\boldsymbol\omega|$) |

実装も `rans_sst_source_d`
([`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu)) 内で完結:

- 渦度の大きさ二乗を既存の速度勾配から組む (新規物理量なし):
  `Om_sq = (dUxdy-dUydx)^2 + (dUxdz-dUzdx)^2 + (dUydz-dUzdy)^2` $= 2\Omega_{ij}\Omega_{ij}=|\boldsymbol\omega|^2$。
  軸対称(無旋回)でもフープは渦度を持たないため**補正不要**(ひずみと非対称)。
- `dilatationCorrection` 適用後の `S_sq` を用い、生産係数を
  $S^2 \to \sqrt{S\_sq}\cdot\sqrt{Om\_sq}$ に置換 (`katoLaunder==1` のとき)。
  $P_k=\mu_t\cdot S\Omega$、$P_\omega=\alpha\rho\cdot S\Omega$ の双方に適用。
- リミタ ($10\beta^*\rho k\omega$) と (B) 等方項は従来どおり後段で適用。
- `katoLaunder == 0` では従来挙動とビット一致。`dilatationCorrection` と直交・併用可。

**検証結果** (case 29 conical 2D, `run_0017`(off) vs `run_0018_conical_kl`(on)):
`katoLaunder:1` で軸中心線 $k$ **17.0→1.93** (μt/μlam 37→10.7), 核 $k$ **5.0→0.47** と偽生産を除去。
一方 **壁 BL は不変** (x=50mm で μt/μlam ピーク 13.0→12.9, $k$ 5374→5369), **推力も不変**
(F=1953.0→1953.1, mdot・λ 一致)。$S\approx\Omega$ の BL を保ち $S\gg\Omega$ の核生産のみ落とす狙い通り。

### 5.4 非軸対称ケースの不変性

`isAxisymmetric == 0` のときは追加項をすべて 0 とし、既存の 2D / 3D SST と
ビット単位で同一の挙動を維持する (子 plan の回帰判定基準)。

## 6. 陰解法との連成 (segregated point-implicit, 2026-06)

平均流 5 式は block DPLUR (5×5)、SST (k/ω) は **segregated point-implicit** で別に解く
（7×7 連成はしない）。block 陰解法 (`timeIntegration==11`) では平均流を commit した後、
同じ擬似時間ステップで k/ω を point-implicit 更新する（旧実装の「凍結」を解除）。

- 消散項のみ陰化（stiff、$\partial D_k/\partial(\rho k)=\beta^\*\omega$, $\partial D_\omega/\partial(\rho\omega)=2\beta\omega$）。
  `ransSource_d.cu` が残差に加えて `src_jac_k`/`src_jac_omega` を出力。
- `update_d.cu` の `applySSTPointImplicit_d` が $D_\phi=V/\Delta\tau+V\,\partial D/\partial(\rho\phi)$、
  $\delta(\rho\phi)=\text{relax}\cdot\text{res}_{\rho\phi}/D_\phi$ を作り $\rho\phi=\max(\rho\phi^N+\delta,\text{floor})$ を適用。
- `main.cpp` `implicitNonlinearUpdate` で `scalarResidualEnabled` のとき `applyBlockImplicitCorrection` 直後に呼ぶ。
- 移流・拡散・生産は残差 lagged。詳細は [`methods/time_integration/`](../time_integration/implementation.md)。

将来、強連成が必要なら `gpu-implicit-plan` の子 plan として 7×7 化を検討する。

## 7. 検証方針

標準の検証起点は [../../procedures/verification/README.md](../../procedures/verification/README.md)
に従う。

推奨順序:

1. `case/23.axi_nozzle/run_0033_slau_axisym_m4_publicrao_full_10k` の複製 run で advection-only の通りを確認
2. laminar / LES case で null regression
3. その後、airfoil 系で壁面境界層と `vis_turb` の立ち上がり確認
4. 最後に 3D case と軸対称子 plan へ展開
## 8. WMLES 代数壁応力モデルの実装

理論は [`theory.md`](theory.md) §10、計画は
[`turbulence-wmles-wall-stress.md`](../../plans/active/turbulence-wmles-wall-stress.md)。

### 8.1 ソース構成

- `cuda_forge/wallLaw_d.cuh` — 壁法則デバイス関数群 (共有ヘッダ)。Reichardt
  `wallLaw_reichardt_uplus` / `wallLaw_reichardt_duplus_dyp` (SST 壁関数 `ransWallFunction_d.cu`
  から昇格、数値経路ビット不変)、WMLES 用 warm start Newton `wallLaw_solve_utau`、Kader
  `wallLaw_kader_tplus`、`viscMethod`/`thermalMethod` 整合の壁面物性 helper。
- `cuda_forge/wmlesWallModel_d.cu` — 壁境界面並列カーネル。取得層 (cell/node の
  マッチング点・壁面状態) → 壁モデル → `bvar_d["utau"/"ypls"/"twall_*"]`、node 用
  `Tau_Wall` / `Qw_Wall` (per-CV) 書き込み。呼び出しは SST 壁関数と同位相
  (`applyBconds` 後・`viscousFlux` 前)。
- `viscousFlux_d.cu` — 壁面カーネルの `wallTreatment` を 3 値化
  (0: なし / 1: SST 壁関数 / 2: WMLES)。2 では cell の接線せん断置換 + 熱流束置換、
  node は AddTauWall / AddQWall (W↔I 再スケール)。

### 8.2 設定

- bcondConfig: `wall` / `wall_isothermal` の `ints:` に `wallModelLES: 1` (壁単位・既定 0)。
- solverConfig `turbulence` セクション: `wmlesNewtonTol` (既定 1e-6) /
  `wmlesNewtonMaxIt` (既定 20) / `wmlesPrt` (既定 0.9)。
- 有効条件: `LESorRANS != 2` (LES/ILES) かつ当該壁の `wallModelLES==1`。SST
  (`LESorRANS==2`) 経路は従来ロジックのままビット不変。
- node は `nodeWallDirichlet=1` (既定) が前提。

### 8.3 診断

- `bvar_d["utau"]` は warm start を兼ねる永続配列 (restart には含めない。初回数 step で回復)。
- Newton 不収束カウンタと $u_\tau/y^+/\tau_w$ の面統計 (min/max/mean) をログ可能にする。
- 単体検証は `tools/verify_wall_law.py` (python 参照実装との突合)。

## 整合オプション (2026-09-08, codex レビュー由来)

codex レビュー (plans/accepted/turbulence-sst-node-corner-heating.md §3.5) が指摘した標準 SST からのずれを、
`turbulence.sst*` キー (procedures/solver-settings.md) で選択できるようにした。既定は 2026-09-08 から壁 k ピン・P_ω 整合・σ ブレンドが ON (plan turbulence-sst-consistency-options §2.1)、等方応力・エネルギー k 源は OFF。

| キー | 内容 | 実装 |
|---|---|---|
| `sstNodeWallKPin` (既定 1) | node 低 Re 壁ノードの k/roK を 0 にピンし残差を 0 化。旧は ghost 反射のみで、node は境界半割面の拡散を skip するため壁 k=0 が効いていなかった | `ransBoundary_d.cu` 低 Re node 分岐, `ransSource_d.cu` 末尾 |
| `sstOmegaProdFromPk` (既定 1, 2026-09-08〜) | $P_\omega = \alpha P_k/\nu_t$ (リミッタ後の $P_k$ と整合)。0 は $P_\omega=\alpha\rho S^2$ (Menter 2003 形)。等方項 $-\tfrac23\rho k\,\nabla\!\cdot\!\mathbf u$ (dilatation 2) はリミッタ $\min(\cdot,10\beta^*\rho k\omega)$ の前に加える (キー非依存) | `ransSource_d.cu` |
| `sstSigmaBlend` (既定 1, 2026-09-08〜) | $\sigma_k = F_1\,0.85 + (1-F_1)\,1.0$, $\sigma_\omega = F_1\,0.5 + (1-F_1)\,0.856$。$F_1$ は同 step の前処理 `ransBlendF1_d_wrapper` (`ransGradient` の直後・`ransTransport` の直前) が `sstF1` に書く。残差評価順は gradient → F1 → transport(拡散) → source | `ransSource_d.cu` (`rans_sst_blend_f1_d`), `scalarTransport_d.cu` 拡散カーネル, `main.cpp` |
| `sstIsotropicStress` | 応力に $-\tfrac23\rho k\,\delta_{ij}$ (内部面のみ) | `viscousFlux_d.cu` 内部面カーネル |
| `sstEnergyKSource` | エネルギー式に $-(P_k - D_k)V$ | `ransSource_d.cu` |

検証と既定値の判断は plan `turbulence-sst-consistency-options.md`。スケール不変性の単体試験は `tools/test_scale_invariance.py`。
