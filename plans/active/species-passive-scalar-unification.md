# 受動スカラ (排気トレーサ・凝縮モーメント) を化学種輸送経路に乗せ、dual-time の物理時間項を揃える

## メタ

- **area**: `convection / diffusion / condensation / time_integration`
- **status**: `in_progress`
- **related_docs**:
  - `methods/thermophysics.md` (実装 §5b 多成分化学種輸送, §5d 拡散, §5 tracer)
  - `methods/condensation.md` (実装 §4 モーメント輸送・§4c 更新クランプ)
  - `methods/convection/theory.md` (多成分 TP の face 組成整合 S2/S3)
  - `methods/time_integration/` (dual-time, scalar-DPLUR)
- **related_plans**: [thermophysics-cea-mole-fraction-species.md](../accepted/thermophysics-cea-mole-fraction-species.md) (トレーサ `roXi` の起源, F-sp1),
  [condensation-source-limiter-steady.md](../accepted/condensation-source-limiter-steady.md) (更新クランプ `cond_moment_update_limited_d`),
  [convection-multispecies-contact-pressure.md](../accepted/convection-multispecies-contact-pressure.md) (S2/S3 face 組成; S3 は cell cfl 4 で発散し experimental 扱い),
  [condensation-followups.md](condensation-followups.md) F-cf8 (モーメントの dual-time 物理時間項) / F-cf9 (旧 mode 0 削除),
  chem ブランチの dual-time 化学種修正 (`feature/chemistry-finite-rate` e296f0d0; [[dualtime-species-frozen-bug]])
- **created**: `2026-09-17`
- **owner**: Claude (session 011spCYH)

## 1. 目的

排気トレーサ `roXi` と凝縮モーメント `rog_s, roQ0_s..roQ2_s` は汎用スカラコア (`scalarTransport_d.cu`: 1 次風上、拡散なし) で、化学種 `roY_s` は化学種輸送経路
(既定は同じ 1 次風上だが、再正規化・入口ピン・陰解法結合・Fick + 乱流拡散を持ち、`speciesFaceReconstruction 2` [S3] で 2 次面移流に切り替え可能) で解かれている。
同じ排気率でも `lumped` の $Y_{EXH}$ と `full`/tracer の $\xi$ が同一 run 内で 1.8e-4 ずれ (case/46 `run_0104`; 原因は未確定, §4.0)、
凝縮 onset は 1 次風上の数値拡散を受ける。ユーザ決定 (2026-09-16): **受動スカラは全て化学種の輸送経路を通す**。
完了時: (a) トレーサ・凝縮モーメントが「熱力学と ΣY から除外した受動種」として化学種と**同じ移流離散化・同じ更新写像・同じ拡散モデル (トレーサ)** で輸送される;
(b) 化学種と受動種の 2 次面移流 (S3) が node の本番 CFL (6 + `implicitRelax 0.7`) で安定に回り、固定点が緩和・擬似 CFL に依らない;
(c) dual-time で化学種が更新されない main のバグを chem ブランチから移植し、受動種にも BDF 物理時間項が入る (F-cf8 を閉じ、`tracer × dualTime` の拒否を解除)。

## 2. スコープ

- **やる**:
  - 化学種輸送に「受動種」の区分: 輸送 (面再構成・移流残差・拡散・陰解法対角・境界・ピン・周期・restart・出力・残差列・dual-time BDF) は化学種と同じ経路、
    **熱力学 (MW/cp/h/R, EOS, 面組成 R_mix/γ, 粘性混合)・ΣY 再正規化・ΣJ=0 補正・エンタルピー拡散・`speciesImplicitCoupling` の予測/commit・入口 X/Y 検証・
    `condGasSpecies`** からは除外。
  - トレーサ `roXi` を受動種として登録 (拡散あり: $D = \mu/(\rho\,Sc) + \mu_t/(\rho\,Sc_t)$、config `Sc`/`Sc_t` を共用、`nSpecies==1` でも動く)。
  - 凝縮モーメント 4 本 (× nCondSpecies) を受動種として移流 (拡散なし)。ソース・更新クランプ・実現可能性クランプ・EOS 結合 (g) は現行のまま。
  - S3 (2 次面移流) の本番化: 残差は常に完全な 2 次 $R_2$ (§4.2)、安定化は増分緩和と擬似 Δτ の制限に限定。node で発散原因を実装前に確認。
  - dual-time: chem e296f0d0 の化学種修正を移植し、受動種にも BDF 項。処理順 (ピン・入口・周期・restart) を明記。
  - 切替 `passiveScalarScheme` (1 = 化学種経路, 0 = 旧汎用スカラ経路 [A/B 用; 旧経路の挙動はビット不変]) と、`speciesFaceReconstruction` (0/1/2) の組で A/B。
- **やらない**: k/ω の 2 次化 (別 plan)、新しいリミッタ関数の設計 (既存 Venkat `limiter_r1_d` をスカラごとに流用)、化学反応、cell 離散化の検証 (ユーザ指示)。
  残差ブレンド型 deferred-correction (固定点を変えるので採らない; codex plan C1)。

## 3. 関連 docs と前提

- 化学種の 2 次面組成は SLAU 経路の `convectiveFlux_slau_d.inc.cuh` (S3 分岐) が `interp_dispatch` + **ψ_ρ (`limiter_ro`)** で作り、ΣY_f=1 に正規化して upwind 側を
  `Yface[ip*nSpecies+s]` に書き、`species_advection_faceY_d` が $\sum \dot m\,Y_f$ で残差を組む (`transport_diag` は 1 次風上の defect-correction)。
  `limiter_Y{s}` (Venkat on ∇Y) は診断用に存在するが、min(ψ_ρ, ψ_Y) を Y に使うと case/28 で発散したため本番は ψ_ρ。**既定 `speciesFaceReconstruction 0` では
  化学種も `scalarTransportResidualMulti_d` (1 次)** — case/44 `run_0170`、case/46 `run_0100`–`0106` は全てこの状態。
- 汎用スカラコアは 1 次風上 + 任意の μ 係数拡散 + 擬似時間 point-implicit (緩和なし)。k/ω はこのまま残す。
- 凝縮モーメントの定常固定点は残差で決まる (更新クランプは収束時に無作用) ので、移流離散化の変更 (S3) は固定点を変える。既存の凝縮回帰 (case/44 `run_0170`, case/16
  `run_0335`, Arthur) は**再取得して差を記録する**対象で、ビット一致は要求しない。
- dual-time の化学種バグ ([[dualtime-species-frozen-bug]]): main の `advanceImplicitDualTime` は化学種を更新しない (k/ω は `update_d.cu` で BDF 項あり)。修正は chem e296f0d0 のみ。
- S3 が発散した実例 `case/28 run_0052_rycl_B_s3_cfl4` は **cell** (discretization 未指定)、`implicitRelax 0.7`、`speciesImplicitCoupling 1` (scalar-DPLUR) の条件で step 433 に NaN
  (codex plan M2)。「未緩和 point-implicit の LHS 不整合」という説明は当てはまらないので、node で原因を実測してから設計する (§4.2)。

## 4. 設計方針

### 4.0 前提と仮説 (調査 2026-09-17, codex plan レビュー反映)

- 本番の化学種移流は 1 次 (§3)。本 plan の「2 次化」= S3 を化学種・受動種の両方で node 本番化すること。
- case/46 `run_0104` の同一 run 内 |ξ − Y_EXH| (最大 1.8e-4, 平均 4.6e-7) の原因は**段階比較で確定 (2026-09-17, `run_0107`/`0108`; §9)**: 移流残差
  (差 1.2e-7 / スケール 2.2e-5)・`transport_diag`・入口ピン・point-implicit の生更新 (差 4.3e-8) は float 精度で同一で、差は**化学種だけの ΣρY=ρ 再正規化
  (`species_renormalize_d`) が全て** (1.40e-4 = 全体; 副次的にトレーサだけの primitive クランプ 3.9e-5)。発生点はカウル後縁の純排気ノード (raw Y_AIR = 0, Σraw/ρ − 1 = 1.4e-4
  → Y_EXH を 1.0 に押し上げる) で、正体は block-DPLUR の δρ と 1 次 point-implicit の δ(ρY) の不整合 (限界サイクル上で残差 ≠ 0 の間は消えない)。
  よって統一経路では受動種を再正規化から除外するのが正しく、lumped の Y_EXH は残差が消えるまでこの不整合を持ち続ける (§6-1 の一致ゲートは再正規化を受けない比較種に限る)。
  `full` と `lumped` は別の流れを解くので機械精度一致は要求しない。

### 4.1 受動種の表現 (化学種経路の再利用)

- 化学種の device ポインタ配列 (`g_roY_dev` 等) と同形の**受動種配列** (roP, roPN/PM, res, transport_diag, src_jac, P, dP/dx.., limiter_P, Pface, P/PP [dual-time]) を
  `nPassive = (tracer ? 1 : 0) + 4·nCondSpecies` 本で持つ。名前は現行のまま (`roXi`/`Xi`, `rog_s`/`g_s`, `roQ2_s`, `roQ1_s`, `roQ0_s`; 順序は `registerCondensation` の
  {g, Q2, Q1, Q0}) → 出力・restart・後処理・`species_meta.yaml`・周期 gather の互換を保つ。
- **熱力学は `cfg.nSpecies` 本のまま**。受動種は `nSpecies==1` (CPG + トレーサ) でも動く (`speciesEnabled` ゲートとは独立)。
- **移流**: `SpeciesArgs` に受動種の再構成入力 (`nPassive, Pd_recon, dPdx/y/z_recon, limiterP_recon, Pface_out`) を足し、SLAU の S3 分岐と同じ `interp_dispatch` で面値を作る。
  **リミッタは受動種ごとの Venkat ψ_P** (`limiter_r1_d` を流用; 化学種の ψ_ρ 流用は「ρ と同次数」の整合が理由だが、受動種は熱力学に入らず ρ との整合は不要で、
  スカラ自身の単調性の方が重要)。**スケーリング (codex plan-2 M1)**: 現行 `limiter_r1_d` は $\Delta_+^2\Delta_-$ 等の 3 次積を float32 で組むため、モーメント
  ($Q_0 \sim 10^{15}$) では Inf/NaN になる。受動種のリミッタは差分をセル局所スケール $\phi_{ref} = \max(|\phi_c|, \max_{nb}|\phi|, \phi_{floor})$ で無次元化し
  ($\tilde\Delta = \Delta/\phi_{ref}$, $\epsilon^2$ も無次元) Venkat 関数は同じ式で評価する (スケール不変なので ψ は変わらない)。単体で g/Q2/Q1/Q0 の実測スケール・
  ゼロ近傍・極値でクリップ前の係数が有限であることを検査する。正規化はせず、面値を下限 0 (トレーサは [0,1]) でクリップ。**面クリップは保存性を壊さない** (同一面の 1 つの流束を両 CV に逆符号で
  加える; codex M3 訂正)。upwind 側を `Pface_out[ip*nPassive+q]` に書き、`species_advection_faceY_d` を受動種ポインタで呼ぶ (`res_<cons>`, 1 次風上 `transport_diag_<prim>`)。
  SLAU 以外、または `speciesFaceReconstruction < 2` のときは受動種も化学種と同じ 1 次経路 (`scalarTransportResidualMulti_d`) に落ちる (常に化学種と同じ次数)。
- **有界性**: 面値の非負・[0,1] クリップは更新後の有界性を保証しない (codex M3 の反例: ステップ状分布で流出面値 > 流入面値 → 負の増分)。有界性は
  (i) ψ_P リミッタで再構成の overshoot を抑え、(ii) 更新後の floor (`scalarTimeIntegration_d` の `max(·, floor)` / モーメントの実現可能性クランプ) で確定し、
  (iii) **floor による保存量の補正量を診断に載せる** (`passiveFloorCorr_<name>`: |Δ(ρφ)| の累積、`condClampCorr` と同じ規約)。試験は `Xi` の表示値でなく生の
  `roXi/ro`・総トレーサ量 ∫ρξ dV・floor 補正量で判定し、モーメントは非負に加え相互の実現可能性 ($q_1/q_0 \le r_{30}$ 等) を確認する。
- **拡散**: トレーサのみ、Fick 形 $D = \mu/(\rho Sc) + \mu_t/(\rho Sc_t)$ (定数 Sc; 粘性 run のみ) を ΣJ=0 補正・エンタルピー項なしで加える `passive_diffusion_d`。
  化学種との一致は**等拡散係数 (定数 Sc、`speciesDiffusionMethod 0`) の条件に限定**し、混合平均 (差動) 拡散では非一致が正しい挙動 (codex M8)。モーメントは拡散なし。
- **境界・ピン**: 化学種の Dirichlet/Neumann/ピン (`species_dirichlet_boundary_d` 系, `species_pin_residual_d`) を受動種ポインタで呼ぶ (トレーサ入口 `Xi`、モーメント入口 0)。
  ピン残差除去は BDF 項追加の**後**にも適用する (§4.4)。
- **周期 node (codex M4)**: 名前ベースの `res_*` gather だけでは不足。処理順を固定する:
  1. 勾配積算は周期半割面を除外して行い (`species_gradient_d` に除外引数を追加)、勾配を gather (既存の周期勾配 gather に受動種・化学種の勾配を登録);
  2. 空間残差 `res_*` と **`transport_diag_*`** を gather (ソース Jacobian `src_jac_*` は体積比で整合);
  3. dual-time では gather 後に合併体積で BDF 項を**一度だけ**追加;
  4. 受動種の更新後に状態 mirror (`periodicMirrorScalarState` に受動種を登録; 平均流更新時の mirror では不十分)。
  5. **化学種の coupling 1/2 も周期整合させる (codex plan-2 M3)**: coupling 2 の予測・EOS クロス項注入は「擬似刻み確定後、BDF 込み・ピン除去済み残差」で行い、
     クロス項は**独立バッファ**に組んでその追加分だけを周期 gather してから流れの RHS に加える (gather 済み残差への直接加算は部分 CV 分、再 gather は重複)。
     coupling 1 (scalar-DPLUR) の近傍補正と各 sweep の `dq` は周期グループで整合させ (流れ側 `main.cpp` の周期補正同期に相当)、化学種の更新後状態も mirror する。
  試験: 周期 node 箱の移流・拡散 (継ぎ目・辺・角で保存量 ∫ρφ dV と更新速度が内部と一致) を受動種だけでなく**非一様組成の化学種 (coupling 0/1/2)** で行い、
  一様過飽和の凝縮ソース (case/09 `run_0064` プロトコル) も再確認する。
- **更新**: 受動種は現行どおり segregated point-implicit (`scalarTimeIntegration_d` / `cond_moment_update_limited_d`)。化学種の結合予測/commit には入れない。
  **上下限 (codex plan-2 M2)**: 更新確定時に**更新済み密度**で $0 \le \rho\xi \le \rho$ (トレーサ)、$\rho\phi \ge 0$ (モーメント) を適用する (流れ更新で ρ が減ると
  非負化だけでは ξ>1 になる)。primitive 段で保存量を書き換えるクランプは撤廃。**補正収支**: 上限・下限それぞれの符号付き補正 $\int \Delta(\rho\phi)\,dV$ と絶対補正
  $\int |\Delta(\rho\phi)|\,dV$ を物理 step 内 (擬似反復の和) と全期間の積算で記録し (`passiveFloorCorr_<name>` は更新ごとの正規化量でなく積算; `condClampCorr` とは
  規約が違うことを明記)、保存誤差は $|\Delta \int\rho\phi\,dV| / \int\rho\phi\,dV$ の相対値で定義、補正量にも合否閾値 (全期間の相対積算 ≤ 1e-4) を置く。
- 切替 `passiveScalarScheme` (0 = 旧汎用スカラ経路 [ビット不変], 1 = 化学種経路)。既定は検証完了後に 1。

### 4.2 S3 の本番化 — 完全な 2 次残差 + 増分緩和 (codex plan C1/M2 反映)

- **残差は常に完全な $R_2$** (面再構成した流束の総和)。残差のブレンドや deferred-correction 係数は使わない (固定点が変わる)。
- **原因確認の結果 (2026-09-17, node case/28 `run_0064`–`0078`; §9)**: S3 の発散は **`speciesImplicitCoupling 0` (segregated point-implicit, LHS = 1 次風上 `transport_diag`,
  緩和なし) に固有**で、組成せん断層で化学種残差が 50–80 step で成長し cfl 4 で step 402 / cfl 6 で step 165 に NaN。**S3 + coupling 1 (scalar-DPLUR [近傍結合の陰解法
  sweep], relax 0.7) は cfl 2/4/6 で S3 固有の発散が無く** (化学種残差は ×0.7 に減少)、baseline の別問題 (下記) で死ぬまで走る。cell の `run_0052` (S3+c1 で発散) とは逆で、
  §10 の contingency は反転する: **脆弱なのは受動種が使う予定だった segregated 更新の方**。
- したがって安定化は次の順で入れる: (i) **受動種も化学種の scalar-DPLUR sweep (`species_dplur_sweep_d`) に乗せる** (`passiveImplicitCoupling 1`, S3 時の既定;
  受動種ポインタ配列でそのまま呼べる。coupling 2 の EOS 結合には入れない); (ii) segregated 更新 (`passiveImplicitCoupling 0`) は `passiveImplicitRelax`
  (既定 = `implicitRelax`) で増分を緩和 (`passiveScalarScheme 1` のみ; 旧経路 0 は緩和なし = ビット不変); **化学種の segregated 更新の緩和は独立キー `speciesImplicitRelax`
  (既定 1.0 = 現行と同じ写像)** (codex plan-2 M5)。無影響試験には `implicitRelax 0.7` の run を含める。(iii) 保険として化学種/受動種だけ擬似 Δτ を絞る `scalarCflMax`
  (既定なし; 物理時間項は変えない)。S3 の本番推奨は「coupling 1 + relax 0.7」を基本にし、coupling 0 + S3 は緩和/上限付きでのみ許す。
- **実装前の原因確認 (node)**: case/28 の S3 発散は cell・`speciesImplicitCoupling 1`・relax 0.7 の条件だった。node の同一 IC/BC で S2/S3 × coupling 0/1 の 4 組を cfl 4 で回し、
  発散の有無・最初の NaN の位置・更新方式との対応を記録してから設計を確定する (発散が coupling 1 [scalar-DPLUR] 固有なら受動種 [coupling 0 相当] は無関係)。
- 固定点不変の検証 (codex plan-2 M4): `run_0471` は組成がほぼ一様 (max−min 3e-8〜7e-7) で組成再構成の寄与が丸め程度なので**無影響回帰にのみ使う**。固定点ゲートには
  **定常でも非ゼロ勾配が残る**ケースを新設する: case/16 `run_0471` プロトコルに `inletProfile` CSV で半径方向に変わる $Y_{H2O}$ (例 0.02→0.06) とトレーサ $\xi$ (0→1) を与えた
  node Euler run (`run_0473` 系; PASS が取れる 12000 step)。$R_2 - R_1$ が反復ノイズより十分大きいことを確認し、`speciesImplicitRelax`/`passiveImplicitRelax` 0.7/1.0、
  `scalarCflMax` 有/無、cfl 2/6 の**交差 restart** で全残差 PASS・収束解が反復ノイズ内・補正無作用 (floor 補正 0) を要求する。固定点不変は「正の緩和率・非特異な更新作用素・
  制約補正が無作用」の条件付きなので、補正無作用を必ず併記する。

### 4.3 凝縮モーメントの移流

- `condensationTransport_d_wrapper` の `scalarTransportResidualMulti_d` を受動種の移流 (4.1) に置換。`res_<cons>`・`transport_diag_<prim>` の名前・ゼロ初期化・順序
  (移流 → ソース → 更新) は現行のまま。ソース・更新クランプ・実現可能性クランプ・EOS 結合は不変。
- 面値の負値クリップは保存的 (4.1)。更新後の floor 補正は `condClampCorr` に既に載る。

### 4.4 dual-time (chem e296f0d0 の移植 + 受動種; codex M6 反映)

- 移植する: `species_shift_levels_d` (`roY{s}P/PP`)、`species_add_unsteady_d` (BDF 残差 $-(V/\Delta t)(a\,\rho Y - b\,\rho Y^P + c\,\rho Y^{PP})$ + 対角 $Va/\Delta t$;
  履歴カウンタはシフト後に増えるので **最初の 1 物理 step が BDF1**)、`g_roPred_dev` (plan-C の δρ 基準 = 予測時点の ρ)、更新後の再正規化・primitive。
  持ち込まない: CMC の再正規化、`CMC7_TFREEZE_BYPASS`、`jacobianMode 2` block LU、`dt_local_sp`。保つ: float32 拡散、node ピン、fused 残差、S3 `Yface`。
- **処理順 (サブ反復ごと)**: assembleResidual (空間残差 + ピン残差除去) → 周期 gather → BDF 項追加 (化学種・受動種; 合併体積で一度だけ) → **ピン残差を再除去**
  (BDF がピン行に残差を戻すため) → 流れ block 更新 → 化学種更新 (coupling 0: point-implicit / 1: scalar-DPLUR [移植元に無い分岐を本ブランチの定常経路から流用] /
  2: EOS 結合 commit) → 再正規化・primitive → 入口 Dirichlet の再適用 → 受動種更新 (BDF 込みの `res`/`transport_diag` で `cond_moment_update_limited_d` /
  `scalarTimeIntegration_d`) → 状態 mirror。
- **履歴契約 (codex plan-2 M6)**: BDF の履歴有効数と係数は流れ・化学種・受動種で**共有する 1 つの状態** (`cfg`/StepContext の `nHistoryValid` と `iStep`; 移植元の
  プロセス内カウンタ `g_speciesLevelShifts` は使わない)。checkpoint は流れ・化学種・受動種の履歴 (`*N/NN`, `*P/PP`)・物理時刻・刻み・履歴有効数をまとめて書き、
  restart は全部揃っているときだけ復元し、1 つでも欠ければ**全系を揃えて BDF1 から**再開する (旧形式 `res_*.h5` はこの経路)。起動時は `P = PP = 現在値`。
- `condLimiterMode 1` の dual-time 自動降格と `tracer × dualTime` 拒否を解除 (検証後)。`scalarCflMax` は物理時間項に触れない。

### 4.5 k/ω

対象外 (§10)。

### 4.6 モーメント更新クランプの非負化は成分ごと (codex result-2 M1; 2026-09-17)

- #13 で入れた「共通 θ による非負化」(いずれかの成分で $N_k=0$ かつ $d_k<0$ なら $\theta=0$ で **4 本全部の増分を停止**) は**撤回**する。
  `run_0257` では 663 ノード (全て $Q_1=0$、蒸気枯渇ではない) で g/Q0 の更新まで止まり、θ≥1e-12 の保証を壊してモーメントの sub-iter 残差床と 1 次相当の
  誤差 (BDF2 で 1.3 次) を作っていた。
- 新しい規則: 共通 θ は $\theta_u$ (閾値 [初回 sub-iter のみ] + 実現可能性 [avail, 蒸発上限]; ≥1e-12) のまま。非負化は**成分ごと**に、
  $N_k+\theta_u d_k<0$ となる成分 $k$ だけ増分を $-N_k$ に切って確定値 0 にする (他の成分は $\theta_u d_k$ のまま)。切った量 $\theta_u|d_k|-N_k$ を
  `passiveLimCorr_<k>` (セル累積) と limStats (絶対量・作動数・符号付き; root のみ) に記録する。`condLim` は $\theta_u$。
- 「停止セルの寄与」の記録: 非負化の作動数と絶対量は monitor の受動種収支行 (成分ごと) に出る。時間次数試験 (§6-6) では floor+lim の収支が総量比 ≤1e-6
  であることを次数評価の前提にする (作動していれば次数の言い訳にせず原因を潰す)。
- 期待: $Q_1=0$ のノードで g/Q0 の更新が止まらなくなり、モーメントの sub-iter 床が流れと同水準 (≥2 桁) に下がる。次数は #12/#18 の 3 水準 × nSub 倍増で再評価し、
  BDF2 2.0±0.3・nSub 倍増差/最小水準差 ≤0.1 を全量 (g, Q0, T, ro) で要求する。

### 4.7 非定常 (dual-time) の受動種の有界化 — 物理 step 末尾の保存的 FCT 補正 (Zalesak 型; codex result-2 M2, plan-3 M1/M2/m2 反映; 2026-09-17)

- **問題**: CV ごとの増分スケーリング $\theta_b$ (#14) は非保存 (縮めた分を隣へ戻さない)。seam ステップ試験 `run_0098` で総量 0.54 % 損失。
- **plan-3 で棄却した案**: sub-iter 内の残差レベルで「陽的低次予測 $\rho\phi^k+\Delta t R_{low}/(aV)$」に Zalesak を掛ける形は、予測が CFL>1 で有界でなく
  (BDF1 初回反復で $[1,0,0,0]\to[-1,2,0,0]$)、収束時も遷移域で $Q^\pm$ が退化して制約が無意味になる。
- **適用範囲**: `passiveScalarScheme 1` かつ S3 (SLAU, `speciesFaceReconstruction ≥ 2`) かつ **dual-time (`timeIntegration 11`, `unsteady 1`, `dualTime 1`)**。
  定常 (局所 Δτ) には掛けない (Δτ 依存の固定点を作る; F-cf7 と同じ罠)。定常は現行の ψ_P 再構成 + 起動時の $\theta_b$ 増分緩和のまま。ただし
  「固定点で無作用」は保証でなく**検証条件**: limCorr の停止と未制限残差 (θ_b が縮めた候補増分の総量) の消滅を別々に確認する。陽的 RK は対象外
  (Chaplin–Colella の RK4 積分流束制限のような先行例はあるが受動種の非定常本番は dual-time; §10)。
- **定式化** (物理 step 末尾、sub-iter 収束後に 1 回。面 $f=(i,j)$、$\dot m>0$ は $i\to j$、$a$ は BDF 係数 [BDF1 1, BDF2 3/2]、hist $=(b\,q^n-c\,q^{n-1})/a$):
  1. **HO 解** $q_H$ = sub-iter が収束した現行の完全陰的 2 次スキームの状態 ($\rho\phi^{n+1}$; 変更なし)。
  2. **低次陰解** $q_L$: 同じ BDF 履歴・同じ質量流束 (最終 sub-iter の `massflux`)・同じソース (最終 sub-iter の $S(q^k)V$ を凍結)・同じ拡散で、移流だけ 1 次風上にした線形系
     $$\Big[\tfrac{aV}{\Delta t}+\textstyle\sum_{out}\tfrac{\dot m}{\rho_i}+D_{ii}\Big]q_{L,i}-\sum_{in}\tfrac{\dot m}{\rho_j}q_{L,j}-\sum_j D_{ij}q_{L,j}=\tfrac{V}{\Delta t}(b\,q^n_i-c\,q^{n-1}_i)+S_iV$$
     を Jacobi sweep で解く (M 行列: 反復ごとに有界 [RHS ≥ 0 なら $q_L\ge0$、上限は流入側の凸結合]; 周期 group は近傍和を gather; ピン行は $q_L=q_H$;
     `passiveFctSweeps` 上限 [既定 100] と相対変化 `passiveFctTol` [既定 1e-6] で停止し、最終の線形残差ノルムを診断に出す)。BE (初回 step) では厳密に有界。
     BDF2 の履歴項 ($-c\,q^{n-1}$) は負の係数なので $q_L$ が僅かに負になり得る (plan-3 の反例: CFL 4 のステップで $-1.8\times10^{-2}$) — これは「保証の穴」として
     残し、最後の floor の収支に記録して試験で定量化する (§6-2)。ソースが負 (蒸発) のときも同様。
  3. **反拡散流束** (面): $A'_f=\dot m_f\,(P_{face,f}-\phi_{L,up})+[J_f(q_H)-J_f(q_L)]$ ($P_{face}$ は最終 sub-iter の ψ_P 再構成面値、$J_f$ はトレーサ Fick 拡散流束)。
     2 つの離散式の差から厳密に $\tfrac{aV}{\Delta t}(q_H-q_L)=\sum_f A'_f$ (ソース・履歴は相殺; 最終 sub-iter の状態と収束状態の差 = sub-iter 誤差だけ残る)。
     内部面 (ip < nNormalPlanes) 以外は $A'_f=0$ (node 境界半割面は自セル値、周期半割面は移流ループに無い)。
  4. **Zalesak**: 前制限 ($A'_f(\phi_L(j)-\phi_L(i))<0\Rightarrow A'_f=0$; `passiveFctPrelimit` 既定 1) → $P^+_i=\sum\max(0,\text{流入})$, $P^-_i=\sum\max(0,\text{流出})$ →
     限界 $\phi_{max/min}$ = {自身, 面隣接} の {$\phi^n=\rho\phi^n/\rho^n$, $\phi_L=q_L/\rho^{n+1}$} の max/min (トレーサはさらに $[0,1]$、モーメントは下限 0) →
     $Q^\pm_i=(\rho^{n+1}_i\phi_{max}-q_{L,i},\ q_{L,i}-\rho^{n+1}_i\phi_{min})\,aV/\Delta t$, $R^\pm=\min(1,Q^\pm/P^\pm)$ ($P=0$ で 1、$Q<0$ で 0) →
     $\alpha_f=\min(R^+_{受け側},R^-_{出し側})$ (受け側は $A'_f$ の符号で決める)。周期 node: $P^\pm$ は和 gather、極値は max/min gather (§4.8 の wrapper)。
  5. **補正**: $q^{n+1}:=q_H-\tfrac{\Delta t}{aV}\sum_f(1-\alpha_f)A'_f$ ($=q_L+\tfrac{\Delta t}{aV}\sum_f\alpha_fA'_f$ を $q_L$ の解誤差に依らず $q_H$ 側から書いた形)。
     面ごとに 1 つの $\alpha_f$ を両 CV に逆符号で加えるので**保存**、$\alpha_f=1$ の面では $q_H$ のまま (時間 2 次の完全陰的スキームは不変)、
     制限が作動した所だけ有界な $q_L$ 側へ寄る。周期 node は補正を独立バッファに組んで gather (和 → broadcast) してから加える。
  6. **後処理の契約**: 入口 Dirichlet の再適用 (ピンへの交換量は境界収支として `passiveFctPinCorr` に記録) → floor (`passive_bounds_d`; `floorCorr` に記録) →
     モーメントの実現可能性クランプ (`condClampCorr`) → primitive → 周期 mirror。$\phi_N\delta\rho$ 項は sub-iter 内の話で、ここでは密度は確定済み。
     `passiveImplicitCoupling` 0/1 のどちらでも同じ補正 (sub-iter の解法に依らない)。
- **診断**: `fctCorr` = $\sum_f(1-\alpha_f)|A'_f|\Delta t/a$ (落とした反拡散の絶対量; 保存量の損失ではない)、作動面数、$q_L$ の線形残差、
  $q^{n+1}$ の限界逸脱量 (= floor が処理した量) を monitor の受動種収支行に併記する。
- **合否は用途で分ける**: 非定常 (dual-time) は「総量保存 ≤1e-6 (float; 周期箱は root-only 総量、非周期は境界流束・ソース込みの収支)」+「`floorCorr + limCorr + fctPinCorr`
  の総量比 ≤1e-6」+「$0\le\xi\le1$, モーメント ≥0」を要求 (1e-12 は double 参照実装の代数試験のみ)。定常は起動時の limCorr を許し固定点で 0。
- キー `passiveFct` (`time.deltaT`; 1 = 上の条件で有効 [既定], 0 = 無効 [A/B]), `passiveFctPrelimit`, `passiveFctSweeps`, `passiveFctTol`。

### 4.8 周期 node のリミッタは group で統合 (codex result-2 M3; 2026-09-17)

- `limiter_r1_scaled_d` (受動種 ψ_P) と `limiter_r1_d` (流れ・化学種 ψ) は各 member の**部分 CV** の内部面だけで極値・スケール・ψ を作るため、周期対で ψ が
  一致しない (`run_0089` の 4532/2982: 0.985 vs 1.0)。状態 mirror が一致しても再構成の作用素が合併 CV と違う。
- 修正: 2 段に分ける。(1) 極値 $Q_{max}/Q_{min}$ (自身 + 内部面隣接) を配列へ → 周期 group で max/min gather (broadcast); (2) 合併極値 (受動種は φ_ref も合併極値から)
  で各 member の面の ψ を評価 → group で min gather (broadcast)。非周期 (root なし) では従来の 1 段 kernel のまま (ビット不変)。流れ・化学種の `limiter_r1_d` も同じ
  2 段化を周期 node で行う (受動種だけ直しても組成で ρ が変わる試験では流れ側の ψ_ρ の不一致が残る)。
- 試験: SLAU + `convMethod 1` + `limiter 2` + `speciesFaceReconstruction 2` (受動種 S3) で `run_0090`/`0091` プロトコル (dual-time coupling 1, ガウス ξ + ガウス組成 [N2,O2])
  の π シフト等価 (内部 vs seam 中心) を追加: max|Δ| が float 床、周期対の `limiter_Xi`/`limiter_ro` の差 = 0。

## 5. 実装ステップ

1. `methods/thermophysics.md` §5b/§5d (受動種、拡散の適用条件)、`methods/condensation.md` 実装 §4 (モーメント移流)、`methods/convection/theory.md` (S3 本番化の安定化)、
   `methods/time_integration/` (dual-time の化学種/受動種 BDF と処理順) を先行更新。
2. **原因確認 (実装前)**: node で S3 の発散再現 (§4.2) と、`run_0104` 差の段階比較 (§4.0) を行い結果を §9 に記録。
3. forge 受動種基盤: `passiveInit_d` (配列・ポインタ)、`SpeciesArgs` 拡張と SLAU S3 分岐、`species_advection_faceY_d`/`species_gradient_d` (周期除外引数)/境界/ピンの
   受動種呼び出し、ψ_P リミッタ、`passive_diffusion_d`、floor 診断、`condensationTransport_d`・`tracerTransport_d` の移流置換、周期 gather/mirror の登録、`passiveScalarScheme`。
4. S3 本番化: 新経路の増分緩和 (`implicitRelax`)、`scalarCflMax`。
5. dual-time: chem 移植 + 受動種 BDF + 処理順 (§4.4) + 拒否/降格の解除。
6. 検証 §6 (node のみ) → docs 同期。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~調査 (コード地図) と §4 の具体化~~ | 済 (2026-09-17): §4.0–4.5 |
| 2 | ~~codex plan レビュー~~ | 1 回目 NO-GO (C1/M7/m1)、2 回目 **GO-with-changes (M6)** を全採用 (§6.1)。実装着手可 (2026-09-17) |
| 3 | ~~docs 先行更新~~ | 済 (2026-09-17, a8745508): thermophysics §5b/§5, condensation §4, convection/theory, time_integration/theory |
| 4 | ~~原因確認 (node)~~ | 済 (2026-09-17): S3 発散は coupling 0 固有 (case/28 node `run_0064`–`0078`)、1.8e-4 は再正規化が全て (case/46 `run_0107`/`0108`) → §4.0/§4.2 に反映。副産物: node TP-SST の case/28 baseline 自体が上境界/軸で ~1500–3000 step 後に発散する別問題 (§10) |
| 5 | ~~受動種基盤 (ステップ 3) + S3 安定化キー (ステップ 4)~~ | 済 (2026-09-17, Phase A; §9): `passiveTransport_d.cuh`/`passiveKernels_d.cuh`/`passiveLimiter_d.cuh` (受動種ポインタ配列、化学種カーネル共用、無次元化 Venkat、`passive_bounds_d` 収支、`passive_diffusion_d`)、SLAU S3 分岐の受動種面値、周期 (勾配除外・`transport_diag` gather・DPLUR dq mirror・coupling 2 クロス項の独立バッファ)、キー `passiveScalarScheme`/`passiveImplicitRelax`/`speciesImplicitRelax`/`scalarCflMax`/`passiveImplicitCoupling`、単体 `test_passive_scalar.cu` ALL PASS。**残**: `speciesImplicitRelax`/`scalarCflMax` の効果 run、S3 長時間 run での上限クランプ無作用の確認 |
| 6 | ~~S3 本番化 (ステップ 4)~~ | 済 (2026-09-17; §9 検証): 非一様組成 PASS ケース case/16 `run_0476` (S3, PASS 全列 3.3–4.7 桁, Y_H2O 0.020–0.060, ξ 1e-4–1.0 の勾配を保持, R2−R1 はノイズの 1850 倍) からの交差 restart (cfl 6 + relax 0.7 `run_0478`, relax 1.0 `run_0479`, scalarCflMax 2 `run_0480`) が反復ノイズ内・floor 補正 0。S3 の本番推奨組合せ = SFR 2 + `speciesImplicitCoupling 1` + pc1 + relax 0.7 (solver-settings に記載)。**既定は SFR 0 のまま** (S3 は凝縮の固定点を動かす; §10) |
| 7 | ~~dual-time 移植 + 受動種 BDF (ステップ 5)~~ | 済 (2026-09-17, Phase B; §9)。残件 3 点も済: coupling 1/2 (case/16 `run_0491`/`0492`: ΣY 7e-8、sub-iter 2.4 桁、coupling 0 との差 1e-5 級 = 反復ノイズ)、ピン行 (case/16 `run_0493`, case/44 `run_0228`: BDF + 再ピン後の res = 0.0、入口値差 0.0)、過渡 IC 凝縮 dual-time (case/44 `run_0226`/`0227`: g/Q2/Q1 ≥2 桁、**Q0 1.4–1.8 桁・流れ/化学種 0.95–1.9 桁は nSub 10 では未達** → nSub 20 の再 run は §6-6 の残)。|
| 9 | ~~(result-1 M1) 周期 DPLUR の近傍寄与を独立バッファで周期群合算~~ | 済 (2026-09-17): `dplurSweepOnce` (近傍寄与 `species_dplur_neighbor_d` → 周期 gather → `species_dplur_solve_d`; 化学種 coupling 1/2 と受動種 DPLUR で共通、非周期はビット不変)。case/09 dual-time 陰解法 (coupling 1 + pc1) の seam 版 `run_0091` vs 内部版 `run_0090`: 周期対 Δ=0、π シフト等価 max|Δξ| 6.5e-7 / |ΔY0| 4.5e-7 (float 床)、sub-iter 4.6–6.0 桁。**残**: coupling 2 の周期試験は CPG 受動組成の箱では案C (NASA 熱力学前提) が成立せず NaN → TP の静止周期箱 (u=0, 組成ガウス, 層流拡散) で再試験 (#17) |
| 10 | ~~(M2) 保存量・補正収支の積分を周期 root のみ合併体積で~~ | 済 (2026-09-17): `passive_bounds_d` に root 引数; `run_0088` のログ総量 148.8301 = 後処理 root-only (差 2e-7; 旧 0067/0068 の保存数は後処理由来で不変); seam のステップ ξ + S3 (`run_0098`) で root-only 総量変化 −4.5376e-1 = 符号付き limCorr −4.5380e-1 (恒等式の不一致 4e-5 = 総量比 5e-7; RK は確定ステージのみ記録) |
| 11 | ~~(M3) checkpoint: dt・生成方式を layout に入れ整合しなければ全系 BDF1~~ | 済 (2026-09-17): layout に dt/bdfOrder/passiveScalarScheme/speciesImplicitCoupling、dt 相対 1e-12 一致が復元条件、scheme 0 は受動種履歴を書かない。非定常 (case/09 ガウス ξ+組成 BDF2): 連続 `run_0094` vs 通常 restart `0095` (差 5e-7 = 復元), dt 変更 `0096` / 履歴削除 `0097` は `history NOT restored` → BDF1 (差 3e-4), scheme 変更 case/44 `run_0233`→`0234`→`0235` は layout mismatch で BDF1 |
| 12 | (M4) 凝縮 dual-time: sub-iter 収束と モーメント込み 3 水準次数 | 調査済 (2026-09-17, case/44 `run_0236`–`0253`): sub-iter の床は擬似 CFL・緩和・pc1・リミッタでなく**流れの 2 次 MUSCL 再構成** (convMethod 0 で全列 median ≥2.2 桁)。convMethod 0 の 3 水準: BDF1 は g/Q0 0.89, T/ro 2.1; BDF2 は T 2.43 / ro 3.27 だが **g/Q0 の次数は未達** (水準差 1e-7 が sub-iter 床に埋没; 大刻み系でも nSub 倍増差が最小水準差の 0.4 倍) → 原因は発達域の更新クランプ θ_u (condLim min ≈0.1) が sub-iter ごとに増分を絞り収束を止めること。**設計決定 (2026-09-17)**: dual-time では θ_u (dg_max/dT_max) を「物理 step あたりの上限」として**各物理 step の初回 sub-iter だけ**に掛け、2 回目以降は θ_u=1 (実現可能性・非負は維持) → 再試験 (#18) |
| 13 | (M5) 受動種の上下限を増分スケーリングで守る | 部分済 (2026-09-17): `passive_limit_increment_d` (θ_b) + モーメントの共通 θ 非負化 (`boundByTheta`), floor と limCorr を別記録。Arthur S3 `run_0108`: floor 積算 ≤1e-20 (旧 5.6e-3; limCorr は起動過渡のみ 7.4e-3)。トレーサ `run_0498` 1.76e-4 / `run_0499` 4.84e-4 は**未達**: 残る floor は全て上限側で、流れ block 更新で ρ が減った後に ρξ_N が ρ_new を超える = 分離更新の δρ 不整合 (θ_b では消えない)。**設計決定 (2026-09-17)**: 受動種の更新に ξ_N·δρ 項を入れる (案C と同じ: δ(ρξ) = z + ξ_N δρ、z は輸送増分; モーメントも同様に φ_N δρ) → 流れ更新と整合し、輸送増分 0 なら ξ 不変。再検証 (#19) |
| 14 | ~~(M6) `implicitRelax` 0.7/1.0 の交差 restart、収束場 restart の判定ツール、case/44 README の表現~~ | 済 (2026-09-17): 緩和キーの実効経路を確認 (coupling 0 → `speciesImplicitRelax`, coupling 1 → `implicitRelax`, 受動 point-implicit/DPLUR → `passiveImplicitRelax`); case/16 `run_0497` (`implicitRelax 1.0` + `passiveImplicitRelax 1.0` の交差 restart) は `run_0476` と反復ノイズ内 (Y0 1.9e-7, ξ 1.5e-5, ro 1.3e-6, floor 補正 0); `check_convergence.py --from-floor REF` (参照 run の末尾 20 % 床に対し peak/tail ≤1.5×) を追加し `run_0478`/`0479`/`0480`/`0497` は **PASS (within 1.5x of reference floor)**; case/44 README の表現を準定常量の比較に限定 |
| 15 | ~~(m1) solver-settings / recommended-settings / thermophysics の旧記述を統一~~ | 済 (2026-09-17, a6892e13) |
| 16 | ~~(m2) §5.1・§10 (ψ 切替)・case/09/16 README の番号同期~~ | 済 (2026-09-17): README の改番前番号を修正 (case/09 時間次数 run_0080–0086, case/16 dual-time run_0486–0490)、全引用 run の実在を確認; ψ_P vs ψ_ρ は §10 で「ψ_P を採用 (無次元化 Venkat; 固定点ケースで補正 0)」として決着 |
| 17 | ~~coupling 2 の周期試験を TP 周期箱で~~ | 済 (2026-09-17): TP [N2, H2O] 周期箱 (組成ガウス + ξ, SLAU 1 次, dual-time, `nStepInner 4` [1 sweep は周期箱で NaN], pc1) **u=10 で coupling 2** (`run_0109` 内部 / `0110` seam): 周期対 Δ=0、π シフト等価 Y1 6.1e-7 / ρ 8.3e-7 / ξ 9.7e-6 (float 床)、∫ρY 6e-7、sub-iter 3–4.8 桁。u=0 では EOS クロス項が作動しない (mdot=0) ので u=10 系列を根拠にする (u=0 系列 `run_0105`–`0108`/`0111`–`0114` は記録のみ) |
| 18 | dual-time の θ_u を物理 step あたりに; 凝縮の 3 水準次数 (g/Q0) | 実装済・**次数は未達** (2026-09-17, case/44 `run_0256`–`0262`; convMethod 0, nSub 40/80, dt 1.6e-5/8e-6/4e-6): sub-iter 低下は流れ・化学種 ≥2.1 桁、rog/Q2/Q1 ≥2.4 桁だが **roQ0 min 1.64** (nSub 80 でも同値; condLim min 0 = 実現可能性クランプ [avail/非負] が常時作動するセルで増分 0); 観測次数 BDF2 g 1.30 / Q0 1.36 / T 1.27 / ro 2.29、BDF1 g 0.95 / Q0 0.70 / T 1.31 / ro 1.93; nSub 倍増差/最小水準差 g 0.43 / Q0 0.37 / T 0.07。**結論**: 流れ・トレーサ・化学種の BDF2 (2.0–2.3) は成立、モーメントは硬い実現可能性制約が作動するセルで時間 1 次相当 (制約整合の時間離散が要る = 本 plan の範囲外; §10) |
| 19 | ~~受動種更新の φ_N·δρ 項; floor ≤1e-4~~ | 済 (2026-09-17): `passive_add_rho_term_d` (ρφ += φ_N δρ, 制限は輸送増分のみ; 単体 z=0 で φ 不変 6e-8)。floor 積算: case/16 `run_0500` (run_0476 型) **3.4e-26** (旧 1.76e-4), `run_0501` (run_0494 型) **0** (旧 4.84e-4), Arthur S3 `run_0109` ≤4e-21 (limCorr は起動過渡のみ 7.2e-3 で per-step 0); 固定点 `0500` vs `run_0476` は反復ノイズ内 (ξ 1.9e-5); scheme 0 無影響 (case/44 `run_0254`/`0255` ノイズ床) |
| 20 | (result-2 M1) モーメントの共通 θ 停止を撤回 (成分ごとの非負化, θ≥1e-12)、停止セルの残差寄与の記録、次数の再評価 (#12/#18 の実体) | forge (必須): 非負化のコードは済 (`cond_moment_update_limited_body`; 切った量 = その成分の未適用残差を limStats に記録)、次数の再評価は #21 の後に S3 でも (§6-6) |
| 21 | (M2 → plan-3 で再設計 §4.7) 物理 step 末尾の保存的 FCT 補正 (低次陰解 $q_L$ を限界に $q_H$ から落とす); dual-time seam ステップで総量保存 1e-6; floor/lim/fct 併記 | forge (必須; 単体 `test_passive_fct.cu` + case/09 dual-time seam + §6-6 の S3 次数) |
| 22 | ~~(M3) ψ_P の周期群での極値・スケール統合と係数共有; SLAU + SFR 2 の seam/内部平行移動試験~~ | 済 (2026-09-17): `limiterPeriodic_d.cuh` (極値 → max/min gather → ψ → min gather; fused5 [流れ] / `limiter_r1_d` [Y] / `limiter_r1_scaled_d` [受動種] の全呼出し), `periodicAtomic_d.cuh` (CAS 比較の float atomicMax/Min, NaN member は no-op), 単体 `test_periodic_limiter.cu` PASS (群 2/4/8 混合符号; 分割鎖と内部鎖でビット一致); case/09 `run_0115`/`0116` (SLAU convMethod 1 + limiter 2 + SFR 2, dual-time coupling 1): 周期対 `limiter_Xi` 差 2.5e-2/1.07e-1 → **0**, π シフト等価 limiter_Xi 1.04e-1 → 7.4e-6, ξ 6.3e-7 (不変); u=−10 (`run_0119`/`0120`) も全 ψ の周期対差 0。非周期 (case/46 `run_0114`–`0116`) は同一バイナリ反復ノイズ内 |
| 23 | (m1) `--from-floor`: 参照の通常 PASS・列対応・ゼロ参照列の非ゼロ化を検査 | tools |
| 24 | (m2) methods/condensation.md §4c (θ, φ_N·δρ, limCorr, dual-time の初回 sub-iter 限定 = 物理 step 上限ではない)、§4/§10/#8、plans/README の同期 | docs |
| 8 | 検証 (§6 1–8, node のみ) と codex result レビュー | §6-1/2/7 (短 run) Phase A; §6-3 拡散 (解析解 0.23 %, 2 次収束; 等拡散一致 [コア 2.7e-6] / 混合平均非一致 [3.7e-3 一定] 済), §6-4 固定点, §6-5 凝縮回帰 (S3 は固定点を動かす), §6-7 24000 step は済 (§9)。**§6-6 のモーメント時間次数と §6-2 の非定常保存は未達** (#20–#22 が必須残作業)。既定 `passiveScalarScheme 1` (2026-09-17)。codex result 2 回 NO-GO |

## 6. 検証 (node のみ; cell はユーザ指示で対象外)

全比較で `speciesFaceReconstruction`・`speciesImplicitCoupling`・`implicitRelax`・実効スカラ CFL (`scalarCflMax`) を明記して固定する。

1. **差の原因の段階比較** (`run_0104` プロトコル, lumped [EXH, AIR] + tracer, node Euler): 凍結流れ (`FORGE_FREEZE_*` 相当または収束場から 1 step) で
   残差 → 更新前後 → 再正規化前後 → クランプ前後の各段で |ξ − Y_EXH| を記録し、差の由来を確定 (**済 2026-09-17**: 再正規化が全て, §9)。**改定ゲート (Phase A 実測)**:
   残差・対角・生更新は float 精度で同一; ξ と再正規化前の化学種 raw 更新の差は「上限クランプ (raw>ρ) が無作用のノード」で平均 ≤1e-6 (max はクランプ履歴の伝播で
   1e-5 級を許容)、クランプ作用ノードは補正収支で説明できること。**厳密一致ゲート (max ≤ 1e-6) は全離散作用素と更新写像を揃えた
   制御試験に限定** (1 次移流 [`speciesFaceReconstruction 0`]・拡散なし・再正規化なし・同じ緩和・同じ floor; codex plan-2 M5)。S3 では化学種 (ψ_ρ + 面正規化) と
   受動種 (ψ_P + 独立クリップ) のリミッタが違うので一致は要求せず差を記録する。
2. **保存・有界**: 一様流中のステップ状 ξ の移流 (node 箱 + 周期, 1 次/S3): ∫ρξ dV の保存 1e-6、floor 補正量の記録、0 ≤ roXi/ρ ≤ 1 (floor 後)、
   S3 が 1 次より鋭いこと。モーメントは一様過飽和の凝縮ソース試験で非負・実現可能性・継ぎ目/辺/角の更新速度が内部と一致 (case/09 `run_0064` プロトコル)。
   **追加 (plan-3 M4/m1, 2026-09-17)**: 非定常の保存は float32 で 1e-6 (double 集計の 1e-12 は `test_passive_fct.cu` の double 参照実装のみ)。dual-time の seam ステップ (`run_0091` プロトコル + `run_0098` のステップ ξ, SLAU S3, FCT 有効) で root-only 総量保存 ≤1e-6、`floorCorr+limCorr+fctPinCorr` の総量比 ≤1e-6、$0\le\xi\le1$; 非周期入口 (case/16 のトレーサ分布 run を dual-time で短時間) は境界流束・ソース込みの収支; 周期の辺・角 (33³ 箱の角を跨ぐガウス) でも同じ収支。定常の θ_b は「limCorr の停止」と「未制限残差 (縮めた候補増分の総量) の消滅」を別々に確認する。
3. **拡散 (F-sp1 を閉じる)**: 解析解付き node 拡散試験 (1 次元ガウス核の拡散; 層流 `Sc`、乱流 `Sc_t` [`vis_turb` を人為的に与える]、無流束壁、入口、周期) で
   誤差が 2 次収束; `nSpecies==1` + トレーサでも動作; モーメントに拡散が入らない; 等拡散係数条件で化学種と一致、混合平均拡散では非一致 (記録のみ)。
4. **S3 本番化 (安定性・固定点)**: node で S2/S3 × coupling 0/1 を cfl 4 (原因確認, ステップ 2); **非一様組成 + トレーサの PASS ケース** (§4.2: `run_0471` プロトコル +
   半径方向 $Y_{H2O}$/ξ 分布, `run_0473` 系) + S3 で `speciesImplicitRelax`/`passiveImplicitRelax` 0.7/1.0、`scalarCflMax` 有/無、cfl 2/6 の交差 restart が全て **PASS**、
   収束解が反復ノイズ内、floor 補正 0 (固定点不変); `run_0471` は無影響回帰; case/44 `run_0170` プロトコルを S3 + cfl 6 + relax 0.7 で完走 (NaN 0, series STEADY) —
   こちらは準定常回帰 (固定点の証拠ではない)。
5. **凝縮回帰 (準定常, 差の記録)**: case/44 `run_0170` プロトコルで `passiveScalarScheme` 0/1 × `speciesFaceReconstruction` 0/2 の onset・出口 g・g max・series STEADY・
   condLim 1・補正 0・floor 補正量; case/16 `run_0335` プロトコル (Wysłouzil) で onset と実験の差の変化; Arthur N2 node (`case/34 run_0106` プロトコル) の onset。
   同一バイナリ反復ノイズ床を併記。
6. **dual-time**: (i) 多成分 dual-time で ΣY=1 と化学種の時間発展 — 「`roY` が動く」だけでは不十分で、**非一様組成の化学種を含めて** 3 水準の時間次数で判定する;
   (ii) 時間精度: 滑らかな解 (ガウス状 ξ・モーメント・組成の移流 + ソース作動時の凝縮) を同一最終時刻で Δt, Δt/2, Δt/4 の 3 水準で比較し、3 水準の差から求めた次数が
   BDF2 で 2.0±0.3 (BDF1 では 1.0±0.3)、サブ反復数倍増で結果の差 ≤ 3 水準最小差の 1/10、各物理 step のサブ反復で全化学種・受動種の BDF 込み残差が ≥2 桁下がること
   (物理 step の `outer_end` だけで判定しない); (iii) coupling 0/1/2 の各分岐でピン行の残差が 0、入口値が保たれる; (iv) **restart**: 連続実行 vs 途中 checkpoint からの
   restart が反復ノイズ内で一致、旧形式入力 (履歴なし) は BDF1 再開で完走、`nSpecies==1` + 受動種でも同じ (codex plan-2 M6)。
7. **無影響**: `passiveScalarScheme 0` + `speciesFaceReconstruction 0` で現行とビット一致 (case/44 `run_0170`, case/46 `run_0100`, **`implicitRelax 0.7` の run** [case/44
   `run_0183`])、化学種のみの run (case/16 `run_0471`) は `passiveScalarScheme` に依らず不変、`speciesImplicitRelax` 省略時 (1.0) は化学種の更新写像が現行と同一。
8. **判定基準**: 上のゲート + `check_convergence.py` / `check_quasisteady.py` VERDICT、NaN 0、step 時間の増分を記録。既定変更 (`passiveScalarScheme 1`,
   `speciesFaceReconstruction 2` の推奨) は 1–7 を全て通した後。
   **追加 (plan-3 M4, 2026-09-17)**: 時間次数試験は **FCT が作動する SLAU + `speciesFaceReconstruction 2`** でも行う (case/44 `run_0256` プロトコルを S3 で; BDF1/BDF2 × `passiveImplicitCoupling` 0/1, 3 水準 + nSub 倍増)。評価量に g, Q0 に加え **Q1, Q2**、T, ρ; 制限作動集合 (α<1, 非負化作動, floor 作動) の未制限残差を記録し、作動集合を除いた次数と全域の次数を併記。実現可能性不等式 $Q_1^2\le Q_0Q_2$, $Q_2^2\le Q_1Q_3$ ($Q_3\propto g$) の違反ノード数を最終場で数える (後処理 `analyze_moment_order.py`)。保存は最終 BDF 残差の積分から見積もる誤差 ≤1e-6。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result (2 回目) | `2026-09-17` | [2026-09-16-species-passive-scalar-unification-result-2.md](../../notes/reviews/2026-09-16-species-passive-scalar-unification-result-2.md) | **NO-GO**, C0/M3/m2 | **全採用 (2026-09-17 反映中, §5.1 #20–#24)**: M1 (モーメント次数未達の真因は #13 で入れた共通 θ 非負化: いずれかのモーメントが N_k=0 かつ d_k<0 なら **全モーメントの増分を θ=0 で停止** [`run_0257` で 663 ノード, Q1=0, 蒸気枯渇ではない, θ≥1e-12 の保証を無効化]; T の次数 1.27・ro の sub-iter 比 0.17・最小刻みの rms_roY1 1.9 桁も未達) → 共通停止を撤回し成分ごとに負値だけを非負化 (θ≥1e-12 維持)、停止セルの残差寄与を記録してから次数を再評価 (#12/#18 は必須のまま); M2 (CV ごとの増分制限は非保存: seam ステップ試験で総量 0.54 % 損失、limCorr が floor から移っただけ) → 非定常の受動種は**面流束の共有制限** (Zalesak 型: 2 次補正流束を両 CV で同じ係数で縮小) で有界化し総量保存 1e-6 を満たす; 定常の起動緩和と非定常の保存を分けて判定、floor と lim を併記; M3 (S3 受動種リミッタ ψ_P が member の部分 CV ステンシルで計算され周期点で不一致 [0.985 vs 1.0]) → 周期群で極値・スケールを統合し係数の最小値を共有; SLAU + SFR 2 の seam/内部平行移動試験を追加; m1 (`--from-floor` が参照の通常 PASS・列対応・ゼロ参照列の非ゼロ化を検査しない) → 検査を追加し不成立は拒否; m2 (condensation.md の θ≥1e-12 記述、φ_N·δρ / limCorr / 初回 sub-iter 限定の未反映、「物理 step あたりの上限」と実装の差、#8 の完了状態、plans/README の draft) → 同期 |
| result (1 回目) | `2026-09-17` | [2026-09-16-species-passive-scalar-unification-result.md](../../notes/reviews/2026-09-16-species-passive-scalar-unification-result.md) | **NO-GO**, C0/M6/m2 | **全採用 (2026-09-17 反映中, §5.1 #9–#16)**: M1 (周期 DPLUR が合併 CV の近傍補正を解かず root の dq をコピーするだけ) → sweep の近傍寄与を独立バッファに作り周期群で合算してから解く + 非一様組成/トレーサが seam・辺・角を横切る陰解法試験 (#9); M2 (周期で補正収支・総量が全メンバーで積分され重複; `run_0067` 総量 166.2 vs root のみ 148.8) → 積分は周期 root のみ合併体積で; seam 局在補正の収支試験 (#10); M3 (checkpoint の履歴有効判定不足: dt 変更で警告のみ、layout に scheme 無し) → dt・生成方式の整合を復元条件に、不一致は全系 BDF1; 非定常の組成/トレーサ/モーメントで連続・通常 restart・刻み変更・方式変更・履歴欠落を試験 (#11); M4 (凝縮 dual-time の sub-iter 桁低下 <2 [Q0 1.58, roY0 1.92, roY1 1.77]、sub-iter 依存 g 3.7e-5 が刻み半減差 9.5e-5 に対し小さくない、モーメントの 3 水準無し) → sub-iter の収束を改善 (擬似 CFL・モーメントの DPLUR 更新) しモーメント込み 3 水準の次数試験; §10 送りにせず §5.1 必須へ (#12); M5 (累積 floor 補正が上限 1e-4 超: `run_0476` ξ 2.2e-4、`run_0494` 5.0e-4、Arthur S3 Q0 5.6e-3) → 上下限を事後クリップでなく増分スケーリング (θ_u 型) で守り floor は最後の保険にする; 再検証 (#13); M6 (緩和比較が実効でない: coupling 1 は `implicitRelax` を使う; 収束場 restart の判定ツール無し; case/44 の「固定点は同一」表現) → `implicitRelax` 0.7/1.0 の交差 restart、収束場 restart の判定 (基準 PASS の床と量の変動幅) をツール化、case/44 は準定常量の比較に表現限定 (#14); m1 (solver-settings/recommended-settings/thermophysics の旧仕様) → 同期 (#15); m2 (§5.1 の済/未達併記、§10 ψ 切替未決、case/09/16 README の改番前番号) → 同期 (#16) |
| plan (3 回目) | `2026-09-17` | [2026-09-16-species-passive-scalar-unification-plan-3.md](../../notes/reviews/2026-09-16-species-passive-scalar-unification-plan-3.md) | **NO-GO**, C0/M4/m2 (§4.6–4.8 の設計) | **全採用 (2026-09-17 反映)**: M1 (sub-iter 内の陽的低次予測は CFL>1 で有界でない; BDF2 の負係数と負ソースで低次厳密解も負になり得る) → §4.7 を「物理 step 末尾の保存的 FCT 補正: 実際に解いた低次陰解 $q_L$ を限界に、$q^{n+1}=q_H-\frac{\Delta t}{aV}\sum(1-\alpha)A'$」に再設計、BDF2/負ソースの穴は floor 収支で定量化; M2 (処理順: ピン復活・境界半割面の stale Pface・φ_N δρ との順序・coupling 0/1) → §4.7-3/6 の契約 (内部面のみ A'、補正後に再ピンと境界収支、密度確定後なので δρ 項なし、両 coupling 共通); M3 (流れ側は `limiter_r1_fused5_d`、順序保存 int 化は負値に使えない) → §4.8 に fused5 を明記、CAS 比較の float atomicMax/Min (全符号・±0・NaN 規約) + 単体試験 (群 2/4/8, 負速度 run); M4 (試験が FCT の失敗を検出しない) → §6-6 に FCT 有効 (SLAU+SFR2) の 3 水準 × nSub 倍増 (BDF1/2, coupling 0/1)、Q1/Q2 の誤差、実現可能性不等式 $Q_1^2\le Q_0Q_2$、圧縮流・非周期入口・周期辺角の収支、BDF 残差積分の保存誤差 ≤1e-6 を追加; m1 (float32 に 1e-12 は無理) → 1e-6 + 丸め測定、1e-12 は double 参照のみ; m2 (前制限は P± の前、受け側は A の符号) → 反映 |
| plan (2 回目) | `2026-09-17` | [2026-09-16-species-passive-scalar-unification-plan-2.md](../../notes/reviews/2026-09-16-species-passive-scalar-unification-plan-2.md) | **GO-with-changes**, C0/M6/m0 | **全採用 (2026-09-17 反映)**: M1 (Venkat の float32 3 次積がモーメント 1e15 で Inf/NaN) → 受動種のリミッタはセル局所スケールで無次元化した差分で評価 + 有限性単体 (§4.1); M2 (非負化だけでは ξ≤1 を保証しない、`condClampCorr` は積算でない) → 更新確定時に更新済み ρ で 0≤ρξ≤ρ、符号付き/絶対補正の体積積分を step 内・全期間で積算、相対保存誤差と補正閾値 (§4.1); M3 (周期で coupling 1/2 の予測・EOS クロス項・dq 同期が未定義) → クロス項は独立バッファで gather、DPLUR dq と更新後状態の周期整合、非一様組成 coupling 0/1/2 の周期試験 (§4.1-5); M4 (`run_0471` は組成ほぼ一様で固定点ゲートが空振り) → 非一様 Y_H2O/ξ 分布の PASS ケース `run_0473` 系を新設し交差 restart + 補正無作用 (§4.2, §6-4); M5 (化学種の緩和を受動種の切替に紐付けると無影響ゲートと矛盾、S3 のリミッタ差で 1e-6 一致は不成立) → `speciesImplicitRelax` (既定 1.0) を独立キーに、無影響試験に relax 0.7 を含め、厳密一致は全作用素を揃えた制御試験に限定 (§4.2, §6-1, §6-7); M6 (BDF 履歴の共有・checkpoint・restart 契約と時間次数の数値基準が無い) → 履歴有効数と係数を全系で共有、checkpoint 一括復元/欠落時 BDF1、連続 vs restart・旧形式・nSpecies==1 の試験、次数 2.0±0.3 等の数値基準 (§4.4, §6-6)。**実装着手可** (原因確認は並行) |
| plan | `2026-09-17` | [2026-09-16-species-passive-scalar-unification-plan.md](../../notes/reviews/2026-09-16-species-passive-scalar-unification-plan.md) | **NO-GO**, C1/M7/m1 | **全採用 (2026-09-17 反映)**: C1 (残差ブレンド deferred-correction は固定点を変える) → 撤回、残差は常に完全な $R_2$、緩和は増分と擬似 Δτ のみ (§4.2); M2 (case/28 の発散例は cell・relax 0.7・coupling 1 で「未緩和 point-implicit」ではない) → node で S2/S3 × coupling 0/1 の原因確認を実装前に (§4.2, §6-4); M3 (面クリップは保存的、有界性は面値では保証されない、`Xi` は輸送に読む原始量) → 保存性の記述を訂正、ψ_P リミッタ + 更新後 floor + floor 補正量の診断、生の roXi/ρ・総量・実現可能性で判定 (§4.1, §6-2); M4 (node 周期は名前ベース gather では不足) → 勾配の周期除外と gather、`transport_diag` gather、更新後 mirror、BDF は gather 後に一度、周期試験必須 (§4.1); M5 (1.8e-4 の原因は未確定) → 仮説に戻し段階比較で確定、一致ゲートは同じ更新方式の制御試験に限定 (§4.0, §6-1); M6 (dual-time の処理順・ピン・coupling 分岐・時間精度ゲート) → §4.4 の処理順、3 水準の時間精度、サブ反復残差、ソース作動試験 (§6-6); M7 (完了条件から S3 安定性・原因切り分けが抜け、`passiveScalarScheme 1` だけでは 1 次のまま) → 検証 1–7 を完了条件に、A/B は SFR/coupling/relax/scalar CFL を固定、固定点不変は PASS ケースで (§6, §8); M8 (F-sp1 を閉じる拡散検証が無い) → 解析解付き拡散試験 (§6-3); m1 (§1/§3 の「2 次」「limiter_Y」記述、`implicitRelax` 既定 1.0、BDF1 は最初の 1 step、緩和は新経路のみ) → 本文修正 |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/{speciesTransport_d.cu,cuh, convection/convectiveFlux_d.cu, convection/convectiveFlux_slau_d.inc.cuh, convection/convectiveFlux_common_d.cuh,
  limiter_d.cu, condensationTransport_d.cu, tracerTransport_d.*, periodicNode_d.cu, scalarTransport_d.cu, update_d.cu}`, `variables.cpp`, `main.cpp` (定常・RK・dual-time),
  `input/solverConfig.*`。
- 既存 run: `passiveScalarScheme` 省略時の既定を 1 にすると (再正規化・クランプの扱いが変わり) トレーサ・モーメントの結果が変わる; S3 推奨化で凝縮 run の固定点が変わる。
  回帰 run の再取得と README 記録が必要。旧経路 0 はビット不変。
- docs: `methods/thermophysics.md`, `methods/condensation.md`, `methods/convection/theory.md`, `methods/time_integration/`, `procedures/solver-settings.md`,
  `procedures/recommended-settings.md`。

## 8. 完了条件

- [ ] 関連 methods を更新済み
- [ ] 実装・§6 の検証 **1–7** (原因切り分け / 保存・有界 / 拡散 / S3 安定性と固定点不変 / 凝縮回帰 / dual-time / 無影響) を満たす
- [ ] codex レビュー (`plan` 2 回目 GO 以降 / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status: done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`、`plans/README.md` 同期、followups F-cf8 / F-sp1 を閉じる

## 9. 変更ログ

- `2026-09-17` — codex plan レビュー 3 回目 **NO-GO (M4/m2)**: §4.7 (sub-iter 内 FCT) を物理 step 末尾の保存的補正に再設計 (§6.1)。§4.8 (周期リミッタの 2 段化) は実装済 (#22)。#20 の成分ごと非負化を実装。
- `2026-09-17` — codex result レビュー 2 回目 **NO-GO (M3/m2)** を全採用 (§6.1, §5.1 #20–#24): モーメント次数未達の真因は #13 の共通 θ 停止 (バグ)、増分制限の非保存、ψ_P の周期不整合。
- `2026-09-17` — #17–#19: φ_N·δρ 項で受動種の更新を流れの密度更新と整合 (floor 補正が 3 case とも ≤1e-4、固定点不変)、TP 周期箱 u=10 で coupling 2 の周期整合を float 床まで確認、dual-time の θ_u を物理 step 初回 sub-iter のみに。**モーメントの BDF2 次数は未達 (1.3)**: 残る sub-iter 床 (roQ0 1.64 桁) は実現可能性クランプが常時作動するセル (condLim 0) 由来で nSub に依らず、流れ・化学種・トレーサ (2.0–2.3) と切り分けた。§6-6 のモーメント次数ゲートは本 plan では満たせない → §10 に理由と後続 (制約整合の時間離散) を記載し、codex result 2 回目で採否を諮る。
- `2026-09-17` — result-1 M1/M2/M3/M5 の solver 修正と M4 調査 (§5.1 #9–#13): 周期 DPLUR の近傍寄与を周期群で合算 (π シフト等価 6.5e-7)、収支の root 化、checkpoint の dt/方式整合 (不一致は全系 BDF1)、増分スケーリング θ_b (Arthur S3 の floor 5.6e-3 → ≤1e-20)。未達 2 件に設計決定: (M5) トレーサの上限側 floor (1.8e-4/4.8e-4) は分離更新の δρ 不整合 → ξ_N·δρ 項を導入; (M4) モーメントの BDF2 次数は θ_u が sub-iter 収束を止めるため → dual-time では θ_u を物理 step の初回 sub-iter のみに。sub-iter 床の主因は流れの 2 次 MUSCL (convMethod 0 で全列 ≥2.2 桁)。
- `2026-09-17` — result-1 M6/m2 対応: 実効緩和の交差 restart `run_0497` (反復ノイズ内)、`check_convergence.py --from-floor` (収束場 restart の判定; `run_0478`/`0479`/`0480`/`0497` PASS)、README 番号・表現の同期。
- `2026-09-17` — codex result レビュー 1 回目 **NO-GO (M6/m2)** を全採用 (§6.1, §5.1 #9–#16)。
- `2026-09-17` — 既定 `passiveScalarScheme 1` で `build-passive` と既定 `build/` をフルリビルド (HEAD 33dc6588 のソルバソースと一致)。単体: `test_passive_scalar` / `test_cond_limiter_steady` / `test_cond_float_device` / `test_species_db_host` (31) / `test_solver_config_species` (15; `tracer × dualTime` は scheme 0 明示時のみエラーに試験を更新) / `test_convert_species_field_fail` (21) / `test_tpgas_lowT` (5) / `test_interp_field_tracer` (3) 全 PASS。既定キー確認 case/44 `run_0231` (キー無し) vs `run_0232` (明示 1): 差は node 反復ノイズ内 (ro 2.7e-6, g 9.6e-5)。**codex result レビューへ**。
- `2026-09-17` — §6-3 等拡散一致/混合平均非一致 (case/16 `run_0494`/`0495`/`0496`, node NS 層流 [MIXDRY, H2O] + トレーサ, 入口ステップ Y_H2O 0.05 ↔ ξ 1, 定数 Sc 0.7 共有, SFR 0, coupling 0, relax 1.0): 定数 Sc では |Y_H2O/0.05 − ξ| が残差とともに単調減少 (24000 step: max 4.8e-4 / 平均 4.8e-5 / コア 2.3e-5 → 72000 step: max 4.6e-5 [入口角の壁列] / 平均 4.7e-6 / **コア 2.7e-6**)、両スカラの残差はスケール後 0.4 % で一致 (rms_roXi/20 = rms_roY1) = 移流・拡散・更新写像が同一で、残る差は化学種だけの再正規化の過渡分; `run_0494` PASS (全列 3.2–4.0 桁)。混合平均拡散 (`speciesDiffusionMethod 1`, `run_0495`): 差は界面で 3.66e-3 に**一定** (減衰しない) = 差動拡散による正しい非一致。
- `2026-09-17` — **検証 §6-3/4/5/7 (node, スナップショット binary)**: (§6-4) 非一様組成 + トレーサの PASS ケース case/16 `run_0476` (`run_0471` プロトコル + inletProfile で Y_H2O 0.06 軸→0.02 壁, ξ 1→0; S3 + coupling 1 + relax 0.7, 24000 step): `check_convergence` PASS 全列 (ro 4.3, roUy 3.3, roY 4.2–4.3, roXi 4.6 桁), Y_H2O 0.020–0.060 / ξ 1.2e-4–0.99983 の勾配を保持, floor 補正は過渡のみ (累積 2.2e-4, step 5600 以降 0); SFR 0 双子 `run_0477` との差 Y0 3.7e-4 / ξ 9.4e-3 = 反復ノイズ (`run_0481`: Y0 2.0e-7, ξ 1.5e-5) の 1850×/630×; 交差 restart `run_0478` (cfl 6 + relax 0.7) / `0479` (relax 1.0) / `0480` (scalarCflMax 2) は `run_0476` と反復ノイズ内 (Y0 ≤2.7e-7, ξ ≤1.5e-5, ro rel ≤3.0e-6), floor 補正累積 0.000 (verdict は収束場からの restart のため drop 基準で NOT CONVERGED; 床 4.7e-9 で不変)。(§6-3) case/09 周期箱の RK4 ガウス拡散 `run_0071`–`0074` (33³/64³): 分散増加が解析 2Dt に対し −0.93 % / −0.23 % (比 4.1 = 2 次), L∞ 6.4e-3→1.5e-3, 保存 1e-7; SST の乱流拡散でトレーサ界面が層流の 2 倍に広がる (case/16 `run_0484`/`0485`); モーメントは拡散なし (Wys の g 壁法線分布が S3/scheme 0 で同一)。**注意**: 受動種/化学種の拡散は `viscMethod 0` では加わらない (定数粘性でなく拡散なし扱い; solver-settings に明記)。(§6-5) case/44 `run_0170` プロトコル 24000 step: scheme 0 `run_0220` / scheme 1 `run_0221` (SFR 0) は run_0170 と g_0 7e-6 / g L1 3e-8 (床内), **S3 `run_0222` (cfl 2) / `run_0223` (cfl 6 + relax 0.7)**: onset 16.024 r_t (+0.18 r_t 下流), x30 18.497 (−0.15), 出口 g 0.5639 % (−3.5 %), 前線が鋭くなり出口近くの径方向前線で T 差 18.6 K; 0222 と 0223 は同一 (g_0 3.5e-6) = S3 固定点は cfl/relax 非依存; 全 run ALL STEADY・condLim 1.0000・補正 0。Wysłouzil (`run_0470/lim1e` プロトコル 48000 step): scheme 0 `run_0483` = lim1e (床内), **S3 `run_0482`: onset 23.23 mm (+0.72 mm 下流; 実験からはさらに遠い)**, 壁偏差 +2.47/−5.0/+5.5/+5.5 %, PASS (全列 3.1–6.9 桁), ALL STEADY。Arthur N2 node (`run_0106` ラベル `passA_s0` / `passA_s1_sfr2`): s0 は 28/28 PASS、S3 は 15/28 (固定点変化: g_0 4.3e-2, T 3.4e-3), Δp/p onset 2.10 in は両方同一, g>1e-4 onset 2.01 vs 2.03 in。(§6-7) `run_0220` vs run_0170 24000 step: 床内 (T max 0.0165 K は壁 16 ノード, 反復値の 2.2 倍)。step 時間: scheme 1 SFR 0 は scheme 0 より速い (2.37 vs 2.91 ms), S3 2.8 ms, Wys S3 3.77 vs 2.91。**決定**: `passiveScalarScheme` の既定を 1 に (1 次経路は scheme 0 と床内一致, dual-time 対応); `speciesFaceReconstruction` の既定は 0 のまま (S3 は凝縮 onset を下流へ動かし実験比較を変える; 使うときは coupling 1 + relax 0.7)。
- `2026-09-17` — 過渡 IC 凝縮 dual-time の sub-iter 数依存 (case/44 `run_0229` nSub 20 / `run_0230` nSub 40): BDF 込み sub-iter の桁低下は nSub 10→20→40 で**頭打ち** (rog/Q2/Q1 ≥2.2 桁、Q0 1.5–1.7、ro/roY0 ≈1.9、roe 1.5、roUy 1.1) = 擬似時間反復の残差床 (cfl_pseudo 2 の point-implicit と発達域の θ_u クランプ [condLim min 0.10]) であり sub-iter 数では解消しない。同一物理時刻の場の差: nSub 10 vs 20 で ro 3.8e-3 (Δt 半減差 3.6e-3 と同程度 = nSub 10 は不足)、nSub 20 vs 40 で ro 1.9e-4 / g 3.7e-5 / Q0 4.1e-4 (モーメントは反復床)。→ この case の時間精度評価は nSub ≥20 が必要; 「全列 ≥2 桁」は擬似 CFL・更新方式側の課題として §10 に残す。
- `2026-09-17` — Phase B 残件: dual-time × coupling 1/2 (case/16 `run_0491`/`0492`) は coupling 0 と反復ノイズ内、ピン行残差・入口値は 0.0 (env `FORGE_PIN_DIAG=1` の host 読み戻し診断; case/16 `run_0493`, case/44 `run_0228`)、凝縮発達中の過渡 IC からの dual-time (case/44 `run_0226` dt 8e-6 / `0227` dt 4e-6; IC = `run_0186/res_4000`): NaN 0、g/Q2/Q1 の sub-iter 桁低下 ≥2.2、Q0 1.4–1.8・流れ/化学種 0.95–1.9 (nSub 10 では 2 桁未達; dt 半減で 1.3→1.8)、Δt vs Δt/2 の同一物理時刻差 ro 3.6e-3 / g 9.5e-5 (過渡変化量の ~5 %)。
- `2026-09-17` — **Phase B 実装 (dual-time)**: chem e296f0d0 を移植 (`g_roPred_dev` 予測時点 ρ、`speciesShift/AddUnsteadyTimeTerm` [BDF 残差 + 対角 Va/Δt; 係数は呼び出し側]、CMC/block LU/`dt_local_sp` は除外)、受動種版 `passiveShift/AddUnsteadyTimeTerm` (scheme 1 のみ; scheme 0 は不変)、`roY{s}P/PP`・`roXiP/PP`・`<cons>P/PP` 登録、共有履歴 `nHistoryValid` (BDF2 は `bdfOrder>=2 && nHistoryValid>=1`)、`advanceImplicitDualTime` を §4.4 の順に (BDF → ピン再除去 → 残差記録 → 予測/クロス項 → block → commit → SST → 化学種 [coupling 0 relax / 1 DPLUR / 2 案C] → 再正規化 → ミラー → primitive → 入口再適用 → 受動種更新 → 入口再適用)、dual-time の DPLUR 非対角 ρ は現在値、`/CHECKPOINT` group (layout・totalTime・dt・nHistoryValid + 全履歴) を `res_*.h5` に書き、全部揃うときだけ復元 (欠ければ P=PP=現在値・BDF1)。`condLimiterMode 1` の降格と `tracer × dualTime` 拒否は scheme 0 のみ。**検証 (node)**: (i) 多成分 smoke case/16 `run_0486` (収束場, ΣY 6.2e-8) / `run_0490` (過渡 IC: BDF 込み sub-iter 残差が全列 median 2.4 桁低下, ΣY 7.0e-8); (ii) **時間精度** case/09 `run_0080`–`0086` (周期箱 KEEP + 移流ゲージ, ガウス ξ と CPG 受動組成 Y0, 最終時刻 0.2 s, dt 4e-3/2e-3/1e-3): **BDF2 観測次数 ξ 2.03 / Y0 2.03、BDF1 0.91 / 0.91**、sub-iter 倍増の差 ≤ 最小水準差の 1/1000、sub-iter 残差 4.7–6.0 桁低下; (iii) 凝縮ソース作動下 case/44 `run_0224`/`0225` (run_0170 収束場, dt 8e-6 / 4e-6 [物理 CFL≈10], 100/200 step): NaN 0、`condLimiterMode 1` 有効、モーメントは定常場に留まる (ドリフト g 8.5e-5〜1.2e-4 = 反復床、onset/g max 不変、floor 補正 ≤3e-12); (iv) restart case/16 `run_0487` (checkpoint) → `run_0488` (`history restored`, 連続 `run_0486` との差 ro 1.2e-6 = 反復差内), 旧形式 `run_0489` (BDF1 再開・完走), CPG + トレーサ case/09 `run_0079` (rms_roXi 毎 step 3.8–4.8 桁低下)。未実施: coupling 1/2 の dual-time 分岐の run、ピン行残差 0 の直接確認、過渡 IC の凝縮 dual-time。(run 番号は検証エージェントと重複したため Phase B 側を改番: case/16 0486–0490, case/44 0224/0225, case/09 0079–0086)
- `2026-09-17` — **Phase A 実装 (受動種の化学種経路化 + S3 安定化キー)**: 新規 `passiveTransport_d.cuh` / `passiveKernels_d.cuh` / `passiveLimiter_d.cuh` / `limiterFunctions_d.cuh`、化学種カーネル (Dirichlet/Neumann/ピン/勾配 [周期半割面除外引数]/DPLUR sweep) を受動種ポインタで共用、SLAU S3 分岐で受動種を ψ_P (セル局所スケールで無次元化した Venkat) 再構成・下限 0/[0,1] クリップ、`passive_bounds_d` (更新済み ρ で上下限 + 符号付き/絶対補正の体積積分を step 内/全期間で積算 `passiveFloorCorr_<prim>`)、`passive_diffusion_d` (トレーサ Fick)、`cond_moment_update_limited_passive_d` (relax / dtScale / DPLUR 増分入力; θ_u と診断は同じ)、周期 gather (`transport_diag_*`, 受動種・化学種勾配)・DPLUR dq mirror・coupling 2 クロス項の独立バッファ。キー: `passiveScalarScheme` (既定 0 = 旧経路ビット不変), `passiveImplicitRelax` (= implicitRelax), `speciesImplicitRelax` (1.0), `scalarCflMax`, `passiveImplicitCoupling` (自動: scheme 1 + SFR≥2 で 1)。単体 `tests/unit/test_passive_scalar.cu` ALL PASS (無次元化 Venkat は Q0 8.7e14 で有限 [旧式は NaN]; 1-D ステップ移流 S3 は L1 が 1 次の 0.155 倍・遷移幅 20→4 セル、S3 + Venkat + 前進 Euler は plateau 端で ξ>1 を 2 % 作り上限クランプで確定 [収支に記録])。**CFD (node, 短 run)**: 無影響 = HEAD 反復ノイズ内 (case/46 `run_0109`–`0112`, case/44 `run_0204`–`0209` [relax 0.7 含む], case/16 `run_0473`–`0475`); 制御試験 case/46 `run_0113` (scheme 1, SFR 0, relax 1.0): 残差・対角・生更新は同一、ξ と Y_EXH の差は再正規化 1.05e-4 + raw>ρ ノード (3053) のクランプで、クランプ無作用域では平均 4.6e-7・max 4.0e-5 (クランプ履歴の伝播) → **≤1e-6 は全域では不成立** (設計どおりの差; §6-1 ゲートを「再正規化・クランプ無作用域で平均 ≤1e-6」に改定); 周期 case/09 `run_0067`/`0068` (トレーサ 1 次/S3): 周期対 ΔXi = 0, ∫ρξ 3.7e-8 / 3.1e-7, S3 解析解誤差 8.9e-3 (1 次 4.0e-2); `run_0069`/`0070` 凝縮ソース seam 比 1.000000; 1 次経路 scheme 1 vs 0 (case/44 `run_0214`/`0215`): g 6.8e-5, onset 同一, floor 補正 0。**S3 smoke (case/44 `run_0170` 入力 2000 step)**: 化学種 `speciesImplicitCoupling 0` では S3 が不健全 (`run_0210`–`0213`: rms_ro 3–4e-2, condLim 0, モーメント floor 補正 1e-2 級; cfl 6 + pc1 は step 1098 で軸上 ρ<0 → NaN); **`speciesImplicitCoupling 1` では 4 本 (`run_0216`–`0219`: cfl 2/6 × pc 0/1) 完走**、流れ残差 1 次と同水準、pc1 が優位 (floor 補正 3e-6〜1e-5 vs pc0 1e-2, rms_rog 3e-6 vs 4e-5)。推奨 = S3 + `speciesImplicitCoupling 1` + `passiveImplicitCoupling 1` (自動既定)。step 時間 2.4→2.8 ms。
- `2026-09-17` — **原因確認 (実装前, node)**: (A) case/28 He/空気 coaxial を node 変換し `run_0065` (cfl 2) の場から S2/S3 × coupling 0/1 × cfl 2/4/6 (`run_0066`–`0078`, relax 0.7): S3 + coupling 0 は組成せん断層で化学種残差が先に成長し cfl 4 step 402 / cfl 6 step 165 で NaN、cfl 2 は緩やか (onset 334); S3 + coupling 1 は cfl 2/4/6 とも S3 固有の発散なし (化学種残差 ×0.7)。全 run は baseline の別問題 (node TP-SST の上境界 `outlet_statPress`/軸で rms_roe 主導、restart 後 ~1500–3000 step で発散; `run_0076` で同一設定の継続でも再現) で終わる。(B) case/46 `run_0104` の res_6000 を roXi:=roY0 で restart し 1/10/100 step を診断出力 (`run_0107`/`0108`, 診断は scratch build の env ゲート出力のみ): 残差・対角・生更新は float 精度で同一、差は化学種の再正規化 1.40e-4 (全体) + トレーサのクランプ 3.9e-5、発生点はカウル後縁の純排気ノード。
- `2026-09-17` — 初稿 (ユーザ決定 2026-09-16: トレーサを化学種カーネルの受動種に、凝縮モーメントも化学種経路、dual-time の化学種修正移植と受動種の BDF 項を一括で)。
- `2026-09-17` — 調査で前提を訂正 (§4.0): 本番の化学種移流も 1 次 (S3 は experimental)。本 plan の 2 次化 = S3 の node 本番化を含む。
- `2026-09-17` — codex plan レビュー 2 回目 **GO-with-changes (M6)** を全採用: リミッタの無次元化、上下限と補正収支、周期の coupling 1/2 整合、非一様組成の固定点ケース、`speciesImplicitRelax` の分離、BDF 履歴契約と時間次数の数値基準。
- `2026-09-17` — codex plan レビュー 1 回目 **NO-GO (C1/M7/m1)** を全採用して改訂: 残差ブレンドを撤回 (完全な $R_2$ + 増分緩和)、S3 発散の原因確認を実装前に、
  面クリップの保存性訂正と有界性の診断、周期 node の処理順、1.8e-4 は仮説に、dual-time の処理順と 3 水準時間精度、拡散の解析解試験、完了条件 1–7。

## 10. 未確定事項

- ~~S3 を既定にするか~~ 決着 (2026-09-17, §9): 既定は SFR 0 のまま。S3 は onset を 0.2 r_t / 0.7 mm 下流に動かす (前線の数値拡散減) が Wysłouzil では実験からさらに遠ざかる。onset の実験差 (~5 mm) の主因は移流次数ではなく核生成モデル側 (plan condensation-followups)。
- k/ω も化学種経路 (2 次) に乗せるか: ユーザは「全部化学種の経路」と述べたが、SST の生産項・壁関数との結合の検証が別途要るので本 plan では見送り、後続とする (要確認)。
- ~~ψ_P (受動種ごとの Venkat) と ψ_ρ のどちらが有界性と安定性で優れるか~~ 決着 (2026-09-17): 受動種は ψ_P (セル局所スケールで無次元化した Venkat) を採用。非一様組成の固定点ケース (case/16 `run_0476` 系) で floor 補正 0・交差 restart が反復ノイズ内、S3 の 1-D ステップで遷移幅 4 セル。ψ_ρ 版は実装していない (化学種との整合は不要)。
- ~~S3 の node 発散が coupling 1 (scalar-DPLUR) 固有だった場合の contingency~~ 決着 (2026-09-17, §9): 逆で coupling 0 固有。受動種も DPLUR sweep に乗せる (§4.2)。
- ~~モーメントの dual-time 時間次数 (§6-6) は制約由来で後続~~ 撤回 (2026-09-17, codex result-2 M1): 真因は #13 で導入した共通 θ 非負化が N_k=0 かつ d_k<0 の成分があると全モーメントの増分を θ=0 で止めていたこと (Q1=0 の 663 ノード; 蒸気枯渇ではない)。#20 で修正して次数を再評価する (必須)。
- 発達中の凝縮 dual-time で sub-iter 残差が 2 桁落ちない列 (Q0, 流れ, 化学種; case/44 `run_0229`/`0230`) は擬似時間反復の床 (cfl_pseudo 2 + θ_u クランプ) で、nSub では解消しない。cfl_pseudo 上げ・`passiveImplicitCoupling 1`・更新クランプの dual-time 挙動の切り分けは後続 (時間精度の判定は nSub ≥20 で行う)。
- node TP-SST の case/28 baseline 発散 (上境界 `outlet_statPress`/軸, rms_roe 主導, restart 後 ~1500–3000 step; `thermoHrefTemp` 無し・絶対基準 h が候補) は本 plan の外 → followups へ登録 (S3 の定量 A/B をこの case でやる前に要解決)。
