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
| 6 | S3 本番化 (ステップ 4) | 増分緩和、`scalarCflMax`、固定点不変 |
| 7 | dual-time 移植 + 受動種 BDF (ステップ 5) | 処理順・履歴・restart |
| 8 | 検証 (§6 1–8, node のみ) と codex result レビュー | 完了条件 §8 |

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

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
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

- `2026-09-17` — **Phase A 実装 (受動種の化学種経路化 + S3 安定化キー)**: 新規 `passiveTransport_d.cuh` / `passiveKernels_d.cuh` / `passiveLimiter_d.cuh` / `limiterFunctions_d.cuh`、化学種カーネル (Dirichlet/Neumann/ピン/勾配 [周期半割面除外引数]/DPLUR sweep) を受動種ポインタで共用、SLAU S3 分岐で受動種を ψ_P (セル局所スケールで無次元化した Venkat) 再構成・下限 0/[0,1] クリップ、`passive_bounds_d` (更新済み ρ で上下限 + 符号付き/絶対補正の体積積分を step 内/全期間で積算 `passiveFloorCorr_<prim>`)、`passive_diffusion_d` (トレーサ Fick)、`cond_moment_update_limited_passive_d` (relax / dtScale / DPLUR 増分入力; θ_u と診断は同じ)、周期 gather (`transport_diag_*`, 受動種・化学種勾配)・DPLUR dq mirror・coupling 2 クロス項の独立バッファ。キー: `passiveScalarScheme` (既定 0 = 旧経路ビット不変), `passiveImplicitRelax` (= implicitRelax), `speciesImplicitRelax` (1.0), `scalarCflMax`, `passiveImplicitCoupling` (自動: scheme 1 + SFR≥2 で 1)。単体 `tests/unit/test_passive_scalar.cu` ALL PASS (無次元化 Venkat は Q0 8.7e14 で有限 [旧式は NaN]; 1-D ステップ移流 S3 は L1 が 1 次の 0.155 倍・遷移幅 20→4 セル、S3 + Venkat + 前進 Euler は plateau 端で ξ>1 を 2 % 作り上限クランプで確定 [収支に記録])。**CFD (node, 短 run)**: 無影響 = HEAD 反復ノイズ内 (case/46 `run_0109`–`0112`, case/44 `run_0204`–`0209` [relax 0.7 含む], case/16 `run_0473`–`0475`); 制御試験 case/46 `run_0113` (scheme 1, SFR 0, relax 1.0): 残差・対角・生更新は同一、ξ と Y_EXH の差は再正規化 1.05e-4 + raw>ρ ノード (3053) のクランプで、クランプ無作用域では平均 4.6e-7・max 4.0e-5 (クランプ履歴の伝播) → **≤1e-6 は全域では不成立** (設計どおりの差; §6-1 ゲートを「再正規化・クランプ無作用域で平均 ≤1e-6」に改定); 周期 case/09 `run_0067`/`0068` (トレーサ 1 次/S3): 周期対 ΔXi = 0, ∫ρξ 3.7e-8 / 3.1e-7, S3 解析解誤差 8.9e-3 (1 次 4.0e-2); `run_0069`/`0070` 凝縮ソース seam 比 1.000000; 1 次経路 scheme 1 vs 0 (case/44 `run_0214`/`0215`): g 6.8e-5, onset 同一, floor 補正 0。**S3 smoke (case/44 `run_0170` 入力 2000 step)**: 化学種 `speciesImplicitCoupling 0` では S3 が不健全 (`run_0210`–`0213`: rms_ro 3–4e-2, condLim 0, モーメント floor 補正 1e-2 級; cfl 6 + pc1 は step 1098 で軸上 ρ<0 → NaN); **`speciesImplicitCoupling 1` では 4 本 (`run_0216`–`0219`: cfl 2/6 × pc 0/1) 完走**、流れ残差 1 次と同水準、pc1 が優位 (floor 補正 3e-6〜1e-5 vs pc0 1e-2, rms_rog 3e-6 vs 4e-5)。推奨 = S3 + `speciesImplicitCoupling 1` + `passiveImplicitCoupling 1` (自動既定)。step 時間 2.4→2.8 ms。
- `2026-09-17` — **原因確認 (実装前, node)**: (A) case/28 He/空気 coaxial を node 変換し `run_0065` (cfl 2) の場から S2/S3 × coupling 0/1 × cfl 2/4/6 (`run_0066`–`0078`, relax 0.7): S3 + coupling 0 は組成せん断層で化学種残差が先に成長し cfl 4 step 402 / cfl 6 step 165 で NaN、cfl 2 は緩やか (onset 334); S3 + coupling 1 は cfl 2/4/6 とも S3 固有の発散なし (化学種残差 ×0.7)。全 run は baseline の別問題 (node TP-SST の上境界 `outlet_statPress`/軸で rms_roe 主導、restart 後 ~1500–3000 step で発散; `run_0076` で同一設定の継続でも再現) で終わる。(B) case/46 `run_0104` の res_6000 を roXi:=roY0 で restart し 1/10/100 step を診断出力 (`run_0107`/`0108`, 診断は scratch build の env ゲート出力のみ): 残差・対角・生更新は float 精度で同一、差は化学種の再正規化 1.40e-4 (全体) + トレーサのクランプ 3.9e-5、発生点はカウル後縁の純排気ノード。
- `2026-09-17` — 初稿 (ユーザ決定 2026-09-16: トレーサを化学種カーネルの受動種に、凝縮モーメントも化学種経路、dual-time の化学種修正移植と受動種の BDF 項を一括で)。
- `2026-09-17` — 調査で前提を訂正 (§4.0): 本番の化学種移流も 1 次 (S3 は experimental)。本 plan の 2 次化 = S3 の node 本番化を含む。
- `2026-09-17` — codex plan レビュー 2 回目 **GO-with-changes (M6)** を全採用: リミッタの無次元化、上下限と補正収支、周期の coupling 1/2 整合、非一様組成の固定点ケース、`speciesImplicitRelax` の分離、BDF 履歴契約と時間次数の数値基準。
- `2026-09-17` — codex plan レビュー 1 回目 **NO-GO (C1/M7/m1)** を全採用して改訂: 残差ブレンドを撤回 (完全な $R_2$ + 増分緩和)、S3 発散の原因確認を実装前に、
  面クリップの保存性訂正と有界性の診断、周期 node の処理順、1.8e-4 は仮説に、dual-time の処理順と 3 水準時間精度、拡散の解析解試験、完了条件 1–7。

## 10. 未確定事項

- k/ω も化学種経路 (2 次) に乗せるか: ユーザは「全部化学種の経路」と述べたが、SST の生産項・壁関数との結合の検証が別途要るので本 plan では見送り、後続とする (要確認)。
- ψ_P (受動種ごとの Venkat) と ψ_ρ のどちらが有界性と安定性で優れるかは §6-2/6-4 で実測して決める (両方を切替可能に実装)。
- ~~S3 の node 発散が coupling 1 (scalar-DPLUR) 固有だった場合の contingency~~ 決着 (2026-09-17, §9): 逆で coupling 0 固有。受動種も DPLUR sweep に乗せる (§4.2)。
- node TP-SST の case/28 baseline 発散 (上境界 `outlet_statPress`/軸, rms_roe 主導, restart 後 ~1500–3000 step; `thermoHrefTemp` 無し・絶対基準 h が候補) は本 plan の外 → followups へ登録 (S3 の定量 A/B をこの case でやる前に要解決)。
