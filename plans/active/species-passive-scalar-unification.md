# 受動スカラ (排気トレーサ・凝縮モーメント) を化学種カーネルに乗せ、dual-time の物理時間項を揃える

## メタ

- **area**: `convection / diffusion / condensation / time_integration`
- **status**: `draft`
- **related_docs**:
  - `methods/thermophysics.md` (実装 §5b 多成分化学種輸送, §5d 拡散, §5 tracer)
  - `methods/condensation.md` (実装 §4 モーメント輸送・§4c 更新クランプ)
  - `methods/time_integration/` (dual-time, scalar-DPLUR)
- **related_plans**: [thermophysics-cea-mole-fraction-species.md](../accepted/thermophysics-cea-mole-fraction-species.md) (トレーサ `roXi` の起源, F-sp1),
  [condensation-source-limiter-steady.md](../accepted/condensation-source-limiter-steady.md) (更新クランプ `cond_moment_update_limited_d`),
  [condensation-followups.md](condensation-followups.md) F-cf8 (モーメントの dual-time 物理時間項) / F-cf9 (旧 mode 0 削除),
  chem ブランチの dual-time 化学種修正 (`feature/chemistry-finite-rate` e296f0d0; [[dualtime-species-frozen-bug]])
- **created**: `2026-09-17`
- **owner**: Claude (session 011spCYH)

## 1. 目的

排気トレーサ `roXi` と凝縮モーメント `rog_s, roQ0_s..roQ2_s` は汎用スカラコア (1 次風上、拡散なし/μ 係数形) で解かれていて、
化学種 `roY_s` (2 次 MUSCL + Venkatakrishnan、Fick + 乱流拡散) と数値の扱いが違う。同じ排気率でも `lumped` の $Y_{EXH}$ と `full` の
$\xi$ が 1.8e-4 ずれ (case/46 `run_0104`)、凝縮 onset は 1 次の数値拡散を受ける。ユーザ決定 (2026-09-16): **受動スカラは全て化学種の経路を通す**。
完了時: (a) トレーサは化学種カーネルの「EOS と ΣY から除外した受動種」として 2 次 + 拡散 (層流 Sc + 乱流 Sc_t) で輸送され、同一 run 内で
$\xi$ と $Y_{EXH}$ が float 精度で一致する; (b) 凝縮モーメントも同じ 2 次面再構成で移流される (ソース・更新クランプ・実現可能性は現行のまま);
(c) dual-time で化学種が更新されない main のバグを chem ブランチから移植し、受動種 (トレーサ・モーメント) にも BDF 物理時間項が入る
(F-cf8 を閉じ、`tracer × dualTime` の拒否を解除)。

## 2. スコープ

- **やる**:
  - 化学種輸送に「受動種」の区分を追加: 輸送 (面再構成・移流残差・拡散・陰解法対角・境界・周期・restart・出力・残差列) は化学種と同じ経路、
    **熱力学 (MW/cp/h/R, EOS, 対流流束の面組成 R_mix/γ)・ΣY 正規化・ΣJ=0 補正・エンタルピー拡散・`speciesImplicitCoupling` の予測/commit・入口の X/Y 検証・
    `condGasSpecies`** からは除外。
  - トレーサ `roXi` を受動種として登録 (拡散あり: $D = \mu/(\rho\,Sc) + \mu_t/(\rho\,Sc_t)$、config `Sc`/`Sc_t` を共用)。汎用スカラ経路の `tracerTransport_d` は撤去
    (残す場合は A/B 用の `passiveScalarScheme: 0` に限定)。
  - 凝縮モーメント 4 本 (× nCondSpecies) を受動種として移流 (拡散なし、面再構成は化学種と同じ MUSCL + リミッタ)。ソース (`condensationSource_d`)・更新クランプ
    (`cond_moment_update_limited_d`)・実現可能性クランプ・EOS 結合 (g) は現行のまま `res_*`/`transport_diag_*` を受け取る。
  - dual-time: chem ブランチの化学種修正 (時間レベル `roY{s}P/PP`、BDF 残差 + 対角、予測→commit、再正規化) を main に移植し、受動種にも同じ BDF 項を付ける。
    `condLimiterMode 1` の dual-time 自動降格と `tracer × dualTime` 拒否を解除 (検証後)。
  - 切替キー `passiveScalarScheme` (1 = 化学種経路 [既定候補], 0 = 旧汎用スカラ経路) を A/B・回帰のために残す (最終的に 0 は F-cf9 と同時に削除)。
- **やらない**: k/ω の 2 次化 (SST の安定性検証が別途要る; 後続 plan)、リミッタの新設 (化学種の Venkat をそのまま使う)、化学反応、cell 離散化の検証 (ユーザ指示で cell は使わない)。

## 3. 関連 docs と前提

- 化学種の面組成は `convectiveFlux` が MUSCL + `limiter_Y{s}` で決め、化学種移流残差は同じ面組成 (`Yface`) で $\sum \dot m\,Y_f$ を組む
  (`methods/thermophysics.md` §5b)。受動種は同じ `Yface` 配列を拡張して使う。
- 汎用スカラコア (`scalarTransport_d.cu`) は 1 次風上 + 任意の μ 係数拡散 + 擬似時間 point-implicit。k/ω はこのまま残す。
- 凝縮モーメントの定常固定点は残差だけで決まる (更新クランプは収束時に無作用) ので、移流離散化の変更は固定点を変える (2 次化で onset の数値拡散が減る)。
  したがって既存の凝縮回帰 (case/44 `run_0170`, case/16 `run_0335`, Arthur) は**再取得して差を記録する**対象で、ビット一致は要求しない。
- dual-time の化学種バグ ([[dualtime-species-frozen-bug]]): main の `advanceImplicitDualTime` は化学種を更新しない。修正は chem ブランチ e296f0d0 のみ。

## 4. 設計方針

### 4.0 前提の訂正 (調査 2026-09-17)

- **本番の化学種移流も 1 次風上**である。2 次の面組成再構成 (`time.deltaT.speciesFaceReconstruction`; S2 = face thermo のみ、S3 = 移流も) は
  [convection-multispecies-contact-pressure](../accepted/convection-multispecies-contact-pressure.md) で実装済みだが既定 0 で、S3 (=2) は
  「2 次流束 (RHS) と 1 次風上 `transport_diag` (LHS) の defect-correction 不整合で cfl 4 で発散 (case/28)」のため cfl≤2 限定の experimental。
  case/44 `run_0170` も case/46 `run_0100`–`0106` も `speciesFaceReconstruction` 未指定 = 化学種は `scalarTransportResidualMulti_d` (1 次)。
- したがって case/46 `run_0104` の「同一 run 内 |ξ − Y_EXH| 1.8e-4 はカーネルの次数差」という記録は誤りで、両者とも 1 次風上。残る差は
  **化学種だけに掛かる ΣρY=ρ 再正規化 (`species_renormalize_d`) と、トレーサだけの primitive 段クランプ (0≤roXi≤ρ を保存量に書き戻す)**、および fused/single
  カーネルの加算順の違いである (本 plan の検証 1 で切り分ける)。README/accepted plan の該当記述は本 plan で訂正する。
- 本 plan の「2 次化」= **S3 を化学種・受動種の両方で本番化する**こと。鍵は LHS 不整合の解消で、これが無いと cfl 6 + `implicitRelax 0.7` (現行推奨) で回らない。

### 4.1 受動種の表現 (化学種カーネルの再利用)

- 化学種の device ポインタ配列 (`g_roY_dev` 等, `speciesInit_d`) と同じ形の**受動種配列** (`g_pas_*`: roP, roPN, res, transport_diag, src_jac, P, dP/dx.., Pface) を
  `nPassive = (tracer ? 1 : 0) + 4·nCondSpecies` 本で持つ。名前は現行のまま (`roXi`/`Xi`, `rog_s`/`g_s`, `roQ2_s`/`Q2_s`, `roQ1_s`, `roQ0_s`; 順序は
  `variables::registerCondensation` の {g, Q2, Q1, Q0}) → 出力・restart・後処理・`species_meta.yaml`・周期 gather (`res_*` 名) の互換を壊さない。
- **熱力学は `cfg.nSpecies` 本のまま** (EOS・面 R_mix/γ・粘性混合・入口 X/Y 検証・`condGasSpecies`・ΣY 正規化・ΣJ=0 補正・エンタルピー拡散・
  `speciesImplicitCoupling` 予測/commit は受動種に触れない)。受動種は `nSpecies==1` (CPG + トレーサ) でも動く (`speciesEnabled` ゲートとは独立)。
- **移流**: `SpeciesArgs` に受動種の再構成入力 (`nPassive, Pd_recon, dPdx/y/z_recon, Pface_out`) を足し、SLAU の S3 分岐 (`convectiveFlux_slau_d.inc.cuh` 267–296)
  と同じ `interp_dispatch` (リミッタは化学種と同じ **ψ_ρ = `limiter_ro`**; 新規リミッタ無し) で面値を作り、**正規化はせず**下限 0 (トレーサは [0,1]) でクリップして
  upwind 側を `Pface_out[ip*nPassive+q]` に書く。残差は `species_advection_faceY_d` を受動種ポインタで呼ぶ (`res_<cons>`, 1 次風上の `transport_diag_<prim>`)。
  SLAU 以外の対流スキーム、または `speciesFaceReconstruction < 2` のときは、受動種も化学種と同じ 1 次経路 (`scalarTransportResidualMulti_d`) に落ちる
  (常に化学種と**同じ次数**)。勾配は `species_gradient_d` (Green–Gauss, 軸対称対応) を受動種ポインタで流用。
- **拡散**: トレーサのみ、化学種の Fick 形 $D = \mu/(\rho Sc) + \mu_t/(\rho Sc_t)$ (`speciesDiffusionMethod 0` 相当の定数 Sc、粘性 run のみ) を ΣJ=0 補正と
  エンタルピー項なしで加える小カーネル `passive_diffusion_d` (species_diffusion_d の D 計算を流用)。モーメントは拡散なし (液滴は拡散しない)。
- **境界・ピン・周期**: 化学種の Dirichlet/Neumann/ピン (`species_dirichlet_boundary_d` 系) を受動種ポインタで呼ぶ (トレーサ入口 `Xi`、モーメント入口 0)。
  周期 gather/mirror は現行の名前ベース登録のまま。
- **更新**: 受動種は現行どおり segregated point-implicit (`scalarTimeIntegration_d` / `cond_moment_update_limited_d`)。化学種の結合予測/commit には入れない。
  トレーサの primitive 段クランプは **保存量を書き換えない** (Xi の表示だけクリップ; 化学種と同じ扱い)。
- 切替 `passiveScalarScheme` (1 = 化学種経路 [検証後の既定], 0 = 旧汎用スカラ経路 [A/B 用])。

### 4.2 S3 の本番化 — LHS/RHS 不整合の解消 (化学種・受動種共通)

defect-correction (RHS 2 次・LHS 1 次) の segregated point-implicit は、増分が $\delta = \Delta\tau R_2/(V + \Delta\tau(\mathrm{diag}_1 + \mathrm{sj}))$ で、
$\Delta\tau\to\infty$ で $\delta \to R_2/\mathrm{diag}_1$ と 2 次流束の反対称部分が減衰せず発散する (case/28 cfl 4 の指紋)。対策を 2 段で入れる:

1. **増分の緩和**: スカラ更新に `implicitRelax` (流れと同じ 0.7 既定) を掛ける (現行スカラ更新は緩和なし)。
2. **2 次流束の deferred-correction 分割**: $R_2 = R_1 + (R_2 - R_1)$ とし、$(R_2-R_1)$ に緩和係数 $\omega_{dc}$ (config `speciesDeferredCorrection`, 既定 0.5) を掛ける。
   固定点は $\omega_{dc}$ に依らず $R_2 = 0$ (定常で $R_2 - R_1$ の項は残差として全量残る) なので**解は変わらない**が、擬似時間の安定領域が広がる
   (標準的な segregated ソルバの手法)。$R_1$ は同じ面ループで `phi_upwind` から追加コストなしに得る。
3. それでも cfl 6 で発散する場合の保険: 化学種/受動種だけ擬似 Δτ を絞る `scalarCflMax` (既定なし)。

検証は case/28 (S3 が cfl 4 で発散した実例; node で再現) と case/44 cfl 6 + relax 0.7 で行う。

### 4.3 凝縮モーメントの移流

- `condensationTransport_d_wrapper` の `scalarTransportResidualMulti_d` を受動種の移流 (4.1) に置換。`res_<cons>` と `transport_diag_<prim>` の名前・ゼロ初期化・
  順序 (移流 → ソース → 更新) は現行のまま。ソース (`condensationSource_d`)、更新クランプ (`cond_moment_update_limited_d`)、実現可能性クランプ、EOS 結合 (g) は不変。
- 面値の非負クリップは局所的に非保存 (面で質量が一致しない) だが、モーメントは元々クランプで非保存なので許容し、量を `condClampCorr` と同様に診断に載せる。

### 4.4 dual-time (chem e296f0d0 の移植 + 受動種)

- 化学種: `species_shift_levels_d` (`roY{s}P/PP`)、`species_add_unsteady_d` (BDF 残差 + 対角 $Va/\Delta t$; 最初の 2 step は BDF1)、`g_roPred_dev` による
  予測→commit (plan-C は予測時点の ρ を δρ 基準に)、再正規化、primitive を `advanceImplicitDualTime` に入れる。chem 固有 (CMC の再正規化、`CMC7_TFREEZE_BYPASS`、
  `jacobianMode 2` の block LU、`dt_local_sp`) は持ち込まない。本ブランチ固有 (float32 拡散、node ピン、fused 残差、S3 の `Yface`) は保つ。
- 受動種: 同じ時間レベル配列 (`<cons>P/PP`) と BDF 残差・対角を `cond_moment_update_limited_d` / `scalarTimeIntegration_d` の入力 (`res_*`, `transport_diag_*`) に
  足す (更新式は不変)。`condLimiterMode 1` の dual-time 自動降格と `tracer × dualTime` の拒否を解除。k/ω は既に BDF 項あり (`update_d.cu`)。

### 4.5 k/ω

対象外 (§10)。同じ受動種経路に乗せられる構造にはなるが、SST の生産項・壁関数との結合の検証は別 plan。

## 5. 実装ステップ

1. `methods/thermophysics.md` §5b/§5d に受動種を追記、`methods/condensation.md` 実装 §4 にモーメント移流の 2 次化、`methods/time_integration/` に dual-time の化学種/受動種 BDF。
2. forge: 受動種の device 配列 (`passiveInit_d`)、`SpeciesArgs` の受動種拡張と SLAU S3 分岐の面値、`species_advection_faceY_d`/`species_gradient_d`/境界/ピンの
   受動種呼び出し、`passive_diffusion_d`、`condensationTransport_d`・`tracerTransport_d` の移流置換、`passiveScalarScheme`。
2b. forge: S3 の本番化 — スカラ更新の `implicitRelax`、deferred-correction 係数 `speciesDeferredCorrection`、保険の `scalarCflMax`。
3. dual-time の移植 (chem e296f0d0) + 受動種の BDF 項 + 拒否/降格の解除。
4. 単体 (§6) → node 回帰 (case/46 `run_0104` 同一 run 比較, case/44 `run_0170` プロトコル, case/16 `run_0335` プロトコル, Arthur) → dual-time 時間精度 → docs 同期。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~調査 (コード地図) と §4 の具体化~~ | 済 (2026-09-17): §4.0–4.5 |
| 2 | codex plan レビュー | §4/§6 確定後、実装前 |
| 3 | docs 先行更新 | ステップ 1 |
| 4 | 受動種の実装 (トレーサ・モーメント) | ステップ 2 |
| 5 | dual-time 移植 + 受動種 BDF | ステップ 3 |
| 6 | 回帰 (node のみ) と codex result レビュー | ステップ 4 |

## 6. 検証

- **単体 / ビルド**: 化学種の既存単体 (`test_species_eos_cross` 等) が不変。受動種の 1 セル/1 次元移流試験: 一様流中のステップ状 ξ が 2 次で輸送され (1 次より鋭い)、
  0 ≤ ξ ≤ 1 を保つ; モーメントの面値非負; dual-time の BDF 項が物理 Δt 半減で誤差 1/4 (BDF2) になる 1 次元移流。
- **検証ケース (node のみ)**:
  1. **トレーサ = 化学種の一致**: case/46 `run_0104` プロトコル (lumped [EXH, AIR] + tracer, m6_on node Euler 6000 step) で同一 run 内 $|\xi - Y_{EXH}|$ が
     float 精度 (max ≤ 1e-5、平均 ≤ 1e-7; 現行 1.8e-4 / 4.6e-7)。`full` (`run_0101` プロトコル) の $\xi$ と lumped の $Y_{EXH}$ も同程度。力 C_T/C_L/C_M は現行と 1e-3 内、
     時系列 0.1 % STEADY。
  2. **モーメント 2 次化の凝縮回帰**: case/44 `run_0170` プロトコル (入口 Tt 分布, cfl 2 → 推奨 cfl 6 + relax 0.7 も) で `passiveScalarScheme` 1 vs 0 の onset・出口 g・g max・
     series STEADY・condLim 1・補正 0 を記録 (差は 2 次化の効果として採否判断; ビット一致は要求しない); case/16 `run_0335` プロトコル (Wysłouzil) で onset と実験の差
     (現行 ~5 mm 下流) がどう動くか; Arthur N2 node (`case/34 run_0106` プロトコル) の onset。同一バイナリ反復ノイズ床を併記。
  3. **dual-time**: (i) 多成分 dual-time で ΣY=1 と化学種の時間発展 (main のバグ修正; case/16 TP dual-time の短い run で `roY` が動くこと)、(ii) モーメント/トレーサの
     BDF: 物理 Δt を半減してサブ反復収束を揃え、$\xi$ とモーメントの時間発展が Δt に依らない (差 ≤ BDF2 の 1/4 スケーリング)。
  4. **無影響**: `passiveScalarScheme: 0` + `speciesFaceReconstruction 0` で現行とビット一致 (case/44 `run_0170`, case/46 `run_0100`)、化学種のみの run
     (case/16 `run_0471`) は `passiveScalarScheme` に依らず不変。
  5. **S3 の本番化 (安定性)**: case/28 (S3 が cfl 4 で発散した実例) を node で `speciesFaceReconstruction 2` + `implicitRelax 0.7` + deferred-correction で cfl 4/6 完走・
     残差床が S2 と同等; case/44 `run_0170` プロトコルを `speciesFaceReconstruction 2` + cfl 6 + relax 0.7 で完走 (NaN 0, series STEADY)。固定点が
     `speciesDeferredCorrection` に依らないこと (0.5 vs 1.0 で場が反復ノイズ内)。
  6. **前提訂正の裏付け**: case/46 `run_0104` プロトコルで、`passiveScalarScheme 1` (受動種も化学種と同じ 1 次経路、再正規化なし・クランプなし) にすると同一 run 内
     |ξ − Y_EXH| が float 精度になることを示し、残差 1.8e-4 の由来 (再正規化/クランプ) を README・accepted plan で訂正する。
- **判定基準**: 上のゲート + `check_convergence.py` / `check_quasisteady.py` VERDICT、NaN 0、step 時間の増分を記録 (受動種 5 本の 2 次化でどれだけ増えるか)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/{speciesTransport_d.cu,cuh, convectiveFlux*_d.cu, condensationTransport_d.cu, tracerTransport_d.*, periodicNode_d.cu, scalarTransport_d.cu}`,
  `variables.cpp`, `main.cpp` (定常・RK・dual-time), `input/solverConfig.*`。
- 既存 run: `passiveScalarScheme` 省略時の既定を 1 にすると凝縮 run の固定点が変わる (2 次化)。回帰 run の再取得と README 記録が必要。
- docs: `methods/thermophysics.md`, `methods/condensation.md`, `methods/time_integration/`, `procedures/solver-settings.md`, `procedures/recommended-settings.md`。

## 8. 完了条件

- [ ] 関連 methods を更新済み
- [ ] 実装・§6 の検証 1–4 を満たす
- [ ] codex レビュー 2 回を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status: done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`、`plans/README.md` 同期、followups F-cf8 / F-sp1 を閉じる

## 9. 変更ログ

- `2026-09-17` — 初稿 (ユーザ決定 2026-09-16: トレーサを化学種カーネルの受動種に、凝縮モーメントも化学種経路、dual-time の化学種修正移植と受動種の BDF 項を一括で)。
- `2026-09-17` — 調査で前提を訂正 (§4.0): 本番の化学種移流も 1 次 (S3 は experimental)。本 plan の 2 次化 = S3 の本番化 (LHS 不整合の解消, §4.2) を含む。

## 10. 未確定事項

- k/ω も化学種経路 (2 次) に乗せるか: ユーザは「全部化学種の経路」と述べたが、SST の安定性 (壁関数・生産項との結合) の検証が別途要るので本 plan では見送り、後続とする (要確認)。
- 受動種の面値に負値クリップを入れると保存性が局所的に崩れる (面での質量が一致しない) — モーメントは元々クランプで非保存なので許容するが、量を診断で記録する。
