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

(§4.1–4.4 は調査結果を受けて具体化する: 触る関数・配列レイアウト・呼び出し順)

### 4.1 受動種の表現

- 化学種配列を `nSpecies` (熱力学に入る) + `nPassive` (受動) の連結にし、輸送カーネルは `nTransport = nSpecies + nPassive` 本を同じループで扱う。
  受動種の index は `nSpecies..nTransport-1`。名前は `Xi`、`g_<s>`, `Q0_<s>`, `Q1_<s>`, `Q2_<s>` (保存量 `roXi`, `rog_<s>`, ...) を**現行のまま**保ち、
  出力・restart・後処理・`species_meta.yaml`・`forge_species.py` の互換を壊さない。
- 熱力学は `nSpecies` 本だけを見る (面組成の R_mix/γ、EOS、`thermo_init_db`)。ΣY 正規化・ΣJ=0 補正・`speciesImplicitCoupling` の予測/commit も `nSpecies` 本に限定。
- 面再構成: `Yface` を `nTransport` 本に拡張し、受動種も MUSCL + Venkat (`limiter_<name>`) で面値を作る。凝縮モーメントは非負の床 (0) を面値でも守る (リミッタ後の負値は 0 にクリップし、
  現行 1 次のときと同じ実現可能性クランプで確定)。
- 拡散: トレーサは化学種と同じ Fick 形 ($D=\mu/(\rho Sc)+\mu_t/(\rho Sc_t)$, 粘性 run のみ)、ΣJ=0 補正とエンタルピー拡散からは除外。モーメントは拡散なし。
- 陰解法: 受動種は現行どおり segregated point-implicit (対角 = 1 次風上の流出流束 + 拡散対角 + `src_jac`)。化学種の結合予測/commit には入れない。
- 境界: 入口 Dirichlet (`Xi` / モーメント 0)、他は Neumann、node ピン、周期 gather/mirror は化学種の経路に受動種を追加する。

### 4.2 凝縮モーメントの移流

- `condensationTransport_d_wrapper` の `scalarTransportResidualMulti_d` 呼び出しを、化学種移流残差カーネル (受動種区間) に置き換える。`res_rog_s` 等と
  `transport_diag_*` の配列は現行のまま (ソース・更新クランプ・実現可能性が読む)。
- 面値の非負性: $g_f, Q_{k,f} \ge 0$ を保証 (リミッタ後クリップ)。$Q_1/Q_0 \le r_{30}$ 等の実現可能性は面値では要求しない (セル更新後のクランプに任せる)。

### 4.3 dual-time

- chem ブランチ e296f0d0 の `speciesShiftDualTimeLevels` / `speciesAddUnsteadyTimeTerm` (BDF1/BDF2, 対角 $V a/\Delta t$) / 予測→commit (δρ 基準は予測時点の ρ) / 再正規化を
  移植。受動種は同じ時間レベル配列 (`<name>P/PP`) と BDF 残差・対角を持ち、正規化は受けない。
- `condLimiterMode 1` は dual-time でも更新クランプを使う (自動降格を撤回)。`tracer × dualTime` の拒否を解除。
- 検証は物理 Δt 半減とサブ反復数変更で時間精度 (F-cf8 の要求)。

### 4.4 切替と互換

- `passiveScalarScheme: 1` (化学種経路) / `0` (旧汎用スカラ経路)。0 は回帰の対照用に残し、検証後に既定を 1 にする。

## 5. 実装ステップ

1. `methods/thermophysics.md` §5b/§5d に受動種を追記、`methods/condensation.md` 実装 §4 にモーメント移流の 2 次化、`methods/time_integration/` に dual-time の化学種/受動種 BDF。
2. forge: 受動種の登録 (`variables.cpp` / `speciesInit_d` の配列を `nTransport` に)、`convectiveFlux` の面組成拡張、化学種残差/拡散/境界/周期の受動種対応、
   `condensationTransport_d` と `tracerTransport_d` の移流を化学種経路へ。切替キー。
3. dual-time の移植 (chem e296f0d0) + 受動種の BDF 項 + 拒否/降格の解除。
4. 単体 (§6) → node 回帰 (case/46 `run_0104` 同一 run 比較, case/44 `run_0170` プロトコル, case/16 `run_0335` プロトコル, Arthur) → dual-time 時間精度 → docs 同期。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | 調査 (コード地図) と §4 の具体化 | 進行中 (2026-09-17) |
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
  4. **無影響**: `passiveScalarScheme: 0` で現行とビット一致 (case/44 `run_0170`, case/46 `run_0100`)、化学種のみの run (case/16 `run_0471`) は両設定で不変。
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

## 10. 未確定事項

- k/ω も化学種経路 (2 次) に乗せるか: ユーザは「全部化学種の経路」と述べたが、SST の安定性 (壁関数・生産項との結合) の検証が別途要るので本 plan では見送り、後続とする (要確認)。
- 受動種の面値に負値クリップを入れると保存性が局所的に崩れる (面での質量が一致しない) — モーメントは元々クランプで非保存なので許容するが、量を診断で記録する。
