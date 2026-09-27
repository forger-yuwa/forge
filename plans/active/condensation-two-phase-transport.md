# 凝縮域の輸送: 物性は気相組成で、拡散は蒸気の勾配で

## メタ

- **area**: `condensation / thermophysics`
- **status**: `draft`
- **related_docs**:
  - [`methods/condensation.md`](../../methods/condensation.md) (二相 EOS・受動スカラー輸送・実現可能性クランプ)
  - [`methods/thermophysics.md`](../../methods/thermophysics.md) (§4 輸送係数, §5 化学種拡散とエネルギー結合)
- **related_plans**:
  - [`thermophysics-solver-owned-species-db.md`](thermophysics-solver-owned-species-db.md) (§10 から本 plan へ移管。lump の輸送物性展開と同じ `gasProperties_d` を触る)
  - [`species-passive-scalar-unification.md`](../accepted/species-passive-scalar-unification.md) (凝縮モーメントの受動種経路・`passive_diffusion_d`)
- **created**: `2026-09-27`
- **owner**: `Claude (主セッション) / ユーザ`

## 1. 目的

凝縮 ON の NS run で、液相を水蒸気として扱っている輸送の近似を正す。現状 (2026-09-27 コード読み取り):

- kinetic 輸送 (`viscMethod: 2`) の組成は `roY` (H2O は蒸気 + 液の総水分) から作り、液相分率 `g` を参照しない (`cuda_forge/gasProperties_d.cu:71-85`)。
  液の質量も水蒸気の分子として μ・λ・拡散係数の混合に入る。
- 化学種拡散は総水分 `Y_w` の勾配を気相の拡散係数 (分子 + 乱流 `μt/(ρSc_t)`) で Fick 拡散させ、蒸気のエンタルピー h_v を運ぶ。
  **液 `rog` (凝縮モーメント) は拡散しない** (受動拡散 `passiveDiffusion_d_wrapper` はトレーサからしか呼ばれない, `tracerTransport_d.cu:176`;
  2026-09-27 訂正: 当初「液は別経路で受動拡散」と書いたのは誤り)。したがって拡散で動くのは実質「蒸気」で (g は不変)、エネルギー勘定は蒸気として整合しているが、
  駆動勾配が `∇Y_w` (蒸気でなく総水分) である。差 (蒸気) の非負性は拡散では保証されず、毎ステップの実現可能性クランプ `0 ≤ rog ≤ roY_w`
  (`condensationRealizability_d.cuh:92-113`, 液を削る = 強制蒸発, 総水分・エネルギーは保存) に頼っている。

完了時には、μ・λ・拡散係数は気相組成で評価され、拡散は「蒸気の勾配で分子拡散・液は分子拡散なし (任意で微小な液 Sc)・乱流は同じ Sc_t」になり、
蒸気の非負性が拡散作用素の構造で保たれる (クランプは保険として残し、働いた量を監視する)。

## 2. スコープ

- **やる**:
  - (A) 輸送物性 (μ・λ・化学種拡散係数) を**気相組成** (蒸気 `Y_w − g`、`1 − g` で正規化) で評価する。液滴の懸濁効果 (粘性増加・有効熱伝導) は無視 (2026-09-27 ユーザ決定)。
  - (B) 拡散: 総水分の分子拡散流束を蒸気の勾配 `−ρD∇(Y_w − g)` に、液 `rog` の分子拡散は既定で 0 (安定化用に**微小な液シュミット数を任意指定**できる選択肢を残す; ユーザ要望)、
    乱流拡散は総水分・液とも同じ `Sc_t`。輸送変数 (`roY_w`, `rog`, モーメント) は分け直さない。
  - (C) 監視: 実現可能性クランプの液の補正量 (累積, 総液量比) を凝縮 run で常にログと判定ツールに出す。
- **やらない**: 液滴の懸濁効果のモデル化、二温度 (液滴温度) の輸送への反映、液滴の慣性・スリップ。

## 3. 関連 docs と前提 (観測事実)

- 物性の誤差の大きさ: `case/16.nozzle_wys/run_0483_passive_wys_s0_sfr0` (SST + 非平衡凝縮) の液最大点 (T 207.6 K, Y_w 0.0110, g 0.0109) で、
  液を蒸気として数えた μ は気相組成の μ より **−0.50 %** (N2/H2O の Chapman–Enskog + Wilke, 当方の Python 検算 2026-09-27)。
- 液滴体積分率 φ = gρ/ρ_l: `run_0483` の液最大点で ρ 0.201 kg/m³・g 0.0109・ρ_l 992 → **φ = 2.2e-6** (液密度最大点でも 2.7e-6)。懸濁の粘性補正 Einstein 2.5φ = **5.5e-6**、有効熱伝導 Maxwell ≤3φ = **6.6e-6** で、組成効果 (μ −0.5 %) の 1/1000 以下 → 「混合物の μ・λ = 気相組成の μ・λ」で確定 (2026-09-27)。液が輸送に効くのは (i) 気相組成 (§4.1) と (ii) 液流束が運ぶエンタルピー (§4.2 の h_l) の 2 つ。
- 実現可能性の実測 (最終スナップショット): `run_0483` で液/総水分 最大 0.991 (0.99 超 166 ノード)、蒸気分率の最小 9.7e-5、最終 step のクランプ量 1.5e-19 (実質 0)。
  `case/45.isobutane_m6_d155/run_0039_ns_final_cond` 最大 0.034、`case/42.isobutane_wt/run_0071_…_ns_cond_evap_default` 最大 0.26 (両者クランプ場の出力なし)。
  **計算途中のクランプ量は多くの run で記録が無い** (`[passive] clampBudget` のログは一部の経路でしか出ない, `speciesTransport_d.cu:1512-1524`)。
- ~~エネルギー項の現状 (旧記述: 乱流で運ばれる液が h_v を運び受け手の温度が高く出る)~~ → **2026-09-27 訂正**: 液は拡散しないので、拡散で動くのは蒸気で h_v は整合。欠けているのは**液と蒸気が逆向きに乱流混合することで運ばれる潜熱流束** `L·ρD_t∇g` (目標モデルでは質量は打ち消し、エンタルピーは L ぶん残る)。以下の見積もりはその量の旧解釈による値:見積もり (次元解析のみ): 誤差 `L·J_g` と乱流熱流束の比 ≈ (L/cp)(Δg/ΔT)(Pr_t/Sc_t) ≈ 2500 K × (0.01/50 K) ≈ 0.5、g が 1 % 動くごとに温度誤差 ~L·Δg/cp ≈ 25 K の規模。境界層に液が入る NS + 凝縮 run で効き得る。
- 移流の ΣY: 既定 (`speciesFaceReconstruction 0`) は 1 次風上で ṁ·Y_up の和が ṁ と整合。S3 (2 次) は種ごとの再構成・リミッタで面の ΣY_f が 1 とは限らず、面で正規化しているかは未確認。いずれも毎 step の `species_renormalize_d` (`speciesTransport_d.cu:152,859`) で ΣρY=ρ に比例配分するので、個々の種の総量は厳密には保存されない。
- **流束の実測 (2026-09-27, `case/16.nozzle_wys/analyze_liquid_diffusion_error.py`, 最終スナップショット)**: `run_0483` (旧経路) / `run_0482` (S3) の g>1e-4 のノードで、
  「液の勾配による余計な蒸気分子流束 / 本来の蒸気分子流束」= 1 (全点) = 入口組成一様で Y_w がほぼ一様なので**現行の水の拡散流束はほぼ 0**。
  欠けている潜熱の乱流流束 L·ρD_t|∇g| と熱流束 (λ + cp μt/Pr_t)|∇T| の比は、乱流域 (μt/μ>10) で p50 0.77/0.75、p95 4.8/5.3。
  D_w (H2O–N2) p50 8.2e-5 m²/s、D_t p50 1.6e-4 m²/s。**これは流束の大きさの比で、温度への影響ではない** (下)。
- **温度への影響は潜熱が相殺する (codex diagnose 2026-09-27, [`notes/reviews/2026-09-27-condensation-diffusion-error-diagnose.md`](../../notes/reviews/2026-09-27-condensation-diffusion-error-diagnose.md))**:
  二相 EOS `e = e_gas + g(R_wT − L)` (`condensationEOS_d.cuh:381`) では、液の乱流輸送 r_g に伴う潜熱流束 −L r_g と EOS の −L r_g が打ち消し、温度に効くのは `−R_wT r_g`
  (T≈210 K で L の約 4 %) だけ (局所定係数の見積もり)。したがって「流束比 ~5」「L·g/cp ≈ 26 K」は温度変化の予測にならない (当方の旧解釈を撤回)。
  **エネルギー式だけに潜熱流束を足し `rog` を拡散させない実装は、この相殺を壊すので禁止**。影響の大きさは §6 の A/B で測る。
- S3 の面組成は `[0,1]` クランプ後に ΣY=1 へ正規化 (`convectiveFlux_slau_d.inc.cuh:364`) され、化学種移流 `ṁ·Y_f` は float 丸めの範囲で連続の式と整合 (確認済み)。
- 経緯: 当初 (B) は「既知の近似として受け入れる」とした (2026-09-27 ユーザ) が、同日「将来的にはこの方針でいきたい、液 Sc は安定化用に残したい」「μ・λ はできれば考慮したい」に更新。

## 4. 設計方針 (未レビュー; 実装前に上位諮問と codex plan 段レビューを通す)

### 4.1 (A) 気相組成の輸送物性

- **2026-09-27 改訂 (輸送物性の作り直し後)**: μ・λ は種ごとの輸送物性 (`viscMethod: 2` + `physProp.transport`; plan `thermophysics-solver-owned-species-db` §4.3c、`feature/species-transport`) になり、セル・壁とも `transport_mix_Y` / 表引き `transport_mix_Y_tab` (`cuda_forge/transportMix_d.cuh`・`transportTables_d.cuh`) が輸送種の Y から組成を作る。旧 Wilke 経路は撤去済み。
- 凝縮 carrier (TP) のとき、この入口に渡す組成を気相組成にする: `Y_s^gas = Y_s/(1−g)` (s ≠ 凝縮種), `Y_w^gas = (Y_w − g)/(1−g)` (液を除き再正規化)。その後は通常どおり輸送種 X → lump 展開 → CEA 形混合。
  変更は `gas_transport_cell_Y` (`gasProperties_d.cu`) と壁の組成を作る箇所に限り、`transport_mix_Y` の式は変えない。
- 化学種拡散係数 (混合平均、`thermo_Dmix_species_f`) も同じ気相組成で評価する。CPG carrier (空気) と `viscMethod 0/1` は組成を見ないので変更なし。

### 4.2 (B) 拡散作用素の統一

- 総水分: 分子 `−ρD_w∇(Y_w − g)` + 乱流 `−(μt/Sc_t)∇Y_w`。液: 分子 `−(μ/Sc_l)∇g` (既定 `Sc_l = ∞` = 0, 任意指定) + 乱流 `−(μt/Sc_t)∇g`。
  → 蒸気 `Y_w − g` の流束は分子 `−ρD_w∇(Y_w−g)` + 乱流 `−(μt/Sc_t)∇(Y_w−g)` (+ 液 Sc を入れた分) となり、蒸気も拡散方程式に従う。
- **補正と乱流を分離する** (codex diagnose M): 分子流束だけを気相内の分率 `Y_k/(1−g)` で補正して `Σ_気相 J^mol* = 0`、乱流流束は液・凝縮モーメント `Q0/Q1/Q2` を含む全相に共通の `μt/Sc_t` で
  (Σ∇Y=0 なので補正不要)。面ごとに `Σ J_gas + J_l = 0`、`J_w = J_v + J_l` を満たす。エネルギーは `Σ h_s J_s(gas) + h_l J_l` (h_l = h_v − L(T), 既存 EOS と同じ基準)。
- **有限の液シュミット数 Sc_l**: 液の分子流束を総水分と全相の収支にも入れないと、蒸気 `J_v = J_w − J_l` に `+(μ/Sc_l)∇g` が残り、蒸気の非負性が構造的に保たれない (codex diagnose M)。
  opt-in として残す方針は維持し、入れる場合は全相の収支に整合させる。初回の A/B では液の分子拡散を 0 に固定する。
- 液シュミット数のキー名・既定・上限は実装前に決める (例 `condensation.condLiquidSchmidt`, 既定 0 = 分子拡散なし)。

### 4.3 (C) 監視

- 凝縮 run では実現可能性クランプの g 成分の累積 (符号付き・絶対・総液量比) を毎ログで出し、`CONVERGENCE_VERDICT` と並ぶ判定 (閾値は実装前に決める, 例: 総液量比 1e-6 超で WARN) にする。

## 5. 実装ステップ

1. (C) 監視を先に入れ、既存の NS + 凝縮 run を再実行して現状のクランプ量を測る (直す前の基準)。
2. (A) 気相組成の輸送物性。
3. (B) 拡散作用素の統一と液 Sc オプション。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4 の上位諮問と codex plan 段レビュー | `cuda_forge` の輸送・拡散を変えるので AGENTS.md エスカレーション 1・6。§4.2 の補正速度・エネルギー項の式と §6 の合否を確定してから | F |
| 1b | 影響の A/B (codex diagnose 2026-09-27) | 共通初期場 = `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1` の最終保存量 (`restart_field.py` で 2 run に複製)。同じバイナリ・S3・BC・CFL。**変える点は拡散モデルの 1 点**: A = 現行、B = 蒸気の分子+乱流拡散 + 液と Q0/Q1/Q2 の同じ Sc_t による乱流輸送 + 整合したエンタルピー流束 (試験用実装が必要)。4000 step (200 step ごと保存) から、両者が `check_convergence` PASS・報告量 STEADY になるまで延長。測る量: 共通初期場での **EOS に投影した温度変化率 δṪ = (∂T/∂U)(R_B − R_A)/V** (総水分・液の残差も EOS に通す)、定常での共通凝縮域の体積重み **p95\|ΔT\|** (事前判別閾値 1 K)、onset・出口 g・壁圧。試験用実装は `cuda_forge` の変更なので #1 の後 | F |
| 2 | (C) クランプ量の監視 | ログ常時出力 + 判定。合格: 既存 NS + 凝縮 run (case/16 `run_0483` 系, case/42 `run_0071` 系) の再実行で数値が出て、閾値判定が働く | O |
| 3 | (A) 気相組成の輸送物性 | §4.1。合格: g=0 のとき現行とビット一致、g>0 で気相組成を渡した独立参照 (`tests/unit/transport_reference.py`) と一致 (float 格納 ≤1e-5)。case/16 NS + 凝縮で μ の変化量を記録 | O |
| 4 | (B) 拡散作用素の統一 + 液 Sc | §4.2。合格は §6 (実装前に確定) | O |
| 5 | docs | `methods/condensation.md`・`methods/thermophysics.md`・`procedures/solver-settings.md` | O |

## 6. 検証 (骨子; 合否の数値は #1 で確定する)

- 単体: g=0 で現行と一致、g>0 で気相組成の μ・λ・D が独立参照計算と一致。
- 実現可能性: (B) 後、クランプなしで蒸気 ≥ 0 が保たれること (クランプ量が 0 または丸め程度)。
- 収支: 総水分・液・エネルギーの収支 (式と許容差は #1 で確定)。
- 回帰: case/16 Wysłouzil NS + 凝縮で onset・g・壁熱流束の変化を記録 (Euler run は不変であること)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `cuda_forge/gasProperties_d.cu`, `cuda_forge/wmlesWallModel_d.cu`, `cuda_forge/speciesTransport_d.cu`, `cuda_forge/condensationTransport_d.cu`, `cuda_forge/condensationRealizability_d.cuh`
- Euler (visc 0) の run は影響なし。NS + 凝縮 run の結果は変わる (変化量を記録)

## 8. 完了条件

- [ ] 関連 `methods/` を更新済み
- [ ] §6 を満たす
- [ ] codex レビュー 2 回を §6.1 に記録
- [ ] status を done にし、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-27` — §4.1 を輸送物性の作り直し後の実装 (`transport_mix_Y`、CEA 形混合、旧 Wilke 撤去) に合わせて改訂。作業は `feature/species-transport` で行う。
- `2026-09-27` — 流束の実測と codex diagnose (`notes/reviews/2026-09-27-condensation-diffusion-error-diagnose.md`) を反映: 液は拡散していない (旧記述を訂正)、潜熱は EOS で相殺するので流束比から温度影響を読まない、補正は分子流束だけ・乱流は全相共通 Sc_t、有限 Sc_l の整合、影響を測る A/B (#1b)。解析スクリプトの拡散係数の定数誤り (1e4 倍小) を修正 (比には影響なし)。
- `2026-09-27` — 起票。ユーザ決定: 懸濁効果は無視、μ・λ は気相組成で考慮したい、拡散は将来「蒸気の勾配で分子拡散・液は分子拡散なし (安定化用の微小な液 Sc は選択肢として残す)・乱流は同じ Sc_t」。
  種 DB plan (`thermophysics-solver-owned-species-db.md`) §10 から移管。
