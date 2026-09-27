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
- 化学種拡散は総水分 `Y_w` の勾配を気相の拡散係数で Fick 拡散させ (液も分子拡散する)、液 `rog` は別経路の受動拡散 `μ/(ρSc)+μt/(ρSc_t)` を持つ。
  総水分と液の拡散が別作用素なので、差 (蒸気) の非負性は拡散では保証されず、毎ステップの実現可能性クランプ `0 ≤ rog ≤ roY_w`
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
- 液滴体積分率 φ = gρ/ρ_l は 1e-6〜1e-5 程度で、懸濁の粘性補正 (Einstein 2.5φ) と有効熱伝導の補正は無視できる → 「混合物の μ・λ = 気相の μ・λ」が妥当な近似。
- 実現可能性の実測 (最終スナップショット): `run_0483` で液/総水分 最大 0.991 (0.99 超 166 ノード)、蒸気分率の最小 9.7e-5、最終 step のクランプ量 1.5e-19 (実質 0)。
  `case/45.isobutane_m6_d155/run_0039_ns_final_cond` 最大 0.034、`case/42.isobutane_wt/run_0071_…_ns_cond_evap_default` 最大 0.26 (両者クランプ場の出力なし)。
  **計算途中のクランプ量は多くの run で記録が無い** (`[passive] clampBudget` のログは一部の経路でしか出ない, `speciesTransport_d.cu:1512-1524`)。
- 経緯: 当初 (B) は「既知の近似として受け入れる」とした (2026-09-27 ユーザ) が、同日「将来的にはこの方針でいきたい、液 Sc は安定化用に残したい」「μ・λ はできれば考慮したい」に更新。

## 4. 設計方針 (未レビュー; 実装前に上位諮問と codex plan 段レビューを通す)

### 4.1 (A) 気相組成の輸送物性

- `gasProperties_d` (と `wmlesWallModel_d` の壁物性) で、凝縮 carrier (TP) のとき組成を `Y_s^gas = Y_s/(1−g)` (s ≠ 凝縮種), `Y_w^gas = (Y_w − g)/(1−g)` として X を作る。
- 化学種拡散係数 (混合平均) も同じ気相組成で評価する。CPG carrier (空気) と `viscMethod 0/1` は組成を見ないので変更なし。

### 4.2 (B) 拡散作用素の統一

- 総水分: 分子 `−ρD_w∇(Y_w − g)` + 乱流 `−(μt/Sc_t)∇Y_w`。液: 分子 `−(μ/Sc_l)∇g` (既定 `Sc_l = ∞` = 0, 任意指定) + 乱流 `−(μt/Sc_t)∇g`。
  → 蒸気 `Y_w − g` の流束は分子 `−ρD_w∇(Y_w−g)` + 乱流 `−(μt/Sc_t)∇(Y_w−g)` (+ 液 Sc を入れた分) となり、蒸気も拡散方程式に従う。
- 補正速度 (ΣJ = 0, `speciesTransport_d.cu:271`) とエネルギー項 `Σh_s J_s` は**気相の流束**で組む (蒸気のエンタルピーは h_v、液の乱流輸送は h_l)。実装前に式を確定する。
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
| 2 | (C) クランプ量の監視 | ログ常時出力 + 判定。合格: 既存 NS + 凝縮 run (case/16 `run_0483` 系, case/42 `run_0071` 系) の再実行で数値が出て、閾値判定が働く | O |
| 3 | (A) 気相組成の輸送物性 | §4.1。合格: 単体試験で g=0 のとき現行とビット一致、g>0 で気相組成の Wilke と一致。case/16 NS + 凝縮で μ の変化量を記録 | O |
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

- `2026-09-27` — 起票。ユーザ決定: 懸濁効果は無視、μ・λ は気相組成で考慮したい、拡散は将来「蒸気の勾配で分子拡散・液は分子拡散なし (安定化用の微小な液 Sc は選択肢として残す)・乱流は同じ Sc_t」。
  種 DB plan (`thermophysics-solver-owned-species-db.md`) §10 から移管。
