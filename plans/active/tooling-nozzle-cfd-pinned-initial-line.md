# CFD ピン: MOC の初期線 (スロート特性線) を CFD の遷音速場から取る

## メタ

- **area**: `tooling / nozzle design`
- **status**: `draft`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (§軸 Mach law・§壁の決め方 — 初期線の出所に CFD ピンを追記する)
- **related_plans**:
  - [`verification-m6-axis-wave-mesh-su2.md`](verification-m6-axis-wave-mesh-su2.md) (§5.1 #13・#15・#16、§9 2026-10-05: 本 plan の動機。V0/V4/V4b の Euler A/B と寄与率)
  - [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (§4.9 C2 recipe — 本 plan が通れば初期線の出所を recipe の段に入れる)
  - 旧 CFD アンカー (cell 時代、`design/forge_design/feedback/cfd_anchor.py` `CFDThroat`、軸の縦線初期値): 今回は「スロート特性線上の状態」を取る別物
- **created**: `2026-10-05`
- **owner**: `Claude (主セッション) / ユーザ`

## 1. 目的

M6 ノズル (case/45) の設計チェーンでは、MOC の初期線を Hall の遷音速解 (R=2 の円弧スロートの 3 次系列) で与えている。
Hall の初期線は MOC の C⁻ 適合条件を積算 0.12〜0.13° 満たさず (codex 4 回目の再計算)、スロート直後の MOC 壁に約 0.1° の角度オフセットを作る。
さらに Euler A/B で、**スロート直上流 1.5 r_t の壁の 1.5e-3 r_t の違いが試験部の圧力の傾きを 0.30 %pt 動かす**ことが分かった (verification-m6 §9, V4b)。
試験部の一様性がスロート近傍の扱いに強く依存するので、初期線を実際の流れ (CFD) に合わせ、Hall の打切り誤差を設計から外す。
完了時には、problem YAML の 1 キーで「Hall」か「CFD の場」かを初期線の出所として選べ、CFD ピンの設計が Hall 版より試験部の一様性で劣らないことを Euler で確かめた状態にする。

ユーザ決定 (2026-10-05): 「CFD ピンは形状を詰める段階の選択肢」→ V4b の結果を受けて前倒し (「やって」)。

## 2. スコープ

- **やる**:
  - CFD 場 (node Euler) からスロート特性線を抽出する provider (`HallThroat` 互換の `throat_characteristic(n)`・`mach(x, 0)`・`axis_anchor(x0)`)。
  - `design_chain` に初期線の出所の切り替え (`geometry.initial_line: hall | cfd`, `geometry.initial_line_run: <run>`) を追加。
  - 固定点反復: 壁 → Euler → 初期線 → 壁 を、初期線の変化が許容差以下になるまで (上限 3 回)。
  - case/45 で Hall 版 (現行) との Euler A/B。
- **やらない**:
  - 縮流部の形の変更 (V0 の縮流部に固定。縮流部を設計変数にする案は diagnostician 2 回目で却下)。
  - NS での CFD ピン (Euler の遷音速場を使う。境界層の効果は従来どおり δ_r で別に扱う)。
  - 壁表現 (補間/当てはめ/V4b) の最終決定 — CFD ピン後の点群で verification-m6 §5.1 #16 を決め直す。

## 3. 関連 docs と前提 (観測事実)

- 実現性 (`case/45.isobutane_m6_d155/cfd_initial_line_probe.py` → `_band_ab/cfd_initial_line_probe_run_00{59,62}.json`、既存 Euler 最終場の後処理のみ):
  Euler 場の中でスロート壁点 (x=0) から軸へ C⁻ (dr/dx = tan(θ−μ)、線形補間、RK2、ds 2e-4) をたどれる。軸への着地 x₀ は CFD 0.509 (V0)・0.511 (補間壁) vs Hall 0.520。
  線上の差 (CFD − Hall) は M 最大 0.008 (壁 1.117 vs 1.121、軸 1.214 vs 1.223)、θ 最大 0.045°。
  **CFD データの C⁻ 適合残差の積算は −0.003° (V0)・−0.019° (補間壁)** — Hall の 0.12〜0.13° より 1 桁以上小さい。
- 壁の違い (V0 vs 補間壁) で抽出線が M で ~0.0015 動く → 固定点反復が要る見込み。
- 軸 law の始点 x_A = x₀ と端条件 (M, M′, M″) は現在 Hall の `axis_anchor(x0)`。CFD ピンではこれも CFD 場の軸分布から取る (既存の `x_reach_cfd` / `axis_curve_node` の考え方)。

## 4. 設計方針 (案 — 諮問・plan 段レビュー前)

1. **provider `CFDThroatCharacteristic`** (`design/forge_design/feedback/cfd_initial_line.py` 新設):
   node Euler run の末尾 3 スナップショット平均の (M, θ) を構造格子 (ni × nj) 上で補間し、壁の x=0 の点から C⁻ を軸まで追跡。
   r 等間隔 n_start 点に再標本化 (軸→壁、`HallThroat.throat_characteristic` と同じ規約)。軸 (r=0) は最寄り点群の偶関数 (M)・奇関数 (θ) 当てはめで外挿。
   壁足の θ は 0 (設計スロートの壁接線、`moc_inverse.py:541` と同じ扱い)。`axis_anchor(x0)` と `mach(x, 0)` は軸分布の平滑化スプライン (`axis_curve_node` の evenfit) から。
2. **ガスの整合**: MOC は設計側の semi-perfect ガス (`pm_nu(M, g)`) で ν を作る。CFD の M はその run のガスの音速で定義される — 同じ組成・同じ熱力学なので一致する前提。抽出時に ν↔M の往復誤差を検査する。
3. **固定点**: pass 0 = Hall 版の壁で Euler → 線 L₁ → 壁 W₁ → Euler → 線 L₂ … 。停止: max|ΔM| ≤ 2e-4 かつ max|Δθ| ≤ 0.005° (線上、隣接 pass 間)、上限 3 pass。各 pass は新しい run ディレクトリ。
4. **壁表現**: 下流は位置+壁角の同時当てはめ。始点条件は V0 型 (r′(0)=0、r″(0)=縮流部の端曲率 1/R) を既定にし、CFD 線で始点のオフセットが消えるなら V0 型のまま点上ゲートを満たすはず (§6 V1)。縮流部は V0 と同一。
5. **problem YAML**: `geometry.initial_line: hall` (既定、従来どおり) / `cfd` + `geometry.initial_line_run`。`prepare_info.json` に出所と抽出元 run・snapshot を記録。

## 5. 実装ステップ

1. provider と単体試験 (合成場: Hall 場を格子に載せて抽出し Hall の線を再現できること)。`design/forge_design/feedback/cfd_initial_line.py`、`design/tests/`。
2. `design_chain` の切り替え (`design/forge_design/evaluate/runner_axismach.py`)。
3. case/45 で固定点 pass (Euler, AWS)。
4. Euler A/B と判定。methods/design/overview.md 更新。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の諮問 | diagnostician に設計方針と検証計画を諮る (ユーザは前回 diagnostician を指定) | F |
| 2 | plan 段 codex レビュー | `codex_review.py <this> --stage plan` | O |
| 3 | provider + 単体試験 | §5 ステップ 1。合格: 合成 Hall 場から抽出した線が Hall の線と M 1e-4・θ 0.002° 以内 | O |
| 4 | design_chain 切り替え | §5 ステップ 2。合格: `initial_line: hall` で既存の壁とビット一致 | O |
| 5 | 固定点 pass と Euler A/B | §6 に事前登録した基準で判定 | F |

## 6. 検証 (案 — 諮問で確定させる)

- **V0 (抽出器の検証)**: 合成場 (Hall 場を case/45 の Euler 格子に載せる) から抽出した線が Hall の線を再現 (M 1e-4・θ 0.002°)。
- **V1 (形状)**: CFD ピンの MOC 点群で、最初の区間の自己不整合 c₀ ≤ 0.005° (Hall 版 0.044°)。V0 型始点 (r′(0)=0, r″(0)=0.5) の当てはめが全点ゲート (|Δr| ≤ 5e-6 r_t・|Δθ| ≤ 0.005°) を満たす。
- **V2 (固定点)**: 3 pass 以内に停止条件 (§4.3) を満たす。
- **V3 (Euler)**: CFD ピンの壁 (3 本 + 延長、#15/#16 と同手順・評価器 `eval_wallfit_euler.py --fixed-coef`) を V0・V4b と比較。
  主量: |試験部 P の傾き| η0/η0.1、オーバーシュート η0/η0.1、軸 M の law からの差 (x ∈ [x₀, x_E])。ガード: M 波・P 波・出口コア M が悪化方向に different なら不採用。合否の数値は諮問後に記入。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `design/forge_design/feedback/cfd_initial_line.py` (新規)、`design/forge_design/evaluate/runner_axismach.py` (design_chain・prepare_info)、`methods/design/overview.md`。
- 既定は Hall のまま (`initial_line: hall`) なので既存 case の結果は変わらない。

## 8. 完了条件

- [ ] `methods/design/overview.md` を更新
- [ ] §6 を満たす
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done`、§9 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-10-05` — 起票。動機・実現性は §1・§3 (verification-m6 §9 の V4b 結果、`cfd_initial_line_probe.py`)。ユーザ決定「CFD ピンを前倒し」(「やって」)。
