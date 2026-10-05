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
  Euler 場の中でスロート壁点 (x=0) から軸へ C⁻ をたどれる。軸着地 x₀ は CFD 0.509 (V0)・0.511 (補間壁) vs Hall 0.520。線上の差 (CFD − Hall) は M 最大 0.008、θ 最大 0.045°。
- **訂正 (diagnostician 2026-10-05)**: 「CFD 線は C⁻ 適合が Hall より 1 桁よい」は CPG での積算だけの話。MOC が使う semi-perfect `g` で評価すると積算 +0.023° (CFD) vs +0.135° (Hall)、区間最大 2.3e-4 vs 3.7e-4 rad、幾何残差 (割線−平均角) は CFD 線 4.5e-4 rad > Hall 9.4e-5 rad。
- 抽出器の試作 `case/45.isobutane_m6_d155/cfd_initial_line.py` ((x, η) 平面の 3 次スプライン補間、壁は 3 次補間、RK4): 合成 Hall 場 (case/45 の Euler 格子に Hall を載せる) から Hall の線を M 1.4e-5・θ 0.00026°・x₀ 1e-5 で再現 (`test_cfd_initial_line.py`)。
- 線は物理的には x>0 の壁に依存しない (依存域は上流)。数値的な結合は壁足 1〜2 セルの角度経由のみ: 補間壁 vs V0 (位置差 ≤ 3e-5、角度差 0.087° が x<0.08) で線は ΔM 1.5e-3・Δθ 0.016°。再実行 3 本間 ΔM ≤ 7e-7。
- 現行壁 (V0) は Hall アンカーの軸 law を Euler で最大 1.0e-3 (x=32)、スロート近傍 ≤ 8e-4 で実現している。CFD 軸の x₀ 局所 4 次フィットの M″ は 0.052〜0.060 (Hall 0.0214 の 2.5〜2.8 倍、壁の 3e-5 差で 15 % 動く)。
- m* (`moc_inverse.py:541`、線上の `_flux_along`) は CFD 線で Hall 比 +4.1e-4 → 下流壁が r·2.0e-4 膨らみ出口半径 0.775 m で +0.16 mm。+4.1e-4 が Hall 打切りか Euler 離散化かは格子細分でしか決まらない。
- CFD 0 step の pass 1 予測 (diagnostician の再計算): 始点オフセット a_θ 0.096° → 0.047° (V0 場) / 0.055° (補間壁場)、c₀ 0.044° → 0.009°。CFD 軸そのものを目標にした再現試験でも a_θ 0.058° が残る → **約 0.05° は MOC の壁足第 1 セルの性質で Hall 線由来ではない**。

## 4. 設計方針 (diagnostician 2026-10-05 を採用)

1. **provider** (`HallThroat` を継承し `throat_characteristic` / `axis_anchor` / `mach` だけ上書き): node Euler 場の (M, θ) を (x, η=r/r_w) 平面で補間 (壁は断面の壁節点の 3 次補間)。
   壁足は (0, 1, θ=0) に厳密に置き、追跡は壁足 1e-6 内側から RK4 で軸まで。壁足の M は内側 10 点の線形外挿。軸端は追跡点の r ≤ 0.05 の 2 次多項式 (偶関数当てはめは使わない — C⁻ は軸を斜めに横切る)。r 等間隔 n_start 点に再標本化。
2. **軸アンカー**: M_A = 線の軸端の M (線と同じ出所)、M′_A = 軸 evenfit の 4 次窓フィット (窓 0.25)、**M″_A は Hall の `axis_anchor` を据え置く** (CFD 局所 M″ は law を実測軸から 8.5e-3 離す)。x_A = 線の軸着地 x₀ (Hall 0.520 → 0.509)。L_c 固定で x_E は x_A に従う (x_F 不変)。
3. **m***: CFD 線から `_flux_along` で (Hall 比 +4.1e-4)。A/B は r_t 固定、生産化では出口半径から r_t を解き直す (別項)。
4. **凍結 + 検証 1 回** (固定点反復はしない): pass 0 = V0 の Euler 場 (run_0062、3 本の差 7e-7) から線・m*・(M_A, M′_A) を取り凍結 → pass 1 の壁 → Euler → 再抽出で D₁ = |L(pass1) − L(pass0)| を測る (§6 V2)。上限 2 pass、超過は諮問、D₁ > D₀ (= |L(V0) − L_Hall|) は棄却。
5. **壁表現**: 縮流部は V0 と同一、下流は V0 型の同時当てはめ (r′(0)=0, r″(0)=1/R)。始点オフセットは CFD ピンで半減するが消えない前提 (点上ゲートは i ≥ 1)。
6. **ガス**: C⁻ 適合残差は MOC の semi-perfect `g` で評価・記録。ν↔M の往復誤差 ≤ 2e-4 を検査。
7. **problem YAML**: `geometry.initial_line: hall` (既定) / `cfd` + `geometry.initial_line_run` + `initial_line_res`。`prepare_info.json` に出所と snapshot を記録。

## 5. 実装ステップ

1. provider と単体試験 (合成場: Hall 場を格子に載せて抽出し Hall の線を再現できること)。`design/forge_design/feedback/cfd_initial_line.py`、`design/tests/`。
2. `design_chain` の切り替え (`design/forge_design/evaluate/runner_axismach.py`)。
3. case/45 で固定点 pass (Euler, AWS)。
4. Euler A/B と判定。methods/design/overview.md 更新。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~§4・§6 の諮問~~ | 完了 2026-10-05 (diagnostician、§6.1)。判断: 凍結 + 検証 1 回、M″_A は Hall、V0 主対照、|傾き| 両側判定 | F |
| 2 | plan 段 codex レビュー | `codex_review.py <this> --stage plan` | O |
| 3 | provider を §4.1/4.2 に合わせて改修 + V0 試験 | 試作 `case/45.isobutane_m6_d155/cfd_initial_line.py` を壁足厳密・1e-6 内側開始・軸端 2 次多項式・アンカー (M_A=線の軸端、M′_A evenfit、M″_A Hall) に。合格 §6 V0 | O |
| 4 | pass 1 の MOC (CFD 0 step) | V0 場 run_0062 から凍結した線で design_chain → §6 V0b・V1。r″ の山が下がらなければ止まる | O |
| 5 | Euler A/B (V3) と D₁ (V2) | AWS。run_0077〜0079 + 延長 0080〜0082 | F |
| 6 | design_chain への組み込み・methods 更新 | `initial_line: hall` でビット一致 | O |
| 7 | 生産化 (r_t の解き直し)・格子細分での m* 帰属 | 別項 | F |

## 6. 検証 (事前登録、diagnostician 2026-10-05 の数値)

- **V0 (抽出器)**: 合成 Hall 場から 線 max|ΔM| ≤ 2e-4・max|Δθ| ≤ 0.004°・|Δx₀| ≤ 2e-4、アンカー |ΔM| ≤ 2e-5・|ΔM′| ≤ 1e-4・|ΔM″| ≤ 1e-3、m* の再現 ≤ 1e-4 (相対)、ν↔M 往復 ≤ 2e-4。
- **V0b (アンカー整合)**: law − 実測軸 (V0 場) ≤ 2e-3 on [x₀+0.3, x_K] (予測 1.2e-3)。
- **V1 (形状、CFD 0 step)**: c₀ ≤ 0.015° (予測 0.009°; Hall 0.044°)、a_θ ≤ 0.06° (予測 0.047°)、点上ゲート i ≥ 1 で |Δr| ≤ 5e-6 r_t・|Δθ| ≤ 0.005°、V0 型当てはめの r″ の山 (x∈[0,0.3]) が V0 の 0.534 より下がる (予測 ≤ 0.52)。併記: 壁の変化 (pass 1 − V0、m* 分 r·2.05e-4 を除いた形状) 予測 −2e-4〜−1e-3 r_t on [0,15)、m* +4.1e-4、r_F +2.0e-4。**r″ の山が下がらなければ Euler に進まない**。
- **V2 (凍結の自己整合)**: D₁ の max|ΔM| ≤ 5e-4・max|Δθ| ≤ 0.01°・|Δm*/m*| ≤ 2e-4 → pass 1 を採用。超過は諮問 (反復しない)、D₁ > D₀ は棄却。
- **V3 (Euler)**: CFD ピン壁 3 本 + 延長 (#15/#16 と同手順、評価器 `eval_wallfit_euler.py --fixed-coef`)。**主対照 V0** (run_0056〜0058 + 0062〜0064; 変数は初期線・m*・M_A だけ)、V4b は副対照 (参考)、補間壁は記録のみ。4 量 (|P 傾き| η0/η0.1、オーバーシュート η0/η0.1) を両側判定: 4 量とも q_pin − q_V0 ≤ −U → 改善で採用 / いずれか ≥ +U → 悪化で不採用 / それ以外 → 差なし (V1 の形状で採否)。U は #15 の式、傾きの Δ_slope = 0.03 %pt。傾きの符号が反転したら |·| が減っても「反転」と記録 (改善と呼ばない)。m* のレベルシフト (試験部 M +O(1e-4)) があるので、オーバーシュートは 100(M/M_exit−1) の最大も併記。出口コア M が下がる方向に different なら m* ピンの符号が逆 → 棄却。予測: |Δ傾き| 0.05〜0.15 %pt (符号未定)、|Δ傾き| < U なら効きは V4b の blend 区間 (x<0) 側と解釈し、CFD ピンの効果は m*・出口 M に限ると書く。
- **不足 (別項)**: 線のスナップショット間変動 (AWS の末尾 5 枚)、Euler 格子細分での線と m* の変化 (m* +4.1e-4 の帰属)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-05` | diagnostician (Fable) — 本 plan §4・§6 の諮問 (ユーザ指定の諮問先) | 方針は半分よい: M″_A を CFD 局所フィットから取るのは不可、固定点の許容差 2e-4 反復は無意味 (抽出床以下・角度利得 0.5)、始点オフセットは消えない (0.047° 残る)、CFD 線の C⁻ 適合の優位は CPG 積算だけ、m* +4.1e-4 で出口 +0.16 mm、抽出器の偶関数外挿は誤り | 全件採用 → §3 訂正・§4・§5.1・§6 書き換え |

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
- `2026-10-05` — diagnostician の判断を採用し §3 (訂正)・§4・§5.1・§6 を書き換え。抽出器の試作は合成 Hall 場で合格 (M 1.4e-5・θ 0.00026°)。
