# スロート直後の壁 r″ を単調にする当てはめ (joint 壁に r‴ ≤ 0 の不等式拘束)

## メタ

- **area**: `tooling`
- **status**: `draft`
- **related_docs**:
  - `methods/design/overview.md` (設計区間の壁表現 `geometry.wall_repr` の `joint`)
  - `case/45.isobutane_m6_d155/README.md`
- **related_plans**:
  - [tooling-nozzle-cfd-pinned-initial-line.md](../accepted/tooling-nozzle-cfd-pinned-initial-line.md) (joint 壁 + CFD ピンの生産化。本 plan はその壁表現の変更)
  - [verification-m6-axis-wave-mesh-su2.md](verification-m6-axis-wave-mesh-su2.md) (§5.1 #15・#16、§9 2026-10-05: 補間壁 / joint V0 / V4 / V4b の Euler A/B と判定の枠組み)
- **created**: `2026-10-06`
- **owner**: `Claude (Opus 5.5)`

## 1. 目的

最終壁 (`problem_d155_ns_finemesh_recal_final.yaml`、run_0117/0118 の壁) では、スロート直後の x = 0.014 r_t で r″ が 0.5102 になり、1/R = 0.5 を超える山がある (幅 0.034 r_t)。
ユーザ (2026-10-06):「流れ場にあまり影響はないと思うが、結構気持ち悪い」。説明ページ (Artifact「スロート r″ の山」) を受けて **案 B (r″ 単調拘束) で進める** と決定した。
joint 壁の当てはめに「[0, 1.5] r_t で r″ が増えない」不等式拘束を加えて山を消す。そのうえで、Euler で現行壁に対して非劣化であることを確かめてから生産の既定にする。

## 2. スコープ

- **やる**:
  - `joint_fit_wall` / `JointFitCFDWall` に単調拘束のオプション `mono_r2=(a, b)` を追加する。未指定ならビット同一。
  - `design_chain` の `geometry.wall_fit_mono_r2` キーから渡す。
  - 形状ゲートを課す。
  - Euler A/B を行う (生産 Euler 格子 G1、各腕 3 回)。
  - 合格したら生産問題 YAML にキーを入れ、`methods/design/overview.md` を更新する。
- **やらない**:
  - MOC 第 1 点が 1/R より強く曲がる原因 (CFD ピン初期線の誤差か離散化か) の追及。必要なら別 plan とする。
  - 上流 Hermite の端条件の変更 (V4 / V4b 型)。これは試験部の P 傾きを 0.3 %pt 動かしたため不採用 (verification-m6 §9)。
  - NS の再計算。壁の変化は最大 0.5 µm で、スロートの δ_r 107 µm の 0.5 %。Euler の結果を見て §6 の条件に当たったときだけ別途判断する。
  - 軸上 M の +0.2 % (verification-m6 §5.1 #17、別セッション)。

## 3. 関連 docs と前提

形状だけの事前検討。スクリプトは `case/45.isobutane_m6_d155/throat_r2_explainer.py`、出力 `throat_r2_explainer.json`。いずれも commit c821b71e。

- **原因**:
  - MOC 壁点の区間平均曲率 (tan θ の差 / Δx) は、第 1 区間 (x 0〜0.025) だけ 0.515 で、以降は 0.484, 0.475, 0.464 … と下がる。
  - 当てはめは x = 0 で r″ = 1/R (上流 Hermite と C²) に固定されているため、第 1 点の流れ角 0.7452° に寄るには 0.5 を超えるしかない。
  - 第 1 点の流れ角は、曲率 1/R の円弧の角度 atan(x₁/R) = 0.7233° より 0.0219° 大きい。
- **山の高さは λ 依存**: λ 1e-10 / 1e-9 / 1e-8 で、r″ の最大は 0.5212 / 0.5102 / 0.5028、第 1 点の流れ角ずれは 0.007 / 0.012 / 0.019° (現行の始点拘束のまま)。
- **試し当てはめ C** (始点拘束は現行と同じ。r‴ の B-spline 係数のうち台が [0, 1.5] にかかるものを ≤ 0 に拘束。有効制約法):
  - 全 λ で r″ ≤ 0.5。
  - 第 1 点の流れ角ずれは 0.0219 / 0.0219 / 0.0220°、半径ずれは 1.04 / 1.05 / 1.11 µm (r_t = 76.654 mm)。
  - 現行との差 (λ 1e-9) は、半径で最大 0.50 µm (x ≈ 0.06)、壁角で最大 0.0105°。差が 0.05 µm を超えるのは x 0.017〜0.12 だけで、x > 0.3 では 1e-3 µm。
  - r″ が平らな区間 (x < 0.03) から下がり始める所で、r‴ が 0 → −0.9 と変わる (x 0.04〜0.06)。λ 1e-8 ではこれが緩い。
  - λ 1e-9 での有効制約は 4 本、解き直しは 5 回。
- **物理壁**: `PhysicalNozzleWall` は joint 壁に δ_r を足す解析経路なので、設計壁の山がそのまま残る (現行の物理壁の r″ 最大は 0.5114、x 0.014)。設計壁の山を消せば物理壁の山も δ_r″ 分だけになる見込み (§6 S6 で確認)。
- **流れへの感度の参考**: V4 では始点の角度を 0.09° 変えたところ (今回の壁角差の約 9 倍)、試験部の P 傾き・オーバーシュートが Δq を超えて動いた (交絡あり)。今回の差がそれより十分小さいとは事前には言えないので、Euler で測る。

## 4. 設計方針

### 4.1 当てはめ

- `geometry/wall_axismach.py::joint_fit_wall(…, mono_r2=None)`。`mono_r2=(a, b)` のとき、次の不等式を課す。
  - 対象は r‴ (2 次スプライン、`BSpline(t, e_i, 5).derivative(3)` の係数を cᵢ の一次式として組む) の係数のうち、台 [t₃ⱼ, t₃ⱼ₊₃] が (a, b) にかかるもの。これを ≤ 0 にする。
  - B-spline は係数の凸包に収まるので、これは r″ が [a, b] で単調非増加であることの**十分条件**になる (必要条件ではない。少し保守的)。
- **解法**: 目的関数・ノット・等式拘束 (始点 r, r′, r″、出口 r, r′) は現行のまま。目的関数の 2 次形式は対角の最大値で正規化する。
  - 有効制約法で解く。違反最大の制約を追加し、乗数が負の制約を外す。
  - 上限は 200 回で、超えたら例外にする。収束判定は違反 ≤ 1e-10 × max|G c|。
  - `fit_diag` に `mono_r2`・有効制約数・反復数を記録する。
- **既定値**: `mono_r2=None` なら現行と同じ連立方程式 1 回で、**ビット同一** (正規化も行わない経路)。
- **区間 [0, 1.5]**: r″ はスロート下流で単調に下がり、x ≈ 2.15 で浅い極小を持つ (現行・試し C とも。形状固有)。その手前で止めるため 1.5 とする。λ は生産値 1e-9 のまま。

### 4.2 設計チェーン

- `runner_axismach.design_chain`: `wall_repr == "joint"` のとき、`geometry.wall_fit_mono_r2` (2 要素のリスト、無ければ None) を `JointFitCFDWall(…, mono_r2=…)` に渡す。
- `design_chain` の戻り値の `wall_fit` に記録する。
- `prepare` / `prepare_ns` の `prepare_info.json` に残る経路を確認する。

### 4.3 判定の流れ

形状ゲート (S, CFD なし) → Euler A/B (E) → 合格なら生産キー化。

- S が FAIL なら Euler に進まず諮問する。
- E の判定は「非劣化」(現行と比べて悪化していないことを不確かさ込みで示す)。改善は主張の条件にしない。山を消すのは形状上の要求 (ユーザ決定) で、流れの改善は期待していない。

## 5. 実装ステップ

1. `design/forge_design/geometry/wall_axismach.py`: `joint_fit_wall` と `JointFitCFDWall` に `mono_r2`。
2. `design/forge_design/evaluate/runner_axismach.py`: `geometry.wall_fit_mono_r2` を渡して記録する。
3. `design/tests/run_joint_fit_mono_tests.py` (新規): 既定のビット同一、単調性、不正入力の拒否。
4. `case/45.isobutane_m6_d155/throat_mono_shape_gate.py` (新規): §6 S1〜S6 → `_band_ab/throat_mono_shape_gate.json`。
5. Euler A/B の入力と評価器を一般化し (§5.1 #5)、AWS で回す。
6. 合格時: 生産問題 YAML にキーを追加し、`methods/design/overview.md`・case README・本 plan を更新する。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の諮問 | codex diagnose (`~/.config/forge/diagnose-backend` = codex)。ブリーフ `notes/reviews/briefs/2026-10-06-throat-monotone-r2.md` | F |
| 2 | codex plan 段レビュー | `codex_review.py <本 plan> --stage plan` (run_in_background) | O |
| 3 | 当てはめの実装 | §5 の 1〜3。合格条件: 新テストが FAIL 0。既存 `design/tests/run_mesh_params_tests.py`・`run_cfd_initial_line_tests.py`・`run_physical_wall_analytic_tests.py` が FAIL 0 | O |
| 4 | 形状ゲート S | §5 の 4。合格条件: §6 S1〜S6 が全 PASS | O |
| 5 | Euler A/B の準備 | 次の 3 点。<br>・問題 `problem_d155_euler_pin_G1_recal_mono.yaml` (G1_recal + `wall_fit_mono_r2: [0.0, 1.5]`)。<br>・`eval_wallfit_euler.py` の一般化: X_E・X_F を各腕の `prepare_info.json` / 設計壁から読む (現状は c2final n2400 の定数)、run 名の glob を 4 桁に。<br>・投入スクリプト `run_throat_mono_ab.sh`。IC は §6 E1 | O |
| 6 | Euler A/B の実行と判定 | §6 E1〜E4。各 run に `check_quasisteady` / `check_convergence`。README の run 一覧に追記 | O (判定の解釈は F) |
| 7 | 生産化 | E 合格時: `problem_d155_ns_finemesh_recal_final*.yaml`・`problem_d155_euler_pin_G1_recal.yaml` 系にキー。`methods/design/overview.md`。報告 (pptx) の壁は再生成しない旨を §9 に | O |
| 8 | codex result 段レビュー | `--stage result` | O |

## 6. 検証

### 6.0 事前登録 (2026-10-06、諮問の前の初稿。諮問で変えたら §9 に残す)

**形状ゲート S** (最終問題の MOC 点群、λ = 1e-9、`mono_r2 = (0, 1.5)`):

- **S1 単調性**: [0, 1.5] で r‴ の B-spline 係数がすべて ≤ 1e-12。加えて 1e-4 刻みの評価で r″(x) ≤ 1/R + 1e-9。
- **S2 接続**: x = 0 での上流 Hermite との跳び。r、r′、r″ が各 ≤ 1e-8。
- **S3 局所性**: x ≥ 0.3 で、現行壁との差が |Δr| ≤ 1e-7 r_t かつ |Δ壁角| ≤ 1e-4°。
- **S4 点上の忠実度**: 第 1 点 (i = 1) の流れ角ずれは、下限 θ₁ − atan(x₁/R) (= 0.0219°) との差 ≤ 0.002°。i ≥ 2 は全点で |Δθ| ≤ 0.02°、|Δr| ≤ 2e-5 r_t。
- **S5 λ 頑健性**: λ ∈ {1e-10, 1e-8} でも S1 が成り立ち、第 1 点の流れ角ずれの幅 ≤ 0.001°。
- **S6 物理壁**: run_0117 の δ_r で `PhysicalNozzleWall` を作り、[0, 0.3] の r″ の最大 ≤ 1/R + max|δ_r″| + 1e-6。`validate()` が空であること。
- **S7 既定のビット同一**: キー無しの設計壁が現行と max|Δr| = 0。

**Euler A/B E** (生産 Euler 格子 G1 = `problem_d155_euler_pin_G1_recal.yaml` の mesh、バイナリ `~/forge-wallfit-bin`、AWS):

- **E1 腕と IC**:
  - 腕 A (現行壁): run_0140〜0142 `euler_wallfit_pinG1_r{1,2,3}`。
  - 腕 B (単調壁): run_0143〜0145 `euler_wallfit_monoG1_r{1,2,3}`。
  - 両腕とも IC は run_0114 の最終場。A は同一メッシュなので `restart_field`。B は節点座標が最大 0.5 µm 動くので `interp_field`。
  - 段 soft (1 次 cfl 0.5、3000) → 本段 2 次 cfl 2・relax 0.7・18000 step・1000 ごと出力 (E2 の 12000 + 延長 6000 を 1 本にしたもの)。
  - B の壁が使われた証拠として、変換後メッシュの壁節点から求めた r″ の山の有無を記録する。
- **E2 量と許容悪化幅 Δq** (verification-m6 #15・CFD ピン plan V3 と同じ定義・評価器 `eval_wallfit_euler.py --fixed-coef`):
  - |P 傾き| η0 / η0.1: 0.03 %pt
  - オーバーシュート η0 / η0.1、出口コア M で規格化したオーバーシュート: 0.003 %pt
  - M 波: 0.001 %pt
  - P 波: 0.010 %pt
  - |出口コア M − 6|: 0.00018
- **E3 判定**:
  - 不確かさ U = max(3R_A, 3R_B, 2T, 2E_fixed, Δq/10)。R は腕内 3 回の幅、T は末尾 5 枚の幅、E は評価刻み感度。
  - **採用**: 全量で (q_B − q_A) + U ≤ Δq。
  - **不採用**: いずれかで q_B − q_A ≥ Δq + U。
  - それ以外は保留として諮問。η0 は T ≤ Δq/4 のときだけ判定に入れる (V3 と同じ)。
- **E4 準定常**: 各 run の評価量に `check_quasisteady --series-csv`。末尾 5 枚 STEADY でない量は判定に入れず保留扱い。`check_convergence` は全残差 (plateau NOT CONVERGED は既知、記録のみ)。
- **E5 予測** (外れたら止めて諮問): 出口コア M の差 |q_B − q_A| ≤ 1e-4 (較正のやり直しは不要)。軸の M・P の差は x < 2 に集中し、x > 10 で 1e-4 未満。

**生産化の条件**: S 全 PASS かつ E3 採用。NS は再計算しない。根拠は、壁変化 ≤ 0.5 µm (δ_r 107 µm の 0.5 %) と Euler 非劣化の 2 つ。これを限界として §9 と README に書く。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `design/forge_design/geometry/wall_axismach.py`、`design/forge_design/evaluate/runner_axismach.py` (オプション追加。既定はビット同一)。
- `case/45.isobutane_m6_d155/` の評価器・問題 YAML・投入スクリプト。
- `methods/design/overview.md` (joint の節)。

## 8. 完了条件

- [ ] `methods/design/overview.md` の joint の節を更新済み
- [ ] §6 S・E を満たす (または判定保留を記録して諮問済み)
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done` にして accepted へ移動し、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-06` — 初稿。ユーザ決定「B (r″ 単調拘束) で進めてよい」。形状の事前検討は commit c821b71e (`throat_r2_explainer.py`)。
