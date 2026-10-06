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
  - (候補形状のあいだは) NS の再計算。**生産化する (既存 NS 性能・凝縮結果を新形状に引き継ぐ) なら NS の再評価が要る** (諮問 2026-10-06 ⑥: NS 出口コア M 5.998887 の下限 5.9988 までの余裕 8.7e-5 は Euler の Δq 1.8e-4 より小さく、形状変位/δ_r の比からは性能感度を上限評価できない)。§6 N で事前登録する。
  - 軸上 M の +0.2 % (verification-m6 §5.1 #17、別セッション)。

## 3. 関連 docs と前提

形状だけの事前検討。スクリプトは `case/45.isobutane_m6_d155/throat_r2_explainer.py`、出力 `throat_r2_explainer.json`。いずれも commit c821b71e。

- **原因**:
  - MOC 壁点の区間平均曲率 (tan θ の差 / Δx) は、第 1 区間 (x 0〜0.025) だけ 0.515 で、以降は 0.484, 0.475, 0.464 … と下がる。
  - 当てはめは x = 0 で r″ = 1/R (上流 Hermite と C²) に固定されているため、第 1 点の流れ角 0.7452° に寄るには 0.5 を超えるしかない。
  - 第 1 点の流れ角は、曲率 1/R 一定の曲線 r = 1 + x²/2R の角度 atan(x₁/R) = 0.7233° より 0.0219° 大きい。
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
- **代償の下限 (諮問 ①)**: r′(0) = 0 と r″ ≤ 1/R を積分すると r′(x₁) ≤ x₁/R なので、第 1 点の流れ角ずれは解析的に**少なくとも θ₁ − atan(x₁/R) = 0.0219°** (atan(x/R) は r″ = 1/R 一定の放物線 r = 1 + x²/2R の接線角)。点の重みを下げても、局所的に λ を強めても、この下限は消えない (局所 λ では単調性も保証できない)。不等式拘束を直接課す。
- **停止条件・許容誤差の統一 (諮問 Minor)**: 2 次形式は対角最大で正規化し、不等式行列 G は行ごとに最大絶対値で正規化する。違反判定は正規化後の max(G c) ≤ 1e-10 · max(1, max|G c|)、乗数は ≥ −1e-12 · max(1, max|μ|)。`fit_diag` に等式残差・正規化後の不等式最大値・乗数最小値・有効制約の台を記録する。S1 は同じ正規化で判定する。
- **平滑化項の積分 (諮問 ②・判別 A/B, 2026-10-06)**: 現行の ∫(r‴)² は全長 8000 点の和 (刻み 0.0119 ≈ 最初のノット間隔 0.0125) で、スロート近傍では基底ごとの重みが厳密値と最大 4.5 倍ずれる。判別 A/B (§9 同日) で単調拘束版は積分方法に鈍感 (D₃ 2.4 %、D₄ 5.1 % ≤ 10 %) と確認したので、**積分は現行のまま (既定のビット同一を保つ)**。現行 V0 のほうが積分に敏感 (山 0.5102 → 0.5082) なことは既知の課題として記録し、本 plan では変えない。

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
| 1 | ~~§4・§6 の諮問~~ 完了 | 判断: 2026-10-06 codex (diagnose) — 単調拘束は採用、Euler 前に積分方法の判別 A/B (→ 実施、積分主因説を棄却)。Major 5・Minor 1 を全件採用し §2・§4.1・§6 に反映 (§6.1) | F |
| 2 | codex plan 段レビュー | `codex_review.py <本 plan> --stage plan` (run_in_background) | O |
| 3 | 当てはめの実装 | §5 の 1〜3。合格条件: 新テストが FAIL 0。既存 `design/tests/run_mesh_params_tests.py`・`run_cfd_initial_line_tests.py`・`run_physical_wall_analytic_tests.py` が FAIL 0 | O |
| 4 | 形状ゲート S | §5 の 4。合格条件: §6 S1〜S8 が全 PASS (S4・S5 は記録項目を含む) | O |
| 5 | Euler A/B の準備 | 次の 3 点。<br>・問題 `problem_d155_euler_pin_G1_recal_mono.yaml` (G1_recal + `wall_fit_mono_r2: [0.0, 1.5]`)。<br>・`eval_wallfit_euler.py` の一般化: X_E・X_F を各腕の `prepare_info.json` / 設計壁から読む (現状は c2final n2400 の定数)、run 名の glob を 4 桁に。<br>・投入スクリプト `run_throat_mono_ab.sh`。IC は §6 E1 (B の IC 写像の記録を含む)。<br>・出口コア M は両腕とも同じ評価器・同じ標本 (§6 E2) | O |
| 6 | Euler A/B の実行と判定 | §6 E1〜E4。各 run に `check_quasisteady` / `check_convergence`。README の run 一覧に追記 | O (判定の解釈は F) |
| 7 | NS 再評価 | E 採用時: §6 N の NS (および凝縮) を回して判定。合格まで生産キーは入れない (候補形状) | O (判定の解釈は F) |
| 7b | 生産化 | N 合格時: 生産問題 YAML にキー、`methods/design/overview.md`、報告 (pptx) を新しい NS で再生成 | O |
| 8 | codex result 段レビュー | `--stage result` | O |

## 6. 検証

### 6.0 事前登録 (初稿 2026-10-06 → 同日諮問で改訂。初稿との差分は §9)

**形状ゲート S** (最終問題の MOC 点群、λ = 1e-9、`mono_r2 = (0, 1.5)`。S4・S5 の数値は試し C を見た後の値なので、**独立の裏付けではなく開発用の受入条件**として扱う。諮問 ③):

- **S1 単調性**: [0, 1.5] で正規化後の r‴ の B-spline 係数がすべて ≤ 1e-10 (§4.1 の停止条件と同じ正規化)。加えて 1e-4 刻みの評価で r″(x) ≤ 1/R + 1e-9。
- **S2 接続**: x = 0 での上流 Hermite との跳び。r、r′、r″ が各 ≤ 1e-8。
- **S3 局所性**: x ≥ 0.3 で、現行壁との差が |Δr| ≤ 1e-7 r_t かつ |Δ壁角| ≤ 1e-4°。
- **S4 点上の忠実度**: 第 1 点の流れ角ずれは、§4.1 の解析的下限 0.0219° との差 ≤ 0.002°。i ≥ 2 は**記録のみ** (要求に伴う忠実度の損失として明示。試し C: 最大 0.0139°・1.37e-5 r_t、現行 0.0087°・7.8e-6 r_t)。
- **S5 λ 感度** (記録): λ ∈ {1e-10, 1e-8} で S1 が成り立つこと (合否)。全域の Δr・Δθ、[0, 0.3] の r″・r‴・r⁗ の差を記録する。
- **S6 物理壁 (生産経路)**: `prepare_ns` と同じ経路 (新しい設計壁から `integral_bl` を再計算 → 平滑化 → `delta_r_from_table` → `PhysicalNozzleWall`、pw_ramp [−11, −6]) で物理壁を作る。合否は次の 3 つ。
  - [0, 1.5] に r″ の増加区間 (r‴ > 1e-9) が無い。または、増加が δ_r″ の寄与だけで説明でき、その大きさを記録する。
  - [0, 3] の r″ の内部極値の数 ≤ 現行の物理壁。
  - `validate()` が空。
  - 加えて、現行の物理壁との全域の形状差 (Δr の最大と位置) を記録する。固定 δ_r (run_0117) での検査は S1 からほぼ自明なので行わない (諮問 S6)。
- **S7 既定のビット同一**: キー無しの設計壁が現行と max|Δr| = 0。
- **S8 滑らかさ** (諮問 ②): [0, 0.3] で、区間ごとに厳密に求めた ∫(r‴)² と max|r⁗| が、どちらも現行壁 (同じ厳密評価) 以下。悪化した量があれば「滑らかさ改善」とは呼ばず、記録して諮問する。判別 A/B 時点の値は、試し C 0.0950 / 58.3、現行 0.1059 / 225。

**Euler A/B E** (生産 Euler 格子 G1 = `problem_d155_euler_pin_G1_recal.yaml` の mesh、バイナリ `~/forge-wallfit-bin`、AWS):

- **E1 腕と IC**:
  - 腕 A (現行壁): run_0140〜0142 `euler_wallfit_pinG1_r{1,2,3}`。
  - 腕 B (単調壁): run_0143〜0145 `euler_wallfit_monoG1_r{1,2,3}`。
  - 両腕とも IC は run_0114 の最終場。A は同一メッシュなので `restart_field`。B は節点座標が最大 0.5 µm 動くので `interp_field` (規則どおり。両腕への interp_field や座標検査の緩和はしない)。
  - **B の IC 写像を記録する**: 対応節点の一致率・最大移動量、保存量の差 (最大・RMS)。IC の影響は消去済みとは扱わず、限界として記録する (諮問 ④)。
  - 段 soft (1 次 cfl 0.5、3000) → 本段 2 次 cfl 2・relax 0.7・18000 step・1000 ごと出力。
  - B の壁が使われた証拠として、変換後メッシュの壁節点の座標と設計壁 (`wall_design.csv`) の差を記録する。r″ の山の有無は壁節点 2〜3 個からは認定しない (解析形で S1 が担う)。
- **E2 量と許容悪化幅 Δq** (verification-m6 #15・CFD ピン plan V3 と同じ定義・評価器 `eval_wallfit_euler.py --fixed-coef`):
  - |P 傾き| η0 / η0.1: 0.03 %pt
  - オーバーシュート η0 / η0.1、出口コア M で規格化したオーバーシュート: 0.003 %pt
  - M 波: 0.001 %pt
  - P 波: 0.010 %pt
  - |出口コア M − 6|: 0.00018
  - **出口コア M の定義** (諮問 ⑤): 判定は両腕とも `eval_wallfit_euler.py` の η 重み付き面積平均で行う。較正側の定義 (`euler_grid_ab.py` の G1 帯内節点の共通標本 M_common) も両腕について併記する。2 つの定義を「較正済み M = 6」として混ぜない。出口の評価誤差として、2 定義の差の腕間変化と、出口断面の標本 (帯 η 0.05〜0.7 の端を ±0.05 動かす) 感度を U に入れる。
- **E3 判定**:
  - 不確かさ U = max(3R_A, 3R_B, 2T, 2E_fixed, E_exit, Δq/10)。R は腕内 3 回の幅、T は末尾 5 枚の幅、E は評価刻み感度、E_exit は E2 の出口評価誤差 (出口 M のみ)。3R は再実行ばらつきの目安で、保証された信頼限界ではない (諮問 ⑤)。
  - **採用 (候補として非劣化)**: 全量で (q_B − q_A) + U ≤ Δq。
  - **不採用**: いずれかで q_B − q_A ≥ Δq + U。
  - それ以外は保留として諮問。η0 は T ≤ Δq/4 のときだけ判定に入れる (V3 と同じ)。
- **E4 準定常**: 各 run の評価量に `check_quasisteady --series-csv`。末尾 5 枚 STEADY でない量は判定に入れず保留扱い。`check_convergence` は全残差 (plateau NOT CONVERGED は既知、記録のみ)。
- **E5 予測** (外れたら止めて諮問): 出口コア M の差 |q_B − q_A| ≤ 1e-4。軸の M・P の差は x < 2 に集中し、x > 10 で 1e-4 未満。

**NS 再評価 N** (E 採用時のみ、生産化の前提。諮問 ⑥):

- 問題は `problem_d155_ns_finemesh_recal_final.yaml` + `wall_fit_mono_r2: [0.0, 1.5]`。r_t・k_f・Md_moc_offset は run_0117 と同じで、較正し直さない。物理壁は生産経路 (S6)。
- run: run_0146 `ns_mono_final`。IC は run_0117 最終場を `interp_field`。段階起動 → cfl 1・60000 step (run_0117 と同じ手順)。
- ゲートは run_0117 と同じ: 出口半径 0.775 ± 0.1 mm、δ_E/δ_C(x_F) 1 ± 0.5 %、出口コア M 6.000 ± 0.02 %、波 η0.1 ≤ 0.01 %、オーバーシュート η0.1 ≤ +0.035 %、壁解像 PASS。各量に check_quasisteady、全残差の check_convergence。
- run_0117 との差 (出口コア M・δ_E・波・オーバーシュート) を表にする。
- 出口コア M が下限を割った場合は、旧壁と同様にユーザに判断を仰ぐ。
- 凝縮 ON (run_0147 `ns_mono_final_cond`、run_0118 と同手順) はゲート合格後に回す。

**生産化の条件**: S 全 PASS かつ E3 採用かつ N のゲート合格 (またはユーザ判断)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-06` | [`notes/reviews/2026-10-06-throat-monotone-r2-diagnose.md`](../../notes/reviews/2026-10-06-throat-monotone-r2-diagnose.md) (ブリーフ [`briefs/2026-10-06-throat-monotone-r2.md`](../../notes/reviews/briefs/2026-10-06-throat-monotone-r2.md)) | 方針採用・C0/M5/m1 + 判別 A/B 1 件 | 全件採用: ① 第 1 点の解析的下限を §4.1 に、② S8 (厳密 ∫(r‴)²・max\|r⁗\|) を追加、③ S4・S5 を開発用受入条件 + 記録に格下げ、④ B の IC 写像を記録し IC 影響は限界として扱う、⑤ 出口 M の定義を両腕同一評価器に固定し評価誤差を U に、⑥ 生産化には NS 再評価 (§6 N)、S6 を生産経路の物理壁に、Minor 停止条件の統一。判別 A/B は実施 (§9) |

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
- `2026-10-06` — **codex (diagnose) に諮った**: `notes/reviews/2026-10-06-throat-monotone-r2-diagnose.md` — 単調拘束は採用、ただし Euler 前に平滑化項の積分方法だけを変える判別 A/B を行う。採否は §6.1 (全件採用)。§6.0 の初稿からの変更: S1 の正規化、S4・S5 を開発用受入条件 + 記録に、S6 を生産経路の物理壁に、S8 (滑らかさ) 追加、E1 に IC 写像の記録、E2 に出口 M の定義固定と評価誤差、E3 の U に E_exit、生産化に N (NS 再評価) を追加。初稿の「NS は再計算しない」は撤回。
- `2026-10-06` — **判別 A/B (CFD 0 step、事前閾値 D₃/D₄ ≤ 0.10)**: `case/45.isobutane_m6_d155/throat_mono_integration_ab.py` → `_band_ab/throat_mono_integration_ab.json`。単調拘束版で、積分を全長 8000 点 (現行) と区間 3 点 Gauss (厳密) に変えた。結果は次のとおり。
  - D₂ 0.0002・D₃ 0.024・D₄ 0.051 で、**積分誤差主因説を棄却**。急変は拘束に伴う形状上の代償として評価する。
  - 全域の差は Δr 1.1e-8 r_t・Δθ 7.6e-5°。
  - 有効制約は 4 本で、台はすべて [0, 0.05] 内。KKT は等式残差 7e-12、不等式最大 3e-16、乗数最小 +2e-9。
  - 厳密評価の [0, 0.3] は、∫(r‴)² 0.0950 (現行 0.1059)、max|r⁗| 58.3 (現行 225)。
  - 参考: 現行 V0 は積分に敏感 (D₃ 0.87・D₄ 1.43、r″ の山 0.5102 → 0.5082)。既知の課題として記録し、本 plan では変えない。
