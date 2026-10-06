# ノズルの物理壁を全域 1 本の x の 5 次 B-spline で表す (δ_r の平滑化を x 空間に)

## メタ

- **area**: `tooling`
- **status**: `draft`
- **related_docs**:
  - `methods/design/overview.md` (joint 壁、物理壁の解析経路、δ_r の平滑化)
- **related_plans**:
  - [tooling-nozzle-throat-monotone-r2.md](../accepted/tooling-nozzle-throat-monotone-r2.md) (現行の生産壁。単調壁の問題と NS・凝縮のゲート N・K)
  - [discretization-moc-axis-limit-and-corrector.md](discretization-moc-axis-limit-and-corrector.md) (並行する設計壁の変更。生産に入れる場合の NS・凝縮の再評価をまとめる)
- **created**: `2026-10-07`
- **owner**: `Claude (Opus 5.5)`

## 1. 目的

今の物理壁は、区間ごとに別の式の和になっている。

- 入口直管の定数。
- 収縮部の 5 次 Hermite H(x)。
- 拡大部の 5 次 B-spline S(x) (係数 259 個)。
- 境界層の排除厚 δ_r(x)。δ_r は 2 段で作っている。
  - 第 1 段: ξ = √(x + 16.5) の 5 次 P-spline、係数 71 個。
  - 第 2 段: その 1500 点の値を通す 5 次補間スプライン (約 1500 区間)。

ユーザ決定 (2026-10-07):「設計側の平滑化そのものを x の B-spline に替える」「できれば共通化させたい」「(i) で、全域 1 本にして進めて」。

物理壁を、入口から出口まで 1 本の x の 5 次 B-spline (ノット列と係数) で表す。式を簡潔にし、CAD にそのまま渡せる形にする。

## 2. スコープ

- **やる**:
  - δ_r の平滑化を、x の 5 次 B-spline (不等間隔ノット、∫(δ_r‴)² dx の罰則) に替える。
  - 物理壁を全域 1 本の 5 次 B-spline に組み立てる。式がそのまま 5 次多項式になる所は厳密に、ならない所 (ランプ区間) は許容誤差つきの当てはめで。
  - メッシュ生成・報告・NS の準備がこの 1 本を使う。ノット・係数を書き出す (CAD 用)。
  - 生産に入れるかどうかの判断。物理壁が µm 級で変わるので、NS・凝縮の再評価が要る。MOC の軸処理の plan も生産に入れる場合は、再評価をまとめて 1 回にする。
- **やらない**:
  - 設計壁 (Euler が使う非粘性の壁) の変更。こちらは MOC の plan が扱う。
  - δ_r の物理モデル (CONTUR の積分法・k_f) の変更。
  - STEP/IGES などの CAD 形式への書き出し。ノット・係数の CSV/JSON までにとどめる。

## 3. 関連 docs と前提

- **今の δ_r の平滑化** (`design/forge_design/metrics/deltastar.py::smooth_delta_quintic`): 横軸 ξ = √(x − x₀ + 4) で等間隔のノット (x₀ = −12.5)、係数 71 個 (中央でノット間隔 ≈ 2 r_t)、係数の 3 階差分の罰則 λ = 1。当てはめの残差は相対 rms 1.9 %、最大 3.5e-4 r_t (run_0147)。
- **第 2 段の補間** (`runner_axismach.delta_r_from_table`): 1500 点 (0.0719 r_t 間隔) を通す 5 次補間スプライン。第 1 段との差は、値で 6e-11 r_t、2 階微分で 3e-7 (2026-10-07 の試算)。
- **物理壁の今の式** (r_t 単位):
  - [−12.5, −12) は 6.485。
  - [−12, −11) は H。[−11, −6) は H + s·δ_r。[−6, 0) は H + δ_r。ランプは s = 10ξ³ − 15ξ⁴ + 6ξ⁵、ξ = (x + 11)/5。
  - [0, x_e] は S + δ_r。
  - 継ぎ目 (−12・−11・−6・0) で 2 階微分まで連続。
- ランプ区間の H + s·δ_r は、区間ごとに 10 次の多項式になるので、5 次 B-spline では厳密に書けない。

## 4. 設計方針

### 4.1 δ_r の平滑化 (x 空間)

- 5 次 B-spline in x。内部ノットは、今の ξ の等間隔ノットを x に戻した位置 (x_k = ξ_k² − 16.5) に置く。今と同じ係数の数で、スロート側で密になる。
- 罰則は λ ∫(δ_r‴)² dx。区間ごとの Gauss 求積で厳密に評価する。不等間隔ノットでは、係数の差分の罰則は意味がずれるため。
- **λ の決め方 (結果を見る前に登録)**: 今の P-spline と有効自由度 (平滑化行列のトレース) が等しくなる λ を選ぶ。等しくできなければ、最も近い値を選んで記録する。
- 断熱壁の「δ_r ≥ 0」の切り落とし (`positive`) は、今と同じく値にだけ掛け、発火したら記録する (発火すると滑らかさが壊れるため。今の生産では最小 0.0011 r_t で発火していない)。

### 4.2 全域 1 本の B-spline の組み立て

- 区間ごとの式 (§3) から、1 本の 5 次 B-spline を組み立てる。
  - 直管 (定数)、[−12, −11) の H、[−6, 0) の H + δ_r、[0, x_e] の S + δ_r は、5 次の区分多項式なので、ノットの和集合の上で厳密に書く (変換誤差は丸め程度)。
  - 継ぎ目では、元の連続性に合わせてノットを重ねる (2 階微分まで連続の継ぎ目は重複度 3)。
  - ランプ [−11, −6) は H + s·δ_r が 10 次なので、区間を細かいノットで分けて最小二乗で当てはめる。値・1〜2 階微分を両端で合わせる拘束つきで、許容誤差は半径で 0.1 µm (1.3e-6 r_t) 以下。満たすまでノットを細かくする。
- 結果は `PhysicalNozzleWall` と同じ呼び出し (`r(x, deriv)`・`x_in`・`x_e`) で使えるクラスにし、メッシュ生成・報告・NS 準備はそれを使う。
- ノット列・係数・次数を run の `prepare_info.json` と、別ファイル (CSV/JSON、m 単位と r_t 単位) に書き出す。

### 4.3 既定値と切り替え

- 問題 YAML の `geometry.physical_wall_repr` (`legacy` | `single_bspline`、既定 `legacy`)。既定は今の壁とビット同一。
- 生産に入れるかは §6 の結果を見てユーザが決める。

## 5. 実装ステップ

1. `methods/design/overview.md`: δ_r の平滑化と、物理壁の表現の節を更新する。
2. `design/forge_design/metrics/deltastar.py`: x 空間の P-spline (§4.1)。今の関数は残す。
3. `design/forge_design/geometry/wall_axismach.py`: 全域 1 本の B-spline の組み立てと壁クラス (§4.2)。
4. `design/forge_design/evaluate/runner_axismach.py`: キーの読み取り、`prepare_ns` の経路 (積分法初期化と `delta_r_csv` の両方)、書き出し。
5. `design/tests/run_wall_single_bspline_tests.py` (新規)。
6. §6 の W1〜W5。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の諮問 | codex diagnose | F |
| 2 | codex plan 段レビュー | `codex_review.py <本 plan> --stage plan` | O |
| 3 | 実装 | §5 の 1〜5。MOC の軸処理の実装が終わってから (同じ `runner_axismach.py` を触る)。合格条件: 新テスト FAIL 0、design/tests の既存テスト FAIL 0 | O |
| 4 | 検証 (CFD 0 step) | §6 W1〜W3・W5 | O (解釈は F) |
| 5 | 生産への反映 (候補にする場合) | W4 (NS・凝縮の再評価。MOC の plan と合わせて 1 回) → codex result 段レビュー → ユーザ判断 | ユーザ・O |

## 6. 検証

### 6.0 事前登録 (初稿 2026-10-07、諮問の前)

- **W0 既定のビット同一**: キー無しで、今の物理壁 (run_0147 の `wall_physical.csv`・メッシュの座標) と完全一致。design/tests の既存テストが FAIL 0。
- **W1 δ_r の新しい平滑化** (単調壁の生産問題、CFD 0 step): 今の δ_r との差 (値・1〜3 階微分の最大と位置)、当てはめの残差 (相対 rms・最大)、有効自由度、`positive` の発火を記録する。
  - 合格: 残差の相対 rms ≤ 今 (1.9 %) × 1.1。物理壁の r″ の高周波 (報告の `r2_highfreq_max_x_gt2`) ≤ 今の値。単調壁の形状ゲート S6 (物理壁の r″ の最大増加 ≤ 0.002) が PASS。出口半径の変化 ≤ 0.01 mm (r_t を解き直さなくてよい目安)。
- **W2 全域 1 本の表現**: 区間ごとの式との差。厳密に書く区間は ≤ 1e-10 r_t、ランプ区間は ≤ 1.3e-6 r_t (0.1 µm)。継ぎ目の 1〜2 階微分の跳び ≤ 1e-8、それ以外の区間内は 4 階微分まで連続。ノット数と係数の数を記録する。
- **W3 生産の入口の再現**: monotone plan の `throat_mono_entry_prep_ab.py` と同じ方法で、`prep_c2pin.py` と `deltastar_loop --init-integral` の両方が同じ 1 本の B-spline と同じメッシュを作ることを確かめる。
- **W4 NS・凝縮の再評価** (生産に入れる場合): monotone plan の §6 N・K と同じ手順・ゲート。MOC の plan も入れる場合は、両方の変更を入れた壁で 1 回にする。
- **W5 下流の道具**: メッシュ品質 (AR ≤ 5000・skew ≤ 0.9) が PASS。報告 (nozzle_report) が新しい壁で動く。書き出したノット・係数から壁を復元すると、元と ≤ 1e-12 r_t で一致する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `design/forge_design/metrics/deltastar.py`・`geometry/wall_axismach.py`・`evaluate/runner_axismach.py` (オプション追加、既定はビット同一)。
- 報告・メッシュ生成は壁クラス経由で使うので、呼び出しは変わらない見込み。
- `methods/design/overview.md`。

## 8. 完了条件

- [ ] `methods/design/overview.md` を更新
- [ ] 完了の区分を §9 に明記する: **機能実装の完了** (W0〜W3・W5) か **生産採用** (加えて W4 と result レビュー、ユーザ判断)
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done` にして accepted へ移動し、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-07` — 初稿。ユーザ決定「設計側の平滑化そのものを x の B-spline に替える」「できれば共通化させたい」「(i) で、全域 1 本にして進めて」。δ_r の第 1 段を直接使う案 (値の差 6e-11 r_t) は、x の区分多項式にならず CAD に渡せないので、x 空間の平滑化に置き換える方針にした。
