# ノズルの物理壁: スロートより上流に δ_r を足さず配管〜スロートを多項式 1 本にする / 寸法をスロート径から決める選択肢

## メタ

- **area**: `tooling`
- **status**: `draft`
- **related_docs**:
  - `methods/design/overview.md` (joint 壁の物理壁の解析経路、ランプ、`solve_rt`)
  - `procedures/nozzle-design-workflow.md` (寸法の決め方、`pw_ramp`)
- **related_plans**:
  - [tooling-nozzle-wall-single-bspline.md](tooling-nozzle-wall-single-bspline.md) (物理壁の表現。本 plan の上流多項式ではランプの当てはめが消え、全域が厳密になる)
  - [discretization-moc-axis-limit-and-corrector.md](discretization-moc-axis-limit-and-corrector.md) (設計壁の変更。本 plan は設計壁を変えない)
  - [tooling-nozzle-throat-monotone-r2.md](../accepted/tooling-nozzle-throat-monotone-r2.md) (現行の生産壁)
- **created**: `2026-10-07`
- **owner**: `Claude (Opus 5.5)`

## 1. 目的

今の物理壁は、スロートより上流の縮流部にも排除厚 δ_r を足している (x ∈ [−11, −6] のランプで 0 から全量へ、x ≥ −6 で全量)。これをやめ、上流は「配管〜スロートの 5 次多項式 1 本」、δ_r はスロートから下流にだけ足す形を選べるようにする。
あわせて、寸法 (設計の r_t) を出口径から決める今の方式に加えて、**物理スロート径から逆算する**方式を選べるようにする。

ユーザ決定 (2026-10-07):
- 「スロート点より上流に排除厚さを足しこむ処理を、やめにしませんか」「上流側は、配管〜スロート点まで、の多項式で形状を作る形」→ 上流の多項式化を標準にする (「a yes」)。
- スロート径を固定したい場合の方式は、下流で δ_r を 0 から徐変させる案ではなく、スロート径から r_t を逆算する案にする (「b 逆算にしようか」)。徐変案は、実効的な壁が設計壁から壁角で 0.12〜0.23° ずれる試算 (§3) のため採らない。

## 2. スコープ

- **やる**:
  - 物理壁の上流の作り方の選択肢 (`ramp` = 今 / `poly` = 本 plan) (§4.1)。
  - 寸法の決め方の選択肢 (出口径から / 物理スロート径から) の関数 (§4.2)。
  - 標準の手順 (`procedures/nozzle-design-workflow.md`) を `poly` にする。
- **やらない**:
  - 設計壁 (MOC・当てはめ・単調拘束) の変更。
  - δ_r の平滑化・積分法・CFD 抽出の変更。
  - 下流で δ_r を徐変させる方式 (§1 の理由で採らない)。
  - case/45 の生産壁の切り替えそのもの (NS・凝縮の再評価を経てユーザが決める。§5.1 #6)。

## 3. 関連 docs と前提

- **今の物理壁** (joint 壁の解析経路、`PhysicalNozzleWall.r`、`wall_axismach.py:497-500`): r_W = r_design + s(x)·δ_r(x)。r_design は直管 r_U ([x_in, −L_U))・上流 Hermite H ([−L_U, 0))・S ([0, x_e])。s はランプ `geometry.pw_ramp` の 5 次 smoothstep。
- **事前試算** (2026-10-07、CFD 0 step、run_0147 の物理壁。`case/45.isobutane_m6_d155/upstream_poly_probe.py` → `_band_ab/upstream_poly_probe.txt`):
  - 多項式 Q: 配管端 x = −L_U で (r_U, 0, 0)、設計スロート x = 0 で下流の物理壁 (S + δ_r) の (r, r′, r″) = (1.0013897, 7.76e-4, 0.501184) に 2 階微分まで連続。
  - 物理スロート (r′ = 0) の位置・半径は今と同じ (x = −0.1188 mm、r = 1.0013891 r_t)。
  - 今の物理壁との差: 縮流部で最大 −0.52 mm (x = −7.1 r_t。多項式の方が細い)。スロートから 1 r_t 以内で ≤ 6 µm、0.2 r_t 以内で ≤ 0.06 µm。今そこに足している δ_r は 0.48〜0.73 mm (x = −11〜−6)。
  - 上流 (x < x_t) で r′ ≤ 0 は成り立つ。r″ の今との差の最大 3.1e-3 /r_t (今のランプのゲート 5e-3 以内)。
  - 徐変案 (採らない): スロートの δ_r(0) = 0.107 mm で実効スロート面積 −0.28 %。徐変区間で実効的な壁 (物理壁 − 排除厚) の設計壁からのずれの傾きは、徐変長さ 1〜20 r_t で 0.12〜0.23°。
- **寸法の今の決め方**: `feedback/deltastar_loop.solve_rt(problem, R_exit_m, prev_run, euler_run)` が R = r_t·[r_F + δ_r(x_F)] を r_t について解く (CFD 前は積分法の δ_r、NS 後は抽出した δ_E を Re^−0.2 で換算)。結果は `spec.r_throat` に手で書く (`procedures/nozzle-design-workflow.md` ④・C2 方式)。
- **Euler の評価は設計壁を使う** (`runner_axismach.prepare` は `d["wall"]` でメッシュを作る)。物理壁の上流の変更は Euler の結果を変えない見込み。

## 4. 設計方針

### 4.1 上流の多項式 (`geometry.pw_upstream`)

- キー `geometry.pw_upstream`: `ramp` (今の作り方) | `poly` (本 plan)。joint 壁の物理壁の解析経路でだけ有効。他の経路で `poly` を指定したら例外。不正値は例外。
- **コードの既定は `ramp` のまま**にする (キーの無い過去の問題 YAML と `rerun_conditions` の再現を壊さないため)。**標準の手順** (`procedures/nozzle-design-workflow.md`・`procedures/recommended-settings.md` の該当行) は `poly` にする。ユーザの「既定に」はこの形で実現する (2026-10-07 の応答で説明)。
- `poly` の物理壁:
  - x ∈ [x_in, −L_U): r_U (直管)。
  - x ∈ [−L_U, 0): 5 次多項式 Q。端条件は x = −L_U で (r_U, 0, 0)、x = 0 で下流の物理壁 (S + δ_r) の x = 0⁺ の (r, r′, r″)。6 条件で一意。
  - x ∈ [0, x_e]: S + δ_r (今と同じ。この区間は今の物理壁とビット同一)。
  - `pw_ramp` は使わない (`poly` で `pw_ramp` も書かれていたら警告して無視するか例外にするかは諮問で決める)。
- **物理スロート**: r′ = 0 の根。Q の中、x = 0 の少し上流にある (case/45 で −0.0015 r_t)。今の囲い込みで探す (直管の r′ ≈ 0 を拾わない)。
- **ゲート** (今のランプのゲートを置き換える。不合格は例外で chain を止め、`prepare_info.json` に記録):
  - 上流の単調性: x ∈ [−L_U, x_t) で r′ ≤ 0 (x = −L_U の端は r′ = 0)。
  - 上流の曲率の変化: [−L_U, 0) で \|Q″ − H″\| ≤ 5e-3 (今のランプのゲートと同じ値)。
  - 継ぎ目 (−L_U・0) の値・1 階・2 階微分の跳び ≤ 1e-8。
- 診断として、今の物理壁 (`ramp`) との差の最大と位置を記録する。

### 4.2 寸法をスロート径から決める (`solve_rt_throat`)

- `feedback/deltastar_loop.solve_rt_throat(problem, R_throat_m, prev_run=None, euler_run=None)` を足す。`solve_rt` と同じ形で、物理スロートの半径 r_t·ρ_t(r_t) = R_throat を r_t について解く。ρ_t は物理壁 (`pw_upstream` に従う) の最小半径 (r_t 単位)。
  - CFD 前: 積分法の δ_r(x; r_t) で物理壁を作り直して不動点反復 (δ_r の r_t 依存は弱い)。
  - NS 後: 抽出した δ_E のスロートの値を Re^−0.2 で換算 (`solve_rt` の NS 後と同じ扱い)。
  - 戻り値に、反復の履歴、物理スロートの位置と半径、そのとき決まる出口半径 R_exit = r_t·[r_F + δ_r(x_F)] を入れる。
- 出口径から決める今の `solve_rt` は変えない。どちらで決めたかは、結果を `spec.r_throat` に書くときに問題 YAML のコメントと `prepare_info.json` (`sizing`) に残す。
- k_f の較正 (C2 方式、出口の δ を NS に合わせる) は寸法の決め方と独立なので、スロート径から決める場合にも使える。

## 5. 実装ステップ

1. `methods/design/overview.md`: joint 壁の物理壁の節に `pw_upstream`、寸法の節に `solve_rt_throat`。
2. `design/forge_design/geometry/wall_axismach.py`: `PhysicalNozzleWall` の解析経路に `upstream="poly"` (§4.1)。
3. `design/forge_design/evaluate/runner_axismach.py`: キーの読み取り・検査・`prepare_info.json` への記録。
4. `design/forge_design/feedback/deltastar_loop.py`: `solve_rt_throat` (§4.2)。
5. `design/tests/run_pw_upstream_poly_tests.py` (新規)。
6. `procedures/nozzle-design-workflow.md`・`procedures/recommended-settings.md`: 標準を `poly` に、スロート径から決める手順を追記。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の諮問 | codex diagnose。特に: コードの既定を `ramp` に残す判断、x = 0 で接続する判断 (物理スロートで接続する案との比較)、NS の再評価の方法 (§6 U4) | F |
| 2 | codex plan 段レビュー | `codex_review.py <本 plan> --stage plan` | O |
| 3 | 実装 | §5 の 1〜6。全域 1 本の B-spline の実装 (同じ `wall_axismach.py`・`runner_axismach.py` を触る) が終わってから。合格条件: 新テスト FAIL 0、design/tests の既存テスト FAIL 0 | O |
| 4 | 検証 (CFD 0 step) | §6 U0〜U3 | O (解釈は F) |
| 5 | codex result 段レビュー | `--stage result` (機能実装の完了でも必須) | O |
| 6 | case/45 の生産への反映 | §6 U4 (NS・凝縮)。MOC の plan の V5′ との順序・まとめ方は諮問で決める。採否はユーザ | ユーザ・O |

## 6. 検証

### 6.0 事前登録 (初稿 2026-10-07。諮問で確定する)

- **U0 既定のビット同一**: キー無しで、今の物理壁・`prepare_info.json` の既存の値・`solve_rt` の結果が変更前と完全一致。design/tests の既存テストが FAIL 0。
- **U1 上流の多項式の形** (case/45 の単調壁の生産問題、CFD 0 step):
  - [0, x_e] の物理壁が `ramp` とビット同一。
  - 物理スロートの差: 位置 ≤ 1e-6 r_t、半径 ≤ 1e-9 r_t (事前試算では 7 桁一致)。
  - §4.1 のゲートが PASS (上流の単調性、\|Q″ − H″\| ≤ 5e-3、継ぎ目の跳び ≤ 1e-8)。
  - 今の物理壁との差の最大と位置を記録 (事前試算 −0.52 mm、x = −7.1)。
  - 一般性: case/45 と違う L_U・r_U (標準の L_U = 3.5、r_U = 2.5) の問題でも、ゲートを通るかを記録する (通らなければ、その形状では `poly` を使えないことを記録)。
- **U2 スロート径から決める**:
  - 自己整合: 目標を、今の r_t (0.0766539 m) で同じ δ_r の経路 (CFD 前は積分法) から作った物理スロート半径にすると、r_t が今の値に 1e-9 m で戻る (`ramp` と `poly` の両方。物理スロートは両者で同じ)。
  - 目標を変えた場合 (例 ±5 %): 解いた r_t で作った物理壁の最小半径が目標に 1e-9 m で一致。反復の収束 (履歴) と、決まる出口半径を記録。
  - NS 後の経路: 抽出結果のある run (run_0147 の保存物) で動き、Re^−0.2 換算が `solve_rt` と同じ式であることをテストで確かめる。
- **U3 全域 1 本の B-spline との関係** (記録): `poly` の物理壁を全域 1 本の B-spline に作り直したときの誤差 (ランプが無いので全区間がノットの和集合で表せる見込み)。
- **U4 case/45 の NS・凝縮の再評価** (生産に入れる場合。方法は諮問で確定): 縮流部の形が最大 0.52 mm 変わるので、dry NS と凝縮 NS を回し直す。初期値の移し方 (前回の番号写像の上限 1 µm を超える) と、試験部の量 (出口 M・波・壁圧など) の比較の合否を、回す前に登録する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `design/forge_design/geometry/wall_axismach.py`・`evaluate/runner_axismach.py`・`feedback/deltastar_loop.py`。既定はビット同一。
- `methods/design/overview.md`・`procedures/nozzle-design-workflow.md`・`procedures/recommended-settings.md`。

## 8. 完了条件

- [ ] `methods/design/overview.md` を更新
- [ ] §6 U0〜U3 の結果を §9 に記録し、完了の区分を明記する: **機能実装の完了** (case/45 の生産は変えない) か **生産採用** (U4 とユーザ判断)
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done` にして accepted へ移動し、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-07` — 初稿。ユーザ決定「スロート点より上流に排除厚さを足しこむ処理を、やめにしませんか」「上流側は、配管〜スロート点まで、の多項式で形状を作る形」「a yes、b 逆算にしようか」。事前試算 (`case/45.isobutane_m6_d155/upstream_poly_probe.py`): 物理スロートは不変、縮流部で最大 −0.52 mm、徐変案は実効的な壁のずれ 0.12〜0.23° で採らない。
