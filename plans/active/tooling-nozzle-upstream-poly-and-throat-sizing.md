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
  - 物理スロート (r′ = 0) の位置・半径の今との差は、位置 −1.46e-9 r_t・半径 −3.77e-13 r_t (諮問の独立計算。x = −0.1188 mm、r = 1.0013891 r_t)。接続点は設計スロート x = 0 で、物理スロートではない。
  - 今の物理壁との差: 縮流部で最大 −0.52 mm (x = −7.1 r_t。多項式の方が細い)。スロートから 1 r_t 以内で ≤ 6 µm、0.2 r_t 以内で ≤ 0.06 µm。今そこに足している δ_r は 0.48〜0.73 mm (x = −11〜−6)。
  - 上流 (x < x_t) で r′ ≤ 0 は成り立つ。Q′ の実根は定義域内で x = −12 (重根、配管との接続) と x = −0.00155 (物理スロート) だけ。
  - \|Q″ − H″\| の最大は 1.184e-3 (x → 0⁻。= δ_r″(0) 相当)。初稿の「3.1e-3」は Q″ と今の物理壁 (H + s·δ_r) の差で、ゲートの量 (Q″ − H″) ではなかった (2026-10-07 諮問の指摘、主セッションで再計算して確認)。
  - 縮流部の幾何断面積は、x = −7.085 で今の物理壁より −0.270 % (諮問の独立計算)。
  - 徐変案 (採らない): スロートの δ_r(0) = 0.107 mm で実効スロート面積 −0.28 %。徐変区間で実効的な壁 (物理壁 − 排除厚) の設計壁からのずれの傾きは、徐変長さ 1〜20 r_t で 0.12〜0.23°。
- **寸法の今の決め方**: `feedback/deltastar_loop.solve_rt(problem, R_exit_m, prev_run, euler_run)` が R = r_t·[r_F + δ_r(x_F)] を r_t について解く (CFD 前は積分法の δ_r、NS 後は抽出した δ_E を Re^−0.2 で換算)。結果は `spec.r_throat` に手で書く (`procedures/nozzle-design-workflow.md` ④・C2 方式)。
- **Euler の評価は設計壁を使う** (`runner_axismach.prepare` は `d["wall"]` でメッシュを作る)。物理壁の上流の変更は Euler の結果を変えない見込み。

## 4. 設計方針

### 4.1 上流の多項式 (`geometry.pw_upstream`)

- キー `geometry.pw_upstream`: `ramp` (今の作り方) | `poly` (本 plan)。joint 壁の物理壁の解析経路でだけ有効。他の経路で `poly` を指定したら例外。不正値は例外。
- **コードの既定を `poly` にする** (ユーザ決定 2026-10-07「現在開発中なんだから過去の開発とか気にしなくていい」「既定はランプとしなくていい」。初稿の「過去の YAML の再現のために既定を `ramp` に残す」は取り下げ)。今の作り方は `pw_upstream: ramp` を明示したときだけ使う。
  - 既定が適用されるのは joint 壁の物理壁の解析経路だけ。他の経路ではキーが無ければ今の振る舞いのまま、`poly` を明示したら例外。
  - キーの無い過去の問題 YAML (case/45 の生産 YAML を含む) は、次に回すと `poly` の壁になる。どちらで作ったかは `prepare_info.json` に実効値として記録する。
  - 標準の手順 (`procedures/nozzle-design-workflow.md`・`procedures/recommended-settings.md`) も `poly` に揃え、`pw_ramp` を旧設定として記す。
- `poly` の物理壁:
  - x ∈ [x_in, −L_U): r_U (直管)。
  - x ∈ [−L_U, 0): 5 次多項式 Q。端条件は x = −L_U で (r_U, 0, 0)、x = 0 で下流の物理壁 (S + δ_r) の x = 0⁺ の (r, r′, r″)。6 条件で一意。
  - x ∈ [0, x_e]: S + δ_r (今と同じ。この区間は今の物理壁とビット同一)。
  - `pw_ramp` は使わない。`poly` (明示・既定とも) で `pw_ramp` も書かれていたら例外にする (諮問の推奨)。
  - 生成する問題 YAML・テンプレートには `pw_upstream` の値を明記し、解決済みの値を `prepare_info.json` に保存する (諮問の推奨)。
- **物理スロート**: r′ = 0 の根。Q の中、x = 0 の少し上流にある (case/45 で −0.0015 r_t)。Q′ (4 次多項式) の全実根と Q の極値から求め、定義域内で最小点が一意であることを確かめる (固定の囲い込みだけでは大域最小を保証しないため。諮問の指摘)。
- **ゲート** (今のランプのゲートを置き換える。不合格は例外で chain を止め、`prepare_info.json` に記録):
  - 正の半径、最小点が一意、その前 (x < x_t) で r′ ≤ 0・後で r′ ≥ 0 (Q′ の全実根で判定。x = −L_U の端は r′ = 0)。
  - 上流の曲率の変化: [−L_U, 0) で \|Q″ − H″\| ≤ 5e-3。これは今のランプのゲートと同じ値の**幾何の変化の上限**で、流れの品質を保証するものではない。
  - 継ぎ目 (−L_U・0) の値・1 階・2 階微分の跳び ≤ 1e-8 (左右の極限で比べる)。
- 診断として、今の物理壁 (`ramp`) との差の最大と位置を記録する。

### 4.1b 全域 1 本の B-spline は代数的に組む (最小二乗なし)

ユーザ指摘 (2026-10-07「ランプが無くなるので、シンプルに表現ができて、最小二乗フィットを省けますよね」) を採用する。`poly` の物理壁は全区間が x の 5 次の区分多項式なので、全域 1 本の B-spline (plan tooling-nozzle-wall-single-bspline) を当てはめでなく厳密に組む。

- 直管は定数の係数、Q は 5 次多項式を B-spline の係数へ直接変換 ([−L_U, 0] の Bézier)。
- [0, x_e] は S と δ_r (`delta_r_from_table` の補間スプライン) のノットを和集合にそろえ (ノット挿入。δ_r は x = 0 で重複度 6 まで挿入して区間を切り出す)、係数を足す。
- 継ぎ目 (−L_U・0) は重複度 3 (2 階微分まで連続)。区間ごとの表現を重複度 6 で並べたものから、連続性の分だけノットを除去する (除去の誤差を記録し、1e-12 r_t を超えたら例外)。
- ノット挿入は係数の凸結合だけなので、誤差は倍精度の丸め程度の見込み (U3 で測る)。全域 1 本の plan の最小二乗の版 (Gauss 点・特異値分解) は `ramp` のときだけに使う。生産壁を `poly` に切り替えた後で、`ramp` の版を残すかを決める (今の生産壁 run_0147 は `ramp` なので、それまでは STEP 化に要る)。

### 4.2 寸法をスロート径から決める (`solve_rt_throat`)

- `feedback/deltastar_loop.solve_rt_throat(problem, R_throat_m, ...)` を足す。方程式は r_t·min_x r_W(x; r_t) = R_throat (r_W は物理壁、r_t 単位)。
- **δ_r の供給経路を生産の壁と共通にする** (2026-10-07 諮問 Major): 反復のたびに、`prepare_ns` と同じ経路 (`integral_delta_r`: k_f の `cf_scale`・熱条件・5 次 P-spline の平滑化・`delta_r_from_table`) で δ_r を作り、同じ壁の構築 (`PhysicalNozzleWall`、`pw_upstream` に従う) で最小半径を求める。今の `solve_rt` の CFD 前の経路 (`integral_bl` を直接・未較正・未平滑化) は真似ない。run_0147 の保存物で、x = 0 の δ_r は直接経路 0.00101353 r_t、生産経路 0.00138968 r_t (実寸差約 28.8 µm、諮問の独立計算)。
- **NS 後の経路**: スロートの δ_E 1 点だけでは、接続の傾き・2 階微分・最小点・出口径が決まらない。NS 後に寸法を決めるときは、次の壁に実際に使う補正関数 (δ_E の全分布から作る δ_r 関数) を明示し、同じ関数で寸法と壁を作る。Re^−0.2 の換算は予測の近似として限定して使う。
- **k_f を較正し直したら寸法を解き直す** (k_f はスロートの δ_r も変えるため)。
- 戻り値に、反復の履歴、使った δ_r の出典と設定、物理スロートの位置と半径、そのとき決まる出口半径を入れる。どちらで寸法を決めたかを `prepare_info.json` (`sizing`) と問題 YAML のコメントに残す。
- スロート径を固定する試験では、出口径を同時に固定条件にしない。
- 出口径から決める今の `solve_rt` の CFD 前の経路にも同じ不一致がある (未較正・未平滑化)。同型の欠陥として、同じ共通経路に揃える (振る舞いが変わる。§7)。NS 後の経路 (C2 方式で使っている) は変えない。

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
| 1 | ~~§4・§6 の諮問~~ 完了 | 判断: 2026-10-07 codex (diagnose) — x = 0 接続は採用、ゲートは全実根で判定、`solve_rt_throat` は生産と共通の δ_r 経路で往復検証、U4 は `poly` を MOC 固定で先に評価し cross-mesh 移送。採否は §6.1 | F |
| 2 | codex plan 段レビュー | `codex_review.py <本 plan> --stage plan` | O |
| 3 | 実装 | §5 の 1〜6。全域 1 本の B-spline の実装 (同じ `wall_axismach.py`・`runner_axismach.py` を触る) が終わってから。合格条件: 新テスト FAIL 0、design/tests の既存テスト FAIL 0 | O |
| 4 | 検証 (CFD 0 step) | §6 U0〜U3 | O (解釈は F) |
| 5 | codex result 段レビュー | `--stage result` (機能実装の完了でも必須) | O |
| 6 | case/45 の生産への反映 | §6 U4 (NS・凝縮)。MOC の plan の V5′ との順序・まとめ方は諮問で決める。採否はユーザ | ユーザ・O |

## 6. 検証

### 6.0 事前登録 (初稿 2026-10-07。諮問で確定する)

- **U0 旧経路のビット同一**: `pw_upstream: ramp` を明示すると、今の物理壁・`prepare_info.json` の既存の値が変更前と完全一致。`solve_rt` (出口径から決める) の結果は、キーによらず変更前と一致する (下流の壁が同じため)。design/tests の既存テストは、ランプの振る舞いを試すものに `pw_upstream: ramp` を明示したうえで FAIL 0。キー無しが `poly` になることもテストで確かめる。
- **U1 上流の多項式の形** (case/45 の単調壁の生産問題、CFD 0 step):
  - [0, x_e] の物理壁が `ramp` とビット同一。
  - 物理スロートの差: 位置 ≤ 1e-6 r_t、半径 ≤ 1e-9 r_t (事前試算では 7 桁一致)。
  - §4.1 のゲートが PASS (上流の単調性、\|Q″ − H″\| ≤ 5e-3、継ぎ目の跳び ≤ 1e-8)。
  - 今の物理壁との差の最大と位置を記録 (事前試算 −0.52 mm、x = −7.1)。
  - 一般性: case/45 と違う L_U・r_U (標準の L_U = 3.5、r_U = 2.5) の問題でも、ゲートを通るかを記録する (通らなければ、その形状では `poly` を使えないことを記録)。
- **U2 スロート径から決める: 往復検証** (CFD 0 step、2026-10-07 諮問で改訂。逆算関数の内部の自己一致では判定しない):
  - 判別 A/B: 変えるのは逆算器の δ_r の供給経路だけ (A = 今の `solve_rt` 型の `integral_bl` 直接、B = `prepare_ns` と共通の経路)。MOC・`pw_upstream`・k_f・熱条件は固定。
  - 目標: 基準 (今の r_t で生産経路から作った物理スロート半径) と ±5 %。各目標で逆算した r_t を問題に入れ、**生産経路 (`prepare_ns` と同じ) で壁を作り直し**、連続な壁の最小半径と目標の差を測る。
  - 合格: B が全目標で \|差\| ≤ 1e-9 m (数値解法の許容差であって、float32 のメッシュや排除厚モデルの精度ではない)。A の結果も記録する。A だけ外れれば共通経路の採用を支持。両方通れば「経路の違いで外れる」仮説を棄却。B も外れれば反復の停止条件・寸法の反映・壁の生成の不一致を調べて止める。
  - 反復の履歴、使った δ_r の出典・設定、最小点の位置、決まる出口半径を記録する。
- **U3 全域 1 本の B-spline (代数的に組む、§4.1b)**: `poly` の物理壁と、ノット挿入で組んだ 1 本の B-spline の差 (半径・1〜2 階微分、各ノット区間の内部の密な評価点と継ぎ目の左右の極限)。合格: 半径 ≤ 1e-12 r_t、1 階 ≤ 1e-10、2 階 ≤ 1e-8 (丸めの範囲。全域 1 本の plan の W1 より厳しい)。最小二乗の版との差と、所要時間も記録する。
- **U4 case/45 の NS・凝縮の再評価** (生産に入れる場合、2026-10-07 諮問で方法を登録):
  - **順序**: 今の MOC (legacy・fixed2) に固定して `poly` だけを評価する。MOC の軸処理 (plan discretization-moc-axis-limit-and-corrector の V5′) はその後に別の変更として評価する (2 つの変更を同時に入れない)。
  - **初期値**: 縮流部の節点が最大約 0.5 mm 動くので、番号写像 (実績は約 0.5 µm の変更まで) は使わない。上限の数値だけを広げることもしない。cross-mesh の移送 (`interp_field.py`) とし、ドナーの対応・壁の層の対応・化学種と EOS の基底・新しい壁距離を検査する。移送直後の場を定常の比較に使わない。
  - **合否** (monotone plan §6 N・K の条件を転記): dry は出口コア M 5.9988〜6.0012、波 η0.1 ≤ 0.01 %、オーバーシュート ≤ 0.035 %、δ_E/δ_C = 1 ± 0.5 %、それぞれ `check_quasisteady` で STEADY。凝縮は §6 K と同じ。
  - **記録する量 (判定なし)**: 質量流量、音速線の位置、スロート付近の壁圧と排除厚の補正量の時系列。許容差を決める根拠が無いので、今回は値と run_0147・0148 との差を記録するだけにする (判定には使わない)。
  - Euler の評価は設計壁を使うので、上流の変更を検出できない。CFD でピン止めしている初期線 (run_0062 の Euler 由来) は固定した比較の基準として保持し、妥当性は物理壁の NS で確かめる。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-07` | [`notes/reviews/2026-10-07-upstream-poly-throat-sizing-diagnose.md`](../../notes/reviews/2026-10-07-upstream-poly-throat-sizing-diagnose.md) (ブリーフ [`briefs/2026-10-07-upstream-poly-throat-sizing.md`](../../notes/reviews/briefs/2026-10-07-upstream-poly-throat-sizing.md)) | M5/m2 | 採用: 上流の多項式化は採用、流れへの影響は Euler で検出できないので NS で確かめる・初期線は固定した基準 (M) / x = 0 接続、接続点は設計スロートと明記 (m) / ゲートの量の訂正 (初稿 3.1e-3 は Q″ − 今の壁、Q″ − H″ は 1.184e-3。主セッションで再計算)、Q′ の全実根で一意な最小点と単調性、5e-3 は幾何の上限 (M) / `solve_rt_throat` は生産と共通の δ_r 経路、U2 を往復検証の判別 A/B に、NS 後は δ_E の全分布から作る補正関数で寸法と壁を作る、k_f 較正後に寸法を解き直す (M×2) / U4 は MOC 固定で `poly` を先に、cross-mesh 移送、N・K の条件を転記、追加量は記録のみ (M) / `poly` と `pw_ramp` の併記は例外、解決済みの値を保存 (m)。不要: 「コードの既定を `ramp` に残す」への回答 (ユーザ決定で既定は `poly`)。追加: 今の `solve_rt` の CFD 前の経路も同型の不一致なので共通経路に揃える |

## 7. 影響範囲

- `design/forge_design/geometry/wall_axismach.py`・`evaluate/runner_axismach.py`・`feedback/deltastar_loop.py`。`solve_rt` の CFD 前の経路の結果が変わる (共通の δ_r 経路に揃えるため)。joint 壁の物理壁の既定が `poly` に変わる (キーの無い問題 YAML の物理壁が縮流部で変わる)。`pw_upstream: ramp` で今の壁を再現できる。
- `methods/design/overview.md`・`procedures/nozzle-design-workflow.md`・`procedures/recommended-settings.md`。

## 8. 完了条件

- [ ] `methods/design/overview.md` を更新
- [ ] §6 U0〜U3 の結果を §9 に記録し、完了の区分を明記する: **機能実装の完了** (case/45 の生産は変えない) か **生産採用** (U4 とユーザ判断)
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done` にして accepted へ移動し、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-07` — 初稿。ユーザ決定「スロート点より上流に排除厚さを足しこむ処理を、やめにしませんか」「上流側は、配管〜スロート点まで、の多項式で形状を作る形」「a yes、b 逆算にしようか」。事前試算 (`case/45.isobutane_m6_d155/upstream_poly_probe.py`): 物理スロートは不変、縮流部で最大 −0.52 mm、徐変案は実効的な壁のずれ 0.12〜0.23° で採らない。
- `2026-10-07` — **ユーザ決定「現在開発中なんだから過去の開発とか気にしなくていいんだが。既定はランプとしなくていい」**: コードの既定を `poly` に変更 (§4.1・§6 U0・§7)。旧経路は `pw_upstream: ramp` の明示で再現する。
- `2026-10-07` — **codex (diagnose) に諮った**: `notes/reviews/2026-10-07-upstream-poly-throat-sizing-diagnose.md` — 上流の多項式化は採用。`solve_rt_throat` は生産と共通の δ_r 経路で、U2 は往復検証の判別 A/B。U4 は MOC を固定して `poly` を先に評価する。全件採用 (既定に関する問いはユーザ決定で不要) し §3・§4.1・§4.2・§5.1・§6 U2・U4・§7 を改訂。
- `2026-10-07` — **ユーザ指摘「ランプが無くなるので、シンプルに表現ができて、最小二乗フィットを省けますよね」**: 採用。`poly` では全域 1 本をノット挿入で代数的に組む (§4.1b)。U3 を合否つきに改訂。最小二乗の版は `ramp` (今の生産壁) の STEP 化のために当面残す。
