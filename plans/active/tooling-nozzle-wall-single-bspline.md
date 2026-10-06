# ノズルの物理壁を全域 1 本の x の 5 次 B-spline で表し、STEP で渡す (設計の中身は変えない)

## メタ

- **area**: `tooling`
- **status**: `draft`
- **related_docs**:
  - `methods/design/overview.md` (joint 壁、物理壁の解析経路)
- **related_plans**:
  - [tooling-nozzle-throat-monotone-r2.md](../accepted/tooling-nozzle-throat-monotone-r2.md) (現行の生産壁)
  - [discretization-moc-axis-limit-and-corrector.md](discretization-moc-axis-limit-and-corrector.md) (並行する設計壁の変更。本 plan は設計を変えないので独立)
- **created**: `2026-10-07`
- **owner**: `Claude (Opus 5.5)`

## 1. 目的

今の物理壁は、区間ごとに別の式の和になっている。入口直管の定数、収縮部の 5 次 Hermite H(x)、拡大部の 5 次 B-spline S(x) (係数 259 個)、境界層の排除厚 δ_r(x) (ξ = √(x + 16.5) の P-spline を 1500 点で補間し直したもの)。
これを、入口から出口まで **1 本の x の 5 次 B-spline** (ノット列と係数) に作り直し、メッシュ生成・報告・CAD のすべてで同じものを使う。CAD には **STEP** (平面の B-spline 曲線) で渡す。

ユーザ決定 (2026-10-07): 「できれば共通化させたい」→ 当初の「δ_r の平滑化を x の B-spline に替える」案 (B) は、事前に決めた条件を満たす λ が見つからず (plan 段レビュー)、物理壁が数 µm 動いて NS の再評価が要るため、**「設計の中身は変えず表現だけ 1 本にする」案 (A)** に切り替えた (「A でいいかな」)。CAD 形式は「STEP でいいよ」。

## 2. スコープ

- **やる**:
  - 今の物理壁を、全域 1 本の 5 次 B-spline に、許容誤差つきで作り直す (§4.1)。設計の中身 (MOC・当てはめ・δ_r の平滑化・ランプ) は変えない。
  - その 1 本を壁クラスとして、NS の準備・初期値・メッシュ生成・報告が使う (§4.2)。
  - STEP の書き出しと読み直しの検査 (§4.3)。
  - 問題 YAML のキーで選ぶ。既定は今の壁とビット同一。
- **やらない**:
  - δ_r の平滑化の変更 (案 B は取り下げ。§9 に経緯)。
  - 外側の形 (肉厚・フランジ) の CAD 化。
  - IGES・DXF の書き出し。

## 3. 関連 docs と前提

- **物理壁の今の式** (r_t 単位): [−12.5, −12) は 6.485。[−12, −11) は H、[−11, −6) は H + s·δ_r (s は 5 次 smoothstep)、[−6, 0) は H + δ_r、[0, x_e] は S + δ_r。継ぎ目 (−12・−11・−6・0) で 2 階微分まで連続。
- δ_r は x の区分多項式ではない (ξ = √(x − x₀ + x_stretch) の多項式)。そのため、拡大部も含めて全域で**当てはめ**が要る (厳密な変換ではない)。
- **メッシュの座標は float32** (`geom_float`、`flowFormat.hpp:7`)。r ≈ 0.775 m で 1 ulp ≈ 6e-8 m (0.06 µm)。壁の差を 0.01 µm (1.3e-7 r_t) 以下にすれば、変換後のメッシュ座標は丸めの範囲でしか変わらない見込み (§6 で確かめる)。
- **STEP の書き出し**: FreeCAD 1.1.1 (OpenCascade) が `/home/sano/opt/squashfs-root/usr/bin/freecadcmd` にあり、5 次の B-spline 曲線を STEP (`B_SPLINE_CURVE_WITH_KNOTS`) に書き出せることを確かめた (2026-10-07)。

## 4. 設計方針

### 4.1 1 本の B-spline への作り直し

- 定義域 [x_in, x_e] (= [−12.5, x_e])。次数 5。両端はノットの重複度 6。
- **継ぎ目のノット**: x = −12・−11・−6・0 は、元の連続性 (2 階微分まで) に合わせて重複度 3 のノットにする。こうしておけば、各区間の中だけで近似すればよく、継ぎ目の値・1〜2 階微分は両側から拘束で合わせられる。
- **区間内のノット**: 初期値は、拡大部では設計壁 S のノットと、δ_r の P-spline のノット (ξ のノットを x に戻した位置) の和集合。収縮部では δ_r のノット。
- **当てはめ**: 各区間で、継ぎ目の値・1〜2 階微分を拘束した最小二乗。評価点は区間ごとに十分密にする。
- **許容誤差** (元の物理壁との差、区間多項式の極値と密な評価点の両方で検査):
  - 半径 \|Δr\| ≤ 1.3e-7 r_t (0.01 µm)。
  - 1 階微分 \|Δr′\| ≤ 1e-7、2 階微分 \|Δr″\| ≤ 1e-5 (半径誤差の振幅 ε・幅 h の振動は r″ に ε/h² で出るため、微分にも別に上限を置く)。
  - 超えた区間は、そこのノットを二分して当てはめ直す。上限の反復で満たさなければ不合格として止める。
- **形のゲートの維持**: 単調壁の形状ゲートのうち、物理壁にかかるもの (S6: 物理壁の r″ の最大増加 ≤ 0.002) と、ランプのゲート (`max|r″ − r″_design| ≤ 0.005`、`max r′ < 0`、`wall_axismach.py:399`) を、新しい 1 本でも満たすこと。
- 係数の数・ノットの数・各区間の最大誤差を記録する。

### 4.2 壁クラスと下流の道具

- 新しい壁クラスは、今の `PhysicalNozzleWall` が下流に出している属性と診断をすべて持つ (plan 段レビューの指摘): `r(x, deriv)`・`theta`・`x_in`・`x_e`・`validate()`・スロート量 (`x_throat`・`r_throat`・`kappa_throat`)・`offset_mode`・`_dstar_hist` 相当など。属性の一覧は実装前に `runner_axismach.py`・`evaluate/ic.py`・`report/nozzle_report.py` の読み出し箇所から洗い出して表にする。
  - スロート量は、1 本の B-spline から求め直す (r′ = 0 の点と、そこでの r″)。元の壁の値との差を記録する。
  - `evaluate/ic.py` は、スロート属性が無いと黙って (0, 1) に戻る (`ic.py:68`)。新しいクラスでは属性を必ず持たせ、欠けたら例外にする。
- 報告 (`nozzle_report.py`) は、CSV の再補間ではなく、run に保存した係数・ノットから壁を復元して評価する (共通の読み込み関数を作る)。
- 対応する組み合わせ: joint 壁 + 物理壁の解析経路 (`offset: radial`)。それ以外は、新しいキーを指定したら例外。
- ノット・係数・次数・定義域・単位 (r_t と m) を run の `prepare_info.json` と JSON ファイルに保存する。

### 4.3 STEP の書き出し

- 平面の B-spline 曲線 C(u) = (x(u), r(u), 0)。助変数 u = x、次数 5、同じノット列、重み 1。
- 制御点は Pᵢ = (x̄ᵢ, cᵢ)。x̄ᵢ はグレビル点 (ノット 5 個の平均)。B-spline は 1 次関数を正確に表すので、x(u) = u が厳密に再現され、曲線は r(x) と一致する。
- 単位は mm (r_t = 76.6539 mm を掛ける)。x 軸を流れ方向 (軸)、原点はスロート。曲線は上半分 (r ≥ 0)。
- 書き出しは `freecadcmd` で `Part.BSplineCurve.buildFromPolesMultsKnots` を使う。STEP には曲線 1 本だけを入れ、入口端・出口端の座標を添え書き (JSON) に記す。
- 回転体 (内面) にするのは CAD 側の作業とする。

### 4.4 既定値と切り替え

- 問題 YAML の `geometry.physical_wall_repr` (`legacy` | `single_bspline`、既定 `legacy`)。既定は今の壁とビット同一。
- 生産に入れるか (case/45 の問題 YAML にキーを入れるか) は、§6 の結果を見てユーザが決める。

## 5. 実装ステップ

1. `methods/design/overview.md`: 物理壁の表現と STEP の節。
2. `design/forge_design/geometry/wall_axismach.py`: 1 本の B-spline への作り直し (§4.1) と壁クラス (§4.2)。
3. `design/forge_design/evaluate/runner_axismach.py`・`evaluate/ic.py`・`report/nozzle_report.py`: キーの読み取り、保存、共通の読み込み関数。
4. `design/forge_design/export/wall_step.py` (新規): STEP の書き出しと読み直しの検査 (`freecadcmd` を呼ぶ)。
5. `design/tests/run_wall_single_bspline_tests.py` (新規)。
6. §6 の W0〜W5。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の諮問 (案 A への書き直し後) | codex diagnose | F |
| 2 | codex plan 段レビュー | `codex_review.py <本 plan> --stage plan` | O |
| 3 | 実装 | §5 の 1〜5。MOC の軸処理の実装が終わってから (同じ `runner_axismach.py` を触る)。合格条件: 新テスト FAIL 0、design/tests の既存テスト FAIL 0 | O |
| 4 | 検証 (CFD 0 step) | §6 W0〜W5 | O (解釈は F) |
| 5 | codex result 段レビュー | `--stage result` (機能実装の完了でも必須) | O |
| 6 | 生産への反映 | case/45 の問題 YAML にキーを入れるか (ユーザ判断)。W3 で メッシュの差が丸めの範囲なら CFD のやり直しは要らない | ユーザ |

## 6. 検証

### 6.0 事前登録 (案 A、2026-10-07)

- **W0 既定のビット同一**: キー無しで、今の物理壁 (run_0147 の `wall_physical.csv`・メッシュの座標) と完全一致。design/tests の既存テストが FAIL 0。
- **W1 作り直しの精度** (単調壁の生産問題、CFD 0 step): §4.1 の許容誤差 (半径 ≤ 1.3e-7 r_t、r′ ≤ 1e-7、r″ ≤ 1e-5、区間多項式の極値と密な評価点の両方)。継ぎ目の 1〜2 階微分の跳び ≤ 1e-8。区間の中は 4 階微分まで連続。形のゲート (S6・ランプのゲート) が PASS。係数・ノットの数を記録。
- **W2 スロート量**: 1 本の B-spline から求めたスロートの位置・半径・曲率と、元の壁の値の差を記録する。合格: 位置 ≤ 1e-6 r_t、半径 ≤ 1.3e-7 r_t、曲率 ≤ 1e-5。
- **W3 メッシュへの影響**: 同じメッシュ設定で、元の壁と新しい 1 本からメッシュを作り、変換後の座標 (float32) を比べる。
  - 合格: 全節点で座標の差 ≤ 2 ulp (float32)、接続・境界対応が同一、メッシュ品質 PASS。
  - 合格なら、物理壁を新しい 1 本に替えても CFD の結果は変わらない (メッシュが丸めの範囲で同じ) とみなし、NS・凝縮のやり直しは要らない。不合格なら止めて諮問する。
- **W4 下流の道具**: `prepare_ns`・初期値 (`ic.py`)・報告 (`nozzle_report.py`) が新しい壁を使って動く。報告が CSV の再補間ではなく保存した係数から評価していることを確かめる。§4.2 の属性の表のすべてが埋まっていること。
- **W5 STEP**:
  - 書き出した STEP を OpenCascade (`freecadcmd`) で読み直し、曲線の次数・ノット・重複度・制御点数が書き出し時と一致すること。
  - 曲線上の点・接線・曲率を、独立の評価器 (scipy の `BSpline` を使わない de Boor の自前実装) と、元の物理壁の関数の両方と比べる。合格: 位置 ≤ 1e-6 mm。接線方向の角度 ≤ 1e-9 rad。曲率は絶対 ≤ 1e-9 /mm または相対 ≤ 1e-9 (零曲率の直管があるため、絶対と相対の両方の基準を置く)。重複ノットの位置では左右の極限で比べる。
  - mm への換算で、n 階微分が r_t^(1−n) 倍になることを確かめる。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-07` | [`notes/reviews/2026-10-07-wall-single-bspline-diagnose.md`](../../notes/reviews/2026-10-07-wall-single-bspline-diagnose.md) (ブリーフ [`briefs/2026-10-07-wall-single-bspline.md`](../../notes/reviews/briefs/2026-10-07-wall-single-bspline.md)) | C0/M6/m1 | 全件採用: ① λ を EDF 同等で決めるのは却下 (残差 2.66 % で不合格、諮問の独立試算)。残差と局所形状の条件を満たす最大の λ に改訂、判別 A/B (λ = 0 対 EDF 同等) を記録 / ② r″ の高周波を新旧の壁関数で直接、区間別に評価 / ③ ランプに Δr′・Δr″ の上限と今のゲートを追加、positive の発火は不合格 / ④ 出口半径の変化だけで較正の要否を決めない。較正器に同じ平滑化を通し、W4 で出口半径と δ_E/δ_C を別々に判定 / ⑤ δ_r を使うすべての入口 (delta_r_csv・δ* ループ・較正器) に同じ平滑化 / ⑥ NS・凝縮の統合評価は、MOC と本件の形状検証を先に終えてから / ⑦ CAD の仕様と独立評価器での検査 |
| plan | `2026-10-07` | [`notes/reviews/2026-10-07-tooling-nozzle-wall-single-bspline-plan.md`](../../notes/reviews/2026-10-07-tooling-nozzle-wall-single-bspline-plan.md) (案 B に対して) | GO-with-changes, C0/M4/m2 | 案 B を取り下げ、案 A に書き直した (ユーザ決定)。M1 (W1 を満たす λ が 242 候補中 0) が取り下げの理由。M3 (壁クラスの属性・ic.py の黙った既定・報告の CSV 再補間) と m5 (量ごとの許容差・零曲率) は案 A の §4.2・§6 W4・W5 に反映。M2 (CSV の意味・二重平滑化) は平滑化を変えないので該当しない。M4 (IC の移送) は W3 でメッシュの差が丸めの範囲なら CFD を回さないので該当しない (不合格なら諮問)。m6 は §5.1 #5 に反映 |

## 7. 影響範囲

- `design/forge_design/geometry/wall_axismach.py`・`evaluate/runner_axismach.py`・`evaluate/ic.py`・`report/nozzle_report.py`・`export/wall_step.py` (新規)。既定はビット同一。
- `methods/design/overview.md`。

## 8. 完了条件

- [ ] `methods/design/overview.md` を更新
- [ ] §6 W0〜W5 の結果を §9 に記録し、完了の区分を明記する: **機能実装の完了** (生産の YAML は変えない) か **生産採用** (case/45 の YAML にキーを入れる、ユーザ判断)
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done` にして accepted へ移動し、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-07` — 初稿。ユーザ決定「設計側の平滑化そのものを x の B-spline に替える」「できれば共通化させたい」「(i) で、全域 1 本にして進めて」。δ_r の第 1 段を直接使う案 (値の差 6e-11 r_t) は、x の区分多項式にならず CAD に渡せないので、x 空間の平滑化に置き換える方針にした。
- `2026-10-07` — **codex (diagnose) に諮った**: `notes/reviews/2026-10-07-wall-single-bspline-diagnose.md` — 全域 1 本の方針は維持し、λ の決め方を「残差と局所形状の条件を満たす λ」に差し戻す。諮問の独立試算 (run_0147 の保存入力): EDF 同等の λ = 17.58 で残差 2.66 % (今 1.90 %)、物理壁の変化は最大約 10 µm (x ≈ −5.4)、出口 +4.3 µm。同じ基底で λ = 0 なら残差 0.68 %。全件採用し §4・§6 を改訂。
- `2026-10-07` — **plan 段レビュー (案 B)**: `notes/reviews/2026-10-07-tooling-nozzle-wall-single-bspline-plan.md` (GO-with-changes、C0/M4/m2)。codex の独立試算で、δ_r を x の B-spline で平滑化すると、残差の条件を満たす λ (173 候補) のどれもランプの r″ の高周波が今を超えた (最良 1.001 倍、x = −11)。物理壁は出口で約 6 µm 動く。事前の取り決めどおり止めてユーザに相談。
- `2026-10-07` — **ユーザ決定「A でいいかな」「STEP でいいよ、それで進めて」**: 設計の中身は変えず、今の物理壁を全域 1 本の 5 次 B-spline に許容誤差つきで作り直す案 A に plan を書き直した。CAD には STEP (平面の B-spline 曲線、グレビル点の制御点、mm) で渡す。FreeCAD 1.1.1 で 5 次の B-spline 曲線の STEP 書き出しを確かめた。
