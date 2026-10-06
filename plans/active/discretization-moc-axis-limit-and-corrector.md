# MOC の軸上ソース項を解析極限 θ_r にし、予測修正を収束まで回す

## メタ

- **area**: `discretization`
- **status**: `draft`
- **related_docs**:
  - `methods/design/overview.md` (逆 MOC の軸対称源項 A11、「今後の課題」の軸上の解析極限)
- **related_plans**:
  - [discretization-moc-axisymmetric-source-term.md](../accepted/discretization-moc-axisymmetric-source-term.md) (A11: 源項を点自身で評価。軸上だけ相手から極限を代用する現行の形を決めた plan)
  - [tooling-nozzle-throat-monotone-r2.md](../accepted/tooling-nozzle-throat-monotone-r2.md) (§3 仮説 H、§5.1 #10、§9 2026-10-06〜07: 始点付近の角度差の切り分けと事前試算)
- **created**: `2026-10-07`
- **owner**: `Claude (Opus 5.5)`

## 1. 目的

逆 MOC の単位過程には、軸の 1 段目に 2 つの誤差がある (monotone plan の試験 (b))。

- 軸上の端点では sinθ/r が 0/0 なので、相手の点の値を借りて代用している。借りる値が組み合わせる相手によって変わる。
- 予測修正の反復回数が 2 回固定 (`n_corr=2`) で、軸の近くでは収束しない (振動の比 ≈ −1/2)。

これを次の 2 点で直し、MOC の始点付近の角度差 (スロート直後の r″ の山の主因の一部) を減らす。

- 軸上の端点のソース項に、軸上の M 分布から求めた解析極限 θ_r を使う。
- 予測修正を収束まで回す。

ユーザ決定 (2026-10-07)。「試験用に直した手順 (予測修正を 20 回まで回す)、総点数 2400 をメインに」→ 軸上のソース項は 0 ではなく θ_r にする案で合意。θ_r は semi-perfect で求める。

## 2. スコープ

- **やる**:
  - `design/forge_design/geometry/moc_kernel.py` (と `moc_inverse.py` の呼び出し) に、軸上の極限 θ_r と、収束判定つきの予測修正を入れる。
  - 問題 YAML のキーで選ぶ (`geometry.moc_axis_limit: analytic`・`geometry.moc_corrector: converge`)。既定は現行どおりで、ビット同一。
  - 厳密解 (放射源流) で精度を確かめる。再出発の自己整合試験、case/45 の形状の試算、所要時間を測る。
  - case/45 の生産に入れるかどうかは §6 の結果で決める。入れる場合は Euler で出口 M を確かめる。
- **やらない**:
  - 他の case の問題 YAML の変更。
  - `AXIS_LIMIT_FRAC` (軸に近い軸外の点の扱い) の変更。
  - 初期線のデータ側の不整合の追及 (monotone plan §5.1 #9)。
  - 単調拘束の撤去。事前試算では、どの条件でも山は 4 割残る。

## 3. 関連 docs と前提

- **θ_r の導出**: 軸の近くでは θ ≈ θ_r · r で、sinθ/r → θ_r。等エントロピー流では dν = √(M² − 1) · d ln u (Prandtl–Meyer の定義)、連続の式は d ln(ρu) = (1 − M²) d ln u。軸の近くの質量保存 2ρu θ_r = −d(ρu)/dx から、θ_r = ½ √(M² − 1) · dν/dx = ½ √(M² − 1) · (dν/dM) · dM/dx。
  - 比熱一定なら dν/dM = √(M² − 1)/(M(1 + (γ − 1)M²/2)) で、θ_r = ½ (M² − 1)/(M(1 + (γ − 1)M²/2)) · M′ に戻る。
  - semi-perfect では、MOC 本体と同じ気体モデルの ν(M) 表を使う (`moc_kernel.dnu_dM` は表の 3 次スプライン微分を返す)。γ の温度依存は ν(M) の中に入る。
  - methods の「今後の課題」にある ∂θ/∂r|r=0 = −½ d ln F/dx と同じ量。
- **事前試算** (monotone plan §9 2026-10-07、`case/45.isobutane_m6_d155/throat_moc_fix_probe*.py`、CFD 0 step):
  - 生産の n_start 41 で、第 1 点 (x = 0.025) の角度差は、生産の手順・dx0 0.03 で 0.022°。
  - 試験用の手順 (修正子 20 回 + 軸上のソース項 0) では 0.018°。修正子だけを収束させた版 (conv) でもほぼ同じ (n_start 321・dx0 0.015 で 0.019° 対 生産 0.022°)。効き目の大半は修正子の収束にある。
  - 軸上のソース項 0 は物理的には誤り (正しくは θ_r ≈ 0.10 rad/r_t、x_A で)。試験では、組み合わせで値が変わる不整合を消すために使っただけ。
  - 所要時間 (design_chain 1 回): 生産 2.3 s、試験用 (修正子 20 回固定) 8.6 s。
- **A11 の記述との関係**: methods は「`AXIS_LIMIT_FRAC` による代用は生産経路で 1 度も発火しない」と書いている。一方、軸上の端点 (r ≤ 1e-9) は常に代用の分岐に入る。前者は軸外の点のしきい値の話と読める。本 plan で methods の記述を正確にする。

## 4. 設計方針

### 4.1 軸上の極限 θ_r

- 軸上の端点 (`r_p ≤ 1e-9`) の sinθ/r を、相手の値ではなく、その点の θ_r にする。
- θ_r = ½ √(M² − 1) · ν_M(M) · M′(x)。M と M′ は軸上の M 分布 (target law) から取る。初期線の軸端 (x_A) では、アンカーの M_A・M′_A を使う。
- M′ は軸則の解析微分を使う (law が微分を持たなければ、軸節点の ν を 3 次スプラインで微分する。どちらを使ったかを記録する)。
- 軸上の点に θ_r を持たせる経路: 軸節点の配列に θ_r の列を足し、`interior_vec` (とスカラー版) が軸上の端点で使う。
- 軸外の点 (r_p > 1e-9) は従来どおり点自身で評価する。`AXIS_LIMIT_FRAC` の分岐は変えない。

### 4.2 予測修正の収束

- 修正子の反復を、θ と ν の変化の最大が許容値 (1e-12 rad) を下回るまで回す。上限は 50 回。
- 上限で収束しない対があれば、その数を記録して警告する。設計は止めない (反証可能にするため数を出す)。
- 使った反復回数の分布 (最大・平均) を `design_chain` の診断に出す。

### 4.3 既定値と切り替え

- 問題 YAML の `geometry.moc_axis_limit` (`legacy` | `analytic`、既定 `legacy`) と `geometry.moc_corrector` (`fixed2` | `converge`、既定 `fixed2`)。
- 既定はビット同一。case/45 の生産に入れるかは §6 V4・V5 の結果で決める (ユーザ判断)。
- 他の case の既定を切り替えるかは、本 plan の範囲外 (別途判断)。

## 5. 実装ステップ

1. `methods/design/overview.md`: 軸上の極限と修正子の記述、A11 の「発火しない」の記述を正確にする。
2. `design/forge_design/geometry/moc_kernel.py`・`moc_inverse.py`: §4.1・§4.2。
3. `design/forge_design/evaluate/runner_axismach.py`: キーを読んで渡し、診断を記録する。
4. `design/tests/run_moc_axis_limit_tests.py` (新規): θ_r の式 (CPG の閉形式との一致、semi-perfect は ν(M) 表と連続の式の数値微分との一致)、既定のビット同一、収束判定。
5. §6 の V2〜V6 を回す。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の諮問 | codex diagnose (`~/.config/forge/diagnose-backend` = codex) | F |
| 2 | codex plan 段レビュー | `codex_review.py <本 plan> --stage plan` | O |
| 3 | 実装 | §5 の 1〜4。合格条件: 新テスト FAIL 0、design/tests の既存テスト FAIL 0 | O |
| 4 | 検証 | §6 の V2〜V6 | O (解釈は F) |
| 5 | case/45 の生産への反映 | V4・V5 の結果を見てユーザが決める | ユーザ |
| 6 | codex result 段レビュー | `--stage result` | O |

## 6. 検証

### 6.0 事前登録 (初稿 2026-10-07、諮問の前)

- **V1 既定のビット同一**: キー無しで、case/45 の設計壁 (単調壁の問題) の当てはめ後の係数とノット、MOC 点群が現行と完全一致。design/tests の既存テストが FAIL 0。
- **V2 厳密解 (放射源流)**: A11 と同じ厳密解の試験で、新しい手順 (analytic + converge) と現行を比べる。誤差の最大と、点数を倍にしたときの収束次数を記録する。合格: 誤差が現行以下、かつ収束次数が現行 (1.7) 以上。
- **V3 再出発の自己整合**: `throat_moc_restart_test.py` と同じ試験 (網の列を初期線に戻す) を新しい手順で回す。合格: 列の 1 段目の作り直しの θ の差 ≤ 1e-6° (現行は −0.05〜−0.11°)。
- **V4 case/45 の形状の試算** (CFD 0 step、単調壁の生産問題、n_start 41・n_axis 2400・dx0 0.03): 第 1 点の角度差、単調拘束なしの r″ の山、設計壁の変化 (最大・位置・出口)、試験部の壁角の変化、出口半径の変化を記録する。
  - case/45 に入れる条件: 壁の変化 ≤ 10 µm、出口での変化 ≤ 0.1 µm、単調壁の形状ゲート S1〜S8 が PASS。外れたら諮問。
- **V5 Euler の確認** (V4 で入れる場合だけ): 生産 Euler 格子 G1 で 1 本 (run_0114 の場から番号写像、run_0113/0114 と同じ段)。出口コア M の変化が |ΔM| ≤ 1e-4 なら出口較正 (`Md_moc_offset`) は据え置き。超えたら較正をやり直す。
- **V6 所要時間**: design_chain 1 回の時間が現行の 5 倍以下。修正子の反復回数の最大・平均と、上限 50 回で収束しなかった対の数を記録する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `design/forge_design/geometry/moc_kernel.py`・`moc_inverse.py`・`design/forge_design/evaluate/runner_axismach.py` (オプション追加、既定はビット同一)。
- 逆 MOC を使うすべての設計 (case/41・42・44・45・46 など) は、キーを足さない限り変わらない。
- `methods/design/overview.md`。

## 8. 完了条件

- [ ] `methods/design/overview.md` を更新
- [ ] §6 V1〜V6 (V5 は入れる場合だけ) の結果を §9 に記録
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done` にして accepted へ移動し、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-07` — 初稿。ユーザ決定「予測修正を 20 回まで回す手順をメインに」「θ_r は semi-perfect で」。軸上のソース項は、試験で使った 0 ではなく解析極限 θ_r にする。
