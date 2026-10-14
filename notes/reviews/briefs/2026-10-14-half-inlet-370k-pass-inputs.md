# 諮問ブリーフ: H (入口の半径を半分) の 370 K の最初の固定壁の pass の入力 (δ_r・積分法・Euler の参照) と事前登録の形 (2026-10-14)

- 依頼者: 主セッション (エスカレーション条件 1: §6 の検証計画を新規に書く、条件 4: 新しい run の前)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.9 (codex の採否まで)・§5.1 #33〜#36
- 前の諮問: `notes/reviews/2026-10-14-half-inlet-370k-design-diagnose.md` (H は物理壁の上流の条件を明示してチェーンに入れる、H の格子は Gref から、main の修正を先に、生成の経路の A/B (CFD 0 step) を最初に、370 K の止め方: stop_rule の状態の直接の確認・緩和の前の max|δ_E − δ_in| ≤ 0.05 mm と物理壁の変化 ≤ 0.05 mm をスロートから出口まで・時間の変動 0.01 mm・`band_select="edge"`・設計要求、名前は「H・B0 による暫定の設計」)
- ユーザ (2026-10-14): 入口の配管の半径を半分にした形 H を使う、370 K は B0 で進める、H の格子は Gref から。

## 観測事実 (#33・#34、commit 6bfbf7a0・01c4a625)

- origin/main (PR #7) をこのブランチへ merge した。design の試験 47 本で新しい失敗なし。
- `geometry.physical_upstream: {r_inlet: 3.2425, L_U: 9.0, L_pipe: 3.5}` を足した (物理壁の直管 r と Q の始点 (−L_U, r_inlet, 0, 0) だけが変わる。設計壁と CFD ピンの照合はそのまま、H の問題でも 16 項目が通る。制約: joint 壁の pw_upstream poly だけ、入口の位置は設計と同じ、μ ≤ 20 (H は 18.06)、δ_r は `delta_r_csv` でだけ渡せる)。キーの無い問題は HEAD とビット同一 (29 項目)。
- 生成の経路の A/B (CFD 0 step、`case/45.isobutane_m6_d155/_band_ab/half_inlet/half_inlet_ab.json`): 結果 A。壁 A − B は位置 3.3e−14・r′ 1.6e−14・r″ 7.3e−12 r_t、x ≥ 0 は C と同じ、第一層の距離の差 最大 0.072 %、品質 AR 4574・skew 0.32 (≤ 5000)。修正後の C は壁 609 → 610 点・終点 +5.353 mm、Gref の表の終端を延ばして ni 4719 のまま、station は旧から x < 0 で ≤ 9.2e−4・試験部で ≤ 0.069 r_t 動く。
- **見つかった穴**:
  1. **δ_r**: AWS の旧い δ_r の表 (`_band_ab/prod_confirm/prep/delta_r_initial.csv`、sha 0bc9da5b) は延びた壁を覆わず、`SingleBSplinePhysicalWall` が例外を出す。A/B では、生産の問題 (断熱) の積分法で修正後の設計から作り直した表 (`delta_r_initial_C.csv`、sha 24251cb6、numpy 2.5.2・scipy 1.18.0) を A・B・C すべてに使った。この表は断熱の積分法で、370 K の冷却壁のものではない。
  2. **積分法 (CONTUR) の初期壁と physical_upstream の併用**: 今は例外。積分法は設計の縮流部を前進積分する (x_virtual = L_pipe·r_t も 0.5 → 3.5 になる)。H の上流で積分するかは設計の判断。
  3. **δ* の抽出の Euler の参照**: core-matched の Euler の参照 (`run_0174` など) は旧い縮流部のまま。Euler の `prepare` は H を作れない (例外にした)。H の δ* ループでスロートの近くの帯局所の参照として許容できるかは未判断。
  4. カードの寸法 (wt_card・成果物) は H でも設計の 6.485 / 12 / 0.5 を出す (例外にしていない)。
  5. H の第一層の表は作ったときの δ_r の r_C/r_H で決まり、δ_r が変わっても作り直されない (影響は 1 % よりずっと小さい)。
- 370 K の C50 (今の形・G0c、V-tw370) の末尾: δ_target の全域の振れ 0.272 mm (x ≈ 1.27)、δ_E 0.588 mm (前の諮問で「受け渡しの条件を満たさない」)。

## 問い

1. 370 K・H の最初の固定壁の pass に使う δ_r: (a) 断熱の積分法で作り直した表 (`delta_r_initial_C.csv`) をそのまま、(b) 370 K の冷却壁の積分法 (CONTUR v2 など) で作り直す、(c) 300 K の NS の δ_E から、(d) ほか。最初の pass は「固定壁」で δ* の更新はしないので、どれを選ぶと何が言えるか。
2. 積分法と H の併用 (例外のままにするか、H の上流で積分する設計判断をするか)。
3. δ* の抽出の Euler の参照: 旧い縮流部の Euler (`run_0174`) を H の δ* ループで使ってよいか。使えないなら、H の Euler をどう作るか (Euler の `prepare` に physical_upstream を通すのは新しい機能)。最初の固定壁の pass でも抽出は要るか。
4. 最初の固定壁の pass の事前登録の形: 何を回して何を記録・判定すれば、#36 の δ* ループに進めるか (止め方の条件は前の諮問のもの)。run の数・予算・判定の項目を最小に。

## 読んでよいファイル

- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.8・§4.9・§5.1 #29〜#36・§6 (V-tw370 の結果)
- `case/45.isobutane_m6_d155/half_inlet_ab.py`・`_band_ab/half_inlet/half_inlet_ab.json`・`problem_d155_ns_prod_coldmesh_tw300_cgref_c.yaml`・`problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml`
- `design/forge_design/geometry/wall_axismach.py` (physical_upstream)、`design/forge_design/evaluate/runner_axismach.py` (`_physical_upstream`・`integral_delta_r`・`prepare_ns`)、`design/forge_design/feedback/deltastar_loop.py`
- `notes/reviews/2026-10-14-half-inlet-370k-design-diagnose.md`
