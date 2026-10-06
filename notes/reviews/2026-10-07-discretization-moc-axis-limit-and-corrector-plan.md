# codex レビュー: discretization-moc-axis-limit-and-corrector (plan)

- **plan**: [`plans/active/discretization-moc-axis-limit-and-corrector.md`](../../plans/active/discretization-moc-axis-limit-and-corrector.md)
- **stage**: `plan`
- **date**: 2026-10-07
- **commit**: `e1041e99` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1
- **extra**: `notes/reviews/2026-10-07-moc-axis-limit-diagnose.md`, `plans/accepted/discretization-moc-axisymmetric-source-term.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
軸上極限の式と課題の同定は妥当です。独立のメモリ上試験でも、放射源流の第1段・壁の精度改善を確認しました。  
ただし、IC 写像、CFD の採否条件、数値失敗の扱い、生産化の順序を実装前に修正してください。

1. **Major — 新しい壁変化に、旧 plan の IC 写像条件をそのまま適用できません。**

   **根拠:** 対象 plan の V4 は壁変化を最大 10 µm まで許容し、V5′ は旧 plan の手順を継承します。一方、旧手順の番号写像は移動量 **1 µm 以下**に限定され、実装も超過を拒否します（[ic_index_map.py:444](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ic_index_map.py:444)）。

   今回、既存の生産問題と保存済み初期線を用いた CFD 0 step の試算では、`analytic+converge` による設計壁の最大変化は **6.181 µm、x/r_t ≈ 0.400**でした。これは設計壁の値で、NS 物理壁・変換後メッシュの移動量は未確認ですが、旧手順を無条件に継承できない具体的な兆候です。

   **対案:** V5・V5′の前に、変換後メッシュの移動量・接続・境界対応・要素反転・品質を検査する工程を追加してください。番号写像は今回の変形に対して条件を検証し直し、保存量の直接転送と IC 依存の確認を行う。**上限値だけを引き上げて通すことは推奨しません。**

2. **Major — V5 は較正の据置条件であり、設計変更の採否条件として不足しています。**

   **根拠:** [対象 plan:124](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:124) は出口 M の差だけを判定し、「比較の不確かさ」の算出方法、対象量の準定常条件、欠損時の扱いを定義していません。出口平均が維持されても、M 波・P 波・オーバーシュートは悪化し得ます。

   また、参照する `run_0143〜0145` の旧採用は、自己相関未補正・事後選択窓を含む**限定的な実務判定**です（[旧 plan:227](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:227)、同:243）。その例外を今回の比較へ自動継承できません。今回、この Euler 系列の実データと VERDICT は取得できておらず、CFD の収束・定常性は独立確認していません。

   **対案:** 出口 M の条件は較正判定として残し、旧 E2 相当の波・傾き・オーバーシュートについても許容悪化幅を事前登録してください。評価窓、標本感度、時間変動・再実行変動からの不確かさ算出を固定し、全残差の `check_convergence` と対象量の `check_quasisteady` の VERDICT を判定の前提にする。**不確かさが大きく判別不能な場合は、較正変更ではなく保留**にしてください。

3. **Major — 「有効な対だけで停止判定」の定義によっては、失敗した対を除外して合格できます。**

   **根拠:** [対象 plan:83](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:83) は未収束数 0 を要求しますが、有効性の定義がありません。現在の `interior_vec` は平行・非有限の対を `ok=False` にし（[moc_kernel.py:435](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_kernel.py:435)、同:450）、`fill_levels` は棄却した結果を NaN に置換します（[moc_inverse.py:142](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:142)）。既存の欠損と、今回の反復が起こした破綻を同じマスクで除くと、失敗が未収束数に残りません。

   最終の幾何・適合式残差も「記録」のみで、合否条件になっていません。

   **対案:** 「入力時点の対象対」「既存の欠損」「幾何的棄却」「反復失敗」「収束」を分離して集計してください。対象対で発生した NaN/Inf・反復上限到達は不合格とし、幾何的棄却は許容領域を明示する。最終状態の残差にも許容値を設け、NaN・平行特性線・上限到達を入れた負例テストを追加してください。

4. **Major — 残作業表では、生産反映が result レビューより先です。**

   **根拠:** [対象 plan:109](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:109) の生産反映条件は V4・V5 のみで、次行に result レビューがあります。V5′は本文に存在するものの、生産反映の必須条件として表に結び付いていません。継承元は N・K と result レビューが揃ってから生産化する条件です（[旧 plan:277](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:277)）。

   **対案:** 順序を **V0 → 実装・回帰 → V4 → V5 → V5′ → result レビュー → 生産採用判断**に統一してください。生産不採用で機能実装だけを完了する場合と、生産採用まで完了する場合も分ける。完了条件には V0 と V5′を明記してください。

5. **Minor — §2に、analytic に対して未確認の「山は4割残る」が残っています。**

   **根拠:** [対象 plan:40](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:40) と、同:51 の「analytic の残る山についての観測ではない」が整合していません。今回の試算でも第1点の角度差は **0.021884° → 0.005441°**となり、K2c の観測を analytic に流用できません。なお、これは角度差であり、r″ の山の残存率ではありません。

   **対案:** 「単調拘束は維持する。analytic 適用後の山の残存率は V4 で測定する」に修正してください。

推奨は、**解析極限と収束修正子の採用方針を維持し、上記 Major を番号順に plan へ反映してから実装すること**です。

この推奨には数値的な裏付けがあります。[既存の放射源流試験:268](/home/sano/work/forge-integ-1005/design/tests/run_inverse_tests.py:268)を使い、既知軸端の源項だけをメモリ上で変更しました。壁誤差は x=1.4〜2.4 の既存21標本での最大相対誤差です。

| n_axis / n_start | A: legacy＋converge | B: analytic＋converge |
|---|---:|---:|
| 140 / 25 | 9.7648e−5 | 7.0373e−5 |
| 280 / 49 | 2.9493e−5 | 1.7434e−5 |
| 560 / 97 | 8.7378e−6 | 4.3688e−6 |
| 1120 / 193 | 2.5262e−6 | 1.0996e−6 |

B の最細区間の観測次数は **1.990**。軸節点2点から生成する第1段の最大 θ 誤差は、全解像度で A の **4.37e−4 倍以下**でした。両腕とも更新量 1e−12 以下に達し、上限到達の未収束対は 0、最大反復数は 35 です。したがって、**V0 の主要な精度条件は独立試算で支持されます**。正式実装の回帰・再出発・全形状ゲートの代用ではありません。

目的は A11 が明示的に先送りした残件で、解決済み機能との重複ではありません。対象は NumPy の float64 MOC であり、1e−12 を CUDA/float32 の制約と混同する必要もありません。FVM の node/cell・周期・block-DPLUR を直接変更しないため、検証は MOC の厳密解と case/45 の node 設計チェーンに集中する方針が妥当です。

ファイル変更なし。**plan 未反映**です。反映先は §2、§4.2、§5.1、§6、§8です。

指摘数: Critical 0 / Major 4 / Minor 1
