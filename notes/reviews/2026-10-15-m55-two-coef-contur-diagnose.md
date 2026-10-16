# codex 諮問 (diagnose): m55-two-coef-contur

- **brief**: [`notes/reviews/briefs/2026-10-15-m55-two-coef-contur.md`](../../notes/reviews/briefs/2026-10-15-m55-two-coef-contur.md)
- **plan**: [`plans/active/tooling-nozzle-card-local-mesh.md`](../../plans/active/tooling-nozzle-card-local-mesh.md)
- **date**: 2026-10-15
- **commit**: `36af1655` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 9.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-m55-two-coef-contur-diagnose.codex.log` (codex: ok)
- **結論**: 2 係数方式は採用し、実装・新規 CFD の前に、#3k の固定帯 `refine=4/32` A/B を全分布で行い、平滑化後の端点から得る係数と 2 係数壁の感度を確かめる。
- **extra**: `procedures/nozzle-design-workflow.md`, `case/71.m55_combustion_wt/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| ① 2 係数方式・pass 0 較正 | **採用。ただし係数確定は要再検証・Major** | 保存 CSV から **c＝0.852073279755、δ₀＝−0.322005099 mm** を再現した。初回の較正元は pass 0 を推す。pass 1 や両者の平均に替えても、壁・精度が異なる結果を混ぜるだけで誤差低減の根拠にならない。ただし「同じ CONTUR 壁」は抽出の正確さを保証しないため、下記 A/B を先に行う。 |
| ① 上流接続・寸法逆算 | **方式は採用。「波を作らない」は却下・Major** | `poly` は下流の値・傾き・曲率を端条件に上流多項式を再生成するため、正しく実装すれば C² 接続は保てる。一方、δ₀ は上流形状を変える。「上流に直接加算しない」と「上流が変わらない」は別である。[wall_axismach.py:567](/home/sano/work/forge-cardmesh/design/forge_design/geometry/wall_axismach.py:567)。寸法は x＝0 の半径ではなく**物理壁の大域最小**で解く既存処理を維持する。[deltastar_loop.py:443](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:443)。半径許容 1e−7 m は直径で 0.0002 mm 相当だが、最終生成壁・メッシュ・STEP での再測定が必要。 |
| ① 実装・カードの指定 | **現記載は修正・Major** | 計画の係数 2 キーだけでは `model` が欠け、`resolve_integral_initializer` が `ValueError` を返すことをメモリ上で再現した。[同:1189](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:1189)。`deltastar_initializer.model: contur` も overrides に指定する。変換は実効 `scale` を使った `c·δ_C＋δ₀/scale` とし、値・導関数・`.spline`・出力表を同じ変換結果から作る。[runner_axismach.py:1366](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:1366)。再較正時の δ_C は**係数適用前の CONTUR**に固定する。新 run の `delta_in_prev` は既に補正済みなので、そのまま基底に使わない。 |
| ② 途中の厚すぎと一様性 | **NS で判定を採用** | 保存 CSV で x＝3 m の差 **＋2.04386 mm** を再現した。端点一致は途中の有効壁の一致を保証しない。壁差の傾き・曲率から注意区間は探せるが、試験部 M の振幅を保証する感度は示されていない。M6 の合格実績だけで V5c 合格を予告しない。 |
| ③ V5b・V5c | **幅は候補として採用、判定方法は修正・Major** | V5b は **5.49725〜5.50275**。η＝r/r_w＝0.1 は既存手順にもある測線だが、軸の基準と同等ではない。`run_0005_ns/res_200000.h5` の再抽出では、最大絶対偏差／半振幅が軸で **0.76357／0.65550 %**、η＝0.1 で **0.53559／0.46504 %**。**この保存場では緩めになる**が、一般的な大小関係ではない。これは測線感度の比較であり、定常値の認定ではない。また既存 `wave_pct` は平滑曲線からの最大差で、V5c の半振幅とは異なる。[nozzle_report.py:216](/home/sano/work/forge-cardmesh/design/forge_design/report/nozzle_report.py:216)。指定窓で `max|100(M/5.5−1)|` と `100(max M−min M)/(2·5.5)` を直接評価し、**両量自身の時系列**を準定常判定に追加する。 |
| ③ V5a・V5d | **現状の合格定義は修正・Major** | [plan:109](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:109) のカード寸法 5 項目だけで、成果物全体を合格扱いにしない。既知の整合ゲート不合格は #3j で訂正・再判定する。停止規則系列には V5c の波高が含まれていないため、それだけの `STEADY` は一様性の定常性を保証しない。[stop_rule.py:37](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:37)。 |
| ④ Euler の使い回し | **回し直しを採用** | r_t の一致桁数による安全な閾値はない。カードは `R`・入口半径・各長さも r_t で割り直すため、r_t の近さだけでは設計・格子・参照補間の同一性を保証できない。[wt_card.py:1199](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:1199)。約 4 分なら新しい理想壁で作り直し、参照の出力と SHA を固定する。 |
| ⑤ #3k の省略 | **取りやめは却下・Major** | 使用する値が 2 点でも、**平滑化後の端点は全分布から解いた係数に依存する**。[deltastar.py:924](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:924)。2 断面だけの再抽出では現在の較正手順を検査できない。全分布の A/B を維持し、今回の進行判断は「その差が 2 係数壁にどれだけ残るか」で行う。 |
| ⑥ 初期場 | **run_0005_ns を採用** | pass 0 基準で比較すると、新候補の補正は pass 1 に使った δ_E から途中約 2 mm、pass 0 入力から出口約 7.62 mm 離れる。新候補に近い初期推定として pass 1 が妥当。ただし速い整定の保証ではなく、波が初期場由来かどうかも短い窓で判定しない。 |

結論: **2 係数方式は採用し、実装・新規 CFD の前に、#3k の固定帯 `refine=4/32` A/B を全分布で行い、平滑化後の端点から得る係数と 2 係数壁の感度を確かめる。**

第 1 仮説: **径方向標本化の誤差が平滑化後の端点にも入り、2 係数壁に無視できない差を残す。** 確度: **中**。影響量は未確認。

  根拠: 帯フィットは採用点の等重み最小二乗、欠損積分は `y <= y_b` の点だけを使い、積分端を追加していない。[deltastar.py:381](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:381)、[同:418](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:418)。その出力を全分布で平滑化してから較正するため、端点だけ使ってもこの経路は残る。

  反証条件: 下記 A/B で、較正端点・投影した 2 係数壁・壁角の差がすべて事前許容内なら、**今回の係数決定を妨げる大きさの標本化依存**という仮説を棄却する。

第 2 仮説: **2 係数では表せない途中の排除厚差が残り、V5b を満たしても V5c を外す。** 確度: **中**。途中差は確認済みだが、新壁の NS 応答は未確認。

判別 A/B: **変更は `band_local_deficit.refine` のみ、A＝4、B＝32、forge は 0 step。**

- 入力は `case/71.m55_combustion_wt/run_0003_ns_cont/res_200000.h5`、参照は `run_0001_euler/res_12000.h5`。保存帯、一定 x への補間、平滑化ノット・λ・重みを固定する。
- 既存の [refine_ab_3k.py](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0005_ns/refine_ab_3k.py:56) を再利用する。ただし現台本は腕ごとに重みを再計算する。hard 判定が変わった場合は単純比較を合格にせず、A の保存 raw 再現も必須とする。
- 全分布を平滑化し、同じ x＝0・出口で較正する。δ_C と設計壁は固定し、各腕の c・δ₀ から補正壁を作る。
- 今回の診断許容は、端点の **max|Δδ_E| ≤ 0.01 mm**、下流全域の **max|Δ(cδ_C＋δ₀)| ≤ 0.01 mm**、**max|Δ壁角| ≤ 0.01°**。元の #3k の途中差も記録する。
- → **いずれか超過**なら「2 点較正なら抽出感度を無視できる」を棄却し、抽出を先に直す。  
  → **すべて以内**なら、この標本化依存を新方式の阻害要因から外し、上表の修正を反映して進む。元の直接 δ_E 更新の原因診断は未決着として残す。

やらない方がよいこと: **#3k を省略して係数を確定すること、途中差だけを理由に係数を増やすこと、出口 M の合格で試験部全体を合格とすること。** 再較正する場合も、V5b 超過の原因が排除厚の端点差だと確認せず、機械的にもう一本回さない。

呼び出し側の前提への異議: **「係数の算術が合う」「端点が比較的落ち着いている」「新設計が合格する」は分ける。**

- 保存記録では両 NS とも本段 0〜199,999 が **`NOT CONVERGED`**。`run_0005_ns` について今回も再実行して確認した。出口コア M の系列 2,000〜200,000 は **`STEADY`** だが、波高の判定ではない。
- pass 0 の保存 δ 系列 160,000〜200,000、21 点を既存 `check_quasisteady.classify_series` で確認すると、スロート直後 x＝0.0997 mm と出口はいずれも **`STEADY`**。全窓の振れはそれぞれ **0.00300 mm、0.000759 mm**。端点の時間変動は比較的小さいが、前者は x＝0 そのものではなく、標本化誤差も検査していない。
- 「M6 は合格した」は方式選択の参考になるが、測定窓・波の定義が違うため、今回の V5c 合格の証拠にはならない。

不足情報: 保存元での #3k の結果、厳密な x＝0 の較正値の時系列、新壁の V5 判定。これらは係数確定・最終合格に必要だが、次の一手を選ぶ判断は妨げない。根拠 run の索引は [case README の計算 run 一覧](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:21)。

証拠: 足りる
