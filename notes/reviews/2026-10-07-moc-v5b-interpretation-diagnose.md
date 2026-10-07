# codex 諮問 (diagnose): moc-v5b-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-07-moc-v5b-interpretation.md`](../../notes/reviews/briefs/2026-10-07-moc-v5b-interpretation.md)
- **plan**: [`plans/active/discretization-moc-axis-limit-and-corrector.md`](../../plans/active/discretization-moc-axis-limit-and-corrector.md)
- **date**: 2026-10-07
- **commit**: `8625cc4e` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.6 min, rc=0
- **結論**: **(b) を推奨し、保存場の熱力学量を独立に復元する 0 step A/B を先に行い、全温超過が保存状態にあるのか出力・後処理の不整合なのかを判別する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | 「同じ壁で IC により別の定常状態になる」 | **要再検証**。延長 7 本とも保存済みの区間判定は `NOT CONVERGED (stalled/plateau)`。ISEN の P 傾きは窓 B 全体・末尾 5 枚とも `DRIFTING` を再現した。対案は「異なる初期化履歴からの未整定状態に差が残る」と記述すること。根拠：[判定ファイル](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/v5b_series/run_0160_euler_icdep_mocG1_isen_ext36k/CONVERGENCE_VERDICT_segment.txt:1)。 |
| **Major** | V5 を閉じ、直ちに NS へ進む案 (a) | **現時点では却下、(b) を採用**。全温超過の生成経路が未確認のままでは、比較に使う Euler 場の健全性が分からない。対案は下記の **0 step の熱力学整合性 A/B**。全温超過の記録：[plan:210](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:210)。 |
| **Major** | MOC 修正自体を不採用にする案 (c) | **却下**。保存された V0 は壁誤差の次数 1.990、V4 は候補条件を満たしている。これは候補維持の根拠になるが、生産性能の証明ではない。対案は候補を維持し、生産採用・既定変更を保留すること。[plan:199](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:199)、[V4:203](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:203)。 |
| **Major** | 今回の結果で単調壁の採用根拠も無効になる | **一括撤回は却下、Euler の因果解釈は要再検証**。E′ はユーザが選んだ限定的な実務判定で、IC 独立性を証明していない。単調壁の限定採用には別途 N・K の記録もある。対案は採用履歴を残し、「Euler の差を定常な形状効果と認定できない」と限定すること。[E′:227](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:227)、[限定採用:378](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:378)。 |
| **Minor** | ブリーフの数値表現 | **訂正を採用**。窓 B の腕 B 平均は **0.208242** %pt。E′ の **37 % は (D＋2SE)/Δq** で、差 D 自体は約 **27 %**。対案は両者を分記すること。[評価 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/moc_v5b_ext_eval.json)、[E′ 数値:327](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:327)。 |

結論: **(b) を推奨し、保存場の熱力学量を独立に復元する 0 step A/B を先に行い、全温超過が保存状態にあるのか出力・後処理の不整合なのかを判別する。**

第 1 仮説: **初期化履歴を残した擬似時間反復の未整定が、P 傾きの差に寄与している。別の定常解への分岐とはまだ言えない。**　確度: **中**

  根拠: 以下は `case/45.isobutane_m6_d155/` の各 run に対応する保存 CSV を再集計した**未収束の窓平均**。

| run | 通算 24000〜36000 | 通算 42000〜54000 | 窓 B の P 傾き VERDICT：全13枚／末尾5枚 |
|---|---:|---:|---|
| `run_0157_euler_wallfit_mocG1_r1_ext36k` | 0.16094 | 0.16498 | `TRANSIENT-UNSETTLED / DRIFTING` |
| `run_0160_euler_icdep_mocG1_isen_ext36k` | 0.19826 | 0.21327 | `DRIFTING / DRIFTING` |

  単位は %pt。`check_quasisteady.py` の `classify_series` で再判定した。ISEN の窓間移動 **＋0.01501** は Δq/10＝0.003 の **5.0 倍**。腕 M 平均との差 **＋0.04867** は観測できるが、漸近的な IC 依存量ではない。

  反証条件: 独立復元で差が消えれば、保存場の未整定を主因とする説明は棄却する。逆に、熱力学整合性を確認した両系列が同じ定量的な定常条件を満たしても差が残るなら、「未整定だけで説明できる」は棄却する。

第 2 仮説: **壁近傍のエネルギー異常が初期場から継承、または反復中に生成され、履歴差を維持している。** 確度: 低、未確認。全温超過は plan の記録のみで、異常領域と η＝0.1 の圧力変化の結び付きは未検証。

第 3 仮説: **保存原始量・全エンタルピー・熱物性の復元経路に不整合がある。** 確度: 低、未確認。ただし通常の評価座標ずれを主因候補には戻さない。既存の座標感度検査で P 傾きの変化は最大 **3.9×10⁻⁹ %pt**。抽出コードも `centCoords` ではなく `MESH/COORD` を使用する。[nozzle_report.py:64](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:64)

判別 A/B: **変更点は「熱力学状態の取得経路」だけ、CFD は 0 step。**

- **A**：保存された `P,T,Y*,h0` を使用し、全温は既存の `total_quantities.total_state` で求める。
- **B**：同じ節点の `ro,roUx,roUy,roUz,roe,roY*` と、その run の解決済み熱物性・エンタルピー基準から、float64 で原始量と全エンタルピーを独立復元する。全温逆算には A と同じ熱物性関数を使う。ファイルへは書かない。

対象は上記 `run_0157`・`run_0160` の窓 A/B 各13枚と親の初期・終了場。座標、η、標本、P 傾きの定義は固定する。P 傾きの A/B 差、全温差、`h0` と `(roe＋P)/ro` の整合、壁からの層別最大・分位値・超過節点数を記録する。**診断閾値を P 傾き 0.003 %pt、全温 1 K と事前登録**する。

→ **両経路が全標本で閾値内、かつ数百 K の超過が再現するなら**、出力・復元経路だけが原因という仮説を退け、保存状態側の診断へ進む。  
→ **閾値を超えて食い違うなら**、現在の復元・出力の整合性を棄却し、その修正を先行する。両者の差を物理的な IC 依存と解釈しない。

やらない方がよいこと: **再延長、窓の移動、2SE による救済、出口較正変更、全温のクリップ**。また「µm 級だから Euler 比較は無意味」と一般化しない。現状では必要な精度の定常比較を成立させられていない、という結論までである。

呼び出し側の前提への異議: **全温超過は今回の判断に関係する。** 全腕にあることは相殺の証拠ではなく、中央値 1600 K も局所異常を免責しない。ただし、それを P 傾きの真因ともまだ呼べない。ユーザには、①V0/V4 の設計側の改善、②未収束・未整定という判定、③全温超過の検証結果、④未完了の NS・凝縮評価を並べるべきである。

不足情報: 対象 run の HDF5、実効 config・BC、解決済み熱物性記録、全温超過を計算した手順がこのツリーにないため、保存場の健全性と全温超過は独立確認できなかった。現在の HEAD は `8625cc4e` だが、V5/V5b 評価器の SHA256 は登録値と一致した。**ファイル変更なし。推奨は plan 未反映で、呼び出し側による反映が必要。**
