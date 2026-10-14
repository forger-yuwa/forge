# codex 諮問 (diagnose): core-grid-gref-result

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-gref-result.md`](../../notes/reviews/briefs/2026-10-11-core-grid-gref-result.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `4dac0ea4` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.1 min, rc=0
- **結論**: **今回の判定不能を維持し、次は Gref の半径方向条件を固定した x 粗化だけの A/B を、新しい資格基準と新しい観測窓で事前登録する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0／Major 5／Minor 1）

| 対象 | 重大度・採否 | 根拠と推奨 |
|---|---|---|
| (a) G1・G1x を候補から外す | **Major・要再検証** | **追加計算の優先対象から外すことは採用、不合格の確定は却下**。G1x の θ_r[40,50] は差 3.496 % に対して U = 0.496 % で約7倍。δ_loc の差は 0.846〜2.037 % で、「境界層量が一律3〜7 %」「揺れより1桁以上」は成立しない。さらに ≤0.14 %/4万 step は **G0・Gref の傾き**で、G1・G1x の長期ドリフトを拘束しない。「登録判定は判定不能。暫定差が大きいため、現候補への追加投資を見送る」と記録する。[plan:459](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:459)、[plan:470](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:470) |
| (a) G0 の結果の表現 | **Major・採用、表現を限定** | 「指定窓の平均差は θ_r +0.193〜+0.271 %、δ_loc −0.403〜+0.017 %、Q_w +0.335 %、軸の指標 +0.182 %pt」とする。**境界層の精度確認・軸誤差の確定とは呼ばない**。軸の B は s 抜きでも 0.189 %pt > τ = 0.02 %pt。G0 を次の生産候補として残す根拠はない。[plan:459](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:459)、[plan:467](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:467) |
| (b) 壁際維持・主流上限変更・x粗化をまとめる | **Major・却下** | **まず Gref の半径方向の生成条件を維持し、x だけ G1x 相当に粗くした `Grefx` を試す**。3086 × 160 = 493,760節点で、G0 比約0.865となり削減目的にも沿う。既存の暫定差から計算すると G1x−G1 は θ_r −1.09〜−1.32、Q_w −2.133 ポイントで、x の効果を無視できない。約41万節点案を一度に試すと相殺と原因の切り分けができない。[plan:310](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:310)、[plan:459](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:459) |
| (c) 資格閾値の再登録 | **Major・採用、条件付き** | τ に基づく**新しい試験の運用基準**として登録する。旧判定は変更しない。`check_quasisteady` の drift は窓全体の線形変化、osc は全振幅であり、d・U とは別物。推奨は共通分母で **drift ≤ τ/8、全振幅 ≤ τ/2** とし、既存の d ≤ τ/8・U ≤ τ/4、末端トレンド・漸近値の条件も維持する。既存 Gref 窓は探索比較には使えるが、新基準による確認には**登録後の新しい継続窓**を使う。[判定コード:280](/home/sano/work/forge-coregrid/solver_density_cuda/tools/check_quasisteady.py:280)、[plan:442](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:442) |
| (d) G0・Gref の s を今測る | **Major・却下、次の候補選別後に実施** | G0 は軸の暫定予算ですでに超過し、節点削減にもならない。今は G0 の point 継続に投資せず、次の比較で候補が残ったら **その候補と Gref** の s を測る。これは既登録の順序と整合する。s 未測定の段階では生産採用を認めない。[plan:444](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:444) |
| M の後処理予算 | **Minor・要再検証** | コードは候補・Gref とも出口 0.0006 %、軸等 0.0014 %pt を使う。一方、plan は Gref の実測値を 0.0001 %、0.0005 %pt と記載し、その後に「どの格子も同じ値」とも書いている。保守側の差だが、**実測値と採用する予算を区別して明記**する。実測値採用なら B は出口で0.0005ポイント、軸等で0.0009 %pt小さくなる。今回の資格不成立は変わらない。[集計コード:130](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_gref_judge.py:130)、[plan:445](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:445) |

結論: **今回の判定不能を維持し、次は Gref の半径方向条件を固定した x 粗化だけの A/B を、新しい資格基準と新しい観測窓で事前登録する。**

第 1 仮説: x 方向の粗化による量の変化は、半径方向を Gref 相当に細かくしても許容差 τ の水準に残る。確度: **中**。  
根拠: `case/45.isobutane_m6_d155/run_0603_cg_g1_k1/` と `run_0605_cg_g1x_k1/` の暫定平均差は θ_r 約1.1〜1.3ポイント、Q_w 約2.13ポイント。変更した x の密度と数値は上記参照先に記録されている。ただし未収束の比較であり、細かい半径方向格子への一般化は未確認。  
反証条件: 新しい両腕が資格を満たし、全指標で E + U_A + U_B + p_A + p_B ≤ τ となれば、「この条件でも x 粗化が比較予算を超える」という仮説を退ける。s を含む最終採用は別に判定する。

第 2 仮説: G1・G1x の大きな境界層量の差には、半径方向分布の変更が主要に寄与する。確度: **中**。G0−Gref の境界層量の差が相対的に小さいことと整合するが、「壁際だけ」が原因とは特定できない。

第 3 仮説: 反復の停滞・残存過渡が格子間差の一部を作っている。確度: **中**。s が未測定なので未除外。Gref の記録された VERDICT は **`NOT CONVERGED`**、160k〜200k の準定常判定は **`DRIFTING 5・OSCILLATING 5`**。G0 も同じ準定常判定である。[plan:452](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:452)

判別 A/B:

- **A**: Gref を同じ設定で継続。
- **B**: 半径方向の生成条件・物理壁・第一層の設計関数を Gref と同じにし、**x の密度だけ** G1x 相当に変更した `Grefx`。
- 両腕の起点は `case/45.isobutane_m6_d155/run_0607_cg_gref_k1/res_200000.h5`。A は同一格子 restart、B は cross-mesh 移送。バイナリ・空間離散化・反復設定・hoop キーは固定する。
- 各 **6万 stepを予算上限**とし、最初の2万を比較窓から外し、後半4万を2500 step間隔で保存する。2万で移送過渡が消えたとは仮定せず、後半の資格で検査する。θ_r・δ_loc・Q_w・出口と軸等の M を同じ定義で評価し、新格子の p も測る。
- **結果A**: 両腕が資格を満たし、全量で上記の s 抜き予算内 → x 粗化の大きな影響という仮説を退け、`Grefx` と Gref の s 測定へ進む。
- **結果B**: 両腕が資格を満たしても、いずれかで E − U_A − U_B − p_A − p_B > τ → 観測窓の揺れ・後処理では説明できない x 感度を支持し、この粗化幅での候補化を見送る。
- 境界域、または資格不成立なら**判定不能**。反復誤差の上限が不明なので、結果Bでも真の格子誤差の超過とは断定しない。

やらない方がよいこと: 旧 Gref 窓を緩めた閾値で合格に読み替えること、G1・G1x を精度不合格と確定すること、主流上限と x を同時に変えて改善を壁際分布へ帰属すること。変換器変更を揺れの主因として再調査することも、記録済み A/B の比 1.05・0.92 を覆す新証拠がない限り勧めない。[plan:415](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:415)

呼び出し側の前提への異議: **小さい U と短窓の傾きは、未収束の偏りの上限ではない**。また、τ/8 を CLI にそのまま渡してはいけない。CLI は各系列の平均を分母とする無次元比である。軸の τ = 0.02 %pt は、判定用の M/6 系列では 0.0002 に相当する。分母換算を含む規則を事前登録する。「12倍厳しい」は θ・δ の比較であり、Q_w や M には同じ倍率を使えない。

不足情報: 正式な τ、新格子の品質・壁解像・後処理予算、候補と Gref の s。依頼の読み取り制限に従い、生の場・残差・判定成果物は再検算しておらず、数値と VERDICT は指定 plan の記録に基づく。**ファイル変更・forge 実行なし。plan 未反映**。採用時の反映先は呼び出し側の §4.20・§5.1・§6。
