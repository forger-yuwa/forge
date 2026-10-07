# codex 諮問 (diagnose): throat-mono-result-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-07-throat-mono-result-interpretation.md`](../../notes/reviews/briefs/2026-10-07-throat-mono-result-interpretation.md)
- **plan**: [`plans/active/tooling-nozzle-throat-monotone-r2.md`](../../plans/active/tooling-nozzle-throat-monotone-r2.md)
- **date**: 2026-10-07
- **commit**: `479807e7` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.2 min, rc=0
- **結論**: **生産キーの反映前に、検証済み NS レシピと実際の生産入口を比べる 0-step の準備 A/B を１件実施する。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 2 / Minor 2）

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 「S・E′・N・K を満たした」 | **採用。ただし検証条件内の判定** | S1〜S8 の保存結果は `PASS`。CSV から `check_quasisteady.py` を再実行し、dry の連結窓 60000〜80000、凝縮の窓 14000〜18000 は各４量とも `STEADY` を再現した。保存済み残差判定は `NOT CONVERGED (stalled/plateau)`、RISING なし。したがって「登録した数値ゲート合格」は支持するが、「収束」「性能認定」とは書かない。 |
| 「Euler・dry NS・凝縮 NS のどれでも旧壁との差が許容幅と時間変動の範囲内」 | **却下／Major** | dry の η0.1 の場の差は最大 **0.002199 %pt**、両腕の末尾変動幅の最大は **0.001675 %pt**。既に後者を超える。さらに [比較コード](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_ns_compare.py:51) の帯は各腕の最大最小幅であり、差の信頼区間ではない。N は絶対性能ゲート、K は準定常ゲートで、旧壁との差の許容幅を設定していない。**「Euler は E′ 合格、dry は N 合格、凝縮は K 合格。新旧差は参考値」**に分ける。 |
| 「生産 YAML にキーを入れれば検証結果を引き継げる」 | **要再検証／Major** | 検証経路は [prep_c2pin.py:11](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/prep_c2pin.py:11) で較正済み `initializer` を明示する。対象 YAML には `deltastar_initializer` がなく、通常の `prepare_ns` で引数も省略すると、[runner_axismach.py:981](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:981) の積分法経路に入らず、`correlation_hist_v1` 側になる。**実際に生産で使う入口が検証済み物理壁を再現することを、0-step で確認する。** |
| NS の番号写像・段階起動省略 | **採用** | `run_0147/IC_MAP.json` は全検査 `OK`、９保存量ビット一致、幾何・`wall_dist` 保持、最大移動 0.499189 µm。旧新 dry、旧新凝縮の実際の `solverConfig.yaml`・`bcondConfig.yaml` はそれぞれ同一。投入前変更なので、それ自体を結果後の条件緩和とは扱わない。ただし IC 非依存性や通常起動経路の検証にはならない。 |
| E′ の説明・plan の整合 | **修正採用／Minor** | 実務判定への変更はユーザ決定として有効。ただし自己相関未補正・窓の事後選択により、`2SE` を保証された信頼限界として扱えない。P 傾きの **差 D は許容幅の27.3%**、**D＋2SE が37.0%**。また [plan:243](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:243) の「検出できなかった」と同277行の「E3採用」は、実際の E′ と検出結果に合わせて訂正する。 |
| 「N 合格後に K」の手順 | **結果は採用、履歴訂正／Minor** | [投入スクリプト:36](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_mono_ns_chain.sh:36) は N の準定常合否を確認せず、`run_0147/res_60000` から K を開始している。その時点のオーバーシュートは `DRIFTING`。K 自身の４量 `STEADY` は再現したので直ちに再計算を要求しないが、「登録順序どおり完了」とは記録しない。 |

結論: **生産キーの反映前に、検証済み NS レシピと実際の生産入口を比べる 0-step の準備 A/B を１件実施する。**

第１仮説: 残る主要な穴は、検証結果の失敗ではなく、較正済み物理壁を作る条件が生産入口に引き継がれる保証の不足である。確度: **高（経路分岐は確認済み、生産入口での発現は未確認）**  
　根拠: `initializer` の明示有無で積分法と相関法が分岐する。一方、検証済み `run_0147` は `cf_scale = 1.054129117086371` の積分法を使用している。  
　反証条件: 実際の生産入口も同じ較正値を解決し、同じ環境で設計 spline、δ_r、物理壁、変換後メッシュを再現する。

第２仮説: 新旧 NS 差には壁変更以外の継続時間・初期化履歴の影響も含まれる。確度: **中、寄与は未分離**。旧壁の波は `DRIFTING` で、同一入力再実行幅も未評価なので、微小差を壁変更だけに帰属できない。

判別 A/B: **変更点は `initializer` の供給経路だけ、CFD は両腕とも 0 step**。同一の単調壁 YAML・環境を使い、A は検証時の較正済み辞書を明示、B は生産入口が通常解決する値を使う。他の引数は固定する。供給値、spline、δ_r、物理壁、float32 化後のメッシュ座標・接続を比較する。  
　→ **一致なら**、準備経路の不一致説を棄却し、表の記述修正と正式 result レビュー後に、当該条件に限定して生産採用する。  
　→ **不一致なら**、YAML キー追加だけで検証条件を引き継げる前提を棄却し、生産レシピを修正する。

やらない方がよいこと: 現段階で NS の同一入力再実行や MOC 修正を一律に追加すること。今回の限定した受入判断には先に上記の準備確認が必要であり、NS 差の統計的非劣化まで主張する場合に初めて再実行評価が必要になる。

呼び出し側の前提への異議: **数値ゲート合格と、差が検出できないことは別である。** 新壁の出口 M は 5.998870574、下限まで **7.06e−5**。登録した評価法では合格だが、既知の格子・標本不確かさ約 ±1e−4 より余裕が小さい点は残る。「旧壁の下限境界問題まで解消した」とは解釈しない。

不足情報: ローカルには新規 NS run の `res_*.h5`・残差 CSV、Euler 各 run の原系列がないため、全保存場の NaN と Euler 統計は独立再計算できなかった。NS 準定常は CSV から再判定し、残差・NaN・壁解像・メッシュ品質は保存記録を確認した。

根拠 run は `case/45.isobutane_m6_d155/run_0147_ns_mono_final/`、`run_0149_ns_mono_final_ext/`、`run_0148_ns_mono_final_cond/`。主要根拠は各 CSV、`CONVERGENCE_VERDICT_segment.txt`、`report/report.json`。品質記録は **`VERDICT: PASS (AR<=5000, skew<=0.90)`**、壁解像は **PASS（超過面積3.5%）**。[run 索引](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:21)。

ファイル変更・forge 起動は行っていない。上記提案は **plan 未反映**。
