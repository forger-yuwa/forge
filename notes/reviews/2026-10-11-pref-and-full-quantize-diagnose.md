# codex 諮問 (diagnose): pref-and-full-quantize

- **brief**: [`notes/reviews/briefs/2026-10-11-pref-and-full-quantize.md`](../../notes/reviews/briefs/2026-10-11-pref-and-full-quantize.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `461f3817` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.6 min, rc=0
- **結論**: **次は同じ FP64 バイナリで全 9 保存量の末尾量子化 G9 対 C だけを 2,000 step で判定し、pRef は今回の事前登録から外す。**
- **extra**: `case/45.isobutane_m6_d155/rq_an.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical なし）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| — | **§6.41 の限定した結論：採用** | 保存済みの `rq_judge.json` から区間を再計算し、e_E = 0.00576、区間 [−0.00679, 0.01831]、e_F = 0.91254、区間 [0.57429, 1.25080] を確認した。「登録した初期温度差の成分は、ρE 単独の末尾量子化では再現しない」は支持される。根拠：[集計結果](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rq_judge.json:98)。 |
| **Major** | **「主な成分は残差演算側」：却下** | 未試験の 8 保存量と相互作用が残る。また F と C は初期化後の保存量も異なり、自分の T₀ を引いても、その初期差が後で増幅する効果は消えない。さらに FP64 ビルドにも `ST=float` の経路があり、「演算は double」も正確でない。対案は全 9 保存量の末尾量子化の十分性だけを先に判定すること。根拠：[初期化差](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rq_judge.json:37)、[plan §6.36](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:898)。 |
| **Major** | **pRef の陽性結果から「積・足し込みが主な通り道」とする解釈：却下** | `pRef` は内部面・境界面・hoop ソースを同時に変える。圧力積の丸め、残差の加算、幾何の閉包誤差、hoop との相殺を分離しない。対案は「この圧力ゲージ変更で登録成分が減った」に限定すること。Cp−C は単なる記録ではなく、事前に許容幅を決めた同等性確認が必要。根拠：[内部面](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:618)、[境界面](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:197)、[hoop](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:52)。 |
| **Major** | **pRef の判定式：却下** | `U_Fp ≤ 0.5` では、例えば区間 [−2, −1] という大きな逆符号の偏りまで「抑制」になる。また 0.5 は固定した参照応答の半分であり、今回の e_F の半減ではない。対案：固定参照に対する小ささなら区間全体を [−0.5, 0.5] に収める。「既存 F から半減」を主張するなら e_F の不確かさも含めて比較する。射影の減少を全場の改善と混同せず、温度差の RMS・最大値も併記する。根拠：[ブリーフ:31–34](/home/sano/work/forge-integ-1005/notes/reviews/briefs/2026-10-11-pref-and-full-quantize.md:31)。 |
| **Minor** | **pRef の誤差見積もり：定性的部分のみ採用** | p̃ を作った後で pRef を引くので、それ以前の丸めは残る。ただし 0.25 Pa は 5.5 MPa 付近の値を一度丸める際の半 ulp で、p̃ の総誤差の上限ではない。提示された圧力範囲から、p−pRef は入口で約 −100〜500 Pa、収縮部全体では約 −182,000〜500 Pa、出口では約 −5,496,800 Pa。「数百 Pa」は入口付近だけに当てはまる。対案：「2〜数倍」という改善倍率の予測は外し、未測定とする。根拠：[圧力範囲](/home/sano/work/forge-integ-1005/notes/reviews/briefs/2026-10-11-pref-and-full-quantize.md:18)、[p̃ の生成](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:564)。 |
| **Minor** | **「前提はすべて満たした」：要再検証** | §6.40 は `check_quasisteady` の記録も要求するが、提示された `rq_an.py` と集計 JSON にその記録がない。対案は記録を補うか、未確認と明記すること。有限期間の十分性試験なので、準定常への到達自体を追加の合格条件にはしない。根拠：[登録した前提](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:984)、[判定器](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rq_an.py:50)。 |

結論: **次は同じ FP64 バイナリで全 9 保存量の末尾量子化 G9 対 C だけを 2,000 step で判定し、pRef は今回の事前登録から外す。**

第 1 仮説: **既存の FP64 演算経路に全 9 保存量の末尾量子化を加えるだけでは、登録した初期温度差の成分を再現できない。** 確度: **低**  
　根拠: `case/45.isobutane_m6_d155/run_0598_rq_e1/`・`run_0599_rq_e2/`・`run_0620_rq_e3/` の集計では、ρE 単独の e_E の上限が 0.01831。これは限定的な傍証であり、全 9 保存量についての証拠ではない。  
　反証条件: G9−C の区間下限が 0.5 以上で、陽性対照と介入の実効性が成立すること。

第 2 仮説: **他の保存量の丸め、またはそれらの相互作用を加えると、末尾量子化だけで登録成分を生じさせられる。** 未確認。  
第 3 仮説: **float 固有の初期化差の増幅が応答の一部を担う。** 未確認。C−F 比較だけでは除外されない。

判別 A/B: **変更は全 9 保存量の末尾量子化 OFF／ON の一点。**

- C は既存の `case/45.isobutane_m6_d155/run_0595_rq_c1/`〜`run_0597_rq_c3/`、G9 は同じ FP64 バイナリ・入力・設定で独立に 3 本、各 2,000 step。`FORGE_DIAG_QUANTIZE_STATE_ONLY` は明示的に解除する。`pRef` は変更しない。
- 評価域・参照 v_s・時点 1,000／1,500／2,000・run ごとの a は §6.40 の定義を固定する。
- 比較を **e_G9 = mean(a_G9)−mean(a_C)** と、既存陽性対照 **e_F = mean(a_F)−mean(a_C)** の二つに明記する。Welch 区間は各比較の t 分位を **1−0.05/4** とする。これは今回の二比較に対する補正であり、過去から続く探索全体の 95 % 保証ではない。
- **結果 A**：L_F ≥ 0.5 かつ L_G9 ≥ 0.5。第 1 仮説を退ける。「この介入だけで登録成分の少なくとも半分を生じさせられる」まで。
- **結果 B**：L_F ≥ 0.5 かつ [L_G9, U_G9] ⊂ [−0.1, 0.1]。第 2 仮説を退ける。「9 保存量の末尾量子化だけでは不十分」まで。残差演算への局在化はしない。
- それ以外は判別不能。延長・閾値変更で救済しない。

　投入前に C の再利用条件として、バイナリ・入力・実効設定と res₀ の 9 保存量の一致を確認する。G9 の自己点検は量名まで検査し、単に「合計 0」の文字列だけで通さない。特に [plan §5.1 #19](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:181) にある **入口ピン後の `roY_N` が float32 で表せない可能性**は未解決である。量子化直後の全 CV の検査と、その後の境界上書きの検査を分け、許容する場所・時点を結果を見る前に固定する。

やらない方がよいこと: **今回 pRef と G9 をまとめて走らせ、陰性結果を「状態ではなく残差」と読むこと。** 二つの試験を同じ文書に登録すること自体は可能だが、提示案の「三比較」は対象が明確でない。仮に F の陽性対照まで三比較補正で再判定すると、既存値から得る下限は **0.49561** で、0.5 を満たさない。比較集合と過去の陽性対照の扱いを後から選んではいけない。

呼び出し側の前提への異議: **観測されたのは有限期間の応答である。** 保存済み集計にある C・E・F 全 9 本の収束判定は **`OVERALL: CHECK FAILURES ABOVE`**。今回確認した結果 B を、定常解・下流 θ_r・全場の一致や原因確定へ広げることは受け入れない。

不足情報: 対象の run 本体はこの作業ツリーに無く、HDF5・残差履歴・実効設定の独立再検査はできなかった。確認できたのはコード、計画、台本、保存済み JSON とその区間計算である。`check_quasisteady` の記録も未確認。**ファイル変更・forge 起動は行っていない。提案は plan 未反映で、反映先は `plans/active/architecture-float-state-double-geometry.md` の §6.41 に続く事前登録節。**
