# codex 諮問 (diagnose): conjugate-A-result-disposition

- **brief**: [`notes/reviews/briefs/2026-09-30-conjugate-A-result-disposition.md`](../../notes/reviews/briefs/2026-09-30-conjugate-A-result-disposition.md)
- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **date**: 2026-09-30
- **commit**: `139e922e` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **(ii) を採用する事後改訂を記録し、延長部の源項だけを切り替える参照解A/Bで、延長感度の原因説明を確かめる。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（コード確認：`139e922e`）

| 指摘 | 採否・重大度 | 根拠と必要な処置 |
|---|---|---|
| M1 | **採用・Major** | [eval_conj.py:279](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:279) は各解自身の総入熱で比を作り、変種間の比の差から `U` を求めている。`q_hole` は [solidFem2d.cpp:334](/home/sano/work/forge-cht/solver_density_cuda/conjugate/solidFem2d.cpp:334) で積分済みの `W/rad`。修正方針は妥当。本番の再評価は未確認。 |
| M2 | **旧評価への指摘は採用、延長差を比較の `U` に含める処置は却下・Major** | **選択肢 (ii) を推奨**。[plan:53](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:53) の検証対象は、指定した流れ場・入口・有限領域の伝熱問題。一方、[eval_conj.py:92](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:92) の延長は入口位置・固体端面位置・源項の存在範囲を変更する。延長差は**採用した外挿条件に依存する問題設定感度**として別掲し、同一有限領域の比較から分離する。§4.6・§6の事後改訂として記録し、旧登録条件の判定履歴は残す。 |
| M3 | **採用。ただし対応完了は要再検証・Major** | 正式ツールへの接続と有限性検査は改善した。しかし [series_conj.py:89](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/series_conj.py:89) は `integral`・`dTs_norm`・`Qax_frac` を熱流束群へ入れ、全て `D=0.004` で判定する。登録の「比較許容の1/5」なら `Q_up/Q_tot` は **0.001**。固体指標も対応する許容幅から `dTs_norm=0.002`、`Qax_frac=0.001` として分けるべき。正式判定関数で、積分比の末尾増分 **0.0018** が現設定では `STEADY`、`D=0.001` では `DRIFTING` になることを再現した。 |
| M4 | **採用、未完了・Major** | 本番6本（`case/64.conjugate_pipe_wall/` の `run_0005/0006/0007/0009/0010/0011`）で、収束ツールは全件 **`NO residual_history.csv`**、総合 **`OVERALL: CHECK FAILURES ABOVE`**。本番の流体・固体保存場と `CONVERGENCE_VERDICT.txt` も無い。一次データの回収、評価コード版・入力ハッシュの記録、CSVと判定出力の同時再生成が必要。 |
| M5 | **採用。ただし対応完了は要再検証・Major** | [eval_conj.py:297](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:297) に効果の `U` と倍率判定は実装されている。ただし系列の閾値には上記M3の問題があり、新しい本番系列も未配置。修正後に両指標の正式VERDICTと倍率を出して完了とする。A2の「1〜2桁」という表現の撤回は妥当。 |
| m6 | **採用・Minor** | [README.md:15](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/README.md:15) は旧判定のまま。再評価後、有限領域・FP64・node・流れ場固定という限定と、事後改訂した基準での結果を同期する。 |

追加の **Minor**：[check_quasisteady.py:282](/home/sano/work/forge-cht/solver_density_cuda/tools/check_quasisteady.py:282) は `abs_scale` の正値・有限性を検査していない。`NaN` または `−1` を渡すと、単調に `0→1` へ変わる系列も `STEADY` になることを再現した。入口で **有限かつ正**を必須にする。今回の固定値 `1` の呼び出しには直接影響しない。

結論: **(ii) を採用する事後改訂を記録し、延長部の源項だけを切り替える参照解A/Bで、延長感度の原因説明を確かめる。**

第 1 仮説: 延長差の主因は、追加区間の `Φ + u·∇p` によって旧入口へ到達する温度分布が変わることである。 **確度: 中**  
　根拠: [eval_conj.py:119](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:119) は入口流れ場を外挿し、源項を計算して、入口温度条件を延長端に置く。`test_eval_conj.py` と同じ合成流れを使い、中間参照格子でメモリ内再計算したところ、延長による共通領域の最大壁温差は **0.1720 K、規格化差2.480%**。ブリーフの約2.46%を概ね再現した。これは本番runの実測ではない。  
　反証条件: 追加区間の源項をゼロにしても、元領域の最大壁温差が源項ありの場合の半分以上残る。

第 2 仮説: 固体端面・入口拘束位置の移動、または接続部の微分・離散化の変化が無視できない。**未確認**。延長時には元入口が内部節点となり、源項の微分も延長格子上で再計算される。

判別 A/B: **同じ延長格子・流れ場・BCで、追加区間だけの `S=Φ+u·∇p` を、A＝現状、B＝ゼロとする。** 元領域の源項は両者で同一に保ち、定常線形系を各1回解く。元領域の基準解に対する壁温差の最大ノルムを `D_A`・`D_B` として測る。  
→ **`D_B < 0.5 D_A` なら追加区間の源項が主因という仮説を支持／それ以外なら主因という説明を棄却し、第2仮説を残す。** このBは原因診断専用で、合格用参照解には使わない。

やらない方がよいこと: 散逸だけを止めて「問題を変えずに延長した」と扱うこと、`U` から延長項を外しただけで旧登録条件のPASSに読み替えること、現状の証拠でAを完了扱いにすること。

呼び出し側の前提への異議: **「粘性散逸が熱を足し続ける」は原因説明として不十分。** 合成場の断面積分では、散逸 **+0.04930 W/m/rad** に対し圧力仕事は **−0.04920 W/m/rad**で、ほぼ相殺した。壁温上昇だけから正味の加熱の累積とは言えず、温度分布の発達と分ける必要がある。また、正常合成ケースでも、比較差がゼロであることは別途設定した `U` の上限を満たす保証にはならない。

不足情報: 本番6本の一次出力、改訂後の全節点・固体指標系列と正式VERDICT、評価版・入力ハッシュの対応表。**ファイル変更なし。判断は plan 未反映であり、呼び出し側による§4.6・§4.7・§6・§6.1への反映が必要。**
