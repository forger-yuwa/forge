# codex 諮問 (diagnose): conjugate-benchmarks-results

- **brief**: [`notes/reviews/briefs/2026-09-30-conjugate-benchmarks-results.md`](../../notes/reviews/briefs/2026-09-30-conjugate-benchmarks-results.md)
- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **date**: 2026-09-30
- **commit**: `f8270e3e` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **まず評価器を登録条件に合わせて既存出力を再評価し、C は (a) の判定不能として保留する；現状で完了承認を求める result レビューには進めない。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0／Major 4／Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **現行評価を登録条件どおりの結果として扱う：却下** | 登録した領域切断の不確かさは [plan:102](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:102) にあるが、[eval_conj.py:9](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:9) は意図的に除外している。同じ領域で比較することと、領域切断感度がゼロであることは別。**登録した感度を評価し、各量の U に加える**。省略するなら事後の設計変更として採否を記録する。 |
| **Major** | **「全長・全節点で確認済み」：却下** | [eval_conj.py:196](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:196) は A の評価を加熱区間と −40R 以降の予熱域に限定。[series_conj.py:67](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/series_conj.py:67) は熱流束が尺度の 5% 未満の節点を除外する。**温度は登録どおり全長、準定常は対象全節点で再評価する**。小さい熱流束は登録尺度に対する変動幅で扱い、除外しない。 |
| **Major** | **A の r64 を代表結果にする：条件付き採用。r16 を単に「判定不能」とする：却下** | `case/64.conjugate_pipe_wall/run_0005_a1_r16/EVAL_CONJ.txt:3` は温度 **FAIL**（差 1.1905%＋U 0.28637%）。`run_0011_a2_r16_df20/EVAL_CONJ.txt:4` は熱流束 **FAIL**（1.7807%＋0.53325%）。[eval_conj.py:230](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:230) が別項目の判定不能を優先表示している。**「不合格項目あり、一部判定不能」と併記する**。r64 主判定・r32 補助・r16 粗格子限界という整理は、再評価後の**事後選定**として明記し、事前登録の総合 PASS に遡及させない。 |
| **Major** | **C：(a) 判定不能として検証成立を主張しない、を採用** | `case/65.conjugate_flat_plate/run_0007_c1_n64/CHT_INTERFACE_VERDICT.txt:3` は `res_abs=677.80 W/m²`、`res_rel=0.073939`、`dTw_max=0.0040952 K` で①〜④すべて不合格。`run_0010_c2_n64` も `607.03 W/m²`。**約 5.8 W/m²という要約は細格子の悪化を隠している**。両 run の収束判定は `NOT CONVERGED`、系列評価は `ALL STEADY`。後者で前者を免除できない。no-slip 化は流れそのものを変えるため、原因確認前の再登録・例外受入れは推奨しない。 |
| **Minor** | **A2 r16 の `Df_scale=20`：条件付き採用** | `case/64.conjugate_pipe_wall/run_0008_a2_r16/CONVERGENCE_CHECK.txt:2` は step 32300 で `NOT CONVERGED`。`run_0011_a2_r16_df20/CONVERGENCE_CHECK.txt:2` は `PASS (converged)`、G-if・G-cons も `VERDICT: PASS`。**変更履歴を残した別 run として使用可能**。ただし同一格子で収束した `Df_scale=5` の対照はなく、「固定点への影響を実測で除外済み」とは言えない。許可された閲覧範囲では、その実装上の不変性も未確認。 |

結論: **まず評価器を登録条件に合わせて既存出力を再評価し、C は (a) の判定不能として保留する；現状で完了承認を求める result レビューには進めない。**

第 1 仮説: **A の加熱・予熱域の粗格子誤差は、空間離散化誤差が主因である。** 確度: 中  
  根拠: `case/64.conjugate_pipe_wall/run_0005_a1_r16`／`run_0006_a1_r32`／`run_0007_a1_r64` の `eval_conj_600000.csv` を共通粗格子節点で比較すると、forge 自身の温度格子間差は **0.07211 → 0.01775 K（比 4.06）**。各 `CONVERGENCE_CHECK.txt` は `PASS (converged)`、`SERIES_CONJ.txt` は `OVERALL: ALL STEADY`。単なる forge–参照差の比より強い証拠である。ただし、現在の参照格子も forge 格子と同時に細分化されているため、誤差の帰属は確定していない。  
  反証条件: forge の保存場を固定し、参照格子だけ細分化すると、r16 の不合格が解消する場合。「forge の粗格子だけが原因」という説明は棄却する。

第 2 仮説: **C の細格子では連成反復が揺れを増幅している可能性がある。** 確度: 中。n32→n64 で G-if の `res_abs` は C1 が **5.7708→677.80**、C2 が **6.3217→607.03 W/m²**。slip 欠陥が起点か、連成が独立に不安定かは未確認。

判別 A/B: **A1 r16 の保存場を固定し、参照格子の最大細分回数だけ `--levels 2 → 3` に変える。** forge の追加計算は **0 step**、参照計算を一回追加する。原 run を保持した複製で、温度・熱流束の差、U の内訳、項目別判定を見る。  
  → **A：U が上限内まで下がっても不合格が残る**なら、「参照解の解像不足だけで説明できる」を棄却する。  
  → **B：判定可能になって合格へ変わる**なら、「r16 の不合格は forge 粗格子だけが原因」を棄却する。  
  U 超過が残れば未決であり、どちらの証拠にも使わない。

やらない方がよいこと: **C の閾値を緩める、窓外だから無害とする、上下流を同時に no-slip 化して原因を確定扱いする、r16 の FAIL を総括の「判定不能」で消すこと。**

呼び出し側の前提への異議: **「約 1/4」は全長の温度誤差には成立していない。** CSV の全長再集計では、A1 r64 の最大差は **0.132%**、A2 r64 は **0.134%**で、双方とも除外されていた下流端にある。現在算入された U との和はそれぞれ **0.229%、0.259%**で温度許容内だが、登録した全 U を含む判定ではない。また「窓内温度の振れ 2e−6 K」は時間変動の小ささであり、定常的な流れの誤りや界面残差の影響が小さい証拠にはならない。

不足情報: **登録どおりの全領域・全節点評価、領域切断 U、forge 側の固体厚さ方向温度差と軸方向伝導効果の評価、不合格となった G-if の節点位置・時系列、`Df_scale` の実効設定差分と固定点不変性の根拠。** これらを揃えてから result レビューへ進む。ファイルは変更しておらず、**plan 未反映**。
