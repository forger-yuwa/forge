# codex レビュー: boundary-cht-conjugate-benchmarks (plan)

- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **stage**: `plan`
- **date**: 2026-09-29
- **commit**: `8410764d` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 6.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **extra**: `case/64.conjugate_pipe_wall/conjugate_ref.py`, `case/64.conjugate_pipe_wall/survey_A.py`, `case/65.conjugate_flat_plate/survey_C.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
検証目的は妥当で、既存 plan と重複しません。  
ただし、参照解の軸上離散化と条件選定の集計に誤りがあり、現在の事前登録値を確定値として実装に進むべきではありません。

既存の検証範囲は、[plans/README.md:77](/home/sano/work/forge-cht/plans/README.md:77) と [methods/boundary.md:350](/home/sano/work/forge-cht/methods/boundary.md:350) に記載された限定付きのものです。「流れがあり、固体抵抗・軸方向伝導も効く条件」を追加する意義はあります。`node`・定常陰解法・`fem2d`・`axisymMethod: 0` は既存構造と整合します。今回の結果を `cell`・FP32・周期境界まで一般化しない前提なら、検証範囲の選択も適切です。

以下の追加数値は、ファイルを変更せず実行した**参照計算・合成試験**です。case/64・65 の forge 検証 run の結果ではありません。

1. **Major — 参照解の軸ノードで、軸方向対流が消えている**

   **根拠:** [conjugate_ref.py:151](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/conjugate_ref.py:151) は軸方向質量流束の面積を `H * rw(ys[j])` としており、`r=0` ではゼロになります。しかし、軸ノードの双対セルは有限断面を持ち、その面積は per rad で
   \[
   \int_0^{\Delta r/2}r\,dr=\Delta r^2/8>0
   \]
   です。伝導・源項側にはこの有限体積が入るため、不整合です。

   合成場 `T=x−x²/2`、`ρ=c_p=u=1`、`v=0`、`k=0.1`、`S=1−x+0.1` を代入すると、`x=0.5` の軸行の体積規格化残差は、半径方向分割数 **8・16・32 のすべてで −0.5**。欠落した対流項そのものです。

   **対案:** 軸方向対流面積も流体部分の `∫r dr` で組み、非一様半径格子を含む製造解試験を追加する。既存自己検査は再実行で `VERDICT: PASS` でしたが、この欠陥を検出できていません。

2. **Major — `Q_up` の積分区間と温度の平均方法が、登録指標と整合していない**

   **根拠:** [survey_A.py:19](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/survey_A.py:19) の `xs < 0` による切り出し後の積分は、最後の負の節点から `x=0` までを落とします。同ファイル27行の壁温上昇は、非一様格子上の節点単純平均です。

   §4.5 の条件を既定参照格子で再計算すると、次の差が出ます。

   | 指標 | 現実装 | 積分区間・重みを修正 |
   |---|---:|---:|
   | A1 `Q_up/Q_total` | 0.139013 | 0.140533 |
   | A2 `Q_up/Q_total` | 0.375419 | 0.376580 |
   | A1 加熱区間の壁温上昇平均 | 5.47032 K | 長さ平均 5.99356 K |

   A1 の積分欠落 **0.001520** は、許容幅 `0.005` の約30%、不確かさ上限 `0.005/3` の約91%に相当します。平均温度は主判定の許容幅そのものを左右します。

   **対案:** 積分端点を含め、平均は長さ・面積重みで定義する。forge 側も対応する境界積分を使い、修正後に §4.5 と3格子の数値を再登録する。

3. **Major — A1/A2 は「同じ加熱条件」での伝導率比較になっていない**

   **根拠:** [plan:62](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:62) は同じ加熱条件を要求しますが、89–90行では `h_o` が **570→5700 W/m²K** に変わります。[survey_A.py:12](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/survey_A.py:12) が `Bi_o` を固定しているためです。固体伝導率と外面熱抵抗の効果が交絡します。

   同じ参照計算で、A2 の `h_o` を A1 と同じ約570に固定すると、`Q_up/Q_total` は **0.375419→0.332474**、節点平均壁温上昇は **9.05087→6.05867 K** に変わります。

   **対案:** **`h_o`・`T_c`・加熱区間を固定し、`k_s` だけを変える**条件に統一する。その条件で固体効果と参照不確かさを再評価することを推奨します。

4. **Major — 主参照に必要な源項・写像の検証と、不確かさの算出手順が不足している**

   **根拠:** [plan:45](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:45) の凍結流れ場による比較は、主張を限定するなら妥当です。しかし、[conjugate_ref.py:274](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/conjugate_ref.py:274) の自己検査はすべて軸対称・`v=0`・源項なしです。平面の横方向輸送、非零の散逸・圧力仕事、保存場からの写像を検査していません。

   また、同じ入力流れ場を参照格子だけ細分化しても、入力場の微分・補間に共通する誤差は検出できません。「3水準の最後の差」を不確かさ上限とする根拠も未登録です。

   前身の [`run_0014_g2_dT0_r32/TEMP_REPRODUCE.txt:5`](/home/sano/work/forge-cht/case/63.graetz_cht/run_0014_g2_dT0_r32/TEMP_REPRODUCE.txt:5) でも、温度差 **0.00144 K** に対し角の影響見積もり約 **0.0027 K**、最終記録は **`VERDICT: 判定不能`** です。小さい観測差だけでは主参照の成立を保証できません。

   **対案:** 有限の固体抵抗を持つ二材料製造解、非零 `v`・源項の試験、軸上応力極限の試験を先行させる。微分・写像方式の感度、領域延長、参照格子細分化を別々に評価し、各判定量の不確かさへ換算する手順を登録する。

5. **Major — 実行条件と合否ゲートが、まだ一意に再現できない**

   **根拠:** [plan:115](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:115) はゲート名を列挙していますが、今回の `G-if`・`G-cons`、準定常窓・閾値、熱流束誤差ノルム、規格化に使う参照側／forge側の選択が確定していません。C の上境界も69行では「slip または自由流」のままです。

   [methods/boundary.md:591](/home/sano/work/forge-cht/methods/boundary.md:591) は drift・振幅を比較許容の1/5以下と要求します。例えば前身の `drift=0.001` を絶対温度300 Kへそのまま適用すると、変動尺度は約0.3 Kとなり、今回の温度比較を支えません。また、[同:606](/home/sano/work/forge-cht/methods/boundary.md:606) の正本は `iface_q_eff` であり、参照側の片側温度差分との抽出差も管理が必要です。

   **対案:** §5.1 の run 前に「実行・評価条件の固定」を追加する。境界条件・定数物性・固体格子・IC・連成開始条件を確定し、以下を明記する。
   
   - `q_i = -iface_q_eff`、比較位置、積分重み、誤差ノルム、規格化尺度。
   - `check_convergence.py` の判定区間、`check_quasisteady.py` に渡す全節点・積分量系列と閾値。
   - `G-if`・`G-cons` の数値基準、連成更新位相を拾う保存間隔。
   - メッシュ品質確認、段階起動、格子変更時の restart、負例試験を本番 run より前に置く順序。

6. **Minor — C2 の軸方向伝導効果は、主判定と同じ量・領域で確認すべき**

   **根拠:** [survey_C.py:60](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/survey_C.py:60) の `Qax_max` は板全体の最大値です。C2 の再計算では **11.65%** の最大位置が **`x/L=0.0465`** で、評価窓 `[0.2,0.9]` の外でした。これを局所熱流束の許容3%と直接比較しても、検出能力の証明にはなりません。

   一方、メモリ上で固体の軸方向伝導だけを除くと、C2 の窓内 `θ_i` 最大差は **0.01959**。C2 を残す根拠はあります。

   **対案:** 同じ評価窓で軸方向伝導あり／なしの `T_i`・`q_i` 差を取り、その差が不確かさと主判定許容を超えることを登録する。`Q_ax` は補助指標とする。

7. **Minor — 目的・完了条件に旧方針が残っている**

   **根拠:** [plan:21](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:21) は一様熱流束加熱ですが、採用条件は有限 Robin です。§4.1 は有限差分と記載しますが、実装は節点中心有限体積です。§1 は局所 `Nu` の一致を完了条件に含める一方、§6 にその合否基準がありません。

   **対案:** 目的を有限 Robin の共役問題へ揃え、離散化名を修正する。今回は `T_i`・`q_i`・積分量を合否対象とし、`Nu` は分母が十分大きい領域での参考量と明記する。

**推奨は、現在の枠組みを維持し、参照解と事前登録を修正してから実装へ進むことです。** 優先順は **①軸上離散化・自己検査 → ②積分・平均・加熱条件の修正と再サーベイ → ③源項・写像を含む不確かさ手順 → ④実行・評価ゲートの確定**です。ソルバ変更や別ソルバ導入へ直ちに広げる必要はありません。

ファイル変更なし。上記の提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
