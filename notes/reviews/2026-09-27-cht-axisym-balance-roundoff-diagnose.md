# codex 諮問 (diagnose): cht-axisym-balance-roundoff

- **brief**: [`notes/reviews/briefs/2026-09-27-cht-axisym-balance-roundoff.md`](../../notes/reviews/briefs/2026-09-27-cht-axisym-balance-roundoff.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **date**: 2026-09-27
- **commit**: `799b58d2` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 2.0 min, rc=0
- **結論**: **閾値と今回の FAIL を維持し、失敗した交互対角 N=64 の C++ 行列で、求解精度だけを上げる A/B を一度行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 2）

| 重大度・論点 | 採否・根拠・対案 |
|---|---|
| **Major：`(e)` を `(c‴)` から外して合格にする** | **却下**。plan の適用範囲には曖昧さがあるが、コードは既に `(c‴)` の収支を判定している（[test_solid_fem2d_axisym.py:600](/home/sano/work/forge-cht/solver_density_cuda/tools/test_solid_fem2d_axisym.py:600)）。置き換え後も収支要件を継承すると明記し、今回の **FAIL を保存、V-ax1 全体は未合格**とする。 |
| **Major：「残差総和との一致」で求解丸めと確定する** | **要再検証**。FAIL の評価対象は **C++ 解と C++ の `q_hole`**（[同:552](/home/sano/work/forge-cht/solver_density_cuda/tools/test_solid_fem2d_axisym.py:552)、[同:568](/home/sano/work/forge-cht/solver_density_cuda/tools/test_solid_fem2d_axisym.py:568)）。提示された残差総和の照合は Python 求解であり、失敗した C++ 経路の原因を直接示していない。同じ C++ 行列・荷重で求解精度だけを変えて判別する。 |

結論: **閾値と今回の FAIL を維持し、失敗した交互対角 N=64 の C++ 行列で、求解精度だけを上げる A/B を一度行う。**

第 1 仮説: FP64 の線形求解誤差が集積し、全体収支の厳しい許容を超えた。 **確度: 中**
  根拠: [plan:175](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:175) は最大相対不釣り合い `2.58e-12` を記録している。登録値から入熱は **20 W/rad**、許容は **2.00e-11 W/rad**、報告された不釣り合いは **5.16e-11 W/rad**。ただし、これらは登録値と記録からの計算であり、今回の独立再測定ではない。
  反証条件: 同じ行列・荷重で高精度に評価した残差総和を十分下げても、同じ `q_hole` 評価による収支が `1e-12` を超えて残る。

第 2 仮説: 組立て後の伝導行列の零和性、または組立てと `q_hole` 評価の丸めの差が支配している。**未確認**。

判別 A/B: **交互対角 N=64・自然端面の一組だけ**を使い、C++ の組立て済み行列、荷重、境界条件、出力評価を固定する。**A = 現行解、B = 同じ行列への高精度残差による反復改良を固定 3 回**。forge は不要。補正後は FP64 の温度で既存 `field` を評価し、符号付き収支、残差総和、最大相対残差、温度補正量を見る。
  → **A の FAIL を再現し、B で残差総和と収支がともに低下、収支 ≤1e-12**なら第 1 仮説を支持する。
  → **残差総和を入熱の 1e-14 以下へ下げても収支 >1e-12**なら「求解誤差だけが原因」を棄却し、第 2 仮説へ進む。
  → 残差自体を下げられなければ判別不能とする。

やらない方がよいこと: **「収支＝残差総和」を新しい合格条件に置き換えること**。これは保存精度の要求を満たさない解でも通り得る。新しい格子を選び直す必要もない。まず同じ失敗入力で原因を切り分ける。**現時点で §5.1 #4 へ進む推奨はしない**。補正が有効なら実際の求解経路へ反映し、元の閾値で V-ax1 を再判定する。旧 FAIL は履歴として残し、修正後の PASS と区別する。

呼び出し側の前提への異議: 「剛性行和がゼロだから全桁一致」は浮動小数点では自明でない。伝導行列を \(S\)、Robin 行列を \(R\)、残差を \(r=(S+R)T-b\) とすると、整合した評価の厳密算術では
\[
Q_{\rm out}-Q_{\rm in}=\mathbf1^Tr-\mathbf1^TST。
\]
零和性の丸め誤差と出力積分の評価誤差を確認せず、求解だけに帰属させられない。また、小さい最大相対残差は、小さい残差**総和**を保証しない。

不足情報: 失敗した **C++ 解**の符号付き収支と残差総和、その正規化、および零和性の誤差。既存試験は一時ファイルを書き出すため、今回は再実行していない。**plan 未反映**。
