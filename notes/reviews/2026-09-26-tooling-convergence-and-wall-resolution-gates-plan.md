# codex レビュー: tooling-convergence-and-wall-resolution-gates (plan)

- **plan**: [`plans/active/tooling-convergence-and-wall-resolution-gates.md`](../../plans/active/tooling-convergence-and-wall-resolution-gates.md)
- **stage**: `plan`
- **date**: 2026-09-26
- **commit**: `d079a14d` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 3.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1
- **focus**: §4.6 (ソルバ側 ypls を node/mode 0 で正しく出す) と §6.2b の合格条件だけを見てほしい。実装は commit d079a14d で完了し、§5.1 12c に実測を書いた: (i) ypls 以外の 22 量が ok・ypls だけ設計どおり NG、(ii) 独立再計算との最大相対差 3.26e-5、(iv) 3.0001e-6 m、(iii) 未実施 (cell の陽性対照)。問うのは 3 点: (a) 定義 y1+ = y1*sqrt(rho_W|tau_t|)/mu_lam,W と Normal_Neighbor (cos<0.5 を未評価) が check_wall_resolution.py と本当に同義か、(b) 未評価を -1 にした設計の副作用 (後処理・ParaView・既存ツールが負値を平均に入れる等)、(c) (iii) を未実施のままにしてよいか、字面の根拠 (三項演算子の isNode==0 側が不変) で足りるか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
式自体は妥当ですが、カーネルとツールは同義ではなく、実装は対象外の mode 1/2 も上書きします。  
平板の回帰結果は支持できます。未評価値の扱いと合否手順を修正するまで、§4.6 の完了は認めません。

1. **Major — mode 0 限定が実装されていない。**

   **根拠:** [viscousFlux_d.cu:1212](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/viscousFlux_d.cu:1212) の起動条件は node・`nodeWallStressEdgeKernel`・`wall_flag` だけです。[同:1243](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/viscousFlux_d.cu:1243) で全 mode に `ypls` を渡し、[同:859](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/viscousFlux_d.cu:859) 以降で上書きします。§4.6 の「mode 1/2 はビット不変」と矛盾します。

   mode 1 の従来値は [ransWallFunction_d.cu:310](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/ransWallFunction_d.cu:310) の `utau*y/nu` です。今回の壁ノード物性・接線応力による式との同値性はありません。

   **対案:** 境界ごとの実効 `wallTreatment` が **0 の場合だけ**新診断を書き、1/2 では `ypls_b=nullptr` を渡す。既存の `twall` 補正は維持し、mode 1/2 の出力不変を固定入力の局所試験で確認してください。

2. **Major — 数式は同じでも、第一内部点と未評価集合の定義が異なる。**

   **根拠:** カーネルは [viscousFlux_d.cu:868](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/viscousFlux_d.cu:868) で内部双対面に限定し、壁ノードを除外して、内向きの符号付き射影を最大化します。一方、[check_wall_resolution.py:120](/home/sano/work/forge-cht/solver_density_cuda/tools/check_wall_resolution.py:120) は同じ DOF の境界法線を集約し、隣接点の壁属性を検査せず、**射影の絶対値**を最大化します。

   メモリ上の合成入力でツールの実関数を実行しました。直角二壁の交点で、二隣接とも別の壁上にある配置では、ツールは `y1=1e-5 m, evaluated=True` を返します。カーネルは両隣接を `wall_flag` で除外するため `ypls=-1` です。**§6.2b(ii) の未評価集合完全一致には反例があります。**

   物性異常時も異なります。[ツール:363](/home/sano/work/forge-cht/solver_density_cuda/tools/check_wall_resolution.py:363) は負密度や非正粘性をクランプしますが、カーネルは粘性だけを判定し、密度・応力の有限性を保証していません。

   **対案:** 内部点の資格、法線の単位、内向き判定、同率候補の選択、物性異常時の扱いを共通仕様にする。ツール側もその仕様に合わせ、角・他壁上の隣接・斜交・`cos≈0.5`・不正物性を含む試験で値と未評価集合を比較してください。平板1001点で未評価0件という結果だけでは、この条件を検証できません。

3. **Major — `-1` は既存の統計に混入する。§5.1 #12d の「≤0 を除外」も不適切。**

   **根拠:** [check_wall_resolution.py:377](/home/sano/work/forge-cht/solver_density_cuda/tools/check_wall_resolution.py:377) は `np.mean(sol)`、[cavity_eval.py:648](/home/sano/work/forge-cht/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:648) は無条件の面積加重平均です。例えば `[0.5,-1]` の平均は `-0.25` になり、有効値を表しません。[output.cpp:375](/home/sano/work/forge-cht/solver_density_cuda/output/output.cpp:375) も番兵を通常の数値フィールドとして出力します。

   また、正しい幾何・物性でも `|τ_t|=0` なら **`y1+=0` は有効値**です。`≤0` 除外はこれを欠損扱いします。

   **対案:** `-1` を維持するなら、利用側を同時に更新し、有限かつ `ypls>=0` を有効とする。平均の分母も有効部分に限定し、未評価面積割合を必ず併記する。ParaView 向けにも有効性フィールドと抽出手順を用意し、全点未評価は判定不能にしてください。相対差判定にはゼロ値用の絶対許容差が必要です。

4. **Major — §6.2b(i) は、記載された CLI では達成できない。**

   **根拠:** [check_field_regress.py:147](/home/sano/work/forge-cht/solver_density_cuda/tools/check_field_regress.py:147) の `--quantities` は体積量の選択です。`--boundary` は固定の `BOUNDARY_Q` を読み、境界量を追加するため、`ypls` を除外できません。

   実データで再実行した結果は **`VERDICT: FAIL (1 量が許容外)`**。`ypls` 以外の22量は許容内、非ゼロ比較の比は約0.73–1.50、`wall_4/ypls` の L2 比は約 `2.578e7` でした。非退行を支持する結果ですが、記載されたゲートの `PASS` ではありません。

   **対案:** 境界量にも適用できる明示的な除外指定を実装し、`wall_4/ypls` だけを除外した正式な `PASS` を記録する。除外した `ypls` は、別ゲートで定義・未評価集合・有限性を検査してください。

5. **Minor — (iii) は免除可能だが、「字面で合格」に置き換えてはいけない。**

   **根拠:** [procedures/verification/README.md:63](/home/sano/work/forge-cht/procedures/verification/README.md:63) は「cell の回帰対照も組まない」と明記しています。したがって、cell run の新規作成を要求する必要はありません。

   ただし、式が同じことは実行結果のビット一致の証明ではありません。さらに、変更した三項演算子は `wallTreatment==0` 内なので、**cell＋壁関数の試験は、その変更箇所を通りません**。

   **対案:** §6.2b(iii) を「現行運用方針により免除、cell は実測未検証」と明記する。代わりに指摘1の mode 限定試験を必須にしてください。「cell 不変を実測した」という主張は引き続き不可です。

実測の確認範囲は、別ワークツリー `/home/sano/work/forge/` の `case/48.flat_plate_cooled_m4/run_0061_ypls_{base_a,base_b,new_a,new_b}/`、主に `res_200.h5` と `res_wall_4_200.h5` です。`new_a` でツールと同じ保存済み `vis_lam` を使って再計算すると、最大相対差は **約1.85e-7**、逆算した距離の中央値は **3.000078e-6 m**。平板における(ii)の値比較と(iv)は支持できます。記載の `3.26e-5` 自体は再現しませんでした。

4 run の収束判定はすべて **`NOT CONVERGED`**。ここで支持するのは固定200 step の回帰・診断式の検証です。定常解や派生量の定常性を確認した結果ではありません。

**推奨は、node/mode 0 限定の修正版として進めることです。** 優先順は、①mode 制限、②ツールとの幾何・未評価仕様統一、③番兵の利用側対応、④実行可能な合否ゲート、⑤cell 試験免除の明文化です。ファイル変更禁止の依頼に従い、これらは **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
