# codex レビュー: convection-slau-wall-normal-chi-usage-rule (plan)

- **plan**: [`plans/active/convection-slau-wall-normal-chi-usage-rule.md`](../../plans/active/convection-slau-wall-normal-chi-usage-rule.md)
- **stage**: `plan`
- **date**: 2026-09-25
- **commit**: `e5ed2312` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1
- **focus**: 3 回目の plan 段。2 回目 (notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-2.md) の M2 で生産許容が 0.05 と判明し、Q1/Q2 を撤回して縮小した。縮小の理由が妥当か、§4.1 規則文・§4.2 差区間の 3 値判定・§4.3 量別閾値・§4.4 衝撃足の再同定条件が測る前に固定され全ケースを覆うか、前 plan (plans/accepted/convection-slau-wall-normal-chi.md) への『膨張角窓』訂正が適切かに集中してほしい。
- **extra**: `notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-2.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
Q1/Q2 の撤回と「膨張角窓」への訂正は妥当です。スカラーの差区間判定も前回の問題を解消しています。  
ただし、衝撃足の判定と適用診断には誤判定経路が残るため、以下を測定開始前に修正してください。

現行許容が `|ΔC_M|≤0.05` であることは [sern-3d:1918](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1918) で確認しました。既定0・2Dは0を維持し、正否や格子独立性を主張しない今回のスコープなら、2Dの格子列・CFL試験を外す判断を支持します。ただし、**末尾平均差が許容内であることと、今回の差区間ゲートに合格することは別**です。後者は未判定のままです。

保存系列を再実行した結果は以下でした。

- `case/46.sern_design/_r3_m1m3/sern2d_conv/flag{0,1}/`：本段末尾 `11999` step、両方 **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。
- `case/46.sern_design/_r3_m1m3/perrun_flag{0,1}.csv`：72点、指定閾値で両方 **`OVERALL: ALL STEADY`**。ただし抽出位置は全点 `−0.000877683 m`、flag0の窓内圧力は約 `99.57→56.60→53.55 kPa`。圧縮足の証拠ではなく、[前 plan:335](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:335) の訂正は適切です。

`ct_flag*.csv` と元の72 dumpはローカルにないため、今回、力係数の新閾値判定と新しい衝撃足抽出は実施していません。

1. **Major — 壁圧の相対 L2 に、スカラーの差区間を適用する定義がない**

   **根拠:** [対象 plan:118](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:118) は「相対 L2 は差区間で」としています。しかし、相対 L2 は既に**二つの圧力分布から作る量**であり、各 run のスカラー `C_f` を引く §4.2 の式を直接使えません。前回 M7 のこの部分は未解消です。

   人工的な定常分布 `p0=(100,102)`、`p1=(102,100)` では、各分布の L2 ノルムは同じでも、分布間の相対 L2 は **1.9801%**。各 run をノルムに縮約して差区間を作る方法では、この差を見落とします。

   [対象 plan:85](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:85) には、共通窓の固定方法、積分重み、正規化分母、「局所間隔」の選び方もありません。

   **対案:** 共通の固定窓・物理座標・積分重みを定義し、末尾の全 dump 対について直接
   \[
   \ell_{ij}=\frac{\|p_{1,j}-p_{0,i}\|_w}{\|\bar p_0\|_w}
   \]
   を計算してください。区間 `[min ℓij,max ℓij]` が全て `≤0.01` なら帯内、全て `>0.01` なら差が残る、それ以外は判定不能とします。位置差には現在の §4.2 を使えます。各 run の準定常判定は、この比較とは別に、固定座標上の圧力系列へ適用してください。

2. **Major — 衝撃足の「未検出」を「カウル衝撃が当たらない」と断定できない**

   **根拠:** [対象 plan:82](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:82) の検出条件は「最大の正勾配＋窓内圧力比1.1以上」です。これだけではカウル由来かを識別できず、弱い圧縮や広く拡散した圧縮を取り逃がします。人工圧力分布で確認すると、明瞭な8%の圧力上昇は最大正勾配を持ちますが、窓内比 `1.0800` のため不検出になります。

   また、片方の flag だけで検出、複数候補、窓がランプ端を越える場合、`m4_off` でも未検出の場合の分岐がありません。

   座標にも注意が必要です。既存 [v3sern_series.py:29](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_series.py:29) は `CELLS/centCoords` を使いますが、node の値の位置は [solver-settings.md:360](/home/sano/work/forge-sern-design/procedures/solver-settings.md:360) に従えばノード座標です。この helper をそのまま再利用すると、壁上の位置・局所間隔の定義を誤ります。

   **対案:** 未検出の結論は「登録した検出条件では同定できない」に限定してください。カウル由来との同定には、唇から続く圧縮構造との対応を確認する手順を加えること。ノード座標、`p前/p後` の採取方法、候補選択、窓不足・片側検出・両作動点未検出の処理を先に固定し、識別できない場合は「判定保留」で閉じてください。

3. **Major — `convMethod: 0` だけでは、診断ツールと実ソルバの質量収支が一致しない**

   **根拠:** [対象 plan:53](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:53) は1次区間を診断の前提としています。しかし [diag_wall_cv_budget.py:79](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:79) が計算する流束には、実カーネルの次の項がありません。

   - [convectiveFlux_slau_d.inc.cuh:620](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:620)：`slauContactFloor` による質量流束の追加。
   - [同:516](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:516)：`sstEnergyIncludesK` 有効時の `p*` 差。

   EOS整合の人工面状態で再計算したところ、`convMethod: 0`、`slauContactFloor: 0.01` では、ツールが返す流束は **`+0.09901`（流出）**、追加項を含む実カーネル式は **`−0.19092`（流入）**でした。面積を `1 m²` とした値です。条件(ii)の符号そのものが逆になります。

   **対案:** 今回は診断対象をツールと一致する実効設定に限定し、`solver`、`slauContactFloor: 0`、`sstEnergyIncludesK: 0`、前処理・面状態変更オプション等を機械的に検査してください。未対応設定は「診断不能」とすること。`ρ_i` に使う内点の選択規則も固定し、同一CV・同一面・同一dumpの変更前後を対応づける手順を §5.1 に追加してください。docstring変更だけでは適用規則を担保できません。

4. **Major — 生産許容を評価する量から `C_T_with_shear` が落ちている**

   **根拠:** [対象 plan:73](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:73) の対象は `C_T,C_L,C_M` だけです。[sern_forces.py:90](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_forces.py:90) の `C_T` は摩擦を含まず、摩擦込み推力は同ファイル98行の別列です。

   許容値の出典である [sern-3d:1867](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1867) 自体が、`C_T` の格子差は許容内でも、`C_T_with_shear` の差は許容外だった事例を記録しています。前 plan の [#10b:974](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:974) でも摩擦込み係数を評価していました。今回の縮小でこの列を外す理由はありません。

   **対案:** 既存CSVの `C_T_with_shear` を準定常・差区間判定へ追加してください。新規runは不要です。`C_T`、`C_L`、`C_M` の圧力／摩擦の定義を表に明記し、どの列の合格をもって「生産許容内」とするか固定してください。

5. **Minor — §4.3 の閾値説明と、実際に固定した数値が一致していない**

   **根拠:** [対象 plan:75](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:75) は「絶対許容の1/5」と説明しますが、`C_T≈0.92942` に対する指定値は絶対量換算で drift **`0.0001859`**、osc **`0.0004647`**。`0.002/5=0.0004` ではありません。

   また [check_quasisteady.py:282](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:282) の `--osc` は半振幅ではなく **`(max−min)/|mean|`** です。超過時は比較禁止ではなく `OSCILLATING` になります。

   **対案:** 数値を測定後に動かさず、「`C_T` は前 plan の閾値を継承、`C_L/C_M` は生産許容から換算」と説明を分けてください。各列に相対閾値、絶対換算、`OSCILLATING` の扱いを併記すれば再現可能です。

推奨は、**縮小案を維持し、実装前に「衝撃足の比較式 → 同定と未検出の分岐 → 適用診断の設定検査 → 摩擦込み推力と閾値表」の順で修正すること**です。Q1/Q2の復活やソルバ改修は不要です。`stage_manifest` の `solver` hardキー追加は妥当で、case/46中心の検証も今回の変更範囲に合っています。

ファイルは変更していません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
