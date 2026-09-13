# codex レビュー: condensation-kantrowitz-carrier (result)

- **plan**: [`plans/active/condensation-kantrowitz-carrier.md`](../../plans/active/condensation-kantrowitz-carrier.md)
- **stage**: `result` (diff base `feature/sern-design`)
- **date**: 2026-09-12
- **commit**: `d01a4118` (feature/condensation-air)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m2
- **focus**: 実装 diff (cond_kantrowitz_theta / CondNucCarrier の種別和 / src_jac への伝播 / sigmaScale / 診断) が plan §4 と codex plan レビューの採用事項に一致するか。§9 の主張 (回帰ノイズ床以内, 局所序列 J0>=J3>=J2>=J1, onset 序列, carrier 形が実験より 8 mm 早い解釈=CNT の J 過大の露出; ユーザ指摘の dr/dt 過大 (α=1) との縮退) が run (case/16 run_0350–0358) の実測で裏付けられるか。accepted にしてよいか
- **extra**: `case/16.nozzle_wys/compare_kantrowitz_carrier.txt`, `case/16.nozzle_wys/compare_kantrowitz_carrier_cell.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

補正式と主要な引数伝播は設計に沿っていますが、検証の合格報告に再現できない箇所があります。  
特に「低温 CNT の J 過大が露出した」という結論は、引用文献の温度依存を逆に読んでおり、現状のまま accepted へ移すべきではありません。

1. **Major — J 過大という原因判定と「実験より8 mm早い」は裏付けられていません。**

   **根拠:** [plan:194](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:194)、[methods:118](/home/sano/work/forge-cond/methods/condensation.md:118) は、Wölk–Strey を根拠に低温 CNT の過大評価を説明しています。しかし補正式は
   \[
   J_{\exp}/J_{\rm CNT}=\exp(-27.56+6500/T)
   \]
   で、230 Kでは **2.02倍**、213.8 Kでは **17.15倍**です。この式は当該温度で CNT を増幅します。低温で過小評価、高温で過大評価という説明も一次論文に明記されています。[Hale 2005、式3と本文](https://mst.elsevierpure.com/ws/portalfiles/portal/41572653/Temperature%20Dependence%20of%20Homogeneous%20Nucleation%20Rates%20for%20Water_.pdf)

   また、実測で再現できた **8.230 mm** は `run_0350` の22.513 mmと `run_0353` の14.283 mmの差です。これは**モデル間差**であり、実験の同一定義による onset 差ではありません。[実験CSV:1](/home/sano/work/forge-cond/case/16.nozzle_wys/wyslouzil_fig3_pp0.csv:1) は壁圧データで、中心線 `g=10⁻³` の位置を測っていません。壁圧偏差の「+15〜18 %」も位置誤差ではありません。

   **対案:** 「carrier 補正で計算上の onset が上流へ移動し、21 mmの壁圧偏差が増えた」までを結論にしてください。`J` と `dr/dt` の寄与は未同定です。[成長則:132](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cuh:132) の `alpha=1` を含む成長モデルとの縮退を解く前に、J 抑制較正を必須と決める根拠はありません。§5.1の成長率検証を維持し、§9・§10の断定を撤回する必要があります。

2. **Major — `verify_theta.py` が別の物性を計算しており、報告された PASS も再現しません。**

   **根拠:** [verify_theta.py:17](/home/sano/work/forge-cond/case/16.nozzle_wys/verify_theta.py:17) のエンタルピー計算は `Rw=461.5`、[C++:197](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:197) は `Ru/M=461.5222959` です。現行スクリプトの再実行結果は次のとおりです。

   - `case/16.nozzle_wys/run_0350_kw1_ref/`: 相対誤差 `9.79e-5`、絶対誤差 **`2.09e-2`、FAIL**。
   - `case/16.nozzle_wys/run_0354_kw3_sig097/`: 相対誤差 **`3.98e-4`、FAIL**。
   - `case/16.nozzle_wys/run_0357_cell_kw1/`: 絶対誤差 **`2.09e-2`、FAIL**。

   潜熱計算の気体定数だけをメモリ上で訂正すると、`run_0350` の相対誤差は **`5.96e-8`** に下がります。「float 保存」が主因という説明は誤りです。

   さらに、[同:41](/home/sano/work/forge-cond/case/16.nozzle_wys/verify_theta.py:41) は出力側の `theta>0` だけを比較対象にするため、誤って全ゼロを出力しても PASS になり得ます。[同:51](/home/sano/work/forge-cond/case/16.nozzle_wys/verify_theta.py:51) の許容値も、plan §6の `≤1e-5` と異なります。

   **対案:** 物性定数を合わせ、期待値から比較対象を選び、ゼロ・非有限値・評価対象欠落を検出してください。診断と保存場の評価時刻も確認し、絶対・相対誤差の基準を明文化して再判定してください。

3. **Major — σ感度の未達条件が、残作業表では完了扱いです。**

   **根拠:** `case/16.nozzle_wys/run_0354_kw3_sig097/` の再判定は、

   - `check_convergence.py`: **`NOT CONVERGED`**
   - `compare_condfix.py --series`: **`DRIFTING (h0err)`**
   - 末尾30000–48000 stepの `h0err` 変動幅: **0.1589 kJ/kg**、許容値0.01

   でした。onset 自体の末尾変動は小さいものの、計画した検証全体は未達です。[残作業表:121](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:121) は run 検証を完了扱いにし、[methods:119](/home/sano/work/forge-cond/methods/condensation.md:119) は未収束注記なしでσ感度を確定値として記載しています。

   [case README:351](/home/sano/work/forge-cond/case/16.nozzle_wys/README.md:351) の「実験一致の議論はσの外挿精度の範囲内」も、一定倍率の感度試験からは導けません。

   **対案:** `run_0354` の未収束・エネルギー変動の解消を§5.1へ戻してください。それまでは負側感度を過渡参考値とし、「σ換算約+10 %」や外挿精度内という結論を確定結果から外してください。

4. **Major — 単体試験の「全項目 PASS」は、実装された試験範囲を超えています。**

   **根拠:** [単体試験:66](/home/sano/work/forge-cond/solver_density_cuda/tests/unit/test_cond_kantrowitz_carrier.cu:66) の device 照合は230 Kの一点で、hostで作った `CondNucCarrier` を渡して `theta` と `J` を比較します。400状態の掃引はhost側です。計画した **float入力経路、実kernelでの種別和、ソース・差分係数、キー省略時の既定値**は検査していません。pure N2 の試験もありません。

   [同:97](/home/sano/work/forge-cond/solver_density_cuda/tests/unit/test_cond_kantrowitz_carrier.cu:97) ではmode 0を二度呼んでいますが、比較用の `Jn/rn` は検査に使われません。

   **対案:** 採用済みの検証項目を実際のテストへ追加してください。特にfloat保存量からの種別和と、本体・温度摂動・モーメント摂動のソース照合を優先します。既存テストの `ALL PASS` と、計画した全検証の達成を区別してください。

5. **Minor — cell 回帰の「ノイズ床の≤1.2倍」は過小報告です。**

   **根拠:** `case/16.nozzle_wys/run_0357_cell_kw1/` と旧 `run_0341_condfix_cell_new/` に指定の `diff_res.py --tolfile noise_cell_48000.json --factor 2` を適用すると **FAIL、9変数が超過**しました。

   - `ro`: 床の **2.36倍**
   - `roUy`: **3.02倍**
   - `Uy`: **2.71倍**
   - `h0`: **2.59倍**

   [plan:190](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:190) の数値とは異なります。またmode 3のonsetはnode **14.283 mm**、cell **14.278 mm**で、小数第2位の丸め値が同じです。cellの収束判定は **`NOT CONVERGED`**、系列判定は **`STEADY`**でした。

   **対案:** 「非ゲートの未収束準定常比較、指定の場差判定はFAIL、onsetの丸め値は同じ」と訂正してください。この結果だけでコード回帰を断定する必要もありません。

6. **Minor — plan段で採用した訂正が文書内に残っています。**

   **根拠:** [plan:47](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:47) は依然「pure N2ではmode 1と同値」、[同:57](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:57) は「混合物性から逆算」、58行は「241.8 KまでIAPWS外挿と一致」と記載しています。いずれも採用済みの設計・訂正と矛盾します。[methods/index.md:27](/home/sano/work/forge-cond/methods/index.md:27) も凝縮を「Phase 1 受動スカラー骨格」としています。

   **対案:** 実装済みの種別和、Feder純蒸気極限、methodsの修正済み文献説明へ統一し、目次の実装状態を更新してください。`plans/README.md` のactive登録は現状と整合しています。

確認できた点も明記します。指定diffでは、`a_v` 形の符号・次元、種別和、本体と `src_jac` への伝播、`sigmaScale` の核生成・成長・蒸発への適用、既定値 **0／1.0** は設計に沿っています。新しい補正式に絶対閾値によるゼロ割回避はなく、体積を掛ける既存ソース構造も維持されています。境界半割面・周期seam・軸の幾何処理への変更はありませんが、周期・軸の実測検証を今回の結果から保証することはできません。

`case/16.nozzle_wys/` の `run_0350–0353/0355` は再実行で **`PASS (converged)`／`VERDICT(series): STEADY`**、`check_quasisteady.py` の `pmax,machmax` も **`STEADY`**でした。node品質記録は **`VERDICT: PASS`**、AR最大724.2、skewness最大0.203。対象9 runの全74保存場で `VALUE/*` のNaN/Infは0です。`run_0350` のnode回帰も29変数の許容値内で再現しました。[run一覧](/home/sano/work/forge-cond/case/16.nozzle_wys/README.md:346)

**推奨は、mode 3の実装を維持してplanをactiveに留め、優先順2→4→3で検証を完了し、1・5・6の記述を訂正してからresultレビューを再度受けることです。** J較正の必要性は現段階で確定させません。

ファイル変更はありません。描画を伴う系列解析は描画初期化のみメモリ上で除外し、計量・判定処理を実行しました。CUDA単体バイナリの再実行は未実施です。本指摘は **plan未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
