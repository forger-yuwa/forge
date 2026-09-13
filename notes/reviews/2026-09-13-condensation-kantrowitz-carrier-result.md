# codex レビュー: condensation-kantrowitz-carrier (result)

- **plan**: [`plans/active/condensation-kantrowitz-carrier.md`](../../plans/active/condensation-kantrowitz-carrier.md)
- **stage**: `result` (diff base `feature/sern-design`)
- **date**: 2026-09-13
- **commit**: `4390f17c` (feature/condensation-air)
- **codex**: effort `high`, 5.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m1
- **focus**: 2026-09-12 の result レビュー 1 回目 (NO-GO) の採用: CNT J 過大の解釈撤回 (Wölk–Strey は低温で CNT を増幅)、verify_theta.py の R_w 修正 (mode 1 6e-8)、σ 0.97 の 96000 step 延長結果 (残差 plateau・報告量頭打ち・h0err のみ OSCILLATING)、cell 場差の記述修正、単体試験拡充。accepted にしてよいか
- **extra**: `case/16.nozzle_wys/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

Feder 補正式・種別和・`src_jac` への引数伝播に、移行を阻む実装誤りは見つかりませんでした。  
前回の物性定数修正と主要結果は再現できましたが、σ×0.97 の完了扱い、推奨設定の誤記、単体試験の達成範囲は修正が必要です。

1. **Major — 推奨設定の正本に、既定値と実験比較の誤記が残っています。**

   **根拠:** [procedures/recommended-settings.md:123](/home/sano/work/forge-cond/procedures/recommended-settings.md:123) は「既定は1」「実験より早い」と記載しています。一方、[solverConfig.cpp:689](/home/sano/work/forge-cond/solver_density_cuda/input/solverConfig.cpp:689) の省略時値は **0** です。約8.2 mmは `run_0350` と `run_0353` の**モデル間差**であり、実験との同一定義の位置差ではありません。

   また、同節の診断名 `condTheta`／`condLim` は、実際の `condTheta_0`／`condLim_0` と異なります。

   **対案:** 「グローバル既定0、Wysłouzil参照設定では1を明示」と訂正し、8.2 mmをモデル間差に限定してください。`extraFields` の例には実在する名前を使ってください。

2. **Major — σ×0.97 は準定常の参考結果として使えますが、計画した検証の完了にはなっていません。**

   **根拠:** [plan:122](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:122) は前回M3を閉じていますが、[plan:148](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:148) の合否条件は依然 `PASS`＋`STEADY` です。

   `case/16.nozzle_wys/run_0354_kw3_sig097/` の再判定は次のとおりです。

   - `check_convergence.py`: **`NOT CONVERGED (stalled/plateau)`**
   - `compare_condfix.py --series`: **`OSCILLATING (h0err)`**
   - `check_quasisteady.py --quantity pmax,machmax`: **`STEADY`**

   60000–96000 stepの4保存場では、onsetは **12.088787 mm、半範囲0.000020 mm**、壁偏差@21 mmは **21.17054 %、半範囲0.000285 %pt**。報告量の頭打ちは確認できます。一方、`h0err` は **平均1.19981 kJ/kg、半範囲0.01512 kJ/kg**です。記載された0.03は最大値−最小値であり、平均±振幅の振幅ではありません。疎な4点の非単調性だけでリミットサイクルの成立までは示せません。

   さらに、[比較表](/home/sano/work/forge-cond/case/16.nozzle_wys/compare_kantrowitz_carrier.txt:7) は旧 `res_48000.h5`、@21 mm **+21.4 %**のままで、現在のREADMEの **+21.2 %**と対応しません。

   **対案:** 本計画の受理範囲を「補正式の実装・検証」に限定し、σ×0.97は**未収束の準定常参考結果**として明示してください。エネルギー変動の未解決項目を後続のactive planへ引き継ぎ、§6・§8・§9・methods・比較成果物を同期するのが妥当です。収束済み感度としての扱いは認めません。

3. **Major — 拡充した試験でも、前回M4の「ソース・差分係数の検査」は未達です。**

   **根拠:** [test_cond_kantrowitz_carrier.cu:105](/home/sano/work/forge-cond/solver_density_cuda/tests/unit/test_cond_kantrowitz_carrier.cu:105) 以降は、モーメントをすべて0にして `S0=J` を照合しています。温度摂動も `pv`／`rho_v` を固定しています。一方、実装の [condensationSource_d.cu:238](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cu:238) は蒸気状態を再評価し、`Sg` の差分から `sj_g` を作ります。`q0>0` の `Q1` 摂動分岐と `sj_Q1` は追加試験を通りません。

   また、[plan:123](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:123) の「省略時既定をdry起動ログで確認」は成立しません。`run_0356_dry_regress/solverConfig.yaml` には凝縮セクション自体がなく、ログにも両キーの値がありません。既定値はコードでは確認できますが、読込試験の証拠ではありません。

   **対案:** 非零モーメントを用い、実装と同じ蒸気状態再評価を含む `Sg`／`SQ1` と差分係数を照合してください。凝縮セクションを有効にして両キーだけ省略する読込試験も追加し、§5.1の完了記録を更新してください。

4. **Minor — `verify_theta` の合格基準と結果記述が統一されていません。**

   **根拠:** [plan:150](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:150) は `≤1e-5`、[verify_theta.py:63](/home/sano/work/forge-cond/case/16.nozzle_wys/verify_theta.py:63) はmode 2/3で `max(1e-5, 5×S再構成誤差)` を使います。再実行した `run_0354/res_96000.h5` は相対誤差 **2.57e-4**、許容 **約1.4e-3**でPASSでした。[plan:184](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:184) の「1.2e-4以内」とも異なります。スクリプトのコメントにある絶対許容1e-4も、実判定の1e-3と不一致です。

   **対案:** 評価対象、絶対・相対許容、保存場から再構成する際の誤差を§6にまとめ、run別の実測誤差を記録してください。mode 2は式が`ln S`に依存しないため、許容の説明も訂正が必要です。

確認できた成果は次のとおりです。

| 対象（`case/16.nozzle_wys/` 配下） | 再検証結果 |
|---|---|
| `run_0350_kw1_ref`～`run_0353_kw3_feder_surf`、`run_0355_kw3_sig103` | 収束 **PASS**、系列 **STEADY**、`pmax,machmax` **STEADY** |
| `run_0350_kw1_ref` 対 `run_0335_condfix_new` | ノイズ床×2で **PASS**。29変数を判定、`gamma`は床未定義 |
| `run_0357_cell_kw1` 対 `run_0341_condfix_cell_new` | **9変数FAIL**。床比 `ro`約2.36、`roUy`約3.02、`h0`約2.59という訂正を確認 |
| 対象9 runの全74保存場 | `VALUE/*` のNaN/Inf **0** |
| nodeメッシュ品質記録 | **PASS**、AR最大724.2、skewness最大0.203 |
| `verify_theta` のmode 1 | 相対誤差 **5.96e-8、PASS** |

種別和の次元、蒸気枯渇極限、σ倍率の核生成・成長・蒸発への伝播は整合しています。確認したmode 2/3の保存場に `qhat<0` はありません。Wölk–Streyの低温での増幅方向への訂正も妥当です。[補正式の一次資料](https://mst.elsevierpure.com/files/41572653/Temperature%20Dependence%20of%20Homogeneous%20Nucleation%20Rates%20for%20Water_.pdf)

指定diffには後続の `condensation-air` の境界・EOS変更も含まれます。本判定はcarrier planに対するもので、その別計画の受理や周期・軸の実測保証を含みません。CUDA単体バイナリは実行を試みましたが、OSがGPUアクセスを拒否したためdevice検証は再現できませんでした。host側の400状態掃引はPASSでした。

**推奨は、mode 3を維持し、移行前に③試験の不足、②受理範囲と残課題の引継ぎ、①推奨設定、④判定記録の順で修正してからacceptedへ移すことです。** 成長率α感度・分圧スイープ・既定値決定も後続active planへリンクしてください。ファイル変更はなく、本レビューは **plan未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
