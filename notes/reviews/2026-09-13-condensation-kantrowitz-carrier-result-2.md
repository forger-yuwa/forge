# codex レビュー: condensation-kantrowitz-carrier (result)

- **plan**: [`plans/active/condensation-kantrowitz-carrier.md`](../../plans/active/condensation-kantrowitz-carrier.md)
- **stage**: `result` (diff base `feature/sern-design`)
- **date**: 2026-09-13
- **commit**: `79131b2d` (feature/condensation-air)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m3
- **focus**: 2026-09-13 result レビュー 2 回目 (GO-with-changes, M1 推奨設定の誤記 / M2 σ0.97 の完了扱い / M3 差分係数試験と読込試験 / m4 verify_theta 許容) の採用 (§6.1, §9 2026-09-13 ②): test_cond_kantrowitz_carrier (b3)、run_0034 読込試験、σ0.97 の参考値化と受理範囲 (§8)、後続の condensation-followups.md への引き継ぎ。accepted にしてよいか
- **extra**: `case/16.nozzle_wys/README.md`, `plans/active/condensation-followups.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

Feder 補正式・種別和・`src_jac` への伝播に、受理を阻む実装誤りは見つかりませんでした。
前回の主要修正は成立していますが、比較成果物の誤数値と、試験・引き継ぎ・現在仕様の記述を移行前に直してください。

1. **Major — σ×0.97 の更新成果物に、誤った比較値と無条件の感度記述が残っています。**

   **根拠:** [compare_kantrowitz_carrier.txt:8](/home/sano/work/forge-cond/case/16.nozzle_wys/compare_kantrowitz_carrier.txt:8) の最終列は「先頭 run 比」なのに、`run_0354/res_96000.h5` の場差がすべて **0** です。`run_0350/res_48000.h5` を基準に再計算した `ro/P/T/Ux/g_0` の差は、それぞれ **0.114 / 0.121 / 0.113 / 0.184 / 0.671** でした。

   また、[recommended-settings.md:125](/home/sano/work/forge-cond/procedures/recommended-settings.md:125) は依然「±3 %で onset ∓2.3 mm」と無条件に記載しています。実際の `case/16.nozzle_wys/run_0354_kw3_sig097/` は、収束判定 **`NOT CONVERGED (stalled/plateau)`**、系列判定 **`OSCILLATING (h0err)`** です。

   **対案:** 比較表を全対象 run・同一基準で再生成してください。推奨設定と単独配布される比較表にも、「σ×1.03：+2.5 mm、σ×0.97：−2.2 mmは未収束参考値」と明記してください。plan §8 の受理範囲限定自体は妥当です。

2. **Minor — `(b3)` のモーメント入力が、記載した単分散状態と単位不整合です。**

   **根拠:** [test_cond_kantrowitz_carrier.cu:127](/home/sano/work/forge-cond/solver_density_cuda/tests/unit/test_cond_kantrowitz_carrier.cu:127) は `g/(4πρ_l r³/3)`、つまり質量当たりの `Q0` を作り、そのまま渡しています。一方、[condensationSource_d.cuh:190](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cuh:190) の入力は体積当たりの保存量 `ρQn` です。

   非零成長項と差分経路の検査にはなっていますが、「g=0.002、半径20 nmの単分散状態」の再現にはなっていません。

   **対案:** `q0=rod*g/(4πρ_l rbar³/3)` として `q1/q2` も揃え、`(b3)` を再実行してください。温度差分の誤差と成長項の優勢度も更新してください。

3. **Minor — 後続計画への引き継ぎで、原因判別に必要な観測量が抜けています。**

   **根拠:** [親plan:129](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:129) は、onsetだけでは `J` と `dr/dt` を区別できないため、圧力上昇幅・出口液滴径・SAXS文献取得を要求しています。しかし、残作業の正本とした [condensation-followups.md:43](/home/sano/work/forge-cond/plans/active/condensation-followups.md:43) は onset・壁圧応答だけです。同表45行の既定値候補も「1と2/3」で、現行既定0の維持が明示されていません。

   **対案:** #1 に圧力上昇幅、液滴径、SAXS資料取得を引き継ぎ、#3 は「現行0を維持するか、変更するか」としてください。Wysłouzil参照設定の1とグローバル既定を区別してください。

4. **Minor — 現在仕様と完了記録に、訂正前の説明が残っています。**

   **根拠:** [methods/condensation.md:93](/home/sano/work/forge-cond/methods/condensation.md:93) は Feder carrier 拡張を「未実装」と記載し、直後の実装済み説明と矛盾します。[plan:123](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:123) にも、前回否定された「run_0356 dry の起動ログで既定値確認」が残っています。

   **対案:** 前者を実装済み節への参照に置換し、後者を `case/34.arthur_n2_nozzle/run_0034_cfg_default_check/forge_run.log` に訂正してください。同ログの実効値 **`condKantrowitz=0 condSigmaScale=1`** は確認できました。

再検証した主要証拠は次のとおりです。run索引は [case/16 README](/home/sano/work/forge-cond/case/16.nozzle_wys/README.md) にあります。

| 対象（`case/16.nozzle_wys/` 配下） | 再検証結果 |
|---|---|
| `run_0350_kw1_ref`～`run_0353_kw3_feder_surf`、`run_0355_kw3_sig103` | `check_convergence`: **PASS**、報告量系列: **STEADY**、`pmax,machmax`: **STEADY** |
| `run_0354_kw3_sig097` | **NOT CONVERGED**、系列 **OSCILLATING (h0err)**。末尾4場の `h0err` 全幅 **0.03024 kJ/kg** |
| `run_0350` 対 `run_0335` | `diff_res --factor 2`: **PASS**。29変数判定、`gamma` は床未定義 |
| `run_0357` 対 `run_0341` | **9変数FAIL**。cellは引き続き非ゲートの未収束参考結果 |
| 対象9 run・全74保存場 | `VALUE/*` と全CSV残差のNaN/Inf **0**。nodeメッシュ品質記録 **PASS** |

`verify_theta.py` は対象全 run で **PASS**、記載された誤差を再現しました。既存単体バイナリの host 側では400状態掃引と `(b3)` の差分検査が通りました。ただしCUDA初期化が失敗したため、device側の **ALL PASSは今回再確認できていません**。系列検査は、キャッシュ書込みを要求する描画初期化を除き、元の集計・判定コードをメモリ上で実行しました。

**推奨は、上記1→2→3→4の順で修正し、§8の限定された受理範囲でacceptedへ移すことです。** σ×0.97の収束達成や、同じdiffに含まれる別の `condensation-air` 計画、周期・軸境界の実測保証まで本判定に含めないでください。ファイル変更なし、本レビューは **plan未反映**です。

指摘数: Critical 0 / Major 1 / Minor 3
