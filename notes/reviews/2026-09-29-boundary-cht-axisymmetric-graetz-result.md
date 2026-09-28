# codex レビュー: boundary-cht-axisymmetric-graetz (result)

- **plan**: [`plans/active/boundary-cht-axisymmetric-graetz.md`](../../plans/active/boundary-cht-axisymmetric-graetz.md)
- **stage**: `result` (diff base `feature/cht-axisym-fem2d`)
- **date**: 2026-09-29
- **commit**: `b41c54cc` (feature/cht-axisym-graetz)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m1
- **extra**: `case/63.graetz_cht/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
`Nu` の数値偏差と延長後の準定常性は保存 CSV から再確認できましたが、登録した `V-g1` は FAIL のままです。  
追加の温度再現 A/B は時間方向・空間方向の検証条件を満たしておらず、完了を裏付ける証拠には使えません。

指定の `git diff feature/cht-axisym-fem2d...HEAD -- solver_density_cuda design methods procedures` は空でした。対象 HEAD は `b41c54cc`。ソルバ・既定値の変更による回帰は今回の差分にはありません。追加された case/63 のコードと成果物をレビューしました。

1. **Major — 登録済みの FAIL が残っており、完了条件を満たしていない。**

   **根拠:** [EVAL_SNAP.txt:13](/home/sano/work/forge-cht/case/63.graetz_cht/run_0016_g2_dT10_r32/EVAL_SNAP.txt:13) は、圧力勾配の偏差 **2.2503% > 2%**、対照の加熱開始断面の温度幅 **0.13839 K > 0.05 K**、総合 `VERDICT: FAIL` です。[plan:144](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:144) 自身も「完了は保留」としています。

   保存 CSV から再計算した `max_W |Nu/Nu_ref−1|` は、`run_0020`／`run_0016`／`run_0024` で **0.8074／0.7727／0.8127%**。この数値主張は確認できました。しかし、運動量収支で圧力勾配の差を説明できても、登録条件への適合や温度差し引きの妥当性まで証明したことにはなりません。

   **対案:** `active` に維持し、現結果を「`Nu` の数値閾値は通過、検証全体は保留」と統一する。登録閾値を維持したまま前提を満たす検証条件を事前登録し、その検証を完了条件にする。

2. **Major — 温度再現 A/B が、1枚の場だけで仮説支持を返している。**

   **根拠:** [temp_reproduce.py:183](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce.py:183) は存在するファイルを最大25枚選ぶだけで、枚数不足を拒否しません。[同:209](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce.py:209) 以降の判定も最終値だけです。

   実際の [TEMP_REPRODUCE.txt:7](/home/sano/work/forge-cht/case/63.graetz_cht/run_0014_g2_dT0_r32/TEMP_REPRODUCE.txt:7) は「末尾 **1枚**」「変動 **0**」の直後に「第1仮説を支持」と出力しています。保存 CSV に `check_quasisteady.py` を再実行した結果は、**`TRANSIENT-UNSETTLED`／`OVERALL: NOT ALL STEADY`** でした。plan §5.1 #6g の「末尾25枚・対象量の準定常判定」に適合しません。

   **対案:** 必要枚数と対象量の準定常判定を仮説判定の前提にし、不足時は `判定不能` を返す。現在の結論は「最終スナップショットで B の分布差が 0.00142 K」に限定し、時間系列を確保して再評価する。

3. **Major — 温度再現 A/B の判定位置で、離散化の不確かさを確認できていない。**

   **根拠:** [temp_reproduce.py:158](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce.py:158) は、合成試験の **`x=0` を合否から除外**しています。一方、実データの判定位置は [同:195](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce.py:195) の **`x=0`** です。

   `selftest` を再実行すると、除外地点の差は **0.0145 K**。記録された実データの角の影響見積もりも **0.0027 K** で、登録した不確かさ上限 **0.005/3 ≈ 0.00167 K** を超えます。この見積もり自体が未収束の参照に基づくため、厳密な誤差上限でもありません。下流の合成試験が PASS しても、判定地点の精度保証にはなりません。

   **対案:** `x=0` の混合境界条件を含む合成試験で参照解・後処理解を細分化し、同地点で不確かさを確認する。それまでは、分布差が 0.005 K 以下でも仮説支持の正式判定を保留する。

4. **Major — 一次データの欠落・削除により、登録ゲートの独立再検証ができない。**

   **根拠:** [README.md:61](/home/sano/work/forge-cht/case/63.graetz_cht/README.md:61) は、中間スナップショットと `conjugate_iface_log_4.csv` の削除を明記しています。ローカルの本番・延長14 run には残差原本と場のスナップショットがありません。`run_0016` に `check_convergence.py` を実行すると、**`NO residual_history.csv`／`OVERALL: CHECK FAILURES ABOVE`** でした。これは発散判定ではなく、入力不足です。

   保存された `CONVERGENCE_CHECK.txt` は14本とも `ALL PASS`、共役11本の G-if/G-cons 記録も PASS でした。しかし、全期間の残差・NaN、場から評価 CSV への抽出、界面ゲートを独立に再検証できません。指摘2の25枚検査にも直接影響しています。

   **対案:** AWS に残る圧縮残差・最終場・再開状態を検証可能な形で保存し、削除済みの必要時系列はバックアップから復元する。保存先、入力・評価器の版、ハッシュ、再実行コマンドを台帳に残す。復元できない検査は未検証として扱う。

5. **Minor — 文書と残作業表が最終実装・最新結果に同期していない。**

   **根拠:** [plan:96](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:96) は混合平均を「台形則」と記載しますが、実装は [eval_graetz.py:121](/home/sano/work/forge-cht/case/63.graetz_cht/eval_graetz.py:121) の Simpson 則です。§5.1 #5 には実装済みの `u_it` 集計が残作業として残り、#6g には最新 A/B 結果と判定不足が未反映です。[methods/boundary.md:350](/home/sano/work/forge-cht/methods/boundary.md:350) の検証範囲も親計画の円環・円板までです。

   **対案:** 残作業表を現状に更新し、指摘1〜4を優先順付きで登録する。`methods/`・索引・README は「case/63 は検証継続中」と明記し、古典 Graetz 検証済みという記述は保留する。

確認できた成果もあります。`graetz_ref.py --selftest` と3水準のメッシュ品質は **`VERDICT: PASS`**。保存 CSV の再判定では、`N_r=16/32` と延長後の `run_0023`／`run_0024` は **`ALL STEADY`**、延長前の `run_0021`／`run_0022` は **`NOT ALL STEADY`** で、plan の区別を再現しました。`CMP_Vg3_ext` の格子間差減少も再計算で **PASS** でした。

**推奨は、`accepted` への移動を見送り、一次データを確保して診断の不足を解消し、登録した検証前提を満たすまで `active` で継続することです。** ファイルは変更しておらず、本レビューの指摘は plan 未反映です。

指摘数: Critical 0 / Major 4 / Minor 1
