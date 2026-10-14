# codex 諮問 (diagnose): tw370-result

- **brief**: [`notes/reviews/briefs/2026-10-13-tw370-result.md`](../../notes/reviews/briefs/2026-10-13-tw370-result.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-13
- **commit**: `dd753798` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **選択肢 (c) として、R0・C1 の保存済み全監視系列を既存の停止規則で再判定し、原判定を残した訂正評価と、欠落している末尾の準定常評価を作る。**
- **extra**: `case/45.isobutane_m6_d155/tw370_run.py`, `case/45.isobutane_m6_d155/tw370_judge.py`, `case/45.isobutane_m6_d155/m9_watch.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/tw370_judge.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**原判定「判定不能」は保存し、終了分類だけを訂正した再評価を別記する方針を採用する。次は追加計算ではなく、保存済み系列の再判定を行う。** Critical に相当する根拠は確認できなかった。

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **原判定を結果 B で上書き：却下。訂正評価の追加：採用** | [plan:346–347](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:346) は未到達を B、見張りの誤作動を判定不能としている。原判定と欠陥を保存したうえで、全監視系列に既存の `spacing_ok`・`done_at` を適用し、正常終了・200,000 step・全期間未到達を確認できれば、**「終了分類の訂正後：結果 B〔予算内未到達〕」を追記する**。結果を見た後の実装修正なので、その事実を明示する。閾値・予算・比較対象を変えない分類訂正と、合格条件の変更は区別できる。 |
| **Major** | **約 1 % を定常解の差／上限 50 の誤差とする：却下** | 保存された `check_convergence --segment` は **3 腕とも `OVERALL: CHECK FAILURES ABOVE`**。C50 の末尾 77,500〜97,500 step は **θ_r 全3断面 `DRIFTING`、代理 Q_w `TRANSIENT-UNSETTLED`**。[判定記録:782](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/tw370_judge.json:782)。約 1 % は異なる反復時点の出力差としてのみ記載する。R0・C1 が C50 の保存値に近づく方向に動いているという観測から、同じ漸近値へ行くとは判断できない。 |
| **Major** | **`Observer.run_forge` の戻り直後へ `RUN_RC` を移せば解決：却下** | [runner.py:279](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner.py:279) の `run_forge` は **`run_case.sh` 全体**を待つ。[run_case.sh:53](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/run_case.sh:53) は forge 終了後に収束判定・残差図を作るため、その戻り直後でも遅い。**実際の forge 終了直後に、本段・今回の起動を識別できる終了コードを原子的に記録する**。見張りはその記録を待ち、期限切れは判定不能にする。S1・S2 の終了コードを本段のものとして読ませない。単なる固定時間の待機追加では競合が残る。 |
| **Major** | **今すぐ延長／参照の作り直し：見送り。既存データの再判定：採用** | 現判定器は `EXEC_ERROR` の R0・C1 を準定常・δ* 時系列の処理から除外する。[tw370_judge.py:722](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tw370_judge.py:722)。まず分類を訂正して既存の末尾9出力を評価すれば、GPU 計算なしで欠けた判断材料を得られる。特に C50 の `δ_target` は3断面で `STEADY` でも、全分布の末尾変動幅は **最大 0.272 mm**。[記録:880](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/tw370_judge.json:880)。壁更新の可否には #33 の寸法許容との照合が必要。 |
| **Minor** | **V-g0c-chain の今回の欠陥による判定取消し：不要** | [保存記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/vgc_judge.json) は `REACHED`、120,000 step、`spacing_ok=true`、再計算 `done_at=["REACH",120000]`。この経路は `RUN_RC` を読む終端分類より先に確定する。[m9_watch.py:205](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:205)。**今回の終了分類欠陥については影響なし**。収束・準定常を保証しないという従来の限定は維持する。 |

結論: **選択肢 (c) として、R0・C1 の保存済み全監視系列を既存の停止規則で再判定し、原判定を残した訂正評価と、欠落している末尾の準定常評価を作る。**

第 1 仮説: **R0・C1 の `EXEC_ERROR` は、forge 終了から終了コード記録までの空白期間による誤分類である。** 確度: **高**。  
根拠: [m9_watch.py:209](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:209) は起動から600秒経過後、forge が不在なら `RUN_RC` 未作成を `"?"` として `EXEC_ERROR` にする。一方、記録は [tw370_run.py:1033](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tw370_run.py:1033) の `run_staged_ns` 終了後である。実コードの終端分岐だけを読み取り実行し、同じ `last=200000` でも **`rc="?" → EXEC_ERROR`、`rc="0" → CENSORED`** を確認した。  
反証条件: 原系列に途中の到達・非有限・欠損がある、または終了コードが当該本段の正常終了を示していない場合、**「分類だけの欠陥で、実態は正常な未到達」**という説明は成立しない。

第 2 仮説: **上限1の遅い過渡が、現在の腕間差に寄与している。** 確度: **中**。R0 の θ_r(40) は100,000→200,000 stepで 0.073598→0.073023 と低下している。ただし、差の全量を説明するかは未確認。  
第 3 仮説: **上限変更に依存する別の長期状態が残る。** 確度: **低・未確認**。現在の未収束系列では除外も支持もできない。

判別 A/B: **同じ保存済み系列に対し、終了コードの取得時点だけを変える。forge の追加実行は0 step。**

- A：現行どおり、当時の未作成状態 `"?"` で終端分類する。
- B：当該本段の正常終了記録を確認してから分類する。停止条件・200,000 stepの予算・2,500 step間隔・連続2回の条件は固定する。

→ **A が `EXEC_ERROR`、B が全系列の再計算でも `CENSORED`**なら、終了分類の競合を支持し、訂正評価は「結果 B：予算内未到達」とする。  
→ **B で途中の到達・非有限・系列欠損が見つかる**なら、単なる終了分類ミスという仮説を棄却し、判定不能を維持する。

同時に、R0・C1 の180,000〜200,000 stepの9出力について、登録済みの θ_r・代理 Q_w・δ_E・δ_target・出口 M の準定常判定を補う。これは追加の数値実験ではなく、欠落した結果評価の回復である。

やらない方がよいこと: **元の200,000 step試験を「止まるまで」延長して救済すること、C50 の停止場を「ほぼ収束した参照」に置き換えること、現在の低下率から到達時刻や漸近値を断定すること。** 今は (a)・(b) とも採らない。分類訂正後も延長が必要なら、元試験とは別の予算と判定を登録する。停止規則への到達だけで壁更新へ進めない。

呼び出し側の前提への異議: **「上限1では G0 と G0c の差はほぼ無い」は範囲が広すぎる。** 言えるのは「報告された200,000 stepの抽出値では差が小さい」までで、格子影響の除外ではない。また、300 K の起点も文書上は未収束の停止場であり、「ほぼ収束した場」と言い換えない。Q_w・q_w は登録上の**代理指標**である。

不足情報: 手元には次の実 run がなく、判定 JSON の監視系列も **80／80／39という行数だけ**に縮約されている。そのため、全期間未到達の独立再計算と、ブリーフの参考比較値の再抽出は未実施。

- R0：`case/45.isobutane_m6_d155/run_0811_tw370_g0_c1/`
- C1：`case/45.isobutane_m6_d155/run_0812_tw370_g0c_c1/`
- C50：`case/45.isobutane_m6_d155/run_0813_tw370_g0c_c50/`

必要なのは原 `m9_watch.json`、本段終了記録、末尾9出力・残差系列。恒久索引は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md)。

**plan 未反映。** 呼び出し側の反映先は `tooling-nozzle-isothermal-wall-chain.md` §6「V-tw370 の結果」・§5.1 #32／#33。ファイル変更・forge 起動は行っていない。
