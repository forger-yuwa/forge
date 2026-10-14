# codex レビュー: tooling-nozzle-mknot-auto (result)

- **plan**: [`plans/active/tooling-nozzle-mknot-auto.md`](../../plans/active/tooling-nozzle-mknot-auto.md)
- **stage**: `result` (diff base `88675245`)
- **date**: 2026-10-14
- **commit**: `16c144d7` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.1 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-tooling-nozzle-mknot-auto-result-2.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M0/m2
- **focus**: 前回の result 段レビュー (notes/reviews/2026-10-14-tooling-nozzle-mknot-auto-result.md、NO-GO M2/m1) の指摘への修正 (e0055023: 境目の探し直しの上限で未確定・窓の外を N と 2N の両方で) と、その後の AWS B の再検証 (16c144d7) を確かめてほしい。対象は design/forge_design/evaluate/mknot_select.py・design/forge_design/intake/wt_card.py の M_knot の部分・試験・case/70.m5_air_wt_demo/mknot_*.py。diff にある他の plan は対象外

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

前回の Major 2件は修正され、反例の再実行でも解消を確認しました。
AWS B の保存結果と現行コードの判定は整合しています。移行前に、検証範囲と残作業表の記述を2点修正してください。

1. **Minor — M6 の「再現」の説明が、実際の検証範囲より広い**

   根拠: [methods/design/overview.md:404](/home/sano/work/forge-deliv/methods/design/overview.md:404) は既定の自動規則の説明に続けて「M6 の L_c 45 → 2.6・50 → 2.7 を再現する」と記載しています。

   実際の合格判定は、[run_mknot_select_tests.py:734](/home/sano/work/forge-deliv/design/tests/run_mknot_select_tests.py:734) の **`dM_safe=0`・`tie_deg=0`** による目的関数の回帰です。今回も **2.6／2.7** を再現しましたが、既定設定では CSV の未収録点を不成立にするため、両ケースとも未確定になります。これは実形状の不成立を意味しません。

   **対案:** 「M6 の記録に対し、安全距離・同点規則を無効化した目的関数の回帰で2.6／2.7を再現」と限定してください。既定規則全体の実データ検証は M5、という範囲を明記すれば十分です。

2. **Minor — 残作業表に修正前の現在形と旧所要時間が残る**

   根拠: [plan:128](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:128) の #3d は、修正・再検証済みの説明と「今は警告だけでゲートも PASS」「今は PASS」が混在しています。[同:129](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:129) の所要短縮も **67〜74分**のままです。

   AWS B の保存値は final **5404秒（90.1分）**、explore **5840秒（97.3分）**です。

   **対案:** #3d を「済」とし、反例の挙動は「修正前」に変更してください。#6 は最新の90〜97分へ更新し、非阻害の残課題として維持してください。

確認した実装・結果は以下のとおりです。

| 対象 | 確認結果 |
|---|---|
| 前回 M1 | 再探索上限で `MknotUndeterminedError`。旧反例のゲートは `settled=False`・`distance=False`。距離 **0.000732422** を合格させない |
| 前回 M2 | 2N の窓外境目を検出。上側の反例は **2.08→2.07**、下側は **1.92→1.93**。拡張窓の例も成功 |
| V1a | 保存評価表を現行ゲートで再判定し **PASS**。選択 **2.069270614**、距離 **0.0975**、比較記録 **37点・不一致0** |
| V6 | 保存判定は `settled` を含め全項目 **PASS**。選択 **2.07／2.07**、距離 **0.098125**。素の最小 **2.07→2.05** という説明も整合 |
| V4 final／explore | 保存成果物の検査と評価表からのゲート再判定がともに **PASS**。選択 **2.069270614**、距離 **0.09625**、寸法・YAML／card_log照合も成功 |
| V3・既存試験 | AWS B のログは数値指定2カードの差0・探索呼び出し0、`run_wt_card_tests.py --slow` は `ALL OK` |

証跡は [V1a・V6](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_select/)、[V4 final](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_card_v4_3d/check_final.txt:1)、[V4 explore](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_card_v4_3d/check_explore.txt:1)、[AWS B ログ](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_card_v4_3d/aws_logs/)です。V4 の N/2N 選択差は **0.01で合格限界**にあり、余裕が大きいとは評価できません。

指定diffを取得し、`e0055023` 以降に対象実装の変更がないことを確認しました。書込みを伴う要約試験を除いた **93チェックを再実行し、全件成功**しています。対象変更に CUDA・node境界半割面・周期処理・float32演算の変更はなく、MOC の seam 除外は参照実装と同じです。既定変更はplan・methods・手順書に明記され、`methods/index.md` の参照も有効です。

V1〜V6 は設計計算の検証で、CFD収束・準定常性を合格根拠にしていません。対象外の `run_*` の VERDICT を流用する必要はありません。今回、重い設計計算そのものは再実行していません。

**推奨:** 指摘1→2の順に文書を修正し、今回のレビュー結果を記録して `accepted` へ移してください。追加の実装修正を要求する欠陥は確認していません。

証拠: 足りる

指摘数: Critical 0 / Major 0 / Minor 2
