# codex レビュー: tooling-nozzle-mknot-auto (result)

- **plan**: [`plans/active/tooling-nozzle-mknot-auto.md`](../../plans/active/tooling-nozzle-mknot-auto.md)
- **stage**: `result` (diff base `88675245`)
- **date**: 2026-10-14
- **commit**: `f8e77f8c` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 5.4 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-tooling-nozzle-mknot-auto-result.codex.log` (codex: ok)
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **focus**: この plan の範囲は design/forge_design/evaluate/mknot_select.py・design/forge_design/intake/wt_card.py の M_knot の部分・design/tests/run_mknot_select_tests.py・run_wt_card_tests.py の M_knot の部分・case/70.m5_air_wt_demo/mknot_*.py。diff にある他の plan (mesh-size-guard・review-fallback-backend・停止規則) は対象外

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

M5 の保存結果は報告と整合し、V4 final/explore の寸法・記録照合も再実行で PASS でした。
しかし、安全距離 0.05 未満の候補を解像度ゲートが PASS にする反例を2件再現しました。現状での `accepted` 移行には反対します。

1. **Major — 再探索の上限で、既に発見した境目を捨てて合格させる**

   根拠: [mknot_select.py:518](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:518)、[同:588](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:588)。3回目にも境目が増えた場合、最新の `nb`／`nbN`／`nb2` を採用せず、古い境目による選択を警告付きで返します。その警告はゲートの不合格条件に入りません。

   **現行コードでの再現:** 成立範囲を M ≥ 2.0037、ただし不成立区間を `(2.053, 2.055)`、`(2.105, 2.107)`、`(2.157, 2.159)`、目的関数を θ = 15 + 10(M − 2.02)² としました。既定の `select_mknot` は **2.16、距離 0.052929688** を返します。しかし、その戻り値の評価表から既存の `find_boundaries` で求め直すと、既に評価済みの境目 **2.159267578** までの距離は **0.000732422** です。同じ評価関数による `resolution_gate` も **PASS、全 checks=True** でした。

   対案: 最新の境目で再選択が確定するまで合格させず、反復上限では `MknotUndeterminedError` にする。この反例を `select_mknot` と `resolution_gate` 両方の回帰試験に追加してください。

2. **Major — 2N にだけ存在する窓外の境目を探索しない**

   根拠: [mknot_select.py:688](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:688)。窓外の追加評価は N だけで、2N は窓内と「N で見つかった境目の両端」しか調べません。[wt_card.py:1928](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1928) の最終評価も N だけを ±0.15 に広げています。

   **現行コードでの再現:** 入力 2.00、窓 `[1.90, 2.10]`、θ = 15 + 100(M − 2.08)²、成立上限を N では 2.30、2N では 2.12 としました。カードと同じ N の ±0.15 の評価表を渡すと、**両選択 2.08、境目なし、距離 Infinity、PASS**。2N の窓外も評価すると境目は `[2.12, 2.120625]`、成立側端までの距離は **0.04 < 0.05** でした。

   対案: 初回・拡張後とも、**N と 2N の両方**で候補窓の外へ ΔM_safe 分を評価し、それぞれ独立に境目を求めてください。これは実装修正に加え、plan §4.4 の「N だけを広げる」方針の修正も必要です。上側・下側双方の反例を追加してください。

3. **Minor — 残作業表に完了済み作業が未了として残っている**

   根拠: [plan:119](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:119) の #1 は「次は plan 段のレビュー」、[同:123](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:123) の #3 は「V4 explore は未了」のままですが、§6.1 と #3c には完了記録があります。

   対案: #1・#3 を現在の状態へ更新し、今回の2件を未解決項目として追加してください。所要短縮は #6 に残っており、記録漏れではありません。

検証結果については、次を確認しました。

| 対象 | 確認結果 |
| --- | --- |
| V1a | 保存 JSON は選択 **2.069270614**、θ_max **15.288322993°**、距離 **0.0975**、比較 **37/37・不一致0**、ゲート PASS |
| V1b・V2 | ファイルを書かない V2(a)〜(l) と V1b を再実行し全件成功。M6 の選択は **2.6／2.7** |
| V3 | 保存ログは数値指定2カードとも **PASS**、YAML・card_log の差0、探索呼び出し0 |
| V4 final | 保存成果物の検査を再実行して **PASS**。拡大部長差 **+0.3680 mm**、物理径差 **−3.531×10⁻⁸ m**、所要 **4449 s** |
| V4 explore | 同じく **PASS**。拡大部長差 **+0.3691 mm**、物理径差 **−2.358×10⁻⁸ m**、所要 **4007 s** |
| V6 | 保存 JSON は選択 **2.07／2.07**、距離 **0.098125**。素の最小 **2.07→2.05** と同点規則による合格、という説明も整合 |

根拠の保存先は [V1a・V6](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_select/)、[V3 ログ](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_card_v4_3b/aws_logs/b2_v3.log)、[V4 final](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_card_v4_3b/check_final.txt)、[V4 explore](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_card_v4_3b/check_explore.txt) です。V4 の N/2N 差は **0.01** で、合格限界にあります。

指定 diff を取得し、AWS 検証時の `b2528f7f` から HEAD まで対象実装に変更がないことも確認しました。既定の自動化は plan に明記され、`methods`・手順書・skill に反映されています。`methods/index.md` の参照も有効です。対象変更には CUDA の境界半割面・周期処理・float32 演算の変更はなく、MOC の seam 除外は参照実装と同じでした。

今回の V1〜V6 は設計計算の検証で、CFD の収束・準定常性を合格根拠としていません。したがって、対象外の `run_*` の収束 VERDICT をこれらの PASS の裏付けには使っていません。

**推奨は、移行を保留し、指摘1→2→3の順に修正することです。** 2件の反例を回帰試験に加え、修正後の境目集合で V1a・V6・V4 final/explore を再検証してから result レビューを更新してください。

証拠: 足りる

指摘数: Critical 0 / Major 2 / Minor 1
