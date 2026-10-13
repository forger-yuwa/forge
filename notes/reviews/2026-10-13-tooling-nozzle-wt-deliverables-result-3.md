# codex レビュー: tooling-nozzle-wt-deliverables (result)

- **plan**: [`plans/active/tooling-nozzle-wt-deliverables.md`](../../plans/active/tooling-nozzle-wt-deliverables.md)
- **stage**: `result` (diff base `origin/main`)
- **date**: 2026-10-13
- **commit**: `ddb75d6a` (feature/nozzle-deliverables)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m2
- **focus**: 2 回目の result レビュー (notes/reviews/2026-10-13-tooling-nozzle-wt-deliverables-result-2.md、NO-GO C0/M3/m1) の各指摘が最新 commit で解消したか、§6.1 の採否が妥当か、accepted にしてよいか。巨大ファイルは grep で該当箇所だけ
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**  
`ddb75d6a` で前回の Major 3件・Minor 1件は解消しています。§6.1 の採用判断も妥当です。  
移動前に、単独 CLI の不具合と現在仕様の記述漏れを修正してください。

1. **Minor — 単独の寸法ゲート CLI が、正常な入力でも合格できません。**

   根拠: [dimension_gate.py:241](/home/sano/work/forge-deliv/design/forge_design/export/dimension_gate.py:241) で `chord_pass` を必須にしましたが、[同:269](/home/sano/work/forge-deliv/design/forge_design/export/dimension_gate.py:269) の `main()` は `chord` を渡していません。

   `case/70.m5_air_wt_demo/run_0001_precfd_final/` の保存済み STEP 測定値と実メッシュを使い、CLI の判定経路を実行すると、**`未判定`・終了コード1・missing = `["chord_pass"]`** でした。同じ入力に弦検査を渡した `evaluate_gate()` は `PASS` です。一括の `deliverables` 経路には影響しません。

   **対案:** `--msh` 指定時に `chord_from_mesh()` を実行して渡し、CLI 経由の正常入力が `PASS` になる試験を追加してください。

2. **Minor — 現在仕様と配布表の説明が、修正前の判定方法のままです。**

   根拠:
   - [methods/design/overview.md:1296](/home/sano/work/forge-deliv/methods/design/overview.md:1296) は準定常判定を「16区間」と記述。
   - [deliverables.py:321](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:321) の xlsx 注記も「総熱量と16区間」。
   - [plan:296](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:296) の V5 正規化は `ΣQ` のまま。実装と §6.1 は `Σ|Q|` です。

   **対案:** 現在仕様と生成する注記を「全節点の最悪判定、16区間は補助」に統一し、V5 の式を `|Σq·ring − ΣQ| / Σ|Q|` に合わせてください。過去の変更ログは履歴として残して構いません。

前回指摘の再確認結果は以下です。

| 前回指摘 | 今回の確認 |
|---|---|
| M1 寸法ゲート | 入口端を10 mmずらす負例は `FAIL`。メッシュ測定・弦検査の欠落は `未判定`。配布 msh と STEP 由来の弦誤差 **1.189989 µm** を確認 |
| M2 積分検査の伝播 | 実 CSV の `ring_area` を1.01倍するメモリ上の負例で、相対差 **0.0100029084**、Excel part・MANIFEST・`thermal_use`・表の status がすべて `FAIL` |
| M3 局所非定常 | 同一区間内の隣接2節点で熱量を交換する負例は、総熱量・全16区間が `STEADY`、対象節点は `DRIFTING`。最新実装は全節点を判定対象にしている |
| m4 配布 README | STL の偏差設定 **0.05 mm／0.05 rad** と、CGNS の境界型を付け直す注意を確認 |

`case/70.m5_air_wt_demo/run_0006_ns_pass0_iface_60k/` の再計算では、ΣQ = **732005.554220 W**、Σq·ring = **732007.662078 W**、積分相対差 = **2.8795657×10⁻⁶**。9点の熱流束形式間の差も **1.35〜2.60%** でした。V5 の再定義の根拠と、旧基準への不合格を残した記録は妥当です。

同 run の評価は、本段 `main` の残差について **`NOT CONVERGED`**、60枚の壁ダンプについて **`DRIFTING`（総熱量21.1%/tail）**、壁解像は **`VERDICT: FAIL`（y₁⁺>1の面積3.8%）**。成果物の `thermal_use: FAIL` と [case README の run 一覧](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/README.md:11) は整合しています。未収束場を収束解として扱う誤記は、確認対象にはありませんでした。

ソルバ側の差分は `gmshReader.hpp` のタグ順検査のみで、計算カーネルの変更はありません。既知の後続事項は §5.1 に残り、`check_plans.py` は `PASS`。`run_0001`・`run_0006` の MANIFEST 掲載ファイル計42件は SHA256 が一致しました。

**推奨は、上記1→2の順に修正・確認し、レビュー記録と移動先リンクを更新して accepted に移すことです。** Fluent 実機確認と V8 の未実施は記録どおり残してください。書き込みを伴う試験一式は再実行しておらず、成功件数は追認していません。ファイル変更なし、plan 未反映です。

指摘数: Critical 0 / Major 0 / Minor 2
