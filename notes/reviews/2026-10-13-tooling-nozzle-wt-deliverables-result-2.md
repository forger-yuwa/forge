# codex レビュー: tooling-nozzle-wt-deliverables (result)

- **plan**: [`plans/active/tooling-nozzle-wt-deliverables.md`](../../plans/active/tooling-nozzle-wt-deliverables.md)
- **stage**: `result` (diff base `origin/main`)
- **date**: 2026-10-13
- **commit**: `e1cd114f` (feature/nozzle-deliverables)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: 前回の result レビュー (notes/reviews/2026-10-13-tooling-nozzle-wt-deliverables-result.md、NO-GO C0/M5/m2) の各指摘が commit e1cd114f で解消したか、§6.1 の採否と §6 V5 の再定義 (結果を見てから許容を変えた記録) が妥当か。巨大ファイルは grep で該当箇所だけ。accepted にしてよいかの判定
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**  
`e1cd114f` で前回の M2・M3・M4・m6 は解消していますが、検査不合格・未実施を `PASS` にできる経路が残っています。  
V5 の再定義を支える実測値は再現できました。ただし、現状を「全件解消」として accepted に移すことには反対です。

1. **Major — 寸法ゲートは、直管長が公差外のメッシュや必須検査の欠落をまだ合格させます。**

   根拠: [dimension_gate.py:193](/home/sano/work/forge-deliv/design/forge_design/export/dimension_gate.py:193) は入口端の x 座標を検査していません。また、[同:208](/home/sano/work/forge-deliv/design/forge_design/export/dimension_gate.py:208) は「結果が存在する検査」だけを必須リストに入れます。

   メモリ上で `case/70.m5_air_wt_demo/run_0001_precfd_final/nozzle.msh` の直管部分を短縮した結果:

   - 入口端: −500 → **−490 mm**
   - メッシュの直管長区間: **189.4966〜190.7824 mm**
   - カード・STEP の直管長: **200 mm**、公差 **1 mm**
   - 構造検査: 1600×97 の規則を維持
   - 寸法ゲート: **PASS**

   別の負例でも、弦誤差が欠落・`未計算` の場合、およびメッシュ測定自体が無い場合に **PASS** でした。前回 M1 は部分修正です。

   **対案:** 入口端と直管長区間を判定に加え、出荷判定に必要な検査項目を固定してください。欠落は `未判定` または `FAIL` とし、弦誤差も固定名の `nozzle.h5` ではなく、配布対象メッシュから測るべきです。

2. **Major — 再定義した V5 の積分検査が失敗しても、総合判定は `PASS` です。**

   根拠: [deliverables.py:305](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:305) は `cfd_integral.pass` を作りますが、`part_excel()` の合否へ伝えません。[同:423](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:423) の総合判定は各 part 直下の `ok`・`pass` しか見ません。

   `run_0006` の CSV をメモリ上で読み、`ring_area_m2` だけを 1.01 倍した負例では:

   - 積分相対差: **0.0100029084**
   - 許容差: **0.0001**
   - `cfd_integral.pass`: **false**
   - MANIFEST: **PASS**、`failed_parts: []`

   **対案:** 積分検査の失敗を Excel part と総合判定へ伝播させ、この負例を回帰試験にしてください。熱設計可否にも整合検査を組み込み、ゼロ総熱量・符号相殺時の正規化方法を明示する必要があります。

3. **Major — 16 区間の熱量だけでは、局所熱流束の非定常を除外できません。**

   根拠: [wall_tables.py:565](/home/sano/work/forge-deliv/design/forge_design/export/wall_tables.py:565) は節点熱量を区間内で足し、[同:670](/home/sano/work/forge-deliv/design/forge_design/export/wall_tables.py:670) はその和だけを準定常判定へ渡します。同一区間内の変化は相殺されます。

   `run_0006` の実際の 1,600 壁節点・最終熱量を使い、同一区間内の隣接 2 節点で同量を逆向きに交換する合成系列を作りました。一方の節点は ±20% 変動します。

   - 総熱量と全16区間: **`STEADY`、変動 0%**
   - 当該節点: **`DRIFTING`、変動幅 40%**

   したがって、収束・壁解像が合格なら、この局所非定常を `thermal_use: PASS` にできる構造です。追加試験は「別区間間の交換」を検出しただけです。

   **対案:** CSV・1次元壁伝導で使う各節点の q を同じ判定区間で評価し、その最悪判定を採用してください。同一区間内で熱量を交換する負例を追加し、16区間の和は補助診断に留めることを推奨します。

4. **Minor — 配布 README に、計画で必須とした使用条件が欠けています。**

   根拠: [_readme():378](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:378) と保存済み README には、V2 の **STL の弦誤差設定値**、§4.4 の **CGNS 読み込み後に境界型を付け直す注意**がありません。後者は [fluent_msh.py:577](/home/sano/work/forge-deliv/design/forge_design/export/fluent_msh.py:577) のコードコメントにしかありません。

   **対案:** MANIFEST の STL 設定値を README に転記し、CGNS の境界型設定が必要であることを配布物に明記してください。

前回指摘の確認結果は次のとおりです。

| 前回指摘 | 今回の確認 |
|---|---|
| M1 寸法ゲート | 部分解消。上記1が残存 |
| M2 CAD 検査 | 必須条件への追加を確認 |
| M3 表への判定伝播 | `run_0006` の JSON・conditions・xlsx に `FAIL` を確認 |
| M4 公差省略 | 同じ既定値が補完され、判定が通ることを再実行で確認 |
| M5 積分検査 | 実測と再定義の説明は確認。上記2が残存 |
| m6 M5 A/B | YAML の差は `L_total_tol` のみ。`run_0001` **+0.065873 mm／PASS**、`run_0007` **−1.907017 mm／FAIL** |
| m7 文書状態 | methods・手順書・一覧の状態更新を確認。`check_plans.py` は **PASS** |

V5 の再定義について、`case/70.m5_air_wt_demo/run_0006_ns_pass0_iface_60k/deliverables/excel/wall_cfd_nodes.csv` から再計算すると、ΣQ = **732005.554220 W**、Σq·ring = **732007.662078 W**、相対差 = **2.8795657e−6** でした。`ring/A−1` の最大絶対値 **2.1119e−5** と Q 加重平均も説明に合います。**旧基準への不合格を残して再定義した記録は妥当**です。ただし、1e−4 は今回採用した工学的許容差であり、一般的な誤差上限が証明されたわけではありません。

同 run の再判定は **`NOT CONVERGED`**、熱量系列は **`DRIFTING`（21.1%/tail）**、壁解像は **`VERDICT: FAIL`（y₁⁺>1 の面積3.8%）**でした。[run 一覧](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/README.md:11)の未収束・熱設計不可という記述は適切です。未収束 snapshot を単位・符号検査に使うこと自体は問題視していません。

既存ソルバへの差分は `gmshReader.hpp` のタグ順検査で、計算カーネルの変更はありません。既知の後続事項も §5.1 に残っています。ただし、今回の問題に対応する #6・#7 の完了扱いは戻す必要があります。

**推奨は、accepted への移動を保留し、上記1→2→3→4の順で修正・負例検証を行い、result レビューを再実施することです。** 書き込みを伴う試験一式は再実行しておらず、成功件数は追認していません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
