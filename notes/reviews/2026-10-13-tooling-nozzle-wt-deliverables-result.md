# codex レビュー: tooling-nozzle-wt-deliverables (result)

- **plan**: [`plans/active/tooling-nozzle-wt-deliverables.md`](../../plans/active/tooling-nozzle-wt-deliverables.md)
- **stage**: `result` (diff base `origin/main`)
- **date**: 2026-10-13
- **commit**: `0785d720` (feature/nozzle-deliverables)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: 実装 diff と §6 の事前登録・§9 の結果の突き合わせ。巨大ファイル (wt_card.py 1400 行・wall_tables.py 860 行・fluent_msh.py 810 行) は grep で該当箇所だけ読む。methods/design/overview.md は「実寸の要求カードと成果物一式」節だけ
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**  
主要な変換経路は動作していますが、検査失敗を `PASS` にする経路と、熱設計不可の判定が配布表から欠落する不具合があります。  
対象は `origin/main...HEAD` の `0785d720`。ファイル変更なし。レビュー中に発生した未コミット変更は受入判定に含めていません。

1. **Major — 寸法ゲートがメッシュ側の不一致を合否に反映していません。**

   根拠: [dimension_gate.py:183](/home/sano/work/forge-deliv/design/forge_design/export/dimension_gate.py:183)。入口半径差・スロート半径差・直管終端の包含を記録しますが、合否に使うのは `x_end_pass` だけです。`chord` も判定に使われず、[deliverables.py:310](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:310) から渡されません。

   実測: `run_0001_precfd_final` の測定辞書をメモリ上で変更し、メッシュの入口半径を **+10 mm**、スロート径を **+20 mm**、直管終端区間を別位置にしても `evaluate_gate()` は **PASS** でした。§4.2a の出荷ゲートとして不十分です。

   対案: 配布メッシュと STEP の各寸法差を、端点の転送許容差・弦誤差・接合点の測定区間に分けて判定へ組み込み、この負例を拒否させてください。

2. **Major — 3D STEP の独立検査が失敗しても CAD 成果物が合格します。**

   根拠: [deliverables.py:107](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:107)。gmsh の読み直し失敗は `surface_gmsh_error` に記録するだけで、125 行の `surface_area_rel_err` の既定値 **0.0** により合格できます。FreeCAD の `readback_valid`、面数、bbox も合否に入りません。

   実行確認: 外部呼び出しをメモリ上で置き換え、`readback_valid=False`・面数 2・gmsh 終了コード 1 を返しても、`part_cad()` は **`pass: True`** を返しました。

   対案: V2 の各検査を明示的な必須条件とし、失敗・未実行・結果欠落を合格にしない実装へ修正してください。

3. **Major — `thermal_use: FAIL` が CSV の付属 JSON と xlsx に伝わりません。**

   根拠: [wall_tables.py:502](/home/sano/work/forge-deliv/design/forge_design/export/wall_tables.py:502) の `.status.json` に `thermal_use` がありません。[deliverables.py:270](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:270) は判定辞書を先頭 500 文字で切るため、末尾にある総合判定を落とします。

   実測: `case/70.m5_air_wt_demo/run_0006_ns_pass0_iface_60k/deliverables/` では、MANIFEST は `thermal_use: FAIL` ですが、CSV の付属 JSON は **`status: valid` のみ**。xlsx の文字列を展開しても **`FAIL`・`NOT CONVERGED` は存在しません**。§4.5 が求める、表を単独で渡した場合の使用可否表示を満たしていません。

   対案: 総合判定・対象 step・各検査の判定を専用フィールドとして JSON と xlsx に保存し、総合判定を文字数制限の対象から外してください。

4. **Major — 正式に許可した「公差省略カード」が成果物生成で失敗します。**

   根拠: [wt_card.py:254](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:254) は公差省略時に既定値を補います。一方、[dimension_gate.py:135](/home/sano/work/forge-deliv/design/forge_design/export/dimension_gate.py:135) は両公差の明示を必須とし、一括コマンドはカードを正規化せず読みます。

   実行確認: M5 カードから `tolerances` だけを除くと、寸法判定が「カードに tolerances.length_m と tolerances.diameter_m が要る」で停止しました。

   対案: 生成器と出荷ゲートで同じカード検証・既定値補完を使い、公差省略カードの一貫した往復試験を追加してください。

5. **Major — V5 の「合格」は、事前登録した積分検査を満たしていません。**

   根拠: [plan:291](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:291) は ΣQ と ∫q·2πr ds の相対差 ≤1e−9 を要求しています。しかし [run_wall_tables_tests.py:402](/home/sano/work/forge-deliv/design/tests/run_wall_tables_tests.py:402) は、`q = Q/A_node` に同じ `A_node` を掛け戻す自己整合を検査しています。

   `run_0006_ns_pass0_iface_60k` の保存済み表では次の値でした。

   | 積分 | 値 |
   |---|---:|
   | ΣQ_node | 732005.5542196009 W |
   | Σq·ring_area | 732007.6620776893 W |
   | 相対差 | **2.8796e−6** |

   事前条件の約 **2880 倍**です。これは定常熱量としての報告ではなく、同一 snapshot の積分整合検査です。`iface_Qf_eff` からの単位・符号変換を否定する結果ではありませんが、§9 の V5 全体の合格宣言は成立しません。

   対案: 保存熱量を正本として維持し、ソルバ面積と幾何面積の差を分離してください。float32 の丸めと境界面の離散化に基づく許容差を根拠付きで再定義し、自己整合とは独立した積分検査を行うべきです。

6. **Minor — V7 の事前登録 A/B が別ケースの検証に置き換わっています。**

   根拠: [plan:297](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:297) は **M5 の同じ手順**で既定 `0.05 r_t` と比較する要求です。結果の [plan:359](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:359) は **case/45 の緩和結果**を代わりに示しています。別形状で FAIL を検出したことは、指定 A/B の実施記録にはなりません。

   対案: M5 の対照を実施して両成果物のパスを記録してください。それまでは V7 の当該項目を未完として §5.1 に残すべきです。

7. **Minor — 完了扱いに必要な文書・残作業の整理が未完です。**

   根拠: [methods/design/overview.md:1244](/home/sano/work/forge-deliv/methods/design/overview.md:1244) は「計画中」、[procedures/nozzle-design-outputs.md:79](/home/sano/work/forge-deliv/procedures/nozzle-design-outputs.md:79) は「整備中」、[plans/README.md:62](/home/sano/work/forge-deliv/plans/README.md:62) は `draft` のままです。plan §7 の変更範囲には実際に変更した `gmshReader.hpp` がありません。`check_plans.py` は **PASS** でしたが、この内容整合までは検査していません。

   対案: 修正・検証後に現在仕様と状態を同期し、§5.1 の完了項目と継続項目を明示してください。#14〜#17・#12 の後続事項は引継ぎ先を確定し、今回の指摘と result レビュー記録も追加してください。V8 未実施と Fluent 未確認は、明記された扱いに従って残せます。

確認できた正常部分もあります。M5 `run_0001_precfd_final` の 20 ファイルと `run_0006_ns_pass0_iface_60k` の 22 ファイルは MANIFEST のハッシュと一致しました。M5 の拡大部長差 **+0.065873 mm**、保存済み品質判定 **`VERDICT: PASS (AR<=5000, skew<=0.90)`**、Fluent の全 **308703 面で c0 が左側**を確認しました。

`run_0006` の再判定は **`NOT CONVERGED`**、熱量系列は **`DRIFTING`（末尾 drift 21.1%）**、壁解像は **`VERDICT: FAIL`（y₁⁺>1 の面積 3.8%）**でした。全残差 CSV と最終 `VALUE/*` の非有限値はありません。[case の run 一覧](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/README.md:11)の未収束・熱設計不可という説明は妥当です。未収束であること自体を、V5 の単位検査の不合格理由にはしていません。

**推奨は、accepted への移動を保留し、上記順に修正・再検証して result レビューをやり直すことです。** 今回は書き込み禁止のため **plan 未反映**です。ファイルを生成する既存試験一式は再実行しておらず、その成功件数を追認したものではありません。

指摘数: Critical 0 / Major 5 / Minor 2
