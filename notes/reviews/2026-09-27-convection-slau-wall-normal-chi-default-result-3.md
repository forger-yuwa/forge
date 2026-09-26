# codex レビュー: convection-slau-wall-normal-chi-default (result)

- **plan**: [`plans/active/convection-slau-wall-normal-chi-default.md`](../../plans/active/convection-slau-wall-normal-chi-default.md)
- **stage**: `result` (diff base `22976398`)
- **date**: 2026-09-27
- **commit**: `17370e30` (feature/sern-design)
- **codex**: effort `high`, 3.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m1
- **focus**: result 3 回目。result-2 の 3 件は §5.1 #13 で対応 ( 復元 + 形状生成の回帰試験、出口角除外の再現可能な判定と全壁点 CSV、表の同期)。chi 関係の変更に絞って accepted にしてよいか判断してほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
result-2 の重大指摘2件は解消しています。chi の実装と保存データの再判定は受入を支持します。  
移行前に必要なのは、最終PASS記録への参照の同期1件です。

1. **Minor — 結果表と run 索引が、出口角除外前の FAIL 記録を参照しています。**

   [plan §6.2:194](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:194) と [case/40 README:113](/home/sano/work/forge-sern-design/case/40.nozzle_design_tool/README.md:113) は、合格の記録として `CHI_DEFAULT_B1a_isoT{,_ext}.txt` を挙げています。しかし、その延長版の結論は `VERDICT B1-a 量: FAIL` です。

   除外後の判定は別の [CHI_DEFAULT_B1a_isoT_ext_cornerexcl.txt:3](/home/sano/work/forge-sern-design/case/40.nozzle_design_tool/CHI_DEFAULT_B1a_isoT_ext_cornerexcl.txt:3) に存在します。全壁点CSVから独立に再計算し、**221点、除外1点（node 17900）、全点1.358338%、除外後0.031287%**を確認しました。証拠不足は解消しています。

   **対案:** 移行前に結果表とrun索引へ、この最終判定ファイルと全壁点CSVへの参照を追加してください。元のFAIL記録は、規則変更前の履歴として残します。

確認結果は次のとおりです。

- **実装:** 指定diffを取得して確認。三値解決・明示1の構成検査・カーネルへの `> 0` 判定は設計どおりです。chiの流束式、符号、単位、境界処理、浮動小数点演算・ゼロ割閾値への変更はありません。
- **回帰:** `_dv` 復元後の形状生成 **4/4 PASS**、manifest **17/17 PASS**、設計DB選別 **23/23 PASS**。書き込みを伴う試験用ファイルはメモリへ置き換えて実行しました。
- **残差:** 保存残差21本を `check_convergence` の判定関数で再評価し、**全本 NOT CONVERGED、NaN/Inf・RISINGなし**。全8組の `check_floor_ratio` は **PASS**、最大末尾平均比 **1.141214倍**でした。[一次判定記録](/home/sano/work/forge-sern-design/case/46.sern_design/chi_evidence/JUDGEMENTS.txt:1) と整合します。これはplanのプラトー許容による受入であり、通常基準での収束ではありません。
- **派生量:** `case/40.nozzle_design_tool/run_0958/0959_chidef_b1a_isoT_*_ext` の保存系列は、`check_quasisteady` で壁圧・熱流束・積分量とも **STEADY** を再現。SERN `run_0965`–`0968` は、m4_offの一部が **OSCILLATING**、m10_onは **STEADY**、5列すべて差区間が帯内でした。
- **文書・残作業:** 既定変更と出口角除外の事後決定は明記されています。起動との厳密な対応付け、3D格子収束、衝撃足評価の未解決事項も委譲先・制限として残っています。plan lintは **PASS** です。

**推奨は、上記1件の参照を同期し、本レビューを§6.1に記録したうえで `accepted` へ移すことです。** 今回のchi既定化について追加CFD計算を要求する根拠はありません。

ファイルは変更していません。本指摘は **plan未反映**です。

指摘数: Critical 0 / Major 0 / Minor 1
