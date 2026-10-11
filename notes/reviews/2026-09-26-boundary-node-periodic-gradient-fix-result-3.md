# codex レビュー: boundary-node-periodic-gradient-fix (result)

- **plan**: [`plans/active/boundary-node-periodic-gradient-fix.md`](../../plans/active/boundary-node-periodic-gradient-fix.md)
- **stage**: `result` (diff base `1266aba1`)
- **date**: 2026-09-26
- **commit**: `11c77f72` (feature/sern-design)
- **codex**: effort `high`, 4.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m3
- **focus**: 3 回目 (2 回目 NO-GO C0/M4/m1 への対応確認)。§8 に固定した accepted 条件 1–6 に対する採点表 (§6.1『§8 採点』行、§6.2 #8a–#8d) の妥当性を採点してほしい。特に (1) #8a 面寄与ダンプの 3 バイナリ比較と、並進周期の基準を測定前に訂正したこと (2) #8b の定義補正 2 点 (角 group 3 格子×1、丸め床) が事後緩和でないか、H' 追試 (3) #8c 定数場 3 区分 (4) #8d の帰結 (旧の restart スパイク = F1 初回上書き) の書き方 (5) 継続課題 (dF1_inf・壁閉包 5ε・軸対称×周期スカラー未合併・ジッタ格子の LSQ 次数・R2 機構分解) の追跡先。§6.1 の運用規則: §8 の外側の新しい要求で中核の符号・適用条件・単位の誤りを示さないものは後続 plan 項目になる。accepted にしてよいか
- **extra**: `case/09.Taylor-Green/_g0_lsq_seam/R8a_facedump.txt`, `case/09.Taylor-Green/_g0_lsq_seam/R8a_R3_dumpoff.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G_channel.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G2p_fixedshape.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G2p_fixedshape_Halt.txt`, `case/39.periodic_hills/_f1split/F1_SPLIT.txt`, `notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-result-2.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
中核実装と追加試験の数値は支持します。前回の Major 4件は実質的に解消しています。  
`accepted` 移行前に、判定規則・機械判定の記録・継続課題の追跡先を同期してください。

指定の diff を HEAD `11c77f72` まで確認しました。合併 LSQ の実変位・重複配分・打ち切り、半割面除外の適用条件、`ransGradient → gather → ransBlendF1` の順序に、新たな符号・単位・条件の誤りは見つかりませんでした。

指摘は優先順に次の3件です。

1. **Minor — §8・残作業表・機械判定が、採用した規則に同期していません。**

   **根拠:** [plan:263](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:263) は、並進周期も旧新全面一致としています。しかし、測定前に訂正した [#8a:131](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:131) の正しい条件は「半割面以外は一致、半割面は新0・旧非0」です。訂正は commit `f6164c99` に残り、実装からも必然なので支持します。

   `#8b` も、結果欄では3格子の角を採用していますが、[保存 VERDICT:48](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/G2p_fixedshape_Halt.txt:48) は依然「判定保留」です。また [plan:240](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:240) の `H′zz=0.3、H′xy=0.2` は、実際には `κ=0.2` を掛けた **0.06、0.04** です。

   **対案:** §8・#8b・スクリプトを採用済みの規則へ揃え、旧判定を履歴として保持したうえで新しい集約 VERDICT を残してください。`H′` は実際の行列を記載します。R3も「全場ビット同一」ではなく、既存の勾配・リミタ判定と atomicAdd のノイズ床を区別してください。

2. **Minor — 継続課題の追跡先に、未作成・対象外の参照があります。**

   **根拠:** [plan:137](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:137) の「後続の LSQ 統一 plan」は、対応ファイルがまだありません。`dF1_inf` の参照先とする本 plan §7にも、その項目がありません。さらに、軸対称×周期の移管先である [回転周期 plan §4.5](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:81) は、軸対称を受付対象外としています。

   **対案:** 後続 plan 作成まで本 plan §5.1 #7を正式な追跡先と明記し、5課題を個別に管理してください。軸対称×周期は回転周期から分け、ジッタ次数とスカラー合併にも完了条件を付けます。解決そのものを今回の移行条件に追加する必要はありません。

3. **Minor — F1修正の適用範囲とrun索引に更新漏れがあります。**

   **根拠:** [methods/gradient.md:157](/home/sano/work/forge-sern-design/methods/gradient.md:157) は後続項目をすべて `periodicSeamMergeActive` 条件付きとしていますが、F1初期化の変更は [variables.cpp:328](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:328) で周期条件に依存しません。[case/39 README:79](/home/sano/work/forge-sern-design/case/39.periodic_hills/README.md:79) も、F1仮説を「識別中」としたままです。

   **対案:** F1修正を周期条件付きの説明から分離し、非周期SSTの初回も変わることを明記してください。run索引には `#8d` の結果を反映します。

§8の採点は次のとおりです。

| 条件 | 採点・根拠 |
|---|---|
| 1 中核 | **成立**。7変種の G0・G2・G1-a・G1-b・定数場の保存判定は PASS。`#8c` の壁∩継ぎ目は閉包参照との差 **0.91ε**。絶対閉包誤差 **5.31ε** は既知制約として残す整理を支持します。 |
| 2 固定形状 | **数値成立、記録同期が必要**。ローカルの実メッシュで面3・辺3・3格子の角3を再評価し、`H′` の全成分・両縮小対で次数 **1.0000** を再現しました。角の個数補正は妥当で、追試は丸め床除外に頼っていません。ジッタ系列の **0.855／判定不能** は別途残します。 |
| 3 同一性 | **訂正した規則で成立**。[面寄与ダンプ](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/R8a_facedump.txt) は化学種・受動種の双方を比較し、誤り版では軸対称周期面の **164/164・82/82** に差を検出しています。 |
| 4 R1・F1 | **登録した限定条件で成立**。保存収束判定は旧新とも `OVERALL: ALL PASS`。準定常ツールの再実行は旧 `ALL STEADY`、新は `dF1_inf: DRIFTING / NOT ALL STEADY`。末尾64点は **0.03334079–0.03344835**。`#8d` の修正限定パッチと **408.5→0.987倍** は、この旧restartのスパイクをF1初回上書きへ帰属する根拠として妥当です。 |
| 5 R2 | **成立**。`run_0173–0176` の保存CSVから4本の保存上限 PASSを再確認しました。機構を断定せず、観測とcommit範囲への帰属に限定する結論を支持します。 |
| 6 文書・追跡 | **上記3件の修正が必要**。 |

R1の連結収束は保存 VERDICT による確認です。AWS上の元HDF5・残差CSVや面寄与ダンプを再集計したものではありません。

**推奨は、上記3件を修正してから `accepted` へ移すことです。** 追加の大規模計算は要求しません。ファイルは変更しておらず、指摘はplan未反映です。

指摘数: Critical 0 / Major 0 / Minor 3
