# codex レビュー: boundary-cht-axisymmetric-fem2d (result)

- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **stage**: `result` (diff base `488975be`)
- **date**: 2026-09-27
- **commit**: `20d1e27b` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 7.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m1
- **focus**: 軸対称 fem2d の V-ax0〜4 が主張の範囲 (§5.1 #12) どおりに裏付けられているか、done・拒否解除 (受理条件のコードは既に入っている) してよいか。流体側の未解決 (一様 IC で非一様格子が偽流れに居座る) は別 plan
- **extra**: `case/61.conjugate_annulus/README.md`, `case/62.conjugate_disk/README.md`, `case/62.conjugate_disk/vax0_evidence/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
§5.1 #12 の限定された連成精度・保存・定常性の主張は、実装と実測が支えています。  
ただし、再開契約の検査漏れと文書の不整合を直してから `accepted` へ移すべきです。

1. **Major — 新形式の再開状態を、属性の部分欠落によって旧形式として受理できる。**

   [conjugateWall.cpp:611](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:611) は、`geometry` が無ければ平面側で契約検査を通過させます。`load_unit`・`state_contract` の検査は `geometry` がある場合だけです。[check_cht_balance.py:96](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:96) も同じ穴を持っています。

   `case/61.conjugate_annulus/run_0005_annulus_r32s16_w20k/conjugate_state_4.h5` の属性をメモリ上で複製し、`geometry` だけを除いた負例で確認しました。`load_unit=W/rad`・`state_contract=2` が残っているのに、平面 run に対する `_state_geometry()` は **`(False, None)`＝拒否理由なし**を返します。C++ 側も同じ条件で通過するコードです。荷重履歴がある状態では、異なる単位の `SOLID/QBUF` を復元し得ます。

   **対案:** 旧形式として受理する条件を「識別属性 **3 個すべてが無い、かつ平面**」に限定する。部分欠落は C++・評価器とも拒否し、各属性の単独欠落と旧平面状態の互換性を V-ax4 の負例へ追加してください。今回の連続・再開温度比較の PASS は、この検査漏れを検出しません。

2. **Minor — 現行仕様・索引・ケース手順に、実装前の拒否条件と未決事項が残っている。**

   [methods/boundary.md:559](/home/sano/work/forge-cht/methods/boundary.md:559) は「軸対称の面内伝導」を起動時拒否対象としており、同文書の新しい受理条件と矛盾します。[plans/README.md:28](/home/sano/work/forge-cht/plans/README.md:28) も「全モードで拒否済み」「上位に諮る前」のままです。

   [case/61 README:79](/home/sano/work/forge-cht/case/61.conjugate_annulus/README.md:79) は `warmup: 200` を案内し、[同:103](/home/sano/work/forge-cht/case/61.conjugate_annulus/README.md:103) は格子感度を「保留」としています。実際の合格構成は `warmup: 20000` で、感度系列も実行済みです。

   **対案:** 現行の受理条件と合格 run の再現手順へ同期し、旧記述は日付付き履歴として区別してください。`accepted` 移動時には索引・リンク・§6.1 のレビュー採否も更新する必要があります。

実測の再確認結果は次のとおりです。指定 diff を取得して確認し、固体の剛性・Robin 行列・集中量・`q_hole` の式、積分済み荷重の受け渡し、流体面積の半径重みと床の対応に、新たな符号・単位の誤りは見つかりませんでした。

| 対象 | 確認結果 |
|---|---|
| V-ax0 | 軸上節点の拒否、method 1 の NaN 化・属性出力はコードと[保存記録](/home/sano/work/forge-cht/case/62.conjugate_disk/vax0_evidence/README.md)で確認。GPU 負例の再実行はしていません。 |
| V-ax1 | 置換試験の Python 再計算で交互対角の次数 **1.968 / 1.984**、同方向対角の最細誤差 **0.0107868 K**を再現。既存円環・円板固体で C++–Python 求解差は最大 **1.04e-11 K**。収支 **2.58e-12 > 1e-12** は引き続き例外受け入れであり、PASS への変更は不可。丸め A/B の原データは今回再検証していません。 |
| V-ax2 | `case/61.conjugate_annulus/run_0005_annulus_r32s16_w20k/`：壁温誤差 **0.037526%**、G-cons **1.532e-6%**を HDF5 から再計算。残差 **PASS**、準定常 **ALL STEADY**、保存済み G-if **PASS**。 |
| V-ax2b | `case/62.conjugate_disk/run_0013_disk_r32u_w20k/` と `run_0019_disk_r32g_condic/`：誤差 **0.047650% / 0.085774%**、系列比 **0.3174 / 0.2978**を再現。双方、残差 **PASS**・準定常 **ALL STEADY**・保存済み G-if **PASS**。非一様側は純伝導 IC 限定。 |
| V-ax3 | `case/58.conjugate_slot/run_0015_vax3_{new_a,new_b,base}/`：新旧差 **1.181460e-3 K ≤ 3.964181e-3 K**を再現。有限反復の回帰判定は PASS。保存済み残差判定は **NOT CONVERGED**で、物理解の収束証明には使えません。 |
| V-ax4 | `case/61.conjugate_annulus/run_0012_vax4_cont_{a,b}/`・`run_0013_vax4_resume/`：再開差 **1.897706e-8 K ≤ 2.970046e-8 K**を再現。登録比較は PASS。再開 run 単独の残差判定は **NOT CONVERGED**。 |

円環・円板の本体／感度用 **9 run** は、全残差履歴の再判定で **PASS**、全界面節点の温度降下・熱流束系列の再判定で **ALL STEADY**でした。判定条件は登録どおり `--tail 0.5 --drift 0.001 --osc 0.001`。メッシュ品質も全件 **VERDICT: PASS**。読み出した場の主要量に非有限値、密度・温度・圧力の非正値はありませんでした。

円環の最後の格子間差 **0.006550 K**、流体半径感度 **0.001757 K**も再現しました。粗格子の円環 `run_0009` の熱流束誤差 **0.763577%**、円板 `run_0011` の壁温誤差 **0.504119%**は、個別評価器に FAIL として残っており、総括による上書きはありません。成果物の索引は [case/61 の run 一覧](/home/sano/work/forge-cht/case/61.conjugate_annulus/README.md:115)と [case/62 README](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md)です。

**推奨は一つです。** 指摘1の契約修正・負例確認を最優先に行い、次に指摘2の文書を同期したうえで、§5.1 #12 の限定と例外を維持して `accepted` へ移してください。一様 IC の未解決事項は[流体側 plan §5.1 #4](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:79)に残っており、本 plan の完了条件へ戻す必要はありません。

ファイル変更なし。本レビューの指摘・推奨は **plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 1
