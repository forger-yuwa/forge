# codex 諮問 (diagnose): graetz-result-disposition

- **brief**: [`notes/reviews/briefs/2026-09-29-graetz-result-disposition.md`](../../notes/reviews/briefs/2026-09-29-graetz-result-disposition.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-graetz.md`](../../plans/active/boundary-cht-axisymmetric-graetz.md)
- **date**: 2026-09-29
- **commit**: `d4ccea53` (feature/cht-axisym-graetz)
- **codex**: effort `high`, 1.9 min, rc=0
- **結論**: **旧判定の撤回と閉鎖文面・再検証限界を同期したうえで、V-g1 FAIL を保持した「ユーザ決定による例外閉鎖」として完了させる。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0／Major 3／Minor 2）

| 指摘 | 採否・根拠・対案 |
|---|---|
| **M1／Major** | **採用。例外閉鎖は可、現文面のままの閉鎖は不可。** [plan:211](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:211) の「散逸・圧力仕事で外れた」は、同じ行の「②は判定不能」と矛盾する。①の診断結果と②の未確定原因を分け、`done` は「例外受け入れによる作業終了」であり、登録検証の合格ではないと明記する。 |
| **M2／Major** | **採用、記録への反映は未完。** [判定コード:215](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce.py:215) は判定不能を返すが、[保存結果:8](/home/sano/work/forge-cht/case/63.graetz_cht/run_0014_g2_dT0_r32/TEMP_REPRODUCE.txt:8) は依然「第 1 仮説を支持」。旧記録に撤回・後継記録への参照を付け、修正版の判定結果と実行条件を保存する。 |
| **M3／Major** | **採用。ただし「48倍」の解釈を限定する。** [不確かさ評価:58–71](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce_x0_uncertainty.py:58) は、合成場の `T32−T128` を評価している。これは実場の源項・境界・抽出を含む総不確かさの保証ではない。「48倍だから物理原因が確定」は不可。「この固定場の補助計算では源項の有無で温度幅が大きく変わった」に留める。 |
| **M4／Major** | **採用、履行確認は一部要再検証。** [台帳:3](/home/sano/work/forge-cht/case/63.graetz_cht/MANIFEST_primary_data.sha256:3) は中間場・界面ログの喪失を明記する。最終場と残差の回収では、元の場から全節点の準定常性・`u_it`・温度診断系列を再抽出できない。G-if に加え、この再検証限界も記録する。181ファイルの照合と14本・3対の再実行は、許可された閲覧範囲では独立確認できていない。 |
| **m5／Minor** | **採用。** [plan:144](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:144) の限定を `methods` にも引き継ぐ。「流れのある軸対称CHTを検証済み」だけでは広すぎる。なお [不確かさ評価:6](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce_x0_uncertainty.py:6) の「N_r=16 の半径位置」は実装と不一致なので訂正する。 |

結論: **旧判定の撤回と閉鎖文面・再検証限界を同期したうえで、V-g1 FAIL を保持した「ユーザ決定による例外閉鎖」として完了させる。**

第 1 仮説: 現在の閉鎖上の問題は、訂正した判定が成果物と原因説明に反映しきれていないこと。 確度: 高  
  根拠: `TEMP_REPRODUCE.txt:8` は「支持」、`temp_reproduce.py:216–222` は「判定不能」。判定部だけを読み取り専用で実行し、1枚では枚数不足、25枚相当でも不確かさ超過により、ともに `VERDICT: 判定不能`、終了コード2を確認した。  
  反証条件: 保存結果に正式な撤回・後継判定が紐付き、閉鎖文面から②の原因断定が除去されていること。

第 2・第 3 仮説: 追加しない。温度幅の物理原因は今回の証拠では確定しない。

判別 A/B: **追加の物理計算は不要。** 同一の比較数値・登録枚数25・不確かさを固定し、判定部へ渡す枚数だけ A=1／B=25 に変更する、0 step の確認を実施済み。Aだけ判定不能なら枚数不足が唯一の障害、Bも判定不能なら不確かさも独立した障害。**実測は後者**であり、枚数を揃えるだけでは診断は成立しない。

やらない方がよいこと: 許容差を変更する、旧「支持」を有効な判定として残す、例外閉鎖を検証合格へ読み替える、FP64の結果をFP32へ一般化する。

呼び出し側の前提への異議: ユーザの例外受け入れは閉鎖の根拠になるが、原因の証明にはならない。「0.82 %以内」は差し引きNuの数値偏差として記載し、絶対熱流束の精度やV-g1成立の保証にしない。

不足情報: 閉鎖記録には、再実行した14本・3対の対象パス、コマンド・コード版、判定区間、保存VERDICTへの参照と、ハッシュ照合結果を紐付けること。今回それらの独立照合は行っていない。plan未反映（依頼どおり変更なし）。

閉鎖文案（plan §1・`methods`）:

> case/63 の登録した node・FP64・M=0.05 条件では、評価窓 \(x^+\in[0.003,0.1]\) における対照差し引き局所Nuの古典Graetz基準からの数値偏差は、保存評価記録上、3格子で最大0.82 %以内だった。前提検査V-g1のFAILは保持し、2026-09-29のユーザ決定により例外として閉鎖した。温度非一様の原因診断は判定不能であり、絶対熱流束の精度は保証しない。中間場・界面ログの喪失により、時系列・界面判定の一部は保存済み記録に依存する。
