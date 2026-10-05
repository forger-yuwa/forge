# codex 諮問 (diagnose): sern-r7a-g4

- **brief**: [`notes/reviews/briefs/2026-10-05-sern-r7a-g4.md`](../../notes/reviews/briefs/2026-10-05-sern-r7a-g4.md)
- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **date**: 2026-10-05
- **commit**: `30a7e7a4` (feature/sern-design)
- **codex**: effort `high`, 3.2 min, rc=0
- **結論**: **旧 g4 判定を「判定不能」として保存し、B_g4 の同一設定＋20000 step継続を、固定長の追加診断として事前登録する。**
- **extra**: `notes/investigations/2026-10-05-sern-r7a/R7A_G4_VERDICT.txt`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 採否 | 根拠と対案 |
|---|---|---|
| **Major** | **(a) 追加延長を別試験として事前登録する案を採用** | 旧 g4 試験は **`VERDICT: 判定不能` のまま確定・保存**する。`run_1048` の `C_M` 前窓差は7.13×10⁻⁴で、登録閾値5×10⁻⁴を超える。追加試験は停止時点・判定窓・必要条件を先に固定し、旧試験の合否を書き換えない。[判定原本:2](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R7A_G4_VERDICT.txt:2) |
| **Major** | **(b) 今回の救済策として `C_M` を外す案を却下** | 現行評価器は4量を必要条件とし、`C_M` は最適化の制約にも使う。推力だけの将来試験を別目的で定義することは可能だが、それでは現行評価器の成立を認定できない。今回は4量の条件を維持する。[受理仕様:878](/home/sano/work/forge-sern-design/methods/design/overview.md:878) |
| **Major** | **g3 の限定認定は採用、「格子依存は10⁻⁷級」は却下** | g3 の保存判定は **`分岐 A: 二水準を識別できる`**。一方、g3 の E＝6.89×10⁻⁶とg4 の U＝6.32×10⁻⁶を同じ幅の扱いで足すと、差の差の比較幅は **1.321×10⁻⁵**。中心値−1.285×10⁻⁷はその約1/100であり、分解できた格子感度とは言えない。対案は「g3で登録基準を満たした。g4の参考中心値も約＋0.002086だが、全量判定は未成立」と記録する。[g3原本:8](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R7A_VERDICT.txt:8)、[g4原本:4](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R7A_G4_VERDICT.txt:4) |
| **Major** | **「初期写像に依存しない」は限定表現へ修正** | 写像対照は推力2量・`C_L`について登録閾値内だが、全体は **`VERDICT: 判定不能`**。さらにg3 Bを写像対照に置き換えるだけで、推力の差の差の中心値は約−1.3×10⁻⁷から＋3.0×10⁻⁷へ変わる。「今回の写像変更による差は指定許容内」までとする。[写像対照:3](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R7A_MAP_VERDICT.txt:3) |
| **Major** | **実機トリム認定を保留し、既存の暫定制約を維持** | コードの `C_M` はノズル面の圧力モーメントで、摩擦・機体面を含まない。ユーザー決定の加重平均 `cm_min: -7.0` は探索用条件として維持できるが、作動点別の実機トリム保証にはならない。定義を明記し、実機窓を捏造しない。[集計:309](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:309)、[決定履歴:437](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:437) |

結論: **旧 g4 判定を「判定不能」として保存し、B_g4 の同一設定＋20000 step継続を、固定長の追加診断として事前登録する。**

第 1 仮説: **B_g4 は緩慢な過渡が残っており、さらに20000 step進めれば、現行定義の `C_M` も登録窓条件を満たす。** 確度: **中・未確認**。  
根拠: ブリーフの前窓差は9.75×10⁻³→7.13×10⁻⁴と縮小し、後者は[判定原本:2](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R7A_G4_VERDICT.txt:2)でも確認できる。ただし、生履歴がないため、減衰形状や漸近値は未検証。  
反証条件: 追加20000 step後も、同じ前後10000 stepの窓条件を満たさないこと。この場合、「もう20000 stepで足りる」という仮説を棄却する。過渡全般を否定したことにはしない。

第 2・第 3 仮説: **追加しない。** 現資料だけでは、振動・写像依存の残存・別の数値機構を順位付けできない。

判別 A/B:

- **A＝既存の累積40000 step結果**：`case/46.sern_design/run_1048_r7a_g4_lsw10_c40k/`。
- **B＝累積60000 step結果**：Aの最終場を、別の新規 `run_NNNN_*` に `restart_field.py` で継承する。変更は反復数だけ。格子・実効設定・バイナリ・出力間隔500 stepを固定する。
- **固定長**：追加20000 step。途中の閾値通過で止めず、追加区間の前半10000／後半10000 stepで判定する。
- **必要条件**：推力2量・`C_L` の前窓差≤5×10⁻⁵、`C_M`≤5×10⁻⁴、4量の `check_quasisteady` が `STEADY`、有限値・床・品質ゲートを満たす。全残差の `check_convergence` VERDICTと判定区間も保存する。
- **結果A：必要条件を満たす** → 追加診断では反復不足説を支持。既存A_g4との推力差とUを再計算し、`|Δ|−U>0.002`なら、**追加試験として**g4の識別成立を記録する。
- **結果B：必要条件を満たさない** → 固定長延長で解消する仮説を棄却し、判定不能を維持する。自動的な再延長や閾値変更は行わない。

やらない方がよいこと: **通るまで延長を繰り返すこと、`C_M`の基準点を変えて窓条件を通すこと、g4参考値の小数点一致を精度保証にすること。** `GATES PASS`と今回の厳しい絶対窓条件は別なので、両者の食い違いだけで判定器の不具合とは言えない。

呼び出し側の前提への異議:

- **`C_L`合格・`C_M`未達は矛盾しない。** 現行基準点はx_ref＝−20H。コードのモーメント式では、`C_L`変化3.57×10⁻⁵だけでも、基準点移動に対応する項は約7.13×10⁻⁴になる。これは長い腕による感度の説明であり、今回の原因を確定するものではない。[モーメント式:259](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:259)
- **R7bへの順序**は、①上記追加診断、②現行`C_M`の面集合・基準点・正規化・摩擦除外を仕様として明記、③`L_sw`を物理stationとして指定・再現できるようにする、④その形状で他作動点の評価成立と設計差の格子・領域感度を確認、⑤3D最適化へ接続、とする。現在は最寄りstationへ丸める実装なので、連続探索前に③が必要。[station選択:148](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:148)
- 実機窓の出典がないことは、暫定制約による探索準備を止める理由にはしない。ただし、その探索結果を実機トリム成立とは表現しない。

不足情報: **対象runの生履歴・個別VERDICT・実効config・起動記録は、このcheckoutにない。** したがって、残差・準定常性・有限値・品質の独立再判定はできていない。追加試験前に、`force_history.csv`、`residual_history.csv`、個別判定、restart来歴、実行した集約スクリプトを保存・照合する必要がある。run索引は[case README:375](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:375)。

**ファイル変更・forge起動なし。plan 未反映。** 提案の反映先は `tooling-nozzle-sern-chain.md` §5.1 R7a／R7b。
