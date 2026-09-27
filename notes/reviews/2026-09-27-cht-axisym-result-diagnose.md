# codex 諮問 (diagnose): cht-axisym-result

- **brief**: [`notes/reviews/briefs/2026-09-27-cht-axisym-result.md`](../../notes/reviews/briefs/2026-09-27-cht-axisym-result.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **date**: 2026-09-27
- **commit**: `5696a53b` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 4.3 min, rc=0
- **結論**: **条件つきで、現時点の `done`・拒否解除は不可；まず V-ax0 の未完了項目を閉じ、総括と仕様文書を訂正してから codex result レビューへ進める。**
- **extra**: `case/62.conjugate_disk/README.md`, `case/61.conjugate_annulus/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **条件つきで、現時点の `done`・拒否解除は不可；まず V-ax0 の未完了項目を閉じ、総括と仕様文書を訂正してから codex result レビューへ進める。**

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 「V-ax0 PASS、解除条件がそろった」 | **却下・Major** | [plan:130](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:130) は軸上辺の負例と method 1 の HDF5 属性を未了と明記。一方、[登録:145](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:145) は両方を要求する。出力実装も [output.cpp:317](/home/sano/work/forge-cht/solver_density_cuda/output/output.cpp:317) では `step`・`step_abs` のみで、未検証属性の書き込みがない。**属性を実装し、未実施の負例を実行・保存するまで V-ax0 は未完了とする。** |
| V-ax1 の例外受け入れ | **採用** | [plan:180](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:180) にユーザ決定がある。旧 FAIL と収支 `2.58e-12 > 1e-12` は保持し、「全項目 PASS」へ変更しない。丸めの説明は**当該行列・演算経路の観測**に限定し、FP64 一般の不可避な精度限界とは主張しない。 |
| V-ax2・V-ax2b の限定合格 | **採用** | 下記の実測と判定を再確認した。非一様格子は純伝導 IC 限定。一様 IC の失敗は [流体側 plan §5.1](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:72) に未解決として残し、CHT 合格を流体側の解決と読み替えない。 |
| V-ax3・V-ax4 の比較結果 | **採用** | 壁ダンプから再計算し、V-ax3 は差 `1.18146e-3 K ≤ 3.96418e-3 K`、V-ax4 は `1.89771e-8 K ≤ 2.97005e-8 K` を再現。**登録された有限反復の回帰・再開比較**として扱い、物理解の収束証明には使わない。 |
| 現在仕様の記述 | **修正を採用・Minor** | [methods/boundary.md:493](/home/sano/work/forge-cht/methods/boundary.md:493) は「計画中」、同 `:500` は全軸対称を拒否と記載。完了時には受理条件、検証範囲、IC 制限、再開契約へ更新する。`q_compact` の受理を今回の `q_eff` 検証で保証せず、既存 opt-in を理由なく削除しない。 |

実測の確認結果:

| 根拠 run | 再確認した数値・VERDICT |
|---|---|
| `case/61.conjugate_annulus/run_0005_annulus_r32s16_w20k/` | 壁温誤差 **0.037526 %**。`check_convergence`: **PASS**、`check_quasisteady`: **ALL STEADY**、保存済み G-if: **PASS** |
| `case/62.conjugate_disk/run_0019_disk_r32g_condic/` | 壁温誤差 **0.085774 %**。同じく **PASS / ALL STEADY / G-if PASS** |

残差判定は各 run の全履歴、最終 step `299999`。準定常は更新系列 `20000–299950` の末尾半分を、全界面節点の温度降下・熱流束について再判定した。両 run の最終 `res_300000.h5` の `VALUE/*` に非有限値はなかった。格子感度 **0.006550 K** と円板系列の誤差比 **0.3174 / 0.2978** も評価器で再現した。成果物の索引は [case/61 README](/home/sano/work/forge-cht/case/61.conjugate_annulus/README.md)・[case/62 README](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md)。

主張の範囲の推奨文:
> 定常・node・FP64・`axisymMethod: 0`・`fem2d`・`q_eff` の登録円環／円板条件で連成精度・保存・定常性を確認した；非一様円板は純伝導 IC に限定し、固体単体の収支 FAIL は明記された例外として受け入れる。

第 1 仮説: **今回の完了判断を阻む主因は、連成精度の不足ではなく、V-ax0 の未完了項目を総括で PASS にしたこと。** 確度: **高**  
　根拠: 上記 `plan:130` と `:145`、`output.cpp:317`。  
　反証条件: 対象コミットに未検証属性の出力実装が存在し、軸上辺を含む全登録負例の再現可能な合格記録が提示されること。

第 2・第 3 仮説: 追加しない。今回の解除判断のために流体側の原因仮説を増やす必要はない。

判別 A/B: **CHT なし・`interfaceDiag: 1` の同一入力で、`mesh.axisymMethod` だけ A=0／B=1、各 1 step。** 出力対象の等温壁で、A は対象診断量が有限、B は `iface_q_eff`・`iface_q_eff_raw`・`iface_Qf_eff` が NaN かつ未検証属性を持つことを確認する。  
→ **この対照が成立すれば** method 1 の無効化・識別契約を支持。**B が NaN のみで属性なしなら**値の無効化だけが実装され、V-ax0 未完了という判断を支持する。軸上辺の負例は、これとは別に既登録の未実施項目として完了させる。

やらない方がよいこと: 先に `done` にすること、例外受け入れを PASS に書き換えること、純伝導 IC の成功で一様 IC の失敗を解決済みにすること、今回の不足を埋めるために長時間の連成 run を追加すること。

呼び出し側の前提への異議: **「解除条件はそろった」は受け入れない。** また、粗格子の個別評価器 FAIL と登録された系列合格は区別する必要がある。円環 `run_0009` の熱流束誤差 **0.76358 %**、円板 `run_0011` の壁温誤差 **0.50412 %** は実際に FAIL と保存されており、系列の合格で上書きしない。

不足情報: V-ax0 の軸上辺試験と method 1 属性の保存済み証拠。V-ax1 の丸め A/B は今回は原データから再検証しておらず、原因説明の採用範囲は plan の記録まで。**本回答は plan 未反映**。
