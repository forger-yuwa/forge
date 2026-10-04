# codex 諮問 (diagnose): farfield-v0b-result

- **brief**: [`notes/reviews/briefs/2026-10-03-farfield-v0b-result.md`](../../notes/reviews/briefs/2026-10-03-farfield-v0b-result.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-10-03
- **commit**: `eaa598c8` (feature/sern-design)
- **codex**: effort `high`, 4.7 min, rc=0
- **結論**: ?
- **extra**: `notes/investigations/2026-10-03-farfield-evidence/c46_V0B_REPRO_20261003.txt`, `notes/investigations/2026-10-03-farfield-evidence/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | FAIL は規則の偏りなので同等と扱う | **却下**。記録は全9保存量で `VERDICT: FAIL`、したがって「登録条件で同等性を示せなかった」が正式な結論。ただし、コード回帰を証明したわけでもない。[判定原本:19](/home/sano/work/forge-sern-design/notes/investigations/2026-10-03-farfield-evidence/c46_V0B_REPRO_20261003.txt:19)。元の FAIL を保持し、原因診断を別試験にする。 |
| **Major** | 初回ダンプ一致で、非 farfield 経路の変更を除外できる | **却下**。ダンプは最初の評価の `massflux` と `ro,Ux,Uy,Uz,P,sonic` だけ。保存量更新・乱流／化学種輸送・陰解法更新を網羅しない。[実装:565](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:565)。対案は、同一入力と固定加算順での完全な1更新比較。 |
| **Major** | バイナリ間の差は farfield と出力専用変更だけ | **要再検証**。候補コミット `2fa3826c..a6ceee0b` には、再構成速度を変更する `FORGE_DIAG_FACE_VEL_CELL` の追加もある。[読込:227](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:227)、[作用箇所:280](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:280)。発動した証拠はないが、実バイナリのソース対応と実行時環境を照合するまで除外できない。 |
| **Major** | 未決のまま限定受理で `accepted` に進む | **却下**。§2 の縮小決定は V2c・V2d-2 を扱っており、新しく判明した V0 後半の FAIL の免除ではない。[plan §2:24](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:24)、[§6 V0:176](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:176)。限定運用を維持し、`in_progress` に留める。 |

**結論:** V0 後半の FAIL を保持し、次は旧・新の同一入力に対する「固定加算順での完全な1更新 A/B」を事前登録する。

**第 1 仮説:** 非決定的な積算による再実行変動と、未校正の最大値比較が FAIL の主因である。**確度: 中**

根拠：

- `case/46.sern_design/run_1032_v0b_{o,n}{1,2,3}/` の集計では、`roUx` の `S_old / S_new / D` は `0.2778 / 0.3838 / 0.5125`、`roOmega` は `4.459e5 / 7.511e5 / 6.332e5`。新旧差と同一バイナリ内の差が同じ桁にある。ただし、これだけで同分布とは言えない。
- 流れの残差だけでなく、スカラーの残差・輸送対角にも `atomicAdd` がある。[流れ:661](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:661)、[スカラー:181](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/scalarTransport_d.cu:181)。
- 規則の問題は「9組対3組」だけではない。**3回の観測範囲を、後続標本の許容限界として使っている。** 単一スカラーを同じ連続分布から独立に3個ずつ取る簡単な場合でも、`D ≤ S_old` が成立するのは6個の最小・最大がともに旧群に入る場合だけで、確率は **20%**。全20割当の列挙でも PASS 4・FAIL 16だった。これは今回の全場試験の FAIL 確率ではないが、**同分布でも頻繁に落ちる規則であることの反例**にはなる。9保存量も互いに独立な証拠ではない。

反証条件：入力・加算順を固定し、同一版の反復がビット一致する条件でも、新旧の面寄与・残差・対角・更新量に再現する差が残れば、「再実行変動と判定規則だけで説明できる」を棄却する。

**第 2 仮説:** 共有するスカラー輸送等の変更、または生成コードの違いが決定的な作用素差を生んだ。未確認。farfield が無い場合に面値ポインタが `nullptr` となる経路は確認したが、それは更新全体の一致証明ではない。

**第 3 仮説:** 入力・環境変数・ビルド条件の交絡。未確認。特に診断介入の無効化と、実バイナリのソース対応が不足している。

**判別 A/B:** 変える要因は旧版／新版だけ。`run_0971` の同じ保存量・メッシュ・BC・物性を使い、元の内反復数を保った**完全な1 outer step**を各2回実行する。

- 両版に同じ診断機構を適用し、更新へ入る積算を同一順序に固定する。残差だけでなく、実際に使う勾配・輸送対角・陰解法組立ても対象。面寄与、段別残差、対角／行列、`ΔQ`、更新後の全9保存量を比較する。
- **前提ゲート:** 旧同士・新同士がビット一致しなければ診断は不成立。差を新旧の作用素差と解釈しない。
- **A＝同一版内も新旧間も一致:** 試験した状態・1更新における決定的な作用素差を棄却する。200 step 全体の同等性や V0 合格には昇格させない。
- **B＝同一版内は一致、新旧間だけ不一致:** 再実行変動だけの説明を棄却し、最初に異なる段へ原因を限定する。

これは既存の `FORGE_DUMP_MASSFLUX` だけでは実施できず、診断実装が必要である。**通常の1 step比較を「決定的な比較」と呼ぶことはできない。**

**やらない方がよいこと:** `S_old` を結果に合わせて `max(S_old,S_new)` や倍率付き閾値へ変更すること、同じ規則のまま本数だけ増やして PASS を待つこと、`atomicAdd` を真因と断定すること。旧3＋旧3の負対照は規則の誤判定を示すには有用だが、旧同士でも FAIL したことは新旧同等性の証明にならないため、今回の推奨一手には選ばない。

**呼び出し側の前提への異議:** 「同等性未証明」と「規則に統計的問題がある」は両立する。どちらか一方を選ぶ問題ではない。また、初回ダンプ一致と多 step 更新の不変性を同一視しない。

**不足情報:** ローカルには対象6 run の実体がなく、最終 HDF5・入力ハッシュ・実効設定・実行時環境・ビルド条件・評価スクリプトを再監査できなかった。`run_0991` の新バイナリ識別子は索引で `778179984023`、今回の新バイナリは `b0240cd7…` と異なるため、初回一致の証拠を今回へそのまま転用できない。[索引:107](/home/sano/work/forge-sern-design/notes/investigations/2026-10-03-farfield-evidence/INDEX.tsv:107)。また、[case の run 台帳](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:345) には今回の `run_1032` 群の追記が必要。

ファイル変更・forge 実行はしていない。提案は **plan 未反映**。
