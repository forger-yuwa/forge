# codex 諮問 (diagnose): g4-admission-jacobian

- **brief**: [`notes/reviews/briefs/2026-10-11-g4-admission-jacobian.md`](../../notes/reviews/briefs/2026-10-11-g4-admission-jacobian.md)
- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-11
- **commit**: `a3695e9b` (feature/sern-design)
- **codex**: effort `high`, 4.9 min, rc=0
- **結論**: **投入を止めたまま、`run_1080`／`run_1081` の最終 node HDF5 で、対応ヘキサごとのヤコビアンを比較する 0-step A/B を行う。**
- **extra**: `plans/active/convection-zero-thickness-edge-reconstruction.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 2 / Minor 1）

| 重大度 | 採否 | 根拠と対案 |
|---|---|---|
| **Major 1** | **提示された読み替えによる投入は却下** | 再生成すると、B は「変形領域に 0.65 未満を新しく作らない」を満たさない。新規超過は **4,180 ヘキサ**。現行判定は全域最小値しか比較しておらず、局所劣化を表示しない。[判定コード:770](/home/sano/work/forge-sern-design/case/46.sern_design/diag/te_wake_grid_check.py:770)。対案は、実入力 HDF5 の対応するヘキサごとに A/B を比較すること。 |
| **Major 2** | **g4 基準格子の品質受入れは要再検証** | A/B とも `MESH_QUALITY: SOFT-PASS` だが、AR 最大 **12865.6**、AR > 5000 が **691 個**ある。[品質記録:4](/home/sano/work/forge-sern-design/notes/investigations/2026-10-11-g4-admission/run_1081_tewake_g4_B10_m10/MESH_QUALITY.txt:4)。手順は SOFT-PASS の外れ値の位置確認を要求する。[手順:382](/home/sano/work/forge-sern-design/procedures/calculation-workflow.md:382)。旧 g4 の実績は参考になるが、現在の m10_on の品質保証にはならない。対案は、基準格子の例外受入れと局所変形の劣化判定を分けて記録すること。 |
| **Minor 1** | **復帰角の現行 PASS は採用** | 未丸め値は **5.999850829970171°** で、実装も表示値ではなく未丸め値を 6° と比較している。[保存値:200](/home/sano/work/forge-sern-design/notes/investigations/2026-10-11-g4-admission/run_1081_tewake_g4_B10_m10/TE_WAKE_PROVENANCE.json:200)、[比較:693](/home/sano/work/forge-sern-design/case/46.sern_design/diag/te_wake_grid_check.py:693)。近いことだけを理由に事後に FAIL へ変更しない。ただし余裕は **0.00014917°** しかなく、この格子に限定した合格とする。 |

結論: **投入を止めたまま、`run_1080`／`run_1081` の最終 node HDF5 で、対応ヘキサごとのヤコビアンを比較する 0-step A/B を行う。**

第 1 仮説: **全域最小値は基準格子の遠方セルに固定されているが、後縁の局所変形は別のセルで品質を悪化させている。** 確度: **高**

根拠: 対象 YAML からメモリ上で両格子を再生成し、メッシャの **10 桁出力 → float32 格納**を再現して、各ヘキサの全 8 頂点を検査した。[出力形式:698](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:698)、[ヤコビアンの定義:391](/home/sano/work/forge-sern-design/case/46.sern_design/diag/te_wake_grid_check.py:391)。

| 検査量 | A：変形なし | B：1.0 H |
|---|---:|---:|
| 全域の最小値 | 0.4523499229 | 0.4523499229 |
| 最小値が 0.65 未満のヘキサ数 | 40,700 | **44,880** |
| 変形節点を含む領域内の同ヘキサ数 | 6,105 | **10,285** |
| 同領域内の最小値 | 0.5264462965 | **0.5174910493** |

さらに、**新しく 0.65 未満になるヘキサは 4,180 個、既に 0.65 未満でさらに悪化するヘキサは 6,105 個**。変形節点を含まないヘキサの差は 0 だった。

具体例は、再生成時の 0 始まりヘキサ index で次のとおり。

- **931802**：0.73435596 → **0.64715481**。B の重心位置は `(x/H, y/H, z/H) = (1.88407, −0.72072, 0.33059)`。
- **931643**：0.64734454 → **0.57541870**。位置は `(1.88407, −1.34764, 0.21258)`。

再生成の移動節点 **127,848 個**、変形領域 **134,750 ヘキサ**、復帰角、全域最小値は保存記録を再現した。B の YAML ハッシュも来歴記録と一致し、使用したメッシャ・runner・判定コードには `d734392d` からの差分がない。ただし、**AWS の実入力 HDF5 自体の再検査ではない**。

反証条件: 来歴を照合した実入力 HDF5 の比較で、新規の 0.65 未満が 0 個、既存の低値セルにも有意な悪化がなく、上記の局所劣化が再現しないこと。その場合は再生成との相違を調べる。

第 2 仮説: **g4 の第一層厚の下流ブレンドと中間線変形の組合せが、下側バンドの局所劣化を増幅している。** 確度: **中・原因分離は未確認**。第一層厚は x 依存で変わり、その値と中間線位置から帯全体の節点配置を再計算している。[ブレンド:241](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:241)、[層配置:267](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:267)。

判別 A/B: **変更点は `te_wake_blend_H` の 0／1.0 だけ、時間積分は 0 step。** 次の既存入力を読み取り比較する。

- A：`case/46.sern_design/run_1080_tewake_g4_A0_m10/sern.h5`
- B：`case/46.sern_design/run_1081_tewake_g4_B10_m10/sern.h5`

接続・番号対応・署名を確認し、「動いた節点を 1 個以上含むヘキサ」を変形領域とする。各ヘキサの最小 scaled Jacobian、新規の 0.65 未満、既存低値セルの悪化、領域外の差を集計する。

→ **結果 A：新規超過または既存低値セルの悪化が再現する。** 「失敗は基準格子だけに由来する」を棄却する。閾値の読み替えでは進まず、g4 の下側バンドの分布と変形の接続を扱う格子設計へ戻す。

→ **結果 B：両方とも 0 で、低品質部分が不変と確認できる。** 今回の局所劣化仮説を棄却する。基準格子の品質を別途受け入れたうえで、旧 FAIL を保存し、適用対象・基準格子署名・局所比較条件を明記した新しい投入規則として改訂を検討できる。旧規則で PASS だったとは扱わない。

やらない方がよいこと: **0.65 を 0.45 付近へ下げること、全域最小値の一致だけで投入すること、L_b を調整して再挑戦すること。** 現在の保存判定は `ADMISSION VERDICT: FAIL`、品質は A/B とも `SOFT-PASS`、双対閉性は `PASS` のまま維持する。

呼び出し側の前提への異議: **40,700 個は再生成では A の値で、B は 44,880 個だった。** また、最悪セルが変形領域外にあることは、低品質セル全体が変形領域外にある証拠ではない。双子の記録も **11,620「組」ではなく、5,810 組・11,620 節点**である。[双子の検査記録](/home/sano/work/forge-sern-design/notes/investigations/2026-10-11-g4-admission/run_1081_tewake_g4_B10_m10/TWIN_IC_VERDICT.txt)。

不足情報: 実入力 HDF5 による上記の対応比較、AR 外れ値 691 個の位置と受入れ根拠。時間積分へ進む際には、別途停止条件となっている g3 の有限化入力による EOS 比較の合格記録も必要。[plan:98](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:98)。

**ファイル変更・forge 起動なし。plan 未反映**。採否と追加検査の反映先は `tooling-sern-te-wake-grid.md` §5.1 #4・§6 とし、呼び出し側で記録する。
