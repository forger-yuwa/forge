# codex 諮問 (diagnose): line-viscous-jprobe

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-jprobe.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-jprobe.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `b6298420` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **正式判定は (d) で閉じ、追加確認は TP 面エンタルピー評価の float／double A/B 一組だけに限定し、その後はマスク5の評価を広げず本線の仕上げ込み総壁時計評価へ戻る。**
- **extra**: `case/45.isobutane_m6_d155/_jprobe/s0p7n_compare.txt`, `case/45.isobutane_m6_d155/_jprobe/s1p7n_compare.txt`, `case/45.isobutane_m6_d155/_jprobe/s1p7r_compare.txt`, `case/45.isobutane_m6_d155/jprobe.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0、Major 4、Minor 2）

| 論点 | 判定 | 根拠と対案 |
|---|---|---|
| S0・S1 とも分岐 (d) | **採用** | JSON を再集計した。ε = 1e-2 の合格行を保持し、不合格行だけ ε = 1e-1 で再試行すると、有効行は各方向とも **S0 3/20、S1 4/24**。半数未満なので結論は変わらない。[判定規則](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:336) |
| H-a「リミッタ等の不連続が原因で接線検証は成立しない」 | **Major・却下** | 「今回の幅では差分が確定しない」は観測だが、原因の説明は未立証。`min/max` は入力が連続なら値の不連続を作らない。[リミッタ本体](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/limiterFunctions_d.cuh:41)と面間の最小値・clamp だけでは「跳び」の説明にならない。**測定不成立と、その原因未同定を分けて記録する**。 |
| H-b を探索的観測として残す | **採用。ただし Major の留保あり** | 数値差は記録してよい。しかし比較対象には意図的な安定化項が含まれる。`lineViscCoupling: 3` は[スカラー対角を加算](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:966)し、[比較コード](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jprobe.py:133)が引くのは時間項だけ。差をそのまま「Jacobian の欠陥」や破綻原因に帰属させず、**安定化項を含む近似作用素と評価残差の差**と呼ぶ。 |
| H-c の追加確認 | **採用・Major** | ブリーフにない **TP 面エンタルピーの float 評価**がある。これはエネルギー行の微小差分を直接変え得る。CPG/TP 比較より、下記の精度だけの A/B を一組行う価値がある。 |
| H-d「float では大きさを説明できない」 | **Major・要再検証** | 面エンタルピー経路が見積もりから抜けている。同一入力の再評価差では、入力依存の量子化は除外できない。mass の跳びについても、該当面の感度・相殺前の流束・分岐の変化を測っていないため、除外済みにはできない。 |
| 集計・指標の表現 | **Minor 2件・修正** | ① S1・p5・ライン12の `momx` の再現差は **0.002514**で、ブリーフの「≤ 1e-3」から外れる。② 熱伝導・仕事・残りは**ベクトルのノルム比**で、足して100%になる寄与率ではない。S1・p5・ライン65では heat **0.0440**、rest **1.0039**。[元出力](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_jprobe/s1p5n_compare.txt:5)。「残りが96〜100%を占める」ではなく、各ノルム比をそのまま示す。 |

結論: **正式判定は (d) で閉じ、追加確認は TP 面エンタルピー評価の float／double A/B 一組だけに限定し、その後はマスク5の評価を広げず本線の仕上げ込み総壁時計評価へ戻る。**

第 1 仮説: **H-c の大きなエネルギー行の差には、微小摂動に対する面エンタルピーの量子化・実質的な凍結が寄与している。** 確度: **中**

根拠:
- [対流流束](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:399)は `thermo_h_mix_f` で面エンタルピーを作る。double ビルド用の[オーバーロード](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/thermo_d.cuh:1079)も、組成と温度を float に変換し、float を返す。**指定された commit `846727be` にも同じ経路があることを確認した。**
- この経路は `thermoFloat` の切替とは独立である。`thermoFloat: 0` だけでは除去できない。
- S0・p7・ε = 1e-6 の最大 |δρ/ρ| は **7.03e-10**。[入力記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_jprobe/s0p7n_fields.json:15)。面温度が同じ float 値に丸められる範囲では、h(T) の温度応答が消え、境界を跨ぐと段差になる可能性がある。
- ライン2183のエネルギー行は、差分再現誤差が S0 **9.31e-6**、S1 **2.46e-5**と小さい一方、|J_a p|/|J_t p| は **0.0284／0.0242**。この再現性は「量子化された評価経路内で再現する」ことを示しても、連続な TP モデルの接線を測れている保証にはならない。

反証条件: 面エンタルピー評価だけを double にしても、有効な差分同士でライン2183のエネルギー応答が **10%未満しか変わらず**、約0.98の不整合が残れば、この異常に対する主要因という仮説を棄却する。**面エンタルピーだけでは、同じ入力状態の mass 行の跳びは説明しない。**

第 2 仮説: **人工スカラー対角、一次 FVS の近似、凍結物性と再構成残差との差が、探索で見えた作用素差の一部を説明する。** 確度: 中。項の存在はコードで確認できるが、寄与の大きさと破綻への因果は未確認。

第 3 仮説: **mass 等の差分不成立には、音速など別の float 経路、または再構成・境界等の切替が寄与する。** 確度: 低。出所は未同定。[γ・音速の float 経路](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/dependentVariables_d.cu:186)は候補に残す。

判別 A/B:
- **変える一点**: TP 対流流束の面エンタルピー評価精度。A は現行 float、B は同じ NASA 係数・datum を用い、温度入力から h の出力まで double で評価する診断分岐。既存設定にこの切替はないため、追加するなら診断専用とする。
- **固定条件**: S0、保存済み p7、元の摂動支持領域、元の D/K。p を再計算しない。EOS、輸送物性、リミッタ、組成、BC は共通。
- **長さ**: 各側で q0・q0再評価・±1e-6・±5e-7を各1 step、計12評価。新しい run に出し、評価対象は事前に **ライン2183の自由なエネルギー行**へ固定する。
- **測定**: 面の Δh、エネルギー行の J_t p、ε/ε/2差、再評価ノイズ、固定した J_a p との差。比較の有効条件は再現誤差 ≤1%、再評価ノイズ／信号 ≤1e-3。
- **結果の分岐**: Bで J_t p が10%以上変われば、面熱力学の精度依存を支持する。さらに作用素との相対差が半減すれば、H-cの大差の説明としても支持する。変化が10%未満なら、この条件での主要因説を棄却する。Bでも差分が成立しなければ判別不能のまま終了し、探索を追加しない。

やらない方がよいこと:
- ε = 1e-6 の結果で、事前登録の (d) を (a) に置き換える。
- 跳ぶ節点だけを事後に削除して「そのラインで微分が成立した」と扱う。探索結果は**状態×方向×ライン×保存量行**ごとに、拘束行以外を削らない既存ノルムで掲載する。
- CPGへの変更、ρEだけの方向、`thermoFloat: 0` を、この面エンタルピー経路の単因子検査とみなす。
- 原因未同定を理由に、マスク5の長期評価や熱伝導Kの列削除を続ける。

呼び出し側の前提への異議:
- `jumps.py` が数えるのは「最小幅の中心差分から10%以上外れた点」であり、数学的な不連続点ではない。[実装](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_jprobe/jumps.py:14)。また表示する k は異なる値の先頭15個だけなので、これだけで「k = 4〜26に集中」とは確定できない。
- p は選択ライン上だけ非零であり、0.7pの試験も**切り出した補正の有限振幅応答**である。全域の平均流・乱流更新を含む実際の1 stepの予測誤差とは区別する。
- 微小差分の跳びが残差床に関係する可能性はあるが、ライン陰解法の破綻との因果は未証明。別の残作業として「残差評価の精度依存・差分成立域の監査」を登録し、**今回の短いA/B以上には広げない**。
- plan §6.14には正式判定(d)、有効行一覧と採用ε、探索結果の限定、量子化経路、未同定事項を残す。§5.1にはこの一組と、本線への復帰を記す。値2・3の保留と既定値は維持する。**本回答ではファイル変更禁止に従い、plan未反映。**

不足情報: ローカルには比較JSON/TXTがある一方、元の `_jprobe/npz`、ライン行列、S0/S1のHDF5、対象runの実効設定・VERDICTがない。したがって集計とコード経路は確認したが、節点別応答、床未到達、AWSバイナリ実体までは独立検証していない。`forge` は起動していない。成果物の索引は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:175)。
