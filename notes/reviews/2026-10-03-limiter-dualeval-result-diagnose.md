# codex 諮問 (diagnose): limiter-dualeval-result

- **brief**: [`notes/reviews/briefs/2026-10-03-limiter-dualeval-result.md`](../../notes/reviews/briefs/2026-10-03-limiter-dualeval-result.md)
- **plan**: [`plans/active/limiter-inlet-column-oscillation.md`](../../plans/active/limiter-inlet-column-oscillation.md)
- **date**: 2026-10-03
- **commit**: `75e61d3a` (feature/species-transport)
- **codex**: effort `high`, 3.9 min, rc=0
- **結論**: **共通の状態・BC・勾配から通常ψ／保存ψへ分岐する二重評価に直し、100標本を取り直すことを唯一の次手として推奨します。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | `A'−A` が状態変化の影響も含む上限：**却下** | [main.cpp:2162](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2162) は復元なしで `B→A'→A` を実行します。後半の変化が小さくても、最初の変化を拘束しません。**同じ組立入力から分岐して比較**してください。 |
| **Major** | 事前登録に従い第1仮説を支持と確定：**要再検証** | [plan:138](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:138) の状態不変という前提が未達です。現結果は「大きな直接感度と整合するが、組立順序の交絡が残る」まで。`roK≈0.94`・`roUy≈0.33` は、前提を満たしても登録基準では保留です。 |
| **Minor** | `S_B/S_A` を追加：**採用**。既存CSVからの算出：**却下** | [main.cpp:2198](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2198) が保存するのは `‖A‖`・`‖B−A‖`・`‖A'−A‖`。方向情報がありません。**`ΣB²` と `ΣA·(B−A)` を追加**してください。 |
| **Major** | 次に長期凍結または入口の数値処理変更へ進む：**却下** | [plan:29](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:29) のスコープと、未成立の同状態比較からは早計です。まず既存の診断を修正。診断修正は現planの§5.1・§6に反映し、入口近傍集合の変更・局所1次化は別の修正planで扱ってください。 |

結論: **共通の状態・BC・勾配から通常ψ／保存ψへ分岐する二重評価に直し、100標本を取り直すことを唯一の次手として推奨します。**

第 1 仮説: **入口2列のψの差し替えは、同じ状態の入口密度残差に直接、大きく寄与する。** 確度: **中**

  根拠: [limiter_d.cu:459](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:459) は実ノードの指定集合だけを上書きし、SLAUの組成再構成にも `limiter_ro` が伝わります（[convectiveFlux_slau_d.inc.cuh:349](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:349)）。入口境界半割面の質量流束は境界状態から作ります（[convectiveFlux_boundary_d.inc.cuh:192](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:192)）。  
  `case/16.nozzle_wys/run_0534_liminlet_dualeval2/` の **D(res_ro)=1.070/1.062** はこの仮説と整合しますが、これはplan記載値であり、今回独立集計した実測値ではありません。

  反証条件: 共通入力・差し替え・再評価誤差の検査を満たした比較で、**D(res_ro)≤0.1が前後両窓で再現**すること。その場合、「この保存ψへの差し替えによる大きな直接寄与」を棄却します。

第 2 仮説: **全組立の反復によるEOS・BC・派生量の変化が、現Dを大きく汚染している。** 確度: **低**。変化する経路は確認済みですが、寄与の大きさは未確認です。

判別 A/B:

- **変えるのは入口集合Sの5個のψだけ**。EOS・壁射影・BC・勾配・渦粘性などの共通準備を1回行い、最初の流束評価直前（[main.cpp:1874](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1874)）から、通常ψと保存ψで残差を評価します。
- 両枝の入力を同一にし、`res_*`・`massflux`・面組成・作業配列への書込みを分離または復元。**保存量だけの復元では不十分**です。時間更新には通常枝の残差と作業状態だけを使います。
- `iStep=100` のψを保存し、**101〜200を含む100標本**、前半・後半各50標本。現行の200外反復では199までになり、99標本でした。
- 主指標は従来のD。追加で `Q=√Σ‖B‖² / √Σ‖A‖²`、全保存量の応答、共通入力のバイト一致、通常ψ同士の再評価誤差を記録します。再評価誤差は差し替え差の1%以下を検査条件として事前記載してください。
- **A：D≥1が両窓** → 第1仮説を支持し、「直接寄与は小さい」を棄却。**B：D≤0.1が両窓** → 第1仮説を棄却。中間・入力不一致・誤差条件未達は保留。収束や床の原因確定には使いません。

やらない方がよいこと: **1〜2 ulpだから無害と扱うこと、D≈1を残差低減と読むこと、入口ノードを極値集合から外して改善しただけで丸めを真因と決めること。** D≈1.06だけから分かるQの範囲は約 **0.06〜2.06** で、低減も増大も可能です。

呼び出し側の前提への異議:

- **内部`roe`の変化源は、壁射影よりEOS再構成が有力です。** TP経路はfloatの `roe/ro−ek` を入力にし（[dependentVariables_d.cu:83](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/dependentVariables_d.cu:83)）、反転した温度から**床の作動有無によらず**`roe`を書き戻します（同:209）。この往復はビット冪等とは限りません。壁射影は壁フラグのある実ノードだけに作用し、運動量をゼロにします（[nodeWallDirichlet_d.cu:24](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:24)）。後段が運動量を書き戻さなければ、次回のKE除去はゼロです。**壁に近いという位置情報だけでは壁射影に帰属できません。**
- **ghostの変化はBC再計算で説明可能です。** 圧力入口はownerの圧力・Machから状態を作り直し、ghostの`ro`・`roUx`を書きます（[boundaryCond_d.cu:1111](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/boundaryCond_d.cu:1111)、同:1193）。組成境界は `roY[ig]=ro[ig]·Yin` です（[speciesTransport_d.cu:128](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:128)）。ただし、報告されたghostがどの境界かは未確認です。
- **保存量比較は組立入力全体の比較ではありません。** EOSの初期推定には前の`T`も使われます（[dependentVariables_d.cu:116](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/dependentVariables_d.cu:116)）。また各組立は`roM`・`roeM`なども更新します。現実装は、前回レビューの「`assembleResidual`をそのまま複数回呼ばない」という指定（[前回レビュー:39](/home/sano/work/forge-species/notes/reviews/2026-10-03-limiter-psi-freeze-diag-diagnose.md:39)）を満たしていません。

不足情報: 対象の`run_0533`・`run_0534`はローカルにありません。**`psi_dualeval.csv`、実効config、変更ノードと境界の対応、実行バイナリの同一性**が必要です。関連runの収束・準定常VERDICTも未確認なので、収束・定常性についての結論は出しません。主要診断・EOS・BCソースは`aef6c1e1`から現HEADまで差分なしを確認しました。ファイル変更・forge起動は行っておらず、提案は **plan未反映**です。
