# codex 諮問 (diagnose): limiter-psi-freeze-diag

- **brief**: [`notes/reviews/briefs/2026-10-03-limiter-psi-freeze-diag.md`](../../notes/reviews/briefs/2026-10-03-limiter-psi-freeze-diag.md)
- **plan**: [`plans/active/limiter-inlet-column-oscillation.md`](../../plans/active/limiter-inlet-column-oscillation.md)
- **date**: 2026-10-03
- **commit**: `51d59632` (feature/species-transport)
- **codex**: effort `high`, 4.4 min, rc=0
- **結論**: **長期の F-in／F-all を投入する前に、同じ状態・BC・勾配で通常 ψ と入口2列の保存 ψ を切り替える、短い残差二重評価 A/B を1つ行ってください。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 相当の欠陥は確認していません）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 凍結を「状態・BC・勾配を変えずに ψ だけを固定」と解釈：**却下** | 上書きした瞬間には成立します。しかし残差→状態更新→次回の BC・勾配再計算という経路があります（[main.cpp:1835](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1835)、同:1847、2083、2130）。長時間の凍結は閉ループ全体への介入です。**直接の残差応答は、同じ状態で通常 ψ／保存 ψ を差し替えた二重評価で測る**べきです。 |
| **Major** | `1/5`・`0.5〜2`・「L0 床の3倍」で原因を確定：**却下** | [plan:136](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:136) の判定は強すぎます。凍結は再構成と固定点を変えます。改善しても丸め・BC・非線形反復のどれが起点かは決まらず、改善しなくても別の残差との相殺変化があり得ます。**比は介入の効果量として使い、因果の一方向性や唯一原因の判定には使わない**。L0 の `1.38e-8` は参考到達値に留めます。 |
| **Major** | `10000 step`・出力 `1000` で床を判定：**却下** | 同じ [plan:136](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:136) では、判定区間の `R_入口` は約6標本しかありません。L0 の緩和時間を凍結作用素へ転用する根拠もありません。**10000 step は途中確認まで**。長期試験なら凍結後を独立区間として全残差の `check_convergence` と対象量の `check_quasisteady` を確認し、下降中なら延長する必要があります。 |
| **Major** | 「化学種・受動種は対象外」：**配列の説明として採用、影響範囲として却下** | SLAU の組成再構成は実際に `limiter_ro` を使います（[convectiveFlux_slau_d.inc.cuh:349](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:349)）。受動種固有の ψ は別ですが、質量流束と上流選択は変わります（同:579、593）。**5配列への介入が、その読者すべてへ伝播する試験**と定義し、`speciesFaceReconstruction`・`multispeciesRhoYCommonLimiter`・`reconT` の実効値を記録してください。 |
| **Major** | `nCells_all` 全体への座標マスク・周期対応：**要再検証** | 通常計算は全配列を1で初期化してから実ノードの ψ を計算します（[limiter_d.cu:451](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:451)、同:475）。周期 node は group の ψ を一致させます（同:493）。**保存容量は `nCells_all`、node の上書き対象は `ic<nCells`**を基本とし、ghost を除外する。周期 group をマスクが分断する場合は拒否するか group 全体へ適用してください。 |
| **Minor** | wrapper の挿入位置：**条件付き採用** | [limiter_d.cu:513](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:513) 直後なら通常・周期両経路の後です。ただし `reconT=1` の `limiter_T` は凍結されません。また `LIMG1` のカウンタは既にカーネル内で計算済みです（同:399）。**後段の印字が CLEAN でも凍結後の有界性は保証しない**ため、適用した ψ に対する検査を別に行ってください。 |
| **Minor** | 診断専用実装の承認履歴：**採用、文書整合は必要** | `#5b` のユーザ了承は確認しました。一方、[plan:30](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:30)、同:151、155 は「凍結なし／ソース変更なし」のままです。**診断限定の例外と変更範囲を記載し、ソース変更なしを理由とするレビュー免除を更新**してください。追加のユーザ承認を求める理由にはしません。 |

結論: **長期の F-in／F-all を投入する前に、同じ状態・BC・勾配で通常 ψ と入口2列の保存 ψ を切り替える、短い残差二重評価 A/B を1つ行ってください。**

第 1 仮説: **入口2列の ψ の変動は、同じ状態における入口密度残差へ直接、大きく寄与している。** 確度: **中**。

  根拠: 内部面の SLAU 再構成は5変数の ψ を読みます。一方、入口境界半割面の質量流束は境界状態から `Sρ_R U_n,R` として作られ、ψ を直接読みません（[convectiveFlux_boundary_d.inc.cuh:192](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:192)）。したがって、同じ BC・状態のまま内部面だけを変える介入は可能です。ただし、寄与の大きさは run 未確認です。

  反証条件: ψ の差替えが実際に効いていることを確認した上で、二重評価による残差差の入口ノルムが、通常残差の **10%以下**に連続窓で留まれば、「直接の大きな寄与」という仮説を棄却します。長期のフィードバック増幅まで棄却する条件ではありません。

第 2 仮説: **ψ と残差の相関の主部分は状態・BC・勾配の変動を共通入力とする応答、または反復を通じた増幅である。** 確度: **低、未確認**。短い同状態比較で直接寄与が小さければ、こちらを残します。

判別 A/B:

- **変更点は入口集合 S の5個の ψ だけ**。通常の軌道を200外反復進め、100回目の ψ を保存。101〜200回目の各状態で、A＝その時点の通常 ψ、B＝S だけ保存 ψ として残差を評価します。**状態更新には A だけを使用**します。
- 通常の `assembleResidual` をそのまま二度呼ぶ設計は避けてください。これは状態射影や BC 更新も行います。分岐は BC・勾配確定後とし、B の `res_*`・`massflux`・面組成などは診断用の別領域へ出す必要があります（[main.cpp:1812](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1812)、[convectiveFlux_d.cu:166](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:166)）。
- 主指標は **D = √ΣₙΣᵢ∈S(Rᴮᵢ,ₙ−Rᴬᵢ,ₙ)² / √ΣₙΣᵢ∈S(Rᴬᵢ,ₙ)²**。補助として通常／差替え後の残差ノルム比、全保存量の応答、適用 ψ と保存値の最大差を記録します。前半・後半の各50標本でも確認します。
- **A：D≥1**が両窓で再現 → 直接寄与が小さいという説明を棄却し、第1仮説を支持。ただし残差を減らす方向とは限らないため、ノルム比も併記します。
- **B：D≤0.1**が両窓で再現 → 第1仮説を棄却し、第2仮説を残します。中間は保留です。これは感度試験であり、収束・定常性の主張には使いません。
- 残差組立てには `atomicAdd` があるため、同じ通常 ψ の再評価差を測り、差替え効果が組立て再現誤差を十分上回ることも確認します（[convectiveFlux_slau_d.inc.cuh:638](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:638)）。

やらない方がよいこと: **F-all が L0 の床へ近づいたことだけで「床は ψ の動的変化そのもの」と結論すること、10000 step で改善しないことを原因否定に使うこと、旧バイナリの `run_0527` を無検証で新実装の対照にすること。**

呼び出し側の前提への異議:

- **呼び出し回数**: 定常陰解法では1外反復につき1回です（[main.cpp:2274](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2274)、同:2083、1866）。`nStepInner=5` はその後の線形 sweep 数です（同:1978）。通常経路で1始まりの `N=100` なら、第100更新前に保存し、最初に値が変わり得るのは第101更新前です。ただし起動時の `twoPhaseAudit` が有効なら追加呼び出しがあります（同:2711）。**call番号と外反復番号を併記**してください。
- **restart**: プロセス内カウンタと保存 ψ は通常の場の restart では引き継がれません。長期試験を再起動で延長すると再凍結になります。保存 ψ を永続化するか、同じプロセスで継続する必要があります。
- **x 閾値**: 提示された3列の座標が実データどおりなら `x<−0.059 m` は妥当です。ただし実ノードID・個数・選択列を実測して固定してください。変更した面流束は両端へ加算されるので、**2列目と3列目の間の面を通じ、3列目の残差にも直接効きます**（上記 SLAU:638、644）。「入口2列だけの残差を変える介入」ではありません。
- **状態固定と介入対象固定は別**です。長期凍結試験自体は有用ですが、答えられる問いは「この凍結で反復の挙動がどう変わるか」です。

不足情報: `case/16.nozzle_wys/run_0524_floor_dry_L1/`、`run_0527_liminlet_ctrl_cfl2/`、`run_0528_liminlet_cfl1/`、`run_0529`〜`0532` はローカルにありません。実効 config、メッシュ、残差履歴、HDF5、収束・準定常 VERDICT が必要です。**既報の `1.93e-7`、床の数値、ALL STEADY は今回は独立検証していません。** 指定主要ソースは `8ce8da3a` との差分なしを確認しましたが、実行バイナリの同一性は未確認です。ファイル変更・forge 起動は行っておらず、提案は **plan 未反映**です。
