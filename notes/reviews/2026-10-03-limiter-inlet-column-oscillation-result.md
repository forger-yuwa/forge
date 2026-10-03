# codex レビュー: limiter-inlet-column-oscillation (result)

- **plan**: [`plans/active/limiter-inlet-column-oscillation.md`](../../plans/active/limiter-inlet-column-oscillation.md)
- **stage**: `result` (diff base `1826a143~1`)
- **date**: 2026-10-03
- **commit**: `5d721fcc` (feature/species-transport)
- **codex**: effort `high`, 5.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **focus**: 調査 plan を Z (受け入れて閉じる) で閉じる前の点検。(1) 変更ログ・§4.6・methods/limiter.md の新節の主張が run の実測 (plan 本文の数値) と矛盾・過大でないか、(2) ソースに残す opt-in (limiterEpsConst、FORGE_DIAG_PSI_DUALEVAL、assembleResidual の前半/後半分割) が既定経路の挙動を変えていないか (diff を確認)、残すべきか削るべきか (limiter-config-simplify はキーを増やしたくない方針)、(3) accepted へ移してよいか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
Z として調査を閉じる方針は支持します。ただし「原因確定・解への無影響」への書き換えは認められません。  
一次成果物の照合、結論の限定、不要な試験キーの撤去を済ませてから `accepted` へ移してください。

指定 diff は取得しました。`assembleResidualPre` と `Post` の本体を連結すると、コメント・空白を除いて変更前と完全一致しました。既定 off の診断は組立を追加せず、`limiterEpsConst=0` の ε 計算も元の式です。境界半割面の除外、周期集約、軸の射影順序を変える差分はありません。**ソース上の既定経路の破壊は見つかりませんでしたが、GPU 実行のビット一致は未検証です。**

1. **Major — 一次成果物を確認できず、result レビューの数値検証が成立していません。**

   [case README:394](/home/sano/work/forge-species/case/16.nozzle_wys/README.md:394) には AWS の run が登録されていますが、この checkout と隣接 worktree に対象 run はありません。実際に `check_convergence.py` を実行した結果、次の４本すべてが `NO residual_history.csv`、総合結果は `OVERALL: CHECK FAILURES ABOVE` でした。

   - `case/16.nozzle_wys/run_0527_liminlet_ctrl_cfl2/`
   - `case/16.nozzle_wys/run_0528_liminlet_cfl1/`
   - `case/16.nozzle_wys/run_0536_limeps_ctrl/`
   - `case/16.nozzle_wys/run_0537_limeps_const2e-6/`

   これは計算の失敗判定ではなく、**レビュー環境での証拠不足**です。`run_0535` の D・Q・復元検査も CSV を直接確認できていません。以下の数値指摘は plan と台帳の照合に基づきます。

   **対案:** config・ビルド識別子・全残差履歴・両判定ツールの出力と実行引数・`wall_series.csv`・`psi_dualeval.csv`・入口残差の両末尾窓集計を、レビュー可能な場所へ揃えてください。準定常判定の既定閾値は drift 5%／変動10%なので、`ALL STEADY` だけでは今回の壁圧0.1%／壁温0.1 Kという比較精度を保証しません。[check_quasisteady.py:542](/home/sano/work/forge-species/solver_density_cuda/tools/check_quasisteady.py:542)

2. **Major — 壁温差を20倍過小に要約し、測った範囲を「解全体」へ拡張しています。**

   [plan:148](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:148) と [methods/limiter.md:137](/home/sano/work/forge-species/methods/limiter.md:137) は CFL 半減も含めて壁温差を `≤1e-4 K` としています。しかし [plan #2r:179](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:179) と台帳の CFL 比較は **−0.002 K** です。

   また [wall_series.py:15](/home/sano/work/forge-species/notes/investigations/2026-10-03-limiter-inlet/wall_series.py:15) が抽出するのは上側壁で、温度は指定測点の平均です。壁温分布全体・入口局所場・保存収支を検証した結果ではありません。

   **対案:** CFL 比較と ε 比較を分けて数値を記載し、「測定した上側壁の４報告量は設定した許容差内」としてください。`0.002 K` でも事前許容差 `0.1 K` は満たします。結論を強めるために数値を縮める必要はありません。「解は不変」「影響が無い」は撤回すべきです。

3. **Major — 原因、ε の有効範囲、収束状態について、観測を超えた断定があります。**

   [methods/limiter.md:133](/home/sano/work/forge-species/methods/limiter.md:133) は丸めによる極値・面選択の切り替えを床の原因として説明します。しかし [plan #5er:187](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:187) の登録判定は **D=0.915/0.914で保留**。Q≈0.704 は特定状態・特定保存ψへの差し替え効果であり、長期凍結後の残差床ではありません。

   [plan #8r:190](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:190) で試した定数 ε̂ は **2e-6 の一点**です。「精度を守れる範囲では消せない」という範囲全体の結論は出せません。同じ行の「入口以外の残差だけが下がった」も、入口ノルム比0.73〜0.80という記載と矛盾します。

   さらに [methods/limiter.md:141](/home/sano/work/forge-species/methods/limiter.md:141) の「床の所在を見れば真の未収束と区別できる」は不適切です。ツールの停滞判定は **`NOT CONVERGED (stalled/plateau …)`** です。[check_convergence.py:408](/home/sano/work/forge-species/solver_density_cuda/tools/check_convergence.py:408)

   **対案:** 「入口近傍のリミッタ活動と残差床の関連を観測したが、原因は未確定。試した ε̂=2e-6 は採用基準未達。未収束の準定常評価として、限定した報告量への感度が小さいため追加調査を打ち切る」と統一してください。残差位置を合格条件にしてはいけません。

4. **Major — 不採用の `limiterEpsConst` を公開キーとして残す根拠が不足しています。**

   [limiter-config-simplify:16](/home/sano/work/forge-species/plans/active/limiter-config-simplify.md:16) の目的は設定表面の縮小です。本試験は不採用で、標準ケース・格子精度検証も未実施。それでも [plan:225](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:225) はキーを残すとしています。`methods/`・`procedures/` にキーの仕様もありません。

   実装にも公開設定としての穴があります。[solverConfig.cpp:642](/home/sano/work/forge-species/solver_density_cuda/input/solverConfig.cpp:642) は非有限値を拒否しません。[limiter_d.cu:60](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:60) の二乗と float32 変換を再現すると、入力 `1e-30` は `-0.0`、`1e20` は `-Inf` になります。また `limiterScaled:1` の検査後に cell 経路では0へ落ちるため、要求した介入が無効になる組み合わせも受理されます。[solverConfig.cpp:681](/home/sano/work/forge-species/solver_density_cuda/input/solverConfig.cpp:681)

   **対案:** **`limiterEpsConst` と定数 ε の分岐を撤去**し、実験再現は commit `65056017` と成果物で保持してください。診断 `FORGE_DIAG_PSI_DUALEVAL` と前半／後半分割は残してよいと判断します。前者は通常作用素を更新に使う診断で、周期 node 拒否・実ノード限定上書きも確認できました。試験キーの恒久化とは用途が異なります。

5. **Minor — 残作業表と完了条件が、Z で閉じる判断に追従していません。**

   [plan:195](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:195) の `#5c` は未完了で、`#2r` の入口ノルムは末尾10%／20%を分けて記録していません。`#3` は終了窓番号が未確定のまま、`#6` は旧 §4.4 参照です。§8 は依然として原因特定または「どれでもない」を要求し、**保留のまま調査を終了する今回の判断を表現していません**。[plan:221](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:221)

   §7 と plan レビュー免除理由も「ソース変更は診断だけ」のままです。構造 lint は `PASS` でしたが、内容の整合性は満たしません。

   **対案:** 各行を完了・中止・未解決として整理し、原因未確定を残してください。`#5c` は既存データの再集計で埋め、§8には「原因特定による完了ではなく、限定的な運用判断による調査終了」と記録します。移動時に `plans/README.md` と参照リンクを同期してください。`methods/index.md` の既存索引は維持されています。

6. **Minor — van Albada の逸脱量に係数２の誤りがあります。**

   [plan:151](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:151) の式 `½(r+r²)/(1+r²)` を `r=1−√2` で評価すると **−0.103553** です。本文の「約10%」は正しい一方、括弧内の最小値 `−0.207 δ_ij` は誤りです。

   **対案:** `−0.103553 δ_ij` に訂正してください。約10%の逸脱という判断自体は変わりません。

推奨は、**原因未確定の限定的な調査終了としてZを採用すること**です。移動前の優先順は、①一次成果物の照合、②数値・因果・収束表現の訂正、③`limiterEpsConst` の撤去、④残作業表・完了条件・リンクの同期です。追加のリミッタ開発や長期凍結試験は、この閉鎖の条件にしません。

ファイル変更は行っていません。本レビューの採否は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
