# codex レビュー: tooling-sern-te-wake-grid (result)

- **plan**: [`plans/accepted/tooling-sern-te-wake-grid.md`](../../plans/accepted/tooling-sern-te-wake-grid.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-11
- **commit**: `f78b7eef` (feature/sern-design)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
局所化格子の冷点解消と、報告された力係数差は保存時系列で裏付けられます。  
ただし、床判定を迂回する学習条件、生産 YAML の実行長不足、設計点ごとの評価方式 ID の分裂があり、生産受入れの完了は認めません。

1. **Major — 床判定が FAIL／欠測の評価を、新しい学習系列へ混入できる**

   **根拠:** [driver_sern.py:95](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:95) と同ファイルの `_learnable` は、保存された `floor_events` の合否・判定区間を検査しません。[runner_sern.py:458](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:458) の評価方式 ID にも床ゲートの要件は入りません。

   メモリ上で次を再現しました。

   - `floor_events: 1 → 0`、判定窓 `10000 → 1` と変更しても、`eval_method_id` は同一。
   - 現行方式 ID・既存の必須フラグを満たす `status: PASS` 行に、`floor_events: {required: false, verdict: FAIL}` を入れても、`_learnable(...) == True`。

   これは成立可能な行です。[sern_gates.py:265](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:265) は、床ゲートが任意なら床判定 FAIL でも総合合格を許します。§4 の「生産系列ではカウンタを必須入力にし、台帳の採否まで接続」を満たしていません。

   **対案:** キャンペーンの床判定要件を学習条件として照合し、各作動点について PASS・必要な区間長・検査更新数を要求すること。欠測、短窓、FAIL の旧行は再判定まで除外し、この３種類の混入拒否試験を追加してください。

2. **Major — 生産 YAML を標準 runner で実行すると、必ず床判定が判定不能になる**

   **根拠:** 生産２ファイルは `floor_events_window_steps: 10000` に対し、`nStepOuter: 4000` のままです。[lswx08.yaml:94](/home/sano/work/forge-sern-design/case/46.sern_design/problem_3d_prod_3op_wallres_lswx08.yaml:94)

   [runner_sern3d.py:418](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:418) は段階起動後にそのまま `collect` し、検証で使った追加 20000 step を実行しません。`check_floor_events.judge` に同じ長さ条件を与え、次を再現しました。

   > INDETERMINATE: 判定区間 10000 step が最後の session の長さ 4000 step を超える

   また、単に反復数を増やすだけでは検証手順と同等になりません。`run_1097_tewake_lsw10_m6_cont20k` は **GATES PASS** ですが、保存された温度判定は **TRANSIENT-UNSETTLED**、C_M の前窓差は **1.9723e−3 > 5e−4** です。今回は手動で延長されています。

   **対案:** 検証済みの「段階起動 → 20000 step → §6 #1〜#3 の判定 → 必要時のみ +20000」を生産の標準実行経路に接続してください。窓を短くして通す修正は不可です。§5.1 #7 の「済」は、その経路を確認してからにすべきです。

3. **Major — `L_sw` を変えた設計点が別の評価方式になり、同時に学習できない**

   **根拠:** [runner_sern.py:434](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:434) は格子パラメータ全体をハッシュ化し、設計由来として除外するのは２種類の角度だけです。`L_sw` の値はそのまま方式 ID に入ります。

   実測記録でも、同じ g3・局所化・同一メッシャ／変換器なのに、

   - `run_1084`、`L_sw=0.8`: `a50c9185…`
   - `run_1096`、`L_sw=1.0`: `013caf37…`

   と分裂しています。両者の方式レコードで異なるのは `mesh_recipe_sha256` だけです。現行コードでも `L_sw` だけを変更して再現しました。ID の完全一致を要求する [driver_sern.py:146](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:146) では、両点を同じ要求方式で採れません。

   **対案:** 設計変数の値と格子生成方式を分離してください。`L_sw` は設計点・来歴に保持し、方式には `L_sw_exact`、解像度、局所化方式などを残す。そのうえで、今回の `L_sw=0.8/1.0` の２点が同じ系列に入り、旧格子は除外される試験を追加してください。

4. **Minor — §5.1 #6b の「各 run は #1〜#3 を満たす」は実測と矛盾する**

   **根拠:** `run_1087_tewake_m6_A0` の末尾判定窓には、節点 517159 の **約102.53 K** が残っています。`ALL_n_lt150=1`、`ALL_n_lt200=2` です。[TE_MONITOR.csv:41](/home/sano/work/forge-sern-design/notes/investigations/2026-10-11-c6/run_1087_tewake_m6_A0/TE_MONITOR.csv:41)

   したがって「T < 150 K なし」という #2 は不合格です。また、本文の「150 K 未満２節点」は、200 K 未満の件数との取り違えです。旧評価を持ち越さない結論自体は変わりません。

   **対案:** A は「床事象ゼロ・準定常だが低温条件不合格」、修正後 B は「#1〜#3 合格」と分け、件数を訂正してください。

5. **Minor — 現行手順と残作業の記述が整理されていない**

   **根拠:** [plan:102](/home/sano/work/forge-sern-design/plans/accepted/tooling-sern-te-wake-grid.md:102) の同じ行に、g4 への `interp_field` を禁止する判断と、それを実行する旧手順が併存しています。§6 の説明も、現在の g3 の起点 `run_1055` に対して `run_1079` のままです。

   さらに、[plan:165](/home/sano/work/forge-sern-design/plans/accepted/tooling-sern-te-wake-grid.md:165) に識別粒度の未決事項と、EOS 床の commit 問題の「別 plan 未起票」が残っています。一方、2D・ghost 初期化の移管先は存在し、`plans/README.md` に登録されています。

   **対案:** §5.1 を現在の完了・未完了・移管先が分かる表に整理し、旧手順は履歴として明示してください。識別粒度の決着と EOS commit 問題の追跡先を残すこと。また、[methods/design/overview.md:894](/home/sano/work/forge-sern-design/methods/design/overview.md:894) の「g3/g4 で冷点が消え」は、「g3 で消失、g4 の別起動履歴でも冷点なし」に限定すべきです。

保存データの再確認結果は、格子修正そのものを支持します。

| 確認対象 | 今回の確認結果 |
|---|---|
| `run_1086 − run_1084` の末尾窓平均差 | C_T −5.26786e−5、C_T_with_shear −3.567995e−4、C_L −1.98803e−4、C_M +6.34811e−3。plan と整合 |
| 修正後の設計差 | ΔJ = **+1.0659299e−3** を再現。E_J 未確定なので順位確定の根拠にはできない |
| 採用側７ run の温度・力 | 保存 CSV に `check_quasisteady.py --tail 0.5` を再実行し、すべて **ALL STEADY** |
| `run_1084` の床記録 | 圧縮 CSV をメモリ上で再判定し、(10000, 20000] は **PASS、T/ρ/P 各０件** |
| 既定挙動 | 小型格子で変更前版と比較し、blend 0 の座標・接続・境界面はバイト一致 |

保存された収束判定は **NOT CONVERGED** です。特に `run_1096/1098/1099` の `rms_roK` は `falling` であり、「全列プラトー」とも呼べません。今回の評価は準定常の限定比較です。AWS の実入力 HDF5・元の残差 CSV は手元になく、収束判定の再実行と GPU 回帰の独立再現はしていません。

**推奨は、局所化方式を維持したまま `accepted` への移動を止め、指摘1→2→3の順に生産採否経路を修正することです。** その後、実際の生産設定と保存評価を用いた接続試験、指摘4・5の文書修正を確認して再レビューしてください。格子形状の再調整を支持する根拠はありません。

ファイルは変更しておらず、本レビューの指摘は plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 2
