# codex 諮問 (diagnose): orchestrator-v3

- **brief**: [`notes/reviews/briefs/2026-10-18-orchestrator-v3.md`](../../notes/reviews/briefs/2026-10-18-orchestrator-v3.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-18
- **commit**: `1ef4b7d2` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 5.4 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-18-orchestrator-v3-diagnose-v3.codex.log` (codex: ok)
- **結論**: まず §6 V3・§5.1 #13 の事前登録を上記の実行順序・判定経路・網羅性・見積もり契約で修正し、13a を先行、13b を費用承認待ちとして分ける。
- **extra**: `design/forge_design/campaign/recipes/wt_card_chain.py`, `design/forge_design/campaign/budget.py`, `design/forge_design/campaign/driver.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| Q1：N1〜N3・S2 の実行順序 | **要再検証・Major** | `approval.preview` は問題 YAML を必要とし、無ければ問題を返す（[approval.py:362](/home/sano/work/forge-cardmesh/design/forge_design/intake/approval.py:362)）。したがって「承認 → 初めて `wt_card`」では成立しない。**承認前に用いる問題 YAML の準備方法を明記**する。固定 fixture を使うなら、そのハッシュと、承認後に生成した問題との `approval.check` を登録する。S2 は同じ正常 fixture を与え、`wall` の出所欠落を理由に承認できず、後続の生成に進まないことを検査する。「問題 YAML が無いから止まった」を成功に数えない。 |
| Q1：S1 の公差未達 | **要再検証・Major** | 実物の `case/90.qs_trial_a/run_0003_ns_p0/` に、現在の条件を使って `AP.check_run` を読み取りのみで実行すると、`MISMATCH`、差分は **`quantities.residual_finite.required (今の値に無い)`**。既存の N15 と同じ結果（[run_approval_tests.py:185](/home/sano/work/forge-cardmesh/design/tests/run_approval_tests.py:185)）。正式な承認モードでは最終結果が `incomplete` になる（[wt_qualify.py:1177](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:1177)）。**S1 は旧データの非承認モードによる公差判定の部品試験と明記**し、承認不一致の負例と分ける。既存 Q2 も非承認モードである（[run_wt_qualify_rv_tests.py:117](/home/sano/work/forge-cardmesh/design/tests/run_wt_qualify_rv_tests.py:117)）。「厳しく」ではなく公差・評価窓・反復番号・上限を数値で固定し、`back_to_correction`／`rejected` と、その原因まで照合する。 |
| Q1：既存試験の集計だけで負例を網羅 | **要再検証・Major** | `run_driver_tests.py` は systemd が使えないと**試験を飛ばして終了コード 0**を返す（[run_driver_tests.py:54](/home/sano/work/forge-cardmesh/design/tests/run_driver_tests.py:54)）。プロセスの成功だけでは負例を実行した証拠にならない。**§6 V3 の各条件 → 個別の検査 ID → 期待する拒否理由**の対応表を固定し、必須検査の未実行・欠落・skip は合格にしない。既存検査の再利用は採用するが、「新しい検査は書かない」は網羅性を確認してから決める。 |
| Q3：各段の見積もりを起動前チェックへ接続 | **採用。ただし数値の一律適用は却下・Major** | 空の見積もりを返す点は事実（[wt_card_chain.py:203](/home/sano/work/forge-cardmesh/design/forge_design/campaign/recipes/wt_card_chain.py:203)）。残予算 0.01 GPU h に対し、`check_start({}, estimate={})` 相当の確認は通り、見積もりを 4/60 GPU h にすると拒否した。既存 driver は検査を試行記録・起動より前に行っている（[driver.py:251](/home/sano/work/forge-cardmesh/design/forge_design/campaign/driver.py:251)）。**依頼・実行ホストに対応する段別見積もりを台帳に固定して渡す**。根拠 run、格子規模、計算条件、出力条件を記録し、ディスクは既存分を除く増分とする。case/71 の時間を別条件へ無条件に流用しない。これは既存の予算契約に沿う（[plan:368](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:368)）。 |
| Q2：13a／13b の分割 | **採用** | 機械検査を先行し、モデル実行を費用承認待ちにする分割は妥当。ただし **13a PASS と V3 完了を分ける**。モデル版・役割・介入条件・各3回・機械と報告それぞれの誤合格集計は V3 の要求として残る（[plan:641](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:641)）。13b 未実施の間は V3 と §8 の完了を宣言しない。 |
| 事前登録の固定条件 | **採用に必要な補足・Minor** | 現案の「厳しく」「既定」「既存データ」だけでは再現条件が確定しない。[plan:609](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:609) に従い、カード・問題・判定器／依存環境・参照場と判定入力のハッシュ、正確な公差と判定区間、外部データの所在を登録する。N3 は格子部品の成功として扱う（[mesh_stage.py:17](/home/sano/work/forge-cardmesh/design/forge_design/campaign/mesh_stage.py:17)）。 |

結論: **まず §6 V3・§5.1 #13 の事前登録を上記の実行順序・判定経路・網羅性・見積もり契約で修正し、13a を先行、13b を費用承認待ちとして分ける。**

第 1 仮説: S3 が予算不足で起動前に止まらない主因は、recipe が段の見積もりを渡していないことである。 確度: 高  
根拠: `wt_card_chain.py:203` は `estimate: {}`。`budget.py:66` は欠けた見積もりを 0 として比較する。読み取りのみの関数確認では、残予算 0.01 GPU h に空の見積もりは `ok: true`、0.0666667 GPU h の見積もりは `ok: false` だった。  
反証条件: 有効な見積もりが残予算を超えている状態でも、実際の driver が Euler の試行を記録・起動する。

第 2 仮説: S1 を現在の承認モードで実施すると、公差未達とは別の承認不一致により、期待した `back_to_correction` に到達しない。 確度: 高。実物 run の `check_run` で前記の不一致を確認した。

判別 A/B: **変えるのは Euler 段の `estimate.gpu_hours` だけ**。同じ mesh 完了状態・残予算 0.01 GPU h から、A は現状の未指定、B は試験用の 4/60 GPU h とする。ディスク条件等は同一にし、各1回、driver の起動判定まで実行する。起動処理を到達回数だけ記録する代替にして、GPU・forge は起動しない。見る量は Euler 試行の追加数、起動到達数、停止理由、既存台帳の保持。  
→ **A は追加1・到達1、B は追加0・到達0・予算理由で停止**なら、見積もり未接続の仮説を支持する。B も到達すれば受け渡し／検査にも欠陥がある。A も止まれば対照不成立。この値は試験入力であり、生産用見積もりの妥当性を証明しない。

やらない方がよいこと: S1 の旧基準を合格目的で書き換えること、例外や任意の停止を負例成功に数えること、case/71 の時間を全依頼に固定すること、13a の結果だけでモデル運用まで検証済みとすること。

呼び出し側の前提への異議: 「GPU なしでは予算の停止を検証できない」は受け入れない。実行中の予算停止も既存 T7 がダミーの `sleep` で検査している（[run_driver_tests.py:195](/home/sano/work/forge-cardmesh/design/tests/run_driver_tests.py:195)）。問題は GPU の有無ではなく、本物の recipe から見積もりが接続されていないこと。また「試行を記録しない」は **Euler 試行を追加しない**という意味に限定し、先行する mesh 試行と予算停止の記録は保持する。

不足情報: 生産用見積もりの根拠 run・ホスト・ディスク増分、および13bの実行モデルと費用上限は登録時に必要。今回の設計判断は妨げない。ファイル変更・forge 起動・試験一式の再実行はしていない。**plan 未反映**。反映先は `tooling-design-orchestrator.md` §4.5.1・§5.1 #13・§6 V3・§8。

証拠: 足りる
