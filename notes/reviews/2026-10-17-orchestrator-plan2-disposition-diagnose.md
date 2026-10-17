# codex 諮問 (diagnose): orchestrator-plan2-disposition

- **brief**: [`notes/reviews/briefs/2026-10-17-orchestrator-plan2-disposition.md`](../../notes/reviews/briefs/2026-10-17-orchestrator-plan2-disposition.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-17
- **commit**: `9436b084` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.4 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-orchestrator-plan2-disposition-diagnose.codex.log` (codex: ok)
- **結論**: 5 件とも採用し、まず §4.5 を「外部の service 管理が attempt 全体の停止を保証する契約」に確定して、GPU 不要の障害 A/B を実装の先頭に置いてください。
- **extra**: `notes/reviews/2026-10-17-tooling-design-orchestrator-plan-2.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 指摘 | 採否 | 根拠と反映の評価・対案 |
|---|---|---|
| **Major M1** 実行担当の障害 | **採用。ただし反映は未完了** | [stop_rule.py:1901](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1901) は別セッションで子を起動し、例外時の停止も同じ監視プロセスが実行します（同:1980）。**監視自身への `SIGKILL` では、その停止処理は動きません。** [plan:250](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:250) の「仕組みは実装で決める」は広すぎます。停止責任を実行担当の外に置き、対象・期限・停止確認・予算精算まで plan で決めてください。 |
| **Major M2** 節点数の承認流用 | **採用** | [hand_mesh.py:354](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_mesh/hand_mesh.py:354) の固定承認文と点数からの上限拡大を確認しました。§4.6・#8・V2c の修正は適切です。ただし [plan:271](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:271) の「そのまま包む」は「格子生成の数値処理を再利用し、承認・上限生成部分は #8 で置換」に直すべきです。Euler・NS の両方を照合対象にしてください。 |
| **Major M3** V1b の期待結果 | **採用** | [WT_QUALIFICATION.json:50](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/run_0014_chain_v2/wt_qualify_provenance/WT_QUALIFICATION.json:50) は実際に `accepted / done`。V1b の修正は妥当です。「補正へ戻す」は**試走 A の `run_0003_ns_p0` に固定**してください。根拠は下記。 |
| **Major M4** 完了条件 | **採用** | [plan:382](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:382) と [plan:410](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:410) により、途中介入なしの通し成功が必須になりました。指摘の意図を満たしています。V4(a) では未較正カードと入力参照を固定し、既存の補正済み壁を途中投入して合格する余地をなくしてください。 |
| **Minor m5** run 索引 | **採用。同期条件を補足** | [plan:302](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:302) は作成・終了時を追加しましたが、[AGENTS.md:62](/home/sano/work/forge-cardmesh/AGENTS.md:62) の**破棄時**も必要です。また実行担当が死亡すると終了時更新は実行できません。再開時に台帳から修復し、複数 campaign が同じ README を更新する場合は case 単位のロックと原子的更新で他行を保持してください。 |

M1 は、**attempt ごとの `systemd` service に実行担当と全子孫を収容する方式**を推奨します。実行担当を main process とし、`ExitType=main`、`KillMode=control-group`、有限の停止期限、`Restart=no` を契約にします。main process の終了で service を停止し、残ったプロセス全体を停止する仕組みです。[公式 service 定義](https://raw.githubusercontent.com/systemd/systemd/main/man/systemd.service.xml)、[公式 KillMode 定義](https://raw.githubusercontent.com/systemd/systemd/main/man/systemd.kill.xml)

この方式なら `start_new_session=True` を維持できます。§4.5 の `nohup` は service 起動に置き換えます。**死んだ実行担当自身に後始末を期待する構成にはしない**でください。停止期限内に全子孫の終了を確認できなければ新 attempt を禁止し、確認できるまで予算予約を解放しない契約も必要です。API の呼び方は実装に任せられますが、この責任分担は plan で決める事項です。

fixture の選択根拠は以下です。どちらも `check_settled.py --end-step … --no-write` を今回再実行しています。

| 用途・run | 判定区間と VERDICT | 確認した値 |
|---|---|---|
| V1b 正例：`case/71.m55_combustion_wt/run_0014_chain_v2/` | 36,000〜46,000、停止観測 46,171。**`VERDICT: SETTLED`** | 出口コア M＝5.502841913。V5c 最大偏差＝0.062305 %、片振幅＝0.057285 %。保存 qualification は `accepted`。 |
| 補正への分岐：`case/90.qs_trial_a/run_0003_ns_p0/` | 32,000〜42,000、停止観測 42,060。**`VERDICT: SETTLED`** | 出口コア M＝5.519563193、要求 5.5 から **+0.355694 %**。V5c 最大偏差＝**0.404298 %**、片振幅＝0.030797 %。公差 0.1 % では V5b と V5c 最大偏差が FAIL。 |

両 run とも全残差履歴の `rms_*` 14 列と停止時の `VALUE/*` 23 dataset の非有限数は **0**。保存メッシュ判定は両方 **`VERDICT: PASS (AR<= 5000, skew<= 0.90)`**。試走 A の保存壁解像判定も **PASS** です。監視対象の出口コア M・線上量について `SETTLED` を根拠にしており、単一場だけから定常性を推定していません。

試走 A は `card_p0.yaml`、`res_42060.h5`、公差、評価窓、反復上限を固定してください。既存の [decide:877](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:877) に確認済み判定を渡す単体確認では、`back_to_correction / escalate` でした。**`next=recalibrate` には較正规則の `apply=true` が別途必要**です（同:908）。qualification 全体を今回再実行したという意味ではありません。

case/71 の初回壁は、[README:26](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:26) に未整定・ディスク満杯、延長も CENSORED と記録されています。「解析は成立しているが公差未達」という分岐の fixture には試走 A の方が適切です。恒久索引は [case/90 README:6](/home/sano/work/forge-cardmesh/case/90.qs_trial_a/README.md:6) と [case/71 README:38](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:38) です。

結論: **5 件とも採用し、まず §4.5 を「外部の service 管理が attempt 全体の停止を保証する契約」に確定して、GPU 不要の障害 A/B を実装の先頭に置いてください。**

第 1 仮説: 現在案に残る最大の欠落は、実行担当死亡後の停止責任が未定義なことです。 確度: **高**  
　根拠: `stop_rule.py:1901`・`:1980`、`plan:250`。同じ `start_new_session=True` の構成をダミー子で確認したところ、親の `SIGKILL` 後も子は状態 `S` で生存しました。試験用の子は後始末済みです。  
　反証条件: 実行担当を `SIGKILL` しても、担当外の仕組みにより別セッションの全子孫が規定期限内に停止し、二重起動防止・予算記録も維持されること。

第 2・第 3 仮説: 追加なし。

判別 A/B: **変更点は実行担当への `SIGKILL` の有無だけ**。提案した service 内で、別セッションのダミー `run_case` と、その子のダミー forge を動かし、100 ms ごとに生存を観測します。A は無介入で 10 秒、B は 2 秒時点で実行担当だけを強制終了。停止期限は試験前に 5 秒へ固定します。→ **A が継続し B の全子孫が期限内に消えれば契約の成立を支持／B が残れば M1 未解決**。A も止まる場合は試験不成立です。forge・GPU は使いません。

やらない方がよいこと: `start_new_session=False` や親の `finally` だけで死亡時停止を保証したことにする、既存 runner を再呼び出す、case/71 初回壁の未整定を公差未達の正例に使うこと。

呼び出し側の前提への異議: **指摘の全件採用と、反映完了は別です。** M1 の仕組みは全面的に実装へ委ねられません。また「補正へ戻す」と「自動較正を続行する」は現行コード上で別判定です。

不足情報: 対象 AWS・手元環境での service 管理の利用可否は未確認です。実装開始時の事前検査にし、利用できなければ `nohup` へ黙って戻さず起動を止めてください。

証拠: 足りる
