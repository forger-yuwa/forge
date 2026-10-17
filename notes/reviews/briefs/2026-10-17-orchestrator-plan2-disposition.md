# 諮問ブリーフ: orchestrator plan の plan 段レビュー (2 回目) の Major 4 件・Minor 1 件の採否と反映 (2026-10-17)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件 5: codex の Major の採否を決める。
- plan: `plans/active/tooling-design-orchestrator.md` (作業ツリーの未 commit の変更に反映済み。`git diff` で見られる)
- レビュー: `notes/reviews/2026-10-17-tooling-design-orchestrator-plan-2.md` (GO-with-changes、C0/M4/m1)
- commit: 9436b084 + 作業ツリーの変更

## 観測事実 (主セッションが確かめた根拠)

- Major 1: `design/forge_design/evaluate/stop_rule.py:1901` は `run_case.sh` を `start_new_session=True` で起動し、`Popen` と監視の状態を呼び出し元が持つ。
  `:1879` は本段の出力が既にあると起動を拒む。
- Major 2: `case/71.m55_combustion_wt/_mesh/hand_mesh.py:354` 付近は節点数から `max_nodes` を切り上げ、理由に固定の文「ユーザ決定「a」2026-10-15: …」を書く。
  `runner_axismach.py` の `node_budget` (790 行付近) は正の整数・空でない理由だけを見る。規則は `plans/accepted/tooling-nozzle-mesh-size-guard.md` §4.3。
- Major 3: `case/71.m55_combustion_wt/run_0014_chain_v2/wt_qualify_provenance/WT_QUALIFICATION.json` は `result: accepted`・`next: done`、
  V5b 5.502842 (範囲 5.4945〜5.5055)、V5c 最大の偏差 0.0623 %・片振幅 0.0573 % (公差 0.1 %)、壁解像 PASS、flow PASS。plan の V1b は「補正へ戻す」が正解と書いていた (誤り)。
- Major 4: 前の §8 は V4 を「詰まる所の一覧と対処」で完了としていた。V1b は前段を固定した 1 pass、V3 は格子まで。
- Minor 5: AGENTS.md「case の run 索引」は run の作成・破棄で case README の run 一覧の同期を必須にしている。

## 採否 (案、全件採用) と反映先

| 指摘 | 採否 | 反映 |
| --- | --- | --- |
| M1 実行器の障害時の監視・予算・再接続 | 採用 | §4.5 に「試行ごとの実行担当」(runner・停止規則の監視・予算の監視を 1 つのプロセスが持つ; CLI は接続だけ; 実行担当が死んだら子も止める; 死んだ試行は再接続せず完了の証拠を確かめて新しい試行へ; 既存の runner を再び呼ばない)。§5.1 #6 の合格に障害の試験 3 つ |
| M2 格子の節点数の承認の流用 | 採用 | §4.6 の契約に「依頼ごとの承認済みの上限と承認元」、必要な点数と上限を分けて格子の前に照合。§4.7 のユーザ欄に上限。§5.1 #8 (#10 を待たない)、§6 V2c (別のカード + 過去の承認の文・承認を超える作り直し の負例) |
| M3 V1b の期待の結果 | 採用 | V1b は `accepted` + 成果物。「補正へ戻す」は別の fixture (候補: 試走 A `case/90.qs_trial_a` か case/71 の 1 回目の壁; 事前登録で決める) |
| M4 完了条件 | 採用 | §6 V4 を (a) 決定的な recipe の通し (成功の条件: 既知の正常のカード 1 件が入口から必須の成果物と最終の `accepted` まで、人もモデルも直さない) と (b) モデルの測定 (記録) に分けた。§8 に V4 (a) を必須に |
| m5 case の run 一覧の同期 | 採用 | §4.8 と §5.1 #6 (台帳からその campaign の行を生成、他の行は残す) |

## 問い

1. 全件採用でよいか。反映 (`git diff plans/active/tooling-design-orchestrator.md`) が指摘の意図を外している所、または新しい矛盾を作った所はあるか。
2. M1 の「実行担当が死んだら子も止める」は、停止規則 v2 が子を別のセッションで起動している今の作りと両立するか。仕組みを実装に任せてよいか、plan で決めておくべきか。
3. M3 の「補正へ戻す」の fixture は、試走 A と case/71 の 1 回目の壁のどちらが良いか (あるいは別のもの)。

## 読んでよいファイル

- `plans/active/tooling-design-orchestrator.md` (作業ツリー)、`git diff` (同)
- `notes/reviews/2026-10-17-tooling-design-orchestrator-plan-2.md`
- `design/forge_design/evaluate/stop_rule.py` (1860〜1990 行)、`case/71.m55_combustion_wt/_mesh/hand_mesh.py`、`design/forge_design/evaluate/runner_axismach.py` (780〜820 行)
- `case/71.m55_combustion_wt/README.md`、`case/71.m55_combustion_wt/run_0014_chain_v2/wt_qualify_provenance/`
- `case/90.qs_trial_a/` (試走 A、git 追跡外)、`case/71.m55_combustion_wt/_runs/qs_trial_a_evidence/`
