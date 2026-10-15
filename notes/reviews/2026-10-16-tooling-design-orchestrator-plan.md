# codex レビュー: tooling-design-orchestrator (plan)

- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **stage**: `plan`
- **date**: 2026-10-16
- **commit**: `633b3fee` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-tooling-design-orchestrator-plan.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **extra**: `plans/active/tooling-design-problem-campaign-recipe.md`, `../forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

既存の campaign・判定ツールを再利用し、判断役を限定された選択に閉じ込める方向は妥当です。
ただし、実装の前提、設計完了の判定、再開の契約、V1a の合格条件に修正が必要です。
以下の Major を優先順に計画へ反映してから実装してください。

1. **Major — 親計画の未実装の契約が、実行器の前提から抜けている。**

   **根拠:** 本計画は親の「ステップ5の実装」としていますが、親の [スキーマ・adapter・`resolve_evaluation`](/home/sano/work/forge-cardmesh/plans/active/tooling-design-problem-campaign-recipe.md:181) は未完了です。コード検索でも `resolve_evaluation`・`ExecutionPlan`・新スキーマの実装はありません。一方、子の [残作業表](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:191) は gate 統合から実行器へ進み、この依存を扱っていません。

   これは名称だけの問題ではありません。親が要求する snapshot の固定に対し、現行の [`prepare_ns`](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2511) は `ic_from` ディレクトリの最大番号の場を選びます。参照 run に出力が増えると、同じ入力指定でも初期場が変わります。

   **対案:** search 無しの最小 campaign、入力参照の解決、pass ごとの `execution_plan.json`、全段の設定照合を実行器より先に実装する項目として追加してください。snapshot・メッシュ・物性の識別を固定し、計画確定後に参照 run が増えても入力が変わらない負例試験を置くべきです。

2. **Major — 計算の成立・受け渡し成功と、最終設計の要求達成が分離されていない。**

   **根拠:** [V1b](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:204) は `REACHED`・関門通過・抽出成功・成果物の実在を確認します。しかし、[`judge_run` の総合判定](/home/sano/work/forge-gate/solver_density_cuda/tools/judge_run.py:1490) に出口 Mach や一様性の要求公差はなく、[`handoff_gate_v2`](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:1388) も整定・保存・壁解像・抽出の成立を扱います。既存の V5b/V5c 判定は [case 側の別処理](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0010_ns_b0/judge_0010.sh:47) です。

   また、[「RUN_VERDICT が SETTLE_VERDICT と矛盾しない」](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:191) だけでは統合の合格条件が曖昧です。`SETTLED` でも壁解像が FAIL なら、総合 FAIL が正しい結果です。

   **対案:** recipe に最終 qualification を明記してください。要求公差・寸法・壁解像・必要成果物をすべて満たした場合だけ `accepted`、評価済みで公差未達なら補正へ戻し、反復上限で `rejected`、必要な判定や成果物が欠ければ `incomplete` とします。公差未達でも整定している例、`SETTLED` でも壁解像 FAIL の例、必須 CAD 欠落の例を試験に加えてください。

3. **Major — `--resume` が中断した処理を安全に扱う契約になっていない。**

   **根拠:** [台帳の再開規則](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:150) は完了段のハッシュ照合だけです。しかし、停止規則は [`start_new_session=True`](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1548) で子プロセスを起動するため、実行器だけが終了しても計算が残り得ます。既存 runner の再呼び出しは [開始・停止記録を削除](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2946) し、段間には [初期場の更新と出力削除](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2981) もあります。単純な再実行では二重起動や証拠の消失につながります。

   **対案:** 台帳の単一書き手ロック、原子的更新、attempt ID、起動プロセスとの結び付け、実行中段の照合手順を定義してください。生存中なら再起動せず監視へ復帰し、停止済みなら完了証拠を確認、不明なら人に返します。再試行は新しい run に分け、回数・累積予算を再開後も保持します。起動直後・結果公開直後の強制終了、同時 `--resume`、古い問いへの回答を使った障害注入試験が必要です。idle 自動停止の未決項目も AWS 対応の前提へ繰り上げてください。

4. **Major — V1a の「GPU 不要・入力一式のバイト同一」は現行実装では成立しない。**

   **根拠:** [V1a](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:203) の比較対象には、生成時刻を含む [`settle_criteria.yaml`](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/settle_presets.py:373)、[`RECIPE.json`](/home/sano/work/forge-gate/solver_density_cuda/tools/recipe_settings.py:859)、[`MAIN_STAGE_STARTED`](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2777) が含まれます。最後の記録は前段計算後に生成され、その初期場のハッシュも含みます。時刻だけを固定しても、GPU 不要で前段の結果まで再生成したことにはなりません。

   **対案:** 決定的な入力の比較と、実行時の来歴照合を分けてください。時刻・配置パスは項目を限定して正規化し、参照先の内容とハッシュの結び付けは別途検査します。GPU 不要試験では固定した前段成果物を使って全段の起動入力を捕捉し、初期場、バイナリ、変換器、数値に効く環境変数も比較対象に含めます。実際の段間引き継ぎと完了記録は V1b で検証してください。

5. **Major — 「ユーザに聞かずに埋めたら止まる」を検査する情報がない。**

   **根拠:** [§4.7 の11](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:173) は聞くべき欄の無断補完を拒否するとしています。しかし、現行の [カード検証](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:277) が受け取るのは値とキーであり、ユーザ指定と AI の推測を識別できません。欄の分類表と必須チェックだけでは、もっともらしい値を埋めたカードが通ります。[V4](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:207) も事後の記録に留まっています。

   **対案:** 要求値と出所を分け、ユーザ入力の参照、recipe の既定、導出、未確認の提案を保存してください。ユーザ指定が必要な項目は、既存の指示との対応または確認済みのカード版に結び付け、未確認なら起動前に止めます。「欠落」だけでなく「値はあるが根拠がない」カードを負例に追加してください。

6. **Minor — 判断役の「全員 digest だけ」が既存 backend の再利用では保証されない。**

   **根拠:** [§4.7 の5・6](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:167) は入力条件を揃えていますが、既存の [`run_codex`](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/codex_review.py:444) はリポジトリを読める実行環境、Claude は [`Read,Grep,Glob`](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/codex_review.py:68) を持ちます。切り替え機構に `local` を足すだけでは比較条件が揃いません。

   **対案:** 外部レビュー用と行動選択用の呼び出しを分け、後者は全 backend で道具なし・同一 digest・同一回答スキーマに限定してください。共用するのは接続・timeout・試行記録の部分に留めます。

7. **Minor — 改訂前の実装指示と影響範囲が残っている。**

   **根拠:** §4.7 は新しい列挙と格子一般化を取り消していますが、[§5](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:178) は両方を実装手順として残しています。[§7](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:219) の「既存の部品は包むだけ」も、`judge_run` の v2 対応や backend 拡張と一致しません。

   **対案:** §2・§4・§5・§7を採用済み方針へ直接書き換え、実装者が優先規則で旧記述を読み替える必要をなくしてください。

**推奨:** search 無しの axismach campaign に絞り、上記1〜5の契約と負例試験を先に確定する進め方を推奨します。case/71 と別カードの case/70 は、この範囲の検証先として妥当です。ソルバ演算を変更しない限り、cell や周期箱の追加計算を本計画の必須条件にする必要はありません。

依頼ファイル内の plan と作業ツリーの plan が同一であることを確認しました。判定は計画と現行コードの照合に基づきます。V6 の完了や計算結果の収束を認定するレビューではありません。ファイル変更・CFD 実行は行っていません。

証拠: 足りる

指摘数: Critical 0 / Major 5 / Minor 2
