"""設計の実行役の部品 (plan tooling-design-orchestrator §4.1.1・§4.5、親 plan tooling-design-problem-campaign-recipe の限定の実装 #6a・#7a)。

- `ic`: 初期場の固定・継承する値の解決・適合性 (§5.1 #4a)
- `plan`: 実行の計画 `execution_plan.json` (準備の前の固定・準備の後の確定・計画からの実行と照合; §5.1 #4b)
- `campaign`: search の無い campaign `campaign.yaml` (`forge.campaign/v1`; §5.1 #4c)
- `ledger`・`budget`・`stop`: 台帳 (driver が唯一の書き手)・予算の勘定・停止の記録 (§5.1 #6a)
- `step`・`driver`: 試行の実行担当 (service の main process、予算の監視) と、段を順に回す driver (記録 → 起動 → 一度だけの取り込み、再開の照合; §5.1 #6b)
"""
