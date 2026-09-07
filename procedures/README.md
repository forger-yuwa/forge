# forge 運用・開発ハンドブック (`procedures/`)

`procedures/` は forge を**使う/開発するための手順・運用ルール**を置く。
理論や実装の仕様 (`methods/`)、設計判断 (`plans/`)、調査ログ (`notes/`) とは役割が異なり、
ここは「**どう動かすか・どう開発するか**」の how-to を担う。

これらの文書は元々 `.github/` 配下にあったが、GitHub プラットフォーム固有ではない
(Copilot を常用しなくなった) ため `procedures/` に移した。`.github/` は GitHub 設定
(`copilot-instructions.md`、将来の `workflows/` 等) のみとする。

## 計算運用

| 文書 | 役割 |
| --- | --- |
| [`calculation-workflow.md`](calculation-workflow.md) | 計算ケース準備・メッシュ生成/変換・`forge` 実行の標準手順 |
| [`divergence-and-startup.md`](divergence-and-startup.md) | 発散 (NaN/残差爆発) の主因と段階起動の手順。**新規計算・発散時に必ず参照** |
| [`recommended-settings.md`](recommended-settings.md) | **推奨解析設定の正本**: 解析種別ごとの現行レシピ (日付付き) と旧設定 (superseded) 一覧。config を組む前に必ず参照 |
| [`solver-settings.md`](solver-settings.md) | `convMethod` / `limiter` などの数値設定リファレンス (運用上の設定の意味) |
| [`su2-cross-check.md`](su2-cross-check.md) | 同一メッシュ・同一 BC で SU2 と比較し forge 固有の問題を切り分ける手順 |

## 開発環境

| 文書 | 役割 |
| --- | --- |
| [`development-environment.md`](development-environment.md) | 開発環境とビルド (Docker / WSL native)、速度評価の基準環境 |
| [`cloud-aws-gpu.md`](cloud-aws-gpu.md) | AWS EC2 GPU インスタンスでの環境構築・計算投入・速度計測・結果回収の手順 |
| [`coding-conventions.md`](coding-conventions.md) | ソース構成・C++/CUDA 命名規約・ビルド/テスト実行手順 |

## 検証運用

| 文書 | 役割 |
| --- | --- |
| [`verification/README.md`](verification/README.md) | コード変更時の検証ケース選定と確認手順 (索引) |
| [`verification/`](verification/) | 各標準検証ケースの個別手順 (`08-bump.md`, `13-nozzle-h.md`, `20-naca-ml.md` 等) |

ルール全体の正本は [`../AGENTS.md`](../AGENTS.md)。設定値の意味は記憶や推測で判断せず、必ず
[`solver-settings.md`](solver-settings.md) を参照すること。
