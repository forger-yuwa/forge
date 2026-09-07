# <タスク名>

<!--
配置: 新規計画は `plans/active/<area>-<short-slug>.md` に作る (例: `gradient-least-squares.md`)。
ファイル名は固定し、ライフサイクルはフォルダ移動で表す:
  active/ (検討中・進行中) → accepted/ (実装・検証済みで現役) / archived/ (superseded・終了)。
-->

## メタ

- **area**: `<gradient | limiter | convection | diffusion | time_integration | boundary | poisson | architecture | その他>`
- **status**: `draft`  <!-- draft / in_progress / done / superseded。done で accepted/、superseded で archived/ へ移動 -->
- **related_docs**:
  - `methods/<area>/` (現在仕様ドキュメント)
- **related_plans**: <!-- 親 / 子 / 関連計画があれば列挙、無ければ削除 -->
- **created**: `YYYY-MM-DD`
- **owner**: `<担当>`

## 1. 目的

このタスクが解決する課題と、完了時に得られる状態を 2〜4 行で記述する。

## 2. スコープ

- **やる**: 本計画で実装する範囲を箇条書き。
- **やらない**: 将来計画 / 別 plan に切り出す範囲を明示。

## 3. 関連 docs と前提

理論的背景・既存実装の参照先 (`methods/<area>/`) と、依存する別計画を列挙する。
理論や実装方針に変更が及ぶ場合は、本計画の実装に入る前に対応する `methods/<area>/`
の現在仕様を先に更新する。

## 4. 設計方針

数式・データ構造・呼び出し関係など、実装に直結する設計判断を記述する。
スキーム名や物理プロセスはここでは `methods/<area>/` の現在仕様へリンクし、
本計画には差分・判断のみ書く。

## 5. 実装ステップ

1. <ステップ 1: 触るファイルと変更概要>
2. <ステップ 2>
3. <ステップ 3>

各ステップで触る主要ファイルを併記すること。

### 5.1 残作業 (優先順)

**残作業の正本はこの表**。`notes/sessions/` の引き継ぎ文書には写しとポインタだけを置く。
方針が決まった時点で (実装の有無を問わず) ここを更新する。

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | <最優先の残作業> | <何をどこで> |
| 2 | <次> | <> |

## 6. 検証

- **単体 / ビルド**: コンパイル・ユニット検証手順。
- **検証ケース**: `case/<case_name>/` のどれを使うか、[`../procedures/verification/README.md`](../procedures/verification/README.md) 参照。
- **判定基準**: 残差・物理量・性能の合否ライン。

## 7. 影響範囲

- 触るモジュール / ファイル一覧
- 既存ケース・実行手順への影響
- ドキュメント (`methods/<area>/`, `methods/index.md`) の更新箇所

## 8. 完了条件

- [ ] 関連 `methods/<area>/` の現在仕様を更新済み
- [ ] 実装・検証完了 (本計画の §6 を満たす)
- [ ] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [ ] ファイルを `plans/active/` → `plans/accepted/` (superseded なら `archived/`) へ移動
- [ ] [`plans/README.md`](README.md) の一覧を同期 (移動元・移動先)

## 9. 変更ログ

- `YYYY-MM-DD` — 初稿。
