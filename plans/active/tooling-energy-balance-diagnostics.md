# エネルギー収支診断の出力 (拘束前残差と境界数値流束)

## メタ

- **area**: `tooling / architecture (output)`
- **status**: `draft`
- **related_docs**:
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md) — 出力登録と時間積分の現在仕様
  - [`procedures/solver-settings.md`](../../procedures/solver-settings.md) — `output.level` / `extraFields`
- **related_plans**:
  - [`case-hypersonic-gap-heating-validation.md`](case-hypersonic-gap-heating-validation.md) — **発注元**。すきま壁の熱量を離散収支で裏取りしたいが、現状は取得経路が無い (あちらの §4.8)
- **created**: `2026-09-19`
- **owner**: `sano`

## 1. 目的

壁熱量や領域熱収支を**離散スキームの言葉で**検算できるようにする。現状、
[`nodeWallDirichlet_d.cu`](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu) が等温壁ノードの
エネルギー残差をゼロ化し、拘束前の `res_roe` とエネルギー数値流束は
[`variables.hpp`](../../solver_density_cuda/variables.hpp) に未登録なので
([`output.cpp`:55](../../solver_density_cuda/output/output.cpp) の仕様上 `extraFields` でも出ない)、
**場から積分した物理収支しか取れない**。これは SLAU の数値散逸・再構成・node の拘束操作を含まないので、
「壁熱量が離散的にも閉じているか」を答えられない。

## 2. スコープ

- **やる**: (a) 拘束前 `res_roe` の出力、(b) CV 境界を横切るエネルギー数値流束の出力、
  (c) 符号・単位・対象 CV の定義、(d) 解除試験 (1D 熱伝導で機械精度で閉じること)
- **やらない**: 収支に基づく自動判定ツール、他方程式 (運動量・化学種) への拡張、cell 離散化での同等機能

## 3. 関連 docs と前提

- 出力登録は `variables.hpp` の `output_cellValNames` に無いと `extraFields` でも出ない (警告のみ)。
- node の壁は Dirichlet 拘束でエネルギー残差を潰すため、**拘束前**の値が要る。

## 4. 設計方針

- 診断は **opt-in** (既定 OFF、`output.extraFields` で拾う形) とし、既定の出力量・性能を変えない。
- 対象 CV 群 (壁列、壁から 1 層内側、開口断面で囲まれた領域) を指定できるようにする。
- 数値流束は**面単位**で集計し、CV 境界の向き (外向き正) と単位 [W] を明示する。

## 5. 実装ステップ

1. `variables.hpp` への登録 (`res_roe_raw`, `eflux_face` 等) と `output.cpp` の経路確認。
2. 拘束前 `res_roe` の退避 (node 拘束の直前)。
3. 面エネルギー流束の集計と出力。
4. 解除試験の追加。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | 要求の確定 | 発注元 ([case-hypersonic-gap-heating-validation](case-hypersonic-gap-heating-validation.md) §4.8) と対象 CV・量・単位をすり合わせる |
| 2 | 出力登録 | `variables.hpp` / `output.cpp` の経路。既定挙動を変えない |
| 3 | 拘束前残差の退避 | `nodeWallDirichlet` 直前の値 |
| 4 | 面流束の集計 | 外向き正・[W]・面単位 |
| 5 | 解除試験 | 1D 熱伝導で収支が機械精度で閉じる。これが通って初めて発注元の離散収支評価を開始できる |
| 6 | codex レビュー | `plan` 段 (§4 が書けた時点)、`result` 段 |

## 6. 検証

- **単体**: 1D 熱伝導 (解析解あり) で、CV 境界の数値流束和 = 内部エネルギー変化率 が機械精度で閉じる。
- **回帰**: 診断 OFF で既存ケースの場・残差・速度が不変。
- **判定基準**: 収支残差が保存量の丸め誤差程度 (double 相当) / float 経路では丸め見積り以内。

### 6.1 レビュー記録 (codex)

<!-- 実装着手前に `--stage plan` を回し、ここに 1 行 (段階/日付/記録パス/判定/対応) を入れる。 -->

## 7. 影響範囲

- `solver_density_cuda/variables.hpp`, `output/output.cpp`, `cuda_forge/nodeWallDirichlet_d.cu` 周辺。
- 既定 OFF なので既存ケースへの影響は無い想定 (回帰で確認)。

## 8. 完了条件

- [ ] 解除試験 (1D 熱伝導) が機械精度で閉じる
- [ ] 診断 OFF の回帰が不変
- [ ] codex レビュー (`plan` / `result`) を §6.1 に記録
- [ ] `status: done` にし `plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-09-19` — 初稿。[case-hypersonic-gap-heating-validation](case-hypersonic-gap-heating-validation.md) の
  codex plan レビュー (2–4 巡目) で「離散収支の取得経路が無い」と指摘されたのを受けて起票。
