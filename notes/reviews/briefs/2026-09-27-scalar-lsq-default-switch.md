# 諮問ブリーフ: node の `mesh.scalarGradient` 既定を gg → lsq に切り替える実装方針 (2026-09-27)

plan: `plans/active/gradient-scalar-lsq-unification.md` §4.4・§5.1 #6 (前提はユーザ決定「B」で差し替え済み)。

## 決めてほしいこと

下の実装案で既定を切り替えてよいか。足りない検証・危険があれば 1 つに絞って。

## 読んでよいファイル

plan 全文、`solver_density_cuda/input/solverConfig.{hpp,cpp}` の `scalarGradient` 周辺 (hpp:545-555、cpp:240-270)、`design/forge_design/evaluate/runner_sern.py:30-60,230-250,930-945`、`design/forge_design/opt/driver_sern.py:40-85,180-240`、`procedures/solver-settings.md` の `mesh.scalarGradient` 節。

## 前提 (事実)

- Phase 1 合格 (2026-09-27): S0/S1、S2 (ユーザ決定で基準を現行 gg 双子に差し替え、lsq/gg 末尾残差 ≤ 1.42 倍・物理量差は上限内)、S3 (lsq/gg 0.997 / 1.018)、FCT smoke の収支差はユーザ決定で揺れの範囲として合格。codex result 2 回は NO-GO だったが指摘は手続き・ツールで、全件対応済み。
- ユーザ決定「B」: 起動記録の結び付け (#2g) を待たない。省略時の警告 + 「既定切り替え日をまたぐ run は途中再開せず最初から回し直す」運用ルール。

## 実装案

1. `solverConfig.cpp`: node で `scalarGradient` 省略 → `lsq` (由来 `default`)。cell は従来どおり gg。省略時は 1 行警告 `[config] mesh.scalarGradient 省略 → 既定 lsq (2026-09-27 から。以前に gg 既定で始めた run を途中から再開しないこと。旧挙動は gg を明記)`。明示 gg/lsq は従来どおり。
2. `runner_sern.py` `FLAG_POLICY` を `2026-09-27` に (前日までの SERN 評価行を学習から外す。chi 既定化と同じ仕組み)。runner は `scalarGradient` を書かない (既定に従う)。
3. 文書: `procedures/solver-settings.md`・`recommended-settings.md` §9 (「2026-09-27 以前の node run は gg。再現は gg 明記」)、`methods/gradient.md`・`discretization.md` §7.3。
4. 検証: (a) node 省略で起動エコー `lsq (default)`、cell 省略で `gg`、明示 gg で `gg (explicit)`。(b) 既定 lsq の 1 step が明示 lsq と同じ配列 (S1 の規則で比較、周期ありの tgv と周期なしの case/48)。(c) 既存の回帰試験 (stage_manifest・ゲート) PASS。

## 懸念

- 同じワークツリー・ブランチで並行している SERN 設計セッションの評価が、既定変更で黙って変わる。`FLAG_POLICY` で DB は分けるが、走行中のキャンペーンは切り替え日をまたぐ。
