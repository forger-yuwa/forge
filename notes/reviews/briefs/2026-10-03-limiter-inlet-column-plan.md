# 諮問ブリーフ: 入口列のリミッタの揺れ — 調査計画 (§4・§6) の点検

エスカレーション条件 1 (plan §4・§6 を新規に書く)。対象 plan: `plans/active/limiter-inlet-column-oscillation.md`。

## 観測事実
plan §3.1 にすべて転記した (run パス・設定・数値)。要点: case/16 (平面 2D node SST 陰解法) の残差の床が
`limiterScaled 1` で `limiterScaled 0` の約 20 倍 (凝縮の有無と無関係)、Σres² の 71 % が入口の最初の 2 列の節点、
入口境界の壁の隣 3 節点で limiter_Ux = 0、ε̂ が丸め以下という仮説は否定 (勾配増分が ε̂ の 300 倍以上)。

## 期待値と出典
`limiterScaled 0` では同じ IC から 3.3 dec 下がって ALL PASS (run_0525)。limiter-config-simplify §4.3 で ψ 凍結は入れないと決定済み。

## 実施済み
condensation-two-phase-transport §5.1 #4i・#4l (2×2)。後処理 `notes/investigations/2026-10-03-twophase-onoff/{floor_analysis,eps_probe}.py`。

## 仮説
plan §4.2 の H1〜H4。

## 問い
1. §4 の調査順序と A/B の選び方は妥当か。抜けている有力仮説 (例: 入口 BC `inlet_Pressure` の node 実装、半 CV の再構成距離、
   SLAU の低マッハ域での ψ 感度、`nodeWallDirichlet` で固定された壁ノードを近傍に含む極値) はあるか。コード (`solver_density_cuda/cuda_forge/limiter_d.cu`,
   `limiterFunctions_d.cuh`, 入口 BC の node 実装) を読んで、入口境界節点の ψ がどの近傍集合・どの評価点で決まるかを確認してほしい。
2. §4.3 の判定基準 (入口 2 列の占有 <10 % かつ床 ≤1/5 で主因、占有 ≥50 % または変化 2 倍以内で否定) は事前登録として妥当か。
3. H1 を config だけで切り分ける方法はあるか (既存の診断オプション・境界の 1 次化キー。ただし `mesh.bndFirstOrder` は使用禁止)。
4. A/B の長さ (24000 step) と IC (L1 の収束床の場) で判定できるか。L0 が L1 の床から 1.5 dec 下がるのに 48000 step かかった点を考慮して。

読んでよいファイル: 対象 plan、`plans/active/limiter-config-simplify.md`、`methods/limiter.md`、`methods/boundary.md`、
`solver_density_cuda/cuda_forge/` 配下、`case/16.nozzle_wys/README.md`、上記後処理スクリプト。
