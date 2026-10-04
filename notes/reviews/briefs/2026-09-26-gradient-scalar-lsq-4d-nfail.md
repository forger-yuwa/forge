# 諮問ブリーフ: #4d 再判定で N 区分の 1 配列が 2 倍規則 FAIL — 解釈と次の一手 (2026-09-26)

plan: `plans/active/gradient-scalar-lsq-unification.md` (§5.1 #4d、§6「ダンプ非干渉と S1 の 1 step 比較の区分」)。
前回諮問: `notes/reviews/2026-09-26-gradient-scalar-lsq-4a-closeout-diagnose.md` (E/P/N 区分を提案、全件採用済み)。

## 決めてほしいこと

1. 下の N 区分 FAIL (tgv lsq res_1 `res_ro`) をどう扱うか: #4d の VERDICT は事前に書いた条件 (E 厳密 ∧ P 集合内) で A だが、§6 の N 区分には「S1 の訂正規則 (2 倍規則) を適用」と書いてある。
2. 親の案 (下の H) で良いか。良ければ判定規則 (反復数・閾値) をここで固定してから回す。

## 読んでよいファイル

- `case/09.Taylor-Green/_g0_lsq_seam/DUMP_RECHECK_tgv.txt`・`DUMP_RECHECK_case48.txt` (今回の結果、全文)
- `case/09.Taylor-Green/_g0_lsq_seam/dump_recheck.py` (判定スクリプト)
- plan §5.1 #4d・§6 の該当 bullet
- `solver_density_cuda/cuda_forge/periodicNode_d.cu:50-80` (`periodicGather1ToRoot_d`)

## 観測事実 (バイナリ sha256 `a398afdc…`、1 step、既存 run の再解析のみ)

- tgv (三重周期 SST + ξ): ダンプ無し `~/sglsq/s1h/s1_tgv_{gg,lsq}_{a,b,c}`、ダンプ有り `~/sglsq/pg/pg_tgv_{gg,lsq}_{a,b,c}` の 12 本。
  - E0: pre-gather NS 18 配列 (init・loop1) がダンプ有り 6 本でビット同一。
  - E: res_0 の状態量 35 配列、NS 18 勾配の group 外節点 (res_0・res_1) が 12 本すべて基準とビット一致。
  - P: 2 member 以上の節点で基準と値が違う (run, 節点) は res_0 dUxdx 368・dUxdz 352、res_1 dUxdx 792・dUxdy 64・dUxdz 92。**集合外は 0**。
  - N: 88 配列が全対ビット一致。残りは 2 倍規則 PASS が 43 行、**FAIL が 1 行**: res_1 `res_ro` (lsq) — 無し同士 n≤11 / max 2.98e-8、有り同士 n≤14 / 2.98e-8、無し–有り **n≤40** / 2.98e-8。最大差は 3 組とも同じ (1 値)。同じ配列の gg 側は無し同士 12、有り同士 **45**、無し–有り 42 で PASS。
- case48 (周期無し SST 平板): E・P・N すべて通過 (N は 17 配列すべて 2 倍以内)。VERDICT A。

## 期待値と出典

- #4d の A/B (事前、plan §5.1 #4d): A = E 厳密 ∧ P 集合内。N は「不一致数と最大差を 3 組とも出す」。
- §6 N 区分: 「S1 の訂正規則 (両側 3 本のノイズ対、最大差 ≤ 2 倍かつ不一致数 ≤ 2 倍) を適用」。

## 仮説

- H: `res_ro` は面ループ atomicAdd の残差で、lsq/gg はこの配列の計算経路に影響しない (スカラー勾配は質量残差に入らない; 1 step 目の質量流束は NS 再構成のみ)。不一致数は 3 本同士の最大で見ているので標本が小さく、gg 側の有り同士 45 が示すように同じ配列のノイズの不一致数は 11〜45 の幅で出る。**3 本ずつでは不一致数の 2 倍規則の検出力が足りない** (前提 plan #6a で廃止した「桁規則」と同じ種類の脆さ)。
- 親の案: (a) 本件は gg/lsq でノイズ源が同じ配列なので、ノイズ対を kind を問わずプールする (無し同士 6 + 有り同士 6 = gg/lsq それぞれ 3 対ずつ) と n_noise = 45 で PASS — ただしこれは結果を見てからの規則変更になるので採らない方がよいと考える。(b) 代わりに、判定を変えずに**反復を増やす**: lsq のダンプ無し・有りを各 +5 本 (計 8 本ずつ、1 step で数分) 回し、同じ 2 倍規則を 28 対 + 28 対のノイズ対 vs 64 対の交差で判定する。この規則 (反復数・閾値) を回す前に固定する。FAIL が残れば B′ 相当 (ダンプが残差のノイズ分布を変えている) として #4a を「診断専用・非干渉は未確認」と記録し、S2 には進む (ダンプは S2 で使わない)。
- 問い: (b) の S2 への扱い (ダンプは S2 で使わないので S2 の前提にしなくてよいか) も含めて判断してほしい。
