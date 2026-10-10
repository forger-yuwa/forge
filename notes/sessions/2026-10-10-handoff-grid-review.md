# 引き継ぎ ②: case/45 の格子の見直し (縮流部の軸方向の密度・境界層と主流のブロック) (2026-10-10)

共通ルールは [`2026-10-10-handoff-common.md`](2026-10-10-handoff-common.md)。発端は `plans/active/time_integration-line-implicit-speed.md` の §5.1 #21。**新しい plan を自分で立てて**進める (例 `plans/active/meshing-nozzle-core-grid.md`)。

## 背景

- 上位の観察 (`notes/reviews/2026-10-10-line-implicit-speed-results-diagnose.md` の「問い 4」と「第 1 仮説」、残差の局在は呼び出し側で再計算済み):
  - 縮流部〜スロートの 2448 列 (全 4719 列の 52 %) が、スロート並みの軸方向の密度 (壁の弧長の間隔 0.0021 r_t) になっている。
  - 半径方向は壁から軸まで 1 本の等比数列 (比 1.073〜1.09)。境界層に 88〜104 節点が入り、主流は 17〜33 節点。軸際の半径方向の幅は 0.14〜0.69 r_t。
  - 縮流部の軸側は dr/dx = 50 (x/r_t −2) で、半径方向に長いセルになっている。
  - 水準に入った時点の残差の床の 85 % が縮流部 (x/r_t −5〜0) の内部にある。
- 速度 plan のふるい (§6.16〜§6.20) では、ライン陰解法の 1 step の時間はラインの長さにほぼ比例した (1 節点で約 0.11 ms)。
- 格子のデータ: `case/45.isobutane_m6_d155/run_0353_m9_L5/res_115000.{h5,xmf}` (手元。ParaView で見られる)。
- 冷却壁の格子の罠 (メモ):
  - `cooled-wall-mesh-precision`: FP64 の変換器・msh 17 桁・近壁は壁法線
  - `mesh-wall-layer-x-blend`
  - `geometry-must-not-follow-mesh-spacing`
  - `mesh-quality-gate`

## やること (案。plan で決める)

1. 縮流部の上流の軸方向の密化を緩める (3〜5 倍粗く) 案と、境界層ブロック (等比、AR 5〜20 まで) と主流ブロック (緩い伸び) に分ける案の、節点数と品質 (`check_mesh_quality.py`) を見積もる。
2. 格子を変えると解が変わる。θ_r・δ_r・Q_w・軸上の M の格子感度を、今の格子と比べて検証する計画を立てる (codex の plan 段を通す)。
3. **生産の格子は、検証が済むまで置き換えない。** 元のセッションは今の格子で float 化を検証している。

- run 番号は case/45 の 06xx を使う。
- 判断の場面 (plan §4・§6 の新規、格子の採否) は AGENTS.md のとおり上位に諮る。
