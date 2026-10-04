# 諮問: SERN 3D カウル後縁の低温スポット (plan tooling-nozzle-sern-3d §5.1 R5h) — 今の形状で掘るか、有限厚後縁 plan に任せるか

## 問い (1 つ)
R5h を「厚さ 0 の後縁形状のまま診断を続ける (§5.1 R5h の (a) 時系列と EOS 前後、(b) 回復温度 vs 対流流束)」か、
「形状起因として [tooling-sern-mesh-blocking.md](../../../plans/active/tooling-sern-mesh-blocking.md) の有限厚後縁 `cowl_base` に委ね、
R5h は監視項目 (力係数への影響の上限を測るだけ) に下げる」か。後者なら、影響の上限を測る最小の A/B を 1 つ指定してほしい。

## 観測事実 (2026-09-27、AWS、バイナリ 2fa3826c、chi 1・lsq・出口 outflow、生産条件 frozen_tp・等温壁 1000 K・SST 低 Re)
| run | 格子 | T 最小 | 位置 (x,y,z) m | wall_dist | T<120 K | T<150 K | GATES |
| --- | --- | --- | --- | --- | --- | --- | --- |
| case/46.sern_design/run_0973_3d_g1_chidef | g1 (first_wall_frac 2.56e-3) | 176.6 (プルーム遠方、後縁ではない) | (1.007,0.125,0.25) | 2.1e-1 | 0 | 0 | PASS |
| run_0971_3d_g3_chidef_cont16k | g3 (1.6e-4) | **105.4** | (0.12000,−0.01053,0.09719) | 3.5e-5 | 1 | 1 | PASS |
| run_0972_3d_g4_chidef_cont40k | g4 (4e-5 + x ブレンド) | **94.7** | (0.12000,−0.01053,0.09719) | 2.9e-5 | 1 | 3 | PASS |
- 全場の T p0.01 は 196–216 K。力係数 4 量は全 run STEADY、g3→g4 は §8 許容内 (plan §4.46)。床張り付き 0。
- 位置は x = L_cowl (0.12 m) の後縁。後縁 station では cowl_in(5)/cowl_out(6) が同一ノード ID を共有 (25 z station 全部、y 差 0 = 厚さ 0 の閉じた後縁)。
  低温節点は後縁ノード (y −0.010499) の直下 (y −0.01053、外部流側) の内部節点。z = 0.09719 は W/2 = 0.1 の 6 z-station 内側 (側壁は L_sw = 0.08 m で終わっており x = L_cowl には側壁が無い)。なぜこの z なのかは未確認。
- 板厚則は `np.interp(xs, [−L_up, 0.8·L_cowl, L_cowl], [t, t, 0])` (cowl_thickness 0.005 H) で後縁で厳密に 0 (plan §4.35)。
- 旧データ: 断熱 CPG run_0415 で 37.6 K (8 節点)、等温壁 run_0418 で 109.6 K (1 節点) (plan §4.35/§4.37.1)。

## 期待値と出典
- 物理的には後縁の外部流側で ~200 K を大きく下回る理由は無い (外部流の静温 ~200 K、プルーム過膨張の最小 176 K)。出典: 本表の g1・場の p0.01。

## 実施済み
- 等温壁化で 37.6 → 109.6 K に改善 (回復温度は寄与したが主因でない、plan §4.37.1)。
- 本日の現行コード (chi・lsq・メッシャ修正・outflow) でも残存し、格子細分で悪化 (上表)。

## 仮説 (未検証)
- H1: 厚さ 0 の後縁 1 station の形状起因 (2772 K 排気と外部流が 1 本の壁ノードを挟む)。有限厚ベースで消える。
- H2: 数値 (対流流束・リミッタ) の局所不整合で、形状を変えても残る。

## 読んでよいファイル
- plans/active/tooling-nozzle-sern-3d.md の §4.35・§4.37.1・§4.46・§5.1 (R5h, R5j)
- plans/active/tooling-sern-mesh-blocking.md の冒頭〜§4 の cowl_base 関連と §5.1 (B1d, B6)
- design/forge_design/meshing/mesh_sern3d.py (後縁・板厚・側端テーパ部分のみ)
巨大な h5・ログは読まないこと。
