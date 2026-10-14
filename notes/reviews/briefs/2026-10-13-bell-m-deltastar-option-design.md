# 諮問ブリーフ: ③ベル (M) の δ* 補正を設計ツールの選択肢にする形 (plan §4.4 の新段落・§5.1 #7c)

AGENTS.md のエスカレーション条件 1 (plan §4 の設計方針を新しく書く)。実装の前に点検してほしい。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `*.npz`, `*.png` を読まない** (巨大・バイナリ)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの
- `sed -n '96,97p' plans/active/tooling-nozzle-bell-two-mode.md` (δ* 補正の扱い [ユーザ決定 10-12] と、今回の新段落「δ* 補正の提供の形」)
- `sed -n '149,150p' plans/active/tooling-nozzle-bell-two-mode.md` (§5.1 #7c・#7d)
- `sed -n '153,153p' plans/active/tooling-nozzle-bell-two-mode.md` (§5.1 #8、V3 の前提)
- `sed -n '177,182p' plans/active/tooling-nozzle-bell-two-mode.md` (§6 V3)
- `sed -n '1,43p' case/68.bell_two_mode/rao_m/check_edge_moc_vs_euler.py` と `sed -n '1,17p' case/68.bell_two_mode/rao_m/check_edge_moc_vs_euler.txt`
- `sed -n '1,155p' case/68.bell_two_mode/rao_m/run_corrected_wall.py` (L100 の補正壁を作った台本)
- `sed -n '302,345p' design/forge_design/metrics/bell_deltastar.py` (`contur_bell`: 今は縁を同じ壁の Euler から取る)
- `sed -n '311,400p' design/forge_design/geometry/wall.py` (`DisplacementOffset`)

## 観測事実
- L100 の補正壁 (V-ext、run_0026) は凍結した判定で成立 (壁圧の誤差 2.15 % → 0.74 %、共通の点では 4.03 % → 0.25 %)。η は +3.5e-5 で
  推力はほぼ変わらない。x < 1.17 の既定は V-ext (L100 の暫定、codex 2026-10-13 result diagnose)。
- ユーザ決定 (2026-10-13): 残る −0.24 % は後回しにし、補正を設計ツールの選択肢として組み込み、3 列の比較 (V3) に進む。
- 縁の M の取り方の比較 (`check_edge_moc_vs_euler.txt`): MOC の壁の M (x ≥ 0) + 凍結源 Euler run_0015 の壁 (x < 0) を縁にした contur_v2 と、
  同じ壁の Euler の縁 (較正で使った取り方) の差は、L100 で最大 0.30 % (x 0.41〜2)・出口 −0.04 %、L70 で最大 0.62 % (x 0.41〜1.17)・
  出口 +0.05 %。x = 0 の M は両者 1.18124 で一致 (ピン線の凍結源が run_0015)。較正後の残差は RMS 1.1 %・最大 4.0 %。
- kδ = 1.0289 は L100 の Euler の縁で較正。MOC の縁にすると C が L100 で x ≥ 1.17 の平均 −0.03〜−0.26 % 変わる (係数は較正し直していない)。

## 期待値と出典
- ユーザ決定 10-12: 補正は設計者の選択肢、方法は理論式 + NS 較正 (plan 96 行)。
- codex 10-12 (calibration diagnose): 較正の NS・Euler 基準・検証・補正形状の再生成はキャンペーン準備費として数える。補正の有無の性能差と U_k を分ける。

## 案 (主セッション)
1. 縁を MOC の壁の M + 凍結源 Euler にして、補正壁ごとの Euler を不要にする (kδ は据え置き)。
2. V3 の補正列は L70〜L100 の 4 本を、この生産の手順 (MOC の縁) で作り直して NS で評価する (L100 も、Euler の縁で作った run_0026 を
   使わず作り直す — 手順をそろえるため。差は δ で 0.3 % 以下の見込み)。
3. 補正後の壁で理論式を計算し直す反復はしない。

## 問い (1 つ)
この提供の形 (特に案 1・2) で実装してよいか。最も重い欠陥があれば 1 つ挙げ、どう変えるかを示してほしい。
