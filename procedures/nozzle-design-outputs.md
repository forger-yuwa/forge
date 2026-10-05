# ノズル設計の標準出力と PowerPoint 報告

ノズル設計 (軸対称 axis-Mach チェーン; `design/forge_design`) の結果は、**毎回同じ図・同じ量・同じ並び**で出し、
PowerPoint にまとめる (2026-10-05 ユーザ指示でルール化)。設計ごとに図の種類や並びが変わると比較できないため。

- 実体: `design/forge_design/report/nozzle_report.py` (図・評価量・条件表 → `<run>/report/report.json` と `fig_*.png`) と
  `design/forge_design/report/build_pptx.py` (リポジトリの `.venv-pptx` の python-pptx で `<run>/report/<run>_report.pptx`)。
- **自動実行**: 設計チェーンの NS pass (`feedback/deltastar_loop.py` の `run_pass` / `run_pass0_integral`) の最後に走る。
  失敗しても chain は止めない (ログに「design report: 失敗」と出る)。止めたいときだけ環境変数 `FORGE_NO_DESIGN_REPORT=1`。
- **手で回す** (既存 run・凝縮 ON run など):

  ```bash
  cd design
  .venv-opt/bin/python -m forge_design.report.nozzle_report ../case/<case>/<run> --euler ../case/<case>/<Euler 参照 run>
  ```

  `--no-pptx` で図と json だけ。`.venv-pptx` が無い環境では自動で図と json だけになる。

## 評価の観点 (ユーザ決定 2026-10-05)

軸中心マッハ数は形状の生成器にすぎない。良否は次で見る:

1. **壁がなめらか** — 2 階微分の高周波 (0.5 r_t 移動平均からの残差) を、壁全体・境界層補正 (r_w − r_inv)・設計壁に分けて出す。
2. **試験部に圧力波が無い** — r/r_w = 0.1 の M/M_d − 1 から 10 r_t の P-spline を引いた残差の最大 (軸ノード r=0 は既知の値の癖があるので併記のみ)。
3. **オーバーシュートが小さい** — 膨張の終わり以降 (x ≥ x_E − 15) の M/M_d − 1 の最大。過膨張は凝縮を強める。
4. 出口コアの M (r/r_w 0.05〜0.7 の平均と最小〜最大)、流量比 ṁ_NS/ṁ_Euler、試験部の傾き (1 次近似の両端差)。
5. 凝縮 ON のとき: 凝縮の始まり (軸で g > 1e-4 になる x)、出口の g (軸・コア平均)、最大過飽和度。

試験部は x ∈ [x_E + 2, x_F − 1] (x_E = 軸 M 則の終点、`prepare_info.json`)。

## スライドの構成 (この順)

| # | スライド | 中身 |
| --- | --- | --- |
| 1 | 表紙 | case / run、出口コア M、圧力波、オーバーシュート、出口半径 (凝縮 ON は凝縮の始まり) |
| 2 | 解析領域と境界 | `nozzle.h5` の格子 (間引き表示) と境界パッチを physID ごとに色分け、スロート近傍の拡大 |
| 3 | 境界条件 | `bcondConfig.yaml` から転記 (physID・境界名・kind・値) |
| 4 | 解析設定と物性 | 離散化・対流スキーム・時間積分・step 数・乱流モデル・格子・初期場・壁の δ_r の出所・凝縮設定 / 気体・粘性・熱伝導・r_t・設計 M |
| 5 | 判定ゲート | `check_convergence` (OVERALL と不合格の列)・`check_mesh_quality`・NaN/Inf・準定常 (末尾 5 枚の変動)・壁解像 y₁⁺ — **不合格・未収束もそのまま載せる** |
| 6 | 評価量 | 上の 1〜5 の表と、試験部の M/M_d − 1 の図 (r=0 と r/r_w=0.1、Euler 参照も重ねる) |
| 7〜 | コンタ | マッハ数・M/M_d − 1 (±0.5 %) / (r_t/p)∂p/∂x・数値シュリーレン log10(\|∇ρ\| r_t/ρ) / 静圧・静温 / 凝縮 ON: 液滴の質量分率 g・過冷却度 T_sat − T・過飽和度 S |
| | 軸に沿った分布 | r=0 と r/r_w=0.1 で M・静圧・静温・密度・動圧 0.5ρv²・0.5v²・全圧 (+ 凝縮量) |
| | 出口断面の分布 | 同じ量を r/r_w に対して + 流れ角 atan(v/u) |
| | 壁面の分布 | 壁圧・壁温・C_f・y₁⁺ |
| | 壁の形と微分 | r・r′・r″ (設計壁と重ねる) と設計壁との差 |

## 図の規約

- カラーマップは全量 turbo (skill `forge-contour`)。M/M_d − 1 は ±0.5 %、圧力勾配は ±0.02 で頭打ち、他は 0.2〜99.8 % で切る (切った値は題に書く)。
- 静圧・密度・動圧の軸方向分布は対数目盛。
- 全圧・全温は `solver_density_cuda/tools/total_quantities.py` の `total_state` で h0 から作る (AGENTS.md「出力と後処理の原則」)。h0 の無い古い res では全圧を省く。
- y₁⁺ の判定値 (ゲート表・`report.json` の `metrics.wall_resolution`) は**正式ツール `solver_density_cuda/tools/check_wall_resolution.py` を全 no-slip 壁で実行した結果** (VERDICT・目標超過の面積割合・最大とその位置 x/r_t) をそのまま載せる (2026-10-05, codex result M4: 速度差の近似を x>0 の節点割合で出していて最大 8.7 / 18 % と正式値 13.4 / 29 % より甘かった)。許容面積は `--wall-over-frac` (既定はツール既定 2 %)。壁面分布の図の y₁⁺ は近似で図示だけ。**ソルバの `ypls` は使わない** (AGENTS.md「壁解像確認」)。
- 各スライドのノートに図のファイル名と出典の run / res を書く。

## 報告の扱い

- pptx は run ディレクトリの `report/` に置く (run の成果物なので commit しない)。結論・判定は対応する plan と case README の run 一覧に書く。
- 設計の比較では、比べる run 同士で同じ Euler 参照・同じバイナリを使う (plan `verification-m6-axis-wave-mesh-su2` §4.1 の教訓)。
