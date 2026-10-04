# 諮問 (設計): SERN 3D R4d — 遠方境界と領域独立性の検証計画 (plan tooling-nozzle-sern-3d §5.1 R4d) を事前に固定したい

## 背景
- R4d (codex 2026-09-19 plan M4 採用): `side_far` の slip 固定 (反射) を見直し、生産 TP・SST で C_L/C_M まで含めた領域独立性を測る。許容値を係数ごとに数値で固定する。
- 旧 R4 (`case/46.sern_design/r4_domain_study.py`、run_0102): smoke 設計・加速点・3D **Euler**・|ΔC_T|<0.002 のみ。成果物は checkout に無く結果も README に残っていない (「実行中」のまま)。
- 2026-09-27 に現行コード (chi 1・lsq・R5r メッシャ修正・出口 outflow) で生産格子列を取り直した: g1 `run_0973` / g3 `run_0971` / g4 `run_0972`、すべて GATES PASS・4 量 STEADY、g3→g4 は §8 許容内 (§4.46)。

## 現在の境界と領域 (生産 problem_3d_prod_m6on_*.yaml、H = 0.1 m)
- inlet_nozzle / inlet_ext: 一様流入 (外部流 M ≈ 6、Ux 1788 m/s、Ps 2851 Pa)。outlet・bottom: outflow (外挿)。
- **top_out: slip** (`top_out_kind: slip`、「3D は top_out を slip 壁 (ランプ延長) にする (run_0027 で確立)」)、**side_far: slip**、sym (z=0): slip、underside_far: slip。
- 寸法: W 2.0 (半幅 W/2)、Z_ext 1.5 (側方外部の幅)、top_depth 2.0、bot_depth 3.0、x_out_extra 2.0 (ノズル出口から下流)、単位は H。
- 生産の力係数はノズル面のみ (ramp 幅内 + cowl + ダクト側壁)、機体面は別枠。

## 設計案 (呼び出し側、未確定)
- 基準格子: **g1** (1.2M セル、段階起動 + 本段 20000 step で約 40 分) を使い、1 因子ずつ変えた変種を回す。理由: 領域・遠方 BC の感度は外部流 (非粘性支配) で、g1→g3 の差は摩擦が主 (§4.40/§4.46)。
- 変種 (各 1 因子):
  V1 side_far: slip → outflow、V2 top_out: slip → outflow、V3 Z_ext 1.5 → 2.25 (nz_out も比例)、V4 top_depth 2.0 → 3.0、V5 bot_depth 3.0 → 4.5、V6 x_out_extra 2.0 → 3.0。
- 判定 (案): base に対する |Δ| が §8 許容の 1/4 以下 (C_T・C_T_with_shear・C_L 0.0005、C_M 0.0125) を「領域独立」、両者 check_quasisteady STEADY の上で平均差 + 双方の振幅で評価。
  1/4 にする理由: 格子 (g3→g4) の差と領域の差が同じ許容を食い合わないように。
- 補助診断: 遠方境界上での自由流からの圧力偏差 (最大 |p/p∞ − 1|) と、そこへ波が当たっているか。
- 費用: 6 変種 × ~40 分。3D 変換は OOM 回避のため 1 本ずつ (AWS 16 GB)。

## 問い
1. 基準格子を g1 にしてよいか (g3 が要るか)。
2. 変種の集合と、一度に 1 因子の設計でよいか。slip→outflow の切り替えは「反射の有無」を測れるか (outflow 外挿は超音速流入側があると不適切か)。
3. 許容値 (§8 の 1/4) の妥当性。事前に固定すべき数値を指定してほしい。
4. 省いてよい変種・足すべき変種。

## 読んでよいファイル
- plans/active/tooling-nozzle-sern-3d.md §4.46・§5.1 (R4d)・§8
- design/forge_design/evaluate/runner_sern3d.py (_bcond_config)、design/forge_design/meshing/mesh_sern3d.py (寸法パラメタ部分のみ)
- case/46.sern_design/problem_3d_prod_m6on_g1.yaml、case/46.sern_design/r4_domain_study.py
巨大な h5・ログは読まない。
