# 諮問: Euler の全温の異常 — E1 の結果と、格子による違いの観測を受けた次の手

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 4・7。ユーザ決定「切り分けを続ける」(2026-10-07)。
plan: `plans/active/verification-case45-euler-total-enthalpy.md` (全文)。関連 `plans/active/discretization-moc-axis-limit-and-corrector.md` (§9 の V5〜V5c)。作業ツリー `/home/sano/work/forge-integ-1005` (HEAD c67fc07d)。

## 観測事実

- **E1** (`case/45.isobutane_m6_d155/euler_t0_stage_ab.py`、出力 `_band_ab/euler_t0_stage_ab.json`、AWS で実行): 等エントロピー IC の prep (起動前) は全温 1599.99〜1600.00 K。soft 段 (1 次・cfl 0.5・3000 step) の後は 1245〜1915 K (max\|T₀ − Tt\| 354.8 K、軸上 x = 52.4 r_t)、+1 K を超える節点 95,983・−1 K を下回る節点 91,475・+100 K を超える節点 707。本段の res_0 の保存量は ro・roU* が soft 段の後とビット一致、roe・roY* は一致しない。
- **格子による違い** (2026-10-07、主セッション、別ツリー `/home/sano/work/forge/case/45.isobutane_m6_d155/` の保存物、`total_quantities.py --Tt 1600` の経路 A):

  | run | 格子 ni × nj | wall_first_frac | cfl | 全温の最大 − Tt | Tt + 1 K を超える節点 |
  | --- | --- | --- | --- | --- | --- |
  | run_0047・0059〜0064・0086 (Euler、壁 interp/fit/pincal) | 1100 × 65 | 0.005 | 6 | +0.07〜0.12 K | 0 |
  | run_0120 (Euler、rerun) | 1100 × 65 | 0.005 | 2 | +0.073 K | 0 |
  | run_0114 (Euler、G1、壁 joint、出口較正に使用) | 2000 × 97 | 1.3e-5 | 2 | +241 K | 3,490 |
  | run_0140〜0145・0150〜0160 (Euler、G1、AWS) | 2000 × 97 | 1.3e-5 | 2 | +334〜373 K | 3,500〜6,400 |

  run_0120 と run_0114 の `solverConfig.yaml` は、`cfl`・`nStepOuter`・`outStepInterval` 以外に差が無い (node、isAxisymmetric 1、axisCentroidShift 1、convMethod 1、limiter 2、SLAU、block-DPLUR)。G1 のメッシュ品質: AR の最大約 4290 (`check_mesh_quality --ar-max 5000` で PASS)、skew 0.44。G1 の壁の第 1 セルは約 0.35 µm (NS 用の壁の細分化をそのまま Euler に使っている)。
- **設計チェーンの入力との関係**: CFD でピン止めする初期線は run_0062 (1100 × 65、全温は正常) から取る。出口較正 (`Md_moc_offset`) は細分格子 G1 の Euler run_0113+0114 で決めた。単調壁の採用 (E′) と MOC の V5 は G1 の Euler の比較。
- 関連の経験 (メモリ、未検証の手がかり): 軸対称の r 重みの閉合 Σ r_f S_f が高 AR の壁 CV で float32 だと崩れた例 (`axisym-rweight-closure-fp32`)、node の双対 CV の重心・体積の float32 の桁落ちを直した例 (`node-yp1-dual-geometry-float32-fix`)。

## 仮説 (未確認)

- H1: 擬似時間の未整定 (残差が下げ止まる) と全温のずれは同じ現象の表れ。
- H2: 壁際の極端に薄い・高 AR の CV (G1 の NS 用の細分化) で、Euler のすべり壁・node の双対 CV・軸対称の扱いのどこかが全エンタルピーを保存しない (幾何の float32 の精度を含む)。
- H3: 格子の解像度 (ni × nj) の違いによる。

## 問い

1. H1〜H3 を判別する、安い順の A/B。候補: (a) G1 と同じ ni × nj で wall_first_frac だけを 0.005 にした Euler (同じ壁・同じ設定、IC は等エントロピーか G0 からの補間)、(b) 1100 × 65 で wall_first_frac だけを 1.3e-5 にした Euler、(c) 既存の場で全温の異常の節点と CV の AR・幾何の閉合誤差の関係を見る (0 step)、(d) 幾何を倍精度にしたビルドでの比較。どれを先に、何を変え、何を固定し、判定をどう書くか。
2. 出口較正 (run_0113+0114 由来の `Md_moc_offset`) と、G1 の Euler の比較に基づく判断 (単調壁の E′、MOC の V5) への影響をどう見積もるか。今の段階で何を「影響あり/なし」と言ってよいか。
3. 次の手の登録として plan §6 に書く条件。

## 読んでよいもの

- 上記 2 つの plan、`case/45.isobutane_m6_d155/{euler_t0_stage_ab.py, README.md}`、`_band_ab/{euler_t0_stage_ab.json, moc_v5c_thermo_ab.json}`
- `design/forge_design/evaluate/runner_axismach.py` (`mesh_params`・`prepare`)、`design/forge_design/meshing/mesh2d.py`
- `solver_density_cuda/tools/{total_quantities.py, check_mesh_quality.py}`、`procedures/solver-settings.md`
- 保存物は別ツリー `/home/sano/work/forge/case/45.isobutane_m6_d155/` (run_0062・run_0114・run_0120 の `solverConfig.yaml`・`prepare_info.json`・`MESH_QUALITY.txt`・`res_6000.h5`)

編集は禁止。
