# 諮問ブリーフ: M5 デモ (冷却壁) の出口マッハ数 −0.19 % を詰める経路 (case/70、plan tooling-nozzle-wt-deliverables §5.1 #15・#18)

AGENTS.md 条件 4 (承認済みの手順に無い修正・再実行の前) と 6 相当 (設計チェーン本体の変更の前)。判断役のスイッチが codex。ユーザ指示 (2026-10-13)
「その前に codex に諮って」。

## 対象と条件

- case/70.m5_air_wt_demo: 空気 (semiperfect)・Pt 3.0 MPa・Tt 1250 K・M_d 5・物理スロート径 50 mm (r_t は物理スロート径から解く)・入口半径 150 mm・
  直管 0.20 m・縮流部 0.30 m・設計スロート→出口 1.20 m・**等温壁 350 K** (断熱壁温 ~1150 K に対し強い冷却)。仮の条件 (実条件ではない)。
- 設計: 要求カード → `wt_card` → axis-Mach (Hall の出発線、`initial_line: hall`、`Md_moc_offset` なし、joint 壁・single_bspline の物理壁)。
- 段取りと結果 (case/70 README の run 一覧):

| run | 内容 | r_t [m] | 結果 |
| --- | --- | --- | --- |
| run_0003_euler_36k | Euler (問題 `problem_m5_air_final_euler.yaml`) | 0.024968985 | NOT CONVERGED (2.6 桁で頭打ち)。出口コア M 4.9993・軸出口 4.9996 (参考) |
| run_0006_ns_pass0_iface_60k | NS pass 0 (CONTUR k_f 1 の物理壁、等エントロピー IC、60000 step) | 0.024968985 (Euler と同じ) | NOT CONVERGED。熱量は 12k 1.51 MW → 60k 0.73 MW と未整定 |
| run_0009_ns_c2_final | C2 1 回目: run_0006 の出口 δ_E に k_f 0.837 を合わせ、カードから作り直し (r_t を物理スロート径から解き直し)、IC = run_0006 (cross-mesh)・60000 step | 0.024973230 | 出口コア M 4.9918 (−0.16 %)、壁解像 FAIL (2.4 %) |
| run_0010_ns_c2b_final | C2 2 回目: run_0009 の δ_E に k_f 0.8246、第 1 層 0.89 倍、IC = run_0009・60000 step | 0.024973561 | `check_convergence` PASS (warm start 直後の過渡から)、`check_quasisteady --quantity machmax,pmax` ALL STEADY、壁解像 PASS。**出口コア M 4.9907 (−0.19 %)**、r/r_w 0.1 の波 0.004 %、ṁ_NS/ṁ_Euler 1.0077 |

## 観測事実 (抽出 `extract_and_merge(.., band_select="edge")` の `delta_r_next.csv` の `delta_E` 列と `delta_r_extract_summary.json` の `massflow`)

| NS run | Euler 参照との r_t の差 | δ_E [r_t] x=0 / 1 / 3 / 10 / 30 / 48 | 平滑化の resid_rel_rms | ṁ_NS/ṁ_E (`massflow_ratio`) | delta_r_throat_eff (= r_throat_phys − √ratio) |
| --- | --- | --- | --- | --- | --- |
| run_0006 | 0 | −0.0009 / +0.0000 / +0.0044 / 0.0356 / 0.1738 / 0.2892 | 26.3 | 1.003407 | −0.00046 (物理壁のスロート 1.001242) |
| run_0009 | +0.017 % | **+0.2403 / +0.5821 / +0.6993** / 0.0411 / 0.1697 / 0.2857 | **1137.7** | 1.009005 | −0.00342 (物理壁のスロート 1.001072) |

- `massflow_ratio` (`design/forge_design/metrics/deltastar.py:93-117`) は実寸の流量の比 (Euler と NS の r_t が違うと面積比 (0.0249732/0.0249690)² = 1.00034 だけ比が上がる)。
- 平滑化の正値制約は断熱壁だけ (`design/forge_design/feedback/deltastar_loop.py:83`)。等温壁の本ケースでは負の δ も表せる。
- run_0006 と run_0009 で NS の流量が 0.56 % 違う (6.592 → 6.629 kg/s、物理スロートの差は 0.017 % 程度)。run_0006 は熱量が未整定だった。
- 出口 M の低下は「スロートを 0.8 % 多く流れる → 実効の面積比が小さい」で桁が合う (A/A* の感度で ΔM/M ≈ −ΔA/A / 4 ≈ −0.2 %)。

## 当方の仮説 (未検証)

H1. run_0009 の抽出がスロート付近で壊れたのは、Euler 参照と NS の寸法 (r_t) が 0.017 % 違うため (run_0006 は同じ r_t で健全)。
H2. 冷却壁でスロートの実効の δ* は負 (−0.0005〜−0.003 r_t)。C2 (出口の δ だけ k_f で合わせる) はスロートを CONTUR (+0.001 r_t) のまま残すので、
    実効スロートが大きすぎて出口 M が下がる。

## 提案している段取り (ユーザに提示済み)

1. Euler を最終形 (run_0010 の問題 `problem_m5_air_c2b_euler.yaml`、r_t 0.024973561) で回し直す。
2. run_0010 の場を新しい Euler 参照で抽出し直し、スロート付近が健全か (流量からの実効値と整合するか) を見る。
3. 健全なら、NS 後の寸法 (`solve_rt_throat --prev --delta-r-out`、δ_E の全分布) で r_t と補正の表を作り、`prepare_ns(delta_r_csv=...)` の壁で最終 NS。
   ただし抽出の表から single_bspline の物理壁を作る段は、表の終点が float32 の格子座標で設計壁の出口より 8e-7 r_t 手前になり
   `SingleBSplinePhysicalWall.applicability` (`design/forge_design/geometry/wall_axismach.py:1450-1460`) で止まる (plan §5.1 #18)。
   直すなら表の端を許容つきで出口まで延ばす (設計チェーンの変更)。
4. 健全でなければ、スロートは流量からの実効値で決める (置く仕組みが無いので変更になる)。
5. あわせて Euler から `initial_line: cfd` と `Md_moc_offset` を入れるか (生産 M6 のレシピ)。

## 問い

**出口コア M を 5.000 ± 0.05 % (M6 の生産の判定基準をデモにも使う) に入れるための経路を 1 つに絞って推奨し、その経路で最初に回す判別 A/B を 1 つ設計してほしい。**
特に: (a) H1 を確かめる最小の試験 (段取り 1・2 で十分か)、(b) スロートの δ* の正本を抽出と流量のどちらにすべきか、(c) #18 の直し方 (表の端の許容つき延長) の安全性、
(d) 段取り 5 を今回入れるべきか (Hall の出発線由来の誤差は Euler で −0.01 % 程度)。候補に縛られなくてよい。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.cgns`, `*.stl`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,40p' case/70.m5_air_wt_demo/README.md`
- `sed -n '60,125p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '280,360p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '93,117p' design/forge_design/metrics/deltastar.py`
- `grep -n "def deltastar_from_core_matched_euler" -A60 design/forge_design/metrics/deltastar.py | head -90`
- `sed -n '1440,1470p' design/forge_design/geometry/wall_axismach.py`
- `sed -n '1290,1330p' design/forge_design/evaluate/runner_axismach.py`
- `grep -n "質量流量実効値\|スロート補正" case/45.isobutane_m6_d155/README.md | head -10`
- `grep -n "initial_line\|Md_moc_offset" procedures/nozzle-design-workflow.md | head -20`
