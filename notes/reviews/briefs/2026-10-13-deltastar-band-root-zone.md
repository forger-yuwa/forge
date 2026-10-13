# 諮問ブリーフ: 帯局所抽出が δ_r = r_w を返す真因と、M5 デモ (冷却壁) の出口 M を閉じる経路 (plan tooling-nozzle-deltastar-band-root-zone)

AGENTS.md 条件 1 (新規 plan の §4・§6) と条件 4 (「原因は○○」と plan に書く前)。判断役のスイッチが codex。
ユーザ決定 (2026-10-13): 抽出の頑健化と #18 を別 plan で起こし、codex に諮ってから実装 → NS の δ_E 全分布で壁と r_t を作り直して出口 M を 4.9975〜5.0025 に。

## 観測事実

- case/70 (空気 M5、Pt 3 MPa、Tt 1250 K、等温壁 350 K、断熱壁温 ~1150 K) の最終 run `run_0010_ns_c2b_final`: `check_convergence` PASS、
  `check_quasisteady --quantity machmax,pmax` ALL STEADY、壁解像 PASS、**出口コア M 4.9907 (−0.19 %)**、ṁ_NS/ṁ_Euler 1.0077。
- 帯局所抽出 (`extract_and_merge(.., band_select="edge")`) の生の δ_r が x ≈ 1.09〜2.83 r_t の 99 断面 (run_0009) / 38 断面 (run_0010) で
  **壁の半径 r_w と一致** (1.15〜1.49 r_t)。理由は soft の `core_frac_sens` だけで hard_ok = 1 なので平滑化に入り、δ_E が x = 0〜3 で +0.05〜+0.77 r_t に広がる。
- 機構 (当方の読み): `band_local_deficit` (`design/forge_design/metrics/deltastar.py:255-283`) の q_ref = (a + b(y − y_b)) q_E は帯の線形フィットを全半径へ外挿する。
  帯の傾きが大きい (run_0009 x=1.09: 帯の幅での比の変化 −5.1 %、`band_slope`) とコアで q_ref < 0、累積流量 F (壁から) が途中で減り末尾が D を下回る。
  `np.interp(D, F_from_wall, r_desc)` (282 行) は x ≥ xp[-1] で fp[-1] (軸 ≈ 0) を返し δ_r = r_w。直前の検査 `D > F_zone` (帯の下端での F) は通る。
  edge の帯は固定 (`max_retry = 0`) なので帯の傾きでの再試行 (240 行付近) も効かない。
- **A/B** (`case/70.m5_air_wt_demo/ab_root_zone.py`、forge 0 step、Euler 参照 run_0003): 探索を壁〜y_b の区間の F に限る腕 zone。

| NS run | 腕 | δ_r = r_w の断面 | max \|raw\| (0≤x≤3) | δ_E x=0 / 1 / 2 / 3 / 10 / 48 [r_t] | 流量からの実効値 δ_eff,ṁ |
| --- | --- | --- | --- | --- | --- |
| run_0009 | orig | 99 | 1.48 | +0.2403 / +0.5821 / +0.7737 / +0.6993 / 0.0411 / 0.2857 | −0.003421 |
| run_0009 | zone | 0 | 0.0034 | −0.003601 / −0.002972 / −0.001804 / +0.000092 / 0.0308 / 0.2857 | 同 |
| run_0010 | orig | 38 | 1.49 | +0.0525 / +0.1959 / +0.3002 / +0.2995 / 0.0352 / 0.2867 | −0.002978 |
| run_0010 | zone | 0 | 0.0029 | −0.003022 / −0.002336 / −0.001150 / +0.000832 / 0.0320 / 0.2867 | 同 |

  zone の腕: 区間の中で F が単調でない呼び出し 0 / 114,323、全半径の F の末尾が D を下回る呼び出し 5,521。δ_eff,ṁ = r_throat_phys − √(ṁ_NS/ṁ_E) (`deltastar.py:93-117`)。
- 手元の抽出の記録 33 個の走査: 同じ症状は case/70 の 3 個と case/45 `run_0089_ns_c2pin_pass1` の縮流部 8 断面 (x −5.08〜−4.45、断熱壁)。
- #18: 設計の出口 x_e 48.031377427 r_t、抽出の表の終点 48.031375885 (短い分 1.54e-6 r_t = 3.9e-8 m)。格子の座標は float32 [m] (x_e·S の 1 ulp 1.19e-7 m = 4.8e-6 r_t)、
  抽出は x を r_t 単位の float32 で持つ (1 ulp 3.8e-6 r_t)。`SingleBSplinePhysicalWall.applicability` (`wall_axismach.py:1452-1455`) が [0, x_e] を覆うことを要求して止める。

## 実施済みの操作と結果

- 尺度の A/B (H1: Euler と NS の r_t の違いが原因) は棄却済み (`_ab_extract_scale/`、D 0.780 → 0.790 r_t)。
- C2 方式 (出口の δ だけ k_f で合わせる) は 2 回回した (run_0009 4.9918、run_0010 4.9907)。スロートは CONTUR の +0.001 r_t のまま。

## 当方の仮説と計画 (plan §4・§6 を読むこと)

- 真因は §3 の `np.interp` の全半径の探索 (冷却壁に固有ではない)。修正は探索を積分域に限り、区間で F が単調でなければ hard。
- #18 は「同じ出口であることを照合したうえで終点を x_e に揃える (ε_f32 = 2 回の丸め)」。
- 出口 M は、修正後の抽出 (Euler 参照を最終の r_t で回し直し) で δ_E 全分布 → `solve_rt_throat --prev` → 最終 NS で閉じる見込み。

## 問い

**(1) 真因の記述 (§3) に反証や見落としがあれば 1 つ挙げ、(2) plan §4・§6 のうち、実装後に「出口 M が 4.9975〜5.0025 に入らない / 入っても理由を取り違える」
形で失敗する可能性が最も高い点を 1 つ挙げて、防ぐ変更を 1 つ推奨してほしい。** 候補 (縛られなくてよい): (a) soft の `core_frac_sens` を平滑化に入れ続けること、
(b) edge の梯子が正の δ だけを使う (`edge_band_positions` の `d > 0`、`deltastar.py:327`) ため冷却壁のスロートで縁が NaN になること、
(c) §4.3 の Re^−0.2 の換算と、r_t 単位の設計を据え置いて実寸の長さが動くこと、(d) V3 (c) の許容が結果を見た後に書かれていること、
(e) V5 の対照 (run_0010) と比べて変えたものが 1 つでないこと (壁・r_t・Euler 参照)。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.cgns`, `*.stl`, `plans/README.md`, `delta_r_equiv_diag.json` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,200p' plans/active/tooling-nozzle-deltastar-band-root-zone.md`
- `sed -n '93,131p;185,290p;310,330p;560,595p' design/forge_design/metrics/deltastar.py`
- `sed -n '63,125p;288,357p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '1126,1150p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '1440,1460p' design/forge_design/geometry/wall_axismach.py`
- `sed -n '1,60p' case/70.m5_air_wt_demo/ab_root_zone.py`
- `sed -n '15,25p' case/70.m5_air_wt_demo/README.md`
- `cat case/70.m5_air_wt_demo/_ab_root_zone/run_0010/zone_root_ab.json`
