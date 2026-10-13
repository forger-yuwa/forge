# 諮問ブリーフ: Euler 参照と入力 δ_E の受入れ、plan レビュー Major 5 件の採否 (plan tooling-nozzle-deltastar-band-root-zone)

AGENTS.md 条件 3 (事前に書いた判定が FAIL・判定不能) と条件 5 (codex の Major の採否)。判断役のスイッチが codex。判断単位でまとめて諮る。

## 1. plan 段のレビュー (`notes/reviews/2026-10-13-tooling-nozzle-deltastar-band-root-zone-plan.md`、GO-with-changes C0/M5/m2)

当方の採否案: **全件採用**。M1 hard/held 点が壁更新に入る → 使用域 [0, x_e] に hard/held があれば NS 後の寸法を止める + 伝播の試験。M2 符号分岐の前に参照流束の有限・正と探索域を検査
(合成の負例 δ_r = −4.98e13)。M3 #18 は両側 |Δ| ≤ ε_f32 (run_0009 は出口を 2.7e-6 r_t 超える)。M4 入力 δ_E の準定常 (下の 3.)。
M5 出口 M の準定常は専用の許容 (末尾の変化と推定の残りの変化 ≤ 2.5e-4、単調なら漸近値も合否範囲内)。m6 V5 は「各設計の受入れた定常解どうし」。m7 methods を実装前に。

## 2. Euler 参照の判別 A/B の結果 (事前登録: plan §5.1 #5)

- A = `case/70.m5_air_wt_demo/run_0012_euler_c2b` (r_t 0.024973561 m、本段 cfl 1.0・relax なし [run_0003 は cfl 6・relax 0.7、段取りの逸脱]、36000 step)。
  `check_convergence --segment`: NOT CONVERGED (全列 2.0〜2.3 桁で下降中)。`qs_series.py`: 出口コア M・無次元の流量 STEADY (既定の許容)。
- B = A の res_36000 から `restart_field` (ビット一致) で同じ設定 12000 step (`run_0012b_euler_c2b_cont12k`)。`check_convergence`: NOT CONVERGED (**stalled/plateau**、全列 0.1 桁で横ばい = 床、ρ 9.6e-8)。
  出口コア M・流量 STEADY。
- run_0010 を A と B に対して修正後の抽出 (`ab_root_zone.py` の zone の腕 = plan §4.1 と同じ差し替え) に通した (`_ab_euler_steps/`):

| 量 | A | B | 差 | 事前の閾値 |
| --- | --- | --- | --- | --- |
| 出口コア M (Euler の末尾 5 枚の平均) | 4.99936 | 4.99937 | ≈1e-5 | 2.5e-4 |
| ṁ_NS/ṁ_Euler | 1.0077795 | 1.0077433 | 3.6e-5 (相対) | 1e-4 |
| max\|Δδ_E\| 全分布 | | | **7.2e-3 r_t (x ≈ −8、縮流部)** | 1e-4 r_t → **FAIL** |
| 同 区間別 [−20,−12) / [−12,−6) / [−6,−2) / [−2,0) / [0,3) / [3,10) / [10,20) / [20,35) / [35,48] | | | 1.7e-3 / 7.2e-3 / 4.9e-3 / 9.4e-5 / 7.6e-5 / 1.5e-4 / 1.5e-4 / 8e-6 / 5e-6 | |
| δ_E x = 0 | −0.003086 | −0.003052 | 3.4e-5 | |
| 流量からの実効値 | −0.002824 | −0.002806 | | |

- 縮流部の生の δ_r は A と B で大きく違う (x −12〜−4 の 5 点: A −0.0068/−0.0076/−0.0032/+0.0002/+0.0003、B −0.0078/−0.0133/−0.0109/−0.0055/−0.0025)。
- **物理壁は x < 0 で δ_r を使わない**: `wall_axismach.py:1721` の式 "x < −L_U: r_U (直管) / −L_U ≤ x < 0: Q(x) (upstream_poly) / 0 ≤ x ≤ x_e: S(x) + δ_r(x)"
  (Q は x = 0 の値と導関数で S + δ_r に接続)。寸法 (`solve_rt_throat`) は物理壁の最小 (x ≈ 0) を使う。
- 判定 (plan の文言): A・B とも収束の受入れ (check_convergence PASS) を満たさず、δ_E の閾値は全分布では FAIL → 「B も不合格なら判別未了」。

## 3. 入力 δ_E の準定常 (レビュー M4 の確認、`case/70.m5_air_wt_demo/deltaE_series.py`・`_deltaE_series_run_0010/`)

run_0010 の 40000〜60000 step の 5 枚を、Euler A を固定して同じ修正済みの抽出・平滑化に通した。

| 区間 | 40k−60k | 45k−60k | 50k−60k | 55k−60k |
| --- | --- | --- | --- | --- |
| [0,2) | 8.5e-5 | 2.0e-5 | 2.8e-5 | 2.0e-5 |
| [2,6) | 4.6e-4 | 4.4e-4 | 4.9e-4 | 5.6e-5 |
| [6,12) | 4.6e-4 | 4.4e-4 | 4.8e-4 | 7.6e-5 |
| [12,25) | 2.8e-4 | 1.0e-4 | 6.6e-5 | 7.6e-5 |
| [25,48] | 2.2e-4 | 1.4e-4 | 9.0e-5 | 5.0e-5 |

- δ_E x = 5: 0.007022 / 0.007054 / 0.006991 / **0.007502** / 0.007446 (50k → 55k で +5.1e-4 の段差、その後は 5.6e-5)。x = 0: −0.003136 / −0.003104 / −0.003083 / −0.003100 / −0.003086。
- `check_quasisteady --series-csv` (末尾 5 枚): x = 0・1・10・20・30・48 は STEADY、x = 3 は DRIFTING (0 交差付近で相対量が大きい、絶対値 8.5e-4)。
- run_0010 は `check_convergence` PASS、出口コア M・流量 STEADY (4.990617〜4.990682)。成果物の壁の総熱量は TRANSIENT-UNSETTLED。

## 当方の仮説 (未検証)

- 縮流部の δ_E の差は、低マッハの縮流部で帯局所の比 q_NS/q_E が Euler の床の雑音に敏感なためで、壁 (x ≥ 0) と寸法には入らない。閾値は使用域 [0, x_e] に限るべきだった (事前の文言の不備)。
- δ_E の 50k → 55k の段差は、edge の帯の位置 (梯子から決まる離散的な選択) の切り替わりで、流れの過渡ではない可能性。

## 問い

**この状態から最終 NS に進むために、(1) Euler 参照を何を根拠に受け入れる (または受け入れない) か、(2) 入力 δ_E を何で決めるか (最終 1 枚・末尾の平均・その他) を、
それぞれ 1 つに絞って推奨し、判定条件を事前登録の形で示してほしい。** 加えて、1. の採否案に反対の点があれば挙げてほしい。
候補 (縛られなくてよい): (a) 閾値を使用域 [0, x_e] に限って A/B を再判定 (事後の変更と明記)、(b) Euler を cfl 6・relax 0.7 (run_0003 の段取り) で回し直す、
(c) δ_E は末尾 2 枚 (55k・60k) の平均、(d) 段差の原因 (edge の帯の切り替わり) を先に確かめる。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.cgns`, `*.stl`, `plans/README.md`, `delta_r_equiv_diag.json` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,220p' plans/active/tooling-nozzle-deltastar-band-root-zone.md`
- `sed -n '15,140p' notes/reviews/2026-10-13-tooling-nozzle-deltastar-band-root-zone-plan.md`
- `sed -n '1715,1725p;1785,1800p' design/forge_design/geometry/wall_axismach.py`
- `sed -n '310,375p' design/forge_design/metrics/deltastar.py`
- `cat case/70.m5_air_wt_demo/_deltaE_series_run_0010/deltaE_series.csv`
- `cat case/70.m5_air_wt_demo/_ab_euler_steps/A/zone_root_ab.json case/70.m5_air_wt_demo/_ab_euler_steps/B/zone_root_ab.json`
- `grep -n "run_0012" case/70.m5_air_wt_demo/README.md`
