# case/45 — イソブタン燃焼ガス M6 風洞・出口径 1.55 m の最短全長スタディ

Pt 5.5 MPa / Tt 1600 K / φ=0.9 燃焼ガス (semi-perfect NASA-9) / M_design 6 /
**出口径 1.55 m は δ\* 込み物理壁で定義**。全長 (スロート→物理出口) の最短化と
粘性壁 (δ\* 物理壁 + SST NS) までの検証。計画:
[`plans/active/design-isobutane-m6-d155.md`](../../plans/active/design-isobutane-m6-d155.md)。

- 最短探索: `search_shortest.py` → `search_shortest.json`
  (基準 = shortest-robust study: margin≥1°/topo≥0.02/hard gate/単峰、n1200 確認)。
  **勝者 R2 / L_c 39.3 / M_K 2.7 → x_F = 95.104 r_t** (R3: 39.7/2.8 → 95.433)。
- r_t = (0.775 − δ\*_exit)/(r_F/r_t), 初期推定 0.0771 m (δ\*_exit ≈ 0.66 r_t)。
  A/A\*(M6, Tt1600) = 88.21。
- 凝縮事前見積り: 出口 T 233.5 K vs Tsat(H₂O) 263.9 K = **−30 K 過冷却** →
  dry 本体 + 採用点の凝縮 ON 再評価 (方針 a)。

> **2026-09-04 更新**: 粘性トリム (Md 6.0144, run_0005–0007) は撤回。新チェーン (積分法初期壁 + 固定 Euler 基準・帯局所抽出) で
> **Md 6.0 トリムなし**のまま ṁ_NS/ṁ_E 1.0002・出口面コア M 6.0002 (`run_0038_ns_final_rt77p02` = **最終形**: r_t 77.02 mm、出口半径 0.7750 m、全長 7.325 m; 点列 `points_d155_final_ns.csv`)。旧 v3 のスロート δ\* は実効値の 8 倍過大だった
> ([調査ノート](../../notes/investigations/nozzle-deltastar-throat-review.md), [plan](../../plans/accepted/tooling-nozzle-deltastar-core-matched-euler.md))。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_euler_shortest_dry` | 最短点 (R2/L_c39.3/M_K2.7) の node Euler TP dry (`problem_d155_R2_Lc39.3_tp.yaml`, split_h2o, 段階起動 12000 step) | **dM_max 0.036 % M_d・軸出口 M 6.00003・overshoot 0.034 %・ε_M_rms 0.0057 %・ε_θ 0.0054°・コア 64/65 面**。品質 PASS (AR OK/skew 0.44)・NaN 0・quasisteady ALL STEADY・軸 M 凍結 8k→12k 5.6e-5。残差は warm 床 (rms_ro 1.1e-6, `check_convergence` NOT CONVERGED 判定 = M6 系列の既知パターン)。`metrics.json`/`residual_history.png` | active (**最短点 Euler 検証の正本**) |
| `run_0002_ns_coarse` | NS 中継 (y+~50, `problem_d155_ns_coarse.yaml`, 物理壁 δ\* v1 相関, IC=run_0001) | 完走 12000 step・NaN 0。残差 2–3 桁降下 (ro/roe は低レベルでプラトー = 中継品質として想定内)。出口壁半径 0.7771 m (v1 相関 δ\* ≈ 0.71 r_t) | active (中継) |
| `run_0003_ns_v1` | NS 本計算 v1 (y+~1.4, `problem_d155_ns.yaml`, 相関 δ\*, 24000 step, IC=run_0002) | 完走・NaN 0・品質 SOFT-PASS。**軸出口 M 6.00098・dM_max 0.327 % M_d (ゲート内)**・ε_M_rms 0.075 %。残差 2 桁降下 warm 床。δ\* v3 抽出元 (`dstar_v3.csv`: 相関は中流過大 median 比 0.909, x14 で 0.669) | active (v1 記録・v3 抽出元) |
| `run_0004_ns_v3` | NS 本計算 v3 (CFD 抽出 δ\* 全域採用, IC=run_0003) | 完走・NaN 0・品質 SOFT-PASS・16k→24k 軸 M 凍結 1.2e-4。**δ\* 固定点確認 (2巡目比 median 1.001 [0.981,1.005])**。出口面軸 M **5.9856 (−0.24 % M_d)**、x_E ディップ 5.956 (law 側残差)、一様区間 5.975–5.996 のうねり。metrics の `M_axis_exit` は x_E 評価値であることに注意 | active (**dry NS 固定点・トリム根拠**) |
| `run_0005_euler_trim` | **粘性トリム版** Euler (`problem_d155_trim_tp.yaml`: Md 6.0144, L_c 39.7/M_K 2.7 [トリム後最短 x_F 96.00 r_t], r_t 0.076511) | dM_max 0.036 % M_d・x_E で M 6.01445 (=目標)・ε_M_rms 0.0056 %・コア 64/65。NaN 0 | active (トリム Euler) |
| `run_0006_ns_trim_v1` | トリム版 NS v1 (相関 δ\*, `problem_d155_trim_ns.yaml`, IC=run_0004 cross-mesh) | 完走・NaN 0・凍結 2.7e-4。出口面軸 M **5.9810** (= 旧 v1 5.9667 + トリム 0.24 % — トリムは線形に作用)。δ\* v3' 抽出元 | active (トリム v1・抽出元) |
| `run_0007_ns_trim_v3` | トリム版 NS v3 (CFD 抽出 δ\*, IC=run_0006) | 完走・NaN 0・凍結 1.2e-4・**δ\* 固定点 (2巡目比 median 0.999)**。**出口面軸 M = 6.00079 (+0.013 % M_d)** — トリム狙い通り。quasisteady machmax/pmax STEADY (shock 指標は無衝撃のため対象外)。出口壁半径 0.77591 m → 最終 r_t 補正 −0.12 % で 0.775 m 合わせ | active (**dry 最終形の正本**) |
| `run_0017_laminar_coarse` | **モデル梯子用 laminar NS** (`model: "none"` 化 + warm interp, coarse y+~50, 18000 step 段階起動) — cfl 上限スイープの probe 基地 | 完走・NaN 0。probe: cfl4 OK / 8@127 / 16@16 (SST・Euler と同一境界) | active (モデル梯子) |
| `run_0018_ns_trim_cfl6_r07` | **生産推奨設定の段階起動検証** (run_0007 と同一入力: IC=run_0006, dstar_v3+blend(−1,−0.5), 24000 step。差分は本段 `cfl 6 + implicitRelax 0.7` のみ) | 完走・NaN 0 (P min 2188 Pa)。**metrics が run_0007 と 5 桁一致** (出口軸 M 5.970306 vs 5.970316)。quasisteady machmax/pmax STEADY。残差は run_0007 の 24000 step 到達水準に 6.4k–11.3k step で到達 (2〜3.8 倍速)・roK/roOmega はさらに 5.8×/3.7× 深い。本段 192.9 s (cfl1 は 184.1 s = relax の追加コスト +5 %) | active (**生産設定検証の正本**) |
| `run_0019_ns_cm_pass1` | **排除厚さ更新計画 V1 pass 1** ([plan](../../plans/active/tooling-nozzle-deltastar-core-matched-euler.md)): pass 0 = `run_0004` の場から固定 Euler (`run_0001`) 基準で δ_r を抽出 (この pass はコア全体 α 方式)、ω=0.5、半径方向オフセット、**Md 6.0 トリムなし**、IC=run_0004 cross-mesh、`problem_d155_ns.yaml`、24000 step | 完走・NaN 0・品質 SOFT-PASS・quasisteady STEADY・残差 2 桁 warm 床 (NOT CONVERGED 判定 = 系列共通)。**ṁ_NS/ṁ_E 1.020 → 1.010** (スロート補正 0.0114 → 0.0065 r_t)、x≤30 の軸残差 −0.4 → −0.1 %。しかし試験部軸 M に +0.8 % (x 55–75) / −1.6 % (x 85–94) の波 = コア全体 α 抽出が上流壁誤差の波を欠損に取り込んだ (x=20–50 で δ +15 %)。→ 抽出器を帯局所参照に変更 (plan §4.3) | active (pass 1 記録・帯局所化の根拠) |
| **`run_0020_ns_cm_pass2`** | **V1 pass 2**: `run_0019` の場から**帯局所参照** (`band_local_deficit`) で δ_r を全域抽出、ω=1.0、半径方向オフセット、Md 6.0 トリムなし、IC=run_0019 | 完走・NaN 0・品質 SOFT-PASS・quasisteady ALL STEADY・残差 2 桁 warm 床 (系列共通)。**ṁ_NS/ṁ_E = 0.9999** (スロート補正 0.0014 r_t = 質量流量実効値)・**出口面コア M 6.0008〜6.0011 (+0.01 %) — トリムなしで目標達成**。軸 M は x=20–75 で Euler 設計 +0.05〜0.17 %、x≥80 の軸上に ±0.5〜1 % の波 (pass1→2 の壁変化 x≈40 が末端軸へ届いたもの、コア平均には出ない)。固定点差: use/in median 1.003、x=40–60 で +4 % → pass 3 で収縮確認 | active (**新チェーンの到達点 (暫定)**) |
| `run_0021_ns_cm_pass3` | **V1 pass 3**: `run_0020` の場から帯局所抽出、ω=0.5、IC=run_0020 (収縮確認) | 完走・NaN 0・SOFT-PASS・ALL STEADY・残差 warm 床。**ṁ_NS/ṁ_E = 1.0000、出口面コア M 6.005/5.996/5.998 (x=80/90/94, ±0.08 %)**、軸波 ±1 % → ±0.5 % に減衰。固定点差は x=20–60 で ±5〜7 % (帯幅依存の下流揺らぎ、帯パラメータ再検討中) | active (pass 3) |
| **`run_0022_ns_ib_pass0`** | **V2 pass 0: 積分法 (CONTUR 運動量積分, 断熱壁) 初期壁のみ、CFD 帰還なし** (`--init-integral`, 半径方向オフセット, Md 6.0 トリムなし, IC=run_0021 cross-mesh) | 完走・NaN 0・SOFT-PASS・ALL STEADY・残差 warm 床。**ṁ_NS/ṁ_E = 0.9991、出口面コア M 6.0026〜6.0037 (+0.04〜0.06 %)、軸 M は x=20〜94 で Euler +0.09〜0.21 % で平坦 (波なし)**。抽出 δ_r vs 積分法 δ_r: x≤60 で ±2 %、x=90 で −1.7 % (固定点にほぼ乗っている)。初期壁だけで主ゲート達成 | active (**積分法初期壁の到達点**) |
| **`run_0023_ns_ib_pass1`** | **V2 pass 1**: `run_0022` (積分法初期壁) の場から帯局所抽出 (適応帯)、ω=1.0、単調性ガード (λ=1.0 追加平滑化)、IC=run_0022 | 完走・NaN 0・SOFT-PASS・ALL STEADY・残差 warm 床。**ṁ_NS/ṁ_E = 0.9999、出口面コア M 6.0008〜6.0013 (+0.01〜0.02 %)、固定点差 use/in p10–p90 0.997〜1.002 (収束)**。軸 M は x=20–94 で Euler −0.14〜+0.45 % (x=70 に小さな山)。V1 pass 3 (相関初期壁から 3 pass) と δ_r が x=20–60 で 3 % 以内で一致 = **初期推定器に依らない固定点** | active (**新チェーンの正本 (積分法初期壁 + 1 pass)**) |
| `run_0024_ns_cm_pass4_q` | **δ_r 平滑化 (5 次 P-spline) の検証, V1 系 pass 4**: `run_0021` の場から帯局所抽出 → 5 次 P-spline (3 階差分ペナルティ, ノット 2 r_t, λ=1) → ω=0.5 (緩和後も再平滑化) | 完走・NaN 0・SOFT-PASS・ALL STEADY。**壁曲率の高周波ノイズ 7.7e-3 → 3.8e-4 [1/r_t]** (積分法初期壁と同水準)。ṁ 0.9997、出口面コア M ±0.1 %、軸 M の試験部最大偏差 0.49 % (run_0021 の 0.40 % と同等 = 軸の波は壁の凸凹ではない) | active (平滑化検証) |
| **`run_0025_ns_ib_pass2_q`** | **同 V2 系 pass 2**: `run_0023` の場から抽出 → P-spline → ω=1.0 | 完走・NaN 0・SOFT-PASS・ALL STEADY。**ṁ 1.0002、出口面コア M +0.01 %、固定点差 p10–p90 0.999〜1.006 (収束)**。軸 M は x=70 に +0.42 % (run_0023 の +0.45 % と同じ) — 滑らかな壁でも残るので、積分法壁→抽出壁の ~1 % の δ 差 (x≈15–25) への軸の集束応答。コア平均は平坦 | active (**生産形の最終 (平滑化込み)**) |
| `run_0026_ns_ib_pass3_nosoft_cfl5` | **起動レシピ A/B (a)**: run_0025 の固定点から `--stages none` (soft/mid なし) + **cfl 5 + implicitRelax 0.7**、24000 step | 完走・NaN 0・SOFT-PASS・ALL STEADY。metrics は run_0025 と同一 (ṁ 1.0002、出口面コア M 6.0000)、roK は 6e-5 と 3 倍深い。**全残差列が run_0025 (cfl1, 30000 step) の最終水準に 5000〜8900 step で到達、出口 M は 8000 step で凍結** → 12000 step で十分 (NS ≈ 95 s)。NS 壁時計 191 s (24000 step) | active (**warm start の生産レシピ根拠**) |
| `run_0027_ns_ib_pass3_ramp_cfl5` | **起動レシピ A/B (b)**: `--stages ramp --ramp 1,2,3.5 --ramp-steps 1000` → 本段 cfl 5 + relax 0.7 | 完走・NaN 0・SOFT-PASS。metrics は (a) と同一 (差 1e-4 以下)。ランプの利得なし (warm start では不要; cold start 用の保険として残す)。本段 step が outStepInterval の倍数でなく最終 res が 16000 止まり → runner 側で丸めるよう修正済み | active (A/B 記録) |
| `run_0028_euler_cfl6_r07` | **Euler の CFL 引き上げ**: run_0001 と同一設計、`--cfl 6 --implicit-relax 0.7 --stages soft` (soft 3000 + 本段 12000; run_0001 は soft+mid+本段 cfl2 18000 step) | 完走・NaN 0・`check_convergence` **ALL PASS**・ALL STEADY。metrics は run_0001 と同一 (dM 0.0363 % vs 0.0361 %、出口軸 M 6.00007 vs 6.00003、ε_M 0.0058 %)。残差は run_0001 最終水準に 4700〜8900 step で到達 → 本段 8000 step で可 (run_0032 で確認) | active (Euler 起動レシピ) |
| **`run_0029_ns_ib_pass0_norelay`** | **中継なし NS pass 0**: 積分法初期壁、IC = Euler `run_0001` を直接 y+~1 メッシュへ cross-mesh (k/ω 貼り + ω 床)、既定 `full` 起動 (soft/mid + 本段 cfl1 24000) | 完走・NaN 0・SOFT-PASS・ALL STEADY。**ṁ 0.9999、出口面コア M +0.00〜+0.04 %、固定点差 0.985〜1.002 = 中継あり (run_0022) と同等 → 中継 (y+~50 の踏み台 run) は不要**。NS 壁時計 514 s は別セッションの GPU 共有中の値 (通常 ~190 s) | active (**中継不要の根拠**) |
| `run_0030_ns_ib_pass0_norelay_cfl5` | 中継なし + soft/mid なし + 本段 cfl 5 (cold start を一気に) | **step 1 で NaN** → cold start には 1 次・低 CFL の soft/mid が必要。ディレクトリ削除済み (記録のみ、`_logs/run_0030.log`) | 削除済 |
| **`run_0031_ns_ib_pass0_norelay_full_cfl5`** | **cold start の短縮版**: 中継なし (IC = Euler run_0001)、`--stages full` (soft/mid 6000) + **本段 cfl 5 + implicitRelax 0.7、12000 step** | 完走・NaN 0・SOFT-PASS。**ṁ 1.0000、出口面コア M +0.00〜+0.03 %、固定点差 0.982〜0.998 = run_0029 (cfl1 24000) と同等**。outStepInterval 8000 が 12000 を割り切らず最終 res が 8000 止まり (prepare_ns 側で修正済み) → quasisteady はスナップショット不足で判定不能、残差は roK 2.8 桁 falling。NS 壁時計 411 s は GPU 共有中の値 | active (**cold start 生産レシピ候補**) |
| `run_0032_euler_cfl6_r07_8k` | Euler `--cfl 6 --implicit-relax 0.7 --stages soft --steps 8000` (run_0028 の短縮版) | 完走・NaN 0・品質 PASS。metrics は run_0001 と同一 (dM 0.0361 %、出口軸 M 6.00005、ε_M 0.0057 %) だが `check_convergence` は未達判定 (12000 step の run_0028 は ALL PASS) → Euler 本段は 12000 を推奨 | active (Euler 短縮の記録) |
| `run_0033_euler_Lpipe10` | 上流配管の影響テスト用 Euler: `problem_d155_R2_Lc39.3_tp_Lpipe10.yaml` (直管 L_pipe 0.5 → 10 r_t)、cfl 6/relax 0.7/soft | 完走・NaN 0 (run_0034 の固定 Euler 基準) | active (基準) |
| `run_0035_euler_rt77p05` / `run_0036_ns_final_rt77p05` | 最終 r_t 補正の 1 回目: r_t 77.05 mm (run_0025 の**生抽出** δ_r(x_F−0.3)=0.684 で解いた値) の Euler + NS (run_0025 から warm start, none/cfl5/12000) | 完走・NaN 0・SOFT-PASS・ALL STEADY。ṁ 1.0002、出口面コア M 6.0000、固定点 0.996〜1.002。**物理出口半径 0.7753 m と 0.3 mm 残った** = 壁に載るのは平滑化後の δ_r,exit 0.688 なので、solve_rt をその値を使うよう修正して r_t 77.02 mm でやり直し (run_0037/0038) | active (補正 1 回目の記録) |
| `run_0034_ns_ib_pass0_Lpipe10` | **上流配管の影響**: 直管 L_pipe 0.5 → 10 r_t (入口境界層が 6 倍厚い)、他は run_0031 と同じ (中継なし, full, cfl5, 12000; ni 1400) | 完走・NaN 0・SOFT-PASS・ALL STEADY。**スロート以降の δ_r は L_pipe 0.5 (run_0031) と 0.3 % 以内で同一 (x=0: 0.0014/0.0014, x=90: 0.659/0.661)、出口面コア M +0.00〜0.03 % も同一**。ṁ 比 0.9992 (スロート実効 δ 0.0018 vs 0.0014 = 0.03 mm)。→ 収縮部の加速で上流の境界層履歴は消える (積分法の θ0 無記憶と整合) | active (**配管非依存の根拠**) |
| `run_0037_euler_rt77p02` | 最終 r_t 77.02 mm の Euler (固定基準, cfl 6/relax 0.7/soft) | 完走・NaN 0・PASS・**ALL PASS**・ALL STEADY | active (最終形の基準) |
| **`run_0038_ns_final_rt77p02`** | **最終形**: r_t 77.02 mm (solve_rt: run_0025 の壁に載せた δ_r,exit 0.688 で 0.775 m 合わせ)、δ_r = run_0025 抽出 (P-spline)、warm start (none/cfl5/relax0.7/12000) | 完走・NaN 0・SOFT-PASS・ALL STEADY (残差 warm 床)。**物理出口半径 0.7750 m (spec 一致)、ṁ_NS/ṁ_E 1.0002、出口面コア M 6.0002 (+0.00 %)、固定点 0.993〜1.002、全長 (スロート→出口) 7.325 m、物理スロート半径 77.14 mm**。点列 `points_d155_final_{ns,euler}.csv` | active (**最終形の正本**) |
| **`run_0039_ns_final_cond`** | **最終壁 (run_0038) の凝縮 ON restart** (`problem_d155_ns_rt77p02_cond.yaml`: Kw+HK condModel1+Kantrowitz, 蒸発 ON; IC=run_0038 同一メッシュ, none/cfl1, 12000 step; 評価 `eval_cond.py`) | 完走・NaN 0・SOFT-PASS・**ALL STEADY (4k/8k/12k で同一)**。軸 onset x≈57.9 r_t、出口 (x=94.1) g 0.27 % (H₂O の 3 %)、S_max 17.4。**出口軸 M 5.918 (−1.36 %)、出口面コア M 5.989 (−0.18 %)** (dry run_0038: 6.0017 / 6.0002)。旧トリム壁 run_0008 (onset 68.6, 軸 −0.99 %, コア −0.17 %) と同じ結論 = Tt 1600 K では凝縮不可避、影響は軸に集束して見えコア平均は −0.2 % | active (**最終形の凝縮評価**) |
| `run_0008_ns_trim_cond` | 凝縮 ON restart (`problem_d155_trim_ns_cond.yaml`: Kw+HK condModel1+Kantrowitz, 蒸発 ON, IC=run_0007, 12000 step) | 完走・NaN 0・**STEADY** (4k/8k/12k で M_exit 差 5e-4)。軸 onset x≈69 r_t、出口 g 0.20 % (H₂O の 2 %)、**出口軸 M 5.9273 (−1.2 %)**・試験区間に M 低下勾配 (x60→96 で 6.00→5.93)。dry の軸は x≈24 r_t (M5.5) で飽和線越え S≈14 (`axis_values.csv` の Tsat_post) | active (**凝縮評価の正本**) |

## SU2 クロスチェック (境界層厚さ, 2026-09-01)

δ\* 抽出の軸対称化 (円環重み) と合わせて、**同一メッシュ・同一 BC・同一ガス (CPG γ\*=1.27354)** で forge node SST と SU2 v8.5 axisym SST の境界層厚さを比較する ([procedures/su2-cross-check.md](../../procedures/su2-cross-check.md) 準拠)。

| run | 内容 | 結果 | 状態 |
| --- | --- | --- | --- |
| `run_0009_cpg_euler` | CPG 版 Euler (`problem_d155_cpg_euler.yaml`) | 完走・NaN 0。CPG は A/A*=159 (semi-perfect 88.2) でノズルが伸びる (出口壁半径 1.03 m) — 比較チェーン専用で製品形状ではない | active (比較チェーン) |
| `run_0010_cpg_ns_coarse` | CPG NS 中継 (y+~50) | 完走・NaN 0 | active (中継) |
| `run_0011_cpg_ns` | **CPG NS 本計算** (wall_first 6.5e-5=y+~2 [4.5e-5 は CPG 長尺で AR1403 FAIL→緩和], 相関 δ\* 物理壁, 24000 step) — 比較の forge 側 | 完走・NaN 0・品質 PASS。残差 2 桁 warm 床 (既知パターン) | active (**比較 forge 側**) |
| `run_0016_chamber_cpg` | **チャンバー A/B** (Codex レビュー準拠): 完全形状 (入口管+収縮+ノズル) + 出口後方 4D_e×外径 3D_e チャンバー+フランジ壁の 4 ブロック transfinite (`gen_chamber_mesh.py` → `mesh_chamber/`)。品質 PASS (AR296)。静止雰囲気 1389 Pa IC + 段階起動 18000 step で plume 確立 (ジェット軸 M5.97, 外周静止)。**probe (cfl4 OK/8@181/16@21) は素メッシュと同一挙動 = 出口境界は無罪**。**line-implicit 併用の 4 象限も全て同一** (coarse×line 8@216/16@15, chamber×line 8@194/16@15; チャンバーのライン被覆 81.9%・噴流域 point fallback) — 壁法線×出口境界の直積消去で streamwise 内部モード律速が確定 | active (**A/B の正本**) |
| `run_0012_su2_sst` | **SU2 v8.5 RANS-SST** (同一メッシュ msh→su2 変換, axisym, ROE+MUSCL, SST V2003m, Sutherland/Pr0.72/Prt0.9, 40k+20k iter) | 収束: rms[Rho] −6.5・rms[RhoE] −2.3・出口 massflux/Mach ドリフト 0.003 %/5k (手順書合格実績以上)。`history.csv` は継続分のみ | active (**比較 SU2 側**) |
| `run_0013_cpg_ns_plainsst` | forge 素 SST 化 A/B (`dilatationCorrection: 0, katoLaunder: 0` — prepare 後に solverConfig 直接編集。YAML への独自キー追記は黙って無視されるため不可), IC=run_0011, 12000 step | 完走・NaN 0。**素 SST にすると SU2 と一致: δ99 ≤3 %・θ ≤0.4 %・δ* ≤1.3 %** | active (**帰属 A/B の正本**) |
| `run_0014_cpg_ns_dilat0_kl1` | 帰属分離 (dilatation OFF / KL ON) | run_0013 と ≤0.7 % 差 → **Kato-Launder は無関係、差は dilatationCorrection 単独** | active (帰属分離) |

### SU2 クロスチェックの結論 (2026-09-01)

- **forge 生産設定 SST は SU2 比で境界層が薄い**: δ99 −17〜21 %・θ −16 %・δ* −5 % (4 ステーション一貫)。
- **帰属は `dilatationCorrection: 2` (圧縮性生産項の正確形: deviatoric trace 除去 + 等方項 −⅔ρk∇·u) 単独**。素 SST 化した forge は SU2 と δ99 ≤3 %/θ ≤0.4 %/δ* ≤1.3 % で一致 = **ソルバ・離散化は無罪、意図した乱流モデル差**。Kato-Launder の寄与 ≤0.7 %。
- ノズルは全域 ∇·u>0 (強膨張) のため等方項が k のシンクとして常時働く。急膨張の乱流抑制は物理的に実在する効果で forge 形は正当だが、**SST のモデル形式差として境界層厚に ±16 % 級の不確かさ**があると認識すること。δ* への影響は ±5 % (出口 M ±0.1〜0.15 %) に留まり、Md トリムは forge 自身の NS で較正しているため製品性能の結論は不変。
- 比較図・数値: `compare_bl_su2.{png,json}` (3 者重ね描き)、40k 時点の記録 `compare_bl_su2_iter40k.json`。

## 陰解法 cfl_pseudo 上限スイープ (2026-09-02)

`run_0015_cflsweep/` (probe 22 本, `sweep_cfl_implicit.py` + `summary.json`。中間 res は削除済み・各 probe の config/residual/NaN ダンプは保持)。warm 場 (run_0007 TP / run_0011 CPG) restart 2000 step + coarse IC 収束レース 4000 step。

- **素の上限 (implicitRelax=1): TP は cfl 2 (3 で発散)、CPG は 3 (4 で発散)**。既定挙動としては単桁前半で頭打ち。
- **原因は低マッハ調整ではない** (仮説棄却): ① 現行 config は `lowMachPrecond: 0` で低マッハ処理は不活性、② 発散箇所は亜音速チェンバでなく**下流端の超音速壁 BL 内** (壁から ~0.003 r_t の対数層)、③ `lowMachPrecond: 2` を入れると逆に悪化 (step 97→13)。
- **終端症状 = 陰的更新の P アンダーシュート → EOS 床洗浄 → NaN** (ω 爆発は増幅表示)。**帰属の最終確定 (2026-09-02, チャンバー A/B 完了)**: 発散種の位置はメッシュ依存 (fine=壁∩出口コーナー / coarse=出口手前コア / 高cfl=スロート軸上) で、①背圧値 100/10 Pa・②逆流 Tt 1500K・③出口近傍局所 dt キャップ・④**チャンバー付加 (出口境界を 8 m 後方へ, run_0016)** の全てに**不感** (発散 cfl・step ほぼ不変: chamber cfl4 OK/8@181/16@21 vs 素 cfl4 OK/8@215/16@17)。⑤ **乱流状態更新の凍結 (`FORGE_FREEZE_TURB=1`) でも発散 step が 1 step 単位で一致** (cpg@98/28/18, tp@53/24)。**確定したのは**: k-ω segregated 更新は発散に不要 / ω 先行 NaN は症状 / 流れブロックだけで発散する。**未消去**: μt は凍結 roK,roOmega から毎 step 再計算される (dependentVariables/turbulent_viscosity) ため SST closure・μt を含む粘性剛性の寄与は残る。⑥ **laminar NS (`model: "none"`, run_0017 で定常化後 probe)**: cfl4 OK / 8@127 / 16@16 — SST と同じ上限。⑦ **Euler (run_0009 base, 粘性・乱流・BL なし)**: cfl4 OK / 8@150 / 16@55 / 32@19 — **やはり同じ 4→8 境界** (種は x 4–12 r_t の近軸内部 = 最急膨張域)。⇒ 判定木の「Euler まで同じ」分岐に着地: **律速は対流の streamwise defect-correction 反復不安定** (SST closure・μt・粘性剛性は全て消去済み)。根治は streamwise ADI / LU-SGS 流下順序が妥当。
- **切り分け**: `nStepInner` 10/20 無効 (DPLUR 内部収束でない)・`implicitRelaxSST` 0.5 無効 (SST 方程式でない) — **`implicitRelax` (流れ方程式の緩和) だけが効く**。
- **緩和込みの上限**: relax0.7 → cfl 8 (12 で発散)、relax0.5 → cfl 16+。積 cfl×relax ≈ 6-8 で飽和。
- **収束レース**: 実効速度も cfl×relax でスケールし **cfl8+relax0.7 が最速** (4000 step で rms_roUx 現行 cfl1 比 1/3、rms_roOmega 1/4)。cfl16+relax0.5 は頭打ち。図 `run_0015_cflsweep/race_residuals.png`。
- **生産推奨**: CPG `cfl 8 + implicitRelax 0.7` / TP `cfl 6 + implicitRelax 0.7`。現行 cfl1 比で 2〜4 倍速。**TP cfl6+r0.7 は実段階起動系列で検証済み** (`run_0018_ns_trim_cfl6_r07`: run_0007 と metrics 5 桁一致・machmax/pmax STEADY・NaN 0。素の TP cfl3 は発散するため **implicitRelax 0.7 とのペアで固定**すること)。停止条件 (十分速く安定に答えが出る) を満たしたため **ADI/LU-SGS の実装は棚上げ** — 速度が再びボトルネックになったら着手 (GPU 可否は検討済み: 第2ライン族 ADI は v1 機構流用で可、古典 LU-SGS は 2D 波面並列が細く不適 = DPLUR の存在理由)。
- **line-implicit 追試 (2026-09-02, probes7–10, [plan](../../plans/accepted/time_integration-line-implicit.md))**: 壁法線ライン block-Thomas (`lineImplicit: 1`) を実装して検証。ライン構築は完璧 (1250 本・被覆 100 %) だが **cfl 上限は不変 (積 ≈6–8)** — cfl8-line の発散は出口域全断面で、**律速は壁法線でなく streamwise lag** と判明 (壁法線剛性は粘性対角が既に吸収)。効果は同一設定の収束 −14 %/step のみ。実装過程で lu5 ピボット罠 (LASWP 先行必須)・device printf 引数上限罠を踏んで記録。
- **追試 (2026-09-02, probes5/6)**: commit 時の正値性ガード (`updateGuardAlpha`, 新実装 opt-in) は**上限を上げない** (発散遅延のみ。α=0.5 は毎 step 半減を許すため床への歩行が続く)。**lowMachPrecond=2+ガード併用も不成立** — P が床に触れないまま SST チャネル (ω 2.8e22) で爆発。⇒ 真因は「defect-correction 反復 (近似 Jacobian + SST 隔離更新) の不安定」で、床 NaN は終端症状。`implicitRelax` が効くのは反復スペクトル半径を縮めるから。詳細 [plan](../../plans/accepted/time_integration-update-positivity-guard.md)

## 結論 (2026-08-31)

- **最短全長: スロート→物理出口 7.337 m** (x_F 96.00 r_t, r_t 76.42 mm)、入口配管端→出口 8.292 m。うち E→F 一様化区間 ≈4.36 m は物理の床。
- **出口軸 M (dry): 6.0008 (+0.013 %)** — Md 6.0144 トリム (NS 固定点の −0.24 % 欠損を打ち消し) で達成。
- **凝縮リスク**: Tt 1600 K では dry 不成立 (S≈14)。凝縮 ON で出口 M 5.927・軸勾配残り。凝縮フリーは Tt ≳ 1820 K。
- 結果ページ: https://claude.ai/code/artifact/f1cfbdf8-2415-4bfc-ae2d-156793569bd0

注: `run_0006_ns_trim` (旧設計 run_0004 の δ\* CSV を新設計に流用するショートカット) は物理壁フィルタ不合格 (スロート下流非単調) で prepare 段階に失敗し**削除済み**。トリム版も正規フロー (v1→抽出→v3) で回す。
