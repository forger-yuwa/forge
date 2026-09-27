# 諮問ブリーフ: Graetz plan の codex plan レビュー (C0/M4/m2) の採否 + 起動時の擬似 CFL (2026-09-27)

AGENTS.md エスカレーション条件 **5** (codex Critical/Major の採否)。plan: `plans/active/boundary-cht-axisymmetric-graetz.md`、
レビュー: `notes/reviews/2026-09-27-boundary-cht-axisymmetric-graetz-plan.md`。前回の諮問 `notes/reviews/2026-09-27-graetz-design-diagnose.md` は全件採用済み。

## 0. 読んでよいファイル
- 上の 2 つ + plan。評価器 `case/63.graetz_cht/eval_graetz.py`、run 生成 `case/63.graetz_cht/make_run.py`、テンプレート `case/63.graetz_cht/template/solverConfig_{iso,cht}.yaml`
- 判定ツール `solver_density_cuda/tools/check_cht_balance.py` (usage 冒頭)、`check_cht_interface.py`
- スモーク run (登録 run でない。ローカル float・未収束): `case/63.graetz_cht/run_000{2..6}_*`

## 1. 採否案 (私の案。全件採用)

| 指摘 | 案 | 具体 |
| --- | --- | --- |
| M1 x=0 断面の温度一様は正しい解でも FAIL | **採用** | (a) 加熱 run の**入口面 x=−L_up** の T(r) max−min ≤ 0.05 K、(b) **対照 run (ΔT=0) の x=0 断面** (壁節点を除く) の T max−min ≤ 0.05 K (予熱が無いので入口 BC の非一様だけを見る。評価器は既にこれ)、(c) 加熱 run の x=0 断面 T(r) は ellip (Pe 720) の場と**参考比較** (合否に使わない) |
| M2 V-g3 がモデル差込みの誤差の単調減少 | **採用** | V-g3 を「共通位置 (x⁺ = 3e-3, 1e-2, 3e-2, 0.1 と窓内の全 N_r=16 壁節点位置) での**格子間差** \|Nu_16−Nu_32\| > \|Nu_32−Nu_64\| (各点と最大)」に変更。観測次数は参考値のみ、漸近域とは扱わない。**固体半径 8 層**の感度: N_r=32・ΔT 10 の共役で 16 層の run を 1 本足し、差し引き Nu の差 ≤ 0.05 % (窓内最大) を登録 |
| M3 欠損・局所ドリフトを見逃す | **採用** | 窓内の期待節点集合をメッシュ (流体の壁節点で x⁺ が窓内) から固定し、壁ダンプに欠損・重複・`iface_ok`=0 があれば **REFUSED** (終了コード 2)。準定常は窓内**全節点**の分子・分母・Nu の系列 (各節点) を `check_quasisteady.py --series-csv` に渡し ALL STEADY、加えて末尾半分の各節点の Nu の (max−min)/mean ≤ 0.03 %。負例に「最大誤差節点を壁ダンプから落とす → REFUSED」「4 観測点を避けた 1 節点の局所ドリフト → FAIL」を足す |
| M4 G-if/G-cons の数値が一意でない | **採用** | §6 に全パラメータを書き、テンプレートをその値にする。**G-if** (`conjugate.gate`、全域): `eps_abs_Wm2` **1.0** (= 窓内で最小の q ≈ k_f·3.657·ΔT_b/D ≈ 209 W/m² (ΔT 10、x⁺ 0.1 で ΔT_b ≈ 2 K) の 0.5 %。親 V-ax2 の「0.5 % of q*」と同じ作り方) / `eps_rel` **1e-3** / `dT_K` **1e-3** (= 1e-4 × ΔT 10) / `tol_solid` **1e-9 W/rad** (窓内の節点荷重 q·Δx·R ≈ 2e-4 W/rad の約 1e-5) / `n_consec` **80**。**対照 run (ΔT=0) も同じ絶対値**を使う (対照の q0 は加熱 run の q から引かれるので、その誤差は加熱 run の熱量尺度で測るべき)。**G-cons** (`check_cht_balance.py --solid-mode fem2d --phys-id 4 --phys-name wall_heat`): `--q-floor` **8.4e-5 W/rad** (= 1e-3 × 加熱 run の総入熱の見積もり ṁ' c_p ΔT (1−θ_b,end) ≈ 1.02e-5·1004.5·10·0.82 ≈ 0.084 W/rad)、`--tol-rel` **1e-3**、`--tol-abs` **8.4e-5 W/rad**、対象は最終スナップショット、対照 run も同じ値 |
| m5 固体の T_c と ΔT の一致 | 採用 | `make_run.py` で solid.h5 の `ROBIN/TC` (実際にソルバが読む値) を読み、T_in+ΔT と不一致なら停止、実値を RUN_INPUTS に記録 |
| m6 主張の範囲 | 採用 | §1 目的・完了時の主張を「登録した node・FP64 条件で、**対照差し引き後の局所 Nu** が許容差内」に。§4.4 の「連成が入れた誤差」を V-g4 と同じ「UWT 近似を含む差」に統一。加熱・対照に共通の加算誤差は消えるので絶対熱流束の精度は保証しない、と明記 |

## 2. 追加の観測事実 (レビュー後に判明、plan 未記載)

スモーク (ローカル float、N_r=16、Poiseuille IC、`cfl_pseudo` 5・`implicitRelax` 0.7):
- `run_0003_smoke_iso0_r16` (ΔT=0、等温壁 300 K): step 14 から x=0・r=0.9375 mm (壁の隣の内部節点、最も細かい軸方向セル 41 µm) で指数的に崩れ、step 100 で P −650 kPa・|u_r| 161 m/s。
- `run_0004_smoke_walls_adiab_r16` (壁を全部断熱): 崩れない (200 step で max|u_r| 1.3e-3 m/s)。`run_0004_smoke_walls_iso_r16` (全部等温 300 K): 同じく崩れる → 断熱/等温の継ぎ目でなく**等温壁 (強制形) そのもの**。
- `run_0005_smoke_iso0_cfl{2,1}_r16`: `cfl_pseudo` 2 / 1 は 400 step まで安定 (max 偏差 0.026 / 0.015 m/s)。
- 既知: 強制形の等温壁は擬似 CFL の上限 ~5 (`plans/active/boundary-weak-isothermal-wall.md`)。
- **案**: 本番の `cfl_pseudo` を **2** (`implicitRelax` 0.7 のまま) に固定して登録する。発散手順の 1 回目の対処で解けたので条件 2 には当たらないと判断。

## 3. 諮りたいこと (推奨を 1 つに絞る)
1. 上の採否案で良いか。数値 (G-if/G-cons、固体 16 層の 0.05 %、Nu 変動 0.03 %) の作り方に穴はないか。
2. `cfl_pseudo` 2 の固定で良いか (N_r=64 は軸方向の最小セルが 10 µm と更に細かい。3 格子で同じ CFL にすべきか)。
3. 見落とし。

**両論併記で逃げず、推奨を 1 つに絞ること。根拠は `ファイル:行` か上の数値で示すこと。**
