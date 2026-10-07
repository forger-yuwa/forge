# 諮問: 物理壁を全域 1 本の x の 5 次 B-spline にする方針 (§4) と検証計画 (§6)

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 1 (plan §4・§6 を新規に書く)。
plan: `plans/active/tooling-nozzle-wall-single-bspline.md` (全文を読むこと)。作業ツリー `/home/sano/work/forge-integ-1005`。
別の実装担当が同じツリーで `design/forge_design/geometry/moc_*.py`・`evaluate/runner_axismach.py` を編集中 (MOC の plan)。読むのは自由だが、未 commit の差分が混ざっていることに注意。

## 観測事実

- **今の物理壁** (r_t 単位、単調壁の生産問題 `case/45.isobutane_m6_d155/problem_d155_ns_finemesh_recal_final_mono.yaml`、run_0147):
  - 6.485 ([−12.5, −12))。
  - H(x) ([−12, −11))、H + s·δ_r ([−11, −6))、H + δ_r ([−6, 0))。H は 5 次 Hermite、s は 5 次 smoothstep。
  - S(x) + δ_r(x) ([0, 95.245])。S は 5 次 B-spline で係数 259・区間 254。
- **δ_r の 2 段**:
  - 第 1 段 (`metrics/deltastar.py::smooth_delta_quintic`): ξ = √(x + 16.5) の等間隔ノットの 5 次 B-spline、係数 71、係数の 3 階差分の罰則 λ = 1。残差は相対 rms 1.9 %、最大 3.5e-4 r_t。
  - 第 2 段 (`runner_axismach.delta_r_from_table`): 第 1 段を 1500 点で評価し、それを通す 5 次補間スプライン。
- **第 1 段を直接使った場合との差** (CFD 0 step、ローカル): 値 6.3e-11 r_t、1 階微分 6.3e-9、2 階微分 3.5e-7、3 階微分 7.7e-6 (差分近似の誤差を含む)。δ_r の最小は 0.0011 r_t (`positive` の切り落としは発火していない)。
- ユーザ決定: 「設計側の平滑化そのものを x の B-spline に替える」「共通化させたい」「(i) 全域 1 本」。

## 方針 (plan §4)

- δ_r を x の 5 次 B-spline で平滑化する。ノットは今の ξ の等間隔ノットを x に戻した位置、罰則は λ∫(δ_r‴)² dx (Gauss で厳密)、λ は今と有効自由度が等しくなる値。
- 物理壁を 1 本の 5 次 B-spline に組み立てる。
  - 5 次の区分多項式になる区間は、ノットの和集合の上で厳密に書く。継ぎ目は重複ノットで 2 階微分まで連続にする。
  - ランプ区間 (10 次) は、両端の値・1〜2 階微分を拘束した最小二乗で、0.1 µm 以下に当てはめる。
- 新しいキー `geometry.physical_wall_repr` で選び、既定はビット同一。

## 問い

1. δ_r の x 空間の平滑化の設計 (ノット・罰則・λ の選び方) は妥当か。今の ξ の P-spline と比べて、何が変わりうるか (δ_r と物理壁の変化の見込みの大きさ)。
2. 全域 1 本の組み立て。厳密区間の和集合のノットと重複度、ランプ区間の当てはめ (拘束・許容誤差 0.1 µm) は妥当か。ランプを当てはめにすることで、継ぎ目や r″ の滑らかさに問題が出ないか。
3. §6 W1 の合格条件 (残差 rms ≤ 今 × 1.1、物理壁の r″ の高周波 ≤ 今、S6、出口半径の変化 ≤ 0.01 mm) は、結果を見る前の登録として十分か。C2 の較正 (k_f・r_t) を解き直す条件は何にすべきか。
4. 生産に入れる場合、NS・凝縮の再評価を MOC の plan とまとめて 1 回にすること (両方の変更を入れた壁で) は妥当か。変更が 2 つ交絡するのは問題にならないか。
5. CAD に渡す書き出し (ノット・係数) の検査として W5 で十分か。

## 読んでよいもの

- 上記 plan、`plans/active/discretization-moc-axis-limit-and-corrector.md`、`plans/accepted/tooling-nozzle-throat-monotone-r2.md`
- `design/forge_design/metrics/deltastar.py` (`smooth_delta_quintic`)、`design/forge_design/geometry/wall_axismach.py` (`PhysicalNozzleWall`・`joint_fit_wall`)、`design/forge_design/evaluate/runner_axismach.py` (`integral_delta_r`・`delta_r_from_table`・`prepare_ns`)
- `methods/design/overview.md`

編集は禁止。
