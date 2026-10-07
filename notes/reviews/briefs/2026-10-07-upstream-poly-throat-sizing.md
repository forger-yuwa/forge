# 諮問: 物理壁の上流を多項式 1 本にする (δ_r はスロートから下流だけ) / 寸法をスロート径から逆算する — plan §4・§6

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 1 (plan §4・§6 を新規に書く)。
plan: `plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md` (全文を読むこと)。作業ツリー `/home/sano/work/forge-integ-1005`。
同じツリーで別の実装担当が `design/forge_design/geometry/wall_axismach.py`・`evaluate/runner_axismach.py`・`report/nozzle_report.py`・`methods/design/overview.md` を編集中 (全域 1 本の B-spline の plan)。未 commit の差分が混ざるので、HEAD の版を基準に読むこと。

## 観測事実

- **今の物理壁** (joint 壁の解析経路、`wall_axismach.py:497-500`): r_W = r_design + s(x)·δ_r(x)。r_design は直管 6.485 ([−12.5, −12))・上流 Hermite H ([−12, 0))・S ([0, 95.245])。s は `geometry.pw_ramp` [−11, −6] の 5 次 smoothstep。ランプのゲートは \|r_W″ − r_design″\| ≤ 5e-3 と r_W′ < 0。
- **δ_r の大きさ** (run_0147、積分法 CONTUR、k_f = 1.054): x = −11 で 0.48 mm、−9 で 0.73 mm、−6 で 0.57 mm、−3 で 0.22 mm、スロート 0.107 mm (0.00139 r_t)、出口 56 mm (0.733 r_t)。
- **事前試算** (`case/45.isobutane_m6_d155/upstream_poly_probe.py` → `_band_ab/upstream_poly_probe.txt`、CFD 0 step):
  - 上流の 5 次多項式 Q (x = −12 で (6.485, 0, 0)、x = 0 で下流の物理壁の (r, r′, r″) = (1.0013897, 7.76e-4, 0.501184))。
  - 物理スロート (r′ = 0) は今と同じ (x = −0.1188 mm、r = 1.0013891 r_t)。
  - 今の物理壁との差: 最大 −0.52 mm (x = −7.1、Q が細い)、1 r_t 以内で ≤ 6 µm。上流の r′ ≤ 0 は成立。r″ の今との差の最大 3.1e-3。
  - 徐変案 (スロートを動かさず下流で δ_r を 0 から全量へ): 実効スロート面積 −0.28 %、実効的な壁 (物理壁 − 排除厚) の設計壁からのずれの傾きは徐変長さ 1〜20 r_t で 0.12〜0.23°。
- **Euler の評価は設計壁を使う** (`runner_axismach.prepare` は `d["wall"]` でメッシュを作る)。
- **寸法の今の決め方**: `deltastar_loop.solve_rt` が出口半径から r_t を解く (CFD 前は積分法、NS 後は抽出した δ_E を Re^−0.2 換算)。結果は `spec.r_throat` に手で書く。
- 生産の NS (run_0147) の段階判定は main 区間で NOT CONVERGED (stalled/plateau)。単調壁の plan はこの条件で量の判定を限定して採用した (`plans/accepted/tooling-nozzle-throat-monotone-r2.md`)。

## ユーザ決定 (2026-10-07)

- 「スロート点より上流に排除厚さを足しこむ処理を、やめにしませんか」「上流側は、配管〜スロート点まで、の多項式で形状を作る形」→ 上流の多項式化を標準にする。
- 「時にはスロート径を固定したい場合もある」→ 徐変案ではなくスロート径から r_t を逆算する案 (「b 逆算にしようか」)。

## 方針 (plan §4)

- §4.1 `geometry.pw_upstream: ramp | poly`。**コードの既定は `ramp`** (過去の YAML と `rerun_conditions` の再現のため)。標準の手順を `poly` にする。`poly` は x = 0 (設計スロート) で下流の物理壁に 2 階微分まで接続する 5 次多項式。ゲートは上流の単調性、\|Q″ − H″\| ≤ 5e-3、継ぎ目の跳び ≤ 1e-8。
- §4.2 `solve_rt_throat`: 物理スロート半径 r_t·ρ_t(r_t) = R_throat を r_t について解く。CFD 前は積分法の不動点反復、NS 後は δ_E のスロート値を Re^−0.2 換算。

## 問い

1. 上流に δ_r を足さないことで、物理的に何を失うか。縮流部の実効的な壁は設計の H から最大約 0.5 mm 外れる (流路が狭い側)。亜音速の縮流部でこれが遷音速域・スロートの音速線・CFD でピン止めしている初期線 (run_0062 の Euler 由来) に与える影響の見込みと、確かめ方。
2. 接続位置を x = 0 (設計スロート) にする案と、物理スロート (r′ = 0) にする案のどちらがよいか。x = 0 で接続すると物理スロートは多項式の中 (−0.0015 r_t) にある。
3. コードの既定を `ramp` に残し、標準の手順だけを `poly` にする判断は妥当か (ユーザの言う「既定」をこの形で満たすことの是非)。`poly` で `pw_ramp` が書かれていたら例外か警告か。
4. §4.1 のゲートの値 (\|Q″ − H″\| ≤ 5e-3 は今のランプのゲートの流用) は妥当か。短い縮流部 (標準の L_U = 3.5、r_U = 2.5) で Q が単調でなくなる・曲率が崩れる可能性と、その時の扱い。
5. `solve_rt_throat` の定式化 (物理スロートの定義・δ のどの値を使うか・NS 後の Re^−0.2 換算) と U2 の合格条件は妥当か。スロートの δ は相関で 3〜12 倍過大になった前例がある (積分法の初期壁; 今は k_f 較正後でスロートの積分法 0.00139 r_t、NS 抽出 0.00144 r_t)。スロート径から決める場合、どの δ を使うべきか。
6. U4 (case/45 の生産に入れる場合の NS・凝縮の再評価) の方法: 初期値の移し方 (縮流部の節点が最大 0.5 mm 動き、前回の番号写像の上限 1 µm を超える)、比較する量と合否の登録。MOC の軸処理の plan (`discretization-moc-axis-limit-and-corrector.md`、V5′ で NS・凝縮を回す予定、生産候補にするかはユーザ判断待ち) との順序・まとめ方 (別々の A/B にするか、片方を先に入れるか)。

## 読んでよいもの

- 上記 plan、`plans/active/tooling-nozzle-wall-single-bspline.md`、`plans/active/discretization-moc-axis-limit-and-corrector.md`、`plans/accepted/tooling-nozzle-throat-monotone-r2.md`
- `design/forge_design/geometry/wall_axismach.py` (`PhysicalNozzleWall`・`default_pw_ramp`・`UpstreamThroatPoly`)、`design/forge_design/evaluate/runner_axismach.py` (`prepare`・`prepare_ns`・`integral_delta_r`)、`design/forge_design/feedback/deltastar_loop.py` (`solve_rt`)
- `methods/design/overview.md` (joint 壁の物理壁・物理スロート A13)、`procedures/nozzle-design-workflow.md`
- `case/45.isobutane_m6_d155/upstream_poly_probe.py`・`_band_ab/upstream_poly_probe.txt`
- run_0147 の保存物は別ツリー `/home/sano/work/forge/case/45.isobutane_m6_d155/run_0147_ns_mono_final/`

編集は禁止。
