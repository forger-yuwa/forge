# 軸対称風洞ノズル (axismach) の設計・計算手順

問題タイプ `wind_tunnel_axisym_axismach` (軸 Mach 則 → 逆 MOC で壁を作る風洞ノズル) を対象に、
よくある 3 つのユースケースの作業手順をまとめる。最新の実例は `case/45.isobutane_m6_d155`
(M6、出口直径 1.55 m、イソブタン燃焼ガス、2026-10 時点の生産レシピ)。

- 関連: 要件インタビュー [`/design-intake`](../.claude/skills/design-intake/SKILL.md)、
  標準出力と pptx [`nozzle-design-outputs.md`](nozzle-design-outputs.md)、
  計算の一般手順 [`calculation-workflow.md`](calculation-workflow.md)、
  発散対処 [`divergence-and-startup.md`](divergence-and-startup.md)、AWS 投入 [`/forge-aws-run`](../.claude/skills/forge-aws-run/SKILL.md)。
- 設計判断の正本: [`plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md`](../plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md)
  (CFD ピン・当てはめ壁・細分格子・Euler 較正、§5.1 #11〜#11f)、
  [`plans/active/tooling-design-problem-campaign-recipe.md`](../plans/active/tooling-design-problem-campaign-recipe.md)
  (再現可能なパイプライン化、§4.9 C2 方式 — **未実装部分が多い**)。
- コマンドは `design/` 直下で `.venv-opt/bin/python` を使う。run は検証 run でも AWS で回す。
- **計算環境 (2026-10-07)**: 積分法の δ_r (scipy の RK45) は numpy・scipy の版で 1e-6 r_t の桁 (case/45 の出口付近で約 0.5 µm) 動く。
  設計として問題になる大きさではないが、壁・メッシュの**ビット同一を要する照合** (保存した壁の作り直し・W3/W4 型の入力の同一性・生産の
  問題の準備が生産の run の入力を再現するかの確認) は、**照合する相手と同じ環境で**行う。
  - 生産の準備と照合は、生産の run を回す AWS の環境で行う (例: `case/45.isobutane_m6_d155/prod_confirm.py`)。
  - 手元では `design/.venv-opt` の Python を使う。システムの `python3` は numpy・scipy の版が違い (numpy 1.26 には `numpy.trapezoid` も無い)、
    手元の過去の成果物とも一致しない。
  - `prepare`・`prepare_ns` は `prepare_info.json` の `environment` (Python・numpy・scipy の版、実行ファイル、ホスト、逆 MOC の
    backend と numba・llvmlite の版) に環境を残す。照合が合わないときは、まず両方の `environment` を比べる。
  - **numba が要る (2026-10-15)**: 逆 MOC の充填 (`InverseMOC.fill_levels`・`fill_arrays`) の 1 段の単位過程は既定で numba の版
    (`geometry/moc_fill_nb.py`。numpy の版とはビット同一でなく、網の位置の差は 1e−14 r_t の桁、受入は 1e−10 r_t; plan
    [tooling-nozzle-mknot-fast](../plans/accepted/tooling-nozzle-mknot-fast.md) §4.4・§4.5・§6 V5)。壁のビット同一を要する照合は
    両方とも同じ backend で行う (`prepare_info.json` の `environment.moc_kernel`)。
    design の venv に `numba 0.68.0`・`llvmlite 0.50.0` を **`--no-deps`** で入れる (`.venv-opt/bin/pip install --no-deps numba==0.68.0
    llvmlite==0.50.0`; numpy を変えると δ_r が動く)。numba が無いと ImportError で止まる (黙って numpy の版に落ちない)。
    `cplus_flux_wall` と網の検査 (`characteristic_topology_diagnostics`) も numba を使う (出力は numpy の書き方と完全一致、ただし網の検査の
    記録の位置の一覧 `pileup_locs_interior` の選び方・並びは段・添字の順に変えた) ので、
    `FORGE_MOC_KERNEL=numpy` でも numba は要る。
    環境変数 `FORGE_MOC_KERNEL=numpy` で numpy の版 (A/B と回帰の基準)、他の値は ValueError。初回はコンパイルに数秒〜十数秒かかり、
    結果は `geometry/__pycache__` にキャッシュされる。試験: `design/tests/run_moc_jit_tests.py [--real]`。

## 0. 共通: problem YAML (契約書)

すべての出発点は `case/NN/problem_*.yaml`。読み込みと検査は
`design/forge_design/probdef.py` の `load_problem()` (未知の type・欠けた spec・組成の不整合で例外)。

| ブロック | 主なキー | 注意 |
| --- | --- | --- |
| `spec` | `Pt`, `Tt`, `p_ambient`, `r_throat` [m], `M_design`、任意で `sizing` | **出口径は直接書けない** (`D_e` は未実装)。出口半径は `r_throat` × MOC の面積比 + δ で決まる。物理出口径に合わせるには `solve_rt`、物理スロート径から決めるには `solve_rt_throat` (§1 ④)。どちらで決めたかを `sizing: {method: exit \| throat, target_m: …}` に書くと `prepare_info.json` の `sizing` に目標との差が出る |
| `gas` | `model: semiperfect`, `species` (組成), `transport` | Tt > 600 K・燃焼ガスは semiperfect。条件を変えたら `geometry.r_inlet` も付け替える |
| `dv` | `L_c.value` | 軸 law の長さ。runner は `min/max` を使わない |
| `geometry` | `R` (スロート曲率半径/r_t)、`r_inlet` (入口配管半径/**r_t 単位**)、`L_U` (縮流部長)、`L_pipe`、`axis_law: knot`、`M_knot`、`Lc_mode`、`start_line: throat_char`、`wall_mode: cplus`、`n_axis_inv` | `Lc_mode`: `explicit` (L_c を使う、許容窓の外は例外) / `max` / `from_length` (`dv.L_total` から L_c を解く)。r_t を変えても r_inlet・L_total は r_t 単位のまま (実寸がずれる) **要求カードでは M_knot は既定で自動** (最大壁角が最小・成立の境目から 0.05 以上、2026-10-14; plan `tooling-nozzle-mknot-auto`)。 |
| `geometry` (生産レシピ) | `initial_line: cfd` + `initial_line_run`/`initial_line_res`、`wall_repr: joint`、`Md_moc_offset`、`pw_upstream: poly` | §2 ⑤ 参照。新規設計の初回は `initial_line: hall`・`wall_repr: interp` (既定) でよい。`pw_ramp` は旧設定 (`pw_upstream: ramp` と組で使う。`poly` と併記すると例外) |
| `mesh` | `ni`, `nj`, `wall_first_frac` (+ `_throat`・ブレンド)、`throat_refine`、`throat_width`、`axis_cap_frac`、`ar_max` | **NS (`prepare_ns`) の格子**。壁解像 PASS の格子を使う (§2 ⑥)。Euler には効かない (2026-10-07 から) |
| `mesh_euler` | `mesh` と同じキー (未知のキーは例外) | **Euler (`prepare`) の格子** (2026-10-07、plan `verification-case45-euler-total-enthalpy` §4)。`mesh` からは補完しない。既定は全断面で `wall_first_frac` 0.005 の等比 (スロートの別指定なし・軸側の cap なし、ni × nj は 321 × 65)。**`mesh_euler` の無い問題は Euler の `prepare` が移行先を示して止まる** (旧格子の記録を再現するなら `mesh` の全キーを写す)。case/45 は ni 2000 × nj 97・軸方向は G1 と同じ (`throat_refine 4`・`throat_width 3`)・全域 0.005 |
| `evaluate` | `nStepOuter` (**Euler の本段だけ**)、`ns_steps` (NS の本段、既定 60,000)、`outStepInterval`, `cfl_main`, `tp_species`, `condensation` | 凝縮は `tp_species: split_h2o` + `condensation`。**NS の本段の step 数は `--steps` > `evaluate.ns_steps` > 既定 60,000** (2026-10-14 ユーザ決定。以前は雛形の `nStepOuter` 12,000 が NS にも効いて短すぎた)。`prepare_info.json` の `nStepOuter_source` に出所。**NS の本段は既定で停止規則つき・上限 20 万 step** (2026-10-14 ユーザ決定; `evaluate.stop_rule`・`stop_rule_cap`、記録は run の `stop_rule.json`、plan `tooling-nozzle-cooled-wall-ns-settling` §4.6) |
| `mesh`・`mesh_euler` の上限 | `max_nodes`, `max_nodes_reason` (組で書く) | **2D の節点数の上限 30 万点** (2026-10-14、plan `tooling-nozzle-mesh-size-guard`)。生成器・`prepare`・`prepare_ns` が格子を作る前に実効の節点数 (`axis_gap_frac` の導出 nj を含む) で判定して止まる。引き上げはユーザの了承を取ってから理由と組で書く (AGENTS.md「メッシュの規模は作る前に示す」)。`prepare_info.json` の mesh 欄に `n_nodes` |

新しい問題は既存の YAML を複製して作る。複製元のコメントに残る古い設定 (`M_design` の試作較正値、
`mesh.bndFirstOrder`、廃止キー) を持ち込まない ([`recommended-settings.md`](recommended-settings.md) §9)。
検査: `cd design && .venv-opt/bin/python -c "from forge_design.probdef import load_problem; load_problem('../case/NN/problem.yaml')"`。

## 1. 新規の非粘性ノズルを設計する (軸長・出口径・入口径・入口条件を指定)

1. **要件を決める** — skill `design-intake` (要件インタビュー → 検証済み YAML・方針メモ。本書は以降の実行手順を担う) で、Pt・Tt・組成・設計マッハ数 `M_design`・入口配管半径・軸長
   (または全長)・出口径を確定する。凝縮の可能性 (H₂O を含み膨張で低温) をここで確認する。
   **実寸で与えられた要求** (物理スロート径・入口半径・直管長・縮流部長・設計スロート→出口の長さ・壁温) は要求カード
   (`forge.wt_card/v1`) に書き、`-m forge_design.intake.wt_card CARD --out problem.yaml` で ② の YAML を作る
   (r_t を解いた後の長さの割り直しを生成器が行う。④ で r_t を解き直したら `--r-throat` で作り直す)。
2. **YAML を作る** — 近い既存ケースの `problem_*.yaml` を複製し、`spec`・`gas`・`geometry.r_inlet`
   (= 入口半径 / r_t)・`dv.L_c` (または `Lc_mode: from_length` + `dv.L_total`) を書き換える。
   `r_throat` は仮の値でよい (④ で出口径から解く)。
3. **MOC だけで壁を作る** — `design_chain(p)` (`design/forge_design/evaluate/runner_axismach.py`) が
   Hall のスロート初期線 → 軸 Mach law → 逆 MOC → 壁の QA を決定的に実行する (CFD なし)。
   ファイルに落とすには Euler run の準備だけを行う:
   ```
   .venv-opt/bin/python -m forge_design.evaluate.runner_axismach ../case/NN/problem.yaml ../case/NN/run_XXXX_moc --prepare-only
   ```
   → `wall_design.csv` (x, r, 壁角, 壁 M)・`target_axis_M.csv`・`nozzle.msh`/`nozzle.h5`・`MESH_QUALITY.txt`。
   設計が成立しないと例外になる (`law.gates()`・`wall_qa`: 半径の単調性・出口角 ≤ 0.2°・r_F の 1D 比 ±3 %・
   スロート曲率、`wall.validate()`: スプラインのうねり・C¹/C² 接続、縮流部 μ ≤ 20)。
   **fold・特性線トポロジ・壁マージン (μ_w − θ_w ≥ 1°) は `design_chain` に入っていない** —
   `design/forge_design/geometry/moc_diagnostics.py` で別に確認する (§2 ①)。
4. **寸法 (r_t) を決める** — 設計は r_t 単位なので、寸法は r_t の 1 変数。**出口径から**決めるか、**物理スロート径から**決めるかを選ぶ
   (同時に両方を固定条件にしない)。解いた r_t で `spec.r_throat` を手で書き換え、`spec.sizing` とコメントにどちらで決めたかを残す
   (plan `tooling-nozzle-upstream-poly-and-throat-sizing` §4.2):
   ```
   # 出口径から (物理出口半径 0.775 m)
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ../case/NN/problem.yaml --euler-ref X --run-dir Y --solve-rt 0.775
   # 物理スロート径から (物理スロート半径 [m])
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ../case/NN/problem.yaml --euler-ref X --run-dir Y --solve-rt-throat 0.07676
   ```
   ```yaml
   spec:
     r_throat: 0.0766539     # [m] solve_rt_throat (物理スロート半径 0.07676 m) で決めた (YYYY-MM-DD)
     sizing: {method: throat, target_m: 0.07676}
   ```
   (`--euler-ref`/`--run-dir` は必須引数だがこのモードでは使われない。) CFD 前の δ は `prepare_ns` と同じ経路 (YAML の
   `deltastar_initializer` の k_f・熱条件・平滑化、2026-10-07 から。以前の `solve_rt` は未較正・未平滑化で、生産の壁と出口半径が 2.4 mm
   ずれた) で、物理壁も `prepare_ns` と同じ構築 (`pw_upstream` に従う)。**k_f を較正し直したら寸法も解き直す** (k_f はスロートの δ も変える)。
   **CFD 前の寸法は初期見積もり**: 残差はスロート 1e-7 m・出口 1e-5 m (`SIZING_TOL_PRE_CFD_THROAT_M`・`SIZING_TOL_PRE_CFD_EXIT_M`、
   引数 `tol_R_m` で明示も可) で止め、30 回で収まらなければ不合格 (例外)。CFD 前の δ (積分法) は r_t に対して滑らかでなく
   (r_t を 1e-12 m 変えるだけで出口の δ が ~1e-5 r_t 揺れる)、局所の散らばりはスロート 7.4e-9 m・出口 2.6e-6 m (2026-10-07 U2c) で、
   精密な寸法・感度の評価には未検証 (plan `tooling-nozzle-upstream-poly-and-throat-sizing` §6 U2d・既知の制約)。NS 後の `solve_rt_throat` は
   `SIZING_TOL_M` (1e-9 m) で止め、NS 後の `solve_rt` (`--prev`) は固定回数の反復で許容差を判定しない。
   NS 後に物理スロート径から決め直すときは `--prev RUN --delta-next RUN/_extract_edge/delta_r_next.csv --delta-r-out FILE`
   (δ_E の全分布で寸法を決め、同じ補正の表を FILE に書く。次の壁は `prepare_ns(delta_r_csv=FILE)` でその表から作る)。
   **実寸の長さを保って作り直すとき (要求カードから; 2026-10-13、plan `tooling-nozzle-deltastar-band-root-zone`)**: r_t が動くと r_t 単位の設計を
   据え置く上の経路では実寸の長さが公差を外れる (case/70: r_t ×1.0041 で拡大部長 +4.4 mm)。カードから作り直す:
   ```
   .venv-opt/bin/python -m forge_design.intake.wt_card CARD --out P.yaml --ns-correction PREV_RUN --delta-next PREV_RUN/<抽出>/delta_r_next.csv
       → P.yaml・P_delta_r.csv (最終の r_t での補正の表)・P.yaml.card_log.json
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem P.yaml --euler-ref EULER_RUN --delta-r-csv P_delta_r.csv \
       --ic-from PREV_RUN --run-dir RUN [--stages none --cfl 1 --steps 60000]
   ```
   `--prev` / `--ns-correction` の前の run は `wall_repr.json` (前の設計の区切りの点) が必須。表の端は前の設計の入口端・出口に float32 の丸めの範囲で揃え、
   区間で正規化した座標 s で今の設計へ写す。使用域 [0, x_e] に hard 不合格 (held) の点があれば止まる。`--delta-r-csv` は表の sha256 をカードの記録と照合する。
   カードの寸法ゲートは**設計壁から測った値**で、最終の合格は成果物 (STEP・メッシュ) の出荷寸法ゲート。
   **冷却壁の注意**: case/70 (T_w/T_aw ≈ 0.3) の NS は、残差が床に着いた後も出口 M・流量が 10 万 step 単位で動き続けた (plan §5.1 #11)。
   作り直しの入力の δ_E と最終の出口 M は、長い窓の系列で動きが止まったか (出口 M は末尾の変化と推定の残りの変化 ≤ 2.5e-4) を確かめてから使う。
   **作り直しの入力の δ_E は、整定した場の 5 枚の窓 (5k ごと) で各枚の抽出 (帯は各枚で選ぶ、既定) の点ごとの平均にする** (2026-10-13、plan `tooling-nozzle-cooled-wall-ns-settling` §4.4)。1 枚だけでは、梯子が平らにならない区間 (冷却壁のスロート付近・直管・縮流部) で帯位置が切り替わり、局所で δ_E が 5e-4 r_t 動く。5 枚の平均でも 1 枚の切り替わりで 1e-4 r_t 程度は残る (既知の限界)。
5. **形を見る・出す** — 最終 NS の run から成果物一式を 1 コマンドで出す (plan
   [`tooling-nozzle-wt-deliverables`](../plans/accepted/tooling-nozzle-wt-deliverables.md)、2026-10-13):
   ```
   .venv-opt/bin/python -m forge_design.export.deliverables ../case/NN/run_XXXX --card ../case/NN/card_<tag>.yaml --problem ../case/NN/problem_<tag>.yaml
   ```
   → `run_XXXX/deliverables/`: STEP (2D の壁曲線・3D の内面)・STL・解析メッシュ (msh4.1/2.2)・STEP から同じ構造のメッシュを
   作り直す gmsh スクリプトと構造ファイル・Fluent 2D 軸対称メッシュ (`.msh`、境界の名前と型つき)・CGNS・Fluent の境界条件表
   (動作圧 0 Pa)・Excel 向けの形状表と壁面量表 (CSV、`xlsxwriter` があれば `.xlsx`)。README の冒頭に**寸法の出荷ゲート**
   (STEP とメッシュから測り直した寸法と要求カードの公差) と、壁面の熱量を熱設計に使えるか (`thermal_use`) が出る。
   物理壁の表現が `single_bspline` の run (最終 NS 用の雛形 `axismach_final_v1`) でないと STEP は出ない。
   壁の 2 階微分 (r″) の連続性は `nozzle_report` の「壁の形と微分」で見る (ユーザが重視する観点)。
   旧来の case 側のスクリプト (`case/42.isobutane_wt/export_wall_m42.py`、`case/44.vitiated_air_wt/export_wall_va.py`、
   `case/45.isobutane_m6_d155/cad_deliverables.py`・`viz_*.py`) は絶対パスを含むので新しい設計には使わない。

## 2. 出口・入口・条件を固定し、軸長を MOC でパラスタ → Euler → NS → 形状調整 → 凝縮

M6 (case/45) で実際に通した順番。各段の「何で判定するか」を先に決めてから回す (事前登録、AGENTS.md)。

1. **MOC パラスタ (CFD なし)** — 汎用ツールは無い。case のスクリプトを複製して使う:
   `case/45.isobutane_m6_d155/search_shortest.py` (R × L_c × M_knot の格子で、ハードゲート +
   `moc_diagnostics` の壁マージン ≥ 1°・特性線トポロジ ≥ 0.02 を満たす最短 x_F を探す)、
   `case/42.isobutane_wt/sweep_axislaw_A_MK.py`、`case/44.vitiated_air_wt/study_design_va.py`。
   定数 (条件・パス) がスクリプト内に書かれているので、問題に合わせて書き換える。
   有望形状の基準: 全長 x_F、壁マージン、特性線の健全さ、壁の r″ のなめらかさ。
2. **Euler (slip) で候補を評価** — 候補ごとに problem YAML を作り
   ```
   .venv-opt/bin/python -m forge_design.evaluate.runner_axismach ../case/NN/problem_X.yaml ../case/NN/run_XXXX_euler_X [--stages soft --cfl 2 --implicit-relax 0.7]
   ```
   → `metrics.json` (`collect()`: 軸 M の偏差・出口 M・オーバーシュート・出口一様性・流量比)。
   形状比較は `nozzle_report` の評価量 (r/r_w = 0.1 のMach 波・オーバーシュート・出口コア M) で行う。
   複数候補は `case/44.../run_study_va.py` のようにバッチで回す。
3. **NS + 境界層補正 (δ* ループ)** — `design/forge_design/feedback/deltastar_loop.py`:
   ```
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ../case/NN/problem_ns.yaml --euler-ref ../case/NN/run_EULER --init-integral --ic-from ../case/NN/run_EULER --run-dir ../case/NN/run_XXXX_ns_pass0
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ... --euler-ref ... --prev ../case/NN/run_XXXX_ns_pass0 --run-dir ../case/NN/run_YYYY_ns_pass1 [--stages none --cfl 5 --implicit-relax 0.7 --steps 12000]
   ```
   物理壁 = 設計壁 + δ_r。δ は CONTUR 積分法で初期化し、NS から E 法 (`band_select="edge"`) で δ_E を測る。
   **注意 (2026-10-13)**: `--ic-from <Euler run>` は、`interp_field.py` が化学種の照合で拒否する (Euler は輸送の記録が無く、NS は
   `transport` を持つので species_hash が違う)。新しい設計の最初の NS は `--ic-from` を付けず、prepare_ns の 1 次元等エントロピー場から
   段階起動 (`--stages full`) で始める (case/70 run_0005)。2 回目以降の NS は前の NS の場から (`--prev`)。
4. **C2 方式 (出口 δ で CONTUR を較正し r_t を解く)** — CONTUR の摩擦係数倍率 k_f を「出口での δ_C = 測った δ_E」
   に合わせ、出口径 0.775 m (例) から r_t を解き直す。**k_f は YAML キーが無く case スクリプト経由**:
   `case/45.isobutane_m6_d155/prep_c2pin.py` (run 準備、k_f を引数で渡す)、`c2pin_solve.py PASS1_RUN EULER_REF
   [--base YAML] [--out NAME] [--cond-steps N] [--cond-out N]` (未緩和 δ_E から k_f・r_t を一発で解き
   `<NAME>.yaml`・`<NAME>_cond.yaml` を書く)、バッチ例 `run_recal_chain.sh`。
   δ は **`delta_r_next.csv` の `delta_E` 列 (未緩和)** を使う (第 2 列 `delta_target` は ω 緩和後 — 取り違えた前例あり)。
5. **形状を整える (生産レシピ、2026-10)** — 探索段階では不要、形状を詰める段階で入れる:
   - `wall_repr: joint` — MOC 点を位置 + 壁角で同時に当てはめた 5 次 B-spline (壁の r″ が連続でなめらか)。
   - `initial_line: cfd` — スロート初期線と軸アンカーを Hall 解でなく Euler の遷音速場から取る
     (`initial_line_run`/`initial_line_res` で凍結源を指定; 凍結源は Hall 初期線の当てはめ壁を Euler で解いた場)。
   - `Md_moc_offset` — Euler の出口コア M を 6 に合わせる 1 係数の較正。**Euler は `mesh_euler` の壁に寄せない配点で決める**
     (2026-10-07 から。case/45 は 2000 × 97・軸方向は NS の G1 と同じ・全断面で壁の第 1 間隔の比 0.005)。NS と同じ壁に寄せた
     G1 の配点 (1.3e-5・スロート 4.5e-6) の Euler は、スロート付近の全温が Tt を数百 K 超えたまま整定しなかった
     (plan `verification-case45-euler-total-enthalpy` §9 E2。全域 0.005 の配点は同じ窓で全領域 |T₀ − Tt| ≤ 0.103 K)。出口コア M は最終断面の
     η ∈ [0.05, 0.7] で、NS の基準格子 G1 の固定の η の列に線形補間した平均 (M_common) で較正し、自格子の平均も残す (同 plan §6 E4)。
     Euler の較正の合格を NS に移せるとは限らない (配点も粘性も違う) ので、**NS の出口 M は NS で判定する**。
     旧手順「生産 NS と同じ格子パラメータの Euler で決める」(plan `tooling-nozzle-cfd-pinned-initial-line` §5.1 #11f。粗い Euler 格子
     1100 × 65 で決めると細分格子の NS で出口 M が約 −0.0008 ずれた経緯) は、G1 の Euler の全温が整定しない (E2) ため取り下げた。
     case/45 の旧値 +3.770e-4 は G1 の Euler (run_0113+0114) 由来で、新しい配点でやり直した (同 plan §6 E4・E4V)。現行値は 6.8825162455159465e-06 (run_0164)。
   - `pw_upstream: poly` (既定。新しい YAML には明記する) — スロートより上流に δ_r を足さず、配管〜設計スロートを 5 次多項式 1 本にする
     (plan `tooling-nozzle-upstream-poly-and-throat-sizing`、2026-10-07)。`pw_ramp` は書かない (併記すると例外)。
     旧設定: `pw_upstream: ramp` + `pw_ramp: [−11, −6]` (縮流部側で δ_r をなめらかに入れる区間)。case/45 の生産 YAML は NS・凝縮の再評価
     (同 plan §6 U4) で採否を決めるまでこの旧設定を明示している。
   - `wall_fit_mono_r2: [0.0, 1.5]` — joint 壁のスロート直後の r″ の山を消す単調拘束 (case/45 の生産、2026-10-07)。
     Euler は実務判定 (全量で 差 + 2SE ≤ 許容幅) に合格、dry NS は登録したゲートに合格、凝縮 NS は 4 量の準定常に合格。新旧の場の差は参考値 (plan `tooling-nozzle-throat-monotone-r2`)。
   - `deltastar_initializer: {model: contur, a_crocco: 1.0, cf_scale: <k_f>, n_scale: 1.0}` — C2 で較正した k_f を YAML に書く。
     `deltastar_loop --init-integral` も読む (書かないと k_f = 1 の物理壁になる)。最終 NS の数値レシピ (CFL 1・60000 step など) は
     YAML ではなく投入スクリプトが持つ (case/45 の手本 `run_mono_ns_chain.sh`、延長 `run_mono_ns_ext.sh`)。
6. **NS の格子** — 壁解像 (y1+ ≤ 1、`check_wall_resolution.py --over-frac 5`) を満たす細分格子を使い、
   格子ゲート (細分前後で δ_E(x_F) の差 ≤ 1 %、三水準目で確認) を通す。case/45 の生産格子 (NS の `mesh`、G1):
   ni 2000 × nj 97、第 1 セル 1.3e-5 (スロート 4.5e-6)、`throat_refine 4`・`throat_width 3`、AR ≤ 5000 (境界層層の例外)。
   この壁に寄せた配点は NS 用で、Euler には使わない (Euler は `mesh_euler`、§0)。
   **細分格子の NS は本段 CFL 5 で発散する** — 段階起動 (soft → mid → 本段) + 本段 cfl 1・60000 step
   (累積 CFL を cfl 5 × 12000 にそろえる)。
7. **最終判定** — 事前登録したゲート (case/45 の例: 出口半径 0.775 m ± 0.1 mm、δ_E/δ_C 1 ± 0.5 %、出口コア M 6.000 ± 0.05 % [2026-10-07 ユーザ決定、旧 ± 0.02 %]、
   r/r_w = 0.1 の波 ≤ 0.01 %・オーバーシュート ≤ 0.035 %、壁解像 PASS) を、各量の時系列に
   `check_quasisteady.py --series-csv` をかけて判定する (`nozzle_report` の末尾幅は VERDICT の代わりにならない)。
   時系列は `case/45.isobutane_m6_d155/exitM_sampling_ab.py` (環境変数 `EULER_REF`・`SOLVE_JSON`・`FINAL_PROBLEM`)。
8. **凝縮 ON** — 最終 NS の場を `convert_species_field.py --mode conserve` で液相付きの種構成に変換して引き継ぐ
   (restart_field は種構成が違うと拒否する):
   ```
   python3 solver_density_cuda/tools/convert_species_field.py LAST_res.h5 RUN_COND/nozzle.h5 --meta RUN_COND/species_meta.yaml --src-run RUN_NS --dst-run RUN_COND
   ```
   problem は `<NAME>_cond.yaml` (`tp_species: split_h2o` + `condensation`)、cfl 1。凝縮 4 量
   (始まり位置・S_max・出口 g・出口コア M) が末尾 5 枚で STEADY になるまで回す。
9. **報告** — `nozzle_report` → pptx ([`nozzle-design-outputs.md`](nozzle-design-outputs.md))。
   case README の run 一覧と plan §9 に記録する。

## 3. 作成済みのノズルで Euler / NS を回す

**まず「形状を固定するか」を決める。** problem YAML の `spec` (Pt・Tt)・`gas` を変えて `runner_axismach` や
`deltastar_loop` に渡すと、`design_chain` が MOC と δ 補正をやり直し、**壁そのものが変わる** (= 作り直し、§1・§2)。
既にあるノズルで条件だけ変えたいときは、YAML を通さず run ディレクトリの config を直接書き換える (3a)。

### 3a. 形状は固定、入口条件・背圧・入口分率だけ変えて回す (よくある用途)

ツール **`solver_density_cuda/tools/rerun_conditions.py`** を使う (plan [`tooling-rerun-conditions.md`](../plans/accepted/tooling-rerun-conditions.md)、仕様は `methods/design/overview.md`「既存 run の条件変更」)。
既存 run の `nozzle.h5` (物理壁の格子) を再利用し、bcond の floats だけを書き換え、参照 run の最終場を初期場に入れ、変更点を `RERUN_CONDITIONS.json` に残す。**forge は起動しない。**

```
python3 solver_density_cuda/tools/rerun_conditions.py REF_RUN NEW_RUN [--Pt P --Ps P | --keep-Ps] [--Tt T] [--Y H2O=0.09] \
    [--Tw T | --keep-Tw] [--scale-ic pt] [--steps N --out-interval N --cfl C] [--dry-run]
```

- 対応する入力: 単一の `inlet_Pressure` (Y{s} 形式、inletProfile なし)・`outlet_statPress`・`wall`/`wall_isothermal`/`slip`・`axis`、`nozzle.h5` を mesh/value に使う run。外れれば作る前に止まる。
  NS / Euler は実効 config で判定する (粘着壁があれば NS; 全壁 `slip` は粘性・熱伝導・種拡散・乱流が無効な config だけ Euler として受理し、輸送が有効なら止まる)。重複した YAML キーを含む config は止まる (ツールとソルバで解釈が違うため)。
- **Pt を変えるときは背圧の指定 (`--Ps` か `--keep-Ps`) が必須**。node の出口は壁際の列が常に背圧を見るので、Pt だけ下げると出口の壁際から崩れる。
  ツールは参照の出口断面の静圧に対する `Ps/(f·P_exit_ref)` を表示する (相似にするなら Ps も同じ比で)。
- 乾き成分 (lump) の組成変更は v1 では止まる (化学種の定義が変わり、どの引き継ぎ経路も拒否する)。H2O 分率 (`--Y H2O=...`) は可。
- 入口の k・ω は変えない (指定時のみ)。等温壁で Tt を変えるときは `--Tw` か `--keep-Tw` を明示する。
- **回し方** (`RERUN_CONDITIONS.json` の `recommended_stages` に従う): 条件の変更なし → `run_staged_ns(run, stages="none")`; 条件を変えた → `stages="full"` (段階起動);
  **Pt だけを変えた → full + 本段 cfl 1 (`--cfl 1.0 --steps 60000`)、`--scale-ic pt` (保存量を Pt 比で一様スケールする初期場変換) を推奨** (2026-10-06 確定: scale なしは規定時間内に準定常に達せず、入口配管の壁際に逆流域が残った)。**Pt と Tt・組成・壁温を同時に変える、または凝縮 run で Pt を変えるときは `--scale-ic pt` が使えないか推奨外で、起動・整定は未検証** (full・本段 cfl 1・60000 step と警告が記録される)。
  粗格子の Pt 0.8 倍は、段階起動でも本段 cfl 5 では発散した (plan §6 (ii′))。Tt・H2O を変えたときは本段の量の整定に時間がかかる (粗格子で本段 36000 step でも流量・δ_E は DRIFTING; 整定長は未確定、plan §5.1 #11)。
- **Euler 参照**: δ_E の抽出や流量比には、同じ条件にした Euler の rerun を対で作る (旧条件の Euler 参照だと δ が数 % 動く)。
- 判定: `check_convergence.py` (`--segment` で本段だけ)・`check_quasisteady.py --series-csv` (報告する量)・NS は `check_wall_resolution.py` (Pt・Tt が変わると Re と y1+ が変わる)。
- 注意: 壁は元の条件で設計したもの。条件を変えると境界層の厚さが変わり、出口 M・試験部の一様性は設計値からずれる (それを見るのがこの計算の目的のことが多い)。ずれを壁に返したいなら §2 の作り直し。

### 3b. 既存の problem YAML から同じ設計を回し直す (作り直しを含む)

- Euler: `.venv-opt/bin/python -m forge_design.evaluate.runner_axismach PROBLEM RUN_DIR [--stages full|soft|none] [--cfl C] [--implicit-relax 0.7] [--steps N] [--ic-from RUN]`
- NS: runner に NS フラグは無い。`deltastar_loop` (`--init-integral` で CONTUR の δ から、`--prev RUN` で前の NS の δ_E から)、
  または C2 で解いた k_f で `prep_c2pin.py RUN K_F --problem YAML --ic RUN --stages full` → `run_staged_ns`
  (`case/45.isobutane_m6_d155/run_recal_chain.sh` が手本)。
- YAML の Pt・Tt・ガスを変えると壁が変わる (上記)。同じ壁を再現したいだけなら YAML を変えない。
- **コードの版も壁を変える**: 2026-10-13 の出口の終点の修正 (de35cb0b) から、終端の壁点を丸めで落としていた設計は同じ YAML でも壁・メッシュ・STEP が 1 点間隔長くなる (case/70 c2b +0.846 mm、case/45 生産 +5.353 mm)。古い run と同じ壁を再現するには、その run の `wall_repr.json`・`wall_physical.csv` (保存された壁) と `RUN_PROVENANCE.txt` のコードの版を確かめる (YAML だけでは再現しない)。

### 3c. 外部の壁座標 (CSV・CAD) しか無い場合

axismach パイプラインに取り込むキーは**未実装** (`wall_csv` 相当なし)。Gmsh `.geo` → `gmsh` → 変換器 → config を手で組む
一般経路 ([`calculation-workflow.md`](calculation-workflow.md)) を使い、config は skill `forge-config` で組む。

### 共通

- 同じ格子の続きは `restart_field.py` (ビット一致を検査)、格子が違えば `interp_field.py`。一様初期場からは始めない。
- 判定: `check_mesh_quality.py` (投入前)、`check_convergence.py` (同一設定の区間)、`check_quasisteady.py` (報告する量そのもの)、
  NS は `check_wall_resolution.py`。残差が plateau (`NOT CONVERGED`) になるのはこの系の常態 — そのまま記録し、量の準定常で判断する。

## 4. 未実装・注意 (2026-10-06 時点)

- MOC パラスタ・作図・凝縮の引き継ぎは **case スクリプト頼み** (汎用 CLI なし)。壁の出力は `forge_design.export.deliverables` (§1 ⑤)、
  C2 の k_f は `deltastar_initializer.cf_scale`。
  統一は campaign-recipe plan の残作業 (#6〜#8・#14・#15)。
- 出口径・スロート径を spec に書けない (r_throat + M_design → `solve_rt` / `solve_rt_throat` → 手で書き換え。`spec.sizing` は記録だけ)。
- `design_chain` には fold・トポロジ・壁マージンの検査が入っていない。
- AWS の共有バイナリは他セッションが作り直すことがある — 自分の worktree でビルドしたバイナリを使い、
  `RUN_PROVENANCE.txt` の sha256 を記録する (case/45 は `~/forge-wallfit-bin`)。
