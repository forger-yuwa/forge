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
- コマンドは `design/` 直下で `.venv-opt/bin/python` (または `python3`) を使う。run は検証 run でも AWS で回す。

## 0. 共通: problem YAML (契約書)

すべての出発点は `case/NN/problem_*.yaml`。読み込みと検査は
`design/forge_design/probdef.py` の `load_problem()` (未知の type・欠けた spec・組成の不整合で例外)。

| ブロック | 主なキー | 注意 |
| --- | --- | --- |
| `spec` | `Pt`, `Tt`, `p_ambient`, `r_throat` [m], `M_design` | **出口径は直接書けない** (`D_e` は未実装)。出口半径は `r_throat` × MOC の面積比 + δ で決まる。物理出口径に合わせるには `solve_rt` (§2 ④) |
| `gas` | `model: semiperfect`, `species` (組成), `transport` | Tt > 600 K・燃焼ガスは semiperfect。条件を変えたら `geometry.r_inlet` も付け替える |
| `dv` | `L_c.value` | 軸 law の長さ。runner は `min/max` を使わない |
| `geometry` | `R` (スロート曲率半径/r_t)、`r_inlet` (入口配管半径/**r_t 単位**)、`L_U` (縮流部長)、`L_pipe`、`axis_law: knot`、`M_knot`、`Lc_mode`、`start_line: throat_char`、`wall_mode: cplus`、`n_axis_inv` | `Lc_mode`: `explicit` (L_c を使う、許容窓の外は例外) / `max` / `from_length` (`dv.L_total` から L_c を解く)。r_t を変えても r_inlet・L_total は r_t 単位のまま (実寸がずれる) |
| `geometry` (生産レシピ) | `initial_line: cfd` + `initial_line_run`/`initial_line_res`、`wall_repr: joint`、`Md_moc_offset`、`pw_ramp` | §2 ⑤ 参照。新規設計の初回は `initial_line: hall`・`wall_repr: interp` (既定) でよい |
| `mesh` | `ni`, `nj`, `wall_first_frac` (+ `_throat`・ブレンド)、`throat_refine`、`throat_width`、`ar_max` | Euler と NS で同じキーが効く (2026-10-06 修正)。NS は壁解像 PASS の格子を使う (§2 ⑥) |
| `evaluate` | `nStepOuter`, `outStepInterval`, `cfl_main`, `tp_species`, `condensation` | 凝縮は `tp_species: split_h2o` + `condensation` |

新しい問題は既存の YAML を複製して作る。複製元のコメントに残る古い設定 (`M_design` の試作較正値、
`mesh.bndFirstOrder`、廃止キー) を持ち込まない ([`recommended-settings.md`](recommended-settings.md) §9)。
検査: `cd design && .venv-opt/bin/python -c "from forge_design.probdef import load_problem; load_problem('../case/NN/problem.yaml')"`。

## 1. 新規の非粘性ノズルを設計する (軸長・出口径・入口径・入口条件を指定)

1. **要件を決める** — skill `design-intake` (要件インタビュー → 検証済み YAML・方針メモ。本書は以降の実行手順を担う) で、Pt・Tt・組成・設計マッハ数 `M_design`・入口配管半径・軸長
   (または全長)・出口径を確定する。凝縮の可能性 (H₂O を含み膨張で低温) をここで確認する。
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
4. **出口径に合わせる** — 非粘性の出口半径は r_t × r_F。物理出口径 (境界層込み) に合わせるときは
   `solve_rt` で r_t を解き、`spec.r_throat` を手で書き換える:
   ```
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ../case/NN/problem.yaml --euler-ref X --run-dir Y --solve-rt 0.775
   ```
   (`--euler-ref`/`--run-dir` は必須引数だがこのモードでは使われない。δ は CONTUR 積分法の見積もり。)
5. **形を見る・出す** — 作図・CAD 用の全輪郭出力は case 側のスクリプトしかない
   (`case/42.isobutane_wt/export_wall_m42.py`、`case/44.vitiated_air_wt/export_wall_va.py`、
   `case/45.isobutane_m6_d155/viz_*.py`)。**中に `/home/sano/work/forge/...` の絶対パスが書かれているものがある**
   ので、使う前に確認する。壁の 2 階微分 (r″) の連続性も見る (ユーザが重視する観点)。

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
   - `Md_moc_offset` — Euler の出口コア M を 6 に合わせる 1 係数の較正。**生産 NS と同じ格子パラメータの Euler で決める**
     (粗い Euler 格子で決めると細分格子の NS で出口 M が約 −0.0008 ずれる; plan §5.1 #11f)。
   - `pw_ramp` — 縮流部側で δ_r をなめらかに入れる区間 (case/45 は [−11, −6])。
6. **NS の格子** — 壁解像 (y1+ ≤ 1、`check_wall_resolution.py --over-frac 5`) を満たす細分格子を使い、
   格子ゲート (細分前後で δ_E(x_F) の差 ≤ 1 %、三水準目で確認) を通す。case/45 の生産格子:
   ni 2000 × nj 97、第 1 セル 1.3e-5 (スロート 4.5e-6)、`throat_refine 4`・`throat_width 3`、AR ≤ 5000 (境界層層の例外)。
   **細分格子の NS は本段 CFL 5 で発散する** — 段階起動 (soft → mid → 本段) + 本段 cfl 1・60000 step
   (累積 CFL を cfl 5 × 12000 にそろえる)。
7. **最終判定** — 事前登録したゲート (case/45 の例: 出口半径 ±0.1 mm、δ_E/δ_C 1 ± 0.5 %、出口コア M 6.000 ± 0.02 %、
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

ツール **`solver_density_cuda/tools/rerun_conditions.py`** を使う (plan [`tooling-rerun-conditions.md`](../plans/active/tooling-rerun-conditions.md)、仕様は `methods/design/overview.md`「既存 run の条件変更」)。
既存 run の `nozzle.h5` (物理壁の格子) を再利用し、bcond の floats だけを書き換え、参照 run の最終場を初期場に入れ、変更点を `RERUN_CONDITIONS.json` に残す。**forge は起動しない。**

```
python3 solver_density_cuda/tools/rerun_conditions.py REF_RUN NEW_RUN [--Pt P --Ps P | --keep-Ps] [--Tt T] [--Y H2O=0.09] \
    [--Tw T | --keep-Tw] [--scale-ic pt] [--steps N --out-interval N --cfl C] [--dry-run]
```

- 対応する入力: 単一の `inlet_Pressure` (Y{s} 形式、inletProfile なし)・`outlet_statPress`・`wall`/`wall_isothermal`/`slip`・`axis`、`nozzle.h5` を mesh/value に使う run。外れれば作る前に止まる。
- **Pt を変えるときは背圧の指定 (`--Ps` か `--keep-Ps`) が必須**。node の出口は壁際の列が常に背圧を見るので、Pt だけ下げると出口の壁際から崩れる。
  ツールは参照の出口断面の静圧に対する `Ps/(f·P_exit_ref)` を表示する (相似にするなら Ps も同じ比で)。
- 乾き成分 (lump) の組成変更は v1 では止まる (化学種の定義が変わり、どの引き継ぎ経路も拒否する)。H2O 分率 (`--Y H2O=...`) は可。
- 入口の k・ω は変えない (指定時のみ)。等温壁で Tt を変えるときは `--Tw` か `--keep-Tw` を明示する。
- **回し方** (`RERUN_CONDITIONS.json` の `recommended_stages` に従う): 条件の変更なし → `run_staged_ns(run, stages="none")`; 条件を変えた → `stages="full"` (段階起動);
  **Pt を変えた → full + 本段 cfl 1 (`--cfl 1.0 --steps 60000`)、`--scale-ic pt` (保存量を Pt 比で一様スケールする初期場変換) を推奨** (2026-10-06 確定: scale なしは規定時間内に準定常に達せず、入口配管の壁際に逆流域が残った)。
  粗格子の Pt 0.8 倍は、段階起動でも本段 cfl 5 では発散した (plan §6 (ii′))。Tt・H2O を変えたときは本段の量の整定に時間がかかる (確認中、plan §6 (iv″))。
- **Euler 参照**: δ_E の抽出や流量比には、同じ条件にした Euler の rerun を対で作る (旧条件の Euler 参照だと δ が数 % 動く)。
- 判定: `check_convergence.py` (`--segment` で本段だけ)・`check_quasisteady.py --series-csv` (報告する量)・NS は `check_wall_resolution.py` (Pt・Tt が変わると Re と y1+ が変わる)。
- 注意: 壁は元の条件で設計したもの。条件を変えると境界層の厚さが変わり、出口 M・試験部の一様性は設計値からずれる (それを見るのがこの計算の目的のことが多い)。ずれを壁に返したいなら §2 の作り直し。

### 3b. 既存の problem YAML から同じ設計を回し直す (作り直しを含む)

- Euler: `.venv-opt/bin/python -m forge_design.evaluate.runner_axismach PROBLEM RUN_DIR [--stages full|soft|none] [--cfl C] [--implicit-relax 0.7] [--steps N] [--ic-from RUN]`
- NS: runner に NS フラグは無い。`deltastar_loop` (`--init-integral` で CONTUR の δ から、`--prev RUN` で前の NS の δ_E から)、
  または C2 で解いた k_f で `prep_c2pin.py RUN K_F --problem YAML --ic RUN --stages full` → `run_staged_ns`
  (`case/45.isobutane_m6_d155/run_recal_chain.sh` が手本)。
- YAML の Pt・Tt・ガスを変えると壁が変わる (上記)。同じ壁を再現したいだけなら YAML を変えない。

### 3c. 外部の壁座標 (CSV・CAD) しか無い場合

axismach パイプラインに取り込むキーは**未実装** (`wall_csv` 相当なし)。Gmsh `.geo` → `gmsh` → 変換器 → config を手で組む
一般経路 ([`calculation-workflow.md`](calculation-workflow.md)) を使い、config は skill `forge-config` で組む。

### 共通

- 同じ格子の続きは `restart_field.py` (ビット一致を検査)、格子が違えば `interp_field.py`。一様初期場からは始めない。
- 判定: `check_mesh_quality.py` (投入前)、`check_convergence.py` (同一設定の区間)、`check_quasisteady.py` (報告する量そのもの)、
  NS は `check_wall_resolution.py`。残差が plateau (`NOT CONVERGED`) になるのはこの系の常態 — そのまま記録し、量の準定常で判断する。

## 4. 未実装・注意 (2026-10-06 時点)

- MOC パラスタ・壁出力・作図・C2 の k_f・凝縮の引き継ぎは **case スクリプト頼み** (汎用 CLI・YAML キーなし)。
  統一は campaign-recipe plan の残作業 (#6〜#8・#14・#15)。
- 出口径を spec に書けない (r_throat + M_design → `solve_rt` → 手で書き換え)。
- `design_chain` には fold・トポロジ・壁マージンの検査が入っていない。
- AWS の共有バイナリは他セッションが作り直すことがある — 自分の worktree でビルドしたバイナリを使い、
  `RUN_PROVENANCE.txt` の sha256 を記録する (case/45 は `~/forge-wallfit-bin`)。
