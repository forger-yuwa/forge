# case/66 — ホストメモリ削減の回帰ハーネス

## ケース概要

plan [architecture-solver-host-memory](../../plans/active/architecture-solver-host-memory.md) §5.1 #3 (基準入力の確定と回帰ハーネス) の実体。
§6 の構成表の各構成を、変更前ビルド (base) と変更後ビルド (new) で**同じ入力・同じ手順で 3 回ずつ**回し、
§6 の (a) step 0 のビット一致・(b) N step 後の差が同ビルド内ばらつきの 2 倍以内・(c) NaN/Inf 無し・出力互換性 (データセット集合・shape・dtype・属性) を
`compare_runs.py` で判定する。

- **base** = `9c9f623c` (AWS `~/bin-hostmem/forge_9c9f623c`、sha256 `c65dbb39f772…`、変換器 `~/bin-hostmem/convert_9c9f623c`)。
- **new** = `93e55957` (R1–R3。AWS `~/bin-hostmem/forge_93e55957` sha256 `7121a8e51e7b…`、変換器 `~/bin-hostmem/convert_93e55957` `144db0f8b20f…`。
  `~/forge-b4` を 93e55957 へ `git checkout -f` し、新しいビルドディレクトリ `solver_density_cuda/build_93e55957` でクリーンビルド。
  base と同じ `RelWithDebInfo`・`CMAKE_CUDA_ARCHITECTURES=86`・`CMAKE_CXX_FLAGS=-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl`)。
  ハーネスはどちらも `.bin/<ビルド名>/forge` (シンボリックリンク) 経由で起動する。
- **run の置き場所は AWS** `~/forge-b4/case/66.hostmem_regression/run_NNNN_<構成>_<ビルド>_r<k>/`。比較も AWS で行い、手元にはテキストだけを持ち帰る
  (ローカルでは forge を回さない。重い h5 を手元へ引かない)。このディレクトリ (ローカル) には入力 config・スクリプト・README だけを置く
  (h5 は `.gitignore` で追跡外)。
- 元入力は別セッションのワークツリー `/home/sano/work/forge/case/...` にあり、**読むだけ**で複製して使う (`prepare_inputs.py`)。
- SERN 3D (g3) は `case/46.sern_design` の base `run_1072`〜`run_1074`・new `run_1075`〜`run_1077` (AWS `~/forge-r8/case/46.sern_design/`、[case/46 README](../46.sern_design/README.md) の run 一覧)。

## ファイル

| ファイル | 役割 |
| --- | --- |
| `matrix_spec.py` | **構成表の正本** (入力の元 run・複製するファイル・設定の修正・種・checkpoint・構成ごとの環境変数) |
| `prepare_inputs.py` | (ローカル) 元 run から `inputs/<入力名>/` を作る。修正は正規表現 + 期待一致数で行い、PyYAML で期待値を検査。`mesh.bndFirstOrder`・`wallTreatmentSST: 1` が残れば失敗 |
| `run_matrix.py` | (AWS) `seed` (restart_field.py を種に掛ける)・`set-ckpt`・`launch` (run を作って 1 本ずつ順に回すワーカーを裏で起動)・`resume`・`status`・`verify` (RUN_PROVENANCE の forge_bin/sha256 照合)・`note` (README の状態・比較から除外)・`table` (下の run 一覧の行) |
| `compare_runs.py` | (AWS) §6 の判定 (登録判定 A、既定 `--metric m`)。`--cfg X` / `--all` (registry.tsv から run を選ぶ)、`--base … --new …` (直接指定)、`--diff2 A.h5 B.h5` (分割と連続など 2 ファイル)。`--metric abs` で追加診断 B (plan §6.2) |
| `fixedwidth_eval.py` | 固定幅の独立 A/B (plan §6.3): `freeze` (既存 base 3 本から T = 2·S0 を凍結) と `eval` (段階 2 の 12 本を凍結した T で評価、§6.3 の判定表) |
| `fixedwidth_c44dual_ckpt100/` | §6.3 の段階 1 の成果物 (`T_frozen.tsv`・`T_frozen.sha256`・`PLAN.txt`・`plan.json`・`notes.txt`)。段階 2 は `run_matrix.py launch-plan fixedwidth_c44dual_ckpt100/plan.json` |
| `test_compare_abs.py` | 追加診断 B の判定関数 `judge_abs`・自己検査 `self_check_abs` の単体試験 (codex の最小再現 2 つを含む)。`python3 test_compare_abs.py` |
| `memlog_summary.py` | (AWS) `FORGE_MEMLOG=1` の工程別 RSS/HWM/GPU と `--memwatch` の 1 s 採取 (`mem_samples.csv`) を表にする |
| `split_vs_cont.py` | (AWS) dual-time の分割 (ckpt100 → 再開 100) と連続 200 の差を、ビルドごとに反復内の差と並べる (判定はしない) |
| `results/<日付>_<base>_vs_<new>_abs/` | 追加診断 B のテキスト (構成ごと・`all_quantities.tsv`・`summary.txt`・`selfcheck.txt`・`test_compare_abs.txt`・`sern_g3/`) |
| `results/<日付>_<base>_vs_<new>/` | AWS から持ち帰った判定のテキスト (`summary.txt`・構成ごとの報告・`sern_g3.txt`・`sern_g3_memlog.txt`・`split_vs_cont_c09.txt`・`registry.tsv`) |
| `inputs/<入力名>/` | 入力テンプレート。`SOURCE.txt` (元・sha256・修正の全記録)、`FILES` (run に写すファイル)、AWS では `SEEDED.txt` (種の適用記録)・`CKPT_FROM.txt` |

## 構成表

N は入力の `nStepOuter`。出力は `outStepInterval = N` で初期出力 `res_0` と最終 `res_N` の 2 回 (c09cont200 だけ 100 ごと)。2D/小格子だけなので N は全部 100〜200。

| 構成 | 入力 (`inputs/`) ← 元 run (`/home/sano/work/forge/case/` 相対) | 直した点 | N | 環境変数 | 確認する成果物 |
| --- | --- | --- | --- | --- | --- |
| `c36node` | `c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 旧キー LESorRANS/LESmodel/RANSmodel → `model: sst`、`wallTreatmentSST` 1→0、probe 3 点を追加 (20 step ごと) | 200 | - | `res_*.h5` (集合・shape・dtype・属性)・残差・壁 h5 `res_wall_4_*`・`point_probe_*.out` |
| `c36node_impdiag` | `c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 旧キー LESorRANS/LESmodel/RANSmodel → `model: sst`、`wallTreatmentSST` 1→0、probe 3 点を追加 (20 step ごと) | 200 | `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `diag.csv` |
| `c36node_psidual` | `c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 旧キー LESorRANS/LESmodel/RANSmodel → `model: sst`、`wallTreatmentSST` 1→0、probe 3 点を追加 (20 step ごと) | 200 | `FORGE_DIAG_PSI_DUALEVAL=2,4` | `psi_dualeval.csv`・log `snapshot: X of Y device arrays` (base は 269 of 269) |
| `c09ckpt100` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | - | `/CHECKPOINT` 14 データセットと属性 (layout・nHistoryValid・totalTime) |
| `c09restart100` | `c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | 無修正。`ic_ckpt.h5` = base r1 の c09ckpt100 の res_100.h5 (全ビルド・全反復で同じファイル) | 100 | - | log の履歴復元行・`/CHECKPOINT`・分割 (本構成の res_100) と連続 (c09cont200 の res_200) の差 |
| `c09cont200` | `c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | 無修正 (連続 200 step、res_100/res_200) | 200 | - | res_100・res_200 |
| `c09ckpt100_outres` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | `FORGE_OUT_RESIDUALS=1`, `FORGE_RESID_SNAP=0` | `res_*`・`dq_*`・`res_*_m`・`dq_*_new` の有無 |
| `c09ckpt100_rawdiag` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | `FORGE_SPECIES_RAW_DIAG=1` | `roYraw0/1` の有無 |
| `c09ckpt100_pindiag` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | `FORGE_PIN_DIAG=1` | log のピン診断行 |
| `c44steady` | `c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 無修正 (定常 active run、res_24000 から restart)。凝縮モーメント rog_0/roQ*_0 は変換直後の nozzle.h5 に無いので 0 で作ってから写す。種の記録が無い場なので restart_field は --force-species、solver は FORGE_ALLOW_UNVERIFIED_SPECIES=1。species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須、現行 solver が拒否する) | 200 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | 種・モーメント (`rog_0`・`roQ*_0`)、`res_outlet_2_*` |
| `c44dual_ckpt100` | `c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須)。run_0376 型 (plan §6) のうち推奨設定に合う run_0468 (run_0376 は dual-time で implicitRelax 0.7・SFR 2 なのに convMethod 0 で check_solver_config が FAIL)。200 step を 100 + 100 に分割 (前半) | 100 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `/CHECKPOINT/<cons>P` (26 データセット) |
| `c44dual_restart100` | `c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | species_db.yaml から H2O を外す。後半。`valueFileName` を ic.h5 → ic_ckpt.h5 (= base r1 の c44dual_ckpt100 の res_100.h5) | 100 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | log の履歴復元行 |
| `c44dual_pindiag` | `c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須)。run_0376 型 (plan §6) のうち推奨設定に合う run_0468 (run_0376 は dual-time で implicitRelax 0.7・SFR 2 なのに convMethod 0 で check_solver_config が FAIL)。200 step を 100 + 100 に分割 (前半) | 100 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1`, `FORGE_PIN_DIAG=1` | log のピン診断行 (入口あり) |
| `c52cht` | `c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | 無修正 (res_8000 から restart。CHT の内部状態は res に無いので界面は再初期化) | 200 | - | 壁 h5 (`iface*`)・`conjugate_*.csv` |
| `c36cell` | `c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | 旧キー → `model: sst`、`wallTreatmentSST` 1→0 | 200 | - | 出力・残差 (cell は CONNE を cells から組む) |
| `c20cell_rk3` | `c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | 無修正 (cell RK3、変換時の初期場から) | 200 | - | 出力・残差・壁 h5 |
| `c20cell_dual` | `c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | 旧キー LESorRANS 1 / LESmodel 1 → `model: wale` | 100 | - | 出力・残差・壁 h5 |
| `c20cell_impdiag` | `c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | 旧キー LESorRANS 0 → `model: none` (cell 陰解法。FORGE_IMPLICIT_DIAG_CSV の cell 版) | 200 | `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `diag.csv` |
| `c57lm` | `c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | 無修正 (mesh.h5 は run_0011 の場を種にしたもの = roGamma を読む経路) | 200 | - | `lm*`、log `[variables] transition … read from input` |
| `c57lm_fromsst` | `c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | run_0013 の LM 設定 + run_0012 (SST だけ) の mesh.h5 に res_20000 を写した場 (roGamma が無い = 初期化する経路) | 200 | - | log `… not in input: will be initialised` |
| `c56lineimp` | `c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | 無修正 (res_100000 から restart) | 200 | - | 出力の集合 (extraFields)、log `[lineImplicit] lines=…` |
| `c56extra` | `c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | 無修正 (extraFields に res_ro 等。元は 1 step → 200 step) | 200 | - | extraFields の `res_*` |
| `c48absorb` | `c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | 無修正 (extraFields に roN.. と dq_block_old_*。元は 1 step → 200 step) | 200 | - | extraFields の `roN..`・`dq_block_old_*` |
| `c26optin` | `c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | `wallTreatmentSST` 1→0 (env による名前の登録・除去は壁処理に依らない: variables.cpp:301-356) | 200 | - | データセット集合、警告 `… 確保されていない変数なので無視する` |
| `c26optin_env` | `c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | `wallTreatmentSST` 1→0 (env による名前の登録・除去は壁処理に依らない: variables.cpp:301-356) | 200 | `FORGE_WI_FORCE_DIAG=1`, `FORGE_WF_CLOSURE_DIAG=1`, `FORGE_OMEGA_BUDGET=1`, `FORGE_WF_REP_DIAG=1` | データセット集合 (wi_*・wf_g・omg_*・rep_*) |
| `v36node` | `v36node` ← `36.passive_pseudoshock_control` | c36node の config (node) で mesh/passive_solid.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v36cell` | `v36cell` ← `36.passive_pseudoshock_control` | c36cell の config (cell) で mesh/passive_solid.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v09` | `v09` ← `09.Taylor-Green` | c09ckpt100 の config (node 周期・2 種・トレーサ) で mesh/Taylor-Green.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v52` | `v52` ← `52.conjugate_slab` | case/52 mesh/ の変換 config のまま | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v44` | `v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | c44steady の config (軸対称・2 種・凝縮) で nozzle.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |

**元から差し替えた入力**: §6 の「`run_0376` 型」は `run_0376` 自体が dual-time で `implicitRelax: 0.7`・`speciesFaceReconstruction 2` なのに `convMethod 0`
(`check_solver_config.py` が FAIL) なので、同じ型 (軸対称 node・[MIXDRY, H2O]・凝縮・scheme 1・level 2・inletProfile) で推奨設定に合う
`run_0468_sweep_cflp12_nsub20_float` (PASS) を使った。

**起動のために複製側だけ直したもの** (全部 `inputs/*/SOURCE.txt` に記録):
- 旧キー `LESorRANS`/`LESmodel`/`RANSmodel` (現行 solver が起動時に拒否) → `turbulence.model` (c36node・c36cell・c20cell_dual・c20cell_impl)。
- `wallTreatmentSST: 1` (使用禁止) → 0 (c36node・c36cell・c26optin)。c26 の env 診断は名前の登録・除去が壁処理に依らない (`variables.cpp:301-356`) ので試験の目的は保たれる。
- case/44: `species_db.yaml` から H2O を外した (凝縮 ON では凝縮する H2O は組込みでなければならず、DB の H2O は MW 0.0180153 vs 組込み 0.01801528 で拒否される)。
  元の場に化学種の記録 (`species_hash`) が無いので、restart は `--force-species`、solver は `FORGE_ALLOW_UNVERIFIED_SPECIES=1` (監査の env 一覧で host `c` の名前集合に効かない側)。
  c44steady は変換直後の `nozzle.h5` に凝縮モーメントが無いので `rog_0`・`roQ*_0` を 0 で作ってから写した (13 量をビット一致で移送)。
- 同一メッシュ restart の種: c36node (res_80000)・c36cell (res_60000)・c44steady (res_24000)・c52cht (res_8000)・c56lineimp (res_100000)・
  c57lm_fromsst (run_0012 res_20000 を run_0012 の mesh.h5 へ)。いずれも `restart_field.py` (9c9f623c) の VERDICT OK (保存量ビット一致)。
- dual-time の再開 (c09restart100・c44dual_restart100) は **base r1 の ckpt100 の `res_100.h5` を全ビルド・全反復で同じ checkpoint として使う**
  (`ic_ckpt.h5`、`CKPT_FROM.txt` に sha256)。c09 は `run_0001` (sha256 `04a0f362…`)、c44 は `run_0029` (`1473d282…`)。

## 判定 (`compare_runs.py`、§6 を実装したもの。結果を見てから変えない)

- **(a)** `residual_history.csv` の最初の step の全行を全列ビット一致 (文字列一致) で比べる。base 3 回で値が割れた列は「new の値が base の観測値のどれか」で可。
  **初期出力** (最初の `res_*.h5`) は保存量・原始量・幾何量 (`T0_VALUE` の正規表現 + `/VALUE` 以外の全データセット) だけをビット一致で比べ、
  残差・補正量・勾配・リミッタ・診断は比べない (情報として不一致の名前だけ出す)。
- **(b)** 最終出力 (最終 `res_N.h5`・境界出力・CSV/probe 出力・残差履歴の各列) で m(A,B) = max|A−B| / max|A|。
  S = 同ビルド内ペア (base 3 対 + new 3 対) の最大、D = base×new 9 対の最大。合格は **D ≤ 2·S** (S = 0 なら D = 0)。new が無いときは S_base を出す。
- **(c)** 各 run の `NANCHECK.txt` (ワーカーが残差 CSV の数値列と run が書いた全 h5 の浮動小数データセットを検査)。
- **出力互換**: 出力ファイル集合、全 h5 のデータセット集合・shape・dtype・属性 (値まで) を全 run で比べる。
- **ログ**: 名前登録・checkpoint 復元・遷移初期化・ψ 退避件数 (`snapshot: X of Y`)・警告の行が全 run で一致すること。

**比較スクリプトの対象の直し (2026-10-07、new の比較で誤 FAIL・比較漏れが出た所。§6 の判定条件そのものは変えていない)**:
- ログは**行の集合**で比べる (順序は見ない)。変更後は環境変数の出力登録を確保の前へ移すので、`[FORGE_OUT_RESIDUALS]` 等の行の位置が変わるのは設計どおり (監査 §2 (i))。
- `[memlog]` 行 (FORGE_MEMLOG の計測値。"lineImplicit" で拾われ、RSS の実測値が run ごとに違う) と、変更後が設計どおり新しく出す
  `[variables] host cell arrays (gpu: 1): N of M registered` 行はログ比較から外す (後者の N は下の「結果」に別表で載せる)。
- 空白区切りの CSV (`conjugate_Tw_*.csv`) を列に分けられず比較から漏れていた → `,` が無い表は空白で切る。
- `mem_samples.csv` (ハーネスの採取で forge の出力ではない) は比較しない。
- 判定には使わない情報を足した: (a) の base 内 ulp 幅と new の最近傍 base 値までの ulp 距離、(b) の D/2S 最大のデータセット、FAIL の量の max|A|・max|A−B|。

上の (a)〜ログを**登録判定 A** と呼ぶ (2026-10-07 の結果は下の「結果」、FAIL のまま記録し書き換えない)。

### 追加診断 B (plan §6.2 で事前登録、2026-10-07。`compare_runs.py --metric abs`)

A の (b) の尺度 m = max|A−B|/max|A| は**分母がペアの片側**で、同ビルド内は `i<j` の片方向・ビルド間は base→new の全組合せなので
S と D の尺度が揃わず、**反復の並び順で判定が変わる** (base [100,2,1]・new [100,2,1] が FAIL、new を [1,2,100] に並べ替えると PASS)。
さらに `inf ≤ 2·inf` が真になり、**比較不能を PASS にする** (base [0,1,1]・new [100,100,100] が PASS)。codex (diagnose)
[`notes/reviews/2026-10-07-hostmem-regression-fail-diagnose.md`](../../notes/reviews/2026-10-07-hostmem-regression-fail-diagnose.md) の指摘、
最小再現は `test_compare_abs.py` に入れた。B は**既存の base 3 本・new 3 本の同じ保存時点・同じ配列・同じ CSV 行で尺度だけを替える** (forge の追加実行なし):

- d(A,B) = max|A − B| (float64、絶対 L∞。ペアに対称)。S_abs = 同ビルド内 6 対 (base 3 対 + new 3 対) の最大、D_abs = ビルド間 9 対の最大。
- 合格は **D_abs ≤ 2·S_abs** (S_abs = 0 なら D_abs = 0)。全入力が有限で shape・列・行キー (`step`・`inner`・`phase`・`Step`・`var`・`physID`) が
  対応すること — **比較不能は FAIL** (inf を PASS にしない)。整数・文字列の配列・列は全 run で厳密一致。データセットの集合・shape・dtype・属性も全 run で一致。
- 対象は A の (b) と同じ (最終出力 h5 の全データセット・境界出力・CSV/probe 出力・残差履歴の全列)。FAIL の 7 量だけでなく全量・全構成と SERN g3 を再評価。
- **自己検査**: 各量で (i) base・new それぞれの全順列 (3!×3! = 36 通り) と (ii) base/new の交換で S_abs・D_abs・判定が変わらないこと
  (d は**順序つき**の組で計算するので、非対称なら検出される)。(iii) 単体試験 `test_compare_abs.py` (codex の最小再現 2 つ、比較不能・欠落・
  shape 違い・S = 0・境界値・整数の厳密一致・乱数での並べ替え不変)。
- 出力: `results/<日付>_<base>_vs_<new>_abs/` に構成ごとの報告 (`<構成>.txt`、全量を D/2S の大きい順、FAIL には各 run の max|x| と最悪ペア)、
  全量の表 `all_quantities.tsv` (構成・ファイル・量・S_abs・D_abs・D/2S・判定・理由)、`summary.txt`、`selfcheck.txt`、`test_compare_abs.txt`。
  SERN g3 は `sern_g3/` (`--base … --new … --tag sern_g3`)。
- **読み方 (plan §6.2)**: B で超過が消えた量は、その A の FAIL を「尺度依存で説明できる (追加診断で反復内差の 2 倍以内)」と記録する。
  B でも超過する量は変更起因の差を候補に戻す。**B の PASS は A の書換えや「非決定性だけだった」証明には使わない**。
## 使い方 (変更後のビルドで回す手順)

前提: AWS が `running` (`bash solver_density_cuda/tools/aws_instance.sh status`、start/stop は自分でしない)、`pgrep -x forge` の cwd で他セッションを確認。

```bash
# 0) AWS 側の case/66 は ~/forge-b4 の未追跡ファイル。case/66 を含むコミットへ checkout するときは
#    同名の追跡ファイルと衝突するので `git checkout -f <commit>` にする (run_*・h5 は .gitignore 済みで残る)。
#    変更後のバイナリは ~/bin-hostmem/forge_<hash> に置く。
cd ~/forge-b4/case/66.hostmem_regression
python3 run_matrix.py status                     # running? (中断) / queued を確認
python3 run_matrix.py note <中断した run ...> --state 破棄予定 --text "AWS 停止で中断" --exclude
python3 run_matrix.py resume                     # queued のまま残った run を登録どおりに回す
# 1) base の不足分を足す (各構成 3 本そろうまで。番号は自動で連番)
python3 run_matrix.py launch <構成 ...> --build base --bin ~/bin-hostmem/forge_9c9f623c \
    --convert-bin ~/bin-hostmem/convert_9c9f623c --reps 2 3
# 2) 変更後 (ビルド名は段階ごとに変える: r1 / r2 / r3 / new 等)
python3 run_matrix.py launch $(python3 -c "import matrix_spec as m; print(' '.join(m.CONFIGS))") \
    --build new --bin ~/bin-hostmem/forge_<hash> --convert-bin ~/bin-hostmem/convert_<hash> --reps 1 2 3
python3 run_matrix.py status --build new; python3 run_matrix.py verify
# 3) 判定 (構成ごとの報告 + summary.txt)
python3 compare_runs.py --all --new-build new --out-dir compare_new
python3 compare_runs.py --diff2 run_<c09restart100>/res_100.h5 run_<c09cont200>/res_200.h5   # 分割と連続
# 4) テキストだけ持ち帰る (除外を先に書く)
rsync -az --exclude='*.h5' --exclude='residual_history.csv' --include='*/' --include='*.txt' --include='registry.tsv' \
    --include='summary.txt' --exclude='*' ubuntu@<IP>:forge-b4/case/66.hostmem_regression/ case/66.hostmem_regression/_aws/
python3 run_matrix.py table   # → 下の「計算 run 一覧」を置き換える
```

SERN g3 (メモリ): `case/46` に `run_1072` と同じ作り方で run を作り (sern.h5 はハードリンクで共有・読み取り専用)、
`python3 run_matrix.py launch-dir <run ...> --bin <bin> --label <ビルド名> --env FORGE_MEMLOG=1 --memwatch`、
`python3 memlog_summary.py <run ...> --items` で工程別の表。比較は `compare_runs.py --base <base 3 本> --new <new 3 本>`。

## 結果 (2026-10-07、base 9c9f623c 3 回 対 new 93e55957 3 回)

原本: [`results/2026-10-07_base9c9f623c_vs_new93e55957/`](results/2026-10-07_base9c9f623c_vs_new93e55957/summary.txt)。全 run `VERIFY: PASS` (RUN_PROVENANCE の forge_bin・sha256 が登録どおり)。
短い回帰 run なので収束の主張には使わない (全 forge run の check_convergence は NOT CONVERGED、回帰差の判定だけに使う)。

| 構成 | (a) step 0 | step 0 の ulp (base 内の幅 / new の最近傍距離) | 初期出力 | (b) D ≤ 2S | (b) D/2S 最大のデータセット | (c) NaN | 出力互換 | ログ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `c36node` | FAIL (9 値) | 44 / 13 | PASS | PASS | `res_wall_4_200.h5:VALUE/utau` 0.67 (D 4.35e-04, S 3.27e-04) | PASS | PASS | PASS |
| `c36node_impdiag` | FAIL (18 値) | 75 / 33 | PASS | PASS | `point_probe_2.out:T` 0.75 (D 3.22e-07, S 2.15e-07) | PASS | PASS | PASS |
| `c36node_psidual` | FAIL (15 値) | 42 / 5 | PASS | PASS | `psi_dualeval.csv:S_ApmA` 0.92 (D 1.10e+00, S 5.97e-01) | PASS | PASS | PASS |
| `c09ckpt100` | FAIL (300 値) | 4476783 / 3664948 | PASS | PASS | `res_100.h5:VALUE/dUydx` 0.60 (D 5.19e-01, S 4.32e-01) | PASS | PASS | PASS |
| `c09restart100` | FAIL (315 値) | 3575316 / 3435028 | PASS | PASS | `residual_history.csv:rms_roUx` 0.55 (D 5.97e-02, S 5.46e-02) | PASS | PASS | PASS |
| `c09cont200` | FAIL (316 値) | 3744452 / 3408329 | PASS | PASS | `residual_history.csv:rms_roXi` 0.75 (D 1.31e-07, S 8.75e-08) | PASS | PASS | PASS |
| `c09ckpt100_outres` | FAIL (312 値) | 3952886 / 2146423 | PASS | PASS | `res_100.h5:VALUE/res_roUy` 0.75 (D 3.00e+00, S 2.00e+00) | PASS | PASS | PASS |
| `c09ckpt100_rawdiag` | FAIL (311 値) | 4208959 / 2727567 | PASS | PASS | `res_100.h5:VALUE/res_roUy` 0.58 (D 1.75e+00, S 1.50e+00) | PASS | PASS | PASS |
| `c09ckpt100_pindiag` | FAIL (315 値) | 3959569 / 2485060 | PASS | PASS | `res_100.h5:VALUE/res_roUy` 0.67 (D 2.00e+00, S 1.50e+00) | PASS | PASS | PASS |
| `c44steady` | FAIL (54 値) | 13561 / 21961 | PASS | PASS | `res_200.h5:VALUE/res_ro` 0.95 (D 2.73e+00, S 1.43e+00) | PASS | PASS | PASS |
| `c44dual_ckpt100` | FAIL (627 値) | 1850076 / 1176704 | PASS | FAIL (2 量) | `res_100.h5:VALUE/condClampCorrQ_0` 2016.58 (D 1.79e+05, S 4.45e+01) | PASS | PASS | PASS |
| `c44dual_restart100` | FAIL (628 値) | 1301738 / 1454616 | PASS | FAIL (1 量) | `res_100.h5:VALUE/condClampCorrQ_0` 2192.98 (D 4.39e+03, S 1.00e+00) | PASS | PASS | PASS |
| `c44dual_pindiag` | FAIL (626 値) | 2482237 / 1193940 | PASS | PASS | `res_100.h5:CHECKPOINT/rog_0_fctH` 0.87 (D 5.23e-04, S 3.02e-04) | PASS | PASS | PASS |
| `c52cht` | FAIL (18 値) | 19 / 29 | PASS | PASS | `res_wall_bot_3_200.h5:VALUE/ifaceRro` 0.72 (D 4.82e-05, S 3.35e-05) | PASS | PASS | PASS |
| `c36cell` | FAIL (21 値) | 28 / 9 | PASS | PASS | `res_200.h5:VALUE/roUx` 0.60 (D 3.95e-04, S 3.28e-04) | PASS | PASS | PASS |
| `c20cell_rk3` | FAIL (14 値) | 18051 / 7617 | PASS | FAIL (1 量) | `residual_history.csv:rms_roUz` 4.00 (D 1.02e+03, S 1.28e+02) | PASS | PASS | PASS |
| `c20cell_dual` | FAIL (172 値) | 1025677 / 984352 | PASS | FAIL (3 量) | `res_100.h5:VALUE/Uz` 1.68 (D 8.87e+01, S 2.65e+01) | PASS | PASS | PASS |
| `c20cell_impdiag` | PASS | 0 / 0 | PASS | PASS | `residual_history.csv:rms_roUx` 0.61 (D 1.69e-06, S 1.39e-06) | PASS | PASS | PASS |
| `c57lm` | FAIL (72 値) | 38642 / 24679 | PASS | PASS | `res_200.h5:VALUE/lmPgamma` 0.75 (D 2.10e-05, S 1.40e-05) | PASS | PASS | PASS |
| `c57lm_fromsst` | FAIL (54 値) | 21209 / 15991 | PASS | PASS | `residual_history.csv:rms_roOmega` 0.65 (D 1.76e-01, S 1.35e-01) | PASS | PASS | PASS |
| `c56lineimp` | FAIL (27 値) | 48 / 40 | PASS | PASS | `res_200.h5:VALUE/roOmega` 0.85 (D 2.35e-05, S 1.38e-05) | PASS | PASS | PASS |
| `c56extra` | FAIL (27 値) | 66 / 33 | PASS | PASS | `res_200.h5:VALUE/roOmega` 0.76 (D 8.69e-06, S 5.72e-06) | PASS | PASS | PASS |
| `c48absorb` | FAIL (54 値) | 4409 / 3672 | PASS | PASS | `residual_history.csv:rms_roK` 0.68 (D 4.41e-03, S 3.24e-03) | PASS | PASS | PASS |
| `c26optin` | FAIL (21 値) | 199 / 239 | PASS | PASS | `res_200.h5:VALUE/k` 0.57 (D 1.83e-05, S 1.60e-05) | PASS | PASS | PASS |
| `c26optin_env` | FAIL (24 値) | 94 / 146 | PASS | PASS | `res_200.h5:VALUE/roUy` 0.74 (D 1.44e-03, S 9.75e-04) | PASS | PASS | PASS |
| `v36node` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v36cell` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v09` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v52` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v44` | - | - | - | PASS | - | PASS | PASS | PASS |
| `SERN g3 (case/46 run_1072–1077)` | FAIL (45 値) | 52 / 15 | PASS | PASS | `res_vehicle_base_18_100.h5:VALUE/twall_z` 1.00 (D 4.11e-07, S 2.05e-07) | PASS | PASS | PASS |

**まとめ**: 初期出力 (保存量・原始量・幾何量のビット一致)・(c) NaN・出力互換 (ファイル集合・データセット集合・shape・dtype・属性値)・ログ行は**全構成 PASS**。
変換器 5 構成は出力 h5 が base 同士・new 同士・base 対 new とも全データセットで差 0。
ψ 二重評価の退避件数は new 3 回とも `snapshot: 269 of 269 device arrays` (base と同じ = R2 の pdeSize の fallback 修正で面配列が漏れていない)。
dual-time の再開 (c09restart100・c44dual_restart100) は new でも `history restored from ic_ckpt.h5 (/CHECKPOINT: 10 levels …)` と FCT 履歴の復元行が出る。

**(a) step 0 は §6 の文言どおりでは c20cell_impdiag 以外の全構成 FAIL**。§6 の前提「既知の 1 ulp の非決定性で、変更前ビルド同士でも同じ 2 値が出る」が
base 3 回の時点で成り立っていない: 最初の組立の残差でも 3 回で 3 値に割れる列があり、幅は 2D node 定常で 19〜75 ulp、SERN g3 で 52 ulp、
dual-time は step 0 に内反復 (22 行) が入るので 10⁶ ulp 級。new の値と最近傍 base 値の距離は大半の構成で base 内の幅以下、
上回るのは c44steady (21961 / 13561)・c52cht (29 / 19)・c26optin (239 / 199)・c26optin_env (146 / 94)・c44dual_restart100 (1454616 / 1301738)。
c52cht の `rms_roUx` は base 3 回・new 3 回がそれぞれ揃って 1 ulp 違ったので、追加で base・new 各 4 回 (run_0186–0193、比較から除外) を回した:
base 7 回中 4 回が new と同じ値を出す (= 系統差でなく 2 値の非決定性)。判定条件の扱い (ulp 幅で読み替えるか) は plan 側の決定事項 (ここでは変えない)。

**(b) の FAIL (4 構成・7 量)** — いずれも名目上ゼロか疎な診断量で、base 自身の反復の間で桁が動いている:

| 構成 | データセット | S (base / new) | D | run ごとの最大絶対値 (base r1–r3 / new r1–r3) |
| --- | --- | --- | --- | --- |
| c44dual_ckpt100 | `res_100:condClampCorrQ_0` (凝縮モーメントのクランプ補正、非零 ~1400 節点) | 44.5 / 1.0 | 1.79e5 | 8.2e17, 8.1e17, 3.6e19 / **1.5e23**, 8.1e17, 1.7e17 (c44dual_pindiag の base r3 も 1.45e23) |
| c44dual_ckpt100 | `res_100:condR30_0` | 1.04 / 1.0 | 3.14 | 1.6e-5, 1.3e-5, 1.4e-5 / 4.2e-5, 2.5e-5, 1.3e-5 |
| c44dual_restart100 | `res_100:condClampCorrQ_0` | 1.0 / 1.0 | 4.39e3 | 1.4e10, 8.0e6, 1.8e3 / 7.9e6, 7.9e6, 3.1e2 |
| c20cell_rk3 | `residual_history.csv:rms_roUz` (擬似 2D の spanwise 残差、`rms_roUx` 6e-3 に対し 1e-9〜3e-7) | 1.0 / 128 | 1.02e3 | (残差履歴の列) |
| c20cell_dual | `res_100:Uz`・`roUz`・`CHECKPOINT/roUzN` (擬似 2D の spanwise 速度) | 26.5 / 1.3 | 88.7 | roUz: 1.7e-8, **4.3e-7**, 5.7e-9 / **4.2e-7**, 5.8e-9, 5.7e-9 (roUx は 156) |

保存量・原始量の主要データセット (ro・roUx・roe・P・T・ρY・凝縮モーメント・k/ω・γ/Re_θt) は全構成で D ≤ 2S。

**dual-time の分割と連続** ([`split_vs_cont_c09.txt`](results/2026-10-07_base9c9f623c_vs_new93e55957/split_vs_cont_c09.txt)): c09 は一様流 (ro・roUx・roe・P・T・ρY は
分割・連続・反復のどれでもビット一致) で、動くのはトレーサ roXi と名目ゼロの roUy・roUz。roXi の m (分割 vs 連続 / 連続内 / 分割内) は
base 3.31e-6 / 2.81e-6 / 6.0e-7、new 2.71e-6 / 2.76e-6 / 6.0e-7 (checkpoint の roXiP も同程度) — 分割と連続の差は反復のばらつきと同じ大きさで、ビルド間で変わらない。

**SERN g3 のメモリ** (192 万節点、本段設定 100 step、`FORGE_MEMLOG=1`、[`sern_g3_memlog.txt`](results/2026-10-07_base9c9f623c_vs_new93e55957/sern_g3_memlog.txt)。3 回とも工程別に同値):

| 工程 | base RSS / HWM [MiB] | new RSS / HWM [MiB] | GPU 使用 [MiB] (両方) |
| --- | --- | --- | --- |
| readMesh 末尾 (平坦な読込配列が生存) | 2614 / 2614 | 2614 / 2614 | — |
| initMatrix の後 | 2721 / 2721 | 2143 / 2614 | 260 |
| allocVariables の後 | 4864 / 4864 | 2344 / 2614 | 2440 |
| setStructuralVariables の後 | 4865 / **5157** | 2338 / **2630** | 2602 |
| 主ループの後 / 終了 | 5001 / 5157 | 2474 / 2630 | 2982 |

- ホスト VmHWM **5157 → 2630 MiB (2816 → 1436 B/節点、−49 %)**。終了時の内訳 (new): host `c` 23/223 本 183 MB (base 1774)、host `p` 0 (357)、`mat_ns` 0 (572)。
  ピークは readMesh 末尾 (2614) のすぐ上の setStructuralVariables (2630) に移った。
- GPU は工程別に base と同値 (プロセス 2974 MiB、初回 memlog から +2724 MB) → §6 の「GPU 傾き ±2 % 以内」を満たす。
- ローカルの縮小格子 s070/s085 からの外挿 (HWM 傾き 2737 → 1371 B/節点、切片 163 → 133 MiB) では g3 で base 5175・new 2644 MiB の見込みで、
  実測 5157・2630 MiB と −0.3 %・−0.5 % で合う。§6 の「ホスト VmHWM の傾き ≤ 1500 B/節点」は g3 単点の平均 1436 B/節点 (切片込み) で下回る
  (傾きの判定は 2 格子が要り、g4 が無いので縮小格子の傾き 1371 を併用)。

**new の host cell 確保数** (`[variables] host cell arrays (gpu: 1): N of M registered`、new r1 の起動ログ):
c36node 19/191・impdiag 20/191・c09ckpt100 98/241 (outres 114/241・rawdiag 100/243・pindiag 102/241)・c44steady 120/303・c44dual 133/303 (pindiag 140/303)・
c52cht・c36cell・c20cell_rk3 19/191・c20cell_dual 26/191・c57lm 92/211・c57lm_fromsst 24/202・c56lineimp 23/191・c56extra 25/191・c48absorb 29/191・
c26optin 20/191 (env 23/212)・SERN g3 23/223。

### 追加診断 B の結果 (2026-10-07、同じ 6 本、forge の追加実行なし)

原本: [`results/2026-10-07_base9c9f623c_vs_new93e55957_abs/`](results/2026-10-07_base9c9f623c_vs_new93e55957_abs/summary.txt)
(AWS で `compare_runs.py --metric abs --all --new-build new` と SERN g3 の直接指定、所要 66 s = 30 構成 35 s + SERN 31 s)。

- **全 30 構成 3147 量・SERN g3 255 量、FAIL 0** (比較不能 0、整数・文字列の厳密一致 185 量はすべて一致)。
- **自己検査**: 全 3402 量 × 並べ替え 36 通り (計 122472 評価) で S_abs・D_abs・判定が変わったもの 0、base/new の交換で変わったもの 0
  ([`selfcheck.txt`](results/2026-10-07_base9c9f623c_vs_new93e55957_abs/selfcheck.txt)・[`sern_g3/selfcheck.txt`](results/2026-10-07_base9c9f623c_vs_new93e55957_abs/sern_g3/selfcheck.txt))。
  単体試験 `test_compare_abs.py` は AWS・ローカルとも ALL PASS (8 試験。A の 2 つの欠陥も再現する)。
- **A で FAIL だった 7 量の B の値** (いずれも B では PASS。A の FAIL は**尺度依存で説明できる**と記録する):

  | 構成 | 量 | S_abs | D_abs | D/2S | S・D を決めている run (各 run の max|x| は「結果」の (b) の表) |
  | --- | --- | --- | --- | --- | --- |
  | c44dual_ckpt100 | `res_100:condClampCorrQ_0` | 1.462e23 | 1.462e23 | 0.50 | new r1 (1.46e23) が S (new 内) と D の両方を決める |
  | c44dual_ckpt100 | `res_100:condR30_0` | 4.203e-5 | 4.203e-5 | 0.50 | new r1 (4.2e-5) |
  | c44dual_restart100 | `res_100:condClampCorrQ_0` | 1.433e10 | 1.433e10 | 0.50 | base r1 (1.4e10) |
  | c20cell_rk3 | `residual_history.csv:rms_roUz` | 3.277e-7 | 3.277e-7 | 0.50 | (残差履歴の列) |
  | c20cell_dual | `res_100:roUz` | 4.315e-7 | 4.314e-7 | 0.50 | base r2 (4.3e-7) と new r1 (4.2e-7) |
  | c20cell_dual | `res_100:Uz` | 3.654e-7 | 3.653e-7 | 0.50 | 同上 |
  | c20cell_dual | `CHECKPOINT/roUzN` | 4.316e-7 | 4.315e-7 | 0.50 | 同上 |

  D/2S = 0.50 は「1 本だけ飛び抜けた run が S と D の両方を決めている」形で、B はこの量について base と new を区別できていない
  (その run が base 側か new 側かに依らず同じ値になる)。plan §6.2 のとおり、B の PASS を「非決定性だけだった」証明には使わない。
- **D/2S の上位 10** (全部 PASS): SERN g3 `res_vehicle_base_18_100:twall_z` **1.000** (S_abs 7.812e-3、D_abs 1.562e-2 = 2·S_abs ちょうど。
  float32 の値の刻みでちょうど 2 倍になっている)、SERN g3 `res_sidewall_in_11_100:utau` 0.957、SERN g3 `res_cowl_in_5_100:ypls` 0.950、
  SERN g3 `res_sidewall_in_11_100:omegab` 0.914、c44dual_pindiag `CHECKPOINT/rog_0_fctH` 0.866、c44dual_pindiag `CHECKPOINT/rog_0_fctHsrc` 0.856、
  c56lineimp `res_200:roOmega` 0.850、c56lineimp `res_gap_6_200:qwall` 0.845、SERN g3 `res_sidewall_out_12_100:twall_y` 0.827、
  c44dual_restart100 `res_100:sonic` 0.820。30 構成だけの上位と SERN の上位は `summary.txt`・`sern_g3/summary.txt`。

### 固定幅の独立 A/B (plan §6.3、構成 c44dual_ckpt100)

- **段階 1 (2026-10-07 済み、forge の実行なし)**: 既存 base 3 本 (run_0029・0038・0068) だけから $S_0 = \max_{i<j}\|B_i - B_j\|_\infty$、
  $T = 2S_0$ を凍結した ([`fixedwidth_c44dual_ckpt100/T_frozen.tsv`](fixedwidth_c44dual_ckpt100/T_frozen.tsv)、sha256 `45168ae5adf4…`)。
  最終出力の量 196 (幅 121、差 0 を要求 72、厳密一致 3) + 行キー 1 + 初期出力 44。new は使っていない。
  段階 2 の投入計画 (run_0194–0205 を b,n 交互、バイナリ・入力の sha256、評価法、判定表) は [`PLAN.txt`](fixedwidth_c44dual_ckpt100/PLAN.txt)・`plan.json`。
  追加記録 (twall_z の最悪差の位置と ULP、condClampCorrQ_0 の最大位置) は [`notes.txt`](fixedwidth_c44dual_ckpt100/notes.txt)。
- **段階 2 (2026-10-07、commit 55aaf586 で凍結した計画どおり)**: **判定 (plan §6.3): 新規 base も超過 → 判定不能**
  (基準 3 本では再現性を捉えられていない。T・本数・順序・対象は変えていない)。原本 [`RESULT.txt`](fixedwidth_c44dual_ckpt100/RESULT.txt)。
  - run_0194–0205 (b,n 交互、各 100 step、83 s で 12 本)。全 rc 0・NaN PASS・verify PASS、構造・行キー・初期出力は全本で既存 B1 と一致。中断・再投入なし。
  - 全量が幅内だったのは新規 base 2/6・new 2/6。幅を超えた量 (本数 b/n): `condClampCorrQ_0` 3/4、`passiveFloorCorr_g_0` 4/3、
    `passiveFloorCorr_Q0_0` 3/3、`passiveLimCorr_Q2_0` 3/1、`condR30_0` 1/1、new 1 本だけ `condClampCorr_0`・`Q0_0`・`roQ0_0`・`CHECKPOINT/roQ1_0_fctH`。
  - 記述: `condClampCorrQ_0` > 1e22 (節点 19954、1.45〜1.46e23) に入った本数は**新規 base 2/6・new 2/6** (T の算定に使った既存 base は 0/3)。
    全本で最大値 × 1e-30 がその節点の `roQ1_0` か `roQ2_0` と相対 4e-8 以下で一致 (分母床 1e-30 と整合)。

## 既知の注意

- **AWS の自動停止 (2026-10-07 01:05 UTC 頃 = JST 10:05、base 投入中に停止した)**: idle 自動停止 (`solver_density_cuda/tools/cloud/idle_autostop.sh`、5 分おき・6 回連続で停止) は
  「GPU 使用率 0・`pgrep -x forge` 0 件・ログイン 0・load < 1」を idle と数える。base は `forge_9c9f623c` という名前で起動していたので
  **`pgrep -x forge` に一致せず**、短い 2D run の合間 (GPU 使用率の瞬間値 0) と準備作業の時間で 30 分が数えられたと推定する (停止の記録は見ていない)。
  対策として `run_matrix.py` は `--bin` を `.bin/<ビルド名>/forge` というシンボリックリンク経由で起動するようにした (comm が `forge` になる)。
  再起動後の最初の投入で `pgrep -x forge` が自分の run に一致する (comm=forge、exe は実ファイル) ことを確認した (2026-10-07)。
- **`~/forge-b4` の `git checkout -f` で `solver_density_cuda/third_party/HighFive` が空になった** (この worktree の HighFive は
  `.git` ファイルがローカル PC のパスを指す複製だった)。同じ pin (`dfc06537`) の `~/forge-sern` から中身を複製して戻した
  (`third_party/HighFive.COPIED_FROM.txt`)。次に checkout するときも空になったら同じ手順で戻す。
- 変換器は終了時の `cudaFree` で GPUassert (invalid argument) を出して exit 1 になる既知の罠がある (出力 h5 は完全)。ワーカーは rc でなく
  出力 h5 が開けて `/VALUE`・`/MESH` を持つかで完了を判定する (`written=1`)。cell モードは `/VIZMESH` を書かない。
- case/36 の c36node は 80223 節点・200 step で約 0.3 s、全 2D 構成で 1 本数秒〜十数秒。ワーカーは 1 本ずつ順に回す。
- `run_case.sh` は `~/forge-b4` のものを使うので `RUN_PROVENANCE.txt` の `git_head` は `~/forge-b4` の HEAD (base では 9c9f623c)。

## 保留

- **SERN g4**: `run_1045`/`run_1046` (と延長 `run_1048`/`run_1049`) の最終場・g4 の `sern.h5` は AWS に残っていない (場の h5 は削除済み、
  g4 メッシュも無い)。g4 を使うには格子の再生成 (`problem_3d_prod_m6on_g4.yaml`) と g3 からの cross-mesh restart か段階起動が要る。
  2 格子での傾きは 2026-10-07 の縮小格子 3 点 (`notes/investigations/2026-10-07-forge-memlog/`) でしか出せない。
- **NaN ダンプ (表の P)**: 入力 (s050 縮小格子) が別セッションの scratch にあるので入れていない。
- case/13 `run_slau` (cell) は入れていない (cell は case/36・case/20 の 3 構成で代表)。
- (a) step 0 の判定条件の扱いと (b) の FAIL 7 量 (名目ゼロ・疎な診断量) の扱いは plan の決定事項 (上の「結果」)。

## 計算 run 一覧

AWS `~/forge-b4/case/66.hostmem_regression/` (2026-10-07)。base = `.bin/base/forge` → `~/bin-hostmem/forge_9c9f623c` (r1 の 34 本だけは実ファイル名のまま起動)、
new = `.bin/new/forge` → `~/bin-hostmem/forge_93e55957`、全部 `FORGE_CUDA_BLOCKSIZE=128`。N step の短い回帰 run なので収束判定は NOT CONVERGED
(回帰差の基準であって収束の主張には使わない)。判定済みの run からは入力 h5 の複製を削除した (`INPUTS_REMOVED.txt`、テンプレートは `inputs/` に残る)。
出力は `res_0`・最終 `res_N`・境界出力・CSV だけなので消していない。表は `python3 run_matrix.py table` の出力。

| run | 目的・設定差分 | 元の入力 | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0001_c09ckpt100_base_r1` | dual-time 100 step (checkpoint 書出し) (base r1) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0002_c44dual_ckpt100_base_r1` | 軸対称・凝縮 dual-time 前半 (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 1, NaN PASS, ?。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0003_c36node_base_r1` | 2D node 標準 (base r1) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0004_c36node_impdiag_base_r1` | 2D node + 陰解法診断 CSV (表 M) (base r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0005_c36node_psidual_base_r1` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r1) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0006_c09cont200_base_r1` | dual-time 連続 200 step (base r1) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0007_c09ckpt100_outres_base_r1` | dual-time + 残差出力 (表 K・L) (base r1) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0008_c09ckpt100_rawdiag_base_r1` | dual-time + roYraw (表 I) (base r1) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0009_c09ckpt100_pindiag_base_r1` | dual-time + ピン診断 (表 N) (base r1) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0010_c44steady_base_r1` | 軸対称・多成分・凝縮 (定常) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS, ?。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0011_c44dual_pindiag_base_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 1, NaN PASS, ?。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0012_c52cht_base_r1` | 共役伝熱 (base r1) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0013_c36cell_base_r1` | cell 定常 SST (base r1) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0014_c20cell_rk3_base_r1` | cell RK3 (base r1) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0015_c20cell_dual_base_r1` | cell LES dual-time (base r1) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0016_c20cell_impdiag_base_r1` | cell 陰解法 + 診断 CSV (base r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0017_c57lm_base_r1` | 遷移 (roGamma を読む) (base r1) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0018_c57lm_fromsst_base_r1` | 遷移 (SST の場から初期化) (base r1) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0019_c56lineimp_base_r1` | line-implicit + extraFields (base r1) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0020_c56extra_base_r1` | extraFields (res_* を含む) + FP64 アキュムレータ (base r1) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0021_c48absorb_base_r1` | extraFields (roN・dq_block_old_*) (base r1) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0022_c26optin_base_r1` | env 診断なし (名前が消えて警告) (base r1) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0023_c26optin_env_base_r1` | env 診断 4 種 ON (表 I) (base r1) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0024_v36node_base_r1` | 変換器 node (base r1) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0025_v36cell_base_r1` | 変換器 cell (base r1) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0026_v09_base_r1` | 変換器 node 周期・種 (base r1) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0027_v52_base_r1` | 変換器 node (CHT 用スラブ) (base r1) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0028_v44_base_r1` | 変換器 軸対称・種・凝縮 (base r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0029_c44dual_ckpt100_base_r1` | 軸対称・凝縮 dual-time 前半 (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0030_c44steady_base_r1` | 軸対称・多成分・凝縮 (定常) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0031_c44dual_pindiag_base_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0032_v44_base_r1` | 変換器 軸対称・種・凝縮 (base r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0033_c09restart100_base_r1` | dual-time 再開 100 step (checkpoint 復元) (base r1) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0034_c44dual_restart100_base_r1` | 軸対称・凝縮 dual-time 後半 (再開) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0035_c09ckpt100_base_r2` | dual-time 100 step (checkpoint 書出し) (base r2) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | done rc=? (中断) nan=-, NOT CONVERGED。AWS 停止 (2026-10-07 01:05 UTC 頃) で中断 | 破棄予定 |
| `run_0036_c09restart100_base_r2` | dual-time 再開 100 step (checkpoint 復元) (base r2) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0037_c09cont200_base_r2` | dual-time 連続 200 step (base r2) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0038_c44dual_ckpt100_base_r2` | 軸対称・凝縮 dual-time 前半 (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0039_c44dual_restart100_base_r2` | 軸対称・凝縮 dual-time 後半 (再開) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0040_c36node_base_r2` | 2D node 標準 (base r2) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0041_c36node_impdiag_base_r2` | 2D node + 陰解法診断 CSV (表 M) (base r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0042_c36node_psidual_base_r2` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r2) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0043_c09ckpt100_outres_base_r2` | dual-time + 残差出力 (表 K・L) (base r2) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0044_c09ckpt100_rawdiag_base_r2` | dual-time + roYraw (表 I) (base r2) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0045_c09ckpt100_pindiag_base_r2` | dual-time + ピン診断 (表 N) (base r2) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0046_c44steady_base_r2` | 軸対称・多成分・凝縮 (定常) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0047_c44dual_pindiag_base_r2` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0048_c52cht_base_r2` | 共役伝熱 (base r2) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0049_c36cell_base_r2` | cell 定常 SST (base r2) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0050_c20cell_rk3_base_r2` | cell RK3 (base r2) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0051_c20cell_dual_base_r2` | cell LES dual-time (base r2) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0052_c20cell_impdiag_base_r2` | cell 陰解法 + 診断 CSV (base r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0053_c57lm_base_r2` | 遷移 (roGamma を読む) (base r2) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0054_c57lm_fromsst_base_r2` | 遷移 (SST の場から初期化) (base r2) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0055_c56lineimp_base_r2` | line-implicit + extraFields (base r2) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0056_c56extra_base_r2` | extraFields (res_* を含む) + FP64 アキュムレータ (base r2) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0057_c48absorb_base_r2` | extraFields (roN・dq_block_old_*) (base r2) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0058_c26optin_base_r2` | env 診断なし (名前が消えて警告) (base r2) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0059_c26optin_env_base_r2` | env 診断 4 種 ON (表 I) (base r2) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0060_v36node_base_r2` | 変換器 node (base r2) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0061_v36cell_base_r2` | 変換器 cell (base r2) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0062_v09_base_r2` | 変換器 node 周期・種 (base r2) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0063_v52_base_r2` | 変換器 node (CHT 用スラブ) (base r2) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0064_v44_base_r2` | 変換器 軸対称・種・凝縮 (base r2) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0065_c09ckpt100_base_r3` | dual-time 100 step (checkpoint 書出し) (base r3) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0066_c09restart100_base_r3` | dual-time 再開 100 step (checkpoint 復元) (base r3) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0067_c09cont200_base_r3` | dual-time 連続 200 step (base r3) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0068_c44dual_ckpt100_base_r3` | 軸対称・凝縮 dual-time 前半 (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0069_c44dual_restart100_base_r3` | 軸対称・凝縮 dual-time 後半 (再開) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0070_c36node_base_r3` | 2D node 標準 (base r3) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0071_c36node_impdiag_base_r3` | 2D node + 陰解法診断 CSV (表 M) (base r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0072_c36node_psidual_base_r3` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r3) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0073_c09ckpt100_outres_base_r3` | dual-time + 残差出力 (表 K・L) (base r3) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0074_c09ckpt100_rawdiag_base_r3` | dual-time + roYraw (表 I) (base r3) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0075_c09ckpt100_pindiag_base_r3` | dual-time + ピン診断 (表 N) (base r3) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0076_c44steady_base_r3` | 軸対称・多成分・凝縮 (定常) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0077_c44dual_pindiag_base_r3` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0078_c52cht_base_r3` | 共役伝熱 (base r3) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0079_c36cell_base_r3` | cell 定常 SST (base r3) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0080_c20cell_rk3_base_r3` | cell RK3 (base r3) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0081_c20cell_dual_base_r3` | cell LES dual-time (base r3) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0082_c20cell_impdiag_base_r3` | cell 陰解法 + 診断 CSV (base r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0083_c57lm_base_r3` | 遷移 (roGamma を読む) (base r3) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0084_c57lm_fromsst_base_r3` | 遷移 (SST の場から初期化) (base r3) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0085_c56lineimp_base_r3` | line-implicit + extraFields (base r3) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0086_c56extra_base_r3` | extraFields (res_* を含む) + FP64 アキュムレータ (base r3) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0087_c48absorb_base_r3` | extraFields (roN・dq_block_old_*) (base r3) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0088_c26optin_base_r3` | env 診断なし (名前が消えて警告) (base r3) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0089_c26optin_env_base_r3` | env 診断 4 種 ON (表 I) (base r3) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0090_v36node_base_r3` | 変換器 node (base r3) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0091_v36cell_base_r3` | 変換器 cell (base r3) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0092_v09_base_r3` | 変換器 node 周期・種 (base r3) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0093_v52_base_r3` | 変換器 node (CHT 用スラブ) (base r3) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0094_v44_base_r3` | 変換器 軸対称・種・凝縮 (base r3) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0095_c09ckpt100_base_r2` | dual-time 100 step (checkpoint 書出し) (base r2) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0096_c36node_new_r1` | 2D node 標準 (new r1) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0097_c36node_impdiag_new_r1` | 2D node + 陰解法診断 CSV (表 M) (new r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0098_c36node_psidual_new_r1` | 2D node + ψ 二重評価 (R2 の pdeSize) (new r1) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0099_c09ckpt100_new_r1` | dual-time 100 step (checkpoint 書出し) (new r1) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0100_c09restart100_new_r1` | dual-time 再開 100 step (checkpoint 復元) (new r1) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0101_c09cont200_new_r1` | dual-time 連続 200 step (new r1) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0102_c09ckpt100_outres_new_r1` | dual-time + 残差出力 (表 K・L) (new r1) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0103_c09ckpt100_rawdiag_new_r1` | dual-time + roYraw (表 I) (new r1) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0104_c09ckpt100_pindiag_new_r1` | dual-time + ピン診断 (表 N) (new r1) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0105_c44steady_new_r1` | 軸対称・多成分・凝縮 (定常) (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0106_c44dual_ckpt100_new_r1` | 軸対称・凝縮 dual-time 前半 (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0107_c44dual_restart100_new_r1` | 軸対称・凝縮 dual-time 後半 (再開) (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0108_c44dual_pindiag_new_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0109_c52cht_new_r1` | 共役伝熱 (new r1) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0110_c36cell_new_r1` | cell 定常 SST (new r1) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0111_c20cell_rk3_new_r1` | cell RK3 (new r1) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0112_c20cell_dual_new_r1` | cell LES dual-time (new r1) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0113_c20cell_impdiag_new_r1` | cell 陰解法 + 診断 CSV (new r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0114_c57lm_new_r1` | 遷移 (roGamma を読む) (new r1) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0115_c57lm_fromsst_new_r1` | 遷移 (SST の場から初期化) (new r1) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0116_c56lineimp_new_r1` | line-implicit + extraFields (new r1) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0117_c56extra_new_r1` | extraFields (res_* を含む) + FP64 アキュムレータ (new r1) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0118_c48absorb_new_r1` | extraFields (roN・dq_block_old_*) (new r1) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0119_c26optin_new_r1` | env 診断なし (名前が消えて警告) (new r1) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0120_c26optin_env_new_r1` | env 診断 4 種 ON (表 I) (new r1) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0121_v36node_new_r1` | 変換器 node (new r1) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0122_v36cell_new_r1` | 変換器 cell (new r1) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0123_v09_new_r1` | 変換器 node 周期・種 (new r1) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0124_v52_new_r1` | 変換器 node (CHT 用スラブ) (new r1) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0125_v44_new_r1` | 変換器 軸対称・種・凝縮 (new r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0126_c36node_new_r2` | 2D node 標準 (new r2) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0127_c36node_impdiag_new_r2` | 2D node + 陰解法診断 CSV (表 M) (new r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0128_c36node_psidual_new_r2` | 2D node + ψ 二重評価 (R2 の pdeSize) (new r2) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0129_c09ckpt100_new_r2` | dual-time 100 step (checkpoint 書出し) (new r2) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0130_c09restart100_new_r2` | dual-time 再開 100 step (checkpoint 復元) (new r2) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0131_c09cont200_new_r2` | dual-time 連続 200 step (new r2) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0132_c09ckpt100_outres_new_r2` | dual-time + 残差出力 (表 K・L) (new r2) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0133_c09ckpt100_rawdiag_new_r2` | dual-time + roYraw (表 I) (new r2) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0134_c09ckpt100_pindiag_new_r2` | dual-time + ピン診断 (表 N) (new r2) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0135_c44steady_new_r2` | 軸対称・多成分・凝縮 (定常) (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0136_c44dual_ckpt100_new_r2` | 軸対称・凝縮 dual-time 前半 (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0137_c44dual_restart100_new_r2` | 軸対称・凝縮 dual-time 後半 (再開) (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0138_c44dual_pindiag_new_r2` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0139_c52cht_new_r2` | 共役伝熱 (new r2) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0140_c36cell_new_r2` | cell 定常 SST (new r2) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0141_c20cell_rk3_new_r2` | cell RK3 (new r2) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0142_c20cell_dual_new_r2` | cell LES dual-time (new r2) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0143_c20cell_impdiag_new_r2` | cell 陰解法 + 診断 CSV (new r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0144_c57lm_new_r2` | 遷移 (roGamma を読む) (new r2) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0145_c57lm_fromsst_new_r2` | 遷移 (SST の場から初期化) (new r2) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0146_c56lineimp_new_r2` | line-implicit + extraFields (new r2) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0147_c56extra_new_r2` | extraFields (res_* を含む) + FP64 アキュムレータ (new r2) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0148_c48absorb_new_r2` | extraFields (roN・dq_block_old_*) (new r2) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0149_c26optin_new_r2` | env 診断なし (名前が消えて警告) (new r2) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0150_c26optin_env_new_r2` | env 診断 4 種 ON (表 I) (new r2) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0151_v36node_new_r2` | 変換器 node (new r2) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0152_v36cell_new_r2` | 変換器 cell (new r2) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0153_v09_new_r2` | 変換器 node 周期・種 (new r2) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0154_v52_new_r2` | 変換器 node (CHT 用スラブ) (new r2) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0155_v44_new_r2` | 変換器 軸対称・種・凝縮 (new r2) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0156_c36node_new_r3` | 2D node 標準 (new r3) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0157_c36node_impdiag_new_r3` | 2D node + 陰解法診断 CSV (表 M) (new r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0158_c36node_psidual_new_r3` | 2D node + ψ 二重評価 (R2 の pdeSize) (new r3) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0159_c09ckpt100_new_r3` | dual-time 100 step (checkpoint 書出し) (new r3) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0160_c09restart100_new_r3` | dual-time 再開 100 step (checkpoint 復元) (new r3) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0161_c09cont200_new_r3` | dual-time 連続 200 step (new r3) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0162_c09ckpt100_outres_new_r3` | dual-time + 残差出力 (表 K・L) (new r3) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0163_c09ckpt100_rawdiag_new_r3` | dual-time + roYraw (表 I) (new r3) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0164_c09ckpt100_pindiag_new_r3` | dual-time + ピン診断 (表 N) (new r3) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0165_c44steady_new_r3` | 軸対称・多成分・凝縮 (定常) (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0166_c44dual_ckpt100_new_r3` | 軸対称・凝縮 dual-time 前半 (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0167_c44dual_restart100_new_r3` | 軸対称・凝縮 dual-time 後半 (再開) (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0168_c44dual_pindiag_new_r3` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0169_c52cht_new_r3` | 共役伝熱 (new r3) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0170_c36cell_new_r3` | cell 定常 SST (new r3) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0171_c20cell_rk3_new_r3` | cell RK3 (new r3) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0172_c20cell_dual_new_r3` | cell LES dual-time (new r3) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0173_c20cell_impdiag_new_r3` | cell 陰解法 + 診断 CSV (new r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0174_c57lm_new_r3` | 遷移 (roGamma を読む) (new r3) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0175_c57lm_fromsst_new_r3` | 遷移 (SST の場から初期化) (new r3) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0176_c56lineimp_new_r3` | line-implicit + extraFields (new r3) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0177_c56extra_new_r3` | extraFields (res_* を含む) + FP64 アキュムレータ (new r3) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0178_c48absorb_new_r3` | extraFields (roN・dq_block_old_*) (new r3) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0179_c26optin_new_r3` | env 診断なし (名前が消えて警告) (new r3) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0180_c26optin_env_new_r3` | env 診断 4 種 ON (表 I) (new r3) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0181_v36node_new_r3` | 変換器 node (new r3) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0182_v36cell_new_r3` | 変換器 cell (new r3) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0183_v09_new_r3` | 変換器 node 周期・種 (new r3) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0184_v52_new_r3` | 変換器 node (CHT 用スラブ) (new r3) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0185_v44_new_r3` | 変換器 軸対称・種・凝縮 (new r3) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0186_c52cht_base_r4` | 共役伝熱 (base r4) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0187_c52cht_base_r5` | 共役伝熱 (base r5) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0188_c52cht_base_r6` | 共役伝熱 (base r6) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0189_c52cht_base_r7` | 共役伝熱 (base r7) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0190_c52cht_new_r4` | 共役伝熱 (new r4) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0191_c52cht_new_r5` | 共役伝熱 (new r5) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0192_c52cht_new_r6` | 共役伝熱 (new r6) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0193_c52cht_new_r7` | 共役伝熱 (new r7) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0194_fw_c44dual_b1` | 固定幅の独立 A/B (plan §6.3) b1: base (9c9f623c)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 4 量、`condClampCorrQ_0` 最大 4.62e18 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0195_fw_c44dual_n1` | 固定幅の独立 A/B (plan §6.3) n1: new (93e55957)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 3 量、`condClampCorrQ_0` 最大 7.75e20 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0196_fw_c44dual_b2` | 固定幅の独立 A/B (plan §6.3) b2: base (9c9f623c)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 2 量、`condClampCorrQ_0` 最大 **1.45e23** ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0197_fw_c44dual_n2` | 固定幅の独立 A/B (plan §6.3) n2: new (93e55957)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 0 量、`condClampCorrQ_0` 最大 2.44e17 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0198_fw_c44dual_b3` | 固定幅の独立 A/B (plan §6.3) b3: base (9c9f623c)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 0 量、`condClampCorrQ_0` 最大 8.14e17 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0199_fw_c44dual_n3` | 固定幅の独立 A/B (plan §6.3) n3: new (93e55957)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 7 量、`condClampCorrQ_0` 最大 **1.45e23** ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0200_fw_c44dual_b4` | 固定幅の独立 A/B (plan §6.3) b4: base (9c9f623c)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 4 量、`condClampCorrQ_0` 最大 **1.46e23** ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0201_fw_c44dual_n4` | 固定幅の独立 A/B (plan §6.3) n4: new (93e55957)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 2 量、`condClampCorrQ_0` 最大 7.75e20 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0202_fw_c44dual_b5` | 固定幅の独立 A/B (plan §6.3) b5: base (9c9f623c)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 4 量、`condClampCorrQ_0` 最大 7.75e20 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0203_fw_c44dual_n5` | 固定幅の独立 A/B (plan §6.3) n5: new (93e55957)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 4 量、`condClampCorrQ_0` 最大 **1.46e23** ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0204_fw_c44dual_b6` | 固定幅の独立 A/B (plan §6.3) b6: base (9c9f623c)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 0 量、`condClampCorrQ_0` 最大 1.32e17 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
| `run_0205_fw_c44dual_n6` | 固定幅の独立 A/B (plan §6.3) n6: new (93e55957)、c44dual_ckpt100 と同じ入力・100 step + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0、NaN PASS、幅超過 0 量、`condClampCorrQ_0` 最大 5.46e16 ([RESULT.txt](fixedwidth_c44dual_ckpt100/RESULT.txt)) | active |
