# case/66 — ホストメモリ削減の回帰ハーネス

## ケース概要

plan [architecture-solver-host-memory](../../plans/active/architecture-solver-host-memory.md) §5.1 #3 (基準入力の確定と回帰ハーネス) の実体。
§6 の構成表の各構成を、変更前ビルド (base) と変更後ビルド (new) で**同じ入力・同じ手順で 3 回ずつ**回し、
§6 の (a) step 0 のビット一致・(b) N step 後の差が同ビルド内ばらつきの 2 倍以内・(c) NaN/Inf 無し・出力互換性 (データセット集合・shape・dtype・属性) を
`compare_runs.py` で判定する。

- **base** = `9c9f623c` (AWS `~/bin-hostmem/forge_9c9f623c`、sha256 `c65dbb39f772…`、変換器 `~/bin-hostmem/convert_9c9f623c`)。
- **run の置き場所は AWS** `~/forge-b4/case/66.hostmem_regression/run_NNNN_<構成>_<ビルド>_r<k>/`。比較も AWS で行い、手元にはテキストだけを持ち帰る
  (ローカルでは forge を回さない。重い h5 を手元へ引かない)。このディレクトリ (ローカル) には入力 config・スクリプト・README だけを置く
  (h5 は `.gitignore` で追跡外)。
- 元入力は別セッションのワークツリー `/home/sano/work/forge/case/...` にあり、**読むだけ**で複製して使う (`prepare_inputs.py`)。
- SERN 3D (g3) は `case/46.sern_design` の `run_1072`〜`run_1074` (AWS `~/forge-r8/case/46.sern_design/`、[case/46 README](../46.sern_design/README.md) の run 一覧)。

## ファイル

| ファイル | 役割 |
| --- | --- |
| `matrix_spec.py` | **構成表の正本** (入力の元 run・複製するファイル・設定の修正・種・checkpoint・構成ごとの環境変数) |
| `prepare_inputs.py` | (ローカル) 元 run から `inputs/<入力名>/` を作る。修正は正規表現 + 期待一致数で行い、PyYAML で期待値を検査。`mesh.bndFirstOrder`・`wallTreatmentSST: 1` が残れば失敗 |
| `run_matrix.py` | (AWS) `seed` (restart_field.py を種に掛ける)・`set-ckpt`・`launch` (run を作って 1 本ずつ順に回すワーカーを裏で起動)・`resume`・`status`・`verify` (RUN_PROVENANCE の forge_bin/sha256 照合)・`note` (README の状態・比較から除外)・`table` (下の run 一覧の行) |
| `compare_runs.py` | (AWS) §6 の判定。`--cfg X` / `--all` (registry.tsv から run を選ぶ)、`--base … --new …` (直接指定)、`--diff2 A.h5 B.h5` (分割と連続など 2 ファイル) |
| `memlog_summary.py` | (AWS) `FORGE_MEMLOG=1` の工程別 RSS/HWM/GPU と `--memwatch` の 1 s 採取 (`mem_samples.csv`) を表にする |
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

## 既知の注意

- **AWS の自動停止 (2026-10-07 01:05 UTC 頃 = JST 10:05、base 投入中に停止した)**: idle 自動停止 (`solver_density_cuda/tools/cloud/idle_autostop.sh`、5 分おき・6 回連続で停止) は
  「GPU 使用率 0・`pgrep -x forge` 0 件・ログイン 0・load < 1」を idle と数える。base は `forge_9c9f623c` という名前で起動していたので
  **`pgrep -x forge` に一致せず**、短い 2D run の合間 (GPU 使用率の瞬間値 0) と準備作業の時間で 30 分が数えられたと推定する (停止の記録は見ていない)。
  対策として `run_matrix.py` は `--bin` を `.bin/<ビルド名>/forge` というシンボリックリンク経由で起動するようにした (comm が `forge` になる)。
  **この対策は AWS 停止後に入れたので未実行**。次回投入時に `pgrep -x forge` で見えることを確認すること。
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
- **base の 3 回がそろっていない** (下の run 一覧): r1 は全構成完走 (NaN PASS)、r2/r3 と再開 r1 は投入直後に AWS が止まったので未確認。
  同ビルド内ばらつき (§6 (b) の基準 S_base) はまだ出ていない。

## 計算 run 一覧

AWS `~/forge-b4/case/66.hostmem_regression/` (2026-10-07)。全部 base (`~/bin-hostmem/forge_9c9f623c`、`FORGE_CUDA_BLOCKSIZE=128`)。
N step の短い回帰 run なので収束判定は NOT CONVERGED (回帰差の基準であって収束の主張には使わない)。

| run | 目的・設定差分 | 元の入力 | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0001_c09ckpt100_base_r1` | dual-time 100 step (checkpoint 書出し) (base r1) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0002_c44dual_ckpt100_base_r1` | 軸対称・凝縮 dual-time 前半 (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0003_c36node_base_r1` | 2D node 標準 (base r1) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0004_c36node_impdiag_base_r1` | 2D node + 陰解法診断 CSV (表 M) (base r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0005_c36node_psidual_base_r1` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r1) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0006_c09cont200_base_r1` | dual-time 連続 200 step (base r1) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0007_c09ckpt100_outres_base_r1` | dual-time + 残差出力 (表 K・L) (base r1) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0008_c09ckpt100_rawdiag_base_r1` | dual-time + roYraw (表 I) (base r1) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0009_c09ckpt100_pindiag_base_r1` | dual-time + ピン診断 (表 N) (base r1) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0010_c44steady_base_r1` | 軸対称・多成分・凝縮 (定常) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0011_c44dual_pindiag_base_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0012_c52cht_base_r1` | 共役伝熱 (base r1) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0013_c36cell_base_r1` | cell 定常 SST (base r1) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0014_c20cell_rk3_base_r1` | cell RK3 (base r1) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0015_c20cell_dual_base_r1` | cell LES dual-time (base r1) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0016_c20cell_impdiag_base_r1` | cell 陰解法 + 診断 CSV (base r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0017_c57lm_base_r1` | 遷移 (roGamma を読む) (base r1) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0018_c57lm_fromsst_base_r1` | 遷移 (SST の場から初期化) (base r1) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0019_c56lineimp_base_r1` | line-implicit + extraFields (base r1) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0020_c56extra_base_r1` | extraFields (res_* を含む) + FP64 アキュムレータ (base r1) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0021_c48absorb_base_r1` | extraFields (roN・dq_block_old_*) (base r1) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0022_c26optin_base_r1` | env 診断なし (名前が消えて警告) (base r1) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0023_c26optin_env_base_r1` | env 診断 4 種 ON (表 I) (base r1) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0024_v36node_base_r1` | 変換器 node (base r1) | `inputs/v36node` ← `36.passive_pseudoshock_control` | 完了 (rc 1 = 終了時 cudaFree の GPUassert、既知の罠。h5 は完全)・NaN PASS | active |
| `run_0025_v36cell_base_r1` | 変換器 cell (base r1) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | 完了 (rc 1 = 終了時 cudaFree の GPUassert、既知の罠。h5 は完全)・NaN PASS | active |
| `run_0026_v09_base_r1` | 変換器 node 周期・種 (base r1) | `inputs/v09` ← `09.Taylor-Green` | 完了 (rc 1 = 終了時 cudaFree の GPUassert、既知の罠。h5 は完全)・NaN PASS | active |
| `run_0027_v52_base_r1` | 変換器 node (CHT 用スラブ) (base r1) | `inputs/v52` ← `52.conjugate_slab` | 完了 (rc 1 = 終了時 cudaFree の GPUassert、既知の罠。h5 は完全)・NaN PASS | active |
| `run_0028_v44_base_r1` | 変換器 軸対称・種・凝縮 (base r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0029_c44dual_ckpt100_base_r1` | 軸対称・凝縮 dual-time 前半 (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0030_c44steady_base_r1` | 軸対称・多成分・凝縮 (定常) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0031_c44dual_pindiag_base_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 完走 rc 0・NaN PASS・NOT CONVERGED (短い回帰 run、収束の主張には使わない) | active |
| `run_0032_v44_base_r1` | 変換器 軸対称・種・凝縮 (base r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 完了 (rc 1 = 終了時 cudaFree の GPUassert、既知の罠。h5 は完全)・NaN PASS | active |
| `run_0033_c09restart100_base_r1` | dual-time 再開 100 step (checkpoint 復元) (base r1) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | 2026-10-07 01:04 投入直後に AWS が停止 (01:05 UTC 頃、自動停止と推定)。完了したか未確認 | 要確認 |
| `run_0034_c44dual_restart100_base_r1` | 軸対称・凝縮 dual-time 後半 (再開) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 2026-10-07 01:04 投入直後に AWS が停止 (01:05 UTC 頃、自動停止と推定)。完了したか未確認 | 要確認 |
| `run_0035_c09ckpt100_base_r2` | dual-time 100 step (checkpoint 書出し) (base r2) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0036_c09restart100_base_r2` | dual-time 再開 100 step (checkpoint 復元) (base r2) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0037_c09cont200_base_r2` | dual-time 連続 200 step (base r2) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0038_c44dual_ckpt100_base_r2` | 軸対称・凝縮 dual-time 前半 (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0039_c44dual_restart100_base_r2` | 軸対称・凝縮 dual-time 後半 (再開) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0040_c36node_base_r2` | 2D node 標準 (base r2) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0041_c36node_impdiag_base_r2` | 2D node + 陰解法診断 CSV (表 M) (base r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0042_c36node_psidual_base_r2` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r2) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0043_c09ckpt100_outres_base_r2` | dual-time + 残差出力 (表 K・L) (base r2) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0044_c09ckpt100_rawdiag_base_r2` | dual-time + roYraw (表 I) (base r2) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0045_c09ckpt100_pindiag_base_r2` | dual-time + ピン診断 (表 N) (base r2) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0046_c44steady_base_r2` | 軸対称・多成分・凝縮 (定常) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0047_c44dual_pindiag_base_r2` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0048_c52cht_base_r2` | 共役伝熱 (base r2) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0049_c36cell_base_r2` | cell 定常 SST (base r2) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0050_c20cell_rk3_base_r2` | cell RK3 (base r2) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0051_c20cell_dual_base_r2` | cell LES dual-time (base r2) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0052_c20cell_impdiag_base_r2` | cell 陰解法 + 診断 CSV (base r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0053_c57lm_base_r2` | 遷移 (roGamma を読む) (base r2) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0054_c57lm_fromsst_base_r2` | 遷移 (SST の場から初期化) (base r2) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0055_c56lineimp_base_r2` | line-implicit + extraFields (base r2) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0056_c56extra_base_r2` | extraFields (res_* を含む) + FP64 アキュムレータ (base r2) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0057_c48absorb_base_r2` | extraFields (roN・dq_block_old_*) (base r2) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0058_c26optin_base_r2` | env 診断なし (名前が消えて警告) (base r2) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0059_c26optin_env_base_r2` | env 診断 4 種 ON (表 I) (base r2) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0060_v36node_base_r2` | 変換器 node (base r2) | `inputs/v36node` ← `36.passive_pseudoshock_control` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0061_v36cell_base_r2` | 変換器 cell (base r2) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0062_v09_base_r2` | 変換器 node 周期・種 (base r2) | `inputs/v09` ← `09.Taylor-Green` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0063_v52_base_r2` | 変換器 node (CHT 用スラブ) (base r2) | `inputs/v52` ← `52.conjugate_slab` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0064_v44_base_r2` | 変換器 軸対称・種・凝縮 (base r2) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0065_c09ckpt100_base_r3` | dual-time 100 step (checkpoint 書出し) (base r3) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0066_c09restart100_base_r3` | dual-time 再開 100 step (checkpoint 復元) (base r3) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0067_c09cont200_base_r3` | dual-time 連続 200 step (base r3) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0068_c44dual_ckpt100_base_r3` | 軸対称・凝縮 dual-time 前半 (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0069_c44dual_restart100_base_r3` | 軸対称・凝縮 dual-time 後半 (再開) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0070_c36node_base_r3` | 2D node 標準 (base r3) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0071_c36node_impdiag_base_r3` | 2D node + 陰解法診断 CSV (表 M) (base r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0072_c36node_psidual_base_r3` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r3) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0073_c09ckpt100_outres_base_r3` | dual-time + 残差出力 (表 K・L) (base r3) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0074_c09ckpt100_rawdiag_base_r3` | dual-time + roYraw (表 I) (base r3) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0075_c09ckpt100_pindiag_base_r3` | dual-time + ピン診断 (表 N) (base r3) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0076_c44steady_base_r3` | 軸対称・多成分・凝縮 (定常) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0077_c44dual_pindiag_base_r3` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0078_c52cht_base_r3` | 共役伝熱 (base r3) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0079_c36cell_base_r3` | cell 定常 SST (base r3) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0080_c20cell_rk3_base_r3` | cell RK3 (base r3) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0081_c20cell_dual_base_r3` | cell LES dual-time (base r3) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0082_c20cell_impdiag_base_r3` | cell 陰解法 + 診断 CSV (base r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0083_c57lm_base_r3` | 遷移 (roGamma を読む) (base r3) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0084_c57lm_fromsst_base_r3` | 遷移 (SST の場から初期化) (base r3) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0085_c56lineimp_base_r3` | line-implicit + extraFields (base r3) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0086_c56extra_base_r3` | extraFields (res_* を含む) + FP64 アキュムレータ (base r3) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0087_c48absorb_base_r3` | extraFields (roN・dq_block_old_*) (base r3) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0088_c26optin_base_r3` | env 診断なし (名前が消えて警告) (base r3) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0089_c26optin_env_base_r3` | env 診断 4 種 ON (表 I) (base r3) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0090_v36node_base_r3` | 変換器 node (base r3) | `inputs/v36node` ← `36.passive_pseudoshock_control` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0091_v36cell_base_r3` | 変換器 cell (base r3) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0092_v09_base_r3` | 変換器 node 周期・種 (base r3) | `inputs/v09` ← `09.Taylor-Green` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0093_v52_base_r3` | 変換器 node (CHT 用スラブ) (base r3) | `inputs/v52` ← `52.conjugate_slab` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
| `run_0094_v44_base_r3` | 変換器 軸対称・種・凝縮 (base r3) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 投入済み・AWS 停止で中断 (未完了の可能性。`status` で running のまま/queued なら除外して再投入) | 要確認 |
