# 諮問ブリーフ: R2c (実 GPU の延長経路) の解釈と、同じ型の欠陥の残り (R2h) の切り分け

- 日付: 2026-10-12、ブランチ `feature/sern-design`、commit `9081e272`
- 関連 plan: [`plans/accepted/tooling-sern-te-wake-grid.md`](../../../plans/accepted/tooling-sern-te-wake-grid.md) §5.1 R2〜R2h・§6 #11
- 前回の諮問: [`notes/reviews/2026-10-12-sern-judgment-path-diagnose.md`](../2026-10-12-sern-judgment-path-diagnose.md) (M1〜M3・Minor を全件採用 → R2a〜R2e)
- 諮る理由: AGENTS.md 条件 7 (result 段の解釈を確定する前)。R2c の担当は「解釈は F」、R2h は「切り分けは F」。

## 1. 観測事実

### R2c: `case/46.sern_design/run_1100_r2c_lsw10_m6_judgment/` (AWS-B、commit 15968c33、写し `notes/investigations/2026-10-12-r2c/`)
- 手順は §6 #11 の事前登録どおり: 生産 YAML `problem_3d_prod_3op_wallres_lswx10.yaml` の m6_on で prepare → 格子署名 5f22ff6e…・双対ハッシュ 475ce484… が run_1094 と一致、本段 config・bcondConfig・species_meta の差なし → `run_1094_tewake_lsw10_m6/res_4000.h5` から `runner_sern.run_judgment` (N 20000・W 10000・延長 1 回、強制の分岐なし)。`FORGE_CUDA_BLOCKSIZE=128`、`FORGE_ALLOW_UNVERIFIED_SPECIES` は未設定 (起点の場は印付き)。
- 1 回目 (`res_4000.h5` から): rc 0、FAIL/JUDGMENT_UNSETTLED。床 PASS (10000, 20000]・更新 10000・事象 0。窓 (a / d、許容): C_T 3.70e-6 / 2.65e-5 OK、**C_T_with_shear 3.32e-6 / 6.92e-5 NG (5e-5)**、C_L 1.96e-5 / 4.37e-5 OK、**C_M 5.23e-4 / 1.92e-3 NG (5e-4)**。温度: R_TE 220.8 K・R_SE 224.7 K・RET 220.6 K は STEADY、**ALL 210.0 K は TRANSIENT-UNSETTLED**、T < 150 K は全領域 0。延長の理由「未定常 (窓・定常性の未達) だけ」。
- 2 回目 (`res_20000.h5` から自動延長): rc 0、**PASS**。床 PASS (10000, 20000]・更新 10000・事象 0。窓: C_T 2.80e-6 / 3.98e-7、C_T_with_shear 2.67e-6 / 3.70e-7、C_L 1.15e-5 / 8.28e-7、C_M 2.55e-4 / 2.38e-5 (すべて許容内・STEADY)。温度 4 領域 STEADY (ALL 209.5 K)、T < 150 K 0。
- 手作業の検証 run との比較 (同じ起点): 1 回目 ↔ run_1097 (C_T_with_shear の d 6.85e-5・C_M の d 1.97e-3 で NG、同じ形)。2 回目 ↔ run_1099 (手作業の +20000、`restart_field.py`): 平均 C_T 0.8926329 / 0.8926327、C_T_with_shear 0.8725507 / 0.8725505、C_L 0.07148224 / 0.0714822、C_M −1.848920 / −1.848921。
- 残差 (`check_convergence.py`、判定区間の各 session): 両回とも NOT CONVERGED (全列 plateau または低下中、rising なし。plan §6 のプラトー受理の方針どおりで「収束」とは呼ばない)。`--from-floor run_1094` は REFUSED (参照 run が通常基準で PASS しないため。情報のみ)。
- `metrics.json` (正本) は 2 回目の判定と一致、1 回目の記録は `*_judgment1.*` に退避。`JUDGMENT.json` に env・plan・required_state。

### R2f・R2g (実装済み、commit 1595e86b・9081e272)
- R2f: 受理ゲートの目的量・量を `evaluate.model` から (SST で列が欠けたら判定不能)。保存済み 3D 12 本・case/46 の 10 本で判定が前後一致 (case/46 の metrics.json 69 本の全走査でも同一)。
- R2g: `field_health`・`floor_gate`・`check_convergence.py` の必須量を乱流モデルから (SST に roK・roOmega)、モデル名は `probdef.sern_turbulence_model` で小文字の sst・euler だけ受理。case/46 の 144 本で判定が前後一致、`check_convergence.py` は CLI 14 本と主チェックアウトの 2866 run の `analyze` の戻りが前後一致。試験はすべて PASS (`run_sern_required_state_tests.py` 102、`test_gate_bad_input.py` 71 ほか 14 本)。

### R2h: 同じ型の欠陥のさらに残り (R2g の実装担当の洗い出し、未修正)
plan §5.1 R2h の (1)〜(8)。要約:
1. `sern_gates.residual_health` が check_convergence の判定不能 (必須列の欠落・末尾窓の代表値 0) を受理ゲートに反映しない (R2g の SST 列の検査が受理に効かない)。
2. 化学種 roY* を h5 にあるときだけ検査 (`field_health`・`floor_gate` の ΣY)。
3. `sern_forces.py:81`・`runner_sern3d.py:398`: ある step で壁の部品の出力が無いと、黙ってその部品抜きの C_T になる。
4. `sern_forces.py:92`・`runner_sern3d.py:432`: twall の無い部品の摩擦を 0 として足し、C_T_with_shear が部分和でも印が無い。
5. 初期場・restart で roK・roOmega をあるときだけ写す (判定区間の段の移行だけは R2d が検査)。
6. `check_quasisteady.py:418`: 適用できない量を黙って落とし、0 量で ALL STEADY (`--quantity asym` を非対称格子の run_0198 に当てて実測)。`:447` `wall_model_tau` も utau の無い壁を飛ばす。
7. `check_floor_events.py:149`: res_*.h5 が無いと最終 step を照合しない。
8. `te_monitor.py:51` の L_cowl 既定 1.2 H (runner の形状照合で拾う)、`:74・:86` 空の領域の行で DictWriter の列がずれる。

## 2. 期待値と出典
- §6 #11 (事前登録 2026-10-12): 合格 = 1 回目が「未定常だけ」で不合格 → 自動で延長 → 2 回目が #1〜#3 PASS・GATES PASS。記録: 2 回目の床の session と区間・#2・#3・残差 VERDICT・`*_judgment1.*` の退避・最終の metrics と JUDGMENT の対応。
- §6 の事前登録: 各 run 20000、窓未達なら +20000 を 1 回、末尾 10000 step で判定、許容 5e-5 (C_T・C_T_with_shear・C_L)・5e-4 (C_M)。

## 3. 問い
- Q1. R2c は §6 #11 の合格とみなしてよいか。1 回目の不合格の理由に温度 (ALL の Tmin が TRANSIENT-UNSETTLED) が含まれていたが、延長の条件 (未達が「未定常」だけ) の中で扱った。これは事前登録の「未定常だけ」に含まれるか。
- Q2. R2h の (1)〜(8) のうち、この plan (生産受入れ) の status done の前に直すべきものはどれか。残りの行き先 (親 plan R7b・新規 plan・共有ツールの plan など) をどう分けるか。実装担当は (1) を最優先と見ている。(3)(4) は C_T・C_T_with_shear の値そのものに効きうる。(6)(7) は全セッションが使う共有ツール。
- Q3. 上を反映した後、codex の `--stage result` の再レビューに回してよいか。

読んでよいファイル: plan 本体、`notes/investigations/2026-10-12-r2c/` (R2C_SUMMARY.txt・JUDGMENT.json・R2C_CONVERGENCE_*・R2C_QUASISTEADY_*)、`notes/investigations/2026-10-11-lsw10/run_1097_*`・`run_1099_*`、`design/forge_design/metrics/sern_gates.py`・`sern_forces.py`、`design/forge_design/evaluate/runner_sern.py`・`runner_sern3d.py`、`solver_density_cuda/tools/check_quasisteady.py`・`check_floor_events.py`・`check_convergence.py`、`case/46.sern_design/diag/te_monitor.py`。
