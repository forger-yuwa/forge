# case/68 — ③ベル 2 モード比較 (TP のロケット排気)

plan [`tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md) の検証・比較 case。
③ベルの壁を (M) 特性曲線で辿る Rao (TP) と (D) 形状を直接変える C2 B-spline の 2 つの設計手順で作り、同じ評価器で比べる。
親は [`tooling-nozzle-design-tool.md`](../../plans/active/tooling-nozzle-design-tool.md) §3 確定事項 8〜11。

## 基準ケース

- ガス: LOX/LH2・燃焼室圧 10 MPa・O/F 6 の**燃焼室平衡組成を凍結した TP** (NASA CEA2 本体で計算)。
  CFD は単一の lump 擬似種 `EXH` (`thermalMethod: 2`・`thermoHrefTemp: 298.15`・`viscMethod: 2` + `physProp.transport`)。
  - 反応物: H2(L) 20.27 K、O2(L) 90.17 K。燃焼器は無限面積 (CEA `rocket` の既定)。凍結位置 = 燃焼室 (`rocket frozen nfz=1`)。
  - 生成物の候補: H2・O2・H2O・OH・H・O (forge の熱物性と輸送データの両方にある種)。全生成物で解いた場合に出る HO2・H2O2・O3 は
    モル分率の和 4.5e-5、燃焼室温度の差 −0.17 K (`exhaust_lox_lh2_pc100_of6.yaml` の `omitted_species`)。
  - 燃焼室温度 3523.96 K、モル分率 H2O 0.67299 / H2 0.24799 / OH 0.042114 / H 0.031081 / O2 0.002975 / O 0.002854 (CEA の印字 5 桁)。
- 面積比 ε 50、真空の設計点 (出口は `outflow`)、断熱壁、SST 低 Re (`wallTreatmentSST: 0`)。スロート半径 0.14 m (plan §8-1 の仮値)。

## ファイル

| ファイル | 内容 |
| --- | --- |
| `exhaust_lox_lh2_pc100_of6.yaml` | 排気ファイル (schema `forge_rocket_exhaust_v1`)。組成・燃焼室温度・燃焼室圧の**単一ソース**と、同じ条件の CEA 凍結膨張の表 (ε 10/25/50/100、1D 理想値の照合の参照)。手で編集しない |
| `problem_bell_tp_top.yaml` | 評価器 (TP 配線) の確認用の問題定義。形状は既存の TOP を ε 50 に合わせた仮の 1 点 (θn 32°・θe 9°・L 18 r_t)。設計の比較ではない |
| `problem_bell_tp_top_2b.yaml` | plan §5.1 #2b (推力の CV 会計とゲート) の確認用。形状・ガス・格子は上と同じで、段階起動 (`evaluate.stages`・`stage_steps`・本段の `nStepOuter`/`outStepInterval`) を足したもの。`python -m forge_design.evaluate.runner <YAML> <run> --staged` (継続は `--continue-from <run>`) |
| `problem_bell_tp_top_2b_wallfine.yaml` | 同上で第一セル厚をスロート近傍だけ詰めた格子 (`wall_first_frac_throat` と前後のブレンド、1201×121、AR ≤ 1000)。引き継ぎ (c) の確認用 |

排気ファイルの作り直し (CEA2 は `FORGE_CEA_DIR`、無ければ主作業ツリーの `.venv-cea/nasa_cea`):

```bash
cd design
python -m forge_design.gas.rocket_exhaust --supar 10,25,50,100 --out ../case/68.bell_two_mode/exhaust_lox_lh2_pc100_of6.yaml
```

1D 理想値 (同じ凍結 TP の等エントロピー膨張) と CEA 凍結計算の照合は `design/tests/run_bell_tp_tests.py`
(ε 50 で出口温度・ṁ_1D・出口速度・出口圧力・F_1D,vac とも相対 2.2e-5 以内)。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_tp_soft_smoke` | plan §5.1 #2a: TP 配線の起動確認。`problem_bell_tp_top.yaml`、soft 段 (`convMethod 0`・`limiter 0`・`cfl_pseudo 0.5`・`nStepInner 10`・`implicitRelax 0.7`) 1000 step、準 1D 等エントロピー IC (冷間)。メッシュ 561×121 (`wall_first_frac 2.5e-5`、品質 PASS: AR 最大 836・skew 0.417)。forge は本作業ツリーの native ビルド (`FORGE_BIN`) | NaN/Inf なし (`residual_history.csv` 全 rms 列、`res_0/250/500/750/1000.h5` の VALUE 全変数)。IC の化学種属性 `verified`。`CONVERGENCE_VERDICT.txt` は NOT CONVERGED (soft 1000 step、収束は要求外)。過渡で縮小部 (x/r_t −4.0〜−0.6、ノードの 6.9 %) の P・T が Pc・Tc を超える (res_1000 で最大 1.068e7 Pa・3567 K、res_250 の 1.159e7 Pa・3622 K から減少中)。2.9 ms/step。`residual_history.png` | active |
| `run_0002_tp_staged_cv` | plan §5.1 #2b: 推力の CV 会計とゲートの全流れ解。`problem_bell_tp_top_2b.yaml` (形状・ガス・格子 561×121 は run_0001 と同じ、品質 PASS: AR 836・skew 0.417)。`runner --staged`: soft 2000 → mid 2000 → 本段 12000 (`convMethod 1`・`limiter 2`・`cfl_pseudo 4`・`implicitRelax 0.7`・`nStepInner 4`) を同一 run・同一メッシュで、段の間は `restart_field.py`、`stage_manifest.json` に段を記録。出力に勾配の `extraFields` | 本段 12000 では未整定: CV 残差 −6.6e-4 (step 1000 の −1.4e-2 から単調に減少中)、ṁ 収支 +9.3e-4。`CONVERGENCE_SEGMENT_VERDICT.txt` (区間 S3_main): NOT CONVERGED (stalled)。NaN/Inf なし。`thrust_history.csv`・`GATES.json`・`residual_history.png` | active (run_0003 の継続元) |
| `run_0003_tp_cv_continue` | run_0002 の本段を同一メッシュで 24000 step 継続 (`runner --continue-from`: `restart_field.py`、段の履歴と manifest を引き継ぎ、判定区間は S3_main → S4_main_p の本段 36000 step) | 壁会計 1.135026e6 N・出口面会計 1.135027e6 N、**CV 残差 −5.2e-7** (F_1D 比)、η (出口面会計) 0.97922、ṁ/ṁ_1D 0.99255、ṁ 収支 +5.4e-5、内壁の摩擦 −1.18e4 N。`CONVERGENCE_SEGMENT_VERDICT.txt`: NOT CONVERGED (stalled; rms_ro 2.5・rms_roUy 2.2・rms_roe 2.5・rms_roOmega 1.6 桁で横ばい)。`QUASISTEADY_THRUST_VERDICT.txt`: 両推力会計・η・ṁ STEADY、cv_residual DRIFTING (0 へ近づく系列の相対ドリフト)。局所 y1+ 最大 5.68 (x/r_t −0.31)、y1+>1 が点の 55.6 %・面積の 25.7 % (`WALL_RESOLUTION_VERDICT.txt` FAIL、`wall_yplus_profile.csv`)。収縮部 (x/r_t −4.0〜−0.6) の P は全点 < Pc、T > Tc が 127 ノード (最大 +0.35 K、壁から 48 µm 以内)。NaN/Inf なし。2.25 ms/step | active |
| `run_0004_tp_staged_wallfine` | 引き継ぎ (c): `problem_bell_tp_top_2b_wallfine.yaml` — 1201×121、`wall_first_frac_throat 7e-6` (スロートで 0.97 µm)・前後のブレンド・`throat_refine 8`・`throat_width 3` (品質 PASS: AR 最大 979.5・skew 0.417)。soft 2000 → mid 2000 → 本段 36000 | CV 残差 −7.3e-5 (減少中)、η 0.97921、ṁ/ṁ_1D 0.99264。区間 S3_main: NOT CONVERGED (stalled; rms_roUy 2.7・rms_roOmega 1.7 桁、他は 3 桁以上)。推力・η・ṁ STEADY、cv_residual DRIFTING。局所 y1+ 最大 2.55 (入口直管の始端 x/r_t −4.13)、スロート (\|x/r_t\| ≤ 0.35) 最大 1.60、y1+>1 が点の 26.0 %・面積の 6.2 % (FAIL)。収縮部の P は全点 < Pc、T > Tc が 1095 ノード (最大 +0.51 K、壁から 57 µm 以内)。NaN/Inf なし。3.74 ms/step | active |
| `v0_moc_gate/` V0-a (CFD なし) | plan §6 V0: 前進 MOC ゲート (`geometry/moc_forward.py`) の A/B。(D) の格子 135 点 (θ_max {25,32.5,40}°×θ_e {5,10,15}°×x_m/L {0.1,0.2,0.3}×L_k 60〜100 %) の成立判定と、成立 8 点・TOP (θn 25.8°・θe 10°・R_u=R_d=1) 5 点のゲート。ガス EXH 凍結 TP、初期線 Hall (R_t 1)、主設定 n 321 (確認 n 161・641、参考の整合初期線 n 321/641/1281) | (D) 成立 8/135・通過 3/8 → **判定不能 (標本の退化)** (diagnostician 2026-10-11、plan §6 V0-a)。TOP 0/5 (C⁻ 族の交差、始まり x 3.9〜5.0 r_t)。Hall n 641 はスロート直後の初期線の不整合で全点不受理 (アーチファクト)。出力は `v0a_*` に改名: `v0a_summary.json`・`v0a_gate.csv`・`v0a_d_feasibility.csv`・`v0a_run_v0.log`、図 (git 管理外) `v0a_gate_map.png`・`v0a_nets.png`。スクリプト `run_v0.py` (§5.1 #3c のガード後は n 641 が拒否されるので再実行は commit `c636e716` で) | ref (V0-b に置き換え) |
| `v0_moc_gate/` V0-b (CFD なし) | plan §6 V0 の回し直し (§5.1 #3c): (D) の x_m を内部決定 (`design_bell_c2`: 0.5〜4.0 r_t を 141 点で走査 → 成立区間で ∫(r''')² 最小を黄金分割)、ゲートのガード 3 点 (Hall R_t ≤ 1 で n > 321 を拒否、前線 1 の反転 = 判定不能 (初期線不整合)、根拠が j ≤ 2 だけ = 判定不能 (軸近傍))。格子 θ_max {20,25,32.5,40}°×θ_e {5,10,15}°×L_k 60〜100 % = 60 点と TOP 5 点。ガス・ε・寸法・初期線は V0-a と同じ、主設定 Hall n 321、確認 n 161 | (D) 成立 31/60 (L60 7・L70 7・L80 6・L90 5・L100 6) → 有効条件を満たす。n 321: 通過 18・不受理 12・判定不能 1 (軸近傍) → 通過率 60.0 % (成立点を分母にすると 58.1 %) → **中間**。marginal 1 点 (L90・32.5°・5°、C⁻ 間隔比 0.056)。健全性の確認は**不一致** (40° 不受理 6/6 は一致、L60・32.5° は 3 点とも通過、32.5°・L ≥ 70 % は 9 点中 4 点が不受理) → 事前登録の文言で「ゲートが帯内の x_m に敏感」。TOP 0/5 (交差の始まり x 3.90〜5.04 r_t・r 0.86〜0.96)。n 161: 通過 20・不受理 11 → 中間、2 点で合否が違う (L100・20°・10° 不受理→通過、L100・25°・10° 判定不能→通過)。走査範囲: 成立する x_m は 60 点とも 1.08〜2.75 r_t (0.01 刻み) で 0.5〜4.0 の外に無し。補足 (`v0b_sliver_check.py`): 成立の判定 (2000 点の標本) が細い区間に依存 — 0.001 刻みでは 10 点が新たに成立、成立 31 点のうち 2 点 (L60・20°・10°、L70・25°・5°) は x_m の下流で θ が θ_max を超える (30.18°・26.55°)。成果物 `v0b_summary.json`・`v0b_gate.csv`・`v0b_d_design.csv`・`v0b_xm_scan_coverage.csv`・`v0b_sliver_check.csv`・`v0b_run.log`・`v0b_sliver_check.log`、図 (git 管理外) `v0b_gate_map.png`。スクリプト `run_v0b.py`・`v0b_sliver_check.py` | ref (V0-c に置き換え) |
| `v0_moc_gate/` V0-c (CFD なし) | plan §6 V0-c (事前登録、§5.1 #3d の後の記録の差し替え): (D) の成立条件 (ii) を「r'' ≥ −tol on [0, x_m] かつ r'' ≤ tol on [x_m, L]」(tol = `RESID_TOL` 1e-10) に置き換え、選んだ成立区間の端を走査点と隣の不成立点の間で二分法 (1e-6 r_t) で詰めてから黄金分割 (J 最小の規則は不変)。格子・TOP・ガス・初期線・閾値は V0-b と同じ (主 Hall n 321、確認 n 161) | (D) 成立 29/60 (L60 6・L70 6・L80 6・L90 5・L100 6) → 有効条件を満たす。V0-b の接触型 2 点 (L60・20°・10°、L70・25°・5°) が不成立になり、残り 29 点は成立のまま (増減なし)。x_m 1.109〜2.354 r_t、29 点中 18 点が詰めた端で確定 (V0-b からの x_m の差 −0.010〜+0.023 r_t)。n 321: 通過 15・不受理 13・判定不能 1 (軸近傍、L100・25°・10°) → 通過率 53.6 % (成立点を分母にすると 51.7 %) → **中間**、marginal 0 点。n 161: 通過 17・不受理 12 → 中間、主と違う点は 2 (V0-b と同じ L100・20°・10° 不受理→通過、L100・25°・10° 判定不能→通過)。V0-b の既存 29 点との合否の反転 1 点 (L90・32.5°・5° 通過→不受理: C⁻ 間隔比の最小 0.056 → −0.008、交差開始 x 20.35 r_t = x/L 0.998)、他の 28 点は verdict 同じ。健全性の確認は不一致 (40° 不受理 6/6 は一致、L60・32.5° は 3 点とも通過、32.5°・L ≥ 70 % は 9 点中 5 点が不受理)。TOP 0/5 (V0-b と同じ verdict)。走査範囲: 成立する x_m は 1.08〜2.75 r_t (0.01 刻み) で 0.5〜4.0 の外に無し、0.01 刻みでだけ成立する点も無し。**事前登録の規則による判定は A** (分類 中間・n 161 との不一致 2 ≤ 2・反転 1 < 4。合否を pass↔reject だけ・verdict の変化で数えても、不一致を D だけで数えても A。解釈は未確定)。189 s。成果物 `v0c_summary.json`・`v0c_gate.csv`・`v0c_d_design.csv`・`v0c_xm_scan_coverage.csv`・`v0c_vs_v0b.csv` (V0-b で成立した 31 点の対照表)・`v0c_run.log`、図 (git 管理外) `v0c_gate_map.png`。スクリプト `run_v0c.py`。注: 実行時、作業ツリーの `evaluate/runner.py`・`geometry/wall.py`・`opt/driver.py`・`meshing/mesh2d.py` に別作業の未 commit の変更があった (`v0c_summary.json` の `design_tree_dirty`)。V0-c が使うのは `wall_bell_c2.py`・`moc_forward.py`・`moc_kernel.py` と `wall.py` の `NozzleWall` (変更なし) | active |
| `run_0005_tp_cv_continue2` | plan §5.1 #2c (事前登録、diagnostician 2026-10-11): run_0003 を同一メッシュ・**設定変更なし**で +24000 継続 (本段計 60000)。修正前のバイナリ | **結果 A**: $\eta_{exit}$ 0.9791941 (予測漸近値 0.979183 との差 1.1e-5 ≤ 1.5e-4)、継続区間の $|cv|$ 最大 4.25e-5 ≤ 1e-4、推力系 4 列 STEADY (1.5e-4)、ṁ STEADY (1e-3)、残差 NOT CONVERGED (stalled、併記)、NaN/Inf 0 | ref |
| `run_0006_tp_hoopfix_staged` | plan §5.1 #2e: 近軸 hoop 修正 (`axisSegmentRWeight` ON、`/PLANES/rSurfVect`) を取り込んだバイナリ・変換器で基準ケースを作り直し、段階起動 (本段 12000)。作業ツリー `forge-thruster-merge` で実行し、merge 後にここへ移した | NaN/Inf 0、残差 NOT CONVERGED (stalled) | active |
| `run_0007_tp_hoopfix_cont48k` | run_0006 を同一メッシュで +48000 継続 (本段計 60000 = run_0005 と同じ長さ) | run_0005 との差: $\eta_{exit}$ +5.1e-6、$\eta_{wall}$ −3.4e-5、ṁ +5.0e-6、**CV 残差 4.25e-5 → 3.7e-6**。推力系 4 列 STEADY (1.5e-4)、NaN/Inf 0 | active |
