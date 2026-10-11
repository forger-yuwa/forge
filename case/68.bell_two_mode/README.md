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
| `v0_moc_gate/` (CFD なし) | plan §6 V0: 前進 MOC ゲート (`geometry/moc_forward.py`) の A/B。(D) の格子 135 点 (θ_max {25,32.5,40}°×θ_e {5,10,15}°×x_m/L {0.1,0.2,0.3}×L_k 60〜100 %) の成立判定と、成立 8 点・TOP (θn 25.8°・θe 10°・R_u=R_d=1) 5 点のゲート。ガス EXH 凍結 TP、初期線 Hall (R_t 1)、主設定 n 321 (確認 n 161・641、参考の整合初期線 n 321/641/1281) | 主: (D) 成立 8/135・通過 3/8 (37.5 %) → 中間 (**暫定**: (D) の族の不成立率と Hall 初期線の不整合を diagnostician に諮問中)。TOP 0/5 (C⁻ 族の交差、始まり x 3.9〜5.0 r_t)。Hall n 641 はスロート直後の初期線の不整合で全点不受理 (アーチファクト)。`v0_summary.json`・`v0_gate.csv`・`v0_d_feasibility.csv`、図 (git 管理外) `v0_gate_map.png`・`v0_nets.png`、再生成 `run_v0.py` | active |
