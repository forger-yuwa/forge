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
