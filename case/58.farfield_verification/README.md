# case/58 遠方境界 `farfield` の検証

node 用の特性型遠方境界 `farfield` (外側状態 = TRRS + 滑らかな超音速の重み、境界面 HLLC) の受入れ試験。
計画・合格条件は [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md) §6、
結果の記録は同 §5.1 #2・#3。run は AWS (`~/forge-pgrad-new/case/58.farfield_verification/`) にあり、手元には入力スクリプトだけを置く。

- `make_box_msh.py`: 構造 hex の直方体 (gmsh 4.1、physID 1–6 = xmin, xmax, ymin, ymax, zmin, zmax、7 = 体積)
- `setup_v1.py RUN VARIANT [--steps N]`: V1 の run を作る (1 × 0.5 × 0.5 m、24 × 12 × 12、全 6 面 farfield、一様 IC を h5 に直接書く)。
  VARIANT a = CPG M0.5 / b = TP 外気 M6 (SERN の EXH/AIR、Y_EXH 0) / c = b を xy 面内で 30° 傾ける / d = b + SST
- `eval_v1.py RUN...`: V1 (a)–(c) の判定 (最終 `res_*.h5` の自由流からの最大相対ずれ ≤ 1e-5)
- `eval_v1d.py RUN`: V1 (d) の判定 (帳簿ダンプ全節点で、対流と k・ω 輸送の残差 / 接する面流束の絶対和 ≤ 1e-5)。
  run は `FORGE_DUMP_LEDGER=ledger.csv FORGE_DUMP_LEDGER_CALLS=1 FORGE_DUMP_LEDGER_NODES=<全節点> FORGE_DUMP_FARFIELD=ffdump` で回す

V1 の `check_convergence.py` は NOT CONVERGED (全列が丸め床で横ばい) になるが、一様 IC から始める自由流保持の試験で
残差は step 0 から float32 の丸め床にあり、低下しようがない。判定は plan §6 V1 の場のずれで行う。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_v1a` | V1(a) CPG M0.5、2000 step、block-DPLUR | 最大ずれ dρ 7.5e-7・dP 3.1e-7・du 6.7e-7 → **PASS**。置換・退避 0 | active |
| `run_0002_v1b` | V1(b) TP M6 (不具合版バイナリ) | dP 4.2e-5 FAIL。側面の全面が置換 (GPU 上で外側状態の組成が 0 になる不具合) | 破棄予定 |
| `run_0003_v1c` | V1(c) 30° 傾き (不具合版バイナリ) | 側面が置換。不具合版 | 破棄予定 |
| `run_0004_v1d` | V1(d) TP M6 + SST (k 479.653、ω 119844)、2 step、帳簿ダンプ全 4225 節点 | 輸送残差 / 規模 最大 5.8e-7 → **PASS** (`V1D_VERDICT.txt`)。置換・退避 0 | active |
| `run_0005_v1b_fix` | V1(b) 修正版バイナリ (`ff_copy`) | dρ 4.9e-6・dP 6.7e-6・du 6.1e-7・dY 0 → **PASS**。step 0 で 1e-6 → 500 step で 6e-6、以後横ばい | active |
| `run_0006_v1c_fix` | V1(c) 修正版バイナリ | dρ 3.5e-6・dP 6.3e-6・du 4.7e-7・dY 0 → **PASS**。置換・退避 0 | active |
| `dbg_b` | run_0002 の 1 step 複製 (不具合の切り分け用) | 修正版で置換 0 を確認 | 破棄予定 |
