# case/58 遠方境界 `farfield` の検証

node 用の特性型遠方境界 `farfield` (外側状態 = TRRS + 滑らかな超音速の重み、境界面 HLLC) の受入れ試験。
計画・合格条件は [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md) §6、
結果の記録は同 §5.1 #2・#3。run は AWS (`~/forge-pgrad-new/case/58.farfield_verification/`) にあり、手元には入力スクリプトだけを置く。

- `make_box_msh.py`: 構造 hex の直方体 (gmsh 4.1、physID 1–6 = xmin, xmax, ymin, ymax, zmin, zmax、7 = 体積)
- `setup_v1.py RUN VARIANT [--steps N]`: V1 の run を作る (1 × 0.5 × 0.5 m、24 × 12 × 12、全 6 面 farfield、一様 IC を h5 に直接書く)。
  VARIANT a = CPG M0.5 / b = TP 外気 M6 (SERN の EXH/AIR、Y_EXH 0) / c = b を xy 面内で 30° 傾ける / d = b + SST
- `eval_v1.py RUN...`: V1 (a)–(c) の判定 (最終 `res_*.h5` の自由流からの最大相対ずれ ≤ 1e-5)
- `v0u_reject.py`: V0u (v) 非対応構成の起動拒否 (7 構成 + 対照)
- `setup_v2a.py` / `eval_v2a.py`: V2a 音響反射 (薄板チャネル、dual-time、プローブ時系列の短−長差)。リミッタ基準値は自由流で固定
- `setup_v2b.py` / `eval_v2b.py`: V2b 保存収支 (帳簿全節点 + 面ダンプ、最終場 restart の 1 評価 `_bal`)・V2f 局所逆流
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
| `v0u_reject_*` (一時) | V0u (v) 起動拒否。実行ごとに作って消す | 7 構成拒否・対照完走 → **PASS** (`V0U_REJECT_VERDICT.txt`) | 破棄済み |
| `run_0010`–`run_0018_v2a_*` | V2a 初回 (dt 5e-6、プローブ 6 桁) | dt/2 との差 7.6 % で時間精度不足。プローブが 1 Pa 刻み | 破棄予定 |
| `run_0020`–`run_0028_v2a_*` | V2a dt 1.25e-6 (リミッタ基準値は自動) | 差の最大 1.1 % が入射通過時刻に出た = 短・長で L_ref が違い離散化が別物 | 破棄予定 |
| `run_0040_v2a_m03_ff_short` / `run_0041_v2a_m03_long` | V2a 本試験 M0.3 SLAU、dt 1.25e-6、nSub 20、リミッタ基準値固定 | 反射率 **0.099 %** → PASS (`V2A_VERDICT.txt`) | active |
| `run_0042_v2a_m03_ff_short_slau2` / `run_0043_v2a_m03_long_slau2` | 同 SLAU2 | 反射率 0.118 % → PASS | active |
| `run_0044_v2a_m03_ff_short_nsub40` / `run_0045_v2a_m03_ff_short_dthalf` | 時間精度 (nSub 40 / dt 6.25e-7) | run_0040 との差 0.15 % / 0.58 % (≤ 1 %) → PASS | active |
| `run_0046_v2a_m0_slip_short` / `run_0047_v2a_m0_long` / `run_0048_v2a_m0_ff_short` | M0: slip 対照 / 長領域 / farfield | slip 91.3 % (≥ 90 %、試験が反射を検出できる)・farfield 0.203 % → PASS | active |
| `run_0030_v2b_x_ek0` / `run_0031_v2b_x_ek1` (+ `_bal`) | V2b 流出配置 (+x M0.5)、sstEnergyIncludesK 0/1、3000 step | 恒等式 ≤ 1e-9・全体収支 ≤ 6e-9・Y 置換 → PASS (`V2B_VERDICT.txt`) | active |
| `run_0032_v2b_obl_ek0` / `run_0033_v2b_obl_ek1` (+ `_bal`) | V2b 流入配置 (3 面から斜め流入) | 同上 PASS。残差 4.4 桁低下で float 床 | active |
| `run_0034_v2b_x_sfr2` / `run_0035_v2b_obl_sfr2` (+ `_bal`) | V2b 化学種 S3 経路 (speciesFaceReconstruction 2) | 同上 PASS | active |
| `run_0036_v2f` (+ `_bal`) | V2f: 自由流 +x M2、x ≥ 0.8 m 帯を U_x −0.95a・Y 0.13 | 評価 1 で xmax 169 面が ṁ<0・外側組成 = 外気、恒等式 3 評価とも一致 → PASS (`V2F_VERDICT.txt`) | active |
