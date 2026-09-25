# G0: node 並進周期の継ぎ目勾配 (plan boundary-node-periodic-gradient-fix §6 G0)

32³ の三重周期 TGV メッシュ (`case/09.Taylor-Green/mesh`、node 変換) に線形場を焼き込み、1 step で勾配を出す。

| ファイル | 内容 |
| --- | --- |
| `bake_linear.py` | 変換済み `Taylor-Green.h5` の VALUE に Ux=10+y、k=1+0.1y、ω=100+2y を焼き込む |
| `aggregate.py` | x 継ぎ目と内部の勾配を集計し plan §6 G0 の閾値で判定 (k/ω の GG は G1 の参考値) |
| `solverConfig.yaml` / `bcondConfig.yaml` | 修正前の測定に使った設定 |
| `solverConfig_sst.yaml` | k/ω 勾配も出す SST 版 (`turbulence.model: sst`、`output.level: 2`) |
| `G0_translational.txt` | 修正前 (`3de26cba`): 継ぎ目 dUx/dy = 2.000000 |
| `G0_translational_m1.txt` | 修正後 (M1 実変位版、`res_1.h5`): LSQ PASS、GG k/ω は継ぎ目法線成分に ~20 ulp·φ/h の丸め |

手順: `python3 bake_linear.py <run>/Taylor-Green.h5` → `FORGE_CUDA_BLOCKSIZE=128 FORGE_CUDA_BLOCKSIZE_SMALL=128 solver_density_cuda/tools/run_case.sh <run>`
(SLAU カーネルは既定ブロック 512 で起動できない) → `python3 aggregate.py <run>/res_1.h5`。
