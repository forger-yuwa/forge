# 諮問ブリーフ: 軸対称 node の静止ガスが非一様格子でだけ動き、圧力が下がる (2026-09-27)

AGENTS.md 条件 **2** (発散対処 2 回で解けない — 同 run は warmup 200 で固体 454.6 K、warmup 20000 で固体 276 K) と **4** (原因を書く前)。

## 再現条件
- `case/62.conjugate_disk/run_0014_disk_r32g_w20k` (AWS `~/forge-cht-axi/case/62.conjugate_disk/`、FP64 ビルド = commit `6f64b873` の flowFormat を double にしたもの、`FORGE_CUDA_BLOCKSIZE=128`)。
- 静止ガス層 x∈[0, H=5 mm] × r∈[5, 20] mm (軸を含まない)、x=0 等温 350 K、x=H 等温 325 K (warmup 中は連成なし)、r 両端 slip。`mesh.isAxisymmetric: 1`・`axisymMethod: 0`・node・`nodeWallDirichlet: 1`、定数物性、p0 1013.25 Pa、一様 IC 325 K・静止。
  格子: x 方向 16 一様 × r 方向 32 **等比 1.1 (r が大きいほど広い)**。
- 対照 `run_0013_disk_r32u_w20k`: 同一で r 方向 32 **一様**。

## 観測 (warmup 終了 step 20000、連成前)
| | max|U| [m/s] | 位置 | P [Pa] | Σρ·V·r (指標、V は mesh.h5 の CELLS/volume) |
| --- | --- | --- | --- | --- |
| 一様 r32 (`run_0013`) | 4.3e-4 | (x 3.44, r 5.00) mm | 1027.282–1027.283 | 1.0184e-8 → 9.947e-9 (step 1e4 以降一定) |
| 非一様 r32 (`run_0014`) | **0.231** | (x 4.69, r 6.19) mm | **728.3–742.4** | 1.0184e-8 → **7.147e-9 → 7.125e-9** |
- 非一様では共役壁 (x=H) の壁熱流束 `iface_Qf_eff`/A が r≈5.9–6.6 mm で**隣接節点ごとに符号反転** (+4087 / −3614 / −3198 / +2962 W/m²、物理的には ≈120 W/m²)。
- 質量指標は V が r 重み済みか未確認なので**絶対値は信用しない**。ただし一様と非一様で同じ指標を使っている。
- 同じ非一様格子系列の他 3 本 (一様 8/16/32) は正常に連成が収束中。円環 (case/61、一様格子) も正常。

## 読んでよいファイル
- `case/62.conjugate_disk/gen_mesh.py`・`template/solverConfig.yaml`・`bcondConfig.yaml`・README
- `solver_density_cuda/variables.cpp`:495–640 (r 重み幾何・hoop 用の A_planar)、軸対称の hoop ソースを入れているカーネル (grep `A_planar` / `uy_over_r` in `solver_density_cuda/cuda_forge/`)
- 関連 memory 由来の既知事項 (私のメモ): 軸対称の pRef ゲージと hoop ソースの不整合を 2026-08 に修正済み、r 重み面積の閉包 Σ r_f S_f = (0, A_planar) が高 AR で float32 桁落ちした件 (FP64 では?)。
- run の生データは AWS にあり、ローカルにはまだ無い (必要なら私が取ってくる)。

## 仮説 (未確認)
1. 非一様格子で node の r 重み幾何の閉包 (Σ_f r_f S_f と hoop 用 A_planar の整合) が崩れ、静止ガスに偽の圧力勾配/ソースが立つ。
2. slip 端 (r=5, 20 mm) の境界処理が非一様格子の端セルで非保存。
3. 格子ではなく低マッハの市松 (odd-even) モード。

## 諮りたいこと
1. 最も確からしい仮説と、それを判別する**最小の A/B を 1 つ** (例: 同じ非一様格子を平面 (`isAxisymmetric: 0`) で回す / 非一様を r 方向の向きだけ反転 / 数百 step の質量の厳密な保存量で見る)。
2. これは本 plan (CHT 軸対称) の範囲か、流体の別 plan に切り出すべきか。V-ax2b (非一様) の判定をどう保留・記録すべきか。
