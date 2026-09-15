# case/48 超音速冷却壁平板 (等温壁 × 低 Re SST)

対象: 等温壁 (`wall_isothermal`, `Ts`) と低 Re SST (`wallTreatmentSST: 0`) の**冷却された超音速乱流境界層**の再現性。
ノズル設計チェーンの等温壁化 ([plan](../../plans/active/tooling-nozzle-isothermal-wall-chain.md) §4.3) の根拠ケース。
壁 BC・SST・粘性流束・node 壁ピンを触ったときに回す。

## 手順

1. メッシュ: `python3 case/48.flat_plate_cooled_m4/gen_mesh.py --y1 3` (平面 2D quad → node 変換 → `check_mesh_quality.py`)。
2. 断熱基準: `python3 gen_runs.py --run run_NNNN_A --mesh fp_y1_3um --wall adiabatic --cfl 2 --main-steps 48000 --out-int 4000`
   (層流暖機 → SST soft/mid → 2 次ランプ 0.5/1/2 → 本段 cfl 2。cfl 4 は前縁で NaN)。
3. 冷却壁: `--wall 300 --stages soft --ramp 0.5,1,2 --ic-from run_NNNN_A` (壁温切替は soft/mid 段を挟む)。
4. 評価: `design/.venv-opt/bin/python tools/cooled_plate_eval.py RUN --ref RUN_A --closure --series --integrals [--su2 DIR] --out RUN/eval`。

## 判定 (2026-09-12 基準値, run_0004 / run_0005)

- 冷却壁 (T_w 300 K, T_w/T_aw 0.26): $C_f$/van Driest II 0.97–0.99、$2St/C_f$ 1.16、δ\*/θ の CONTUR 比 +4 % / 0 %、エネルギー閉合 1.017、series STEADY。
- 断熱: 回復係数 0.88、$C_f$/VD-II 0.87–0.92 (低 $Re_\theta$)。
- SU2 (素 SST 同士): $C_f$・θ・$q_w$ ≤1 %、δ\* ≤4 %。
- これらから 3 % 以上ずれたら回帰とみなす (残差床 rms_ro ≈1e-7 / roe ≈0.15 は前縁特異点由来の既知の床)。
