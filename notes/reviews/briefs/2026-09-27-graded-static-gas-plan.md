# 諮問ブリーフ: 非一様格子の静止ガス偽流れ — 新 plan の §4 / §6 (2026-09-27)

AGENTS.md 条件 **1** (新規 plan の §4 設計方針・§6 検証計画)。

## 読んでよいファイル
- 新 plan: `plans/active/axisymmetric-graded-grid-static-gas.md` (全文、短い。§3 が観測事実)
- 既知欠陥: `notes/investigations/node-slip-tangential-density-spurious-flow.md`
- 前回の諮問: `notes/reviews/2026-09-27-axisym-nonuniform-static-gas-diagnose.md` (閉包 A/B の登録)
- 構成: `case/62.conjugate_disk/gen_mesh.py`・`template/solverConfig.yaml`・`bcondConfig.yaml`、`case/62.conjugate_disk/README.md` の run 一覧
- 出力 (ローカルにある): `case/62.conjugate_disk/run_0015_hoop_A_hoop0/` (500 step、10 step ごとの res_*.h5 と壁ダンプ)、`run_0014_disk_r32g_w20k/res_20000.h5`、`run_0013_disk_r32u_w20k/res_20000.h5`。**h5 は必要な配列だけ numpy で**
- 必要なら `solver_density_cuda/cuda_forge/boundaryCond_d.cu` の slip (node)、`solver_density_cuda/variables.cpp`:495-640

## 諮りたいこと (推奨を 1 つに)
1. §4 の H1/H2/H3 のどれが最も確からしいか、手元の出力 (特に `run_0015_hoop_A_hoop0` の 10 step ごとの場: 偽流れが最初に立つ位置と時刻 — slip 境界上か、壁際か、格子の粗い側か) から言えることがあれば。
2. **最初に回す判別 A/B を 1 つ** (§4 の (i) 加熱なし / (ii) 平面 / (iii) slip→no-slip or periodic、またはその他)。変える 1 変数・step 数・見る量・判定を事前登録できる形で。
3. §6 の静止保持の閾値 (max|U|、P) を何で決めるべきか。
