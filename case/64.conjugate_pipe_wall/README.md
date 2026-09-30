# 64. 厚肉管の共役 Graetz (A1/A2) — 固体の軸方向伝導が効く共役熱伝達

計画: [`plans/accepted/boundary-cht-conjugate-benchmarks.md`](../../plans/accepted/boundary-cht-conjugate-benchmarks.md) (§4.2・§4.5・§4.7)。判定基準は plan の登録文に §4.6 の事後改訂 (上流延長の感度を U から外して別記) を加えたもの。

- 流体: 管 R 1 mm、M 0.05・Re 1000・Pr 0.72 (case/63 と同じ定数物性・放物入口)、上流 80R・加熱 10R・下流 10R。壁は全長が共役 (physID 3)。
- 固体: R〜2R の殻 (全長)、k_s/k_f 10 (A1) / 100 (A2)、外面 Robin h_o 570 W/m²K・T_c 310 K を加熱区間の辺だけ、他は断熱。
- 参照解 `conjugate_ref.py` (forge と独立、流体 + 固体の節点中心有限体積) と評価器 `eval_conj.py` (forge の流れ場を固定した主参照との比較)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_r16` / `run_0002_dry_r32` / `run_0003_dry_r64` | 乾式 1 step (壁 300 K 等温)。壁ダンプ = 固体界面節点の出所 (ローカル float) | 壁が r=R、全長、`iface_ok` 全 1 | ref (固体の入力) |
| `run_0004_smoke_a1_r16` | 評価器の通し確認 (ローカル float、A1、6000 step、未収束) | 評価器が動くこと | 破棄予定 (登録外) |
| `run_0005_a1_r16` / `run_0006_a1_r32` / `run_0007_a1_r64` | **本番 A1** (AWS FP64、600000 step) | 前提ゲート全 PASS。再評価 (`EVAL_CONJ_L2.txt`、評価版 70b155cd): r16 FAIL (一部判定不能) / r32・r64 全項目 PASS。準定常 `SERIES_CONJ_ALLNODES.txt` PASS | active |
| `run_0008_a2_r16` | 本番 A2 r16 | step 32300 で連成の安全停止 (max\|dTw\| の 10 更新連続増) | ref (失敗の記録) |
| `run_0009_a2_r32` / `run_0010_a2_r64` / `run_0011_a2_r16_df20` | 本番 A2 (r16 は Df_scale 20 で再走) | 前提ゲート全 PASS。再評価: r16 (`run_0011`) FAIL (一部判定不能、参照 4 水準でも FAIL — forge 粗格子誤差の寄与を示唆、参照側は排除していない) / r32・r64 全項目 PASS。準定常 PASS | active |


## 一次データ

- 台帳 `MANIFEST_primary_data_A.sha256` (AWS `~/forge-graetz/case/64.conjugate_pipe_wall` で 2026-09-30 作成)。手元には最終 step (600000) の場・壁/固体ダンプ・`residual_history.csv.gz`・系列 CSV・評価出力のみ。中間スナップショットは AWS 側。
- 延長区間の源項 A/B: `ab_ext/`。参照 4 水準 (`L3`) での r16 再判定: `ab_levels3/`。
