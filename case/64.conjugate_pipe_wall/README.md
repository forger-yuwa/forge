# 64. 厚肉管の共役 Graetz (A1/A2) — 固体の軸方向伝導が効く共役熱伝達

計画: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md) (§4.2・§4.5・§4.7)。判定基準は plan の登録文そのまま。

- 流体: 管 R 1 mm、M 0.05・Re 1000・Pr 0.72 (case/63 と同じ定数物性・放物入口)、上流 80R・加熱 10R・下流 10R。壁は全長が共役 (physID 3)。
- 固体: R〜2R の殻 (全長)、k_s/k_f 10 (A1) / 100 (A2)、外面 Robin h_o 570 W/m²K・T_c 310 K を加熱区間の辺だけ、他は断熱。
- 参照解 `conjugate_ref.py` (forge と独立、流体 + 固体の節点中心有限体積) と評価器 `eval_conj.py` (forge の流れ場を固定した主参照との比較)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_r16` / `run_0002_dry_r32` / `run_0003_dry_r64` | 乾式 1 step (壁 300 K 等温)。壁ダンプ = 固体界面節点の出所 (ローカル float) | 壁が r=R、全長、`iface_ok` 全 1 | ref (固体の入力) |
| `run_0004_smoke_a1_r16` | 評価器の通し確認 (ローカル float、A1、6000 step、未収束) | 評価器が動くこと | 破棄予定 (登録外) |
