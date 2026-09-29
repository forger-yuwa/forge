# 65. 共役平板 (C1/C2) — 厚さ方向の抵抗 (C1) と板の軸方向伝導 (C2) が効く共役熱伝達

計画: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md) (§4.3・§4.5・§4.7)。判定基準は plan の登録文そのまま。

- 流体: 一様流 M 0.1、Re_L 1e4 (L 10 mm)、定数物性 (Pr 0.72)、平面 2D。前縁の上流 L/2 と後縁の下流 L/2 は slip、上境界 (6δ) は slip。
- 固体: 板 0 ≤ x ≤ L、−b ≤ y ≤ 0 (b/L 0.2)、k_s/k_f 10 (C1) / 100 (C2)、下面 Robin h 1e8・T_h 310 K、両端面は断熱。
- 評価器は case/64 の `eval_conj.py C <run>`、条件選定は `survey_C.py`。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_n16` / `run_0002_dry_n32` / `run_0003_dry_n64` | 乾式 1 step (板 300 K 等温)。壁ダンプ = 板の界面節点の出所 (ローカル float) | 板の壁が y=0、0 ≤ x ≤ L | ref (固体の入力) |
| `run_0004_smoke_c1_n16` | 評価器の通し確認 (ローカル float、C1、6000 step、未収束) | 評価器の C 経路が動くこと | 破棄予定 (登録外) |
