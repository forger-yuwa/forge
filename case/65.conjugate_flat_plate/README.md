# 65. 共役平板 (C1/C2) — 厚さ方向の抵抗 (C1) と板の軸方向伝導 (C2) が効く共役熱伝達

計画: [`plans/accepted/boundary-cht-conjugate-benchmarks.md`](../../plans/accepted/boundary-cht-conjugate-benchmarks.md) (§4.3・§4.5・§4.7)。原因切り分けは後継 [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)。判定基準は plan の登録文そのまま。

- 流体: 一様流 M 0.1、Re_L 1e4 (L 10 mm)、定数物性 (Pr 0.72)、平面 2D。前縁の上流 L/2 と後縁の下流 L/2 は slip、上境界 (6δ) は slip。
- 固体: 板 0 ≤ x ≤ L、−b ≤ y ≤ 0 (b/L 0.2)、k_s/k_f 10 (C1) / 100 (C2)、下面 Robin h 1e8・T_h 310 K、両端面は断熱。
- 評価器は case/64 の `eval_conj.py C <run>`、条件選定は `survey_C.py`。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_n16` / `run_0002_dry_n32` / `run_0003_dry_n64` | 乾式 1 step (板 300 K 等温)。壁ダンプ = 板の界面節点の出所 (ローカル float) | 板の壁が y=0、0 ≤ x ≤ L | ref (固体の入力) |
| `run_0004_smoke_c1_n16` | 評価器の通し確認 (ローカル float、C1、6000 step、未収束) | 評価器の C 経路が動くこと | 破棄予定 (登録外) |
| `run_0005_c1_n16` / `run_0006_c1_n32` / `run_0007_c1_n64` / `run_0008_c2_n16` / `run_0009_c2_n32` / `run_0010_c2_n64` | **本番 C1・C2** (AWS FP64、600000 step) | 主判定 6 本 PASS。**前提ゲートの流体収束・G-if は NOT CONVERGED** (slip 境界の揺れ: 後縁下流・前縁上流) | active |
| `run_0013_pre_c1_n64_coupled` / `run_0014_pre_c1_n64_fixed_tw` | B1 予備確認 (後継 plan §4.1.1): `run_0007` の step 600000 から 200 step・毎 step 出力。A は連成継続 (warmup 5000)、B は壁温固定 (warmup 900000) | 受入条件 PASS・採取間隔の検定 FAIL (`AB_PRE_CHECK_run0013_0014.txt`) → 毎 step FP64 抜き出しへ改訂 | ref (予備確認) |
| `run_0011_ab_c1_n64_coupled` / `run_0012_ab_c1_n64_fixed_tw` | **B1 本段** (後継 plan §4.1・§4.1.1b): 同じ起点から各 60000 step。変えるのは `conjugate.warmup` だけ (A 5000 / B 900000)。毎 step 出力を `ab_extract.py` が `ab_series.npz` に抜き出す | `AB_JUDGE_run0011_0012.txt`: 「動的な連成だけが停滞を維持する」を棄却。壁温を固定しても残差・流体変動・前縁の界面熱流束の変動は同じ (B/A ≈ 1)。両 run とも NOT CONVERGED (stalled)。系列は各 run の `ab_series.npz` | active |
| `run_0015_cfl_c1_n64_cfl2` / `run_0016_cfl_c1_n64_cfl05` | **B2′** (後継 plan §4.3): `run_0007` の step 600000・壁温固定から `cfl_pseudo` 2 (12000 step) / 0.5 (48000 step)。CFL × step を揃える。毎 step 抜き出しに `Ux/ro/dt_local/limiter_*` を追加 | 実行中 (2026-10-01 投入)。判定器 `ab_judge_cfl.py` | active |
