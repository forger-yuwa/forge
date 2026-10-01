# 65. 共役平板 (C1/C2) — 厚さ方向の抵抗 (C1) と板の軸方向伝導 (C2) が効く共役熱伝達

計画: [`plans/accepted/boundary-cht-conjugate-benchmarks.md`](../../plans/accepted/boundary-cht-conjugate-benchmarks.md) (§4.3・§4.5・§4.7)。原因切り分けは後継 [`plans/accepted/boundary-cht-conjugate-flat-plate.md`](../../plans/accepted/boundary-cht-conjugate-flat-plate.md)。判定基準は発注元 plan の登録文に、後継 plan §4.6・§4.6.1・§4.6.2 の事後改訂 (limiter 0・cfl 0.5・同一格子 restart・④ tol_solid の訂正・C1 の軸方向熱量を参考表示・C2 の固体効果を登録量 Δθ_ax で判定・n16 の参照 4 水準) を加えたもの。**結論は n32・n64 の限定結果** (後継 plan §1)。

- 流体: 一様流 M 0.1、Re_L 1e4 (L 10 mm)、定数物性 (Pr 0.72)、平面 2D。前縁の上流 L/2 と後縁の下流 L/2 は slip、上境界 (6δ) は slip。
- 固体: 板 0 ≤ x ≤ L、−b ≤ y ≤ 0 (b/L 0.2)、k_s/k_f 10 (C1) / 100 (C2)、下面 Robin h 1e8・T_h 310 K、両端面は断熱。
- 評価器は case/64 の `eval_conj.py C <run>`、条件選定は `survey_C.py`。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_n16` / `run_0002_dry_n32` / `run_0003_dry_n64` | 乾式 1 step (板 300 K 等温)。壁ダンプ = 板の界面節点の出所 (ローカル float) | 板の壁が y=0、0 ≤ x ≤ L | ref (固体の入力) |
| `run_0004_smoke_c1_n16` | 評価器の通し確認 (ローカル float、C1、6000 step、未収束) | 評価器の C 経路が動くこと | 破棄予定 (登録外) |
| `run_0005_c1_n16` / `run_0006_c1_n32` / `run_0007_c1_n64` / `run_0008_c2_n16` / `run_0009_c2_n32` / `run_0010_c2_n64` | **本番 C1・C2** (AWS FP64、600000 step、登録設定 limiter 2・cfl 2) | 主判定 6 本 PASS だが、**前提ゲートの流体収束・G-if が NOT CONVERGED のため判定不能** (保持)。再走 `run_0022〜0027` の restart 元 | ref (判定不能) |
| `run_0013_pre_c1_n64_coupled` / `run_0014_pre_c1_n64_fixed_tw` | B1 予備確認 (後継 plan §4.1.1): `run_0007` の step 600000 から 200 step・毎 step 出力。A は連成継続 (warmup 5000)、B は壁温固定 (warmup 900000) | 受入条件 PASS・採取間隔の検定 FAIL (`AB_PRE_CHECK_run0013_0014.txt`) → 毎 step FP64 抜き出しへ改訂。AWS の毎 step 全場は 2026-10-01 に削除 (`PRUNED_2026-10-01.txt`、res_0/res_200 は保持) | ref (予備確認) |
| `run_0011_ab_c1_n64_coupled` / `run_0012_ab_c1_n64_fixed_tw` | **B1 本段** (後継 plan §4.1・§4.1.1b): 同じ起点から各 60000 step。変えるのは `conjugate.warmup` だけ (A 5000 / B 900000)。毎 step 出力を `ab_extract.py` が `ab_series.npz` に抜き出す | `AB_JUDGE_run0011_0012.txt`: 「動的な連成だけが停滞を維持する」を棄却。壁温を固定しても残差・流体変動・前縁の界面熱流束の変動は同じ (B/A ≈ 1)。両 run とも NOT CONVERGED (stalled)。系列は各 run の `ab_series.npz` | active |
| `run_0015_cfl_c1_n64_cfl2` / `run_0016_cfl_c1_n64_cfl05` | **B2′** (後継 plan §4.3): `run_0007` の step 600000・壁温固定から `cfl_pseudo` 2 (12000 step) / 0.5 (48000 step)。CFL × step を揃える。毎 step 抜き出しに `Ux/ro/dt_local/limiter_*` を追加 | `AB_JUDGE_CFL_run0015_0016.txt`: 中間 (判定不能)。cfl 0.5 で前縁の界面 q と上流 slip の振動は消えた。残差は 0.43〜0.51 倍で横ばい、後縁直後の揺れは残った (0.84)。両 run とも NOT CONVERGED | active |
| `run_0017_lim_c1_n64_lim2` / `run_0018_lim_c1_n64_lim0` | **B3** (後継 plan §4.4): `run_0016` の最終場 (cfl 0.5・壁温固定) から `space.limiter` 2 / 0 (0 = 無制限の 2 次再構成)、各 12000 step。残差の領域別の二乗和とリミッタの節点ごとの統計も採る | `AB_JUDGE_LIM_run0017_0018.txt`: 中間 (判定不能)。limiter 0 は 2.0〜2.7 桁下降中 (still converging)、rms_roUx は過渡中。limiter 2 の残差は後縁帯 (x/L 1.004〜1.017) に 95〜99 % 局在し、そこでリミッタ係数が時間変動 | active |
| `run_0019_lim0ext_c1_n64` | **B3 の延長** (後継 plan §4.5): `run_0018` の最終場から同じ設定 (limiter 0・cfl 0.5・壁温固定) で 36000 step (累計 48000) | `AB_JUDGE_LIMEXT_run0017_0019.txt`: 判定不能。rms_roUx だけ A の 0.21 倍でゆっくり下降中、他の列と下流の P は A の 1/200〜1/300。連結 check_convergence は NOT CONVERGED (rms_roUx still converging、他は 3.7 桁) | active |
| `run_0020_lim0_c1_n64_cfl2` | **cfl 判別** (後継 plan §4.6): `run_0019` の最終場から limiter 0・`cfl_pseudo` 2、6000 step。本番の cfl を決める | `AB_JUDGE_CFL2LIM0.txt`: cfl 0.5 採用 (条件 4 のパワー比 0.624 で FAIL。窓内の下降トレンドを拾った誤検出で、判定は登録どおり) | ref |
| `run_0022_c1_n64_lim0` 〜 `run_0027_c2_n32_lim0` | **C の本番再走** (後継 plan §4.6、事後改訂設定 limiter 0): `run_0005〜0010` の step 600000 から同一格子 restart、150000 step、連成あり、`tol_solid` を登録式の値に。順は 0022 C1 n64 → 0023 C2 n64 → 0024 C1 n16 → 0025 C2 n16 → 0026 C1 n32 → 0027 C2 n32。still converging は 1 回だけ `<run>_ext` で延長 | **C の結論の根拠** (後継 plan §1・§6.2)。収束 6 本とも PASS (3.9〜7.8 桁、全列下降中)。準定常・G-if (登録式 tol_solid)・G-cons 全 PASS。主判定は C2 n16 だけ判定不能 (q の U 1.03 % > 1 %) → 参照 4 水準 (`ab_levels3/`) で **FAIL** (2.35 % + 0.98 % > 3 %)、他 5 本 PASS。一次データ台帳 `MANIFEST_primary_data_C.sha256`。各 run の `CONVERGENCE_CHECK.txt`・`CHT_INTERFACE_VERDICT*.txt`・`CHT_BALANCE_VERDICT.txt`・`SERIES_CONJ.txt`・`EVAL_CONJ.txt`・`residual_history.png` | active |
