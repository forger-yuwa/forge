# 諮問: 共役平板 B3 (limiter 2 vs 0) が「中間」— 解釈と次の一手

AGENTS.md 条件 3・4・7。plan `plans/active/boundary-cht-conjugate-flat-plate.md` §4.4・§5.1 #1e (commit a949af93)。判定出力 `case/65.conjugate_flat_plate/AB_JUDGE_LIM_run0017_0018.txt`、判定器 `ab_judge_lim.py`、抜き出し `ab_extract.py --region-stats`。

## 観測事実 (AWS FP64 `86115cb1`、起点 run_0016 res_48000・cfl_pseudo 0.5・壁温固定)
- A `run_0017_lim_c1_n64_lim2` (12000 step): run_0016 を ±1.6 % で再現、check_convergence NOT CONVERGED (stalled 0.2 桁)。
- B `run_0018_lim_c1_n64_lim0` (12000 step、limiter 0 = 無制限 2 次再構成): check_convergence NOT CONVERGED (**still converging**)、rms_ro 7.36e-8→1.22e-10 (2.7 桁)、rms_roUx 3.25e-6 (peak 5.22e-6)→5.41e-8 (2.0 桁)、rms_roUy 2.5 桁、rms_roe 2.7 桁、なお下降中。最終場 NaN/Inf なし、Pmin 101270.68、Tmin 299.9754。
- 登録判定 B/A (A 最終区間基準、区間 1/2/3): rms_ro 0.491/0.064/0.044、**rms_roUx 13.457/2.214/1.408**、rms_roUy 0.826/0.109/0.073、rms_roe 0.463/0.062/0.043、下流 slip P 0.496/0.054/0.034 → rms_roUx が 1/10 条件を満たさず **中間 (判定不能)**。切り替え直後の過渡で rms_roUx が跳ね、下降途中。
- 残差の局在 (時間積算の残差二乗和、全域 RMS は CSV と同一 step で 8.9e-16 一致):
  - A: 後縁帯 (x/L 0.98〜1.1、y≤0.05L、1240/15665 節点) に res_ro 97.9 %、res_roUx 94.9 %、res_roUy 98.7 %、res_roe 97.9 %。上位節点は x/L 1.004〜1.017 の y=0 と第 1 層 (y 2.5e-6 m)。同じ場所で limiter_Uy の時間標準偏差 0.38〜0.41 (x/L 1.017〜1.042、y=0)、limiter_P 0.23〜0.29 (x/L 1.004〜1.011)、limiter_ro 0.02 (x/L 0.983)。
  - B: 後縁帯 0.4〜0.8 %、前縁帯 (x/L −0.1〜0.02) 91〜93 % (下降中の過渡)。
- 先行 A/B: 壁温の動的更新を止めても停滞は同じ (B1)。cfl 2→0.5 で前縁・上流の振動は消えたが残差は半分で横ばい (B2′)。
- 設定: node・SLAU・Venkatakrishnan (limiter 2、limiterScaled 既定 1・venkatK 既定 0.05)・blockDPLUR・定常。後縁 x/L=1 は no-slip 板の端で下流は slip (y=0)。

## 仮説
H1: 後縁直後 (no-slip→slip の接合と後流の始まり) でリミッタ係数が step ごとに切り替わり、有界振動 (停滞) を維持している。limiter 0 で停滞が解ける。
H2: limiter 0 は後縁で解が振動せずに済むだけで、無制限再構成が本番設定として妥当かは別問題 (M0.1・衝撃なしの流れなので、リミッタの本来の役目は薄い)。

## 質問
1. 登録判定は「中間」のまま記録する。書いてよい解釈の範囲は (H1 の支持の強さ)。
2. 次の一手: (a) B を延長して収束 (check_convergence PASS) まで回す (rms_roUx の過渡が抜けるか)、(b) リミッタの切り替わりを抑える既存キー (venkatK を上げる等) の A/B、(c) 本番 6 本を limiter 0 (と cfl 0.5) で回し直して C の完了条件を判定する。どの順で、何を事前登録すべきか。
3. limiter 0 を C の検証設定として使ってよいか (登録条件の変更になる。発注元 plan §4.5 のレシピとの関係)。
