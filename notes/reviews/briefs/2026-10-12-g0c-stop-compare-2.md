# 諮問ブリーフ (続き): 試験 G の前提「両格子の check_convergence 合格」が G0 で満たせない件 (2026-10-12)

- 前回の諮問: `notes/reviews/2026-10-12-g0c-stop-compare-diagnose.md` (結論: G0 と G0c を同じ FP64・上限 1 で `run_0754_t2_d2` res_142500 から各 60,000 step、両腕の `check_convergence` 合格と各量の準定常を比較の前提に)
- plan: `plans/active/architecture-float-state-double-geometry.md` §5.1 #25

## 観測事実 (新しい)
- `check_convergence.py run_0754_t2_d2`: 再開 (rms_ro 2.2e-7) から 14.45 万 step、全列が「STALLED (plateau — needs scheme change, not more steps)」(rms_ro 2.20e-7 → 3.02e-7、ρu 4.7e-5 → 7.7e-5、ρE 0.21 → 0.34、ρk 2.0e-4 → 1.9e-4、ρω 6.6e2 → 1.4e3、低下 0.0〜0.2 桁)。OVERALL: CHECK FAILURES。
- `--from-floor run_0722_c1_d2` は REFUSED (参照の `run_0722` は Q32 の起点から 30,000 step で 2.1 桁 (ρ) ・0.4〜0.6 桁 (ρk・ρω) しか下がらず、通常の基準 (--drop 3.0) を満たさない)。
- G0 のこれまでのどの計算 (B0・上限 1・point) も、残差は約 2 桁下がったところで横ばいになり、`check_convergence` は合格していない。
- G0c の B0 の FP64 は、Q32 の起点から 30,000 step で rms_ro が 3.6〜3.8 桁下がった (§6.52、`run_0725`・`run_0726`、check_convergence は CHECK FAILURES、どの列が落ちたかは未確認)。

## 問い
1. 前回の前提 (両腕の check_convergence 合格) は G0 では構成上満たせない。試験 G を「停止規則 (§6.11) での比較 + 限界の明記」に組み直してよいか。それとも別の前提 (例: 収束の判定を残差でなく量の窓の変化で置き換える) がよいか。
2. 組み直しの案 (下) に穴はないか。

## 組み直しの案
- 3 本を B で同時に (数値の応答の収集、速度は別に占有で測る): (G1) G0c・上限 1、(G2a・G2b) G0c・上限 50 (B0) の 2 本 (再実行の幅)。どれも FP64 (`~/forge-fgeom7-fp64`)、起点は G0c の Q32 の起点 (`_g0c/q32_init.h5`、= G0 の B0 の停止場を移したもの)。停止規則 (§6.11、上限 20 万 step) まで。
- 比べ方: (i) G1 対 G0・上限 1 (T2 の FP64 の停止値) = 格子の効果。(ii) G2 対 G1 = G0c での解き方 (上限) の効果と、G0c の B0 が停滞しないか。基準は前回の暫定の基準 (θ_r・δ_loc・Q_w ≤ 0.1 %、q_w の局所 ≤ 1 % 全区間、出口のマッハ数の分布、収支) を停止時の値に当てる。各量の時系列と `check_quasisteady`、`check_convergence` の VERDICT と落ちた列を記録し、「収束した解の比較ではない」と書く。
- 物性: `reduce_fields` の CP 1360・Pr 0.72・Sutherland は case/45 の実効の物性 (TP、`thermalMethod 2`、輸送は CEA) と違う。比較は両格子に同じ式を当てた相対の差として使い、絶対の熱負荷とは呼ばない。符号つきのせん断は足す。
