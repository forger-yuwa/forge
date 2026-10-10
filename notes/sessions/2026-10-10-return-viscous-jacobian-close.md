# 引き継ぎ ④ の戻り: 粘性ヤコビアン plan を閉じた (2026-10-10)

元のセッション (速度 plan・float 化 plan の担当) への申し送り。正本は [`plans/accepted/time_integration-line-viscous-jacobian.md`](../../plans/accepted/time_integration-line-viscous-jacobian.md) §6.17 と §6.1。このメモは写しとポインタ。

- ブランチ: `feature/faceh-audit-viscjac-close` (作業ツリー `/home/sano/work/forge-faceh`)。元のブランチへのマージは区切りで相談する。
- 結果: status `done`、結論は「本線不採用・診断として残置」で `plans/accepted/` へ移した。コードは変えていない。
  - codex result 段: [`notes/reviews/2026-10-10-time_integration-line-viscous-jacobian-result.md`](../reviews/2026-10-10-time_integration-line-viscous-jacobian-result.md) (GO-with-changes、C0/M3/m1、全件採用)。
  - 訂正したこと: U-J の PASS は行列全体の最大値で正規化した判定に限る (登録した列ごとの判定は未確認)、U2 の「収束」は温度の解析解誤差だけが根拠、E1 は比較無効、密度の列は「除いても壊れる・同じ機構に必要かは未確定」。case/45・case/52 の README の該当行も直した。

## 元のセッションに反映してほしいこと (このブランチでは触っていない)

1. **速度 plan のリンク**: `plans/active/time_integration-line-implicit-speed.md` §5.1 #6 の `[time_integration-line-viscous-jacobian](time_integration-line-viscous-jacobian.md)` は、移動でリンクが切れる。`../accepted/time_integration-line-viscous-jacobian.md` に直す。
2. **U0 の扱い**: U0 は未実施のまま閉じた。値 0 については速度 plan §6.0 の B1 の短期の一致 (`run_0276`〜`run_0281_abB_*`、`_band_ab/cold_pair/AB_abB.json`、4d394a71 と 64ed2cd6 の比較、値 0・キー 5・方向別・上限 50) を補助証拠として引用した。
   値 1・point・ISP 1 は確かめていない。`procedures/solver-settings.md` に「値 1 は 5ab83056 以降の不変を確かめていない」と書いた。
3. **float 化 plan への適用範囲**: 値 2・3 の動作は FP64 ビルドでだけ確かめた (FP32 は起動も未確認)。float の状態のビルドで値 2・3 を使うなら、再開の条件 (§6.17) の最初の試験 (FP32 での起動と U2・U3) を先に行う。
4. 総時間の引用 (§6.17): L5 1.18 h vs L0 1.16 h (速度 plan §6.15) は判別不能として本 plan の結論に使った。留保として「専有単価の推定・NOT CONVERGED・TRANSIENT-UNSETTLED・B0 は到達条件が違う」を添えた。

## ③ (面エンタルピーの A/B) の結果

正本は [`plans/active/time_integration-implicit-thermal-jacobian.md`](../../plans/active/time_integration-implicit-thermal-jacobian.md) §6.2 (事前登録)・§6.3 (結果)。§5.1 #5 は「探索試験の判定済み」として閉じた。
codex 諮問: [`notes/reviews/2026-10-10-faceh-floor-result-diagnose.md`](../reviews/2026-10-10-faceh-floor-result-diagnose.md)。

**float 化の担当への申し送り** (諮問の文面のとおり):

> `case/45.isobutane_m6_d155/run_0354_m9_L5cut/res_40000.h5` を出発点に、同一の `lineM_fp64` (sha256 `05ad8bdf…`)・point 設定で、面エンタルピーの評価の精度だけを変えた各 2000 step・3 本の探索試験を行った。
> 共通の評価器 d のエネルギーの残差は B/A = 0.990128、登録の規則では `NOT_SUPPORT`。この条件・期間では、面エンタルピーの double の評価を必要とする 10 % 以上の改善を支持しなかった。
> 面エンタルピーの精度の影響がゼロであること、長期の残差の床、FP32 の状態での妥当性、物理の精度・速度については未判定である。他の混合精度の経路も除外していない。

- 補足: B の水準は 3 本とも約 1 % 低かった (共通の評価でも範囲が分かれた) が、登録外の観測で、切替による再現可能な改善かは未確定。k・ω の残差は B のほうが約 1〜1.7 % 高い。
- FP32 のビルドでは、診断の切替の戻り値が `flow_float` (float) に丸められるので、今回の介入と同じにはならない。今回の結果を FP32 の精度の配分にそのまま使わない。
- run: AWS `~/forge-wallfit/case/45.isobutane_m6_d155/` の `run_0500`〜`run_0505` (軌道)・`run_0510`〜`run_0537` (評価)、判定 `_band_ab/cold_pair/fh_floor_judge.json`。
