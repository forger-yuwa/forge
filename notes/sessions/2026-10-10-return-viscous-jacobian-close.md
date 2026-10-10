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

## ③ (面エンタルピーの A/B) の状況

[`plans/active/time_integration-implicit-thermal-jacobian.md`](../../plans/active/time_integration-implicit-thermal-jacobian.md) §6.2 に事前登録した。run は case/45 の 05xx (`run_0500`〜`run_0537`)、AWS の `~/forge-wallfit/case/45.isobutane_m6_d155/` で回す。結果は同 plan と case/45 の README に書く。
