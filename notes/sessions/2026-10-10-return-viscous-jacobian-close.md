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

## ④ の続き: 粘性ヤコビアン plan の再開 (plan `time_integration-line-viscous-jacobian-faceh`) からの申し送り (2026-10-10)

ユーザの方針「速度でなく筋のいい手法」で、全部入り (値 3・マスク 7) が case/45 で壊れる理由を切り分けている。float 化の段 ③ に関わる観測があるので共有する
(codex 諮問 [`notes/reviews/2026-10-10-lvcaudit-result-diagnose.md`](../reviews/2026-10-10-lvcaudit-result-diagnose.md) の申し送り案を元にした文面)。

> FP64 ビルドでも ISP 0 の block DPLUR は、座標を float にしてから差を取ります (段 ③ の前の経路)。case/45 の監査の記録 (`run_0183` の res_100000、5 本の壁法線のライン、
> `run_0550_lvcaudit`、plan §6.8) では、壁法線のライン面の β・κ に、double の座標から計算した係数との差が最大 10.8 % ありました (壁から 0〜15 番目の節点で 3〜11 %、52 番目まで 1e-3 超、294 面)。
> 例: 壁の隣の dcc が double 2.679e-8 に対し float の座標の差で 3.003e-8 (y の float の ulp 7.5e-9)。ライン外の面は 1e-3 以内。
> 同じ幾何は B0 のスカラーの対角と、キー 5 の熱伝導の対角にも使われます。元の plan §4.2a・段 ③ の対象には既に含まれています。
> 今回の発散への因果と B0 の収束の挙動への影響は未確定です。段 ② と段 ③ の FP64 のバイナリで、LHS の座標の差だけを変えた値 3・マスク 7 の A/B を回します (plan §6.9、run は case/45 の 0560〜0566)。
> 腕の前に、同じ凍結入力で段 ② と段 ③ の残差 (全節点の 6〜8 場と `rhs_s0`) が変わらないことと、係数の変化が double の差 e で説明できることを事前のゲート (`lvcgeom_pregate.py`) で確かめます。
> 訂正: そちらの `run_0419` のビット一致は段 ③ の中の旧式・新式の比較で、段 ② と段 ③ のバイナリの間の比較ではないと理解しました (codex plan-5 レビュー M1)。
> 結果は共有します。残差の固定状態での不変性と、LHS・更新の履歴の変化を分けて検証してください。

- 段 ② と段 ③ のバイナリ (`~/forge-fgeom2-fp64`・`~/forge-fgeom3-fp64`) は読むだけで、再ビルドしない。run は AWS の自分の作業ツリー `~/forge-faceh-audit/case/45.isobutane_m6_d155/` に置く (`~/forge-wallfit` は触らない、`run_0183` はリンク)。
- 監査用のビルド (`-DFORGE_LINE_AUDIT`、`~/forge-faceh-audit`) と対照の通常のビルド (`~/forge-faceh-ctrl`) は、このブランチの 2d58b457 + typedef double。
- メッシュの座標の精度: `run_0183` の `nozzle.h5` の `MESH/COORD` は float64 で、書き出した 605 節点の座標は float32 に丸めると変わる (x 0/605、y 5/605 が一致) = double の精度を持っている。
  今回のずれは LHS で float に落としてから引くことから来ていて、メッシュの側ではない。
