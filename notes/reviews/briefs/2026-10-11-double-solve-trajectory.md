# 諮問ブリーフ: §6.30 の 1 step の診断が判定不能 — 次に「float のビルドで線形の解きだけ double」の軌跡の腕を足してよいか (2026-10-11)

- 依頼者: 主セッション (AGENTS.md のエスカレーション条件 3: §6.30 が判定不能、4: 承認済みの手順に無い run の前、1: §6 を足す)
- plan: `plans/active/architecture-float-state-double-geometry.md` (§6.26〜§6.31)

## 0. 問い
1. §6.31 の結果 (1 step の精度の差は再実行の揺れと同じ桁で段を分けられない) を受けて、次に下の A/B を回してよいか。組み方・判定に穴はないか。
2. §6.27 の FP64 の腕 (`run_0487`) と float の腕 (`run_0488`・`0489`) を、新しい腕の比べる相手としてそのまま使ってよいか (同じバイナリ・同じ入力・同じ環境変数で、時刻だけ違う)。

## 1. 観測事実
- §6.27 (確定): 共通の Q32 の起点から、FP64 は θ_r が ≤ 0.12 % しか動かず、float の 2 本は +1.0〜1.9 % 動く (精度依存を支持)。
- §6.29: SST の輸送の更新を止めても float は凍結した FP64 から δ +0.62〜+1.17 ポイント離れる (判別不能)。
- §6.31: 共通の入力からの最初の 1 回の更新の精度の差は、再実行の揺れの 1.5〜3 倍 (10 倍に届かない)。差の内訳はほぼ要求の更新 d の側 (‖Δbb + Δd‖/‖Δa‖ 0.83〜1.00) だが揺れの範囲内。軸の近くでは残差の精度の間の差が R 自身の大きさと同じ桁 (‖ΔR‖/‖R‖ ≈ 1.1〜1.2、R の再実行の差は未測定)。初期化で ρE が精度ごとに組み直され、Q0 は ρE で丸めの水準でずれる。
- 既存の切り替え `time.deltaT.implicitSolvePrecision: 1`: float のビルドのまま、block DPLUR の線形の解き (`implicit_defect_correction_block_d<double>`、ライン陰解法の line_prev・line_next もこのカーネル) を double で行う (`timeIntegration_d.cu:1710-1770`、`solverConfig.cpp:1390-1414`)。blockDPLUR 1・lowMachPrecond 0/1・timeIntegration 11 が条件で、B0 は満たす。過去の経緯 (メモリ `mixed-precision-axisym-refuted`): 軸対称の近軸の固着で double の解きが効いたが、真因は粘性の対角の幾何で、double の解きは悪条件の左辺を押し切る対症療法だった。

## 2. 呼び出し側の案 (検証していない)
- 腕: C = float + `implicitSolvePrecision: 1`、C′ = C の再実行。§6.26 と同じ共通の入力 `_tr/q32_init.h5`・格子・キー 1・環境変数 (`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500`)・バイナリ `~/forge-fgeom7-f32`、各 30,000 step・500 ごと、インスタンス B。違いは設定の 1 行だけ。比べる相手は §6.27 の A (`run_0487`、FP64) と B・B′ (`run_0488`・`0489`、float)。
- 判定 (θ_r の 3 断面、W1 20〜25k・W2 25〜30k、Δ は §6.26 と同じ): 
  - 「解きの精度でずれが消える」: 両窓・3 断面で |mean(ΔC, ΔC′) − ΔA| ≤ 0.1|D|、かつ |mean(ΔB, ΔB′) − mean(ΔC, ΔC′)| ≥ 0.5|D|。
  - 「解きの精度ではずれは消えない」: 両窓・3 断面で mean(ΔC, ΔC′) − ΔA が D の向きに 0.5|D| 以上、かつ |ΔC − ΔC′| の 10 倍を超える。
  - それ以外は判別不能。
- 解釈の範囲: 消えれば「線形の解き (左辺・ライン Thomas) の float の丸めが、ずれの通り道に要る」。消えなければ「残差の評価・状態の更新の側だけでずれが育つ」(定常の固定点の偏り)。どちらでも演算の特定ではない。

## 3. 読んでよいファイル
- float plan §6.26〜§6.31
- `solver_density_cuda/cuda_forge/timeIntegration_d.cu` (1700〜1780)、`solver_density_cuda/input/solverConfig.cpp` (1385〜1415)
- `case/45.isobutane_m6_d155/tr.sh`・`tr_an.py`
