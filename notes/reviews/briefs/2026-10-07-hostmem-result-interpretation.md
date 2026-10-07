# 諮問: ホストメモリ削減 R1–R3 — 追加診断 B が全量 PASS。受け入れてよいか、残る弱点をどう扱うか (2026-10-07)

関連 plan: `plans/active/architecture-solver-host-memory.md` (§6 登録条件、§6.2 登録結果・追加診断 B の事前登録と判断)。前回諮問: `notes/reviews/2026-10-07-hostmem-regression-fail-diagnose.md`。
エスカレーション条件 7 (result 段の解釈を確定する前、codex `--stage result` の前)。

## 観測事実
- 登録判定 A (§6、相対尺度 m): (c) NaN・初期出力・出力互換・ログ・変換器は全 PASS。(a) step 0 は 1 構成を除き FAIL (前提の 1 ulp 2 値が base 自身で崩れていた)。(b) は 4 構成 7 量で FAIL。**A は FAIL のまま記録**。
- 比較器の欠陥 (非対称な分母・inf を PASS) は呼び出し側でも確認。
- **追加診断 B** (§6.2 で事前登録、絶対 L∞、追加の forge 実行なし、既存 base 3 本・new 3 本): 30 構成の 3147 量 + SERN g3 の 255 量で **FAIL 0・比較不能 0**。整数・文字列 185 量は全 run で厳密一致。
  自己検査: 全 3402 量で base・new の全順列 (36 通り)、計 122472 回評価して S・D・判定の変化 0、base/new 交換でも変化 0。単体試験 8 件 (codex の最小再現 2 つを含む) PASS。
  原本: `case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957_abs/` (`summary.txt`・`all_quantities.tsv`・`selfcheck.txt`)、README の「追加診断 B」。
- A で FAIL の 7 量は B で全部 D/2S = 0.50 (S と D を同じ 1 本の外れ run が決めている): c44dual_ckpt100 `condClampCorrQ_0` (S=D=1.462e23、new r1) と `condR30_0` (4.203e-5、new r1)、c44dual_restart100 `condClampCorrQ_0` (1.433e10、base r1)、c20cell_rk3 `rms_roUz` (3.277e-7)、c20cell_dual `roUz`・`Uz`・`CHECKPOINT/roUzN` (≈4.3e-7、base r2 と new r1)。→ **この 7 量では B は base と new を区別できていない**。
- D/2S の上位: SERN g3 `res_vehicle_base_18_100:twall_z` **1.000** (S_abs 2⁻⁷、D_abs 2⁻⁶ = float32 の刻みで 2 倍、境界ちょうど)、SERN `utau` 0.957・`ypls` 0.950・`omegab` 0.914、c44dual_pindiag の FCT 履歴 0.866/0.856、c56lineimp `roOmega` 0.850。
- メモリ: SERN g3 ホスト VmHWM 5157 → 2630 MiB (2816 → 1436 B/節点)、GPU 不変、ローカル傾き 2737 → 1371 B/節点。§6 の ≤ 1500 は満たす (g4 は入力消失で g3 単点 + ローカル傾き)。
- 実装の他の確認 (ローカル): 初期出力・dual-time checkpoint のビット一致、負例 2 つ (H から `P`・`roN` を外す) の名前つき停止、変換器 398 データセット一致。
- 設計: デバイス側の確保・転送・カーネルは変えていない (diff は `main.cpp` の初期化順序・`hostCell` 経由のホスト参照、`variables.cpp` の条件付き確保、`output.cpp` の名前を共通関数から取る形)。

## 問い
1. R1–R3 を受け入れて (`done` へ) よいか。受け入れるなら、A の FAIL・B の区別できない 7 量・twall_z の境界を、どう記録すべきか (限定事項として)。
2. 受け入れの前に追加で要る確認はあるか (例: 区別できない 7 量 [凝縮クランプ補正量・名目ゼロの roUz] について、ゼロ近傍の扱いや、反復数を事前に決めて増やす試験。twall_z は float32 の刻みの問題として扱ってよいか)。追加するなら合格条件を事前に。
3. 段階 2 (R4+R5 `readMesh` の平坦化、R6 GPU の未使用変数) に進む前に、この結果から決めておくべきこと。
