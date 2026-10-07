# 諮問: ホストメモリ削減 (R1–R3) の回帰で、事前登録した (a) step 0 一致と (b) D ≤ 2S が FAIL — 解釈と、判別のための追加試験 (2026-10-07)

関連 plan: `plans/active/architecture-solver-host-memory.md` (§4 設計、§6 合格条件 [事前登録]、§5.1 #3–#5)。エスカレーション条件 3 (事前登録の比較が FAIL)。
結果の原本: `case/66.hostmem_regression/README.md` (判定・run 一覧)、`case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957/*.txt` (構成ごとの比較)、比較スクリプト `case/66.hostmem_regression/compare_runs.py`。
コード: base = 9c9f623c (計測ログのみ追加、挙動は変更前と同一)、new = 93e55957 (R1: `mat_ns` を作らない、R2: ホスト面変数 `p` を確保しない、R3: ホストのセル変数を H だけ確保 + 共通関数 + `hostCell` アクセサ + 初期化の順序変更)。設計上、**デバイス側の確保・転送・カーネルは変えていない**。

## 観測事実
- 30 構成 (forge 25・変換器 5) + SERN g3、base 3 回・new 3 回 (AWS g5、`FORGE_CUDA_BLOCKSIZE=128`、同一入力)。全 run rc 0・NaN/Inf なし。
- **PASS**: (c) NaN、初期出力 (保存量・原始量・幾何量のビット一致、level 2 の 100 データセットを含む)、出力互換 (ファイル集合・データセット集合・shape・dtype・属性)、ログ行 (行の集合で比較)、変換器 5 構成 (base 内・new 内・base 対 new の全データセットで差 0)、dual-time の分割 (100 + checkpoint 再開 100) 対連続 200 (`split_vs_cont_c09.txt`: 保存量はビット一致、トレーサ roXi の m は base 3.31e-6/2.81e-6/6.0e-7、new 2.71e-6/2.76e-6/6.0e-7 = 反復ばらつきと同程度)。再開直後の `/CHECKPOINT` は入力とビット一致 (ローカル確認)。ψ 退避 `snapshot: 269 of 269`。
- **(a) step 0 FAIL** (c20cell_impdiag 以外の全構成): §6 は「既知の 1 ulp の 2 値のどちらか」を前提にしたが、**base 自身が 3 回で 3 値に割れる列がある** (幅: 2D 定常 19–75 ulp、SERN 52 ulp、dual-time は step 0 に内反復 22 行が入り 1e6 ulp 級)。new と最も近い base 値との距離は大半の構成で base 内の幅以下。上回ったのは c44steady (21961 / 13561 ulp)、c52cht (29 / 19)、c26optin (239 / 199)、c26optin_env (146 / 94)、c44dual_restart100 (1.45M / 1.30M)。c52cht は追加で各 4 回回すと base 7 回中 4 回が new と同じ値を出した (`rms_roUx`)。
- **(b) D ≤ 2S FAIL (4 構成・7 量)** — 保存量・原始量 (ro・roUx・roe・P・T・ρY・凝縮モーメント・k/ω・γ) は全構成で D ≤ 2S (D/2S 最大 0.55–0.95、SERN g3 の最悪は `res_vehicle_base_18_100:twall_z` で 1.00 = 境界ちょうど)。FAIL はいずれも名目ゼロか疎な診断量で、base 自身の run 間で桁が動く:
  - c44dual_ckpt100 `condClampCorrQ_0` (S 44.5 / 1.0、D 1.79e5): run ごとの最大 base 8.2e17・8.1e17・3.6e19、new 1.5e23・8.1e17・1.7e17 (c44dual_pindiag の base r3 も 1.45e23)。`condR30_0` D 3.14・S 1.04。
  - c44dual_restart100 `condClampCorrQ_0`: base 1.4e10・8.0e6・1.8e3、new 7.9e6・7.9e6・3.1e2 (D 4.39e3、S 1.0)。
  - c20cell_rk3 `rms_roUz` (擬似 2D の spanwise 残差 1e-9〜3e-7、rms_roUx は 6e-3): D/2S 4.0。
  - c20cell_dual `Uz`・`roUz`・`CHECKPOINT/roUzN`: roUz の最大 base 1.7e-8・4.3e-7・5.7e-9、new 4.2e-7・5.8e-9・5.7e-9 (roUx は 156)。D/2S 1.68。
- メモリ: SERN g3 ホスト VmHWM 5157 → 2630 MiB (2816 → 1436 B/節点、−49 %)、GPU 不変。ローカル縮小格子 2 点の傾き 2737 → 1371 B/節点 (§6 の ≤ 1500 は満たす。g4 は入力消失で g3 単点 + ローカル傾きで判定)。

## 期待値と出典
- plan §6 (2026-10-07 事前登録): (a) step 0 の残差行は全列ビット一致 (1 ulp の 2 値の列は例外)、(b) 各データセットで「base 3 回・new 3 回の同ビルド内ペア差 (計 6 対) の最大 S の 2 倍以内」(D = base 対 new のペア差の最大)。[[repeat-range-as-limit-rule-trap]] を意識して係数 2・両側プール。

## 仮説 (呼び出し側)
- 差は全て run 間の非決定性 (atomicAdd 集計) で、変更起因の系統差は無い。根拠: デバイス側を変えていない設計、初期出力と checkpoint のビット一致、保存量・原始量は全構成で D ≤ 2S、FAIL 量は base 内で桁が動く重い裾の診断量。
- ただし (a) の前提 (1 ulp の 2 値) は事前に誤っていたので、(a) は判別力を持たなかった。

## 問い
1. この FAIL をどう扱うか (登録判定は FAIL のまま記録し、追加の判別試験で結論を出す、が呼び出し側の案)。
2. 判別試験の案と、その事前登録すべき合格条件:
   (i) **初期化直後の全デバイス配列のハッシュ**: env で有効な小さな診断 (例 `FORGE_DEVICE_HASH=1`: 初期化の終わりに全 `c_d`/`p_d`/`bvar_d`/格子マップを D2H して名前ごとのハッシュを出力) を base コミットと new コミットの両方に載せた 2 バイナリで、全構成のハッシュが一致すること (計算は決定的なので一致すべき)。一致すれば、以後の差は非決定性だけ、と言えるか。
   (ii) FAIL した構成・量だけ反復を増やす (例 各 10 回) 統計比較 (順位和検定や max の分布比較)。必要か、(i) で十分か。
   (iii) 1 step 目の組立直後のデバイス状態のハッシュ (atomicAdd の前) は取れるか / 要るか。
3. (a) の基準をどう作り直すべきだったか (今後の回帰の作法として plan と procedures に残すため)。
