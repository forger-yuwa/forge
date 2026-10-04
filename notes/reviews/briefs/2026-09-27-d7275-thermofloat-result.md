# 諮問ブリーフ: D-7275 試験 26 の thermoFloat A/B (T4-0b-TF) の解釈と、自由流域の残差床の次の一手 (AGENTS 条件 7・4)

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` §4.12、§6 G15、§5.1 #61 の末尾
- 事前登録と結果: `case/60.flatplate_d7275_m7/acceptance.json` の **T4-0b-TF** (`outcome` に数値を全部記載)、T4-0b-H
- 台帳: `case/60.flatplate_d7275_m7/README.md` (run_0015〜0018)
- 前回の諮問: `notes/reviews/2026-09-27-d7275-height-ab-result-diagnose.md` (この A/B はその提案どおり)
- 集計: `case/60.flatplate_d7275_m7/tools/tf_ab.py`
- 変換器: `solver_density_cuda/` の convertGmshToForge (メッシュ h5 の書き出し型)、FP64 ビルドは `flowFormat.hpp` の typedef (flow_float・geom_float とも double) だけを変えたもの

## 1. 観測事実 (AWS。数値は acceptance.json の T4-0b-TF outcome)

- A = thermoFloat 1 (`run_0015_t26_tf1` → 延長 `run_0017_t26_tf1_ext`)、B = thermoFloat 0 (`run_0016_t26_tf0` → `run_0018_t26_tf0_ext`)。run_0014 最終場から restart_field (12 量ビット一致)、各 5k + 延長 5k。延長は 5k で k/ω の窓間中央値が 5 % 超動いたため (事前登録どおり)。
- 10k 時点の窓間中央値変化: 流れ・化学種 ≤ 1.1 %、B の roOmega −7.7 %・roK −6.9 % (まだ動く)。
- **自由流域 (y 0.03–0.79) の √Σres² B/A = 0.91–1.05 (全成分)**。
- **第一内部列 (y = 3 µm)**: roUy 0.028、roe 0.042、ro 0.10、roY 0.10–0.20、roUx 0.48、roOmega 0.46、roK 0.84。
- ω 収支 (A の res_roOmega 上位 5 点): trans ≈ dest ≈ 6.0–6.5e8、和 = res_roOmega が A ±0.8–1.3e2 → B ±1.5e-3〜3.1。
- 全域 RMS 履歴 B/A (末尾 1k): ro 0.96、roUx 0.99、roUy 1.05、roe 0.98、roY 0.96、roK 0.77、roOmega 0.39。両腕 `NOT CONVERGED (stalled/plateau)`。
- 固定節点の振幅 (100 step 毎 11 枚): 壁列の roOmega は一部の点で B が 1/20–1/30、他は同程度。自由流 3 点は ρ・P・k・ω とも両腕 ~1e-7。
- 比較量: 比較域 q_w の B/A−1 最大 9.2e-7、R_A 平均 1.313。
- **新しい観測**: B の最終場で、自由流域の |res_ro| / (ρ U √vol) は中央値 1.2e-8、90 % 点 5.3e-8、最大 3.3e-7。**変換済みメッシュ h5 の MESH/COORD と CELLS/volume は float32**。

## 2. 期待値と出典

- 事前登録 (T4-0b-TF decision): B の対象残差が A の 1/10 以下かつ局所振幅減少 → 支持 / 0.8–1.2 倍 → 棄却 / 流れだけ改善して ω が残れば帰属を流れに限定 / 中間は保留。診断条件で G15 の代替にしない。
- 過去の同型の事例 (memory): 軸対称の r 重み閉包が高 AR 壁 CV で float32 桁落ち、双対 CV 重心・体積の float32 桁落ち (いずれも別ケース・別修正)。

## 5. 仮説 (未確認)

- 呼び出し側の読み: 自由流域の床について「float 熱力学経路が主因」説は棄却。第一内部列の流れ (roUy・roe) には寄与が大きい。ω は判別保留。
- H-geo: 自由流域の床 (相対 1e-8) は float32 で保存された幾何 (面ベクトル・体積) の閉包誤差による。一様流で Σ_f S_f ≠ 0 が float32 の丸めで残り、残差 ~ ε_f32 × 流束。未確認 — 面ベクトルの閉包を数えていない。

## 6. 聞きたいこと

1. 上の読み (自由流域は棄却・第一内部列の流れは寄与大・ω は保留) を確定してよいか。
2. 自由流域の床の次の一手を 1 つ。H-geo なら、どう判別するか (例: メッシュ h5 の幾何を double で書き出して同じ場から継続する A/B — 変換器の変更が要るか、後処理で double の幾何を再計算して差し替えられるか)。
3. G15 の扱い: 比較量は 8 回の変更で max|Δq/q| ≤ 5e-6 しか動かない。どこまで追うべきか、どこで「残差床は比較量に効かない既知の数値床」として G15 を条件付きで扱う提案をユーザに出すべきか。
