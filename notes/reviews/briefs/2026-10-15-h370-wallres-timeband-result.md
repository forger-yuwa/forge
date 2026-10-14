# 諮問ブリーフ: V-h370-diag の壁の解像 (#36o) と帯の時間方向の A/B (#36p) の結果の解釈と次の一手 (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 7: 結果の解釈の確定の前、条件 1・4: 次の修正・登録の前)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6「帯の時間方向の扱いの A/B の事前登録」とその結果、「壁の解像の切り分けの結果」、§5.1 #36n〜#36p
- 前の諮問: `notes/reviews/2026-10-14-h370-fixedband-result-diagnose.md` (#36n の解釈、#36p の登録、#36o の先行)
- commit: 7584da7d (#36o)、6f49b928 (#36p)

## 観測事実

結果の数値は plan §6 の 2 つの結果の項に書いた (記録 `_band_ab/half_inlet/h370_diag_time_band.json`・`h370_wallres_split.json`)。要点:

- #36p: D_time = 0.004316 mm (≤ 0.01、結果 B)。腕 A の再現・先頭の時刻の再現はビット同一。2 窓の最大−最小は腕 A の δ_E で 0.00816・0.00817 mm、窓の平均の差 0.00054 mm。`check_quasisteady --series-csv` は両腕 ALL STEADY (4 列)。腕 A の 0〜10 r_t の時間の変化 (≤ 0.0082 mm) は帯の選び直しが運んでいる (腕 B では 1e−5 以下)。梯子は 3〜10 r_t で交差の区間の変化・飛びがある。
- #36o:
  - 受け渡しの関門 `deltastar_loop.wall_resolution_verdict` (`deltastar_loop.py:652`) は `--groups` を渡さず、生成器が出口の境界の出力を書く (`runner_wt.py:284`) ので、関門の y₁⁺ の判定は**常に判定不能**。壁だけで判定すると FAIL (1 超えの面積 2.9 % > 2 %)。
  - y₁⁺ > 1 は x −12.5〜−5.36 r_t の 1 区間 (直管と縮流部の入口側、最大 2.77 は直管の中 x −10.53)。実際の y₁ は設計の第一層と一致。下流の y₁⁺ ≤ 1 の隣と比べ y₁ ×2.10・τ_t ×0.74。
  - 形 C の 370 K の run (格子 g0c、未収束) の同じ区間は τ_t が 1/11.4。
- Euler の資格は予備の資格 FAIL のまま (残差 2.0〜2.5 桁で停滞、CFL で床が動く: plan §6 の CFL の A/B の結果)。

## 期待値と出典

- §6 V-h370: 受け渡しには Euler の資格・全使用域の σ_ref,use ≤ 0.01 mm・2 窓・局所 y₁⁺ ≤ 1 (超過の面積の許容 2 %) などが要る。
- 格子の第一層の配点は `mesh.wall_first_up_x0/x1` などで決まる (`wt_card.py` の「第 1 層の配点」の候補の選び方、`design/forge_design/geometry/mesh2d.py` の `_first_frac_profile`)。

## 呼び出し側の仮説・案

1. #36p は登録どおり結果 B と読む。#36n (帯を共通にすれば参照の入れ替えの差 0.0011 mm) と合わせ、「δ_E の参照・時間の変化の大半は帯の選び直しが運ぶが、この窓ではどちらも 0.01 mm 以内」まで言える。
2. 関門の欠陥 (常に判定不能) は 1 ファイルのバグ修正として直す: `wall_resolution_verdict` に壁の群だけを渡す (`nozzle_report` と同じ選び方)。負例の試験 (出口の群だけ・壁の群に応力が無い) を付ける。
3. 格子: 入口側 (x < −5 r_t) の第一層を下流と同じ比で薄くする (配点の上流側の値を変える)。y₁⁺ の目安の見積もり (`wt_card` の第 1 層の配点の候補) が H の入口側の大きなせん断を見積もれていないかも確かめる。
4. Euler の資格: 独立に収束を確かめた参照 (別登録) は、格子の修正の後で。

## 問い

1. 解釈 1 の言える範囲。
2. 案 2 (関門の修正) の採否と、これで過去の受け渡しの判定に影響があるか (今まで「判定不能」で止まっていた / 通っていた経路はあるか)。
3. 案 3 の格子の修正の前に何を確かめるか (第一層の配点の設計の目安と実際の τ_t の差、H と C の違い)。修正後に回す最小の run と合格条件。
4. 次の一手の順序 (関門の修正・格子・Euler の資格・帯の頑健化)。

## 読んでよいファイル

- plan §5.1 #36a〜#36p、§6 V-h370 系の登録と結果
- `case/45.isobutane_m6_d155/h370_wallres_split.py`・`h370_judge.py`、`case/45.isobutane_m6_d155/_band_ab/half_inlet/h370_wallres_split.json`・`h370_diag_time_band.json`
- `design/forge_design/feedback/deltastar_loop.py` (`wall_resolution_verdict`)、`design/forge_design/report/nozzle_report.py` (`wall_resolution`)、`solver_density_cuda/tools/check_wall_resolution.py`
- `design/forge_design/geometry/mesh2d.py`、`design/forge_design/intake/wt_card.py` (第 1 層の配点)、`case/45.isobutane_m6_d155/` の問題 YAML (370 K・H)
- `notes/reviews/2026-10-14-h370-fixedband-result-diagnose.md`
