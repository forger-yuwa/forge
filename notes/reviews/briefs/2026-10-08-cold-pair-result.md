# 諮問: 冷却壁の NS の対 (V-c45) の 1 回目 — ゲート不成立の扱いと、CONTUR が冷却の効果を外す件

日付 2026-10-08。諮問先 codex (diagnose)。エスカレーション条件 3 (事前登録の判定が判定不能) と 7 (結果の解釈を確定する前)。
plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` の §4.7、§5.1 #10〜#20、§6 V-c45 (「延長の決め方」と「1 回目の結果」を含む。全文を読むこと)。
作業ツリー `/home/sano/work/forge-integ-1005` (commit 75a8f28b 以降)。判定のコード `case/45.isobutane_m6_d155/cold_pair.py` (`judge`)、出力 `_band_ab/cold_pair/V_c45.json`・`gates_aws.json`。

## 観測事実

- run: `case/45.isobutane_m6_d155/run_0181_ns_coldmesh_ad` (断熱) と `run_0182_ns_coldmesh_tw300` (300 K)。
  - 格子は冷却壁用 (ni 4719 × nj 121、近壁は壁法線、msh 17 桁)、FP64 のビルド、生産の壁。
  - IC は run_0179 の cross-mesh。段階起動 full → 本段 100000 step。
  - 投入の前提はすべて成立 (品質 PASS、第一層の誤差 0、高 AR のスキューセル 0、壁の差 0)。
- **ゲート 1 (NaN)**: 合格 (全段の残差と 21 枚の場)。
- **ゲート 2 (収束)**: 両 run とも `NOT CONVERGED (stalled/plateau)`、RISING なし。残差の床 (本段の最後の 5000 step の中央値) は生産の run_0179 (FP32・別格子 19.4 万節点) を上回る。
  - 断熱: rms_ro 7.0e-6 (run_0179 は 3.4e-7)、roUx 3.3e-3 (2.8e-4)、roe 10.2 (0.54)、roK 2.9e-5 (6.5e-5、これだけ下)、roOmega 573 (29.8)。
  - 300 K: rms_ro 9.1e-6、roUx 3.3e-3、roe 10.5、roK 8.5e-4、roOmega 1230。
  - 末尾 40k step の推移: 断熱は ro・roUx・roe が横ばい (−1 %)、roK・roOmega は下がり続ける。
    300 K は ro 1.11e-5 → 9.1e-6、roK 1.6e-2 → 8.5e-4 (1/19)、roOmega 2736 → 771 → 1230 (下がってから上がる)。
  - 同じ生産の格子では FP64 の方が FP32 より床が低い (400 step の比較で rms_ro 3.4e-7 → 7.4e-8)。今回の床の高さは格子の違い (節点 2.9 倍、第一層が最大 1/13) から来ている可能性があるが、確かめていない。
- **ゲート 3 (準定常、classify、5 枚全部、drift・osc 0.1 %)**:
  - 断熱の δ_E と壁温は全点 STEADY。
  - 300 K の δ_E は全点 DRIFTING。5000 step ごとに単調に増え、増分は前の 0.65〜0.85 倍。x = 94 で 43.80 → 44.08 mm、等比外挿の残りは約 +0.17 mm。
  - Q_w は −12.34 → −11.49 MW で DRIFTING (増分の比 0.91〜0.92、外挿の残りは約 2 MW)。
- **ゲート 4 (壁解像、面積、300 K)**: 全壁 1.8 %、試験部 0 %、[−1, 40) 0 % は合格。縮流部 (入口の角を除く) は y1+ > 1 が 24 % で上限 10 % を超える (最大 1.23、x ≈ −11.4)。断熱は全域 0 %。
- **記録 (判定に使わない)**: R_NS = δ_E(300 K) / δ_E(断熱) は x = 40 で 0.725、60 で 0.758、80 で 0.779、94 で 0.787。冷却で δ_E が 21〜27 % 薄くなる。
  - CONTUR の予測 (同じ壁、生産の k_f) は、温度形 1.049 → 1.012、エンタルピー形 1.020 → 0.986、contur_v2 (熱閉包 + 粘性 + 加速の項) 1.012 → 0.977。
  - e (|R_m/R_NS − 1| の最大) は温度形 0.45、エンタルピー形 0.40 (k_f 1 でも同じ)。区間の判定の形式上は「エンタルピー形を支持」だが、両方とも冷却の効果をほぼ再現しない。
  - 不確かさ U は 0.8〜1.5 %。格子と精度を替えた感度 (この格子の断熱の δ_E / 生産の δ_E) は試験部の平均 +0.32 %、最大 0.52 %。
- **参考**:
  - case/44 (空気 TP、Tt 1060 K) では断熱 → 300 K で出口コア M が +0.35〜0.48 % 上がり、「積分法の初期壁は冷却効果を過小評価する」と記録している (plan §5.1 #4)。
  - case/48 平板 (M4.19、300 K) では、CONTUR の δ* が NS より 3.5〜6 % 小さく、θ は ±0 % だった。
  - CONTUR は冷却 (300 K) で θ が約 2 倍、H が 13.9 → 7.2 になり、δ* がほぼ変わらない (θ と H の変化が打ち消し合う)。NS の θ・H・C_f は今回まだ抽出していない。

## 実施済みの操作

- 判定は事前登録どおり「判定不能 (ゲート不成立)」とした。
- 事前登録の「延長の決め方」に従い、DRIFTING の 300 K の腕だけを延長中: `run_0183_ns_coldmesh_tw300_ext` (run_0182 の res_100000 から restart_field、ビット一致、100000 step、5000 ごと)。判定窓は延長の最後の 5 枚。
- 原因は plan に書いていない。

## 仮説 (私の案)

1. ゲート 2 の「残差の床が生産以下」は、格子が違う run どうしの絶対値の比較で、基準として不適切だった。
   - 代案: 同じ格子の中で、各保存量の床が末尾 40k step で横ばいか下がっていること。RISING がないこと。
   - それに加え、比べる量 (δ_E・壁温・Q_w) が STEADY であること。
   - ただし、結果を見た後に基準を替えることになるので、替えてよいか、どう記録するかを問いたい。
2. ゲート 4 の縮流部の超過 (最大 1.23) は、試験部の δ_E の比較への影響は小さいと見るが、未確認。
3. CONTUR が冷却の効果を外す原因の候補 (確かめていない、順不同):
   - (a) 冷却で θ が 2 倍になる摩擦の閉包 (C_f の圧縮性変換が冷却壁で過大)
   - (b) 形状係数 H の関係 (速度分布 N と温度分布 Walz/Crocco の組み合わせ)
   - (c) 縮流部・スロートの強い冷却と加速の履歴
   - (d) NS 側の δ_E の抽出 (等温壁の符号付き、core-matched Euler の基準) が冷却壁で別のものを測っている
   - 切り分けには NS の θ・H・C_f を両方の壁で抽出して CONTUR と比べるのが最初の一手だと考える。

## 問い

1. ゲート 2・4 の不成立をどう扱うべきか。基準の書き直しは許されるか (許されるなら、どう記録するか)。延長の run だけで判定し直してよいか。
2. 延長 (run_0183) で 300 K の δ_E が STEADY になったとき、V-c45 の「熱閉包の形を選ぶ」という問い自体に意味は残るか。両方とも e ≈ 0.4 で外れているので、判定の問いを「CONTUR は冷却の効果を当てるか」に替えるべきか。
3. CONTUR の冷却の効果の外れの切り分けで、最小の A/B は何か (NS 側の θ・H・C_f の抽出、δ_E の抽出の別定義、CONTUR の閉包の差し替えのうち何を先に)。
4. 設計の手順への含意: 壁温を変えたら NS の δ* 反復が必須、と書いてよいか。contur_v2 を既定にする価値は残るか。

## 読んでよいもの

- 上記 plan、`case/45.isobutane_m6_d155/cold_pair.py`・`cold_pair_mesh.py`・`delta_contur_compare.py` (predict・hform)、`design/forge_design/feedback/deltastar_integral.py`
- `case/45.isobutane_m6_d155/_band_ab/cold_pair/V_c45.json`・`gates_aws.json`、`_band_ab/delta_contur/cooling_ratio_predictions.json`
- `design/forge_design/metrics/deltastar.py` (`deltastar_from_core_matched_euler`)、`design/forge_design/feedback/deltastar_loop.py` (`extract_and_merge`)
- `notes/reviews/2026-10-08-contur-property-temperature-diagnose.md`、`notes/reviews/2026-10-08-tooling-nozzle-isothermal-wall-chain-plan.md`
