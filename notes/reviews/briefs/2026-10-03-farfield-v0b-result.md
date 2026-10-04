# 諮問: farfield V0 後半 (多 step 再現性) が事前規則で FAIL — 解釈と次の一手 (2026-10-03)

関連 plan: `plans/active/boundary-node-farfield-characteristic.md` §5.1 #5a (事前登録)、§6 V0。result 段レビュー: `notes/reviews/2026-10-03-boundary-node-farfield-characteristic-result.md` (GO-with-changes、M3 の一次証拠回収で V0 後半の欠落が判明)。
証拠: `notes/investigations/2026-10-03-farfield-evidence/c46_V0B_REPRO_20261003.txt`、索引 `INDEX.tsv`。

## 観測事実
- V0 前半 (済 2026-09-29): farfield を含まない SERN 3D 構成で、同一初期状態からの初回 massflux・状態ダンプが変更前バイナリとビット一致 (`run_0991_ff_v0_old/new`)。
- V0 後半 (本日、事前登録どおり): run_0971 (g3、farfield なし、SST・2 種 TP) の設定と最終場から、旧 = farfield 実装前 `~/sglsq/forge_2fa3826c` (sha dcd17eb3) × 3、新 = farfield 実装後 build-ff (sha b0240cd7、種 DB 取り込み前) × 3、各 200 step。保存量の最終場で変数ごとに S_old = 旧同士 3 組の max_node|ΔQ|、D = 新×旧 9 組の max_node|ΔQ|。
- 結果: 全変数で D > S_old (FAIL)。D/S_old: ro 1.08、roUx 1.84、roUy 1.14、roUz 1.30、roe 1.08、roK 1.24、roOmega 1.42、roY0 1.08、roY1 1.82。新同士の S_new も同じ桁 (roUx S_old 0.278 / S_new 0.384 / D 0.513、roY1 1.54e-4 / 2.04e-4 / 2.80e-4、roOmega 4.46e5 / 7.51e5 / 6.33e5)。全値有限。
- 規則の性質: D は 9 組の最大、S_old は 3 組の最大。同じ分布から取っても D の方が大きく出やすい。
- 旧 → 新の間のコード差は farfield の追加 (非 farfield 面は既存経路、初回評価はビット一致) と、その後の V2 用の小修正 (評価器・プローブ桁など、出力専用) のはず。node 残差は atomicAdd で、同一バイナリでも 1 step でビット再現しない (既知)。

## 期待値と出典
- §6 V0: 「更新後保存量は旧バイナリ 3 回反復の再現性幅以内」。事前規則は本日 §5.1 #5a に「全変数で D ≤ S_old」と具体化した (結果を見る前)。

## 仮説
1. 新旧の違いは無く、FAIL は組数の違う max 比較の統計的な偏り (新同士も旧同士を上回っている)。
2. 非 farfield 経路に、初回評価では現れない小さな差 (更新・陰解法・別カーネル) が入っている。

## 問い
1. この FAIL を「同等性が示されない」と読むべきか、規則の不備 (組数の偏り) と読むべきか。結果を見た後で規則を変えずに判別するには、どんな追加試験が要るか (例: 同じ 6 本を旧 3 + 旧 3 の 2 群に分けた対照で D/S の分布を作る、本数を増やす、初回以降の各 step でビット比較を続ける — 1 step 目から atomicAdd で割れるか)。
2. 2 の仮説を潰す最小の試験は何か (例: 非 farfield 構成で旧・新の 1 step 後の res_* を FORGE_DUMP_MASSFLUX で照合、決定的な比較)。
3. これが未決のまま、§2 の限定受理で accepted に進めてよいか。
