# 諮問: 冷却壁の NS の対で、float32 の格子座標では第一セルが表せない — FP64 のビルドに替えてよいか

日付 2026-10-08。諮問先 codex (diagnose、`~/.config/forge/diagnose-backend` = codex)。
エスカレーション条件 4 (事前登録した手順 plan §6 V-c45 の「forge の sha256 が run_0167 と同じ」から外れる再実行の前)。
plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.7・§5.1 #13・#14・§6 V-c45 (全文を読むこと)。作業ツリー `/home/sano/work/forge-integ-1005` (commit 656bd37f)。
ユーザ: 「(第 1 層の実装と冷却壁の NS を) やってって。残念ながら冷却壁は y+ 上がっちゃうから気をつけてね」(2026-10-08)。

## 観測事実

1. **冷却の y1+ の倍率**: 生産の断熱 NS (`case/45.isobutane_m6_d155/run_0179_ns_n012_N2_ext`、step 20000) の壁の y1+ 分布
   (`check_wall_resolution.py --profile-csv`、`_band_ab/cold_pair/y1p_ad_wall.csv`) に、局所の倍率
   √((T_w,ad/300)·C_f 比)·μ(T_w,ad)/μ(300) (C_f 比は CONTUR のエンタルピー形、μ は NS と同じ CEA 混合気) を掛けると 8.5〜9.4 倍。
   断熱の y1+ は平均 0.277・最大 1.49 なので、同じ格子の 300 K は入口直管で最大 14、スロートで 10、x 5〜20 で 2.8、y1+ > 1 の面積 42 % の見込み。
2. **格子の設計** (`cold_pair_mesh.py`、`mesh2d` に表の格子を追加、既存格子はビット同一を確認):
   第一セル厚 / 局所半径を y1+_cold ≤ 0.8 になるように表にすると 3.44e-7 (スロート) 〜 1.53e-5 (出口)。ni 4496 × nj 121 (54 万節点)。
3. **投入の前提で止まった**: `run_cold_pair.sh` の準備で `check_mesh_quality --ar-max 5000` が FAIL (AR 最大 6014、5000 超 701 セル 0.13 %)。
   msh から倍精度で測った AR は最大 5488 で、場所は縮流部 (x −9〜−1、壁の傾きで壁沿いの辺が長くなる分を設計で見ていなかった)。
   作りかけの run dir は消した (forge は起動していない)。
4. **格子座標は float32**: 変換後の `nozzle.h5` の `MESH/COORD`・`CELLS/centCoords`・`CELLS/volume`・`PLANES/surfVect`・`VALUE/wall_dist` が float32
   (`flowFormat.hpp` の `geom_float = float`、生産のバイナリ sha256 6b47811b…)。品質ツールの AR (6014) と倍精度の AR (5488) の差もこれと合う。
   第一セル厚 / 局所半径 3.44e-7 は、半径座標の float32 の相対精度 (6e-8〜1.2e-7) で 3〜6 ulp。生産の格子 (4.5e-6) は約 40 ulp。
   寸法を変えても比は変わらない (相対精度の問題)。
5. **FP64 のビルド**: 生産のバイナリと同じソース (`~/forge-wallfit-bin`、HEAD e2696d8f0) を複製し、`flowFormat.hpp` の 4 つの typedef だけを double にして
   ビルド中 (`~/forge-wallfit-bin-fp64`、同じ cmake の設定)。別の古い commit の全域 double ビルド (`~/forge56-double`) の実績はある。
   FP32 の生産の格子は 5000 step 52 s (19.4 万節点、3 本並列のとき)。FP64 の速度は未計測。

## 期待値と出典

- 冷却壁 (T_w/T_aw 0.2) の低 Re SST は y1+ ≤ 1 が要る (case/48 平板: δ*/θ は y1+ ≤ 3.6 で 2 % 内、q_w は y1+ ≤ 1 で 1.5 % 内; plan §5.1 #2・#9)。
- V-c45 で判別したい差は、温度形とエンタルピー形の 300 K/断熱の δ_r 比の差で最大 2.8 % (plan §5.1 #10a)。

## 選択肢

- **(A) FP64 ビルドで 2 本とも回す**: 同じ格子・同じ FP64 バイナリで断熱と 300 K を回す。判定は 2 本の比なので、生産の FP32 との違いは比の判定に効かない
  (生産との差は「格子の感度」と同じく記録のみ)。事前登録の「forge の sha256 が run_0167 と同じ」を「FP64 ビルドの sha256 (ソースは生産と同じ、typedef 4 行だけ違う)」に替える。
  格子は縮流部の AR を見込みに入れて作り直す (壁沿いの辺 = dx·√(1 + r′²)、AR の目標を下げる)。
- **(B) FP32 のまま第一セル厚 / 局所半径 ≥ 5e-6 に制限**: スロートの y1+ は 12 前後、x 5〜20 で 1.7〜3、試験部で ≤ 1 の見込み。判別したい差 (≤ 2.8 %) より
  格子の誤差が大きくなりうる。
- **(C) 座標だけ double (geom_float = double、flow_float = float)**: コンパイルできるか未確認。

## 問い

1. (A) に替えてよいか。替えるなら、事前登録 (§6 V-c45) をどう書き直すべきか (バイナリの照合、断熱の FP64 と生産の FP32 の比較を何に使うか)。
2. FP64 で判定の前に確かめるべきこと (例: 同じ格子の断熱を FP32 と FP64 で比べる、座標の量子化の影響の見積もり)。
3. 格子の設計の直し方 (縮流部の AR、目標の y1+ 0.8 は妥当か、入口直管の角をどう扱うか)。
4. (B) や (C) を選ぶべき理由があるか。
5. 生産の断熱の格子 (第一セル 4.5e-6 で約 40 ulp) の float32 の量子化は、これまでの生産の判定 (δ_E/δ_C、壁解像) に効いていないと言えるか (記録として)。

## 読んでよいもの

- 上記 plan、`case/45.isobutane_m6_d155/cold_pair_mesh.py`・`cold_pair.py`・`run_cold_pair.sh`、`design/forge_design/meshing/mesh2d.py`
- `solver_density_cuda/flowFormat.hpp`、`solver_density_cuda/tools/check_mesh_quality.py`
- `notes/reviews/2026-10-08-contur-property-temperature-diagnose.md`
