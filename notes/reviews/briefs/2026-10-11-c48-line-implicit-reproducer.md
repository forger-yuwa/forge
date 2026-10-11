# 諮問ブリーフ: case/45 の float の偏りを小さい case/48 で再現させる試験 (line implicit を足す) (2026-10-11)

- 依頼者: 主セッション (エスカレーション条件 1・4: §6 を足す、承認済みの手順に無い run の前)
- plan: `plans/active/architecture-float-state-double-geometry.md` (§6.21〜§6.24 case/48、§6.26〜§6.36 case/45)
- ユーザ決定 (2026-10-11): case/45 は 1 本 1〜2 時間で切り分けが遅すぎる。偏りが再現する小さいケースを先に作る。§6.36 (case/45 で 10 本) は保留。

## 0. 問い
1. 下の試験の組み方・判定に穴はないか。「再現した」の基準をどう置くか (結果を見る前に決める)。
2. 再現しなかったときに、次に足す要素 (軸対称・TP 多成分) へ進む順序の考え方。

## 1. 観測事実
- case/45 (偏りあり): 共通の Q32 の起点から 30,000 step で、FP64 は θ_r が ≤ 0.12 % しか動かず、float は +0.9〜1.9 % 動く (§6.27)。SST の更新の停止 (§6.29)・block の系の組み立ての double 化 (§6.33)・平均流の FP64 での累積 (§6.35) は、どれも判別不能で、大幅な抑制は観測されなかった。設定: node、**軸対称**、SLAU、**TP 多成分** (thermalMethod 2、2 種)、SST、定常陰解法 blockDPLUR + **lineImplicit 1 + lineDtDirectional 1 (上限 50) + implicitThermalJacobian 5**、**CFL 4 / cfl_pseudo 4**、implicitRelax 0.7、nStepInner 5、limiter 2、57 万節点。
- case/48 (平板 M4・冷却壁 300 K、偏りなし): §6.21〜§6.24 で float ×2・FP64 ×2 を同じ入力から 48,000 step。窓 (24,000〜48,000) の平均は float と FP64 で 7e-5 以内、float の反復の揺れ 1e-4 (§6.22〜§6.24、座標の訂正後)。設定: node、2D 平面、SLAU、CPG 単成分、SST、blockDPLUR 1 (**line implicit なし**)、CFL 2 / cfl_pseudo 2、implicitRelax 0.7、nStepInner 5、limiter 2、scalarGradient lsq、8.9 万節点。
- つまり case/45 と case/48 は「line implicit 系の設定・CFL」「軸対称」「TP 多成分」の 3 つが重なって違う。
- line implicit の線は実行時に壁の bcond から作る (`mesh/mesh.cpp:1200〜1312`) ので case/48 にも使える。

## 2. 呼び出し側の案 (未検証)
- §6.21 と同じバイナリ (`~/forge-fgeom7-{f32,fp64}`)・格子 (`case/45/_conv614/c48/new64/m.h5`)・入力 (`_prec/mesh_init.h5`)・設定 (`run_0057` の設定) に、case/45 の陰解法の設定 `lineImplicit: 1, lineDtDirectional: 1, lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5` と `cfl: 4.0, cfl_pseudo: 4.0` を足す。腕は float ×2 (`run_0061`・`0062`)・FP64 ×2 (`run_0063`・`0064`)、48,000 step・2,000 ごと、AWS A で 4 本同時。
- 共通の入力を精度で揃えるため、§6.26 と同じく保存量を float32 に丸めて double に戻した入力にする案 (§6.21 の入力が既に float32 で表せるかは未確認)。
- 判定は §6.21 の道具 (`c48_prec_series.py`・`c48_prec_judge.py`、窓 24,000〜48,000、τ = 0.0005 (θ・δ*)・0.001 (Cf・q_w・CD・HF)、D − E > τ で「許容外の精度差を検出」) をそのまま使う。「再現した」= どれかの量で D − E > τ、かつ「0.5 % 以上の明瞭な差」(D − E ≥ 0.005) が θ か q_w にある。
- 再現しなければ、line implicit 系の設定・CFL は単独では偏りを生まないと読み、次は case/45 の粗い格子 (`problem_d155_ns_coarse.yaml`) か、軸対称・TP 多成分を足した小さいケース。

## 3. 読んでよいファイル
- float plan §6.21〜§6.24、§6.26〜§6.36
- `case/48.flat_plate_cooled_m4/c48_prec.sh`・`tools/c48_prec_judge.py`・`tools/c48_prec_series.py`
- `solver_density_cuda/mesh/mesh.cpp` (1200〜1320)
