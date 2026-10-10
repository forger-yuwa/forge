# 諮問ブリーフ: case/45 ライン陰解法の総時間 — ここまでの結果の解釈と速度向上策の優先順 (2026-10-10)

- 依頼者: ユーザ「ここまでの結果と、速度向上策について、diagnostician に諮ってくれないか」
- plan: `/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md` (§4・§5.1・§6.6〜§6.16)
- 関連 plan: `/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md`、`/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md`
- コードの版: ブランチ `feature/nozzle-wall-fit-and-pipeline` (worktree `/home/sano/work/forge-integ-1005`)。実行バイナリ lineL = f9be0c4f + `typedef double` (sha256 e10a195d…、AWS `~/forge-linevisc-fp64`)
- run の実体は AWS (`~/forge-wallfit/case/45.isobutane_m6_d155/`) にある。手元にあるのは判定の記録と、ParaView 用の `case/45.isobutane_m6_d155/run_0353_m9_L5/res_115000.{h5,xmf}` だけ

## 0. 問い (判断してほしいこと)

1. **ここまでの結果の解釈に誤り・飛躍はないか**。特に次の 4 点。
   - (a) ラインと point の止まる値の差 (θ_r 約 0.1 %)
   - (b) 緩和 > 0.7 で揺れが大きくなること
   - (c) 終わりを「水準」だけで決める基準の妥当性
   - (d) 粘性入りのライン (L5) と値 0 のライン (L0) が総時間で区別できないこと
2. **速度向上策の優先順と、総時間で比べる試験の設計**。候補は §4。1 本の総時間 run が専有換算で約 1.1 h かかる。何を組み合わせ、何本回し、どこから出発すべきか。
3. **総時間の評価を安くする方法**はあるか (途中の場から出発する、短い代理の指標を使う、など)。ユーザは判断に要らない長い run を嫌う (「3h は長すぎ。なんでその計算かける？」)。
4. **格子** (壁→軸を 1 本の等比数列で張っている) について、速度・収束の観点から指摘すべきことはあるか。

## 1. 観測事実

### 1.1 ケース

- case/45 イソブタン燃焼ガス M6 ノズル (d155)
- 2D 軸対称、node、定常 (`unsteady 0`)
- 等温壁 300 K (`wall_isothermal Ts 300`)
- 入口 `inlet_Pressure` Pt 5.5 MPa・Tt 1600 K、出口 `outlet_statPress` Ps 2237 Pa
- TP 2 成分 (`thermalMethod 2`)、SST 低 Re (`wallTreatmentSST 0`、`dilatationCorrection 2`、`katoLaunder 1`)
- `convMethod 1`、`limiter 2` (参照値つき)、`blockDPLUR 1`、`nStepInner 5`、`cfl_pseudo 4`、`implicitRelax 0.7`
- 格子: 4719 列 × 121 節点 = 570999 節点、品質 PASS (AR ≤ 5000 緩和、最大 4573、スキュー 0.457)
- 出発点: `run_0183_ns_coldmesh_tw300_ext` の res_100000。この run 自体は NOT CONVERGED (stalled/plateau、`rms_roOmega` RISING)

### 1.2 ライン陰解法の構成と 1 step の内訳

ライン:
- 壁の節点を種に、壁法線へ伸ばして軸まで。4719 本 × 121 節点で全 CV を覆う
- 1 ライン 1 スレッドの block Thomas (5×5)

1 step の流れ:
- 分解 1 回 → 代入 5 回。sweep 間はライン外と Jacobi で結合する

単価 (`m9_unit.json`、専有、1000 step × 3 本):

| 構成 | 1 step [ms] |
| --- | --- |
| point | 17.83 |
| L0 (値 0 + 方向別 Δτ + 上限 50 + キー 5、LU) | 34.79 |
| L5 (粘性 値 3・マスク 5、上限なし) | 36.90 |
| LAYOUT2 (L0 の Thomas の配列を並べ替え + 連鎖の短縮。ビット一致、§6.12) | 32.14 |

ncu の内訳 (§6.6・§6.8・§6.12):
- 従来の Thomas は DRAM 74 %、読み込み命令 1 回あたり約 31 セクタ、SM あたり 2 warp (ライン数 4719 が少ない)
- LAYOUT2 で分解 4.58 ms・代入 1.60 ms × 5、Thomas の合計 12.6 ms/step

不採用 (§6.3・§6.5・§6.8):
- 逆行列の保存: +19 %
- Thomas を float に: η 1.9e-5 > 上限 1e-5
- LAYOUT1: −0.47 ms

### 1.3 総時間 (終わりまでの step 数 × 単価)

判定器 `m9_judge.py`、結果 `_band_ab/cold_pair/m9_judge.json`、plan §6.14:

| 腕 | 段 | E2 (参照 R からの差) を含む終わり | 総時間 |
| --- | --- | --- | --- |
| P (point だけ、事後の集計、古いバイナリの系列) | point 61.0 万 step | E2 −0.075/−0.064/−0.059/+0.085 % | 3.03 ± 0.05 h |
| L0 (事後の集計) | ライン 13.5 万 + point 3.5 万 | E2 +0.034/+0.032/+0.031/−0.048 % | 1.49 ± 0.08 h |
| L5 (事前登録、`run_0353_m9_L5` → `run_0354_m9_L5cut`) | ライン 11.5 万 + point 4.0 万 | E2 +0.021/+0.020/+0.019/−0.044 % | 1.38 ± 0.08 h |

- `check_convergence` は L5 の両段とも NOT CONVERGED。
  - ライン段: `rms_roUx`・`rms_roe`・`rms_roK` が plateau。
  - point 段: 全列が plateau、低下 0〜0.5 桁。
- 到達の窓の `check_quasisteady`: θ_r は TRANSIENT-UNSETTLED、Q_w は STEADY。

### 1.4 ユーザ決定 (§6.15、2026-10-10)

- 参照 R (`run_0263` の最後の出力、point 通算 68 万 step) は point が止まった場所であって、R(Q) = 0 の解ではない。
  - R 自体も NOT CONVERGED で、θ_r は +0.015 %/2 万 step で動き続けている。
  - したがって E2 を終わりの基準から外し、**水準だけ** で終える。水準 = |流量欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step。入口流量は約 98.8 kg/s。
- P は測り直さない。
- 水準だけで数え直した結果 (`levelreach.py`): L0 ライン 12.0 万 step (≈ 1.16 h)、L5 ライン 11.5 万 step (≈ 1.18 h) → 判別不能。
- 本番の手順の案: 値 0 のラインで水準まで (point の仕上げなし)。

### 1.5 止まった後の変動 (非公式の集計、4 万 step の窓、5000 step ごと 9 点)

θ_r (単位は %):

| 腕 | 揺れ (直線を引いた残りの標準偏差) | 窓内の幅 | 傾き (%/2 万 step) |
| --- | --- | --- | --- |
| point (`run_0263`) | 0.001 | 0.03〜0.04 | +0.014〜+0.018 |
| L0 ライン (`run_0252`) | 0.007〜0.013 | 0.03〜0.05 | +0.009〜+0.015 |
| L5 ライン (`run_0363_m9_L5line`、水準の後に追加 2.5〜6.5 万 step) | 0.007〜0.010 | 0.025〜0.04 | −0.005〜+0.009 |
| ライン → point (`run_0262`・`run_0354`) | 0.006〜0.016 | 0.07〜0.12 | −0.04〜−0.07 |
| 緩和 0.85 / 1.0 のライン (`run_0364`・`0366`・`0368`) | 0.08〜0.6 | 0.3〜1.6 | 振れる |

残差 (最後の 2 万 step の中央値、括弧内は log10(p95/p5)):

| 腕 | `rms_ro` | `rms_roe` | `rms_roK` | `rms_roUx` |
| --- | --- | --- | --- | --- |
| point | 3.7e-6 (0.02) | 6.5 (0.02) | 3.3e-4 | 3.1e-3 |
| ライン 緩和 0.7 | 8.8〜9.2e-6 (0.2〜0.3) | 8.1〜8.2 (0.1) | 5.7e-4 (0.3〜0.6) | 3.1e-3 |
| ライン 緩和 0.85 / 1.0 | 3.9〜5.3e-5 | 27〜48 | 4〜10e-3 | 6.5e-3〜1.9e-2 |

- `rms_roUx` は point とラインで同じ 3.1e-3 にある (床)。
- ラインの止まった値は、R から θ_r で約 +0.1 % 離れている。これはラインの揺れ (0.01 %) の 10 倍。
- point は R からラインの値の側へ +0.015 %/2 万 step で動いている。ラインから point に切り替えた腕は −0.04〜−0.07 %/2 万 step で下がっている途中。

### 1.6 緩和を上げる腕 (§6.13・§6.14)

- `run_0370_m9r_L5_r100` (L5・緩和 1.0) は step 159 で `ro` が非有限になり、detectNaN で停止。
- 他の 3 腕は 14〜15.5 万 step まで水準に入らず、ユーザ判断で打ち切った。緩和は 0.7 のまま。
- 初期の欠損の減りは速かった (L0 r0.85 の 4 万 step で 0.154 kg/s、r0.7 では 0.282)。

### 1.7 格子の壁法線方向の並び (`run_0353_m9_L5/res_115000.h5` から)

| x/r_t | r_w [r_t] | 第 1 層 | 比 Δy_{k+1}/Δy_k | δ99 (u ≥ 0.99 u_max) の節点 | δ99/r_w | AR > 100 / > 20 / > 5 の節点 |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 3.67 | 0.76 µm | 1.090 (全節点で一定) | 88 | 6 % | 45 / 64 / 80 |
| 40 | 8.38 | 5.6 µm | 1.079 (同) | 104 | 30 % | 29 / 50 / 68 |
| 94 | 10.10 | 11.8 µm | 1.073 (同) | 96 | 18 % | 20 / 43 / 63 |

- 壁から軸まで 1 本の等比数列。境界層の外の主流は 17〜33 節点しかない。
- x = 94 では節点 100→120 の 20 節点で半径の 75 % を覆い、軸際のセルの半径方向の幅は約 0.7 r_t。
- 軸方向の間隔はおよそ 0.02 r_t の見積もり (未確認)。したがって軸の近くでは、半径方向に長いセルになっている。

### 1.8 1 step の時間のふるい (§6.16、完了、`scr.sh`、判定 `_band_ab/cold_pair/scr_judge.json`)

- 照合 `run_0372_ml64_cmp` (ライン長 64 の部分被覆で LAYOUT2 と LU のビット一致): **PASS** (dq の差 0、η 1.6e-16)。
- 1 step の時間: B0 (値 0 + 上限 50 + LAYOUT2) と各案を、B0 で挟んで交互に回す。巡ごとに順を変えて 3 巡。
  - 33 本とも確認が OK (GPU 専有・モード・被覆・`nStepInner`・残差と最終場の有限性)。
  - 雑音 (隣り合う B0 の差の最大) は 0.066 ms。

| 案 | 1 step [ms] (3 巡) | 前後の B0 [ms] | B0 との差 (中央値) | 判定 |
| --- | --- | --- | --- | --- |
| `FORGE_LINE_MAXLEN` 100 | 29.95 / 29.92 / 29.93 | 32.14 / 32.11 / 32.12 | −6.8 % | 総時間へ進む |
| 64 | 25.93 / 25.93 / 25.89 | 32.13 / 32.11 / 32.16 | −19.3 % | 総時間へ進む |
| 48 | 24.23 / 24.22 / 24.20 | 32.12 / 32.15 / 32.14 | −24.7 % | 総時間へ進む |
| `nStepInner` 3 | 27.34 / 27.36 / 27.36 | 32.12 / 32.14 / 32.14 | −14.9 % | 総時間へ進む |
| `nStepInner` 2 | 24.95 / 24.94 / 24.94 | 32.12 / 32.13 / 32.13 | −22.4 % | 総時間へ進む |

- この 1000 step では収束の速さは測っていない (step 数への影響は未知)。5 案すべてが「総時間へ進む」なので、どれを組み合わせて総時間を回すかが問い 2 になる。

## 2. 期待値と出典

- 定常の解は R(Q) = 0 で左辺に依らない。ラインも point も残差が下がりきっていないので、どちらの止まる値も解ではない。
- 水準の値 (0.1 kg/s、0.05 %/2 万 step) は §6.11 で事前に登録した。E2 は §6.15 でユーザ決定により外した。
- θ_r の SU2 との照合の許容差は 3 % (`case/45.isobutane_m6_d155/cold_xcheck.py:286`)。
- ふるいの合格は「組の判定が速い かつ B0 比 −5 % 以下 (M 系は照合 PASS も)」(§6.16)。5 % は計算予算の優先の基準である。

## 3. 再現条件

- 準備: `python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext <run> --steps N --out 5000 --cfl 4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --itj 5 --cap 50 [--inner n] [--lvc 3]`
- 実行: `cold_cfl.py run <run>`。環境変数は `FORGE_LINE_LAYOUT=2`・`FORGE_LINE_MAXLEN=N`・`FORGE_LVC_TERMS=5` (L5 のとき)。
- `FORGE_CUDA_BLOCKSIZE=128`。
- 見張りは `m9_watch.py <run> --phase line --budget 200000` (水準で止める)。

## 4. 速度向上策の候補 (呼び出し側の整理)

| # | 案 | 1 step への効果 | 数値への影響 | 準備 |
| --- | --- | --- | --- | --- |
| #13 | LAYOUT2 を既定にする | −7.6 % (測定済み) | なし (ビット一致)。部分被覆も PASS | 残りの経路 (可変長・周期・拘束の行) の確認 |
| #20 | ライン長の上限 (壁から N 節点、先は point) | −7〜−25 % (ふるい) | あり (ライン外の節点は point の Δτ と Jacobi) | 設定だけ |
| #16 | sweep の回数を減らす | −15〜−22 % (ふるい) | あり (ライン間の結合が弱まる) | 設定だけ |
| #17 | 分解を step をまたいで使い回す | 最大 −4.6 ms (−14 %) 程度 | あり | コードの変更 (定常では `lineKFreeze` が効かない、`main.cpp:2070`) |
| #18 | block カーネルが D・K をラインの並びで書く | −1〜2 ms の見積もり | なし | コードの変更 |
| — | 格子を変える (主流側を粗く・細かく、境界層と主流で比を分ける) | 節点数次第 | 解が変わる | 格子の作り直しと別の検証 |

- #14 (壁際の Δτ の上限を上げる) はやめた。上限なしの L5 と上限 50 の L0 で、step 数がほぼ同じだったため。

## 5. 呼び出し側の仮説 (検証していない)

- **H1**: ラインと point の止まる値の差 (θ_r 約 0.1 %) は、どちらも残差床にいることによる反復誤差である。設計上は無視できる (格子・モデルの差は % の桁)。
- **H2**: 水準までの step 数を決めているのは、壁法線方向の剛性ではない。流れ方向 (軸方向 4719 列) へ過渡が抜ける速さである。上限を外しても step 数が同じ (L5 ≈ L0) ことと合う。
- **H3**: ライン長を AR > 20 の範囲 (≈ 64) に縮めても、step 数はほとんど増えない。1 step は −19 %。
  - 根拠: 軸の近くは半径方向に長いセルなので、半径方向のラインが要らない。
  - 一方で、ライン外の節点は point の Δτ になり、ラインの終わりで Δτ が段差になる。
- **H4**: sweep を 2〜3 回に減らすと、ライン間 (軸方向) の結合が弱まる。H2 が正しければ、step 数が増えて相殺される可能性がある。
- **H5**: 緩和 > 0.7 の不安定は、方向別 Δτ (上限 50) と緩和の積が、point で知られた安定の境目 (cfl × 緩和 ≈ 6〜8) を、局所的に超えることによる。

## 6. 読んでよいファイル

- plan: `/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md` (全体。特に §4・§5.1・§6.6〜§6.16)
- `/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md` の run 一覧のうち run_0183・0223〜0263・0353〜0375 の行 (ファイルは大きいので grep で)
- コード:
  - `/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.cpp:1040-1160` (ラインの構築・`FORGE_LINE_MAXLEN`)
  - `/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2048-2100` (sweep の流れ)
  - `/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu` (block カーネル・方向別 Δτ・LP カーネル。8000 行を超えるので grep で該当箇所だけ)
- 判定器: `/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/{m9_judge.py,m9_watch.py,levelreach.py,scr.sh,scr_judge.py,scr_check.py}`
- 過去の諮問とレビュー: `/home/sano/work/forge-integ-1005/notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md`、`/home/sano/work/forge-integ-1005/notes/reviews/2026-10-10-time_integration-line-implicit-speed-plan-6.md`
- 巨大なファイル (`residual_history.csv`・`forge_run.log`・`*.log` のレビューログ) は読まないこと。AWS には入れない。
