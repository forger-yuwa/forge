# codex 諮問 (diagnose): float-geometry-design

- **brief**: [`notes/reviews/briefs/2026-10-10-float-geometry-design.md`](../../notes/reviews/briefs/2026-10-10-float-geometry-design.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-10
- **commit**: `6fece0c7` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.5 min, rc=0
- **結論**: **案 A を目標にするが、まず幾何係数の符号・ガード・境界での定義を修正し、保存場を固定した粘性・拡散の A/B を一度だけ行って、座標差の修復が実際に誤差を減らすか判別する。**
- **extra**: `notes/investigations/2026-10-10-geometry-precision-inventory.md`, `notes/reviews/2026-10-10-m64-and-float-diagnose.md`, `methods/architecture/overview.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **案 A を目標にするが、まず幾何係数の符号・ガード・境界での定義を修正し、保存場を固定した粘性・拡散の A/B を一度だけ行って、座標差の修復が実際に誤差を減らすか判別する。**

確認したコードは `6fece0c7`。対象 run は作業ツリーに存在せず、速度・第一層誤差・Δq/ULP の数値はブリーフと記録によるものです。今回、それらを独立に実測確認したとは扱いません。

**採否表**

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 単一の `ediff = |S|²/|e·S|`：**却下** | [viscousFlux_d.cu:151](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:151) は絶対値＋`1e-30`、[scalarTransport_d.cu:110](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/scalarTransport_d.cu:110) は**符号付き**＋`1e-6·|e|·ss` の相対ガード。[passiveFct_d.cuh:32](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/passiveFct_d.cuh:32) は符号付き絶対ガード。共通係数への置換は精度変更だけでは済まない。**利用側ごとの既存式を維持した係数**を作り、統一は別変更とする。`ss²` を `S·S` に替えることも同一性確認の対象。 |
| **Major** | 「面積・体積は相対精度で十分」、closure は double の半径で修復：**要再検証** | 面ベクトルは読み込みで丸まり、さらに float の半径を掛けている。[mesh.cpp:254](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.cpp:254)、[variables.cpp:599](/home/sano/work/forge-integ-1005/solver_density_cuda/variables.cpp:599)。現行 closure は**流束に渡す面ベクトルの和**。[variables.cpp:637](/home/sano/work/forge-integ-1005/solver_density_cuda/variables.cpp:637)。生の double 幾何だけで closure を作ると、float 面流束との閉包がずれる可能性がある。**座標だけでなく元の `surfVect/surfArea` も double で保持し、幾何精度と離散閉包を別々に検査**する。 |
| **Major** | 壁関数の y を壁節点ごとに保存：**却下** | [ransWallFunction_d.cu:127](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransWallFunction_d.cu:127) は境界面 `ib` ごとの法線で代表内点を選ぶ。角部では同じ節点でも面により代表点・距離が異なる。**境界面ごとの `(irep, y)` を一緒に事前計算**する。WMLES も同様。`conjugateWall` は BC 内で合算した法線を使うため、そちらの定義を壁関数へ流用しない。[conjugateWall.cpp:87](/home/sano/work/forge-integ-1005/solver_density_cuda/conjugateWall.cpp:87) |
| **Major** | §4.2 の対象一覧は十分：**要再検証** | 周期対応付けは丸めた面重心を使用し、ライン構築も丸めた節点を double に広げている。[mesh.cpp:568](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.cpp:568)、[mesh.cpp:1067](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.cpp:1067)。**周期 partner・LSQ の同値類・ライン接続の生成も対象へ追加**する。距離だけ改善して接続が FP64 と違えば、同じ B0 の比較にならない。 |
| **Minor** | node の再構成は `±0.5e`：**限定して採用** | SLAU の条件は `g_reconEdgeMid == 1 && ip < nNormalPlanes`。[convectiveFlux_slau_d.inc.cuh:148](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:148)。**内部辺には正しいが、全境界面へ一般化しない**。非中点経路の `pc−cc` と境界の `fx` は残す。 |
| **Major** | V2 の「残差相対差 ≤1e-12」：**却下** | 相殺後の残差を分母にすると、零近傍で指標が破綻する。さらに前述のガード差は丸め順序ではない。[plan:117](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:117)。**係数・面流束・残差・更新量を分け、相殺前の量で正規化**する。 |
| **Major** | V4 の R3 閾値をそのまま同等性の合格条件にする：**却下** | 元の R3 は「診断なしに採用しない」ためのフラグで、収束保証ではない。[line-implicit-speed:502](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:502)。**0.05 %・0.1 % は工学量の比較限界として維持**し、判定窓・不確かさ・VERDICT を追加する。 |
| **Major** | V6 の `case/08.bump`：**却下** | 現行手順は入口 BC の問題と、修正後も 2 次で残差が頭打ちになることを理由に、回帰基準には使えないと明記。[verification/README.md:23](/home/sano/work/forge-integ-1005/procedures/verification/README.md:23)。下記の node ケースへ置き換える。 |
| **Major** | H2「場全体が止まる心配は小さい」：**要再検証** | FP64 の 9→10 step の変化は、float 計算で要求された更新と実際に反映された更新の測定ではない。[plan:39](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:39)。**`dq_requested` と `Q_after−Q_before` を別々に記録**する。`qAccumulatorFP64` は別 plan とする判断を採用するが、停滞原因の除外には使わない。 |

案 A は**検証後の最終設計として採用**します。新しい常設 opt-in は不要です。比較には固定した旧コミットを使えます。ただし、今回の修正を理由に既存の opt-in を削除することは含めません。

実装は一括変更せず、① double 入力・幾何係数の生成、② LSQ・周期・壁代表点・ライン接続などの前処理、③粘性・拡散・陰解法、④再構成・リミタ、の順に分け、各段で V2 を通してください。①では利用側を切り替えず比較できるようにします。

**第 1 仮説:** 絶対座標を float に丸めてから差を取る処理が、壁際の粘性・拡散評価誤差の主要な寄与である。**確度: 中**

- 根拠: [viscousFlux_d.cu:113](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:113) は float の絶対座標差から距離・拡散係数を算出する。スカラー拡散も同じ構造。ただし「残差約 10 倍の主因」まで実証されてはいない。
- 反証条件: 下記の固定場比較で、正しい係数生成を確認しても壁際の粘性・拡散残差誤差が半減しなければ、**この直接経路だけで誤差の過半を説明する仮説**を棄却する。

**第 2 仮説:** LSQ の変位誤差、面ベクトルの丸めと closure の不整合が、距離修復後にも残る誤差を支配する。**確度: 中、寄与未確認。** LSQ も丸めた座標に依存する。[calcGradient_d.cu:716](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/calcGradient_d.cu:716)

**第 3 仮説:** 幾何を修復しても、残差演算と commit の丸めが θ_r の遅い変化を妨げる。**確度: 中、case/45 での支配性は未確認。** 軸対称の Qacc 制限には実装上の理由があり、単に解除できない。[main.cpp:3371](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:3371)

**判別 A/B:** **変えるものは「粘性・拡散が使う変位と派生幾何係数の生成精度」だけ。時間発展は 0 step、固定場の演算評価を各 3 回。**

- 同じ保存場から状態・勾配・物性・面ベクトルを一度だけ float にそろえ、両腕で固定する。
- 旧腕は現在の絶対座標差、新腕は元の double 座標差から作った float の係数を使用。符号・ガード・面の選択は変えない。
- 参照は、**同じ固定された float の状態・勾配・物性を double に正確に広げた演算**と、double 座標から作る幾何。元の FP64 状態を直接使って入力丸めを混ぜない。
- 見る量は面ごとの粘性応力・熱流束・k/ω 拡散と、その残差への寄与。第一層、壁から内部 3 層、全域を別集計する。ρ残差はこの試験の主判定に使わない。
- **結果 A:** 熱・運動量・k/ω の各対象で新旧誤差比 ≤0.5、かつ再実行差を十分上回る改善  
  → 第 1 仮説を支持し、「直接の座標差誤差は小さい」を棄却する。
- **結果 B:** 誤差比 >0.5  
  → 対応する量について第 1 仮説の「過半を説明」を棄却する。第 2・第 3 仮説を確定したとは扱わない。
- 旧誤差が既に再実行差程度なら判別不能。これは因果の切り分けであり、到達性能の合格試験ではありません。

**やらない方がよいこと:**

- 全カーネルを変更してから、1000 step の残差だけで成功判定する。
- node 境界半割面の `e≈0` に対し、一律に `ediff` を評価する。現行 block-DPLUR はその粘性対角を明示的に除外している。[timeIntegration_d.cu:946](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:946)
- double の closure が精密だから、float の流束とも整合すると仮定する。
- V4 不合格を即座に commit のせいにする、または Qacc の軸対称拒否だけを外す。

**呼び出し側の前提への異議:**

§6 は次の形で事前登録してください。以下の数値は**提案する判定基準**であり、実測で裏付け済みの許容差ではありません。

1. **V1：係数だけでなく定義・接続も検査する。**  
   非退化面の `|e|` と各用途の拡散係数は FP64 比 ≤1e-6。ベクトル成分の零割を避け、`‖e_new−e_ref‖/‖e_ref‖` で判定する。退化面は相対比較せず、有限値・利用禁止の扱いを確認する。周期 partner・LSQ 同値類・ライン接続は参照と一致を要求する。

   closure は `Σ±S_device − A_closure` を、局所の `A_planar` と相殺前の `Σ|S_device|` の**両方**で正規化して記録する。一様圧試験には `P=pRef` だけでなく `P≠pRef` も含める。前者だけではゲージが誤差を隠す。[axisymmetricSource_d.cu:51](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:51)

2. **V2：FP64 の同一性。**  
   同じコミット系列・コンパイラ・フラグ・入力・実効設定で比較する。各保存量について、セル残差差を「参照の各面流束・源項の絶対値和」で正規化し、L1 と最大値の両方を評価する。閾値は `max(1e-12, 旧版の再実行差の10倍)`。面流束、残差、`dq` を別に保存し、20 step 後の状態は固定した物理スケールで正規化して ≤1e-10 を補助条件とする。分岐・ガードの変更があれば丸め差として許容しない。

3. **V3：「同じ状態」を厳密にする。**  
   全保存量・組成・乱流量を共通の `Q32` にそろえ、FP64 側にも `double(Q32)` を渡す。比較は **commit 前の残差場**。全域と壁際で、各保存量の正規化誤差について  
   `E_new ≤ 1.1 E_old + 再実行差の10倍`  
   を非劣化条件とする。旧誤差が再実行差の10倍を超える熱・運動量・乱流の対象には、少なくとも半減を改善条件として課す。全保存量を一つのノルムへ混ぜない。

4. **V4：0.05 %・0.1 % は維持するが、到達点一枚で判定しない。**  
   両者の同一設定区間で、末尾 2 万 step・2500 step 間隔の少なくとも 9 点を比較する。θ_r は各断面別、Q_w は符号を保持し、参照窓平均を分母とする。比較する両窓のドリフト・振幅を閾値の各 1/5 以下に抑え、平均差に両者の時間変動幅を加えても比較限界内であることを要求する。単調なら漸近値も併記する。

   θ_r・Q_w の正規の抽出系列を `check_quasisteady.py --series-csv` に渡す。組み込み `theta` は平板専用であり、そのままノズルへ使えない。既定の drift=0.05 は **5 %**で、この比較には緩すぎる。[check_quasisteady.py:519](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:519)

   `check_convergence` が `NOT CONVERGED` なら結論は「指定水準への到達比較」に限定する。`PASS`・対象量の `STEADY`・判定区間がそろわない限り「FP64 と同等に収束」と書かない。

5. **V5：単価と到達時間を分ける。**  
   専有 native の交互測定で 0.70 倍以下という目標は採用する。ただし採用判断には V4 までの**総時間**も含める。20 万 step を許すと、単価だけ速くても総時間では負け得る。

6. **V6：node の対象を具体化する。**

   | 対象 | 主な物差し |
   |---|---|
   | `case/05.sod_shock_tube` | 同じ物理時刻の厳密解に対する ρ・u・p の L1、逸脱・正値性 |
   | `case/09.Taylor-Green` | 周期の質量・運動量・エネルギー収支、KEEP のエネルギー履歴、受動種収支 |
   | `case/48.flat_plate_cooled_m4` | SST、壁熱流束・摩擦・θ の分布、全残差 |
   | `case/44.vitiated_air_wt` | 軸対称・TP・凝縮の成分収支、出口量、凝縮量 |
   | `case/16.nozzle_wys` の既存 3D node 構成 | 角部・壁代表点・三方向変位、熱流束と全残差 |

   時系列または準定常窓で、新 float の参照からの誤差を `旧 float の1.05倍＋再実行差の10倍` 以下とする。Sod は FP64 より厳密解を優先する。その他の FP64 は**比較参照であって真値ではない**。未実行の壁関数・WMLES 等は、係数検査だけで実行回帰済みと扱わない。

**不足情報:**

- `run_0383*`・`run_0387_dqulp_B0`・B0 参照 run の実体、実効設定、バイナリ由来、判定区間付きの両 VERDICT。
- B0 の到達時系列。元 plan 自身が到達時の HDF5 は削除済みと記載しており、別時点の場では V2/V3 の起点を代用できても、V4 の到達点確認は代用できません。[line-implicit-speed:504](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:504)
- V6 の基準 run ID・固定する物理時刻／窓・再実行差。対象ケースだけでは完全な事前登録になりません。

ファイル変更・forge 起動は行っていません。**plan 未反映**。呼び出し側で `architecture-float-state-double-geometry.md` §4.2・§4.3・§5・§6 に採否と判定条件を反映してください。
