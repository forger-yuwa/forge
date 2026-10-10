# codex 諮問 (diagnose): float-cause-step2-design

- **brief**: [`notes/reviews/briefs/2026-10-11-float-cause-step2-design.md`](../../notes/reviews/briefs/2026-10-11-float-cause-step2-design.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `65094ac4` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.1 min, rc=0
- **結論**: **共通の Q32 初期状態から、全領域・固定30,000 step の float／FP64 軌跡 A/B を、float の再実行1本付きで行う。部分 double 化と領域短縮は先行させない。**
- **extra**: `plans/active/axisymmetric-freestream-hoop-gauge.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 「SST の状態・commit が第一原因」 | **要再検証**。`S_lostfrac` は Σ｜要求更新−実更新｜/Σ｜要求更新｜で、切り上げも含む。内部／境界の集計だけで、符号・収縮部への局在も分からない。[commitLossDiag.cu:61](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/commitLossDiag.cu:61)。**SST は重点観測対象にするが、commit に原因を絞らない。** |
| Major | 「V3 PASS なので残差側は一致」 | **却下**。PASS は旧実装からの相対的改善。記録上、新実装でも ω の残差誤差は参照ノルム比 0.37、壁際の半径運動量は 1.42。[plan §6.10](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:451)。**commit 前の作用素差を候補に残す。** |
| Major | (a) FP64 保存場からの軌跡診断 | **条件付き採用**。元の保存場を両ビルドで読むだけでは、float 側だけ初期状態が丸められる。**保存量を一度だけ Q32 にそろえ、FP64 側にも double(Q32) を渡す。** 既存の [q32.py:10](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/q32.py:10) は `VALUE` の全浮動小数点配列を丸めるので、そのまま使って `wall_dist` まで変更しない。 |
| Major | (b) k・ω の状態／commit の double 化を先行 | **現段階では却下**。状態の精度を上げると、残差への入力も変わり、commit 単独の介入にならない。また既存 `qAccumulatorFP64` は流れの保存量 5 本だけが対象で、k・ω は対象外。[main.cpp:3420](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:3420)。**まず既存の診断で発生域と再現性を絞る。** |
| Minor | 「B で約15分」 | **修正**。提示された単価では、30,000 step の float＋FP64 は合計約26分の GPU 処理量。下記の float 再実行込みでは約37分で、出力費用は別。根拠は [plan §6.17](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:578)。同時投入で専有時の単価が維持されるとは見積もらない。 |

結論: **共通の Q32 初期状態から、全領域・固定30,000 step の float／FP64 軌跡 A/B を、float の再実行1本付きで行う。部分 double 化と領域短縮は先行させない。**

第 1 仮説: **float に依存する作用素・更新の差が、収縮部の遅い結合モードを偏らせている。** 確度: **中**。SST の残差形成を優先して観測するが、発生源を SST に限定しない。

- 根拠: 計画に記録された `case/45.isobutane_m6_d155/run_0460_ab_f32cont/`・`run_0461_ab_f32cont2/`・`run_0462_ab_fp64switch/` では、同じ float 保存場から FP64 に切り替えた腕だけ θ_r が約1.4〜2.2%動いた。ただし総合判定は **「判別不能」**、3本とも **`CHECK FAILURES (NOT CONVERGED)`**。[plan §6.19](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:606)。
- commit 前にも精度差がある。さらに旧格子・固定 Q32 の ω 拡散診断は `r_all = 1.004` で、**その状態・壁際の残差差の過半を拡散面流束の差だけでは説明しなかった**。[plan §6.13](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:515)。これはソース・対流・集積を調べる理由になるが、新格子の長期差全体から拡散を除外する証拠ではない。
- 反証条件: 共通 Q32 から両精度が同方向に大きく動き、精度間の差が下記の小差条件に収まるなら、**この起点・30,000 step では精度差が支配的**という仮説を棄却する。

第 2 仮説: **SST の commit 丸めが偏りを維持する主要経路。** 確度: **低、未確認**。`roK += dk`・`roOmega += dw` は実在するが、残差・対角・更新量の計算も `flow_float` であり、最後の加算だけの責任とは分からない。[update_d.cu:403](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:403)。

第 3 仮説: **共通の長い過渡、restart 時の再構築、再実行の揺れが、停止時点の差を増幅している。** 確度: **低〜中**。新格子の両腕も記録上は **`CHECK FAILURES`／`NOT ALL STEADY`** で、停止 step が120,000と145,000に異なる。[plan §6.25](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:689)。

判別 A/B:

1. **入力と腕**
   - 起点は `case/45.isobutane_m6_d155/run_0483_hp7_k1/res_145000.h5`。同一格子へ `restart_field.py` で9保存量を移し、保存量だけ Q32 に丸めた共通入力を作る。幾何・壁距離は保持する。
   - **A = FP64、B = float、B′ = B の再実行**。B′は別の介入ではなく、再実行差の対照。
   - ソース版、格子、BC、B0 の実効設定、キー1、環境変数、診断条件を固定する。初期保存量と起動後の拘束・再構築による差も記録する。

2. **期間・出力**
   - 各 **30,000 step 固定、500 step ごと**。到達判定で途中停止せず、結果を見て40,000へ延長しない。
   - 主判定窓を **20,000〜25,000／25,000〜30,000 step** の二つに固定する。
   - `output.level: 1` に必要な残差を追加し、全出力を保持する。既存の `FORGE_OMEGA_BUDGET=1` と `omg_trans/prod/dest/cross` も全腕で有効化する。この配列は環境変数なしでは削除される。[variables.cpp:349](/home/sano/work/forge-integ-1005/solver_density_cuda/variables.cpp:349)。
   - `FORGE_DIAG_COMMIT_LOSS=500` は補助記録に使う。既存の集計から局所的・符号付きの commit 偏りは推定しない。

3. **物差し**
   - 主量は θ_r(40/70/94)。Q_w は独立した副判定として残す。旧試験の Q_w の不合格を、今回の主判定に合わせて書き換えない。
   - 場は **ρk・ρω、k・ω、μt、ρ・P・T・速度**を、収縮部全体、その内部 `x/r_t ∈ [−5,−1), j=20〜60`、軸近傍、壁際、下流に分けて追う。
   - 正規化尺度と領域は起点から固定する。量ごとの領域 L2、差の分位点・符号・位置を使い、微小な μt を分母にした「最大30倍」だけで順位を決めない。
   - 局所差の立ち上がりは、同じ物差しで **B−A が B−B′の10倍を超える状態が3出力連続**した最初の区間として記録する。500 step 未満の先後は不明とする。**先に離れた変数は、原因の確定ではない。**

4. **事前の分岐**
   - 既存の停止時差 `D = (+1.892, +1.647, +1.431)%` は、真値ではなく今回の効果量の尺度にだけ使う。
   - **精度依存を支持する結果**: 両窓・3断面で、BとB′がAから既存 float 側へ **0.5D以上**離れ、その差が各窓のB−B′の平均差の10倍を超える。A自身の起点からの移動は **0.1D以下**。→ 共通の過渡だけという説明を退け、第1仮説を支持する。**SST commit の証明にはしない。**
   - **共通過渡を支持する結果**: 両窓・3断面で、精度間の差が **0.1D以下**で、両精度が起点から **0.5D以上同方向に動く**。→ 第1仮説の「この期間で支配的」を棄却し、第3仮説を優先する。
   - 両方とも動かない、再実行差が大きい、窓間で反転する場合は **判別不能**。上記係数は暫定の診断閾値で、信頼区間ではない。
   - 全期間の非有限検査と、区間を明示した `check_convergence`・対象系列の `check_quasisteady` の VERDICT を併記する。未収束でも軌跡の応答は評価できるが、定常解の一致とは呼ばない。

やらない方がよいこと:

- `S_lostfrac = 4.9%`、`rms_roK` 比2.83、μt の最大比を足し合わせて「SST commit が原因」と結論すること。μt は速度勾配と壁距離にも依存する。[turbulent_viscosity_d.cu:205](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu:205)。
- `qAccumulatorFP64` の軸対称拒否だけを外すこと。k・ω は対象外で、軸の基準状態の射影にも未対応である。
- **今回の A/B と同時に x/r_t≈3 で領域を切ること。** 差の局在は下流からの影響がない証拠ではない。短縮には、切断面全体の法線 Mach 数・逆流・壁際の亜音速域、粘性／SST 拡散と出口条件、保持領域の幾何・壁距離・接続の維持を確認する必要がある。さらに保持領域の指標が全領域の場合から十分小さくしか変わらないことが必要。列を61%残すなら、節点数比例の理想化でも速度向上は約1.64倍に留まる。

呼び出し側の前提への異議: **観測は「未収束 run の停止時差」と「精度切り替えへの応答」であり、「異なる定常解」や「SST の commit が主因」ではない。** case/48 は候補を除外しない。FP64 ビルドでの面エンタルピー・LHS の試験も float の候補を除外しない。一方、case/45 の抽出座標は `res_*.h5` ではなく `nozzle.h5` を読むため、case/48 の座標抽出問題をそのまま再提示する根拠もない。[cold_xcheck.py:107](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:107)。

不足情報: 対象 run、判定 JSON、残差 CSV、保存場はこの作業ツリーに存在せず、**今回、実測値と VERDICT の独立再集計はできていない**。上記の数値は計画の記録に基づく。投入前に呼び出し側で原本、バイナリ・入力のハッシュ、実効設定を照合する必要がある。ファイル変更・forge 起動は行っていない。**plan 未反映**であり、反映先は `architecture-float-state-double-geometry.md` §5.1 #15 と新しい検証節とする。
