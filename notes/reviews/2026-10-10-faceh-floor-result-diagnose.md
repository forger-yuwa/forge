# codex 諮問 (diagnose): faceh-floor-result

- **brief**: [`notes/reviews/briefs/2026-10-10-faceh-floor-result.md`](../../notes/reviews/briefs/2026-10-10-faceh-floor-result.md)
- **plan**: [`plans/active/time_integration-implicit-thermal-jacobian.md`](../../plans/active/time_integration-implicit-thermal-jacobian.md)
- **date**: 2026-10-10
- **commit**: `c51f9fc4` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 2.6 min, rc=0
- **結論**: **追加runは行わず、上記の限定と証拠の保存先をplanに記録して、§5.1 #5を「探索試験の判定済み」として閉じることを推奨する。**
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/fh_floor_judge.json`, `case/45.isobutane_m6_d155/fh_floor_judge.py`, `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| ① 登録判定 `NOT_SUPPORT` | **採用** | 保存JSONの28評価値から独立に再計算すると、評価器dで低下率 **0.9872%**、u **4.3077%**、合計 **5.2949% < 10%**。評価器fも同判定。[判定値](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/fh_floor_judge.json:2593)、[判定実装](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor_judge.py:219)。ただし、登録した探索規則による判定であり、10%以上の真の効果を統計的に棄却した意味ではない。 |
| H1の「同じ状態のRMS差が小さいから主要因ではない」という根拠 | **Major・却下** | 測定したのは ‖R_d‖/‖R_f‖ で、‖R_d−R_f‖/‖R_f‖ ではない。前者が1に近くても残差ベクトルの差や反復への影響は小さいとは限らない。**H1の根拠は2000 stepの共通評価による登録判定に置き換える。** また、同じ状態の差は一定ではなく、B1の2000 stepでは E_d/E_f＝0.9998538。出発Sの約2e−5だけで代表させない。 |
| 「ドリフト許容内」 | **Minor・文言を限定して採用** | 6本とも保存された傾き・窓差は登録閾値±0.01桁以内。[規則](/home/sano/work/forge-faceh/plans/active/time_integration-implicit-thermal-jacobian.md:299)。記録は「**rms_roeについて、末尾500 stepで登録閾値を超えるドリフトを検出しなかった**」とする。「有意」は統計検定の結果ではなく、全残差・目的量の定常性や残差床への到達も示さない。 |
| ② H2の約1%差 | **観測の記録は採用、因果・改善の確定は要再検証** | 共通評価dでもB平均はAより0.9872%低く、B最大6.47689 < A最小6.48765。したがって、保存された終状態の差は表示器の切替だけでは説明できない。ただし3本ずつの探索試験で、Bの幅は1.4359%。**「今回の標本ではBが低かった。切替による再現可能な改善かは未確定」**と記録する。k・ωはそれぞれ約0.96%・1.65%高く、全残差の改善とも呼べない。[全列の比](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/fh_floor_judge.json:2602)。追加runは今回の判定に不要。 |
| ③ H3のfloat化への含意 | **限定付き採用。FP32への一般化はMajor・却下** | 既定経路はFP64ビルドでもT・Yをfloatへ変換する。一方、診断経路の戻り値は `flow_float` に変換され、FP32ビルドではfloatへ丸められる。今回とFP32では介入の内容が同一ではない。[入力変換](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/thermo_d.cuh:1079)、[切替実装](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:399)。申し送りは下記の限定文とする。 |
| ④ §5.1 #5を閉じる | **条件付き採用** | **「登録した2000 step探索試験の判定済み」**として閉じてよい。§6.2に結果と適用範囲、§6.1に今回の採否、§5.1 #5と§9に完了記録を残す。元CSV・設定・`COLD_PAIR.json`・収束判定・残存HDF5の保存場所も明記する。残差床の原因解明、FP32検証、plan全体の完了とは分ける。[現行の残作業](/home/sano/work/forge-faceh/plans/active/time_integration-implicit-thermal-jacobian.md:111)。 |

float化担当への申し送り文：

> `case/45.isobutane_m6_d155/run_0354_m9_L5cut/res_40000.h5` を出発点に、同一の `lineM_fp64`（sha256 `05ad8bdf…`）・point設定で、面エンタルピー評価精度だけを変えた各2000 step・3本の探索試験を行った。共通評価器dのエネルギー残差はB/A＝0.990128、登録規則では `NOT_SUPPORT`。この条件・期間では、面エンタルピーのdouble評価を必要とする10%以上の改善を支持しなかった。面エンタルピー精度の影響ゼロ、長期の残差床、FP32状態での妥当性、物理精度・速度については未判定である。他の混合精度経路も除外していない。

結論: **追加runは行わず、上記の限定と証拠の保存先をplanに記録して、§5.1 #5を「探索試験の判定済み」として閉じることを推奨する。**

第 1 仮説: 面エンタルピーのfloat評価は、今回の条件で2000 step以内に10%以上改善できる停滞の主要因ではない。確度: **中**。登録判定の算術的再現性は高いが、原因としての一般化はできない。  
根拠: 軌道 `case/45.isobutane_m6_d155/run_0500_fh_a1/`〜`run_0505_fh_b3/`、評価 `case/45.isobutane_m6_d155/run_0510_fhe_s0_f/`〜`run_0537_fhe_b1k2000_dr/` の保存集計で、低下＋u＝5.2949%。保存された `check_convergence --segment` の **VERDICTは6本とも `NOT CONVERGED`**、判定区間は `main`、step 0〜1999。[集計](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/fh_floor_judge.json:2583)、[run索引](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/README.md:189)。  
反証条件: 同条件の有効な独立確認で、登録規則の「低下−u ≥ 10%」かつ「B最大 ≤ 0.9×A最小」を満たすこと。

第 2 仮説: 精度切替が小さな軌道差に寄与した可能性はある。確度: **低**。共通評価での差は観測されているが、再実行のばらつきとの因果的な分離は未完了。

判別 A/B: **今回は実施済みの一組で区切る。** A＝既定、B＝`FORGE_DIAG_FACE_H_DOUBLE=1`のみ変更、各2000 step・3本、主量は共通評価dの終状態 `rms_roe`。登録の支持条件なら短期の寄与を支持、低下＋u < 10%なら今回の主要因説を支持しない。今回は後者である。

やらない方がよいこと: 約1%差を採用基準へ事後変更する、範囲の分離を統計的有意差と呼ぶ、RMS比から残差ベクトル誤差を推定する、今回だけでdouble既定化やFP32の精度配分を決める、自動延長する。

呼び出し側の前提への異議: **観測されたRMS比の小ささと、反復を駆動する残差ベクトル差の小ささは別である。** H1の限定付き結論は採用するが、その補助根拠は採用しない。H2は観測として残し、改善の確定には使わない。

不足情報: この作業ツリーには対象runのCSV・HDF5・実行バイナリ・収束判定原本がない。プロンプト内のplan・判定コード・JSONと実ファイルの一致、JSON内評価値からの判定再計算は確認したが、**315ゲートの実測再確認、ドリフトのCSVからの再計算、実行バイナリの同一性確認はできていない**。`forge`は起動せず、ファイルも変更していない。**plan未反映。**
