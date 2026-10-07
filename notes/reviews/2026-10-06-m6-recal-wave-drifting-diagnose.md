# codex 諮問 (diagnose): m6-recal-wave-drifting

- **brief**: [`notes/reviews/briefs/2026-10-06-m6-recal-wave-drifting.md`](../../notes/reviews/briefs/2026-10-06-m6-recal-wave-drifting.md)
- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **date**: 2026-10-06
- **commit**: `60c22a29` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.1 min, rc=0
- **結論**: **③未達・④保留を維持し、既存スナップショットだけで「最大値探索の空間標本間隔」を変える後処理 A/B を一度行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 0）

| 重大度 | 論点 | 採否・根拠・対案 |
|---|---|---|
| Major | Q1：閾値内なので③を合格とする | **却下。** 登録は振幅と準定常性の両条件。[plan:92](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:92)。提示記録の `wave01: DRIFTING` は未達として保持し、抽出感度を調べる。 |
| Major | Q2：③合格として④へ進む | **却下。** ④は③の全ゲート判定後という登録がある。現状は「振幅は閾値内、準定常性未達、④保留」。`NOT CONVERGED` も残す。 |
| Major | Q3：位置固定・上位 N 点平均で最大値を置き換える | **代替としては却下、補助診断として採用。** 現指標は試験窓全体の最大絶対残差。[nozzle_report.py:211](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:211)。位置固定は移動先のピークを、平均化は局所超過を見逃す。最大値を維持し、最大位置・波形の時間変化を追加する。 |

結論: **③未達・④保留を維持し、既存スナップショットだけで「最大値探索の空間標本間隔」を変える後処理 A/B を一度行う。**

第 1 仮説: 抽出された波形に時間変化が残っており、最大位置の切り替わりだけでは `DRIFTING` を説明できない。確度: **中**。
  根拠: `case/45.isobutane_m6_d155/run_0117_ns_recal_final_ext/` の提示記録は、波の末尾平均 0.00651 %、幅 0.00053 %pt、`DRIFTING`。`run_0116_ns_recal_final/` も `DRIFTING`。両 run の残差判定は `NOT CONVERGED`。ただし、これらは実データ再検算ではなく提示記録による。
  
  [check_quasisteady.py:280](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:280) は末尾平均で規格化した線形トレンドを判定する。報告された 5.4 %/tail は絶対変化約 **0.000352 %pt/末尾区間**に相当し、ゼロ近傍の分母だけで説明する状況ではない。
  
  反証条件: 同じ場・同じ平滑化曲線について最大値探索だけを細かくすると、このトレンドの大半が消え、同じ判定条件で `STEADY` になること。

第 2 仮説: 移動するピークの空間的な取り逃しが、最大値の時系列を変動させている。確度: **低・未確認**。現行抽出はノズル全長を 2401 点で標本化する。[exitM_sampling_ab.py:27](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/exitM_sampling_ab.py:27)。この標本間隔への感度は提示されていない。

判別 A/B: **CFD 0 step。変更は最大値探索の標本間隔だけ。**

- 対象は `run_0116`・`run_0117` の保存済み各12枚。Aは現行標本点、Bは同じ試験窓の各標本間隔を4分割する。
- η=0.1、Md、補間方法、試験窓を固定する。**P-spline は各時刻の現行標本から当てはめた同じ曲線を両腕で使い、Bの追加点で再フィットしない。**
- 各時刻の最大値・最大位置・B−Aを記録し、各 run の末尾5枚（40000〜60000 step）について符号付きトレンドと `check_quasisteady` の VERDICT を比較する。
- **A/B差が全枚で ≤ 5×10⁻⁵ %pt**なら、観測幅の約1割以下なので、探索の粗さを主因とする説を退ける。
- **Bが `STEADY` となり、末尾の絶対トレンドがAから80%以上減る**なら、第1仮説の「最大振幅自体の持続的ドリフト」を退け、抽出感度を支持する。
- 中間結果は判別不能として保留する。Bは診断値であり、その場で登録ゲートを置き換えない。

やらない方がよいこと: 無条件の再延長、`--drift` の事後緩和、合格する固定位置や N の選択、④の自動投入。将来の登録には、最大値を残したまま、要求精度に基づく**絶対変動許容幅と複数窓の確認**を加えるのが妥当だが、今回の合格への読み替えには使わない。

呼び出し側の前提への異議: **「最大位置が跳ぶから値が段で変わる」は原因説明として不十分。** 同じ標本集合の残差ベクトル r に対して、最大絶対値の変化は ‖r₂−r₁‖∞ 以下である。位置の切り替わりだけでは値の変動は生まれず、残差波形・平滑化・標本化のどこが変わったかを示す必要がある。また `extremum at tail-end` は**波指標の時系列の極値**を指し、空間的な最大位置の移動を検出してはいない。[check_quasisteady.py:290](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:290)

不足情報: 対象 run はローカルに存在せず、`quantities_series.csv`、実際の判定コマンド、全残差の VERDICT、最大位置・残差波形を未確認。特に δ_E/δ_C の再検算には、旧参照が既定値として残る `EULER_REF`・`SOLVE_JSON`・`FINAL_PROBLEM` の実効指定も必要。**ファイル変更なし・plan 未反映**。採用時の反映先は当該 plan §5.1・§6.1・§9。
