# codex 諮問 (diagnose): d7275-outlet-corner-plateau

- **brief**: [`notes/reviews/briefs/2026-09-27-d7275-outlet-corner-plateau.md`](../../notes/reviews/briefs/2026-09-27-d7275-outlet-corner-plateau.md)
- **date**: 2026-09-27
- **commit**: `e1e1a26a` (feature/gap-heating-precision)
- **codex**: effort `high`, 4.3 min, rc=0
- **結論**: **(a) を選び、G15 を未達のまま維持して、同一初期場から `time.deltaT.implicitRelax: 1.0 / 0.5` だけを変える A/B を行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| (a) 原因を追う | **採用** | 次は境界条件を変えず、反復の緩和だけを変える。空間離散化の問題と反復不安定を先に分ける。 |
| (b) 比較域だけで収束判定 | **却下・Major** | G15 は全残差の同一設定区間 PASS を要求する（[plan:1030](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1030)）。局在率と熱流束の変更感度は、その代替にならない。**ゲート未達を維持**する。 |
| (c) x=2.6 m で切断 | **現段階では却下・Major** | 出口∩上端の角は残り、下端には出口∩等温壁の角ができる。比較域下流端 2.46 m との距離も 0.34→0.14 m に縮む（[mesh:3](/home/sano/work/forge/case/60.flatplate_d7275_m7/mesh/fp_d7275_y3.geo:3)、[README:35](/home/sano/work/forge/case/60.flatplate_d7275_m7/README.md:35)）。「角を遠ざける／なくす」試験にならない。対案は同一メッシュの下記 A/B。 |

結論: **(a) を選び、G15 を未達のまま維持して、同一初期場から `time.deltaT.implicitRelax: 1.0 / 0.5` だけを変える A/B を行う。**

第 1 仮説: **出口近傍のプラトーは、無緩和の流れ・SST 連成更新が維持する数値振動であり、定常解の不存在ではない。** 確度: **低**
  
  根拠: `case/60.flatplate_d7275_m7/run_0004_t26_wdA` の台帳値は、末尾中央値 `rms_ro=1.01e−6`、`rms_roK=5.44e−2`、`rms_roOmega=0.628`、`NOT CONVERGED (stalled/plateau)`（[acceptance.json:26](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:26)）。生成器は緩和を 1.0 に固定している（[gen_runs.py:71](/home/sano/work/forge/case/60.flatplate_d7275_m7/gen_runs.py:71)）。実装では `implicitRelax` が流れの block 補正と SST の更新量の両方に効く（[timeIntegration_d.cu:1108](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1108)、[update_d.cu:389](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:389)）。**振動の符号付き時系列は未確認**なので、これは原因の認定ではなく、安く判別できる第一候補。
  
  反証条件: 緩和が実際に効いたことを確認し、両腕の過渡が落ち着いた後も、未達残差列と出口近傍の状態振幅がほぼ同水準なら、**「0.5 への緩和で解消する反復振動」説を棄却**する。

第 2 仮説: **node の出口外挿と slip 閉包に共通する空間離散化上の問題。** 確度: 低。境界共有ノードの矛盾・定常解不存在まで示す証拠はない。

第 3 仮説: **AWS の実行精度・コード世代・入力幾何に由来する残留誤差。** 未確認。ローカル入力の座標・面積ベクトル・体積は `float32`。ただし独立に計算した双対面閉性誤差 ‖ΣS‖/Σ‖S‖ の最大は全域 2.23e−10、下流域 4.43e−11 で、**大きな幾何閉性破れを積極的に支持しない**。z 座標も 1 種類で、ローカルメッシュは押し出し 2 ノード問題に該当しない。

判別 A/B:

- **唯一の変更**: A=`implicitRelax: 1.0`、B=`0.5`。同一の検証済みバイナリ・メッシュ・保存量初期場を使い、CFL、内反復数、出口種別、壁距離、空間スキームを固定する。新規 run を作り、`restart_field.py --keep-src-dtype` で初期保存量の一致を確認する。
- **長さ**: 各 20,000 step。全残差履歴を残し、空間残差は 2,000 step ごと。出口上部・下部の代表節点では、開始・終了付近の連続反復について符号付き残差と `ρ,P,k,ω` の振幅も取る。`res_roK/res_roOmega` を診断対象に追加する。
- **事前判別**: 末尾 4,000 step と直前 4,000 step の中央値差が各未達列で 5% 以内になったことを確認する。
  - **B の未達残差列がすべて A の 1/10 以下となり、局所状態振幅も減衰** → 第1仮説を支持。「境界条件が矛盾して定常解を持たない」という説明は退ける。
  - **全未達列と局所状態振幅が A の 0.8–1.2 倍に留まる** → 上記の第1仮説を棄却。共通の空間閉包を候補として残すが、確定しない。
  - 中間的な改善、列ごとの食い違い、過渡継続は判別保留。**1/10 に届かないだけで「効かない」としない。**
- 両腕で `check_convergence` と固定10点の St 系列に対する `check_quasisteady` の VERDICT を残す。この判別閾値は **G15 の PASS を代替しない**。

やらない方がよいこと: **Major — 角を除外して残差ゲートを通すこと、直ちに領域を切ること、`mesh.bndFirstOrder` を使うこと。** 比較域への変更感度が小さいことを、比較量の誤差上限へ読み替えない。対案は上記の、作用素を固定した反復緩和 A/B。

呼び出し側の前提への異議:

- **Major — 出口種別 A/B は、出口閉包を除外していない。** `outlet_statPress` の超音速流出分岐は内部状態を外挿し、`outflow` も内部状態をコピーする（[boundaryCond_d.cu:610](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:610)、[同:1593](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:1593)）。したがって、該当分岐では変更が実質的に効かない。対案は出口各面の局所 Mach・法線速度・分岐を確認し、「出口の種類を替えても改善しなかった」に結論を限定する。
- **Major — 壁距離 A/B の読みは事前登録より強すぎる。** 登録は「残差水準が変わらないなら棄却」だが、実測報告は 27–47% の低下（[acceptance.json:18](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:18)、[同:27](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:27)）。対案は **「壁距離変更だけでは解消しないが、寄与はある」**。同じ A/B の再実施は不要。
- **Major — 「比較域 0.00%」は収束の証明ではない。** 出力した4成分には SST 残差が含まれず、表示丸めもある（[README:35](/home/sano/work/forge/case/60.flatplate_d7275_m7/README.md:35)）。また RMS は保存残差の二乗和で、局所体積による正規化ではない（[residualMonitor_d.cu:132](/home/sano/work/forge/solver_density_cuda/cuda_forge/residualMonitor_d.cu:132)）。対案は成分別の局在率に加え、比較域の絶対 RMS・最大値・体積当たり残差を時系列で示す。
- **Major — FP64 という申告だけでは丸め由来を除外できない。** ローカル入力幾何は `float32` で、その値を読み込む実装である（[mesh.cpp:257](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:257)）。対案は AWS の実バイナリ・型定義・入力ファイルを確認する。今回のローカル閉性検査では、この候補の優先度は低い。

不足情報: AWS の実効 YAML、バイナリ識別情報、全残差履歴、判定区間、両判定ツールの VERDICT 原本、残差局在の集計コード、出口の局所状態時系列。**ローカルに対象 `run_*` はなく、run 数値と VERDICT は台帳の引用であり独立再判定ではない。** 独立確認したのはコード・登録条件・入力メッシュの統計。ファイル変更・`forge` 起動なし。**plan 未反映。**
