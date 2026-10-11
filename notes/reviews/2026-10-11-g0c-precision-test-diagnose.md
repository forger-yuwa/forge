# codex 諮問 (diagnose): g0c-precision-test

- **brief**: [`notes/reviews/briefs/2026-10-11-g0c-precision-test.md`](../../notes/reviews/briefs/2026-10-11-g0c-precision-test.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `871788cb` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.2 min, rc=0
- **結論**: **上限 50 対 1 の2回目を完了・集計し、上限 50 の陽性対照と G0c 入力の検査を確認してから、B0 固定の G0c で FP64 ×2・float ×2 の30,000 step 比較を行う。**
- **extra**: `case/45.isobutane_m6_d155/make_g0c.py`, `case/45.isobutane_m6_d155/g0c_prep.sh`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 採否 | 指摘・根拠・対案 |
|---|---|---|
| **Major** | 要再検証 | **旧 G0 の結果だけでは陽性対照が不足する。** §6.51 でも旧結果に依存せず、上限 50 の対照を取り直す方針になっている。[plan:1154](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1154)。今回も `run_0717_c50_d1`・`run_0721_c50_d2`・`run_0719_c50_f1`・`run_0723_c50_f2` を使い、両窓・3 断面で大差が再現したことを前提にする。G0c の直接の比較相手は **G0c の FP64** とする。 |
| **Major** | 要再検証 | **「準備済み」を台本の完了だけでは保証できない。** `g0c_prep.sh` は閉性・品質・補間・Q32 化の失敗後も続行でき、最後に `.done` を作る。[g0c_prep.sh:19](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/g0c_prep.sh:19)、[同:32](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/g0c_prep.sh:32)。各必須処理の失敗で停止し、成功した入力だけを採用する。既存成果物についても各検査の成功記録を確認する。これは今回の準備が失敗したという指摘ではない。 |
| **Major** | 要再検証 | **格子の同一性の照合先が足りない。** `make_g0c.py` が比較しているのは、その場で再生成した G0 と G0c であり、実際の `_tr/q32_init.h5` の G0 ではない。[make_g0c.py:24](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/make_g0c.py:24)。実入力との座標・節点対応を照合し、区間外の9保存量の移植一致、変更区間の最近傍対応距離・壁法線方向の対応を確認する。「Q32 化で変化 0」は、正しく移植できた証明にはならない。 |
| **Major** | 修正して採用 | **起点からの増分差と、窓内の値そのものの差は別。** 共通 Q32 入力でも、point 試験の初期化後には精度間で `roe` 最大相対差 `5.52e−4` が記録されている。[plan:1125](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1125)。各精度内の `res_0` 一致と精度間の差を別々に確認する。Δ の比較だけなら結論は「θ_r の**増分の差**が小さい」に限定する。値自体の小差を述べるには、起点差を戻した窓平均の差も確認する。 |
| **Minor** | 修正して採用 | **D の出典と正規化を固定する。** `(+1.892, +1.647, +1.431)` は §6.27 の窓別実測差ではなく、§6.25 から持ち込んだ尺度である。[plan:721](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:721)。G0c の起点を分母にするなら、旧尺度も同じ分母へ換算する。N は2反復で観測した幅であり、信頼区間とは呼ばない。 |
| **Minor** | 修正して採用 | **符号付き条件は「同方向の大差」の検出条件。** 逆符号の大差は現案では判別不能になる。§6.51 は絶対値でも大差を判定している。[plan:1168](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1168)。今回も絶対値で大差を判定し、符号は別記する。 |

結論: **上限 50 対 1 の2回目を完了・集計し、上限 50 の陽性対照と G0c 入力の検査を確認してから、B0 固定の G0c で FP64 ×2・float ×2 の30,000 step 比較を行う。**

第 1 仮説: 縮流部の格子と B0 の反復特性の相互作用が精度差の増幅に関わり、G0c では登録した窓の θ_r の増分差が小さくなる。確度: **中**。  
根拠: [plan §5.1 #22](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:184) にある Grefc の改善は動機になる。ただし別格子の報告であり、G0c の証明ではない。確認できた point 試験の判定 JSON では、W4 の `max|ΔF−ΔD|+N` は `0.0319〜0.0336` ポイントだったが、4 本とも次の判定だった。

- `check_convergence`: **OVERALL: CHECK FAILURES ABOVE**
- `check_quasisteady`: **DRIFTING**（100,000〜120,000 step）

根拠 run は `case/45.isobutane_m6_d155/run_0713_pt_d1`・`run_0714_pt_d2`・`run_0715_pt_f1`・`run_0716_pt_f2`。確認した記録は [rpt_judge_s2.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rpt_judge_s2.json)。「point で精度差を除外した」とは扱わない。  
反証条件: 陽性対照が成立し、G0c でも両窓・3 断面で、反復幅を差し引いた精度間の増分差が旧尺度の半分以上残ること。

第 2 仮説: 格子変更は差の発達速度を変えるだけで、30,000 step より後に大差が育つ。未確認。  
第 3 仮説: 格子移植・初期化と精度の相互作用が応答を変える。未確認。G0c の初期状態の照合が必要。

判別 A/B: **G0c・B0・共通入力・設定・環境変数を固定し、ビルド精度だけを変える。** 各2反復、30,000 step、W1＝20,000〜25,000、W2＝25,000〜30,000。主指標は θ_r(40/70/94)、Q_w は副指標とする。

式は次のように事前登録する。

- S_c＝指定した G0c FP64 腕の step 0 の |θ_r|。全腕に共通で使う。
- Δ＝100 ×（窓平均 − 各腕の step 0）/ S_c。
- D_c＝D_old × S_old / S_c。S_old は旧判定で使った分母。
- eᵢⱼ＝ΔFᵢ − ΔFP64ⱼ、N＝max（精度ごとの2反復のΔの差）。

**結果 A:** 両窓・3 断面で `minᵢⱼ |eᵢⱼ| − N ≥ 0.5|D_c|`。  
→ 第1仮説の有限期間の予測を退ける。符号も報告する。

**結果 B:** 両窓・3 断面で `maxᵢⱼ |eᵢⱼ| + N ≤ 0.1|D_c|`。  
→ 第1仮説の有限期間の予測を支持する。結論は「登録した窓で θ_r の増分差が小さい」まで。窓平均値自体の差も同じ小差条件を満たした場合に限り、「θ_r の値の差が小さい」と述べる。

中間・断面間の混在・前提不成立は判別不能。収束・準定常の VERDICT を併記し、未収束なら有限期間の応答として報告する。

やらない方がよいこと: **同じ AWS B で C1 の4本に G0c の4本を追加しない。** 既存試験の同時実行条件を変え、利用予定の陽性対照にも影響する。準備・判定器の確認は並行でよいが、投入は待つ。待つ理由は上限1の結論が理論上必須だからではなく、対照の確認と実行条件の維持である。また、G0c の小差を「遅いモードを同定した」「double は不要」と解釈しない。

呼び出し側の前提への異議: Gref の「遅いモード」は G0 の精度差の原因としては仮説である。格子変更は離散作用素・局所 dt・初期場の移植も変えるため、今回だけで寄与を分離できない。診断限定の AR ≤ 10,000 はユーザ決定として受け入れる。

不足情報: ローカルには G0c の HDF5、対象 run の生データ、`rcap_judge.json` がなく、準備の実測値と陽性対照を独立検証できなかった。また `rcap_an.py` の設定比較は各精度の1回目だけなので、対照を流用する前に2回目も照合する必要がある。[rcap_an.py:49](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rcap_an.py:49)。G0c に対応した断面抽出の確認も必要。

ファイル変更・forge 起動は行っていない。**plan 未反映**。上記は呼び出し側が §6.52 等へ事前登録するためのレビュー結果である。
