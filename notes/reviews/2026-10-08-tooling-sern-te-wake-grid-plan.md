# codex レビュー: tooling-sern-te-wake-grid (plan)

- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **stage**: `plan`
- **date**: 2026-10-08
- **commit**: `11c7ee0e` (feature/sern-design)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
局所格子変形を生産候補として検証する方針は妥当です。  
ただし、床計測の導入順序・生産ゲートへの接続・学習系列の識別・2D 検証を、実装前に具体化する必要があります。

目的・スコープについて、`plans/README.md` と関連する `accepted/` を確認しました。本件の生産受入れが解決済みという記録はありません。`w` の受入れと有限厚メッシュ開発から分離する判断も妥当です。既存 Hermite 曲線は端点の位置・勾配を満たし、途中で傾きが旧線より急になる性質も[実装に明記](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:425)されています。復帰区間まで幾何・温度を監視する設計を支持します。

実測の再確認では、`case/46.sern_design/run_1078_tewake_A0_m10/`・`run_1079_tewake_B10_m10/` の[保存 CSV](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-te-wake-ab/)を公式 `check_quasisteady.py` の関数で再判定しました。区間は **10000 < step ≤ 20000**、温度4量・力4量は双方 **`STEADY`**。B の R_TE 最低温度は **227.3638 K**、B−A′ は計画記載の係数差を再現しました。ただし、残差の保存 VERDICT は双方 **`NOT CONVERGED (stalled/plateau)`**。メッシュの保存判定は **`PASS`**、B は **`ADMISSION VERDICT: PASS`** です。入力・結果 HDF5 と残差 CSV はローカルにないため、この部分は保存記録への依存です。run 索引は [case README](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:391)です。

1. **Major — 床カウンタより先に行う g4 計算では、登録した受入れ条件を判定できない。**

   **根拠:** [§5・§5.1](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:54)は g4 → カウンタの順ですが、[§6 #4](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:83)は両格子で毎更新の床補正ゼロを要求しています。既存 B の CSV にあるのは保存時点の `n_floor_plus1` で、過去の更新事象を復元できません。g3 の再実行を[未確定事項](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:109)に残すのも、この条件と矛盾します。

   **対案:** カウンタ実装・検証を先行させ、その共通バイナリで g3 と g4 の判定区間を取得してください。既存 B は初期場として使用できます。更新単位、実節点数、計測区間、restart 時のリセット、欠測・切り詰め・整数オーバーフローの検出を記録形式として固定し、「記録がない」をゼロ件として扱わないことを試験します。

2. **Major — カウンタを出力するだけでは、生産評価への不良データ混入を防げない。**

   **根拠:** 現在の [`floor_gate()`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:200)は最後の HDF5 だけを検査し、[`evaluate_gates()`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:279)もその結果を採用します。計画の実装対象には、この受理経路の変更がありません。判定区間中に床補正され、その後回復した評価は、現行経路では通り得ます。

   **対案:** 新しい生産系列では毎更新カウンタを必須入力にし、`sern_gates` → `metrics.json` → 台帳の採否まで接続してください。少なくとも「途中で1件・最終場正常」「記録欠落」「restart 後の区間欠落」を拒否する試験を追加します。旧評価の持ち越しは、§6 #6・#7 に基づく明示的な例外として記録します。

3. **Major — 実格子の同一性と、サロゲートが学習する評価系列の互換性が区別されていない。**

   **根拠:** [§4](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:42)は実格子署名・双対ハッシュ等を識別に含め、「識別が違えば別系列」としています。しかし [`driver_sern` の設計変数](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:41)には `L_cowl` や角度があり、正常な設計探索でも座標・接続署名は変わります。この署名の一致を学習条件にすると設計点ごとに系列が分裂し、記録するだけなら旧方式の混入を防げません。現在の[学習可否判定](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:91)にも、今回の互換性契約はありません。

   **対案:** 次の二つを分離してください。
   - `stage_key` と個別評価の来歴には、実入力の格子署名・双対ハッシュを保存する。
   - 学習可否には、曲線版、blend 指定、格子生成レシピ、変換器の版、離散化等からなる評価方式の識別を使う。

   「同じ方式の異なる設計点は一緒に学習可能」「同じ形状でも旧方式・不明来歴は除外」を両方試験してください。持ち越し承認の対象と根拠も別項目にします。

4. **Major — 2D の検証は、既存ツールを領域の読み替えだけで実行できない。**

   **根拠:** 格子検査の [`_read_h5()`](/home/sano/work/forge-sern-design/case/46.sern_design/diag/te_wake_grid_check.py:146)はヘキサ接続を要求します。温度監視も[側端を z = 0.1 m に固定](/home/sano/work/forge-sern-design/case/46.sern_design/diag/te_monitor.py:27)しており、H = 0.1 m の平面 z = 0 にそのまま適用すると、半径 0.2H の R_SE は空です。

   また、[§6 #10](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:87)は A/B のどちらに #1〜#3 を要求するか不明です。症状再現を要求する A に低温排除条件も課すと、試験が成立しません。

   **対案:** 2D 用の格子品質・双対閉性・実辺角度・層厚検査と、形状から監視領域を生成する処理を実装項目に追加してください。存在しない側端領域は明示的に対象外とし、必要領域が空なら判定不能にします。判別条件は **A＝症状再現、B＝#1〜#3 合格**と明記し、初期場の作成方法、2D runner のオプション受け渡し、既定0の不変試験まで登録してください。

5. **Minor — 「残差・保存量がビット一致」は、現在の回帰手順と整合しない。**

   **根拠:** [§5.1 #3](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:67)は残差のビット一致を要求しますが、[検証手順](/home/sano/work/forge-sern-design/procedures/verification/README.md:78)は同一バイナリでも `atomicAdd` の加算順により残差が変わるため、決定的状態と時間発展を分けています。「無効化不要」と「無効時」の関係も未定義です。

   **対案:** 共通入力に対する EOS 更新前後の保存量など、決定的比較ではビット一致を要求してください。組立残差・時間発展には変更前だけで校正・凍結した許容幅を使い、カウンタ有効時を検証対象に明記します。周期・軸対称・凝縮への適用範囲も指定してください。`cell` の実行回帰は[現行規則](/home/sano/work/forge-sern-design/procedures/verification/README.md:63)どおり追加せず、未検証範囲として残す扱いで構いません。

**推奨は、`te_wake_blend_H: 1.0`・`w` 無効を維持し、計測と受理経路を先に完成させてから生産受入れ試験へ進むことです。** 実装前の修正優先順は、①カウンタ・生産ゲートの契約、②評価方式と格子来歴の分離、③共通バイナリでの g3/g4 取得順序、④2D の検証仕様、⑤回帰基準です。g4 の0格子は投入条件の基準用と確定してよく、今回の目的には、その流れ計算まで必須にする必要はありません。

ファイル変更・`forge` 起動は行っていません。指摘と推奨は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
