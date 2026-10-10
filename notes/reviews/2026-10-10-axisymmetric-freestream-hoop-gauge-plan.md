# codex レビュー: axisymmetric-freestream-hoop-gauge (plan)

- **plan**: [`plans/active/axisymmetric-freestream-hoop-gauge.md`](../../plans/active/axisymmetric-freestream-hoop-gauge.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `78638540` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **focus**: §4.5 の設計 (区間ごとの r 重みの面ベクトル、ユーザ決定) と §4.6 の合格条件に集中。§4.3・§4.4 の幾何の A/B の結果も踏まえて。実装の前の点検

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
案 (a) の W_f = Σ r_k S_k は、集約後に半径を掛ける誤差を直す適切な方針です。  
ただし、§4.6 は壁で残る欠損・float32・非静止場の退行を十分に検出できず、このままの既定化は承認できません。

対象の `run_0183`・`run_0252`・`run_0401`・`run_0452` は参照可能な作業ツリーに存在せず、実測値と VERDICT は独立に追認できませんでした。以下では plan の数値を**報告値**と明記し、コード検分と書き込みなしの CPU 演算試験を根拠にしています。

1. **Major — §4.4 の不合格原因を未解決のまま、§4.6 で合格基準を置き換えている。**

   **根拠:** [plan:86](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:86) の報告値では、B は全域条件に x 2,664 CV・y 1,534 CV が不合格です。座標丸めという説明は未確定ですが、[同:113](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:113) は全域 1e−6 に切り替えています。この変更だけでは、残存欠損が丸めなのか境界構築の不整合なのかを判別できません。

   コードにも検査すべき違いがあります。内部区間は丸めた中点から `G−M` を作り、境界半割ベクトルは `0.5·(B−A)` で作ります（[gmshReader.hpp:1493](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:1493)、[同:1657](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:1657)）。浮動小数点では、後者と `M−A` は同一とは限りません。

   **対案:** 不合格 CV と隣接区間を抽出し、同じ入力座標から高精度で W・双対面積を再計算してください。端点生成・区間積・集約の誤差を分離し、座標丸めを含む許容誤差を説明してから基準を改訂するべきです。旧 FAIL は維持し、「厳密に釣り合う」は厳密演算上の恒等式に限定してください。

2. **Major — 静止場の「全域残差」試験は、問題が残った壁 CV を検査できない。**

   **根拠:** [main.cpp:2015](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2015) は軸・壁の残差を射影し、[nodeWallDirichlet_d.cu:56](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:56) は壁の `res_roUx`・`res_roUy` を明示的にゼロにします。参照する既存集計器は出力済み `VALUE/res_roUy` を読むだけです（[fg2_an.py:126](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fg2_an.py:126)）。

   したがって、§4.4 で不合格の大半を占める壁について、§4.6 の残差・速度試験が通っても幾何閉性の証明になりません。また、変更する W_x の初期残差基準がありません。

   **対案:** デバイスへ渡した最終 W から、**境界拘束とは独立に** E_x・E_y を全 CV で評価する試験を必須にしてください。流れの試験は射影前後を分け、自由 DOF の `res_roUx`・`res_roUy` を測定します。一定 U の質量残差試験も残してください。軸では hoop source 自体を省略するため、射影前の半径残差を通常 CV と同じゼロ条件にしてはいけません。

3. **Major — float32 を既定対象にするのに、float32 の自由流保持と軸面の合格条件がない。**

   **根拠:** [plan:110–115 相当](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:110) は float のビルド・幾何転送を確認しますが、静止場と変更経路の生産試験は FP64 です。W を double で正しく作っても、各面を float に丸めれば閉性は再び崩れます。

   CPU の合成試験でも確認しました。半径 1、厚さ 1e−7、軸長 1e−4 の長方形 CV（AR≈1000）で、面ごとの W_y を float32 に丸めると、**double で集計しても**閉性欠損 / A は約 **0.455** でした。これは実 run の値ではありませんが、「double 生成で float の閉性も保証される」という前提への反例です。

   軸床も要注意です。W=(0, −1e−24, 0) のノルムは、float の二乗和では **0**、double でノルムを取ってから float にすれば **1e−24** でした。下流は `ss` で除算します（[setDT_d.cu:64](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:64)）。

   **対案:** float32 のキー 0/1・再実行対照を追加し、領域別の閉性欠損、静止場残差、10 step の速度に定量的な非劣化条件を置いてください。軸面判定は「W がほぼゼロ」ではなく幾何条件で定義し、床を適用したベクトルと `ss` の計算順序・有限性・正規化条件を明記してください。

4. **Major — 非静止場で悪化しても合格できるため、「既定 1」の採用ゲートになっていない。**

   **根拠:** [plan:104](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:104) は既定 1ですが、[同:115](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:115) は解への影響を「合否に使わない」としています。

   W は旧 S と一般に平行でなく、変更は圧力項だけに閉じません。粘性の `|S|²/|e·S|`（[viscousFlux_d.cu:169](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:169)）、block-DPLUR の法線と面積（[timeIntegration_d.cu:883](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:883)）も変わります。静止場ではこれらの影響を十分に試せません。

   **対案:** §4.6 #6 に最低限、非有限・非物理値・発散・規定 step 上限への未到達を不合格とする条件を追加してください。同一バイナリ・同一メッシュ・同一保存量からキー 0/1 を比較し、既存 `run_0452` を流用するなら入力と旧経路の同一性を確認します。量の差は修正効果も含むため、旧解との一致を機械的に要求せず、許容差超過時は既定化を保留する設計が適切です。

   なお、参照先の停止規則は「水準への到達」であり収束判定ではありません。§4.1 #4 の「収束まで」と整合させ、`PASS`・対象量の `STEADY` がなければ定常解比較とは呼ばないでください。

5. **Minor — 新データセットが存在する場合の適用除外と、境界集約の検証が不足している。**

   **根拠:** §4.6 の平面・3D 試験では新データセットを生成しないため、「データセットはあるが `axisymMethod=1`／`axisRFloor>0` なので無視する」分岐を検査できません。また、基準スクリプトは quad 専用です（[dual_segment_closure.py:39](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/dual_segment_closure.py:39)）。変換器には入口半割面を壁へ移す所有変更もあります（[gmshReader.hpp:1668](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:1668)）。

   **対案:** 適用条件の組合せを読み込み時の配列比較で検証し、無効分岐の旧経路一致を確認してください。小さな tri・曲がった境界・複数マーカ角の幾何試験も追加します。cell の流れ計算を新たに要求する必要はありません。

6. **Minor — 現在の課題と、過去に解決済みの課題の記述が混在している。**

   **根拠:** §1・§2 は既実装の `pRef`／closure 修正を中心に説明し、[methods/axisymmetric/implementation.md:399](/home/sano/work/forge-integ-1005/methods/axisymmetric/implementation.md:399) も自由流閉性を float32 起因として説明しています。一方、今回直すのは FP64 にも残る集約誤差です。`precision-mixed-axisym.md` は `accepted` ではなく `archived` にあり、float 精度を真因とする説明も撤回済みです（[同文書:20](/home/sano/work/forge-integ-1005/plans/archived/precision-mixed-axisym.md:20)）。

   **対案:** 冒頭に今回の追加スコープを明示し、関連リンクと現在仕様を更新してください。既存の double 幾何化・粘性対角修正とは別問題であり、今回の修正自体は重複ではありません。残作業表は「残存欠損の説明 → データ生成・読込 → 精度別／適用分岐の検証 → 既定化判断」の順に具体化してください。

**推奨は案 (a) の継続です。** 実装前に、優先順で **①全域 FAIL の原因確認、②射影に隠されない x/y 検査、③float32・軸床の仕様と試験、④非静止場の採用ゲート**を plan に追加してください。再変換時のメッシュ品質確認と、同一メッシュへの保存量の厳密な restart コピーも検証手順に含めます。`hoopAreaFromClosure` の既定化は混ぜない方針を維持してください。

ファイル変更・solver 起動は行っていません。レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
