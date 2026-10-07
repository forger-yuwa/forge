# codex レビュー: tooling-nozzle-wall-single-bspline (result)

- **plan**: [`plans/active/tooling-nozzle-wall-single-bspline.md`](../../plans/active/tooling-nozzle-wall-single-bspline.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-07
- **commit**: `10f8bb36` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m2
- **focus**: W0〜W5 の結果 (§9) が事前登録 (§6.0) を満たすか。特に W4 の統合検査の対象を run_0147 から run_0166 に変えた改訂と、判定器の場所ラベルの修正が登録の範囲内か。生産への反映 (§5.1 #6) はユーザ判断
- **extra**: `design/forge_design/geometry/wall_axismach.py`, `design/forge_design/export/wall_step.py`, `case/45.isobutane_m6_d155/wsb_w4_integration.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
ノット挿入による壁表現と保存 STEP の精度は、独立計算でも許容内でした。  
W4 の対象変更・場所ラベルの修正は妥当です。ただし、STEP 判定器の検査漏れと完了記録を修正してから accepted へ移してください。

1. **Major — STEP 判定器が、異なる立体位置・途中で切れた辺を PASS にする**

   [wall_step.py:199](/home/sano/work/forge-integ-1005/design/forge_design/export/wall_step.py:199) は位置誤差を x・y 成分だけで計算し、`poles_max_abs_z` は記録するだけです。また、[freecad_wall_step_job.py:46](/home/sano/work/forge-integ-1005/design/forge_design/export/freecad_wall_step_job.py:46) が取得する `edge_first`・`edge_last` を判定に使わず、評価するのも辺の基底曲線です。

   保存 STEP を基にした読み戻しデータへの負例試験で、次を再現しました。

   - 全体を z 方向に **1 mm** 移動：`pass=True`。
   - `edge_last` を定義域の中央に変更：`pass=True`。

   **対案:** 位置誤差を3次元で判定し、平面性、辺の始終端パラメータ、実際の端点を検査してください。両負例を回帰試験に追加します。併せて、[CLI の終了判定](/home/sano/work/forge-integ-1005/design/forge_design/export/wall_step.py:328)にも回転面の妥当性を含めてください。今回の保存 STEP 自体に、この異常があるという指摘ではありません。

2. **Minor — W4 の参照成果物を固定するハッシュ記録が不足している**

   [plan §6.0 W4](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:132) は参照結果の取得元とハッシュを要求しています。しかし、[W4 の記録処理](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/wsb_w4_integration.py:69)が保存するハッシュは問題 YAML であり、参照した結果 HDF5・両側の `report.json`・壁ファイルのハッシュはありません。結果へのシンボリックリンクだけでは、後日の更新後に比較対象を特定できません。

   **対案:** `run_0166_ns_n012_N1/` と Euler 参照について、実際に読んだ結果ファイル、入力・壁ファイル、比較した報告 JSON のハッシュ一覧を保存し、W4 の記録から参照してください。

3. **Minor — 現行仕様・完了状況・残作業の記述が揃っていない**

   [plan §4.1](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:48)には、ランプ両端のノット、最小二乗、ランプ開始からの入力範囲が現行仕様の形で残っています。[plans/README.md:28](/home/sano/work/forge-integ-1005/plans/README.md:28)も「係数1747・実装待ち」です。現行実装は `poly` 専用・係数1588の検証済み方式です。

   また、[変更ログの未解決事項](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:182)である `run_0147` との δ_r 再現差は、§5.1 の残作業表にありません。

   **対案:** 現行仕様を本文にまとめ、旧方式は履歴と明示してください。W0 の「既定不変」は上流方式などを固定した表現変更の範囲に限定し、別 plan による `poly` 既定化と区別します。δ_r 再現差は残作業表に残すか、引継ぎ先を明記してください。

W4 の変更については、`ba2d8a84` で対象変更が記録されており、`poly` 専用化に伴う変更理由も成立しています。[初版の結果](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/wsb/W4_integration_v1_pathlabels.json)の差分は `wall_resolution.cmd` と `conditions.case` の2か所です。場所ラベルの正規化は、登録された数値評価量の条件緩和には当たりません。

独立確認では、保存壁の最大差は半径 **1.07e−14 r_t**、r′ **1.08e−13**、r″ **3.09e−11**。保存 STEP の7,895点で位置 **4.98e−10 mm**、角度 **2.58e−11 rad**、曲率 **7.39e−12 /mm**でした。構築・一般性・拒否・IC の読み取り専用テスト部分も FAIL 0 です。W3・W4 の HDF5 本体はローカルにないため、ビット同一性は保存判定記録までの確認です。W4 の PASS は報告処理の回帰検証であり、流れの収束を証明するものではありません。

**推奨は、上記1→2→3を修正し、W5を再判定して「機能実装の完了」として accepted へ移すことです。** 生産採用は §5.1 #6 のユーザ判断として残してください。ファイルは変更しておらず、本レビューは plan 未反映です。

指摘数: Critical 0 / Major 1 / Minor 2
