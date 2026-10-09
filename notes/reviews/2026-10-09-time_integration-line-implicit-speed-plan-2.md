# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-09
- **commit**: `3c2c3c18` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: §6.2 (逆行列の保存 #7 の事前登録) と、その実装 (solver_density_cuda/cuda_forge/timeIntegration_d.cu の lineThomasFactor_d<INV>・lineThomasSolve_d<INV>・lineInvEnabled・比較の経路) だけを点検してほしい。判定の基準・尺度・run の組み方・分岐に穴がないか、実装が §6.2 の記述どおりか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

逆行列の生成・W の LU 代入・前進の行列ベクトル積は、§6.2 と整合しています。ただし、比較・到達判定・計測に誤判定の穴があり、現状の事前登録では採否を決められません。  
予定の `run_0324〜0336` は手元に存在せず、逆行列版の実測合格は未確認です。以下はコード監査と既存 JSON・比較処理の再現結果です。

1. **Major — 比較コードが NaN を「最大差 0」に落とす。**

   **根拠:** [`timeIntegration_d.cu:2434`](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2434) は差のビット列を `atomicMax` で集約し、同ファイル:2451 で `std::max(md, d)` に戻します。`d=NaN` なら比較が偽となり、初期値の `md=0` が残ります。この処理順を host で再現すると、`NaN vs 1 → 最大差 0、不一致 1`、`Inf vs Inf → 最大差 0、不一致 0` でした。有限の差が混在していても、NaN が集約値を占めれば最大差を失います。

   §6.2(1) は dq の不一致件数を失格条件にしていないため、誤合格し得ます。また、[`lu5_factor:1843`](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1843) の失敗判定は微小ピボットだけで、NaN は失敗扱いになりません。`fail` の一致だけでは不足です。

   **対案:** LU/INV 双方の係数・因子・補正について非有限件数を明示的に数え、1 件でも比較失格にする。分解失敗は「不一致件数」に加えて各腕の絶対件数・ライン番号を記録する。NaN/Inf を注入した比較器の確認を、実 run の前に通してください。

2. **Major — 「全ライン・全 sweep」の後退誤差を、予定の dump では検証できない。**

   **根拠:** [`plan:139`](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:139) は factor の第 1・20 回、6 ラインだけを書き出す一方、合格条件を「全ライン・全 sweep」としています。全 4719 本の残りと中間の factor は未評価です。(1) の LU との近さは、両者が共有する誤差を検出できず、後退誤差の代わりにはなりません。

   さらに、[`timeIntegration_d.cu:1514`](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1514) の `FORGE_LINE_DUMP_NODES` は**節点番号**を受け取り、その節点を含むラインを選びます。plan の「ライン 12・33・…」をそのまま渡すと、別の対象になります。

   **対案:** 比較 run の全 solve で、保存済み D/K/rhs と緩和前補正から η を計算し、最大値・ライン・step・sweep・評価件数を残す。6 本の dump は独立した CPU 解との照合に使い、ライン番号と指定節点の対応を登録する。LU も同じ入力で評価し、両方が失格なら「逆行列の精度不足」と断定せず、共通系の問題として判定不能にしてください。

3. **Major — 到達条件が初期出力を誤合格させ、負の大きな収支誤差も通す。**

   **根拠:** [`plan:144`](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:144) は「欠損 2πΣres_ro が初めて 1.0 以下」とだけ定義しています。既存の [`run_0223 の時系列 JSON:5`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/series_run_0223_ns_coldmesh_tw300_linedir_tj5_cap50.json:5) には次の記録があります。

   | step | 保存された `deficit` |
   |---:|---:|
   | 0 | 0.0 |
   | 2500 | 1.243481931 |
   | 5000 | 1.080046174 |
   | 7500 | 0.983260038 |

   登録条件を文字どおり適用した再計算では、到達 step は **0** です。初期出力を除いた隣接点の補間では約 **7068 step** になります。これは到達条件の監査値であり、収束判定ではありません。[`cold_series.py:35`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_series.py:35) は符号付き総和なので、例えば −10 kg/s も現在の条件を通ります。

   **対案:** 未評価の初期残差を除外し、`|2πΣres_ro| ≤ 1.0 kg/s` とする。複数回連続で条件を満たすことを要求し、最初の横断区間と補間値を保存する。1000 step 間隔の補間誤差が ±10% の判定を左右する場合は、閾値付近の出力間隔を細かくしてください。

4. **Major — 「同じ品質」の条件と、未到達時の採否分岐が閉じていない。**

   **根拠:** [`plan:144–149 の開始位置`](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:144) では、質量残差の総和だけで品質を代表させ、`NOT CONVERGED` を許しています。これは途中の到達速度の指標としては使えますが、運動量・エネルギー・SST や θ_r・Q_w の品質を保証しません。§6 本文の [`plan:102`](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:102) にある `check_quasisteady` が、§6.2 の採否条件から抜けています。

   また、以下が未定義です。
   - 片腕だけが閾値に到達する場合。
   - 到達 step 差は ±10% 以内だが、推定総時間が短くない場合。
   - 両腕未到達の場合の「末尾の欠損 × 時間」を、何の基準で採否に変換するか。

   最後の積は、同じ品質までの時間ではありません。

   **対案:** (4) は「欠損閾値への途中到達時間」と明記して予備選別に限定する。片腕・両腕未到達は上限打切りとして記録し、採用候補の判定を保留する。品質の採否には、全残差の `check_convergence` と対象量の `check_quasisteady`、腕間の許容差を事前登録する。必要な判定と監査が終わるまで中間 HDF5 を保持してください。

5. **Major — 計測する `Time` が、登録した「出力を除く時間」ではない。**

   **根拠:** [`plan:103`](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:103) は出力を除く区間を要求しますが、§6.2(3) は最終出力付き run の `Time` を使います。[`main.cpp:3446`](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:3446) で開始した時計は、ループ内の [`writeStepOutputs:2536`](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2536) を経て、:3479 で読まれます。最終 HDF5 出力も計測対象です。

   さらに「中央値が少しでも小さい」を合格にすると、計測のばらつき以下の差でも候補になります。最初の 1000 step の単価を長期の到達 step に掛けた値も、実測総時間ではなく推定です。

   **対案:** ウォームアップ後、出力を含まない同期済み区間を揃えて測り、比較・dump・profiler の環境変数を無効化する。品質 run では到達区間の累積壁時計を直接記録する。反復ペアの差とばらつきを使い、差を識別できない場合は判定不能とする。GPU 専有確認も `pgrep -x forge` だけで済ませず、GPU 上の実プロセスを確認してください。

6. **Minor — `FORGE_LINE_INV` が黙って無効になる併用経路がある。**

   **根拠:** [`timeIntegration_d.cu:2497`](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2497) は `lineInvEnabled()` を呼ぶ前に、ラインなし・`FORGE_LINE_MONO`・`FORGE_LINE_DEBUG_POINT`・`FORGE_LINE_NOOP` で return します。これらが残った環境では、INV を指定しても逆行列経路を通らず、PAR との併用拒否も実行されない場合があります。「起動時に止める」という記述は厳密には成立していません。

   **対案:** 起動時にスイッチの整合を一度検査し、INV とこれらの診断経路の併用を拒否する。各 run に実効モード、factor/solve の実行件数、バイナリ SHA256、精度と環境変数を残してください。

**推奨は、double の逆行列保存を opt-in の限定実験として継続することです。**  
ただし、試験再開前に **①非有限検出 → ②全対象の後退誤差 → ③到達条件 → ④品質・採否分岐 → ⑤時計 → ⑥スイッチ検査** の順で修正してください。カーネルの数式を変える必要は見つかっておらず、先に判定系を直すのが妥当です。

ファイルは変更していません。上記修正案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
