# codex レビュー: config-key-pruning (plan)

- **plan**: [`plans/active/config-key-pruning.md`](../../plans/active/config-key-pruning.md)
- **stage**: `plan`
- **date**: 2026-09-18
- **commit**: `2493caf2` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: 3 回目。2 回目 GO-with-changes (C0/M7/m1) を全採用済み: 性能スイッチ 3 件をコードごと復元、棚卸しを入れ子探索 4030 本・コメント除去・節/拒否専用/必須の分離 (値キー 167)・別名解決・絶対許容差撤去へ修正、未知キー検出を完全修飾パス化、wmlesPrt は引数撤去、文書を最終対象 (削除 2 パス + 定数化 5 + 撤去 1、残置 9) に同期、§6.2' の基準を事前確定し §6.4' で 4 経路 + 拒否 8 + 受理 9 + opt-in 3 を全 PASS。問い: (A) 残置決定の読み違い・棚卸しの誤判定はまだ残るか。(B) §6.2'/§6.4' は非退行の根拠として十分か。(C) 第 2 陣に進んでよいか。§4 の基準そのものの再指摘は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

残置した性能スイッチ3件の復元は確認できましたが、棚卸しの誤判断は残っています。  
§6.4' は壁量と一部の分岐を検証できておらず、「全 PASS」を第2陣の着手根拠にはできません。

対象コード `2493caf2` と、レビュー中に更新された文書 `770c9de6` を確認しました。ファイル変更はしていません。

1. **Major — 第2陣で削除予定の `time.deltaT.detectNaN` は、使用ゼロではありません。**

   貼付版後に追加された、実ファイルの [plan §5.3 m](/home/sano/work/forge/plans/active/config-key-pruning.md:170) は、`time.deltaT` 側を「0 run・移送不要」としています。しかし棚卸し再実行では、`time.deltaT.detectNaN` は **記載2,879件、非既定2,781件**。ゼロなのは `time.deltaT.detectNaNInterval` です。

   実例も [case/44 の `run_0117` の config](/home/sano/work/forge/case/44.vitiated_air_wt/run_0117_va_ns_ar5k_iso300_samewall/solverConfig.yaml:12) にあります。[パーサ](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:399) はこの値を実際に読みます。削除すると、多数の既存入力が拒否されるか、無視する実装なら NaN 検知が無効になります。

   **対案:** §5.3 m の削除決定を撤回し、2パスを別々に集計すること。`detectNaN` は互換読みを維持し、将来統合する場合はトップレベル優先の規則を保って入力を移送してください。

2. **Major — 回帰判定ツールは NaN を含む候補を `PASS` にします。**

   [check_field_regress.py:105](/home/sano/work/forge/solver_density_cuda/tools/check_field_regress.py:105) は、非有限値を検査せず `max(cand2, n2)` を取ります。`max(0.0, NaN)` が `0.0` になるため、異常が消えます。

   ファイルを書かない再現試験で、反復3本を `[1, 2]`、候補を `[1, NaN]` とすると、**全7量の差が0、`VERDICT: PASS`、終了コード0**になりました。なお、今回の保存済み16 run の `VALUE/*` と残差196,400行では、非有限値は検出していません。問題は判定器の保証能力です。

   **対案:** 全入力の必須量・形状・有限性を数値比較より先に検査し、不備は即失敗にすること。[既存の `perf_regress.py`](/home/sano/work/forge/solver_density_cuda/tools/perf_regress.py:117) に同種の検査があります。候補側と反復側の双方に NaN/Inf を入れた退行試験を追加してください。

3. **Major — §6.2' で必須とした壁量が、§6.4' の比較から脱落しています。**

   [plan:215](/home/sano/work/forge/plans/active/config-key-pruning.md:215) は `Tau_Wall` / `Qw_Wall` を必須としていますが、[ツール:27](/home/sano/work/forge/solver_density_cuda/tools/check_field_regress.py:27) は任意量に分類し、[81行目](/home/sano/work/forge/solver_density_cuda/tools/check_field_regress.py:81) で全 run に存在する量だけを比較します。読み込むのも `res_<step>.h5` だけです。

   実際に必須量を明示して再実行すると、SSTの `case/26.flat_plate_sst/run_0081_prune_regress_sst_new_r1`〜`run_0084_prune_regress_sst_base` と、line-implicit の `case/39.periodic_hills/run_0028_prune_lineimp_ddes_new_r1`〜`run_0031_prune_lineimp_ddes_base` は、**`FAIL: 必須の比較量が出力に無い: ['Tau_Wall', 'Qw_Wall']`** でした。WMLESの「8量 PASS」にも壁量は含まれていません。

   **対案:** 経路別の必須量を明示し、欠落を失敗にすること。既存の境界ファイルには `twall_*` / `qwall` が保存されているので、要求量との対応を確認して境界出力の比較を追加できます。planと各case READMEの「全 PASS」は、その再判定後に更新してください。

4. **Major — 受理9件の完走は、対象分岐への到達を証明していません。**

   `case/26.flat_plate_sst/run_0085_prune_keytest/logs/solverConfig.base.yaml` は、`thermalMethod: 0`、単成分、`viscMethod: 0`、`unsteady: 0`、`dualTime: 0` です。この設定では少なくとも次が不活性です。

   - `multispeciesRhoYCommonLimiter`: [再構成側の条件](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:233) を満たさない。
   - `passiveFctTolAbs`: [FCTの有効条件](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1692) が要求する dual-time ではない。
   - `turbulentSchmidt`: [化学種拡散](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:813) が `viscMethod != 0` を要求する。

   **対案:** 現記録を「パーサ受理9件」と明記し、この3件には各機能が作動する入力で到達確認を追加すること。値の読み込みログと、消費先の実行確認を分けてください。

5. **Major — 共有コードを変更しているのに、cell の回帰がありません。**

   §6.4' の4経路はすべて `discretization: node` です。しかし SST緩和の変更は [共有wrapper](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:356)、WMLESの引数変更は [cell/nodeで取得処理が分岐するカーネル](/home/sano/work/forge/solver_density_cuda/cuda_forge/wmlesWallModel_d.cu:114) に及びます。[検証手順:40](/home/sano/work/forge/procedures/verification/README.md:40) の両離散化検証を満たしていません。

   **対案:** cell の SST陰解法とWMLESについて、旧新比較を追加してください。node専用変更までcellで再検証する必要はありませんが、共有変更の未検証側は残せません。

6. **Minor — 棚卸しの「既定値のみ」集計に、既定値不明・必須キーが混入しています。**

   [config_key_inventory.py:271](/home/sano/work/forge/solver_density_cuda/tools/config_key_inventory.py:271) は `default_known` を条件にしていません。今回の出力「既定値のみ38件」には、既定を判定できない35件も含まれ、`gpu`、`solver`、`physProp.species` などが混入します。既定値を解決できたものに限定すると3件です。

   **対案:** 非既定数を不明時には `null` とし、「既定のみ」「必須」「既定不明」を分離してください。第2陣の分類表には、`mesh.gradLSQ` のnode時固定など、条件付きの実効既定値も記載してください。

**推奨は、第1陣の変更を維持したまま、検証と分類表を先に完成させ、第2陣の削除実装を保留することです。** 優先順は、①誤った `detectNaN` 削除判断の撤回、②判定器の修正、③壁量・不活性分岐・cell の検証補完、④167パスの分類確定です。レビュー提案は **plan未反映**です。

SSTの再実行判定は旧新とも `check_convergence: NOT CONVERGED`、`check_quasisteady: ALL STEADY` でした。指定時点の回帰として扱う方針は妥当ですが、収束解不変の根拠にはなりません。

指摘数: Critical 0 / Major 5 / Minor 1
