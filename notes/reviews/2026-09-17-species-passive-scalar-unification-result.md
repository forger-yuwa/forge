# codex レビュー: species-passive-scalar-unification (result)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-17
- **commit**: `d67bd0d0` (feature/sern-design)
- **codex**: effort `high`, 7.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: 3 回目 (1・2 回目 NO-GO と plan 3–12 回目はすべて全採用済)。主眼は §6-6/§6-2 の結果: (a) モーメント時間次数 1.3–1.6 の真因を「float32 の丸め床 (核生成率の T 感度で増幅) + 核生成 onset 前線の非平滑性」と結論し倍精度 + 核生成抑止の滑らかブロブで BDF2 (1.92–1.97, nSub 倍増差 0.000) を実証; (b) S3 系列の sub-iter 依存の真因は implicitRelax 単独; (c) FCT 流束形履歴の restart 修正 (case/09 run_0164); (d) ゲート 3 種の VERDICT (PASS/FAIL/PARTIAL/INDETERMINATE) に合否を委ねてよいか。run の実測は case/44.vitiated_air_wt/README.md (run_0330–0381) と case/09.Taylor-Green/README.md (run_0156–0164) を自分で読むこと。問い: 因果の主張に穴はないか、done にして accepted へ移してよいか。設計 §4 の再設計提案は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

モーメントの観測次数と restart 改善は再現しました。しかし、正式系列は収支・全保存量次数の両ゲートで未合格です。
さらに「S3」の検証が実効一次風上になっており、因果の結論と完了宣言を支持できません。ファイルは変更していません。

1. **Major — 実効一次風上の試験で、S3 の検証完了と原因確定を宣言している**

   `run_0372`–`0381` は `convMethod: 0` です。受動種も同じ `conv_scheme` を `interp_dispatch` に渡し、`0` ならセル値を返します。`speciesFaceReconstruction: 2` と FCT の起動ログがあっても、二次面再構成は働いていません。根拠: [run_0379 の設定:30](/home/sano/work/forge/case/44.vitiated_air_wt/run_0379_passiveG_fctorder_double_bdf2_dt8e-6_nsub40/solverConfig.yaml:30)、[受動種の再構成:509](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:509)、[一次風上への分岐:163](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_common_d.cuh:163)。

   再集計した `ρ` の nSub 40/80 差は、無緩和 Venkat `9.044e-6`、無緩和 Barth `8.979e-6`、緩和 0.7 `5.651e-4`。数値は記載どおりですが、この試験ではリミッタが再構成に使われません。元の `run_0322` は `convMethod: 1` なので、「二次 MUSCL 系列の原因は `implicitRelax` 単独」「limiter は無関係」への一般化は未成立です。

   **対案:** `convMethod: 1` を固定した同一 IC・同一精度の緩和 A/B と、実際に二次面再構成が働く次数系列を取得する。現在の結果は「一次風上＋FCT の時間次数」「一次風上での緩和依存」と記録する。

2. **Major — 正式な次数系列が、評価の前提である補正収支も満たしていない**

   `case/44.vitiated_air_wt/run_0378_passiveG_fctorder_double_bdf2_dt1.6e-5_nsub40/` に対する `check_passive_budget.py --double` は **`VERDICT: FAIL`**。クランプ絶対補正の総量比は `Q1=6.345134e-5`、`Q2=8.842479e-5` で、要求 `1e-6` の約63倍・88倍です。`|H_rem|` もそれぞれ `7.94e-5`、`1.21e-4`。根拠: [forge_run.log:206](/home/sano/work/forge/case/44.vitiated_air_wt/run_0378_passiveG_fctorder_double_bdf2_dt1.6e-5_nsub40/forge_run.log:206)。

   `run_0378`–`0381` の次数ゲートも **`VERDICT: FAIL`** を再現しました。モーメントは `1.915 / 1.923 / 1.942 / 1.965`、全列の sub-iter 低下は `≥3.20` 桁ですが、`roUx=1.623`、`roe=1.653` が未達です。誤差上位1%の事後除外は、[§6-6 の全保存量ゲート](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:332)の代替になりません。`run_0379`–`0381` の収支は `PASS` でも、最粗水準の失格は残ります。

   **対案:** 補正収支・場・全保存量次数を同一系列ですべて通す。§5.1 #20/#21/#8 の「達成」「残りはレビューのみ」を未完了へ戻す。前線への誤差集中は再現できましたが、「非平滑性が唯一の原因」という断定も診断上の仮説へ限定する。

3. **Major — `passiveScalarScheme: 0` の旧経路不変契約が破れている**

   新しい実現可能性射影の有効条件に `passiveScalarScheme` がありません。既定 `condRealizProject=1` により、旧経路でも定常・RK の primitive 処理と dual-time 末尾で `Q1/Q2` を変更します。根拠: [condensationTransport_d.cu:172](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:172)、[同:204](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:204)。

   例えば `Q0>0, g>0, Q1=Q2=0` は新射影で単分散へ再初期化されます。[実装:21](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:21)が明示する状態変更であり、丸め差ではありません。

   **対案:** 新射影を scheme 1 に限定し、射影が実際に作動する入力で scheme 0 の旧実装との不変性を検証する。

4. **Major — 収支ゲートが、固定閾値の独立照合まで緩和して偽合格する**

   plan §6-2 は独立照合を `1e-6` に固定していますが、[check_passive_budget.py:157](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:157) は step 数で拡大した `tol_float` を使います。

   メモリ上の実書式ログで、100 step、開始総量 `1`、終了総量 `1.000002`、増分・境界流束・ソース・残りがすべて `0` の反例を通すと、`evaluate()` は **`True`**。独立照合誤差約 `2e-6` を受理しました。また [同:242](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:242) は `--double` にも float32 の step 比例許容を適用します。

   **対案:** 独立照合は `tol` に戻す。丸め許容は対象項と演算精度を限定し、この反例を回帰試験へ追加する。既存 `test_check_passive_budget.py` は `ALL PASS` でしたが、この経路を捕捉していません。

5. **Major — 場ゲートの絶対閾値が、実現不能なモーメントを丸ごと検査対象外にする**

   [passive_gate_common.py:256](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:256) は `rog≤1e-30` を実現可能性検査から除外します。solver の射影条件は `g>0` であり、同じ条件ではありません。

   メモリ上の float32 HDF5 に `ρ=1, T=300, rog=1e-31, roQ0=1, roQ1=roQ2=1e8` を与えると、`check_field()` は **`(True, [])`**。巨大な不整合モーメントまで無視されました。`1e-31` は float32 の非正規化数でもありません。

   **対案:** 正の状態はスケーリングして検査する。塵を除外するなら全モーメントの寄与を評価・記録し、判定できない状態を `PASS` にしない。この反例も追加する。

6. **Minor — 現在仕様と残作業表に、撤回済み説明が残っている**

   [methods/condensation.md:841](/home/sano/work/forge/methods/condensation.md:841) は共通 θ 停止を次数低下の原因と断定したままです。[同:848](/home/sano/work/forge/methods/condensation.md:848) の退化閾値 `x≤1e-3 / y≤1e-6` も、実装の `≤1e-30` と異なります。§5.1 #8/#21 は未達・達成が併記され、今回確認した収支失格が残作業にありません。

   **対案:** 現在の仕様・確定した事実・未検証の仮説を分け、§5.1 を現在の未完了項目に同期する。

restart 修正自体は支持します。[`case/09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed/`](/home/sano/work/forge/case/09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed) は収支 **`VERDICT: PASS`**、連続実行との差 `max|ΔroXi|=3.159046e-6` を再現しました。ただし「流れはビット一致」は厳密には誤りで、`roUy/roUz` に最大 `7.29e-11 / 1.62e-10` の差があります。丸め床の説明も、倍精度化によるモーメント nSub 比の約 `0.4` から `0.006–0.007` への低下が支持します。

**推奨は、`in_progress` のまま残すことです。** 優先順は、偽合格するゲートと旧経路回帰の修正 → 実効二次再構成での因果・次数・収支の再検証 → 文書同期です。修正後も `FAIL / PARTIAL / INDETERMINATE` は完了扱いにせず、必要なゲートすべての `PASS` を移行条件にしてください。本レビューの指摘は、read-only 指示に従い **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
