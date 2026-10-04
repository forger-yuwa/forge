# codex レビュー: condensation-two-phase-default (plan)

- **plan**: [`plans/active/condensation-two-phase-default.md`](../../plans/active/condensation-two-phase-default.md)
- **stage**: `plan`
- **date**: 2026-10-04
- **commit**: `c93a6f93` (feature/species-transport)
- **codex**: effort `high`, 7.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **extra**: `plans/active/condensation-two-phase-transport.md`, `methods/condensation.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

段階既定化と `Solver 1 / NonnegLimit 0` の採用方針は妥当です。  
ただし、既定 ON の適用範囲と検証範囲が一致せず、収支・数値許容にも不足があります。S2 の FCT 案は再設計が必要です。

ファイルは変更していません。参照 run 本体はローカルにないため、CFD の再実行は行わず、コード・保存済み VERDICT・台帳と、書き込みを伴わない代数検算で確認しました。

1. **Major — 「検証済み包絡」が、実際の検証より広い**

   **根拠:** [対象 plan:47](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:47) は周期・軸対称・cell を除外せず、`speciesImplicitCoupling 0` や `speciesFaceReconstruction 0` も通します。一方、G2 は case/16 の平面 node・coupling 1・SFR 2 だけです。[親 plan:136](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:136) が本 plan に移管した「ON の周期・軸対称検証」が脱落しています。G0 の Euler は二相拡散が不活性なので代替になりません。

   また、G2 の `nStepInner 5` は、現行生産レシピの [4](/home/sano/work/forge-species/procedures/recommended-settings.md:36) と異なります。

   **対案:** S1 の対応表を先に固定し、周期・軸対称では二相拡散が実際に作動する NS 試験を追加してください。通常の省略 config と、解決後の値を明示した config の比較も必要です。cell は現行方針どおり実行検証を要求せず、未検証として明示してください。実施できない構成を「検証済み」として既定化しないことが条件です。

2. **Major — G3 は補正監視であり、閉じた収支の検証になっていない**

   **根拠:** [対象 plan:85](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:85) は commit・floor・clamp を判定しますが、再正規化と射影は記録のみです。ON では再正規化係数を液・Q にも掛けるため、OFF と `max|f−1|` が同程度でも、液・Q の収支への影響は同じではありません。[speciesTransport_d.cu:2157](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2157)

   この懸念に対して親 plan が要求した [#4k(2) の閉じた収支](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:131) が、新しいゲートに入りません。保存済みの OFF 基準も `VERDICT: NOT CONVERGED (stalled/plateau)` です。[保存済み判定](/home/sano/work/forge-species/notes/investigations/2026-10-03-limiter-inlet/evidence/run_0520_twophase_off_16k/CONVERGENCE_VERDICT.txt:2)

   **対案:** G3 に、同一制御体積での移流・拡散境界流束、相変化ソース、残差、数値補正の収支を追加してください。親 #4k の「不整合が観測された液流束差の 10% 未満」を引き継ぎ、再正規化・射影も成分別に照合します。既に了承されたプラトー受入を撤回する必要はありませんが、準定常性だけで輸送収支を代替できません。

3. **Major — S2 の「液だけ FCT、蒸気は射影」は、親 plan の却下案に戻っている**

   **根拠:** [対象 plan:59](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:59) は総水分を FCT に通さず、液だけを補正します。既存 FCT の制限係数は成分別で、液の上限 `ρY_w` を参照しません。[passiveFct_d.cuh:297](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/passiveFct_d.cuh:297)

   単純な代数反例として、総水分 `[0.10, 0.20]`、液 `[0.09, 0.11]` に保存的な液の補正 `[+0.02, −0.02]` を加えると、蒸気は `[-0.01, 0.11]` になります。既存の上限クランプは液を削るので、この例では総液量を 5% 失います。[condensationRealizability_d.cuh:125](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:125)  
   これは CFD 実測ではなく、成分別保存と蒸気非負が両立するとは限らないことの反例です。

   さらに既存 FCT は SFR ≥2・SLAU 系などの条件でしか作動しません。[speciesTransport_d.cu:1781](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1781) S2 の包絡には、この前提もありません。

   **対案:** S2 は後継 plan に分離し、当面の起動拒否を維持してください。蒸気・液・他の気相種を連成して制限し、総水分を蒸気＋液から戻す方式について、エネルギー流束と BDF 履歴まで定義してから再レビューすることを推奨します。補正量の記録だけでは非負・保存の保証になりません。

4. **Major — G1 の「8ε₃₂ で一致」と過去の p50 は、そのままでは合否基準にできない**

   **根拠:** [対象 plan:83](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:83) は誤差尺度を定めていません。本番は `ρg/ρ` の差から流束を作るので、小さい勾配では差し引きの丸めが支配します。[twoPhaseDiffusion_d.cuh:114](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:114)

   同じ float32 格納値 `ρ=[1,1.000001]`、`ρg=[0.01,0.01000001]` を使った検算では、流束係数を除いた差は float32 で `9.313226e−10`、float64 で `7.078047e−10`、相対差は約 `0.316` でした。流束自身を分母にした `8ε₃₂` 判定なら、正しい演算でも落ちます。

   また、過去の `0.77` は [解析スクリプト:49](/home/sano/work/forge-species/case/16.nozzle_wys/analyze_liquid_diffusion_error.py:49) の節点回復勾配と固定の潜熱・比熱による統計です。本番の面二点差分・面幾何による流束とは評価量が異なります。[speciesTransport_d.cu:2060](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2060)

   **対案:** 同一面・同一格納入力・同一係数で独立 double 参照を作り、差し引き前の項の大きさを含む絶対誤差尺度を事前定義してください。面積・符号・評価領域も固定します。`0.77` は参考値とし、合否には同じ入力と離散化で作った参照値を使用してください。

5. **Major — G2 は、有意差未確認の量まで ON−OFF 差を許容の分母にしている**

   **根拠:** [対象 plan:84](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:84) は全 7 量に「レシピ差 ≤ ON−OFF 差の 10%」を要求します。しかし [親 plan:130](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:130) は onset 差 `−0.0012 mm` を「有意差未確認」としています。これを使うと許容は `0.00012 mm` です。差がゼロなら許容もゼロになり、反復ノイズとレシピ依存を区別できません。

   **対案:** 同条件反復から量ごとのノイズを先に測り、有意な ON−OFF 差がある量にだけ 10% 条件を適用してください。無影響量には別の絶対許容を登録し、判別能力不足は「レシピ依存」と断定せず判定不能にします。比較する 7 量の抽出領域・基準マスクも数値照合が必要です。

6. **Minor — S0 の分類確定と、検査ツールの判定順を先に定める必要がある**

   **根拠:** [対象 plan:76](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:76) は平衡凝縮の分類を S2 より後に置いていますが、分類は S0 の実装に必要です。コード上は `condEquilibrium 1` が液輸送＋緩和ソース、`2` が EOS による液量決定です。[condensationSourceKernels_d.cuh:101](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:101)

   また、G5 の「dual-time＋ON→FAIL」を生のキー値で実装すると、G0 が通す予定の不活性 Euler まで拒否します。

   **対案:** 分類確定を S0 の前提に移し、省略・明示 0・明示 1と、active・inactive・unsupported の判定表を作ってください。ソルバ、検査ツール、HDF5、manifest が同じ実効判定を使い、G5 を既定変更前の必須ゲートにします。

7. **Minor — G0 の参照 run が取り違えられている**

   **根拠:** [対象 plan:82](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:82) の case/44 `run_0487` は、台帳では `run_0487_reach_rhoylim` という再構成分岐の到達確認です。[case/44 README:777](/home/sano/work/forge-species/case/44.vitiated_air_wt/README.md:777) dual-time checkpoint の `run_0487_passiveB_dt_ckpt100` は [case/16 README:403](/home/sano/work/forge-species/case/16.nozzle_wys/README.md:403) にあります。

   **対案:** 意図した凝縮 ON・dual-time・`viscMethod 0` の入力を特定し、完全な run パス、config、IC を固定してください。番号だけの参照では誤った不活性試験を実施する危険があります。

推奨は、**S0・S1 に絞った段階既定化**です。実装前の修正優先順は、①適用範囲と分類表の確定、②閉じた収支の追加、③G1・G2 の許容と参照入力の修正、④S2 の分離です。実効作用素を記録する方針は維持してください。現行 manifest が OFF/ON を同一キーにすることは検算でも確認でき、ここを直す価値は明確です。以上はレビュー提案であり、**plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
