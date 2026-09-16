# codex レビュー: species-passive-scalar-unification (result)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `result` (diff base `3b62903a`)
- **date**: 2026-09-16
- **commit**: `45747ae0` (feature/sern-design)
- **codex**: effort `high`, 7.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **focus**: 2 回目 (result-1 の M1–M6/m1–m2 反映後; §6.1 採否表と §5.1 #9–#19、§9 の最新項)。前回の再指摘は不要。実物で確認: M1 周期 DPLUR (dplurSweepOnce, case/09 run_0090/0091 [coupling 1], run_0109/0110 [TP 周期箱 coupling 2 u=10]), M2 root-only 収支 (run_0088/0089/0098), M3 checkpoint 整合 (case/09 run_0094–0097, case/44 run_0233–0235), M5 φ_N·δρ 項と増分スケーリング (case/16 run_0500/0501, case/34 run_0109 の floor ≤1e-4, 固定点不変), M6 実効緩和 (case/16 run_0497) と check_convergence --from-floor (run_0478/0479/0480/0497 PASS), M4 凝縮 dual-time (case/44 run_0236–0262): sub-iter 床が流れの 2 次 MUSCL 由来で convMethod 0 では流れ・化学種 ≥2.1 桁、しかし**モーメントの BDF2 次数は 1.3 で未達** (実現可能性クランプが常時作動するセルの床; §10 に理由と後続を記載) — この未達を「範囲外の制約整合の時間離散」として accepted に移してよいか、それとも必須残作業か、根拠つきで判定してほしい。既定 passiveScalarScheme 1 / SFR 0 のまま。cell は対象外。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

凝縮 dual-time の次数未達は、本 plan の必須残作業です。「制約セルだけ時間 1 次」とする根拠は不足しています。  
さらに、増分制限後の保存誤差と、S3 リミッタの周期整合に未解決事項があります。

`3b62903a...45747ae0` の diff と指定 run を確認しました。前回修正された周期 DPLUR、root-only 集計、checkpoint の不整合時 BDF1 再開は確認できたため、再指摘しません。`run_0478/0479/0480/0497` はツール再実行で **`PASS (within 1.5x of reference floor …)`**、`run_0500/0501` は **`PASS (converged)`** でした。

1. **Major — 新しい共通増分制限が停止点を作っており、次数未達を「範囲外の物理制約」として処理できない。**

   [condensationUpdateLimiter_d.cuh:88](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:88) は、いずれかのモーメントについて `N_k=0, d_k<0` なら共通 `theta=0` とし、**他のモーメントの増分も停止**させます。直前の `theta≥1e-12` 保証を、新設された処理が無効化しています。

   実測でも、`case/44.vitiated_air_wt/run_0257_passiveD_order_bdf2_dt8e-6_nsub40/res_200.h5` の `condLim=0` は **663 ノード**。全て `Q1=0`、一方で `min(Y1−g)=0.0376953` なので、この集合は蒸気枯渇ではありません。`condClampCorr`／`condClampCorrQ` も全てゼロです。「硬い実現可能性クランプが原因」という説明では、今回追加した更新停止と従来の状態補正が区別されていません。

   最新 `run_0256–0259_passiveD_order_*` の最終 HDF5 から、同じ非重み付き L2 差で再計算すると：

   | 量 | BDF2 観測次数 | sub-iter 倍増差／最小水準差 |
   |---|---:|---:|
   | `g_0` | 1.303 | 0.431 |
   | `Q0_0` | 1.360 | 0.373 |
   | `T` | 1.274 | 0.073 |
   | `ro` | 2.290 | **0.170** |

   温度も次数未達、密度も反復誤差の `≤0.1` ゲート未達です。さらに最小刻み `run_0258_passiveD_order_bdf2_dt4e-6_nsub40/residual_history.csv` の `rms_roY1` は sub-iter 低下 **最小 1.867 桁、中央値 1.937 桁**。「流れ・化学種 ≥2.1 桁」は系列全体には成立しません。

   **対案:** §5.1 #12/#18 を必須のまま残す。制限理由別に、停止セルの未制限残差と誤差への寄与を記録し、境界上で全モーメントを停止させない制約付き更新を検証してください。反復誤差を十分小さくした後に次数を再評価すべきです。現在の数値から「作動セルは1次、非作動セルは2次」とは判定できず、制約セルの残差を除外するだけでは解決になりません。

2. **Major — `floorCorr` の減少をもって保存ゲート合格とはできない。増分制限による質量損失が残っている。**

   [passiveKernels_d.cuh:112](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:112) は CV ごとの増分を独立に縮小します。補正を隣接 CV に逆符号で戻す処理がなく、保存的な面流束制限ではありません。

   `case/09.Taylor-Green/run_0098_passiveC_m2_step_seam_s3_final/` の `res_0.h5` と `res_500.h5` を周期 root 相当の代表点のみ、合併体積で積分すると：

   - 総トレーサ量：`83.71696370 → 83.26320518`
   - 相対損失：**`5.42015e-3`（0.542%）**
   - ログの符号付き `limCorr`：`−4.538e-1`、`floorCorr=0`

   収支の説明は正しくなりましたが、[plan §6-2:197](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:197) の保存誤差 `1e-6` は満たしません。また最新 `run_0500` の累積 `limCorr` 比は `1.728e-4`、Arthur の `case/34.arthur_n2_nozzle/run_0109_passiveD_rhoterm_arthur_s1_sfr2/` では `Q0` が `7.208e-3` です。`floor` だけを判定対象にすると、保存量補正が別の診断欄へ移った分を見落とします。

   **対案:** 周期の非定常保存試験には、両 CV で共有する面流束の制限を採用し、総量保存を確認してください。`floorCorr` と `limCorr` を併記し、定常計算の起動時緩和と、非定常計算の保存誤差を分けて合否判定する必要があります。

3. **Major — 周期 DPLUR は修正されたが、S3 の受動種リミッタは部分 CV のステンシルのまま。**

   [passiveLimiter_d.cuh:35](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveLimiter_d.cuh:35) は各 member の内部面だけから極値・局所スケール・リミッタを計算します。[limiter_d.cu:315](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:315) に周期群の極値統合やリミッタ同期はありません。

   `case/09.Taylor-Green/run_0089_passiveC_m2_rerun0068/res_500.h5` では、同一周期点のノード `4532`／`2982` で：

   - `Xi=0.828983128`、`dXidx=0.371410549` は同値
   - `limiter_Xi=0.985149026`／`1.0`、差 **`0.014850974`**

   状態 mirror が一致しても、再構成は合併 CV と同じ作用素になっていません。今回の `run_0090/0091` は KEEP、`run_0109/0110` は SLAU `convMethod: 0` なので、確認できた π シフト等価性はこの S3 分岐を検証していません。

   **対案:** 周期群で隣接値の最大・最小とスケールを統合し、全所属面に対する係数の最小値を共有してください。その後、**SLAU・SFR 2・2次再構成**で seam／内部の平行移動試験を追加する必要があります。

4. **Minor — `--from-floor` が前提不成立を検出せず PASS を返す。**

   [check_convergence.py:120](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:120) は参照 run の通常判定を確認しません。また同ファイル149行では、参照床のない列を合格扱いでスキップします。参照でゼロだった成分が対象で非ゼロになっても検出できません。

   実行確認では `run_0088_passiveC_m2_rerun0067` は通常 **`NOT CONVERGED`**、同じ run を `--from-floor` の参照にすると **`PASS`、exit 0** でした。今回の正しい参照 `run_0476` による4本の PASS を否定する指摘ではありません。

   **対案:** 参照の通常 PASS、必要列の対応、ゼロ参照列の非ゼロ化を検査し、不成立なら判定を拒否してください。

5. **Minor — 最新の設計変更が現在仕様と残作業表に反映しきれていない。**

   [methods/condensation.md:834](/home/sano/work/forge/methods/condensation.md:834) は依然 `theta≥1e-12`・停止なしと説明していますが、指摘1の実装は `theta=0` を許します。最新の `φ_N·δρ`、`limCorr`、初回 sub-iter 限定処理も現在仕様に十分反映されていません。

   また [plan:177](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:177) の「物理 step あたりの上限」に対し、実装は [condensationTransport_d.cu:297](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:297) の初回 sub-iter だけの制限です。以後の増分を累積制限しないため、物理 step 全体の上限を保証しません。§5.1 #18 は未達なのに #8 は §6-6 済としており、完了状態も矛盾しています。

   **対案:** §4・§5.1・§10 と `methods/` を実装の契約に合わせて更新し、未達ゲートを #8 から参照してください。`plans/README.md` の `draft` 表記も同期が必要です。`methods/index.md` は今回ファイル追加・移動がなく、目次変更の欠落は認めません。

**推奨は、`in_progress` を維持し、この plan 内で指摘1〜3を解消してから result レビューを再実施することです。** 凝縮 dual-time は当初から目的・検証・完了条件に含まれており、現段階での後続課題化による accepted 移行には反対します。レビュー提案は **plan 未反映**です。指定どおりファイルは変更していません。

指摘数: Critical 0 / Major 3 / Minor 2
