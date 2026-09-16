# codex レビュー: species-passive-scalar-unification (result)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `result` (diff base `3b62903a`)
- **date**: 2026-09-16
- **commit**: `9549f490` (feature/sern-design)
- **codex**: effort `high`, 8.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m2
- **focus**: 実装 diff (受動種の化学種経路化 Phase A: passiveTransport_d.cuh / passiveKernels_d.cuh / passiveLimiter_d.cuh / speciesTransport_d.cu / convectiveFlux_slau_d.inc.cuh / condensationTransport_d.cu / tracerTransport_d.cu / periodicNode_d.cu / limiter_d.cu / variables.cpp / main.cpp / solverConfig; Phase B: dual-time の化学種移植 + 受動種 BDF + 共有履歴契約 + /CHECKPOINT) が plan §4.1–4.4 (codex plan レビュー 1/2 回目の採否表 §6.1 の対応を含む) と一致するか。検証 (§9 の 2026-09-17 の各項; case/16 README run_0473–0496, case/44 README run_0204–0232, case/09 README run_0067–0086, case/28 README node 表, case/34 README passA ラベル, case/46 README run_0107–0113) が §6 1–7 のゲートを満たすか: 原因確認 (S3 発散は coupling 0 固有 / 1.8e-4 は再正規化)、制御試験の改定ゲート、非一様組成 PASS ケースの固定点不変 (交差 restart は収束場からなので drop 基準では NOT CONVERGED)、解析解拡散 2 次、等拡散一致 (コア 2.7e-6)・混合平均非一致、凝縮回帰で S3 が固定点を動かす (onset 下流) こと、dual-time の 3 水準時間次数 (BDF2 2.03 / BDF1 0.91)・restart・ピン行・過渡凝縮 (nSub 頭打ちの床)、無影響、既定変更 (passiveScalarScheme 1, SFR 0 のまま)。cell はユーザ指示で対象外。accepted に移してよいか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

基本経路の統一と BDF 移植は確認できましたが、周期 DPLUR・補正収支・checkpoint の実装に不整合があります。凝縮の時間精度と累積補正量も合格条件を満たしておらず、現状で `accepted` へ移すことには反対します。

1. **Major — 周期 DPLUR が、合併 CV の近傍補正を解いていません。**

   **根拠:** [speciesTransport_d.cu:348](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:348) は各ノード自身の `cell_planes` だけで `neighbor` を組みます。その後は [同:851](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:851) と [同:1273](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1273) で root の `dq` をコピーするだけです。メッシュ側も体積を合併する一方、面リストは各部分 CV のままです（[mesh.cpp:710](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:710)、[同:839](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:839)）。

   したがって、残差・対角は合併済みなのに、非対角項は root 側だけになります。状態の周期一致は保証できても、内部と同じ更新作用素にはなりません。有限サブ反復での精度・安定性に影響します。

   `case/09.Taylor-Green/run_0067–0068` は RK、`run_0080–0086` は化学種 coupling 0 です。case/16 の coupling 1/2 試験は非周期で、要求した非一様組成の周期 coupling 1/2 検証を代替できません。

   **対案:** 各 sweep の近傍寄与を独立バッファに作り、周期群で合算してから解く。非一様組成・トレーサが seam・辺・角を横切る陰解法試験で、内部との更新速度と保存収支を確認してください。

2. **Major — 周期境界の補正収支が重複計上されています。**

   **根拠:** [passiveKernels_d.cuh:44](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:44) は「二重計上されるが相対量は不変」としていますが、これは誤りです。[speciesTransport_d.cu:1287](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1287) が合併体積を渡し、全メンバーで積分しています。内部・面・辺・角で重複倍率が異なるため、局在する補正と総量の比では相殺しません。

   実測でも `case/09.Taylor-Green/run_0067_passiveA_tracer_seam_1st/res_500.h5` の総トレーサ量は、ログの **166.246** に対し、周期重複を除いた積分では **148.830** でした。

   **対案:** 保存量と補正量の積分は、周期 root のみを合併体積で数える方式に統一する。補正を seam に局在させた試験で符号付き・絶対収支を検証し、既存の保存・補正判定を再集計してください。

3. **Major — checkpoint の「履歴が有効」という判定が不足しています。**

   **根拠:** [main.cpp:1024](/home/sano/work/forge/solver_density_cuda/main.cpp:1024) は `dt` が違っても警告だけで履歴を復元し、[同:1705](/home/sano/work/forge/solver_density_cuda/main.cpp:1705) は固定刻み BDF2 係数を使います。刻み変更 restart は、その最初の物理 step で誤った時間微分になります。

   また、[output.cpp:163](/home/sano/work/forge/solver_density_cuda/output/output.cpp:163) は scheme 0 の凝縮でも履歴を書きますが、[speciesTransport_d.cu:1425](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1425) は scheme 0 の受動種履歴をシフトしません。`layout` に scheme の識別がないため、scheme 0 の checkpoint を scheme 1 で読むと、古い初期値を有効な物理履歴として受け入れます。

   **対案:** 刻み・履歴を生成した方式の整合を復元条件に加え、不一致時は全系を BDF1 に戻す。非定常の組成・トレーサ・モーメントについて、連続実行、通常 restart、刻み変更、scheme 変更、履歴欠落を検証してください。現在のほぼ定常・化学種のみの restart 試験では検出できません。

4. **Major — 凝縮 dual-time は §6-6 の時間精度ゲートを通っていません。**

   **根拠:** CSV を再集計すると、`case/44.vitiated_air_wt/run_0230_passiveB_dt_cond_transient_nsub40/` のサブ反復残差低下中央値は、`roQ0` **1.584 桁**、`roY0` **1.922 桁**、`roY1` **1.770 桁**です。要求する「各物理 step・全化学種／受動種で ≥2 桁」を満たしません。

   HDF5 の再比較でも、`run_0229` と `run_0230` の差は `g` **3.74e-5**、`Q0` **4.08e-4**。一方、`run_0226` と `run_0227` の刻み半減差は **9.52e-5 / 4.64e-4**です。サブ反復依存が十分小さいとはいえず、凝縮について第3の刻み水準もありません。

   ガウス移流 `case/09.Taylor-Green/run_0080–0086` の **BDF2 2.033 次、BDF1 0.915 次**は再現できました。ただしこれは凝縮モーメント・ソースの時間精度の証明ではありません。「頭打ち」は確認できても、その原因を CFL・θクランプに確定する比較も不足しています。

   **対案:** 凝縮過渡でサブ反復誤差を分離し、モーメントを含む3水準の次数試験を完了する。§10 の後続課題に移して完了扱いにせず、§5.1 の必須残作業へ戻してください。

5. **Major — 累積 floor 補正が、採用済み上限 `1e-4` を超えています。**

   **根拠:** [plan:104](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:104) の全期間相対積算上限に対し、実ログは次の値です。

   | run | 累積 `abs/total` |
   |---|---:|
   | `case/16.nozzle_wys/run_0476_tp5_node_euler_profile_base/` | `roXi`: **2.216e-4** |
   | `case/16.nozzle_wys/run_0494_passive_eqdiff_constSc/` | `roXi`: **5.029e-4** |
   | `case/34.arthur_n2_nozzle/run_0106_limiter_regress_n2_node/passA_s1_sfr2/` | `Q2`: **1.241e-3**、`Q1`: **2.688e-3**、`Q0`: **5.609e-3** |

   これらは周期重複とは別の超過です。「過渡だけ作用し、後半はゼロ」は固定点近傍の説明にはなりますが、全期間ゲートの合格理由にはなりません。

   **対案:** 起動・更新条件を改めて上限内に収め、補正込みの収支を再検証する。過渡補正を許容する別基準を採るなら、結果への影響を定量化して plan の設計判断から再レビューすべきです。

6. **Major — 固定点試験で、化学種の緩和変更が実際には試されていません。**

   **根拠:** `case/16.nozzle_wys/run_0479_tp5_node_euler_profile_xr_relax10/solverConfig.yaml:21` は `speciesImplicitCoupling: 1`、`implicitRelax: 0.7` のままです。変更した `speciesImplicitRelax: 1.0` は coupling 0 専用で、coupling 1 は [speciesTransport_d.cu:829](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:829) の `cfg.implicitRelax` を使います。基準でも `speciesImplicitRelax` は既定 1.0 なので、化学種側は比較になっていません。

   再実行した判定は、基準 `run_0476` が **`PASS (converged)`**、交差 restart `run_0478` が **`NOT CONVERGED (stalled/plateau)`**でした。restart が既に低残差から始まる説明は妥当ですが、§6-4 の「全て PASS」とは未整合です。また [case/44 README:721](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:721) は `NOT CONVERGED` の準定常回帰を「固定点は同一」と表現しています。

   **対案:** 実際に効く緩和係数を変更した試験を追加する。収束場 restart は、基準 PASS の全残差・場の変動幅を使う判定をツール化して plan に明記する。case/44 は、確認済みの **`OVERALL: ALL STEADY`** に基づく「準定常量の比較」に表現を限定してください。

7. **Minor — 現在仕様と運用手順に旧仕様が残っています。**

   **根拠:** [solver-settings.md:198](/home/sano/work/forge/procedures/solver-settings.md:198) は tracer を「拡散なし・dual-time 禁止」、[recommended-settings.md:124](/home/sano/work/forge/procedures/recommended-settings.md:124) は凝縮 dual-time を自動降格と説明しています。[thermophysics.md:278](/home/sano/work/forge/methods/thermophysics.md:278) も受動種更新を point-implicit のままと記載しています。

   **対案:** scheme 0/1 の条件、DPLUR 自動選択、dual-time 対応を現在仕様として統一する。既定変更自体は実装・plan に記載されていますが、追記だけでは矛盾が解消していません。

8. **Minor — 残作業表と run 索引が、実際の検証状態に同期していません。**

   **根拠:** [plan:173](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:173) は「済」と未達を併記し、次行では §6-6 を完了扱いにしています。§10 の ψ 切替比較も未決のままです。さらに [case/09 README:40](/home/sano/work/forge/case/09.Taylor-Green/README.md:40) は時間次数試験を改番前の `run_0072–0078` と記載し、拡散試験の番号と衝突しています。case/16 の restart 説明にも改番前番号が残っています。

   **対案:** 未達ゲートを優先順付きで §5.1 に戻し、索引を実在する run と成果物に統一する。`check_plans.py` は今回 **PASS** でしたが、これは構造検査であり検証完了の裏付けにはなりません。

**推奨は、`active` に留めて修正・再検証することです。** 優先順は、周期更新・収支と checkpoint の修正 → 凝縮時間精度・補正上限の達成 → 有効な緩和比較 → 文書・残作業表の同期です。既定変更を「全ゲート達成済み」として承認する根拠はありません。

ファイルは変更していません。本レビューの指摘は plan 未反映です。

指摘数: Critical 0 / Major 6 / Minor 2
