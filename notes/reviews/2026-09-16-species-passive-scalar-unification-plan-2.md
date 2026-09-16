# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `e434b1f7` (feature/sern-design)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m0
- **focus**: 2 回目 (1 回目 C1/M7/m1 反映後; §6.1 の採否表)。前回指摘の再指摘は不要で、改訂後の §4.0–4.4 と §6 1–8 が実装可能かつ十分か: (1) 残差は完全な R2、緩和は増分と擬似 Δτ のみ (§4.2) で固定点不変が成立するか; (2) 受動種ごとの Venkat ψ_P + 更新後 floor + floor 補正診断で有界性の扱いが妥当か; (3) 周期 node の処理順 (勾配の周期除外・transport_diag gather・BDF 一度・更新後 mirror); (4) dual-time の処理順 (BDF 後のピン再除去、coupling 0/1/2 分岐、履歴初期化・restart); (5) 検証ゲート (段階比較、拡散の解析解、3 水準時間精度、固定点不変は PASS ケース、凝縮回帰は準定常) の実効性と抜け。cell は対象外。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **GO-with-changes**

完全な `R2` と増分緩和を採る方向は妥当です。ただし、float32 のモーメント再構成、上限保証、周期の結合更新、固定点検証に実装前の修正が必要です。  
以下は改訂後に残る問題です。ファイル変更・新規計算は行っていません。

## 指摘一覧

### 1. **Major — 既存 Venkat 式をモーメントへ直接流用するとオーバーフローする**

**根拠:** [limiter_d.cu:34](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:34) は、分子・分母に `delta_p² * delta_m` などの三次積を float32 で計算します。

`case/44.vitiated_air_wt/run_0170_va3_M4.19_Lc8_noneq_inletTt_lim1e_cfl2/res_24000.h5` を読み直すと、`Q0_0` の最大値は **8.670266e14**。このスケールの差分を扱うための対策が §4.1 にありません。

現行式をホストの float32 で再現すると、`delta_p=1e14, delta_m=1e13, eps²=1e-5` で **分子・分母とも Inf、商は NaN** になります。同じ式の float64 評価は有限です。CUDA の後段 `min/max` が NaN を隠しても、正常なリミッタ評価にはなりません。

**対案:** Venkat 関数自体は維持し、差分と `eps²` を整合してスケーリングした、オーバーフローしない評価式にすること。`g/Q2/Q1/Q0` の実測スケール、ゼロ近傍、極値で、クリップ前の係数まで有限性を検査する試験を追加してください。

### 2. **Major — 更新後 floor だけではトレーサの上限を保証できない**

**根拠:** [plan:93](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:93) は更新後を非負化のみにしますが、[plan:158](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:158) は `roXi/ro ≤ 1` を要求しています。実際の [scalarTransport_d.cu:290](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:290) は下限しか処理しません。

例えば更新前に `roXi=ro=1`、流れ更新後に `ro=0.9`、スカラ増分がゼロなら、非負化後も `roXi/ro=1.111…` です。面クリップ・`ψ_P` はこの密度更新との不整合を防ぎません。化学種には再正規化がありますが、受動種にはありません。

診断にも相違があります。[condensationUpdateLimiter_d.cuh:78](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:78) の `condClampCorr` は更新ごとにリセットされる正規化量であり、時間を通じた保存量補正の積算ではありません。

**対案:** 新経路の更新確定時に、更新済み密度を使って **`0 ≤ roXi ≤ ro`** を適用する仕様にしてください。上下限それぞれの符号付き補正と絶対補正を体積積分し、物理ステップ内・全期間で記録する。保存誤差 `1e-6` は相対値として定義し、補正量にも合否閾値を設けてください。

### 3. **Major — 周期処理の改訂に coupling 1/2 の補正経路が含まれていない**

**根拠:** §4.4 の順序には coupling 2 の **予測・EOS クロス項注入**が明示されていません。現行の [main.cpp:1388](/home/sano/work/forge/solver_density_cuda/main.cpp:1388) は流れ block 更新前に予測を行い、[speciesTransport_d.cu:445](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:445) が `res_roe` に面ごとの補正を加えます。

これを既に周期 gather 済みの残差へ直接加えると、その追加分は部分 CV のままです。残差全体を再 gather すると、今度は既に合算済みの空間残差・BDF 項を重複加算します。

また、化学種 DPLUR の [speciesTransport_d.cu:846](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:846) は sweep 後に swap するだけです。流れ側の [main.cpp:1347](/home/sano/work/forge/solver_density_cuda/main.cpp:1347) にある周期補正同期に相当する処理がありません。更新後の受動種 mirror だけでは、この経路を扱えません。

**対案:** 処理順に次を追加してください。

- 現サブ反復の擬似刻み確定後、BDF 込み・ピン除去済み残差で coupling 2 を予測する。
- EOS クロス項は独立バッファに組み、その追加分だけ周期 gather して流れ RHS に加える。
- coupling 1/2 の近傍補正と各 sweep の `dq` を周期グループで整合させ、化学種の更新後状態も同期する。

周期試験は受動種だけでなく、**非一様組成を持つ coupling 0/1/2** を対象にしてください。

### 4. **Major — `run_0471` のままでは固定点検証が実質的に空振りする**

**根拠:** [plan:162](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:162) が指定する `case/16.nozzle_wys/run_0471_tp5_node_euler_va3/` は、トレーサ・凝縮なしのほぼ一様組成です。

`res_12000.h5` の再集計では、全域の `max(Ys)-min(Ys)` は **2.89e-8～6.56e-7**。S3 に切り替えても、検証したい組成再構成の寄与が丸め程度に留まります。

今回のツール再判定は以下です。

- `check_convergence.py`: **OVERALL: ALL PASS**
- `check_quasisteady.py --quantity machmax,pmax`: **OVERALL: ALL STEADY**

この判定は既存ケースの妥当性を支持しますが、非一様な受動スカラの固定点不変を裏付けません。索引は [case/16 README:373](/home/sano/work/forge/case/16.nozzle_wys/README.md:373) です。

**対案:** このケースは無影響回帰に残し、固定点ゲートには、**定常でも非ゼロ勾配が残る化学種・トレーサ・モーメント輸送ケース**を追加してください。`R2−R1` が反復ノイズより十分大きいことを確認し、緩和・擬似 CFL の交差 restart、全残差 `PASS`、補正無作用を要求するべきです。

`Q←Q+ωM⁻¹R2` の固定点不変は、正の緩和率と非特異な更新作用素、制約補正が無作用という条件付きです。floor や再正規化込みの更新停止だけでは証明できません。

### 5. **Major — 緩和の切替仕様と無影響・一致ゲートが両立していない**

**根拠:** [plan:99](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:99) は `passiveScalarScheme=1` のとき化学種の segregated 更新にも `implicitRelax` を掛けます。一方、[plan:171](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:171) は化学種のみの run が同フラグに依存しないことを要求しています。

現行 coupling 0 は [speciesTransport_d.cu:774](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:774) から緩和なしの更新へ進むため、`implicitRelax=0.7` では変更前後の写像が異なります。既存 `run_0471` の設定が `1.0` なので、そのケースだけでは矛盾を検出できません。

さらに S3 では、化学種は `ψ_ρ` と面正規化、受動種は `ψ_P` と独立クリップです。同じ更新方式・等拡散係数だけでは §6-1 の `1e-6` 一致条件になりません。

**対案:** 化学種の segregated 緩和を受動種の切替から分離し、独立した設定で明示してください。無影響試験には `implicitRelax=0.7` を含める。化学種対トレーサの厳密一致ゲートは、**一次移流を含め全離散作用素と更新写像を揃えた制御試験**に限定し、S3 の異なるリミッタ間には適用しないでください。

### 6. **Major — restart の BDF 履歴契約と検証が未完了**

**根拠:** [plan:120](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:120) は履歴欠落時の初期化を定義していますが、履歴がある場合の復元条件、物理時刻・刻み・履歴有効数の対応を定義していません。

移植元 `e296f0d0` の `speciesTransport_d.cu:827,846` はプロセス内の `g_speciesLevelShifts` で BDF1/BDF2 を決めます。一方、流れは [main.cpp:1560](/home/sano/work/forge/solver_density_cuda/main.cpp:1560) の `iStep` で判定します。配列名の互換だけでは、全保存量が同じ時間離散を使う保証になりません。

§6-6 の三水準比較には restart がなく、化学種については「`roY` が動く」確認に留まっています。誤った BDF 項でも `roY` は動きます。

**対案:** 全保存量で共有する履歴有効数と BDF 係数を定義してください。完全な checkpoint は流れ・化学種・受動種の履歴、時刻、刻みをまとめて復元し、不完全なら全系を揃えて BDF1 から再開する。

検証には、連続実行対 restart、旧形式入力、`nSpecies==1`＋受動種を追加してください。時間次数は非一様な化学種も含めて評価し、三水準の差から求める次数とサブ反復誤差の許容値を数値で定義する必要があります。

## 推奨

**現行の統合方針を維持し、上記を実装前に plan へ反映することを推奨します。** 目的は未解決の F-sp1/F-cf8 に対応しており、重複実装ではありません。受動種を熱力学・化学種正規化から分離する構成も妥当です。

修正優先順は **① float32 再構成 → ②上下限・補正収支 → ③周期結合処理 → ⑤切替契約 → ⑥履歴契約 → ④有効な固定点ゲート**。その後、予定どおり node の原因確認を先行させてください。

読み取り専用レビューのため、**plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 0
