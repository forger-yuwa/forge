# codex レビュー: condensation-source-limiter-steady (result)

- **plan**: [`plans/active/condensation-source-limiter-steady.md`](../../plans/active/condensation-source-limiter-steady.md)
- **stage**: `result` (diff base `f1e40f41`)
- **date**: 2026-09-15
- **commit**: `357905c1` (feature/sern-design)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m1
- **focus**: 実装 diff (θ の残差撤去・蒸発の瞬間速度形 cond_evap_source_rate・更新クランプ kernel condensationUpdateLimiter_d.cuh・RK/dual-time 自動降格) が plan §4 と一致するか、§6 の検証 run (case/44 run_0132–0136, case/34 run_0105–0107, case/16 run_0470; 各 case README の該当行と plan §9) の実測が変更ログの主張を裏付けるか、accepted にしてよいか。単体テストは solver_density_cuda/tests/unit/test_cond_limiter_steady.cu

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

# 判定: **NO-GO**

残差から θ を外す主要変更は確認でき、報告された最終場の差分値も概ね再現しました。  
ただし、検証ゲート未達、補正量診断の欠落、plan と異なる蒸発式、平衡モデルへの意図しない適用があります。  
現状で「定常解の一致を検証完了」として `accepted` に移すことは認められません。

## 指摘一覧

### 1. **Major — 「検証完了」を裏付ける収束・定常性ゲートが成立していない**

**根拠:** 対象全 run に `check_convergence.py` を再実行し、保存済み判定と同じ結果を確認しました。

| 対象 | 再確認した結果 |
|---|---|
| `case/44.vitiated_air_wt/run_0132`–`0136` の対象 run | 全て **NOT CONVERGED** |
| `case/34.arthur_n2_nozzle/run_0105`–`0107` の `lim0` / `lim1`、`0105/lim1_r2` | 全て **NOT CONVERGED** |
| `case/16.nozzle_wys/run_0470_limiter_regress_wys/lim0` / `lim1` | 両方 **NOT CONVERGED** |

特に、[run_0132 の判定](/home/sano/work/forge/case/44.vitiated_air_wt/run_0132_va3_M4.19_Lc8_noneq_inletTt_lim1_cfl2/CONVERGENCE_VERDICT.txt:1)では `rms_roQ1_0` が **2.9 桁**、`rms_roQ0_0` が **2.7 桁**の低下です。`run_0135` も同様で、[plan §6](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:127)の代替ゲート「全モーメント残差が 3 桁以上低下」も満たしません。

Wysłouzil は保存場が step 0 と 48000 の **2 枚だけ**です。`check_quasisteady.py` の再実行結果は両モードとも **TRANSIENT-UNSETTLED**。報告量の定常性を確認できません。

数値そのものの再現性はあります。

- `run_0136` 対 `run_0132`: g 相対 L1 **0.002099**、最大温度差 **0.22853 K**、最大 Mach 差 **0.0022006**。
- Wysłouzil `lim1` 対 `lim0`: 正規化最大差は `ro` **1.91e-6**、g **9.89e-6**。
- Arthur の `perf_regress.py` は記載どおり。ただし `run_0105/lim1` は **FAIL（27/28）**、`ro` 差 **6.18e-4**、ノイズ床 **2.90e-4**です。これは収束判定とは別の比較です。

**対案:** 「最終保存場の差が小さい」と「定常解が一致した」を区別し、§5.1 の検証項目を未完了へ戻してください。全残差ゲートと、報告する全派生量の **STEADY** を満たしてから再判定するべきです。`lim1_r2` の成功だけで最初の FAIL を合格扱いせず、事前に定めた反復比較基準で評価してください。

### 2. **Major — `condClampCorr=0` では硬クランプ無作用を証明できない**

**根拠:** [condensationUpdateLimiter_d.cuh:73](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:73)では全モーメントを非負へ補正しますが、78 行で記録するのは **負の `rog` を 0 に戻した量だけ**です。

以下の補正は記録されません。

- `roQ0/1/2` の負値補正。
- 後段の `rog ≤ roY_w` / `0.99ρ`。
- 液滴消滅による 4 モーメントの削除。

後二者は [condensationTransport_d.cu:83](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:83)、[同:114](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:114)で実行され、診断配列を受け取ってすらいません。

したがって、実測された `condClampCorr=0` は、plan が要求する「補正量の総和がゼロ」を意味しません。残差が残っていても、未記録の射影が更新を打ち消す可能性を排除できません。

**対案:** 候補更新から最終確定までの**全補正**を、モーメント別の単位を保って記録してください。上限超過・各モーメントの負値・消滅を個別に試験し、その診断で検証ゲートを再評価してください。

### 3. **Major — 更新クランプの評価状態が plan と異なる**

**根拠:** [plan §4.2](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:66)は、密度変化を含む

`Δg = (δ_g − g_old Δρ) / ρ_new`

と更新後の温度を要求しています。一方、[kernel:43](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:43)は `g_old=N_g/ρ_new`、`Δg=δ_g/ρ_new` です。

例えば `ρ_old=1, g_old=0.01, ρ_new=0.5, δ_g=0` なら、実際の g は **0.01 増加**しますが、この kernel は `Δg=0, θ=1` と判断します。

さらに、[main.cpp:1360](/home/sano/work/forge/solver_density_cuda/main.cpp:1360)で流れの保存量を更新してから、[同:1410](/home/sano/work/forge/solver_density_cuda/main.cpp:1410)で limiter を呼ぶまで温度・比熱の再計算がありません。新密度・新組成と旧温度・旧比熱が混在します。

**対案:** 実際の `g_new(θ)−g_old` と、整合した EOS 評価状態で制限してください。モーメント増分を縮めるだけでは制限を満たせない密度変更も扱う必要があります。密度・組成が同時に変化する試験を追加し、評価状態の定義を plan と実装で一致させてください。

### 4. **Major — 蒸発ソースは plan の瞬間速度式を実装していない**

**根拠:** [plan:59](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:59)は、保存量表記へ合わせると

- `S_Q1 = q0 drdt`
- `S_Q2 = 2 q1 drdt`
- `S_g = 4πρ_l q2 drdt`

です。しかし [condensationSource_d.cuh:335](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cuh:335)は `a=drdt/r30` として `a q1`, `2a q2`, `3a ρg` を使います。float 版も同じです。

これは旧 λ スケーリングの微分であり、plan の式と一致するのは単分散の場合だけです。例えば半径 `r` と `2r` の液滴を同数含む分布では、両式の `S_g` は約 **9 %**異なります。Δτ 不変性だけでは、このモデル差を検出できません。

**対案:** 本 plan では §4.2 の式に揃えてください。多分散の実現可能なモーメントを使い、ソース値そのものを検証する試験と、蒸発域の回帰が必要です。

### 5. **Major — 据え置き予定の `condEquilibrium=1` に新クランプが適用される**

**根拠:** [condensationTransport_d.cu:342](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:342)の新経路条件は `timeIntegration==11 && condLimiterMode==1` だけで、`condEquilibrium` を除外していません。

平衡緩和ソースは従来の `condEqDTmax` などで制限済みですが、その後さらに新しい **`condDTmaxStep=1 K`** の更新クランプを受けます。従来と更新量・反復経路が変わり、[methods:820](/home/sano/work/forge/methods/condensation.md:820)の「据え置き」と矛盾します。平衡形でも更新演算が float から double 中間演算へ変わります。

**対案:** 新経路を **`condEquilibrium==0` に限定**し、平衡形は従来の更新へ戻してください。設定省略時も含む平衡形の回帰試験を追加してください。

### 6. **Major — 周期 node seam で凝縮ソースを二重計上する既存欠陥が残る**

**根拠:** ソース wrapper は両精度で [condensationSource_d.cu:64](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cu:64)、85 行の `volume` を渡します。しかし周期 node の `volume` は [mesh.cpp:710](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:710)で**合併体積**に置換済みです。

各 member が `S·Vmerged` を加算した後、[periodicNode_d.cu:90](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:90)で凝縮残差も合算します。2 member の seam では、本来 `S·Vmerged` のところ **`2S·Vmerged`**になります。

これは今回の差分で新設されたバグではありません。ただし、依頼された境界監査では未解決の重大欠陥であり、現在の非周期ノズル検証では検出できません。

**対案:** ソース積分には `volumePartial_d`、更新分母には合併体積を使ってください。周期 seam と内部で同一状態・同一ソースになる試験を追加するべきです。

### 7. **Major — 必須試験の欠落を残作業表が完了扱いしている**

**根拠:** [plan:107](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:107)が要求する「輸送とソースが非ゼロで釣り合う 1 セル固定点試験」は、[test_cond_limiter_steady.cu](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:176)にありません。密度・組成の前後変化を扱う試験もありません。

CFD も、case/44 の `condFloat=0` は CFL 2 の一本だけです。node/cell × 両精度それぞれの CFL 比較にはなっていません。Arthur cell の旧新比較は、大型ノズルで limiter が強く作動する条件の CFL 不変性試験を代替しません。

**対案:** §5.1 の単体・回帰を再開し、未実施項目を具体的に残してください。特に固定点試験は実際の硬クランプまで通す必要があります。

### 8. **Minor — 文書・残作業の同期が未完了**

**根拠:**

- [plan §4.4](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:94)は、降格する RK にも更新クランプが効くと記載。
- [plan §7](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:154)と [case/44 README](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:584)は、§6 で参考扱いした `run_0131` を固定点・正本と記載。
- [recommended-settings.md:117](/home/sano/work/forge/procedures/recommended-settings.md:117)は `condEquilibrium=2` を既定と記載していますが、コード既定は 0。
- F-cf8 は [followups:60](/home/sano/work/forge/plans/active/condensation-followups.md:60)に登録済みですが、対象 plan の残作業表は未完了のまま。旧 mode 0 の削除課題は followups に未登録です。

**対案:** 実装・検証状況を確定してから、plan、現在仕様、推奨設定、case 台帳を同期してください。`methods/index.md` は今回ファイル追加がないため、変更不要という判断で問題ありません。

## 推奨

**`in_progress` を維持し、実装と診断を修正したうえで、未加工残差・全補正量・派生量時系列による検証をやり直してください。** 最終場の近さは確認できましたが、現時点では受入条件を満たしません。

レビュー対象は `f1e40f41...HEAD`、実装 commit は `357905c1`。ファイルは変更していません。CUDA 単体テストの再実行は driver/runtime 不整合で実施不能でした。これはコードの試験失敗とは区別します。本レビューの指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
