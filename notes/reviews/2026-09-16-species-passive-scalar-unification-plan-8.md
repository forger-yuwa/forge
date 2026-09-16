# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `208be829` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: 8 回目 (plan-7 NO-GO M1–M6/m1 の反映; §4.7 v6, 実装・単体済 [commit 後])。変更: (1) check_passive_budget.py を総量比の合計 ≤ tol・有限性・FCT 記録必須 ([passiveFct] active 行)・総量 0 の補正は rel=1 で FAIL に (失敗系 tests/unit/test_check_passive_budget.py); (2) 境界の実現流束 Σ_∂ G Δt (符号付き) と落とした量、H = H_src + H_rem (物理ソース履歴 vs 残差・クランプ由来) の分離積算で収支を閉じる; (3) 実現可能性: 厳密 0/アンダーフローの特異不整合を不等式検査より先に単分散再初期化 ((0,0) も)、δ 撤廃、小さい正の状態は最近点 (連続); (4) EOS 後は Q1/Q2 だけの射影カーネル (g・消滅は触らない); (5) 線形残差と r_H を成分ごとに正規化し全成分で受入、区間最大を log; (6) E_ρ = M(ρ^{n+1}−ρ^n) − Bṁ^eff と Lρ − f を実作用素で評価し fctDensity に出す (逸脱は基点逸脱として FAIL に接続); (7) 正本式 A^raw = (1/a)F_H + (c/a)G^n − F_L(ṁ^eff, q_L), M(q_H − q_L) = B A^raw + r_L − r_H/a。問い: 残る論理の穴、収支の閉じ方 (非周期: Δ総量 = −Σ_∂ G Δt + Σ H_src Δt + Σ H_rem Δt) の妥当性、§6-2/§6-6 のゲート (check_passive_budget.py FAIL, analyze_moment_order.py) で失敗を検出できるか。§4.1–4.6, §4.8 の再指摘は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
流束形 BDF2 と `H_src/H_rem` 分離の方針は妥当です。  
しかし v6 でも、収支ゲートの偽 PASS、密度診断の境界項の誤り、最近点射影と時間次数ゲートの欠陥が残っています。

対象は HEAD `208be829`。目的・スコープは妥当で、既存 accepted plan の未解決事項 F-sp1／F-cf8 を引き継いでおり、重複ではありません（[F-sp1](/home/sano/work/forge/plans/accepted/thermophysics-cea-mole-fraction-species.md:186)、[F-cf8](/home/sano/work/forge/plans/accepted/condensation-source-limiter-steady.md:90)）。以下は指定された §4.7 と検証ゲートに限定します。

非周期収支
\[
\Delta Q=-\sum_{\partial}G\Delta t+\sum H_{\rm src}\Delta t+\sum H_{\rm rem}\Delta t
\]
は符号も含めて正しいです。ただし、`H := MΔq−BG` から作るため、**閉合そのものは解の正しさを保証しません**。独立した総量差との照合と、非物理的な `H_rem` の制限が必要です。ピン・クランプを `H_rem` に含めた場合、閉合式にそれらを再加算してはいけません。

1. **Major — 実ログの FCT 行を読み取れず、最終区間も検査していない**

   **根拠:** [check_passive_budget.py:21](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:21) は種名の後に空白1個を要求しますが、出力は [speciesTransport_d.cu:1462](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1462) の `%-8s` です。`roXi` の後には複数の空白が入ります。

   実際に `case/09.Taylor-Green/run_0134_passiveF_fct_step_seam_venkat100/forge_run.log` を読むと、floor 記録は1成分、**FCT 記録は0成分**。ツールは `FAIL(NO_FCT_RECORD,SUM>tol)` を返しました。ログには FCT 行が存在します。

   また、この100 step run の最終収支記録は **step 91** です。[main.cpp:2047](/home/sano/work/forge/solver_density_cuda/main.cpp:2047) は monitor 間隔でしか出力せず、ツールも完了 step との一致を検査しません。末尾9 step の補正を見落とします。

   **対案:** 空白を `\s+` で解析し、終了時に全成分の最終収支を必ず出力してください。完了 step・成分集合・記録時点を照合し、不完全な記録は拒否すること。現在の単体試験は `evaluate()` への辞書入力だけなので、**実際の出力書式から parse→判定まで通す試験**を追加してください。

2. **Major — 収支・有限性・補正総和の拒否条件が未完成で、偽 PASS が残る**

   **根拠:** [check_passive_budget.py:47](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:47) 以降は、境界流束・ソース・総増分を解析しても、閉合誤差を計算しません。有限性検査も一部の相対値に限定され、クランプ収支は他の補正と別枠で判定されます。

   現行 `evaluate()` を直接実行し、次の全例が **PASS** になることを確認しました。

   | 入力した異常 | 判定 |
   |---|---|
   | 総増分 `1`、境界・ソース・残りはすべて `0` | PASS |
   | 境界流束が `NaN` | PASS |
   | `qL` 相対残差が `1` | PASS |
   | 同じ成分の基点補正 `6e-7` ＋クランプ補正 `6e-7` | PASS |

   さらに [speciesTransport_d.cu:1486](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1486) は、**クランプ収支だけ総量0なら補正が非ゼロでも相対値0**を出します。「総量0の非ゼロ補正は拒否」の修正が全経路に届いていません。

   **対案:** 独立集計した初期・最終総量から閉合誤差を検査し、全生値の有限性、全期間の線形解受入、成分ごとの全補正合計を判定してください。クランプも同じ成分の合計に組み込み、総量0・非ゼロ補正は明示的に FAIL にするべきです。区間最大残差も最終区間だけでなく全期間を検査してください。

3. **Major — `fctDensity` は境界付きの実作用素と一致せず、FAIL にも接続されていない**

   **根拠:** [speciesTransport_d.cu:1771](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1771) は `nbRho` の計算に `qP=roN` を渡します。[passiveFct_d.cuh:138](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:138) の境界流入は、その結果 **`|ṁeff|×1`** になります。一方、実際のスカラ RHS は `|ṁeff|φⁿ` です。これを [同:437](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:437) で差し引くため、記録値は要求された `Lρ−f` ではありません。

   1 CV、`M=ρ=1`、流入・流出とも1、`φⁿ=0.5`、局所有効履歴項 `0.75` とすると、

   - 実際の低次解：`qL=(0.5+0.75+0.5)/2=0.875`
   - 正しい上限余裕：`Lρ−f=+0.25`
   - 現行診断：`−0.25`

   となり、許容解を逸脱と数えます。また、診断量は [speciesTransport_d.cu:1782](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1782) の別バッファに積算されるだけで、`baseViol` に加わらず、収支ツールにも parser がありません。

   **対案:** 内部近傍作用と境界の定数 RHS を分離し、同じ `L` と完全な `f` から評価してください。そのうえで正規化した逸脱量を拒否条件へ接続し、上記の許容例と不許容例を単体試験に追加してください。

4. **Major — 小さい正の状態に対する「最近点射影」が最近点になっていない**

   **根拠:** [condensationRealizability_d.cuh:25](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:25) は、常に区間 `[0,1]` を黄金分割で60回探索します。探索分解能が約 `1e-13` なので、計画が明示的に扱う小さい正の状態では精度が足りません。

   同じ演算を CPU で再現すると、

   ```
   入力: x=1e-14, y=1e-28×(1−2e-6)
   出力: x=1.444480187e-13, y=2.086523011e-26
   ```

   です。`x` は14.4倍、`y` は208.7倍になります。しかし既知の許容点 `(1e-14,1e-28)` までの距離は **約 `2e-34`**。現行出力までの距離は **約 `1.34e-13`** です。ゼロ特異点の修正は確認できましたが、「小さい正の状態には連続な最小補正」という契約は満たしていません。

   **対案:** 境界上の距離の停留条件を、局所スケールに応じた区間と相対停止条件で解いてください。単純な許容候補への距離を上界として保持し、それより悪い候補を採用しない保証も必要です。対数的に小さくした入力で、補正距離・成分相対変化を試験してください。

5. **Major — `analyze_moment_order.py` は §6-6 の合否ゲートを実装していない**

   **根拠:** [analyze_moment_order.py:44](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:44) 以降は次数を表示するだけで、BDF1/BDF2 の閾値判定がありません。実現可能性不等式、同一最終時刻、刻み比、全 step のサブ反復低下も検査しません。`FAIL` を表示した場合も非ゼロ終了しません。

   人工的な3水準データで全量の次数を **1.0**、サブ反復差を0にすると、`VERDICT: see orders above` で正常終了しました。

   実データでも、`case/44.vitiated_air_wt/run_0256_passiveD_order_bdf2_dt1.6e-5_nsub40/` を含む指定の3水準実行は、[同:16](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:16) が `res_wall_*.h5` まで拾い、`AttributeError` で停止します。

   境界出力を除いて既存 `run_0256`–`0259` を独立集計すると、次数は `g=1.303、Q1=1.470、Q2=1.448`、サブ反復差比はそれぞれ `0.431、0.377、0.397`。これは**既知の旧試験の未達値の再確認**であり、v6 の結果ではありません。こうしたデータを確実に拒否する回帰試験が必要です。

   **対案:** `res_<整数>.h5` のみ選び、BDF 次数・時刻・刻み比・必要全成分・有限性・実現可能性・サブ反復条件を検査し、明確な PASS/FAIL と終了コードを返してください。ケース選定は周期・凝縮という目的に合っていますが、現在のツールでは計画の条件を強制できません。

6. **Minor — 基点の式に `1/a` の脱落が残る**

   **根拠:** [plan:197](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:197) の正本式は修正済みですが、直後の198行と [passiveFct_d.cuh:13](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:13) は、依然として `qL+M⁻¹(rL−rH)` です。

   **対案:** 両方を `qL+M⁻¹(rL−rH/a)` に統一してください。実装の直接差分とは別に、残差から基点誤差を評価する際の誤りを防ぐ必要があります。

**推奨は、流束形 BDF2＋post-step FCT を維持し、検証ゲートを先に完成させてから統合検証を再開することです。** 優先順は、①ログ読込・最終記録と収支拒否、②境界付き密度診断、③小スケール射影、④時間次数・実現可能性ゲートです。これらを §5.1 #21 と §6 の必須条件にしてください。

既存 Python 単体試験は `ALL PASS` でしたが、上記の反例を含んでいません。再実行した `check_convergence.py` は `run_0134` と `run_0257` とも **NOT CONVERGED**。この通常判定を、非定常計算の各物理 step 内の収束判定には代用していません。[case/09 索引](/home/sano/work/forge/case/09.Taylor-Green/README.md)、[case/44 索引](/home/sano/work/forge/case/44.vitiated_air_wt/README.md)。

ファイル変更・新規 CFD 計算・GPU 単体試験の再実行は行っていません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
