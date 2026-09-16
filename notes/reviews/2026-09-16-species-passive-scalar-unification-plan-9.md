# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `dfdc102b` (feature/sern-design)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m1
- **focus**: 9 回目 (plan-8 NO-GO M1–M5/m1 の反映 + 検証 1 巡目の不具合 5 件; §4.7 v8, 実装・単体済, commit 済)。変更: (1) check_passive_budget.py: \s+ 解析・終了時の収支出力・最終 step 照合・閉合 (増分 = −境界流束 + ソース履歴 + 残り) と独立照合 (最終−初期 = 増分)・全生値の有限性・全期間最大の低次/HO 残差・同一成分のクランプ合算・総量 0 は rel=1、実書式 parse→判定の失敗系試験; (2) 上限診断を内部近傍作用素 Lρ と境界定数 RHS 込みの f_full で評価し upperViol として合計へ (単体 (f): 許容 +0.25 / 不許容 0.2); (3) 射影は単純クランプ候補を上界に局所スケール区間・相対停止 (x=1e-14 で補正 2e-6 級); (4) analyze_moment_order.py を res_<整数> のみ・時刻/刻み比/有限性/次数/nSub 比/sub-iter 低下/実現可能性の PASS/FAIL に (旧 run_0256–0259 は FAIL); (5) 検証 1 巡目の不具合: 射影の r=0 で x=∞ → NaN (非有限ガード), restart layout が読み側で passiveFct=0 (設定だけの passiveFctConfigured), KEEP で S3/FCT が黙って非作動 (WARNING), θ_b が FCT と併用で −0.7 % (FCT 有効時は θ_b とトレーサ sub-iter floor を切る), 射影が sub-iter 内で数百万回 (step 末尾のみに; condRealizProject)。問い: (a) FCT 有効時に θ_b とトレーサの sub-iter floor を切る判断 (モーメントの floor は残す) の妥当性; (b) 射影を step 末尾だけにする判断 (定常は毎 step); (c) 収支ゲートの閉合・独立照合の定義; (d) 残る穴と §6-2/§6-6 のゲート。§4.1–4.6, §4.8 の再指摘は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
`θ_b`／トレーサ floor の sub-iter 内解除と、射影の step 末尾への移動は妥当です。  
しかし、事後拒否を担う収支・時間次数ゲートに偽 PASS と誤判定が残り、現状のゲートを根拠に先へ進めません。

対象は HEAD `dfdc102b`。目的・スコープは妥当で、既存の [F-sp1](/home/sano/work/forge/plans/accepted/thermophysics-cea-mole-fraction-species.md:186)／[F-cf8](/home/sano/work/forge/plans/active/condensation-followups.md:60) の未解決事項を引き継いでいます。周期箱・凝縮ノズルの選択も検証目的に合っています。指定どおり、以下は §4.7 と §6-2／§6-6 に限定します。

1. **Major — 独立収支照合の初期・最終状態が、増分積算と揃っていない**

   **根拠:** [speciesTransport_d.cu:1463](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1463) は、最初の monitor 時点、つまり **step 1 後**の総量を `initialTotal` にします。一方、`increment` は最初の物理 step から積算します。

   `case/44.vitiated_air_wt/run_0278_passiveG_smoke_order3/` を実測すると、`rog_0` は次の状態です。

   | 量 | 値 |
   |---|---:|
   | `res_0.h5` の初期総量 | `1.6118299086e-4` |
   | log の `initialTotal` | `1.611792e-4` |
   | HDF5 の最終−初期 | `−8.60769362e-8` |
   | log の積算 `increment` | `−8.607694e-8` |

   **正しい初期値なら独立照合誤差は総量比 `2.36e-11`**ですが、現行ツールは `TOTAL_VS_INCREMENT(2.3e-05)` と誤判定します。

   また、総量は [main.cpp:1855](/home/sano/work/forge/solver_density_cuda/main.cpp:1855) の floor 時点で集計され、その後の消滅・射影を含む確定状態では再集計されません。

   **対案:** 初期総量は積分開始前、最終総量は全後処理後の確定保存量から、周期 root のみで独立集計してください。収支の始点・終点を完全に揃え、十分な出力桁数で保存するべきです。

2. **Major — 収支ゲートに、記録欠落・非有限値・保存誤差の偽 PASS が残る**

   **根拠:** [check_passive_budget.py:36](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:36) は、FCT／クランプ行を step と結び付けず、最後に読めた辞書を保持します。[同:77](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:77) の有限性検査も floor の生値を含みません。クランプ記録は欠けていても省略されます。

   現行の実書式テスト用入力を使い、parse→evaluate で以下を再現しました。

   | 異常入力 | 現行判定 |
   |---|---|
   | floor の絶対補正が `NaN`、相対値は `0` | PASS |
   | モーメントの `clampBudget` が欠落 | PASS |
   | step 91 の FCT 行だけ存在し、step 100 は floor 行のみ | PASS |
   | 境界・ソース・増分がすべて `0`、総量だけ `5e-6` 増加 | PASS |

   最後の例は [同:98](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:98) が独立照合だけ `1e-5` を許すためで、§6-2 の保存 `1e-6` を強制できていません。

   solver 側でも [speciesTransport_d.cu:1776](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1776) の `std::max` に NaN を渡すと、有限な既存最大値が残り得ます。run-max が有限でも、全期間有限だった証拠にはなりません。

   **対案:** 必須成分集合と step を持つ完結した収支レコードにし、欠落・解析失敗・全生値の非有限を拒否してください。solver 側にも非有限発生を保持するフラグが必要です。独立照合は出力精度を改善して `1e-6` に統一し、上記の失敗例を試験へ追加してください。

3. **Major — 時間次数ゲートが、必要な sub-iter 検証を実施せず PASS を返す**

   **根拠:** [analyze_moment_order.py:41](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:41) は `inner_iter` 行がない CSV に `{}` を返し、[同:127](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:127) はその場合に失敗を登録しません。初回残差がゼロの列も無条件に除外されます。`--nsub` は任意で、実際に反復数が倍になったかも確認しません。

   ファイルを作らない模擬入力で、次数だけを `2.0` にすると、**sub-iter 記録なし・nSub 比較なしで `VERDICT: PASS`, exit 0** を再現しました。`totalTime=NaN` でも PASS でした。

   既存 `run_0256`–`0259` が今回は `VERDICT: FAIL` になることは確認しました。しかし、それだけではこの欠落経路を検証できていません。

   **対案:** §6-6 の正式ゲートでは、全物理 step・必要全残差列・初回と最終 sub-iter の存在と有限性を必須にしてください。ゼロ始点から非ゼロになった列は除外せず拒否します。nSub 倍増比較を必須とし、実設定の BDF、FCT 作動、刻み、反復数も照合してください。

4. **Major — 実現可能性ゲートが solver と異なる物性・許容条件を検査している**

   **根拠:** [analyze_moment_order.py:137](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:137) は定数 `ρ_l=1000` を使いますが、solver の射影は [condensationRealizability_d.cuh:235](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:235) で EOS 後の温度に対応する物性を使います。H2O の定義は [condensationProperties_d.cuh:215](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:215) です。

   `run_0278/res_6.h5` の `Q2²≤Q1Q3` 違反数は、ツール既定の定数密度なら **221 ノード**、solver の温度依存式なら **31 ノード**でした。現在の違反数を、そのまま射影の失敗数として扱えません。

   さらに、検査は種 `_0` の二不等式だけです。`Q0>0,g>0,Q1=Q2=0` は両式を通りますが、§4.7 が明示的に修復する特異不整合です。負値・他の凝縮種も網羅していません。

   **対案:** 全凝縮種について、確定保存量と EOS 後の温度から solver と同じ物性・無次元化・許容条件で判定してください。非負性、特異不整合、非有限値も検査対象に含めます。float 書き戻し後の値に対して判定し、残る違反を物性差や丸めと区別してください。

5. **Minor — 基点の式に `1/a` の脱落が残っている**

   **根拠:** [plan:198](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:198) は依然として  
   `q_B = q_L + M⁻¹(r_L − r_H)`  
   です。直前の正本式、および修正済みの [passiveFct_d.cuh:13](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:13) と一致しません。

   **対案:** `q_B = q_L + M⁻¹(r_L − r_H/a)` に統一してください。

重点の判断は次のとおりです。

- **(a) 制限解除:** FCT 有効時に `θ_b` とトレーサの sub-iter floor を切る判断を支持します。有限な中間逸脱を許し、最後の確定状態で制限する構成は妥当です。モーメントの floor を残すことも支持しますが、その補正収支を免除してはいけません。
- **(b) 射影の時機:** dual-time は EOS 更新後の step 末尾に一回、定常は毎 step という判断を支持します。ただし「射影回数が減った」ことは精度の証拠にならず、補正量と時間次数の両ゲートが必要です。
- **(c) 閉合:** `Δ総量 = −境界流束 + ソース履歴 + H_rem` は正しい定義です。ただし `H` を確定増分から作るため、閉合は主に集計の整合確認です。独立総量差との照合と `|H_rem|` の制限が不可欠で、閉合式に floor／ピンを再加算してはいけません。

**推奨は一つです。流束形 BDF2＋step 末尾 FCT を維持し、①収支の時点整合、②欠落・非有限・保存ゲート、③時間次数と実現可能性ゲートの順に修正してから、統合検証を再開してください。**

最新スモーク `case/44.vitiated_air_wt/run_0278_passiveG_smoke_order3/` は `check_passive_budget.py: VERDICT: FAIL`。初期総量の誤りを除いても、`rog_0` の lim 補正だけで `5.89e-6` と閾値を超えます。通常の `check_convergence.py` も **NOT CONVERGED** でしたが、これは非定常の sub-iter 判定とは分けて扱っています。[run 索引](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:727)。

Python 単体試験は `ALL PASS`。ファイル変更・新規 CFD 実行はしていません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
