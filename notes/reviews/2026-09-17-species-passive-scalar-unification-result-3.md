# codex レビュー: species-passive-scalar-unification (result)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-17
- **commit**: `98b11dc1` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **focus**: 5 回目 (4 回目 NO-GO C0/M3/m1 を全採用; §5.1 #28・§6.1・§9 に反映)。M1: 半径 (Q3/Q0)^{1/3} の underflow を指数分離 (cond_moment_xy) で解消し solver は「半径が存在し得ない塵」を Q1=Q2=0 に修復、ゲートは判定不能セルを FAIL・許容をファイルの精度から作る (codex の反例と実セル run_0421 を確認、回帰試験追加)。M2: FCT の次数緩和は BDF2 の受動種だけに限定 (流れ・化学種 [1.7,2.3]、BDF1 [0.7,1.3])。M3: §6-8 の検証行列を FCT 有効条件で取得 — BDF2 x coupling 0/1 は次数・収支とも PASS (モーメント 1.740/1.924/1.875/1.381、流れ・化学種 1.93-1.95)、BDF1 x coupling 0/1 は漸近水準 run_0453-0456 (dt 1e-6/0.5e-6/0.25e-6) で PASS (モーメント 0.965-0.974)、FCT 無効の対照 run_0445-0448 は 1.987-2.013 で PASS。m1: 推奨設定 (定常/dual-time の implicitRelax)・残作業表・run 索引を整合。問い: 残る偽合格・因果の穴はあるか、この plan を status done にして accepted へ移してよいか。run の実測は case/44.vitiated_air_wt/README.md (run_0429-0456) と case/09.Taylor-Green/README.md を自分で読むこと。設計 §4 の再設計提案は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
最新の BDF1/BDF2 × coupling 0/1 の次数・収支 `PASS` は再現でき、前回 M2/M3 の対応は支持します。  
ただし、場ゲートの偽合格と、指数分離による修復の未完了箇所が残っています。

1. **Major — 丸め許容によって、明確に実現不能な float32 モーメントが合格する**

   根拠: [passive_gate_common.py:277](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:277)。`rog` の相対丸めから作る共通 `eps` が1を超えると、281行の下限制約 `y < x*x*(1-eps)` が事実上無効になります。

   メモリ上の HDF5 に次の **float32 保存値**を入れ、実際の `check_field()` に通しました。`ro=1, T=250 K` です。

   ```text
   rog_0  = 1.401298464324817e-45
   roQ0_0 = 1
   roQ1_0 = 3.4747253568784376e-17
   roQ2_0 = 4.8294863261832316e-36

   Q1²/(Q0·Q2) = 250.000010
   check_field() → (True, [])
   ```

   `Q1² ≤ Q0·Q2` は250倍破れています。この不等式には `rog` が含まれず、正常数である `Q0/Q1/Q2` の丸めでは説明できません。「ファイルの精度から許容を作る」だけでは偽合格を閉じられていません。

   **対案:** 不等式ごとに関係する保存量の丸め幅を評価する。特に第1不等式を `rog` 由来の許容で緩めないこと。保証できないセルは既定方針どおり `FAIL` とし、この反例を回帰試験に追加してください。

2. **Major — 指数分離が書き戻しまで貫徹されず、射影後も実現不能状態が残る**

   根拠: [condensationRealizability_d.cuh:281](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:281)、同ファイル294行。`x,y` は対数で作りますが、書き戻しには依然として `rr=cbrt(q3/q0)` を使用します。途中の除算が非正規化数へ丸められても、`rr>0` なら対数経路を使いません。

   当該演算を double で再評価した反例です。

   ```text
   q3=1e-310, q0=3e13, q1=q2=0
   q3/q0 → 4.9406564584124654e-324
   rr    → 1.7031839360032601e-108
   対数で求める半径 → 1.4938015821857023e-108

   単分散 (x,y)=(1,1) へ射影して書き戻した後:
   x=1.1401674468, y=1.2999818068
   ```

   修復直後に `x≤1` を約14%破ります。同じ書き戻しが125行・216行にもあります。既存試験は無次元値の射影までを検査しており、保存量への往復を検査していません。

   また、[同ファイル284行](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:284)の修復失敗分岐は `Q1/Q2` をゼロに変更して290行で戻るため、306行以降の成分別収支に変更量が入りません。

   **対案:** 保存量への書き戻しも指数分離した共通処理に統一し、全修復分岐を収支記録へ通す。試験は「射影→保存量へ書き戻し→再判定」と符号付き・絶対補正量まで確認してください。

3. **Minor — 現在の判定と旧記述が混在し、残作業の完了状態が一意でない**

   根拠と修正対象:

   - [case/44 README:746](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:746)の旧系列には「場ゲートとも全水準 PASS」が残りますが、`run_0421_fctorder5noevap_double_dt4e-6_nsub40` の再判定は **`VERDICT: FAIL`、判定不能6セル**でした。
   - [plan:295](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:295)の #20 は「再取得中」、#21 は「検証中」「流れはビット一致」を残しています。後者は [case/09 README:59](/home/sano/work/forge/case/09.Taylor-Green/README.md:59)で既に訂正されています。
   - [solver-settings.md:169](/home/sano/work/forge/procedures/solver-settings.md:169)の `implicitRelax 0.7` 推奨には定常限定がなく、修正済みの `recommended-settings.md` と適用範囲が揃っていません。

   **対案:** 旧結果を削除せず「旧ゲート判定・現在は棄却」と明示し、残作業表は最新状態と後継 run への参照に整理してください。

再検証した主要結果は次のとおりです。番号は [case/44 の run 一覧](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:744)に対応します。

| 条件 | `case/44.vitiated_air_wt/` 配下の run | 次数ゲート | 収支ゲート |
|---|---|---|---|
| BDF2・coupling 1 | `run_0429`–`0432` | PASS | 全4本 PASS |
| BDF2・coupling 0 | `run_0433`–`0436` | PASS | 全4本 PASS |
| BDF1・coupling 1 | `run_0439,0449,0453,0454` | PASS | 全4本 PASS |
| BDF1・coupling 0 | `run_0443,0451,0455,0456` | PASS | 全4本 PASS |
| BDF2・FCT 無効対照 | `run_0445`–`0448` | PASS | INDETERMINATE |

FCT 有効の BDF2 モーメント次数 **1.740/1.924/1.875/1.381**、BDF1 **0.965/0.972/0.974/0.971**、無効対照 **1.987–2.013**を再現しました。次数緩和が BDF2 の受動種だけに限定されたことも確認しました。

`case/09.Taylor-Green/run_0156/0157/0159/0164` の収支も再判定で `PASS`。定常回帰の `case/44.vitiated_air_wt/run_0222` は **`NOT CONVERGED` / 準定常量 `ALL STEADY`**であり、固定点一致の証拠には使えません。

**推奨は、`in_progress` を維持し、指摘1→2→3の順に修正してから再レビューすることです。** §4の再設計や検証行列の全面再取得は現時点では求めません。修正したゲートで正式系列を再判定し、修正カーネルの反例試験と影響する回帰を確認してください。

ファイル変更・新規CFD実行はしていません。今回の指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 1
