# codex レビュー: species-passive-scalar-unification (result)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-17
- **commit**: `b89da23f` (feature/sern-design)
- **codex**: effort `high`, 5.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: 4 回目 (3 回目 NO-GO C0/M5/m1 を全採用して再検証; §5.1 #27 と §9 に反映)。M1: convMethod 1 で全系列を取り直し (case/44 run_0384-0428) — sub-iter 依存の原因は implicitRelax 単独が 2 次でも成立 (run_0399-0404: S1 1.40e-5 / S3+FCT 1.32e-5 / +relax 6.32e-4)。M2: 空間スキーム整合 IC (run_0394) + 0.64 ms + dt 4e-6/2e-6/1e-6 の正式系列 run_0417-0420 は収支・場ゲート全水準 PASS、流れ・化学種 1.93-1.95、g 1.74 / Q2 1.92 / Q1 1.88、roQ0 1.38 の原因は FCT の非線形リミッタと確定 (対照 run_0425-0428 [passiveFct 0] は全モーメント 1.99-2.01 で PASS; 射影・蒸発・丸め・緩和は個別に排除) → 次数ゲートは「FCT 作動時は下限 1.3、スキームの 2 次は FCT 無効の対照で示す」に。M3: 射影と物理 step 末尾の後処理 (EOS 再更新含む) を passiveScalarScheme 1 限定にし旧経路の無影響を回復 (run_0382/0383/0393 で反復ノイズと同等)。M4/M5: 独立照合は固定 tol・--double は許容なし、場ゲートは g>0 かつ Q0>0 の全セルを判定 (許容はセルごとの丸め)、反例 2 件を回帰試験に。副産物: 射影が極端な (x,y) で候補距離 overflow により違反を残す欠陥を修正。問い: 残る偽合格・因果の穴はあるか、done にして accepted へ移してよいか。run の実測は case/44.vitiated_air_wt/README.md と case/09.Taylor-Green/README.md を自分で読むこと。設計 §4 の再設計提案は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
正式系列の次数・収支 PASS と、`implicitRelax` による sub-iter 差の増大は再現しました。  
ただし、場・次数ゲートに偽合格が残り、必須の検証行列も未完です。現状で `accepted` へは移せません。

1. **Major — 場ゲートが倍精度の実現不能状態を合格させる**

   根拠: [passive_gate_common.py:263](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:263) は半径 `rr` が計算できない wet セルを検査対象から外し、同ファイル:265–274 は倍精度データにも float32 の丸め許容を適用します。`--double` でも変わりません。

   ファイルを変更せず、メモリ上で `rog=1e-60, roQ0=1, (x,y)=(2,0.1)` を与えたところ、明白な不等式違反を `check_field=(True, [])` と判定しました。

   **実 run にも見逃しがあります。**  
   `case/44.vitiated_air_wt/run_0421_fctorder5noevap_double_dt4e-6_nsub40/res_160.h5` の配列 index `12430` は、指数を分離して計算すると `x≈1.21e20, y≈6.09e37`。しかし `check_passive_field.py --double` は **`VERDICT: PASS`** です。solver 側も [condensationRealizability_d.cuh:242](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:242) の `q3` の underflow → `rr=0` → return で修復を飛ばします。

   **対案:** 保存精度に対応した許容を使い、`wet && !okr` と非有限な `x,y` を不合格にする。solver・ゲートとも指数分離などで半径を安全に計算し、上記反例と実セルを回帰試験に追加する。なお、正式系列 `run_0417–0420` は独立した厳しい再検査でも違反ゼロでした。

2. **Major — FCT 用の次数緩和が流れ・化学種まで緩め、BDF1 判定も壊す**

   根拠: [analyze_moment_order.py:59](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:59) で決めた下限 `1.3` を、同ファイル:130 で**全評価量**に適用しています。

   正式系列の読込結果をメモリ上だけで差し替え、`ro` の観測次数を **1.400**、sub-iter 差をゼロにした反例が **`VERDICT: PASS`, exit 0** になりました。FCT によるモーメントの次数低下という根拠では、流れの低次数化まで許容できません。

   また、`--bdf 1 --expect-fct` は許容区間が **`[1.3,1.3]`** になります。実際に `run_0325–0328` の再判定で、`roQ0=0.869` などをこの区間で評価していました。

   **対案:** 閾値を BDF 次数・評価成分ごとに分ける。緩和対象の受動種だけに BDF2 下限 `1.3` を適用し、流れ・化学種は `[1.7,2.3]`、BDF1 は `[0.7,1.3]` を維持する。両反例を CLI 終了コードまで検査する。

3. **Major — 完了条件の BDF1/BDF2 × coupling 0/1 が揃っていない**

   根拠: [plan §6-8:342](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:342) は、FCT 有効条件で **BDF1/BDF2 × `passiveImplicitCoupling` 0/1、それぞれ3水準＋nSub倍増**を要求しています。

   再取得した `run_0417–0420` は BDF2・coupling 1 のみです。case/44 の config を走査しても、SLAU・`convMethod≥1`・SFR≥2・FCT 有効の **coupling 0 run はありません**。BDF1 の旧系列 `run_0325–0328` は、ゲートの区間誤りを除いても `ro` 次数 **3.328**、nSub 比 **0.546** で未達です。

   **対案:** 未完の3組を §5.1 に明記し、現在の滑らかな IC・修正済みゲートで取得する。BDF2・coupling 1 の成功を、検証行列全体の完了として扱わない。

4. **Minor — 現行の推奨・残作業・run 索引が最新判断と不整合**

   根拠:
   - [recommended-settings.md:128](/home/sano/work/forge/procedures/recommended-settings.md:128) は S3 に `implicitRelax 0.7` を一律推奨し、plan §5.1 #26 の「dual-time 時間精度試験では使わない」を反映していません。
   - [plan:303](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:303) の #8 は依然「未達」「再取得中」で、#27 の「済」と衝突します。
   - [case/44 README:740](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:740) は棄却された `run_0378–0381` を現在も「正式」「正本」と呼んでいます。

   **対案:** 定常と dual-time の推奨を分け、残作業表を今回の未完事項へ更新する。旧系列は履歴として残し、正式系列への置換を明記する。

再確認できた点も明確にします。[case/44 の run 索引](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:744) の正式系列は、次数 **g 1.740 / Q2 1.924 / Q1 1.875 / Q0 1.381**、収支4 run とも **`VERDICT: PASS`**。FCT 無効対照はモーメント **1.987–2.013、次数 `VERDICT: PASS`** でした。[case/09 の run 索引](/home/sano/work/forge/case/09.Taylor-Green/README.md) の `run_0156/0157/0159/0164` も収支 PASS を再現しています。旧経路 `run_0382/0383/0393` の差が反復ノイズと同程度という記録も裏付けられました。

`run_0417` の定常収束ツールは **`NOT CONVERGED (stalled/plateau)`** です。これは非定常の次数試験なので、それ自体を不合格理由にはしていません。

**推奨は、`in_progress` を維持し、1 → 2 → 3 → 4 の順に解消してから再レビューすることです。** §4 の再設計は不要です。ファイル変更は行っておらず、今回の指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
