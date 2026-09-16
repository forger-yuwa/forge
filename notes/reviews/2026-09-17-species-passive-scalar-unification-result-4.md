# codex レビュー: species-passive-scalar-unification (result)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-17
- **commit**: `de472488` (feature/sern-design)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m1
- **focus**: 6 回目 (5 回目 NO-GO C0/M2/m1 を全採用; §5.1 #29・§6.1・§9 に反映)。M1: 場ゲートの判定を保存量の 3 不等式 (Q1^3<=Q0^2 Q3 / Q1^2<=Q0 Q2 / Q2^2<=Q1 Q3) に分け、許容はその不等式に入る量の丸めだけから作る (rog の非正規化数が Q1^2<=Q0 Q2 を緩めない); 厳密に成り立つ不等式には許容を使わず、破れが許容を超えかつ許容が 1e-2 超なら判定不能 (FAIL)。codex の反例 (Q1^2/(Q0 Q2)=250) は FAIL を確認・回帰試験に追加。M2: 射影の書き戻しを cond_moment_writeback (q3/q0 が非正規化数なら対数経路) に統一し、修復不能分岐も収支記録へ; 往復 (射影→保存量→再判定) の単体試験を追加 (codex の q3=1e-310, q0=3e13 の例は x=1.14 → 1.000000)。m1: 旧系列の位置づけ・残作業表・solver-settings.md の implicitRelax (定常限定) を整合。修正後の binary で検証行列を 1 本再実行し旧結果との差が反復ノイズと同等 (run_0457/0458) であることを確認したので行列 (BDF2/BDF1 x coupling 0/1 + FCT 無効対照, run_0429-0456) は据え置き。問い: 残る偽合格・因果の穴はあるか、この plan を status done にして accepted へ移してよいか。設計 §4 の再設計提案は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
前回の場ゲート・書き戻しの反例は解消され、正式検証行列の PASS も再現しました。  
ただし、収支ゲートに別の偽合格が残ります。以下を修正してから `accepted` へ移すことを推奨します。

1. **Major — `|H_rem|` が補正収支の合計から漏れている**

   根拠: [check_passive_budget.py:148](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:148)。`base/pin/upper` は合計へ加えますが、`rem_rel` は別枠で比較するだけです。これは [plan §4.6:166](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:166) の「`|H_rem|` を含む総量比の合計 ≤1e-6」と異なります。

   実書式のログをメモリ上で生成し、実際の `parse_lines()` → `evaluate()` を通して再現しました。倍精度相当の固定許容でも、次が合格します。

   ```text
   floor絶対補正 / 総量 = 6e-7
   |H_rem|       / 総量 = 6e-7
   その他の補正        = 0
   閉合誤差・独立照合差 = 0

   要求される合計 = 1.2e-6 > 1e-6
   実際の判定     = True（表示 sum は 6e-7）
   ```

   **対案:** 倍精度では `rem_rel` も同じ合計判定へ入れる。float32 の step 比例許容は、その適用項目と合算方法を明文化し、他の補正まで緩めないこと。この反例を回帰試験へ追加してください。

   **今回の実測への影響:** 正式行列の FCT 有効16本について、漏れた項を含めて再集計した最大値は **2.951e-7**（`run_0433`、`roQ1_0`）でした。`run_0457/0458` も最大 **1.831e-7**。今回の行列がこの穴で誤合格した証拠はありません。

2. **Minor — 完了判断に使う記録と残作業の所在が未同期**

   根拠と修正対象:

   - [plan §9:406](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:406) に、`run_0222/0223` から「S3 固定点は cfl/relax 非依存」とする記述が残ります。両者の収束判定は **`NOT CONVERGED`**、準定常ツールの再実行は **`OVERALL: ALL STEADY`** です。支持できるのは準定常量の比較までです。
   - 最新修正の根拠となる `run_0457/0458` が、[case/44 の run 一覧](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:744) にありません。
   - [plan §10:420](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:420) の周期一様流発散は「F-cf12 候補」のままで、[followups の残作業表](/home/sano/work/forge/plans/active/condensation-followups.md:63)には登録されていません。

   **対案:** 固定点の旧主張を訂正済み履歴として明示し、最新2本を run 索引へ追加する。周期発散は後続 plan の残作業表へ登録し、本 plan §5.1 から引き継ぎ先をリンクしてください。

再検証結果は以下のとおりです。run 番号はすべて `case/44.vitiated_air_wt/` 配下で、主要成果物は各 run の最終 `res_*.h5`、`residual_history.csv`、`forge_run.log` です。

| 条件 | run | 次数 VERDICT | 収支 VERDICT |
|---|---|---|---|
| BDF2・coupling 1 | `0429–0432` | PASS | 全4本 PASS |
| BDF2・coupling 0 | `0433–0436` | PASS | 全4本 PASS |
| BDF1・coupling 1 | `0439,0449,0453,0454` | PASS | 全4本 PASS |
| BDF1・coupling 0 | `0443,0451,0455,0456` | PASS | 全4本 PASS |
| BDF2・FCT 無効対照 | `0445–0448` | PASS | INDETERMINATE |

FCT 有効 BDF2 のモーメント次数 **1.740/1.924/1.875/1.381**、無効対照 **1.987–2.013**を再現しました。前回の250倍違反は現行 `check_field()` で拒否され、書き戻しの反例も現行式では **x=y=1.000000** になります。`run_0430→0457` の密度最大差 **3.383e-7** は、再実行 `0457→0458` の **3.660e-7** と同程度でした。

**推奨は、指摘1のゲート修正・反例試験・正式系列の再判定、次に指摘2の記録同期を済ませてから `done/accepted` にすることです。**

ファイル変更・新規CFD実行はしていません。今回の指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 1
