# codex 諮問 (diagnose): conjugate-benchmarks-closure

- **brief**: [`notes/reviews/briefs/2026-09-30-conjugate-benchmarks-closure.md`](../../notes/reviews/briefs/2026-09-30-conjugate-benchmarks-closure.md)
- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **date**: 2026-09-30
- **commit**: `70372d63` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 2.5 min, rc=0
- **結論**: **A/C を分割し、上記の限定を入れた A の部分結果を result レビューに提出し、C は原因未確定・未完了として active に残す。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | A の結論文：**限定して採用** | `case/64.conjugate_pipe_wall/run_0006_a1_r32`・`run_0007_a1_r64`・`run_0009_a2_r32` の `EVAL_CONJ.txt:9` は `VERDICT: PASS`。`run_0010_a2_r64/EVAL_CONJ.txt:8` は `Q_up/Q_tot` のみ判定不能。**FP64・node・forge の流れ場を固定した伝熱問題**への限定を明記する。float32、流れ場自体の正しさ、共役伝熱全般の検証へ拡張しない。 |
| **Major** | 「r16 の不合格は forge の粗格子誤差」：**A2 は要再検証** | 追加参照格子の A/B は **A1 だけ**。A2 の `run_0011_a2_r16_df20/EVAL_CONJ.txt:4` は差 **1.7807 %**、U **0.54178 %**、許容 **2 %**。差単独は許容内で、U を加えて FAIL になっている。対案は「A1 は粗格子誤差が主因と強く示唆される。A2 も不合格項目があるが原因の確定は保留」。 |
| **Major** | C を保留したまま全体を閉じる：**却下**。A/C 分割：**採用** | plan §1 は両ケースの成立を完了条件とし、§5.1 #6b（`plans/active/boundary-cht-conjugate-benchmarks.md:132`）も C の前提ゲート失敗を明記する。**C を未完了の後継 plan に移し、元の登録条件・失敗結果・完了条件を引き継ぐ**。A は限定された結果として result レビューに進める。分割を事前登録の総合 PASS に読み替えない。 |
| **Minor** | 「全節点 PASS」「格子倍増で約 1/4」：**記述修正を採用** | C1 n64 は `EVAL_CONJ.txt:1` の壁節点 **161** に対し、`SERIES_CONJ.txt:1` の対象は窓内 **77**。「評価窓内の全節点」と書く。また A1 の加熱区間温度差は **1.1905→0.29668→0.075254 %** だが、全長差は **1.2209→0.34192→0.13160 %**。約 1/4 は指標・領域を限定する。 |

結論: **A/C を分割し、上記の限定を入れた A の部分結果を result レビューに提出し、C は原因未確定・未完了として active に残す。**

第 1 仮説: A の粗格子差には forge の空間離散化誤差が支配的に寄与しているが、A2 r16 の FAIL まで同じ原因で確定する証拠は不足している。確度: **中**
  
  根拠: `case/64.conjugate_pipe_wall/ab_levels3/EVAL_CONJ_a1_r16_levels3.txt:3` では、A1 の参照を4水準にしても温度差 **1.2486 % + U 0.12003 % > 許容1 %**、`VERDICT: FAIL`。一方、A2 r16 は上表のとおり U が合否を左右する。A1 の試験で棄却できたのは「参照格子の解像不足だけで説明できる」である。
  
  反証条件: A2 の forge 保存場を変えず参照だけを細分化して熱流束が PASS になれば、「A2 r16 の現在の FAIL は forge 側の許容超過を示す」という解釈は撤回する。

第 2・第 3 仮説: C の node slip 欠陥起因説は**未確認**。端部への残差集中と窓内の小変動だけでは、境界実装・連成反復・抽出のどれが原因か区別できない。

判別 A/B: **A2 r16 の参照格子水準数だけを3→4にする**。対象は `case/64.conjugate_pipe_wall/run_0011_a2_r16_df20` の同じ step 600000 保存場。forge の追加計算は **0 step**、参照評価を1回行う。加熱区間・予熱域の熱流束について差と U を別々に比較する。  
→ **A：U が上限内で FAIL が残る**なら、「参照を1段細分すれば不合格が解消する」を棄却し、粗格子誤差説を補強する。  
→ **B：両領域が PASS**なら、元の FAIL を forge の許容超過の確証とする解釈を棄却する。判定不能なら原因帰属を保留する。この A/B は A の限定的な結果レビューを妨げる条件にはしない。

やらない方がよいこと: C の窓外を事後的にゲートから外して完了扱いにすること、slip を真因と決めて修正に入ること、A2 r64 の `Q_up/Q_tot` を小さい観測差だけで PASS とすること。

呼び出し側の前提への異議: **観測された項目別 PASS と、原因帰属・計画完了を分ける必要がある。** A の各 `SERIES_CONJ.txt:4` は `VERDICT: PASS`、入力系列は step 12010–600000 の50枚だが、それだけで流体収束や全前提ゲートの成立を独立確認したことにはならない。

不足情報: result レビューには、対象 run ごとの収束・G-if・G-cons・メッシュ品質の原出力と判定区間、`check_quasisteady.py` の実行条件・VERDICT、および §6 が要求する forge 側の **ΔT_s・Q_ax と不確かさの比較**を添付すること。今回の指定閲覧範囲ではこれらの原証拠と評価器コードを確認できず、完了承認までは出せない。**ファイル変更なし・plan 未反映。**
