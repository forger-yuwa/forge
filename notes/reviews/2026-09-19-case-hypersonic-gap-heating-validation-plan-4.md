# codex レビュー: case-hypersonic-gap-heating-validation (plan)

- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `0070d54d` (feature/sern-design)
- **codex**: effort `high`, 3.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m2
- **focus**: 3 巡目 (NO-GO, C0/M5/m1) を全件採用した改訂稿。M1 区間判定器の先行依存 (§4.7-1), M2 系列別温度比台帳と G13 の温度比1.00固定 (§4.4b), M3 T2-G0 の定量化と事前登録の順序 (§4.2b/§5/§5.1), M4 判定不能=不合格 (§4.7-3/G9), M5 離散収支の別 plan 依存 (§4.8/§5.1 #16), m1 抽出器の必須出力 (§4.8) が実質的に閉じたかを検証してほしい。まだ実装に進めない穴があれば指摘し、無ければ GO 判定とその条件を明示してほしい。原報 PDF は papers/gap_heating/
- **extra**: `notes/investigations/hypersonic-gap-cavity-heating-survey.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

T1/T2の構成と検証範囲は妥当で、3巡目の主要な設計問題は概ね解消しています。  
ただし、G10の壁解像判定には実測で確認できる誤合格経路があり、離散収支の依存関係にも記述の不整合が残ります。  
ファイルは変更していません。以下は **plan 未反映**です。

`plans/README.md`・`plans/accepted/` を確認した範囲で、本検証を代替する完了計画はありません。平面node・SLAU・block-DPLURの採用は現在仕様に整合します。共通ソルバを変更しない今回の範囲では、cell・周期・軸対称の追加回帰を要求する理由はありません。

3巡目への対応は次のように評価します。

| 前回指摘 | 今回の評価 |
|---|---|
| M1 区間判定器 | **計画上は解消**。三つの変更が同一キーになる欠陥を再現しましたが、段別runで単独判定する代替経路は成立します。修正前の`--segment`使用は禁止のままです。 |
| M2 系列別温度比 | **解消**。ローカル原報Fig.8(a)と三系列の値が一致し、G13の温度比1.00固定も適切です。 |
| M3 G0・事前登録 | **方針は解消**。等温基準、比較尺度、事前登録の順序が本文に入りました。ただし数値許容の確定は未完了で、T0/T2-G0投入の前提です。 |
| M4 判定不能＝不合格 | **計画上は解消**。G9も更新済みです。既存`check_quasisteady.py`だけでは残余ゲートを実装できないため、別途その判定実装が必要です。 |
| M5 離散収支 | **部分解消**。必要量と解除試験は明確になりましたが、依存先の実体と残作業表が未整合です。下記2。 |
| m1 必須出力 | **解消**。`thermCond`は登録済みで明示出力可能、`cp`の再計算と欠損時失敗も妥当です。[出力実装:46](/home/sano/work/forge/solver_density_cuda/output/output.cpp:46) |

1. **Major — G10の「局所 \(y_1^+\le1\)」を、指定ツールの既定PASSでは保証できない**

   **根拠:** [plan:193](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:193)は専用ツールによる局所判定を要求します。しかし、[判定実装:285](/home/sano/work/forge/solver_density_cuda/tools/check_wall_resolution.py:285)は「面積割合」を**点数割合**で計算し、評価不能点も最大10%まで許します。さらに既定では、目標超過を2%まで許してPASSにします。

   既存の`case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/`を再判定すると、

   - 最大 \(y_1^+=3.915\)、超過点割合0.4%。
   - 既定実行：`VERDICT: PASS (壁解像)`
   - `--target 1 --over-frac 0`：`VERDICT: FAIL`

   となりました。これは定常解の妥当性を主張する数値ではなく、**G10とツール判定の不一致の実証**です。非一様に細分する本計画では、点数割合を面積割合と扱うことも不適切です。

   **対案:** G10の呼出契約を、**全対象壁の出力確認・全評価点で有限値・`--target 1 --over-frac 0`**に固定してください。既存ツールが許す未評価点も、本計画側では不合格にします。面積割合の報告は実際の境界面積で重み付けし、ツール修正は既存の判定ツールplanへ依存登録する。これでソルバ改修なしに誤合格を防げます。

2. **Minor — 離散収支の対象範囲と、優先順の正本がまだ矛盾している**

   **根拠:** [plan:258](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:258)は診断planを「本計画の先行依存」としながら、直後には離散収支を「本計画では実施しない」としています。`plans/`全体を検索しても`tooling-energy-balance-diagnostics.md`は存在せず、[残作業#7:333](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:333)には依然「離散収支」が残っています。

   また「優先順」表では事前登録#12がT0/T1の#8より後です。本文の先行条件によって誤実行は禁止されていますが、正本として一本化されていません。

   **対案:** 診断planを実際に起票・リンクし、依存を**「離散収支評価の開始条件」**と明記してください。本計画の抽出器とG8は物理収支診断に限定し、#7から離散収支を外す。残作業表も「条件・幾何・分母 → 事前登録 → 抽出器検証 → T0/T2-G0 → 本体」の順に並べ直してください。現状のスコープなら、診断ソルバの完成を文献整理や物理熱流束抽出器の着手条件にする必要はありません。

3. **Minor — 文献基盤に旧結論が残り、確定した比較条件を再び取り違え得る**

   **根拠:** [survey:214](/home/sano/work/forge/notes/investigations/hypersonic-gap-cavity-heating-survey.md:214)はTH76を「同一位置の平滑模型実測」の分母に分類していますが、[同文書:92](/home/sano/work/forge/notes/investigations/hypersonic-gap-cavity-heating-survey.md:92)と改訂planは相関式を正本としています。また[同文書:107](/home/sano/work/forge/notes/investigations/hypersonic-gap-cavity-heating-survey.md:107)のcase/49深部への断定は、plan §4.11の「適用不確かさは未定量」より強い主張です。

   **対案:** surveyの要約・比較表も、**TH76の分母は相関式、温度比は系列別、case/49への傾向移植は未検証の仮説**に統一してください。

なお、case/48の上記runは本段単独で再実行した`check_convergence.py`が **`NOT CONVERGED (stalled/plateau)`**、`rms_roe`最終値 \(1.42\times10^{-1}\)、低下0.1桁でした。これを無条件の合格基準にしない方針は正しいです。[run索引](/home/sano/work/forge/case/48.flat_plate_cooled_m4/README.md)

**推奨は、T1/T2の構成を維持し、上記修正後に段階的に実装へ進むことです。** 実装前の優先順は、①G10の誤合格防止を仕様化、②依存関係と残作業表を一本化、③surveyを同期、です。その後、数値許容の事前登録と抽出器・残余判定の検証を終えてからT0/T2-G0を投入してください。計画の全面再設計は不要ですが、現状のままの無条件GOにはしません。

指摘数: Critical 0 / Major 1 / Minor 2
