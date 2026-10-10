# codex レビュー: time_integration-line-viscous-jacobian-faceh (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `6dee5bb4` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m0
- **focus**: §6.11 (2026-10-10 に追加した、値 3 のマスク 7/5 の比較の事前登録) と、台本 lvc75.sh・事前のゲート lvc75_pregate.py・本判定 lvc75_judge.py (作業ツリーの未 commit の新規ファイル) だけを点検する。§6.9 の lvcgeom.sh・lvcgeom_pregate.py・lvcgeom_judge.py (plan-5 レビューを採用済み、実行して PASS・棄却) を元にした。腕の違いは環境変数 FORGE_LVC_TERMS (7/5) だけ、バイナリは段 ③ の FP64 で共通。事前のゲートの V で、D と K の行 0〜3 の不変と、K の行 4 の変化が独立に計算した熱伝導の近傍 K の項と一致することを確かめる設計が十分か、分岐と言える範囲が妥当か、抜けている確認は何か。
- **extra**: `case/45.isobutane_m6_d155/lvc75.sh`, `case/45.isobutane_m6_d155/lvc75_pregate.py`, `case/45.isobutane_m6_d155/lvc75_judge.py`, `notes/reviews/2026-10-10-lvcgeom-result-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

同一バイナリでマスク7/5だけを変える方針と、結論を「この条件・2000 step以内の非有限化の回避」に限定する分岐は妥当です。  
ただし、事前ゲートには誤って `PASS` を返す経路があります。以下3点を本試験前に修正してください。

1. **Major — 非有限の参照値や残差が、事前ゲートを通過する。**

   根拠: [lvc75_pregate.py:176](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvc75_pregate.py:176) で残差の有限性を `I` に追加しますが、最終判定は `R` と `V` だけを集計します。同じ NaN が3本に入る人工入力で、**`I=False` にもかかわらず `PASS`** を確認しました。本判定では拒否されますが、事前停止には間に合いません。

   さらに、[同ファイル:245](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvc75_pregate.py:245) の `worst = max(worst, r)` は `r=NaN` を見逃します。比較コードに104面の人工入力を与え、1面の生の熱伝導率を NaN にすると、**解像面103・最大誤差0・V合格**になりました。この経路では不合格項目自体が残りません。

   **対案:** 使用する監査入力、κ・h・許容・誤差比をすべて有限性検査し、不正なら即 `INVALID`。最終判定でも全 `I` 項目の合格を必須にする。NaN・Infを注入する回帰試験を追加してください。

2. **Major — Kの行4の照合が、必要な面を網羅したことを保証していない。**

   根拠: [lvc75_pregate.py:219](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvc75_pregate.py:219) は監査記録の `fr[5]`・`fr[45]` に従って比較対象を選び、対象面の欠落・重複を検査しません。「解像面100以上」は網羅性の代わりになりません。

   Vの実コードによる人工入力試験では、**1面を対象外にして、その面のK差へ `1e12` を入れても合格**しました。また、同じ面の記録を重複させると、解像面数が104から208に増えて合格します。

   **対案:** `node_line` の順序と `line_prev/next` から期待する有向の接続集合を作り、各接続に面がちょうど1つあること、往復の面番号・向き・両端が整合することを先に検査する。拘束行も現在の壁フラグから独立に照合してください。既存の [lvcaudit_judge.py:173](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:173) の構造検査を再利用できます。そのうえで全対象の行4を比較し、100面条件は信号の解像条件として残すのが妥当です。

3. **Major — 事前ゲートの合格が、本試験の条件と保存された証拠に結び付いていない。**

   根拠: [lvc75_pregate.py:108](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvc75_pregate.py:108) は設定の一部だけを検査し、[lvc75_judge.py:284](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvc75_judge.py:284) は本試験4本どうしを比較します。**書き出しと本試験の全設定の比較がありません。** 例えば、書き出しの `turbulentPrandtl=0.85`、本試験4本が共通して `1.7` でも、この設定検査は通ります。事前に確認したκと本試験のκが異なる状態です。

   また、[同ファイル:308](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvc75_judge.py:308) が確認する証拠は保存JSONと `D.f64` の存在だけです。人工入力で、**K・状態・残差の書き出しが存在しなくても `intervention()` の全項目が合格**しました。`arms` モードも保存JSONの `PASS` だけで進みます。

   **対案:** 書き出しと本試験の設定差を、登録した step数・出力間隔・`extraFields` に限定して全項目照合する。事前ゲートのJSONには、使用した設定・監査記録・D/K・状態・残差のハッシュと判定器の版を保存し、`arms` 起動時と本判定時に照合してください。

**推奨:** 段③の同一FP64バイナリによる7/5比較を維持し、上記1→2→3を修正・§6.11へ反映してから実施してください。熱伝導項の式と符号は共通関数と整合しています。分岐の人工入力試験も登録どおりであり、今回の診断にcell・周期・速度評価を追加する必要はありません。

シェル・Python構文確認とメモリ上の人工入力試験を実施しました。AWSの実バイナリ・生データは再測定していません。ファイル変更・forge起動なし。**plan未反映**です。

指摘数: Critical 0 / Major 3 / Minor 0
