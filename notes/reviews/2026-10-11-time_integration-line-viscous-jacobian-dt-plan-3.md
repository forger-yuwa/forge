# codex レビュー: time_integration-line-viscous-jacobian-dt (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-dt.md`](../../plans/active/time_integration-line-viscous-jacobian-dt.md)
- **stage**: `plan`
- **date**: 2026-10-11
- **commit**: `1cc9fd49` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **focus**: §6.5 (2026-10-11 に追加した、全部入り L7 の再実行 1 本の事前登録) と、台本 lvcmp_rerun.sh・判定 lvcmp_rerun_judge.py だけを点検する。§6.3 の lvcmp.sh・lvcmp_judge.py (plan-2 レビュー採用済み・実行済み) を流用している。再現の閾値 (4 量 ≥ +0.3 % かつ窓の中央値の床の比 ≥ 2、再現せず ≤ 0.1 % かつ < 1.5) と言える範囲が妥当か、判定器・台本に抜けや誤通過の穴が無いか。設計の諮問は diagnostician (notes/reviews/2026-10-11-lvcmp-rerun-design-diagnose.md)。
- **extra**: `case/45.isobutane_m6_d155/lvcmp_rerun.sh`, `case/45.isobutane_m6_d155/lvcmp_rerun_judge.py`, `case/45.isobutane_m6_d155/lvcmp_judge.py`, `notes/reviews/2026-10-11-lvcmp-rerun-design-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

L7 を 1 本追加する目的と、閾値を固定した探索的な比較は妥当です。ただし現状の判定器は正常入力でも必ず `INVALID` になり、その修正後にも誤通過の穴が残ります。「止まる状態の再現」という結論も、保存済みの準定常判定と整合しません。

確認範囲は §6.5、台本・判定器とその依存処理です。元の `run_0553`〜`0555` は手元にないため、実測の確認は保存済み `lvcmp_judge.json` までです。コードの反例はファイルを書き換えず、実関数をメモリ上へ抽出して確認しました。

1. **Major — `L7r` は入力ゲートを必ず通過できない。**

   **根拠:** [lvcmp_rerun_judge.py:24](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_rerun_judge.py:24) は `J.ARMS` だけを拡張します。一方、[lvcmp_judge.py:145](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_judge.py:145) は `ARM_CFG[name]`、同ファイル153行は固定辞書 `[name]` を参照します。どちらにも `L7r` がありません。両箇所で `KeyError: 'L7r'` を再現しました。例外はゲート不合格になるため、計算が正常でも判定不能です。

   **対案:** run の識別名と構成種別を分離し、`L7r` の期待設定・ログ検証には `L7` の構成を使ってください。`ARM_CFG` の追加だけでは153行の問題が残ります。投入前に、正常な L7r 入力が全ゲートを通る試験を必須にします。

2. **Major — 床の中央値が、登録した 2 万 step の窓を保証していない。**

   **根拠:** [lvcmp_rerun_judge.py:29](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_rerun_judge.py:29) の検査は「1000 行以上・有限」だけです。人工入力で以下をすべて受理しました。

   - `[90000,110000]` の指定に対して `90000〜90999` の1000行だけ。
   - `step=100000` を1000回重複。
   - 残差が全行 `−1`。

   欠損や重複によって中央値と床比を任意に偏らせられます。既存の収束判定テキストを読むだけの `segverdict()` は、この CSV の完全性を保証しません。

   **対案:** `step`・`phase` を列名で読み、出力仕様に対応する窓の全 step が重複なく存在すること、残差が有限・非負であること、分母の中央値が正であることを確認してください。窓・行数・CSV のハッシュを記録し、旧 L7 にも同じ中央値と閾値を適用してから「2 本で再現」とします。

3. **Major — 非有限な派生量が照合を通り、「再現」に入りうる。**

   **根拠:** [lvcmp_judge.py:215](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_judge.py:215) の相対誤差比較には有限性検査がありません。再計算値が有限、記録値が `Inf` の場合、比較は `Inf > Inf` となって偽です。実際の `assess()` に人工入力を与え、**有限なスナップショット再計算値に対して `Q_w=Inf` の記録が `REACHED` を返す**ことを確認しました。

   その後、[lvcmp_rerun_judge.py:80](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_rerun_judge.py:80) では `Inf >= 0.3` が成立します。他の3量を `+0.4 %`、床比を `2.1` とすると再現条件を満たします。

   **対案:** 記録値・再計算値・参照値・相対差・床比の全段階で `isfinite` を必須にしてください。分母の有効性も確認し、不正値は `INVALID` にします。JSON 出力も `allow_nan=False` にしてください。

4. **Major — 判定できるのは「停止規則で選んだ到達時の差」であり、「止まる状態の再現」ではない。**

   **根拠:** [plan:300](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:300) は「止まる状態のずれと高い床」の再現を主張します。しかし保存済みの [lvcmp_judge.json:201](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcmp_judge.json:201) では、

   - `run_0553_lvcmp_l7`: `NOT CONVERGED`、`Q_w: DRIFTING`。
   - `run_0554_lvcmp_l0`: `NOT CONVERGED`、`theta_r_70/94: TRANSIENT-UNSETTLED`。

   L7r の準定常判定も記録だけで、再現判定の条件ではありません。

   **閾値の評価:** `+0.3 % / ±0.1 %` と床比 `2 / 1.5` は、**探索的な運用上の判定幅としては採用可能**です。ただし `0.7 %・4〜5倍` という数値自体の再現や、統計的な同等性は示しません。L5 は別構成なので、L0 の再実行変動の代わりにもなりません。また床比 `<1.5` は `0.1` も含むため、判定器の「床とも生産と区別できない」は不適切です。

   **対案:** 結論を「同じ停止規則による到達時に、登録幅以上の正の差と高い残差中央値が再度観測された」に統一してください。「再現せず」は「今回、登録した差のパターンを再現しなかった」とします。準定常 VERDICT を併記し、「速度の利点が無い」も「今回の比較では登録幅以上の到達優位を検出しなかった」に限定します。

5. **Major — 初回保存前の発散を、登録した `DIVERGED` 分岐で扱えない。**

   **根拠:** [lvcmp_judge.py:191](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_judge.py:191) は系列が空なら即 `INVALID` を返し、223行以降の発散証拠確認に到達しません。一方、[m9_watch.py:194](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/m9_watch.py:194) は `res_nan_*` があれば初回2500 step出力前でも `DIVERGED` を記録します。`fail_step=100, rows=[]` の人工入力が `INVALID` になることを確認しました。

   **対案:** `DIVERGED` は発散証拠を先に確認し、発散前に存在すべき保存出力だけを照合してください。初回保存前なら空系列を許容し、`res_nan`・ログ・失敗 step の整合を要求します。

6. **Minor — 他セッション待ちが、同じ作業ツリーの別 run を見落とす。**

   **根拠:** [lvcmp_rerun.sh:119](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_rerun.sh:119) は `~/forge-faceh-audit/*` 配下の全 `forge` を除外します。同じツリーの別 case にある `forge` を模擬すると、`others() = 0` になりました。§6.5 の「他のセッションの forge が無いとき」と一致しません。

   **対案:** 起動前は自分の `forge` がまだ存在しないため、ツリー単位の除外をなくしてください。除外が必要なら、明示的に所有する PID に限定します。

7. **Minor — 約束した「全列の窓の中央値」が出力されない。**

   **根拠:** [plan:303](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:303) は全列を記録するとしていますが、[lvcmp_rerun_judge.py:104](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_rerun_judge.py:104) が保存するのは `rms_ro` の中央値だけです。区間付き収束判定は、この到達窓の全列中央値を代替しません。

   **対案:** 同じ窓で存在する保存量残差の全列を集計し、列別の中央値を JSON に残してください。主判定を `rms_ro` に限定する登録は維持できます。

**推奨は、上記を修正してから L7r を 1 本だけ実行することです。** 優先順は、①入力ゲートの成立、②床・派生量の不正入力拒否、③早期発散の分類、④結論の限定、⑤待機条件と記録の補完です。L0 再実行や point 解き直しは今回追加せず、正常入力・欠損・重複・非有限・早期発散の人工入力試験を投入条件にしてください。ファイルは変更しておらず、提案は plan 未反映です。

指摘数: Critical 0 / Major 5 / Minor 2
