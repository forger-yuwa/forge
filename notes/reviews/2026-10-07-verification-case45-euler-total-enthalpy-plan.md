# codex レビュー: verification-case45-euler-total-enthalpy (plan)

- **plan**: [`plans/active/verification-case45-euler-total-enthalpy.md`](../../plans/active/verification-case45-euler-total-enthalpy.md)
- **stage**: `plan`
- **date**: 2026-10-07
- **commit**: `08fcf648` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m2
- **focus**: 遅延の plan 段レビュー: E3 (mesh_euler の実装、commit 3a9bc9f8) と E4/E4V (出口較正のやり直し) は実施済み。方針 (Euler は壁に寄せない配点、NS は mesh のまま) と登録・判定の穴、§8 完了条件 (原因の記録・影響の見積もり) を満たすのに足りないものを指摘してほしい
- **extra**: `design/forge_design/evaluate/runner_axismach.py`, `procedures/nozzle-design-workflow.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
Euler と NS の格子設定を分離する方針は妥当です。  
ただし、出口較正の合格を参照場全体の適格性に広げる条件と、§8 の原因・影響評価を完了する手順が不足しています。

コード、関連 plan・仕様・検証手順、手元の評価 JSON を確認しました。格子設定の分離・旧入力の拒否・E4V の判定分岐について、書き込みを伴わない6項目の確認は通りました。ソルバ本体を変更しない今回の範囲では、node の case/45 を中心に検証する選択は妥当です。

ただし、手元の `e4_recal_eval.json` は評価器 `39722e1a…` による段1までの記録で、最終状態は「保留・段2未実施」です。`run_0164` の実データもありません。**E4V 合格の数値は plan の報告として扱い、今回独立に再確認した結果とは区別します。**

1. **Major — 出口較正の合格条件では、δ_E・C2 用の参照場の適格性を保証できない。**

   **根拠:** [plan:107](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:107) の前提は残差・全温・出口 `M_common` ですが、直後の109行では合格場を δ_E 抽出・C2 にも採用しています。一方、[deltastar.py:419](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:419) は Euler の **ρUx の空間分布**を読み、各断面へ補間して欠損を計算します。出口 Mach と全温が条件内でも、この分布の時間変動や参照変更への感度は未評価です。

   **対案:** 出口較正の採用と δ_E・C2 参照の採用を別のゲートにしてください。後者には、登録した時系列での質量流量・使用する断面の ρUx、または固定した NS 保存場に対して Euler のスナップショットだけを替えた δ_E の感度評価を追加します。許容幅は下流の δ_E・C2 の誤差予算から事前に決め、対応量の `check_quasisteady` の VERDICT を残してください。

2. **Major — 凍結初期線 `run_0062` への影響評価が、全温の正側最大値だけで止まっている。**

   **根拠:** [plan:150](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:150) は「正常」の根拠を最大超過 `+0.073 K` としています。しかし、これは負側偏差の上限を与えません。また、初期線の実装が読むのは [cfd_initial_line.py:215](/home/sano/work/forge-integ-1005/design/forge_design/feedback/cfd_initial_line.py:215) の `Ux・Uy・sonic` から得る Mach・流れ角です。全温の健全性だけでは、これらの入力への影響を見積もったことにはなりません。

   **対案:** 凍結スナップショットと既存の末尾系列について、全温の両側偏差、初期線上の Mach・角度、軸アンカー・質量流束の変動と判定をまとめてください。既存の `cfdpin_line_stability_run_0062.json` は再利用できます。凍結源を変更する必要はありませんが、§8 には「今回の異常の影響を評価した範囲」と「格子誤差など未評価の範囲」を分けて記録する必要があります。

3. **Major — §8 を閉じる作業が残作業表にない。**

   **根拠:** [plan:54](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:54) の未完了項目はレビュー2件だけですが、[§8:137](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:137) は発生段・原因候補・設計入力への影響の記録を要求しています。E1 では `roe・roY*` の B→C 差の大きさが未記録（148行）、E2 は判別不能です。E3/E4 の完了だけでは、この診断を閉じられません。

   **対案:** §5.1 に、①発生段と未確定原因の整理、②初期線・出口較正への影響表、③下流参照の採用条件、を追加してください。保存物があれば B→C の同一節点差を定量化し、なければ未確定として残します。較正値の変更量 `−3.7011748e−4` は記録できますが、それを全温異常だけが引き起こした出口 Mach 誤差とは扱わないでください。§8 は原因候補と未確定範囲での完了を認めており、根本原因の確定まで無制限に計算を続ける必要はありません。

4. **Minor — E3 の「完了」と登録したテスト条件が一致していない。**

   **根拠:** [plan:101](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:101) は `design/tests` の「FAIL 0」を要求しますが、[166行](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:166) は既存の3件が失敗したと記録しています。

   **対案:** 3件の試験名・変更前後の commit・失敗内容を対応づけ、「登録条件は未達、既存失敗として免除」と明記してください。新規失敗ゼロという回帰評価は妥当ですが、FAIL 0 達成とは別です。

5. **Minor — 現行スコープと採用済み設定の記述が古い。**

   **根拠:** [plan:24](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:24) は出口較正の調整を対象外とする一方、E4/E4V は調整と採用を実施済みです。[methods/design/overview.md:360](/home/sano/work/forge-integ-1005/methods/design/overview.md:360) も旧値 `+3.770e−4` を現行値として扱っています。

   **対案:** §2 をユーザ決定後の範囲に更新し、仕様・手順・索引には新値と適用範囲を同期してください。旧判断は履歴として保持します。result レビューには、E4V の修正後評価 JSON、量別 VERDICT、メッシュ品質記録とハッシュをそろえてください。

**推奨は、`mesh_euler` の分離を維持し、上記1→2→3の順で登録と証拠を補ってから下流への参照採用・plan 完了を確定することです。** E2 と E4 段1の「判別不能」、E4V の独立検証という区別、0.1 K の閾値は維持してください。

ファイルは変更していません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
