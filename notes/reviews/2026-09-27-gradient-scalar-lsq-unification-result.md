# codex レビュー: gradient-scalar-lsq-unification (result)

- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **stage**: `result` (diff base `1cfd5623`)
- **date**: 2026-09-27
- **commit**: `03934178` (feature/sern-design)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m1
- **focus**: Phase 2 (node の既定を lsq に切り替え) の result。前提 #2g はユーザ決定『B』で警告 + 運用ルールに差し替え済み (§5.1 #6)。Phase 1 は合格記録済み (§6.1)。既定切り替えの実装 (24365d72)・設計 DB ゲート・確認 (DEFSW_*.txt) で plan を done / accepted にしてよいかを判断してほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**  
node の既定切り替えと LSQ 中核実装は支持します。ユーザ決定「B」と Phase 1 の合格判断を覆す根拠はありません。  
ただし、設計 DB の選別漏れ 2 件と文書の不整合を修正してから `accepted` に移してください。

1. **Major — 最終起動の実効値が不明でも、過去の `lsq` を引き継いで学習に採用する。**

   根拠: [runner_sern.py:44](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:44) の `_last_launch_value` は、キー欠落・JSON 解析失敗で `continue` するため、「最後の起動」ではなく「最後にキーを読めた起動」を返します。

   実関数への読み取り専用試験で、起動記録を `lsq → キー欠落`、`lsq → JSON破損` とすると、**どちらも実効値 `lsq`、`_learnable=True`** になりました。明示 `gg` なら拒否されます。実キャンペーンでの発生は未確認ですが、不明値を拒否する契約は破れています。

   **対案:** 最後の非空レコードだけを評価し、欠落・破損・未知値は `None` として拒否する。起動記録→実効値転記→学習採否まで試験する。現行の [5 例の試験:38](/home/sano/work/forge-sern-design/design/tests/run_sern_scalar_gradient_gate_tests.py:38) は完成済みの行を `_XF` に渡すだけなので、この欠陥を検出できません。後回しにした起動 ID 対応の実装までは不要です。

2. **Major — 学習から除外した `gg`・旧方針の評価が、最終 Pareto に混入する。**

   根拠: `_XF` は `_learnable` を使いますが、[driver_sern.py:289](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:289) の `summary` は `status == PASS` だけで選別します。さらに出力時に `flag_policy`・`scalar_gradient_effective` を落とすため、混在が成果物から判別できません。

   同一長さで `C_T_w=0.95` の現行 `lsq`、`0.99` の現行 `gg`、`1.00` の旧方針という合成 3 行を実関数へ渡すと、**学習は 1 行なのに、要約は `n_pass=3`、Pareto は旧方針の行だけ**になりました。既存の要約経路に残る統合漏れです。

   **対案:** Pareto・HV の母集団にも `_learnable` を適用し、出力に方針と実効値を残す。`rejudge` の行にも必要な情報を転記し、学習対象と最終選択対象が一致する回帰試験を追加する。

3. **Minor — 現行仕様と完了条件に、差し替え前の記述が残っている。**

   根拠:
   - [plan §4.4:67](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:67)、§5、§5.1 の #2g・#5j・#6、S4 は、なお「#2g 完了が既定化の前提」と記載。
   - [plans/README.md:27](/home/sano/work/forge-sern-design/plans/README.md:27) も「Phase 2 は #2g 待ち」。
   - [methods/gradient.md:172](/home/sano/work/forge-sern-design/methods/gradient.md:172) はスカラーを GG として説明し、「LSQ 化で解消する見込み」「既定 `gg`」が残存。

   **対案:** 旧判断は履歴として残し、現行条件を「B による警告＋運用ルール」に統一する。境界寄与の説明は既定 LSQ と明示 GG を区別する。§5.1 は未完了・移管済み・差し替え済みを明示し、今回の指摘の処置を記録してから閉じる。`methods/index.md` のリンクには問題を認めませんでした。

確認結果として、指定の `git diff 1cfd5623...HEAD` を取得し、関連する現行カーネルも確認しました。共有 `cInt`、差分形、境界 incidence の skip、並進周期だけの合併、スカラーの二重 gather 回避は設計に沿っています。Phase 2 による符号・単位・除算ガードの新たな欠陥は見つかりませんでした。

回収済み残差を現行 checker の判定関数で再計算すると、5 組とも通常判定は **`NOT CONVERGED`**、NaN/Inf・RISING はなく、変更後の床比判定は **`PASS`** でした。最大比は `case/48.flat_plate_cooled_m4/run_0953_sglsq_s2_lsq_ext/` の `rms_roe`、対 `run_0952_sglsq_s2_gg_ext/` で **1.418695**。記録と整合しますが、通常基準で「収束した」という意味ではありません。[run 索引](/home/sano/work/forge-sern-design/case/48.flat_plate_cooled_m4/README.md:55)

case/39 と凝縮双子の回収系列では `check_quasisteady.py` の **`OVERALL: ALL STEADY`** を再現しました。`DEFSW` は両ケースとも **E/P=A・N=PASS**、S3 記録は **0.9974／1.0184** で上限内です。GPU 計算と AWS 上の HDF5 からの量抽出は再実行していません。

**推奨は、既定 `lsq` を維持し、移動前に ①実効値取得の誤採用を修正、②Pareto まで選別を統一、③文書・残作業表を同期することです。** 今回の指摘の解消には Python の回帰試験で足り、S2/S3 の再計算は不要です。ファイルは変更しておらず、指摘は plan 未反映です。

指摘数: Critical 0 / Major 2 / Minor 1
