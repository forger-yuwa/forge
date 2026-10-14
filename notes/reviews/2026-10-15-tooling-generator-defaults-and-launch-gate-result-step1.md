# codex レビュー: tooling-generator-defaults-and-launch-gate (result)

- **plan**: [`plans/active/tooling-generator-defaults-and-launch-gate.md`](../../plans/active/tooling-generator-defaults-and-launch-gate.md)
- **stage**: `result` (diff base `origin/main`)
- **date**: 2026-10-15
- **commit**: `e9fd8811` (feature/generator-defaults-launch-gate)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 6.7 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-generator-defaults-and-launch-gate-result-step1.codex.log` (codex: ok)
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: 段取り 1 (宣言・関門を移行期の WARN で入れ、生成器の数値は変えない) を main に merge してよいか。実装 diff (origin/main..HEAD) と §5.1 の各項目の合格条件・試験の結果を突き合わせ、merge で他のセッションの run・台本が壊れる経路 (run_case.sh の関門、restart_field.py の厳格化、生成器の config へのキーの挿入、migrate) を洗ってほしい。段取り 2・3 (#12・#15) は未着手で、この review の対象外

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

現 HEAD は段取り 1 のマージ条件を満たしていません。case/45 の既存台本を止める回帰が残り、宣言と実際の計算条件の照合にも未実装部分があります。
段取り 2・3（#12・#15）は判定対象から除外しました。

1. **Major — case/45 の互換性修正が未完了で、明記されたマージ条件を満たさない**

   新しい生成器は B0 の `soft`・`mid` に `lineImplicit: 0`、`lineDtDirectional: 0`、`implicitThermalJacobian: 0` を挿入します。一方、[h370_run.py:559](/home/sano/work/forge-gate/case/45.isobutane_m6_d155/h370_run.py:559) と [h370_problems.py:507](/home/sano/work/forge-gate/case/45.isobutane_m6_d155/h370_problems.py:507) は、これらが **0 でも存在するだけで失敗**にします。保存済み fixture に現在の `finish_config` を適用し、この変化を再現しました。

   [plan #26:285](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:285) に別ブランチの修正記録がありますが、現 HEAD には未統合で、その修正にも `h370_problems.py` の `prep-check` が残っています。#11a の「merge は #26 の後」を満たしていません。

   **対案:** #26 を `prep-check` まで完了・統合し、旧 config／新 config の両方で登録設定の検査と `test_h370_launch.py` を通す。単なるキー存在判定を、明示化前後の実効値を確認する比較に置き換えてください。

2. **Major — プロファイルの `dimension`・`mesh` が照合されない**

   [recipe_settings.py:531](/home/sano/work/forge-gate/solver_density_cuda/tools/recipe_settings.py:531) の `derive_applies` は次元と格子種別を返しません。[同:677](/home/sano/work/forge-gate/solver_density_cuda/tools/recipe_settings.py:677) は取得できた属性だけを比較するため、この 2 条件は黙って無視されます。[launch_gate.py:445](/home/sano/work/forge-gate/solver_density_cuda/tools/launch_gate.py:445) にも、メッシュの実測条件とプロファイルを結び付ける照合がありません。

   SERN 2D の正常 config に対し、プロファイルの適用条件だけを `dimension: 3d`、さらに存在しない次元・格子種別へ変更しても、すべて **FAIL 0／WARN 0** でした。プロファイルに適用条件を書いても、実際には保証されません。

   **対案:** 変換成果物から得た次元・格子種別を宣言の `applies` と照合する。判定不能を一致扱いせず、2D／3D の取り違えと構造／非構造の取り違えを負例に追加してください。

3. **Major — 採用した起動後の実効値照合が未実装**

   [plan §4.3:162](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:162) は、未知のバイナリを拒否しない代わりに、起動エコーがある量を起動後に照合すると定めています。しかし [run_case.sh:66](/home/sano/work/forge-gate/solver_density_cuda/tools/run_case.sh:66) は `slauWallNormalChi`・`scalarGradient` の行を転記するだけで、要求値との比較はありません。

   #10 の記録にも未実装とあります。現在ソースの既定値を明示したことと、旧ビルド・FP64 ビルドがその値を実際に採用したことは別です。#26 に残る旧 FP64 バイナリの互換性未確認も、この穴では検出できません。

   **対案:** 起動エコーと要求値を比較し、不一致・確認不能を結果記録へ明示する処理を追加する。既定値が異なる偽バイナリで検出を確認し、この作業を独立した残件と合格条件として記載してください。

4. **Minor — 現在仕様の説明が実装・残作業表と矛盾する**

   [recommended-settings.md:55](/home/sano/work/forge-gate/procedures/recommended-settings.md:55) は「投入の関門はまだ無い」としていますが、実際には起動前に実行され、一部の失敗は移行期でも停止します。また [methods/design/overview.md:62](/home/sano/work/forge-gate/methods/design/overview.md:62) の「全段の config は自分のプロファイルに PASS」は、#27 に残る未結線経路・`ramp`・`warm_adapt` 等を含む記述として過大です。

   **対案:** 移行期の WARN／FAIL 条件を現実装に合わせ、検証済み範囲を「代表 11 機種・30 段」と明記する。#27 の残件と今回対象外の #12・#15 を保持し、段取り 1 の完了と計画全体の完了を分けてください。

独立確認では、`test_migrate_solver_config.py`、`test_check_solver_config.py`、レシピから生成した文書との一致確認は通りました。`case/` の移行 dry-run も **475 本中、移行対象 150 本・手動対応 2 本・失敗 0 本**で記録どおりです。保存済み代表 fixture に対する仕上げ処理では、**30 段・319 キーの挿入**とプロファイル照合を確認しました。ただし、これは全生成経路の再実行ではありません。

メッシュ品質の主張も、主作業ツリー `/home/sano/work/forge` の既存成果物を読み取り専用で再判定し、裏付けました。

| 成果物（主作業ツリーからの相対パス） | 再判定 |
| --- | --- |
| `case/45.isobutane_m6_d155/run_0139_rerun_pt08_noscale_cfl1_ext/nozzle.h5` | `VERDICT: FAIL`。緩和できない高 AR セル 85、最大 AR 約 1209.6、最大 φ 約 39.63° |
| `case/48.flat_plate_cooled_m4/run_0062_ypls_modefix/mesh.h5` | `VERDICT: FAIL`。壁層外の高 AR セル 82、最大 AR 約 2598.2 |

`test_launch_gate.py` は、一時ディレクトリを作れない読み取り専用環境のため開始できませんでした。統合試験一式を再実行したとは扱っていません。CFD の収束・準定常性について新たな合格判断はしていません。

**推奨:** 現 HEAD のマージと `accepted` 移動を保留し、上記の互換性回帰と照合漏れを修正してから、書き込み可能な検証環境で段取り 1 の統合試験を再実行してください。

証拠: 足りる

指摘数: Critical 0 / Major 3 / Minor 1
