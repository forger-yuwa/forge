# codex レビュー: tooling-convergence-and-wall-resolution-gates (plan)

- **plan**: [`plans/active/tooling-convergence-and-wall-resolution-gates.md`](../../plans/active/tooling-convergence-and-wall-resolution-gates.md)
- **stage**: `plan`
- **date**: 2026-10-15
- **commit**: `20f6347b` (feature/generator-defaults-launch-gate)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.6 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-convergence-and-wall-resolution-gates-plan-judge.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M6/m2
- **focus**: §4.6 (計算の後の判定を 1 本の道具 judge_run.py にまとめ run_case.sh・設計チェーン・Stop フックから自動で回す) だけを見てほしい。他の節は 2026-09-19 のレビュー済み。判定の対象 (特に準定常の量) の決め方、判定不能の扱い、既存の判定ツール (check_convergence --segment・check_quasisteady・check_wall_resolution・停止規則) との重複と責務、run_case.sh に足すことの影響 (関門・起動の後のエコーの照合と同じ枝)

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

§4.6 の自動化は必要ですが、現在の呼び出し順では本段の記録が確定する前に判定してしまいます。
既存ツールを束ねるだけでは欠損入力による誤合格も残ります。以下を計画に反映してから実装してください。
対象は §4.6 のみです。ファイルは変更していません。

1. **Major — 判定の実行時点が、段情報・停止記録の確定より早い。**

   **根拠:** [`runner_axismach.py:2414`](/home/sano/work/forge-gate/design/forge_design/evaluate/runner_axismach.py:2414) は `run_case.sh` の終了後に段履歴と manifest を記録します。本段の `_record` と `STOP_RULE_STATE.json` の作成も同様です（[2467行](/home/sano/work/forge-gate/design/forge_design/evaluate/runner_axismach.py:2467)）。したがってラッパー内で `--segment` を呼ぶと、前段までの manifest を判定する可能性があります。

   停止規則もラッパー終了を待ってから、停止位置より後の出力を削除し、最終記録を確定します（[`stop_rule.py:953`](/home/sano/work/forge-gate/design/forge_design/evaluate/stop_rule.py:953)、[998行](/home/sano/work/forge-gate/design/forge_design/evaluate/stop_rule.py:998)）。ラッパー内の全出力監査は、この整理前の出力を読んでしまいます。

   **対案:** 単段と段階起動を明示的に分ける。段階起動ではラッパーが段単位の実行結果・健全性を保存し、親が履歴保存・manifest 更新・停止記録確定・出力整理を終えてから、総合判定を一度だけ確定する。停止規則側の代替後処理 `postprocess()` も同時に移行対象へ追加してください。

2. **Major — 機種だけでは、判定すべき量・基準・用途を決められない。**

   **根拠:** 停止規則が記録する M は「出口コア M」であり、計画の「出口の軸の M」と異なります。また `REACHED` は設計量の定常性を保証しないと明記されています（[`stop_rule.py:35`](/home/sano/work/forge-gate/design/forge_design/evaluate/stop_rule.py:35)、[46行](/home/sano/work/forge-gate/design/forge_design/evaluate/stop_rule.py:46)）。

   既存の設計ゲートには、質量収支の**絶対的な平均・振幅条件**と、δ_E・δ_target の**連続する二窓の条件**があります（[`deltastar_loop.py:669`](/home/sano/work/forge-gate/design/forge_design/feedback/deltastar_loop.py:669)、[698行](/home/sano/work/forge-gate/design/forge_design/feedback/deltastar_loop.py:698)）。一定の質量収支誤差は準定常でも保存則を満たしません。SERN はさらに `C_T_with_shear`・`C_L`・`C_M` を判定します（[`sern_gates.py:128`](/home/sano/work/forge-gate/design/forge_design/metrics/sern_gates.py:128)）。

   **対案:** 機種の既定を基に、**run ごとの判定仕様を実行前に確定・保存**する。量の定義、単位、抽出器、対象領域、区間、窓、閾値、必須性を含める。停止規則の系列は定義・出典・区間が一致するときだけ再利用する。既存の設計受理条件は各設計ゲートに残し、`judge_run.py` は共通検査の結果を提供する責務に限定してください。

   同じ仕様に、継続 run の `--from-floor` と非定常解析の適用範囲も必要です。現案の「manifest がなければ全体」だけでは、既存の床参照運用を自動化できません（[`calculation-workflow.md:367`](/home/sano/work/forge-gate/procedures/calculation-workflow.md:367)）。

3. **Major — 既存ツールの正常終了を、そのまま完全な検査の証拠にできない。**

   **根拠と再現:** ファイル I/O をメモリ上に置き換え、既存関数を実行しました。

   - 同一区間の二段のうち本段履歴を欠損させると、前段だけを連結して `convergence_ok=True`。欠損段を `continue` する実装です（[`check_convergence.py:234`](/home/sano/work/forge-gate/solver_density_cuda/tools/check_convergence.py:234)）。
   - `asym` だけを要求し、その抽出器が適用不能だと、評価量ゼロでも `STEADY`・severity `0`（[`check_quasisteady.py:417`](/home/sano/work/forge-gate/solver_density_cuda/tools/check_quasisteady.py:417)、[507行](/home/sano/work/forge-gate/solver_density_cuda/tools/check_quasisteady.py:507)）。
   - `step=[1,1,1,1]`、値が一定の系列も `STEADY`。時間方向の進展を確認できていません。

   壁ツールにも、指定 step がなければ別 step を使う経路があります（[`check_wall_resolution.py:327`](/home/sano/work/forge-gate/solver_density_cuda/tools/check_wall_resolution.py:327)）。計画の事前照合は必須です。

   **対案:** §5.1 #12 の前提に、**要求した対象と実際に評価した対象の完全一致**を追加する。段・列・量・壁・時刻の欠損、空集合、時刻の重複・逆行、必要窓の不足は `INDETERMINATE`。共通ツール側で修正し、集約側にも対象数と被覆の照合を持たせてください。manifest 不在時の単段扱いも、単段と確認できる場合に限定すべきです。

4. **Major — `viscMethod ≠ 0` を壁解像検査の条件にするのは誤り。**

   **根拠:** `viscMethod=0` は定数粘性です。ソルバは `vis_lam=visc_lam` を設定しています（[`gasProperties_d.cu:91`](/home/sano/work/forge-gate/solver_density_cuda/cuda_forge/gasProperties_d.cu:91)）。

   実在する node run、[`case/26.flat_plate_sst/run_0029_node_outletic_subsonic/solverConfig.yaml:16`](/home/sano/work/forge-gate/case/26.flat_plate_sst/run_0029_node_outletic_subsonic/solverConfig.yaml:16) は `viscMethod: 0`、`visc: 1.8e-5 Pa·s` で、[`bcondConfig.yaml:32`](/home/sano/work/forge-gate/case/26.flat_plate_sst/run_0029_node_outletic_subsonic/bcondConfig.yaml:32) に no-slip 壁があります。現案では検査対象から落ちます。

   **対案:** 実効粘性と BC から適用条件を決め、定数粘性も含める。必要な壁出力・`vis_lam` 等の有無を投入前に確認する。この run は壁の `outputHDFflg: 0` なので、終了後だけの確認では再計算が必要になる点も計画へ反映してください。

5. **Major — 総合状態と終了コードの対応が未定義。**

   **根拠:** 停止規則の実状態は `REACHED / CENSORED / DIVERGED / EXEC_ERROR / DATA_ERROR / UNAVAILABLE / NOT_USED` です。計画の `NOT_REACHED` では失敗と判定不能を区別できません（[`stop_rule.py:59`](/home/sano/work/forge-gate/design/forge_design/evaluate/stop_rule.py:59)）。

   また、収束・準定常 CLI は異なる失敗理由を同じ非ゼロ終了へまとめます（[`check_convergence.py:424`](/home/sano/work/forge-gate/solver_density_cuda/tools/check_convergence.py:424)、[`check_quasisteady.py:562`](/home/sano/work/forge-gate/solver_density_cuda/tools/check_quasisteady.py:562)）。`run_case.sh` は現在、solver の終了コードを基本とし、strict 指定時だけ収束判定を伝播します（[`run_case.sh:73`](/home/sano/work/forge-gate/solver_density_cuda/tools/run_case.sh:73)）。

   **対案:** 原判定を保持した構造化結果と、総合状態の真理値表を先に定義する。既知の `FAIL` は判定不能と併存しても隠さず、適用外は理由付き `NOT_APPLICABLE` とする。solver 異常終了、正常な監視停止、判定器例外を区別する。`CONVERGENCE_VERDICT.txt` には収束サブ判定を残し、既存 strict 変数の意味を無断で総合判定へ拡張しないでください。起動関門・`LAUNCH_ECHO`・残差プロットも維持する必要があります。

6. **Major — Stop フックが読む判定の鮮度と、何を block するかが不足している。**

   **根拠:** 現フックはセッションの claim を尊重し、残差 CSV より判定ファイルが古い場合に block する仕組みです。合格そのものを要求してはいません（[`hook_convergence_gate.py:44`](/home/sano/work/forge-gate/solver_density_cuda/tools/hook_convergence_gate.py:44)）。

   総合判定には config、manifest、準定常系列、壁出力、停止記録、`judge.yaml` も影響します。計画にある「道具の sha」だけでは、これらの変更後に古い PASS を再利用することを防げません。

   **対案:** 判定へ実行識別子・確定状態・入力の識別情報を保存し、原子的に公開する。フックは同じセッション範囲で「判定済みか・現入力に対応するか」を確認する。FAIL や判定不能を正しく報告してターンを終えることは許し、PASS になるまで終了を阻止する設計にはしないでください。

7. **Minor — 全スナップショットの再監査は現行ルールと衝突し、性能条件も曖昧。**

   **根拠:** [`AGENTS.md:74`](/home/sano/work/forge-gate/AGENTS.md:74) は、通常の終了時監査を「残差全履歴＋最終場」とし、途中の全場を開く必要はないと明記しています。§4.6 は全場・全境界を毎回読み直します。停止規則にも既に体積場の非有限監査があります（[`stop_rule.py:908`](/home/sano/work/forge-gate/design/forge_design/evaluate/stop_rule.py:908)）。

   **対案:** 通常監査は現行ルールに合わせ、完全監査を明示オプションにする。準定常抽出で読む場は同時に監査し、対応する既存監査記録も再利用する。「2分以内」は出力量・保存媒体・キャッシュ条件を固定した性能試験として登録してください。

8. **Minor — §4.6 用の統合試験と影響範囲が不足。**

   **根拠:** [§6.2](/home/sano/work/forge-gate/plans/active/tooling-convergence-and-wall-resolution-gates.md:181) は旧ツールの単体条件が中心で、[§7](/home/sano/work/forge-gate/plans/active/tooling-convergence-and-wall-resolution-gates.md:190) は依然「tools の3ファイル」です。

   **対案:** 上記の欠損反例に加え、段階起動、正常な監視停止、solver 異常終了、判定器例外、古い判定、定数粘性壁、他セッションの除外を統合試験へ追加する。起動関門失敗時の未起動、エコー照合・プロットの維持、終了コードの互換性も確認する。変更対象には runner・watcher・設計ゲート・フックを列挙してください。

**推奨は一つです。** `judge_run.py` を、実行前に確定した判定仕様と実行後の確定記録を使う、既存判定ツールの共通集約層として実装してください。実装前の修正順は **①確定タイミング → ②判定仕様と責務 → ③入力完全性 → ④壁の適用条件 → ⑤状態・終了コード・鮮度 → ⑥統合試験**です。

実データ確認では、保存済み `run_0528_liminlet_cfl1/wall_series.csv` の20点・step 2400〜48000を既存ツールで再判定し、指定4量の `OVERALL: ALL STEADY` を再現しました。一方、同 run の保存済み収束判定は `NOT CONVERGED` です（[準定常記録](/home/sano/work/forge-gate/notes/investigations/2026-10-03-limiter-inlet/evidence/run_0528_liminlet_cfl1/QUASISTEADY_wall.txt:1)、[収束記録](/home/sano/work/forge-gate/notes/investigations/2026-10-03-limiter-inlet/evidence/run_0528_liminlet_cfl1/CONVERGENCE_VERDICT.txt:1)）。これらを独立した判定として保持する設計が必要です。未実装の総合判定・性能の合格は主張していません。

証拠: 足りる

指摘数: Critical 0 / Major 6 / Minor 2
