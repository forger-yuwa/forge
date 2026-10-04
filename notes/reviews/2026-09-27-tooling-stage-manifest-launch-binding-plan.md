# codex レビュー: tooling-stage-manifest-launch-binding (plan)

- **plan**: [`plans/active/tooling-stage-manifest-launch-binding.md`](../../plans/active/tooling-stage-manifest-launch-binding.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `b3384298` (feature/sern-design)
- **codex**: effort `high`, 3.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
課題の同定は正しく、起動 ID による直接対応と SHA-256 による区間分離を支持します。  
ただし、履歴の生成・退避、成果物への ID 付与、DB の選別、収束判定への接続が不足しています。以下を実装前に計画へ反映してください。

1. **Major — `StageManifest.add` で起動 ID を読む設計は、既存の呼び出し順では成立しない**

   **根拠:** [case/57 の生成器:144](/home/sano/work/forge-sern-design/case/57.transition_flat_plate/gen_runs.py:144) と [case/49 の生成器:325](/home/sano/work/forge-sern-design/case/49.plate_annular_cavity_m5/gen_runs.py:325) は、**起動前**に `add`・`write` を実行します。この時点では対象起動の sidecar はありません。また、[main.cpp:277](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:277) は CSV を毎起動 `trunc` し、[main.cpp:2205](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2205) は step を 0 から始めます。「再試行は一つの履歴へ追記」という T4(c) は現行動作と異なります。

   履歴の退避も CSV 単独です（[runner_sern.py:774](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:774)、[r1_driver.sh:32](/home/sano/work/forge-sern-design/case/39.periodic_hills/r1_driver.sh:32)）。`StageManifest.add` だけを変更すると、正常な段でも ID を保存できません。

   **対案:** 段の予定登録と、起動後の履歴確定を分けてください。現行の CSV 上書きを維持し、CSV と sidecar を一組として退避・確定する共通処理を導入します。sidecar には履歴の世代・行範囲・内容照合情報を持たせ、旧 sidecar の残留、途中停止、同じ step 範囲での再試行を識別してください。生成器と退避処理を §5・§7 に追加し、別 run の履歴を参照する [r1_concat_manifest.py:32](/home/sano/work/forge-sern-design/case/39.periodic_hills/r1_concat_manifest.py:32) も試験対象に含めるべきです。

2. **Major — 「評価成果物の起動」を識別する記録先が設計にない**

   **根拠:** [plan:40](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:40) の ID 出力先は起動ログ・実行ログ・残差 sidecar だけですが、[plan:45](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:45) は `res` の生成元を要求しています。現行の場出力はローカル step のファイル名で上書きし（[output.cpp:97](/home/sano/work/forge-sern-design/solver_density_cuda/output/output.cpp:97)）、境界出力も同様です（[output.cpp:324](/home/sano/work/forge-sern-design/solver_density_cuda/output/output.cpp:324)）。

   さらに SERN 評価は単一の最終場だけを読みません。[sern_forces.py:119](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_forces.py:119) は境界 HDF5 を step 順に集めます。別設定の短い失敗起動が一部を上書きすると、旧起動の後半と新起動の前半が混在します。最後のファイルに正しい ID を付けても、定常性を判定した時系列全体の由来は保証できません。

   **対案:** `output/output.cpp` をスコープに追加し、場・境界 HDF5 の双方へ `launch_id` を属性として保存してください。評価に使った全ファイルを起動記録へ照合し、異なる hard キーの成果物や由来不明の成果物を混ぜた評価を拒否します。T7 は「失敗起動が出力前に落ちる場合」に加え、**一部の場・境界出力を上書きしてから落ちる場合**を必須にしてください。

3. **Major — `runner_sern.py` の列追加だけでは設計 DB と選別に反映されない**

   **根拠:** DB 行への転記は [driver_sern.py:163](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:163) の `_op_summary` が列を列挙して行います。`runner_sern.py` に追加した `launch_id`・SHA・`scalarGradient` は、このままでは転記されません。

   また、[driver_sern.py:190](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:190) は行の `flag_policy` に現在の定数を付け、[driver_sern.py:236](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:236) はその一致だけで学習対象を選びます。これは成果物の実効設定に基づく選別ではありません。§7 に同ファイルがないため、現在の実装ステップでは T7 を満たせません。

   **対案:** `driver_sern.py` の転記・台帳保存・再読込・`_XF` を変更対象へ追加してください。作動点ごとの成果物 provenance を保存し、必要な全作動点で SHA・実効値がキャンペーンの受入条件を満たす場合だけ選別します。`launch_id` は追跡用とし、起動ごとに異なる ID 自体を互換性条件にはしないこと。T7 は `metrics.json` だけでなく、**台帳を再読込した後の学習対象集合**まで検査すべきです。

4. **Major — 起動の対応が正しくても、収束判定側が欠落した履歴を飛ばして別区間を判定できる**

   **根拠:** [check_convergence.py:229](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:229) は履歴が欠けると読み飛ばし、区間が一段なら無条件に `residual_history.csv` を代用します。空の履歴も読み飛ばします。これらには、選択した段の起動 ID と実際に読む CSV の照合がありません。

   ファイル I/O をメモリへ置換して現行関数を実行したところ、同一区間の `s1→s2` で **終段 `s2.csv` を欠落させても、`s1` の2行だけから判定用 CSV が生成されました**。これは収束の実測ではなく、判定対象の欠落を受理する挙動の再現です。

   **対案:** `check_convergence.py` を明示的な変更対象にしてください。選択区間の全履歴について、存在・非空・sidecar の対応・対象行の被覆を確認し、不足時は理由付きの判定不能と非ゼロ終了にします。代替 CSV も ID の一致を必須にしてください。T4/T6 に CLI 経由の試験を追加し、欠落・`unknown`・SHA 不明から PASS が出ないことを検査します。「途中に `unknown` があるが、その後の終区間は完全に既知」の扱いも明文化が必要です。

5. **Minor — 却下済みの対応方式が実装指示に残っている**

   **根拠:** §4 は ID による直接対応ですが、[plan:26](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:26) は「設定ハッシュ＋起動順」、[plan:50](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:50) は「`_normalize` を単調対応に」と指示しています。[plans/README.md:28](/home/sano/work/forge-sern-design/plans/README.md:28) も旧案です。§7 の「起動記録に2項目追加」も、sidecar と成果物まで扱う変更規模を表していません。

   **対案:** §2・§5・§5.1・§7 と索引を §4 に統一してください。実装順は「記録形式と失敗時の扱い → ソルバ出力 → 段生成・退避 → 区間判定 → DB → 統合試験」とするのが妥当です。

6. **Minor — 既存テストの期待値と新仕様の整合が未整理**

   **根拠:** [chi テスト:36](/home/sano/work/forge-sern-design/solver_density_cuda/tools/test_stage_manifest_wall_normal_chi.py:36) は、起動 ID・SHA・`scalarGradient` の証拠がない旧 manifest を自動連結できることを要求します。これは §4.5 の `unknown` 方針と両立しません。「異なる `cfg_fnv` の fixture を直す」だけでは不十分です。

   また、T3 が完全一致を要求する `value_origin` は、現行コードでは cell の強制 GG が `explicit(cell->gg)`（[solverConfig.cpp:261](/home/sano/work/forge-sern-design/solver_density_cuda/input/solverConfig.cpp:261)）です。§4.4 の三値への変換規則がありません。

   **対案:** 旧テストを「証拠付き legacy」と「証拠なし unknown」に分け、新仕様に合わせて期待値を更新してください。由来は正規化した種別と元の理由文字列を別々に保存します。実バイナリ試験は使用ケース・段数・保存間隔を指定し、T3 の cell は設定解決の試験と位置付けてください。通常の計算回帰は [verification/README.md](/home/sano/work/forge-sern-design/procedures/verification/README.md) の node 方針に合わせるべきです。

目的・数値面については、既存の accepted plan が扱った YAML キーによる分離を補完する課題であり、解決済みの重複ではありません。同一 `cfg_fnv` の起動実効値 `[0,1]` が段では `[1,1]`・1区間になることも再現しました。SHA が異なれば区間を分ける方針は保守的ですが、この用途では妥当です。数値作用素を変更しないため、float32・block-DPLUR・周期・軸対称の数値検証を一式やり直す必要はありません。

**推奨は、起動 ID の直接対応案を維持し、指摘1〜4の経路を実装前に具体化することです。** 優先順は履歴の確定・退避、成果物への ID 保存、収束判定の欠落拒否、DB の保存・選別です。その後、指摘5〜6の文書と試験を同期してください。

既存の manifest テスト2本は実行して PASS でしたが、上記の欠陥は検出できていません。`test_gate_bad_input.py` は一時ファイルを作成するため今回は未実行です。ファイル変更・CFD run の新規実行はしておらず、レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
