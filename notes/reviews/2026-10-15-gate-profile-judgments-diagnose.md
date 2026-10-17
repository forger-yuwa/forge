# codex 諮問 (diagnose): gate-profile-judgments

- **brief**: [`notes/reviews/briefs/2026-10-15-gate-profile-judgments.md`](../../notes/reviews/briefs/2026-10-15-gate-profile-judgments.md)
- **plan**: [`plans/active/tooling-generator-defaults-and-launch-gate.md`](../../plans/active/tooling-generator-defaults-and-launch-gate.md)
- **date**: 2026-10-15
- **commit**: `ebe0d6b6` (feature/generator-defaults-launch-gate)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-gate-profile-judgments-diagnose.codex.log` (codex: ok)
- **結論**: 次は plan §4.3・§6.2 を「出所によらない明示値と理由つき逸脱」に統一し、省略を拒否する契約を forge 0 step の負例で固定する。
- **extra**: `notes/investigations/2026-10-14-launch-path-inventory.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（対象 HEAD: `ebe0d6b6`。Critical は確認なし）

| 案 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 1：生成器だけ `fixed` の省略を FAIL | **却下・Major** | 省略値は手書きでも実行バイナリに依存する。現状は作業ツリーのソースから既定を抽出し、省略を許可する（[recipe_settings.py:154](/home/sano/work/forge-gate/solver_density_cuda/tools/recipe_settings.py:154)、[同:475](/home/sano/work/forge-gate/solver_density_cuda/tools/recipe_settings.py:475)）。**出所によらず、通常のプロファイル宣言では省略を FAIL** にする。省略自体を試す場合は、既存の `key=~` と理由つき逸脱を使い、実効値未確定を記録する。`nodeInletCornerWall` を変換来歴で確認する部分は採用。 |
| 2：CFL を `advisory` に移し、自動で逸脱記録 | **却下・Major** | `advisory` は範囲内の差を記録せず、範囲外でも無宣言で通る（[recipe_settings.py:481](/home/sano/work/forge-gate/solver_density_cuda/tools/recipe_settings.py:481)）。メモリ上で SERN 本段の `cfl_pseudo` を 0.5→0.75、仮の範囲を `[0.1,1.0]` としたところ、無宣言で FAIL/WARN とも0件、`applied=[]`。逸脱を渡しても「使われていない」となった。**機種・段の CFL は `fixed` に残し、明示した問題 YAML の上書きだけを理由つき逸脱へ変換する**。パス・キーは出所として自動記録し、変更理由は1行を共用できればよい。 |
| 3：代表 YAML を fixture に固定 | **採用・Minor** | 現状は可変の case YAML を直接読む（[regen_profiles.py:405](/home/sano/work/forge-gate/solver_density_cuda/tools/recipes/regen_profiles.py:405)）。固定コピー・元パス・sha256、元との差を別検査で WARN、でよい。参照ファイルの解決先も維持する。ただし fixture の更新は、推奨値の変更承認とは分ける。 |
| 4(i)：継続チャンクの restart 修正、SERN の段終了ゲート | **採用・Major** | `euler_loop` は格子をコピーした直後に補間している（[euler_loop.py:344](/home/sano/work/forge-gate/design/forge_design/feedback/euler_loop.py:344)）。SERN は欠落・長さ不一致の量をスキップし得る（[runner_sern.py:749](/home/sano/work/forge-gate/design/forge_design/evaluate/runner_sern.py:749)）。既存の段終了ゲートと保存量コピーを共用する。層流→SST では、前段の必須量検査と次段で追加する乱流量の初期化を区別する。 |
| 4(ii)：再生成格子を座標一致だけで同一と認定 | **却下・Major** | [restart_field.py:74](/home/sano/work/forge-gate/solver_density_cuda/tools/restart_field.py:74) は接続を検査せず、座標欠落時も続行する。さらに float64 座標を float32 に丸めて比較する。実コードの照合部分をメモリ上で実行すると、元入力から座標を約 `1.0e−8` 変えた別入力も「同一メッシュ」と判定した。**CPG の同一形状継続は元の変換済み格子を再利用する**。再生成が必要なら、元入力と新入力の座標・接続・CV順序・境界所属・幾何の同一性を来歴に結び付けて確認する。格子が実際に異なる場合だけ cross-mesh とする。 |
| 5：`passiveFct` の誤報を本体で修正 | **採用・Minor** | 省略を既定1として扱い、定常でも不作動 WARN を出す（[check_solver_config.py:146](/home/sano/work/forge-gate/solver_density_cuda/tools/check_solver_config.py:146)、[同:179](/home/sano/work/forge-gate/solver_density_cuda/tools/check_solver_config.py:179)）。**省略された定常・RK の通常動作では警告しない**。明示的な有効化要求の不成立は警告を残す。現在の `recipe_settings` にある「省略なら一律削除」は撤去し、判定を本体に集約する。 |

#10 の前には、次の2点も決める必要がある。

- **Major：強制開始の条件が未達。** 棚卸しには WT NS v1 のプロファイル欠落、axis-Mach の `ramp`、SERN の適応段の未収録が残る（[棚卸し:60](/home/sano/work/forge-gate/notes/investigations/2026-10-14-launch-path-inventory.md:60)）。SERN の試験も `restart_by_index` を差し替えており、段間移植の安全性を検証していない（[run_launch_path_tests.py:203](/home/sano/work/forge-gate/design/tests/run_launch_path_tests.py:203)）。**各経路の宣言方法と段間の負例試験が揃うことを、強制開始条件にする。** 関門の実装着手と全経路への強制開始は分ける。
- **Major：「省略禁止」と「config のバイト同一」は両立しない。** plan は明示を要求する一方、#11 はバイト同一を要求する（[plan:158](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:158)、[同:244](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:244)）。**追加を許すキーと旧実効値を列挙し、その明示化だけを例外にする。** また、書き写しは出力済みキーだけを採るため、SERN の `scalarGradient`・`implicitRelax` など省略された依存値は `fixed` 検査の対象にもならない（[regen_profiles.py:391](/home/sano/work/forge-gate/solver_density_cuda/tools/recipes/regen_profiles.py:391)）。必要な明示キーの一覧を先に確定する。

結論: 次は plan §4.3・§6.2 を「出所によらない明示値と理由つき逸脱」に統一し、省略を拒否する契約を forge 0 step の負例で固定する。

第 1 仮説: 関門の主要な抜け穴は、宣言した実効設定の確認を、現在のソース既定との一致で代用していることである。確度: 高  
  根拠: `sern2d_finite_base_frozen_tp_sst/main` の `space.limiterScaled: 1` を省略しても、実測で FAIL/WARN とも0件。config を変えず、メモリ上の既定値表だけを1→0にすると同キーで FAIL になった。`recipe_settings.py:475` の分岐と一致する。  
  反証条件: 通常宣言では出所・ソース既定によらず省略が拒否され、理由つき省略だけが明示的な逸脱として扱われること。

第 2 仮説: 再生成格子への保存量コピーでは、座標照合を同一性の証明に使うと異なる格子を受理し得る。確度: 高。`restart_field.py:74` の照合部分で再現済み。

判別 A/B: 同じ SERN 本段のプロファイル・config を使い、A は `space.limiterScaled: 1` を明記、B はその1キーだけを削除する。forge 0 step、見る量は FAIL の対象キー・WARN・逸脱記録。**Aが通りBが同キーで拒否されれば省略禁止の契約が成立、両方通れば第1仮説を支持する。** 現 HEAD のメモリ上の検査結果は両方とも通過だった。

やらない方がよいこと: 検査を通すために CFL を `advisory` へ移すこと、手書きという出所だけで省略を許すこと、座標照合成功を根拠に再生成格子へ無条件で index コピーすること。

呼び出し側の前提への異議: 「生成器×代表問題の出力を写した」は確認できるが、「機種の生産レシピを検証した」とは別である。偽の forge の PASS は書き写し用の仮データ（`regen_profiles.py:63`）であり、収束の証拠ではない。また、格子生成の決定性と、今回の新旧格子の同一性は別の条件である。

不足情報: #12 の生産既定の採否には、機種ごとの根拠 run・同一設定区間・収束／準定常 VERDICT が必要。本回答はその数値既定の妥当性を認定していない。ファイル変更・forge 起動は行っていない。

証拠: 足りる
