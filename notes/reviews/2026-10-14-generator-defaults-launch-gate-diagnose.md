# codex 諮問 (diagnose): generator-defaults-launch-gate

- **brief**: [`notes/reviews/briefs/2026-10-14-generator-defaults-launch-gate.md`](../../notes/reviews/briefs/2026-10-14-generator-defaults-launch-gate.md)
- **plan**: [`plans/active/tooling-generator-defaults-and-launch-gate.md`](../../plans/active/tooling-generator-defaults-and-launch-gate.md)
- **date**: 2026-10-14
- **commit**: `f66113ee` (feature/generator-defaults-launch-gate)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **数値既定を一括変更する前に、plan §4.3 を「明示した profile・stage と、変換済みメッシュの実効条件に束縛した関門」へ修正し、CFD なしの負例試験を先行させる。**
- **extra**: `procedures/recommended-settings.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Major 8 件、Minor 1 件。Critical と断定できる証拠は無し）

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| **H1：生産レシピ別プロファイル** | **採用。ただし根拠の整理は要再検証／Major** | SERN の引用表は **2D 有限ベース限定**（`plans/active/tooling-nozzle-sern-3d.md:591`）。「上限 3」の判断は同ファイル `:627` で保留され、後続の緩和比較も `:761` では両者 **NOT CONVERGED** と記録されている。さらに壁解像 3D には CFL 0.5 で発散した記録がある（`:1258`）。**SERN 一括の検証済み設定にしない**。次元・物理・壁処理・メッシュ適用条件・段を持つプロファイルに分ける。Euler の `nStepInner: 5→4` も、case/16 SST の根拠だけで検証済みとは扱わない。 |
| **CFL の範囲** | **共通範囲外を無条件 FAIL は却下／Major** | `procedures/recommended-settings.md:61` の CFL は「目安」。一方、plan `:112` はこれを必須範囲にしている。**推奨範囲は WARN、選択した生産レシピからの無宣言変更は FAIL** と分ける。通常運用の例外は専用プロファイル、個別試験は理由付き逸脱にする。 |
| **`space.pRef` 必須化** | **要再検証／Major** | `procedures/recommended-settings.md:251` は「動作静圧」としか定義していない。入口全圧・出口静圧・平均静圧のどれを採る根拠も不足。**今回の共通必須化から外す**。先に圧力流束・境界・hoop 源での差し引きの整合、一定圧自由流の残差、代表ノズルの全残差床と定常目的量への影響を確認する。 |
| **段階起動の一律置換** | **要再検証／Minor** | `procedures/solver-settings.md:29` は起動時 `convMethod: 0, limiter: 2`、`recommended-settings.md:119` は `0,0` で、正本同士も不一致。**limiter は当面 2 を維持**し、段の数値を明示プロファイルへ移す。2000/3000 step や mid の次数は別の挙動変更として検証する。指定された資料だけでは、全関連経路で limiter の実効差が無いとは断言できない。同一メッシュの `restart_field.py` 化は採用。 |
| **壁直交性と AR 緩和** | **機械検査は採用、角度の即決は却下／Major** | `check_mesh_quality.py:273` は AR・skew のみ。既存関数を実行すると、**壁傾斜 40°、AR 4000、skew 0.444444** は上限 5000/0.9 を通った。一方、node は LSQ 勾配・エッジ中点再構成（`procedures/solver-settings.md:471`）なので、primal の角度だけから精度限界 5°は導けない。第一層の全対象辺と、**AR 緩和を受ける全セル・全層**を確認し、角部と通常壁を区別する。φ・内角許容は壁応力／熱流束等の誤差予算から定める。条件を証明できないセルには AR 緩和を認めない。 |
| **H2：`run_case.sh` が全経路を保証** | **前提は却下、設計チェーンの設置場所としては採用／Major** | `hook_forge_guard.py:12` 自身が変数展開による呼び出しを検出できないと明記し、`:32` はコマンド内の文字列でも除外する。**全直接起動を防ぐ仕組みではない**。今回の対象を管理された設計チェーンに限定し、各段の起動をラッパーへ統一する。forge 本体の証明書必須化は今回行わず、「全起動を保証」という記述を直す。 |
| **H3：config からプロファイルを推定** | **却下／Major** | `runner_axismach.py:1584` の前段生成は次数・CFL・反復数を変更するが、提案された推定項目である物理モデル・軸対称・壁種は本段と同じ。**物理の分類から実行目的・段は復元できない**。手書きも `profile・version・stage` の宣言を必須にし、config は宣言の適用条件と照合する。多重継承の衝突も暗黙上書きせず拒否する。 |
| **起動時 YAML と mesh SHA だけの照合** | **却下／Major** | `nodeInletCornerWall` は変換時に焼き込まれ、起動時には効かない（`procedures/solver-settings.md:429`）。また plan `:139` のキャッシュキーには、壁判定に使う BC・閾値・検査器版が無い。**変換時の実効設定・BC・変換器と生成 HDF5 の対応を記録・検証**し、品質キャッシュも BC・検査条件・検査器版に束縛する。未知／誤配置キーは現状 WARN・終了 0（`check_solver_config.py:117`, `:242`）なので、管理対象の起動では FAIL にする。YAML 重複キー、検査不能、証跡不明も通さない。 |
| **§6.3：短い起動試験だけで既定変更を承認** | **却下／Major** | plan `:222` の条件では、安定して別の解へ進む変更を検出できない。SERN の記録自身が「30 step の順位と4000 step の順位が異なる」としている（`tooling-nozzle-sern-3d.md:754`）。**起動健全性と生産受入を分ける**。前者は全段の実効設定・有限性・正値性、後者は同一作用素区間の収束 VERDICT と目的量の準定常 VERDICT、事前登録した許容差で判定する。ベルと SERN 3D を除外するには実効設定・経路の同等性の証明が必要。現資料から十分な固定 step 数・最小本数は決められない。 |

結論: **数値既定を一括変更する前に、plan §4.3 を「明示した profile・stage と、変換済みメッシュの実効条件に束縛した関門」へ修正し、CFD なしの負例試験を先行させる。**

第 1 仮説: 現案の関門は、設定表への適合を確認しても、実際に使うメッシュの不適合を通してしまう。 **確度: 高**
  根拠: `procedures/solver-settings.md:429` の変換時専用設定と、`plans/active/tooling-generator-defaults-and-launch-gate.md:139` の mesh SHA のみのキャッシュ。既存品質関数でも、上記の傾斜 40°・AR 4000 の反例を再現した。
  反証条件: 起動 YAML だけを修正した旧変換メッシュを拒否し、同じ mesh SHA でも BC・判定条件が変われば再判定する実装と負例試験が示されること。

第 2 仮説: config からの推定は、前段と本段を誤分類する。 **確度: 高**。`runner_axismach.py:1584` の段生成で、提案された推定項目が変わらない。

第 3 仮説: 生産設定の不整合の一部は、適用範囲と撤回履歴を落として表だけ転記したことに由来する。 **確度: 中**。SERN の 2D 有限ベース表、撤回、3D 壁解像の記録が一括扱いされている。現行 run による再確認は未実施。

判別 A/B: **壁直交性キャッシュの単体試験を 1 本**行う。固定した傾斜 40°・AR 100 の合成メッシュで、A は対象境界を `slip` として判定・キャッシュし、B は **その境界の `kind` だけを `wall` に変更**して再判定する。forge は **0 step**。試験用 φ_max は 5°に固定するが、生産閾値の根拠とは扱わない。
  → **A が対象外・通常品質 PASS、B が壁直交性 FAIL**なら、BC を見落とすキャッシュ仮説は反証される。
  → **B が A の判定を再利用して PASS**なら、mesh SHA だけの契約は不十分。見る量は対象壁数・最大 φ・キャッシュ使用有無・VERDICT。関門未実装のため、この A/B は未実施。

やらない方がよいこと: `nStepInner`・CFL・緩和・`pRef`・壁格子を同時に変更し、短時間 NaN が出ないことだけで新しい生産既定を承認すること。`FORGE_LAUNCH_GATE=warn` で禁止キーや無効メッシュまで一括解除すること。

呼び出し側の前提への異議: **生成器の直書きと値のずれは再現できた。しかし「検証済み生産設定」「全経路の通り道」「短い起動試験で既定変更を承認可能」は受け入れない。** さらに `check_mesh_quality.py:283` は少数の閾値超過を `SOFT-PASS`・終了 0 にするため、終了コードだけでは「全セルが緩和条件を満たす」を保証できない。新関門では扱いを明記する必要がある。

不足情報: このツリーの case/45・46・70 には `residual_history.csv` と `CONVERGENCE_VERDICT.txt` が無く、生産結果を独立再検証できなかった。必要なのは、採用予定プロファイルごとの具体的 run・実効設定・判定区間・両 VERDICT、`pRef` と limiter の関連実装、壁直交性と目的量誤差を結ぶデータ。**ファイル変更・forge 起動なし。plan 未反映（依頼どおり呼び出し側で反映）。**
