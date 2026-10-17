# codex レビュー: tooling-generator-defaults-and-launch-gate (plan)

- **plan**: [`plans/active/tooling-generator-defaults-and-launch-gate.md`](../../plans/active/tooling-generator-defaults-and-launch-gate.md)
- **stage**: `plan`
- **date**: 2026-10-14
- **commit**: `7221fb9b` (feature/generator-defaults-launch-gate)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m1
- **extra**: `procedures/recommended-settings.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

数値既定の変更を分離し、先に生成・投入経路を固める方針は妥当です。
ただし、現案には正常入力の誤拒否と、変換済みメッシュ・実効設定の取り違えを防げない箇所があります。
以下の Major 6 件を実装前に修正してください。ファイルは変更しておらず、**plan 未反映**です。

課題自体は未解決です。`run_case.sh` に投入前検査はありません。一方、重複・merge キー拒否と YAML の位置による編集は、accepted の `tooling-rerun-conditions` で既に実装されています。新規実装ではなく [`yaml_strict.py`](/home/sano/work/forge-gate/solver_density_cuda/tools/yaml_strict.py:57) を共用すべきです。`campaign-recipe` から数値設定部分を切り出す範囲分担も妥当ですが、起動記録との結び付けを後回しにするなら、下記 3 の制約が必要です。

1. **Major — 既存チェッカーの WARN をそのまま FAIL にすると、正しい設定まで拒否する。**

   **根拠:** 現行 `unknown_keys()` を実行すると、次の誤検出を再現しました。

   - `conjugate.mode: fem2d` →「正しくは `cj.mode`」
   - `conjugate.gate: {eps_rel: 0.01}` →「solverConfig が読まないキー」

   ソルバは両方を明示的に受理します（[solverConfig.cpp:541](/home/sano/work/forge-gate/solver_density_cuda/input/solverConfig.cpp:541)、[同:558](/home/sano/work/forge-gate/solver_density_cuda/input/solverConfig.cpp:558)）。チェッカーは C++ の正規表現抽出を受理キー一覧として使っており、その一覧が不完全です（[check_solver_config.py:77](/home/sano/work/forge-gate/solver_density_cuda/tools/check_solver_config.py:77)）。

   **対案:** §5 の厳格化より先に、受理キー一覧の完全性を修正・検証する工程を追加してください。未知キーの負例だけでなく、正当な CHT・TP・凝縮・周期・非定常設定を含む正例を必須にします。警告を一律に失敗へ変えるだけでは不可です。

2. **Major — メッシュの変換来歴に必要な入力と、成果物との結び付けが不足する。**

   **根拠:** [plan:141](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:141) の記録項目には、少なくとも次が欠けています。

   - `wallDistExtraPhysIDs`：変換時の `wall_dist` を変える（[calcWallDistance_kdtree.cpp:128](/home/sano/work/forge-gate/solver_density_cuda/input/calcWallDistance_kdtree.cpp:128)）。
   - 入口境界の種類・physID：壁だけでなく、`inlet_*` 判定も半割面の所属変更を決める（[gmshReader.hpp:1601](/home/sano/work/forge-gate/solver_density_cuda/mesh/gmshReader.hpp:1601)、[同:1625](/home/sano/work/forge-gate/solver_density_cuda/mesh/gmshReader.hpp:1625)）。
   - 記録が対応する変換成果物のハッシュ：品質キャッシュのハッシュだけでは、別メッシュの sidecar を取り違えたことを検出できません。

   また、入力 HDF5 は段間 restart で `/VALUE` を更新します（[restart_field.py:96](/home/sano/work/forge-gate/solver_density_cuda/tools/restart_field.py:96)）。変換時のファイル全体 SHA を固定すると正常な restart まで拒否します。

   **対案:** 変換が依存する設定・境界分類を棚卸しし、**形状・接続・境界所属・`wall_dist` の不変部分のハッシュ**と結び付けてください。初期場のハッシュは別に記録します。負例には sidecar の取り違え、入口分類だけの変更、`wallDistExtraPhysIDs` だけの変更を追加し、保存量だけを更新した正常 restart は通す試験を置きます。

3. **Major — config のバイト一致では、プロファイルの実効値を固定できない。**

   **根拠:** §5.1 #9 はバイト一致を要求しますが、省略値は実行バイナリで決まります。例えば `scalarGradient` は node で省略時 `lsq` に解決されます（[solverConfig.cpp:249](/home/sano/work/forge-gate/solver_density_cuda/input/solverConfig.cpp:249)）。`run_case.sh` は `FORGE_BIN` による別ビルドの実行も認めます（[run_case.sh:28](/home/sano/work/forge-gate/solver_density_cuda/tools/run_case.sh:28)）。

   現案の `RECIPE.json` はプロファイルの版と config SHA を持ちますが、対応バイナリ・省略値の解決規則との照合を定めていません。この問題は後回しの [launch-binding plan:31](/home/sano/work/forge-gate/plans/active/tooling-stage-manifest-launch-binding.md:31) と直接重なります。

   **対案:** 第 1 段は、プロファイルに対応バイナリ SHA、プロファイル内容ハッシュ、省略キーの期待実効値を記録し、対応不明のバイナリを拒否する方式に絞ってください。バイト一致は「同一バイナリでの生成器互換性」の試験として残します。全起動記録の改修を先行させる必要はありませんが、バイナリを無条件に自由化したまま実効値保証を名乗るべきではありません。

4. **Major — 全起動へ強制する範囲に対し、プロファイルの収録範囲と移行経路が不足する。**

   **根拠:** 宣言必須は全 `run_case.sh` 利用者へ及びますが、正例・棚卸しは主にノズル生成器と case/45・46・70 です（[plan:130](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:130)、[同:227](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:227)）。標準検証には衝撃管・KEEP 周期箱・凝縮などもあります（[verification/README.md:11](/home/sano/work/forge-gate/procedures/verification/README.md:11)）。

   さらに、古いメッシュについて「いつ FAIL にするか」はありますが、変換来歴を裏付けられないメッシュをどう移行するかはありません。待機期間だけでは記録は復元されません。

   **対案:** 強制開始条件を日付だけで決めず、対象経路の棚卸し完了にしてください。標準検証・手書き実験用の宣言方式も定義し、推奨適合と実験的逸脱を区別します。旧メッシュは元入力と変換条件を検証できる場合だけ移行し、根拠がなければ再変換を要求してください。現在の config から過去の来歴を推定して埋める移行は禁止します。

5. **Major — 段階起動の実体と、宣言・検証の単位が揃っていない。**

   **根拠:** 計画の段は `soft/mid/main` ですが、実コードには SERN の層流暖機・暖機 ramp・作動点変更後の適応段があります（[runner_sern.py:910](/home/sano/work/forge-gate/design/forge_design/evaluate/runner_sern.py:910)、[同:926](/home/sano/work/forge-gate/design/forge_design/evaluate/runner_sern.py:926)）。axis-Mach にも本段次数のまま進む ramp があります（[runner_axismach.py:1593](/home/sano/work/forge-gate/design/forge_design/evaluate/runner_axismach.py:1593)）。

   また、SERN・WT・walldriven には現在 `StageManifest` の登録がありません。`RECIPE.json` を書くだけでは、§6.4 の `check_convergence.py --segment` は成立しません（[check_convergence.py:219](/home/sano/work/forge-gate/solver_density_cuda/tools/check_convergence.py:219)）。

   `interp_field` → `restart_field` も実際に初期場を変える修正です。config の一致試験では検証できず、ビット一致だけでは負密度・Inf を排除できません。axis-Mach には既に別の段終了ゲートがあります（[runner_axismach.py:1575](/home/sano/work/forge-gate/design/forge_design/evaluate/runner_axismach.py:1575)）。

   **対案:** 実在する全起動経路について、段のプロファイル、最終 config、`RECIPE`、履歴退避、`StageManifest` を一組で定義してください。第 1 段の試験も「基本 config 5 個」から「全段の起動列」へ拡張します。restart は保存量・種属性・`wall_dist` 保持を検査し、欠損保存量・Inf・負密度なら次段を起動しない負例を追加してください。

6. **Major — §6.3 のキャッシュ試験は、判定仕様と矛盾し、失敗原因も分離できない。**

   **根拠:** §4.4 は AR ≤ 1000 の大きい φ を WARN としますが、§6.3 は AR 100・φ 40° を FAIL と期待します（[plan:158](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:158)、[同:239](/home/sano/work/forge-gate/plans/active/tooling-generator-defaults-and-launch-gate.md:239)）。

   既存 `cell_metrics()` に合成四角形を与えた実測は **AR = 100、skew = 0.444444、φ = 40°**。第 1 段の規則なら品質判定は WARN 相当です。また、起動関門全体を通す試験なら、`slip` → `wall` は変換来歴の不一致でも止まるため、FAIL だけではキャッシュ無効化の成功を証明できません。

   **対案:** この試験は品質チェッカー単体で行い、A は対象壁 0、B は対象壁あり・φ 40°・WARN、かつキャッシュ再利用なしを期待値にしてください。来歴不一致の拒否は別試験にします。キャッシュキーには境界種類だけでなく **physID と壁面集合の対応**を明記し、種類の集合が同じまま壁にする physID を入れ替える負例も追加してください。

7. **Minor — 第 2 段の基準 run を「生産設定」という名称だけで受入基準にしない。**

   **根拠:** 別作業ツリー `/home/sano/work/forge` の既存データを読み取り、現行ツールで再判定しました。

   | run（同ツリーのリポジトリ相対パス） | 再判定 |
   |---|---|
   | `case/46.sern_design/run_0352_long_cfl1ref/` | `NOT CONVERGED`。末尾 `rms_roe = 2.44`、低下 0.4 桁 |
   | `case/46.sern_design/run_0353_long_relax05/` | `NOT CONVERGED`。末尾 `rms_roe = 81.7`、低下 0.0 桁 |

   前者のローカルデータに対する `check_quasisteady.py --quantity shock,machmax,pmax` は **`TRANSIENT-UNSETTLED`**（対象スナップショット 2 枚、必要数 4）でした。このデータだけから「収束解が緩和で変わった」とは判断できません。run の索引は [case/46 README](/home/sano/work/forge/case/46.sern_design/README.md) です。

   **対案:** §5.1 #12 では、基準 run の区間・収束 VERDICT・目的量の準定常 VERDICT・保存データの所在まで事前登録してください。壁直交性の受入には `check_wall_resolution.py` と局所壁応力・熱流束の比較も含めます。現時点の承認範囲は第 1 段までとし、第 2 段は登録済み基準を確認してから進めるべきです。

**推奨は、第 1 段を「現行挙動を固定した投入契約」として先に完成させることです。** 実装前の修正順は、①受理キー一覧の修正、②メッシュ来歴の完全化、③バイナリとプロファイルの固定、④全起動経路の移行設計、⑤全段の記録・restart 検証、⑥キャッシュ試験の修正。数値既定と壁法線格子の既定化は、その後の受入に分離してください。

指摘数: Critical 0 / Major 6 / Minor 1
