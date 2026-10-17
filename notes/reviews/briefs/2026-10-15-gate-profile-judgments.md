# 諮問ブリーフ: 機種のプロファイルと関門の実装で出た判断 4 件 (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1: §4 の方針の詳細を決める前、条件 4: plan §5.1 に無い修正の前)
- plan: `plans/active/tooling-generator-defaults-and-launch-gate.md` (§4.0〜§4.7、§5.1 #6・#9・#10・#11・#12、§6.2・§6.3)
- 検証対象: ブランチ `feature/generator-defaults-launch-gate`、HEAD ebe0d6b6 (#9 = c9b28d25、#6 = 5637d090、#5 = 6429e86f、#8c の取り込み = ebe0d6b6)
- 前の諮問: `notes/reviews/2026-10-14-generator-defaults-launch-gate-diagnose.md` (宣言を必須に、生産の受入と起動の健全性を分ける)、`notes/reviews/2026-10-14-wall-orthogonality-judgments-diagnose.md`

## 観測事実 (実装の結果)

- #6: `solver_density_cuda/tools/recipes/recommended_settings.yaml` (版 1) に機種 11 のプロファイル。値は各生成器を代表の問題 YAML で偽の forge に通し、全段の config から書き写した (`recipes/regen_profiles.py`)。照合は `recipe_settings.py check`: `fixed` の不一致とプロファイルに無いキーは FAIL (理由つきの逸脱で通す)。`problem_keys` (物性・組成・本段の長さ・出力の間隔・pRef・リミッタの基準値・kInit/omegaInit・extraFields・condensation) は照合しない。
- #6 の照合の規則で、**`fixed` のキーを config に書いていなくても、solver の既定 (`config_key_inventory.py` がソースから抜いた値) が同じなら PASS** にした。`nodeInletCornerWall` は変換器が読むので既定を抜けず、差の表では「不明」。
- #6 の書き写しは「生成器 × 代表の問題」で決まる。SERN の cfl の梯子 (2D 0.1/0.1/0.5、3D 0.2/0.2/0.25) と `cfl_main` は問題 YAML の値を写した。同じ機種で別の問題 YAML を使う (例えば cfl_main を変える) と、今の規則では `fixed` の不一致 = 逸脱の宣言が要る。
- 代表の問題は追跡済みの YAML (case/40・41・45・46・68) を直接読む。これを誰かが編集すると `design/tests/run_recipe_profiles_tests.py` の (a) (台本の出力 = YAML) が落ちる。case/45 の YAML は他のセッションが編集中のことがある。
- #9 の棚卸し (`notes/investigations/2026-10-14-launch-path-inventory.md` §3) で、同じ格子への `interp_field.py` が 2 系統残っている: (i) `euler_loop` の継続チャンク `cont_<j>` (同じ nozzle.h5)、(ii) `opt/driver`・`opt/polish` の CPG のベル (同じ幾何から作り直した格子; 生産のベルは TP の経路)。SERN の段間 `restart_by_index` には段終了ゲート (必要保存量・有限・ρ>0) が無い (化学種の照合と座標一致ノードの扱いはある)。
- `check_solver_config.py` の passiveFct の WARN が、そのキーを書いていない定常の config すべてに出る (誤報)。今は `recipe_settings` の側で落としている。
- 観察: SERN 2D の生産 YAML は `node_inlet_corner_wall` を書いていない (runner のコメントは「生産 YAML は 1」)。

## 期待値と出典

- ユーザ決定 (2026-10-14): 機種ごとの推奨 = 最新の生産レシピ、生成器の既定 = 機種の推奨。推奨から外すのは理由だけで可 (根拠は任意、「試し計算は手軽に」)。棚卸しの完了後は宣言の無い手書きの run を `run_case.sh` で止める (それまで WARN)。
- plan §4.3 (4a): 「省略に頼らない」(前の plan 段レビュー M3 の採用)。
- AGENTS.md: 同一メッシュの restart は `restart_field.py`、`interp_field.py` は使わない。

## 呼び出し側の仮説 (案)

1. **`fixed` の省略**: 生成器が作った config (RECIPE.json の出所 = generator) では省略を FAIL (生成器に書かせる)、手書きの run では WARN。solver の既定が将来変わると、省略した config の実効値が黙って変わるため。`nodeInletCornerWall` はメッシュの来歴 (#7) の実効条件で照合する。
2. **問題 YAML の数値の上書き** (cfl の梯子・cfl_main): プロファイルでは `fixed` から `advisory` (範囲つき) に移し、範囲外は WARN、プロファイルとの差は生成器が RECIPE.json に「問題 YAML <path> の <key>」を理由として自動で逸脱を書く。`problem_keys` には入れない (照合から消えるため)。
3. **代表の問題 YAML**: `design/tests/fixtures/recipe_profiles/` に写しを置き、sha256 と元のパスを記録。試験は写しを読む。元の YAML と写しの違いは別の検査で WARN (書き写しの更新の契機)。
4. **同一格子の `interp_field.py` の残り**: (i) euler_loop の継続チャンクは #11 で `restart_field.py` + #9 の段終了ゲートに置き換える (AGENTS.md の規則違反なので必須)。(ii) CPG のベルは、格子の生成が決定的なら `restart_field.py` の座標一致の検査で置き換え、決定的でなければ `interp_field.py` を理由つきで残す。SERN の `restart_by_index` に段終了ゲートを足すのも #11。
5. passiveFct の誤報は #10 で `check_solver_config.py` 本体を直す。

## 問い

1. 案 1〜5 の採否。特に案 1 (生成器だけ省略を FAIL) と案 2 (`advisory` に移す + 自動の逸脱の記録) が、「省略に頼らない」「宣言を必須に」と「試し計算は手軽に」を両立するか。もっと単純な形があるか。
2. 案 4 の (ii): 作り直した格子の決定性を `restart_field.py` の座標一致で確かめるのは十分か。
3. #10 (関門) の前に、他に決めておくべきことはあるか。

## 読んでよいファイル

- plan 全体、`notes/investigations/2026-10-14-launch-path-inventory.md`
- `solver_density_cuda/tools/recipe_settings.py`、`solver_density_cuda/tools/recipes/recommended_settings.yaml`、`solver_density_cuda/tools/recipes/regen_profiles.py`
- `solver_density_cuda/tools/check_solver_config.py`、`solver_density_cuda/tools/config_key_inventory.py`
- `design/forge_design/evaluate/runner_wt.py`・`runner_walldriven.py`・`runner_sern.py`・`euler_loop.py`、`design/forge_design/opt/driver.py`・`polish.py`
- `design/tests/run_recipe_profiles_tests.py`・`run_launch_path_tests.py`、`solver_density_cuda/tests/unit/test_recipe_settings.py`
- `procedures/recommended-settings.md`
