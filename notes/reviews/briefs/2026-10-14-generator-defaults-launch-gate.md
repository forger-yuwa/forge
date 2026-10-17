# 諮問ブリーフ: 設計チェーンの生成器の既定を推奨レシピから組み、投入前の関門で設定とメッシュを止める (2026-10-14)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1: plan §4・§6 を新規に書く)
- plan: `plans/active/tooling-generator-defaults-and-launch-gate.md` (draft。§4 は起案、§4.8 が問い)
- ユーザ要求 (2026-10-14、原文): 「なんか、汎用ランナーが、古い解析設定を生成するみたい。推奨設定を、あほなaiでも実施してくれるように改修できる?」「あと、壁面メッシュが、壁に対して直交しないものが、あほなaiで生成されてる。対処法ある？」
- 作業ツリー: `forge-gate` / `feature/generator-defaults-launch-gate` (`origin/main` f66113ee から)。

## 観測事実

すべて読み取りの調査 (2026-10-14)。行番号は `feature/nozzle-wall-fit-and-pipeline` の作業ツリーのもので、`origin/main` ではずれる。plan §3 に表がある。

1. 生成器は 5 つ: `runner.py` `_solver_config` (ベル)、`runner_wt.py` `_config_euler` (cell 固定)・`_config_euler_node`・`_config_sst_node` (axis-Mach の生産 NS もこれ)、`runner_sern.py` `_solver_config` (SERN 2D・3D)。値はすべて直書き。`procedures/recommended-settings.md` を読むコードは無い。
2. 推奨とのずれ (plan §3.2): `nStepInner` 全生成器 5 (推奨 4、書かないと solver 既定 1)。`implicitRelax` は指定時のみ (推奨 0.7、solver 既定 1.0)。NS 本段 CFL の既定 1.0 (axis-Mach `prepare_ns`・intake レシピ・SERN 3D、推奨 4〜6)。`nodeInletCornerWall` は R・WS で欠落 (solver 既定 0、推奨 1)。`turbulentPrandtl` は R・S で欠落 (solver 既定 0.85、推奨 0.9)。軸対称の `pRef` はどこも書かない (solver 既定 0)。WE と SERN の既定は cell (ユーザは 2026-09-16 に cell を使わないと決定)。段階起動は 3000 step・soft の limiter 2・Euler と SERN の mid は 2 次 (推奨 §1.2 は 2000 step・limiter 0・1 次)。runner_wt と walldriven は同一メッシュの段間に `interp_field.py` を使う (AGENTS.md は `restart_field.py` を要求)。
3. 推奨の側の不整合 (plan §3.3):
   - SERN の生産設定は `cfl_main 1.0`・`implicitRelax 1.0`・`nStepInner 5` (`plans/active/tooling-nozzle-sern-3d.md` の生産設定の表、「上限 3、4 で発散」)。推奨 §1 は CFL 4〜6・緩和 0.7・`nStepInner` 4。
   - 推奨 §4 (Euler 設計評価、2026-09-04) は `nStepInner 5`。§1 の改定 (2026-09-12) より古い。§4 は「`runner.py` が生成する config が正本」と書くが、`runner.py` は Euler を作れない。
   - case/45 の生産 NS (B0) は `main_line` で CFL 4・緩和 0.7 を明示、`nStepInner` は書かない (= 5)。
   - `rerun_conditions.py` の `recommended_stages` は NS の Pt 変更で CFL 1.0・60000 step を勧める。
   - `pRef` の「動作静圧」の決め方は書かれていない。
4. 既存の照合: `check_solver_config.py` は削除キー・`bndFirstOrder`・dual-time の緩和などを見るが推奨値は見ず、生成器からも `run_case.sh` からも呼ばれない。`run_case.sh` は来歴と収束判定だけ。design の `run_forge` は `run_case.sh` を呼ぶので、`run_case.sh` が全経路の通り道 (エージェントの Bash は `hook_forge_guard.py` が `run_case.sh` 経由を強制)。`run_case.sh` は実行中に編集すると走行中の run を壊す (ファイル冒頭の注記、2026-09-20 の事故)。
5. メッシュ: `mesh2d.wall_normal_layer` (近壁を壁法線に、2026-10-08) は `feature/nozzle-core-grid` 系にあり `origin/main` に無い。既定 `None` で j 方向の線は半径方向。縮流部の壁は最大 40° 傾くので近壁セルの内角は 50°、skew 0.44 で `check_mesh_quality.py` (skew ≤ 0.9) を通る。AGENTS.md の AR ≤ 5000 の緩和は「壁法線の構造層 (内角 90° 近傍)」に限るが、`--ar-max 5000` は条件を検査しない。

## 期待値と出典

- 推奨設定の正本は `procedures/recommended-settings.md` (AGENTS.md「計算・実行ルール」)。
- ユーザ決定 (2026-10-07): 標準と決めた方式はコードの既定にし、過去 YAML の再現のために旧既定を残さない。
- AGENTS.md: 同一メッシュの restart は `restart_field.py`、メッシュ品質は AR ≤ 1000・skew ≤ 0.9、緩和は壁法線の構造層だけ AR ≤ 5000。

## 再現条件

CFD は不要。生成器を呼んで出力 YAML を読めば再現する (`design/tests/` の既存の試験が `nStepInner` 5 を固定している: `run_main_line_tests.py`・`run_stop_rule_tests.py`)。

## 実施済みの操作と結果

- 生成器と推奨の突き合わせ (読み取りのみ)。結果は plan §3。
- plan の起案 (§4.1〜§4.7)。実装・run はしていない。

## 仮説 (呼び出し側の案。誤りなら指摘してほしい)

- H1: 推奨の値は「共通 1 つ」ではなく「生産レシピごとのプロファイル + 根拠」で持つべき。§1 の共通値を全生成器に当てると、SERN のように検証済みの生産設定を壊す。
- H2: 関門は `run_case.sh` → `launch_gate.py` で足りる (全経路の通り道)。forge 本体に関門の記録の照合を入れるのは過剰。
- H3: 手書き config のプロファイルは config から推定できる。推定できないときだけ FAIL にすればよい。

## 問い (plan §4.8)

1. **正本の優先順**: §1 の共通値と生産レシピ (SERN の CFL 1.0・緩和なし、B0 の明示値) が食い違うとき、各生成器の既定をどちらにするか。プロファイルを生産レシピごとに分ける案 (plan §4.1) でよいか。§4 (Euler) の `nStepInner` 5 を §1 の 4 に揃えるか。
2. **範囲の扱い**: CFL「4〜6」のような範囲を関門の FAIL の境界にするか WARN にするか。範囲の外が要る生産設定はプロファイルで表すか逸脱 (理由つき) で表すか。
3. **`space.pRef`**: 軸対称ノズルで何を入れるべきか。生成器は書いておらず、case/45 の生産 run も `pRef` 0。必須化すると生産の結果が変わる。必須にする前に何を確かめるべきか。
4. **段階起動の正本**: §1.2 (2000 step・limiter 0・mid 1 次) と現行の生成器 (3000 step・limiter 2・mid 2 次) のどちらに揃えるか。convMethod 0 のときの limiter の値に実効差があるか。
5. **壁直交性の判定**: 近壁の格子線と壁法線の角 φ の上限、見る層の数、AR 緩和で認める内角の幅。node (median-dual) の離散化で近壁の非直交が誤差に効く機構 (勾配の再構成・壁の双対面・`wall_dist`) と照らして、根拠のある値の決め方。
6. **関門の置き場所**: H2 でよいか。直接起動を塞ぐべきか (試験・性能計測の起動との両立)。
7. **プロファイルの推定**: H3 でよいか。誤推定 (例: 軸対称の Euler を NS と取る、段の config を本段と取る) をどう防ぐか。宣言を必須にする案と比べてどちらが「判断力の低いエージェント」に対して堅いか。
8. **既定の変更の確かめ方** (plan §6.3): 既定が変わる生成器ごとに、何を何 step 回して何を合格にすれば足りるか。最小の本数で。

## 読んでよいファイル

- `plans/active/tooling-generator-defaults-and-launch-gate.md` (全体)
- `procedures/recommended-settings.md` (§1・§1.2・§2・§4・§5・§8・§9)
- `design/forge_design/evaluate/runner.py`・`runner_wt.py`・`runner_axismach.py`・`runner_sern.py`・`runner_sern3d.py`・`runner_walldriven.py`
- `solver_density_cuda/tools/check_solver_config.py`・`check_mesh_quality.py`・`run_case.sh`・`hook_forge_guard.py`・`rerun_conditions.py`
- `plans/active/tooling-nozzle-sern-3d.md` (生産設定の表)、`plans/active/tooling-design-problem-campaign-recipe.md` (§4)
- `AGENTS.md`
