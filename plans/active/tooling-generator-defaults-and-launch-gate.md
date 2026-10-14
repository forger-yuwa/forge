# 生成器の既定を推奨レシピから組み、投入前の関門で止める (設定とメッシュ)

## メタ

- **area**: `tooling` (`design/forge_design/evaluate/runner*.py`、`design/forge_design/meshing/mesh2d.py`、`solver_density_cuda/tools/check_solver_config.py`・`check_mesh_quality.py`・`run_case.sh`)
- **status**: `draft` (2026-10-14 codex diagnose 反映済み、plan 段レビュー前)
- **related_docs**:
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) (推奨設定の正本。本 plan で値の正本を機械可読のファイルへ移し、本文はそこへの解説にする)
  - [`procedures/calculation-workflow.md`](../../procedures/calculation-workflow.md) 「メッシュ品質チェック」
  - [`methods/design/overview.md`](../../methods/design/overview.md) (設計チェーンが forge の入力を作る箇所)
- **related_plans**:
  - [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (draft。recipe が「推奨設定のどの節に従うか」を持つ構想。本 plan はその数値設定の部分だけを先に実体化する。§4.7)
  - [`tooling-stage-manifest-launch-binding.md`](tooling-stage-manifest-launch-binding.md) (draft・後回し。起動記録の扱いが重なる)
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) §5.1 #15 (`mesh2d.wall_normal_layer` の導入元)
  - [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (SERN の生産設定 `cfl_main 1.0`・緩和なし・`nStepInner 5` の根拠)
- **created**: `2026-10-14`
- **owner**: `M6ノズルの設計` セッション (作業ツリー `forge-gate` / `feature/generator-defaults-launch-gate`、`origin/main` f66113ee から)

## 1. 目的

設計チェーンの生成器 (`runner*.py`) は推奨設定の正本を読まず、値を各ファイルに直書きしている。正本を改定しても生成器には伝わらない (例: `nStepInner` は 2026-09-12 に 5 → 4 と改定されたが、全生成器が 5 のまま)。投入の前に推奨値を照合する関門も無い。メッシュも同様で、傾斜壁で壁に直交しない近壁格子が判定を通って計算に入る。

ユーザの要求 (2026-10-14): 「推奨設定を、あほなaiでも実施してくれるように」「壁面メッシュが壁に対して直交しないものが生成されている」。手順書を読む・読まないに依存せず、**生成器は推奨値しか作れず、手書きの入力は投入時に止まる**状態にする。

## 2. スコープ

- **やる** (第 1 段、数値の既定は変えない):
  1. 生産レシピごとのプロファイル (適用条件・段・根拠・推奨との差) を機械可読のファイル 1 つに置き、値の正本にする (§4.1)。
  2. 生成器の直書きをやめ、共通の組み立て関数から config を作り、プロファイルの宣言を書く (§4.2)。
  3. `run_case.sh` の起動前に、宣言したプロファイル・段と変換済みメッシュの実効条件に照らして設定とメッシュを照合し、FAIL なら forge を起動しない (§4.3)。
  4. メッシュ判定に壁直交性を足し、AR の緩和条件 (壁法線の構造層) を機械判定にする (§4.4)。
  5. 同一メッシュの段間の引き継ぎを `restart_field.py` に揃える (§4.5)。
- **やる** (第 2 段): プロファイルごとに受入の試験を事前登録して回し、合格したものから既定を推奨の値へ切り替える。`mesh2d` の壁法線の近壁層の既定化もここ (§4.0・§6.4)。
- **やらない**:
  - ソルバ側の既定値の変更 (`solverConfig.cpp`)。書かないと推奨と違う値になるキーは、生成器と関門で書かせる。
  - campaign / recipe の契約全体 (`tooling-design-problem-campaign-recipe` の範囲)。
  - 過去の run の config の書き換え。関門は起動時だけに効く。
  - 3D 非構造 (四面体) の壁直交性。三角形・四面体の近壁は既存の AR ≤ 1000・skew ≤ 0.9 だけで見る。

## 3. 関連 docs と前提 (観測事実、2026-10-14 調査)

調査は `feature/nozzle-wall-fit-and-pipeline` の作業ツリーで行った (読み取りのみ)。`origin/main` の同じファイルは行番号がずれるが、内容は同じ系統である。

### 3.1 生成器

| 略号 | 関数 | 使う経路 |
| --- | --- | --- |
| R | `runner.py` `_solver_config`・`_bcond_config` | ベル型ノズル (`opt/driver.py`) |
| WE | `runner_wt.py` `_config_euler` (**cell 固定**) | `runner_wt.prepare(euler=True)`、`euler_loop` (cell 指定時) |
| WEn | `runner_wt.py` `_config_euler_node` | axis-Mach の Euler 本番、`runner_walldriven`、`euler_loop` |
| WS | `runner_wt.py` `_config_sst_node` | `runner_wt.prepare`、axis-Mach の生産 NS (`prepare_ns`) |
| S | `runner_sern.py` `_solver_config` | SERN 2D・3D |

この上に axis-Mach 側の書き換え (`_apply_gas_to_config`・CFL の置換・`_apply_main_line`・段の設定) が重なる。段の書き換えのうち SERN・walldriven・runner_wt・driver は正規表現か文字列置換で、書式が変わると黙って当たらなくなる。

### 3.2 推奨 (`recommended-settings.md`) とのずれ

| 項目 | 推奨 | 生成器 | 書かないときの solver 既定 |
| --- | --- | --- | --- |
| `time.nStepInner` | 4 (§1、2026-09-12) | 全生成器で 5。段の論理も 5 を前提 (axis-Mach は「5 のときだけ 10 に上げる」) | 1 |
| `time.deltaT.implicitRelax` | 0.7 (§1) | 指定したときだけ書く | 1.0 |
| NS 本段の CFL | 4〜6 (§1) | axis-Mach の生産 NS・intake のレシピ・SERN 3D の既定は 1.0 | — |
| `mesh.nodeInletCornerWall` | 1 (§1) | R・WS (生産 NS を含む) は書かない。S は既定 0 | 0 |
| `turbulence.turbulentPrandtl` | 0.9 (§2) | R・S は書かない | 0.85 |
| `space.pRef` (軸対称) | 動作静圧 (§5) | 軸対称の生成器はどれも書かない | 0 |
| `mesh.discretization` | node (§1、cell は使わない: 2026-09-16 ユーザ指示) | WE が cell 固定、S の既定が cell | cell |
| 段階起動 | 2000 step、soft・mid は 1 次・limiter 0 (§1.2) | 3000 step、limiter 2 のまま、Euler と SERN の mid は 2 次 | — |
| 同一メッシュの段間 | `restart_field.py` (AGENTS.md) | runner_wt・walldriven は `interp_field.py` | — |

廃止キー (§9) と `bndFirstOrder` を書く生成器は無い。`dt`・`dt_min`・`dt_max` の値は推奨と違うが、定常では効かない。

### 3.3 正本の側の不整合 (照合を機械化すると表に出るもの)

1. **推奨が 1 つでない**: §1 は CFL 4〜6 + 緩和 0.7・`nStepInner` 4。SERN の生産設定は検証の結果 `cfl_main 1.0`・緩和なし・`nStepInner 5` (CFL 4 で発散、`tooling-nozzle-sern-3d.md` の生産設定の表)。メモリ「SERN: CFL を上げる / implicit_relax は使わない」は緩和 0.7 が解を 0.05 % 動かすとする。§4 (Euler 設計評価) は `nStepInner` 5 のまま (2026-09-04、§1 の改定より古い)。case/45 の生産 NS (B0) は `main_line` で CFL 4・緩和 0.7 を明示し、`nStepInner` は書かない (= 生成器の 5)。
2. **範囲で書かれた値**: CFL は「4〜6」「0.5〜2 (TP)」「≤ 50 (易しい流れ)」など。
3. **決め方の無い値**: `space.pRef` の「動作静圧」(出口・外気・入口のどれか)。
4. **古い記述**: §4 は「`runner.py` が生成する config が正本」とするが、`runner.py` は Euler を作れない。§2 は「生成器は `evaluate.wall_treatment_sst` を読む」とするが、読むのは SERN だけ。§3 の化学種の例は生成器の書き方と違う。
5. **生成器ごとの独自推奨**: `tools/rerun_conditions.py` の `recommended_stages` は NS の Pt 変更で CFL 1.0・60000 step、Euler で CFL 2 を勧める。

### 3.4 既存の照合

- `check_solver_config.py`: 削除・改名キーを FAIL、誤った節・存在しないキーを WARN、`bndFirstOrder` を FAIL、dual-time の緩和などを検査する。**推奨値は見ない**。design の生成器からも `run_case.sh` からも呼ばれていない (skill `forge-config` の手順にあるだけ)。
- solver の起動時検査: 削除キー・node の旧キー・旧乱流キー・`limiterScaled 2` の拒否、`wallTreatmentSST 1` の警告。
- `run_case.sh`: 来歴 (`RUN_PROVENANCE.txt`) と収束判定だけ。設定・メッシュの照合はしない。design の `run_forge` は `run_case.sh` を呼び、エージェントの Bash は `hook_forge_guard.py` が `run_case.sh` 経由を強制する。ただし全起動の保証ではない (同フックは変数展開による呼び出しを検出できないと自ら書いている。diagnose 2026-10-14)。

### 3.5 メッシュ

- `mesh2d.wall_normal_layer` (近壁を壁法線に置く、2026-10-08) は `feature/nozzle-core-grid` 系 (c096d66c) にあり、`origin/main` には無い。既定は `None` (j 方向の線は半径方向)。縮流部の壁は最大 40° 傾くので、近壁セルの内角は 50° になる。
- `check_mesh_quality.py` の skew 閾値は 0.9。内角 50° のセルの skew は 0.44 で合格する。
- AGENTS.md の AR 緩和 (≤ 5000) は「内角 90° 近傍でスキューの無い壁法線の構造層」に限るが、`--ar-max 5000` はその条件を検査しない。

## 4. 設計方針 (2026-10-14 codex diagnose を反映、[記録](../../notes/reviews/2026-10-14-generator-defaults-launch-gate-diagnose.md))

### 4.0 二段に分ける (diagnose の結論)

- **第 1 段 (本 plan の中心)**: 値の正本 (プロファイル)・共通の組み立て関数・宣言・投入の関門・メッシュの壁直交性の機械判定。
  **プロファイルの値は、各生成器が今使っている生産の値をそのまま書き写す**。§1 の推奨と違う所は、プロファイルの中に
  「推奨との差 (未受入)」として記録し、関門の出力にも出す。数値の既定は変えない。例外は AGENTS.md の既存の規則の適用だけ
  (同一メッシュの段間を `restart_field.py` にする、§4.5)。
- **第 2 段 (プロファイルごとの受入)**: 推奨への切り替え (`nStepInner` 4、緩和 0.7、CFL、`nodeInletCornerWall`、`turbulentPrandtl`、
  cell → node、段の長さ) は、プロファイルごとに**受入の試験**を事前登録して回し、合格したものから既定を切り替える (§6.4)。
  起動の健全性 (NaN が無い) だけでは承認しない。安定したまま別の解へ進む変更を検出できないからである
  (SERN の記録に、30 step と 4000 step で順位が入れ替わった例がある: `tooling-nozzle-sern-3d.md` の該当節)。
- 理由: diagnose は「§1 の値を検証済みの生産設定に当てると壊す」(SERN の CFL 4 は 2D 有限ベースで発散) と「短い起動試験で既定を承認しない」を
  指摘した。判断力の低いエージェントへの堅さは、値そのものより「宣言と関門で止まること」から来るので、第 1 段だけで要求の大半を満たす。

### 4.1 値の正本: プロファイルの機械可読ファイル

- `solver_density_cuda/tools/recipes/recommended_settings.yaml` を置く (design からも tools からも読める場所)。読み込みと照合は `solver_density_cuda/tools/recipe_settings.py`。
- プロファイルは**適用条件**を持つ: 次元 (2D 平面・2D 軸対称・3D)、物理 (Euler / 層流 / SST、CPG / TP)、壁の扱い (slip / no-slip 断熱 / 等温)、
  メッシュの種類 (構造・非構造)、**段** (`soft`・`mid`・`main`)。生産レシピごとに分ける (例: `axismach_euler_node/main`、`axismach_ns_sst_b0/main`、
  `sern2d_finite_base_ns/main`。SERN を 1 つにまとめない。2D 有限ベースの設定を 3D や壁解像の 3D に流用しない)。
- 各キーは 3 種のどれかで持つ:
  - `fixed`: 宣言したプロファイルの値。違えば FAIL (理由つきの逸脱を宣言すれば通す)。
  - `advisory`: 推奨の範囲 (例: CFL 4〜6)。外れたら WARN。
  - `forbid`: 書けば FAIL (`bndFirstOrder`、§9 の削除キー、dual-time の `implicitRelax` など)。
- 各プロファイルに `source` (plan・節・日付) と、§1 との差の一覧 `differs_from_recommended` (キー・§1 の値・この値・受入の状態) を持たせる。
- 継承は 1 本の鎖に限る。2 つの親が同じキーに別の値を与えたら、黙って上書きせず読み込み時に例外にする。
- `recommended-settings.md` の値の表は `recipe_settings.py --doc` で生成した区間に置き換え、試験で YAML との一致を確かめる。

### 4.2 生成器: 共通の組み立て関数

- `design/forge_design/evaluate/config_builder.py` に `build_solver_config(profile, stage, problem, deviations=None) -> (text, record)` を置く。5 つの生成器はこれを呼ぶ。
- 値の出所を 3 つに分けて記録する: プロファイル、問題 YAML の物理条件 (気体・壁温・入口条件)、明示の逸脱 (キー・値・理由・根拠)。理由の無い逸脱は例外。
- 生成器は run ディレクトリに `RECIPE.json` (プロファイル・版・段・逸脱・生成した config の sha256) を書く。段ごとに config を作る生成器は段ごとに書く。
- 段の設定は YAML の位置で書き換える方式 (axis-Mach の `_cfg_set`/`_cfg_check`) に揃え、正規表現・文字列置換による書き換えを廃止する。

### 4.3 投入の関門: 宣言に束縛した照合

- **宣言を必須にする** (diagnose: config からの推定は前段と本段を区別できない。`runner_axismach.py` の前段は次数・CFL・反復数だけを変え、物理・軸対称・壁種は本段と同じ)。
  生成器は `RECIPE.json` を書く。手で作る run は `recipe_settings.py declare <run_dir> --profile P --stage S` で宣言する
  (宣言の時点で config を照合し、通らなければ書かない)。宣言の無い run は関門で FAIL (移行の扱いは §4.7)。
- `run_case.sh` は forge の起動の直前に `solver_density_cuda/tools/launch_gate.py <run_dir>` を呼び、終了コードが 0 でなければ forge を起動しない。
  `run_case.sh` への変更はこの呼び出しの数行だけにし、照合の中身は Python 側に置く (`run_case.sh` は実行中に編集すると走行中の run を壊す)。
- **範囲の記述**: `run_case.sh` は設計チェーンと、フックが強制するエージェントの Bash の起動の通り道だが、全起動の保証ではない
  (`hook_forge_guard.py` は変数展開による呼び出しを検出できないと自ら書いている)。本 plan はこの範囲を対象にし、forge 本体での照合はしない。
- 照合の中身:
  1. **設定**: 既存の `check_solver_config.py` の検査に加え、宣言したプロファイルとの照合。**未知のキー・誤った節のキー・YAML の重複キー・merge key は FAIL**
     (`check_solver_config.py` は今は WARN・終了 0)。検査できない入力も FAIL。
  2. **メッシュの実効条件**: 変換時にだけ効く設定がある (`nodeInletCornerWall` は変換時に焼き込まれ、起動時の値は効かない: `procedures/solver-settings.md` の該当節)。
     変換時の実効設定 (`discretization`・`nodeInletCornerWall`・壁の bcond の種類と physID・変換器の sha256) を `<mesh>.provenance.json` に書き、
     関門は run の config・bcond と照合する。記録の無いメッシュは FAIL (移行の扱いは §4.7)。
  3. **メッシュの品質**: §4.4 の判定。結果の再利用は (メッシュの sha256、壁の bcond の種類、閾値、検査器の版) が全部一致したときだけ。
  4. 結果を `LAUNCH_GATE.txt` (VERDICT・プロファイル・段・逸脱・WARN・メッシュの sha) に書き、`RUN_PROVENANCE.txt` にも 1 行写す。
- `FORGE_LAUNCH_GATE=warn` は `fixed` の不一致だけを WARN に落とす。`forbid`、宣言の欠落、メッシュの FAIL、検査不能は落とさない。

### 4.4 メッシュ: 壁直交性と AR の緩和条件

- `check_mesh_quality.py` に壁直交性を足す。壁の bcond (`wall`・`wall_isothermal`) の各境界辺 (3D は境界面) について、隣接セルの
  「壁上の節点から壁の外の節点へ向かう辺」と壁法線のなす角 φ を測る。四角形・六面体・プリズムを対象にし、三角形・四面体は対象外と出す。
- **AR の緩和の機械化**: AR が 1000 を超えるセルは、壁に接する構造層の中にあり、層の全セル・全層で φ と内角が条件を満たすときだけ 5000 まで認める。
  条件を確かめられないセルには緩和を認めない (`--ar-max` の指定によらない)。角部 (2 つの壁が会う所) と通常の壁を分けて数える。
- **閾値は二つに分ける**:
  - 構成の判定 (この層が壁法線に沿って作られているか): 半径方向の配置 (傾斜 40° で φ = 40°) と壁法線の配置 (φ ≈ 0) を分ける値。
    暫定値は試験で決め、精度の限界としては扱わない。
  - 精度の予算 (どこまで傾いてよいか): 壁応力・熱流束・δ\* の誤差から決める。今は根拠が無い (diagnose: node は LSQ 勾配と辺の中点の再構成なので、
    primal の角度だけから 5° のような限界は導けない)。第 2 段の A/B (§6.4) で決める。
- AR ≤ 1000 の近壁で φ が大きいセルは、精度の予算が決まるまで WARN にする (数と位置と最大 φ を出す)。
- 今の `SOFT-PASS` (0.1 % 未満の外れ値で終了 0) は、関門では緩和の対象セルに当てない。
- `mesh2d`: `wall_normal_layer` を既定で有効にするのは第 2 段 (格子が変わるので、結果が変わる)。`origin/main` に無いので、
  `feature/nozzle-core-grid` の main への取り込みを待つか、`mesh2d.py` の該当部分だけを cherry-pick する (§5.1)。

### 4.5 段階起動

- 段 (`soft`・`mid`・`main`) はプロファイルの段として持ち、値は今の生成器のものを書き写す (limiter は当面 2 を維持。§1.2 の 0 と
  `procedures/solver-settings.md` の 2 で正本同士が食い違うので、第 2 段で決める)。段の長さ・mid の次数も第 2 段。
- 同一メッシュの段間の引き継ぎは `restart_field.py` に揃える (runner_wt・walldriven の `interp_field.py` を置き換える)。AGENTS.md の既存の規則の適用。

### 4.6 推奨の側の直し

- `recommended-settings.md` の古い記述 (§4 の「`runner.py` が正本」、§2 の `evaluate.wall_treatment_sst`、§3 の化学種の例) を直す。
- `pRef` は今回必須にしない (diagnose: 「動作静圧」の決め方の根拠が無い)。圧力流束・境界・hoop 源での差し引きの整合、一定圧の自由流の残差、
  代表ノズルの残差の床と目的量への影響を確かめる項目を第 2 段に置く。
- `rerun_conditions.py` の独自の推奨 (`recommended_stages`) はプロファイルを参照する形にする。

### 4.7 移行と進行中の作業

- 第 1 段は数値の既定を変えないので、進行中の比較 (case/45 の V-h370 など) は影響を受けない (生成器の出力は同じ。宣言が付くだけ)。
- 宣言の無い run・記録の無いメッシュの扱い (全セッションの `run_case.sh` の起動に効く) はユーザの判断 (§5.1 #2)。案: main に入れた日から
  2 週間は WARN、以後は FAIL。移行の期間でも `forbid` とメッシュの FAIL は止める。
- `run_case.sh` を含む merge は、その作業ツリーで run が走っていないときに行う。

### 4.8 `tooling-design-problem-campaign-recipe` との関係

あちらの recipe は「形状・メッシュ・段・δ\* 反復・停止条件・数値設定 (推奨設定のどの節に従うか)・検査する VERDICT」を持つ構想 (draft)。
本 plan のプロファイルはその「数値設定」の部分を先に実体化したものとして作り、あちらの `resolve_evaluation` が実装されたら recipe から参照する。

## 5. 実装ステップ

1. 文書: `recommended-settings.md` の古い記述 (§4.6) を直し、値の正本が YAML に移ることを冒頭に書く。`methods/design/overview.md` に生成器がプロファイルから組むことを書く。
2. `recipes/recommended_settings.yaml` (今の生成器の値の書き写し + 推奨との差) と `recipe_settings.py` (読み込み・継承・照合・宣言・文書の生成)。
3. `check_solver_config.py` の厳格化 (未知・誤配置・重複・merge を FAIL にするモード) とプロファイルとの照合。
4. 変換器の呼び出し側でメッシュの実効条件の記録 (`<mesh>.provenance.json`)。design の `run_converter` と手順書の変換コマンド。
5. `check_mesh_quality.py` に壁直交性・AR の緩和条件・再利用の鍵。
6. `launch_gate.py` と `run_case.sh` の呼び出し。
7. `config_builder.py` と 5 つの生成器の置き換え (値は不変、宣言を書く)。段の書き換えを YAML の位置の方式へ。`interp_field` → `restart_field`。
8. AGENTS.md・skill `forge-config`・`procedures/calculation-workflow.md` に宣言と関門を書く。調整板で知らせる。
9. 第 2 段: プロファイルごとの受入 (§6.4) を事前登録して回し、合格したものから既定を推奨へ切り替える。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~§4 の論点を上位に諮る~~ | 判断: 2026-10-14 codex diagnose ([記録](../../notes/reviews/2026-10-14-generator-defaults-launch-gate-diagnose.md)) — 「数値既定の一括変更の前に、関門を明示したプロファイル・段と変換済みメッシュの実効条件に束縛し、CFD なしの負例試験を先行」。全件採用 (§4.0〜§4.7、§6.1) | F |
| 2 | 移行の扱いのユーザ判断 | 宣言の無い run・記録の無いメッシュを、main に入れた後いつから FAIL にするか (§4.7) | ユーザ |
| 3 | codex plan 段レビュー | `codex_review.py <この plan> --stage plan`。Critical / Major の採否を §6.1 と本表へ | F |
| 4 | 文書の直し (§5 の 1) | 合格: `check_plans.py` PASS、`recommended-settings.md` の古い記述 3 件が直っている | O |
| 5 | プロファイルと照合・宣言 (§5 の 2・3) | 合格: §6.2 の設定の負例がすべて対象キーを名指しで FAIL、正例が PASS、継承の衝突が例外 | O |
| 6 | メッシュの実効条件の記録 (§5 の 4) | 合格: 記録の無いメッシュ・`nodeInletCornerWall` が変換時と違う run が FAIL | O |
| 7 | 壁直交性 (§5 の 5) | 合格: §6.2 のメッシュの例と、再利用の鍵の A/B (§6.3) | O |
| 8 | 関門 (§5 の 6) | 合格: 負例で forge が起動しない、正例で起動する、`FORGE_LAUNCH_GATE=warn` が `forbid` を落とさない | O |
| 9 | 生成器の置き換え (§5 の 7) | 合格: 全生成器 × 代表の問題で、**置き換えの前後で生成した config がバイト同一** (`interp_field` → `restart_field` の段間を除く)、関門を通る | O |
| 10 | `wall_normal_layer` を main 系へ | forge-25 の main 取り込みの予定を確かめ、待つか cherry-pick するか | O |
| 11 | 周知 | AGENTS.md・skill・調整板 | O |
| 12 | 第 2 段の受入の事前登録 (§6.4) | プロファイルごとに受入の run・判定区間・量・許容差を書いて上位に諮る。`pRef` の確認、limiter の段の値、壁直交性の精度の予算を含む | F |

## 6. 検証

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | 2026-10-14 | [2026-10-14-generator-defaults-launch-gate-diagnose.md](../../notes/reviews/2026-10-14-generator-defaults-launch-gate-diagnose.md) | C0/M8/m1 | 全件採用。生産レシピの根拠は適用条件つきで分ける (SERN を一括にしない)、CFL の範囲は WARN・宣言したプロファイルからの変更は FAIL、`pRef` の必須化は外す、limiter は当面 2、壁直交性は機械判定するが角度の限界は誤差の予算から、`run_case.sh` は全起動の保証ではない、推定を廃して宣言を必須、変換時の実効条件と再利用の鍵を束縛、起動の健全性と生産の受入を分ける (§4.0〜§4.7) |

### 6.2 単体 (CFD なし)

- **生成器の正例**: R・WE・WEn・WS・S に代表の問題 (ベル・axis-Mach Euler・axis-Mach NS・SERN 2D・SERN 3D) を与え、置き換えの前後で config がバイト同一、関門を通る。
- **設定の負例** (それぞれ該当キーを名指しで FAIL): 宣言したプロファイルと違う `nStepInner` (逸脱の宣言なし)、`bndFirstOrder`、
  dual-time に `implicitRelax`、未知のキー、誤った節のキー、YAML の重複キー、理由の無い逸脱、宣言の欠落、
  宣言した段と config の段の不一致 (soft の宣言で 2 次の config)。advisory の範囲外の CFL は WARN で終了 0。
- **メッシュの実効条件の負例**: 記録の無いメッシュ、変換時の `nodeInletCornerWall` が 0 で run の config が 1。
- **壁直交性**: 40° 傾斜壁の半径方向の格子 (AR 4000) が緩和を受けられず FAIL、`wall_normal_layer` の格子 (AR 4000) が PASS、
  AR 800 の半径方向の格子は WARN (φ の最大と位置を出す)、四角形と三角形の混在で三角形を対象外と出す。
- **文書の同期**: `recommended-settings.md` の生成区間と YAML の一致。
- **既存の生産入力の棚卸し**: case/45 (B0)・case/46 (SERN)・case/70 の直近の生産 run の config に宣言を付けて関門に通し、`differs_from_recommended` の一覧を記録する。

### 6.3 再利用の鍵の A/B (diagnose の判別試験、forge 0 step)

固定した傾斜 40°・AR 100 の合成メッシュで、A は対象の境界を `slip` として判定して結果を保存、B は**その境界の種類だけを `wall` に変えて**再判定する。
試験用の構成の判定の閾値は 5° に固定する (生産の閾値の根拠とは扱わない)。

- A が「壁直交性は対象外・品質 PASS」、B が「壁直交性 FAIL」→ 再利用の鍵が bcond を含むことの確認。
- B が A の結果を再利用して PASS → 鍵が不十分。直してから先へ進む。
- 記録する量: 対象の壁の数・最大 φ・再利用の有無・VERDICT。

### 6.4 第 2 段の受入 (プロファイルごと、事前登録してから回す)

- 起動の健全性と生産の受入を分ける。健全性: 全段の実効設定・有限性・正値性。受入: 同じ作用素の区間の収束 VERDICT
  (`check_convergence.py --segment`)、目的量の準定常 VERDICT (`check_quasisteady.py`)、事前登録した許容差での旧既定との比較。
- 対象・量・許容差・step の上限はプロファイルごとに §5.1 #12 で事前登録する。結果を見てから条件を変えない。

## 7. 影響範囲

- design: `runner.py`・`runner_wt.py`・`runner_axismach.py`・`runner_sern.py`・`runner_sern3d.py`・`runner_walldriven.py`・`opt/driver.py`・`feedback/euler_loop.py`・`meshing/mesh2d.py`、新規 `evaluate/config_builder.py`。design の試験のうち `nStepInner` 5 を固定しているもの (`run_main_line_tests.py`・`run_stop_rule_tests.py` など)。
- tools: `check_solver_config.py`・`check_mesh_quality.py`・`run_case.sh`・`rerun_conditions.py` (独自の推奨をレシピに寄せる)、新規 `recipe_settings.py`・`launch_gate.py`・`recipes/recommended_settings.yaml`。
- 他のセッション: main に入ると、全セッションの `run_case.sh` の起動に関門が掛かる。進行中の作業は逸脱の宣言が要る。`run_case.sh` を含む merge は、その作業ツリーで run が走っていないときに行う (実行中の編集で run が壊れる)。
- 文書: `recommended-settings.md`・`calculation-workflow.md`・`methods/design/overview.md`・AGENTS.md・skill `forge-config`。

## 8. 完了条件

- [ ] `recommended-settings.md`・`methods/design/overview.md` を更新済み
- [ ] §6 を満たす
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status: done`、§9 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-14` — codex diagnose を全件採用して §4〜§6 を改訂: 二段に分け (第 1 段は値を変えずに宣言と関門、第 2 段はプロファイルごとの受入で既定を切り替え)、推定を廃して宣言を必須、変換時の実効条件・再利用の鍵の束縛、`pRef` の必須化を外す、壁直交性は構成の判定と精度の予算を分ける。
- `2026-10-14` — 初稿 (ユーザ要求「推奨設定を、あほなaiでも実施してくれるように」「壁に直交しない壁面メッシュ」)。§3 は生成器と推奨設定の突き合わせの調査 (読み取りのみ)。§4 は起案で、§4.8 の論点を上位に諮る前。
