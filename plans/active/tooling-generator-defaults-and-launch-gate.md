# 生成器の既定を推奨レシピから組み、投入前の関門で止める (設定とメッシュ)

## メタ

- **area**: `tooling` (`design/forge_design/evaluate/runner*.py`、`design/forge_design/meshing/mesh2d.py`、`solver_density_cuda/tools/check_solver_config.py`・`check_mesh_quality.py`・`run_case.sh`)
- **status**: `draft` (2026-10-14 codex diagnose・plan 段レビュー反映済み。実装前)
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

- **やる** (§4.0 の段取り 1・2):
  1. 生産レシピごとのプロファイル (適用条件・段・根拠・推奨との差) を機械可読のファイル 1 つに置き、値の正本にする (§4.1)。
  2. 生成器の直書きをやめ、共通の組み立て関数から config を作り、プロファイルの宣言を書く (§4.2)。
  3. `run_case.sh` の起動前に、宣言したプロファイル・段と変換済みメッシュの実効条件に照らして設定とメッシュを照合し、FAIL なら forge を起動しない (§4.3)。
  4. メッシュ判定に壁直交性を足し、AR の緩和条件 (壁法線の構造層) を機械判定にする (§4.4)。
  5. 同一メッシュの段間の引き継ぎを `restart_field.py` に揃える (§4.5)。
  6. 生成器の既定を各機種の最新の生産レシピに揃える (§4.0 の段取り 2)。
- **やる** (§4.0 の段取り 3): 機種の推奨を §1 に寄せるキーは、機種ごとに受入の計算を事前登録して回してから切り替える。`mesh2d` の壁法線の近壁層の既定化もここ (§6.4)。
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

### 4.0 機種ごとの推奨設定 (ユーザ決定 2026-10-14、codex diagnose を反映)

**ユーザ決定 (2026-10-14、原文)**: 「機種ごとに推奨設定を変えたらいいだけでは？」(二段で進める案への回答)。

- 推奨設定は**機種と段ごと**に持つ (風洞軸対称 axis-Mach の Euler・NS、ベル、デュアルベル、SERN 2D 有限ベース、SERN 3D、walldriven など。§4.1 の「プロファイル」はこの 1 セットを指す)。
  値は**その機種の最新の生産レシピ** (検証済みの run がある値) とし、根拠の plan・run を付ける。`recommended-settings.md` §1 は機種に依らない共通部分と、
  機種の決まっていない手書きの run の既定にする。
- **生成器の既定 = その機種の推奨設定**にする。問題 YAML に書かなくても最新の生産レシピで回るようにする。今は axis-Mach NS の既定が CFL 1.0・緩和なしで、
  生産の B0 は問題 YAML に `main_line` を明示して初めて入る。判断力の低いエージェントが新しい問題 YAML を書くと古い既定で回るのはこの形である。
- 段取り (diagnose「数値既定を一括で変える前に関門を固める」を採用):
  1. **関門と宣言**: 機種ごとの推奨設定を今の生成器の値で書き写し、宣言・共通の組み立て関数・関門・壁直交性の判定を入れる。この時点では生成器の出力は変えない。
  2. **生成器の既定を機種の生産レシピに揃える**: 生成器の今の既定と機種の生産レシピが違うキーは、生産レシピ側に揃える。根拠の run があるので新たな受入の計算は要らない
     (ただし揃えた後の生成器の出力が、根拠の run の実効設定と一致することを確かめる)。どの run を根拠にするかは機種ごとに表にして上位に確かめる (§5.1)。
  3. **§1 との差**: 機種の生産レシピに根拠が無く、生成器の値が §1 と違うキー (例: ベルの `turbulentPrandtl` 0.85、`nodeInletCornerWall` の欠落) は、
     今の値をその機種の推奨として残し、§1 との差を記録する。§1 に揃えるときだけ、機種ごとに受入の計算を事前登録して回す (§6.4)。
     起動の健全性 (NaN が無い) だけでは承認しない (SERN の記録に、30 step と 4000 step で順位が入れ替わった例がある)。
- diagnose の指摘「§1 の値を検証済みの生産設定に当てると壊す」(SERN の CFL 4 は 2D 有限ベースで発散) は、機種ごとに持つことで起きない。

### 4.1 値の正本: プロファイルの機械可読ファイル

- `solver_density_cuda/tools/recipes/recommended_settings.yaml` を置く (design からも tools からも読める場所)。読み込みと照合は `solver_density_cuda/tools/recipe_settings.py`。
- プロファイル (= 機種と段ごとの推奨設定) は**適用条件**を持つ: 次元 (2D 平面・2D 軸対称・3D)、物理 (Euler / 層流 / SST、CPG / TP)、壁の扱い (slip / no-slip 断熱 / 等温)、
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
- 値の出所を 3 つに分けて記録する: プロファイル、問題 YAML の物理条件 (気体・壁温・入口条件)、明示の逸脱 (キー・値・理由、任意で根拠)。**理由の無い逸脱は例外。根拠 (plan の節・run パス) は任意** (ユーザ決定 2026-10-14「試し計算は手軽にやりたいですね」: 根拠の必須化は試し計算の手間を増やすので採らない)。
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
     **ただし FAIL にする前に受理キーの一覧を直す** (plan 段 M1): 今の一覧は C++ の正規表現抽出で不完全で、ソルバが受理する
     `conjugate.mode`・`conjugate.gate` を誤検出する (`check_solver_config.py:77`、`solverConfig.cpp:541`・`:558`)。一覧をソルバの実装と突き合わせて直し、
     CHT・TP・凝縮・周期・非定常の正当な config を正例にした試験で誤拒否が無いことを確かめてから FAIL にする。重複キー・merge key の検査は
     既存の `yaml_strict.py` を共用する (新規に書かない)。
  2. **メッシュの実効条件**: 変換時にだけ効く設定がある (`nodeInletCornerWall` は変換時に焼き込まれ、起動時の値は効かない: `procedures/solver-settings.md` の該当節)。
     変換時の実効設定を `<mesh>.provenance.json` に書き、関門は run の config・bcond と照合する。記録の項目 (plan 段 M2 で追加):
     `discretization`・`nodeInletCornerWall`・`wallDistExtraPhysIDs` (`wall_dist` を変える)・**全境界**の種類と physID
     (壁だけでなく `inlet_*` も半割面の所属を決める: `gmshReader.hpp:1601`・`:1625`)・変換器の sha256。
     記録は**変換成果物の不変部分のハッシュ** (座標・接続・境界の所属・`wall_dist`) に結び付ける。段間の restart は `/VALUE` を書き換えるので
     ファイル全体の sha は使わず、初期場のハッシュは別に記録する。記録の無いメッシュは FAIL (移行の扱いは §4.7)。
  3. **メッシュの品質**: §4.4 の判定。結果の再利用は (メッシュの不変部分のハッシュ、**physID と壁の面の集合の対応**、閾値、検査器の版) が全部一致したときだけ (plan 段 M6)。
  4. **バイナリと省略値** (plan 段 M3、一部採用): 省略したキーの値は実行するバイナリで決まる (例 `scalarGradient` は node で省略時 `lsq`、
     `solverConfig.cpp:249`)。`run_case.sh` は `FORGE_BIN` で別ビルドも回す。そこで (a) プロファイルが値に依存するキーは**省略せず書く**
     (省略に頼るキーを作らない)、(b) `LAUNCH_GATE.txt` にバイナリの sha256 とプロファイルの内容のハッシュを記録する、
     (c) forge の起動時の実効値の出力 (`'... effective'` の行) がある量は起動後に照合する。**対応の分からないバイナリの拒否は採らない**
     (理由: 旧ビルド・倍精度ビルドとの A/B は日常の検証で使う。拒否すると検証の経路を塞ぐ。代わりに記録と WARN)。
  5. 結果を `LAUNCH_GATE.txt` (VERDICT・プロファイル・段・逸脱・WARN・メッシュのハッシュ・バイナリ) に書き、`RUN_PROVENANCE.txt` にも 1 行写す。
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
    primal の角度だけから 5° のような限界は導けない)。段取り 3 の A/B (§6.4) で決める。
- AR ≤ 1000 の近壁で φ が大きいセルは、精度の予算が決まるまで WARN にする (数と位置と最大 φ を出す)。
- 今の `SOFT-PASS` (0.1 % 未満の外れ値で終了 0) は、関門では緩和の対象セルに当てない。
- `mesh2d`: `wall_normal_layer` を既定で有効にするのは段取り 3 (格子が変わるので結果が変わる。受入の計算を回してから)。`origin/main` に無いので、
  `feature/nozzle-core-grid` の main への取り込みを待つか、`mesh2d.py` の該当部分だけを cherry-pick する (§5.1)。

### 4.5 段階起動

- 段 (`soft`・`mid`・`main`) はプロファイルの段として持ち、値は今の生成器のものを書き写す (limiter は当面 2 を維持。§1.2 の 0 と
  `procedures/solver-settings.md` の 2 で正本同士が食い違うので、段取り 2・3 で機種ごとに決める)。段の長さ・mid の次数も同じ。
- 同一メッシュの段間の引き継ぎは `restart_field.py` に揃える (runner_wt・walldriven の `interp_field.py` を置き換える)。AGENTS.md の既存の規則の適用。
  初期場が変わる修正なので、config の一致ではなく restart の検査 (保存量・化学種の属性・`wall_dist` の保持、欠けた保存量・Inf・負の密度なら次の段を起動しない) を負例つきで置く (plan 段 M5)。
- **実在する起動の経路をすべて棚卸しする** (plan 段 M5): `soft`・`mid`・`main` のほかに、SERN の層流暖機・暖機 ramp・作動点を変えた後の適応段、
  axis-Mach の本段の次数のままの ramp がある。経路ごとに「段のプロファイル・最終 config・`RECIPE.json`・履歴の退避・`StageManifest`」を一組で定義する。
  SERN・WT・walldriven は今 `StageManifest` を書いていないので、書かせないと段取り 3 の `check_convergence.py --segment` が成立しない。

### 4.6 推奨の側の直し

- `recommended-settings.md` の古い記述 (§4 の「`runner.py` が正本」、§2 の `evaluate.wall_treatment_sst`、§3 の化学種の例) を直す。
- `pRef` は今回必須にしない (diagnose: 「動作静圧」の決め方の根拠が無い)。圧力流束・境界・hoop 源での差し引きの整合、一定圧の自由流の残差、
  代表ノズルの残差の床と目的量への影響を確かめる項目を段取り 3 に置く。
- `rerun_conditions.py` の独自の推奨 (`recommended_stages`) はプロファイルを参照する形にする。

### 4.7 移行と進行中の作業

- 段取り 1 は生成器の出力を変えない (宣言が付くだけ)。段取り 2 で既定が変わる機種の進行中の比較 (case/45 の V-h370 など) は、問題 YAML に値を明示してあるので影響を受けない (明示した値は既定より優先)。揃えた後に、進行中の問題 YAML の生成結果が変わらないことを試験で確かめる。
- 宣言の無い run・記録の無いメッシュの扱い (全セッションの `run_case.sh` の起動に効く): **ユーザ決定 2026-10-14「aで」** = 棚卸しの完了の後は、宣言の無い手書きの run も止める (それまでは WARN。設計チェーンの run は生成器が宣言を書く)。
  plan 段 M4 を採用し、**強制の開始は日付でなく、対象の経路の棚卸しが終わったこと**を条件にする: 標準検証のケース (`procedures/verification/README.md` の衝撃管・
  KEEP 周期箱・凝縮など) と設計チェーンの全経路にプロファイルか宣言の方法があること。手書きの実験の run には「実験」の宣言 (推奨への適合ではなく、外すキーと理由を書く宣言。**理由は 1 行の自由記述で足り、plan や比較 run の指定は要らない**: ユーザ決定 2026-10-14) を用意し、推奨に従う run と区別する。
- 旧いメッシュは、元の入力と変換の条件を確かめられるときだけ記録を作って移す。確かめられなければ再変換を求める。
  今の config から過去の変換条件を推定して記録を埋めることはしない。
- 移行の期間でも `forbid` とメッシュの FAIL は止める。
- `run_case.sh` を含む merge は、その作業ツリーで run が走っていないときに行う。

### 4.8 `tooling-design-problem-campaign-recipe` との関係

あちらの recipe は「形状・メッシュ・段・δ\* 反復・停止条件・数値設定 (推奨設定のどの節に従うか)・検査する VERDICT」を持つ構想 (draft)。
本 plan のプロファイルはその「数値設定」の部分を先に実体化したものとして作り、あちらの `resolve_evaluation` が実装されたら recipe から参照する。

## 5. 実装ステップ

1. 文書: `recommended-settings.md` を「機種ごとの推奨設定 (値は YAML が正本) + §1 共通」の形に直し、古い記述 (§4.6) を直す。`methods/design/overview.md` に生成器が機種の推奨設定から組むことを書く。
2. 受理キーの一覧を直す (`check_solver_config.py`、plan 段 M1)。正当な config の正例で誤拒否が無いこと。
3. `recipes/recommended_settings.yaml` (機種と段ごとの推奨設定。値は今の生成器の書き写し + §1 との差 + 根拠) と `recipe_settings.py` (読み込み・継承・照合・宣言・文書の生成)。
4. 変換時の実効条件の記録 (`<mesh>.provenance.json`、不変部分のハッシュ、plan 段 M2)。design の `run_converter` と手順書の変換コマンド。
5. `check_mesh_quality.py` に壁直交性・AR の緩和条件・再利用の鍵 (physID と壁の面の集合)。
6. 起動の経路の棚卸し (全段・`StageManifest`・restart の検査、plan 段 M5)。
7. `launch_gate.py` と `run_case.sh` の呼び出し (バイナリと省略値の記録、plan 段 M3)。
8. `config_builder.py` と 5 つの生成器の置き換え (出力は不変、宣言を書く)。段の書き換えを YAML の位置の方式へ (`yaml_strict.py` 共用)。`interp_field` → `restart_field`。
9. 段取り 2: 機種ごとに「生成器の既定 ↔ 最新の生産レシピ (根拠の run)」の表を作り、上位に確かめてから既定を揃える。
10. AGENTS.md・skill `forge-config`・`procedures/calculation-workflow.md` に宣言と関門を書く。調整板で知らせる。
11. 段取り 3: §1 に寄せる機種・キーごとに受入の計算 (§6.4) を事前登録して回す。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~§4 の論点を上位に諮る~~ | 判断: 2026-10-14 codex diagnose ([記録](../../notes/reviews/2026-10-14-generator-defaults-launch-gate-diagnose.md)) — 「数値既定の一括変更の前に、関門を明示したプロファイル・段と変換済みメッシュの実効条件に束縛し、CFD なしの負例試験を先行」。全件採用 (§4.0〜§4.7、§6.1) | F |
| 2 | ~~codex plan 段レビュー~~ | 判断: 2026-10-14 ([記録](../../notes/reviews/2026-10-14-tooling-generator-defaults-and-launch-gate-plan.md)、GO-with-changes C0/M6/m1) — M1・M2・M4・M5・M6・m7 採用、M3 は一部採用 (対応の分からないバイナリの拒否は却下: 旧・倍精度ビルドとの A/B を塞ぐため。記録と WARN にする)。§4.3〜§4.7・§6 に反映 | F |
| 3 | ~~強制の開始のユーザ判断~~ | 宣言の無い手書きの run・記録の無いメッシュを、棚卸しの完了後に止めるか (§4.7)。設計チェーンの run は生成器が宣言を書くので最初から止めてよいか。**2026-10-14 の対話**: 案 a (棚卸しの完了後に手書きの run も止める) の説明に対しユーザは「推奨から変えることはアドリブで実施可能」を確認し、「試し計算は手軽にやりたい」と決定 (逸脱は理由だけで通す、根拠は任意: §4.2・§4.7)。**判断: 2026-10-14 ユーザ決定「aで」** — 棚卸しの完了 (§4.7 の条件) の後は、宣言の無い手書きの run も `run_case.sh` で止める。それまでは WARN | ユーザ |
| 4 | 受理キーの一覧の修正 (§5 の 2) | 合格: CHT (`conjugate.mode`・`conjugate.gate`)・TP・凝縮・周期・非定常の正当な config が誤検出されない、未知・誤配置のキーの負例が名指しで出る。**実装 (2026-10-14、implementer、未 commit)**: 取りこぼしの原因は `config_key_inventory.py` の手書きの変数→節の表 `VAR2PATH` (`cj`・`g` が無い) と、既定値に括弧を含むヘルパ呼び出し (`std::string("local1d")`) で第 4 引数 (節名) を落とす正規表現。束縛 (`auto cj = config["conjugate"];`) をソースから導き、引数を括弧の対応で分ける方式に替えた (`extract_keys_with_report`)。一覧 235 → 236 パス (`cj.mode/flux/back` → `conjugate.*`、`conjugate.gate` を追加、`space.limiterMatchRecon` を拒否へ)。静的に入れられない書き方: 列の要素への添字 7 (値 `physProp.species` の中身)・`mode: fem2d` のときだけの拒否 5 キー (`solverConfig.cpp:593`)。実在の正例 8 本と負例を `tests/unit/test_check_solver_config.py` に追加し ALL PASS。WARN → FAIL は #10 | O |
| 5 | 文書の直し (§5 の 1) | 合格: `check_plans.py` PASS、`recommended-settings.md` の古い記述 3 件が直り、機種ごとの節の形になっている | O |
| 6 | 機種ごとの推奨設定と照合・宣言 (§5 の 3) | 合格: §6.2 の設定の負例がすべて対象キーを名指しで FAIL、正例が PASS、継承の衝突が例外、「実験」の宣言は理由が無いと書けない | O |
| 7 | メッシュの実効条件の記録 (§5 の 4) | 合格: §6.2 のメッシュの来歴の負例 (記録なし・記録の取り違え・入口の分類だけの変更・`wallDistExtraPhysIDs` だけの変更・`nodeInletCornerWall` の不一致) が FAIL、保存量だけを更新した正常な restart は PASS | O |
| 8 | 壁直交性 (§5 の 5) | 合格: §6.2 のメッシュの例と §6.3 の再利用の鍵の試験。**実装 (2026-10-14、implementer、未 commit)**: `check_mesh_quality.py --bcond` (横断辺の φ・壁からの層の行進・AR 緩和の条件・角部の別集計・`<h5>.quality.json` の再利用)。暫定 `--phi-construct 5°`・`--theta-tol 5°` (構成の判定): case/45 の生産の格子 (`..._tw370_cgref_h_b0.yaml`、float64) で層 1 の φ 最大 0.06°、AR > 1000 の層 1..19 で φ・内角のずれ最大 0.07°、同じ壁の半径方向の線で 26.6°、合成の 40° 傾斜で 40.0°。判定は格納された座標で行い、float32 の丸めで条件を外れる層は「丸めの範囲」として緩和しない (同じ生産の格子の float32 版は FAIL)。`tests/unit/test_check_mesh_quality_wall.py` (§6.2・§6.3・3D 六面体・角部・cell 変換の経路) ALL PASS、`--bcond` なしの VERDICT・終了コードは 18/18 で従来と同じ | O |
| 8b | 壁直交性の判定で出た判断 (#8 の実装から) | (1) 判定は格納された座標で行う (丸めは理由の区別だけ)。float32 の h5 の薄い層は φ が丸めで数度〜十数度揺れる: case/45 の生産の格子は float64 で φ 0.07° だが、float32 で格納すると φ 11.96°・AR 5609 で FAIL。今 `--ar-max 5000` で通っている float32 のメッシュが関門で落ちうる (対策は倍精度の幾何を書く変換器 = M6 のブランチの cb6a0d8b、1 回目の統合で main に入る)。この方式でよいか。(2) AR 緩和の層の条件を壁の physID ごとの帯全体で判定する (§4.4 の字義どおり。1 か所の傾きでその壁の全セルの緩和が消える) か、列ごとにするか。(3) 角部の扱い (physID が 2 つ会う節点か法線の折れ ≥ 45° を角部とし、帯の条件から外して個別に判定) でよいか。(4) `--bcond` なしの `--ar-max 5000` は旧い判定のまま (回帰を守るため)。関門は必ず `--bcond` を渡す (#10)。暫定の閾値 φ 5°・θ_tol 5° は構成の判定の値 (生産の格子の 0.07° の約 70 倍、40° の半径方向の線を拒む) で精度の限界ではない | F |
| 9 | 起動の経路の棚卸し (§5 の 6) | 合格: 全経路の段の一覧 (段のプロファイル・最終 config・`RECIPE.json`・履歴・`StageManifest`) が表になり、SERN・WT・walldriven が `StageManifest` を書く。restart の負例 (欠けた保存量・Inf・負の密度) で次の段を起動しない | O |
| 10 | 関門 (§5 の 7) | 合格: 負例で forge が起動しない、正例で起動する、`FORGE_LAUNCH_GATE=warn` が `forbid` を落とさない、バイナリの sha とプロファイルのハッシュが記録される | O |
| 11 | 生成器の置き換え (§5 の 8) | 合格: 同一バイナリで、全生成器 × 代表の問題の**全段の起動列**の config が置き換えの前後でバイト同一 (`interp_field` → `restart_field` の段間を除く)、関門を通る | O |
| 12 | 機種ごとの「生成器の既定 ↔ 生産レシピ」の表 (§5 の 9) | 機種ごとに根拠の run (パス・区間・収束 VERDICT・準定常 VERDICT・保存データの所在、plan 段 m7) を表にし、揃えるキーを決める。上位に確かめてから既定を揃える | F |
| 13 | ~~`wall_normal_layer` を main 系へ~~ | 済み (2026-10-14): PR #8 (81e92491) で main に入った (c096d66c、M6 のブランチ経由)。この作業ツリーには PR #9 の取り込み (7187cf87) で入っている | O |
| 14 | 周知 | AGENTS.md・skill・調整板 | O |
| 15 | §1 に寄せる受入の事前登録 (§6.4) | 機種・キーごとに受入の run・判定区間・量・許容差を書いて上位に諮る。`pRef` の確認、limiter の段の値、壁直交性の精度の予算 (`check_wall_resolution.py` と局所の τ_w・q_w の比較を含む) | F |

## 6. 検証

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-14 | [2026-10-14-tooling-generator-defaults-and-launch-gate-plan.md](../../notes/reviews/2026-10-14-tooling-generator-defaults-and-launch-gate-plan.md) | GO-with-changes, C0/M6/m1 | M1 採用→§4.3・#4 (受理キーの一覧を先に直す、`yaml_strict.py` 共用)、M2 採用→§4.3・#7 (変換来歴の項目と不変部分のハッシュ)、M3 一部採用→§4.3 (省略に頼らない・バイナリの記録と起動後の照合。対応の分からないバイナリの拒否は却下: 旧・倍精度ビルドとの A/B を塞ぐ)、M4 採用→§4.7・#3 (強制の開始は棚卸しの完了、「実験」の宣言、旧メッシュは再変換)、M5 採用→§4.5・#9・#11 (全起動経路の段と `StageManifest`、restart の検査)、M6 採用→§6.3 (チェッカー単体、WARN と再利用なし、physID の入れ替えの負例)、m7 採用→#12・#15 (根拠の run の VERDICT と所在を事前登録) |
| diagnose | 2026-10-14 | [2026-10-14-generator-defaults-launch-gate-diagnose.md](../../notes/reviews/2026-10-14-generator-defaults-launch-gate-diagnose.md) | C0/M8/m1 | 全件採用。生産レシピの根拠は適用条件つきで分ける (SERN を一括にしない)、CFL の範囲は WARN・宣言したプロファイルからの変更は FAIL、`pRef` の必須化は外す、limiter は当面 2、壁直交性は機械判定するが角度の限界は誤差の予算から、`run_case.sh` は全起動の保証ではない、推定を廃して宣言を必須、変換時の実効条件と再利用の鍵を束縛、起動の健全性と生産の受入を分ける (§4.0〜§4.7) |

### 6.2 単体 (CFD なし)

- **生成器の正例**: R・WE・WEn・WS・S に代表の問題 (ベル・axis-Mach Euler・axis-Mach NS・SERN 2D・SERN 3D) を与え、**全段の起動列**の config が置き換えの前後でバイト同一 (同一バイナリ)、関門を通る。
- **restart の負例**: 欠けた保存量・Inf・負の密度の場では次の段を起動しない。
- **設定の負例** (それぞれ該当キーを名指しで FAIL): 宣言したプロファイルと違う `nStepInner` (逸脱の宣言なし)、`bndFirstOrder`、
  dual-time に `implicitRelax`、未知のキー、誤った節のキー、YAML の重複キー、理由の無い逸脱、宣言の欠落、
  宣言した段と config の段の不一致 (soft の宣言で 2 次の config)。advisory の範囲外の CFL は WARN で終了 0。
- **受理キーの正例**: CHT (`conjugate.mode`・`conjugate.gate`)・TP・凝縮・周期・非定常の正当な config が未知・誤配置と判定されない。
- **メッシュの実効条件の負例**: 記録の無いメッシュ、記録の取り違え (別メッシュの記録)、入口の分類だけの変更、`wallDistExtraPhysIDs` だけの変更、
  変換時の `nodeInletCornerWall` が 0 で run の config が 1。**正例**: 保存量だけを更新した restart の後のメッシュは通る。
- **壁直交性**: 40° 傾斜壁の半径方向の格子 (AR 4000) が緩和を受けられず FAIL、`wall_normal_layer` の格子 (AR 4000) が PASS、
  AR 800 の半径方向の格子は WARN (φ の最大と位置を出す)、四角形と三角形の混在で三角形を対象外と出す。
- **文書の同期**: `recommended-settings.md` の生成区間と YAML の一致。
- **既存の生産入力の棚卸し**: case/45 (B0)・case/46 (SERN)・case/70 の直近の生産 run の config に宣言を付けて関門に通し、`differs_from_recommended` の一覧を記録する。

### 6.3 再利用の鍵の試験 (品質チェッカー単体、forge 0 step。diagnose の判別試験を plan 段 M6 で直した形)

固定した傾斜 40°・AR 100 の合成メッシュ (既存の `cell_metrics()` で AR 100・skew 0.444・φ 40°)。

- A: 対象の境界を `slip` として判定して結果を保存 → 期待: 対象の壁 0、品質 PASS。
- B: **その境界の種類だけを `wall` に変えて**再判定 → 期待: 対象の壁あり・φ 40° で **WARN** (AR ≤ 1000 なので FAIL ではない)、**A の結果を再利用しない**。
- C: 種類の集合は同じまま、壁にする physID を入れ替える → 期待: 再利用しない。
- B か C が A の結果を再利用したら鍵が不十分。直してから先へ進む。来歴の不一致の拒否は §6.2 の別の試験で見る。
- 記録する量: 対象の壁の数・最大 φ・再利用の有無・VERDICT。試験用の構成の判定の閾値は 5° (生産の閾値の根拠とは扱わない)。

### 6.4 §1 に寄せる受入 (機種・キーごと、事前登録してから回す)

- 起動の健全性と生産の受入を分ける。健全性: 全段の実効設定・有限性・正値性。受入: 同じ作用素の区間の収束 VERDICT
  (`check_convergence.py --segment`)、目的量の準定常 VERDICT (`check_quasisteady.py`)、事前登録した許容差での旧既定との比較。
- 対象・量・許容差・step の上限・基準の run (区間・両 VERDICT・保存データの所在) は機種・キーごとに §5.1 #15 で事前登録する。壁直交性の精度の予算には `check_wall_resolution.py` と局所の τ_w・q_w の比較を含める。結果を見てから条件を変えない。

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

- `2026-10-14` — main (PR #8・#9) を取り込んだ (7187cf87)。#13 は PR #8 で済み。ユーザ指示「④も進めて」で #9 → #6 → #5 → #7 → #10 → #11 の順に着手、#8b・#12 は上位に諮る。
- `2026-10-14` — ユーザ決定「aで」: 棚卸しの完了の後は宣言の無い手書きの run も止める (§4.7・§5.1 #3)。
- `2026-10-14` — ユーザ決定「試し計算は手軽にやりたいですね」: 逸脱・「実験」の宣言は理由 (1 行の自由記述) だけで通し、根拠 (plan の節・run パス) は任意にする (§4.2・§4.7)。
- `2026-10-14` — ユーザ決定「機種ごとに推奨設定を変えたらいいだけでは？」を §4.0 に反映 (推奨設定は機種と段ごと、値は最新の生産レシピ、生成器の既定 = 機種の推奨)。codex plan 段レビュー (GO-with-changes C0/M6/m1) の採否を §4.3〜§4.7・§5.1・§6 に反映 (M3 は一部却下)。
- `2026-10-14` — codex diagnose を全件採用して §4〜§6 を改訂: 二段に分け (第 1 段は値を変えずに宣言と関門、第 2 段はプロファイルごとの受入で既定を切り替え)、推定を廃して宣言を必須、変換時の実効条件・再利用の鍵の束縛、`pRef` の必須化を外す、壁直交性は構成の判定と精度の予算を分ける。
- `2026-10-14` — 初稿 (ユーザ要求「推奨設定を、あほなaiでも実施してくれるように」「壁に直交しない壁面メッシュ」)。§3 は生成器と推奨設定の突き合わせの調査 (読み取りのみ)。§4 は起案で、§4.8 の論点を上位に諮る前。
