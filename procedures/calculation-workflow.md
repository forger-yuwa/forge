# Forge Calculation Workflow

この文書は、このリポジトリで計算ケースを準備して `forge` を実行する際の標準手順をまとめたものです。計算方法を案内するときは、まずこの文書の流れを基準にし、ケース側に専用スクリプトや README がある場合はそちらを優先します。

## 対象範囲

`case/` 配下の多くのケースは、概ね次の構成を取ります。

- `mesh/`: `.geo` や生成スクリプトなど、メッシュ作成用の入力
- `run_*`: `run_slau`、`run_roe`、`run_keep` など、ソルバー実行用ディレクトリ

質問が特定ケースに限定されていない場合も、基本的にはこの構成を前提に案内します。

## 標準手順

計算を実行するときは、既存の `run_*` ディレクトリをそのまま使い回さず、必ず新しい `run_*` ディレクトリを作ってから実行する。

理由:

- 既存の結果やログを上書きしないため
- 実行時の設定差分を追跡しやすくするため
- `residual_history.csv`、`res_*.h5`、実行ログなどを run ごとに保存するため
- 実行後に residual plot PNG を run ごとに残すため

## run ディレクトリ命名規則

`run_*` が増えたときに新旧や派生関係が分からなくならないよう、run 名は次の形式で統一する。

```text
run_<NNNN>_<scheme>_<short-slug>
```

例:

```text
run_0001_slau_baseline
run_0002_slau_lowbp_1st
run_0003_slau_lowbp_muscl_restart
run_0004_roe_baseline
```

ルール:

- `<NNNN>` は `0000` から `9999` までの 4 桁ゼロ埋め連番とする。
- `<scheme>` は `slau`, `roe`, `keep`, `ausm`, `ausmUP` などスキーム識別子を入れる。
- 連番はケースディレクトリごとに単調増加で使い、スキームを変えても振り直さない。
- 欠番は気にせず、既存最大番号の次を使う。
- `<short-slug>` は 2〜4 語程度の短い ASCII ラベルにし、変更意図を一目で分かるようにする。
- 日付は原則 run 名に入れず、時系列は連番で管理する。
- 既存 run を複製して派生 run を作るときも、新しい連番を必ず振り直す。

推奨運用:

- まず `ls` や `find` でケース内の既存 `run_` を確認する。
- その中の最大番号に `+1` した名前で新しい run を作る。
- `short-slug` は `baseline`, `lowbp`, `lowbp_1st`, `muscl_restart`, `axisym_m2` のように目的中心で付ける。

非推奨:

- `run_slau_test`, `run_new_slau`, `run_latest_slau` のような意味が時間で変わる名前
- `run_20260531_slau_xxx` のように日付だけで順序管理する名前
- `run_0007_slau_copy`, `run_0008_slau_final2` のような派生関係が読めない名前

1. 対象のケースディレクトリを `case/` 配下から特定する。
2. 基準にする既存の `run_*` を複製して、新しい `run_*` ディレクトリを作る。
3. 必要なら新しい `run_*` 内の `solverConfig.yaml` などを調整する。
4. `mesh/` に入り、Gmsh で `.geo` から `.msh` を生成する。
5. 生成した `.msh` を新しい `run_*` ディレクトリへコピーする。
6. その新しい `run_*` ディレクトリで `convertGmshToForge` を実行し、HDF5 に変換する。
7. 実行する `forge` バイナリが最新ソースで再ビルドされているか確認する (下記「実行前のバイナリ鮮度チェック」)。
8. 同じ新しい `run_*` ディレクトリで `forge` を実行する。
9. 実行後、その新しい `run_*` ディレクトリで `residual_history.csv` から `residual_history.png` を生成する。

重要なのは、`solverConfig.yaml` がある新しい `run_*` ディレクトリを起点に変換と実行を行うことです。

## 実行前のバイナリ鮮度チェック (重要)

`forge` を実行するコマンドは `solver_density_cuda/build/forge` を直接呼ぶだけで、**ソースの再ビルドを自動では行わない**。`cuda_forge/*.cu` などを編集した後に古い `build/forge` のまま計算すると、現行ソースと無関係な **stale バイナリの結果が silently に出る** (例: ソース側で実装した乱流モデルが反映されず、結果が壊れて見える)。

`forge` を実行する前に、必ずバイナリがソースより新しいことを確認すること。

```bash
cd /home/sano/work/forge/solver_density_cuda
# build/forge より新しいソースがあれば stale (= 再ビルドが必要)
find cuda_forge *.cpp *.cu *.hpp -newer build/forge 2>/dev/null
```

上のコマンドが何か出力したら、計算前に再ビルドする。Docker dev image 内で `tools/build.sh` を使う:

```bash
cd /home/sano/work/forge
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/workspace" \
  -e HDF5_INC=/usr/include/hdf5/serial \
  -e HDF5_LIBDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial \
  forge-solver:cuda-dev \
  bash -c "cd /workspace/solver_density_cuda && bash tools/build.sh"
```

設定変更 (例: `convMethod`) が解に与える影響を診断するときは、まず stale でないバイナリで「変更前設定」を再現する control run を取り、スキーム差とコード世代の差を取り違えないこと。

## 標準コマンド例

```bash
# 既存 run を複製して新しい run を作る
cp -r run_0003_slau_baseline run_0004_slau_implicit_check

# <case>/mesh で実行
gmsh -3 input.geo -o case.msh -format msh4

# 生成メッシュを run ディレクトリへコピー
cp case.msh ../run_0004_slau_implicit_check/

# <case>/run_0004_slau_implicit_check で実行
convertGmshToForge case.msh case.h5
# forge 実行前に build/forge が stale でないか確認 (「実行前のバイナリ鮮度チェック」参照)
forge
python3 /home/sano/work/forge/solver_density_cuda/tools/plot_implicit_residuals.py \
	--input residual_history.csv \
	--output residual_history.png
```

Docker 経由で上記を手動実行する場合は、必ず `--user "$(id -u):$(id -g)"` を付けること (詳細は `procedures/development-environment.md` のファイル所有者ルールを参照)。

読み替えポイント:

- `input.geo`: 実際のジオメトリ入力ファイル名
- `case.msh`: 生成するメッシュファイル名
- `run_0004_slau_implicit_check`: 実際に使う新しい `run_*` ディレクトリ名
- `case.h5`: `solverConfig.yaml` の `mesh.meshFileName` と `mesh.valueFileName` に一致する HDF5 名

## Docker を使う場合の流れ

`solver_density_cuda` の Docker 環境を使う場合も、計算の論理的な流れは同じです。

1. 必要なら Docker イメージをビルドする。
2. 基準にする既存の `run_*` を複製して、新しい `run_*` を作る。
3. ケースの `mesh/` で `.msh` を生成する。
4. `.msh` を新しい `run_*` にコピーする。
5. その新しい `run_*` ディレクトリで `convertGmshToForge` を実行する。
6. 同じ新しい `run_*` ディレクトリで `forge` を実行する。
7. 実行後、その新しい `run_*` ディレクトリで `residual_history.csv` から `residual_history.png` を生成する。

代表例:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/docker_build.sh dev      # forge-solver:cuda-dev

cd /home/sano/work/forge/case/20.naca_ml/001.test
cp -r run_slau run_slau_20260510_implicit_check

cd /home/sano/work/forge/case/20.naca_ml/001.test/run_slau_20260510_implicit_check
./run_case.sh
```

ケースに `run_case.sh` のような専用ヘルパーがある場合は、手順を手で再構成するより、そのヘルパーを優先します。ただし、その場合も既存の `run_*` を直接使わず、複製した新しい `run_*` 側で実行する。

Residual plot は PNG で残せば十分とし、既定では `residual_history.csv` を入力に `residual_history.png` を生成する。

## Gmsh / ParaView の起動

Docker 環境経由で Gmsh や ParaView を開くときは、`solver_density_cuda/tools/` 配下のホスト側ラッパースクリプトを使います。

Gmsh:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_gmsh_gui.sh
```

ParaView:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_paraview_gui.sh
```

ファイルを直接開く例:

```bash
./tools/run_gmsh_gui.sh /home/sano/work/forge/case/08.bump/mesh/bump.geo
./tools/run_paraview_gui.sh /home/sano/work/forge/case/08.bump/run_slau/res_0.xmf
```

補足:

- これらのスクリプトはコンテナ内ではなくホスト側から実行する。
- WSL2 では WSLg を前提に GUI 転送する。
- ParaView のアイコン表示が崩れる場合は `solver_density_cuda/tools/docker_build.sh dev` で Docker イメージを再ビルドする。

## ParaView プラグイン (マッハ数 / シュリーレン / Q 値 / ヘリシティ)

`res_*.xmf` には保存量と原始変数しか入っていないため、可視化用の派生量は
`solver_density_cuda/tools/paraview/forge_filters.py` の Python プラグインで計算する。
フィルタ名は **Forge Derived Quantities** (Filters > Alphabetical)。

読み込みは Tools > Manage Plugins > Load New... で `forge_filters.py` を選ぶ
("Auto Load" を有効にすれば次回起動から自動で入る)。

**WSL native の ParaView は起動方法に注意**。Ubuntu 24.04 の paraview 5.11.2 は Python 3.12 上で動くが、
ParaView 同梱の `paraview/detail/pythonalgorithm.py` が Python 3.11 で削除された `inspect.getargspec` を
import するため、素で `paraview` と打つと Python プラグインの読み込みが必ず次のエラーで失敗する。

```
WARN| Failed to load Python plugin: Failed to import `paraview.detail.pythonalgorithm`.
ImportError: cannot import name 'getargspec' from 'inspect'
ERR | .../forge_filters.py: invalid ELF header
```

native ではラッパー経由で起動すれば、この shim (`tools/paraview/pvshim/sitecustomize.py` を PYTHONPATH に追加)
が入るので Manage Plugins も auto load もそのまま通る。

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_paraview_native.sh                       # native (WSL)
./tools/run_paraview_native.sh ../case/08.bump/run_0006_slau_loM_imp_cfl5/res_8000.xmf
```

Docker 経由 (`./tools/run_paraview_gui.sh`) の ParaView は Python 3.10 なので shim 不要でそのまま読める。
どうしても素の `paraview` から使いたいときは、Macros > Add new macro... で
`tools/paraview/macro_load_forge_filters.py` を登録し、そのボタンから読み込む
(マクロ側で同じ shim を入れてからプラグインを読む)。

出力される配列:

| 配列 | 定義 |
| --- | --- |
| `U` | 速度ベクトル (`Ux`,`Uy`,`Uz` を束ねたもの) |
| `Mach`, `sound_speed` | $M=\lvert U\rvert/a$。$a$ は配列 `sonic`、無ければ $a=\sqrt{\gamma P/\rho}$ |
| `grad_ro` | $\nabla\rho$ (ベクトル) |
| `schlieren_mag` | $\lvert\nabla\rho\rvert/\max\lvert\nabla\rho\rvert$ (0–1) |
| `schlieren_dir` | $(\hat n\cdot\nabla\rho)/\max\lvert\hat n\cdot\nabla\rho\rvert$ (−1–1)。$\hat n$ は Schlieren Direction (既定 x=流れ方向) |
| `schlieren` | $\exp(-k\,\lvert\nabla\rho\rvert/\max\lvert\nabla\rho\rvert)$ の疑似シュリーレン画像 (白黒反転して表示すると実験像に近い) |
| `vorticity`, `vorticity_mag` | $\omega=\nabla\times U$ |
| `Q`, `Q_norm` | $Q=\tfrac12(\lVert\Omega\rVert^2-\lVert S\rVert^2)$、`Q_norm` はその $\tfrac12(\lVert\Omega\rVert^2+\lVert S\rVert^2)$ 正規化 (−1–1) |
| `lambda2` | $S^2+\Omega^2$ の中間固有値 (既定 OFF) |
| `helicity`, `helicity_norm` | $U\cdot\omega$ と $U\cdot\omega/(\lvert U\rvert\lvert\omega\rvert)$ |

使い方のポイント:

- **勾配の出所**: `Use Solver Gradients` (既定 ON) ならソルバが出力した `dUxdx…`/`drodx…` をそのまま使い、
  無ければ `vtkGradientFilter` でメッシュから計算する。ソルバの離散勾配と可視化を一致させたいときは ON、
  出力に勾配が無い run や外部データでは OFF でよい。
- **cell / node 両対応**: `ro` が Cell Data にあれば Cell、Point Data にあれば Point として処理する。
  cell モードの結果を滑らかに描きたいときは後段に `Cell Data to Point Data` を挟む。
- **シュリーレンの正規化**は「その時刻・全ブロックの最大値」で行うため、時系列でスケールが変わる。
  アニメーションでレンジを固定したいときは色マップの Rescale を切る。`Schlieren Exponent` $k$ は
  大きいほど弱い密度勾配を強調する (既定 15)。
- **Q 値の等値面**は `Q` (次元 1/s²) に閾値を決めにくいので、まず `Q_norm > 0.1` 程度で当たりを付けるとよい。
- 検算: ABC (Beltrami) 流 $\omega=U$ の解析場で `vorticity`・`grad_ro`・`Q` が 2 次精度で一致し、
  `helicity_norm`=1.000000 になることを確認済み。

## ケースごとの注意

- `run_case.sh` などのケース専用ヘルパーがあるなら、それを優先する。
- 実行前に、既存の `run_*` を複製して新しい `run_*` を作る。既存ディレクトリへそのまま再実行しない。
- HDF5 変換後のファイル名は、必ず `solverConfig.yaml` の `mesh.meshFileName` と `mesh.valueFileName` に合わせる。
- `forge` は必ず `solverConfig.yaml` がある新しい `run_*` ディレクトリで実行する。
- `forge` 実行後は、同じ新しい `run_*` ディレクトリで `residual_history.csv` から `residual_history.png` を生成する。
- `run_*` が複数ある場合は、ユーザー指定のもの、または使いたいスキームに対応するものを選ぶ。

## Fluent (CFF HDF5) メッシュの取り込み

世の中で配布されている ANSYS Fluent の HDF5 メッシュ (`*.msh.h5` / `*.cas.h5`, いわゆる CFF) を forge 入力 HDF5 に変換するには `solver_density_cuda/tools/fluent_h5_to_forge.py` を使う。Gmsh を経由せず、Fluent の面 (face) 情報から forge の plane/cell/bcond を直接生成する (体積は発散定理で計算するので多面体セルでも解は回る。ただし多面体は forge の結果可視化では正しく描画されない)。

手順:

```bash
# 1) 構造とゾーン一覧を確認 (名前・zoneType・境界/内部・単位)
python3 solver_density_cuda/tools/fluent_h5_to_forge.py dump mesh.cas.h5

# 2) ゾーン -> physID 対応表と bcondConfig.yaml の雛型を出力
python3 solver_density_cuda/tools/fluent_h5_to_forge.py template mesh.cas.h5 \
        -o mapping.yaml --bcond bcondConfig.yaml

# 3) mapping.yaml / bcondConfig.yaml を編集後、変換 (単位換算は --scale)
python3 solver_density_cuda/tools/fluent_h5_to_forge.py convert mesh.cas.h5 mesh.h5 \
        --map mapping.yaml --bcond bcondConfig.yaml --scale 0.0254 \
        --ro 1.0 --ux 100 --p 101325
```

要点:

- 境界は **physID (整数)** で実行時に `bcondConfig.yaml` と突き合わされる (`boundaryCond.cpp::readBcondConfig`)。Fluent ゾーン名 -> physID の対応 (`mapping.yaml`) と、同じ physID をもつ `bcondConfig.yaml` を必ず用意する。`kind`・境界値は `bcondConfig.yaml` が正本 (H5 の `bcondKind` 属性は上書きされる)。
- 出力 `mesh.h5` 名は `solverConfig.yaml` の `mesh.meshFileName` / `mesh.valueFileName` に合わせる。
- 座標単位が m でない場合 (Fluent は inch/mm が多い) は `--scale` で SI へ換算する (`dump` が単位を警告する)。`--double` は double ビルド向けに float64 で書く。
- 初期値 (`--ro/--ux/--uy/--uz/--p`/`--t`, RANS は `--rok/--roomega`) は一様値で `/VALUE/*` に焼き込まれる。`wall_dist` は `kind` が `wall` で始まる境界面から自動計算する。
- CFF のレイアウトは Fluent バージョンで揺れる。読めない場合は `dump` で構造を確認し、必要なら `--nodes-path` 等で補う。

### 変換メッシュを安定に回す実行レシピ (重要)

非構造/非直交の Fluent メッシュを発散させずに回すには、実行設定が肝になる (検証済み: fan(hex) 完全収束、StaticMixer(tet+prism) 健全完走)。

1. **1 次風上で投入する** (`space.convMethod: 0`)。歪んだセルで MUSCL(2)は overshoot しやすい。収束後に必要なら 2 次へ上げる。
2. **初期場の流速を入口流れとおおむね一致させる**。静止場に速度入口を当てると危険 (起動衝撃)。入口が `Ux=50` なら初期場も `U=(50,0,0)` 程度にする。
3. **定常は陰解法 (`timeIntegration: 11`, `blockDPLUR: 1`, `cfl_pseudo`≈1, `nStepInner`≈10)**。**陽解法 (`timeIntegration: 3`) + 定常 (`unsteady: 0`) は局所時間刻みで不安定**になりやすい。陽解法を使うなら `unsteady: 1` (時間精度)。
4. **非直交メッシュは `space.pRef` に動作静圧を入れる** (free-stream 保存)。既定 0 で従来挙動。
5. **出口は `outlet_statPress`** (静圧指定)。`outflow` は超音速用なので亜音速には使わない。
6. **`outlet_statPress` には逆流用の `Pt`/`Tt` も必ず設定する**。逆流 (出口で内向き) が生じた面では `Pt`/`Tt` から逆流ガスの静的状態を作るため、`Ps` だけだと `Pt=0`→密度 0→NaN になる。混合・剥離など出口で逆流が出うる流れでは必須。

## forge は run_case.sh 経由で実行する (収束チェック強制)

forge を直接 `build/forge` で起動せず、**`solver_density_cuda/tools/run_case.sh <run_dir>`** 経由で実行する。ラッパーが実行後に必ず `check_convergence.py` を走らせ、VERDICT を `<run_dir>/CONVERGENCE_VERDICT.txt` に残し残差プロットも生成する。長時間 run は `nohup solver_density_cuda/tools/run_case.sh <run_dir> &`。

`.claude/settings.json` のフックでこれを強制する: **PreToolUse** が直接 `build/forge` 実行を deny、**Stop** が「最近 forge を回したのに収束チェックしていない run」があるとターン終了を block する (実行方法に依らずファイル状態で検査)。これにより「衝撃波位置の安定や `rms_ro` 単独で収束と判断する」近道を構造的に防ぐ。**収束/一致を主張する応答では VERDICT 行を引用する**。

## メッシュ変更後の restart (cross-mesh interpolation)

メッシュを変えた (quad↔tri↔構造化, 解像度変更, スロット数変更など) ときは、**uniform 初期値から計算を始めない**。超音速/衝撃波/SST は uniform IC から数ステップで発散する (実例: tri 化したメッシュを uniform から起動し step 4 で `res_nan` ダンプ)。

過去の収束済み場を新メッシュへ移植してから計算する:

```bash
# 新メッシュを convert した直後の入力 h5 に、過去 run の res を最近傍interpolateして貼る
python3 solver_density_cuda/tools/interp_field.py <past_run>/res_NNNN.h5 <new_run>/mesh.h5
```

- 保存量 (ro,roUx,roUy,roUz,roe)・乱流 (roK,roOmega)・スカラー輸送 (roY*) を移植。
- **wall_dist は移植しない** (新メッシュで convert 時に計算された値を使う; 別メッシュの距離は不整合)。
- SRC は res (primitives) でも入力 h5 (conserved) でも可。scipy `cKDTree` 最近傍。
- 同一メッシュ内での restart (背圧変更など) は同ケースの `restart_field.py` を使う。

## メッシュ品質チェック (計算前・必須)

メッシュを `convertGmshToForge` で HDF5 化したら、**計算を投入する前に必ずメッシュ品質を確認する**。歪んだ/極端に細長いセルは発散・精度劣化・非物理解の原因になる。専用ツールで判定する:

```bash
python3 solver_density_cuda/tools/check_mesh_quality.py <run_dir>/mesh.h5
```

- **アスペクト比 (AR) ≤ 1000 を目標**。最長辺/最短辺。境界層クラスタリングで薄いセルを作るときに監視する。
- **スキューネス (equiangle skew) ≤ 0.9 を目標**。四角形の内角の直交からのずれ (0=直交, 1=退化)。
- ツールは AR・skew の max / p99 / 違反セル数を出し `VERDICT: PASS / SOFT-PASS / FAIL` を返す。**FAIL なら計算を投入しない**。`SOFT-PASS` (違反<0.1%) は局所外れ値として許容しうるが、場所を確認する。
- 近壁を細分化 (wall-resolved, 第一セル数μm) すると AR が増えやすい。**接線方向セルを細かくしすぎず、AR が 1000 を超えないよう第一セル厚と接線長のバランスを取る** (高 Re では y+~1 と AR≤1000 は両立しないことがあり、その場合は y+~30-80 + `wallTreatmentSST=1` を選ぶ)。
- 「メッシュできた」「収束した」と報告する応答には、本ツールの品質 VERDICT も根拠として併記する。
- なお `check_mesh_quality.py` は **primal (cell) 変換の h5 専用**で、median-dual (node) 変換した
  h5 を渡すと `CONNE が NumberOfElements より短い` で落ちる (ツール側の制約)。node メッシュの
  品質ゲートは、**同じ `.msh` を `discretization: "cell"` で変換した h5** に対して実行する
  (primal 形状は同一なので AR/skew の判定として妥当)。

## node モードのメッシュは平面 2D で作る (2D 問題・必須)

**node-centered (median-dual) で 2D 問題を解くときは、cell で常用する「1 層押し出し + spanwise slip」の
メッシュを使ってはならない。** 押し出しメッシュは node では spanwise に**ノードを 2 枚**作り、
2 ノードしかない方向では **2 次精度 MUSCL の左右再構成が厳密に一致して上流差分の散逸が完全に消える**。
その結果 spanwise 市松モードだけが無減衰になり、丸め誤差から指数成長して計算を壊す。
**リミッタ (barth/venkata) も CFL 低減も LSQ 勾配も効かない** (機構と実測は
[`methods/discretization.md`](../methods/discretization.md) §5.1)。

- **2D 問題**: `.geo` で `Extrude` せず、`Physical Curve` で境界タグを付けた**平面 2D メッシュ**を作る
  (例: [`case/26.flat_plate_sst/mesh/flat_plate_yp30_planar.geo`](../case/26.flat_plate_sst/mesh/flat_plate_yp30_planar.geo),
  [`flat_plate_planar.geo`](../case/26.flat_plate_sst/mesh/flat_plate_planar.geo))。
  `bcondConfig.yaml` から spanwise の `side1`/`side2` エントリを落とすこと。
- **3D で均質方向を持つ場合**: その方向に**最低 3 ノード (2 層以上)** 確保する。
- **症状の見分け方**: 発散したとき、まず**同一 (x,y) にある spanwise ノード間の値のばらつき**を測る。
  疑似 2D なら厳密に 0 のはずで、これが成長していれば本件。局所的な「前縁の暴走」等に見えても
  最終スナップショットだけ見ると発生源を取り違える (case/26 の実例)。
- 既存の node run を複製するときは、**種メッシュが押し出し版でないか**を確認する
  (`CELLS/centCoords` の z ユニーク値が 2 個なら該当)。

## 結果ディレクトリの明示 (投入時・まとめ時)

検証計算の所在が後から分からなくならないよう、計算の **投入時** と **結果まとめ時** の両方で、出力先 `run_*` をリポジトリルートからの相対パスで明示する (正本ルールは `AGENTS.md` の「計算結果ディレクトリの明示」)。

- **投入時**: 起動する応答内で各 run の出力先パスを列挙し、run ごとの設定差分 (スキーム・`convMethod`・物性など) を 1 行で添える。
  - 例: `case/05.sod_shock_tube/run_0004_slau_tracer_n2n2/` — SLAU, `species:[N2,N2]` 識別子トレーサ。
- **まとめ時**: 結論・比較表とあわせて、根拠となった `run_*` パスと主要成果物 (`res_*.h5`, `residual_history.png`, 後処理図) を併記する。複数 run 比較では「どの数値/図がどの run か」を取り違えないよう run パスを明記する。
- 一時 run・破棄予定 run もパスと用途を明記し、所在不明にしない。

## NaN / 発散チェック (投入直後・まとめ前)

計算が exit 0 で完走しても妥当とは限らず、**初期ステップから発散している**ことが多い。結果を使う前に必ず NaN/Inf を確認する (正本ルールは `AGENTS.md` の「NaN / 発散チェック」)。

- 投入直後に序盤 (数十ステップ) で `residual_history.csv` の `rms_*` が NaN になっていないか早期確認する。`step 0`/`step 1` から NaN なら初手発散。
- まとめ前に最終・中間 `res_*.h5` の `VALUE/*` (`ro`,`P`,`T`,`roe`,勾配など) の NaN/Inf と物理妥当性 (静圧≤全圧, `ro>0`, `T>0`) を確認する。
- `max cfl` ログが有限でも残差は NaN のことがある。CFL だけで安定と判断しない。
- 発散時は「収束」と報告せず、最初に NaN が出たステップ・境界/セルを切り分けて原因 (BC のゼロ割・IC と BC の不整合・CFL 過大・メッシュ不良) を特定してから対処する。

## 運用上の位置づけ

Copilot が計算手順を説明するときは、この文書を標準の参照先とします。個別ケースの検証方針や、Docker と native 実行の使い分けは、別の関連文書を優先して参照します。