# Docker イメージの多段分割とビルド時間短縮

## メタ

- **area**: `architecture` (開発環境・ツーリング)
- **status**: `done`
- **related_docs**:
  - [`solver_density_cuda/README_docker.md`](../../solver_density_cuda/README_docker.md)
  - [`procedures/development-environment.md`](../../procedures/development-environment.md)
- **related_plans**: [`tooling-cloud-gpu-env.md`](../active/tooling-cloud-gpu-env.md)
- **created**: `2026-09-06`
- **owner**: ユーザー + Claude

## 1. 目的

Docker イメージのビルドが遅く、gmsh の取得に失敗すると apt からやり直しになる状態を解消する。
あわせてコンテナ内の `forge` ビルド時間を詰める。完了時には、イメージを段階的に作れて
失敗時の再実行範囲が局所化され、ソルバのフルビルドが従来の 1/2 (キャッシュ温なら 1/17) になる。

## 2. スコープ

- **やる**: `Dockerfile.cuda.dev` / `Dockerfile.cuda.cloud` を多段 `Dockerfile.cuda` に統合、
  レイヤ分割、`.dockerignore` の追加、`tools/build.sh` への CUDA arch 絞り込みと ccache 導入、
  未宣言だった `libboost-dev` の明示、`boost::scoped_array` → `std::unique_ptr` 置換。
- **やる (2026-09-06 追加)**: `mesh/gmshReader.hpp` の `boost::split` 置換と `libboost-dev` の削除。
- **やらない**: CUDA base イメージ (11.4 GB) 自体の縮小。CPU 専用イメージ。`nsys` の同梱。
  ビルド対象外の `solvePoisson_amgcl.*` (boost::property_tree + amgcl) の復活。

## 3. 関連 docs と前提

- イメージ構成と使い方の正本は [`solver_density_cuda/README_docker.md`](../../solver_density_cuda/README_docker.md)。
- クラウド (AWS) 側の運用は [`procedures/cloud-aws-gpu.md`](../../procedures/cloud-aws-gpu.md) と
  [`tooling-cloud-gpu-env.md`](../active/tooling-cloud-gpu-env.md) が引き続き所管する。
  本計画は cloud イメージを「独立ファイル」から「`Dockerfile.cuda` の `cloud` ステージ」へ移した。

## 4. 設計方針

### 4.1 3 段構成

| stage | tag | 内容 |
| --- | --- | --- |
| `base` | `forge-solver:cuda-base` | コンパイラ + HDF5 / yaml-cpp / METIS / boost(headers) + ccache |
| `cloud` | `forge-solver:cuda-cloud` | base + Python (numpy/h5py/matplotlib) + gmsh (ヘッドレス) |
| `dev` | `forge-solver:cuda-dev` | cloud + ParaView / GUI 依存 / gdb |

`forge` のビルドと計算実行に必要なのは `base` だけ。GUI は `dev` にしか無い。

### 4.2 レイヤ分割の原則

「変わりにくいもの → 変わりやすいもの」の順に 1 工程 = 1 `RUN` で並べる。
最も壊れやすい gmsh の取得を `cloud` の最後に置く。従来は apt・gmsh・ParaView が
1 つの `RUN` にまとまっていたため、gmsh の失敗で apt (~2 GB) から丸ごとやり直しになっていた。

BuildKit のキャッシュマウント (`--mount=type=cache`) は採用しない。この環境には buildx が
入っておらず legacy builder では通らないため、レイヤ分割だけで目的 (失敗時の再実行範囲の局所化)
を達成する。

### 4.3 boost の扱い

`mesh/gmshReader.hpp` が `boost::split` / `boost::is_space` を使うため boost ヘッダは必須。
にもかかわらず両 Dockerfile に `libboost-dev` の記述が無く、**ParaView (dev) と
python3-matplotlib (cloud) が引きずり込む依存にたまたま乗っていただけ**だった。
GUI を外した `base` ステージで即座にビルドが落ちて発覚した。`base` に明示宣言する。

`boost::scoped_array` の 2 箇所 (`probe/point_probes`, `input/calcWallDistance_kdtree`) は
`std::unique_ptr<T[]>` へ置換した (どちらも `HAVE_KDTREE` ガード内で、現環境では未コンパイル)。

**`boost::split` も置換して boost 依存を完全に外した (2026-09-06)**。`common/stringUtil.hpp` に
`splitOnSpace(out, in)` を追加し、`gmshReader.hpp` の 14 箇所を差し替えた。
`boost::split(out, in, boost::is_space())` の規約を厳密に保つ必要がある:

- 連続する区切りを圧縮しない (`token_compress_off`) → **空トークンを残す**
- 入力が空でも要素数 1 (空文字列) を返す
- 先頭・末尾の区切りも空トークンを生む (CRLF の `\r` は区切り扱い)

`.msh` は列位置でインデックスアクセスするため、ここを崩すとメッシュ読み込みが静かに壊れる。
等価性は boost 実装との差分テストで確認した (§6)。

`solvePoisson_amgcl.cpp` / `amgcl_cuda/solvePoisson_amgcl_cuda.cu` は `boost::property_tree` を
使ったままだが、どの CMake ターゲットにも入っていない (amgcl 自体も vendoring されていない)。
これらを復活させるときは `libboost-dev` と amgcl を `base` ステージへ戻すこと。

### 4.4 ソルバのビルド時間

- **CUDA アーキテクチャ**: `CMakeLists.txt` の既定 `75;86;89;90` は全 `.cu` を 4 回コード生成する。
  `tools/build.sh` が `nvidia-smi` から compute capability を読んで 1 つに絞る
  (`FORGE_CUDA_ARCHITECTURES` で上書き可)。他機配布用 fatbin が要るときだけ複数指定する。
- **ccache**: コンパイラランチャに使う。`CCACHE_DIR=/workspace/.ccache` (= ホスト側
  `solver_density_cuda/.ccache`) なのでコンテナを捨てても残る。`FORGE_NO_CCACHE=1` で無効。

### 4.5 ビルドコンテキスト

`Dockerfile.cuda` は `COPY` / `ADD` を持たない (ソースはバインドマウント) ので、
`.dockerignore` に `*` を書いてコンテキストを空にする。従来は `build*/` や `run_*/` を含む
約 1 GB がビルドのたびに docker デーモンへ転送されていた (実測 1.058 GB → 6.6 kB)。

## 5. 実装ステップ

1. `Dockerfile.cuda` を新設 (base / cloud / dev の 3 ステージ)。`Dockerfile.cuda.dev` と
   `Dockerfile.cuda.cloud` を削除。
2. `.dockerignore` を追加。
3. `tools/docker_build.sh` を新設 (`base|cloud|dev|all`)。
4. `tools/build.sh` に CUDA arch 自動判定と ccache を追加。
5. boost 依存の除去 (`probe/point_probes.{cu,cuh,cpp,hpp}`, `input/calcWallDistance_kdtree.{cpp,hpp}`,
   `common/stringUtil.hpp` の `splitOnSpace` 追加、`mesh/gmshReader.hpp` の 14 箇所差し替え)。
6. 参照箇所の更新 (`compose.yml`, `tools/run_{gmsh,paraview}_gui.sh`,
   `tools/run_{nsys,ncu}_profile.sh`, `tools/cloud/setup_instance.sh`,
   `case/20.naca_ml/001.test/smoke_node_euler/run_case.sh`, 各 `procedures/*.md`, `README_docker.md`)。
7. `.gitignore` に `solver_density_cuda/.ccache/` を追加。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | buildx の導入 | **要 sudo のためエージェントからは実行不可**。`sudo apt-get install -y docker-buildx` (noble universe, 0.30.1)。入れれば `--mount=type=cache` で apt/pip の再ダウンロードも避けられる。legacy builder は deprecated なのでいずれ必要。導入後に `Dockerfile.cuda` の apt/pip 行へキャッシュマウントを足す |
| 2 | AWS 側での確認 | `tools/cloud/setup_instance.sh` の変更 (`docker_build.sh cloud` 呼び出し) を実機で 1 回通す |
| 3 | native ビルドへの ccache | `tools/build_native_wsl.sh` は未対応。native には ccache 未インストール |
| 4 | msh 2.2 リーダの実体化 (別件) | `readMeshFormat` は "2.2" を受理するが、以降は 4.1 のレイアウト (`$Entities` 前提) で読むため実ファイルでは `stoi` 例外で落ちる。本件の調査中に判明した**既存の不具合**で、boost 除去とは無関係 (boost 版でも同一に落ちる) |

## 6. 検証

- **ビルド**: `./tools/docker_build.sh all` が 3 ステージとも成功 (実測 7 分 54 秒, 全段コールド)。
- **ソルバビルド**: `forge-solver:cuda-base` (GUI も matplotlib も無い最小構成) で
  `./tools/build.sh` が 47/47 完走し `forge` / `convertGmshToForge` がリンクされること。
- **`splitOnSpace` の等価性**: boost 実装と新実装を同一入力に流して結果ベクタを完全一致で比較する
  (エッジケース + ランダム文字列 + 実 `.msh` の全行)。
- **メッシュ変換の回帰**: boost 版バイナリと新バイナリで同一 `.msh` を `convertGmshToForge` にかけ、
  出力 HDF5 を**データセット単位で**比較する (HDF5 のファイルレイアウトは同一バイナリでも
  非決定的なのでバイト比較は使えない)。
- **回帰**: `case/05.sod_shock_tube/run_0012_aws_repro_cell` の入力を複製して 200 step 実行し、
  `res_200.h5` を同 run の既存結果と比較する。
- **判定基準**: 主要変数の相対差が cell モードの atomicAdd ノイズ床 (~2e-6) 以内、`residual_history.csv` に NaN 無し。

## 7. 影響範囲

- 追加: `solver_density_cuda/Dockerfile.cuda`, `.dockerignore`, `tools/docker_build.sh`
- 削除: `solver_density_cuda/Dockerfile.cuda.dev`, `Dockerfile.cuda.cloud`
- 変更: `tools/build.sh`, `compose.yml`, `tools/run_{gmsh,paraview}_gui.sh`,
  `tools/run_{nsys,ncu}_profile.sh`, `tools/cloud/setup_instance.sh`,
  `probe/point_probes.{cu,cuh,cpp,hpp}`, `input/calcWallDistance_kdtree.{cpp,hpp}`,
  `common/stringUtil.hpp`, `mesh/gmshReader.hpp`,
  `.gitignore`, `case/20.naca_ml/001.test/smoke_node_euler/run_case.sh`
- docs: `README_docker.md`, `procedures/development-environment.md`,
  `procedures/calculation-workflow.md`, `procedures/coding-conventions.md`, `procedures/cloud-aws-gpu.md`

## 8. 完了条件

- [x] 関連 docs (`README_docker.md`, `procedures/*`) を更新済み
- [x] 実装・検証完了 (§6)
- [x] `status` を `done` に変更し §9 に変更ログを記載
- [x] `plans/accepted/` に配置
- [x] [`plans/README.md`](../README.md) の一覧を同期

## 9. 変更ログ

- `2026-09-06` — 初稿 + 実装 + 検証。実測値は次の通り (RTX 3060 / 12 並列)。
  - ビルドコンテキスト: 1.058 GB → 6.6 kB
  - イメージ: `base` 11.9 GB / `cloud` 12.9 GB / `dev` 14.0 GB (従来の単一 dev は 14.1 GB)。
    サイズはほぼ CUDA devel ベース (11.4 GB) が占めるため縮小効果は小さい。効くのは**再ビルド範囲**の方。
  - `forge` フルビルド: 4 arch・ccache 無し 2 分 25 秒 → 1 arch 1 分 07 秒 → `build/` 削除 + ccache 温 3.9 秒
  - sod 回帰: `res_200.h5` の `ro`/`P`/`T`/`Ux`/`roe` が参照 run と相対 ~2e-6 以内で一致
    (`Uy` は 1D 問題の数値ノイズ ±3e-4 で無意味)、`residual_history.csv` に NaN 無し。
  - **判明した罠**: `libboost-dev` はどちらの Dockerfile にも書かれておらず、
    ParaView / python3-matplotlib の依存に相乗りしていただけだった (§4.3)。
- `2026-09-06` — **boost 依存を完全に除去** (§4.3)。`common/stringUtil.hpp` に `splitOnSpace` を追加し
  `mesh/gmshReader.hpp` の `boost::split` 14 箇所を差し替え、`Dockerfile.cuda` の `base` から
  `libboost-dev` を落とした。あわせて未使用の `boost/foreach.hpp` include も削除。
  - **等価性**: boost 実装との差分テスト **532,628 ケース全一致** (エッジケース 21 +
    ランダム文字列 200,000 + 実 `.msh` 3 ファイル 332,607 行)。
  - **メッシュ変換の回帰**: boost 版バイナリと boost なしイメージで作った新バイナリで
    4 メッシュを変換し、出力 HDF5 が**データセット単位で完全一致** (435 データセット)。
    内訳: 2D 軸対称 cell (case/40 run_0032, 92) / 2D 軸対称 node (case/40 run_0034, 101) /
    3D 周期 node (case/35 run_0002, 133) / 3D RANS cell 15 MB (case/29 run_3d_conical_rans, 109)。
    **HDF5 のバイト列は同一バイナリを 2 回走らせても一致しない** (ファイルレイアウトが非決定的)
    ので、バイト比較ではなくデータセット比較で判定した。
  - **sod 回帰**: boost なしイメージでビルドした `forge` で 200 step 実行し、参照 run と
    `ro` 2.8e-6 / `P` 9.5e-7 / `T` 2.1e-6 / `Ux` 3.8e-6 / `roe` 1.4e-6、NaN 無し。
  - `base` イメージから boost が消えたことを確認 (`/usr/include/boost` 無し、boost パッケージ 0)。
  - 副産物: **msh 2.2 リーダは実質未実装**と判明 (§5.1 #4)。boost 版でも同じ位置で落ちるため本件とは無関係。
