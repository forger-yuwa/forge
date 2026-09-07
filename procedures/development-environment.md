# Forge Development Environment

この文書は、Forge の開発時に Docker と native 環境をどう使い分けるかを整理したものです。

## 基本方針

- 通常のコード修正、ビルド、ケース実行は Docker コンテナを基本とする。
- `solver_density_cuda` の開発環境は、`Dockerfile.cuda` の `dev` ステージとそのコンテナ実行を既定の流れとして扱う。
  イメージは 1 つの多段 Dockerfile から `base` (ビルド最小) / `cloud` (+Python・gmsh) / `dev` (+ParaView・GUI) を切り出す。
  ビルドは `solver_density_cuda/tools/docker_build.sh {base|cloud|dev|all}`。詳細は [`solver_density_cuda/README_docker.md`](../solver_density_cuda/README_docker.md)。
- コンテナ内でビルドした成果物を使ってケースを実行する運用を標準とする。

## Docker を既定とする理由

- リポジトリのビルド手順や依存関係が Docker 前提で揃っている。
- 開発者ごとの差分を減らしやすい。
- `run_case.sh` のようなケース専用スクリプトも、Docker ベースの実行を前提にしているものがある。

## 例外: NVIDIA の計測ツールは native を優先する

NVIDIA 提供の GPU 速度計測・プロファイリングツールを使う場合は、Docker 内でうまくいかないことがあるため、WSL native または Linux native の環境を優先する。

対象として想定するツール:

- Nsight Compute (`ncu`)
- Nsight Systems (`nsys`)
- その他、NVIDIA 提供の GPU 計測ツール

通常の solver 実行は Docker を基準にしてよいが、これらのツールで計測したいときは、まず native 実行へ切り替える。

## 速度評価・ベンチマークは native を既定とする (ルール)

forge の **計算速度を評価・比較・プロファイルする作業は、原則 WSL native / Linux native で行う**。Docker での計測は、計測ツールがイメージに同梱されていなかったり、計測値がコンテナ層の影響で不安定になるため、速度評価の基準環境としない。

対象となる作業:

- 1 ステップあたりの所要時間や `nStep` あたりの total 時間の計測
- スキーム・設定変更前後の速度比較 (control run を含む)
- カーネル単位のボトルネック特定 (`ncu` など)
- 内蔵プロファイラ (`FORGE_PROFILE`) を使ったセクション別計時

native で計測するときの標準手順:

1. `solver_density_cuda/tools/build_native_wsl.sh` で native バイナリ (`.build-native/relwithdebinfo/forge`) を最新ソースから再ビルドする (Docker の `build/forge` とは別物。stale に注意)。
2. 計測対象の新しい `run_*` ディレクトリで native バイナリを実行する。
3. セクション別の計時は `FORGE_PROFILE=1` (詳細は `FORGE_PROFILE_VERBOSE=1`) を付けて起動し、正常終了時に出る "Runtime Profile Summary" を読む (summary は正常終了時のみ出力されるため、計測 run は完走する step 数にする)。
4. カーネル単位の内訳が要るときは `ncu`、タイムラインが要るときは `nsys` を native で使う。

この環境 (WSL2 + Ubuntu, RTX 3060 / CC 8.6) での計測ツールの実情:

- `ncu` (Nsight Compute) は単体で動作する。ウォームアップ後の代表区間に `--launch-skip` / `--launch-count` で絞り、`--metrics gpu__time_duration.sum,sm__throughput.avg.pct_of_peak_sustained_elapsed,gpu__compute_memory_throughput.avg.pct_of_peak_sustained_elapsed` でカーネル別の実時間と compute/memory バウンド傾向を取るのが軽い。
- `nsys` は target パッケージのみ入っており host 側 importer (`QdstrmImporter`) が無いため、native では `.qdstrm` を収集できても `.nsys-rep` / `.sqlite` への変換・`nsys stats` ができない。タイムラインが必要なときは Windows 側 Nsight Systems GUI で `.qdstrm` を開く (`tools/check_wsl_profile_target.sh` 参照)。当面のカーネル計測は `ncu` と `FORGE_PROFILE` を優先する。
- `nvprof` は CC 8.0 以上 (RTX 3060 = 8.6) を kernel profiling 対象にできないため使わない。

## native 実行時の目安

- WSL native / Linux native では `solver_density_cuda/tools/build_native_wsl.sh` をビルドの起点として扱う。
- Ubuntu 24.04 + CUDA toolkit 12.0 の実績では、CUDA host compiler に `gcc-12` / `g++-12` を使う構成が安定している。
- GPU アーキテクチャは環境に応じて `FORGE_CUDA_ARCHITECTURES` で上書きする。

代表例:

```bash
cd /home/sano/work/forge/solver_density_cuda
FORGE_BUILD_JOBS=1 FORGE_CUDA_ARCHITECTURES=86 ./tools/build_native_wsl.sh
```

この例は native 側のビルド確認用であり、通常開発の既定経路を Docker から native に置き換える意図ではない。

## Docker 実行時のファイル所有者ルール

`docker run` を含むすべてのスクリプト・コマンドには、必ず `--user "$(id -u):$(id -g)"` を付加すること。

```bash
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$repo_root":/workspace \
  ...
```

**理由:** Docker コンテナはデフォルトで root として実行される。`--user` を省略すると、コンテナ内で生成したファイル (`*.h5`, `*.msh`, `*.csv`, `*.png` など) がホスト側でも root 所有になり、通常ユーザーから読み書きできなくなる。

**適用範囲:**
- `case/` 配下の `run_case.sh`
- `solver_density_cuda/tools/` 配下のラッパースクリプト (既適用)
- エージェントや手動で書く `docker run` 呼び出しすべて

コンテナ内でツールが `$HOME` を参照する場合は、`-e HOME=/tmp/forge-home` も合わせて指定する (ParaView など)。

## 使い分けの実務ルール

- ビルドと通常のケース確認: まず Docker で行う。
- Gmsh / ParaView の GUI 起動: `solver_density_cuda/tools/` 配下のラッパースクリプトをホスト側から使う。
- NVIDIA の計測ツールを試す: Docker 内に固執せず、WSL native または Linux native に切り替える。
- Docker で計測が失敗した場合: 例外扱いではなく、native へ切り替えるのが既定対応。

## クラウド (AWS EC2) での実行

AWS の GPU インスタンス (基準: `g5.xlarge`, us-east-1) で Docker ビルド・計算投入・速度計測を行う手順は
[`cloud-aws-gpu.md`](cloud-aws-gpu.md) を参照。クラウドの Linux native ホストでは本文書の
「速度評価は native」ルールをそのまま適用でき、WSL で不能だった `nsys` の解析までインスタンス内で完結する。
クラウド用イメージは `Dockerfile.cuda` の `cloud` ステージ (ParaView GUI 無し、gmsh はヘッドレス同梱)。`./tools/docker_build.sh cloud` で作る。
メッシングは gmsh スクリプトで、可視化は pvserver + SSH トンネルでクラウド完結できる (手順は同文書)。

## 関連文書

- 計算準備と `forge` 実行手順: `procedures/calculation-workflow.md`
- 標準検証ケース: `procedures/verification/README.md`
- Docker 開発環境の補足: `solver_density_cuda/README_docker.md`
- native ビルドスクリプト: `solver_density_cuda/tools/build_native_wsl.sh`
