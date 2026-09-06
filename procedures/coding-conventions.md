# Forge Coding Conventions

`solver_density_cuda` (密度ベース CUDA CFD ソルバ) の開発で従うソース構成・命名規約・
ビルド/テスト手順をまとめる。本文の運用ルール全体の正本は [AGENTS.md](../AGENTS.md)。
理論・実装解説は [`methods/`](../methods/index.md) を参照する。

## ソース構成

```
solver_density_cuda/
├── main.cpp              # エントリポイント (設定読込→メッシュ→初期化→時間ループ→出力)
├── *.cpp / *.hpp         # CPU 側 (variables, output など)
├── cuda_forge/           # GPU カーネル群 (静的ライブラリ libcuda_forge)
│   ├── <feature>_d.cu    # デバイス実装
│   └── <feature>_d.cuh   # 対応ヘッダ
├── probe/                # 点プローブ (静的ライブラリ libprobe)
├── input/                # 設定パース (solverConfig)、壁距離計算
├── mesh/                 # メッシュ I/O・変換 (convertGmshToForge)・分割 (mesh_part)
├── output/               # 結果出力 (HDF5 / XMF)
├── common/               # 汎用ユーティリティヘッダ
├── third_party/HighFive/ # vendored (git submodule)
├── tools/                # ビルド/プロファイル/GUI 補助スクリプト
└── tests/                # 回帰テスト (tests/README.md 参照)
```

ビルド成果物: `forge` (ソルバ)、`convertGmshToForge` (Gmsh→HDF5)、`mesh_part` (分割)。

## 命名規約 (既存コードに準拠)

- **ファイル**: GPU 実装は `<feature>_d.cu` と `<feature>_d.cuh` の対で置く。CPU は `.cpp`/`.hpp`。
- **デバイス/ホスト**: デバイス側のスコープには接尾辞 `_d` を付ける (例: カーネル `SLAU_d`、device ポインタ `c_d` / `p_d`)。ホストから呼ぶ起動関数は `<feature>_d_wrapper` とする。
- **カーネル/デバイス関数**: `__global__ void <Name>_d(...)`、補助は `__device__` / `__device__ __forceinline__`。
- **浮動小数型**: 実数は `flow_float` (typedef) を使い、`float`/`double` を直書きしない。
- **流れ場変数**: セル/プレーン量は文字列キーの map (`variables::c`, `c_d` 等) で保持し、キーは `"ro"`, `"roUx"`, `"roe"`, `"roK"`, `"roOmega"` のような保存量名を用いる。時間レベルは接尾辞 `N` / `NN` で表す。
- **スキーム識別子**: `SLAU`, `ROE`, `AUSM+`, `AUSM+UP`, `KEEP_SLAU` を設定値・run 名で統一して使う (run 名規則は [calculation-workflow](calculation-workflow.md))。

## コメント・言語

- コメント本文は日本語。識別子・関数名・スキーム名は原語のまま。
- 理論的背景はコメントに長文で書かず、対応する [`methods/<area>/theory.md`](../methods/index.md) を参照する。

## ビルド

通常開発は Docker を基本とし、GPU プロファイリング時のみ native を優先する (方針は
[forge-development-environment.md](development-environment.md))。

### Docker

```bash
cd solver_density_cuda
# イメージ (初回 / Dockerfile 変更時)
./tools/docker_build.sh dev      # forge-solver:cuda-dev (ビルドだけなら base で足りる)
# コンテナ内でビルド (build/ に出力)
./tools/build.sh
```
`compose.yml` の `dev-cuda` サービス、または `.devcontainer/` でも同イメージを使える。

### WSL / Linux native

```bash
cd solver_density_cuda
./tools/build_native_wsl.sh        # .build-native/relwithdebinfo/ に出力
```
主な環境変数:

| 変数 | 用途 |
| --- | --- |
| `FORGE_CUDA_ARCHITECTURES` | 対象 CUDA arch (既定 `75;86;89;90`) |
| `FORGE_BUILD_JOBS` | 並列ジョブ数 (既定 `nproc`) |
| `FORGE_NATIVE_BUILD_TYPE` | ビルドタイプ (既定 `RelWithDebInfo`) |
| `FORGE_NATIVE_BUILD_DIR` | 出力ディレクトリ上書き |

依存: HDF5 / yaml-cpp / METIS (必須)、HighFive (vendored)。

## テスト / 回帰

GPU 必須。詳細は [`solver_density_cuda/tests/README.md`](../solver_density_cuda/tests/README.md)。

```bash
cd solver_density_cuda
python3 tests/regression/run_regression.py --smoke   # 短ステップ・有限性のみ
python3 tests/regression/run_regression.py           # baseline と比較 (PASS/FAIL)
python3 tests/regression/run_regression.py --update-baseline  # 数値を変える変更の後
```

数値スキームや時間積分を意図的に変えた場合は、`--update-baseline` で golden を更新し、
差分の根拠を commit / plan の変更ログに記載する (AGENTS.md の開発フロー参照)。
回帰は AGENTS.md の複製ルールに従い、ハーネスが一時 run を作って実行する。

## commit / PR

- 1 commit は 1 つの論理変更にまとめ、メッセージ本文は日本語可・命令形の要約から始める。
- 数値・設計を変える変更は、対応する [`methods/`](../methods/index.md) と [`.github/plans/`](../plans/README.md) を先に更新してから実装する (AGENTS.md「開発フロー」)。
- PR は対象・非対象、検証に使ったケースと結果 (回帰 PASS/FAIL、必要なら residual PNG) を記載する。

## 今後の課題 (backlog)

- GPU CI (self-hosted runner) もしくは build-only の CPU CI。
- 回帰指標の拡充 (壁面分布・line profile 比較)。
- 既存検証ケースの設定が旧キー (`nStep` / `nInnerLoop`) のまま残っているものがあり、現行ソルバは `nStepOuter` / `nStepInner` を要求する。順次更新する (回帰ハーネスは複製した run 内で自動正規化して回避している)。
- [`methods/`](../methods/index.md) の薄い領域 (poisson / diffusion) と未整備領域 (GPU メモリ管理 / I/O / 入力スキーマ) の拡充。
