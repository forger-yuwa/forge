# 諮問ブリーフ: 変換器が常に double で幾何を書く設計 (2026-10-10)

- 依頼者: 主セッション (AGENTS.md のエスカレーション条件 1: plan §4・§6 を新しく書く)
- plan: `plans/active/architecture-float-state-double-geometry.md` の **§4.6 (新しい設計)・§6.14 (合格条件)**、§5.1 #13
- ユーザの決定 (2026-10-10): 「変換器が常に double で幾何を書くようにする」「手順で FP64 の変換器を使うことを必須にする」の両方をやる。
- コード: ブランチ `feature/nozzle-wall-fit-and-pipeline` の HEAD (段 ①〜④ と commit の診断まで commit 済み)

## 0. 問い

1. §4.6 の方針 (変換器はビルドの型に依らず幾何を double で計算・出力、ソルバは float32 の格子に警告、手順で確認) に抜けや危険はあるか。
2. 実装のやり方 (a) 変換器だけを double の型で別にコンパイル / (b) 変換器の幾何の経路の型を double にする — のどちらを勧めるか。下の §1 の制約 (変換器が `cuda_forge` を link している) を踏まえて。
3. §6.14 の合格条件は事前登録として十分か。特に「float の変換器の出力 = FP64 の変換器の出力 (ビット単位)」は成り立つと見込めるか (成り立たないなら、どの量の、どの程度の差を許すべきか)。

## 1. 観測事実

- 変換器 `convertGmshToForge` のソース (`solver_density_cuda/CMakeLists.txt:114-126`): `mesh/convertGmshToForge.cpp`・`mesh/mesh.cpp`・`variables.cpp`・`boundaryCond.cpp`・`input/solverConfig.cpp`・`input/calcWallDistance_kdtree.cpp`。`cuda_forge` (float のビルドでは float の型でコンパイルされた静的ライブラリ) を link している。
- `flowFormat.hpp` は `flow_float`・`geom_float` を typedef する。FP64 のビルドは typedef を double に替えた全体のビルド。
- 変換器の幾何の計算は `mesh/gmshReader.hpp` (2627 行、`geom_float` 119 か所) と `input/calcWallDistance_kdtree.cpp` (`geom_float` 30 か所)。双対の面ベクトル・面積・面重心・体積・重心・境界の半割面は `std::vector<geom_float>` の配列 (`gmshReader.hpp:45-56`)。座標は `geom_float` で読む (`:505-514`、double のときだけ `stod`)。出力は `writeInputH5` (`:2322-`) が `/MESH/COORD`・`/PLANES/*`・`/CELLS/*`・`/DUAL/*` を書く。
- ソルバは段 ① から、同じデータセットを double の写し (`coord64`・`surfVect64`・`planeCent64`・`cc64` など) にも読む (`mesh/mesh.cpp:236-260`)。float32 のデータセットなら写しは float の値そのもの。
- case/45 の格子は FP64 のビルドの変換器 (`~/forge-wallfit-bin-fp64`) で作っていて double。ほかの標準ケースの格子は float32 (V6 の下調べ)。
- 変換器の HDF5 は同じバイナリでもバイト単位では再現しない (メモ converter-h5-not-byte-reproducible)。データセットの値で比べる。
- 変換器は終了時に `GPUassert: invalid argument cudaWrapper.cu 81` を出すことがある (AWS で既知、h5 は書けている)。変換器が CUDA の初期化をしている。

## 2. 期待値と出典

- §4.2a の取り決め: 差を取る幾何は double の値から 1 回だけ丸める。格子の HDF5 が float32 だと、その「double の値」が存在しない。
- FP64 のビルドの変換器の出力は今と同じ値であるべき (既存の FP64 の格子と run の再現性)。

## 3. 呼び出し側の案 (検証していない)

- (b) を第一候補と考えている: `gmshReader.hpp` と `calcWallDistance_kdtree.cpp` の幾何の型を専用の `conv_float = double` にし、`writeInputH5` で double のまま書く。`mesh` の構造体 (ソルバと共有) に入れる経路があれば、そこだけ変換器側で double の配列を持つ。
- (a) は `cuda_forge` を double でもう一度コンパイルするか、変換器から `cuda_forge` への依存を外す必要があり、ビルドの時間が倍になるか構造の変更が大きいと見ている。

## 4. 読んでよいファイル

- `plans/active/architecture-float-state-double-geometry.md` (§4.2a・§4.6・§5.1 #13・§6.14)
- `solver_density_cuda/CMakeLists.txt`、`solver_density_cuda/flowFormat.hpp`
- `solver_density_cuda/mesh/gmshReader.hpp` (grep で該当箇所だけ)、`solver_density_cuda/mesh/convertGmshToForge.cpp`、`solver_density_cuda/input/calcWallDistance_kdtree.cpp`
- `solver_density_cuda/mesh/mesh.cpp` (`:230-380`)
- 巨大なファイル (`*.h5`・ログ) は読まないこと。
