# solver_density_cuda のコード構成

## 1. この文書の目的

この文書は、`solver_density_cuda` ディレクトリにある解析コードの現状構成を、新しくコードを読む開発者向けに整理したものである。主な目的は次の 3 点である。

- どのファイル群が何を担当しているかを把握できるようにする
- `main.cpp` を起点に、初期化から 1 ステップの更新までの流れを追えるようにする
- 今後、流束法・境界条件・出力項目などを追加するときに、どこを見ればよいかを分かるようにする

数値スキームの理論導出そのものではなく、あくまでコード構成と責務分担の説明に重点を置く。

## 2. 対象範囲

対象は `solver_density_cuda` 配下のソルバ本体である。

- 対象に含むもの: ソルバ本体、入力読込、メッシュ読込、変数管理、境界条件、CUDA カーネル群、出力、プローブ、補助ユーティリティ
- 対象に含まないもの: `case/` 配下の個別ケース説明、`post_tool/` の後処理ツール詳細、Docker 運用手順の詳細

`case/` は実行時に `solverConfig.yaml` や `bcondConfig.yaml`、変換済みメッシュ、出力ファイルが置かれる場所としてだけ最小限に参照する。

## 3. 全体像

このコードは、非構造格子上の圧縮性流れを GPU で計算するためのソルバとして構成されている。大まかな流れは次の通りである。

1. `solverConfig.yaml` を読み込んで計算条件を構築する
2. HDF5 形式のメッシュを読み込む
3. セル値・面値・GPU 側バッファを確保する
4. `bcondConfig.yaml` を読み込んで境界条件オブジェクトを初期化する
5. 構造量や従属変数を初期化する
6. CUDA 側で勾配、リミタ、流束、時間積分を繰り返し計算する
7. 必要なタイミングで HDF5/XDMF とプローブ出力を書き出す

制御の中心は `main.cpp` にあり、CPU 側は処理順の制御と I/O を担当し、重い数値計算は `cuda_forge/` の wrapper 経由で GPU 側に渡している。

## 4. ディレクトリごとの役割

| パス | 主な役割 |
| --- | --- |
| `main.cpp` | 全体の処理順を決めるエントリポイント |
| `input/` | `solverConfig.yaml` の読込、初期値設定、壁距離計算など |
| `mesh/` | メッシュ構造の保持、HDF5 メッシュ読込、前処理ユーティリティ |
| `variables.*` | セル値・面値・GPU 側配列の確保、転送、初期値読込 |
| `boundaryCond.*` | 境界条件設定の読込と CPU 側の境界条件処理 |
| `dependentVariables.*` | 保存変数から圧力、温度、音速などの従属変数を更新 |
| `setStructualVariables.*` | 体積や重心などの構造量を変数側へセット |
| `gradient.*` | 勾配計算の CPU 側入口 |
| `convectiveFlux.*` | 対流流束計算の CPU 側入口 |
| `implicitCorrection.*` | 陽解法ループ中の補正処理 |
| `update.*` | CPU/GPU 間の更新補助 |
| `output/` | HDF5/XDMF 出力、境界出力 |
| `probe/` | 点プローブ出力 |
| `cuda_forge/` | GPU 側の実計算カーネルと wrapper |
| `common/` | 文字列やベクタ処理の小さな共通ユーティリティ |
| `third_party/` | HighFive などの外部依存 |

## 5. `main.cpp` から見た実行フロー

### 5.1 初期化フェーズ

`main.cpp` の先頭では、計算開始前に次の処理が順番に実行される。

1. `solverConfig cfg; cfg.read("solverConfig.yaml");`
   - 計算条件、時間積分条件、物性値、出力条件などを YAML から読み込む
2. `mesh msh; msh.readMesh(cfg.meshFileName);`
   - HDF5 形式のメッシュを読み込む
3. `matrix mat_ns; mat_ns.initMatrix(msh);`
   - 疎行列構造を初期化する
4. `readBcondConfig(cfg, msh.bconds);`
   - `bcondConfig.yaml` を読み、境界面に対して境界条件種別とパラメータを対応づける
5. `variables var; var.allocVariables(cfg.gpu, msh);`
   - セル値・面値・デバイス側バッファを確保する
6. `var.readValueHDF5(cfg.valueFileName, msh);`
   - 初期値を HDF5 から読み込む
7. `cudaConfig cuda_cfg = cudaConfig(msh); msh.setMeshMap_d();`
   - GPU 側で使うメッシュ接続マップを構築する
8. `msh.setPeriodicPartner();`
   - 周期境界の対応を設定する
9. `var.setStructuralVariables(cfg, cuda_cfg, msh);`
   - 体積、重心などの構造量を変数に格納する
10. `dependentVariables(...)`
    - 保存変数から圧力、温度、速度、全エンタルピーなどを更新する
11. `gasProperties_d_wrapper(...)`
    - GPU 側で比熱比や比熱などの物性値を更新する
12. `fluct_variables` の初期化
    - 変動流入境界用の補助量を確保・初期化する
13. `applyBconds(...)`
    - 境界条件を適用する
14. `calcGradient_d_wrapper(...)`
    - 初期状態の勾配を GPU 側で計算する
15. `updateVariablesOuter(...)`
    - 必要な値を GPU から CPU 側へ戻す
16. `setDT_d_wrapper(...)`
    - 初期の時間刻みを設定する
17. `point_probes::init(...)`
    - 点プローブを初期化する
18. `outputH5_XDMF(...)`
    - 初期状態を出力する

ここで重要なのは、初期化の段階ですでに「従属変数更新 → 境界条件適用 → 勾配計算」という、時間発展ループ内と似た並びが一度走ることである。

### 5.2 時間発展ループ

時間発展は `cfg.mainLoopCount()` で決まる主ループで回る。明示法では `nStepOuter` のループ、その内側に `cfg.nStage` 回のステージループが入る。陰解法・非定常では `nStepOuter` の外側ループと `nStepInner` の内側反復、陰解法・定常では `nStepInner` を主ループとして回す。1 回の更新で行われる処理順序は次の通りである。

1. `updateVariablesInner(...)`
   - ステージ計算の前に必要な値を GPU 側で使える形にそろえる
2. `dependentVariables(...)`
   - 保存変数から圧力、温度、速度などを再計算する
3. `gasProperties_d_wrapper(...)`
   - 物性値を GPU 側で更新する
4. `applyBconds(...)`
   - 境界条件を適用する
5. `calcGradient_d_wrapper(...)`
   - 勾配を計算する
6. `limiter_d_wrapper(...)`
   - リミタを計算する
7. `ducrosSensor_d_wrapper(...)`
   - センサ量を計算する
8. `turbulent_viscosity_d_wrapper(...)`
   - 乱流粘性を更新する
9. `convectiveFlux_d_wrapper(...)`
   - 対流流束を計算する
10. `viscousFlux_d_wrapper(...)`
    - 粘性流束を計算する
11. `implicitCorrection(...)`
    - 必要な補正を加える
12. `timeIntegration_d_wrapper(...)`
    - ステージの時間積分を進める

この順序は文書を読む上で最も重要である。多くの拡張は、この 12 個のどこに新しい処理を差し込むかという形で考えると整理しやすい。

### 5.3 ステップ後処理

各ステップの最後では次が実行される。

1. `updateVariablesOuter(...)`
2. `outputH5_XDMF(...)`
3. `outputBconds_H5_XDMF(...)`
4. `pprobes.outputProbes(...)`
5. `setDT_d_wrapper(...)`
6. `cfg.totalTime += cfg.dt;`

すなわち、数値計算の中心は GPU に置きつつ、出力と時間管理は CPU 側でまとめて制御している。

## 6. 主要データ構造

### 6.1 `mesh`

`mesh/mesh.hpp` の `mesh` 構造体は、節点 `node`、面 `plane`、セル `cell`、境界条件 `bcond` の集合を保持する。

- `nodes`: 節点座標と接続情報
- `planes`: 面の節点、隣接セル、面法線ベクトル、面積、面重心
- `cells`: セルを構成する節点、面、体積、セル重心、要素型
- `bconds`: 境界条件ごとの面集合、境界値、境界補助データ

さらに GPU 用の接続配列として、`map_plane_cells_d`、`map_cell_planes_index_d`、`map_cell_planes_d` などを持つ。これらは `setMeshMap_d()` によって準備される。

### 6.2 `variables`

`variables.hpp` の `variables` クラスは、計算中に使う流体変数を一元管理する。

- `c`: ホスト側のセル変数
- `p`: ホスト側の面変数
- `c_d`: デバイス側のセル変数
- `p_d`: デバイス側の面変数

変数名は文字列で管理されており、たとえば次が含まれる。

- 保存変数: `ro`, `roUx`, `roUy`, `roUz`, `roe`
- 従属変数: `Ux`, `Uy`, `Uz`, `T`, `P`, `Ht`, `sonic`
- 勾配: `dUxdx`, `dPdx`, `drodx` など
- 数値安定化用: `limiter_*`, `ducros`, `dt_local`
- 乱流・物性: `vis_lam`, `vis_turb`, `wall_dist`, `gamma`, `cp`

`output_cellValNames` に含まれる変数だけが標準出力の対象になる。

### 6.3 `bcond`

`bcond` は境界条件 1 種類分の情報を持つ。

- 境界に属する面とセルのインデックス
- `inputInts`、`inputFloats` による設定値
- `bvar` による境界面上の実値
- `bint` による周期境界などの整数補助情報
- GPU 側の境界変数バッファ `bvar_d`、`bint_d`

境界条件は単に「名前だけ」を持つのではなく、対象面集合、入力値、出力対象フラグまでまとめて 1 つのオブジェクトにして扱っている。

## 7. 設定ファイルの読み方

### 7.1 `solverConfig.yaml`

`input/solverConfig.cpp` の `solverConfig::read()` が読み込む。主に次の情報を担当する。

- メッシュ形式とファイル名
- 初期値ファイル名
- ソルバ種別
- GPU 使用有無
- 時間積分条件
- CFL や `dt` などの時間刻み制御
- 空間離散化設定
- 乱流モデル設定
- 圧縮性、粘性、熱物性関連の設定
- 出力間隔

`getValidatedValue()` を通して読み込んでいるため、キーの欠落や型不整合は比較的早い段階で検出される。

### 7.2 `bcondConfig.yaml`

`boundaryCond.cpp` の `readBcondConfig()` が読み込む。各境界名に対して次を受け取る。

- `physID`
- `kind`
- `ints`
- `floats`
- `outputHDFflg`

読み込んだ結果は `physID` をキーに `msh.bconds` に対応づけられる。つまり、メッシュ側の物理境界 ID と YAML 設定がこの段階で結び付けられる。

## 8. 出力系の構成

### 8.1 メイン出力

`output/output.cpp` の `outputH5_XDMF()` が、セル中心量を HDF5 と XDMF で書き出す。

- HDF5 側にはメッシュ座標 `MESH/COORD`、接続 `MESH/CONNE`、各変数 `VALUE/...` を書き込む
- XDMF 側には ParaView などから読めるトポロジ、座標、属性の定義を書く

GPU 計算中でも、出力前には必要なセル変数を `copyVariables_cell_D2H()` でホスト側に戻してから書き出している。

### 8.2 境界出力

`outputBconds_H5_XDMF()` は、`outputHDFflg == 1` の境界条件だけを個別ファイルとして出力する。壁面量や流入境界量を後から確認したいときに使う。

### 8.3 プローブ出力

`probe/` の `point_probes` は、指定点に対応する情報を初期化し、各ステップで `outputProbes()` により時系列出力を行う。

### 8.4 残差モニタと host-device 同期

`residual_history.csv` への残差 RMS 出力は `main.cpp` の `ResidualCsvLogger` が担い、reduction 本体は
`cuda_forge/residualMonitor_d.cu` にある。GPU 経路では残差二乗和を **fused 1 カーネルで全保存量一括に取り、
device バッファに常駐**させる (`DeviceResidualReducer` / `residualSumSq_d`)。host へは `monitorInterval`
ステップごとに 1 回だけまとめて D2H 転送し CSV へ書き出す。これにより毎ステップの値は保ったまま、
変数ごとの `thrust::transform_reduce` (host スカラ返り = `cudaStreamSynchronize`) による per-step 同期を除く。
`monitorInterval` は console のモニタ行 (§8.5) の出力頻度も兼ねる (モニタリング出力の共通間隔)。
dt 適応 (`setDT` の `thrust::max_element`→`cfg.dt`) はこれとは独立に制御する (定常では間引き可・解に不影響、
explicit/dual-time は必要に応じ毎ステップ適応; 表示のみ `monitorInterval` で間引く)。

> 背景: 定常 implicit (timeIntegration=11) は 1 step の GPU 実働が wall の約半分しかなく、残り半分は
> launch/host 同期だった。毎 step の残差 RMS と detectNaN がその主因。`detectNaN` も fused な device int
> フラグ化し、`detectNaNInterval` ステップごとにのみフラグを host 読み出しする。`monitorInterval` /
> `detectNaNInterval` の既定は 1 で従来挙動。設計詳細は
> [`.github/plans/architecture-residual-monitor-async.md`](../../plans/accepted/architecture-residual-monitor-async.md)。

GPU 計算結果をホストで使う一般則は変わらない: 出力・監視のために値を host で読む箇所は同期点になるため、
頻度を必要最小限にするか device 常駐＋まとめ転送にするのが基本方針。

### 8.5 console モニタ行 (stdout / `forge_run.log`)

`main.cpp` の `StepMonitor` が、`monitorInterval` ステップごとに **1 行**の要約を stdout へ出す。残差の正本は
`residual_history.csv` (§8.4) で変わらないが、`tail -f forge_run.log` だけで「速いか・収束しているか・壊れ始めて
いないか」を読めるようにするのが目的。旧形式 (`Step :` / `max cfl :` / `dt :` の 3〜4 行/step と RK/inner の
`Stage : n` 行) は廃止した。

**出力する量はモードで変える**。判定軸は `time.unsteady` であり、`isImplicit` ではない (下表)。

| モード | `unsteady` | `isImplicit` | 時間刻みの意味 | モニタ行に出す時間量 |
| --- | --- | --- | --- | --- |
| 定常 陰解法 (block-DPLUR) | 0 | 1 | `dt_local = cfl_pseudo·V/λ` (局所擬似時間)。`cfg.dt` は打ち消されて無意味 | 出さない (`cfl_pseudo`/`implicitRelax` を起動時に 1 回) |
| 定常 陽解法 (局所 dt) | 0 | 0 | 同上 (`setDTlocal_pseudo_cell_d`) | 出さない (同上) |
| 非定常 陽解法 (RK) | 1 | 0 | `cfg.dt` が一様物理 dt (dtControl=1 なら CFL 適応) | `t`, `dt`, `maxCFL` (物理 CFL) |
| 非定常 dual-time 陰解法 | 1 | 1 | `cfg.dt` が固定物理 dt、サブ反復は `cfl_pseudo` | `t`, `dt`, `maxCFL` (物理 CFL; 経験的に ≲12 が安定域) |

`max cfl` は `setDT_d_wrapper` が host 読みした値を `cfg.monitorCflMax` に格納するだけにし、印字はモニタ行が担う。
CFL を評価した刻みは `cfg.monitorCflDt` に対で保存する: `dtControl: 1` (CFL 適応) の陽解法では `setDT` が直後に
`cfg.dt` を次 step 用に更新するので、モニタ行は「この step の前進に使った `dt`」と「その dt での `maxCFL`」を対で出し、
更新後の値は `dt_next` として別に出す。定常モードでは `printCfl=false` を渡すので**表示目的**の host 読みは無くなる
(`dtControl: 1` のときの dt 適応目的の読み (`adaptDt`) は従来どおり残る。定常では適応結果は解に効かない)。

モニタ行の構成 (1 行, 空白区切り, `|` で区分):

```
step     200 | 3.96 ms/step elapsed 0.8 s eta 0 s | rms_ro 7.14e-06 (-0.3) worst rms_roOmega 6.16e+03 (+2.4)
step     200 | t 5.0000e-05 dt 2.50e-07 maxCFL 0.07 | 4.06 ms/step elapsed 0.6 s eta 0 s | rms_ro 8.16e-03 (-0.1) worst rms_roUx 6.43e+00 (-0.0) from0 rms_roUy 3.7e-06
```
(1 行目: 定常陰解法 case/16 run_0330、2 行目: 非定常陽解法 case/05 run_0013。`(+2.4)` は増大 = 初期比 10^2.4 倍)

- **時間量** (unsteady のみ): 物理時間 `t`、物理 `dt`、`maxCFL`。
- **速度**: 直前 monitor 区間の壁時計 (`steady_clock`) 平均 `ms/step`、開始からの経過秒、`nStepOuter` までの予想残り秒。
  終了時の `Time = ... s` も CPU 時間 (`clock()`) から壁時計へ変えた (GPU 待ちが含まれる値になる)。
- **残差**: その step の `outer_end` 行から、`rms_ro` と「初期値 (step 0 `outer_begin`) からの低下桁数が最小の列」
  (= worst) を値と `log10(現在/初期)` (負 = 低下, 正 = 増大) で示す。初期値が 0 の列は桁数が定義できないので、
  現在非 0 なら最大のものを `from0` として別表示する。全列の詳細は CSV。非有限値があれば末尾に `*** NaN ***` を付ける。
  値は `ResidualCsvLogger` が flush 時に保持した最新 `outer_end` 行 (`latestSummary()`) を使うので、CSV と同じ数値・
  追加の reduction 無し。

同期コストは monitor ステップでの CSV flush (D2H 1 回) と、unsteady の `max cfl` 読み (従来どおり) だけ。
`monitorInterval: 1` (既定) なら毎ステップ 1 行。

## 9. CUDA 側実装の見取り図

`cuda_forge/` には、GPU 側の主要計算がまとめられている。役割ごとに見ると次のように整理できる。

### 9.1 構造量・従属変数・更新

- `calcStructualVariables_d.*`
- `dependentVariables_d.*`
- `update_d.*`

保存変数から従属変数を組み立てたり、GPU 側の更新を補助したりする層である。

### 9.2 境界条件

- `boundaryCond_d.*`
- `fluct_variables_d.*`

境界面やゴーストセルに対して、GPU 側で境界条件を適用する。

### 9.3 空間離散化と安定化

- `calcGradient_d.*`
- `limiter_d.*`
- `ducrosSensor_d.*`

勾配、リミタ、センサ量を計算し、後段の流束計算に必要な量を準備する。

### 9.4 流束と乱流・物性

- `convectiveFlux_d.*`
- `viscousFlux_d.*`
- `turbulent_viscosity_d.*`
- `gasProperties_d.*`

対流流束、粘性流束、乱流粘性、物性値の更新を担当する。

### 9.5 時間積分

- `implicitCorrection_d.*`
- `timeIntegration_d.*`
- `setDT_d.*`

ステージ更新、補正、時間刻み計算を担当する。

### 9.6 GPU 化の土台

GPU 側の計算を読むときは、個別カーネルだけでなく次の 2 点も重要である。

1. `variables::copyVariables_*`
   - ホスト側とデバイス側の変数転送を受け持つ
2. `cudaConfig` と `mesh::setMeshMap_d()`
   - メッシュ接続情報を GPU 側で引きやすい形へ展開する

この 2 つがあるため、`main.cpp` 側は比較的単純な wrapper 呼び出しで GPU カーネル群を扱えるようになっている。

## 10. ビルド対象と補助実行ファイル

`CMakeLists.txt` では、ソルバ本体以外にも次の補助実行ファイルをビルドしている。

- `forge`
  - メインの解析実行ファイル
- `convertGmshToForge`
  - Gmsh メッシュを forge 用 HDF5 形式へ変換する
- `mesh_part`
  - METIS を使ってメッシュ分割を行う

通常の解析実行では `forge` を直接見るのが中心になるが、入力データの準備まで含めて理解したい場合は `convertGmshToForge` もあわせて読むと全体像がつかみやすい。

## 11. 読み進める順番の推奨

初めてコードを読む場合は、次の順で読むと迷いにくい。

1. `main.cpp`
2. `input/solverConfig.*`
3. `mesh/mesh.*`
4. `variables.*`
5. `boundaryCond.*`
6. `dependentVariables.*`
7. `output/output.*`
8. `cuda_forge/` の wrapper と主要カーネル

最初から個別の CUDA カーネルに入るより、まず CPU 側の呼び出し順とデータ構造を押さえてから GPU 側に降りた方が理解しやすい。

## 12. 拡張ポイント

### 12.1 新しい流束法を追加したい場合

最初に見るべき場所は `convectiveFlux.*` と `cuda_forge/convectiveFlux_d.*` である。設定値 `convMethod` がどこで分岐に使われているかを追うと、追加位置を決めやすい。

### 12.2 新しい境界条件を追加したい場合

`boundaryCond.cpp` と `cuda_forge/boundaryCond_d.*` を対応で見る。YAML 側の `kind`、必要な `ints`/`floats`、`bcond` への格納方法、GPU 適用部の順で追うと整理しやすい。

### 12.3 新しい出力項目を追加したい場合

`variables.hpp` の `output_cellValNames` と `output/output.cpp` の出力処理を確認する。GPU 計算結果を出したい場合は、出力前にホスト側へ戻す必要があるかどうかもあわせて確認する。

## 13. まとめ

`solver_density_cuda` の構成は、大きく分けると「CPU 側の制御と I/O」と「GPU 側の数値計算」に分かれている。`main.cpp` はその接続点であり、`mesh` と `variables` がデータの土台、`boundaryCond` と `output` が周辺機能、`cuda_forge/` が数値計算の本体を担っている。

実装を追うときは、まず `main.cpp` の処理順を固定し、その後で各モジュールの責務とデータの持ち方を確認するのが最も効率的である。