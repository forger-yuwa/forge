# クラウド GPU 環境構築 (Docker + 計算投入 + 速度計測)

## メタ

- **area**: `その他 (tooling / infra)`
- **status**: `in_progress`
- **related_docs**:
  - [`procedures/development-environment.md`](../../procedures/development-environment.md) (現行の Docker / native 使い分けルール)
  - [`solver_density_cuda/Dockerfile.cuda.dev`](../../solver_density_cuda/Dockerfile.cuda.dev) (現行の開発イメージ)
- **related_plans**:
  - (将来) `architecture-mpi-multigpu` — MPI 化の実装計画は別 plan として切り出す。本計画はその**実行基盤**を先に用意する
- **created**: `2026-07-20`
- **owner**: `sano`

## 1. 目的

ローカル (WSL2 + RTX 3060) に限られている forge の実行環境をクラウド GPU に拡張し、
(1) Docker で再現可能な環境構築、(2) 計算ケースの投入・回収、(3) 安定した速度計測 (プロファイル含む) をクラウド上で行えるようにする。
将来の速度最適化の試行錯誤と MPI 化 (single-node multi-GPU → multi-node) の受け皿となるプラットフォームを確定させる。

## 2. スコープ

- **やる**:
  - クラウドサービスの選定と、その根拠の記録 (本計画 §4)
  - 単一 GPU インスタンスでの環境構築手順の確立 (Docker build / ケース実行 / 速度計測 / ncu プロファイル)
  - メッシング・可視化のクラウド完結 (ヘッドレス gmsh 同梱 / pvserver + SSH トンネル / 必要時のみ NICE DCV)
  - 結果成果物 (`run_*`) の回収経路 (S3) とコスト制御 (spot / auto-stop / billing alert)
  - multi-GPU インスタンス (4 GPU) の起動確認まで (MPI 実装の前提確認)
  - 手順の `procedures/` への文書化
- **やらない**:
  - MPI 化そのもの (ソルバのドメイン分割・halo 交換) → 別 plan `architecture-mpi-multigpu`
  - multi-node (EFA / ParallelCluster) の構築 → MPI single-node が動いてから拡張 (§4 に方向性のみ記す)
  - CI 常駐化・自動回帰実行 (必要になったら別 plan)

## 3. 関連 docs と前提

- 現行ルール ([`procedures/development-environment.md`](../../procedures/development-environment.md)):
  - 通常のビルド・ケース実行は Docker 基準 (`Dockerfile.cuda.dev`)
  - **速度評価・プロファイルは native 基準** — これは WSL2 の Docker 層が計測を不安定にするための規定。クラウドの Linux native ホストではこの問題は軽く、`ncu` / `nsys` (host importer 含む) がフルに動くのがクラウド化の利点の一つ
- リポジトリは GitHub (`forger-yuwa/forge`, private) にあり、クラウド側は deploy key (read-only) で clone できる
- `case/` 配下の run 成果物は commit しない運用のため、クラウド上の結果はリポジトリ外 (S3) で回収する
- ローカル GPU は RTX 3060 (Ampere, CC 8.6)。ビルドは `FORGE_CUDA_ARCHITECTURES` で対象アーキテクチャを指定する

## 4. 設計方針 (サービス選定と構成)

### 4.1 サービス選定: AWS EC2 (spot) を採用

要件は「① 自前 Docker が動く root 付き VM ② `ncu`/`nsys` が動く (= GPU パフォーマンスカウンタへのアクセス) ③ multi-GPU → 将来 multi-node MPI への一本道 ④ 使う時だけ課金」。この 4 点で比較した。

| サービス | 形態 | ①Docker | ②ncu/nsys | ③MPI への道 | ④コスト | 評価 |
| --- | --- | --- | --- | --- | --- | --- |
| **AWS EC2 (G 系)** | VM (root) | ○ | ○ (root で可) | ○ 4GPU インスタンス → EFA/ParallelCluster まで同一基盤 | spot で on-demand 比 3〜7 割引 | **採用** |
| GCP (g2, L4) | VM (root) | ○ | ○ | △ multi-node HPC は AWS より手薄 | 同等 | 対抗馬。AWS 経験・EFA の分だけ AWS 優先 |
| Lambda Labs | VM (root) | ○ | ○ | △ 1x/8x に二極化し中間帯 (4GPU) が薄い、spot 無し | A10 $1.29/h 程度 | 単発利用なら簡単だが伸びしろ不足 |
| RunPod / Vast.ai | **コンテナ貸し** | × (docker-in-docker 不可) | × (権限制約) | × | 最安 (RTX4090 ~$0.3-0.6/h) | 本目的には**不適** (自前 Dockerfile・プロファイル・MPI がいずれも制約) |
| Azure (NC/NV) | VM (root) | ○ | ○ | ○ | 高価 | コスト劣位 |

RunPod/Vast 系は価格最安だが「Docker で環境構築」自体ができない (借りるのがコンテナ) ため除外。VM 系の中で、spot・4GPU 中間帯インスタンス・EFA (multi-node HPC ネットワーク) を全部持つ AWS を採用する。

### 4.2 インスタンスの使い分け (段階)

| 段階 | インスタンス | GPU | CC | 参考価格 (us-east-1, on-demand)* | 用途 |
| --- | --- | --- | --- | --- | --- |
| P1: 環境構築・回帰・速度計測基準 | `g5.xlarge` | 1×A10G 24GB | **8.6** | ~$1.01/h | 入口。**ローカル RTX 3060 と同じ Ampere/CC 8.6** でカーネル特性を直接比較できる基準機。環境構築〜スモーク回帰〜ベンチ基準まで一気にここで行う (ユーザーは AWS 経験ありのため最安機 g4dn での手順練習フェーズは省略、2026-07 決定) |
| P1': 新アーキ試験 (任意) | `g6.xlarge` | 1×L4 24GB | 8.9 | ~$0.80/h | Ada 世代での挙動確認 (g5 より安い) |
| P2: MPI single-node | `g6.12xlarge` (または `g4dn.12xlarge` で安価に) | 4×L4 / 4×T4 | 8.9 / 7.5 | ~$4.6/h / ~$3.9/h | MPI 化のデバッグ・スケーリング測定 |
| P3: multi-node (スコープ外) | EFA 対応サイズ + ParallelCluster | — | — | — | MPI plan 側で扱う |

\* 2026-07 時点の概数。spot は概ね 3〜7 割引。投入前に [instances.vantage.sh](https://instances.vantage.sh/) で実勢を確認する。**リージョンは `us-east-1` に決定** (spot 在庫・価格が最良、2026-07 決定。東京は 2〜3 割高)。

### 4.3 構成要素

- **AMI**: AWS Deep Learning Base GPU AMI (Ubuntu 22.04) — NVIDIA driver / Docker / nvidia-container-toolkit 導入済み。素の Ubuntu + 手動導入はしない (手順が増えるだけ)
- **ストレージ**: EBS gp3 200GB (インスタンス stop 中も保持、課金は ~$0.08/GB·月)。長期成果物は S3 バケット `forge-runs` (仮) に `aws s3 sync`
- **リポジトリ取得**: GitHub deploy key (read-only) で clone。`.git` が ~900MB あるので `git clone --depth 1 --single-branch` を既定にする
- **イメージ**: インスタンス上で **クラウド用イメージ `Dockerfile.cuda.cloud`** を直接 build (単一ノードなら ECR 不要)。dev イメージから ParaView GUI/Qt/X サーバ層を除き、**gmsh はヘッドレス (バイナリ + python API) で同梱**してメッシュ生成をクラウドで完結させる。multi-node 化の際に ECR へ移行
- **可視化**: 定型 PNG 後処理はイメージ内 matplotlib で完結。対話 3D は **pvserver (egl 版) + SSH トンネル** でローカル ParaView クライアントから接続 (`res_*.h5` を転送せず GPU 側レンダリング)。gmsh GUI 等のフルデスクトップが要るときのみ NICE DCV (EC2 上ライセンス無料) を立てる
- **ビルドアーキテクチャ**: T4=`FORGE_CUDA_ARCHITECTURES=75`, A10G=`86`, L4=`89`
- **速度計測**: リポジトリのルール通り native (ホスト直) ビルドで行う。`tools/build_native_wsl.sh` を Linux native でも使えるよう汎用化 (or `build_native_linux.sh` 追加)。`ncu` は root 実行または `NVreg_RestrictProfilingToAdminUsers=0` を設定。WSL で不能だった `nsys` の `.nsys-rep` 変換・`nsys stats` がクラウドでは完結する
- **コスト制御**:
  - 通常は on-demand + **stop 運用** (使い終わったら stop、EBS のみ課金)。長時間バッチは spot (中断は `res_*.h5` restart で許容)
  - GPU 使用率 0% が 30 分続いたら self-shutdown する cron をインスタンスに仕込む (消し忘れ保険)
  - AWS Budgets で月額アラート **$50** を設定 (2026-07 決定。g5.xlarge on-demand で月 ~50h 相当)
- **セキュリティ**: inbound は SSH (自宅 IP 制限) のみ、または SSM Session Manager でポーレス運用

## 5. 実装ステップ

1. **AWS 側準備 (手作業)**: 既存アカウントを使用。`us-east-1` で G 系 vCPU quota を on-demand / spot とも 48 vCPU へ引き上げ申請 → **結果: 両方とも 8 vCPU で承認 (2026-07-22)**。P1 (`g5.xlarge`=4 vCPU) には十分。**P2 の `*.12xlarge` (48 vCPU) は現枠では不可**のため、P2 着手前にケース再開でユースケース (MPI 開発で 4 GPU 機を短時間利用) を添えて再申請する。ほか S3 バケット、Budgets アラート ($50/月)
2. **ブートストラップ整備**: `solver_density_cuda/tools/cloud/` を新設し、(a) launch 用 user-data スクリプト (clone → docker build)、(b) `Dockerfile.cuda.cloud`、(c) `sync_results_s3.sh`、(d) idle auto-stop cron を置く
3. **native ビルドスクリプトの汎用化**: `build_native_wsl.sh` を確認した結果、WSL 固有処理は無く Linux native でそのまま動く (2026-07-20 確認、変更不要。エラーメッセージの文言に WSL とあるのみ)
4. **P1 スモーク回帰**: `g5.xlarge` で Docker build 後、[`procedures/verification/README.md`](../../procedures/verification/README.md) の標準に沿って次を実行する (詳細は §6)。run は新規 `run_*` ディレクトリ + 各 case README の run 一覧表更新 (通常ルール通り)
5. **P1 ベンチ基準確立**: 同じ `g5.xlarge` (CC 8.6) で native ビルド + `FORGE_PROFILE=1` 計測、同一ケースを 3 回実行して再現性 (ばらつき) を確認。ローカル RTX 3060 との比較表を作り、以後の速度最適化のベースラインとする。`ncu` / `nsys stats` の動作確認もここで行う
6. **P2 前提確認**: 4 GPU インスタンスを短時間起動し、`nvidia-smi` で 4 GPU 認識・CUDA-aware OpenMPI 入りイメージのビルドまで確認 (MPI 実装は別 plan)
7. **文書化**: `procedures/cloud-aws-gpu.md` を新規作成 (起動→実行→回収→停止の一連手順)、`procedures/development-environment.md` と `procedures/README.md` からリンク

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | 4GPU 起動確認 (§5 ステップ 6) | `*.12xlarge` を短時間起動し 4 GPU 認識と CUDA-aware OpenMPI 入りイメージのビルドまで。**現 quota (8 vCPU) では起動できないため再申請が前提**。これが済めば本計画は accepted へ移せる |
| 2 | `setup_instance.sh` の実機再確認 | 2026-09-06 の Docker 多段化に伴い `docker_build.sh cloud` 呼び出しへ変更済み。P1 実機で 1 回通す |

## 6. 検証

- **単体 / ビルド**: EC2 上で Docker build 成功、native build 成功 (対象 CC で)
- **検証ケース (P0 スモーク回帰)**: [`procedures/verification/README.md`](../../procedures/verification/README.md) の標準構成に従い、順に:
  1. `case/05.sod_shock_tube` — 最短の 1D スモーク (厳密解あり)。環境が根本的に壊れていないかの初手確認
  2. `case/20.naca_ml/001.test/run_slau` — 既定の検証先。メッシュ生成 → `forge` 実行のパイプライン一式がクラウドで通ることの確認
  3. `case/08.bump/verify/run_verification.sh` — **定量回帰の本命**。低マッハ/高マッハ × 陽解法/陰解法を流し、y≈0.3 line profile の一致と収束度合いの非劣化を自動判定する既存スイート
  - あわせてクラウド完結ワークフローの確認: naca_ml のメッシュ生成をコンテナ内ヘッドレス gmsh で実施 (2 の run_case.sh が自然に通ることで確認)、pvserver + SSH トンネルでローカル ParaView からの接続・描画を 1 回確認
  - 1〜3 は軽いケースなので **cell / node 両方** (`discretization` 2 run) で実行する (軽い計算は両系統・LES/DES 級の重いケースは node のみ、という運用方針に従う)。node 未対応構成はその旨明記のうえ cell のみ
- **判定基準**:
  - 回帰: 各 run で `check_convergence.py` PASS (VERDICT 添付) + `run_verification.sh` の自動判定 PASS。ローカル (RTX 3060) の同設定結果との突合は、**node は run-to-run ほぼ決定的 (~1e-6) なので厳しく、cell は atomicAdd 非決定性のノイズ床 (~6e-4) を許容**して判断する
  - 速度計測: 同一 run 3 回の total 時間ばらつきが数 % 以内 (計測基準環境として成立)
  - コスト: P1 の作業合計が想定内 (目安 $20〜30) に収まることを記録

## 7. 影響範囲

- 追加: `solver_density_cuda/tools/cloud/` (スクリプト群), `solver_density_cuda/Dockerfile.cuda.cloud`, `procedures/cloud-aws-gpu.md`
- 変更: `procedures/development-environment.md` / `procedures/README.md` (リンク追記)。`build_native_wsl.sh` は変更不要と確認済み
- ソルバ本体 (`solver_density_cuda/cuda_forge/`) には触れない

## 8. 完了条件

- [ ] AWS 上で Docker build + 回帰ケース (cell/node) がローカルと整合
- [ ] `g5.xlarge` で速度計測ベースライン (対 RTX 3060 比較表 + `ncu` 動作) を記録
- [ ] 4 GPU インスタンスの起動・認識確認 (MPI plan への引き継ぎ条件)
- [ ] `procedures/cloud-aws-gpu.md` 整備、コスト制御 (auto-stop / Budgets) 稼働
- [ ] 本計画の `status` を `done` にし `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-07-20` — 初稿。サービス比較 (AWS/GCP/Lambda/RunPod/Azure) の上で AWS EC2 spot を採用、段階構成を策定。
- `2026-07-20` — ユーザー決定を反映: リージョン `us-east-1` / 入口インスタンスは `g5.xlarge` (AWS 経験ありのため g4dn での P0 練習フェーズを省略し P1 に統合) / Budgets $50/月。quota 申請は on-demand/spot とも 48 vCPU を一括で出す。回帰は「軽いケースは cell/node 両方、LES/DES は node のみ」の運用に確定。
- `2026-07-20` — §5 ステップ 2〜3 実装: `Dockerfile.cuda.cloud` (計算専用軽量イメージ)、`tools/cloud/{setup_instance,idle_autostop,sync_results_s3}.sh`、`procedures/cloud-aws-gpu.md` (起動→実行→計測→回収→停止の手順) を追加、`procedures/README.md` / `development-environment.md` にリンク。`build_native_wsl.sh` は WSL 依存無しでそのまま Linux native 可と確認 (変更不要)。quota 申請 (48 vCPU × on-demand/spot) はユーザーが提出済み。残: P1 実機検証 (§5 ステップ 4〜6)。
- `2026-07-20` — スコープ追加 (ユーザー要望): メッシング・可視化もクラウド完結に。`Dockerfile.cuda.cloud` にヘッドレス gmsh を同梱、対話 3D は pvserver + SSH トンネル (クライアント/サーバ同一版必須、egl 版)、フル GUI は NICE DCV を任意手順として `procedures/cloud-aws-gpu.md` に追記。
- `2026-08-30` — **P1 実機立ち上げ + sod スモーク完了** (`g5.xlarge`, us-east-1, Deep Learning Base AMI Ubuntu 24.04, A10G driver 595.91.07)。立ち上げで直した問題 3 件: ① `setup_instance.sh` の GitHub 認証チェックが pipefail 下で常に失敗 (GitHub は成功でも ssh exit 1) → 出力受けに修正。② gmsh.info ダウンロード不通 (サイト障害) → gmsh は PyPI wheel のみに変更。③ コンテナ内ビルドが `CMAKE_CUDA_ARCHITECTURES=52` 既定で double atomicAdd コンパイル不能 (`project(CUDA)` が検出時に定義するため CMakeLists のガード無効) → ホスト GPU CC を `CUDAARCHS` で明示。既知事象 2 件: `convertGmshToForge` が終了時 `cudaFree` GPUassert (invalid argument) で exit≠0 (出力 h5 は完全、driver CUDA 13.2 × コンテナ 12.4 で顕在化とみられる。要フォロー)。組み込み `initial:"sod"` は非次元で使用不可 (既知バグ再確認) → `case/05.sod_shock_tube/gen_sod3d_ic.py` (node/cell 両対応の物理 IC 生成) を追加。**検証結果**: sod 3D periodic を node/cell 両方でクラウド完結パイプライン (コンテナ gmsh → convert → 実行) により完走、NaN なし・mesh 品質 PASS。同一コード・同一入力のローカル RTX 3060 と**最大相対差 ~2e-6 / L2 ~1e-7 で一致** (run_0009/0010=EC2, run_0011/0012=ローカル再現)。残: naca_ml (§6-2)・bump 回帰スイート (§6-3)・速度ベースライン + ncu/nsys (§5-5)・4GPU 起動確認 (§5-6, quota 再申請待ち)。
- `2026-08-30` — **P1 検証メニュー完了 (naca_ml / bump / 速度・プロファイル)**。
  - **naca_ml スモーク (§6-2) PASS**: コンテナ gmsh → convert → 4000 step 実行。mesh 品質 PASS。同一入力のローカル再現 (`case/20.naca_ml/001.test/run_aws_p1_repro`) と最大相対差 ~2.5e-5 / L2 ~1e-6 (EC2 側 run は `run_aws_p1_slau`)。check_convergence は両環境とも同一表示の「主要3成分 falling (収束途中)・2D の roUz 機械ゼロ床がプラトー扱い」。
  - **bump 回帰 (§6-3)**: loM exp/imp は **PASS** (ΔP≤0.2%・ΔM≤0.43%・収束フロア=ベースライン同等)。hiM exp/imp は committed baseline (2026-06-11, ローカル作成) 比で FAIL 表示 (ΔP 相対 10.4%) だが、**真因は GPU 個体系統差**と切り分け完了: 差は衝撃波下流 x≈0.83–0.91 の 8 点のみ (ΔP/Pt 絶対 0.0197)、ローカルは baseline と 1e-4 一致、**EC2 の run-to-run も 1.1e-4 で再現**、CUDA 12.4 コンテナ版/13.2 native 版・陽/陰でも同値 → cell 離散化の atomicAdd 集約順が GPU (A10G vs RTX 3060, 同 sm_86) で異なり、臨界セル上の衝撃波ロック位置が 1 セルずれる。既知の「cell cross-run ~0.5% / node 決定的」と整合し **AWS 環境の欠陥ではない**。cell 衝撃ケースの跨マシン回帰にはマシン別 baseline か衝撃許容メトリクスが要る (node 主体方針なら実害小)。残差は 2.0–2.8 桁低下後プラトーで baseline_convergence.csv のフロアと同一 (スイートの設計通り)。
  - **速度ベースライン (§5-5)**: bump hiM imp 3000 step ×3、native ビルド・`FORGE_PROFILE=1`。**A10G 6.75/6.74/6.79 s (ばらつき 0.7% → 計測基準として成立)** vs RTX 3060 (WSL native) 15.51/14.12/15.64 s (ばらつき ~10%) → **A10G ≈2.2 倍速**。`ncu` OK (要 sudo)。`nsys` OK (スタンドアロン版導入・手順は procedures 参照)。カーネル内訳 (nsys, hiM imp): block-DPLUR 59% / limiter 10% / SLAU 8.8%。
  - 立ち上げで追加した修正: idle_autostop に CPU load 条件 (native ビルド中の誤 shutdown 実績への対処)、bump loM 正準 run の probe.yaml を追跡化、native ビルドの CUDA13/CCCL・boost・pip h5py 対応 (procedures に反映)。
  - **残タスク: 4GPU 起動確認 (§5-6) のみ** (ユーザー側 quota 状況の確認待ち)。完了後に本計画を accepted へ移動する。
- `2026-09-06` — **クラウド用イメージの所在が変わった**。`Dockerfile.cuda.cloud` は廃止し、多段
  [`Dockerfile.cuda`](../../solver_density_cuda/Dockerfile.cuda) の `cloud` ステージに統合した。
  ビルドは `solver_density_cuda/tools/docker_build.sh cloud` で、`tools/cloud/setup_instance.sh`
  も同様に更新済み (実機での再確認は未実施)。経緯と設計判断は
  [tooling-docker-image-split.md](../accepted/tooling-docker-image-split.md) を参照。
