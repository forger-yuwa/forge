# AWS GPU インスタンスでの forge 実行手順

クラウド (AWS EC2, `us-east-1`) で forge を Docker ビルド・計算投入・速度計測するための手順。
選定根拠と段階構成は [`plans/active/tooling-cloud-gpu-env.md`](../plans/active/tooling-cloud-gpu-env.md) を参照。

## 構成サマリ

| 項目 | 決定値 |
| --- | --- |
| リージョン | `us-east-1` |
| 基準インスタンス | `g5.xlarge` (1×A10G 24GB, CC 8.6 = ローカル RTX 3060 と同一 CC) |
| MPI 用 (将来) | `g6.12xlarge` (4×L4) または `g4dn.12xlarge` (4×T4) |
| AMI | Deep Learning Base AMI with Single CUDA (Ubuntu 24.04) — 旧 "Base OSS Nvidia Driver GPU AMI (22.04)" の後継 |
| ストレージ | EBS gp3 200GB (stop 中も保持、~$0.08/GB·月) |
| 成果物回収 | S3 (`sync_results_s3.sh`) |
| コスト保険 | idle 30 分で自動 stop + AWS Budgets $50/月 |

## 事前準備 (一度だけ)

1. **vCPU quota**: `us-east-1` の「Running On-Demand G and VT instances」「All G and VT Spot Instance Requests」を各 48 vCPU へ引き上げ (承認済みであること)。
2. **S3 バケット**: バケット名はグローバル一意なのでサフィックスを付ける (例 `forge-runs-<name>`)。

   ```bash
   aws s3 mb s3://forge-runs-<name> --region us-east-1
   ```

3. **Budgets**: 月 $50 のコストアラートを作成 (Billing → Budgets → Cost budget)。
4. **キーペア / セキュリティグループ**: SSH (22) は自宅 IP 制限で 1 つ作成。inbound は SSH のみ。
5. **GitHub deploy key**: 計算機からの clone 用に read-only deploy key を `forger-yuwa/forge` に登録する。鍵はローカルで `ssh-keygen -t ed25519 -f forge-deploy` して公開鍵を GitHub → リポジトリ Settings → Deploy keys へ。秘密鍵は起動したインスタンスの `~/.ssh/` に `scp` で配置し、`~/.ssh/config` に:

   ```
   Host github.com
     IdentityFile ~/.ssh/forge-deploy
   ```

## インスタンス起動

コンソールなら: EC2 → インスタンスを起動 → AMI 検索で「**Deep Learning Base AMI with Single CUDA (Ubuntu 24.04)**」→ `g5.xlarge` → ストレージ 200GB gp3 → 上記キーペア/SG。

- 「Base」= ドライバ + Docker + nvidia-container-toolkit 入り・フレームワーク無し (これで十分)。
- 「Single CUDA」で足りる (計算は CUDA 同梱のコンテナ内で閉じるため。Multi CUDA は不要)。
- ホスト Ubuntu 24.04 の既定 gcc は 13。native ビルドが gcc 非対応でこけたら `sudo apt install gcc-12 g++-12` (`build_native_wsl.sh` は gcc-12 があれば自動で優先する)。

CLI なら (AMI ID は SSM パラメータで常に最新を引ける):

```bash
ami=$(aws ssm get-parameter --region us-east-1 \
  --name /aws/service/deeplearning/ami/x86_64/base-oss-nvidia-driver-gpu-ubuntu-24.04/latest/ami-id \
  --query Parameter.Value --output text)
aws ec2 run-instances --region us-east-1 \
  --image-id "$ami" --instance-type g5.xlarge \
  --key-name <keypair> --security-group-ids <sg-id> \
  --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":200,"VolumeType":"gp3"}}]' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=forge-g5}]'
```

## 初回セットアップ

SSH ログイン (ユーザーは `ubuntu`) し、deploy key を配置してから:

```bash
git clone --depth 1 git@github.com:forger-yuwa/forge.git ~/forge
~/forge/solver_density_cuda/tools/cloud/setup_instance.sh
```

スクリプトが行うこと: GPU/Docker 確認 → HighFive submodule 取得 → クラウド用イメージ `forge-solver:cuda-cloud` のビルド ([`Dockerfile.cuda`](../solver_density_cuda/Dockerfile.cuda) の `cloud` ステージ: dev から ParaView GUI/Qt を除き、gmsh はヘッドレスで同梱) → コンテナ内 Release ビルド → idle auto-stop cron の登録。

### 移設で必ず踏む 4 点 (2026-09-05, case/46 SERN チェーンで実測)

1. **submodule**: `git clone` 直後は HighFive が無く cmake が `HighFive headers not found` で止まる →
   `GIT_SSH_COMMAND="ssh -i ~/.ssh/forge-deploy" git submodule update --init --recursive`。
2. **CUDA 13 の thrust/cub**: CUDA 13.x は thrust/cub を `$CUDA_HOME/include/cccl/` に移したため
   `fatal error: thrust/extrema.h: No such file or directory` になる → configure に
   `-DCMAKE_CUDA_FLAGS=-I/usr/local/cuda/include/cccl -DCMAKE_CXX_FLAGS=-I/usr/local/cuda/include/cccl` を足す。
   `CUDAARCHS` も必須 (g5 の A10G は `86`)。
3. **converter の終了 assert**: `convertGmshToForge` は h5 を書き切ってから終了時に
   `GPUassert: invalid argument (cudaWrapper.cu 73)` で非零終了する (既知・無害)。
   成否は **exit code ではなく出力ファイルの存在とサイズ**で判定する
   (`evaluate/runner_sern.py::convert_mesh` がその形。`check=True` にすると全 run が落ちる)。
4. **design/ の python 依存**: `numpy scipy h5py pyyaml matplotlib` だけでは `opt/` が動かない。
   MOO には **`pymoo` / `scikit-learn` / `smt`** も要る (ローカルと同版を推奨):

   ```bash
   python3 -m venv design/.venv-opt
   design/.venv-opt/bin/pip install numpy scipy h5py pyyaml matplotlib \
       "pymoo==0.6.2" "scikit-learn==1.9.0" "smt==2.14.1"
   ```

`evaluate/` の runner はリポジトリ位置から forge のパスを導く (`FORGE_ROOT` 環境変数で上書き可)
ので、`~/forge` 以外に置いても動く。

## メッシングと可視化 (クラウドで完結させる)

### メッシュ生成 — コンテナ内 gmsh (ヘッドレス)

リポジトリのメッシュ生成はスクリプト駆動 (`case/*/mesh/make_*.py`, `.geo`) なので GUI 不要。
イメージに gmsh (バイナリ + python API) が入っており、ローカルと同じ手順がそのまま通る:

```bash
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp/forge-home \
  -v ~/forge:/repo -w /repo/case/<case>/mesh \
  forge-solver:cuda-cloud python3 make_XXXX.py
```

`.geo` を対話編集したいときだけ後述の NICE DCV を使う (通常は不要)。
メッシュ品質チェック (`check_mesh_quality.py`) と HDF5 変換もクラウド側で通常ルール通り実施する。

### 定型の後処理図 — matplotlib (同梱)

`residual_history.png`・line profile・断面図などのスクリプト後処理はイメージ内の
python3-h5py/numpy/matplotlib で完結する (`plot_residual.py` 等をコンテナ内で実行)。

### 対話的な 3D 可視化 — pvserver + SSH トンネル (既定)

`res_*.h5` (数 GB になり得る) をダウンロードせず、クラウド GPU 側でレンダリングして
ローカルの ParaView クライアントから操作する。

1. **バージョン一致が必須** (クライアントとサーバは同一版)。ローカルで使う ParaView の版を確認し、
   同じ版の Linux サーバ用バイナリ (`egl` 版 = ヘッドレス GPU レンダリング対応) を
   [paraview.org/download](https://www.paraview.org/download/) からインスタンスへ取得・展開する。
2. インスタンス側で起動:

   ```bash
   ~/paraview/bin/pvserver --server-port=11111
   ```

3. ローカル側でトンネルを張り、ParaView → File → Connect → `localhost:11111`:

   ```bash
   ssh -L 11111:localhost:11111 ubuntu@<instance>
   ```

セキュリティグループは追加開放不要 (SSH トンネル経由のため 22 のみで済む)。

### フル GUI が欲しいとき — NICE DCV (任意)

gmsh GUI や ParaView をデスクトップごと使いたい場合は、EC2 上でライセンス無料の
NICE DCV (GPU アクセラレートされたリモートデスクトップ) をインストールする。
手順は AWS 公式 (DCV Server + ubuntu-desktop) に従い、接続は DCV クライアントか
ブラウザ。ポートは SSH トンネル (`-L 8443:localhost:8443`) で通せば SG 開放不要。
帯域を食うので常用は pvserver、DCV は対話編集が必要なときだけ立てる。

## 計算実行

ローカルの Docker 手順と同じ。イメージ名だけ `forge-solver:cuda-cloud` に読み替える:

```bash
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v ~/forge:/repo -w /repo/case/<case>/run_XXXX_<slug> \
  forge-solver:cuda-cloud /repo/solver_density_cuda/build/forge solverConfig.yaml
```

`run_*` の複製・NaN/収束チェック・README の run 一覧表更新など通常ルール ([AGENTS.md](../AGENTS.md)) はクラウドでも同一に適用する。

## コード改修のイテレーション

2 モードを使い分ける:

- **モード A: push → pull (既定)** — ローカルで編集して feature ブランチへ push し、クラウドで pull → 差分ビルド。shallow clone のため初回のみ明示 fetch が要る:

  ```bash
  git fetch --depth 1 origin feature/<name>:feature/<name>
  git checkout feature/<name>
  docker run --rm --gpus all --user "$(id -u):$(id -g)" \
    -v ~/forge/solver_density_cuda:/workspace -w /workspace \
    forge-solver:cuda-cloud ./tools/build.sh
  ```

- **モード B: VS Code Remote-SSH で直編集** — カーネル最適化などの試行錯誤ループ用。インスタンスの `~/forge` を Remote-SSH で直接編集し、編集→ビルド→`ncu` を commit を挟まず回す。成果が出た変更のみ commit/push する。インスタンスから push するには deploy key を書き込み可で登録する (read-only のままなら確定 diff をローカルへ持ち帰って push)。

**stale build の罠 (既知)**: `solverConfig.hpp` などヘッダ変更後の差分ビルドは CUDA obj の
取りこぼしで step0 dt=0/NaN 凍結を起こした前歴がある。ヘッダを触ったら
`rm -rf build` してフルビルドすること。

## 速度計測・プロファイル

速度評価は Docker でなく **native ビルド**で行う (ローカルと同じルール。EC2 の Linux native では WSL と違い `nsys` の解析まで完結する):

```bash
# 依存 (初回のみ)。apt の python3-h5py は numpy と ABI 不整合で壊れるので pip で入れる
sudo apt-get install -y cmake ninja-build libyaml-cpp-dev libhdf5-dev libmetis-dev \
  libboost-dev python3-numpy python3-matplotlib
sudo python3 -m pip install --break-system-packages h5py

cd ~/forge/solver_density_cuda
# ホスト CUDA が 13.x の場合、g++ 単体コンパイルに Thrust (CCCL) の include が要る
CXXFLAGS="-I/usr/local/cuda/include/cccl" \
FORGE_CUDA_ARCHITECTURES=86 ./tools/build_native_wsl.sh   # A10G=86, L4=89, T4=75
```

- 粗い計時: `FORGE_PROFILE=1` (詳細 `FORGE_PROFILE_VERBOSE=1`) で実行し Runtime Profile Summary を読む。
- `ncu`: `ERR_NVGPUCTRPERM` が出たら `sudo` で実行するか、`/etc/modprobe.d/` に `options nvidia NVreg_RestrictProfilingToAdminUsers=0` を書いて再起動。
- `nsys`: Base DLAMI には入っていない (CUDA 同梱の `/usr/local/cuda/bin/nsys` はスタブでエラーになる)。スタンドアロン版を入れて直接叩く:

  ```bash
  sudo apt-get install -y nsight-systems-2026.1.3   # apt-cache search ^nsight-systems で最新を確認
  sudo chmod 755 /opt/nvidia/nsight-systems         # パッケージが 0700 で入る場合がある
  mkdir -p ~/tmp                                    # /tmp/nvidia が root 専有のことがある
  TMPDIR=$HOME/tmp /opt/nvidia/nsight-systems/2026.1.3/bin/nsys profile -o rep ./forge
  TMPDIR=$HOME/tmp /opt/nvidia/nsight-systems/2026.1.3/bin/nsys stats rep.nsys-rep --report cuda_gpu_kern_sum
  ```

  `nsys profile` → `nsys stats` までインスタンス内で完結する (WSL では不可能だった部分)。
- 計測はインスタンスタイプ (= GPU) をローカルの記録と併記し、同一 run 3 回でばらつきを確認してから比較に使う。

## 結果回収と停止

```bash
# 成果物を S3 へ (リポジトリ相対パスのままキー化される)
FORGE_S3_BUCKET=s3://forge-runs-<name> \
  ~/forge/solver_density_cuda/tools/cloud/sync_results_s3.sh case/<case>/run_XXXX_<slug>

# 停止 (EBS のみ課金になる)
sudo shutdown -h now        # または コンソール/CLI から stop-instances
```

- **stop** は EBS が残り再開できる (通常はこれ)。**terminate** はディスクごと消える (撤収時のみ)。
- 消し忘れ保険: `idle_autostop.sh` が root cron で 5 分おきに監視し、「GPU 使用率 0 + forge プロセス無し + ログイン無し」が 30 分続くと自動 shutdown する。長時間バッチは forge プロセスが生きている限り止まらない。
- ローカルへの取り込みは `aws s3 sync s3://forge-runs-<name>/case/... case/...` で逆方向に。

## spot で長時間バッチを回す場合

- 中断リスクを `outStepInterval` を適切に設定した `res_*.h5` からの restart で受ける (`restart_field.py` / 同一メッシュ restart)。
- ただし SST の restart は乱流場を完全再現しない既知の非忠実性がある (restart 直後の `vis_turb`・壁 omega が再現されず、擬似衝撃波位置など敏感な量が動く)。そうした量を扱う run は on-demand で通すのが安全。

## コスト早見 (2026-07, us-east-1, 参考値)

| 項目 | 単価 | 備考 |
| --- | --- | --- |
| g5.xlarge on-demand | ~$1.01/h | 基準。spot で ~$0.3-0.5/h |
| g6.12xlarge on-demand | ~$4.6/h | MPI 検証時のみ短時間 |
| EBS gp3 200GB | ~$16/月 | stop 中も課金される。長期離脱時は AMI 化 or terminate |
| S3 | ~$0.023/GB·月 | 成果物置き場 |
