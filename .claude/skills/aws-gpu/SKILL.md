---
name: aws-gpu
description: AWS の GPU インスタンス (g5.xlarge) を自分で起動・停止して run を回す。aws_instance.sh の使い方、認証情報の置き場所 (/mnt/c に置くと chmod が効かない)、idle 自動停止が正常であること、落ちた run の見分け方と再投入、ローカルと AWS の使い分け。「AWS で回して」「インスタンス起動して」と言われたとき、またはローカルの資源では足りない run を投入するときに使う。
---

# /aws-gpu — GPU インスタンスを自分で起こして回す

手順の正本は **[`procedures/cloud-aws-gpu.md`](../../../procedures/cloud-aws-gpu.md)**。
本 skill はそこに書かれていない **運用の癖と失敗パターン**だけを足す。重複は書かない。

## 1. 起こす・止める

```bash
export AWS_SHARED_CREDENTIALS_FILE=/home/sano/.aws-wsl/credentials \
       AWS_CONFIG_FILE=/home/sano/.aws-wsl/config AWS_PROFILE=forge
solver_density_cuda/tools/aws_instance.sh status   # state / IP / type
IP=$(solver_density_cuda/tools/aws_instance.sh start)   # running + SSH 可まで待つ (実測 26 秒)
solver_density_cuda/tools/aws_instance.sh stop
```

- **新しいシェルでは `.bashrc` が効くので export は不要**。効かない環境 (`NoCredentials` / `NoRegion`) では
  上の 3 変数を明示する。`NoRegion` は `config` を `[forge]` と書いたときに出る — **`config` だけ
  `[profile forge]`** が正しい (`credentials` は `[forge]`)。
- **IP は毎回変わる**。変数に取って使い回す。ハードコードした IP を次のターンで再利用しない。

## 2. 止まっているのは異常ではない

インスタンス側の `idle_autostop.sh` (root cron) が **「GPU 0 % + forge プロセス無し + ログイン無し」が
30 分**続くと shutdown する。2026-09-24 に 3 回落ちたが、いずれも run 終了後にローカルで解析・plan 更新を
していた 30 分だった。**これはコスト保護が正しく働いた結果**である。

- **回避しないこと**。偽のプロセスを走らせて保護を無効化しない。起動が 26 秒で済む以上、起こす方が安い。
- 長時間バッチは forge が生きている限り止まらない。**だから run は `nohup` でインスタンス側に残す**
  (SSH セッションが切れても走り続ける)。

## 3. 落ちた run を必ず見分ける

**停止をまたいだ run は「走った」と仮定しない。** 実際に `run_0448` を投入したつもりで dump が 0 個だった
(起動前に停止していた)。再接続したら**必ず成果物を数える**:

```bash
ls RUN/res_[0-9]*.h5 | sed 's/.*res_//;s/\.h5//' | sort -n | tail -1   # **sort -n**。文字列ソートは 9000 > 27000 になる
pgrep -f relwithdebinfo/forge >/dev/null && echo RUNNING || echo DONE
```

`ls | tail` の文字列ソートで**古い dump を最終値と誤認した事例がある**。必ず `sort -n`。

## 4. SSH の待ち方

`nohup ... &` を含む ssh は**手元の timeout で切られても向こうは走り続ける**。切られたら失敗と決めつけず、
`pgrep` で確認する。完了待ちは until ループを `run_in_background` で回す:

```bash
until ! ssh ... 'pgrep -f relwithdebinfo/forge >/dev/null'; do sleep 45; done
```

## 5. ローカルと AWS の使い分け

- **ローカル (RTX 3060 12 GB / RAM 11 GB / ディスクに余裕なし)**: 数万節点・1 step 級の切り分け、
  面流束ダンプの照合、メッシュ生成と品質検査。
- **AWS (A10G 23 GB)**: 100 万節点以上、数万 step、格子感度の対。
- **速度評価・プロファイルは native を既定** (AGENTS.md)。Docker で測らない。
- 転送は `rsync -az`。**同時に 2 本走らせると帯域が半分になる** (実際に踏んだ)。

## 6. 終わったら止める

自分の run が終わり、結果を取り込んだら `stop` する。idle 保護はあるが 30 分ぶん課金される。
**他セッションの run が走っていないことを確認してから止める** (`pgrep -f relwithdebinfo/forge`)。
