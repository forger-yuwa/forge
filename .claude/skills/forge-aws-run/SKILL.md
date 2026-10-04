---
name: forge-aws-run
description: 共有 AWS GPU インスタンスで forge の run を投入・継続・待ち合わせ・後片付けするときの手順と、踏んだ罠の一覧。起動確認・ssh 越しの完了判定 (pgrep の自己一致)・インスタンス停止の検知・ディスク逼迫・restart の型と化学種・実行中スクリプトの編集禁止を扱う。「AWS で回して」「延長して」「終わったら知らせて」「ディスク空けて」と言われたとき、または AWS 上で run を起動・監視するときに使う。
---

# /forge-aws-run — 共有 AWS で run を回す

前提: ローカル PC はユーザ作業で使うので、検証 run は小さくても AWS で回す。インスタンスは 1 台を**他セッションと共有**している。
接続・起動は `solver_density_cuda/tools/aws_instance.sh {status|ip|start|stop|ssh}` (別ブランチにしか無いときは
`git show <commit>:solver_density_cuda/tools/aws_instance.sh` で取り出す)。鍵は `~/.ssh/test.pem`、認証情報は `~/.aws-wsl/`。

## 1. 起動・停止

- **start は毎回ユーザに確認する** (課金が始まる)。stop は手動でしない — idle 自動停止に任せる (他セッションの準備中の作業を切らない)。
- 起動は `InsufficientInstanceCapacity` で失敗することがある。5 分おきの再試行を `run_in_background` で回し、成功で抜ける。
- 起動のたびに **public IP が変わる**。`aws_instance.sh status` の 2 列目を毎回読み直す。
- 触る前に `nvidia-smi` と `pgrep -x forge` の cwd で他セッションの利用を見る。GPU は共有してよいが、他人の run・ディレクトリには触らない。

## 2. 投入

- `nohup` で起動し、**stdin を `/dev/null` にする** (`(... &)` + `< /dev/null`)。そうしないと ssh がプロセス終了まで戻らない。
- `FORGE_CUDA_BLOCKSIZE=128` を必ず付ける (既定 512 は node SLAU のレジスタ上限を超えて起動できない。値は結果に効く)。
- 起動後に `RUN_PROVENANCE.txt` の `forge_bin` を確認する (倍精度ビルドなど別バイナリを使うとき)。
- 別バイナリ (例: 全域 FP64) は同じ commit の worktree を作り、`flowFormat.hpp` の typedef だけを変えて native build する。
  サブモジュール (HighFive 等) は元の checkout のものをシンボリックリンクで流用する。
- 既存のビルド済みツリーを複製して別バイナリを作るときは `build/` を持ち込まず新規に cmake する (CMakeCache が元の絶対パスを指す)。
  AWS では `-DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_CXX_FLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl"` が要る
  (CXX_FLAGS が無いと `speciesDB.cpp` 等で `vector_types.h: No such file` になる)。元の CMakeCache と CMAKE_CXX_FLAGS を突き合わせて確認する。

## 3. 継続・引き継ぎ (restart)

- **同一メッシュは `restart_field.py`、別メッシュは `interp_field.py`**。同一メッシュに interp_field を使わない (保存量が丸めでずれる)。
- **倍精度の場を種にするときは `restart_field.py --keep-src-dtype`**。付けないと「型の縮小で丸めあり」になる。
- **VERDICT の移した量の数を必ず見る**。新しく変換したメッシュには化学種 `roY*` が無いので、多成分 run では 12 量のはずが 7 量しか移らない
  → 事前に src の `roY*` データセットを dst に作っておく (値は restart_field が上書きする)。
- 継続先には `species_db.yaml` も複製する (忘れると forge が `failed to read speciesDBFile` で起動しない)。
- 継続は新しい run ディレクトリで行い、既存 run を上書きしない。収束判定は同一設定なら接続区間で行う (作用素が変わったら接続しない)。
- **保存量がビット一致でも再開直後に残差が跳ねる** (2026-09-27 case/60: rms_ro 6.8e-12 → 最大 5.5e-9 → 数千 step で戻る)。保存量以外の状態は移らない。
  延長 run の窓間変化・低下桁数にはこの戻りが入るので、A/B の両腕を同じ回数だけ再開して揃え、延長区間だけの判定では跳ねから測っていることを明記する。
- 上に層を積んだだけのメッシュ (節点座標が完全一致) へ場を移すときは interp_field でなく、座標一致でビット単位コピーする (case/60 `tools/stack_init.py`)。

## 4. 待ち合わせ (完了通知を取りこぼさない)

- **完了は成果物で判定する**: 残差 CSV の最終 step が `nStepOuter−1`、ビルドなら `build/forge` の mtime、ログの終了行。
- **ssh 越しの `pgrep -f "パターン"` は使わない** — `bash -c '<コマンド全文>'` 自体がパターンを含むので自分に一致する。
  待ち合わせでは「終わっても実行中」と出続け、`kill $(pgrep -f ...)` では自分のシェルを殺して exit 255 になる (2026-09-26/27 に 2 回)。
  プロセスを見る・止めるときは `pgrep -x forge` + `readlink /proc/PID/cwd` で run ディレクトリを絞る。
- ローカルで `&` 起動したプロセスを待つときは、外側シェルでなく**実プロセスの PID** を控えて `kill -0` で待つ。
- **ssh の失敗 (接続不可・インスタンス停止) はそれ自体を事象として即座に抜けて知らせる**。空の結果で回り続けると、
  idle 自動停止に気づかず 2 時間待つ (2026-09-27)。
- ssh の終了コードで接続失敗を判定するときは、リモート側のコマンド列の末尾に `; true` を付ける (`grep -c` は 0 件で終了コード 1 を返し、接続失敗と区別できなくなる)。ssh 自体の失敗は 255。
- **ssh の終了コードで接続失敗を判定するときは、リモート側のコマンドの終了コードと混ぜない**。`grep -c` は 0 件で終了コード 1 を返すので、`ssh ... "grep -c ERROR log"` を `|| { echo SSH FAIL; break; }` で受けると、エラー 0 件 = 正常なのに「接続失敗」でループを抜ける (2026-09-27: 実行中の run を完了扱いで判定した)。リモート側は `if grep -q ...; then echo ERR; else echo WAIT; fi` のように**常に 0 で終わる**形にして状態を文字列で返す。
- ループには上限と「タイムアウトした」出力を付ける。1 本終わるごとに抜ける形にして、終わった run から順に判定する。
- 実行中のシェルスクリプト (バッチ) を編集しない — bash は逐次読みなので再開位置がずれる。変えたいなら次の投入で。

## 5. ディスク

- 共有ディスクは他セッションが使うので急に減る。投入前と長時間 run の途中で `df -h ~` を見る。
- **満杯になると forge は `HighFive::DataSetException ... Write failed` で落ち、その後 idle 自動停止する** (2026-09-27: 自分の run の段出力 `_<段名>_res_*` と
  中間スナップショットで case/59 が 11 GB になり、A/B 2 本が step 1900 で死んだ)。待ち合わせは空き容量と `Write failed` も監視する。
- 段階起動の生成器が残す段出力 (`_<段名>_res_*`) と、継続 run の `res_0` (引き継ぎ元の最終場の複製) も判定後に消してよい。
- 自分の run は**判定が済んだら中間の全場スナップショット `res_<n>.h5` を消す**。残すもの: `res_0`・最終場・壁/表面出力
  (`res_wall_*`/`res_plate_*`/`res_gap_*`、準定常の時系列に使う)・残差履歴・VERDICT・時系列 CSV。消したことは case README に書く。
- 他セッションのディレクトリの削除はユーザに判断を仰ぐ (大きい順に候補を挙げるだけ)。

## 6. 判定

- 判定ツール (`check_convergence.py`・`check_quasisteady.py`・壁解像) は AWS 上で回し、結果と run パスを case README の run 一覧に書く。
- 事前登録した比較 (acceptance.json) は、run を回す前に commit しておく。
