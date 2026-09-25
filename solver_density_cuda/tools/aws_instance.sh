#!/usr/bin/env bash
# forge の GPU インスタンス (g5.xlarge) を操作する。**鍵はこのファイルに含まない**。
#
#   tools/aws_instance.sh status|ip|busy|start|stop|ssh
#
# **インスタンスは 1 台を複数セッションで共有している**。別セッションが起動して run を回して
# いることがある。`stop` は forge プロセス・他ログイン・GPU 使用があれば拒否する
# (上書きは FORGE_AWS_FORCE_STOP=1。ユーザが明示したときだけ使う)。`start` は既に running なら
# 使用状況を表示する — 自分が起こしたと思い込まず、他の run と GPU を取り合わないか確認すること。
#
# 認証情報は WSL ネイティブ側に置くこと (`/mnt/c` 配下は drvfs でパーミッションが効かず
# chmod 600 が 777 になる):
#   AWS_SHARED_CREDENTIALS_FILE=/home/sano/.aws-wsl/credentials   ([forge])
#   AWS_CONFIG_FILE=/home/sano/.aws-wsl/config                    ([profile forge])
# IAM は当該インスタンスの Start/Stop + DescribeInstances だけ (Terminate/RunInstances は付けない)。
#
# **止まっているのは正常**: インスタンス側の `idle_autostop.sh` が root cron で
# 「GPU 0 % + forge プロセス無し + ログイン無し」30 分で shutdown する
# (procedures/cloud-aws-gpu.md)。長時間バッチは forge が生きている限り止まらない。
set -euo pipefail
IID=${FORGE_AWS_INSTANCE:-i-0b1a5e0b8dc152f00}
REGION=${AWS_REGION:-us-east-1}
KEY=${FORGE_AWS_KEY:-$HOME/.ssh/test.pem}
# 使用状況: forge プロセス数 / 対話ログイン数 (コマンド指定の ssh は who に出ないので自分は数えない) / GPU 使用率 (%) を 1 行で返す
busy() {
  ssh -o StrictHostKeyChecking=no -o ConnectTimeout=10 -i "$KEY" "ubuntu@$(q PublicIpAddress)" \
    'printf "forge=%s logins=%s gpu=%s\n" "$(pgrep -xc forge || true)" \
       "$(who | wc -l)" \
       "$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | head -1)"'
}
q() { aws ec2 describe-instances --region "$REGION" --instance-ids "$IID" \
        --query "Reservations[0].Instances[0].$1" --output text; }

case "${1:-status}" in
  status) printf '%s\t%s\t%s\n' "$(q State.Name)" "$(q PublicIpAddress)" "$(q InstanceType)" ;;
  ip)     q PublicIpAddress ;;
  start)
    # **stopping / pending からは start できない** (IncorrectInstanceState)。終端状態まで待つ。
    for _ in $(seq 60); do
      st=$(q State.Name)
      case "$st" in running|stopped) break ;; esac
      sleep 5
    done
    if [ "$st" = running ]; then
      echo "既に running (別セッションが使用中の可能性): $(busy 2>/dev/null || echo '使用状況取得失敗')" >&2
    fi
    [ "$st" = running ] || aws ec2 start-instances --region "$REGION" --instance-ids "$IID" >/dev/null
    for _ in $(seq 60); do
      [ "$(q State.Name)" = running ] && break; sleep 5
    done
    ip=$(q PublicIpAddress)
    # SSH が上がるまで待つ (起動直後は接続を拒否する)
    for _ in $(seq 60); do
      ssh -o StrictHostKeyChecking=no -o ConnectTimeout=5 -i "$KEY" "ubuntu@$ip" true 2>/dev/null && break
      sleep 5
    done
    echo "$ip" ;;
  busy)   [ "$(q State.Name)" = running ] && busy || echo "not running" ;;
  stop)
    if [ "$(q State.Name)" = running ] && [ "${FORGE_AWS_FORCE_STOP:-0}" != 1 ]; then
      u=$(busy) || { echo "使用状況を取得できないので止めない (FORGE_AWS_FORCE_STOP=1 で強制)" >&2; exit 3; }
      case "$u" in
        "forge=0 logins=0 gpu=0") ;;
        *) echo "使用中なので止めない: $u (別セッションの run の可能性。FORGE_AWS_FORCE_STOP=1 で強制)" >&2; exit 3 ;;
      esac
    fi
    aws ec2 stop-instances --region "$REGION" --instance-ids "$IID" \
            --query 'StoppingInstances[0].CurrentState.Name' --output text ;;
  ssh)    exec ssh -o StrictHostKeyChecking=no -i "$KEY" "ubuntu@$(q PublicIpAddress)" ;;
  *) echo "usage: $0 status|ip|busy|start|stop|ssh" >&2; exit 2 ;;
esac
