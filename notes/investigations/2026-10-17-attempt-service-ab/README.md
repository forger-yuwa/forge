# 試行ごとの service の契約: 実行担当への SIGKILL の A/B と、ssh を切ったときの生存 (2026-10-17)

plan [`tooling-design-orchestrator`](../../../plans/active/tooling-design-orchestrator.md) §4.5「停止の責任は実行担当の外に置く」・§5.1 #15 の実測。
GPU・forge は使っていない (ダミーのプロセスだけ)。

## 何を確かめたか

設計の実行役は、試行ごとに systemd の user service を 1 つ作り、実行担当とその子孫 (`run_case.sh`・forge) を収める。
停止規則 v2 は `run_case.sh` を別のセッション (`start_new_session=True`、`design/forge_design/evaluate/stop_rule.py:1901`) で起動し、
例外のときの forge の停止も同じ監視のプロセスが行う。監視のプロセスが `SIGKILL` で死ぬと、その停止の処理は動かない。
service の `KillMode=control-group` なら、実行担当が死んだ時点で systemd が cgroup ごと止める、という契約を確かめた。

- 試験の台本: [`design/tests/run_attempt_service_ab.py`](../../../design/tests/run_attempt_service_ab.py)
- service の設定: `ExitType=main`・`KillMode=control-group`・`TimeoutStopSec=5`・`Restart=no` (`systemd-run --user --collect`)
- 判定 (plan §4.5 に事前に書いた): A (介入なし 10 秒) で 3 つ (実行担当・ダミーの run_case・ダミーの forge) が全観測で生き、
  B (3 つがそろって 2 秒の時点で実行担当だけを `SIGKILL`) で全子孫が停止期限 5 秒の内に消えれば成立。
- 対照 (記録だけ): service に入れず同じ形で起動して実行担当を `SIGKILL` する。

## 結果

| 環境 | systemd | A | B (SIGKILL → 全子孫が消えるまで) | 対照 (service なし) | 記録 |
| --- | --- | --- | --- | --- | --- |
| 手元 (WSL2) | 255 | PASS (101 回の観測で 3 つとも生存) | PASS (0.101 秒) | run_case・forge が残った | [`local_wsl.json`](local_wsl.json) |
| AWS B (g5.xlarge) | 255 | PASS (101 回) | PASS (0.100 秒) | run_case・forge が残った | [`aws_b.json`](aws_b.json) |

- 両環境とも、3 つのプロセスはその試験の unit の cgroup に入り (`/proc/PID/cgroup` が `…/forge-attempt-ab-<arm>-<id>.service`)、
  実行担当と run_case はセッションが分かれていた (`sessions_differ: true`)。停止規則 v2 と同じ形でも、cgroup で止められる。
- 対照では、実行担当を `SIGKILL` した後もダミーの run_case と forge が生き残った (codex の観測 2026-10-17 を再現)。この試験は「残る」を検出できる。
- 後始末の後、両環境とも試験のプロセスは残っていない。

## ssh を切ったときの生存 (Linger)

AWS B の既定は `Linger=no`。この状態で user service を起動して ssh を切ると、約 12 秒後に systemd が service を止めた
(journal: 10:23:59 Started → 10:24:11 Stopping/Stopped)。`sudo loginctl enable-linger ubuntu` で `Linger=yes` にすると、
ssh を切って 25 秒後も生きていた (10:24:41 Started、25 秒後に ALIVE)。確かめた後は `disable-linger` で `Linger=no` に戻した。

- AWS で実行役を回すには、そのインスタンスの `Linger=yes` が要る (1 回だけの設定)。実行役の事前検査で `Linger=no` なら起動を止める。
- 手元 (WSL2) も `Linger=no`。手元で長い計算を回すときも同じ。

## 制約

- 試験したのは停止の仕組みだけで、実行担当の実装・予算の監視・台帳はまだ無い (plan §5.1 #6)。
- `TimeoutStopSec` で止まらない (D 状態など) プロセスの扱いは試していない。plan §4.5 の「停止期限の内に終了を確かめられなければ新しい試行を起動しない・予算を解放しない」で受ける。
