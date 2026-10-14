#!/bin/bash
# plan architecture-float-state-double-geometry §6.61: 速さの比較 S を B で逐次に回す (G0c → G0)、終わったら判定。
set -uo pipefail
C45=~/forge-wallfit/case/45.isobutane_m6_d155; cd $C45
bash sp_arm.sh run_0804_sp_g0c g0c
bash sp_arm.sh run_0805_sp_g0 g0
python3 sp_judge.py >> sp.log 2>&1
touch sp.done
