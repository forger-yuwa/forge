#!/usr/bin/env bash
# forge run キュー 1 (順次): A 断熱 / B 300 K / C 700 K, 全て y1=3um 生産 SST
cd "$(dirname "$0")"
python3 gen_runs.py --run run_0001_A_ad_y3   --mesh fp_y1_3um --wall adiabatic > _queue1_A.log 2>&1
python3 gen_runs.py --run run_0002_B_tw300_y3 --mesh fp_y1_3um --wall 300      > _queue1_B.log 2>&1
python3 gen_runs.py --run run_0003_C_tw700_y3 --mesh fp_y1_3um --wall 700      > _queue1_C.log 2>&1
echo QUEUE1_DONE >> _queue1_C.log
