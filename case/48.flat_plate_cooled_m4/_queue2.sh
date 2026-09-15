#!/usr/bin/env bash
# キュー 2 (順次): C 700K (A から) / B の y1+ 掃引 (B から cross-mesh) / A-plain, B-plain (同一メッシュ, 素 SST)
cd "$(dirname "$0")"
while ! grep -q "main rc" _B.log 2>/dev/null; do sleep 10; done
python3 gen_runs.py --run run_0006_C_tw700_y3 --mesh fp_y1_3um --wall 700 --stages soft --ramp "0.5,1,2" --cfl 2 --main-steps 48000 --out-int 4000 --ic-from run_0004_A_ad_y3_cont > _C.log 2>&1
for y in 6 12 24; do
  python3 gen_runs.py --run run_$(printf %04d $((y==6?7:(y==12?8:9))))_B_tw300_y${y} --mesh fp_y1_${y}um --wall 300 --stages soft --ramp "0.5,1,2" --cfl 2 --main-steps 36000 --out-int 4000 --ic-from run_0005_B_tw300_y3 > _B_y${y}.log 2>&1
done
python3 gen_runs.py --run run_0010_Aplain_ad_y3 --mesh fp_y1_3um --wall adiabatic --plain --stages soft --ramp "1,2" --cfl 2 --main-steps 36000 --out-int 4000 --ic-from run_0004_A_ad_y3_cont > _Aplain.log 2>&1
python3 gen_runs.py --run run_0011_Bplain_tw300_y3 --mesh fp_y1_3um --wall 300 --plain --stages soft --ramp "1,2" --cfl 2 --main-steps 36000 --out-int 4000 --ic-from run_0005_B_tw300_y3 > _Bplain.log 2>&1
echo QUEUE2_DONE >> _Bplain.log
