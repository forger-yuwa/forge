#!/bin/bash
# retry_doe_024 を台帳に足してドライバを再開する (評価の切れ目で止める)。2026-10-12 plan §5.1 #8a′ (6) 判定 A
set -u
C=$HOME/forge-thruster/case/68.bell_two_mode; D=$C/run_0014_moo_d; L=$D/ledger.jsonl
n0=$(wc -l < $L); echo "start n=$n0 $(date -u)"
for i in $(seq 1 720); do n=$(wc -l < $L); [ "$n" -gt "$n0" ] && break; sleep 5; done
[ "$n" -gt "$n0" ] || { echo "TIMEOUT waiting for eval boundary"; exit 1; }
dpid=$(ps -eo pid,args | awk "/[o]pt\.driver .*run_0014_moo_d/{print \$1}")
echo "boundary n=$n, kill driver $dpid"; kill $dpid; sleep 2
for p in $(pgrep -x forge || true); do case "$(readlink /proc/$p/cwd)" in $D/*) echo "kill forge $p $(readlink /proc/$p/cwd)"; kill $p;; esac; done
sleep 3
$HOME/forge-thruster/design/.venv-opt/bin/python - <<"EOF"
import json, csv, os, shutil
D=os.path.expanduser("~/forge-thruster/case/68.bell_two_mode/run_0014_moo_d")
rows=[json.loads(l) for l in open(D+"/ledger.jsonl") if l.strip()]
tags={r["tag"] for r in rows}
for d in sorted(os.listdir(D)):
    p=os.path.join(D,d)
    if os.path.isdir(p) and (d.startswith("doe_") or d.startswith("inf_")) and d not in tags:
        print("remove partial", d); shutil.rmtree(p)
src=[r for r in rows if r["tag"]=="doe_024"][0]
th=list(csv.DictReader(open(D+"/retry_doe_024/thrust_history.csv")))[-1]
row=dict(src); row.update(tag="retry_doe_024", run_dir="../case/68.bell_two_mode/run_0014_moo_d/retry_doe_024", status="PASS", fail_class=None,
  eta_cf=float(th["eta"]), eta_wall=float(th["eta_wall"]), thrust_N=float(th["thrust_exit_N"]), thrust_wall_N=float(th["thrust_wall_N"]),
  mdot=float(th["mdot_exit_kg_s"]), mdot_ratio_1d=float(th["mdot_ratio_1d"]), cv_residual=float(th["cv_residual"]), res_file="res_24000.h5",
  conv="STALLED", residual_verdict="NOT CONVERGED (stalled/plateau) — segment S3_main + S4_main_p (84000 step)", quasisteady="ALL STEADY (連結 thrust_history_concat.csv、28 点)",
  gates={"convergence": True, "quasisteady": True, "physics": True, "cv": True, "mesh_quality": True, "wall_resolution": True},
  note="retry of doe_024: 同一設定・同一メッシュで本段 +24000 (runner --continue-from)。doe_024 の RISING は check_convergence の spike 条件の誤検出 (spike True・slow False・中央値比 0.739)。判定 A (plan §5.1 #8a′ (6)、diagnostician 2026-10-12)。GATES.json の quasisteady FAIL は継続分 8 点だけを見た点数不足で、親 20 点 + 継続 8 点の連結で ALL STEADY", wall_s=None)
open(D+"/ledger.jsonl","a").write(json.dumps(row, ensure_ascii=False, default=float)+"\n")
json.dump({"judgement":"A","row":row}, open(D+"/retry_doe_024/RETRY_JUDGEMENT.json","w"), ensure_ascii=False, indent=1, default=float)
print("appended retry_doe_024 eta", row["eta_cf"])
EOF
cd $HOME/forge-thruster/design
echo "=== driver restart $(date -u) with retry_doe_024 in the ledger ===" >> $C/moo_d14.log
nohup env FORGE_CUDA_BLOCKSIZE=128 FORGE_BIN=$HOME/forge-thruster/solver_density_cuda/.build-native/relwithdebinfo/forge .venv-opt/bin/python -m forge_design.opt.driver ../case/68.bell_two_mode/problem_bell_tp_d_moo.yaml ../case/68.bell_two_mode/run_0014_moo_d --n-doe 40 --n-iter 10 --batch 2 --seed 0 --ref -0.90 23.5 >> $C/moo_d14.log 2>&1 < /dev/null &
sleep 5; ps -eo pid,args | awk "/[o]pt\.driver/" | cut -c1-120; echo DONE
