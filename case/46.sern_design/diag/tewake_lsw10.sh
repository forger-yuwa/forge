#!/bin/bash
# plan tooling-sern-te-wake-grid §5.1 #6c: L_sw 1.0 の 3 作動点を修正後の格子 (局所化) で段階起動 → 20000 step の判定区間 (AWS-B)。
set -u
cd /home/ubuntu/forge-r8
export FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_CUDA_BLOCKSIZE=128
C=case/46.sern_design; PA=problem_3d_prod_3op_wallres_lswx10_tewake_g3_A0.yaml; PL=problem_3d_prod_3op_wallres_lswx10_tewake_g3_L10.yaml
A=run_1092_tewake_lsw10_A0_ref
declare -A RUN=([m10_on]=run_1093_tewake_lsw10_m10 [m6_on]=run_1094_tewake_lsw10_m6 [m4_off]=run_1095_tewake_lsw10_m4)
declare -A CONT=([m10_on]=run_1096_tewake_lsw10_m10_cont20k [m6_on]=run_1097_tewake_lsw10_m6_cont20k [m4_off]=run_1098_tewake_lsw10_m4_cont20k)
prep() { prob=$1; r=$2; op=$3
  [ -e $C/$r ] && { echo "EXISTS $r"; exit 1; }
  (cd design && python3 -m forge_design.evaluate.runner_sern3d ../$C/$prob ../$C/$r --op $op --prepare-only > ../$C/$r.prepare.log 2>&1 < /dev/null); echo "$r prepare rc=$?"
  python3 solver_density_cuda/tools/check_dual_closure.py $C/$r/sern.h5 --tol 1e-5 > $C/$r/DUAL_CLOSURE.txt 2>&1; echo "$r dual rc=$?"
  rm -f $C/$r/sern.msh; }
prep $PA $A m10_on
ok=1
for op in m10_on m6_on m4_off; do r=${RUN[$op]}
  prep $PL $r $op
  python3 $C/diag/te_wake_grid_check.py $C/$A $C/$r --admission --dual-closure-b $C/$r/DUAL_CLOSURE.txt --dual-closure-a $C/$A/DUAL_CLOSURE.txt --code-commit-a "$(git rev-parse HEAD)" --code-commit-b "$(git rev-parse HEAD)" --provenance-out $C/$r/TE_WAKE_PROVENANCE.json > $C/$r/TE_WAKE_ADMISSION.txt 2>&1
  python3 $C/diag/te_wake_hex_compare.py $C/$A $C/$r --ar-max 5000 --json $C/$r/HEX_COMPARE.json > $C/$r/HEX_COMPARE.txt 2>&1
  python3 $C/diag/twin_ic_check.py --g3 $C/run_1084_tewake_g3_L10_m10 --src-res $C/run_1084_tewake_g3_L10_m10/res_20000.h5 --g4 $C/$r --json $C/$r/TWIN_IC.json > $C/$r/TWIN_IC_VERDICT.txt 2>&1
  echo "== $op $r"; grep "^行ごと" $C/$r/TE_WAKE_ADMISSION.txt; grep "TE_WAKE_LOCAL VERDICT" $C/$r/HEX_COMPARE.txt; tail -1 $C/$r/TWIN_IC_VERDICT.txt
  grep -q "TE_WAKE_LOCAL VERDICT: PASS" $C/$r/HEX_COMPARE.txt || ok=0
  grep -q "TWIN_IC VERDICT: B_OK" $C/$r/TWIN_IC_VERDICT.txt || ok=0
  # 投入条件の行: 全域 0.65 の行 (primal の品質) 以外が PASS であること (ユーザ判断 2026-10-11)
  grep -E "^(FAIL|判定不能)" $C/$r/TE_WAKE_ADMISSION.txt | grep -v "スケール済みヤコビアンの最小" | grep -q . && ok=0
done
[ $ok = 1 ] || { echo "GRID/IC CHECK NOT OK -> stop before time integration"; echo LSW10_DONE_STOP $(date); exit 0; }
staged() { op=$1; r=${RUN[$op]}
  (cd design && python3 - <<PY > ../$C/$r.staged.log 2>&1
from forge_design.evaluate.runner_sern3d import run_staged, collect
from forge_design.evaluate.runner_sern import load_problem
o = (load_problem("../$C/$PL").raw.get("opt") or {})
rc = run_staged("../$C/$r", "full", int(o.get("soft_steps", 2000)), soft_cfl=float(o.get("soft_cfl", 0.5)), soft_conv=int(o.get("soft_conv", 0)),
                warm_lam_steps=int(o.get("warm_lam_steps", 0)), warm_lam_cfl=float(o.get("warm_lam_cfl", 0.2)), mid_steps=int(o.get("mid_steps", 0)))
out = collect("../$C/$PL", "../$C/$r", rc=rc)
g = out["gates"]; print("STAGED rc", rc, "GATES", g["verdict"], g["fail_class"], g["reasons"])
PY
  )
  echo "$r $(tail -1 $C/$r.staged.log)"
  [ -n "$(ls $C/$r | grep -E '^res_[0-9]+\.h5$' | head -1)" ] && $C/diag/tewake_ext.sh $r ${CONT[$op]} $PL > $C/${CONT[$op]}.ext.log 2>&1
  echo "${CONT[$op]} $(tail -1 $C/${CONT[$op]}.ext.log)"; }
staged m10_on & sleep 120; staged m6_on & sleep 120; staged m4_off & wait
echo LSW10_DONE $(date)
