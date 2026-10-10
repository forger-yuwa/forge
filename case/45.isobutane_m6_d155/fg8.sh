#!/bin/bash
# plan architecture-float-state-double-geometry §6.14 (事前登録) の確認 (2026-10-10): 変換器は常に double で幾何を書く。
#   1: 新 FP64 対 旧 FP64 (段 ④ + 診断 = ~/forge-fgeom5-fp64、56b1a6fe) の変換の出力が全件ビット一致 (case/45・48・39)
#   2: 新 float 対 新 FP64 が --geometry-only --require-float64 で一致
#   3: ソルバの警告 (float32 の格子で出る、新しい double の格子で出ない、混在の格子で出る)
#   4: 新 float のソルバが新しい double の格子を読み double の写しを持つ (段 ② のダンプ、FP64 のダンプとの比較)
#   5: case/56 (float32 の格子) の 20 step の互換性 (新 = fgeom6-f32、旧 = fgeom5-f32 と再実行)
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C56=~/forge-wallfit/case/56.gap_tp1187; C48=~/forge-wallfit/case/48.flat_plate_cooled_m4
W=$C45/_conv614; LOG=$C45/fg8.log; RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
CMP=~/forge-fgeom6-f32/solver_density_cuda/tools/compare_mesh_h5.py
V_OLD64=~/forge-fgeom5-fp64/solver_density_cuda/build/convertGmshToForge
V_NEW64=~/forge-fgeom6-fp64/solver_density_cuda/build/convertGmshToForge
V_NEW32=~/forge-fgeom6-f32/solver_density_cuda/build/convertGmshToForge
F_NEW32=~/forge-fgeom6-f32/solver_density_cuda/build/forge; F_NEW64=~/forge-fgeom6-fp64/solver_density_cuda/build/forge
F_OLD32=~/forge-fgeom5-f32/solver_density_cuda/build/forge
echo "== 開始 $(date -Is)" >> $LOG
for b in V_OLD64 V_NEW64 V_NEW32 F_NEW32 F_NEW64 F_OLD32; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
conv() {   # conv <case> <arm> <変換器> <msh> <出力名>
  local d=$W/$1/$2; mkdir -p $d; cp $W/inp/$1/*.yaml $d/
  ( cd $d && $3 $W/inp/$1/$4 $5 > conv.log 2>&1 ); local rc=$?
  if [ $rc -ne 0 ]; then
    if grep -q "Write Input HDF5" $d/conv.log && grep -q "GPUassert: invalid argument" $d/conv.log; then rc=0; fi   # 終了時の GPUassert は既知 (conv_tolerant.sh と同じ扱い)
  fi
  echo "conv $1 $2 rc=$rc" >> $LOG
}
for c in c45:nozzle.msh:nozzle.h5 c48:fp_y1_3um.msh:m.h5 c39:hill_des_80x50x30.msh:hill_des_80x50x30.h5; do
  IFS=: read cs msh out <<< "$c"
  conv $cs old64 $V_OLD64 $msh $out; conv $cs new64 $V_NEW64 $msh $out; conv $cs new32 $V_NEW32 $msh $out
  echo "== 1 $cs 新 FP64 対 旧 FP64 (全件)" >> $LOG; python3 $CMP $W/$cs/old64/$out $W/$cs/new64/$out >> $LOG 2>&1; echo "  rc=$?" >> $LOG
  echo "== 2 $cs 新 float 対 新 FP64 (--geometry-only --require-float64)" >> $LOG; python3 $CMP $W/$cs/new32/$out $W/$cs/new64/$out --geometry-only --require-float64 >> $LOG 2>&1; echo "  rc=$?" >> $LOG
done
# ---- 3・4: ソルバ (case/48 の格子) ----
mk48() {   # mk48 <run> <格子の h5>
  mkdir $C48/$1 && cp $C48/run_0046_fg2_dump32/{solverConfig.yaml,bcondConfig.yaml,probe.yaml} $C48/$1/ && cp $2 $C48/$1/mesh.h5
}
mk48 run_0048_fg6_dump32 $W/c48/new32/m.h5; mk48 run_0049_fg6_dump64 $W/c48/new64/m.h5
mk48 run_0050_fg6_warn_mixed $W/c48/new32/m.h5
python3 - "$C48/run_0050_fg6_warn_mixed/mesh.h5" <<'PY' >> $LOG 2>&1
import h5py, numpy as np, sys
with h5py.File(sys.argv[1], "r+") as h:
    a = np.asarray(h["/PLANES/surfVect"][:]); del h["/PLANES/surfVect"]; h.create_dataset("/PLANES/surfVect", data=a.astype(np.float32))
print("[fg8] 混在の格子: /PLANES/surfVect だけ float32 にした")
PY
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; }
run $C48/run_0048_fg6_dump32 $F_NEW32 FORGE_DIAG_GEOM_STAGE2_DUMP=$C48/run_0048_fg6_dump32/stage2.h5
run $C48/run_0049_fg6_dump64 $F_NEW64 FORGE_DIAG_GEOM_STAGE2_DUMP=$C48/run_0049_fg6_dump64/stage2.h5
run $C48/run_0050_fg6_warn_mixed $F_NEW32 FORGE_DIAG_GEOM_STAGE2_DUMP=$C48/run_0050_fg6_warn_mixed/stage2.h5
for r in run_0048_fg6_dump32 run_0049_fg6_dump64 run_0050_fg6_warn_mixed; do echo "== 3 警告 $r: $(grep -h 'geometry datasets stored as float32' $C48/$r/forge_run.log $C48/$r/run_case_stdout.log 2>/dev/null | sort -u | head -1 | cut -c1-200)" >> $LOG; done
# ---- 5: case/56 (float32 の格子) の 20 step、新旧と旧の再実行 ----
for r in run_0090_fg6_id32_new run_0091_fg6_id32_old run_0092_fg6_id32_old2; do mkdir $C56/$r && cp $C56/run_0073_lay2cmp_varlen/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json,solverConfig.yaml} $C56/$r/; done
run $C56/run_0090_fg6_id32_new $F_NEW32; run $C56/run_0091_fg6_id32_old $F_OLD32; run $C56/run_0092_fg6_id32_old2 $F_OLD32
echo "== 3 警告 case/56 (float32 の格子): $(grep -h 'geometry datasets stored as float32' $C56/run_0090_fg6_id32_new/forge_run.log $C56/run_0090_fg6_id32_new/run_case_stdout.log 2>/dev/null | sort -u | head -1 | cut -c1-200)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch fg8.done
