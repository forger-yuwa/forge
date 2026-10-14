#!/bin/bash
# plan time_integration-implicit-thermal-jacobian §6 V0 の追加 (2026-10-09、事後): V0 (run_0225〜0232) は 4 組ともビット不一致だった。
# forge は面の流束を atomicAdd で足すので同じバイナリでも再実行でビット一致しない可能性がある → 同じバイナリの再実行の差を測る。
#   (a) 20 step: pt0・dir0 を旧・新それぞれ 2 本追加 (run_0225〜0232 と合わせて各 3 本)
#   (b) 1 step: pt0・dir0 を旧・新それぞれ 2 本 (1 step なら加算順の丸めの差だけが残るはず)
# forge は cold_pair.py / cold_cfl.py の run (= run_case.sh) 経由で回す。AWS の case dir で実行する。
cd "$(dirname "$0")"
OLD=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/forge
NEW=$HOME/forge-thermjac-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES
n=233
one() {  # $1 = run 名, $2 = pt|dir, $3 = old|new, $4 = step 数
  local r=$1 opt="--steps $4 --cfl 4 --out $4 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext"
  [ $2 = dir ] && opt="$opt --line dir"
  if [ $3 = old ]; then
    export FORGE_BIN=$OLD; unset COLD_ALT_BINARY
    python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r $opt > /dev/null && python3 cold_pair.py run $r > $r/cold_pair_run_stdout.log 2>&1
  else
    export FORGE_BIN=$NEW COLD_ALT_BINARY=thermjac_cap_fp64
    python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r $opt > /dev/null && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
  fi
  echo "$r rc=$? $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> v0_repeat.log
  rm -f $r/nozzle.msh
}
for cfg in pt dir; do for b in old new; do for k in b c; do
  one $(printf "run_%04d_v0_%s0_%s_%s" $n $cfg $b $k) $cfg $b 20; n=$((n+1))
done; done; done
for cfg in pt dir; do for b in old new; do for k in a b; do
  one $(printf "run_%04d_v0s1_%s0_%s_%s" $n $cfg $b $k) $cfg $b 1; n=$((n+1))
done; done; done
touch v0_repeat.done
