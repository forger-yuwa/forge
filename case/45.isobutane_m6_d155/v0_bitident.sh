#!/bin/bash
# plan time_integration-implicit-thermal-jacobian §6 V0: キー無し (既定) で旧バイナリ (65be5e28…) と新バイナリ (35e498b1…) の res_20 がビット一致するか。
# 設定: point cfl 4 と directional cfl 4 (run_0203 と同じ) × implicitSolvePrecision 0/1。run_0183 の res_100000 から 20 step。
# forge は cold_pair.py / cold_cfl.py の run (= run_case.sh) 経由で回す。AWS の case dir で実行する。
cd "$(dirname "$0")"
OLD=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/forge
NEW=$HOME/forge-thermjac-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CUDA_BLOCKSIZE=128 FORGE_CONVERTER=$PWD/conv_tolerant.sh
unset FORGE_ALLOW_UNVERIFIED_SPECIES
n=225
for cfg in pt dir; do for isp in 0 1; do for b in old new; do
  r=$(printf "run_%04d_v0_%s%d_%s" $n $cfg $isp $b); n=$((n+1))
  opt="--steps 20 --cfl 4 --out 20 --limiter-ref-from run_0183_ns_coldmesh_tw300_ext"
  [ $cfg = dir ] && opt="$opt --line dir"
  [ $isp = 1 ] && opt="$opt --isp 1"
  if [ $b = old ]; then
    export FORGE_BIN=$OLD; unset COLD_ALT_BINARY
    python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r $opt > /dev/null && python3 cold_pair.py run $r > $r/cold_pair_run_stdout.log 2>&1
  else
    export FORGE_BIN=$NEW COLD_ALT_BINARY=thermjac_cap_fp64
    python3 cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext $r $opt > /dev/null && python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1
  fi
  echo "$r rc=$? $(grep forge_sha256 $r/RUN_PROVENANCE.txt | cut -c15-30)" >> v0_bitident.log
done; done; done
python3 - <<'PY' >> v0_bitident.log
import glob
import json

import h5py
import numpy as np

out = {}
for cfg in ("pt", "dir"):
    for isp in (0, 1):
        o = glob.glob(f"run_*_v0_{cfg}{isp}_old/res_20.h5"); nw = glob.glob(f"run_*_v0_{cfg}{isp}_new/res_20.h5")
        if not o or not nw:
            out[f"{cfg}{isp}"] = "res_20 が無い"
            continue
        A = h5py.File(o[0])["VALUE"]; B = h5py.File(nw[0])["VALUE"]
        keys = [k for k in A.keys() if k in B]
        diff = {k: float(np.max(np.abs(A[k][:] - B[k][:]))) for k in keys}
        out[f"{cfg}{isp}"] = {"bit_identical": all(np.array_equal(A[k][:], B[k][:]) for k in keys), "n_fields": len(keys),
                              "max_abs_diff_nonzero": {k: v for k, v in diff.items() if v != 0.0}}
json.dump(out, open("_band_ab/cold_pair/V0_bitident.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(out, ensure_ascii=False))
PY
touch v0_bitident.done
