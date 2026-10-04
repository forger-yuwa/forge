#!/bin/bash
# w/d 掃引: 実験と理論が一致する広い側 (0.383) で forge を検証する
set -e
cd "$(dirname "$0")"
for spec in "0.211:B1:t1_wd0211_v2:run_0008_T1_wd0211" "0.383:C1:t1_wd0383_v2:run_0009_T1_wd0383"; do
  IFS=: read wd series mesh run <<< "$spec"
  echo "=== $run (w/d=$wd, series $series) ==="
  python3 tools/make_case.py --run "$run" --mesh "$mesh" --series "$series" --main-steps 20000 --out-int 5000
  # 本段を継続して評価用の時系列を取る
  long="${run}_long"
  mkdir -p "$long"
  cp "$run"/{mesh.h5,bcondConfig.yaml,solverConfig.yaml,species_db.yaml,probe.yaml,case_setup.json} "$long/"
  echo "$run (main 段の収束場)" > "$long/CONTINUED_FROM"
  python3 - "$long" <<'PY'
import sys,pathlib,re
p=pathlib.Path(sys.argv[1])/'solverConfig.yaml'; s=p.read_text()
s=re.sub(r'nStepOuter: \d+','nStepOuter: 60000',s); s=re.sub(r'outStepInterval: \d+','outStepInterval: 5000',s)
p.write_text(s)
PY
  LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial ../../solver_density_cuda/tools/run_case.sh "$long" >/dev/null 2>&1 || true
  python3 tools/cavity_eval.py "$long" --series | tail -6
done
