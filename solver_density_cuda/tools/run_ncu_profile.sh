#!/usr/bin/env bash
set -euo pipefail

image_name="${FORGE_IMAGE_NAME:-forge-solver:cuda-dev}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
solver_dir="$(cd -- "$script_dir/.." && pwd)"
repo_root="$(cd -- "$solver_dir/.." && pwd)"
host_uid="$(id -u)"
host_gid="$(id -g)"
default_kernel_regex='SLAU_d|AUSMp_d|AUSMp_UP_d|ROE_d|KEEP_SLAU_d|KEEP_FVS_d|convectiveFlux_boundary_d|viscousFlux_d|viscousFlux_wall_d|calcGradient_1_d|calcGradient_2_d|limiter_r1_d|ducrosSensor_d|runge_kutta_exp_d|runge_kutta_exp_4th_d|gasProperties_d|WALE_d'

usage() {
  cat <<'EOF'
Usage:
  ./tools/run_ncu_profile.sh <case-dir>
  ./tools/run_ncu_profile.sh <case-dir> <kernel-regex>
  ./tools/run_ncu_profile.sh <case-dir> <kernel-regex> <report-prefix>

Examples:
  ./tools/run_ncu_profile.sh /home/sano/work/forge/case/08.bump/run_slau_mach1.65_4pctHt
  ./tools/run_ncu_profile.sh case/08.bump/run_slau_mach1.65_4pctHt 'SLAU_d|viscousFlux_d' bump-flux

What it does:
  - Runs forge inside the CUDA Docker image under Nsight Compute.
  - Saves an .ncu-rep file and stdout log under <case-dir>/profiles/ncu/.
  - Uses a temporary profiling case directory under profiles/ncu/ so the source case is not modified.
  - Collects sections aimed at GPU efficiency: LaunchStats, Occupancy, SchedulerStats,
    SpeedOfLight, and MemoryWorkloadAnalysis.

Notes:
  - Run this on the host, not inside the container.
  - <case-dir> must be inside the forge repository and contain solverConfig.yaml.
  - The default kernel regex targets the main compute-heavy kernels in this solver.
  - Override FORGE_NCU_NSTEP and FORGE_NCU_OUTSTEP_INTERVAL to shorten profiling runs.
  - Override NCU_TARGET_PROCESSES, NCU_LAUNCH_COUNT, or NCU_IMPORT_SOURCE via env vars if needed.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 || $# -gt 3 ]]; then
  usage
  if [[ $# -lt 1 || $# -gt 3 ]]; then
    exit 1
  fi
  exit 0
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "[run_ncu_profile] ERROR: docker command not found" >&2
  exit 1
fi

if ! docker image inspect "$image_name" >/dev/null 2>&1; then
  echo "[run_ncu_profile] ERROR: docker image '$image_name' not found" >&2
  echo "[run_ncu_profile] Build it first from solver_density_cuda:" >&2
  echo "  ./tools/docker_build.sh dev" >&2
  exit 2
fi

case_dir_input="$1"
case_dir_host="$(readlink -f -- "$case_dir_input")"
repo_root_real="$(readlink -f -- "$repo_root")"

if [[ ! -d "$case_dir_host" ]]; then
  echo "[run_ncu_profile] ERROR: case directory '$case_dir_input' not found" >&2
  exit 3
fi

case "$case_dir_host" in
  "$repo_root_real"/*) ;;
  *)
    echo "[run_ncu_profile] ERROR: case directory must be inside $repo_root_real" >&2
    exit 4
    ;;
esac

if [[ ! -f "$case_dir_host/solverConfig.yaml" ]]; then
  echo "[run_ncu_profile] ERROR: solverConfig.yaml not found in $case_dir_host" >&2
  exit 5
fi

kernel_regex="${2:-$default_kernel_regex}"
rel_case_dir="${case_dir_host#"$repo_root_real"/}"
profile_dir_host="$case_dir_host/profiles/ncu"
mkdir -p "$profile_dir_host"

timestamp="$(date +%Y%m%d_%H%M%S)"
default_prefix="$(basename -- "$case_dir_host")_$timestamp"
report_prefix_name="${3:-$default_prefix}"
report_prefix_host="$profile_dir_host/$report_prefix_name"
temp_case_host="$profile_dir_host/${report_prefix_name}.case"
rm -rf "$temp_case_host"
mkdir -p "$temp_case_host"

while IFS= read -r source_file; do
  base_name="$(basename -- "$source_file")"
  if [[ "$base_name" == "solverConfig.yaml" ]]; then
    cp "$source_file" "$temp_case_host/solverConfig.yaml"
    continue
  fi

  link_target="$(python3 - "$temp_case_host" "$source_file" <<'PY'
from pathlib import Path
import os
import sys

temp_case = Path(sys.argv[1])
source = Path(sys.argv[2])
print(os.path.relpath(source, temp_case))
PY
)"
  ln -s "$link_target" "$temp_case_host/$base_name"
done < <(find "$case_dir_host" -maxdepth 1 -type f \
  ! -name 'res_*.h5' \
  ! -name 'res_*.xmf' \
  ! -name '*.ncu-rep' \
  ! -name '*.nsys-rep' \
  ! -name '*.sqlite' \
  ! -name '*.stdout.log' \
  -print | sort)

if [[ ! -f "$temp_case_host/solverConfig.yaml" ]]; then
  echo "[run_ncu_profile] ERROR: failed to stage solverConfig.yaml into $temp_case_host" >&2
  exit 6
fi

if [[ -n "${FORGE_NCU_NSTEP:-}" ]]; then
  python3 - "$temp_case_host/solverConfig.yaml" "$FORGE_NCU_NSTEP" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
nstep = sys.argv[2]
text = path.read_text()
updated, count = re.subn(r'(^\s*nStepOuter:\s*)\S+', rf'\g<1>{nstep}', text, count=1, flags=re.MULTILINE)
if count != 1:
  raise SystemExit('[run_ncu_profile] ERROR: could not find nStepOuter in solverConfig.yaml')
path.write_text(updated)
PY
fi

if [[ -n "${FORGE_NCU_OUTSTEP_INTERVAL:-}" ]]; then
  python3 - "$temp_case_host/solverConfig.yaml" "$FORGE_NCU_OUTSTEP_INTERVAL" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
interval = sys.argv[2]
text = path.read_text()
updated, count = re.subn(r'(^\s*outStepInterval:\s*)\S+', rf'\g<1>{interval}', text, count=1, flags=re.MULTILINE)
if count != 1:
    raise SystemExit('[run_ncu_profile] ERROR: could not find outStepInterval in solverConfig.yaml')
path.write_text(updated)
PY
fi

temp_case_rel="${temp_case_host#"$repo_root_real"/}"
case_dir_container="/workspace/$temp_case_rel"
report_prefix_container="${case_dir_container}/$report_prefix_name"
stdout_log_container="${report_prefix_container}.stdout.log"

docker_args=(
  run --rm --gpus all
  --user "$host_uid:$host_gid"
  -v "$repo_root:/workspace"
  -w "$case_dir_container"
  -e HOME=/tmp/forge-home
  -e XDG_CONFIG_HOME=/tmp/forge-home/.config
)

echo "[run_ncu_profile] Image:         $image_name" >&2
echo "[run_ncu_profile] Source case:   $case_dir_host" >&2
echo "[run_ncu_profile] Temp case:     $temp_case_host" >&2
echo "[run_ncu_profile] Kernel regex:  $kernel_regex" >&2
echo "[run_ncu_profile] Report prefix: $report_prefix_host" >&2
if [[ -n "${FORGE_NCU_NSTEP:-}" ]]; then
  echo "[run_ncu_profile] Override nStepOuter: ${FORGE_NCU_NSTEP}" >&2
fi
if [[ -n "${FORGE_NCU_OUTSTEP_INTERVAL:-}" ]]; then
  echo "[run_ncu_profile] Override outStepInterval: ${FORGE_NCU_OUTSTEP_INTERVAL}" >&2
fi

exec docker "${docker_args[@]}" "$image_name" bash -lc '
  set -euo pipefail
  if ! command -v ncu >/dev/null 2>&1; then
    echo "[run_ncu_profile] ERROR: ncu is not installed in the container image" >&2
    exit 10
  fi
  mkdir -p "$HOME" "$XDG_CONFIG_HOME"
  blocked_files="$(find . -maxdepth 1 -type f \( -name "res_*.h5" -o -name "res_*.xmf" \) ! -writable -print)"
  if [[ -n "$blocked_files" ]]; then
    echo "[run_ncu_profile] ERROR: existing result files are not writable in $(pwd)" >&2
    echo "$blocked_files" >&2
    echo "[run_ncu_profile] Remove or chown those files before profiling." >&2
    exit 11
  fi
  mkdir -p "$(dirname -- "$1")"
  set +e
  ncu \
    --target-processes "${NCU_TARGET_PROCESSES:-all}" \
    --force-overwrite \
    --kernel-name-base demangled \
    --kernel-name "regex:$2" \
    --launch-count "${NCU_LAUNCH_COUNT:-1}" \
    --import-source "${NCU_IMPORT_SOURCE:-no}" \
    --section LaunchStats \
    --section Occupancy \
    --section SchedulerStats \
    --section SpeedOfLight \
    --section MemoryWorkloadAnalysis \
    --export "$1" \
    /workspace/solver_density_cuda/build/forge 2>&1 | tee "$3"
  ncu_status=${PIPESTATUS[0]}
  set -e

  if grep -q "ERR_NVGPUCTRPERM" "$3"; then
    echo "[run_ncu_profile] ERROR: NVIDIA performance counters are restricted on the host." >&2
    echo "[run_ncu_profile] Ask the host admin to enable GPU profiling counters, or set NVreg_RestrictProfilingToAdminUsers=0 for the NVIDIA kernel module." >&2
    echo "[run_ncu_profile] See https://developer.nvidia.com/ERR_NVGPUCTRPERM" >&2
    exit 12
  fi

  exit "$ncu_status"
' _ "$report_prefix_container" "$kernel_regex" "$stdout_log_container"