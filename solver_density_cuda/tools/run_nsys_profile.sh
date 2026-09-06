#!/usr/bin/env bash
set -euo pipefail

image_name="${FORGE_IMAGE_NAME:-forge-solver:cuda-dev}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
solver_dir="$(cd -- "$script_dir/.." && pwd)"
repo_root="$(cd -- "$solver_dir/.." && pwd)"
host_uid="$(id -u)"
host_gid="$(id -g)"
profile_env_default="${FORGE_PROFILE:-1}"

usage() {
  cat <<'EOF'
Usage:
  ./tools/run_nsys_profile.sh <case-dir>
  ./tools/run_nsys_profile.sh <case-dir> <report-prefix>

Examples:
  ./tools/run_nsys_profile.sh /home/sano/work/forge/case/08.bump/run_slau_mach1.65_4pctHt
  ./tools/run_nsys_profile.sh case/08.bump/run_slau_mach1.65_4pctHt bump-short

What it does:
  - Runs forge inside the CUDA Docker image under Nsight Systems.
  - Saves .nsys-rep, .sqlite, and stdout log files under <case-dir>/profiles/nsys/.
  - Enables FORGE_PROFILE=1 by default so the runtime section summary is captured too.

Notes:
  - Run this on the host, not inside the container.
  - <case-dir> must be inside the forge repository and contain solverConfig.yaml.
  - To disable the runtime summary, set FORGE_PROFILE=0 when invoking this script.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 || $# -gt 2 ]]; then
  usage
  if [[ $# -lt 1 || $# -gt 2 ]]; then
    exit 1
  fi
  exit 0
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "[run_nsys_profile] ERROR: docker command not found" >&2
  exit 1
fi

if ! docker image inspect "$image_name" >/dev/null 2>&1; then
  echo "[run_nsys_profile] ERROR: docker image '$image_name' not found" >&2
  echo "[run_nsys_profile] Build it first from solver_density_cuda:" >&2
  echo "  ./tools/docker_build.sh dev" >&2
  exit 2
fi

case_dir_input="$1"
case_dir_host="$(readlink -f -- "$case_dir_input")"
repo_root_real="$(readlink -f -- "$repo_root")"

if [[ ! -d "$case_dir_host" ]]; then
  echo "[run_nsys_profile] ERROR: case directory '$case_dir_input' not found" >&2
  exit 3
fi

case "$case_dir_host" in
  "$repo_root_real"/*) ;;
  *)
    echo "[run_nsys_profile] ERROR: case directory must be inside $repo_root_real" >&2
    exit 4
    ;;
esac

if [[ ! -f "$case_dir_host/solverConfig.yaml" ]]; then
  echo "[run_nsys_profile] ERROR: solverConfig.yaml not found in $case_dir_host" >&2
  exit 5
fi

rel_case_dir="${case_dir_host#"$repo_root_real"/}"
case_dir_container="/workspace/$rel_case_dir"
profile_dir_host="$case_dir_host/profiles/nsys"
mkdir -p "$profile_dir_host"

timestamp="$(date +%Y%m%d_%H%M%S)"
default_prefix="$(basename -- "$case_dir_host")_$timestamp"
report_prefix_name="${2:-$default_prefix}"
report_prefix_host="$profile_dir_host/$report_prefix_name"
report_prefix_container="${case_dir_container}/profiles/nsys/$report_prefix_name"
stdout_log_container="${report_prefix_container}.stdout.log"

docker_args=(
  run --rm --gpus all
  --user "$host_uid:$host_gid"
  -v "$repo_root:/workspace"
  -w "$case_dir_container"
  -e FORGE_PROFILE="$profile_env_default"
  -e FORGE_PROFILE_VERBOSE="${FORGE_PROFILE_VERBOSE:-0}"
  -e HOME=/tmp/forge-home
  -e XDG_CONFIG_HOME=/tmp/forge-home/.config
)

echo "[run_nsys_profile] Image:         $image_name" >&2
echo "[run_nsys_profile] Case dir:      $case_dir_host" >&2
echo "[run_nsys_profile] Report prefix: $report_prefix_host" >&2

exec docker "${docker_args[@]}" "$image_name" bash -lc '
  set -euo pipefail
  if ! command -v nsys >/dev/null 2>&1; then
    echo "[run_nsys_profile] ERROR: nsys is not installed in the container image" >&2
    echo "[run_nsys_profile] Use run_ncu_profile.sh for kernel metrics, or extend Dockerfile.cuda (dev ステージ) with Nsight Systems if you need timelines." >&2
    exit 10
  fi
  mkdir -p "$HOME" "$XDG_CONFIG_HOME"
  blocked_files="$(find . -maxdepth 1 -type f \( -name "res_*.h5" -o -name "res_*.xmf" \) ! -writable -print)"
  if [[ -n "$blocked_files" ]]; then
    echo "[run_nsys_profile] ERROR: existing result files are not writable in $(pwd)" >&2
    echo "$blocked_files" >&2
    echo "[run_nsys_profile] Remove or chown those files before profiling." >&2
    exit 11
  fi
  mkdir -p "$(dirname -- "$1")"
  nsys profile \
    --force-overwrite true \
    --sample=none \
    --cpuctxsw=none \
    --trace=cuda,nvtx,osrt \
    --stats=true \
    --output "$1" \
    /workspace/solver_density_cuda/build/forge 2>&1 | tee "$2"
' _ "$report_prefix_container" "$stdout_log_container"