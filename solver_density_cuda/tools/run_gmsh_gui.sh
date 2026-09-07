#!/usr/bin/env bash
set -euo pipefail

image_name="${FORGE_IMAGE_NAME:-forge-solver:cuda-dev}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
solver_dir="$(cd -- "$script_dir/.." && pwd)"
repo_root="$(cd -- "$solver_dir/.." && pwd)"
host_uid="$(id -u)"
host_gid="$(id -g)"
wayland_display_name="${WAYLAND_DISPLAY:-wayland-0}"
use_wslg_runtime=0

usage() {
  cat <<'EOF'
Usage:
  ./tools/run_gmsh_gui.sh
  ./tools/run_gmsh_gui.sh <path-to-.geo-or-.msh>

Notes:
  - Run this on the host, not inside the container.
  - The target file must be inside the forge repository.
  - GUI display requires WSLg or X11/Wayland forwarding on the host.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "[run_gmsh_gui] ERROR: docker command not found" >&2
  exit 1
fi

if ! docker image inspect "$image_name" >/dev/null 2>&1; then
  echo "[run_gmsh_gui] ERROR: docker image '$image_name' not found" >&2
  echo "[run_gmsh_gui] Build it first from solver_density_cuda:" >&2
  echo "  ./tools/docker_build.sh dev" >&2
  exit 2
fi

target_in_container=""
workdir_in_container="/workspace/solver_density_cuda"

if [[ $# -ge 1 ]]; then
  target_path="$(readlink -f -- "$1")"
  repo_root_real="$(readlink -f -- "$repo_root")"
  case "$target_path" in
    "$repo_root_real"/*) ;;
    *)
      echo "[run_gmsh_gui] ERROR: target must be inside $repo_root_real" >&2
      exit 3
      ;;
  esac

  rel_path="${target_path#"$repo_root_real"/}"
  target_in_container="/workspace/$rel_path"
  workdir_in_container="$(dirname -- "$target_in_container")"
fi

docker_args=(
  run --rm -it --gpus all
  --user "$host_uid:$host_gid"
  -v "$repo_root:/workspace"
  -w "$workdir_in_container"
  -e DISPLAY
  -e WAYLAND_DISPLAY
  -e PULSE_SERVER
)

if [[ -S "/mnt/wslg/runtime-dir/$wayland_display_name" ]]; then
  docker_args+=(
    --tmpfs "/tmp/xdg-runtime:mode=700,uid=$host_uid,gid=$host_gid"
    -e XDG_RUNTIME_DIR=/tmp/xdg-runtime
    -e WAYLAND_DISPLAY="$wayland_display_name"
    -v "/mnt/wslg/runtime-dir/$wayland_display_name:/tmp/xdg-runtime/$wayland_display_name"
  )

  if [[ -S /mnt/wslg/PulseServer ]]; then
    docker_args+=(
      -e PULSE_SERVER=unix:/tmp/xdg-runtime/PulseServer
      -v /mnt/wslg/PulseServer:/tmp/xdg-runtime/PulseServer
    )
  fi

  use_wslg_runtime=1
elif [[ -n "${XDG_RUNTIME_DIR:-}" && -d "$XDG_RUNTIME_DIR" ]]; then
  docker_args+=( -e XDG_RUNTIME_DIR -v "$XDG_RUNTIME_DIR:$XDG_RUNTIME_DIR" )
fi

if [[ -d /tmp/.X11-unix ]]; then
  docker_args+=( -v /tmp/.X11-unix:/tmp/.X11-unix )
fi

if [[ -d /mnt/wslg && "$use_wslg_runtime" -eq 0 ]]; then
  docker_args+=( -v /mnt/wslg:/mnt/wslg )
fi

echo "[run_gmsh_gui] Image:    $image_name" >&2
echo "[run_gmsh_gui] Repo:     $repo_root" >&2
echo "[run_gmsh_gui] Workdir:  $workdir_in_container" >&2
if [[ -n "$target_in_container" ]]; then
  echo "[run_gmsh_gui] Opening:  $target_in_container" >&2
fi

if [[ -n "$target_in_container" ]]; then
  exec docker "${docker_args[@]}" "$image_name" gmsh "$target_in_container"
fi

exec docker "${docker_args[@]}" "$image_name" gmsh