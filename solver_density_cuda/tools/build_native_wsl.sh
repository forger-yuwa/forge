#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
solver_dir="$(cd -- "$script_dir/.." && pwd)"
# Keep the native WSL build outside the Docker-created build/ tree, which is
# often root-owned when the project has been built inside containers.
build_dir="${FORGE_NATIVE_BUILD_DIR:-$solver_dir/.build-native/relwithdebinfo}"
build_type="${FORGE_NATIVE_BUILD_TYPE:-RelWithDebInfo}"
nvcc_path="${CUDACXX:-$(command -v nvcc || true)}"
cuda_architectures="${FORGE_CUDA_ARCHITECTURES:-}"
default_cc="${CC:-}"
default_cxx="${CXX:-}"

if [[ -z "$default_cc" ]]; then
  if command -v gcc-12 >/dev/null 2>&1; then
    default_cc="$(command -v gcc-12)"
  else
    default_cc="/usr/bin/gcc"
  fi
fi

if [[ -z "$default_cxx" ]]; then
  if command -v g++-12 >/dev/null 2>&1; then
    default_cxx="$(command -v g++-12)"
  else
    default_cxx="/usr/bin/g++"
  fi
fi

if [[ -z "$nvcc_path" ]]; then
  echo "[build_native_wsl] ERROR: nvcc not found. Install the CUDA toolkit for WSL before using this script." >&2
  exit 1
fi

if [[ ! -f "$solver_dir/third_party/HighFive/include/highfive/highfive.hpp" ]]; then
  echo "[build_native_wsl] ERROR: HighFive headers not found under third_party/HighFive." >&2
  echo "[build_native_wsl] Run: git submodule update --init --recursive" >&2
  exit 2
fi

mkdir -p "$build_dir"

echo "[build_native_wsl] solver dir: $solver_dir"
echo "[build_native_wsl] build dir:  $build_dir"
echo "[build_native_wsl] nvcc:       $nvcc_path"
echo "[build_native_wsl] cc:         $default_cc"
echo "[build_native_wsl] cxx:        $default_cxx"

cmake_args=(
  -S "$solver_dir"
  -B "$build_dir"
  -G Ninja
  -DCMAKE_BUILD_TYPE="$build_type"
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
  -DCMAKE_C_COMPILER="$default_cc"
  -DCMAKE_CXX_COMPILER="$default_cxx"
  -DCMAKE_CUDA_COMPILER="$nvcc_path"
  -DCMAKE_CUDA_HOST_COMPILER="$default_cxx"
)

if [[ -n "$cuda_architectures" ]]; then
  echo "[build_native_wsl] cuda arch:  $cuda_architectures"
  cmake_args+=(-DCMAKE_CUDA_ARCHITECTURES="$cuda_architectures")
fi

# CUDA 13.x は thrust・cub を include/cccl/ に移した (procedures/cloud-aws-gpu.md「移設で必ず踏む 4 点」の 2)。
# .cpp も vector_types.h を読むので、toolkit の include を C++ 側にも渡す (runfile 版の CUDA は /usr/include に無い)。
cuda_root="$(cd -- "$(dirname -- "$(readlink -f "$nvcc_path")")/.." && pwd)"
cuda_major="$("$nvcc_path" --version | sed -n 's/.*release \([0-9]*\)\..*/\1/p' | head -n 1)"
if [[ -n "$cuda_major" && "$cuda_major" -ge 13 && -d "$cuda_root/include/cccl" ]]; then
  echo "[build_native_wsl] CUDA $cuda_major: add $cuda_root/include/cccl"
  cmake_args+=(-DCMAKE_CUDA_FLAGS="-I$cuda_root/include/cccl"
               -DCMAKE_CXX_FLAGS="-I$cuda_root/include -I$cuda_root/include/cccl")
fi

cmake "${cmake_args[@]}"

cmake --build "$build_dir" -j "${FORGE_BUILD_JOBS:-$(nproc)}"

echo "[build_native_wsl] Built: $build_dir/forge"
