#!/usr/bin/env bash
set -euo pipefail
# コンテナ内で呼ぶビルドヘルパ。
# HDF5 の include/lib を CMake に渡すほか、以下でビルド時間を詰める:
#   - CUDA アーキテクチャを実機 1 つに絞る (既定は 75;86;89;90 の 4 種 = 約 2 倍の時間)
#   - ccache があればコンパイラランチャに使う (CCACHE_DIR はイメージ側で /workspace/.ccache)
#
# 環境変数:
#   FORGE_CUDA_ARCHITECTURES  例 "86"。未指定なら nvidia-smi から自動判定、取れなければ CMake 既定
#   FORGE_BUILD_DIR           既定 build
#   FORGE_NO_CCACHE=1         ccache を使わない

build_dir="${FORGE_BUILD_DIR:-build}"

if [[ -f "$build_dir/CMakeCache.txt" ]] && ! grep -Fxq "CMAKE_HOME_DIRECTORY:INTERNAL=$PWD" "$build_dir/CMakeCache.txt"; then
	rm -f "$build_dir/CMakeCache.txt"
	rm -rf "$build_dir/CMakeFiles"
fi

cmake_args=(
	-G Ninja
	-DCMAKE_BUILD_TYPE=Release
	-Dhdf5_inc="$HDF5_INC"
	-Dhdf5_libdir="$HDF5_LIBDIR"
)

# --- CUDA アーキテクチャ: 実機の compute capability 1 つに絞る ---
cuda_arch="${FORGE_CUDA_ARCHITECTURES:-}"
if [[ -z "$cuda_arch" ]] && command -v nvidia-smi >/dev/null 2>&1; then
	cuda_arch="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null \
		| head -1 | tr -d ' .' || true)"
	[[ "$cuda_arch" =~ ^[0-9]+$ ]] || cuda_arch=""
fi
if [[ -n "$cuda_arch" ]]; then
	echo "[build] CUDA architectures: $cuda_arch"
	cmake_args+=(-DCMAKE_CUDA_ARCHITECTURES="$cuda_arch")
else
	echo "[build] CUDA architectures: (CMake default; set FORGE_CUDA_ARCHITECTURES to speed up)"
fi

# --- ccache ---
if [[ "${FORGE_NO_CCACHE:-0}" != "1" ]] && command -v ccache >/dev/null 2>&1; then
	echo "[build] ccache: $(command -v ccache) (dir=${CCACHE_DIR:-default})"
	cmake_args+=(
		-DCMAKE_C_COMPILER_LAUNCHER=ccache
		-DCMAKE_CXX_COMPILER_LAUNCHER=ccache
		-DCMAKE_CUDA_COMPILER_LAUNCHER=ccache
	)
fi

mkdir -p "$build_dir"
cd "$build_dir"
build_jobs="${CMAKE_BUILD_PARALLEL_LEVEL:-$(nproc)}"
cmake "${cmake_args[@]}" ..
cmake --build . -j"$build_jobs"
