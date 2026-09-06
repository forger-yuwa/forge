#!/usr/bin/env bash
set -euo pipefail
# Dockerfile.cuda の各ステージをビルドするヘルパ。
#
#   ./tools/docker_build.sh base    -> forge-solver:cuda-base   (ビルド最小構成)
#   ./tools/docker_build.sh cloud   -> forge-solver:cuda-cloud  (+ Python / gmsh)
#   ./tools/docker_build.sh dev     -> forge-solver:cuda-dev    (+ ParaView / GUI)  [既定]
#   ./tools/docker_build.sh all     -> 上記 3 つを base から順に
#
# 環境変数:
#   FORGE_IMAGE_PREFIX  タグの接頭辞 (既定: forge-solver)
#   DOCKER_BUILD_ARGS   docker build に追加で渡す引数 (例: --no-cache)
#
# Dockerfile.cuda は BuildKit 固有の記法を使っていないので legacy builder でもそのまま通る。

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
solver_dir="$(cd -- "$script_dir/.." && pwd)"
prefix="${FORGE_IMAGE_PREFIX:-forge-solver}"
target="${1:-dev}"

dockerfile="$solver_dir/Dockerfile.cuda"

build_stage() {
	local stage="$1"
	echo "[docker_build] target=${stage} tag=${prefix}:cuda-${stage}"
	docker build \
		-f "$dockerfile" \
		--target "$stage" \
		-t "${prefix}:cuda-${stage}" \
		${DOCKER_BUILD_ARGS:-} \
		"$solver_dir"
}

case "$target" in
	base|cloud|dev) build_stage "$target" ;;
	all)            build_stage base && build_stage cloud && build_stage dev ;;
	*) echo "usage: $0 [base|cloud|dev|all]" >&2; exit 2 ;;
esac

echo "[docker_build] done."
