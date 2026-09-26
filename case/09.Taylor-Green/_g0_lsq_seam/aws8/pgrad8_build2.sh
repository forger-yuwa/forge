#!/usr/bin/env bash
# #8a R3 の比較相手: bd22376d (f67fe877 との solver 差分はダンプのみ)。build.sh の終了を待ってからクリーンビルド
set -uo pipefail
while ! grep -q ALL_DONE ~/pgrad8/build.out; do sleep 20; done
dir=$HOME/forge-pgrad-bd; sha=bd22376da521335470208b856a4d3ae341adf96b
t0=$(date +%s)
[ -d "$dir/.git" ] || { mkdir -p "$dir" && cd "$dir" && git init -q && git remote add origin git@github.com:forger-yuwa/forge.git; }
cd "$dir"
git fetch -q --depth 1 origin "$sha" && git checkout -q --force FETCH_HEAD || { echo FETCH_FAIL; exit 1; }
git submodule update --init --recursive --depth 1 -q || { echo SUBMOD_FAIL; exit 1; }
echo "== $dir head $(git rev-parse HEAD)"
cd solver_density_cuda && rm -rf .build-native build
CXXFLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" FORGE_CUDA_ARCHITECTURES=86 ./tools/build_native_wsl.sh > ~/pgrad8/build_forge-pgrad-bd.log 2>&1
echo "BUILD_EXIT $dir $? ($(( $(date +%s) - t0 )) s)"
ln -sfn .build-native/relwithdebinfo build
sha256sum build/forge
df -h ~ | tail -1
echo ALL_DONE2
