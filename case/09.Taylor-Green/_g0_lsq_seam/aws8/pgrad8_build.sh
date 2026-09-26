#!/usr/bin/env bash
# plan boundary-node-periodic-gradient-fix §5.1 #8a/#8d 用のビルド (全部クリーン)。
#   ~/forge-pgrad-new      : f67fe877 (新 + ダンプ)
#   ~/forge-pgrad-olddump  : 1266aba1 + olddump.patch (旧 + 同じダンプ)
#   ~/forge-pgrad-wrong    : f67fe877 + wrong.patch (除外条件を periodicNodeActive に戻した誤り版 + ダンプ)
#   ~/forge-pgrad-m3       : 1266aba1 + m3only.patch (F1 初期値の移動だけ)
set -uo pipefail
NEW=f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
OLD=1266aba1a820a53b96832f88f69ec330d9b1c12d
build_one() {
  local dir=$1 sha=$2 patch=${3:-}
  local t0=$(date +%s)
  if [ ! -d "$dir/.git" ]; then
    mkdir -p "$dir" && cd "$dir" && git init -q && git remote add origin git@github.com:forger-yuwa/forge.git
  fi
  cd "$dir"
  git fetch -q --depth 1 origin "$sha" && git checkout -q --force FETCH_HEAD || { echo "FETCH_FAIL $dir"; return 1; }
  git submodule update --init --recursive --depth 1 -q || { echo "SUBMOD_FAIL $dir"; return 1; }
  if [ -n "$patch" ]; then git apply --check "$patch" && git apply "$patch" || { echo "PATCH_FAIL $dir"; return 1; }; fi
  echo "== $dir head $(git rev-parse HEAD) patch=${patch:-none}"; git diff --stat | cat
  cd solver_density_cuda
  rm -rf .build-native build
  CXXFLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" FORGE_CUDA_ARCHITECTURES=86 ./tools/build_native_wsl.sh > ~/pgrad8/build_$(basename $dir).log 2>&1
  echo "BUILD_EXIT $dir $? ($(( $(date +%s) - t0 )) s)"
  ln -sfn .build-native/relwithdebinfo build
  sha256sum build/forge build/convertGmshToForge
}
build_one $HOME/forge-pgrad-new     $NEW
build_one $HOME/forge-pgrad-olddump $OLD $HOME/pgrad8/olddump.patch
build_one $HOME/forge-pgrad-wrong   $NEW $HOME/pgrad8/wrong.patch
build_one $HOME/forge-pgrad-m3      $OLD $HOME/pgrad8/m3only.patch
sha256sum $HOME/forge-pgrad-old/solver_density_cuda/build/forge
df -h ~ | tail -1
echo ALL_DONE
