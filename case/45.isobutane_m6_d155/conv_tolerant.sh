#!/bin/bash
# convertGmshToForge の終了時 GPUassert (cudaFree invalid argument, AWS g5 で既知: 出力 h5 は完全) だけを許容するラッパー。
# 変換が "Write Input HDF5" まで進み、出力 h5 が開けて /MESH/COORD を持つときに限り終了コード 0 にする。それ以外は元の終了コードを返す。
# usage: FORGE_CONVERTER=<このファイル> (REAL_CONVERTER で本体を指定。既定は ~/forge-integ の統合ビルド)
REAL=${REAL_CONVERTER:-$HOME/forge-integ/solver_density_cuda/build/convertGmshToForge}
log=$(mktemp)
"$REAL" "$@" 2>&1 | tee "$log"; rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ] && grep -q "Write Input HDF5" "$log" && grep -q "GPUassert: invalid argument" "$log"; then
  out="${2%.h5}.h5"
  if python3 -c "import h5py,sys; f=h5py.File(sys.argv[1]); assert f['/MESH/COORD'].size>0" "$out" 2>/dev/null; then
    echo "[conv_tolerant] 終了時 GPUassert を許容 (出力 $out は /MESH/COORD を持つ)"; rm -f "$log"; exit 0
  fi
fi
rm -f "$log"; exit "$rc"
