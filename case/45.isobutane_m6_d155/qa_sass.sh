#!/bin/bash
# plan architecture-float-state-double-geometry §6.34 の前提 1 の SASS の照合 (qa_build.sh から呼ぶ。単独でも回せる)。
# 使い方: bash qa_sass.sh <commit> <単体試験の合否 True|False>。基準 ~/forge-fgeom7-f32、新しい ~/forge-fgeom9-f32。
# cuobjdump は A に無い (CUDA 13.2 の導入に含まれない) ので ~/cuobj_venv の pip の nvidia-cuda-cuobjdump・nvidia-cuda-nvdisasm 13.2 を使う。
set -u
SHA=$1; ut=$2
C45=~/forge-wallfit/case/45.isobutane_m6_d155; PRE=$C45/_band_ab/cold_pair/qa_pre.json; LOG=$C45/qa_build.log
S=~/forge-fgeom9-f32/solver_density_cuda; OLD=~/forge-fgeom7-f32/solver_density_cuda
W=~/qa_sass; rm -rf $W/variants; mkdir -p $W/variants
CUOBJ=$(find ~/cuobj_venv -name cuobjdump -type f | head -1); [ -x "$CUOBJ" ] || { echo "cuobjdump が無い"; exit 1; }
export PATH=$(dirname $CUOBJ):$PATH   # cuobjdump -sass は nvdisasm を PATH から呼ぶ (同じ venv の nvidia-cuda-nvdisasm 13.2)
echo "cuobjdump: $CUOBJ $($CUOBJ --version | tail -1)"
$CUOBJ -sass $OLD/build/forge > $W/base.sass; $CUOBJ -sass $S/build/forge > $W/new.sass
python3 -c "import json; json.dump({'commit': '$SHA', 'unit_test_pass': $ut, 'forge_sha256': '$(sha256sum $S/build/forge | cut -c1-64)', 'base_forge_sha256': '$(sha256sum $OLD/build/forge | cut -c1-64)'}, open('$PRE', 'w'), ensure_ascii=False, indent=1)"
python3 $C45/qa_sass_cmp.py $W/base.sass $W/new.sass $PRE > $W/cmp1.json
# 一致しない関数があれば、その関数を含む基準のソースの .cu を同じフラグで 3 回再コンパイルして揺れの候補を作り、照合し直す
FN=$(python3 -c "import json; d=json.load(open('$W/cmp1.json')); print(' '.join(d['differ']))")
if [ -n "$FN" ]; then
  FM=$OLD/build/cuda_forge/CMakeFiles/cuda_forge.dir/flags.make
  DEF=$(sed -n 's/^CUDA_DEFINES = //p' $FM); INC=$(sed -n 's/^CUDA_INCLUDES = //p' $FM); FL=$(sed -n 's/^CUDA_FLAGS = //p' $FM)
  for f in $FN; do
    nm=$(python3 -c "import re,sys; f=sys.argv[1]; m=re.search(r'_(\d+)_([A-Za-z0-9_]+?)_cu_', f); n=re.match(r'_Z(\d+)', f)
print(m.group(2) if m else f[n.end():n.end()+int(n.group(1))])" $f)
    src=$(grep -lE "__global__ void $nm\b|^$nm\.cu" $OLD/cuda_forge/*.cu | head -1); [ -z "$src" ] && [ -f $OLD/cuda_forge/$nm.cu ] && src=$OLD/cuda_forge/$nm.cu
    echo "揺れの確認: $nm → $src" >> $LOG
    [ -n "$src" ] || continue
    for i in 1 2 3; do
      (cd $OLD/build/cuda_forge && eval /usr/local/cuda/bin/nvcc $DEF $INC $FL -x cu -c $src -o $W/variants/${nm}_base_$i.o) >> $LOG 2>&1
      $CUOBJ -sass $W/variants/${nm}_base_$i.o > $W/variants/${nm}_base_$i.sass
    done
  done
  python3 $C45/qa_sass_cmp.py $W/base.sass $W/new.sass $PRE --variants $W/variants > $W/cmp2.json
fi
python3 -c "import json; d=json.load(open('$PRE')); print('単体試験', d['unit_test_pass'], '照合で残った不一致', d['sass_mismatch_except_axis'], '揺れとみなした', list(d.get('differ_resolved_as_nvcc_variant', {})), '関数の数', d['n_base'], d['n_new'])"
