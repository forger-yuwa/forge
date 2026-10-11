#!/bin/bash
# plan architecture-float-state-double-geometry §5.1 #22・§6.52 の準備: G0c (make_g0c.py が作った .msh) を B で変換し、検査し、共通の起点を作る。
# 1. 変換: ~/forge-fgeom7-fp64 の convertGmshToForge を、run と同じ設定一式 (_tr/*.yaml) を置いた _g0c/ で回す (G0 の _conv9/c45/new64 と同じ変換器)。
# 2. 検査: hoop_verify_conv.py check --rw (rSurfVect の独立な作り直しとの照合)、check_mesh_quality.py --ar-max 10000 (ユーザ決定: 診断に限り AR ≤ 10,000)。
# 3. 初期値: G0 の共通の起点 _tr/q32_init.h5 の場を interp_field.py で G0c に移し、保存量 9 つを float32 に丸めて double に戻す → _g0c/q32_init.h5。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/g0c_prep.log; TL=~/forge-wallfit/solver_density_cuda/tools
CONV=~/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge; F64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45
echo "== 開始 $(date -Is) 変換器 $(sha256sum $CONV | cut -c1-16) msh $(sha256sum _g0c/nozzle.msh | cut -c1-16)" > $LOG
cp _tr/*.yaml _g0c/
( cd _g0c && $CONV nozzle.msh nozzle.h5 > conv.log 2>&1 ); rc=$?
if [ $rc -ne 0 ] && grep -q "Write Input HDF5" _g0c/conv.log && grep -q "GPUassert: invalid argument" _g0c/conv.log; then rc=0; fi
echo "変換 rc=$rc $(sha256sum _g0c/nozzle.h5 | cut -c1-16)" >> $LOG
[ $rc = 0 ] || exit 1
echo "== 閉性 (hoop_verify_conv.py check --rw)" >> $LOG; python3 hoop_verify_conv.py check _g0c/nozzle.h5 --rw >> $LOG 2>&1; echo "  rc=$?" >> $LOG
echo "== 格子の品質 (--ar-max 10000)" >> $LOG; python3 $TL/check_mesh_quality.py _g0c/nozzle.h5 --ar-max 10000 >> $LOG 2>&1; echo "  rc=$?" >> $LOG
echo "== 格子の品質 (既定・参考)" >> $LOG; python3 $TL/check_mesh_quality.py _g0c/nozzle.h5 2>&1 | grep -E "VERDICT|aspect|AR" | head -5 >> $LOG
# 初期値
cp _g0c/nozzle.h5 _g0c/q32_init.h5
python3 - _tr/q32_init.h5 _g0c/q32_init.h5 >> $LOG 2>&1 <<'EOF'
import sys, re, h5py
with h5py.File(sys.argv[1], "r") as s, h5py.File(sys.argv[2], "r+") as d:
    ks = sorted(k for k in s["VALUE"] if re.fullmatch(r"roY\d+", k))
    for k in ks:
        if k not in d["VALUE"]: d["VALUE"].create_dataset(k, shape=(d["VALUE/ro"].shape[0],), dtype=s["VALUE/" + k].dtype)
    print("[g0c] 宛先に roY* を作った:", ks)
EOF
python3 $TL/interp_field.py _tr/q32_init.h5 _g0c/q32_init.h5 --dst-run _g0c --forge $F64 > _g0c/interp.log 2>&1; echo "interp_field rc=$? $(grep -h VERDICT _g0c/interp.log | tail -1)" >> $LOG
python3 - _g0c/q32_init.h5 >> $LOG 2>&1 <<'EOF'
import sys, h5py, numpy as np
Q = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
with h5py.File(sys.argv[1], "r+") as h:
    V = h["VALUE"]; n = 0
    for k in Q:
        a = np.asarray(V[k][:])
        if a.dtype != np.float64: raise SystemExit(f"{k} が float64 でない ({a.dtype})")
        if not np.all(np.isfinite(a)): raise SystemExit(f"{k} に非有限")
        b = a.astype(np.float32).astype(np.float64); n += int(np.count_nonzero(a != b)); V[k][...] = b
    print(f"[Q32] G0c の保存量 9 つを float32 に丸めた (変わった値 {n} 個)")
EOF
echo "共通の入力 _g0c/q32_init.h5 sha $(sha256sum _g0c/q32_init.h5 | cut -c1-16)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG; touch g0c_prep.done
