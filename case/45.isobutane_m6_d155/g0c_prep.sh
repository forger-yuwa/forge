#!/bin/bash
# plan architecture-float-state-double-geometry §5.1 #22・§6.52 の準備: G0c (make_g0c.py が作った .msh) を B で変換し、検査し、共通の起点を作る。
# 1. 変換: ~/forge-fgeom7-fp64 の convertGmshToForge を、run と同じ設定一式 (_tr/*.yaml) を置いた _g0c/ で回す (G0 の _conv9/c45/new64 と同じ変換器)。
# 2. 検査: hoop_verify_conv.py check --rw (rSurfVect の独立な作り直しとの照合)、check_mesh_quality.py --ar-max 10000 (ユーザ決定: 診断に限り AR ≤ 10,000)。
# 3. 初期値: G0 の共通の起点 _tr/q32_init.h5 の場を interp_field.py で G0c に移し、保存量 9 つを float32 に丸めて double に戻す → _g0c/q32_init.h5。
# 4. 照合: g0c_check.py (実際の G0 の入力と、区間の外の座標・9 保存量のビット一致、区間の中の最近傍の対応)。
# 2026-10-11 (codex diagnose 指摘 2): 各段で失敗したら止める。最後まで通ったときだけ g0c_prep.done を作る。
# 使い方: bash g0c_prep.sh <nozzle.msh の sha256 の先頭 16 桁>  (_g0c/nozzle.msh を置き、_g0c/nozzle.h5 が無い状態で回す)
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/g0c_prep.log; TL=~/forge-wallfit/solver_density_cuda/tools
CONV=~/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge; F64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45
die() { echo "中止: $*" >> $LOG; echo "中止: $*"; exit 1; }
rm -f g0c_prep.done
echo "== 開始 $(date -Is) 変換器 $(sha256sum $CONV | cut -c1-16) msh $(sha256sum _g0c/nozzle.msh | cut -c1-16)" > $LOG
[ "$(sha256sum _g0c/nozzle.msh | cut -c1-16)" = "${1:-none}" ] || die ".msh の sha が引数と違う"
[ -e _g0c/nozzle.h5 ] && die "_g0c/nozzle.h5 が既にある (古い変換を消してから回す)"
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || die "G0 の共通の入力が §6.26 と違う"
cp _tr/*.yaml _g0c/ || die "設定一式の複製"
( cd _g0c && $CONV nozzle.msh nozzle.h5 > conv.log 2>&1 ); rc=$?
# 変換器は HDF5 を書いた後の GPU の後片付けで落ちることがある (既知、forge-aws-run)。書き終えていれば通す
if [ $rc -ne 0 ] && grep -q "Write Input HDF5" _g0c/conv.log && grep -q "GPUassert: invalid argument" _g0c/conv.log; then rc=0; fi
echo "変換 rc=$rc $(sha256sum _g0c/nozzle.h5 | cut -c1-16)" >> $LOG
[ $rc = 0 ] || die "変換"
echo "== 閉性 (hoop_verify_conv.py check --rw)" >> $LOG; python3 hoop_verify_conv.py check _g0c/nozzle.h5 --rw >> $LOG 2>&1; rc=$?; echo "  rc=$rc" >> $LOG
[ $rc = 0 ] || die "閉性の検査"
echo "== 格子の品質 (--ar-max 10000)" >> $LOG; python3 $TL/check_mesh_quality.py _g0c/nozzle.h5 --ar-max 10000 > _g0c/quality.log 2>&1; rc=$?; cat _g0c/quality.log >> $LOG; echo "  rc=$rc" >> $LOG
[ $rc = 0 ] && grep -q "^VERDICT: PASS" _g0c/quality.log || die "格子の品質 (AR ≤ 10,000)"
echo "== 格子の品質 (既定・参考)" >> $LOG; python3 $TL/check_mesh_quality.py _g0c/nozzle.h5 2>&1 | grep -E "VERDICT|aspect|AR" | head -5 >> $LOG
# 初期値
cp _g0c/nozzle.h5 _g0c/q32_init.h5 || die "初期値の複製"
python3 - _tr/q32_init.h5 _g0c/q32_init.h5 >> $LOG 2>&1 <<'PY' || die "roY* の作成"
import sys, re, h5py
with h5py.File(sys.argv[1], "r") as s, h5py.File(sys.argv[2], "r+") as d:
    ks = sorted(k for k in s["VALUE"] if re.fullmatch(r"roY\d+", k))
    if ks != ["roY0", "roY1"]: raise SystemExit(f"roY* が想定と違う {ks}")
    for k in ks:
        if k not in d["VALUE"]: d["VALUE"].create_dataset(k, shape=(d["VALUE/ro"].shape[0],), dtype=s["VALUE/" + k].dtype)
    print("[g0c] 宛先に roY* を作った:", ks)
PY
python3 $TL/interp_field.py _tr/q32_init.h5 _g0c/q32_init.h5 --dst-run _g0c --forge $F64 > _g0c/interp.log 2>&1; rc=$?
echo "interp_field rc=$rc $(grep -h 'moved \[' _g0c/interp.log | tail -1)" >> $LOG; grep -h "最近傍距離" _g0c/interp.log >> $LOG
[ $rc = 0 ] && grep -q "moved \['ro', 'roUx', 'roUy', 'roUz', 'roe', 'roK', 'roOmega', 'roY0', 'roY1'\]" _g0c/interp.log || die "interp_field (9 量が移っていない)"
python3 - _g0c/q32_init.h5 >> $LOG 2>&1 <<'PY' || die "Q32 の丸め"
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
PY
echo "== 照合 (g0c_check.py: 実際の G0 の入力と)" >> $LOG; python3 g0c_check.py > _g0c/check.log 2>&1; rc=$?; grep -E "bit_equal|not_bit|seg_|VERDICT" _g0c/check.log >> $LOG; echo "  rc=$rc" >> $LOG
[ $rc = 0 ] || die "実際の G0 の入力との照合"
echo "共通の入力 _g0c/q32_init.h5 sha $(sha256sum _g0c/q32_init.h5 | cut -c1-16)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG; touch g0c_prep.done
