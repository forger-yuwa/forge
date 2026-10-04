#!/usr/bin/env python3
"""check_lm_kernel.py の偽合格の回帰試験 (codex result-2 M1)。

実 run の level 2 出力を 1 つ読み (既定 case/57 の run_0014)、メモリ上で壊して判定器に掛ける:
  正常 -> PASS / 1 節点 NaN -> FAIL / 診断場の欠落 -> FAIL / 対角の欠落 -> FAIL / 反復数の欠落 -> FAIL /
  カーネルの評価マスクが広く 0 に落ちた -> FAIL / ソースを 1 節点だけ 1 % ずらした -> FAIL
usage: python3 test_check_lm_kernel.py [RES_H5]
"""
import contextlib, io, os, sys
import h5py, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import check_lm_kernel as ck

f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "../../case/57.transition_flat_plate/run_0014_t3a_lm_unitcheck/res_50.h5")
with h5py.File(f, "r") as h:
    V0 = {k: np.array(h["VALUE"][k], dtype=float) for k in h["VALUE"]}


def run(V):
    with contextlib.redirect_stdout(io.StringIO()):
        W = ck.check_inputs(dict(V))
        return 1 if W is None else ck.judge(W, "(test)", 1e-3, 1e-5)


ok = True
def chk(name, cond):
    global ok; ok &= bool(cond); print("  [%s] %s" % ("OK " if cond else "NG ", name))

chk("正常入力 -> PASS", run(V0) == 0)
i = int(np.argmax(V0["lmPgamma"]))
V = dict(V0); V["lmPgamma"] = V0["lmPgamma"].copy(); V["lmPgamma"][i] = np.nan; chk("lmPgamma の 1 節点が NaN -> FAIL", run(V) == 1)
for k in ("lmFonset", "src_jac_gamma", "src_jac_reth", "lmCorrIter", "gammaEff"):
    V = {n: v for n, v in V0.items() if n != k}; chk(f"{k} が無い -> FAIL", run(V) == 1)
V = dict(V0); V["lmFlength"] = V0["lmFlength"].copy(); V["lmFlength"][: len(V["lmFlength"]) // 2] = 0.0; chk("カーネルの評価マスクが半分 0 -> FAIL", run(V) == 1)
V = dict(V0); V["lmPgamma"] = V0["lmPgamma"].copy(); V["lmPgamma"][i] *= 1.01; chk("P_gamma を 1 節点だけ 1 % ずらす -> FAIL", run(V) == 1)
V = dict(V0); V["gammaEff"] = V0["gammaEff"].copy(); V["gammaEff"][i] += 0.01; chk("gammaEff を 1 節点だけ 0.01 ずらす -> FAIL", run(V) == 1)
V = dict(V0); V["lmCorrIter"] = V0["lmCorrIter"].copy(); V["lmCorrIter"][i] = 100; chk("相関反復が上限 100 に達した節点 -> FAIL", run(V) == 1)
print("VERDICT:", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
