"""plan tooling-nozzle-core-grid §4.90 P-g の線形性の確認: Q0 + δ★/2 の状態を forge の入力 h5 に書く。
Q0 = 凍結の状態 (run_0693 の res_20000)、Q1★ = 800 sweep の run の res_1 (更新後の保存量)。保存量 9 量 (ro・roUx・roUy・roUz・roe・roK・roOmega・roY0・roY1)
それぞれに Q0 + ½ (Q1★ − Q0) を書く。書く前に、入力 h5 の値が Q0 とビット一致であること (restart_field の直後) を確かめる。
usage: python3 cg_pg_half.py <Q0 の h5> <Q1★ の h5> <書き込む入力 h5>"""
import sys
import h5py, numpy as np
q0p, q1p, dstp = sys.argv[1:4]
NAMES = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
with h5py.File(q0p, "r") as a, h5py.File(q1p, "r") as b, h5py.File(dstp, "r+") as d:
    if not np.array_equal(np.asarray(a["MESH/COORD"][:]), np.asarray(d["MESH/COORD"][:])) or not np.array_equal(np.asarray(b["MESH/COORD"][:]), np.asarray(d["MESH/COORD"][:])):
        sys.exit("MESH/COORD が一致しない")
    rel = {}
    for k in NAMES:
        q0 = np.asarray(a["VALUE/" + k][:], float); q1 = np.asarray(b["VALUE/" + k][:], float); cur = np.asarray(d["VALUE/" + k][:], float)
        if not np.array_equal(cur, q0): sys.exit(f"{k}: 入力 h5 の値が Q0 と一致しない (restart の直後でない)")
        half = q0 + 0.5 * (q1 - q0)
        d["VALUE/" + k][...] = half.astype(d["VALUE/" + k].dtype, copy=False)
        den = np.max(np.abs(q0)); rel[k] = float(np.max(np.abs(q1 - q0)) / den) if den > 0 else 0.0
print("Q0 + δ★/2 を書いた。max|δ★|/max|Q0|: " + " ".join(f"{k} {v:.2e}" for k, v in rel.items()))
