"""二相拡散 既定化 plan §5.1 #4s: G1(ii) の指数スケール A/B。
export: tp_faces.h5 の on_in/* を binary に書く。judge: A (k=0) と B (k=64) の Jl を事前登録の判定で評価する。
  python3 scale_ab.py export tp_faces.h5 faces.bin
  python3 scale_ab.py judge tp_faces.h5 A.bin B.bin 64 faces.bin.idx.npy
"""
import sys, numpy as np, h5py

def export(h5p, out):
    f = h5py.File(h5p); g = lambda k: np.array(f[k][...])
    skip = g("face/skip"); ev = np.where(skip == 0)[0]
    n = g("on_in/rY0").shape[1]; iw = int(f.attrs.get("iw", n - 1)) if "iw" in f.attrs else n - 1
    cols = [g("on_in/rho0"), g("on_in/rho1"), g("on_in/rg0"), g("on_in/rg1"), g("on_in/f"), g("on_in/geo"), g("on_in/geo_abs"), g("on_in/ct"), g("on_in/L")]
    for key in ("rY0", "rY1", "D", "h"):
        a = g("on_in/" + key)
        cols += [a[:, s] for s in range(n)]
    for key in ("rQ0", "rQ1"):
        a = g("on_in/" + key)
        cols += [a[:, m] for m in range(3)]
    X = np.stack([c[ev].astype(np.float32) for c in cols], axis=1)
    with open(out, "wb") as fo:
        np.array([len(ev), n, iw, X.shape[1]], dtype=np.int32).tofile(fo); X.tofile(fo)
    np.save(out + ".idx.npy", ev)
    print("exported", len(ev), "faces, n", n, "iw", iw, "ncol", X.shape[1])

def read(p, N):
    with open(p, "rb") as fi:
        Jf = np.fromfile(fi, dtype=np.float32, count=N); Jd = np.fromfile(fi, dtype=np.float64, count=N)
    return Jf.astype(np.float64), Jd

def judge(h5p, pa, pb, k, idxp):
    f = h5py.File(h5p); g = lambda key: np.array(f[key][...])
    idx = np.load(idxp)
    N = len(idx); eps = 2.0 ** -23; scale = 2.0 ** k
    JfA, JdA = read(pa, N); JfB, JdB = read(pb, N)
    stored = g("on_f/Jl")[idx].astype(np.float64)
    rho0 = g("on_in/rho0")[idx].astype(float); rho1 = g("on_in/rho1")[idx].astype(float)
    rg0 = g("on_in/rg0")[idx].astype(float); rg1 = g("on_in/rg1")[idx].astype(float)
    ctg = (g("on_in/ct")[idx].astype(float) * g("on_in/geo")[idx].astype(float))
    A_A = np.abs(ctg) * (np.abs(rg0 / rho0) + np.abs(rg1 / rho1)); A_B = A_A * scale
    bitA = np.array_equal(JfA.astype(np.float32).view(np.int32), stored.astype(np.float32).view(np.int32))
    print("precondition: A Jl_f bit-identical to the stored diagnostic value:", bitA, "(mismatches", int(np.sum(JfA != stored)), ")")
    with np.errstate(divide="ignore", invalid="ignore"):
        rA = np.abs(JfA - JdA) / (8 * eps * A_A); rB = np.abs(JfB - JdB) / (8 * eps * A_B)
    act = A_A > 0
    failA = act & (rA > 1)
    tiny = np.finfo(np.float32).tiny
    normal = act & (np.abs(JdA) >= tiny)
    ctrl = np.where(normal)[0][np.argsort(-np.abs(JdA[normal]))[: int(failA.sum())]]
    print(f"A: active faces {int(act.sum())}, over old scale {int(failA.sum())}, max ratio {np.nanmax(rA[act]):.3e}")
    print(f"B (2^{k}): over old scale on the A-failing faces {int(np.sum(rB[failA] > 1))} / {int(failA.sum())}, max ratio {np.nanmax(rB[failA]) if failA.any() else 0:.3e}; "
          f"control {len(ctrl)} faces max ratio {np.nanmax(rB[ctrl]) if len(ctrl) else 0:.3e} over {int(np.sum(rB[ctrl] > 1))}; all active faces over {int(np.sum(rB[act] > 1))} max ratio {np.nanmax(rB[act]):.3e}")
    nf = int(np.sum(~np.isfinite(JfB[act]))) + int(np.sum(~np.isfinite(JdB[act])))
    print("B non-finite:", nf, "| B results normal:", int(np.sum(np.abs(JfB[failA]) >= tiny)), "/", int(failA.sum()))
    ok = bitA and failA.any() and np.sum(rB[act] > 1) == 0 and nf == 0
    print("VERDICT #4s:", "SUPPORT (A fails, B all within the old scale) -> derive underflow-aware scale" if ok else "REJECT or precondition not met -> stop relaxing the scale")

if __name__ == "__main__":
    if sys.argv[1] == "export": export(sys.argv[2], sys.argv[3])
    else: judge(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6])
