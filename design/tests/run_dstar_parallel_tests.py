"""抽出の高速化 (plan tooling-nozzle-card-local-mesh §5.1 #3q (2)) の試験: 列ごとの前計算 (`band_local_pre`) と列の並列
(`_map_columns`、FORGE_DSTAR_WORKERS) が、逐次・前計算なしと**ビット同一**の結果を返すこと。
(1) `band_local_deficit(..., _pre=band_local_pre(...))` = `_pre` なし (合成の断面、帯の位置を変えて)。
(2) 合成の①の run で edge の抽出: 並列 4 = 逐次 1 (全配列)。
(3) 他のスレッドが動いているときは逐次に落ちる (fork の子のロックの残りを避ける)。
使い方: design/.venv-opt/bin/python tests/run_dstar_parallel_tests.py
"""
import json, os, sys, tempfile, threading
from pathlib import Path
import h5py
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.metrics import deltastar as ds  # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


# (1) 1 断面
r = np.linspace(0.0, 1.0, 61) ** 0.7
qe = lambda rr: 2.0 + 0.05 * rr                                      # noqa: E731
qn = qe(r) * np.where(1.0 - r < 0.08, ((1.0 - r) / 0.08) ** (1 / 7), 1.0)
same = True
for yb in (0.02, 0.05, 0.09, 0.15, 0.3):
    a = ds.band_local_deficit(r, qn, qe, 1.0, y_b_fixed=yb)
    b = ds.band_local_deficit(r, qn, qe, 1.0, y_b_fixed=yb, _pre=ds.band_local_pre(r, qn, qe, 1.0, 4, 0.30, 1))
    same &= (a is None and b is None) or (a is not None and b is not None and all(
        (np.array_equal(np.asarray(a[k]), np.asarray(b[k]), equal_nan=True)) for k in a))
check("(1) _pre あり = なし (5 つの帯の位置、全キー)", same)

# (2) 合成の run
NI, NJ, S = 101, 33, 0.1
XI = np.linspace(-2.0, 10.0, NI)
rw = lambda x: 1.0 + 0.08 * np.clip(x, 0, None) + np.where(x < 0, 0.05 * x * x, 0.0)   # noqa: E731
Q = lambda x: 2.0 / (1.0 + 0.1 * np.clip(x, 0, None))                                 # noqa: E731


def qns(x, rr):
    dl = 0.03 + 0.006 * np.clip(x, 0, None); s = np.clip((rw(x) - rr) / dl, 0, None)
    return Q(x) * np.where(s < 1, s ** (1 / 7), 1.0)


def write_run(rd, field):
    rd.mkdir(parents=True)
    t = 1.0 - np.arange(NJ) / (NJ - 1)
    X = np.repeat(XI[:, None], NJ, axis=1); R = rw(XI)[:, None] * (1 - t[None, :] ** 1.6)
    co = np.zeros((NI * NJ, 3)); co[:, 0] = X.ravel() * S; co[:, 1] = R.ravel() * S
    with h5py.File(rd / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = co.ravel()
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": S, "x_E": 9.0}))
    with h5py.File(rd / "res_100.h5", "w") as f:
        f["/VALUE/ro"] = np.ones(NI * NJ); f["/VALUE/Ux"] = np.asarray(field(X, R), float).ravel()
    return rd


T = Path(tempfile.mkdtemp(prefix="dstar_par_"))
ns = write_run(T / "ns", qns); eu = write_run(T / "eu", lambda X, R: Q(X))
kw = dict(band_select="edge", remap_constant_x=True, remap_edge="missing")
os.environ["FORGE_DSTAR_WORKERS"] = "1"; o1 = ds.deltastar_from_core_matched_euler(ns, eu, **kw)
os.environ["FORGE_DSTAR_WORKERS"] = "4"; o4 = ds.deltastar_from_core_matched_euler(ns, eu, **kw)
keys = [k for k in o1 if isinstance(o1[k], np.ndarray)]
diff = [k for k in keys if not np.array_equal(o1[k], o4[k], equal_nan=(o1[k].dtype.kind == "f"))]
check("(2) edge の抽出: 並列 4 = 逐次 1 (全配列)", not diff and len(keys) > 5, f"配列 {len(keys)} 個、違い {diff}")

# (3) 他のスレッドが動いていれば逐次
ev = threading.Event(); th = threading.Thread(target=ev.wait); th.start()
try:
    pid_seen = ds._map_columns(lambda i: os.getpid(), range(64))
finally:
    ev.set(); th.join()
check("(3) 他のスレッドがあるときは逐次 (子を作らない)", set(pid_seen) == {os.getpid()})
os.environ.pop("FORGE_DSTAR_WORKERS", None)
print("\nALL PASS" if not FAIL else f"\nFAIL {len(FAIL)}: {FAIL}")
sys.exit(1 if FAIL else 0)
