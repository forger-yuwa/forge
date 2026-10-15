"""流量の帳簿を一定 x の断面で積分する (plan tooling-nozzle-card-local-mesh §5.1 #3i) の試験。

(1) 列が一定 x の格子では、一定 x にそろえても流量の帳簿はビット同一 (既存の①の run の値は変わらない)。
(2) 列が傾いた格子 (壁法線の近壁層) では、抽出の帳簿 (`deltastar_from_core_matched_euler` の `massflow`) が抽出と同じ断面
    (remap_constant_x) で測る: True なら一定 x の値、False (旧い呼び方) なら今までどおり列のままの値。
(3) `runner_axismach.collect` が一定 x の帳簿を呼ぶ。
使い方: design/.venv-opt/bin/python tests/run_massflow_constant_x_tests.py
"""
import inspect
import json
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.metrics import deltastar as ds  # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


NI, NJ, P, S, SIN = 121, 41, 1.6, 0.1, 0.3
XI = np.linspace(-3.0, 10.0, NI)


def rw(x):
    return 1.0 + 0.08 * np.clip(x, 0.0, None) + np.where(x < 0, 0.05 * x * x, 0.0)


def Q(x):
    return 2.0 / (1.0 + 0.1 * np.clip(x, 0.0, None))


def q_ns(x, r):
    dl = 0.03 + 0.006 * np.clip(x, 0.0, None)
    s = np.clip((rw(x) - r) / dl, 0.0, None)
    return Q(x) * np.where(s < 1.0, s ** (1.0 / 7.0), 1.0)


def write_run(rd: Path, slant: float, field):
    rd.mkdir(parents=True, exist_ok=True)
    t = 1.0 - np.arange(NJ) / (NJ - 1)
    X = np.empty((NI, NJ)); R = np.empty((NI, NJ))
    for i, x0 in enumerate(XI):
        y = rw(x0) * t ** P
        R[i] = rw(x0) - y
        X[i] = x0 + slant * y * SIN * (1.0 - t)
    co = np.zeros((NI * NJ, 3)); co[:, 0] = X.ravel() * S; co[:, 1] = R.ravel() * S
    with h5py.File(rd / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = co.ravel()
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": S, "x_E": 9.0}))
    with h5py.File(rd / "res_100.h5", "w") as f:
        f["/VALUE/ro"] = np.ones(NI * NJ); f["/VALUE/Ux"] = np.asarray(field(X, R), dtype=float).ravel()
    return rd


TMP = Path(tempfile.mkdtemp(prefix="mdot_cx_test_"))
eu = write_run(TMP / "eu", 0.0, lambda X, R: Q(X))
ns0 = write_run(TMP / "ns_straight", 0.0, lambda X, R: q_ns(X, R))
ns3 = write_run(TMP / "ns_slant", 3.0, lambda X, R: q_ns(X, R))
KEYS = ("mdot_ns", "mdot_ns_spread", "mdot_euler", "mdot_euler_spread", "mdot_ratio", "cd_rel", "delta_r_throat_eff")

# (1) 列が一定 x の格子はビット同一
a = ds.massflow_ratio(ns0, eu); b = ds.massflow_ratio(ns0, eu, remap_constant_x=True)
check("(1) 列が一定 x の格子では一定 x にそろえても帳簿はビット同一", all(a[k] == b[k] for k in KEYS),
      " ".join(f"{k} {a[k]!r}/{b[k]!r}" for k in ("mdot_ns", "mdot_ratio")))

# (2) 列が傾いた格子: 抽出の帳簿は抽出と同じ断面
col = ds.massflow_ratio(ns3, eu); cx = ds.massflow_ratio(ns3, eu, remap_constant_x=True)
check("(2) 傾いた格子では列のままと一定 x の帳簿が違う (試験の前提)", abs(col["mdot_ratio"] / cx["mdot_ratio"] - 1.0) > 1e-3,
      f"列 {col['mdot_ratio']:.6f} / 一定 x {cx['mdot_ratio']:.6f}; 散らばり 列 {col['mdot_ns_spread']:.2e} 一定 x {cx['mdot_ns_spread']:.2e}")
# 真値 = 一定 x の断面の厳密な比 (合成の場は x で流量が保存しないので散らばりでなく比で見る)
xs_all = XI[(XI > -2.5) & (XI < 9.0)]
rt_ = []
for x in xs_all:
    r = np.linspace(0.0, rw(x), 400001)
    rt_.append(np.trapezoid(q_ns(x, r) * r, r) / np.trapezoid(Q(x) * r, r))
m_true = float(np.median(rt_))
check("(2) 一定 x の帳簿の比は列のままより真値に近い", abs(cx["mdot_ratio"] / m_true - 1) < 0.2 * abs(col["mdot_ratio"] / m_true - 1),
      f"真値 {m_true:.6f}、一定 x {cx['mdot_ratio'] / m_true - 1:+.2e}、列 {col['mdot_ratio'] / m_true - 1:+.2e}")
o_t = ds.deltastar_from_core_matched_euler(ns3, eu, remap_constant_x=True, remap_edge="missing")
o_f = ds.deltastar_from_core_matched_euler(ns3, eu)
check("(2) 抽出 (remap_constant_x=True) の帳簿 = 一定 x の帳簿", all(o_t["massflow"][k] == cx[k] for k in KEYS))
check("(2) 抽出 (旧い呼び方、remap なし) の帳簿 = 列のままの帳簿 (振る舞いを変えない)", all(o_f["massflow"][k] == col[k] for k in KEYS))

# (3) collect は一定 x の帳簿を呼ぶ
src = inspect.getsource(RA.collect)
check("(3) runner_axismach.collect は massflow_ratio(..., remap_constant_x=True) を呼ぶ",
      "massflow_ratio(run_dir, info[\"euler_ref\"], remap_constant_x=True)" in src)

print("\nALL PASS" if not FAIL else f"\nFAIL {len(FAIL)}: {FAIL}")
sys.exit(1 if FAIL else 0)
