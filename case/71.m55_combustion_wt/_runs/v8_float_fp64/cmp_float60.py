# card plan §6 V8: float (run_0015・0016) と FP64 (run_0012・0013) を step 60,000 の全場で比べる (AWS B、道具は最新の ~/forge-cardmesh/design)
import json, os, sys, hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0, os.path.expanduser("~/forge-cardmesh/design"))
from forge_design.report.nozzle_report import load_field, eta_line
from forge_design.evaluate.stop_rule import extract_delta_e
import h5py
C = Path(os.path.expanduser("~/forge-cardmesh/case/71.m55_combustion_wt")); MD = 5.5; K = 60000
OUT = Path(os.path.expanduser("~/m55_p5/cmp_float60")); OUT.mkdir(parents=True, exist_ok=True)
eref = {"run": str(C / "run_0009_euler"), "res": "res_12000.h5", "sha256": hashlib.sha256((C / "run_0009_euler/res_12000.h5").read_bytes()).hexdigest()}
runs = [("FP64", "run_0012_ns_free"), ("FP64", "run_0013_ns_free260"), ("float", "run_0015_ns_float_a"), ("float", "run_0016_ns_float_b")]
res = {}
for kind, r in runs:
    F = load_field(C / r, f"res_{K}.h5"); S = float(F["info"]["scale_m"])
    M = F["V"]["M"]; eta = F["R"][-1] / F["R"][-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    xq = np.linspace(4.0 / S, 7.95 / S, 2001)
    L03 = eta_line(F, "M", 0.03, xq); L01 = eta_line(F, "M", 0.1, xq)
    d03 = 100 * (L03 / MD - 1); d01 = 100 * (L01 / MD - 1)
    with h5py.File(C / r / f"res_wall_3_{K}.h5", "r") as h:
        Q = np.asarray(h["VALUE/iface_Qf_eff"][:], float) if "VALUE/iface_Qf_eff" in h else None
    de_dir = OUT / f"delta_E_{r}"
    if not (de_dir / "delta_r_next.csv").exists():
        rec = extract_delta_e(C / r, K, eref, de_dir)
    a = np.genfromtxt(de_dir / "delta_r_next.csv", delimiter=",", names=True)
    res[r] = dict(kind=kind, coreM=float(M[-1][core].mean()), d03=d03, d01=d01, Q=Q, de=a,
                  maxabs03=float(np.abs(d03).max()), half03=float(100 * (L03.max() - L03.min()) / (2 * MD)), mean03=float(d03.mean()),
                  maxabs01=float(np.abs(d01).max()), half01=float(100 * (L01.max() - L01.min()) / (2 * MD)), S=S)
    print(f"{kind:5s} {r:22s} 出口コア M {res[r]['coreM']:.7f}  η0.03: 最大偏差 {res[r]['maxabs03']:.4f}% 片振幅 {res[r]['half03']:.4f}% 平均 {res[r]['mean03']:+.4f}%  η0.1: {res[r]['maxabs01']:.4f}/{res[r]['half01']:.4f}%  壁の総熱量 {np.nansum(Q) if Q is not None else None:.6g}", flush=True)
names = res[runs[0][1]]["de"].dtype.names
col = "delta_E"
print("δ_E の表の列:", names, "使う列:", col)
def de_mm(r):
    a = res[r]["de"]; S = res[r]["S"]
    return a["x_rt"] * S, a[col] * S * 1e3   # x [m]、δ_E [mm] (表は r_t 単位)
pairs = [("run_0012_ns_free", "run_0013_ns_free260", "FP64 どうし"), ("run_0015_ns_float_a", "run_0016_ns_float_b", "float どうし"),
         ("run_0012_ns_free", "run_0015_ns_float_a", "FP64−float"), ("run_0013_ns_free260", "run_0016_ns_float_b", "FP64−float")]
print("\n差 (A − B):")
for a, b, lab in pairs:
    A, B = res[a], res[b]
    dM = 100 * (A["coreM"] / B["coreM"] - 1)
    dpos = float(np.abs(A["d03"] - B["d03"]).max())
    dQ = 100 * (np.nansum(A["Q"]) / np.nansum(B["Q"]) - 1) if A["Q"] is not None and B["Q"] is not None else None
    dQloc = float(np.nanmax(np.abs(A["Q"] - B["Q"]) / np.maximum(np.abs(B["Q"]), 1e-30))) if A["Q"] is not None else None
    xa, va = de_mm(a); xb, vb = de_mm(b)
    sel = xa >= 0
    dde = float(np.nanmax(np.abs(va[sel] - np.interp(xa[sel], xb, vb)))) if col else None
    print(f"  {lab:10s} {a} − {b}: 出口コア M {dM:+.5f}%、η0.03 の位置ごとの差の最大 {dpos:.5f}%pt、最大偏差の差 {A['maxabs03'] - B['maxabs03']:+.5f}%pt、片振幅の差 {A['half03'] - B['half03']:+.5f}%pt、δ_E (x≥0) の差の最大 {dde} mm、壁の総熱量 {dQ:+.4f}%、壁の節点ごとの相対差の最大 {dQloc:.3g}")
