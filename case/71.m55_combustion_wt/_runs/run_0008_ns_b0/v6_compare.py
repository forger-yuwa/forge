#!/usr/bin/env python3
"""plan tooling-nozzle-card-local-mesh §6 V6 (事前登録 2026-10-15): point (run_0007_ns res_200000) と B0 (run_0008_ns_b0 res_90000) の止まった場を
同じものさしで比べる。格子は同じ (MESH/COORD の sha256 一致を確認済み)。θ_r は帯が run ごとに違うので、両方の帯 (A・B) で両方の場を測る。
δ_E は同じ Euler の参照・同じ固定の抽出の条件 (handoff_extract_settings) で両方の場から run の外に抽出する。
使い方 (AWS A): PYTHONPATH=~/forge-cardmesh/design python v6_compare.py <case71> <b0 dir (res_90000.h5・stop_rule.json)> <out dir>
"""
import json, os, sys
from pathlib import Path
import numpy as np
from forge_design.evaluate import stop_rule as SR
from forge_design.feedback import deltastar_loop as DL

C, B, O = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); O.mkdir(parents=True, exist_ok=True)
RA_ = C / "run_0007_ns"
geo = SR.RunGeometry(RA_)
VA = SR.read_values(RA_ / "res_200000.h5"); VB = SR.read_values(B / "res_90000.h5")
bands = {"bandA(run_0007)": json.loads((RA_ / "stop_rule.json").read_text())["band"],
         "bandB(run_0008)": json.loads((B / "stop_rule.json").read_text())["band"]}
xf = json.loads((RA_ / "stop_rule.json").read_text())["settings"].get("x_fracs") or SR.X_FRACS
out = {"theta": {}}
for bn, band in bands.items():
    ra = SR.snapshot(geo, 200000, band, x_fracs=xf, values=VA); rb = SR.snapshot(geo, 90000, band, x_fracs=xf, values=VB)
    out["theta"][bn] = {k: {"point": ra[k], "b0": rb[k], "rel_diff_pct": 100.0 * (rb[k] / ra[k] - 1.0)} for k in ra if k.startswith("theta_r_")}
    out["exit_core_M"] = {"point": ra["exit_core_M"], "b0": rb["exit_core_M"], "rel_diff_pct": 100.0 * (rb["exit_core_M"] / ra["exit_core_M"] - 1.0)}
print(json.dumps(out, indent=1, ensure_ascii=False), flush=True)
# δ_E: B0 の場は run_0007 の格子・prepare_info を使った仮の run dir で抽出する (格子が同じなので同じ)
pb = O / "b0run"; pb.mkdir(exist_ok=True)
for f in ("nozzle.h5", "prepare_info.json"):
    if not (pb / f).exists():
        os.symlink(RA_ / f, pb / f)
if not (pb / "res_90000.h5").exists():
    os.symlink(B / "res_90000.h5", pb / "res_90000.h5")
info = json.loads((RA_ / "prepare_info.json").read_text())
eref = SR.euler_ref_spec(info["stop_rule"]["euler_ref"]) if isinstance(info.get("stop_rule"), dict) else None
st = DL.handoff_extract_settings()
S = float(info["scale_m"]); mm = S * 1e3
dE = {}
for tag, run, res in (("point", RA_, "res_200000.h5"), ("b0", pb, "res_90000.h5")):
    od = O / f"extract_{tag}"
    DL.extract_and_merge(run, eref, out_dir=od, ns_res=res, **st)
    d = DL.read_delta_r_next(od / "delta_r_next.csv"); dE[tag] = (np.asarray(d["x_rt"]), np.asarray(d["delta_E"]))
    print(tag, "抽出 済み", flush=True)
x = dE["point"][0]; a = dE["point"][1]; b = np.interp(x, dE["b0"][0], dE["b0"][1])
m = x >= 0.0
i = int(np.argmax(np.abs(b - a)[m]))
out["delta_E"] = {"throat_mm": {"point": float(np.interp(0.0, x, a) * mm), "b0": float(np.interp(0.0, x, b) * mm)},
                  "exit_mm": {"point": float(a[-1] * mm), "b0": float(b[-1] * mm)},
                  "max_abs_diff_mm_x_ge_0": float(np.abs(b - a)[m][i] * mm), "x_m_at": float(x[m][i] * S)}
(O / "v6_compare.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out["delta_E"], indent=1, ensure_ascii=False))
