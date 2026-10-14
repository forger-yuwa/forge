"""case/70: 帯局所抽出の等価半径の探索の A/B (plan tooling-nozzle-deltastar-band-root-zone §3、読むだけ・forge 0 step)。

仮説の測定: band_local_deficit の r_eff 探索を壁〜y_b の区間に限ると x≈1–3 の δ_r ≈ r_w が消えるか。
腕 orig = 現行 (np.interp を全半径の F で)、腕 zone = 壁〜y_b の区間の F だけで探索。
usage (design/.venv-opt): python ab_root_zone.py CASE_DIR NS_RUN EULER_RUN OUT_DIR"""
import inspect, json, sys
from pathlib import Path
import numpy as np
C = Path(sys.argv[1]).resolve(); sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS
from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
ns, eu, out = C / sys.argv[2], C / sys.argv[3], Path(sys.argv[4])
orig = DS.band_local_deficit
src = inspect.getsource(orig)
old = "            r_eff = float(np.interp(D, F_from_wall, r_desc))\n"
new = ("            kz = int(np.searchsorted(rwN - r_desc, y_b, side='right'))\n"
       "            Fz, rz = F_from_wall[:kz + 1], r_desc[:kz + 1]\n"
       "            _NONMONO.append(bool(np.any(np.diff(Fz) < 0)))\n"
       "            _GLOBAL_TAIL.append(float(F_from_wall[-1]) < D)\n"
       "            r_eff = float(np.interp(D, Fz, rz))\n")
assert src.count(old) == 1
ns_ = dict(DS.__dict__); ns_["_NONMONO"] = []; ns_["_GLOBAL_TAIL"] = []
exec(compile(src.replace(old, new), "patched", "exec"), ns_)
patched = ns_["band_local_deficit"]
res = {}
for arm, fn in (("orig", orig), ("zone", patched)):
    DS.band_local_deficit = fn
    ns_["band_local_deficit"] = fn
    od = out / arm; od.mkdir(parents=True, exist_ok=True)
    ns_["_NONMONO"].clear(); ns_["_GLOBAL_TAIL"].clear()
    summ = extract_and_merge(ns, eu, band_select="edge", out_dir=od)
    d = json.loads((od / "delta_r_equiv_diag.json").read_text())
    x = np.array(d["x"]); raw = np.array(d["delta_r_raw"], float); rwN = np.array(d["r_wall_ns"], float)
    nx = read_delta_r_next(od / "delta_r_next.csv"); xn, dE = np.asarray(nx["x_rt"]), np.asarray(nx["delta_E"])
    m = (x >= 0) & (x <= 3)
    res[arm] = dict(n_eq_rw=int(np.sum(np.abs(raw - rwN) < 1e-3 * rwN)),
                    max_abs_raw_0_3=float(np.nanmax(np.abs(raw[m]))),
                    raw_at={str(s): float(raw[int(np.argmin(abs(x - s)))]) for s in (0, 0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 30, 48)},
                    D_max_abs_deltaE_0_3=float(np.max(np.abs(dE[(xn >= 0) & (xn <= 3)]))),
                    dE_at={str(s): float(np.interp(s, xn, dE)) for s in (0, 1, 2, 3, 10, 30, 48)},
                    massflow_ratio=(summ.get("massflow") or {}).get("ratio"),
                    calls_zone_nonmono=(int(sum(ns_["_NONMONO"])) if arm == "zone" else None),
                    calls_global_tail_below_D=(int(sum(ns_["_GLOBAL_TAIL"])) if arm == "zone" else None),
                    calls=(len(ns_["_NONMONO"]) if arm == "zone" else None))
DS.band_local_deficit = orig
(out / "zone_root_ab.json").write_text(json.dumps(res, indent=1, default=float))
print(json.dumps(res, indent=1, default=float))
