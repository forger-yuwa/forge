#!/usr/bin/env python3
"""plan tooling-nozzle-card-local-mesh (codex diagnose 2026-10-15 m55-settle-and-monitor の判別 A/B): run_0010_ns_b0 の res_86000 を固定し、
δ_E の抽出の帯の位置だけを変える (A = 選ばれた y_b、B = 1.25 y_b に固定、梯子は探し直さない)。forge は 0 step。
固定: Euler の参照 run_0009_euler/res_12000 (sha256 は prepare_info の記録)、一定 x の補間 (missing)、refine、平滑化 λ 1・ノット 2 r_t・重み。
事前の許容: x 0.5〜7.9 m で max|δ_E,B − δ_E,A| ≤ 0.05 mm、かつ (δ_in − δ_E) の符号が x 4〜6 m で負・7〜7.9 m で正のまま。
A は保存済みの δ_E (extract/delta_r_next.csv) を再現すること。
使い方: PYTHONPATH=<design> python band_ab.py <run> <保存済み delta_r_next.csv> <out.json> [<Euler の参照の run の場所>]
"""
import json, sys
from pathlib import Path
import numpy as np
from forge_design.evaluate import stop_rule as SR
from forge_design.feedback import deltastar_loop as DL
from forge_design.metrics.deltastar import deltastar_from_core_matched_euler, smooth_delta_quintic

run, saved, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
info = json.loads((run / "prepare_info.json").read_text())
er = dict(info["stop_rule"]["euler_ref"])
if len(sys.argv) > 4:                                   # Euler の参照の run の場所だけ差し替える (sha256 は記録のまま照合)
    er["run"] = sys.argv[4]
eref = SR.euler_ref_spec(er)
st = DL.handoff_extract_settings()
kw = {k: v for k, v in st.items() if k in ("band_select", "remap_constant_x", "remap_edge")}
o = deltastar_from_core_matched_euler(run, eref, ns_res="res_86000.h5", yb_sens_factors=[1.25], **kw)
x = o["x"]; S = float(info["scale_m"]); mm = S * 1e3
w = (o["hard_ok"] & np.isfinite(o["delta_r_raw"])).astype(float)
fA, _ = smooth_delta_quintic(x, o["delta_r_raw"], weights=w, knot_spacing=st["knot_spacing"], lam=st["smooth_lam"], positive=False)
rawB = o["delta_r_sens_yb"][:, 0]
wB = (o["hard_ok"] & np.isfinite(rawB)).astype(float)
fB, _ = smooth_delta_quintic(x, rawB, weights=wB, knot_spacing=st["knot_spacing"], lam=st["smooth_lam"], positive=False)
A, B = fA(x), fB(x)
nx = DL.read_delta_r_next(saved)
rep = {"A_reproduces_saved_max_abs_mm": float(np.max(np.abs(A - np.asarray(nx["delta_E"]))) * mm),
       "n_nonfinite_B": int(np.sum(~np.isfinite(rawB))), "n_hard_fail": int(np.sum(~o["hard_ok"]))}
xm = x * S; m = (xm >= 0.5) & (xm <= 7.9)
d = (B - A) * mm; i = int(np.argmax(np.abs(d[m])))
din = np.asarray(nx["delta_in_prev"])
eA = (din - A) * mm; eB = (din - B) * mm
neg = (xm >= 4.0) & (xm <= 6.0); pos = (xm >= 7.0) & (xm <= 7.9)
rep.update({"max_abs_dB_minus_A_mm": float(np.abs(d[m][i])), "x_m_at": float(xm[m][i]),
            "sign_A": {"neg_4_6": bool(np.all(eA[neg] < 0)), "pos_7_79": bool(np.all(eA[pos] > 0))},
            "sign_B": {"neg_4_6": bool(np.all(eB[neg] < 0)), "pos_7_79": bool(np.all(eB[pos] > 0))},
            "err_at_mm": {f"{v}": {"A": float(np.interp(v, xm, eA)), "B": float(np.interp(v, xm, eB))} for v in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.9)}})
ok = (rep["max_abs_dB_minus_A_mm"] <= 0.05 and rep["sign_B"]["neg_4_6"] and rep["sign_B"]["pos_7_79"] and rep["n_nonfinite_B"] == 0)
rep["result"] = ("比較不成立 (A が保存を再現しない)" if rep["A_reproduces_saved_max_abs_mm"] > 1e-9 else
                 ("許容内 (帯の位置依存を主因から外す; 形の残差の仮説を残す)" if ok else "許容外 (今の δ_E を形の誤差の根拠にしない; 4 係数化を止める)"))
out.write_text(json.dumps(rep, indent=1, ensure_ascii=False)); print(json.dumps(rep, indent=1, ensure_ascii=False))
