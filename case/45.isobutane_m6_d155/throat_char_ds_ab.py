"""スロート特性線の追跡刻み ds の A/B (形状のみ, CFD 0 step)。plan verification-m6-axis-wave-mesh-su2 §5.1 #13 ① (事前登録 commit 4bd81610)。
A: ds 2e-4 (既定) / B: ds 2e-5。n_start 41・n_axis_inv 2400・ガス・拘束・当てはめは固定。
見る量: (1) 初期線の壁側 C⁻ 適合残差 (Hall γ の CPG 式, 台形), (2) MOC 壁の最初の区間の割線−平均角度差, (3) 同時当てはめ B (λ=1e-9) の全点誤差。
判定 (事前登録): (1)(2) がともに 80 % 以上減る → 追跡誤差が主因。減らない → 追跡刻み主因の仮説を棄却。
usage: design/.venv-opt/bin/python throat_char_ds_ab.py → _band_ab/throat_char_ds_ab.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.geometry import transonic
from forge_design.geometry.moc_kernel import pm_nu
from forge_design.evaluate.runner_axismach import design_chain, load_problem
import moc_wall_fit_ab as F  # joint_fit を流用 (import で moc_wall_fit_ab の測定も 1 回走る)

C = Path(__file__).resolve().parent
orig = transonic.HallThroat.throat_characteristic


def cminus_resid(x, r, M, th, g):
    nu = np.array([pm_nu(float(m), g) for m in M]); mu = np.arcsin(1.0 / M)
    out = []
    for i in range(len(x) - 1):
        if min(r[i], r[i + 1]) < 0.05:
            out.append(np.nan); continue
        ang = 0.5 * (th[i] + th[i + 1]) - 0.5 * (mu[i] + mu[i + 1])
        f = 0.5 * (np.sin(mu[i]) * np.sin(th[i]) / r[i] + np.sin(mu[i + 1]) * np.sin(th[i + 1]) / r[i + 1])
        out.append((th[i + 1] + nu[i + 1]) - (th[i] + nu[i]) - f / np.cos(ang) * (x[i + 1] - x[i]))
    return np.degrees(np.array(out))


res = {}
for name, ds in (("A_ds2e-4", 2e-4), ("B_ds2e-5", 2e-5)):
    transonic.HallThroat.throat_characteristic = (lambda self, n=61, ds=ds, max_steps=4000000: orig(self, n=n, ds=ds, max_steps=max_steps))
    p = load_problem(C / "problem_d155_ns_c2final.yaml"); p.geometry["n_axis_inv"] = 2400; p.geometry["n_start"] = 41
    d = design_chain(p); tb = d["wall_inv"]; g = float(d["gamma_hall"])
    ht = transonic.HallThroat(R=d["R"], gamma=g)
    xl, rl, Ml, tl = ht.throat_characteristic(n=41)
    cr = cminus_resid(np.asarray(xl), np.asarray(rl), np.asarray(Ml), np.asarray(tl), g)
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    seg = np.degrees(np.arctan(np.diff(r[:4]) / np.diff(x[:4])) - np.arctan(np.tan(0.5 * (th[1:4] + th[:3]))))
    s, nc = F.joint_fit(tb, d["R"], 1e-9)
    dr = s(x) - r; dth = np.degrees(np.arctan(s(x, 1)) - th)
    res[name] = dict(ds=ds, x0_axis=float(xl[0]), cminus_wall_side_deg=float(cr[-1]), cminus_max_abs_deg=float(np.nanmax(np.abs(cr))),
                     first_seg_deg=float(seg[0]), next_segs_deg=[float(v) for v in seg[1:]],
                     fitB_max_dr=float(np.abs(dr).max()), fitB_x_max_dr=float(x[np.argmax(np.abs(dr))]),
                     fitB_max_dth_deg=float(np.abs(dth).max()), fitB_x_max_dth=float(x[np.argmax(np.abs(dth))]),
                     x_wall_end=float(x[-1]), n_pts=len(tb))
transonic.HallThroat.throat_characteristic = orig
a, b = res["A_ds2e-4"], res["B_ds2e-5"]
red = {k: 1 - abs(b[k]) / abs(a[k]) for k in ("cminus_wall_side_deg", "first_seg_deg")}
res["reduction"] = red
res["verdict"] = ("追跡誤差が主因 (両方 ≥80 % 減)" if min(red.values()) >= 0.8 else "追跡刻み主因の仮説を棄却 (少なくとも一方が 80 % 未満)")
(C / "_band_ab/throat_char_ds_ab.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print(json.dumps(res, indent=1, ensure_ascii=False))
