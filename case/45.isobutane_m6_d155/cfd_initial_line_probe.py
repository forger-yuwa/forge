"""CFD ピンの実現性: Euler 解の中でスロート壁点 (x_w, 1) から軸へ C⁻ (dr/dx = tan(θ−μ)) をたどり、その線上の M, θ を Hall 解と比べる。
CFD データの C⁻ 適合残差 (moc_kernel と同じ台形式・設計側のガス) も出す。形状・既存 run の後処理のみ (CFD 0 step)。
plan: (起票予定) CFD ピン。usage: design/.venv-opt/bin/python cfd_initial_line_probe.py RUN [RES] → _band_ab/cfd_initial_line_probe_<run>.json
"""
import json, sys
from pathlib import Path
import numpy as np, h5py
from scipy.interpolate import LinearNDInterpolator
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from forge_design.geometry.moc_kernel import pm_nu, _sin_over_r_vec  # noqa: E402

run = C / sys.argv[1]; resn = sys.argv[2] if len(sys.argv) > 2 else sorted(run.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1].name
info = json.loads((run / "prepare_info.json").read_text()); S = float(info["scale_m"])
with h5py.File(run / "nozzle.h5") as f:
    nc = f["/MESH/COORD"][:].reshape(-1, 3)
with h5py.File(run / resn) as f:
    Ux, Uy, son = (f["/VALUE/" + k][:] for k in ("Ux", "Uy", "sonic"))
x, r = nc[:, 0] / S, nc[:, 1] / S
m = (x > -1.0) & (x < 2.0)
M = np.hypot(Ux, Uy) / son; th = np.arctan2(Uy, Ux)
itp = LinearNDInterpolator(np.c_[x[m], r[m]], np.c_[M[m], th[m]])
# 壁: 各 x station の r 最大 (構造格子) → スロート壁点 = r 最小の station
ni = int(info["mesh"]["ni"]); X = x.reshape(ni, -1); Rr = r.reshape(ni, -1)
xw, rw = X[:, -1], Rr[:, -1]
p0 = np.array([0.0, float(np.interp(0.0, xw, rw))])   # 設計スロート (x=0) の壁点。Hall の特性線と同じ始点
# C⁻ を壁の少し内側から軸へ (RK2, ds 2e-4)。壁直上は補間が壁外に出るので 1e-4 内側から
pts = [p0 - np.r_[0.0, 1e-4]]; ds = 2e-4
while pts[-1][1] > 1e-4 and len(pts) < 200000:
    q = pts[-1]
    def ang(qq):
        Mv, tv = itp(qq[0], qq[1])[0] if np.ndim(itp(qq[0], qq[1])) > 1 else itp(qq[0], qq[1])
        return tv - np.arcsin(1.0 / max(Mv, 1.0 + 1e-9))
    a1 = ang(q); qm = q + 0.5 * ds * np.r_[np.cos(a1), np.sin(a1)]; a2 = ang(qm)
    if not np.isfinite(a2):
        break
    pts.append(q + ds * np.r_[np.cos(a2), np.sin(a2)])
P = np.array(pts); MV = np.array([itp(*p).ravel() for p in P])
ok = np.all(np.isfinite(MV), axis=1) & (P[:, 1] > 0); P, MV = P[ok], MV[ok]
# 41 点に r 等間隔で再標本化 (Hall と同じ規約: 軸→壁)
o = np.argsort(P[:, 1]); rq = np.linspace(0.0, P[0, 1], 41)
# 軸 (r=0) は最寄り 20 点の偶関数 (M) / 奇関数 (θ) 当てはめで外挿
k20 = o[:20]; cM = np.polyfit(P[k20, 1] ** 2, MV[k20, 0], 1); cx = np.polyfit(P[k20, 1], P[k20, 0], 2)
Mq = np.where(rq < P[o[0], 1], np.polyval(cM, rq ** 2), np.interp(rq, P[o, 1], MV[o, 0]))
tq = np.where(rq < P[o[0], 1], MV[o[0], 1] * rq / P[o[0], 1], np.interp(rq, P[o, 1], MV[o, 1]))
xq = np.where(rq < P[o[0], 1], np.polyval(cx, rq), np.interp(rq, P[o, 1], P[o, 0]))
p = load_problem(C / "problem_d155_euler_c2final_n2400.yaml"); d = design_chain(p); g = d["gamma_hall"]
ht = HallThroat(R=d["R"], gamma=g); xh, rh, Mh, thh = ht.throat_characteristic(n=41)
# CFD データの C⁻ 適合残差 (Hall γ の CPG 式)
nu = np.array([float(pm_nu(float(v), g)) for v in Mq]); mu = np.arcsin(1 / Mq); E = []
for k in range(1, 41):
    a = 0.5 * (tq[k - 1] + tq[k] - mu[k - 1] - mu[k])
    fB = 0.5 * (np.sin(mu[k]) * float(_sin_over_r_vec(np.r_[rq[k]], np.r_[tq[k]], np.r_[rq[k - 1]], np.r_[tq[k - 1]])[0])
                + np.sin(mu[k - 1]) * float(_sin_over_r_vec(np.r_[rq[k - 1]], np.r_[tq[k - 1]], np.r_[rq[k]], np.r_[tq[k]])[0]))
    E.append((tq[k - 1] + nu[k - 1]) - (tq[k] + nu[k]) - fB / np.cos(a) * (xq[k - 1] - xq[k]))
E = np.degrees(np.array(E))
out = dict(run=sys.argv[1], res=resn, throat_wall=[float(p0[0]), float(p0[1])], cfd_axis_x=float(xq[0]), trace_end_r=float(P[o[0], 1]), hall_axis_x=float(xh[0]),
           M_wall=[float(Mq[-1]), float(Mh[-1])], M_axis=[float(Mq[0]), float(Mh[0])],
           dM_max=float(np.abs(Mq - np.interp(rq, rh, Mh)).max()), dtheta_max_deg=float(np.degrees(np.abs(tq - np.interp(rq, rh, thh)).max())),
           cfd_cminus_resid_sum_deg=float(E.sum()), cfd_cminus_resid_max_deg=float(np.abs(E).max()),
           line=dict(r=rq.tolist(), x=xq.tolist(), M=Mq.tolist(), theta=tq.tolist()))
(C / f"_band_ab/cfd_initial_line_probe_{sys.argv[1][:8]}.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "line"}, indent=1))
