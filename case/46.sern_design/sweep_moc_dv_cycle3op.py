"""cycle3op の dv 箱で MOC の L_ramp / C_T / C_L / C_M 分布を見る (opt.cm_min の根拠, plan §5.1-1b)。CFD 不要。
実行: .venv-chem/bin/python case/46.sern_design/sweep_moc_dv_cycle3op.py > case/46.sern_design/sweep_moc_dv_cycle3op.log"""
import sys, itertools, numpy as np
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec, wall_forces
M_IN, P_IN, GAM, P_EXT = 1.6745, 101027.0, 1.1828, 2851.0
ratio, XR, LMAX = P_EXT / P_IN, -20.0, 20.0
grid = dict(M_c=[2.4, 2.7, 3.0, 3.3], f=[0.35, 0.475, 0.6], tr=[10.0, 16.0, 22.0], tc=[2.0, 5.0, 8.0], Lc=[0.6, 1.55, 2.5])
rows = []
for M_c, f, tr, tc, Lc in itertools.product(*grid.values()):
    spec = SernKernelSpec(M_in=M_IN, theta_r0=np.deg2rad(tr), theta_c0=np.deg2rad(tc), L_cowl=Lc,
                          gamma=GAM, p_ext_over_p_in=ratio, x_max=14.0, nj=151, dx=4e-3)
    try:
        k = PlanarMOC(spec).march(stop_at=(f, M_c)); d = k.design_ramp(M_c=M_c, f=f, ds=2e-3)
        if d.info["warnings"] or d.L_ramp > LMAX: continue
        fr = wall_forces(d, M_IN, GAM, pa_over_pin=ratio, x_ref=XR, y_ref=0.0)
        rows.append((M_c, f, tr, tc, Lc, d.L_ramp, fr["C_T"], fr["C_L"], fr["C_M"]))
        print(f"ok M_c {M_c} f {f} tr {tr} tc {tc} Lc {Lc} L {d.L_ramp:.2f} C_T {fr['C_T']:.4f} C_M {fr['C_M']:+.3f}", flush=True)
    except Exception as e:
        print(f"-- M_c {M_c} f {f} tr {tr} tc {tc} Lc {Lc}: {str(e)[:40]}", flush=True)
a = np.array([r[5:] for r in rows]); n_all = int(np.prod([len(v) for v in grid.values()]))
print(f"\n成立 {len(rows)} / {n_all} 点")
for i, nm in enumerate(("L_ramp", "C_T", "C_L", "C_M")):
    q = np.percentile(a[:, i], [0, 10, 25, 50, 75, 90, 100])
    print(f"{nm:7} min {q[0]:8.3f} p10 {q[1]:8.3f} p25 {q[2]:8.3f} med {q[3]:8.3f} p75 {q[4]:8.3f} p90 {q[5]:8.3f} max {q[6]:8.3f}")
print("\n-- C_T 上位 8 --")
for r in sorted(rows, key=lambda r: -r[6])[:8]:
    print(f"M_c {r[0]:.2f} f {r[1]:.3f} tr {r[2]:4.1f} tc {r[3]:3.1f} Lc {r[4]:4.2f} | L {r[5]:6.2f} C_T {r[6]:.4f} C_L {r[7]:+.4f} C_M {r[8]:+8.3f}")
