import sys, time, warnings
sys.path.insert(0, "/home/sano/work/forge-thruster/design")
import numpy as np
from forge_design.gas.rocket_exhaust import RocketExhaust
from forge_design.geometry import rao_axisym as ra
EXH = "/home/sano/work/forge-thruster/case/68.bell_two_mode/exhaust_lox_lh2_pc100_of6.yaml"
ig = ra.as_isentrope(RocketExhaust.from_file(EXH).moc_gas())
line = ra._hall_line(1.0, ig, 161)
L = 0.7 * 22.657533530077966; L_march = ra.default_L_march(L); R_t = 1.0
dr = 3.2 / 161
def sB_of(c, thB):
    if c == 0.0: return thB
    disc = 1.0 + 2.0 * c * thB
    return None if disc < 0 else (-1.0 + np.sqrt(disc)) / c
def kside(net, bd, eta):
    d_eta = bd.dspl(eta); d_xi = bd.dxi(eta)
    return float((d_xi[2] + d_eta[2]) / np.hypot(d_xi[0] + d_eta[0], d_xi[1] + d_eta[1]))
cases = {0.0: (34.0, 36.09, 38.0), -0.35: (34.0, 36.5, 39.0), -0.7: (35.0, 37.5, 40.0)}
for cR, thBs in cases.items():
    c = cR / R_t**2
    rng = ra.sB_range(R_t, c)
    t0 = time.time()
    try:
        net = ra.kernel_field(ig, R_t, c, line, rng[1], L_march=L_march)
    except Exception as e:
        print("c", cR, "kernel_field failed", e); continue
    print(f"=== c*Rt^2 {cR:+.2f}: sB_range {rng} s_need {net.s_need:.4f} s_lin {net.wall.s_lin:.4f} K {net.K} ({time.time()-t0:.0f}s)")
    for thBd in thBs:
        sB = sB_of(c, np.deg2rad(thBd))
        if sB is None or sB > net.s_need or sB < net.s_wall[2]:
            print(f"  theta_B {thBd}: s_B {sB} outside net range [{net.s_wall[2]:.4f}, {net.s_need:.4f}]"); continue
        kB = 1.0 / R_t + c * sB
        t0 = time.time()
        bd = ra.trace_bd(net, sB)
        sd = ra.solve_D(bd, eps=50.0, p_a=0.0, dr=dr, n_scan=36)
        if not sd["roots"]:
            print(f"  theta_B {thBd:.3f} s_B {sB:.4f} kappa_B {kB:+.4f}: no inner root ({sd['why']}) jumps {len(sd['jumps'])} ({time.time()-t0:.0f}s)"); continue
        for r in sd["roots"]:
            D = r["D"]; E = r["E"]; rl = r["rl"]
            de = rl.dtheta_dsigma_D(); k = kside(net, bd, r["eta"])
            print(f"  theta_B {thBd:.3f} s_B {sB:.4f} kappa_B {kB:+.4f} | eta_D {r['eta']:.2f} det_D {rl.det_D:+.3f} D(x {D['x']:.3f} r {D['r']:.3f} th {np.rad2deg(D['th']):.2f} M {D['M']:.3f}) x_E-L {E['x']-L:+.3f} th_E {np.rad2deg(E['th']):.2f} | de {de:+.5f} kernel {k:+.5f} g~ {de-k:+.5f} ({time.time()-t0:.0f}s)")
