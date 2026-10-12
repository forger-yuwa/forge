#!/usr/bin/env python3
# ③ベル (M): 理論式 contur_v2 の縁を MOC の壁の M (x ≥ 0) + 凍結源 Euler run_0015 の壁 (x < 0、収縮部は全 (M) 共通) にしたときと、
# 同じ壁の Euler の縁 (較正で使った取り方) のときの δ_r の差 — plan tooling-nozzle-bell-two-mode §5.1 #7c の案 (2026-10-13)。
# 実行: /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/check_edge_moc_vs_euler.py
# MOC の壁の M (x ≥ 0) + 凍結源 Euler run_0015 の壁 (x < 0) を縁にした contur_v2 と、同じ壁の Euler の縁の contur_v2 の比較
import sys, numpy as np
sys.path.insert(0, '/home/sano/work/forge-thruster/design')
from pathlib import Path
from forge_design.metrics import bell_deltastar as bd
from forge_design.evaluate.runner import build_wall, rao_bell_spec
from forge_design.evaluate.runner_axismach import contur_mu_fn
from forge_design.feedback.deltastar_integral import integral_bl
from forge_design.probdef import load_problem
C68 = Path('/home/sano/work/forge-thruster/case/68.bell_two_mode')
def run(problem, wall_npz, eu_run, src_run='run_0015_tp_m_euler_l100'):
    p = load_problem(C68 / problem); wall, rec = build_wall(p); gas = p.gas_model; rt = float(p.spec['r_throat'])
    up = bd.euler_wall_edge(C68 / src_run, gas)
    z = np.load(C68 / wall_npz); W = z['wall']
    m_up = up['x'] < 0.0
    xM = np.r_[up['x'][m_up], W[:, 0]]; MM = np.r_[up['M'][m_up], W[:, 3]]
    o = np.argsort(xM); xM, MM = xM[o], MM[o]; keep = np.r_[True, np.diff(xM) > 1e-9]; xM, MM = xM[keep], MM[keep]
    xe_ = float(wall.x_e); Me_ = float(np.interp(xe_, xM, MM)); k_ = xM < xe_ - 1e-9; xM = np.r_[xM[k_], xe_]; MM = np.r_[MM[k_], Me_]
    xF = xe_
    rwt = np.asarray(wall.r(np.clip(xM, wall.x_in, wall.x_e)))
    wall_tbl = np.c_[xM, rwt, np.zeros(len(xM)), MM]
    x0 = float(wall.x_in) + 1e-6
    ref = np.load(C68 / 'rao_m' / ({'run_0018_tp_m_euler_l100_pin': 'contur_l100_euler_run0018.npz', 'run_0023_tp_m_euler_l70_pin': 'contur_l70_euler_run0023.npz'}[eu_run]), allow_pickle=True)
    xo = np.clip(ref['x'][(ref['x'] >= x0 - 2e-6) & (ref['x'] <= xF)], x0, xF)
    res = integral_bl(wall, wall_tbl, gas, p.cp, float(p.spec['Pt']), float(p.spec['Tt']), rt, thermal_bc=p.wall_thermal_bc_integral,
                      x_virtual_m=float(rao_bell_spec(p).L_pipe) * rt, x_out=xo, rtol=1e-6, closure_version='contur_v2',
                      mu_fn=contur_mu_fn(p), edge_n_grid=4000, edge_M_full=np.c_[xM, MM])
    rw = np.asarray(wall.r(res['x']))
    Cm = bd.delta_r_annular(res['dstar_n'], res['cos_phi'], rw)
    Ce = np.interp(res['x'], ref['x'], ref['delta_r_ann'])
    r = Cm / Ce - 1
    print(problem, 'M jump at 0: euler_up', float(np.interp(-1e-6, up['x'], up['M'])), 'moc', float(W[0, 3]))
    for lo, hi in [(-4, 0), (0, 0.41), (0.41, 1.17), (1.17, 2), (2, 8), (8, 16), (16, 99)]:
        mm = (res["x"] >= lo) & (res["x"] < hi)
        if not mm.any(): continue
        print(f'  [{lo},{hi}) C_moc/C_euler-1 mean {np.mean(r[mm])*100:+.3f}% max|{np.max(abs(r[mm]))*100:.3f}%')
    print('  exit', Cm[-1], Ce[-1], (Cm[-1] / Ce[-1] - 1) * 100, '%')
run('problem_bell_tp_m_ns_l100.yaml', 'rao_m/rao_pin_parts/L100_pin_n321_wall.npz', 'run_0018_tp_m_euler_l100_pin')
run('problem_bell_tp_m_ns_l70.yaml', 'rao_m/rao_pin_parts/L70_pin_n321_wall.npz', 'run_0023_tp_m_euler_l70_pin')
