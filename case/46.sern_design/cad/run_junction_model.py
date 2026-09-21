#!/usr/bin/env python3
"""接続模型 (`hex_junction_model.py` の msh) を、生産 runner の設定・気体状態・段階起動で回す (plan tooling-sern-mesh-blocking B5 の前段)。

形状は MOC 設計ではなく試作の放物線ランプなので、**力の係数は生産値と比べられない**。見るのは
「有限厚の板と接続部で計算が成立するか」= 発散しないか、密度・圧力の床に張り付く節点が無いか、壁解像。
設定は問題 YAML (既定は生産 g1) を正本にし、格子と境界条件のタグだけを差し替える。

usage:  python3 run_junction_model.py PROBLEM.yaml JUNCTION.msh RUN_DIR [--prepare-only] [--stages full|none|main]
"""
import sys, json, subprocess, shutil, argparse
from pathlib import Path
import numpy as np
import h5py

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner_sern as R2          # noqa: E402
from forge_design.probdef import load_problem               # noqa: E402

MESH = R2.MESH
WALLS = ("ramp", "vehicle", "sidewall_in", "sidewall_out", "sidewall_end", "cowl_in", "cowl_out", "cowl_side", "cowl_base")
# 接続模型の形状 (hex_junction_model.py の P0 と同じ。IC の領域分けに使う)
H, ZW, TSW, TC = 0.1, 0.1, 0.005, 0.002


def phys_ids(msh):
    """msh4.1 の $PhysicalNames から 2D のタグ名 -> physID"""
    out = {}
    with open(msh) as f:
        for line in f:
            if line.startswith("$PhysicalNames"):
                for _ in range(int(next(f))):
                    d, t, n = next(f).split(maxsplit=2)
                    if d == "2": out[n.strip().strip('"')] = int(t)
                break
    return out


def bcond(p, st, ids):
    model = p.evaluate.get("model", "euler"); ex, en = st["exhaust"], st["ext"]
    okind = str(p.evaluate.get("outlet_kind", "statPress"))
    inlet = lambda n, s: (f"{n}: {{physID: {ids[n]}, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {s['ro']:.6g}, Ux: {s['u']:.6g}, "
                          f"Uy: 0.0, Uz: 0.0, Ps: {s['P']:.6g}, k: {s['k']:.6g}, omega: {s['omega']:.6g}{R2.inlet_species_floats(s)}}}}}\n")
    def outlet(n):
        if okind == "outflow": return f"{n}: {{physID: {ids[n]}, kind: outflow, outputHDFflg: 0, ints: , floats: }}\n"
        return f"{n}: {{physID: {ids[n]}, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {en['P']:.6g}, Pt: {en['P']:.6g}, Tt: {en['T']:.6g}}}}}\n"
    slip = lambda n: f"{n}: {{physID: {ids[n]}, kind: slip, outputHDFflg: 0, ints: , floats: }}\n"
    s = inlet("inlet", ex) + inlet("ext_in", en) + outlet("outlet") + outlet("far_bottom") + slip("far_side") + slip("symmetry")
    for w in WALLS:
        s += f"{w}: {p.wall_bcond_line(model == 'euler', phys_id=ids[w], output=1)}\n"
    return s


def prepare(problem, msh, run_dir, op=None):
    p = load_problem(problem); run_dir = Path(run_dir); run_dir.mkdir(parents=True, exist_ok=False)
    R2.select_operating_point(p, op); st = R2.gas_states(p)
    shutil.copy(msh, run_dir / "sern.msh"); ids = phys_ids(run_dir / "sern.msh")
    missing = [n for n in WALLS + ("inlet", "ext_in", "outlet", "far_bottom", "far_side", "symmetry") if n not in ids]
    if missing: raise RuntimeError(f"msh にタグが無い: {missing}")
    n = int(p.evaluate.get("nStepOuter", 4000)); out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    cfg = R2._solver_config(p, n, out_int, float(p.evaluate.get("cfl_main", 1.0)), st["ext"]["P"])
    (run_dir / "bcondConfig.yaml").write_text(bcond(p, st, ids)); (run_dir / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    R2.write_species_db(p, run_dir, R2.frozen_gases(p))
    (run_dir / "solverConfig.yaml").write_text(cfg); (run_dir / "solverConfig_main.yaml").write_text(cfg)
    R2.convert_mesh(run_dir, "sern.msh", MESH)
    tools = ROOT / "solver_density_cuda" / "tools"; txt = ""
    for cmd in (["check_mesh_quality.py", MESH, "--mode", "3d"], ["check_dual_closure.py", MESH]):
        q = subprocess.run([sys.executable, str(tools / cmd[0])] + cmd[1:], cwd=run_dir, capture_output=True, text=True); txt += q.stdout + q.stderr
        if cmd[0] == "check_mesh_quality.py" and q.returncode != 0:
            (run_dir / "MESH_QUALITY.txt").write_text(txt); raise RuntimeError("メッシュ品質 FAIL:\n" + q.stdout)
    (run_dir / "MESH_QUALITY.txt").write_text(txt)
    with h5py.File(run_dir / MESH, "r+") as f:             # 領域別一様 IC: 排気 = カウル上面より上 かつ 側壁内面より内側 (板厚の中央で分ける)
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        upper = (cc[:, 1] > -0.5 * TC) & (cc[:, 2] < ZW + 0.5 * TSW)
        R2.write_ic_arrays(f["/VALUE"], R2.region_ic_arrays(upper, st, p.gamma))
    info = {"problem": str(problem), "msh": str(msh), "run_dir": str(run_dir), "H_m": H, "states": st, "gas_model": st["gas_model"], "phys_ids": ids,
            "model": p.evaluate.get("model", "euler"), "dim": 3, "nodes": int(len(cc)), "n_exhaust_ic": int(upper.sum())}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1)); return info


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("problem"); ap.add_argument("msh"); ap.add_argument("run_dir")
    ap.add_argument("--prepare-only", action="store_true"); ap.add_argument("--stages", default="full", choices=["full", "none", "main"]); ap.add_argument("--op", default=None)
    a = ap.parse_args()
    if a.stages != "main":
        info = prepare(a.problem, a.msh, a.run_dir, a.op); print(json.dumps({k: info[k] for k in ("nodes", "n_exhaust_ic", "phys_ids", "gas_model")}, indent=1))
    if a.prepare_only: return 0
    o = (load_problem(a.problem).raw.get("opt") or {})
    return R2.run_staged(a.run_dir, a.stages, int(o.get("soft_steps", 2000)), soft_cfl=float(o.get("soft_cfl", 0.5)), soft_conv=int(o.get("soft_conv", 0)),
                         warm_lam_steps=int(o.get("warm_lam_steps", 0)), warm_lam_cfl=float(o.get("warm_lam_cfl", 0.2)), mid_steps=int(o.get("mid_steps", 0)))


if __name__ == "__main__":
    raise SystemExit(main())
