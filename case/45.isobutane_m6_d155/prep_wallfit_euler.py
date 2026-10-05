"""補間壁 vs 位置+壁角の同時当てはめ壁の Euler A/B の run 準備 (solver は回さない)。plan verification-m6-axis-wave-mesh-su2 §5.1 #15。
problem_d155_euler_c2final_n2400.yaml (MOC 2400 点) から、arm=interp は現行の AxisMachCFDWall、arm=fit は同じ MOC 点群の
joint_fit (λ=1e-9, h0 0.0125, h1 0.5; moc_wall_fit_ab.py) で設計区間を置き換えた壁で、runner_axismach.prepare (Euler, slip, 段階起動 soft) を作る。
usage: python3 prep_wallfit_euler.py <run_dir> {interp|fit} --ic RUN     (準備のみ: 本段 cfl 6・implicitRelax 0.7・12000 step)
       python3 prep_wallfit_euler.py --run <run_dir>                     (準備済み run を soft 3000 → 本段で回す)
"""
import json, sys
from pathlib import Path
C = Path(__file__).resolve().parent
sys.path.insert(0, str(C)); sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402

if sys.argv[1] == "--run":
    rd = Path(sys.argv[2]); rc = RA.run_staged(rd, cfl_main=6.0, mid_stage=False, stages="soft")
    print(f"forge exit={rc}"); raise SystemExit(rc)
run_dir, arm = Path(sys.argv[1]), sys.argv[2]
ic = sys.argv[sys.argv.index("--ic") + 1] if "--ic" in sys.argv else None
if arm not in ("interp", "fit"):
    raise SystemExit("arm は interp か fit")
if arm == "fit":
    from moc_wall_fit_ab import joint_fit  # noqa: E402  (import で形状 A/B の測定も 1 回走る)
    _dc = RA.design_chain

    def design_chain_fit(p):
        d = _dc(p)
        s, nc = joint_fit(d["wall_inv"], d["R"], 1e-9)
        d["wall"]._spl = s
        d["wall_fit"] = {"kind": "joint_fit", "lam": 1e-9, "h0": 0.0125, "h1": 0.5, "n_cp": nc}
        msgs = d["wall"].validate()
        if msgs:
            print("validate:", msgs)
        return d
    RA.design_chain = design_chain_fit
info = RA.prepare(C / "problem_d155_euler_c2final_n2400.yaml", run_dir, nsteps=12000, ic_from=ic, cfl_main=6.0, implicit_relax=0.7)
info["stages"] = "soft"
info["wall_arm"] = arm
(run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
print(arm, json.dumps(info.get("mesh")), (run_dir / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1] if (run_dir / "MESH_QUALITY.txt").exists() else "no MESH_QUALITY")
