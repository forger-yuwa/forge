#!/usr/bin/env python3
"""#13-3 V5(i) 用: runner_axismach --prepare-only と同じ準備をしたうえで、メッシュを V0 (run_0509/0522 の nozzle.msh) に差し替える。

  FORGE_BIN=<新バイナリ> python3 prepare_v5i_v0mesh.py PROBLEM.yaml RUN_DIR V0_NOZZLE_MSH [--cfl 6 --implicit-relax 0.7]

背景: 段 3 で内蔵 Ar の第 2 区間が CEA に変わり、設計 (MOC) 壁が 1e-12 相対だけ動く。メッシュ生成で 1 節点が float32 の 1 ulp
(6e-8 m) 動き nozzle.msh の md5 が V0 と変わる。V5(i) は V0 と md5 一致のメッシュで回す (plan #13-3 合格条件 (4))。
手順: (1) runner_axismach.prepare (新データで config・bcond・IC・記録を作る) → (2) nozzle.msh を V0 のものに置き換え md5 を確認 →
(3) 同じ変換器で nozzle.h5 を作り直し → (4) runner と同じ関数で等エントロピー IC を貼り、化学種属性を付け直す。
runner のソースは変えない (design/forge_design/evaluate/runner_axismach.py の関数を呼ぶだけ)。
"""
import argparse, hashlib, shutil, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve()
ap = argparse.ArgumentParser()
ap.add_argument("problem"); ap.add_argument("run_dir"); ap.add_argument("v0_msh")
ap.add_argument("--repo", default="/home/sano/work/forge-species", help="design/ と solver_density_cuda/ のあるワークツリー")
ap.add_argument("--cfl", type=float, default=6.0); ap.add_argument("--implicit-relax", type=float, default=0.7)
a = ap.parse_args()
sys.path.insert(0, str(Path(a.repo) / "design"))
from forge_design.evaluate import runner_axismach as RA   # noqa: E402
from forge_design.evaluate.runner import _ENV, converter_path   # noqa: E402
from forge_design.probdef import load_problem   # noqa: E402


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


run = Path(a.run_dir)
RA.prepare(a.problem, run, cfl_main=a.cfl, implicit_relax=a.implicit_relax)
gen_md5 = md5(run / "nozzle.msh")
shutil.copyfile(a.v0_msh, run / "nozzle.msh")
assert md5(run / "nozzle.msh") == md5(a.v0_msh)
(run / "nozzle.h5").unlink()
subprocess.run([str(converter_path()), "nozzle.msh", "nozzle.h5"], cwd=run, env=_ENV, check=True, capture_output=True, text=True)
p = load_problem(a.problem)
d = RA.design_chain(p)
RA.paste_isentropic_ic(run / "nozzle.h5", d["wall"], float(p.spec["r_throat"]),
                       float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                       gas=(None if str(p.evaluate.get('cfd_gas', 'same')) == 'cpg' else p.gas_model),
                       h_ref_T=float(p.evaluate.get('thermo_href_temp', 298.15)),
                       species_Y=RA._tp_species_Y(p))
RA._stamp_ic_species(p, run)
(run / "V0_MESH.txt").write_text(f"nozzle.msh replaced by {a.v0_msh}\nmd5 V0 {md5(a.v0_msh)}\nmd5 generated (discarded) {gen_md5}\n"
                                 f"md5 in run {md5(run / 'nozzle.msh')}\n")
print((run / "V0_MESH.txt").read_text())
