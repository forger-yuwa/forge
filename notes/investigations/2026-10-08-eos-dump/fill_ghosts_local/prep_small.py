"""小さい 3D SERN (g3 の構造・farfield の side_far・m10_on) を runner_sern3d.prepare で作る (eos_ab_dir の fill-ghosts の局所試験用)。
粗くする値は design/tests/run_sern_twin_ic_tests.py の (b) g3c と同じ。"""
import sys
from pathlib import Path
import yaml
ROOT = Path("/home/sano/work/forge-sern-design")
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
here = Path(__file__).resolve().parent
prob = ROOT / "case" / "46.sern_design" / "problem_3d_prod_3op_wallres_lswx08_tewake_g3_B10.yaml"
y = yaml.safe_load(open(prob))
y["mesh3d"].update({"first_wake_frac": 0.004, "first_z_frac": 0.02, "nj_vside": 5,
                    "first_wall_frac": 2.56e-3, "nz_in": 4, "nz_out": 3, "ni_noz": 30})
y["mesh"].update({"nj_ext_top": 9, "first_top_frac": 0.004})
pth = here / "prob_small_g3c.yaml"
pth.write_text(yaml.safe_dump(y, allow_unicode=True))
info = R3.prepare(pth, here / "run_small", op="m10_on")
print("nodes", info["mesh"]["nodes"], "evaluate side_far_kind", y["evaluate"]["side_far_kind"])
