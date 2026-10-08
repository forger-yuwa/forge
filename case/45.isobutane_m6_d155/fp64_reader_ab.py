"""FP64 の変換器の座標読み込みの A/B (plan tooling-nozzle-isothermal-wall-chain §5.1 #16、2026-10-08 事前登録、CFD 0 step、AWS)。
冷却壁の格子 (problem_d155_ns_prod_coldmesh.yaml の mesh、生産の物理壁 run_0167 の wall_repr.json) を mesh2d で倍精度に作り、msh (17 桁) に書き、
A = stof の FP64 変換器、B = stod の FP64 変換器で変換して、生成時の倍精度座標と比べる。
判定: B の第一層厚 (壁の各 station で壁節点と第一内部節点の距離) の相対誤差 ≤ 1e-6・非正の層厚なし → 採用。A は記録。
品質は check_mesh_quality.py --ar-max 5000 の VERDICT (厳密に PASS) を両方で記録する。
usage (AWS の case dir): python3 fp64_reader_ab.py <変換器 A> <変換器 B>  → _band_ab/cold_pair/fp64_reader_ab.json
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design"))
from forge_design.geometry.wall_axismach import load_wall_file  # noqa: E402
from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d  # noqa: E402

TOL = 1e-6


def main(conv_a: str, conv_b: str):
    W = load_wall_file(HERE / "run_0167_ns_n012_N2")
    ph = W["physical"]; d0, d1 = (float(v) for v in W["domain"]); rt = float(W["scale_m"])

    class Wall:
        x_in, x_e = d0, d1

        def r(self, x, d=0):
            return ph.r(np.asarray(x), d) if d else ph.r(np.asarray(x))
    m = yaml.safe_load(open(HERE / "problem_d155_ns_prod_coldmesh.yaml"))["mesh"]
    prm = Mesh2DParams(ni=m["ni"], nj=m["nj"], wall_first_frac=m["wall_first_frac"], throat_refine=m["throat_refine"],
                       throat_width=m["throat_width"], wall_first_frac_table=m["wall_first_frac_table"],
                       x_density_table=m["x_density_table"], wall_normal_layer=m["wall_normal_layer"], scale=rt)
    coords, quads, bed = generate_axisym_mesh(Wall(), prm)
    ni, nj = prm.ni, prm.nj
    G = coords[:, :2].reshape(ni, nj, 2)
    y1_gen = np.linalg.norm(G[:, -1] - G[:, -2], axis=1)
    out = {"ni": ni, "nj": nj, "tol": TOL, "msh_digits": int(m["msh_digits"]), "arms": {}}
    with tempfile.TemporaryDirectory(prefix="fp64ab_") as td:
        td = Path(td)
        write_msh41_2d(td / "nozzle.msh", coords, quads, bed, digits=int(m["msh_digits"]))
        for arm, conv in (("A_stof", conv_a), ("B_stod", conv_b)):
            wd = td / arm; wd.mkdir()
            shutil.copy2(td / "nozzle.msh", wd / "nozzle.msh")
            for f in (HERE / "run_0179_ns_n012_N2_ext").glob("*.yaml"):
                shutil.copy2(f, wd / f.name)
            r = subprocess.run([conv, "nozzle.msh", "nozzle.h5"], cwd=wd, capture_output=True, text=True)
            rec = {"converter": conv, "rc": r.returncode}
            if r.returncode != 0 or not (wd / "nozzle.h5").is_file():
                rec["error"] = (r.stdout + r.stderr)[-1500:]; out["arms"][arm] = rec; continue
            with h5py.File(wd / "nozzle.h5", "r") as h:
                C = np.array(h["MESH/COORD"]).reshape(-1, 3); rec["coord_dtype"] = str(h["MESH/COORD"].dtype)
                rec["volume_dtype"] = str(h["CELLS/volume"].dtype)
            same_order = C.shape[0] == coords.shape[0] and np.allclose(C[:, :2], coords[:, :2], rtol=0, atol=1e-6)
            rec["same_node_order"] = bool(same_order)
            if same_order:
                H = C[:, :2].reshape(ni, nj, 2)
                y1 = np.linalg.norm(H[:, -1] - H[:, -2], axis=1)
                rel = np.abs(y1 / y1_gen - 1.0)
                rec.update(coord_max_abs_err_m=float(np.abs(C[:, :2] - coords[:, :2]).max()),
                           first_layer_rel_err_max=float(rel.max()), first_layer_rel_err_median=float(np.median(rel)),
                           x_at_max=float(H[np.argmax(rel), -1, 0] / rt), nonpositive_layers=int(np.sum(np.diff(H[:, :, 1], axis=1) <= 0)))
            q = subprocess.run([sys.executable, str(HERE.parents[1] / "solver_density_cuda/tools/check_mesh_quality.py"), "nozzle.h5", "--ar-max", "5000"],
                               cwd=wd, capture_output=True, text=True)
            rec["mesh_quality"] = [l for l in q.stdout.splitlines() if l.startswith(("aspect", "skewness", "VERDICT"))]
            rec["ok"] = bool(same_order and rec.get("first_layer_rel_err_max", 1.0) <= TOL and rec.get("nonpositive_layers", 1) == 0
                             and any(l.strip() == "VERDICT: PASS (AR<= 5000, skew<= 0.90)" or l.startswith("VERDICT: PASS ") for l in rec["mesh_quality"]))
            out["arms"][arm] = rec
    out["verdict"] = "B 採用 (第一層厚を保持)" if out["arms"].get("B_stod", {}).get("ok") else "B 不合格 (冷却 NS は投入しない)"
    dst = HERE / "_band_ab" / "cold_pair"; dst.mkdir(parents=True, exist_ok=True)
    (dst / "fp64_reader_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
