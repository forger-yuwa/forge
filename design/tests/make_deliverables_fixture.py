"""成果物の試験 (run_{mesh_structure,fluent_msh,wall_tables}_tests.py) の fixture を作る (plan tooling-nozzle-wt-deliverables §6)。

壁ファイルのある準備ディレクトリ (wall_repr.json・prepare_info.json・bcondConfig.yaml。例: case/45 の生産の確認
`_band_ab/prod_confirm/prep`) から、
- `nozzle.msh` を prepare_info.json の格子パラメータで mesh2d から作り直し、記録された生成器と msh の sha256 と一致することを確かめる
- `step/wall_physical.step` を wall_step で書く (FreeCAD が要る)
を OUT に置く。試験は環境変数 FORGE_DELIV_FIXTURE=OUT で使う。

usage (design/ を cwd に):
  .venv-opt/bin/python tests/make_deliverables_fixture.py PREP_DIR OUT_DIR
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.export import wall_step  # noqa: E402
from forge_design.geometry.wall_axismach import load_wall_file  # noqa: E402
from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d  # noqa: E402


def _upd(h, name, a):
    a = np.ascontiguousarray(a)
    h.update(name.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())


def main(prep, out) -> int:
    prep, out = Path(prep), Path(out)
    (out / "step").mkdir(parents=True, exist_ok=True)
    for n in ("wall_repr.json", "prepare_info.json", "bcondConfig.yaml", "wall_physical.csv", "wall_design.csv"):
        if (prep / n).exists():
            shutil.copy2(prep / n, out / n)
    info = json.loads((prep / "prepare_info.json").read_text())
    m = info["mesh"]
    W = load_wall_file(prep)
    coords, quads, bedges = generate_axisym_mesh(W["physical"], Mesh2DParams(**m["params"]))
    hg = hashlib.sha256()
    _upd(hg, "coords", np.asarray(coords, dtype=np.float64))
    _upd(hg, "quads", np.asarray(quads, dtype=np.int64))
    for g in ("inlet", "outlet", "wall", "axis"):
        _upd(hg, g, np.asarray(bedges[g], dtype=np.int64))
    write_msh41_2d(out / "nozzle.msh", coords, quads, bedges)
    ok_gen = hg.hexdigest() == m["hashes"]["generator_sha256"]
    ok_msh = hashlib.sha256((out / "nozzle.msh").read_bytes()).hexdigest() == m["hashes"]["msh_sha256"]
    print(f"generator sha256 一致: {ok_gen} / msh sha256 一致: {ok_msh}")
    data = wall_step.step_curve_data(W["record"])
    wall_step.write_step(data, out / "step" / "wall_physical.step")
    print(f"fixture: {out}")
    return 0 if (ok_gen and ok_msh) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
