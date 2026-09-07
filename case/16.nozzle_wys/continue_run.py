#!/usr/bin/env python3
"""同一メッシュ継続 run の生成: SRC_RUN の最終 res を新 run の mesh h5 に **index コピー**して nStepOuter=N で本段のみ実行。
(旧版は interp_field で貼っていたが、(x,y) 最近傍は 3D で z 列を混同するため index コピーに変更 2026-09-08)
usage: python3 continue_run.py SRC_RUN DST_RUN N"""
import sys, shutil, re, h5py, yaml
from pathlib import Path
CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "design"))
from forge_design.evaluate.runner import run_forge
src, dst, n = CASE / sys.argv[1], CASE / sys.argv[2], int(sys.argv[3])
dst.mkdir(exist_ok=False)
mesh = yaml.safe_load((src / "solverConfig.yaml").read_text())["mesh"]["meshFileName"]
for f in ["solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "probe.yaml", mesh]:
    if (src / f).exists():
        shutil.copy(src / f, dst / f)
res = sorted(src.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
with h5py.File(res, "r") as s, h5py.File(dst / mesh, "r+") as d:
    nn = len(d["VALUE/ro"])
    for k in s["VALUE"]:
        if k != "wall_dist" and k in d["VALUE"] and len(s["VALUE"][k]) == nn:
            d["VALUE"][k][:] = s["VALUE"][k][:]
c = (dst / "solverConfig.yaml").read_text()
c = re.sub(r"nStepOuter: \d+", f"nStepOuter: {n}", c); c = re.sub(r"outStepInterval: \d+", f"outStepInterval: {n//4}", c)
(dst / "solverConfig.yaml").write_text(c)
(dst / "CONTINUED_FROM.txt").write_text(f"{src.name}/{res.name}\n")
print("main rc", run_forge(dst))
print((dst / "CONVERGENCE_VERDICT.txt").read_text().strip().splitlines()[-1])
