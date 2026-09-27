#!/usr/bin/env python3
r"""case/62 の run ディレクトリを template から組む (実行はしない)。

    python3 case/62.conjugate_disk/make_run.py dry  case/62.conjugate_disk/run_0001_dry_r8u  --mesh case/62.conjugate_disk/mesh/disk_r8_u.h5
    python3 case/62.conjugate_disk/make_run.py cht  case/62.conjugate_disk/run_0005_disk_r8u \
        --mesh case/62.conjugate_disk/mesh/disk_r8_u.h5 --solid case/62.conjugate_disk/mesh/solid_disk_r8_u.h5
    python3 case/62.conjugate_disk/make_run.py axisprobe case/62.conjugate_disk/run_0009_axisprobe \
        --mesh case/62.conjugate_disk/mesh/disk_axisprobe_r16.h5

実行は run ディレクトリで `FORGE_CUDA_BLOCKSIZE=128 <forge>` (本番は AWS の全域 FP64 ビルド)。
`axisprobe` は template/*_axisprobe.yaml (共役なし 1 step、axisRFloor 0.5 mm)。
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "61.conjugate_annulus"))
import axcht  # noqa: E402

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["dry", "cht", "axisprobe"])
    ap.add_argument("run")
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--solid", default=None)
    ap.add_argument("--nstep", type=int, default=None, help="ローカルの起動確認だけに使う (nStepOuter を差し替え)")
    ap.add_argument("--note", default="")
    a = ap.parse_args()
    tdir = HERE / "template"
    if a.kind == "axisprobe":
        # dry と同じ組み方で、config だけ軸負例のものに差し替える
        axcht.make_run(tdir, Path(a.run), Path(a.mesh), None, dry=True, nstep=a.nstep, note=a.note or "軸負例")
        shutil.copy(tdir / "solverConfig_axisprobe.yaml", Path(a.run) / "solverConfig.yaml")
        shutil.copy(tdir / "bcondConfig_axisprobe.yaml", Path(a.run) / "bcondConfig.yaml")
        with open(Path(a.run) / "RUN_INPUTS.txt", "a") as f:
            f.write("config: template/solverConfig_axisprobe.yaml / bcondConfig_axisprobe.yaml\n")
        return
    axcht.make_run(tdir, Path(a.run), Path(a.mesh), Path(a.solid) if a.solid else None,
                   dry=(a.kind == "dry"), nstep=a.nstep, note=a.note)


if __name__ == "__main__":
    main()
