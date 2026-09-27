#!/usr/bin/env python3
r"""case/61 の run ディレクトリを template から組む (実行はしない)。

    # 乾式 1 step (界面節点の確定用。共役なし)
    python3 case/61.conjugate_annulus/make_run.py dry  case/61.conjugate_annulus/run_0001_dry_r32 \
        --mesh case/61.conjugate_annulus/mesh/annulus_r32_x4.h5
    # 共役の本 run
    python3 case/61.conjugate_annulus/make_run.py cht  case/61.conjugate_annulus/run_0002_annulus \
        --mesh case/61.conjugate_annulus/mesh/annulus_r32_x4.h5 \
        --solid case/61.conjugate_annulus/mesh/solid_s16_x4.h5

実行は run ディレクトリで `FORGE_CUDA_BLOCKSIZE=128 <forge>` (本番は AWS の全域 FP64 ビルド)。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import axcht


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["dry", "cht"])
    ap.add_argument("run")
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--solid", default=None)
    ap.add_argument("--nstep", type=int, default=None, help="ローカルの起動確認だけに使う (nStepOuter を差し替え)")
    ap.add_argument("--note", default="")
    a = ap.parse_args()
    axcht.make_run(axcht.HERE / "template", Path(a.run), Path(a.mesh),
                   Path(a.solid) if a.solid else None, dry=(a.kind == "dry"), nstep=a.nstep, note=a.note)


if __name__ == "__main__":
    main()
