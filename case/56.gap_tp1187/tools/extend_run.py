#!/usr/bin/env python3
"""終わった 3D run を**新しい run ディレクトリ**へ index コピーで継続する。

すきま内部の発達は残差ノルムにほとんど乗らない (自由流に埋もれる) ので、
残差がプラトーでも派生量が `DRIFTING` のことがある。事前登録 `T4-3D-AB` の gating は
その場合「比較せず run を伸ばす」と定めている。これはそのための道具。

    python3 tools/extend_run.py --src run_0049_ab3d_wall --dst run_0051_ab3d_wall_long \\
        --steps 150000 [--out-int 10000] [--dry]

**同一バイナリで回すこと** (A/B の両腕を別ビルドで回すと比較が交絡する)。
`FORGE_CUDA_BLOCKSIZE` も元 run の `RUN_PROVENANCE.txt` から引き継ぐ。
"""
import argparse, glob, os, re, shutil, subprocess, sys
from pathlib import Path

CASE = Path(__file__).resolve().parents[1]
ROOT = CASE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"


def latest_res(run):
    fs = [(int(re.search(r"res_(\d+)\.h5$", f).group(1)), f)
          for f in glob.glob(str(run / "res_[0-9]*.h5"))]
    if not fs:
        raise SystemExit(f"{run}: res_*.h5 が無い")
    return Path(sorted(fs)[-1][1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--dst", required=True)
    ap.add_argument("--steps", type=int, required=True)
    ap.add_argument("--out-int", type=int, default=10000)
    ap.add_argument("--cfl", type=float, default=None,
                    help=("cfl と cfl_pseudo を差し替える。**変えたら収束先が同じか"
                          "確認すること**。推奨は node NS 定常本段で 4-6 "
                          "(implicitRelax 1.0 なら上限 6)"))
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    src, dst = CASE / a.src, CASE / a.dst
    if dst.exists():
        raise SystemExit(f"{dst} が既にある (連番が衝突している)")
    seed = latest_res(src)
    print(f"[継続] {a.src} の {seed.name} から {a.dst} へ")

    dst.mkdir()
    for f in ("mesh.h5", "solverConfig.yaml", "bcondConfig.yaml", "probe.yaml",
              "species_db.yaml", "inlet_profile_1.csv", "_inlet_table.txt",
              "case_setup.json", "stages.json"):
        p = src / f
        if p.exists():
            shutil.copy(p, dst / f)
    # **種を読んでから**古い出力を消す (混在事故の防止)
    subprocess.run([sys.executable, str(TOOLS / "restart_field.py"),
                    str(seed), str(dst / "mesh.h5")], check=True)
    for f in list(dst.glob("res_*")) + list(dst.glob("*.log")) + list(dst.glob("residual_history*")):
        f.unlink()

    cfg = (dst / "solverConfig.yaml").read_text(encoding="utf-8")
    cfg = re.sub(r"nStepOuter:\s*\d+", f"nStepOuter: {a.steps}", cfg)
    cfg = re.sub(r"outStepInterval:\s*\d+", f"outStepInterval: {a.out_int}", cfg)
    if a.cfl is not None:
        cfg = re.sub(r"cfl: [0-9.]+, cfl_pseudo: [0-9.]+",
                     f"cfl: {a.cfl:g}, cfl_pseudo: {a.cfl:g}", cfg)
        print(f"  cfl = cfl_pseudo = {a.cfl:g} に差し替え")
    (dst / "solverConfig.yaml").write_text(cfg, encoding="utf-8")
    m = re.search(r"nStepOuter:\s*(\d+)", cfg)
    print(f"  nStepOuter={m.group(1)}  outStepInterval={a.out_int}")

    bs = None
    prov = src / "RUN_PROVENANCE.txt"
    if prov.exists():
        mm = re.search(r"FORGE_CUDA_BLOCKSIZE=(\d+)", prov.read_text())
        bs = mm.group(1) if mm else None
    if bs:
        os.environ["FORGE_CUDA_BLOCKSIZE"] = bs
        print(f"  FORGE_CUDA_BLOCKSIZE={bs} (元 run から引き継ぎ)")
    else:
        print("  警告: 元 run の FORGE_CUDA_BLOCKSIZE が分からない")

    if a.dry:
        print("dry: 投入しない"); return
    rc = subprocess.run([str(TOOLS / "run_case.sh"), str(dst)]).returncode
    print(f"  run_case rc={rc}")


if __name__ == "__main__":
    main()
