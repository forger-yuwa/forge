#!/usr/bin/env python3
"""内側円柱の偏心 x_off を振って、CAD→メッシュ→計算→評価を一気に回す (ユーザ要件 2026-09-19)。

usage:
  python3 sweep_offset.py --offsets 0 0.5 1.0 1.5 [--stage D] [--run-prefix run_01] [--dry]
  python3 sweep_offset.py --offsets 0.5 --skip-mesh          # 既存メッシュで計算だけ

各点で:
  case.json の geometry.x_off を書き換え -> setup.py --resolve -> manifest.json
    -> cad/build_geom.py (FreeCAD)  -> cavity_fluid_half.step
    -> cad/mesh_salome.py (SALOME)  -> med -> med_to_msh41 -> forge.msh
    -> gen_runs.py convert          -> mesh/off<NNN>.h5   (品質 VERDICT を残す)
    -> gen_runs.py run              -> run_01NN_off<NNN>  (段階起動)
    -> tools/cavity_eval.py         -> cavity_eval.json

**偏心は x 方向のみ** (y=0 対称が保たれるので半割のまま)。すきまは周方向に
Ro-Ri-|x_off| 〜 Ro-Ri+|x_off| になり、VL 層数は最小すきまから自動で引き直される。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FREECAD = "/home/sano/opt/squashfs-root/usr/bin/freecadcmd"
SALOME = "/home/sano/opt/salome/SALOME-9.14.0-native-UB24.04-SRC/salome"
SALOME_LIB = "/home/sano/opt/salome/syslibs/usr/lib/x86_64-linux-gnu"
VENV_MESH = ROOT / ".venv-mesh" / "bin" / "python3"


def sh(cmd, cwd=None, env=None, log=None):
    e = dict(os.environ)
    if env:
        e.update(env)
    r = subprocess.run(cmd, cwd=cwd, env=e, capture_output=True, text=True, shell=isinstance(cmd, str))
    if log:
        Path(log).write_text(r.stdout + r.stderr)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offsets", type=float, nargs="+", required=True, help="x_off [mm]")
    ap.add_argument("--stage", default=None, help="mesh.stage (既定は case.json のまま)")
    ap.add_argument("--run-prefix", default="run_01")
    ap.add_argument("--main-steps", type=int, default=24000)
    ap.add_argument("--cfl", type=float, default=2.0)
    ap.add_argument("--skip-mesh", action="store_true")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    base = json.loads((HERE / "case.json").read_text())
    for i, off in enumerate(a.offsets):
        tag = "off%03d" % round(off * 100)
        print("\n================ x_off = %.2f mm  (%s) ================" % (off, tag), flush=True)
        cfg = json.loads(json.dumps(base))
        cfg["geometry"]["x_off"] = off
        if a.stage:
            cfg["mesh"]["stage"] = a.stage
        (HERE / "case.json").write_text(json.dumps(cfg, indent=2, ensure_ascii=False))
        r = sh([sys.executable, "setup.py", "--resolve"], cwd=HERE)
        print([l for l in r.stdout.splitlines() if l.startswith(("gap", "VL"))])
        if a.dry:
            continue

        mesh_h5 = HERE / "mesh" / ("%s.h5" % tag)
        if not a.skip_mesh or not mesh_h5.exists():
            r = sh([FREECAD, "build_geom.py"], cwd=HERE / "cad",
                   env={"QT_QPA_PLATFORM": "offscreen"}, log=HERE / "cad" / ("_geom_%s.log" % tag))
            if "wrote" not in r.stdout:
                print(r.stdout[-1500:]); sys.exit("build_geom 失敗 (%s)" % tag)
            med = "cavity_%s.med" % tag
            r = sh([SALOME, "--keep-paths", "-t", "mesh_salome.py",
                    "args:cavity_fluid_half.step,%s" % med], cwd=HERE / "cad",
                   env={"LD_LIBRARY_PATH": SALOME_LIB + ":" + os.environ.get("LD_LIBRARY_PATH", "")},
                   log=HERE / "cad" / ("_mesh_%s.log" % tag))
            nl = [l for l in r.stdout.splitlines() if l.startswith("nodes:")]
            print("  ", nl[0] if nl else "メッシュ失敗")
            if not nl or " tet 0 " in nl[0]:
                sys.exit("mesh 失敗 (tet 0 = VL 総厚 > 壁面サイズ の疑い): %s" % tag)
            r = sh([str(VENV_MESH), "med_to_msh41.py", med, "forge_%s.msh" % tag], cwd=HERE / "cad",
                   log=HERE / "cad" / ("_conv_%s.log" % tag))
            r = sh([sys.executable, "gen_runs.py", "convert", "--msh", "cad/forge_%s.msh" % tag,
                    "--out", "mesh/%s.h5" % tag], cwd=HERE)
            print("  ", [l for l in r.stdout.splitlines() if "VERDICT" in l])

        rd = HERE / ("%s%02d_%s" % (a.run_prefix, i, tag))
        if rd.exists():
            shutil.rmtree(rd)
        r = sh([sys.executable, "gen_runs.py", "run", "--run", rd.name, "--mesh", "mesh/%s.h5" % tag,
                "--inlet-csv", "mesh/inlet_profile.csv", "--main-steps", str(a.main_steps),
                "--cfl", str(a.cfl)], cwd=HERE, log=HERE / ("_sweep_%s.log" % tag))
        print("  run:", [l for l in r.stdout.splitlines() if l.startswith(("  [", "main rc"))][-3:])
        r = sh([sys.executable, "tools/cavity_eval.py", rd.name], cwd=HERE)
        for l in r.stdout.splitlines():
            if any(k in l for k in ("合計 Q", "dT_mouth", "dT_mid", "dT_floor", "dT_up", "dT_dn", "zpen(25")):
                print("   ", l.strip())

    (HERE / "case.json").write_text(json.dumps(base, indent=2, ensure_ascii=False))
    print("\ncase.json を元に戻しました")


if __name__ == "__main__":
    main()
