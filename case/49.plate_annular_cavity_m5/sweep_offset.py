#!/usr/bin/env python3
"""内側円柱の偏心 x_off を振って、メッシュ生成→変換→計算→評価を一気に回す (ユーザ要件 2026-09-19)。

usage:
  python3 sweep_offset.py --offsets 0 0.5 1.0 1.5 [--scale 1.0] [--run-prefix run_02] [--dry]
  python3 sweep_offset.py --offsets 0.5 --skip-mesh          # 既存メッシュで計算だけ
  python3 sweep_offset.py --offsets 0.5 --mesher salome      # 旧 tet+prism 経路 (ローカル専用)

各点で (既定 = ヘキサ経路):
  case.json の geometry.x_off を書き換え -> setup.py --resolve -> manifest.json
    -> cad/build_hex_mesh.py --scale S   -> cad/hex_off<NNN>.msh   (gmsh API のみ)
    -> gen_runs.py convert               -> mesh/off<NNN>.h5       (品質 VERDICT を残す)
    -> tools/check_mesh_extra.py         -> 追加ゲート
    -> gen_runs.py run                   -> run_<prefix><NN>_off<NNN>  (段階起動)
    -> tools/cavity_eval.py              -> cavity_eval.json

**偏心は x 方向のみ** (y=0 対称が保たれるので半割のまま)。すきまは周方向に
Ro-Ri-|x_off| 〜 Ro-Ri+|x_off| になる。ヘキサ経路は外筒軸からのレイが内円と交わる
半径 r_i(φ) を節点座標に使うだけなので、偏心でもブロッキングは変わらない
(tet 経路は VL 層数を最小すきまから引き直すため、偏心ごとに格子族が変わってしまう)。

**ヘキサ経路は FreeCAD も SALOME も要らない** → AWS 上で完結できる (`pip install gmsh`)。
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


def venv_python():
    """gmsh が入っている python (ローカル .venv-mesh / AWS ~/venv-mesh / 現行 python)。"""
    for p in (ROOT / ".venv-mesh" / "bin" / "python3",
              Path.home() / "venv-mesh" / "bin" / "python3"):
        if p.exists():
            return str(p)
    return sys.executable


def sh(cmd, cwd=None, env=None, log=None):
    e = dict(os.environ)
    if env:
        e.update(env)
    r = subprocess.run(cmd, cwd=cwd, env=e, capture_output=True, text=True, shell=isinstance(cmd, str))
    if log:
        Path(log).write_text(r.stdout + r.stderr)
    return r


def grep(out, *keys):
    return [l.strip() for l in out.splitlines() if any(k in l for k in keys)]


def build_hex(tag, scale, log):
    """gmsh API で全ヘキサ・メッシュを作る。戻り値 = .msh のパス (失敗なら None)。"""
    msh = HERE / "cad" / ("hex_%s.msh" % tag)
    r = sh([venv_python(), "build_hex_mesh.py", "--out", str(msh), "--scale", str(scale)],
           cwd=HERE / "cad", log=log)
    nl = grep(r.stdout, "節点")
    print("  ", nl[-1] if nl else "メッシュ失敗 (log: %s)" % log)
    if not msh.exists():
        print(r.stdout[-1500:])
        return None
    return msh


def build_salome(tag, log_dir):
    """旧 tet+prism 経路 (FreeCAD -> SALOME -> med -> msh4.1)。ローカル専用。"""
    r = sh([FREECAD, "build_geom.py"], cwd=HERE / "cad",
           env={"QT_QPA_PLATFORM": "offscreen"}, log=log_dir / ("_geom_%s.log" % tag))
    if "wrote" not in r.stdout:
        print(r.stdout[-1500:])
        return None
    med = "cavity_%s.med" % tag
    r = sh([SALOME, "--keep-paths", "-t", "mesh_salome.py",
            "args:cavity_fluid_half.step,%s" % med], cwd=HERE / "cad",
           env={"LD_LIBRARY_PATH": SALOME_LIB + ":" + os.environ.get("LD_LIBRARY_PATH", "")},
           log=log_dir / ("_mesh_%s.log" % tag))
    nl = grep(r.stdout, "nodes:")
    print("  ", nl[0] if nl else "メッシュ失敗")
    if not nl or " tet 0 " in nl[0]:
        print("   tet 0 = VL 総厚 > 壁面サイズ の疑い")
        return None
    msh = HERE / "cad" / ("forge_%s.msh" % tag)
    sh([venv_python(), "med_to_msh41.py", med, msh.name], cwd=HERE / "cad",
       log=log_dir / ("_conv_%s.log" % tag))
    return msh if msh.exists() else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offsets", type=float, nargs="+", required=True, help="x_off [mm]")
    ap.add_argument("--mesher", choices=("hex", "salome"), default="hex")
    ap.add_argument("--scale", type=float, default=1.0, help="ヘキサの細分率 (build_hex_mesh --scale)")
    ap.add_argument("--stage", default=None, help="mesh.stage (tet 経路のみ。既定は case.json のまま)")
    ap.add_argument("--run-prefix", default="run_02")
    ap.add_argument("--main-steps", type=int, default=24000)
    ap.add_argument("--cfl", type=float, default=2.0)
    ap.add_argument("--gas", default="CPG")
    ap.add_argument("--skip-mesh", action="store_true")
    ap.add_argument("--skip-run", action="store_true", help="評価だけやり直す")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    base = json.loads((HERE / "case.json").read_text())
    summary = []
    try:
        for i, off in enumerate(a.offsets):
            tag = "off%03d" % round(off * 100)
            print("\n================ x_off = %.2f mm  (%s) ================" % (off, tag), flush=True)
            cfg = json.loads(json.dumps(base))
            cfg["geometry"]["x_off"] = off
            if a.stage:
                cfg["mesh"]["stage"] = a.stage
            (HERE / "case.json").write_text(json.dumps(cfg, indent=2, ensure_ascii=False))
            r = sh([sys.executable, "setup.py", "--resolve"], cwd=HERE)
            print("  ", grep(r.stdout, "gap", "VL "))
            if a.dry:
                continue

            mesh_h5 = HERE / "mesh" / ("%s.h5" % tag)
            if not (a.skip_mesh and mesh_h5.exists()):
                if a.mesher == "hex":
                    msh = build_hex(tag, a.scale, HERE / "cad" / ("_hex_%s.log" % tag))
                else:
                    msh = build_salome(tag, HERE / "cad")
                if msh is None:
                    print("   メッシュ失敗 -> この偏心点をスキップ")
                    continue
                r = sh([sys.executable, "gen_runs.py", "convert", "--msh",
                        str(msh.relative_to(HERE)), "--out", "mesh/%s.h5" % tag], cwd=HERE,
                       log=HERE / ("_conv_%s.log" % tag))
                print("  ", grep(r.stdout, "VERDICT"))
                r = sh([sys.executable, "tools/check_mesh_extra.py", "mesh/%s.h5" % tag], cwd=HERE)
                print("  ", grep(r.stdout, "VERDICT"))

            rd = HERE / ("%s%02d_%s" % (a.run_prefix, i, tag))
            if not a.skip_run:
                if rd.exists():
                    shutil.rmtree(rd)
                r = sh([sys.executable, "gen_runs.py", "run", "--run", rd.name,
                        "--mesh", "mesh/%s.h5" % tag, "--inlet-csv", "mesh/inlet_profile.csv",
                        "--gas", a.gas, "--main-steps", str(a.main_steps), "--cfl", str(a.cfl)],
                       cwd=HERE, log=HERE / ("_sweep_%s.log" % tag))
                print("  run:", grep(r.stdout, "  [", "main rc")[-3:])
            r = sh([sys.executable, "tools/cavity_eval.py", rd.name], cwd=HERE)
            keep = grep(r.stdout, "合計 Q", "dT_mouth", "dT_mid", "dT_floor", "dT_up", "dT_dn",
                        "zpen(25", "imbalance", "CV 収支")
            for l in keep:
                print("   ", l)
            summary.append((off, tag, rd.name, keep))
    finally:
        (HERE / "case.json").write_text(json.dumps(base, indent=2, ensure_ascii=False))
        print("\ncase.json を元に戻しました")

    if summary:
        print("\n================ 偏心スイープ まとめ ================")
        for off, tag, rn, keep in summary:
            print("x_off %.2f mm  [%s]" % (off, rn))
            for l in keep:
                print("   ", l)


if __name__ == "__main__":
    main()
