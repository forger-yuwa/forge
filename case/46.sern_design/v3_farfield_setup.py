#!/usr/bin/env python3
"""farfield plan (plans/accepted/boundary-node-farfield-characteristic.md §6 V3) の SERN run を作る (AWS で実行)。
R4d と同じ g3・同じ設定で、変えるのは side_far の種別 (slip / farfield) と遠方面の幅だけ。全 run 共通で:
  - 新バイナリ (build-ff、run_case.sh に FORGE_BIN で渡す)
  - リミッタ基準値を run_0986 の自動値に固定 (幅ごとの自動決定の交絡を除く、plan §5.1 #3 V2a の経緯)
  - 20000 step、500 step 出力 (R4d と同じ ε・窓条件で評価)
farfield の自由流は各 run の inlet_ext 行と同じ状態 (外気)。

  python3 v3_farfield_setup.py DST SRC_RUN --kind slip|farfield [--ic RES_H5 | --ic-from-mesh]
    --ic RES_H5      : 同一格子の最終場を restart_field.py で写す
    --ic-from-mesh   : SRC_RUN/sern.h5 の VALUE をそのまま使う (run_0989 = 3.42 H 最終場の共通部 index コピー + 追加層の最近傍)
"""
import argparse, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "..", "..", "solver_density_cuda", "tools")
LIM = "limiterRefLength: 2.136214053, limiterRoRef: 0.04544068206, limiterPRef: 4411.497933, limiterARef: 368.6412254"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dst"); ap.add_argument("src"); ap.add_argument("--kind", choices=("slip", "farfield"), required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--ic"); g.add_argument("--ic-from-mesh", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.dst)
    for f in ("solverConfig.yaml", "bcondConfig.yaml", "sern.h5", "probe.yaml", "species_db.yaml", "species_meta.yaml",
              "prepare_info.json", "cowl_contour.csv", "ramp_contour.csv", "MESH_QUALITY.txt"):
        if os.path.exists(os.path.join(a.src, f)):
            shutil.copy(os.path.join(a.src, f), a.dst)
    s = open(os.path.join(a.dst, "solverConfig.yaml")).read()
    s = re.sub(r"nStepOuter: \d+", "nStepOuter: 20000", s)
    s = re.sub(r"outStepInterval: \d+", "outStepInterval: 500", s)
    if "limiterRefLength" not in s:   # 延長 (V3 の run 自身が元) なら既に入っている
        s = s.replace("venkatK: 0.05}", f"venkatK: 0.05, {LIM}}}")
    assert LIM in s, "space 行の書き換えに失敗"
    open(os.path.join(a.dst, "solverConfig.yaml"), "w").write(s)
    b = open(os.path.join(a.dst, "bcondConfig.yaml")).read().splitlines()
    ext = [l for l in b if l.startswith("inlet_ext:")][0]
    side = [i for i, l in enumerate(b) if l.startswith("side_far:")][0]
    pid = re.search(r"physID: (\d+)", b[side]).group(1)
    if a.kind == "farfield":
        b[side] = re.sub(r"^inlet_ext: \{physID: \d+, kind: inlet_uniformVelocity", f"side_far: {{physID: {pid}, kind: farfield", ext)
        assert "kind: farfield" in b[side]
    open(os.path.join(a.dst, "bcondConfig.yaml"), "w").write("\n".join(b) + "\n")
    if a.ic:
        r = subprocess.run([sys.executable, os.path.join(TOOLS, "restart_field.py"), a.ic, os.path.join(a.dst, "sern.h5")], capture_output=True, text=True)
        print(r.stdout[-400:], r.stderr[-400:])
        if r.returncode != 0:
            raise SystemExit("restart_field に失敗")
        ic = f"restart_field.py {a.ic} -> sern.h5"
    else:
        ic = f"{a.src}/sern.h5 の VALUE をそのまま (R4d で作った初期場)"
    open(os.path.join(a.dst, "IC_FROM.txt"), "w").write(
        f"IC: {ic}。farfield plan §6 V3: side_far {a.kind}、設定は {a.src} と同一 + リミッタ基準値固定 ({LIM})、nStepOuter 20000\n")
    print(f"prepared {a.dst}: side_far {a.kind}\n  {b[side]}")


if __name__ == "__main__":
    main()
