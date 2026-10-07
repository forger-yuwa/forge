"""plan tooling-nozzle-throat-monotone-r2 の result 段レビュー (2026-10-07 M1) の対応: 実際の生産の入口を通す 0 step の準備 A/B (AWS で回す)。
基準 A = run_0147_ns_mono_final (検証に使った NS の入力)。単調壁の生産問題 YAML (problem_d155_ns_finemesh_recal_final_mono.yaml) を、
  B1: 検証の NS レシピ = prep_c2pin.py RUN K_F --problem YAML --stages none + run_mono_ns_chain.sh と同じ CFL 1・60000 step・5000 ごとへの書き換え
  B2: 設計チェーンの入口 = python -m forge_design.feedback.deltastar_loop --init-integral --prepare-only (k_f は YAML の deltastar_initializer)
で準備し (forge は起動しない)、A と比べる:
  実効 config (solverConfig.yaml・bcondConfig.yaml を辞書で)、prepare_info の当てはめ後 spline (t・c・k・mono_r2)、
  delta_r_initial.csv、wall_physical.csv、nozzle.h5 の MESH/COORD・MESH/CONNE (VALUE は IC なので比べない)。
判定 (事前に決める): B1 は config・幾何とも A とビット一致 → 「検証レシピは YAML + run_mono_ns_chain.sh で再現する」。
  B2 は幾何 (spline・δ_r・物理壁・メッシュ) が A とビット一致 → 「設計チェーンの入口も同じ物理壁を作る」。B2 の config (段・CFL・step) は
  pass 0 の設定で検証レシピとは違うので、違いを列挙するだけ (NS の数値レシピは run_mono_ns_chain.sh が持つ)。
usage (AWS, case dir): python3 throat_mono_entry_prep_ab.py [WORK_DIR (既定 _entry_prep_ab)] → _band_ab/throat_mono_entry_prep_ab.json
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

C = Path(__file__).resolve().parent
A = C / "run_0147_ns_mono_final"
PN = "problem_d155_ns_finemesh_recal_final_mono.yaml"
EU = "run_0143_euler_wallfit_monoG1_r1"
WORK = C / (sys.argv[1] if len(sys.argv) > 1 else "_entry_prep_ab")
KF = json.loads((C / "c2pin_solve_recal.json").read_text())["k_f"]


def sh(cmd, cwd=C, env=None):
    q = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=env)
    if q.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))} rc={q.returncode}\n{q.stdout[-1500:]}\n{q.stderr[-1500:]}")
    return q.stdout


def chain_edit(rd: Path):
    """run_mono_ns_chain.sh と同じ書き換え (sed と同じ正規表現)。"""
    f = rd / "solverConfig.yaml"; s = f.read_text()
    s = re.sub(r"cfl: [0-9.]+, cfl_pseudo: [0-9.]+", "cfl: 1.0, cfl_pseudo: 1.0", s)
    s = re.sub(r"nStepOuter: [0-9]+", "nStepOuter: 60000", s)
    s = re.sub(r"outStepInterval: [0-9]+", "outStepInterval: 5000", s)
    f.write_text(s)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def dict_diff(a, b, path=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: {'無' if k not in a else a[k]!r} ≠ {'無' if k not in b else b[k]!r}")
            else:
                out += dict_diff(a[k], b[k], f"{path}/{k}")
    elif a != b:
        out.append(f"{path}: {a!r} ≠ {b!r}")
    return out


def compare(rd: Path) -> dict:
    res = {}
    for f in ("solverConfig.yaml", "bcondConfig.yaml"):
        da, db = yaml.safe_load((A / f).read_text()), yaml.safe_load((rd / f).read_text())
        res[f] = dict_diff(da, db) or "identical"
    ia, ib = json.loads((A / "prepare_info.json").read_text()), json.loads((rd / "prepare_info.json").read_text())
    sa, sb = (ia.get("wall_fit") or {}).get("spline"), (ib.get("wall_fit") or {}).get("spline")
    res["wall_fit_spline"] = ("identical" if sa == sb and sa is not None else
                              {"A_has": sa is not None, "B_has": sb is not None,
                               "max_abs_dc": (float(np.abs(np.asarray(sa["c"]) - np.asarray(sb["c"])).max())
                                              if sa and sb and len(sa["c"]) == len(sb["c"]) else None)})
    res["mono_r2"] = [(ia.get("wall_fit") or {}).get("mono_r2"), (ib.get("wall_fit") or {}).get("mono_r2")]
    res["dstar_source"] = [ia.get("dstar_source"), ib.get("dstar_source")]
    res["initializer"] = [ia.get("initializer"), ib.get("initializer")]
    for f in ("delta_r_initial.csv", "wall_physical.csv"):
        if (A / f).is_file() and (rd / f).is_file():
            xa = np.loadtxt(A / f, delimiter=",", skiprows=1); xb = np.loadtxt(rd / f, delimiter=",", skiprows=1)
            res[f] = ("identical" if sha(A / f) == sha(rd / f) else
                      {"shape": [list(xa.shape), list(xb.shape)],
                       "max_abs_diff": float(np.abs(xa - xb).max()) if xa.shape == xb.shape else None})
        else:
            res[f] = {"A_exists": (A / f).is_file(), "B_exists": (rd / f).is_file()}
    with h5py.File(A / "nozzle.h5") as fa, h5py.File(rd / "nozzle.h5") as fb:
        for k in ("MESH/COORD", "MESH/CONNE"):
            va, vb = np.asarray(fa[k]), np.asarray(fb[k])
            res[k] = ("identical" if va.shape == vb.shape and np.array_equal(va, vb) else
                      {"shape": [list(va.shape), list(vb.shape)],
                       "max_abs_diff": float(np.abs(va.astype(float) - vb.astype(float)).max()) if va.shape == vb.shape else None})
    return res


def main():
    if WORK.exists():
        raise SystemExit(f"{WORK} が既にある — 消してから回す (既存の結果を上書きしない)")
    WORK.mkdir()
    env = dict(os.environ)
    b1, b2 = WORK / "B1_prep_c2pin", WORK / "B2_deltastar_loop"
    sh([sys.executable, "prep_c2pin.py", str(b1), repr(KF), "--problem", PN, "--stages", "none"], env=env)
    chain_edit(b1)
    sh([sys.executable, "-m", "forge_design.feedback.deltastar_loop", "--problem", str(C / PN), "--euler-ref", str(C / EU),
        "--init-integral", "--prepare-only", "--run-dir", str(b2)], cwd=C.parents[1] / "design", env=env)
    out = {"A": A.name, "problem": PN, "k_f": KF, "B1_prep_c2pin_plus_chain_edit": compare(b1), "B2_deltastar_loop_init_integral": compare(b2)}
    geo = ("wall_fit_spline", "delta_r_initial.csv", "wall_physical.csv", "MESH/COORD", "MESH/CONNE")
    out["verdict"] = {
        "B1": "config・幾何とも A とビット一致" if all(out["B1_prep_c2pin_plus_chain_edit"][k] == "identical"
                                                for k in geo + ("solverConfig.yaml", "bcondConfig.yaml")) else "A と不一致 (上の差を見る)",
        "B2_geometry": "幾何は A とビット一致" if all(out["B2_deltastar_loop_init_integral"][k] == "identical" for k in geo) else "幾何が A と不一致",
    }
    (C / "_band_ab").mkdir(exist_ok=True)
    (C / "_band_ab/throat_mono_entry_prep_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print(json.dumps(out, indent=1, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
