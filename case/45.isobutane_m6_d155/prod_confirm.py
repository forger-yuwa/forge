"""生産の問題 (problem_d155_ns_prod.yaml) の準備が、生産の run (N2 の run_0167_ns_n012_N2) の入力を再現することの確認 (CFD 0 step)。
2026-10-07 ユーザ決定「生産採用: 上流を多項式に・MOC の新方式・全域 1 本の B スプライン壁」(plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md §9)。

  prepare : prepare_ns を ns_n012.py prep-dry と同じ引数で (IC は貼らない) → _band_ab/prod_confirm/prep
  compare : run_0167 と比べる → _band_ab/prod_confirm/PROD_CONFIRM.json
            合格: nozzle.h5 の幾何 (/MESH・/CELLS・/PLANES・/BCONDS の VALUE 以外) がビット同一、wall_design.csv・delta_r_initial.csv・
            bcondConfig.yaml がバイト同一、壁ファイル (wall_repr.json、single_bspline) がある。wall_physical.csv と throat_physical の差は記録
            (生産は 1 本の B-spline から評価、run_0167 は区分表現から評価。全域 1 本の plan の W1 で差は 1e-14 r_t 級)。

usage (AWS の case dir で): FORGE_BIN=... REAL_CONVERTER=... FORGE_CONVERTER=$(pwd)/conv_tolerant.sh python3 prod_confirm.py prepare
                           python3 prod_confirm.py compare
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ns_n012 as NS  # noqa: E402
import wsb_verify as WV  # noqa: E402

OUT = HERE / "_band_ab" / "prod_confirm"
PROB = HERE / "problem_d155_ns_prod.yaml"
REF = HERE / "run_0167_ns_n012_N2"


def prepare():
    import os
    if not os.environ.get("FORGE_CONVERTER"):
        raise SystemExit("FORGE_CONVERTER を指定する")
    d = OUT / "prep"
    if d.exists():
        raise SystemExit(f"{d} が既にある")
    OUT.mkdir(parents=True, exist_ok=True)
    kf = NS.k_f_checked(PROB)
    NS._prepare_ns(PROB, d, nsteps=12000, ic_from=None, initializer=NS.initializer(kf), cfl_main=5.0, implicit_relax=NS.RELAX)
    rec = {"problem": PROB.name, "problem_sha256": NS.sha256_file(PROB), "k_f": kf, "git_head": NS.git_head(), "created": NS.now()}
    NS.jdump(OUT / "prepare_record.json", rec)
    print(f"[prod] prepared {d}")


def compare():
    from forge_design.geometry.wall_axismach import WALL_FILE, load_wall_file
    d = OUT / "prep"
    cmp = WV.h5_compare(d / "nozzle.h5", REF / "nozzle.h5")
    geo = {k: v["identical"] for k, v in cmp["datasets"].items()
           if k.startswith(("MESH/", "CELLS/", "PLANES/")) or (k.startswith("BCONDS/") and "/VALUE/" not in k)}
    files = {n: NS.sha256_file(d / n) == NS.sha256_file(REF / n) for n in ("wall_design.csv", "delta_r_initial.csv", "bcondConfig.yaml")}
    wa = np.loadtxt(d / "wall_physical.csv", delimiter=",", skiprows=1)
    wb = np.loadtxt(REF / "wall_physical.csv", delimiter=",", skiprows=1)
    ia, ib = NS.jload(d / "prepare_info.json"), NS.jload(REF / "prepare_info.json")
    W = load_wall_file(d / WALL_FILE) if (d / WALL_FILE).is_file() else None
    out = {"item": "生産の問題の準備が run_0167 の入力を再現するか (CFD 0 step)", "problem": PROB.name, "ref": REF.name,
           "nozzle_h5_geometry": {"n": len(geo), "all_identical": bool(geo) and all(geo.values()), "diff": [k for k, v in geo.items() if not v]},
           "files_identical": files,
           "wall_physical_csv": {"same_shape": wa.shape == wb.shape, "max_abs_diff_m": (float(np.abs(wa - wb).max()) if wa.shape == wb.shape else None),
                                 "x_identical": bool(wa.shape == wb.shape and np.array_equal(wa[:, 0], wb[:, 0]))},
           "throat_physical": {"prod": ia.get("throat_physical"), "run_0167": ib.get("throat_physical")},
           "physical_wall_repr": (ia.get("physical_wall") or {}).get("repr"), "wall_file": bool(W),
           "moc": {k: (ia.get("moc") or {}).get(k) for k in ("axis_limit", "corrector")}, "moc_gate_pass": ((ia.get("moc") or {}).get("gate") or {}).get("pass"),
           "pw_upstream": (ia.get("pw_upstream") or {}).get("value"), "Md_moc_offset": ia.get("Md_moc_offset"), "scale_m": ia.get("scale_m")}
    ok = (out["nozzle_h5_geometry"]["all_identical"] and all(files.values()) and out["physical_wall_repr"] == "single_bspline" and out["wall_file"]
          and out["moc_gate_pass"] is True and ia.get("scale_m") == ib.get("scale_m") and ia.get("Md_moc_offset") == ib.get("Md_moc_offset"))
    out["verdict"] = "PASS (生産の問題の準備は run_0167 の入力を再現)" if ok else "FAIL (差の経路を特定して止める)"
    NS.jdump(OUT / "PROD_CONFIRM.json", out)
    print(json.dumps({k: out[k] for k in ("verdict", "nozzle_h5_geometry", "files_identical", "wall_physical_csv", "physical_wall_repr", "moc_gate_pass")},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    {"prepare": prepare, "compare": compare}[sys.argv[1]]()
