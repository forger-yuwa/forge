"""plan tooling-nozzle-isothermal-wall-chain §5.1 #32・§6「V-tw370 の事前登録」(b5df7416、訂正 b6aeb263) の問題 YAML 3 本の設定の差の検査。
読み込んだ辞書 (yaml.safe_load) の差のパスで、登録の「3 本の設定の差」を確かめる (name は設定に数えない):
  R0 (g0_c1) と C1 (g0c_c1)        : 格子だけ (mesh.ni・mesh.x_respace・mesh.ar_max)
  C1 と C50 (g0c_c50)              : evaluate.main_line.cap だけ (1 → 50)
  R0 と tw300 (G0、V-c45)          : spec.wall_thermal.Tw (300 → 370) と evaluate だけ (cfl_main を外し、implicit_relax・main_line・nStepOuter・
                                     outStepInterval) で、evaluate は C1 の evaluate と同じ
  C1 と tw300_g0c_chain (V-g0c-chain): spec.wall_thermal.Tw と evaluate.main_line.cap (50 → 1) だけ
  C50 と tw300_g0c_chain            : spec.wall_thermal.Tw だけ
--design <作業ツリー/design> を付けると、forge_design の load_problem・_main_line_spec を 3 本に通し、解決した main_line (cap・itj・cfl・relax・
基準値・出力の量)・壁の熱境界条件 (isothermal 370 K)・bcond の壁の行 (Ts: 370.0) も確かめる。
usage: python3 tw370_yaml_check.py [--design DIR] [--json OUT]   (合格で rc 0、不成立で rc 1)
"""
import argparse
import json
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
Y = {"R0": "problem_d155_ns_prod_coldmesh_tw370_g0_c1.yaml", "C1": "problem_d155_ns_prod_coldmesh_tw370_g0c_c1.yaml",
     "C50": "problem_d155_ns_prod_coldmesh_tw370_g0c_c50.yaml",
     "tw300": "problem_d155_ns_prod_coldmesh_tw300.yaml", "tw300_g0c_chain": "problem_d155_ns_prod_coldmesh_tw300_g0c_chain.yaml"}
MISSING = "<無し>"
ML_WANT = {"line": "dir", "itj": 5, "cfl": 4.0, "relax": 0.7,
           "limiter_ref": {"limiterRoRef": 0.8739869154, "limiterPRef": 387491.0573, "limiterARef": 359.7712185},
           "extra_fields": ["res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega", "omg_prod", "omg_dest", "omg_cross",
                            "omg_trans", "omg_axisym", "dt_local"]}
CAP = {"R0": 1.0, "C1": 1.0, "C50": 50.0}
TW = 370.0


def diff_paths(a, b, path=()):
    """2 つの辞書の差 {パス: (a の値, b の値)}。片側にしか無いキーは MISSING (型が違うのも差)。"""
    out = {}
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            out.update(diff_paths(a.get(k, MISSING), b.get(k, MISSING), path + (k,)))
        return out
    if type(a) is not type(b) or a != b:
        out[path] = (a, b)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    D = {k: yaml.safe_load((HERE / v).read_text()) for k, v in Y.items()}
    fails, rec = [], {}

    def pair(x, y, allowed, want=None):
        d = {p: v for p, v in diff_paths(D[x], D[y]).items() if p != ("name",)}
        keys = {".".join(map(str, p)) for p in d}
        rec[f"{x} vs {y}"] = {".".join(map(str, p)): [str(v[0])[:80], str(v[1])[:80]] for p, v in d.items()}
        extra = sorted(keys - set(allowed))
        lack = sorted(set(allowed) - keys)
        if extra or lack:
            fails.append(f"{x} と {y} の差が登録と違う: 余分 {extra}、足りない {lack}")
        for p, (va, vb) in (want or {}).items():
            got = d.get(tuple(p.split(".")))
            if got != (va, vb):
                fails.append(f"{x} と {y} の {p} = {got} ({(va, vb)} であること)")

    pair("R0", "C1", ["mesh.ni", "mesh.x_respace", "mesh.ar_max"],
         {"mesh.ni": (4719, 4467), "mesh.ar_max": (5000, 10000),
          "mesh.x_respace": (MISSING, [{"x0": -8.0, "x1": -2.0, "blend": 0.5, "factor": 2.0, "ni_ref": 4719}])})
    pair("C1", "C50", ["evaluate.main_line.cap"], {"evaluate.main_line.cap": (1, 50)})
    pair("R0", "tw300", ["spec.wall_thermal.Tw", "evaluate.nStepOuter", "evaluate.outStepInterval", "evaluate.cfl_main",
                         "evaluate.implicit_relax", "evaluate.main_line"],
         {"spec.wall_thermal.Tw": (370.0, 300.0), "evaluate.cfl_main": (MISSING, 1.0)})
    pair("C1", "tw300_g0c_chain", ["spec.wall_thermal.Tw", "evaluate.main_line.cap"],
         {"spec.wall_thermal.Tw": (370.0, 300.0), "evaluate.main_line.cap": (1, 50)})
    pair("C50", "tw300_g0c_chain", ["spec.wall_thermal.Tw"], {"spec.wall_thermal.Tw": (370.0, 300.0)})
    if D["R0"]["evaluate"] != D["C1"]["evaluate"]:
        fails.append("R0 と C1 の evaluate が違う")
    rec["names"] = {k: D[k]["name"] for k in ("R0", "C1", "C50")}
    if a.design:
        sys.path.insert(0, a.design)
        from forge_design.evaluate import runner_axismach as RA
        from forge_design.evaluate.runner_wt import _bcond
        des = {}
        for k in ("R0", "C1", "C50"):
            p = RA.load_problem(HERE / Y[k])
            ml = RA._main_line_spec(p)
            wt = p.wall_thermal
            bc = _bcond(p, euler=False)
            wall = [l for l in bc.splitlines() if l.startswith("wall:")]
            des[k] = {"main_line": ml, "wall_thermal": wt, "bcond_wall": wall, "nStepOuter": p.evaluate["nStepOuter"],
                      "outStepInterval": p.evaluate["outStepInterval"], "implicit_relax": p.evaluate.get("implicit_relax"),
                      "cfl_main": p.evaluate.get("cfl_main", MISSING)}
            bad = [n for n, v in ML_WANT.items() if ml.get(n) != v]
            if ml.get("cap") != CAP[k] or bad or ml.get("prestage_relax") != 0.7:
                fails.append(f"{k}: _main_line_spec が登録と違う (cap {ml.get('cap')}、違うキー {bad}、prestage_relax {ml.get('prestage_relax')})")
            if wt != {"mode": "isothermal", "Tw": TW}:
                fails.append(f"{k}: wall_thermal {wt}")
            if len(wall) != 1 or "kind: wall_isothermal" not in wall[0] or "Ts: 370.0}" not in wall[0]:
                fails.append(f"{k}: bcond の壁の行 {wall}")
            if (p.evaluate["nStepOuter"], p.evaluate["outStepInterval"], p.evaluate.get("implicit_relax"), "cfl_main" in p.evaluate) != (200000, 2500, 0.7, False):
                fails.append(f"{k}: evaluate の step 数・出力の間隔・implicit_relax・cfl_main が登録と違う")
        rec["design"] = des
    out = {"ok": not fails, "failures": fails, "diffs": rec}
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print("VERDICT:", "PASS (3 本の設定の差は登録どおり)" if not fails else "FAIL")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
