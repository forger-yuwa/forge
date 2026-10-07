"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U0 の `solve_rt` の記録 (CFD 0 step、forge は起動しない)。

`solve_rt` の **NS 後の経路** (prev_run あり、C2 方式で使う、変えない) と **CFD 前の経路** (prev_run なし、共通の δ_r 経路に揃える
ので変わる) を、DESIGN_DIR の design/ で回して JSON に書く。変更前 (HEAD の写し) と変更後の 2 回回し、upoly_verify.py u0 で比べる。

usage: DESIGN_DIR=<design/> python upoly_solve_rt_record.py PROBLEM.yaml PREV_RUN OUT.json [R_EXIT_M]
"""
import json
import os
import sys
import time
from pathlib import Path

C = Path(__file__).resolve().parent
DESIGN = Path(os.environ.get("DESIGN_DIR", str(C.parents[1] / "design"))).resolve()
sys.path.insert(0, str(DESIGN))


def main(argv):
    prob, prev, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    R = float(argv[3]) if len(argv) > 3 else 0.775
    from forge_design.feedback.deltastar_loop import solve_rt
    rec = {"design_dir": str(DESIGN), "problem": str(prob), "prev_run": str(prev), "R_exit_m": R}
    t0 = time.time()
    rec["ns_after"] = solve_rt(prob, R, prev_run=prev)
    rec["ns_after_sec"] = time.time() - t0
    t0 = time.time()
    try:
        rec["cfd_before"] = solve_rt(prob, R)
    except Exception as e:  # noqa: BLE001  (変更後は反復の上限で SizingNotConverged になりうる — そのまま記録する)
        rec["cfd_before"] = {"error": f"{type(e).__name__}: {e}", "history": getattr(e, "history", None)}
    rec["cfd_before_sec"] = time.time() - t0
    out.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    print(json.dumps({k: (v.get("r_t_m") if isinstance(v, dict) else v) for k, v in rec.items()}, indent=1, default=str))


if __name__ == "__main__":
    main(sys.argv[1:])
