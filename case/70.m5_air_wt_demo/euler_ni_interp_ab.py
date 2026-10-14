"""case/70: Euler の ni が δ* の抽出 (δ_E・流量の比) に与える影響を測る (plan tooling-nozzle-mesh-size-guard §6 V1a・V1b)。

- stride (V1a): Euler の解は同じまま、抽出が使う列だけを間引く。測れるのは `euler_q_at` の列間の線形補間の誤差だけ。
- runs (V1b): ni の違う Euler の run (同じ段取りで回したもの) を同じ NS に対して抽出し、最初の run との差を測る (Euler の解そのものの ni 依存)。
forge は回さない。

usage (case dir、design/.venv-opt):
  python euler_ni_interp_ab.py NS_RUN EULER_RUN X_E OUT_DIR [STRIDE ...]          # V1a、既定 STRIDE 1 2 3
  python euler_ni_interp_ab.py --runs NS_RUN X_E OUT_DIR EULER_RUN_REF EULER_RUN ...  # V1b
"""
import json
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next  # noqa: E402

ORIG = DS._load_structured


def run_stride(ns_run, eu_run, stride, out):
    eu = eu_run.resolve()

    def load(run_dir):
        d = ORIG(run_dir)
        if Path(run_dir).resolve() != eu or stride == 1:
            return d
        ni = d["x"].shape[0]
        keep = np.unique(np.r_[np.arange(0, ni, stride), ni - 1])      # 端の列は必ず残す
        return dict(d, x=d["x"][keep], r=d["r"][keep], q=d["q"][keep])

    DS._load_structured = load
    try:
        od = out / f"s{stride}"
        summ = extract_and_merge(ns_run, eu_run, band_select="edge", out_dir=od)
    finally:
        DS._load_structured = ORIG
    nx = read_delta_r_next(od / "delta_r_next.csv")
    return np.asarray(nx["x_rt"]), np.asarray(nx["delta_E"]), summ.get("massflow")


def run_runs(ns_run, eu_runs, x_e, out):
    """V1b: 各 Euler の run で抽出し、最初の run との差 (使用域 [0, x_e] の max|Δδ_E|・流量の比の差) を出す。"""
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    for eu in eu_runs:
        od = out / eu.name
        summ = extract_and_merge(ns_run, eu, band_select="edge", out_dir=od)
        nx = read_delta_r_next(od / "delta_r_next.csv")
        res[eu.name] = (np.asarray(nx["x_rt"]), np.asarray(nx["delta_E"]), summ.get("massflow"))
    ref = eu_runs[0].name
    x0, d0, m0 = res[ref]
    use = (x0 >= 0.0) & (x0 <= x_e)
    summ = {"ns_run": ns_run.name, "x_e": x_e, "ref": ref, "rows": []}
    for eu in eu_runs:
        x, d, m = res[eu.name]
        assert np.array_equal(x, x0)
        diff = np.abs(d - d0)[use]
        k = int(np.argmax(diff))
        summ["rows"].append({"euler_run": eu.name, "max_abs_dE_use_rt": float(diff.max()), "x_at_max_rt": float(x0[use][k]),
                             "mdot_ratio": m["mdot_ratio"], "d_mdot_ratio": m["mdot_ratio"] - m0["mdot_ratio"], "massflow": m})
    (out / "euler_ni_runs_ab.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False, default=float))
    print(json.dumps(summ, indent=1, ensure_ascii=False, default=float))


def main():
    if sys.argv[1] == "--runs":
        run_runs(C / sys.argv[2], [C / a for a in sys.argv[5:]], float(sys.argv[3]), Path(sys.argv[4]))
        return
    ns_run, eu_run, x_e, out = C / sys.argv[1], C / sys.argv[2], float(sys.argv[3]), Path(sys.argv[4])
    strides = [int(s) for s in sys.argv[5:]] or [1, 2, 3]
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    for s in strides:
        res[s] = run_stride(ns_run, eu_run, s, out)
    x0, d0, m0 = res[strides[0]]
    use = (x0 >= 0.0) & (x0 <= x_e)
    summ = {"ns_run": str(ns_run.name), "euler_run": str(eu_run.name), "x_e": x_e, "ref_stride": strides[0], "rows": []}
    for s in strides:
        x, d, m = res[s]
        assert np.array_equal(x, x0)
        diff = np.abs(d - d0)[use]
        k = int(np.argmax(diff))
        summ["rows"].append({"stride": s, "max_abs_dE_use_rt": float(diff.max()), "x_at_max_rt": float(x0[use][k]),
                             "massflow": m})
    (out / "euler_ni_interp_ab.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False, default=float))
    print(json.dumps(summ, indent=1, ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()
