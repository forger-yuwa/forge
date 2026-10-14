"""case/70: 冷却壁の NS の遅いモードの切り分け (plan tooling-nozzle-cooled-wall-ns-settling §5.1 #1、forge 0 step、読むだけ)。

run の各スナップショットで、どの量がどの速さで動いているかを測る:
  - 出口コア M (η 0.05〜0.7 の平均)、無次元の流量 (x ∈ (−2.5, x_E) の中央値; qs_series.py と同じ定義)
  - スロートの断面 (x = 0 に最も近い列) の NS だけの排除厚さの近似: コア (r < 0.5 r_w) の平均の ρU を縁の値として
    δ_m = r_w − sqrt(2∫ρU r dr / ρU_core)  [r_t]
  - 壁の熱流束 qwall (壁ダンプ `res_wall_3_<step>.h5`) の区間ごとの平均 (直管・縮流部・スロート |x| < 1 r_t・拡大部の前半・後半) と
    総熱量の近似 Q = ∫ qwall 2π r ds (節点の台形)
連続した run (延長) は通算の step で並べる。増分の比 (10k ごと) と、等比で外挿した残りの変化を出す。

usage (case dir、design/.venv-opt):  python slow_mode_probe.py OUT.json RUN[:OFFSET] [RUN[:OFFSET] ...] [--every K]
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import _res_files, load_field  # noqa: E402

argv = list(sys.argv[1:])
every = 5
if "--every" in argv:
    k = argv.index("--every"); every = int(argv[k + 1]); del argv[k:k + 2]
args = argv
out = Path(args[0])
rows = []
for spec in args[1:]:
    name, off = (spec.split(":") + ["0"])[:2]
    run = (C / name).resolve()
    S = float(json.loads((run / "prepare_info.json").read_text())["scale_m"])
    for r in _res_files(run)[::every]:
        F = load_field(run, r)
        step = int("".join(ch for ch in Path(str(r)).stem.split("_")[-1] if ch.isdigit()))
        X, R, V = F["X"], F["R"], F["V"]
        q = V["ro"] * V["Ux"]
        eta = R[-1] / R[-1, -1]
        core = (eta >= 0.05) & (eta <= 0.7)
        xs = X[:, 0]
        xE = float(F["info"].get("x_E", 40.0))
        m = 2 * np.pi * np.trapezoid(q * R, R, axis=1)
        i0 = int(np.argmin(np.abs(xs)))
        rr, qq = R[i0], q[i0]
        rw = float(rr[-1])
        qc = float(qq[rr < 0.5 * rw].mean())
        reff = np.sqrt(2 * np.trapezoid(qq * rr, rr) / qc)
        row = {"run": name, "step": step + int(off), "exit_core_M": float(V["M"][-1][core].mean()),
               "mdot_rt": float(np.median(m[(xs > -2.5) & (xs < xE)])), "delta_m_throat_rt": float(rw - reff)}
        wf = run / f"res_wall_3_{step}.h5"
        if wf.exists():
            with h5py.File(wf) as f:
                xyz = f["MESH/COORD"][:].reshape(-1, 3)
                qw = f["VALUE/qwall"][:]
            xw, rwall = xyz[:, 0] / S, xyz[:, 1]
            o = np.argsort(xw)
            xw, rwall, qw = xw[o], rwall[o], qw[o]
            ds = np.hypot(np.diff(xw * S), np.diff(rwall))
            g = qw * 2 * np.pi * rwall
            row["Q_total_W"] = float(np.sum(0.5 * (g[1:] + g[:-1]) * ds))
            for tag, lo, hi in (("pipe", -1e9, -12.0), ("contraction", -12.0, -1.0), ("throat", -1.0, 1.0),
                                ("div_front", 1.0, 20.0), ("div_rear", 20.0, 1e9)):
                mm = (xw >= lo) & (xw < hi)
                row[f"qwall_mean_{tag}"] = float(qw[mm].mean()) if mm.any() else None
        rows.append(row)
        print(json.dumps(row), flush=True)
rows.sort(key=lambda d: d["step"])
out.write_text(json.dumps(rows, indent=1))
