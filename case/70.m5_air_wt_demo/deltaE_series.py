"""case/70: NS の末尾スナップショットごとに同じ修正済みの抽出 (等価半径の探索を積分域に限る; ab_root_zone.py の zone の腕と同じ差し替え) と
平滑化で δ_E(x) を作る (plan tooling-nozzle-deltastar-band-root-zone、codex plan レビュー 2026-10-13 M4: 再設計の入力 δ_E の準定常確認)。
Euler 参照は固定。forge は回さない。

usage (case dir、design/.venv-opt):
  python deltaE_series.py one NS_RUN EULER_RUN STEP OUT_DIR [EDGE_CSV]   # 1 枚 (res_<STEP>.h5) を抽出して OUT_DIR/<STEP>/ に書く
                                                                         # EDGE_CSV: 帯位置 y_b(x) をその表に固定 (診断用、§5.1 #5b の腕 B)
  python deltaE_series.py summary OUT_DIR X_E                     # 各枚の δ_E を集め、使用域 [0, X_E] の変化と check_quasisteady の判定を出す
"""
import csv
import inspect
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next  # noqa: E402

STATIONS = (0.0, 1.0, 3.0, 10.0, 20.0, 30.0, 48.0)


def patch_zone():
    """band_local_deficit の等価半径の探索を壁〜y_b の区間に限る (plan §4.1 と同じ; 実装前の測定用の差し替え)。"""
    src = inspect.getsource(DS.band_local_deficit)
    old = "            r_eff = float(np.interp(D, F_from_wall, r_desc))\n"
    new = ("            kz = int(np.searchsorted(rwN - r_desc, y_b, side='right'))\n"
           "            r_eff = float(np.interp(D, F_from_wall[:kz + 1], r_desc[:kz + 1]))\n")
    assert src.count(old) == 1
    ns = dict(DS.__dict__)
    exec(compile(src.replace(old, new), "patched", "exec"), ns)
    DS.band_local_deficit = ns["band_local_deficit"]
    ns["band_local_deficit"] = ns["band_local_deficit"]


def fix_band(edge_csv):
    """診断用 (plan §5.1 #5b の腕 B): edge の帯位置 y_b(x) を別の抽出の `delta_r_edge.csv` の値に固定する
    (梯子・縁の推定は計算するが y_b だけ差し替える; 流束のフィット・積分・平滑化は各時刻で行う)。"""
    tb = np.loadtxt(edge_csv, delimiter=",", skiprows=1)
    xf, ybf = tb[:, 0], tb[:, 5]
    orig = DS.edge_band_positions

    def fixed(x, *a, **k):
        out = orig(x, *a, **k)
        xx = np.asarray(x, dtype=float)
        assert xx.min() >= xf.min() - 1e-6 and xx.max() <= xf.max() + 1e-6, "固定する帯位置の表が抽出の範囲を覆わない"
        out = dict(out)
        out["y_b_reselected"] = out["y_b"]
        out["y_b"] = np.interp(xx, xf, ybf)
        return out

    DS.edge_band_positions = fixed


def one(ns_run, eu_run, step, out):
    ns_run, eu_run = (C / ns_run).resolve(), (C / eu_run).resolve()
    orig = DS._load_structured

    def load(run_dir):
        run_dir = Path(run_dir)
        if run_dir.resolve() != ns_run:
            return orig(run_dir)
        import h5py
        info = json.loads((run_dir / "prepare_info.json").read_text())
        ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
        S = float(info["scale_m"])
        with h5py.File(run_dir / "nozzle.h5") as f:
            nc = f["/MESH/COORD"][:].reshape(-1, 3)
        res = run_dir / f"res_{step}.h5"
        with h5py.File(res) as f:
            q = f["/VALUE/ro"][:] * f["/VALUE/Ux"][:]
        return dict(info=info, S=S, res=res.name, x=nc[:, 0].reshape(ni, nj) / S, r=nc[:, 1].reshape(ni, nj) / S,
                    q=q.reshape(ni, nj))

    DS._load_structured = load
    patch_zone()
    od = Path(out) / str(step)
    summ = extract_and_merge(ns_run, eu_run, band_select="edge", out_dir=od)
    (od / "massflow.json").write_text(json.dumps(summ.get("massflow"), default=float))


def summary(out, x_e):
    out = Path(out)
    steps = sorted(int(p.name) for p in out.iterdir() if p.is_dir() and p.name.isdigit())
    tabs = {}
    for s in steps:
        nx = read_delta_r_next(out / str(s) / "delta_r_next.csv")
        tabs[s] = (np.asarray(nx["x_rt"]), np.asarray(nx["delta_E"]))
    x0 = tabs[steps[-1]][0]
    use = (x0 >= 0.0) & (x0 <= x_e)
    rows, res = [], {"steps": steps}
    for s in steps:
        x, d = tabs[s]
        assert np.array_equal(x, x0)
        rows.append({"step": s, **{f"dE_{st:g}": float(np.interp(st, x, d)) for st in STATIONS}})
    fin = tabs[steps[-1]][1]
    res["max_abs_change_vs_last_use_domain_rt"] = {str(s): float(np.max(np.abs(tabs[s][1] - fin)[use])) for s in steps}
    res["max_abs_change_last_two_use_domain_rt"] = float(np.max(np.abs(tabs[steps[-1]][1] - tabs[steps[-2]][1])[use]))
    p = out / "deltaE_series.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    cols = ",".join(k for k in rows[0] if k != "step")
    q = subprocess.run([sys.executable, str(C.parents[1] / "solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(p),
                        "--series-cols", cols, "--tail", f"{min(1.0, 4.5 / len(rows)):.6f}"], capture_output=True, text=True)
    (out / "quasisteady_deltaE.txt").write_text(q.stdout + q.stderr)
    res["quasisteady"] = (q.stdout + q.stderr).strip().splitlines()[-12:]
    res["rows"] = rows
    (out / "deltaE_series.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps({k: res[k] for k in ("max_abs_change_vs_last_use_domain_rt", "max_abs_change_last_two_use_domain_rt", "quasisteady")},
                     indent=1, ensure_ascii=False))


if __name__ == "__main__":
    if sys.argv[1] == "one":
        if len(sys.argv) > 6:                       # 6 番目の引数: 帯位置を固定する delta_r_edge.csv (腕 B)
            fix_band(sys.argv[6])
        one(sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
    else:
        summary(sys.argv[2], float(sys.argv[3]))
