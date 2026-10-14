"""plan tooling-nozzle-core-grid §4.15 (codex 諮問 2026-10-11 の指摘): 残っている末尾の出力で、登録した指標を集計し直す。
腕: G0′ = run_0483_hp7_k1 (別セッション、読むだけ。到達までの末尾)、G1・G1x のキー 1・0 (run_0603〜0606、20 万 step までの末尾)。
- 指標 (§4.5): θ_r・δ_loc = PCHIP でつないだ x の窓 [40,50]・[65,75]・[84,94] の平均、Q_w、出口の主流 M の線平均、軸と η 0.1 の試験部の平均 100(M/6−1)。
- u = 窓の中の最後の値からのずれの最大 (相対 %、M の量は %pt)。定常性は check_quasisteady.py の系列モード (θ・δ は drift/osc 1e-4、Q_w は 2e-4、--tail 1.0)。
- 収束は各 run の CONVERGENCE_VERDICT.txt (無ければ check_convergence.py を回す)。
これは「到達後 4 万 step の窓」(§4.4) ではなく、残っている末尾 (約 2 万 step) での暫定の集計である。
usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): python3 cg_judge_tail.py → _band_ab/core_grid/judge_tail.json"""
import json, re, subprocess, sys
from pathlib import Path
import h5py
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "design"))
import cold_xcheck as XC
from forge_design.report import nozzle_report as NR
TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
ARMS = {"G0p_k1": Path.home() / "forge-wallfit/case/45.isobutane_m6_d155/run_0483_hp7_k1",
        "G1_k1": HERE / "run_0603_cg_g1_k1", "G1_k0": HERE / "run_0604_cg_g1_k0", "G1x_k1": HERE / "run_0605_cg_g1x_k1", "G1x_k0": HERE / "run_0606_cg_g1x_k0"}
WINDOWS = ((40.0, 50.0), (65.0, 75.0), (84.0, 94.0))
yb_x, yb = XC.common_yb()
def snap(run, res):
    with h5py.File(run / res, "r") as h:
        xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
        v = [np.asarray(h["VALUE/" + k][:], float) for k in ("ro", "Ux", "Uy", "T", "k")]
    ni, nj, S = XC.mesh_info(len(v[0]), run)
    rec = {}
    for prof in ("pchip", "linear"):
        o = XC.reduce_fields(xy, *v, False, ni, nj, S, yb_x, yb, profile=prof)
        for lo, hi in WINDOWS:
            w = (o["x"] >= lo) & (o["x"] <= hi)
            for q in ("theta_r", "delta_loc"):
                rec[f"{q}_w{int(lo)}_{int(hi)}{'' if prof == 'pchip' else '_lin'}"] = float(np.trapezoid(o[q][w], o["x"][w]) / (o["x"][w][-1] - o["x"][w][0]))
        if prof == "pchip": rec["Q_w"] = float(o["Q_w"])
    F = NR.load_field(run, res); info = F["info"]; Md = float(info.get("Md", 6.0)); xE = float(info["x_E"]); xF = float(F["X"][-1, 0])
    eta_last = F["R"][-1] / F["R"][-1, -1]; rec["exit_M_line"] = NR.exit_core_line_mean(eta_last, F["V"]["M"][-1])
    xq = np.linspace(float(F["X"][0, 0]), xF, 2401); w = (xq >= xE + 2) & (xq <= xF - 1)
    for et in (0.0, 0.1):
        rec[f"eta{et}_mean_pct"] = float((100 * (NR.eta_line(F, "M", et, xq) / Md - 1))[w].mean())
    return rec
out = {}
for name, run in ARMS.items():
    st = json.loads((run / "m9_watch.json").read_text()); end = st.get("reach_step") or st["rows"][-1]["step"]
    res = sorted((int(m.group(1)), f.name) for f in run.iterdir() if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)) and int(m.group(1)) > 0)
    res = [(n, f) for n, f in res if end - 20000 <= n <= end]
    rows = []
    for n, f in res:
        r = snap(run, f); r["step"] = n; rows.append(r); print(f"[{name}] {f}", flush=True)
    keys = [k for k in rows[0] if k != "step"]
    a = {"steps": [r["step"] for r in rows], "n": len(rows)}
    for k in keys:
        v = np.array([r[k] for r in rows]); last = v[-1]
        if k.endswith("_pct"):
            a[k] = {"mean": float(v.mean()), "last": float(last), "u_pctpt": float(np.abs(v - last).max())}
        else:
            a[k] = {"mean": float(v.mean()), "last": float(last), "u_pct": float(100 * np.abs(v - last).max() / abs(last))}
    sc = HERE / "_band_ab" / "core_grid" / f"tail_series_{name}.csv"
    with open(sc, "w") as fh:
        fh.write("step," + ",".join(keys) + "\n")
        for r in rows: fh.write(f"{r['step']}," + ",".join(f"{r[k]:.10g}" for k in keys) + "\n")
    qs = {}
    for cols, dr in ((",".join(k for k in keys if k.startswith(("theta_r_w", "delta_loc_w")) and not k.endswith("_lin")), "0.0001"), ("Q_w", "0.0002")):
        p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", cols,
                            "--tail", "1.0", "--drift", dr, "--osc", dr], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        qs[cols.split(",")[0] + ("…" if "," in cols else "")] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
    a["check_quasisteady"] = qs
    cv = run / "CONVERGENCE_VERDICT.txt"
    a["check_convergence"] = ([l for l in cv.read_text().splitlines() if "VERDICT" in l or "OVERALL" in l] or ["(行なし)"])[-1] if cv.is_file() else "(ファイルなし)"
    out[name] = a
def rel(a, b, k):
    if k.endswith("_pct"): return out[a][k]["mean"] - out[b][k]["mean"]
    return 100 * (out[a][k]["mean"] / out[b][k]["mean"] - 1)
keys = [k for k in out["G1_k1"] if isinstance(out["G1_k1"][k], dict) and "mean" in out["G1_k1"][k]]
out["grid_diff_vs_G0p"] = {a: {k: rel(a, "G0p_k1", k) for k in keys} for a in ("G1_k1", "G1x_k1")}
out["key_diff_k1_vs_k0"] = {g: {k: rel(g + "_k1", g + "_k0", k) for k in keys} for g in ("G1", "G1x")}
(HERE / "_band_ab" / "core_grid" / "judge_tail.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for a in ("G1_k1", "G1x_k1"):
    print(a, "対 G0′:", {k: round(v, 3) for k, v in out["grid_diff_vs_G0p"][a].items()})
for g in ("G1", "G1x"):
    print(g, "キー 1 − 0:", {k: round(v, 3) for k, v in out["key_diff_k1_vs_k0"][g].items()})
for n in ARMS:
    print(n, "u:", {k: round(out[n][k].get("u_pct", out[n][k].get("u_pctpt")), 3) for k in keys}, out[n]["check_quasisteady"], out[n]["check_convergence"])
