"""plan axisymmetric-freestream-hoop-gauge §4.6 の 7 (事前登録、既定化の条件・非静止場) の判定。hp7.sh の後に AWS の case/45 で回す。
- 腕: キー 0 = run_0482_hp7_k0、キー 1 = run_0483_hp7_k1 (FP64、同じバイナリ・格子・保存量、同じ停止規則)。
- 比較が成り立つ条件: 到達を保存した系列 (m9_watch.json) から計算し直して確かめる (出力の間隔 2500・末尾 2 万 step に 9 点以上)。満たさなければ判定不能。
- 不合格 (腕ごと): 非有限・非物理 (ρ ≤ 0・T ≤ 0・P ≤ 0、残した全スナップショット)・発散・上限 (20 万 step) までに水準に入らない。
- 既定化を保留する条件: 到達時の θ_r(40/70/94) のキー 1 とキー 0 の差が 0.5 % を超える、または Q_w の差が 1 % を超える (分母はキー 0 の到達時の値の絶対値)。
  超えたら原因を見てから決める (旧の解との一致は求めない)。
- 記録: 軸の上のマッハ数の分布の差、到達の step、check_convergence と check_quasisteady (系列、θ_r は drift/osc 1e-4、Q_w は 2e-4、--tail 1.0)。
  停止規則は水準への到達であって収束の判定ではないので、check_convergence の PASS と対象量の STEADY が無ければ「定常解の比較」とは呼ばない。
結果は標準出力と _band_ab/cold_pair/hp7_judge.json。"""
import json, subprocess, sys, h5py, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
ARMS = {"k0": HERE / "run_0482_hp7_k0", "k1": HERE / "run_0483_hp7_k1"}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94")
NJ, RT = 121, 0.076807
OUT = {"arms": {}, "fail": [], "hold": [], "undecidable": []}
def drift(rows, key):
    s_ = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000: return None
    j = int(np.searchsorted(s_, s_[-1] - 20000)); return float(100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j]))
def reach(rows):   # m9_watch の --phase line と同じ水準 (v4_judge.py と同じ式)
    run_ok = 0
    for n in range(len(rows)):
        r = rows[n]
        if r.get("nonfinite") or r.get("nonfinite_all"): return ("DIV", r["step"])
        dr = [drift(rows[: n + 1], k) for k in KEYS]
        win = sum(1 for x in rows[: n + 1] if r["step"] - 20000 <= x["step"] <= r["step"])
        ok = r.get("deficit") is not None and None not in dr and abs(r["deficit"]) <= 0.1 and all(abs(d) <= 0.05 for d in dr) and win >= 9
        run_ok = run_ok + 1 if ok else 0
        if run_ok >= 2: return ("REACH", r["step"])
    return None
def axis_mach(run, n):
    with h5py.File(run / f"res_{n}.h5", "r") as h:
        V = h["VALUE"]; X = np.asarray(h["MESH/COORD"][:], np.float64).reshape(-1, 3)
        M = np.hypot(np.asarray(V["Ux"][:], np.float64), np.asarray(V["Uy"][:], np.float64)) / np.asarray(V["sonic"][:], np.float64)
    j = np.arange(M.size) % NJ
    return X[j == 0, 0] / RT, M[j == 0]
for name, run in ARMS.items():
    st = json.loads((run / "m9_watch.json").read_text())
    rows = sorted([r for r in st["rows"] if r["step"] > 0], key=lambda r: r["step"])
    steps = [r["step"] for r in rows]
    a = {"status": st["status"], "reach_step_watch": st.get("reach_step"), "n_rows": len(rows)}
    if steps != [2500 * (i + 1) for i in range(len(steps))]: OUT["undecidable"].append(f"{name}: 出力の間隔が 2500 でない")
    d = reach(rows); a["reach_recomputed"] = d
    # 非物理: 残した全スナップショット
    unphys = []
    for p in sorted(run.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1])):
        with h5py.File(p, "r") as h:
            V = h["VALUE"]; bad = {k: int(np.count_nonzero(~(np.asarray(V[k][:], np.float64) > 0))) for k in ("ro", "T", "P")}
            nf = sum(int(np.count_nonzero(~np.isfinite(np.asarray(V[k][:])))) for k in V.keys())
        if any(bad.values()) or nf: unphys.append(f"{p.name}: ≤0 {bad}・非有限 {nf}")
    a["unphysical"] = unphys
    if unphys: OUT["fail"].append(f"{name}: 非物理・非有限 {unphys[:3]}")
    if d is None: OUT["fail"].append(f"{name}: 上限までに水準に入らない (最後の出力 {steps[-1] if steps else None})")
    elif d[0] != "REACH": OUT["fail"].append(f"{name}: 発散 ({d})")
    else:
        n = d[1]; a["reach"] = n
        if st.get("reach_step") != n: OUT["undecidable"].append(f"{name}: 見張りの到達 {st.get('reach_step')} と計算し直した到達 {n} が違う")
        R = next(r for r in rows if r["step"] == n); a["at_reach"] = {k: R[k] for k in KEYS + ("Q_w", "deficit")}
        a["window_points"] = sum(1 for x in rows if n - 20000 <= x["step"] <= n)
        a["drift_at_reach"] = [drift([x for x in rows if x["step"] <= n], k) for k in KEYS]
    OUT["arms"][name] = a
    print(f"== {name}: 見張り {a['status']} 到達 {a['reach_step_watch']}、計算し直し {d}、出力 {a['n_rows']} 点、非物理 {len(unphys)}")
k0, k1 = OUT["arms"]["k0"], OUT["arms"]["k1"]
if OUT["undecidable"]:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"])
elif OUT["fail"]:
    verdict = "不合格: " + "; ".join(OUT["fail"])
else:
    for k in KEYS + ("Q_w",):
        dv = 100 * (k1["at_reach"][k] - k0["at_reach"][k]) / abs(k0["at_reach"][k]); lim = 1.0 if k == "Q_w" else 0.5
        print(f"  {k:11s} キー 0 {k0['at_reach'][k]:.6g}  キー 1 {k1['at_reach'][k]:.6g}  差 {dv:+.4f} %  (保留の閾値 {lim} %)")
        OUT.setdefault("diff_pct", {})[k] = dv
        if abs(dv) > lim: OUT["hold"].append(f"{k} {dv:+.4f} %")
    verdict = "既定化の条件を満たす (到達時の差が保留の閾値以下)" if not OUT["hold"] else "既定化を保留 (原因を見てから決める): " + "; ".join(OUT["hold"])
print(f"== VERDICT §4.6 の 7: {verdict}")
OUT["verdict"] = verdict
print("== 記録: 軸の上のマッハ数 (到達時、キー 1 − キー 0)")
if k0.get("reach") and k1.get("reach"):
    x0, m0 = axis_mach(ARMS["k0"], k0["reach"]); x1, m1 = axis_mach(ARMS["k1"], k1["reach"])
    if not np.array_equal(x0, x1): print("  軸の節点の座標が一致しない (記録を省く)")
    else:
        o = np.argsort(x0); x0, m0, m1 = x0[o], m0[o], m1[o]
        rel = (m1 - m0) / np.maximum(np.abs(m0), 1e-3)
        sup = x0 >= 0
        print(f"  全域: 最大 |ΔM| {np.max(np.abs(m1 - m0)):.4e} (x/r_t {x0[np.argmax(np.abs(m1 - m0))]:.2f})、超音速部 (x ≥ 0) の最大 |ΔM/M| {np.max(np.abs(rel[sup])):.4e}")
        tab = []
        for xr in (-5, -2, -1, 0, 1, 2, 5, 10, 20, 40, 70, 94):
            i = int(np.argmin(np.abs(x0 - xr))); tab.append((float(x0[i]), float(m0[i]), float(m1[i])))
            print(f"    x/r_t {x0[i]:7.2f}: M キー 0 {m0[i]:.5f}  キー 1 {m1[i]:.5f}  差 {m1[i] - m0[i]:+.2e}")
        OUT["axis_mach"] = {"max_abs_dM": float(np.max(np.abs(m1 - m0))), "max_rel_sup": float(np.max(np.abs(rel[sup]))), "table": tab}
print(f"   到達の step: キー 0 {k0.get('reach_recomputed')}  キー 1 {k1.get('reach_recomputed')}")
print("== 記録: 判定ツール")
for name, run in ARMS.items():
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    OUT["arms"][name]["check_convergence"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
    print(f"  {name} check_convergence: {OUT['arms'][name]['check_convergence']}")
    a = OUT["arms"][name]
    if a.get("reach"):
        n = a["reach"]; st = json.loads((run / "m9_watch.json").read_text())
        rows = [x for x in st["rows"] if n - 20000 <= x["step"] <= n]
        sc = run / "hp7_series_tail.csv"
        with open(sc, "w") as fh:
            fh.write("step," + ",".join(KEYS + ("Q_w",)) + "\n")
            for x in sorted(rows, key=lambda x: x["step"]): fh.write(f"{x['step']}," + ",".join(f"{x[k]:.10g}" for k in KEYS + ("Q_w",)) + "\n")
        for cols, dr in ((",".join(KEYS), "0.0001"), ("Q_w", "0.0002")):
            p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", cols,
                                "--tail", "1.0", "--drift", dr, "--osc", dr], capture_output=True, text=True)
            ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
            a[f"check_quasisteady_{cols}"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
            print(f"  {name} check_quasisteady ({cols}): {a[f'check_quasisteady_{cols}']}")
o = HERE / "_band_ab/cold_pair/hp7_judge.json"; o.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", o)
