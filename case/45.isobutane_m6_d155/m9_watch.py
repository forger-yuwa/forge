"""§5.1 #9 の見張り (plan time_integration-line-implicit-speed §6.9、2026-10-10)。
走っている run の res_<N>.h5 (N > 0、書き終わって 60 s 以上) を 1 つずつ `cold_series.snapshot` で 1 行にして系列 (`_band_ab/cold_pair/series_<run>.json`) に足し、
  - 非有限 → DIVERGED を書いて forge を止める
  - 水準 (|欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step、窓はこの run の中だけ) → REACHED を書いて forge を止め、それより後の出力を消す
  - 5 万 step ごとの出力と最新の出力以外の res_<N>.h5 (と同じ N の res_wall・res_outlet) を消す (ディスクのため)
forge が終わって (自分の run の forge が無い) 未処理の出力が無くなったら抜ける。止めるのは cwd がこの run の forge だけ (pgrep -x forge + /proc/PID/cwd)。
usage: python3 m9_watch.py <run> → <run>/m9_watch.json (状態)、標準出力にログ"""
import json, os, re, signal, subprocess, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
run = HERE / sys.argv[1]
OUT = XC.OUTD / f"series_{run.name}.json"
state = {"run": run.name, "status": "running", "reach_step": None, "processed": []}
rows = json.loads(OUT.read_text())["rows"] if OUT.exists() else []
done = {r["step"] for r in rows}
yb_x, yb = XC.common_yb()

def forge_pids():
    out = subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split()
    pids = []
    for p in out:
        try:
            if os.path.realpath(f"/proc/{p}/cwd") == str(run.resolve()): pids.append(int(p))
        except OSError:
            pass
    return pids

def stop(reason):
    for pid in forge_pids():
        try: os.kill(pid, signal.SIGTERM)
        except OSError: pass
    print(f"[m9_watch] {run.name}: {reason} — forge を止めた", flush=True)

def drift(rows, key):
    st = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if st[-1] - st[0] < 20000: return None
    j = int(np.searchsorted(st, st[-1] - 20000))
    return float(100 * (v[-1] / v[j] - 1) * 20000 / (st[-1] - st[j]))

def save():
    OUT.write_text(json.dumps({"run": run.name, "rows": rows, "level_rule": "|欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step (窓はこの run の中)"}, indent=1, ensure_ascii=False))
    (run / "m9_watch.json").write_text(json.dumps(state, indent=1, ensure_ascii=False))

def res_steps():
    return sorted(int(m.group(1)) for f in run.glob("res_*.h5") if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)))

while True:
    alive = bool(forge_pids())
    new = [n for n in res_steps() if n > 0 and n not in done and (time.time() - (run / f"res_{n}.h5").stat().st_mtime > 60 or not alive)]
    for n in new:
        try:
            rec = CS.snapshot(run, n, yb_x, yb)
        except Exception as e:  # 書きかけ等
            print(f"[m9_watch] {run.name} {n}: 読めない ({e}) — あとで", flush=True); break
        rows.append(rec); rows.sort(key=lambda r: r["step"]); done.add(n); state["processed"].append(n)
        dr = [drift(rows, k) for k in ("theta_r_40", "theta_r_70", "theta_r_94")]
        print(f"[m9_watch] {run.name} {n}: 欠損 {rec['deficit']:.4f}、θ_r(70) {rec['theta_r_70']:.6f}、Q_w {rec['Q_w']/1e6:.4f} MW、ドリフト {dr}、非有限 {rec['nonfinite']}", flush=True)
        if rec["nonfinite"]:
            state.update(status="DIVERGED", fail_step=n); save(); stop(f"非有限 (step {n})"); break
        if rec["deficit"] is not None and abs(rec["deficit"]) <= 0.1 and None not in dr and all(abs(x) <= 0.05 for x in dr):
            state.update(status="REACHED", reach_step=n, drift=dr); save(); stop(f"水準に到達 (step {n})")
            time.sleep(20)
            for m in res_steps():
                if m > n:
                    for f in run.glob(f"res*_{m}.h5"): f.unlink()
                    for f in run.glob(f"res*_{m}.xmf"): f.unlink()
            break
        # 消す: 5 万 step ごとと最新の処理済み以外
        latest = max(done)
        for m in sorted(done):
            if m != latest and m % 50000 != 0:
                for f in list(run.glob(f"res_{m}.h5")) + list(run.glob(f"res_*_{m}.h5")) + list(run.glob(f"res*_{m}.xmf")):
                    f.unlink()
        save()
    if state["status"] != "running":
        break
    if not alive and not [n for n in res_steps() if n > 0 and n not in done]:
        state["status"] = "ENDED"; save(); print(f"[m9_watch] {run.name}: forge が終わった (水準に届かず)", flush=True); break
    time.sleep(30)
save()
