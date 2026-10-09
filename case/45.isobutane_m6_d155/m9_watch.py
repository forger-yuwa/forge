"""§5.1 #9 の見張り v2 (plan time_integration-line-implicit-speed §6.11、2026-10-10、codex plan-2 の指摘の反映)。
usage: python3 m9_watch.py <run> --phase line|point --budget <最大 step> → <run>/m9_watch.json (状態と系列をまとめて原子的に保存)
- 書き終わった res_<N>.h5 (N > 0、60 s 以上前に書かれたもの、または forge が終わった後) を 1 つずつ `cold_series.snapshot` で 1 行にする。
  VALUE の全ての浮動小数点の量の非有限も数える (cold_series は ρ・u・T・k だけ)。読めないものは 3 回まで再試行、それでも読めなければ DATA_ERROR。
- 終わりの判定 (step 0 は使わない): 水準 = |欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step (窓はこの run の中)。
  --phase point では加えて θ_r(40/70/94)・Q_w が参照 R (m9_ref.json) から相対 0.1 % 以内 (E2)。満たした最初の出力で REACHED。
- 状態: REACHED / CENSORED (forge が上限まで正常に回って届かず) / DIVERGED (非有限) / EXEC_ERROR (forge が上限より前に異常終了) / DATA_ERROR。
- REACHED・DIVERGED では cwd がこの run の forge に SIGTERM、180 s 待って残れば SIGKILL、終了を確かめてからそれより後の出力を消し、到達の出力の sha256 を記録する。
- 消す: 5 万 step ごと・最新 5 つ (準定常の判定の窓)・到達の出力以外の res_<N>.h5 と同じ N の res_* (状態を保存した後に消す)。
- 再開: 既存の m9_watch.json を読み、最終の状態なら何もしない。"""
import argparse, hashlib, json, os, re, signal, subprocess, sys, time
from pathlib import Path
import h5py
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--phase", choices=("line", "point", "line_e2"), required=True)   # line_e2 = ラインのまま E (水準 + E2) まで (§6.11 の追加の腕); ap.add_argument("--budget", type=int, required=True)
a = ap.parse_args()
run = HERE / a.run
SF = run / "m9_watch.json"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")
REF = None
if a.phase in ("point", "line_e2"):                    # E2 の参照は数値 4 つだけ (説明の文字列は使わない、codex plan-3 M2)
    _r = json.loads((HERE / "m9_ref.json").read_text())
    REF = {k: float(_r[k]) for k in KEYS}
    if not all(np.isfinite(v) and v > 0 for v in REF.values()): print("[m9_watch] 参照が有限・正でない — 止める"); sys.exit(2)
FINAL = ("REACHED", "CENSORED", "DIVERGED", "EXEC_ERROR", "DATA_ERROR")
st = json.loads(SF.read_text()) if SF.exists() else {"run": run.name, "phase": a.phase, "budget": a.budget, "status": "running", "rows": [], "tries": {}}
if st["status"] in FINAL:
    print(f"[m9_watch] {run.name}: 既に {st['status']}"); sys.exit(0)
if st.get("phase") != a.phase or st.get("budget") != a.budget:
    print(f"[m9_watch] {run.name}: 状態の phase/budget が違う — 止める"); sys.exit(2)
yb_x, yb = XC.common_yb()

def save():
    tmp = SF.with_suffix(".tmp"); tmp.write_text(json.dumps(st, indent=1, ensure_ascii=False)); os.replace(tmp, SF)
def pids():
    out = subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split(); r = []
    for p in out:
        try:
            if os.path.realpath(f"/proc/{p}/cwd") == str(run.resolve()): r.append(int(p))
        except OSError: pass
    return r
def stop():
    ps = pids()
    for p in ps:
        try: os.kill(p, signal.SIGTERM)
        except OSError: pass
    t0 = time.time()
    while pids() and time.time() - t0 < 180: time.sleep(5)
    for p in pids():
        try: os.kill(p, signal.SIGKILL)
        except OSError: pass
    t0 = time.time()
    while pids() and time.time() - t0 < 60: time.sleep(2)
    return not pids()
def steps():
    return sorted(int(m.group(1)) for f in run.glob("res_*.h5") if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)))
def drift(rows, key):
    s_ = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000: return None
    j = int(np.searchsorted(s_, s_[-1] - 20000)); return float(100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j]))
def nonfinite_all(n):
    c = 0
    with h5py.File(run / f"res_{n}.h5", "r") as h:
        for k, d in h["VALUE"].items():
            x = np.asarray(d[:])
            if x.dtype.kind == "f": c += int(np.count_nonzero(~np.isfinite(x)))
    return c
def finish(status, **kw):
    st.update(status=status, **kw); save(); print(f"[m9_watch] {run.name}: {status} {kw}", flush=True)
def cleanup(keep_extra=()):
    done = sorted(r["step"] for r in st["rows"])
    keep = set(done[-5:]) | {d for d in done if d % 50000 == 0} | set(keep_extra)
    for m in done:
        if m in keep: continue
        for f in list(run.glob(f"res_{m}.h5")) + list(run.glob(f"res_*_{m}.h5")) + list(run.glob(f"res*_{m}.xmf")): f.unlink()

def done_at(rows):
    """保存済みの系列の最初の終わりの出力 (水準 [+ E2]) を返す。非有限があればその step を負で返す。"""
    for n_ in range(len(rows)):
        r_ = rows[n_]
        if r_.get("nonfinite") or r_.get("nonfinite_all"): return -r_["step"]
        dr_ = [drift(rows[: n_ + 1], k) for k in KEYS[:3]]
        e2_ = {k: 100 * (r_[k] / REF[k] - 1) for k in KEYS} if REF else None
        if r_.get("deficit") is not None and None not in dr_ and abs(r_["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr_) \
           and (e2_ is None or all(abs(v) <= 0.1 for v in e2_.values())):
            return r_["step"]
    return None
def conclude(n):
    """終わりの出力 n で止める (何度呼んでも同じ結果): forge を止め、後の出力を消し、sha256 を記録して REACHED。"""
    ok = stop()
    if not ok: finish("EXEC_ERROR", note="forge を止められない", step=n); sys.exit(1)
    for m in steps():
        if m > n:
            for f in list(run.glob(f"res_{m}.h5")) + list(run.glob(f"res_*_{m}.h5")) + list(run.glob(f"res*_{m}.xmf")): f.unlink()
    if not (run / f"res_{n}.h5").exists(): finish("DATA_ERROR", note=f"到達の出力 res_{n}.h5 が無い", step=n); sys.exit(1)
    sha = hashlib.sha256((run / f"res_{n}.h5").read_bytes()).hexdigest()
    finish("REACHED", reach_step=n, reach_sha256=sha); cleanup(keep_extra=(n,)); sys.exit(0)
# 再開: 保存済みの系列で先に判定し直す (codex plan-3 M4)
_d = done_at(st["rows"])
if _d is not None:
    if _d < 0: stop(); finish("DIVERGED", fail_step=-_d); sys.exit(0)
    conclude(_d)

def main_loop():
  while True:
      alive = bool(pids())
      have = {r["step"] for r in st["rows"]}
      new = [n for n in steps() if n > 0 and n not in have and (not alive or time.time() - (run / f"res_{n}.h5").stat().st_mtime > 60)]
      progressed = False
      for n in new:
          try:
              rec = CS.snapshot(run, n, yb_x, yb); rec["nonfinite_all"] = nonfinite_all(n)
          except Exception as e:
              st["tries"][str(n)] = st["tries"].get(str(n), 0) + 1; save()
              print(f"[m9_watch] {run.name} {n}: 読めない ({e}) — {st['tries'][str(n)]} 回目", flush=True)
              if st["tries"][str(n)] >= 3: finish("DATA_ERROR", step=n); stop(); sys.exit(1)
              break
          st["rows"].append(rec); st["rows"].sort(key=lambda r: r["step"]); progressed = True
          rows = st["rows"]
          dr = [drift(rows, k) for k in ("theta_r_40", "theta_r_70", "theta_r_94")]
          e2 = {k: 100 * (rec[k] / REF[k] - 1) for k in KEYS} if REF else None
          rec["drift"] = dr; rec["e2"] = e2
          print(f"[m9_watch] {run.name} {n}: 欠損 {rec['deficit']:.4f}、ドリフト {dr}、E2 {e2}、非有限 {rec['nonfinite']}/{rec['nonfinite_all']}", flush=True)
          if rec["nonfinite"] or rec["nonfinite_all"]:
              save(); ok = stop(); finish("DIVERGED", fail_step=n, stopped=ok); sys.exit(0)
          lvl = rec["deficit"] is not None and None not in dr and abs(rec["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr)
          if lvl and (e2 is None or all(abs(v) <= 0.1 for v in e2.values())):
              save(); conclude(n)
          save(); cleanup()
      if not progressed and not alive and not [n for n in steps() if n > 0 and n not in {r["step"] for r in st["rows"]}]:
          rcf = run / "RUN_RC"; rc = rcf.read_text().strip() if rcf.exists() else "?"
          last = max((r["step"] for r in st["rows"]), default=0)
          if rc == "0" and last >= a.budget: finish("CENSORED", last_step=last)
          else: finish("EXEC_ERROR", rc=rc, last_step=last)
          sys.exit(0)
      time.sleep(30)

try:
    main_loop()
except SystemExit:
    raise
except Exception as e:                                   # 見張りの例外でも forge を止めて状態を残す (codex plan-3 M2)
    import traceback; traceback.print_exc()
    ok = stop(); finish("EXEC_ERROR", note=f"見張りの例外: {e!r}", stopped=ok); sys.exit(1)
