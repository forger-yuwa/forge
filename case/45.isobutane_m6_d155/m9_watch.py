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
def _emergency_stop(run_name, why):
    """初期化の失敗でも、cwd がその run の forge を止めて状態を残す (codex plan-4 M1)。"""
    import signal as _sg, subprocess as _sp, time as _t
    rd = os.path.realpath(str(HERE / run_name))
    def _p():
        r_ = []
        for q in _sp.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split():
            try:
                if os.path.realpath(f"/proc/{q}/cwd") == rd: r_.append(int(q))
            except OSError: pass
        return r_
    for sig_, w_ in ((_sg.SIGTERM, 180), (_sg.SIGKILL, 60)):
        for q in _p():
            try: os.kill(q, sig_)
            except OSError: pass
        t0_ = _t.time()
        while _p() and _t.time() - t0_ < w_: _t.sleep(3)
    try:
        sf_ = HERE / run_name / "m9_watch.json"
        st_ = json.loads(sf_.read_text()) if sf_.exists() else {"run": run_name, "rows": []}
        st_.update(status="EXEC_ERROR", note=f"見張りの初期化の失敗: {why}")
        sf_.write_text(json.dumps(st_, indent=1, ensure_ascii=False))
    except Exception: pass
    print(f"[m9_watch] {run_name}: 初期化の失敗 ({why}) — forge を止めた", flush=True); sys.exit(1)

ap = argparse.ArgumentParser(); ap.add_argument("run")
ap.add_argument("--phase", choices=("line", "point", "line_e2"), required=True)   # line_e2 = ラインのまま E (水準 + E2) まで (§6.11 の追加の腕 L5L)
ap.add_argument("--budget", type=int, required=True)
ap.add_argument("--consec", type=int, default=1, help="水準 (point と line_e2 は E2 も) を何出力続けて満たしたら到達とするか (既定 1 = 従来。§6.18 は 2、diagnostician 2026-10-10)")
ap.add_argument("--inherit", default=None, help="line_e2: 同じ構成の分岐元の run。その到達の出力までの系列を step − 到達 (≤ 0) で引き継ぐ (codex plan-4 M2)")
try:
    a = ap.parse_args()
except SystemExit:
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"): _emergency_stop(sys.argv[1], "引数の解釈")
    raise
try:
    run = HERE / a.run
    SF = run / "m9_watch.json"
    KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")
    REF = None
    if a.phase in ("point", "line_e2"):                    # E2 の参照は数値 4 つだけ (説明の文字列は使わない、codex plan-3 M2)
        _r = json.loads((HERE / "m9_ref.json").read_text())
        REF = {k: float(_r[k]) for k in KEYS}
        if not all(np.isfinite(v) and v > 0 for v in REF.values()): raise ValueError("参照が有限・正でない")
    if (a.phase == "line_e2") != (a.inherit is not None): raise ValueError("--inherit は line_e2 のときだけ、必ず付ける")
    FINAL = ("REACHED", "CENSORED", "DIVERGED", "EXEC_ERROR", "DATA_ERROR")
    if SF.exists():
        st = json.loads(SF.read_text())
    else:
        st = {"run": run.name, "phase": a.phase, "budget": a.budget, "consec": a.consec, "status": "running", "rows": [], "tries": {}}
        if a.inherit:
            src = json.loads((HERE / a.inherit / "m9_watch.json").read_text())
            if src.get("status") != "REACHED": raise ValueError(f"分岐元 {a.inherit} が REACHED でない")
            nb = int(src["reach_step"])
            st["inherit"] = {"run": a.inherit, "reach_step": nb, "reach_sha256": src.get("reach_sha256")}
            st["rows"] = [dict(r, step=r["step"] - nb, inherited=True) for r in src["rows"] if 0 < r["step"] <= nb]
    if a.consec < 1: raise ValueError("--consec は 1 以上")
    if st.get("phase") != a.phase or st.get("budget") != a.budget or st.get("consec", 1) != a.consec:   # 終端の状態でも先に照合する (codex plan-7 M3)
        raise ValueError("状態の phase/budget/consec が違う")
    if st["status"] in FINAL:
        print(f"[m9_watch] {run.name}: 既に {st['status']}"); sys.exit(0)
except SystemExit:
    raise
except Exception as e:
    _emergency_stop(a.run, repr(e))
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
    done = sorted(r["step"] for r in st["rows"] if r["step"] > 0 and not r.get("inherited"))
    keep = set(done[-1:]) | set(keep_extra)               # 最新と到達の出力だけ残す (準定常は系列から判定するので不要、ディスクのため。2026-10-10)
    for m in done:
        if m in keep: continue
        for f in list(run.glob(f"res_{m}.h5")) + list(run.glob(f"res_*_{m}.h5")) + list(run.glob(f"res*_{m}.xmf")): f.unlink()

def done_at(rows):
    """保存済みの系列の最初の終わりを ("REACH", step) か ("DIV", step) で返す (引き継いだ行の step は ≤ 0)。無ければ None。
    REACH は水準 (と E2) を --consec 個の連続した出力で満たした、その最後の出力 (§6.18)。"""
    run_ok = 0
    for n_ in range(len(rows)):
        r_ = rows[n_]
        if r_.get("nonfinite") or r_.get("nonfinite_all"): return ("DIV", r_["step"])
        dr_ = [drift(rows[: n_ + 1], k) for k in KEYS[:3]]
        e2_ = {k: 100 * (r_[k] / REF[k] - 1) for k in KEYS} if REF else None
        ok_ = r_.get("deficit") is not None and None not in dr_ and abs(r_["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr_) \
              and (e2_ is None or all(abs(v) <= 0.1 for v in e2_.values()))
        run_ok = run_ok + 1 if ok_ else 0
        if run_ok >= a.consec:
            return ("REACH", r_["step"])
    return None
def conclude(n):
    """終わりの出力 n で止める (何度呼んでも同じ結果): forge を止め、後の出力を消し、sha256 を記録して REACHED。"""
    ok = stop()
    if not ok: finish("EXEC_ERROR", note="forge を止められない", step=n); sys.exit(1)
    for m in steps():
        if m > n:
            for f in list(run.glob(f"res_{m}.h5")) + list(run.glob(f"res_*_{m}.h5")) + list(run.glob(f"res*_{m}.xmf")): f.unlink()
    if n <= 0:                                           # 分岐点 (引き継いだ行) で既に E を満たす: 追加の step は 0
        finish("REACHED", reach_step=0, note="分岐点で既に E を満たす (追加 0 step)"); sys.exit(0)
    if not (run / f"res_{n}.h5").exists(): finish("DATA_ERROR", note=f"到達の出力 res_{n}.h5 が無い", step=n); sys.exit(1)
    sha = hashlib.sha256((run / f"res_{n}.h5").read_bytes()).hexdigest()
    finish("REACHED", reach_step=n, reach_sha256=sha); cleanup(keep_extra=(n,)); sys.exit(0)
# 再開: 保存済みの系列で先に判定し直す (codex plan-3 M4)
_d = done_at(st["rows"])
if _d is not None:
    if _d[0] == "DIV": stop(); finish("DIVERGED", fail_step=_d[1]); sys.exit(0)
    conclude(_d[1])

T0 = time.time()
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
          _dn = done_at(st["rows"])                       # 連続の数え方は再開時と同じ関数で (§6.18 の --consec)
          if _dn is not None and _dn[0] == "REACH":
              save(); conclude(_dn[1])
          save(); cleanup()
      started = (run / "RUN_RC").exists() or time.time() - T0 > 600     # forge の起動前に見張りが始まっても待つ (起動の猶予 10 分)
      if started and not progressed and not alive and not [n for n in steps() if n > 0 and n not in {r["step"] for r in st["rows"]}]:
          rcf = run / "RUN_RC"; rc = rcf.read_text().strip() if rcf.exists() else "?"
          last = max((r["step"] for r in st["rows"]), default=0)
          nan = sorted(run.glob("res_nan_*.h5"))
          if rc == "0" and last >= a.budget: finish("CENSORED", last_step=last)
          elif nan: finish("DIVERGED", fail_step=int(nan[0].stem.split("_")[-1]), note="forge の detectNaN が止めた")
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
