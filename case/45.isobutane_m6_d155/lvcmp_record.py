#!/usr/bin/env python3
"""見張り (m9_watch.py) が途中の出力を消す前に、各出力の sha256・全保存量の有限性・ρ・P・T の最小値を記録する
(plan time_integration-line-viscous-jacobian-dt §6.3、codex plan-2 の m5)。見張りは速度 plan のものなので変えず、横で並べて動かす。
usage: python3 lvcmp_record.py <run>   → <run>/lvcmp_outputs.json ({"rows": {step: {...}}}、原子的に保存)
- res_<N>.h5 (N > 0) を、最終更新から 20 s 経ったものから 1 つずつ読む (見張りは 60 s 経ってから読むので、それより先に読める)。
- forge が終わり、記録していない出力が残っていなければ終わる。読めない出力は 3 回まで再試行し、だめなら "error" を記録する。
"""
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

run = Path(sys.argv[1]).resolve()
SF = run / "lvcmp_outputs.json"
st = json.loads(SF.read_text()) if SF.exists() else {"run": run.name, "rows": {}, "tries": {}}


def save():
    tmp = SF.with_suffix(".tmp"); tmp.write_text(json.dumps(st, indent=1, ensure_ascii=False)); os.replace(tmp, SF)


def alive():
    for p in subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split():
        try:
            if os.path.realpath(f"/proc/{p}/cwd") == str(run):
                return True
        except OSError:
            pass
    return False


def outputs():
    return sorted(int(m.group(1)) for f in run.glob("res_*.h5") if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)) and int(m.group(1)) > 0)


t0 = time.time()
while True:
    a = alive()
    pending = [n for n in outputs() if str(n) not in st["rows"]]
    did = False
    for n in pending:
        p = run / f"res_{n}.h5"
        try:
            if a and time.time() - p.stat().st_mtime < 20:
                continue
            sha = hashlib.sha256(p.read_bytes()).hexdigest()
            nonfin = 0
            mins = {}
            with h5py.File(p, "r") as h:
                for k, d in h["VALUE"].items():
                    x = np.asarray(d[:])
                    if x.dtype.kind == "f":
                        nonfin += int(np.count_nonzero(~np.isfinite(x)))
                for k in ("ro", "P", "T"):
                    x = np.asarray(h["VALUE"][k][:], dtype=np.float64)
                    mins[k] = float(np.nanmin(x)) if x.size else None
            st["rows"][str(n)] = {"sha256": sha, "nonfinite_all": nonfin, "min": mins,
                                  "positive": all(v is not None and v > 0 for v in mins.values())}
            did = True
        except FileNotFoundError:
            st["rows"].setdefault(str(n), {"error": "記録の前に消えた"})
        except Exception as e:
            st["tries"][str(n)] = st["tries"].get(str(n), 0) + 1
            if st["tries"][str(n)] >= 3:
                st["rows"][str(n)] = {"error": repr(e)}
        save()
    if not a and time.time() - t0 > 600 and not [n for n in outputs() if str(n) not in st["rows"]]:
        st["finished"] = True; save(); break
    if not did:
        time.sleep(10)
