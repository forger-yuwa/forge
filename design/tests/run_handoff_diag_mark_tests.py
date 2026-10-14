"""受け渡しの関門が診断専用の run (印のある run) を拒むことの試験 (plan tooling-nozzle-isothermal-wall-chain §6 V-h370-diag、2026-10-14)。

印 (`deltastar_loop.DIAG_ONLY_MARKS`) のある run に `handoff_gate` を呼ぶと、上書きの有無によらず HandoffBlocked で止まり、
run には何も書かない。印の無い run では今までどおり先へ進む (ここでは印の検査を通り越すことだけを見る)。
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

fails = 0


def check(name, ok):
    global fails
    print(("ok   " if ok else "FAIL ") + name)
    fails += 0 if ok else 1


for mark in DL.DIAG_ONLY_MARKS:
    for kw in ({}, {"allow_not_reached": "試験", "allow_unconverged": "試験", "allow_unsettled": "試験"}):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / mark).write_text("{}")
            before = sorted(p.name for p in d.iterdir())
            try:
                DL.handoff_gate(d, **kw)
                blocked = False
            except DL.HandoffBlocked as e:
                blocked = "診断専用" in str(e)
            after = sorted(p.name for p in d.iterdir())
            check(f"印 {mark} の run は止まる (上書き {'あり' if kw else 'なし'})", blocked)
            check(f"印 {mark} の run には何も書かない (上書き {'あり' if kw else 'なし'})", before == after)

with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    try:
        r = DL.handoff_gate(d, write=False, extract=False)
        passed_mark = True
        check("印の無い run は印の検査を通り越す (記録の無い run として止まる)", r.get("proceed") is False)
    except DL.HandoffBlocked as e:
        check("印の無い run は印の検査を通り越す", "診断専用" not in str(e))

print("ALL PASS" if fails == 0 else f"FAIL {fails}")
sys.exit(1 if fails else 0)
