"""plan tooling-nozzle-core-grid §4.104 P-hs の見張りの判定 (出力ごとに cg_hs.sh が呼ぶ。最新の出力 s で止めるかを決める)。
  腕 b0 (H・B0): 停止規則に初めて届いた出力 s₁ を記録し (後から書き換えない)、s − 40k ≥ s₁ かつ [s − 40k, s] の全出力で停止規則を満たし、
                 同じ窓で 8 量の資格 (cg_hs_lib.qualify) を満たした最初の s で止める。
  腕 point (H・point): [s − 40k, s] で 8 量の資格を満たし、s で停止規則を満たした最初の s で止める。
  上限 (台本の nStepOuter) まで満たさなければ止めない (判定は cg_hs_judge.py)。
usage: python3 cg_hs_check.py <run> <b0|point> → 標準出力に STOP か CONTINUE、_band_ab/core_grid/hs_status_<run>.json を書く"""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cg_hs_lib as L  # noqa: E402
run, arm = sys.argv[1], sys.argv[2]
stp = L.OUTD / f"hs_status_{run}.json"
st = json.loads(stp.read_text()) if stp.exists() else {"run": run, "arm": arm, "first_reach": None, "checks": []}
rows, jd = L.stop_judgements(run)
ser = L.series(run)
dec = "CONTINUE"
if rows and ser:
    s = min(rows[-1]["step"], ser[-1]["step"])
    okmap = {j["step"]: bool(j.get("ok")) for j in jd}
    if st["first_reach"] is None:
        fr = [j["step"] for j in jd if j.get("ok")]
        if fr:
            st["first_reach"] = fr[0]
    rec = {"s": s, "stop_ok_at_s": okmap.get(s)}
    if s - L.CONF >= 0 and len([r for r in ser if s - L.CONF <= r["step"] <= s]) == 17:
        held = all(okmap.get(k) for k in range(s - L.CONF, s + 1, 2500))
        rec["stop_held_40k"] = held
        cond_pre = (arm == "point" and okmap.get(s)) or (arm == "b0" and st["first_reach"] is not None and s - L.CONF >= st["first_reach"] and held)
        if cond_pre:
            q = L.qualify(ser, s)
            rec["qual_ok"] = q["ok"]; rec["qual_fail"] = [k for k, v in q["q"].items() if not v["ok"]]
            if q["ok"]:
                dec = "STOP"
    rec["decision"] = dec
    st["checks"].append(rec)
    if dec == "STOP":
        st["stop_step"] = s
stp.write_text(json.dumps(st, ensure_ascii=False, indent=1))
print(dec)
