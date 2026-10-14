"""§6.41・§6.42 の記録 (判定ではない): C・E・F・G9 の 12 本の check_quasisteady (run のディレクトリの既定の量) の VERDICT を集める。
有限期間 (2,000 step) の応答の試験なので、準定常は合格条件にしない。→ _band_ab/cold_pair/qs_record.json"""
import json, os, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
RUNS = ("run_0595_rq_c1", "run_0596_rq_c2", "run_0597_rq_c3", "run_0598_rq_e1", "run_0599_rq_e2", "run_0620_rq_e3",
        "run_0621_rq_f1", "run_0622_rq_f2", "run_0623_rq_f3", "run_0706_rp_g9a", "run_0707_rp_g9b", "run_0708_rp_g9c")
out = {}
for r in RUNS:
    if not (HERE / r / "RUN_RC").exists(): out[r] = "run が無い"; continue
    p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / r)], capture_output=True, text=True)
    lines = [l.strip() for l in (p.stdout + p.stderr).splitlines() if any(k in l for k in ("STEADY", "DRIFT", "OSCILL", "TRANSIENT", "OVERALL"))]
    out[r] = lines[:10]
    print(r, (lines[-1] if lines else "VERDICT の行が無い")[:120])
(HERE / "_band_ab/cold_pair/qs_record.json").write_text(json.dumps(out, indent=1, ensure_ascii=False)); print("→ _band_ab/cold_pair/qs_record.json")
