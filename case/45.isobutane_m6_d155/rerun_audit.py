"""plan tooling-rerun-conditions §6 の検証 run を監査できる形にする (codex result 2026-10-06 Major 4)。
usage: python3 rerun_audit.py RUN [RUN ...]   (case/45 の run ディレクトリがある場所で)
各 RUN に QS_VERDICT.txt (§6 の判定規約で check_quasisteady をかけたコマンドと出力) を書き、末尾 5 枚の平均・幅・直前窓差と、
最終場があれば逆流節点数 (Ux < 0) を標準出力に出す。"""
import csv, subprocess, sys, datetime
from pathlib import Path
import numpy as np
TOOLS = Path(__file__).resolve().parents[2] / "solver_density_cuda" / "tools"
THR = {"delta_E": ("5e-5", "1e-4"), "exitM_A": ("3e-6", "3e-6"), "mdot": ("1e-5", "2e-5")}
for name in sys.argv[1:]:
    run = Path(name); qs = run / "quantities_series.csv"
    if not qs.exists():
        print(f"{name}: quantities_series.csv なし"); continue
    rows = list(csv.DictReader(open(qs)))
    lines = [f"# {datetime.datetime.now().isoformat(timespec='seconds')} plan tooling-rerun-conditions §6 の判定規約 (--tail 0.4 --min-snaps 10、量ごとの閾値)"]
    summ = []
    for c, (dr, osc) in THR.items():
        if c not in rows[0]:
            continue
        cmd = [sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(qs), "--series-cols", c,
               "--tail", "0.4", "--min-snaps", "10", "--drift", dr, "--osc", osc]
        out = subprocess.run(cmd, capture_output=True, text=True).stdout
        lines.append("$ " + " ".join(cmd[1:]))
        lines += [l for l in out.splitlines() if "csv" in l or c in l]
        v = np.array([float(r[c]) for r in rows]); f5, p5 = v[-5:], v[-10:-5]
        vd = next((w for w in ("STEADY", "DRIFTING", "OSCILLATING", "TRANSIENT-UNSETTLED", "NONFINITE") if any(w in l and c in l for l in out.splitlines())), "?")
        summ.append(f"{c} {f5.mean():.7g} (幅 {np.ptp(f5):.3g}, 直前窓差 {f5.mean() - p5.mean() if len(p5) == 5 else float('nan'):.3g}) {vd}")
    (run / "QS_VERDICT.txt").write_text("\n".join(lines) + "\n")
    rev = ""
    res = sorted(run.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    if res and res[-1].stem != "res_0":
        import h5py
        with h5py.File(res[-1]) as f:
            rev = f"; {res[-1].name}: Ux<0 {int((f['VALUE/Ux'][:] < 0).sum())} 節点"
    print(f"{name}: " + " | ".join(summ) + rev)
