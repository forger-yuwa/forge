"""plan tooling-nozzle-throat-monotone-r2 §6 N・K: NS run の量別の準定常判定 (check_quasisteady --series-csv、末尾 5 枚) を
nozzle_report --verdicts 用の JSON にする。比較のため、旧壁の run_0117 / run_0118 にも同じコマンドをかける (run_0117 の既存 verdicts と同じ窓)。
  dry (quantities_series.csv、5000 step ごと 12 枚): exitM_A → exit_core_M、overshoot01 → overshoot_eta0.1、wave01 → wave_eta0.1、dE_over_dC
  cond (cond_series.csv、1000 step ごと 18 枚): onset_x_axis → cond_onset_x_axis、S_max → cond_S_max、exit_g_core → cond_exit_g_core、exit_core_M
判定は check_quasisteady の VERDICT をそのまま使う (末尾 5 枚 = --tail 4.5/n)。値は末尾 5 枚の平均と幅を note に書く。
usage: python3 throat_mono_ns_verdicts.py RUN {dry|cond} [--note KEY=TEXT ...]  → _band_ab/verdicts_<RUN 名の run_NNNN>.json と標準出力
"""
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

C = Path(__file__).resolve().parent
TOOLS = C.parents[1] / "solver_density_cuda/tools"
MAP = {"dry": ("quantities_series.csv", {"exitM_A": "exit_core_M", "overshoot01": "overshoot_eta0.1", "wave01": "wave_eta0.1",
                                         "dE_over_dC": "dE_over_dC"}),
       "cond": ("cond_series.csv", {"onset_x_axis": "cond_onset_x_axis", "S_max": "cond_S_max", "exit_g_core": "cond_exit_g_core",
                                    "exit_core_M": "exit_core_M"})}


def main():
    run, mode = C / sys.argv[1], sys.argv[2]
    notes = dict(a.split("=", 1) for a in sys.argv[3:] if "=" in a)
    fname, cols = MAP[mode]
    f = run / fname
    rows = list(csv.DictReader(open(f)))
    n = len(rows)
    if n < 6:
        raise SystemExit(f"{f}: 行が {n} しかない (末尾 5 枚の判定には 6 枚以上)")
    tail = 4.5 / n
    steps = [int(float(r["step"])) for r in rows][-5:]
    out, table = {}, []
    for col, key in cols.items():
        q = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(f),
                            "--series-cols", col, "--tail", f"{tail:.6f}"], capture_output=True, text=True)
        txt = q.stdout + q.stderr
        m = re.findall(r"\b(STEADY|OSCILLATING|TRANSIENT-UNSETTLED|DRIFTING|NONFINITE)\b", txt.split(col, 1)[-1] if col in txt else txt)
        if not m:
            raise SystemExit(f"{col}: check_quasisteady の VERDICT を読めない:\n{txt[-800:]}")
        v = [float(r[col]) for r in rows][-5:]
        mean, rng = sum(v) / 5, max(v) - min(v)
        note = f"末尾 5 枚の平均 {mean:.7g}・幅 {rng:.3g}"
        if key in notes:
            note += "。" + notes[key]
        out[key] = {"verdict": m[0], "window": f"{run.name} res_{steps[0]}–{steps[-1]} ({fname})", "note": note}
        table.append((key, m[0], mean, rng))
    tag = re.match(r"(run_\d{4})", run.name).group(1)
    (C / "_band_ab").mkdir(exist_ok=True)
    dst = C / "_band_ab" / f"verdicts_{tag}_monoeval.json"
    dst.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for key, vd, mean, rng in table:
        print(f"{run.name} {key:20s} {vd:20s} 平均 {mean:.7g}  幅 {rng:.3g}")
    print("->", dst)


if __name__ == "__main__":
    main()
