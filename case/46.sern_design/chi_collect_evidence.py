#!/usr/bin/env python3
"""chi 既定化 (plan convection-slau-wall-normal-chi-default) の B0/B1 一次記録をレビュー可能な形で集める (codex result 2026-09-27 M1)。

各 run について `<OUT>/<case>/<run 名から run_ を除く>/` に: 設定 (solverConfig/bcondConfig/IC_FROM/RUN_PROVENANCE/
CONVERGENCE_VERDICT)、起動記録 `forge_launches.jsonl`、残差の outer_end 行 (`residual_outer_end.csv.gz`)、派生量 CSV
(ct.csv / chidef_*.csv / point_probe_*.out)。`<OUT>/INDEX.md` に元ファイルと入力 h5 の sha256、`<OUT>/JUDGEMENTS.txt` に
元ファイルへ判定ツール (check_convergence・check_floor_ratio --factor 2) を当てた全出力を書く。
リポジトリ直下で実行: python3 case/46.sern_design/chi_collect_evidence.py OUT_DIR
"""
import csv
import glob
import gzip
import hashlib
import os
import shutil
import subprocess
import sys

T = "solver_density_cuda/tools"
PAIRS = [  # (case, 明示 0, 省略)
    ("case/46.sern_design", "run_0965_chidef_2d_m4_off_flag0", "run_0966_chidef_2d_m4_off_omit"),
    ("case/46.sern_design", "run_0967_chidef_2d_m10_on_flag0", "run_0968_chidef_2d_m10_on_omit"),
    ("case/16.nozzle_wys", "run_0965_chidef_sst_flag0", "run_0966_chidef_sst_omit"),
    ("case/39.periodic_hills", "run_0952_chidef_b1p_flag0", "run_0953_chidef_b1p_omit"),
    ("case/40.nozzle_design_tool", "run_0954_chidef_b1a_flag0", "run_0955_chidef_b1a_omit"),
    ("case/40.nozzle_design_tool", "run_0956_chidef_b1a_isoT_flag0", "run_0957_chidef_b1a_isoT_omit"),
    ("case/40.nozzle_design_tool", "run_0958_chidef_b1a_isoT_flag0_ext", "run_0959_chidef_b1a_isoT_omit_ext"),
    ("case/16.nozzle_wys", "run_0963_chidef_b1c_flag0", "run_0964_chidef_b1c_omit"),
]
SINGLES = [("case/46.sern_design", "run_0961_chidef_junction_omit"), ("case/46.sern_design", "run_0960_chidef_junction_warm"),
           ("case/46.sern_design", "run_0963_chidef_2d_m4_off_base"), ("case/46.sern_design", "run_0964_chidef_2d_m10_on_base"),
           ("case/16.nozzle_wys", "run_0962_chidef_sst_base")]
KEEP = ("solverConfig.yaml", "bcondConfig.yaml", "IC_FROM.txt", "RUN_PROVENANCE.txt", "CONVERGENCE_VERDICT.txt", "forge_launches.jsonl",
        "prepare_info.json", "probe.yaml")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return f"$ {cmd}\n(rc={r.returncode})\n{r.stdout}{r.stderr}\n"


def collect(out, case, run, index):
    src = os.path.join(case, run)
    dst = os.path.join(out, case.split("/")[1].split(".")[0], run[4:])
    os.makedirs(dst, exist_ok=True)
    for f in KEEP:
        if os.path.exists(os.path.join(src, f)):
            shutil.copy2(os.path.join(src, f), dst)
    for pat in ("ct.csv", "chidef_*.csv", "point_probe_*.out", "force_history.csv", "cond_series.csv"):
        for f in glob.glob(os.path.join(src, pat)):
            shutil.copy2(f, dst)
    rh = os.path.join(src, "residual_history.csv")
    if os.path.exists(rh):
        rows = list(csv.reader(open(rh)))
        ph = rows[0].index("phase") if "phase" in rows[0] else None
        with gzip.open(os.path.join(dst, "residual_outer_end.csv.gz"), "wt", newline="") as g:
            w = csv.writer(g); w.writerow(rows[0]); w.writerows(x for x in rows[1:] if ph is None or x[ph] == "outer_end")
        index.append(f"| `{case}/{run}` | residual `{sha(rh)[:16]}` |")
    for h5 in sorted(glob.glob(os.path.join(src, "*.h5"))):
        if not os.path.basename(h5).startswith("res_"):
            index.append(f"| `{case}/{run}` | 入力 `{os.path.basename(h5)}` `{sha(h5)[:16]}` |")


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    index = ["# chi 既定化 B0/B1 一次記録 (元は AWS `~/forge-pgrad-new/<case>/<run>/`)", "", "| run | sha256 (先頭 16) |", "| --- | --- |"]
    J = ["# 判定ツールを元ファイルに当てた全出力 (plan convection-slau-wall-normal-chi-default §6、codex result 2026-09-27 M1)\n"]
    for case, a, b in PAIRS:
        collect(out, case, a, index); collect(out, case, b, index)
        J.append(f"## {case}: 明示 0 {a} / 省略 {b}\n")
        J.append(sh(f"python3 {T}/check_convergence.py {case}/{a} {case}/{b}"))
        J.append(sh(f"python3 {T}/check_floor_ratio.py --start {case}/{a} --factor 2 {case}/{b}"))
    for case, r in SINGLES:
        collect(out, case, r, index)
        J.append(sh(f"python3 {T}/check_convergence.py {case}/{r}"))
    for b in ("forge_2fa3826c", "forge_22976398"):
        p = os.path.expanduser(f"~/sglsq/{b}")
        if os.path.exists(p):
            index.append(f"| バイナリ `{b}` | `{sha(p)[:16]}` |")
    open(os.path.join(out, "INDEX.md"), "w").write("\n".join(index) + "\n")
    open(os.path.join(out, "JUDGEMENTS.txt"), "w").write("\n".join(J))
    print("ok", out)


if __name__ == "__main__":
    main()
