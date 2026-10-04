#!/usr/bin/env python3
"""S2 の一次記録を、レビュー可能な小さい形で集める (plan gradient-scalar-lsq-unification、codex result-2 M2)。

各 run について `s2_evidence/<case>/<run>/` に:
  residual_outer_end.csv.gz  残差履歴の outer_end 行だけ (判定ツールが読む行。元ファイルの sha256 を README に)
  solverConfig.yaml / bcondConfig.yaml / IC_FROM.txt / RUN_PROVENANCE.txt / CONVERGENCE_VERDICT.txt (あれば)
  passive_log.txt            forge_run.log の `[passive]` 行 (FCT の丸め前の収支値を含む、あれば)
  series.csv                 cond_series.csv / r1_series.csv (あれば)
さらに `s2_evidence/JUDGEMENTS.txt` に、元ファイルに対して判定ツールを実行したコマンドと全出力を書く。

  python3 s2_collect_evidence.py OUT_DIR   (リポジトリ直下で実行)
"""
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import sys

T = "solver_density_cuda/tools"
S = os.path.expanduser("~/sglsq/s2start")
RUNS = {
    "case/48.flat_plate_cooled_m4": ["run_0950_sglsq_s2_gg", "run_0951_sglsq_s2_lsq", "run_0952_sglsq_s2_gg_ext", "run_0953_sglsq_s2_lsq_ext",
                                     "run_0956_sglsq_s2_gg_chi0", "run_0957_sglsq_s2_base36d8", "run_0958_sglsq_s2_base36d8_ext"],
    "case/40.nozzle_design_tool": ["run_0950_sglsq_s2_gg", "run_0951_sglsq_s2_lsq", "run_0952_sglsq_s2_gg_ext", "run_0953_sglsq_s2_lsq_ext"],
    "case/39.periodic_hills": ["run_0950_sglsq_s2_gg", "run_0951_sglsq_s2_lsq"],
    "case/16.nozzle_wys": ["run_0950_sglsq_s2_0476_gg", "run_0951_sglsq_s2_0476_lsq", "run_0956_sglsq_s2_0476_gg_cont", "run_0957_sglsq_s2_0476_lsq_cont",
                           "run_0952_sglsq_s2_0482_gg", "run_0953_sglsq_s2_0482_lsq", "run_0954_sglsq_s2_fct_gg", "run_0955_sglsq_s2_fct_lsq",
                           "run_0958_sglsq_s2_fct_gg_nsub30", "run_0959_sglsq_s2_fct_lsq_nsub30", "run_0960_sglsq_s2_fct_gg_base36d8",
                           "run_0961_sglsq_s2_fct_gg_nokey"],
}
# 差し替え後の S2 判定 (plan §6 の差し替え規則): gg 双子 → lsq 双子
PAIRS = [("case/48.flat_plate_cooled_m4", "run_0952_sglsq_s2_gg_ext", "run_0953_sglsq_s2_lsq_ext"),
         ("case/40.nozzle_design_tool", "run_0952_sglsq_s2_gg_ext", "run_0953_sglsq_s2_lsq_ext"),
         ("case/39.periodic_hills", "run_0950_sglsq_s2_gg", "run_0951_sglsq_s2_lsq"),
         ("case/16.nozzle_wys", "run_0956_sglsq_s2_0476_gg_cont", "run_0957_sglsq_s2_0476_lsq_cont"),
         ("case/16.nozzle_wys", "run_0952_sglsq_s2_0482_gg", "run_0953_sglsq_s2_0482_lsq")]
FCT = ["run_0954_sglsq_s2_fct_gg", "run_0955_sglsq_s2_fct_lsq", "run_0958_sglsq_s2_fct_gg_nsub30",
       "run_0959_sglsq_s2_fct_lsq_nsub30", "run_0960_sglsq_s2_fct_gg_base36d8", "run_0961_sglsq_s2_fct_gg_nokey"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return f"$ {cmd}\n(rc={r.returncode})\n{r.stdout}{r.stderr}\n"


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    index = ["# S2 一次記録 (元ファイルは AWS `~/forge-pgrad-new/<case>/<run>/`)", "",
             "| run | residual_history.csv sha256 | 行数 (全体 / outer_end) |", "| --- | --- | --- |"]
    for case, runs in RUNS.items():
        for r in runs:
            src = os.path.join(case, r)
            dst = os.path.join(out, case.split("/")[1].split(".")[0], r)
            os.makedirs(dst, exist_ok=True)
            rh = os.path.join(src, "residual_history.csv")
            rows = list(csv.reader(open(rh)))
            head, body = rows[0], rows[1:]
            ph = head.index("phase") if "phase" in head else None
            keep = [x for x in body if ph is None or x[ph] == "outer_end"]
            with gzip.open(os.path.join(dst, "residual_outer_end.csv.gz"), "wt", newline="") as g:
                w = csv.writer(g); w.writerow(head); w.writerows(keep)
            index.append(f"| `{case}/{r}` | `{sha(rh)}` | {len(body)} / {len(keep)} |")
            for f in ("solverConfig.yaml", "bcondConfig.yaml", "IC_FROM.txt", "RUN_PROVENANCE.txt", "CONVERGENCE_VERDICT.txt"):
                if os.path.exists(os.path.join(src, f)):
                    shutil.copy2(os.path.join(src, f), dst)
            for f in ("cond_series.csv", "r1_series.csv"):
                if os.path.exists(os.path.join(src, f)):
                    shutil.copy2(os.path.join(src, f), os.path.join(dst, "series.csv"))
            log = os.path.join(src, "forge_run.log")
            pl = [l for l in open(log, errors="replace") if l.startswith("[passive")]
            if pl:
                open(os.path.join(dst, "passive_log.txt"), "w").writelines(pl)
    open(os.path.join(out, "INDEX.md"), "w").write("\n".join(index) + "\n")
    J = ["# 判定ツールを元ファイルに対して実行した記録 (plan §6 S2 の差し替え規則、codex result-2 M2)\n"]
    for case, gg, lsq in PAIRS:
        J.append(f"## {case}: gg {gg} → lsq {lsq}\n")
        J.append(sh(f"python3 {T}/check_convergence.py {case}/{gg} {case}/{lsq}"))
        J.append(sh(f"python3 {T}/check_floor_ratio.py --start {case}/{gg} {case}/{lsq}"))
    J.append("## case/16 系列の準定常 (plan §6 表の閾値)\n")
    for r in ("run_0952_sglsq_s2_0482_gg", "run_0953_sglsq_s2_0482_lsq"):
        J.append(sh(f"python3 {T}/check_quasisteady.py --series-csv case/16.nozzle_wys/{r}/cond_series.csv --series-cols onset_mm,dev_pct,g_exit,M_exit"))
    for r in ("run_0956_sglsq_s2_0476_gg_cont", "run_0957_sglsq_s2_0476_lsq_cont"):
        J.append(sh(f"python3 {T}/check_quasisteady.py --quantity machmax,pmax case/16.nozzle_wys/{r}"))
    for r in ("run_0950_sglsq_s2_gg", "run_0951_sglsq_s2_lsq"):
        J.append(sh(f"python3 {T}/check_quasisteady.py --series-csv case/39.periodic_hills/{r}/r1_series.csv --series-cols Cf_x05,Cf_x2,Cf_x6,xr_h,r_gradu,r_gradk,r_gradw,dF1_inf --drift 0.002 --osc 0.005 --tail 0.4"))
    J.append("## FCT smoke の保存収支 (checker 全出力 + 丸め前の remainder)\n")
    for r in FCT:
        J.append(sh(f"python3 {T}/check_passive_budget.py --mode fct case/16.nozzle_wys/{r}"))
        J.append(sh(f"grep -h 'fctCorr' case/16.nozzle_wys/{r}/forge_run.log | tail -1"))
    open(os.path.join(out, "JUDGEMENTS.txt"), "w").write("\n".join(J))
    print("ok", out)


if __name__ == "__main__":
    main()
