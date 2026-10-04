#!/usr/bin/env python3
"""A/B (縦すきま上流端: 壁 vs 上流横すきま) を事前登録の基準で判定する。

閾値は**ここに書かない** — `acceptance.json` の `T4-3D-AB` から読む
(結果を見てから閾値を動かさないため)。両腕の `check_quasisteady` の VERDICT が
`STEADY` であることを先に確かめ、揃っていなければ判定を拒否する。

    python3 tools/ab_compare.py --a run_0049_ab3d_wall --b run_0050_ab3d_open
"""
import argparse, json, re, subprocess, sys
from pathlib import Path
import numpy as np

CASE = Path(__file__).resolve().parents[1]
ROOT = CASE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
KEYS = ["q_tc92", "q_tc91", "q_tc90", "p_fwd_max"]


def series(run):
    out = CASE / run / "_ab_series.csv"
    subprocess.run([sys.executable, str(CASE / "tools" / "gap3d_eval.py"),
                    "--run", run, "--csv", str(out)], check=True,
                   stdout=subprocess.DEVNULL)
    raw = np.loadtxt(out, delimiter=",", skiprows=1)
    cols = out.read_text().splitlines()[0].split(",")
    if raw.ndim == 1:
        raw = raw[None, :]
    return {c: raw[:, i] for i, c in enumerate(cols)}, out


def quasisteady(csv):
    r = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"),
                        "--series-csv", str(csv), "--series-cols", ",".join(KEYS)],
                       capture_output=True, text=True)
    return r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--force", action="store_true",
                    help="定常でなくても数値を出す (判定は出さない)")
    a = ap.parse_args()

    acc = json.loads((CASE / "acceptance.json").read_text(encoding="utf-8"))
    gate = {g["id"]: g for g in acc["gates"]}["T4-3D-AB"]
    thr = None
    for c in gate["criteria"]:
        m = re.search(r"(\d+)\s*%", c["criterion"])
        if m:
            thr = float(m.group(1)) / 100.0
            break
    if thr is None:
        sys.exit("acceptance.json から閾値を読めなかった")
    print(f"事前登録 {gate['id']} ({gate['status']}): 閾値 |A-B|/max < {thr*100:.0f} %")

    res, ok = {}, True
    for tag, run in (("A", a.a), ("B", a.b)):
        s, csv = series(run)
        res[tag] = s
        qs = quasisteady(csv)
        v = [l for l in qs.splitlines() if "VERDICT" in l or "STEADY" in l or "DRIFT" in l]
        print(f"\n--- {tag} = {run}  ({len(s['step'])} スナップショット)")
        for l in v[:8]:
            print("   ", l.strip())
        if "STEADY" not in qs or "DRIFTING" in qs or "TRANSIENT" in qs:
            ok = False

    print("\n量          腕 A            腕 B          |A-B|/max")
    verdict = True
    for k in KEYS:
        va, vb = float(res["A"][k][-1]), float(res["B"][k][-1])
        rel = abs(va - vb) / max(abs(va), abs(vb), 1e-30)
        flag = "OK" if rel < thr else "**超過**"
        if rel >= thr:
            verdict = False
        print(f"{k:11s} {va:13.4g} {vb:13.4g}   {rel*100:8.2f} %  {flag}")

    if not ok and not a.force:
        print("\nVERDICT: 判定しない — どちらかの腕が STEADY でない (--force で数値だけ出せる)")
        sys.exit(2)
    print(f"\nVERDICT: {'簡略 2 は使える (腕 A を生産構成にする)' if verdict else '簡略 2 は使えない (腕 B を生産構成にする)'}")


if __name__ == "__main__":
    main()
