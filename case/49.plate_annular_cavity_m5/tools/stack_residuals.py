#!/usr/bin/env python3
"""段階起動の残差履歴を 1 本に繋いで、起動全体の収束 VERDICT を出す。

段階起動 (S0 slip -> S1 laminar -> S2 等温 -> S3/S4 SST -> S5 2次ランプ -> S6 本段) では、
各段が前段の収束場から始まるので**本段だけを `check_convergence.py` にかけると残差の低下桁数が
小さく出て必ず "STALLED (plateau)" になる** (絶対値は既に床まで落ちている)。起動全体で
何桁落ちたかを見るため、`residual_history_<段>.csv` と本段の `residual_history.csv` を
step をずらして連結し、それに対して正本ツール `check_convergence.py` を回す。

usage:
  python3 tools/stack_residuals.py <run_dir> [...] [--keep] [--plot]

出力:
  <run_dir>/residual_history_all.csv    全段連結 (**起動診断用**)
  <run_dir>/residual_history_turb.csv   乱流を解いた段以降
  <run_dir>/residual_history_same.csv   **同じ方程式・BC・空間離散化の区間だけ** (判定用)

**判定に使うのは `_same` (既定 S5 ランプ以降 = 2 次・SST・等温で固定) である** (2026-09-19
codex Major)。全段連結は slip / 層流 / 断熱 / 1 次 の大きな過渡を含み、`check_convergence.py` が
系列全体の最大値を低下桁数の基準に取るため、**別の境界条件の過渡を本段の基準にしてしまう**。
CFL だけが違う区間の連結は可、残差の定義と正規化が同じことを前提にする。
収束済み場からの継続は `check_convergence.py --from-floor <参照 run>` を使う。
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CHECK = ROOT / "solver_density_cuda" / "tools" / "check_convergence.py"
ORDER = ["S0_slip", "S1_lam", "S2_iso", "S3_sst_soft", "S4_sst_mid"]


def stage_files(rd):
    """起動順に (ラベル, csv) を返す。S5 ランプは番号順、最後に本段。"""
    out = []
    for tag in ORDER:
        p = rd / ("residual_history_%s.csv" % tag)
        if p.exists():
            out.append((tag, p))
    ramps = sorted(rd.glob("residual_history_S5_ramp*.csv"),
                   key=lambda p: int(p.stem.split("ramp")[1].split("_")[0]))
    out += [(p.stem.replace("residual_history_", ""), p) for p in ramps]
    main = rd / "residual_history.csv"
    if main.exists():
        out.append(("S6_main", main))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--keep", action="store_true", help="連結 csv を消さない (既定でも残す)")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--same-from", default="S5",
                    help="判定区間の開始段プレフィックス (既定 S5 = 2 次ランプ以降)")
    a = ap.parse_args()
    rc = 0
    for run in a.runs:
        rd = Path(run)
        files = stage_files(rd)
        print("\n================ %s ================" % run)
        if not files:
            print("  残差履歴が無い"); rc = 1; continue
        # 段ごとに列が違う (S0-S2 は乱流無しなので rms_roK/rms_roOmega が無い) ので
        # **列名でそろえる**。欠けている列は空欄にする (check_convergence 側が無視する)。
        staged, hdr, off, bounds = [], [], 0, []
        for tag, p in files:
            with open(p) as f:
                r = list(csv.reader(f))
            if len(r) < 2:
                continue
            h = [c.strip() for c in r[0]]
            for c in h:
                if c not in hdr:
                    hdr.append(c)
            body = [x for x in r[1:] if x and x[0].strip().lstrip("-").isdigit()]
            if not body:
                continue
            n0, n1 = int(body[0][0]), int(body[-1][0])
            staged.append((tag, h, body, off, n0, n1))
            bounds.append((tag, off + n0, off + n1, len(body)))
            off += n1 + 1
        def write(path, cols, subset):
            """cols だけを列に持つ連結 csv を書く (subset の段だけ使う)。"""
            rows = []
            for tag, h, body, off0, n0, n1 in subset:
                idx = {c: i for i, c in enumerate(h)}
                for x in body:
                    row = [x[idx[c]] for c in cols]
                    row[0] = str(int(x[0]) + off0)
                    rows.append(row)
            with open(path, "w", newline="") as f:
                w = csv.writer(f); w.writerow(cols); w.writerows(rows)
            return len(rows)

        for tag, s0, s1, n in bounds:
            print("  %-18s step %7d .. %7d  (%d 行)" % (tag, s0, s1, n))

        # **全段に共通する列だけ**を連結する (S0-S2 は乱流方程式を解いていないので
        # rms_roK/rms_roOmega の行が存在しない。空欄で埋めると check_convergence が落ちる)。
        common = [c for c in hdr if all(c in h for _, h, _, _, _, _ in staged)]
        write(rd / "residual_history_all.csv", common, staged)
        turb = [st for st in staged if "rms_roK" in st[1]]
        outs = [("起動全体 (共通列) — **起動診断用。合否判定に使わない**",
                 rd / "residual_history_all.csv")]
        if turb and len(turb) < len(staged):
            tcols = [c for c in hdr if all(c in h for _, h, _, _, _, _ in turb)]
            write(rd / "residual_history_turb.csv", tcols, turb)
            outs.append(("乱流を解いた段以降 (%s 〜) — 参考" % turb[0][0],
                         rd / "residual_history_turb.csv"))
        # **判定区間**: 数値設定が最終形で固定されている段だけ (既定 S5 ランプ以降)
        same = [st for st in staged if st[0].startswith(a.same_from) or st[0] == "S6_main"]
        if same:
            scols = [c for c in hdr if all(c in h for _, h, _, _, _, _ in same)]
            write(rd / "residual_history_same.csv", scols, same)
            outs.append(("**判定区間** (%s 〜 = 同じ方程式・BC・空間離散化)" % same[0][0],
                         rd / "residual_history_same.csv"))
        for label, path in outs:
            print("  --- %s ---" % label)
            cmd = [sys.executable, str(CHECK), str(path)]
            if a.plot:
                cmd.append("--plot")
            r = subprocess.run(cmd, capture_output=True, text=True)
            print(r.stdout.strip() or r.stderr.strip()[-800:])
            rc = max(rc, r.returncode)
    return rc


if __name__ == "__main__":
    sys.exit(main())
