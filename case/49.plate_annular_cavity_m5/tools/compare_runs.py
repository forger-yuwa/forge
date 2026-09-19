#!/usr/bin/env python3
"""複数 run (格子・条件違い) の結論量を並べて比較する。格子収束ゲート (plan §6.4) 用。

usage:
  python3 tools/compare_runs.py run_0002_stageA_qwall run_0003_stageB_cpg run_0004_stageC_cpg [--plot]

各 run の `cavity_eval.json` (無ければその場で評価) から:
  キャビティ 3 壁の総入熱 Q / 面平均 q'' / h_ref / 開口ガス温度 / 深さ中央・底の温度 / 侵入深さ
を取り、**最も細かい run を基準にした相対差**を出す。節点数は mesh h5 から読む。
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent

KEYS = [("Q_cav", "キャビティ3壁 Q [W]"), ("qpp_cav", "面平均 q'' [W/m2]"),
        ("h_ref", "h_ref [W/m2K]"), ("T_mouth", "開口ガス T [K]"),
        ("T_mid", "中央 T [K]"), ("T_floor", "底 T [K]"), ("zpen25", "侵入(25K) [mm]"),
        ("mdot_in", "開口流入 [kg/s]")]


def nodes_of(run):
    m = Path(run) / "mesh.h5"
    if not m.exists():
        return float("nan")
    with h5py.File(m, "r") as f:
        return int(f["VALUE/ro"].shape[0])


def load(run, force=False):
    j = Path(run) / "cavity_eval.json"
    if force or not j.exists():
        subprocess.run([sys.executable, str(HERE / "cavity_eval.py"), str(run)],
                       check=True, capture_output=True, text=True)
    d = json.loads(j.read_text())
    f, w = d["field"], d["wall"]
    cav = [g for g in ("cav_outer", "cyl_side", "cav_floor") if g in w]
    Q = sum(w[g]["Q_W"] for g in cav)
    A = sum(w[g]["area_m2"] for g in cav)
    hr = sum(w[g].get("h_ref", float("nan")) * w[g]["area_m2"] for g in cav) / max(A, 1e-30)
    return {"Q_cav": Q, "qpp_cav": Q / max(A, 1e-30), "h_ref": hr,
            "T_mouth": 500.0 + f["dT_mouth"], "T_mid": 500.0 + f["dT_mid"],
            "T_floor": 500.0 + f["dT_floor"], "zpen25": f["zpen_25"] * 1e3,
            "mdot_in": f["mdot_in"], "nodes": nodes_of(run)}


def richardson(rows):
    r"""ASME (Celik 2008) の格子収束手順。rows は節点数の**昇順**で渡す。

    代表寸法は $h \propto N^{-1/3}$。細分比が一定でない場合も扱える一般形を使う
    (ここは scale 0.7/1.0/1.4 なので r ≈ 1.43 / 1.40 でわずかに違う)。
      p を  p = |ln|e32/e21| + q(p)| / ln r21,   q(p) = ln((r21^p - s)/(r32^p - s))
    で不動点反復し、外挿値 f_ext = (r21^p f1 - f2)/(r21^p - 1)、
    GCI_fine = 1.25 |(f1-f2)/f1| / (r21^p - 1) を出す。
    **p が 0 以下や極端に大きい (>4) 場合は単調収束していない**ので値を信じない。
    """
    f_ = [d for _, d in rows]
    h = [d["nodes"] ** (-1.0 / 3.0) for d in f_]
    # h1 = 最細
    i1, i2, i3 = len(rows) - 1, len(rows) - 2, len(rows) - 3
    r21, r32 = h[i2] / h[i1], h[i3] / h[i2]
    print("\n--- Richardson 外挿 (ASME/Celik) ---")
    print("  最細 %s (%d 節点) / 中 %s / 粗 %s,  r21 = %.4f, r32 = %.4f"
          % (Path(rows[i1][0]).name, f_[i1]["nodes"], Path(rows[i2][0]).name,
             Path(rows[i3][0]).name, r21, r32))
    print("  %-16s %12s %12s %12s %8s %14s %10s %9s"
          % ("量", "粗", "中", "最細", "p", "外挿値", "GCI[%]", "判定"))
    for k, lab in KEYS:
        f1, f2, f3 = f_[i1][k], f_[i2][k], f_[i3][k]
        e21, e32 = f2 - f1, f3 - f2
        if abs(e21) < 1e-30 or abs(e32) < 1e-30:
            print("  %-16s %12.5g %12.5g %12.5g  (差がゼロ)" % (lab, f3, f2, f1)); continue
        sgn = 1.0 if (e32 / e21) > 0 else -1.0
        pp = 2.0
        for _ in range(200):
            q = np.log((r21 ** pp - sgn) / (r32 ** pp - sgn))
            pn = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
            if abs(pn - pp) < 1e-10:
                pp = pn; break
            pp = 0.5 * pp + 0.5 * pn
        fext = (r21 ** pp * f1 - f2) / (r21 ** pp - 1.0)
        gci = 100.0 * 1.25 * abs((f1 - f2) / f1) / (r21 ** pp - 1.0)
        note = "OK" if (0.5 <= pp <= 4.0 and sgn > 0) else ("振動" if sgn < 0 else "次数外")
        print("  %-16s %12.5g %12.5g %12.5g %8.2f %14.5g %10.2f %9s"
              % (lab, f3, f2, f1, pp, fext, gci, note))
    print("  注: 単調 (sgn>0) かつ 0.5<=p<=4 でないものは外挿値を信じない。GCI は最細格子の数値不確かさ。")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--force", action="store_true", help="cavity_eval を回し直す")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--richardson", action="store_true",
                    help="3 点以上で観測次数 p・外挿値・GCI を出す (ASME/Celik 手順)")
    a = ap.parse_args()
    rows = [(r, load(CASE / r if not Path(r).is_absolute() else r, a.force)) for r in a.runs]
    rows.sort(key=lambda t: t[1]["nodes"])
    fin = rows[-1][1]

    print("%-26s %9s " % ("run", "nodes") + " ".join("%14s" % lab for _, lab in KEYS))
    for r, d in rows:
        print("%-26s %9d " % (Path(r).name, d["nodes"])
              + " ".join("%14.5g" % d[k] for k, _ in KEYS))
    print("\n--- 最細 run (%s, %d 節点) からの相対差 [%%] ---" % (Path(rows[-1][0]).name, fin["nodes"]))
    print("%-26s %9s " % ("run", "nodes") + " ".join("%14s" % lab for _, lab in KEYS))
    for r, d in rows[:-1]:
        print("%-26s %9d " % (Path(r).name, d["nodes"])
              + " ".join("%14.2f" % (100.0 * (d[k] / fin[k] - 1.0)) for k, _ in KEYS))
    if len(rows) >= 3:
        pen = [abs(rows[-2][1][k] / fin[k] - 1.0) * 100 for k, _ in KEYS]
        print("\n最細 2 点の差の最大 = %.2f %%  (plan §6.4 の格子感度ゲートは 5 %%)" % max(pen))

    if a.richardson and len(rows) >= 3:
        richardson(rows)

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import font_manager
        for fp in Path.home().joinpath(".fonts").glob("NotoSansCJKjp-Regular.otf"):
            font_manager.fontManager.addfont(str(fp))
            matplotlib.rcParams["font.family"] = font_manager.FontProperties(fname=str(fp)).get_name()
        h = np.array([d["nodes"] ** (-1 / 3) for _, d in rows])
        fig, ax = plt.subplots(2, 4, figsize=(15.5, 7))
        for i, (k, lab) in enumerate(KEYS):
            b = ax[i // 4][i % 4]
            b.plot(h, [d[k] for _, d in rows], "o-", lw=2, color="#1d4ed8")
            b.set_xlabel("$N^{-1/3}$ (格子の細かさ →)")
            b.set_title(lab, fontsize=11)
            b.grid(alpha=.3)
            b.set_xlim(0, max(h) * 1.1)
        fig.suptitle("case/49 格子収束 (" + " / ".join("%s %dk" % (Path(r).name.split("_")[2], d["nodes"] // 1000)
                                                      for r, d in rows) + ")", fontsize=13)
        fig.tight_layout()
        out = CASE / "grid_convergence.png"
        fig.savefig(out, dpi=130, bbox_inches="tight")
        print("wrote", out)


if __name__ == "__main__":
    main()
