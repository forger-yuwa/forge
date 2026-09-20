#!/usr/bin/env python3
r"""NASA CR-168015 (Hylton et al. 1983) から**翼型座標と試験条件**を抽出する。

一次資料: `papers/cht/NASA-CR-168015_Hylton_1983.pdf` (NTRS 19830020105, **git 追跡外**)。
スキャン PDF の OCR なので、**数値は必ず検証してから書き出す**:

1. 表 II (Mark II) / III (C3X) は **cm と inch の両方**を載せているので、
   `cm == 2.54 * in` (相対 0.5 %) を全点で課す。**通らない点は捨てて REPAIRS から補う**。
2. `REPAIRS` は**ページ画像を目視して読んだ値**で、各行を同じ 2.54 則で検算してから入れてある。
3. 点数が揃わなければ**書き出さずに落ちる** (欠測のまま形状を作らない)。

出力:
  case/53.c3x_vane_cht/ref/vane_c3x.csv       78 点 (x_cm, y_cm)
  case/54.markii_vane_cht/ref/vane_markii.csv 60 点
  case/53.c3x_vane_cht/ref/test_conditions.csv 表 VIII/IX の主要行

usage: python3 case/53.c3x_vane_cht/tools/extract_vane_data.py [--pdf <path>] [--plot]
"""
import argparse
import math
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

# ---- 目視で読んだ補修値 (OCR が壊した行だけ)。各行 cm == 2.54*in を満たすことを確認済み ----
REPAIRS = {
    "markii": {   # page 11 (PDF 17)
        7:  (2.9812, 11.0833),   # (1.1737) (4.3635)
        22: (5.7269, 3.9444),    # (2.2547) (1.5529)
        24: (6.0295, 2.9741),    # (2.3738) (1.1709)
        27: (6.4554, 1.5128),    # (2.5415) (0.5956)
        57: (1.6337, 8.9891),    # (0.6432) (3.5390)
        59: (1.0208, 9.5456),    # (0.4019) (3.7581)
    },
    "c3x": {      # page 12 (PDF 18)
        # **報告自身が不整合**: 29 行は "0.4115 (0.0162)" と印字されているが 0.0162 in = 0.0411 cm。
        # 幾何で決める: 28 (7.6624, 0.6391) と 30 (7.8161, -0.0053) の間隔は ~0.66 cm で、
        # 0.4115 なら 28 から 0.28 / 30 から 0.42 と等間隔に並ぶ。0.0411 だと 30 から 0.047 cm しか離れず
        # 点列が重なる。→ **cm 値 0.4115 を採り、inch 側 (0.1620 であるべき) の誤植とみなす**。
        3:  (0.7658, 12.6764),   # (0.3015) (4.9907)  ← ヘッダの R_LE/R_TE に潰されていた行
        29: (7.8115, 0.4115),
        11: (4.2885, 10.5766),   # (1.6884) (4.1640)
        65: (1.5519, 9.3932),    # (0.6110) (3.6981)
        70: (0.2822, 10.4094),   # (0.1111) (4.0982)
        75: (0.0064, 11.0802),   # (0.0025) (4.3623)
    },
}
NPTS = {"markii": 60, "c3x": 78}
PAGES = {"markii": 17, "c3x": 18}


def norm(t):
    t = t.replace('"', '.').replace('·', '.')
    t = re.sub(r'(?<=\d),(?=\d)', '.', t)
    t = re.sub(r'\bi0\b', '10', t)
    t = t.replace('II.', '11.').replace('(I.', '(1.').replace('(i.', '(1.')
    return t


def extract_points(doc, name):
    body = norm(doc[PAGES[name]].get_text())
    # **表ヘッダを捨てる**: "R_LE = 1.168 cm (0.460 in.) R_TE = 0.173 cm (0.068 in.)" は
    # cm == 2.54*in を満たすので、そのままだとデータ行として通ってしまう (実際に点 3 を潰した)。
    # 列見出し ("Position number" が 2 回) の後ろだけを読む。
    k = body.rfind("number")
    if k > 0:
        body = body[k:]
    nums = re.findall(r'-?\d+\.\d+|-?\d+', body)
    # 位置番号は 1..N で重複なし、かつ cm == 2.54*in を課す。番号の逐次チェックまでやると
    # OCR が 1 行壊したときに同期を失うので、ここは**重複禁止**までに留め、欠測は REPAIRS で埋める。
    pts, i = {}, 0
    while i + 4 < len(nums):
        try:
            pos = int(float(nums[i])); xc = float(nums[i + 1]); xi = float(nums[i + 2])
            yc = float(nums[i + 3]); yi = float(nums[i + 4])
        except ValueError:
            i += 1; continue
        okx = abs(xc - 2.54 * xi) <= max(0.01, 0.005 * abs(xc))
        oky = abs(yc - 2.54 * yi) <= max(0.01, 0.005 * abs(yc))
        if 1 <= pos <= NPTS[name] and okx and oky and pos not in pts:
            pts[pos] = (xc, yc); i += 5
        else:
            i += 1
    n_ocr = len(pts)
    for pos, xy in REPAIRS[name].items():
        pts.setdefault(pos, xy)
    missing = sorted(set(range(1, NPTS[name] + 1)) - set(pts))
    print(f"[{name}] OCR validated {n_ocr}/{NPTS[name]}, repaired {len(REPAIRS[name])}, missing {missing}")
    if missing:
        sys.exit(f"[{name}] {len(missing)} points missing — ページ画像を読んで REPAIRS に足すこと (推測で埋めない)")
    return [pts[i] for i in range(1, NPTS[name] + 1)]


def write_csv(path, pts, header):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(header)
        f.write("i,x_cm,y_cm\n")
        for i, (x, y) in enumerate(pts, 1):
            f.write(f"{i},{x:.4f},{y:.4f}\n")
    print(f"  -> {path.relative_to(ROOT)} ({len(pts)} points)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", default=str(ROOT / "papers/cht/NASA-CR-168015_Hylton_1983.pdf"))
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    import fitz
    doc = fitz.open(a.pdf)

    meta = {
        "c3x": ("# C3X vane coordinates — NASA CR-168015 Table III (report p.12)\n"
                "# R_LE = 1.168 cm (0.460 in), R_TE = 0.173 cm (0.068 in)\n"
                "# 検証: 全点で cm == 2.54*in (相対 0.5%)。点 29 は報告の inch 値が誤植 (抽出スクリプト参照)\n"),
        "markii": ("# Mark II vane coordinates — NASA CR-168015 Table II (report p.11)\n"
                   "# R_LE = 1.280 cm (0.504 in), R_TE = 0.000 (blunt)\n"
                   "# 検証: 全点で cm == 2.54*in (相対 0.5%)\n"),
    }
    out = {"c3x": ROOT / "case/53.c3x_vane_cht/ref/vane_c3x.csv",
           "markii": ROOT / "case/54.markii_vane_cht/ref/vane_markii.csv"}
    pts = {}
    for name in ("c3x", "markii"):
        pts[name] = extract_points(doc, name)
        write_csv(out[name], pts[name], meta[name])

    # 試験条件 (表 VIII/IX, report p.30 = PDF 36) — **ページ画像から読んだ値**。
    # 注意: 表の SI 列 ("PT1--Pa") は psia 列と 51.7 倍ずれており内部矛盾している。
    #       psia 側が文献の引用値 (C3X 4411 で ~3.2 atm) と合うので **psia を正**とする。
    cond = ROOT / "case/53.c3x_vane_cht/ref/test_conditions.csv"
    with open(cond, "w") as f:
        f.write("# NASA CR-168015 Table VIII (Mark II) / IX (C3X), report p.30 — ページ画像から読取\n"
                "# PT1 は psia 列を正とする (表の SI 列は psia と 51.7 倍不整合)\n"
                "vane,code,run,PT1_psia,PT1_kPa,TT1_K,M1,Re1_e6,M2,Re2_e6,Tu_pct,Tw_over_Tg\n"
                "markii,5411,42,48.89,337.1,788,0.19,0.56,1.04,2.01,6.5,0.68\n"
                "markii,4411,43,49.64,342.3,784,0.18,0.57,0.89,1.98,6.5,0.69\n"
                "c3x,4411,108,46.34,319.5,786,0.17,0.52,0.90,1.99,6.5,0.73\n"
                "c3x,5411,107,45.24,311.9,798,0.17,0.51,1.05,1.97,6.5,0.72\n")
    print(f"  -> {cond.relative_to(ROOT)}")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 5))
        for k, name in enumerate(("c3x", "markii")):
            p = pts[name]
            ax[k].plot([q[0] for q in p] + [p[0][0]], [q[1] for q in p] + [p[0][1]], "-o", ms=2.5)
            ax[k].set_title(f"{name} ({len(p)} pts)"); ax[k].set_aspect("equal"); ax[k].grid(alpha=.3)
            ax[k].set_xlabel("x [cm]"); ax[k].set_ylabel("y [cm]")
        fig.tight_layout()
        png = ROOT / "case/53.c3x_vane_cht/ref/vane_shapes.png"
        fig.savefig(png, dpi=130)
        print(f"  -> {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
