#!/usr/bin/env python3
"""基準温度の取り方を変えた熱伝達率の比較 (plan §4.14)。

同じ壁面熱流束 q'' に対し、基準温度 T_ref の取り方だけを変えて h = q''/(T_ref - Tw) を出す。

  (A) すきま中央面の局所 T0   … `cavity_eval.py` の `h_eff`。面上の局所 ΔT で重み付けした
                                h_eff = Σq·w / Σ(ΔT·w)。分布を使う一番素直な定義。
  (B) キャビティ内の体積平均 T … ざっくり「中の平均ガス温度」を 1 点で代表させる定義。
  (C) 主流総温 Tt_inf         … 外部流の総温。設計で最も保守的 (ΔT が最大 -> h が最小)。
  (D) 断熱壁温 Taw           … 教科書的な定義 (回復温度基準)。cavity_eval.py の `h_aw`。

(B) の平均は **primal セル体積で重み付け**する。node 方式では近壁にノードが密集するので、
ノードの単純平均を取ると壁温側に大きく引きずられて「平均ガス温度」にならない。

usage: python3 tools/href_defs.py RUN [RUN ...]
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import geom_common as gc  # noqa: E402

# XDMF コード -> (節点数, 四面体分割)。全ヘキサ生産メッシュを主対象にする。
TETRA_OF = {
    6: (4, [(0, 1, 2, 3)]),                                                    # tetra
    9: (8, [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6), (1, 4, 5, 6), (3, 4, 6, 7)]),  # hexa
    8: (6, [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]),                        # wedge
    7: (5, [(0, 1, 2, 4), (0, 2, 3, 4)]),                                      # pyramid
}
CAV_WALLS = [("cav_outer", "外筒内壁"), ("cyl_side", "円柱側面"), ("cav_floor", "底面")]


def cells_uniform(conne, ncells):
    """全セルが同一型なら (idx[ncells,nn], code) を返す。違えば None (混在は未対応)。"""
    code = int(conne[0])
    if code not in TETRA_OF:
        return None
    nn = TETRA_OF[code][0]
    if len(conne) != ncells * (nn + 1):
        return None
    a = conne.reshape(ncells, nn + 1)
    if not np.all(a[:, 0] == code):
        return None
    return a[:, 1:].astype(np.int64), code


def tet_volumes(p, idx, tets):
    """四面体分割による**絶対値**体積 (向きは問わない)。"""
    vol = np.zeros(len(idx))
    for t in tets:
        a, b, c, d = (p[idx[:, t[k]]] for k in range(4))
        vol += np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0
    return vol


def cavity_mean_T(res, man):
    """キャビティ内 (z<0 のすきま) のセル体積重み平均 T と、比較用の単純ノード平均。"""
    with h5py.File(res, "r") as f:
        coord = np.array(f["MESH/COORD"]).reshape(-1, 3).astype(np.float64)
        T = np.array(f["VALUE/T"]).astype(np.float64)
        conne = np.array(f["MESH/CONNE"])
    nn_guess = TETRA_OF[int(conne[0])][0] + 1
    ncells = len(conne) // nn_guess
    got = cells_uniform(conne, ncells)
    if got is None:
        raise SystemExit("セル型が混在しているメッシュには未対応 (全ヘキサ想定)")
    idx, code = got
    tets = TETRA_OF[code][1]
    cent = coord[idx].mean(axis=1)
    m = gc.cavity_mask(cent[:, 0], cent[:, 1], cent[:, 2], man)
    if not m.any():
        raise SystemExit("キャビティ内のセルが 0 個 (manifest と res の幾何が不一致)")
    idx = idx[m]
    vol = tet_volumes(coord, idx, tets)
    Tc = T[idx].mean(axis=1)
    nm = gc.cavity_mask(coord[:, 0], coord[:, 1], coord[:, 2], man)
    return (float((vol * Tc).sum() / vol.sum()), float(vol.sum()),
            int(len(idx)), float(T[nm].mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()

    rows = []
    for r in a.runs:
        run = Path(r)
        man = gc.load_manifest(run=run)
        D = json.loads((run / "conditions.json").read_text())
        ev = json.loads((run / "cavity_eval.json").read_text())
        snaps = sorted(run.glob("res_[0-9]*.h5"),
                       key=lambda p: int(p.stem.split("_")[1]))
        Tcav, vol, ncav, Tnode = cavity_mean_T(snaps[-1], man)
        tp = D.get("gas_used") == "TP"
        Tt = float(D["Tt_tp"] if tp else D["Tt_cpg"])
        Taw = float(D["Taw_tp"] if tp else D["Taw_cpg"])
        rows.append(dict(run=str(run), Tw=float(D["wall_T"]), Tcav=Tcav, Tnode=Tnode,
                         vol=vol, ncav=ncav, Tt=Tt, Taw=Taw, ev=ev, step=snaps[-1].stem))

    for w in rows:
        print(f"{w['run']}  ({w['step']})")
        print(f"  壁温 Tw              {w['Tw']:9.2f} K ({w['Tw'] - 273.15:.0f} degC)")
        print(f"  キャビティ体積平均 T   {w['Tcav']:9.2f} K   "
              f"(セル {w['ncav']}, 体積 {w['vol'] * 1e6:.1f} cm3)")
        print(f"   参考: ノード単純平均  {w['Tnode']:9.2f} K   <- 近壁密集で壁温寄り。使わない")
        print(f"  主流総温 Tt_inf       {w['Tt']:9.2f} K   断熱壁温 Taw {w['Taw']:.2f} K")

    hdr = "".join("%12s" % (f"{r['Tw'] - 273.15:.0f} degC") for r in rows)
    print("\n" + " " * 30 + hdr)
    print("-" * (30 + 12 * len(rows)))
    for key, jp in CAV_WALLS:
        print(f"[{jp}]  q''_mean [kW/m2]".ljust(30)
              + "".join("%12.4g" % (r["ev"]["wall"][key]["qpp_mean"] / 1e3) for r in rows))
        for tag, dT in (("(A) すきま中央面 局所T0", None),
                        ("(B) キャビティ平均 T", "cav"),
                        ("(C) 主流総温 Tt", "tt"),
                        ("(D) 断熱壁温 Taw", "aw")):
            vals = []
            for r in rows:
                W = r["ev"]["wall"][key]
                if dT is None:
                    vals.append(W["h_eff"])
                elif dT == "aw":
                    vals.append(W["h_aw"])
                else:
                    d = (r["Tcav"] if dT == "cav" else r["Tt"]) - r["Tw"]
                    vals.append(W["qpp_mean"] / d if abs(d) > 1e-9 else float("nan"))
            print(("   h " + tag).ljust(30) + "".join("%12.4g" % v for v in vals))
        print()

    print("基準温度差 ΔT = T_ref - Tw [K]".ljust(30)
          + "".join("%12s" % "" for _ in rows))
    for tag, f in (("(B) キャビティ平均", lambda r: r["Tcav"] - r["Tw"]),
                   ("(C) 主流総温", lambda r: r["Tt"] - r["Tw"]),
                   ("(D) 断熱壁温", lambda r: r["Taw"] - r["Tw"])):
        print(("   " + tag).ljust(30) + "".join("%12.4g" % f(r) for r in rows))

    print("\n壁温に対する h の振れ幅 (最大/最小 - 1)  ** 小さいほど外挿に使える **")
    for key, jp in CAV_WALLS:
        out = []
        for tag, g in (("A", lambda W, r: W["h_eff"]),
                       ("B", lambda W, r: W["qpp_mean"] / (r["Tcav"] - r["Tw"])),
                       ("C", lambda W, r: W["qpp_mean"] / (r["Tt"] - r["Tw"])),
                       ("D", lambda W, r: W["h_aw"])):
            v = [g(r["ev"]["wall"][key], r) for r in rows]
            out.append("%s %+6.0f %%" % (tag, 100.0 * (max(v) / min(v) - 1.0)))
        print("   " + jp.ljust(10) + "   ".join(out))

    if a.plot:
        plot(rows)


def plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in Path.home().glob(".fonts/NotoSansCJK*"):
        try:
            font_manager.fontManager.addfont(str(f))
        except Exception:
            pass
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    tw = [r["Tw"] - 273.15 for r in rows]
    fig, ax = plt.subplots(1, 3, figsize=(14.0, 4.4), sharex=True)
    defs = [("(A) すきま中央面 局所T0", lambda W, r: W["h_eff"], "tab:blue", "o"),
            ("(B) キャビティ体積平均 T", lambda W, r: W["qpp_mean"] / (r["Tcav"] - r["Tw"]),
             "tab:orange", "s"),
            ("(C) 主流総温 Tt", lambda W, r: W["qpp_mean"] / (r["Tt"] - r["Tw"]),
             "tab:red", "^"),
            ("(D) 断熱壁温 Taw", lambda W, r: W["h_aw"], "tab:green", "v")]
    for b, (key, jp) in zip(ax, CAV_WALLS):
        for lab, g, col, mk in defs:
            b.plot(tw, [g(r["ev"]["wall"][key], r) for r in rows], mk + "-",
                   color=col, label=lab)
        b.set_yscale("log")
        b.set_title(jp, fontsize=12)
        b.set_xlabel("壁温 [°C]")
        b.grid(alpha=0.3, which="both")
    ax[0].set_ylabel("熱伝達率 h [W/m²K]")
    ax[0].legend(fontsize=9)
    fig.suptitle("case/49 形状2 同心 M9 — 基準温度の取り方による熱伝達率の違い "
                 "(同じ壁面熱流束を割る温度差だけを変えたもの)", fontsize=12.5)
    fig.tight_layout()
    out = Path(rows[0]["run"]).parent / "href_defs.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
