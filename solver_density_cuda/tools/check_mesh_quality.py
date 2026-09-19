#!/usr/bin/env python3
"""forge 入力 H5 メッシュの品質指標 (アスペクト比・スキューネス) を計算・判定する。

使い方:
    python3 check_mesh_quality.py mesh.h5 [--ar-max 1000] [--skew-max 0.9] [--mode auto|2d|3d]

- アスペクト比 (AR): 各セルの 最長辺/最短辺。薄い境界層セルの監視用。
- スキューネス (equiangle skew): 内角 θ から max[(θmax-ideal)/(180-ideal),(ideal-θmin)/ideal]
  (0=直交/正三角、1=退化)。四角形 ideal=90、三角形 ideal=60。
- 既定しきい値: AR ≤ 1000, skew ≤ 0.9 (本リポジトリの計算前メッシュ品質ルール、
  AGENTS.md / forge-calculation-workflow.md)。VERDICT を PASS/FAIL で返す。

2D / 3D の扱い (重要):
- **2D モード** (準2D・z 押し出し・軸対称): セルを x-y 平面に射影し多角形として測る (従来挙動)。
- **3D モード** (完全非構造 3D: tet/prism/hex/pyramid が z 方向にも構造を持つ): x-y 射影は
  無意味 (z 平行エッジが一点に潰れ AR=inf 等の偽値) なため、**各セルの面ごとに equiangle skew
  を測りその最大**を skew、**全辺の最長/最短**を AR とする (真の 3D 形状で判定)。
- `--mode auto` (既定): 3D セル型を含み z 平面が 3 枚以上なら 3D、それ以外は 2D と自動判定。

`solver_density_cuda/tools/res_h5_to_vtu.py` の parse_conne を流用して CONNE を読む。
node (median-dual) 変換の h5 は `/MESH/CONNE` が双対 CV (多角形/多面体) なので、`/VIZMESH/CONNE` (primal
セル) があればそちらを測る (品質は primal の性質。cell 変換をやり直す必要はない)。
計算は VTK セル型ごとに numpy でベクトル化 (2026-09-08: 2M セル hex で 15 分 → 数秒)。`--loop` で旧の
セル毎ループ (検証用) を使う。
"""
import sys, os, argparse
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from res_h5_to_vtu import parse_conne

# VTK cell type -> 面 (セル局所節点 index)。parse_conne の節点順 (gmsh/XDMF/VTK) 準拠。
# TETRA=10, HEXAHEDRON=12, WEDGE(prism)=13, PYRAMID=14
FACES_3D = {
    10: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
    12: [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
    13: [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    14: [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}
VTK_3D = set(FACES_3D)


def poly_skew(pts):
    """多角形 (n,2 or n,3) の equiangle skew を返す (cell_metrics の skew 部と同一定義)。"""
    n = len(pts)
    angs = []
    for i in range(n):
        a = pts[(i - 1) % n] - pts[i]
        b = pts[(i + 1) % n] - pts[i]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na == 0 or nb == 0:
            angs.append(0.0); continue
        c = np.clip(np.dot(a, b) / (na * nb), -1.0, 1.0)
        angs.append(np.degrees(np.arccos(c)))
    angs = np.array(angs)
    ideal = 90.0 if n == 4 else (60.0 if n == 3 else 180.0 * (n - 2) / n)
    return max((angs.max() - ideal) / (180.0 - ideal), (ideal - angs.min()) / ideal)


def cell_metrics(pts):
    """pts: (n,2) 多角形頂点 (CCW or CW)。(aspect_ratio, skew) を返す (2D モード用)。"""
    n = len(pts)
    edges = [np.linalg.norm(pts[(i + 1) % n] - pts[i]) for i in range(n)]
    emin, emax = min(edges), max(edges)
    ar = emax / emin if emin > 0 else float("inf")
    return ar, poly_skew(pts)


# VTK セル型ごとの四面体分割 (符号付き体積の合計を取るため)
TETRA_OF = {
    10: [(0, 1, 2, 3)],                                              # tetra
    12: [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6), (1, 4, 5, 6), (3, 4, 6, 7)],  # hexa
    13: [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)],                  # wedge (prism)
    14: [(0, 1, 2, 4), (0, 2, 3, 4)],                                # pyramid
}


def cell_volumes(coord3, conn, offs, types):
    """セルの符号付き体積 (四面体分割の和)。型が未知なら NaN。**成立性検査専用**。"""
    vols = np.full(len(offs), np.nan)
    s = 0
    for i, o in enumerate(offs):
        idx = conn[s:o]; s = o
        tets = TETRA_OF.get(int(types[i]))
        if tets is None or len(idx) < 4:
            continue
        pts = coord3[idx]
        v = 0.0
        okc = True
        for t in tets:
            if max(t) >= len(pts):
                okc = False
                break
            p0, p1, p2, p3 = pts[t[0]], pts[t[1]], pts[t[2]], pts[t[3]]
            v += np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0
        vols[i] = v if okc else np.nan
    return vols


def cell_metrics_3d(pts, vtk):
    """pts: (m,3) セル節点。面ごとの equiangle skew の最大と、全辺の 最長/最短 を返す。"""
    faces = FACES_3D[vtk]
    skew = 0.0
    edges = []
    for fc in faces:
        fp = pts[list(fc)]
        skew = max(skew, poly_skew(fp))
        for i in range(len(fc)):
            edges.append(np.linalg.norm(fp[(i + 1) % len(fc)] - fp[i]))
    emin, emax = min(edges), max(edges)
    ar = emax / emin if emin > 0 else float("inf")
    return ar, skew


def _poly_skew_batch(P):
    """P: (nc, n, d) 多角形頂点列。各多角形の equiangle skew (nc,) をベクトル計算 (poly_skew と同一定義)。"""
    n = P.shape[1]
    a = np.roll(P, 1, axis=1) - P          # pts[i-1] - pts[i]
    b = np.roll(P, -1, axis=1) - P         # pts[i+1] - pts[i]
    na = np.linalg.norm(a, axis=2); nb = np.linalg.norm(b, axis=2)
    den = na * nb
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(den > 0, np.einsum("ijk,ijk->ij", a, b) / np.where(den > 0, den, 1.0), 1.0)
    ang = np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
    ang = np.where(den > 0, ang, 0.0)      # 退化辺 (長さ 0) は角 0 (poly_skew と同じ)
    ideal = 90.0 if n == 4 else (60.0 if n == 3 else 180.0 * (n - 2) / n)
    return np.maximum((ang.max(axis=1) - ideal) / (180.0 - ideal), (ideal - ang.min(axis=1)) / ideal)


def _edge_ar_batch(P):
    """P: (nc, n, d)。各多角形の辺長 (nc, n) を返す。"""
    return np.linalg.norm(np.roll(P, -1, axis=1) - P, axis=2)


def metrics_vectorized(coord, conn, offs, types, is3d):
    """セル型ごとにまとめて (ars, skews) を計算する (cell_metrics / cell_metrics_3d と同一定義)。"""
    ncells = len(offs)
    ars = np.empty(ncells); skews = np.empty(ncells)
    starts = np.concatenate(([0], offs[:-1]))
    sizes = offs - starts
    for t in np.unique(types):
        sel = np.where(types == t)[0]
        nn = np.unique(sizes[sel])
        for m in nn:                                   # 同型でも節点数が違えば分ける (安全側)
            idx = sel[sizes[sel] == m]
            gather = starts[idx][:, None] + np.arange(m)[None, :]
            P = coord[conn[gather]]                     # (nc, m, d)
            if not is3d:
                e = _edge_ar_batch(P)
                emin = e.min(axis=1); emax = e.max(axis=1)
                ars[idx] = np.where(emin > 0, emax / np.where(emin > 0, emin, 1.0), np.inf)
                skews[idx] = _poly_skew_batch(P)
            else:
                faces = FACES_3D.get(int(t))
                if faces is None:
                    raise SystemExit(f"unsupported 3D VTK type {t}")
                sk = np.zeros(len(idx)); emin = np.full(len(idx), np.inf); emax = np.zeros(len(idx))
                for fc in faces:
                    fp = P[:, list(fc), :]
                    sk = np.maximum(sk, _poly_skew_batch(fp))
                    e = _edge_ar_batch(fp)
                    emin = np.minimum(emin, e.min(axis=1)); emax = np.maximum(emax, e.max(axis=1))
                ars[idx] = np.where(emin > 0, emax / np.where(emin > 0, emin, 1.0), np.inf)
                skews[idx] = sk
    return ars, skews


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--ar-max", type=float, default=1000.0)
    ap.add_argument("--skew-max", type=float, default=0.9)
    ap.add_argument("--mode", choices=["auto", "2d", "3d"], default="auto")
    ap.add_argument("--loop", action="store_true", help="旧のセル毎ループで計算 (検証用, 遅い)")
    ap.add_argument("--no-vizmesh", action="store_true", help="/VIZMESH があっても /MESH/CONNE を測る")
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        coord3 = np.array(f["MESH/COORD"]).reshape(-1, 3)
        if "VIZMESH/CONNE" in f and not a.no_vizmesh:
            # node (median-dual) 変換: primal セル (可視化トポロジ) を測る
            conne = np.array(f["VIZMESH/CONNE"]); ncells = int(f["VIZMESH"].attrs["nVizCells"])
            print(f"source: /VIZMESH/CONNE (primal cells of a node-converted mesh, {ncells} cells)")
        else:
            conne = np.array(f["MESH/CONNE"])
            ncells = f["VALUE/ro"].shape[0] if "VALUE/ro" in f else None
    if ncells is None:
        # count from CONNE
        from res_h5_to_vtu import _count_cells
        ncells = _count_cells(conne)
    conn, offs, types = parse_conne(conne, ncells)

    # モード判定: 3D セル型を含み z 平面が 3 枚以上なら 3D (x-y 射影が破綻するため)
    if a.mode == "auto":
        nz = len(np.unique(np.round(coord3[:, 2], 6)))
        is3d = bool(np.isin(types, list(VTK_3D)).any()) and nz > 2
    else:
        is3d = (a.mode == "3d")
    print(f"mode: {'3D (per-face skew / edge AR)' if is3d else '2D (x-y projection)'}")

    coord = coord3 if is3d else coord3[:, :2]
    if a.loop:
        ars = np.empty(ncells); skews = np.empty(ncells)
        s = 0
        for i, o in enumerate(offs):
            pts = coord[conn[s:o]]; s = o
            ars[i], skews[i] = (cell_metrics_3d(pts, types[i]) if is3d else cell_metrics(pts))
    else:
        ars, skews = metrics_vectorized(coord, np.asarray(conn), np.asarray(offs), np.asarray(types), is3d)

    # --- メッシュの**成立性**を、形状品質の外れ値許容より前に無条件で検査する ---
    # (2026-09-19 codex: NaN 座標が 1 セルあっても `SOFT-PASS` exit 0、体積ゼロの四面体は
    #  AR 1.732 / skew 0.500 で通常合格していた。形状の良し悪し以前の問題は許容しない。)
    fatal = []
    nonfinite = int(np.count_nonzero(~np.isfinite(coord3)))
    if nonfinite:
        bad_nodes = int(np.count_nonzero(~np.isfinite(coord3).all(axis=1)))
        fatal.append("非有限座標: 成分 %d 個 / 節点 %d 個" % (nonfinite, bad_nodes))
    if is3d:
        vols = cell_volumes(coord3, np.asarray(conn), np.asarray(offs), np.asarray(types))
        nv_bad = int(np.count_nonzero(~np.isfinite(vols)))
        nz_bad = int(np.count_nonzero(np.abs(vols) <= 0.0))
        neg = int(np.count_nonzero(vols < 0.0))
        pos = int(np.count_nonzero(vols > 0.0))
        if nv_bad:
            fatal.append("体積が非有限のセル %d 個" % nv_bad)
        if nz_bad:
            fatal.append("体積ゼロのセル %d 個" % nz_bad)
        if neg and pos:
            fatal.append("体積の符号が混在 (正 %d / 負 %d) = 向きの不整合" % (pos, neg))
        good = vols[np.isfinite(vols) & (vols != 0.0)]
        if good.size:
            print("cell volume  : |V| min=%.3e  median=%.3e  max=%.3e  (0 or 非有限: %d)"
                  % (np.abs(good).min(), np.median(np.abs(good)), np.abs(good).max(),
                     nv_bad + nz_bad))
    if fatal:
        print("FATAL: " + " / ".join(fatal))
        print("VERDICT: FAIL (メッシュとして成立していない。形状品質の外れ値許容は適用しない)")
        sys.exit(1)

    ar_bad = int((ars > a.ar_max).sum())
    sk_bad = int((skews > a.skew_max).sum())
    print(f"cells: {ncells}")
    print(f"aspect ratio : max={ars.max():.1f}  p99={np.percentile(ars,99):.1f}  "
          f"mean={ars.mean():.1f}  >{a.ar_max:.0f}: {ar_bad} ({100*ar_bad/ncells:.2f}%)")
    print(f"skewness     : max={skews.max():.3f}  p99={np.percentile(skews,99):.3f}  "
          f"mean={skews.mean():.3f}  >{a.skew_max:.2f}: {sk_bad} ({100*sk_bad/ncells:.2f}%)")
    ok = (ars.max() <= a.ar_max) and (skews.max() <= a.skew_max)
    # 少数の外れ値は許容するソフト判定も併記
    soft = (ar_bad <= 0.001 * ncells) and (sk_bad <= 0.001 * ncells)
    print(f"VERDICT: {'PASS' if ok else ('SOFT-PASS (<0.1% outliers)' if soft else 'FAIL')} "
          f"(AR<= {a.ar_max:.0f}, skew<= {a.skew_max:.2f})")
    sys.exit(0 if (ok or soft) else 1)


if __name__ == "__main__":
    main()
