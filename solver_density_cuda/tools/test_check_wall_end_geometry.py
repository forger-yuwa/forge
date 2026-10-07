#!/usr/bin/env python3
"""check_wall_end_geometry.py (壁終端の格子の幾何の点検) の単体試験。GPU 不要。

    python3 solver_density_cuda/tools/test_check_wall_end_geometry.py

合成の小さい node 変換 h5 (MESH/COORD・PLANES/STRUCT・VIZMESH/CONNE・BCONDS の vizBface) で:
  (1) 厚さ 0 の平板が途中で終わり、下流の格子線を角度 θ (0° と 40°) で折った格子 — 3D hex (span 3 節点)・
      疑似 2D (1 層押し出し、span 2 節点)・2D quad の 3 種。期待値は格子の作り方から解析的に出す:
        kink = θ、len = (Δx_d / cos θ) / Δx_u、disp/h = ∓Δx_d tan θ / (y_{k+1} − y_k) (上側 −、下側 +)、
        LSQ の固有値 = 既知の近傍の方向から直接組んだ M の固有値。全層 (0〜3) と全 span の節点で確かめる。
      端の分類 (knife が 1 本だけ、前縁・span 端の辺は端にしない)、--pair の共有数、扇のセル数 4。
  (2) 段差 (凸の角 270°、壁 2 つ) → convex。床の面は下流の辺がまっすぐ (kink 0)、ベースの面は格子線が上へ延びる。
  (3) 壁が同じ面内で出口の境界に続く → open。壁でない面は解析しない。
  (4) 入力の拒否: タグ誤記・共有なしの --pair・未対応のセル型・node 変換でない h5、CLI の終了コードと CSV、--help の定義と限界。
  (5) 板の後縁の線と露出した側端の線の角 (SERN のカウル後縁と側端の角と同じ接続: 角の節点で上下の面が 3 節点を共有し
      扇が閉じる)。2 本の knife の線に分かれ、角の節点も含めて flag なしで後縁は θ、側端は 0° の折れと期待の辺長比。
  (6) ランプの膨張の角のような 45° 未満の壁の折れ (φ 約 200°) は既定では終端にしない (--convex-deg を下げれば出る)。
"""
import math
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_wall_end_geometry as cw  # noqa: E402

TOOL = os.path.join(HERE, "check_wall_end_geometry.py")
fail = 0


def check(name, ok, detail=""):
    global fail
    print(("ok   " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    fail += 0 if ok else 1


# ---------------------------------------------------------------------------
# 合成の node 変換 h5 (境界面はセルの面のうち 1 回しか現れない面、タグは tagger で決める)
# ---------------------------------------------------------------------------
HEX_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))


def write_mesh(path, coord, cells, tagger, kinds, coord_dtype=np.float64):
    """cells: (nc, 8) hex または (nc, 4) quad。tagger(面の重心, セルの重心) → タグ名。kinds = {名前: (physID, kind)}。"""
    coord = np.asarray(coord, dtype=float)
    cells = np.asarray(cells, dtype=np.int64)
    hexa = cells.shape[1] == 8
    faces = {}
    if hexa:
        for c, row in enumerate(cells):
            for fl in HEX_FACES:
                f = tuple(int(row[i]) for i in fl)
                faces.setdefault(tuple(sorted(f)), []).append((f, c))
        edge_loc = cw.HEX_EDGES
    else:
        for c, row in enumerate(cells):
            for i in range(4):
                f = (int(row[i]), int(row[(i + 1) % 4]))
                faces.setdefault(tuple(sorted(f)), []).append((f, c))
        edge_loc = ((0, 1), (1, 2), (2, 3), (3, 0))
    bf = {nm: [] for nm in kinds}
    for occ in faces.values():
        if len(occ) != 1:
            continue
        f, c = occ[0]
        nm = tagger(coord[list(f)].mean(axis=0), coord[cells[c]].mean(axis=0))
        bf[nm].append(f)
    edges = set()
    for row in cells:
        for a, b in edge_loc:
            edges.add((min(row[a], row[b]), max(row[a], row[b])))
    edges = sorted(edges)
    n = coord.shape[0]
    with h5py.File(path, "w") as f:
        g = f.create_group("MESH")
        g.attrs["nNodes"] = n; g.attrs["nCells"] = n
        g.attrs["nNormalPlanes"] = len(edges); g.attrs["nPlanes"] = len(edges)
        g.attrs["nBconds"] = len(kinds); g.attrs["nBPlanes"] = 0
        g.create_dataset("COORD", data=coord.astype(coord_dtype).reshape(-1))
        f.create_dataset("PLANES/STRUCT", data=np.array([[2, a, b, 2, a, b] for a, b in edges], dtype=np.int64).reshape(-1))
        code = 9 if hexa else 5
        conne = np.concatenate([np.full((cells.shape[0], 1), code), cells], axis=1).reshape(-1)
        v = f.create_group("VIZMESH"); v.attrs["nVizCells"] = cells.shape[0]; v.attrs["vizCONNE_dim"] = conne.size
        v.create_dataset("CONNE", data=conne)
        for nm, (pid, kind) in kinds.items():
            b = f.create_group(f"BCONDS/{pid}")
            b.attrs["bcondKind"] = kind
            fs = bf[nm]
            b.create_dataset("vizBfaceSizes", data=np.array([len(x) for x in fs], dtype=np.int64))
            b.create_dataset("vizBfaceNodes", data=np.array([i for x in fs for i in x], dtype=np.int64))
    return bf


def write_bcond(path, kinds):
    with open(path, "w") as fh:
        for nm, (pid, kind) in kinds.items():
            fh.write(f"{nm}: {{physID: {pid}, kind: {kind}, outputHDFflg: 0, ints: , floats: }}\n")


def args_for(h5, bc=None, **kw):
    a = cw.build_parser().parse_args([h5] + (["--bcond-config", bc] if bc else []) + ["--quiet", "--no-csv"])
    for k, v in kw.items():
        setattr(a, k, v)
    return a


# ---------------------------------------------------------------------------
# (1) 厚さ 0 の平板の後縁
# ---------------------------------------------------------------------------
DXU, DXD, NU, ND = 0.1, 0.25, 3, 3
OFFS = [0.0, 0.01, 0.025, 0.05, 0.1, 0.2]          # 板からの層の位置 (上下対称)
DZ = 0.05
PLATE_KINDS = {"inlet": (1, "inlet_uniformVelocity"), "outlet": (2, "outflow"), "bottom": (3, "slip"), "top": (4, "slip"),
               "plate_up": (5, "wall_isothermal"), "plate_lo": (6, "wall"), "zmin": (7, "slip"), "zmax": (8, "slip")}


def plate(theta_deg, nz):
    """nz = 0 は 2D quad、nz ≥ 2 は span nz 節点の hex。後縁 x = 0、板は y = 0 (x ≤ 0)、下流 (x > 0) は全行を
    y を −x tan θ だけ剛体的にずらす (旧メッシャの中間線の折れと同じ形)。板の上面 (x < 0) は座標一致の別 ID。"""
    tn = math.tan(math.radians(theta_deg))
    xs = np.r_[-DXU * np.arange(NU, 0, -1), 0.0, DXD * np.arange(1, ND + 1)]
    ys = np.r_[-np.array(OFFS[::-1]), np.array(OFFS[1:])]
    JP = len(OFFS) - 1; ITE = NU
    NX, NY = len(xs), len(ys)
    NZ = max(nz, 1)
    zs = DZ * np.arange(NZ)

    def base(i, j, k):
        return (i * NY + j) * NZ + k
    coord = []
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                coord.append((xs[i], ys[j] - (xs[i] * tn if xs[i] > 0 else 0.0), zs[k]))
    dup = {}
    for i in range(ITE):
        for k in range(NZ):
            dup[(i, k)] = len(coord); coord.append((xs[i], 0.0, zs[k]))

    def nd(i, j, k, upper):
        return dup[(i, k)] if (upper and j == JP and (i, k) in dup) else base(i, j, k)
    cells = []
    for i in range(NX - 1):
        for j in range(NY - 1):
            up = j >= JP
            if nz == 0:
                cells.append((nd(i, j, 0, up), nd(i + 1, j, 0, up), nd(i + 1, j + 1, 0, up), nd(i, j + 1, 0, up)))
                continue
            for k in range(NZ - 1):
                cells.append((nd(i, j, k, up), nd(i + 1, j, k, up), nd(i + 1, j + 1, k, up), nd(i, j + 1, k, up),
                              nd(i, j, k + 1, up), nd(i + 1, j, k + 1, up), nd(i + 1, j + 1, k + 1, up), nd(i, j + 1, k + 1, up)))
    x0, x1 = xs[0], xs[-1]

    def tagger(fc, cc):
        if nz and abs(fc[2]) < 1e-12:
            return "zmin"
        if nz and abs(fc[2] - zs[-1]) < 1e-12:
            return "zmax"
        if abs(fc[0] - x0) < 1e-12:
            return "inlet"
        if abs(fc[0] - x1) < 1e-12:
            return "outlet"
        if fc[0] < 0 and abs(fc[1]) < 1e-12:
            return "plate_up" if cc[1] > 0 else "plate_lo"
        return "bottom" if fc[1] < cc[1] else "top"          # 面がセルの下なら下の境界
    info = {"base": base, "dup": dup, "JP": JP, "ITE": ITE, "NZ": NZ, "xs": xs, "tn": tn}
    return np.array(coord), cells, tagger, info


def plate_expect(theta_deg, side_up, k):
    """後縁の列の第 k 層の期待値 (kink, len, disp/h)。"""
    tn = math.tan(math.radians(theta_deg))
    h = OFFS[k + 1] - OFFS[k]
    disp = DXD * tn * (-1.0 if side_up else 1.0)
    return theta_deg, (DXD / math.cos(math.radians(theta_deg))) / DXU, disp / h


def plate_lsq_expect(theta_deg, nz, k_span, dim2):
    """後縁の節点 (層 0) の LSQ 近傍の方向 (格子の作り方から直接)。"""
    tn = math.tan(math.radians(theta_deg))
    vecs = [(-DXU, 0, 0), (-DXU, 0, 0), (0, OFFS[1], 0), (0, -OFFS[1], 0), (DXD, -DXD * tn, 0)]
    if not dim2:
        if k_span > 0:
            vecs.append((0, 0, -DZ))
        if k_span < nz - 1:
            vecs.append((0, 0, DZ))
    d = np.array(vecs, dtype=float)
    if dim2:
        d = d[:, :2]
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return np.linalg.eigvalsh(d.T @ d)


with tempfile.TemporaryDirectory() as td:
    P = lambda name: os.path.join(td, name)
    bc_plate = P("bcond_plate.yaml")
    write_bcond(bc_plate, PLATE_KINDS)

    for nz, label in ((3, "3D hex"), (2, "疑似 2D (1 層押し出し)"), (0, "2D quad")):
        dim2 = nz == 0
        for theta in (0.0, 40.0):
            c, cells, tagger, info = plate(theta, nz)
            h5 = P(f"plate_{nz}_{int(theta)}.h5")
            write_mesh(h5, c, cells, tagger, PLATE_KINDS)
            res = cw.run(args_for(h5, bc_plate, pair=["plate_up,plate_lo"]))
            tag = f"(1) {label} θ={theta:g}°"
            lines = res["lines"]
            check(f"{tag}: 終端線は knife 1 本 (plate_lo|plate_up)、前縁・span 端は端にしない",
                  len(lines) == 1 and lines[0]["cls"] == "knife" and lines[0]["tags"] == ("plate_lo", "plate_up"),
                  str([(L["cls"], L["tags"]) for L in lines]))
            te = sorted(info["base"](info["ITE"], info["JP"], k) for k in range(info["NZ"]))
            check(f"{tag}: 終端線の節点 = 後縁の列 (全 span)", lines and sorted(lines[0]["nodes"]) == te,
                  f"{lines[0]['nodes'] if lines else None} vs {te}")
            check(f"{tag}: 流体角 φ = 360° (厚さ 0)", lines and abs(lines[0]["phi"][0] - 360.0) < 1e-9 and abs(lines[0]["phi"][1] - 360.0) < 1e-9,
                  str(lines[0]["phi"] if lines else None))
            pr = res["pairs"][0]
            n_sh = 1 if dim2 else info["NZ"] - 1
            check(f"{tag}: --pair plate_up,plate_lo の共有 {n_sh}・すべて knife", pr["n_shared"] == n_sh and pr["n_knife"] == n_sh, str(pr))
            rows = res["rows"]
            check(f"{tag}: 行数 = 節点 × 2 面 × 4 層", len(rows) == len(te) * 2 * 4, str(len(rows)))
            bad = []
            for r in rows:
                up = r["side"] == "plate_up"
                ek, el, ed = plate_expect(theta, up, r["layer"])
                ok = (r["flag"] == "" and abs(r["kink_deg"] - ek) < 1e-9 and abs(r["len_ratio"] - el) < 1e-9
                      and abs(r["disp_ratio"] - ed) < 1e-6 * max(1.0, abs(ed)) and r["fan_cells"] == 4)
                if not ok:
                    bad.append((r["side"], r["layer"], r["node"], r.get("kink_deg"), r.get("len_ratio"), r.get("disp_ratio"), r["flag"]))
            check(f"{tag}: 全節点・両面・層 0〜3 で kink・len・disp/h が期待値", not bad, str(bad[:3]))
            # 層の節点が格子の作り方どおり (下側の第 k 層 = 板から k 行下、上側 = k 行上)
            ok_nodes = True
            for r in rows:
                kk = int(round(r["z"] / DZ)) if not dim2 else 0
                j = info["JP"] + (r["layer"] if r["side"] == "plate_up" else -r["layer"])
                ok_nodes &= r["node"] == info["base"](info["ITE"], j, kk)
            check(f"{tag}: 第 k 層の節点 = 後縁の列の板から k 行目", ok_nodes)
            # LSQ (層 0 の後縁の節点)
            bad = []
            for r in rows:
                if r["layer"] != 0:
                    continue
                kk = 0 if dim2 else int(round(r["z"] / DZ))
                lam_e = plate_lsq_expect(theta, info["NZ"], kk, dim2)
                if len(r["lam"]) != lam_e.size or np.max(np.abs(np.array(r["lam"]) - lam_e)) > 1e-12 or r["n_lsq"] != (5 if dim2 else 5 + (kk > 0) + (kk < info["NZ"] - 1)):
                    bad.append((kk, r["lam"], lam_e.tolist(), r["n_lsq"]))
            check(f"{tag}: 後縁の節点の LSQ 近傍の数と M の固有値", not bad, str(bad[:2]))
            if theta == 0.0 and nz == 3:
                r0 = [r for r in rows if r["layer"] == 0 and abs(r["z"] - DZ) < 1e-12][0]
                check(f"{tag}: 内部の span の後縁の M = diag(3, 2, 2) (u が上下 2 つ・z±)", np.allclose(r0["lam"], [2, 2, 3]), str(r0["lam"]))
            # 要約
            sm = [s for s in res["summary"] if s["layer"] == 0]
            check(f"{tag}: 要約は 2 面 (上下を分ける) × 層、kink 最大 = θ", len(res["summary"]) == 8 and len(sm) == 2
                  and all(abs(s["kink_max"] - theta) < 1e-9 for s in sm), str([(s["side"], s["kink_max"]) for s in sm]))

    # float32 座標 (変換器の既定) でも 40° を 1e-4° で拾う
    c, cells, tagger, info = plate(40.0, 3)
    write_mesh(P("plate_f32.h5"), c, cells, tagger, PLATE_KINDS, coord_dtype=np.float32)
    res = cw.run(args_for(P("plate_f32.h5"), bc_plate))
    check("(1) float32 座標: kink 40° (1e-4° 以内)", all(abs(r["kink_deg"] - 40.0) < 1e-4 for r in res["rows"]))

    # ---------------------------------------------------------------------------
    # (2) 段差 (凸の角): 床 (y = 0, x < 0) とベース (x = 0, y < 0) の 2 つの壁
    # ---------------------------------------------------------------------------
    STEP_KINDS = {"inlet": (1, "inlet_uniformVelocity"), "outlet": (2, "outflow"), "bottom": (3, "slip"), "top": (4, "slip"),
                  "floor": (5, "wall"), "base": (6, "wall"), "zmin": (7, "slip"), "zmax": (8, "slip")}
    xs = np.array([-0.3, -0.2, -0.1, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25])
    ys = np.array([-0.3, -0.15, -0.05, 0.0, 0.01, 0.03, 0.1, 0.3])
    zs = np.array([0.0, 0.05, 0.1])
    ids = -np.ones((len(xs), len(ys), len(zs)), dtype=np.int64)
    coord = []
    for i in range(len(xs)):
        for j in range(len(ys)):
            if xs[i] < 0 and ys[j] < 0:
                continue
            for k in range(len(zs)):
                ids[i, j, k] = len(coord); coord.append((xs[i], ys[j], zs[k]))
    cells = []
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            if xs[i] < 0 and ys[j] < 0:
                continue                       # 固体 (段の下)
            for k in range(len(zs) - 1):
                cells.append((ids[i, j, k], ids[i + 1, j, k], ids[i + 1, j + 1, k], ids[i, j + 1, k],
                              ids[i, j, k + 1], ids[i + 1, j, k + 1], ids[i + 1, j + 1, k + 1], ids[i, j + 1, k + 1]))

    def step_tag(fc, cc):
        if abs(fc[2]) < 1e-12: return "zmin"
        if abs(fc[2] - zs[-1]) < 1e-12: return "zmax"
        if abs(fc[0] - xs[0]) < 1e-12: return "inlet"
        if abs(fc[0] - xs[-1]) < 1e-12: return "outlet"
        if abs(fc[1] - ys[-1]) < 1e-12: return "top"
        if abs(fc[1] - ys[0]) < 1e-12: return "bottom"
        if abs(fc[1]) < 1e-12 and fc[0] < 0: return "floor"
        if abs(fc[0]) < 1e-12 and fc[1] < 0: return "base"
        raise AssertionError(f"未分類の境界面 {fc}")
    write_mesh(P("step.h5"), np.array(coord), cells, step_tag, STEP_KINDS)
    bc_step = P("bcond_step.yaml"); write_bcond(bc_step, STEP_KINDS)
    res = cw.run(args_for(P("step.h5"), bc_step))
    L = res["lines"]
    check("(2) 段差: 終端線は convex 1 本 (base|floor)、φ = 270°", len(L) == 1 and L[0]["cls"] == "convex" and L[0]["tags"] == ("base", "floor")
          and abs(L[0]["phi"][0] - 270) < 1e-9 and abs(L[0]["phi"][1] - 270) < 1e-9, str([(x["cls"], x["tags"], x["phi"]) for x in L]))
    fl = [r for r in res["rows"] if r["side"] == "floor"]
    bs = [r for r in res["rows"] if r["side"] == "base"]
    check("(2) 床の面: 下流の辺は x 方向 (kink 0、len 2.5、disp 0)、第 1 層 = 0.01 上",
          all(abs(r["kink_deg"]) < 1e-9 and abs(r["disp"]) < 1e-12 for r in fl)
          and all(abs(r["len_ratio"] - 2.5) < 1e-9 for r in fl if r["layer"] == 0)
          and all(abs(r["h_layer"] - 0.01) < 1e-12 for r in fl if r["layer"] == 0) and all(r["fan_cells"] == 3 and r["flag"] == "" for r in fl),
          str([(r["layer"], r.get("kink_deg"), r.get("len_ratio"), r.get("h_layer"), r["flag"]) for r in fl[:4]]))
    check("(2) ベースの面: 格子線は x = 0 を上へ (kink 0、len 0.01/0.05 = 0.2、第 1 層 = 下流の 0.25)",
          all(abs(r["kink_deg"]) < 1e-9 for r in bs) and all(abs(r["len_ratio"] - 0.2) < 1e-9 and abs(r["h_layer"] - 0.25) < 1e-12
                                                            for r in bs if r["layer"] == 0) and all(r["flag"] == "" for r in bs),
          str([(r["layer"], r.get("kink_deg"), r.get("len_ratio"), r.get("h_layer"), r["flag"]) for r in bs[:4]]))

    # ---------------------------------------------------------------------------
    # (3) 壁 → 出口の境界 (open)
    # ---------------------------------------------------------------------------
    OPEN_KINDS = {"inlet": (1, "inlet_uniformVelocity"), "outlet": (2, "outflow"), "top": (4, "slip"),
                  "floor": (5, "wall"), "floor_out": (6, "outflow"), "zmin": (7, "slip"), "zmax": (8, "slip")}
    xs = np.array([-0.3, -0.2, -0.1, 0.0, 0.25, 0.5]); ys = np.array([0.0, 0.01, 0.03, 0.1, 0.3]); zs = np.array([0.0, 0.05])
    ids = np.arange(len(xs) * len(ys) * len(zs)).reshape(len(xs), len(ys), len(zs))
    coord = [(x, y, z) for x in xs for y in ys for z in zs]
    cells = [(ids[i, j, 0], ids[i + 1, j, 0], ids[i + 1, j + 1, 0], ids[i, j + 1, 0], ids[i, j, 1], ids[i + 1, j, 1], ids[i + 1, j + 1, 1], ids[i, j + 1, 1])
             for i in range(len(xs) - 1) for j in range(len(ys) - 1)]

    def open_tag(fc, cc):
        if abs(fc[2]) < 1e-12: return "zmin"
        if abs(fc[2] - zs[-1]) < 1e-12: return "zmax"
        if abs(fc[0] - xs[0]) < 1e-12: return "inlet"
        if abs(fc[0] - xs[-1]) < 1e-12: return "outlet"
        if abs(fc[1] - ys[-1]) < 1e-12: return "top"
        return "floor" if fc[0] < 0 else "floor_out"
    write_mesh(P("open.h5"), np.array(coord), cells, open_tag, OPEN_KINDS)
    bc_open = P("bcond_open.yaml"); write_bcond(bc_open, OPEN_KINDS)
    res = cw.run(args_for(P("open.h5"), bc_open))
    L = res["lines"]
    check("(3) 壁 → 出口: open 1 本 (floor|floor_out)、φ = 180°、面は壁 (floor) だけ",
          len(L) == 1 and L[0]["cls"] == "open" and L[0]["tags"] == ("floor", "floor_out") and abs(L[0]["phi"][0] - 180) < 1e-9
          and {r["side"] for r in res["rows"]} == {"floor"}, str([(x["cls"], x["tags"], x["phi"]) for x in L]))
    check("(3) 壁 → 出口: kink 0、len 2.5、扇 2 セル", all(abs(r["kink_deg"]) < 1e-9 and r["fan_cells"] == 2 and r["flag"] == "" for r in res["rows"])
          and all(abs(r["len_ratio"] - 2.5) < 1e-9 for r in res["rows"] if r["layer"] == 0))

    # ---------------------------------------------------------------------------
    # (4) 入力の拒否・CLI
    # ---------------------------------------------------------------------------
    h5 = P("plate_3_40.h5")

    def err(**kw):
        try:
            cw.run(args_for(h5, bc_plate, **kw))
        except cw.WallEndError as e:
            return str(e)
        return None
    m = err(pair=["plate_up,plate_lw"])
    check("(4) --pair のタグ誤記 → エラー", m is not None and "plate_lw" in m, m or "")
    m = err(pair=["plate_up,top"])
    check("(4) 共有しない --pair → エラー", m is not None and "共有" in m, m or "")
    m = err(walls="plate_up,plate_low")
    check("(4) --walls のタグ誤記 → エラー", m is not None and "plate_low" in m, m or "")
    m = err(wall=["plate_up=99"])
    check("(4) --wall の physID が h5 に無い → エラー", m is not None and "99" in m, m or "")
    # --wall で明示 (bcondConfig なし): 名前は引数のもの
    res = cw.run(args_for(h5, None, wall=["up=5", "lo=6"]))
    check("(4) --wall NAME=ID で明示", len(res["lines"]) == 1 and res["lines"][0]["tags"] == ("lo", "up"), str([x["tags"] for x in res["lines"]]))
    # 片面だけ壁にすると knife は 1 面だけ解析
    res = cw.run(args_for(h5, bc_plate, walls="plate_lo"))
    check("(4) --walls plate_lo: 面は plate_lo だけ", {r["side"] for r in res["rows"]} == {"plate_lo"}, str({r["side"] for r in res["rows"]}))
    # 未対応のセル型
    with h5py.File(P("prism.h5"), "w") as f, h5py.File(h5, "r") as g:
        for k in ("MESH", "PLANES", "BCONDS"):
            g.copy(g[k], f, k)
        v = f.create_group("VIZMESH"); v.attrs["nVizCells"] = 1
        v.create_dataset("CONNE", data=np.array([9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5], dtype=np.int64))
    try:
        cw.run(args_for(P("prism.h5"), bc_plate)); m = None
    except cw.WallEndError as e:
        m = str(e)
    check("(4) hex 以外のセル (prism) → エラー", m is not None and "prism" in m, m or "")
    with h5py.File(P("cellish.h5"), "w") as f, h5py.File(h5, "r") as g:
        for k in ("MESH", "PLANES", "BCONDS"):
            g.copy(g[k], f, k)
    try:
        cw.run(args_for(P("cellish.h5"), bc_plate)); m = None
    except cw.WallEndError as e:
        m = str(e)
    check("(4) node 変換でない h5 (VIZMESH なし) → エラー", m is not None and "node" in m, m or "")
    # CLI
    r = subprocess.run([sys.executable, TOOL, h5, "--bcond-config", bc_plate, "--pair", "plate_up,plate_lo", "--top", "5",
                        "--out-prefix", P("out")], capture_output=True, text=True)
    ok_csv = os.path.exists(P("out_lines.csv")) and os.path.exists(P("out_nodes.csv"))
    n_lines = sum(1 for _ in open(P("out_nodes.csv"))) - 1 if ok_csv else -1
    ok = r.returncode == 0 and "終端線" in r.stdout and "上位 5 節点" in r.stdout and ok_csv and n_lines == 24
    check("(4) CLI: 終了 0、表と CSV (全行 = 3 節点 × 2 面 × 4 層)", ok, "" if ok else (r.stdout + r.stderr)[-300:])
    r = subprocess.run([sys.executable, TOOL, h5, "--bcond-config", bc_plate, "--pair", "plate_up,top", "--no-csv"], capture_output=True, text=True)
    check("(4) CLI: 共有しない --pair は終了 2", r.returncode == 2 and "ERROR" in r.stderr, r.stderr.strip()[-120:])
    r = subprocess.run([sys.executable, TOOL, "--help"], capture_output=True, text=True)
    ok = r.returncode == 0 and "定義:" in r.stdout and "限界:" in r.stdout and "knife" in r.stdout and "disp/h" in r.stdout
    check("(4) --help に定義と限界", ok, "" if ok else r.stdout[:100])

    # ---------------------------------------------------------------------------
    # (5) 後縁の線と露出した側端の線の角
    # ---------------------------------------------------------------------------
    theta = 40.0; tn = math.tan(math.radians(theta))
    xs = np.r_[-DXU * np.arange(NU, 0, -1), 0.0, DXD * np.arange(1, ND + 1)]
    ys = np.r_[-np.array(OFFS[::-1]), np.array(OFFS[1:])]
    zs = np.array([0.0, 0.05, 0.1, 0.2, 0.3]); KS = 2                  # 板は z ≤ 0.1 (k ≤ KS)、その外は板なし
    JP = len(OFFS) - 1; ITE = NU; NX, NY, NZ = len(xs), len(ys), len(zs)

    def base(i, j, k):
        return (i * NY + j) * NZ + k
    coord = [(xs[i], ys[j] - (xs[i] * tn if xs[i] > 0 else 0.0), zs[k]) for i in range(NX) for j in range(NY) for k in range(NZ)]
    dup = {}
    for i in range(ITE):                       # 上面の複製は板の内部 (i < ITE かつ k < KS) だけ。後縁と側端は共有
        for k in range(KS):
            dup[(i, k)] = len(coord); coord.append((xs[i], 0.0, zs[k]))

    def nd(i, j, k, upper):
        return dup[(i, k)] if (upper and j == JP and (i, k) in dup) else base(i, j, k)
    cells = []
    for i in range(NX - 1):
        for j in range(NY - 1):
            for k in range(NZ - 1):
                up = j >= JP
                cells.append((nd(i, j, k, up), nd(i + 1, j, k, up), nd(i + 1, j + 1, k, up), nd(i, j + 1, k, up),
                              nd(i, j, k + 1, up), nd(i + 1, j, k + 1, up), nd(i + 1, j + 1, k + 1, up), nd(i, j + 1, k + 1, up)))

    def corner_tag(fc, cc):
        if abs(fc[2]) < 1e-12: return "zmin"
        if abs(fc[2] - zs[-1]) < 1e-12: return "zmax"
        if abs(fc[0] - xs[0]) < 1e-12: return "inlet"
        if abs(fc[0] - xs[-1]) < 1e-12: return "outlet"
        if fc[0] < 0 and abs(fc[1]) < 1e-12 and fc[2] < zs[KS]:
            return "plate_up" if cc[1] > 0 else "plate_lo"
        return "bottom" if fc[1] < cc[1] else "top"
    write_mesh(P("corner.h5"), np.array(coord), cells, corner_tag, PLATE_KINDS)
    res = cw.run(args_for(P("corner.h5"), bc_plate, pair=["plate_up,plate_lo"]))
    L = res["lines"]
    te_nodes = sorted(base(ITE, JP, k) for k in range(KS + 1))
    side_nodes = sorted(base(i, JP, KS) for i in range(ITE + 1))
    got = sorted((sorted(x["nodes"]), x["cls"]) for x in L)
    check("(5) 角: knife 2 本 (後縁の線・側端の線)、角の節点は両方に入る",
          got == sorted([(te_nodes, "knife"), (side_nodes, "knife")]), str(got))
    corner = base(ITE, JP, KS)
    lid_te = [x["id"] for x in L if sorted(x["nodes"]) == te_nodes]
    lid_sd = [x["id"] for x in L if sorted(x["nodes"]) == side_nodes]
    rows = res["rows"]
    check("(5) 角: 全行 flag なし・扇 4 セル (角の節点も)", all(r["flag"] == "" and r["fan_cells"] == 4 for r in rows),
          str([(r["line"], r["side"], r["layer"], r["node"], r["flag"]) for r in rows if r["flag"]][:4]))
    bad = []
    for r in rows:
        if r["line"] in lid_te:
            ek, el, ed = plate_expect(theta, r["side"] == "plate_up", r["layer"])
        else:
            ek, el, ed = 0.0, (zs[KS + 1] - zs[KS]) / (zs[KS] - zs[KS - 1]), 0.0
        if abs(r["kink_deg"] - ek) > 1e-9 or abs(r["len_ratio"] - el) > 1e-9 or abs(r["disp_ratio"] - ed) > 1e-6 * max(1, abs(ed)):
            bad.append((r["line"], r["side"], r["layer"], r["node"], r["kink_deg"], r["len_ratio"], r["disp_ratio"]))
    check("(5) 角: 後縁の線は kink 40°・len・disp/h が期待値、側端の線は kink 0・len 2・disp 0 (角の節点を含む)", not bad, str(bad[:3]))
    rc = [r for r in rows if r["node"] == corner and r["layer"] == 0]
    check("(5) 角の節点: 後縁の線の下流 = x の次の点、側端の線の下流 = z の次の点",
          {(r["line"] in lid_te, r["down_node"]) for r in rc} == {(True, base(ITE + 1, JP, KS)), (False, base(ITE, JP, KS + 1))}, str(rc[:2]))
    check("(5) 角: --pair の共有辺 = 後縁 KS 本 + 側端 ITE 本、すべて knife",
          res["pairs"][0]["n_shared"] == KS + ITE and res["pairs"][0]["n_knife"] == KS + ITE, str(res["pairs"]))

    # ---------------------------------------------------------------------------
    # (6) 45° 未満の壁の折れ (ランプの膨張の角) は終端にしない
    # ---------------------------------------------------------------------------
    BEND_KINDS = {"inlet": (1, "inlet_uniformVelocity"), "outlet": (2, "outflow"), "top": (4, "slip"),
                  "ramp": (5, "wall"), "zmin": (7, "slip"), "zmax": (8, "slip")}
    xs = np.array([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]); ys = np.array([0.0, 0.01, 0.03, 0.1, 0.3]); zs = np.array([0.0, 0.05])
    bend = math.tan(math.radians(20.0))
    ids = np.arange(len(xs) * len(ys) * len(zs)).reshape(len(xs), len(ys), len(zs))
    coord = [(x, y - (x * bend if x > 0 else 0.0), z) for x in xs for y in ys for z in zs]   # 床が x = 0 で 20° 下へ折れる
    cells = [(ids[i, j, 0], ids[i + 1, j, 0], ids[i + 1, j + 1, 0], ids[i, j + 1, 0], ids[i, j, 1], ids[i + 1, j, 1], ids[i + 1, j + 1, 1], ids[i, j + 1, 1])
             for i in range(len(xs) - 1) for j in range(len(ys) - 1)]

    def bend_tag(fc, cc):
        if abs(fc[2]) < 1e-12: return "zmin"
        if abs(fc[2] - zs[-1]) < 1e-12: return "zmax"
        if abs(fc[0] - xs[0]) < 1e-12: return "inlet"
        if abs(fc[0] - xs[-1]) < 1e-12: return "outlet"
        return "ramp" if fc[1] < cc[1] else "top"
    write_mesh(P("bend.h5"), np.array(coord), cells, bend_tag, BEND_KINDS)
    bc_bend = P("bcond_bend.yaml"); write_bcond(bc_bend, BEND_KINDS)
    res = cw.run(args_for(P("bend.h5"), bc_bend))
    check("(6) 20° の壁の折れ (φ 200°) は既定 (225°) では終端にしない", res["lines"] == [], str([(x["cls"], x["phi"]) for x in res["lines"]]))
    res = cw.run(args_for(P("bend.h5"), bc_bend, convex_deg=190.0))
    check("(6) --convex-deg 190 なら convex 1 本 (φ 200°)・kink 20°",
          len(res["lines"]) == 1 and abs(res["lines"][0]["phi"][0] - 200.0) < 1e-9
          and all(abs(r["kink_deg"] - 20.0) < 1e-9 for r in res["rows"] if r["layer"] == 0), str([(x["cls"], x["phi"]) for x in res["lines"]]))

print("RESULT:", "PASS" if fail == 0 else f"FAIL ({fail})")
sys.exit(1 if fail else 0)
