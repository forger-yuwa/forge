"""Fluent 書き出し器・CGNS・境界条件表の試験 (plans/accepted/tooling-nozzle-wt-deliverables.md §5.1 #5、§6 V4)。

(a) 合成メッシュ (`generate_axisym_mesh` ni 11・nj 5、節点とセルを並べ替えた写し、六角形を 3 枚の四角形に割った非構造の例) で
    書いて `check_fluent_msh_2d` が全部合格。
(b) 実物との照合: TGrid 2D の elbow.msh (環境変数 `FLUENT_2D_SAMPLE`、無ければ飛ばす。GPL の配布物なのでリポジトリに入れない) で
    `measure_orientation_2d` が全 1,454 面で c0 を左 (内部 1,300・境界 154) と測る。
(c) 負例: n0, n1 を入れ替えたファイルは右 100 % で検査が不合格。時計回りの四角形・境界に属さない外周辺・内部の辺を境界に
    入れた・型の無い境界名は書き出し器が例外。
(d) M6 fixture (環境変数 `FORGE_DELIV_FIXTURE` = fixture のディレクトリ) で書き出しと検査が通る。所要時間を出す。
(e) ParaView (`pvpython` があれば) の `vtkFLUENTReader` で読み直し: fluid の点数・セル数・ゾーン名と型、全セルの符号付き面積が正、
    点が元の座標の float32 丸めと一致 (この読み込み器は点を float32 で持つので 1e-15 の照合はできない)。
(f) CGNS: gmsh で書いて読み直した節点数・要素数・物理名 (と物理グループの要素数) が msh と一致。
(g) 境界条件表 (V4(d)): 動作圧 0 Pa の行があり、bcond の全圧・背圧などの圧力と表の絶対値・入力欄の値 (絶対値 − 0) が一致。
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.export.fluent_msh import (  # noqa: E402
    OPERATING_PRESSURE_PA, check_fluent_msh_2d, measure_orientation_2d, write_cgns_via_gmsh, write_fluent_bc_table,
    write_fluent_msh_2d)
from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
SCRATCH = "/tmp/claude-1000/-home-sano-work-forge/bd18adac-85f4-4627-811b-8f3186867995/scratchpad"
ELBOW = Path(os.environ.get("FLUENT_2D_SAMPLE", f"{SCRATCH}/dl_elbow/elbow.msh"))
ELBOW_SHA256 = "a7a3c88110d19f5534b5a225eb99102f1416d5113331c61ee47de1a416c299e3"
FIXTURE = Path(os.environ.get("FORGE_DELIV_FIXTURE", f"{SCRATCH}/fixtures/m6_prod"))
GMSH_PY = Path(os.environ.get("FORGE_GMSH_PYTHON", "/home/sano/work/forge/.venv-mesh/bin/python"))
PVPYTHON = Path(os.environ.get("FORGE_PVPYTHON",
                               str(Path.home() / "opt/ParaView-6.1.1-MPI-Linux-Python3.12-x86_64/bin/pvpython")))
CASE45_PROBLEM = REPO / "case/45.isobutane_m6_d155/problem_d155_ns_prod.yaml"
fails = 0


def check(name, ok, info=""):
    global fails
    print(("ok   " if ok else "FAIL ") + name + (f" ({info})" if info else "")); fails += (not ok)


def skip(name, why):
    print(f"skip {name} ({why})")


def raises(fn, *a, **k):
    try:
        fn(*a, **k)
    except ValueError as e:
        return str(e)
    return None


class _Wall:
    """縮小・拡大のある簡単な壁 (r_t 単位)。"""
    x_in, x_e = -2.0, 4.0

    def r(self, x):
        x = np.asarray(x, dtype=float)
        return 1.0 + 0.1 * x ** 2 / (1.0 + 0.05 * x ** 2)


def small_mesh():
    return generate_axisym_mesh(_Wall(), Mesh2DParams(ni=11, nj=5, scale=0.05))


def write_and_check(td, tag, coords, quads, bedges, **kw):
    p = Path(td) / f"{tag}.msh"
    s = write_fluent_msh_2d(p, coords, quads, bedges, **kw)
    r = check_fluent_msh_2d(p, coords, quads, bedges, zone_types=kw.get("zone_types"))
    return p, s, r


def swap_n0_n1(src, dst):
    """面の節 (13) の各行の n0, n1 を入れ替えて書く (向きを逆にした負例)。"""
    out, in_face = [], False
    for line in Path(src).read_text().splitlines():
        if re.match(r"\(13 \([0-9a-f]+ [0-9a-f]+ [0-9a-f]+ [0-9a-f]+ 2\)\($", line):
            in_face = True
        elif in_face and line.startswith(")"):
            in_face = False
        elif in_face:
            t = line.split()
            line = " ".join([t[1], t[0]] + t[2:])
        out.append(line)
    Path(dst).write_text("\n".join(out) + "\n")


PV_SCRIPT = r'''
import json, sys, vtk
import numpy as np
from vtk.util import numpy_support as ns
r = vtk.vtkFLUENTReader(); r.SetFileName(sys.argv[1]); r.Update()
d = r.GetOutput()
names = [d.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME()) for i in range(d.GetNumberOfBlocks())]
b = d.GetBlock(0)
P32 = ns.vtk_to_numpy(b.GetPoints().GetData())
P = P32.astype(float)
ca = b.GetCells()
off = ns.vtk_to_numpy(ca.GetOffsetsArray()).astype(np.int64)
con = ns.vtk_to_numpy(ca.GetConnectivityArray()).astype(np.int64)
npc = np.diff(off)
cell = np.repeat(np.arange(len(npc)), npc)
nxt = np.arange(len(con)) + 1
nxt[off[1:] - 1] = off[:-1]
pa, pb = P[con, :2], P[con[nxt], :2]
ref = np.stack([np.bincount(cell, weights=pa[:, k]) / npc for k in range(2)], 1)
pa, pb = pa - ref[cell], pb - ref[cell]
area = 0.5 * np.bincount(cell, weights=pa[:, 0] * pb[:, 1] - pb[:, 0] * pa[:, 1])
np.save(sys.argv[2], P[:, :2])
print("FORGE_JSON " + json.dumps({"names": names, "points": int(b.GetNumberOfPoints()), "cells": int(b.GetNumberOfCells()),
      "nodes_per_cell": sorted(set(npc.tolist())), "n_nonpos": int(np.sum(~(area > 0))), "area_min": float(area.min()),
      "area_sum": float(area.sum()), "dtype": str(P32.dtype)}))
'''


def paraview_check(td, tag, msh, coords, quads, bedges, zones):
    cas = Path(td) / f"{tag}.cas"                      # vtkFLUENTReader は拡張子 .cas で読む
    shutil.copyfile(msh, cas)
    scr = Path(td) / "pv_check.py"
    scr.write_text(PV_SCRIPT)
    npy = Path(td) / f"{tag}_pv.npy"
    t0 = time.perf_counter()
    p = subprocess.run([str(PVPYTHON), str(scr), str(cas), str(npy)], capture_output=True, text=True, timeout=600)
    dt = time.perf_counter() - t0
    line = [s for s in p.stdout.splitlines() if s.startswith("FORGE_JSON ")]
    if p.returncode != 0 or not line:
        check(f"(e) {tag}: ParaView の読み直し", False, (p.stderr or p.stdout)[-500:])
        return
    r = json.loads(line[-1][len("FORGE_JSON "):])
    want = sorted(["fluid:fluid"] + [f"{n}:{v['type']}" for n, v in zones.items()])
    P = np.load(npy)
    xy = np.asarray(coords, float)[:, :2]
    same32 = P.shape == xy.shape and bool(np.all(P == xy.astype(np.float32).astype(float)))
    rel = float(np.max(np.abs(P - xy) / np.maximum(np.abs(xy), 1e-300))) if P.shape == xy.shape else float("nan")
    check(f"(e) {tag}: ParaView の fluid の点数・セル数・四角形・ゾーン名と型",
          r["names"][0] == "fluid:fluid" and r["points"] == len(xy) and r["cells"] == len(quads)
          and r["nodes_per_cell"] == [4] and sorted(r["names"]) == want,
          f"点 {r['points']}・セル {r['cells']}・ゾーン {r['names']}・{dt:.1f} s")
    check(f"(e) {tag}: ParaView で全セルの符号付き面積が正", r["n_nonpos"] == 0,
          f"正でない {r['n_nonpos']}・最小 {r['area_min']:.3e}・総和 {r['area_sum']:.9g}")
    check(f"(e) {tag}: ParaView の点 = 元の座標の float32 丸め (点は {r['dtype']})", same32, f"元との相対差の最大 {rel:.2e}")


def parse_md_table(path, header0):
    """`| header0 | ...` で始まる Markdown の表を [{列名: 値}] で返す。"""
    rows, hdr = [], None
    for line in Path(path).read_text().splitlines():
        if not line.startswith("|"):
            if hdr is not None and rows:
                break
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if hdr is None:
            if cells[0] == header0:
                hdr = cells
            continue
        if set(cells[0]) <= {"-", " "}:
            continue
        rows.append(dict(zip(hdr, cells)))
    return rows


def bc_table_checks(td, tag, problem, bcond, prep, expect_guess, **kw):
    p = Path(td) / f"bc_{tag}.md"
    write_fluent_bc_table(p, problem, bcond, prep, **kw)
    txt = p.read_text()
    setting = parse_md_table(p, "項目")
    op = [r for r in setting if r.get("項目") == "Operating Pressure"]
    check(f"(g) {tag}: Operating Pressure = 0 Pa の行", len(op) == 1 and re.match(r"^0(\.0+)? Pa", op[0]["値"]) is not None,
          op[0]["値"] if op else "行が無い")
    rows = parse_md_table(p, "ゾーン名")
    bc = bcond if isinstance(bcond, dict) else __import__("yaml").safe_load(Path(bcond).read_text())
    n_p, bad = 0, []
    for zone, ent in bc.items():
        for key, v in (ent.get("floats") or {}).items():
            if not str(key).startswith("P"):
                continue
            n_p += 1
            hit = [r for r in rows if re.search(rf"bcondConfig\.yaml` {zone}\.floats\.{key}$", r["出典"])]
            if len(hit) != 1:
                bad.append(f"{zone}.{key}: 行 {len(hit)}")
                continue
            r = hit[0]
            if not (r["単位"] == "Pa" and float(r["絶対値"]) == float(v)
                    and float(r["Fluent の入力欄の値"]) == float(v) - OPERATING_PRESSURE_PA):
                bad.append(f"{zone}.{key}: 表 {r['絶対値']}/{r['Fluent の入力欄の値']} ≠ bcond {v}")
    check(f"(g) {tag}: bcond の全ての圧力 ({n_p} 個) と表の絶対値・入力欄の値 (絶対値 − 動作圧) が一致", n_p > 0 and not bad, bad)
    guess = [r for r in rows if "初期ゲージ静圧の目安" in r["項目 (Fluent の入力欄)"]]
    if expect_guess:
        pt = float(bc["inlet"]["floats"]["Pt"])
        ok = len(guess) == 1 and guess[0]["絶対値"] not in ("目安なし", "不明") and 0.9 * pt < float(guess[0]["絶対値"]) < pt
        check(f"(g) {tag}: 入口の初期ゲージ静圧の目安が等エントロピーの値 (0.9 Pt < p < Pt)", ok,
              guess[0]["絶対値"] if guess else "行が無い")
    else:
        check(f"(g) {tag}: 寸法・γ が無ければ初期ゲージ静圧は「目安なし」", len(guess) == 1 and guess[0]["絶対値"] == "目安なし",
              guess[0]["出典"] if guess else "行が無い")
    return p, rows, txt


with tempfile.TemporaryDirectory() as td:
    # ---------------- (a) 合成メッシュ ----------------
    c, q, b = small_mesh()
    p_small, s, r = write_and_check(td, "small", c, q, b)
    check("(a) 小メッシュ (11×5): check_fluent_msh_2d が全部合格", r["ok"], r["failures"])
    check("(a) ゾーン ID は wall 3・inlet 4・outlet 5・axis 6 (名前順でない)",
          [(n, v["id"]) for n, v in s["zones"].items()] == [("interior-fluid", 2), ("wall", 3), ("inlet", 4), ("outlet", 5), ("axis", 6)],
          s["zones"])
    ori = r.get("orientation", {})
    check("(a) 幾何から測って c0 が左 100 %・内部面の c1 が反対側", ori.get("left_fraction") == 1.0
          and ori["interior"]["c1_opposite"] == ori["interior"]["left"], ori)
    check("(a) 座標の相対差 0 (17 有効桁)・面積の総和・境界の閉路",
          r.get("coord_max_rel") == 0.0 and r.get("area_total_rel", 1) <= 1e-12 and r.get("boundary_single_loop"),
          {k: r.get(k) for k in ("coord_max_rel", "area_total_rel", "boundary_loop_len", "boundary_loop_area_rel")})
    txt = p_small.read_text()
    check("(a) 節の書式 (2 2)・10・12 (型 1・要素 3)・13 (0 番の宣言)・45", all(t in txt for t in (
        "(2 2)", f"(10 (0 1 {len(c):x} 0 2))", f"(12 (1 1 {len(q):x} 1 3))", "(45 (1 fluid fluid)())",
        "(45 (6 axis axis)())", "(13 (6 ")) and txt.isascii())
    # 並べ替えた写し: 節点番号・セルの順・各四角形の始点をずらす (一般の四角形メッシュとして扱えるか)
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(c))
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(c))
    c2 = c[perm]
    q2 = inv[q][rng.permutation(len(q))]
    q2 = np.array([np.roll(row, int(k)) for row, k in zip(q2, rng.integers(0, 4, len(q2)))])
    b2 = {n: [(int(inv[e1]), int(inv[e0])) for e0, e1 in v] for n, v in b.items()}   # 境界辺の向きも逆にして渡す
    _, _, r2 = write_and_check(td, "perm", c2, q2, b2)
    check("(a) 節点・セル・始点を並べ替え境界辺を逆向きに渡しても合格", r2["ok"], r2["failures"])
    # 非構造: 正六角形を中心節点で 3 枚の四角形に割る (価数 3 の節点)
    ang = np.deg2rad(np.arange(6) * 60.0)
    c3 = np.vstack([[0.0, 0.0], np.c_[1.0 + np.cos(ang), 1.5 + np.sin(ang)]])
    c3[0] = [1.0, 1.5]
    q3 = np.array([[0, 1, 2, 3], [0, 3, 4, 5], [0, 5, 6, 1]])
    b3 = {"wall": [(1, 2), (2, 3), (3, 4)], "inlet": [(4, 5), (5, 6), (6, 1)]}
    _, s3, r3 = write_and_check(td, "hex", c3, q3, b3)
    check("(a) 非構造 (六角形を 3 枚の四角形に割った) も合格", r3["ok"] and s3["zones"]["interior-fluid"]["faces"] == 3, r3["failures"])

    # ---------------- (b) 実物の 2D Fluent メッシュ ----------------
    if ELBOW.is_file():
        sha = hashlib.sha256(ELBOW.read_bytes()).hexdigest()
        check("(b) elbow.msh の sha256 が測定したものと同じ", sha == ELBOW_SHA256, sha)
        m = measure_orientation_2d(ELBOW)
        check("(b) 実物 (TGrid 2D elbow.msh) で c0 が左: 全 1,454 面 (内部 1,300・境界 154)",
              m["faces"] == 1454 and m["interior"]["left"] == 1300 and m["interior"]["right"] == 0
              and m["boundary"]["left"] == 154 and m["boundary"]["right"] == 0 and m["interior"]["c1_opposite"] == 1300,
              f"節点 {m['nodes']}・セル {m['cells']}・内部 {m['interior']}・境界 {m['boundary']}")
    else:
        skip("(b) 実物の 2D Fluent メッシュとの照合", f"{ELBOW} が無い (FLUENT_2D_SAMPLE で指定)")

    # ---------------- (c) 負例 ----------------
    p_rev = Path(td) / "small_reversed.msh"
    swap_n0_n1(p_small, p_rev)
    mr = measure_orientation_2d(p_rev)
    rr = check_fluent_msh_2d(p_rev, c, q, b)
    check("(c) n0, n1 を入れ替えたファイルは右 100 %", mr["interior"]["right"] + mr["boundary"]["right"] == mr["faces"]
          and mr["boundary"]["right"] > 0, mr)
    check("(c) n0, n1 を入れ替えたファイルは check_fluent_msh_2d が不合格", not rr["ok"], rr["failures"][:2])
    e = raises(write_fluent_msh_2d, Path(td) / "cw.msh", c, q[:, ::-1], b)
    check("(c) 時計回りの四角形は書き出し器が例外", e is not None and "時計回り" in e, e)
    bm = dict(b, wall=b["wall"][:-1])
    e = raises(write_fluent_msh_2d, Path(td) / "miss.msh", c, q, bm)
    check("(c) どの境界にも属さない外周辺があれば例外", e is not None and "どの境界にも属さない" in e, e)
    bi = dict(b, wall=list(b["wall"]) + [(int(q[0, 1]), int(q[0, 2]))])
    e = raises(write_fluent_msh_2d, Path(td) / "inner.msh", c, q, bi)
    check("(c) 内部の辺を境界に入れたら例外", e is not None and "内部の辺" in e, e)
    bx = dict(b, wall=list(b["wall"]) + [(0, len(c) - 1)])
    e = raises(write_fluent_msh_2d, Path(td) / "nonedge.msh", c, q, bx)
    check("(c) 四角形の辺でない境界辺は例外", e is not None and "一致しない" in e, e)
    e = raises(write_fluent_msh_2d, Path(td) / "empty.msh", c, q, dict(b, extra=[]), zone_types={"extra": 3})
    check("(c) 辺の無い境界は例外 (空の面ゾーン)", e is not None and "辺が無い" in e, e)
    e = raises(write_fluent_msh_2d, Path(td) / "unk.msh", c, q, {("farfield" if k == "outlet" else k): v for k, v in b.items()})
    check("(c) 型の無い境界名は例外 (zone_types で与える)", e is not None and "型が無い" in e, e)
    _, s4, r4 = write_and_check(td, "ztype", c, q, {("farfield" if k == "outlet" else k): v for k, v in b.items()},
                                zone_types={"farfield": 9})
    check("(c') zone_types で型を与えれば書けて合格 (farfield = pressure-far-field、ID は既定 4 名の後)",
          r4["ok"] and s4["zones"]["farfield"] == {"id": 6, "type": "pressure-far-field", "faces": 4}, s4["zones"])

    # ---------------- (d) M6 fixture ----------------
    fx_dir = FIXTURE if FIXTURE.is_dir() else FIXTURE.parent
    fx_msh = FIXTURE if FIXTURE.is_file() else FIXTURE / "nozzle.msh"
    m6 = None
    if fx_msh.is_file():
        t0 = time.perf_counter()
        try:
            from forge_design.export.mesh_structure import read_msh41_mesh2d
            cm, qm, bm6, _ = read_msh41_mesh2d(fx_msh)
            src = "read_msh41_mesh2d"
        except Exception as ex:                       # 読み込み器が無い/作りかけ: 壁ファイルと格子パラメータから作り直す
            from forge_design.geometry.wall_axismach import load_wall_file
            prep = json.loads((fx_dir / "prepare_info.json").read_text())
            wall = load_wall_file(fx_dir / "wall_repr.json")["physical"]
            cm, qm, bm6 = generate_axisym_mesh(wall, Mesh2DParams(**prep["mesh"]["params"]))
            src = f"generate_axisym_mesh (read_msh41_mesh2d が使えない: {type(ex).__name__})"
        t1 = time.perf_counter()
        p_m6 = Path(td) / "m6.msh"
        s6 = write_fluent_msh_2d(p_m6, cm, qm, bm6)
        t2 = time.perf_counter()
        r6 = check_fluent_msh_2d(p_m6, cm, qm, bm6)
        t3 = time.perf_counter()
        m6 = (p_m6, cm, qm, bm6, s6)
        print(f"     M6: 入力 {src} {t1 - t0:.2f} s・書き出し {t2 - t1:.2f} s・検査 {t3 - t2:.2f} s・"
              f"{p_m6.stat().st_size / 1e6:.1f} MB")
        print(f"     M6: 節点 {r6['nodes']}・セル {r6['cells']}・面 {s6['faces']} {r6['zone_faces']}")
        print(f"     M6: 座標の相対差 {r6.get('coord_max_rel')}・面積の総和 {r6.get('area_total')} m² (元との相対差 "
              f"{r6.get('area_total_rel')})・最小セル面積 {r6.get('cell_area_min')}・境界の閉路 {r6.get('boundary_loop_len')} 面 "
              f"(囲む面積の相対差 {r6.get('boundary_loop_area_rel')})・向き {r6.get('orientation', {}).get('interior')} "
              f"{r6.get('orientation', {}).get('boundary')}")
        check("(d) M6 (194,000 節点・191,904 セル) の書き出しと検査が全部合格",
              r6["ok"] and r6["nodes"] == 194000 and r6["cells"] == 191904, r6["failures"])
        check("(d) M6 の書き出し + 検査が 1 分以内 (目安)", t3 - t1 < 60.0, f"{t3 - t1:.1f} s")
    else:
        skip("(d) M6 fixture", f"{fx_msh} が無い (FORGE_DELIV_FIXTURE で fixture のディレクトリを指定)")

    # ---------------- (e) ParaView ----------------
    if PVPYTHON.is_file():
        paraview_check(td, "small", p_small, c, q, b, s["zones"])
        if m6 is not None:
            paraview_check(td, "m6", m6[0], m6[1], m6[2], m6[3], m6[4]["zones"])
    else:
        skip("(e) ParaView の読み直し", f"{PVPYTHON} が無い (FORGE_PVPYTHON で指定)")

    # ---------------- (f) CGNS ----------------
    if GMSH_PY.is_file():
        msh_small = Path(td) / "small_gmsh.msh"
        write_msh41_2d(msh_small, c, q, b)
        cases = [("small", msh_small, len(c), len(q), {k: len(v) for k, v in b.items()})]
        if fx_msh.is_file() and m6 is not None:
            cases.append(("m6", fx_msh, len(m6[1]), len(m6[2]), {k: len(v) for k, v in m6[3].items()}))
        for tag, msh, nn, nq, nb in cases:
            cg = write_cgns_via_gmsh(msh, Path(td) / f"{tag}.cgns", GMSH_PY)
            a, g = cg["msh"], cg["cgns"]
            want_phys = sorted([[1, n, k] for n, k in nb.items()] + [[2, "fluid", nq]])
            got_phys = sorted([[d, n, k] for d, _, n, k in g["physical"]])
            check(f"(f) {tag}: CGNS を読み直して節点数・要素数・物理名 (と要素数) が msh と一致",
                  cg["match"] and g["nodes"] == nn and g["elements"].get("Quadrilateral 4") == nq
                  and g["elements"].get("Line 2") == sum(nb.values()) and got_phys == want_phys,
                  f"節点 {g['nodes']}・要素 {g['elements']}・物理 {got_phys}・{cg['seconds']:.1f} s"
                  + ("" if cg["tags_preserved"] else "・物理タグの番号は振り直される"))
    else:
        skip("(f) CGNS", f"{GMSH_PY} が無い (FORGE_GMSH_PYTHON で指定)")

    # ---------------- (g) 境界条件表 ----------------
    bc_yaml, prep_json = fx_dir / "bcondConfig.yaml", fx_dir / "prepare_info.json"
    if bc_yaml.is_file() and prep_json.is_file():
        problem = CASE45_PROBLEM if CASE45_PROBLEM.is_file() else None
        _, rows, txt = bc_table_checks(td, "m6", problem, bc_yaml, prep_json, expect_guess=problem is not None)
        if problem is not None:
            check("(g) m6: ガス組成 (質量分率)・気体モデル・輸送物性・節点中心とセル中心の注記",
                  all(t in txt for t in ("組成 CO2 (質量分率) | 0.1677", "| 気体モデル (forge) | semiperfect",
                                         "H2O: custom:h2o_iapws_cea_v1", "節点中心", "セル中心", "約半分",
                                         "Fluent 側の物性の選び方")))
        wall = [r for r in rows if r["ゾーン名"] == "wall"]
        check("(g) m6: 壁 (kind wall) は断熱 (熱流束 0)", len(wall) == 1 and "断熱" in wall[0]["項目 (Fluent の入力欄)"]
              and wall[0]["絶対値"] == "0.0", wall)
        bc_table_checks(td, "m6_noproblem", None, bc_yaml, None, expect_guess=False)
        # physID から msh の物理名を引く: bcond のキー名を変えてもゾーン名は msh の名前になる
        bc_ren = __import__("yaml").safe_load(bc_yaml.read_text())
        bc_ren = {("nozzle_wall" if k == "wall" else k): v for k, v in bc_ren.items()}
        phys = {"inlet": 1, "outlet": 2, "wall": 3, "axis": 4, "fluid": 5}
        _, rows_r, txt_r = bc_table_checks(td, "m6_renamed", problem, bc_ren, prep_json, expect_guess=problem is not None,
                                           msh_physical=phys)
        check("(g) m6_renamed: physID 3 の bcond `nozzle_wall` は msh のゾーン名 wall で書き、注記する",
              [r["ゾーン名"] for r in rows_r if "断熱" in r["項目 (Fluent の入力欄)"]] == ["wall"]
              and "`nozzle_wall` (physID 3) は msh のゾーン `wall`" in txt_r)
    else:
        skip("(g) fixture の境界条件表", f"{bc_yaml} / {prep_json} が無い")
    # 等温壁・未知の kind (辞書で渡す)
    bc_iso = {"inlet": {"physID": 1, "kind": "inlet_Pressure", "floats": {"Pt": 1.0e6, "Tt": 400.0}},
              "outlet": {"physID": 2, "kind": "outlet_statPress", "floats": {"Ps": 1.0e4}},
              "wall": {"physID": 3, "kind": "wall_isothermal", "floats": {"Ux": 0.0, "Uy": 0.0, "Uz": 0.0, "Ts": 350.0}},
              "axis": {"physID": 4, "kind": "axis", "floats": None},
              "extra": {"physID": 5, "kind": "inlet_supersonic", "floats": {}}}
    prob_iso = {"gas": {"model": "cpg", "gamma": 1.4, "cp": 1005.0}, "geometry": {"r_inlet": 3.0},
                "spec": {"wall_thermal": {"mode": "isothermal", "Tw": 350.0}}}
    _, rows, txt = bc_table_checks(td, "isothermal", prob_iso, bc_iso, {"throat_physical": {"r": 1.0}}, expect_guess=True)
    wall = [r for r in rows if r["ゾーン名"] == "wall"]
    check("(g) 等温壁 (wall_isothermal) は Ts を写す", len(wall) == 1 and "等温壁" in wall[0]["項目 (Fluent の入力欄)"]
          and wall[0]["絶対値"] == "350.0" and wall[0]["単位"] == "K", wall)
    ex = [r for r in rows if r["ゾーン名"] == "extra"]
    check("(g) 対応を決めていない kind は「不明」", len(ex) == 1 and ex[0]["絶対値"] == "不明", ex)
    gs = [r for r in rows if "初期ゲージ静圧" in r["項目 (Fluent の入力欄)"]]
    # A/A* = 9、γ 1.4 の亜音速解 (独立に: M ≈ 0.06483, p/Pt ≈ 0.99707)
    def _area_ratio_g14(m):                            # γ = 1.4 の A/A* (指数 (γ+1)/(2(γ−1)) = 3)
        return (1.0 / m) * ((2.0 / 2.4) * (1.0 + 0.2 * m * m)) ** 3

    mref = 0.0648
    for _ in range(100):                              # Newton 法 (微分は中心差分)
        h = 1e-7 * mref
        mref -= (_area_ratio_g14(mref) - 9.0) / ((_area_ratio_g14(mref + h) - _area_ratio_g14(mref - h)) / (2 * h))
    pref = 1.0e6 * (1 + 0.2 * mref * mref) ** -3.5
    check("(g) 初期ゲージ静圧の目安 = 等エントロピーの静圧 (A/A* 9・γ 1.4、Newton 法の参照と相対 1e-9)",
          len(gs) == 1 and abs(float(gs[0]["絶対値"]) - pref) <= 1e-9 * pref, (gs[0]["絶対値"] if gs else None, pref))
    check("(g) cpg のとき semiperfect の注記を書かない", "NASA-9" not in txt)

print(f"FAIL 件数: {fails}"); sys.exit(1 if fails else 0)
