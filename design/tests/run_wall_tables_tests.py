#!/usr/bin/env python3
"""形状表・3D STEP/STL・CFD 壁面量表・xlsx の試験 (plans/active/tooling-nozzle-wt-deliverables.md §5.1 #3・#6、§6 V2・V5)。

1. 形状表 (fixture の生産の壁 `wall_repr.json`): 列名 ASCII・ヘッダ 1 行・両端を含む・刻み (1 mm / スロート付近 0.1 mm)・
   r が `load_wall_file` の評価と 1e-9 mm・曲率半径の符号と inf・累積側面積の最終値が ∫2πr ds (20001 点) と相対 1e-6
2. 3D STEP/STL (FreeCAD): 読み直しで面 1 枚以上・妥当・面積正・ノット区間ごとの面積が ∫2πr ds と相対 1e-6・bbox の x 範囲が定義域と 1e-6 mm
3. gmsh (OCC、別経路) の読み直し (`cad_check_gmsh.py`): 面積が ∫2πr ds と相対 1e-6・bbox の x 範囲が定義域と 1e-6 mm。
   形状表の累積側面積の最終値と FreeCAD・gmsh の面積が相対 1e-6 (V2・V5)
4. CFD 表 (case/45 run_0051、断熱壁・iface なし): `missing` と列名の切り替え (`allow_diagnostic`)・ring_area の総和 = 壁の多角形の
   側面積 (メッシュの辺から別に計算) 相対 1e-12・y1+ が `check_wall_resolution.py` の CLI の要約 (平均・最大) と一致
5. CFD 表 (case/63 軸対称・等温): 小さな run で `valid`、収束した ΔT 10 K の run で表の q と −qwall・iface_q_compact が内部の 2 箇所以上で
   5 % 以内・同じ符号、ΣQ_node = Σ q A_node (1e-9)、断熱壁で `not_applicable`
6. 合成データの単位試験: 半径 0.01・0.05 m の 2 本の円筒で `iface_Qf_eff` = q r Δs [W/rad] から q が相対 1e-12 で戻る・旧定義
   (`iface_q_eff` = r q) は r 倍ずれる (負例)・`iface_q_eff` だけでは `missing` (置き換えない)・一部 NaN で `nonfinite`・全 NaN で
   `not_applicable`・異なる step の混入で例外・平面の run・y1+ の式・接線応力
7. `wall_heat_series`・`thermal_use` (合成 run で PASS / 判定材料なしで未評価、未収束の case/45 で FAIL)
8. xlsx (`.venv-pptx` が無ければ飛ばしたことを返す)
9. 読んだ他セッションの run ディレクトリに何も書いていない (一覧・大きさ・時刻が前後で同じ)

出力はテンポラリ (環境変数 `WALL_TABLES_TEST_OUT` を指定するとそこに残す)。fixture は `FORGE_DELIV_FIXTURE` で上書きできる。
usage (cwd = リポジトリ): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_wall_tables_tests.py
"""
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.export import wall_tables as WT  # noqa: E402
from forge_design.export.wall_step import FREECADCMD, write_surface  # noqa: E402
from forge_design.geometry.wall_axismach import load_wall_file  # noqa: E402

FIX = Path(os.environ.get("FORGE_DELIV_FIXTURE",
                          "/tmp/claude-1000/-home-sano-work-forge/bd18adac-85f4-4627-811b-8f3186867995/scratchpad/fixtures/m6_prod"))
MESH_PY = os.environ.get("FORGE_MESH_PY", "/home/sano/work/forge/.venv-mesh/bin/python")
GMSH_CHECK = Path(WT.__file__).with_name("cad_check_gmsh.py")
RUN45 = Path("/home/sano/work/forge/case/45.isobutane_m6_d155/run_0051_ns_final_c2")
C63 = Path("/home/sano/work/forge-cht/case/63.graetz_cht")
R63_SMOKE = C63 / "run_0004_smoke_walls_iso_r16"     # 小さな等温 run (8 step の壁ダンプ、ΔT 0、未収束)
R63_CONV = C63 / "run_0013_g1_dT10_r32"              # 収束した等温 run (ΔT 10 K、600000 step、FP64)
R63_DRY = C63 / "run_0001_dry_r16"
CWR = WT.TOOLS / "check_wall_resolution.py"
FAIL = 0
SKIPPED = []
INFO = {}


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def raises(fn, exc=Exception):
    try:
        fn()
    except exc as e:  # noqa: PERF203
        return f"{type(e).__name__}: {e}"
    return None


def snap(d: Path):
    """ディレクトリ直下の一覧・大きさ・時刻 (他セッションの run に書いていないことの確認)。"""
    return sorted((p.name, p.stat().st_size, p.stat().st_mtime_ns) for p in d.iterdir()) if d.is_dir() else None


def read_csv(p):
    with open(p) as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def col(hdr, rows, name):
    j = hdr.index(name)
    return np.array([float(r[j]) if r[j] != "" else np.nan for r in rows])


def cli_y1plus(run_dir, group, step):
    q = subprocess.run([sys.executable, str(CWR), str(run_dir), "--groups", group, "--step", str(step)],
                       capture_output=True, text=True, timeout=600)
    m = re.search(r"y1\+ 平均\s+([0-9.]+) / p99\s+([0-9.]+) / 最大\s+([0-9.]+)", q.stdout)
    return (float(m.group(1)), float(m.group(3))) if m else None, q.stdout


def lateral_ref(W, n=20001):
    """∫2πr ds [mm²] (wall_repr の曲線を n 点の Simpson で)。"""
    from scipy.integrate import simpson
    s = 1000.0 * W["scale_m"]
    x0, x1 = W["domain"]
    x = np.linspace(x0, x1, n)
    f = 2.0 * np.pi * W["physical"].r(x) * s * np.sqrt(1.0 + W["physical"].r(x, 1) ** 2) * s
    return float(simpson(f, x=x))


# ---------------------------------------------------------------------------------------------------------------------
# 合成データ (半径の違う 2 本の円筒)
# ---------------------------------------------------------------------------------------------------------------------

SYN_WALLS = {3: ("wall_a", 0.01, [0.0, 0.01, 0.025, 0.04, 0.06, 0.08, 0.1], 0.001),
             4: ("wall_b", 0.05, [0.2, 0.22, 0.25, 0.27, 0.3], 0.002)}
SYN_VISC, SYN_RO = 1.8e-5, 1.2


def syn_q(pid, x):
    return 1000.0 * (1.0 + 10.0 * x) if pid == 3 else -3000.0 * (1.0 + x)


def make_synth(td: Path, axi=True, steps=(100,), qf="valid", tw=5.0, residual=False, only_q_eff=False):
    """node 方式・2D の合成 run: 壁 2 本 (円筒 r = 0.01 / 0.05 m)、各壁節点の内側に 1 節点。壁ダンプは iCells と違う順。
    qf: valid / nan_one / nan_all。iface_Qf_eff = q r_f |S_f| (軸対称、1 rad あたり) または q |S_f| (平面)。"""
    td.mkdir(parents=True, exist_ok=True)
    coords, wall_ids, inner_ids, planes = [], {}, {}, []
    for pid, (_, r, xs, dr) in SYN_WALLS.items():
        w0 = len(coords)
        coords += [[x, r, 0.0] for x in xs]
        wall_ids[pid] = list(range(w0, w0 + len(xs)))
        i0 = len(coords)
        coords += [[x, r - dr, 0.0] for x in xs]
        inner_ids[pid] = list(range(i0, i0 + len(xs)))
    coords = np.array(coords)
    sa, sv, pc, struct, bc = [], [], [], [], {}
    for pid, (_, r, xs, _) in SYN_WALLS.items():
        ip, ds = [], []
        for i, x in enumerate(xs):
            left = 0.5 * (x - xs[i - 1]) if i > 0 else 0.0
            right = 0.5 * (xs[i + 1] - x) if i < len(xs) - 1 else 0.0
            d = left + right
            xc = (left * (x - 0.5 * left) + right * (x + 0.5 * right)) / d
            ip.append(len(sa))
            sa.append(d); sv.append([0.0, d, 0.0]); pc.append([xc, r, 0.0])
            struct += [1, wall_ids[pid][i], 1, wall_ids[pid][i]]
            ds.append(d)
        bc[pid] = (ip, ds)
    for pid in SYN_WALLS:
        for w, n in zip(wall_ids[pid], inner_ids[pid]):
            sa.append(1.0); sv.append([1.0, 0.0, 0.0]); pc.append(list(0.5 * (coords[w] + coords[n])))
            struct += [2, w, n, 2, w, n]
    with h5py.File(td / "mesh.h5", "w") as f:
        f["MESH/COORD"] = coords.ravel()
        f["CELLS/centCoords"] = coords.ravel()
        f["CELLS/volume"] = np.ones(len(coords))
        f["PLANES/surfArea"] = np.array(sa)
        f["PLANES/surfVect"] = np.array(sv).ravel()
        f["PLANES/centCoords"] = np.array(pc).ravel()
        f["PLANES/STRUCT"] = np.array(struct, dtype=np.int32)
        for pid in SYN_WALLS:
            ip, _ = bc[pid]
            order = list(range(len(ip)))[::-1]                      # iCells を逆順に (順序を仮定しない)
            f[f"BCONDS/{pid}/iCells"] = np.array([wall_ids[pid][i] for i in order], dtype=np.int32)
            f[f"BCONDS/{pid}/iPlanes"] = np.array([ip[i] for i in order], dtype=np.int32)
            f[f"BCONDS/{pid}/iBPlanes"] = np.arange(len(ip), dtype=np.int32)
            w = wall_ids[pid]
            f[f"BCONDS/{pid}/vizBfaceNodes"] = np.array([[w[i + 1], w[i]] for i in range(len(w) - 1)], dtype=np.int32).ravel()
            f[f"BCONDS/{pid}/vizBfaceSizes"] = np.full(len(w) - 1, 2, dtype=np.int32)
    (td / "solverConfig.yaml").write_text(
        "mesh: {discretization: \"node\", isAxisymmetric: %d, axisymMethod: 0, meshFileName: \"mesh.h5\", valueFileName: \"mesh.h5\"}\n"
        "physProp: {viscMethod: 0, visc: %r}\ntime: {unsteady: 0}\noutput: {level: 1, interfaceDiag: 1}\n" % (1 if axi else 0, SYN_VISC))
    (td / "bcondConfig.yaml").write_text("".join(
        "%s: {physID: %d, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ts: 300.0}}\n" % (nm, pid)
        for pid, (nm, *_) in SYN_WALLS.items()) + "axis: {physID: 6, kind: axis, outputHDFflg: 0, ints: , floats: }\n")
    truth = {}
    for pid, (nm, r, xs, _) in SYN_WALLS.items():
        ip, ds = bc[pid]
        xs_, ds_ = np.array(xs), np.array(ds)
        q = syn_q(pid, xs_)
        Qf = q * r * ds_ if axi else q * ds_
        q_eff_old = Qf / ds_                                    # main のソルバの iface_q_eff (r を掛けない面積で割る) = r q
        if qf == "nan_one":
            Qf = Qf.copy(); Qf[2] = np.nan
        elif qf == "nan_all":
            Qf = np.full_like(Qf, np.nan)
        perm = [2, 0, 4, 1, 3] + list(range(5, len(xs)))            # 壁ダンプの順も別に
        for st in steps:
            with h5py.File(td / f"res_{nm}_{pid}_{st}.h5", "w") as f:
                f.attrs["step"] = np.int32(st)
                f["MESH/COORD"] = coords[[wall_ids[pid][i] for i in perm]].ravel()
                vals = {"ro": np.full(len(xs), SYN_RO), "Ps": np.full(len(xs), 1.0e5), "Ts": np.full(len(xs), 300.0),
                        "twall_x": np.full(len(xs), tw), "twall_y": np.zeros(len(xs)), "twall_z": np.zeros(len(xs)),
                        "qwall": -1.01 * q, "iface_q_eff": q_eff_old}
                if not only_q_eff:
                    vals["iface_Qf_eff"] = Qf
                for k, v in vals.items():
                    f["VALUE/" + k] = np.asarray(v, dtype=float)[perm]
        truth[pid] = {"x": xs_, "r": r, "ds": ds_, "q": q, "Qf": Qf, "q_eff_old": q_eff_old}
    if residual:
        n = 500
        with open(td / "residual_history.csv", "w") as f:
            f.write("step,inner,phase,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n")
            for i in range(n):
                v = 10.0 ** (-6.0 * i / (n - 1))
                f.write(f"{i},-1,outer_end,{v},{v},{v},0,{v}\n")
    return truth


def main():
    keep = os.environ.get("WALL_TABLES_TEST_OUT")
    out = Path(keep) if keep else Path(tempfile.mkdtemp(prefix="wall_tables_tests_"))
    out.mkdir(parents=True, exist_ok=True)
    watched = [RUN45, R63_SMOKE, R63_CONV, R63_DRY]
    before = {str(d): snap(d) for d in watched}
    try:
        run(out)
    finally:
        after = {str(d): snap(d) for d in watched}
        check("読んだ他セッションの run ディレクトリに何も書いていない (一覧・大きさ・時刻が同じ): "
              + ", ".join(Path(d).name for d in before), before == after)
        if not keep:
            shutil.rmtree(out, ignore_errors=True)
    print("V2・V5 の数値:", json.dumps(INFO, ensure_ascii=False, default=float))
    print(f"FAIL 件数: {FAIL}")
    if SKIPPED:
        print("飛ばした項目 (合格扱いにしない):", SKIPPED)
    sys.exit(1 if FAIL else (2 if SKIPPED else 0))


def run(out: Path):
    # ---------------- 1. 形状表 ----------------
    W = load_wall_file(FIX)
    s_mm = 1000.0 * W["scale_m"]
    x0, x1 = W["domain"]
    ref_area = lateral_ref(W)
    INFO["lateral_area_ref_mm2 (∫2πr ds、20001 点)"] = ref_area
    shp = WT.shape_table(FIX, out / "excel" / "wall_shape_mm.csv")
    raw = (out / "excel" / "wall_shape_mm.csv").read_text().splitlines()
    hdr, rows = read_csv(out / "excel" / "wall_shape_mm.csv")
    check(f"形状表: 列名 {hdr} が仕様どおり・ASCII", hdr == list(WT.SHAPE_COLUMNS) and all(c.isascii() for c in hdr))
    check("形状表: ヘッダは列名だけの 1 行 (コメント行なし、2 行目以降はすべて数値か inf)",
          raw[0] == ",".join(WT.SHAPE_COLUMNS) and all(not ln.startswith("#") for ln in raw)
          and all(np.isfinite(float(v)) or v == "inf" for r in rows for v in r))
    x = col(hdr, rows, "x_mm")
    r = col(hdr, rows, "r_mm")
    check(f"形状表: 両端を含む (x {float(x[0])!r} .. {float(x[-1])!r} mm = 定義域 × {s_mm:.6f})",
          x[0] == x0 * s_mm and x[-1] == x1 * s_mm and np.all(np.diff(x) > 0))
    hw = max(50.0, 2.0 * s_mm)
    dx = np.diff(x)
    inz = (x[:-1] >= -hw) & (x[1:] <= hw)
    check(f"形状表: 刻み 1 mm 以下・スロート付近 |x| ≤ {hw:.3f} mm は 0.1 mm 以下 (最大 {dx.max():.6f} / {dx[inz].max():.6f})",
          dx.max() <= 1.0 + 1e-9 and dx[inz].max() <= 0.1 + 1e-9 and np.any(np.isclose(x, 0.1)))
    rr = W["physical"].r(np.clip(x / s_mm, x0, x1)) * s_mm
    check(f"形状表: r が load_wall_file の評価と一致 (最大差 {np.abs(r - rr).max():.2e} mm ≤ 1e-9)", np.abs(r - rr).max() <= 1e-9)
    rc = col(hdr, rows, "curvature_radius_mm")
    it = int(np.argmin(np.abs(x - W["throat"]["x"] * s_mm)))
    check(f"形状表: 曲率半径の符号 (物理スロート付近 {rc[it]:.3f} mm > 0 = 流路側から凸)・直管は inf ({int(np.isinf(rc).sum())} 点)",
          rc[it] > 0 and np.isinf(rc[0]) and np.isinf(rc).sum() > 0)
    per = col(hdr, rows, "perimeter_mm")
    check("形状表: 周長 = 2πr", np.allclose(per, 2 * np.pi * r, rtol=1e-15, atol=0))
    A_tab = col(hdr, rows, "lateral_area_cum_mm2")[-1]
    s_tab = col(hdr, rows, "s_mm")[-1]
    INFO["shape_lateral_area_final_mm2"] = A_tab
    INFO["shape_area_err_est_rel"] = shp["integration"]["area_err_est_rel"]
    check(f"形状表: 累積側面積の最終値 {A_tab:.6f} mm² と ∫2πr ds (20001 点) {ref_area:.6f} が相対 {abs(A_tab / ref_area - 1):.1e} ≤ 1e-6",
          abs(A_tab / ref_area - 1) <= 1e-6)
    check(f"形状表: 弧長の最終値 {s_tab:.6f} mm が戻り値と一致・積分の精度の目安 (面積 相対 {shp['integration']['area_err_est_rel']:.1e}) を返す",
          s_tab == shp["integration"]["s_total_mm"] and shp["integration"]["area_err_est_rel"] < 1e-8)
    check("形状表: 負例 — n_sub が奇数なら例外", raises(lambda: WT.shape_table(FIX, out / "bad.csv", n_sub=15), ValueError) is not None)

    # ---------------- 2. 3D STEP/STL (FreeCAD) ----------------
    step2d = FIX / "step" / "wall_physical.step"
    step3d, stl = out / "cad" / "wall_surface.step", out / "cad" / "wall_surface.stl"
    step3d.parent.mkdir(parents=True, exist_ok=True)
    fc = None
    if Path(FREECADCMD).is_file() and step2d.is_file():
        fc = write_surface(step2d, step3d, stl)
        a_fc = fc["readback_area_knot_split_mm2"]
        INFO["freecad_readback_area_knot_split_mm2"] = a_fc
        INFO["freecad_shape_area_mm2 (OCC GProp, 参考)"] = fc["readback_area_mm2"]
        check(f"3D STEP (FreeCAD): 読み直しで面 {fc['readback_faces']} 枚・妥当 {fc['readback_valid']}・面積 > 0",
              fc["readback_faces"] >= 1 and fc["readback_valid"] and fc["readback_area_mm2"] > 0 and fc["is_valid"])
        check(f"3D STEP (FreeCAD): ノット区間ごとの面積 {a_fc:.6f} mm² が ∫2πr ds と相対 {abs(a_fc / ref_area - 1):.1e} ≤ 1e-6",
              a_fc is not None and abs(a_fc / ref_area - 1) <= 1e-6)
        print(f"     参考: FreeCAD の Shape.Area (1 枚の回転面の GProp) は {fc['readback_area_mm2']:.3f} mm² (相対 {fc['readback_area_mm2'] / ref_area - 1:+.2e})")
        bb = fc["readback_bbox_mm"]
        ex = max(abs(bb[0] - x0 * s_mm), abs(bb[3] - x1 * s_mm))
        INFO["freecad_bbox_x_err_mm"] = ex
        check(f"3D STEP (FreeCAD): bbox の x 範囲が定義域と {ex:.1e} mm ≤ 1e-6", ex <= 1e-6)
        check(f"STL: 書けて三角形 {fc['stl_triangles']} 枚 (弦誤差の設定 {fc['stl_linear_deflection_mm']} mm / {fc['stl_angular_deflection_rad']} rad)",
              stl.is_file() and stl.stat().st_size > 0 and fc["stl_triangles"] > 0)
        check(f"3D STEP: 相対 1e-6 で形状表の累積側面積の最終値と FreeCAD の面積が一致 ({abs(A_tab / a_fc - 1):.1e})",
              abs(A_tab / a_fc - 1) <= 1e-6)
    else:
        SKIPPED.append(f"3D STEP/STL (freecadcmd {FREECADCMD} か {step2d} が無い)")
    # ---------------- 3. gmsh の読み直し ----------------
    if fc is not None and Path(MESH_PY).exists():
        q = subprocess.run([MESH_PY, "-I", str(GMSH_CHECK), str(step3d), "--out", str(out / "cad" / "gmsh_check.json")],
                           capture_output=True, text=True, timeout=600)
        g = json.loads((out / "cad" / "gmsh_check.json").read_text()) if (out / "cad" / "gmsh_check.json").exists() else {}
        check(f"gmsh: 読めて面 {g.get('n_surfaces')} 枚・面積の求積が収束 (exit {q.returncode})",
              q.returncode == 0 and g.get("n_surfaces", 0) >= 1 and g.get("area_quadrature_converged"))
        a_g = float(g.get("area_mm2", float("nan")))
        INFO["gmsh_area_mm2"] = a_g
        INFO["gmsh_occ_mass_mm2 (参考)"] = g.get("area_occ_mass_mm2")
        check(f"gmsh: 面積 {a_g:.6f} mm² が ∫2πr ds と相対 {abs(a_g / ref_area - 1):.1e} ≤ 1e-6 (V2)", abs(a_g / ref_area - 1) <= 1e-6)
        print(f"     参考: gmsh の occ.getMass は {g.get('area_occ_mass_mm2'):.3f} mm² (相対 {g.get('area_occ_mass_mm2') / ref_area - 1:+.2e})")
        bb = g.get("bbox_mm") or [np.nan] * 6
        ex = max(abs(bb[0] - x0 * s_mm), abs(bb[3] - x1 * s_mm))
        INFO["gmsh_bbox_x_err_mm"] = ex
        check(f"gmsh: bbox の x 範囲が定義域と {ex:.1e} mm ≤ 1e-6 (V2)", ex <= 1e-6)
        check(f"形状表の累積側面積の最終値と gmsh の面積が相対 {abs(A_tab / a_g - 1):.1e} ≤ 1e-6 (V5)", abs(A_tab / a_g - 1) <= 1e-6)
    else:
        SKIPPED.append(f"gmsh の読み直し ({MESH_PY} か 3D STEP が無い)")

    # ---------------- 4. CFD 表 case/45 (断熱壁) ----------------
    if RUN45.is_dir():
        c45 = out / "c45"
        st = WT.cfd_wall_table(RUN45, c45 / "wall_cfd_nodes.csv")
        hdr, rows = read_csv(c45 / "wall_cfd_nodes.csv")
        check(f"CFD 表 (case/45 断熱壁): 状態 {st['status']} = missing (iface_* なし)・step {st['step']}",
              st["status"] == "missing" and st["step"] == 12000 and st["walls"]["3"]["heat"]["status"] == "missing")
        check("CFD 表: 正本の列名 q_into_wall_W_m2・Q_node_W で中身は空欄 (診断値に置き換えない)",
              "q_into_wall_W_m2" in hdr and "q_into_wall_diag_W_m2" not in hdr
              and all(r[hdr.index("q_into_wall_W_m2")] == "" and r[hdr.index("Q_node_W")] == "" for r in rows))
        check(f"CFD 表: 列名が ASCII・x 昇順 ({len(rows)} 行)・.status.json がある",
              all(c.isascii() for c in hdr) and np.all(np.diff(col(hdr, rows, "x_m")) > 0)
              and (c45 / "wall_cfd_nodes.status.json").exists())
        sd = WT.cfd_wall_table(RUN45, c45 / "wall_cfd_nodes_diag.csv", allow_diagnostic=True)
        hd, rd_ = read_csv(c45 / "wall_cfd_nodes_diag.csv")
        with h5py.File(RUN45 / "res_wall_3_12000.h5") as f:
            qw = np.array(f["VALUE/qwall"], dtype=float)
            wx = np.array(f["MESH/COORD"], dtype=float).reshape(-1, 3)
        qd = col(hd, rd_, "q_into_wall_diag_W_m2")
        check(f"CFD 表 (allow_diagnostic): 状態 {sd['status']}・列名 q_into_wall_diag_W_m2・Q_node_diag_W に切り替わり、値 = −qwall",
              sd["status"] == "diagnostic" and "q_into_wall_diag_W_m2" in hd and "Q_node_diag_W" in hd and "q_into_wall_W_m2" not in hd
              and np.array_equal(qd, -qw[np.argsort(wx[:, 0])]))
        ring = col(hdr, rows, "ring_area_m2")
        with h5py.File(RUN45 / "nozzle.h5") as f:
            nod = np.array(f["MESH/COORD"], dtype=float).reshape(-1, 3)
            e = np.array(f["BCONDS/3/vizBfaceNodes"]).reshape(-1, 2)
        L = np.hypot(*(nod[e[:, 1], :2] - nod[e[:, 0], :2]).T)
        poly = float(np.sum(np.pi * (nod[e[:, 0], 1] + nod[e[:, 1], 1]) * L))
        INFO["c45_ring_area_sum_m2"] = float(ring.sum())
        INFO["c45_polygon_lateral_area_m2"] = poly
        INFO["c45_ring_vs_polygon_rel"] = abs(ring.sum() / poly - 1)
        check(f"CFD 表: ring_area の総和 {ring.sum():.15g} m² と壁の多角形の側面積 (メッシュの辺から) {poly:.15g} が相対 {abs(ring.sum() / poly - 1):.1e} ≤ 1e-12 (V5)",
              abs(ring.sum() / poly - 1) <= 1e-12)
        y1 = col(hdr, rows, "y1plus")
        cl, txt = cli_y1plus(RUN45, "wall", 12000)
        INFO["c45_y1plus_mean_max_table"] = [float(np.nanmean(y1)), float(np.nanmax(y1))]
        INFO["c45_y1plus_mean_max_cli"] = cl
        check(f"CFD 表: y1+ (平均 {np.nanmean(y1):.4f}・最大 {np.nanmax(y1):.4f}) が check_wall_resolution.py の CLI ({cl}) と一致 (表示の桁)",
              cl is not None and abs(np.nanmean(y1) - cl[0]) <= 5e-4 + 1e-9 and abs(np.nanmax(y1) - cl[1]) <= 5e-4 + 1e-9)
        tau = col(hdr, rows, "tau_w_Pa")
        check(f"CFD 表: tau_w は接線成分で付着流で正 (正の割合 {np.mean(tau > 0):.3f})", np.mean(tau > 0) > 0.99)
        th = WT.thermal_use(RUN45)
        INFO["c45_thermal_use"] = th["verdict"]
        check(f"thermal_use (未収束・断熱の case/45): {th['verdict']} ≠ PASS (収束 {th['convergence']['result']}: {th['convergence'].get('verdict_line')})",
              th["verdict"] == "FAIL" and th["convergence"]["result"] == "FAIL" and th["heat_status"] == "missing")
    else:
        SKIPPED.append(f"case/45 の run ({RUN45}) が無い")

    # ---------------- 5. CFD 表 case/63 (軸対称・等温) ----------------
    if R63_SMOKE.is_dir():
        s4 = WT.cfd_wall_table(R63_SMOKE, out / "c63s" / "wall_cfd_nodes.csv", step=8, wall_ids=[3])
        h = s4["walls"]["3"]["heat"]
        check(f"CFD 表 (case/63 run_0004 smoke、wall_up step 8): 状態 {s4['status']} = valid (iface_Qf_eff {h.get('n_finite')}/{h.get('n')} 有限)",
              s4["status"] == "valid" and h.get("n_finite") == h.get("n") == 25)
        hdr, rows = read_csv(out / "c63s" / "wall_cfd_nodes.csv")
        q = col(hdr, rows, "q_into_wall_W_m2")
        with h5py.File(R63_SMOKE / "res_wall_up_3_8.h5") as f:
            qw = -np.array(f["VALUE/qwall"], dtype=float)[np.argsort(np.array(f["MESH/COORD"]).reshape(-1, 3)[:, 0])]
        INFO["c63_smoke_q_median_vs_minus_qwall_median"] = [float(np.median(q)), float(np.median(qw))]
        print(f"     参考 (ΔT 0・8 step で未収束、物理の熱流束は ~0): 表の q の中央値 {np.median(q):.3g} W/m²、−qwall の中央値 {np.median(qw):.3g} W/m²、"
              f"符号の一致 {np.mean(np.sign(q) == np.sign(qw)):.2f}")
    else:
        SKIPPED.append(f"case/63 の smoke run ({R63_SMOKE}) が無い")
    if R63_CONV.is_dir():
        sc = WT.cfd_wall_table(R63_CONV, out / "c63" / "wall_cfd_nodes.csv", wall_ids=[4])
        hdr, rows = read_csv(out / "c63" / "wall_cfd_nodes.csv")
        xq, q, Q, A = (col(hdr, rows, c) for c in ("x_m", "q_into_wall_W_m2", "Q_node_W", "A_node_m2"))
        with h5py.File(R63_CONV / "res_wall_heat_4_600000.h5") as f:
            o = np.argsort(np.array(f["MESH/COORD"]).reshape(-1, 3)[:, 0])
            qw = -np.array(f["VALUE/qwall"], dtype=float)[o]
            qc = np.array(f["VALUE/iface_q_compact"], dtype=float)[o]
        check(f"CFD 表 (case/63 run_0013 収束・ΔT 10 K、wall_heat): 状態 {sc['status']} = valid", sc["status"] == "valid")
        sel = [int(np.argmin(np.abs(xq - xs))) for xs in (0.01, 0.05, 0.15)]
        rel_w = [abs(q[i] / qw[i] - 1) for i in sel]
        rel_c = [abs(q[i] / qc[i] - 1) for i in sel]
        INFO["c63_q_vs_minus_qwall_rel_at_x(0.01,0.05,0.15)"] = rel_w
        INFO["c63_q_vs_q_compact_rel"] = rel_c
        check(f"CFD 表: q と −qwall・iface_q_compact が内部 3 箇所 (x = 0.01/0.05/0.15 m) で 5 % 以内 ({[f'{v:.3f}' for v in rel_w]})",
              max(rel_w) <= 0.05 and max(rel_c) <= 0.05)
        check(f"CFD 表: 符号が −qwall と全節点で一致 (流体→壁が正; 加熱壁なので負 {np.mean(q < 0):.2f})",
              np.all(np.sign(q) == np.sign(qw)))
        tot = float(np.sum(Q))
        check(f"CFD 表: ΣQ_node {tot:.6g} W と Σ q A_node が相対 {abs(np.sum(q * A) / tot - 1):.1e} ≤ 1e-9", abs(np.sum(q * A) / tot - 1) <= 1e-9)
        ends = [(float(xq[i]), float(q[i] / qw[i])) for i in (0, len(q) - 1)]
        INFO["c63_end_nodes_q_over_minus_qwall"] = ends
        print(f"     参考: 端節点 (隣の断熱壁と共有) の q / (−qwall) = {ends}")
        y1 = col(hdr, rows, "y1plus")
        cl, _ = cli_y1plus(R63_CONV, "wall_heat", 600000)
        check(f"CFD 表 (case/63): y1+ (平均 {np.nanmean(y1):.4f}・最大 {np.nanmax(y1):.4f}) が CLI ({cl}) と一致",
              cl is not None and abs(np.nanmean(y1) - cl[0]) <= 5e-4 + 1e-9 and abs(np.nanmax(y1) - cl[1]) <= 5e-4 + 1e-9)
        sa = WT.cfd_wall_table(R63_CONV, out / "c63" / "wall_up.csv", wall_ids=[3])
        check(f"CFD 表 (case/63 run_0013 の断熱壁 wall_up): 状態 {sa['status']} = not_applicable (データセットはあるが全 NaN)",
              sa["status"] == "not_applicable")
        xl = WT.write_xlsx(out / "excel" / "nozzle_tables.xlsx", out / "excel" / "wall_shape_mm.csv", out / "c63" / "wall_cfd_nodes.csv",
                           {"Pt_Pa": 5.5e6, "gas": {"preset": "air"}}, ["試験の注記"])
        if xl.get("skipped"):
            SKIPPED.append(f"xlsx ({xl.get('reason')})")
            check("xlsx: 飛ばしたことを返す (written False・理由つき)", xl.get("written") is False and bool(xl.get("reason")))
        else:
            import zipfile
            with zipfile.ZipFile(out / "excel" / "nozzle_tables.xlsx") as z:
                wb = z.read("xl/workbook.xml").decode()
                names = re.findall(r'sheet name="([^"]+)"', wb)
                s5 = z.read("xl/worksheets/sheet5.xml").decode()
            check(f"xlsx: シート {names}", names == ["shape", "cfd_wall", "conditions", "notes", "wall_1d"])
            qi = hdr.index("q_into_wall_W_m2")
            qL = chr(ord("A") + qi)
            check(f"xlsx: wall_1d の式が cfd_wall の q 列 ({qL}) を参照し T_back = T_w − q t / k",
                  f"cfd_wall!{qL}2" in s5 and "C8-D8*$B$1/$B$2" in s5 and xl["info"]["wall_1d"]["q_column"] == "q_into_wall_W_m2")
        xs = WT.write_xlsx(out / "excel" / "x.xlsx", out / "excel" / "wall_shape_mm.csv", None, {}, [], python=str(out / "no_python"))
        check("xlsx: python が無ければ飛ばしたことを返す (例外にしない)", xs.get("skipped") is True and xs.get("written") is False)
    else:
        SKIPPED.append(f"case/63 の収束 run ({R63_CONV}) が無い")
    if R63_DRY.is_dir():
        sdry = WT.cfd_wall_table(R63_DRY, out / "c63d" / "wall_up.csv", wall_ids=[3])
        check(f"CFD 表 (case/63 run_0001_dry の断熱壁 wall_up): 状態 {sdry['status']} = not_applicable", sdry["status"] == "not_applicable")

    # ---------------- 6. 合成データ ----------------
    syn = out / "syn"
    tr = make_synth(syn / "axi")
    s = WT.cfd_wall_table(syn / "axi", syn / "axi_out" / "t.csv")
    hdr, rows = read_csv(syn / "axi_out" / "t.csv")
    wid = col(hdr, rows, "wall_id")
    q, Q, A, ring, tau, y1 = (col(hdr, rows, c) for c in ("q_into_wall_W_m2", "Q_node_W", "A_node_m2", "ring_area_m2", "tau_w_Pa", "y1plus"))
    xx = col(hdr, rows, "x_m")
    ok_q, ok_Q, ok_old, ok_ring = [], [], [], []
    for pid, t in tr.items():
        m = wid == pid
        o = np.argsort(xx[m])
        qq, QQ, AA = q[m][o], Q[m][o], A[m][o]
        ok_q.append(float(np.max(np.abs(qq / t["q"] - 1))))
        ok_Q.append(float(np.max(np.abs(QQ / (2 * np.pi * t["Qf"]) - 1))))
        ok_old.append(float(np.max(np.abs(t["q_eff_old"] / qq - t["r"]))))
        ok_ring.append(float(np.max(np.abs(ring[m][o] / AA - 1))))
    INFO["synthetic_q_rel_err"] = max(ok_q)
    check(f"合成 (r = 0.01・0.05 m): 状態 {s['status']}・表の q が与えた q に相対 {max(ok_q):.1e} ≤ 1e-12 で戻る", s["status"] == "valid" and max(ok_q) <= 1e-12)
    check(f"合成: Q_node = 2π iface_Qf_eff (相対 {max(ok_Q):.1e})・円筒では ring_area = A_node ({max(ok_ring):.1e})",
          max(ok_Q) <= 1e-12 and max(ok_ring) <= 1e-12)
    check(f"合成 (負例): 旧定義 iface_q_eff (= r q) をそのまま使うと r 倍ずれる (iface_q_eff / q − r の最大 {max(ok_old):.1e})、"
          f"表の q とは違う", max(ok_old) <= 1e-12 and not np.allclose(tr[3]["q_eff_old"], tr[3]["q"], rtol=0.5))
    yref = SYN_WALLS[3][3] * math.sqrt(SYN_RO * 5.0) / SYN_VISC
    check(f"合成: tau_w = 5 Pa (接線、+x)・y1+ = y1 √(ρ τ)/μ = {yref:.6g} (壁 a)",
          np.allclose(tau, 5.0, rtol=1e-15) and abs(y1[wid == 3][0] / yref - 1) <= 1e-12)
    check("合成: x 昇順 (壁ごと、壁は始点の x の順)", np.all(np.diff(xx[wid == 3]) > 0) and np.all(np.diff(xx[wid == 4]) > 0) and xx[wid == 3].max() < xx[wid == 4].min())
    s_one = WT.cfd_wall_table(syn / "axi", syn / "axi_out" / "one.csv", wall_ids=[4])
    check("合成: 壁を絞れる (wall_ids)", s_one["status"] == "valid" and list(s_one["walls"]) == ["4"])
    make_synth(syn / "nan1", qf="nan_one")
    s = WT.cfd_wall_table(syn / "nan1", syn / "nan1_out" / "t.csv")
    hdr, rows = read_csv(syn / "nan1_out" / "t.csv")
    check(f"合成 (一部 NaN): 状態 {s['status']} = nonfinite・q 列は空欄",
          s["status"] == "nonfinite" and all(r[hdr.index("q_into_wall_W_m2")] == "" for r in rows))
    make_synth(syn / "nanall", qf="nan_all")
    check("合成 (全 NaN): not_applicable", WT.cfd_wall_table(syn / "nanall", syn / "nanall_out" / "t.csv")["status"] == "not_applicable")
    make_synth(syn / "qeff", only_q_eff=True)
    s = WT.cfd_wall_table(syn / "qeff", syn / "qeff_out" / "t.csv")
    sd = WT.cfd_wall_table(syn / "qeff", syn / "qeff_out" / "d.csv", allow_diagnostic=True)
    hd, rd_ = read_csv(syn / "qeff_out" / "d.csv")
    check(f"合成 (iface_q_eff だけ): 状態 {s['status']} = missing (iface_q_eff に置き換えない)・診断を許すと diagnostic で −qwall",
          s["status"] == "missing" and sd["status"] == "diagnostic"
          and np.allclose(col(hd, rd_, "q_into_wall_diag_W_m2")[col(hd, rd_, "wall_id") == 3], 1.01 * tr[3]["q"], rtol=1e-12))
    make_synth(syn / "mix", steps=(100, 200))
    (syn / "mix" / "res_wall_a_3_200.h5").unlink()
    e1 = raises(lambda: WT.cfd_wall_table(syn / "mix", syn / "mix_out" / "t.csv"), ValueError)
    e2 = raises(lambda: WT.cfd_wall_table(syn / "mix", syn / "mix_out" / "t.csv", step=200), FileNotFoundError)
    check(f"合成 (異なる step の混入): 最終 step が壁で違えば例外 ({e1})・欠けた step の指定も例外", e1 is not None and e2 is not None)
    check("合成: step を揃えて指定すれば読める", WT.cfd_wall_table(syn / "mix", syn / "mix_out" / "t.csv", step=100)["step"] == 100)
    shutil.copy2(syn / "mix" / "res_wall_b_4_100.h5", syn / "mix" / "res_wall_b_4_300.h5")
    shutil.copy2(syn / "mix" / "res_wall_a_3_100.h5", syn / "mix" / "res_wall_a_3_300.h5")
    e3 = raises(lambda: WT.cfd_wall_table(syn / "mix", syn / "mix_out" / "t.csv", step=300), ValueError)
    check(f"合成: ファイル名の step と属性 step が違えば例外 ({e3})", e3 is not None and "属性 step" in e3)
    trp = make_synth(syn / "planar", axi=False)
    s = WT.cfd_wall_table(syn / "planar", syn / "planar_out" / "t.csv")
    hdr, rows = read_csv(syn / "planar_out" / "t.csv")
    m = col(hdr, rows, "wall_id") == 3
    o = np.argsort(col(hdr, rows, "x_m")[m])
    qp = col(hdr, rows, "q_into_wall_W_m2")[m][o]
    Qp = col(hdr, rows, "Q_node_W")[m][o]
    check(f"合成 (平面): Q_node = iface_Qf_eff・A = |S| で q が戻る ({np.max(np.abs(qp / trp[3]['q'] - 1)):.1e})・状態 {s['status']}",
          s["status"] == "valid" and np.max(np.abs(qp / trp[3]["q"] - 1)) <= 1e-12 and np.max(np.abs(Qp / trp[3]["Qf"] - 1)) <= 1e-12)

    # ---------------- 7. 時系列と thermal_use (合成) ----------------
    make_synth(syn / "pass", steps=(100, 200, 300, 400, 500), tw=1e-5, residual=True)
    ser = WT.wall_heat_series(syn / "pass", syn / "pass_out" / "series.csv")
    hs, rs = read_csv(syn / "pass_out" / "series.csv")
    tot_ref = sum(2 * np.pi * np.sum(t["Qf"]) for t in tr.values())
    check(f"wall_heat_series: 列 {hs}・{ser['n_steps']} step・総熱量 = Σ 2π iface_Qf_eff (相対 {abs(float(rs[0][1]) / tot_ref - 1):.1e})",
          hs[:2] == ["step", "Q_total_W"] and ser["n_steps"] == 5 and abs(float(rs[0][1]) / tot_ref - 1) <= 1e-12)
    th = WT.thermal_use(syn / "pass")
    INFO["synthetic_thermal_use_pass"] = {k: th[k]["result"] for k in ("convergence", "quasisteady", "wall_resolution")}
    check(f"thermal_use (合成: 収束・定常・y1+ < 1・valid): {th['verdict']} = PASS ({INFO['synthetic_thermal_use_pass']})",
          th["verdict"] == "PASS" and th["quasisteady"]["verdict_line"] is not None and th["wall_resolution"]["verdict_line"].startswith("VERDICT: PASS"))
    th2 = WT.thermal_use(syn / "pass", step=300)
    check(f"thermal_use: 最終でない step は収束の区間を切り出さず未評価 ({th2['verdict']}: {th2['convergence'].get('reason')})",
          th2["verdict"] == "未評価" and th2["convergence"]["result"] == "未評価")
    make_synth(syn / "nomat", tw=1e-5)
    th3 = WT.thermal_use(syn / "nomat")
    check(f"thermal_use (合成: 残差履歴なし・ダンプ 1 枚): {th3['verdict']} = 未評価 (判定材料なし; "
          f"収束 {th3['convergence']['result']}・準定常 {th3['quasisteady']['result']}・壁解像 {th3['wall_resolution']['result']})",
          th3["verdict"] == "未評価" and th3["convergence"]["result"] == "未評価" and th3["quasisteady"]["result"] == "未評価")
    th5 = WT.thermal_use(syn / "axi")
    check(f"thermal_use (合成: y1+ 136 で壁解像 FAIL): {th5['verdict']} = FAIL ({th5['wall_resolution']['verdict_line']})",
          th5["verdict"] == "FAIL" and th5["wall_resolution"]["result"] == "FAIL")
    make_synth(syn / "pass_nan", steps=(100, 200, 300, 400, 500), tw=1e-5, residual=True, qf="nan_one")
    th4 = WT.thermal_use(syn / "pass_nan")
    check(f"thermal_use (合成: 熱量 nonfinite): {th4['verdict']} = FAIL", th4["verdict"] == "FAIL")


if __name__ == "__main__":
    main()
