"""構造ファイル `export/mesh_structure.py` と STEP→gmsh 再構築 `export/rebuild_from_step.py` の試験
(plan tooling-nozzle-wt-deliverables §5.1 #4・§6 V3, 2026-10-13)。

- 合成メッシュ (ni 21・nj 9): msh4.1 の読み込み → 構造 → 往復 (構造から座標を作り直して一致)。
- 負例: 1 節点の x を 1e-6 L ずらした msh (縦線)・節点番号を入れ替えた msh (接続)・軸から浮いた節点・書式違い。
- M6 fixture (case/45 の生産の NS 格子 2000×97 と同じ壁の STEP): 再構築の比較 JSON が PASS。
  助変数を付け直した STEP (点列から occ.addBSpline、x(u) ≠ u) でも節点が作れ、曲線を密に評価した r(x) と壁節点が ≤ 1e-5 mm。
  折り返す曲線・station を含まない曲線は拒否。station より長い曲線は切って使う。
- 出力の形: `OUT_PREFIX.msh` は forge の変換器が読める形 (write_msh41_2d と同じ書式; 10 桁ならバイト一致)、`_cad.msh` は違う。
- V3(d): 変換器 (FORGE_CONVERTER、既定は solver_density_cuda/.build-native/relwithdebinfo/convertGmshToForge) があれば、
  元の msh と作り直した `.msh` を cell の config (FORGE_QC_CONFIG の solverConfig.yaml を cell に直したもの) で変換し、
  `check_mesh_quality.py --ar-max 5000` がどちらも PASS。作業はテンポラリ (run ディレクトリには書かない)。

fixture は環境変数 FORGE_DELIV_FIXTURE (既定は下の FIXTURE_DEFAULT)、gmsh の python は FORGE_MESH_PYTHON
(既定は .venv-mesh)。無ければ該当の試験を飛ばしたと出力する。作業フォルダは終わりに消す (FORGE_KEEP_TMP=1 で残す)。
usage: design/.venv-opt/bin/python design/tests/run_mesh_structure_tests.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from forge_design.export.mesh_structure import (  # noqa: E402
    StructureError, extract_structure, load_structure_json, read_msh41_mesh2d, sha256_file, structure_coords,
    structure_from_msh, structure_s, write_structure_json)
from forge_design.meshing.mesh2d import PHYS, Mesh2DParams, _radial_fracs, generate_axisym_mesh, write_msh41_2d  # noqa: E402
# 単体スクリプトだが numpy しか import しないので書き出し器と形の検査を直接試せる
from forge_design.export.rebuild_from_step import forge_format_check, write_msh41_forge  # noqa: E402

FIXTURE_DEFAULT = "/tmp/claude-1000/-home-sano-work-forge/bd18adac-85f4-4627-811b-8f3186867995/scratchpad/fixtures/m6_prod"
REBUILD = HERE.parents[1] / "forge_design" / "export" / "rebuild_from_step.py"
CONVERTER_DEFAULT = HERE.parents[2] / "solver_density_cuda" / ".build-native" / "relwithdebinfo" / "convertGmshToForge"
QC_TOOL = HERE.parents[2] / "solver_density_cuda" / "tools" / "check_mesh_quality.py"
QC_CONFIG_DEFAULT = "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0051_ns_final_c2/solverConfig.yaml"
fails = 0


def check(name, ok, info=""):
    global fails
    print(("ok   " if ok else "FAIL ") + name + (f" ({info})" if info != "" else "")); fails += (not ok)


def skip(name, why):
    print(f"skip {name} ({why})")


def mesh_python():
    cands = [os.environ.get("FORGE_MESH_PYTHON"), HERE.parents[2] / ".venv-mesh" / "bin" / "python",
             "/home/sano/work/forge/.venv-mesh/bin/python"]
    for c in cands:
        if c and Path(c).exists():
            return str(c)
    return None


def raises(fn, exc, word):
    try:
        fn()
    except exc as e:
        return word in str(e), str(e)
    return False, "例外なし"


class Wall:
    """合成の壁 r(x) = 1 + 0.1 x² [r*] (スロート x = 0)。"""
    x_in, x_e = -2.0, 6.0

    @staticmethod
    def r(x):
        return 1.0 + 0.1 * np.asarray(x, dtype=float) ** 2


def as_arrays(bedges):
    return {k: np.asarray(v, dtype=np.int64) for k, v in bedges.items()}


tmp = Path(tempfile.mkdtemp(prefix="mesh_structure_tests_"))
print(f"作業フォルダ: {tmp}")

# ------------------------------------------------------------------------------------------
# 1. 合成メッシュ: 読み込み → 構造 → 往復
# ------------------------------------------------------------------------------------------
prm = Mesh2DParams(ni=21, nj=9, wall_first_frac=0.05, throat_refine=2.0, scale=0.05)
coords, quads, bedges = generate_axisym_mesh(Wall, prm)
syn = tmp / "syn.msh"
write_msh41_2d(syn, coords, quads, bedges)
c1, q1, b1, p1 = read_msh41_mesh2d(syn)
check("(1a) 読み込み: 四角形・境界辺・物理名が書いたものと一致",
      np.array_equal(q1, quads) and all(np.array_equal(b1[k], np.asarray(v)) for k, v in bedges.items())
      and list(b1) == ["inlet", "outlet", "wall", "axis"] and p1 == PHYS)
check("(1a) 読み込み: 座標が 10 有効桁の丸めの範囲で一致", c1.shape == (21 * 9, 3) and c1.dtype == np.float64
      and np.abs(c1 - coords).max() <= 1e-9 * np.abs(coords).max(), f"{np.abs(c1 - coords).max():.2e} m")
st = structure_from_msh(syn)
check("(1b) ni・nj の推定 (壁辺数 + 1・入口辺数 + 1)", (st["ni"], st["nj"]) == (21, 9))
check("(1b) format・connectivity・物理名とタグ・sha256",
      st["format"] == "forge.mesh_structure/v1" and st["connectivity"] == "mesh2d_ij_v1"
      and {k: v["tag"] for k, v in st["physical"].items()} == PHYS and st["physical"]["fluid"]["dim"] == 2
      and st["source"]["sha256"] == sha256_file(syn))
check("(1b) x_mm・壁半径・第 1 層厚の基準値",
      np.allclose(st["x_mm"], c1[::9, 0] * 1e3, rtol=0, atol=1e-12)
      and np.allclose(st["wall_r_mm"], c1[8::9, 1] * 1e3, rtol=0, atol=1e-12)
      and np.allclose(st["first_layer_mm"], (c1[8::9, 1] - c1[7::9, 1]) * 1e3, rtol=0, atol=1e-12))
check("(1b) msh の丸めで行が揃わないので s は station ごと (ni×nj)", "s" in st and "s_common" not in st)
js = tmp / "syn_structure.json"
write_structure_json(st, js)
st2 = load_structure_json(js)
rt = structure_coords(st2)
L = c1[-9, 0] - c1[0, 0]
check("(1c) 往復: JSON から作り直した座標が読んだ msh と一致 (≤ 1e-14 L)", np.abs(rt - c1).max() <= 1e-14 * L,
      f"{np.abs(rt - c1).max() / L:.2e} L")
st_d = extract_structure(coords, quads, bedges)
check("(1d) 丸めのない配列では s_common 1 行 (= _radial_fracs)", "s_common" in st_d
      and np.abs(np.asarray(st_d["s_common"]) - _radial_fracs(9, 0.05)).max() <= 1e-15)
check("(1d) s_common の往復", np.abs(structure_coords(st_d) - coords).max() <= 1e-15 * 10)
prm_v = Mesh2DParams(ni=21, nj=9, wall_first_frac=0.05, wall_first_frac_throat=0.01, wall_first_blend_x0=0.5,
                     wall_first_blend_x1=4.0, scale=0.05)
cv, qv, bv = generate_axisym_mesh(Wall, prm_v)
st_v = extract_structure(cv, qv, bv)
check("(1e) 第 1 層厚が x で変わる格子は s が ni×nj で往復一致", "s" in st_v
      and structure_s(st_v).shape == (21, 9) and np.abs(structure_coords(st_v) - cv).max() <= 1e-15 * 10)
ok, msg = raises(lambda: extract_structure(c1, q1, b1, ni=20), StructureError, "ni")
check("(1f) 与えた ni が推定と違えば拒否", ok, msg)

# ------------------------------------------------------------------------------------------
# 2. 負例
# ------------------------------------------------------------------------------------------
cx = coords.copy()
Lm = coords[-1, 0] - coords[0, 0]
cx[10 * 9 + 4, 0] += 1e-6 * Lm
write_msh41_2d(tmp / "neg_x.msh", cx, quads, bedges)
ok, msg = raises(lambda: structure_from_msh(tmp / "neg_x.msh"), StructureError, "縦線")
check("(2a) 1 節点の x を 1e-6 L ずらした msh は縦線の検査で拒否", ok, msg[:90])
cs = coords.copy()
cs[10 * 9 + 4, 0] += 0.5e-9 * Lm
check("(2a') 許容 (1e-9 L) 内のずれは通す (配列で)", extract_structure(cs, quads, bedges)["ni"] == 21)
cs[10 * 9 + 4, 0] += 1.0e-9 * Lm
ok, msg = raises(lambda: extract_structure(cs, quads, bedges), StructureError, "縦線")
check("(2a') 許容を超える 1.5e-9 L は拒否 (配列で)", ok)
rng = np.random.default_rng(0)
perm = np.arange(len(coords))
perm[[30, 31]] = perm[[31, 30]]                                     # 2 節点の番号を入れ替える (形は同じ)
cp = np.empty_like(coords)
cp[perm] = coords
bp = {k: perm[np.asarray(v)] for k, v in bedges.items()}
write_msh41_2d(tmp / "neg_perm.msh", cp, perm[quads], bp)
ok, msg = raises(lambda: structure_from_msh(tmp / "neg_perm.msh"), StructureError, "接続")
check("(2b) 節点番号を入れ替えた msh (形は同じ) は接続の検査で拒否", ok, msg[:90])
perm2 = rng.permutation(len(coords))
cp2 = np.empty_like(coords)
cp2[perm2] = coords
write_msh41_2d(tmp / "neg_perm2.msh", cp2, perm2[quads], {k: perm2[np.asarray(v)] for k, v in bedges.items()})
ok, msg = raises(lambda: structure_from_msh(tmp / "neg_perm2.msh"), StructureError, "接続")
check("(2b') 番号を乱数で並べ替えた msh も拒否", ok)
ca = coords.copy()
ca[5 * 9 + 0, 1] = 1e-6
write_msh41_2d(tmp / "neg_axis.msh", ca, quads, bedges)
ok, msg = raises(lambda: structure_from_msh(tmp / "neg_axis.msh"), StructureError, "軸")
check("(2c) 軸の節点が r = 0 でない msh は拒否", ok, msg[:90])
txt = syn.read_text()
(tmp / "neg_fmt.msh").write_text(txt.replace("4.1 0 8", "2.2 0 8", 1))
ok, msg = raises(lambda: read_msh41_mesh2d(tmp / "neg_fmt.msh"), ValueError, "4.1")
check("(2d) 書式が 4.1 でない msh は ValueError", ok, msg)
ok, msg = raises(lambda: extract_structure(c1, q1, {k: v for k, v in b1.items() if k != "axis"}), StructureError, "境界辺")
check("(2e) 境界辺が 4 つ揃わなければ拒否", ok)

# 書き出し器 (rebuild_from_step.write_msh41_forge) が write_msh41_2d と同じ並び・書式か
ORDER = ["inlet", "outlet", "wall", "axis"]
write_msh41_forge(tmp / "w10.msh", coords, quads, bedges, PHYS, ORDER, "fluid", digits=10)
check("(2f) forge の形の書き出し器: 同じ配列・10 桁なら write_msh41_2d とバイト一致", (tmp / "w10.msh").read_bytes() == syn.read_bytes())
write_msh41_forge(tmp / "w17.msh", coords, quads, bedges, PHYS, ORDER, "fluid")
c17, q17, b17, p17 = read_msh41_mesh2d(tmp / "w17.msh")
t10 = syn.read_text().splitlines()
t17 = (tmp / "w17.msh").read_text().splitlines()
k0, k1 = t10.index("$Nodes"), t10.index("$EndNodes")
check("(2f) 17 桁: 座標の行以外は 10 桁版と同じ行・座標は元の配列と一致・接続と物理名も一致",
      t10[:k0 + 3 + 21 * 9] == t17[:k0 + 3 + 21 * 9] and t10[k1:] == t17[k1:] and np.array_equal(c17, coords)
      and np.array_equal(q17, quads) and all(np.array_equal(b17[k], np.asarray(v)) for k, v in bedges.items()) and p17 == PHYS)
ff = forge_format_check(tmp / "w17.msh")
check("(2f) forge の形の検査: 節点ブロック 1 つ・節点番号と要素番号がファイル順に 1..N", ff["ok"], ff)

# ------------------------------------------------------------------------------------------
# 3. gmsh を使う試験 (.venv-mesh)
# ------------------------------------------------------------------------------------------
HELPER = r'''
import sys, json, numpy as np, gmsh
mode = sys.argv[1]
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
def curve(step):
    gmsh.model.occ.importShapes(step); gmsh.model.occ.synchronize()
    t = gmsh.model.getEntities(1)[0][1]
    lo, hi = (float(v[0]) for v in gmsh.model.getParametrizationBounds(1, t))
    return t, lo, hi
if mode == "reparam":            # 元の STEP の曲線を不等間隔の点列で標本化し、点列を制御点にした B-spline で書き直す
    src, out, n = sys.argv[2], sys.argv[3], int(sys.argv[4])
    t, lo, hi = curve(src)
    s = np.linspace(0.0, 1.0, n)
    s = s + 0.3 * np.sin(2 * np.pi * s) / (2 * np.pi)       # 単調な不等間隔 (x(u) を非線形に)
    P = np.asarray(gmsh.model.getValue(1, t, (lo + (hi - lo) * s).tolist())).reshape(-1, 3)
    gmsh.clear(); gmsh.model.add("re")
    pts = [gmsh.model.occ.addPoint(*p) for p in P]
    c = gmsh.model.occ.addBSpline(pts)
    gmsh.model.occ.synchronize()
    gmsh.write(out)
elif mode == "spline":           # 点列 (x, r) [mm] を通る補間スプラインを STEP に
    pts_xy, out = json.loads(open(sys.argv[2]).read()), sys.argv[3]
    pts = [gmsh.model.occ.addPoint(x, r, 0.0) for x, r in pts_xy]
    gmsh.model.occ.addSpline(pts); gmsh.model.occ.synchronize()
    gmsh.write(out)
elif mode == "dense":            # STEP の曲線を密に評価
    step, n, out = sys.argv[2], int(sys.argv[3]), sys.argv[4]
    t, lo, hi = curve(step)
    u = np.linspace(lo, hi, n)
    P = np.asarray(gmsh.model.getValue(1, t, u.tolist())).reshape(-1, 3)
    np.savez(out, u=u, x=P[:, 0], r=P[:, 1], z=P[:, 2], lo=lo, hi=hi)
elif mode == "eval":             # STEP の曲線を与えた助変数で評価
    step, uj, out = sys.argv[2], sys.argv[3], sys.argv[4]
    t, lo, hi = curve(step)
    u = np.asarray(json.loads(open(uj).read()), dtype=float)
    P = np.asarray(gmsh.model.getValue(1, t, np.clip(u, lo, hi).tolist())).reshape(-1, 3)
    np.savez(out, x=P[:, 0], r=P[:, 1])
elif mode == "nodes":            # msh の節点 (番号順) と壁の曲線の助変数範囲
    msh, out = sys.argv[2], sys.argv[3]
    gmsh.open(msh)
    tg, xyz, _ = gmsh.model.mesh.getNodes()
    o = np.argsort(tg)
    np.savez(out, tags=np.asarray(tg)[o], xyz=np.asarray(xyz).reshape(-1, 3)[o])
gmsh.finalize()
'''
MPY = mesh_python()
fx = Path(os.environ.get("FORGE_DELIV_FIXTURE", FIXTURE_DEFAULT))
has_fx = (fx / "nozzle.msh").exists() and (fx / "step" / "wall_physical.step").exists()


def helper(*args):
    hp = tmp / "gmsh_helper.py"
    if not hp.exists():
        hp.write_text(HELPER)
    r = subprocess.run([MPY, str(hp), *map(str, args)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"helper {args[0]} 失敗: {r.stderr[-2000:]}")


def rebuild(step, structure, prefix, compare=None):
    cmd = [MPY, str(REBUILD), str(step), str(structure), str(prefix)] + (["--compare", str(compare)] if compare else [])
    return subprocess.run(cmd, capture_output=True, text=True)


def wall_dense_check(step, prefix, struct, tag):
    """作り直した壁節点の r と、その STEP の曲線を密に評価した r(x) (3 次の局所補間) の差 [mm]。x は station のまま。"""
    from scipy.interpolate import CubicSpline
    helper("dense", step, 400001, tmp / f"{tag}_dense.npz")
    helper("nodes", f"{prefix}.msh", tmp / f"{tag}_nodes.npz")
    d = np.load(tmp / f"{tag}_dense.npz")
    nd = np.load(tmp / f"{tag}_nodes.npz")
    ni, nj = struct["ni"], struct["nj"]
    xyz = nd["xyz"] * 1e3
    xw, rw = xyz[nj - 1::nj, 0], xyz[nj - 1::nj, 1]
    xs, rs = d["x"], d["r"]
    if xs[-1] < xs[0]:
        xs, rs = xs[::-1], rs[::-1]
    r_ref = CubicSpline(xs, rs)(np.clip(xw, xs[0], xs[-1]))
    return {"max_wall_r_minus_curve_mm": float(np.abs(rw - r_ref).max()),
            "max_station_x_mm": float(np.abs(xw - np.asarray(struct["x_mm"])).max()),
            "n_nodes": int(len(nd["tags"])), "x_minus_u_max_mm": float(np.abs(d["x"] - d["u"]).max()),
            "param": [float(d["lo"]), float(d["hi"])], "r_wall_mm": rw}


if MPY is None:
    skip("(3)〜(5) gmsh を使う試験", "gmsh の python (.venv-mesh) が見つからない; FORGE_MESH_PYTHON で指定")
else:
    # 3. 合成の構造で、曲線の検査と切り詰め
    xs_mm = np.asarray(st["x_mm"])

    def rmm(x):
        return 50.0 + 0.002 * np.asarray(x) ** 2

    def spline_step(name, xpts):
        jp = tmp / f"{name}_pts.json"
        jp.write_text(json.dumps([[float(x), float(rmm(x))] for x in xpts]))
        helper("spline", jp, tmp / f"{name}.step")
        return tmp / f"{name}.step"

    xf = list(np.linspace(-100.0, 100.0, 21)) + [110.0, 104.0, 120.0] + list(np.linspace(130.0, 300.0, 18))
    r = rebuild(spline_step("fold", xf), js, tmp / "fold_out")
    check("(3a) 折り返す曲線 (x が途中で減る) は拒否 (終了コード 2)", r.returncode == 2 and "狭義単調" in r.stderr,
          r.stderr.strip()[:120])
    r = rebuild(spline_step("short", np.linspace(-90.0, 300.0, 40)), js, tmp / "short_out")
    check("(3b) station の範囲を含まない曲線は拒否", r.returncode == 2 and "含まない" in r.stderr, r.stderr.strip()[:120])
    ext = spline_step("ext", np.linspace(-130.0, 330.0, 93))
    r = rebuild(ext, js, tmp / "ext_out", compare=syn)
    cj = json.loads((tmp / "ext_out_compare.json").read_text()) if (tmp / "ext_out_compare.json").exists() else {}
    check("(3c) station より長い曲線は両端を切って作り直す", r.returncode in (0, 1) and cj.get("stations", {}).get("trimmed") == ["inlet", "outlet"],
          (r.returncode, cj.get("stations", {}).get("trimmed"), r.stderr.strip()[:200]))
    if cj:
        wd = wall_dense_check(ext, tmp / "ext_out", st, "ext")
        check("(3c) 切った曲線でも壁節点が曲線上 (≤ 1e-5 mm)・x は station のまま・節点数",
              wd["max_wall_r_minus_curve_mm"] <= 1e-5 and wd["max_station_x_mm"] <= 1e-12 and wd["n_nodes"] == 21 * 9,
              f"{wd['max_wall_r_minus_curve_mm']:.2e} mm")
        check("(3c) 四角形の接続・物理グループは元と一致 (距離は補間スプラインの近似誤差なので判定外)",
              cj["a_quads"]["match"] and cj["c_physical_groups"]["match"],
              f"最大距離 {cj['b_distance']['max_dist_mm']:.2e} mm")
        rev = spline_step("rev", np.linspace(330.0, -130.0, 93))     # 助変数の向きが x と逆
        r = rebuild(rev, js, tmp / "rev_out", compare=syn)
        cr = json.loads((tmp / "rev_out_compare.json").read_text()) if (tmp / "rev_out_compare.json").exists() else {}
        if cr:
            wr = wall_dense_check(rev, tmp / "rev_out", st, "rev")
            check("(3e) x が u で減る向きの曲線でも作り直せる (壁節点が曲線上・接続と物理グループが一致)",
                  cr["curve"]["direction"] == "x 減少" and wr["max_wall_r_minus_curve_mm"] <= 1e-5
                  and cr["a_quads"]["match"] and cr["c_physical_groups"]["match"],
                  f"{wr['max_wall_r_minus_curve_mm']:.2e} mm・最大距離 {cr['b_distance']['max_dist_mm']:.2e} mm")
        else:
            check("(3e) x が u で減る向きの曲線でも作り直せる", False, r.stderr.strip()[-300:])
        st_rb = structure_from_msh(tmp / "ext_out.msh")
        check("(3d) forge の形の .msh は mesh2d の規則どおり (読み直した構造の ni・nj・station が元と一致)",
              (st_rb["ni"], st_rb["nj"]) == (21, 9) and np.abs(np.asarray(st_rb["x_mm"]) - xs_mm).max() <= 1e-12,
              f"{np.abs(np.asarray(st_rb['x_mm']) - xs_mm).max():.1e} mm")
        try:
            read_msh41_mesh2d(tmp / "ext_out_cad.msh")
            check("(3d) gmsh が書いた _cad.msh (節点ブロックが実体ごと) は read_msh41_mesh2d が ValueError", False)
        except ValueError as e:
            check("(3d) gmsh が書いた _cad.msh (節点ブロックが実体ごと) は read_msh41_mesh2d が ValueError", True, str(e)[:60])
        fc = forge_format_check(tmp / "ext_out_cad.msh")
        check("(3d) forge の形の検査: _cad.msh は節点ブロックが 1 つでなく番号もファイル順でない → ok False",
              not fc["ok"] and fc["n_node_blocks"] > 1 and not fc["node_tags_in_file_order"], fc["n_node_blocks"])

    # 4. M6 fixture
    if not has_fx:
        skip("(4) M6 fixture の試験", f"fixture が無い: {fx} (FORGE_DELIV_FIXTURE で指定)")
    else:
        msh6, step6 = fx / "nozzle.msh", fx / "step" / "wall_physical.step"
        st6 = structure_from_msh(msh6)
        check("(4a) M6: 構造 ni 2000・nj 97・縦線 0・軸 r = 0", (st6["ni"], st6["nj"]) == (2000, 97)
              and st6["checks"]["vertical_max_dev_over_L"] <= 1e-9 and st6["checks"]["axis_max_abs_r_m"] == 0.0,
              f"縦線 {st6['checks']['vertical_max_dev_over_L']:.1e} L")
        js6 = tmp / "m6_structure.json"
        write_structure_json(st6, js6)
        c6, q6, b6, p6 = read_msh41_mesh2d(msh6)
        write_msh41_forge(tmp / "m6_w10.msh", c6, q6, b6, p6, list(b6), "fluid", digits=10)
        check("(4a) M6: 書き出し器に読んだ配列を 10 桁で渡すと元の nozzle.msh とバイト一致",
              (tmp / "m6_w10.msh").read_bytes() == msh6.read_bytes())
        check("(4a) M6: 構造の往復 (≤ 1e-14 L)", np.abs(structure_coords(load_structure_json(js6)) - c6).max()
              <= 1e-14 * st6["length_mm"] * 1e-3)
        r = rebuild(step6, js6, tmp / "m6_out", compare=msh6)
        cj6 = json.loads((tmp / "m6_out_compare.json").read_text()) if (tmp / "m6_out_compare.json").exists() else {}
        bd = cj6.get("b_distance", {})
        check("(4b) M6: 再構築の比較 JSON が PASS (終了コード 0)", r.returncode == 0 and cj6.get("verdict") == "PASS",
              f"最大距離 {bd.get('max_dist_mm', float('nan')):.3e} mm・第 1 層厚比 {bd.get('max_dist_over_first_layer', float('nan')):.3e}"
              f"・第 1 層厚の最小 {bd.get('first_layer_min_mm', float('nan')):.3e} mm" + ("" if r.returncode == 0 else f" / {r.stderr[-300:]}"))
        if cj6:
            check("(4b) M6: 四角形 191904 が全要素で一致・物理グループの名前・タグ・要素数が一致・境界辺も同一",
                  cj6["a_quads"]["match"] and cj6["a_quads"]["n_orig"] == 191904 and cj6["c_physical_groups"]["match"]
                  and cj6["info_group_elements_identical"])
            check("(4b) M6: msh2.2 は番号を振り直すが四角形の角の座標と物理グループの要素数は同じ",
                  cj6["info_v22_readback"]["quad_corner_coords_identical"] and cj6["info_v22_readback"]["group_counts_match"])
            check("(4b) M6: 元の msh の sha256 が構造ファイルと一致", cj6["orig_sha256_matches_structure"] is True)
            check("(4b) M6: (d) .msh が forge の形 (節点ブロック 1 つ・番号がファイル順に 1..N)", cj6["d_forge_format"]["ok"])
            cad = cj6["cad_msh"]
            check("(4b) M6: _cad.msh も (a)(b)(c) を満たす・forge の形ではない (変換器に入れない)",
                  cad["a_quads"]["match"] and cad["c_physical_groups"]["match"]
                  and cad["b_distance"]["max_dist_mm"] <= 1e-5 and cad["b_distance"]["max_dist_over_first_layer"] <= 0.01
                  and not cad["info_forge_format"]["ok"],
                  f"最大距離 {cad['b_distance']['max_dist_mm']:.3e} mm・節点ブロック {cad['info_forge_format']['n_node_blocks']}")
            check("(4b) M6: CGNS を書いた", (tmp / "m6_out.cgns").exists() and (tmp / "m6_out.cgns").stat().st_size > 0)
            print("     比較 JSON:", json.dumps({"verdict": cj6["verdict"], "b_distance": {k: v for k, v in bd.items() if k != "tol"},
                                                  "stations": {k: cj6["stations"][k] for k in ("max_solve_residual_mm", "max_end_residual_mm", "n_iter", "trimmed")},
                                                  "x_minus_u_max_abs_mm": cj6["curve"]["x_minus_u_max_abs_mm"]}, ensure_ascii=False))
        # 助変数の付け直し (正例): 不等間隔の 2001 点を制御点にした B-spline (x(u) ≠ u)
        rep = tmp / "m6_reparam.step"
        helper("reparam", step6, rep, 2001)
        r = rebuild(rep, js6, tmp / "m6_re_out")
        check("(4c) 助変数を付け直した STEP で作り直せる (終了コード 0)", r.returncode == 0, r.stderr.strip()[-300:])
        if r.returncode == 0:
            wd = wall_dense_check(rep, tmp / "m6_re_out", st6, "m6re")
            check("(4c) 付け直した曲線は x(u) ≠ u (助変数の範囲と |x(u) − u|)", wd["x_minus_u_max_mm"] > 1.0,
                  f"助変数 {wd['param'][0]:.3g}〜{wd['param'][1]:.4g}・|x(u) − u| 最大 {wd['x_minus_u_max_mm']:.4g} mm")
            check("(4c) 壁節点の r と曲線を密に評価した r(x) の差 ≤ 1e-5 mm・x は station のまま",
                  wd["max_wall_r_minus_curve_mm"] <= 1e-5 and wd["max_station_x_mm"] <= 1e-9,
                  f"{wd['max_wall_r_minus_curve_mm']:.2e} mm")
            uj = tmp / "m6_x.json"
            uj.write_text(json.dumps(st6["x_mm"]))
            helper("eval", step6, uj, tmp / "m6_orig_at_x.npz")             # 元の STEP は u = x
            ro = np.load(tmp / "m6_orig_at_x.npz")["r"]
            print(f"     (記録) 点列近似そのものの誤差: 付け直した曲線と元の STEP の壁半径の差 最大 "
                  f"{np.abs(wd['r_wall_mm'] - ro).max():.3e} mm (2001 点の制御点による近似; 判定外)")

    # 5. V3(d): forge の変換器 (cell の config) で元と作り直し (forge の形) を変換し、check_mesh_quality --ar-max 5000
    conv = Path(os.environ.get("FORGE_CONVERTER", CONVERTER_DEFAULT))
    cfg_src = Path(os.environ.get("FORGE_QC_CONFIG", QC_CONFIG_DEFAULT))
    if not has_fx:
        skip("(5) V3(d) の品質", "M6 fixture が無い")
    elif not conv.exists():
        skip("(5) V3(d) の品質", f"変換器が無い: {conv} (FORGE_CONVERTER で指定)")
    elif not cfg_src.exists() or not QC_TOOL.exists():
        skip("(5) V3(d) の品質", f"config の元 {cfg_src} か {QC_TOOL} が無い (FORGE_QC_CONFIG で指定)")
    elif not (tmp / "m6_out.msh").exists():
        check("(5) V3(d) の品質: 作り直した M6 の .msh がある", False)
    else:
        import re
        cfg = cfg_src.read_text(encoding="utf-8")
        cfg = cfg.replace('discretization: "node"', 'discretization: "cell"')
        cfg = re.sub(r",\s*(axisCentroidShift|nodeWallDirichlet):\s*\d+", "", cfg)
        cfg = re.sub(r'(meshFileName|valueFileName):\s*"[^"]*"', r'\1: "nozzle_qc.h5"', cfg)
        if 'discretization: "cell"' not in cfg or cfg.count("nozzle_qc.h5") != 2:
            check("(5) V3(d) の品質: cell の config を作れる", False, cfg.splitlines()[0][:150])
        else:
            for nm, src in (("元の msh", msh6), ("作り直した .msh (forge の形)", tmp / "m6_out.msh")):
                d = tmp / ("qc_" + ("orig" if src == msh6 else "rebuilt"))
                d.mkdir()
                shutil.copyfile(src, d / "nozzle.msh")
                (d / "solverConfig.yaml").write_text(cfg, encoding="utf-8")
                shutil.copyfile(fx / "bcondConfig.yaml", d / "bcondConfig.yaml")
                rc = subprocess.run([str(conv), "nozzle.msh", "nozzle_qc.h5"], cwd=d, capture_output=True, text=True)
                if rc.returncode != 0 or not (d / "nozzle_qc.h5").exists():
                    check(f"(5) V3(d) {nm}: 変換器が通る", False, (rc.stdout + rc.stderr)[-300:])
                    continue
                q = subprocess.run([sys.executable, str(QC_TOOL), "nozzle_qc.h5", "--ar-max", "5000"], cwd=d,
                                   capture_output=True, text=True)
                lines = [ln for ln in q.stdout.splitlines() if ln.startswith(("aspect", "skew", "VERDICT"))]
                check(f"(5) V3(d) {nm}: cell 変換 → check_mesh_quality --ar-max 5000 が PASS",
                      q.returncode == 0 and "VERDICT: PASS" in q.stdout, " / ".join(lines) or q.stderr[-300:])
                (d / "nozzle_qc.h5").unlink(missing_ok=True)

if os.environ.get("FORGE_KEEP_TMP"):
    print(f"作業フォルダを残した: {tmp}")
else:
    shutil.rmtree(tmp, ignore_errors=True)
print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
