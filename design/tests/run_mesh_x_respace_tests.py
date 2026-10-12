"""区間の局所の並べ直し `mesh.x_respace` の試験 (plan tooling-nozzle-isothermal-wall-chain §4.8-1・§5.1 #29, 2026-10-13)。
(a) 旧方式との一致: forge-25 の G0c の作り方 (commit 96197821 の case/45.isobutane_m6_d155/core_grid_mesh.py:451-473
    `grefc_stations`、下の _ref_grefc_stations にそのまま写した) と、station・生成座標・接続がビット一致する
    (case/45 の冷却壁の mesh ブロック、壁は解析的な仮の壁で domain だけ実物に合わせる)。区間の端は −7.9904141995・−2.0011935537、
    ni 4467 (codex diagnose 2026-10-13 の再計算)。端の station から blend を取る別の解釈では 4468 列になる (採らないことの確認)。
(a') case/45 の物理壁 (_band_ab/prod_confirm/prep、git 管理外) があれば、problem_d155_ns_prod_coldmesh_tw300_g0c.yaml から
    _mesh_params_from → generate_axisym_mesh → write_msh41_2d (17 桁) で作った .msh の sha256 が make_g0c.py 96197821 の出力
    (285f0fa4ecc0ed1c…) と一致する。無ければ SKIP と表示する (合格に数えない)。
(b) x_respace の無い問題は変更前と同じ: 変更前の _x_stations (commit 09d1d39f の写し) と station がビット一致 (ガウスの細分・局所の細分・表)。
(c) 区間の外の列は元の配列 (ni_ref 列) とビット一致し、区間の中の列だけが減る。
(d) 不正な指定・ni の不一致・複数の区間は ValueError で止まる。Euler の mesh_euler には書けない。
(e) _mesh_params_from / mesh_params が x_respace を Mesh2DParams に渡す (キーが無ければ None)。
2026-10-13 の結果は plan §5.1 #29 の行を参照。"""
import hashlib
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.meshing import mesh2d  # noqa: E402
from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
C45 = REPO / "case/45.isobutane_m6_d155"
G0_YAML = C45 / "problem_d155_ns_prod_coldmesh_tw300.yaml"
G0C_YAML = C45 / "problem_d155_ns_prod_coldmesh_tw300_g0c.yaml"
WALL_REF = C45 / "_band_ab" / "prod_confirm" / "prep"
G0C_MSH_SHA256 = "285f0fa4ecc0ed1c04cb771d2da6d437d18509b0699326f78670f29a137ded6b"   # make_g0c.py 96197821 (2026-10-13 に作り直して照合)
fails = 0


def check(name, ok, info=""):
    global fails
    print(("ok   " if ok else "FAIL ") + name + (f" ({info})" if info else ""))
    fails += (not ok)


def raises(fn):
    try:
        fn()
    except Exception as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"
    return None


def _old_x_stations(x0, x1, ni, refine, width, local_center=0.0, local_refine=1.0, local_width=0.75, density_table=None):
    """変更前 (commit 09d1d39f) の mesh2d._x_stations の写し。"""
    xs = np.linspace(x0, x1, 8001)
    dens = 1.0 + (refine - 1.0) * np.exp(-((xs / width) ** 2))
    if local_refine > 1.0:
        dens = dens * (1.0 + (local_refine - 1.0)
                       * np.exp(-(((xs - local_center) / local_width) ** 2)))
    if density_table is not None:
        t = mesh2d._check_table("x_density_table", density_table, x0, x1)
        dens = dens * np.interp(xs, t[:, 0], t[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xs))])
    cum /= cum[-1]
    return np.interp(np.linspace(0.0, 1.0, ni), cum, xs)


GREFC_SEG = (-8.0, -7.5, -2.5, -2.0)


def _ref_grefc_stations(x0, x1, refine, width, local_center, local_refine, local_width, table, ni_ref=4719, seg=GREFC_SEG, alt_ends=False):
    """forge-25 (commit 96197821) の core_grid_mesh.py:455-468 `grefc_stations` の写し (_ORIG_X_STATIONS = 変更前の _x_stations)。
    alt_ends=True は端の station から blend を取る別の解釈 (codex が 4468 列になると示したもの; 採らないことの確認用)。"""
    xs = _old_x_stations(x0, x1, ni_ref, refine, width, local_center, local_refine, local_width, table)
    t = np.asarray(table, float); a0, a1, b1, b0 = seg
    ia = int(np.argmin(np.abs(xs - a0))); ib = int(np.argmin(np.abs(xs - b0))); xa, xb = float(xs[ia]), float(xs[ib])
    if alt_ends:
        a1, b1 = xa + (seg[1] - seg[0]), xb - (seg[3] - seg[2])
    xf = np.linspace(xa, xb, 20001); d = np.interp(xf, t[:, 0], t[:, 1])
    w = np.clip((xf - xa) / (a1 - xa), 0, 1) * np.clip((xb - xf) / (xb - b1), 0, 1)
    w = w * w * (3 - 2 * w)
    dn = d * (1 - w * (1 - np.maximum(1.0, d / 2) / d))
    cr = np.concatenate([[0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(xf))]); cn = np.concatenate([[0], np.cumsum(0.5 * (dn[1:] + dn[:-1]) * np.diff(xf))])
    nseg = int(round((ib - ia) * cn[-1] / cr[-1]))
    s = np.interp(np.linspace(0, cn[-1], nseg + 1), cn, xf); s[0], s[-1] = xa, xb
    return np.concatenate([xs[:ia], s, xs[ib + 1:]]), (ia, ib, nseg)


class Wall:
    """解析的な仮の壁 (r_t 単位)。domain は case/45 の物理壁 (wall_repr.json) と同じ。"""
    x_in, x_e = -12.5, 95.22757793472597

    def r(self, x, d=0):
        x = np.asarray(x, dtype=float)
        if d == 1:
            return np.where(x < 0, 0.08 * x, 0.09 - 0.0006 * x)
        return np.where(x < 0, 1.0 + 0.04 * x * x, 1.0 + 0.09 * x - 0.0003 * x * x)


W = Wall()
XR = [{"x0": -8.0, "x1": -2.0, "blend": 0.5, "factor": 2.0, "ni_ref": 4719}]
m_g0 = yaml.safe_load(G0_YAML.read_text())["mesh"]
m_g0c = yaml.safe_load(G0C_YAML.read_text())["mesh"]


def same_mesh(a, b):
    ca, qa, ba = a
    cb, qb, bb = b
    return (ca.shape == cb.shape and np.array_equal(ca, cb) and np.array_equal(qa, qb)
            and all(np.array_equal(np.asarray(ba[k]), np.asarray(bb[k])) for k in ("inlet", "outlet", "wall", "axis")))


# --- (a) 旧方式との一致 ---------------------------------------------------------------------------------------------------
check("(a) G0c の YAML は x_respace・ni 4467・ar_max 10000 で、ほかの mesh のキーは G0 と同じ",
      m_g0c["x_respace"] == XR and m_g0c["ni"] == 4467 and m_g0c["ar_max"] == 10000
      and {k: v for k, v in m_g0c.items() if k not in ("ni", "x_respace", "ar_max")} == {k: v for k, v in m_g0.items() if k not in ("ni", "ar_max")})
prm_c = RA._mesh_params_from(m_g0c, 1.0, 561, 97, 4.5e-5)
xs_new, info = mesh2d._x_stations_respaced(W.x_in, W.x_e, prm_c)
xs_ref, (ia_r, ib_r, nseg_r) = _ref_grefc_stations(W.x_in, W.x_e, prm_c.throat_refine, prm_c.throat_width, prm_c.local_center,
                                                   prm_c.local_refine, prm_c.local_width, prm_c.x_density_table)
check("(a) station が forge-25 の grefc_stations とビット一致", xs_new.shape == xs_ref.shape and np.array_equal(xs_new, xs_ref),
      (len(xs_new), len(xs_ref)))
check("(a) 列の数 4467・区間の列 589 → 337・区間の端 (−7.9904141995, −2.0011935537)",
      len(xs_new) == 4467 and (info["ia"], info["ib"], info["n_new"]) == (ia_r, ib_r, nseg_r) and info["n_old"] == 589 and info["n_new"] == 337
      and abs(info["xa"] + 7.9904141995) < 1e-10 and abs(info["xb"] + 2.0011935537) < 1e-10, info)
xs_alt, _ = _ref_grefc_stations(W.x_in, W.x_e, prm_c.throat_refine, prm_c.throat_width, prm_c.local_center,
                                prm_c.local_refine, prm_c.local_width, prm_c.x_density_table, alt_ends=True)
check("(a) 端の station から blend を取る別の解釈は 4468 列 (採らない; ブレンドの内側の端は指定の x0 + blend・x1 − blend)",
      len(xs_alt) == 4468, len(xs_alt))
_orig = mesh2d._x_stations
try:
    mesh2d._x_stations = lambda *a, **k: xs_ref                    # forge-25 の make_g0c.py と同じ差し替えで作る基準
    prm_ref = RA._mesh_params_from(dict(m_g0, ni=4467), 1.0, 561, 97, 4.5e-5)
    mesh_ref = generate_axisym_mesh(W, prm_ref)
finally:
    mesh2d._x_stations = _orig
mesh_new = generate_axisym_mesh(W, prm_c)
check("(a) 生成座標・接続・境界辺が差し替え方式 (make_g0c.py) とビット一致 (壁法線の近壁層・第一層の表を含む)", same_mesh(mesh_new, mesh_ref))

# --- (a') case/45 の物理壁での .msh のバイト一致 (壁のファイルが無ければ SKIP) ------------------------------------------------
if (WALL_REF / "wall_repr.json").exists():
    from forge_design.geometry.wall_axismach import load_wall_file

    class PhysWall:
        """make_g0c.py が使う forge-25 core_grid_mesh.Wall と同じ (物理壁・domain・scale_m)。"""
        def __init__(self, ref):
            Wf = load_wall_file(ref)
            self.ph = Wf["physical"]; self.x_in, self.x_e = (float(v) for v in Wf["domain"]); self.scale = float(Wf["scale_m"])

        def r(self, x, d=0):
            return self.ph.r(np.asarray(x), d) if d else self.ph.r(np.asarray(x))

    pw = PhysWall(WALL_REF)
    prm_p = RA._mesh_params_from(m_g0c, pw.scale, int(m_g0c["ni"]), int(m_g0c["nj"]), float(m_g0c["wall_first_frac"]))
    c, q, b = generate_axisym_mesh(pw, prm_p)
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "nozzle.msh"
        write_msh41_2d(str(f), c, q, b, digits=int(m_g0c["msh_digits"]))
        h = hashlib.sha256(f.read_bytes()).hexdigest()
    check("(a') case/45 の物理壁: G0c の YAML から作った .msh が make_g0c.py 96197821 の出力とバイト一致", h == G0C_MSH_SHA256, h[:16])
else:
    print(f"SKIP (a') case/45 の物理壁 {WALL_REF} が無い (git 管理外) — .msh の sha256 の照合をしていない (合格に数えない)")

# --- (b) x_respace の無い問題は変更前と同じ ------------------------------------------------------------------------------
cases = [(-12.0, 95.0, 2000, 4.0, 3.0, 0.0, 1.0, 0.75, None), (-12.0, 95.0, 561, 3.0, 1.5, 1.2, 3.0, 0.5, None),
         (W.x_in, W.x_e, 4719, 1.0, 1.0, 0.0, 1.0, 0.75, m_g0["x_density_table"]),
         (W.x_in, W.x_e, 1500, 2.5, 2.0, 0.0, 1.0, 0.75, m_g0["x_density_table"])]
check("(b) _x_stations が変更前とビット一致 (ガウス・局所の細分・表)",
      all(np.array_equal(mesh2d._x_stations(*a), _old_x_stations(*a)) for a in cases))
prm_0 = RA._mesh_params_from(m_g0, 1.0, 561, 97, 4.5e-5)
check("(b) x_respace の無い mesh ブロックは Mesh2DParams.x_respace = None", prm_0.x_respace is None and Mesh2DParams().x_respace is None)
try:
    mesh2d._x_stations = _old_x_stations
    mesh_old = generate_axisym_mesh(W, prm_0)
finally:
    mesh2d._x_stations = _orig
mesh_0 = generate_axisym_mesh(W, prm_0)
check("(b) x_respace の無い G0 の生成座標・接続が変更前の _x_stations で作ったものとビット一致", same_mesh(mesh_0, mesh_old))

# --- (c) 区間の外の列はビット一致 -----------------------------------------------------------------------------------------
ni0, nj = prm_0.ni, mesh_0[0].shape[0] // prm_0.ni
C0 = mesh_0[0].reshape(ni0, nj, 3); C1 = mesh_new[0].reshape(len(xs_new), nj, 3)
ia, ib, nn = info["ia"], info["ib"], info["n_new"]
check("(c) 区間の前の列 (ia 列) と後ろの列が G0 とビット一致",
      np.array_equal(C0[:ia], C1[:ia]) and np.array_equal(C0[ib + 1:], C1[ia + nn + 1:]), (ia, ni0 - ib - 1))
check("(c) 区間の端の列も G0 とビット一致 (端は元の列に固定)", np.array_equal(C0[ia], C1[ia]) and np.array_equal(C0[ib], C1[ia + nn]))
dx0 = np.diff(xs_ref[ia:ia + nn + 1]); dxg = np.diff(_old_x_stations(W.x_in, W.x_e, 4719, 1.0, 1.0, 0.0, 1.0, 0.75, m_g0["x_density_table"])[ia:ib + 1])
check("(c) 区間の中は粗くなる (最大の間隔の比 ≤ factor)", 1.0 < dx0.max() / dxg.max() <= 2.0 + 1e-12, float(dx0.max() / dxg.max()))

# --- (d) 止まる --------------------------------------------------------------------------------------------------------------
def gen(xr, ni=4467, **kw):
    return generate_axisym_mesh(W, RA._mesh_params_from(dict(m_g0, ni=ni, x_respace=xr, **kw), 1.0, 561, 97, 4.5e-5))


e = raises(lambda: gen(XR, ni=4468))
check("(d) ni が並べ直した列の数と違えば止まる (列の数を示す)", e is not None and e.startswith("ValueError") and "4467" in e, (e or "")[:80])
two = [dict(XR[0]), {"x0": 10.0, "x1": 20.0, "blend": 0.5, "factor": 2.0, "ni_ref": 4719}]
for bad, why in ((two, "2 区間"), ([], "空の配列"), (XR[0], "辞書 (配列でない)"),
                 ([{k: v for k, v in XR[0].items() if k != "blend"}], "キーの不足"), ([dict(XR[0], facter=2.0)], "未知のキー"),
                 ([dict(XR[0], x0=-1.0)], "x0 ≥ x1"), ([dict(XR[0], x0=-20.0)], "x0 が x_in の外"), ([dict(XR[0], blend=3.5)], "2·blend > 幅"),
                 ([dict(XR[0], blend=0.0)], "blend ≤ 0"), ([dict(XR[0], factor=-1.0)], "factor ≤ 0"), ([dict(XR[0], ni_ref=4719.5)], "ni_ref が整数でない"),
                 ([dict(XR[0], ni_ref=30)], "ni_ref が粗すぎて端の列がブレンドの内側 (ni の照合の前に)")):
    e = raises(lambda bad=bad: gen(bad))
    check(f"(d) {why}は止まる", e is not None and e.startswith("ValueError")
          and ("区間の端の列" in e if "ni の照合の前" in why else True), (e or "")[:90])


class PE:
    def __init__(self, mesh_euler):
        self.mesh = {}; self.raw = {"mesh": {}, "mesh_euler": mesh_euler}; self.mesh_euler = mesh_euler


e = raises(lambda: RA.mesh_params_euler(PE({"ni": 2000, "x_respace": XR}), 1.0))
check("(d) Euler の mesh_euler には書けない (MESH_EULER_KEYS に無い)", e is not None and "x_respace" in e and "x_respace" not in RA.MESH_EULER_KEYS,
      (e or "")[:60])

# --- (e) 問題からの読み込み ---------------------------------------------------------------------------------------------------
p = load_problem(G0C_YAML)
mp = RA.mesh_params(p, 0.0766, 561, 97, 4.5e-5)
check("(e) mesh_params (prepare_ns の経路) が x_respace・ni 4467 を Mesh2DParams に渡す", mp.x_respace == XR and mp.ni == 4467, (mp.ni, mp.x_respace))
check("(e) G0 の問題は x_respace = None", RA.mesh_params(load_problem(G0_YAML), 0.0766, 561, 97, 4.5e-5).x_respace is None)

print(f"FAIL 件数: {fails}"); sys.exit(1 if fails else 0)
