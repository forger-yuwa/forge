"""plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #34: 入口の配管の半径を半分にした形 H を、明示の物理壁の上流の条件
(`geometry.physical_upstream`) でチェーンに入れたときの生成の経路の A/B (codex の判別 A/B、CFD 0 step、変換しない)。

- 修正後の C: Gref (`problem_d155_ns_prod_coldmesh_tw300_cgref.yaml`、core-grid の Gref の写し、ni 4719・nj 160・axis_cap_frac 0.015) に
  main の出口の終点の修正 (de35cb0b、壁の終点 +0.0698 r_t) を入れた設計。Gref の x の密度の表・第一層の表は終端が延びた壁を覆わない
  (mesh2d._check_table が例外) ので、両方の表の最後に [X_END_NEW, 最後の値] を足す (既存の節点はそのまま、終端の先は端の値で一定)。
  ni は Gref のまま。x の station は全域で正規化し直されるので旧と同じにはならない (codex の指摘どおり、受入れの条件にしない)。
  → `problem_d155_ns_prod_coldmesh_tw300_cgref_c.yaml`
- H (チェーンの経路 B の問題): C に `geometry.physical_upstream: {r_inlet: 3.2425, L_U: 9.0, L_pipe: 3.5}` を足し、第一層の表の x < 0 の
  節点の比を r_C(x)/r_H(x) 倍する (旧格子の第一層の実寸を保つ。節点で掛けて有効数字 6 桁に丸める; codex の Minor)。
  → `problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml`
- δ_r: 旧の δ_r の表 (prod_confirm/prep/delta_r_initial.csv、AWS) は旧の終端 95.2276 で終わり、延びた壁を覆わないので 1 本の B-spline の
  物理壁を作れない。生産の問題 (problem_d155_ns_prod.yaml、断熱) の積分法の初期壁 (prepare_ns と同じ `integral_delta_r`) を修正後の設計で
  作り直し、`delta_r_initial_C.csv` (prepare_ns と同じ列) に書いて、A・B・C のすべてがそれを `delta_r_from_table` で読む。
- A (参照の経路): C の物理壁 (チェーン) の x < 0 を、既往の H の定義 (forge-25 の cg_hi_prep.py の HalfInletWall) の 5 次 Hermite
  (−9, 3.2425, 0, 0) → x = 0 の C の (r, r′, r″) で直接置き換え、第一層の比を各 station で (C の比) × r_C/r_A にする。
- B (チェーンの経路): H の問題で `design_chain` → `build_physical_wall` → `mesh_params` → `generate_axisym_mesh` (prepare_ns と同じ関数)。
- 比べるもの: 壁の位置・1 階・2 階微分 (A − B ≤ 1e-8 r_t)、x ≥ 0 の壁と格子が C と同じ、第一層の距離 (A と B の差 ≤ 1 %)、品質
  (check_mesh_quality.metrics_vectorized を生成の座標に当てて AR ≤ 5000・skew ≤ 0.9)。旧 H (forge-25、出口の修正の前、
  prod_confirm/prep/wall_repr.json の壁から同じ定義で作る) との差・両方の表の終端・生成の ni・実際の区間の端も記録する。

使い方 (design の venv で): /home/sano/work/forge/design/.venv-opt/bin/python half_inlet_ab.py <出力ディレクトリ> [--pin <run_0062 のディレクトリ>]
  --pin: CFD ピンの凍結元 (geometry.initial_line_run)。無ければ問題 YAML の相対パス (この case の run_0062_euler_wallfit_fit_r1_ext6k)。
  照合は消さない (問題 YAML の写しの initial_line_run だけを書き換えて読む)。
出力: <出力>/half_inlet_ab.json、<出力>/delta_r_initial_C.csv、<出力>/wall_repr_H.json (B の壁ファイル)。問題 YAML 2 つはこの case に書く
  (既にあれば生成した中身と一致することを確かめ、違えば止める)。"""
import argparse, hashlib, json, math, sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design"))
sys.path.insert(0, str(HERE.parents[1] / "solver_density_cuda" / "tools"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.geometry.wall_axismach import check_required_attrs, load_wall_file, save_wall_file  # noqa: E402
from forge_design.geometry.wall_walldriven import _hermite_quintic, _poly_eval  # noqa: E402
from forge_design.meshing import mesh2d  # noqa: E402
from check_mesh_quality import metrics_vectorized  # noqa: E402

GREF = HERE / "problem_d155_ns_prod_coldmesh_tw300_cgref.yaml"        # core-grid の Gref の写し (sha256 85b13b3d…)
GREF_SHA = "85b13b3ddc75e25b4abf3e2865113906a6360f833308d2267f137fc37f90fa62"
PROD = HERE / "problem_d155_ns_prod.yaml"                              # δ_r の作り直し (断熱、Gref の壁と同じ規約)
C_YAML = HERE / "problem_d155_ns_prod_coldmesh_tw300_cgref_c.yaml"
H_YAML = HERE / "problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml"
OLD_WALL = HERE / "_band_ab" / "prod_confirm" / "prep"                # 出口の修正の前の生産の壁 (forge-25 の H の元の壁)
R_U_NEW, L_U_NEW, L_PIPE_NEW = 3.2425, 9.0, 3.5                         # H (core-grid plan §4.98)
X_U_NEW = -L_U_NEW
PU_LINE = f"  physical_upstream: {{r_inlet: {R_U_NEW}, L_U: {L_U_NEW}, L_pipe: {L_PIPE_NEW}}}"
NI, NJ = 4719, 160
TOL_WALL, TOL_FIRST = 1e-8, 0.01
FORGE25_JOIN = {"r": 1.0013889174618704, "r1": 0.0007743057045309643, "r2": 0.5011830155628192}   # half_inlet_wall.json (forge-25)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def table_line(key, tbl):
    return f"  {key}: [" + ", ".join(f"[{float(x)!r}, {float(v)!r}]" for x, v in tbl) + "]"


def replace_line(lines, prefix, new):
    idx = [i for i, ln in enumerate(lines) if ln.startswith(prefix)]
    if len(idx) != 1:
        raise SystemExit(f"{prefix!r} で始まる行が {len(idx)} 行 — 止める")
    old = lines[idx[0]]
    lines[idx[0]] = new
    return old


def load_pinned(txt: str, name: str, pin: Path, out: Path):
    """問題 YAML の写し (initial_line_run だけ凍結元の絶対パスに書き換え) を out に置いて読む。CFD ピンの照合は消さない。"""
    raw = yaml.safe_load(txt)
    raw["geometry"]["initial_line_run"] = str(pin)
    q = out / f"_pin_{name}"
    q.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    return RA.load_problem(q)


def write_or_check(path: Path, text: str):
    if path.exists():
        if path.read_text() != text:
            raise SystemExit(f"{path.name} が既にあり、生成した中身と違う — 止める (消してから回すか、差を確かめる)")
        return "unchanged"
    path.write_text(text)
    return "written"


class HalfInletWall:
    """参照の経路 A: 元の物理壁の x < 0 だけを既往の H の定義 (forge-25 の cg_hi_prep.py) で置き換えた壁。"""

    def __init__(self, base):
        self._b = base
        self.x_in, self.x_e = float(base.x_in), float(base.x_e)
        self.join = [float(np.asarray(base.r(np.array([0.0]), d))[0]) for d in (0, 1, 2)]
        self._L = -X_U_NEW
        self._c = _hermite_quintic(X_U_NEW, 0.0, R_U_NEW, 0.0, 0.0, *self.join)

    def r(self, x, deriv: int = 0):
        xa = np.atleast_1d(np.asarray(x, dtype=float))
        out = np.array(self._b.r(xa, deriv), dtype=float, copy=True)
        m = xa < 0.0
        if np.any(m):
            xi = (np.clip(xa[m], X_U_NEW, 0.0) - X_U_NEW) / self._L
            v = _poly_eval(self._c, xi, deriv) / self._L ** deriv
            out[m] = np.where(xa[m] < X_U_NEW, (R_U_NEW if deriv == 0 else 0.0), v)
        return out


class RestoredWall:
    """壁ファイルの物理壁 (load_wall_file) を mesh2d の壁として使う。"""

    def __init__(self, w):
        self._p = w["physical"]
        self.x_in, self.x_e = w["domain"]

    def r(self, x, deriv: int = 0):
        return self._p.r(np.asarray(x, dtype=float), deriv)


def mesh_with_first(wall, prm, first_fn=None):
    """generate_axisym_mesh。first_fn があれば第一層の比を差し替える (A の各 station の補正; forge-25 と同じ差し替え)。"""
    orig = mesh2d._first_frac_profile
    if first_fn is not None:
        mesh2d._first_frac_profile = lambda xs, p_: first_fn(np.asarray(xs, dtype=float), orig(xs, p_))
    try:
        c, q, b = mesh2d.generate_axisym_mesh(wall, prm)
    finally:
        mesh2d._first_frac_profile = orig
    return c, q, b


def quality(coords, quads):
    nq = len(quads)
    ars, sk = metrics_vectorized(coords[:, :2], quads.ravel(), np.arange(1, nq + 1) * 4, np.full(nq, 9), False)
    i = int(np.argmax(ars))
    return {"cells": nq, "ar_max": float(ars.max()), "ar_p99": float(np.percentile(ars, 99)), "skew_max": float(sk.max()),
            "ar_gt_5000": int((ars > 5000).sum()), "skew_gt_0.9": int((sk > 0.9).sum()),
            "pass": bool(ars.max() <= 5000.0 and sk.max() <= 0.9), "quad_at_ar_max": i}


def first_layer(coords, scale):
    P = coords[:, :2].reshape(NI, NJ, 2) / scale
    return np.hypot(*(P[:, NJ - 1] - P[:, NJ - 2]).T), P


def maxdiff(fa, fb, xs):
    out = {}
    for n in range(3):
        d = np.abs(np.asarray(fa(xs, n)) - np.asarray(fb(xs, n)))
        i = int(np.argmax(d))
        out[f"d{n}"] = {"max": float(d[i]), "x": float(xs[i])}
    return out


REGIONS = {"pipe_contraction_x_lt_0": (-np.inf, 0.0), "throat_0_1": (0.0, 1.0), "x_1_40": (1.0, 40.0), "test_40_94": (40.0, 94.0),
           "end_94_old_end": (94.0, np.inf)}


def region_maxdiff(fa, fb, xs):
    """区間ごとの |fa − fb| の最大 (位置・1 階・2 階)。"""
    out = {}
    for k, (lo, hi) in REGIONS.items():
        m = (xs >= lo) & (xs < hi)
        if m.any():
            out[k] = {f"d{n}": float(np.max(np.abs(np.asarray(fa(xs[m], n)) - np.asarray(fb(xs[m], n))))) for n in range(3)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--pin", default=None)
    a = ap.parse_args()
    out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    if sha(GREF) != GREF_SHA:
        raise SystemExit(f"{GREF.name} の sha256 が core-grid の Gref と違う — 止める")
    rec = {"plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §4.9・§5.1 #34", "tool": "half_inlet_ab.py",
           "environment": RA.environment_record(), "gref_sha256": GREF_SHA}

    # --- 1. 修正後の設計 (CFD ピンの照合つき) と δ_r の作り直し ---
    g_txt = GREF.read_text()
    g = yaml.safe_load(g_txt)
    pin = Path(a.pin).resolve() if a.pin else (HERE / g["geometry"]["initial_line_run"]).resolve()
    rec["cfd_pin_dir"] = str(pin)
    pp = load_pinned(PROD.read_text(), PROD.name, pin, out)
    d = RA.design_chain(pp)
    w_d = d["wall"]
    x_in, x_e = float(w_d.x_in), float(w_d.x_e)
    res, drx_fn, info = RA.integral_delta_r(pp, d, pp.raw["deltastar_initializer"])
    csv = out / "delta_r_initial_C.csv"
    np.savetxt(csv, np.c_[res["x"], res["delta_r"], res["dstar_n"], res["theta"] * float(pp.spec["r_throat"]), res["H"], res["N"],
                          res["Cf"], res["M"], res["Tw"], res["Taw"]],
               delimiter=",", comments="", header="x_rt,delta_r_rt,dstar_n_rt,theta_m,H,N,Cf,M_e,Tw,Taw")
    tbl_r = np.loadtxt(csv, delimiter=",", skiprows=1)
    drx = RA.delta_r_from_table(tbl_r[:, 0], tbl_r[:, 1])               # prepare_ns の delta_r_csv の読み方
    rec["design"] = {"n_wall_inv": int(len(d["wall_inv"])), "x_in": x_in, "x_e": x_e, "x_F": float(d["wall_inv"][-1, 0]),
                     "L_c": float(d["L_c"]), "cfd_pin": d["initial_line"]}
    rec["delta_r"] = {"source": "integral_delta_r (problem_d155_ns_prod.yaml、断熱、修正後の設計)", "csv": str(csv), "sha256": sha(csv),
                      "x_range": [float(tbl_r[0, 0]), float(tbl_r[-1, 0])], "throat": info["delta_r_throat"], "exit": info["delta_r_exit"]}

    # --- 2. 修正後の C の問題 YAML (表の終端を延ばす) ---
    x_end_old = float(g["mesh"]["x_density_table"][-1][0])
    if float(g["mesh"]["wall_first_frac_table"][-1][0]) != x_end_old:
        raise SystemExit("Gref の 2 つの表の終端が違う — 止める")
    x_end_new = math.ceil(x_e * 1e4) / 1e4
    lines = g_txt.splitlines()
    for k in ("x_density_table", "wall_first_frac_table"):
        if table_line(k, g["mesh"][k]) not in lines:
            raise SystemExit(f"Gref の {k} の行を同じ書式で書き直せない — 止める")
    xd_new = [list(r) for r in g["mesh"]["x_density_table"]] + [[x_end_new, float(g["mesh"]["x_density_table"][-1][1])]]
    ff_new = [list(r) for r in g["mesh"]["wall_first_frac_table"]] + [[x_end_new, float(g["mesh"]["wall_first_frac_table"][-1][1])]]
    replace_line(lines, "name: ", "name: problem_d155_ns_prod_coldmesh_tw300_cgref_c")
    replace_line(lines, "  x_density_table:", table_line("x_density_table", xd_new))
    replace_line(lines, "  wall_first_frac_table:", table_line("wall_first_frac_table", ff_new))
    note_end = (f"# 表の終端: Gref の {x_end_old!r} は出口の終点の修正で延びた壁の終点 {x_e!r} r_t を覆わないので、x_density_table・\n"
                f"# wall_first_frac_table の最後に [{x_end_new!r}, 最後の値] を足した (既存の節点はそのまま。ni は Gref の {NI} のまま、x の station は全域で正規化し直される)。\n")
    hdr_c = ("# plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #34: 修正後の今の形 C (Gref に main の出口の終点の修正 de35cb0b を入れた設計)。\n"
             "# half_inlet_ab.py が problem_d155_ns_prod_coldmesh_tw300_cgref.yaml (core-grid の Gref の写し) から作る (手で編集しない)。\n"
             "# Gref との差は name と 2 つの表の終端だけ。\n" + note_end)
    body = list(lines)
    c_txt = hdr_c + "\n".join(body) + "\n"
    c = yaml.safe_load(c_txt)
    g2, c2 = dict(g), dict(c)
    for k in ("name", "mesh"):
        g2.pop(k); c2.pop(k)
    gm, cm = dict(g["mesh"]), dict(c["mesh"])
    for k in ("x_density_table", "wall_first_frac_table"):
        gm.pop(k); cm.pop(k)
    if g2 != c2 or gm != cm:
        raise SystemExit("C の問題 YAML の差が name と 2 つの表以外にもある — 止める")
    rec["yaml_C"] = {"file": C_YAML.name, "status": write_or_check(C_YAML, c_txt)}

    scale = float(pp.spec["r_throat"])
    pc = load_pinned(c_txt, C_YAML.name, pin, out)
    dc = RA.design_chain(pc)
    WC = RA.build_physical_wall(pc, dc, scale, delta_r_x=drx, offset="radial")
    check_required_attrs(WC)
    if WC.validate():
        raise SystemExit(f"C の物理壁が不合格: {WC.validate()}")
    mpC = RA.mesh_params(pc, scale, ni=561, nj=97, wall_first_frac=4.5e-5)
    # Gref の表のまま (終端 x_end_old) では修正後の壁を覆わず mesh2d の表の検査で止まる (codex の読み取りの再現)
    pg = load_pinned(g_txt, GREF.name, pin, out)
    try:
        mesh2d.generate_axisym_mesh(WC, RA.mesh_params(pg, scale, ni=561, nj=97, wall_first_frac=4.5e-5))
        rec["gref_tables_on_fixed_wall"] = "例外にならなかった (想定外)"
    except ValueError as e:
        rec["gref_tables_on_fixed_wall"] = f"例外 (想定どおり): {e}"

    # --- 3. H の問題 YAML (physical_upstream + 第一層の表の x < 0 の節点を r_C/r_H 倍) ---
    h0 = yaml.safe_load(c_txt)
    h0["geometry"]["physical_upstream"] = {"r_inlet": R_U_NEW, "L_U": L_U_NEW, "L_pipe": L_PIPE_NEW}
    ph0 = load_pinned(yaml.safe_dump(h0, allow_unicode=True, sort_keys=False), "h_noratio.yaml", pin, out)
    WH0 = RA.build_physical_wall(ph0, RA.design_chain(ph0), scale, delta_r_x=drx, offset="radial")
    ff_h = []
    for x, v in ff_new:
        if x < 0.0:
            xq = np.array([max(float(x), x_in)])
            ratio = float(WC.r(xq)[0] / WH0.r(xq)[0])
            ff_h.append([float(x), float(f"{v * ratio:.6g}")])
        else:
            ff_h.append([float(x), float(v)])
    hl = list(body)
    replace_line(hl, "name: ", "name: problem_d155_ns_prod_coldmesh_tw300_cgref_h")
    replace_line(hl, "  wall_first_frac_table:", table_line("wall_first_frac_table", ff_h))
    i_pw = [i for i, ln in enumerate(hl) if ln.startswith("  pw_upstream: poly")]
    if len(i_pw) != 1:
        raise SystemExit("pw_upstream: poly の行が 1 行でない — 止める")
    hl.insert(i_pw[0] + 1, PU_LINE + "   # 物理壁だけの上流 (配管 r 3.2425・縮流部 x ∈ [−9, 0] の 5 次多項式)。r_inlet・L_U・L_pipe (上) は設計壁と CFD ピンの照合の値")
    hdr_h = ("# plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #34: 入口の配管の半径を半分にした形 H (core-grid plan §4.98 の定義) を\n"
             "# 明示の物理壁の上流の条件 geometry.physical_upstream で表した問題。half_inlet_ab.py が ..._cgref_c.yaml から作る (手で編集しない)。\n"
             "# C との差は name・geometry.physical_upstream・wall_first_frac_table の x < 0 の節点の比 (r_C(x)/r_H(x) 倍、有効数字 6 桁; 旧格子の第一層の実寸)。\n"
             "# r_C・r_H は同じ δ_r (half_inlet_ab.py の delta_r_initial_C.csv) の C・H の物理壁の節点での値。δ_r は delta_r_csv で渡す (積分法の初期壁とは併用しない)。\n"
             + note_end)
    h_txt = hdr_h + "\n".join(hl) + "\n"
    h = yaml.safe_load(h_txt)
    c3, h3 = dict(c), dict(h)
    for k in ("name", "mesh", "geometry"):
        c3.pop(k); h3.pop(k)
    cg, hg = dict(c["geometry"]), dict(h["geometry"])
    hg.pop("physical_upstream")
    cm, hm = dict(c["mesh"]), dict(h["mesh"])
    cm.pop("wall_first_frac_table"); hm.pop("wall_first_frac_table")
    if c3 != h3 or cg != hg or cm != hm:
        raise SystemExit("H の問題 YAML の差が name・physical_upstream・第一層の表以外にもある — 止める")
    rec["yaml_H"] = {"file": H_YAML.name, "status": write_or_check(H_YAML, h_txt)}

    # --- 4. B (チェーンの経路) ---
    ph = load_pinned(h_txt, H_YAML.name, pin, out)
    dh = RA.design_chain(ph)                      # CFD ピンの照合は geometry の r_inlet・L_U・L_pipe (凍結した設計) で通る
    WB = RA.build_physical_wall(ph, dh, scale, delta_r_x=drx, offset="radial")
    check_required_attrs(WB)
    vB = WB.validate()
    if vB:
        raise SystemExit(f"B の物理壁が不合格: {vB}")
    wpath, wrec, wsha = save_wall_file(out, WB, scale, "single_bspline")
    wpath.rename(out / "wall_repr_H.json")
    mpB = RA.mesh_params(ph, scale, ni=561, nj=97, wall_first_frac=4.5e-5)
    if not (mpB.ni == mpC.ni == NI and mpB.nj == mpC.nj == NJ and mpB.axis_cap_frac == mpC.axis_cap_frac == 0.015):
        raise SystemExit("格子の構成が Gref (ni 4719・nj 160・axis_cap_frac 0.015) でない — 止める")

    # --- 5. A (参照の経路) ---
    WA = HalfInletWall(WC)

    # --- 6. 格子 ---
    cC, qC, _ = mesh2d.generate_axisym_mesh(WC, mpC)
    # 比を先に作ってから掛ける (f·r_C/r_A の順に割ると x ≥ 0 で比 1 でも丸めで 1 ulp ずれる。forge-25 の _ff はその順で 1.4e-14 r_t ずれていた)
    cA, qA, _ = mesh_with_first(WA, mpC, lambda xs, f: f * (np.asarray(WC.r(xs)) / np.asarray(WA.r(xs))))
    cB, qB, _ = mesh2d.generate_axisym_mesh(WB, mpB)
    dC, PC = first_layer(cC, scale)
    dA, PA = first_layer(cA, scale)
    dB, PB = first_layer(cB, scale)
    xs_st = PC[:, NJ - 1, 0]
    same_st = bool(np.array_equal(PA[:, NJ - 1, 0], xs_st) and np.array_equal(PB[:, NJ - 1, 0], xs_st))
    dn = xs_st >= 0.0
    cols_dn = np.where(dn)[0]

    # --- 7. 比べる ---
    xs = np.unique(np.r_[np.linspace(x_in, x_e, 400001), xs_st, [X_U_NEW, 0.0], np.linspace(-1e-3, 1e-3, 2001)])
    xs = xs[(xs >= x_in) & (xs <= x_e)]
    ab = maxdiff(WA.r, WB.r, xs)
    xdn = xs[xs >= 0.0]
    ab_C_dn = {"A": maxdiff(WA.r, WC.r, xdn), "B": maxdiff(WB.r, WC.r, xdn)}
    relA_B = np.abs(dB / dA - 1.0)
    i_r = int(np.argmax(relA_B))
    rec["AB"] = {
        "wall_A_minus_B": ab,
        "wall_pass": bool(all(v["max"] <= TOL_WALL for v in ab.values())),
        "x_ge_0_vs_C": ab_C_dn,
        "x_ge_0_pass": bool(all(v["max"] <= TOL_WALL for s in ab_C_dn.values() for v in s.values())),
        "mesh_x_ge_0_cols_vs_C": {"A_bit_equal": bool(np.array_equal(cA.reshape(NI, NJ, 3)[cols_dn], cC.reshape(NI, NJ, 3)[cols_dn])),
                                  "B_bit_equal": bool(np.array_equal(cB.reshape(NI, NJ, 3)[cols_dn], cC.reshape(NI, NJ, 3)[cols_dn])),
                                  "B_max_abs_m": float(np.max(np.abs(cB.reshape(NI, NJ, 3)[cols_dn] - cC.reshape(NI, NJ, 3)[cols_dn]))),
                                  "n_cols": int(dn.sum())},
        "stations_same_A_B_C": same_st,
        "first_layer": {"B_over_A_max_rel": float(relA_B[i_r]), "x_at_max": float(xs_st[i_r]), "pass": bool(relA_B[i_r] <= TOL_FIRST),
                        "A_over_C_range": [float(np.min(dA / dC)), float(np.max(dA / dC))],
                        "B_over_C_range": [float(np.min(dB / dC)), float(np.max(dB / dC))]},
        "quality": {"C": quality(cC, qC), "A": quality(cA, qA), "B": quality(cB, qB)},
        "tol": {"wall_rt": TOL_WALL, "first_layer_rel": TOL_FIRST, "ar_max": 5000, "skew_max": 0.9},
    }
    rec["AB"]["pass"] = bool(rec["AB"]["wall_pass"] and rec["AB"]["x_ge_0_pass"] and rec["AB"]["first_layer"]["pass"] and same_st
                             and all(v["pass"] for v in rec["AB"]["quality"].values()))
    rec["H_wall"] = {"cfd_pin_match_H_problem": dh["initial_line"].get("match"), "physical_upstream": WB.physical_upstream, "throat": {"x": WB.x_throat, "r": WB.r_throat, "kappa": WB.kappa_throat},
                     "upstream_gate": {k: WB.upstream_gate.get(k) for k in ("pass", "max_abs_d2_change_vs_H", "max_seam_jump", "H_ref")},
                     "join_at_x0_A": WA.join, "wall_file": str(out / "wall_repr_H.json"), "wall_file_sha256": wsha,
                     "C_throat": {"x": WC.x_throat, "r": WC.r_throat, "kappa": WC.kappa_throat}}

    # --- 8. 旧 (forge-25 の H・Gref、出口の修正の前) との差 ---
    wold = load_wall_file(OLD_WALL)
    Cold = RestoredWall(wold)
    Hold = HalfInletWall(Cold)
    xo_e = float(wold["domain"][1])
    xcm = xs[xs <= xo_e]
    st_old = mesh2d._x_stations(float(wold["domain"][0]), xo_e, NI, mpC.throat_refine, mpC.throat_width, mpC.local_center,
                                mpC.local_refine, mpC.local_width, g["mesh"]["x_density_table"])
    dst = xs_st - st_old
    j = int(np.argmax(np.abs(dst)))
    dens_int = {}
    for lab, (x0_, x1_, t_) in {"old": (float(wold["domain"][0]), xo_e, g["mesh"]["x_density_table"]),
                                "new": (x_in, x_e, xd_new)}.items():
        xx = np.linspace(x0_, x1_, 8001)
        de = mesh2d._x_density(xx, mpC.throat_refine, mpC.throat_width, mpC.local_center, mpC.local_refine, mpC.local_width,
                               np.asarray(t_, dtype=float))
        dens_int[lab] = float(np.sum(0.5 * (de[1:] + de[:-1]) * np.diff(xx)))
    rec["vs_old"] = {
        "old_wall": str(OLD_WALL / "wall_repr.json"), "old_wall_sha256": sha(OLD_WALL / "wall_repr.json"),
        "domain_old": wold["domain"], "domain_new": [x_in, x_e], "x_e_shift_rt": x_e - xo_e, "x_e_shift_mm": (x_e - xo_e) * scale * 1e3,
        "table_end_old": x_end_old, "table_end_new": x_end_new, "ni_old": NI, "ni_new": NI,
        "ni_to_keep_mean_spacing": 1 + (NI - 1) * dens_int["new"] / dens_int["old"],
        "station_shift": {"max_abs_rt": float(abs(dst[j])), "x_at_max": float(xs_st[j]), "at_x0_col": float(dst[int(np.argmin(np.abs(xs_st)))]),
                          "first": float(dst[0]), "last": float(dst[-1])},
        "join_old_H_vs_forge25_record": {k: Hold.join[i] - FORGE25_JOIN[k] for i, k in enumerate(("r", "r1", "r2"))},
        "station_shift_by_region": {k: float(np.max(np.abs(dst[(xs_st >= lo) & (xs_st < hi)]))) for k, (lo, hi) in REGIONS.items()},
        "C_new_minus_C_old_by_region": region_maxdiff(WC.r, Cold.r, xcm),
        # 内訳: 設計壁 S (joint の当てはめは全点の同時の当てはめなので、終点の 1 点の追加で全域が少し動く) と δ_r (作り直し)
        "design_S_new_minus_old_by_region": region_maxdiff(w_d.r, wold["design"].r, xcm[xcm >= 0.0]),
        "delta_r_new_minus_old_by_region": region_maxdiff(lambda x, n: WC.r(x, n) - w_d.r(x, n),
                                                          lambda x, n: Cold.r(x, n) - wold["design"].r(x, n), xcm[xcm >= 0.0]),
        "H_new_B_minus_H_old_by_region": region_maxdiff(WB.r, Hold.r, xcm),
        "C_new_minus_C_old": maxdiff(WC.r, Cold.r, xcm),
        "H_new_A_minus_H_old": maxdiff(WA.r, Hold.r, xcm),
        "H_new_B_minus_H_old": maxdiff(WB.r, Hold.r, xcm),
        "join_new_minus_old": [WA.join[i] - Hold.join[i] for i in range(3)],
    }
    # 旧 δ_r の表の範囲 (旧の終端で切れた表) では延びた壁を 1 本の B-spline にできない (C を作り直す理由の確かめ)
    m_old = tbl_r[:, 0] <= xo_e
    try:
        RA.build_physical_wall(pc, dc, scale, delta_r_x=RA.delta_r_from_table(tbl_r[m_old, 0], tbl_r[m_old, 1]), offset="radial")
        rec["vs_old"]["old_table_end_check"] = "例外にならなかった (想定外)"
    except ValueError as e:
        rec["vs_old"]["old_table_end_check"] = f"例外 (想定どおり): {str(e)[:160]}"
    (out / "half_inlet_ab.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    print(json.dumps({"AB_pass": rec["AB"]["pass"], "wall_A_minus_B": ab, "first_layer": rec["AB"]["first_layer"],
                      "quality": {k: (v["ar_max"], v["skew_max"], v["pass"]) for k, v in rec["AB"]["quality"].items()},
                      "x_ge_0": rec["AB"]["mesh_x_ge_0_cols_vs_C"], "vs_old": {k: rec["vs_old"][k] for k in
                      ("x_e_shift_mm", "table_end_old", "table_end_new", "ni_to_keep_mean_spacing", "station_shift", "join_new_minus_old")}},
                     ensure_ascii=False, indent=1, default=str))
    if not rec["AB"]["pass"]:
        sys.exit("A/B 不合格")


if __name__ == "__main__":
    main()
