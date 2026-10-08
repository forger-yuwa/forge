"""Euler 専用の格子の設定 `mesh_euler` の試験 (plan verification-case45-euler-total-enthalpy §5.1 #5a・§6 E3, 2026-10-07)。
(a) mesh だけで mesh_euler の無い問題は、Euler の prepare が run dir を作る前に移行先を示して止まる。
(b) mesh_euler の既定は全断面で wall_first_frac 0.005 の等比 (スロートの別指定なし・軸側の cap なし、ni・nj は従来の Euler の 321 × 65)。
(c) mesh_euler は mesh から補完しない (mesh に書いた壁際の細分・cap は Euler に効かない)。未知のキー・node でない discretization・
    辞書でない値は止まる。他の type に mesh_euler を書くと load_problem が止まる。
(d) NS (mesh_params) は mesh だけを読み、変更前の関数 (下の _old_mesh_params の写し) と Mesh2DParams も生成座標もビット一致する
    (case/45 の NS の問題すべて)。
(e) case/45 の Euler の問題: 全域 0.005 へ移した問題は 2000 × 97・0.005・スロートの別指定なし・cap なし。記録の格子を明示した問題
    (G1・旧 1100 × 65) は、変更前の Euler 経路 (mesh を読む) と Mesh2DParams・生成座標がビット一致する。E2 の腕 B の問題は、
    mesh_euler (標準の書き方) と変更前の mesh (両 0.005 + ブレンドの鍵) で生成座標がビット一致する (E3 の合格条件 (1) の生成器側)。
生成座標は解析的な仮の壁 (x_in・x_e・r(x)) で作る (座標は Mesh2DParams と壁の r(x) だけで決まる)。"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
C45 = REPO / "case/45.isobutane_m6_d155"
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


def _old_mesh_params(p, scale, ni, nj, wall_first_frac):
    """変更前 (commit 99431498) の runner_axismach.mesh_params の写し (Euler も NS も p.mesh を読んでいた)。"""
    m = p.mesh
    opt = lambda k: None if m.get(k) is None else float(m[k])  # noqa: E731
    return Mesh2DParams(ni=int(m.get("ni", ni)), nj=int(m.get("nj", nj)),
                        wall_first_frac=float(m.get("wall_first_frac", wall_first_frac)),
                        throat_refine=float(m.get("throat_refine", 3.0)),
                        throat_width=float(m.get("throat_width", 1.5)),
                        wall_first_frac_throat=opt("wall_first_frac_throat"),
                        wall_first_blend_x0=float(m.get("wall_first_blend_x0", 0.5)),
                        wall_first_blend_x1=float(m.get("wall_first_blend_x1", 6.0)),
                        wall_first_up_x0=opt("wall_first_up_x0"), wall_first_up_x1=opt("wall_first_up_x1"),
                        axis_gap_frac=opt("axis_gap_frac"), axis_cap_frac=opt("axis_cap_frac"), scale=scale)


class Wall:
    """解析的な仮の壁 (r_t 単位): 入口 x_in = −12、出口 x_e = 95、スロート x = 0 で r = 1。"""
    x_in, x_e = -12.0, 95.0

    def r(self, x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0, 1.0 + 0.04 * x * x, 1.0 + 0.09 * x - 0.0003 * x * x)


W = Wall()


def coords_of(mp):
    c, q, b = generate_axisym_mesh(W, mp)
    return c, q, b


def same_mesh(a, b):
    ca, qa, ba = a
    cb, qb, bb = b
    return (ca.shape == cb.shape and np.array_equal(ca, cb) and np.array_equal(qa, qb)
            and all(np.array_equal(np.asarray(ba[k]), np.asarray(bb[k])) for k in ("inlet", "outlet", "wall", "axis")))


class P:
    """mesh_params / mesh_euler_block が見る属性だけの仮の問題。"""
    def __init__(self, mesh, mesh_euler=None, has_euler=True):
        self.mesh = mesh
        self.raw = {"mesh": mesh}
        if has_euler:
            self.raw["mesh_euler"] = mesh_euler
        self.mesh_euler = mesh_euler


TMP = Path(tempfile.mkdtemp(prefix="mesh_euler_tests_"))


def write_yaml(src: Path, name: str, edit) -> Path:
    doc = yaml.safe_load(src.read_text())
    edit(doc)
    out = TMP / name
    out.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False))
    return out


# --- (a) mesh だけの Euler の問題は止まる -----------------------------------------------------------------------------
y_old = write_yaml(C45 / "problem_d155_euler_pin_G1_recal_mono.yaml", "mesh_only.yaml", lambda d: d.pop("mesh_euler"))
rd = TMP / "never_euler"
e = raises(lambda: RA.prepare(y_old, rd))
check("(a) mesh だけの問題で Euler の prepare が止まる (run dir を作る前)", e is not None and "mesh_euler" in e and not rd.exists(), (e or "")[:60])
check("(a) 止まるときに移行先 (新しい既定と旧格子の写し) を示す",
      e is not None and "wall_first_frac: 0.005" in e and "mesh の全キーを mesh_euler に写す" in e and "wall_first_frac_throat: 4.5e-06" in e)
e2 = raises(lambda: RA.mesh_euler_block(P({"ni": 321}, has_euler=False)))
check("(a) mesh_euler_block: mesh_euler の無い問題は ValueError", e2 is not None and e2.startswith("ValueError"))

# --- (b) 既定 ------------------------------------------------------------------------------------------------------------
md = RA.mesh_params_euler(P({}, {}), 0.0768075)
want = Mesh2DParams(ni=321, nj=65, wall_first_frac=0.005, throat_refine=3.0, throat_width=1.5, scale=0.0768075)
check("(b) mesh_euler: {} の既定 = 321 × 65・0.005・スロートの別指定なし・cap なし", md == want and md.wall_first_frac_throat is None
      and md.axis_cap_frac is None and md.axis_gap_frac is None, md)
mp = RA.mesh_params_euler(P({}, {"ni": 400, "nj": 97}), 1.0)
c, _, _ = coords_of(mp)
R = c[:, 1].reshape(mp.ni, mp.nj)
gw = (R[:, -1] - R[:, -2]) / R[:, -1]
ga = R[:, 1] / R[:, -1]
check("(b) 全断面で壁の第 1 間隔の比が 0.005 (等比、断面によらない)", np.allclose(gw, 0.005, rtol=1e-12, atol=0) and np.ptp(ga) < 1e-15,
      (float(gw.min()), float(gw.max()), float(ga.min())))

# --- (c) mesh から補完しない・不正な指定は止まる ----------------------------------------------------------------------
g1 = {"ni": 2000, "nj": 97, "wall_first_frac": 1.3e-5, "wall_first_frac_throat": 4.5e-6, "wall_first_up_x0": -9.0, "wall_first_up_x1": -4.0,
      "wall_first_blend_x0": 1.0, "wall_first_blend_x1": 17.0, "throat_refine": 4.0, "throat_width": 3.0, "axis_cap_frac": 0.02, "ar_max": 5000}
mx = RA.mesh_params_euler(P(g1, {"ni": 2000, "nj": 97}), 1.0)
check("(c) mesh の壁際の細分・cap・軸方向の鍵は mesh_euler に継承されない",
      (mx.wall_first_frac, mx.wall_first_frac_throat, mx.axis_cap_frac, mx.throat_refine, mx.throat_width) == (0.005, None, None, 3.0, 1.5), mx)
for bad, why in (({"ni": 2000, "wall_frist_frac": 0.01}, "未知のキー"), ({"discretization": "cell"}, "node でない discretization"),
                 ([2000, 97], "辞書でない")):
    e = raises(lambda bad=bad: RA.mesh_params_euler(P({}, bad), 1.0))
    check(f"(c) mesh_euler の{why}は止まる", e is not None and e.startswith("ValueError"), (e or "")[:70])
check("(c) discretization: node は受ける", RA.mesh_params_euler(P({}, {"discretization": "node"}), 1.0).ni == 321)
y_wt = write_yaml(REPO / "case/41.wind_tunnel_design/problem_m4_axismach.yaml", "other_type.yaml",
                  lambda d: (d.__setitem__("type", "wind_tunnel_axisym"), d.__setitem__("mesh_euler", {"ni": 321})))
e = raises(lambda: load_problem(y_wt))
check("(c) 他の type に mesh_euler を書くと load_problem が止まる", e is not None and "mesh_euler" in e, (e or "")[:60])
y_nd = write_yaml(C45 / "problem_d155_euler_pin_G1_recal_mono.yaml", "nondict.yaml", lambda d: d.__setitem__("mesh_euler", [2000, 97]))
e = raises(lambda: load_problem(y_nd))
check("(c) mesh_euler が辞書でないと load_problem が止まる", e is not None and "mesh_euler" in e, (e or "")[:60])
pp = load_problem(C45 / "problem_d155_euler_pin_G1_recal_mono.yaml")
check("(c) Problem.mesh_euler に YAML の mesh_euler が入る", pp.mesh_euler == pp.raw["mesh_euler"] and pp.mesh_euler["ni"] == 2000)

# --- (d) NS は mesh だけを読み、変更前と同じ ----------------------------------------------------------------------------
ns_files = sorted(f for f in C45.glob("problem_d155*.yaml") if "_ns" in f.name)
# 2026-10-08 に足した opt-in の格子のキー (表の第一層・x 密度・壁法線の近壁層) を使う問題は「変更前と同一」の対象外 (変更前の関数はこれらを知らない)。
# 代わりに、そのキーが Mesh2DParams に入ることを別に確かめる (plan tooling-nozzle-isothermal-wall-chain §5.1 #13・#15)
NEW_MESH_KEYS = ("wall_first_frac_table", "x_density_table", "wall_normal_layer")
new_key_files = [f for f in ns_files if any(k in (load_problem(f).mesh or {}) for k in NEW_MESH_KEYS)]
ns_files = [f for f in ns_files if f not in new_key_files]
for f in new_key_files:
    p = load_problem(f); mp = RA.mesh_params(p, 0.0768075, 561, 97, 4.5e-5)
    check(f"(d) 新しい格子のキーが Mesh2DParams に入る ({f.name})",
          all((getattr(mp, k) is not None) == (k in p.mesh) for k in NEW_MESH_KEYS))
nbad = []
for f in ns_files:
    p = load_problem(f)
    new = RA.mesh_params(p, 0.0768075, 561, 97, 4.5e-5)
    old = _old_mesh_params(p, 0.0768075, 561, 97, 4.5e-5)
    if new != old:
        nbad.append(f.name)
check(f"(d) NS の Mesh2DParams が変更前と同一 (case/45 の NS の問題 {len(ns_files)} 本)", not nbad and len(ns_files) >= 20, nbad[:3])
for f in ("problem_d155_ns_finemesh_recal_final_mono.yaml", "problem_d155_ns_recal_final_cap020.yaml", "problem_d155_ns_rt77p02_axgap025.yaml",
          "problem_d155_ns.yaml"):
    p = load_problem(C45 / f)
    check(f"(d) NS の生成座標・接続が変更前とビット一致: {f}",
          same_mesh(coords_of(RA.mesh_params(p, 0.0768075, 561, 97, 4.5e-5)), coords_of(_old_mesh_params(p, 0.0768075, 561, 97, 4.5e-5))))
pb = P(g1, {"ni": 2000, "nj": 97})
check("(d) 両ブロックがあっても NS は mesh だけを読む", RA.mesh_params(pb, 1.0, 561, 97, 4.5e-5) == _old_mesh_params(pb, 1.0, 561, 97, 4.5e-5))

# --- (e) case/45 の Euler の問題の移行 --------------------------------------------------------------------------------------
STD = ("problem_d155_euler_pin_G1_recal.yaml", "problem_d155_euler_pin_G1_recal_mono.yaml",
       "problem_d155_euler_pin_G1_recal_mono_moc.yaml", "problem_d155_euler_t0cluster_u5em3.yaml",
       "problem_d155_euler_e4_recal_d0.yaml")
STD = STD + tuple(sorted(f.name for f in C45.glob("problem_d155_euler_e4_recal_d[1-9].yaml")))   # E4 の段 2 (make-d1 が作る)
KEEP = ("problem_d155_euler_t0cluster_g1.yaml", "problem_d155_euler_pin_G1.yaml", "problem_d155_euler_pin_G0.yaml",
        "problem_d155_euler_c2final_n2400.yaml", "problem_d155_euler_c2final_n2400_pincal.yaml", "problem_d155_R2_Lc39.3_tp.yaml",
        "problem_d155_R2_Lc39.3_tp_Lpipe10.yaml", "problem_d155_R2_Lc39.3_tp_rt77p02.yaml", "problem_d155_cpg_euler.yaml",
        "problem_d155_trim_tp.yaml")
for f in STD:
    p = load_problem(C45 / f)
    m = RA.mesh_params_euler(p, 0.0768075)
    ok = ((m.ni, m.nj, m.wall_first_frac, m.wall_first_frac_throat, m.axis_cap_frac, m.axis_gap_frac, m.throat_refine, m.throat_width)
          == (2000, 97, 0.005, None, None, None, 4.0, 3.0)) and (m.throat_refine, m.throat_width) == (float(p.mesh["throat_refine"]), float(p.mesh["throat_width"]))
    check(f"(e) 全域 0.005 へ移した: {f} (軸方向の鍵は mesh と同じ)", ok, m)
for f in KEEP:
    p = load_problem(C45 / f)
    m = RA.mesh_params_euler(p, 0.0768075)
    o = _old_mesh_params(p, 0.0768075, 321, 65, 5.0e-3)      # 変更前の Euler 経路 (mesh を読む)
    check(f"(e) 記録の格子を明示: {f} の Mesh2DParams が変更前の Euler 経路と同一", m == o, (m.ni, m.nj, m.wall_first_frac))
check("(e) 記録の格子を明示した問題の生成座標が変更前とビット一致 (G1: t0cluster_g1)",
      same_mesh(coords_of(RA.mesh_params_euler(load_problem(C45 / KEEP[0]), 0.0768075)),
                coords_of(_old_mesh_params(load_problem(C45 / KEEP[0]), 0.0768075, 321, 65, 5.0e-3))))
pB = load_problem(C45 / "problem_d155_euler_t0cluster_u5em3.yaml")
mB_new, mB_old = RA.mesh_params_euler(pB, 0.0768075), _old_mesh_params(pB, 0.0768075, 321, 65, 5.0e-3)
check("(e) E2 の腕 B: mesh_euler (標準) と変更前の mesh (両 0.005 + ブレンド) は Mesh2DParams は違うが生成座標・接続がビット一致",
      mB_new != mB_old and same_mesh(coords_of(mB_new), coords_of(mB_old)))
p44 = load_problem(REPO / "case/44.vitiated_air_wt/problem_va3_M4.19_Lc8_dry_lumpX.yaml")   # run_species_attrs_ic_tests の Euler prepare
check("(e) case/44 va3 dry lumpX (試験が使う Euler の問題) は記録の格子を明示し、変更前の Euler 経路と同一",
      RA.mesh_params_euler(p44, 0.05) == _old_mesh_params(p44, 0.05, 321, 65, 5.0e-3))
eul = [f.name for f in C45.glob("problem_d155*.yaml") if "_ns" not in f.name]
miss = [f for f in eul if "mesh_euler" not in yaml.safe_load((C45 / f).read_text())]
check(f"(e) case/45 の Euler の問題 ({len(eul)} 本) はすべて mesh_euler を持つ", not miss and set(eul) == set(STD) | set(KEEP), miss)
ns_with = [f.name for f in ns_files if "mesh_euler" in yaml.safe_load(f.read_text())]
check("(e) case/45 の NS の問題には mesh_euler を足していない", not ns_with, ns_with)

print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
