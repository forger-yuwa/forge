#!/usr/bin/env python3
"""格子の来歴と評価方式の識別の試験 (plan tooling-sern-te-wake-grid §4「格子の来歴と評価方式の識別を分ける」・§5.1 #3、
codex plan レビュー 2026-10-08 M3)。forge は回さない (変換器があれば粗い格子を runner の prepare で作る)。

合格 (測る前に固定):
  (a) 同じ評価方式の異なる設計点 (L_cowl が違い格子署名が違う) は一緒に学習できる: runner の prepare で L_cowl だけ違う粗い格子を
      作ると (2D: 1.0 / 1.1、3D: L_b 1.0 H で 2 設計点)、評価方式の id が同じ (= キャンペーンの要求 eval_method_required) で
      格子署名・双対幾何のハッシュは違う。collect → 台帳の作動点要約 → 修正後の系列のキャンペーンで両方の行が学習に入る。
      id は設計 (MOC の θ_b・θ_e) に依らず、格子パラメータ・明示の角度・格子変形の L_b で変わる
  (b) 同じ形状でも旧方式・来歴不明は混ぜない: 修正後の系列 (L_b 1.0) で、格子変形なし (識別の導入前の行・off)・変換器の版違い・
      メッシャの版違い・レシピ違い・不明 (None)・prepare の後に格子が変わった評価を外す。修正前の系列 (off) は格子変形を入れた
      評価・不明を外し、識別の導入前の行は従来どおり採る
  (c) stage_key: 後縁下流の格子変形の OFF→ON・双対幾何だけの変更 (座標・接続は同じ = 格子署名が同じ)・曲線版・L_b で区間が
      分かれる。/VALUE・/AUX の書き換えでは分かれない。キーが無い config・問題 YAML (格子に宣言の属性が無い・0) の stage_key と、
      te_wake_blend_H の有無によらず生成 config は変更前 (commit 9c930e0a) とバイト一致
  (d) 旧評価の持ち越しの承認 (opt.eval_method_carry_over: 作動点・旧評価の評価方式・plan §6 #6・#7 の判定・根拠) は識別の一致と
      独立: 承認が無ければ旧評価は外す、承認 (両判定 PASS) があれば名指しの作動点の該当する旧評価だけ採る、PASS でない承認は効かない、
      承認は処置以外の採否を緩めない、識別が一致する評価に承認は要らない、不正な記録はキャンペーンの作成で止まる
  (e) 来歴: metrics.json に格子署名・双対幾何のハッシュ (いま再計算した値と prepare の記録) と評価方式の識別。双対幾何のハッシュは
      格納の型・圧縮・属性・/VALUE に依らない。te_wake_grid_check の --provenance-out と --same-mesh に双対幾何のハッシュ
      (署名が同じで双対だけ違えば NO)

  design/.venv-opt/bin/python design/tests/run_sern_mesh_provenance_tests.py
"""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
CHECK = ROOT / "case" / "46.sern_design" / "diag" / "te_wake_grid_check.py"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(TOOLS))
try:
    import pymoo  # noqa: F401
except ImportError:   # 選別関数だけを試すので、最適化モジュールは空のスタブで足りる
    _m = types.ModuleType("forge_design.opt.moo"); _m.propose_infill = None
    sys.modules["forge_design.opt.moo"] = _m
try:
    import smt  # noqa: F401
except ImportError:
    _s = types.ModuleType("forge_design.opt.surrogate"); _s.KrigingSet = object
    sys.modules["forge_design.opt.surrogate"] = _s

import mark_zero_thickness_edges as MZ  # noqa: E402
import stage_manifest as SM  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.meshing.mesh_sern3d import TE_WAKE_CURVE_VERSION  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
from forge_design.probdef import Problem, load_problem  # noqa: E402

REF_COMMIT = "9c930e0a"     # 本変更の直前の HEAD
fails = 0
skipped = []
tmp = Path(tempfile.mkdtemp(prefix="mesh_prov_"))


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""))


def raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    except Exception as e:  # noqa: BLE001
        print(f"    (別の例外 {type(e).__name__}: {e})")
        return False
    return False


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_ref(rel, modname, package=None):
    """変更前のソースを git から読んで別名のモジュールにする (相対 import は package で解決)。読めなければ None。"""
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{REF_COMMIT}:{rel}"], capture_output=True, text=True)
    if r.returncode != 0:
        return None
    d = tmp / "ref" / "forge_design" / "evaluate"      # 旧ソースの FORGE_ROOT の既定 (parents[3]) が評価できる深さに置く
    d.mkdir(parents=True, exist_ok=True)
    f = d / (modname.replace(".", "_") + ".py")
    f.write_text(r.stdout)
    spec = importlib.util.spec_from_file_location(modname, f)
    m = importlib.util.module_from_spec(spec)
    if package:
        m.__package__ = package
    saved = os.environ.get("FORGE_ROOT")
    os.environ["FORGE_ROOT"] = str(ROOT)
    try:
        spec.loader.exec_module(m)
    finally:
        if saved is None:
            os.environ.pop("FORGE_ROOT", None)
        else:
            os.environ["FORGE_ROOT"] = saved
    return m


def prob(evaluate=None, mesh=None, mesh3d=None):
    raw = {"gas": {}}
    if mesh3d is not None:
        raw["mesh3d"] = dict(mesh3d)
    return Problem(name="t", type="sern_2d", gamma=1.4, cp=1004.5, spec={"H_m": 0.1}, dv={}, geometry={},
                   mesh=dict(mesh or {}), evaluate=dict(evaluate or {}), raw=raw)


# ------------------------------------------------------------------------------------------ 合成の node 格子
COORD = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]], float)
EDGES = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]


def make_node_h5(path, coord=COORD, dtype=np.float32, compress=False, surf_bump=0.0, drop=(), attrs=None, val=1.0):
    """mark_zero_thickness_edges の格子署名が読む構造 (MESH・PLANES/STRUCT・BCONDS・VIZMESH) と、双対幾何 (PLANES/surfVect・
    surfArea・centCoords、CELLS/volume・centCoords) を持つ 2D の小さな node 格子。surf_bump は面 0 の面積ベクトルだけを
    相対で変える (座標・接続はそのまま = 格子署名は同じ)。"""
    n, nf = len(coord), len(EDGES)
    kw = {"compression": "gzip"} if compress else {}
    sv = np.array([[0.0, 1.0, 0.0] if abs(coord[a][1] - coord[b][1]) < 1e-12 else [1.0, 0.0, 0.0] for a, b in EDGES])
    sv[0] *= 1.0 + surf_bump
    geo = {"PLANES/surfVect": sv.ravel(), "PLANES/surfArea": np.linalg.norm(sv, axis=1),
           "PLANES/centCoords": np.array([0.5 * (coord[a] + coord[b]) for a, b in EDGES]).ravel(),
           "CELLS/volume": np.full(n, 0.25), "CELLS/centCoords": np.asarray(coord, float).ravel()}
    with h5py.File(path, "w") as f:
        g = f.create_group("MESH")
        for k_, v_ in (("nNodes", n), ("nCells", n), ("nNormalPlanes", nf), ("nPlanes", nf), ("nBconds", 2), ("nBPlanes", 0)):
            g.attrs[k_] = v_
        g.create_dataset("COORD", data=np.asarray(coord, np.float32).ravel())
        f.create_dataset("PLANES/STRUCT", data=np.array([x for a_, b_ in EDGES for x in (2, a_, b_, 2, a_, b_)], np.int32))
        for k_, v_ in geo.items():
            if k_ not in drop:
                f.create_dataset(k_, data=np.asarray(v_, dtype), **kw)
        v = f.create_group("VIZMESH"); v.attrs["nVizCells"] = 1; v.attrs["vizCONNE_dim"] = 5
        v.create_dataset("CONNE", data=np.array([5, 0, 1, 4, 3], np.int32))
        for pid, faces in ((5, [[0, 1]]), (6, [[1, 2]])):
            b_ = f.create_group(f"BCONDS/{pid}"); b_.attrs["bcondKind"] = "wall"
            b_.create_dataset("vizBfaceSizes", data=np.array([len(x) for x in faces], np.int32))
            b_.create_dataset("vizBfaceNodes", data=np.array([i for x in faces for i in x], np.int32))
        for k_ in ("ro", "roUx", "roUy", "roUz", "roe"):
            f.create_dataset(f"VALUE/{k_}", data=np.full(n, val, np.float32))
        for k_, v_ in (attrs or {}).items():
            f.attrs[k_] = v_


ON_ATTRS = {SM.TE_WAKE_ATTR: np.float64(1.0), SM.TE_WAKE_CURVE_ATTR: TE_WAKE_CURVE_VERSION}
BC = "ramp: {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: }\n"


def run_dir_with(name, cfg, h5=True, **kw):
    d = tmp / name
    d.mkdir()
    (d / "solverConfig.yaml").write_text(cfg); (d / "solverConfig_main.yaml").write_text(cfg)
    (d / "bcondConfig.yaml").write_text(BC)
    if h5 is True:
        make_node_h5(d / "sern.h5", **kw)
    elif isinstance(h5, (bytes, bytearray)):
        (d / "sern.h5").write_bytes(h5)
    return d


# ======================================================================== (e) 双対幾何のハッシュ (単体)
a0 = tmp / "dual_a0.h5"; make_node_h5(a0)
dh = SM.dual_geometry_hash
for i_, (label, kw) in enumerate((("格納の型 float64", {"dtype": np.float64}), ("gzip 圧縮", {"compress": True}),
                                  ("ルート属性 (格子変形の宣言)", {"attrs": ON_ATTRS}), ("/VALUE の書き換え", {"val": 3.0}))):
    pth = tmp / f"dual_var{i_}.h5"; make_node_h5(pth, **kw)
    check(f"(e) 双対幾何のハッシュは{label}で変わらない (格子署名も同じ)", dh(pth) == dh(a0) and MZ.mesh_signature(str(pth)) == MZ.mesh_signature(str(a0)))
a_bump = tmp / "dual_bump.h5"; make_node_h5(a_bump, surf_bump=1e-6)
check("(e) 面積ベクトル 1 面だけの変更 (座標・接続は同じ): 格子署名は同じ・双対幾何のハッシュは違う",
      MZ.mesh_signature(str(a_bump)) == MZ.mesh_signature(str(a0)) and dh(a_bump) != dh(a0))
for ds in ("PLANES/surfArea", "PLANES/centCoords", "CELLS/volume", "CELLS/centCoords"):
    pth = tmp / f"dual_drop_{ds.replace('/', '_')}.h5"; make_node_h5(pth, drop=(ds,))
    check(f"(e) {ds} の有無でハッシュが違う", dh(pth) != dh(a0))
a_cv = tmp / "dual_cv.h5"; make_node_h5(a_cv)
with h5py.File(a_cv, "r+") as f:
    f["CELLS/volume"][2] = np.float32(0.26)
check("(e) CV 体積 1 つの変更でハッシュが違う", dh(a_cv) != dh(a0))
mi = SM.mesh_identity(a0)
check("(e) mesh_identity: 格子署名 (道具の再計算)・双対幾何のハッシュ・版・離散化",
      mi["mesh_signature"] == MZ.mesh_signature(str(a0)) and mi["mesh_signature_version"] == MZ.SIG_VERSION
      and mi["dual_hash"] == dh(a0) and mi["dual_hash_version"] == SM.DUAL_HASH_VERSION and mi["discretization"] == "node")
check("(e) mesh_identity: ファイルが無ければ None", SM.mesh_identity(tmp / "no_such.h5") is None)

# ======================================================================== (c) stage_key
REFSM = load_ref("solver_density_cuda/tools/stage_manifest.py", "_stage_manifest_ref")
REF2 = load_ref("design/forge_design/evaluate/runner_sern.py", "forge_design.evaluate._runner_sern_ref", "forge_design.evaluate")
VARIANTS = [("euler cell (既定)", {}, {}), ("euler node", {}, {"discretization": "node"}),
            ("sst node 生産風", {"model": "sst", "limiter_ref": {"length": 1.0, "ro": 0.1, "p": 1e4, "a": 500.0}, "implicit_relax": 0.7,
                               "p_min": 1.0, "wall_treatment_sst": 0, "floor_events": 1}, {"discretization": "node", "node_inlet_corner_wall": 1}),
            ("sst node 処置あり", {"model": "sst", "zero_thickness_edge_velocity": {"tags": ["cowl_in", "cowl_out"], "rings": 2}},
             {"discretization": "node"})]
if REF2 is None or REFSM is None:
    skipped.append("(c) git に変更前の版が無い")
    print("SKIP (c) バイト一致: git show が使えない")
else:
    for name, ev, ms in VARIANTS:
        old = REF2._solver_config(prob(ev, ms), 4000, 500, 2.5, 2851.0)
        same = True
        for lb_label, m3, m2x in (("キー無し", None, {}), ("mesh3d 0", {"te_wake_blend_H": 0.0}, {}), ("mesh3d 1.0", {"te_wake_blend_H": 1.0}, {}),
                                  ("mesh 1.0 (2D)", None, {"te_wake_blend_H": 1.0})):
            p_ = prob(ev, {**ms, **m2x}, m3)
            same &= (R2._solver_config(p_, 4000, 500, 2.5, 2851.0) == old and R3._solver_config(p_, 4000, 500, 2.5, 2851.0) == old)
        check(f"(c) {name}: 生成 config (2D/3D 本段) は te_wake_blend_H の有無・値によらず変更前とバイト一致", same)
        disc = ms.get("discretization", "cell")
        check(f"(c) {name}: 品質確認用 cell の config も変更前とバイト一致", R2.qc_cell_config(old, disc) == REF2.qc_cell_config(old, disc))
    CFGS = [R2._solver_config(prob(ev, ms), 4000, 500, 2.5, 2851.0) for _, ev, ms in VARIANTS] + [
        "turbulence:\n  model: sst\nspace:\n  convMethod: 1\n  limiter: 2\nmesh:\n  scalarGradient: lsq\n  meshFileName: sern.h5\n",
        "conjugate: {mode: 1, flux: q}\nspace: {convMethod: 1}\n", "space: {convMethod: 1\n"]
    DIRS = [("run_dir None", None)]
    for label, kw in (("格子に宣言なし", {}), ("宣言 0", {"attrs": {SM.TE_WAKE_ATTR: np.float64(0.0)}}),
                      ("宣言 0 (整数)", {"attrs": {SM.TE_WAKE_ATTR: 0}})):
        DIRS.append((label, run_dir_with("sk_off_" + str(len(DIRS)), CFGS[2], **kw)))
    DIRS.append(("HDF5 でない sern.h5", run_dir_with("sk_bad", CFGS[2], h5=b"not an hdf5 file")))
    DIRS.append(("sern.h5 無し", run_dir_with("sk_none", CFGS[2], h5=False)))
    def _key(mod, c, rdir):
        try:
            return json.dumps(mod.stage_key(c, BC, rdir), sort_keys=True)
        except Exception as e:  # noqa: BLE001   (処置ありの段で h5 が読めない等は、変更前と同じ例外になることを確かめる)
            return "raise " + type(e).__name__
    for i, c in enumerate(CFGS):
        for label, rdir in DIRS:
            check(f"(c) キーが無い段 #{i} ({label}): stage_key が変更前とバイト一致", _key(SM, c, rdir) == _key(REFSM, c, rdir))

cfg = R2._solver_config(prob({"model": "sst"}, {"discretization": "node"}), 4000, 500, 2.5, 2851.0)
MK = lambda k: {a: b for a, b in k.items() if a.startswith(SM.MESH_PREFIX)}
d_off = run_dir_with("sk2_off", cfg)
d_on = run_dir_with("sk2_on", cfg, attrs=ON_ATTRS)
k_off, k_on = SM.stage_key(cfg, BC, d_off), SM.stage_key(cfg, BC, d_on)
check("(c) 宣言の無い段に格子の来歴のキーは足さない", not MK(k_off))
check("(c) OFF → ON (格子に L_b 1.0 の宣言) で別キー", k_off != k_on)
check("(c) ON の段のキー = 宣言 (L_b・曲線版) と再計算した格子署名・双対幾何のハッシュ",
      MK(k_on) == {SM.MESH_PREFIX + ".te_wake_blend_H": "1.0", SM.MESH_PREFIX + ".te_wake_curve_version": TE_WAKE_CURVE_VERSION,
                   SM.MESH_PREFIX + ".mesh_signature": MZ.mesh_signature(str(d_on / "sern.h5")),
                   SM.MESH_PREFIX + ".dual_hash": dh(d_on / "sern.h5")}, MK(k_on))
seg = SM.segments({"stages": [{"tag": t, "key": k} for t, k in (("off1", k_off), ("off2", k_off), ("on1", k_on), ("on2", k_on))]})
check("(c) OFF,OFF,ON,ON → 2 区間 (ON 同士は連結)", [[s["tag"] for s in sg] for sg in seg] == [["off1", "off2"], ["on1", "on2"]])
d_dual = run_dir_with("sk2_dual", cfg, attrs=ON_ATTRS, surf_bump=1e-6)
k_dual = SM.stage_key(cfg, BC, d_dual)
check("(c) 双対幾何だけの変更 (座標・接続は同じ = 格子署名が同じ) で別キー",
      MK(k_dual)[SM.MESH_PREFIX + ".mesh_signature"] == MK(k_on)[SM.MESH_PREFIX + ".mesh_signature"] and k_dual != k_on
      and {a for a in k_on if k_on[a] != k_dual.get(a)} == {SM.MESH_PREFIX + ".dual_hash"})
check("(c) 双対幾何の変更は区間を分ける", len(SM.segments({"stages": [{"tag": "a", "key": k_on}, {"tag": "b", "key": k_dual}]})) == 2)
for label, at in (("曲線版の違い", {SM.TE_WAKE_ATTR: np.float64(1.0), SM.TE_WAKE_CURVE_ATTR: "hermite3-lower-tangent-v0"}),
                  ("L_b 0.5", {SM.TE_WAKE_ATTR: np.float64(0.5), SM.TE_WAKE_CURVE_ATTR: TE_WAKE_CURVE_VERSION})):
    d_ = run_dir_with("sk2_" + label[:2] + str(len(label)), cfg, attrs=at)
    check(f"(c) {label} で別キー", SM.stage_key(cfg, BC, d_) != k_on)
for label, at in (("負", -1.0), ("文字列", "1.0"), ("非有限", float("nan"))):
    d_ = run_dir_with("sk2_bad_" + str(len(label)) + label[:1], cfg, attrs={SM.TE_WAKE_ATTR: at})
    k_ = SM.stage_key(cfg, BC, d_)
    check(f"(c) 読めない宣言 ({label}) は無効と連結しない (invalid)",
          k_ != k_off and str(k_.get(SM.MESH_PREFIX + ".te_wake_blend_H")).startswith("invalid"), k_.get(SM.MESH_PREFIX + ".te_wake_blend_H"))
d_nc = run_dir_with("sk2_nocurve", cfg, attrs={SM.TE_WAKE_ATTR: np.float64(1.0)})
check("(c) 曲線版の属性が無い宣言は missing (ON の別キー)", MK(SM.stage_key(cfg, BC, d_nc)).get(SM.MESH_PREFIX + ".te_wake_curve_version") == "missing")
with h5py.File(d_on / "sern.h5", "r+") as f:
    f["VALUE/ro"][...] = 7.0
    f.create_dataset("AUX/w_recon_vel", data=np.ones(len(COORD), np.float32))
check("(c) /VALUE・/AUX の書き換え (restart_field・処置の w) では同じキー", SM.stage_key(cfg, BC, d_on) == k_on)
d_rep = run_dir_with("sk2_replace", cfg)
man = SM.StageManifest(d_rep); man.add("soft", cfg, BC); man.add("mid", cfg, BC)
make_node_h5(d_rep / "sern.h5", attrs=ON_ATTRS)          # 同じ run で格子を L_b 1.0 の格子に差し替えて続ける
man.add("main", cfg, BC); man.write()
check("(c) StageManifest: 同じ run の中で格子を ON に差し替えた段は別区間",
      [[s["tag"] for s in sg] for sg in SM.segments(json.load(open(d_rep / "stage_manifest.json")))] == [["soft", "mid"], ["main"]])
_w0 = tmp / "w0c.h5"; make_node_h5(_w0); _b0 = _w0.read_bytes()
check("(c) 宣言の書き込み: minfo にキーが無い (2D) / 0 は書かない", R2.write_te_wake_attrs(_w0, {}) is None and R2.write_te_wake_attrs(_w0, {"te_wake_blend_H": 0.0}) is None
      and _w0.read_bytes() == _b0)
_w1 = tmp / "w1.h5"; make_node_h5(_w1)
R2.write_te_wake_attrs(_w1, {"te_wake_blend_H": 1.0, "te_wake_curve_version": TE_WAKE_CURVE_VERSION})
check("(c) 宣言の書き込み: L_b 1.0 は stage_manifest が読める宣言になる",
      SM.te_wake_declared(_w1) == {"te_wake_blend_H": "1.0", "te_wake_curve_version": TE_WAKE_CURVE_VERSION})

# ======================================================================== (a)(b) 評価方式の識別の単体 (設計に依らない・格子の作り方で変わる)
p3 = load_problem(ROOT / "case" / "46.sern_design" / "problem_3d_prod_m6on_g1.yaml")
CONV = {"sha256": "c" * 64}


def eid3(p_, theta_b=0.0, theta_e=0.0, conv=CONV):
    prm = R3.sern_mesh3d_params(p_, theta_b, theta_e)
    L_b, curve = R2.te_wake_spec(p_.raw, 3)
    return R2.eval_method_of(3, "node", R2.mesh_recipe(prm, p_.mesh), R2.mesher_version(3), L_b, curve, conv["sha256"])[1]


def with_raw(p_, mesh=None, mesh3d=None):
    raw = json.loads(json.dumps(p_.raw))
    raw["mesh"].update(mesh or {}); raw.setdefault("mesh3d", {}).update(mesh3d or {})
    return Problem(name=p_.name, type=p_.type, gamma=p_.gamma, cp=p_.cp, spec=raw["spec"], dv=raw["dv"], geometry=raw["geometry"],
                   mesh=raw["mesh"], evaluate=raw["evaluate"], raw=raw)


base_id = eid3(p3, -0.70, 0.05)
check("(a) 評価方式の id は設計の角度 (θ_b・θ_e) に依らない", eid3(p3, -0.50, 0.10) == base_id)
check("(a) 評価方式の id は dv (L_cowl) に依らない", eid3(with_raw(p3)) == base_id
      and (lambda r: eid3(Problem(name="x", type="sern_2d", gamma=1.4, cp=1004.5, spec=r["spec"], dv={**r["dv"], "L_cowl": 1.37},
                                  geometry=r["geometry"], mesh=r["mesh"], evaluate=r["evaluate"], raw=r)))(json.loads(json.dumps(p3.raw))) == base_id)
for label, kw in (("格子パラメータ (ni_noz)", {"mesh3d": {"ni_noz": 91}}), ("第一層厚", {"mesh3d": {"first_wall_frac": 2.0e-3}}),
                  ("格子変形 L_b 1.0", {"mesh3d": {"te_wake_blend_H": 1.0}}), ("明示の interface_angle_rad", {"mesh": {"interface_angle_rad": -0.7}})):
    check(f"(b) {label} で評価方式の id が変わる", eid3(with_raw(p3, **kw)) != base_id)
check("(b) 明示の L_b 0 はキー無しと同じ格子 (id は同じ)", eid3(with_raw(p3, mesh3d={"te_wake_blend_H": 0.0})) == base_id)
check("(b) 変換器の版違いで id が変わる、版が不明なら id は None",
      eid3(p3, conv={"sha256": "d" * 64}) != base_id and R2.eval_method_of(3, "node", {}, R2.mesher_version(3), 0.0, None, None)[1] is None)
check("(b) 2D と 3D は別の評価方式 (メッシャの版も別)", R2.mesher_version(2) != R2.mesher_version(3))
check("(b) 格子変形の識別子", R2.te_wake_effective_of(0.0, None) == "off"
      and R2.te_wake_effective_of(1.0, TE_WAKE_CURVE_VERSION) == f"L_b=1.0;curve={TE_WAKE_CURVE_VERSION}"
      and R2.te_wake_spec({"mesh3d": {"te_wake_blend_H": 1.0}}, 3) == (1.0, TE_WAKE_CURVE_VERSION)
      and R2.te_wake_spec({"mesh": {"te_wake_blend_H": 1.0}}, 2) == (1.0, None) and R2.te_wake_spec({}, 3) == (0.0, None))

# ======================================================================== 実際の runner (変換器があるとき): 3D・2D の粗い格子
class Fake:
    """キャンペーンの偽物 (選別に要る属性だけ)。"""

    def __init__(self, rows, te="off", req=None, co=()):
        self.rows, self.te_wake_required, self.eval_method_required, self.eval_method_carry_over = rows, te, req, tuple(co)
        self.ref, self.ops = (-0.9, 20.0), []


HAVE_CONV = R2.CONVERTER.exists()
runs = {}


def stub_collect(mod, rd, prob_path, cfg_main):
    """力の履歴とゲートを差し替えて collect (metrics.json) を回す。起動記録は本段の config で 1 行。"""
    with open(rd / "forge_launches.jsonl", "a") as fh:
        fh.write(json.dumps({"time": 0, "cfg_fnv": SM.fnv1a64(cfg_main), "slauWallNormalChi": 1, "scalarGradient": "lsq"}) + "\n")
    gates = lambda *a, **k: {"steadiness": {"series": {}}, "objective": "C_T", "verdict": "PASS", "fail_class": None, "reasons": []}
    saved = (R2.force_history, R2.evaluate_gates, R3.evaluate_gates)
    R2.force_history = lambda *a, **k: []; R2.evaluate_gates = gates; R3.evaluate_gates = gates
    try:
        out = mod.collect(prob_path, rd, rc=0)
    finally:
        R2.force_history, R2.evaluate_gates, R3.evaluate_gates = saved
    return out


if not HAVE_CONV:
    skipped.append("変換器が無い: 実際の prepare の試験")
    print(f"SKIP 実際の prepare: {R2.CONVERTER} が無い")
else:
    # ---- 3D: 生産の g1 を粗くした格子。A = キー無し、B1 = L_b 1.0、B2 = L_b 1.0 で L_cowl を変えた設計点
    src3 = yaml.safe_load(open(ROOT / "case" / "46.sern_design" / "problem_3d_prod_m6on_g1.yaml"))
    cache = {}
    _saved_design = R2.design_from_problem

    def _cached_design(p_, design=None):
        key = json.dumps(p_.dv, sort_keys=True, default=str)
        if key not in cache:
            cache[key] = _saved_design(p_, design=design)
        return cache[key]
    R2.design_from_problem = _cached_design
    try:
        for lab, lb, dLc in (("A", None, 0.0), ("B1", 1.0, 0.0), ("B2", 1.0, 0.05)):
            y = json.loads(json.dumps(src3))
            y["mesh3d"].update({"nz_in": 4, "nz_out": 3, "ni_noz": 30, "nj_vside": 5, "first_wake_frac": 0.004, "first_z_frac": 0.02})
            y["mesh"].update({"nj_ext_top": 9, "first_top_frac": 0.004})
            if lb is not None:
                y["mesh3d"]["te_wake_blend_H"] = lb
            y["dv"]["L_cowl"] = float(y["dv"]["L_cowl"]) + dLc
            pth = tmp / f"p3_{lab}.yaml"; pth.write_text(yaml.safe_dump(y, allow_unicode=True))
            info = R3.prepare(pth, tmp / f"r3_{lab}")
            runs[lab] = (tmp / f"r3_{lab}", pth, info)
    finally:
        R2.design_from_problem = _saved_design
    (rA, pA, iA), (rB1, pB1, iB1), (rB2, pB2, iB2) = runs["A"], runs["B1"], runs["B2"]
    print(f"--- 3D 粗い格子: A {iA['mesh']['nodes']} 節点、B1 {iB1['mesh']['nodes']}、B2 {iB2['mesh']['nodes']}")
    with h5py.File(rA / "sern.h5", "r") as f:
        a_attrs = dict(f.attrs)
    check("(c) 3D キー無しの問題 YAML: 格子に格子変形の宣言を書かない", SM.TE_WAKE_ATTR not in a_attrs and SM.TE_WAKE_CURVE_ATTR not in a_attrs)
    check("(c) 3D L_b 1.0: 格子に L_b と曲線版の宣言を書く",
          SM.te_wake_declared(rB1 / "sern.h5") == {"te_wake_blend_H": "1.0", "te_wake_curve_version": TE_WAKE_CURVE_VERSION}
          and iB1["mesh_provenance"]["te_wake_attrs"] == {SM.TE_WAKE_ATTR: 1.0, SM.TE_WAKE_CURVE_ATTR: TE_WAKE_CURVE_VERSION})
    cA, cB1 = (rA / "solverConfig_main.yaml").read_text(), (rB1 / "solverConfig_main.yaml").read_text()
    bA, bB1 = (rA / "bcondConfig.yaml").read_text(), (rB1 / "bcondConfig.yaml").read_text()
    if REFSM is not None:
        check("(c) 3D キー無しの問題 YAML の run: stage_key が変更前とバイト一致",
              json.dumps(SM.stage_key(cA, bA, rA), sort_keys=True) == json.dumps(REFSM.stage_key(cA, bA, rA), sort_keys=True))
    kA, kB1 = SM.stage_key(cA, bA, rA), SM.stage_key(cB1, bB1, rB1)
    check("(c) 3D: config・bcond は同じで、OFF → ON の stage_key は別 (来歴のキーだけ増える)",
          cA == cB1 and bA == bB1 and kA != kB1 and {a: b for a, b in kB1.items() if not a.startswith(SM.MESH_PREFIX)} == kA
          and MK(kB1)[SM.MESH_PREFIX + ".dual_hash"] == iB1["mesh_provenance"]["dual_hash"])
    eA, eB1, eB2 = (i["eval_method"] for i in (iA, iB1, iB2))
    reqA = R2.eval_method_required(load_problem(pA), 3)
    reqB = R2.eval_method_required(load_problem(pB1), 3)
    check("(a) 3D: L_b 1.0 の 2 設計点 (L_cowl 違い) は評価方式の id が同じ・格子署名と双対幾何のハッシュは違う",
          eB1[R2.EVAL_METHOD_ID] == eB2[R2.EVAL_METHOD_ID] is not None
          and iB1["mesh_provenance"]["mesh_signature"] != iB2["mesh_provenance"]["mesh_signature"]
          and iB1["mesh_provenance"]["dual_hash"] != iB2["mesh_provenance"]["dual_hash"])
    check("(a) 3D: prepare の id = キャンペーンの要求 (問題 YAML・現在のコード・現在の変換器、MOC を解かない)",
          reqB[R2.EVAL_METHOD_ID] == eB1[R2.EVAL_METHOD_ID] and reqA[R2.EVAL_METHOD_ID] == eA[R2.EVAL_METHOD_ID])
    check("(b) 3D: 同じ設計の旧方式 (キー無し) は別の id・格子変形の実効値 off / L_b 1.0",
          eA[R2.EVAL_METHOD_ID] != eB1[R2.EVAL_METHOD_ID] and eA[R2.TE_WAKE_EFFECTIVE] == "off"
          and eB1[R2.TE_WAKE_EFFECTIVE] == f"L_b=1.0;curve={TE_WAKE_CURVE_VERSION}")
    check("(e) 3D: prepare_info に変換器の識別 (実行したバイナリの sha256)",
          eB1["converter"]["sha256"] == R2.file_sha256(R2.CONVERTER) == eB1["eval_method"]["converter_sha256"])
    out3 = {lab: stub_collect(R3, runs[lab][0], runs[lab][1], (runs[lab][0] / "solverConfig_main.yaml").read_text()) for lab in runs}
    m3 = json.loads((rB1 / "metrics.json").read_text())
    check("(e) 3D collect: metrics.json に格子の来歴 (いまの再計算 = prepare の記録) と評価方式の id",
          m3["mesh_provenance"]["same_as_prepare"] is True and m3["mesh_provenance"]["now"]["dual_hash"] == iB1["mesh_provenance"]["dual_hash"]
          and m3[R2.EVAL_METHOD_ID] == eB1[R2.EVAL_METHOD_ID] and m3[R2.TE_WAKE_EFFECTIVE] == eB1[R2.TE_WAKE_EFFECTIVE])
    # 双対幾何だけを変えた格子 (同じ run で h5 を差し替えた想定)
    rBd = tmp / "r3_B1_dual"; shutil.copytree(rB1, rBd)
    with h5py.File(rBd / "sern.h5", "r+") as f:
        sv = f["PLANES/surfVect"]; a_ = np.asarray(sv[...]); j_ = int(np.argmax(np.abs(a_)))   # 0 でない成分を 1 つだけ変える
        a_.flat[j_] *= 1.0 + 1e-3; sv[...] = a_
    check("(c) 3D: 双対幾何だけの変更は格子署名が同じで stage_key が別",
          MZ.mesh_signature(str(rBd / "sern.h5")) == MZ.mesh_signature(str(rB1 / "sern.h5")) and SM.stage_key(cB1, bB1, rBd) != kB1)
    emd = R2.eval_method_provenance(rBd)
    check("(b) 3D: prepare の後に双対幾何が変わった評価は評価方式が不明 (None)", emd[R2.EVAL_METHOD_ID] is None and "変わった" in emd["reason"], emd["reason"])
    TC = load_module(CHECK, "te_wake_grid_check_prov")
    e_ = TC.provenance_entry(str(rB1))
    check("(e) te_wake_grid_check の来歴に双対幾何のハッシュ (再計算)", e_.get("dual_hash") == iB1["mesh_provenance"]["dual_hash"]
          and e_.get("dual_hash_version") == SM.DUAL_HASH_VERSION)
    r_ = subprocess.run([sys.executable, str(CHECK), "--same-mesh", str(rB1), str(rBd)], capture_output=True, text=True)
    check("(e) --same-mesh: 署名が同じで双対幾何だけ違う → SAME MESH: NO (終了コード 1)",
          r_.returncode == 1 and "SAME MESH: NO" in r_.stdout and "双対幾何だけ" in r_.stdout, r_.stdout[-400:] + r_.stderr[-400:])
    r_ = subprocess.run([sys.executable, str(CHECK), "--same-mesh", str(rB1), str(rB1 / "sern.h5")], capture_output=True, text=True)
    check("(e) --same-mesh: 同じ格子 → YES", r_.returncode == 0 and "SAME MESH: YES" in r_.stdout, r_.stdout[-300:])
    # 識別の導入前の run (prepare_info に記録が無い) の来歴
    rLeg = tmp / "r3_A_legacy"; shutil.copytree(rA, rLeg)
    li = json.loads((rLeg / "prepare_info.json").read_text()); li.pop("eval_method"); li.pop("mesh_provenance")
    li["mesh"]["te_wake_blend_H"] = 1.0; li["mesh"]["te_wake_curve_version"] = TE_WAKE_CURVE_VERSION
    (rLeg / "prepare_info.json").write_text(json.dumps(li))
    eml = R2.eval_method_provenance(rLeg)
    check("(b) 識別の導入前の run: id は None (来歴不明)、格子変形は info の mesh から (L_b 1.0)",
          eml[R2.EVAL_METHOD_ID] is None and eml[R2.TE_WAKE_EFFECTIVE] == f"L_b=1.0;curve={TE_WAKE_CURVE_VERSION}" and "導入前" in eml["reason"])

    # ---- 2D: driver の評価経路。L_cowl 1.0 / 1.1 の粗い node 格子
    src2 = yaml.safe_load(open(ROOT / "case" / "46.sern_design" / "problem_moo_sst_node_3op.yaml"))
    for Lc in (1.0, 1.1):
        y = json.loads(json.dumps(src2))
        y["mesh"].update({"ni_noz": 40, "ni_plume": 60, "nj_top": 31, "nj_bot": 21, "first_wall_frac": 0.01, "ar_max": 5000})
        y["dv"]["L_cowl"]["value"] = Lc
        pth = tmp / f"p2_{Lc}.yaml"; pth.write_text(yaml.safe_dump(y, allow_unicode=True))
        runs[f"2d_{Lc}"] = (tmp / f"r2_{Lc}", pth, R2.prepare(pth, tmp / f"r2_{Lc}", op="cruise"))
    (r21, p21, i21), (r22, p22, i22) = runs["2d_1.0"], runs["2d_1.1"]
    req2 = R2.eval_method_required(load_problem(p21), 2)
    check("(a) 2D: L_cowl 1.0 / 1.1 は評価方式の id が同じ (= キャンペーンの要求)・格子署名は違う",
          i21["eval_method"][R2.EVAL_METHOD_ID] == i22["eval_method"][R2.EVAL_METHOD_ID] == req2[R2.EVAL_METHOD_ID] is not None
          and i21["mesh_provenance"]["mesh_signature"] != i22["mesh_provenance"]["mesh_signature"])
    check("(c) 2D: 格子に宣言を書かない・stage_key は変更前とバイト一致", SM.te_wake_declared(r21 / "sern.h5") is None and (REFSM is None or json.dumps(
        SM.stage_key((r21 / "solverConfig_main.yaml").read_text(), (r21 / "bcondConfig.yaml").read_text(), r21), sort_keys=True)
        == json.dumps(REFSM.stage_key((r21 / "solverConfig_main.yaml").read_text(), (r21 / "bcondConfig.yaml").read_text(), r21), sort_keys=True)))
    o21 = stub_collect(R2, r21, p21, (r21 / "solverConfig_main.yaml").read_text())
    o22 = stub_collect(R2, r22, p22, (r22 / "solverConfig_main.yaml").read_text())
    s21, s22 = D.SernCampaign._op_summary(o21, r21), D.SernCampaign._op_summary(o22, r22)
    check("(a) 2D: 台帳の作動点要約に評価方式の id・格子変形の実効値・格子署名・双対幾何のハッシュ",
          s21[D.EVAL_METHOD_ID] == req2[R2.EVAL_METHOD_ID] and s21[D.TE_WAKE_EFFECTIVE] == "off"
          and s21["mesh_signature"] == i21["mesh_provenance"]["mesh_signature"] and s21["dual_hash"] == i21["mesh_provenance"]["dual_hash"])
    rows2 = [{"status": "PASS", "flag_policy": D.FLAG_POLICY, "x": [float(Lc)], "C_T_w": 0.9, "L_ramp": 1.0, "C_M_w": 0.0, "ops": {"m6_on": s}}
             for Lc, s in ((1.0, s21), (1.1, s22))]
    check("(a) 2D: 同じ評価方式の 2 設計点は、その id を要求するキャンペーンで一緒に学習に入る",
          len(D.SernCampaign._XF(Fake(rows2, "L_b=1.0;curve=x", req2[R2.EVAL_METHOD_ID]))[0]) == 2)
    check("(b) 2D: 同じ 2 行も、別の id (変換器の版違い) を要求するキャンペーンでは学習に入らない",
          len(D.SernCampaign._XF(Fake(rows2, "L_b=1.0;curve=x", "e" * 64))[0]) == 0)
    check("(b) 2D: 修正前の系列 (off) では識別を記録した off の行も従来どおり学習に入る", len(D.SernCampaign._XF(Fake(rows2))[0]) == 2)
    rowsB = [{"status": "PASS", "flag_policy": D.FLAG_POLICY, "x": [0.0], "C_T_w": 0.9, "L_ramp": 1.0, "C_M_w": 0.0,
              "ops": {"m6_on": {**D.SernCampaign._op_summary(out3[lab], runs[lab][0]), "scalar_gradient_effective": "lsq",
                                "slau_wall_normal_chi_effective": 1}}} for lab in ("B1", "B2")]
    check("(a) 3D: L_b 1.0 の 2 設計点の評価は、L_b 1.0 のキャンペーンで一緒に学習に入る",
          len(D.SernCampaign._XF(Fake(rowsB, eB1[R2.TE_WAKE_EFFECTIVE], reqB[R2.EVAL_METHOD_ID]))[0]) == 2)
    check("(b) 3D: L_b 1.0 の評価は修正前の系列 (off) では学習に入らない", len(D.SernCampaign._XF(Fake(rowsB))[0]) == 0)

# ======================================================================== (b)(d) 設計 DB の選別 (作動点要約を組み立てる)
ID_NEW, ID_OLD, ID_CONV2 = "a" * 64, "b" * 64, "f" * 64
ON = f"L_b=1.0;curve={TE_WAKE_CURVE_VERSION}"
OPS2 = ("m6_on", "m10_on")


def op_(mid="new", te=None, sg="lsq"):
    o = {"scalar_gradient_effective": sg, "slau_wall_normal_chi_effective": 1, "mesh_signature": os.urandom(8).hex()}
    if mid != "legacy":
        o[D.EVAL_METHOD_ID] = {"new": ID_NEW, "old": ID_OLD, "conv2": ID_CONV2, "none": None}[mid]
        o[D.TE_WAKE_EFFECTIVE] = te if te is not None else (ON if mid == "new" else "off")
    return o


def row_(ops, x=0.0):
    return {"status": "PASS", "flag_policy": D.FLAG_POLICY, "x": [x], "C_T_w": 0.9, "L_ramp": 1.0, "C_M_w": 0.0, "ops": ops}



def n_learn(rows, te="off", req=None, co=()):
    return len(D.SernCampaign._XF(Fake(rows, te, req, co))[0])


A6 = {"op": "m6_on", "from_eval_method": "legacy", "check_6": "PASS", "check_7": "PASS", "evidence": "plans/active/tooling-sern-te-wake-grid.md §6 #6・#7"}
for label, rows, te, req, co, want in (
        ("(a) 修正後の系列: 同じ id の 2 設計点 (格子署名は違う) は両方採る", [row_({k: op_() for k in OPS2}, x) for x in (0.1, 0.2)], ON, ID_NEW, (), 2),
        ("(b) 修正後の系列: 識別の導入前の行 (格子変形なし) は外す", [row_({k: op_("legacy") for k in OPS2})], ON, ID_NEW, (), 0),
        ("(b) 修正後の系列: 格子変形なし (off) の新しい行は外す", [row_({k: op_("old") for k in OPS2})], ON, ID_NEW, (), 0),
        ("(b) 修正後の系列: 変換器の版違い (L_b 1.0 でも別 id) は外す", [row_({k: op_("conv2", te=ON) for k in OPS2})], ON, ID_NEW, (), 0),
        ("(b) 修正後の系列: 不明 (None) は外す", [row_({"m6_on": op_(), "m10_on": op_("none")})], ON, ID_NEW, (), 0),
        ("(b) 修正後の系列: 片方の作動点だけ旧方式でも外す", [row_({"m6_on": op_(), "m10_on": op_("old")})], ON, ID_NEW, (), 0),
        ("(b) 修正後の系列で要求の id が不明: 承認の無い行は採らない", [row_({k: op_() for k in OPS2})], ON, None, (), 0),
        ("(b) 修正前の系列: 識別の導入前の行は従来どおり採る", [row_({k: op_("legacy") for k in OPS2})], "off", None, (), 1),
        ("(b) 修正前の系列: off の新しい行は採る (変換器の版は照合しない)", [row_({k: op_("conv2", te="off") for k in OPS2})], "off", None, (), 1),
        ("(b) 修正前の系列: 格子変形を入れた評価は外す", [row_({k: op_() for k in OPS2})], "off", None, (), 0),
        ("(b) 修正前の系列: 格子変形の実効値が不明 (None) は外す", [row_({k: {**op_("old"), D.TE_WAKE_EFFECTIVE: None} for k in OPS2})], "off", None, (), 0),
        ("(d) 承認なし: 旧評価 (識別の導入前) は外す", [row_({"m6_on": op_("legacy"), "m10_on": op_()})], ON, ID_NEW, (), 0),
        ("(d) 承認 [m6_on, legacy, PASS/PASS]: m6_on の旧評価 + m10_on の新方式を採る", [row_({"m6_on": op_("legacy"), "m10_on": op_()})], ON, ID_NEW, (A6,), 1),
        ("(d) 同じ承認でも名指しの無い m10_on の旧評価は外す", [row_({"m6_on": op_("legacy"), "m10_on": op_("legacy")})], ON, ID_NEW, (A6,), 0),
        ("(d) legacy の承認は識別のある旧方式 (off の id) を持ち越さない", [row_({"m6_on": op_("old"), "m10_on": op_()})], ON, ID_NEW, (A6,), 0),
        ("(d) 承認 [m6_on, id 指定]: その id の旧評価だけ採る", [row_({"m6_on": op_("old"), "m10_on": op_()})], ON, ID_NEW, ({**A6, "from_eval_method": ID_OLD},), 1),
        ("(d) id 指定の承認は別の id (変換器の版違い) を持ち越さない", [row_({"m6_on": op_("conv2", te="off"), "m10_on": op_()})], ON, ID_NEW, ({**A6, "from_eval_method": ID_OLD},), 0),
        ("(d) 承認は不明 (None) を持ち越さない", [row_({"m6_on": op_("none"), "m10_on": op_()})], ON, ID_NEW, (A6,), 0),
        ("(d) 承認は処置以外の採否を緩めない (scalarGradient gg)", [row_({"m6_on": op_("legacy", sg="gg"), "m10_on": op_()})], ON, ID_NEW, (A6,), 0),
        ("(d) 識別が一致する評価に承認は要らない (承認があっても同じ)", [row_({k: op_() for k in OPS2})], ON, ID_NEW, (A6,), 1)):
    got = n_learn(rows, te, req, co)
    check(f"{label}: 採用 {got} 件 (期待 {want})", got == want)
check("(b) 既存の試験の偽物 (評価方式の属性なし) は要求なしとして動く", D._method_policy(object()) == ("off", None, ())
      and D._learnable(row_({k: op_("legacy") for k in OPS2})))

# 承認の検査 (PASS でない承認は効かない・不正な記録は止まる)
ITEMS_OK = [A6, {**A6, "op": "m10_on", "check_6": "FAIL"}, {**A6, "op": "m10_on", "check_7": "UNDECIDABLE", "from_eval_method": ID_OLD}]
allit, appr = D.parse_eval_method_carry_over(ITEMS_OK, set(OPS2), ON)
check("(d) 判定が両方 PASS の項目だけが承認になる (FAIL・UNDECIDABLE は記録だけ)", len(allit) == 3 and [a["op"] for a in appr] == ["m6_on"])
check("(d) PASS でない承認は効かない", n_learn([row_({"m6_on": op_(), "m10_on": op_("legacy")})], ON, ID_NEW, appr) == 0)
for label, items, te in (("無い作動点", [{**A6, "op": "m4_off"}], ON), ("id の形でない", [{**A6, "from_eval_method": "abc"}], ON),
                         ("判定の値", [{**A6, "check_6": "OK"}], ON), ("根拠が空", [{**A6, "evidence": " "}], ON),
                         ("未知キー", [{**A6, "note": "x"}], ON), ("キーの欠け", [{k: v for k, v in A6.items() if k != "check_7"}], ON),
                         ("リストでない", A6, ON), ("修正前の系列 (off) で承認を書く", [A6], "off")):
    check(f"(d) 不正な承認の記録は止まる: {label}", raises(lambda: D.parse_eval_method_carry_over(items, set(OPS2), te)))

# 実際のキャンペーン: 基準 YAML の L_b で要求が決まる (driver は 2D runner)
BASE2 = yaml.safe_load(open(ROOT / "case" / "46.sern_design" / "problem_moo_sst_node_3op.yaml"))


def campaign(name, lb=None, co=None):
    d = tmp / name; d.mkdir()
    raw = json.loads(json.dumps(BASE2))
    if lb is not None:
        raw["mesh"]["te_wake_blend_H"] = lb
    if co is not None:
        raw.setdefault("opt", {})[D.EVAL_METHOD_CARRY_OVER_KEY] = co
    (d / "problem.yaml").write_text(yaml.safe_dump(raw, allow_unicode=True))
    return D.SernCampaign(d / "problem.yaml", d / "camp")


c0 = campaign("c_off")
check("(b) 基準 YAML にキーが無いキャンペーン: 要求なし (off、id は照合しない)", c0.te_wake_required == "off" and c0.eval_method_required is None)
A6c = {**A6, "op": "cruise"}
c1 = campaign("c_on", 1.0, [A6c])
check("(b) L_b 1.0 のキャンペーン: 格子変形と評価方式の id (変換器があれば 64 桁) を要求",
      c1.te_wake_required == "L_b=1.0;curve=None" and (c1.eval_method_required is None) == (not HAVE_CONV)
      and (not HAVE_CONV or len(c1.eval_method_required) == 64))
check("(d) キャンペーンが承認の記録を持ち、要約に残す", [a["op"] for a in c1.eval_method_carry_over] == ["cruise"]
      and c1.summary()[D.EVAL_METHOD_CARRY_OVER_KEY] == [A6c] and c1.summary()["required_te_wake"] == "L_b=1.0;curve=None")
check("(d) 修正前の系列 (キー無し) で承認を書くとキャンペーンの作成で止まる", raises(lambda: campaign("c_bad", None, [A6c])))
check("(d) 基準 YAML に無い作動点の承認はキャンペーンの作成で止まる", raises(lambda: campaign("c_bad2", 1.0, [A6])))
_led = tmp / "c_resume"; _led.mkdir(); (_led / "camp").mkdir()
(_led / "camp" / "ledger.jsonl").write_text("\n".join(json.dumps(dict(r, tag=f"t{i}")) for i, r in enumerate(
    [row_({k: op_("legacy") for k in ("cruise", "accel", "lownpr")})])) + "\n")
_raw = json.loads(json.dumps(BASE2)); _raw["mesh"]["te_wake_blend_H"] = 1.0
(_led / "problem.yaml").write_text(yaml.safe_dump(_raw, allow_unicode=True))
check("(b) 識別の導入前の台帳の上で L_b 1.0 のキャンペーンを再開: 旧評価は学習に入らない (台帳には残る)",
      D.SernCampaign(_led / "problem.yaml", _led / "camp").summary()["n_pass"] == 0
      and len((_led / "camp" / "ledger.jsonl").read_text().splitlines()) == 1)
(_led / "problem.yaml").write_text(yaml.safe_dump({**_raw, "mesh": {k: v for k, v in _raw["mesh"].items() if k != "te_wake_blend_H"}}, allow_unicode=True))
check("(b) 同じ台帳をキー無しのキャンペーンで読むと従来どおり採る", D.SernCampaign(_led / "problem.yaml", _led / "camp").summary()["n_pass"] == 1)

shutil.rmtree(tmp, ignore_errors=True)
if skipped:
    print("省略: " + "; ".join(skipped))
print("VERDICT:", ("PASS" if not skipped else "PASS (一部省略)") if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
