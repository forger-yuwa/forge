#!/usr/bin/env python3
"""厚さ 0 の板の自由端の速度再構成の処置: 設計チェーンへの配線の試験
(plan convection-zero-thickness-edge-reconstruction §4「設計チェーンへの配線」・§5 の 2、codex plan レビュー 2026-10-08 M3、
codex diagnose 2026-10-08 edge-weight-preprocessing)。forge は回さない。

合格 (測る前に固定):
  (a) 未指定で生成される config が変更前 (commit cac3c9a4) とバイト一致 (2D/3D の本段・品質確認用 cell)
  (b) 指定で space に `zeroThicknessEdgeVelocity: {field: w_recon_vel}` が出る (差はその 1 か所だけ)・品質確認用 cell には出ない・
      段階起動の全段に同じ値が出て、各段がソルバの起動条件 (node・SLAU・非周期・非軸対称・w あり) を満たす・不正な指定は止まる
  (c) stage_key: OFF→ON・rings 2/3・w の変更・格子の変更・フィールド名の変更で区間が分かれる、明示の無効はキー無しと同じ、
      キー無しの段は変更前と同じ値、ハッシュは属性から転記せず再計算する
  (d) 来歴 (metrics) に有効状態・タグ・rings・格子署名・w のハッシュが残り、起動の記録がそろわない有効は不明 (None)
  (e) 設計 DB: 処置ありと無しの評価が別扱い、処置なしの旧評価の持ち越しは作動点を名指しした明示の指定があるときだけ
  (f) 同一格子 restart (restart_field / restart_by_index) は宛先の /AUX を保持して元の /AUX を写さない、
      cross-mesh (interp_field) は元の /AUX を移植しない
前処理の道具 (mark_zero_thickness_edges.py) は未実装なので、その mesh_signature / field_hash は偽物を差し込む (配線の試験)。
"""
import hashlib
import importlib.util
import json
import os
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

# 前処理の道具の偽物 (格子署名 = MESH/COORD のハッシュ、w のハッシュ = 値のハッシュ、どちらも 16 進 64 桁)。
# sys.modules が sys.path より優先されるので、実際の道具があってもこちらが使われる (実際の道具は (g) で試す)
_fake = types.ModuleType("mark_zero_thickness_edges")
_fake.mesh_signature = lambda h5path: hashlib.sha256(
    b"sig" + np.ascontiguousarray(h5py.File(h5path, "r")["MESH/COORD"][...]).tobytes()).hexdigest()
_fake.field_hash = lambda w: hashlib.sha256(np.ascontiguousarray(np.asarray(w, "<f4")).tobytes()).hexdigest()
sys.modules["mark_zero_thickness_edges"] = _fake

import stage_manifest as SM  # noqa: E402
from forge_design.evaluate import runner_sern as RS  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
from forge_design.probdef import Problem  # noqa: E402

REF_COMMIT = "cac3c9a4"     # 配線を入れる前の HEAD
fails = 0
skipped = []
tmp = Path(tempfile.mkdtemp(prefix="zte_wiring_"))


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


def load_ref(rel, modname, package=None):
    """変更前のソースを git から読んで別名のモジュールにする (相対 import は package で解決)。読めなければ None。"""
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{REF_COMMIT}:{rel}"], capture_output=True, text=True)
    if r.returncode != 0:
        return None
    # 旧ソースは `os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3])` で既定値を先に評価するので、
    # 写しは親が 3 段以上ある場所に置く (/tmp 直下の一時ディレクトリでは parents[3] が IndexError になる)
    d = tmp / "ref" / "forge_design" / "evaluate"
    d.mkdir(parents=True, exist_ok=True)
    f = d / (modname.replace(".", "_") + ".py")
    f.write_text(r.stdout)
    spec = importlib.util.spec_from_file_location(modname, f)
    m = importlib.util.module_from_spec(spec)
    if package:
        m.__package__ = package
    saved = os.environ.get("FORGE_ROOT")
    os.environ["FORGE_ROOT"] = str(ROOT)          # 一時ファイルの位置からリポジトリを導けないので明示する
    try:
        spec.loader.exec_module(m)
    finally:
        if saved is None:
            os.environ.pop("FORGE_ROOT", None)
        else:
            os.environ["FORGE_ROOT"] = saved
    return m


def prob(evaluate=None, mesh=None):
    return Problem(name="t", type="sern_2d", gamma=1.4, cp=1004.5, spec={}, dv={}, geometry={},
                   mesh=dict(mesh or {}), evaluate=dict(evaluate or {}), raw={"gas": {}})


ZTE = {"tags": ["cowl_out", "cowl_in"], "rings": 2}
COORD = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.5, 0.5, 0.0]])


def make_h5(path, coords=COORD, w="default", tags="cowl_in,cowl_out", rings=2, val=1.0, sha_attr="deadbeef"):
    """合成の入力 h5: MESH/COORD・VALUE (保存量)・/AUX/w_recon_vel (属性 tags/rings と偽の sha256)。w=None なら /AUX 無し。"""
    n = len(coords)
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=np.asarray(coords, float).ravel())
        for k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"):
            f.create_dataset(f"VALUE/{k}", data=np.full(n, val * (2.0 if k == "roe" else 1.0), np.float32))
        if w is not None:
            ww = np.ones(n) if isinstance(w, str) else np.asarray(w, float)
            if isinstance(w, str):
                ww[-1] = 0.0
            ds = f.create_dataset(f"AUX/{RS.ZTE_FIELD}", data=ww)
            ds.attrs["tags"] = tags; ds.attrs["rings"] = rings; ds.attrs["sha256"] = sha_attr
            ds.attrs["generator"] = "fake"; ds.attrs["n_edge"] = 3; ds.attrs["n_zero"] = int((ww == 0).sum())


def run_dir_with(name, cfg, bcond="ramp: {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: }\n", **h5kw):
    d = tmp / name; d.mkdir()
    (d / "solverConfig.yaml").write_text(cfg); (d / "solverConfig_main.yaml").write_text(cfg)
    (d / "bcondConfig.yaml").write_text(bcond)
    make_h5(d / "sern.h5", **h5kw)
    return d


# ======================================================================== (a) 未指定は変更前とバイト一致
REF2 = load_ref("design/forge_design/evaluate/runner_sern.py", "forge_design.evaluate._runner_sern_ref", "forge_design.evaluate")
VARIANTS = [
    ("euler cell (既定)", {}, {}),
    ("euler node", {}, {"discretization": "node"}),
    ("sst node 生産風", {"model": "sst", "limiter_ref": {"length": 1.0, "ro": 0.1, "p": 1e4, "a": 500.0}, "implicit_relax": 0.7,
                       "p_min": 1.0, "wall_treatment_sst": 0}, {"discretization": "node", "node_inlet_corner_wall": 1}),
    ("sst node chi 0・旧リミッタ", {"model": "sst", "limiter_scaled": 0, "limiter": 1}, {"discretization": "node", "slau_wall_normal_chi": 0}),
    ("sst cell", {"model": "sst"}, {"discretization": "cell"}),
]
if REF2 is None:
    skipped.append("(a) git に変更前の版が無い")
    print("SKIP (a): git show が使えない")
else:
    for name, ev, ms in VARIANTS:
        p = prob(ev, ms)
        new, old = RS._solver_config(p, 4000, 500, 2.5, 2851.0), REF2._solver_config(p, 4000, 500, 2.5, 2851.0)
        check(f"(a) 2D {name}: 本段 config がバイト一致", new == old)
        check(f"(a) 3D {name}: 本段 config がバイト一致", R3._solver_config(p, 4000, 500, 2.5, 2851.0) == old)
        disc = p.mesh.get("discretization", "cell")
        legacy_qc = (old.replace(f'discretization: "{disc}"', 'discretization: "cell"')
                     .replace(", nodeWallDirichlet: 1", "").replace(", nodeInletCornerWall: 1", ""))
        check(f"(a) {name}: 品質確認用 cell の config がバイト一致", RS.qc_cell_config(new, disc) == legacy_qc)

# ======================================================================== (b) 指定で space にキーが出る
p_on = prob({"model": "sst", "zero_thickness_edge_velocity": ZTE}, {"discretization": "node"})
p_off = prob({"model": "sst"}, {"discretization": "node"})
cfg_on, cfg_off = RS._solver_config(p_on, 4000, 500, 2.5, 2851.0), RS._solver_config(p_off, 4000, 500, 2.5, 2851.0)
check("(b) space に {field: w_recon_vel}", (yaml.safe_load(cfg_on)["space"] or {}).get("zeroThicknessEdgeVelocity") == {"field": "w_recon_vel"})
check("(b) 無効との差はキーの 1 か所だけ", cfg_on.count(RS.ZTE_SPACE_FRAGMENT) == 1 and cfg_on.replace(RS.ZTE_SPACE_FRAGMENT, "") == cfg_off)
check("(b) 3D も同じ", R3._solver_config(p_on, 4000, 500, 2.5, 2851.0) == cfg_on)
check("(b) 品質確認用 cell の config には出ない", RS.qc_cell_config(cfg_on, "node") == RS.qc_cell_config(cfg_off, "node")
      and "zeroThickness" not in RS.qc_cell_config(cfg_on, "node"))
check("(b) タグは並べ替え", RS.zte_spec(p_on) == {"tags": ["cowl_in", "cowl_out"], "rings": 2})
check("(b) 識別子", RS.zte_signature_of(RS.zte_spec(p_on)) == "tags=cowl_in,cowl_out;rings=2" and RS.zte_signature_of(None) == "off")
check("(b) cell + 指定は止まる", raises(lambda: RS._solver_config(prob({"zero_thickness_edge_velocity": ZTE}, {"discretization": "cell"}), 1, 1, 1.0, 1.0)))
T2 = ["cowl_in", "cowl_out"]
for bad_name, bad in (("tags 空", {"tags": [], "rings": 2}), ("tags が文字列", {"tags": "cowl_in", "rings": 2}),
                      ("tags 1 つ", {"tags": ["cowl_in"], "rings": 2}), ("tags 3 つ", {"tags": T2 + ["ramp"], "rings": 2}),
                      ("tags 重複", {"tags": ["cowl_in", "cowl_in"], "rings": 2}), ("rings 0", {"tags": T2, "rings": 0}),
                      ("rings 負", {"tags": T2, "rings": -1}), ("rings 真偽値", {"tags": T2, "rings": True}),
                      ("rings 小数", {"tags": T2, "rings": 2.5}), ("未知キー", {"tags": T2, "rings": 2, "field": "x"}),
                      ("タグに空白", {"tags": ["cowl in", "cowl_out"], "rings": 2}), ("map でない", True)):
    check(f"(b) 不正な指定は止まる: {bad_name}", raises(lambda: RS.zte_spec(prob({"zero_thickness_edge_velocity": bad}))))
from forge_design.meshing.mesh_sern import PHYS_SERN  # noqa: E402
from forge_design.meshing.mesh_sern3d import PHYS_SERN3D  # noqa: E402
check("(b) 格子に無いタグは止まる (2D/3D)", raises(lambda: RS.zte_spec(prob({"zero_thickness_edge_velocity": {"tags": ["cowl_x", "cowl_in"], "rings": 2}}), PHYS_SERN.keys()))
      and raises(lambda: RS.zte_spec(prob({"zero_thickness_edge_velocity": {"tags": ["cowl_x", "cowl_in"], "rings": 2}}), PHYS_SERN3D.keys())))
check("(b) 生産のタグは 2D/3D とも通る", RS.zte_spec(p_on, PHYS_SERN.keys()) is not None and RS.zte_spec(p_on, PHYS_SERN3D.keys()) is not None)
check("(b) null / false は無効", RS.zte_spec(prob({"zero_thickness_edge_velocity": None})) is None
      and RS.zte_spec(prob({"zero_thickness_edge_velocity": False})) is None)


# 段階起動の全段 (層流暖機 ramp → soft ramp → mid → 本段、warm start の適応段) の config を捕まえる
def staged(cfg, name, **kw):
    rd = run_dir_with(name, cfg)
    seen, k = [], [0]

    def fake_run(cmd, **_kw):
        d = Path(cmd[1]); seen.append((d / "solverConfig.yaml").read_text())
        k[0] += 1; (d / f"res_{k[0]}.h5").write_bytes(b"")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")
    saved = (RS.subprocess, RS.restart_by_index, RS.warm_from_run)
    RS.subprocess = types.SimpleNamespace(run=fake_run)
    RS.restart_by_index = lambda *a, **k_: None
    RS.warm_from_run = lambda *a, **k_: {}
    try:
        rc = RS.run_staged(rd, "full", **kw)
    finally:
        RS.subprocess, RS.restart_by_index, RS.warm_from_run = saved
    return rd, rc, seen


ZKEY = lambda k: {a: b for a, b in k.items() if a.startswith(SM.ZTE_PREFIX)}
for label, kw, n_stage in (("層流暖機 ramp + soft ramp + mid", dict(soft_steps=20, warm_lam_steps=20, warm_lam_ramp=[0.1, 0.2], soft_ramp=[0.3, 0.5], mid_steps=20), 7),
                           ("warm start (適応段) + mid", dict(soft_steps=20, mid_steps=20, warm_src=tmp), 3)):
    rd, rc, seen = staged(cfg_on, "staged_" + str(len(label)), **kw)
    bc = (rd / "bcondConfig.yaml").read_text()
    kinds = ["1 次" if "convMethod: 0" in c else "2 次" for c in seen]
    check(f"(b) {label}: {len(seen)} 段 ({','.join(kinds)}) すべてに同じ値", len(seen) == n_stage and all(c.count(RS.ZTE_SPACE_FRAGMENT) == 1 for c in seen))
    ok = True
    for c in seen:
        try:
            RS.check_zte_stage(c, bc, rd)
        except ValueError as e:
            ok = False; print("    ", e)
    check(f"(b) {label}: 全段がソルバの起動条件 (node・SLAU・非周期・非軸対称・w あり) を満たす", ok)
    zk = [json.dumps(ZKEY(SM.stage_key(c, bc, rd)), sort_keys=True) for c in seen]
    check(f"(b) {label}: 全段の stage_key の処置の部分が同じ", len(set(zk)) == 1 and '"1"' in zk[0])
    _zl = rd / RS.ZTE_LAUNCHES
    recs = [json.loads(l) for l in _zl.read_text().splitlines()] if _zl.exists() else []
    check(f"(b) {label}: 起動前の記録が段ごとに 1 行 (cfg_fnv が各段の config と対応、識別量は再計算値)",
          [r["cfg_fnv"] for r in recs] == [SM.fnv1a64(c) for c in seen]
          and {(r["mesh_signature"], r["field_hash"]) for r in recs} == {(_fake.mesh_signature(rd / "sern.h5"),
                                                                          _fake.field_hash(h5py.File(rd / "sern.h5", "r")["AUX/w_recon_vel"][...]))})
rd, rc, seen = staged(cfg_off, "staged_off", soft_steps=20, warm_lam_steps=20, mid_steps=20)
check("(b) 無効の段階起動: キーは出ず、起動前の記録も作らない", all("zeroThickness" not in c for c in seen) and not (rd / RS.ZTE_LAUNCHES).exists())

# 起動直前の検査: 条件違反は forge を起動する前に ValueError
for label, cfg_bad, bcond in (("cell", cfg_on.replace('discretization: "node"', 'discretization: "cell"'), None),
                              ("ROE", cfg_on.replace('solver: "SLAU"', 'solver: "ROE"'), None),
                              ("軸対称", cfg_on.replace("isAxisymmetric: 0", "isAxisymmetric: 1"), None),
                              ("gpu 0", cfg_on.replace("gpu: 1", "gpu: 0"), None),
                              ("space の値が null", cfg_on.replace("{field: w_recon_vel}", "null"), None),
                              ("周期境界", cfg_on, "side: {physID: 9, kind: periodic, outputHDFflg: 0, ints: , floats: }\n")):
    d = run_dir_with("bad_" + label, cfg_bad, **({"bcond": bcond} if bcond else {}))
    called = []
    saved = RS.subprocess
    RS.subprocess = types.SimpleNamespace(run=lambda *a, **k: called.append(1))
    try:
        r = raises(lambda: RS.run_forge(d))
    finally:
        RS.subprocess = saved
    check(f"(b) 起動直前の検査で止まる: {label} (forge は起動しない)", r and not called)
d = run_dir_with("bad_no_w", cfg_on, w=None)
check("(b) 起動直前の検査で止まる: meshFileName に w が無い", raises(lambda: RS._record_zte_launch(d)))
d = run_dir_with("ok_slau2", cfg_on.replace('solver: "SLAU"', 'solver: "SLAU2"'))
check("(b) SLAU2 は通る (ソルバの契約と同じ)", RS.check_zte_stage((d / "solverConfig.yaml").read_text(), "", d) is not None)

# prepare が使う照合: w を書く関数の差し替えと要求との一致
d = run_dir_with("verify", cfg_on)
check("(b) verify_zte_field: 要求と一致すれば識別量を返す", len(RS.verify_zte_field(d, RS.zte_spec(p_on))["field_hash"]) == 64)
check("(b) verify_zte_field: rings が要求と違えば止まる", raises(lambda: RS.verify_zte_field(d, {"tags": ["cowl_in", "cowl_out"], "rings": 3})))
d = run_dir_with("verify_len", cfg_on, w=np.ones(3))
check("(b) verify_zte_field: w の長さ ≠ 節点数 (別格子の w) は止まる", raises(lambda: RS.verify_zte_field(d, RS.zte_spec(p_on))))
_saved_tool_path = RS.ZTE_TOOL
RS.ZTE_TOOL = tmp / "no_such_tool.py"
try:
    check("(b) mark_zte_field: 道具が無ければ止まる", raises(lambda: RS.mark_zte_field(d, RS.zte_spec(p_on)), RuntimeError))
finally:
    RS.ZTE_TOOL = _saved_tool_path

# ======================================================================== (c) stage_key
REFSM = load_ref("solver_density_cuda/tools/stage_manifest.py", "_stage_manifest_ref")
BC = "ramp: {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: }\n"
d_on = run_dir_with("sk_on", cfg_on)
CFGS_OFF = [cfg_off, RS._solver_config(prob({}, {}), 1, 1, 1.0, 1.0),
            "turbulence:\n  model: sst\nspace:\n  convMethod: 1\n  limiter: 2\nmesh:\n  scalarGradient: lsq\n",
            "conjugate: {mode: 1, flux: q}\nspace: {convMethod: 1}\n", "space: {convMethod: 1\n"]
if REFSM is None:
    skipped.append("(c) git に変更前の stage_manifest が無い")
else:
    for i, c in enumerate(CFGS_OFF):
        for rdir in (None, d_on):
            check(f"(c) キー無しの段 #{i} (run_dir={'あり' if rdir else 'None'}): stage_key が変更前と同じ",
                  SM.stage_key(c, BC, rdir) == REFSM.stage_key(c, BC, rdir))
check("(c) キー無しの段に処置のキーは足さない", not ZKEY(SM.stage_key(cfg_off, BC, d_on)))
k_off, k_on = SM.stage_key(cfg_off, BC, d_on), SM.stage_key(cfg_on, BC, d_on)
check("(c) OFF → ON で別キー", k_off != k_on)
seg = SM.segments({"stages": [{"tag": t, "key": k} for t, k in (("off1", k_off), ("off2", k_off), ("on1", k_on), ("on2", k_on))]})
check("(c) OFF,OFF,ON,ON → 2 区間 (ON 同士は連結)", [[s["tag"] for s in sg] for sg in seg] == [["off1", "off2"], ["on1", "on2"]])
_with = lambda v: cfg_off.replace("venkatK: 0.05}", "venkatK: 0.05, zeroThicknessEdgeVelocity: %s}" % v)
for label, v in (("enabled: 0", "{enabled: 0}"), ("enabled: false", "{enabled: false}"), ("enabled: 0 + field", "{enabled: 0, field: w_recon_vel}")):
    check(f"(c) 明示の無効 ({label}) はキー無しと同じ", _with(v) != cfg_off and SM.stage_key(_with(v), BC, d_on) == k_off)
check("(c) {enabled: 1} は {field: w_recon_vel} と同じキー (field の既定)", SM.stage_key(_with("{enabled: 1}"), BC, d_on) == SM.stage_key(cfg_on, BC, d_on))
for label, v in (("null", "null"), ("false", "false"), ("未知キー", "{enabled: 1, rings: 2}"), ("enabled 2", "{enabled: 2}")):
    k_bad = SM.stage_key(_with(v), BC, d_on)
    check(f"(c) ソルバが拒否する値 ({label}) は無効と連結しない (invalid)", k_bad != k_off and str(k_bad.get(SM.ZTE_PREFIX + ".field")).startswith("invalid"))
d_r3 = run_dir_with("sk_r3", cfg_on, rings=3)
check("(c) rings 2 → 3 で別キー", SM.stage_key(cfg_on, BC, d_on) != SM.stage_key(cfg_on, BC, d_r3))
d_w = run_dir_with("sk_w", cfg_on, w=np.array([1.0, 1.0, 1.0, 0.0, 0.0]))
check("(c) w の変更 (同じ属性) で別キー", SM.stage_key(cfg_on, BC, d_on) != SM.stage_key(cfg_on, BC, d_w))
d_m = run_dir_with("sk_mesh", cfg_on, coords=COORD * 2.0)
check("(c) 格子の変更 (同じ w・属性) で別キー", SM.stage_key(cfg_on, BC, d_on) != SM.stage_key(cfg_on, BC, d_m))
d_t = run_dir_with("sk_tags", cfg_on, tags="cowl_in")
check("(c) タグの変更で別キー", SM.stage_key(cfg_on, BC, d_on) != SM.stage_key(cfg_on, BC, d_t))
check("(c) フィールド名の変更で別キー", SM.stage_key(cfg_on, BC, d_on) != SM.stage_key(cfg_on.replace("field: w_recon_vel", "field: w_other"), BC, d_on))
d_sha = run_dir_with("sk_sha", cfg_on, sha_attr="0123456789abcdef")
check("(c) 属性の sha256 だけ違っても同じキー (ハッシュは転記せず再計算)", SM.stage_key(cfg_on, BC, d_on) == SM.stage_key(cfg_on, BC, d_sha))
kk = ZKEY(SM.stage_key(cfg_on, BC, d_on))
check("(c) 有効の段のキー = 有効・フィールド名・タグ (並べ替え)・rings・格子署名・w のハッシュ",
      kk == {SM.ZTE_PREFIX + ".enabled": "1", SM.ZTE_PREFIX + ".field": "w_recon_vel", SM.ZTE_PREFIX + ".tags": "cowl_in,cowl_out",
             SM.ZTE_PREFIX + ".rings": "2", SM.ZTE_PREFIX + ".mesh_signature": _fake.mesh_signature(d_on / "sern.h5"),
             SM.ZTE_PREFIX + ".field_hash": _fake.field_hash(h5py.File(d_on / "sern.h5", "r")["AUX/w_recon_vel"][...])}, kk)
d_vf = tmp / "sk_valuefile"; d_vf.mkdir(); make_h5(d_vf / "sern.h5", w=None); make_h5(d_vf / "other.h5")
cfg_vf = cfg_on.replace('valueFileName: "sern.h5"', 'valueFileName: "other.h5"')
check("(c) w は meshFileName からだけ読む (valueFileName にあっても missing)",
      ZKEY(SM.stage_key(cfg_vf, BC, d_vf)).get(SM.ZTE_PREFIX + ".field_hash") == "missing")
check("(c) run_dir 無し → unknown (OFF と連結しない)", ZKEY(SM.stage_key(cfg_on, BC, None)).get(SM.ZTE_PREFIX + ".field_hash") == "unknown"
      and SM.stage_key(cfg_on, BC, None) != SM.stage_key(cfg_off, BC, None))
sm_man = SM.StageManifest(d_on); sm_man.add("mid", cfg_on, BC); sm_man.add("main", cfg_on, BC); sm_man.write()
check("(c) StageManifest も同じキーで連結", len(SM.segments(json.load(open(d_on / "stage_manifest.json")))) == 1)
_saved_mod = SM.ZTE_TOOL_MODULE
SM.ZTE_TOOL_MODULE = "no_such_zte_tool_module"
try:
    _stop = False
    try:
        SM.stage_key(cfg_on, BC, d_on)
    except SystemExit:
        _stop = True
    check("(c) 道具を import できなければ止める (識別できないまま連結しない)", _stop)
    check("(c) 同じとき runner の起動直前の検査は ValueError (driver で ERROR に分類、プロセスは止めない)",
          raises(lambda: RS.check_zte_stage(cfg_on, BC, d_on)))
    _pv = RS.zte_provenance(d_on)
    check("(c) 同じとき来歴は不明 (collect を止めない)", _pv["effective"] is None and "import" in _pv["reason"])
finally:
    SM.ZTE_TOOL_MODULE = _saved_mod

# ======================================================================== (d) 来歴
def launch(d, cfg, zte=None, flat=False):
    rec = {"time": 0, "cfg_fnv": SM.fnv1a64(cfg), "slauWallNormalChi": 1, "scalarGradient": "lsq"}
    if zte is not None:
        if flat:
            rec.update({"zeroThicknessEdgeVelocity": 1, **{f"zeroThicknessEdgeVelocity_{k}": v for k, v in zte.items()}})
        else:
            rec["zeroThicknessEdgeVelocity"] = zte
    with open(d / "forge_launches.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def ident(d):
    a = SM.zte_identity(d / "sern.h5", "w_recon_vel")
    return {"mesh_signature": a["mesh_signature"], "field_hash": a["field_hash"]}


d = run_dir_with("prov_off", cfg_off); launch(d, cfg_off)
pv = RS.zte_provenance(d)
check("(d) 無効の run → effective off (有効状態 False)", pv["effective"] == "off" and pv["enabled"] is False)
d = run_dir_with("prov_off_nolaunch", cfg_off)
check("(d) 無効・起動記録なし (旧 run) → off", RS.zte_provenance(d)["effective"] == "off")
d = run_dir_with("prov_off_but_on", cfg_off); launch(d, cfg_off, {"enabled": 1, "field": "w_recon_vel"})
check("(d) 設定は無効なのに起動が有効を記録 → 不明", RS.zte_provenance(d)["effective"] is None)


def prov_on(name, solver=True, flat=False, tamper=None, **h5kw):
    d = run_dir_with(name, cfg_on, **h5kw)
    RS._record_zte_launch(d)                      # runner の起動前の記録 (本段)
    z = {"enabled": 1, "field": "w_recon_vel", **ident(d), "n_E": 7, "n_S2": 41} if solver else None
    if tamper:
        z = tamper(z) if z else z
    launch(d, cfg_on, z, flat=flat)
    return d


d = prov_on("prov_on")
pv = RS.zte_provenance(d)
check("(d) 有効の run → effective = タグ・rings の識別子", pv["effective"] == "tags=cowl_in,cowl_out;rings=2", pv["reason"])
check("(d) 来歴にタグ・rings・格子署名・w のハッシュ・節点数・属性の数・ソルバの数",
      pv["tags"] == ["cowl_in", "cowl_out"] and pv["rings"] == 2 and pv["mesh_signature"] == ident(d)["mesh_signature"]
      and pv["field_hash"] == ident(d)["field_hash"] and pv["n_nodes"] == 5 and pv["attrs"]["n_edge"] == 3
      and pv["launch"]["n_S2"] == 41 and pv["pre_launch"]["field_hash"] == ident(d)["field_hash"])
check("(d) 来歴の w のハッシュは属性の sha256 を転記していない", pv["field_hash"] != "deadbeef" and pv["attrs"]["sha256"] == "deadbeef")
check("(d) 起動記録が平坦な形でも読む", RS.zte_provenance(prov_on("prov_flat", flat=True))["effective"] == "tags=cowl_in,cowl_out;rings=2")
check("(d) ソルバの起動記録に処置の記録が無い (旧バイナリ・未記録) → 不明", RS.zte_provenance(prov_on("prov_nosolver", solver=False))["effective"] is None)
check("(d) ソルバの記録した w のハッシュが再計算と違う → 不明",
      RS.zte_provenance(prov_on("prov_badhash", tamper=lambda z: {**z, "field_hash": "fh-other"}))["effective"] is None)
check("(d) ソルバが無効を記録 → 不明", RS.zte_provenance(prov_on("prov_soloff", tamper=lambda z: {**z, "enabled": 0}))["effective"] is None)
d = prov_on("prov_changed")
with h5py.File(d / "sern.h5", "r+") as f:
    f["AUX/w_recon_vel"][0] = 0.5
check("(d) 起動後に w が変わった → 不明", RS.zte_provenance(d)["effective"] is None)
d = prov_on("prov_notmain"); launch(d, cfg_on.replace("cfl: 2.5", "cfl: 0.5"), {"enabled": 1})
check("(d) 最後の起動が本段でない → 不明", RS.zte_provenance(d)["effective"] is None)
d = run_dir_with("prov_nopre", cfg_on); launch(d, cfg_on, {"enabled": 1, "field": "w_recon_vel"})
check("(d) runner の起動前の記録が無い (手で起動) → 不明", RS.zte_provenance(d)["effective"] is None)
d = prov_on("prov_bad_kw", tamper=lambda z: {**z, "enabled": "x"})
check("(d) 起動記録の壊れた値 → 不明 (例外にしない)", RS.zte_provenance(d)["effective"] is None)


# ソルバの起動記録に処置のキーが無いときは起動ログ forge_run.log の effective 行 (input/zeroThicknessEdge.cpp の書式) で確かめる
def prov_log(name, eff=1, fh=None, ms=None, full=True):
    d = prov_on(name, solver=False)
    i = ident(d); fh = fh or i["field_hash"]; ms = ms or i["mesh_signature"]
    lines = ["Init Thermo DB", f"'zeroThicknessEdgeVelocity' effective: {eff}"
             + (f" (field /AUX/w_recon_vel, w<1 nodes 1, field_sha256 {fh[:16]}, mesh_signature {ms[:16]})" if eff else " (off; bit-identical)")]
    if eff and full:
        lines += ["[zte] 生成元: mark_zero_thickness_edges.py 1 mode=edges tags=cowl_in,cowl_out tag_physids=cowl_in:5,cowl_out:6 rings=2 "
                  "n_edge_nodes(E)=3 n_marked_nodes(S)=1 n_nodes=5",
                  f"[zte] 検証: 値 [0,1]・有限、field_sha256 一致 ({fh})、mesh_signature 一致 ({ms}, zte-mesh-sig-v1)"]
    (d / "forge_run.log").write_text("\n".join(lines) + "\n")
    return RS.zte_provenance(d)


ON2_ = "tags=cowl_in,cowl_out;rings=2"
pv = prov_log("log_ok")
check("(d) 起動ログで有効を確認 (全桁のハッシュ・属性の数も拾う)", pv["effective"] == ON2_, pv["reason"])
check("(d) 起動ログの記録が来歴に残る", pv["launch"]["source"] == "forge_run.log" and pv["launch"]["n_marked_nodes"] == 1
      and pv["launch"]["n_edge_nodes"] == 3 and len(pv["launch"]["field_hash"]) == 64)
check("(d) 起動ログ (先頭 16 桁だけ) でも確認できる", prov_log("log_prefix", full=False)["effective"] == ON2_)
check("(d) 起動ログのハッシュが再計算と違う → 不明", prov_log("log_badhash", fh="0" * 64)["effective"] is None)
check("(d) 起動ログが effective 0 → 不明", prov_log("log_off", eff=0)["effective"] is None)


# collect が metrics.json に書く (2D/3D)。力の履歴とゲートは差し替える
def collect_case(mod, name, cfg, solver=True):
    d = prov_on(name) if cfg is cfg_on else run_dir_with(name, cfg)
    if cfg is not cfg_on:
        launch(d, cfg)
    prob_yaml = d / "problem.yaml"
    prob_yaml.write_text(yaml.safe_dump({"name": "t", "type": "sern_2d", "gas": {}, "spec": {}, "dv": {}, "geometry": {}, "mesh": {}, "evaluate": {}}))
    info = {"states": {"exhaust": {"ro": 1.0, "u": 1.0, "P": 1.0}, "ext": {"P": 1.0}}, "H_m": 1.0, "F_ideal_N_per_m": 1.0,
            "design": {"L_ramp": 1.0}, "moc_forces": {}, "discretization": "node", "half_W_m": 1.0}
    (d / "prepare_info.json").write_text(json.dumps(info))
    gates = lambda *a, **k: {"steadiness": {"series": {}}, "objective": "C_T", "verdict": "PASS"}
    saved = (RS.force_history, RS.evaluate_gates, R3.evaluate_gates)
    RS.force_history = lambda *a, **k: []; RS.evaluate_gates = gates; R3.evaluate_gates = gates
    try:
        mod.collect(prob_yaml, d, rc=0)
    finally:
        RS.force_history, RS.evaluate_gates, R3.evaluate_gates = saved
    return json.loads((d / "metrics.json").read_text())


for mod, dim in ((RS, "2D"), (R3, "3D")):
    m_on, m_off = collect_case(mod, f"col_on_{dim}", cfg_on), collect_case(mod, f"col_off_{dim}", cfg_off)
    check(f"(d) {dim} collect: metrics に処置の来歴と識別子",
          m_on[RS.ZTE_EFFECTIVE] == "tags=cowl_in,cowl_out;rings=2" and m_on["zero_thickness_edge_velocity"]["field_hash"]
          and m_off[RS.ZTE_EFFECTIVE] == "off" and m_off["zero_thickness_edge_velocity"]["enabled"] is False)

# ======================================================================== (e) 設計 DB
ON2, ON3 = "tags=cowl_in,cowl_out;rings=2", "tags=cowl_in,cowl_out;rings=3"


def row(effs, legacy=False):
    ops = {}
    for op, e in effs.items():
        o = {"scalar_gradient_effective": "lsq", "slau_wall_normal_chi_effective": 1}
        if not legacy:
            o[RS.ZTE_EFFECTIVE] = e
        ops[op] = o
    return {"status": "PASS", "flag_policy": D.FLAG_POLICY, "x": [0.0], "C_T_w": 0.9, "L_ramp": 1.0, "C_M_w": 0.0, "ops": ops}


class Fake:
    def __init__(self, rows, req="off", co=()):
        self.rows, self.zte_required, self.zte_carry_over, self.ref, self.ops = rows, req, tuple(co), (-0.9, 20.0), []


def n_learn(rows, req="off", co=()):
    return len(D.SernCampaign._XF(Fake(rows, req, co))[0])


OPS2 = ("m6_on", "m10_on")
for label, rows, req, co, want in (
        ("処置なしのキャンペーン: 旧行 (キー無し)", [row(dict.fromkeys(OPS2, "off"), legacy=True)], "off", (), 1),
        ("処置なしのキャンペーン: 処置なしの行", [row(dict.fromkeys(OPS2, "off"))], "off", (), 1),
        ("処置なしのキャンペーン: 処置ありの行は除外", [row(dict.fromkeys(OPS2, ON2))], "off", (), 0),
        ("処置なしのキャンペーン: 片方だけ処置あり", [row({"m6_on": "off", "m10_on": ON2})], "off", (), 0),
        ("不明 (None) は除外", [row(dict.fromkeys(OPS2, None))], "off", (), 0),
        ("処置ありのキャンペーン: 処置ありの行", [row(dict.fromkeys(OPS2, ON2))], ON2, (), 1),
        ("処置ありのキャンペーン: 旧行 (処置なし) は持ち越さない", [row(dict.fromkeys(OPS2, "off"), legacy=True)], ON2, (), 0),
        ("処置ありのキャンペーン: 処置なしの行は持ち越さない", [row(dict.fromkeys(OPS2, "off"))], ON2, (), 0),
        ("処置ありのキャンペーン: rings 3 の行は除外", [row(dict.fromkeys(OPS2, ON3))], ON2, (), 0),
        ("処置ありのキャンペーン: 不明は除外", [row({"m6_on": ON2, "m10_on": None})], ON2, ("m10_on",), 0),
        ("明示の持ち越し [m6_on]: m6_on 処置なし + m10_on 処置あり → 採用", [row({"m6_on": "off", "m10_on": ON2})], ON2, ("m6_on",), 1),
        ("明示の持ち越し [m6_on]: 両方処置なし → m10_on は名指し無しで除外", [row(dict.fromkeys(OPS2, "off"))], ON2, ("m6_on",), 0),
        ("明示の持ち越し [m6_on, m10_on]: 両方処置なし (旧行) → 採用", [row(dict.fromkeys(OPS2, "off"), legacy=True)], ON2, OPS2, 1),
        ("明示の持ち越しでも rings 3 は採用しない", [row({"m6_on": ON3, "m10_on": ON2})], ON2, ("m6_on",), 0)):
    got = n_learn(rows, req, co)
    check(f"(e) {label}: 採用 {got} 件 (期待 {want})", got == want)
s_ = D.SernCampaign.summary(Fake([dict(row(dict.fromkeys(OPS2, ON2)), tag="a"), dict(row(dict.fromkeys(OPS2, "off")), tag="b")], ON2))
check("(e) Pareto の母集団も同じ選別・要約に要求と作動点ごとの実効値",
      s_["n_pass"] == 1 and [p_["tag"] for p_ in s_["pareto"]] == ["a"] and s_["required_zero_thickness_edge_velocity"] == ON2
      and s_["pareto"][0]["ops"]["m6_on"][RS.ZTE_EFFECTIVE] == ON2)
check("(e) 既存の試験の偽物 (属性なし) は処置なしとして動く", n_learn([row(dict.fromkeys(OPS2, "off"), legacy=True)]) == 1)

# キャンペーン: 基準 YAML の要求・持ち越し指定の検査・評価行に識別子が残る
BASE = {"type": "sern_2d", "dv": {k: {"min": 0.0, "max": 1.0} for k in D.DV_ORDER},
        "spec": {"inflow": {"M_in": 2.5, "p_in": 1e5}, "external": {"p_inf": 5e3}, "H_m": 0.1,
                 "operating_points": [{"name": "m6_on", "weight": 0.6, "external": {}}, {"name": "m10_on", "weight": 0.4, "external": {}}]},
        "gas": {"gamma": 1.4, "cp": 1004.5}, "geometry": {"L_ramp_max": 12.0}, "evaluate": {}, "opt": {}}


def campaign(name, zte=None, co=None):
    d = tmp / name; d.mkdir()
    raw = json.loads(json.dumps(BASE))
    if zte is not None:
        raw["evaluate"]["zero_thickness_edge_velocity"] = zte
    if co is not None:
        raw["opt"][D.ZTE_CARRY_OVER_KEY] = co
    (d / "problem.yaml").write_text(yaml.safe_dump(raw))
    return D.SernCampaign(d / "problem.yaml", d / "camp")


check("(e) キャンペーンの要求: 未指定 → off、指定 → 識別子", campaign("c_off").zte_required == "off" and campaign("c_on", ZTE).zte_required == ON2)
check("(e) 持ち越し指定: 無い作動点は止まる", raises(lambda: campaign("c_bad1", ZTE, ["m4_off"])))
check("(e) 持ち越し指定: 処置なしのキャンペーンでは止まる", raises(lambda: campaign("c_bad2", None, ["m6_on"])))
check("(e) 持ち越し指定: リストでなければ止まる", raises(lambda: campaign("c_bad3", ZTE, "m6_on")))


def fake_runner(eff):
    m = types.SimpleNamespace()

    def prepare(prob_, rd, op=None, **kw):
        rd = Path(rd); rd.mkdir(parents=True)
        info = {"design": {"L_ramp": 10.0, "warnings": []}}
        (rd / "prepare_info.json").write_text(json.dumps(info)); return info

    def run_staged(rd, stages, **kw):
        (Path(rd) / "res_6000.h5").write_bytes(b""); return 0

    def collect(prob_, rd, out_dir=None, rc=None, require_residual_pass=False):
        st = {k: {"verdict": "STEADY"} for k in ("C_T_with_shear", "C_T", "C_L", "C_M")}
        return {"forge_rc": 0, "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0, "step": 6000,
                "scalar_gradient_effective": "lsq", "slau_wall_normal_chi_effective": 1, "flag_policy": D.FLAG_POLICY,
                RS.ZTE_EFFECTIVE: eff, "zero_thickness_edge_velocity": {"mesh_signature": "msig-x", "field_hash": "fh-y"},
                "gates": {"verdict": "PASS", "fail_class": None, "reasons": [], "objective": "C_T_with_shear",
                          "residual": {"verdict": "NOT CONVERGED"}, "steadiness": {"series": st}}}
    m.prepare, m.run_staged, m.collect = prepare, run_staged, collect
    return m


_R = D.R
try:
    c_on = campaign("e_on", ZTE); D.R = fake_runner(ON2)
    r_on = c_on.evaluate([0.5] * 5, "doe_000")
    check("(e) 評価行: 要求と作動点ごとの実効値・格子署名・w のハッシュが台帳に残る",
          r_on["zero_thickness_edge_velocity"] == ON2 and r_on["ops"]["m6_on"][RS.ZTE_EFFECTIVE] == ON2
          and r_on["ops"]["m6_on"]["zero_thickness_edge_velocity_field_hash"] == "fh-y")
    check("(e) 処置ありのキャンペーンで処置ありの評価は学習に入る", c_on.summary()["n_pass"] == 1)
    D.R = _R; c_off = campaign("e_off"); D.R = fake_runner(ON2)
    c_off.evaluate([0.5] * 5, "doe_000")
    s_off = c_off.summary()
    check("(e) 処置なしのキャンペーンでは処置ありの評価 (PASS) を学習に入れない", s_off["n_pass"] == 0 and s_off["n_pass_excluded_by_policy"] == 1)
    # 処置ありのキャンペーンを、処置なしの台帳の上で再開する (旧評価がある状態) → 持ち越さない
    D.R = _R
    cmp_dir = tmp / "e_resume"; cmp_dir.mkdir()
    (cmp_dir / "camp").mkdir()
    (cmp_dir / "camp" / "ledger.jsonl").write_text(json.dumps(dict(row(dict.fromkeys(OPS2, "off"), legacy=True), tag="old")) + "\n")
    raw = json.loads(json.dumps(BASE)); raw["evaluate"]["zero_thickness_edge_velocity"] = ZTE
    (cmp_dir / "problem.yaml").write_text(yaml.safe_dump(raw))
    check("(e) 処置なしの旧台帳の上で処置ありを再開: 旧評価は学習に入らない", D.SernCampaign(cmp_dir / "problem.yaml", cmp_dir / "camp").summary()["n_pass"] == 0)
    raw["opt"][D.ZTE_CARRY_OVER_KEY] = ["m6_on", "m10_on"]
    (cmp_dir / "problem.yaml").write_text(yaml.safe_dump(raw))
    check("(e) 同じ台帳で作動点を名指しした持ち越し指定があれば採用", D.SernCampaign(cmp_dir / "problem.yaml", cmp_dir / "camp").summary()["n_pass"] == 1)
finally:
    D.R = _R

# ======================================================================== (f) restart の道具は /AUX を写さない
W_SRC, W_DST = np.array([1.0, 0.0, 1.0, 0.0, 1.0]), np.array([0.0, 1.0, 1.0, 1.0, 0.0])
CPG_CFG = "physProp: {thermalMethod: 0}\nmesh: {meshFileName: \"dst.h5\"}\n"


def aux(path):
    with h5py.File(path, "r") as f:
        return np.asarray(f["AUX/w_recon_vel"][...]) if "AUX/w_recon_vel" in f else None


def run_tool(script, *args):
    return subprocess.run([sys.executable, str(TOOLS / script), *map(str, args)], capture_output=True, text=True)


for label, w_dst in (("宛先に w あり", W_DST), ("宛先に w なし", None)):
    d = tmp / f"rf_{w_dst is None}"; d.mkdir(); (d / "solverConfig.yaml").write_text(CPG_CFG)
    make_h5(d / "src.h5", w=W_SRC, val=3.0); make_h5(d / "dst.h5", w=w_dst, val=1.0)
    r = run_tool("restart_field.py", d / "src.h5", d / "dst.h5", "--force-species")
    got = aux(d / "dst.h5")
    moved = float(h5py.File(d / "dst.h5", "r")["VALUE/ro"][0]) == 3.0
    check(f"(f) restart_field (同一格子, {label}): 保存量は写り、宛先の /AUX は保持・元の /AUX は写らない",
          r.returncode == 0 and moved and ((got is None) if w_dst is None else np.array_equal(got, w_dst)), r.stdout[-300:] + r.stderr[-300:])
    d2 = tmp / f"rbi_{w_dst is None}"; d2.mkdir(); (d2 / "solverConfig.yaml").write_text(CPG_CFG); (d2 / "solverConfig_main.yaml").write_text(CPG_CFG)
    s2 = tmp / f"rbi_src_{w_dst is None}"; s2.mkdir(); (s2 / "solverConfig.yaml").write_text(CPG_CFG); (s2 / "solverConfig_main.yaml").write_text(CPG_CFG)
    make_h5(s2 / "res_100.h5", w=W_SRC, val=3.0); make_h5(d2 / "sern.h5", w=w_dst, val=1.0)
    try:
        RS.restart_by_index(s2 / "res_100.h5", d2 / "sern.h5")
        err = ""
    except Exception as e:  # noqa: BLE001
        err = f"{type(e).__name__}: {e}"
    got = aux(d2 / "sern.h5")
    check(f"(f) runner の restart_by_index (同一格子, {label}): 宛先の /AUX は保持・元の /AUX は写らない",
          not err and float(h5py.File(d2 / "sern.h5", "r")["VALUE/ro"][0]) == 3.0
          and ((got is None) if w_dst is None else np.array_equal(got, w_dst)), err)
    d3 = tmp / f"if_{w_dst is None}"; d3.mkdir(); (d3 / "solverConfig.yaml").write_text(CPG_CFG)
    coords_new = np.vstack([COORD * 1.01, [[0.2, 0.7, 0.0]]])
    make_h5(d3 / "src.h5", w=W_SRC, val=3.0)
    make_h5(d3 / "dst.h5", coords=coords_new, w=(np.append(W_DST, 1.0) if w_dst is not None else None), val=1.0)
    r = run_tool("interp_field.py", d3 / "src.h5", d3 / "dst.h5", "--force-species")
    got = aux(d3 / "dst.h5")
    check(f"(f) interp_field (別格子, {label}): 元の /AUX を移植しない (宛先のものは不変)",
          r.returncode == 0 and ((got is None) if w_dst is None else np.array_equal(got, np.append(W_DST, 1.0))), r.stdout[-300:] + r.stderr[-300:])

# ======================================================================== (g) 実際の道具で端から端まで (道具があるときだけ)
# 小さな node の 2D 格子 (6 節点、内部面 7、cowl_in = 線 [0,1]、cowl_out = 線 [1,2] が端点 1 を共有) に runner の
# mark_zte_field で道具の CLI を回し、再計算した識別量が道具の書いた属性 (= ソルバが照合する値) と一致することを確かめる
def make_node_h5(path):
    coord = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]], float); n = len(coord)
    edges = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
    with h5py.File(path, "w") as f:
        g = f.create_group("MESH")
        for k_, v_ in (("nNodes", n), ("nCells", n), ("nNormalPlanes", len(edges)), ("nPlanes", len(edges)), ("nBconds", 2), ("nBPlanes", 0)):
            g.attrs[k_] = v_
        g.create_dataset("COORD", data=coord.astype(np.float32).ravel())
        f.create_dataset("PLANES/STRUCT", data=np.array([x for a_, b_ in edges for x in (2, a_, b_, 2, a_, b_)], np.int32))
        f.create_dataset("CELLS/volume", data=np.ones(n, np.float32))
        v = f.create_group("VIZMESH"); v.attrs["nVizCells"] = 1; v.attrs["vizCONNE_dim"] = 5
        v.create_dataset("CONNE", data=np.array([5, 0, 1, 4, 3], np.int32))
        for pid, faces in ((5, [[0, 1]]), (6, [[1, 2]])):
            b_ = f.create_group(f"BCONDS/{pid}"); b_.attrs["bcondKind"] = "wall"
            b_.create_dataset("vizBfaceSizes", data=np.array([len(x) for x in faces], np.int32))
            b_.create_dataset("vizBfaceNodes", data=np.array([i for x in faces for i in x], np.int32))
        for k_ in ("ro", "roUx", "roUy", "roUz", "roe"):
            f.create_dataset(f"VALUE/{k_}", data=np.ones(n, np.float32))


if not (TOOLS / "mark_zero_thickness_edges.py").exists():
    skipped.append("(g) 道具が無い")
    print("SKIP (g): mark_zero_thickness_edges.py が無い")
else:
    _fake_mod = sys.modules.pop("mark_zero_thickness_edges")
    try:
        real = importlib.import_module("mark_zero_thickness_edges")
        d = tmp / "e2e"; d.mkdir()
        (d / "bcondConfig.yaml").write_text("cowl_in: {physID: 5, kind: wall}\ncowl_out: {physID: 6, kind: wall}\n")
        (d / "solverConfig.yaml").write_text(cfg_on); (d / "solverConfig_main.yaml").write_text(cfg_on)
        make_node_h5(d / "sern.h5")
        spec = {"tags": ["cowl_in", "cowl_out"], "rings": 1}
        try:
            a_ = RS.mark_zte_field(d, spec); err = ""
        except Exception as e:  # noqa: BLE001
            a_, err = None, f"{type(e).__name__}: {e}"
        check("(g) runner の mark_zte_field で道具の CLI が通る", a_ is not None, err + (d / "ZTE_MARK.txt").read_text()[-400:] if (d / "ZTE_MARK.txt").exists() else err)
        if a_ is not None:
            check("(g) 再計算した w のハッシュ・格子署名 = 道具が書いた属性 (ソルバが照合する値)",
                  a_["field_hash"] == a_["attrs"]["field_sha256"] and a_["mesh_signature"] == a_["attrs"]["mesh_signature"])
            check("(g) w は端点 1 から 1 リング (節点 0,1,2,4) が 0", np.array_equal(np.asarray(h5py.File(d / "sern.h5", "r")["AUX/w_recon_vel"][...]),
                                                                            np.array([0, 0, 0, 1, 0, 1], np.float32)))
            k1 = ZKEY(SM.stage_key(cfg_on, "", d))
            check("(g) stage_key に道具の値が入る", k1.get(SM.ZTE_PREFIX + ".field_hash") == a_["attrs"]["field_sha256"]
                  and k1.get(SM.ZTE_PREFIX + ".mesh_signature") == a_["attrs"]["mesh_signature"] and k1.get(SM.ZTE_PREFIX + ".rings") == "1")
            src = d / "res_100.h5"; make_node_h5(src)
            with h5py.File(src, "r+") as f:
                f["VALUE/ro"][...] = 3.0
            r = run_tool("restart_field.py", src, d / "sern.h5", "--force-species")
            check("(g) 同一格子の restart_field の後も stage_key の処置の部分は同じ (格子署名に /VALUE を含まない)",
                  r.returncode == 0 and ZKEY(SM.stage_key(cfg_on, "", d)) == k1, r.stdout[-300:] + r.stderr[-300:])
            RS._record_zte_launch(d); launch(d, cfg_on)
            (d / "forge_run.log").write_text(f"'zeroThicknessEdgeVelocity' effective: 1 (field /AUX/w_recon_vel, w<1 nodes 4, "
                                             f"field_sha256 {a_['field_hash'][:16]}, mesh_signature {a_['mesh_signature'][:16]})\n")
            check("(g) 来歴: 実際の道具の w で有効が確定", RS.zte_provenance(d)["effective"] == "tags=cowl_in,cowl_out;rings=1")
    finally:
        sys.modules["mark_zero_thickness_edges"] = _fake_mod

print(f"\n一時ディレクトリ: {tmp}")
if skipped:
    print("省略: " + "; ".join(skipped))
print("VERDICT:", ("PASS" if not skipped else "PASS (一部省略)") if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
