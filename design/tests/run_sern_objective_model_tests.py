#!/usr/bin/env python3
"""受理ゲート (4) の目的量・定常性の量を問題設定の乱流モデルから決める試験 (plan tooling-sern-te-wake-grid §5.1 R2f、
R2a の横展開)。forge は回さない。

旧実装の `sern_gates.objective_key` は「最後の履歴に C_T_with_shear の列があるか」で目的量と定常性の量の数を決めていたので、
判定区間の段を持たない SST キャンペーンでは列が欠けると C_T の 3 量で通った (R2a で judgment_gate について直した欠陥と同じ型)。

合格 (測る前に固定):
  (f1) 保存済みの記録で判定が変わらない: 3D の検証 run の写し 12 本 (notes/investigations/2026-10-11-{lsw10,g3-local,g4-local,c6})
       と case/46.sern_design にローカルで残っている run 10 本 (2D の Euler 3・SST 5、3D の SST 2。PASS と FAIL を含む) の
       force_history.csv・metrics.json で、変更前 (commit 15968c33 の sern_gates) と変更後の objective・定常性 (量の並び・各量の
       判定と詳細・未定常の量・理由)・evaluate_gates の verdict・fail_class が一致し、metrics.json の保存値 (objective・verdict・
       fail_class・各量の判定) とも一致する。他のゲート (場・残差・最終場の床・床事象・残差の桁) は metrics.json の保存値で差し替える
  (f2) SST で C_T_with_shear の列を消した履歴 (全標本から・最後の標本だけから) は判定不能 (FAIL/JUDGMENT_UNVERIFIABLE、延長の
       対象外。driver_sern は verdict が PASS でない作動点を台帳 FAIL にするので学習に入らない)。未定常が重なっても判定不能が先。
       変更前の実装では同じ履歴が C_T の 3 量で PASS していた (欠陥の再現)
  (f3) Euler 指定では 3 量 (目的量 C_T) で判定できる: 列を消した履歴も、列がある履歴も (C_T_with_shear を見ない)
  (f4) モデルの受け渡し: evaluate_gates・steadiness_gate は model を省略できない (TypeError)。None・未知のモデルは ValueError
       (列から推測する経路も、None を既定に読み替える経路も残さない)。旧引数 obj は廃止 (TypeError)
  (f5) collect (2D・3D) が問題設定のモデル (evaluate.model。無ければ runner の _solver_config と同じ既定 euler) を
       evaluate_gates に渡す

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_objective_model_tests.py
"""
import csv
import json
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
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

from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

fails = 0
n_checks = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_objective_model_"))
C46 = ROOT / "case" / "46.sern_design"
NOTES = ROOT / "notes" / "investigations"
PRE = "15968c33"          # R2f の変更前 (列の有無で目的量を決めていた sern_gates)


def check(name, ok, detail=""):
    global fails, n_checks
    n_checks += 1
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""))


def raises(fn, exc):
    try:
        fn()
    except exc:
        return True
    except Exception as e:  # noqa: BLE001
        print(f"    (別の例外 {type(e).__name__}: {e})")
        return False
    return False


def stop(msg):
    print("止める: " + msg)
    sys.exit(2)


# --- 変更前の sern_gates (commit PRE) を別モジュールとして読む ------------------------------------------------
_src = subprocess.run(["git", "-C", str(ROOT), "show", f"{PRE}:design/forge_design/metrics/sern_gates.py"], capture_output=True, text=True)
if _src.returncode != 0:
    stop(f"変更前の sern_gates ({PRE}) を git から読めない: {_src.stderr.strip()[:200]}")
OLD = types.ModuleType("forge_design.metrics._sern_gates_pre_r2f")
OLD.__package__ = "forge_design.metrics"
OLD.__file__ = G.__file__
exec(compile(_src.stdout, f"{PRE}:sern_gates.py", "exec"), OLD.__dict__)
if "model" in OLD.evaluate_gates.__code__.co_varnames or not hasattr(OLD, "objective_key"):
    stop(f"{PRE} の sern_gates が変更前の形 (列の有無で目的量を決める) でない")


# --- 保存済みの記録 -------------------------------------------------------------------------------------------
def load_hist(path, drop=(), drop_last=()):
    rows = list(csv.DictReader(open(path)))
    out = [{k: (int(v) if k == "step" else float(v)) for k, v in r.items() if v not in ("", None) and k not in drop} for r in rows]
    for k in drop_last:
        out[-1].pop(k, None)
    return out


def model_of_problem(name):
    return load_problem(C46 / name).evaluate.get("model", "euler")


# 3D の検証 run の写し: 問題 YAML は case/46.sern_design/diag/tewake_*.sh の記録 (全て SST の生産系列)
Y08_G3L, Y08_G3A = "problem_3d_prod_3op_wallres_lswx08_tewake_g3_L10.yaml", "problem_3d_prod_3op_wallres_lswx08_tewake_g3_A0.yaml"
Y08_G4L, Y10_G3L = "problem_3d_prod_3op_wallres_lswx08_tewake_g4_L10.yaml", "problem_3d_prod_3op_wallres_lswx10_tewake_g3_L10.yaml"
REC3D = [(NOTES / "2026-10-11-g3-local", Y08_G3L),
         (NOTES / "2026-10-11-g4-local" / "run_1085_tewake_g4_L10_cont20k", Y08_G4L),
         (NOTES / "2026-10-11-g4-local" / "run_1086_tewake_g4_L10_cont40k", Y08_G4L),
         (NOTES / "2026-10-11-c6" / "run_1087_tewake_m6_A0", Y08_G3A),
         (NOTES / "2026-10-11-c6" / "run_1088_tewake_m6_L10", Y08_G3L),
         (NOTES / "2026-10-11-c6" / "run_1089_tewake_m4_A0", Y08_G3A),
         (NOTES / "2026-10-11-c6" / "run_1090_tewake_m4_L10", Y08_G3L),
         (NOTES / "2026-10-11-c6" / "run_1091_tewake_m6_L10_cont40k", Y08_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1096_tewake_lsw10_m10_cont20k", Y10_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1097_tewake_lsw10_m6_cont20k", Y10_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1098_tewake_lsw10_m4_cont20k", Y10_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1099_tewake_lsw10_m6_cont40k", Y10_G3L)]
# case/46.sern_design にローカルで残っている run (run_* は追跡外なので、linked worktree では主チェックアウトの case/ を見る)
RUNS_LOCAL = ["run_0100_species_regress_euler_lumped", "run_0101_species_regress_euler_full", "run_0102_species_regress_euler_lumped_24k",
              "run_0123_cfl_A_base", "run_0128_cflR1_0p5", "run_0163_rankB_2000", "run_0198_rb_prop_A", "run_0308_wake_soft05",
              "run_0121_r5_vehicleside_lsw08", "run_0122_r5_outflow"]


def _local_case_dirs():
    yield C46
    r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--path-format=absolute", "--git-common-dir"], capture_output=True, text=True)
    if r.returncode == 0 and r.stdout.strip():
        yield Path(r.stdout.strip()).parent / "case" / "46.sern_design"


def find_local(name):
    for d in _local_case_dirs():
        if (d / name / "force_history.csv").exists() and (d / name / "metrics.json").exists():
            return d / name
    return None


records = []          # (ラベル, 履歴の CSV, metrics.json の dict, model)
for d, yml in REC3D:
    if not ((d / "force_history.csv").exists() and (d / "metrics.json").exists()):
        stop(f"保存済みの記録 {d} の force_history.csv / metrics.json が無い")
    records.append((d.name, d / "force_history.csv", json.loads((d / "metrics.json").read_text()), model_of_problem(yml)))
for name in RUNS_LOCAL:
    d = find_local(name)
    if d is None:
        stop(f"case/46.sern_design の {name} (force_history.csv・metrics.json) が見つからない ({[str(x) for x in _local_case_dirs()]})")
    info = json.loads((d / "prepare_info.json").read_text())
    if info.get("model") not in ("sst", "euler"):
        stop(f"{d}/prepare_info.json に model が無い ({info.get('model')!r})")
    records.append((f"{name} ({info.get('dim', 2)}D)", d / "force_history.csv", json.loads((d / "metrics.json").read_text()), info["model"]))
models = sorted({r[3] for r in records})
check(f"(f1) 記録 {len(records)} 本 (3D の写し {len(REC3D)}・case/46 {len(RUNS_LOCAL)})、モデル {models}、PASS と FAIL を含む",
      len(records) == 22 and models == ["euler", "sst"]
      and {r[2]["gates"]["verdict"] for r in records} == {"PASS", "FAIL"})

# 記録の時点でまだ無かったゲート (最終場の床・床事象・残差の桁) は合格の中立値で置く (当時の判定に無かったものを足さない)
NEUTRAL = {"floors": {"ok": True, "counts": {}, "reasons": []},
           "floor_events": {"ok": True, "required": False, "verdict": None, "fail_class": None, "reasons": [], "notes": []},
           "residual_scale": {"ok": True, "skipped": "保存値なし (このゲートの導入前の記録)"}}


def gates_from_saved(mod, saved, hist, **kw):
    """evaluate_gates を回す。保存場・残差・床・床事象・残差の桁は metrics.json の保存値 (run の HDF5 は写しに無い)。"""
    vals = {"field_health": saved["field"], "residual_health": saved["residual"], "floor_gate": saved.get("floors") or NEUTRAL["floors"],
            "floor_event_gate": saved.get("floor_events") or NEUTRAL["floor_events"],
            "residual_scale_gate": saved.get("residual_scale") or NEUTRAL["residual_scale"]}
    keep = {k: getattr(mod, k) for k in vals}
    for k, v in vals.items():
        setattr(mod, k, lambda *a, _v=v, **k_: _v)
    try:
        return mod.evaluate_gates(tmp / "saved", hist, saved.get("rc"), require_residual_pass=bool(saved.get("require_residual_pass")),
                                  require_residual_plateau=bool(saved.get("require_residual_plateau")), **kw)
    finally:
        for k, v in keep.items():
            setattr(mod, k, v)


def series_view(stead, with_detail=True):
    return [(k, v["verdict"]) + ((v["detail"],) if with_detail else ()) for k, v in stead["series"].items()]


# ======================================================================== (f1) 保存済みの記録で判定が変わらない
for label, path, met, model in records:
    hist = load_hist(path)
    saved = met["gates"]
    g_new = gates_from_saved(G, saved, hist, model=model)
    g_old = gates_from_saved(OLD, saved, hist)
    s_new, s_old = g_new["steadiness"], g_old["steadiness"]
    same_prepost = (g_new["verdict"] == g_old["verdict"] and g_new["fail_class"] == g_old["fail_class"]
                    and g_new["objective"] == g_old["objective"] and series_view(s_new) == series_view(s_old)
                    and s_new["unsteady"] == s_old["unsteady"] and s_new["ok"] == s_old["ok"] and g_new["reasons"] == g_old["reasons"]
                    and not s_new["missing"])
    same_saved = (g_new["verdict"] == saved["verdict"] and g_new["fail_class"] == saved["fail_class"]
                  and g_new["objective"] == saved["objective"] == met.get("objective")
                  and series_view(s_new, False) == series_view(saved["steadiness"], False))
    check(f"(f1) {label} [{model}]: 変更前後で objective {g_new['objective']}・定常性 {len(s_new['series'])} 量・"
          f"verdict {g_new['verdict']}/{g_new['fail_class']} が一致し、metrics.json の保存値とも一致",
          same_prepost and same_saved,
          {"new": (g_new["verdict"], g_new["fail_class"], g_new["objective"], series_view(s_new, False)),
           "old": (g_old["verdict"], g_old["fail_class"], g_old["objective"], series_view(s_old, False)),
           "saved": (saved["verdict"], saved["fail_class"], saved["objective"], series_view(saved["steadiness"], False))})

# ======================================================================== (f2) SST で列が欠けると判定不能
R1096 = NOTES / "2026-10-11-lsw10" / "run_1096_tewake_lsw10_m10_cont20k"
saved96 = json.loads((R1096 / "metrics.json").read_text())["gates"]
hA = load_hist(R1096 / "force_history.csv")
for tag, hX in (("全標本から", load_hist(R1096 / "force_history.csv", drop=("C_T_with_shear",))),
                ("最後の標本だけから", load_hist(R1096 / "force_history.csv", drop_last=("C_T_with_shear",)))):
    g_new = gates_from_saved(G, saved96, hX, model="sst")
    g_old = gates_from_saved(OLD, saved96, hX)
    ext, why = G.judgment_extendable(g_new)
    check(f"(f2) run_1096 の履歴の C_T_with_shear を{tag}消す [sst]: 判定不能 (FAIL/JUDGMENT_UNVERIFIABLE)・延長の対象外",
          g_new["verdict"] == "FAIL" and g_new["fail_class"] == "JUDGMENT_UNVERIFIABLE" and g_new["objective"] == "C_T_with_shear"
          and g_new["steadiness"]["missing"] == ["C_T_with_shear"] and not g_new["steadiness"]["ok"]
          and any("C_T_with_shear" in r and "列が力の履歴に無い" in r for r in g_new["reasons"]) and not ext and "必須量の列" in why,
          (g_new["verdict"], g_new["fail_class"], g_new["reasons"][:2], why))
    check(f"(f2) 同じ履歴 ({tag}) は変更前の実装では C_T の 3 量で PASS していた (欠陥の再現)",
          g_old["verdict"] == "PASS" and g_old["objective"] == "C_T" and list(g_old["steadiness"]["series"]) == ["C_T", "C_L", "C_M"],
          (g_old["verdict"], g_old["objective"], list(g_old["steadiness"]["series"])))
# 2D の SST でも同じ
d123 = find_local("run_0123_cfl_A_base")
saved123 = json.loads((d123 / "metrics.json").read_text())["gates"]
h123B = load_hist(d123 / "force_history.csv", drop=("C_T_with_shear",))
g123 = gates_from_saved(G, saved123, h123B, model="sst")
check("(f2) 2D の SST (run_0123) でも列を消すと判定不能 (変更前は C_T で PASS)",
      g123["fail_class"] == "JUDGMENT_UNVERIFIABLE" and gates_from_saved(OLD, saved123, h123B)["objective"] == "C_T"
      and gates_from_saved(OLD, saved123, h123B)["verdict"] == "PASS", (g123["fail_class"], g123["reasons"][:2]))
check("(f2) 判定不能の理由は列の欠落だけ (他の 3 量の定常性は判定して載せる)",
      sorted(gates_from_saved(G, saved96, load_hist(R1096 / "force_history.csv", drop=("C_T_with_shear",)), model="sst")
             ["steadiness"]["series"]) == ["C_L", "C_M", "C_T"])
# 欠落と未定常が同時なら判定不能が先 (延長しない)
hBd = load_hist(R1096 / "force_history.csv", drop=("C_T_with_shear",))
for i, h in enumerate(hBd):
    h["C_L"] = h["C_L"] + 1.0e-3 * i
gBd = gates_from_saved(G, saved96, hBd, model="sst")
check("(f2) 列の欠落 + C_L のドリフト → fail_class は判定不能 (UNSTEADY にしない = 延長の対象にしない)",
      gBd["fail_class"] == "JUDGMENT_UNVERIFIABLE" and "C_L" in gBd["steadiness"]["unsteady"] and not G.judgment_extendable(gBd)[0],
      (gBd["fail_class"], gBd["steadiness"]["unsteady"]))

# ======================================================================== (f3) Euler 指定は 3 量
hB = load_hist(R1096 / "force_history.csv", drop=("C_T_with_shear",))
for tag, hX in (("列を消した履歴", hB), ("列がある履歴", hA)):
    gE = gates_from_saved(G, saved96, hX, model="euler")
    check(f"(f3) Euler 指定・{tag}: 3 量 (C_T・C_L・C_M、目的量 C_T) で PASS (C_T_with_shear を見ない)",
          gE["verdict"] == "PASS" and gE["objective"] == "C_T" and list(gE["steadiness"]["series"]) == ["C_T", "C_L", "C_M"]
          and gE["steadiness"]["missing"] == [], (gE["verdict"], gE["objective"], list(gE["steadiness"]["series"])))
check("(f3) objective_key: sst → C_T_with_shear、euler → C_T (R2a の表 JUDGMENT_QUANTITIES と同じ大小文字の扱い)",
      G.objective_key("sst") == "C_T_with_shear" and G.objective_key("euler") == "C_T" and G.objective_key(" SST ") == "C_T_with_shear")
st_sst = G.steadiness_gate(hA, "sst")
check("(f3) steadiness_gate の量は必須量の表から: SST は目的量を先頭に 4 量 (変更前と同じ並び)",
      list(st_sst["series"]) == ["C_T_with_shear", "C_T", "C_L", "C_M"] and st_sst["model"] == "sst"
      and set(st_sst["series"]) == set(G.JUDGMENT_QUANTITIES["sst"]), list(st_sst["series"]))

# ======================================================================== (f4) モデルを省略・推測しない
run_any = tmp / "nothing"
check("(f4) evaluate_gates は model を省略できない (TypeError。列から推測しない)", raises(lambda: G.evaluate_gates(run_any, hA, 0), TypeError))
check("(f4) steadiness_gate も model を省略できない (TypeError)", raises(lambda: G.steadiness_gate(hA), TypeError))
check("(f4) model=None は拒否 (ValueError。既定 euler に読み替えない — collect が問題設定から解決して渡す)",
      raises(lambda: G.evaluate_gates(run_any, hA, 0, model=None), ValueError) and raises(lambda: G.steadiness_gate(hA, None), ValueError)
      and raises(lambda: G.objective_key(None), ValueError))
check("(f4) 未知のモデル (rans・空) は拒否 (ValueError)",
      raises(lambda: G.evaluate_gates(run_any, hA, 0, model="rans"), ValueError) and raises(lambda: G.objective_key(""), ValueError))
check("(f4) 履歴が無くてもモデルは先に検査する (未知のモデルを NO_FORCES に紛れさせない)",
      raises(lambda: G.evaluate_gates(run_any, [], 0, model="rans"), ValueError))
check("(f4) 旧引数 obj は廃止 (TypeError)", raises(lambda: G.evaluate_gates(run_any, hA, 0, obj="C_T", model="sst"), TypeError))
check("(f4) 旧呼び出し steadiness_gate(hist, 目的量) は未知のモデルとして拒否 (ValueError)",
      raises(lambda: G.steadiness_gate(hA, "C_T_with_shear"), ValueError))
g_empty = gates_from_saved(G, saved96, [], model="sst")
check("(f4) 履歴なしは従来どおり NO_FORCES (目的量はモデルから C_T_with_shear、量の欠落とは数えない)",
      g_empty["fail_class"] == "NO_FORCES" and g_empty["objective"] == "C_T_with_shear" and g_empty["steadiness"]["missing"] == []
      and g_empty["reasons"] == ["no force history (no wall output)"]
      and gates_from_saved(OLD, saved96, [])["fail_class"] == "NO_FORCES", (g_empty["fail_class"], g_empty["reasons"]))

# ======================================================================== (f5) collect がモデルを渡す
seen = []


def _capture(*a, **k):
    seen.append(k.get("model", "<なし>"))
    return {"steadiness": {"series": {}}, "objective": "C_T", "verdict": "FAIL", "fail_class": "NO_FORCES", "reasons": []}


INFO = {"states": {"exhaust": {"ro": 1.0, "u": 1.0, "P": 1.0}, "ext": {"P": 1.0}}, "H_m": 1.0, "F_ideal_N_per_m": 1.0,
        "design": {"L_ramp": 1.0}, "moc_forces": {}, "discretization": "node", "half_W_m": 1.0}
saved_fns = (R2.force_history, R2.evaluate_gates, R3.evaluate_gates)
try:
    R2.force_history = lambda *a, **k: []
    R2.evaluate_gates = _capture
    R3.evaluate_gates = _capture
    for mod, dim in ((R2, "2D"), (R3, "3D")):
        for ev, want in (({"model": "sst"}, "sst"), ({"model": "euler"}, "euler"), ({}, "euler")):
            d = tmp / f"collect_{dim}_{want}_{len(ev)}"; d.mkdir()
            (d / "prepare_info.json").write_text(json.dumps(INFO))
            py = d / "problem.yaml"
            py.write_text(yaml.safe_dump({"name": "t", "type": "sern_2d", "gas": {}, "spec": {}, "dv": {}, "geometry": {}, "mesh": {},
                                          "evaluate": ev}))
            n0 = len(seen)
            mod.collect(py, d, rc=0)
            check(f"(f5) collect ({dim}, evaluate {ev or '{}'}): evaluate_gates に model={want!r} を渡す",
                  len(seen) == n0 + 1 and seen[-1] == want, seen[n0:])
finally:
    R2.force_history, R2.evaluate_gates, R3.evaluate_gates = saved_fns

print(f"\n{n_checks - fails}/{n_checks} PASS" + ("" if not fails else f"  ({fails} FAIL)"))
print(f"(作業ディレクトリ {tmp})")
sys.exit(1 if fails else 0)
