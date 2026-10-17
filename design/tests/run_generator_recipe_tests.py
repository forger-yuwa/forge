#!/usr/bin/env python3
"""生成器の置き換えの試験 (plan tooling-generator-defaults-and-launch-gate §5.1 #11、§6.2 の「生成器の正例」)。

forge・変換器は使わない (段の config は `recipes/regen_profiles.py` の偽の forge で各生成器の全段の起動列を取る; 代表の問題は
`recipes/fixtures` の写し)。

(a) バイト同一: 全生成器 × 代表の問題の全段の起動列の config を、置き換えの前 (`fixtures/generator_launch_before_11.json`、
    ab9767fe で同じ台本が取ったもの) と比べる。例外は #16 の例外表 (238d2706 の `explicit`: 明示すべきキーと旧い実効値) の
    明示化だけ: (1) 変化は文字の挿入だけ (difflib の opcodes が equal・insert だけ)、(2) 足したキー = その段の例外表で、値も
    旧い実効値と同じ (実効値は変わらない)、(3) 前からあるキーの値は不変。
(b) 段の宣言: 各起動の直前に RECIPE.json があり、プロファイル = 機種、段 = その段、config の sha256 = 起動した config、逸脱なし、
    `explicit_added` = その段の例外表、段の写し RECIPE_<tag>.json = RECIPE.json。風洞 B10 の帰還 (euler_loop) は台本が起動しないので、
    パスの config (`euler_loop.pass_config`) を run に置いて宣言を確かめる。
(c) 置き換えた後の全段の config が自分のプロファイルに PASS (FAIL・WARN なし)。推奨設定のファイルは版を変えず、各段の `fixed` =
    置き換えの前の `fixed` ∪ `explicit` (値は不変)、`explicit` は空。
(d) 問題 YAML の上書きの逸脱への変換 (plan §4.3 の 4b): 代表の問題の写しの数値を変えて回すと、違う段・キーだけが理由つきの逸脱
    (出所 = 問題 YAML のパスとキーと値) になり、照合 (関門と同じ check_config) が通る。引数による上書き (axis-Mach の `--cfl`) は
    出所 = 引数。
(e) 関門: 本物の run_case.sh + 偽の forge (FORGE_BIN) で、生成器が書いた run (axis-Mach Euler の本段・ベル TP NS の本段) が起動し、
    LAUNCH_GATE.txt の設定・宣言に FAIL・WARN が無い (メッシュは合成の h5 と来歴の記録)。
(f) 負例: 出所の分からない違い (手で変えた CFL) は宣言を書かず RECIPE_<tag>.FAILED.json に残す (関門は宣言の無い run として扱う)、
    合うプロファイルの無い config (風洞 NS v1) は変えず、宣言もしない。前の段の宣言 (config が違う) は消す。
(h) 生成器の既定 (plan §5.1 #27): 問題 YAML に書いていないキーの生成器の既定 (SERN 2D の opt.soft_steps 3000・cfl_main 4.0、
    gas.transport 無しの viscMethod 1・thermCondMethod 1、driver_sern の opt.soft_steps 1500) が出所 generator_default つきの逸脱になり
    宣言が書ける (関門と同じ照合が通る)。既定は生成器ごと (driver_sern の 1500 は runner_sern.main の段としては FAILED)。

実行: design/.venv-opt/bin/python design/tests/run_generator_recipe_tests.py   (約 1 分)
"""
import contextlib
import difflib
import hashlib
import io
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TOOLS / "recipes"))
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(ROOT / "design"))
import recipe_settings as RS  # noqa: E402
import regen_profiles as RG  # noqa: E402
import yaml_strict as YS  # noqa: E402

from forge_design.evaluate import config_builder as CB  # noqa: E402

FAIL = 0
TMP = Path(tempfile.mkdtemp(prefix="gen_recipe_"))
BEFORE = json.loads((Path(__file__).resolve().parent / "fixtures" / "generator_launch_before_11.json").read_text(encoding="utf-8"))


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail and not cond else ""), flush=True)
    if not cond:
        FAIL += 1


def sha(text):
    return hashlib.sha256(text.encode() if isinstance(text, str) else text).hexdigest()


# ------------------------------------------------------------------------------------------------ 起動ごとの宣言を控える
LAUNCHES = []                                  # (run_dir の名前, config, RECIPE.json の中身 or None, RECIPE_<tag> の写し)
_orig_call = RG.Capture.__call__


def _capture(self, run_dir):
    rd = Path(run_dir)
    cfg = (rd / "solverConfig.yaml").read_text()
    decl = json.loads((rd / "RECIPE.json").read_text()) if (rd / "RECIPE.json").exists() else None
    tag = ((decl or {}).get("generator") or {}).get("tag")
    copy = (rd / f"RECIPE_{tag}.json").read_bytes() if tag and (rd / f"RECIPE_{tag}.json").exists() else None
    LAUNCHES.append((rd, cfg, decl, copy == (rd / "RECIPE.json").read_bytes() if decl else None))
    return _orig_call(self, run_dir)


RG.Capture.__call__ = _capture

text = RG.YAML_PATH.read_text(encoding="utf-8")
recipes = RS.load_recipes(text=text)
by_machine = {}
for sp in RG.SPECS:
    n0 = len(LAUNCHES)
    work = TMP / sp["name"]
    work.mkdir()
    problem = RG.materialize(sp, work)
    stages, bcond = sp["run"](problem, work)
    by_machine[sp["name"]] = {"stages": stages, "bcond": bcond, "launches": LAUNCHES[n0:], "work": work, "problem": problem}

# ------------------------------------------------------------------------------------------------ (a) バイト同一
print("== (a) 置き換えの前後のバイト同一 (例外 = 例外表のキーの挿入)")
n_stage = n_keys = 0
for name, m in by_machine.items():
    old = BEFORE["machines"][name]
    check(f"(a) {name}: 起動列の段の並びが前と同じ", [s for s, _ in m["stages"]] == [s for s, _ in old["stages"]],
          f"{[s for s, _ in m['stages']]} vs {[s for s, _ in old['stages']]}")
    for (st, new), (_, prev) in zip(m["stages"], old["stages"]):
        ops = difflib.SequenceMatcher(None, prev, new, autojunk=False).get_opcodes()
        only_insert = all(op in ("equal", "insert") for op, *_ in ops)
        fn, fo = RS.flatten(YS.load(new)), RS.flatten(YS.load(prev))
        added = {k: v for k, v in fn.items() if k not in fo}
        want = old["explicit"][st]
        same_vals = all(RS.same(fn.get(k), v) for k, v in fo.items()) and set(fo) <= set(fn)
        check(f"(a) {name}/{st}: 変化は文字の挿入だけ ({sum(1 for op, *_ in ops if op == 'insert')} か所)", only_insert,
              [op for op in ops if op[0] not in ("equal", "insert")][:3])
        check(f"(a) {name}/{st}: 足したキー = 例外表 ({len(want)} キー、値 = 旧い実効値)",
              set(added) == set(want) and all(RS.same(added[k], v) and type(added[k]) is type(v) for k, v in want.items()),
              f"足した {sorted(set(added) - set(want))[:4]}、足していない {sorted(set(want) - set(added))[:4]}、"
              f"値の違い {[(k, added.get(k), v) for k, v in want.items() if not RS.same(added.get(k), v)][:3]}")
        check(f"(a) {name}/{st}: 前からあるキーの値は不変", same_vals,
              [(k, fo[k], fn.get(k)) for k in fo if not RS.same(fn.get(k), fo[k])][:3])
        n_stage += 1
        n_keys += len(added)
print(f"     比べた段 {n_stage}、足したキーの延べ {n_keys} (例外表の延べ "
      f"{sum(len(e) for m in BEFORE['machines'].values() for e in m['explicit'].values())})")

# ------------------------------------------------------------------------------------------------ (b) 段の宣言
print("== (b) 段の宣言 (起動の直前の RECIPE.json)")
for name, m in by_machine.items():
    if name == "wt_feedback_euler_node":
        continue
    stages = [s for s, _ in m["stages"]]
    L = m["launches"]
    check(f"(b) {name}: 起動 {len(L)} 回 = 段 {len(stages)}", len(L) == len(stages))
    for (rd, cfg, decl, copy_ok), st in zip(L, stages):
        ok = (decl is not None and decl.get("profile") == name and decl.get("stage") == st and decl.get("kind") == "profile"
              and decl["config"]["sha256"] == sha(cfg) and decl.get("deviations") == []
              and decl.get("explicit_added") is not None and copy_ok)
        ex = BEFORE["machines"][name]["explicit"][st]
        check(f"(b) {name}/{st}: RECIPE.json (プロファイル・段・config の sha256 = 起動した config・逸脱なし・段の写し)", ok,
              json.dumps({k: (decl or {}).get(k) for k in ("profile", "stage", "deviations")}, ensure_ascii=False)[:200])
        check(f"(b) {name}/{st}: 宣言の explicit_added = 例外表", decl is not None and
              set((decl.get("explicit_added") or {})) == set(ex) and all(RS.same(decl["explicit_added"][k], v) for k, v in ex.items()))
# euler_loop: パスの config を run に置いて宣言
from forge_design.evaluate import runner_wt as W  # noqa: E402
from forge_design.feedback import euler_loop as EL  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402
m = by_machine["wt_feedback_euler_node"]
p = load_problem(m["problem"])
rd = TMP / "euler_loop_pass"
rd.mkdir()
cfg, fin = EL.pass_config(p, max(24000, 12000))
(rd / "solverConfig.yaml").write_text(cfg)
(rd / "bcondConfig.yaml").write_text(W._bcond(p, euler=True))
r = CB.declare_stage(rd, "main", "wt_feedback", finished=fin, problem=m["problem"], generator="euler_loop.run_loop")
decl = json.loads((rd / "RECIPE.json").read_text()) if (rd / "RECIPE.json").exists() else {}
check("(b) wt_feedback_euler_node/main: euler_loop のパスの config (pass_config) = 台本の本段、宣言 (逸脱なし)",
      cfg == m["stages"][0][1] and r["declared"] and decl.get("profile") == "wt_feedback_euler_node" and not decl.get("deviations")
      and decl["config"]["sha256"] == sha(cfg), str(r)[:200])

# ------------------------------------------------------------------------------------------------ (c) 自分のプロファイルに PASS
print("== (c) 置き換えた後の全段の config が自分のプロファイルに PASS")
for name, m in by_machine.items():
    for st, cfg in m["stages"]:
        r = RS.check_config(recipes, cfg, name, st, bcond=m["bcond"])
        check(f"(c) {name}/{st}: PASS (FAIL・WARN なし)", not r["fails"] and not r["warns"], f"{r['fails'][:3]} {r['warns'][:3]}")
old_rc = RS.load_recipes(text=subprocess.run(["git", "-C", str(ROOT), "show",
                                               f"{BEFORE['exception_table']['commit']}:{BEFORE['exception_table']['file']}"],
                                              capture_output=True, text=True).stdout) \
    if subprocess.run(["git", "-C", str(ROOT), "cat-file", "-e", BEFORE["exception_table"]["commit"]],
                      capture_output=True).returncode == 0 else None
if old_rc is None:
    print("     (238d2706 が無いので、版・fixed の比較は飛ばした)")
else:
    check("(c) 推奨設定のファイルの版は変えない", recipes.get("version") == old_rc.get("version"),
          f"{recipes.get('version')} vs {old_rc.get('version')}")
    for name in by_machine:
        a, b = RS.resolve(old_rc, name), RS.resolve(recipes, name)
        for st in a["stages"]:
            want = {**a["stages"][st]["fixed"], **a["stages"][st]["explicit"]}
            got = b["stages"][st]["fixed"]
            check(f"(c) {name}/{st}: fixed = 前の fixed ∪ explicit (値は不変)、explicit は空",
                  set(got) == set(want) and all(RS.same(got[k], v) for k, v in want.items()) and not b["stages"][st]["explicit"]
                  and a["stages"][st]["advisory"] == b["stages"][st]["advisory"],
                  f"{sorted(set(got) ^ set(want))[:4]}")


# ------------------------------------------------------------------------------------------------ (d) 問題 YAML の上書き → 逸脱
print("== (d) 問題 YAML の上書きの逸脱への変換 (plan §4.3 の 4b)")


def rerun(name, edit, tag, argv_extra=None):
    """代表の問題の写しを edit で変えて回す。戻り (run_dir, 変えた問題 YAML のパス)。"""
    sp = next(s for s in RG.SPECS if s["name"] == name)
    work = TMP / f"d_{name}_{tag}"
    work.mkdir()
    prob = RG.materialize(sp, work)
    raw = yaml.safe_load(prob.read_text(encoding="utf-8"))
    edit(raw)
    prob.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False), encoding="utf-8")
    if argv_extra is None:
        sp["run"](prob, work)
    else:
        from forge_design.evaluate import runner_axismach as A
        cap = RG.Capture()
        rd0 = work / "run"

        def fake_prepare(problem_path, run_dir, nsteps=None, ic_from=None, cfl_main=None, implicit_relax=None):
            p = load_problem(problem_path)
            run_dir = Path(run_dir)
            run_dir.mkdir(parents=True, exist_ok=False)
            n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
            out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
            (run_dir / "bcondConfig.yaml").write_text(A._bcond_with_species(W._bcond(p, euler=True), A._tp_species_Y(p)))
            (run_dir / "solverConfig.yaml").write_text(A._euler_config(p, run_dir, n, out_int, cfl_main, implicit_relax))
            return {}
        with RG.patched([(A, "prepare", fake_prepare), (A, "run_forge", cap), (A, "collect", lambda *a: {})] + RG._common_patches()):
            RG._quiet(A.main, [str(prob), str(rd0)] + argv_extra)
    rd = work / "run" if (work / "run").exists() else work / "campaign" / "regen"
    return rd, prob


def decls(rd):
    """run の段の写し RECIPE_<tag>.json → {tag: 宣言}。"""
    return {f.name[len("RECIPE_"):-len(".json")]: json.loads(f.read_text()) for f in sorted(rd.glob("RECIPE_*.json"))
            if not f.name.endswith(".FAILED.json") and f.name != "RECIPE_CONTEXT.json"}


def devs_of(d):
    return {x["key"]: x for x in d.get("deviations") or []}


CFL = ("time.deltaT.cfl", "time.deltaT.cfl_pseudo")


def expect(rd, prob, want):
    """want {tag: {キー: (値, {問題 YAML のキー: 値} or None, 引数の名前 or None)}}。段の宣言の逸脱が want とちょうど同じ。"""
    ds = decls(rd)
    rel = CB._rel(prob)
    for tag, d in ds.items():
        dv = devs_of(d)
        w = want.get(tag, {})
        good = set(dv) == set(w)
        for k, (val, keys, arg) in w.items():
            x = dv.get(k) or {}
            src = x.get("source") or {}
            good &= RS.same(x.get("value"), val) and bool(str(x.get("reason") or "").strip())
            if keys is not None:
                good &= src.get("problem") == rel and src.get("keys") == keys
            if arg is not None:
                good &= arg in json.dumps(src.get("argument"), ensure_ascii=False)
        check(f"(d) {rd.parent.name}/{tag}: 逸脱 = {sorted(w) or 'なし'} (出所つき)", good,
              json.dumps({k: {kk: x.get(kk) for kk in ("value", "source")} for k, x in dv.items()}, ensure_ascii=False)[:300])
    check(f"(d) {rd.parent.name}: 全段の宣言がある ({sorted(ds)})", set(want) <= set(ds) and not list(rd.glob("RECIPE_*.FAILED.json")),
          sorted(ds))
    return ds


def recheck(rd, name, ds, cfgs_by_tag):
    """段の宣言の逸脱で、関門と同じ check_config が段の config に通る。"""
    bc = YS.load((rd / "bcondConfig.yaml").read_text())
    for tag, d in ds.items():
        cfg = cfgs_by_tag.get(tag)
        if cfg is None:
            continue
        r = RS.check_config(recipes, cfg, name, d["stage"], {x["key"]: x["value"] for x in d["deviations"]}, bc)
        check(f"(d) {rd.parent.name}/{tag}: 宣言の逸脱で照合が通る (関門と同じ check_config)",
              not r["fails"] and sha(cfg) == d["config"]["sha256"], r["fails"][:3])


# SERN 2D: 本段の cfl_main と soft の CFL (mid も soft_cfl を使う)
LAUNCHES.clear()
rd, prob = rerun("sern2d_finite_base_frozen_tp_sst",
                 lambda r: (r["evaluate"].__setitem__("cfl_main", 0.75), r["opt"].__setitem__("soft_cfl", 0.2)), "cfl")
ds = expect(rd, prob, {"main": {k: (0.75, {"evaluate.cfl_main": 0.75}, None) for k in CFL},
                       "soft_cfl0.2": {k: (0.2, {"opt.soft_cfl": 0.2}, None) for k in CFL},
                       "mid_cfl0.2": {k: (0.2, {"opt.soft_cfl": 0.2}, None) for k in CFL},
                       "warm_cfl0.1": {}})
recheck(rd, "sern2d_finite_base_frozen_tp_sst", ds, {d["generator"]["tag"]: c for (_, c, d, _) in LAUNCHES if d})
# ベル TP NS: 本段の cfl_main と soft の段の長さ (evaluate.stage_steps)
LAUNCHES.clear()
rd, prob = rerun("bell_tp_ns_sst", lambda r: (r["evaluate"].__setitem__("cfl_main", 5.0),
                                              r["evaluate"].__setitem__("stage_steps", {"soft": 1500})), "cfl")
ds = expect(rd, prob, {"S3_main": {k: (5.0, {"evaluate.cfl_main": 5.0}, None) for k in CFL},
                       "S1_soft": {"time.last.nStepOuter": (1500, {"evaluate.stage_steps.soft": 1500}, None)},
                       "S2_mid": {}})
recheck(rd, "bell_tp_ns_sst", ds, {d["generator"]["tag"]: c for (_, c, d, _) in LAUNCHES if d})
# axis-Mach NS (main_line なし): implicit_relax は全段に効く
LAUNCHES.clear()
rd, prob = rerun("axismach_ns_sst_tp", lambda r: r["evaluate"].__setitem__("implicit_relax", 0.7), "relax")
ds = expect(rd, prob, {t: {"time.deltaT.implicitRelax": (0.7, {"evaluate.implicit_relax": 0.7}, None)}
                       for t in ("S1_soft", "S2_mid", "main")})
recheck(rd, "axismach_ns_sst_tp", ds, {d["generator"]["tag"]: c for (_, c, d, _) in LAUNCHES if d})
# axis-Mach Euler: 問題 YAML の cfl_main と、引数 --cfl (問題 YAML より優先)
LAUNCHES.clear()
rd, prob = rerun("axismach_euler_tp", lambda r: r["evaluate"].__setitem__("cfl_main", 3.0), "cfl")
ds = expect(rd, prob, {"main": {k: (3.0, {"evaluate.cfl_main": 3.0}, None) for k in CFL}, "soft": {}, "mid": {}})
recheck(rd, "axismach_euler_tp", ds, {d["generator"]["tag"]: c for (_, c, d, _) in LAUNCHES if d})
LAUNCHES.clear()
rd, prob = rerun("axismach_euler_tp", lambda r: None, "arg", argv_extra=["--cfl", "3.5"])
ds = expect(rd, prob, {"main": {k: (3.5, None, "--cfl") for k in CFL}, "soft": {}, "mid": {}})
recheck(rd, "axismach_euler_tp", ds, {d["generator"]["tag"]: c for (_, c, d, _) in LAUNCHES if d})

# ------------------------------------------------------------------------------------------------ (h) 生成器の既定
print("== (h) 生成器の既定を逸脱の出所に (plan §5.1 #27): 問題 YAML に書いていないキーの生成器の既定")


def gd_of(d, k):
    return ((devs_of(d).get(k) or {}).get("source") or {}).get("generator_default") or {}


def expect_gd(rd, want, label):
    """want {tag: {キー: (値, 関数, 問題 YAML のキー)}}。段の宣言の逸脱が want とちょうど同じで、出所が generator_default。"""
    ds = decls(rd)
    for tag, w in want.items():
        d = ds.get(tag) or {}
        dv = devs_of(d)
        good = bool(d) and set(dv) == set(w)
        for k, (val, fn, key) in w.items():
            g = gd_of(d, k)
            good &= RS.same((dv.get(k) or {}).get("value"), val) and g.get("function") == fn and g.get("key") == key
        check(f"(h) {label}/{tag}: 逸脱 = {sorted(w) or 'なし'} (出所 generator_default)", good,
              json.dumps({k: {kk: x.get(kk) for kk in ("value", "source")} for k, x in dv.items()}, ensure_ascii=False)[:300])
    check(f"(h) {label}: FAILED の宣言が無い", not list(rd.glob("RECIPE_*.FAILED.json")), [f.name for f in rd.glob("RECIPE_*.FAILED.json")])
    return ds


# runner_sern.main: opt.soft_steps を書かない → 既定 3000 (代表の問題は 4000)。evaluate.cfl_main を書かない → runner_sern.prepare の 4.0
LAUNCHES.clear()
rd, prob = rerun("sern2d_finite_base_frozen_tp_sst",
                 lambda r: (r["opt"].pop("soft_steps"), r["evaluate"].pop("cfl_main")), "gd_main")
ds = expect_gd(rd, {"soft_cfl0.1": {"time.last.nStepOuter": (3000, "runner_sern.main", "opt.soft_steps")},
                    "main": {k: (4.0, "runner_sern.prepare", "evaluate.cfl_main") for k in CFL},
                    "warm_cfl0.1": {}, "mid_cfl0.1": {}}, "SERN 2D runner_sern.main")
recheck(rd, "sern2d_finite_base_frozen_tp_sst", ds, {d["generator"]["tag"]: c for (_, c, d, _) in LAUNCHES if d})
# gas.transport を書かない → runner_sern._solver_config は viscMethod 1 (空気の Sutherland) + thermCondMethod 1。以前は viscMethod が
# 段の fixed (2) と違うだけでどのプロファイルにも合わず、宣言も FAILED の記録も書かなかった
LAUNCHES.clear()
rd, prob = rerun("sern2d_finite_base_frozen_tp_sst", lambda r: r["gas"].pop("transport"), "gd_transport")
_tr = {"physProp.viscMethod": (1, "runner_sern._solver_config", "gas.transport"),
       "physProp.thermCondMethod": (1, "runner_sern._solver_config", "gas.transport")}
ds = expect_gd(rd, {t: dict(_tr) for t in ("warm_cfl0.1", "soft_cfl0.1", "mid_cfl0.1", "main")}, "SERN 2D gas.transport なし")
recheck(rd, "sern2d_finite_base_frozen_tp_sst", ds, {d["generator"]["tag"]: c for (_, c, d, _) in LAUNCHES if d})


def run_driver_sern(edit, tag, generator=None):
    """opt/driver_sern.SernCampaign._eval_op の標準レシピ (prepare・forge・collect・設計の検査を差し替え)。generator を与えると、
    同じ段の値で runner_sern.run_staged を宣言の生成器の名前だけ変えて回す (既定が生成器ごとであることを見る)。"""
    from forge_design.evaluate import runner_sern as S
    from forge_design.opt import driver_sern as DS
    sp = next(x for x in RG.SPECS if x["name"] == "sern2d_finite_base_frozen_tp_sst")
    work = TMP / f"h_{tag}"
    work.mkdir()
    prob = RG.materialize(sp, work)
    raw = yaml.safe_load(prob.read_text(encoding="utf-8"))
    edit(raw)
    prob.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False), encoding="utf-8")
    cap = RG.Capture()
    pats = [(S, "prepare", RG._sern_fake_prepare(S, False)), (S, "run_forge", cap), (S, "restart_by_index", lambda res, mesh, **k: None),
            (S, "collect", lambda *a, **k: {"gates": {"verdict": "PASS", "fail_class": None, "reasons": []}}),
            (DS.SernCampaign, "_check_design", lambda self, info: None)]
    with RG.patched(pats):
        camp = DS.SernCampaign(prob, work / "campaign")
        op = camp.ops[0]["name"]
        if generator is None:
            rd_, rc_, out_, deg_ = RG._quiet(camp._eval_op, prob, "e1", op, {})
        else:                                       # _eval_op の標準レシピと同じ段の値 (driver_sern.py の ladder[0])
            oc = camp.optcfg
            rd_ = work / "campaign" / f"e1_{op}"
            S.prepare(prob, rd_, op=op)
            RG._quiet(S.run_staged, rd_, "full", soft_steps=int(oc.get("soft_steps", 1500)), soft_cfl=float(oc.get("soft_cfl", 0.5)),
                      warm_lam_steps=int(oc.get("warm_lam_steps", 0)), warm_lam_cfl=float(oc.get("warm_lam_cfl", 0.2)),
                      mid_steps=int(oc.get("mid_steps", 0)), warm_adapt_steps=int(oc.get("warm_adapt_steps", 500)),
                      recipe={"family": "sern2d", "problem": prob, "generator": generator})
    return Path(rd_)


rd = run_driver_sern(lambda r: r["opt"].pop("soft_steps"), "drv")
expect_gd(rd, {"soft_cfl0.1": {"time.last.nStepOuter": (1500, "opt.driver_sern._eval_op", "opt.soft_steps")},
               "warm_cfl0.1": {}, "mid_cfl0.1": {}, "main": {}}, "driver_sern の標準レシピ")
rd = run_driver_sern(lambda r: r["opt"].pop("soft_steps"), "drv_as_main", generator="runner_sern.main")
fj = rd / "RECIPE_soft_cfl0.1.FAILED.json"
check("(h) 既定は生成器ごと: driver_sern の 1500 を runner_sern.main (既定 3000) の段として宣言しようとすると FAILED",
      fj.exists() and {f["key"] for f in json.loads(fj.read_text())["fails"]} == {"time.last.nStepOuter"},
      fj.read_text()[:300] if fj.exists() else sorted(f.name for f in rd.glob("RECIPE_*")))
# 問題 YAML に書いてあれば既定は使わない (書いた値が出所): 代表の問題の値どおりなら逸脱なし
rd = run_driver_sern(lambda r: None, "drv_written")
expect_gd(rd, {"soft_cfl0.1": {}, "warm_cfl0.1": {}, "mid_cfl0.1": {}, "main": {}}, "driver_sern (opt.soft_steps を書いた問題)")

# ------------------------------------------------------------------------------------------------ (f) 負例
print("== (f) 負例: 出所の分からない違い・合うプロファイルの無い config")
src = by_machine["sern2d_finite_base_frozen_tp_sst"]
rdn = TMP / "f_hand_cfl"
shutil.copytree(src["work"] / "run", rdn)
cfg = (rdn / "solverConfig.yaml").read_text()
(rdn / "solverConfig.yaml").write_text(YS.replace_scalars(cfg, {("time", "deltaT", "cfl"): "0.9", ("time", "deltaT", "cfl_pseudo"): "0.9"}))
r = CB.declare_stage(rdn, "main", "sern2d", problem=src["problem"], generator="試験")
fj = rdn / "RECIPE_main.FAILED.json"
failed = json.loads(fj.read_text()) if fj.exists() else {}
check("(f) 手で変えた本段の CFL (問題 YAML に出所が無い) は宣言を書かず、FAILED に名指しで残す",
      not r["declared"] and not (rdn / "RECIPE.json").exists() and not (rdn / "RECIPE_main.json").exists()
      and {f["key"] for f in failed.get("fails", [])} == set(CFL), str(r)[:200])
# 関門は宣言の無い run として扱う (移行期は WARN)
sys.path.insert(0, str(TOOLS))
import launch_gate as LG  # noqa: E402
out = LG.Findings()
info = {}
raw = (rdn / "solverConfig.yaml").read_bytes()
doc, pf = RS.parse_config(raw.decode())
LG.check_settings(str(rdn), raw, doc, pf, YS.load((rdn / "bcondConfig.yaml").read_text()), {"enforce": False, "warn": False}, out, info)
check("(f) 関門: 宣言の無い run (移行期は WARN、強制なら FAIL)", info.get("declaration") == "なし"
      and any(a == "宣言" and lv == "WARN" for lv, a, _, _ in out.items))
# 合うプロファイルの無い config (風洞 NS v1 の本段) は変えず、宣言もしない。前の段の宣言 (別の config) は消す
P = type("P", (), {"cp": 1004.5, "gamma": 1.4})()
cfg_ns = W._config_sst_node(P, 12000, 4000, 4.0)
bc_slip = "inlet: {physID: 1, kind: inlet_Pressure}\nwall: {physID: 3, kind: wall}\n"
txt, fin = CB.finish_config(cfg_ns, bc_slip, "wt_modef")
check("(f) 風洞 NS v1 (node SST) は wt_modef のプロファイルに合わず、config は変えない", txt == cfg_ns and fin["profile"] is None)
rdv = TMP / "f_ns_v1"
rdv.mkdir()
(rdv / "solverConfig.yaml").write_text(cfg_ns)
(rdv / "bcondConfig.yaml").write_text(bc_slip)
(rdv / "RECIPE.json").write_text(json.dumps({"config": {"file": "solverConfig.yaml", "sha256": "0" * 64}}))
r = CB.declare_stage(rdv, "main", "wt_modef")
check("(f) 合うプロファイルの無い config は宣言せず、前の段の宣言 (config の sha256 が違う) を消す",
      not r["declared"] and not (rdv / "RECIPE.json").exists() and not list(rdv.glob("RECIPE_*")))

# ------------------------------------------------------------------------------------------------ (e) 関門 (run_case.sh + 偽の forge)
print("== (e) 関門: 本物の run_case.sh + 偽の forge で生成器の run が起動する")
import mesh_provenance as MP  # noqa: E402

FAKE = TMP / "fake_forge.sh"
FAKE.write_text("#!/bin/sh\n# 偽の forge (試験用): 起動されたら印を書くだけ\necho started > FORGE_STARTED\nexit 0\n")
FAKE.chmod(FAKE.stat().st_mode | stat.S_IXUSR)


def write_mesh(path, kinds, ni=21, nj=6, L=1.0, H=0.2):
    """node 変換と同じ置き方の合成の h5 (test_launch_gate.py と同じ形)。physID 1 = 左 (入口)、2 = 右 (出口)、3 = 上 (壁)、4 = 下 (軸)。"""
    xs, ys = np.linspace(0.0, L, ni), np.linspace(0.0, H, nj)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    coord = np.c_[X.ravel(), Y.ravel(), np.zeros(ni * nj)]
    nid = lambda i, j: i * nj + j                                                       # noqa: E731
    quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    edges = {1: [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)], 2: [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)],
             3: [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)], 4: [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)]}
    with h5py.File(path, "w") as f:
        g = f.create_group("MESH")
        g.attrs["nCells"] = np.int64(len(coord))
        f["MESH/COORD"] = coord.ravel()
        f["MESH/CONNE"] = np.array([x for n in range(len(coord)) for x in (1, n)], dtype=np.int64)
        f["VIZMESH/CONNE"] = np.array([x for q in quads for x in [5] + q], dtype=np.int64)
        f["VIZMESH"].attrs["nVizCells"] = np.int64(len(quads))
        f["CELLS/volume"] = np.ones(len(coord))
        for k in ("ro", "roUx", "roUy", "roUz", "roe"):
            f["VALUE/" + k] = np.ones(len(coord), np.float32)
        f["VALUE/wall_dist"] = (H - coord[:, 1])
        for pid, es in edges.items():
            b = f.create_group(f"BCONDS/{pid}")
            b["iCells"] = np.array(sorted({n for e in es for n in e}), dtype=np.int64)
            b["vizBfaceNodes"] = np.array([n for e in es for n in e], dtype=np.int64)
            b["vizBfaceSizes"] = np.full(len(es), 2, dtype=np.int64)
            b.attrs["bcondKind"] = kinds[pid]


def run_case(d):
    env = dict(os.environ, FORGE_BIN=str(FAKE))
    for k in ("FORGE_LAUNCH_GATE", "FORGE_STRICT_CONVERGENCE"):
        env.pop(k, None)
    p = subprocess.run(["bash", str(TOOLS / "run_case.sh"), str(d)], capture_output=True, text=True, env=env)
    gate = (d / "LAUNCH_GATE.txt").read_text(encoding="utf-8") if (d / "LAUNCH_GATE.txt").exists() else ""
    return p.returncode, p.stdout + p.stderr, (d / "FORGE_STARTED").exists(), gate


for name in ("axismach_euler_tp", "bell_tp_ns_sst"):
    m = by_machine[name]
    d = TMP / f"e_{name}"
    shutil.copytree(m["work"] / "run", d, ignore=shutil.ignore_patterns("res_*", "residual_*", "CONVERGENCE_VERDICT.txt"))
    cfg_t = (d / "solverConfig.yaml").read_text()
    bc_t = (d / "bcondConfig.yaml").read_text()
    bc = YS.load(bc_t)
    kinds = {int(v["physID"]): v["kind"] for v in bc.values()}
    mesh = d / YS.load(cfg_t)["mesh"]["meshFileName"]
    mesh.unlink(missing_ok=True)
    write_mesh(mesh, kinds)
    MP.write_record(str(mesh), MP.build_record(str(mesh), cfg_t, bc_t))
    rc, out_txt, started, gate = run_case(d)
    gl = [ln for ln in gate.splitlines() if ln.startswith(("FAIL", "WARN"))]
    bad = [ln for ln in gl if "[設定]" in ln or "[宣言]" in ln]
    check(f"(e) {name}: run_case.sh が関門を通って起動 (rc 0)", rc == 0 and started and "[run_case] forge exit=0" in out_txt,
          f"rc {rc} {gate[-800:]}")
    check(f"(e) {name}: LAUNCH_GATE.txt の宣言 = 生成器の宣言 (プロファイル・段)、設定・宣言の FAIL・WARN なし",
          f"profile     : {name}" in gate and "stage       : main" in gate and not bad, bad[:3])
    print(f"     {name}: LAUNCH_GATE の FAIL・WARN = {gl}")

# ------------------------------------------------------------------------------------------------ (g) yaml_strict.insert_scalars
print("== (g) yaml_strict.insert_scalars (明示すべきキーの挿入: 文字の挿入だけ)")
T = ('mesh: {discretization: "node", meshFileName: "nozzle.h5"}\n'
     'physProp: {thermalMethod: 2,\n           cp: 1.0, transport: {A: cea}}\n'
     'time:\n  deltaT: {cfl: 4.0,\n           detectNaN: 1}\n  nStepInner: 5\n'
     'space: {}\n'
     '# コメント\n')
out, done = YS.insert_scalars(T, {("mesh", "scalarGradient"): '"lsq"', ("physProp", "thermoFloat"): "1",
                                  ("time", "deltaT", "implicitRelax"): "1.0", ("time", "extra"): "3", ("space", "venkatK"): "0.05"})
ops = difflib.SequenceMatcher(None, T, out, autojunk=False).get_opcodes()
check("(g) flow (1 行・複数行・入れ子・空 {})・block に足す: 挿入だけ、値はトークンどおり、元のコメントは残る",
      all(op in ("equal", "insert") for op, *_ in ops) and YS.load(out)["mesh"]["scalarGradient"] == "lsq"
      and YS.load(out)["physProp"]["transport"] == {"A": "cea"} and YS.load(out)["time"]["deltaT"]["implicitRelax"] == 1.0
      and YS.load(out)["time"]["extra"] == 3 and YS.load(out)["space"] == {"venkatK": 0.05} and "# コメント" in out
      and len(done) == 5, out)
for what, bad in (("既にあるキー", {("mesh", "discretization"): '"cell"'}), ("親が無い", {("nope", "x"): "1"}),
                  ("スカラーでないトークン", {("mesh", "x"): "[1, 2]"}), ("複数行のトークン", {("mesh", "x"): "1\n2"})):
    try:
        YS.insert_scalars(T, bad)
        check(f"(g) {what} は ValueError", False)
    except ValueError:
        check(f"(g) {what} は ValueError", True)
try:
    YS.insert_scalars("a: &x {p: 1}\nb: *x\n", {("b", "q"): "1"})
    check("(g) alias の親 (共有実体) は ValueError", False)
except ValueError:
    check("(g) alias の親 (共有実体) は ValueError", True)

if FAIL == 0:
    shutil.rmtree(TMP, ignore_errors=True)
print(f"\n{'ALL PASS' if FAIL == 0 else f'FAIL {FAIL}  (一時ディレクトリ {TMP} を残した)'}")
sys.exit(1 if FAIL else 0)
