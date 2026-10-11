#!/usr/bin/env python3
"""必須の量・列を乱流モデル (設定) から決める試験 (plan tooling-sern-te-wake-grid §5.1 R2g、R2f の洗い出しの残り)。forge は回さない。

R2a・R2f で直した欠陥 (「ファイル・列にあれば検査する」作りで、量が欠けると判定対象が減って通る) の残り 5 箇所:
`sern_gates.field_health` (roK・roOmega)、`sern_gates.floor_gate` (roOmega の床)、`check_convergence.py` (rms_roK・rms_roOmega)、
乱流モデル名の大小文字 (runner は完全一致、ゲートは小文字化)、`driver_sern._op_summary` の C_T 欄。

合格 (測る前に固定):
  (g1) field_health: SST で最終 h5 に roK・roOmega のどちらかが無ければ判定不能 (ok False・missing・理由に「判定不能」)。
       evaluate_gates は FAIL/JUDGMENT_UNVERIFIABLE で延長の対象外。流体の量 (ro・roUx・roUy・roUz・roe・P・T) の欠落は
       どちらのモデルでも判定不能。Euler では roK・roOmega が無くても従来どおり通る。model の省略は TypeError、None・未知は ValueError
  (g2) floor_gate: SST で roOmega が無ければ不合格 (counts「roOmega 不在」= -1。P・T・ro が無いときの既存の扱いと同じ)。
       Euler で roOmega が無ければ従来どおり見ない。roOmega があれば両モデルとも床を数える
  (g3) check_convergence: run の solverConfig.yaml の乱流モデルが SST (sst・sst-ddes・sst-iddes、旧書式 LESorRANS 2 + RANSmodel 1)
       なら rms_roK・rms_roOmega を必須にし、欠けたら判定不能 (analyze・--from-floor・CLI)。乱流なし (none・旧書式 LESorRANS 0) は
       従来どおり 5 列で判定。solverConfig.yaml が無い・読めない・モデルが未知なら SST の 2 列も必須 (最も厳しい集合。揃って
       いれば判定し、欠けていれば判定不能) — 黙って 5 列に緩めない
  (g4) 乱流モデル名: 問題 YAML の evaluate.model は sst・euler (小文字) だけ。Euler・SST・前後の空白・null・未知は load_problem で拒否し、
       judgment_plan (driver_sern の学習の要件も) も同じ関数 (probdef.sern_turbulence_model) で拒否する。受理した値では runner の
       solver config (turbulence.model) と受理ゲートの量が同じモデルを指す
  (g5) driver_sern._op_summary: 目的量の値が out に無ければ C_T 欄は None (欠測。圧力の C_T を入れない)。目的量が無いときも None
  (g6) 保存済みの記録で判定が変わらない: case/46.sern_design の run (RUNS_LOCAL の 10 本。--all で metrics.json と res_*.h5 の
       ある全 run) で、変更前 (commit PRE の sern_gates) と変更後の field_health・floor_gate・evaluate_gates (verdict・fail_class・
       reasons) が一致する。3D の写し 12 本 (notes/investigations/2026-10-11-*) は h5 が無いので保存値で確かめる: 保存された
       floors の counts に roOmega の行があり (SST で roOmega があった)、field が ok
  (g7) check_convergence の CLI: SST 8 本・乱流なし 6 本 (cell・node、旧書式、--segment、--from-floor を含む、case/46 以外も) で
       変更前 (commit PRE) と変更後の VERDICT 行・OVERALL・終了コードが一致する (--segment は必要なファイルを作業ディレクトリに写して回す)

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_required_state_tests.py [--all] [--no-cli]
"""
import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import h5py
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

import check_convergence as CC  # noqa: E402
from forge_design import probdef as PD  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--all", action="store_true", help="(g6) case/46.sern_design の metrics.json と res_*.h5 がある全 run も比べる")
ap.add_argument("--no-cli", action="store_true", help="(g7) の CLI 比較を省く")
args = ap.parse_args()

fails = 0
n_checks = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_required_state_"))
C46 = ROOT / "case" / "46.sern_design"
NOTES = ROOT / "notes" / "investigations"
PRE = "1595e86b"          # R2g の変更前 (R2f の実装。roK・roOmega・rms_roK・rms_roOmega を「あれば検査する」)


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


def git_show(path):
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{PRE}:{path}"], capture_output=True, text=True)
    if r.returncode != 0:
        stop(f"変更前の {path} ({PRE}) を git から読めない: {r.stderr.strip()[:200]}")
    return r.stdout


# --- 変更前の sern_gates (commit PRE) を別モジュールとして読む ------------------------------------------------
OLD = types.ModuleType("forge_design.metrics._sern_gates_pre_r2g")
OLD.__package__ = "forge_design.metrics"
OLD.__file__ = G.__file__
exec(compile(git_show("design/forge_design/metrics/sern_gates.py"), f"{PRE}:sern_gates.py", "exec"), OLD.__dict__)
if "model" in OLD.field_health.__code__.co_varnames or "model" not in OLD.evaluate_gates.__code__.co_varnames:
    stop(f"{PRE} の sern_gates が変更前の形 (evaluate_gates は model 必須、field_health は model なし) でない")


def _local_case_dirs():
    """case/ の run_* は追跡外なので、linked worktree では主チェックアウトの case/ も見る。"""
    yield ROOT / "case"
    r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--path-format=absolute", "--git-common-dir"], capture_output=True, text=True)
    if r.returncode == 0 and r.stdout.strip():
        yield Path(r.stdout.strip()).parent / "case"


def find_run(rel, need=("residual_history.csv",)):
    for d in _local_case_dirs():
        if all((d / rel / n).exists() for n in need):
            return d / rel
    return None


def load_hist(path):
    rows = list(csv.DictReader(open(path)))
    return [{k: (int(v) if k == "step" else float(v)) for k, v in r.items() if v not in ("", None)} for r in rows]


def _last_res(d):
    return sorted(Path(d).glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]


# ======================================================================== 合成の run (SST 2D node の実 run を写して量を消す)
SRC = find_run("46.sern_design/run_0123_cfl_A_base", ("residual_history.csv", "force_history.csv", "solverConfig.yaml", "metrics.json"))
if SRC is None:
    stop("case/46.sern_design/run_0123_cfl_A_base (SST 2D node) が見つからない")


def synth_run(name, drop=(), cfg=None):
    """run_0123 の最終 h5・残差・設定を作業ディレクトリに写し、VALUE から drop の量を消す (run_0123 自体には書かない)。"""
    d = tmp / name
    d.mkdir()
    shutil.copy2(_last_res(SRC), d / "res_12000.h5")
    for n in ("residual_history.csv", "solverConfig.yaml"):
        shutil.copy2(SRC / n, d / n)
    if cfg is not None:
        (d / "solverConfig.yaml").write_text(cfg)
    if drop:
        with h5py.File(d / "res_12000.h5", "r+") as f:
            for k in drop:
                del f[f"VALUE/{k}"]
    return d


HIST = load_hist(SRC / "force_history.csv")
base = synth_run("full")
with h5py.File(base / "res_12000.h5", "r") as f:
    keys0 = set(f["VALUE"].keys())
check("(準備) 写した SST の最終 h5 に流体の量と roK・roOmega がそろう",
      set(G.FIELD_REQUIRED["sst"]) <= keys0 and set(G.FIELD_REQUIRED["euler"]) < set(G.FIELD_REQUIRED["sst"]), sorted(keys0))

# ======================================================================== (g1) field_health
fh_full = {m: G.field_health(base, model=m) for m in ("sst", "euler")}
check("(g1) 量がそろった h5 は SST・Euler とも従来どおり (変更前と ok・理由・非有限・非正・最小値が一致)",
      all(fh_full[m]["ok"] and fh_full[m]["missing"] == [] for m in fh_full)
      and all({k: fh_full[m][k] for k in ("ok", "reasons", "nonfinite", "nonpositive", "min", "species_checked", "file")}
              == {k: OLD.field_health(base)[k] for k in ("ok", "reasons", "nonfinite", "nonpositive", "min", "species_checked", "file")}
              for m in fh_full), {m: (v["ok"], v["reasons"]) for m, v in fh_full.items()})
for drop in (("roK",), ("roOmega",), ("roK", "roOmega")):
    d = synth_run("no_" + "_".join(drop), drop)
    fs, fe, fo = G.field_health(d, model="sst"), G.field_health(d, model="euler"), OLD.field_health(d)
    check(f"(g1) SST で {'・'.join(drop)} が無い → 判定不能 (変更前は素通り)",
          not fs["ok"] and fs["missing"] == list(drop) and any("判定不能" in r and drop[0] in r for r in fs["reasons"]) and fo["ok"],
          (fs["ok"], fs["missing"], fs["reasons"], "old", fo["ok"]))
    check(f"(g1) Euler で {'・'.join(drop)} が無い → 従来どおり通る", fe["ok"] and fe["missing"] == [], (fe["ok"], fe["reasons"]))
    g = G.evaluate_gates(d, HIST, 0, model="sst")
    ext, why = G.judgment_extendable(g)
    check(f"(g1) SST で {'・'.join(drop)} が無い → evaluate_gates は FAIL/JUDGMENT_UNVERIFIABLE・延長の対象外",
          g["verdict"] == "FAIL" and g["fail_class"] == "JUDGMENT_UNVERIFIABLE" and not ext and "保存場" in why
          and any("判定不能" in r for r in g["reasons"]), (g["verdict"], g["fail_class"], g["reasons"][:2], why))
for v in ("T", "roUz", "roe"):
    d = synth_run("no_" + v, (v,))
    check(f"(g1) 流体の量 {v} が無い → SST・Euler とも判定不能 (変更前は素通り)",
          all(not G.field_health(d, model=m)["ok"] and G.field_health(d, model=m)["missing"] == [v] for m in ("sst", "euler"))
          and OLD.field_health(d)["ok"])
d = synth_run("nonfinite_and_missing", ("roK",))
with h5py.File(d / "res_12000.h5", "r+") as f:
    a = f["VALUE/T"][:]; a[0] = float("nan"); f["VALUE/T"][:] = a
g = G.evaluate_gates(d, HIST, 0, model="sst")
check("(g1) 非有限の値と量の欠落が重なれば DIVERGED が先 (欠落だけのときに限り判定不能)",
      g["fail_class"] == "DIVERGED" and g["field"]["missing"] == ["roK"] and g["field"]["nonfinite"], (g["fail_class"], g["field"]["reasons"]))
check("(g1) field_health は model を省略できない (TypeError)", raises(lambda: G.field_health(base), TypeError))
check("(g1) model=None・未知のモデルは ValueError (既定に読み替えない)",
      raises(lambda: G.field_health(base, model=None), ValueError) and raises(lambda: G.field_health(base, model="rans"), ValueError))
check("(g1) res_*.h5 が無い run は従来どおり (field の不合格は DIVERGED のまま)",
      not G.field_health(tmp / "nothing", model="sst")["ok"] and G.field_health(tmp / "nothing", model="sst")["missing"] == [])

# ======================================================================== (g2) floor_gate
fg_full = {m: G.floor_gate(base, 20.0, model=m) for m in ("sst", "euler")}
fg_old = OLD.floor_gate(base, 20.0)
check("(g2) roOmega がある h5 は SST・Euler とも変更前と同じ (roOmega の床を数える)",
      all(fg_full[m] == fg_old for m in fg_full) and "roOmega<=1e-20" in fg_old["counts"], (fg_full["sst"], fg_old))
d = synth_run("floor_no_omega", ("roOmega",))
fs, fe, fo = G.floor_gate(d, 20.0, model="sst"), G.floor_gate(d, 20.0, model="euler"), OLD.floor_gate(d, 20.0)
check("(g2) SST で roOmega が無い → 不合格 (roOmega 不在 = -1。変更前は roOmega の床を見ずに通った)",
      not fs["ok"] and fs["counts"].get("roOmega 不在") == -1 and fo["ok"] and "roOmega 不在" not in fo["counts"], (fs, fo["ok"]))
check("(g2) Euler で roOmega が無い → 従来どおり見ない (変更前と同じ結果)", fe == fo and fe["ok"], (fe, fo))
check("(g2) floor_gate も model を省略できない・None・未知は ValueError",
      raises(lambda: G.floor_gate(base, 20.0), TypeError) and raises(lambda: G.floor_gate(base, 20.0, model=None), ValueError)
      and raises(lambda: G.floor_gate(base, 20.0, model="SST-X"), ValueError))
check("(g2) res_*.h5 が無い run は従来どおり判定不能 (モデルの検査はその前)",
      not G.floor_gate(tmp / "nothing", model="sst")["ok"] and raises(lambda: G.floor_gate(tmp / "nothing", model=None), ValueError))

# ======================================================================== (g3) check_convergence の必須列
COLS5 = ["step", "phase"] + ["rms_" + k for k in ("ro", "roUx", "roUy", "roUz", "roe")]
SER = [10.0 ** (-i / 4.0) for i in range(40)]


def conv_dir(name, cfg, cols):
    d = tmp / ("cc_" + name)
    d.mkdir()
    if cfg is not None:
        (d / "solverConfig.yaml").write_text(cfg)
    with open(d / "residual_history.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols); w.writerows([[i, "outer_end"] + [v] * (len(cols) - 2) for i, v in enumerate(SER)])
    return d


CFGS = {"sst flow": 'turbulence: {model: "sst", wallTreatmentSST: 0}\n',
        "sst block": "turbulence:\n  model: sst\n  scalarDiffusion: 1\n",
        "sst-ddes": "turbulence: {model: sst-ddes}\n",
        "旧書式 SST (LESorRANS 2 + RANSmodel 1)": "turbulence: {LESorRANS: 2, LESmodel: 1, RANSmodel: 1}\n"}
NONE_CFGS = {"none": 'turbulence: {model: "none"}\n', "wale": "turbulence: {model: wale}\n",
             "旧書式 乱流なし (LESorRANS 0)": "turbulence: {LESorRANS: 0, LESmodel: 1}\n"}
for i, (lab, cfg) in enumerate(CFGS.items()):
    for miss in (("rms_roK", "rms_roOmega"), ("rms_roK",), ("rms_roOmega",)):
        cols = COLS5 + [c for c in ("rms_roK", "rms_roOmega") if c not in miss]
        d = conv_dir(f"sst{i}_{len(miss)}_{miss[0]}", cfg, cols)
        r = CC.analyze(str(d / "residual_history.csv"), 3.0, 0.2)
        check(f"(g3) {lab} で {'・'.join(miss)} が無い → 判定不能",
              r is not None and r[2] is False and "(入力)" in r[1] and all(c in r[1]["(入力)"][0] for c in miss)
              and "判定不能" in r[1]["(入力)"][0], None if r is None else (r[2], r[1].get("(入力)")))
    d = conv_dir(f"sst{i}_full", cfg, COLS5 + ["rms_roK", "rms_roOmega"])
    r = CC.analyze(str(d / "residual_history.csv"), 3.0, 0.2)
    check(f"(g3) {lab} で 7 列がそろう → 従来どおり合格", r is not None and r[2] is True and "(入力)" not in r[1], None if r is None else r[1])
for i, (lab, cfg) in enumerate(NONE_CFGS.items()):
    d = conv_dir(f"none{i}", cfg, COLS5)
    r = CC.analyze(str(d / "residual_history.csv"), 3.0, 0.2)
    check(f"(g3) {lab} の 5 列 → 従来どおり合格 (SST の列を要求しない)", r is not None and r[2] is True, None if r is None else r[1])
for lab, cfg in (("solverConfig.yaml が無い", None), ("YAML として読めない", "turbulence: {model: [\n"),
                 ("turbulence.model が未知の綴り (SST)", 'turbulence: {model: "SST"}\n'), ("turbulence が無い", "gpu: 1\n")):
    d5 = conv_dir("unk5_" + str(abs(hash(lab))), cfg, COLS5)
    r5 = CC.analyze(str(d5 / "residual_history.csv"), 3.0, 0.2)
    d7 = conv_dir("unk7_" + str(abs(hash(lab))), cfg, COLS5 + ["rms_roK", "rms_roOmega"])
    r7 = CC.analyze(str(d7 / "residual_history.csv"), 3.0, 0.2)
    check(f"(g3) {lab}: 5 列は判定不能 (SST の 2 列も必須、理由を書く)・7 列は判定する (合格、注記あり)",
          r5 is not None and r5[2] is False and "rms_roK" in r5[1]["(入力)"][0] and "乱流モデルを決められない" in r5[1]["(入力)"][0]
          and r7 is not None and r7[2] is True and "(設定)" in r7[1] and r7[1]["(設定)"][1] is True,
          (None if r5 is None else r5[1].get("(入力)"), None if r7 is None else (r7[2], r7[1].get("(設定)"))))
# --from-floor: 対象が SST で 2 列欠落 (参照は 7 列、床は 5 列だけ渡す) → 判定不能
d = conv_dir("ff_sst5", CFGS["sst flow"], COLS5)
fl = {c: 1.0 for c in COLS5[2:]}
r = CC.analyze_from_floor(str(d / "residual_history.csv"), fl, set(), 1e30, 0.2)
check("(g3) --from-floor: SST で rms_roK・rms_roOmega が無い対象 → 判定不能",
      r is not None and r[2] is False and "rms_roK" in r[1].get("(入力)", ("",))[0], None if r is None else r[1].get("(入力)"))
dn = conv_dir("ff_none5", NONE_CFGS["none"], COLS5)
r = CC.analyze_from_floor(str(dn / "residual_history.csv"), fl, set(), 1e30, 0.2)
check("(g3) --from-floor: 乱流なしの 5 列 → 従来どおり合格", r is not None and r[2] is True, None if r is None else r[1])
# CLI (system の python3 でも動くこと: 他のセッションはこれで呼ぶ)
for py in (sys.executable, "/usr/bin/python3"):
    if not Path(py).exists():
        continue
    rr = subprocess.run([py, str(TOOLS / "check_convergence.py"), str(conv_dir(f"cli_{Path(py).name}_{abs(hash(py))}", CFGS["sst block"], COLS5))],
                        capture_output=True, text=True)
    out = rr.stdout + rr.stderr
    check(f"(g3) CLI ({py}): SST で 2 列欠落 → 非ゼロ終了・判定不能・NOT CONVERGED", rr.returncode == 1 and "判定不能" in out
          and "Traceback" not in out and "NOT CONVERGED" in out, out[-300:])
check("(g3) turbulence_kind: sst 系・旧書式は sst、none 系は other、無い・読めない・未知は None",
      CC.turbulence_kind(str(conv_dir("k_sst", CFGS["sst-ddes"], COLS5) / "residual_history.csv"))[0] == "sst"
      and CC.turbulence_kind(str(conv_dir("k_old", CFGS["旧書式 SST (LESorRANS 2 + RANSmodel 1)"], COLS5) / "residual_history.csv"))[0] == "sst"
      and CC.turbulence_kind(str(conv_dir("k_none", NONE_CFGS["旧書式 乱流なし (LESorRANS 0)"], COLS5) / "residual_history.csv"))[0] == "other"
      and CC.turbulence_kind(str(conv_dir("k_nocfg", None, COLS5) / "residual_history.csv"))[0] is None)

# ======================================================================== (g4) 乱流モデル名
check("(g4) sern_turbulence_model: 無ければ euler、sst・euler はそのまま",
      PD.sern_turbulence_model({}) == "euler" and PD.sern_turbulence_model(None) == "euler"
      and PD.sern_turbulence_model({"model": "sst"}) == "sst" and PD.sern_turbulence_model({"model": "euler"}) == "euler")
for bad in ("Euler", "SST", " sst", "sst ", "EULER", None, "rans", "none", 1):
    check(f"(g4) evaluate.model {bad!r} は拒否 (ValueError。読み替えない)", raises(lambda b=bad: PD.sern_turbulence_model({"model": b}), ValueError))
SMOKE = yaml.safe_load(open(C46 / "problem_smoke_euler.yaml"))
for val, ok in (("euler", True), ("sst", True), ("Euler", False), ("SST", False), (" sst", False), (None, False)):
    doc = json.loads(json.dumps(SMOKE)); doc["evaluate"]["model"] = val
    py = tmp / f"prob_{abs(hash(repr(val)))}.yaml"; py.write_text(yaml.safe_dump(doc, allow_unicode=True))
    if ok:
        p = PD.load_problem(py)
        cfg = R2._solver_config(p, 100, 50, 0.5, 0.0)
        tm = yaml.safe_load(cfg)["turbulence"]["model"]
        req = G.judgment_required_quantities(PD.sern_turbulence_model(p.evaluate))
        check(f"(g4) load_problem (evaluate.model {val!r}): solver config の turbulence.model {tm!r} とゲートの量 {len(req)} が同じモデル",
              (tm == "sst") == (val == "sst") and (len(req) == 4) == (val == "sst") and (tm == "none") == (val == "euler"), (tm, req))
    else:
        check(f"(g4) load_problem (evaluate.model {val!r}) は拒否 (ソルバは SST・ゲートは Euler の組を作らない)",
              raises(lambda q=py: PD.load_problem(q), ValueError))
doc = json.loads(json.dumps(SMOKE)); doc["evaluate"].pop("model", None)
py = tmp / "prob_nomodel.yaml"; py.write_text(yaml.safe_dump(doc, allow_unicode=True))
check("(g4) evaluate.model が無い問題 YAML は従来どおり euler (受理)",
      yaml.safe_load(R2._solver_config(PD.load_problem(py), 100, 50, 0.5, 0.0))["turbulence"]["model"] == "none")
PROD_EV = yaml.safe_load(open(C46 / "problem_3d_prod_3op_wallres_lswx10.yaml"))["evaluate"]
check("(g4) 生産 YAML (lswx10、sst) の judgment_plan は 4 量 (変更前と同じ)",
      R2.judgment_plan(PROD_EV)["required"] == G.JUDGMENT_QUANTITIES["sst"] and PROD_EV.get("judgment_steps"))
for bad in ("SST", "Euler", " sst"):
    ev = dict(PROD_EV, model=bad)
    check(f"(g4) judgment_plan・driver の学習の要件 (evaluate.model {bad!r}) は拒否",
          raises(lambda e=ev: R2.judgment_plan(e), ValueError) and raises(lambda e=ev: D.judgment_requirement_of(e), ValueError))

# ======================================================================== (g5) _op_summary の C_T 欄
S_ = D.SernCampaign._op_summary
c1 = S_({"gates": {"objective": "C_T_with_shear", "verdict": "FAIL"}, "C_T": 0.91}, "rd")
c2 = S_({"gates": {"objective": "C_T_with_shear", "verdict": "PASS"}, "C_T": 0.91, "C_T_with_shear": 0.88}, "rd")
c3 = S_({"gates": {"objective": "C_T", "verdict": "PASS"}, "C_T": 0.91}, "rd")
c4 = S_({"gates": {"verdict": "FAIL", "fail_class": "DIVERGED"}, "C_T": 0.91}, "rd")
c5 = S_({"gates": {"verdict": "FAIL", "fail_class": "DIVERGED"}}, "rd")
check("(g5) 目的量 (C_T_with_shear) の値が無い → C_T 欄は None (圧力の C_T を入れない)、C_T_p は圧力の C_T",
      c1["C_T"] is None and c1["C_T_p"] == 0.91 and c1["objective"] == "C_T_with_shear", c1["C_T"])
check("(g5) 目的量の値があればその値 (SST は C_T_with_shear、Euler は C_T)", c2["C_T"] == 0.88 and c3["C_T"] == 0.91, (c2["C_T"], c3["C_T"]))
check("(g5) 目的量が無い (起動段で落ちた等) → C_T 欄は None", c4["C_T"] is None and c5["C_T"] is None, (c4["C_T"], c5["C_T"]))

# ======================================================================== (g6) 保存済みの記録で判定が変わらない
RUNS_LOCAL = ["run_0100_species_regress_euler_lumped", "run_0101_species_regress_euler_full", "run_0102_species_regress_euler_lumped_24k",
              "run_0123_cfl_A_base", "run_0128_cflR1_0p5", "run_0163_rankB_2000", "run_0198_rb_prop_A", "run_0308_wake_soft05",
              "run_0121_r5_vehicleside_lsw08", "run_0122_r5_outflow"]
FIELD_KEYS = ("ok", "reasons", "nonfinite", "nonpositive", "min", "species_checked", "file", "nan_files")


def compare_run(d):
    info = json.loads((d / "prepare_info.json").read_text())
    model = info.get("model")
    if model not in ("sst", "euler"):
        return None, f"prepare_info.json の model {model!r}"
    met = json.loads((d / "metrics.json").read_text())
    saved = met.get("gates") or {}
    fn, fo = G.field_health(d, model=model), OLD.field_health(d)
    # collect が渡した p_min (問題 YAML の evaluate.p_min) は保存された床の counts の見出し "P<=X" から (無ければ config の pMin)
    pk = [k for k in ((saved.get("floors") or {}).get("counts") or {}) if k.startswith("P<=")]
    pmin = float(pk[0][3:]) if pk else None
    gn, go = G.floor_gate(d, pmin, model=model), OLD.floor_gate(d, pmin)
    hist = load_hist(d / "force_history.csv") if (d / "force_history.csv").exists() else []
    kw = dict(require_residual_pass=bool(saved.get("require_residual_pass")), require_residual_plateau=bool(saved.get("require_residual_plateau")),
              p_min=pmin)
    rc = met.get("forge_rc", saved.get("rc"))
    en, eo = G.evaluate_gates(d, hist, rc, model=model, **kw), OLD.evaluate_gates(d, hist, rc, model=model, **kw)
    same = ({k: fn.get(k) for k in FIELD_KEYS} == {k: fo.get(k) for k in FIELD_KEYS} and fn["missing"] == [] and gn == go
            and (en["verdict"], en["fail_class"], en["reasons"]) == (eo["verdict"], eo["fail_class"], eo["reasons"]))
    sv = saved.get("verdict")
    return (same, model, info.get("dim", 2), fn["ok"], gn["ok"], en["verdict"], en["fail_class"], sv, saved.get("fail_class")), None


names = list(RUNS_LOCAL)
if args.all:
    seen = set(names)
    for base_dir in _local_case_dirs():
        for m in sorted((base_dir / "46.sern_design").glob("run_*/metrics.json")):
            if m.parent.name not in seen and list(m.parent.glob("res_[0-9]*.h5")) and (m.parent / "prepare_info.json").exists():
                seen.add(m.parent.name); names.append(m.parent.name)
summary = {}
n_local = 0
# この比較は R2g の変更前後 (PRE) を見る。後から入れた残差の時系列の完全性 (§5.1 R2j) は外す: 発散で切り詰められた保存 run
# (run_0121 など) は R2j で理由が 1 つ増えるので、R2g の比較に混ぜると R2g の不変性を確かめられない。R2j の保存 run への影響は
# run_sern_residual_series_tests.py の (s) で別に数える
_series_check = G.residual_series_check
G.residual_series_check = lambda run_dir: ([], {})
for name in names:
    d = find_run("46.sern_design/" + name, ("metrics.json", "prepare_info.json"))
    if d is None or not list(d.glob("res_[0-9]*.h5")):
        if name in RUNS_LOCAL:
            stop(f"case/46.sern_design/{name} (metrics.json・res_*.h5) が見つからない")
        continue
    res, why = compare_run(d)
    if res is None:
        print(f"    (飛ばす {name}: {why})")
        continue
    same, model, dim, fok, gok, v, fc, sv, sfc = res
    n_local += 1
    summary[(model, dim, v)] = summary.get((model, dim, v), 0) + 1
    if name in RUNS_LOCAL or not same:
        check(f"(g6) case/46 {name} [{model} {dim}D]: 変更前後で field_health・floor_gate・evaluate_gates が一致 "
              f"(field {fok}・floors {gok}・verdict {v}/{fc}; 保存値 {sv}/{sfc})", same, res)
G.residual_series_check = _series_check
if args.all:
    check(f"(g6) --all: case/46 の {n_local} 本で一致 (内訳 {sorted(summary.items())})", n_local >= len(RUNS_LOCAL))
# 3D の写し (h5 は無い): 保存された床の counts に roOmega の行がある (SST で roOmega があった) と field が ok
REC3D = [NOTES / "2026-10-11-g3-local", NOTES / "2026-10-11-g4-local" / "run_1085_tewake_g4_L10_cont20k",
         NOTES / "2026-10-11-g4-local" / "run_1086_tewake_g4_L10_cont40k", NOTES / "2026-10-11-c6" / "run_1087_tewake_m6_A0",
         NOTES / "2026-10-11-c6" / "run_1088_tewake_m6_L10", NOTES / "2026-10-11-c6" / "run_1089_tewake_m4_A0",
         NOTES / "2026-10-11-c6" / "run_1090_tewake_m4_L10", NOTES / "2026-10-11-c6" / "run_1091_tewake_m6_L10_cont40k",
         NOTES / "2026-10-11-lsw10" / "run_1096_tewake_lsw10_m10_cont20k", NOTES / "2026-10-11-lsw10" / "run_1097_tewake_lsw10_m6_cont20k",
         NOTES / "2026-10-11-lsw10" / "run_1098_tewake_lsw10_m4_cont20k", NOTES / "2026-10-11-lsw10" / "run_1099_tewake_lsw10_m6_cont40k"]
ok3 = []
for d in REC3D:
    if not (d / "metrics.json").exists():
        stop(f"保存済みの記録 {d}/metrics.json が無い")
    g = json.loads((d / "metrics.json").read_text())["gates"]
    ok3.append((d.name, "roOmega<=1e-20" in (g.get("floors") or {}).get("counts", {}), g["field"]["ok"], g["verdict"]))
check(f"(g6) 3D の写し {len(REC3D)} 本: 保存された floors に roOmega の行があり (SST で roOmega があった)・field が ok "
      "(h5 は写しに無いので field_health の再計算はできない。verdict は R2f の試験と同じく保存値で一致)",
      all(x[1] and x[2] for x in ok3), [x for x in ok3 if not (x[1] and x[2])])

# ======================================================================== (g7) check_convergence の CLI: 変更前後の VERDICT
if not args.no_cli:
    old_tools = tmp / "old_tools"
    old_tools.mkdir()
    for f in ("check_convergence.py", "stage_manifest.py"):
        (old_tools / f).write_text(git_show(f"solver_density_cuda/tools/{f}"))
    CLI = [  # (ラベル, run, 追加引数)。--from-floor の参照も run (case/ からの相対)
        ("S1 sst node (現行の判定で PASS)", "40.nozzle_design_tool/run_0081_rao80_rd08_stgA", []),
        ("S2 sst cell (現行の判定で PASS)", "40.nozzle_design_tool/run_0079_rao80_cell_stgA", []),
        ("S3 sst (旧書式) cell", "36.passive_pseudoshock_control/run_sym_G_sst_bp170", []),
        ("S4 sst node (case/46)", "46.sern_design/run_0123_cfl_A_base", []),
        ("S5 sst node --segment", "45.isobutane_m6_d155/run_0119_rerun_ctrl", ["--segment"]),
        ("S6 sst + 遷移 node --segment", "57.transition_flat_plate/run_0004_t3a_lm", ["--segment"]),
        ("S7 sst node --from-floor (参照は PASS、自己参照)", "40.nozzle_design_tool/run_0081_rao80_rd08_stgA", ["--from-floor", "40.nozzle_design_tool/run_0081_rao80_rd08_stgA"]),
        ("S8 sst --from-floor (参照が PASS しない → REFUSED)", "48.flat_plate_cooled_m4/run_0019_cht_insolver_cont", ["--from-floor", "48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3"]),
        ("N1 none node (現行の判定で PASS)", "08.bump/run_0056_pl_nolim", []),
        ("N2 none cell (現行の判定で PASS)", "40.nozzle_design_tool/run_0085_forge_euler_cell_stgA", []),
        ("N3 none (旧書式) cell (現行の判定で PASS)", "35.uniform_periodic_box/run_0014_cell_keep_cbd_pure", []),
        ("N4 none node (case/46)", "46.sern_design/run_0100_species_regress_euler_lumped", []),
        ("N5 none node --segment", "44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX", ["--segment"]),
        ("N6 none node --from-floor", "16.nozzle_wys/run_0478_tp5_node_euler_profile_xr_cfl6_relax07", ["--from-floor", "16.nozzle_wys/run_0476_tp5_node_euler_profile_base"]),
    ]

    def staged_copy(src, tag):
        """--segment は run ディレクトリに residual_history_segment.csv を書くので、区間の判定に要るファイルだけ写す。"""
        dst = tmp / "seg" / tag
        dst.mkdir(parents=True)
        for pat in ("stage_manifest.json", "forge_launches.jsonl", "residual_history*.csv", "solverConfig.yaml"):
            for f in glob.glob(str(src / pat)):
                if os.path.basename(f) != "residual_history_segment.csv":
                    shutil.copy2(f, dst)
        return dst

    def verdict_lines(out):
        return [ln for ln in out.splitlines() if " -> " in ln or ln.startswith(("OVERALL", "REFUSED")) or "判定不能" in ln]

    kinds = {"sst": 0, "none": 0}
    for label, rel, extra in CLI:
        src = find_run(rel)
        refs = [x for x in extra if x != "--segment"]
        ref_dirs = [find_run(x) for x in refs[1:]] if refs else []
        if src is None or any(r is None for r in ref_dirs):
            stop(f"(g7) case/{rel} {refs} の residual_history.csv が見つからない")
        res = {}
        for name, tool in (("old", old_tools / "check_convergence.py"), ("new", TOOLS / "check_convergence.py")):
            tgt = staged_copy(src, f"{name}_{src.name}") if "--segment" in extra else src
            argv = [sys.executable, str(tool), str(tgt)] + (["--segment"] if "--segment" in extra else []) \
                + (["--from-floor", str(ref_dirs[0])] if refs else [])
            r = subprocess.run(argv, capture_output=True, text=True)
            out = (r.stdout + r.stderr).replace(str(tgt), "<RUN>")
            res[name] = (r.returncode, verdict_lines(out))
        kinds[label.split()[1]] = kinds.get(label.split()[1], 0) + 1
        check(f"(g7) {label} (case/{rel}{' ' + ' '.join(refs) if refs else ''}): 変更前後で VERDICT 行・終了コードが一致 "
              f"[{(res['new'][1] or ['?'])[0].strip()[:90]}]", res["old"] == res["new"], (res["old"], res["new"]))
    check(f"(g7) SST {kinds['sst']} 本・乱流なし {kinds['none']} 本 (各 3 本以上)", kinds["sst"] >= 3 and kinds["none"] >= 3)

print(f"\n{n_checks - fails}/{n_checks} PASS" + ("" if not fails else f"  ({fails} FAIL)"))
print(f"(作業ディレクトリ {tmp})")
sys.exit(1 if fails else 0)
