#!/usr/bin/env python3
"""受理検査の版と残差の入力の検査済み要約を学習の採否に (plan tooling-sern-te-wake-grid §5.1 R2k、codex result 再レビュー 3 回目
2026-10-12 `notes/reviews/2026-10-12-tooling-sern-te-wake-grid-result-2.md` の M1)。forge は回さない (0 step)。

合格 (測る前に固定):
  (a) 記録: evaluate_gates は gates に受理検査の版 (sern_gates.ACCEPTANCE_VERSION) を書き、residual_health は要求した残差列・欠けた列と
      時系列の要約を返す。2D・3D の collect (力の履歴だけ差し替え、ゲートは本物) の metrics.json の gates に版と残差の要約が残る。
      _op_summary は版と残差の入力の検査済み要約 (residual_evidence) を作動点要約に載せ、JSON の往復で変わらない。
      gates に版が無い (旧い判定) なら作動点要約の版は None (現行の版を付けない)
  (b) codex の表: R2c (run_1100) の 1 回目の生の残差 (無ければ fixture を展開) の 7 種 (完全・rms_roY0 欠落・化学種 2 列欠落・先頭の 1 行だけ・
      末尾の 1 行だけ・途中の 1 step の欠測・1 step の重複) を、残差以外の証拠を 2 回目の保存値に固定して**旧実装 (9e5f0ffd) の sern_gates**
      で判定 → 7 種とも GATES PASS (codex の再現) → 台帳行 (R2k 前 [d2f3ad68] の _op_summary = 保存済みの台帳の形と、現行の _op_summary の
      2 通り) → R2k 前の _learnable はキャンペーンの要件 (評価方式・床事象・判定区間) をすべて渡して 7 種とも採る (codex の再現)。
      現行の _learnable・学習配列 (_XF)・Pareto (summary) は 7 種とも外す
  (c) 同じ入力を現行ゲートで再判定: 完全は PASS で、その行だけが学習配列と Pareto に入る。不完全な 6 種は FAIL/JUDGMENT_UNVERIFIABLE で、
      status を PASS に書き換えた行も (版は現行でも) 残差の入力の要約の不合格で外す
  (d) 版の無い旧行・版の古い行・要約の欠けた行 (要約なし・入力の欄の欠け・必須列の欠け・欠けた列あり・時系列の欄の欠けと不一致・最終場の不一致)
      は外す。作動点の 1 つだけの欠陥でも外す。旧評価の持ち越しの承認 (評価方式・処置) は緩めない
  (e) 一律: 床事象・判定区間を要求しない 2D のキャンペーン (problem_moo_sst_node_3op.yaml) でも版・要約の無い行を外す (台帳には残る)。
      既存の再判定経路 (SernCampaign.rejudge、collect は現行ゲートを通す偽物) で判定し直した行は版と要約を持ち、完全な残差の行だけが
      pareto_rejudged の学習の母集団に戻る。Pareto の要約に要求する版を残す

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_learn_acceptance_version_tests.py [--r2c-raw DIR | --fixture]
"""
import argparse
import csv
import gzip
import hashlib
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

from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

DEFAULT_RAW = "/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad/r2c_raw"
RAW_NAME, RAW_SHA256 = "residual_judgment1.csv", "1679d4c919d83e211e095a833d3bd6664c76be1c78ffe47fbf08786f1f2c7155"   # 1 回目 (12,964,240 B)
ap = argparse.ArgumentParser()
ap.add_argument("--r2c-raw", default=os.environ.get("SERN_R2C_RAW", DEFAULT_RAW),
                help="R2c の生の残差 (residual_judgment1.csv = 1 回目) を置いたディレクトリ")
ap.add_argument("--fixture", action="store_true", help="生の CSV なしに、間引いた fixture (展開して使う) だけで回す")
args = ap.parse_args()

fails = 0
n_checks = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_learn_acceptance_"))
C46 = ROOT / "case" / "46.sern_design"
R2C = ROOT / "notes" / "investigations" / "2026-10-12-r2c" / "run_1100_r2c_lsw10_m6_judgment"
FIXTURE = R2C / "residual_judgment1_stride40.csv.gz"
PROD10 = C46 / "problem_3d_prod_3op_wallres_lswx10.yaml"
PRE_GATES = "9e5f0ffd"     # R2i・R2j の前 (codex が旧判定に使った実装)
PRE_DRIVER = "d2f3ad68"    # R2k の前 (driver_sern は 9e5f0ffd から変わっていない)
N_J = 20000                # 判定区間の段の nStepOuter


def check(name, ok, detail=""):
    global fails, n_checks
    n_checks += 1
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def stop(msg):
    print("止める: " + msg)
    sys.exit(2)


def old_module(rev, path, name, package):
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{rev}:{path}"], capture_output=True, text=True)
    if r.returncode != 0:
        stop(f"変更前の {path} ({rev}) を git から読めない: {r.stderr.strip()[:200]}")
    m = types.ModuleType(name)
    m.__package__ = package
    m.__file__ = str(ROOT / path)
    exec(compile(r.stdout, f"{rev}:{path}", "exec"), m.__dict__)
    return m


def jrt(x):
    """台帳 (ledger.jsonl) の JSON の往復。"""
    return json.loads(json.dumps(x, default=str))


# --- 旧実装: sern_gates (9e5f0ffd、その時点の check_convergence を使わせる) と driver_sern (R2k の前) -----------------------------
OLD_CC = old_module(PRE_GATES, "solver_density_cuda/tools/check_convergence.py", "_check_convergence_pre_r2i", None)
OLD_G = old_module(PRE_GATES, "design/forge_design/metrics/sern_gates.py", "forge_design.metrics._sern_gates_pre_r2i", "forge_design.metrics")
OLD_G._cc = OLD_CC
OLD_D = old_module(PRE_DRIVER, "design/forge_design/opt/driver_sern.py", "forge_design.opt._driver_sern_pre_r2k", "forge_design.opt")
if hasattr(OLD_G, "residual_series_check") or hasattr(OLD_G, "ACCEPTANCE_VERSION") or hasattr(OLD_D, "_acceptance_row_ok"):
    stop("旧実装のモジュールが変更前の形でない")


def load_hist(path):
    rows = list(csv.DictReader(open(path)))
    return [{k: (int(v) if k == "step" else float(v)) for k, v in r.items() if v not in ("", None)} for r in rows]


# --- 残差以外の証拠: R2c の 2 回目 (PASS) の保存値 --------------------------------------------------------------------------
P10 = load_problem(PROD10)
JP = R2.judgment_plan(P10.evaluate)
M2 = json.loads((R2C / "metrics.json").read_text())
CTX = {"saved": M2["gates"], "hist": load_hist(R2C / "force_history.csv"), "mon": G.read_monitor_csv(R2C / "TE_MONITOR.csv")}
check("前提: R2c の 2 回目の保存された判定は PASS で、受理検査の版を持たない (R2k の前の記録)",
      CTX["saved"]["verdict"] == "PASS" and G.ACCEPTANCE_VERSION_KEY not in CTX["saved"])
# キャンペーンの要件 (codex の再現と同じく評価方式・床事象・判定区間をすべて渡す): 生産 YAML (lswx10) と run_1100 の評価方式の id
TE, ID = M2[D.TE_WAKE_EFFECTIVE], M2[D.EVAL_METHOD_ID]
FR, JR = D.floor_requirement_of(P10.evaluate), D.judgment_requirement_of(P10.evaluate)
REQ = {"method": (TE, ID, ()), "floor": FR, "judgment": JR}
check("前提: キャンペーンの要件 (評価方式の id・床事象 10000・判定区間 10000) がそろう",
      isinstance(ID, str) and len(ID) == 64 and TE != "off" and FR and FR["window_steps"] == 10000 and JR and JR["window_steps"] == 10000)
PROV = {k: M2.get(k) for k in ("slau_wall_normal_chi_effective", "scalar_gradient_effective", D.ZTE_EFFECTIVE, "zero_thickness_edge_velocity",
                               D.EVAL_METHOD_ID, D.TE_WAKE_EFFECTIVE, "eval_method", "mesh_provenance", "flag_policy")}


def gates(mod, d):
    """残差 (residual_health・residual_scale_gate) だけを実際に評価し、場・最終場の床・床事象は 2 回目の保存値に固定、力の履歴と
    判定区間の §6 #2・#3 は 2 回目の保存系列から計算して evaluate_gates を通す (run_sern_residual_series_tests と同じ固定の仕方)。"""
    jg = mod.judgment_gate(CTX["hist"], JP["window_steps"], JP["required"], JP["planned_steps"], monitor_rows=CTX["mon"])
    vals = {"field_health": CTX["saved"]["field"], "floor_gate": CTX["saved"]["floors"], "floor_event_gate": CTX["saved"]["floor_events"]}
    orig = {n: getattr(mod, n) for n in vals}
    for n, v in vals.items():
        setattr(mod, n, lambda *a, _v=v, **kw: _v)
    try:
        return mod.evaluate_gates(d, CTX["hist"], 0, p_min=float(P10.evaluate.get("p_min", 1.0)), require_floor_events=True,
                                  floor_window_steps=JP["window_steps"], judgment=jg, model="sst")
    finally:
        for n, v in orig.items():
            setattr(mod, n, v)


def mkrun(name, hdr, body, drop=()):
    d = tmp / name
    d.mkdir()
    (d / "solverConfig.yaml").write_text((R2C / "solverConfig_judgment.yaml").read_text())
    keep = [i for i, c in enumerate(hdr) if c not in drop]
    with open(d / "residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([hdr[i] for i in keep])
        for r in body:
            w.writerow([r[i] for i in keep])
    (d / f"res_{N_J}.h5").write_bytes(b"")
    return d


def expand_stride(body, n_steps=N_J):
    """40 行おきに間引いた fixture を、各行の値を次の行の step の手前まで保持して outer_end の step 0..N−1 に戻す
    (run_sern_residual_series_tests と同じ。試験だけの仕掛け)。"""
    out = []
    for k, r in enumerate(body):
        s1 = int(body[k + 1][0]) if k + 1 < len(body) else n_steps
        out += [[str(s)] + r[1:] for s in range(int(r[0]), s1)]
    return out


def make_out(g):
    """collect の戻り (metrics.json) の形: ゲートだけを差し替え、力と来歴は 2 回目の保存値。"""
    return {**PROV, "gates": g, "forge_rc": 0, "step": M2.get("step"),
            **{k: M2[k] for k in ("C_T", "C_T_with_shear", "C_L", "C_M")}}


def row_(ops, tag, x, ct=0.9, L=1.0):
    return {"tag": tag, "status": "PASS", "fail_class": None, "flag_policy": D.FLAG_POLICY, "x": [float(x)] * 5, "C_T_w": ct,
            "L_ramp": L, "C_M_w": M2["C_M"], "degraded": False, "degraded_ops": [], "ops": ops}


class Fake:
    """キャンペーンの偽物 (選別と Pareto の要約に要る属性だけ)。要件は codex の再現と同じ (評価方式・床事象・判定区間)。"""

    def __init__(self, rows, method=(TE, ID, ()), floor=FR, judgment=JR, zte="off", zco=()):
        self.rows, self.ref, self.ops = rows, (-0.90, 20.0), []
        self.te_wake_required, self.eval_method_required, self.eval_method_carry_over = method[0], method[1], tuple(method[2])
        self.floor_required, self.judgment_required, self.zte_required, self.zte_carry_over = floor, judgment, zte, tuple(zco)


# --- 残差の入力: 生の CSV (無ければ展開した fixture) ------------------------------------------------------------------------
if not args.fixture:
    p = Path(args.r2c_raw) / RAW_NAME
    if not p.exists():
        stop(f"R2c の生の残差 {p} が無い (12 MB の CSV はリポジトリに置かない)。--r2c-raw で場所を渡すか、--fixture で fixture だけで回す")
    if hashlib.sha256(p.read_bytes()).hexdigest() != RAW_SHA256:
        stop(f"{p} の sha256 が登録した R2c の残差と違う")
    rows = list(csv.reader(open(p)))
    HDR, BODY = rows[0], rows[1:]
    SRC = "生の CSV (1 回目、60000 行)"
else:
    with gzip.open(FIXTURE, "rt", newline="") as fh:
        fx = list(csv.reader(fh))
    HDR, BODY = fx[0], expand_stride(fx[1:])
    SRC = f"fixture ({FIXTURE.name} を step-hold で {len(BODY)} 行に展開)"
print(f"    残差の入力: {SRC}")
OE = [r for r in BODY if r[2] == "outer_end"]


def steps_of(rows, lo, hi):
    return [r for r in rows if lo <= int(r[0]) < hi]


VARIANTS = [
    ("完全", {}),
    ("rms_roY0 欠落", {"drop": ("rms_roY0",)}),
    ("化学種 2 列欠落", {"drop": ("rms_roY0", "rms_roY1")}),
    ("先頭の 1 行だけ", {"body": [OE[0]]}),
    ("末尾の 1 行だけ", {"body": [OE[-1]]}),
    ("途中の 1 step の欠測", {"body": steps_of(BODY, 0, 10000) + steps_of(BODY, 10001, N_J)}),
    ("1 step の重複", {"body": steps_of(BODY, 0, 10001) + steps_of(BODY, 10000, 10001) + steps_of(BODY, 10001, N_J)}),
]
SC = {}
for i, (lab, kw) in enumerate(VARIANTS):
    d = mkrun(f"v{i}", HDR, kw.get("body", BODY), drop=kw.get("drop", ()))
    g_old, g_new = gates(OLD_G, d), gates(G, d)
    SC[lab] = {"d": d, "g_old": g_old, "g_new": g_new,
               "s_old_pre": jrt(OLD_D.SernCampaign._op_summary(make_out(g_old), d)),     # 保存済みの台帳の行の形 (R2k の前)
               "s_old_cur": jrt(D.SernCampaign._op_summary(make_out(g_old), d)),          # 旧い判定を現行の _op_summary に通す
               "s_new": jrt(D.SernCampaign._op_summary(make_out(g_new), d))}             # 現行ゲートで再判定

# ======================================================================== (a) 記録
s_full = SC["完全"]["s_new"]
g_full = SC["完全"]["g_new"]
check("(a) evaluate_gates は gates に現行の受理検査の版を書く",
      g_full.get(G.ACCEPTANCE_VERSION_KEY) == G.ACCEPTANCE_VERSION and isinstance(G.ACCEPTANCE_VERSION, str) and G.ACCEPTANCE_VERSION)
rh = g_full["residual"]
want_cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
check("(a) residual_health は要求した残差列 (流体 5・SST 2・化学種 2) と欠けた列 (なし)・時系列の要約を返す",
      sorted(rh.get("required_columns") or []) == sorted(want_cols) and rh.get("missing_columns") == []
      and rh["series"]["n_steps"] == N_J and rh["series"]["n_outer_end"] == N_J, (rh.get("required_columns"), rh.get("missing_columns")))
ev = s_full.get("residual_evidence") or {}
check("(a) _op_summary: 版と残差の入力の検査済み要約 (入力の妥当性・必須列・時系列) が作動点要約に載る (JSON の往復後)",
      s_full.get(G.ACCEPTANCE_VERSION_KEY) == G.ACCEPTANCE_VERSION and ev.get("input_ok") is True and ev.get("unverifiable") is False
      and sorted(ev.get("required_columns") or []) == sorted(want_cols) and ev.get("missing_columns") == []
      and ev.get("series", {}).get("final_field") == f"res_{N_J}.h5" and D._acceptance_op_ok(s_full), ev)
ev_d = SC["rms_roY0 欠落"]["s_new"]["residual_evidence"]
check("(a) 現行ゲートで rms_roY0 が欠けた残差: 要約に入力の不合格と欠けた列が残る",
      ev_d["input_ok"] is False and ev_d["unverifiable"] is True and ev_d["missing_columns"] == ["rms_roY0"], ev_d)
check("(a) 旧い判定 (9e5f0ffd) を現行の _op_summary に通しても現行の版は付かない (None)・要約の時系列と必須列は None",
      all(SC[lab]["s_old_cur"].get(G.ACCEPTANCE_VERSION_KEY) is None
          and SC[lab]["s_old_cur"]["residual_evidence"]["series"] is None
          and SC[lab]["s_old_cur"]["residual_evidence"]["required_columns"] is None for lab, _ in VARIANTS))
check("(a) 保存済みの台帳の形 (R2k 前の _op_summary) には版も要約も無い",
      all(G.ACCEPTANCE_VERSION_KEY not in SC[lab]["s_old_pre"] and "residual_evidence" not in SC[lab]["s_old_pre"] for lab, _ in VARIANTS))


# 2D・3D の collect: 力の履歴だけを差し替え、ゲート (evaluate_gates) は本物。metrics.json の gates に版と残差の要約が残る
def synth_run(d: Path, steps=400, model="sst", drop=(), truncate=0):
    """forge と同じ時系列の小さな run: solverConfig (nStepOuter・乱流モデル・physProp)・最終場 res_N (有限・正の場)・残差 (outer_end 0..N−1)。"""
    d.mkdir(parents=True, exist_ok=True)
    n = 40
    fields = [("ro", 1.0), ("roUx", 2.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", 3.0e5), ("P", 1.0e5), ("T", 300.0)]
    if model == "sst":
        fields += [("roK", 1.0), ("roOmega", 1.0e3)]
    with h5py.File(d / f"res_{steps}.h5", "w") as f:
        for k, val in fields:
            f.create_dataset(f"VALUE/{k}", data=np.full(n, val))
    turb = "sst" if model == "sst" else "none"
    (d / "solverConfig.yaml").write_text(f'physProp: {{thermalMethod: 0}}\nturbulence: {{model: "{turb}"}}\ntime:\n  last: {{nStepOuter: {steps}}}\n')
    cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe"] + (["rms_roK", "rms_roOmega"] if model == "sst" else [])
    cols = [c for c in cols if c not in drop]
    with open(d / "residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "inner", "phase"] + cols)
        for s in range(steps - truncate):
            v = 10 ** (-1 - 1.5 * min(s / steps * 3, 1))
            w.writerow([s, -1, "outer_end"] + [v * (1 + i) for i in range(len(cols))])
    (d / "run_case_stdout.log").write_text("[run_case] forge exit=0\n")
    return d


HIST = [{"step": 10 * (i + 1), "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0} for i in range(40)]


def real_collect(mod, name):
    d = synth_run(tmp / name, model="euler")
    prob = d / "problem.yaml"
    prob.write_text(yaml.safe_dump({"name": "t", "type": "sern_2d", "gas": {}, "spec": {}, "dv": {}, "geometry": {}, "mesh": {}, "evaluate": {}}))
    info = {"states": {"exhaust": {"ro": 1.0, "u": 1.0, "P": 1.0}, "ext": {"P": 1.0}}, "H_m": 1.0, "F_ideal_N_per_m": 1.0,
            "design": {"L_ramp": 1.0}, "moc_forces": {"C_T": 0.0, "C_L": 0.0, "C_M": 0.0}, "discretization": "node", "half_W_m": 1.0}
    (d / "prepare_info.json").write_text(json.dumps(info))
    hist3 = [{k: v for k, v in h.items() if k != "C_T_with_shear"} for h in HIST]
    saved = (R2.force_history, R2.write_force_history_csv, R3.force_parts_plan, R3._part_steps, R3.forces3d, R3.write_force_history_csv)
    R2.force_history = lambda *a, **k: [dict(h, C_T_wall=0.0, T_wall=0.0, L=0.0, M_noseup=0.0) for h in hist3]
    R3.force_parts_plan = lambda *a, **k: {"required": ["ramp"], "friction": [], "faces_source": "試験", "reasons": []}
    R3._part_steps = lambda *a, **k: [h["step"] for h in hist3]
    R3.forces3d = lambda run_dir, s, *a, **k: dict(next(h for h in hist3 if h["step"] == s))
    R2.write_force_history_csv = R3.write_force_history_csv = lambda *a, **k: None
    try:
        mod.collect(prob, d, rc=0)
    finally:
        (R2.force_history, R2.write_force_history_csv, R3.force_parts_plan, R3._part_steps, R3.forces3d, R3.write_force_history_csv) = saved
    return json.loads((d / "metrics.json").read_text())


for mod, dim in ((R2, "2D"), (R3, "3D")):
    try:
        m = real_collect(mod, f"collect_{dim}")
        gm = m["gates"]
        ok = (gm.get(G.ACCEPTANCE_VERSION_KEY) == G.ACCEPTANCE_VERSION and gm["verdict"] == "PASS"
              and gm["residual"]["input"]["ok"] is True and gm["residual"]["missing_columns"] == []
              and gm["residual"]["series"]["n_outer_end"] == 400 and D._acceptance_op_ok(jrt(D.SernCampaign._op_summary(m, "run"))))
        check(f"(a) {dim} の collect: metrics.json の gates に受理検査の版・残差の入力 (input)・必須列・時系列が残り、台帳の要約が合格",
              ok, (gm.get(G.ACCEPTANCE_VERSION_KEY), gm["verdict"], gm["reasons"][:3]))
    except Exception as e:  # noqa: BLE001
        check(f"(a) {dim} の collect: metrics.json の gates に受理検査の版", False, f"{type(e).__name__}: {e}")

# ======================================================================== (b) codex の表 (旧判定の PASS 行)
print("\n    codex の表 (残差以外の証拠 = R2c の 2 回目の保存値、キャンペーンの要件 = 評価方式・床事象・判定区間):")
print(f"      {'残差の入力':22s} {'旧ゲート':10s} {'現行ゲート':28s} {'R2k 前の学習':12s} {'現行の学習 (旧行)':16s} {'現行の学習 (再判定)':16s}")
for lab, _ in VARIANTS:
    sc = SC[lab]
    r_pre, r_cur, r_new = row_({"m6_on": sc["s_old_pre"]}, "x", 0), row_({"m6_on": sc["s_old_cur"]}, "x", 0), row_({"m6_on": sc["s_new"]}, "x", 0)
    sc["learn"] = (OLD_D._learnable(r_pre, **REQ), D._learnable(r_pre, **REQ), D._learnable(r_cur, **REQ), D._learnable(r_new, **REQ))
    gn = sc["g_new"]
    print(f"      {lab:22s} {sc['g_old']['verdict']:10s} {(gn['verdict'] + '/' + str(gn['fail_class'] or '')):28s} "
          f"{'採用' if sc['learn'][0] else '除外':12s} {'採用' if (sc['learn'][1] or sc['learn'][2]) else '除外':16s} "
          f"{'採用' if sc['learn'][3] else '除外':16s}")
for lab, _ in VARIANTS:
    sc = SC[lab]
    check(f"(b) {lab}: 旧実装 (9e5f0ffd) のゲートは GATES PASS、R2k 前の _learnable は要件をすべて渡しても採る (codex の再現)",
          sc["g_old"]["verdict"] == "PASS" and sc["learn"][0] is True, (sc["g_old"]["verdict"], sc["g_old"]["fail_class"], sc["learn"]))
    check(f"(b) {lab}: 旧判定の PASS 行は現行の _learnable で外す (保存済みの台帳の形・現行の _op_summary の両方)",
          sc["learn"][1] is False and sc["learn"][2] is False, sc["learn"])

# ======================================================================== (c) 現行ゲートで再判定した行
check("(c) 完全な残差を現行ゲートで再判定: GATES PASS・学習に採る",
      SC["完全"]["g_new"]["verdict"] == "PASS" and SC["完全"]["learn"][3] is True, (SC["完全"]["g_new"]["fail_class"], SC["完全"]["learn"]))
for lab, _ in VARIANTS[1:]:
    sc = SC[lab]
    check(f"(c) {lab}: 現行ゲートは FAIL/JUDGMENT_UNVERIFIABLE、status を PASS に書き換えた行も残差の要約の不合格で外す (版は現行)",
          (sc["g_new"]["verdict"], sc["g_new"]["fail_class"]) == ("FAIL", "JUDGMENT_UNVERIFIABLE") and sc["learn"][3] is False
          and sc["s_new"].get(G.ACCEPTANCE_VERSION_KEY) == G.ACCEPTANCE_VERSION and not D._acceptance_op_ok(sc["s_new"]),
          (sc["g_new"]["verdict"], sc["g_new"]["fail_class"], sc["learn"]))

# 学習配列と Pareto: 旧判定の行 (2 つの形 × 7 種)・現行ゲートで不合格の書き換え行 (6 種) を優位な位置 (高い C_T_w・短い L_ramp) に置き、
# 再判定した完全な行だけを劣位に置く。外れていれば学習配列は 1 行・Pareto はその 1 点だけ
rows = []
for i, (lab, _) in enumerate(VARIANTS):
    rows.append(row_({"m6_on": SC[lab]["s_old_pre"]}, f"old_pre_{i}", 0.10 + 0.01 * i, ct=0.95 + 0.001 * i, L=1.0 - 0.01 * i))
    rows.append(row_({"m6_on": SC[lab]["s_old_cur"]}, f"old_cur_{i}", 0.20 + 0.01 * i, ct=0.95 + 0.001 * i, L=1.0 - 0.01 * i))
    if i:
        rows.append(row_({"m6_on": SC[lab]["s_new"]}, f"new_bad_{i}", 0.30 + 0.01 * i, ct=0.95 + 0.001 * i, L=1.0 - 0.01 * i))
rows.append(row_({"m6_on": SC["完全"]["s_new"]}, "rejudged", 0.5, ct=0.90, L=2.0))
rows = jrt(rows)
fk = Fake(rows)
X, F = D.SernCampaign._XF(fk)
check("(c) 学習配列 (_XF): 現行ゲートで再判定した完全な行だけ (1 行)",
      X.shape[0] == 1 and np.allclose(X[0], [0.5] * 5) and np.allclose(F[0], [-0.90, 2.0]), X.tolist())
sm = D.SernCampaign.summary(fk)
check("(c) Pareto (summary): 母集団はその 1 行・Pareto はその 1 点、要求する版を要約に残す",
      sm["n_pass"] == 1 and [p_["tag"] for p_ in sm["pareto"]] == ["rejudged"] and sm["n_pass_excluded_by_policy"] == len(rows) - 1
      and sm["required_acceptance_version"] == G.ACCEPTANCE_VERSION
      and sm["pareto"][0]["ops"]["m6_on"][G.ACCEPTANCE_VERSION_KEY] == G.ACCEPTANCE_VERSION,
      (sm["n_pass"], [p_["tag"] for p_ in sm["pareto"]]))
Xo, _ = OLD_D.SernCampaign._XF(fk)
smo = OLD_D.SernCampaign.summary(fk)
check(f"(c) 対照: R2k 前の学習配列・Pareto は旧判定の行 (14)・書き換え行 (6) まで採り ({Xo.shape[0]} 行)、Pareto は旧判定の行になる",
      Xo.shape[0] == len(rows) and "rejudged" not in [p_["tag"] for p_ in smo["pareto"]]
      and any(p_["tag"].startswith("old_") for p_ in smo["pareto"]), (Xo.shape[0], [p_["tag"] for p_ in smo["pareto"]]))

# ======================================================================== (d) 版・要約の欠け
OK = SC["完全"]["s_new"]
EV = OK["residual_evidence"]


def with_ev(**kw):
    return {**OK, "residual_evidence": {**EV, **kw}}


def with_series(**kw):
    return with_ev(series={**EV["series"], **kw})


def drop_keys(o, *keys):
    return {k: v for k, v in o.items() if k not in keys}


N = EV["series"]["n_steps"]
D_CASES = [
    ("対照: 現行の版と完全な要約", OK, True),
    ("版の無い旧行 (版・要約のキーが無い)", drop_keys(OK, G.ACCEPTANCE_VERSION_KEY, "residual_evidence"), False),
    ("版のキーが無い (要約はある)", drop_keys(OK, G.ACCEPTANCE_VERSION_KEY), False),
    ("版が None", {**OK, G.ACCEPTANCE_VERSION_KEY: None}, False),
    ("版の古い行 (sern-acceptance-v0)", {**OK, G.ACCEPTANCE_VERSION_KEY: "sern-acceptance-v0"}, False),
    ("要約のキーが無い", drop_keys(OK, "residual_evidence"), False),
    ("要約が None", {**OK, "residual_evidence": None}, False),
    ("入力の妥当性の欄が None (旧い residual_health)", with_ev(input_ok=None), False),
    ("入力が妥当でない", with_ev(input_ok=False), False),
    ("判定不能の欄が None", with_ev(unverifiable=None), False),
    ("判定不能", with_ev(unverifiable=True), False),
    ("要求した列が None", with_ev(required_columns=None), False),
    ("要求した列が空", with_ev(required_columns=[]), False),
    ("要求した列に流体の列 (rms_roe) が無い", with_ev(required_columns=[c for c in EV["required_columns"] if c != "rms_roe"]), False),
    ("欠けた列がある", with_ev(missing_columns=["rms_roY0"]), False),
    ("欠けた列の欄が None", with_ev(missing_columns=None), False),
    ("時系列の要約が None", with_ev(series=None), False),
    ("時系列: outer_end が 1 行足りない", with_series(n_outer_end=N - 1), False),
    ("時系列: 先頭が step 1", with_series(first=1), False),
    ("時系列: 末尾が N−2", with_series(last=N - 2), False),
    ("時系列: 最終場が res_{N−1}", with_series(final_field=f"res_{N - 1}.h5"), False),
    ("時系列: 最終場の欄が None", with_series(final_field=None), False),
    ("時系列: nStepOuter が None", with_series(n_steps=None), False),
    ("時系列: nStepOuter が真偽値", with_series(n_steps=True, n_outer_end=1, first=0, last=0, final_field="res_1.h5"), False),
]
for lab, o, want in D_CASES:
    got = D._learnable(jrt(row_({"m6_on": o}, "d", 0)), **REQ)
    check(f"(d) {lab}: {'採る' if want else '外す'}", got is want, got)
OPS3 = ("m6_on", "m10_on", "m4_off")
check("(d) 3 作動点のうち 1 つだけ版が無い行は外す (他の 2 つは現行)",
      D._learnable(jrt(row_({k: OK for k in OPS3}, "d", 0)), **REQ) is True
      and D._learnable(jrt(row_({"m6_on": OK, "m10_on": drop_keys(OK, G.ACCEPTANCE_VERSION_KEY), "m4_off": OK}, "d", 0)), **REQ) is False)
check("(d) 3 作動点のうち 1 つだけ要約の時系列が不完全な行は外す",
      D._learnable(jrt(row_({"m6_on": OK, "m10_on": OK, "m4_off": with_series(n_outer_end=1)}, "d", 0)), **REQ) is False)
# 旧評価の持ち越しの承認 (評価方式・処置) は受理検査の版の要件を緩めない
A6 = {"op": "m6_on", "from_eval_method": "legacy", "check_6": "PASS", "check_7": "PASS", "evidence": "試験"}
leg_ok = drop_keys(OK, D.EVAL_METHOD_ID, D.TE_WAKE_EFFECTIVE)
leg_old = drop_keys(leg_ok, G.ACCEPTANCE_VERSION_KEY, "residual_evidence")
check("(d) 評価方式の持ち越しの承認 [m6_on, legacy]: 版と要約を持つ旧方式の行は採り、版の無い旧行は承認があっても外す",
      D._learnable(jrt(row_({"m6_on": leg_ok, "m10_on": OK}, "d", 0)), method=(TE, ID, (A6,)), floor=FR, judgment=JR) is True
      and D._learnable(jrt(row_({"m6_on": leg_old, "m10_on": OK}, "d", 0)), method=(TE, ID, (A6,)), floor=FR, judgment=JR) is False)
z_old = {**drop_keys(OK, G.ACCEPTANCE_VERSION_KEY), D.ZTE_EFFECTIVE: "off"}
check("(d) 処置の持ち越し指定 [m6_on]: 版の無い行は持ち越し指定があっても外す",
      len(D.SernCampaign._XF(Fake([row_({"m6_on": z_old}, "z", 0)], zte="tags=cowl_in;rings=2", zco=("m6_on",)))[0]) == 0
      and len(D.SernCampaign._XF(Fake([row_({"m6_on": {**OK, D.ZTE_EFFECTIVE: "off"}}, "z", 0)], zte="tags=cowl_in;rings=2", zco=("m6_on",)))[0]) == 1)
check("(d) 要件を渡さない呼び出し (評価方式・床・判定区間なし) でも版・要約は要求する (一律)",
      D._learnable(jrt(row_({"m6_on": drop_keys(OK, G.ACCEPTANCE_VERSION_KEY, D.EVAL_METHOD_ID, D.TE_WAKE_EFFECTIVE)}, "d", 0))) is False
      and D._learnable(jrt(row_({"m6_on": drop_keys(OK, D.EVAL_METHOD_ID, D.TE_WAKE_EFFECTIVE)}, "d", 0))) is True)

# ======================================================================== (e) 床事象・判定区間を要求しない 2D のキャンペーンと再判定の経路
BASE2 = yaml.safe_load(open(C46 / "problem_moo_sst_node_3op.yaml"))
OPS2D = [o["name"] for o in BASE2["spec"]["operating_points"]]


def campaign_dir(name, rows):
    d = tmp / name
    (d / "camp").mkdir(parents=True)
    raw = json.loads(json.dumps(BASE2))
    (d / "problem.yaml").write_text(yaml.safe_dump(raw, allow_unicode=True))
    (d / "camp" / "ledger.jsonl").write_text("\n".join(json.dumps(r, ensure_ascii=False) for r in rows) + "\n")
    return d


def old_op():
    """R2k の前の台帳の作動点要約 (処置以外の採否は満たす。版・要約なし)。"""
    return {"scalar_gradient_effective": "lsq", "slau_wall_normal_chi_effective": 1, "gate": "PASS", "residual": "NOT CONVERGED (stalled/plateau)"}


def acc_op():
    return {**old_op(), G.ACCEPTANCE_VERSION_KEY: G.ACCEPTANCE_VERSION, "residual_evidence": EV}


cd = campaign_dir("c2d", [row_({k: old_op() for k in OPS2D}, "t_old", 0.1, ct=0.95, L=1.0),
                          row_({k: acc_op() for k in OPS2D}, "t_cur", 0.2, ct=0.90, L=2.0)])
c2 = D.SernCampaign(cd / "problem.yaml", cd / "camp")
s2 = c2.summary()
check("(e) 2D のキャンペーン (床事象・判定区間の要件なし): 版の無い旧行は外し、現行の版と要約の行だけ採る (台帳には 2 行とも残る)",
      c2.floor_required is None and c2.judgment_required is None and s2["n_pass"] == 1 and [p_["tag"] for p_ in s2["pareto"]] == ["t_cur"]
      and s2["n_eval"] == 2 and len((cd / "camp" / "ledger.jsonl").read_text().splitlines()) == 2
      and s2["required_acceptance_version"] == G.ACCEPTANCE_VERSION, (s2["n_pass"], s2["n_eval"]))

# 再判定の経路: 台帳の 2 評価 (旧行。残差が完全な評価と末尾 3 step が切り詰められた評価) の run を、collect を現行ゲートに通す偽物で再判定する
rj = campaign_dir("rejudge", [row_({k: old_op() for k in OPS2D}, "t_full", 0.3), row_({k: old_op() for k in OPS2D}, "t_trunc", 0.4)])
camp = rj / "camp"
for tag in ("t_full", "t_trunc"):
    (camp / f"{tag}.yaml").write_text((rj / "problem.yaml").read_text())
    for op in OPS2D:
        rd = synth_run(camp / f"{tag}_{op}", truncate=(3 if tag == "t_trunc" else 0))
        (rd / "prepare_info.json").write_text(json.dumps({"design": {"L_ramp": 1.0, "warnings": []}}))


def fake_collect(prob, rd, out_dir=None, rc=None, require_residual_pass=False):
    g = G.evaluate_gates(rd, HIST, 0, require_residual_pass=require_residual_pass, model="sst")
    out = {"gates": g, "forge_rc": 0, "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0, "step": HIST[-1]["step"],
           "scalar_gradient_effective": "lsq", "slau_wall_normal_chi_effective": 1, "flag_policy": D.FLAG_POLICY, D.ZTE_EFFECTIVE: "off"}
    if out_dir is not None:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        (Path(out_dir) / "metrics.json").write_text(json.dumps(out, default=str))
    return out


c_rj = D.SernCampaign(rj / "problem.yaml", camp)
before = c_rj.summary()
_R = D.R
D.R = types.SimpleNamespace(collect=fake_collect)
try:
    after = c_rj.rejudge(tmp / "rejudge_out")
finally:
    D.R = _R
new_rows = [json.loads(line) for line in (tmp / "rejudge_out" / "ledger_rejudged.jsonl").read_text().splitlines() if line.strip()]
by_tag = {r["tag"]: r for r in new_rows}
check("(e) 再判定の前: 版の無い旧行 2 行はどちらも学習・Pareto に入らない", before["n_pass"] == 0 and before["n_eval"] == 2)
check("(e) 再判定 (--rejudge の経路): 完全な残差の評価は現行ゲートで PASS・版と要約を持ち、学習の母集団 (pareto_rejudged) に戻る",
      by_tag["t_full"]["status"] == "PASS" and all(o.get(G.ACCEPTANCE_VERSION_KEY) == G.ACCEPTANCE_VERSION and D._acceptance_op_ok(o)
                                                    for o in by_tag["t_full"]["ops"].values())
      and after["n_pass"] == 1 and [p_["tag"] for p_ in after["pareto"]] == ["t_full"],
      (by_tag["t_full"]["status"], by_tag["t_full"].get("note"), "base", D._learnable_base(by_tag["t_full"]),
       "acceptance", D._acceptance_row_ok(by_tag["t_full"]), "zte", D._zte_row_ok(by_tag["t_full"]), after["n_pass"]))
check("(e) 再判定: 末尾 3 step が切り詰められた残差の評価は FAIL/JUDGMENT_UNVERIFIABLE で戻らない",
      by_tag["t_trunc"]["status"] == "FAIL" and by_tag["t_trunc"]["fail_class"] == "JUDGMENT_UNVERIFIABLE", (by_tag["t_trunc"]["status"], by_tag["t_trunc"]["fail_class"]))
pj = json.loads((tmp / "rejudge_out" / "pareto_rejudged.json").read_text())
check("(e) 再判定の結果の Pareto の要約に要求する版が残り、元の台帳は書き換えない (2 行のまま・版なし)",
      pj["required_acceptance_version"] == G.ACCEPTANCE_VERSION
      and all(G.ACCEPTANCE_VERSION_KEY not in o for line in (camp / "ledger.jsonl").read_text().splitlines() if line.strip()
              for o in json.loads(line)["ops"].values()))

print(f"\n{n_checks - fails}/{n_checks} PASS" + ("" if not fails else f"  ({fails} FAIL)"))
print(f"(作業ディレクトリ {tmp})")
shutil.rmtree(tmp, ignore_errors=True)
sys.exit(1 if fails else 0)
