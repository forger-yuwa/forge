#!/usr/bin/env python3
"""学習の採否に床事象の判定を接続する試験 (plan tooling-sern-te-wake-grid §5.1 R1、codex result レビュー 2026-10-11 M1)。
forge は回さない (design/.venv-opt の Python で実行)。

合格 (測る前に固定):
  (1) キャンペーンの床判定の要件は基準 YAML の evaluate.floor_events / floor_events_window_steps から作る:
      floor_events 1 → 0 で要件が消え、窓 10000 → 1 で要件が変わる (値の検査は sern_gates と同じ)
  (2) 要件のあるキャンペーンは、各作動点の保存された判定が PASS・区間長 ≥ 要件・検査した更新数 ≥ 要件の行だけを採る:
      欠測 (要約なし・None・区間なし・INDETERMINATE)・短い区間・FAIL の 3 種の混入を拒否し、PASS (要件以上の区間) は採る。
      評価時にゲートが任意 (required: false) で FAIL の行 (レビューの再現) も拒否する。片方の作動点だけの欠陥でも拒否
  (3) 要件なし (floor_events 0・キー無し) のキャンペーンは従来どおり (床の要約を見ない)
  (4) 実際のゲート (sern_gates.evaluate_gates → 作動点要約 → JSON の往復) の PASS / FAIL / 判定不能 (段階起動の 4000 step に
      10000 の窓 = レビューの再現) が (2) のとおりに選別される
  (5) 評価方式の承認 (eval_method_carry_over) は床の要件を緩めない。Pareto の要約に要件を残す

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_floor_learn_tests.py
"""
import csv
import json
import shutil
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np
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

from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
import test_gate_bad_input as TGB  # noqa: E402  (合成の floor_events.csv)
import sern_acceptance_fixture as AF  # noqa: E402  受理検査の版と残差の入力の要約 (§5.1 R2k)

fails = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_floor_learn_"))


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


# ======================================================================== (1) 要件の作り方
R10k = D.floor_requirement_of({"floor_events": 1, "floor_events_window_steps": 10000})
check("(1) floor_events 1・窓 10000 → 必須・区間長 ≥ 10000・更新数 ≥ 10000",
      R10k == {"required": True, "window_steps": 10000, "min_updates": 10000}, R10k)
check("(1) floor_events 1 → 0 で要件が消える (窓のキーが残っていても)",
      D.floor_requirement_of({"floor_events": 0, "floor_events_window_steps": 10000}) is None)
R1 = D.floor_requirement_of({"floor_events": 1, "floor_events_window_steps": 1})
check("(1) 窓 10000 → 1 で要件が変わる", R1 is not None and R1 != R10k and R1["window_steps"] == 1, R1)
check("(1) キー無し・evaluate なしは要件なし", D.floor_requirement_of({}) is None and D.floor_requirement_of(None) is None)
RNW = D.floor_requirement_of({"floor_events": 1})
check("(1) 窓の指定なし: 必須・区間は空でないこと (長さ・更新数の下限 1)",
      RNW == {"required": True, "window_steps": None, "min_updates": 1}, RNW)
for label, ev in (("floor_events 2", {"floor_events": 2}), ("窓 0", {"floor_events": 1, "floor_events_window_steps": 0}),
                  ("窓 文字列", {"floor_events": 1, "floor_events_window_steps": "10000"})):
    check(f"(1) 不正な値は止まる: {label}", raises(lambda: D.floor_requirement_of(ev)))

# ======================================================================== (2)(3) 作動点要約の選別
ID = "a" * 64
ON = "L_b=1.0;curve=x"
OPS = ("m6_on", "m10_on", "m4_off")


def fe_(verdict="PASS", window=(10000, 20000), n=None, required=True):
    """sern_gates.floor_events_summary の形 (JSON の往復後と同じ: window は list)。"""
    if verdict is None:
        return None
    d = {"required": required, "verdict": verdict, "fail_class": None if verdict == "PASS" else "FLOOR_EVENT",
         "window": None if window is None else list(window),
         "n_updates_checked": (n if n is not None else (0 if window is None else window[1] - window[0])),
         "n_events": {"T": 0, "Rho": 0, "P": 0}, "first_events": [], "reasons": []}
    return d


def op_(fe="pass", mid=ID, te="off"):
    # 受理検査の版と残差の入力の要約は現行のゲートで判定した形 (2026-10-12 から学習の採否が一律に要求する、§5.1 R2k。この試験は床の要件だけを差にする)
    o = {"scalar_gradient_effective": "lsq", "slau_wall_normal_chi_effective": 1, D.EVAL_METHOD_ID: mid, D.TE_WAKE_EFFECTIVE: te,
         **AF.accepted_op_fields()}
    if fe != "absent":
        o["floor_events"] = fe_() if fe == "pass" else fe
    return o


def row_(ops, x=0.0):
    return {"status": "PASS", "flag_policy": D.FLAG_POLICY, "x": [x], "C_T_w": 0.9, "L_ramp": 1.0, "C_M_w": 0.0, "ops": ops}


def all_ops(**kw):
    return {k: op_(**kw) for k in OPS}


class Fake:
    """キャンペーンの偽物 (選別に要る属性だけ)。floor=… で床の要件、te/req で評価方式の要求。"""

    def __init__(self, rows, floor=None, te="off", req=None, co=()):
        self.rows, self.floor_required = rows, floor
        self.te_wake_required, self.eval_method_required, self.eval_method_carry_over = te, req, tuple(co)
        self.ref, self.ops = (-0.9, 20.0), []


def n_learn(rows, floor=None, te="off", req=None, co=()):
    return len(D.SernCampaign._XF(Fake(rows, floor, te, req, co))[0])


CASES = (
    # (ラベル, 行, 要件ありの期待, 要件なしの期待)
    ("PASS (区間 10000・更新 10000)", [row_(all_ops())], 1, 1),
    ("PASS (より長い区間 20000)", [row_(all_ops(fe=fe_(window=(20000, 40000))))], 1, 1),
    ("欠測: 要約のキーが無い (識別の導入前・ゲート導入前の行)", [row_(all_ops(fe="absent"))], 0, 1),
    ("欠測: 要約が None", [row_(all_ops(fe=None))], 0, 1),
    ("欠測: 判定不能 INDETERMINATE (区間なし・記録の欠け)", [row_(all_ops(fe=fe_("INDETERMINATE", window=None)))], 0, 1),
    ("欠測: PASS でも区間が読めない (None)", [row_(all_ops(fe={**fe_(), "window": None}))], 0, 1),
    ("欠測: PASS でも検査した更新数が無い", [row_(all_ops(fe={**fe_(), "n_updates_checked": None}))], 0, 1),
    ("欠測: 区間が整数でない (小数・真偽値)", [row_(all_ops(fe={**fe_(), "window": [10000.0, 20000.0]})),
                                       row_(all_ops(fe={**fe_(), "window": [True, 20000]}), 0.1)], 0, 2),
    ("短い区間: 段階起動の本段 4000 step の末尾 2000", [row_(all_ops(fe=fe_(window=(2000, 4000))))], 0, 1),
    ("短い区間: 区間長は足りても検査した更新数が足りない", [row_(all_ops(fe=fe_(window=(10000, 20000), n=9999)))], 0, 1),
    ("FAIL: ゲートが任意 (required: false) で FAIL (レビューの再現)", [row_(all_ops(fe=fe_("FAIL", required=False)))], 0, 1),
    ("FAIL: 必須で FAIL", [row_(all_ops(fe=fe_("FAIL")))], 0, 1),
    ("片方の作動点だけ FAIL", [row_({"m6_on": op_(), "m10_on": op_(fe=fe_("FAIL", required=False)), "m4_off": op_()})], 0, 1),
    ("片方の作動点だけ欠測", [row_({"m6_on": op_(), "m10_on": op_(fe="absent"), "m4_off": op_()})], 0, 1),
    ("PASS 2 行と FAIL 1 行: PASS の 2 行だけ", [row_(all_ops(), 0.1), row_(all_ops(), 0.2),
                                          row_(all_ops(fe=fe_("FAIL", required=False)), 0.3)], 2, 3),
)
for label, rows, want_req, want_free in CASES:
    got = n_learn(rows, R10k)
    check(f"(2) 要件あり (窓 10000): {label} → 採用 {got} (期待 {want_req})", got == want_req)
    got = n_learn(rows, None)
    check(f"(3) 要件なし: {label} → 採用 {got} (期待 {want_free}、従来どおり)", got == want_free)
check("(2) 窓 1 の要件なら短い区間 (2000) の PASS も採る (要件が窓に追従する)",
      n_learn([row_(all_ops(fe=fe_(window=(2000, 4000))))], R1) == 1)
check("(2) 窓の指定なしの要件: 空でない区間の PASS は採る・区間 0 は採らない",
      n_learn([row_(all_ops(fe=fe_(window=(2000, 4000))))], RNW) == 1
      and n_learn([row_(all_ops(fe=fe_(window=(4000, 4000))))], RNW) == 0)
# レビューの再現: 現行方式 id・既存の必須フラグを満たす PASS 行に floor_events {required: false, verdict: FAIL}
rev = row_(all_ops(fe=fe_("FAIL", required=False), te=ON))
check("(2) レビューの再現: 評価方式の要求を満たす行でも床 FAIL なら _learnable は False (要件ありのキャンペーン)",
      D._learnable(rev, method=(ON, ID, ())) is True and D._learnable(rev, method=(ON, ID, ()), floor=R10k) is False)

# ======================================================================== (4) 実際のゲート → 作動点要約 → JSON の往復
def make_run(d: Path, steps):
    d.mkdir(parents=True)
    n = 40
    with h5py.File(d / f"res_{steps}.h5", "w") as f:
        # SST の場 (roK・roOmega も。2026-10-12 から SST では必須: plan tooling-sern-te-wake-grid §5.1 R2g)
        for k, val in (("ro", 1.0), ("roUx", 2.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", 3.0e5), ("P", 1.0e5), ("T", 300.0),
                       ("roK", 1.0), ("roOmega", 1.0e3)):
            f.create_dataset(f"VALUE/{k}", data=np.full(n, val))
    # 残差は流体 5 列 + SST の 2 列、run の solverConfig.yaml (SST) も置く: 旧 fixture の 3 列 (rms_ro・rms_roUx・rms_roe) で設定なしの
    # CSV は check_convergence の判定不能 (必須列の欠落) で、2026-10-12 から受理ゲートも判定不能にする (plan tooling-sern-te-wake-grid §5.1 R2h1)
    # forge と同じ時系列 (outer_end の step 0..N−1 を 1 行ずつ、最終場 res_N、config に nStepOuter N と physProp) にする: 2026-10-12 から
    # 受理ゲートは残差の時系列の欠測・終了位置の不一致 (旧 fixture は 500 step おき・step N まで) と、化学種の列を決める physProp の
    # 無い config を判定不能にする (plan tooling-sern-te-wake-grid §5.1 R2i・R2j)
    (d / "solverConfig.yaml").write_text(f'physProp: {{thermalMethod: 0}}\nturbulence: {{model: "sst"}}\ntime:\n  last: {{nStepOuter: {steps}}}\n')
    with open(d / "residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"])
        for s in range(steps):
            v = 10 ** (-1 - 1.5 * min(s / steps * 3, 1))
            w.writerow([s, -1, "outer_end", v, v * 10, v * 10, v * 10, v * 100, v * 10, v * 1000])
    (d / "run_case_stdout.log").write_text("[run_case] forge exit=0\n")
    return d


hist = [{"step": 500 * (i + 1), "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0} for i in range(40)]


def summary_after_json(run, required, window):
    g = G.evaluate_gates(run, hist, 0, require_floor_events=required, floor_window_steps=window, model="sst")
    out = {"gates": g, "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0, "forge_rc": 0,
           "slau_wall_normal_chi_effective": 1, "scalar_gradient_effective": "lsq", D.EVAL_METHOD_ID: ID, D.TE_WAKE_EFFECTIVE: "off", D.ZTE_EFFECTIVE: "off"}
    s = D.SernCampaign._op_summary(out, run)
    return json.loads(json.dumps(s)), g


r_ok = make_run(tmp / "g_ok", 20000); TGB.write_floor_csv(str(r_ok), 20000, sessions=[2000, 2000, 2000, 4000])
s_ok, g_ok = summary_after_json(r_ok, True, 10000)
check("(4) 実ゲート: 判定区間の段 (20000 step) の最後の session に窓 10000 → PASS、要約は JSON の往復後も採る",
      g_ok["floor_events"]["verdict"] == "PASS" and s_ok["floor_events"]["window"] == [10000, 20000]
      and D._floor_op_ok(s_ok, R10k), s_ok["floor_events"])
r_st = make_run(tmp / "g_staged", 4000); TGB.write_floor_csv(str(r_st), 4000, sessions=[2000, 2000, 2000])
s_st, g_st = summary_after_json(r_st, False, 10000)
check("(4) 実ゲート: 段階起動の本段 4000 step に窓 10000 (レビューの再現) → 判定不能、ゲートが任意で総合 PASS でも採らない",
      g_st["verdict"] == "PASS" and g_st["floor_events"]["verdict"] == "INDETERMINATE" and not D._floor_op_ok(s_st, R10k)
      and n_learn([row_({k: s_st for k in OPS})], R10k) == 0 and n_learn([row_({k: s_st for k in OPS})], None) == 1)
r_ev = make_run(tmp / "g_event", 20000); TGB.write_floor_csv(str(r_ev), 20000, events={15000: (1, 0, 0)})
s_ev, g_ev = summary_after_json(r_ev, False, 10000)
check("(4) 実ゲート: 区間内に 1 件 (ゲート任意で総合 PASS) → 要約 FAIL、採らない",
      g_ev["verdict"] == "PASS" and s_ev["floor_events"]["verdict"] == "FAIL" and not D._floor_op_ok(s_ev, R10k))
r_no = make_run(tmp / "g_none", 20000)
s_no, _ = summary_after_json(r_no, False, 10000)
check("(4) 実ゲート: 記録なし → 判定不能、採らない", s_no["floor_events"]["verdict"] == "INDETERMINATE" and not D._floor_op_ok(s_no, R10k))
s_short, _ = summary_after_json(r_ok, True, 5000)
check("(4) 実ゲート: 窓 5000 で判定した PASS は窓 10000 の要件では採らない (短窓)", s_short["floor_events"]["verdict"] == "PASS"
      and not D._floor_op_ok(s_short, R10k) and D._floor_op_ok(s_short, D.floor_requirement_of({"floor_events": 1, "floor_events_window_steps": 5000})))

# ======================================================================== (5) 承認は床の要件を緩めない・要約・実際のキャンペーン
A6 = {"op": "m6_on", "from_eval_method": "legacy", "check_6": "PASS", "check_7": "PASS", "evidence": "plan §6 #6・#7"}
leg = {k: v for k, v in op_(fe="absent").items() if k not in (D.EVAL_METHOD_ID, D.TE_WAKE_EFFECTIVE)}
rows_co = [row_({"m6_on": leg, "m10_on": op_(te=ON), "m4_off": op_(te=ON)})]
check("(5) 承認済みの持ち越し (m6_on の旧評価) は評価方式の要求では採るが、床の要件ありでは欠測として外す",
      n_learn(rows_co, None, ON, ID, (A6,)) == 1 and n_learn(rows_co, R10k, ON, ID, (A6,)) == 0)

BASE = yaml.safe_load(open(ROOT / "case" / "46.sern_design" / "problem_moo_sst_node_3op.yaml"))


def campaign(name, ev_extra, rows=()):
    d = tmp / name; d.mkdir()
    raw = json.loads(json.dumps(BASE))
    ev = raw.setdefault("evaluate", {})
    for k in ("floor_events", "floor_events_window_steps"):
        ev.pop(k, None)
    ev.update(ev_extra)
    (d / "problem.yaml").write_text(yaml.safe_dump(raw, allow_unicode=True))
    if rows:
        (d / "camp").mkdir()
        (d / "camp" / "ledger.jsonl").write_text("\n".join(json.dumps(dict(r, tag=f"t{i}")) for i, r in enumerate(rows)) + "\n")
    return D.SernCampaign(d / "problem.yaml", d / "camp")


opnames = [o["name"] for o in BASE["spec"]["operating_points"]]
led = [row_({k: op_() for k in opnames}, 0.1), row_({k: op_(fe=fe_("FAIL", required=False)) for k in opnames}, 0.2),
       row_({k: op_(fe="absent") for k in opnames}, 0.3), row_({k: op_(fe=fe_(window=(2000, 4000))) for k in opnames}, 0.4)]
c_on = campaign("c_on", {"floor_events": 1, "floor_events_window_steps": 10000}, led)
c_off = campaign("c_off", {"floor_events": 0, "floor_events_window_steps": 10000}, led)
c_w1 = campaign("c_w1", {"floor_events": 1, "floor_events_window_steps": 1}, led)
c_nokey = campaign("c_nokey", {}, led)
check("(1) 実際のキャンペーン: 基準 YAML の floor_events 1 → 0・窓 10000 → 1 で要件が変わる",
      c_on.floor_required == R10k and c_off.floor_required is None and c_w1.floor_required == R1 and c_nokey.floor_required is None)
s_on, s_off, s_w1, s_nk = (c.summary() for c in (c_on, c_off, c_w1, c_nokey))
check("(2)(3) 実際のキャンペーン (台帳 4 行: PASS・FAIL・欠測・短窓): 要件 10000 は 1 行、窓 1 は 2 行 (短窓も可)、要件なしは 4 行",
      (s_on["n_pass"], s_w1["n_pass"], s_off["n_pass"], s_nk["n_pass"]) == (1, 2, 4, 4),
      (s_on["n_pass"], s_w1["n_pass"], s_off["n_pass"], s_nk["n_pass"]))
check("(5) Pareto の要約に要件を残し、外した行は台帳に残る (n_eval 4)",
      s_on["required_floor_events"] == R10k and s_off["required_floor_events"] is None and s_on["n_eval"] == 4
      and s_on["n_pass_excluded_by_policy"] == 3)
check("(3) 既存の試験の偽物 (床の属性なし) は要件なしとして動く", D._floor_policy(object()) is None
      and D._learnable(row_(all_ops(fe="absent"))))

shutil.rmtree(tmp, ignore_errors=True)
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
