#!/usr/bin/env python3
"""判定区間の経路の欠測・未達を拒否する試験 (plan tooling-sern-te-wake-grid §5.1 R2a・R2b・R2d・R2e、
codex diagnose 2026-10-12 [notes/reviews/2026-10-12-sern-judgment-path-diagnose.md] の M1・M2・Minor・穴 (e))。forge は回さない。

合格 (測る前に固定):
  R2a 必須量と標本の格子を設定から固定
    (a1) diagnose の判別 A/B を受入れ試験に: run_1096 (L_sw 1.0・m10_on) の保存済みの履歴 A は PASS、同じ履歴から
         C_T_with_shear の列だけを消した B は判定不能 (JUDGMENT_UNVERIFIABLE)。必須量は問題設定 (生産 YAML の model: sst) から
    (a2) 標本の格子: 各 10000 step 窓の後半 5000 step だけに 20 標本を置いた入力 (標本数は 20/20)・step の重複・終了 step の
         不一致 (不足・超過) を判定不能にする。温度の監視も同じ予定 step 列と照合する (末尾窓の重複・予定外)
    (a3) Euler 指定 (model: euler) では 3 量で判定できる (B の履歴が PASS)。必須量・予定 step 列の作り方の検査
         (未知のモデル・N < 2W・割り切れない間隔・窓あたり 20 未満・outStepInterval の明示なしは止まる)
    (a4) collect (3D) が問題設定から作った必須量・予定 step 列を judgment_gate に渡す
  R2b 判定区間の判定を学習の採否に接続
    (b1) run_1097 (L_sw 1.0・m6_on) の保存記録 (旧 GATES PASS・床 PASS) から作った行は、方式 ID をキャンペーンの要求に合わせても、
         判定区間を要求するキャンペーンでは拒否する: 判定の要約なし・judgment: UNSETTLED の明示・保存履歴の再判定
         (UNSETTLED、C_M の d 1.97e-3 > 5e-4) のどれでも。判定 PASS の行 (run_1096 の再判定) は採る
    (b2) 要約の欠け・食い違い (要約なしの PASS・要約の窓が短い・必須量の欠け [Euler で判定した行を SST のキャンペーンで]・
         温度の領域の欠け・片方の作動点だけ未達) を拒否。旧評価の持ち越しの承認 (eval_method_carry_over) では解除しない
    (b3) 判定区間を要求しないキャンペーンは従来どおり (run_1097 の行も採る)。実際のキャンペーンの要件と Pareto の要約
  R2d 段の移行の検査
    (d1) 必須の保存量の集合を config から: 生産 config (TP・SST・2 種) は流体 5 + roK・roOmega + roY0・roY1、Euler CPG は 5、
         トレーサ付きは roXi も。未知の乱流モデルは止まる
    (d2) 実際の restart_by_index (CPG・SST) の index コピーは通す。転送先に roK・roOmega が無い・−0 を +0 にした
         (np.array_equal は等しいと言う)・dtype 違い・shape 違い・継続元に無い量を拒否
    (d3) run_judgment が config の必須量で検査する: roK・roOmega を写さない継続は止まる
  R2e 起動の環境の記録
    (e1) 結果に効く環境変数の実際の値が runner_launches.jsonl (run_forge・監視つきの起動) と JUDGMENT.json (全体・各段) に残る。
         未設定は null + unset に名前。その他の FORGE_* も残す。runner は FORGE_ALLOW_UNVERIFIED_SPECIES を既定で付けない

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_judgment_path_tests.py
"""
import csv
import json
import os
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

import stage_manifest as SM  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402
import sern_acceptance_fixture as AF  # noqa: E402  受理検査の版と残差の入力の要約 (§5.1 R2k)

fails = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_judgment_path_"))
C46 = ROOT / "case" / "46.sern_design"
PROD10 = C46 / "problem_3d_prod_3op_wallres_lswx10.yaml"
LSW10 = ROOT / "notes" / "investigations" / "2026-10-11-lsw10"
R1096 = LSW10 / "run_1096_tewake_lsw10_m10_cont20k"      # case/46.sern_design/run_1096_… の保存済みの写し (AWS の run の記録)
R1097 = LSW10 / "run_1097_tewake_lsw10_m6_cont20k"


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


for d in (R1096, R1097):
    for f in ("force_history.csv", "TE_MONITOR.csv", "metrics.json"):
        if not (d / f).exists():
            print(f"止める: 保存済みの履歴 {d / f} が無い (判別 A/B の入力)")
            sys.exit(2)


def load_hist(path, drop=()):
    rows = list(csv.DictReader(open(path)))
    return [{k: (int(v) if k == "step" else float(v)) for k, v in r.items() if v not in ("", None) and k not in drop} for r in rows]


# ======================================================================== R2a
P10 = load_problem(PROD10)
JP = R2.judgment_plan(P10.evaluate)          # 生産 YAML (model: sst、N 20000・W 10000・outStepInterval 500)
check("(a3) 生産 YAML の予定: SST の 4 量・予定 step 列 500..20000 の 40 標本・窓 10000",
      tuple(JP["required"]) == ("C_T", "C_T_with_shear", "C_L", "C_M") and JP["planned_steps"] == list(range(500, 20001, 500))
      and JP["window_steps"] == 10000 and JP["out_interval"] == 500, {k: JP[k] for k in ("required", "window_steps", "out_interval")})
mon96 = G.read_monitor_csv(R1096 / "TE_MONITOR.csv")
hA = load_hist(R1096 / "force_history.csv")
hB = load_hist(R1096 / "force_history.csv", drop=("C_T_with_shear",))
gA = G.judgment_gate(hA, JP["window_steps"], JP["required"], JP["planned_steps"], monitor_rows=mon96)
gB = G.judgment_gate(hB, JP["window_steps"], JP["required"], JP["planned_steps"], monitor_rows=mon96)
check("(a1) 判別 A/B の A: run_1096 の保存済みの履歴 (SST 4 量) は PASS", gA["verdict"] == "PASS" and gA["ok"]
      and sorted(gA["force_window"]["quantities"]) == sorted(JP["required"]), gA["reasons"])
check("(a1) 判別 A/B の B: 同じ履歴から C_T_with_shear の列だけを消すと判定不能 (3 量に減らして PASS しない)",
      gB["verdict"] == "INDETERMINATE" and gB["fail_class"] == "JUDGMENT_UNVERIFIABLE" and not gB["extendable"]
      and any("C_T_with_shear" in r and "列が力の履歴に無い" in r for r in gB["reasons"]), gB["reasons"])
check("(a1) B の温度の判定は A と同じ (力の列だけが原因)", gB["temperature"]["verdict"] == gA["temperature"]["verdict"] == "PASS")
EULER = G.judgment_required_quantities("euler")
gBe = G.judgment_gate(hB, 10000, EULER, JP["planned_steps"], monitor_rows=mon96)
check("(a3) Euler 指定 (model: euler) では 3 量で判定できる: B の履歴が PASS (C_T_with_shear を見ない)",
      gBe["verdict"] == "PASS" and sorted(gBe["force_window"]["quantities"]) == sorted(EULER), gBe["reasons"])


def synth(steps, **kw):
    return [{"step": s, "C_T": 0.9, "C_T_with_shear": 0.88, "C_L": 0.07, "C_M": -1.8, **kw} for s in steps]


def synth_mon(steps):
    rows = []
    for s in steps:
        r = {"step": float(s), "nonfinite_T": 0.0}
        for R in G.TE_REGIONS:
            r.update({f"{R}_n": 100.0, f"{R}_Tmin": 220.0, f"{R}_Tmin_node": 7.0, f"{R}_Tmin_x": 0.1, f"{R}_Tmin_y": 0.0,
                      f"{R}_Tmin_z": 0.0, f"{R}_n_lt150": 0.0, f"{R}_n_lt200": 0.0, f"{R}_vol_lt200": 0.0})
        rows.append(r)
    return rows


PL = JP["planned_steps"]
SST4 = JP["required"]
ok0 = G.judgment_gate(synth(PL), 10000, SST4, PL, monitor_rows=synth_mon(PL))
check("(a2) 合成: 予定どおりの 40 標本 (定常) は PASS", ok0["verdict"] == "PASS", ok0["reasons"])
# 標本の偏り: 各 10000 step 窓の後半 5000 step だけに 250 step 間隔で 20 標本 (標本数は 20/20 で数だけなら足りる)
biased = [s for w0 in (0, 10000) for s in range(w0 + 5250, w0 + 10001, 250)]
gbias = G.judgment_gate(synth(biased), 10000, SST4, PL, monitor_rows=synth_mon(PL))
fw = gbias["force_window"]
check("(a2) 標本の偏り (各窓の後半 5000 step だけに 20 標本): 判定不能 (予定の標本の欠け・予定外の step)",
      len(biased) == 40 and sum(s > 10000 for s in biased) == 20 and gbias["fail_class"] == "JUDGMENT_UNVERIFIABLE"
      and any("欠けている" in r for r in fw["reasons"]) and any("予定外" in r for r in fw["reasons"]), fw["reasons"])
dup = synth(PL) + synth([15000])
dup.sort(key=lambda h: h["step"])
gdup = G.judgment_gate(dup, 10000, SST4, PL, monitor_rows=synth_mon(PL))
check("(a2) step の重複 (15000 が 2 標本): 判定不能", gdup["fail_class"] == "JUDGMENT_UNVERIFIABLE"
      and any("重複" in r for r in gdup["force_window"]["reasons"]), gdup["force_window"]["reasons"])
short = synth(PL[:-1])
gshort = G.judgment_gate(short, 10000, SST4, PL, monitor_rows=synth_mon(PL))
check("(a2) 終了 step の不一致 (19500 で終わる = 最後の標本の欠け): 判定不能",
      gshort["fail_class"] == "JUDGMENT_UNVERIFIABLE" and any("終了 step" in r for r in gshort["force_window"]["reasons"]),
      gshort["force_window"]["reasons"])
over = synth(PL + [20500])
gover = G.judgment_gate(over, 10000, SST4, PL, monitor_rows=synth_mon(PL + [20500]))
check("(a2) 終了 step の不一致 (予定より長く 20500 まで): 判定不能 (力・温度とも)",
      gover["fail_class"] == "JUDGMENT_UNVERIFIABLE" and any("終了 step" in r for r in gover["force_window"]["reasons"])
      and gover["temperature"]["verdict"] == "INDETERMINATE", (gover["force_window"]["reasons"], gover["temperature"]["reasons"]))
gmdup = G.judgment_gate(synth(PL), 10000, SST4, PL, monitor_rows=sorted(synth_mon(PL) + synth_mon([15000]), key=lambda r: r["step"]))
check("(a2) 温度の監視の末尾窓に重複: 判定不能", gmdup["fail_class"] == "JUDGMENT_UNVERIFIABLE"
      and any("重複" in r for r in gmdup["temperature"]["reasons"]), gmdup["temperature"]["reasons"])
gmx = G.judgment_gate(synth(PL), 10000, SST4, PL, monitor_rows=sorted(synth_mon(PL) + synth_mon([15250]), key=lambda r: r["step"]))
check("(a2) 温度の監視の末尾窓に予定外の step: 判定不能", gmx["fail_class"] == "JUDGMENT_UNVERIFIABLE"
      and any("予定外" in r for r in gmx["temperature"]["reasons"]), gmx["temperature"]["reasons"])
check("(a2) 予定 step 列・必須量が無い・壊れた呼び出しは判定不能 (合格にしない)",
      G.judgment_gate(synth(PL), 10000, SST4, [], monitor_rows=synth_mon(PL))["fail_class"] == "JUDGMENT_UNVERIFIABLE"
      and G.judgment_gate(synth(PL), 10000, (), PL, monitor_rows=synth_mon(PL))["fail_class"] == "JUDGMENT_UNVERIFIABLE"
      and G.judgment_gate(synth(PL), 10000, ("C_T", "C_L"), PL, monitor_rows=synth_mon(PL))["fail_class"] == "JUDGMENT_UNVERIFIABLE"
      and G.judgment_gate(synth(PL), 10000, SST4, [1000, 500], monitor_rows=synth_mon(PL))["fail_class"] == "JUDGMENT_UNVERIFIABLE")

check("(a3) 必須量: model なし = euler (runner の既定と同じ)・大文字の SST も sst・未知のモデルは止まる",
      G.judgment_required_quantities(None) == EULER and G.judgment_required_quantities("SST") == SST4
      and raises(lambda: G.judgment_required_quantities("kw")))
for label, args in (("N < 2W", (19000, 10000, 500)), ("間隔が W を割り切らない", (20000, 10000, 300)),
                    ("窓あたり 10 標本", (20000, 10000, 1000)), ("間隔 0", (20000, 10000, 0)), ("真偽値", (20000, 10000, True))):
    check(f"(a3) 予定 step 列の作成は止まる: {label}", raises(lambda: G.judgment_planned_steps(*args)))
ev = dict(P10.evaluate)
check("(a3) judgment_plan: outStepInterval の明示が無い・小数・判定区間の段なしの扱い",
      raises(lambda: R2.judgment_plan({k: v for k, v in ev.items() if k != "outStepInterval"}))
      and raises(lambda: R2.judgment_plan({**ev, "outStepInterval": 500.0}))
      and R2.judgment_plan({k: v for k, v in ev.items() if k not in ("judgment_steps", "judgment_extend_once")}) is None)
check("(a3) judgment_plan: model: euler なら 3 量", tuple(R2.judgment_plan({**ev, "model": "euler"})["required"]) == EULER)

# (a4) collect (3D) の経路: 問題設定の必須量・予定 step 列が judgment_gate に渡る
rdc = tmp / "collect3d"; rdc.mkdir()
(rdc / "prepare_info.json").write_text(json.dumps({"states": {"exhaust": {"ro": 1.0, "u": 1000.0, "P": 1e5}, "ext": {"P": 1e3}}, "H_m": 0.1,
                                                   "F_ideal_N_per_m": 1e4, "half_W_m": 0.1, "discretization": "node", "moc_forces": {}}))
seen = []
_saved = (R3.judgment_gate, R3.evaluate_gates)
try:
    R3.judgment_gate = lambda hist, W, req, planned, **kw: (seen.append((W, tuple(req), list(planned))),
                                                            _saved[0](hist, W, req, planned, **kw))[1]
    R3.evaluate_gates = lambda *a, **k: {"steadiness": {"series": {}}, "objective": "C_T", "verdict": "FAIL", "fail_class": "NO_FORCES",
                                         "reasons": [], "judgment": k.get("judgment")}
    for yml, want in ((PROD10, SST4), (None, EULER)):
        if yml is None:
            y = yaml.safe_load(PROD10.read_text()); y["evaluate"]["model"] = "euler"
            yml = tmp / "p_euler.yaml"; yml.write_text(yaml.safe_dump(y, allow_unicode=True))
        out = R3.collect(yml, rdc, rc=0)
        check(f"(a4) collect (3D, {Path(yml).name}): judgment_gate に必須量 {len(want)} 量・予定 step 列 (500..20000)・窓 10000 を渡す",
              seen and seen[-1] == (10000, tuple(want), PL) and out["gates"]["judgment"]["fail_class"] == "JUDGMENT_UNVERIFIABLE", seen[-1:])
finally:
    R3.judgment_gate, R3.evaluate_gates = _saved

# ======================================================================== R2b
m96, m97 = json.loads((R1096 / "metrics.json").read_text()), json.loads((R1097 / "metrics.json").read_text())
check("(b1) run_1097 の保存記録: 旧 GATES PASS・床 PASS (10000, 20000]・判定区間の判定なし (段の導入前)",
      m97["gates"]["verdict"] == "PASS" and m97["floor_events"]["verdict"] == "PASS" and m97["floor_events"]["window"] == [10000, 20000]
      and "judgment" not in m97["gates"])
JR = D.judgment_requirement_of(P10.evaluate)
FR = D.floor_requirement_of(P10.evaluate)
check("(b3) 判定区間の要件 (生産 YAML): 窓 10000・SST 4 量・温度の 4 領域・標本 20、#1 の床の要件を含む (キャンペーンの床の要件と同じ)",
      JR == {"window_steps": 10000, "quantities": list(SST4), "regions": list(G.TE_REGIONS), "min_samples": 20, "floor": FR}
      and FR == {"required": True, "window_steps": 10000, "min_updates": 10000}, (JR, FR))
ID, TE = m97[R2.EVAL_METHOD_ID], m97[R2.TE_WAKE_EFFECTIVE]
check("(b1) run_1096・run_1097 の方式 ID は同じ (L_sw は方式に入らない、R3)", m96[R2.EVAL_METHOD_ID] == ID and ID)


def op_from(m, jg=None):
    """保存された metrics.json (collect の出力) から作動点要約を作る。jg を渡すと判定区間の判定を付け直す (= --rejudge と同じ形)。
    受理検査の版と残差の入力の要約は現行のゲートで再判定した形にそろえる (2026-10-12 から学習の採否は一律に現行の版を要求する
    [§5.1 R2k] が、保存された metrics は版の導入前。この試験は判定区間の要件だけを差にする)。"""
    out = dict(m)
    if jg is not None:
        out["gates"] = {**m["gates"], "judgment": jg}
    return json.loads(json.dumps({**D.SernCampaign._op_summary(out, "run"), **AF.accepted_op_fields()}))


def row_(ops, x=0.0):
    return {"status": "PASS", "flag_policy": D.FLAG_POLICY, "x": [x], "C_T_w": 0.9, "L_ramp": 1.0, "C_M_w": 0.0, "ops": ops}


class Fake:
    """キャンペーンの偽物 (選別に要る属性だけ)。"""

    def __init__(self, rows, judgment=JR, floor=FR, co=()):
        self.rows, self.judgment_required, self.floor_required = rows, judgment, floor
        self.te_wake_required, self.eval_method_required, self.eval_method_carry_over = TE, ID, tuple(co)
        self.ref, self.ops = (-0.9, 20.0), []


def n_learn(rows, **kw):
    return len(D.SernCampaign._XF(Fake(rows, **kw))[0])


j97 = G.judgment_gate(load_hist(R1097 / "force_history.csv"), 10000, SST4, PL, monitor_rows=G.read_monitor_csv(R1097 / "TE_MONITOR.csv"))
j96 = gA
q97 = j97["force_window"]["quantities"]["C_M"]
check("(b1) run_1097 の保存履歴の再判定: UNSETTLED、C_M の d 1.9722766e-3 > 5e-4 (diagnose の再計算と同じ)",
      j97["verdict"] == "UNSETTLED" and abs(q97["d"] - 1.9722766e-3) < 1e-9 and q97["d"] > q97["tol"] == 5e-4, (j97["verdict"], q97))
o97_none, o97_rj, o96_rj = op_from(m97), op_from(m97, j97), op_from(m96, j96)
o97_str = {**o97_none, "judgment": "UNSETTLED"}
check("(b1) 作動点要約: 保存記録のままは判定なし、再判定すると judgment と要約が付く",
      "judgment" not in o97_none and o97_rj["judgment"] == "UNSETTLED" and o97_rj["judgment_summary"]["verdict"] == "UNSETTLED"
      and o96_rj["judgment"] == "PASS" and o96_rj["judgment_summary"]["force_window"]["required"] == list(SST4))
check("(b1) 方式 ID・床・既存の必須フラグはキャンペーンの要求を満たす (判定区間の要件だけが差)",
      all(D._learnable(row_({"m6_on": o}), method=(TE, ID, ()), floor=FR) for o in (o97_none, o97_rj, o97_str, o96_rj)))
for label, o in (("判定の要約なし (旧 GATES PASS のまま)", o97_none), ("judgment: UNSETTLED を明示 (diagnose の再現)", o97_str),
                 ("保存履歴の再判定 UNSETTLED (C_M の d 超過)", o97_rj)):
    check(f"(b1) run_1097 の行 [{label}]: 方式 ID を合わせても判定区間を要求するキャンペーンでは採らない",
          n_learn([row_({"m6_on": o})]) == 0 and D._learnable(row_({"m6_on": o}), method=(TE, ID, ()), floor=FR, judgment=JR) is False)
check("(b1) 判定 PASS の行 (run_1096 の保存履歴を再判定) は採る", n_learn([row_({"m10_on": o96_rj})]) == 1)


def with_js(o, **upd):
    """要約の一部を書き換えた作動点要約 (upd: 'force_window.window' のような点区切りのキー)。"""
    o = json.loads(json.dumps(o))
    for k, v in upd.items():
        d = o["judgment_summary"]; *path, last = k.split(".")
        for p_ in path:
            d = d[p_]
        d[last] = v
    return o


fw96 = o96_rj["judgment_summary"]["force_window"]
CASES = (
    ("judgment PASS だが要約が無い", {**o96_rj, "judgment_summary": None}),
    ("judgment PASS だが要約の判定が UNSETTLED (食い違い)", with_js(o96_rj, verdict="UNSETTLED")),
    ("要約の窓の長さ 5000 < 10000", with_js(o96_rj, window_steps=5000)),
    ("力の窓の区間が 5000 step", with_js(o96_rj, **{"force_window.window": [15000, 20000]})),
    ("力の窓の標本 19", with_js(o96_rj, **{"force_window.n_prev": 19})),
    ("必須量 C_T_with_shear の判定が無い (Euler で判定した行)", with_js(o96_rj, **{"force_window.ok": {q: True for q in EULER}})),
    ("必須量の 1 つが不合格", with_js(o96_rj, **{"force_window.ok": {**fw96["ok"], "C_M": False}})),
    ("温度の領域 R_SE が無い", with_js(o96_rj, **{"temperature.regions": ["ALL", "RET", "R_TE"]})),
    ("温度の窓の判定が INDETERMINATE", with_js(o96_rj, **{"temperature.verdict": "INDETERMINATE"})),
)
for label, o in CASES:
    check(f"(b2) {label} → 採らない", n_learn([row_({"m10_on": o})]) == 0)
check("(b2) 片方の作動点だけ判定なし → 採らない", n_learn([row_({"m10_on": o96_rj, "m6_on": o97_none})]) == 0)
o96_floorx = {**o96_rj, "floor_events": {**o96_rj["floor_events"], "verdict": "INDETERMINATE"}}
check("(b2) 判定 PASS でも床の要件 (#1) を満たさなければ採らない (床の要約が判定不能)。キャンペーンの床の要件を外しても判定区間の要件が拒否する",
      n_learn([row_({"m10_on": o96_floorx})]) == 0 and n_learn([row_({"m10_on": o96_floorx})], floor=None) == 0
      and n_learn([row_({"m10_on": o96_rj})], floor=None) == 1)
# 旧評価の持ち越しの承認 (旧 ID の行を名指しで承認) は判定区間の要件を緩めない
OLD = "f" * 64
A6 = {"op": "m6_on", "from_eval_method": OLD, "check_6": "PASS", "check_7": "PASS", "evidence": "plan §6 #6・#7"}
o97_old = {**o97_none, R2.EVAL_METHOD_ID: OLD}
check("(b2) 旧 ID の承認済みの持ち越し: 判定区間の要件なしでは採る (承認の効果) が、要件ありでは判定なしとして外す",
      n_learn([row_({"m6_on": o97_old})], judgment=None, co=(A6,)) == 1 and n_learn([row_({"m6_on": o97_old})], co=(A6,)) == 0
      and n_learn([row_({"m6_on": {**o97_rj, R2.EVAL_METHOD_ID: OLD}})], co=(A6,)) == 0)
check("(b3) 判定区間を要求しないキャンペーンは従来どおり: run_1097 の行 (判定なし・UNSETTLED とも) も採る",
      n_learn([row_({"m6_on": o97_none})], judgment=None) == 1 and n_learn([row_({"m6_on": o97_rj})], judgment=None) == 1)
check("(b3) 既存の試験の偽物 (判定区間の属性なし) は要件なしとして動く", D._judgment_policy(object()) is None)

# 実際のキャンペーン (基準 YAML の evaluate に判定区間の段): 要件・台帳の選別・Pareto の要約
BASE = yaml.safe_load(open(C46 / "problem_moo_sst_node_3op.yaml"))
opnames = [o["name"] for o in BASE["spec"]["operating_points"]]


def campaign(name, ev_extra, rows=(), drop=()):
    d = tmp / name; d.mkdir()
    raw = json.loads(json.dumps(BASE))
    e = raw.setdefault("evaluate", {})
    for k in ("floor_events", "floor_events_window_steps", "judgment_steps", "judgment_extend_once") + tuple(drop):
        e.pop(k, None)
    e.update(ev_extra)
    (d / "problem.yaml").write_text(yaml.safe_dump(raw, allow_unicode=True))
    if rows:
        (d / "camp").mkdir()
        (d / "camp" / "ledger.jsonl").write_text("\n".join(json.dumps(dict(r, tag=f"t{i}")) for i, r in enumerate(rows)) + "\n")
    return D.SernCampaign(d / "problem.yaml", d / "camp")


def base_op(o):
    """キャンペーン (2D、格子変形の要求なし = off) の方式の要求に合わせた作動点要約 (方式のキーは識別の導入前と同じく無し)。"""
    return {k: v for k, v in o.items() if k not in (R2.EVAL_METHOD_ID, R2.TE_WAKE_EFFECTIVE)}


EVJ = {"floor_events": 1, "floor_events_window_steps": 10000, "judgment_steps": 20000, "judgment_extend_once": 1,
       "outStepInterval": 500, "model": "sst"}
led = [row_({k: base_op(o96_rj) for k in opnames}, 0.1), row_({k: base_op(o97_rj) for k in opnames}, 0.2),
       row_({k: base_op(o97_none) for k in opnames}, 0.3)]
c_j = campaign("c_j", EVJ, led)
c_nj = campaign("c_nj", {"floor_events": 1, "floor_events_window_steps": 10000}, led)
s_j, s_nj = c_j.summary(), c_nj.summary()
check("(b3) 実際のキャンペーン: 判定区間の段ありは要件が付き、台帳 3 行 (PASS・UNSETTLED・判定なし) のうち 1 行だけ採る。"
      "段なしは従来どおり 3 行", c_j.judgment_required == JR and c_nj.judgment_required is None
      and (s_j["n_pass"], s_nj["n_pass"]) == (1, 3) and s_j["required_judgment"] == JR and s_nj["required_judgment"] is None
      and s_j["n_eval"] == 3, (s_j["n_pass"], s_nj["n_pass"]))
check("(b3) 判定区間の段の設定が壊れた基準 YAML はキャンペーンを作れない (outStepInterval なし)",
      raises(lambda: campaign("c_bad", {k: v for k, v in EVJ.items() if k != "outStepInterval"}, drop=("outStepInterval",))))

# ======================================================================== R2d
p10m6 = load_problem(PROD10); R2.select_operating_point(p10m6, "m6_on")
cfg_prod = R3._solver_config(p10m6, 4000, 500, 0.25, 2851.0)
k_prod = R2.required_state_keys(cfg_prod)
check("(d1) 生産 config (TP・SST・2 種): 流体 5 + roK・roOmega + roY0・roY1", k_prod == ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"],
      k_prod)
cfg_eu = 'physProp: {thermalMethod: 0, viscMethod: 0}\nturbulence: {model: "none"}\n'
cfg_tr = cfg_prod.replace("thermoHrefTemp: 298.15", "thermoHrefTemp: 298.15, tracer: exhaust")
check("(d1) Euler CPG は 5 量、トレーサ付き TP は roXi も、未知の乱流モデルは止まる",
      R2.required_state_keys(cfg_eu) == ["ro", "roUx", "roUy", "roUz", "roe"] and R2.required_state_keys(cfg_tr)[-1] == "roXi"
      and raises(lambda: R2.required_state_keys(cfg_eu.replace('"none"', '"kw"'))))

SST7 = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"]
rng = np.random.default_rng(3)


def write_h5(path, keys, n=50, dtype=np.float64, neg0=None, scale=1.0):
    with h5py.File(path, "w") as f:
        for i, k in enumerate(keys):
            a = (rng.standard_normal(n) * scale + i).astype(dtype)
            if neg0 is not None and k == neg0:
                a[3] = -0.0
            f.create_dataset(f"VALUE/{k}", data=a)
        f.create_dataset("VALUE/P", data=np.ones(n))


def run_dir_cpg(name, dst_keys, **kw):
    rd = tmp / name; rd.mkdir()
    (rd / "solverConfig.yaml").write_text('physProp: {thermalMethod: 0}\nturbulence: {model: "sst"}\n')
    write_h5(rd / "res_100.h5", SST7, **kw)
    write_h5(rd / "sern.h5", dst_keys)
    return rd


rd = run_dir_cpg("d_ok", SST7, neg0="roUy")
R2.restart_by_index(rd / "res_100.h5", rd / "sern.h5")
check("(d2) 実際の restart_by_index (CPG・SST) の index コピーは通す (−0 も写る)", R2._restart_mismatch(rd / "res_100.h5", rd / "sern.h5", SST7) == [])
rd = run_dir_cpg("d_noturb", SST7[:5])
R2.restart_by_index(rd / "res_100.h5", rd / "sern.h5")
bad = R2._restart_mismatch(rd / "res_100.h5", rd / "sern.h5", SST7)
check("(d2) 転送先に roK・roOmega が無い (restart_by_index は写さずに通る): 拒否 (旧検査は共通の量だけ見て空だった)",
      bad == ["roK が継続先に無い", "roOmega が継続先に無い"], bad)
rd = run_dir_cpg("d_neg0", SST7, neg0="roUy")
R2.restart_by_index(rd / "res_100.h5", rd / "sern.h5")
with h5py.File(rd / "sern.h5", "r+") as f:
    a = f["VALUE/roUy"][:]; a = a + 0.0; f["VALUE/roUy"][:] = a          # −0 + 0 = +0 (符号だけを失った継続)
with h5py.File(rd / "res_100.h5", "r") as s, h5py.File(rd / "sern.h5", "r") as d_:
    eq = np.array_equal(s["VALUE/roUy"][:], d_["VALUE/roUy"][:]) and np.signbit(s["VALUE/roUy"][3]) and not np.signbit(d_["VALUE/roUy"][3])
bad = R2._restart_mismatch(rd / "res_100.h5", rd / "sern.h5", SST7)
check("(d2) −0 を +0 にした継続 (np.array_equal は等しいと言う): バイト比較で拒否", eq and bad == ["roUy のバイト列が違う"], (eq, bad))
rd = run_dir_cpg("d_dtype", SST7)
with h5py.File(rd / "sern.h5", "r+") as f:
    del f["VALUE/roe"]; f.create_dataset("VALUE/roe", data=np.zeros(50, np.float32))
with h5py.File(rd / "res_100.h5", "r+") as f:
    f["VALUE/roe"][:] = np.arange(50, dtype=np.float64)                       # float32 で厳密に表せる値 (数値は一致する)
R2.restart_by_index(rd / "res_100.h5", rd / "sern.h5")
bad = R2._restart_mismatch(rd / "res_100.h5", rd / "sern.h5", SST7)
check("(d2) dtype 違い (float64 → float32、値は一致): 拒否", bad == ["roe の dtype が違う (float64 / float32)"], bad)
rd = run_dir_cpg("d_shape", SST7)
with h5py.File(rd / "sern.h5", "r+") as f:
    del f["VALUE/roK"]; f.create_dataset("VALUE/roK", data=np.zeros(49))
bad = R2._restart_mismatch(rd / "res_100.h5", rd / "sern.h5", SST7)
check("(d2) shape 違い: 拒否", any(b.startswith("roK の shape が違う") for b in bad), bad)
bad = R2._restart_mismatch(rd / "res_100.h5", rd / "sern.h5", SST7 + ["roY0"])
check("(d2) 必須の化学種が継続元に無い: 拒否", "roY0 が継続元に無い" in bad, bad)

# (d3) run_judgment の経路: config の必須量 (SST) で検査し、roK・roOmega を写さない継続を止める
PYAML = tmp / "p_j.yaml"
PYAML.write_text(yaml.safe_dump({**yaml.safe_load(open(PROD10)),
                                 "evaluate": {"floor_events": 1, "floor_events_window_steps": 40, "judgment_steps": 80, "judgment_extend_once": 1,
                                              "outStepInterval": 2, "model": "sst"}}, allow_unicode=True))
CFG_SST = 'turbulence: {model: "sst"}\ntime:\n  last: {nStepOuter: 20}\n  outStepInterval: 2\noutput: {floorEvents: 1}\n'


def make_rd(name, cfg=CFG_SST):
    rd = tmp / name; rd.mkdir()
    (rd / "solverConfig_main.yaml").write_text(cfg); (rd / "solverConfig.yaml").write_text(cfg)
    write_h5(rd / "sern.h5", SST7, dtype=np.float32); write_h5(rd / "res_20.h5", SST7, dtype=np.float32)
    (rd / "residual_history.csv").write_text("step,phase\n20,outer_end\n")
    return rd


def copy_keys(keys):
    def _c(res, mesh):
        with h5py.File(res, "r") as s, h5py.File(mesh, "r+") as d_:
            for k in keys:
                d_["VALUE"][k][:] = s["VALUE"][k][:]
    return _c


launched = []


def fake_forge(rd, monitor, poll_s):
    launched.append(rd)
    write_h5(rd / "res_80.h5", SST7, dtype=np.float32)
    return 0, {}


def fake_collect(problem_path, rd, rc=None, require_residual_pass=False):
    g = {"verdict": "PASS", "fail_class": None, "reasons": [], "rc": 0, "steadiness": {"series": {}}}
    (Path(rd) / "metrics.json").write_text("{}")
    return {"gates": g, "step": 80}


_saved_fns = (R2._run_forge_monitored, R2.restart_by_index, R2._ENV)
try:
    R2._run_forge_monitored = fake_forge
    R2.restart_by_index = copy_keys(SST7[:5])           # roK・roOmega を写さない継続
    rd = make_rd("j_noturb")
    check("(d3) run_judgment: roK・roOmega を写さない継続は時間積分の前に止まる (config の SST から必須)",
          raises(lambda: R2.run_judgment(PYAML, rd, fake_collect, 0), RuntimeError) and not launched)
    R2.restart_by_index = copy_keys(SST7)
    # (e1) JUDGMENT.json に起動の環境 (forge に渡す _ENV の実際の値)
    R2._ENV = {"PATH": "/usr/bin", "FORGE_CUDA_BLOCKSIZE": "128", "FORGE_KERNEL_SYNC": "1"}
    rd = make_rd("j_ok")
    R2.run_judgment(PYAML, rd, fake_collect, 0)
    j = json.loads((rd / R2.JUDGMENT_FILE).read_text())
    check("(d3) 全部写す継続は通り、JUDGMENT.json に必須の保存量の集合を残す",
          j["verdict"] == "PASS" and j["required_state"] == SST7 and len(launched) == 1, j.get("reason"))
    e = j["env"]
    check("(e1) JUDGMENT.json の env: FORGE_CUDA_BLOCKSIZE は値 '128'、FORGE_ALLOW_UNVERIFIED_SPECIES は null で unset に名前、"
          "その他の FORGE_* (FORGE_KERNEL_SYNC) も残す、各段 (attempts[].env) にも同じ値",
          e["named"]["FORGE_CUDA_BLOCKSIZE"] == "128" and e["named"]["FORGE_ALLOW_UNVERIFIED_SPECIES"] is None
          and "FORGE_ALLOW_UNVERIFIED_SPECIES" in e["unset"] and "FORGE_CUDA_BLOCKSIZE" not in e["unset"]
          and e["other_forge_env"] == {"FORGE_KERNEL_SYNC": "1"} and j["attempts"][0]["env"] == e, e)
    check("(e1) JUDGMENT.json の plan: 必須量・出力間隔・予定の終了 step",
          j["plan"] == {"required": list(SST4), "out_interval": 2, "planned_end": 80, "n_planned": 40}, j.get("plan"))
finally:
    R2._run_forge_monitored, R2.restart_by_index, R2._ENV = _saved_fns

# ======================================================================== R2e (実際の run_forge・監視つきの起動の記録)
fake_tools = tmp / "fake_tools"; fake_tools.mkdir()
(fake_tools / "run_case.sh").write_text("#!/bin/bash\necho fake run_case\nexit 0\n"); os.chmod(fake_tools / "run_case.sh", 0o755)
_saved_env = (R2.FORGE_TOOLS, R2._ENV)
try:
    R2.FORGE_TOOLS = fake_tools
    R2._ENV = {"PATH": os.environ.get("PATH", ""), "FORGE_CUDA_BLOCKSIZE": "128", "FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}
    rd = tmp / "e_launch"; rd.mkdir()
    (rd / "solverConfig.yaml").write_text("time: {}\n")
    rc1 = R2.run_forge(rd)
    rc2, _st = R2._run_forge_monitored(rd, monitor=lambda *a, **k: True, poll_s=0.01)
    recs = [json.loads(x) for x in (rd / R2.RUNNER_LAUNCHES).read_text().splitlines() if x.strip()]
    check("(e1) run_forge・監視つきの起動とも runner_launches.jsonl に 1 行ずつ: config の FNV-1a (forge_launches.jsonl の cfg_fnv と同じ)・"
          "環境の実際の値 (その実行で許可した FORGE_ALLOW_UNVERIFIED_SPECIES=1 も値として残る)",
          rc1 == 0 and rc2 == 0 and len(recs) == 2 and all(r_["cfg_fnv"] == SM.fnv1a64((rd / "solverConfig.yaml").read_bytes()) for r_ in recs)
          and all(r_["env"]["named"]["FORGE_ALLOW_UNVERIFIED_SPECIES"] == "1" and r_["env"]["named"]["FORGE_CUDA_BLOCKSIZE"] == "128"
                  and r_["env"]["unset"] == ["FORGE_CUDA_BLOCKSIZE_SMALL", "FORGE_BIN"] for r_ in recs), recs)
finally:
    R2.FORGE_TOOLS, R2._ENV = _saved_env
env_noallow = {k: v for k, v in os.environ.items() if k != "FORGE_ALLOW_UNVERIFIED_SPECIES"}
check("(e1) runner は FORGE_ALLOW_UNVERIFIED_SPECIES を既定で付けない (呼び出し側の環境に無ければ記録も null)",
      R2.launch_env(env_noallow)["named"]["FORGE_ALLOW_UNVERIFIED_SPECIES"] is None
      and ("FORGE_ALLOW_UNVERIFIED_SPECIES" in R2._ENV) == ("FORGE_ALLOW_UNVERIFIED_SPECIES" in os.environ))

shutil.rmtree(tmp, ignore_errors=True)
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
