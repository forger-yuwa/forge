#!/usr/bin/env python3
"""生産の標準実行経路の判定区間の段の試験 (plan tooling-sern-te-wake-grid §5.1 R2、codex result レビュー 2026-10-11 M2)。
forge は回さない (GPU の接続試験は run_sern_judgment_connection_test.py)。

合格 (測る前に固定):
  (1) 設定: evaluate.judgment_steps / judgment_extend_once の検査 (床事象のゲートと窓の長さが必須、N ≥ 2 × 窓、
      窓を短くして通さない)。生産 YAML (lswx08・lswx10) は N 20000・延長 1 回・窓 10000。判定区間の config は本段の
      config と nStepOuter だけが違い、outStepInterval が窓を割り切り 20 標本以上 (足りない設定は時間積分の前に止まる)
  (2) 判定 (§6 #2・#3) は検証の手作業と同じ: 保存された検証 run 12 本 (notes/investigations/2026-10-11-*) の force_history.csv・
      TE_MONITOR.csv を通すと、a・d が v3_farfield_eval の出力と一致し、STEADY の判定が check_quasisteady --tail 0.5 の出力と
      一致する。plan の記録どおり run_1087 (A m6_on) は低温 FAIL、run_1085・run_1088・run_1097 は未定常 (延長の対象)、
      run_1084・run_1086・run_1091・run_1099 などは PASS
  (3) 判定不能の入力を合格にしない: 力の履歴が 2 窓に足りない (段階起動の本段 4000 step に窓 10000 = レビューの再現)・
      窓の標本の不足・非有限・監視の欠け・空の領域・監視が力と別の場で終わる・監視が使えない (2D・形状が合わない)
  (4) 延長の可否: 未定常 (既存の UNSTEADY と JUDGMENT_UNSETTLED) だけなら延長、低温・床事象・判定不能・発散・残差は延長しない
  (5) 段の運び (偽の forge): PASS で 1 回、未定常 → 延長 → PASS、未定常が 2 回で INDETERMINATE、低温・床の判定不能は延長せず、
      延長しない設定、起動前の判定不能 (precheck)・段階起動の失敗で段を回さない。同じ run ディレクトリ・config の差・
      記録の退避 (*_judgment1.*)・JUDGMENT.json・継続のビット一致の検査・スナップショットの監視と削除 (最後の 1 枚を残す)

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_judgment_tests.py
"""
import csv
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))

from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

fails = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_judgment_"))
C46 = ROOT / "case" / "46.sern_design"
INV = ROOT / "notes" / "investigations"


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


# ======================================================================== (1) 設定
EV = {"floor_events": 1, "floor_events_window_steps": 10000}
for prod in ("problem_3d_prod_3op_wallres_lswx08.yaml", "problem_3d_prod_3op_wallres_lswx10.yaml"):
    jo = R2.judgment_options(load_problem(C46 / prod).evaluate)
    check(f"(1) 生産 YAML {prod}: N 20000・延長 1 回・窓 10000", jo == {"steps": 20000, "extend_once": True, "window_steps": 10000}, jo)
check("(1) キー無しは段なし (従来の経路)", R2.judgment_options(EV) is None and R2.judgment_options({}) is None)
check("(1) 延長なしの既定", R2.judgment_options({**EV, "judgment_steps": 20000})["extend_once"] is False)
for label, ev in (("床事象のゲートなし", {"judgment_steps": 20000, "floor_events_window_steps": 10000}),
                  ("窓の長さなし", {"floor_events": 1, "judgment_steps": 20000}),
                  ("N < 2 × 窓 (窓を短くして通さない)", {**EV, "judgment_steps": 19999}),
                  ("延長だけ", {**EV, "judgment_extend_once": 1}), ("N が文字列", {**EV, "judgment_steps": "20000"}),
                  ("N 0", {**EV, "judgment_steps": 0}), ("延長 2", {**EV, "judgment_steps": 20000, "judgment_extend_once": 2})):
    check(f"(1) 不正な設定は止まる: {label}", raises(lambda: R2.judgment_options(ev)))
p08 = load_problem(C46 / "problem_3d_prod_3op_wallres_lswx08.yaml")
R2.select_operating_point(p08, "m10_on")        # 作動点の組成を入れる (prepare と同じ順)
cfg_main = R3._solver_config(p08, 4000, 500, 0.25, 2851.0)
cfg_j = R2.judgment_config(cfg_main, R2.judgment_options(p08.evaluate))
dl = [(a, b) for a, b in zip(cfg_main.splitlines(), cfg_j.splitlines()) if a != b]
check("(1) 判定区間の config は本段と nStepOuter だけが違う (同一設定)",
      len(cfg_main.splitlines()) == len(cfg_j.splitlines()) and len(dl) == 1 and "nStepOuter: 4000" in dl[0][0]
      and "nStepOuter: 20000" in dl[0][1] and "output: {floorEvents: 1}" in cfg_j, dl)
jo1 = R2.judgment_options(p08.evaluate)
for oi, why in ((300, "窓を割り切らない"), (1000, "窓あたり 10 標本"), (0, "出力なし")):
    check(f"(1) outStepInterval {oi} ({why}) は時間積分の前に止まる",
          raises(lambda: R2.judgment_config(re.sub(r"outStepInterval: \d+", f"outStepInterval: {oi}", cfg_main), jo1)))
check("(1) nStepOuter の置換に失敗したら止まる", raises(lambda: R2.judgment_config("time: {}\noutStepInterval: 500\n", jo1), RuntimeError))


# ======================================================================== (2) 保存された検証 run の再判定
# 必須量と予定 step 列は問題設定から (§5.1 R2a): 検証 run は SST・N 20000・窓 10000・出力 500 step (生産 YAML と同じ)
SST4 = G.judgment_required_quantities("sst")
PLAN = G.judgment_planned_steps(20000, 10000, 500)


def load_hist(path):
    rows = list(csv.DictReader(open(path)))
    return [{k: (int(v) if k == "step" else float(v)) for k, v in r.items() if v not in ("", None)} for r in rows]


def parse_v3(path):
    out = {}
    for line in open(path):
        m = re.match(r"\s+(\S+)\s+平均 (\S+)\s+a (\S+)\s+前窓差 (\S+)", line)
        if m:
            out[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return out


def parse_qs(path):
    out = {}
    for line in open(path):
        m = re.match(r"\s+(\S+)\s+:.*\s(STEADY|OSCILLATING|TRANSIENT-UNSETTLED|DRIFTING|NONFINITE)\s*$", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


saved = sorted({p.parent for p in INV.glob("2026-10-1*/**/force_history.csv")} | {p.parent for p in INV.glob("2026-10-1*/force_history.csv")})
saved = [d for d in saved if (d / "TE_MONITOR.csv").exists() and (d / "V3_EVAL.txt").exists()]
check("(2) 保存された検証 run が 12 本ある (notes/investigations/2026-10-11-*)", len(saved) == 12, [d.name for d in saved])
verdicts = {}
for d in saved:
    name = d.name if d.name.startswith("run_") else "run_1084_tewake_g3_L10_m10"
    hist = load_hist(d / "force_history.csv")
    jg = G.judgment_gate(hist, 10000, SST4, PLAN, monitor_rows=G.read_monitor_csv(d / "TE_MONITOR.csv"))
    verdicts[name[:8]] = jg
    v3, qf, qt = parse_v3(d / "V3_EVAL.txt"), parse_qs(d / "QUASISTEADY_forces.txt"), parse_qs(d / "QUASISTEADY_Tmin.txt")
    q = jg["force_window"]["quantities"]
    same_ad = all(np.isclose(q[k]["a"], v3[k][1], rtol=6e-3, atol=0) and np.isclose(q[k]["d"], v3[k][2], rtol=6e-3, atol=0)
                  and np.isclose(q[k]["mean"], v3[k][0], rtol=0, atol=1e-6) for k in v3)
    check(f"(2) {name[:8]}: a・d・平均が v3_farfield_eval の出力と一致", same_ad and set(v3) == set(q),
          {k: (q[k]["a"], q[k]["d"], v3[k]) for k in v3})
    check(f"(2) {name[:8]}: 力の STEADY が check_quasisteady --tail 0.5 の出力と一致",
          {k: q[k]["steady"] for k in q} == {k: qf[k] for k in q}, ({k: q[k]["steady"] for k in q}, qf))
    tr = jg["temperature"]["regions"]
    check(f"(2) {name[:8]}: 最低温度の STEADY が check_quasisteady --tail 0.5 の出力と一致",
          {f"{R}_Tmin": tr[R]["steady"] for R in tr} == qt, ({R: tr[R]["steady"] for R in tr}, qt))
    rows = G.read_monitor_csv(d / "TE_MONITOR.csv")
    lt150 = sum(int(r[f"{R}_n_lt150"]) for r in rows if r["step"] > max(x["step"] for x in rows) - 10000 for R in G.TE_REGIONS)
    want_ok = (all(max(v3[k][1], v3[k][2]) <= G.JUDGMENT_FORCE_TOL[k] for k in v3) and all(v == "STEADY" for v in qf.values())
               and all(v == "STEADY" for v in qt.values()) and lt150 == 0)
    check(f"(2) {name[:8]}: 合否が検証の判定 (窓 max(a, d) ≤ 許容・STEADY・T < 150 K なし) と一致 → {jg['verdict']}",
          jg["ok"] == want_ok, jg["reasons"][:3])
EXPECT = {"run_1087": ("FAIL", "COLD_SPOT"), "run_1085": ("UNSETTLED", "JUDGMENT_UNSETTLED"), "run_1088": ("UNSETTLED", "JUDGMENT_UNSETTLED"),
          "run_1097": ("UNSETTLED", "JUDGMENT_UNSETTLED"), "run_1084": ("PASS", None), "run_1086": ("PASS", None), "run_1091": ("PASS", None),
          "run_1096": ("PASS", None), "run_1098": ("PASS", None), "run_1099": ("PASS", None)}
for k, (v, fc) in EXPECT.items():
    jg = verdicts.get(k)
    check(f"(2) plan の記録どおり {k}: {v} ({fc})", jg is not None and jg["verdict"] == v and jg["fail_class"] == fc,
          None if jg is None else (jg["verdict"], jg["fail_class"], jg["reasons"][:2]))
j87 = verdicts.get("run_1087")
check("(2) run_1087 の低温は節点 517159 の 1 節点 (約 102.5 K、判定窓の 20 標本すべて。result レビュー m1 の訂正どおり)",
      j87 is not None and j87["temperature"]["regions"]["ALL"]["n_lt150_window_sum"] == 20
      and abs(j87["temperature"]["regions"]["ALL"]["Tmin_window"] - 102.53) < 0.01
      and j87["temperature"]["regions"]["ALL"]["Tmin_node"] == 517159, None if j87 is None else j87["temperature"]["regions"]["ALL"])
j97 = verdicts.get("run_1097")
check("(2) run_1097 の未達は C_T_with_shear・C_M の窓と ALL の最低温度 (延長 run_1099 の理由と同じ)",
      j97 is not None and {k for k, v in j97["force_window"]["quantities"].items() if not v["ok"]} == {"C_T_with_shear", "C_M"}
      and j97["temperature"]["verdict"] == "UNSETTLED" and j97["extendable"])

# ======================================================================== (3) 判定不能の入力
def synth_hist(N=20000, dt=500, drift=0.0, start=500):
    return [{"step": s, "C_T": 0.9 + drift * s, "C_T_with_shear": 0.88, "C_L": 0.07, "C_M": -1.8} for s in range(start, N + 1, dt)]


def synth_mon(N=20000, dt=500, cold_at=None, region_empty=None, start=500, drop=()):
    rows = []
    for s in range(start, N + 1, dt):
        if s in drop:
            continue
        r = {"step": float(s), "nonfinite_T": 0.0, "n_floor_plus1": 0.0}
        for R in G.TE_REGIONS:
            r.update({f"{R}_n": 0.0 if R == region_empty else 100.0, f"{R}_Tmin": 220.0, f"{R}_Tmin_node": 7.0, f"{R}_Tmin_x": 0.1,
                      f"{R}_Tmin_y": 0.0, f"{R}_Tmin_z": 0.0, f"{R}_n_lt150": 1.0 if (cold_at == s and R == "ALL") else 0.0,
                      f"{R}_n_lt200": 0.0, f"{R}_vol_lt200": 0.0})
        rows.append(r)
    return rows


ok_g = G.judgment_gate(synth_hist(), 10000, SST4, PLAN, monitor_rows=synth_mon())
check("(3) 合成: 定常な力・温度 → PASS", ok_g["ok"] and ok_g["verdict"] == "PASS" and ok_g["force_window"]["n_tail"] == 20
      and ok_g["force_window"]["n_prev"] == 20 and ok_g["temperature"]["n_window"] == 20, ok_g["reasons"])
for label, kw, want in (
        ("段階起動の本段 4000 step に窓 10000 (レビューの再現)", {"hist": synth_hist(4000), "mon": synth_mon(4000)}, "JUDGMENT_UNVERIFIABLE"),
        ("出力間隔 1000 で窓の標本 10", {"hist": synth_hist(dt=1000, start=1000), "mon": synth_mon(dt=1000, start=1000)}, "JUDGMENT_UNVERIFIABLE"),
        ("力に NaN", {"hist": [dict(h, C_L=float("nan")) if h["step"] == 15000 else h for h in synth_hist()], "mon": synth_mon()},
         "JUDGMENT_UNVERIFIABLE"),
        ("監視の記録が無い", {"hist": synth_hist(), "mon": None}, "JUDGMENT_UNVERIFIABLE"),
        ("末尾窓の監視が 1 枚欠ける", {"hist": synth_hist(), "mon": synth_mon(drop=(15000,))}, "JUDGMENT_UNVERIFIABLE"),
        ("監視が力より前の場で終わる", {"hist": synth_hist(), "mon": synth_mon(N=19500)}, "JUDGMENT_UNVERIFIABLE"),
        ("領域 R_SE が空", {"hist": synth_hist(), "mon": synth_mon(region_empty="R_SE")}, "JUDGMENT_UNVERIFIABLE"),
        ("監視が使えない (2D)", {"hist": synth_hist(), "mon": synth_mon(), "why": R2.MONITOR_UNAVAILABLE_2D}, "JUDGMENT_UNVERIFIABLE"),
        ("末尾窓で ALL に T < 150 K が 1 節点", {"hist": synth_hist(), "mon": synth_mon(cold_at=12000)}, "COLD_SPOT"),
        ("低温は前窓だけ (判定の窓の外) なら合格", {"hist": synth_hist(), "mon": synth_mon(cold_at=8000)}, None),
        ("C_T のずれ (d 1e-4 > 5e-5)", {"hist": synth_hist(drift=1e-8), "mon": synth_mon()}, "JUDGMENT_UNSETTLED")):
    g = G.judgment_gate(kw["hist"], 10000, SST4, PLAN, monitor_rows=kw["mon"], monitor_unavailable=kw.get("why"))
    check(f"(3) {label} → {g['fail_class']} (期待 {want})", g["fail_class"] == want and g["ok"] == (want is None), g["reasons"][:2])
nf = synth_mon(); nf[-1]["nonfinite_T"] = 3.0
check("(3) 末尾窓の非有限の温度は低温と同じく FAIL (延長しない)", G.judgment_gate(synth_hist(), 10000, SST4, PLAN, monitor_rows=nf)["fail_class"] == "COLD_SPOT")

# ======================================================================== (4) 延長の可否
def gates_(fail=None, **over):
    base = {"verdict": "PASS" if fail is None else "FAIL", "fail_class": fail, "rc": 0, "field": {"ok": True}, "floors": {"ok": True},
            "floor_events": {"ok": True}, "residual_scale": {"ok": True}, "residual": {"ok": True, "converged": False},
            "steadiness": {"ok": True, "series": {"C_T": {}}}, "judgment": {"ok": True, "extendable": False, "fail_class": None},
            "require_residual_pass": False, "require_residual_plateau": False, "plateau": {"falling": []}}
    base.update(over)
    return base


UNS = {"ok": False, "extendable": True, "fail_class": "JUDGMENT_UNSETTLED"}
for label, g, want in (("合格", gates_(), False),
                       ("判定区間の窓・定常性だけ未達", gates_("JUDGMENT_UNSETTLED", judgment=UNS), True),
                       ("既存の UNSTEADY だけ", gates_("UNSTEADY", steadiness={"ok": False, "series": {"C_T": {}}}), True),
                       ("UNSTEADY + 窓", gates_("UNSTEADY", steadiness={"ok": False, "series": {"C_T": {}}}, judgment=UNS), True),
                       ("低温", gates_("COLD_SPOT", judgment={"ok": False, "extendable": False, "fail_class": "COLD_SPOT"}), False),
                       ("判定不能", gates_("JUDGMENT_UNVERIFIABLE", judgment={"ok": False, "extendable": False, "fail_class": "JUDGMENT_UNVERIFIABLE"}), False),
                       ("床事象 + 窓", gates_("FLOOR_EVENT", floor_events={"ok": False}, judgment=UNS), False),
                       ("rc ≠ 0 + 窓", gates_("DIVERGED", rc=1, judgment=UNS), False),
                       ("残差の上昇 + 窓", gates_("RESIDUAL_RISING", residual={"ok": False}, judgment=UNS), False),
                       ("最終場の床 + 窓", gates_("FLOOR_STUCK", floors={"ok": False}, judgment=UNS), False),
                       ("力の履歴なし", gates_("NO_FORCES", steadiness={"ok": False, "series": {}}), False),
                       ("require_residual_pass で未収束 + 窓", gates_("NOT_CONVERGED", require_residual_pass=True, judgment=UNS), False)):
    got, why = G.judgment_extendable(g)
    check(f"(4) 延長の可否: {label} → {got} (期待 {want})", got == want, why)
# evaluate_gates への接続 (判定区間の未達は FAIL、渡さなければ従来どおり)
rd0 = tmp / "eg"; rd0.mkdir()
with h5py.File(rd0 / "res_20000.h5", "w") as f:
    # SST の場 (roK・roOmega も。2026-10-12 から SST では必須: plan tooling-sern-te-wake-grid §5.1 R2g)
    for k, v in (("ro", 1.0), ("roUx", 2.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", 3e5), ("P", 1e5), ("T", 300.0),
                 ("roK", 1.0), ("roOmega", 1.0e3)):
        f.create_dataset(f"VALUE/{k}", data=np.full(30, v))
# 残差は流体 5 列 + SST の 2 列、run の solverConfig.yaml (SST) も置く: 旧 fixture の 3 列 (rms_ro・rms_roUx・rms_roe) で設定なしの
# CSV は check_convergence の判定不能 (必須列の欠落) で、2026-10-12 から受理ゲートも判定不能にする (plan tooling-sern-te-wake-grid §5.1 R2h1)
# 時系列は forge と同じ (outer_end の step 0..19999 を 1 行ずつ、最終場 res_20000、config に nStepOuter 20000 と physProp): 2026-10-12 から
# 受理ゲートは残差の時系列の欠測・終了位置の不一致 (旧 fixture は 500 step おき・step 20000 まで) と、化学種の列を決める physProp の
# 無い config を判定不能にする (plan tooling-sern-te-wake-grid §5.1 R2i・R2j)
(rd0 / "solverConfig.yaml").write_text('physProp: {thermalMethod: 0}\nturbulence: {model: "sst"}\ntime:\n  last: {nStepOuter: 20000}\n')
with open(rd0 / "residual_history.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"])
    for s in range(20000):
        v = 10 ** (-1 - 1.5 * min(s / 20000 * 3, 1)); w.writerow([s, -1, "outer_end", v, v * 10, v * 10, v * 10, v * 100, v * 10, v * 1000])
h0 = synth_hist()
g_none = G.evaluate_gates(rd0, h0, 0, model="sst")
g_pass = G.evaluate_gates(rd0, h0, 0, judgment=ok_g, model="sst")
g_uns = G.evaluate_gates(rd0, h0, 0, judgment=G.judgment_gate(synth_hist(drift=1e-8), 10000, SST4, PLAN, monitor_rows=synth_mon()), model="sst")
check("(4) evaluate_gates: judgment を渡さなければキーも出さず従来どおり、PASS の判定は合格のまま",
      "judgment" not in g_none and g_none["verdict"] == "PASS" and g_pass["verdict"] == "PASS" and g_pass["judgment"]["ok"])
check("(4) evaluate_gates: 判定区間の未達は FAIL/JUDGMENT_UNSETTLED で延長の対象",
      g_uns["verdict"] == "FAIL" and g_uns["fail_class"] == "JUDGMENT_UNSETTLED" and G.judgment_extendable(g_uns)[0])


# ======================================================================== (5) 段の運び (偽の forge)
PYAML = tmp / "p_j.yaml"
PYAML.write_text(yaml.safe_dump({**yaml.safe_load(open(C46 / "problem_3d_prod_3op_wallres_lswx08.yaml")),
                                 "evaluate": {"floor_events": 1, "floor_events_window_steps": 40, "judgment_steps": 80, "judgment_extend_once": 1,
                                              "outStepInterval": 2, "model": "sst"}}, allow_unicode=True))
CFG_MAIN = "time:\n  last: {nStepOuter: 20}\n  outStepInterval: 2\noutput: {floorEvents: 1}\n"


def write_state(path, val):
    with h5py.File(path, "w") as f:
        for k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "P"):
            f.create_dataset(f"VALUE/{k}", data=np.full(12, val + len(k), np.float32))


def make_rd(name):
    rd = tmp / name; rd.mkdir()
    (rd / "solverConfig_main.yaml").write_text(CFG_MAIN); (rd / "solverConfig.yaml").write_text(CFG_MAIN)
    write_state(rd / "sern.h5", 0.0); write_state(rd / "res_20.h5", 5.0)
    (rd / "residual_history.csv").write_text("step,phase\n20,outer_end\n")
    return rd


def simple_restart(res, mesh):
    with h5py.File(res, "r") as s, h5py.File(mesh, "r+") as d:
        for k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"):
            d["VALUE"][k][:] = s["VALUE"][k][:]


class FakeForge:
    """_run_forge_monitored の代わり: 継続した場が sern.h5 に入っていることを確かめ、res_80.h5 と記録を書く。"""

    def __init__(self):
        self.calls = []

    def __call__(self, rd, monitor, poll_s):
        with h5py.File(rd / "sern.h5", "r") as f:
            ro0 = float(f["VALUE/ro"][0])
        self.calls.append({"cfg": (rd / "solverConfig.yaml").read_text(), "ro0": ro0, "res_left": sorted(x.name for x in rd.glob("res_*")),
                           "monitor_csv_at_launch": (rd / "TE_MONITOR.csv").exists()})
        write_state(rd / "res_80.h5", 10.0 * len(self.calls))
        (rd / "residual_history.csv").write_text(f"step,phase\n80,outer_end\n# call {len(self.calls)}\n")
        (rd / "TE_MONITOR.csv").write_text(f"step\n80\n# call {len(self.calls)}\n")
        return 0, {"monitored": 40, "failed": 0, "removed": 39}


def fake_collect_seq(seq):
    calls = []

    def _c(problem_path, rd, rc=None, require_residual_pass=False):
        g = seq[min(len(calls), len(seq) - 1)]
        calls.append(rc)
        out = {"gates": {**gates_(), **g}, "step": 80}
        (Path(rd) / "metrics.json").write_text(json.dumps({"call": len(calls), "fail_class": out["gates"]["fail_class"]}))
        return out
    _c.calls = calls
    return _c


PASS_G = {}
UNS_G = {"verdict": "FAIL", "fail_class": "JUDGMENT_UNSETTLED", "judgment": UNS}
COLD_G = {"verdict": "FAIL", "fail_class": "COLD_SPOT", "judgment": {"ok": False, "extendable": False, "fail_class": "COLD_SPOT"}}
FLOORX_G = {"verdict": "FAIL", "fail_class": "FLOOR_UNVERIFIABLE", "floor_events": {"ok": False}, "judgment": UNS}
saved_fns = (R2._run_forge_monitored, R2.restart_by_index)
try:
    R2.restart_by_index = simple_restart
    for label, seq, n_forge, verdict in (("1 回目で PASS", [PASS_G], 1, "PASS"), ("未定常 → 延長 → PASS", [UNS_G, PASS_G], 2, "PASS"),
                                         ("未定常が 2 回 → INDETERMINATE", [UNS_G, UNS_G], 2, "INDETERMINATE"),
                                         ("低温 → 延長しない FAIL", [COLD_G], 1, "FAIL"),
                                         ("床の判定不能 (+窓) → 延長しない INDETERMINATE", [FLOORX_G], 1, "INDETERMINATE")):
        rd = make_rd("seq_" + str(abs(hash(label)) % 10 ** 6))
        ff = FakeForge(); R2._run_forge_monitored = ff
        col = fake_collect_seq(seq)
        rc, out = R2.run_judgment(PYAML, rd, col, 0, monitor=lambda *a, **k: True)
        j = json.loads((rd / R2.JUDGMENT_FILE).read_text())
        check(f"(5) {label}: forge {len(ff.calls)} 回 (期待 {n_forge})・JUDGMENT {j['verdict']} (期待 {verdict})",
              len(ff.calls) == n_forge and j["verdict"] == verdict and len(j["attempts"]) == n_forge, j["reason"])
        c0 = ff.calls[0]["cfg"]
        check(f"(5) {label}: 判定区間の config = 本段の nStepOuter だけ 80、solverConfig_judgment.yaml に記録",
              c0 == CFG_MAIN.replace("nStepOuter: 20", "nStepOuter: 80") and (rd / R2.JUDGMENT_CONFIG).read_text() == c0)
        check(f"(5) {label}: 本段の最終場 (res_20) から継続し、起動前に res_* を消す",
              ff.calls[0]["ro0"] == 7.0 and ff.calls[0]["res_left"] == [] and (rd / "residual_main.csv").exists())
        if n_forge == 2:
            check(f"(5) {label}: 延長は 1 回目の最終場 (res_80) から継続し、1 回目の記録を *_judgment1.* に退避",
                  ff.calls[1]["ro0"] == 12.0 and (rd / "metrics_judgment1.json").exists() and (rd / "TE_MONITOR_judgment1.csv").exists()
                  and json.loads((rd / "metrics_judgment1.json").read_text())["call"] == 1
                  and json.loads((rd / "metrics.json").read_text())["call"] == 2
                  and (rd / "residual_judgment1.csv").exists() and (rd / "residual_judgment2.csv").exists()
                  and "archived" in j["attempts"][0] and "extend_reason" in j["attempts"][0])
        check(f"(5) {label}: 各段の起動時に監視の CSV は空 (前の段の記録は退避済み)", not any(c["monitor_csv_at_launch"] for c in ff.calls))
    # 段より前の監視の記録が残っている run (再実行・手作業の残り): 退避してから回す (監視済みとみなしてスナップショットを消さない)
    rd = make_rd("stale_monitor"); (rd / "TE_MONITOR.csv").write_text("step\n2\n4\n"); (rd / "TE_MONITOR.log").write_text("old\n")
    ff = FakeForge(); R2._run_forge_monitored = ff
    R2.run_judgment(PYAML, rd, fake_collect_seq([PASS_G]), 0, monitor=lambda *a, **k: True)
    j = json.loads((rd / R2.JUDGMENT_FILE).read_text())
    check("(5) 段より前の TE_MONITOR.csv・.log は *_before_judgment1.* に退避してから回し、JUDGMENT.json に記録",
          not ff.calls[0]["monitor_csv_at_launch"] and (rd / "TE_MONITOR_before_judgment1.csv").read_text() == "step\n2\n4\n"
          and sorted(j["attempts"][0].get("stale_monitor_archived", [])) == ["TE_MONITOR_before_judgment1.csv", "TE_MONITOR_before_judgment1.log"])
    # 延長しない設定
    p_noext = tmp / "p_noext.yaml"
    y = yaml.safe_load(PYAML.read_text()); y["evaluate"]["judgment_extend_once"] = 0; p_noext.write_text(yaml.safe_dump(y, allow_unicode=True))
    rd = make_rd("noext"); ff = FakeForge(); R2._run_forge_monitored = ff
    R2.run_judgment(p_noext, rd, fake_collect_seq([UNS_G]), 0)
    j = json.loads((rd / R2.JUDGMENT_FILE).read_text())
    check("(5) 延長しない設定: 未定常 1 回で INDETERMINATE (延長しない)", len(ff.calls) == 1 and j["verdict"] == "INDETERMINATE" and "延長しない設定" in j["reason"])
    # 起動前の判定不能・段階起動の失敗
    rd = make_rd("precheck"); ff = FakeForge(); R2._run_forge_monitored = ff
    col = fake_collect_seq([PASS_G])
    R2.run_judgment(PYAML, rd, col, 0, precheck=lambda _rd: R2.MONITOR_UNAVAILABLE_2D)
    j = json.loads((rd / R2.JUDGMENT_FILE).read_text())
    check("(5) 起動前に §6 #2 を判定できない (2D): 段を回さず INDETERMINATE、collect は 1 回 (metrics.json は残る)",
          len(ff.calls) == 0 and j["verdict"] == "INDETERMINATE" and len(col.calls) == 1 and (rd / "metrics.json").exists())
    rd = make_rd("startup_fail"); ff = FakeForge(); R2._run_forge_monitored = ff
    R2.run_judgment(PYAML, rd, fake_collect_seq([PASS_G]), 1)
    check("(5) 段階起動の rc ≠ 0: 段を回さず FAIL", len(ff.calls) == 0 and json.loads((rd / R2.JUDGMENT_FILE).read_text())["verdict"] == "FAIL")
    # 継続のビット一致の検査
    R2._run_forge_monitored = FakeForge()
    R2.restart_by_index = lambda res, mesh: None          # 写さない継続 (検査が拾うこと)
    rd = make_rd("bits")
    check("(5) 継続元とビット一致しない継続は止まる", raises(lambda: R2.run_judgment(PYAML, rd, fake_collect_seq([PASS_G]), 0), RuntimeError))
finally:
    R2._run_forge_monitored, R2.restart_by_index = saved_fns
check("(5) 段なしの問題 YAML に run_judgment は使えない", raises(lambda: R2.run_judgment(C46 / "problem_3d_prod_3op_wallres_lswx08_tewake_g3_L10.yaml",
                                                                         tmp, fake_collect_seq([PASS_G]), 0)))

# スナップショットの監視と削除 (最後の 1 枚を残す、実行中は最新に触らない)
rd = tmp / "drain"; rd.mkdir()
for s in (0, 2, 4, 6):
    write_state(rd / f"res_{s}.h5", s); (rd / f"res_{s}.xmf").write_text("x")
(rd / "res_ramp_4_2.h5").write_text("wall"); (rd / "res_nan_3.h5").write_text("nan")
seen = []


def mon(rdir, res, step0=False):
    seen.append((Path(res).name, step0))
    if not step0:
        with open(Path(rdir) / R2.MONITOR_CSV, "a") as fh:
            fh.write(("step,x\n" if not (Path(rdir) / R2.MONITOR_CSV).exists() or not (Path(rdir) / R2.MONITOR_CSV).read_text() else "")
                     + f"{int(Path(res).stem.split('_')[1])},1\n")
    return Path(res).name != "res_4.h5"      # 1 枚だけ監視に失敗させる


st1 = R2._drain_snapshots(rd, mon, final=False)
left1 = sorted(x.name for x in rd.glob("res_*"))
check("(5) 実行中: 最新 (res_6) に触らず、step 0 は別枠、他は監視して消す (壁面の出力・発散ダンプは触らない)",
      seen == [("res_0.h5", True), ("res_2.h5", False), ("res_4.h5", False)] and left1 == ["res_6.h5", "res_6.xmf", "res_nan_3.h5", "res_ramp_4_2.h5"]
      and st1 == {"monitored": 2, "failed": 1, "removed": 3}, (seen, left1, st1))
st2 = R2._drain_snapshots(rd, mon, final=True)
check("(5) 終了後: 最後の 1 枚 (res_6) を監視して残す・同じ step を二重に記録しない",
      seen[-1] == ("res_6.h5", False) and (rd / "res_6.h5").exists() and st2 == {"monitored": 1, "failed": 0, "removed": 0}
      and R2._drain_snapshots(rd, mon, final=True)["monitored"] == 0)

shutil.rmtree(tmp, ignore_errors=True)
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
