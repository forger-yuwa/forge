#!/usr/bin/env python3
"""残差の化学種・トレーサの列 (plan tooling-sern-te-wake-grid §5.1 R2i) と残差の時系列の完全性 (§5.1 R2j) の受理試験
(codex result 再レビュー 2026-10-12 `notes/reviews/2026-10-12-tooling-sern-te-wake-grid-result.md` の M1・M2)。forge は回さない (0 step)。

合格 (測る前に固定):
  (i)  R2c (run_1100) の 1 回目の session の残差 (生の CSV。無ければ 40 行おきの fixture を試験の中で step-hold で 0..19999 に展開) で、
       残差以外の証拠を保存値に固定し residual_health → evaluate_gates → judgment_extendable を通す:
         1 回目相当 (他の証拠 = 1 回目の保存値。窓が未定常): 完全な残差 → FAIL/JUDGMENT_UNSETTLED・延長可 (保存された判定と同じ)。
           rms_roY0 だけ・rms_roY0 と rms_roY1 の両方を消す → FAIL/JUDGMENT_UNVERIFIABLE・延長不可
         2 回目相当 (他の証拠 = 2 回目の保存値。PASS): 完全な残差 → PASS。消した 2 種 → FAIL/JUDGMENT_UNVERIFIABLE
       変更前 (PRE) の実装では消した 2 種とも 1 回目相当は延長可・2 回目相当は PASS (codex の再現)。
       共有の check_convergence (既定 = CLI) は 1 本だけ消えた CSV を判定不能にし、全部消えた CSV は旧版の forge と区別できないので
       従来どおり (require_species=True の SERN だけが塞ぐ)
  (j)  同じ残差の時系列だけを崩す: 先頭 step 0 の 1 行だけ・末尾 step 19999 の 1 行だけ・途中の 1 step の欠測・1 step の重複・末尾の
       切り詰め・最終場の step の不一致 (res_19999・res_40000)・最終場が無い・別の session の連結 (step が戻る)・範囲外の step の行・
       未知の phase・config に nStepOuter が無い → 1 回目相当・2 回目相当とも FAIL/JUDGMENT_UNVERIFIABLE・延長不可。
       変更前は 1 行だけの 2 種が 1 回目相当で延長可・2 回目相当で PASS (codex の再現)。完全な残差は (i) の完全と同じ
  (s)  保存済みの run: case/46.sern_design の metrics.json・res_*.h5・prepare_info.json のある run (R2h の (s) --all と同じ選び方) で、
       変更前後の residual_health と evaluate_gates (verdict・fail_class・reasons) を比べる (力の履歴は保存されたものを両方に渡す)。
       時系列が完全な run (この試験が CSV・config・res の名前から独立に判定) は全部一致。不完全な run は新しい実装で残差の入力が
       判定不能になり、PASS/FAIL は変わらない (fail_class が変わる run は一覧に出す)。R2c の生の残差 (2 session) は完全で、
       新しい実装でも保存された判定 (1 回目 JUDGMENT_UNSETTLED・延長可、2 回目 PASS) と一致し、analyze の各列は保存された VERDICT と一致

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_residual_series_tests.py [--r2c-raw DIR | --fixture]
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(TOOLS))

import check_convergence as CC  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

DEFAULT_RAW = "/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad/r2c_raw"
RAW_SHA256 = {"residual_judgment1.csv": "1679d4c919d83e211e095a833d3bd6664c76be1c78ffe47fbf08786f1f2c7155",    # 1 回目 (12,964,240 B)
              "residual_history.csv": "36875370a8de5f01416d2197a65e2dd7f6937358c25c1acf3dd746980ed1ca4c"}      # 2 回目 (12,993,076 B)
ap = argparse.ArgumentParser()
ap.add_argument("--r2c-raw", default=os.environ.get("SERN_R2C_RAW", DEFAULT_RAW),
                help="R2c の生の残差 (residual_judgment1.csv = 1 回目、residual_history.csv = 2 回目) を置いたディレクトリ")
ap.add_argument("--fixture", action="store_true", help="生の CSV なしに、間引いた fixture (展開して使う) だけで回す")
args = ap.parse_args()

fails = 0
n_checks = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_residual_series_"))
C46 = ROOT / "case" / "46.sern_design"
R2C = ROOT / "notes" / "investigations" / "2026-10-12-r2c" / "run_1100_r2c_lsw10_m6_judgment"
FIXTURE = R2C / "residual_judgment1_stride40.csv.gz"
PROD10 = C46 / "problem_3d_prod_3op_wallres_lswx10.yaml"
PRE = "9e5f0ffd"          # R2i・R2j の変更前 (R2h の後)
N_J = 20000               # 判定区間の段の nStepOuter (solverConfig_judgment.yaml)


def check(name, ok, detail=""):
    global fails, n_checks
    n_checks += 1
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def stop(msg):
    print("止める: " + msg)
    sys.exit(2)


def git_show(path):
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{PRE}:{path}"], capture_output=True, text=True)
    if r.returncode != 0:
        stop(f"変更前の {path} ({PRE}) を git から読めない: {r.stderr.strip()[:200]}")
    return r.stdout


def old_module(path, name, package):
    m = types.ModuleType(name)
    m.__package__ = package
    m.__file__ = str(ROOT / path)
    exec(compile(git_show(path), f"{PRE}:{path}", "exec"), m.__dict__)
    return m


# --- 変更前 (commit PRE) のモジュール。sern_gates は変更前の check_convergence を使わせる (変更後の analyze を通すと再現にならない)
OLD_CC = old_module("solver_density_cuda/tools/check_convergence.py", "_check_convergence_pre_r2i", None)
OLD = old_module("design/forge_design/metrics/sern_gates.py", "forge_design.metrics._sern_gates_pre_r2i", "forge_design.metrics")
OLD._cc = OLD_CC
if hasattr(OLD, "residual_series_check") or hasattr(OLD_CC, "species_required"):
    stop(f"{PRE} のモジュールが変更前の形でない")


def load_hist(path):
    rows = list(csv.DictReader(open(path)))
    return [{k: (int(v) if k == "step" else float(v)) for k, v in r.items() if v not in ("", None)} for r in rows]


# --- 残差以外の証拠: 1 回目 (未定常で延長した回) と 2 回目 (PASS) の保存値 ---------------------------------------------
P10 = load_problem(PROD10)
JP = R2.judgment_plan(P10.evaluate)
CTX = {}
for k, met, fh, mon in ((1, "metrics_judgment1.json", "force_history_judgment1.csv", "TE_MONITOR_judgment1.csv"),
                        (2, "metrics.json", "force_history.csv", "TE_MONITOR.csv")):
    saved = json.loads((R2C / met).read_text())["gates"]
    hist = load_hist(R2C / fh)
    CTX[k] = {"saved": saved, "hist": hist, "mon": G.read_monitor_csv(R2C / mon)}
check("前提: 保存された判定は 1 回目 FAIL/JUDGMENT_UNSETTLED、2 回目 PASS (§6 #11 の合格の形)",
      (CTX[1]["saved"]["verdict"], CTX[1]["saved"]["fail_class"], CTX[2]["saved"]["verdict"]) == ("FAIL", "JUDGMENT_UNSETTLED", "PASS"))


def gates(mod, d, k):
    """残差 (residual_health・residual_scale_gate) だけを実際に評価し、場・最終場の床・床事象は k 回目の保存値に固定、力の履歴と
    判定区間の §6 #2・#3 は k 回目の保存系列から計算して evaluate_gates → judgment_extendable を通す。戻り = (gates, (延長可否, 理由))。"""
    c = CTX[k]
    jg = mod.judgment_gate(c["hist"], JP["window_steps"], JP["required"], JP["planned_steps"], monitor_rows=c["mon"])
    vals = {"field_health": c["saved"]["field"], "floor_gate": c["saved"]["floors"], "floor_event_gate": c["saved"]["floor_events"]}
    orig = {n: getattr(mod, n) for n in vals}
    for n, v in vals.items():
        setattr(mod, n, lambda *a, _v=v, **kw: _v)
    try:
        g = mod.evaluate_gates(d, c["hist"], 0, p_min=float(P10.evaluate.get("p_min", 1.0)), require_floor_events=True,
                               floor_window_steps=JP["window_steps"], judgment=jg, model="sst")
    finally:
        for n, v in orig.items():
            setattr(mod, n, v)
    return g, mod.judgment_extendable(g)


def mkrun(name, hdr, body, drop=(), res=(N_J,), cfg_edit=None):
    """判定区間の段の run ディレクトリ: solverConfig.yaml (= solverConfig_judgment.yaml、cfg_edit で書き換え)、残差 (drop の列を除く)、
    最終場 res_<n>.h5 (中身は使わない: 場のゲートは保存値に固定し、時系列の照合は名前だけを見る)。"""
    d = tmp / name
    d.mkdir()
    cfg = (R2C / "solverConfig_judgment.yaml").read_text()
    (d / "solverConfig.yaml").write_text(cfg_edit(cfg) if cfg_edit else cfg)
    keep = [i for i, c in enumerate(hdr) if c not in drop]
    with open(d / "residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([hdr[i] for i in keep])
        for r in body:
            w.writerow([r[i] for i in keep])
    for n in res:
        (d / f"res_{n}.h5").write_bytes(b"")
    return d


def expand_stride(body, n_steps=N_J):
    """40 行おきに間引いた fixture を、各行の値を次の行の step の手前まで保持して outer_end の step 0..N−1 に戻す (試験だけの仕掛け。
    受理ゲートの時系列の照合には手を入れない)。列ごとの判定が生の CSV と同じことは、生の CSV があるとき下で確かめる"""
    out = []
    for k, r in enumerate(body):
        s1 = int(body[k + 1][0]) if k + 1 < len(body) else n_steps
        out += [[str(s)] + r[1:] for s in range(int(r[0]), s1)]
    return out


def col_classes(report):
    """analyze の列ごとの (trend, 合否) と全体の旗 (ok, nan, stalled, converging)。"""
    cols = {c: ((re.search(r"\b(rising|falling|flat)\b", m) or [None, None])[1] if "drop=" in m else m.split("<--")[0].strip(), ok)
            for c, (m, ok) in report[1].items()}
    return cols, tuple(bool(x) for x in report[2:])


# --- 残差の入力: 生の CSV (無ければ展開した fixture) -----------------------------------------------------------------
raw_dir = Path(args.r2c_raw)
raw = {}
if not args.fixture:
    for fn, sha in RAW_SHA256.items():
        p = raw_dir / fn
        if not p.exists():
            stop(f"R2c の生の残差 {p} が無い (12 MB の CSV はリポジトリに置かない)。--r2c-raw で場所を渡すか、--fixture で fixture だけで回す")
        if hashlib.sha256(p.read_bytes()).hexdigest() != sha:
            stop(f"{p} の sha256 が登録した R2c の残差と違う")
        rows = list(csv.reader(open(p)))
        raw[fn] = (rows[0], rows[1:])
with gzip.open(FIXTURE, "rt", newline="") as fh:
    fx_rows = list(csv.reader(fh))
FX = (fx_rows[0], expand_stride(fx_rows[1:]))
HDR, BODY = raw["residual_judgment1.csv"] if raw else FX
SRC = "生の CSV (1 回目、60000 行)" if raw else f"fixture ({FIXTURE.name} を step-hold で {len(FX[1])} 行に展開)"
print(f"    残差の入力: {SRC}")
if raw:
    da, db = mkrun("fx_vs_raw_fx", *FX), mkrun("fx_vs_raw_raw", HDR, BODY)
    ca, cb = col_classes(CC.analyze(str(da / "residual_history.csv"), 3.0, 0.2)), col_classes(CC.analyze(str(db / "residual_history.csv"), 3.0, 0.2))
    check("展開した fixture の列ごとの trend・合否と全体の旗は生の CSV と同じ (fixture だけで回すときの代わりになる)", ca == cb, (ca, cb))
    ga, gb = gates(G, da, 1), gates(G, db, 1)
    check("展開した fixture の 1 回目相当の判定 (verdict・fail_class・延長可否) は生の CSV と同じ",
          (ga[0]["verdict"], ga[0]["fail_class"], ga[1][0]) == (gb[0]["verdict"], gb[0]["fail_class"], gb[1][0]))

# ======================================================================== (i) 化学種の残差列 (R2i)
SC = {}       # シナリオ → {(実装, 回): (verdict, fail_class, 延長可)}


def scenario(label, d, pre=True):
    out = {}
    for name, mod in (("new", G),) + ((("pre", OLD),) if pre else ()):
        for k in (1, 2):
            g, e = gates(mod, d, k)
            out[(name, k)] = (g["verdict"], g["fail_class"], e[0], g)
    SC[label] = out
    return out


full = scenario("完全", mkrun("i_full", HDR, BODY))
g1, g2 = full[("new", 1)][3], full[("new", 2)][3]
check("(i) 完全な残差: 1 回目相当は FAIL/JUDGMENT_UNSETTLED・延長可 (保存された 1 回目と同じ)、2 回目相当は PASS。入力は妥当・プラトー受理は不変",
      full[("new", 1)][:3] == ("FAIL", "JUDGMENT_UNSETTLED", True) and full[("new", 2)][:3] == ("PASS", None, False)
      and g1["residual"]["input"]["ok"] and g1["residual"]["verdict"] == "NOT CONVERGED (stalled/plateau)"
      and g1["residual"]["series"]["n_outer_end"] == N_J and g1["residual"]["series"]["final_field"] == f"res_{N_J}.h5",
      {k: v[:3] for k, v in full.items()})
check("(i) 完全な残差: 変更前と同じ判定 (verdict・fail_class・理由・延長可否)",
      all((full[("new", k)][3]["verdict"], full[("new", k)][3]["fail_class"], full[("new", k)][3]["reasons"], full[("new", k)][2])
          == (full[("pre", k)][3]["verdict"], full[("pre", k)][3]["fail_class"], full[("pre", k)][3]["reasons"], full[("pre", k)][2]) for k in (1, 2)))
for lab, drop in (("rms_roY0 だけ消す", ("rms_roY0",)), ("rms_roY0・rms_roY1 を消す", ("rms_roY0", "rms_roY1"))):
    d = mkrun("i_" + "_".join(drop), HDR, BODY, drop=drop)
    r = scenario(lab, d)
    rh = r[("new", 1)][3]["residual"]
    check(f"(i) {lab}: 1 回目相当・2 回目相当とも FAIL/JUDGMENT_UNVERIFIABLE・延長不可 (残差の入力に {'・'.join(drop)})",
          r[("new", 1)][:3] == ("FAIL", "JUDGMENT_UNVERIFIABLE", False) and r[("new", 2)][:3] == ("FAIL", "JUDGMENT_UNVERIFIABLE", False)
          and rh["unverifiable"] and all(c in " ".join(rh["input"]["reasons"]) for c in drop)
          and G.judgment_extendable(r[("new", 1)][3])[1] == "未定常以外の不合格: 残差の入力 (判定不能)", {k: v[:3] for k, v in r.items()})
    check(f"(i) {lab}: 変更前の実装では 1 回目相当は延長可・2 回目相当は GATES PASS だった (codex の再現)",
          r[("pre", 1)][:3] == ("FAIL", "JUDGMENT_UNSETTLED", True) and r[("pre", 2)][:3] == ("PASS", None, False), {k: v[:3] for k, v in r.items()})
    rep = CC.analyze(str(d / "residual_history.csv"), 3.0, 0.2)          # 共有ツールの既定 (CLI)
    if len(drop) == 1:
        check(f"(i) {lab}: 共有の check_convergence (既定) も判定不能 (rms_roY1 が残る = 列を書く版の forge なので設定の 2 列を要求)",
              rep[2] is False and "rms_roY0" in rep[1].get("(入力)", ("",))[0])
    else:
        check(f"(i) {lab}: 共有の check_convergence (既定) は化学種の列が 1 本も無い CSV を旧版の forge と区別できず従来どおり "
              f"(require_species=True では判定不能 = SERN の受理が塞ぐ)",
              "(入力)" not in rep[1] and CC.analyze(str(d / "residual_history.csv"), 3.0, 0.2, require_species=True)[2] is False)

# ======================================================================== (j) 残差の時系列 (R2j)
OE = [r for r in BODY if r[2] == "outer_end"]


def steps_of(rows, lo, hi):
    return [r for r in rows if lo <= int(r[0]) < hi]


blk = steps_of(BODY, 10000, 10001)
VARIANTS = [
    ("先頭 step 0 の outer_end 1 行だけ", dict(body=[OE[0]]), True),
    ("末尾 step 19999 の outer_end 1 行だけ", dict(body=[OE[-1]]), True),
    ("途中の 1 step (10000) の欠測", dict(body=steps_of(BODY, 0, 10000) + steps_of(BODY, 10001, N_J)), False),
    ("1 step (10000) の重複", dict(body=steps_of(BODY, 0, 10001) + blk + steps_of(BODY, 10001, N_J)), False),
    ("末尾 3 step の切り詰め", dict(body=steps_of(BODY, 0, N_J - 3)), False),
    ("最終場が res_19999 (残差の session と合わない)", dict(body=BODY, res=(N_J - 1,)), False),
    ("最終場が res_40000 (後の session の場)", dict(body=BODY, res=(N_J, 2 * N_J)), False),
    ("最終場が無い", dict(body=BODY, res=()), False),
    ("別の session の連結 (前の session の末尾 5000 step の後に本 session)", dict(body=steps_of(BODY, 15000, N_J) + BODY), False),
    ("session の範囲外の step (20000) の行", dict(body=BODY + [[str(N_J)] + OE[-1][1:2] + ["outer_begin"] + OE[-1][3:]]), False),
    ("未知の phase の行", dict(body=BODY[:100] + [BODY[100][:2] + ["outer_mid"] + BODY[100][3:]] + BODY[100:]), False),
    ("config に time.last.nStepOuter が無い", dict(body=BODY, cfg_edit=lambda t: re.sub(r"last: \{nStepOuter: \d+\}", "last: {}", t)), False),
]
for lab, kw, codex in VARIANTS:
    body = kw.pop("body")
    d = mkrun("j_" + str(len(SC)), HDR, body, **kw)
    r = scenario(lab, d, pre=True)
    rh = r[("new", 1)][3]["residual"]
    check(f"(j) {lab}: 1 回目相当・2 回目相当とも FAIL/JUDGMENT_UNVERIFIABLE・延長不可 (残差の入力が判定不能)",
          r[("new", 1)][:3] == ("FAIL", "JUDGMENT_UNVERIFIABLE", False) and r[("new", 2)][:3] == ("FAIL", "JUDGMENT_UNVERIFIABLE", False)
          and rh["unverifiable"] and not rh["input"]["ok"] and rh["input"]["reasons"]
          and G.judgment_extendable(r[("new", 1)][3])[1] == "未定常以外の不合格: 残差の入力 (判定不能)",
          ({k: v[:3] for k, v in r.items()}, rh["input"]["reasons"][:3]))
    if codex:
        check(f"(j) {lab}: 変更前の実装では NOT CONVERGED (stalled/plateau) として 1 回目相当は延長可・2 回目相当は GATES PASS だった (codex の再現)",
              r[("pre", 1)][:3] == ("FAIL", "JUDGMENT_UNSETTLED", True) and r[("pre", 2)][:3] == ("PASS", None, False)
              and r[("pre", 1)][3]["residual"]["verdict"] == "NOT CONVERGED (stalled/plateau)", {k: v[:3] for k, v in r.items()})
check("(j) 完全な残差 (raw は 3 phase × 20000 step、fixture は outer_end だけ) は時系列の照合を通る (上の (i) の完全と同じ判定)",
      g1["residual"]["series"]["first"] == 0 and g1["residual"]["series"]["last"] == N_J - 1 and not G.residual_series_check(tmp / "i_full")[0])

print("\n    codex の再現シナリオの判定 (新 / 変更前、1 回目相当 = 延長、2 回目相当 = 受理):")
for lab in ("完全", "rms_roY0 だけ消す", "rms_roY0・rms_roY1 を消す", "先頭 step 0 の outer_end 1 行だけ", "末尾 step 19999 の outer_end 1 行だけ"):
    r = SC[lab]
    print(f"      {lab:34s} 新: 1 回目 {r[('new', 1)][1]}・延長 {'可' if r[('new', 1)][2] else '不可'} / 2 回目 {r[('new', 2)][0]} {r[('new', 2)][1] or ''}"
          f"   変更前: 1 回目 {r[('pre', 1)][1]}・延長 {'可' if r[('pre', 1)][2] else '不可'} / 2 回目 {r[('pre', 2)][0]} {r[('pre', 2)][1] or ''}")

# ======================================================================== (s) 保存済みの run
if raw:
    for fn, k, verd in (("residual_judgment1.csv", 1, "CONVERGENCE_VERDICT_judgment1.txt"), ("residual_history.csv", 2, "CONVERGENCE_VERDICT.txt")):
        d = mkrun(f"s_r2c_{k}", *raw[fn])
        (g, e), (go, eo) = gates(G, d, k), gates(OLD, d, k)
        rep = CC.analyze(str(d / "residual_history.csv"), 3.0, 0.2)
        txt = (R2C / verd).read_text()
        same_txt = all(m.split("<--")[0].strip() in txt for c, (m, _) in rep[1].items())
        sv = CTX[k]["saved"]
        check(f"(s) R2c の生の残差 ({k} 回目、{fn}): 時系列は完全、判定は保存値 ({sv['verdict']}/{sv['fail_class']}) と変更前の実装に一致し、"
              f"analyze の各列は保存された VERDICT ({verd}) の行と一致",
              not G.residual_series_check(d)[0] and (g["verdict"], g["fail_class"]) == (sv["verdict"], sv["fail_class"])
              and (g["verdict"], g["fail_class"], g["reasons"], e) == (go["verdict"], go["fail_class"], go["reasons"], eo) and same_txt,
              (g["verdict"], g["fail_class"], g["residual"]["input"]))
else:
    print("    (生の CSV が無いので R2c の 2 回目の生の残差の照合は省いた)")


def _main_case_dirs():
    yield ROOT / "case"
    r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--path-format=absolute", "--git-common-dir"], capture_output=True, text=True)
    if r.returncode == 0 and r.stdout.strip():
        yield Path(r.stdout.strip()).parent / "case"


def series_complete(d):
    """この試験の独立な判定: outer_end の step 列 = 0..N−1 (N = config の time.last.nStepOuter)、全行の step が 0..N−1 で昇順、
    最後の res_<n>.h5 の n = N。戻り = (完全か, 説明)。"""
    try:
        N = yaml.safe_load((d / "solverConfig.yaml").read_text())["time"]["last"]["nStepOuter"]
    except Exception as e:  # noqa: BLE001
        return False, f"nStepOuter を読めない ({type(e).__name__})"
    oe, allsteps = [], []
    with open(d / "residual_history.csv", newline="") as fh:
        rd = csv.reader(fh)
        h = next(rd)
        si, pi = h.index("step"), h.index("phase")
        for row in rd:
            if row:
                allsteps.append(int(row[si]))
                if row[pi] == "outer_end":
                    oe.append(int(row[si]))
    res = sorted(int(m.group(1)) for f in d.glob("res_*.h5") for m in [re.fullmatch(r"res_(\d+)\.h5", f.name)] if m)
    ok = oe == list(range(N)) and allsteps == sorted(allsteps) and all(0 <= s < N for s in allsteps) and bool(res) and res[-1] == N
    return ok, f"N {N}・outer_end {len(oe)} 行 (最後 {oe[-1] if oe else None})・最後の場 {res[-1] if res else None}"


RESID_KEYS = ("ok", "converged", "verdict", "nan", "rising", "stalled", "converging", "columns", "last_step", "reasons", "unverifiable", "input")
names, seen_n = [], set()
for b in _main_case_dirs():
    for m in sorted((b / "46.sern_design").glob("run_*/metrics.json")):
        d = m.parent
        if d.name not in seen_n and list(d.glob("res_[0-9]*.h5")) and (d / "prepare_info.json").exists():
            seen_n.add(d.name); names.append(d)
n_same = n_inc = 0
bad_same, bad_inc, fc_change, verdict_change, inc_list = [], [], [], [], []
for d in names:
    info = json.loads((d / "prepare_info.json").read_text())
    model = info.get("model")
    if model not in ("sst", "euler"):
        print(f"    (飛ばす {d.name}: prepare_info.json の model {model!r})")
        continue
    met = json.loads((d / "metrics.json").read_text())
    saved = met.get("gates") or {}
    hist = met.get("history") or (load_hist(d / "force_history.csv") if (d / "force_history.csv").exists() else [])
    pk = [k for k in ((saved.get("floors") or {}).get("counts") or {}) if k.startswith("P<=")]
    kw = dict(require_residual_pass=bool(saved.get("require_residual_pass")), require_residual_plateau=bool(saved.get("require_residual_plateau")),
              p_min=float(pk[0][3:]) if pk else None)
    rc = met.get("forge_rc", saved.get("rc"))
    rn, ro = G.residual_health(d), OLD.residual_health(d)
    en, eo = G.evaluate_gates(d, hist, rc, model=model, **kw), OLD.evaluate_gates(d, hist, rc, model=model, **kw)
    complete, why = series_complete(d)
    same_r = json.dumps({k: rn.get(k) for k in RESID_KEYS}, sort_keys=True, default=str) == json.dumps({k: ro.get(k) for k in RESID_KEYS}, sort_keys=True, default=str)
    same_g = (en["verdict"], en["fail_class"], en["reasons"]) == (eo["verdict"], eo["fail_class"], eo["reasons"])
    if complete:
        n_same += 1
        if not (same_r and same_g and not G.residual_series_check(d)[0]):
            bad_same.append((d.name, why, en["verdict"], en["fail_class"], eo["verdict"], eo["fail_class"]))
    else:
        n_inc += 1
        inc_list.append((d.name, why, f"{eo['verdict']}/{eo['fail_class']} → {en['verdict']}/{en['fail_class']}", rc))
        if not (rn["unverifiable"] and any("判定不能" in r for r in rn["input"]["reasons"])):
            bad_inc.append((d.name, why, rn["input"]))
        if en["verdict"] != eo["verdict"]:
            verdict_change.append((d.name, eo["verdict"], en["verdict"]))
        if en["fail_class"] != eo["fail_class"]:
            fc_change.append((d.name, why, f"{eo['verdict']}/{eo['fail_class']} → {en['verdict']}/{en['fail_class']}"))
check(f"(s) case/46 の保存 run {n_same + n_inc} 本: 時系列が完全な {n_same} 本は変更前後で residual_health・evaluate_gates が一致",
      n_same + n_inc >= 100 and not bad_same, bad_same[:5])
check(f"(s) 時系列が不完全な {n_inc} 本 (発散で切り詰め・最終場の不一致): 新しい実装で残差の入力が判定不能", not bad_inc, bad_inc[:5])
check(f"(s) 不完全な {n_inc} 本でも PASS/FAIL は変わらない (変わった run: {verdict_change})", not verdict_change, verdict_change)
print(f"    (s) 時系列が不完全な保存 run ({len(inc_list)} 本。判定 = 変更前 → 変更後、rc = 保存された forge の終了コード):")
for x in inc_list:
    print(f"        {x[0]}: {x[2]}  rc {x[3]}  ({x[1]})")
print(f"    (s) fail_class が変わった保存 run ({len(fc_change)} 本、PASS/FAIL は不変。要承認として報告する):")
for x in fc_change:
    print(f"        {x[0]}: {x[2]}  ({x[1]})")

print(f"\n{n_checks - fails}/{n_checks} PASS" + ("" if not fails else f"  ({fails} FAIL)"))
print(f"(作業ディレクトリ {tmp})")
shutil.rmtree(tmp, ignore_errors=True)
sys.exit(1 if fails else 0)
