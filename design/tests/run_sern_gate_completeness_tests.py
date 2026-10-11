#!/usr/bin/env python3
"""入力の欠落を判定不能 (または継承の停止) にする試験 (plan tooling-sern-te-wake-grid §5.1 R2h1〜R2h5、codex diagnose
`notes/reviews/2026-10-12-sern-r2c-and-sweep-diagnose.md` の R2h(1)〜(5))。forge は回さない (0 step)。

R2a・R2f・R2g で直した「ファイル・列にあれば検査する」型の残り 5 箇所:
  R2h1 sern_gates.residual_health が check_convergence.analyze の判定不能を捨てて ok を組み直す
  R2h2 sern_gates.field_health・floor_gate が h5 にある roY* だけを見る
  R2h3 runner_sern3d.forces3d が出力の無い部品を足さずに通常の係数を返す
  R2h4 runner_sern3d.forces3d が 1 部品でも摩擦があれば他の部品の摩擦を 0 として C_T_with_shear を作る
  R2h5 runner_sern.restart_by_index が必須の状態量の欠落・寸法の不一致を飛ばす (標準起動の段の移行と初期化)

合格 (測る前に固定):
  (h1) R2c の 1 回目の session (run_1100、未定常だけで不合格 → 延長した回) の残差で、config・力・温度・床事象・保存場・最終場の床を
       固定し residual_health → evaluate_gates → judgment_extendable を通す。A (残差完全) は FAIL/JUDGMENT_UNSETTLED・延長可
       (保存された判定と同じ)、B (rms_roOmega 列だけ除く) は FAIL/JUDGMENT_UNVERIFIABLE・延長不可。変更前 (commit PRE) の
       実装では B も延長可だった (欠陥の再現)。A まで拒否したら方針 (プラトー受理) を誤って変えている。
       生の CSV (12 MB) は --r2c-raw (既定 = この試験を作ったセッションの scratchpad) から読み、無ければ止める (exit 2)。
       リポジトリには 40 行おきに間引いた fixture (notes の R2c の写し、residual_judgment1_stride40.csv.gz) を置き、試験の中で
       step-hold で 0..19999 に展開して使う (残差の時系列の完全性 §5.1 R2j を通すため。仕掛けは試験だけ)。生の CSV が
       あるときは列ごとの判定 (trend・合否) と A/B の結果が生の CSV と同じことを確かめる。--fixture で生の CSV なしに fixture だけで回す。
       合成の残差: CSV が無い・空・必須列の欠落・末尾窓の代表値 0 は判定不能 (延長不可)、NaN が重なれば DIVERGED が先、
       完全な列のプラトーは従来どおり合格 (require_residual_pass は使わない)
  (h2) 実効 config (solverConfig.yaml) の TP の化学種 roY0..roY{n−1} (トレーサは roXi) が最終 h5 に無い・ro と寸法が違う・config を
       読めないときは field_health が判定不能 (evaluate_gates は JUDGMENT_UNVERIFIABLE、延長不可)、floor_gate も不合格。
       全種の欠落 (変更前は組成の検査自体が消えて PASS)・欠けた種の質量分率が 0 (ΣρY = ρ は成り立つ) も拒否。CPG は要求しない。
       規則は runner_sern.required_state_keys と同じ関数 (config_species_keys)
  (h3) 部品の要求は形状 (prepare_info の wall_faces、無い旧 run は mesh = メッシャの info の条件) と bcondConfig から固定する:
       面の無い部品 (外側空間なしの側壁外・R4e の vehicle) は要求しない、面があり出力を要求する部品が 1 step でも欠ければ
       その標本は力係数を持たず incomplete → evaluate_gates は JUDGMENT_UNVERIFIABLE・延長不可 (変更前は部品を減らした係数で通る)。
       メッシャの条件 (wall_faces_from_minfo) はメッシャの境界 quad (6 構成) と手元の 3D run 6 本の sern.h5 の BCONDS に一致
  (h4) 摩擦を要求する部品 (面があり出力を要求するノズルの粘性壁) に twall が 1 部品でも無ければ C_T_with_shear・C_T_friction を作らず
       incomplete (変更前は他の部品の摩擦を 0 として足した)。slip 壁 (Euler) は摩擦を要求しない (従来どおり)
  (h5) restart_by_index(required=...) は継承で必須の量が継続元・継続先に無い・節点数が違えば**書き込みの前に**止める (継続先は不変)。
       inherit_state_keys: SST→SST は roK・roOmega を含み、層流 (none)→SST は含まない (継続先の初期場の乱流量を保つ)。
       run_staged の段の移行 (層流暖機 → soft → 本段) に配線され、SST の soft 段の res に roK が無ければ本段の前に止まり、層流段の res に
       roK が無いのは通る。check_initial_state は段の config が解く量が初期場に無ければ止める。warm_from_same_mesh は変更しない
  (s)  保存済みの run の判定が変わらない: case/46.sern_design の run (RUNS_LOCAL 10 本、--all で metrics.json と res_*.h5 のある
       全 run) で変更前後の field_health・floor_gate・residual_health・evaluate_gates (verdict・fail_class・reasons) が一致。3D の run
       (手元で壁面の出力がある 5 本) は力の履歴を新旧の経路で再計算して一致 (標本の step・全係数、incomplete なし)。
       R2c の保存系列 (2 session × 40 標本) と 3D の写し 12 本の保存系列は新しい部品・摩擦の検査を通り、保存された C_T_friction を
       保存された部品の twall から再計算した値と一致する (0 で補った標本が無い)

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_gate_completeness_tests.py [--r2c-raw DIR | --fixture] [--all]
"""
import argparse
import csv
import gzip
import inspect
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

import check_convergence as CC  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.meshing.mesh_sern3d import PHYS_SERN3D, SernMesh3DParams, generate_sern_mesh3d  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

DEFAULT_RAW = "/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad/r2c_raw"
RAW_SHA256 = "1679d4c919d83e211e095a833d3bd6664c76be1c78ffe47fbf08786f1f2c7155"     # 生の residual_judgment1.csv (12,964,240 バイト)
ap = argparse.ArgumentParser()
ap.add_argument("--r2c-raw", default=os.environ.get("SERN_R2C_RAW", DEFAULT_RAW),
                help="R2c の 1 回目の session の生の残差 residual_judgment1.csv を置いたディレクトリ")
ap.add_argument("--fixture", action="store_true", help="生の CSV なしに、間引いた fixture だけで (h1) を回す")
ap.add_argument("--all", action="store_true", help="(s) case/46.sern_design の metrics.json と res_*.h5 がある全 run も比べる")
args = ap.parse_args()

fails = 0
n_checks = 0
tmp = Path(tempfile.mkdtemp(prefix="sern_gate_completeness_"))
C46 = ROOT / "case" / "46.sern_design"
NOTES = ROOT / "notes" / "investigations"
R2C = NOTES / "2026-10-12-r2c" / "run_1100_r2c_lsw10_m6_judgment"
FIXTURE = R2C / "residual_judgment1_stride40.csv.gz"
PROD10 = C46 / "problem_3d_prod_3op_wallres_lswx10.yaml"
PRE = "b296ee8c"          # R2h の変更前 (R2g の後)


def check(name, ok, detail=""):
    global fails, n_checks
    n_checks += 1
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def raises(fn, exc=Exception):
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


def old_module(path, name, package):
    m = types.ModuleType(name)
    m.__package__ = package
    m.__file__ = str(ROOT / path)
    exec(compile(git_show(path), f"{PRE}:{path}", "exec"), m.__dict__)
    return m


# --- 変更前 (commit PRE) のモジュール -------------------------------------------------------------------------
OLD = old_module("design/forge_design/metrics/sern_gates.py", "forge_design.metrics._sern_gates_pre_r2h", "forge_design.metrics")
OLD3 = old_module("design/forge_design/evaluate/runner_sern3d.py", "forge_design.evaluate._runner_sern3d_pre_r2h", "forge_design.evaluate")
OLD2 = old_module("design/forge_design/evaluate/runner_sern.py", "forge_design.evaluate._runner_sern_pre_r2h", "forge_design.evaluate")
if hasattr(OLD, "run_species_keys") or "required" in inspect.signature(OLD2.restart_by_index).parameters \
        or "required_parts" in inspect.signature(OLD3.forces3d).parameters:
    stop(f"{PRE} のモジュールが変更前の形でない")


def load_hist(path):
    rows = list(csv.DictReader(open(path)))
    return [{k: (int(v) if k == "step" else float(v)) for k, v in r.items() if v not in ("", None)} for r in rows]


def stub_gates(mod, saved, keep=("residual_health", "residual_scale_gate")):
    """evaluate_gates の下のゲートのうち keep 以外を保存値で差し替える (戻り = 元に戻す関数)。"""
    vals = {"field_health": saved["field"], "residual_health": saved["residual"], "floor_gate": saved["floors"],
            "floor_event_gate": saved["floor_events"], "residual_scale_gate": saved["residual_scale"]}
    orig = {k: getattr(mod, k) for k in vals if k not in keep}
    for k in orig:
        setattr(mod, k, lambda *a, _v=vals[k], **k_: _v)

    def restore():
        for k, v in orig.items():
            setattr(mod, k, v)
    return restore


COLS7 = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"]
# 合成の run は forge と同じ時系列にする (outer_end の step 0..N−1 を 1 行ずつ、config に nStepOuter N と physProp、最終場 res_N):
# 2026-10-12 から受理ゲートは残差の時系列の欠測・終了位置の不一致と、化学種の列を決める physProp の無い config を判定不能にする
# (plan tooling-sern-te-wake-grid §5.1 R2i・R2j。旧 fixture は 100 step おき・nStepOuter・physProp なし)。試験の意図は変えない
NSTEP = 60
CFG_SST = f'physProp: {{thermalMethod: 0}}\nturbulence: {{model: "sst"}}\ntime:\n  last: {{nStepOuter: {NSTEP}}}\n'


def write_resid(d, cols=COLS7, n=NSTEP, pattern="plateau", zero_col=None, nan_col=None):
    with open(Path(d) / "residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["step", "inner", "phase"] + list(cols))
        for i in range(n):
            t = i / (n - 1)
            v = 10 ** (-1 - 1.5 * min(t * 3, 1)) if pattern == "plateau" else 10 ** (-1 - 2 * min(t * 3, 1)) * (1 + 40 * max(t - 0.8, 0))
            row = [i, -1, "outer_end"]
            for c in cols:
                x = v * (10 if c != "rms_ro" else 1)
                if c == zero_col and t > 0.5:
                    x = 0.0
                if c == nan_col and t > 0.9:
                    x = float("nan")
                row.append(x)
            w.writerow(row)


# ======================================================================== (h1) R2c の残差 A/B
SAVED1 = json.loads((R2C / "metrics_judgment1.json").read_text())
G1 = SAVED1["gates"]
P10 = load_problem(PROD10)
JP = R2.judgment_plan(P10.evaluate)
H1 = load_hist(R2C / "force_history_judgment1.csv")
MON1 = G.read_monitor_csv(R2C / "TE_MONITOR_judgment1.csv")


def r2c_ab(resid_rows, hdr, label):
    """生 (または fixture) の残差の行から A (完全)・B (rms_roOmega を除く) の run ディレクトリを作り、ゲートを通す。
    最終場は判定区間の段の終了 step の res_20000.h5 (中身は使わない: 場のゲートは保存値に固定し、残差の時系列の照合 [§5.1 R2j] が名前だけを見る)。"""
    out = {}
    for tag in ("A", "B"):
        d = tmp / f"r2c_{label}_{tag}"; d.mkdir()
        shutil.copy2(R2C / "solverConfig_judgment.yaml", d / "solverConfig.yaml")
        (d / "res_20000.h5").write_bytes(b"")
        keep = [i for i, c in enumerate(hdr) if not (tag == "B" and c == "rms_roOmega")]
        with open(d / "residual_history.csv", "w", newline="") as fh:
            w = csv.writer(fh); w.writerow([hdr[i] for i in keep])
            for r in resid_rows:
                w.writerow([r[i] for i in keep])
        res = {}
        for name, mod in (("new", G), ("old", OLD)):
            jg = mod.judgment_gate(H1, JP["window_steps"], JP["required"], JP["planned_steps"], monitor_rows=MON1)
            restore = stub_gates(mod, G1)
            try:
                g = mod.evaluate_gates(d, H1, 0, p_min=float(P10.evaluate.get("p_min", 1.0)), require_floor_events=True,
                                       floor_window_steps=JP["window_steps"], judgment=jg, model="sst")
            finally:
                restore()
            res[name] = (g, mod.judgment_extendable(g), jg)
        res["resid"] = G.residual_health(d)
        res["report"] = CC.analyze(str(d / "residual_history.csv"), 3.0, 0.2)
        out[tag] = res
    return out


def col_classes(report):
    """analyze の列ごとの (trend, 合否) と全体の旗 (ok, nan, stalled, converging)。"""
    import re
    cols = {c: ((re.search(r"\b(rising|falling|flat)\b", m) or [None, None])[1] if "drop=" in m else m.split("<--")[0].strip(), ok)
            for c, (m, ok) in report[1].items()}
    return cols, tuple(bool(x) for x in report[2:])


raw = Path(args.r2c_raw) / "residual_judgment1.csv"
with gzip.open(FIXTURE, "rt", newline="") as fh:
    fx_rows = list(csv.reader(fh))
fx_hdr, fx_body = fx_rows[0], fx_rows[1:]


def expand_stride(body, n_steps=20000):
    """40 行おきに間引いた fixture を、各行の値を次の行の step の手前まで保持して outer_end の step 0..N−1 に戻す (試験だけの仕掛け)。
    受理ゲートは残差の時系列の欠測を判定不能にする (§5.1 R2j) ので、間引いたままでは A (完全) の試験にならない。列ごとの判定 (trend・合否)
    が生の CSV と同じことは、生の CSV があるとき下で確かめる"""
    out = []
    for k, r in enumerate(body):
        s1 = int(body[k + 1][0]) if k + 1 < len(body) else n_steps
        out += [[str(s)] + r[1:] for s in range(int(r[0]), s1)]
    return out
runs_h1 = []
if raw.exists():
    import hashlib
    sha = hashlib.sha256(raw.read_bytes()).hexdigest()
    if sha != RAW_SHA256:
        stop(f"{raw} の sha256 {sha} が登録した R2c の 1 回目の残差 ({RAW_SHA256}) と違う")
    rows = list(csv.reader(open(raw)))
    runs_h1.append(("生の CSV (60000 行)", rows[0], rows[1:]))
elif not args.fixture:
    stop(f"R2c の 1 回目の生の残差 {raw} が無い (12 MB の CSV はリポジトリに置かない)。--r2c-raw で場所を渡すか、"
         f"--fixture で間引いた fixture ({FIXTURE.name}) だけで回す")
runs_h1.append((f"fixture ({FIXTURE.name}, 40 行おき {len(fx_body)} 行を step-hold で 20000 行に展開)", fx_hdr, expand_stride(fx_body)))
ab = {}
for label, hdr, body in runs_h1:
    key = "raw" if label.startswith("生") else "fx"
    ab[key] = r2c_ab(body, hdr, key)
    A, B = ab[key]["A"], ab[key]["B"]
    (gA, eA, jA), (gB, eB, _) = A["new"], B["new"]
    check(f"(h1) {label}: 判定区間の §6 #2・#3 は保存された 1 回目の判定と同じ (UNSETTLED / JUDGMENT_UNSETTLED)",
          jA["verdict"] == G1["judgment"]["verdict"] == "UNSETTLED" and jA["fail_class"] == G1["judgment"]["fail_class"],
          (jA["verdict"], jA["fail_class"]))
    check(f"(h1) {label} A (残差完全): 入力は妥当・プラトー (rms_roK だけ低下中) で ok、evaluate_gates は FAIL/JUDGMENT_UNSETTLED・"
          f"延長可 (保存された 1 回目と同じ。プラトー受理を変えていない)",
          A["resid"]["ok"] and A["resid"]["input"]["ok"] and not A["resid"]["unverifiable"]
          and A["resid"]["verdict"] == "NOT CONVERGED (stalled/plateau)" and A["resid"]["converging"] == ["rms_roK"]
          and gA["verdict"] == "FAIL" and gA["fail_class"] == "JUDGMENT_UNSETTLED" == G1["fail_class"] and eA[0] is True,
          (A["resid"]["verdict"], A["resid"]["input"], gA["fail_class"], eA))
    check(f"(h1) {label} B (rms_roOmega 列だけ除く): 判定不能 — evaluate_gates は FAIL/JUDGMENT_UNVERIFIABLE、延長不可 (理由に残差の入力)",
          not B["resid"]["ok"] and B["resid"]["unverifiable"] and any("rms_roOmega" in r for r in B["resid"]["input"]["reasons"])
          and gB["verdict"] == "FAIL" and gB["fail_class"] == "JUDGMENT_UNVERIFIABLE" and eB[0] is False and "残差の入力" in eB[1],
          (B["resid"]["input"], gB["fail_class"], eB))
    (gBo, eBo, _), (gAo, eAo, _) = B["old"], A["old"]
    check(f"(h1) {label}: 変更前 ({PRE}) の実装では B も A と同じく JUDGMENT_UNSETTLED・延長可だった (欠陥の再現)",
          gBo["fail_class"] == "JUDGMENT_UNSETTLED" and eBo[0] is True and gAo["fail_class"] == "JUDGMENT_UNSETTLED" and eAo[0] is True,
          (gBo["fail_class"], eBo))
    check(f"(h1) {label}: analyze は B を「必須の保存量残差列が無い: rms_roOmega」とし、他の列の trend は A と同じ (収束状態は別に判定)",
          "rms_roOmega" in B["report"][1].get("(入力)", ("",))[0]
          and {c: v for c, v in col_classes(B["report"])[0].items() if c != "(入力)"}
          == {c: v for c, v in col_classes(A["report"])[0].items() if c != "rms_roOmega"})
if "raw" in ab:
    rawA = ab["raw"]["A"]["report"]
    saved_txt = (R2C / "CONVERGENCE_VERDICT_judgment1.txt").read_text()
    same_saved = all(m.split("<--")[0].strip() in saved_txt for c, (m, _) in rawA[1].items())
    check("(h1) 生の CSV の analyze は保存された 1 回目の VERDICT (CONVERGENCE_VERDICT_judgment1.txt) の各列の行と一致", same_saved,
          [m for c, (m, _) in rawA[1].items() if m.split("<--")[0].strip() not in saved_txt])
    for tag in ("A", "B"):
        cr, cf = col_classes(ab["raw"][tag]["report"]), col_classes(ab["fx"][tag]["report"])
        gr, gf = ab["raw"][tag]["new"], ab["fx"][tag]["new"]
        check(f"(h1) fixture の {tag}: 列ごとの trend・合否と全体の旗、evaluate_gates の verdict・fail_class・延長可否が生の CSV と同じ",
              cr == cf and (gr[0]["verdict"], gr[0]["fail_class"], gr[1][0]) == (gf[0]["verdict"], gf[0]["fail_class"], gf[1][0]),
              (cr, cf))
else:
    print("    (生の CSV が無いので fixture だけで判定した。fixture が生の CSV と同じ判定になることは作成時に確かめた)")

# 合成の残差: 入力の不備の種類と優先順位
OK_SAVED = {"field": {"ok": True, "nonfinite": {}, "nonpositive": {}, "nan_files": [], "reasons": []}, "residual": None,
            "floors": {"ok": True, "counts": {}, "reasons": []},
            "floor_events": {"ok": True, "required": False, "verdict": None, "fail_class": None, "reasons": [], "notes": []},
            "residual_scale": {"ok": True, "reasons": []}}
STEADY_UNS = [{"step": 500 * (i + 1), "C_T": 0.95, "C_T_with_shear": 0.94 + 0.01 * i, "C_L": 0.1, "C_M": -1.0} for i in range(8)]


def gates_resid_only(mod, d, hist=None):
    """残差 (residual_health) だけを実際に評価し、他のゲートは合格に固定する。力の履歴は既定で未定常 (延長の対象になりうる形)。"""
    restore = stub_gates(mod, OK_SAVED)
    try:
        return mod.evaluate_gates(d, STEADY_UNS if hist is None else hist, 0, model="sst")
    finally:
        restore()


def resid_dir(name, cols=COLS7, cfg=CFG_SST, **kw):
    d = tmp / ("rh_" + name); d.mkdir()
    (d / f"res_{NSTEP}.h5").write_bytes(b"")      # 最終場 (中身は使わない。場のゲートは合格に固定)
    if cfg is not None:
        (d / "solverConfig.yaml").write_text(cfg)
    if cols is not None:
        write_resid(d, cols, **kw)
    return d


rh_ok = G.residual_health(resid_dir("ok"))
check("(h1) 合成: 完全な 7 列のプラトーは従来どおり ok (input.ok・stalled)", rh_ok["ok"] and rh_ok["input"]["ok"] and rh_ok["stalled"]
      and rh_ok["verdict"] == "NOT CONVERGED (stalled/plateau)", rh_ok)
for lab, d in (("CSV が無い", resid_dir("nocsv", cols=None)), ("空の CSV", None),
               ("rms_roUz の欠落", resid_dir("nouz", cols=[c for c in COLS7 if c != "rms_roUz"])),
               ("SST で rms_roK・rms_roOmega の欠落", resid_dir("nokw", cols=COLS7[:5])),
               ("設定が無く SST の列も無い (最も厳しい集合)", resid_dir("nocfg", cols=COLS7[:5], cfg=None)),
               ("列の末尾窓の代表値 0", resid_dir("zero", zero_col="rms_roe"))):
    if d is None:
        d = tmp / "rh_empty"; d.mkdir(); (d / "solverConfig.yaml").write_text(CFG_SST)
        (d / "residual_history.csv").write_text("step,inner,phase," + ",".join(COLS7) + "\n")
    rh, rho = G.residual_health(d), OLD.residual_health(d)
    g = gates_resid_only(G, d)
    check(f"(h1) 合成: {lab} → residual_health は判定不能 (input.ok False・unverifiable)、evaluate_gates (他のゲートは合格に固定) は "
          f"FAIL/JUDGMENT_UNVERIFIABLE、延長不可 (残差の入力)",
          not rh["ok"] and rh["unverifiable"] and not rh["input"]["ok"] and rh["input"]["reasons"]
          and g["fail_class"] == "JUDGMENT_UNVERIFIABLE" and G.judgment_extendable(g) == (False, "未定常以外の不合格: 残差の入力 (判定不能)"),
          (rh["verdict"], rh["input"], g["fail_class"], G.judgment_extendable(g)))
    if lab in ("rms_roUz の欠落", "SST で rms_roK・rms_roOmega の欠落", "列の末尾窓の代表値 0"):
        check(f"(h1) 合成: {lab} は変更前の実装では ok (欠陥の再現)", rho["ok"], rho)
d_nan = resid_dir("nan_and_missing", cols=COLS7[:5], nan_col="rms_ro")
g_nan = gates_resid_only(G, d_nan)
check("(h1) 合成: NaN と列の欠落が重なれば DIVERGED が先 (入力の不備だけのときに限り判定不能)",
      G.residual_health(d_nan)["nan"] and G.residual_health(d_nan)["unverifiable"] and g_nan["fail_class"] == "DIVERGED", g_nan["fail_class"])
d_rise = resid_dir("rise_and_missing", cols=COLS7[:5], pattern="rising")
rr = G.residual_health(d_rise)
check("(h1) 合成: rising と列の欠落が重なれば判定不能 (収束状態の rising は別に記録)",
      rr["rising"] and rr["unverifiable"] and rr["verdict"] == "UNVERIFIABLE (input)"
      and gates_resid_only(G, d_rise)["fail_class"] == "JUDGMENT_UNVERIFIABLE", rr["verdict"])
d_rise7 = resid_dir("rise_full", pattern="rising")
check("(h1) 合成: 完全な列の rising は従来どおり RESIDUAL_RISING (変更前と同じ判定)",
      gates_resid_only(G, d_rise7)["fail_class"] == "RESIDUAL_RISING" == gates_resid_only(OLD, d_rise7)["fail_class"]
      and not G.residual_health(d_rise7)["unverifiable"])
g_pl, g_plo = gates_resid_only(G, resid_dir("plateau_uns")), gates_resid_only(OLD, tmp / "rh_plateau_uns")
check("(h1) 合成: 完全な列のプラトー + 力の未定常は従来どおり UNSTEADY・延長可 (変更前と同じ)",
      g_pl["fail_class"] == "UNSTEADY" == g_plo["fail_class"] and G.judgment_extendable(g_pl)[0] and g_pl["reasons"] == g_plo["reasons"],
      (g_pl["fail_class"], G.judgment_extendable(g_pl)))
check("(h1) 合成: 欠落列はあっても require_residual_pass で隠さない (既定 False のまま判定不能)",
      inspect.signature(G.evaluate_gates).parameters["require_residual_pass"].default is False)
check("(h1) judgment_extendable: 保存済みの旧 gates (residual に input キーが無い) は従来どおり",
      G.judgment_extendable({"verdict": "FAIL", "rc": 0, "field": {"ok": True}, "floors": {"ok": True}, "floor_events": {"ok": True},
                             "residual_scale": {"ok": True}, "residual": {"ok": True, "verdict": "NOT CONVERGED (stalled/plateau)"},
                             "steadiness": {"ok": False, "series": {"C_T": {}}}})[0] is True)

# ======================================================================== (h2) 化学種の必須集合
CFG_TP = ('physProp: {thermalMethod: 2, species: [{name: "EXH", lump: {"N2": 1.0}, basis: mole}, '
          '{name: "AMB", lump: {"N2": 0.79, "O2": 0.21}, basis: mole}]}\nturbulence: {model: "sst"}\n')
CFG_TP_TR = CFG_TP.replace("basis: mole}]}", "basis: mole}], tracer: exhaust}")
CFG_CPG = 'physProp: {thermalMethod: 0}\nturbulence: {model: "sst"}\n'
N = 40


def tp_run(name, cfg=CFG_TP, drop=(), y=(0.3, 0.7), shape_bad=None, extra=None):
    d = tmp / ("tp_" + name); d.mkdir()
    if cfg is not None:      # 最終場 res_100.h5 と同じ session の長さ (§5.1 R2j)
        (d / "solverConfig.yaml").write_text(cfg + "time:\n  last: {nStepOuter: 100}\n")
    # 残差は forge が書く列: TP の 2 種は rms_roY0・rms_roY1、トレーサは rms_roXi (§5.1 R2i。場の欠落の試験で残差は完全にしておく)
    write_resid(d, COLS7 + (["rms_roY0", "rms_roY1"] if cfg in (CFG_TP, CFG_TP_TR) else []) + (["rms_roXi"] if cfg == CFG_TP_TR else []),
                n=100)
    with h5py.File(d / "res_100.h5", "w") as f:
        for k, v in (("ro", 1.0), ("roUx", 2.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", 3e5), ("P", 1e5), ("T", 300.0),
                     ("roK", 1.0), ("roOmega", 1e3), ("roY0", y[0]), ("roY1", y[1]), ("roXi", 0.5)):
            if k in drop:
                continue
            n = N - 1 if k == shape_bad else N
            f.create_dataset(f"VALUE/{k}", data=np.full(n, v))
        for k, v in (extra or {}).items():
            f.create_dataset(f"VALUE/{k}", data=v)
    return d


STEADY = [{"step": 500 * (i + 1), "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0} for i in range(8)]
d = tp_run("full")
fh, fg, g = G.field_health(d, model="sst"), G.floor_gate(d, model="sst"), G.evaluate_gates(d, STEADY, 0, model="sst")
check("(h2) TP の 2 種がそろう場は従来どおり (field ok・floor ok・PASS、変更前と field・floor・判定が同じ)",
      fh["ok"] and fh["species_required"] == ["roY0", "roY1"] and fg["ok"] and g["verdict"] == "PASS"
      and fg == OLD.floor_gate(d, model="sst") and fh["reasons"] == OLD.field_health(d, model="sst")["reasons"]
      and (g["verdict"], g["fail_class"], g["reasons"]) == tuple(OLD.evaluate_gates(d, STEADY, 0, model="sst")[k] for k in ("verdict", "fail_class", "reasons")),
      (fh, fg, g["reasons"]))
check("(h2) 規則は段の移行と同じ関数: required_state_keys(config) の化学種 = config_species_keys",
      R2.required_state_keys(CFG_TP) == ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"]
      and R2.required_state_keys(CFG_TP_TR)[-3:] == ["roY0", "roY1", "roXi"]
      and inspect.getsource(R2.required_state_keys).count("config_species_keys(doc)") == 1)
cases = [("roY1 の欠落", dict(drop=("roY1",)), ["roY1"]),
         ("全種の欠落 (変更前は組成の検査自体が消える)", dict(drop=("roY0", "roY1")), ["roY0", "roY1"]),
         ("欠けた種の質量分率が 0 (ΣρY = ρ は成り立つ)", dict(drop=("roY1",), y=(1.0, 0.0)), ["roY1"])]
for lab, kw, miss in cases:
    d = tp_run("m" + str(len(lab)) + miss[-1], **kw)
    fh, fg, g = G.field_health(d, model="sst"), G.floor_gate(d, model="sst"), G.evaluate_gates(d, STEADY, 0, model="sst")
    ext = G.judgment_extendable(g)
    check(f"(h2) {lab}: field_health は判定不能 (missing {miss})、floor_gate も不合格 (不在 -1)、evaluate_gates は "
          f"FAIL/JUDGMENT_UNVERIFIABLE・延長不可",
          not fh["ok"] and fh["missing"] == miss and not fg["ok"] and all(fg["counts"].get(f"{k} 不在") == -1 for k in miss)
          and g["fail_class"] == "JUDGMENT_UNVERIFIABLE" and not ext[0] and any("判定不能" in r for r in fg["reasons"]),
          (fh["missing"], fg, g["fail_class"], ext))
    if lab != "roY1 の欠落":
        go = OLD.evaluate_gates(d, STEADY, 0, model="sst")
        check(f"(h2) {lab}: 変更前の実装では PASS (欠陥の再現)", go["verdict"] == "PASS" and OLD.floor_gate(d, model="sst")["ok"], go["reasons"])
d = tp_run("shape", shape_bad="roY1")
fh, fg = G.field_health(d, model="sst"), G.floor_gate(d, model="sst")
check("(h2) roY1 の配列寸法が ro と違う: field_health は判定不能 (shape_mismatch)、floor_gate は寸法の違いを -1 で不合格 (例外にしない)",
      not fh["ok"] and fh["shape_mismatch"] == {"roY1": [N - 1]} and not fg["ok"] and any("配列寸法" in k for k in fg["counts"])
      and G.evaluate_gates(d, STEADY, 0, model="sst")["fail_class"] == "JUDGMENT_UNVERIFIABLE", (fh["reasons"], fg))
d = tp_run("nocfg", cfg=None)
fh, fg = G.field_health(d, model="sst"), G.floor_gate(d, model="sst")
check("(h2) 実効 config が無い: 必須の化学種を決められないので判定不能 (field・floor とも)",
      not fh["ok"] and fh["species_required"] is None and fh["unverifiable"] and not fg["ok"]
      and G.evaluate_gates(d, STEADY, 0, model="sst")["fail_class"] == "JUDGMENT_UNVERIFIABLE", (fh["reasons"], fg["reasons"]))
d = tp_run("tracer", cfg=CFG_TP_TR, drop=("roXi",))
check("(h2) トレーサ付き TP で roXi が無い → 判定不能", G.field_health(d, model="sst")["missing"] == ["roXi"])
d = tp_run("cpg", cfg=CFG_CPG, drop=("roY0", "roY1", "roXi"))
check("(h2) CPG (thermalMethod 0) は化学種を要求しない (従来どおり)",
      G.field_health(d, model="sst")["ok"] and G.floor_gate(d, model="sst")["ok"] and G.field_health(d, model="sst")["species_required"] == [])
d = tp_run("cpg_nan", cfg=CFG_CPG, extra={"roY5": np.full(N, np.nan)})
check("(h2) CPG でも h5 にある化学種は従来どおり検査する (NaN は DIVERGED)",
      G.field_health(d, model="sst")["nonfinite"].get("roY5") == N and G.evaluate_gates(d, STEADY, 0, model="sst")["fail_class"] == "DIVERGED")

# ======================================================================== (h3)(h4) 3D の力の部品と摩擦
P = PHYS_SERN3D


def sq(z0, z1, y=1.0):
    return [(0, y, z0), (1, y, z0), (1, y, z1), (0, y, z1)]


def quad_file(path, quads, ps, tau=None):
    xyz, conn = [], []
    for q in quads:
        base = len(xyz); xyz.extend(q); conn.extend([5, base, base + 1, base + 2, base + 3])
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=np.asarray(xyz, float).ravel()); f.create_dataset("MESH/CONNE", data=np.array(conn, np.int64))
        f.create_dataset("VALUE/Ps", data=np.repeat(np.asarray(ps, float), 4))
        if tau is not None:
            for c, v in zip("xyz", tau):
                f.create_dataset(f"VALUE/twall_{c}", data=np.full(4 * len(quads), v))


GEOM = {"ramp": ([sq(0, 0.5), sq(0.5, 1)], [3, 3]), "cowl_in": ([sq(0, 1, 0.0)], [4]), "cowl_out": ([sq(0, 1, -0.1)], [2]),
        "sidewall_in": ([[(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]], [2.5]), "sidewall_out": ([[(0, 0, 1.01), (1, 0, 1.01), (1, 1, 1.01), (0, 1, 1.01)]], [1.5]),
        "vehicle_top": ([sq(0, 1, 2.0)], [1.2]), "vehicle_side": ([[(0, 1, 1), (1, 1, 1), (1, 2, 1), (0, 2, 1)]], [1.1]),
        "vehicle_base": ([[(1, 1, 0), (1, 2, 0), (1, 2, 1), (1, 1, 1)]], [1.3])}
BC_LINE = "{name}: {{physID: {pid}, kind: {kind}, outputHDFflg: {out}, ints: , floats: }}\n"


def force_run(name, faces, kinds=None, outs=None, steps=(500, 1000), skip=(), no_tau=(), wall_faces=True, minfo=None):
    """合成の 3D run: prepare_info (wall_faces か mesh)、bcondConfig、面のある部品の壁面出力。skip = 出さない (部品, step)。"""
    d = tmp / ("f_" + name); d.mkdir()
    info = {"mesh": minfo or {}}
    if wall_faces:
        info[R3.WALL_FACES_KEY] = {n: faces.get(n, 0) for n in R3._FORCE_PARTS}
    (d / "prepare_info.json").write_text(json.dumps(info))
    bc = ""
    for n, (pid, _) in R3._FORCE_PARTS.items():
        kind = (kinds or {}).get(n, "wall")
        bc += BC_LINE.format(name=n, pid=pid, kind=kind, out=(outs or {}).get(n, 1))
    (d / "bcondConfig.yaml").write_text(bc)
    for s in steps:
        for n, cnt in faces.items():
            if cnt and (n, s) not in skip and (n, "*") not in skip:
                q, ps = GEOM[n]
                tau = None if (n in no_tau or (kinds or {}).get(n) == "slip") else (0.01, 0.0, 0.0)
                quad_file(d / f"res_{n}_{P[n]}_{s}.h5", q, ps, tau)
    return d


KW3 = dict(p_a=1.0, F_ideal_per_m=10.0, half_W=1.0, H=1.0)
FACES8 = {"ramp": 2, "cowl_in": 1, "cowl_out": 1, "sidewall_in": 1, "sidewall_out": 1, "vehicle": 0, "vehicle_top": 1,
          "vehicle_side": 1, "vehicle_base": 1}


def hist_new(d):
    fp = R3.force_parts_plan(d)
    steps = R3._part_steps(d, fp["required"] or ["ramp"])
    return fp, [R3.forces3d(d, s, **KW3, required_parts=fp["required"], friction_parts=fp["friction"]) for s in steps]


def hist_old(d):
    steps = sorted({int(f.stem.split("_")[-1]) for f in Path(d).glob("res_ramp_*_[0-9]*.h5")})
    return [OLD3.forces3d(d, s, **KW3) for s in steps]


def coeffs(h):
    return {k: v for k, v in h.items() if k.startswith("C_") or k in ("step", "T_wall", "L", "M_noseup", "F_gross")}


d = force_run("full", FACES8)
fp, hn = hist_new(d)
ho = hist_old(d)
check("(h3) 8 部品の構成: 要求 = 面があり出力を要求する 8 部品 (vehicle は 0 面で要求しない)、摩擦はノズルの粘性壁 5 部品",
      fp["required"] == ["ramp", "cowl_in", "cowl_out", "sidewall_in", "sidewall_out", "vehicle_top", "vehicle_side", "vehicle_base"]
      and fp["friction"] == ["ramp", "cowl_in", "cowl_out", "sidewall_in", "sidewall_out"] and not fp["reasons"], fp)
check("(h3) 全部品がそろう標本は変更前と同じ係数 (incomplete なし)", [coeffs(h) for h in hn] == [coeffs(h) for h in ho]
      and not any("incomplete" in h for h in hn) and all("C_T_with_shear" in h for h in hn))
d = force_run("miss_cowl", FACES8, skip={("cowl_out", 1000)})
fp, hn = hist_new(d)
ho = hist_old(d)
g = G.evaluate_gates(tmp / "nothing", hn, 0, model="sst")
check("(h3) step 1000 で cowl_out の出力が無い → その標本は力係数を持たず incomplete、evaluate_gates は JUDGMENT_UNVERIFIABLE・"
      "延長不可 (力の部品の出力)",
      [("C_T" in h, bool(h.get("incomplete"))) for h in hn] == [(True, False), (False, True)]
      and g["force_parts"]["n_incomplete"] == 1 and not g["force_parts"]["ok"]
      and "力の部品の出力" in G.judgment_extendable(dict(g, rc=0))[1] and not G.judgment_extendable(dict(g, rc=0))[0],
      (hn[-1].get("incomplete"), g["force_parts"]))
check("(h3) 同じ入力で変更前の forces3d は cowl_out を足さない係数を返していた (欠陥の再現)",
      len(ho) == 2 and "C_T" in ho[1] and ho[1]["C_L"] != ho[0]["C_L"] and "cowl_out" not in ho[1]["parts"])
g_full = G.evaluate_gates(tmp / "nothing", hist_new(force_run("full2", FACES8))[1], 0, model="sst")
check("(h3) 全部品がそろう履歴は力の部品のゲートを通る (force_parts ok)", g_full["force_parts"]["ok"] and g_full["force_parts"]["n_incomplete"] == 0)
d = force_run("miss_ramp", FACES8, skip={("ramp", 1000)})
fp, hn = hist_new(d)
check("(h3) ランプの出力が無い step も落とさない (標本の step は要求する全部品の出力の和集合。変更前はランプの出力の step だけ)",
      [h["step"] for h in hn] == [500, 1000] and bool(hn[1].get("incomplete")) and [h["step"] for h in hist_old(d)] == [500])
F_NOOUT = dict(FACES8, sidewall_out=0, vehicle=0, vehicle_top=0, vehicle_side=0, vehicle_base=0)
d = force_run("noouter", F_NOOUT)
fp, hn = hist_new(d)
check("(h3) 面の無い部品 (外側空間なしの側壁外・機体面) は bcond が出力を要求しても要求しない (存在しない面と出力の欠落を区別)",
      fp["required"] == ["ramp", "cowl_in", "cowl_out", "sidewall_in"] and not any(h.get("incomplete") for h in hn), fp)
d = force_run("nozzle_noout", FACES8, outs={"cowl_out": 0})
fp = R3.force_parts_plan(d)
check("(h3) ノズルの部品 (cowl_out) に面があるのに bcond の出力が無い → 判定不能の理由 (全標本の incomplete)",
      any("cowl_out" in r and "判定不能" in r for r in fp["reasons"]), fp["reasons"])
d = force_run("veh_noout", FACES8, outs={"vehicle_top": 0})
check("(h3) 機体の面 (vehicle_top) の出力が無い設定は要求しない (ノズル力の帳簿の外。従来どおり)",
      "vehicle_top" not in R3.force_parts_plan(d)["required"] and not R3.force_parts_plan(d)["reasons"])
d = force_run("pid_bad", FACES8)
(d / "bcondConfig.yaml").write_text((d / "bcondConfig.yaml").read_text().replace("cowl_in: {physID: 5,", "cowl_in: {physID: 6,"))
check("(h3) bcond の physID が部品のタグと違う → 判定不能", any("physID" in r for r in R3.force_parts_plan(d)["reasons"]))
d = force_run("nobc", FACES8); (d / "bcondConfig.yaml").unlink()
check("(h3) bcondConfig.yaml が無い → 判定不能", any("bcondConfig" in r for r in R3.force_parts_plan(d)["reasons"]))
MINFO_OLD = {"nz": 39, "k_sw": 24}         # 機体面の記録の無い旧メッシャの info
d = force_run("minfo_old", {k: v for k, v in FACES8.items() if k in R3._NOZZLE_FACES}, wall_faces=False, minfo=MINFO_OLD,
              outs={"vehicle": 1})
check("(h3) wall_faces の無い旧 run: メッシャの info から面を決め、面の有無を決められない機体面に出力を要求していれば判定不能",
      R3.force_parts_plan(d)["faces_source"].startswith("prepare_info.mesh")
      and any("vehicle" in r and "決められない" in r for r in R3.force_parts_plan(d)["reasons"]), R3.force_parts_plan(d))
d = force_run("tau", FACES8, no_tau=("cowl_in",), steps=tuple(500 * (i + 1) for i in range(8)))     # 定常性の判定に 4 標本以上
fp, hn = hist_new(d)
ho = hist_old(d)
check("(h4) 粘性壁の cowl_in に twall が無い → C_T_with_shear・C_T_friction を作らず incomplete (圧力の係数は作る)",
      all("C_T_with_shear" not in h and "C_T_friction" not in h and "C_T" in h and h.get("incomplete") for h in hn), hn[0].get("incomplete"))
check("(h4) 同じ入力で変更前の forces3d は cowl_in の摩擦を 0 として C_T_with_shear を作っていた (欠陥の再現)",
      all("C_T_with_shear" in h for h in ho))
g = gates_resid_only(G, tmp / "rh_ok", hist=hn)
check("(h4) evaluate_gates (SST、他のゲートは合格に固定): FAIL/JUDGMENT_UNVERIFIABLE・延長不可 (力の部品の出力・C_T_with_shear の列)",
      g["fail_class"] == "JUDGMENT_UNVERIFIABLE" and not G.judgment_extendable(g)[0] and "力の部品の出力" in G.judgment_extendable(g)[1]
      and g["steadiness"]["missing"] == ["C_T_with_shear"], (g["fail_class"], g["reasons"][:3], G.judgment_extendable(g)))
g_o = gates_resid_only(OLD, tmp / "rh_ok", hist=hist_old(d))
check("(h4) 同じ入力で変更前のゲートは C_T_with_shear (摩擦 0 補完) で PASS していた (欠陥の再現)", g_o["verdict"] == "PASS", g_o["reasons"][:3])
d = force_run("tau_veh", FACES8, no_tau=("vehicle_top",))
fp, hn = hist_new(d)
check("(h4) 機体の面だけ twall が無い → C_T_with_shear は作る (目的量に寄与しない面は要求しない)",
      all("C_T_with_shear" in h and not h.get("incomplete") for h in hn))
d = force_run("slip", FACES8, kinds={n: "slip" for n in R3._FORCE_PARTS})
fp, hn = hist_new(d)
check("(h4) slip 壁 (Euler) は摩擦を要求せず、twall が無いので C_T_with_shear も作らない (従来どおり)",
      fp["friction"] == [] and all("C_T_with_shear" not in h and not h.get("incomplete") for h in hn)
      and [coeffs(h) for h in hn] == [coeffs(h) for h in hist_old(d)], fp)

# collect が要求を forces3d に渡し、要求を固定できない run は全標本に理由を付ける
seen = []


def _capture(run_dir, hist, rc, **kw):
    seen.append(hist)
    return {"steadiness": {"series": {}}, "objective": "C_T_with_shear", "verdict": "FAIL", "fail_class": "X", "reasons": []}


INFO3 = {"states": {"exhaust": {"ro": 1.0, "u": 1.0, "P": 1.0}, "ext": {"P": 1.0}}, "H_m": 1.0, "F_ideal_N_per_m": 10.0,
         "half_W_m": 1.0, "moc_forces": {}, "discretization": "node", "mesh": {}}
PY3 = tmp / "p3.yaml"
PY3.write_text(yaml.safe_dump({"name": "t", "type": "sern_2d", "gas": {}, "spec": {}, "dv": {}, "geometry": {}, "mesh": {},
                               "evaluate": {"model": "sst"}}))
_saved = R3.evaluate_gates
try:
    R3.evaluate_gates = _capture
    for lab, kw in (("全部品", {}), ("cowl_out 欠落", {"skip": {("cowl_out", 1000)}})):
        d = force_run("collect_" + str(len(lab)), FACES8, **kw)
        info = json.loads((d / "prepare_info.json").read_text()); info.update(INFO3)
        (d / "prepare_info.json").write_text(json.dumps(info))
        out = R3.collect(PY3, d, out_dir=tmp / ("co_" + d.name), rc=0)
        h = seen[-1]
        fp_full = R3.force_parts_plan(d)["required"]
        check(f"(h3) collect ({lab}): 部品の要求を forces3d に渡し、metrics の force_parts に要求を残す",
              bool(fp_full) and out["force_parts"]["required"] == fp_full
              and [bool(x.get("incomplete")) for x in h] == ([False, False] if not kw else [False, True]), (out.get("force_parts"), [x.get("incomplete") for x in h]))
    d = force_run("collect_nobc", FACES8); (d / "bcondConfig.yaml").unlink()
    info = json.loads((d / "prepare_info.json").read_text()); info.update(INFO3); (d / "prepare_info.json").write_text(json.dumps(info))
    R3.collect(PY3, d, out_dir=tmp / "co_nobc", rc=0)
    check("(h3) collect: 要求を固定できない run (bcond なし) は全標本に判定不能の理由", seen[-1] and all(x.get("incomplete") for x in seen[-1]))
finally:
    R3.evaluate_gates = _saved

# メッシャの条件 (wall_faces_from_minfo) がメッシャの境界 quad と一致する
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec  # noqa: E402
k_ = PlanarMOC(SernKernelSpec(M_in=2.5, theta_r0=np.deg2rad(15.0), theta_c0=np.deg2rad(5.0), L_cowl=1.0,
                              x_max=6.0, nj=201, dx=4e-3, p_ext_over_p_in=0.05)).march()
dz = k_.design_ramp(M_c=3.3, f=0.45)
BASEP = dict(ni_up=6, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=6, W=2.0, Z_ext=1.5,
             interface_angle=float(k_.TH[-1, 0]), top_ext_angle=dz.info["theta_e"])
EXT = dict(ext_top=True, nj_ext_top=9, top_depth=1.5, vehicle_taper=0.3, vehicle_wedge_deg=3.0)
for lab, kw in (("既定 (機体面なし・W_vehicle なし)", {}), ("W_vehicle 3.0", {"W_vehicle": 3.0}), ("外側空間なし", {"nz_out": 0}),
                ("ext_top + 機体側面", dict(EXT, W_vehicle=3.0)), ("ext_top・機体側面なし", dict(EXT, W_vehicle=3.0, vehicle_side=False)),
                ("ext_top・外側空間なし", dict(EXT, nz_out=0))):
    _, _, B_, info_, _ = generate_sern_mesh3d(dz, SernMesh3DParams(**{**BASEP, **kw}))
    truth = {n: len(B_.get(n, [])) > 0 for n in R3._FORCE_PARTS}
    rule = R3.wall_faces_from_minfo(info_)
    check(f"(h3) メッシャ ({lab}): wall_faces_from_minfo = 境界 quad の有無 {sorted(n for n, v in truth.items() if v)}", rule == truth,
          {n: (rule[n], truth[n]) for n in truth if rule[n] != truth[n]})


def _main_case_dirs():
    yield ROOT / "case"
    r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--path-format=absolute", "--git-common-dir"], capture_output=True, text=True)
    if r.returncode == 0 and r.stdout.strip():
        yield Path(r.stdout.strip()).parent / "case"


def find_run(rel, need=("prepare_info.json",)):
    for b in _main_case_dirs():
        if all((b / rel / n).exists() for n in need):
            return b / rel
    return None


RUNS3D = ["run_0023_smoke_sst_node_3d_noouter", "run_0027_diag3d_euler_accel_slip", "run_0082_3d_sst_lsw08", "run_0121_r5_vehicleside_lsw08",
          "run_0122_r5_outflow", "run_0401_3d_wf0"]
for name in RUNS3D:
    rd = find_run("46.sern_design/" + name, ("prepare_info.json", "sern.h5", "bcondConfig.yaml"))
    if rd is None:
        stop(f"case/46.sern_design/{name} (prepare_info.json・sern.h5・bcondConfig.yaml) が見つからない")
    info = json.loads((rd / "prepare_info.json").read_text())
    rule = R3.wall_faces_from_minfo(info["mesh"])
    with h5py.File(rd / "sern.h5", "r") as f:
        have = {int(x) for x in f["BCONDS"].keys()}
    bad = {n: (rule[n], R3._FORCE_PARTS[n][0] in have) for n in R3._FORCE_PARTS if rule[n] is not None and rule[n] != (R3._FORCE_PARTS[n][0] in have)}
    check(f"(h3) {name}: wall_faces_from_minfo (決められた部品 {sum(v is not None for v in rule.values())}/9) が sern.h5 の BCONDS と一致",
          not bad, bad)

# ======================================================================== (h5) 初期化と段の移行
SST7 = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"]
CFG_NONE = 'physProp: {thermalMethod: 0}\nturbulence: {model: "none"}\n'
CFG_SSTc = 'physProp: {thermalMethod: 0}\nturbulence: {model: "sst"}\n'
check("(h5) inherit_state_keys: SST→SST は roK・roOmega を含み、層流 (none)→SST は含まない、TP は化学種を含む",
      R2.inherit_state_keys(CFG_SSTc, CFG_SSTc) == SST7 and R2.inherit_state_keys(CFG_NONE, CFG_SSTc) == SST7[:5]
      and R2.inherit_state_keys(CFG_TP, CFG_TP)[-2:] == ["roY0", "roY1"] and R2.inherit_state_keys(CFG_SSTc, CFG_NONE) == SST7[:5])
rng = np.random.default_rng(5)


def h5_state(path, keys, n=30, dtype=np.float64):
    with h5py.File(path, "w") as f:
        for i, k in enumerate(keys):
            f.create_dataset(f"VALUE/{k}", data=(rng.standard_normal(n if not isinstance(n, dict) else n.get(k, 30)) + i).astype(dtype))
        f.create_dataset("VALUE/wall_dist", data=np.ones(30))


def rdir(name, src_keys, dst_keys=SST7, cfg=CFG_SSTc, n=30):
    d = tmp / ("rs_" + name); d.mkdir()
    (d / "solverConfig.yaml").write_text(cfg)
    h5_state(d / "res_100.h5", src_keys, n=n)
    h5_state(d / "sern.h5", dst_keys)
    return d


def digest(path):
    with h5py.File(path, "r") as f:
        return {k: f["VALUE"][k][()].tobytes() for k in f["VALUE"]}


d = rdir("sst_nok", SST7[:5] + ["roOmega"])
before = digest(d / "sern.h5")
ok_stop = raises(lambda: R2.restart_by_index(d / "res_100.h5", d / "sern.h5", required=R2.inherit_state_keys(CFG_SSTc, CFG_SSTc)), ValueError)
check("(h5) SST→SST で継続元に roK が無い → 書き込みの前に止まる (継続先はビット不変)", ok_stop and digest(d / "sern.h5") == before)
d2 = rdir("sst_nok_old", SST7[:5] + ["roOmega"])
before2 = digest(d2 / "sern.h5")
OLD2.restart_by_index(d2 / "res_100.h5", d2 / "sern.h5")
check("(h5) 同じ入力で変更前の restart_by_index は roK を飛ばして通っていた (欠陥の再現)",
      digest(d2 / "sern.h5")["roK"] == before2["roK"] and digest(d2 / "sern.h5")["ro"] != before2["ro"])
d = rdir("sst_short", SST7, n={"roK": 29})
before = digest(d / "sern.h5")
check("(h5) 継続元の roK の節点数が違う → 止まる (変更前は黙って飛ばした)",
      raises(lambda: R2.restart_by_index(d / "res_100.h5", d / "sern.h5", required=SST7), ValueError) and digest(d / "sern.h5") == before)
d = rdir("ro_short", SST7, n={k: 29 for k in SST7})
before = digest(d / "sern.h5")
check("(h5) required を渡さない呼び出しでも流体 5 量の節点数違い (全量を黙って飛ばしていた) は止まる",
      raises(lambda: R2.restart_by_index(d / "res_100.h5", d / "sern.h5"), ValueError) and digest(d / "sern.h5") == before)
d = rdir("dst_nok", SST7, dst_keys=SST7[:5])
check("(h5) 継続先に roK が無い (SST→SST) → 止まる", raises(lambda: R2.restart_by_index(d / "res_100.h5", d / "sern.h5", required=SST7), ValueError))
d = rdir("lam", SST7[:5], cfg=CFG_NONE)
before = digest(d / "sern.h5")
R2.restart_by_index(d / "res_100.h5", d / "sern.h5", required=R2.inherit_state_keys(CFG_NONE, CFG_SSTc))
after = digest(d / "sern.h5")
check("(h5) 層流 (none)→SST: 継続元に roK・roOmega が無くても通り、継続先の初期場の乱流量はビット不変、流体 5 量は写る",
      after["roK"] == before["roK"] and after["roOmega"] == before["roOmega"] and all(R2._restart_mismatch(d / "res_100.h5", d / "sern.h5", SST7[:5]) == [] for _ in [0]))
d = rdir("full", SST7)
R2.restart_by_index(d / "res_100.h5", d / "sern.h5", required=SST7)
check("(h5) SST→SST で全部そろう → 従来どおり写る (段の移行の検査 _restart_mismatch も空)", R2._restart_mismatch(d / "res_100.h5", d / "sern.h5", SST7) == [])
d = rdir("init_ok", SST7)
check("(h5) check_initial_state: SST の config で初期場に 7 量がそろう → 通る", R2.check_initial_state(d / "sern.h5", CFG_SSTc) == SST7)
d = rdir("init_nok", SST7, dst_keys=SST7[:5])
check("(h5) check_initial_state: SST の config で初期場に roK・roOmega が無い (write_ic_arrays が飛ばす) → 止まる",
      raises(lambda: R2.check_initial_state(d / "sern.h5", CFG_SSTc), ValueError)
      and R2.check_initial_state(d / "sern.h5", CFG_NONE) == SST7[:5])
d = rdir("init_tp", SST7, cfg=CFG_TP)
check("(h5) check_initial_state: TP の config で初期場に roY0・roY1 が無い → 止まる", raises(lambda: R2.check_initial_state(d / "sern.h5", CFG_TP), ValueError))
import ast  # noqa: E402


def fn_ast(src, name):
    return next(ast.dump(n) for n in ast.walk(ast.parse(src)) if isinstance(n, ast.FunctionDef) and n.name == name)


check("(h5) warm_from_same_mesh (keep_turb=True、Euler からの暖機) は変更しない (関数の構文木が変更前と同じ)",
      fn_ast((ROOT / "design/forge_design/evaluate/runner_sern3d.py").read_text(), "warm_from_same_mesh")
      == fn_ast(git_show("design/forge_design/evaluate/runner_sern3d.py"), "warm_from_same_mesh"))

# run_staged の配線: 層流暖機 → soft → 本段 (forge は偽物: 段の config の乱流モデルに応じた res を書く)
CFG_MAIN = ('physProp: {thermalMethod: 0}\nturbulence: {model: "sst", scalarDiffusion: 1}\nspace: {convMethod: 1}\n'
            'time:\n  last: {nStepOuter: 40}\n  deltaT: {cfl: 0.25, cfl_pseudo: 0.25}\n  outStepInterval: 10\n')


def staged_dir(name):
    d = tmp / ("st_" + name); d.mkdir()
    (d / "solverConfig.yaml").write_text(CFG_MAIN); (d / "solverConfig_main.yaml").write_text(CFG_MAIN)
    h5_state(d / "sern.h5", SST7)
    return d


def fake_forge_factory(drop_k_in_sst_soft):
    calls = []

    def fake(run_dir):
        rd = Path(run_dir); cfg = (rd / "solverConfig.yaml").read_text(); calls.append(cfg)
        lam = 'model: "none"' in cfg
        keys = SST7[:5] if (lam or (drop_k_in_sst_soft and "convMethod: 0" in cfg)) else SST7
        h5_state(rd / f"res_{10 * len(calls)}.h5", keys)
        return 0
    return fake, calls


_saved_rf = R2.run_forge
try:
    fake, calls = fake_forge_factory(False)
    R2.run_forge = fake
    d = staged_dir("ok")
    rc = R2.run_staged(d, "full", soft_steps=20, warm_lam_steps=20)
    check("(h5) run_staged: 層流暖機 (res に roK なし) → soft (SST) → 本段は通る (層流 → SST は乱流量を継承しない)",
          rc == 0 and len(calls) == 3 and 'model: "none"' in calls[0] and 'model: "sst"' in calls[1], [c[:60] for c in calls])
    fake, calls = fake_forge_factory(True)
    R2.run_forge = fake
    d = staged_dir("nok")
    ok_stop = raises(lambda: R2.run_staged(d, "full", soft_steps=20, warm_lam_steps=20), ValueError)
    check("(h5) run_staged: SST の soft 段の res に roK が無い → 本段の前に止まる (SST→SST の継承は必須)",
          ok_stop and len(calls) == 2, len(calls))
finally:
    R2.run_forge = _saved_rf

# ======================================================================== (s) 保存済みの run の判定が変わらない
RUNS_LOCAL = ["run_0100_species_regress_euler_lumped", "run_0101_species_regress_euler_full", "run_0102_species_regress_euler_lumped_24k",
              "run_0123_cfl_A_base", "run_0128_cflR1_0p5", "run_0163_rankB_2000", "run_0198_rb_prop_A", "run_0308_wake_soft05",
              "run_0121_r5_vehicleside_lsw08", "run_0122_r5_outflow"]
FIELD_KEYS = ("ok", "reasons", "nonfinite", "nonpositive", "min", "species_checked", "file", "nan_files", "missing")
RESID_KEYS = ("ok", "converged", "verdict", "nan", "rising", "stalled", "converging", "columns", "last_step", "reasons")


def compare_run(d):
    info = json.loads((d / "prepare_info.json").read_text())
    model = info.get("model")
    if model not in ("sst", "euler"):
        return None, f"prepare_info.json の model {model!r}"
    met = json.loads((d / "metrics.json").read_text())
    saved = met.get("gates") or {}
    pk = [k for k in ((saved.get("floors") or {}).get("counts") or {}) if k.startswith("P<=")]
    pmin = float(pk[0][3:]) if pk else None
    fn, fo = G.field_health(d, model=model), OLD.field_health(d, model=model)
    gn, go = G.floor_gate(d, pmin, model=model), OLD.floor_gate(d, pmin, model=model)
    rn, ro = G.residual_health(d), OLD.residual_health(d)
    if int(info.get("dim", 2)) == 3:
        fp = R3.force_parts_plan(d, info)
        hist = [R3.forces3d(d, s, ((info["states"])["ext"])["P"], info["F_ideal_N_per_m"], info["half_W_m"], info["H_m"],
                            mdot_u_in=info["states"]["exhaust"]["ro"] * info["states"]["exhaust"]["u"] ** 2 * info["H_m"],
                            p_in=info["states"]["exhaust"]["P"], twall_on_fluid=(info["discretization"] == "cell"),
                            half_W_vehicle=info.get("half_W_vehicle_m"), required_parts=fp["required"], friction_parts=fp["friction"])
                for s in R3._part_steps(d, fp["required"] or ["ramp"])]
        steps_old = sorted({int(f.stem.split("_")[-1]) for f in d.glob("res_ramp_*_[0-9]*.h5")})
        hist_o = [OLD3.forces3d(d, s, info["states"]["ext"]["P"], info["F_ideal_N_per_m"], info["half_W_m"], info["H_m"],
                                mdot_u_in=info["states"]["exhaust"]["ro"] * info["states"]["exhaust"]["u"] ** 2 * info["H_m"],
                                p_in=info["states"]["exhaust"]["P"], twall_on_fluid=(info["discretization"] == "cell"),
                                half_W_vehicle=info.get("half_W_vehicle_m")) for s in steps_old]
        same_hist = [coeffs(h) for h in hist] == [coeffs(h) for h in hist_o] and not any(h.get("incomplete") for h in hist) and not fp["reasons"]
    else:
        hist = hist_o = load_hist(d / "force_history.csv") if (d / "force_history.csv").exists() else []
        same_hist, fp = True, None
    kw = dict(require_residual_pass=bool(saved.get("require_residual_pass")), require_residual_plateau=bool(saved.get("require_residual_plateau")),
              p_min=pmin)
    rc = met.get("forge_rc", saved.get("rc"))
    en, eo = G.evaluate_gates(d, hist, rc, model=model, **kw), OLD.evaluate_gates(d, hist_o, rc, model=model, **kw)
    parts = {"field": {k: fn.get(k) for k in FIELD_KEYS} == {k: fo.get(k) for k in FIELD_KEYS}, "floor": gn == go,
             # columns の drop は判定できない列で NaN (NaN != NaN) なので JSON の文字列で比べる
             "residual": json.dumps({k: rn.get(k) for k in RESID_KEYS}, sort_keys=True) == json.dumps({k: ro.get(k) for k in RESID_KEYS}, sort_keys=True)
             and rn["input"]["ok"],
             "gates": (en["verdict"], en["fail_class"], en["reasons"]) == (eo["verdict"], eo["fail_class"], eo["reasons"]), "hist": same_hist}
    same = all(parts.values())
    return (same, model, info.get("dim", 2), en["verdict"], en["fail_class"], saved.get("verdict"), saved.get("fail_class"),
            len(hist), (fp or {}).get("required"), [k for k, v in parts.items() if not v]), None


names = list(RUNS_LOCAL) + [n for n in RUNS3D if n not in RUNS_LOCAL]
if args.all:
    seen_n = set(names)
    for b in _main_case_dirs():
        for m in sorted((b / "46.sern_design").glob("run_*/metrics.json")):
            if m.parent.name not in seen_n and list(m.parent.glob("res_[0-9]*.h5")) and (m.parent / "prepare_info.json").exists():
                seen_n.add(m.parent.name); names.append(m.parent.name)
summary, n_cmp, diff = {}, 0, []
# この比較は R2h の変更前後 (PRE) を見る。後から入れた残差の時系列の完全性 (§5.1 R2j) は外す: 発散で切り詰められた保存 run
# (run_0121 など) は R2j で理由が 1 つ増えるので、R2h の比較に混ぜると R2h の不変性を確かめられない。R2j の保存 run への影響は
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
    n_cmp += 1
    summary[(res[1], res[2], res[3])] = summary.get((res[1], res[2], res[3]), 0) + 1
    if not res[0]:
        diff.append((name, res))
    if name in RUNS_LOCAL or name in RUNS3D or not res[0]:
        check(f"(s) case/46 {name} [{res[1]} {res[2]}D]: 変更前後で field・floor・residual・evaluate_gates が一致"
              + (f"、力の履歴 {res[7]} 標本 (要求 {len(res[8] or [])} 部品) も新旧の経路で一致" if res[2] == 3 else "")
              + f" (verdict {res[3]}/{res[4]}; 保存値 {res[5]}/{res[6]})", res[0], res)
G.residual_series_check = _series_check
if args.all:
    check(f"(s) --all: case/46 の {n_cmp} 本で一致 (内訳 {sorted(summary.items())})", not diff, diff[:5])


# R2c の保存系列と 3D の写し: 新しい部品・摩擦の検査を通り、保存された C_T_friction を保存された部品から再計算できる
def series_ok(hist, req, fric, label):
    bad_parts = [(h["step"], sorted(set(req) - set(h["parts"]))) for h in hist if set(req) - set(h["parts"])]
    bad_tau = [(h["step"], [n for n in fric if "Fx_tau" not in h["parts"][n]]) for h in hist if any("Fx_tau" not in h["parts"].get(n, {}) for n in fric)]
    rel = []
    for h in hist:
        ft = -sum(v.get("Fx_tau", 0.0) for n, v in h["parts"].items() if n not in R3._VEHICLE_FACES) / h["F_ideal"]   # node: sgn = −1
        rel.append(abs(ft - h["C_T_friction"]) / max(abs(h["C_T_friction"]), 1e-30))
    return bad_parts, bad_tau, max(rel) if rel else None


fp_r2c = R3.force_parts_plan(R2C)
REQ8 = ["ramp", "cowl_in", "cowl_out", "sidewall_in", "sidewall_out", "vehicle_top", "vehicle_side", "vehicle_base"]
check("(s) R2c の写し (prepare_info・bcondConfig): 要求 8 部品 (vehicle は 0 面)・摩擦 5 部品、面の有無はメッシャの info から",
      fp_r2c["required"] == REQ8 and fp_r2c["friction"] == REQ8[:5] and not fp_r2c["reasons"]
      and fp_r2c["faces_source"].startswith("prepare_info.mesh"), fp_r2c)
for lab, mpath in (("1 回目", R2C / "metrics_judgment1.json"), ("2 回目", R2C / "metrics.json")):
    hist = json.loads(mpath.read_text())["history"]
    bp, bt, rel = series_ok(hist, fp_r2c["required"], fp_r2c["friction"], lab)
    check(f"(s) R2c {lab} の保存系列 ({len(hist)} 標本): 全標本に要求する 8 部品と摩擦 5 部品がそろい (新しい検査を通る)、"
          f"C_T_friction を部品の twall から再計算した値と一致 (相対差 {rel:.1e})",
          len(hist) == 40 and not bp and not bt and rel < 1e-12, (bp[:3], bt[:3], rel))
# 3D の写し 12 本: prepare_info・bcondConfig が写しに無いので、問題 YAML の境界の構成が R2c と同じことを確かめて同じ要求を当てる
Y08_G3L, Y08_G3A = "problem_3d_prod_3op_wallres_lswx08_tewake_g3_L10.yaml", "problem_3d_prod_3op_wallres_lswx08_tewake_g3_A0.yaml"
Y08_G4L, Y10_G3L = "problem_3d_prod_3op_wallres_lswx08_tewake_g4_L10.yaml", "problem_3d_prod_3op_wallres_lswx10_tewake_g3_L10.yaml"
REC3D = [(NOTES / "2026-10-11-g3-local", Y08_G3L),
         (NOTES / "2026-10-11-g4-local" / "run_1085_tewake_g4_L10_cont20k", Y08_G4L),
         (NOTES / "2026-10-11-g4-local" / "run_1086_tewake_g4_L10_cont40k", Y08_G4L),
         (NOTES / "2026-10-11-c6" / "run_1087_tewake_m6_A0", Y08_G3A), (NOTES / "2026-10-11-c6" / "run_1088_tewake_m6_L10", Y08_G3L),
         (NOTES / "2026-10-11-c6" / "run_1089_tewake_m4_A0", Y08_G3A), (NOTES / "2026-10-11-c6" / "run_1090_tewake_m4_L10", Y08_G3L),
         (NOTES / "2026-10-11-c6" / "run_1091_tewake_m6_L10_cont40k", Y08_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1096_tewake_lsw10_m10_cont20k", Y10_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1097_tewake_lsw10_m6_cont20k", Y10_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1098_tewake_lsw10_m4_cont20k", Y10_G3L),
         (NOTES / "2026-10-11-lsw10" / "run_1099_tewake_lsw10_m6_cont40k", Y10_G3L)]


def boundary_layout(yml):
    raw_ = yaml.safe_load(open(C46 / yml))
    m3, m, ev = raw_.get("mesh3d") or {}, raw_.get("mesh") or {}, raw_.get("evaluate") or {}
    return (int(m3.get("nz_out", 17)) > 0, bool(m3.get("vehicle_side", True)), int(m3.get("ext_top", m.get("ext_top", 0))) or int(m.get("ext_top", 0)),
            m3.get("W_vehicle"), str(ev.get("vehicle_top_kind", "wall")), ev.get("model"), (raw_.get("spec") or {}).get("wall_thermal", {}).get("mode"))


lay_r2c = boundary_layout(PROD10.name)
n_rec = 0
for d, yml in REC3D:
    lay = boundary_layout(yml)
    hist = json.loads((d / "metrics.json").read_text())["history"]
    bp, bt, rel = series_ok(hist, REQ8, REQ8[:5], d.name)
    n_rec += len(hist)
    check(f"(s) 3D の写し {d.name} ({len(hist)} 標本): 境界の構成が R2c と同じ ({yml})、全標本に 8 部品・摩擦 5 部品、"
          f"C_T_friction の再計算と一致 (相対差 {rel:.1e})", lay == lay_r2c and not bp and not bt and rel < 1e-12, (lay, lay_r2c, bp[:2], bt[:2], rel))

print(f"\n{n_checks - fails}/{n_checks} PASS" + ("" if not fails else f"  ({fails} FAIL)"))
print(f"(作業ディレクトリ {tmp})")
sys.exit(1 if fails else 0)
