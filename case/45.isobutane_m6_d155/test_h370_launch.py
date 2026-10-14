#!/usr/bin/env python3
"""V-h370 の投入の台本 (h370_prep.sh・h370_run.sh・h370_run.py・h370_judge.py; plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36f・
§6 V-h370) の手元の試験。forge・変換器・AWS は使わない (場は合成、設定は本物の生成器で作る)。

  A. 登録の値の写し: h370_run.py・h370_judge.py の定数が h370_problems.py・問題 YAML (停止規則・B0・Euler の本段・格子・H)・
     V-tw370 (tw370_run.py のバイナリ・環境変数・numpy) と同じ。台本の構文 (bash -n)、h370_run.sh の環境変数が tw370_run.sh と同じ。
  B. 壁の出力の検査 (wall_output_check) と NS の段の設定の検査 (ns_config_checks) を、問題 YAML から本物の生成器 (_bcond・_ns_configs) で
     作った bcond・設定に当てる: 合格 (壁 wall_isothermal・outputHDFflg 1・Ts 370、本段の output.level 1、wallTreatmentSST 0、間隔 2,500)。
     負例 (outputHDFflg 0・output.level 0・wallTreatmentSST 1・prepare_info の sha の不一致) で不合格。
  C. Euler の段の設定の導出 (euler_stage_configs) が runner_axismach.run_staged が実際に起動する設定 (soft・mid・本段) と sha で一致。
     Observer が導いた段の外の設定を拒む (forge を起動しない)。
  D. Euler の予備の資格 (h370_judge.py --euler-prequalify、§5.1 #36g (1)) を合成の Euler の run で: 合格 → 固定 (res_60000・sha)・σ_ref の
     代わりの参照 (末尾の窓の代表 40,000・45,000・50,000・55,000 と sha)・末尾の窓の全出力・スクリーニングの出力、2 回目は同じ固定を使う。
     負例: 出口 M・流量のドリフト (DRIFTING)、場の NaN、残差の停滞 (check_convergence NOT CONVERGED)、既にある違う固定、run の欠け、
     代わりの参照の出力が消された (判定不能)。
  E. 抽出の参照の A/B の分類 (classify_ab)・参照の最終の資格 (reference_final)・1 出力の比較 (ab_pair) の単体の試験 (結果 A・B・判定不能の
     境界、§5.1 #36g (3)〜(5))。codex の判定の基準の諮問の合成の反例 3 つ (中央 1 列の欠測・停止の 1 行だけ・途中の 1 出力だけの参照の感度
     0.030 mm) が「結果 A」にならず理由つきで拒まれる。
  F. NS の判定 (h370_judge.py) を合成の NS・Euler の run (本物の抽出 stop_rule.delta_snapshot・extract_and_merge、関門と次の壁は偽物) で:
     結果 A (腕 A = 範囲の短い同じ参照: 比べる列は座標と範囲から先に固定、範囲の外の列は除く; σ_ref は判定窓の全出力 × 代わりの参照 4 つ)、
     結果 B (腕 A の参照の形を変える)、判定不能 (σ_ref が大きい → 参照の最終の資格なし → 更新に進めない; 代わりの参照の sha が違う →
     σ_ref を測れない)、行が欠けても res があれば抽出し直して補う、反例 3 つの組み込みの版、比べる列の固定に出口の幾何の許容が効く、
     関門で止まる、本段の直前の照合の記録が無い → 判定不能。
  G. h370_run.py の段取り (--dry-run・本番) を偽の準備 (本物の生成器の設定・偽の格子) で: 投入前の検査の合格・不成立・壁の出力が無い
     (終了コード 4)・既存の run で止まる・dry-run と不成立で forge を回さない・資格未達で NS に進まない・環境変数と固定の欠けで止まる・
     h370_problems.py check を作業ツリーの版で回す。
usage: /home/sano/work/forge/design/.venv-opt/bin/python case/45.isobutane_m6_d155/test_h370_launch.py   (所要 約 3〜5 分; 終了コード 0 = ALL PASS)
"""
from __future__ import annotations

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

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[1]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(CASE))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")

import h370_judge as HJ  # noqa: E402
import h370_run as HR  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
from forge_design.metrics import deltastar as DS  # noqa: E402

TMP = Path(tempfile.mkdtemp(prefix="h370_launch_"))
FAIL = 0
TOOLS = ROOT / "solver_density_cuda" / "tools"


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def raises(fn, exc=Exception):
    try:
        fn()
    except exc as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"
    return None


class Patch:
    """属性の差し替え (with で戻す)。"""

    def __init__(self, *items):
        self.items, self.old = items, []

    def __enter__(self):
        for obj, name, val in self.items:
            self.old.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)
        return self

    def __exit__(self, *a):
        for obj, name, val in reversed(self.old):
            setattr(obj, name, val)


# ====================================================================================== A. 登録の値の写し
def test_constants():
    print("== A. 登録の値の写し・台本の構文", flush=True)
    import h370_problems as HP
    import tw370_run as TR
    check("δ_r の CSV の sha256・停止規則の閾値と断面・Euler の上限と throat_refine・壁温が h370_problems.py と同じ",
          HR.CSV_SHA == HP.CSV_SHA and HR.THRESHOLDS == HP.STOP_RULE_THRESHOLDS and HR.X_FRACS == HP.STOP_RULE_X_FRACS
          and HR.EU_STEPS == HP.EU_NSTEP and HR.EU_MESH["throat_refine"] == HP.EU_THROAT_REFINE and HR.TW == HP.TW_NEW
          and HR.NS_STAGE_WANT["main"]["nStepOuter"] == 200000)
    pn, pe = RA.load_problem(CASE / HR.NS_PROBLEM), RA.load_problem(CASE / HR.EU_PROBLEM)
    ml = RA._main_line_spec(pn)
    ev = pn.evaluate
    check("NS の問題 YAML: 停止規則 (上限・保持・閾値)・出力の間隔・B0・格子 (Gref 4719 × 160)・H の上流・壁温が台本の定数と同じ",
          ev["stop_rule"] is True and ev["stop_rule_cap"] == HR.NS_CAP and ev["stop_rule_hold_steps"] == HR.NS_HOLD
          and ev["outStepInterval"] == HR.NS_OUT and ev["nStepOuter"] == HR.NS_CAP and ev["stop_rule_thresholds"] == HR.THRESHOLDS
          and {k: ml[k] for k in HR.MAIN_LINE} == HR.MAIN_LINE and (pn.mesh["ni"], pn.mesh["nj"]) == (HR.NS_MESH["ni"], HR.NS_MESH["nj"])
          and HR.phys_up_ok({"physical_upstream": {"requested": RA._physical_upstream(pn.geometry)}})
          and pn.wall_thermal == {"mode": "isothermal", "Tw": HR.TW},
          str({k: ml[k] for k in HR.MAIN_LINE}))
    me = RA.mesh_euler_block(pe)
    check("Euler の問題 YAML: 本段 60,000・出力 1,000・cfl_main 2.0・mid_stage・格子 2000 × 97・throat_refine 5・H の上流",
          (pe.evaluate["nStepOuter"], pe.evaluate["outStepInterval"], pe.evaluate["cfl_main"], pe.evaluate["mid_stage"])
          == (HR.EU_STEPS, HR.EU_OUT, HR.EU_CFL, True) and (me["ni"], me["nj"], me["throat_refine"]) == tuple(HR.EU_MESH.values())
          and HR.phys_up_ok({"physical_upstream": {"requested": RA._physical_upstream(pe.geometry)}})
          and pe.geometry["initial_line_run"] == HR.CFD_PIN_RUN)
    check("physical_upstream の照合: 解決済みの値 (x_in・mu・frozen_design 付き) の 3 値で比べる; 違う半径は不合格",
          not HR.phys_up_ok({"physical_upstream": {"requested": {**RA._physical_upstream(pe.geometry), "r_inlet": 6.485}}})
          and not HR.phys_up_ok({}))
    check("バイナリ (forge 7e95e5bbd8f50c9e・変換器 711f4438)・FORGE_* の環境変数・numpy 2.5.2 が V-tw370 (tw370_run.py) と同じ",
          (HR.BIN, HR.BIN_SHA, HR.CONV_SHA, HR.ENV_WANT, HR.NUMPY_WANT) == (TR.BIN, TR.BIN_SHA, TR.CONV_SHA, TR.ENV_WANT, TR.NUMPY_WANT))
    check("h370_judge.py の run の名前・固定と資格の記録の場所・バイナリ・A/B の許容 (関門の 2 窓の予算 0.01 mm) が h370_run.py と同じ",
          (HJ.EU_RUN, HJ.NS_RUN, HJ.NS_PROBLEM, HJ.OLD_EULER, HJ.PIN, HJ.QUAL, HJ.EU_STEPS, HJ.BIN_SHA)
          == (HR.EU_RUN, HR.NS_RUN, HR.NS_PROBLEM, HR.OLD_EULER, HR.PIN, HR.QUAL, HR.EU_STEPS, HR.BIN_SHA)
          and HJ.TOL_MM == DL.HANDOFF_WINDOW_TOL_MM)
    sh = {n: subprocess.run(["bash", "-n", str(CASE / n)], capture_output=True, text=True) for n in ("h370_prep.sh", "h370_run.sh")}
    check("h370_prep.sh・h370_run.sh の構文 (bash -n)", all(r.returncode == 0 for r in sh.values()), str({k: v.stderr for k, v in sh.items()}))
    env = lambda p: [ln.strip() for ln in (CASE / p).read_text().splitlines() if ln.strip().startswith(("export FORGE", "for v in $(compgen"))]  # noqa: E731
    check("h370_run.sh の FORGE_* の設定 (export・unset) が tw370_run.sh と同じ行", env("h370_run.sh") == env("tw370_run.sh") and len(env("h370_run.sh")) == 4,
          str(env("h370_run.sh")))


# ====================================================================================== B. 壁の出力・NS の段の設定
def gen_ns_configs(rd: Path):
    """問題 YAML から prepare_ns と同じ生成器で bcond・前の段の元の設定・本段の設定を作る (格子は偽物)。"""
    p = RA.load_problem(CASE / HR.NS_PROBLEM)
    ml = RA._main_line_spec(p)
    n, out_int, _ = RA._ns_steps(p, None)
    rd.mkdir(parents=True)
    (rd / "bcondConfig.yaml").write_text(RA._bcond_with_species(RA._bcond(p, euler=False, profile_ok=True), RA._tp_species_Y(p)))
    cfg_ns, cfg_main = RA._ns_configs(p, rd, n, out_int, None, None, ml)
    (rd / "solverConfig.yaml").write_text(cfg_main)
    (rd / RA.PRESTAGE_CONFIG).write_text(cfg_ns)
    with h5py.File(rd / "nozzle.h5", "w") as f:              # 偽の格子 (mesh_identity が読む最小の中身)
        f["MESH/COORD"] = np.arange(12.0)
        f["VALUE/wall_dist"] = np.ones(4)
        f["VALUE/ro"] = np.ones(4)
    (rd / "nozzle.msh").write_text("$MeshFormat 偽物\n")
    (rd / "probe.yaml").write_text("")
    info = {"ns_stage_configs": RA.ns_stage_sha256(cfg_main, cfg_ns), "main_line": {**ml, "source": "evaluate.main_line"},
            "delta_r_csv": str(HR.HERE / HR.CSV), "physical_upstream": {"requested": RA._physical_upstream(p.geometry), "wall": {}},
            "wall_thermal": dict(p.wall_thermal), "mesh": {"ni": 4719, "nj": 160, "hashes": {"generator_sha256": "x" * 64}},
            "initial_line": {"match": {"checked": ["R", "L_U"], "unchecked": []}}}
    (rd / "prepare_info.json").write_text(json.dumps(info))
    return p, info


def test_wall_output_and_configs():
    print("== B. 壁の出力の検査・NS の段の設定の検査 (本物の生成器の設定)", flush=True)
    rd = TMP / "cfg_ok"
    gen_ns_configs(rd)
    w = HR.wall_output_check(rd)
    check(f"壁の出力: 合格 ({w.get('wall_dump_name')}・output.level {w.get('output_level')}・間隔 {w.get('outStepInterval')}・wallTreatmentSST "
          f"{w.get('wallTreatmentSST')})", w["ok"] and w.get("wall_dump_name") == "res_wall_3_<step>.h5" and w.get("output_level") == 1
          and w.get("outStepInterval") == 2500 and w.get("wallTreatmentSST") == 0, str(w.get("failures")))
    c = HR.ns_config_checks(RA, rd)
    check("NS の段の設定: 合格 (S1 0.5・S2 1.0 は 1 次・3,000 step・B0 のキーなし、本段 cfl 4・200,000・2,500・B0、3 段の implicitRelax 0.7、"
          "prepare_info の sha と一致)", c["ok"] and c["stages"]["main"]["lineDtDirectionalCap"] == 50.0 and c["stages"]["S1_soft"]["lineImplicit"] is None
          and all(c["prepare_info_match"].values()), str(c.get("failures")))
    for tag, fn, want in (("outputHDFflg 0", lambda d: (d / "bcondConfig.yaml").write_text((d / "bcondConfig.yaml").read_text().replace(
                              "kind: wall_isothermal,  outputHDFflg: 1", "kind: wall_isothermal,  outputHDFflg: 0")), "outputHDFflg"),
                          ("output.level 0", lambda d: (d / "solverConfig.yaml").write_text((d / "solverConfig.yaml").read_text().replace(
                              "output: {level: 1,", "output: {level: 0,")), "output.level"),
                          ("wallTreatmentSST 1", lambda d: (d / "solverConfig.yaml").write_text((d / "solverConfig.yaml").read_text().replace(
                              "wallTreatmentSST: 0", "wallTreatmentSST: 1")), "wallTreatmentSST")):
        d = TMP / f"cfg_bad_{tag.split()[0]}"
        shutil.copytree(rd, d)
        before = (d / "bcondConfig.yaml").read_text() + (d / "solverConfig.yaml").read_text()
        fn(d)
        changed = before != (d / "bcondConfig.yaml").read_text() + (d / "solverConfig.yaml").read_text()
        r = HR.wall_output_check(d)
        check(f"壁の出力: 負例 {tag} → 不合格", changed and not r["ok"] and any(want in f for f in r["failures"]), str(r.get("failures")))
    d = TMP / "cfg_bad_sha"
    shutil.copytree(rd, d)
    info = json.loads((d / "prepare_info.json").read_text())
    info["ns_stage_configs"]["main"] = "0" * 64
    (d / "prepare_info.json").write_text(json.dumps(info))
    r = HR.ns_config_checks(RA, d)
    check("NS の段の設定: prepare_info の本段の sha が違う → 不合格", not r["ok"] and not r["prepare_info_match"]["main"])
    r = HR.ns_prepare_info_checks(json.loads((rd / "prepare_info.json").read_text()), None, prelaunch=True)
    check("NS の prepare_info (投入前の準備): 合格 (δ_r の CSV・H・壁温・B0・格子・CFD ピンの照合)", r["ok"], str(r["failures"]))
    bad = {**json.loads((rd / "prepare_info.json").read_text()), "wall_thermal": {"mode": "isothermal", "Tw": 300.0}}
    r = HR.ns_prepare_info_checks(bad, None, prelaunch=True)
    check("NS の prepare_info: 壁温 300 K → 不合格", not r["ok"] and any("370" in f for f in r["failures"]))


# ====================================================================================== C. Euler の段の設定
def test_euler_stage_configs():
    print("== C. Euler の段の設定の導出 = run_staged が起動する設定", flush=True)
    pe = RA.load_problem(CASE / HR.EU_PROBLEM)
    rd = TMP / "eu_cfg"
    rd.mkdir()
    cfg = RA._apply_gas_to_config(RA._config_euler_node(pe, HR.EU_STEPS, HR.EU_OUT, 4.0, 1), pe, rd)
    (rd / "solverConfig.yaml").write_text(cfg)
    (rd / "bcondConfig.yaml").write_text(RA._bcond_with_species(RA._bcond(pe, euler=True), RA._tp_species_Y(pe)))
    want = HR.euler_stage_configs(RA, cfg, HR.EU_CFL, True)
    got = []

    def fake_forge(run_dir):
        t = (Path(run_dir) / "solverConfig.yaml").read_text()
        got.append(HR.sha256_text(t))
        (Path(run_dir) / f"res_{RA._cfg_value(t, 'nStepOuter')}.h5").write_bytes(b"")
        return 0
    stub_res = lambda *a, **k: {"ok": True, "failures": []}  # noqa: E731
    with Patch((RA, "run_forge", fake_forge), (RA, "stage_gate", lambda *a: []), (RA, "_restart_same_mesh", lambda *a: None), (HR, "resources", stub_res)):
        tags = {HR.sha256_text(c): t for t, c in want.items()}
        obs = HR.Observer(RA, rd, tags, "main")
        obs.install()
        try:
            rc = RA.run_staged(rd, cfl_main=HR.EU_CFL, mid_stage=True, stages="full")
        finally:
            obs.uninstall()
        check("run_staged が起動した soft・mid・本段の設定の sha = euler_stage_configs の導出 (Observer は 3 段とも通す)",
              rc == 0 and got == [HR.sha256_text(want[t]) for t in ("soft", "mid", "main")], f"{len(got)} 段")
        ev = [json.loads(ln) for ln in (rd / "h370_stages.jsonl").read_text().splitlines()]
        check("Observer の記録: 段の start・end と本段の前の資源", [e.get("tag") for e in ev if e["event"] == "start"] == ["soft", "mid", "main"]
              and any(e["event"] == "resources" for e in ev))
        rd2 = TMP / "eu_cfg2"
        shutil.copytree(rd, rd2, ignore=shutil.ignore_patterns("res_*", "h370_stages.jsonl", "stage_manifest.json", "residual_history_*"))
        (rd2 / "solverConfig.yaml").write_text(cfg)
        got.clear()
        obs = HR.Observer(RA, rd2, {k: v for k, v in tags.items() if v != "mid"}, "main")
        obs.install()
        try:
            e = raises(lambda: RA.run_staged(rd2, cfl_main=HR.EU_CFL, mid_stage=True, stages="full"), RuntimeError)
        finally:
            obs.uninstall()
        check(f"Observer: 導いた段の外の設定 (mid を外した表) → forge を起動しない ({(e or '')[:40]}…)", e is not None and len(got) == 1)
    st = {t: {k: RA._cfg_value(c, k) for k in ("convMethod", "cfl", "nStepOuter", "outStepInterval")} for t, c in want.items()}
    check("Euler の段: soft は 1 次・cfl 0.5・3,000、mid は cfl 1.0・3,000、本段は cfl 2.0・60,000・出力 1,000",
          st["soft"]["convMethod"] == 0 and st["soft"]["cfl"] == 0.5 and st["soft"]["nStepOuter"] == 3000 and st["mid"]["cfl"] == 1.0
          and st["mid"]["nStepOuter"] == 3000 and (st["main"]["cfl"], st["main"]["nStepOuter"], st["main"]["outStepInterval"]) == (2.0, 60000, 1000),
          str(st))


# ====================================================================================== D. Euler の資格
EU_NI, EU_NJ, EU_S = 80, 9, 0.0767


def make_eu_run(case: Path, drift: float = 0.0, nan_step=None, flat_resid: bool = False, steps=None):
    """合成の Euler の run (構造格子 80 × 9、本段の出力 1,000〜60,000、段の残差の履歴と stage_manifest、RUN_RC 0、投入前の検査の合格)。"""
    from stage_manifest import StageManifest
    run = case / HJ.EU_RUN
    run.mkdir(parents=True)
    x = np.linspace(-12.5, 95.3, EU_NI)
    rw = np.where(x < 0, 1.0 + 0.017 * x * x, 1.0 + 0.09 * x)
    eta = np.linspace(0.0, 1.0, EU_NJ)
    X, R = np.repeat(x[:, None], EU_NJ, 1), rw[:, None] * eta[None, :]
    with h5py.File(run / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = np.c_[X.ravel() * EU_S, R.ravel() * EU_S, np.zeros(X.size)].ravel()
    M0 = np.where(x[:, None] < 0, 0.2 + 0.8 * np.exp(x[:, None] / 3.0), 1.0 + 5.0 * (1.0 - np.exp(-np.maximum(x[:, None], 0.0) / 20.0))) \
        * (1.0 - 0.002 * eta[None, :] ** 2)
    for st in (steps or range(1000, 60001, 1000)):
        k = 1.0 + 1e-3 * np.exp(-st / 5000.0) + drift * st / 60000.0
        M = M0 * k
        T = 1600.0 / (1.0 + 0.135 * M ** 2)
        a = np.sqrt(1.27 * 287.0 * T)
        P = 5.5e6 * (T / 1600.0) ** (1.27 / 0.27)
        ro = P / (287.0 * T)
        with h5py.File(run / f"res_{st}.h5", "w") as f:
            for n_, v in (("ro", ro), ("Ux", M * a), ("Uy", 0.0 * M), ("Uz", 0.0 * M), ("P", P), ("T", T), ("sonic", a), ("roe", P / 0.27)):
                v = np.asarray(v, dtype=float).ravel().copy()
                if nan_step == st and n_ == "Ux":
                    v[5] = np.nan
                f[f"/VALUE/{n_}"] = v
    (run / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": EU_NI, "nj": EU_NJ}, "scale_m": EU_S, "x_E": 40.0, "Md": 6.0}))
    cfgs = {"soft": "space: {convMethod: 0, limiter: 2}\ntime:\n  deltaT: {cfl: 0.5, cfl_pseudo: 0.5}\n",
            "mid": "space: {convMethod: 1, limiter: 2}\ntime:\n  deltaT: {cfl: 1.0, cfl_pseudo: 1.0}\n",
            "main": "space: {convMethod: 1, limiter: 2}\ntime:\n  deltaT: {cfl: 2.0, cfl_pseudo: 2.0}\n"}
    sm = StageManifest(run)
    for tag, n in (("soft", 3000), ("mid", 3000), ("main", 60000)):
        s = np.arange(1, n + 1, 10)
        r = np.full(len(s), 1e-3) if (flat_resid and tag == "main") else 1e-2 * 10.0 ** (-5.0 * s / n)
        with open(run / f"residual_history_{tag}.csv", "w") as fh:
            fh.write("step,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n")
            for si, ri in zip(s, r):
                fh.write(f"{si},{ri:.6e},{2 * ri:.6e},{0.5 * ri:.6e},0,{3 * ri:.6e}\n")
        sm.add(tag, cfgs[tag], "wall: {kind: slip}\n", history=f"residual_history_{tag}.csv")
    sm.write()
    (run / "RUN_RC").write_text("0\n")
    (run / "H370_PRECHECK.json").write_text(json.dumps({"ok": True}))
    return run


def qualify(case: Path) -> int:
    return HJ.main(["--euler-prequalify", "--case-dir", str(case), "--design", str(ROOT / "design"), "--tools", str(TOOLS)])


def test_euler_qualify():
    print("== D. Euler の予備の資格 (合成の run)", flush=True)
    case = TMP / "q_ok"
    run = make_eu_run(case)
    rc = qualify(case)
    q, pin = json.loads((case / HJ.QUAL).read_text()), HJ.jload(case / HJ.PIN)
    v = {c: d["verdict"] for c, d in q["quasisteady"]["verdicts"].items()}
    check(f"合格: 予備の資格あり (rc {rc})、資格の量 {len(v)} 列すべて STEADY、check_convergence --segment PASS、記録の名前は予備の資格 (prequalified)",
          rc == 0 and q["prequalified"] and "qualified" not in q and set(v) == set(HJ.QUAL_COLS) and set(v.values()) == {"STEADY"}
          and q["convergence"]["rc"] == 0 and HJ.QUAL.endswith("h370_euler_prequalify.json"), str(q["failures"] + q["undecidable"]))
    check("固定: 本段の最後の出力 res_60000.h5 とその sha256 (h370_euler_ref_pin.json)",
          pin is not None and pin["res"] == "res_60000.h5" and pin["sha256"] == HJ.sha256_file(run / "res_60000.h5") and pin["run"] == HJ.EU_RUN)
    sr = q["sigma_ref_snapshots"]
    check("σ_ref の代わりの参照: 末尾の窓 [40,000, 60,000] の 5,000 step ごとの代表 40,000・45,000・50,000・55,000 (固定を除く)、sha つき; "
          "末尾の窓の全 21 出力の sha も記録",
          [o["step"] for o in sr] == [40000, 45000, 50000, 55000] and all(o["sha256"] == HJ.sha256_file(run / o["res"]) for o in sr)
          and [o["step"] for o in q["tail_outputs"]] == list(range(40000, 60001, 1000))
          and q["sigma_ref_selection"]["tail_window"] == [40000, 60000] and "340 抽出" in q["sigma_ref_selection"]["why"], str([o["step"] for o in sr]))
    td = q.get("tail_field_dev") or {}
    check("代表の選び方の裏付け (記録だけ): 末尾の窓の 20 出力の場の差と、代表の中・外の最大", len(td.get("per_output", [])) == 20
          and td.get("max_rel_representatives") is not None and td.get("max_rel_others") is not None, str(td)[:200])
    check("スクリーニングの出力: check_quasisteady の末尾の窓 (既定 0.4 → 24 出力) の始まり 37,000 と中ほど 49,000、sha つき",
          [s_["step"] for s_ in q["screening_snapshots"]] == [37000, 49000] and q["quasisteady_tail_window"] == [37000, 60000]
          and all(len(s_["sha256"]) == 64 for s_ in q["screening_snapshots"]))
    check("記録だけの量 (nozzle_report.metrics の η 0.1 の線の wave・overshoot・slope・range) の VERDICT も残る",
          set(q["quasisteady_record"]["record_cols (既定の閾値)"]["verdicts"]) == set(HJ.REC_COLS))
    rc2 = qualify(case)
    check("2 回目: 同じ固定をそのまま使う (rc 0、固定のファイルは変えない)", rc2 == 0 and HJ.jload(case / HJ.PIN) == pin
          and json.loads((case / HJ.QUAL).read_text()).get("pin_note") == "既にある同じ固定をそのまま使う")
    for name, kw, want_rc, want in (("drift", {"drift": 0.2}, 1, "STEADY"), ("nan", {"nan_step": 45000}, 1, "非有限"),
                                    ("flat", {"flat_resid": True}, 1, "check_convergence")):
        c = TMP / f"q_{name}"
        make_eu_run(c, **kw)
        rc = qualify(c)
        q = json.loads((c / HJ.QUAL).read_text())
        check(f"負例 {name}: 未達 (rc {rc})・固定を書かない", rc == want_rc and not q["prequalified"] and not (c / HJ.PIN).exists()
              and any(want in f for f in q["failures"]), str(q["failures"])[:300])
    c = TMP / "q_conflict"
    shutil.copytree(case, c)
    p = json.loads((c / HJ.PIN).read_text())
    p["sha256"] = "0" * 64
    (c / HJ.PIN).write_text(json.dumps(p))
    rc = qualify(c)
    check("既にある違う固定 → 判定不能 (rc 2)・固定を消さない", rc == 2 and json.loads((c / HJ.PIN).read_text())["sha256"] == "0" * 64)
    c = TMP / "q_missing"
    c.mkdir()
    check("run が無い → 判定不能 (rc 2)", qualify(c) == 2)
    c = TMP / "q_short"
    make_eu_run(c, steps=range(1000, 50001, 1000))
    rc = qualify(c)
    check("本段が上限 60,000 まで回っていない → 未達 (rc 1)", rc == 1 and any("上限" in f for f in json.loads((c / HJ.QUAL).read_text())["failures"]))
    c = TMP / "q_gap"
    make_eu_run(c, steps=[s_ for s_ in range(1000, 60001, 1000) if s_ != 45000])
    rc = qualify(c)
    q = json.loads((c / HJ.QUAL).read_text())
    check("σ_ref の代わりの参照の出力 (45,000) が消されている → 判定不能 (rc 2)・固定を書かない", rc == 2 and not (c / HJ.PIN).exists()
          and any("45000" in u for u in q["undecidable"]), str(q["undecidable"])[:200])


# ====================================================================================== E. A/B の分類
def test_classify_ab():
    print("== E. 抽出の参照の A/B の分類・参照の最終の資格・1 出力の比較 (単体)", flush=True)
    req2 = [0, 1]
    per = lambda Ds, pf=None: [{"step": i, "D_mm": d, "pf_A": True if pf is None else pf[i], "pf_B": True} for i, d in enumerate(Ds)]  # noqa: E731

    def sig(mm, alts_pf=None, status="ok"):
        """σ_ref の記録の形 (per_step の alts の pf_F は代わりの参照の時刻の腕 B の合否)。"""
        ps = None if alts_pf is None else {k: {"alts": [{"res": f"res_{e}.h5", "pf_F": v} for e, v in enumerate(alts_pf)]} for k in req2}
        return {"status": status, "mm": mm, "per_step": ps, "why": "試験"}
    cases = [("全出力 0.003 mm・σ 0.001 → 結果 A", per([0.003, 0.002]), 0.001, "結果 A", "A"),
             ("全出力 0.02 mm・σ 0.001 → 結果 B", per([0.02, 0.03]), 0.001, "結果 B", "B_D"),
             ("窓の中で 0.01 mm をまたぐ (0.008・0.012) → 判定不能 (時間の変動)", per([0.008, 0.012]), 0.0, "判定不能", "U_time"),
             ("D 0.008・σ 0.005 → 判定不能 (参照の未整定)", per([0.008, 0.008]), 0.005, "判定不能", "U_sigma_straddle"),
             ("D 0.012・σ 0.005 → 判定不能 (参照の未整定)", per([0.012, 0.013]), 0.005, "判定不能", "U_sigma_straddle"),
             ("合否が変わる出力・代わりの参照の全部で違いが保たれる → 結果 B (D が小さくても)",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": True}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              sig(0.0, [True, True]), "結果 B", "B_pf"),
             ("合否が変わる出力・代わりの参照で確かめられない (per_step が無い) → 判定不能",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": True}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              0.0, "判定不能", "U_pf_unconfirmed"),
             ("合否が変わる出力・代わりの参照の 1 つで違いが消える → 判定不能",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": True}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              sig(0.0, [True, False]), "判定不能", "U_pf_unconfirmed"),
             ("両腕とも比べる列の中に欠測 (合否は同じ) → 判定不能 (比較を合格にしない)",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": False}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              0.0, "判定不能", "U_fixed_cols"),
             ("ちょうど 0.01 mm (max + σ) → 結果 A (≤)", per([0.009, 0.009]), 0.001, "結果 A", "A"),
             ("抽出の失敗 → 判定不能", per([0.001]) + [{"step": 1, "error_A": "x"}], 0.0, "判定不能", "U_error"),
             ("σ を測れない → 判定不能", per([0.001, 0.001]), None, "判定不能", "U_sigma_unmeasurable"),
             ("σ が判定不能 (代わりの参照で合否が変わる) → 判定不能", per([0.001, 0.001]), sig(0.001, status="undecidable"), "判定不能",
              "U_sigma_undecidable"),
             ("出力が無い → 判定不能", [], 0.0, "判定不能", "U_empty")]
    for name, p_, sg, want, code in cases:
        c = HJ.classify_ab(p_, sg, req2)
        check(name, c["result"] == want and c["code"] == code, f"{c['result']} {c['code']}: {c['why'][:120]}")
    c = HJ.classify_ab(per([0.001, 0.001]), 0.0, None)
    check("必要な出力の集合が無い → 判定不能", c["result"] == "判定不能" and c["code"] == "U_required")
    # --- codex の判定の基準の諮問の合成の反例 (2026-10-14-h370-judge-criteria-diagnose.md) ---
    req17 = list(range(22500, 62501, 2500))
    check("必要な出力の列挙: 停止 62,500・窓 20,000・間隔 2,500 の 2 窓 → 17 出力 (22,500〜62,500)",
          SR.required_window_steps(62500, 20000, 2500) == req17 and len(req17) == 17)
    c = HJ.classify_ab([{"step": 62500, "D_mm": 0.0, "pf_A": True, "pf_B": True}], 0.0, req17)
    check(f"反例 2 (停止時の 1 行だけ): 必要な 17 出力のうち 16 を比べていない → 判定不能 ({c['why'][:60]}…)", c["result"] == "判定不能"
          and c["code"] == "U_incomplete")
    per17 = [{"step": k, "D_mm": 0.005, "pf_A": True, "pf_B": True} for k in req17]
    ps = {k: {"alts": [{"res": "res_45000.h5", "pf_F": True, "max_mm_F": (0.030 if k == 42500 else 0.0)}]} for k in req17}
    sg3 = {"status": "ok", "mm": max(a["max_mm_F"] for v in ps.values() for a in v["alts"]), "per_step": ps}
    c = HJ.classify_ab(per17, sg3, req17)
    rf = HJ.reference_final({"classify": c, "sigma_ref": sg3}, True)
    check(f"反例 3 (17 出力で D 0.005 mm、停止の出力の感度 0・途中の 1 出力だけ代わりの参照で 0.030 mm): σ_ref = 全出力の最大 0.030 mm → "
          f"判定不能 (参照の未整定)、参照の最終の資格なし ({c['code']})",
          sg3["mm"] == 0.030 and c["result"] == "判定不能" and c["code"] == "U_sigma_straddle" and rf["qualified"] is False)
    # 反例 1 (中央 1 列の欠測): 比べる列 F を先に固定 → 腕 A の F の中の欠測は合否の違い (結果 A にしない)
    x3 = np.array([10.0, 11.0, 12.0])
    F3 = np.ones(3, dtype=bool)
    dB = {"delta_E": [1e-3] * 3, "delta_target": [2e-3] * 3, "held": [0] * 3, "missing": []}
    dA = {"delta_E": [1e-3, float("nan"), 1e-3], "delta_target": [2e-3, float("nan"), 2e-3], "held": [0] * 3, "missing": [1]}
    r1 = HJ.ab_pair(62500, dA, dB, 76.6, "記録", x3, F3)
    check(f"反例 1 (中央 1 列だけ腕 A が欠測): ab_pair は比べる列 3 を保ち、腕 A は合否 否 (欠測 1・x {r1.get('x_bad_A')})・D は出さない",
          r1["n_fixed"] == 3 and r1["pf_A"] is False and r1["pf_B"] is True and r1["D_mm"] is None and r1["check_A"]["n_missing"] == 1
          and r1.get("x_bad_A") == [11.0], str(r1))
    c1 = HJ.classify_ab([r1], {"status": "ok", "mm": 0.0, "per_step": {62500: {"alts": [{"res": "r", "pf_F": True}]}}}, [62500])
    c1b = HJ.classify_ab([r1], 0.0, [62500])
    check(f"反例 1: 代わりの参照でも腕 B が合 → 結果 B (合否の違い) / 代わりの参照で確かめられなければ判定不能 — どちらも結果 A ではない "
          f"({c1['code']}・{c1b['code']})", c1["result"] == "結果 B" and c1["code"] == "B_pf" and c1b["result"] == "判定不能")
    r_both = HJ.ab_pair(62500, dA, dA, 76.6, "記録", x3, F3)
    check("両腕とも同じ列が欠測 → 合否は同じ (否)・D を出さない → 判定不能 (F の中の欠測)", r_both["D_mm"] is None
          and HJ.classify_ab([r_both], 0.0, [62500])["code"] == "U_fixed_cols")
    r_out = HJ.ab_pair(62500, dA, dB, 76.6, "記録", x3, np.array([True, False, True]))
    check("欠測の列が座標で先に決めた比べる列の外 (範囲の外) なら比べない: 合否は同じ・D 0", r_out["pf_A"] and r_out["D_mm"] == 0.0)
    # 参照の最終の資格 (reference_final)
    mk = lambda code: {"classify": {"code": code, "result": "x", "why": "試験"}, "sigma_ref": {"mm": 0.001, "status": "ok"}}  # noqa: E731
    check("参照の最終の資格: A・B・時間の変動の判定不能 → 資格あり / σ を測れない・判定不能・参照の未整定 → 資格なし / 出力の欠け・抽出の失敗・"
          "F の中の欠測 → 確かめられない (None) / --no-ab → 確かめられない",
          all(HJ.reference_final(mk(k), True)["qualified"] is True for k in ("A", "B_D", "B_pf", "U_time"))
          and all(HJ.reference_final(mk(k), True)["qualified"] is False for k in ("U_sigma_unmeasurable", "U_sigma_undecidable", "U_sigma_straddle"))
          and all(HJ.reference_final(mk(k), True)["qualified"] is None for k in ("U_incomplete", "U_error", "U_fixed_cols", "U_pf_unconfirmed"))
          and HJ.reference_final(None, False)["qualified"] is None)


# ====================================================================================== F. NS の判定
NS_S, NI, NJ, X0, X1 = 0.0766, 56, 48, -2.0, 20.0
S_STOP, W_TEST = 62500, 5000


def _rwE(x):
    return 1.0 + 0.002 * x * x


def _dl(x):
    return 0.01 + 0.002 * (x + 2.0)


def _core(x, r, rw):
    return 100.0 * (1.0 + 0.02 * np.tanh(x / 10.0)) * (1.0 - 0.01 * (r / rw) ** 2)


def _radial(n):
    return 1.0 - (1.0 - np.linspace(0.0, 1.0, n)) ** 2


def write_field_run(rd: Path, X, R, Q, info: dict, steps, perturb=None):
    rd.mkdir(parents=True, exist_ok=True)
    ni, nj = X.shape
    with h5py.File(rd / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = np.c_[X.ravel() * NS_S, R.ravel() * NS_S, np.zeros(ni * nj)].ravel()
    for st in steps:
        k = 1.0 if perturb is None else perturb(st)
        with h5py.File(rd / f"res_{st}.h5", "w") as f:
            f["/VALUE/ro"] = np.ones(ni * nj)
            f["/VALUE/Ux"] = (Q * k).ravel()
            f["/VALUE/Uy"] = np.zeros(ni * nj)
            f["/VALUE/sonic"] = np.full(ni * nj, 300.0)
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": ni, "nj": nj}, "scale_m": NS_S, "x_E": 15.0, **info}))


def synth_ns():
    x = np.linspace(X0, X1, NI)
    rw = _rwE(x) + _dl(x)
    R = rw[:, None] * _radial(NJ)[None, :]
    d = rw[:, None] - R
    slope = (0.004 * x + 0.002) / np.sqrt(1.0 + (0.004 * x + 0.002) ** 2)
    beta = np.clip(1.0 - d / (0.1 * rw[:, None]), 0.0, 1.0)
    X = x[:, None] + beta * d * slope[:, None]
    y = rw[:, None] - R
    return X, R, _core(X, R, rw[:, None]) * np.minimum(1.0, np.maximum(y / (5.0 * _dl(X)), 0.0) ** (1.0 / 7.0))


def synth_eu(x1: float, shape: float = 0.0):
    n = int(round(69 * (x1 - X0 + 0.5) / (X1 - X0 + 1.0))) + 1
    x = np.linspace(X0 - 0.5, x1, n)
    rw = _rwE(x)
    R = rw[:, None] * _radial(NJ)[None, :]
    X = np.repeat(x[:, None], NJ, axis=1)
    return X, R, _core(X, R, rw[:, None]) * (1.0 + shape * (R / rw[:, None]) ** 8)


def write_eu_shape(rd: Path, step: int, shape: float):
    X, R, Q = synth_eu(X1 + 0.5, shape)
    with h5py.File(rd / f"res_{step}.h5", "w") as f:
        f["/VALUE/ro"] = np.ones(Q.size)
        f["/VALUE/Ux"] = Q.ravel()
        f["/VALUE/Uy"] = np.zeros(Q.size)
        f["/VALUE/sonic"] = np.full(Q.size, 300.0)


ALT_SHAPES = {40000: 1e-4, 45000: 8e-5, 50000: 5e-5, 55000: 2e-5}     # σ_ref の代わりの参照 (末尾の窓の代表) の形の小さな違い
SCREEN_SHAPES = {37000: 1e-4, 49000: 5e-5}                              # スクリーニングの出力


def make_ns_case(case: Path):
    """合成の NS の判定の case dir: Euler B (run_0950、固定 res_60000、σ_ref の代わりの参照 res_40000〜55000、スクリーニング res_37000・
    res_49000)・腕 A (run_0174、範囲 x ≤ 19 の同じ形)・NS (run_0951、出力 52,500〜62,500、停止の出力 62,500、窓 5,000 の 2 窓)・記録一式。"""
    case.mkdir(parents=True)
    eb = case / HJ.EU_RUN
    X, R, Q = synth_eu(X1 + 0.5)
    write_field_run(eb, X, R, Q, {}, (60000,))
    for st, sh in {**ALT_SHAPES, **SCREEN_SHAPES}.items():
        write_eu_shape(eb, st, sh)
    pin = {"run": HJ.EU_RUN, "res": "res_60000.h5", "sha256": HJ.sha256_file(eb / "res_60000.h5"), "step": 60000}
    (case / "_band_ab/half_inlet").mkdir(parents=True)
    (case / HJ.PIN).write_text(json.dumps(pin))
    snap = lambda steps: [{"res": f"res_{s_}.h5", "step": s_, "sha256": HJ.sha256_file(eb / f"res_{s_}.h5")} for s_ in steps]  # noqa: E731
    qual = {"prequalified": True, "pin": pin, "sigma_ref_snapshots": snap(sorted(ALT_SHAPES)), "screening_snapshots": snap(sorted(SCREEN_SHAPES)),
            "sigma_ref_selection": HJ.sigma_ref_selection(sorted([*ALT_SHAPES, 60000]))}
    (case / HJ.QUAL).write_text(json.dumps(qual))
    ea = case / HJ.OLD_EULER
    X, R, Q = synth_eu(19.0)
    write_field_run(ea, X, R, Q, {}, (12000,))
    spec_b = {"run": str(eb), "res": "res_60000.h5", "sha256": pin["sha256"]}
    ns = case / HJ.NS_RUN
    sr = SR.resolve_settings({"stop_rule": True, "stop_rule_thresholds": {"window": W_TEST, "min_points": 2}}, nsteps_arg=None, n_main=200000,
                             n_main_source="試験", out_int=2500, euler_ref=spec_b)
    steps = list(range(S_STOP - 2 * W_TEST, S_STOP + 1, 2500))
    X, R, Q = synth_ns()
    write_field_run(ns, X, R, Q, {"wall_thermal": {"mode": "isothermal", "Tw": 370.0}, "physical_wall": {"domain": [X0, X1]},
                                  "euler_ref": spec_b, "stop_rule": sr}, steps, perturb=lambda st: 1.0 + 1e-6 * (st - S_STOP) / 2500.0)
    uc = SR.use_columns(ns)
    ext = {k: v for k, v in sr["extract"].items() if k != "source"}
    rows = []
    for st in steps:
        row = {"step": st, "nonfinite_all": 0, "r": 1e-5 * np.sin(st), "exit_core_M": 5.99 + 1e-6 * np.cos(st), "mdot_rt": 10.0,
               SR.theta_key(0.5): 0.01}
        if st != steps[1]:                                      # 1 出力は記録なし (判定が抽出し直す)
            row["delta"] = SR.delta_snapshot(ns, st, spec_b, ext, uc["x"])
        rows.append(row)
    rec = {"version": 1, "settings": {**sr, "x_fracs": [0.5]}, "status": "REACHED", "rows": rows, "stop_step": S_STOP,
           "sections": [{"key": SR.theta_key(0.5), "x_frac": 0.5, "x_requested": 10.0}],
           "delta_series": {"x": uc["x"].tolist(), "use_domain": [0.0, uc["x_w"]], "scale_m": uc["scale_m"], "extract": sr["extract"],
                            "euler_ref": sr["euler_ref"], "span_steps": 2 * W_TEST}}
    SR.save_record(ns, rec)
    SR.write_state(ns, {"status": "REACHED", "stop_step": S_STOP, "forge_rc": 143, "stop_sha256": HJ.sha256_file(ns / f"res_{S_STOP}.h5")})
    (ns / f"res_wall_3_{S_STOP}.h5").write_bytes(b"")
    (ns / "H370_PRECHECK.json").write_text(json.dumps({"ok": True}))
    (ns / "H370_MAIN_START.json").write_text(json.dumps({"ok": True}))
    (ns / "RUN_RC").write_text("0\n")
    (ns / "RUN_PROVENANCE.txt").write_text("forge_sha256: 7e95e5bbd8f50c9e0000\nenv         : FORGE_CUDA_BLOCKSIZE=128\n")
    for t in ("S1_soft", "S2_mid", "main"):
        (ns / f"residual_history_{t}.csv").write_text("step,rms_ro\n1,1e-3\n2,1e-4\n")
    (case / HJ.NS_PROBLEM).write_text("# 試験の偽物 (関門・次の壁は偽物で切り替える)\n")
    return case


def stub_gate(ok=True):
    def g(run, **kw):
        return {"proceed": ok, "blocked_by": [] if ok else ["収束の判定 NOT_CONVERGED (試験の偽物)"],
                "stop_rule": {"status": "REACHED"}, "convergence": {"kind": "PASS" if ok else "NOT_CONVERGED"},
                "checks": {"two_windows": {"pass": ok}, "wall_resolution": {"pass": ok, "kind": "PASS"}}, "extract": {"made": []},
                "override": {"allow_not_reached": False}, "_kw": {k: (str(v) if isinstance(v, Path) else v) for k, v in kw.items()}}
    return g


def stub_next_wall(problem, table):
    a = np.loadtxt(table, delimiter=",", skiprows=1, ndmin=2)
    return {"pass": True, "throat": {"x": 0.0, "r": 1.0 + float(np.interp(0.0, a[:, 0], a[:, 1])), "kappa": 1.0}}


def judge(case: Path, extra=()) -> tuple:
    rc = HJ.main(["--case-dir", str(case), "--design", str(ROOT / "design"), "--tools", str(TOOLS), *extra])
    return rc, HJ.jload(case / HJ.OUT)


def snap_hook(fn):
    """本物の stop_rule.delta_snapshot の結果を fn(run, step, spec, d) で書き換える試験の差し替え (反例を組み込みで作る)。"""
    def g(run, step, spec, ext, x_use):
        return fn(run, step, spec, SR.delta_snapshot(run, step, spec, ext, x_use))
    return g


def mid_col(R):
    """比べる列 F の中ほどの列の番号 (delta_series.x の番号)。"""
    x = np.asarray(json.loads((R["_case"] / HJ.NS_RUN / SR.RECORD).read_text())["delta_series"]["x"])
    lo, hi = R["ab"]["columns"]["x_range"]
    return int(np.argmin(np.abs(x - 0.5 * (lo + hi))))


def test_ns_judge():
    print("== F. NS の判定 (合成の NS・Euler の run、本物の抽出、関門と次の壁は偽物)", flush=True)
    base = make_ns_case(TMP / "ns_base")
    hooks = {"handoff_gate": stub_gate(True), "next_wall_check": stub_next_wall}
    req = SR.required_window_steps(S_STOP, W_TEST, 2500)
    with Patch((HJ, "HOOKS", dict(hooks))):
        case = TMP / "ns_A"
        shutil.copytree(base, case)
        rc, R = judge(case)
        R["_case"] = case
        ab = R["ab"]
        srcs = {r.get("B_source", "")[:2] for r in ab["per_output"]}
        sg = ab["sigma_ref"]
        check(f"結果 A: 腕 A = 範囲の短い同じ参照 (D max {ab.get('D_mm', {}).get('max', float('nan')):.4f} mm、σ_ref {sg['mm']:.2e} mm "
              f"(使用域 {sg['mm_use']:.2e})、{ab['classify']['code']})", rc == 0 and ab["result"] == "結果 A" and ab["D_mm"]["max"] < HJ.TOL_MM
              and sg["status"] == "ok", ab["why"])
        cols = ab["columns"]
        check(f"比べる列は抽出の前に座標と範囲から固定: 使用域 {cols['n_use']} 列のうち F {cols['n']} 列 (腕 A の範囲の外 {cols['excluded']['A']['n']} 列 "
              f"x {cols['excluded']['A']['x_range']} は比べない)", cols["excluded"]["A"]["n"] == 3 and cols["excluded"]["A"]["x_range"][0] > 19.0
              and cols["excluded"]["B"]["n"] == 0 and cols["n"] == cols["n_use"] - 3 and cols["x_range"][1] <= 19.0 and cols["edge_tol"])
        check(f"判定窓の完全性: 入るべき出力 {ab['completeness']['n_required']} (出力の間隔 2,500 と 2 窓から列挙) をすべて比べた; 腕 B の値は記録 "
              f"(sha が同じ) と、記録の無い出力は抽出し直し", ab["completeness"]["required"] == req and [r["step"] for r in ab["per_output"]] == req
              and srcs == {"記録", "抽出"} and ab["completeness"]["reextracted_B"] == [req[1]], str(srcs))
        check(f"σ_ref: 判定窓の全 {len(req)} 出力 × 代わりの参照 4 つ (40,000・45,000・50,000・55,000) で抽出し直した (合否の変化なし)",
              sorted(int(k) for k in sg["per_step"]) == req and all(len(v["alts"]) == 4 and all(a_["pf_F"] and a_["pf_use"] and not a_["change"] for a_ in v["alts"])
                                                    for v in sg["per_step"].values()) and sg["alts"] == [f"res_{k}.h5" for k in sorted(ALT_SHAPES)])
        check("スクリーニング (停止の出力 × 37,000・49,000) は記録だけ", ab["screening"]["label"].startswith("スクリーニング")
              and len(ab["screening"]["per_snapshot"]) == 2)
        rf = R["reference_final"]
        check(f"参照の最終の資格: あり ({rf['why'][:50]}…)", rf["qualified"] is True and rf["ab_code"] == "A")
        check(f"判定の全体: 「更新に進める条件をすべて満たす / 抽出の参照の A/B: 結果 A」({R['verdict'][:60]}…)",
              R["verdict"].startswith("更新に進める条件をすべて満たす") and "結果 A" in R["verdict"] and R["update_conditions_met"])
        st = ab["at_stop"]
        check("停止の出力の両腕の全体の抽出: 共通の列・帯の位置 (y_b・y*)・感度・held・δ_target・未緩和の固定点の残差・次の壁の真のスロート",
              st["compare"]["n_common_use"] > 0 and "band_y_b_max_diff_mm" in st["compare"] and "band_sensitivity_rel" in st["compare"]
              and "B_fixed_point_residual_mm" in st["compare"] and "throat_diff" in st["next_wall"]
              and all((case / HJ.NS_RUN / HJ.AB_DIR / d / "delta_r_next.csv").is_file() for d in ("A", "B", "A_on_common")), str(st.get("compare"))[:300])
        q = R.get("quasisteady") or {}
        check("check_quasisteady の記録 (2 窓の r・θ_r・出口 M・流量・δ_E・δ_target)", q.get("n") == 5 and "dE_40" in q["default"]["verdicts"]
              and "exit_core_M" in q["default"]["verdicts"], str(q.get("default", {}).get("missing")))
        check("関門は固定の Euler の参照・ω 0.5・上書きなしで呼ばれる",
              R["gate"]["proceed"] and R["stop"]["stop_sha256_match"] and R["stop"]["wall_dump_at_stop"] == [f"res_wall_3_{S_STOP}.h5"])
        j_mid = mid_col(R)
        mmpr = NS_S * 1e3
        # 結果 B: 腕 A の参照の形を変える
        case = TMP / "ns_B"
        shutil.copytree(base, case)
        shutil.rmtree(case / HJ.OLD_EULER)
        X, Rr, Q = synth_eu(19.0, 0.01)
        write_field_run(case / HJ.OLD_EULER, X, Rr, Q, {}, (12000,))
        rc, R = judge(case)
        check(f"結果 B: 腕 A の参照の形を変える (D min {R['ab'].get('D_mm', {}).get('min', float('nan')):.4f} mm > 0.01 mm)",
              rc == 0 and R["ab"]["result"] == "結果 B" and R["ab"]["classify"]["code"] == "B_D" and R["reference_final"]["qualified"] is True,
              R["ab"]["why"])
        # 判定不能: 代わりの参照 (45,000) の形が大きく違う → σ_ref が大きい → 参照の最終の資格なし → 更新に進めない
        case = TMP / "ns_sigma"
        shutil.copytree(base, case)
        eb = case / HJ.EU_RUN
        write_eu_shape(eb, 45000, 0.003)
        q = json.loads((case / HJ.QUAL).read_text())
        q["sigma_ref_snapshots"] = [{**s_, "sha256": HJ.sha256_file(eb / s_["res"])} for s_ in q["sigma_ref_snapshots"]]
        (case / HJ.QUAL).write_text(json.dumps(q))
        rc, R = judge(case)
        check(f"判定不能: σ_ref {R['ab']['sigma_ref']['mm']:.4f} mm で 0.01 mm をまたぐ (参照の未整定) → 参照の最終の資格なし → "
              "関門が進めても更新に進める条件を満たさない",
              R["ab"]["result"] == "判定不能" and R["ab"]["classify"]["code"] == "U_sigma_straddle" and R["reference_final"]["qualified"] is False
              and R["gate"]["proceed"] and not R["update_conditions_met"] and R["verdict"].startswith("更新に進める条件を満たさない")
              and "参照の最終の資格" in R["verdict"], R["ab"]["why"])
        case = TMP / "ns_sigma_sha"
        shutil.copytree(base, case)
        q = json.loads((case / HJ.QUAL).read_text())
        q["sigma_ref_snapshots"][1]["sha256"] = "0" * 64
        (case / HJ.QUAL).write_text(json.dumps(q))
        rc, R = judge(case, ["--no-ab"])
        check("--no-ab: A/B は未実施・参照の最終の資格を確かめていない → 更新に進める条件を満たさない",
              R["ab"]["result"].startswith("未実施") and R["reference_final"]["qualified"] is None and not R["update_conditions_met"])
        rc, R = judge(case)
        check("代わりの参照の sha が違う → σ_ref を測れない → 判定不能・参照の最終の資格なし", R["ab"]["result"] == "判定不能"
              and R["ab"]["sigma_ref"]["status"] == "unmeasurable" and R["ab"]["sigma_ref"]["mm"] is None
              and R["reference_final"]["qualified"] is False and not R["update_conditions_met"])
        # 行が欠けても res があれば抽出し直して補う (結果 A のまま)
        case = TMP / "ns_rows"
        shutil.copytree(base, case)
        rp = case / HJ.NS_RUN / SR.RECORD
        rec = json.loads(rp.read_text())
        rec["rows"] = [r for r in rec["rows"] if r["step"] == S_STOP]
        rp.write_text(json.dumps(rec))
        rc, R = judge(case)
        check("stop_rule.json に停止の行だけ (res は全部ある) → 入るべき出力を列挙して抽出し直して補う → 結果 A",
              R["ab"]["result"] == "結果 A" and R["ab"]["completeness"]["rows_missing"] == req[:-1]
              and R["ab"]["completeness"]["reextracted_B"] == req[:-1], R["ab"]["why"])
        # 反例 2 (組み込み): 停止の行だけ・他の出力の res も無い → 判定不能
        for k in req[:-1]:
            (case / HJ.NS_RUN / f"res_{k}.h5").unlink()
        rc, R = judge(case)
        check(f"反例 2 (組み込み): 停止の 1 行だけで他の res も無い (再抽出で補えない) → 判定不能 ({R['ab']['classify']['code']})",
              R["ab"]["result"] == "判定不能" and R["ab"]["classify"]["code"] == "U_incomplete" and R["ab"]["completeness"]["res_missing"] == req[:-1]
              and not R["update_conditions_met"], R["ab"]["why"])
    # 反例 1 (組み込み): 腕 A の中ほどの 1 列だけ欠測 (F の中) → 結果 A にしない
    def hole(run, step, spec, d):
        if str(Path(spec["run"]).resolve()).endswith(HJ.OLD_EULER):
            d = {**d, "delta_E": list(d["delta_E"]), "missing": list(d["missing"]) + [j_mid]}
            d["delta_E"][j_mid] = float("nan")
        return d
    with Patch((HJ, "HOOKS", {**hooks, "delta_snapshot": snap_hook(hole)})):
        case = TMP / "ns_cx1"
        shutil.copytree(base, case)
        rc, R = judge(case)
        per = R["ab"]["per_output"]
        check(f"反例 1 (組み込み): 腕 A の比べる列の中ほど (x {per[0].get('x_bad_A')}) だけ欠測 → 比べる列は減らさず腕 A は合否 否 → "
              f"代わりの参照の全部でも腕 B は合 → 結果 B ({R['ab']['classify']['code']}、結果 A ではない)",
              R["ab"]["result"] == "結果 B" and R["ab"]["classify"]["code"] == "B_pf" and all(r["pf_A"] is False and r["n_fixed"] == R["ab"]["columns"]["n"]
                                                                                             for r in per), R["ab"]["why"][:200])
    # 反例 3 (組み込み): 途中の 1 出力だけ代わりの参照 (50,000) で δ_E が 0.030 mm 動く → 判定不能・参照の最終の資格なし
    k_mid = req[2]

    def bump(run, step, spec, d):
        if spec.get("res") == "res_50000.h5" and step == k_mid:
            d = {**d, "delta_E": list(d["delta_E"])}
            d["delta_E"][j_mid] += 0.030 / mmpr
        return d
    with Patch((HJ, "HOOKS", {**hooks, "delta_snapshot": snap_hook(bump)})):
        case = TMP / "ns_cx3"
        shutil.copytree(base, case)
        rc, R = judge(case)
        sg = R["ab"]["sigma_ref"]
        at = [a_["max_mm_F"] for a_ in sg["per_step"][str(k_mid)]["alts"] if a_["res"] == "res_50000.h5"]
        check(f"反例 3 (組み込み): 停止の出力 {S_STOP} の感度は小さく、途中の出力 {k_mid} × res_50000 だけ 0.030 mm → σ_ref {sg['mm']:.4f} mm → "
              f"判定不能 (参照の未整定)・参照の最終の資格なし・更新に進めない",
              0.029 < sg["mm"] < 0.031 and at and abs(at[0] - sg["mm"]) < 1e-12 and R["ab"]["classify"]["code"] == "U_sigma_straddle"
              and R["reference_final"]["qualified"] is False and not R["update_conditions_met"], R["ab"]["why"][:200])
    # 比べる列の固定に出口の幾何の許容 (Euler の端の外側だけ、端の列の間隔 × 1e−4) が効く
    xu = np.asarray(json.loads((base / HJ.NS_RUN / SR.RECORD).read_text())["delta_series"]["x"])
    for name, short, want in (("1.68e−9 r_t 短い (変換の丸め)", 1.68e-9, True), ("0.07 r_t 短い (範囲の不足)", 0.07, False)):
        rd = TMP / f"fx_{name[:4]}"
        X, Rr, Q = synth_eu(X1)
        X = X.copy()
        X[-1, :] = xu[-1] - short
        write_field_run(rd, X, Rr, Q, {}, (60000,))
        fx = HJ.fixed_columns(xu, {"B": rd}, edge_tol=True)
        fx0 = HJ.fixed_columns(xu, {"B": rd}, edge_tol=False)
        check(f"比べる列の固定: Euler の端が NS の最後の列より {name} → 最後の列は {'F に入る' if want else '入らない (欠測のまま)'}"
              f" (許容 {fx['cover']['B']['tol_hi']:.2e} r_t; 許容なし (±1e−9) では入らない)",
              bool(fx["mask"][-1]) == want and not fx0["mask"][-1] and fx["mask"][:-2].all())
    with Patch((HJ, "HOOKS", {"handoff_gate": stub_gate(False), "next_wall_check": stub_next_wall})):
        case = TMP / "ns_blocked"
        shutil.copytree(base, case)
        rc, R = judge(case, ["--no-ab"])
        check("関門で止まる → 「更新に進める条件を満たさない: 理由」", rc == 0 and R["verdict"].startswith("更新に進める条件を満たさない")
              and "NOT_CONVERGED" in R["verdict"] and not R["update_conditions_met"])
        case = TMP / "ns_nomain"
        shutil.copytree(base, case)
        (case / HJ.NS_RUN / "H370_MAIN_START.json").unlink()
        (case / HJ.NS_RUN / "DISK_STOP").write_text("試験")
        rc, R = judge(case, ["--no-ab"])
        check("本段の直前の照合の記録が無い・DISK_STOP → 判定不能 (rc 2)", rc == 2 and R["verdict"].startswith("判定不能")
              and any("H370_MAIN_START" in u for u in R["undecidable"]) and any("DISK_STOP" in u for u in R["undecidable"]))
        case = TMP / "ns_nowall"
        shutil.copytree(base, case)
        (case / HJ.NS_RUN / f"res_wall_3_{S_STOP}.h5").unlink()
        rc, R = judge(case, ["--no-ab"])
        check("停止の出力の壁の出力 (res_wall_3_<s>.h5) が無い → 判定不能 (局所 y₁⁺ を判定できない)", rc == 2
              and any("res_wall_3_" in u for u in R["undecidable"]))
        case = TMP / "ns_prov"
        shutil.copytree(base, case)
        (case / HJ.NS_RUN / "RUN_PROVENANCE.txt").write_text("forge_sha256: deadbeef\n")
        rc, R = judge(case, ["--no-ab"])
        check("RUN_PROVENANCE.txt のバイナリ・環境変数が違う → 判定不能", rc == 2 and any("RUN_PROVENANCE" in u for u in R["undecidable"]))


# ====================================================================================== G. h370_run.py の段取り
def fake_prepare(problem, run):
    run = Path(run)
    if run.exists():
        raise FileExistsError(run)
    run.mkdir(parents=True)
    pe = RA.load_problem(problem)
    cfg = RA._apply_gas_to_config(RA._config_euler_node(pe, HR.EU_STEPS, HR.EU_OUT, 4.0, 1), pe, run)
    (run / "solverConfig.yaml").write_text(cfg)
    (run / "bcondConfig.yaml").write_text("wall: {kind: slip}\n")
    for f in ("nozzle.h5", "nozzle.msh"):
        (run / f).write_bytes(b"")
    (run / "MESH_QUALITY.txt").write_text("VERDICT: PASS (AR<= 1000, skew<= 0.90)\n")
    info = {"mesh": {"ni": 2000, "nj": 97, "params": {"throat_refine": 5.0}, "hashes": {"generator_sha256": "g" * 64}},
            "physical_upstream": {"requested": RA._physical_upstream(pe.geometry), "wall": {}}, "scale_m": 0.0766, "nStepOuter": HR.EU_STEPS,
            "euler_wall": {"kind": "design_with_physical_upstream", "delta_r": 0.0, "gate": {"pass": True}, "throat": {"x": 0.0}},
            "initial_line": {"match": {"checked": ["R"], "unchecked": []}}}
    (run / "prepare_info.json").write_text(json.dumps(info))
    return info


def make_fake_prepare_ns(mutate=None):
    def f(problem, run, delta_r_csv=None, offset=None, euler_ref=None, stop_rule=None, **kw):
        run = Path(run)
        if run.exists():
            raise FileExistsError(run)
        p, info = gen_ns_configs(run)
        info.update(delta_r_csv=str(delta_r_csv), stop_rule={"requested": False, "use": False})
        (run / "prepare_info.json").write_text(json.dumps(info))
        (run / "MESH_QUALITY.txt").write_text("VERDICT: PASS (AR<= 5000, skew<= 0.90)\n")
        if mutate:
            mutate(run)
        return info
    return f


def ok_item(*a, **k):
    return {"ok": True, "failures": [], "judged": {"max_abs_mm": 7.96e-4}, "check_mesh_quality_nozzle_h5": "VERDICT: PASS (試験の偽物)"}


def bad_item(*a, **k):
    return {"ok": False, "failures": ["δ_in の照合が不合格 (試験の偽物)"], "judged": {"max_abs_mm": 2e-3}}


def fake_case() -> Path:
    """偽の B の case dir (TMP/repo/case/45x): 問題 YAML・δ_r の CSV・CFD ピンの偽物・h370_judge.py。repo/design・solver_density_cuda は本物へ
    symlink (h370_run.py の HERE.parents[1] の代わり)。"""
    repo = TMP / "repo"
    if not repo.exists():
        (repo / "case").mkdir(parents=True)
        (repo / "design").symlink_to(ROOT / "design")
        (repo / "solver_density_cuda").symlink_to(ROOT / "solver_density_cuda")
    c = repo / "case" / f"c{len(list((repo / 'case').iterdir()))}"
    (c / "_band_ab/half_inlet").mkdir(parents=True)
    for n in (HR.EU_PROBLEM, HR.NS_PROBLEM, HR.CSV, "h370_judge.py"):
        shutil.copy(CASE / n, c / n)
    (c / HR.CFD_PIN_RUN).mkdir()
    for n in HR.CFD_PIN_FILES:
        (c / HR.CFD_PIN_RUN / n).write_text("偽物")
    return c


def args(**kw):
    return types.SimpleNamespace(**{"commit": "abc1234", "dry_run": None, "no_judge": True, "phase": "euler", **kw})


def test_run_orchestration():
    print("== G. h370_run.py の段取り (偽の準備・本物の生成器の設定)", flush=True)
    calls = []
    no_forge = lambda *a, **k: calls.append("forge") or (_ for _ in ()).throw(AssertionError("forge を回した"))  # noqa: E731
    base = [(HR, "resources", lambda *a, **k: {"ok": True, "failures": []}), (HR, "WT", TMP / "no_worktree"),
            (HR, "quality_check", ok_item), (HR, "euler_wall_restore_check", ok_item), (HR, "delta_in_check", ok_item),
            (RA, "prepare", fake_prepare), (RA, "prepare_ns", make_fake_prepare_ns()), (RA, "run_staged", no_forge), (RA, "run_staged_ns", no_forge)]
    c = fake_case()
    with Patch(*base, (HR, "HERE", c)):
        rc = HR.phase_euler(args(dry_run="t1"), RA, {})
        pc = json.loads((c / "_h370_dry_t1_euler" / "H370_PRECHECK.json").read_text())
        check(f"--dry-run: Euler と NS の投入前の検査が合格 (rc {rc})・forge を回さない・run_0950 を作らない",
              rc == 0 and pc["ok"] and not calls and not (c / HR.EU_RUN).exists() and (c / "_h370_dry_t1_ns").is_dir()
              and set(pc["items"]) >= {"E1_wall_restore", "E2_quality", "E3_prepare", "N1_quality", "N2_delta_in", "N3_wall_output", "N4_config", "N5_prepare"},
              str(pc["failures"]))
        check("--dry-run の記録: Euler の段の sha (soft・mid・本段)・NS の段の sha・壁の出力の名前",
              set(pc["items"]["E3_prepare"]["sha256"]) == {"soft", "mid", "main"} and set(pc["items"]["N4_config"]["sha256"]) >= {"S1_soft", "S2_mid", "main"}
              and pc["items"]["N3_wall_output"]["wall_dump_name"] == "res_wall_3_<step>.h5")
        e = raises(lambda: HR.phase_euler(args(dry_run="t1"), RA, {}), SystemExit)
        check("同じ名前の dry-run の準備が既にある → 止める (消さない)", e is not None)
    c = fake_case()
    with Patch(*base, (HR, "HERE", c), (HR, "delta_in_check", bad_item)):
        rc = HR.phase_euler(args(), RA, {})
        check("本番: δ_in の照合が不合格 → 前提の検査の不成立 (rc 2)・PRECHECK_FAILED・forge を回さない",
              rc == 2 and (c / HR.EU_RUN / "PRECHECK_FAILED").is_file() and not calls)
    c = fake_case()
    mut = lambda run: (run / "bcondConfig.yaml").write_text((run / "bcondConfig.yaml").read_text().replace("outputHDFflg: 1, ints: , floats: {Ux", "outputHDFflg: 0, ints: , floats: {Ux"))  # noqa: E731
    with Patch(*base, (HR, "HERE", c), (RA, "prepare_ns", make_fake_prepare_ns(mut))):
        rc = HR.phase_euler(args(), RA, {})
        check("本番: NS が停止の出力で壁の出力を書かない設定 → 終了コード 4 (親に返す)・forge を回さない", rc == HR.EXIT_WALL_OUTPUT and not calls)
    # 合格 → Euler を回す (偽物) → 資格の判定 (本物の h370_judge.py を子プロセスで) が判定不能 → NS に進まない
    c = fake_case()

    def fake_staged(run, **kw):
        calls.append(("run_staged", kw))
        return 0
    with Patch(*base, (HR, "HERE", c), (RA, "run_staged", fake_staged), (RA, "collect", lambda *a: {"偽物": True})):
        rc = HR.phase_euler(args(phase="all"), RA, {})
        rr = json.loads((c / HR.EU_RUN / "RUN_END.json").read_text())
        check(f"本番: 検査の合格 → run_staged を cfl_main 2.0・mid_stage・full で呼ぶ → 資格が無い (偽の Euler) → rc 2 で NS に進まない",
              rc == 2 and calls and calls[-1][1] == {"cfl_main": 2.0, "mid_stage": True, "stages": "full"} and rr["run_rc"] == "0"
              and (c / HR.QUAL).is_file() and not json.loads((c / HR.QUAL).read_text())["prequalified"] and not (c / HR.PIN).exists(), str(calls[-1:]))
        calls.clear()
        e = raises(lambda: HR.phase_ns(args(phase="ns"), RA, {}), SystemExit)
        check("固定が無い → NS を準備しない (止める)", e is not None and not (c / HR.NS_RUN).exists())
    # 環境・prep.done の検査
    with Patch((HR, "HERE", c)):
        old = {k: v for k, v in os.environ.items() if k.startswith("FORGE_")}
        try:
            for k in old:
                os.environ.pop(k)
            os.environ.update(HR.ENV_WANT)
            os.environ["FORGE_EXTRA"] = "1"
            e1 = raises(lambda: HR.code_checks(args()), SystemExit)
            os.environ.pop("FORGE_EXTRA")
            e2 = raises(lambda: HR.code_checks(args()), SystemExit)
        finally:
            for k in [k for k in os.environ if k.startswith("FORGE_")]:
                os.environ.pop(k)
            os.environ.update(old)
        check("FORGE_* が登録と違う → 止める; h370_prep.done が無い → 止める", e1 is not None and e2 is not None)
    # h370_problems.py check を作業ツリーの版で
    c = fake_case()
    (c / "_band_ab").mkdir(exist_ok=True)
    with Patch((HR, "HERE", c), (HR, "WT", ROOT), (HR, "load_forge_design", lambda commit: (RA, {}))):
        e = raises(lambda: HR.check_problems(args(phase="check-problems")), SystemExit)
        check("check-problems: V5d の起動の記録 (moc_v5d_run.log) が B に無い → 止める", e is not None)
        shutil.copy(CASE / "_band_ab" / "moc_v5d_run.log", c / "_band_ab" / "moc_v5d_run.log")
        rc = HR.check_problems(args(phase="check-problems"))
        out = HJ.jload(c / "_band_ab/half_inlet/h370_problems_B.json")
        check("check-problems: 作業ツリーの h370_problems.py check が PASS、出力は h370_problems_B.json (commit した json を上書きしない)",
              rc == 0 and out and out["VERDICT"] == "PASS")


# ====================================================================================== H. NS の段取り・部品
def test_ns_phase_and_parts():
    print("== H. NS の段取り (偽の準備・偽の run_staged_ns で本段の直前の照合)・検査の部品", flush=True)
    from forge_design.feedback.deltastar_loop import handoff_extract_settings
    c = fake_case()
    eb = c / HR.EU_RUN
    eb.mkdir()
    (eb / "res_60000.h5").write_bytes(b"euler")
    sha = HR.sha256_file(eb / "res_60000.h5")
    for k in (40000, 45000, 50000, 55000):
        (eb / f"res_{k}.h5").write_bytes(f"euler {k}".encode())
    (c / HR.PIN).write_text(json.dumps({"run": HR.EU_RUN, "res": "res_60000.h5", "sha256": sha}))
    (c / HR.QUAL).write_text(json.dumps({"prequalified": True, "pin": {"sha256": sha},
                                         "sigma_ref_snapshots": [{"res": f"res_{k}.h5", "step": k, "sha256": HR.sha256_file(eb / f"res_{k}.h5")}
                                                                 for k in (40000, 45000, 50000, 55000)]}))
    fake_ns = make_fake_prepare_ns()

    def fake_prepare_ns_real(problem, run, delta_r_csv=None, offset=None, euler_ref=None, stop_rule=None, **kw):
        info = fake_ns(problem, run, delta_r_csv=delta_r_csv)
        sr = {"requested": True, "use": True, "cap": HR.NS_CAP, "hold_steps": HR.NS_HOLD, "out_step_interval": HR.NS_OUT, **HR.THRESHOLDS,
              "extract": {**handoff_extract_settings(), "source": "試験"}, "euler_ref": euler_ref}
        info.update(stop_rule=sr, euler_ref=euler_ref)
        (Path(run) / "prepare_info.json").write_text(json.dumps(info))
        return info
    seen = {}

    def fake_staged_ns(run, stages="full", **kw):
        """run_staged_ns の段の導出 → 本段の通知 (Observer が照合) だけをまねる (forge は回さない)。"""
        run = Path(run)
        cfg_main = (run / "solverConfig.yaml").read_text()
        cfg_pre = (run / RA.PRESTAGE_CONFIG).read_text()
        pre, main, _ = RA.ns_stage_configs(cfg_main, stages, cfg_pre=cfg_pre)
        done = [{"tag": t, "nsteps": n, "config_sha256": HR.sha256_text(cfg)} for t, cfg, n, _ in pre]
        main = RA._cfg_set(main, {"nStepOuter": HR.NS_CAP})
        (run / "solverConfig.yaml").write_text(main)
        seen["rec"] = RA._write_main_stage_started(run, main, stages, seen.get("done_override", done), cfg_pre)
        return type("RC", (int,), {"stop_rule": {"status": "REACHED"}})(0)
    judge_calls = []
    real_run = subprocess.run

    def spy_run(cmd, *a, **k):
        if any("h370_judge.py" in str(x) for x in cmd):
            judge_calls.append(cmd)
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")
        return real_run(cmd, *a, **k)
    pin_rec = {"run": str(c / HR.EU_RUN), "res": "res_60000.h5", "sha256": sha}
    with Patch((HR, "resources", lambda *a, **k: {"ok": True, "failures": []}), (HR, "WT", TMP / "no_worktree"), (HR, "HERE", c),
               (HR, "quality_check", ok_item), (HR, "delta_in_check", ok_item), (RA, "prepare_ns", fake_prepare_ns_real),
               (RA, "run_staged_ns", fake_staged_ns), (subprocess, "run", spy_run)):
        rc = HR.phase_ns(args(phase="ns", no_judge=False), RA, {})
    run = c / HR.NS_RUN
    pc = json.loads((run / "H370_PRECHECK.json").read_text())
    ms = HJ.jload(run / "H370_MAIN_START.json") or {}
    n7 = pc["items"]["N7_same_as_prelaunch"]
    other = [f for f in pc["failures"] if not f.startswith("(i) の投入前の準備")]
    check(f"NS: (i) の投入前の準備 {HR.NS_PRELAUNCH} が無い → 照合できない → 前提の検査の不成立 (rc {rc})・forge を回さない",
          rc == 2 and not n7["ok"] and not other and (run / "PRECHECK_FAILED").is_file() and "rec" not in seen and not judge_calls, str(other)[:300])
    # (i) の準備があるとき: 合格 → 本段の直前の照合 → 判定を呼ぶ
    c2 = fake_case()
    shutil.copytree(eb, c2 / HR.EU_RUN)
    for f in (HR.PIN, HR.QUAL):
        shutil.copy(c / f, c2 / f)
    pin_rec = {"run": str(c2 / HR.EU_RUN), "res": "res_60000.h5", "sha256": sha}
    with Patch((HR, "HERE", c2)):
        fake_ns(c2 / HR.NS_PROBLEM, c2 / HR.NS_PRELAUNCH, delta_r_csv=str(c2 / HR.CSV))
    base = [(HR, "resources", lambda *a, **k: {"ok": True, "failures": []}), (HR, "WT", TMP / "no_worktree"), (HR, "HERE", c2),
            (HR, "quality_check", ok_item), (HR, "delta_in_check", ok_item), (RA, "prepare_ns", fake_prepare_ns_real),
            (RA, "run_staged_ns", fake_staged_ns), (subprocess, "run", spy_run)]
    with Patch(*base):
        rc = HR.phase_ns(args(phase="ns", no_judge=False), RA, {})
    run = c2 / HR.NS_RUN
    pc = json.loads((run / "H370_PRECHECK.json").read_text())
    ms = HJ.jload(run / "H370_MAIN_START.json") or {}
    check(f"NS: 投入前の検査 (停止規則・Euler の参照の固定・設定・壁の出力・(i) の準備とビット同一) → 本段の直前の照合が合格 (rc {rc})・判定を呼ぶ",
          rc == 0 and pc["ok"] and ms.get("ok") and pc["items"]["N5_prepare"]["euler_ref"] == pin_rec and pc["items"]["N7_same_as_prelaunch"]["ok"]
          and judge_calls, str(pc["failures"] + ms.get("failures", []))[:300])
    check("本段の直前の照合の記録: 本段の設定の sha = 準備の時点で導いた本段 (上限 200,000)・前の段 S1・S2 の sha",
          ms.get("main_stage_started", {}).get("main_config", {}).get("nStepOuter") == HR.NS_CAP
          and [s["tag"] for s in ms["main_stage_started"]["previous_stages"]] == ["S1_soft", "S2_mid"])
    # 負例: 前の段の記録が導いた S1・S2 と違う → WATCH_ABORT・本段を起動させない
    shutil.rmtree(run)
    seen["done_override"] = [{"tag": "S1_soft", "nsteps": 3000, "config_sha256": "0" * 64}]
    with Patch(*base):
        rc = HR.phase_ns(args(phase="ns"), RA, {})
    ms = HJ.jload(run / "H370_MAIN_START.json") or {}
    rr = json.loads((run / "RUN_END.json").read_text())
    check("本段の直前の照合: 前の段の記録が導いた S1・S2 と違う → WATCH_ABORT・例外 (本段の forge を起動させない)",
          not ms.get("ok") and (run / "WATCH_ABORT").is_file() and rr["run_rc"] == "EXC" and "WATCH_ABORT" in rr.get("exception", ""))
    seen.pop("done_override")
    # 部品
    with Patch((HR, "HERE", c2)):
        r = HR.read_pin()
        check("read_pin: 固定 {run, res, sha256} を run の絶対パスで返す (σ_ref の代わりの参照 4 つも在って sha が合う)", r["spec"] == pin_rec)
        (c2 / HR.EU_RUN / "res_45000.h5").unlink()
        check("read_pin: σ_ref の代わりの参照の出力が消されている → NS を回す前に止める", raises(HR.read_pin, SystemExit) is not None)
        (c2 / HR.EU_RUN / "res_45000.h5").write_bytes(b"euler 45000")
        q_ = json.loads((c2 / HR.QUAL).read_text())
        (c2 / HR.QUAL).write_text(json.dumps({**q_, "sigma_ref_snapshots": []}))
        check("read_pin: 予備の資格の記録に σ_ref の代わりの参照が無い → 止める", raises(HR.read_pin, SystemExit) is not None)
        (c2 / HR.QUAL).write_text(json.dumps(q_))
        check("read_pin: 元に戻せば通る", HR.read_pin()["spec"] == pin_rec)
        (c2 / HR.EU_RUN / "res_60000.h5").write_bytes(b"changed")
        check("read_pin: 固定した出力の sha が違う → 止める", raises(HR.read_pin, SystemExit) is not None)
    a, b = TMP / "mi_a.h5", TMP / "mi_b.h5"
    for p_, ro in ((a, 1.0), (b, 2.0)):
        with h5py.File(p_, "w") as f:
            f["MESH/COORD"] = np.arange(6.0)
            f["VALUE/wall_dist"] = np.ones(2)
            f["VALUE/ro"] = np.full(2, ro)
    r = HR.mesh_identity(a, b)
    check("mesh_identity: 格子・wall_dist が同じで初期の場だけ違う → 合格 (初期の場の違いは記録)", r["ok"] and r["initial_field_differ (記録だけ)"] == ["VALUE/ro"])
    with h5py.File(b, "r+") as f:
        f["MESH/COORD"][0] = 1e-12
    check("mesh_identity: 座標が 1 か所でも違う → 不成立", not HR.mesh_identity(a, b)["ok"])
    rd = TMP / "timing_run"
    rd.mkdir()
    (rd / "forge_run.log").write_text("step 1900\nHighFive::DataSetException: Write failed\n")
    HR.Timing(rd).run_once()
    check("Timing: forge_run.log の Write failed → WRITE_FAILED を残す・timing.tsv を書く", (rd / "WRITE_FAILED").is_file() and (rd / "timing.tsv").is_file())
    fake_tools = TMP / "fake_tools"
    fake_tools.mkdir()
    rq = TMP / "q_run"
    rq.mkdir()
    (rq / "MESH_QUALITY.txt").write_text("# cell\nVERDICT: PASS (AR<= 5000, skew<= 0.90)\n")
    for verdict, want in (("PASS (AR<= 5000, skew<= 0.90)", True), ("SOFT-PASS (<0.1% outliers) (AR<= 5000, skew<= 0.90)", False)):
        (fake_tools / "check_mesh_quality.py").write_text(f"print('VERDICT: {verdict}')\n")
        r = HR.quality_check(types.SimpleNamespace(FORGE_TOOLS=fake_tools), rq, 5000)
        check(f"quality_check: check_mesh_quality の VERDICT {verdict.split()[0]} → {'合格' if want else '不成立 (SOFT-PASS は通さない)'}", r["ok"] == want)
    (rq / "MESH_QUALITY.txt").write_text("VERDICT: PASS (AR<= 1000, skew<= 0.90)\n")
    (fake_tools / "check_mesh_quality.py").write_text("print('VERDICT: PASS (AR<= 5000, skew<= 0.90)')\n")
    check("quality_check: prepare の MESH_QUALITY.txt の --ar-max が違う (1000) → 不成立",
          not HR.quality_check(types.SimpleNamespace(FORGE_TOOLS=fake_tools), rq, 5000)["ok"])
    # δ_in の照合 (本物の check_runs、合成の構造格子)
    tbl = np.loadtxt(CASE / HR.CSV, delimiter=",", skiprows=1)
    drx = RA.delta_r_from_table(tbl[:, 0], tbl[:, 1])
    x = np.linspace(-12.5, 95.29, 300)
    rE = 1.0 + 0.0009 * x * x
    for name, shift_mm, want in (("ok", 0.0, True), ("bad", 0.002, False)):
        for kind, rw in (("eu", rE), ("ns", rE + drx(x) + shift_mm * 1e-3 / 0.0767)):
            d = TMP / f"din_{name}_{kind}"
            d.mkdir()
            R2 = rw[:, None] * np.linspace(0.0, 1.0, 5)[None, :]
            with h5py.File(d / "nozzle.h5", "w") as f:
                f["/MESH/COORD"] = np.c_[np.repeat(x, 5) * 0.0767, R2.ravel() * 0.0767, np.zeros(R2.size)].ravel()
            (d / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": 300, "nj": 5}, "scale_m": 0.0767}))
        with Patch((HR, "HERE", CASE)):
            r = HR.delta_in_check(TMP / f"din_{name}_ns", TMP / f"din_{name}_eu")
        check(f"delta_in_check (本物の check_runs): NS の壁 = Euler の壁 + δ_r{'' if want else ' + 0.002 mm'} → {'合格' if want else '不合格'} "
              f"(max {r.get('judged', {}).get('max_abs_mm')} mm)", r["ok"] == want and r.get("csv_sha256") == HR.CSV_SHA, str(r.get("failures")))


def find_pin():
    for c in (CASE / HR.CFD_PIN_RUN, Path("/home/sano/work/forge/case/45.isobutane_m6_d155") / HR.CFD_PIN_RUN):
        if (c / "prepare_info.json").is_file():
            return c.resolve()
    return None


def test_euler_wall_restore():
    print("== I. 正式の Euler の壁の保存・復元の照合 (本物の設計・CFD ピン)", flush=True)
    pin = find_pin()
    if pin is None:
        check("CFD ピンの凍結元 run_0062 が手元に無い — I を回せない", False)
        return
    import h370_problems as HP
    from forge_design.geometry.wall_axismach import save_wall_file
    from forge_design.meshing.mesh2d import generate_axisym_mesh, write_msh41_2d
    td = TMP / "pinned"
    td.mkdir()
    pe = HP.load_pinned(HP.EU, pin, td)
    de = RA.design_chain(pe)
    scale = float(pe.spec["r_throat"])
    wE, er = RA.euler_reference_wall(pe, de, scale)
    c = fake_case()
    run = c / HR.EU_RUN
    run.mkdir()
    mp = RA.mesh_params_euler(pe, scale)
    co, q, b = generate_axisym_mesh(wE, mp)
    write_msh41_2d(run / "nozzle.msh", co, q, b)
    wpath, wrec, wsha = save_wall_file(run, wE, scale, er["repr"])
    info = {"scale_m": scale, "euler_wall": {**er, "wall_file": {"repr": er["repr"], "file": wpath.name, "sha256": wsha}},
            "mesh": {"hashes": {"generator_sha256": HR.generator_sha256(co, q, b)}}}
    (run / "prepare_info.json").write_text(json.dumps(info, default=str))
    shutil.copy(CASE / HR.PREP_CHECK_REC, c / HR.PREP_CHECK_REC)
    with Patch((HR, "HERE", c)):
        r = HR.euler_wall_restore_check(RA, run)
        loc = r.get("local_prep_check (記録だけ)", {})
        check(f"壁ファイル (design_upstream) → 復元した壁で生成した格子の msh が nozzle.msh とバイト一致・生成器の sha が prepare_info と同じ "
              f"(手元の prep-check の記録と同じ: 壁 {loc.get('wall_file_sha256_same')}・格子 {loc.get('generator_sha256_same')})",
              r["ok"] and r["msh_bytes_identical"] and loc.get("wall_file_sha256_same") and loc.get("generator_sha256_same"), str(r.get("failures")))
        m = (run / "nozzle.msh").read_text().splitlines()
        m[len(m) // 2] = m[len(m) // 2] + " "
        (run / "nozzle.msh").write_text("\n".join(m) + "\n")
        r = HR.euler_wall_restore_check(RA, run)
        check("負例: nozzle.msh が 1 文字違う → 不成立", not r["ok"] and not r["msh_bytes_identical"])
        info["euler_wall"]["wall_file"]["sha256"] = "0" * 64
        (run / "prepare_info.json").write_text(json.dumps(info, default=str))
        r = HR.euler_wall_restore_check(RA, run)
        check("負例: 壁ファイルの sha が prepare_info の記録と違う → 不成立", not r["ok"])


def main() -> int:
    for t in (test_constants, test_wall_output_and_configs, test_euler_stage_configs, test_euler_qualify, test_classify_ab, test_ns_judge,
              test_run_orchestration, test_ns_phase_and_parts, test_euler_wall_restore):
        try:
            t()
        except Exception:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            check(f"{t.__name__} が例外", False)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'FAIL {FAIL}'}  (一時ファイル {TMP})")
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
