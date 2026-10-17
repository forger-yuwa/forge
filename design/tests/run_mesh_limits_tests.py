#!/usr/bin/env python3
"""格子の節点数の上限の承認と格子の段の契約の試験 (plan tooling-design-orchestrator §5.1 #8・§6 V2a〜V2c の事前登録 2026-10-18;
2026-10-18 codex 諮問 mesh-limits)。GPU・forge なし (変換器と格子の生成・品質の検査は回す)。

  V2a 手組みの格子を包んだ段の変更の前後の A/B: A = 0161e51c の `hand_mesh.py` (台本だけを一時に写す)、B = 今の `hand_mesh.py build --approval`。
      同じ入力 (case/71 の版 2 のカード + mesh_limits {443000, 467000}、`problem_m55_p3_card.yaml`、--axis-cap-frac 0.03、同じ変換器)。
      問題 YAML の mesh・mesh_euler (max_nodes_reason 以外)・Euler の写し、`check` の NS・Euler の nozzle.msh (バイト)、NS の nozzle.h5 の
      MESH/COORD・MESH/CONNE・BCONDS/* が一致。参考: case/71 の problem_m55_p3.yaml の mesh 節との差を出す (合否にしない)。V6 の再現とは呼ばない。
  V2b 雛形の格子の段 (`campaign.mesh_stage.check_template`): case/70 の版 2 のカード (mesh_limits なし = 既定 30 万) → 生成 → 承認 → NS・Euler の
      品質 PASS・節点数 ≤ 承認の上限・記録に承認元と問題 YAML の sha256・格子のハッシュ。部品の検証。
  V2c 負例 (生成の関数の呼び出し 0 回・出力の問題 YAML 0 件) と正例:
      (a) 版 2 のカードの overrides に max_nodes と過去の理由 → カードの検査で誤り
      (b) case/71 のカードで mesh_limits なし (既定 30 万) の承認 → build が問題 YAML を書かずに止まる
      (c) mesh_limits {400000, 467000} の承認 → 止まる
      (d) 上限 30 万の承認に max_nodes 443000 と過去の理由の問題 YAML → approval.check MISMATCH
      (e) case/70 の承認を case/71 のカードに → 止まる (最初の全格子の生成の前)
      (f) 判別の A/B: 試行のカードの conditions.Pt だけを変える。A (承認の値) は estimate_mesh に到達、B (2 倍) は到達 0 回で止まる
      (g) Euler だけが上限を超える (mesh_euler.axis_gap_frac) → 問題 YAML を書かずに止まる
      (h) 格子の段の後に Euler の写しの max_nodes を上げる → recipe の Euler の段が準備の前に止まる (runner を呼ばない)
      正例: 較正の係数だけを変えた試行のカードは止まらない (承認の照合 OK)

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_mesh_limits_tests.py [--skip-v2b]
所要: 10 分程度 (V2a の build と check を 2 回ずつ、V2b の生成と準備)。
"""
import copy
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DESIGN = REPO / "design"
# 準備の経路 (V2b) は初期場に化学種の属性を書くのに --resolve-species を持つ forge が要る (起動はしない)。runner の起動の環境は import の
# 時点の写しなので import の前に置く
_FB = REPO / "solver_density_cuda" / ".build-native" / "card4a" / "forge"
if not os.environ.get("FORGE_BIN") and _FB.is_file():
    os.environ["FORGE_BIN"] = str(_FB)
sys.path.insert(0, str(DESIGN))
from forge_design.campaign import mesh_stage as MS  # noqa: E402
from forge_design.campaign.recipes import wt_card_chain as WC  # noqa: E402
from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

C71 = REPO / "case" / "71.m55_combustion_wt"
C70 = REPO / "case" / "70.m5_air_wt_demo"
HM_PATH = C71 / "_mesh" / "hand_mesh.py"
RATIO = C71 / "_mesh" / "c45" / "h_over_est_ratio_summary.json"
CARD71, PROB71_IN, PROB71_P5 = C71 / "card_m55_3coef.yaml", C71 / "problem_m55_p3_card.yaml", C71 / "problem_m55_p5_settle_v2.yaml"
CARD70 = C70 / "card_m5_air_demo.yaml"
A_COMMIT = "0161e51c"                      # #8 の前 (事前登録のとおり; af510413 も hand_mesh.py は同じ)
CONV = os.environ.get("FORGE_CONVERTER") or str(REPO / "solver_density_cuda" / ".build-native" / "card4a" / "convertGmshToForge")
PY = sys.executable
FAIL = 0
SKIPPED = []
USER_SECTIONS = ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall", "tolerances", "p_ambient")


def check(name, ok, detail=""):
    global FAIL
    print(("ok  " if ok else "FAIL") + " " + name + (f"  [{str(detail)[:400]}]" if detail else ""), flush=True)
    if not ok:
        FAIL += 1


def v2_card(src: Path, dst: Path, mesh_limits=None, ref="試験の仮の承認", edit=None) -> Path:
    d = yaml.safe_load(src.read_text())
    d["format"] = W.CARD_FORMAT_V2
    d["budget"] = {"gpu_hours": 8.0, "disk_gb": 50.0}
    src_ = {k: {"kind": "user", "ref": ref} for k in USER_SECTIONS + ("budget",) if k in d or k == "budget"}
    if mesh_limits is not None:
        d["mesh_limits"] = mesh_limits
        src_["mesh_limits"] = {"kind": "user", "ref": ref + " (節点数の上限)"}
    d["sources"] = src_
    if edit:
        edit(d)
    dst.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    return dst


CLEAN = {}


def clean_repo() -> Path:
    """判定器の版の記録に使う、きれいな一時の git リポジトリ (この試験は作業ツリーの状態に依らない; run_approval_tests.py と同じ)。"""
    if "r" not in CLEAN:
        r = td / "clean_repo"
        r.mkdir()
        for c in (["init", "-q"], ["-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q", "--allow-empty", "-m", "t"]):
            subprocess.run(["git", "-C", str(r), *c], check=True, capture_output=True)
        CLEAN["r"] = r
    return CLEAN["r"]


def approve(card: Path, prob: Path, out: Path) -> Path:
    pv = AP.preview(card, prob, repo_dir=clean_repo(), judge_env=False)
    if pv["problems"]:
        raise RuntimeError(f"承認の preview に問題: {pv['problems']}")
    pp = out.with_suffix(".preview.json")
    pp.write_text(json.dumps(pv, ensure_ascii=False))
    AP.write(pp, "試験の仮の承認 (run_mesh_limits_tests)", out)
    return out


def trial_card(card: Path, approval: Path, dst: Path, edit=None) -> Path:
    d = yaml.safe_load(card.read_text())
    d["approval"] = {"file": str(approval), "sha256": AP.sha256_file(approval)}
    if edit:
        edit(d)
    dst.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    return dst


def build_cmd(script: Path, card: Path, prob: Path, out: Path, rec: Path, approval: Path | None = None):
    c = [PY, str(script), "build", "--design", str(DESIGN), "--card", str(card), "--problem", str(prob), "--out-yaml", str(out),
         "--record", str(rec), "--axis-cap-frac", "0.03", "--ratio-summary", str(RATIO), "--plan-item", "V2a"]
    return c + (["--approval", str(approval)] if approval is not None else [])


def load_hm():
    spec = importlib.util.spec_from_file_location("hand_mesh_under_test", HM_PATH)
    hm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hm)
    hm._imports(str(DESIGN))
    return hm


class Reached(Exception):
    pass


def build_inproc(card: Path, prob: Path, out: Path, approval: Path | None, *, stub_estimate=False, force_full_grid=None,
                 stub_generate=False):
    """hand_mesh.build を同じプロセスで呼び、最初の全格子の生成 (wt_card.estimate_mesh) と格子の生成 (mesh2d.generate_axisym_mesh) の
    呼び出しを数える。stub_estimate なら estimate_mesh を到達の記録だけにする (Reached で抜ける)。戻り (例外 or None, 呼び出しの数)。"""
    hm = load_hm()
    calls = {"estimate_mesh": 0, "generate_axisym_mesh": 0}
    orig_est, orig_gen = hm.W.estimate_mesh, hm.M2.generate_axisym_mesh

    def est(*a, **k):
        calls["estimate_mesh"] += 1
        if stub_estimate:
            raise Reached()
        if force_full_grid is not None:
            k["full_grid"] = force_full_grid
        return orig_est(*a, **k)

    def gen(*a, **k):
        calls["generate_axisym_mesh"] += 1
        if stub_generate:
            raise Reached()
        return orig_gen(*a, **k)
    hm.W.estimate_mesh, hm.M2.generate_axisym_mesh = est, gen
    import argparse
    a = argparse.Namespace(mode="build", design=str(DESIGN), card=str(card), problem=str(prob), ratio_summary=str(RATIO),
                           out_yaml=str(out), record=str(out.with_suffix(".rec.json")), run_dir=None, axis_cap_frac=0.03,
                           radial_ratio_max=1.2, plan_item="V2c", approval=None if approval is None else str(approval), sanity=None)
    err = None
    try:
        hm.build(a)
    except Exception as e:                       # noqa: BLE001
        err = e
    finally:
        hm.W.estimate_mesh, hm.M2.generate_axisym_mesh = orig_est, orig_gen
    return err, calls


def h5_sets(p: Path) -> dict:
    out = {}
    with h5py.File(p, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and (name.startswith("MESH/") or name.startswith("BCONDS/")):
                out[name] = np.asarray(obj[()])
        f.visititems(visit)
    return out


td = Path(tempfile.mkdtemp(prefix="mesh_limits_tests_"))
try:
    if not (PROB71_IN.is_file() and Path(CONV).is_file()):
        print(f"飛ばした: {PROB71_IN} か変換器 {CONV} が無い")
        sys.exit(0)
    os.environ["FORGE_CONVERTER"] = CONV
    card71 = v2_card(CARD71, td / "card71.yaml", {"max_nodes_ns": 443000, "max_nodes_euler": 467000},
                     ref="case/71 のユーザ決定「a」2026-10-15")
    ap71 = approve(card71, PROB71_P5, td / "approval71.json")

    print("== V2a 手組みの格子を包んだ段の変更の前後の A/B")
    A = td / "A"
    B = td / "B"
    A.mkdir()
    B.mkdir()
    a_script = A / "hand_mesh.py"
    a_script.write_text(subprocess.run(["git", "-C", str(REPO), "show", f"{A_COMMIT}:case/71.m55_combustion_wt/_mesh/hand_mesh.py"],
                                       capture_output=True, text=True, check=True).stdout)
    ra = subprocess.run(build_cmd(a_script, card71, PROB71_IN, A / "problem.yaml", A / "build_record.json"), capture_output=True, text=True)
    rb = subprocess.run(build_cmd(HM_PATH, card71, PROB71_IN, B / "problem.yaml", B / "build_record.json", ap71), capture_output=True, text=True)
    check("V2a 前提: A・B の build が通る", ra.returncode == 0 and rb.returncode == 0,
          f"A rc {ra.returncode} {ra.stderr[-300:]} / B rc {rb.returncode} {rb.stderr[-300:]}")
    ya, yb = (yaml.safe_load((x / "problem.yaml").read_text()) for x in (A, B))

    def strip(y):
        y = copy.deepcopy(y)
        for blk in ("mesh", "mesh_euler"):
            y[blk].pop("max_nodes_reason", None)
        return y
    check("V2a 問題 YAML: max_nodes_reason 以外の全キーが一致 (max_nodes も 443000・467000)",
          strip(ya) == strip(yb) and yb["mesh"]["max_nodes"] == 443000 and yb["mesh_euler"]["max_nodes"] == 467000,
          sorted(k for k in set(ya) | set(yb) if strip(ya).get(k) != strip(yb).get(k)))
    check("V2a B の理由は承認元 (承認の記録の名前と sha256)、A の理由は固定の文",
          "承認の記録 approval71.json" in yb["mesh"]["max_nodes_reason"] and "ユーザ決定「a」" in ya["mesh"]["max_nodes_reason"],
          yb["mesh"]["max_nodes_reason"][:160])
    ea, eb = (x / "problem_euler.yaml" for x in (A, B))
    check("V2a Euler の写し: max_nodes_reason 以外が一致", ea.is_file() and eb.is_file()
          and strip(yaml.safe_load(ea.read_text())) == strip(yaml.safe_load(eb.read_text())))
    ca = subprocess.run([PY, str(a_script), "check", "--design", str(DESIGN), "--card", str(card71), "--problem", str(A / "problem.yaml"),
                         "--record", str(A / "check_record.json"), "--run-dir", str(A / "chk")], capture_output=True, text=True)
    cb = subprocess.run([PY, str(HM_PATH), "check", "--design", str(DESIGN), "--card", str(card71), "--problem", str(B / "problem.yaml"),
                         "--record", str(B / "check_record.json"), "--run-dir", str(B / "chk")], capture_output=True, text=True)
    check("V2a 前提: A・B の check が通る", ca.returncode == 0 and cb.returncode == 0, f"{ca.stderr[-300:]} / {cb.stderr[-300:]}")
    same_msh = all((A / "chk" / t / "nozzle.msh").read_bytes() == (B / "chk" / t / "nozzle.msh").read_bytes() for t in ("ns", "euler"))
    check("V2a NS・Euler の nozzle.msh (座標・接続・境界の辺) がバイトで一致", same_msh)
    for tag, fn in (("NS", "ns/nozzle.h5"), ("NS (cell の写し)", "ns/nozzle_qc.h5"), ("Euler (cell の写し)", "euler/nozzle_qc.h5")):
        fa, fb = A / "chk" / fn, B / "chk" / fn
        ha, hb = (h5_sets(fa), h5_sets(fb)) if fa.is_file() and fb.is_file() else ({}, {"無い": 1})
        check(f"V2a {tag} の変換後の {Path(fn).name} の MESH/*・BCONDS/* のデータセットが一致",
              bool(ha) and ha.keys() == hb.keys() and all(np.array_equal(ha[k], hb[k]) for k in ha), f"{len(ha)} 個")
    hm = load_hm()
    arrs = {}
    for tag, x in (("A", A), ("B", B)):
        pr, dd, _wd, wp, _ini = hm.walls(str(x / "problem.yaml"))
        rt_ = float(pr.spec["r_throat"])
        ns = hm.M2.generate_axisym_mesh(wp, hm.RA.mesh_params(pr, rt_, ni=561, nj=97, wall_first_frac=4.5e-5))
        we, _ew = hm.RA.euler_reference_wall(pr, dd, rt_, pwu=hm.RA._pw_upstream(pr.geometry))
        eu = hm.M2.generate_axisym_mesh(we, hm.RA.mesh_params_euler(pr, rt_))
        arrs[tag] = (ns, eu)

    def same_gen(a, b):
        return all(np.array_equal(np.asarray(x), np.asarray(y)) for x, y in zip(a[:2], b[:2])) and repr(a[2]) == repr(b[2])
    check("V2a 生成した直後の NS・Euler の格子 (座標・接続・境界の辺) がビットで一致 (check と同じ構築)",
          same_gen(arrs["A"][0], arrs["B"][0]) and same_gen(arrs["A"][1], arrs["B"][1]),
          f"NS {arrs['A'][0][0].shape} Euler {arrs['A'][1][0].shape}")
    ref = yaml.safe_load((C71 / "problem_m55_p3.yaml").read_text())
    diffs = sorted(f"{blk}.{k}" for blk in ("mesh", "mesh_euler") for k in set(ref[blk]) | set(yb[blk])
                   if ref[blk].get(k) != yb[blk].get(k))
    print(f"参考 (合否にしない): case/71 problem_m55_p3.yaml (10-15 の版) の mesh 節との差 {diffs}", flush=True)

    print("== V2c 負例と正例")
    d70 = v2_card(CARD70, td / "card70_bad.yaml", edit=lambda d: d.update(
        overrides={"mesh.max_nodes": 443000, "mesh.max_nodes_reason": "ユーザ決定「a」2026-10-15: M5.5 の場所ごとの第一層"}))
    try:
        W.load_card(d70)
        ok, why = False, "通った"
    except ValueError as e:
        ok, why = "書けない" in str(e), str(e)
    check("(a) 版 2 のカードの overrides に max_nodes と過去の理由 → カードの検査で誤り", ok, why[:200])

    card71_def = v2_card(CARD71, td / "card71_default.yaml", None, ref="case/71 (上限なし)")
    ap_def = approve(card71_def, PROB71_P5, td / "approval71_default.json")
    for key, card_, apv in (("b", card71_def, ap_def),):
        out = td / f"{key}" / "problem.yaml"
        out.parent.mkdir()
        err, calls = build_inproc(card_, PROB71_IN, out, apv)
        check("(b) mesh_limits なし (既定 30 万) の承認 → 問題 YAML を書かずに止まる (格子の生成 0 回)",
              err is not None and "書かずに止めた" in str(err) and not out.exists() and not any(out.parent.glob("*.tmp*"))
              and calls["generate_axisym_mesh"] == 0, f"{type(err).__name__}: {str(err)[:200]} {calls}")
    card71_400 = v2_card(CARD71, td / "card71_400.yaml", {"max_nodes_ns": 400000, "max_nodes_euler": 467000}, ref="case/71 (小さい上限)")
    # 承認の preview は問題 YAML の物理・評価を記録するので、照合に使う p5 で作る (上限は check でしか見ない)
    ap400 = approve(card71_400, PROB71_P5, td / "approval71_400.json")
    out = td / "c" / "problem.yaml"
    out.parent.mkdir()
    err, calls = build_inproc(card71_400, PROB71_IN, out, ap400)
    check("(c) mesh_limits {400000, 467000} の承認 → 止まる (承認した上限を超える作り直し)",
          err is not None and "書かずに止めた" in str(err) and not out.exists() and calls["generate_axisym_mesh"] == 0,
          f"{str(err)[:200]} {calls}")
    oA, oB = td / "c_fgA" / "problem.yaml", td / "c_fgB" / "problem.yaml"
    oA.parent.mkdir()
    oB.parent.mkdir()
    eA, cA = build_inproc(card71_400, PROB71_IN, oA, ap400, force_full_grid=True, stub_generate=True)
    eB, cB = build_inproc(card71_400, PROB71_IN, oB, ap400, force_full_grid=False, stub_generate=True)
    check("判別の A/B (codex v2c-fail): (c) の入力で full_grid だけ変える。A (True) は生成に 1 回到達、B (False) は生成 0 回・"
          "442,520 > 400,000 の上限の誤り・問題 YAML と一時のファイル 0 件",
          isinstance(eA, Reached) and cA["generate_axisym_mesh"] == 1 and cB["generate_axisym_mesh"] == 0 and eB is not None
          and "442520" in str(eB) and "400000" in str(eB) and not oB.exists() and not any(oB.parent.glob("*.tmp*")),
          f"A {type(eA).__name__} {cA} / B {type(eB).__name__} {cB} {str(eB)[:160]}")
    r = AP.check(ap_def, card71_def, PROB71_P5)
    check("(d) 上限 30 万の承認に max_nodes 443000 と過去の理由の問題 YAML → approval.check MISMATCH",
          r["status"] == "MISMATCH" and any("mesh_limits.max_nodes_ns" in p for p in r["problems"]), r["problems"][:2])
    card70 = v2_card(CARD70, td / "card70.yaml", None, ref="case/70 の要求")
    pb70 = td / "problem70.yaml"
    gen70 = subprocess.run([PY, "-m", "forge_design.intake.wt_card", str(card70), "--out", str(pb70)], cwd=DESIGN, capture_output=True, text=True)
    check("V2b・(e) 前提: case/70 の版 2 のカードから問題 YAML を生成", gen70.returncode == 0 and pb70.is_file(), gen70.stderr[-300:])
    ap70 = approve(card70, pb70, td / "approval70.json")
    out = td / "e" / "problem.yaml"
    out.parent.mkdir()
    err, calls = build_inproc(card71, PROB71_IN, out, ap70, stub_estimate=True)
    check("(e) case/70 の承認を case/71 のカードに → 止まる (最初の全格子の生成 0 回)",
          isinstance(err, AP.ApprovalError) and calls["estimate_mesh"] == 0 and not out.exists(), f"{type(err).__name__}: {str(err)[:160]}")
    tA = trial_card(card71, ap71, td / "trial_A.yaml")
    tB = trial_card(card71, ap71, td / "trial_B.yaml", edit=lambda d: d["conditions"].update(Pt=2 * float(d["conditions"]["Pt"])))
    outA, outB = td / "fA" / "problem.yaml", td / "fB" / "problem.yaml"
    outA.parent.mkdir()
    outB.parent.mkdir()
    eA, cA = build_inproc(tA, PROB71_IN, outA, ap71, stub_estimate=True)
    eB, cB = build_inproc(tB, PROB71_IN, outB, ap71, stub_estimate=True)
    check("(f) 判別の A/B: A (Pt は承認の値) は estimate_mesh に到達、B (Pt × 2) は到達 0 回・問題 YAML 0 件で止まる",
          isinstance(eA, Reached) and cA["estimate_mesh"] == 1 and isinstance(eB, AP.ApprovalError) and cB["estimate_mesh"] == 0
          and "conditions.Pt" in str(eB) and not outB.exists(), f"A {type(eA).__name__} {cA} / B {type(eB).__name__} {cB} {str(eB)[:160]}")
    tC = trial_card(card71, ap71, td / "trial_C.yaml",
                    edit=lambda d: d["overrides"].update({"deltastar_initializer.delta_scale": 0.70}))
    ck = AP.check_card(ap71, tC)
    check("正例: 較正の係数だけを変えた試行のカードは承認の照合 OK (止まらない)", ck["status"] == "OK" and ck["kind"] == "trial", ck["problems"])
    pg = td / "problem_gap.yaml"
    dg = yaml.safe_load(PROB71_IN.read_text())
    dg["mesh_euler"]["axis_gap_frac"] = 0.005
    pg.write_text(yaml.safe_dump(dg, sort_keys=False, allow_unicode=True))
    shutil.copy(C71 / "problem_m55_p3_card_euler.yaml", td / "problem_gap_euler.yaml")
    out = td / "g" / "problem.yaml"
    out.parent.mkdir()
    err, calls = build_inproc(card71, pg, out, ap71)
    check("(g) Euler だけが上限を超える (mesh_euler.axis_gap_frac 0.005) → 問題 YAML を書かずに止まる (格子の生成 0 回)",
          err is not None and "書かずに止めた" in str(err) and "mesh_euler" in str(err) and not out.exists()
          and calls["generate_axisym_mesh"] == 0, f"{str(err)[:240]} {calls}")
    eu_bad = td / "h_problem_euler.yaml"
    de = yaml.safe_load((B / "problem_euler.yaml").read_text())
    de["mesh_euler"]["max_nodes"] = 900000
    eu_bad.write_text(yaml.safe_dump(de, sort_keys=False, allow_unicode=True))
    runs = []
    orig_run = WC._run
    WC._run = lambda *a, **k: runs.append(a) or subprocess.CompletedProcess(a, 0, "", "")
    try:
        res = WC.exec_step({"step": "euler_p1", "k": 1, "case_dir": str(td), "repo": str(REPO), "approval": str(ap71), "card": str(card71),
                            "problem_euler": str(eu_bad), "run_dir": str(td / "h_run"), "out": str(td / "h.json")})
    finally:
        WC._run = orig_run
    check("(h) 格子の段の後に Euler の写しの max_nodes を上げる → recipe の Euler の段が準備の前に止まる (runner を呼ばない)",
          res["ok"] is False and not runs and "節点数の上限" in (res["detail"] or ""), res["detail"])

    if "--skip-v2b" in sys.argv[1:]:
        SKIPPED.append("V2b (--skip-v2b)")
    else:
        print("== V2b 雛形の格子の段 (case/70)")
        rec = MS.check_template(pb70, td / "ms70", approval=ap70, card=card70)
        check("V2b NS・Euler とも品質 PASS、節点数 ≤ 承認の上限、記録に承認元・問題 YAML の sha256・格子のハッシュ",
              all((rec[t]["quality"] or "").startswith("VERDICT: PASS") for t in ("ns", "euler"))
              and rec["nodes"]["ns"] <= rec["limits"]["max_nodes_ns"] and rec["nodes"]["euler"] <= rec["limits"]["max_nodes_euler"]
              and rec["ns"]["n_nodes"] == rec["nodes"]["ns"] and rec["euler"]["n_nodes"] == rec["nodes"]["euler"]
              and rec["limits"]["approval_sha256"] == AP.sha256_file(ap70) and rec["problem_sha256"] and rec["ns"]["hashes"]
              and (td / "ms70" / MS.RECORD).is_file(),
              {t: (rec[t]["quality"], rec[t]["n_nodes"]) for t in ("ns", "euler")})
finally:
    shutil.rmtree(td, ignore_errors=True)

for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
