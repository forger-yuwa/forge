#!/usr/bin/env python3
r"""要求カードの雛形の NS の既定 (本段 B0・監視用の出力・判定の基準) とリミッタの基準値の記録・継承の試験
(plan tooling-settle-convergence-and-monitor §4.5・§4.6・§4.8 の 3 と 5、§5.1 #6、§6 V5)。

    design/.venv-opt/bin/python design/tests/run_settle_defaults_tests.py

(a) 雛形 axismach_final_v1・axismach_explore_v1 の evaluate: main_line (B0、limiter_ref 無し = 自動)・cfl_euler 1.0・ns_output
    (2,000 / 20,000・4 境界・流量)・settle {preset: wt_nozzle_v3} (2026-10-16 ユーザ決定で v1 → v2 (plan §4.9 の 5) → v3 (残差を止める判定から
    外す; plan §4.4・§5.1 #8))、cfl_main 無し。Euler の nStepOuter 12,000・outStepInterval 1,000 は今のまま。
(b) **B0 と出力は NS の本段だけ**: case/71 の problem_m55_p3.yaml (旧い雛形の出力、point・cfl_main 1.0) の evaluate を新しい雛形の evaluate に
    置き換えた写し (= 新しい雛形で生成器が書く evaluate) で、本段の設定に B0 と出力のキー (outStepInterval 2,000・fieldOutStepInterval 20,000・
    output.boundaryMassflux・nodeProbeFile・interfaceDiag) が入り、前の段の元の設定 (solverConfig_prestage.yaml) には無く、S1・S2 は旧い問題の
    S1・S2 と**文字列で一致**、Euler の本段 (と soft・mid) の設定は旧い問題 (cfl_main 1.0 を run_staged で当てる) と文字列で一致。
    bcondConfig は 4 境界の outputHDFflg だけが 1 に変わる。
(c) settle_criteria.yaml と settle_nodes.txt (`write_settle_files`): 手元の case/71 run_0010_ns_b0 の prepare_info と格子の座標で書き、
    `settle.validate_criteria`・`settle.load_criteria` を通る・前提がすべて満たされる・節点の表が forge の書式で範囲内・重複なし・
    grid.segment が本段の開始の通知を指す (手元に run が無ければ skip)。
(d) point の明示 (main_line: point + cfl_main) と不正な main_line。
(e) 拒否と記録: ns_output・settle の不正値、settle: false (外した理由の記録)、前提を満たせない設定 (凝縮) は例外。
(f) リミッタの基準値 (§4.6): from_run の継承 (forge の書式の limiter_ref.json) で 4 つの値がビットで同じ・本段の設定を読み直しても同じ、
    明示の L、記録が無い・別の起動・基準値を使わない起動・stage_manifest と食い違う記録は ValueError。
(g) 容量 (§4.8 の 5): case/71 の規模・上限 20 万 step・間隔 2,000 / 20,000 で、全場の枚数 ≤ ⌈上限/間隔⌉ + 2、全場以外の出力
    (境界・節点・流量) ≤ 300 MB。1 回分の大きさは 2026-10-16 の手元の短い run (下の MEASURED) を FP64 のビルドに換算した値。
(h) 本段の開始の通知 MAIN_STAGE_STARTED (`_write_main_stage_started`): main_config.fnv・settle_criteria の sha256・前の段の limiter_ref.json を消す。
design_chain・変換器・forge は使わない。
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate import settle_presets as SP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402
import settle  # noqa: E402
import stage_manifest as SM  # noqa: E402

CASE = ROOT / "case/71.m55_combustion_wt"
P3 = CASE / "problem_m55_p3.yaml"
P3E = CASE / "problem_m55_p3_euler.yaml"
RUN10 = CASE / "run_0010_ns_b0"
TMP = Path(tempfile.mkdtemp(prefix="settle_defaults_"))
FAIL = 0
SKIPPED = []
B0_ONLY = [("time", "deltaT", k) for k in ("lineImplicit", "lineDtDirectional", "lineDtDirectionalCap", "implicitThermalJacobian")] \
    + [("space", k) for k in ("limiterRoRef", "limiterPRef", "limiterARef", "limiterRefLength")]
OUT_ONLY = [("time", "fieldOutStepInterval"), ("output",)]


def present(cfg, paths):
    """paths のうち cfg に書いてあるもの。#11 の明示すべきキーとして B0 を「切った値」で書いたもの (`_MAIN_LINE_OFF`; 省略と同じ実効値)
    は数えない (plan tooling-generator-defaults-and-launch-gate §5.1 #11・#16、run_main_line_tests.py の b0_left と同じ)。"""
    return [p for p in paths if get(cfg, p) is not RA._CFG_MISSING
            and not (p in RA._MAIN_LINE_OFF and get(cfg, p) == RA._MAIN_LINE_OFF[p])]


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        FAIL += 1


def rejects(name, fn, must=""):
    try:
        fn()
    except ValueError as e:
        check(f"{name} を拒否", must in str(e), str(e)[:140])
        return
    check(f"{name} を拒否 (通ってしまった)", False)


def get(text, path):
    return RA._cfg_path_get(RA._cfg_load(text), path)


def problem(tag, src, edit):
    raw = yaml.safe_load(src.read_text())
    edit(raw)
    path = TMP / f"{tag}.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
    return load_problem(path)


def rd(tag):
    d = TMP / f"rd_{tag}"
    d.mkdir()
    return d


# ------------------------------------------------------------------------------------------------ (a) 雛形
print("== (a) 雛形の evaluate")
REC = {}
for name in ("axismach_final_v1", "axismach_explore_v1"):
    rec, _, _ = W.load_recipe(name)
    ev = rec["evaluate"]
    REC[name] = ev
    check(f"(a) {name}: main_line = B0 (dir・cap 50・itj 5・cfl 4・relax 0.7)、limiter_ref 無し (= 自動)",
          ev.get("main_line") == {"line": "dir", "cap": 50, "itj": 5, "cfl": 4, "relax": 0.7}, ev.get("main_line"))
    check(f"(a) {name}: cfl_main 無し・cfl_euler 1.0 (Euler の cfl を別に保持)", "cfl_main" not in ev and ev.get("cfl_euler") == 1.0)
    check(f"(a) {name}: ns_output (2,000・20,000・4 境界・流量・界面の診断 (成果物の壁の熱流のため))",
          ev.get("ns_output") == {"outStepInterval": 2000, "fieldOutStepInterval": 20000,
                                  "boundary_h5": ["inlet", "outlet", "wall", "axis"], "boundaryMassflux": 1, "interfaceDiag": 1},
          ev.get("ns_output"))
    check(f"(a) {name}: settle = {{preset: wt_nozzle_v3}} (2026-10-16 ユーザ決定「残さは監視から外そう」で v2 から; v1・v2 はプリセットに残る)",
          ev.get("settle") == {"preset": "wt_nozzle_v3"} and "wt_nozzle_v1" in SP.PRESETS and "wt_nozzle_v2" in SP.PRESETS)
    check(f"(a) {name}: Euler の nStepOuter 12000・outStepInterval 1000・mid_stage は今のまま",
          ev.get("nStepOuter") == 12000 and ev.get("outStepInterval") == 1000 and ev.get("mid_stage") is True)


def new_eval(raw, recipe="axismach_final_v1"):
    """生成器が新しい雛形で書く evaluate (雛形の evaluate + 組成に H2O があるときの tp_species)。"""
    tp = raw["evaluate"].get("tp_species")
    raw["evaluate"] = copy.deepcopy(REC[recipe])
    if tp is not None:
        raw["evaluate"]["tp_species"] = tp


# ------------------------------------------------------------------------------------------------ (b) 本段だけ
if not (P3.is_file() and P3E.is_file()):
    SKIPPED.append(f"(b)〜(e) ({P3.relative_to(ROOT)} が無い)")
else:
    print("== (b) B0 と出力は NS の本段だけ")
    p_old = load_problem(P3)
    p_new = problem("new", P3, new_eval)
    check("(b) 旧い問題は point・cfl_main 1.0 (比べる相手の前提)", RA._main_line_spec(p_old) is None and p_old.evaluate.get("cfl_main") == 1.0)
    ml_src = {}
    ml = RA._main_line_spec(p_new, source_out=ml_src)
    check("(b) 新しい問題: B0 (cap 50・itj 5・cfl 4・relax 0.7)・limiter_ref 自動",
          ml is not None and ml["cap"] == 50.0 and ml["itj"] == 5 and ml["cfl"] == 4.0 and ml["relax"] == 0.7
          and ml["limiter_ref"] == {} and ml_src["limiter_ref"]["mode"] == "auto", ml)
    check("(b) _main_line_spec の戻りのキーは今までと同じ (case/45 の台本が登録の dict と比べる)",
          set(ml) == {"line", "cap", "itj", "cfl", "relax", "prestage_relax", "limiter_ref", "extra_fields"}, sorted(ml))
    n_o, oi_o, _ = RA._ns_steps(p_old)
    n_n, oi_n, _ = RA._ns_steps(p_new)
    check("(b) NS の本段の長さは同じ 60000・出力の間隔は 1000 → 2000 (ns_output)", n_o == n_n == 60000 and oi_o == 1000 and oi_n == 2000)
    base_o, main_o = RA._ns_configs(p_old, rd("old"), n_o, oi_o)
    rec_n = {}                     # _ns_configs が足した明示すべきキー {family, prestage, main} (plan generator-defaults §5.1 #11)
    base_n, main_b0 = RA._ns_configs(p_new, rd("new"), n_n, oi_n, main_line=ml, record=rec_n)
    es, srec = RA._settle_spec(p_new)
    no = RA._ns_output_spec(p_new)
    bc_new = RA._bcond_set_output(RA._bcond_with_species(RA._bcond(p_new, euler=False, profile_ok=True), RA._tp_species_Y(p_new)),
                                  no["boundary_h5"])
    pi_min = {"Md": 5.5, "scale_m": float(p_new.spec["r_throat"]), "mesh": {"ni": 4810, "nj": 92}, "x_E": 40.63,
              "physical_wall": {"domain": [-9.98, 79.81]}}       # run_0010_ns_b0 の prepare_info の値 (丸め)
    needs = RA._settle_needs(SP.resolve(es, pi_min, RA._cfg_load(main_b0), yaml.safe_load(bc_new)))
    check("(b) 判定の基準 (v3) が要る出力: 節点・流量 (壁の熱流は判定から外したので界面の診断は要らない; 成果物のために ns_output で出す)",
          needs == {"nodeProbeFile": True, "boundaryMassflux": True, "interfaceDiag": False}, needs)
    check("(b) 雛形の ns_output に interfaceDiag 1 (成果物の壁の熱流のため)", no.get("interfaceDiag") == 1, no)
    main_n = RA._apply_ns_outputs(main_b0, no, needs)
    check("(b) 本段: B0 のキー (ライン陰解法・上限 50・ITJ 5・cfl 4・relax 0.7)",
          get(main_n, ("time", "deltaT", "lineImplicit")) == 1 and get(main_n, ("time", "deltaT", "lineDtDirectional")) == 1
          and get(main_n, ("time", "deltaT", "lineDtDirectionalCap")) == 50.0
          and get(main_n, ("time", "deltaT", "implicitThermalJacobian")) == 5 and get(main_n, ("time", "deltaT", "cfl")) == 4.0
          and get(main_n, ("time", "deltaT", "implicitRelax")) == 0.7)
    check("(b) 本段: リミッタの基準値を書かない (自動)", all(get(main_n, ("space", k)) is RA._CFG_MISSING
                                               for k in ("limiterRoRef", "limiterPRef", "limiterARef", "limiterRefLength")))
    check("(b) 本段: 出力のキー (outStepInterval 2000・fieldOutStepInterval 20000・boundaryMassflux 1・nodeProbeFile・interfaceDiag 1)",
          get(main_n, ("time", "outStepInterval")) == 2000 and get(main_n, ("time", "fieldOutStepInterval")) == 20000
          and get(main_n, ("output",)) == {"boundaryMassflux": 1, "nodeProbeFile": RA.SETTLE_NODES_FILE, "interfaceDiag": 1},
          get(main_n, ("output",)))
    left = [".".join(p) for p in present(base_n, B0_ONLY + OUT_ONLY)]
    check("(b) 前の段の元の設定 (solverConfig_prestage.yaml) に B0・本段の出力のキーが無い", not left, left)
    check("(b) 前の段の元の設定と旧い問題の設定の差は出力の間隔だけ (1000 → 2000)",
          set(RA._cfg_diff_paths(RA._cfg_load(base_o), RA._cfg_load(base_n))) == {("time", "outStepInterval")})
    pre_n, mn, _ = RA.ns_stage_configs(main_n, "full", cfg_pre=base_n)
    pre_o, mo, _ = RA.ns_stage_configs(main_o, "full")
    check("(b) S1・S2 は旧い問題の S1・S2 と文字列で一致 (point・cfl 0.5 / 1.0・1 次・B0 も出力も無し)",
          [c for _, c, *_ in pre_n] == [c for _, c, *_ in pre_o] and len(pre_n) == 2)
    check("(b) 本段 = 出力と B0 を当てた設定", mn == main_n)
    pre_bad, _, _ = RA.ns_stage_configs(main_n, "full")
    check("(b) 対照: 前の段の元の設定を渡さないと S1 に全場の間隔・出力のキー・B0 が残る",
          get(pre_bad[0][1], ("time", "fieldOutStepInterval")) == 20000 and get(pre_bad[0][1], ("output",)) is not RA._CFG_MISSING
          and get(pre_bad[0][1], ("time", "deltaT", "lineImplicit")) == 1)
    rejects("(b) 前の段の元の設定と本段の差が B0・出力・段の長さ以外 (convMethod) なら段の設定を作らない",
            lambda: RA.ns_stage_configs(main_n.replace("convMethod: 1", "convMethod: 2"), "full", cfg_pre=base_n), "convMethod")
    bo = yaml.safe_load(RA._bcond_with_species(RA._bcond(p_old, euler=False, profile_ok=True), RA._tp_species_Y(p_old)))
    bn = yaml.safe_load(bc_new)
    check("(b) bcondConfig: 4 境界の outputHDFflg 1、他は同じ",
          all(bn[k]["outputHDFflg"] == 1 for k in RA.NS_BOUNDARY_NAMES)
          and all({kk: vv for kk, vv in bn[k].items() if kk != "outputHDFflg"} == {kk: vv for kk, vv in bo[k].items() if kk != "outputHDFflg"}
                  for k in bo) and set(bn) == set(bo) and bo["inlet"]["outputHDFflg"] == 0 and bo["axis"]["outputHDFflg"] == 0)
    # Euler: 旧い問題 (設定 4.0 → run_staged の cfl_main 1.0) と新しい問題 (prepare が cfl_euler 1.0 を書き、run_staged は None) の本段
    pe_old = load_problem(P3E)
    pe_new = problem("new_euler", P3E, new_eval)
    check("(b) Euler の cfl: 旧 cfl_main 1.0 / 新 cfl_euler 1.0 (euler_cfl)", RA.euler_cfl(pe_old) == 1.0 == RA.euler_cfl(pe_new)
          and pe_new.evaluate.get("cfl_main") is None)
    ne, oe = int(pe_new.evaluate["nStepOuter"]), int(pe_new.evaluate["outStepInterval"])
    cfg4_o = RA._apply_gas_to_config(RA._config_euler_node(pe_old, ne, oe, 4.0, 1), pe_old, rd("eu_old"))
    cfg4_n = RA._apply_gas_to_config(RA._config_euler_node(pe_new, ne, oe, 4.0, 1), pe_new, rd("eu_new"))
    main_eo = RA._cfg_set(cfg4_o, {"cfl": 1.0, "cfl_pseudo": 1.0})                      # run_staged(cfl_main=1.0)
    import re as _re
    main_en = _re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {RA.euler_cfl(pe_new)}, cfl_pseudo: {RA.euler_cfl(pe_new)}", cfg4_n)
    check("(b) Euler の本段の設定は旧い問題と文字列で一致 (cfl 1.0; 4.0 に変わらない)", main_en == main_eo and get(main_en, ("time", "deltaT", "cfl")) == 1.0)
    soft_o = RA._first_order(RA._cfg_set(main_eo, {"cfl": 0.5, "cfl_pseudo": 0.5}))
    soft_n = RA._first_order(RA._cfg_set(main_en, {"cfl": 0.5, "cfl_pseudo": 0.5}))
    check("(b) Euler の soft・mid の段も一致 (B0・出力のキー無し)", soft_o == soft_n
          and RA._cfg_set(main_eo, {"cfl": 1.0, "cfl_pseudo": 1.0}) == RA._cfg_set(main_en, {"cfl": 1.0, "cfl_pseudo": 1.0})
          and not [p for p in B0_ONLY + OUT_ONLY if get(main_en, p) is not RA._CFG_MISSING])

    # ------------------------------------------------------------------------------------------- (d) point の明示
    print("== (d) point の明示")

    def point_eval(raw):
        new_eval(raw)
        raw["evaluate"]["main_line"] = "point"
        raw["evaluate"]["cfl_main"] = 1.0
    p_pt = problem("point", P3, point_eval)
    check("(d) main_line: point → B0 を使わない (_main_line_spec None)・記録 point (明示)",
          RA._main_line_spec(p_pt) is None and RA._main_line_mode(p_pt) == {"value": "point", "source": "evaluate.main_line: point (明示)"})
    base_p, main_p = RA._ns_configs(p_pt, rd("point"), n_n, oi_n)
    main_p = RA._apply_ns_outputs(main_p, RA._ns_output_spec(p_pt), needs)
    pre_p, _, _ = RA.ns_stage_configs(main_p, "full", cfg_pre=base_p)
    check("(d) point の本段: B0 のキー無し・cfl 1.0・出力のキーあり、S1・S2 は旧い問題と一致",
          not present(main_p, B0_ONLY) and get(main_p, ("time", "deltaT", "cfl")) == 1.0
          and get(main_p, ("time", "fieldOutStepInterval")) == 20000 and [c for _, c, *_ in pre_p] == [c for _, c, *_ in pre_o])

    def bad_ml(raw):
        new_eval(raw)
        raw["evaluate"]["main_line"] = "line"
    rejects("(d) main_line: line (mapping でも point でもない)", lambda: RA._main_line_spec(problem("bad_ml", P3, bad_ml)), "point")

    # ------------------------------------------------------------------------------------------- (e) 拒否と記録
    print("== (e) 拒否と記録")

    def with_ev(**kv):
        def ed(raw):
            new_eval(raw)
            for k, v in kv.items():
                if v is None:
                    raw["evaluate"].pop(k, None)
                else:
                    raw["evaluate"][k] = v
        return ed
    NO = REC["axismach_final_v1"]["ns_output"]
    for i, (nm, v, must) in enumerate((("ns_output の知らないキー", {**NO, "every": 1}, "知らないキー"),
                                       ("fieldOutStepInterval が outStepInterval の倍数でない", {**NO, "fieldOutStepInterval": 3000}, "倍数"),
                                       ("boundary_h5 の知らない名前", {**NO, "boundary_h5": ["inlet", "nozzle"]}, "boundary_h5"),
                                       ("outStepInterval 0", {**NO, "outStepInterval": 0}, "1 以上"))):
        rejects(f"(e) {nm}", lambda v=v, i=i: RA._ns_output_spec(problem(f"bad_no{i}", P3, with_ev(ns_output=v))), must)
    rejects("(e) NS の本段 (ns_steps 61000) が ns_output.outStepInterval の倍数でない",
            lambda: RA._ns_steps(problem("bad_steps", P3, with_ev(ns_steps=61000))), "倍数")
    rejects("(e) settle: true", lambda: RA._settle_spec(problem("st_true", P3, with_ev(settle=True))), "mapping")
    rejects("(e) settle の知らないキー", lambda: RA._settle_spec(problem("st_key", P3, with_ev(settle={"preset": "wt_nozzle_v1", "foo": 1}))),
            "不正")
    rejects("(e) settle の知らないプリセット", lambda: RA._settle_spec(problem("st_pre", P3, with_ev(settle={"preset": "x"}))), "不正")
    es_off, rec_off = RA._settle_spec(problem("st_off", P3, with_ev(settle=False)))
    check("(e) settle: false → 外す・理由の記録", es_off is None and rec_off["status"] == "off" and "明示" in rec_off["reason"])
    es_un, rec_un = RA._settle_spec(problem("st_un", P3, with_ev(settle=None)))
    check("(e) settle が無い → 今までどおり (未指定の記録)", es_un is None and rec_un["status"] == "unset")
    rejects("(e) settle: false + 全場の間隔 20,000 + 停止規則 (判定の基準が無いので旧方式 = 出力ごとに全場を読む) は run dir を作る前に止める",
            lambda: RA.prepare_ns(problem("st_off2", P3, with_ev(settle=False)).path, TMP / "run_st_off"), "旧方式")
    check("(e) その拒否は run dir を作る前", not (TMP / "run_st_off").exists())
    cond_doc = RA._cfg_load(main_b0)
    cond_doc["condensation"] = {"condensation": 1}
    rejects("(e) 前提を満たせない設定 (凝縮のモーメントの列) は残差を判定するプリセット (wt_nozzle_v2) の settle_presets.resolve が例外 (黙って外さない)",
            lambda: SP.resolve({**es, "preset": "wt_nozzle_v2"}, pi_min, cond_doc, yaml.safe_load(bc_new)), "凝縮")
    c3c = SP.resolve(es, pi_min, cond_doc, yaml.safe_load(bc_new))
    check("(e) 雛形の既定 wt_nozzle_v3 は残差の列名を使わないので凝縮でも解決できる (非有限の見張りは全数値列; plan §5.1 #8)",
          es.get("preset") == "wt_nozzle_v3" and "residual_finite" in [q["kind"] for q in c3c["quantities"]]
          and "residual" not in [q["kind"] for q in c3c["quantities"]])
    rejects("(e) 本段の出力を重ねて書けない (既に fieldOutStepInterval がある)", lambda: RA._apply_ns_outputs(main_n, no, needs), "既に")
    check("(e) キー無しの問題は本段の設定を変えない (文字列一致)", RA._apply_ns_outputs(main_o, None, None) == main_o)

# ------------------------------------------------------------------------------------------------ (c) 判定の基準のファイル
print("== (c) settle_criteria.yaml と settle_nodes.txt")
if not (P3.is_file() and (RUN10 / "prepare_info.json").is_file() and (RUN10 / "nozzle.h5").is_file()):
    SKIPPED.append(f"(c) (手元に {RUN10.relative_to(ROOT)} の prepare_info.json・nozzle.h5 が無い)")
else:
    run = rd("crit")
    info10 = json.loads((RUN10 / "prepare_info.json").read_text())
    with h5py.File(RUN10 / "nozzle.h5", "r") as f, h5py.File(run / "nozzle.h5", "w") as g:
        g.create_dataset("/MESH/COORD", data=f["/MESH/COORD"][:])          # 座標だけの写し (元の run には書かない)
    (run / "solverConfig.yaml").write_text(main_n)
    (run / "bcondConfig.yaml").write_text(bc_new)
    es71 = {"preset": "wt_nozzle_v3", "x_range_m": [4.0, 7.95]}   # 雛形の既定 (v3) で作った本段の設定に合わせる (2026-10-16)
    rec = RA.write_settle_files(run, es71, info10, main_n, bc_new)
    text = (run / settle.CRITERIA_FILE).read_text()
    obj = settle._yaml_load_strict(text)
    v = settle.validate_criteria(obj)
    lc = settle.load_criteria(run)
    check("(c) settle_criteria.yaml が settle.validate_criteria・load_criteria を通る", lc["status"] == "ok" and v["grid"]["out_step_interval"] == 2000,
          lc.get("message"))
    check("(c) v3 の量 5 つ (残差の非有限の見張り・流量・出口コア・線・流れ角)、外した量は残差と壁の熱流 (理由つき)",
          rec["quantities"] == ["residual_finite", "mass_balance", "exit_core_mach", "line_mach", "exit_flow_angle"]
          and sorted(e["kind"] for e in obj["excluded"]) == ["residual", "wall_heat_nodes"], (rec["quantities"], obj.get("excluded")))
    check("(c) 前提 (流量・節点の表・界面の診断・壁の h5) がすべて満たされる", all(c["satisfied"] for c in obj["preconditions"]))
    check("(c) grid.segment = 本段の開始の通知 (main_stage_segment)、start_step 0",
          obj["grid"]["segment"] == SP.main_stage_segment() and obj["grid"]["start_step"] == 0)
    ids = SP.read_node_probe_file(run / RA.SETTLE_NODES_FILE)
    ni, nj = int(info10["mesh"]["ni"]), int(info10["mesh"]["nj"])
    check("(c) settle_nodes.txt: 範囲内・重複なし・出口の列の全節点を含む・記録の数と sha256",
          len(ids) == rec["nodes"]["n"] and min(ids) >= 0 and max(ids) < ni * nj and len(set(ids)) == len(ids)
          and set((ni - 1) * nj + j for j in range(nj)) <= set(ids)
          and rec["nodes"]["sha256"] == hashlib.sha256((run / RA.SETTLE_NODES_FILE).read_bytes()).hexdigest(), (len(ids), rec["nodes"]))
    check("(c) 線の x 範囲は問題 YAML の 4.0〜7.95 m", [q for q in obj["quantities"] if q["kind"] == "line_mach"][0]["x_range_m"] == [4.0, 7.95])
    rejects("(c) 本段の設定に節点の表が無い (出力のキーを書き忘れた) → 例外",
            lambda: RA.write_settle_files(rd("crit_bad"), es71, info10, main_b0, bc_new), "nodeProbeFile")

# ------------------------------------------------------------------------------------------------ (f) リミッタの基準値
print("== (f) リミッタの基準値の継承")
VALS = {"limiterRoRef": 0.59298852861870699, "limiterPRef": 192416.41762690683, "limiterARef": 327.13657637494532,
        "limiterRefLength": 9.0296406421583377}


def src_run(tag, cfg_text, computed=True, fnv=None, auto=True):
    """forge の writeLimiterRefRecord と同じ書式の limiter_ref.json と solverConfig.yaml を持つ run (auto = 自動で決めた値か)。"""
    d = rd(tag)
    (d / "solverConfig.yaml").write_text(cfg_text)
    lines = ['{', '  "schema": "forge.limiter_ref/v1",', '  "time": 1,', '  "cfg_fnv": "%s",' % (fnv or SM.fnv1a64(cfg_text)),
             '  "limiterScaled": 1,', '  "limiterDiag": 0,', '  "computed": %s,' % ("true" if computed else "false")]
    if computed:
        lines += ['  "%s": %s,' % (k, "%.17g" % v) for k, v in VALS.items()]
        tf = "true" if auto else "false"
        lines += ['  "auto": {%s},' % ", ".join('"%s": %s' % (k, tf) for k in VALS)]
    else:
        lines += ['  "%s": null,' % k for k in VALS] + ['  "auto": null,']
    lines += ['  "venkatK": 0.050000000000000003,', '  "h_i": "sqrt(A_planar)"', '}']
    (d / "limiter_ref.json").write_text("\n".join(lines) + "\n")
    return d


def strip_explicit(text, added):
    """#11 で足した明示すべきキーの断片 (flow 形式の `, key: 値`) を除く (足す前の設定; run_main_line_tests.py の同名の関数と同じ)。
    B0 の当て方 (`_apply_main_line`) の単体は足す前の設定で見る (`_ns_configs` は B0 を足す前の設定に当ててから明示すべきキーを足す)。"""
    from forge_design.evaluate.config_builder import _token
    for k, v in (added or {}).items():
        frag = f", {k.split('.')[-1]}: {_token(v)}"
        assert text.count(frag) == 1, (k, frag)
        text = text.replace(frag, "", 1)
    return text


if not P3.is_file():
    SKIPPED.append("(f) (problem_m55_p3.yaml が無い)")
else:
    srcA = src_run("srcA", main_b0)
    vals, lrec = RA.limiter_ref_from_run(srcA)
    check("(f) from_run: 4 つの値が limiter_ref.json とビットで同じ・由来の記録 (sha256・cfg_fnv・元は自動)",
          vals == VALS and lrec["mode"] == "from_run" and lrec["cfg_fnv"] == SM.fnv1a64(main_b0)
          and lrec["sha256"] == hashlib.sha256((srcA / "limiter_ref.json").read_bytes()).hexdigest()
          and all(lrec["auto_in_source"].values()), lrec)

    def ml_from(src):
        def ed(raw):
            new_eval(raw)
            raw["evaluate"]["main_line"] = {**raw["evaluate"]["main_line"], "limiter_ref": {"from_run": str(src)}}
        return ed
    src_f = {}
    ml_f = RA._main_line_spec(problem("from_run", P3, ml_from(srcA)), source_out=src_f)
    main_f = RA._apply_main_line(strip_explicit(base_n, rec_n["prestage"]), ml_f)
    check("(f) 継承した本段の設定を読み直すと 4 つの値がビットで同じ (repr → YAML → float の往復)",
          {k: get(main_f, ("space", k)) for k in VALS} == VALS and src_f["limiter_ref"]["mode"] == "from_run"
          and src_f["limiter_ref"]["values"] == VALS)
    check("(f) 継承の本段と自動の本段の差は基準値 4 つだけ",
          set(RA._cfg_diff_paths(RA._cfg_load(main_f), RA._cfg_load(strip_explicit(main_b0, rec_n["main"])))) == {("space", k) for k in VALS})
    # 区間の鍵: 自動の段 (srcA の記録) と継承の段 (固定の同じ値) は同じ区間 (forge の記録を写す形; 詳しくは tools の試験)
    runM = rd("manifest")
    sm = SM.StageManifest(runM)
    # 両方とも明示すべきキーを足す前の設定で比べる (main_f は足す前の前の段の設定に B0 を当てたもの)
    for tag, c, au in (("main_auto", strip_explicit(main_b0, rec_n["main"]), True), ("main_inherit", main_f, False)):
        shutil.copy(src_run(f"src_{tag}", c, auto=au) / "limiter_ref.json", runM / "limiter_ref.json")
        sm.add(tag, c, bc_new if P3.is_file() else "")
    check("(f) 区間の鍵: 自動の段とその値を継承した段は 1 区間 (forge の記録の値が同じ)", len(SM.segments({"stages": sm.stages})) == 1)
    rel = problem("from_run_rel", P3, ml_from(Path("rd_srcA")))
    check("(f) from_run の相対パスは問題 YAML の場所から", RA._main_line_spec(rel)["limiter_ref"] == VALS)
    rejects("(f) from_run の記録が無い", lambda: RA.limiter_ref_from_run(rd("empty")), "")
    srcB = src_run("srcB", main_b0, fnv="deadbeef")
    rejects("(f) from_run の記録が別の起動 (cfg_fnv が solverConfig.yaml と違う)", lambda: RA.limiter_ref_from_run(srcB), "別の起動")
    srcC = src_run("srcC", main_b0, computed=False)
    rejects("(f) from_run の記録が基準値を使わない起動", lambda: RA.limiter_ref_from_run(srcC), "使わない")
    srcD = src_run("srcD", main_b0)
    (srcD / "stage_manifest.json").write_text(json.dumps({"stages": [{"tag": "main", "limiter_ref": {"values": {**VALS, "limiterRoRef": 0.5}}}]}))
    rejects("(f) from_run の記録と stage_manifest の最後の段が食い違う", lambda: RA.limiter_ref_from_run(srcD), "違う")

    def ml_L(raw):
        new_eval(raw)
        raw["evaluate"]["main_line"] = {**raw["evaluate"]["main_line"], "limiter_ref": {"ro": 0.58, "p": 186000.0, "a": 341.0, "L": 9.0}}
    ml_l = RA._main_line_spec(problem("ml_L", P3, ml_L))
    check("(f) 明示の {ro, p, a, L} は limiterRefLength も固定", ml_l["limiter_ref"] == {"limiterRoRef": 0.58, "limiterPRef": 186000.0,
                                                                                 "limiterARef": 341.0, "limiterRefLength": 9.0}
          and get(RA._apply_main_line(strip_explicit(base_n, rec_n["prestage"]), ml_l), ("space", "limiterRefLength")) == 9.0)

    def ml_bad(raw):
        new_eval(raw)
        raw["evaluate"]["main_line"] = {**raw["evaluate"]["main_line"], "limiter_ref": {"ro": 0.58, "p": 186000.0}}
    rejects("(f) limiter_ref の a の欠け", lambda: RA._main_line_spec(problem("ml_bad", P3, ml_bad)), "limiter_ref")

# ------------------------------------------------------------------------------------------------ (g) 容量
print("== (g) 容量 (§4.8 の 5)")
# 1 回分の出力の中身 (2026-10-16、手元の float のビルドで case/71 run_0010_ns_b0 の格子 (442,520 節点・ni 4810・nj 92) に
# この試験 (b) と同じ本段の設定 (prepare_ns が書いた solverConfig.yaml の間隔だけ短くした写し) で 40 step 回した run の h5 を数えた):
#   ファイル: (float の要素数 (FP64 のビルドでは 8 byte)、int32 の要素数、h5 の固定の上乗せ [byte]、float のビルドでの大きさ [byte])
MEASURED = {"res_wall_3": (202020, 19236, 22528, 907552), "res_axis_4": (81770, 19236, 14336, 418360),
            "res_inlet_1": (1748, 364, 14336, 22784), "res_outlet_2": (1564, 364, 14336, 22048),
            "res_nodes": (38320, 3832, 10240, 178848)}
XMF_BYTES = 5658 + 3403 + 3769 + 3401          # 境界の xmf 4 つ (文字列; ビルドに依らない)
CSV_ROW = 200                                   # boundary_massflux.csv の 1 行 (実測 約 150 byte)
FIELD_FP64 = 100_814_972                        # case/71 run_0010_ns_b0/res_86000.h5 (FP64 のビルド)
ok_meas = all(fl * 4 + it * 4 + ov == sz for fl, it, ov, sz in MEASURED.values())
check("(g) 測った大きさ = float 4 byte × 要素 + int32 × 要素 + 固定の上乗せ (換算の式が測定と合う)", ok_meas)
per_out = sum(fl * 8 + it * 4 + ov for fl, it, ov, _ in MEASURED.values()) + XMF_BYTES + CSV_ROW
cap, oi, fi = 200_000, 2000, 20_000
n_out = cap // oi + 1                            # 2,000〜200,000 の 100 回 + 協調停止の step 1 回 (forge は step 0 に境界・節点を書かない)
fields = sorted({0} | set(range(fi, cap + 1, fi)) | {cap} | {cap - 1})   # 最初・間隔・最後・格子でない停止の step
check(f"(g) 全場の枚数 {len(fields)} ≤ ⌈上限/間隔⌉ + 2 = {math.ceil(cap / fi) + 2} (FP64 で {len(fields) * FIELD_FP64 / 1e9:.2f} GB)",
      len(fields) <= math.ceil(cap / fi) + 2)
check(f"(g) 全場以外の出力 {n_out} 回 × {per_out / 1e6:.3f} MB = {n_out * per_out / 1e6:.1f} MB ≤ 300 MB (FP64 に換算)",
      n_out * per_out <= 300e6, n_out * per_out)
print(f"     内訳 (FP64、1 回分): " + ", ".join(f"{k} {(fl * 8 + it * 4 + ov) / 1e6:.3f} MB" for k, (fl, it, ov, _) in MEASURED.items())
      + f", xmf {XMF_BYTES / 1e6:.3f} MB; 余裕 {(300e6 - n_out * per_out) / 1e6:.1f} MB")

# ------------------------------------------------------------------------------------------------ (h) 本段の開始の通知
print("== (h) 本段の開始の通知 MAIN_STAGE_STARTED")
if not P3.is_file():
    SKIPPED.append("(h) (problem_m55_p3.yaml が無い)")
else:
    runH = rd("started")
    (runH / "solverConfig.yaml").write_text(main_n)
    (runH / RA.PRESTAGE_CONFIG).write_text(base_n)
    (runH / "nozzle.h5").write_bytes(b"initial field")
    (runH / settle.CRITERIA_FILE).write_text("schema: forge.settle_criteria/v1\n")
    (runH / "limiter_ref.json").write_text('{"stale": "S2"}')
    recH = RA._write_main_stage_started(runH, main_n, "full", [{"tag": "S1_soft"}, {"tag": "S2_mid"}], base_n)
    note = json.loads((runH / RA.MAIN_STAGE_STARTED).read_text())
    check("(h) 通知に本段の設定の FNV-1a 64 (forge_launches・limiter_ref.json・stage_manifest と同じ定義)",
          note["main_config"]["fnv"] == SM.fnv1a64(main_n) and note == recH)
    check("(h) 通知に判定の基準の sha256・step 0 から", note["settle_criteria"] == {
        "file": settle.CRITERIA_FILE, "sha256": hashlib.sha256((runH / settle.CRITERIA_FILE).read_bytes()).hexdigest()}
        and note["start_step"] == 0)
    check("(h) 前の段の limiter_ref.json を消した (本段の forge が書き直す)",
          not (runH / "limiter_ref.json").exists() and note["limiter_ref"]["previous_removed"] is True)
    runH2 = rd("started_nocrit")
    (runH2 / "solverConfig.yaml").write_text(main_o)
    (runH2 / "nozzle.h5").write_bytes(b"x")
    n2 = RA._write_main_stage_started(runH2, main_o, "full", [], None)
    check("(h) 判定の基準の無い run は settle_criteria null・limiter_ref.previous_removed false",
          n2["settle_criteria"] is None and n2["limiter_ref"]["previous_removed"] is False)

shutil.rmtree(TMP, ignore_errors=True)
for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
