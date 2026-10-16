#!/usr/bin/env python3
r"""系列のリミッタの基準値の自動の継承の試験 (plan tooling-settle-convergence-and-monitor §4.6・§5.1 #6b、2026-10-16)。

    design/.venv-opt/bin/python design/tests/run_limiter_inherit_tests.py

`evaluate.main_line.limiter_ref` の省略は、prepare_ns の初期場の元 (ic_from) が NS の run ならその基準値の継承 (`limiter_ref_inherit_from_ic`)、
Euler の run・元が無いときは自動。`auto` の明示は継承しない自動。明示 ({ro, p, a[, L]})・from_run は初期場の元より優先。B0 の本段だけ。
(a) 継承: 元の run の limiter_ref.json (cfg_fnv が本段の設定と一致・computed) の 4 つ — 本段の設定に書いて読み直すとビットで同じ・出所の記録
    (ファイル・sha256・cfg_fnv・値)。記録が無い (固定値だけ; run_0010 の形) → 本段の設定の明示 ρ・P・a を継承し limiterRefLength は書かない
    (L も明示なら 4 つ)。記録の cfg_fnv が本段の設定と違う → 固定値に回して理由を記録。stage_manifest は tag main の段と照らす (最後の段でない)。
(b) 明示の優先: {ro, p, a}・from_run は、基準値の分からない NS の元でも止まらず明示の値。
(c) 止まり方 (ValueError、黙って自動に戻さない): 記録も ρ・P・a の明示も無い NS の元、記録が壊れている、computed false、limiterScaled 0 の
    明示、本段の設定の明示と記録の食い違い、MAIN_STAGE_STARTED の本段の設定と solverConfig.yaml の食い違い (sha256・fnv)、stage_manifest の
    main の段との食い違い、NS/Euler の判定の食い違い・決められない元、S1・S2 で落ちた run (MAIN_STAGE_STARTED が無く stage_manifest の
    最後の段が main でない)。MAIN_STAGE_STARTED も stage_manifest も無い旧い run は solverConfig.yaml を本段とみなし「本段の確認なし」と記録。
(d) `auto` の明示: 基準値の分からない NS の元でも自動 (継承しない) と記録。
(e) Euler からの初期場・元が無い・point の本段は自動 (継承を見ない)。
(f) 本段の開始の通知 MAIN_STAGE_STARTED の `limiter_ref.source` に出所を写し、本段の設定の値との一致を書く。
(g) 手元の case/71 run_0010_ns_b0 (固定値 ρ・P・a、L は自動、記録なし) と run_0009_euler (あれば)。
design_chain・変換器・forge は使わない (prepare_ns の試行は別に scratchpad で)。
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate import runner_wt as RW  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402
import stage_manifest as SM  # noqa: E402

CASE = ROOT / "case/71.m55_combustion_wt"
P3 = CASE / "problem_m55_p3.yaml"
RUN10 = CASE / "run_0010_ns_b0"
RUN09 = CASE / "run_0009_euler"
TMP = Path(tempfile.mkdtemp(prefix="limiter_inherit_"))
FAIL = 0
SKIPPED = []
ML = {"line": "dir", "cap": 50, "itj": 5, "cfl": 4, "relax": 0.7}
# 記録の値 (forge の %.17g の書式) と、run_0008 の固定値 (case/71 problem_m55_p3_b0.yaml の limiter_ref)
VALS = {"limiterRoRef": 0.59298852861870699, "limiterPRef": 192416.41762690683, "limiterARef": 327.13657637494532,
        "limiterRefLength": 9.0296406421583377}
FIX = {"limiterRoRef": 0.5810042584, "limiterPRef": 186732.713, "limiterARef": 341.8295597}
NAMES = SM.LIMITER_REF_NAMES


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        FAIL += 1


def rejects(name, fn, must=""):
    try:
        fn()
    except ValueError as e:
        check(f"{name} を拒否", must in str(e), str(e)[:200])
        return
    check(f"{name} を拒否 (通ってしまった)", False)


def get(text, path):
    return RA._cfg_path_get(RA._cfg_load(text), path)


def problem(tag, ml):
    """P3 の写しで evaluate を B0 (main_line = ml; None なら main_line を消す) にした問題。"""
    raw = yaml.safe_load(P3.read_text())
    ev = raw["evaluate"]
    ev.pop("cfl_main", None)
    ev.pop("implicit_relax", None)
    if ml is None:
        ev.pop("main_line", None)
    else:
        ev["main_line"] = ml
    path = TMP / f"{tag}.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
    return load_problem(path)


def with_space(cfg, vals):
    """設定の space の先頭に基準値を書く (_apply_main_line と同じ書き方)。"""
    if not vals:
        return cfg
    return cfg.replace("space: {", "space: {" + ", ".join(f"{k}: {float(v)!r}" for k, v in vals.items()) + ", ", 1)


def record_text(cfg_text, computed=True, fnv=None, auto=True, vals=None):
    """forge の writeLimiterRefRecord と同じ書式の limiter_ref.json。"""
    vals = VALS if vals is None else vals
    lines = ['{', '  "schema": "forge.limiter_ref/v1",', '  "time": 1,', '  "cfg_fnv": "%s",' % (fnv or SM.fnv1a64(cfg_text)),
             '  "limiterScaled": %d,' % (1 if computed else 0), '  "limiterDiag": 0,', '  "computed": %s,' % ("true" if computed else "false")]
    if computed:
        lines += ['  "%s": %s,' % (k, "%.17g" % vals[k]) for k in NAMES]
        tf = "true" if auto else "false"
        lines += ['  "auto": {%s},' % ", ".join('"%s": %s' % (k, tf) for k in NAMES)]
    else:
        lines += ['  "%s": null,' % k for k in NAMES] + ['  "auto": null,']
    lines += ['  "venkatK": 0.050000000000000003,', '  "h_i": "sqrt(A_planar)"', '}']
    return "\n".join(lines) + "\n"


def src_run(tag, cfg_text, chain="axismach_ns", record=None, started="match", manifest=None):
    """初期場の元の run: solverConfig.yaml (本段の設定)・prepare_info.json (chain)・limiter_ref.json (record: 文字列)・
    MAIN_STAGE_STARTED (started: match = 本段の設定の sha256 と fnv / sha_only = #6 より前の形 / bad_sha / bad_fnv / None = 書かない)・
    stage_manifest.json (manifest: 段の list)。res_*.h5 は置かない (継承の判定は読まない)。"""
    d = TMP / f"src_{tag}"
    d.mkdir()
    (d / "solverConfig.yaml").write_text(cfg_text)
    if chain is not None:
        (d / "prepare_info.json").write_text(json.dumps({"chain": chain}))
    if record is not None:
        (d / "limiter_ref.json").write_text(record)
    sha = hashlib.sha256(cfg_text.encode()).hexdigest()
    mc = {"match": {"sha256": sha, "fnv": SM.fnv1a64(cfg_text)}, "sha_only": {"sha256": sha},
          "bad_sha": {"sha256": "0" * 64}, "bad_fnv": {"sha256": sha, "fnv": "deadbeef"}}.get(started)
    if mc is not None:
        (d / RA.MAIN_STAGE_STARTED).write_text(json.dumps({"event": "main_stage_started", "main_config": {"file": "solverConfig.yaml", **mc}}))
    if manifest is not None:
        (d / "stage_manifest.json").write_text(json.dumps({"manifest_version": 3, "stages": manifest}))
    return d


def spec(ml_lr, ic, tag, drop=False):
    """main_line (limiter_ref = ml_lr; drop なら書かない) の _main_line_spec(ic_from=ic) と出所の記録。"""
    ml = dict(ML)
    if not drop:
        ml["limiter_ref"] = ml_lr
    src = {}
    out = RA._main_line_spec(problem(tag, ml), source_out=src, ic_from=ic)
    return out, src.get("limiter_ref")


if not P3.is_file():
    SKIPPED.append(f"(a)〜(f) ({P3.relative_to(ROOT)} が無い)")
else:
    p0 = load_problem(P3)
    BASE = RW._config_sst_node(p0, 60000, 2000, 1.0)       # B0 を当てる前の NS の設定 (point の本段と同じ形)
    NS_AUTO = RA._apply_main_line(BASE, {**ML, "limiter_ref": {}, "extra_fields": None})   # 元の run の本段: 基準値は自動
    NS_FIX = RA._apply_main_line(BASE, {**ML, "limiter_ref": FIX, "extra_fields": None})   # 元の run の本段: ρ・P・a の固定 (run_0010 の形)
    NS_FIX4 = RA._apply_main_line(BASE, {**ML, "limiter_ref": VALS, "extra_fields": None})  # 4 つの固定
    EULER = RW._config_euler_node(p0, 12000, 1000, 1.0)

    # ------------------------------------------------------------------------------------------------ (a) 継承
    print("== (a) 継承")
    sA = src_run("rec", NS_AUTO, record=record_text(NS_AUTO))
    mlA, lrA = spec(None, sA, "a_rec", drop=True)
    mainA = RA._apply_main_line(BASE, mlA)
    check("(a) 記録あり: 4 つを継承 (limiter_ref.json の %.17g とビットで同じ)", mlA["limiter_ref"] == VALS, mlA["limiter_ref"])
    check("(a) 記録あり: 本段の設定を読み直しても 4 つがビットで同じ", {k: get(mainA, ("space", k)) for k in NAMES} == VALS)
    check("(a) 記録あり: 出所 (mode inherit・basis limiter_ref.json・元の run・sha256・cfg_fnv・値・元は自動・L 固定)",
          lrA["mode"] == "inherit" and lrA["basis"] == "limiter_ref.json" and lrA["run"] == str(sA)
          and lrA["sha256"] == hashlib.sha256((sA / "limiter_ref.json").read_bytes()).hexdigest()
          and lrA["cfg_fnv"] == SM.fnv1a64(NS_AUTO) and lrA["values"] == VALS and all(lrA["auto_in_source"].values())
          and lrA["L_ref"].startswith("fixed") and lrA["ic_kind"]["kind"] == "ns"
          and lrA["main_config"]["main_stage_started"] == {"sha256_checked": True, "fnv_checked": True}, lrA)
    sF = src_run("fix", NS_FIX, started="sha_only")
    mlF, lrF = spec(None, sF, "a_fix", drop=True)
    mainF = RA._apply_main_line(BASE, mlF)
    check("(a) 固定値だけ (記録なし・run_0010 の形): ρ・P・a を継承", mlF["limiter_ref"] == FIX, mlF["limiter_ref"])
    check("(a) 固定値だけ: 本段の設定に limiterRefLength を書かない (forge が新しい run の格子から決める)",
          get(mainF, ("space", "limiterRefLength")) is RA._CFG_MISSING and {k: get(mainF, ("space", k)) for k in FIX} == FIX)
    check("(a) 固定値だけ: 出所 (basis solverConfig.yaml・sha256・記録が無い理由・L auto・MAIN_STAGE_STARTED の sha256 だけ照合)",
          lrF["mode"] == "inherit" and lrF["basis"] == "solverConfig.yaml"
          and lrF["sha256"] == hashlib.sha256(NS_FIX.encode()).hexdigest() and "無い" in lrF["record_skipped"]
          and lrF["L_ref"].startswith("auto") and lrF["main_config"]["main_stage_started"] == {"sha256_checked": True, "fnv_checked": False}
          and lrF["main_confirmed"] == f"{RA.MAIN_STAGE_STARTED} の main_config", lrF)
    sF4 = src_run("fix4", NS_FIX4, started=None)
    mlF4, lrF4 = spec(None, sF4, "a_fix4", drop=True)
    check("(a) 4 つの固定 (MAIN_STAGE_STARTED も stage_manifest も無い旧い run): 4 つを継承・L 固定・出所に「本段の確認なし」",
          mlF4["limiter_ref"] == VALS and lrF4["L_ref"].startswith("fixed") and lrF4["main_config"]["main_stage_started"] is None
          and lrF4["main_confirmed"].startswith(RA.IC_MAIN_UNCONFIRMED) and lrF4["main_config"]["main_confirmed"] == lrF4["main_confirmed"], lrF4)
    sO = src_run("other", NS_FIX, record=record_text(NS_FIX, fnv="0123456789abcdef"))
    mlO, lrO = spec(None, sO, "a_other", drop=True)
    check("(a) 記録が別の起動 (cfg_fnv が本段の設定と違う): 固定値に回し、理由と記録の sha256 を残す",
          mlO["limiter_ref"] == FIX and lrO["basis"] == "solverConfig.yaml" and "起動のものでない" in lrO["record_skipped"]
          and hashlib.sha256((sO / "limiter_ref.json").read_bytes()).hexdigest() in lrO["record_skipped"], lrO)
    m_ok = [{"tag": "S2_mid", "cfg_fnv": "x", "limiter_ref": {"status": "forge", "values": {**VALS, "limiterRoRef": 0.5}}},
            {"tag": "main", "order": 2, "cfg_fnv": SM.fnv1a64(NS_AUTO), "limiter_ref": {"status": "forge", "values": VALS}},
            {"tag": "post", "cfg_fnv": "y", "limiter_ref": {"status": "forge", "values": {**VALS, "limiterPRef": 1.0}}}]
    sM = src_run("man_ok", NS_AUTO, record=record_text(NS_AUTO), manifest=m_ok)
    mlM, lrM = spec(None, sM, "a_man", drop=True)
    check("(a) stage_manifest は tag main の段と照らす (前後の段の値が違っても通る)",
          mlM["limiter_ref"] == VALS and lrM["stage_manifest"] == {"checked": True, "tag": "main", "order": 2, "status": "forge"}, lrM)
    sM2 = src_run("man_nomain", NS_AUTO, record=record_text(NS_AUTO), manifest=m_ok[:1])
    check("(a) stage_manifest に main の段が無ければ照らさず理由を記録 (MAIN_STAGE_STARTED があれば本段は確かめ済み; run_0010 の形)",
          spec(None, sM2, "a_man2", drop=True)[1]["stage_manifest"]["checked"] is False)
    m_main = [{"tag": "S1_soft", "cfg_fnv": "a"}, {"tag": "S2_mid", "cfg_fnv": "b"},
              {"tag": "main", "order": 2, "cfg_fnv": SM.fnv1a64(NS_FIX), "limiter_ref": {"status": "unknown", "values": None}}]
    sMM = src_run("nostart_main", NS_FIX, started=None, manifest=m_main)
    mlMM, lrMM = spec(None, sMM, "a_nostart_main", drop=True)
    check("(a) MAIN_STAGE_STARTED が無く stage_manifest の最後の段が main (cfg_fnv 一致) → 本段と確かめて固定値を継承",
          mlMM["limiter_ref"] == FIX and lrMM["main_confirmed"].startswith("stage_manifest") and lrMM["stage_manifest"]["checked"] is False, lrMM)

    # ------------------------------------------------------------------------------------------------ (b) 明示の優先
    print("== (b) 明示の優先")
    sU = src_run("unknown", NS_AUTO)                         # NS・記録なし・明示なし (継承すれば止まる元)
    mlE, lrE = spec({"ro": 0.6, "p": 190000.0, "a": 340.0}, sU, "b_expl")
    check("(b) 明示の {ro, p, a} は初期場の元を見ない (止まらない・mode explicit)",
          mlE["limiter_ref"] == {"limiterRoRef": 0.6, "limiterPRef": 190000.0, "limiterARef": 340.0} and lrE["mode"] == "explicit", lrE)
    mlR, lrR = spec({"from_run": str(sA)}, sU, "b_from")
    check("(b) from_run は初期場の元より優先 (from_run の記録の 4 つ・mode from_run)",
          mlR["limiter_ref"] == VALS and lrR["mode"] == "from_run" and lrR["run"] == str(sA), lrR)
    mlR2, _ = spec({"from_run": str(sA)}, sF, "b_from2")
    check("(b) from_run は継承できる別の元 (固定値) があってもそちら", mlR2["limiter_ref"] == VALS)

    # ------------------------------------------------------------------------------------------------ (c) 止まり方
    print("== (c) 元の基準値が分からないときは止まる")
    rejects("(c) NS の元に記録も ρ・P・a の明示も無い", lambda: spec(None, sU, "c_unknown", drop=True), "基準値が分からない")
    rejects("(c) 明示が ρ・P だけ (a が欠ける)",
            lambda: spec(None, src_run("pa_only", RA._apply_main_line(BASE, {**ML, "extra_fields": None,
                                                                               "limiter_ref": {k: FIX[k] for k in ("limiterRoRef", "limiterPRef")}})),
                         "c_pa", drop=True), "limiterARef")
    rejects("(c) 記録が壊れている (JSON でない) — 固定値があっても落とさない",
            lambda: spec(None, src_run("broken", NS_FIX, record="{not json"), "c_broken", drop=True), "壊れている")
    rejects("(c) 記録の schema が違う",
            lambda: spec(None, src_run("schema", NS_FIX, record=json.dumps({"schema": "x"})), "c_schema", drop=True), "壊れている")
    rejects("(c) 記録が computed false (元の本段が基準値を使わない)",
            lambda: spec(None, src_run("nocomp", NS_AUTO, record=record_text(NS_AUTO, computed=False)), "c_nocomp", drop=True), "使わない")
    rejects("(c) 記録なしで本段の設定が limiterScaled 0",
            lambda: spec(None, src_run("scaled0", NS_AUTO.replace("space: {", "space: {limiterScaled: 0, ", 1)), "c_sc0", drop=True), "limiterScaled 0")
    rejects("(c) 本段の設定の明示と記録の値の食い違い",
            lambda: spec(None, src_run("clash", NS_FIX, record=record_text(NS_FIX)), "c_clash", drop=True), "食い違う")
    rejects("(c) MAIN_STAGE_STARTED の本段の設定の sha256 が solverConfig.yaml と違う",
            lambda: spec(None, src_run("bad_sha", NS_FIX, started="bad_sha"), "c_bsha", drop=True), "MAIN_STAGE_STARTED")
    rejects("(c) MAIN_STAGE_STARTED の本段の設定の fnv が solverConfig.yaml と違う",
            lambda: spec(None, src_run("bad_fnv", NS_FIX, started="bad_fnv"), "c_bfnv", drop=True), "fnv")
    m_bad = [{"tag": "main", "cfg_fnv": SM.fnv1a64(NS_AUTO), "limiter_ref": {"status": "forge", "values": {**VALS, "limiterARef": 300.0}}}]
    rejects("(c) stage_manifest の main の段の値と継承する値が違う",
            lambda: spec(None, src_run("man_bad", NS_AUTO, record=record_text(NS_AUTO), manifest=m_bad), "c_man", drop=True), "stage_manifest")
    m_bad2 = [{"tag": "main", "cfg_fnv": "deadbeef", "limiter_ref": {"status": "forge", "values": VALS}}]
    rejects("(c) stage_manifest の main の段の cfg_fnv が本段の設定と違う",
            lambda: spec(None, src_run("man_badfnv", NS_AUTO, record=record_text(NS_AUTO), manifest=m_bad2), "c_man2", drop=True), "cfg_fnv")
    rejects("(c) NS/Euler の判定の食い違い (設定は SST・prepare_info は Euler の chain)",
            lambda: spec(None, src_run("kind_clash", NS_FIX, chain="axismach"), "c_kind", drop=True), "食い違う")
    rejects("(c) NS か Euler かを決められない (solverConfig.yaml も prepare_info.json も無い)",
            lambda: spec(None, (TMP / "bare").mkdir() or (TMP / "bare"), "c_bare", drop=True), "決められない")
    rejects("(c) prepare_info が NS の chain で solverConfig.yaml が無い",
            lambda: spec(None, (lambda d: (d.mkdir(), (d / "prepare_info.json").write_text('{"chain": "axismach_ns"}'), d)[2])(TMP / "pi_only"),
                         "c_pionly", drop=True), "solverConfig.yaml が無い")
    rejects("(c) limiter_ref の不正な文字列 (auto 以外)", lambda: spec("inherit", sU, "c_str"), "auto")
    # S1・S2 で落ちた run: MAIN_STAGE_STARTED が無く (run_staged_ns は最初に消す)、solverConfig.yaml は前の段 (point・基準値は自動) の設定のまま、
    # limiter_ref.json はその段の起動のもの (cfg_fnv が一致・computed) — 以前はこの値を本段のものとして継承できてしまった
    m_s1 = [{"tag": "S1_soft", "order": 0, "cfg_fnv": SM.fnv1a64(BASE), "limiter_ref": {"status": "forge", "values": VALS}}]
    sS1 = src_run("s1_failed", BASE, started=None, record=record_text(BASE), manifest=m_s1)
    rejects("(c) S1 で落ちた run (MAIN_STAGE_STARTED なし・stage_manifest の最後の段 S1_soft・記録は S1 の起動)",
            lambda: spec(None, sS1, "c_s1", drop=True), "main でない")
    sS2 = src_run("s2_failed", BASE, started=None, record=record_text(BASE),
                  manifest=m_s1 + [{"tag": "S2_mid", "order": 1, "cfg_fnv": SM.fnv1a64(BASE), "limiter_ref": {"status": "forge", "values": VALS}}])
    rejects("(c) S2 で落ちた run (最後の段 S2_mid)", lambda: spec(None, sS2, "c_s2", drop=True), "main でない")
    sMB = src_run("nostart_mainfnv", NS_FIX, started=None, manifest=[{"tag": "main", "cfg_fnv": "deadbeef"}])
    rejects("(c) MAIN_STAGE_STARTED が無く stage_manifest の最後の段 main の cfg_fnv が solverConfig.yaml と違う",
            lambda: spec(None, sMB, "c_mainfnv", drop=True), "cfg_fnv")

    # ------------------------------------------------------------------------------------------------ (d) auto の明示
    print("== (d) auto の明示")
    mlD, lrD = spec("auto", sU, "d_auto")
    check("(d) auto: 基準値の分からない NS の元でも止まらず自動 (limiter_ref {}・mode auto・明示の理由・元の run)",
          mlD["limiter_ref"] == {} and lrD["mode"] == "auto" and "明示" in lrD["reason"] and lrD["ic_from"] == str(sU), lrD)
    mlD2, _ = spec("auto", sA, "d_auto_rec")
    mainD = RA._apply_main_line(BASE, mlD2)
    check("(d) auto: 記録のある NS の元でも継承せず、本段の設定に基準値を書かない",
          mlD2["limiter_ref"] == {} and all(get(mainD, ("space", k)) is RA._CFG_MISSING for k in NAMES))
    mlD3, lrD3 = spec("auto", sS1, "d_auto_s1")
    check("(d) auto: S1 で落ちた元でも止まらず自動", mlD3["limiter_ref"] == {} and lrD3["mode"] == "auto", lrD3)
    mlD4, lrD4 = spec({"ro": 0.6, "p": 190000.0, "a": 340.0}, sS1, "d_expl_s1")
    check("(d) 明示: S1 で落ちた元でも明示の値 (初期場の元を見ない)", lrD4["mode"] == "explicit" and mlD4["limiter_ref"]["limiterRoRef"] == 0.6, lrD4)

    # ------------------------------------------------------------------------------------------------ (e) 自動のまま
    print("== (e) Euler の元・元なし・point は自動")
    sE = src_run("euler", EULER, chain="axismach")
    mlEu, lrEu = spec(None, sE, "e_euler", drop=True)
    check("(e) Euler からの初期場は自動 (limiter_ref {}・mode auto・理由 Euler・判定の根拠)",
          mlEu["limiter_ref"] == {} and lrEu["mode"] == "auto" and "Euler" in lrEu["reason"] and lrEu["ic_kind"]["kind"] == "euler"
          and set(lrEu["ic_kind"]["by"]) == {"solverConfig.yaml", "prepare_info.json"}, lrEu)
    sE2 = src_run("euler_nopi", EULER, chain=None)
    check("(e) Euler の判定は solverConfig.yaml だけでも (prepare_info.json が無い)", spec(None, sE2, "e_euler2", drop=True)[1]["mode"] == "auto")
    sE3 = src_run("euler_otherchain", EULER, chain="wind_tunnel")
    check("(e) 知らない chain は判定に使わない (solverConfig.yaml で Euler)", spec(None, sE3, "e_euler3", drop=True)[1]["mode"] == "auto")
    mlN, lrN = spec(None, None, "e_none", drop=True)
    check("(e) 初期場の元が無い (新しい問題) は自動", mlN["limiter_ref"] == {} and lrN["mode"] == "auto" and "無い" in lrN["reason"], lrN)
    p_pt = problem("e_point", "point")
    check("(e) point の本段は初期場の元を見ない (基準値の分からない NS の元でも None)",
          RA._main_line_spec(p_pt, cfl_main=1.0, ic_from=sU) is None)
    check("(e) main_line の無い問題も None", RA._main_line_spec(problem("e_noml", None), ic_from=sU) is None)
    src_k = {}
    RA._main_line_spec(problem("e_compat", ML), source_out=src_k)
    check("(e) ic_from を渡さない呼び出し (case/45 の台本の照合) は今までどおり自動", src_k["limiter_ref"]["mode"] == "auto")

    # ------------------------------------------------------------------------------------------------ (f) 通知
    print("== (f) MAIN_STAGE_STARTED の出所")
    for tag, ml_lr, src_dir, main_cfg_fn, want in (
            ("rec", None, sA, lambda ml: RA._apply_main_line(BASE, ml), True),
            ("fix", None, sF, lambda ml: RA._apply_main_line(BASE, ml), True),
            ("auto", "auto", sU, lambda ml: RA._apply_main_line(BASE, ml), True),
            ("edited", None, sF, lambda ml: RA._apply_main_line(BASE, ml).replace("limiterRoRef: 0.5810042584", "limiterRoRef: 0.58", 1), False)):
        mlx, lrx = spec(ml_lr, src_dir, f"f_{tag}", drop=(ml_lr is None))
        run = TMP / f"run_f_{tag}"
        run.mkdir()
        cfg = main_cfg_fn(mlx)
        (run / "solverConfig.yaml").write_text(cfg)
        (run / "nozzle.h5").write_bytes(b"x")
        (run / "prepare_info.json").write_text(json.dumps({"main_line": {"limiter_ref_source": lrx}}, default=str))
        note = RA._write_main_stage_started(run, cfg, "full", [], None)
        on_disk = json.loads((run / RA.MAIN_STAGE_STARTED).read_text())
        check(f"(f) {tag}: limiter_ref.source = prepare_info の出所、本段の設定との一致 {want}",
              on_disk == note and note["limiter_ref"]["source"] == json.loads(json.dumps(lrx, default=str))
              and note["limiter_ref"]["source_matches_config"] is want, note["limiter_ref"])
    runP = TMP / "run_f_nopi"
    runP.mkdir()
    (runP / "nozzle.h5").write_bytes(b"x")
    (runP / "solverConfig.yaml").write_text(BASE)
    nP = RA._write_main_stage_started(runP, BASE, "full", [], None)
    check("(f) prepare_info の無い run は source null・一致 null", nP["limiter_ref"]["source"] is None
          and nP["limiter_ref"]["source_matches_config"] is None and nP["limiter_ref"]["config_values"] == {})

# ------------------------------------------------------------------------------------------------ (g) 手元の case/71
print("== (g) 手元の case/71 の run")
if not (RUN10 / "solverConfig.yaml").is_file():
    SKIPPED.append(f"(g) run_0010 ({RUN10.relative_to(ROOT)} が無い)")
else:
    v10, r10 = RA.limiter_ref_inherit_from_ic(RUN10)
    check("(g) run_0010_ns_b0: ρ・P・a = run_0008 の固定値、L は自動 (記録が無い・#6 より前)",
          v10 == FIX and r10["basis"] == "solverConfig.yaml" and r10["L_ref"].startswith("auto") and "無い" in r10["record_skipped"], (v10, r10))
    check("(g) run_0010_ns_b0: MAIN_STAGE_STARTED の sha256 で本段の設定を確かめた (fnv の欄は #6 より前で無い)・NS の判定は 2 つの根拠で一致",
          r10["main_config"]["main_stage_started"] == {"sha256_checked": True, "fnv_checked": False}
          and r10["main_confirmed"] == f"{RA.MAIN_STAGE_STARTED} の main_config"
          and r10["ic_kind"]["kind"] == "ns" and {v["kind"] for v in r10["ic_kind"]["by"].values()} == {"ns"}, r10)
    check("(g) run_0010_ns_b0: stage_manifest に main の段が無い (S1・S2 だけ) ので照らさない", r10["stage_manifest"]["checked"] is False, r10)
if not (RUN09 / "solverConfig.yaml").is_file():
    SKIPPED.append(f"(g) run_0009 ({RUN09.relative_to(ROOT)} が無い)")
else:
    v09, r09 = RA.limiter_ref_inherit_from_ic(RUN09)
    check("(g) run_0009_euler: 自動 (Euler)", v09 == {} and r09["mode"] == "auto" and r09["ic_kind"]["kind"] == "euler", r09)

shutil.rmtree(TMP, ignore_errors=True)
for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
