#!/usr/bin/env python3
"""方向別 dt の上限の A/B の判定 (plan time_integration-line-viscous-jacobian-dt §6、2026-10-11 事前登録)。
faceh の 7/5 の判定器と同じ分類・ゲートで、腕の違いを方向別 dt の上限 (A = 上限なし、B = 上限 50、C = 方向別なし) だけにしたもの。
バイナリは全部段 ③ の FP64 (~/forge-fgeom3-fp64)、値 3・マスク 7 (FORGE_LVC_TERMS=7)。

lvcdt.sh の run_0580〜0583 (A・B の腕)、あれば run_0588・0589 (条件付きの C)、事前のゲートの記録 (lvcdt_pregate.json、run_0584〜0587 の 1 step の書き出しから) を読み、
  (1) 入力のゲート (バイナリ・出発の場・格子の実体・入力ファイル・腕ごとの期待どおりの設定・起動ログ・マスクの表示)、
  (2) 証拠のゲート (初期場・破綻前の場・res_nan・帳簿・FINITE 側の区間付き収束 VERDICT)、
  (3) 各 run の分類 (DIVERGED / FINITE / INVALID)
  (4) 事前のゲートが PASS で、記録した証拠のハッシュと腕の設定 (同じ水準の書き出しと、step 数・出力の間隔・extraFields 以外で同じ) が合う
  (5) 記録だけ: 序盤 200 呼び出しの帳簿の第一内部節点の |Δρ/ρ| の最大、各残差の列が出発の 10 倍に達した step
を行い、§6 の分岐どおりに判定する (主判定は「上限 50 で 2000 step 以内の非有限化を回避したか」に限る。plan-dt レビュー M3)。
--ab-only: C を除いて A/B だけを判定し、_band_ab/cold_pair/lvcdt_judge_ab.json に書く。条件付きの C の起動の直前に台本が呼ぶ
(ゲートが全部合格・INVALID なし・主判定が「回避せず」のときだけ終了コード 0。plan-dt レビュー M5)。出力は _band_ab/cold_pair/lvcdt_judge.json (異常でも必ず書く)。
ゲートの不合格か INVALID の run があれば判定せず、終了コード 1。規則を変えるときは plan §6 を先に改訂する。
"""
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import traceback
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
SHA_NEW = "129de3f4e7f67aa3a80dbd30d5f5cb75df1582d3974e702d76b61c8998598cec"   # 段 ③ の FP64 (~/forge-fgeom3-fp64)
# (run, 面エンタルピーの切替 (常に 0)、期待のバイナリ)
RUNS = {"a1": ("run_0580_lvcdt_a1", 0, SHA_NEW), "b1": ("run_0581_lvcdt_b1", 0, SHA_NEW),
        "a2": ("run_0582_lvcdt_a2", 0, SHA_NEW), "b2": ("run_0583_lvcdt_b2", 0, SHA_NEW)}
RUNS_C = {"c1": ("run_0588_lvcdt_c1", 0, SHA_NEW), "c2": ("run_0589_lvcdt_c2", 0, SHA_NEW)}   # 条件付き (分岐 2 のときだけ回す)
AB_ONLY = sys.argv[1:] == ["--ab-only"]
if not AB_ONLY and any((Path(__file__).resolve().parent / r).exists() for r, _, _ in RUNS_C.values()):
    RUNS.update(RUNS_C)                       # C を回したら 2 本とも判定に入れる (片方だけなら INVALID になる)
ARM = {"a1": "a", "a2": "a", "b1": "b", "b2": "b", "c1": "c", "c2": "c"}
ARM_CFG = {"a": {}, "b": {"time/deltaT/lineDtDirectionalCap": 50.0}, "c": {"time/deltaT/lineDtDirectional": None}}
DT_KEYS = {"time/deltaT/lineDtDirectional", "time/deltaT/lineDtDirectionalCap"}
MASK = {k: 7 for k in ARM}                    # FORGE_LVC_TERMS (台本が run ごとに LVC_TERMS.txt に書く)
LEDGER_FIRST = [1571, 4112, 7984, 4959, 6169]   # 記録だけ: 帳簿の第一内部節点 (前の 3 つは書き出しの r の節点と同じ列、後の 2 つは §6.10 で壊れた列 35〜52 をはさむ列 40・50)
DUMPS = ["run_0584_lvcdt_a_dump", "run_0585_lvcdt_b_dump", "run_0586_lvcdt_a_dump2", "run_0587_lvcdt_c_dump"]   # 事前のゲートが読む 1 step の書き出し
SRC, SRC_RES = "run_0183_ns_coldmesh_tw300_ext", "res_100000.h5"
SRC_SHA16 = "207d39f0e7f4aa03"
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv",
              "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
REQ_COLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
CONS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1", "P", "T"]
NSTEP, OUT, LEDGER_CALLS = 2000, 100, 200
LEDGER_NODES = [int(x) for x in (
    "1572,1571,1570,1569,1568,1567,1566,1565,1564,1563,1562,1561,1560,1559,1558,1557,4113,4112,4111,4110,4109,4108,4107,4106,"
    "4105,4104,4103,4102,4101,4100,4099,4098,4960,4959,4958,4957,4956,4955,4954,4953,4952,4951,4950,4949,4948,4947,4946,4945,"
    "6170,6169,6168,6167,6166,6165,6164,6163,6162,6161,6160,6159,6158,6157,6156,6155,7985,7984,7983,7982,7981,7980,7979,7978,"
    "7977,7976,7975,7974,7973,7972,7971,7970,198560,198559,198558,198557,198556,198555,198554,198553,198552,198551,198550,"
    "198549,198548,198547,198546,198545,264263,264262,264261,264260,264259,264258,264257,264256,264255,264254,264253,264252,"
    "264251,264250,264249,264248").split(",")]
# §6 の期待の設定 (平坦化したキー → 値)。None は「書いていない」ことを要求する
EXPECT = {"time/deltaT/cfl": 4.0, "time/deltaT/cfl_pseudo": 4.0, "time/deltaT/implicitRelax": 0.7,
          "time/deltaT/implicitThermalJacobian": 5, "time/deltaT/lineDtDirectional": 1, "time/deltaT/lineImplicit": 1,
          "time/deltaT/lineViscCoupling": 3, "time/deltaT/blockDPLUR": 1, "time/deltaT/detectNaN": 1,
          "time/timeIntegration": 11, "time/nStepInner": 5, "time/last/nStepOuter": NSTEP, "time/outStepInterval": OUT,
          "time/deltaT/lineDtDirectionalCap": None, "time/deltaT/implicitSolvePrecision": None,
          "output/extraFields": ["res_ro", "volume"],
          "time/deltaT/axisTimestepBeta": None, "time/deltaT/lineViscousDtRelief": None}
DELAY = 2.0      # 「遅らせる」の記録の倍率 (判定には使わない)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def meshsha(p):
    h = hashlib.sha256()

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            h.update(name.encode())
            h.update(obj[...].tobytes())
    with h5py.File(p, "r") as f:
        f["MESH"].visititems(visit)
    return h.hexdigest()


def flat(d, p=""):
    out = {}
    for k, v in d.items():
        kk = f"{p}/{k}" if p else str(k)
        out.update(flat(v, kk) if isinstance(v, dict) else {kk: v})
    return out


def cfg(run):
    return flat(yaml.safe_load((HERE / run / "solverConfig.yaml").read_text()))


def read_csv(run):
    """残差の CSV を厳密に読む。問題は reason に書く (読めた範囲の値も返す)。"""
    p = HERE / run / "residual_history.csv"
    if not p.is_file():
        return None, "残差の CSV が無い"
    with open(p, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return None, "残差の CSV が空"
    hdr = rows[0]
    need = ["step", "inner", "phase"] + REQ_COLS
    if len(set(hdr)) != len(hdr):
        return None, "ヘッダーに重複がある"
    miss = [c for c in need if c not in hdr]
    if miss:
        return None, f"必須の列が無い {miss}"
    ix = {c: hdr.index(c) for c in hdr}
    cols = [c for c in hdr if c.startswith("rms_") and not c.startswith("rms_dq")]
    outer, dup, broken, nonfin, neg = {}, [], [], None, None
    for n, r in enumerate(rows[1:], start=2):
        if len(r) != len(hdr):
            broken.append(n)
            continue
        try:
            s = int(r[ix["step"]])
            v = {c: float(r[ix[c]]) for c in cols}
        except ValueError:
            broken.append(n)
            continue
        if nonfin is None and any(not math.isfinite(x) for x in v.values()):
            nonfin = s
        if neg is None and any(math.isfinite(x) and x < 0 for x in v.values()):
            neg = s
        if r[ix["phase"]] == "outer_begin":
            if s in outer:
                dup.append(s)
            outer[s] = v
    return {"outer": outer, "dup": dup, "broken": broken, "first_nonfinite": nonfin, "first_negative": neg, "cols": cols}, None


def field_ok(p, require_finite=True):
    """場を読めるか、保存量・P・T が有限で ρ・P・T が正か。問題を文字列で返す (問題なしなら None)。"""
    try:
        with h5py.File(p, "r") as h:
            bad = {}
            for v in CONS:
                a = h["VALUE"][v][...]
                if require_finite:
                    n = int(np.sum(~np.isfinite(a))) + (int(np.sum(a <= 0)) if v in ("ro", "P", "T") else 0)
                    if n:
                        bad[v] = n
    except (OSError, KeyError) as e:
        return f"{p.name} が読めない ({e})"
    return f"{p.name} に非有限・非正 {bad}" if bad else None


def ledger_check(run, last_call):
    """帳簿のヘッダー・節点の集合・呼び出しの連続 (1〜last_call) を確かめる。問題を返す。"""
    p = HERE / run / "ledger.csv"
    if not p.is_file():
        return "帳簿が無い"
    nodes, calls = set(), set()
    with open(p, newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd, None)
        if hdr != ["call", "tag", "node", "field", "value"]:
            return f"帳簿のヘッダーが違う {hdr}"
        for r in rd:
            if len(r) != 5:
                return "帳簿に列数の違う行がある"
            calls.add(int(r[0]))
            nodes.add(int(r[2]))
    if nodes != set(LEDGER_NODES):
        return f"帳簿の節点が違う (記録 {len(nodes)}、期待 {len(LEDGER_NODES)})"
    if calls != set(range(1, last_call + 1)):
        return f"帳簿の呼び出しが 1〜{last_call} でない (記録 {min(calls, default=None)}〜{max(calls, default=None)}、{len(calls)} 回)"
    return None


def assess(run, fh):
    """1 本の run の分類と証拠の確認。例外は INVALID に変える。"""
    d = HERE / run
    info = {"run": run}
    try:
        log = (d / "forge_run.log").read_text(errors="replace")
        m = re.search(r"Non-finite value detected in '([A-Za-z0-9_]+)' at step (\d+)", log)
        rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").is_file() else None
        info.update(run_rc=rc, detectNaN=[m.group(1), int(m.group(2))] if m else None)
        c, why = read_csv(run)
        if c is None:
            return "INVALID", {**info, "why": why}
        o = c["outer"]
        info.update(last_step=max(o) if o else None, csv_first_nonfinite=c["first_nonfinite"])
        if c["broken"] or c["dup"] or c["first_negative"] is not None:
            return "INVALID", {**info, "why": f"読めない行 {c['broken'][:3]}・重複 {c['dup'][:3]}・負の値 (step {c['first_negative']})"}
        if sorted(o) != list(range(len(o))):
            return "INVALID", {**info, "why": "outer_begin の行が 0 から連続していない"}
        ev = []                                   # 証拠のゲート
        q = field_ok(d / "res_0.h5")
        if q:
            ev.append(f"初期場: {q}")
        if not (d / "nozzle.h5").is_file():
            ev.append("格子 (nozzle.h5) が無い")
        if m or c["first_nonfinite"] is not None:
            steps = [s for s in (int(m.group(2)) if m else None, c["first_nonfinite"]) if s is not None]
            n = min(steps)
            info.update(diverged_step=n, steps_log_csv=[int(m.group(2)) if m else None, c["first_nonfinite"]])
            for k in range(OUT, n, OUT):          # 破綻より前に保存されるはずの場
                q = field_ok(d / f"res_{k}.h5")
                if q:
                    ev.append(f"破綻前の場: {q}")
            if m:
                q = field_ok(d / f"res_nan_{int(m.group(2))}.h5", require_finite=False)
                if q:
                    ev.append(f"res_nan: {q}")
            q = ledger_check(run, min(LEDGER_CALLS, int(m.group(2)) if m else n))
            if q:
                ev.append(q)
            info["evidence_problems"] = ev
            return ("INVALID", {**info, "why": "証拠の欠け"}) if ev else ("DIVERGED", info)
        if sorted(o) != list(range(NSTEP)) or rc != "0":
            return "INVALID", {**info, "why": f"2000 step に届いていないのに非有限が無い (行 {len(o)}、RUN_RC {rc})"}
        for k in range(OUT, NSTEP + 1, OUT):
            q = field_ok(d / f"res_{k}.h5")
            if q:
                ev.append(q)
        q = ledger_check(run, LEDGER_CALLS)
        if q:
            ev.append(q)
        seg = d / "CONVERGENCE_SEGMENT.txt"
        st = seg.read_text(errors="replace") if seg.is_file() else ""
        mm = re.search(r"->\s*([A-Z][A-Z ]+[A-Z])", st)
        if not (mm and "判定区間" in st):
            ev.append("区間付きの収束 VERDICT が無い")
        else:
            info["segment_verdict"] = mm.group(1)
        info["evidence_problems"] = ev
        return ("INVALID", {**info, "why": "証拠の欠け"}) if ev else ("FINITE", info)
    except Exception as e:                        # 想定外の入力も理由付きの INVALID にする
        return "INVALID", {**info, "why": f"例外 {type(e).__name__}: {e}", "trace": traceback.format_exc(limit=3)}


def records(run):
    c, _ = read_csv(run)
    if c is None:
        return {}
    o = c["outer"]
    steps = sorted(s for s in o if all(math.isfinite(x) for x in o[s].values()))
    rec = {}
    for col in c["cols"]:
        s0 = o[steps[0]][col] if steps else float("nan")
        if not (s0 > 0):
            continue
        vals = [o[s][col] for s in steps]
        rec[col] = {"max_over_start": max(vals) / s0,
                    "first_step_over_3x": next((s for s in steps if o[s][col] > 3 * s0), None),
                    "first_step_over_10x": next((s for s in steps if o[s][col] > 10 * s0), None)}
        if len(steps) == NSTEP:
            x = np.arange(NSTEP - 500, NSTEP, dtype=float)
            rec[col]["tail_slope500"] = float(np.polyfit(x, np.log10([o[int(k)][col] for k in x]), 1)[0] * 500)
    return rec


def ledger_drho(run, last_call):
    """帳簿の第一内部節点の |ρ(call) − ρ(call 1 の入口)| / ρ(call 1 の入口) の最大 (call 1〜last_call、tag after_eos_bc)。
    対象の (call, node, tag) が欠け・重複・非有限・基準 ≤ 0 なら ValueError (判定を INVALID にする。plan-dt レビュー M1)"""
    want = {(1, "entry", n) for n in LEDGER_FIRST} | {(c, "after_eos_bc", n) for c in range(1, last_call + 1) for n in LEDGER_FIRST}
    got = {}
    with open(HERE / run / "ledger.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["field"] != "ro":
                continue
            key = (int(r["call"]), r["tag"], int(r["node"]))
            if key in want:
                if key in got:
                    raise ValueError(f"{run}: 帳簿の行が重複 {key}")
                got[key] = float(r["value"])
    miss = sorted(want - set(got))
    if miss:
        raise ValueError(f"{run}: 帳簿に必要な行が無い ({len(miss)} 件、最初 {miss[0]})")
    if not all(math.isfinite(v) for v in got.values()):
        raise ValueError(f"{run}: 帳簿の値に非有限")
    base = {n: got[(1, "entry", n)] for n in LEDGER_FIRST}
    if not all(v > 0 for v in base.values()):
        raise ValueError(f"{run}: 帳簿の基準の ρ が正でない")
    per = {n: max(abs(got[(c, "after_eos_bc", n)] - base[n]) / base[n] for c in range(1, last_call + 1)) for n in LEDGER_FIRST}
    return {"max": max(per.values()), "per_node": per, "calls": last_call}


def decide(cls):
    A = {k: cls[k] for k in ("a1", "a2")}
    B = {k: cls[k] for k in ("b1", "b2")}
    C = {k: cls[k] for k in ("c1", "c2") if k in cls}
    ca, cb, cc = [v[0] for v in A.values()], [v[0] for v in B.values()], [v[0] for v in C.values()]
    if any(x == "INVALID" for x in ca + cb + cc):
        return "INVALID", "INVALID の run がある", {}
    # 記録だけ: 帳簿の第一層の |Δρ/ρ| (窓は 1〜200 呼び出し。破綻した run は 1〜min(200, 破綻の step − 2)。plan §6 に登録) — 欠けは INVALID
    led = {}
    for k in list(A) + list(B) + list(C):
        ds = cls[k][1].get("diverged_step")
        led[k] = ledger_drho(RUNS[k][0], min(LEDGER_CALLS, ds - 2) if ds is not None else LEDGER_CALLS)
    if any(x == "FINITE" for x in ca):
        return "帰属不能", "A (上限なし) が既知の破綻を再現しない", led
    if all(x == "FINITE" for x in cb):
        return "回避", "A は全部 DIVERGED、B (上限 50) は全部 FINITE", led
    if all(x == "DIVERGED" for x in cb):
        sb = [v[1]["diverged_step"] for v in B.values()]
        if not C:
            return "回避せず", f"B (上限 50) でも全部 DIVERGED (破綻 step {sb})。C (方向別なし) を 2 本回す", led
        if len(C) != 2:
            return "INVALID", "C が 2 本そろっていない", led
        if all(x == "FINITE" for x in cc):
            return "回避せず・C 有限", f"B でも全部 DIVERGED (step {sb})、C (point の dt) は全部 FINITE", led
        if all(x == "DIVERGED" for x in cc):
            return "回避せず・C 発散", f"B でも全部 DIVERGED (step {sb})、C も全部 DIVERGED", led
        return "判別不能", "C の run が分かれた", led
    return "判別不能", "B の run が分かれた", led


def input_gates():
    gates = []
    src_full = sha(HERE / SRC / SRC_RES)
    gates.append(("出発の場の sha256 が事前に固定した値", src_full.startswith(SRC_SHA16), src_full[:16]))
    ref_mesh = meshsha(HERE / SRC / "nozzle.h5")
    first = RUNS["a1"][0]
    for name, (run, fh, binsha) in RUNS.items():
        d = HERE / run
        try:
            prov = (d / "RUN_PROVENANCE.txt").read_text(errors="replace")
            gates.append((f"{run}: forge の sha256 が腕のバイナリ", binsha in prov, binsha[:16]))
            cp = json.loads((d / "COLD_PAIR.json").read_text())
            gates.append((f"{run}: 出発の場が {SRC}/{SRC_RES}", cp["parent"] == SRC and cp["parent_res"] == SRC_RES and cp["parent_res_sha256"] == src_full,
                          [cp["parent"], cp["parent_res"], cp["parent_res_sha256"][:16]]))
            m = meshsha(d / "nozzle.h5")
            gates.append((f"{run}: 格子の実体が出発 run と同じ", m == ref_mesh and (d / "MESH_SHA.txt").read_text().strip() == m, m[:16]))
            diff_in = [fn for fn in SAME_FILES if not (d / fn).is_file() or not (HERE / SRC / fn).is_file() or sha(HERE / SRC / fn) != sha(d / fn)]
            gates.append((f"{run}: 入力ファイルが出発 run と同じ", diff_in == [], diff_in))
            mq = (d / "MESH_QUALITY.txt").read_text(errors="replace")
            gates.append((f"{run}: メッシュ品質の VERDICT が PASS", "VERDICT: PASS" in mq, None))
            c = cfg(run)
            exp_ = {**EXPECT, **ARM_CFG[ARM[name]]}
            bad = {k: c.get(k) for k, v in exp_.items() if (k in c if v is None else c.get(k) != v)}
            gates.append((f"{run}: 設定が §6 の腕 {ARM[name].upper()} の期待どおり", bad == {}, bad))
            dd = sorted(k for k in set(c) | set(cfg(first)) if k not in DT_KEYS and c.get(k) != cfg(first).get(k))
            gates.append((f"{run}: 設定が a1 と dt の 2 キー以外で一致", dd == [], dd))
            log = (d / "forge_run.log").read_text(errors="replace")
            gates.append((f"{run}: 起動ログの lineViscCoupling = 3", "'lineViscCoupling' in 'time.deltaT': 3" in log, None))
            gates.append((f"{run}: 切替の表示 = {fh}", ("FORGE_DIAG_FACE_H_DOUBLE" in log) == bool(fh), None))
            mk = MASK[name]
            shown = ("FORGE_LVC_TERMS=" not in log) if mk == 7 else (f"診断のマスク FORGE_LVC_TERMS={mk} " in log)
            gates.append((f"{run}: マスクの表示が {mk} の腕と合う (7 は表示なし)", shown, None))
            lt = (d / "LVC_TERMS.txt").read_text().strip() if (d / "LVC_TERMS.txt").is_file() else None
            gates.append((f"{run}: 台本が渡した FORGE_LVC_TERMS = {mk}", lt == str(mk), lt))
            gates.append((f"{run}: 実効の並びが LAYOUT2", "Thomas の配列の並び: LAYOUT2" in log, None))
        except Exception as e:
            gates.append((f"{run}: 入力のゲートを評価できない", False, f"{type(e).__name__}: {e}"))
    return gates


def intervention(rec):
    """介入の成立と残差の不変は事前のゲート (lvcdt_pregate.py、plan §6) が判定する。ここではその記録が
    PASS で、記録した証拠のハッシュが今のファイルと同じで、腕の設定が同じ水準の書き出しと合うことを確かめ、要点を写す。"""
    g = []
    try:
        vr = subprocess.run([sys.executable, str(HERE / "lvcdt_pregate.py"), "--verify"], capture_output=True, text=True)
        g.append(("事前のゲートの記録が PASS で、記録した証拠のハッシュ (書き出し・監査の記録・ゲート自身) が今のファイルと同じ", vr.returncode == 0,
                  (vr.stdout + vr.stderr).strip()[-300:]))
        pg = json.loads((HERE / "_band_ab" / "cold_pair" / "lvcdt_pregate.json").read_text())
        # 本試験 4 本の設定が、事前のゲートの書き出しの設定と、登録した差 (step 数・出力の間隔・extraFields) 以外で全項目同じ (plan-6 レビュー M3)
        dcs = pg.get("dump_config") or {}
        allowed = {"time/last/nStepOuter", "time/outStepInterval", "output/extraFields"}
        for name, (run, fh, binsha) in RUNS.items():
            c = cfg(run)
            dc = dcs.get(ARM[name]) or {}
            dd = sorted(k for k in set(c) | set(dc) if k not in allowed and c.get(k) != dc.get(k))
            g.append((f"{run}: 設定が事前のゲートの書き出しと同じ (差は step 数・出力の間隔・extraFields だけ)", bool(dc) and dd == [], dd[:10]))
        rec["pregate"] = {k: pg.get(k) for k in ("VERDICT", "residual", "struct_rhs", "intervention")}
        g.append(("事前のゲートの VERDICT が PASS", str(pg.get("VERDICT", "")).startswith("PASS"), pg.get("VERDICT")))
        g.append(("事前のゲートの項目がすべて合格", bool(pg.get("checks")) and all(c["ok"] for c in pg["checks"]), sum(not c["ok"] for c in pg.get("checks", []))))
        for run in DUMPS:
            g.append((f"{run}: 書き出しが残っている", (HERE / run / "linedump" / "D.f64").is_file(), None))
    except Exception as e:
        g.append(("事前のゲートの記録を読めない", False, f"{type(e).__name__}: {e}"))
    return g


def main():
    rec = {"plan": "time_integration-line-viscous-jacobian-dt §6", "ab_only": AB_ONLY}
    out = HERE / "_band_ab" / "cold_pair" / ("lvcdt_judge_ab.json" if AB_ONLY else "lvcdt_judge.json")
    try:
        gates = input_gates()
    except Exception as e:
        gates = [("入力のゲートを評価できない", False, f"{type(e).__name__}: {e}")]
    gates += intervention(rec)
    cls = {name: assess(run, fh) for name, (run, fh, binsha) in RUNS.items()}
    rec["gates"] = [{"check": g[0], "ok": bool(g[1]), "value": g[2]} for g in gates]
    rec["gates_ok"] = all(g[1] for g in gates)
    rec["runs"] = {}
    for name, (run, fh, binsha) in RUNS.items():
        r = {"class": cls[name][0], **cls[name][1]}
        if cls[name][0] != "INVALID":
            r["records"] = records(run)
        rec["runs"][name] = r
    invalid = [n for n in RUNS if cls[n][0] == "INVALID"]
    if not rec["gates_ok"] or invalid:
        rec["VERDICT"] = f"INVALID (ゲート不合格 {sum(not g[1] for g in gates)} 件、INVALID の run {invalid}): 判定しない"
    else:
        try:
            v3 = decide(cls)
        except Exception as e:
            v3 = ("INVALID", f"判定で例外 {type(e).__name__}: {e}", {})
        txt = {"回避": "この条件・2000 step では、上限 50 で全部入りの非有限化を避けた (全部入りが成り立った・収束した・H1 の機構が正しいとは言わない)",
               "回避せず": "上限 50 でも全部入りの非有限化を避けられない。C を回す",
               "回避せず・C 有限": "point の dt では有限 (H1 の支持にはしない)。次の調査の順は H3 (反復の写像をスケーリングして見る) から",
               "回避せず・C 発散": "point の dt でも非有限。次の調査の順は H2 (実残差との整合) から"}.get(v3[0], v3[0])
        rec["main_v3"] = {"verdict": v3[0], "detail": v3[1], "ledger_first_layer_drho_record_only": v3[2]}
        rec["VERDICT"] = f"主 (値 3・マスク 7): {v3[0]} — {txt}; {v3[1]}"
        if v3[0] == "INVALID":
            invalid = invalid or ["decide"]
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for g in gates:
        if not g[1]:
            print("ゲート不合格:", g[0], g[2])
    for name in RUNS:
        r = rec["runs"][name]
        print(f"{name} {RUNS[name][0]}: {r['class']} step={r.get('diverged_step', r.get('last_step'))} "
              f"log/csv={r.get('steps_log_csv')} nan={r.get('detectNaN')} {r.get('why', '')} {r.get('evidence_problems', '')}")
    print(rec["VERDICT"])
    if AB_ONLY:                                 # C の起動の条件: ゲート合格・INVALID なし・主判定が「回避せず」
        return 0 if rec["gates_ok"] and not invalid and rec.get("main_v3", {}).get("verdict") == "回避せず" else 1
    return 0 if rec["gates_ok"] and not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
