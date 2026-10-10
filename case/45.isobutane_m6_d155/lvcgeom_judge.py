#!/usr/bin/env python3
"""値 3・マスク 7 の破綻と LHS の座標の差 (float の座標の差 / double の座標の差) の A/B の判定
(plan time_integration-line-viscous-jacobian-faceh §6.9、2026-10-10 事前登録)。§6.2 の lvcfh_judge.py と同じ分類・ゲートに、腕ごとのバイナリと介入の成立の確認を足したもの。
A = 元のセッションの float 化の段 ② の FP64 (LHS は ST(ccx[o]) − ST(ccx[ic]))、B = 段 ③ の FP64 (LHS は ST(±ge_x))。

lvcgeom.sh の run_0560〜0563 (腕) と事前のゲートの記録 (lvcgeom_pregate.json、run_0564〜0566 の 1 step の書き出しから) を読み、
  (1) 入力のゲート (バイナリ・出発の場・格子の実体・入力ファイル・§6 の期待どおりの設定・起動ログ)、
  (2) 証拠のゲート (初期場・破綻前の場・res_nan・帳簿・FINITE 側の区間付き収束 VERDICT)、
  (3) 各 run の分類 (DIVERGED / FINITE / INVALID)
  (4) 事前のゲートが PASS (入力の一致・残差の不変・係数の水準での介入の成立。本判定はその記録を写すだけ)
を行い、§6.9 の分岐どおりに判定する。出力は _band_ab/cold_pair/lvcgeom_judge.json (異常でも必ず書く)。
ゲートの不合格か INVALID の run があれば判定せず、終了コード 1。規則を変えるときは plan §6.9 を先に改訂する。
"""
import csv
import hashlib
import json
import math
import re
import sys
import traceback
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
SHA_OLD = "1b8590e84ece63ae50a05cd2d6934440be5f90cbf7f7723a40cd6776ee7285a1"   # 段 ② の FP64 (~/forge-fgeom2-fp64)
SHA_NEW = "129de3f4e7f67aa3a80dbd30d5f5cb75df1582d3974e702d76b61c8998598cec"   # 段 ③ の FP64 (~/forge-fgeom3-fp64)
# (run, 面エンタルピーの切替 (常に 0)、期待のバイナリ)
RUNS = {"a1": ("run_0560_lvcgeom_old_a1", 0, SHA_OLD), "b1": ("run_0561_lvcgeom_new_b1", 0, SHA_NEW),
        "a2": ("run_0562_lvcgeom_old_a2", 0, SHA_OLD), "b2": ("run_0563_lvcgeom_new_b2", 0, SHA_NEW)}
DUMPS = ["run_0564_lvcgeom_old_dump", "run_0565_lvcgeom_new_dump", "run_0566_lvcgeom_old_dump2"]   # 事前のゲートが読む 1 step の書き出し
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
          "output/extraFields": ["res_ro", "volume"]}
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
                    "first_step_over_3x": next((s for s in steps if o[s][col] > 3 * s0), None)}
        if len(steps) == NSTEP:
            x = np.arange(NSTEP - 500, NSTEP, dtype=float)
            rec[col]["tail_slope500"] = float(np.polyfit(x, np.log10([o[int(k)][col] for k in x]), 1)[0] * 500)
    return rec


def pair(A, B):
    ca, cb = [A[k][0] for k in A], [B[k][0] for k in B]
    if any(x == "INVALID" for x in ca + cb):
        return "INVALID", "INVALID の run がある"
    if any(x == "FINITE" for x in ca):
        return "帰属不能", "旧 (段 ②) が既知の破綻を再現しない"
    if all(x == "FINITE" for x in cb):
        return "支持", "旧 (段 ②) は全部 DIVERGED、新 (段 ③) は全部 FINITE"
    if all(x == "DIVERGED" for x in cb):
        sa = max(A[k][1]["diverged_step"] for k in A)
        sb = [B[k][1]["diverged_step"] for k in B]
        note = "遅らせる (記録だけ)" if all(s > DELAY * sa for s in sb) else ""
        return "棄却", f"新でも全部 DIVERGED。新の破綻 step {sb}、旧の最大 {sa} {note}".strip()
    return "判別不能", "新の run が分かれた"


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
            bad = {k: c.get(k) for k, v in EXPECT.items() if (k in c if v is None else c.get(k) != v)}
            gates.append((f"{run}: 設定が §6 の期待どおり", bad == {}, bad))
            dd = sorted(k for k in set(c) | set(cfg(first)) if c.get(k) != cfg(first).get(k))
            gates.append((f"{run}: 設定が a1 と一致", dd == [], dd))
            log = (d / "forge_run.log").read_text(errors="replace")
            gates.append((f"{run}: 起動ログの lineViscCoupling = 3", "'lineViscCoupling' in 'time.deltaT': 3" in log, None))
            gates.append((f"{run}: 切替の表示 = {fh}", ("FORGE_DIAG_FACE_H_DOUBLE" in log) == bool(fh), None))
            gates.append((f"{run}: マスクの表示が無い (= 既定 7)", "FORGE_LVC_TERMS=" not in log, None))
            gates.append((f"{run}: 実効の並びが LAYOUT2", "Thomas の配列の並び: LAYOUT2" in log, None))
        except Exception as e:
            gates.append((f"{run}: 入力のゲートを評価できない", False, f"{type(e).__name__}: {e}"))
    return gates


def intervention(rec):
    """介入の成立と残差の不変は事前のゲート (lvcgeom_pregate.py、plan §6.9 の M1・M2・M3) が判定する。ここではその記録が
    PASS で、ゲートが読んだ書き出し (run_0564・0565・0566) が今もあることだけを確かめ、要点を写す。"""
    g = []
    try:
        pg = json.loads((HERE / "_band_ab" / "cold_pair" / "lvcgeom_pregate.json").read_text())
        rec["pregate"] = {k: pg.get(k) for k in ("VERDICT", "residual", "intervention_new_coef_max_rel_vs_double", "recompute_old_formula_max_ulp",
                                                 "intervention_K", "intervention_D", "near_wall_K_change_max_rel")}
        g.append(("事前のゲートの VERDICT が PASS", str(pg.get("VERDICT", "")).startswith("PASS"), pg.get("VERDICT")))
        g.append(("事前のゲートの項目がすべて合格", bool(pg.get("checks")) and all(c["ok"] for c in pg["checks"]), sum(not c["ok"] for c in pg.get("checks", []))))
        for run in DUMPS:
            g.append((f"{run}: 書き出しが残っている", (HERE / run / "linedump" / "D.f64").is_file(), None))
    except Exception as e:
        g.append(("事前のゲートの記録を読めない", False, f"{type(e).__name__}: {e}"))
    return g


def main():
    rec = {"plan": "time_integration-line-viscous-jacobian-faceh §6.9"}
    out = HERE / "_band_ab" / "cold_pair" / "lvcgeom_judge.json"
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
        v3 = pair({k: cls[k] for k in ("a1", "a2")}, {k: cls[k] for k in ("b1", "b2")})
        txt = {"支持": "この条件・期間では、LHS の座標の差の丸めを直すだけで非有限化を回避できることを支持 (収束・安定・速さ・B0 への影響は言わない)",
               "棄却": "この修復だけで十分という仮説を棄却 (精度の不整合の寄与をゼロとは言わない)"}.get(v3[0], v3[0])
        rec["main_v3"] = {"verdict": v3[0], "detail": v3[1]}
        rec["VERDICT"] = f"主 (値 3): {v3[0]} — {txt}; {v3[1]}"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for g in gates:
        if not g[1]:
            print("ゲート不合格:", g[0], g[2])
    for name in RUNS:
        r = rec["runs"][name]
        print(f"{name} {RUNS[name][0]}: {r['class']} step={r.get('diverged_step', r.get('last_step'))} "
              f"log/csv={r.get('steps_log_csv')} nan={r.get('detectNaN')} {r.get('why', '')} {r.get('evidence_problems', '')}")
    print(rec["VERDICT"])
    return 0 if rec["gates_ok"] and not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
