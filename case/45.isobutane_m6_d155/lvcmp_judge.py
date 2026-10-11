#!/usr/bin/env python3
"""上限 50 での全部入り・粘性入り・生産の比較の判定 (plan time_integration-line-viscous-jacobian-dt §6.3、2026-10-11 事前登録)。

lvcmp.sh の run_0553 (L7 = 値 3・マスク 7)・run_0554 (L0 = 値 0、生産)・run_0555 (L5 = 値 3・マスク 5) と、事前のゲートの記録
(lvcmp_pregate.json、run_0556〜0559 の 1 step の書き出しから) を読み、
  (1) 入力のゲート: バイナリ・出発の場・格子の実体・入力ファイル・腕ごとの期待どおりの設定・起動ログ (値・マスク・LAYOUT2)、
      事前のゲートの --verify、腕の設定が同じ腕の書き出しの設定と step 数・出力の間隔・extraFields 以外で同じ
  (2) 腕ごと (codex plan-2 M1・M2・m5): 見張りの記録のスキーマ (run・phase line・budget 150000・consec 1) と、2500 step ごとの系列が終端まで欠け・重複なく
      そろうこと、記録係 (lvcmp_outputs.json) が系列の全出力について sha256・全保存量の有限・ρ/P/T > 0 を記録していること。
      REACHED: 到達の step を系列から水準の規則で計算し直して一致、到達の出力の sha256 が見張りの記録と一致し、その出力から量を計算し直して系列の行と一致。
      CENSORED: RUN_RC 0・最後の出力が 150000・系列は水準に入っていない。DIVERGED: 非有限の出力の行か、res_nan と detectNaN のログ。
      区間付きの収束: 判定不能を含む・必須の残差の列が無い・VERDICT が無い・DIVERGED は INVALID (REACHED・CENSORED)。EXEC_ERROR・DATA_ERROR・記録なしは INVALID
  (3) 比較: 2 本とも REACHED のとき、到達の step の差が 15000 以上なら「今回の実行では速い/遅い」、未満は「差を検出しなかった」。どちらかが届かなければ比較しない。
      各腕 1 本なので探索的な比較で、再実行の揺れを除いた差とは言わない (codex plan-2 M3)
  (4) 記録だけ: 到達の窓 [N − 2 万, N] の θ_r・Q_w の check_quasisteady (--tail 1.0 --drift 0.0005 --min-snaps 5)、各列の区間付きの収束
を行い、§6.3 の分岐どおりに判定する。出力は _band_ab/cold_pair/lvcmp_judge.json (異常でも必ず書く)。ゲートの不合格か INVALID があれば判定せず、終了コード 1。
"""
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
D = HERE / "_band_ab" / "cold_pair"
D.mkdir(parents=True, exist_ok=True)
SHA_NEW = "129de3f4e7f67aa3a80dbd30d5f5cb75df1582d3974e702d76b61c8998598cec"
SRC, SRC_RES, SRC_SHA16 = "run_0183_ns_coldmesh_tw300_ext", "res_100000.h5", "207d39f0e7f4aa03"
ARMS = {"L7": ("run_0553_lvcmp_l7", "7", "m7"), "L0": ("run_0554_lvcmp_l0", "-", "m0"), "L5": ("run_0555_lvcmp_l5", "5", "m5")}
DUMPS = ["run_0556_lvcmp_l7_dump", "run_0557_lvcmp_l0_dump", "run_0558_lvcmp_l5_dump", "run_0559_lvcmp_l7_dump2"]
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv", "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
BUDGET, OUT, THR = 150000, 2500, 15000
EXPECT = {"time/deltaT/cfl": 4.0, "time/deltaT/cfl_pseudo": 4.0, "time/deltaT/implicitRelax": 0.7, "time/deltaT/implicitThermalJacobian": 5,
          "time/deltaT/lineDtDirectional": 1, "time/deltaT/lineImplicit": 1, "time/deltaT/lineDtDirectionalCap": 50.0, "time/deltaT/blockDPLUR": 1,
          "time/deltaT/detectNaN": 1, "time/timeIntegration": 11, "time/nStepInner": 5, "time/last/nStepOuter": BUDGET, "time/outStepInterval": OUT,
          "time/deltaT/implicitSolvePrecision": None, "time/deltaT/axisTimestepBeta": None, "time/deltaT/lineViscousDtRelief": None,
          "output/extraFields": ["res_ro"]}
ARM_CFG = {"L7": {"time/deltaT/lineViscCoupling": 3}, "L5": {"time/deltaT/lineViscCoupling": 3}, "L0": {"time/deltaT/lineViscCoupling": None}}
KIND = {"L7": "L7", "L5": "L5", "L0": "L0"}   # 腕の名前 → 構成の種類 (§6.5 の再実行 L7r は L7 の構成。codex plan-3 M1)
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def meshsha(p):
    h = hashlib.sha256()

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            h.update(name.encode()); h.update(obj[...].tobytes())
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


def drift(rows, key):
    s_ = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000:
        return None
    j = int(np.searchsorted(s_, s_[-1] - 20000))
    return 100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j])


def reach_level(rows):
    """水準 (速度 plan §6.15、窓は run の中・step 0 は使わない) に初めて入った出力。非有限なら diverged_at"""
    rows = sorted([x for x in rows if x["step"] > 0], key=lambda x: x["step"])
    for n in range(len(rows)):
        pre = rows[: n + 1]; r = pre[-1]
        if r.get("nonfinite") or r.get("nonfinite_all"):
            return {"diverged_at": r["step"]}
        dr = [drift(pre, k) for k in KEYS[:3]]
        if r.get("deficit") is not None and None not in dr and abs(r["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr):
            return {"step": r["step"], "drift": dr, "deficit": r["deficit"], "window": [x for x in pre if x["step"] >= r["step"] - 20000]}
    return None


def qs(window, tag):
    p = D / f"lvcmp_qs_{tag}.csv"
    with open(p, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step"] + list(KEYS))
        for x in window:
            w.writerow([x["step"]] + [x[k] for k in KEYS])
    o = subprocess.run([sys.executable, str(HERE / "../../solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(p), "--series-cols", ",".join(KEYS),
                        "--tail", "1.0", "--drift", "0.0005", "--min-snaps", "5"], capture_output=True, text=True)
    txt = (o.stdout + o.stderr).strip()
    ov = [l for l in txt.splitlines() if "OVERALL" in l]
    # 実行の失敗を正常な判定と区別する (result 段 codex M3): rc は 0/1、OVERALL の行があり、4 量それぞれに分類の行が 1 つずつある。
    # ツールは入力の誤り (列の欠け・空の CSV) でも OVERALL を出して rc 1 で終わるので、OVERALL だけでは区別できない
    cls = {}
    for k in KEYS:
        m = [l.split()[-1] for l in txt.splitlines() if re.match(rf"^\s*{k}\s*:", l)]
        cls[k] = m[0] if len(m) == 1 and m[0] in QS_CLASSES else None
    ok = o.returncode in (0, 1) and bool(ov) and None not in cls.values()
    return {"rc": o.returncode, "overall": ov[-1] if ov else None, "classes": cls, "ok": ok, "tail": txt.splitlines()[-6:]}


QS_CLASSES = ("STEADY", "DRIFTING", "OSCILLATING", "TRANSIENT-UNSETTLED")
REQ_RES = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega")


def segverdict(run):
    """区間付きの収束 VERDICT。判定不能・必須の列の欠け・VERDICT なしは (None, 理由, 照合の記録) を返す (codex plan-2 M2)。
    保存のテキストがその run の残差の履歴から作られたことを lvc_segcheck.recheck で照合する (result 段 codex M3)"""
    import lvc_segcheck
    p = HERE / run / "CONVERGENCE_SEGMENT.txt"
    t = p.read_text(errors="replace") if p.is_file() else ""
    m = re.search(r"->\s*([A-Z][A-Z ]+[A-Z])", t)
    if not (m and "判定区間" in t):
        return None, "区間付きの VERDICT が無い", None
    if "判定不能" in t:
        return None, "判定不能を含む", None
    miss = [c for c in REQ_RES if not re.search(rf"^\s*{c}\s*:", t, re.M)]
    if miss:
        return None, f"必須の残差の列が無い {miss}", None
    ok, rec = lvc_segcheck.recheck(HERE / run)
    if not ok:
        return None, f"保存の VERDICT を残差の履歴と照合できない: {rec.get('why')}", rec
    return m.group(1), None, rec


def input_gates():
    g = []
    src_full = sha(HERE / SRC / SRC_RES)
    g.append(("出発の場の sha256", src_full.startswith(SRC_SHA16), src_full[:16]))
    ref_mesh = meshsha(HERE / SRC / "nozzle.h5")
    for name, (run, mk, _) in ARMS.items():
        d = HERE / run
        try:
            prov = (d / "RUN_PROVENANCE.txt").read_text(errors="replace")
            g.append((f"{run}: forge の sha256", SHA_NEW in prov, SHA_NEW[:16]))
            cp = json.loads((d / "COLD_PAIR.json").read_text())
            g.append((f"{run}: 出発の場", cp["parent"] == SRC and cp["parent_res"] == SRC_RES and cp["parent_res_sha256"] == src_full, cp["parent_res_sha256"][:16]))
            m = meshsha(d / "nozzle.h5")
            g.append((f"{run}: 格子の実体", m == ref_mesh and (d / "MESH_SHA.txt").read_text().strip() == m, m[:16]))
            diff_in = [fn for fn in SAME_FILES if not (d / fn).is_file() or sha(HERE / SRC / fn) != sha(d / fn)]
            g.append((f"{run}: 入力ファイル", diff_in == [], diff_in))
            c = cfg(run)
            exp_ = {**EXPECT, **ARM_CFG[KIND[name]]}
            bad = {k: c.get(k) for k, v in exp_.items() if (k in c if v is None else c.get(k) != v)}
            g.append((f"{run}: 設定が腕 {name} の期待どおり", bad == {}, bad))
            log = (d / "forge_run.log").read_text(errors="replace")
            g.append((f"{run}: LAYOUT2 の分解・代入", "[line] factor 1 回目: モード LAYOUT2" in log and "[line] solve 1 回目: モード LAYOUT2" in log, None))
            g.append((f"{run}: 面エンタルピーの切替なし", "FORGE_DIAG_FACE_H_DOUBLE" not in log, None))
            vshow = "'lineViscCoupling' in 'time.deltaT': 3" in log
            mshow = "診断のマスク FORGE_LVC_TERMS=5 " in log
            okv = {"L7": vshow and "FORGE_LVC_TERMS=" not in log, "L5": vshow and mshow, "L0": ("'lineViscCoupling' in 'time.deltaT'" not in log) and "FORGE_LVC_TERMS=" not in log}[KIND[name]]
            g.append((f"{run}: 値とマスクの表示が腕 {name} と合う", okv, None))
            lt = (d / "LVC_TERMS.txt").read_text().strip() if (d / "LVC_TERMS.txt").is_file() else None
            g.append((f"{run}: 台本が渡した FORGE_LVC_TERMS = {mk}", lt == mk, lt))
        except Exception as e:
            g.append((f"{run}: 入力のゲートを評価できない", False, f"{type(e).__name__}: {e}"))
    try:
        vr = subprocess.run([sys.executable, str(HERE / "lvcmp_pregate.py"), "--verify"], capture_output=True, text=True)
        g.append(("事前のゲートの記録が PASS で、証拠のハッシュが今のファイルと同じ", vr.returncode == 0, (vr.stdout + vr.stderr).strip()[-300:]))
        pg = json.loads((D / "lvcmp_pregate.json").read_text())
        dcs = pg.get("dump_config") or {}
        allowed = {"time/last/nStepOuter", "time/outStepInterval", "output/extraFields"}
        for name, (run, _, dk) in ARMS.items():
            c = cfg(run); dc = dcs.get(dk) or {}
            dd = sorted(k for k in set(c) | set(dc) if k not in allowed and c.get(k) != dc.get(k))
            g.append((f"{run}: 設定が同じ腕の書き出しと同じ (差は step 数・出力の間隔・extraFields だけ)", bool(dc) and dd == [], dd[:10]))
        for run in DUMPS:
            g.append((f"{run}: 書き出しが残っている", (HERE / run / "linedump" / "D.f64").is_file(), None))
    except Exception as e:
        g.append(("事前のゲートの記録を読めない", False, f"{type(e).__name__}: {e}"))
    return g


def assess(name):
    run = ARMS[name][0]
    info = {"run": run}
    p = HERE / run / "m9_watch.json"
    if not p.is_file():
        return "INVALID", {**info, "why": "見張りの記録が無い"}
    w = json.loads(p.read_text())
    st = w.get("status"); info["watch_status"] = st
    if st in ("EXEC_ERROR", "DATA_ERROR") or st is None:
        return "INVALID", {**info, "why": f"見張りの状態 {st}", "note": w.get("note")}
    if not (w.get("run") == run and w.get("phase") == "line" and w.get("budget") == BUDGET and w.get("consec", 1) == 1 and not w.get("inherit")):
        return "INVALID", {**info, "why": "見張りの記録のスキーマ (run・phase・budget・consec) が登録と違う"}
    rows = sorted([x for x in (w.get("rows") or []) if x["step"] > 0], key=lambda x: x["step"])
    steps = [x["step"] for x in rows]
    info["n_rows"] = len(rows)
    if st == "DIVERGED":                           # 初回の出力 (2500 step) の前の発散は系列が空 (codex plan-3 M5)
        nf_rows = [x["step"] for x in rows if x.get("nonfinite") or x.get("nonfinite_all")]
        nan = sorted((HERE / run).glob("res_nan_*.h5"))
        logp = HERE / run / "forge_run.log"
        log = logp.read_text(errors="replace") if logp.is_file() else ""
        fs = w.get("fail_step")
        if not (nf_rows or (nan and "Non-finite value detected" in log)):
            return "INVALID", {**info, "why": "DIVERGED の証拠 (非有限の出力の行、res_nan と detectNaN のログ) が無い"}
        fin = [n for n in steps if n not in nf_rows]
        if fin and fin != list(range(OUT, fin[-1] + 1, OUT)):
            return "INVALID", {**info, "why": "発散の前の系列が 2500 step ごとにそろっていない"}
        if not fin and not (fs is not None and int(fs) <= OUT + OUT):
            return "INVALID", {**info, "why": f"系列が空なのに失敗の step ({fs}) が初回の出力の前でない"}
        info["diverged_at"] = fs
        return "DIVERGED", info
    if not steps or len(set(steps)) != len(steps) or steps != list(range(OUT, steps[-1] + 1, OUT)):
        return "INVALID", {**info, "why": f"系列が {OUT} step ごとに欠け・重複なくそろっていない (行 {len(steps)}、最後 {steps[-1] if steps else None})"}
    # 全行の量の存在と有限性 (到達の行だけでなく途中の行も。result 段 codex M3)。派生の列は要求しない
    import math
    badq = [x["step"] for x in rows if not all(x.get(k) is not None and math.isfinite(float(x[k])) for k in ("deficit",) + KEYS)]
    if badq:
        return "INVALID", {**info, "why": f"見張りの行の量が欠けているか非有限 ({len(badq)} 行、最初 step {badq[0]})"}
    # 記録係: 系列の全出力について sha256・全保存量の有限・ρ/P/T > 0 (見張りが消した出力の代わりの証拠、codex plan-2 m5)
    rp = HERE / run / "lvcmp_outputs.json"
    recs = json.loads(rp.read_text()).get("rows", {}) if rp.is_file() else {}
    last_ok = steps[-1] if st != "DIVERGED" else max([x["step"] for x in rows if not (x.get("nonfinite") or x.get("nonfinite_all"))], default=0)
    bad = [n for n in steps if n <= last_ok and not (str(n) in recs and recs[str(n)].get("nonfinite_all") == 0 and recs[str(n)].get("positive"))]
    if bad:
        return "INVALID", {**info, "why": f"記録係の記録が無いか、非有限・非正の出力がある ({len(bad)} 件、最初 {bad[0]})"}
    seg, why, segrec = segverdict(run)
    info["segment_verdict"] = seg
    info["segment_recheck"] = segrec
    rl = reach_level(rows)
    if st == "REACHED":
        n = w.get("reach_step")
        info["reach_step_watch"] = n
        if rl is None or "step" not in rl or int(rl["step"]) != int(n) or steps[-1] != int(n):
            return "INVALID", {**info, "why": f"到達の step が系列からの再計算と合わない ({rl and rl.get('step')} ≠ {n}、系列の最後 {steps[-1]})"}
        f = HERE / run / f"res_{n}.h5"
        if not f.is_file() or sha(f) != w.get("reach_sha256") or recs.get(str(n), {}).get("sha256") != w.get("reach_sha256"):
            return "INVALID", {**info, "why": "到達の出力が無いか、sha256 が見張り・記録係の記録と違う"}
        import cold_series as CS, cold_xcheck as XC
        yb_x, yb = XC.common_yb()
        again = CS.snapshot(HERE / run, int(n), yb_x, yb)
        row = rows[-1]
        import math
        vals = [(again.get(k), row.get(k)) for k in ("deficit",) + KEYS]
        if not all(a_ is not None and b_ is not None and math.isfinite(float(a_)) and math.isfinite(float(b_)) for a_, b_ in vals):
            return "INVALID", {**info, "why": "到達の出力の量 (記録・計算し直し) に欠けか非有限がある (codex plan-3 M3)"}
        diff = {k: (again.get(k), row.get(k)) for k in ("deficit",) + KEYS
                if abs(float(again[k]) - float(row[k])) > 1e-12 * max(abs(float(again[k])), abs(float(row[k])), 1e-300)}
        if diff:
            return "INVALID", {**info, "why": f"到達の出力から計算し直した量が系列の行と違う {diff}"}
        if seg is None or seg.startswith("DIVERGED"):
            return "INVALID", {**info, "why": f"区間付きの収束: {why or seg}"}
        if segrec["stored"]["last_step"] < int(n):
            return "INVALID", {**info, "why": f"区間付きの収束の last step {segrec['stored']['last_step']} が到達の step {n} より前"}
        q = qs(rl["window"], name)
        if not q["ok"]:
            return "INVALID", {**info, "quasisteady": q, "why": "準定常の判定を実行できなかった (rc・OVERALL・4 量の分類の行のどれかが無い)"}
        info.update(reach_step=int(n), drift=rl["drift"], deficit=rl["deficit"], quasisteady=q)
        return "REACHED", info
    if st == "DIVERGED":
        nf_rows = [x["step"] for x in rows if x.get("nonfinite") or x.get("nonfinite_all")]
        nan = sorted((HERE / run).glob("res_nan_*.h5"))
        log = (HERE / run / "forge_run.log").read_text(errors="replace") if (HERE / run / "forge_run.log").is_file() else ""
        if not (nf_rows or (nan and "Non-finite value detected" in log)):
            return "INVALID", {**info, "why": "DIVERGED の証拠 (非有限の出力の行、res_nan と detectNaN のログ) が無い"}
        info["diverged_at"] = w.get("fail_step")
        return "DIVERGED", info
    if st == "CENSORED":
        rc = (HERE / run / "RUN_RC").read_text().strip() if (HERE / run / "RUN_RC").is_file() else None
        if rc != "0" or steps[-1] != BUDGET or rl is not None:
            return "INVALID", {**info, "why": f"CENSORED の条件 (RUN_RC 0・最後の出力 {BUDGET}・系列が水準に入っていない) を満たさない (rc {rc}、最後 {steps[-1]})"}
        if seg is None or seg.startswith("DIVERGED"):
            return "INVALID", {**info, "why": f"区間付きの収束: {why or seg}"}
        return "CENSORED", info
    return "INVALID", {**info, "why": f"未知の状態 {st}"}


def cmp(cls, x, y):
    cx, cy = cls[x], cls[y]
    if cx[0] == "REACHED" and cy[0] == "REACHED":
        d = cx[1]["reach_step"] - cy[1]["reach_step"]
        v = (f"{x} が速い" if d < 0 else f"{x} が遅い") if abs(d) >= THR else "差を検出しなかった"
        return {"verdict": v, "diff_steps": d, "x": cx[1]["reach_step"], "y": cy[1]["reach_step"]}
    return {"verdict": "比較しない", "x": cx[0], "y": cy[0]}


def main():
    rec = {"plan": "time_integration-line-viscous-jacobian-dt §6.3"}
    out = D / "lvcmp_judge.json"
    try:
        gates = input_gates()
    except Exception as e:
        gates = [("入力のゲートを評価できない", False, f"{type(e).__name__}: {e}")]
    cls = {}
    for n in ARMS:
        try:
            cls[n] = assess(n)
        except Exception as e:
            cls[n] = ("INVALID", {"why": f"例外 {type(e).__name__}: {e}"})
    rec["gates"] = [{"check": g[0], "ok": bool(g[1]), "value": g[2]} for g in gates]
    rec["gates_ok"] = all(g[1] for g in gates)
    rec["arms"] = {n: {"class": c[0], **c[1]} for n, c in cls.items()}
    invalid = [n for n in ARMS if cls[n][0] == "INVALID"]
    if not rec["gates_ok"] or invalid:
        rec["VERDICT"] = f"INVALID (ゲート不合格 {sum(not g[1] for g in gates)} 件、INVALID の腕 {invalid}): 判定しない"
    else:
        c70, c75, c50 = cmp(cls, "L7", "L0"), cmp(cls, "L7", "L5"), cmp(cls, "L5", "L0")
        rec["comparisons"] = {"L7_vs_L0": c70, "L7_vs_L5": c75, "L5_vs_L0": c50}
        s7, s0 = cls["L7"][0], cls["L0"][0]
        # 各腕 1 本の探索的な比較 (codex plan-2 M3): 再現性を含む結論 (手法の利点・熱伝導の K の効き) はここでは言わない
        if s0 != "REACHED":
            main_v = "判別不能 (生産 L0 が今回の実行では水準に届かなかった)"
        elif s7 != "REACHED":
            main_v = f"NO-GO (全部入り L7 は今回の実行では 15 万 step 以内に水準に届かなかった: {s7})"
        elif c70["verdict"] == "L7 が速い":
            main_v = "今回の実行では、全部入りが生産より登録幅 (15000 step) 以上早く水準に届いた (探索的。再現は別に確かめる)"
        elif c70["verdict"] == "L7 が遅い":
            main_v = "NO-GO (今回の実行では、全部入りが生産より登録幅以上遅かった)"
        else:
            main_v = "今回の実行では、全部入りと生産の到達の差を検出しなかった (登録幅 15000 step 未満)"
        rec["main"] = main_v
        rec["VERDICT"] = f"主 (L7 と L0): {main_v}; L7 と L5: {c75['verdict']}; L5 と L0: {c50['verdict']}"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for g in gates:
        if not g[1]:
            print("ゲート不合格:", g[0], g[2])
    for n, c in cls.items():
        print(n, c[0], {k: c[1].get(k) for k in ("run", "reach_step", "diverged_at", "segment_verdict", "why")})
    print(rec["VERDICT"])
    return 0 if rec["gates_ok"] and not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
