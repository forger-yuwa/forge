#!/usr/bin/env python3
"""上限 50 での全部入り・粘性入り・生産の比較の判定 (plan time_integration-line-viscous-jacobian-dt §6.3、2026-10-11 事前登録)。

lvcmp.sh の run_0553 (L7 = 値 3・マスク 7)・run_0554 (L0 = 値 0、生産)・run_0555 (L5 = 値 3・マスク 5) と、事前のゲートの記録
(lvcmp_pregate.json、run_0556〜0559 の 1 step の書き出しから) を読み、
  (1) 入力のゲート: バイナリ・出発の場・格子の実体・入力ファイル・腕ごとの期待どおりの設定・起動ログ (値・マスク・LAYOUT2)、
      事前のゲートの --verify、腕の設定が同じ腕の書き出しの設定と step 数・出力の間隔・extraFields 以外で同じ
  (2) 腕ごと: 見張り (m9_watch.json) の状態。REACHED なら、到達の step を系列から水準の規則で計算し直して一致すること、
      区間付きの収束 VERDICT があって DIVERGED でないこと。EXEC_ERROR・DATA_ERROR・記録なしは INVALID
  (3) 比較: 2 本とも REACHED のとき、到達の step の差が 15000 以上なら速い/遅い、未満は判別不能。どちらかが届かなければ比較しない
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
D = HERE / "_band_ab" / "cold_pair"
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
    return {"rc": o.returncode, "overall": ov[-1] if ov else None, "tail": txt.splitlines()[-6:]}


def segverdict(run):
    p = HERE / run / "CONVERGENCE_SEGMENT.txt"
    t = p.read_text(errors="replace") if p.is_file() else ""
    m = re.search(r"->\s*([A-Z][A-Z ]+[A-Z])", t)
    return (m.group(1) if (m and "判定区間" in t) else None), t


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
            exp_ = {**EXPECT, **ARM_CFG[name]}
            bad = {k: c.get(k) for k, v in exp_.items() if (k in c if v is None else c.get(k) != v)}
            g.append((f"{run}: 設定が腕 {name} の期待どおり", bad == {}, bad))
            log = (d / "forge_run.log").read_text(errors="replace")
            g.append((f"{run}: LAYOUT2 の分解・代入", "[line] factor 1 回目: モード LAYOUT2" in log and "[line] solve 1 回目: モード LAYOUT2" in log, None))
            g.append((f"{run}: 面エンタルピーの切替なし", "FORGE_DIAG_FACE_H_DOUBLE" not in log, None))
            vshow = "'lineViscCoupling' in 'time.deltaT': 3" in log
            mshow = "診断のマスク FORGE_LVC_TERMS=5 " in log
            okv = {"L7": vshow and "FORGE_LVC_TERMS=" not in log, "L5": vshow and mshow, "L0": ("'lineViscCoupling' in 'time.deltaT'" not in log) and "FORGE_LVC_TERMS=" not in log}[name]
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
    rows = w.get("rows") or []
    info["n_rows"] = len(rows)
    seg, _ = segverdict(run)
    info["segment_verdict"] = seg
    if st in ("EXEC_ERROR", "DATA_ERROR") or st is None:
        return "INVALID", {**info, "why": f"見張りの状態 {st}", "note": w.get("note")}
    rl = reach_level(rows)
    if st == "REACHED":
        info["reach_step_watch"] = w.get("reach_step")
        if rl is None or "step" not in rl or int(rl["step"]) != int(w.get("reach_step", -1)):
            return "INVALID", {**info, "why": f"到達の step が系列からの再計算と合わない ({rl and rl.get('step')} ≠ {w.get('reach_step')})"}
        if seg is None or seg.startswith("DIVERGED"):
            return "INVALID", {**info, "why": f"区間付きの収束 VERDICT が無いか DIVERGED ({seg})"}
        info.update(reach_step=int(rl["step"]), drift=rl["drift"], deficit=rl["deficit"], quasisteady=qs(rl["window"], name))
        return "REACHED", info
    if st == "DIVERGED":
        info["diverged_at"] = (rl or {}).get("diverged_at")
        return "DIVERGED", info
    if st == "CENSORED":
        if rl is not None:
            return "INVALID", {**info, "why": "CENSORED なのに系列は水準に入っている"}
        return "CENSORED", info
    return "INVALID", {**info, "why": f"未知の状態 {st}"}


def cmp(cls, x, y):
    cx, cy = cls[x], cls[y]
    if cx[0] == "REACHED" and cy[0] == "REACHED":
        d = cx[1]["reach_step"] - cy[1]["reach_step"]
        v = (f"{x} が速い" if d < 0 else f"{x} が遅い") if abs(d) >= THR else "判別不能"
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
        if s0 != "REACHED":
            main_v = "判別不能 (生産 L0 が水準に届かない)"
        elif s7 != "REACHED":
            main_v = f"NO-GO (全部入り L7 は上限 50 で水準に届かない: {s7})"
        elif c70["verdict"] == "L7 が速い":
            main_v = "全部入りが生産より早く水準に届く (この条件)"
        elif c70["verdict"] == "L7 が遅い":
            main_v = "NO-GO (全部入りは生産より遅い)"
        else:
            main_v = "全部入りと生産は区別できない"
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
