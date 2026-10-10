#!/usr/bin/env python3
"""point 仕上げでの面エンタルピーの精度の A/B の判定 (plan time_integration-implicit-thermal-jacobian §6.2、2026-10-10 事前登録)。

fh_floor.sh の run (軌道 run_0500〜0505、評価 run_0510〜0537) を読み、§6.2 の規則どおりに判定して
_band_ab/cold_pair/fh_floor_judge.json に書く。規則を変えるときは plan §6.2 を先に改訂する (結果を見てから変えない)。

  主判定の量: 評価器 d (面エンタルピーを double で評価) の rms_roe を、A・B の 2000 step の状態で比べる。
  E_X(状態) = 評価 run の step 0 の outer_begin の行 (更新前の残差)。X = f (既定) / d (切替あり)。
  ゲートが 1 つでも外れたら量を判定せず INVALID を出す。
"""
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
COLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
Q = "rms_roe"
S0, S0_RES = "run_0354_m9_L5cut", "res_40000.h5"
TRAJ = {"a1": ("run_0500_fh_a1", 0), "b1": ("run_0501_fh_b1", 1), "a2": ("run_0502_fh_a2", 0),
        "b2": ("run_0503_fh_b2", 1), "a3": ("run_0504_fh_a3", 0), "b3": ("run_0505_fh_b3", 1)}
STATES = ["s0", "a1k1000", "a1k1500", "a1k2000", "b1k1000", "b1k1500", "b1k2000", "a2k2000", "a3k2000", "b2k2000", "b3k2000"]
REPEAT = ["s0", "a1k2000", "b1k2000"]
# 入力の同一性を run_0354 と比べるファイル (prep が run_0183 から複製するもの + 解決済みの化学種)
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv",
              "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
CONS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1", "P", "T"]
# 事前登録の閾値 (§6.2)
DROP = 0.10          # 検出の対象にする低下 (10 %)
NOISE_MULT = 3.0     # 不確かさの幅 u = 3 × noise (統計的な 3σ ではない。探索試験の目安)
DRIFT = 0.01         # 傾き × 500 と窓の水準差の許容 [桁] (0.01 桁 ≈ 2.3 %)
GATE_SELF = 1e-6     # 軌道の step 0 の表示と、同じ状態・同じ評価器の評価の相対差の上限
TAIL, PREV = (1500, 1999), (1000, 1499)


def eval_runs():
    """(状態, 評価器, 回) → run 名。fh_floor.sh と同じ並び。"""
    out, k = {}, 510
    for sd in STATES:
        for x in ("f", "d"):
            out[(sd, x, 1)] = f"run_{k:04d}_fhe_{sd}_{x}"
            k += 1
    for sd in REPEAT:
        for x in ("f", "d"):
            out[(sd, x, 2)] = f"run_{k:04d}_fhe_{sd}_{x}r"
            k += 1
    return out


def read_csv(run):
    """全行 (phase ごと) と outer_begin の行 (step → 値)、行の重複・非有限・負を返す。"""
    allrows, outer, dup, bad = [], {}, [], []
    with open(HERE / run / "residual_history.csv") as f:
        rd = csv.DictReader(f)
        rcols = [c for c in rd.fieldnames if c.startswith("rms_")]
        for r in rd:
            vals = {c: float(r[c]) for c in rcols}
            if any((not math.isfinite(v)) or v < 0 for v in vals.values()):
                bad.append((r["step"], r["phase"]))
            allrows.append(r)
            if r["phase"] == "outer_begin":
                s = int(r["step"])
                if s in outer:
                    dup.append(s)
                outer[s] = vals
    return outer, dup, bad, len(allrows)


def flat(d, p=""):
    out = {}
    for k, v in d.items():
        kk = f"{p}/{k}" if p else str(k)
        if isinstance(v, dict):
            out.update(flat(v, kk))
        else:
            out[kk] = v
    return out


def cfg_diff(run, ref):
    a = flat(yaml.safe_load((HERE / run / "solverConfig.yaml").read_text()))
    b = flat(yaml.safe_load((HERE / ref / "solverConfig.yaml").read_text()))
    skip = ("time/last/nStepOuter", "time/outStepInterval")
    return sorted(k for k in set(a) | set(b) if k not in skip and a.get(k) != b.get(k))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def input_diff(run, ref):
    out = [fn for fn in SAME_FILES if (HERE / ref / fn).is_file() != (HERE / run / fn).is_file()
           or ((HERE / ref / fn).is_file() and sha(HERE / ref / fn) != sha(HERE / run / fn))]
    rs = lambda d: {p.name: sha(p) for p in sorted((HERE / d).glob("resolved_species_*.yaml"))}
    if rs(run) != rs(ref):
        out.append("resolved_species_*.yaml")
    return out


def switch_on(run):
    return "FORGE_DIAG_FACE_H_DOUBLE" in (HERE / run / "forge_run.log").read_text(errors="replace")


def rc_ok(run):
    p = HERE / run / "RUN_RC"
    return p.is_file() and p.read_text().strip() == "0"


def rel(a, b):
    return abs(a - b) / (0.5 * (abs(a) + abs(b)))


def seg_verdict(run):
    p = HERE / run / "CONVERGENCE_SEGMENT.txt"
    if not p.is_file():
        return "判定不能 (CONVERGENCE_SEGMENT.txt が無い)"
    m = re.search(r"->\s*([A-Z][A-Z ]+[A-Z])", p.read_text(errors="replace"))
    return m.group(1) if m else "判定不能 (VERDICT の行が無い)"


def drift(rows):
    """末尾 500 step の log10(rms_roe) の傾き × 500 と、末尾の窓とその前の窓の中央値の差 [桁]、振動の大きさ。"""
    s = np.arange(TAIL[0], TAIL[1] + 1, dtype=float)
    y = np.log10([rows[int(k)][Q] for k in s])
    c = np.polyfit(s, y, 1)
    slope = float(c[0] * 500.0)
    lvl = float(np.log10(np.median([rows[k][Q] for k in range(*TAIL)] + [rows[TAIL[1]][Q]]) /
                         np.median([rows[k][Q] for k in range(*PREV)] + [rows[PREV[1]][Q]])))
    osc = float(np.std(y - np.polyval(c, s)) * math.log(10))   # 傾きを除いた残りの相対の揺れ (step ごと)
    if slope <= -DRIFT and lvl <= -DRIFT:
        cls = "減衰"
    elif slope >= DRIFT and lvl >= DRIFT:
        cls = "増加"
    elif abs(slope) < DRIFT and abs(lvl) < DRIFT:
        cls = "許容内"
    else:
        cls = "混在"
    return {"slope500": slope, "level_tail_vs_prev": lvl, "osc_rel": osc, "class": cls}


def main():
    gates, rec = [], {"plan": "time_integration-implicit-thermal-jacobian §6.2", "quantity": Q}
    ev = eval_runs()
    st = HERE / "_fh_states"
    mesh_ref = (st / "MESH_SHA_REF.txt").read_text().strip()
    s0sha = (st / "s0" / "SHA256").read_text().strip()
    gates.append(("S の sha256 が run_0354/res_40000.h5 と一致", s0sha == sha(HERE / S0 / S0_RES), s0sha[:16]))

    def common(run, state_name):
        gates.append((f"{run}: RUN_RC = 0", rc_ok(run), (HERE / run / "RUN_RC").read_text().strip() if (HERE / run / "RUN_RC").is_file() else None))
        d = cfg_diff(run, S0)
        gates.append((f"{run}: solverConfig が run_0354 と同じ (出力の間隔・step 数を除く)", d == [], d))
        d = input_diff(run, S0)
        gates.append((f"{run}: 境界条件・化学種・壁などの入力が run_0354 と同じ", d == [], d))
        m = (HERE / run / "MESH_SHA.txt").read_text().strip() if (HERE / run / "MESH_SHA.txt").is_file() else ""
        gates.append((f"{run}: 格子 (MESH) が run_0354 と同じ", m == mesh_ref, m[:16]))
        cp = json.loads((HERE / run / "COLD_PAIR.json").read_text())
        want = (st / state_name / "SHA256").read_text().strip()
        gates.append((f"{run}: 出発の場が状態 {state_name}", cp["field_from"] == state_name and cp["parent_res_sha256"] == want,
                      [cp["field_from"], cp["parent_res"], cp["parent_res_sha256"][:16]]))

    # --- ゲート: 軌道 ---
    traj_rows = {}
    for name, (run, fh) in TRAJ.items():
        common(run, "s0")
        rows, dup, bad, n = read_csv(run)
        traj_rows[name] = rows
        gates.append((f"{run}: outer_begin の行が step 0〜1999 で一意・連続", sorted(rows) == list(range(2000)) and not dup,
                      [min(rows, default=None), max(rows, default=None), len(rows), dup[:3]]))
        gates.append((f"{run}: 全行の rms_* が有限・非負", not bad, bad[:3]))
        gates.append((f"{run}: 切替の表示 = {fh}", switch_on(run) == bool(fh), switch_on(run)))
        with h5py.File(HERE / run / "res_2000.h5", "r") as h:
            nonfin = {v: int(np.sum(~np.isfinite(h["VALUE"][v][...]))) for v in CONS}
            nonpos = {v: int(np.sum(h["VALUE"][v][...] <= 0)) for v in ("ro", "P", "T")}
        gates.append((f"{run}: res_2000 の保存量・P・T が有限", sum(nonfin.values()) == 0, nonfin))
        gates.append((f"{run}: res_2000 の ρ・P・T が正", sum(nonpos.values()) == 0, nonpos))
    # --- ゲート: 評価 ---
    E = {}
    for (sd, x, n), run in ev.items():
        common(run, sd)
        rows, dup, bad, nrow = read_csv(run)
        gates.append((f"{run}: step 0 の outer_begin が 1 行", 0 in rows and not dup, [sorted(rows)[:3], dup]))
        gates.append((f"{run}: 全行の rms_* が有限・非負", not bad, bad[:3]))
        gates.append((f"{run}: 切替の表示 = {x}", switch_on(run) == (x == "d"), switch_on(run)))
        e = rows.get(0)
        gates.append((f"{run}: 判定に使う値が正", e is not None and all(e[c] > 0 for c in COLS if c != "rms_roUz"), e and e[Q]))
        E[(sd, x, n)] = e
    ok = all(g[1] for g in gates)
    if ok:   # 自己一致 (保守的な停止の基準。外れても「別の残差を読んだ」とは断定しない)
        self_f = rel(traj_rows["a1"][0][Q], E[("s0", "f", 1)][Q])
        self_d = rel(traj_rows["b1"][0][Q], E[("s0", "d", 1)][Q])
        gates.append(("A1 の step 0 の表示と E_f(S) の相対差 ≤ 1e-6", self_f <= GATE_SELF, self_f))
        gates.append(("B1 の step 0 の表示と E_d(S) の相対差 ≤ 1e-6", self_d <= GATE_SELF, self_d))
        ok = all(g[1] for g in gates)
    rec["gates"] = [{"check": g[0], "ok": bool(g[1]), "value": g[2]} for g in gates]
    rec["gates_ok"] = bool(ok)
    if not ok:
        rec["VERDICT"] = "INVALID (ゲート不合格: 量を判定しない)"
        (HERE / "_band_ab" / "cold_pair" / "fh_floor_judge.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
        for g in gates:
            if not g[1]:
                print("ゲート不合格:", g[0], g[2])
        print(rec["VERDICT"])
        return 1
    # --- 量 ---
    rec["E"] = {f"{sd}|{x}|{n}": E[(sd, x, n)] for (sd, x, n) in E}
    sig = {x: max(rel(E[(sd, x, 1)][Q], E[(sd, x, 2)][Q]) for sd in REPEAT) for x in ("f", "d")}
    rec["sigma_eval"] = sig
    rec["eval_effect_d_over_f"] = {sd: {c: E[(sd, "d", 1)][c] / E[(sd, "f", 1)][c] for c in COLS if E[(sd, "f", 1)][c] > 0} for sd in STATES}
    per = {}
    for x in ("f", "d"):
        A = [E[(f"a{i}k2000", x, 1)][Q] for i in (1, 2, 3)]
        Bv = [E[(f"b{i}k2000", x, 1)][Q] for i in (1, 2, 3)]
        mA, mB = float(np.mean(A)), float(np.mean(Bv))
        wA, wB = (max(A) - min(A)) / mA, (max(Bv) - min(Bv)) / mB
        noise = max(sig[x], wA, wB)
        u = NOISE_MULT * noise
        r = mB / mA
        dec = 1.0 - r
        if dec - u >= DROP and max(Bv) <= (1 - DROP) * min(A):
            v = "SUPPORT"
        elif dec + u < DROP:
            v = "NOT_SUPPORT"
        else:
            v = "INDETERMINATE"
        per[x] = {"A": A, "B": Bv, "mean_A": mA, "mean_B": mB, "ratio": r, "decrease": dec, "w_A": wA, "w_B": wB,
                  "noise": noise, "u": u, "verdict": v,
                  "all_cols_ratio": {c: float(np.mean([E[(f"b{i}k2000", x, 1)][c] for i in (1, 2, 3)]) /
                                         np.mean([E[(f"a{i}k2000", x, 1)][c] for i in (1, 2, 3)]))
                                     for c in COLS if all(E[(f"a{i}k2000", x, 1)][c] > 0 for i in (1, 2, 3))}}
    rec["by_evaluator"] = per
    rec["f_d_concordant"] = per["f"]["verdict"] == per["d"]["verdict"]
    # 自方式の表示 (各腕の CSV) と、反復ごとのドリフト
    disp = {}
    for name, rows in traj_rows.items():
        disp[name] = {"tail_median": {c: float(np.median([rows[k][c] for k in range(TAIL[0], TAIL[1] + 1)])) for c in COLS},
                      "step0": rows[0], "drift": drift(rows), "segment_verdict": seg_verdict(TRAJ[name][0])}
    rec["display"] = disp
    dA = float(np.mean([disp[f"a{i}"]["tail_median"][Q] for i in (1, 2, 3)]))
    dB = float(np.mean([disp[f"b{i}"]["tail_median"][Q] for i in (1, 2, 3)]))
    rec["display_ratio_B_over_A"] = dB / dA
    ee_A = float(np.mean([E[(f"a{i}k2000", "d", 1)][Q] / E[(f"a{i}k2000", "f", 1)][Q] for i in (1, 2, 3)]))
    rec["eval_effect_at_A2000"] = ee_A
    floor = {}
    for arm in ("a", "b"):
        cl = [disp[f"{arm}{i}"]["drift"]["class"] for i in (1, 2, 3)]
        c = cl[0] if len(set(cl)) == 1 else "反復間不一致"
        floor[arm] = {"classes": cl, "class": c,
                      "text": {"許容内": "この窓で有意なドリフトを検出せず", "減衰": "減衰中 = 床の判定は不能 (自動で延長しない)",
                               "増加": "増加 (頭打ちではない)"}.get(c, "判定不能 (" + c + ")")}
    rec["drift"] = floor
    rec["series_vs_S"] = {f"{arm}1": {x: {k: E[(f"{arm}1k{k}", x, 1)][Q] / E[("s0", x, 1)][Q] for k in (1000, 1500, 2000)}
                                      for x in ("f", "d")} for arm in ("a", "b")}
    # 主判定 (評価器 d)
    v = per["d"]["verdict"]
    txt = {"SUPPORT": "この期間の停滞への寄与を支持 (d で評価した離散残差の改善に限る)",
           "NOT_SUPPORT": "この期間の主要因説を支持しない",
           "INDETERMINATE": "判別不能 (10 % の境界を不確かさの幅がまたぐ、または範囲が分かれない)"}[v]
    if v != "SUPPORT" and rec["display_ratio_B_over_A"] <= 1 - DROP:
        u = per["d"]["u"]
        only = abs(1 - per["d"]["ratio"]) <= u and abs(rec["display_ratio_B_over_A"] - ee_A) <= u
        txt += ("; 自方式の表示では 10 % 以上下がったが共通の評価 (d) では基準に届かない" +
                (" — 同じ状態の f/d の差で表示の差を説明できるので、表示の定義の違いだけ" if only else ""))
    rec["VERDICT"] = f"{v}: {txt}; ドリフト: A {floor['a']['text']}・B {floor['b']['text']}; f と d の判定の一致: {rec['f_d_concordant']}"
    out = HERE / "_band_ab" / "cold_pair" / "fh_floor_judge.json"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=float))
    print(rec["VERDICT"])
    print("sigma_eval", sig, "表示の比 B/A", round(rec["display_ratio_B_over_A"], 5), "A の 2000 での d/f", round(ee_A, 5))
    for x in ("f", "d"):
        p = per[x]
        print(f"評価器 {x}: A {['%.5g' % a for a in p['A']]} B {['%.5g' % b for b in p['B']]} 比 {p['ratio']:.4f} u {p['u']:.3e} → {p['verdict']}")
    for name in TRAJ:
        dd = disp[name]["drift"]
        print(f"{name}: 傾き×500 {dd['slope500']:+.4f} 窓の差 {dd['level_tail_vs_prev']:+.4f} 揺れ {dd['osc_rel']:.2e} {dd['class']} / {disp[name]['segment_verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
