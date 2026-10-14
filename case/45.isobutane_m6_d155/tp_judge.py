"""plan architecture-float-state-double-geometry §6.62 (試験 P、事前登録・2026-10-13 の資格の訂正の後) の段 1 の判定: 打ち切りの規則。
P0 = run_0806_tp_p0 (lowMachPrecond 0)、P3 = run_0807_tp_p3 (lowMachPrecond 3)、各 20,000 step・出力 500 ごと、B で占有して逐次。
資格: _band_ab/cold_pair/tp_check.json (tp_check.py)。不成立なら「資格不成立」、判定不能が残れば「判定不能」で止める (窓の判定をしない)。
窓: W1 = step 15,000 以上 17,500 未満、W2 = 17,500 以上 20,000 以下。各指標・各腕・各窓で有限の点が 5 点以上 (足りなければ判定不能)。
  - 全域の指標: residual_history.csv の phase = outer_end の行の全 rms_* 列 (rms_dq_*・rms_roUz を含む)。step は 0 始まり
    (nStepOuter 20,000 なら最後の行は 19,999) なので、W2 は行 17,500〜19,999。
  - 局所の指標: tp_series.jsonl (tp_region.py) の区間 x/r_t ∈ [−8, −2) の全内部節点の Σ|res_q| (q = ro・roUx・roUy・roe・roK・roOmega)。
    res_N は N step 後の場なので、W1 は出力 15,000〜17,000 の 5 点、W2 は 17,500〜20,000 の 6 点。全行の節点の集合 (数・sha) が
    固定したもの (tp_region_set.json) と同じこと (違えば判定不能)。
持続する悪化 (指標ごと): 両窓とも 中央値(P3)/中央値(P0) ≥ 2、かつ 中央値(P3, W2)/中央値(P3, W1) ≥ 0.9。
  中央値は有限の値だけで取り、非有限の数は別に書く (比の判定から分ける)。P0 の中央値が 0 の窓がある指標は比を作らず別に書き、それだけでは打ち切らない。
結果 A (どれかの指標で持続する悪化あり): この設定の長い試験を打ち切る (「この条件では継続の対象から外した」。「低マッハの前処理は効かない」とは言わない)。
結果 B (どの指標にも無い): 段 2 へ (別に登録の追記をしてから。「P3 が優れる」とは言わない)。
記録: 全指標の窓の中央値と比、1 step の時間 (timing.tsv、占有 = 起動から終わりまで forge の数がすべて 1 以下)、step 20,000 の θ_r (x = 40・70・94)・
  Q_w・出口 M (sp_judge.py の fields と同じ式、未収束の値として)、全残差列の check_convergence。
→ _band_ab/cold_pair/tp_judge.json"""
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402

OUTD = HERE / "_band_ab" / "cold_pair"
ARMS = {"P0": "run_0806_tp_p0", "P3": "run_0807_tp_p3"}
WIN = {"W1": (15000, 17500, False), "W2": (17500, 20000, True)}   # (下限, 上限, 上限を含むか)
MIN_PTS = 5
RATIO = 2.0
PERSIST = 0.9
LOCAL_Q = ("ro", "roUx", "roUy", "roe", "roK", "roOmega")
TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
END_STEP = 20000


class Missing(Exception):
    pass


def in_win(s, w):
    lo, hi, inc = WIN[w]
    return lo <= s and (s <= hi if inc else s < hi)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_global(run):
    p = HERE / run / "residual_history.csv"
    if not p.exists():
        raise Missing(f"{run}: residual_history.csv が無い")
    with open(p) as fh:
        rd = csv.DictReader(fh)
        head = rd.fieldnames or []
        if "step" not in head or "phase" not in head:
            raise Missing(f"{run}: residual_history.csv に step・phase の列が無い ({head[:5]})")
        cols = [c for c in head if c.startswith("rms_")]
        if not cols:
            raise Missing(f"{run}: rms_* の列が無い")
        rows = [r for r in rd if r.get("phase") == "outer_end"]
    steps = [int(r["step"]) for r in rows]
    if len(steps) != len(set(steps)):
        raise Missing(f"{run}: outer_end の行の step が重複している")
    return cols, {c: [(int(r["step"]), num(r[c])) for r in rows] for c in cols}, (max(steps) if steps else None)


def load_local(run, fixed):
    p = HERE / run / "tp_series.jsonl"
    if not p.exists():
        raise Missing(f"{run}: tp_series.jsonl が無い")
    rows = {}
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        o = json.loads(line)
        if o["step"] in rows:
            raise Missing(f"{run}: tp_series.jsonl の step {o['step']} が重複している")
        rows[o["step"]] = o
    ids = {(o.get("set_n"), o.get("set_sha16")) for o in rows.values()}
    if ids != {(fixed["n"], fixed["sha16"])}:
        raise Missing(f"{run}: 節点の集合が固定したものと違う行がある ({ids} 対 {(fixed['n'], fixed['sha16'])})")
    ser = {q: [(s, num(o["sum_abs_res"].get(q)) if o.get("sum_abs_res", {}).get(q) is not None else float("nan")) for s, o in sorted(rows.items())]
           for q in LOCAL_Q}
    return ser, sorted(rows)


def win_stats(series):
    """series [(step, 値)] → 窓ごとの {n_finite, n_nonfinite, median}。"""
    out = {}
    for w in WIN:
        v = [x for s, x in series if in_win(s, w)]
        f = [x for x in v if math.isfinite(x)]
        out[w] = {"n_finite": len(f), "n_nonfinite": len(v) - len(f), "median": float(np.median(f)) if f else None}
    return out


def judge_indicator(sp0, sp3):
    """1 つの指標の判定: 点数不足 (None) / P0 の中央値 0 / 持続する悪化 (True/False)。"""
    rec = {"P0": sp0, "P3": sp3}
    short = [f"{a} {w}" for a, st in (("P0", sp0), ("P3", sp3)) for w in WIN if st[w]["n_finite"] < MIN_PTS]
    if short:
        rec["short"] = short
        rec["worsening"] = None
        return rec
    if any(sp0[w]["median"] == 0.0 for w in WIN):
        rec["p0_median_zero"] = True
        rec["worsening"] = False
        return rec
    r1 = sp3["W1"]["median"] / sp0["W1"]["median"]; r2 = sp3["W2"]["median"] / sp0["W2"]["median"]
    persist = sp3["W2"]["median"] / sp3["W1"]["median"] if sp3["W1"]["median"] > 0 else None
    rec.update(ratio_W1=r1, ratio_W2=r2, P3_W2_over_W1=persist)
    rec["worsening"] = bool(r1 >= RATIO and r2 >= RATIO and persist is not None and persist >= PERSIST)
    return rec


def timing(run):
    d = HERE / run
    try:
        t0 = float((d / "START_EPOCH").read_text().strip())
        rows = [r for r in csv.reader(open(d / "timing.tsv"), delimiter="\t") if len(r) >= 3]
    except (OSError, ValueError) as e:
        return {"error": f"時刻の記録を読めない ({e})"}
    nf = [int(r[2]) for r in rows]
    smp = [(float(r[0]), int(r[1])) for r in rows if r[1].strip()]
    out = {"samples": len(rows), "exclusive": bool(rows) and all(n <= 1 for n in nf), "max_forge": max(nf) if nf else None,
           "others": sorted({r[3].strip() for r in rows if len(r) > 3 and r[3].strip()})}
    if len(smp) >= 2:
        # 1 step の時間 = step 2,500 以上の最初の標本から、最後の step に初めて達した標本までの傾き (起動を除く。forge の終了後も
        # RUN_RC ができるまで標本が続くので、最後の標本でなく最後の step に最初に達した標本を終点にする)
        j = next((k for k, (_, s) in enumerate(smp) if s >= 2500), None)
        s_end = max(s for _, s in smp)
        t_end = next(t for t, s in smp if s == s_end)
        if j is not None and s_end > smp[j][1]:
            out["per_step_ms"] = 1e3 * (t_end - smp[j][0]) / (s_end - smp[j][1])
            out["per_step_span"] = [smp[j][1], s_end]
        out["last_step"] = s_end; out["last_step_reached_s"] = t_end - t0
    fl = d / "forge_run.log"
    if fl.exists():
        out["T_forge_end_s (forge_run.log の更新時刻 − 起動)"] = fl.stat().st_mtime - t0
    return out


def fields(run, step):   # sp_judge.py の fields と同じ (§6.58 の gg_judge.py と同じ式)
    yb_x, yb = XC.common_yb()
    xy, ro, ux, uy, T, k, h0 = XC.load_forge(HERE / run, step, with_h0=True)
    ni, nj, S = XC.mesh_info(len(ro))
    o = XC.reduce_fields(xy, ro, ux, uy, T, k, False, ni, nj, S, yb_x, yb, h0=h0)
    with h5py.File(HERE / run / f"res_{step}.h5", "r") as h:
        M = (np.hypot(np.asarray(h["VALUE/Ux"][:], float), np.asarray(h["VALUE/Uy"][:], float)) / np.asarray(h["VALUE/sonic"][:], float)).reshape(ni, nj)
    r = (xy[:, 1] / S).reshape(ni, nj)
    rr = r[-1] / r[-1][-1]
    m80 = rr <= 0.8; m80[0] = False; m95 = rr <= 0.95; m95[0] = False
    rec = {"Q_w": float(o["Q_w"]), "M_exit_mean_r_le_0.8 (軸の行を除く)": float(np.mean(M[-1][m80])),
           "M_exit_max_r_le_0.95 (軸の行を除く)": float(np.max(M[-1][m95])), "M_exit_axis": float(M[-1][0])}
    for xs in XC.XS:
        i = int(np.argmin(np.abs(o["x"] - xs)))
        rec[f"theta_r_{int(xs)}"] = float(o["theta_r"][i]); rec[f"delta_loc_{int(xs)}"] = float(o["delta_loc"][i])
    return rec


def conv(run):
    tool = TOOLS / "check_convergence.py"
    if not tool.exists():
        return f"check_convergence.py が無い ({tool})"
    p = subprocess.run([sys.executable, str(tool), str(HERE / run)], capture_output=True, text=True)
    lines = (p.stdout + p.stderr).splitlines()
    return {"overall": next((l for l in lines if "OVERALL" in l), "(行なし)"), "failed": [l.strip()[:120] for l in lines if "<--" in l]}


def main():
    OUT = {"arms": ARMS, "windows": {w: list(v) for w, v in WIN.items()}, "undecidable": [], "global": {}, "local": {}, "records": {}}
    und = OUT["undecidable"]
    # 資格
    cp = OUTD / "tp_check.json"
    chk = json.loads(cp.read_text()) if cp.exists() else None
    OUT["qualification"] = None if chk is None else {"status": chk.get("status"), "passes": chk.get("passes"), "failed": chk.get("failed"),
                                                     "undecidable": chk.get("undecidable")}
    for nm, run in ARMS.items():
        OUT["records"][f"{nm} timing"] = timing(run)
        p = HERE / run / "RUN_RC"
        OUT["records"][f"{nm} RUN_RC"] = p.read_text().strip() if p.exists() else None
    verdict = None
    if chk is None:
        verdict = "判定不能: tp_check.json が無い (資格を確かめていない)"
    elif chk.get("status") == "不成立":
        verdict = f"資格不成立 (速さの試験の結果にしない): 不合格 {chk.get('failed')}" + (f"、判定不能 {chk.get('undecidable')}" if chk.get("undecidable") else "")
    elif chk.get("status") != "成立":
        verdict = f"判定不能 (資格を確かめられない): 判定不能 {chk.get('undecidable')}"
    if verdict is None:
        # 全域の指標
        try:
            g = {nm: load_global(run) for nm, run in ARMS.items()}
            c0, c3 = g["P0"][0], g["P3"][0]
            if c0 != c3:
                und.append(f"残差の列が腕で違う (P0 {c0}、P3 {c3})")
            for nm in ARMS:
                if g[nm][2] is None or g[nm][2] + 1 < END_STEP:
                    und.append(f"{nm}: 残差の履歴が {END_STEP} step に届かない (最後の step {g[nm][2]})")
            for c in c0:
                if c in c3:
                    OUT["global"][c] = judge_indicator(win_stats(g["P0"][1][c]), win_stats(g["P3"][1][c]))
        except Missing as e:
            und.append(str(e))
        # 局所の指標
        try:
            fixed = json.loads((OUTD / "tp_region_set.json").read_text())
            loc = {nm: load_local(run, fixed) for nm, run in ARMS.items()}
            for q in LOCAL_Q:
                OUT["local"][f"Σ|res_{q}| (x/r_t ∈ [−8, −2))"] = judge_indicator(win_stats(loc["P0"][0][q]), win_stats(loc["P3"][0][q]))
            OUT["records"]["局所の出力の step (P0・P3 の数)"] = [len(loc["P0"][1]), len(loc["P3"][1])]
        except (Missing, OSError, KeyError) as e:
            und.append(str(e))
        allind = {f"全域 {k}": v for k, v in OUT["global"].items()} | {f"局所 {k}": v for k, v in OUT["local"].items()}
        short = [f"{k} ({', '.join(v['short'])})" for k, v in allind.items() if v["worsening"] is None]
        if short:
            und.append("窓の点数が足りない (有限の点 5 点未満): " + "; ".join(short[:6]) + (f" ほか {len(short) - 6}" if len(short) > 6 else ""))
        OUT["p0_median_zero (比を作らず別記)"] = [k for k, v in allind.items() if v.get("p0_median_zero")]
        OUT["nonfinite_in_windows (比の判定から分けた)"] = {k: {a: {w: v[a][w]["n_nonfinite"] for w in WIN} for a in ("P0", "P3")}
                                                    for k, v in allind.items() if any(v[a][w]["n_nonfinite"] for a in ("P0", "P3") for w in WIN)}
        OUT["worsening"] = [k for k, v in allind.items() if v["worsening"]]
        if und:
            verdict = "判定不能: " + "; ".join(und[:4])
        elif OUT["worsening"]:
            verdict = ("結果 A (持続する悪化あり: " + ", ".join(OUT["worsening"]) + ") → この設定の長い試験を打ち切る"
                       " (この条件では継続の対象から外した。低マッハの前処理は効かないとは言わない)")
        else:
            verdict = "結果 B (どの指標にも持続する悪化なし) → 段 2 へ (別に登録の追記をしてから。P3 が優れるとは言わない)"
    # 記録: step 20,000 の量 (未収束の値として) と check_convergence
    F = {}
    for nm, run in ARMS.items():
        try:
            F[nm] = fields(run, END_STEP); OUT["records"][f"{nm} step {END_STEP} (未収束の値)"] = F[nm]
        except (Exception, SystemExit) as e:  # noqa: BLE001  記録だけ (判定に使わない) なので読めなければ理由を書いて先へ
            OUT["records"][f"{nm} step {END_STEP} (未収束の値)"] = f"読めない ({type(e).__name__}: {e})"
        try:
            OUT["records"][f"{nm} check_convergence"] = conv(run)
        except OSError as e:
            OUT["records"][f"{nm} check_convergence"] = f"実行できない ({e})"
    if len(F) == 2:
        OUT["records"][f"P3 − P0 の相対差 [%] (step {END_STEP}、記録)"] = {k: 100 * (F["P3"][k] - F["P0"][k]) / abs(F["P0"][k])
                                                                        for k in F["P0"] if k.startswith(("theta_r", "Q_w")) and F["P0"][k]}
    OUT["verdict"] = verdict
    for k, v in list(OUT["global"].items()) + list(OUT["local"].items()):
        if "ratio_W1" in v:
            print(f"{k:40s} 比 W1 {v['ratio_W1']:.3g}  W2 {v['ratio_W2']:.3g}  P3 W2/W1 {v['P3_W2_over_W1']}  悪化 {v['worsening']}")
        else:
            print(f"{k:40s} {'P0 の中央値 0 (比を作らない)' if v.get('p0_median_zero') else '点数不足 ' + str(v.get('short'))}")
    print("時間", json.dumps({k: v for k, v in OUT["records"].items() if "timing" in k}, ensure_ascii=False, default=float))
    print(f"== VERDICT §6.62 段 1: {verdict}")
    OUTD.mkdir(parents=True, exist_ok=True)
    (OUTD / "tp_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float))
    print("→ _band_ab/cold_pair/tp_judge.json")
    return OUT


if __name__ == "__main__":
    main()
