#!/usr/bin/env python3
"""plan architecture-float-state-double-geometry §6.21 (事前登録、codex diagnose 2026-10-10) の判定: case/48 の float と FP64 の比較。
c48_prec.sh の後に AWS の case/48 で回す。結果は標準出力と case/48/c48_prec_judge.json。

- 腕: float = run_0057_pv_f32a・run_0058_pv_f32b、FP64 = run_0059_pv_f64a・run_0060_pv_f64b (同じソース・同じ double の幾何・同じ初期保存量・同じ実効設定)。
- 前提 (満たさなければ判定不能): 4 本とも 48,000 step を完走し、残差の全期間・全列が有限、prec_series.csv (c48_prec_series.py) が 24 点そろって有限。
  正式ツールが実行できない (VERDICT の行が無い) のも判定不能。check_convergence の NOT CONVERGED は記録する (判定不能にはしない)。
- 窓: 24,000〜48,000 step の 13 点。各 run の窓を check_quasisteady --series-csv --tail 1.0 で判定 (厚さ θ・δ* は drift/osc 1e-4、Cf・q_w・CD・HF は 2e-4)。
  4 本すべてで全量が STEADY でなければ判別不能。
- 判定: 各量で μ32・μ64 = 2 本・全窓点の平均、D = |μ32 − μ64|/|μ64|、b = 各精度内で平均からの最大偏差、E = (b32 + b64)/|μ64|。
  τ = 0.0005 (θ・δ*)、0.001 (Cf・q_w・CD・HF)。全量で D + E ≤ τ → 「登録窓の対象量は許容内」。いずれかで D − E > τ → 「許容外の精度差を検出」。
  それ以外は判別不能。D − E ≥ 0.005 は「0.5 % 以上の明瞭な差」と追加表示する。観測した 4 本の変動幅による判定で、統計的信頼区間ではない。"""
import csv, json, math, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parents[1]
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
ARMS = {"f32": ("run_0057_pv_f32a", "run_0058_pv_f32b"), "f64": ("run_0059_pv_f64a", "run_0060_pv_f64b")}
XS = (0.3, 0.6, 0.9)
THICK = [f"{q}_{x}" for x in XS for q in ("theta", "dstar")]
OTHER = [f"{q}_{x}" for x in XS for q in ("Cf", "qw")] + ["CD", "HF"]
TAU = {**{k: 5e-4 for k in THICK}, **{k: 1e-3 for k in OTHER}}
NEED = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega")
W0, W1 = 24000, 48000
OUT = {"undecidable": [], "runs": {}}
def und(msg): OUT["undecidable"].append(msg); print("  [判定不能]", msg)
series = {}
for arm, runs in ARMS.items():
    for r in runs:
        run = HERE / r; rec = OUT["runs"].setdefault(r, {})
        # 完走と残差の全期間
        rh = run / "residual_history.csv"
        if not rh.exists(): und(f"{r}: residual_history.csv が無い"); continue
        hr = [x for x in csv.DictReader(open(rh)) if x.get("phase", "outer_end") == "outer_end"]
        cols = [c for c in (hr[0] if hr else {}) if c.startswith("rms_") and not c.startswith("rms_dq")]
        miss = [c for c in NEED if c not in cols]
        if miss: und(f"{r}: 残差の列が無い {miss}")
        last = int(hr[-1]["step"]) if hr else -1; rec["last_step"] = last
        if last + 1 < W1: und(f"{r}: 完走していない (最後の step {last})")
        nf = [(x["step"], c) for x in hr for c in cols if not math.isfinite(float(x[c]))]
        if nf: und(f"{r}: 残差に非有限 (最初 {nf[0]})")
        # 判定量の系列
        sc = run / "prec_series.csv"
        if not sc.exists(): und(f"{r}: prec_series.csv が無い"); continue
        S = list(csv.DictReader(open(sc)))
        steps = [int(x["step"]) for x in S]
        if steps != list(range(2000, W1 + 1, 2000)): und(f"{r}: 系列の step がそろわない ({len(steps)} 点)"); continue
        win = [x for x in S if W0 <= int(x["step"]) <= W1]
        vals = {k: np.array([float(x[k]) for x in win]) for k in THICK + OTHER}
        if any(not np.all(np.isfinite(v)) for v in vals.values()): und(f"{r}: 判定量に非有限"); continue
        series[r] = vals
        rec["je_range"] = {f"je_{x}": [int(min(float(w[f'je_{x}']) for w in win)), int(max(float(w[f'je_{x}']) for w in win))] for x in XS}
        wc = run / "prec_window.csv"
        with open(wc, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["step"] + THICK + OTHER); w.writeheader()
            for x in win: w.writerow({k: x[k] for k in ["step"] + THICK + OTHER})
        # 準定常 (正式ツール。run_dir を渡すと組み込みの量まで判定されるので渡さない)
        steady = True
        for grp, tol in ((THICK, "0.0001"), (OTHER, "0.0002")):
            p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(wc), "--series-cols", ",".join(grp),
                                "--tail", "1.0", "--drift", tol, "--osc", tol], capture_output=True, text=True)
            ov = [l for l in (p.stdout + p.stderr).splitlines() if l.startswith("OVERALL")]
            if not ov or p.returncode not in (0, 1): und(f"{r}: check_quasisteady が実行できない (rc {p.returncode})"); steady = False; continue
            rec[f"quasisteady_{tol}"] = ov[-1]; rec[f"quasisteady_{tol}_detail"] = [l for l in p.stdout.splitlines() if "STEADY" in l or "DRIFT" in l or "OSC" in l or "単調" in l][:20]
            if p.returncode != 0: steady = False
        rec["all_steady"] = steady
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        rec["check_convergence"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
        log = (run / "forge_run.log").read_text(errors="replace") if (run / "forge_run.log").exists() else ""
        rec["settings_echo"] = sorted({l.strip()[:160] for l in log.splitlines() if any(k in l for k in ("scalarGradient", "slauWallNormalChi", "limiterScaled", "venkatK", "implicitSolvePrecision"))})[:12]
        print(f"== {r}: 最後の step {last}、準定常 {'ALL STEADY' if steady else 'NOT ALL STEADY'}、check_convergence: {rec['check_convergence']}")
all_present = all(r in series for runs in ARMS.values() for r in runs)
if not all_present and not OUT["undecidable"]: und("系列のそろわない run がある")
res = {}
if all_present:
    for k in THICK + OTHER:
        a32 = np.concatenate([series[r][k] for r in ARMS["f32"]]); a64 = np.concatenate([series[r][k] for r in ARMS["f64"]])
        m32, m64 = a32.mean(), a64.mean(); b32, b64 = np.max(np.abs(a32 - m32)), np.max(np.abs(a64 - m64))
        D = abs(m32 - m64) / abs(m64); Ev = (b32 + b64) / abs(m64)
        res[k] = dict(mu32=float(m32), mu64=float(m64), D=float(D), E=float(Ev), tau=TAU[k],
                      cls="許容内" if D + Ev <= TAU[k] else ("許容外" if D - Ev > TAU[k] else "判別不能"))
        print(f"  {k:10s} μ32 {m32:.6e}  μ64 {m64:.6e}  D {D:.2e}  E {Ev:.2e}  τ {TAU[k]:.0e}  → {res[k]['cls']}{'  (0.5 % 以上の明瞭な差)' if D - Ev >= 0.005 else ''}")
OUT["quantities"] = res
steady_all = all(OUT["runs"].get(r, {}).get("all_steady") for runs in ARMS.values() for r in runs)
if OUT["undecidable"]:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"])
elif not steady_all:
    verdict = "判別不能 (窓の対象量が 4 本すべてで STEADY でない)"
elif any(v["cls"] == "許容外" for v in res.values()):
    verdict = "許容外の精度差を検出: " + ", ".join(f"{k} (D {v['D']:.2e}, E {v['E']:.2e})" for k, v in res.items() if v["cls"] == "許容外")
    if any(v["D"] - v["E"] >= 0.005 for v in res.values()): verdict += " — 0.5 % 以上の明瞭な差あり"
elif all(v["cls"] == "許容内" for v in res.values()):
    verdict = "登録窓の対象量は許容内"
else:
    verdict = "判別不能: " + ", ".join(k for k, v in res.items() if v["cls"] == "判別不能")
OUT["verdict"] = verdict
print(f"== VERDICT §6.21: {verdict}")
(HERE / "c48_prec_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float))
print("→", HERE / "c48_prec_judge.json")
