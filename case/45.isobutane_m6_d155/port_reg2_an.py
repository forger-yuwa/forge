"""plan architecture-float-state-double-geometry §6.60 (事前登録、codex diagnose `2026-10-13-port-reg-redo-diagnose.md`) の判定: 移植 (§5.1 #26) の回帰 R1・R3 の回し直し。
各組 (R1: case/45 の point・B0 × float・FP64、R3: case/48 の平板 × float・FP64) で、移植の前 3 本 (A 群) と後 3 本 (B 群)、同じ入力から 10 step。A で GPU を占有し、前・後を交互に。
量: res_10.h5 の VALUE の全浮動小数点の量 (節点ごと) と residual_history.csv の全 rms_* 列 (step ごと)。
  D_q = max |mean_B − mean_A|、E_A,q = max (A 群の中の最大 − 最小)、E_B,q = 同じく B 群 (どれも節点ごと・step ごとの最大)。
許容 τ_q (結果を見る前に固定、精度の理屈から。左辺は両方のビルドで単精度で組まれるので両ビルドに同じ値):
  状態・原始量 (下の DIAG 以外の VALUE の量): τ = 1e-4 × max|mean_A| + 1e-30。
  打ち消し合いの残りの量 (res_*・omg_*・dt_local): τ = 1e-1 × max|mean_A| + 1e-30。rms_* 列: step ごとに τ = 1e-2 × |mean_A| + 1e-30 (その最小値で判定)。
入力の検査 (欠ければ判定不能): 全 run の RUN_RC 0、res_10.h5・residual_history.csv がある、全量が有限、6 本で VALUE の量の集合・形・rms 列・step 列が一致、空でない。
結果 A (この入力・10 step・登録の許容の範囲で回帰合格): 全組・全量で D + E_A + E_B ≤ τ。
結果 B (基準だけでは説明できない、実装の差を調べる): いずれかの量で D > τ + E_A + E_B。それ以外は判別不能。
→ _band_ab/cold_pair/port_reg2_judge.json"""
import csv, json, math, re, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent; C48 = HERE.parent / "48.flat_plate_cooled_m4"
DIAG = re.compile(r"^(res_|omg_|dt_local)")
OUT = {"groups": {}, "undecidable": []}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
def load(d):
    rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    if rc != "0": raise ValueError(f"{d.name}: RUN_RC {rc}")
    p = d / "res_10.h5"
    if not p.exists(): raise ValueError(f"{d.name}: res_10.h5 が無い")
    with h5py.File(p, "r") as h: V = {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in h["VALUE"] if h["VALUE"][k].dtype.kind == "f"}
    rows = list(csv.DictReader(open(d / "residual_history.csv")))
    if not rows: raise ValueError(f"{d.name}: 残差の履歴が空")
    cols = sorted(c for c in rows[0] if c.startswith("rms_"))
    R = {c: np.array([float(x[c]) for x in rows]) for c in cols}
    steps = np.array([int(x["step"]) for x in rows])
    for k, a in list(V.items()) + list(R.items()):
        if a.size == 0: raise ValueError(f"{d.name}: {k} が空")
        if not np.all(np.isfinite(a)): raise ValueError(f"{d.name}: {k} に非有限")
    return V, R, steps
def judge(base, pre, post):
    try:
        A = [load(base / r) for r in pre]; B = [load(base / r) for r in post]
    except (ValueError, OSError, KeyError) as e:
        return None, str(e)
    keysV = set(A[0][0]); keysR = set(A[0][1])
    for V, R, st in A + B:
        if set(V) != keysV or set(R) != keysR: return None, "量の集合が run の間で違う"
        if any(V[k].shape != A[0][0][k].shape for k in keysV) or not np.array_equal(st, A[0][2]): return None, "形か step の列が run の間で違う"
    res = {"items": {}, "A_ok": True, "B_hit": []}
    for src, keys, idx in (("場", sorted(keysV), 0), ("残差", sorted(keysR), 1)):
        for k in keys:
            a = np.stack([x[idx][k] for x in A]); b = np.stack([x[idx][k] for x in B])
            mA, mB = a.mean(0), b.mean(0)
            D = float(np.max(np.abs(mB - mA))); EA = float(np.max(a.max(0) - a.min(0))); EB = float(np.max(b.max(0) - b.min(0)))
            if src == "場":
                rel = 1e-1 if DIAG.match(k) else 1e-4; tau = rel * float(np.max(np.abs(mA))) + 1e-30
            else:
                tau = float(np.min(1e-2 * np.abs(mA))) + 1e-30; rel = 1e-2
            res["items"][f"{src}|{k}"] = dict(D=D, EA=EA, EB=EB, tau=tau, rel=rel)
            if D + EA + EB > tau: res["A_ok"] = False
            if D > tau + EA + EB: res["B_hit"].append(f"{src} {k}: D {D:.3e} > τ {tau:.3e} + E {EA + EB:.3e}")
    return res, None
groups = []
for tag in ("pt", "b0"):
    for p in ("f32", "fp64"):
        groups.append((f"R1 {tag} {p}", HERE, sorted(x.name for x in HERE.glob(f"run_0*_q1_{tag}_{p}_pre*") if x.is_dir()), sorted(x.name for x in HERE.glob(f"run_0*_q1_{tag}_{p}_post*") if x.is_dir())))
for p in ("f32", "fp64"):
    groups.append((f"R3 {p}", C48, sorted(x.name for x in C48.glob(f"run_0*_q3_{p}_pre*") if x.is_dir()), sorted(x.name for x in C48.glob(f"run_0*_q3_{p}_post*") if x.is_dir())))
verdicts = []
for nm, base, pre, post in groups:
    if len(pre) != 3 or len(post) != 3: und(f"{nm}: run がそろわない (前 {len(pre)}・後 {len(post)})"); continue
    r, e = judge(base, pre, post)
    if e: und(f"{nm}: {e}"); continue
    OUT["groups"][nm] = {"pre": pre, "post": post, **r}
    worst = sorted(r["items"].items(), key=lambda kv: -(kv[1]["D"] + kv[1]["EA"] + kv[1]["EB"]) / kv[1]["tau"])[:3]
    print(nm, "A" if r["A_ok"] else ("B" if r["B_hit"] else "判別不能"), [(k, f"(D+E)/τ {(v['D'] + v['EA'] + v['EB']) / v['tau']:.2e}") for k, v in worst])
    verdicts.append("A" if r["A_ok"] else ("B" if r["B_hit"] else "中間"))
if OUT["undecidable"]: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:3])
elif all(v == "A" for v in verdicts): verdict = "結果 A (全組・全量で D + E_A + E_B ≤ τ、この入力・10 step・登録の許容の範囲で回帰合格)"
elif any(v == "B" for v in verdicts): verdict = "結果 B (D > τ + E の量がある、実装の差を調べる)"
else: verdict = "判別不能"
OUT["verdict"] = verdict
print("== VERDICT §6.60:", verdict)
(HERE / "_band_ab/cold_pair/port_reg2_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/port_reg2_judge.json")
