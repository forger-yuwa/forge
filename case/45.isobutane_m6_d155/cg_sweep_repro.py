"""plan tooling-nozzle-core-grid §4.44 M-d の再試験 (事前登録、codex 2026-10-11 dtedge): 上限 50/1 を同じ条件で各 1 回だけ回し直し、
同じ腕の中の rhs_s0 の不一致と u5 の再現性を確かめる。
入力: 4 run (cap50・cap1 = 既存 run_0628・run_0629、cap50b・cap1b = 再試験) の linedump と、cg_sweep_an.py の sweep_<tag>.json。
見る量: 状態・節点の一致 (sha256)、rhs_s0 の成分別の差 ‖Δ(rhs/V)‖/‖rhs/V‖ (体積重み付き、拘束の節点を除き、ライン全体と島で別々)。
  組: 同じ腕 (cap50 対 cap50b、cap1 対 cap1b) と腕の間 (cap50 対 cap1、cap50b 対 cap1b)。
許容差 (今回の差 (最大 7.0e−9) を見た後に決めた基準であることを記録する): 成分別に 1e−8。
判定: R-A — 同じ腕でも 0 でない rhs_s0 の差が再現し、同じ腕・腕の間のすべての差が許容差内、再試験の 2 本で資格 (2)(3) と島の u5 ≤ 0.1
       (ρ・ρE・全 5 成分) が再現 → 「設定の差が無ければ rhs_s0 は一致する」を退け、非決定的な数値の差という説明を支持。
       §4.42 の再試験の結果 B (この右辺に対する大きな内部の解き残しを退ける) として記録する。
      R-B — 許容差を超える、または u5・資格 (2)(3) が再現しない → 「小さな入力の差だから判定に影響しない」という扱いを退け、入力・書き出しの経路の診断を先にする。
      同じ腕の 2 本がビット一致しただけ (同じ腕の差が 0) → 非決定性は確かめられていない (判定不能、原因を二択で断定しない)。
usage (AWS): python3 cg_sweep_repro.py <run50> <run1> <run50b> <run1b> → _band_ab/core_grid/sweep_repro.json"""
import hashlib, json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
runs = dict(zip(("cap50", "cap1", "cap50b", "cap1b"), sys.argv[1:5]))
TOL = 1e-8
def get(run, name):
    meta = {}
    for ln in (HERE / run / "linedump" / "meta.txt").read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        n, r, c = ln.split(); meta[n] = (int(r), int(c))
    r, c = meta[name]
    return np.fromfile(HERE / run / "linedump" / f"{name}.f64", dtype=np.float64).reshape(r, c)
r0 = runs["cap50"]
nl = get(r0, "node_line"); nodes = nl[:, 0].astype(np.int64)
info = json.loads((HERE / r0 / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
with h5py.File(HERE / r0 / "res_0.h5", "r") as h:
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
nj = len(xy) // ni; x = xy[nodes, 0] / S; jw = nj - 1 - (nodes % nj)
V = get(r0, "dt_vol")[:, 1]; con = np.any(get(r0, "flags_wall_iso_axis") != 0, axis=1)
reg = {"lines": ~con, "island": (x >= -6.0) & (x <= -3.25) & (jw >= 60) & (jw <= 100) & ~con}
D = {k: {"rhs": get(r, "rhs_s0"), "state": hashlib.sha256(get(r, "state_ro_roU_roe_cp_gamma").tobytes()).hexdigest(),
         "nodes": hashlib.sha256(get(r, "node_line")[:, 0].tobytes()).hexdigest()} for k, r in runs.items()}
def diff(a, b):
    o = {}
    for rn, m in reg.items():
        o[rn] = {}
        for c in range(5):
            ref = np.sqrt(np.sum(V[m] * (a[m, c] / V[m]) ** 2))
            if ref == 0:
                continue
            o[rn][str(c)] = float(np.sqrt(np.sum(V[m] * ((a[m, c] - b[m, c]) / V[m]) ** 2)) / ref)
    o["rows_differ"] = int(np.any(a != b, axis=1).sum())
    return o
pairs = {"same_cap50": ("cap50", "cap50b"), "same_cap1": ("cap1", "cap1b"), "cross_orig": ("cap50", "cap1"), "cross_retest": ("cap50b", "cap1b")}
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.44", "runs": runs, "tol": TOL, "tol_note": "今回の差 (最大 7.0e-9) を見た後に決めた基準",
       "same_state_nodes": all(D[k]["state"] == D["cap50"]["state"] and D[k]["nodes"] == D["cap50"]["nodes"] for k in D), "pairs": {}}
for pn, (a, b) in pairs.items():
    out["pairs"][pn] = diff(D[a]["rhs"], D[b]["rhs"])
mx = max(v for p in out["pairs"].values() for rn in ("lines", "island") for v in p[rn].values())
same_nonzero = out["pairs"]["same_cap50"]["rows_differ"] > 0 or out["pairs"]["same_cap1"]["rows_differ"] > 0
J = {t: json.loads((OUT / f"sweep_{t}.json").read_text()) for t in ("cap50b", "cap1b")}
qual = all(X["nonfinite"] == 0 and X["backward_error_max"] <= 1e-8 and all(
    X["region"][r][c]["g40_over_b0"] <= 1e-3 and X["region"][r][c]["conv_40_30"] <= 1e-3 for r in ("island", "lines") for c in ("ro", "roe", "all5"))
    for X in J.values())
u5ok = all(X["region"]["island"][c]["u5"] <= 0.1 for X in J.values() for c in ("ro", "roe", "all5"))
out.update({"max_rel_diff": mx, "same_arm_nonzero": same_nonzero, "retest_qual_2_3": qual, "retest_u5_le_0.1": u5ok,
            "retest_u5_island": {t: {c: X["region"]["island"][c]["u5"] for c in ("ro", "roe", "all5")} for t, X in J.items()}})
if not out["same_state_nodes"] or mx > TOL or not qual or not u5ok:
    res = "R-B: 許容差を超える、または資格 (2)(3)・u5 が再現しない → 入力・書き出しの経路の診断を先にする"
elif not same_nonzero:
    res = "判定不能: 同じ腕の 2 本がビット一致 (非決定性は確かめられていない)"
else:
    res = "R-A: 同じ腕でも rhs_s0 の丸めの差が再現し、すべて許容差内、資格と u5 も再現 → §4.42 の再試験の結果 B (大きな内部の解き残しを退ける)"
out["result"] = res
(OUT / "sweep_repro.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for pn, p in out["pairs"].items():
    print(pn, "行", p["rows_differ"], {rn: {c: f"{v:.1e}" for c, v in p[rn].items()} for rn in ("lines", "island")})
print("状態・節点の一致", out["same_state_nodes"], " 資格 (2)(3)", qual, " u5 ≤ 0.1", u5ok, " 再試験の島の u5", out["retest_u5_island"])
print("結果:", res)
