"""plan tooling-nozzle-core-grid §4.86 の材料 (探索的・登録外、判定に使わない): point の前処理の腕 (=3・=2) で縮流部の残差がどこで増えているか。
各 run の res_20000 (C・P3・P2) から:
 (a) 単位体積あたりの残差 R/V (ρ・ρE) の二乗の体積重み付き和 Σ V (R/V)² の、x/r_t と壁からの節点番号 jw の区間ごとの割合と大きさ。
 (b) 寄与の上位 0.1 % の節点の x/r_t・jw の中央値と範囲。
 (c) dt_local の比 (腕/C) の区間ごとの中央値 (前処理の Δτ′ の拡大がどこで効いているか)。
 (d) 場の差 (腕 − C) の成分: 区間ごとの体積重み付き RMS で、エントロピーの部分 δρ − δp/c² と圧力の部分 δp/c² (c は C の音速) の比。
usage: python3 cg_pp_where.py [C P3 P2 の run] → _band_ab/core_grid/pp_where.json"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
RUNS = sys.argv[1:4] if len(sys.argv) >= 4 else ["run_0688_pp_c", "run_0689_pp_p3", "run_0690_pp_p2"]; N = 20000
XB = [(-99, -8), (-8, -6), (-6, -4), (-4, -2), (-2, 0), (0, 99)]; JB = [(0, 10), (10, 40), (40, 60), (60, 100), (100, 160)]

def load(run):
    info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
    with h5py.File(HERE / run / f"res_{N}.h5", "r") as h:
        v = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("res_ro", "res_roe", "volume", "dt_local", "ro", "P", "sonic", "T")}
        x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
    nj = x.size // ni; idx = np.arange(x.size); jw = nj - 1 - idx % nj
    return v, x, jw

D = {r: load(r) for r in RUNS}
x, jw = D[RUNS[0]][1], D[RUNS[0]][2]
for r in RUNS[1:]:
    if not np.array_equal(D[r][1], x): raise SystemExit("座標が違う")
out = {"runs": RUNS, "step": N, "bins": {}}
for r in RUNS:
    v = D[r][0]; rec = {}
    for k in ("res_ro", "res_roe"):
        e = v["volume"] * (v[k] / v["volume"]) ** 2; tot = e.sum()
        rec[k] = {"total": float(tot), "share": {f"x{a}..{b}|jw{c}..{d}": float(e[(x >= a) & (x < b) & (jw >= c) & (jw < d)].sum() / tot) for a, b in XB for c, d in JB}}
        top = np.argsort(-e)[: max(1, e.size // 1000)]
        rec[k]["top0.1%"] = {"x_median": float(np.median(x[top])), "x_p10_p90": [float(np.percentile(x[top], 10)), float(np.percentile(x[top], 90))],
                             "jw_median": float(np.median(jw[top])), "jw_p10_p90": [float(np.percentile(jw[top], 10)), float(np.percentile(jw[top], 90))]}
    out["bins"][r] = rec
c0 = D[RUNS[0]][0]
out["dt_ratio"] = {}; out["diff_components"] = {}
for r in RUNS[1:]:
    v = D[r][0]; rat = v["dt_local"] / c0["dt_local"]; dp = v["P"] - c0["P"]; dro = v["ro"] - c0["ro"]; c2 = c0["sonic"] ** 2; V = c0["volume"]
    out["dt_ratio"][r] = {}; out["diff_components"][r] = {}
    for a, b in XB:
        for c, d in JB:
            m = (x >= a) & (x < b) & (jw >= c) & (jw < d)
            if not m.any(): continue
            key = f"x{a}..{b}|jw{c}..{d}"
            out["dt_ratio"][r][key] = float(np.median(rat[m]))
            rms = lambda z: float(np.sqrt(np.sum(V[m] * z[m] ** 2) / np.sum(V[m])))
            ent, prs = rms(dro - dp / c2), rms(dp / c2)
            out["diff_components"][r][key] = {"entropy": ent, "pressure": prs, "A": prs / ent if ent > 0 else None, "rms_dT": rms(v["T"] - c0["T"])}
(OUTD / "pp_where.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for r in RUNS:
    for k in ("res_ro", "res_roe"):
        b = out["bins"][r][k]; big = sorted(b["share"].items(), key=lambda t: -t[1])[:5]
        print(f"{r} {k} 合計 {b['total']:.3e}  上位の区間: " + ", ".join(f"{kk} {vv:.2f}" for kk, vv in big))
        print(f"    上位 0.1 % の節点: x 中央値 {b['top0.1%']['x_median']:.2f} ({b['top0.1%']['x_p10_p90'][0]:.2f}〜{b['top0.1%']['x_p10_p90'][1]:.2f})  jw 中央値 {b['top0.1%']['jw_median']:.0f} ({b['top0.1%']['jw_p10_p90'][0]:.0f}〜{b['top0.1%']['jw_p10_p90'][1]:.0f})")
for r in RUNS[1:]:
    print(f"{r} dt 比 (腕/C) 区間ごと: " + ", ".join(f"{k} {v:.2f}" for k, v in out["dt_ratio"][r].items() if "x-8..-6" in k or "x-6..-4" in k or "x-4..-2" in k))
    print(f"{r} 場の差の成分 (縮流部): " + ", ".join(f"{k} A {v['A']:.3g} dT {v['rms_dT']:.2g}" for k, v in out["diff_components"][r].items() if "x-6..-4" in k or "x-4..-2" in k))
