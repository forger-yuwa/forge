"""plan tooling-nozzle-core-grid §4.78 (事前登録): 停滞帯の圧力の鋸歯と双子の差の成分。run は回さず既存の場を読む。
線: 構造格子の添字 (列 i、壁からの節点番号 jw = nj−1−jj、cg_twin_diff.py と同じ)。島 = 壁の列の x/r_t ∈ [−6, −3.25]・jw 60〜100。
  壁法線の線 = 島の各列の jw 58〜102、流れ方向の線 = jw 60〜100 の各行の島の列 ± 2 列 (4 階差分の中心が島に収まる)。
指標 (線ごと): a = RMS(δ⁴φ)/16、b = RMS(δ¹φ)、η = a/b、f = δ²φ の符号が隣どうしで反転する割合。P は a/q (q = ½ρ|u|²、線の中央値) も。
資格: 非有限なし、線の弧長 s(k) の η の中央値 ≤ 0.01。
読み方 (圧力の鋸歯): S = Gref の B0 の 2 枚とも、少なくとも 1 方向で P の f の中央値 ≥ 0.7・η の中央値 ≥ 0.1・|u| の η ≤ P の η/3。
  N = Gref の B0 の 2 枚とも両方向で P の η の中央値 ≤ 0.01。それ以外は判定不能。上限 1 と Grefc は記述だけ。
双子の差の成分: 島の体積重み付き RMS で A = ‖δp/c²‖ / ‖δρ − δp/c²‖ (c は基準側の sonic)。A ≤ 0.1 → エントロピーの成分が主、A ≥ 1 → 圧力の成分が主、その間は混在。
usage: python3 cg_lowmach_saw.py → _band_ab/core_grid/lowmach_saw.json"""
import json
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
FIELDS = {"gref_b0_0608": ("run_0608_cg_gref_cont_k1", 60000), "gref_b0_0619": ("run_0619_cg_gref_b0_ctl", 60000),
          "gref_cap1_0616": ("run_0616_cg_gref_cap1", 20000), "grefc_b0_0625": ("run_0625_cg_grefc", 60000)}
TWINS = {"gref_twin_3000": ("run_0679_tw50r_base", "run_0680_tw50r_pert", 3000)}
for n in (1000, 2000, 3000):
    TWINS[f"grefc_ZR_{n}"] = ("run_0683_rzp_r", "run_0684_rzp_z", n); TWINS[f"grefc_PR_{n}"] = ("run_0683_rzp_r", "run_0685_rzp_p", n)
KEYS = ("P", "ro", "T", "Ux", "Uy", "sonic", "volume", "res_ro")

def load(run, n):
    info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        v = {k: np.asarray(h["VALUE/" + k][:], float) for k in KEYS if "VALUE/" + k in h}
        xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
    nj = v["P"].size // ni
    g = {k: a.reshape(ni, nj)[:, ::-1] for k, a in v.items()}   # [i, jw] (jw = 0 が壁)
    g["umag"] = np.hypot(g["Ux"], g["Uy"]); g["q"] = 0.5 * g["ro"] * g["umag"] ** 2
    X = (xy[:, 0] / S).reshape(ni, nj)[:, ::-1]; Y = (xy[:, 1] / S).reshape(ni, nj)[:, ::-1]
    return g, X, Y, ni, nj

def line_stats(phi):
    d1 = np.diff(phi); d2 = phi[2:] - 2 * phi[1:-1] + phi[:-2]
    d4 = phi[4:] - 4 * phi[3:-1] + 6 * phi[2:-2] - 4 * phi[1:-3] + phi[:-4]
    a = np.sqrt(np.mean(d4 ** 2)) / 16.0; b = np.sqrt(np.mean(d1 ** 2))
    c2 = d2[1:-1]   # 4 階差分の中心と同じ範囲の δ²
    f = float(np.mean(np.sign(c2[1:]) * np.sign(c2[:-1]) < 0)) if c2.size > 1 else float("nan")
    return a, b, (a / b if b > 0 else float("inf")), f

def lines(g, X, Y, ni):
    xc = X[:, 0]; cols = np.where((xc >= -6.0) & (xc <= -3.25))[0]
    out = {"wallnormal": [], "streamwise": []}
    for i in cols:   # 壁法線: jw 58..102 (中心 60..100)
        out["wallnormal"].append((i, slice(58, 103)))
    lo, hi = int(cols[0]) - 2, int(cols[-1]) + 3
    for jw in range(60, 101):
        out["streamwise"].append((slice(lo, hi), jw))
    return out, cols

def field_stats(name, run, n):
    g, X, Y, ni, nj = load(run, n)
    nonfinite = int(sum(np.count_nonzero(~np.isfinite(g[k])) for k in ("P", "ro", "T", "Ux", "Uy")))
    L, cols = lines(g, X, Y, ni); rec = {"run": run, "step": n, "ni": ni, "nj": nj, "n_cols": int(cols.size), "nonfinite": nonfinite}
    for d, ls in L.items():
        acc = {}
        for sel in ls:
            xs, ys = X[sel], Y[sel]; s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xs), np.diff(ys)))])
            acc.setdefault("eta_s", []).append(line_stats(s)[2])
            for k in ("P", "ro", "T", "umag", "res_ro"):
                if k not in g: continue
                a, b, eta, f = line_stats(g[k][sel]); acc.setdefault(f"eta_{k}", []).append(eta); acc.setdefault(f"f_{k}", []).append(f)
                if k == "P": acc.setdefault("aP_over_q", []).append(a / float(np.median(g["q"][sel])))
        rec[d] = {k: float(np.median(v)) for k, v in acc.items()}
        rec[d]["P_eta_p90"] = float(np.percentile(acc["eta_P"], 90)); rec[d]["P_f_p90"] = float(np.percentile(acc["f_P"], 90))
    return rec

def twin_stats(name, rb, rp, n):
    gb, X, Y, ni, nj = load(rb, n); gp, Xp, Yp, _, _ = load(rp, n)
    if not (np.array_equal(X, Xp) and np.array_equal(Y, Yp)): raise SystemExit(f"{name}: 座標が違う")
    xc = X[:, 0]; isl = np.zeros_like(X, bool); cols = (xc >= -6.0) & (xc <= -3.25); isl[np.ix_(cols, np.arange(60, 101))] = True
    V = gb["volume"][isl]; c2 = gb["sonic"][isl] ** 2
    dp = (gp["P"] - gb["P"])[isl]; dro = (gp["ro"] - gb["ro"])[isl]
    du = np.hypot((gp["Ux"] - gb["Ux"])[isl], (gp["Uy"] - gb["Uy"])[isl])
    rms = lambda z: float(np.sqrt(np.sum(V * z * z) / np.sum(V)))
    pp, ss = rms(dp / c2), rms(dro - dp / c2)
    rec = {"base": rb, "pert": rp, "step": n, "rms_dp_over_c2": pp, "rms_entropy_part": ss, "A": pp / ss if ss > 0 else None,
           "dp_over_rhocdu": rms(dp) / rms(gb["ro"][isl] * gb["sonic"][isl] * du) if rms(du) > 0 else None,
           "rms_dT": rms((gp["T"] - gb["T"])[isl]), "rms_dp": rms(dp), "M_island_median": float(np.median(gb["umag"][isl] / gb["sonic"][isl]))}
    L, _ = lines(gb, X, Y, ni)
    for d, ls in L.items():   # 差の格子間隔の成分 (記述)
        for k in ("P", "T"):
            e = [line_stats(gp[k][sel] - gb[k][sel])[2:] for sel in ls]
            rec[f"{d}_eta_d{k}"] = float(np.median([x[0] for x in e])); rec[f"{d}_f_d{k}"] = float(np.median([x[1] for x in e]))
    return rec

F = {k: field_stats(k, *v) for k, v in FIELDS.items()}
T = {k: twin_stats(k, *v) for k, v in TWINS.items()}
qual = {"nonfinite": all(r["nonfinite"] == 0 for r in F.values()),
        "grid_lines_smooth": all(r[d]["eta_s"] <= 0.01 for r in F.values() for d in ("wallnormal", "streamwise"))}
b0 = [F["gref_b0_0608"], F["gref_b0_0619"]]
def present(r):
    return any(r[d]["f_P"] >= 0.7 and r[d]["eta_P"] >= 0.1 and r[d]["eta_umag"] <= r[d]["eta_P"] / 3 for d in ("wallnormal", "streamwise"))
if not all(qual.values()):
    saw = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif all(present(r) for r in b0):
    saw = "S: 停滞帯に圧力の鋸歯がある (速度より格子間隔の上下が強い)"
elif all(r[d]["eta_P"] <= 0.01 for r in b0 for d in ("wallnormal", "streamwise")):
    saw = "N: 停滞帯の圧力に格子間隔の鋸歯は無い"
else:
    saw = "判定不能 (中間)"
def comp(A):
    return None if A is None else ("エントロピーの成分が主" if A <= 0.1 else "圧力の成分が主" if A >= 1 else "混在")
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.78", "qualification": qual, "fields": F, "twins": T,
       "result_sawtooth": saw, "result_twin_components": {k: comp(r["A"]) for k, r in T.items()}}
(OUT / "lowmach_saw.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for k, r in F.items():
    for d in ("wallnormal", "streamwise"):
        q = r[d]
        print(f"{k:16s} {d:10s} η_s {q['eta_s']:.1e}  P η {q['eta_P']:.3g} f {q['f_P']:.2f} a/q {q['aP_over_q']:.2e}  ρ η {q['eta_ro']:.3g} f {q['f_ro']:.2f}  "
              f"T η {q['eta_T']:.3g} f {q['f_T']:.2f}  |u| η {q['eta_umag']:.3g} f {q['f_umag']:.2f}  res_ro η {q.get('eta_res_ro', float('nan')):.3g} f {q.get('f_res_ro', float('nan')):.2f}")
for k, r in T.items():
    print(f"{k:16s} A {r['A']:.3g}  δp/(ρc|δu|) {r['dp_over_rhocdu']:.3g}  RMS δT {r['rms_dT']:.2e} K  RMS δp {r['rms_dp']:.2e} Pa  M {r['M_island_median']:.3f}  "
          f"差の η: 法線 P {r['wallnormal_eta_dP']:.3g} T {r['wallnormal_eta_dT']:.3g} / 流れ P {r['streamwise_eta_dP']:.3g} T {r['streamwise_eta_dT']:.3g}")
print("資格:", qual); print("圧力の鋸歯:", saw); print("双子の差の成分:", out["result_twin_components"])
