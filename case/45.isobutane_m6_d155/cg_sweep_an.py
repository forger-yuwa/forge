"""plan tooling-nozzle-core-grid §4.42 M-d: FORGE_LINE_DUMP の書き出し (1 外反復、nStepInner 41) から、内部反復の解き残しを測る。
書き出し: <run>/linedump/meta.txt (名前・行数・列数) と <名前>.f64 (float64、(節点, 列) の行優先)。節点はラインの順。
  node_line (節点番号・ライン番号)、D (5×5、時間項 V/Δτ 込みの対角)・Kprev・Knext (5×5)、dt_vol (dt_local・体積)、
  flags_wall_iso_axis、各 sweep k の rhs_s{k}・dqnew_s{k}・dqold_s{k}。
実装の反復 (codex 2026-10-11 sweep の定義): ライン行列 L (対角 D、下 −Kprev、上 −Knext)、ω = implicitRelax、b_k = rhs_s{k} として
  L δ_{k+1} = ω b_k、δ_{k+1} = dqnew_s{k} (δ_0 = 0)。δ₅ = dqnew_s4、δ₃₀ = dqnew_s29、δ₄₀ = dqnew_s39、b₄₀ = rhs_s40。
指標 (拘束の節点 = 壁・等温壁・軸のフラグがある節点は除き、別に数える。ノルムは体積重み付きの RMS、成分は ρ_ref・ρ_ref a_ref・ρ_ref a_ref² で無次元化):
  後退誤差 be_k = max_行 |L δ_{k+1} − ω b_k| / (Σ_j |L_ij||δ_j| + |ω b_k|)  (k = 0..40、分母 0 の行は数えて除く)
  線形残差 g₄₀ = (L/ω) δ₄₀ − b₄₀、比 ‖g₄₀/V‖ / ‖b₀/V‖
  conv = ‖δ₄₀ − δ₃₀‖ / ‖δ₄₀‖
  u5 = ‖δ₅ − δ₄₀‖ / ‖δ₄₀‖ (本番の 5 回の解き残し)、u の推移 ‖δ_k − δ₄₀‖ / ‖δ₄₀‖
  κ_k = ‖δ_{k+1} − δ_k‖ / ‖δ_k − δ_{k−1}‖ (k = 10..39) の中央値・最大値
  時間項の割合 (V/Δτ) / D₀₀ と dt_local の中央値
領域: 島 = x ∈ [−6, −3.25] r_t かつ壁から 60〜100 節点目、lines = 書き出したライン全体 (x ∈ [−8, −2] の列)。
usage: python3 cg_sweep_an.py <run> <tag> → _band_ab/core_grid/sweep_<tag>.json (と rhs_s0 の sha256)"""
import hashlib, json, re, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run, tag = sys.argv[1], sys.argv[2]
DD = HERE / run / "linedump"
meta = {}
for ln in (DD / "meta.txt").read_text().splitlines():
    if ln.startswith("#") or not ln.strip():
        continue
    n, r, c = ln.split(); meta[n] = (int(r), int(c))
def get(name):
    r, c = meta[name]
    return np.fromfile(DD / f"{name}.f64", dtype=np.float64).reshape(r, c)
K = max(int(n[len("dqnew_s"):]) for n in meta if n.startswith("dqnew_s")) + 1      # sweep の数 (41 を想定)
if K < 41:
    raise SystemExit(f"sweep が {K} 回しかない (41 回が要る)")
cfg = (HERE / run / "solverConfig.yaml").read_text()
key = lambda k: float(re.search(k + r":\s*([0-9.eE+-]+)", cfg).group(1))
roR, aR, om = key("limiterRoRef"), key("limiterARef"), key("implicitRelax")
info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
nl = get("node_line"); nodes = nl[:, 0].astype(np.int64); line = nl[:, 1].astype(np.int64)
with h5py.File(HERE / run / "res_0.h5", "r") as h:
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
nj = len(xy) // ni
x = xy[nodes, 0] / S; jw = nj - 1 - (nodes % nj)
N = len(nodes)
prv = np.arange(N) - 1; nxt = np.arange(N) + 1
prv[np.r_[True, line[1:] != line[:-1]]] = -1                      # ラインの先頭は下側なし
nxt[np.r_[line[1:] != line[:-1], True]] = -1                      # ラインの末尾は上側なし
Dm = get("D").reshape(N, 5, 5); Kp = get("Kprev").reshape(N, 5, 5); Kn = get("Knext").reshape(N, 5, 5)
dv = get("dt_vol"); dt, V = dv[:, 0], dv[:, 1]
fl = get("flags_wall_iso_axis"); con = np.any(fl != 0, axis=1)
def Lmul(d, absval=False):
    """(L δ) と、absval なら (|L| |δ|)"""
    A, P, Q = (np.abs(Dm), np.abs(Kp), np.abs(Kn)) if absval else (Dm, -Kp, -Kn)
    dd = np.abs(d) if absval else d
    out = np.einsum("nij,nj->ni", A, dd)
    m = prv >= 0; out[m] += np.einsum("nij,nj->ni", P[m], dd[prv[m]])
    m = nxt >= 0; out[m] += np.einsum("nij,nj->ni", Q[m], dd[nxt[m]])
    return out
dq = [np.zeros((N, 5))] + [get(f"dqnew_s{k}") for k in range(K)]     # dq[k] = δ_k (k = 0..K)
rhs = [get(f"rhs_s{k}") for k in range(K)]                           # rhs[k] = b_k (k = 0..K−1)
be, bez = [], []
for k in range(K):
    num = np.abs(Lmul(dq[k + 1]) - om * rhs[k]); den = Lmul(dq[k + 1], True) + np.abs(om * rhs[k])
    ok = den > 0; bez.append(int((~ok).sum())); be.append(float(np.max(num[ok] / den[ok])))
k40 = 40
g40 = Lmul(dq[k40]) / om - rhs[k40]
sc = np.array([roR, roR * aR, roR * aR, roR * aR, roR * aR * aR])
reg = {"island": (x >= -6.0) & (x <= -3.25) & (jw >= 60) & (jw <= 100) & ~con, "lines": ~con}
def vn(a, m, ci):
    """体積重み付き RMS (a は (N,5)、成分 ci)"""
    w = V[m]; b = a[m][:, ci] if isinstance(ci, int) else a[m]
    b2 = b * b if isinstance(ci, int) else np.sum(b * b, axis=1)
    return float(np.sqrt(np.sum(w * b2) / np.sum(w)))
out = {"run": run, "tag": tag, "K": K, "omega": om, "nodes": N, "lines": int(len(np.unique(line))), "constrained_nodes": int(con.sum()),
       "rhs_s0_sha256": hashlib.sha256(rhs[0].tobytes()).hexdigest(), "state_sha256": hashlib.sha256(get("state_ro_roU_roe_cp_gamma").tobytes()).hexdigest(),
       "nodes_sha256": hashlib.sha256(nodes.tobytes()).hexdigest(),
       "nonfinite": int(sum(np.count_nonzero(~np.isfinite(a)) for a in dq + rhs)),
       "backward_error_max": max(be), "backward_error": be, "backward_error_zero_den_rows": bez, "region": {}}
tt = V / dt / Dm[:, 0, 0]
d = [a / sc for a in dq]; d40 = d[k40]
bV0 = rhs[0] / V[:, None] / sc; gV = g40 / V[:, None] / sc
for rn, m in reg.items():
    o = {"n": int(m.sum()), "tt_median": float(np.median(tt[m])), "dt_median": float(np.median(dt[m]))}
    for cn, ci in (("ro", 0), ("roe", 4), ("all5", None)):
        n40 = vn(d40, m, ci)
        kap = [vn(d[k + 1] - d[k], m, ci) / vn(d[k] - d[k - 1], m, ci) for k in range(10, k40)]
        o[cn] = {"u5": vn(d[5] - d40, m, ci) / n40, "u1": vn(d[1] - d40, m, ci) / n40,
                 "conv_40_30": vn(d40 - d[30], m, ci) / n40,
                 "g40_over_b0": vn(gV, m, ci) / vn(bV0, m, ci),
                 "kappa_median": float(np.median(kap)), "kappa_max": float(np.max(kap)),
                 "u_hist": [vn(d[k] - d40, m, ci) / n40 for k in range(1, k40 + 1)]}
    out["region"][rn] = o
(HERE / "_band_ab" / "core_grid" / f"sweep_{tag}.json").write_text(json.dumps(out, indent=1))
print(f"{tag}: 節点 {N} (ライン {out['lines']}、拘束 {out['constrained_nodes']})、非有限 {out['nonfinite']}、後退誤差の最大 {out['backward_error_max']:.2e}")
for rn, o in out["region"].items():
    print(f"  {rn:7s} n {o['n']:6d} tt {o['tt_median']:.3f} " + "  ".join(
        f"{c}: u5 {o[c]['u5']:.3f} g40/b0 {o[c]['g40_over_b0']:.1e} conv {o[c]['conv_40_30']:.1e} κ {o[c]['kappa_median']:.3f}/{o[c]['kappa_max']:.3f}"
        for c in ("ro", "roe", "all5")))
