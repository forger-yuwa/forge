"""plan tooling-nozzle-core-grid §4.34 (事前登録、codex 諮問 2026-10-11 の 2 回目): B0 の停滞の反復の側の診断 D1。
1 本の run の res_2000〜res_2032 (毎 step の書き出し) から、毎 step の更新 ΔQ_k = Q_k − Q_{k−1} (k = 2001〜2032) を節点ごとに作り、
  1. 交互の更新: 変数ごとに窓 W1 = 2001〜2016・W2 = 2017〜2032 で、節点ごとの ΔQ の lag-1・lag-2 の相関 (Pearson) と RMS の振幅。
     「持続する交互更新」= 両窓で lag-1 ≤ −0.8・lag-2 ≥ 0.8・振幅比 W2/W1 が 0.8〜1.25 の節点。丸めに埋もれる節点
     (窓の RMS(ΔQ) < 1e-12·max|Q|) は外す。その節点の組が各変数の残差二乗和 (窓の平均) の何割を担うか。場所 (壁からの節点番号・x) の分布。
  2. 提案の増分と実際の増分: 平均流は dq_block_old_* (res_k の値、implicitRelax を掛けた最終の補正) と ΔQ_k。SST は
     implicitRelax·R/(V/dt_local + V·src_jac + transport_diag) (res_k の値) と ΔQ_k。内部の節点と境界の節点 (軸 j=0・壁 j=nj−1・入口・出口の列) を分けて集計。
     ずれがあっても、それだけで床の作動とは言わない (FORGE_DIAG_COMMIT_LOSS の commit_loss.csv の n_clip と並べる)。
usage (AWS): python3 cg_b0_diag_an.py <run> → _band_ab/core_grid/d1_<run>.json と d1_<run>.npz (節点ごとの統計、比較用)"""
import csv, json, sys
from pathlib import Path
import h5py
import numpy as np

HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
run = sys.argv[1]; R0 = HERE / run
K = list(range(2000, 2033)); RELAX = 0.7
QN = ["ro", "roUx", "roUy", "roe", "roK", "roOmega"]
RN = {q: "res_" + q for q in QN}
DQN = {"ro": "dq_block_old_0", "roUx": "dq_block_old_1", "roUy": "dq_block_old_2", "roe": "dq_block_old_4"}
ni = int(json.loads((R0 / "prepare_info.json").read_text())["mesh"]["ni"]); S = float(json.loads((R0 / "prepare_info.json").read_text())["scale_m"])

def load(k, names):
    with h5py.File(R0 / f"res_{k}.h5", "r") as h:
        return {n: np.asarray(h["VALUE/" + n][:], float) for n in names}

with h5py.File(R0 / f"res_{K[0]}.h5", "r") as h:
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
n = xy.shape[0]; nj = n // ni
jj = np.tile(np.arange(nj), ni)                                   # 0 = 軸、nj−1 = 壁 (列ごとに軸から壁へ)
xx = xy[:, 0] / S
bnd = (jj == 0) | (jj == nj - 1) | (np.repeat(np.arange(ni), nj) == 0) | (np.repeat(np.arange(ni), nj) == ni - 1)

prev = load(K[0], QN)
DQ = {q: np.empty((len(K) - 1, n)) for q in QN}
Rsq = {q: [np.zeros(n), np.zeros(n)] for q in QN}
Qmax = {q: float(np.abs(prev[q]).max()) for q in QN}
prop = {q: {"interior": [0.0, 0.0, 0], "boundary": [0.0, 0.0, 0]} for q in QN}   # Σ|req−act|, Σ|req|, 節点数 (|req|>0)
for t, k in enumerate(K[1:]):
    cur = load(k, QN + [RN[q] for q in QN] + list(DQN.values()) + ["dt_local", "src_jac_k", "src_jac_omega", "transport_diag_k", "transport_diag_omega", "volume"])
    w = 0 if k <= 2016 else 1
    for q in QN:
        d = cur[q] - prev[q]; DQ[q][t] = d
        Rsq[q][w] += cur[RN[q]] ** 2
        if q in DQN:
            req = cur[DQN[q]]
        elif q == "roK":
            V = cur["volume"]; req = RELAX * cur["res_roK"] / np.maximum(V / np.maximum(cur["dt_local"], 1e-30) + V * cur["src_jac_k"] + cur["transport_diag_k"], 1e-30)
        else:
            V = cur["volume"]; req = RELAX * cur["res_roOmega"] / np.maximum(V / np.maximum(cur["dt_local"], 1e-30) + V * cur["src_jac_omega"] + cur["transport_diag_omega"], 1e-30)
        for reg, m in (("interior", ~bnd), ("boundary", bnd)):
            mm = m & (req != 0)
            prop[q][reg][0] += float(np.abs(req[mm] - d[mm]).sum()); prop[q][reg][1] += float(np.abs(req[mm]).sum()); prop[q][reg][2] += int(mm.sum())
    prev = {q: cur[q] for q in QN}

def lagcorr(X, lag):
    a, b = X[:-lag], X[lag:]
    a = a - a.mean(0); b = b - b.mean(0)
    den = np.sqrt((a * a).sum(0) * (b * b).sum(0))
    return np.where(den > 0, (a * b).sum(0) / np.where(den > 0, den, 1), np.nan)

out = {"run": run, "ni": ni, "nj": nj, "n": n, "steps": [K[0], K[-1]], "vars": {}}
npz = {}
for q in QN:
    X = DQ[q]; W = [X[:16], X[16:]]
    l1 = [lagcorr(Wi, 1) for Wi in W]; l2 = [lagcorr(Wi, 2) for Wi in W]; amp = [np.sqrt((Wi ** 2).mean(0)) for Wi in W]
    live = (amp[0] > 1e-12 * Qmax[q]) & (amp[1] > 1e-12 * Qmax[q])
    ratio = np.where(amp[0] > 0, amp[1] / np.where(amp[0] > 0, amp[0], 1), np.nan)
    alt = live & (l1[0] <= -0.8) & (l1[1] <= -0.8) & (l2[0] >= 0.8) & (l2[1] >= 0.8) & (ratio >= 0.8) & (ratio <= 1.25)
    rs = Rsq[q][0] + Rsq[q][1]; tot = float(rs.sum())
    top = np.argsort(rs)[::-1][: max(1, n // 1000)]                # 残差二乗和の上位 0.1 % の節点
    out["vars"][q] = {
        "live_nodes": int(live.sum()), "alt_nodes": int(alt.sum()),
        "alt_residual_share": float(rs[alt].sum() / tot) if tot > 0 else None,
        "alt_amp_rms": float(np.sqrt((amp[0][alt] ** 2 + amp[1][alt] ** 2).mean() / 2)) if alt.any() else 0.0,
        "lag1_median_live": [float(np.nanmedian(l1[0][live])), float(np.nanmedian(l1[1][live]))] if live.any() else None,
        "lag1_median_top_residual": [float(np.nanmedian(l1[0][top])), float(np.nanmedian(l1[1][top]))],
        "top_residual_share": float(rs[top].sum() / tot) if tot > 0 else None,
        "alt_where": {"j_from_wall_hist": np.histogram(nj - 1 - jj[alt], bins=[0, 1, 2, 5, 10, 20, 40, 80, nj])[0].tolist(),
                      "x_hist": np.histogram(xx[alt], bins=[-13, -2, 2, 12, 40, 96])[0].tolist()} if alt.any() else None,
        "top_residual_where": {"j_from_wall_hist": np.histogram(nj - 1 - jj[top], bins=[0, 1, 2, 5, 10, 20, 40, 80, nj])[0].tolist(),
                               "x_hist": np.histogram(xx[top], bins=[-13, -2, 2, 12, 40, 96])[0].tolist()},
        "proposed_vs_actual": {reg: {"lost_frac": (v[0] / v[1]) if v[1] > 0 else None, "nodes": v[2]} for reg, v in prop[q].items()},
    }
    npz[f"{q}_alt"] = alt; npz[f"{q}_amp1"] = amp[0].astype(np.float32); npz[f"{q}_amp2"] = amp[1].astype(np.float32)
    npz[f"{q}_l1"] = np.nanmin(np.vstack(l1), 0).astype(np.float32); npz[f"{q}_rsq"] = rs.astype(np.float32)
out["bins"] = {"j_from_wall": [0, 1, 2, 5, 10, 20, 40, 80, nj], "x_over_rt": [-13, -2, 2, 12, 40, 96]}
cl = R0 / "commit_loss.csv"
if cl.is_file():
    rows = [r for r in csv.DictReader(open(cl)) if 2001 <= int(float(r["step"])) <= 2032]
    agg = {}
    for r in rows:
        a = agg.setdefault((r["qty"], r["region"]), {"n_req": 0, "n_lost": 0, "n_clip": 0, "S_req": 0.0, "S_act": 0.0})
        for c in ("n_req", "n_lost", "n_clip"): a[c] += int(float(r[c]))
        a["S_req"] += float(r["S_req"]); a["S_act"] += float(r["S_act"])
    out["commit_loss"] = {f"{q}|{reg}": v for (q, reg), v in agg.items()}
(OUTD / f"d1_{run}.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
np.savez_compressed(OUTD / f"d1_{run}.npz", **npz)
for q in QN:
    v = out["vars"][q]
    print(f"{q:7s} 生きた節点 {v['live_nodes']} 交互 {v['alt_nodes']} 残差の割合 {v['alt_residual_share']} 上位0.1%の lag1 {v['lag1_median_top_residual']} "
          f"提案と実際のずれ 内部 {v['proposed_vs_actual']['interior']['lost_frac']} 境界 {v['proposed_vs_actual']['boundary']['lost_frac']}")
print(f"[d1] {OUTD / f'd1_{run}.json'}")
