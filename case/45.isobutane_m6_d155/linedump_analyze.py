"""実ライン行列の host 再解 A/B (plan time_integration-line-viscous-jacobian §6.2、codex 諮問 2026-10-09)。
forge の FORGE_LINE_DUMP_DIR が書いた D・Kprev・Knext・各 sweep の rhs・dq から、ラインごとに block 三重対角
  D_k ΔQ_k − Kprev_k ΔQ_{k−1} − Knext_k ΔQ_{k+1} = rhs_k
を密行列で組み、(A) σ = 0 (そのまま) と (B) σ = 1 (ライン面のスカラー 2ν_eff δ/dcc の和を拘束でない行の対角へ戻す) で解く。
  - (A) の解 × implicitRelax と CUDA の dq_new の相対差 (実装経路の確認)
  - 補正の物理量: δu = (δ(ρu) − u δρ)/ρ、δT = {δ(ρE) − u·δ(ρu) + (½|u|² − e) δρ}/(ρ c_v)
  - 無次元化 (S = diag(ρ_i, ρ_i c, ρ_i c, ρ_i c, ρ_i c²)、c = 300 m/s で S⁻¹ A S) した行列の最小特異値と右特異ベクトル (弱いモード) の密度の割合、
    rhs のそのモードへの射影、補正の二乗ノルムに占める密度の成分の割合
usage: python3 linedump_analyze.py <dump_dir> [--sweeps 5] → <dump_dir>/analysis.json と標準出力
"""
import argparse
import json
from pathlib import Path

import numpy as np


def load(d: Path):
    meta = {}
    for line in (d / "meta.txt").read_text().splitlines():
        if line.startswith("#"):
            meta["_head"] = line
            continue
        name, r, c = line.split()
        meta[name] = (int(r), int(c))
    arr = {k: np.fromfile(d / f"{k}.f64", dtype=np.float64).reshape(v) for k, v in meta.items() if k != "_head"}
    relax = float(meta["_head"].split("implicitRelax")[1].split()[0])
    return arr, relax


def main(d: Path, nsweep: int):
    A, relax = load(d)
    nl = A["node_line"]; nodes = nl[:, 0].astype(int); lines = nl[:, 1].astype(int)
    D = A["D"].reshape(-1, 5, 5); Kp = A["Kprev"].reshape(-1, 5, 5); Kn = A["Knext"].reshape(-1, 5, 5)
    sv = A["scalar_visc_line"][:, 0]
    fl = A["flags_wall_iso_axis"]; st = A["state_ro_roU_roe_cp_gamma"]
    dtv = A["dt_vol"]
    out = {"dir": str(d), "implicitRelax": relax, "lines": []}
    for L in np.unique(lines):
        idx = np.where(lines == L)[0]; n = len(idx)
        rho = st[idx, 0]; u = st[idx, 1:4] / rho[:, None]; E = st[idx, 4] / rho; cp = st[idx, 5]; gam = st[idx, 6]
        cv = cp / gam; q2 = (u ** 2).sum(1); e = E - 0.5 * q2
        c = np.full(n, 300.0)   # 無次元化の速度の尺度 (TP は e に基準のずれがあり e から音速を作れないので一定値。尺度にだけ使う)
        dec = np.zeros((n, 5), bool)
        dec[:, 1:4] |= (fl[idx, 0] == 1)[:, None]; dec[:, 4] |= (fl[idx, 1] == 1); dec[:, 2] |= (fl[idx, 2] == 1)   # 拘束の行 (壁の運動量・等温壁のエネルギー・軸の半径運動量)
        rec = {"line": int(L), "n": n, "node_first": int(nodes[idx[0]]), "node_last": int(nodes[idx[-1]]),
               "VdtOverScalar_median": float(np.median((dtv[idx, 1] / np.maximum(dtv[idx, 0], 1e-300)) / np.maximum(sv[idx], 1e-300)))}
        for sigma in (0, 1):
            M = np.zeros((5 * n, 5 * n))
            for k in range(n):
                Dk = D[idx[k]].copy()
                if sigma:
                    for rr in range(5):
                        if not dec[k, rr]:
                            Dk[rr, rr] += sv[idx[k]]
                M[5 * k:5 * k + 5, 5 * k:5 * k + 5] = Dk
                if k > 0:
                    M[5 * k:5 * k + 5, 5 * (k - 1):5 * k] = -Kp[idx[k]]
                if k < n - 1:
                    M[5 * k:5 * k + 5, 5 * (k + 1):5 * k + 10] = -Kn[idx[k]]
            # 無次元化 S⁻¹ M S、S = diag(ρ, ρc, ρc, ρc, ρc²)
            sc = np.stack([rho, rho * c, rho * c, rho * c, rho * c * c], 1).reshape(-1)
            Mn = (M * sc[None, :]) / sc[:, None]
            U_, s_, Vt = np.linalg.svd(Mn)
            vmin = Vt[-1].reshape(n, 5)
            dens_share_mode = float((vmin[:, 0] ** 2).sum() / (vmin ** 2).sum())
            tot = max((vmin ** 2).sum(), 1e-300)
            res = {"sigma_min": float(s_[-1]), "sigma_max": float(s_[0]), "cond": float(s_[0] / s_[-1]), "weak_mode_density_share": dens_share_mode,
                   "weak_mode_component_share": {k: float((vmin[:, j] ** 2).sum() / tot) for j, k in enumerate(("rho", "rhou", "rhov", "rhow", "rhoE"))},
                   "weak_mode_peak_node_index": int(np.argmax((vmin ** 2).sum(1))),
                   "sigma_smallest5": [float(v) for v in s_[-5:]]}
            sweeps = []
            for s in range(nsweep):
                if f"rhs_s{s}" not in A:
                    break
                b = A[f"rhs_s{s}"][idx].reshape(-1)
                x = np.linalg.solve(M, b)
                xn = (x / sc).reshape(n, 5)
                proj = float(abs(U_[:, -1] @ (b / sc)) / max(np.linalg.norm(b / sc), 1e-300))
                X = x.reshape(n, 5)
                du = (X[:, 1:4] - u * X[:, :1]) / rho[:, None]
                dT = (X[:, 4] - (u * X[:, 1:4]).sum(1) + (0.5 * q2 - e) * X[:, 0]) / (rho * cv)
                tn = max((xn ** 2).sum(), 1e-300)
                sw = {"sweep": s, "rhs_proj_on_weak_mode": proj,
                      "corr_component_share": {k: float((xn[:, j] ** 2).sum() / tn) for j, k in enumerate(("rho", "rhou", "rhov", "rhow", "rhoE"))},
                      "corr_peak_node_index": int(np.argmax((xn ** 2).sum(1))),
                      "corr_density_share": float((xn[:, 0] ** 2).sum() / max((xn ** 2).sum(), 1e-300)),
                      "max_abs_drho_over_rho": float(np.max(np.abs(X[:, 0] / rho))), "max_abs_du": float(np.max(np.abs(du))), "max_abs_dT": float(np.max(np.abs(dT))),
                      "norm_scaled_corr": float(np.linalg.norm(xn))}
                if sigma == 0 and f"dqnew_s{s}" in A:
                    cu = A[f"dqnew_s{s}"][idx].reshape(-1)
                    sw["cuda_vs_host_rel"] = float(np.max(np.abs(cu - relax * x)) / max(np.max(np.abs(relax * x)), 1e-300))
                sweeps.append(sw)
            res["sweeps"] = sweeps
            rec[f"sigma{sigma}"] = res
        a0 = rec["sigma0"]["sweeps"]; a1 = rec["sigma1"]["sweeps"]
        if a0 and a1:
            rec["density_corr_ratio_B_over_A_sweep0"] = a1[0]["max_abs_drho_over_rho"] / max(a0[0]["max_abs_drho_over_rho"], 1e-300)
        out["lines"].append(rec)
        s0 = rec["sigma0"]; s1 = rec["sigma1"]
        print(f"ライン {L} (節点 {rec['node_first']}〜{rec['node_last']}、{n} 点、V/Δτ ÷ スカラー粘性 中央値 {rec['VdtOverScalar_median']:.3g}):")
        for tag, r_ in (("A σ=0", s0), ("B σ=1", s1)):
            sw0 = r_["sweeps"][0] if r_["sweeps"] else {}
            print(f"  {tag}: 最小特異値 {r_['sigma_min']:.3e}・条件数 {r_['cond']:.3e}・弱いモードの密度の割合 {r_['weak_mode_density_share']:.3f}; "
                  f"sweep 0: |δρ/ρ| 最大 {sw0.get('max_abs_drho_over_rho', float('nan')):.3e}・|δu| {sw0.get('max_abs_du', float('nan')):.3e}・|δT| {sw0.get('max_abs_dT', float('nan')):.3e}・"
                  f"補正の密度の割合 {sw0.get('corr_density_share', float('nan')):.3f}・rhs の弱いモードへの射影 {sw0.get('rhs_proj_on_weak_mode', float('nan')):.3e}"
                  + (f"・CUDA との相対差 {sw0['cuda_vs_host_rel']:.2e}" if 'cuda_vs_host_rel' in sw0 else ""))
    (d / "analysis.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--sweeps", type=int, default=5)
    a = ap.parse_args()
    main(Path(a.dir), a.sweeps)
