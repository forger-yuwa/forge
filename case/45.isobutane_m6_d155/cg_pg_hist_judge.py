"""plan tooling-nozzle-core-grid §4.92 P-h の判定 (事前登録、codex 諮問 2026-10-13 pg-interp): 評価の履歴の差。
R_A = 継続中の評価 (800 sweep の run の res_2)、R_B = 同じ状態を新しい run に移した再起動直後の評価 (run_0704・0705 の res_1)。
射影 (圧力残差 R_p/c²、r と c² は Q0 = run_0693 の res_20000)・体積重み・領域は §4.90 の判定 (cg_pg_judge.py) と同じ。
D = R_B − R_A、N_old = ΔR(½) − ½ΔR_A(1)、N_new = N_old − D/2 = ΔR(½) − ½ΔR_B(1) (ΔR(s) = R(Q0 + sδ★) − R(Q0)、R(Q0) は各腕の 5 sweep の run の res_1)。
資格: 2 本の restart が res_1 とビット一致・スイッチの行・非有限なし・RUN_RC 0。
履歴の差で説明できる: 両腕で ‖N_new‖/‖½ΔR_B(1)‖ ≤ 0.05 (縮流部)。履歴の差では説明できない: 両腕で ‖D/2‖ ≤ 0.1 ‖N_old‖ (縮流部)。それ以外は判定不能。
usage: python3 cg_pg_hist_judge.py → _band_ab/core_grid/pg_hist_judge.json"""
import json
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
BASE = "run_0693_pdt_off"; SEED = HERE / BASE / "res_20000.h5"
ARMS = {"GI": {"b0": "run_0696_pg_gi5", "full": "run_0698_pg_gi800", "half": "run_0702_pg_gih", "B": "run_0704_pg_gib"},
        "GB": {"b0": "run_0699_pg_gb5", "full": "run_0701_pg_gb800", "half": "run_0703_pg_gbh", "B": "run_0705_pg_gbb"}}
RES = ("res_ro", "res_roUx", "res_roUy", "res_roe")
info = json.loads((HERE / BASE / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
def rd(path, keys):
    with h5py.File(path, "r") as h: return {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}
with h5py.File(SEED, "r") as h: x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
Q0 = rd(SEED, ("ro", "P", "sonic", "h0", "Ux", "Uy", "volume"))
nj = x.size // ni; jw = nj - 1 - np.arange(x.size) % nj; V = Q0["volume"]
REG = {"contraction": (x >= -8) & (x <= -2), "island": (x >= -6) & (x <= -3.25) & (jw >= 60) & (jw <= 100),
       "wall_interior": (x >= -2) & (x <= 0) & (jw >= 1) & (jw <= 10), "wall_node": (x >= -2) & (x <= 0) & (jw == 0)}
c2 = Q0["sonic"] ** 2; gam = c2 * Q0["ro"] / Q0["P"]; kap = gam - 1; ek = 0.5 * (Q0["Ux"] ** 2 + Q0["Uy"] ** 2); chi = c2 - kap * (Q0["h0"] - ek)
def Rp(R): return ((chi + kap * ek) * R["res_ro"] - kap * Q0["Ux"] * R["res_roUx"] - kap * Q0["Uy"] * R["res_roUy"] + kap * R["res_roe"]) / c2 / V
def nrm(z, m): return float(np.sqrt(np.sum(V[m] * z[m] ** 2) / np.sum(V[m])))
logt = lambda run: "".join((HERE / run / f).read_text(errors="replace") for f in ("forge_run.log", "run_case_stdout.log") if (HERE / run / f).exists())
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.92", "arms": {}}
qual = {"restart_bitwise_from_res1": True, "switch_line": True, "finite_completed": True}
for arm, a in ARMS.items():
    rl = (HERE / a["B"] / "restart.log").read_text()
    qual["restart_bitwise_from_res1"] &= "VERDICT: OK (9 量を移した、SRC とビット一致)" in rl and f"{a['full']}/res_1.h5" in (HERE / "cg_pg_hist.log").read_text()
    qual["switch_line"] &= logt(a["B"]).count("FORGE_DIAG_PRECOND_NO_DT_SCALE=1:") >= 1
    b0 = Rp(rd(HERE / a["b0"] / "res_1.h5", RES)); RA = Rp(rd(HERE / a["full"] / "res_2.h5", RES)); RH = Rp(rd(HERE / a["half"] / "res_1.h5", RES))
    RBraw = rd(HERE / a["B"] / "res_1.h5", RES); RB = Rp(RBraw)
    qual["finite_completed"] &= all(np.all(np.isfinite(v)) for v in RBraw.values()) and (HERE / a["B"] / "RUN_RC").read_text().strip() == "0"
    D = RB - RA; Nold = (RH - b0) - 0.5 * (RA - b0); Nnew = (RH - b0) - 0.5 * (RB - b0)
    out["arms"][arm] = {rg: {"norm_D_half": nrm(0.5 * D, m), "norm_N_old": nrm(Nold, m), "norm_N_new": nrm(Nnew, m),
                             "lin_old": nrm(Nold, m) / max(nrm(0.5 * (RA - b0), m), 1e-300), "lin_new": nrm(Nnew, m) / max(nrm(0.5 * (RB - b0), m), 1e-300),
                             "norm_dR_A": nrm(RA - b0, m), "norm_dR_B": nrm(RB - b0, m), "norm_dR_half": nrm(RH - b0, m)} for rg, m in REG.items()}
c = {arm: out["arms"][arm]["contraction"] for arm in ARMS}
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif all(c[a]["lin_new"] <= 0.05 for a in ARMS):
    res = "履歴の差で説明できる: 評価の履歴をそろえると両腕で線形性の誤差 ≤ 0.05 (第 1 仮説を支持)"
elif all(c[a]["norm_D_half"] <= 0.1 * c[a]["norm_N_old"] for a in ARMS):
    res = "履歴の差では説明できない: 両腕で ‖D/2‖ ≤ 0.1 ‖N_old‖ (第 1 仮説を主因として退け、第 2 仮説を優先)"
else:
    res = "判定不能 (中間・腕で割れる)"
out["qualification"] = qual; out["result"] = res
(OUTD / "pg_hist_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for arm in ARMS:
    for rg, v in out["arms"][arm].items():
        print(f"{arm} {rg:13s} ‖D/2‖ {v['norm_D_half']:.3e} ‖N_old‖ {v['norm_N_old']:.3e} ‖N_new‖ {v['norm_N_new']:.3e}  線形性 旧 {v['lin_old']:.3f} 新 {v['lin_new']:.3f}  ‖ΔR_A‖ {v['norm_dR_A']:.3e} ‖ΔR_B‖ {v['norm_dR_B']:.3e} ‖ΔR(½)‖ {v['norm_dR_half']:.3e}")
print("資格:", qual); print("結果:", res)
