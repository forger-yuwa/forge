"""plan tooling-nozzle-core-grid §4.94 P-tf の判定 (事前登録、ユーザ決定・codex 諮問 2026-10-13 thermofloat): 熱物性の評価の経路 (thermoFloat 1 → 0)。
経路 m (1 = 旧、0 = 新) ごとに、圧力残差 R_p/c² (r と c² は Q0 = run_0693 の res_20000、§4.92 と同じ射影) で
N_m = R_m(Q½) − ½R_m(Q1) − ½R_m(Q0)、ℓ_m = ‖N_m‖/‖½[R_m(Q1) − R_m(Q0)]‖。R は各 run の res_1 の残差 (更新前、立ち上げ直後の評価)。
旧: Q0 = run_0696 (GI)・run_0699 (GB)、Q½ = run_0702・run_0703、Q1 = run_0704・run_0705。新: Q0 = run_0900 (共用)、Q½ = run_0901・run_0902、Q1 = run_0903・run_0904。
資格: 新 5 本の thermoFloat が 0、Q0・Q1 の restart がビット一致、Q½ の入力が旧の入力と同じ (sha)、旧の Q0 の GI と GB の縮流部の R_p が相対 1e−6 以内、
  スイッチの行・非有限なし・RUN_RC 0。
A: 両腕で ℓ₀ ≤ 0.05 かつ ‖N₀‖ ≤ 0.1 ‖N₁‖。B: 両腕で ‖N₀‖ ≥ 0.9 ‖N₁‖。それ以外・信号不足・資格不足は判定不能。
補助: T・P・音速の比例の誤差 (縮流部と、縮流部の壁から 1〜10 節点目) を旧・新で併記。
usage: python3 cg_tf_judge.py → _band_ab/core_grid/tf_judge.json"""
import hashlib, json, re
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
BASE = "run_0693_pdt_off"; SEED = HERE / BASE / "res_20000.h5"
OLD = {"GI": ("run_0696_pg_gi5", "run_0702_pg_gih", "run_0704_pg_gib"), "GB": ("run_0699_pg_gb5", "run_0703_pg_gbh", "run_0705_pg_gbb")}
NEW = {"GI": ("run_0900_tf_q0", "run_0901_tf_gih", "run_0903_tf_gib"), "GB": ("run_0900_tf_q0", "run_0902_tf_gbh", "run_0904_tf_gbb")}
RES = ("res_ro", "res_roUx", "res_roUy", "res_roe"); PRIM = ("T", "P", "sonic")
info = json.loads((HERE / BASE / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
def rd(path, keys):
    with h5py.File(path, "r") as h: return {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}
with h5py.File(SEED, "r") as h: x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
Q0 = rd(SEED, ("ro", "P", "sonic", "h0", "Ux", "Uy", "volume"))
nj = x.size // ni; jw = nj - 1 - np.arange(x.size) % nj; V = Q0["volume"]
REG = {"contraction": (x >= -8) & (x <= -2), "contraction_wall1_10": (x >= -8) & (x <= -2) & (jw >= 1) & (jw <= 10), "island": (x >= -6) & (x <= -3.25) & (jw >= 60) & (jw <= 100),
       "wall_interior": (x >= -2) & (x <= 0) & (jw >= 1) & (jw <= 10)}
c2 = Q0["sonic"] ** 2; gam = c2 * Q0["ro"] / Q0["P"]; kap = gam - 1; ek = 0.5 * (Q0["Ux"] ** 2 + Q0["Uy"] ** 2); chi = c2 - kap * (Q0["h0"] - ek)
def Rp(R): return ((chi + kap * ek) * R["res_ro"] - kap * Q0["Ux"] * R["res_roUx"] - kap * Q0["Uy"] * R["res_roUy"] + kap * R["res_roe"]) / c2 / V
def nrm(z, m): return float(np.sqrt(np.sum(V[m] * z[m] ** 2) / np.sum(V[m])))
def cfg(run, key):
    m = re.search(r"'%s' in '[^']*': ([^\s]+)" % key, (HERE / run / "forge_run.log").read_text(errors="replace")); return m.group(1) if m else None
logt = lambda run: "".join((HERE / run / f).read_text(errors="replace") for f in ("forge_run.log", "run_case_stdout.log") if (HERE / run / f).exists())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
newruns = sorted({r for t in NEW.values() for r in t})
qual = {"thermoFloat0": all(cfg(r, "thermoFloat") == "0" for r in newruns),
        "restart_bitwise": all("VERDICT: OK (9 量を移した、SRC とビット一致)" in (HERE / r / "restart.log").read_text() for r in ("run_0900_tf_q0", "run_0903_tf_gib", "run_0904_tf_gbb")),
        "half_input_same": sha(HERE / "run_0901_tf_gih" / "nozzle.h5") == sha(HERE / "run_0702_pg_gih" / "nozzle.h5") and sha(HERE / "run_0902_tf_gbh" / "nozzle.h5") == sha(HERE / "run_0703_pg_gbh" / "nozzle.h5"),
        "switch_line": all(logt(r).count("FORGE_DIAG_PRECOND_NO_DT_SCALE=1:") >= 1 for r in newruns),
        "finite_completed": True}
for r in newruns:
    R = rd(HERE / r / "res_1.h5", RES)
    qual["finite_completed"] &= all(np.all(np.isfinite(v)) for v in R.values()) and (HERE / r / "RUN_RC").read_text().strip() == "0"
b_gi, b_gb = Rp(rd(HERE / OLD["GI"][0] / "res_1.h5", RES)), Rp(rd(HERE / OLD["GB"][0] / "res_1.h5", RES))
m = REG["contraction"]; q0share = nrm(b_gi - b_gb, m) / nrm(b_gi, m); qual["old_Q0_gamma_independent"] = q0share <= 1e-6
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.94", "old_Q0_GI_vs_GB_rel": q0share, "arms": {}}
for arm in ("GI", "GB"):
    out["arms"][arm] = {}
    for path, runs in (("old", OLD[arm]), ("new", NEW[arm])):
        R0, Rh, R1 = (Rp(rd(HERE / r / "res_1.h5", RES)) for r in runs)
        P0, Ph, P1 = (rd(HERE / r / "res_1.h5", PRIM) for r in runs)
        N = Rh - 0.5 * R1 - 0.5 * R0
        rec = {rg: {"normN": nrm(N, mm), "ell": nrm(N, mm) / max(nrm(0.5 * (R1 - R0), mm), 1e-300), "norm_halfdR": nrm(0.5 * (R1 - R0), mm)} for rg, mm in REG.items()}
        for k in PRIM:
            d1 = P1[k] - P0[k]; dh = Ph[k] - P0[k]
            for rg in ("contraction", "contraction_wall1_10"):
                mm = REG[rg]; rec[rg][f"prop_err_{k}"] = nrm(dh - 0.5 * d1, mm) / max(nrm(0.5 * d1, mm), 1e-300)
        out["arms"][arm][path] = rec
c = {a: (out["arms"][a]["old"]["contraction"], out["arms"][a]["new"]["contraction"]) for a in ("GI", "GB")}
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif any(c[a][0]["norm_halfdR"] <= 0 for a in c):
    res = "判定不能 (信号不足)"
elif all(c[a][1]["ell"] <= 0.05 and c[a][1]["normN"] <= 0.1 * c[a][0]["normN"] for a in c):
    res = "A: 熱物性の評価の経路の切り替え (thermoFloat 1 → 0) で比例しない成分の大半が消えた (この状態・補正の方向で)"
elif all(c[a][1]["normN"] >= 0.9 * c[a][0]["normN"] for a in c):
    res = "B: thermoFloat 0 でも比例しない成分が 9 割以上残る (この切り替えで除ける誤差が主因という説を退ける)"
else:
    res = "判定不能 (中間・腕で割れる)"
out["qualification"] = qual; out["result"] = res
(OUTD / "tf_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for a in ("GI", "GB"):
    for path in ("old", "new"):
        for rg, v in out["arms"][a][path].items():
            extra = "  ".join(f"{k[9:]} {v[k]:.2e}" for k in v if k.startswith("prop_err_"))
            print(f"{a} {path:3s} {rg:22s} ‖N‖ {v['normN']:.3e}  ℓ {v['ell']:.3f}  ‖½ΔR‖ {v['norm_halfdR']:.3e}  {extra}")
print("旧の Q0 の GI と GB の差:", q0share); print("資格:", qual); print("結果:", res)
