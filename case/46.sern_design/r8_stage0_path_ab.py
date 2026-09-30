#!/usr/bin/env python3
"""SERN R8 段 (i) の照合 (0) の判別 A/B (codex diagnose 2026-09-30、plan tooling-nozzle-sern-chain R8)。CFD は回さない。
実 run の保存済み物性を、合成の演算経路だけを変えて再計算し再現できるかを見る。
  A = 旧経路: 旧コード (~/forge-pgrad-new/design) の composition.lump_entry を m6_on の EXH に適用 → 旧 run の species_db.yaml と照合
  B = 新経路: C++ speciesDB.cpp synthesizeLump の逐次演算 (x = config の lump 内モル分率 / Σ、MW = Σ x M、a = Σ x a) を
      新 run の config と構成種の生データ (内蔵 forge_species_v1.yaml + species_db_external.yaml) で再計算 → 新 run の
      resolved_species_*.yaml と照合
補助: 保存値同士で EXH の cp の最大相対差と、顕内部エネルギー差 / (cp·max(T, 298.15)) を 50–7000 K の 1007 点で (≤ 1e-12)。
  python3 r8_stage0_path_ab.py OLD_RUN NEW_RUN      (AWS で、OLD は ~/forge-pgrad-new の run_1002、NEW は ~/forge-r8 の run_1005)
"""
import glob, os, subprocess, sys, json
import numpy as np, yaml

OLD, NEW = sys.argv[1], sys.argv[2]
old_db = yaml.safe_load(open(os.path.join(OLD, "species_db.yaml")))
rec_f = sorted(glob.glob(os.path.join(NEW, "resolved_species_*.yaml")))
if not rec_f:
    sys.exit(f"{NEW} に resolved_species_*.yaml が無い")
rec = yaml.safe_load(open(rec_f[0]))
print("新 run の記録:", os.path.basename(rec_f[0]), " keys:", list(rec.keys())[:8])


def rec_species(rec, name):
    for key in ("species", "resolved"):
        if key in rec:
            v = rec[key]
            if isinstance(v, dict) and name in v:
                return v[name]
            if isinstance(v, list):
                for e in v:
                    if e.get("name") == name:
                        return e
    raise KeyError(name)


def coeffs(e):
    lo = e.get("nasa9_low", e.get("low")); hi = e.get("nasa9_high", e.get("high"))
    return float(e["MW"]), [float(v) for v in lo], [float(v) for v in hi]


# ---- A: 旧経路 (旧コードの lump_entry を旧 run と同じ入力で) ----
code = r'''
import sys, json
sys.path.insert(0, "/home/ubuntu/forge-pgrad-new/design")
from forge_design.evaluate import runner_sern as R2
from forge_design.probdef import load_problem
p = load_problem("/home/ubuntu/forge-pgrad-new/case/46.sern_design/problem_moo_frozen_tp_cycle3op.yaml")
R2.design_snapshot(p); R2.select_operating_point(p, "m6_on")
e = R2.frozen_gases(p)["layout"].entries["EXH"]
print(json.dumps({"MW": e.MW, "low": list(e.low), "high": list(e.high)}))
'''
A = json.loads(subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.strip().splitlines()[-1])
oMW, olo, ohi = coeffs(old_db["EXH"])
okA = A["MW"] == oMW and A["low"] == olo and A["high"] == ohi
print(f"[A 旧経路] 再計算 vs 旧 run species_db.yaml: {'ビット一致' if okA else '不一致'}")

# ---- B: 新経路 (C++ の逐次演算を再現) ----
cfg = yaml.safe_load(open(os.path.join(NEW, "solverConfig.yaml")))
lump = [s for s in cfg["physProp"]["species"] if isinstance(s, dict) and s["name"] == "EXH"][0]
members = list(lump["lump"].keys()); fr = [float(lump["lump"][k]) for k in members]
builtin = yaml.safe_load(open("/home/ubuntu/forge-r8/solver_density_cuda/data/species/forge_species_v1.yaml"))
bmap = {}
for e in builtin["species"]:
    for nm in [e["id"]] + list(e.get("aliases") or []):
        bmap[str(nm).upper()] = e
ext = yaml.safe_load(open(os.path.join(NEW, "species_db_external.yaml"))) if os.path.exists(os.path.join(NEW, "species_db_external.yaml")) else {}


def member(k):
    if k in ext:
        e = ext[k]; return float(e["MW"]), [float(v) for v in e["nasa9_low"]], [float(v) for v in e["nasa9_high"]], (e["Tlo"], e["Tmid"], e["Thi"])
    e = bmap[k.upper()]; iv = e["intervals"]
    return float(e["MW"]), [float(v) for v in iv[0]["coeffs"]], [float(v) for v in iv[1]["coeffs"]], (iv[0]["Tlo"], iv[0]["Thi"], iv[1]["Thi"])


ms = [member(k) for k in members]
s = 0.0
for v in fr:
    s += v
x = [v / s for v in fr]
MW = 0.0
for k in range(len(ms)):
    MW += x[k] * ms[k][0]
lo = [0.0] * 9; hi = [0.0] * 9
for k in range(len(ms)):
    for i in range(9):
        lo[i] += x[k] * ms[k][1][i]; hi[i] += x[k] * ms[k][2][i]
nMW, nlo, nhi = coeffs(rec_species(rec, "EXH"))
okB = MW == nMW and lo == nlo and hi == nhi
print(f"[B 新経路] 再計算 vs 新 run の解決済み記録: {'ビット一致' if okB else '不一致'}"
      + ("" if okB else f"  (MW {MW!r} vs {nMW!r}; low 差 {[a - b for a, b in zip(lo, nlo) if a != b][:3]})"))
print("   構成種の温度区切り:", sorted(set(m[3] for m in ms)))

# ---- 保存値同士の差 (A 保存 vs B 保存) ----
ulp = [abs(a - b) / np.spacing(abs(b)) for a, b in zip(nlo + nhi, olo + ohi) if a != b]
print(f"[保存値の差] MW 相対 {abs(nMW - oMW) / oMW:.2e}、係数の差 {len(ulp)} / 18 個、最大 {max(ulp):.1f} ulp")

# ---- 補助: 物性の差 ----
RU = 8.314462618


def cp_R(a, T):
    return a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4


def h_RT(a, T):
    return -a[0] / T**2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3 + a[5] * T**3 / 4 + a[6] * T**4 / 5 + a[7] / T


T = np.linspace(50.0, 7000.0, 1007)
def props(MWx, l, h):
    a = [np.where(T < 1000.0, l[i], h[i]) for i in range(9)]
    cp = cp_R(a, T) * RU / MWx
    hs = h_RT(a, T) * RU * T / MWx - (h_RT([l[i] for i in range(9)], 298.15) * RU * 298.15 / MWx)
    e = hs - RU / MWx * T
    return cp, e
cpo, eo = props(oMW, olo, ohi); cpn, en = props(nMW, nlo, nhi)
dcp = float(np.max(np.abs(cpn - cpo) / np.abs(cpo)))
de = float(np.max(np.abs(en - eo) / (np.abs(cpo) * np.maximum(T, 298.15))))
print(f"[補助] cp 最大相対差 {dcp:.2e}、顕内部エネルギー差 / (cp·max(T,298.15)) 最大 {de:.2e} (≤ 1e-12: {'OK' if max(dcp, de) <= 1e-12 else 'NG'})")
print("VERDICT:", "差を再現 (A・B ともビット一致) → 演算経路の丸めとして説明できる" if (okA and okB) else "再現できない → 入力・DB・バイナリの対応を調べる")
