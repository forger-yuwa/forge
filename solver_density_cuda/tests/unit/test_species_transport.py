#!/usr/bin/env python3
"""種ごとの輸送物性の出所 (physProp.transport) の段 1 試験 — CPU resolver・単成分・独立参照・記録 (GPU 不要)。
plans/active/thermophysics-solver-owned-species-db.md §5.1 #5t2 段 1 (合格条件は実装前に固定したもの)。

  python3 solver_density_cuda/tests/unit/test_species_transport.py [--forge BIN] [--base-forge BASE_BIN] [--case44 DIR] [--keep]

  (G)  生成データ forge_transport_v1.yaml が trans.inp からの再生成と一字一句同じ (生成器 --check)。trans.inp の区間が連続
       (参照側の区間の選び方 = mixing_ab.py 流と実装側 = CEA の kt が同じになる条件)。
  (F)  全 CEA 指定の N2–H2O 16 状態 (T 400/600/1000/2000 K × X_H2O 0/0.1/0.5/1) が FCEA2 (only N2 H2O; frozen) と μ・λ ≤0.1 %。
       FCEA2 はこの場で notes/investigations/2026-09-27-cea-vs-forge-properties/fcea_mix/*.inp を実行する。
  (R)  選択モデル (kinetic・custom:h2o_iapws_cea_v1・fit・混在・lump の重複実種・mass basis・別名) が独立参照
       (transport_reference.py; 実装とコードを共有しない) と double で相対 ≤1e-12 (μ, λ, 実種ごとの μ_i, λ_i, η_ij, 展開後 X)。
  (AB) 400 K 純 H2O で出所だけ cea / custom:h2o_iapws_cea_v1 に変え、それぞれ CEA 値 (μ 1.32788714104e−5・λ 0.0270414202095)・
       IAPWS 値 (μ 1.33545407126e−5・λ 0.0264314431570) と一致 (codex の独立計算値; 12 桁で与えられているので許容は
       max(1e-12, 与えられた桁の丸め幅) で、独立参照とは ≤1e-12)。両者が同値なら分岐が効いていない不良。
  (J)  H2O custom の 500 K・700 K・150 K で、左右極限の値の相対差 ≤1e-12、無次元勾配 T·d ln f/dT の差 ≤1e-10 (μ・λ)。
       極限値は nextafter の両側、勾配は各側の 7 点片側差分 (h: 150 K で 0.05 K、500/700 K で 1 K)。
  (N)  拒否: 未指定種 (単独・lump 構成種) / lump 名をキーに / 未知の種 / 同じ実種に 2 回 (AR と Ar) / 未知モデル (大文字 CEA 含む) /
       未知の custom / custom を H2O 以外に / CEA データの無い種で cea / 組成指定の無い builtin で fit / transport_fit の不正
       (区間 4 つ・Tlo≥Thi・C 欠落・未知キー) / LJ を明示しない外部 DB 種で kinetic。拒否漏れ 0 件。
  (C)  (--forge) forge --resolve-species: 記録が schema v2 + transport_compat を持ち Python load_record が互換性ハッシュを再計算できる、
       輸送指定だけ違う記録の差 (_record_diff) が transport 行を示す、計算の起動は止まる、thermalMethod≠2 と非 mapping は拒否。
  (H)  (--forge と --base-forge) physProp.transport の無い config のハッシュと記録ファイルがバイト不変:
       case/44 run_0509 (4378b7d78339ba27) と内蔵のみ 3 構成。
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, json, math, os, re, shutil, subprocess, sys, tempfile

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.normpath(os.path.join(HERE, "..", ".."))
REPO = os.path.normpath(os.path.join(SOLVER, ".."))
TOOLS = os.path.join(SOLVER, "tools")
sys.path.insert(0, HERE)
sys.path.insert(0, TOOLS)
import forge_species as fs          # noqa: E402
from transport_reference import Reference, TR, iapws_mu, iapws_lam, _loglog_slope, H_IAPWS, L_IAPWS  # noqa: E402

FCEA_DIR = os.path.join(REPO, "notes", "investigations", "2026-09-27-cea-vs-forge-properties", "fcea_mix")
CEA = os.path.join(REPO, ".venv-cea", "nasa_cea")
FAIL = 0


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what, flush=True)
    if not ok:
        FAIL += 1


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


class Driver:
    def __init__(self, work):
        gen = os.path.join(work, "gen")
        subprocess.run(["cmake", f"-DIN={os.path.join(SOLVER, 'data', 'species', 'forge_species_v1.yaml')}",
                        f"-DOUT={os.path.join(gen, 'forge_species_data.hpp')}", "-P",
                        os.path.join(SOLVER, "cmake", "embed_species_data.cmake")], check=True, capture_output=True)
        self.exe = os.path.join(work, "transport_eval_host")
        subprocess.run(["g++", "-O1", "-std=c++17", "-Wno-unknown-pragmas", "-I", SOLVER, "-I", gen,
                        os.path.join(HERE, "transport_eval_host.cpp"), os.path.join(SOLVER, "input", "speciesDB.cpp"),
                        os.path.join(SOLVER, "input", "speciesTransportDB.cpp"), "-lyaml-cpp", "-o", self.exe], check=True)
        self.work = work
        self.n = 0

    def run(self, spec, db=None):
        self.n += 1
        d = os.path.join(self.work, f"d{self.n:03d}")
        os.makedirs(d)
        spec = dict(spec)
        if db is not None:
            with open(os.path.join(d, "db.yaml"), "w") as f:
                yaml.safe_dump(db, f, sort_keys=False)
            spec["speciesDBFile"] = os.path.join(d, "db.yaml")
        p = os.path.join(d, "spec.yaml")
        with open(p, "w") as f:
            yaml.safe_dump(spec, f, sort_keys=False)
        r = subprocess.run([self.exe, p], capture_output=True, text=True)
        return json.loads(r.stdout)


def n2_db_entry(**extra):
    """外部 DB の試験種 (熱物性は内蔵 N2 と同じ係数)。"""
    n2 = [e for e in yaml.safe_load(open(os.path.join(SOLVER, "data", "species", "forge_species_v1.yaml")))["species"] if e["id"] == "N2"][0]
    e = {"MW": n2["MW"], "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
         "nasa9_low": n2["intervals"][0]["coeffs"], "nasa9_high": n2["intervals"][1]["coeffs"]}
    e.update(extra)
    return e


# ---------------------------------------------------------------- (G)
def test_generator():
    r = subprocess.run([sys.executable, os.path.join(TOOLS, "cea_trans_to_forge_transport.py"), "--check"], capture_output=True, text=True)
    check(r.returncode == 0, "G forge_transport_v1.yaml == regeneration from trans.inp (" + r.stdout.strip()[-60:] + ")")
    gaps = []
    for k, rec in TR.items():
        for kind in ("V", "C"):
            iv = rec[kind]
            for a, b in zip(iv, iv[1:]):
                if a[1] != b[0]:
                    gaps.append((k, kind, a[1], b[0]))
    check(not gaps, f"G trans.inp intervals are contiguous (reference and CEA interval selection coincide): gaps {gaps[:3]}")


# ---------------------------------------------------------------- (F)
def run_fcea(work):
    out = {}
    d = os.path.join(work, "fcea")
    os.makedirs(d)
    for f in ("thermo.lib", "trans.lib"):
        shutil.copy(os.path.join(CEA, f), d)
    for xw, name in ((0.0, "mix_0"), (0.1, "mix_0p1"), (0.5, "mix_0p5"), (1.0, "mix_1")):
        shutil.copy(os.path.join(FCEA_DIR, name + ".inp"), d)
        subprocess.run([os.path.join(CEA, "FCEA2")], input=name + "\n", cwd=d, capture_output=True, text=True, check=True)
        t = open(os.path.join(d, name + ".out")).read()
        T = [float(x) for x in re.search(r"T, K\s+([\d. ]+)", t).group(1).split()]
        V = [float(x) for x in re.search(r"VISC,MILLIPOISE\s+([\d. ]+)", t).group(1).split()]
        C = [float(x) for x in re.findall(r"^ CONDUCTIVITY\s+([\d. ]+)$", t, re.M)[-1].split()]   # 最後 = FROZEN 節
        out[xw] = list(zip(T, V, C))
    return out


def test_fcea(D, work):
    if not os.path.exists(os.path.join(CEA, "FCEA2")):
        check(False, "F FCEA2 not found at .venv-cea/nasa_cea/FCEA2")
        return
    ref = run_fcea(work)
    worst_m = worst_l = 0.0
    rows = []
    for xw, pts in ref.items():
        spec = {"species": ["N2", "H2O"], "transport": {"N2": "cea", "H2O": "cea"},
                "states": [{"T": T, "X": [1.0 - xw, xw]} for T, _, _ in pts]}
        d = D.run(spec)
        for (T, v, c), st in zip(pts, d["states"]):
            em = st["mu"] / (v * 1e-4) - 1.0            # millipoise -> Pa s
            el = st["lam"] / (c * 1e-1) - 1.0           # mW/(cm K) -> W/(m K)
            worst_m, worst_l = max(worst_m, abs(em)), max(worst_l, abs(el))
            rows.append(f"X_H2O {xw:4.2f} T {T:6.0f}  mu {100 * em:+.4f} %  lam {100 * el:+.4f} %")
    for r in rows:
        print("       " + r)
    check(len(rows) == 16 and worst_m <= 1e-3 and worst_l <= 1e-3,
          f"F all-CEA N2-H2O 16 states vs FCEA2: max |mu| {100 * worst_m:.4f} %, max |lam| {100 * worst_l:.4f} % (<= 0.1 %)")


# ---------------------------------------------------------------- (R)
T_GRID = [100.0, 149.9, 150.0, 200.0, 253.15, 300.0, 373.2, 400.0, 499.99, 500.0, 500.01, 600.0, 699.99, 700.0, 1000.0,
          1073.2, 2000.0, 5000.0, 6000.0, 12000.0, 20000.0]


def test_reference(D):
    fit_db = {"N2FIT": n2_db_entry(LJ_sigma=3.621, LJ_eps_kB=97.53,
                                   transport_fit={"V": [[150.0, 1000.0, 0.62526577, -31.779652, -1640.7983, 1.7454992],
                                                        [1000.0, 7000.0, 0.87395209, 561.52222, -173948.09, -0.39335958]],
                                                  "C": [[150.0, 1000.0, 0.85439436, 105.73224, -12347.848, 0.47793128]],
                                                  "reference": "test: CEA N2 low interval"}),
              "N2LJ": n2_db_entry(LJ_sigma=3.7, LJ_eps_kB=90.0, LJ_dipole=0.5)}
    cases = [
        ("pure H2O cea", {"species": ["H2O"], "transport": {"H2O": "cea"}}, None, [[1.0]]),
        ("pure H2O custom", {"species": ["H2O"], "transport": {"H2O": "custom:h2o_iapws_cea_v1"}}, None, [[1.0]]),
        ("all kinetic N2/O2/H2O (polar H2O, CE pairs)",
         {"species": ["N2", "O2", "H2O"], "transport": {"N2": "kinetic", "O2": "kinetic", "H2O": "kinetic"}}, None, [[0.7, 0.2, 0.1]]),
        ("N2 kinetic + H2O custom (CEA interaction)",
         {"species": ["N2", "H2O"], "transport": {"N2": "kinetic", "H2O": "custom:h2o_iapws_cea_v1"}}, None, [[0.9, 0.1], [0.5, 0.5], [0.0, 1.0]]),
        ("N2 cea + H2O custom",
         {"species": ["N2", "H2O"], "transport": {"N2": "cea", "H2O": "custom:h2o_iapws_cea_v1"}}, None, [[0.98, 0.02], [0.5, 0.5]]),
        ("AIR kinetic + H2O cea + Ar cea (rigid sphere pairs; alias AR, WATER)",
         {"species": ["AIR", "WATER", "Ar"], "transport": {"AIR": "kinetic", "H2O": "cea", "AR": "cea"}}, None, [[0.8, 0.1, 0.1]]),
        ("two lumps sharing N2 (+ mass basis) + H2O custom",
         {"species": [{"name": "MIX1", "lump": {"N2": 0.8, "O2": 0.2}, "basis": "mole"},
                      {"name": "MIX2", "lump": {"N2": 0.3, "AR": 0.5, "CO2": 0.2}, "basis": "mass"}, "H2O"],
          "transport": {"N2": "cea", "O2": "kinetic", "AR": "cea", "CO2": "kinetic", "H2O": "custom:h2o_iapws_cea_v1"}},
         None, [[0.5, 0.4, 0.1], [0.0, 0.9, 0.1]]),
        ("lump + the same real species listed alone (N2 twice)",
         {"species": [{"name": "MIXDRY", "lump": {"N2": 0.78, "O2": 0.21, "AR": 0.01}, "basis": "mole"}, "N2", "H2O"],
          "transport": {"N2": "kinetic", "O2": "kinetic", "AR": "kinetic", "H2O": "cea"}}, None, [[0.6, 0.3, 0.1]]),
        ("user fit + kinetic polar DB species + cea",
         {"species": ["N2FIT", "N2LJ", "CO2"], "transport": {"N2FIT": "fit", "N2LJ": "kinetic", "CO2": "cea"}}, fit_db,
         [[0.4, 0.4, 0.2]]),
        ("He cea + CO2 kinetic + O2 cea",
         {"species": ["He", "CO2", "O2"], "transport": {"He": "cea", "CO2": "kinetic", "O2": "cea"}}, None, [[0.3, 0.3, 0.4]]),
    ]
    for tag, spec, db, Xs in cases:
        d = D.run({**spec, "states": [{"T": T, "X": X} for X in Xs for T in T_GRID]}, db=db)
        if not d.get("ok"):
            check(False, f"R {tag}: driver failed: {d.get('error')}")
            continue
        ref = Reference(spec["species"], spec["transport"], db)
        worst, where, k = 0.0, "", 0
        for X in Xs:
            for T in T_GRID:
                st, r = d["states"][k], ref.state(T, X)
                k += 1
                for key in ("mu", "lam"):
                    e = rel(st[key], r[key])
                    if e > worst:
                        worst, where = e, f"{key} T={T} X={X}"
                for key in ("mu_i", "lam_i", "eta_ij", "Xreal"):
                    for i, (a, b) in enumerate(zip(st[key], r[key])):
                        e = rel(a, b) if (a or b) else 0.0
                        if e > worst:
                            worst, where = e, f"{key}[{i}] T={T}"
                    if len(st[key]) != len(r[key]):
                        worst, where = 1.0, f"{key} length {len(st[key])} vs {len(r[key])}"
        check(worst <= 1e-12, f"R {tag}: vs independent reference, {k} states, max rel {worst:.2e} ({where}); pairs {d.get('pair_source')}")


# ---------------------------------------------------------------- (AB)
def test_ab(D):
    lit = {"cea": (1.32788714104e-5, 0.0270414202095), "custom:h2o_iapws_cea_v1": (1.33545407126e-5, 0.0264314431570)}
    got = {}
    for m, (mu0, la0) in lit.items():
        d = D.run({"species": ["H2O"], "transport": {"H2O": m}, "states": [{"T": 400.0, "X": [1.0]}]})
        st = d["states"][0]
        r = Reference(["H2O"], {"H2O": m}).state(400.0, [1.0])
        got[m] = (st["mu"], st["lam"])
        # 12 有効桁で与えられた値の丸め幅 (最後の桁の半分) を許容に含める
        tol = lambda v: max(1e-12, 0.5 * 10 ** (math.floor(math.log10(abs(v))) - 11) / abs(v))
        ok_lit = rel(st["mu"], mu0) <= tol(mu0) and rel(st["lam"], la0) <= tol(la0)
        ok_ref = rel(st["mu"], r["mu"]) <= 1e-12 and rel(st["lam"], r["lam"]) <= 1e-12
        check(ok_lit and ok_ref and d["model"] == [m],
              f"AB 400 K pure H2O {m}: mu {st['mu']!r} (codex {mu0}, rel {rel(st['mu'], mu0):.2e} tol {tol(mu0):.1e}; ref rel {rel(st['mu'], r['mu']):.1e}), "
              f"lam {st['lam']!r} (codex {la0}, rel {rel(st['lam'], la0):.2e} tol {tol(la0):.1e}; ref rel {rel(st['lam'], r['lam']):.1e})")
    a, b = got["cea"], got["custom:h2o_iapws_cea_v1"]
    check(rel(a[0], b[0]) > 1e-3 and rel(a[1], b[1]) > 1e-2,
          f"AB branch is effective: custom/cea - 1 = mu {100 * (b[0] / a[0] - 1):+.4f} %, lam {100 * (b[1] / a[1] - 1):+.4f} % (expected +0.57 %, -2.26 %)")


# ---------------------------------------------------------------- (J)
def _weights(m):
    """前進片側差分 x0 + k h (k=0..m) の 1 階微分の重み (m 次精度; Vandermonde を有理数で解く)。"""
    from fractions import Fraction as Fr
    A = [[Fr(k) ** j for k in range(m + 1)] for j in range(m + 1)]
    b = [Fr(0)] * (m + 1)
    b[1] = Fr(1)
    n = m + 1
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [float(M[i][n] / M[i][i]) for i in range(n)]


def test_joins(D):
    W = _weights(6)

    def vals(Ts):
        d = D.run({"species": ["H2O"], "transport": {"H2O": "custom:h2o_iapws_cea_v1"}, "states": [{"T": T, "X": [1.0]} for T in Ts]})
        return [(s["mu"], s["lam"]) for s in d["states"]]

    ceaV, ceaC = TR[("H2O", "")]["V"], TR[("H2O", "")]["C"]

    def cea_slope(rows, T):
        A, B, C, _ = next(r for r in rows if r[0] <= T <= r[1])[2]
        return A - B / T - 2 * C / T**2

    expect = {150.0: (_loglog_slope(H_IAPWS, 150.0), _loglog_slope(L_IAPWS, 150.0)),
              500.0: (_loglog_slope(H_IAPWS, 500.0), _loglog_slope(L_IAPWS, 500.0)),
              700.0: (cea_slope(ceaV, 700.0), cea_slope(ceaC, 700.0))}
    for T0, h in ((150.0, 0.05), (500.0, 1.0), (700.0, 1.0)):
        lo, hi = math.nextafter(T0, 0.0), math.nextafter(T0, math.inf)
        (ml, ll), (mr, lr) = vals([lo, hi])
        dv = max(rel(ml, mr), rel(ll, lr))
        L = vals([T0 - k * h for k in range(7)])
        R = vals([T0 + k * h for k in range(7)])
        res = []
        for q, nm in ((0, "mu"), (1, "lam")):
            gl = -sum(w * math.log(v[q]) for w, v in zip(W, L)) / h * T0
            gr = sum(w * math.log(v[q]) for w, v in zip(W, R)) / h * T0
            res.append((nm, gl, gr, abs(gr - gl), abs(gl - expect[T0][q])))
        ok = dv <= 1e-12 and all(r[3] <= 1e-10 for r in res)
        check(ok, f"J {T0:.0f} K: value left/right rel {dv:.1e} (<=1e-12); "
                  + "; ".join(f"{nm} T dlnf/dT left {gl:.12f} right {gr:.12f} |diff| {dg:.1e} (<=1e-10), vs analytic {da:.1e}"
                              for nm, gl, gr, dg, da in res))


# ---------------------------------------------------------------- (N)
def test_negative(D):
    base_db = {"NOLJ": n2_db_entry(), "N2FIT": n2_db_entry(LJ_sigma=3.621, LJ_eps_kB=97.53)}

    def fitdb(tf):
        return {"N2FIT": n2_db_entry(LJ_sigma=3.621, LJ_eps_kB=97.53, transport_fit=tf)}
    good_rows = [[200.0, 1000.0, 0.6, -30.0, -1600.0, 1.7]]
    lump = {"name": "MIXDRY", "lump": {"N2": 0.79, "O2": 0.21}, "basis": "mole"}
    cases = [
        ("missing species", ["N2", "H2O"], {"N2": "cea"}, None, "no transport model for species H2O"),
        ("missing lump constituent", [lump, "H2O"], {"N2": "cea", "H2O": "cea"}, None, "no transport model for species O2"),
        ("lump name as key", [lump, "H2O"], {"MIXDRY": "cea", "H2O": "cea"}, None, "is a lump"),
        ("unknown species key", ["N2"], {"N2": "cea", "CO2": "cea"}, None, "is not a species of this run"),
        ("same species twice (AR and Ar)", ["Ar"], {"AR": "cea", "Ar": "kinetic"}, None, "is given twice"),
        ("unknown model", ["N2"], {"N2": "sutherland"}, None, "unknown transport model"),
        ("model name is case-sensitive (CEA)", ["N2"], {"N2": "CEA"}, None, "unknown transport model"),
        ("unknown custom", ["H2O"], {"H2O": "custom:h2o_iapws_cea_v2"}, None, "unknown custom model"),
        ("custom on N2", ["N2", "H2O"], {"N2": "custom:h2o_iapws_cea_v1", "H2O": "cea"}, None, "defined only for H2O"),
        ("cea without CEA data (AIR)", ["AIR"], {"AIR": "cea"}, None, "has no entry in CEA trans.inp"),
        ("fit on a built-in species", ["N2"], {"N2": "fit"}, None, "needs 'transport_fit"),
        ("fit without transport_fit", ["N2FIT"], {"N2FIT": "fit"}, base_db, "needs 'transport_fit"),
        ("fit with 4 intervals", ["N2FIT"], {"N2FIT": "fit"},
         fitdb({"V": good_rows * 1 + [[1000.0, 2000.0, 0, 0, 0, 1], [2000.0, 3000.0, 0, 0, 0, 1], [3000.0, 4000.0, 0, 0, 0, 1]], "C": good_rows}),
         "exceed TRANSPORT_MAX_FIT_INTERVALS"),
        ("fit with Tlo >= Thi", ["N2FIT"], {"N2FIT": "fit"}, fitdb({"V": [[1000.0, 200.0, 0, 0, 0, 1]], "C": good_rows}), "needs 0 < Tlo < Thi"),
        ("fit without C", ["N2FIT"], {"N2FIT": "fit"}, fitdb({"V": good_rows}), "no intervals"),
        ("fit with unknown key", ["N2FIT"], {"N2FIT": "fit"}, fitdb({"V": good_rows, "C": good_rows, "units": "SI"}), "unknown key 'units'"),
        ("fit with 5 numbers per interval", ["N2FIT"], {"N2FIT": "fit"}, fitdb({"V": [[200.0, 1000.0, 0.6, -30.0, 1.7]], "C": good_rows}),
         "each interval must be"),
        ("kinetic without explicit LJ (external DB)", ["NOLJ"], {"NOLJ": "kinetic"}, base_db, "no explicit LJ_sigma/LJ_eps_kB"),
    ]
    missed = 0
    for tag, species, tr, db, needle in cases:
        d = D.run({"species": species, "transport": tr, "states": []}, db=db)
        ok = (not d.get("ok")) and needle in d.get("error", "")
        missed += 0 if ok else 1
        check(ok, f"N {tag}: rejected ('{needle}')" + ("" if ok else f" got {d}"))
    check(missed == 0, f"N rejection misses: {missed} / {len(cases)}")
    # 正例: 書かなければ従来どおり (輸送ブロックなし・schema v1)
    d = D.run({"species": ["N2", "H2O"], "states": []})
    check(d["ok"] and not d["transport_enabled"] and '"forge_resolved_species_v1"' in d["record"] and "transport_compat" not in d["record"],
          "N without physProp.transport: no transport block, record schema v1")


# ---------------------------------------------------------------- (C)/(H)
class Cfg:
    def __init__(self, root, seed):
        self.root, self.seed, self.n = root, seed, 0

    def make(self, species, transport=None, db=None, extra=None, keep_db=False):
        self.n += 1
        d = os.path.join(self.root, f"c{self.n:02d}")
        os.makedirs(d)
        t = open(os.path.join(self.seed, "solverConfig.yaml")).read()
        sp = "species: " + json.dumps(species) + ("" if transport is None else ", transport: " + (
            transport if isinstance(transport, str) else json.dumps(transport)))
        t, k = re.subn(r"species: *\[[^\]]*\]", sp, t, count=1)
        assert k == 1
        if not keep_db:
            t = re.sub(r",? *speciesDBFile: *\"?[^,}\"]*\"?", "", t, count=1)
        else:
            shutil.copy(os.path.join(self.seed, "species_db.yaml"), d)
        if db is not None:
            t = t.replace("thermoHrefTemp:", "speciesDBFile: \"test_db.yaml\", thermoHrefTemp:", 1)
            with open(os.path.join(d, "test_db.yaml"), "w") as f:
                yaml.safe_dump(db, f, sort_keys=False)
        for pat, rep in (extra or []):
            t = re.sub(pat, rep, t)
        with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
            f.write(t)
        return d

    def copy(self):
        self.n += 1
        d = os.path.join(self.root, f"c{self.n:02d}")
        os.makedirs(d)
        shutil.copy(os.path.join(self.seed, "solverConfig.yaml"), d)
        shutil.copy(os.path.join(self.seed, "species_db.yaml"), d)
        return d


def resolve(exe, d):
    p = subprocess.run([exe, "--resolve-species"], cwd=d, capture_output=True, text=True)
    lines = p.stdout.strip().splitlines()
    h = lines[-1].strip() if lines else ""
    m = re.search(r"\[species\] record (\S+) \(sha256 ([0-9a-f]{64})\)", p.stderr)
    path = os.path.join(d, m.group(1)) if (p.returncode == 0 and m) else None
    return p.returncode, h, path, p.stderr


def test_forge(a, root):
    C = Cfg(root, a.case44)
    tr = {"N2": "cea", "O2": "cea", "AR": "kinetic", "CO2": "cea", "H2O": "custom:h2o_iapws_cea_v1"}
    lump = {"name": "MIXDRY", "lump": {"N2": 0.708873, "O2": 0.230376, "AR": 0.00850387, "CO2": 0.0522474}, "basis": "mole"}
    d1 = C.make([lump, "H2O"], tr)
    rc, h, path, err = resolve(a.forge, d1)
    rec = fs.load_record(path) if path else None
    ok = (rc == 0 and rec is not None and rec["consistent"] and rec["compat_recomputed"] == h
          and rec["schema"] == "forge_resolved_species_v2" and any(l.startswith("transport.expand[0]") for l in rec["transport_compat"]))
    check(ok, f"C --resolve-species with physProp.transport: schema v2, transport_compat {len(rec['transport_compat']) if rec else 0} lines, "
              f"Python load_record recomputes the hash ({h[:16]})" + ("" if ok else err[-1500:]))
    check("[species]   transport (physProp.transport" in err and "custom:h2o_iapws_cea_v1" in err and "outside their formal range" in err,
          "C startup log shows the transport table and the IAPWS range note")
    # 出所だけ違う (H2O: cea) → ハッシュが違い、差は transport 行
    d2 = C.make([lump, "H2O"], dict(tr, H2O="cea"))
    rc2, h2, path2, _ = resolve(a.forge, d2)
    rec2 = fs.load_record(path2) if path2 else None
    diff = fs._record_diff(rec, rec2) if (rec and rec2) else ["(no record)"]
    check(rc2 == 0 and h2 != h and len(diff) == 1 and diff[0].startswith("transport:") and "h2o_iapws" in diff[0],
          f"C only the H2O model differs -> hash differs ({h[:16]} vs {h2[:16]}) and the diff names the transport line: {diff}")
    # 記録の輸送行を書き換えると自己整合が崩れる
    if path:
        t = open(path).read().replace("model=custom:h2o_iapws_cea_v1", "model=cea", 1)
        bad = os.path.join(root, "resolved_species_edit.yaml")
        open(bad, "w").write(t)
        check(not fs.load_record(bad)["consistent"], "C edited transport line in a record -> load_record reports compat mismatch")
    # 計算の起動は止める (段 2 まで)
    p = subprocess.run([a.forge], cwd=d1, capture_output=True, text=True, timeout=120)
    check(p.returncode != 0 and "GPU transport path is not connected yet" in (p.stderr + p.stdout),
          "C running the solver with physProp.transport stops before GPU init (stage 1)")
    # config の拒否
    d3 = C.make(["N2"], {"N2": "cea"}, extra=[(r"thermalMethod: *2", "thermalMethod: 1")])
    rc3, _, _, err3 = resolve(a.forge, d3)
    p3 = subprocess.run([a.forge], cwd=d3, capture_output=True, text=True, timeout=120)
    check(p3.returncode != 0 and "physProp.transport requires thermalMethod: 2" in (p3.stderr + p3.stdout),
          "C physProp.transport with thermalMethod 1 is rejected")
    d4 = C.make(["N2"], "[cea]")
    rc4, _, _, err4 = resolve(a.forge, d4)
    check(rc4 != 0 and "physProp.transport must be a non-empty mapping" in err4, "C physProp.transport as a list is rejected")
    d5 = C.make(["N2", "H2O"], {"N2": "cea"})
    rc5, _, _, err5 = resolve(a.forge, d5)
    check(rc5 != 0 and "no transport model for species H2O" in err5, "C --resolve-species rejects an unspecified species")

    # (H) physProp.transport の無い config はバイト不変
    cfgs = [("case44 run_0509 (external DB MIXDRY/H2O)", C.copy())]
    for sp in (["N2"], ["N2", "H2O"], ["H2O", "N2", "O2", "AR", "CO2"]):
        cfgs.append((f"builtin {sp}", C.make(sp)))
    for tag, dd in cfgs:
        rcn, hn, pn, _ = resolve(a.forge, dd)
        if tag.startswith("case44"):
            check(rcn == 0 and hn.startswith("4378b7d78339ba27"), f"H {tag}: hash {hn[:16]} == 4378b7d78339ba27")
        if a.base_forge:
            db_ = dd + "_base"
            shutil.copytree(dd, db_, ignore=shutil.ignore_patterns("resolved_species_*"))
            rcb, hb, pb, _ = resolve(a.base_forge, db_)
            same = rcn == 0 and rcb == 0 and hn == hb and pn and pb and open(pn, "rb").read() == open(pb, "rb").read()
            check(same, f"H {tag}: hash {hn[:16]} and record bytes identical to base binary ({hb[:16]})")
        else:
            print(f"[SKIP] H {tag}: --base-forge not given")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", default=os.environ.get("FORGE_BIN"))
    ap.add_argument("--base-forge", default=None)
    ap.add_argument("--case44", default=os.path.join(REPO, "case", "44.vitiated_air_wt", "run_0509_va3_M4.19_Lc8_dry_lumpX"),
                    help="solverConfig.yaml と species_db.yaml の雛形 (読むだけ)")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    root = tempfile.mkdtemp(prefix="forge_species_transport_")
    try:
        D = Driver(root)
        test_generator()
        test_fcea(D, root)
        test_reference(D)
        test_ab(D)
        test_joins(D)
        test_negative(D)
        if a.forge:
            a.forge = os.path.abspath(a.forge)
            if a.base_forge:
                a.base_forge = os.path.abspath(a.base_forge)
            if not os.path.exists(os.path.join(a.case44, "solverConfig.yaml")):
                check(False, f"C/H seed config not found: {a.case44}")
            else:
                test_forge(a, root)
        else:
            print("[SKIP] C/H: --forge not given")
    finally:
        if a.keep:
            print("kept", root)
        else:
            shutil.rmtree(root, ignore_errors=True)
    print("ALL PASS" if FAIL == 0 else f"FAILED: {FAIL}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
