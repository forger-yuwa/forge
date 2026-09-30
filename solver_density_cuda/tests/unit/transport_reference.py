#!/usr/bin/env python3
"""physProp.transport の独立参照計算 (plan thermophysics-solver-owned-species-db §5.1 #5t2 段 1)。

実装 (input/speciesTransportDB.cpp・cuda_forge/transportMix_d.cuh・生成データ forge_transport_v1.yaml) とコードを共有しない:
  - CEA 係数は `trans.inp` を直接読む (生成器・埋め込みデータを経由しない; 読み方は notes/investigations/2026-09-27-cea-vs-forge-properties/mixing_ab.py と同じ正規表現)
  - 区間の選び方は mixing_ab.py と同じ「lo ≤ T ≤ hi の最初の区間、下に外れたら最初・上に外れたら最後」
    (CEA の kt と連続な区間で同じ; 試験で区間の連続性を確認する)
  - 混合則は CEA の式を g/mol で書き直したもの (cea2.f TRANP)、剛体球近似は cea2.f 5565–5570
  - IAPWS 希薄気体 (2008 粘性 μ₀・2011 熱伝導 λ₀) の係数は h2o_blend_and_lambda.py と同じ値、
    253.15 K の対数勾配は S'(Tb) を別の形で微分して求める
MW・NASA-9 (修正 Eucken の c_p)・LJ・双極子は実装と同じ共通データ (forge_species_v1.yaml) か、試験の外部 DB から取る
(MW は実装と同じ値にそろえる: H2O 18.0153 g/mol; codex diagnose 2026-09-27 の指摘)。

  from transport_reference import Reference
  ref = Reference(species=[...], transport={...}, db={...外部 DB (dict) か None})
  out = ref.state(T, X)   # {"mu", "lam", "Xreal", "mu_i", "lam_i", "eta_ij"} (SI)
  out = ref.state_Y(T, Y) # 輸送種の質量分率から (段 2 #5t2-2: 負値は 0 に切り、X_s = (Y_s/M_s)/Σ(Y/M) を作ってから state)
輸送種の分子量 (ref.mw[s]) は構成実種の MW と lump 内モル分率から M = Σ x_k M_k として独自に作る (実装の合成 MW を読まない)。
"""
import math
import os
import re

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
TRANS_INP = os.path.join(REPO, ".venv-cea", "nasa_cea", "trans.inp")
SPECIES_DATA = os.path.join(REPO, "solver_density_cuda", "data", "species", "forge_species_v1.yaml")
RU = 8.314462618
KB_CGS = 1.380649e-16


def load_trans(path=TRANS_INP):
    lines = open(path, encoding="latin-1").read().splitlines()
    num = re.compile(r"-?\d\.\d+E[ +-]\d+")
    tr = {}
    i = 1
    while i < len(lines):
        h = lines[i]
        m = re.search(r"V(\d)C(\d)", h)
        if not m:
            i += 1
            continue
        a, b = h[0:15].strip(), h[15:30].strip()
        nv, nc = int(m.group(1)), int(m.group(2))
        rec = {"V": [], "C": []}
        for k in range(nv + nc):
            L = lines[i + 1 + k]
            rec[L[1]].append((float(L[2:11]), float(L[11:20]), [float(x.replace("E ", "E+")) for x in num.findall(L[20:])]))
        tr[(a, b)] = rec
        i += 1 + nv + nc
    return tr


TR = load_trans()


def fit(ivs, T):
    lo = min(ivs, key=lambda v: v[0])
    hi = max(ivs, key=lambda v: v[1])
    use = next((v for v in ivs if v[0] <= T <= v[1]), lo if T < lo[0] else hi)
    A, B, C, D = use[2]
    return math.exp(A * math.log(T) + B / T + C / T**2 + D)


def omega22(Tstar):
    t = min(max(Tstar, 0.3), 100.0)
    return 1.16145 / t**0.14874 + 0.52487 * math.exp(-0.77320 * t) + 2.16178 * math.exp(-2.43787 * t), t


H_IAPWS = [1.67752, 2.20462, 0.6366564, -0.241605]
L_IAPWS = [2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4]
TC = 647.096


def iapws_mu(T):
    Tb = T / TC
    return 100 * math.sqrt(Tb) / sum(H_IAPWS[k] / Tb**k for k in range(4)) * 1e-6


def iapws_lam(T):
    Tb = T / TC
    return math.sqrt(Tb) / sum(L_IAPWS[k] / Tb**k for k in range(5)) * 1e-3


def _loglog_slope(coef, T):
    """f = √Tb / S(Tb), S = Σ c_k Tb^-k → d ln f/d ln T = 1/2 − Tb S'(Tb)/S(Tb)。"""
    Tb = T / TC
    S = sum(c * Tb**(-k) for k, c in enumerate(coef))
    dS = sum(-k * c * Tb**(-k - 1) for k, c in enumerate(coef))
    return 0.5 - Tb * dS / S


def h2o_iapws_cea_v1(T):
    """custom:h2o_iapws_cea_v1 [SI]。"""
    cm = lambda t: fit(TR[("H2O", "")]["V"], t) * 1e-7
    cl = lambda t: fit(TR[("H2O", "")]["C"], t) * 1e-4
    if T >= 700.0:
        return cm(T), cl(T)
    if T < 253.15:
        nm, nl = _loglog_slope(H_IAPWS, 253.15), _loglog_slope(L_IAPWS, 253.15)
        return iapws_mu(253.15) * (T / 253.15)**nm, iapws_lam(253.15) * (T / 253.15)**nl
    if T <= 500.0:
        return iapws_mu(T), iapws_lam(T)
    s = (T - 500.0) / 200.0
    w = 3 * s * s - 2 * s**3
    return (math.exp((1 - w) * math.log(iapws_mu(T)) + w * math.log(cm(T))),
            math.exp((1 - w) * math.log(iapws_lam(T)) + w * math.log(cl(T))))


def _builtin():
    """共通データの全気相種 (ソルバの内蔵と同じ集合; #13-2)。NASA-9 は全区間 (1〜3 区間; #13-3/#13-5(c)) を持つ。
    LJ: null の種は sigma/eps を None にする (kinetic に使うと参照側でも失敗する)。"""
    raw = yaml.safe_load(open(SPECIES_DATA, encoding="utf-8"))
    out, alias = {}, {}
    for e in raw["species"]:
        if e.get("phase") != "gas":
            continue
        iv = e["intervals"]
        lj = e.get("LJ") or {}
        out[e["id"]] = {"MW": float(e["MW"]),
                        "bounds": [float(iv[0]["Tlo"])] + [float(x["Thi"]) for x in iv],
                        "coefs": [[float(x) for x in v["coeffs"]] for v in iv],
                        "sigma": (float(lj["sigma"]) if lj else None), "eps": (float(lj["eps_kB"]) if lj else None),
                        "dipole": float(lj.get("dipole", 0.0)), "file": False}
        alias[e["id"]] = e["id"]
        for a in e.get("aliases") or []:
            alias[a] = e["id"]
    return out, alias


BUILTIN, ALIAS = _builtin()


def cp_mass(sp, T):
    """修正 Eucken の c_p (質量あたり)。範囲外は端でクランプ (cp 一定)。区間は「区切りちょうどは上の区間」
    (ソルバ `thermo_d.cuh` thermo_interval と同じ規約; 6000 K 超の 3 区間種は第 3 区間を使う)。"""
    b = sp["bounds"]
    Tc = min(max(T, b[0]), b[-1])
    k = 0
    for j in range(1, len(b) - 1):
        if not Tc < b[j]:
            k = j
    a = sp["coefs"][k]
    return RU / sp["MW"] * (a[0] / Tc**2 + a[1] / Tc + a[2] + a[3] * Tc + a[4] * Tc**2 + a[5] * Tc**3 + a[6] * Tc**4)


class Reference:
    def __init__(self, species, transport, db=None):
        """species: physProp.species (str か {name, lump, basis})、transport: {キー: モデル}、db: 外部 DB (dict)。"""
        db = db or {}
        self.reals = []      # [{key, sp, model, trans}]
        self.expand = []     # 種 s → [(r, x)]

        def lookup(name):
            if name in db:
                e = db[name]
                return name.upper(), {"MW": float(e["MW"]),
                                      "bounds": [float(e.get("Tlo", 200.0)), float(e.get("Tmid", 1000.0)), float(e.get("Thi", 6000.0))],
                                      "coefs": [[float(x) for x in e["nasa9_low"]], [float(x) for x in e["nasa9_high"]]],
                                      "sigma": float(e.get("LJ_sigma", 3.6)),
                                      "eps": float(e.get("LJ_eps_kB", 97.0)), "dipole": float(e.get("LJ_dipole", 0.0)),
                                      "file": True, "fit": e.get("transport_fit"), "dbkey": name}
            cid = ALIAS[name]
            return cid, dict(BUILTIN[cid], dbkey=cid)

        def add(name, x):
            key, sp = lookup(name)
            for r, e in enumerate(self.reals):
                if e["key"] == key:
                    return r, x
            self.reals.append({"key": key, "sp": sp, "name": name})
            return len(self.reals) - 1, x

        self.mw = []         # 輸送種 s の分子量 [kg/mol] (lump は Σ x_k M_k)
        for s in species:
            if isinstance(s, str):
                self.expand.append([add(s, 1.0)])
                self.mw.append(lookup(s)[1]["MW"])
            else:
                mem = list(s["lump"].items())
                if s["basis"] == "mole":
                    tot = sum(v for _, v in mem)
                    xs = [v / tot for _, v in mem]
                else:
                    tot = sum(v for _, v in mem)
                    mw = [lookup(k)[1]["MW"] for k, _ in mem]
                    y = [v / tot / m for (_, v), m in zip(mem, mw)]
                    xs = [v / sum(y) for v in y]
                self.expand.append([add(k, x) for (k, _), x in zip(mem, xs)])
                self.mw.append(sum(x * lookup(k)[1]["MW"] for (k, _), x in zip(mem, xs)))
        for e in self.reals:
            model = None
            for k, v in transport.items():
                if k == e["name"] or ALIAS.get(k, k.upper()) == e["key"]:
                    model = v
            assert model is not None, e
            e["model"] = model
            nm = e["sp"]["dbkey"]
            e["trans"] = nm if (nm, "") in TR else next((a for (a, b) in TR if not b and a.upper() == nm.upper()), None)

    # ---- 単成分 ----
    def species(self, r, T):
        e = self.reals[r]
        sp, model = e["sp"], e["model"]
        if model == "cea":
            return fit(TR[(e["trans"], "")]["V"], T) * 1e-7, fit(TR[(e["trans"], "")]["C"], T) * 1e-4
        if model == "custom:h2o_iapws_cea_v1":
            return h2o_iapws_cea_v1(T)
        if model == "fit":
            f = sp["fit"]
            rows = lambda k: [(r_[0], r_[1], r_[2:]) for r_ in f[k]]
            return fit(rows("V"), T) * 1e-7, fit(rows("C"), T) * 1e-4
        # kinetic
        om, Ts = omega22(T / sp["eps"])
        if sp["dipole"]:
            d = (sp["dipole"] * 1e-18)**2 / (2 * sp["eps"] * KB_CGS * (sp["sigma"] * 1e-8)**3)
            om += 0.2 * d * d / Ts
        mu = 2.6693e-6 * math.sqrt(sp["MW"] * 1e3 * T) / (sp["sigma"]**2 * om)
        return mu, mu * (cp_mass(sp, T) + 1.25 * RU / sp["MW"])

    def eta_pair(self, a, b, T, eta):
        A, B = self.reals[a], self.reals[b]
        if A["model"] == "kinetic" and B["model"] == "kinetic":
            sig = 0.5 * (A["sp"]["sigma"] + B["sp"]["sigma"])
            eps = math.sqrt(A["sp"]["eps"] * B["sp"]["eps"])
            ma, mb = A["sp"]["MW"] * 1e3, B["sp"]["MW"] * 1e3
            return 2.6693e-6 * math.sqrt(2 * ma * mb / (ma + mb) * T) / (sig**2 * omega22(T / eps)[0]), "ce"
        ta, tb = A["trans"], B["trans"]
        if ta and tb:
            for k in ((ta, tb), (tb, ta)):
                if k in TR and TR[k]["V"]:
                    return fit(TR[k]["V"], T) * 1e-7, "cea"
        # 剛体球 (cea2.f 5565–5570, i=a j=b)
        Wi, Wj = A["sp"]["MW"] * 1e3, B["sp"]["MW"] * 1e3
        ratio = math.sqrt(Wj / Wi)
        e = 5.656854 * eta[a] * math.sqrt(Wj / (Wi + Wj))
        return e / (1 + math.sqrt(ratio * eta[a] / eta[b]))**2, "rigid"

    def state(self, T, Xs):
        n = len(self.reals)
        X = [0.0] * n
        for s, ex in enumerate(self.expand):
            for r, x in ex:
                X[r] += Xs[s] * x
        eta, con = zip(*[self.species(r, T) for r in range(n)])
        W = [e["sp"]["MW"] * 1e3 for e in self.reals]
        E = {}
        eij = []
        for a in range(n):
            for b in range(a + 1, n):
                E[(a, b)] = E[(b, a)] = self.eta_pair(a, b, T, eta)[0]
                eij.append(E[(a, b)])
        mu = lam = 0.0
        for i in range(n):
            if X[i] <= 0:
                continue
            sv = sc = 0.0
            for j in range(n):
                if i == j:
                    phi = psi = 1.0
                else:
                    phi = 2.0 / (E[(i, j)] * (W[i] + W[j])) * W[j] * eta[i]
                    psi = phi * (1 + 2.41 * (W[i] - W[j]) * (W[i] - 0.142 * W[j]) / (W[i] + W[j])**2)
                sv += phi * X[j]
                sc += psi * X[j]
            mu += eta[i] * X[i] / sv
            lam += con[i] * X[i] / sc
        return {"mu": mu, "lam": lam, "Xreal": X, "mu_i": list(eta), "lam_i": list(con), "eta_ij": eij}

    def X_from_Y(self, Y):
        """輸送種の質量分率 → モル分率 (負値は 0 に切る; 和で正規化)。単成分は [1]。"""
        if len(self.mw) == 1:
            return [1.0]
        y = [max(float(v), 0.0) / m for v, m in zip(Y, self.mw)]
        t = sum(y)
        return [v / t for v in y]

    def state_Y(self, T, Y):
        return self.state(T, self.X_from_Y(Y))
