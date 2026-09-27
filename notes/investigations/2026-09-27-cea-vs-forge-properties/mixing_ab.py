"""#5t ① 混合則の A/B (CFD 0 step)。種ごとの μ_i・λ_i は両方とも CEA trans.inp、変えるのは混合則だけ。
A = forge 現行 (Wilke φ を μ と λ で共用; cuda_forge/thermo_d.cuh thermo_mu_mix / thermo_lambda_mix)
B = CEA frozen 混合則 (cea2.f:5599-5634; ηij は trans.inp の相互作用データ、無ければ剛体球近似 cea2.f:5565)
正解 = FCEA2 の同一 T・同一組成 (only で組成凍結) の VISC と frozen CONDUCTIVITY。合格: 全 16 点で相対 0.1 % 以内。
usage (リポジトリルートから): python3 notes/investigations/2026-09-27-cea-vs-forge-properties/mixing_ab.py <FCEA2 の out があるディレクトリ>
"""
import math, re, sys, os
lines = open(".venv-cea/nasa_cea/trans.inp").read().splitlines()
num = re.compile(r"-?\d\.\d+E[ +-]\d+")
TR = {}; i = 1
while i < len(lines):
    h = lines[i]; m = re.search(r"V(\d)C(\d)", h)
    if not m: i += 1; continue
    a, b = h[0:15].strip(), h[15:30].strip(); nv, nc = int(m.group(1)), int(m.group(2)); rec = {"V": [], "C": []}
    for k in range(nv + nc):
        L = lines[i + 1 + k]; rec[L[1]].append((float(L[2:11]), float(L[11:20]), [float(x.replace("E ", "E+")) for x in num.findall(L[20:])]))
    TR[(a, b)] = rec; i += 1 + nv + nc
def fit(ivs, T):
    lo = min(ivs, key=lambda v: v[0]); hi = max(ivs, key=lambda v: v[1])
    use = next((v for v in ivs if v[0] <= T <= v[1]), lo if T < lo[0] else hi)
    A, B, C, D = use[2]; return math.exp(A * math.log(T) + B / T + C / T**2 + D)
MW = {"N2": 28.0134, "H2O": 18.01528}
def eta(a, b, T):
    if a == b: return fit(TR[(a, "")]["V"], T)
    k = (a, b) if (a, b) in TR else (b, a) if (b, a) in TR else None
    if k: return fit(TR[k]["V"], T)
    ei, ej = eta(a, a, T), eta(b, b, T); ratio = math.sqrt(MW[b] / MW[a])
    return 5.656854 * ei * math.sqrt(MW[b] / (MW[a] + MW[b])) / (1 + math.sqrt(ratio * ei / ej))**2
def con(a, T): return fit(TR[(a, "")]["C"], T)
def mixB(X, T):
    sp = list(X); mu = lam = 0.0
    for i in sp:
        sv = sc = 0.0
        for j in sp:
            if i == j: phi = psi = 1.0
            else:
                phi = 2.0 / (eta(i, j, T) * (MW[i] + MW[j])) * MW[j] * eta(i, i, T)
                psi = phi * (1 + 2.41 * (MW[i] - MW[j]) * (MW[i] - 0.142 * MW[j]) / (MW[i] + MW[j])**2)
            sv += phi * X[j]; sc += psi * X[j]
        mu += eta(i, i, T) * X[i] / sv; lam += con(i, T) * X[i] / sc
    return mu, lam      # μP, μW/(cm K)
def mixA(X, T):
    sp = list(X); mu = lam = 0.0
    wphi = lambda i, j: (1 + math.sqrt(eta(i, i, T) / eta(j, j, T)) * (MW[j] / MW[i])**0.25)**2 / math.sqrt(8 * (1 + MW[i] / MW[j]))
    for i in sp:
        d = sum(X[j] * wphi(i, j) for j in sp)
        mu += X[i] * eta(i, i, T) / d; lam += X[i] * con(i, T) / d
    return mu, lam
def read_out(path):
    t = open(path).read()
    T = [float(x) for x in re.search(r"T, K\s+([\d. ]+)", t).group(1).split()]
    V = [float(x) for x in re.search(r"VISC,MILLIPOISE\s+([\d. ]+)", t).group(1).split()]
    Cs = re.findall(r"^ CONDUCTIVITY\s+([\d. ]+)$", t, re.M)
    C = [float(x) for x in Cs[-1].split()]   # 最後の CONDUCTIVITY 行 = FROZEN 節
    return T, V, C
d = sys.argv[1]
rows = []; worst = {"A": [0, 0], "B": [0, 0]}
for xw, name in ((0.0, "mix_0"), (0.1, "mix_0p1"), (0.5, "mix_0p5"), (1.0, "mix_1")):
    X = {k: v for k, v in (("N2", 1 - xw), ("H2O", xw)) if v > 0}
    T, V, C = read_out(os.path.join(d, name + ".out"))
    for t, v, c in zip(T, V, C):
        muc, lac = v * 1e3, c * 1e3      # millipoise -> μP, mW/(cm K) -> μW/(cm K)
        (mA, lA), (mB, lB) = mixA(X, t), mixB(X, t)
        r = (xw, t, 100 * (mA / muc - 1), 100 * (lA / lac - 1), 100 * (mB / muc - 1), 100 * (lB / lac - 1))
        rows.append(r)
        for key, a, b in (("A", r[2], r[3]), ("B", r[4], r[5])):
            worst[key][0] = max(worst[key][0], abs(a)); worst[key][1] = max(worst[key][1], abs(b))
print("X_H2O   T[K]   A:μ[%]   A:λ[%]   B:μ[%]   B:λ[%]   (FCEA2 比; 合格 |差| ≤ 0.1 %)")
for r in rows: print(f"{r[0]:5.2f} {r[1]:7.1f} {r[2]:+8.3f} {r[3]:+8.3f} {r[4]:+8.3f} {r[5]:+8.3f}")
print(f"最大 |差|: A μ {worst['A'][0]:.3f} % λ {worst['A'][1]:.3f} % / B μ {worst['B'][0]:.3f} % λ {worst['B'][1]:.3f} %")
print("判定:", "A 全点合格" if max(worst['A']) <= 0.1 else "A 不合格", "/", "B 全点合格" if max(worst['B']) <= 0.1 else "B 不合格")
