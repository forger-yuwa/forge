#!/usr/bin/env python3
"""有限速度化学の要否スクリーニング: CEA rocket (equilibrium / frozen nfz=1 / frozen nfz=2) の出口 T, M, 組成を比較する。
実行: cd notes/investigations/scripts/cea_kinetics_screen && python3 run_screen.py
(thermo.lib/trans.lib は case/44/cea からコピー済。FCEA2 は .venv-cea/nasa_cea/FCEA2)"""
import subprocess, re, sys, pathlib, itertools
HERE = pathlib.Path(__file__).parent
FCEA2 = HERE / "../../../../.venv-cea/nasa_cea/FCEA2"
SUPAR = [15.07, 53.2]          # M≈4.2 (case/44) と M≈6 (case/45 級)
P_BAR = 11.39

def problem(case, mode, reac_lines):
    if mode == "eq":   head = f"problem case={case} rocket equilibrium p,bar={P_BAR} supar={','.join(map(str,SUPAR))}"
    elif mode == "fz1": head = f"problem case={case} rocket frozen nfz=1 p,bar={P_BAR} supar={','.join(map(str,SUPAR))}"
    else:               head = f"problem case={case} rocket frozen nfz=2 p,bar={P_BAR} supar={','.join(map(str,SUPAR))}"
    return head + "\nreac\n" + "\n".join(reac_lines) + "\noutput siunits\nend\n"

# (A) case/44 va3 組成を温度だけ上げて与える (解離平衡の効果だけを見る)
VA3 = [("H2O",0.0609135),("N2",0.664860),("O2",0.216072),("Ar",0.00797588),("CO2",0.0490034)]
def reac_va3(Tt): return [f"  name={n} moles={m} t,k={Tt}" for n,m in VA3]
# (B) H2 燃焼加熱器 (air + H2 + O2 補充で残留 O2 = 21 mol%): 反応物 300 K、燃焼温度は CEA が決める
def reac_heater(h):
    r = 0.21*(h+0.7902)/0.79           # 残留 O2 [mol / mol air]
    o2add = r + 0.5*h - 0.2095
    return [f"  oxid=Air moles=1.0 t,k=300", f"  fuel=H2 moles={h:.4f} t,k=300", f"  oxid=O2 moles={o2add:.4f} t,k=300"]

cases = []
for Tt in (1161, 1600, 2000, 2500, 3000): cases.append((f"A{Tt}", reac_va3(Tt)))
for h in (0.10, 0.20, 0.30, 0.40, 0.55, 0.70): cases.append((f"B{int(h*100):02d}", reac_heater(h)))

inp = ""
for name, reac in cases:
    for mode in ("eq","fz1","fz2"):
        inp += problem(f"{name}_{mode}", mode, reac)
(HERE/"screen.inp").write_text(inp)
subprocess.run([str(FCEA2)], input="screen\n", text=True, cwd=HERE, stdout=subprocess.DEVNULL)
txt = (HERE/"screen.out").read_text().splitlines()

def num(t):
    m = re.fullmatch(r"([0-9.]+)([+-]\d+)", t)
    return float(m.group(1))*10**int(m.group(2)) if m else float(t)

def parse(case):
    """CASE ブロックから T, M, MOLE FRACTIONS の列 (CHAMBER, THROAT, EXIT...) を返す"""
    res = {}
    for i,L in enumerate(txt):
        if L.strip().startswith("CASE =") and L.split("=")[1].strip() == case:
            j = i
            while j < len(txt) and not (txt[j].lstrip().startswith("Pinf/P") or txt[j].lstrip().startswith("Pin")): j += 1
            k = j; T = M = None; X = {}
            while k < len(txt) and not txt[k].strip().startswith("CASE ="):
                s = txt[k].lstrip()
                if s.startswith("T, K"): T = [num(t) for t in s.split()[2:]]
                if s.startswith("MACH NUMBER"): M = [num(t) for t in s.split()[2:]]
                if s.startswith("MOLE FRACTIONS") or s.startswith("MASS FRACTIONS"):
                    k += 2
                    while k < len(txt) and txt[k].strip():
                        toks = txt[k].split()
                        # 行内に複数種 (name v v v ...) が並ぶ場合は数値の個数で区切る
                        t = 0
                        while t < len(toks):
                            nm = toks[t].lstrip("*"); t += 1; vals = []
                            while t < len(toks) and re.fullmatch(r"[0-9.]+([+-]\d+)?", toks[t]): vals.append(num(toks[t])); t += 1
                            X[nm] = vals
                        k += 1
                    break
                k += 1
            return T, M, X
    return None

lines = ["| case | mode | T_c [K] | T_exit(M4.2) | M_exit(M4.2) | T_exit(M6) | M_exit(M6) | X_OH,c | X_H,c | X_O,c | X_NO,c | X_H2O,c | X_OH,exit(M4.2) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for name, _ in cases:
    for mode in ("eq","fz1","fz2"):
        r = parse(f"{name}_{mode}")
        if r is None: lines.append(f"| {name} | {mode} | (parse fail) |"); continue
        T, M, X = r
        def x(sp, col=0):
            v = X.get(sp); 
            if not v: return 0.0
            return v[col] if col < len(v) else v[-1]
        def g(v, i): return v[i] if i < len(v) else float('nan')   # CEA が途中で停止した列は NaN
        lines.append(f"| {name} | {mode} | {T[0]:.0f} | {g(T,2):.1f} | {g(M,2):.3f} | {g(T,3):.1f} | {g(M,3):.3f} | {x('OH'):.2e} | {x('H'):.2e} | {x('O'):.2e} | {x('NO'):.2e} | {x('H2O'):.4f} | {x('OH',2):.2e} |")
(HERE/"screen_table.md").write_text("\n".join(lines)+"\n")
print("\n".join(lines))
