#!/usr/bin/env python3
"""NASA CEA `thermo.inp` から forge の `speciesDBFile` (species_db.yaml) を生成する。

  python3 cea_thermo_to_species_db.py thermo.inp --species H2 O2 H O OH HO2 H2O2 H2O N2 [--out species_db.yaml]

- NASA-9 係数 (2 温度域: 200–1000 / 1000–6000 K) を固定幅 16 文字フィールドで読む (Fortran D 指数)。
- 3 温度域以上ある種 (6000–20000 K など) は最初の 2 域だけ使う (forge は 2 域固定)。
- Lennard-Jones パラメータ (σ [Å], ε/kB [K]) は内蔵表 (Cantera h2o2.yaml / gri30.yaml の transport 値) から与える。
  表に無い種は N2 相当 (3.621, 97.53) を入れて警告する。
- 種名キーは引用符付きで書く (`"NO"`, `"N"`: PyYAML など YAML 1.1 実装で真偽値に化けるため)。
- 生成される値は `thermo_d.cu` 内蔵 DB (N2/O2/AR/CO2/HE/H2O) と同一ソース (CEA) なので、既存種を指定しても値は一致する
  (`--check` で内蔵 N2 と突き合わせる自己検証を行う)。
"""
import argparse, re, sys

LJ = {  # (sigma [A], eps/kB [K]) — Cantera h2o2.yaml / gri30.yaml transport データ
    "H2": (2.920, 38.00), "H": (2.050, 145.00), "O": (2.750, 80.00), "O2": (3.458, 107.40),
    "OH": (2.750, 80.00), "H2O": (2.605, 572.40), "HO2": (3.458, 107.40), "H2O2": (3.458, 107.40),
    "N2": (3.621, 97.53), "N": (3.298, 71.40), "NO": (3.621, 97.53), "NO2": (3.500, 200.00),
    "N2O": (3.828, 232.40), "HNO": (3.492, 116.70), "AR": (3.330, 136.50), "CO": (3.650, 98.10),
    "CO2": (3.763, 244.00), "NH": (2.650, 80.00), "NH2": (2.650, 80.00), "NH3": (2.920, 481.00),
    "NNH": (3.798, 71.40), "HCO": (3.590, 498.00), "HE": (2.551, 10.22),
}
ALIAS = {"Ar": "AR", "He": "HE", "AR": "AR"}


def f16(s):
    """Fortran D 指数の 16 文字フィールド → float。空白は None。"""
    s = s.strip()
    if not s:
        return None
    return float(s.replace("D", "E").replace("d", "e"))


def parse_thermo_inp(path):
    """thermo.inp を {name: {"MW":..., "Hf298":..., "intervals":[(Tlo,Thi,[a0..a8]), ...]}} に読む。"""
    lines = open(path, encoding="latin-1").read().splitlines()
    db = {}
    i = 0
    # 先頭の thermo ヘッダ (コメント '!' と "thermo" 行 + 温度域行) を飛ばす
    while i < len(lines) and not lines[i].lower().startswith("thermo"):
        i += 1
    i += 2  # "thermo" 行と温度域行
    while i < len(lines):
        L = lines[i]
        if L.startswith("END") or L.startswith("end"):
            break
        if not L.strip() or L.startswith("!"):
            i += 1
            continue
        name = L[:18].strip().split()[0] if L[:18].strip() else None
        if name is None:
            i += 1
            continue
        hdr = lines[i + 1]
        n_int = int(hdr[0:2])
        MW = float(hdr[52:65])
        Hf = float(hdr[65:80])
        i += 2
        intervals = []
        if n_int == 0:  # 反応物専用 (係数無し): 1 行だけ
            i += 1
        for _ in range(n_int):
            rng = lines[i]
            Tlo, Thi = float(rng[0:11]), float(rng[11:22])
            c1 = lines[i + 1]
            c2 = lines[i + 2]
            a = [f16(c1[k * 16:(k + 1) * 16]) for k in range(5)]
            a += [f16(c2[0:16]), f16(c2[16:32])]
            # c2[32:48] は空 (8 番目の係数は未使用), b1,b2 が c2[48:64], c2[64:80]
            b1, b2 = f16(c2[48:64]), f16(c2[64:80])
            a = [v if v is not None else 0.0 for v in a] + [b1 or 0.0, b2 or 0.0]
            intervals.append((Tlo, Thi, a))
            i += 3
        db[name] = {"MW": MW, "Hf298": Hf, "intervals": intervals}
    return db


def to_entry(name, rec):
    iv = rec["intervals"]
    if len(iv) < 2:
        raise SystemExit(f"{name}: forge は 2 温度域が必要 (thermo.inp は {len(iv)} 域)")
    lo, hi = iv[0], iv[1]
    key = ALIAS.get(name, name.upper() if name.upper() in LJ else name)
    if key in LJ:
        sig, eps = LJ[key]
    else:
        sig, eps = 3.621, 97.53
        print(f"[warn] {name}: LJ パラメータが内蔵表に無い → N2 相当 (3.621, 97.53) を仮置き", file=sys.stderr)
    return {
        "MW": rec["MW"] * 1e-3,           # g/mol → kg/mol
        "LJ_sigma": sig, "LJ_eps_kB": eps,
        "Tlo": lo[0], "Tmid": lo[1], "Thi": hi[1],
        "nasa9_low": lo[2], "nasa9_high": hi[2],
        "_Hf298_J_per_mol": rec["Hf298"],  # 参考 (forge は読まない)
    }


def dump_yaml(entries, out):
    w = out.write
    for name, e in entries.items():
        w(f"\"{name}\":\n")   # NO / N / Y は YAML 1.1 で真偽値に化けるので必ず引用符
        for k in ("MW", "LJ_sigma", "LJ_eps_kB", "Tlo", "Tmid", "Thi"):
            w(f"  {k}: {e[k]!r}\n")
        for k in ("nasa9_low", "nasa9_high"):
            w(f"  {k}:\n")
            for v in e[k]:
                w(f"  - {v!r}\n")
        w(f"  # Hf(298.15) = {e['_Hf298_J_per_mol']} J/mol (CEA thermo.inp)\n")


BUILTIN_N2_LOW = [2.210371497e+04, -3.818461820e+02, 6.082738360e+00, -8.530914410e-03, 1.384646189e-05,
                  -9.625793620e-09, 2.519705809e-12, 7.108460860e+02, -1.076003744e+01]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("thermo_inp")
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--out", default="species_db.yaml")
    ap.add_argument("--check", action="store_true", help="内蔵 N2 係数と照合する")
    a = ap.parse_args()
    db = parse_thermo_inp(a.thermo_inp)
    if a.check:
        n2 = db["N2"]["intervals"][0][2]
        err = max(abs(x - y) / max(abs(y), 1e-30) for x, y in zip(n2, BUILTIN_N2_LOW))
        print(f"[check] N2 low-range coeffs vs thermo_d.cu builtin: max rel diff = {err:.2e}")
    entries = {}
    for s in a.species:
        key = s if s in db else {"AR": "Ar", "HE": "He"}.get(s, s)
        if key not in db:
            raise SystemExit(f"species '{s}' not found in {a.thermo_inp}")
        entries[s] = to_entry(key, db[key])
    with open(a.out, "w") as f:
        dump_yaml(entries, f)
    print(f"wrote {a.out}: {', '.join(entries)}")


if __name__ == "__main__":
    main()
