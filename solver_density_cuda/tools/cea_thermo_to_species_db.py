#!/usr/bin/env python3
"""NASA CEA `thermo.inp` から forge の `speciesDBFile` (species_db.yaml) を生成する。

  python3 cea_thermo_to_species_db.py thermo.inp --species H2 O2 H O OH HO2 H2O2 H2O N2 [--out species_db.yaml]

- NASA-9 係数 (2 温度域: 200–1000 / 1000–6000 K) を固定幅 16 文字フィールドで読む (Fortran D 指数)。
- 3 温度域以上ある種 (6000–20000 K など) は最初の 2 域だけ使う (forge は 2 域固定)。
- Lennard-Jones パラメータ (σ [Å], ε/kB [K]) は内蔵表 (Cantera h2o2.yaml / gri30.yaml の transport 値) から与える。
  表に無い種は N2 相当 (3.621, 97.53) を入れて警告する。
- 種名キーは引用符付きで書く (`"NO"`, `"N"`: PyYAML など YAML 1.1 実装で真偽値に化けるため)。
- ヘッダ 2 行目の元素欄 (5 組 × 8 文字: 元素記号 2 文字 + 原子数 6 文字) を `atoms: {H: 2.0, O: 1.0}` として出力する
  (species_meta.yaml の原子組成・元素混合分率診断が使う)。
- `--check`: 要求した種のうち設計側の転記表 `design/forge_design/gas/semiperfect.py` `SPECIES_NASA9` にある種**すべて**について
  MW と両温度域の係数を照合し、最大相対差の表を出す。係数の相対差 > `--check-tol` (既定 1e-6) か MW の相対差 > 1e-5 が
  1 つでもあれば**非ゼロ終了**する。**既知の不一致** (CEA `thermo.inp` と転記が違う; codex M7 実測): H2O の MW
  (thermo.inp 18.01528 vs 転記 0.0180153, 相対差 1.1e-6 = 既定 MW 許容 1e-5 の内側なので ok 表示) と AR の高温域 a0
  (thermo.inp 20.10538 vs 転記 0, 相対差 1 → **AR を含む `--check` は失敗するのが期待動作**)。許容するなら
  `--check-tol` を上げる (MW は `--check-mw-tol`)。
"""
import argparse, importlib.util, os, re, sys

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
        atoms = parse_atoms(hdr)
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
        db[name] = {"MW": MW, "Hf298": Hf, "intervals": intervals, "atoms": atoms}
    return db


def parse_atoms(hdr):
    """ヘッダ 2 行目 hdr[10:50] = 5 組 × (元素記号 2 文字 + 原子数 6 文字)。例 'H   2.00O   1.00' → {'H': 2.0, 'O': 1.0}。"""
    atoms = {}
    for k in range(5):
        fld = hdr[10 + 8 * k: 18 + 8 * k]
        sym = fld[:2].strip()
        cnt = fld[2:].strip()
        if not sym:
            continue
        try:
            n = float(cnt) if cnt else 0.0
        except ValueError:
            continue
        if n != 0.0:
            atoms[sym.upper() if len(sym) == 1 else sym[0].upper() + sym[1:].lower()] = n
    # CEA は 'AR', 'HE' のように 2 文字を大文字で書く種があるので、単原子希ガスは記号をそのまま大文字で保つ
    return {(k.upper() if k.upper() in ("AR", "HE", "NE", "KR", "XE") else k): v for k, v in atoms.items()}


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
        "atoms": rec.get("atoms", {}),      # 元素組成 (forge は読まない; species_meta / 元素診断用)
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
        if e.get("atoms"):
            w("  atoms: {" + ", ".join(f"{a}: {float(n)!r}" for a, n in e["atoms"].items()) + "}\n")
        w(f"  # Hf(298.15) = {e['_Hf298_J_per_mol']} J/mol (CEA thermo.inp)\n")


DESIGN_SEMIPERFECT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                   "..", "..", "design", "forge_design", "gas", "semiperfect.py"))


def load_builtin_table(path=DESIGN_SEMIPERFECT):
    """設計側の転記表 SPECIES_NASA9 (name -> dict(MW, low, high)) をパスで import する (パッケージ import に依存しない)。"""
    spec = importlib.util.spec_from_file_location("forge_semiperfect_for_check", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.SPECIES_NASA9


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-30)


def check_against_builtin(entries, tol, mw_tol, table_path=DESIGN_SEMIPERFECT):
    """要求種のうち SPECIES_NASA9 にある種すべてを照合。表を出し、許容超えがあれば False。"""
    try:
        table = load_builtin_table(table_path)
    except Exception as e:  # noqa: BLE001
        print(f"[check] SPECIES_NASA9 を読めない ({table_path}): {e}", file=sys.stderr)
        return False
    ok = True
    print(f"[check] CEA thermo.inp vs SPECIES_NASA9 ({table_path}); tol coeff {tol:.1e}, MW {mw_tol:.1e}")
    print(f"[check] {'species':8s} {'MW rel':>10s} {'low max':>10s} {'high max':>10s}  worst")
    for name, e in entries.items():
        key = name if name in table else name.upper()
        if key not in table:
            print(f"[check] {name:8s} {'-':>10s} {'-':>10s} {'-':>10s}  (not in SPECIES_NASA9; skipped)")
            continue
        ref = table[key]
        dmw = rel(e["MW"], float(ref["MW"]))
        dlo = [rel(x, float(y)) for x, y in zip(e["nasa9_low"], ref["low"])]
        dhi = [rel(x, float(y)) for x, y in zip(e["nasa9_high"], ref["high"])]
        worst = ""
        if max(dlo) > tol:
            k = max(range(9), key=lambda i: dlo[i]); worst = f"low a{k}: {e['nasa9_low'][k]!r} vs {ref['low'][k]!r}"
        if max(dhi) > tol and (not worst or max(dhi) > max(dlo)):
            k = max(range(9), key=lambda i: dhi[i]); worst = f"high a{k}: {e['nasa9_high'][k]!r} vs {ref['high'][k]!r}"
        if dmw > mw_tol:
            worst = f"MW: {e['MW']!r} vs {ref['MW']!r}" + (f"; {worst}" if worst else "")
        bad = dmw > mw_tol or max(dlo) > tol or max(dhi) > tol
        ok = ok and not bad
        print(f"[check] {name:8s} {dmw:10.2e} {max(dlo):10.2e} {max(dhi):10.2e}  {'FAIL ' if bad else 'ok   '}{worst}")
    print(f"[check] {'PASS' if ok else 'FAIL'} (known transcription differences: H2O MW, AR high-range a0)")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("thermo_inp")
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--out", default="species_db.yaml")
    ap.add_argument("--check", action="store_true",
                    help="要求種のうち SPECIES_NASA9 (design/forge_design/gas/semiperfect.py) にある種すべての MW・両温度域係数を照合し、"
                         "許容超えがあれば非ゼロ終了")
    ap.add_argument("--check-tol", type=float, default=1e-6, help="係数の相対差の許容 (既定 1e-6)")
    ap.add_argument("--check-mw-tol", type=float, default=1e-5, help="MW の相対差の許容 (既定 1e-5)")
    ap.add_argument("--check-table", default=DESIGN_SEMIPERFECT, help="照合する SPECIES_NASA9 を持つ .py (既定: design 側)")
    a = ap.parse_args()
    db = parse_thermo_inp(a.thermo_inp)
    entries = {}
    for s in a.species:
        key = s if s in db else {"AR": "Ar", "HE": "He"}.get(s, s)
        if key not in db:
            raise SystemExit(f"species '{s}' not found in {a.thermo_inp}")
        entries[s] = to_entry(key, db[key])
    with open(a.out, "w") as f:
        dump_yaml(entries, f)
    print(f"wrote {a.out}: {', '.join(entries)}")
    if a.check:
        if not check_against_builtin(entries, a.check_tol, a.check_mw_tol, a.check_table):
            sys.exit(1)


if __name__ == "__main__":
    main()
