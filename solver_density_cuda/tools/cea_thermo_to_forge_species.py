#!/usr/bin/env python3
"""CEA `thermo.inp` → forge 共通データ形式 (forge_species_data) の種エントリを作る生成器。

plans/active/thermophysics-solver-owned-species-db.md §4.9・§5.1 #13-0 (監査)・#13-2 (生成と内蔵拡大)。仕様 methods/thermophysics.md。

  python3 solver_density_cuda/tools/cea_thermo_to_forge_species.py --thermo PATH/thermo.inp --audit [--audit-out FILE.md]
      [--sern-root /home/sano/work/forge-sern-design]
  python3 solver_density_cuda/tools/cea_thermo_to_forge_species.py --thermo PATH/thermo.inp --write    # 生成ブロックを書き換える
  python3 solver_density_cuda/tools/cea_thermo_to_forge_species.py --thermo PATH/thermo.inp --check    # 往復のビット一致を検査

- 段 0 (#13-0) の `--audit` は監査表 (Markdown) を出すだけ。共通データは書き換えない。
- 段 2 (#13-2) の `--write` は共通データ (`data/species/forge_species_v1.yaml`) の**生成ブロック** (`GENERATED_BEGIN`〜`GENERATED_END`
  の行の間) だけを置き換える。対象は輸送データ (trans.inp 由来 66 種) と thermo.inp に同名 (大小文字区別) がある種のうち、
  生成ブロックの外 (手保守のエントリ: 既存 7 種・`H2O(L)`・CO/H2/OH/H/NO/O) に無いもの。`e-` は除外 (区切り 298.15 K; §4.9)、
  CEA の `Air` (reactants 節) は輸送データに無いので対象外。区間は CEA の区間数・区切りのまま (1〜3 区間)。
  段 3 (#13-3) から `--write` は手保守のエントリも同期する: thermo.inp に同名の気相がある種 (N2/O2/CO2/H2O/Ar/He/
  CO/H2/OH/H/NO/O) の `MW:` 行と `intervals:` 節を thermo.inp そのもの (区間数・区切り・係数・MW) に置き換え、凝縮相
  (`H2O(L)`) は係数を変えず `MW:` だけを気液ペアの気相種と同じ値にする。手保守の欄 (aliases・LJ・legacy_builtin・source・
  deviations・pair_of・extension・コメント) は書き換えない。CEA に同名の無い `AIR` (擬似種) はそのまま。
- LJ は CEA に無い。共通データは出典別の集合 `LJ_sets` を持つ (plan §4.10, #14)。`legacy_v1` (#14 前の値の凍結) は
  `cea_thermo_to_species_db.py` の `LJ` 表 (#13-2 で埋めた値) をそのまま使い、種ごとの実際の出典を `LEGACY_V1_SOURCE` に書く
  (2026-10-01 の監査で表の値が GRI-Mech 3.0 と Svehla 1962 の混在と分かった)。表に無い種は `LJ_sets: {}`
  (輸送データなし; kinetic 輸送・LJ 混合平均拡散に使うとソルバが起動時に拒否する)。
  #14-L1 から `--write` は全気相種 (手保守・生成ブロックとも) の `gri30`・`svehla1962` 集合を書く:
    gri30       Cantera 同梱 `gri30.yaml` (GRI-Mech 3.0 の transport; ck2yaml が gri30_tran.dat から変換) の diameter・well-depth。
                名前の対応は完全一致 → 大小文字無視で一意 → 明示の表 (`GRI30_NAME`: AR→Ar, C2H2→C2H2,acetylene)。
                YAML 1.1 の真偽値 (NO) を避けるため BaseLoader で読む。双極子 (dipole) は種レベルの `dipole` に書く
                (手保守の種は手の値と一致することだけを検査する)。
    svehla1962  Svehla 1962 (NASA TR R-132) Table I(a) の 2 本独立の転記 CSV (notes/investigations/2026-10-01-svehla1962-lj/)。
                生成時に 2 本の全値の一致と、全分子の頁参照 (svehla1962_table1a_pages.csv) を検査する。希ガスは粘性フィット行
                (転記 CSV の行)。名前の対応は gri30 と同じ規則 (明示の表 `SVEHLA_NAME`: Air→AIR, C2H2→C2H2,acetylene)。
  手保守エントリでは `LJ_sets:` の中の `gri30:`・`svehla1962:` の行と、provenance の `lj_sets.gri30`・`lj_sets.svehla1962` の行だけを
  書き換える (`legacy_v1:` と `dipole:` は手保守)。
- `--check` は (a) 輸送データと thermo.inp に同名の全種 (`e-` を含む) を `to_forge_entry` → `dump_entries_yaml` → YAML 読み戻しで
  区間・全係数・MW が thermo.inp のパース値とビット一致、(b) 共通データの生成ブロックが今の生成結果と文字列一致、
  (c) 生成ブロックの各エントリの区間・係数・MW が thermo.inp とビット一致、(d) 手保守のエントリのうち thermo.inp に同名の
  気相がある種の区間・係数・MW が thermo.inp とビット一致、凝縮相の MW が気液ペアの気相種とビット一致、手保守テキストが
  同期結果と文字列一致、を見て、不一致があれば非ゼロ終了する。
- 入力の読み方 (thermo.inp の NASA-9 固定桁; McBride, Zehe, Gordon NASA/TP-2002-211556 App. A):
    記録 1 行目   種名 (1–24 桁; 空白を含まない) + 注記
    記録 2 行目   区間数 [0:2]・参照日付コード [3:9]・元素 5 組 × (記号 2 + 個数 6) [10:50]・相 [50:52] (0 = 気相)・
                  MW [g/mol] [52:65]・Hf(298.15) [J/mol] (気相) / H(298.15) (凝縮相) [65:80]
    区間ごとに 3 行: Tlo [0:11]・Thi [11:22]・…/ a0..a4 (16 桁 × 5)/ a5, a6, (空), b1, b2 (16 桁 × 5; Fortran D 指数)
    区間数 0 (反応物専用) は 3 行目に温度と H(298.15) の 1 行だけ。
  forge の係数並び a0..a8 = CEA の a1..a7, b1, b2 (既存 `cea_thermo_to_species_db.py` と同じ)。
- 数値の変換: 係数は 16 桁欄を D→E にして Python float (最近接丸め)。MW は欄の文字列に `e-3` を付けて float にする
  (g/mol → kg/mol を 10 進のまま行い、`x * 1e-3` の二重丸めを避ける)。
- 名前は大小文字を区別する (`CO` と `Co` は別種; plan §4.9)。
"""
import argparse
import hashlib
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.normpath(os.path.join(HERE, ".."))
REPO = os.path.normpath(os.path.join(SOLVER, ".."))
DEFAULT_THERMO = os.path.join(REPO, ".venv-cea", "nasa_cea", "thermo.inp")
DEFAULT_SPECIES = os.path.join(SOLVER, "data", "species", "forge_species_v1.yaml")
DEFAULT_TRANSPORT = os.path.join(SOLVER, "data", "species", "forge_transport_v1.yaml")
SVEHLA_DIR = os.path.join(REPO, "notes", "investigations", "2026-10-01-svehla1962-lj")
DEFAULT_SVEHLA = (os.path.join(SVEHLA_DIR, "svehla1962_table1a_transcriptA.csv"),
                  os.path.join(SVEHLA_DIR, "svehla1962_table1a_transcriptB.csv"))
DEFAULT_SVEHLA_PAGES = os.path.join(SVEHLA_DIR, "svehla1962_table1a_pages.csv")


def default_gri30():
    """Cantera 同梱 gri30.yaml (このワークツリーの .venv-chem、無ければ隣の forge ワークツリーの .venv-chem)。"""
    import glob
    for root in (REPO, os.path.join(os.path.dirname(REPO), "forge")):
        hit = sorted(glob.glob(os.path.join(root, ".venv-chem", "lib", "python3*", "site-packages", "cantera", "data", "gri30.yaml")))
        if hit:
            return hit[-1]
    return None
RU = 8.314462618            # J/(mol K); cuda_forge/thermo_d.cuh THERMO_RU と同じ
STD_BOUNDS = (200.0, 1000.0, 6000.0, 20000.0)


# ---------------------------------------------------------------------------------------------------------------
# パーサ
# ---------------------------------------------------------------------------------------------------------------
def _f16(s):
    """Fortran D 指数の 16 桁欄 → float。空欄は 0.0。"""
    t = s.strip()
    if not t:
        return 0.0
    return float(t.replace("D", "E").replace("d", "e"))


def _atoms(hdr):
    """元素 5 組 (記号 2 桁 + 個数 6 桁)。記号は thermo.inp の表記のまま (例 'AR', 'CO' = コバルト)。個数 0 は除く。"""
    out = {}
    for k in range(5):
        fld = hdr[10 + 8 * k:18 + 8 * k]
        sym, cnt = fld[:2].strip(), fld[2:].strip()
        if not sym:
            continue
        try:
            n = float(cnt) if cnt else 0.0
        except ValueError:
            continue
        if n != 0.0:
            out[sym] = n
    return out


def parse_thermo_inp(path):
    """thermo.inp を記録の列 (ファイル順) に読む。各記録は dict:
    name, section ('products'|'reactants'), phase (int; 0 = 気相), MW_str (g/mol の欄文字列), MW_gmol, Hf, atoms,
    date, comment, intervals [(Tlo, Thi, [a0..a8]), ...]。"""
    lines = open(path, encoding="latin-1").read().splitlines()
    i = 0
    while i < len(lines) and not lines[i].lower().startswith("thermo"):
        i += 1
    i += 2                                   # "thermo" 行と全体の温度域行
    recs = []
    section = "products"
    while i < len(lines):
        L = lines[i]
        if L.upper().startswith("END PRODUCTS"):
            section = "reactants"
            i += 1
            continue
        if L.upper().startswith("END REACTANTS"):
            break
        if not L.strip() or L.startswith("!"):
            i += 1
            continue
        name = L[:24].split()[0]
        comment = L[18:].strip()
        hdr = lines[i + 1]
        n_int = int(hdr[0:2])
        mw_str = hdr[52:65].strip()
        rec = {"name": name, "section": section, "phase": int(hdr[50:52]), "MW_str": mw_str,
               "MW_gmol": float(mw_str), "Hf": float(hdr[65:80]), "atoms": _atoms(hdr), "date": hdr[3:9].strip(),
               "comment": comment, "intervals": []}
        i += 2
        if n_int == 0:
            i += 1
        for _ in range(n_int):
            rng, c1, c2 = lines[i], lines[i + 1], lines[i + 2]
            Tlo, Thi = float(rng[0:11]), float(rng[11:22])
            a = [_f16(c1[k * 16:(k + 1) * 16]) for k in range(5)]
            a += [_f16(c2[0:16]), _f16(c2[16:32]), _f16(c2[48:64]), _f16(c2[64:80])]
            rec["intervals"].append((Tlo, Thi, a))
            i += 3
        recs.append(rec)
    return recs


def sha256_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------
# forge 共通データ形式のエントリ
# ---------------------------------------------------------------------------------------------------------------
def mw_kg(rec):
    """g/mol の欄文字列を 10 進のまま kg/mol に (x * 1e-3 の二重丸めを避ける)。"""
    return float(rec["MW_str"] + "e-3")


def to_forge_entry(rec, lj=None, dipole=None):
    """1 記録 → forge_species_data のエントリ dict (aliases・pair_of などの手保守欄は付けない)。区間は CEA のまま。
    lj: LJ_sets {集合名: {sigma, eps_kB, source}} (空 / None は輸送データなし)。"""
    return {
        "id": rec["name"],
        "aliases": [],
        "phase": "gas" if rec["phase"] == 0 else "condensed",
        "MW": mw_kg(rec),
        "intervals": [{"Tlo": lo, "Thi": hi, "coeffs": list(a)} for lo, hi, a in rec["intervals"]],
        "LJ_sets": dict(lj or {}),
        "dipole": dipole,
        "atoms": dict(rec["atoms"]) or None,
        "source": {"thermo": f"CEA thermo.inp ({rec['date']}; {rec['comment']})"},
    }


def _num(v):
    """double を YAML の数値として書く。repr (往復が一致する最短表記) に、指数表記で小数点が無いときだけ '.0' を足す
    (PyYAML は YAML 1.1 の規則で '1e-05' を文字列に読むため。'1.0e-05' は同じ double)。"""
    r = repr(float(v))
    if "e" in r and "." not in r:
        m, e = r.split("e")
        r = f"{m}.0e{e}"
    return r


def _yq(s):
    """YAML の二重引用符文字列 (\\ と " だけをエスケープ)。"""
    return '"' + str(s).replace("\\", "\\\\").replace('"', '\\"') + '"'


def dump_entries_yaml(entries, out):
    """エントリ列を forge_species_data の `species:` 以下の形で書く (数値は `_num` = double の往復が一致する最短表記)。"""
    w = out.write
    for e in entries:
        w(f'  - id: {_yq(e["id"])}\n')
        w("    aliases: [" + ", ".join(_yq(a) for a in e["aliases"]) + "]\n")
        w(f"    phase: {e['phase']}\n")
        w(f"    MW: {_num(e['MW'])}\n")
        w("    intervals:\n")
        for iv in e["intervals"]:
            w(f"      - Tlo: {_num(iv['Tlo'])}\n        Thi: {_num(iv['Thi'])}\n")
            w("        coeffs: [" + ", ".join(_num(v) for v in iv["coeffs"]) + "]\n")
        sets = e.get("LJ_sets") or {}
        if sets:
            w("    LJ_sets:\n")
            for k, lj in sets.items():
                w(f"      {k}: {{sigma: {_num(lj['sigma'])}, eps_kB: {_num(lj['eps_kB'])}, source: {_yq(lj['source'])}}}\n")
        else:
            w("    LJ_sets: {}\n")
        if e.get("dipole"):
            w(_dipole_line(e["dipole"]["value"], e["dipole"]["source"]))
        if e["atoms"]:
            w("    atoms: {" + ", ".join(f"{k}: {_num(v)}" for k, v in e["atoms"].items()) + "}\n")
        else:
            w("    atoms: null\n")
        w(f'    source: {{thermo: {_yq(e["source"]["thermo"])}}}\n')


# ---------------------------------------------------------------------------------------------------------------
# 段 2 (#13-2): 共通データの生成ブロック
# ---------------------------------------------------------------------------------------------------------------
GENERATED_BEGIN = "  # ---- BEGIN generated by tools/cea_thermo_to_forge_species.py --write"
GENERATED_END = "  # ---- END generated by tools/cea_thermo_to_forge_species.py"
EXCLUDED = ("e-",)          # 区切り 298.15 K (非標準) で内蔵から除外 (plan §4.9)
_GRI = "GRI-Mech 3.0 transport (gri30_tran.dat; Cantera 同梱 gri30.yaml と同値)"
_SV = "Svehla 1962 (NASA TR R-132) Table I(a)"
# legacy_v1 (#14 前の値の凍結) の種ごとの実際の出典 (2026-10-01 監査: gri30.yaml と Svehla 1962 の転記
# notes/investigations/2026-10-01-svehla1962-lj/ と再照合)。値は cea_thermo_to_species_db.py の LJ 表のまま。
# 表から値を取る種でここに無いものは --write/--check で止める (出典を書かずに値を入れない)。
LEGACY_V1_SOURCE = {
    "N": f"{_SV} p.23 (gri30.yaml も同値)。#13-2 で cea_thermo_to_species_db.py の LJ 表から",
    "NH3": f"{_GRI} (Svehla 1962 は 2.900/558.3)。#13-2 で cea_thermo_to_species_db.py の LJ 表から",
    "NO2": f"{_GRI} (Svehla 1962 Table I(a) に無い)。#13-2 で cea_thermo_to_species_db.py の LJ 表から",
    "N2O": f"{_SV} p.23 (gri30.yaml も同値)。#13-2 で cea_thermo_to_species_db.py の LJ 表から",
}


def _lj_table():
    """LJ の出典 (1 つに固定): cea_thermo_to_species_db.py の LJ 表と別名表をそのまま使う (写さない)。"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_cea_thermo_to_species_db", os.path.join(HERE, "cea_thermo_to_species_db.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.LJ, m.ALIAS


def lj_for(name, table, alias):
    """種名 → LJ_sets ({"legacy_v1": {sigma, eps_kB, source}} または {})。引き方は cea_thermo_to_species_db.py の to_entry と同じ。"""
    key = alias.get(name, name.upper() if name.upper() in table else name)
    if key not in table:
        return {}
    if name not in LEGACY_V1_SOURCE:
        raise SystemExit(f"{name}: LJ 表に値があるが LEGACY_V1_SOURCE に出典が無い (出典を書かずに値を入れない)")
    sig, eps = table[key]
    return {"legacy_v1": {"sigma": float(sig), "eps_kB": float(eps), "source": LEGACY_V1_SOURCE[name]}}


# ---------------------------------------------------------------------------------------------------------------
# #14-L1: LJ の集合 gri30 / svehla1962 (plan §4.10)
# ---------------------------------------------------------------------------------------------------------------
GRI30_NAME = {"AR": "Ar", "C2H2": "C2H2,acetylene"}       # gri30.yaml の名前 → 共通データの id (完全一致・大小文字無視で決まらないもの)
SVEHLA_NAME = {"Air": "AIR", "C2H2": "C2H2,acetylene"}    # Svehla の分子名 → 共通データの id (同上)
LJ_SET_ORDER = ("legacy_v1", "gri30", "svehla1962")        # LJ_sets の書き出し順
# 集合のポテンシャル形 (plan §4.10 「双極子の適用規則」, #14-L1b)。共通データのトップレベル lj_sets (手保守) がこれと一致することを --check が見る
LJ_SET_POTENTIAL = {"legacy_v1": "stockmayer", "gri30": "stockmayer", "svehla1962": "lj12-6"}


def map_name(name, ids, explicit):
    """外部の名前 → 共通データの id (明示の表 → 完全一致 → 大小文字無視で一意)。無ければ None。"""
    if name in explicit:
        return explicit[name]
    if name in ids:
        return name
    hit = [i for i in ids if i.upper() == name.upper()]
    if len(hit) > 1:
        raise SystemExit(f"{name}: 大小文字無視で共通データの id が 1 つに決まらない {hit}")
    return hit[0] if hit else None


def load_svehla(paths=DEFAULT_SVEHLA, pages_path=DEFAULT_SVEHLA_PAGES):
    """転記 CSV 2 本を読み、全分子・全値 (文字列) の一致を検査して {分子: (sigma, eps_kB, 頁)} と来歴を返す。"""
    import csv
    rows = []
    for p in paths:
        with open(p, newline="", encoding="utf-8") as f:
            rows.append({r["molecule"]: (r["sigma_A"], r["eps_over_k_K"]) for r in csv.DictReader(f)})
    a, b = rows
    if a != b:
        diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
        raise SystemExit(f"Svehla 転記 2 本が一致しない ({len(diff)} 分子: {diff[:10]})")
    with open(pages_path, newline="", encoding="utf-8") as f:
        pages = {r["molecule"]: int(r["page"]) for r in csv.DictReader(f)}
    if set(pages) != set(a):
        raise SystemExit(f"Svehla の頁参照が転記と一致しない: {sorted(set(pages) ^ set(a))}")
    out = {m: (float(s), float(e), pages[m]) for m, (s, e) in a.items()}
    prov = (f"Svehla 1962 (NASA TR R-132, NTRS 19630012982) Table I(a), 2 本独立の転記 {len(a)} 分子が全値一致 "
            f"(A sha256 {sha256_file(paths[0])[:16]}, B {sha256_file(paths[1])[:16]}, 頁 {sha256_file(pages_path)[:16]}; "
            f"notes/investigations/2026-10-01-svehla1962-lj/)。希ガスは粘性フィット行")
    return out, prov


def load_gri30(path):
    """Cantera gri30.yaml の transport → {名前: (diameter, well-depth, dipole または None)} と来歴。BaseLoader (NO を真偽値にしない)。"""
    import yaml
    with open(path, encoding="utf-8") as f:
        raw = yaml.load(f, Loader=yaml.BaseLoader)
    out = {}
    for s in raw["species"]:
        tr = s.get("transport") or {}
        if "diameter" not in tr or "well-depth" not in tr:
            continue
        out[str(s["name"])] = (float(tr["diameter"]), float(tr["well-depth"]), float(tr["dipole"]) if "dipole" in tr else None)
    desc = str(raw.get("description", "")).splitlines()[0].strip()
    prov = (f"{desc} の transport (Cantera 同梱 gri30.yaml; ck2yaml (Cantera {raw.get('cantera-version', '?')}) が "
            f"{', '.join(raw.get('input-files') or [])} から変換; sha256 {sha256_file(path)[:16]})")
    return out, prov


def lj_sets_for_ids(ids, gri30_path, svehla_paths=DEFAULT_SVEHLA, pages_path=DEFAULT_SVEHLA_PAGES):
    """共通データの id 列 → ({id: {"gri30": {...}, "svehla1962": {...}}}, {id: gri30 の dipole}, 来歴 {集合: 文字列}, 照合の情報)。"""
    sv, sv_prov = load_svehla(svehla_paths, pages_path)
    gr, gr_prov = load_gri30(gri30_path)
    sets, dip = {i: {} for i in ids}, {}
    info = {"svehla_unmapped": [], "gri30_unmapped": []}
    for m, (s, e, p) in sv.items():
        i = map_name(m, ids, SVEHLA_NAME)
        if i is None:
            info["svehla_unmapped"].append(m)
            continue
        note = f"Svehla 1962 Table I(a) p.{p}" + ("" if m == i else f" ({m})")
        sets[i]["svehla1962"] = {"sigma": s, "eps_kB": e, "source": note}
    for n, (s, e, d) in gr.items():
        i = map_name(n, ids, GRI30_NAME)
        if i is None:
            info["gri30_unmapped"].append(n)
            continue
        sets[i]["gri30"] = {"sigma": s, "eps_kB": e, "source": "GRI-Mech 3.0 (Cantera gri30.yaml)" + ("" if n == i else f" ({n})")}
        if d is not None:
            dip[i] = d
    return sets, dip, {"gri30": gr_prov, "svehla1962": sv_prov}, info


def _lj_line(k, lj):
    return f"      {k}: {{sigma: {_num(lj['sigma'])}, eps_kB: {_num(lj['eps_kB'])}, source: {_yq(lj['source'])}}}\n"


def _dipole_line(value, source):
    return f"    dipole: {{value: {_num(value)}, source: {_yq(source)}}}\n"


DIPOLE_SOURCE_GRI30 = "GRI-Mech 3.0 (Cantera gri30.yaml)。physProp.transport の kinetic (Brokaw の極性補正) だけが読む (#5t2)"


def sync_lj_text(pre, sets, prov):
    """手保守のテキストの各エントリの `LJ_sets:` の gri30・svehla1962 行と、provenance の lj_sets の gri30・svehla1962 行を
    生成結果に置き換える (legacy_v1 の行・dipole・他の欄は一字も変えない)。戻り値 (新テキスト, 書き換えたエントリ)。"""
    lines, blocks = _manual_blocks(pre)
    done = []
    for sid, k0, k1 in reversed(blocks):
        blk = lines[k0:k1]
        ls = [n for n, L in enumerate(blk) if L.startswith("    LJ_sets:")]
        if len(ls) != 1:
            raise SystemExit(f"手保守エントリ {sid}: LJ_sets の行が 1 つでない")
        n0 = ls[0]
        n1 = n0 + 1
        while n1 < len(blk) and blk[n1].startswith("      "):
            n1 += 1
        keep = [L for L in blk[n0 + 1:n1] if L.startswith("      legacy_v1:")]
        others = [L for L in blk[n0 + 1:n1] if not L.startswith(("      legacy_v1:", "      gri30:", "      svehla1962:"))]
        if others:
            raise SystemExit(f"手保守エントリ {sid}: LJ_sets に未知の行 {others}")
        new = keep + [_lj_line(k, sets.get(sid, {})[k]) for k in LJ_SET_ORDER[1:] if k in sets.get(sid, {})]
        blk = blk[:n0] + (["    LJ_sets:\n"] + new if new else ["    LJ_sets: {}\n"]) + blk[n1:]
        done.append(sid)
        lines[k0:k1] = blk
    text = "".join(lines)
    # provenance の lj_sets (legacy_v1 は手保守)
    pl = text.splitlines(keepends=True)
    k = [n for n, L in enumerate(pl) if L.startswith("  lj_sets:")]
    if len(k) != 1:
        raise SystemExit("provenance に lj_sets の行が 1 つでない")
    n1 = k[0] + 1
    while n1 < len(pl) and pl[n1].startswith("    "):
        n1 += 1
    keep = [L for L in pl[k[0] + 1:n1] if L.startswith("    legacy_v1:")]
    pl = pl[:k[0] + 1] + keep + [f"    {s}: {_yq(prov[s])}\n" for s in LJ_SET_ORDER[1:]] + pl[n1:]
    return "".join(pl), list(reversed(done))


def _load_yaml_text(text):
    import yaml
    return yaml.safe_load(text) or {}


def target_names(recs, transport_path):
    """輸送データの種のうち thermo.inp に同名 (大小文字区別) がある種 (輸送データの順) と、名前 → 記録 (最初の記録)。"""
    by = {}
    for r in recs:
        by.setdefault(r["name"], r)
    tr_ids = [str(s["id"]) for s in _load_yaml_text(open(transport_path, encoding="utf-8").read())["species"]]
    return [s for s in tr_ids if s in by], by


def split_generated(text):
    """共通データのテキスト → (前, 生成ブロックの中身, 後)。マーカ行は前・後に含める。マーカが無い・重複は例外。"""
    lines = text.splitlines(keepends=True)
    b = [k for k, L in enumerate(lines) if L.startswith(GENERATED_BEGIN)]
    e = [k for k, L in enumerate(lines) if L.startswith(GENERATED_END)]
    if len(b) != 1 or len(e) != 1 or not b[0] < e[0]:
        raise SystemExit(f"共通データに生成ブロックのマーカ (BEGIN {len(b)} 個, END {len(e)} 個) が 1 組ない")
    return "".join(lines[:b[0] + 1]), "".join(lines[b[0] + 1:e[0]]), "".join(lines[e[0]:])


def lj_context(thermo_path, species_path, transport_path, gri30_path):
    """全気相種 (手保守 + 生成ブロック) の id と、その LJ 集合 (gri30・svehla1962)・gri30 の双極子・来歴・照合の情報。"""
    recs = parse_thermo_inp(thermo_path)
    names, _ = target_names(recs, transport_path)
    pre, _, post = split_generated(open(species_path, encoding="utf-8").read())
    man = _load_yaml_text(pre + post).get("species") or []
    manual = {str(e["id"]) for e in man}
    ids = [str(e["id"]) for e in man if e.get("phase") == "gas"] + [s for s in names if s not in manual and s not in EXCLUDED]
    if gri30_path is None or not os.path.exists(gri30_path):
        raise SystemExit(f"gri30.yaml が無い: {gri30_path} (--gri30 で Cantera 同梱の gri30.yaml を指定)")
    sets, dip, prov, info = lj_sets_for_ids(ids, gri30_path)
    return ids, sets, dip, prov, info


def generated_entries(thermo_path, species_path, transport_path, gri30_path):
    """生成ブロックに入れるエントリ列 (手保守のエントリに無い種だけ) とその YAML テキスト、対象名、名前 → 記録。
    LJ_sets は legacy_v1 (LJ 表) + gri30 + svehla1962、gri30 に双極子があれば種レベルの dipole (#14-L1)。"""
    import io
    recs = parse_thermo_inp(thermo_path)
    names, by = target_names(recs, transport_path)
    pre, _, post = split_generated(open(species_path, encoding="utf-8").read())
    manual = {str(e["id"]) for e in (_load_yaml_text(pre + post).get("species") or [])}
    table, alias = _lj_table()
    _, sets, dip, _, _ = lj_context(thermo_path, species_path, transport_path, gri30_path)
    ents = []
    for s in names:
        if s in manual or s in EXCLUDED:
            continue
        lj = dict(lj_for(s, table, alias))
        lj.update({k: sets[s][k] for k in LJ_SET_ORDER[1:] if k in sets[s]})
        d = {"value": dip[s], "source": DIPOLE_SOURCE_GRI30} if s in dip else None
        ents.append(to_forge_entry(by[s], lj, d))
    buf = io.StringIO()
    dump_entries_yaml(ents, buf)
    return ents, buf.getvalue(), names, by


# ---------------------------------------------------------------------------------------------------------------
# 段 3 (#13-3): 手保守のエントリの MW・区間を CEA そのものへ (手保守の欄は残す)
# ---------------------------------------------------------------------------------------------------------------
def _gas_records(recs):
    """名前 → products 節・気相 (phase 0) の最初の記録。"""
    by = {}
    for r in recs:
        if r["section"] == "products" and r["phase"] == 0:
            by.setdefault(r["name"], r)
    return by


def _manual_blocks(text):
    """手保守のテキスト → (行の列, [(id, 開始行, 終了行 [排他])]) (`  - id: "X"` から次のエントリの直前まで)。"""
    lines = text.splitlines(keepends=True)
    starts = [k for k, L in enumerate(lines) if L.startswith("  - id: ")]
    out = []
    for n, k in enumerate(starts):
        end = starts[n + 1] if n + 1 < len(starts) else len(lines)
        out.append((str(_load_yaml_text(lines[k].strip()[2:])["id"]), k, end))
    return lines, out


def _intervals_lines(intervals):
    """区間列を手保守エントリの書式 (6 字下げ) で書く。数値は生成ブロックと同じ `_num`。"""
    out = ["    intervals:\n"]
    for lo, hi, a in intervals:
        out.append(f"      - Tlo: {_num(lo)}\n        Thi: {_num(hi)}\n")
        out.append("        coeffs: [" + ", ".join(_num(v) for v in a) + "]\n")
    return out


def sync_manual_text(pre, by):
    """手保守のテキスト (生成ブロックの BEGIN マーカまで) の、CEA に同名 (気相) がある種の `MW:` 行と `intervals:` 節を
    thermo.inp の値そのもの (区間数・区切り・係数・MW) に置き換える。それ以外の欄 (aliases・LJ・legacy_builtin・source・
    deviations・コメント) は一字も変えない。凝縮相 (`pair_of` を持つ) は係数を変えず、`MW:` だけを気液ペアの気相種と同じ値にする
    (気液ペアは同じ MW で質量換算する契約; ソルバが起動時に一致を検査する)。戻り値 (新テキスト, 置き換えた種の列)。"""
    lines, blocks = _manual_blocks(pre)
    ents = {str(e["id"]): e for e in (_load_yaml_text(pre + "\n").get("species") or [])}
    done = []
    for sid, k0, k1 in reversed(blocks):
        e = ents[sid]
        blk = lines[k0:k1]
        mw = [n for n, L in enumerate(blk) if L.startswith("    MW: ")]
        if e.get("phase") == "gas" and sid in by:
            r = by[sid]
            iv = [n for n, L in enumerate(blk) if L.startswith("    intervals:")]
            if len(mw) != 1 or len(iv) != 1 or not mw[0] < iv[0]:
                raise SystemExit(f"手保守エントリ {sid}: MW / intervals の行が 1 つずつ (MW が先) でない")
            n1 = iv[0] + 1
            while n1 < len(blk) and blk[n1].startswith("      "):
                n1 += 1
            blk = blk[:iv[0]] + _intervals_lines(r["intervals"]) + blk[n1:]
            blk[mw[0]] = f"    MW: {_num(mw_kg(r))}\n"
            done.append(sid)
        elif e.get("phase") == "condensed" and e.get("pair_of") in by:
            if len(mw) != 1:
                raise SystemExit(f"手保守エントリ {sid}: MW の行が 1 つでない")
            blk[mw[0]] = f"    MW: {_num(mw_kg(by[e['pair_of']]))}\n"
            done.append(sid)
        lines[k0:k1] = blk
    return "".join(lines), list(reversed(done))


def write_generated(thermo_path, species_path, transport_path, gri30_path):
    ents, body, _, _ = generated_entries(thermo_path, species_path, transport_path, gri30_path)
    _, sets, _, prov, _ = lj_context(thermo_path, species_path, transport_path, gri30_path)
    pre, _, post = split_generated(open(species_path, encoding="utf-8").read())
    pre, synced = sync_manual_text(pre, _gas_records(parse_thermo_inp(thermo_path)))
    pre, lj_synced = sync_lj_text(pre, sets, prov)
    with open(species_path, "w", encoding="utf-8") as f:
        f.write(pre + body + post)
    print(f"wrote {len(ents)} generated entries to {species_path}; hand-maintained entries synced to CEA: {synced}; "
          f"LJ sets gri30/svehla1962 synced: {lj_synced}")


def _iv_of(e):
    return [(iv["Tlo"], iv["Thi"], list(iv["coeffs"])) for iv in e["intervals"]]


def _bits_equal(a, b):
    """float の列・入れ子をビットで比較 (型も float であること; YAML が文字列・整数に読んだら不一致)。"""
    if isinstance(a, (list, tuple)):
        return isinstance(b, (list, tuple)) and len(a) == len(b) and all(_bits_equal(x, y) for x, y in zip(a, b))
    return isinstance(a, float) and isinstance(b, float) and a.hex() == b.hex()


def check_generated(thermo_path, species_path, transport_path, gri30_path):
    """往復のビット一致 (a)(b)(c) (docstring)。失敗数を返す。"""
    import io
    fail = 0

    def ck(ok, what):
        nonlocal fail
        print(("[PASS] " if ok else "[FAIL] ") + what)
        fail += 0 if ok else 1

    ents, body, names, by = generated_entries(thermo_path, species_path, transport_path, gri30_path)

    def same_as_cea(e, r):
        return (r is not None and _bits_equal(e["MW"], mw_kg(r))
                and _bits_equal(_iv_of(e), [(lo, hi, a) for lo, hi, a in r["intervals"]]))

    # (a) 同名の全種 (e- を含む) を書いて読み戻す
    buf = io.StringIO()
    dump_entries_yaml([to_forge_entry(by[s]) for s in names], buf)   # (a) は熱物性の往復だけを見る (LJ は (b) が見る)
    back = _load_yaml_text("species:\n" + buf.getvalue())["species"]
    bad = [s for s, e in zip(names, back) if str(e["id"]) != s or not same_as_cea(e, by[s])]
    nint = [len(by[s]["intervals"]) for s in names]
    ck(len(back) == len(names) and not bad,
       f"(a) 往復: 輸送データ ∩ thermo.inp の {len(names)} 種 (2 区間 {nint.count(2)}, 3 区間 {nint.count(3)}, 他 "
       f"{len(nint) - nint.count(2) - nint.count(3)}) の区間・全係数 ({sum(9 * n for n in nint)} 個)・MW が thermo.inp とビット一致"
       + (f" — 不一致 {bad}" if bad else ""))
    # (b) 共通データの生成ブロック = 今の生成結果
    _, cur_body, _ = split_generated(open(species_path, encoding="utf-8").read())
    ck(cur_body == body, f"(b) 共通データの生成ブロック ({len(ents)} 種) が生成結果と文字列一致")
    # (c) 共通データの生成ブロックを YAML として読み、thermo.inp と比較
    got = _load_yaml_text("species:\n" + cur_body).get("species") or []
    bad = [str(e["id"]) for e in got if not same_as_cea(e, by.get(str(e["id"])))]
    ck(not bad and len(got) == len(ents),
       f"(c) 共通データの生成ブロック {len(got)} 種の区間・全係数・MW が thermo.inp とビット一致" + (f" — 不一致 {bad}" if bad else ""))
    # (d) 手保守のエントリ (段 3): CEA に同名 (気相) がある種は区間・全係数・MW が thermo.inp とビット一致、
    #     凝縮相は MW が気液ペアの気相種とビット一致、手保守テキストが同期結果と文字列一致 (--write 済み)
    pre, _, _ = split_generated(open(species_path, encoding="utf-8").read())
    gas = _gas_records(parse_thermo_inp(thermo_path))
    man = _load_yaml_text(pre + "\n").get("species") or []
    mid = {str(e["id"]): e for e in man}
    same = [str(e["id"]) for e in man if e.get("phase") == "gas" and str(e["id"]) in gas]
    bad = [s for s in same if not same_as_cea(mid[s], gas[s])]
    nint = [len(gas[s]["intervals"]) for s in same]
    ck(not bad, f"(d) 手保守 {len(man)} 種のうち CEA と同名の気相 {len(same)} 種 {same} (2 区間 {nint.count(2)}, 3 区間 "
       f"{nint.count(3)}) の区間・全係数 ({sum(9 * n for n in nint)} 個)・MW が thermo.inp とビット一致" + (f" — 不一致 {bad}" if bad else ""))
    liq = [e for e in man if e.get("phase") == "condensed"]
    bad = [str(e["id"]) for e in liq if not (e.get("pair_of") in mid and _bits_equal(e["MW"], mid[e["pair_of"]]["MW"]))]
    ck(not bad, f"(d) 凝縮相 {[str(e['id']) for e in liq]} の MW が気液ペアの気相種とビット一致" + (f" — 不一致 {bad}" if bad else ""))
    ck(sync_manual_text(pre, gas)[0] == pre, "(d) 手保守のテキストが CEA 同期の結果と文字列一致 (--write 済み)")
    # (e) #14-L1: LJ 集合 gri30・svehla1962 (Svehla 転記 2 本の全値一致と頁参照は load_svehla が検査、不一致なら SystemExit)
    ids, sets, dip, prov, info = lj_context(thermo_path, species_path, transport_path, gri30_path)
    ck(sync_lj_text(pre, sets, prov)[0] == pre, "(e) 手保守の LJ_sets の gri30・svehla1962 行と provenance が生成結果と文字列一致 (--write 済み)")
    allent = {str(e["id"]): e for e in (_load_yaml_text(pre + "\n").get("species") or []) + got}
    bad = []
    for i in ids:
        cur = allent[i].get("LJ_sets") or {}
        for k in LJ_SET_ORDER[1:]:
            want = sets[i].get(k)
            have = cur.get(k)
            if (want is None) != (have is None) or (want and not _bits_equal([float(have["sigma"]), float(have["eps_kB"])],
                                                                             [want["sigma"], want["eps_kB"]])):
                bad.append(f"{i}.{k}")
    defs = _load_yaml_text(open(species_path, encoding="utf-8").read()).get("lj_sets") or {}
    pot = {k: (v or {}).get("potential") for k, v in defs.items()}
    ck(pot == LJ_SET_POTENTIAL, f"(e) トップレベル lj_sets の potential {pot} = {LJ_SET_POTENTIAL} (#14-L1b)")
    ng, ns = sum(1 for i in ids if "gri30" in sets[i]), sum(1 for i in ids if "svehla1962" in sets[i])
    ck(not bad, f"(e) 全気相 {len(ids)} 種の LJ_sets.gri30 ({ng} 種)・svehla1962 ({ns} 種) が gri30.yaml・Svehla 転記とビット一致"
       + (f" — 不一致 {bad}" if bad else ""))
    # 双極子: gri30 の値 = 共通データの種レベルの dipole (手保守の種も含む)、gri30 に無い種に dipole があれば手保守として表示
    dbad = [i for i in ids if i in dip and float((allent[i].get("dipole") or {}).get("value", -1.0)) != dip[i]]
    ck(not dbad, f"(e) gri30 の双極子 {', '.join(f'{i} {v}' for i, v in dip.items())} D が共通データの種レベル dipole と一致"
       + (f" — 不一致 {dbad}" if dbad else ""))
    extra = [i for i in ids if allent[i].get("dipole") and i not in dip]
    none = [i for i in ids if not (allent[i].get("LJ_sets") or {})]
    only_legacy = [i for i in ids if set(allent[i].get("LJ_sets") or {}) == {"legacy_v1"}]
    print(f"[INFO] LJ: どの集合にも無い種 {none}; legacy_v1 だけの種 {only_legacy}; gri30 に無い dipole {extra}; "
          f"共通データに無い Svehla 分子 {info['svehla_unmapped']}; 共通データに無い gri30 種 {len(info['gri30_unmapped'])}")
    other = [str(e["id"]) for e in man if str(e["id"]) not in same and e.get("phase") == "gas"]
    print(f"[INFO] 手保守で CEA に同名の気相が無い種 (係数はそのまま): {other}")
    lj = [str(e["id"]) for e in got if (e.get("LJ_sets") or {}).get("legacy_v1")]
    print(f"[INFO] 生成ブロック: legacy_v1 あり {len(lj)} 種 ({', '.join(lj)}), なし {len(got) - len(lj)} 種; "
          f"除外 {list(EXCLUDED)}; 手保守 (生成ブロック外) で CEA と同名 {[s for s in names if s not in EXCLUDED and s not in {str(e['id']) for e in got}]}")
    return fail


# ---------------------------------------------------------------------------------------------------------------
# NASA-9 評価 (double; ソルバ不使用。区間の選び方は T < 区間上端 の最初の区間、範囲外は端の区間で外挿)
# ---------------------------------------------------------------------------------------------------------------
def _pick(intervals, T):
    for lo, hi, a in intervals:
        if T < hi:
            return a
    return intervals[-1][2]


def cp_R(intervals, T):
    a = _pick(intervals, T)
    return a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4


def h_RT(intervals, T):
    a = _pick(intervals, T)
    return (-a[0] / T**2 + a[1] * math.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3 + a[5] * T**3 / 4
            + a[6] * T**4 / 5 + a[7] / T)


def s_R(intervals, T):
    a = _pick(intervals, T)
    return (-a[0] / (2 * T**2) - a[1] / T + a[2] * math.log(T) + a[3] * T + a[4] * T**2 / 2 + a[5] * T**3 / 3
            + a[6] * T**4 / 4 + a[8])


# ---------------------------------------------------------------------------------------------------------------
# 監査
# ---------------------------------------------------------------------------------------------------------------
def _load_yaml(path):
    import yaml
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _cur_intervals(e):
    return [(float(iv["Tlo"]), float(iv["Thi"]), [float(v) for v in iv["coeffs"]]) for iv in e["intervals"]]


def _rel(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def _fmt_T(x):
    return f"{x:.10g}"


def compare_coeffs(cur_iv, cea_iv):
    """現行区間 (n 区間) と CEA の先頭 n 区間を区間ごとに比較。返り値: [(区間番号, 境界一致, ビット一致数/9, 最大相対差, 最大の係数)]。"""
    out = []
    for k, (lo, hi, a) in enumerate(cur_iv):
        if k >= len(cea_iv):
            out.append((k, False, 0, float("nan"), "-"))
            continue
        clo, chi, ca = cea_iv[k]
        nbit = sum(1 for x, y in zip(a, ca) if x == y)
        rels = [_rel(x, y) for x, y in zip(a, ca)]
        j = max(range(9), key=lambda m: rels[m])
        out.append((k, (lo, hi) == (clo, chi), nbit, rels[j], f"a{j}: {a[j]!r} vs {ca[j]!r}" if rels[j] > 0 else ""))
    return out


def delta_props(cur_iv, cur_mw, cea_iv, cea_mw, T_lo=200.0, T_hi=6000.0, dT=1.0):
    """200 ≤ T < 6000 K を dT 刻み (+ 区切り直下/直上) で走査し、max |Δcp|/cp、max |Δh| [J/kg]、max |Δs°| [J/(kg K)] を返す。
    h, s° は各側の MW で質量あたりに換算 (MW 差の効果を含む)。区切り 1000 K は両側 (999.999, 1000) を含む。
    T = 6000 K ちょうどは含めない: 現行 (2 区間) は上端を閉区間で高温区間、CEA 3 区間は 6000 K が第 3 区間の下端になり、
    どちらを取るかは段 1 の区間選択の規約で決まる (その段差は `boundary_jumps` で別に出す)。"""
    Ts = [T_lo + k * dT for k in range(int((T_hi - T_lo) / dT))]      # T_hi 自体は含めない (区間端の規約に依存するため別表)
    Ts += [999.999, 1000.0, 1000.001, 5999.999]
    best = {"cp": (0.0, None), "h": (0.0, None), "s": (0.0, None), "h_sens": (0.0, None)}
    Rc, Rn = RU / cur_mw, RU / cea_mw
    hc298, hn298 = h_RT(cur_iv, 298.15) * 298.15 * Rc, h_RT(cea_iv, 298.15) * 298.15 * Rn
    for T in Ts:
        cpc, cpn = cp_R(cur_iv, T) * Rc, cp_R(cea_iv, T) * Rn
        hc, hn = h_RT(cur_iv, T) * T * Rc, h_RT(cea_iv, T) * T * Rn
        sc, sn = s_R(cur_iv, T) * Rc, s_R(cea_iv, T) * Rn
        vals = {"cp": abs(cpn - cpc) / abs(cpc), "h": abs(hn - hc), "s": abs(sn - sc),
                "h_sens": abs((hn - hn298) - (hc - hc298))}
        for k, v in vals.items():
            if v > best[k][0]:
                best[k] = (v, T)
    return best


def boundary_jumps(intervals):
    """区間境界ごとの (T, cp/R 段差, h/RT 段差, s/R 段差) = 上側区間 − 下側区間 (同じ T で評価)。"""
    out = []
    for k in range(len(intervals) - 1):
        T = intervals[k][1]
        lo, hi = [intervals[k]], [intervals[k + 1]]
        out.append((T, cp_R(hi, T) - cp_R(lo, T), h_RT(hi, T) - h_RT(lo, T), s_R(hi, T) - s_R(lo, T)))
    return out


def audit(thermo_path, species_path, transport_path, sern_root=None):
    recs = parse_thermo_inp(thermo_path)
    by_name = {}
    dup = []
    for r in recs:
        if r["name"] in by_name:
            dup.append(r["name"])
        by_name.setdefault(r["name"], r)
    prod = [r for r in recs if r["section"] == "products"]
    tr = _load_yaml(transport_path)
    sd = _load_yaml(species_path)
    tr_ids = [str(s["id"]) for s in tr["species"]]
    cur = {str(e["id"]): e for e in sd["species"]}

    W = []
    p = W.append
    p(f"- thermo.inp: `{thermo_path}` (SHA-256 `{sha256_file(thermo_path)}`), 記録 {len(recs)} "
      f"(products {len(prod)}, reactants {len(recs) - len(prod)})、完全同名の重複記録 {len(dup)} {dup if dup else ''}")
    p(f"- 輸送データ: `{os.path.relpath(transport_path, REPO)}` {len(tr_ids)} 種、共通データ: `{os.path.relpath(species_path, REPO)}` "
      f"{len(cur)} 種")
    p("- 区間の選び方 (§6 の Δ): `T < Thi` の最初の区間 (ソルバの `Tc < Tmid ? low : high` と同じく境界ちょうどは上側)。"
      f"R_u = {RU} J/(mol K) (`THERMO_RU`)。ソルバは使わない (Python double)")
    p("")

    # ---- 1. trans 66 種と thermo.inp の同名
    have = [s for s in tr_ids if s in by_name]
    miss = [s for s in tr_ids if s not in by_name]
    p("## 1. 輸送 66 種と thermo.inp の同名 (大小文字区別)")
    p("")
    p(f"- 輸送データの種 {len(tr_ids)}、うち thermo.inp に同名あり **{len(have)}**、無し **{len(miss)}**: "
      + ", ".join(f"`{s}`" for s in miss))
    for s in miss:
        ci = [r["name"] for r in recs if r["name"].upper() == s.upper()]
        p(f"  - `{s}`: 大小文字を無視した一致 {ci if ci else 'なし'}")
    sec = [s for s in have if by_name[s]["section"] != "products"]
    cond = [s for s in have if by_name[s]["phase"] != 0]
    p(f"- 同名の記録が reactants 節にしか無い種: {sec if sec else 'なし'}。凝縮相 (phase ≠ 0) の同名: {cond if cond else 'なし'}")
    p("")

    # ---- 2. 区間数・区切り・MW
    p(f"## 2. 同名 {len(have)} 種の区間数・区切り温度・MW")
    p("")
    from collections import Counter
    cnt = Counter(len(by_name[s]["intervals"]) for s in have)
    p("- 区間数の内訳: " + ", ".join(f"{k} 区間 {v} 種" for k, v in sorted(cnt.items())))
    nonstd = []
    for s in have:
        b = [by_name[s]["intervals"][0][0]] + [iv[1] for iv in by_name[s]["intervals"]]
        if any(abs(x - STD_BOUNDS[k]) > 1e-9 for k, x in enumerate(b) if k < len(STD_BOUNDS)) or len(b) > 4:
            nonstd.append((s, b))
        # 区間の連続性
        ivs = by_name[s]["intervals"]
        for k in range(len(ivs) - 1):
            if ivs[k][1] != ivs[k + 1][0]:
                nonstd.append((s, f"不連続 {ivs[k][1]} → {ivs[k + 1][0]}"))
    p("- 境界が 200/1000/6000/20000 の前方一致でない種: "
      + ("; ".join(f"`{s}` {b}" for s, b in nonstd) if nonstd else "なし"))
    p("")
    p("- 右 2 列は区間境界での段差 (上側区間 − 下側区間を同じ T で評価): `Δcp/R` と `Δh/RT`、境界ごとに `/` で区切る。")
    p("")
    p("| # | 種 | 区間数 | 区切り温度 [K] | MW [g/mol] (欄) | 相 | 日付コード | 境界の Δcp/R | 境界の Δh/RT |")
    p("|---|---|---|---|---|---|---|---|---|")
    for k, s in enumerate(have, 1):
        r = by_name[s]
        b = [r["intervals"][0][0]] + [iv[1] for iv in r["intervals"]]
        j = boundary_jumps(r["intervals"])
        p(f"| {k} | `{s}` | {len(r['intervals'])} | {'–'.join(_fmt_T(x) for x in b)} | {r['MW_str']} | {r['phase']} | {r['date']} "
          f"| {' / '.join(f'{x[1]:.1e}' for x in j)} | {' / '.join(f'{x[2]:.1e}' for x in j)} |")
    p("")

    # ---- 3. 既存 7 種 (AIR 除く) + H2O(L)
    solver7 = [i for i, e in cur.items() if "solver" in (e.get("legacy_builtin") or [])]
    targets3 = [i for i in solver7 if i != "AIR"] + ["H2O(L)"]
    p("## 3. 内蔵 (`legacy_builtin: solver`) 種と `H2O(L)`: 現行共通データ vs thermo.inp")
    p("")
    p(f"- `legacy_builtin: solver` の種: {solver7}。`AIR` は CEA に同名が無く (CEA は `Air`, reactants 節) 対象外。")
    p("- 係数は区間ごとに 9 個のビット一致数と最大相対差。現行の区間数 vs CEA の区間数も示す (現行は先頭 n 区間と比べる)。")
    p("")
    p("| 種 | 現行区間数 / CEA | MW 現行 [kg/mol] | MW CEA (欄 e-3) | MW 相対差 | 区間 | 境界一致 | ビット一致 | 最大相対差 | 最大の係数 |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for s in targets3:
        e, r = cur[s], by_name[s]
        civ = _cur_intervals(e)
        mwc, mwn = float(e["MW"]), mw_kg(r)
        for (k, bnd, nbit, mr, worst) in compare_coeffs(civ, r["intervals"]):
            head = (f"`{s}` | {len(civ)} / {len(r['intervals'])} | {mwc!r} | {mwn!r} | {_rel(mwc, mwn):.2e}"
                    if k == 0 else " | | | | ")
            p(f"| {head} | {k} ({_fmt_T(civ[k][0])}–{_fmt_T(civ[k][1])}) | {'yes' if bnd else 'NO'} | {nbit}/9 | {mr:.3e} | {worst} |")
    p("")
    # Ar 高温区間の詳細
    ar = by_name["Ar"]
    p("### 3a. `Ar` 高温区間 (1000–6000 K) の係数")
    p("")
    p("| 係数 | 現行 | thermo.inp |")
    p("|---|---|---|")
    for j in range(9):
        p(f"| a{j} | {_cur_intervals(cur['Ar'])[1][2][j]!r} | {ar['intervals'][1][2][j]!r} |")
    p("")
    if len(ar["intervals"]) > 2:
        lo, hi, a = ar["intervals"][2]
        p(f"- thermo.inp の第 3 区間 {lo:g}–{hi:g} K (現行データに無い): a = {a}")
        p("")

    # ---- 4. 6 種と SERN 外部 DB
    six = ["CO", "H2", "OH", "H", "NO", "O"]
    p("## 4. `CO/H2/OH/H/NO/O`: 現行共通データ・SERN 外部 DB vs thermo.inp")
    p("")
    p("| 種 | 現行区間数 / CEA | MW 現行 | MW CEA | MW 相対差 | 区間 | 境界一致 | ビット一致 | 最大相対差 | 最大の係数 |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for s in six:
        e, r = cur[s], by_name[s]
        civ = _cur_intervals(e)
        mwc, mwn = float(e["MW"]), mw_kg(r)
        for (k, bnd, nbit, mr, worst) in compare_coeffs(civ, r["intervals"]):
            head = (f"`{s}` | {len(civ)} / {len(r['intervals'])} | {mwc!r} | {mwn!r} | {_rel(mwc, mwn):.2e}"
                    if k == 0 else " | | | | ")
            p(f"| {head} | {k} ({_fmt_T(civ[k][0])}–{_fmt_T(civ[k][1])}) | {'yes' if bnd else 'NO'} | {nbit}/9 | {mr:.3e} | {worst} |")
    p("")
    if sern_root:
        p("### 4a. SERN の外部 DB (`species_db_external.yaml`) の中身")
        p("")
        p(sern_external_section(sern_root, six, by_name, cur))
        p("")

    # ---- 5. 大小文字だけ違う名前
    p("## 5. 大小文字だけ違う名前の組 (thermo.inp 全体)")
    p("")
    groups = {}
    for r in recs:
        groups.setdefault(r["name"].upper(), [])
        if r["name"] not in groups[r["name"].upper()]:
            groups[r["name"].upper()].append(r["name"])
    coll = {k: v for k, v in groups.items() if len(v) > 1}
    p(f"- 組の数: {len(coll)}")
    p("")
    p("| 大文字化 | thermo.inp の名前 (節・相) | 輸送 66 種に含まれる名前 |")
    p("|---|---|---|")
    for k, v in sorted(coll.items()):
        desc = ", ".join(f"`{n}` ({by_name[n]['section'][0]}, {by_name[n]['phase']})" for n in v)
        p(f"| {k} | {desc} | {', '.join(f'`{n}`' for n in v if n in tr_ids) or '—'} |")
    inset = [v for v in coll.values() if sum(1 for n in v if n in have) > 1]
    p("")
    p(f"- 同名 {len(have)} 種の**範囲内**で大小文字だけ違う組: {inset if inset else 'なし'}")
    touch = [v for v in coll.values() if any(n in have for n in v)]
    p(f"- {len(have)} 種のどれかを含む組 (相手は範囲外): {touch if touch else 'なし'}")
    # 現行の別名と thermo.inp 名の衝突
    ali = []
    for i, e in cur.items():
        for a in (e.get("aliases") or []):
            hit = [r["name"] for r in recs if r["name"].upper() == str(a).upper()]
            if hit:
                ali.append(f"`{i}` の別名 `{a}` ↔ thermo.inp {hit}")
    p("- 現行共通データの別名と thermo.inp の名前の大小文字無視一致: " + ("; ".join(ali) if ali else "なし"))
    p("")

    # ---- 6. Δ 概算
    p("## 6. 段 3 の予測の根拠: 既存 (`legacy_builtin: solver`, `AIR` 除く) を CEA 化したときの差 (200–6000 K)")
    p("")
    p("- 走査は 200 ≤ T < 6000 K (1 K 刻み + 999.999/1000/1000.001/5999.999 K)。6000 K ちょうどは区間選択の規約次第なので下の別表。")
    p("- 現行 (区間・MW) と CEA (全区間・MW 欄) を質量あたりで比較。`Δh` は絶対 (a7 込み)、`Δ(h−h298)` は 298.15 K 基準の顕熱差、"
      "`Δs°` は J/(kg K)。括弧内は最大を与えた T [K]。")
    p("")
    p("| 種 | max \\|Δcp\\|/cp | max \\|Δh\\| [J/kg] | max \\|Δ(h−h298)\\| [J/kg] | max \\|Δs°\\| [J/(kg K)] |")
    p("|---|---|---|---|---|")
    for s in [i for i in solver7 if i != "AIR"]:
        b = delta_props(_cur_intervals(cur[s]), float(cur[s]["MW"]), by_name[s]["intervals"], mw_kg(by_name[s]))
        p(f"| `{s}` | {b['cp'][0]:.3e} ({_fmt_T(b['cp'][1]) if b['cp'][1] else '-'}) "
          f"| {b['h'][0]:.3e} ({_fmt_T(b['h'][1]) if b['h'][1] else '-'}) "
          f"| {b['h_sens'][0]:.3e} ({_fmt_T(b['h_sens'][1]) if b['h_sens'][1] else '-'}) "
          f"| {b['s'][0]:.3e} ({_fmt_T(b['s'][1]) if b['s'][1] else '-'}) |")
    p("")
    p("- 6000 K ちょうど: 現行は高温区間 (閉区間)。CEA 3 区間の種で第 3 区間を選んだときの段差 (第 3 − 第 2 区間):")
    p("")
    p("| 種 | Δcp/R @6000 | Δh/RT @6000 | Δs°/R @6000 | Δh [J/kg] @6000 |")
    p("|---|---|---|---|---|")
    for s in [i for i in solver7 if i != "AIR"]:
        ivs = by_name[s]["intervals"]
        if len(ivs) < 3:
            p(f"| `{s}` | (CEA も 2 区間: 段差なし) | | | |")
            continue
        T, dcp, dh, ds = boundary_jumps(ivs)[1]
        p(f"| `{s}` | {dcp:.2e} | {dh:.2e} | {ds:.2e} | {dh * T * RU / mw_kg(by_name[s]):.3e} |")
    # 参考: MW だけ変えた場合 (係数は現行) — H2O / He の MW 効果の切り分け
    p("")
    p("- 参考 (MW だけ CEA、係数は現行): " + "; ".join(
        f"`{s}` Δcp/cp {delta_props(_cur_intervals(cur[s]), float(cur[s]['MW']), _cur_intervals(cur[s]), mw_kg(by_name[s]))['cp'][0]:.3e}"
        for s in ["H2O", "He"]))
    # 参考: H2O(L) (273.15–373.15 K) の MW 効果
    e = cur["H2O(L)"]
    r = by_name["H2O(L)"]
    mwc, mwn = float(e["MW"]), mw_kg(r)
    p(f"- 参考 `H2O(L)`: MW {mwc!r} → {mwn!r} (相対 {_rel(mwc, mwn):.2e}); 係数は上の 3 の表のとおり")
    return "\n".join(W) + "\n"


def sern_external_section(sern_root, six, by_name, cur):
    """SERN 設計ワークツリーの生成元 (design/forge_design/gas) を import し、ソルバ内蔵で解決できない 6 種について
    runner が書く `species_db_external.yaml` のテキストを生成元の関数そのもの (`species_db_raw_yaml`) で作って照合する。"""
    import importlib
    import yaml
    out = []
    design = os.path.join(sern_root, "design")
    sys.path.insert(0, design)
    try:
        comp = importlib.import_module("forge_design.gas.composition")
        sp = importlib.import_module("forge_design.gas.semiperfect")
    finally:
        sys.path.pop(0)
    builtin = comp.solver_builtin_names()
    db = comp.ResolvedSpeciesDB.builtin()
    ext = {s: db.entries[s] for s in six if s in db.entries and s.upper() not in builtin}
    text = comp.species_db_raw_yaml(ext)
    raw = yaml.safe_load(text)
    out.append(f"- 生成元: `{sern_root}` (`design/forge_design/gas/composition.py` `species_db_raw_yaml`、"
               f"データは `{os.path.relpath(str(sp.SPECIES_DATA_FILE), sern_root)}` の `legacy_builtin: design`; "
               f"SHA-256 `{sha256_file(str(sp.SPECIES_DATA_FILE))}`)。実 run の `species_db_external.yaml` はワークツリー内に無かった"
               "ので、runner と同じ関数でテキストを作って読み戻した。")
    out.append(f"- ソルバ内蔵で解決できないとして外部 DB に出る種: {list(ext)}")
    out.append("")
    out.append("| 種 | 外部 DB vs 現行共通データ (MW・18 係数・区切り) | 外部 DB vs thermo.inp 先頭 2 区間 (ビット一致 / 最大相対差) | LJ 外部 DB (σ, ε) |")
    out.append("|---|---|---|---|")
    for s in ext:
        e = raw[s]
        ext_iv = [(float(e["Tlo"]), float(e["Tmid"]), [float(v) for v in e["nasa9_low"]]),
                  (float(e["Tmid"]), float(e["Thi"]), [float(v) for v in e["nasa9_high"]])]
        civ = _cur_intervals(cur[s])
        same_cur = (float(e["MW"]) == float(cur[s]["MW"]) and ext_iv == civ)
        cmp = compare_coeffs(ext_iv, by_name[s]["intervals"])
        nb = sum(c[2] for c in cmp)
        mr = max(c[3] for c in cmp)
        out.append(f"| `{s}` | {'ビット一致' if same_cur else '不一致'} | {nb}/18, {mr:.2e}; MW 相対 {_rel(float(e['MW']), mw_kg(by_name[s])):.2e} "
                   f"| ({e['LJ_sigma']}, {e['LJ_eps_kB']}) |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--thermo", default=DEFAULT_THERMO, help="CEA thermo.inp (既定: <repo>/.venv-cea/nasa_cea/thermo.inp)")
    ap.add_argument("--species-data", default=DEFAULT_SPECIES, help="共通データ (--audit/--check は読むだけ、--write は生成ブロックと手保守エントリの MW・区間を書き換える)")
    ap.add_argument("--transport-data", default=DEFAULT_TRANSPORT, help="輸送データ (種の一覧)")
    ap.add_argument("--audit", action="store_true", help="監査表 (Markdown) を出す (段 0)")
    ap.add_argument("--audit-out", help="監査表の出力先 (既定: 標準出力)")
    ap.add_argument("--write", action="store_true", help="共通データの生成ブロックを書き換える (段 2)")
    ap.add_argument("--check", action="store_true", help="往復のビット一致と生成ブロック・手保守エントリの一致を検査する (段 2・3)")
    ap.add_argument("--sern-root", help="SERN 設計ワークツリー (外部 DB の生成元を import して照合; 読むだけ)")
    ap.add_argument("--gri30", default=default_gri30(), help="Cantera 同梱 gri30.yaml (LJ 集合 gri30 の入力; --write/--check)")
    a = ap.parse_args()
    if not os.path.exists(a.thermo):
        raise SystemExit(f"thermo.inp が無い: {a.thermo} (--thermo で指定)")
    if a.write:
        write_generated(a.thermo, a.species_data, a.transport_data, a.gri30)
        return
    if a.check:
        fail = check_generated(a.thermo, a.species_data, a.transport_data, a.gri30)
        print("ALL PASS" if fail == 0 else f"FAIL ({fail})")
        sys.exit(1 if fail else 0)
    if not a.audit:
        raise SystemExit("--audit / --write / --check のどれかを指定する")
    text = audit(a.thermo, a.species_data, a.transport_data, a.sern_root)
    if a.audit_out:
        with open(a.audit_out, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"wrote {a.audit_out}")
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
