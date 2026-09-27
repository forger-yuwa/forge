#!/usr/bin/env python3
"""run ディレクトリの化学種配置 (名前 → index / MW / 凝縮種 / トレーサ) を 1 か所で解決する (importable + CLI)。

後処理 (ParaView `Forge Saturation`, axis_csv_va.py, 種変換 restart) が `Y1` などの index を決め打ちせず、
`solverConfig.yaml` (`physProp.species`, `condensation`) と `species_db.yaml` (MW) から**名前で**引くための共通関数。
`species_meta.yaml` (設計チェーンが書く機械可読メタ) があれば凝縮種名・lump 展開行列もそこから取る。
plans/active/thermophysics-cea-mole-fraction-species.md §4.4。

  python3 forge_species.py RUN_DIR            # 表示
  python3 forge_species.py RUN_DIR --json     # JSON

Python:
  from forge_species import species_info
  info = species_info(run_dir)
  info["names"]            # ['MIXDRY', 'H2O'] (physProp.species の順 = index)
  info["index"]["H2O"]     # 1
  info["MW"]["H2O"]        # 0.0180153 (species_db.yaml; 無ければ内蔵表)
  info["condensing"]       # 'H2O' | None (凝縮 ON のときだけ),  info["condensing_index"]  # 1 | None
  info["h2o_index"]        # H2O (または condensationSpecies で指定した種) の index。凝縮 ON/OFF と無関係。無ければ None
  info["h2o_array"]        # 'Y{h2o_index}' | None
  info["vapor_array"]      # 凝縮 ON なら凝縮種の配列、OFF でも H2O が種にあればその配列 (後処理の分圧は常に名前で解決; codex M9)
  info["tracer"]           # 'exhaust' | None
  species_signature(run_dir)   # restart 照合用: config + species_db.yaml (名前/順序, MW, NASA-9 両域係数, Tlo/Tmid/Thi, LJ, datum, tracer)。
                               # 内蔵種は係数を持たない (None) ので、比較では「照合不能」になる (現在の内蔵表から作り直さない)。
                               # species_meta.yaml の species が config と矛盾すれば例外
  compare_signatures(a, b)     # 不一致・照合不能の説明 list (空なら一致; 係数は rel 1e-12)。source (builtin/file) では比較しない
  required_conserved(sig)      # restart に必須の VALUE/ データセット名

  解決済み記録 (ソルバ出力 resolved_species_<互換16>[_<完全性16>].yaml; plans/active/thermophysics-solver-owned-species-db.md §4.3):
  load_record(path)            # 記録を読み、完全性ハッシュ (全文 SHA-256) と互換性ハッシュ (中身から再計算) を検証した dict
  compat_text(...)             # 互換性ハッシュの正規化テキスト (C++ speciesDB_compatText と一字一句同じ)
  field_species_attrs(h5)      # res_*.h5 / 入力 h5 のルート属性 species_hash / species_record_sha256 / species_input_unverified
  find_record(h5, dirs)        # 属性の完全性ハッシュに一致する記録を h5 の隣 (と dirs) から探す
  signature_from_record(rec)   # 記録から compare_signatures 用の署名 (全種の係数込み)

YAML の注意: 種名 `NO` / `N` / `Y` は PyYAML の既定では真偽値になる (design チェーンの古い config は無引用)。
本モジュールの `load_yaml_str` は真偽値の暗黙解決を外した SafeLoader で読むので `NO` は文字列のまま。
"""
import argparse, hashlib, json, os, re, sys
import yaml


class _StrSafeLoader(yaml.SafeLoader):
    """真偽値の暗黙解決を true/false (大文字小文字違い含む) だけに絞った SafeLoader。YAML 1.1 の yes/no/y/n/on/off は
    文字列のまま (種名 NO / N / Y を壊さない)。"""


_BOOL_TAG = "tag:yaml.org,2002:bool"
_BOOL_RX = re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$")
_StrSafeLoader.yaml_implicit_resolvers = {
    ch: [(tag, rx) for (tag, rx) in rs if tag != _BOOL_TAG] + ([(_BOOL_TAG, _BOOL_RX)] if ch in "tTfF" else [])
    for ch, rs in yaml.SafeLoader.yaml_implicit_resolvers.items()
}


def load_yaml_str(path):
    """真偽値を文字列のまま保つ YAML 読込 (solverConfig.yaml / species_db.yaml / species_meta.yaml 用)。"""
    with open(path) as f:
        return yaml.load(f, Loader=_StrSafeLoader)

# thermo_d.cu / speciesDB.cpp 内蔵 DB の MW [kg/mol] (speciesDBFile が無い run 用)。別名込み。
BUILTIN_MW = {"N2": 0.0280134, "O2": 0.0319988, "AR": 0.039948, "CO2": 0.0440095, "HE": 0.0040026,
              "H2O": 0.0180153, "WATER": 0.0180153, "AIR": 0.0289647}


def _find_ci(d, name):
    """dict d を大文字小文字無視で引く (完全一致優先)。"""
    if name in d:
        return d[name]
    up = str(name).upper()
    for k, v in d.items():
        if str(k).upper() == up:
            return v
    return None


def species_info(run_dir):
    """run_dir の solverConfig.yaml / species_db.yaml / species_meta.yaml から種配置を返す。"""
    run_dir = os.path.abspath(run_dir)
    cfg_path = os.path.join(run_dir, "solverConfig.yaml")
    if not os.path.exists(cfg_path):
        raise FileNotFoundError(f"{cfg_path} が無い")
    cfg = load_yaml_str(cfg_path)
    pp = cfg.get("physProp") or {}
    tm = int(pp.get("thermalMethod", 0))
    names = [str(s) for s in (pp.get("species") or [])]
    if tm == 2 and not names:
        names = ["N2"]                      # solverConfig.cpp の既定 (単成分 N2)
    db = None
    db_file = pp.get("speciesDBFile")
    if db_file:
        p = db_file if os.path.isabs(db_file) else os.path.join(run_dir, db_file)
        if os.path.exists(p):
            db = load_yaml_str(p) or {}
        else:
            raise FileNotFoundError(f"speciesDBFile {p} が無い")
    MW = {}
    for n in names:
        e = _find_ci(db, n) if db else None
        if e is not None and "MW" in e:
            MW[n] = float(e["MW"])
        else:
            mw = _find_ci(BUILTIN_MW, n)
            if mw is None:
                raise KeyError(f"species {n}: MW が species_db.yaml にも内蔵表にも無い")
            MW[n] = float(mw)
    index = {n: i for i, n in enumerate(names)}
    meta_path = os.path.join(run_dir, "species_meta.yaml")
    meta = load_yaml_str(meta_path) if os.path.exists(meta_path) else None
    if meta and meta.get("species") and [str(s) for s in meta["species"]] != names:
        print(f"[forge_species] warning: species_meta.yaml の species {meta['species']} と solverConfig の {names} が違う"
              " (solverConfig を正とする)", file=sys.stderr)

    cond = cfg.get("condensation") or {}
    condensing = None
    if int(cond.get("condensation", 0)) == 1:
        cname = cond.get("condensationSpecies")
        cidx = cond.get("condGasSpecies")
        if cname is not None:
            hit = [n for n in names if n.upper() == str(cname).upper()]
            if not hit:
                raise KeyError(f"condensationSpecies {cname} が physProp.species {names} に無い")
            condensing = hit[0]
            if cidx is not None and int(cidx) != index[condensing]:
                raise ValueError(f"condGasSpecies={cidx} と condensationSpecies={cname} (index {index[condensing]}) が食い違う")
        elif cidx is not None and int(cidx) >= 0:
            if int(cidx) >= len(names):
                raise IndexError(f"condGasSpecies={cidx} が physProp.species {names} の範囲外")
            condensing = names[int(cidx)]
        elif meta and meta.get("condensing_species"):
            condensing = str(meta["condensing_species"])
    tracer = str(pp.get("tracer") or "").strip() or None
    if tracer in ("none", "0"):
        tracer = None
    # H2O (または condensationSpecies で名指しされた種) の index は凝縮 ON/OFF と無関係に名前で解決する (codex 2026-09-16 M9:
    # 凝縮 OFF run の分圧を固定組成に戻すと ~24 % 誤る)。
    h2o_name = None
    cname = cond.get("condensationSpecies")
    for cand in ([str(cname)] if cname else []) + ["H2O", "WATER"]:
        hit = [n for n in names if n.upper() == cand.upper()]
        if hit:
            h2o_name = hit[0]; break
    h2o_index = index.get(h2o_name) if h2o_name else None
    vapor_array = (f"Y{index[condensing]}" if condensing else (f"Y{h2o_index}" if h2o_index is not None else None))
    return {
        "run_dir": run_dir,
        "thermalMethod": tm,
        "names": names,
        "index": index,
        "MW": MW,
        "condensing": condensing,
        "condensing_index": index.get(condensing) if condensing else None,
        "h2o_index": h2o_index,
        "h2o_array": (f"Y{h2o_index}" if h2o_index is not None else None),
        "vapor_array": vapor_array,
        "tracer": tracer,
        "condensation": int(cond.get("condensation", 0)) == 1,
        "condModel": int(cond.get("condModel", 0)),
        "thermoHrefTemp": float(pp.get("thermoHrefTemp", 0.0)),
        "speciesDBFile": db_file,
        "meta": meta,
    }


def species_signature(run_dir):
    """restart 照合用の署名 (codex 2026-09-16 result-2 M2): 実際の solverConfig.yaml と解決済み species_db.yaml から作る。
    {thermalMethod, names (順序込み), tracer, thermoHrefTemp, species: {name: {MW, Tlo, Tmid, Thi, nasa9_low, nasa9_high, source}}}。
    speciesDBFile に無い種 (内蔵 DB) は係数を持たないので source="builtin" (MW は内蔵表) とし、両側 builtin なら同一とみなす。
    species_meta.yaml があれば species の名前/順序と tracer.enabled が config と一致することを要求し、矛盾は例外 (refuse)。
    解決できなければ例外 (呼び手が既定でエラーにする)。CPG (thermalMethod≠2) は names=[]。"""
    info = species_info(run_dir)
    run_dir = info["run_dir"]
    meta = info["meta"]
    if meta is not None and meta.get("species") is not None:
        mnames = [str(x) for x in meta["species"]]
        if mnames != info["names"]:
            raise ValueError(f"{run_dir}: species_meta.yaml の species {mnames} と solverConfig.yaml の physProp.species {info['names']} が矛盾する")
        mt = (meta.get("tracer") or {}).get("enabled")
        if mt is not None and bool(mt) != bool(info["tracer"]):
            # トレーサの有無は保存量 roXi の要否を決めるので、meta と config の矛盾は拒否する (codex 2026-09-16 result-3 M1)。
            raise ValueError(f"{run_dir}: species_meta.yaml tracer.enabled={mt} と solverConfig.yaml physProp.tracer={info['tracer']} が矛盾する"
                             " (species_meta.yaml を config に合わせて直すこと)")
    db = {}
    if info["speciesDBFile"]:
        p = info["speciesDBFile"] if os.path.isabs(info["speciesDBFile"]) else os.path.join(run_dir, info["speciesDBFile"])
        db = {str(k): v for k, v in (load_yaml_str(p) or {}).items()}
    species = {}
    for n in info["names"]:
        e = _find_ci(db, n)
        if e is not None and "nasa9_low" in e and "nasa9_high" in e:
            lo, hi = [float(x) for x in e["nasa9_low"]], [float(x) for x in e["nasa9_high"]]
            if len(lo) != 9 or len(hi) != 9:
                raise ValueError(f"{run_dir}: species_db.yaml {n} の nasa9 係数が 9 個でない")
            # LJ の既定値は C++ speciesDB_resolve と同じ (3.6 Å / 97 K)
            species[n] = {"MW": float(e["MW"]), "Tlo": float(e.get("Tlo", 200.0)), "Tmid": float(e.get("Tmid", 1000.0)),
                          "Thi": float(e.get("Thi", 6000.0)), "LJ_sigma": float(e.get("LJ_sigma", 3.6)),
                          "LJ_eps_kB": float(e.get("LJ_eps_kB", 97.0)), "nasa9_low": lo, "nasa9_high": hi, "source": "file"}
        else:
            # 内蔵種: 係数を Python 側で作り直さない (plan §4.3「過去 run の署名を現在の内蔵表から再生成しない」)。
            # 比較は「照合不能」になる。内容で照合するにはソルバの解決済み記録 (signature_from_record) を使う。
            species[n] = {"MW": info["MW"][n], "Tlo": None, "Tmid": None, "Thi": None, "LJ_sigma": None, "LJ_eps_kB": None,
                          "nasa9_low": None, "nasa9_high": None, "source": "builtin"}
    return {"run_dir": run_dir, "thermalMethod": info["thermalMethod"], "names": list(info["names"]),
            "MW": [info["MW"][n] for n in info["names"]], "thermoHrefTemp": info["thermoHrefTemp"],
            "tracer": info["tracer"], "species": species, "speciesDBFile": info["speciesDBFile"]}


# ---- 解決済み記録 (C++ input/speciesDB.{hpp,cpp} と同じ正規化; 変えるときは両方を同時に変えてスキーマ版を上げる) ----
SPECIES_RECORD_SCHEMA = "forge_resolved_species_v1"
SPECIES_RECORD_EXTRAPOLATION = "nasa9_2interval; low if T<Tmid; cp clamped at Tlo/Thi; h linear with end cp outside [Tlo,Thi]"
SPECIES_RECORD_DATUM = ("coefficients are absolute (before datum); runtime adds -h_abs(Tref)/Ru to a7 of every interval "
                        "when thermoHrefTemp>0")


def _g17(x):
    """C の printf("%.17g") と同じ文字列 (double は往復でビット一致)。"""
    return "%.17g" % float(x)


def compat_text(schema, datum, thermoHrefTemp, extrapolation, species):
    """互換性ハッシュの正規化テキスト。species は [{name, phase, MW, Tlo, Tmid, Thi, LJ_sigma, LJ_eps_kB, nasa9_low, nasa9_high}]。
    source・来歴は入れない。C++ speciesDB.cpp compatTextRaw と一字一句同じにすること。"""
    out = [f"schema: {schema}", f"datum: {datum}", f"thermoHrefTemp: {_g17(thermoHrefTemp)}",
           f"extrapolation: {extrapolation}", f"nSpecies: {len(species)}"]
    for i, e in enumerate(species):
        out.append(f"species[{i}]: name={e['name']} phase={e['phase']}")
        out.append(f"species[{i}].MW: {_g17(e['MW'])}")
        out.append(f"species[{i}].T: {_g17(e['Tlo'])} {_g17(e['Tmid'])} {_g17(e['Thi'])}")
        out.append(f"species[{i}].LJ: {_g17(e['LJ_sigma'])} {_g17(e['LJ_eps_kB'])}")
        out.append(f"species[{i}].low: " + " ".join(_g17(x) for x in e["nasa9_low"]))
        out.append(f"species[{i}].high: " + " ".join(_g17(x) for x in e["nasa9_high"]))
    return "\n".join(out) + "\n"


def load_record(path):
    """解決済み記録を読み、自己整合を検証して返す。
    返り値: {path, integrity (全文 SHA-256), compat_hash (記録内の値), compat_recomputed (中身から再計算), consistent (bool),
             problems (list), schema, datum, extrapolation, thermoHrefTemp, species [..., source], provenance}。
    ファイル名の互換 16 桁と中身の互換ハッシュが違う・記録内の値と再計算が違う、は problems に入る (取り違え・改竄)。"""
    with open(path, "rb") as f:
        raw = f.read()
    rec = yaml.load(raw.decode("utf-8"), Loader=_StrSafeLoader) or {}
    species = []
    for n in rec.get("species") or []:
        species.append({"name": str(n["name"]), "phase": str(n["phase"]), "source": str(n.get("source", "")),
                        "MW": float(n["MW"]), "Tlo": float(n["Tlo"]), "Tmid": float(n["Tmid"]), "Thi": float(n["Thi"]),
                        "LJ_sigma": float(n["LJ_sigma"]), "LJ_eps_kB": float(n["LJ_eps_kB"]),
                        "nasa9_low": [float(x) for x in n["nasa9_low"]], "nasa9_high": [float(x) for x in n["nasa9_high"]]})
    out = {"path": os.path.abspath(path), "integrity": hashlib.sha256(raw).hexdigest(),
           "compat_hash": str(rec.get("compat_hash", "")), "schema": str(rec.get("schema", "")),
           "datum": str(rec.get("datum", "")), "extrapolation": str(rec.get("extrapolation", "")),
           "thermoHrefTemp": float(rec.get("thermoHrefTemp", 0.0)), "species": species,
           "provenance": rec.get("provenance") or {}}
    out["compat_recomputed"] = hashlib.sha256(
        compat_text(out["schema"], out["datum"], out["thermoHrefTemp"], out["extrapolation"], species).encode()).hexdigest()
    problems = []
    if out["compat_recomputed"] != out["compat_hash"]:
        problems.append(f"{path}: compat_hash in file {out['compat_hash'][:16]} != recomputed from content {out['compat_recomputed'][:16]}"
                        " (edited or corrupted record)")
    m = re.match(r"resolved_species_([0-9a-f]{16})(?:_([0-9a-f]{16}))?\.yaml$", os.path.basename(path))
    if m and m.group(1) != out["compat_recomputed"][:16]:
        problems.append(f"{path}: file name compat {m.group(1)} != content compat {out['compat_recomputed'][:16]} (mixed-up record)")
    if m and m.group(2) and m.group(2) != out["integrity"][:16]:
        problems.append(f"{path}: file name integrity {m.group(2)} != sha256 of file {out['integrity'][:16]}")
    if out["schema"] != SPECIES_RECORD_SCHEMA:
        problems.append(f"{path}: schema {out['schema']!r} (this tool knows {SPECIES_RECORD_SCHEMA!r})")
    out["problems"] = problems
    out["consistent"] = not problems
    return out


def field_species_attrs(h5path):
    """h5 のルート属性 (species_hash, species_record_sha256, species_record_file, species_input_unverified)。無いものは None。"""
    import h5py

    def _s(v):
        return v.decode() if isinstance(v, bytes) else (None if v is None else str(v))
    with h5py.File(h5path, "r") as f:
        a = f.attrs
        unv = a.get("species_input_unverified")
        return {"species_hash": _s(a.get("species_hash")), "species_record_sha256": _s(a.get("species_record_sha256")),
                "species_record_file": _s(a.get("species_record_file")),
                "species_input_unverified": None if unv is None else int(unv)}


def find_record(h5path, dirs=None):
    """h5 の属性に対応する記録を探す。完全性ハッシュが属性と一致するものだけを返す。
    返り値: (record dict | None, 説明文字列)。属性が無ければ (None, 'unverifiable: ...')。"""
    at = field_species_attrs(h5path)
    if not at["species_hash"]:
        return None, f"unverifiable: {h5path} has no species_hash attribute"
    cands = []
    for d in [os.path.dirname(os.path.abspath(h5path))] + list(dirs or []):
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.startswith("resolved_species_" + at["species_hash"][:16]) and fn.endswith(".yaml"):
                cands.append(os.path.join(d, fn))
    for p in cands:
        rec = load_record(p)
        if rec["integrity"] == at["species_record_sha256"]:
            if rec["compat_recomputed"] != at["species_hash"]:
                return None, f"record {p} matches the integrity hash but its content compat {rec['compat_recomputed'][:16]} != attribute {at['species_hash'][:16]}"
            return rec, f"record {p}"
    if cands:
        return None, (f"integrity mismatch: record(s) {cands} exist but none has sha256 {str(at['species_record_sha256'])[:16]}"
                      " (mixed-up or edited record); coefficient differences cannot be identified")
    return None, f"no record resolved_species_{at['species_hash'][:16]}*.yaml found; coefficient differences cannot be identified"


def signature_from_record(rec):
    """記録から compare_signatures 用の署名を作る (全種の係数込み; tracer は記録に無いので比較しない)。"""
    return {"run_dir": os.path.dirname(rec["path"]), "thermalMethod": 2, "names": [e["name"] for e in rec["species"]],
            "MW": [e["MW"] for e in rec["species"]], "thermoHrefTemp": rec["thermoHrefTemp"],
            "species": {e["name"]: dict(e) for e in rec["species"]}, "compat_hash": rec["compat_recomputed"],
            "speciesDBFile": (rec["provenance"] or {}).get("speciesDBFile")}


def required_conserved(sig):
    """署名から restart に必須の保存量データセット名 (VALUE/ 以下) を返す。"""
    req = ["ro", "roUx", "roUy", "roUz", "roe"]
    if sig["thermalMethod"] == 2 and len(sig["names"]) >= 2:
        req += [f"roY{s}" for s in range(len(sig["names"]))]
    if sig["tracer"]:
        req.append("roXi")
    return req


def compare_signatures(a, b, mw_rtol=1e-9, coef_rtol=1e-12):
    """2 つの署名の不一致を説明文字列の list で返す (空 = 一致)。**内容で比較する** (plan §4.3):
    種名・順序、MW、NASA-9 両温度域の全係数 (rel 1e-12)、Tlo/Tmid/Thi、LJ、thermoHrefTemp、tracer (両方が持つとき)。
    source (builtin/file) の違いでは拒否しない (同一係数なら一致)。どちらかの係数が不明 (記録の無い内蔵種) なら
    「照合不能」を返す (以前は両側 builtin を一致とみなして省略していた = 内蔵係数の変更を検出できなかった)。"""
    bad = []
    if a["thermalMethod"] != b["thermalMethod"]:
        bad.append(f"thermalMethod {a['thermalMethod']} vs {b['thermalMethod']}")
    if a["names"] != b["names"]:
        bad.append(f"physProp.species {a['names']} vs {b['names']}")
    else:
        for n in a["names"]:
            sa, sb = a["species"][n], b["species"][n]
            if abs(sa["MW"] - sb["MW"]) > mw_rtol * max(abs(sa["MW"]), abs(sb["MW"]), 1e-300):
                bad.append(f"MW[{n}] {sa['MW']!r} vs {sb['MW']!r}")
            if sa.get("nasa9_low") is None or sb.get("nasa9_low") is None:
                side = " / ".join(t for t, x in (("A", sa), ("B", sb)) if x.get("nasa9_low") is None)
                bad.append(f"{n}: unverifiable (照合不能) — coefficients unknown on side {side} "
                           "(built-in species without a solver record resolved_species_*.yaml)")
                continue
            for k in ("Tlo", "Tmid", "Thi", "LJ_sigma", "LJ_eps_kB"):
                if sa.get(k) is None or sb.get(k) is None:
                    continue
                if sa[k] != sb[k]:
                    bad.append(f"{n}.{k} {sa[k]!r} vs {sb[k]!r}")
            for k in ("nasa9_low", "nasa9_high"):
                for i, (x, y) in enumerate(zip(sa[k], sb[k])):
                    if abs(x - y) > coef_rtol * max(abs(x), abs(y), 1e-300):
                        bad.append(f"{n}.{k}[{i}] {x!r} vs {y!r}")
    if a["thermoHrefTemp"] != b["thermoHrefTemp"]:
        bad.append(f"thermoHrefTemp {a['thermoHrefTemp']} vs {b['thermoHrefTemp']}")
    if "tracer" in a and "tracer" in b and (a["tracer"] or None) != (b["tracer"] or None):
        bad.append(f"physProp.tracer {a['tracer']} vs {b['tracer']}")
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    info = species_info(a.run_dir)
    if a.json:
        out = {k: v for k, v in info.items() if k != "meta"}
        print(json.dumps(out, indent=2))
        return
    print(f"run: {info['run_dir']}  thermalMethod={info['thermalMethod']}")
    for n in info["names"]:
        tag = "  <- condensing" if n == info["condensing"] else ""
        print(f"  {info['index'][n]:2d}  {n:10s} MW={info['MW'][n]:.10g}{tag}")
    print(f"  condensing: {info['condensing']} (index {info['condensing_index']}); H2O: {info['h2o_array']}; vapor_array: {info['vapor_array']}")
    print(f"  tracer: {info['tracer']}")


if __name__ == "__main__":
    main()
