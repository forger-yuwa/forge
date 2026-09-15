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
  species_signature(run_dir)   # restart 照合用: config + 解決済み DB (名前/順序, MW, NASA-9 両域係数, Tlo/Tmid/Thi, datum, tracer)。
                               # species_meta.yaml の species が config と矛盾すれば例外
  compare_signatures(a, b)     # 不一致の説明 list (空なら一致; 係数は rel 1e-12)
  required_conserved(sig)      # restart に必須の VALUE/ データセット名

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
            species[n] = {"MW": float(e["MW"]), "Tlo": float(e.get("Tlo", 200.0)), "Tmid": float(e.get("Tmid", 1000.0)),
                          "Thi": float(e.get("Thi", 6000.0)), "nasa9_low": lo, "nasa9_high": hi, "source": "file"}
        else:
            species[n] = {"MW": info["MW"][n], "Tlo": None, "Tmid": None, "Thi": None, "nasa9_low": None, "nasa9_high": None, "source": "builtin"}
    return {"run_dir": run_dir, "thermalMethod": info["thermalMethod"], "names": list(info["names"]),
            "MW": [info["MW"][n] for n in info["names"]], "thermoHrefTemp": info["thermoHrefTemp"],
            "tracer": info["tracer"], "species": species, "speciesDBFile": info["speciesDBFile"]}


def required_conserved(sig):
    """署名から restart に必須の保存量データセット名 (VALUE/ 以下) を返す。"""
    req = ["ro", "roUx", "roUy", "roUz", "roe"]
    if sig["thermalMethod"] == 2 and len(sig["names"]) >= 2:
        req += [f"roY{s}" for s in range(len(sig["names"]))]
    if sig["tracer"]:
        req.append("roXi")
    return req


def compare_signatures(a, b, mw_rtol=1e-9, coef_rtol=1e-12):
    """2 つの署名の不一致を説明文字列の list で返す (空 = 一致)。種名・順序、MW、NASA-9 両温度域の全係数 (rel 1e-12)、
    Tlo/Tmid/Thi、thermoHrefTemp、tracer を比較する。"""
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
            if sa["source"] != sb["source"]:
                bad.append(f"{n}: DB source {sa['source']} vs {sb['source']}")
                continue
            if sa["source"] == "builtin":
                continue
            for k in ("Tlo", "Tmid", "Thi"):
                if sa[k] != sb[k]:
                    bad.append(f"{n}.{k} {sa[k]} vs {sb[k]}")
            for k in ("nasa9_low", "nasa9_high"):
                for i, (x, y) in enumerate(zip(sa[k], sb[k])):
                    if abs(x - y) > coef_rtol * max(abs(x), abs(y), 1e-300):
                        bad.append(f"{n}.{k}[{i}] {x!r} vs {y!r}")
                        break
    if a["thermoHrefTemp"] != b["thermoHrefTemp"]:
        bad.append(f"thermoHrefTemp {a['thermoHrefTemp']} vs {b['thermoHrefTemp']}")
    if (a["tracer"] or None) != (b["tracer"] or None):
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
