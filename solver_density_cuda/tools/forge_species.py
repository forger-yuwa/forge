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
  info["condensing"]       # 'H2O' | None,  info["condensing_index"]  # 1 | None
  info["tracer"]           # 'exhaust' | None
  info["vapor_array"]      # 'Y1' (凝縮種の Y 配列名; 凝縮種が無ければ None)
"""
import argparse, json, os, sys
import yaml

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
    cfg = yaml.safe_load(open(cfg_path))
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
            db = yaml.safe_load(open(p)) or {}
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
    meta = yaml.safe_load(open(meta_path)) if os.path.exists(meta_path) else None
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
    return {
        "run_dir": run_dir,
        "thermalMethod": tm,
        "names": names,
        "index": index,
        "MW": MW,
        "condensing": condensing,
        "condensing_index": index.get(condensing) if condensing else None,
        "vapor_array": (f"Y{index[condensing]}" if condensing else None),
        "tracer": tracer,
        "condensation": int(cond.get("condensation", 0)) == 1,
        "condModel": int(cond.get("condModel", 0)),
        "thermoHrefTemp": float(pp.get("thermoHrefTemp", 0.0)),
        "speciesDBFile": db_file,
        "meta": meta,
    }


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
    print(f"  condensing: {info['condensing']} (index {info['condensing_index']}, array {info['vapor_array']})")
    print(f"  tracer: {info['tracer']}")


if __name__ == "__main__":
    main()
