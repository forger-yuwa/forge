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

  熱物性の読み出し (#8/#9; species_db.yaml を前提にしない):
  run_thermo(run_dir, res_path)   # 種ごとの係数・MW・区間・LJ・datum・lump 構成。res の属性が指す記録 > run dir の記録 >
                                  # 従来の speciesDBFile > forge --resolve-species の順 (source="record"|"speciesDBFile" で限定)
  thermo_gas(th, names)           # その物性の凍結組成混合 _TPGas (total_quantities.py と同じ範囲外規約)
  lump_mass_expansion(th, name)   # lump の lump 内質量分率 {構成種: y}

  入口での属性の扱い (#3b; 属性はその保存量を生成した処理が付ける):
  resolve_species(run_dir)     # `forge --resolve-species` (FORGE_BIN / --forge; GPU 不要) で宛先の互換性ハッシュと記録を得る
  stamp_new_field(h5, run_dir, names, MW, h_ref_T, mixtures)   # 新規初期場: IC の datum・順序・MW・e(T) が記録と一致したら付与
  plan_inherit(src_h5, dst_run_dir) / commit_inherit(dst_h5, plan)   # restart・補間: SRC の記録を検証し宛先と一致なら継承
  plan_convert(src_h5, dst_run_dir)   # 種変換: 入力を検証し、変換後に変換先のハッシュを付ける
  未検証の SRC (属性なし / species_input_unverified=1)・宛先を解決できない (旧バイナリ・solverConfig.yaml なし) は
  **既定で停止** (ソルバと同じ規約, #3c)。許可はその実行だけの FORGE_ALLOW_UNVERIFIED_SPECIES=1 か --force-species で、
  許可して通したときは DST に属性を付けない (宛先のハッシュで埋めない; ソルバ側で未検証として扱われる)。

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
    names, lumps = [], {}
    for s in (pp.get("species") or []):
        if isinstance(s, dict):
            # lump (擬似種, plan #6a): {name, lump: {構成種: 分率}, basis: mole|mass}。係数はソルバが起動時に合成する。
            names.append(str(s["name"]))
            lumps[str(s["name"])] = {"basis": str(s.get("basis", "")),
                                     "members": {str(k): float(v) for k, v in (s.get("lump") or {}).items()}}
        else:
            names.append(str(s))
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
    def _mw(n):
        e = _find_ci(db, n) if db else None
        if e is not None and "MW" in e:
            return float(e["MW"])
        mw = _find_ci(BUILTIN_MW, n)
        if mw is None:
            raise KeyError(f"species {n}: MW が species_db.yaml にも内蔵表にも無い")
        return float(mw)

    MW = {}
    for n in names:
        if n in lumps:
            # lump の MW = Σ x_k M_k (ソルバの合成と同じ式; x は lump 内で正規化したモル分率)。正本はソルバの記録 (resolved_species_*.yaml)。
            lp = lumps[n]
            mws = {k: _mw(k) for k in lp["members"]}
            tot = sum(lp["members"].values())
            if lp["basis"] == "mass":
                nmol = {k: v / tot / mws[k] for k, v in lp["members"].items()}
                den = sum(nmol.values())
                x = {k: v / den for k, v in nmol.items()}
            else:
                x = {k: v / tot for k, v in lp["members"].items()}
            MW[n] = sum(x[k] * mws[k] for k in x)
        else:
            MW[n] = _mw(n)
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
        "lumps": lumps,                     # lump 名 → {basis, members: {構成種: config の分率}} (plan #6a)
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
            # lump (plan #6a) も同じ: 合成係数はソルバの記録にだけある。
            species[n] = {"MW": info["MW"][n], "Tlo": None, "Tmid": None, "Thi": None, "LJ_sigma": None, "LJ_eps_kB": None,
                          "nasa9_low": None, "nasa9_high": None, "source": "lump" if n in info["lumps"] else "builtin"}
    return {"run_dir": run_dir, "thermalMethod": info["thermalMethod"], "names": list(info["names"]),
            "MW": [info["MW"][n] for n in info["names"]], "thermoHrefTemp": info["thermoHrefTemp"],
            "tracer": info["tracer"], "species": species, "speciesDBFile": info["speciesDBFile"]}


# ---- 解決済み記録 (C++ input/speciesDB.{hpp,cpp} と同じ正規化; 変えるときは両方を同時に変えてスキーマ版を上げる) ----
SPECIES_RECORD_SCHEMA = "forge_resolved_species_v1"
# physProp.transport を書いた run の記録 (輸送ブロック transport_compat を持つ; plan thermophysics-solver-owned-species-db #5t2)
SPECIES_RECORD_SCHEMA_TRANSPORT = "forge_resolved_species_v2"
SPECIES_RECORD_EXTRAPOLATION = "nasa9_2interval; low if T<Tmid; cp clamped at Tlo/Thi; h linear with end cp outside [Tlo,Thi]"
SPECIES_RECORD_DATUM = ("coefficients are absolute (before datum); runtime adds -h_abs(Tref)/Ru to a7 of every interval "
                        "when thermoHrefTemp>0")
# 凝縮種の液相 (気液ペア; plan thermophysics-solver-owned-species-db #10)。C++ speciesDB.hpp SPECIES_CONDENSED_* と一字一句同じ。
SPECIES_CONDENSED_EXTENSION = ("h_l: NASA-9 on [Tlo,Thi]; T<Tlo: h(Tlo)-cp_l*(Tlo-T) with cp_l=(h(Tlo+0.5)-h(Tlo-0.5+1e-9))/1; "
                               "T>Thi: h(Thi); mass basis = paired gas MW")
SPECIES_CONDENSED_LATENT = ("L=h_v-h_l clamped to [1.5e6,3.5e6] J/kg; h_v = paired gas species (same coefficients, intervals, "
                            "extrapolation and datum as the species DB)")
SPECIES_CONDENSED_DATUM = ("liquid h gets the same constant Ru*da7/MW as the paired gas (da7=-h_abs,gas(Tref)/Ru) when thermoHrefTemp>0 "
                           "(coefficients recorded before datum)")


def _g17(x):
    """C の printf("%.17g") と同じ文字列 (double は往復でビット一致)。"""
    return "%.17g" % float(x)


def condensed_compat_lines(c):
    """凝縮種の液相 (記録の condensed[0]) の互換性行。C++ speciesDB.cpp condensedCompatLines と一字一句同じ。c が None なら空。"""
    if not c:
        return []
    t = "condensed[0]"
    return [f"{t}: name={c['name']} phase=condensed pair_of={c['pair_of']} gas_index={int(c['gas_index'])}",
            f"{t}.MW: {_g17(c['MW'])}",
            f"{t}.T: {_g17(c['Tlo'])} {_g17(c['Thi'])}",
            f"{t}.coeffs: " + " ".join(_g17(x) for x in c["nasa9"]),
            f"{t}.extension: below={c['extension']['below']} above={c['extension']['above']}",
            f"{t}.rule: {c['rule']}", f"{t}.latent: {c['latent']}", f"{t}.datum: {c['datum']}"]


def compat_text(schema, datum, thermoHrefTemp, extrapolation, species, transport_lines=None, condensed=None):
    """互換性ハッシュの正規化テキスト。species は [{name, phase, MW, Tlo, Tmid, Thi, LJ_sigma, LJ_eps_kB, nasa9_low, nasa9_high}]。
    source・来歴は入れない。C++ speciesDB.cpp compatTextRaw と一字一句同じにすること。
    transport_lines: 輸送ブロック (記録の transport_compat; physProp.transport を書いた run だけ, #5t2)。C++ が作った行をそのまま末尾に足す。
    condensed: 凝縮種の液相 (記録の condensed[0]; 凝縮 ON・H2O の run だけ, #10)。種の後・輸送の前に condensed_compat_lines を足す。"""
    out = [f"schema: {schema}", f"datum: {datum}", f"thermoHrefTemp: {_g17(thermoHrefTemp)}",
           f"extrapolation: {extrapolation}", f"nSpecies: {len(species)}"]
    def _coeff_lines(tag, e):
        out.append(f"{tag}.MW: {_g17(e['MW'])}")
        out.append(f"{tag}.T: {_g17(e['Tlo'])} {_g17(e['Tmid'])} {_g17(e['Thi'])}")
        out.append(f"{tag}.LJ: {_g17(e['LJ_sigma'])} {_g17(e['LJ_eps_kB'])}")
        out.append(f"{tag}.low: " + " ".join(_g17(x) for x in e["nasa9_low"]))
        out.append(f"{tag}.high: " + " ".join(_g17(x) for x in e["nasa9_high"]))

    for i, e in enumerate(species):
        tag = f"species[{i}]"
        out.append(f"{tag}: name={e['name']} phase={e['phase']}")
        _coeff_lines(tag, e)
        # lump (plan #6a) だけ追記する: 合成規約・構成種の名前・lump 内モル分率 x・構成種の係数。basis・入力の分率・source は入れない。
        lump = e.get("lump")
        if lump:
            out.append(f"{tag}.lump: n={len(lump['members'])} synthesis={lump['synthesis']}")
            for k, m in enumerate(lump["members"]):
                mt = f"{tag}.lump[{k}]"
                out.append(f"{mt}: name={m['name']} x={_g17(m['x'])}")
                _coeff_lines(mt, m)
    out.extend(condensed_compat_lines(condensed))
    out.extend(transport_lines or [])
    return "\n".join(out) + "\n"


def _record_condensed(n):
    """記録の condensed[0] (凝縮種の液相; #10)。"""
    ex = n.get("extension") or {}
    return {"name": str(n["name"]), "phase": str(n.get("phase", "condensed")), "source": str(n.get("source", "")),
            "pair_of": str(n["pair_of"]), "gas_index": int(n["gas_index"]), "gas_name": str(n.get("gas_name", "")),
            "MW": float(n["MW"]), "Tlo": float(n["Tlo"]), "Thi": float(n["Thi"]),
            "nasa9": [float(x) for x in n["nasa9"]],
            "extension": {"below": str(ex.get("below", "")), "above": str(ex.get("above", ""))},
            "rule": str(n["rule"]), "latent": str(n["latent"]), "datum": str(n["datum"])}


def _record_coeffs(n):
    """記録の 1 種 (または lump 構成種) の係数ブロック。"""
    return {"MW": float(n["MW"]), "Tlo": float(n["Tlo"]), "Tmid": float(n["Tmid"]), "Thi": float(n["Thi"]),
            "LJ_sigma": float(n["LJ_sigma"]), "LJ_eps_kB": float(n["LJ_eps_kB"]),
            "nasa9_low": [float(x) for x in n["nasa9_low"]], "nasa9_high": [float(x) for x in n["nasa9_high"]]}


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
        e = {"name": str(n["name"]), "phase": str(n["phase"]), "source": str(n.get("source", ""))}
        e.update(_record_coeffs(n))
        if n.get("lump"):
            # lump (plan #6a): 合成後の係数 (上) に加えて構成 (basis・入力の分率・lump 内モル分率 x) と構成種の係数
            lp = n["lump"]
            members = []
            for m in lp.get("members") or []:
                mm = {"name": str(m["name"]), "source": str(m.get("source", "")),
                      "fraction_input": float(m["fraction_input"]), "x": float(m["x"])}
                mm.update(_record_coeffs(m))
                members.append(mm)
            e["lump"] = {"basis": str(lp.get("basis", "")), "synthesis": str(lp.get("synthesis", "")), "members": members}
        species.append(e)
    out = {"path": os.path.abspath(path), "integrity": hashlib.sha256(raw).hexdigest(),
           "compat_hash": str(rec.get("compat_hash", "")), "schema": str(rec.get("schema", "")),
           "datum": str(rec.get("datum", "")), "extrapolation": str(rec.get("extrapolation", "")),
           "thermoHrefTemp": float(rec.get("thermoHrefTemp", 0.0)), "species": species,
           "transport_compat": [str(x) for x in (rec.get("transport_compat") or [])],
           "provenance": rec.get("provenance") or {}}
    problems = []
    # 凝縮種の液相 (#10)。無い記録は凝縮 OFF か #10 以前 (潜熱モデルが記録されていない)。
    out["condensed"] = None
    cl = rec.get("condensed")
    if cl:
        try:
            if not isinstance(cl, list) or len(cl) != 1:
                raise ValueError("condensed must be a list of one entry")
            out["condensed"] = _record_condensed(cl[0])
        except (KeyError, TypeError, ValueError) as e:
            problems.append(f"{path}: malformed condensed block ({e})")
    out["compat_recomputed"] = hashlib.sha256(
        compat_text(out["schema"], out["datum"], out["thermoHrefTemp"], out["extrapolation"], species,
                    out["transport_compat"], out["condensed"]).encode()).hexdigest()
    if out["compat_recomputed"] != out["compat_hash"]:
        problems.append(f"{path}: compat_hash in file {out['compat_hash'][:16]} != recomputed from content {out['compat_recomputed'][:16]}"
                        " (edited or corrupted record)")
    m = re.match(r"resolved_species_([0-9a-f]{16})(?:_([0-9a-f]{16}))?\.yaml$", os.path.basename(path))
    if m and m.group(1) != out["compat_recomputed"][:16]:
        problems.append(f"{path}: file name compat {m.group(1)} != content compat {out['compat_recomputed'][:16]} (mixed-up record)")
    if m and m.group(2) and m.group(2) != out["integrity"][:16]:
        problems.append(f"{path}: file name integrity {m.group(2)} != sha256 of file {out['integrity'][:16]}")
    if out["schema"] not in (SPECIES_RECORD_SCHEMA, SPECIES_RECORD_SCHEMA_TRANSPORT):
        problems.append(f"{path}: schema {out['schema']!r} (this tool knows {SPECIES_RECORD_SCHEMA!r}, {SPECIES_RECORD_SCHEMA_TRANSPORT!r})")
    elif (out["schema"] == SPECIES_RECORD_SCHEMA_TRANSPORT) != bool(out["transport_compat"]):
        problems.append(f"{path}: schema {out['schema']!r} and the transport block (transport_compat) do not go together")
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


def _thermo_entry(e):
    """記録の種 / speciesDBFile のエントリ → run_thermo の種エントリ (絶対基準の係数; datum は thermoHrefTemp で別に持つ)。"""
    out = {"MW": float(e["MW"]), "Tlo": float(e.get("Tlo", 200.0)), "Tmid": float(e.get("Tmid", 1000.0)),
           "Thi": float(e.get("Thi", 6000.0)),
           # LJ の既定値は C++ speciesDB_resolve の speciesDBFile 読込と同じ (3.6 Å / 97 K)
           "LJ_sigma": float(e.get("LJ_sigma", 3.6)), "LJ_eps_kB": float(e.get("LJ_eps_kB", 97.0)),
           "nasa9_low": [float(x) for x in e["nasa9_low"]], "nasa9_high": [float(x) for x in e["nasa9_high"]],
           "lump": e.get("lump")}
    if len(out["nasa9_low"]) != 9 or len(out["nasa9_high"]) != 9:
        raise ValueError("nasa9 係数が 9 個でない")
    return out


def _thermo_from_record(rec, how):
    return {"source": "record", "how": how, "path": rec["path"], "names": [e["name"] for e in rec["species"]],
            "thermoHrefTemp": rec["thermoHrefTemp"], "compat_hash": rec["compat_recomputed"],
            "species": {e["name"]: _thermo_entry(e) for e in rec["species"]},
            # 凝縮種の液相 (#10; 凝縮 ON・H2O の記録だけ。無ければ None = 凝縮 OFF か #10 以前で潜熱モデルが記録されていない)
            "condensed": rec.get("condensed")}


def run_thermo(run_dir, res_path=None, source="auto", forge=None):
    """run の熱物性を 1 か所で読む (plans/active/thermophysics-solver-owned-species-db.md §4.6, #8)。後処理・種変換・入口分布・
    設計 runner の署名はこれを使い、`species_db.yaml` の存在を前提にしない。
    返り値: {source: record|speciesDBFile|resolve, how, path, names (physProp.species の順), thermoHrefTemp, compat_hash (記録のときだけ),
             species: {name: {MW, Tlo, Tmid, Thi, LJ_sigma, LJ_eps_kB, nasa9_low, nasa9_high (絶対基準), lump (記録の lump 構成 | None)}},
             condensed: 凝縮種の液相 (記録の condensed[0]; #10) | None (凝縮 OFF・#10 以前の記録・speciesDBFile)}。
    source="auto" の優先順:
      (1) res_path の属性が指すソルバの解決済み記録 (完全性ハッシュを検証; 場を作った物性そのもの)
      (2) run_dir の記録 resolved_species_*.yaml のうち、種名・順序・thermoHrefTemp が solverConfig.yaml と同じもの (互換性ハッシュが 1 通りのとき)
      (3) 従来の speciesDBFile (記録の無い旧 run。config の全種が DB に必要; lump 記法は読めない)
      (4) `forge --resolve-species` (一時ディレクトリで解決; 記録がまだ無い prepare 直後の run)
    source="record" は (1)(2)(4) だけ、source="speciesDBFile" は (3) だけ (旧経路との照合用)。
    CPG (thermalMethod≠2) は None。どれでも解決できなければ ValueError。"""
    run_dir = os.path.abspath(run_dir)
    if source not in ("auto", "record", "speciesDBFile"):
        raise ValueError(f"run_thermo: source {source!r} は未知 (auto | record | speciesDBFile)")
    info = species_info_config(run_dir)
    if info["thermalMethod"] != 2:
        return None
    names = info["names"]
    up = [n.upper() for n in names]
    notes = []
    if source in ("auto", "record"):
        if res_path is not None:
            st = source_species_state(res_path)
            if st["state"] == "verified":
                return _thermo_from_record(st["record"], f"attributes of {os.path.basename(res_path)}")
            if st["state"] == "broken":
                print(f"[forge_species] warning: {res_path}: species record cannot be verified ({st['why']}); "
                      "falling back to the run directory", file=sys.stderr)
            notes.append(f"{os.path.basename(res_path)}: {st['state']}")
        recs = {}
        for fn in sorted(os.listdir(run_dir)):
            if not re.match(r"resolved_species_[0-9a-f]{16}(?:_[0-9a-f]{16})?\.yaml$", fn):
                continue
            rec = load_record(os.path.join(run_dir, fn))
            if not rec["consistent"]:
                continue
            if [e["name"].upper() for e in rec["species"]] != up or rec["thermoHrefTemp"] != info["thermoHrefTemp"]:
                continue
            recs.setdefault(rec["compat_recomputed"], rec)
        if len(recs) == 1:
            return _thermo_from_record(next(iter(recs.values())), "the only matching record in the run directory")
        notes.append(f"{len(recs)} matching record(s) in {run_dir}")
    if source in ("auto", "speciesDBFile") and info["speciesDBFile"]:
        p = info["speciesDBFile"] if os.path.isabs(info["speciesDBFile"]) else os.path.join(run_dir, info["speciesDBFile"])
        db = {str(k): v for k, v in (load_yaml_str(p) or {}).items()}
        missing = [n for n in names if _find_ci(db, n) is None]
        if not missing:
            return {"source": "speciesDBFile", "how": "speciesDBFile (no solver record)", "path": p, "names": list(names),
                    "thermoHrefTemp": info["thermoHrefTemp"], "compat_hash": None,
                    "species": {n: _thermo_entry(_find_ci(db, n)) for n in names},
                    "condensed": None}   # speciesDBFile は液相 (潜熱モデル) を持たない (#10)
        notes.append(f"speciesDBFile {os.path.basename(p)} lacks {missing}")
    elif source == "speciesDBFile":
        raise ValueError(f"{run_dir}: physProp.speciesDBFile が無い (source=speciesDBFile)")
    if source in ("auto", "record"):
        try:
            r = resolve_species(run_dir, forge, inplace=False)
        except (SpeciesResolveUnavailable, SpeciesCheckError) as e:
            notes.append(f"resolve-only unavailable: {e}")
        else:
            if r is not None:
                th = _thermo_from_record(r["record"], "forge --resolve-species (temporary)")
                th["source"] = "resolve"; th["path"] = None
                return th
    raise ValueError(f"{run_dir}: 熱物性を解決できない (species {names}): " + "; ".join(notes))


def species_info_config(run_dir):
    """solverConfig.yaml だけから種の名前・lump・thermalMethod・thermoHrefTemp・speciesDBFile を返す (DB は読まない)。"""
    cfg = load_yaml_str(os.path.join(run_dir, "solverConfig.yaml")) or {}
    pp = cfg.get("physProp") or {}
    tm = int(pp.get("thermalMethod", 0))
    names = [str(s["name"]) if isinstance(s, dict) else str(s) for s in (pp.get("species") or [])]
    if tm == 2 and not names:
        names = ["N2"]
    return {"thermalMethod": tm, "names": names, "thermoHrefTemp": float(pp.get("thermoHrefTemp", 0.0)),
            "speciesDBFile": pp.get("speciesDBFile")}


def thermo_gas(th, names=None):
    """run_thermo の結果から、ソルバと同じ範囲外規約の凍結組成混合 _TPGas (total_quantities.py) を作る。names 省略時は th の順序。
    名前は大小文字を無視して引く (ソルバの解決と同じ)。"""
    from total_quantities import _TPGas
    names = list(names) if names is not None else list(th["names"])
    db = {}
    for n in names:
        e = _find_ci(th["species"], n)
        if e is None:
            raise KeyError(f"species {n} が解決済み熱物性 ({th['source']}: {th['names']}) に無い")
        db[n] = e
    return _TPGas(db, names, th["thermoHrefTemp"])


def lump_mass_expansion(th, name):
    """lump の構成を lump 内質量分率 {構成種: y} で返す (記録の lump 内モル分率 x と構成種 MW から y = x M / Σ x M)。lump でなければ None。"""
    e = _find_ci(th["species"], name)
    lp = (e or {}).get("lump")
    if not lp:
        return None
    den = sum(m["x"] * m["MW"] for m in lp["members"])
    return {m["name"]: m["x"] * m["MW"] / den for m in lp["members"]}


def required_conserved(sig):
    """署名から restart に必須の保存量データセット名 (VALUE/ 以下) を返す。"""
    req = ["ro", "roUx", "roUy", "roUz", "roe"]
    if sig["thermalMethod"] == 2 and len(sig["names"]) >= 2:
        req += [f"roY{s}" for s in range(len(sig["names"]))]
    if sig["tracer"]:
        req.append("roXi")
    return req


def _compare_lumps(n, la, lb, coef_rtol):
    """lump の構成の差 (有無・合成規約・構成種の名前と順序・モル分率 x・構成種の MW/区間/LJ/係数)。記録由来の署名だけが持つ。"""
    if not la and not lb:
        return []
    if bool(la) != bool(lb):
        return [f"{n}: lump on side {'A' if la else 'B'} only"]
    bad = []
    if la["synthesis"] != lb["synthesis"]:
        bad.append(f"{n}.lump.synthesis {la['synthesis']!r} vs {lb['synthesis']!r}")
    na, nb = [m["name"] for m in la["members"]], [m["name"] for m in lb["members"]]
    if na != nb:
        return bad + [f"{n}.lump members {na} vs {nb}"]
    for ma, mb in zip(la["members"], lb["members"]):
        t = f"{n}.lump.{ma['name']}"
        for k in ("x", "MW", "Tlo", "Tmid", "Thi", "LJ_sigma", "LJ_eps_kB"):
            if ma[k] != mb[k]:
                bad.append(f"{t}.{k} {ma[k]!r} vs {mb[k]!r}")
        for k in ("nasa9_low", "nasa9_high"):
            for i, (x, y) in enumerate(zip(ma[k], mb[k])):
                if abs(x - y) > coef_rtol * max(abs(x), abs(y), 1e-300):
                    bad.append(f"{t}.{k}[{i}] {x!r} vs {y!r}")
    return bad


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
            bad += _compare_lumps(n, sa.get("lump"), sb.get("lump"), coef_rtol)
    if a["thermoHrefTemp"] != b["thermoHrefTemp"]:
        bad.append(f"thermoHrefTemp {a['thermoHrefTemp']} vs {b['thermoHrefTemp']}")
    if "tracer" in a and "tracer" in b and (a["tracer"] or None) != (b["tracer"] or None):
        bad.append(f"physProp.tracer {a['tracer']} vs {b['tracer']}")
    return bad


# ---- 入口での属性の付与・継承 (plans/active/thermophysics-solver-owned-species-db.md §4.3「ハッシュは誰が付けるか」, #3b) ----
# 属性は**その保存量を実際に生成した処理**が付ける:
#   新規初期場 (IC 生成)        : stamp_new_field — 宛先を --resolve-species で解決し、IC に使った datum・種順序・MW・
#                                 エネルギー式が記録と一致したときだけ付ける (一致しなければ属性を付けずに SpeciesCheckError)
#   コピー・restart・補間        : plan_inherit — SRC の属性と記録 (完全性ハッシュ再計算) を検証し、宛先の互換性ハッシュと
#                                 一致したときだけ継承。SRC が未検証・宛先を解決できないときは既定で停止し (#3c)、
#                                 許可 (FORGE_ALLOW_UNVERIFIED_SPECIES=1 / --force-species) したときは DST の属性を消す (宛先のハッシュで埋めない)
#   種変換                      : plan_convert — 入力を検証し、変換の成功後に変換先のハッシュを付ける (入力が未検証なら既定で停止、
#                                 許可したときは未検証のまま = 属性なし)
# CPG (thermalMethod≠2) は対象外 (属性を付けない)。

SPECIES_ATTRS = ("species_hash", "species_record_sha256", "species_record_file", "species_input_unverified")
_RESOLVE_FLAG = b"--resolve-species"


class SpeciesCheckError(Exception):
    """化学種の照合で止めるべき状態 (不一致・記録の欠落/改竄・IC の物性が宛先と違う)。"""


class SpeciesResolveUnavailable(Exception):
    """宛先の物性を解決できない (--resolve-species 対応の forge が無い・solverConfig.yaml が無い)。"""


ALLOW_UNVERIFIED_ENV = "FORGE_ALLOW_UNVERIFIED_SPECIES"


def allow_unverified_species():
    """FORGE_ALLOW_UNVERIFIED_SPECIES=1 (その実行だけの許可; ソルバと同じ環境変数, #3c)。
    未検証の SRC・宛先を解決できない状態を停止でなく警告にする (不一致・記録の破損は通さない)。"""
    return os.environ.get(ALLOW_UNVERIFIED_ENV, "") == "1"


UNVERIFIED_GUIDANCE = (
    "  Either:\n"
    "  (1) regenerate the initial field with a generator that stamps species attributes (it resolves the run with "
    "forge --resolve-species; give a new binary by --forge or FORGE_BIN), or copy/restart from a stamped field, or\n"
    "  (2) if you have checked that the field was produced with the same species/datum as the destination, allow it for THIS "
    "invocation only:  FORGE_ALLOW_UNVERIFIED_SPECIES=1 <tool>   (the destination is written WITHOUT species attributes, "
    "so the solver run on it also needs FORGE_ALLOW_UNVERIFIED_SPECIES=1).\n"
    "  If the species set/DB changed, convert the field with tools/convert_species_field.py instead.")


def refuse_unverified(tool, msg, force=False):
    """未検証・解決不能の状態: 既定は SpeciesCheckError (案内 2 通り付き)。FORGE_ALLOW_UNVERIFIED_SPECIES=1 か
    force (--force-species) なら警告して戻る (呼び出し側は属性を付けずに通す)。"""
    if allow_unverified_species() or force:
        how = "FORGE_ALLOW_UNVERIFIED_SPECIES=1" if allow_unverified_species() else "--force-species"
        print(f"[{tool}] WARNING: {msg}. Allowed for this invocation by {how}; destination written WITHOUT species attributes "
              "(the solver refuses it unless FORGE_ALLOW_UNVERIFIED_SPECIES=1 is set for that invocation, "
              "which marks outputs species_input_unverified=1)")
        return
    raise SpeciesCheckError(f"{msg}: UNVERIFIED. Refusing (nothing written).\n" + UNVERIFIED_GUIDANCE)


def _repo_root():
    return os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def forge_supports_resolve(path):
    """forge バイナリが --resolve-species を持つか (旧バイナリは引数を無視して**計算を始めてしまう**ので起動前に確かめる)。"""
    try:
        with open(path, "rb") as f:
            prev = b""
            while True:
                chunk = f.read(1 << 22)
                if not chunk:
                    return False
                if _RESOLVE_FLAG in prev[-len(_RESOLVE_FLAG):] + chunk:
                    return True
                prev = chunk
    except OSError:
        return False


def find_forge(forge=None):
    """--resolve-species を持つ forge を返す: 明示指定 > 環境変数 FORGE_BIN > solver_density_cuda/build/forge。
    明示指定・FORGE_BIN が対応していなければ例外 (黙って別のバイナリに落とさない)。既定のバイナリが旧版なら None。"""
    for label, cand in (("--forge", forge), ("FORGE_BIN", os.environ.get("FORGE_BIN"))):
        if cand:
            cand = os.path.abspath(cand)
            if not os.path.exists(cand):
                raise SpeciesResolveUnavailable(f"{label}={cand} が無い")
            if not forge_supports_resolve(cand):
                raise SpeciesResolveUnavailable(f"{label}={cand} は --resolve-species を持たない旧バイナリ")
            return cand
    default = os.path.join(_repo_root(), "solver_density_cuda", "build", "forge")
    return default if (os.path.exists(default) and forge_supports_resolve(default)) else None


def _config_thermal_method(run_dir):
    p = os.path.join(run_dir, "solverConfig.yaml")
    if not os.path.exists(p):
        return None
    return int(((load_yaml_str(p) or {}).get("physProp") or {}).get("thermalMethod", 0))


def resolve_species(run_dir, forge=None, inplace=True):
    """run_dir の solverConfig.yaml を `forge --resolve-species` で解決する (GPU 不要)。
    inplace=True: run_dir に記録を書く (宛先 run 用。属性の species_record_file がその隣で引けるように)。
    inplace=False: solverConfig.yaml と speciesDBFile を一時ディレクトリへ複製して解決する (SRC 側・dry-run 用; 元 run を書き換えない)。
    返り値: None (CPG) | {hash, record_file, record_path, record (load_record), run_dir}。
    解決できないときは SpeciesResolveUnavailable、解決が失敗したときは SpeciesCheckError。"""
    import shutil, subprocess, tempfile
    run_dir = os.path.abspath(run_dir)
    tm = _config_thermal_method(run_dir)
    if tm is None:
        raise SpeciesResolveUnavailable(f"{run_dir}/solverConfig.yaml が無い (宛先の物性を解決できない)")
    if tm != 2:
        return None
    exe = find_forge(forge)
    if exe is None:
        raise SpeciesResolveUnavailable("--resolve-species を持つ forge が無い (--forge か FORGE_BIN で新しいバイナリを指定する; "
                                        "既定の solver_density_cuda/build/forge は旧版)")
    tmp = None
    cwd = run_dir
    if not inplace:
        tmp = tempfile.mkdtemp(prefix="forge_resolve_")
        cfg = load_yaml_str(os.path.join(run_dir, "solverConfig.yaml")) or {}
        dbf = (cfg.get("physProp") or {}).get("speciesDBFile")
        depth = 0
        if dbf and not os.path.isabs(dbf):
            parts = os.path.normpath(dbf).split(os.sep)
            depth = sum(1 for x in parts if x == "..")
        cwd = os.path.join(tmp, *(["d"] * depth), "run")     # 相対パス ../x.yaml も一時ディレクトリ内で解決させる
        os.makedirs(cwd)
        shutil.copy(os.path.join(run_dir, "solverConfig.yaml"), cwd)
        if dbf:
            src = dbf if os.path.isabs(dbf) else os.path.join(run_dir, dbf)
            if os.path.exists(src) and not os.path.isabs(dbf):
                dst = os.path.normpath(os.path.join(cwd, dbf))
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy(src, dst)
    try:
        p = subprocess.run([exe, "--resolve-species"], cwd=cwd, capture_output=True, text=True)
        if p.returncode == 2 and "thermalMethod != 2" in p.stderr:
            return None
        lines = p.stdout.strip().splitlines()
        h = lines[-1].strip() if lines else ""
        m = re.search(r"\[species\] record (\S+) \(sha256 ([0-9a-f]{64})\)", p.stderr)
        if p.returncode != 0 or not re.fullmatch(r"[0-9a-f]{64}", h) or not m:
            raise SpeciesCheckError(f"forge --resolve-species failed in {run_dir} (rc={p.returncode}): {p.stderr.strip()[-800:]}")
        rec = load_record(os.path.join(cwd, m.group(1)))
        if not rec["consistent"] or rec["compat_recomputed"] != h or rec["integrity"] != m.group(2):
            raise SpeciesCheckError(f"resolve-only record for {run_dir} is inconsistent: {rec['problems']}")
        out = {"hash": h, "record_file": m.group(1), "record": rec, "run_dir": run_dir,
               "record_path": os.path.join(cwd, m.group(1)) if inplace else None, "forge": exe}
        return out
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)


def write_species_attrs(h5path, attrs):
    """h5 のルート属性 4 つを消し、attrs (dict) があれば書く。attrs=None は「未検証のまま」(ソルバが警告/停止で扱う)。
    h5path は開いた h5py.File (書き込み可) でもよい。"""
    import contextlib
    import h5py
    cm = contextlib.nullcontext(h5path) if isinstance(h5path, h5py.File) else h5py.File(h5path, "r+")
    with cm as f:
        for k in SPECIES_ATTRS:
            if k in f.attrs:
                del f.attrs[k]
        if attrs:
            for k in SPECIES_ATTRS[:3]:
                f.attrs[k] = str(attrs[k])
            f.attrs.create("species_input_unverified", int(attrs.get("species_input_unverified", 0)), dtype="int32")


def source_species_state(h5path, dirs=None):
    """保存場の属性と記録の状態: {state: none|unverified|verified|broken, attrs, record, why}。
    none = 属性なし (旧場・CPG・未認証の初期場)、unverified = species_input_unverified=1、
    broken = 属性はあるが記録が無い/完全性ハッシュが合わない (照合不能: 取り違え・改竄を区別できない)。"""
    at = field_species_attrs(h5path)
    if not at["species_hash"]:
        return {"state": "none", "attrs": at, "record": None, "why": f"{h5path} has no species_hash attribute"}
    if at["species_input_unverified"] == 1:
        return {"state": "unverified", "attrs": at, "record": None,
                "why": f"{h5path} descends from an unverified start (species_input_unverified=1)"}
    rec, why = find_record(h5path, dirs)
    if rec is None or not rec["consistent"]:
        return {"state": "broken", "attrs": at, "record": rec,
                "why": (why if rec is None else "; ".join(rec["problems"]))}
    return {"state": "verified", "attrs": at, "record": rec, "why": why}


def place_record(rec, dst_dir):
    """SRC の記録を宛先ディレクトリへ複製する (上書きしない)。同名で中身が違えば `<互換16>_<完全性16>` の名前にする。返り値 = ファイル名。"""
    import shutil
    base = os.path.basename(rec["path"])
    alt = f"resolved_species_{rec['compat_recomputed'][:16]}_{rec['integrity'][:16]}.yaml"
    for name in (base, alt):
        p = os.path.join(dst_dir, name)
        if os.path.abspath(p) == os.path.abspath(rec["path"]):
            return name
        if os.path.exists(p):
            with open(p, "rb") as f:
                if hashlib.sha256(f.read()).hexdigest() == rec["integrity"]:
                    return name
            continue
        shutil.copyfile(rec["path"], p)
        return name
    raise SpeciesCheckError(f"cannot place record {rec['path']} in {dst_dir}: {base} and {alt} exist with different content")


def _record_diff(src_rec, dst_rec):
    """2 つの記録の差 (種・係数・datum)。"""
    a, b = signature_from_record(src_rec), signature_from_record(dst_rec)
    bad = compare_signatures(a, b)
    # 輸送ブロック (physProp.transport; #5t2): 行単位で最初の差を示す
    ta, tb = src_rec.get("transport_compat") or [], dst_rec.get("transport_compat") or []
    if ta != tb:
        k = 0
        while k < min(len(ta), len(tb)) and ta[k] == tb[k]:
            k += 1
        bad.append(f"transport: SRC {ta[k] if k < len(ta) else '(none)'!r} vs destination {tb[k] if k < len(tb) else '(none)'!r}")
    # 凝縮種の液相 (#10): 正規化行で最初の差を示す (液相だけ違う restart も差として出す)
    ca, cb = condensed_compat_lines(src_rec.get("condensed")), condensed_compat_lines(dst_rec.get("condensed"))
    if ca != cb:
        k = 0
        while k < min(len(ca), len(cb)) and ca[k] == cb[k]:
            k += 1
        bad.append(f"condensed (liquid phase / latent heat): SRC {ca[k] if k < len(ca) else '(none)'!r} vs destination "
                   f"{cb[k] if k < len(cb) else '(none)'!r}")
    return bad


def _check_unverified_src(tool, src_h5, dst_run_dir, st, force=False):
    """SRC が未検証 (属性なし / species_input_unverified=1) のとき: 宛先が TP (thermalMethod 2) か、solverConfig.yaml が無く
    CPG か判定できないなら既定で停止 (refuse_unverified)。宛先が CPG なら属性の対象外なので黙って通す。"""
    tm = _config_thermal_method(dst_run_dir) if dst_run_dir else None
    if tm is not None and tm != 2:
        return
    where = (f"destination {dst_run_dir} is thermally perfect (thermalMethod 2)" if tm == 2 else
             f"destination {dst_run_dir} has no solverConfig.yaml (cannot tell CPG from TP; give --dst-run)")
    refuse_unverified(tool, f"SRC {src_h5} is unverified ({st['why']}) and {where}", force)


def plan_inherit(src_h5, dst_run_dir, forge=None, force=False, tool="restart", inplace=True):
    """コピー・restart・補間の継承判定 (書き込み前に呼ぶ)。返り値 = 書き込み後に DST へ付ける属性 (dict) か None (未検証のまま)。
    不一致・記録の欠落は SpeciesCheckError (force=True なら警告して None = 属性を付けずに通す)。
    SRC が未検証で宛先が TP (または判定不能)・宛先を解決できないときも SpeciesCheckError
    (FORGE_ALLOW_UNVERIFIED_SPECIES=1 か force=True なら警告して None)。"""
    dst_run_dir = os.path.abspath(dst_run_dir)
    st = source_species_state(src_h5)
    if st["state"] in ("none", "unverified"):
        _check_unverified_src(tool, src_h5, dst_run_dir, st, force)
        return None
    if st["state"] == "broken":
        msg = (f"SRC {src_h5} carries species_hash {st['attrs']['species_hash'][:16]} but its record cannot be verified: {st['why']}. "
               "UNVERIFIABLE (照合不能): the coefficient differences cannot be identified. Keep the record resolved_species_*.yaml "
               "next to the field, or use --force-species (copies without species attributes)")
        if force:
            print(f"[{tool}] WARNING (--force-species): {msg}")
            return None
        raise SpeciesCheckError(msg)
    rec = st["record"]
    try:
        dst = resolve_species(dst_run_dir, forge, inplace=inplace)
    except SpeciesResolveUnavailable as e:
        refuse_unverified(tool, f"cannot resolve the destination species ({e}); SRC is verified but the destination "
                                "cannot be checked", force)
        return None
    if dst is None:
        msg = f"destination {dst_run_dir} is calorically perfect (thermalMethod != 2) but SRC {src_h5} is a TP field (species_hash {rec['compat_recomputed'][:16]})"
        if force:
            print(f"[{tool}] WARNING (--force-species): {msg}")
            return None
        raise SpeciesCheckError(msg)
    if dst["hash"] != rec["compat_recomputed"]:
        diff = _record_diff(rec, dst["record"])
        msg = (f"species mismatch: SRC {src_h5} species_hash {rec['compat_recomputed'][:16]} (record {os.path.basename(rec['path'])}) "
               f"!= destination {dst_run_dir} species_hash {dst['hash'][:16]} (forge --resolve-species)\n"
               + "".join(f"    {x}\n" for x in (diff or ["(no coefficient difference found; schema/datum text differs)"]))
               + "  Use the matching species DB / thermoHrefTemp, or tools/convert_species_field.py")
        if force:
            print(f"[{tool}] WARNING (--force-species): {msg}\n  -> copying WITHOUT species attributes")
            return None
        raise SpeciesCheckError(msg)
    print(f"[{tool}] species: SRC record {os.path.basename(rec['path'])} verified (integrity {rec['integrity'][:16]}), "
          f"destination species_hash {dst['hash'][:16]} (forge --resolve-species) matches -> inherit")
    return {"species_hash": rec["compat_recomputed"], "species_record_sha256": rec["integrity"],
            "species_record_file": None, "species_input_unverified": 0, "_record": rec, "_dst_dir": dst_run_dir}


def commit_inherit(dst_h5, plan):
    """plan_inherit の結果を書き込み後に DST へ付ける (plan=None なら属性を消す)。記録は DST の隣へ複製する。"""
    if plan is None:
        write_species_attrs(dst_h5, None)
        return None
    path = dst_h5.filename if hasattr(dst_h5, "filename") else dst_h5
    name = place_record(plan["_record"], os.path.dirname(os.path.abspath(path)))
    attrs = {k: plan[k] for k in SPECIES_ATTRS}
    attrs["species_record_file"] = name
    write_species_attrs(dst_h5, attrs)
    return attrs


def record_energy_gas(rec):
    """記録の係数 (絶対基準) と datum でソルバと同じ e(T) を評価する _TPGas (total_quantities.py) を返す。"""
    from total_quantities import _TPGas
    db = {e["name"]: {"MW": e["MW"], "Tlo": e["Tlo"], "Tmid": e["Tmid"], "Thi": e["Thi"],
                      "nasa9_low": e["nasa9_low"], "nasa9_high": e["nasa9_high"]} for e in rec["species"]}
    return _TPGas(db, [e["name"] for e in rec["species"]], rec["thermoHrefTemp"])


def check_ic_against_record(rec, h_ref_T, names, MW, mixtures, e_rtol=1e-9):
    """IC 生成に使った熱物性が宛先の記録と一致するか。返り値 = 問題の list (空なら一致)。
    h_ref_T: IC の datum (None/≤0 は絶対基準 0)、names/MW: IC が組成を作るのに使った輸送種の順序と MW、
    mixtures: [(label, Y (輸送種順の質量分率), R_ic, e_fn)] — e_fn(T ndarray) は IC が roe に使った内部エネルギー [J/kg]。"""
    import numpy as np
    probs = []
    href = float(h_ref_T) if (h_ref_T is not None and float(h_ref_T) > 0.0) else 0.0
    if href != rec["thermoHrefTemp"]:
        probs.append(f"datum: IC h_ref_T {href!r} != destination thermoHrefTemp {rec['thermoHrefTemp']!r}")
    rnames = [e["name"] for e in rec["species"]]
    if [str(n).upper() for n in names] != [n.upper() for n in rnames]:   # ソルバは現状大小文字を同一視 (speciesDB.cpp; #4 で canonical ID 化)
        probs.append(f"species order: IC {list(names)} != destination {rnames}")
        return probs
    for n, m, e in zip(names, MW, rec["species"]):
        if abs(float(m) - e["MW"]) > 1e-12 * e["MW"]:
            probs.append(f"MW[{n}]: IC {float(m)!r} != destination {e['MW']!r}")
    if probs:
        return probs
    gas = record_energy_gas(rec)
    Tlo = max(e["Tlo"] for e in rec["species"]); Thi = min(e["Thi"] for e in rec["species"])
    T = np.array([t for t in (210.0, 300.0, 600.0, 999.0, 1001.0, 1500.0, 3000.0, 5000.0) if Tlo <= t <= Thi])
    for label, Y, R_ic, e_fn in mixtures:
        Y = [float(y) for y in Y]
        if len(Y) != len(names):
            probs.append(f"{label}: composition length {len(Y)} != species {len(names)}")
            continue
        R_rec = gas.Rmix(Y)
        if abs(float(R_ic) - R_rec) > 1e-10 * R_rec:
            probs.append(f"{label}: gas constant IC {float(R_ic)!r} != destination {R_rec!r}")
        e_rec = gas.h(Y, T) - R_rec * T
        e_ic = np.asarray(e_fn(T), dtype=float).reshape(-1)
        tol = e_rtol * (np.abs(e_rec) + R_rec * T + 1.0e3)
        bad = np.abs(e_ic - e_rec) > tol
        if bad.any():
            i = int(np.argmax(np.abs(e_ic - e_rec)))
            probs.append(f"{label}: internal energy IC {e_ic[i]:.9e} != destination {e_rec[i]:.9e} J/kg at T={T[i]:g} K "
                         f"(diff {e_ic[i] - e_rec[i]:.3e})")
    return probs


def stamp_new_field(h5path, run_dir, names, MW, h_ref_T, mixtures, forge=None, tool="IC"):
    """新規初期場の属性付与 (IC 生成の直後に呼ぶ)。宛先 run_dir を --resolve-species で解決し (記録は run_dir に書かれる)、
    check_ic_against_record が空のときだけ属性 (species_input_unverified=0) を付ける。一致しなければ属性を消して SpeciesCheckError。
    CPG は属性なし。解決できない (旧バイナリ・solverConfig.yaml なし) ときは属性を消して SpeciesCheckError
    (FORGE_ALLOW_UNVERIFIED_SPECIES=1 なら警告して属性なしで "unverified")。返り値 = 状態文字列。"""
    run_dir = os.path.abspath(run_dir)
    try:
        dst = resolve_species(run_dir, forge, inplace=True)
    except SpeciesResolveUnavailable as e:
        write_species_attrs(h5path, None)
        refuse_unverified(tool, f"cannot stamp the initial field {h5path}: {e}")
        return "unverified"
    if dst is None:
        write_species_attrs(h5path, None)
        return "cpg"
    probs = check_ic_against_record(dst["record"], h_ref_T, names, MW, mixtures)
    if probs:
        write_species_attrs(h5path, None)
        raise SpeciesCheckError(f"[{tool}] the initial field was generated with thermophysics that differ from the destination "
                                f"{run_dir} (species_hash {dst['hash'][:16]}); no species attributes attached:\n"
                                + "".join(f"    {x}\n" for x in probs))
    write_species_attrs(h5path, {"species_hash": dst["hash"], "species_record_sha256": dst["record"]["integrity"],
                                 "species_record_file": dst["record_file"], "species_input_unverified": 0})
    print(f"[{tool}] species: initial field stamped with species_hash {dst['hash'][:16]} (record {dst['record_file']}; "
          f"datum {dst['record']['thermoHrefTemp']:g} K, species {list(names)})")
    return "verified"


def plan_convert(src_h5, dst_run_dir, src_run_dir=None, forge=None, force=False, tool="convert", inplace=True, legacy_latent=False):
    """種変換の判定 (書き込み前)。入力を検証し (属性・記録の完全性、SRC config を解決したハッシュ = 場の属性)、
    変換先を解決する。返り値 = {"attrs": 付ける属性 | None, "dst": resolve 結果 | None}。入力が未検証・解決できないときは
    既定で SpeciesCheckError、FORGE_ALLOW_UNVERIFIED_SPECIES=1 か force=True なら attrs=None (属性なし)。
    legacy_latent (convert_species_field.py --src-latent legacy-v0; plan #10 の移行手順): SRC の記録が液相 (condensed) を持たない
    #10 以前の H2O 凝縮 run で、SRC config を今の forge で解決した差が**液相の追加だけ** (記録 + 今の液相行 = 今のハッシュ) なら、
    潜熱モデル以外は検証済みとして通す (潜熱は呼び手が旧モデルで読む)。"""
    st = source_species_state(src_h5)
    if st["state"] in ("none", "unverified"):
        _check_unverified_src(tool, src_h5, os.path.abspath(dst_run_dir), st, force)
        return {"attrs": None, "dst": None}
    if st["state"] == "broken":
        msg = f"SRC {src_h5}: species record cannot be verified: {st['why']} (UNVERIFIABLE)"
        if force:
            print(f"[{tool}] WARNING (--force-species): {msg}")
            return {"attrs": None, "dst": None}
        raise SpeciesCheckError(msg)
    rec = st["record"]
    src_run_dir = os.path.abspath(src_run_dir or os.path.dirname(os.path.abspath(src_h5)))
    try:
        srcr = resolve_species(src_run_dir, forge, inplace=False)
        dst = resolve_species(dst_run_dir, forge, inplace=inplace)
    except SpeciesResolveUnavailable as e:
        refuse_unverified(tool, f"cannot resolve species ({e}); the converted field cannot be verified", force)
        return {"attrs": None, "dst": None}
    if srcr is not None and srcr["hash"] != rec["compat_recomputed"] and legacy_latent and rec.get("condensed") is None \
            and (srcr["record"] or {}).get("condensed") is not None:
        h2 = hashlib.sha256(compat_text(rec["schema"], rec["datum"], rec["thermoHrefTemp"], rec["extrapolation"], rec["species"],
                                        rec["transport_compat"], srcr["record"]["condensed"]).encode()).hexdigest()
        if h2 == srcr["hash"]:
            print(f"[{tool}] species: SRC record {rec['compat_recomputed'][:16]} verified except the latent-heat model "
                  f"(record written before plan #10 has no liquid phase; the source config now resolves to {srcr['hash'][:16]} "
                  "only by adding it). The source latent heat is read with the legacy model (--src-latent legacy-v0).")
            rec = dict(rec, compat_recomputed=srcr["hash"])
    if srcr is None or srcr["hash"] != rec["compat_recomputed"]:
        diff = _record_diff(rec, srcr["record"]) if srcr else ["source run is calorically perfect"]
        if srcr is not None and rec.get("condensed") is None and (srcr["record"] or {}).get("condensed") is not None:
            diff.append("the SRC record has no liquid phase (written before plan #10): if the field was produced by a forge "
                        "before #10, convert it with --src-latent legacy-v0 (migration; the old h2o_latent is used for the source)")
        msg = (f"SRC {src_h5} species_hash {rec['compat_recomputed'][:16]} != the source run config {src_run_dir} "
               f"({srcr['hash'][:16] if srcr else 'CPG'}); the converter would read the field with other properties:\n"
               + "".join(f"    {x}\n" for x in diff))
        if force:
            print(f"[{tool}] WARNING (--force-species): {msg}")
            return {"attrs": None, "dst": None}
        raise SpeciesCheckError(msg)
    if dst is None:
        return {"attrs": None, "dst": None}
    print(f"[{tool}] species: SRC record verified ({rec['compat_recomputed'][:16]} = source config); destination resolved {dst['hash'][:16]}")
    return {"attrs": {"species_hash": dst["hash"], "species_record_sha256": dst["record"]["integrity"],
                      "species_record_file": dst["record_file"], "species_input_unverified": 0}, "dst": dst}


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
