r"""組成の単一ソース — モル/質量分率の換算、解決済み種 DB、統一 `tp_species` スキーマ (full | lumped) の解決。

plan `plans/active/thermophysics-cea-mole-fraction-species.md` §4.1–4.5 の実体。

- **ResolvedSpeciesDB**: 内蔵 `SPECIES_NASA9` (CEA thermo.inp 転記) に `gas.species_db` (CEA 直読み DB,
  `cea_thermo_to_species_db.py` の出力) を上書きした 1 つの DB。名前・MW・2 温度域係数・温度区切り・LJ・原子組成・出典を持ち、
  換算・MOC 熱力学・擬似種生成・IC・`species_db.yaml` 出力の**すべて**がこれを使う (codex M1)。
- **mole_to_mass / mass_to_mole**: $Y_k = X_k M_k/\sum_j X_j M_j$。SERN の `frozen.mole_to_mass` はここへ委譲。
- **SpeciesLayout / resolve_species_layout**: `evaluate.tp_species: {mode, lumps, keep}` (旧 `pseudo` / `split_h2o` / `[EXH, AIR]` は
  別名変換) を、流れ (ノズルは 1 流れ、SERN は排気/外気) ごとの質量配分で輸送種順序・入口ベクトル・lump の展開行列に解決する。
  未配分・二重配分・空 lump・名前衝突・凝縮種が keep に無い等は入力段階で拒否 (codex M2/M3)。
- **species_db_yaml / species_meta**: forge `speciesDBFile` (由来コメント付き) と機械可読メタ (`species_meta.yaml`, codex M5)。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .semiperfect import LJ_PARAMS, RU, SPECIES_NASA9, T_MID

# 内蔵 11 種の原子組成 (CEA thermo.inp の元素欄と同じ)。外部 DB は `atoms` キー (cea_thermo_to_species_db.py が書く) を使う。
BUILTIN_ATOMS = {
    "N2": {"N": 2}, "O2": {"O": 2}, "CO2": {"C": 1, "O": 2}, "H2O": {"H": 2, "O": 1}, "AR": {"AR": 1},
    "H2": {"H": 2}, "OH": {"O": 1, "H": 1}, "H": {"H": 1}, "NO": {"N": 1, "O": 1}, "O": {"O": 1}, "CO": {"C": 1, "O": 1},
}
# 元素の原子量 [kg/mol] (元素質量分率の診断用)
ATOMIC_MW = {"H": 1.00794e-3, "C": 12.0107e-3, "N": 14.0067e-3, "O": 15.9994e-3, "AR": 39.948e-3, "HE": 4.002602e-3}

BUILTIN_SOURCE = "CEA thermo.inp (McBride-Gordon 2002), transcribed in forge_design.gas.semiperfect"
STD_RANGES = (200.0, 1000.0, 6000.0)


def _check_name_key(k):
    if isinstance(k, bool):
        raise ValueError("組成のキーに真偽値がある: YAML で NO/ON/OFF 等がクォート無しだと bool に読まれる ('NO': ... と書く)")
    return str(k).upper()


@dataclass
class SpeciesEntry:
    """1 実種 (または擬似種) の NASA-9 データ。forge `speciesDBFile` の 1 エントリと 1:1。"""
    name: str
    MW: float                       # kg/mol
    low: list                       # a0..a8 (Tlo–Tmid)
    high: list                      # a0..a8 (Tmid–Thi)
    Tlo: float = 200.0
    Tmid: float = 1000.0
    Thi: float = 6000.0
    LJ_sigma: float = 3.621
    LJ_eps_kB: float = 97.53
    atoms: dict = field(default_factory=dict)   # {元素: 個数} (擬似種は構成種の加重平均で実数)
    source: str = BUILTIN_SOURCE
    lump_of: dict | None = None     # 擬似種: {構成種: lump 内モル分率}
    lump_mass: dict | None = None   # 擬似種: {構成種: lump 内質量分率}
    Hf298: float | None = None      # J/mol (参考)

    def to_db_dict(self) -> dict:
        d = {"MW": float(self.MW), "LJ_sigma": float(self.LJ_sigma), "LJ_eps_kB": float(self.LJ_eps_kB),
             "Tlo": float(self.Tlo), "Tmid": float(self.Tmid), "Thi": float(self.Thi),
             "nasa9_low": [float(v) for v in self.low], "nasa9_high": [float(v) for v in self.high]}
        return d

    def standard_ranges(self) -> bool:
        return (abs(self.Tlo - STD_RANGES[0]) < 1e-9 and abs(self.Tmid - STD_RANGES[1]) < 1e-9
                and abs(self.Thi - STD_RANGES[2]) < 1e-9)


class ResolvedSpeciesDB:
    """名前 → SpeciesEntry。`builtin()` は呼び出し時点の `SPECIES_NASA9` を読む (テストの monkeypatch を許す)。"""

    def __init__(self, entries: dict):
        self.entries = dict(entries)

    @classmethod
    def builtin(cls) -> "ResolvedSpeciesDB":
        ents = {}
        for k, sp in SPECIES_NASA9.items():
            lj = LJ_PARAMS.get(k, (3.621, 97.53))
            ents[k] = SpeciesEntry(k, float(sp["MW"]), [float(v) for v in sp["low"]], [float(v) for v in sp["high"]],
                                   LJ_sigma=float(lj[0]), LJ_eps_kB=float(lj[1]), atoms=dict(BUILTIN_ATOMS.get(k, {})),
                                   source=BUILTIN_SOURCE)
        return cls(ents)

    @classmethod
    def from_file(cls, path, base: "ResolvedSpeciesDB | None" = None) -> "ResolvedSpeciesDB":
        """外部 DB (forge `speciesDBFile` 形式の yaml) を base (既定: 内蔵) に**上書き**して返す。
        キーは大文字化 (`Ar`→`AR`)。必須: MW>0, nasa9_low/high 各 9 個, Tlo<Tmid<Thi。"""
        import yaml
        base = base if base is not None else cls.builtin()
        raw = yaml.safe_load(Path(path).read_text())
        if not isinstance(raw, dict):
            raise ValueError(f"species_db {path}: mapping ではない")
        ents = dict(base.entries)
        for k, e in raw.items():
            name = _check_name_key(k)
            if not isinstance(e, dict):
                raise ValueError(f"species_db {path}: {name} が mapping ではない")
            try:
                MW = float(e["MW"]); low = [float(v) for v in e["nasa9_low"]]; high = [float(v) for v in e["nasa9_high"]]
            except (KeyError, TypeError, ValueError) as ex:
                raise ValueError(f"species_db {path}: {name} の MW/nasa9_low/nasa9_high が不正 ({ex})")
            if not MW > 0.0:
                raise ValueError(f"species_db {path}: {name} の MW {MW} は正でない")
            if len(low) != 9 or len(high) != 9:
                raise ValueError(f"species_db {path}: {name} の係数は 9 個ずつ必要 (low {len(low)}, high {len(high)})")
            Tlo, Tmid, Thi = (float(e.get("Tlo", 200.0)), float(e.get("Tmid", 1000.0)), float(e.get("Thi", 6000.0)))
            if not (Tlo < Tmid < Thi):
                raise ValueError(f"species_db {path}: {name} の温度区切り {Tlo}/{Tmid}/{Thi} が単調でない")
            lj0 = LJ_PARAMS.get(name, (3.621, 97.53))
            atoms = {str(a).upper(): float(n) for a, n in (e.get("atoms") or {}).items()} or dict(BUILTIN_ATOMS.get(name, {}))
            ents[name] = SpeciesEntry(name, MW, low, high, Tlo, Tmid, Thi,
                                      float(e.get("LJ_sigma", lj0[0])), float(e.get("LJ_eps_kB", lj0[1])),
                                      atoms=atoms, source=f"species_db file {Path(path).name}",
                                      Hf298=(float(e["_Hf298_J_per_mol"]) if "_Hf298_J_per_mol" in e else None))
        return cls(ents)

    def __contains__(self, name) -> bool:
        return str(name).upper() in self.entries

    def __getitem__(self, name) -> SpeciesEntry:
        k = str(name).upper()
        if k not in self.entries:
            raise KeyError(f"species {k} は DB に無い ({sorted(self.entries)})")
        return self.entries[k]

    @property
    def names(self) -> list:
        return list(self.entries)

    def MW(self, name) -> float:
        return self[name].MW

    def require(self, names) -> None:
        missing = [str(n).upper() for n in names if str(n).upper() not in self.entries]
        if missing:
            raise KeyError(f"species {missing} は DB に無い ({sorted(self.entries)})")

    def coef(self, name, T):
        """(low|high) 係数を T で選ぶ (semiperfect._coef 互換)。"""
        e = self[name]
        return np.asarray(e.low if T < e.Tmid else e.high, dtype=float)

    def cp_mass(self, Y: dict, T):
        """混合 cp [J/kg/K] (質量分率 Y, T_FLOOR 凍結は semiperfect と同じ)。"""
        from .semiperfect import _cp_R
        T = np.atleast_1d(np.asarray(T, dtype=float)); out = np.zeros_like(T)
        for k, y in Y.items():
            e = self[k]
            cpR = np.where(T < e.Tmid, _cp_R(np.asarray(e.low), T), _cp_R(np.asarray(e.high), T))
            out += y * cpR * RU / e.MW
        return out

    def h_mass(self, Y: dict, T):
        """混合 h [J/kg] (絶対基準)。"""
        from .semiperfect import _h_RT
        T = np.atleast_1d(np.asarray(T, dtype=float)); out = np.zeros_like(T)
        for k, y in Y.items():
            e = self[k]
            hRT = np.where(T < e.Tmid, _h_RT(np.asarray(e.low), T), _h_RT(np.asarray(e.high), T))
            out += y * hRT * RU * T / e.MW
        return out


# ---------------------------------------------------------------- 換算 (単一ソース)

def normalize_fractions(x: dict, what: str = "組成") -> tuple:
    """キー大文字化・非負検査・正規化。戻り (正規化 dict, 元の総和)。総和 0 は拒否。"""
    x = {_check_name_key(k): float(v) for k, v in x.items()}
    for k, v in x.items():
        if not np.isfinite(v) or v < 0.0:
            raise ValueError(f"{what}: {k} = {v} は負または非有限")
    tot = sum(x.values())
    if not tot > 0.0:
        raise ValueError(f"{what}: 総和が 0")
    return {k: v / tot for k, v in x.items()}, tot


def mole_to_mass(x: dict, db: ResolvedSpeciesDB | None = None) -> dict:
    """モル分率 → 質量分率 $Y_k = X_k M_k/\\sum_j X_j M_j$ (入力は正規化される)。"""
    db = db if db is not None else ResolvedSpeciesDB.builtin()
    x, _ = normalize_fractions(x, "モル分率")
    db.require(x)
    m = {k: v * db.MW(k) for k, v in x.items()}
    mt = sum(m.values())
    return {k: v / mt for k, v in m.items()}


def mass_to_mole(y: dict, db: ResolvedSpeciesDB | None = None) -> dict:
    """質量分率 → モル分率 $X_k = (Y_k/M_k)/\\sum_j (Y_j/M_j)$。"""
    db = db if db is not None else ResolvedSpeciesDB.builtin()
    y, _ = normalize_fractions(y, "質量分率")
    db.require(y)
    n = {k: v / db.MW(k) for k, v in y.items()}
    nt = sum(n.values())
    return {k: v / nt for k, v in n.items()}


def composition_to_mass(species: dict, basis: str, db: ResolvedSpeciesDB | None = None) -> tuple:
    """問題 YAML の `gas.species` + `composition_basis` → (Y 正規化, X 正規化, 元の総和)。"""
    basis = str(basis or "mass").lower()
    if basis not in ("mass", "mole"):
        raise ValueError(f"gas.composition_basis '{basis}' は未知 (mass | mole)")
    db = db if db is not None else ResolvedSpeciesDB.builtin()
    norm, tot = normalize_fractions(species, "gas.species")
    if basis == "mole":
        X = norm; Y = mole_to_mass(norm, db)
    else:
        Y = norm; X = mass_to_mole(norm, db)
    return Y, X, tot


def element_mass_fractions(Y: dict, db: ResolvedSpeciesDB) -> dict:
    """元素質量分率 $Z_e = \\sum_k Y_k\\, n_{e,k} M_e / M_k$ (診断用)。"""
    Z = {}
    for k, y in Y.items():
        e = db[k]
        for el, n in e.atoms.items():
            Z[el] = Z.get(el, 0.0) + y * n * ATOMIC_MW.get(el, 0.0) / e.MW
    return Z


# ---------------------------------------------------------------- 擬似種 (lump)

def lump_entry(name: str, Y_members: dict, db: ResolvedSpeciesDB) -> SpeciesEntry:
    """構成種 (質量分率, lump 内で正規化) を NASA-9 の質量分率線形混合で 1 擬似種にする (厳密)。
    温度区切りが標準 (200/1000/6000) でない種は畳めない (codex M1 の検証条件)。"""
    Y, _ = normalize_fractions(Y_members, f"lump {name}")
    db.require(Y)
    bad = [k for k in Y if not db[k].standard_ranges()]
    if bad:
        raise ValueError(f"lump {name}: 温度区切りが 200/1000/6000 K でない種 {bad} は畳めない (keep に入れる)")
    MW_mix = 1.0 / sum(y / db.MW(k) for k, y in Y.items())
    low = np.zeros(9); high = np.zeros(9); atoms = {}
    for k, y in Y.items():
        e = db[k]; w = y * MW_mix / e.MW
        low += w * np.asarray(e.low, float); high += w * np.asarray(e.high, float)
        for el, n in e.atoms.items():
            atoms[el] = atoms.get(el, 0.0) + w * n
    sig = sum(y * db[k].LJ_sigma for k, y in Y.items()); eps = sum(y * db[k].LJ_eps_kB for k, y in Y.items())
    X = mass_to_mole(Y, db)
    return SpeciesEntry(name, float(MW_mix), [float(v) for v in low], [float(v) for v in high],
                        LJ_sigma=float(sig), LJ_eps_kB=float(eps), atoms=atoms,
                        source=f"lumped from {list(Y)} (mass-fraction linear mixing of NASA-9)",
                        lump_of={k: float(v) for k, v in X.items()}, lump_mass={k: float(v) for k, v in Y.items()})


# ---------------------------------------------------------------- 統一スキーマ

def parse_tp_species(evaluate: dict) -> dict:
    """`evaluate.tp_species` を mapping `{mode, lumps, keep}` に正規化する。旧形式は別名:
    - 省略 / `pseudo` → `{mode: lumped, lumps: {MIX: {from: composition}}, keep: []}`
    - `split_h2o` → `{mode: lumped, lumps: {MIXDRY: {from: composition, exclude: [keep...]}}, keep: [tp_keep_species|H2O]}`
    - `full` → `{mode: full, keep: []}`
    - `[EXH, AIR]` (SERN) → `{mode: lumped, lumps: {EXH: {from: stream, stream: inflow}, AIR: {from: stream, stream: external}}}`
    文字列/リスト形式と `tp_lump` の併用、mapping と `tp_keep_species` の併用は競合として拒否。"""
    ev = evaluate or {}
    ts = ev.get("tp_species", "pseudo")
    keep_old = ev.get("tp_keep_species")
    tl = ev.get("tp_lump")
    if isinstance(ts, dict):
        if tl is not None or keep_old is not None:
            raise ValueError("evaluate.tp_species (mapping) と tp_lump / tp_keep_species は併用不可")
        mode = str(ts.get("mode", "lumped")).lower()
        if mode not in ("full", "lumped"):
            raise ValueError(f"evaluate.tp_species.mode '{mode}' は未知 (full | lumped)")
        lumps = ts.get("lumps") or {}
        keep = [ _check_name_key(k) for k in (ts.get("keep") or []) ]
        if mode == "full" and lumps:
            raise ValueError("evaluate.tp_species.mode: full では lumps は指定不可")
        out = {"mode": mode, "lumps": {}, "keep": keep, "species": [_check_name_key(k) for k in (ts.get("species") or [])]}
        for name, spec in lumps.items():
            lname = _check_name_key(name)
            if not isinstance(spec, dict) or "from" not in spec:
                raise ValueError(f"evaluate.tp_species.lumps.{lname}: {{from: composition|stream, ...}} が必要")
            frm = str(spec["from"]).lower()
            if frm == "composition":
                out["lumps"][lname] = {"from": "composition", "stream": spec.get("stream"),
                                       "exclude": [_check_name_key(k) for k in (spec.get("exclude") or [])]}
            elif frm == "stream":
                if "stream" not in spec:
                    raise ValueError(f"evaluate.tp_species.lumps.{lname}: from: stream には stream: <名> が必要")
                out["lumps"][lname] = {"from": "stream", "stream": str(spec["stream"]), "exclude": []}
            else:
                raise ValueError(f"evaluate.tp_species.lumps.{lname}.from '{frm}' は未知 (composition | stream)")
        return out
    if isinstance(ts, (list, tuple)):
        names = [_check_name_key(k) for k in ts]
        if names != ["EXH", "AIR"]:
            raise ValueError(f"evaluate.tp_species のリスト形式は [EXH, AIR] のみ ({names})")
        if tl is not None or keep_old is not None:
            raise ValueError("evaluate.tp_species: [EXH, AIR] と tp_lump / tp_keep_species は併用不可")
        return {"mode": "lumped", "keep": [],
                "lumps": {"EXH": {"from": "stream", "stream": "inflow", "exclude": []},
                          "AIR": {"from": "stream", "stream": "external", "exclude": []}}}
    ts = str(ts).lower()
    if ts == "pseudo":
        if tl is not None:
            raise ValueError("evaluate.tp_species: pseudo と tp_lump は併用不可")
        return {"mode": "lumped", "lumps": {"MIX": {"from": "composition", "stream": None, "exclude": []}}, "keep": []}
    if ts == "full":
        if tl is not None:
            raise ValueError("evaluate.tp_species: full と tp_lump は併用不可")
        return {"mode": "full", "lumps": {}, "keep": []}
    if ts in ("split_h2o", "lumped"):
        name = "MIXDRY"; keep = [str(keep_old or "H2O").upper()]
        if tl is not None:
            if not isinstance(tl, dict):
                raise ValueError("evaluate.tp_lump は {name, keep} の mapping")
            if keep_old is not None and "keep" in tl:
                raise ValueError("evaluate.tp_lump.keep と tp_keep_species は併用不可")
            name = _check_name_key(tl.get("name", name)); keep = [_check_name_key(k) for k in tl.get("keep", keep)]
        return {"mode": "lumped", "lumps": {name: {"from": "composition", "stream": None, "exclude": list(keep)}}, "keep": keep}
    raise ValueError(f"evaluate.tp_species '{ts}' は未知 (full | lumped | pseudo | split_h2o | [EXH, AIR] | mapping)")


@dataclass
class SpeciesLayout:
    """解決済みの輸送種配置。`species` が forge `physProp.species` の順序 (index = list の位置)。"""
    mode: str
    species: list                     # 輸送種名 (実種 + 擬似種)
    keep: list
    lumps: dict                       # {擬似種名: {"from", "stream", "members": {実種: lump 内モル分率}, "mass": {実種: lump 内質量分率}}}
    streams: dict                     # {流れ名: {"X": {..}, "Y": {..}, "Y_transport": [..] (species 順), "sum_input": float}}
    entries: dict                     # {輸送種名: SpeciesEntry}
    db: ResolvedSpeciesDB
    condensing_species: str | None = None
    tracer: bool = False              # full × 複数流れ: 排気率 ξ を受動スカラ roXi で輸送

    @property
    def n(self) -> int:
        return len(self.species)

    def index(self, name) -> int:
        k = str(name).upper()
        if k not in self.species:
            raise KeyError(f"species {k} は輸送種に無い ({self.species})")
        return self.species.index(k)

    @property
    def cond_index(self) -> int | None:
        """凝縮種の輸送 index。lump に畳まれて輸送種に無い (凝縮 OFF の pseudo 等) なら None。"""
        if self.condensing_species is None or self.condensing_species not in self.species:
            return None
        return self.index(self.condensing_species)

    def Y_transport(self, stream: str) -> list:
        return list(self.streams[stream]["Y_transport"])

    def expansion_matrix(self) -> dict:
        """lump → 実種の展開 (質量分率): {輸送種: {実種: 重み}} (実種は恒等)。"""
        M = {}
        for s in self.species:
            M[s] = dict(self.lumps[s]["mass"]) if s in self.lumps else {s: 1.0}
        return M

    def db_dict(self) -> dict:
        """forge `speciesDBFile` の内容 (species 順)。"""
        return {s: self.entries[s].to_db_dict() for s in self.species}


def resolve_species_layout(tp: dict, streams: dict, db: ResolvedSpeciesDB,
                           condensing_species: str | None = None, condensation: bool = False) -> SpeciesLayout:
    """`tp` = `parse_tp_species()` の結果、`streams` = {流れ名: 質量分率 dict (実種)} (ノズル: {"inflow": Y}、
    SERN: {"inflow": Y_exh, "external": Y_air})。流れごとに keep を先に取り出し、残りをその流れの lump に配分する。"""
    mode = tp["mode"]; keep = list(tp["keep"])
    if not streams:
        raise ValueError("resolve_species_layout: 流れが無い")
    Ys = {}
    for sname, Y in streams.items():
        Yn, _ = normalize_fractions(Y, f"stream {sname}")
        db.require(Yn)
        Ys[sname] = {k: v for k, v in Yn.items() if v > 0.0}
    all_real = []
    for Y in Ys.values():
        for k in Y:
            if k not in all_real:
                all_real.append(k)
    cond = str(condensing_species).upper() if condensing_species else None
    if cond is not None and cond not in all_real:
        raise ValueError(f"凝縮種 {cond} が組成に無い ({all_real})")
    for k in keep:
        if k not in all_real:
            # 作動点によっては keep の種が組成に無い (SERN m4_off に H2O 無し等)。DB にあれば Y=0 の輸送種として残し、
            # 作動点間で輸送種の配置を揃える
            db.require([k])
            all_real.append(k)
    explicit = [ _check_name_key(k) for k in (tp.get("species") or []) ]
    if explicit:
        if mode != "full":
            raise ValueError("tp_species.species (明示の輸送種順序) は mode: full でのみ指定可")
        db.require(explicit)
        missing = [k for k in all_real if k not in explicit]
        if missing:
            raise ValueError(f"tp_species.species に組成の種 {missing} が無い")
        all_real = list(explicit)
    if condensation and cond is None:
        raise ValueError("凝縮 ON なのに gas.condensing_species が無い")

    if mode == "full":
        species = list(all_real)
        entries = {s: db[s] for s in species}
        lumps = {}
        keep_eff = list(species)
    else:
        lumps_in = tp["lumps"]
        if not lumps_in:
            raise ValueError("tp_species.mode: lumped には lumps が必要")
        if cond is not None and condensation and cond not in keep:
            raise ValueError(f"凝縮種 {cond} は tp_species.keep に入れる (lump に畳むと凝縮モデルが指せない)")
        for lname in lumps_in:
            if lname in all_real:
                raise ValueError(f"擬似種名 {lname} が実種名と衝突")
            if lname in keep:
                raise ValueError(f"擬似種名 {lname} が keep と衝突")
        # 流れ→lump の配分
        stream_names = list(Ys)
        members = {lname: {} for lname in lumps_in}   # lump → {実種: 質量 (流れごとの Y を合算した重み)}
        alloc = {}                                      # (流れ, 実種) → lump (二重配分検査)
        for lname, spec in lumps_in.items():
            if spec["from"] == "stream":
                if spec["stream"] not in Ys:
                    raise ValueError(f"lump {lname}: 流れ '{spec['stream']}' が無い ({stream_names})")
                targets = [spec["stream"]]
            else:
                st = spec.get("stream")
                if st is None:
                    if len(stream_names) != 1:
                        raise ValueError(f"lump {lname}: from: composition は流れが 1 つのときだけ省略可 (流れ {stream_names}; stream: を指定)")
                    st = stream_names[0]
                if st not in Ys:
                    raise ValueError(f"lump {lname}: 流れ '{st}' が無い ({stream_names})")
                targets = [st]
            for st in targets:
                for k, y in Ys[st].items():
                    if k in keep or k in spec["exclude"]:
                        continue
                    if (st, k) in alloc:
                        raise ValueError(f"流れ {st} の種 {k} が lump {alloc[(st, k)]} と {lname} に二重配分")
                    alloc[(st, k)] = lname
                    members[lname][k] = members[lname].get(k, 0.0) + y
        for lname, m in members.items():
            if not m:
                raise ValueError(f"lump {lname} が空 (畳む種が無い)")
        # 未配分の検査 (keep でも lump でもない種)
        for st, Y in Ys.items():
            for k in Y:
                if k not in keep and (st, k) not in alloc:
                    raise ValueError(f"流れ {st} の種 {k} が keep にも lump にも配分されていない")
        # 順序: lumps の宣言順 → keep の宣言順 (旧 split_h2o の [MIXDRY, H2O] と互換)
        species = list(lumps_in) + [k for k in keep if k not in lumps_in]
        entries = {}
        lumps = {}
        for lname, m in members.items():
            e = lump_entry(lname, m, db)
            entries[lname] = e
            lumps[lname] = {"from": lumps_in[lname]["from"], "stream": lumps_in[lname]["stream"],
                            "members": dict(e.lump_of), "mass": dict(e.lump_mass)}
        for k in keep:
            entries[k] = db[k]
        keep_eff = keep
    # 流れごとの入口ベクトル
    streams_out = {}
    for st, Y in Ys.items():
        vec = []
        for s in species:
            if s in lumps:
                spec = tp["lumps"][s]
                own = (spec["stream"] == st) if spec["from"] == "stream" else ((spec.get("stream") or st) == st)
                vec.append(sum(y for k, y in Y.items() if own and k in lumps[s]["mass"] and k not in keep_eff) if own else 0.0)
            else:
                vec.append(float(Y.get(s, 0.0)))
        tot = sum(vec)
        if abs(tot - 1.0) > 1e-9:
            raise ValueError(f"流れ {st} の輸送種ベクトルの和 {tot:.12f} != 1 (配分漏れ)")
        streams_out[st] = {"Y": dict(Y), "X": mass_to_mole(Y, db), "Y_transport": [float(v) for v in vec],
                           "sum_input": float(sum(streams[st].values()))}
    tracer = (mode == "full" and len(Ys) >= 2)
    return SpeciesLayout(mode, species, list(keep_eff if mode == "lumped" else []), lumps, streams_out, entries, db,
                         condensing_species=cond, tracer=tracer)


# ---------------------------------------------------------------- 出力

def _fmt(v) -> str:
    return repr(float(v))


def species_db_yaml(layout: SpeciesLayout) -> str:
    """forge `speciesDBFile` のテキスト (species 順, 由来コメント付き)。種名は引用符付き (NO/N/Y の真偽値化を防ぐ)。"""
    out = []
    for s in layout.species:
        e = layout.entries[s]
        out.append(f'"{s}":')
        for k in ("MW", "LJ_sigma", "LJ_eps_kB", "Tlo", "Tmid", "Thi"):
            out.append(f"  {k}: {_fmt(getattr(e, k))}")
        for k, arr in (("nasa9_low", e.low), ("nasa9_high", e.high)):
            out.append(f"  {k}:")
            out += [f"  - {_fmt(v)}" for v in arr]
        if e.atoms:
            out.append("  atoms: {" + ", ".join(f"{a}: {n:g}" for a, n in e.atoms.items()) + "}")
        if e.lump_of:
            out.append("  # lumped: {" + ", ".join(f"{k}: {v:.6g}" for k, v in e.lump_of.items()) + "} (mole fractions within the lump)")
            out.append("  # lumped_mass: {" + ", ".join(f"{k}: {v:.6g}" for k, v in e.lump_mass.items()) + "} (mass fractions within the lump)")
        out.append(f"  # source: {e.source}")
    return "\n".join(out) + "\n"


def species_meta(layout: SpeciesLayout) -> dict:
    """`species_meta.yaml` (run dir) — 後処理・restart が問題 YAML を再解釈せずに使う機械可読メタ。"""
    return {
        "mode": layout.mode,
        "species": list(layout.species),
        "keep": list(layout.keep),
        "condensing_species": layout.condensing_species,
        "condensing_index": layout.cond_index,
        "tracer": {"enabled": bool(layout.tracer), "name": "roXi" if layout.tracer else None,
                   "definition": "exhaust fraction: 1 at stream inflow, 0 at stream external" if layout.tracer else None},
        "lumps": {n: {"from": l["from"], "stream": l["stream"], "mole_fractions": dict(l["members"]), "mass_fractions": dict(l["mass"])}
                  for n, l in layout.lumps.items()},
        "expansion": layout.expansion_matrix(),
        "streams": {st: {"X": dict(v["X"]), "Y": dict(v["Y"]), "Y_transport": list(v["Y_transport"]), "sum_input": v["sum_input"]}
                    for st, v in layout.streams.items()},
        "atoms": {s: dict(layout.entries[s].atoms) for s in layout.species},
        "MW": {s: float(layout.entries[s].MW) for s in layout.species},
        "source": {s: layout.entries[s].source for s in layout.species},
    }


def write_species_files(layout: SpeciesLayout, run_dir) -> None:
    import yaml
    rd = Path(run_dir)
    (rd / "species_db.yaml").write_text(species_db_yaml(layout))
    (rd / "species_meta.yaml").write_text(yaml.safe_dump(species_meta(layout), sort_keys=False, allow_unicode=True))


def load_species_meta(run_dir) -> dict | None:
    import yaml
    p = Path(run_dir) / "species_meta.yaml"
    return yaml.safe_load(p.read_text()) if p.exists() else None
