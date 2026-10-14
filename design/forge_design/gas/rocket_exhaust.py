r"""ロケット排気 — 燃焼室の平衡組成を凍結した TP (plan tooling-nozzle-bell-two-mode §4.1, §5.1 #2a)。

③ベルの評価器 (両モード共通) のガス。親計画 tooling-nozzle-design-tool §3 確定事項 8・10 (全機種 TP 基準、
基準ケース = LOX/LH2・燃焼室圧 10 MPa・O/F 6、CEA の燃焼室平衡組成を凍結した TP)。

1. **排気ファイル** (`resolve_exhaust` → `write_exhaust_file`、schema `forge_rocket_exhaust_v1`):
   NASA CEA2 本体 (`FCEA2`) で燃焼室の平衡組成・全温を求め、反応物の状態 (名前・温度)・燃焼室圧・O/F・
   生成物の候補 (`only`)・凍結位置 (燃焼室、`rocket frozen nfz=1`) と、同じ条件の CEA 凍結膨張の表 (照合の参照値) を
   1 つの YAML に書く。**組成・全温・全圧の単一ソース**で、問題 YAML は `gas.exhaust` でこのファイルを指す。
2. **契約** `RocketExhaust`: 排気ファイル・解決済み種 DB (`ResolvedSpeciesDB`)・エンタルピー基準 (`thermoHrefTemp`) の
   3 つから、MOC 側の `GasSemiPerfect` (`moc_gas`)、CFD 側の輸送種配置 (lump 擬似種 `EXH`; `layout`)、
   IC・理想値の `FrozenGas` (`frozen_gas`) を**同じ入力から**作る。どれも同じ質量分率 $Y$・同じ NASA-9・同じ $T_c$ を見る。
3. **1D 理想値** `ideal_1d`: 凍結 TP の等エントロピー膨張 ($h_0=h(T_c)$、$s(T,p)=s(T_c,p_c)$)。
   スロートは凍結音速 $u=a_f$ ($a_f^2=\gamma_f R T$、$\gamma_f=c_p/(c_p-R)$)、$\dot m_{1D}=\rho^*a^*A^*$。
   出口は超音速の枝で $\rho_e u_e=\rho^*a^*/\varepsilon$、$F_{1D}(\varepsilon)=\dot m u_e+(p_e-p_a)A_e$。

生成物の候補は既定で forge の熱物性 (共通 species データ) と輸送データ (CEA trans.inp) の両方にある 6 種
(H2・O2・H2O・OH・H・O) に限る。HO2・H2O2・O3 は共通データに無いので、CEA でもこれを除いた平衡を解き
(組成と全温が forge の DB と矛盾しない)、全生成物で解いた場合との差を `omitted_species` に記録する。

CLI (排気ファイルを作る):
    python -m forge_design.gas.rocket_exhaust --out case/68.bell_two_mode/exhaust_lox_lh2_pc100_of6.yaml
"""
from __future__ import annotations

import hashlib
import os
import re
import subprocess
import tempfile
from dataclasses import asdict, dataclass, field
from decimal import Decimal
from pathlib import Path

import numpy as np

SCHEMA = "forge_rocket_exhaust_v1"
EXH_NAME = "EXH"
# forge の熱物性 (forge_species_v1.yaml) と輸送 (forge_transport_v1.yaml = CEA trans.inp) の両方にある H/O 系の生成物
DEFAULT_PRODUCTS = ("H2", "O2", "H2O", "OH", "H", "O")
# CEA の出力のケース名 (1 入力ファイルに 3 問題)
_CASE_EQ, _CASE_FZ, _CASE_FULL = "bell_eq", "bell_fz", "bell_full"


# ---------------------------------------------------------------- CEA の所在

def cea_dir() -> Path:
    """FCEA2・thermo.lib・trans.lib の置き場所。環境変数 FORGE_CEA_DIR > <FORGE_ROOT>/.venv-cea/nasa_cea >
    git の主作業ツリーの .venv-cea/nasa_cea (worktree は共通 git dir の親に置かれた CEA を使う)。"""
    root = Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3]))
    cand = []
    env = os.environ.get("FORGE_CEA_DIR")
    if env:
        cand.append(Path(env))
    cand.append(root / ".venv-cea" / "nasa_cea")
    try:
        common = subprocess.run(["git", "-C", str(root), "rev-parse", "--git-common-dir"],
                                capture_output=True, text=True, check=True).stdout.strip()
        cp = Path(common)
        if not cp.is_absolute():
            cp = root / cp
        cand.append(cp.resolve().parent / ".venv-cea" / "nasa_cea")
    except (OSError, subprocess.CalledProcessError):
        pass
    for c in cand:
        if all((c / f).exists() for f in ("FCEA2", "thermo.lib", "trans.lib")):
            return c
    raise FileNotFoundError("NASA CEA2 (FCEA2・thermo.lib・trans.lib) が見つからない。FORGE_CEA_DIR で指定する "
                            f"(探した場所: {[str(c) for c in cand]})")


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------- 入力

@dataclass(frozen=True)
class ChamberSpec:
    """燃焼室の条件と CEA の入力。反応物の温度と凍結位置 (燃焼室) は固定して記録する (plan §4.1)。"""
    fuel: str = "H2(L)"              # CEA の反応物名 (thermo.lib にある凝縮相)
    fuel_T: float = 20.27            # [K] (NBP の LH2; CEA の H2(L) の定義温度)
    oxid: str = "O2(L)"
    oxid_T: float = 90.17            # [K] (NBP の LOX; CEA の O2(L) の定義温度)
    pc: float = 1.0e7                # 燃焼室圧 [Pa]
    of_ratio: float = 6.0            # 酸化剤/燃料 (質量比)
    products: tuple = DEFAULT_PRODUCTS   # 平衡の生成物の候補 (CEA の only)

    def check(self) -> None:
        if not (self.pc > 0.0 and self.of_ratio > 0.0 and self.fuel_T > 0.0 and self.oxid_T > 0.0):
            raise ValueError(f"ChamberSpec: 正でない値がある ({self})")
        if not self.products:
            raise ValueError("ChamberSpec.products が空")


def _reac_block(s: ChamberSpec) -> str:
    return (f"reac\n  fuel = {s.fuel}  wt%=100.  t,k={s.fuel_T!r}\n"
            f"  oxid = {s.oxid}  wt%=100.  t,k={s.oxid_T!r}\n")


def cea_input_text(spec: ChamberSpec, supar=(50.0,)) -> str:
    """CEA の入力 (3 問題)。
    - bell_eq  : rocket equilibrium + only (燃焼室の平衡組成。trace で微量種も 5 桁の指数表記)
    - bell_fz  : rocket frozen nfz=1 + only (燃焼室で凍結した等エントロピー膨張 = 1D 理想値の照合の参照、輸送物性つき)
    - bell_full: rocket equilibrium、生成物の制限なし (除いた種の量と全温の差の診断)"""
    spec.check()
    sup = ",".join(repr(float(v)) for v in supar)
    pbar = repr(spec.pc / 1.0e5)
    head = f"o/f={spec.of_ratio!r} p,bar={pbar}\n   supar={sup}\n"
    only = "only " + " ".join(spec.products) + "\n"
    return (f"problem case={_CASE_EQ} rocket equilibrium {head}{_reac_block(spec)}{only}"
            "output siunits trace=1.e-15\nend\n"
            f"problem case={_CASE_FZ} rocket frozen nfz=1 {head}{_reac_block(spec)}{only}"
            "output siunits transport\nend\n"
            f"problem case={_CASE_FULL} rocket equilibrium {head}{_reac_block(spec)}"
            "output siunits trace=1.e-15\nend\n")


def run_cea(text: str, workdir=None, name: str = "bell") -> str:
    """FCEA2 を一時ディレクトリ (または workdir) で回して .out の本文を返す。"""
    cdir = cea_dir()
    tmp = None
    if workdir is None:
        tmp = tempfile.TemporaryDirectory(prefix="forge_cea_")
        workdir = tmp.name
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    try:
        for f in ("thermo.lib", "trans.lib"):
            dst = wd / f
            if not dst.exists():
                dst.symlink_to(cdir / f)
        (wd / f"{name}.inp").write_text(text)
        r = subprocess.run([str(cdir / "FCEA2")], input=f"{name}\n", cwd=wd, capture_output=True, text=True)
        out = wd / f"{name}.out"
        if r.returncode != 0 or not out.exists():
            raise RuntimeError(f"FCEA2 が失敗 (rc={r.returncode}):\n{r.stdout[-2000:]}{r.stderr[-2000:]}")
        return out.read_text(encoding="latin-1")
    finally:
        if tmp is not None:
            tmp.cleanup()


# ---------------------------------------------------------------- 出力の解析

_MANT_EXP = re.compile(r"^(-?\d+\.\d*)([+-]\d+)$")      # "2.1329-2" / "1.883-12"
_MANT = re.compile(r"^-?\d+\.\d*$")
_EXP = re.compile(r"^[+-]?\d+$")


def cea_numbers(body: str) -> list:
    """CEA の数値欄 ("4.6118 0" / "2.1329-2" / "2.854 -3" / "-9493.24") を float の列にする。
    小数点を持つ仮数の直後の整数だけの語は指数とみなす (この表の値は常に小数点を持つ)。"""
    toks = body.split()
    out = []
    i = 0
    while i < len(toks):
        t = toks[i]
        m = _MANT_EXP.match(t)
        if m:                           # 10 進の文字列のまま指数を付けて読む (掛け算の丸めを入れない)
            out.append(float(f"{m.group(1)}e{int(m.group(2))}"))
            i += 1
        elif _MANT.match(t) and i + 1 < len(toks) and _EXP.match(toks[i + 1]):
            out.append(float(f"{t}e{int(toks[i + 1])}"))
            i += 2
        else:
            out.append(float(t))
            i += 1
    return out


def _case_block(out: str, case: str) -> list:
    """`CASE = <case>` から次の `CASE =` (または末尾) までの行。"""
    lines = out.splitlines()
    starts = [i for i, L in enumerate(lines) if L.strip().startswith("CASE =")]
    for k, i in enumerate(starts):
        if lines[i].split("=", 1)[1].strip() == case:
            j = starts[k + 1] if k + 1 < len(starts) else len(lines)
            return lines[i:j]
    raise ValueError(f"CEA 出力に CASE = {case} が無い")


# 行ラベル → キー (先頭一致。ラベルの後ろが数値欄)
_ROWS = [("Pinf/P", "Pinf_P"), ("P, BAR", "P_bar"), ("T, K", "T"), ("RHO, KG/CU M", "rho"), ("H, KJ/KG", "h_kJ"),
         ("S, KJ/(KG)(K)", "s_kJ"), ("M, (1/n)", "MW_g"), ("Cp, KJ/(KG)(K)", "cp_kJ"), ("GAMMAs", "gamma_s"),
         ("SON VEL,M/SEC", "a"), ("MACH NUMBER", "M"), ("VISC,MILLIPOISE", "visc_mP"), ("CONDUCTIVITY", "cond_mW_cmK"),
         ("Ae/At", "Ae_At"), ("CSTAR, M/SEC", "cstar"), ("CF", "CF"), ("Ivac, M/SEC", "Ivac"), ("Isp, M/SEC", "Isp")]
_PERF = ("Ae_At", "cstar", "CF", "Ivac", "Isp")          # 性能欄は CHAMBER 列が無い


def parse_rocket_table(out: str, case: str) -> dict:
    """rocket 問題の表 {キー: [CHAMBER, THROAT, EXIT...]} (SI: P は Pa、h/s/cp は J/kg 系、visc は Pa s、cond は W/m/K)。
    性能欄 (Ae/At, CSTAR, CF, Ivac, Isp) は CHAMBER 列を None で埋める。同じラベルは最初の出現だけを使う。"""
    blk = _case_block(out, case)
    rows = {}
    for L in blk:
        s = L.strip()
        for lab, key in _ROWS:
            if key in rows or not s.startswith(lab) or (key == "CF" and s.split()[0] != "CF"):
                continue
            try:
                vals = cea_numbers(s[len(lab):])
            except ValueError:          # 見出し行 ("CONDUCTIVITY IN UNITS OF ...") は数値欄でない
                vals = []
            if vals:
                rows[key] = vals
            break
    if "T" not in rows or "P_bar" not in rows:
        raise ValueError(f"CEA 出力 CASE = {case}: 表 (T, P) が読めない")
    n = len(rows["T"])
    for k in _PERF:
        if k in rows:
            rows[k] = [None] * (n - len(rows[k])) + rows[k]
    # 単位の換算は 10 のべきなので、印字の 10 進表現に指数を足して読み直す (掛け算の丸めを入れない)
    conv = {"P_bar": ("P", 5), "h_kJ": ("h", 3), "s_kJ": ("s", 3), "cp_kJ": ("cp", 3),
            "visc_mP": ("visc", -4), "cond_mW_cmK": ("cond", -1)}
    outd = {}
    for k, v in rows.items():
        if k in conv:
            nk, e = conv[k]
            outd[nk] = [None if x is None else float(Decimal(repr(x)).scaleb(e)) for x in v]
        else:
            outd[k] = list(v)
    outd["n_columns"] = n
    return outd


def parse_mole_fractions(out: str, case: str, column: int = 0) -> dict:
    """`MOLE FRACTIONS` 欄 (trace 指定の 1 種 1 行の形) から列 `column` のモル分率 {種: X} (印字の値そのまま)。
    種名の先頭の `*` (20000 K までの熱物性の印) は落とす。"""
    blk = _case_block(out, case)
    try:
        i0 = next(i for i, L in enumerate(blk) if L.strip() == "MOLE FRACTIONS")
    except StopIteration:
        raise ValueError(f"CEA 出力 CASE = {case}: MOLE FRACTIONS が無い")
    X = {}
    for L in blk[i0 + 1:]:
        s = L.strip()
        if not s:
            if X:
                break
            continue
        if s.startswith("*") and "THERMODYNAMIC" in s:
            break
        name, rest = s.split(None, 1)
        try:
            vals = cea_numbers(rest)
        except ValueError:
            raise ValueError(f"CEA 出力 CASE = {case}: モル分率の行が 1 種 1 行の形でない ({s!r}; output に trace を指定する)")
        X[name.lstrip("*")] = vals[column]
    if not X:
        raise ValueError(f"CEA 出力 CASE = {case}: モル分率が読めない")
    return X


def _program_line(out: str) -> str:
    for L in out.splitlines():
        if "CHEMICAL EQUILIBRIUM PROGRAM" in L:
            return " ".join(L.split())
    return "unknown"


def _thermo_update(out: str) -> str:
    m = re.search(r"LAST thermo\.inp UPDATE:\s*(\S+)", out)
    return m.group(1) if m else "unknown"


# ---------------------------------------------------------------- 排気ファイル

def resolve_exhaust(spec: ChamberSpec | None = None, supar=(50.0,), workdir=None) -> dict:
    """CEA を回して排気ファイルの内容 (dict) を作る。組成は `bell_eq` の燃焼室列 (only の生成物、5 桁)、
    全温・全圧は同じ列、参照表は `bell_fz` (燃焼室で凍結)。"""
    spec = spec or ChamberSpec()
    text = cea_input_text(spec, supar)
    out = run_cea(text, workdir)
    eq = parse_rocket_table(out, _CASE_EQ)
    fz = parse_rocket_table(out, _CASE_FZ)
    full = parse_rocket_table(out, _CASE_FULL)
    X = parse_mole_fractions(out, _CASE_EQ, 0)
    X_full = parse_mole_fractions(out, _CASE_FULL, 0)
    prods = [p.upper() for p in spec.products]
    unknown = [k for k in X if k.upper() not in prods]
    if unknown:
        raise ValueError(f"CEA の生成物 {unknown} が only {spec.products} に無い")
    if abs(fz["T"][0] - eq["T"][0]) > 0.011:
        raise ValueError(f"燃焼室温度が平衡 {eq['T'][0]} と凍結 {fz['T'][0]} の問題で違う")
    omitted = {k: v for k, v in X_full.items() if k.upper() not in prods and v > 0.0}
    cdir = cea_dir()
    cols = ["CHAMBER", "THROAT"] + [f"EXIT_{i + 1}" for i in range(fz["n_columns"] - 2)]
    keep = ("Pinf_P", "P", "T", "rho", "h", "s", "MW_g", "cp", "gamma_s", "a", "M", "visc", "cond",
            "Ae_At", "cstar", "CF", "Ivac", "Isp")
    return {
        "schema": SCHEMA,
        "cea_input": {**{k: (list(v) if isinstance(v, tuple) else v) for k, v in asdict(spec).items()},
                      "pc_unit": "Pa", "of_ratio_basis": "mass", "fuel_T_unit": "K", "oxid_T_unit": "K",
                      "freeze": "chamber (CEA rocket frozen nfz=1; 組成は燃焼室の平衡を凍結)",
                      "combustor": "infinite-area (CEA rocket iac 既定)", "supar": [float(v) for v in supar]},
        "cea_provenance": {"program": _program_line(out), "thermo_inp_update": _thermo_update(out),
                           "fcea2": str(cdir / "FCEA2"), "thermo_lib_sha256": _sha256(cdir / "thermo.lib"),
                           "trans_lib_sha256": _sha256(cdir / "trans.lib")},
        "chamber": {"Pc_Pa": float(eq["P"][0]), "Tc_K": float(eq["T"][0]),
                    "mole_fractions": {k: float(v) for k, v in X.items()},
                    "mole_fraction_sum_printed": float(sum(X.values())),
                    "MW_g_per_mol_cea": float(eq["MW_g"][0]),
                    "visc_Pa_s_cea_frozen": float(fz["visc"][0]) if fz.get("visc") else None,
                    "cond_W_mK_cea_frozen": float(fz["cond"][0]) if fz.get("cond") else None},
        "omitted_species": {"products": "all (bell_full)", "mole_fractions": omitted,
                            "sum": float(sum(omitted.values())), "Tc_full_K": float(full["T"][0]),
                            "dTc_K": float(Decimal(repr(full["T"][0])) - Decimal(repr(eq["T"][0])))},
        "cea_frozen_reference": {"columns": cols[:fz["n_columns"]],
                                 **{k: fz[k] for k in keep if k in fz}},
        "cea_input_text": text,
    }


def write_exhaust_file(path, data: dict) -> None:
    import yaml
    head = ("# ロケット排気 (燃焼室の平衡組成を凍結した TP) — forge_design.gas.rocket_exhaust が NASA CEA2 から生成。手で編集しない。\n"
            "# 組成・全温・全圧の単一ソース: 問題 YAML の gas.exhaust がこのファイルを指す (plan tooling-nozzle-bell-two-mode §4.1)。\n")
    Path(path).write_text(head + yaml.safe_dump(data, sort_keys=False, allow_unicode=True, width=200))


def load_exhaust_file(path) -> dict:
    import yaml
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema") != SCHEMA:
        raise ValueError(f"{path}: schema が {SCHEMA} でない")
    ch = data.get("chamber") or {}
    for k in ("Pc_Pa", "Tc_K", "mole_fractions"):
        if k not in ch:
            raise ValueError(f"{path}: chamber.{k} が無い")
    return data


# ---------------------------------------------------------------- 1D 理想値 (凍結 TP の等エントロピー膨張)

def ideal_1d(gas, T0: float, p0: float, eps: float, p_a: float = 0.0, A_t: float | None = None) -> dict:
    r"""燃焼室 (よどみ: $T_0$, $p_0$, $u=0$) から凍結 TP で等エントロピー膨張した 1D の値。`gas` は `FrozenGas`
    (`h_mass` 絶対基準・`entropy(T, p)`・`cp_mass`・`R`)。
    スロート: $2(h_0-h(T^*))=\gamma_f(T^*)RT^*$。出口 (超音速の枝): $\rho_e u_e=\rho^*a^*/\varepsilon$。
    戻り値の `*_per_At` は単位スロート面積あたり、`A_t` [m²] を与えると `mdot`・`F_vac`・`F` も返す。"""
    from scipy.optimize import brentq

    T0 = float(T0); p0 = float(p0); eps = float(eps)
    if not (eps >= 1.0):
        raise ValueError(f"ideal_1d: 面積比 {eps} < 1")
    R = float(gas.R)
    h0 = float(gas.h_mass(T0)[0])
    s_t = float(gas.entropy(T0, p0)[0])
    s0m = lambda T: float(gas.s0_mass(T)[0])
    h = lambda T: float(gas.h_mass(T)[0])
    p_of = lambda T: 1.0e5 * float(np.exp((s0m(T) - s_t) / R))          # s0(T) − R ln(p/1 bar) = s_t
    u_of = lambda T: float(np.sqrt(max(2.0 * (h0 - h(T)), 0.0)))
    g_of = lambda T: float(gas.gamma(T)[0])
    a_of = lambda T: float(np.sqrt(g_of(T) * R * T))
    G_of = lambda T: p_of(T) / (R * T) * u_of(T)
    xt = dict(xtol=1e-13 * T0, rtol=4 * np.finfo(float).eps, maxiter=500)
    T_s = brentq(lambda T: 2.0 * (h0 - h(T)) - g_of(T) * R * T, 0.3 * T0, T0 * (1.0 - 1e-12), **xt)
    p_s = p_of(T_s); rho_s = p_s / (R * T_s); a_s = a_of(T_s); G_s = rho_s * u_of(T_s)
    out = {"T0": T0, "p0": p0, "h0_abs": h0, "R": R, "eps": eps, "p_a": float(p_a),
           "T_star": T_s, "p_star": p_s, "rho_star": rho_s, "a_star": a_s, "u_star": u_of(T_s),
           "mdot_per_At": G_s, "cstar": p0 / G_s}
    if eps == 1.0:
        T_e = T_s
    else:
        target = G_s / eps
        T_lo = 20.0
        if G_of(T_lo) >= target:
            raise ValueError(f"ideal_1d: 面積比 {eps} の出口温度が {T_lo} K より低い")
        T_e = brentq(lambda T: G_of(T) - target, T_lo, T_s, **xt)
    p_e = p_of(T_e); u_e = u_of(T_e); a_e = a_of(T_e)
    F_vac = G_s * u_e + p_e * eps
    out.update({"T_e": T_e, "p_e": p_e, "rho_e": p_e / (R * T_e), "u_e": u_e, "a_e": a_e, "M_e": u_e / a_e,
                "gamma_e": g_of(T_e), "F_vac_per_At": F_vac, "F_per_At": F_vac - float(p_a) * eps,
                "Ivac": F_vac / G_s, "CF_vac": F_vac / p0})
    if A_t is not None:
        A_t = float(A_t)
        out.update({"A_t": A_t, "mdot": G_s * A_t, "F_vac": F_vac * A_t, "F": (F_vac - float(p_a) * eps) * A_t})
    return out


def cea_reference_at(data: dict, eps: float) -> dict:
    """排気ファイルの CEA 凍結表から面積比 `eps` の出口列を取り出し、1D 理想値と同じ量に換算する。
    ṁ/A* = p_c/c*、u_e = Isp (出口圧まで膨張した比推力 = 出口速度)、F_vac/A* = (p_c/c*)·Ivac。"""
    ref = data["cea_frozen_reference"]
    ae = ref["Ae_At"]
    idx = [i for i, v in enumerate(ae) if v is not None and abs(v - eps) <= 1e-6 * eps]
    if not idx:
        raise ValueError(f"CEA 参照表に Ae/At = {eps} の列が無い ({ae})")
    i = idx[0]
    pc = float(data["chamber"]["Pc_Pa"])
    G = pc / float(ref["cstar"][i])
    return {"column": ref["columns"][i], "T_e": float(ref["T"][i]), "p_e": float(ref["P"][i]), "u_e": float(ref["Isp"][i]),
            "M_e": float(ref["M"][i]), "mdot_per_At": G, "F_vac_per_At": G * float(ref["Ivac"][i]),
            "T_star": float(ref["T"][1]), "p_star": float(ref["P"][1]), "cstar": float(ref["cstar"][i]),
            "Ivac": float(ref["Ivac"][i])}


# ---------------------------------------------------------------- 契約: 1 か所で解決した組成・DB・基準から MOC と CFD を作る

@dataclass
class RocketExhaust:
    """排気ファイル + 解決済み種 DB + エンタルピー基準 → MOC (`moc_gas`)・CFD (`layout`)・IC/理想値 (`frozen_gas`)。

    3 つとも同じ質量分率 `Y` (排気ファイルのモル分率を `db` の MW で換算)・同じ `db` の NASA-9・同じ `Tc`/`Pc` を使う。
    `href_T` は forge の `thermoHrefTemp` (IC の内部エネルギーと CFD の顕エンタルピーの基準) で、MOC は比の量しか
    使わないので基準に依らない (`GasSemiPerfect` は絶対基準の h の差だけを使う)。"""
    data: dict
    db: object                         # gas.composition.ResolvedSpeciesDB
    href_T: float = 298.15
    source: str | None = None          # 排気ファイルのパス (記録用)
    name: str = EXH_NAME
    X: dict = field(init=False)
    Y: dict = field(init=False)
    Tc: float = field(init=False)
    Pc: float = field(init=False)

    def __post_init__(self):
        from .composition import mole_to_mass, normalize_fractions
        if self.data.get("schema") != SCHEMA:
            raise ValueError(f"排気データの schema が {SCHEMA} でない")
        ch = self.data["chamber"]
        X, _ = normalize_fractions(ch["mole_fractions"], "gas.exhaust の燃焼室組成")
        self.db.require(X)
        self.X = X
        self.Y = mole_to_mass(X, self.db)
        self.Tc = float(ch["Tc_K"])
        self.Pc = float(ch["Pc_Pa"])
        if not (self.Tc > 0.0 and self.Pc > 0.0):
            raise ValueError(f"gas.exhaust: Tc {self.Tc} / Pc {self.Pc} が正でない")

    @classmethod
    def from_file(cls, path, db=None, href_T: float = 298.15) -> "RocketExhaust":
        from .composition import ResolvedSpeciesDB
        return cls(load_exhaust_file(path), db if db is not None else ResolvedSpeciesDB.builtin(), float(href_T), str(path))

    # --- 3 つの気体 (同じ Y・db) ---
    def frozen_gas(self):
        from .frozen import FrozenGas
        return FrozenGas(self.Y, self.name, self.href_T, self.db)

    def moc_gas(self, **kw):
        from .semiperfect import GasSemiPerfect
        return GasSemiPerfect(self.Y, Tt=self.Tc, db=self.db, **kw)

    @staticmethod
    def tp_species() -> dict:
        """CFD の輸送種配置 (統一スキーマ `evaluate.tp_species`): 単一の lump 擬似種 EXH (plan §4.1)。"""
        return {"mode": "lumped", "lumps": {EXH_NAME: {"from": "composition"}}}

    def layout(self, tp_species: dict | None = None):
        from .composition import parse_tp_species, resolve_species_layout
        tp = parse_tp_species({"tp_species": tp_species if tp_species is not None else self.tp_species()})
        return resolve_species_layout(tp, {"inflow": self.Y}, self.db)

    # --- 1D 理想値 ---
    def ideal_1d(self, eps: float, p_a: float = 0.0, A_t: float | None = None) -> dict:
        return ideal_1d(self.frozen_gas(), self.Tc, self.Pc, eps, p_a=p_a, A_t=A_t)

    def visc_estimate(self) -> float | None:
        """陰解法対角と局所 dt の剛性見積り `physProp.visc` に使う μ (領域の最大 ≈ 燃焼室; CEA の凍結輸送物性)。"""
        v = self.data["chamber"].get("visc_Pa_s_cea_frozen")
        return float(v) if v is not None else None

    def summary(self) -> dict:
        ci = self.data.get("cea_input", {})
        return {"source": self.source, "name": self.name, "Tc_K": self.Tc, "Pc_Pa": self.Pc, "X": self.X, "Y": self.Y,
                "href_T": self.href_T,
                "reactants": {"fuel": ci.get("fuel"), "fuel_T_K": ci.get("fuel_T"),
                              "oxid": ci.get("oxid"), "oxid_T_K": ci.get("oxid_T"), "of_ratio": ci.get("of_ratio")},
                "freeze": ci.get("freeze"), "products": ci.get("products"),
                "omitted_species_sum": (self.data.get("omitted_species") or {}).get("sum"),
                "cea_program": (self.data.get("cea_provenance") or {}).get("program")}


# ---------------------------------------------------------------- CLI

def main(argv=None) -> int:
    import argparse
    import json

    ap = argparse.ArgumentParser(description="CEA2 で燃焼室の平衡組成を求め、凍結 TP の排気ファイルを書く")
    d = ChamberSpec()
    ap.add_argument("--fuel", default=d.fuel)
    ap.add_argument("--fuel-T", type=float, default=d.fuel_T)
    ap.add_argument("--oxid", default=d.oxid)
    ap.add_argument("--oxid-T", type=float, default=d.oxid_T)
    ap.add_argument("--pc", type=float, default=d.pc, help="燃焼室圧 [Pa]")
    ap.add_argument("--of", type=float, default=d.of_ratio, help="O/F (質量比)")
    ap.add_argument("--products", default=" ".join(d.products), help="生成物の候補 (空白区切り; CEA の only)")
    ap.add_argument("--supar", default="50", help="参照表の面積比 (カンマ区切り)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    spec = ChamberSpec(a.fuel, a.fuel_T, a.oxid, a.oxid_T, a.pc, a.of, tuple(a.products.split()))
    data = resolve_exhaust(spec, [float(v) for v in a.supar.split(",")])
    write_exhaust_file(a.out, data)
    from .composition import ResolvedSpeciesDB
    ex = RocketExhaust(data, ResolvedSpeciesDB.builtin(), source=str(a.out))
    print(json.dumps({"out": a.out, **{k: v for k, v in ex.summary().items() if k not in ("Y",)}}, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
