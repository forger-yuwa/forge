r"""風洞ノズルの実寸の要求カード `forge.wt_card/v1` → 問題 YAML (`wind_tunnel_axisym_axismach`)。

plan: plans/accepted/tooling-nozzle-wt-deliverables.md §4.1 (1〜6)、検証 §6 V1。

カードは実寸 [m, Pa, K] で書く YAML (未知のキー・欠けたキー・型違い・負値は拒否)。生成器は

1. **長さの写し**: r_t で割って `geometry.r_inlet`・`L_U`・`L_pipe`・`dv.L_total` (+ `Lc_mode: from_length`) にする。
   `geometry.L_total_tol` = カードの長さ公差 / 2 / r_t。`from_length` の固定点反復が改善しなくなって止まったら (離散化のノイズ床とは
   区別できない) `n_axis_inv` を 2 倍にして 1 回やり直し、それでも止まれば内部許容差を 0.05 r_t に**明示して**緩める (カードの公差は
   変えず、寸法の合否はカードの公差で測る; codex 諮問 2026-10-13)。x_F が逆 MOC の場の終点に張り付いたら `x_end_margin` を上げる。
2. **r_t の固定点** (物理スロート径のとき): r_t で YAML を書く → `deltastar_loop.solve_rt_throat` (CFD 前の δ の経路) → r_t を
   更新して書き直す、を $|\Delta r_t|/r_t \le$ rel_tol まで繰り返す (max_outer 回で止まらなければ例外)。設計スロート径のときは 1 回。
   `--r-throat` で外から r_t を与えたとき (NS 後の解き直し) は解かずに書く。
3. **縮流部の検査**: $R_{in} > r_t$、$\mu = L_U^2/(R\,\Delta r) \le 20$、`L_pipe` ≥ 0、`M_exit` > 1。外れたら例外。
4. **壁解像の初期見積もり** (`estimate_mesh`): 1 次元の等エントロピー + Eckert の参照温度 + 平板乱流の $C_f$ で $y_1^+ = 0.8$ の第 1 層を
   断面ごとに見積もり、`mesh.wall_first_frac`・`wall_first_frac_throat`・`ni` を埋める。**初期値であって合格の根拠ではない**。
   壁温分布のカードでは断面ごとの局所の $T_w(x)$ を使う (plan tooling-nozzle-wall-temperature-profile §4.3)。
5. **空気の液化の事前警告** (`liquefaction_check`): 出口の静温・静圧と N₂ の飽和蒸気圧 (Antoine) を比べる (止めない)。
6. **雛形** (`recipes/*.yaml`) の構成キーを写し、カードの `overrides` を最後に当てる (上書きは `_card` 節とコメントに残す)。
7. **knot の M_knot** (plan tooling-nozzle-mknot-auto §4.3・§4.4): 既定は自動 (カードに書かない、または `overrides: {geometry.M_knot: auto}`)。
   最初の r_t で `mknot_select.select_mknot` → `resolution_gate` (規則 min_theta_max_v1) の `chosen_N` を選び、r_t の固定点の間は固定する。
   収束した最終の組で詰めと解像度のゲートを回し直し、`chosen_N` が固定点に使った M_K から 0.01 を超えて動いたら r_t の固定点をやり直す
   (最大 2 回)。ゲートで選び直した値が詰めの窓の端なら、窓を半幅 0.1 → 0.2 → 0.3 に広げて選び直す (上限でも端なら例外; §5.1 #3b)。
   記録は card_log の `M_knot_rule` (掃引の表を含む全部)・YAML の `_card.M_knot_rule` (要約だけ; card_log の記録の sha256 で結ぶ)・
   YAML のコメント。数値の `overrides: {geometry.M_knot: x}` は今と同じ (探索しない)。

`_card` 節は問題 YAML の最上位に置く記録 (ソルバ・設計チェーンは読まない。`probdef.load_problem` は未知の最上位キーを拒否しない)。

usage (design/ 直下):
  .venv-opt/bin/python -m forge_design.intake.wt_card CARD --out problem.yaml [--r-throat R] [--no-solve]
    → problem.yaml と problem.yaml.card_log.json (カードの値 → YAML のキーの対応表・警告・r_t の反復・格子の見積もり・寸法)
  NS 後の補正 (2026-10-13, plan tooling-nozzle-deltastar-band-root-zone §4.3 3):
  .venv-opt/bin/python -m forge_design.intake.wt_card CARD --out problem.yaml --ns-correction PREV_RUN [--delta-next CSV]
    → 上に加えて problem_delta_r.csv (最終の r_t の補正の表; NS は deltastar_loop --delta-r-csv で立てる)
"""
from __future__ import annotations

import argparse
import copy
import datetime
import hashlib
import json
import math
import re
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import yaml

CARD_FORMAT = "forge.wt_card/v1"
GENERATOR = "forge_design.intake.wt_card"
RECIPE_DIR = Path(__file__).resolve().parent / "recipes"
DEFAULT_RECIPE = "axismach_explore_v1"

# gas.preset: air = 乾き空気の質量分率 (plan §4.1)。和は 1 (0.7552 + 0.2314 + 0.0129 + 0.0005)
GAS_PRESETS = {"air": {"N2": 0.7552, "O2": 0.2314, "AR": 0.0129, "CO2": 0.0005}}
SPECIES_SUM_TOL = 1e-9
# 出荷ゲート (§4.2a) の公差の既定。**仮の値** (ユーザの実用公差が決まったら置き換える; plan §4.1)
DEFAULT_TOLERANCES = {"length_m": 1.0e-3, "diameter_m": 2.0e-5}
MU_MAX = 20.0                        # 縮流部の単調収縮の条件 μ = L_U²/(R Δr) ≤ 20 (UpstreamThroatPoly.mu、手順書 §0)
# スロートの壁の曲率 (設計壁の曲率半径) のカードのキー: 実寸 [m] (r_t を解き直すたびに割り直す) か r_t 比 (= geometry.R)
THROAT_CURVATURE_KEYS = ("curvature_radius_m", "curvature_ratio")
R_HALL_MIN = 0.5                     # Hall–KL 遷音速解の精度が落ちる R の下限 (geometry/transonic.py HallThroat の警告と同じ値)

# 壁解像の初期見積もり (§4.1 4.)
Y1PLUS_TARGET = 0.8                  # 第 1 層の y1+ の目標
RECOVERY_FACTOR = 0.89               # 断熱壁の回復係数 (乱流 Pr^(1/3)、Pr 0.71)
X_DOWNSTREAM_MIN = 6.0               # wall_first_frac はスロート下流 x > 6 [r_t] の最小
X_THROAT_HALF = 2.0                  # wall_first_frac_throat は |x| ≤ 2 [r_t] の最小
NI_CAP = 6000                        # ni の上限 (超えたら上限にして警告)
AR_MARGIN = 0.97                     # ni の見積もりは AR ≤ AR_MARGIN × ar_max で行う: 見積もりは設計壁で測るが、メッシュは δ_r を足した
                                     # 物理壁で作るので AR が少し上がる (case/70 M5 デモ: 設計壁で 4981 → 物理壁で 5002.4 > 5000、2026-10-13)
FLOAT32_FRAC_WARN = 1e-6             # 第 1 層 / 局所半径がこれ未満なら float32 の座標の警告 (止めない)
SUTHERLAND_AIR = (1.716e-5, 273.15, 110.4)   # μ_ref [Pa s]・T_ref [K]・S [K] (空気)

# N₂ の飽和蒸気圧 (Antoine: log10(P/bar) = A − B/(T + C)、63〜126 K; NIST WebBook の Edejer–Thodos の係数)
N2_ANTOINE = (3.7362, 264.651, -6.788)
N2_ANTOINE_RANGE_K = (63.0, 126.0)

# YAML 1.2 の 10 進の浮動小数 (PyYAML が文字列と読む `3.0e6` を数値として受けるため)
_YAML12_FLOAT = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

# knot の M_knot の自動の選択 (plan tooling-nozzle-mknot-auto §4.3・§4.4; 規則は evaluate/mknot_select.py の min_theta_max_v1)
MKNOT_AUTO = "auto"                  # overrides の geometry.M_knot に書ける「自動」(書かないのと同じ)
MKNOT_MOVE_TOL = 0.01                # 最終の組のゲートの chosen_N と r_t の固定点に使った M_K の差の上限 (§4.3; ゲートの max_change と同じ値)
MKNOT_REDO_MAX = 2                   # その差が上限を超えたときの r_t の固定点のやり直しの上限 (§4.3)
MKNOT_EVAL_MEM_GB_2400 = 1.3         # evaluate_mknot の 1 評価のピークのメモリ [GB] (case/70 M5・n_axis_inv 2400 で 1.28 GB、4800 で 3.86 GB;
MKNOT_EVAL_MEM_EXP = 1.6             # 2026-10-14 実測)。N の (N/2400)^1.6 で見積もる
# M_knot の評価の並列の数 (N, 2N)。2026-10-14: 手元の PC は他の計算と共用で、N を 2 並列にした生成の最中にメモリ不足で PC が落ちた
# → N・2N とも 1 (取り決め 2026-10-14)
MKNOT_WORKERS = (1, 1)

# カードのキー
_TOP_REQUIRED = ("format", "name", "gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall")
_TOP_OPTIONAL = ("p_ambient", "tolerances", "recipe", "overrides")
# overrides で上書きできる節と、カードの値から決まるので上書きできないキー (カードの該当キーを直す)
_OVERRIDE_SECTIONS = ("gas", "geometry", "mesh", "mesh_euler", "evaluate", "deltastar_initializer")
_OVERRIDE_FORBIDDEN = ("gas.model", "gas.species", "gas.composition_basis",
                       "geometry.r_inlet", "geometry.L_U", "geometry.L_pipe", "geometry.Lc_mode", "geometry.L_total_tol")
# mesh のうち estimate_mesh が埋めるキー (overrides で与えたらそちらを使う)
_MESH_ESTIMATED = ("ni", "wall_first_frac", "wall_first_frac_throat")
# 等温壁の壁温の与え方: 一様・表・CSV の厳密に 1 つ (plan tooling-nozzle-wall-temperature-profile §4.1・§4.3)
WALL_TW_KEYS = ("Tw", "Tw_profile", "Tw_profile_csv")
CSV_HEADERS = (("s", "Tw_K"), ("x_m", "Tw_K"))    # Tw_profile_csv のヘッダ行 (カンマ区切り)
_COORD_TEXT = {"s": "s = 区間で正規化した座標 (入口直管 [−2,−1]・縮流部 [−1,0]・拡大部 [0,1]、スロート 0・出口 1)",
               "x_m": "x_m = 実寸 [m]、原点は設計スロート"}
_PROBLEM_ORDER = ("name", "type", "gas", "deltastar_initializer", "spec", "dv", "geometry", "mesh", "mesh_euler", "evaluate", "_card")


# ---------------------------------------------------------------------------------------------------- カードの読み込みと検査

def _num(errs: list, where: str, v, gt: float | None = None, ge: float | None = None):
    """有限の数値 (bool は不可) を float で返す。不正なら errs に足して None。

    PyYAML (YAML 1.1) は指数に符号の無い `3.0e6` を**文字列**と読む (リポジトリの問題 YAML も `Pt: 5.5e6` と書き、probdef は
    float() で読んでいる)。YAML 1.2 では数値なので、その形の文字列 (_YAML12_FLOAT) だけは数値として受ける。他の文字列は型違い。"""
    if isinstance(v, str) and _YAML12_FLOAT.fullmatch(v.strip()):
        v = float(v)
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(float(v)):
        errs.append(f"{where} は有限の数値 (受け取った値: {v!r})")
        return None
    v = float(v)
    if gt is not None and not v > gt:
        errs.append(f"{where} は {gt:g} より大きい値 (受け取った値: {v!r})")
        return None
    if ge is not None and not v >= ge:
        errs.append(f"{where} は {ge:g} 以上 (受け取った値: {v!r})")
        return None
    return v


def _mapping(errs: list, where: str, d, required: tuple, optional: tuple = ()) -> bool:
    """d が mapping で、未知のキー・欠けたキーが無いか。"""
    if not isinstance(d, dict):
        errs.append(f"{where} は mapping (受け取った値: {d!r})")
        return False
    unknown = [str(k) for k in d if k not in required + optional]
    missing = [k for k in required if k not in d]
    if unknown:
        errs.append(f"{where} の未知のキー {unknown} (書けるのは {list(required + optional)})")
    if missing:
        errs.append(f"{where} のキー {missing} が無い")
    return not unknown and not missing


def _flatten_overrides(errs: list, ov) -> dict:
    """overrides を {"節.キー[.下位]": 値} にする。最上位が節名 (点を含まない) で値が mapping なら 1 段だけ展開する
    ({geometry: {R: 2.5}} → geometry.R)。点を含むキーはそのまま ({"geometry.wall_fit_mono_r2": [0, 1.5]})。"""
    if ov is None:
        return {}
    if not isinstance(ov, dict):
        errs.append(f"overrides は mapping (受け取った値: {ov!r})")
        return {}
    flat = {}
    for k, v in ov.items():
        if not isinstance(k, str) or not k:
            errs.append(f"overrides のキーは文字列 (受け取った値: {k!r})")
            continue
        items = ([(f"{k}.{kk}", vv) for kk, vv in v.items()] if ("." not in k and isinstance(v, dict)) else [(k, v)])
        for path, val in items:
            parts = path.split(".")
            if len(parts) < 2 or any(not p for p in parts):
                errs.append(f"overrides.{path}: 「節.キー」の形で書く (例 geometry.R)")
                continue
            if parts[0] not in _OVERRIDE_SECTIONS:
                errs.append(f"overrides.{path}: 節 {parts[0]!r} は上書きできない (書けるのは {list(_OVERRIDE_SECTIONS)}; "
                            "spec・dv・name はカードの値から決まる)")
                continue
            if any(path == f or path.startswith(f + ".") for f in _OVERRIDE_FORBIDDEN):
                errs.append(f"overrides.{path}: カードの値から決まるキーは上書きできない (カードの該当キーを直す)")
                continue
            if path in flat:
                errs.append(f"overrides.{path} が 2 回指定されている")
                continue
            flat[path] = copy.deepcopy(val)
    return flat


def _yaml12_floats(obj):
    """壁温の表 (mapping の値の列) の要素のうち、PyYAML が文字列と読む YAML 1.2 の浮動小数 (`5.0e2` など) だけを数値にする
    (カードの他の数値と同じ扱い; `_num`)。他の値はそのまま返し、検査は `wall_temperature.parse_profile` に任せる。"""
    if not isinstance(obj, dict):
        return obj
    return {k: ([float(e) if isinstance(e, str) and _YAML12_FLOAT.fullmatch(e.strip()) else e for e in v]
                if isinstance(v, list) else v) for k, v in obj.items()}


def read_profile_csv(path_given, card_dir=None) -> tuple:
    """カードの `wall.Tw_profile_csv` を読む → (表 `{s | x_m: [...], Tw_K: [...]}`, 出典 `{path, abs_path, sha256}`)。

    カンマ区切り、最初の空でない行がヘッダ `s,Tw_K` か `x_m,Tw_K` (各欄の前後の空白・空行・UTF-8 の BOM は無視)、以降は 1 行 2 列の数値。
    **相対パスはカードのファイルの所在 (card_dir) が基準** (cwd ではない; codex plan レビュー m5)。card_dir が無い (カードを dict で
    渡した) ときの相対パスは拒否。表の検査 (単調・正・有限・範囲) は呼ぶ側が `wall_temperature.parse_profile` で行う。"""
    if not isinstance(path_given, str) or not path_given.strip():
        raise ValueError(f"CSV のパス (文字列) を書く (受け取った値: {path_given!r})")
    p = Path(path_given)
    if not p.is_absolute():
        if card_dir is None:
            raise ValueError(f"相対パス {path_given!r} はカードのファイルの所在が基準 — カードをファイルから読む (load_card) か、"
                             "絶対パスで書く")
        p = Path(card_dir) / p
    p = p.resolve()
    if not p.is_file():
        raise ValueError(f"CSV {str(p)!r} が無い (カードの値 {path_given!r}; 相対パスはカードのファイルの所在が基準)")
    data = p.read_bytes()
    rows = [(i + 1, ln.strip()) for i, ln in enumerate(data.decode("utf-8-sig").splitlines()) if ln.strip()]
    if not rows:
        raise ValueError(f"CSV {str(p)!r} が空")
    hdr = tuple(c.strip() for c in rows[0][1].split(","))
    if hdr not in CSV_HEADERS:
        raise ValueError(f"CSV {str(p)!r} の列名 (ヘッダ行) は {' か '.join(','.join(h) for h in CSV_HEADERS)} "
                         f"(受け取った値: {rows[0][1]!r}; カンマ区切り)")
    xs, ts = [], []
    for ln_no, ln in rows[1:]:
        cells = [c.strip() for c in ln.split(",")]
        if len(cells) != 2:
            raise ValueError(f"CSV {str(p)!r} の {ln_no} 行目は 2 列 (受け取った値: {ln!r}; カンマ区切り)")
        try:
            xs.append(float(cells[0]))
            ts.append(float(cells[1]))
        except ValueError:
            raise ValueError(f"CSV {str(p)!r} の {ln_no} 行目に数値でない値 ({ln!r})") from None
    return {hdr[0]: xs, "Tw_K": ts}, {"path": path_given, "abs_path": str(p), "sha256": hashlib.sha256(data).hexdigest()}


def _wall_isothermal(errs: list, w: dict, card_dir):
    """`wall: {thermal: isothermal, ...}` の検査。`Tw` (一様)・`Tw_profile` (表)・`Tw_profile_csv` (CSV) の厳密に 1 つ。
    表は `wall_temperature.parse_profile` に通し (書き直さない)、`{"thermal", "Tw_profile": {座標: [...], "Tw_K": [...]}}` の形に揃える。
    CSV は読んだ表を `Tw_profile` に、出典 (カードの値・絶対パス・sha256) を `Tw_profile_csv` に残す (card_log・`_card.card` の記録)。
    不正なら errs に足して None。"""
    from ..wall_temperature import parse_profile
    unknown = [str(k) for k in w if k not in ("thermal",) + WALL_TW_KEYS]
    if unknown:
        errs.append(f"wall (isothermal) の未知のキー {unknown} (書けるのは thermal と {' / '.join(WALL_TW_KEYS)} のどれか 1 つ)")
        return None
    given = [k for k in WALL_TW_KEYS if k in w]
    if len(given) != 1:
        errs.append("wall (isothermal) は Tw (一様 [K])・Tw_profile (表 {s|x_m: [...], Tw_K: [...]})・Tw_profile_csv (CSV のパス) の"
                    f"どれか 1 つ ({'どれも無い' if not given else f'併用している: {given}'})")
        return None
    k = given[0]
    if k == "Tw":
        return {"thermal": "isothermal", "Tw": _num(errs, "wall.Tw", w["Tw"], gt=0.0)}
    try:
        if k == "Tw_profile":
            src = None
            prof = parse_profile(_yaml12_floats(w["Tw_profile"]))
        else:
            tbl, src = read_profile_csv(w["Tw_profile_csv"], card_dir)
            prof = parse_profile(tbl)
    except ValueError as ex:
        errs.append(f"wall.{k}: {ex}")
        return None
    out = {"thermal": "isothermal"}
    if src is not None:
        out["Tw_profile_csv"] = src
    out["Tw_profile"] = {prof["coord"]: list(prof["coords"]), "Tw_K": list(prof["Tw_K"])}
    return out


def validate_card(card, card_dir=None) -> dict:
    """カード (dict) を検査し、既定値を埋めた写しを返す。不正なら全部の理由をまとめて ValueError。

    戻りの形はカードと同じキー (既定値を埋める: p_ambient = "design"、tolerances = DEFAULT_TOLERANCES、recipe = DEFAULT_RECIPE、
    overrides = {})。gas.species の名前は大文字。overrides は「節.キー」に平らにした mapping。
    card_dir: カードのファイルの所在 (`wall.Tw_profile_csv` の相対パスの基準; `load_card` が渡す)。"""
    errs: list = []
    if not isinstance(card, dict):
        raise ValueError(f"要求カード: 最上位は mapping (受け取った値: {type(card).__name__})")
    _mapping(errs, "要求カード", card, _TOP_REQUIRED, _TOP_OPTIONAL)
    out: dict = {}
    fmt = card.get("format")
    if "format" in card and fmt != CARD_FORMAT:
        errs.append(f"format は {CARD_FORMAT!r} (受け取った値: {fmt!r})")
    out["format"] = CARD_FORMAT
    name = card.get("name")
    if "name" in card and (not isinstance(name, str) or not name.strip()):
        errs.append(f"name は空でない文字列 (受け取った値: {name!r})")
    out["name"] = name
    # gas: preset (air) か species ({実種: 質量分率}) のどちらか 1 つ
    g = card.get("gas")
    if "gas" in card:
        if not isinstance(g, dict):
            errs.append(f"gas は mapping (受け取った値: {g!r})")
        else:
            unknown = [str(k) for k in g if k not in ("preset", "species")]
            if unknown:
                errs.append(f"gas の未知のキー {unknown} (preset か species のどちらか 1 つ)")
            if ("preset" in g) == ("species" in g):
                errs.append("gas は preset か species のどちらか 1 つを書く")
            elif "preset" in g:
                if g["preset"] not in GAS_PRESETS:
                    errs.append(f"gas.preset {g['preset']!r} は未知 (書けるのは {list(GAS_PRESETS)})")
                out["gas"] = {"preset": g["preset"]}
            else:
                sp = g["species"]
                if not isinstance(sp, dict) or not sp:
                    errs.append(f"gas.species は空でない mapping {{実種: 質量分率}} (受け取った値: {sp!r})")
                else:
                    from ..gas.composition import ResolvedSpeciesDB
                    known = set(ResolvedSpeciesDB.builtin().names)
                    ys = {}
                    for k, v in sp.items():
                        if isinstance(k, bool):
                            errs.append("gas.species のキーに真偽値がある (YAML で NO/ON はクォートする)")
                            continue
                        kk = str(k).upper()
                        if kk in ys:
                            errs.append(f"gas.species: 種 {kk} が 2 回 (大小文字違いも同じ種)")
                            continue
                        if kk not in known:
                            errs.append(f"gas.species: 種 {kk} は設計側の内蔵 DB に無い (書けるのは {sorted(known)})")
                            continue
                        y = _num(errs, f"gas.species.{kk}", v, gt=0.0)
                        if y is not None:
                            ys[kk] = y
                    if ys and len(ys) == len(sp):
                        tot = sum(ys.values())
                        if abs(tot - 1.0) > SPECIES_SUM_TOL:
                            errs.append(f"gas.species の質量分率の和 {tot!r} が 1 ± {SPECIES_SUM_TOL:g} でない")
                    out["gas"] = {"species": ys}
    c = card.get("conditions")
    if "conditions" in card and _mapping(errs, "conditions", c, ("Pt", "Tt")):
        out["conditions"] = {"Pt": _num(errs, "conditions.Pt", c["Pt"], gt=0.0),
                             "Tt": _num(errs, "conditions.Tt", c["Tt"], gt=0.0)}
    if "M_exit" in card:
        out["M_exit"] = _num(errs, "M_exit", card["M_exit"], gt=1.0)
    t = card.get("throat")
    if "throat" in card:
        if not isinstance(t, dict):
            errs.append(f"throat は mapping (受け取った値: {t!r})")
        else:
            unknown = [str(k) for k in t if k not in ("diameter_physical_m", "diameter_design_m") + THROAT_CURVATURE_KEYS]
            if unknown:
                errs.append(f"throat の未知のキー {unknown} (径は diameter_physical_m か diameter_design_m のどちらか 1 つ、"
                            f"曲率は {' か '.join(THROAT_CURVATURE_KEYS)} を任意で 1 つ)")
            given = [k for k in ("diameter_physical_m", "diameter_design_m") if k in t]
            if len(given) != 1:
                errs.append(f"throat は diameter_physical_m (物理スロート径、δ 込み) か diameter_design_m (設計 = 非粘性) の"
                            f"どちらか 1 つ (受け取ったキー: {given})")
            else:
                out["throat"] = {given[0]: _num(errs, f"throat.{given[0]}", t[given[0]], gt=0.0)}
            # スロートの壁の曲率 (設計壁; 2026-10-13 ユーザ「スロート曲率も指定できるようにして」): 実寸の曲率半径か r_t 比のどちらか
            cg = [k for k in THROAT_CURVATURE_KEYS if k in t]
            if len(cg) > 1:
                errs.append(f"throat の曲率は {' か '.join(THROAT_CURVATURE_KEYS)} のどちらか 1 つ (受け取ったキー: {cg})")
            elif cg and "throat" in out:
                out["throat"][cg[0]] = _num(errs, f"throat.{cg[0]}", t[cg[0]], gt=0.0)
    i = card.get("inlet")
    if "inlet" in card and _mapping(errs, "inlet", i, ("radius_m", "pipe_length_m", "contraction_length_m")):
        out["inlet"] = {"radius_m": _num(errs, "inlet.radius_m", i["radius_m"], gt=0.0),
                        "pipe_length_m": _num(errs, "inlet.pipe_length_m", i["pipe_length_m"], ge=0.0),
                        "contraction_length_m": _num(errs, "inlet.contraction_length_m", i["contraction_length_m"], gt=0.0)}
    if "divergent_length_m" in card:
        out["divergent_length_m"] = _num(errs, "divergent_length_m", card["divergent_length_m"], gt=0.0)
    w = card.get("wall")
    if "wall" in card:
        if not isinstance(w, dict) or w.get("thermal") not in ("isothermal", "adiabatic"):
            errs.append(f"wall は {{thermal: isothermal, Tw: K}}・{{thermal: isothermal, Tw_profile: {{s|x_m: [...], Tw_K: [...]}}}}・"
                        f"{{thermal: isothermal, Tw_profile_csv: <path>}} か {{thermal: adiabatic}} (受け取った値: {w!r})")
        elif w["thermal"] == "isothermal":
            wv = _wall_isothermal(errs, w, card_dir)
            if wv is not None:
                out["wall"] = wv
                # s の表は入口直管 [−2, −1] を区間で写すので直管長 > 0 が要る (codex result レビュー 2026-10-13 m4)
                prof = wv.get("Tw_profile")
                if isinstance(prof, dict) and "s" in prof and "inlet" in out \
                        and isinstance(out["inlet"].get("pipe_length_m"), float) and not out["inlet"]["pipe_length_m"] > 0.0:
                    errs.append("wall.Tw_profile の s の表は入口直管の長さ > 0 が要る (s ∈ [−2, −1] を直管に写す)。"
                                "直管なし (inlet.pipe_length_m = 0) なら x_m の表 (実寸、原点は設計スロート) を使う")
        elif _mapping(errs, "wall (adiabatic)", w, ("thermal",)):
            out["wall"] = {"thermal": "adiabatic"}
    pa = card.get("p_ambient", "design")
    if pa != "design":
        pa = _num(errs, "p_ambient", pa, gt=0.0)
    out["p_ambient"] = pa
    tol = card.get("tolerances")
    if tol is None:
        out["tolerances"] = dict(DEFAULT_TOLERANCES)
        out["tolerances_default"] = True
    elif _mapping(errs, "tolerances", tol, ("length_m", "diameter_m")):
        out["tolerances"] = {"length_m": _num(errs, "tolerances.length_m", tol["length_m"], gt=0.0),
                             "diameter_m": _num(errs, "tolerances.diameter_m", tol["diameter_m"], gt=0.0)}
        out["tolerances_default"] = False
    rec = card.get("recipe", DEFAULT_RECIPE)
    if not isinstance(rec, str) or not (RECIPE_DIR / f"{rec}.yaml").is_file():
        avail = sorted(p.stem for p in RECIPE_DIR.glob("*.yaml"))
        errs.append(f"recipe {rec!r} は無い (書けるのは {avail})")
    out["recipe"] = rec
    out["overrides"] = _flatten_overrides(errs, card.get("overrides"))
    # M_knot (plan tooling-nozzle-mknot-auto §4.4): 既定は自動 (書かない、または overrides の geometry.M_knot: auto; auto は overrides
    # から外して M_knot_mode に置く)。数値の overrides.geometry.M_knot は今と同じに効く (探索しない; 値はそのまま当てる)
    mk = out["overrides"].get("geometry.M_knot", _MISSING)
    if mk is _MISSING:
        out["M_knot_mode"], out["M_knot_auto_source"] = "auto", "カードに geometry.M_knot を書いていない (既定は自動)"
    elif isinstance(mk, str) and mk.strip() == MKNOT_AUTO:
        del out["overrides"]["geometry.M_knot"]
        out["M_knot_mode"], out["M_knot_auto_source"] = "auto", f"カードの overrides の geometry.M_knot: {MKNOT_AUTO}"
    elif isinstance(mk, bool) or not (isinstance(mk, (int, float)) or (isinstance(mk, str) and _YAML12_FLOAT.fullmatch(mk.strip()))):
        errs.append(f"overrides.geometry.M_knot は数値か {MKNOT_AUTO!r} (自動; 書かないのと同じ) (受け取った値: {mk!r})")
    else:
        out["M_knot_mode"] = "fixed"
    if errs:
        raise ValueError("要求カード: " + "; ".join(errs))
    return out


def load_card(path) -> dict:
    """カードのファイルを読んで `validate_card` を通す (`wall.Tw_profile_csv` の相対パスはこのファイルの所在が基準)。"""
    with open(path) as f:
        raw = yaml.safe_load(f)
    return validate_card(raw, card_dir=Path(path).resolve().parent)


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def throat_curvature_ratio(card: dict, rt: float):
    """カードのスロート曲率を R = 曲率半径 / r_t に (指定が無ければ None)。"""
    t = card["throat"]
    if "curvature_ratio" in t:
        return float(t["curvature_ratio"])
    if "curvature_radius_m" in t:
        return float(t["curvature_radius_m"]) / float(rt)
    return None


def throat_kind(card: dict) -> tuple:
    """("physical" | "design", 直径 [m])。"""
    t = card["throat"]
    if "diameter_physical_m" in t:
        return "physical", float(t["diameter_physical_m"])
    return "design", float(t["diameter_design_m"])


def card_species(card: dict) -> dict:
    """質量分率 {大文字の実種: Y} (preset は GAS_PRESETS の値)。"""
    g = card["gas"]
    if "preset" in g:
        return dict(GAS_PRESETS[g["preset"]])
    return dict(g["species"])


# ---------------------------------------------------------------------------------------------------- 雛形

def load_recipe(name: str) -> tuple:
    """雛形 (dict, パス, sha256)。`_recipe` (生成器への指示) の欠けは例外。"""
    path = RECIPE_DIR / f"{name}.yaml"
    if not path.is_file():
        raise ValueError(f"雛形 {name!r} が無い ({path})")
    with open(path) as f:
        rec = yaml.safe_load(f)
    meta = rec.get("_recipe") or {}
    for k in ("name", "M_knot_factor", "gas_transport", "tp_species_if_h2o"):
        if k not in meta:
            raise ValueError(f"雛形 {path}: _recipe.{k} が無い")
    if meta["name"] != name:
        raise ValueError(f"雛形 {path}: _recipe.name {meta['name']!r} がファイル名と違う")
    for sec in ("geometry", "mesh", "mesh_euler", "evaluate"):
        if not isinstance(rec.get(sec), dict):
            raise ValueError(f"雛形 {path}: 節 {sec} が無い")
    return rec, path, _sha256(path)


# ---------------------------------------------------------------------------------------------------- ガス・等エントロピー

def _gas(card: dict):
    from ..gas import GasSemiPerfect
    from ..gas.composition import ResolvedSpeciesDB
    return GasSemiPerfect(card_species(card), Tt=float(card["conditions"]["Tt"]), db=ResolvedSpeciesDB.builtin())


def _check_mach_in_table(gas, M: float) -> None:
    # GasSemiPerfect の超音速の表は T_min = max(Tt/14, 60 K) まで。表の外は np.interp が端の値を返す (黙って間違える) ので止める
    Mmax = float(gas._M[-1])
    if not M < Mmax:
        raise ValueError(f"M_exit {M} は semiperfect の表の上限 M {Mmax:.4f} (T_min = {float(gas._T[-1]):.1f} K) を超える "
                         "(GasSemiPerfect の T_min の既定 max(Tt/14, 60 K))")


def exit_state(card: dict) -> dict:
    """出口の等エントロピー状態 (semiperfect、M_exit): T [K]・p [Pa]。p/Pt は `evaluate.ic._isentropic_pressure_ratio_gas`
    (IC の 1 次元分布と同じ式: ln(p/Pt) = ∫_Tt^T c_p/(R T') dT')、T は `GasSemiPerfect.T_of_M`。"""
    from ..evaluate.ic import _isentropic_pressure_ratio_gas
    gas = _gas(card)
    M = float(card["M_exit"])
    _check_mach_in_table(gas, M)
    T = float(gas.T_of_M(M))
    pr = float(_isentropic_pressure_ratio_gas(np.array([T]), gas)[0])
    return {"M": M, "T_K": T, "p_Pa": float(card["conditions"]["Pt"]) * pr, "p_over_Pt": pr,
            "method": "semiperfect 等エントロピー (GasSemiPerfect.T_of_M + evaluate.ic._isentropic_pressure_ratio_gas)"}


def _mu_sutherland(T):
    mu0, T0, S = SUTHERLAND_AIR
    T = np.asarray(T, dtype=float)
    return mu0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)


# ---------------------------------------------------------------------------------------------------- カード → 問題 (純関数)

def _set_path(d: dict, path: str, value):
    """d の「a.b.c」に値を置き、前の値 (無ければ _MISSING) を返す。途中が mapping でなければ ValueError。"""
    parts = path.split(".")
    cur = d
    for p in parts[:-1]:
        if p not in cur or cur[p] is None:
            cur[p] = {}
        if not isinstance(cur[p], dict):
            raise ValueError(f"overrides.{path}: {p} が mapping でない (値 {cur[p]!r})")
        cur = cur[p]
    old = cur.get(parts[-1], _MISSING)
    cur[parts[-1]] = value
    return old


class _Missing:
    def __repr__(self):
        return "(無し)"


_MISSING = _Missing()


def _sig(v: float, n: int = 6) -> float:
    """有効 n 桁に丸める (参照値の表示用)。"""
    return float(f"{v:.{n - 1}e}")


def card_to_problem(card, r_t_m: float, M_knot: float | None = None) -> dict:
    """要求カードと r_t [m] → 問題 YAML の dict (純関数; 雛形のファイルを読むだけ)。§4.1 の 1・3・6。

    M_knot (plan tooling-nozzle-mknot-auto §4.4): 自動のカード (M_knot_mode auto) では `geometry.M_knot` に M_knot を置く
    (値は `build_problem` が規則 min_theta_max_v1 で選ぶ; 雛形の M_knot_factor は使わない)。knot の軸則で M_knot が無ければ
    ValueError (黙って 0.45 × M_exit に戻さない)。knot でない軸則では M_knot を使わないので `geometry.M_knot` を書かない。
    数値の overrides.geometry.M_knot のカードは今と同じ (雛形の M_knot_factor × M_exit を置いてから overrides で上書き)。

    - 長さは r_t で割る: r_inlet = R_in/r_t、L_U、L_pipe、dv.L_total = L_div/r_t、Lc_mode from_length、
      L_total_tol = 長さ公差/2/r_t。spec は Pt・Tt・p_ambient (design なら semiperfect の等エントロピー出口静圧)・r_throat・M_design・
      wall_thermal (壁温分布は `{mode: isothermal, Tw_profile: {s|x_m: [...], Tw_K: [...]}}`)・(物理スロート径なら) sizing {method: throat, target_m: 径/2}。
    - gas は semiperfect、species はカードの質量分率、transport は雛形 (実種ごとに cea、H2O は custom)、gamma・cp は参照値 (スロート)。
    - mesh の ni・wall_first_frac・wall_first_frac_throat と mesh_euler.ni は None (`build_problem` が `estimate_mesh` で埋める。
      None のまま `prepare` に渡すと例外になる — 黙って Mesh2DParams の既定に落とさない)。overrides で与えたらそちら。
    - overrides を最後に当て、上書きの記録を `_card.overrides` に置く。
    - 検査 (例外): R_in > r_t、μ = L_U²/(R Δr) ≤ 20 (R は雛形か overrides の geometry.R)、L_pipe ≥ 0、M_exit > 1、
      節点数の上限の許可 (mesh・mesh_euler の max_nodes と max_nodes_reason は組)、明示の mesh_euler.ni の節点数 (plan
      tooling-nozzle-mesh-size-guard §4.2)。"""
    card = validate_card(card) if "tolerances_default" not in card else card
    rt = float(r_t_m)
    if not (math.isfinite(rt) and rt > 0.0):
        raise ValueError(f"r_t は正の有限の数値 [m] (受け取った値: {r_t_m!r})")
    mk_auto = card.get("M_knot_mode") == "auto"
    if M_knot is not None and not mk_auto:
        raise ValueError("要求カード → 問題: M_knot= は M_knot が自動のカードだけに渡す (このカードは overrides の geometry.M_knot が数値)")
    rec, rec_path, rec_sha = load_recipe(card["recipe"])
    meta = rec["_recipe"]
    Y = card_species(card)
    gas = _gas(card)
    Pt, Tt = float(card["conditions"]["Pt"]), float(card["conditions"]["Tt"])
    M = float(card["M_exit"])
    if not M > 1.0:
        raise ValueError(f"M_exit {M} は 1 より大きい値")
    _check_mach_in_table(gas, M)
    if card["p_ambient"] == "design":
        es = exit_state(card)
        p_amb, p_amb_src = es["p_Pa"], f"design: 設計出口静圧 ({es['method']}、M {M}、p/Pt {es['p_over_Pt']:.6e})"
    else:
        p_amb, p_amb_src = float(card["p_ambient"]), "カードの値"
    kind, D = throat_kind(card)
    R_in = float(card["inlet"]["radius_m"])
    L_pipe_m = float(card["inlet"]["pipe_length_m"])
    L_U_m = float(card["inlet"]["contraction_length_m"])
    L_div_m = float(card["divergent_length_m"])
    tol_len = float(card["tolerances"]["length_m"])
    wall = card["wall"]
    if wall["thermal"] == "isothermal" and "Tw_profile" in wall:
        # 壁温分布: 表をそのまま写す (s は無次元・x_m は実寸なので r_t で割り直さない; CSV から来た表も埋め込み、YAML 単体で再現できる)
        wt = {"mode": "isothermal", "Tw_profile": copy.deepcopy(wall["Tw_profile"])}
    else:
        wt = {"mode": "isothermal", "Tw": float(wall["Tw"])} if wall["thermal"] == "isothermal" else {"mode": "adiabatic"}
    tr_rule = meta["gas_transport"]
    transport = {k: tr_rule.get(k, tr_rule["default"]) for k in Y}
    spec = {"Pt": Pt, "Tt": Tt, "p_ambient": p_amb, "r_throat": rt, "M_design": M, "wall_thermal": wt}
    if kind == "physical":
        spec["sizing"] = {"method": "throat", "target_m": D / 2.0}
    R_card = throat_curvature_ratio(card, rt)
    geom = {"R": (R_card if R_card is not None else float(rec["geometry"]["R"])), "r_inlet": R_in / rt, "L_U": L_U_m / rt, "L_pipe": L_pipe_m / rt,
            "Lc_mode": "from_length", "L_total_tol": tol_len / 2.0 / rt,
            "M_knot": ((None if M_knot is None else float(M_knot)) if mk_auto else float(meta["M_knot_factor"]) * M)}
    for k, v in rec["geometry"].items():
        if k not in geom:
            geom[k] = copy.deepcopy(v)
    mesh = copy.deepcopy(rec["mesh"])
    for k in _MESH_ESTIMATED:
        mesh[k] = None
    mesh_euler = copy.deepcopy(rec["mesh_euler"])
    mesh_euler["ni"] = None
    evaluate = copy.deepcopy(rec["evaluate"])
    if "H2O" in Y:
        evaluate["tp_species"] = meta["tp_species_if_h2o"]
    prob = {
        "name": str(card["name"]),
        "type": "wind_tunnel_axisym_axismach",
        "gas": {"model": "semiperfect", "transport": transport, "species": dict(Y),
                "gamma": _sig(float(gas.gamma_star)), "cp": _sig(float(gas.cp_mass(np.array([gas.T_star]))[0]))},
        "spec": spec,
        "dv": {"L_total": {"value": L_div_m / rt}},
        "geometry": geom,
        "mesh": mesh,
        "mesh_euler": mesh_euler,
        "evaluate": evaluate,
    }
    if R_card is not None and "geometry.R" in card["overrides"]:
        raise ValueError("要求カード → 問題: スロート曲率を throat と overrides.geometry.R の両方で与えている (どちらか 1 つ)")
    ov_log = []
    for path, val in card["overrides"].items():
        old = _set_path(prob, path, copy.deepcopy(val))
        ov_log.append({"key": path, "before": (None if old is _MISSING else old), "before_missing": old is _MISSING,
                       "value": val})
    # 検査 (§4.1 3.)。R は overrides 後の値
    g = prob["geometry"]
    R = float(g["R"])
    errs = []
    if not R_in > rt:
        errs.append(f"入口半径 {R_in} m が r_t {rt!r} m 以下 (R_in > r_t が必要)")
    else:
        dr = R_in / rt - 1.0
        mu = (L_U_m / rt) ** 2 / (R * dr)
        if mu > MU_MAX:
            errs.append(f"縮流部 μ = L_U²/(R Δr) = {mu:.4g} > {MU_MAX:g} (非単調収縮; L_U = {L_U_m / rt:.6g} r_t、R = {R:g}、"
                        f"Δr = {dr:.6g} r_t。縮流部長 ≤ {math.sqrt(MU_MAX * R * dr) * rt:.6g} m か R を大きく)")
    if not L_pipe_m >= 0.0:
        errs.append(f"直管長 {L_pipe_m} m が負")
    if R < R_HALL_MIN:
        errs.append(f"スロート曲率半径 / r_t = R {R:.4g} < {R_HALL_MIN} (遷音速解 Hall–KL の精度が落ちる範囲; 曲率半径 ≥ {R_HALL_MIN * rt:.6g} m)")
    # 節点数の上限の許可 (plan tooling-nozzle-mesh-size-guard §4.2): max_nodes と max_nodes_reason の組を重い処理の前に検査する。
    # 明示の mesh_euler.ni はここで節点数も判定する (mesh.ni の写しは build_problem が estimate_mesh の後に判定する)
    from ..evaluate.runner_axismach import MESH_EULER_DEFAULTS, _mesh_params_from, check_node_budget, node_budget
    for blk in ("mesh", "mesh_euler"):
        try:
            node_budget(prob.get(blk) or {}, blk)
        except ValueError as ex:
            errs.append(str(ex))
    if not errs and prob["mesh_euler"].get("ni") is not None:
        mp_e = _mesh_params_from(prob["mesh_euler"], 1.0, **MESH_EULER_DEFAULTS)
        try:
            check_node_budget(mp_e, prob["mesh_euler"], "mesh_euler", allow_hint=_card_allow_hint("mesh_euler"))
        except ValueError as ex:
            errs.append(f"明示の mesh_euler.ni {mp_e.ni}: {ex}")
    if errs:
        raise ValueError("要求カード → 問題: " + "; ".join(errs))
    if mk_auto:
        if str(g.get("axis_law", "quintic")) != "knot":
            g.pop("M_knot")                 # knot でない軸則は M_knot を使わない (探索もしない)
        elif M_knot is None:
            raise ValueError("要求カード → 問題: geometry.M_knot は自動 (規則 min_theta_max_v1; plan tooling-nozzle-mknot-auto) — 値は "
                             "build_problem が選ぶ。card_to_problem には選んだ値を M_knot= で渡す (固定するならカードの overrides に "
                             "geometry.M_knot: <数値>)")
    dr = R_in / rt - 1.0
    prob["_card"] = {
        "format": CARD_FORMAT, "generator": GENERATOR, "name": str(card["name"]),
        "recipe": {"name": card["recipe"], "path": str(rec_path), "sha256": rec_sha},
        "r_t_m": rt, "throat": {"kind": kind, "diameter_m": D,
                                "curvature_ratio_R": R, "curvature_radius_design_m": R * rt,
                                "curvature_source": ("card throat." + next(k for k in THROAT_CURVATURE_KEYS if k in card["throat"])
                                                     if R_card is not None else
                                                     ("overrides geometry.R" if "geometry.R" in card["overrides"] else "雛形"))},
        "p_ambient_source": p_amb_src,
        "mu_contraction": (L_U_m / rt) ** 2 / (R * dr),
        "tolerances": dict(card["tolerances"]), "tolerances_default": bool(card["tolerances_default"]),
        "overrides": ov_log,
    }
    return prob


def mapping_table(card: dict, prob: dict) -> list:
    """「カードの値 → YAML のキー」の対応 (標準出力と card_log.json に出す)。"""
    rt = float(prob["spec"]["r_throat"])
    kind, D = throat_kind(card)
    g, s = prob["geometry"], prob["spec"]
    w = card["wall"]
    if "Tw_profile" in w:
        # 壁温分布: mode の行と表の行 (表はそのまま写す。CSV なら読んだ表を埋め込む)
        tbl = s["wall_thermal"]["Tw_profile"]
        coord = next(k for k in ("s", "x_m") if k in tbl)
        how = f"表をそのまま写す ({_COORD_TEXT[coord]}; r_t で割り直さない。実寸の x への写しは呼ぶたびにその場の寸法で)"
        wall_rows = [("wall.thermal", w["thermal"], "spec.wall_thermal.mode", s["wall_thermal"]["mode"], "thermal → mode"),
                     (("wall.Tw_profile_csv", w["Tw_profile_csv"], "spec.wall_thermal.Tw_profile", tbl,
                       "CSV (相対パスはカードの所在が基準; 絶対パスと sha256 を記録) を読み、" + how + "。YAML 単体で再現できる")
                      if "Tw_profile_csv" in w else
                      ("wall.Tw_profile", w["Tw_profile"], "spec.wall_thermal.Tw_profile", tbl, how))]
    else:
        wall_rows = [("wall", w, "spec.wall_thermal", s["wall_thermal"], "thermal → mode")]
    rows = [
        ("conditions.Pt", card["conditions"]["Pt"], "spec.Pt", s["Pt"], "そのまま [Pa]"),
        ("conditions.Tt", card["conditions"]["Tt"], "spec.Tt", s["Tt"], "そのまま [K]"),
        ("M_exit", card["M_exit"], "spec.M_design", s["M_design"], "そのまま"),
        (f"throat.diameter_{kind}_m", D, "spec.r_throat", rt,
         ("物理スロート径: NS 後の補正の solve_rt_throat の固定点 (spec.sizing.target_m = 径/2; _card.ns_correction)"
          if prob["_card"].get("r_t_method") == "physical_throat_ns_corrected" else
          "物理スロート径: CFD 前の solve_rt_throat の固定点 (spec.sizing.target_m = 径/2)") if kind == "physical"
         else "設計スロート径 / 2"),
        (("throat." + next((k for k in THROAT_CURVATURE_KEYS if k in card["throat"]), "curvature (雛形)")),
         card["throat"].get("curvature_radius_m", card["throat"].get("curvature_ratio", "(無し)")), "geometry.R", g["R"],
         "曲率半径 / r_t (設計壁のスロートの曲率; 物理壁は境界層の分だけわずかに違う)"),
        ("inlet.radius_m", card["inlet"]["radius_m"], "geometry.r_inlet", g["r_inlet"], "/ r_t"),
        ("inlet.pipe_length_m", card["inlet"]["pipe_length_m"], "geometry.L_pipe", g["L_pipe"], "/ r_t"),
        ("inlet.contraction_length_m", card["inlet"]["contraction_length_m"], "geometry.L_U", g["L_U"], "/ r_t"),
        ("divergent_length_m", card["divergent_length_m"], "dv.L_total", prob["dv"]["L_total"]["value"],
         "/ r_t (Lc_mode: from_length)"),
        ("tolerances.length_m", card["tolerances"]["length_m"], "geometry.L_total_tol", g["L_total_tol"], "/ 2 / r_t"),
        *wall_rows,
        ("p_ambient", card["p_ambient"], "spec.p_ambient", s["p_ambient"], prob["_card"]["p_ambient_source"]),
        ("gas", card["gas"], "gas.species", prob["gas"]["species"], "質量分率 (model semiperfect、transport は雛形)"),
        ("recipe", card["recipe"], "geometry・mesh・mesh_euler・evaluate", "(雛形の構成キー)",
         (f"M_knot は自動 (雛形 {prob['_card']['recipe']['name']} の M_knot_factor は使わない)" if card.get("M_knot_mode") == "auto" else
          f"M_knot = {prob['_card']['recipe']['name']} の係数 × M_exit")),
    ]
    mkr = prob["_card"].get("M_knot_rule")
    if mkr is not None:
        rows.append(("M_knot (自動)", MKNOT_AUTO, "geometry.M_knot", g.get("M_knot"), _mknot_rule_text(mkr)))
    rows += [(f"overrides.{o['key']}", o["value"], o["key"], o["value"], f"上書き (雛形/生成値: {o['before']!r})")
             for o in prob["_card"]["overrides"]]
    return [{"card_key": a, "card_value": b, "yaml_key": c, "yaml_value": d, "rule": e} for a, b, c, d, e in rows]


# ---------------------------------------------------------------------------------------------------- 壁解像の初期見積もり

def _wall_first_profile(xs, mesh: dict, frac: float, frac_t: float | None) -> np.ndarray:
    """mesh2d.generate_axisym_mesh と同じ規則の断面ごとの第 1 層比 (blend・上流の戻し)。"""
    xs = np.asarray(xs, dtype=float)
    if frac_t is None:
        return np.full(xs.shape, float(frac))
    b0, b1 = float(mesh.get("wall_first_blend_x0", 0.5)), float(mesh.get("wall_first_blend_x1", 6.0))
    t = np.clip((xs - b0) / max(b1 - b0, 1e-9), 0.0, 1.0)
    t = t * t * (3.0 - 2.0 * t)
    u0, u1 = mesh.get("wall_first_up_x0"), mesh.get("wall_first_up_x1")
    if u0 is not None and u1 is not None:
        tu = np.clip((float(u1) - xs) / max(float(u1) - float(u0), 1e-9), 0.0, 1.0)
        tu = tu * tu * (3.0 - 2.0 * tu)
        t = np.maximum(t, tu)
    return frac_t + (frac - frac_t) * t


def _floor_sig2(v: float) -> float:
    """有効 2 桁に切り下げ (第 1 層は薄い側に丸める)。"""
    e = math.floor(math.log10(v))
    return float(f"{math.floor(v / 10.0 ** (e - 1) + 1e-9)}e{e - 1}")


def _wall_ar_fast(wall, ni: int, mesh: dict, frac: float, frac_t) -> tuple:
    """station を mesh2d._x_stations で置き、壁際のセルの AR = max(辺)/min(辺) (辺: 壁の弦・両側の第 1 層) の最大と位置 x。"""
    from ..meshing.mesh2d import _x_stations
    xs = _x_stations(wall.x_in, wall.x_e, int(ni), float(mesh.get("throat_refine", 3.0)), float(mesh.get("throat_width", 1.5)))
    rw = wall.r(xs)
    y1 = _wall_first_profile(xs, mesh, frac, frac_t) * rw
    chord = np.hypot(np.diff(xs), np.diff(rw))
    e = np.stack([chord, y1[:-1], y1[1:]])
    ar = e.max(axis=0) / e.min(axis=0)
    k = int(np.argmax(ar))
    return float(ar[k]), float(0.5 * (xs[k] + xs[k + 1]))


def _ar_mesh_params(ni: int, mesh: dict, frac: float, frac_t):
    """`_mesh_ar_full` が作る全格子の Mesh2DParams (節点数の事前の判定もこれで数える)。"""
    from ..meshing.mesh2d import Mesh2DParams
    opt = lambda k: None if mesh.get(k) is None else float(mesh[k])  # noqa: E731
    return Mesh2DParams(ni=int(ni), nj=int(mesh.get("nj", 81)), wall_first_frac=float(frac),
                        throat_refine=float(mesh.get("throat_refine", 3.0)), throat_width=float(mesh.get("throat_width", 1.5)),
                        wall_first_frac_throat=(None if frac_t is None else float(frac_t)),
                        wall_first_blend_x0=float(mesh.get("wall_first_blend_x0", 0.5)),
                        wall_first_blend_x1=float(mesh.get("wall_first_blend_x1", 6.0)),
                        wall_first_up_x0=opt("wall_first_up_x0"), wall_first_up_x1=opt("wall_first_up_x1"),
                        axis_gap_frac=opt("axis_gap_frac"), axis_cap_frac=opt("axis_cap_frac"), scale=1.0)


def _card_allow_hint(block: str) -> str:
    """要求カードで節点数の上限を許す書き方 (plan tooling-nozzle-mesh-size-guard §4.2)。{n_nodes} は check_node_budget が実効の節点数に置き換える。"""
    return (f"カードの overrides に {block}.max_nodes と {block}.max_nodes_reason を組で書く "
            f"(例 overrides: {{{block}.max_nodes: {{n_nodes}}, {block}.max_nodes_reason: \"…\"}})")


def _mesh_ar_full(wall, ni: int, mesh: dict, frac: float, frac_t) -> dict:
    """実際の mesh2d.generate_axisym_mesh で格子を作り、全セルの AR (check_mesh_quality.py と同じ max 辺/min 辺) を測る。"""
    from ..meshing.mesh2d import generate_axisym_mesh
    prm = _ar_mesh_params(ni, mesh, frac, frac_t)
    coords, quads, _ = generate_axisym_mesh(wall, prm)
    P = coords[:, :2][quads]
    e = np.linalg.norm(np.roll(P, -1, axis=1) - P, axis=2)
    ar = e.max(axis=1) / e.min(axis=1)
    k = int(np.argmax(ar))
    return {"ar_max": float(ar[k]), "x_at_max": float(P[k, :, 0].mean()), "n_cells": int(len(quads)), "nj": prm.nj}


def _ni_for_ar(wall, mesh: dict, frac: float, frac_t, ar_max: float) -> int:
    """壁際のセルの AR (`_wall_ar_fast`) ≤ ar_max となる最小の ni (AR(ni) は ni にほぼ反比例で単調減。倍々で上限を作ってから二分法)。"""
    lo, hi = 2, 64
    while _wall_ar_fast(wall, hi, mesh, frac, frac_t)[0] > ar_max and hi < 64 * NI_CAP:
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if _wall_ar_fast(wall, mid, mesh, frac, frac_t)[0] > ar_max:
            lo = mid
        else:
            hi = mid
    return hi


def _y1plus_at_stations(wall, ni: int, mesh: dict, frac: float, frac_t, x_dense, frac_req) -> tuple:
    """station (`mesh2d._x_stations`) ごとの y1+ の目安 = 0.8 × (配点の第 1 層比) / (y1+ 0.8 に要る第 1 層比)。入口の節点 (平板の起点で
    要る第 1 層が 0 に落ちる) は除く。戻り (最大, その x [r_t])。"""
    from ..meshing.mesh2d import _x_stations
    xs = _x_stations(wall.x_in, wall.x_e, int(ni), float(mesh.get("throat_refine", 3.0)), float(mesh.get("throat_width", 1.5)))[1:]
    y1p = Y1PLUS_TARGET * _wall_first_profile(xs, mesh, frac, frac_t) / np.interp(xs, x_dense, frac_req)
    k = int(np.argmax(y1p))
    return float(y1p[k]), float(xs[k])


def estimate_mesh(card, problem_dict: dict, chain: dict | None = None) -> dict:
    r"""壁解像の第 1 層と ni の**初期見積もり** (§4.1 4.)。合格の根拠ではない — NS 後に `check_wall_resolution.py` の VERDICT で
    判定し、不合格なら細分する (手順書 §2 ⑥)。

    1 次元の等エントロピー (面積比は `design_chain` の設計壁の $r_w(x)^2$、亜音速枝 x < 0・超音速枝 x ≥ 0、semiperfect の表; IC と同じ式)
    で $M,T,p,u$ を作り、Eckert の参照温度 $T^*=T_e(1+0.032M^2+0.58(T_w/T_e-1))$ と平板乱流 $C_f^*=0.0592\,Re_x^{*-0.2}$
    ($x$ = 入口端からの壁の弧長、μ は Sutherland の空気) で $\tau_w = C_f^*\,\rho^* u_e^2/2$、$u_\tau=\sqrt{\tau_w/\rho_w}$、
    $y_1 = 0.8\,\mu_w/(\rho_w u_\tau)$。断熱壁の $T_w$ は回復温度 $T_e + 0.89\,(T_t - T_e)$。
    - **壁温分布** (`wall.Tw_profile`; plan tooling-nozzle-wall-temperature-profile §4.3): 断面ごとの局所の $T_w(x)$ を使う。x は実寸
      (断面の x [r_t] × r_t)。s の表の区切りの点は**その設計の壁** (`geometry.L_pipe`・`L_U`・設計壁の終点 `chain["wall"].x_e` に r_t を
      掛けたもの; `wall_temperature.anchors_from_design`、prepare_ns と同じ写し) から取る (カードの長さではない: 設計壁の終点はカードの
      拡大部長と公差の内でずれ、超えることもある; plan §4.3 [2026-10-13 訂正])。戻り値の `wall_temperature` に区切りの点と出典を残す。
      表が見積もりの範囲 (設計壁の入口端〜終点) を覆うかを `check_covers` で検査し、補間は `eval_Tw` (端値で延長しない)。
    - 戻り値の `sections` (断面ごとの x・局所の $T_w$・丸め前の $y_1$ と $y_1/r_w$) は呼ぶ側の検査用 (card_log には書かない)。
    - `wall_first_frac` = スロート下流 x > 6 の $y_1/r_w$ の最小、`wall_first_frac_throat` = |x| ≤ 2 の最小 (有効 2 桁に切り下げ)。
    - **縮流部側のブレンド** (`wall_first_up_x0/x1`): 雛形の値の配点で、station ごとの y1+ の目安 (入口の節点を除く) が 1 を超えたら、
      問題の形から決め直す。候補を順に試し、目安の最大 ≤ 1 になった最初のものを使う:
      (B) up_x0 = 入口端 x_in・up_x1 = 縮流部の入口 −L_U (縮流部の全体をスロートの値、直管は入口の wall_first_frac からブレンド)、
      (C) null (mesh2d の規則で上流は全域スロートの値)。どれも 1 を超えたら C にして警告。overrides で与えたブレンドは変えない (警告だけ)。
    - `ni` = 壁際のセルの AR (max 辺/min 辺) ≤ `mesh.ar_max` (既定 5000) となる最小 (station は `mesh2d._x_stations` で実際に置き、
      第 1 層は mesh2d と同じブレンドの規則)。10 の倍数に切り上げ、実際の `generate_axisym_mesh` の全セルの AR で確かめる。
      上限 NI_CAP (6000)、超えたら上限にして警告。
    - **節点数の上限** (plan tooling-nozzle-mesh-size-guard §4.2): 全格子 (`_mesh_ar_full`) を作る前に、同じ Mesh2DParams の実効の節点数
      (`mesh2d.node_count`) を `runner_axismach.check_node_budget` で判定し、上限 (既定 30 万、`mesh.max_nodes` と理由で許可) を超えたら
      全格子を作らずに ValueError。戻り値の `nodes` に節点数と上限。
    chain: `design_chain` の結果 (無ければ problem_dict から作る)。"""
    card = validate_card(card) if "tolerances_default" not in card else card
    if chain is None:
        chain = _design_chain_of(problem_dict)
    from ..evaluate.ic import _invert_area_ratio_gas, _isentropic_pressure_ratio_gas
    wall = chain["wall"]
    rt = float(problem_dict["spec"]["r_throat"])
    mesh = dict(problem_dict["mesh"])
    gas = _gas(card)
    Pt, Tt = float(card["conditions"]["Pt"]), float(card["conditions"]["Tt"])
    x = np.linspace(float(wall.x_in), float(wall.x_e), 8001)
    rw = wall.r(x)
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(rw)))]) * rt      # 入口端からの弧長 [m]
    AR = np.maximum(rw, 1.0) ** 2
    Mx = _invert_area_ratio_gas(AR, x >= 0.0, gas)
    T = np.where(Mx >= 1.0, gas.T_of_M(np.maximum(Mx, 1.0)), np.interp(np.minimum(Mx, 1.0), gas._Mu, gas._Tu))
    P = Pt * _isentropic_pressure_ratio_gas(T, gas)
    Rg = float(gas.R)
    u = Mx * np.sqrt(gas.gamma(T) * Rg * T)
    tw_rec = None
    if card["wall"]["thermal"] == "isothermal" and "Tw_profile" in card["wall"]:
        import types
        from ..wall_temperature import anchors_from_design, check_covers, eval_Tw, parse_profile, to_x_m
        prof = parse_profile(card["wall"]["Tw_profile"])
        # s の区切りの点はその設計の壁から (geometry.L_pipe・L_U と設計壁の終点 x_e に r_t を掛ける; prepare_ns と同じ写し、
        # plan §4.3 [2026-10-13 訂正])。anchors_from_design は p.geometry だけを読むので、問題 dict の geometry を渡す
        anc = anchors_from_design(types.SimpleNamespace(geometry=problem_dict["geometry"]), chain, rt)
        x_tab, T_tab = to_x_m(prof, anc)
        xm = x * rt                                        # 断面の x [m] (原点は設計スロート)
        check_covers(x_tab, float(xm[0]), float(xm[-1]), "格子の見積もりの範囲 (設計壁の入口端〜終点)")
        Tw = eval_Tw(x_tab, T_tab, xm)
        tw_note = (f"等温壁 分布 ({prof['coord']} の表 {len(prof['coords'])} 点、T_w {min(prof['Tw_K']):g}〜{max(prof['Tw_K']):g} K)")
        tw_rec = {"coord": prof["coord"], "n_points": len(prof["coords"]),
                  "anchors_m": ({k: float(v) for k, v in anc.items()} if prof["coord"] == "s" else None),
                  "anchors_source": ("design wall: geometry.L_pipe·r_t・geometry.L_U·r_t・設計壁の終点 x_e·r_t "
                                     "(prepare_ns の anchors_from_design と同じ写し; plan §4.3)"
                                     if prof["coord"] == "s" else "x_m の表 (区切りの点を使わない)"),
                  "x_m_table": [float(v) for v in x_tab], "Tw_K_table": [float(v) for v in T_tab],
                  "range_m": [float(xm[0]), float(xm[-1])],
                  "wall_end_minus_card_exit_m": float(xm[-1]) - float(card["divergent_length_m"])}   # 記録 (写しには使わない)
    elif card["wall"]["thermal"] == "isothermal":
        Tw = np.full_like(T, float(card["wall"]["Tw"]))
        tw_note = f"等温壁 Tw {float(card['wall']['Tw'])} K"
    else:
        Tw = T + RECOVERY_FACTOR * (Tt - T)
        tw_note = f"断熱壁: 回復温度 T_e + {RECOVERY_FACTOR}·(Tt − T_e)"
    Ts = T * (1.0 + 0.032 * Mx ** 2 + 0.58 * (Tw / T - 1.0))
    rho_s, mu_s = P / (Rg * Ts), _mu_sutherland(Ts)
    rho_w, mu_w = P / (Rg * Tw), _mu_sutherland(Tw)
    with np.errstate(divide="ignore", invalid="ignore"):
        Re = rho_s * u * s / mu_s
        Cf = 0.0592 * Re ** -0.2
        tau = Cf * 0.5 * rho_s * u * u
        utau = np.sqrt(tau / rho_w)
        y1 = Y1PLUS_TARGET * mu_w / (rho_w * utau)
    frac_req = y1 / (rw * rt)
    m_dn = x > X_DOWNSTREAM_MIN
    m_th = np.abs(x) <= X_THROAT_HALF
    if not m_dn.any() or not m_th.any():
        raise ValueError(f"estimate_mesh: 設計壁 [{wall.x_in:.4g}, {wall.x_e:.4g}] r_t に x > {X_DOWNSTREAM_MIN} か |x| ≤ "
                         f"{X_THROAT_HALF} の断面が無い")
    k_dn = int(np.flatnonzero(m_dn)[np.argmin(frac_req[m_dn])])
    k_th = int(np.flatnonzero(m_th)[np.argmin(frac_req[m_th])])
    warnings = []
    ov = (problem_dict.get("_card") or {}).get("overrides") or []
    ov_keys = {o["key"] for o in ov}
    frac = float(mesh["wall_first_frac"]) if ("mesh.wall_first_frac" in ov_keys) else _floor_sig2(float(frac_req[k_dn]))
    frac_t = (None if mesh.get("wall_first_frac_throat") is None and "mesh.wall_first_frac_throat" in ov_keys else
              float(mesh["wall_first_frac_throat"]) if "mesh.wall_first_frac_throat" in ov_keys else
              _floor_sig2(float(frac_req[k_th])))
    ar_max = float(mesh.get("ar_max", 1000))
    ni_fixed = int(mesh["ni"]) if "mesh.ni" in ov_keys else None
    # 縮流部側のブレンドの候補 (雛形 → B → C)。overrides で与えたブレンドはそのまま
    blend_fixed = bool({"mesh.wall_first_up_x0", "mesh.wall_first_up_x1"} & ov_keys)
    x_in, xU = float(wall.x_in), -float(wall.up.L_U)
    cands = [("recipe", "雛形 (または問題) の値", mesh.get("wall_first_up_x0"), mesh.get("wall_first_up_x1"))]
    if not blend_fixed:
        cands += [("B", "縮流部の全体をスロートの値 (up_x0 = 入口端 x_in、up_x1 = 縮流部の入口 −L_U)", x_in, xU),
                  ("C", "上流は全域スロートの値 (up_x0・up_x1 = null; mesh2d の規則)", None, None)]
    tried = []
    chosen = None
    for label, text, u0, u1 in cands:
        mc = dict(mesh, wall_first_up_x0=u0, wall_first_up_x1=u1)
        ni_req_c = _ni_for_ar(wall, mc, frac, frac_t, AR_MARGIN * ar_max)
        ni_c = ni_fixed if ni_fixed is not None else min(NI_CAP, int(math.ceil(ni_req_c / 10.0) * 10))
        yp, xp = _y1plus_at_stations(wall, ni_c, mc, frac, frac_t, x, frac_req)
        tried.append({"label": label, "text": text, "wall_first_up_x0": u0, "wall_first_up_x1": u1, "ni_required_for_ar": ni_req_c,
                      "ni": ni_c, "y1plus_est_max": yp, "x_y1plus_est_max_rt": xp})
        if yp <= 1.0:
            chosen = tried[-1]
            break
    if chosen is None:
        chosen = tried[-1]
        warnings.append(f"格子: 第 1 層の配点のどの候補でも y1+ の目安の最大が 1 を超える (選んだ {chosen['label']}: "
                        f"{chosen['y1plus_est_max']:.2f} @ x = {chosen['x_y1plus_est_max_rt']:.2f} r_t) — NS 後の check_wall_resolution.py で確かめる")
    mesh_c = dict(mesh, wall_first_up_x0=chosen["wall_first_up_x0"], wall_first_up_x1=chosen["wall_first_up_x1"])
    ni_req = chosen["ni_required_for_ar"]
    ni = chosen["ni"]
    if ni_fixed is None and ni_req > NI_CAP:
        warnings.append(f"格子: 壁際の AR ≤ {ar_max:g} に要る ni {ni_req} が上限 {NI_CAP} を超える — ni を {NI_CAP} にした "
                        f"(AR {_wall_ar_fast(wall, NI_CAP, mesh_c, frac, frac_t)[0]:.0f}; 第 1 層・ar_max・throat_refine を見直す)")
    from ..evaluate.runner_axismach import check_node_budget

    def guard_nodes(ni_: int) -> dict:
        # 節点数の上限 (plan tooling-nozzle-mesh-size-guard §4.2): 全格子 (_mesh_ar_full) を作る前に、同じ Mesh2DParams の実効の節点数で判定する
        prm = _ar_mesh_params(ni_, mesh_c, frac, frac_t)
        try:
            return check_node_budget(prm, mesh, "mesh", allow_hint=_card_allow_hint("mesh"))
        except ValueError as ex:
            src = ("overrides の mesh.ni" if ni_fixed is not None else
                   f"見積もり (壁際の AR ≤ {AR_MARGIN} × {ar_max:g} に要る ni {ni_req}、上限 NI_CAP {NI_CAP}、第 1 層 {frac:g}・スロート {frac_t})")
            raise ValueError(f"要求カード → 問題 (estimate_mesh、NS の格子 ni {ni_} は {src}): {ex}") from ex

    nodes = guard_nodes(ni)
    full = _mesh_ar_full(wall, ni, mesh_c, frac, frac_t)
    n_bump = 0
    while full["ar_max"] > ar_max and ni < NI_CAP and ni_fixed is None and n_bump < 50:
        ni = min(NI_CAP, ni + 10)
        nodes = guard_nodes(ni)
        full = _mesh_ar_full(wall, ni, mesh_c, frac, frac_t)
        n_bump += 1
    if full["ar_max"] > ar_max:
        warnings.append(f"格子: 実際の格子 (ni {ni}) の AR の最大 {full['ar_max']:.0f} が ar_max {ar_max:g} を超える "
                        f"(x = {full['x_at_max']:.3f} r_t)")
    yp, xp = _y1plus_at_stations(wall, ni, mesh_c, frac, frac_t, x, frac_req)
    if chosen["label"] != "recipe":
        warnings.append(f"格子: 雛形のブレンド (wall_first_up_x0/x1 {mesh.get('wall_first_up_x0')}/{mesh.get('wall_first_up_x1')}) では "
                        f"y1+ の目安の最大 {tried[0]['y1plus_est_max']:.2f} (x = {tried[0]['x_y1plus_est_max_rt']:.2f} r_t) > 1 — "
                        f"{chosen['label']} ({chosen['text']}) に決め直した (目安の最大 {yp:.2f} @ x = {xp:.2f} r_t)")
    elif blend_fixed and yp > 1.0:
        warnings.append(f"格子: overrides のブレンドで y1+ の目安の最大 {yp:.2f} > 1 (x = {xp:.2f} r_t)")
    # float32 の座標では第 1 層 / 局所半径 ~3e-6 で ~27 ulp (誤差 ±数 %)。1e-6 を下回る第 1 層は FP64 のビルドと msh の 17 桁が要る
    # (止めない; 冷却壁の第 1 層の精度の注意、2026-10-11)
    for nm, fv in (("wall_first_frac", frac), ("wall_first_frac_throat", frac_t)):
        if fv is not None and fv < FLOAT32_FRAC_WARN:
            warnings.append(f"格子: {nm} {fv:.2g} < {FLOAT32_FRAC_WARN:g} — float32 の座標では第 1 層が局所半径の {fv:.2g} 倍で "
                            f"約 {fv / np.finfo(np.float32).eps:.0f} ulp しかない。FP64 のビルドと msh の 17 桁の座標が要る "
                            "(mesh2d.write_msh41_2d の既定は 10 桁)")
    from ..meshing.mesh2d import _radial_fracs
    nj = int(mesh.get("nj", 81))

    def q_of(fr):
        g = np.diff(_radial_fracs(nj, fr))
        return float(np.max(g[:-1] / g[1:]))
    res = {
        "note": "初期値 (合格の根拠ではない)。NS 後に check_wall_resolution.py の VERDICT で判定し、不合格なら細分する (手順書 §2 ⑥)",
        "wall_first_frac": frac, "wall_first_frac_throat": frac_t, "ni": int(ni),
        "wall_first_up_x0": chosen["wall_first_up_x0"], "wall_first_up_x1": chosen["wall_first_up_x1"],
        "blend_choice": chosen["label"], "blend_text": chosen["text"], "blend_candidates": tried,
        "wall_first_frac_raw": float(frac_req[k_dn]), "x_wall_first_frac_rt": float(x[k_dn]),
        "wall_first_frac_throat_raw": float(frac_req[k_th]), "x_wall_first_frac_throat_rt": float(x[k_th]),
        "y1_m_throat": float(y1[k_th]), "y1_m_downstream_min": float(y1[k_dn]),
        "ni_required_for_ar": int(ni_req), "ar_max_limit": ar_max,
        "ar_max_wall_fast": _wall_ar_fast(wall, ni, mesh_c, frac, frac_t)[0], "mesh_check": full,
        "nodes": nodes,              # NS の格子の実効の節点数と上限 (check_node_budget; plan tooling-nozzle-mesh-size-guard §4.2)
        "y1plus_est_max_with_blend": yp, "x_y1plus_est_max_rt": xp,
        "radial_growth_max": {"downstream": q_of(frac), "throat": (q_of(frac_t) if frac_t is not None else None)},
        "inputs": {"y1plus_target": Y1PLUS_TARGET, "wall": tw_note, "mu": "Sutherland (空気: 1.716e-5 Pa s, 273.15 K, S 110.4 K)",
                   "Cf": "0.0592 Re_x*^-0.2 (x = 入口端からの壁の弧長)", "T_ref": "Eckert T* = T_e(1 + 0.032 M² + 0.58 (T_w/T_e − 1))",
                   "flow_1d": "設計壁の面積比 r_w² の 1 次元等エントロピー (semiperfect)",
                   "y1plus_check": "mesh2d._x_stations の station ごと (入口の節点を除く)、配点は mesh2d と同じブレンド"},
        "warnings": warnings,
    }
    if tw_rec is not None:
        # 壁温分布: 表・区切りの点・第 1 層を決めた断面の局所の T_w (一様・断熱の戻り値・card_log は変えない)
        tw_rec.update(Tw_K_at_wall_first_frac=float(Tw[k_dn]), Tw_K_at_wall_first_frac_throat=float(Tw[k_th]))
        res["wall_temperature"] = tw_rec
    # 断面ごとの値 (丸め前; 検査用の配列。build_problem は card_log に書かない)
    res["sections"] = {"x_rt": x, "x_m": x * rt, "Tw_K": Tw, "y1_m": y1, "y1_over_rw": frac_req}
    return res


# ---------------------------------------------------------------------------------------------------- 液化の事前警告

def liquefaction_check(card) -> dict:
    """出口の静温・静圧 (等エントロピー) と N₂ の飽和蒸気圧 (Antoine、63〜126 K) を比べる (§4.1 5.)。N₂ 分圧が飽和を超えたら
    警告を返す (止めない。判断はユーザ)。範囲外・N₂ を含まないときはそう書く。O₂・CO₂ などは見ない。"""
    card = validate_card(card) if "tolerances_default" not in card else card
    from ..gas.composition import ResolvedSpeciesDB, mass_to_mole
    es = exit_state(card)
    T, p = es["T_K"], es["p_Pa"]
    X = mass_to_mole(card_species(card), ResolvedSpeciesDB.builtin())
    out = {"T_exit_K": T, "p_exit_Pa": p, "M_exit": es["M"], "method": es["method"],
           "antoine_N2": {"A": N2_ANTOINE[0], "B": N2_ANTOINE[1], "C": N2_ANTOINE[2], "range_K": list(N2_ANTOINE_RANGE_K),
                          "form": "log10(P/bar) = A − B/(T + C)"},
           "x_N2": None, "p_N2_Pa": None, "p_sat_N2_Pa": None, "supersaturated": None, "warning": None}
    if "N2" not in X:
        out["verdict"] = "N2 を含まないので判定しない"
        return out
    xN2 = float(X["N2"])
    pN2 = xN2 * p
    out.update(x_N2=xN2, p_N2_Pa=pN2)
    lo, hi = N2_ANTOINE_RANGE_K
    if T > hi:
        out["verdict"] = f"出口静温 {T:.2f} K が Antoine の範囲 ({lo:g}〜{hi:g} K) より高い (N₂ の臨界温度 126.2 K 付近より上) — 液化しない"
        out["supersaturated"] = False
        return out
    if T < lo:
        out["verdict"] = f"出口静温 {T:.2f} K が Antoine の範囲 ({lo:g}〜{hi:g} K) より低い — 判定できない (外挿しない)"
        out["warning"] = (f"液化: 出口静温 {T:.2f} K が N₂ の Antoine の範囲の下限 {lo:g} K (三重点 63.15 K 付近) より低い — "
                          "N₂ の凝縮・固化の可能性を判定できない。Tt を上げるか、凝縮を評価する")
        return out
    A, B, C = N2_ANTOINE
    psat = 1.0e5 * 10.0 ** (A - B / (T + C))
    out["p_sat_N2_Pa"] = psat
    out["supersaturated"] = bool(pN2 > psat)
    ratio = pN2 / psat
    out["S_N2"] = ratio
    if pN2 > psat:
        out["verdict"] = f"過飽和: N₂ 分圧 {pN2:.4g} Pa > 飽和 {psat:.4g} Pa (S = {ratio:.3g}、出口 {T:.2f} K)"
        out["warning"] = (f"液化: 出口 (M {es['M']}、{T:.2f} K、{p:.4g} Pa) で N₂ 分圧 {pN2:.4g} Pa が飽和蒸気圧 {psat:.4g} Pa を超える "
                          f"(S = {ratio:.3g}) — 空気の凝縮の恐れ。Tt を上げるか凝縮を評価する (判断はユーザ)")
    else:
        out["verdict"] = f"未飽和: N₂ 分圧 {pN2:.4g} Pa < 飽和 {psat:.4g} Pa (S = {ratio:.3g}、出口 {T:.2f} K)"
    return out


# ---------------------------------------------------------------------------------------------------- 設計壁から測る寸法

def measure_design_dimensions(card, chain: dict, r_t_m: float, sizing: dict | None = None) -> dict:
    """`design_chain` の設計壁から寸法を測り、カードの値と公差で比べる (§4.2a の出荷ゲートの前段; 最終の STEP・メッシュのゲートは
    成果物の書き出しで別に行う)。入口半径 = r(x_in)、直管長 = (−L_U) − x_in (設計壁の接合位置; 区間で r が一定かも確かめる)、
    縮流部長 = 接合位置 → 設計スロート (r の最小の位置)、スロート→出口長 = 設計壁の終点 x_e (設計スロートが原点)。
    スロート径: 物理スロート径のカードは `solve_rt_throat` の物理スロート半径 × 2 (sizing)、設計スロート径のカードは設計壁の最小半径 × 2。
    MOC の出口 x_F (終端特性線と壁の交点、Lc_mode: from_length が合わせる量) と x_e の差も記録する。"""
    from scipy.optimize import minimize_scalar
    card = validate_card(card) if "tolerances_default" not in card else card
    W = chain["wall"]
    rt = float(r_t_m)
    tl, td = float(card["tolerances"]["length_m"]), float(card["tolerances"]["diameter_m"])
    x_in, xU, x_e = float(W.x_in), -float(W.up.L_U), float(W.x_e)
    xp = np.linspace(x_in, xU, 2001)
    r_in = float(W.r(np.array([x_in]))[0])
    pipe_dev = float(np.max(np.abs(W.r(xp) - r_in)))
    mres = minimize_scalar(lambda v: float(W.r(np.array([v]))[0]), bounds=(-0.5, 0.5), method="bounded",
                           options={"xatol": 1e-12})
    x_t, r_t_min = float(mres.x), float(mres.fun)
    kind, D = throat_kind(card)
    rows = []

    def row(name, card_v, meas, tol, how):
        rows.append({"name": name, "card_m": card_v, "measured_m": meas, "diff_m": meas - card_v, "tol_m": tol,
                     "pass": bool(abs(meas - card_v) <= tol), "how": how})

    row("inlet_radius", float(card["inlet"]["radius_m"]), r_in * rt, td / 2.0, "設計壁の r(x_in) × r_t")
    row("pipe_length", float(card["inlet"]["pipe_length_m"]), (xU - x_in) * rt, tl,
        f"設計壁の直管/縮流部の接合位置 − 入口端 (直管区間の |r − r_in| の最大 {pipe_dev:.2e} r_t)")
    row("contraction_length", float(card["inlet"]["contraction_length_m"]), (x_t - xU) * rt, tl,
        f"接合位置 → 設計スロート (r の最小 x = {x_t:.3e} r_t)")
    row("divergent_length", float(card["divergent_length_m"]), x_e * rt, tl, "設計壁の終点 x_e × r_t (設計スロートが原点)")
    if kind == "physical":
        if sizing is not None:
            ns_after = str(sizing.get("delta_r_route", "")).startswith("ns_after")
            row("throat_diameter_physical", D, 2.0 * float(sizing["physical_throat"]["r_m"]), td,
                f"solve_rt_throat の物理スロート半径 × 2 ({'NS 後の補正の δ_E' if ns_after else 'CFD 前の δ'})")
    else:
        row("throat_diameter_design", D, 2.0 * r_t_min * rt, td, "設計壁の最小半径 × 2")
    x_F = (chain.get("exit") or {}).get("x_F")
    return {"rows": rows, "pass": all(r["pass"] for r in rows), "r_t_m": rt,
            "x_e_rt": x_e, "x_F_moc_rt": (None if x_F is None else float(x_F)),
            "x_F_minus_x_e_m": (None if x_F is None else (float(x_F) - x_e) * rt),
            "Lc_solve": chain.get("Lc_solve"), "x_throat_design_rt": x_t, "r_throat_design_rt": r_t_min,
            "note": "設計壁から測った寸法 (§4.2a の前段)。最終の合格は成果物 (STEP・メッシュ) から測る出荷ゲートで判定する"}


# ---------------------------------------------------------------------------------------------------- 書き出し

def _flow(v) -> str:
    s = yaml.safe_dump(v, default_flow_style=True, allow_unicode=True, width=1 << 30, sort_keys=False).strip()
    if s.endswith("..."):
        s = s[:-3].strip()
    return s


def _one_line(s: str) -> str:
    return " ".join(str(s).split())


def render_problem_yaml(prob: dict, header: list, comments: dict) -> str:
    """問題 dict を YAML の文字列にする (各キーは flow 形式 1 行、コメントつき)。safe_load で読み直すと prob と同じになる。
    comments: {"節.キー" (最上位は "キー"): コメント}。"""
    out = [f"# {_one_line(h)}" if h else "#" for h in header]
    for top in _PROBLEM_ORDER:
        if top not in prob:
            continue
        v = prob[top]
        c = comments.get(top)
        if isinstance(v, dict):
            out.append("")
            out.append(f"{top}:" + (f"   # {_one_line(c)}" if c else ""))
            for k, vv in v.items():
                cc = comments.get(f"{top}.{k}")
                out.append(f"  {k}: {_flow(vv)}" + (f"   # {_one_line(cc)}" if cc else ""))
        else:
            out.append(f"{top}: {_flow(v)}" + (f"   # {_one_line(c)}" if c else ""))
    extra = [k for k in prob if k not in _PROBLEM_ORDER]
    if extra:
        raise ValueError(f"render_problem_yaml: 想定外の最上位キー {extra}")
    text = "\n".join(out) + "\n"
    back = yaml.safe_load(text)
    if back != prob:
        raise RuntimeError("render_problem_yaml: 書いた YAML を読み直すと元の dict と一致しない")
    return text


def _wall_comment(card: dict) -> str:
    """spec.wall_thermal のコメント (一様・断熱は従来どおり「カード wall」)。"""
    w = card["wall"]
    if "Tw_profile" not in w:
        return "カード wall"
    tbl = w["Tw_profile"]
    coord = next(k for k in ("s", "x_m") if k in tbl)
    desc = f"{coord} の表 {len(tbl[coord])} 点、T_w {min(tbl['Tw_K']):g}〜{max(tbl['Tw_K']):g} K; {_COORD_TEXT[coord]}"
    if "Tw_profile_csv" in w:
        src = w["Tw_profile_csv"]
        return f"カード wall.Tw_profile_csv {src['abs_path']} (sha256 {src['sha256'][:16]}…) の表を埋め込んだ ({desc}; YAML 単体で再現できる)"
    return f"カード wall.Tw_profile をそのまま写した ({desc})"


def _key_comments(card: dict, prob: dict) -> dict:
    """生成した各キーのコメント (出典・換算)。"""
    rt = float(prob["spec"]["r_throat"])
    rec = prob["_card"]["recipe"]["name"]
    kind, D = throat_kind(card)
    c = {
        "name": "要求カードの name",
        "gas": "カードの gas (model は semiperfect 固定)",
        "gas.model": "生成器 (カードの gas は常に semiperfect)",
        "gas.transport": f"雛形 {rec}: 実種ごとに cea (H2O は custom:h2o_iapws_cea_v1)",
        "gas.species": f"質量分率 (カードの gas: {_flow(card['gas'])})",
        "gas.gamma": "参照値: スロート T* の γ (semiperfect では NASA-9 が優先)",
        "gas.cp": "参照値: c_p(T*) [J/(kg K)]",
        "spec.Pt": "カード conditions.Pt [Pa]",
        "spec.Tt": "カード conditions.Tt [K]",
        "spec.p_ambient": f"[Pa] {prob['_card']['p_ambient_source']}",
        "spec.M_design": "カード M_exit",
        "spec.wall_thermal": _wall_comment(card),
        "spec.sizing": f"カード throat.diameter_physical_m {D} m の半径 (物理スロート) で寸法を決めた",
        "dv": "Lc_mode: from_length の設計変数 (r_t 単位)",
        "dv.L_total": f"= divergent_length_m {card['divergent_length_m']} m / r_t (設計スロート x=0 → 物理出口)",
        "geometry.R": (f"= throat.curvature_radius_m {card['throat']['curvature_radius_m']} m / r_t (設計壁のスロートの曲率半径)"
                       if "curvature_radius_m" in card["throat"] else
                       f"= throat.curvature_ratio (設計壁のスロートの曲率半径 / r_t)" if "curvature_ratio" in card["throat"] else
                       f"雛形 {rec} (スロート曲率半径 / r_t = {prob['geometry']['R']:g}; 曲率半径 {prob['geometry']['R'] * prob['spec']['r_throat']:.6g} m)"),
        "geometry.r_inlet": f"= inlet.radius_m {card['inlet']['radius_m']} m / r_t",
        "geometry.L_U": f"= inlet.contraction_length_m {card['inlet']['contraction_length_m']} m / r_t "
                        f"(μ = L_U²/(R Δr) = {prob['_card']['mu_contraction']:.4g} ≤ {MU_MAX:g})",
        "geometry.L_pipe": f"= inlet.pipe_length_m {card['inlet']['pipe_length_m']} m / r_t",
        "geometry.Lc_mode": "生成器: dv.L_total から L_c を解く",
        "geometry.L_total_tol": f"= 長さ公差 {card['tolerances']['length_m']} m / 2 / r_t",
        "geometry.M_knot": f"雛形 {rec}: M_knot_factor × M_design",
        "mesh": "NS の格子 (prepare_ns)。ni・wall_first_frac・wall_first_frac_throat は estimate_mesh の初期値",
        "mesh_euler": "Euler の格子 (prepare)。軸方向 (ni・throat_refine・throat_width) は mesh と同じ",
        "_card": "要求カードの記録 (ソルバ・設計チェーンは読まない)",
    }
    if kind == "physical":
        c["spec.r_throat"] = f"[m] 物理スロート径 {D} m から {prob['_card'].get('r_t_method_text', '')}"
    else:
        c["spec.r_throat"] = f"[m] = 設計スロート径 {D} m / 2 (カード throat.diameter_design_m)"
    for k in prob["geometry"]:
        c.setdefault(f"geometry.{k}", f"雛形 {rec}")
    for sec in ("mesh", "mesh_euler", "evaluate"):
        for k in prob.get(sec, {}):
            c.setdefault(f"{sec}.{k}", f"雛形 {rec}")
    est = prob["_card"].get("mesh_estimate")
    if est:
        c["mesh.ni"] = f"estimate_mesh の初期値: 壁際の AR ≤ {est['ar_max_limit']:g} の下限 (要 {est['ni_required_for_ar']})"
        c["mesh.wall_first_frac"] = (f"estimate_mesh の初期値: x > {X_DOWNSTREAM_MIN:g} の y1+ {Y1PLUS_TARGET} の y1/r_w の最小 "
                                     f"{est['wall_first_frac_raw']:.3e} を 2 桁に切り下げ")
        c["mesh.wall_first_frac_throat"] = (f"estimate_mesh の初期値: |x| ≤ {X_THROAT_HALF:g} の最小 "
                                            f"{est['wall_first_frac_throat_raw']:.3e} を 2 桁に切り下げ")
        c["mesh_euler.ni"] = "= mesh.ni (軸方向は mesh と同じ)"
        if est.get("blend_choice") not in (None, "recipe"):
            for k in ("wall_first_up_x0", "wall_first_up_x1"):
                c[f"mesh.{k}"] = (f"estimate_mesh が決め直した ({est['blend_choice']}; null = 上流は全域スロートの値)。"
                                  "雛形の値では y1+ の目安 > 1")
    mkr = prob["_card"].get("M_knot_rule")
    if mkr is not None and mkr.get("mode") == "auto":
        c["geometry.M_knot"] = _mknot_rule_text(mkr) + "。要約は _card.M_knot_rule、掃引の表と全部の記録は card_log の M_knot_rule"
    for a in prob["_card"].get("auto_adjust") or []:
        c[a["key"]] = f"自動調整: {a['from']!r} → {a['to']!r} ({_one_line(a['reason'])[:160]})"
    rl = prob["_card"].get("L_total_tol_relaxed")
    if rl:
        c["geometry.L_total_tol"] = (f"**内部許容差を緩めた**: カードの長さ公差 {card['tolerances']['length_m']} m / 2 / r_t = {rl['from']:.4g} → "
                                     f"{rl['to']} r_t ({STALL_TEXT})。カードの公差は変えていない — 寸法の合否はカードの公差で決める")
    if "tp_species" in prob.get("evaluate", {}):
        c["evaluate.tp_species"] = f"雛形 {rec}: 組成に H2O があるので split_h2o (凝縮の restart 用)"
    for o in prob["_card"]["overrides"]:
        parts = o["key"].split(".")
        key = ".".join(parts[:2])
        before = "(無し)" if o["before_missing"] else _flow(o["before"])
        c[key] = f"overrides で上書き (雛形/生成値: {before})"
    return c


def _write_problem(card: dict, prob: dict, out: Path, header: list) -> str:
    text = render_problem_yaml(prob, header, _key_comments(card, prob))
    out.write_text(text)
    return text


# from_length の停止の言い方 (2026-10-13 codex 諮問 notes/reviews/2026-10-13-wt-length-accuracy-diagnose.md の採用): design_chain の
# from_length は最初の非改善で反復を打ち切るので、到達の限界 (離散化のノイズ床) とは区別できない
STALL_TEXT = ("from_length の固定点反復が改善しなくなって止まった (runner_axismach.py:606-616; 離散化のノイズ床とは区別できない)")
L_TOTAL_TOL_RELAXED = 0.05           # from_length の内部許容差を緩めるときの値 [r_t] (runner_axismach の既定と同じ)
X_END_MARGIN_STEP = 0.5              # x_F が MOC の場の終点に張り付いたときの geometry.x_end_margin の引き上げ幅
X_END_MARGIN_MAX = 5.0               # その上限 (超えても解消しなければ例外)
_MAX_ATTEMPTS = 20                   # 自動調整のやり直しの上限 (安全弁)


def _is_from_length_stall(ex: Exception) -> bool:
    """design_chain の Lc_mode: from_length が L_total_tol に届かずに止まった例外か。"""
    return isinstance(ex, ValueError) and "Lc_mode: from_length 未収束" in str(ex)


def _is_field_short(ex: Exception) -> bool:
    """逆 MOC の場 (軸の延長 x_end_margin) が足りず、終端特性線が壁に届かない例外か。"""
    return (isinstance(ex, RuntimeError) and "終端特性線が壁に到達しない" in str(ex)) or \
        (isinstance(ex, ValueError) and "x_F 追跡に失敗" in str(ex))


def moc_call_record(res: dict, kw: dict) -> dict:
    """逆 MOC (`inverse_design`) の 1 回の結果の出口と場の終点 (r_t 単位)。`x_e` = 切った壁 (= 設計壁) の終点、`x_wall_full_end` =
    切る前の壁 (場で決まった壁) の終点、`spacing_end` = その最後の点間隔。`pinned` = x_F が場の終点から 1 点間隔以内
    (場の終点に張り付いている疑い。cplus の閉包は壁の範囲の外の交点を拾わないが、streamline の terminal_exit は壁の外を端の半径で
    外挿する; moc_inverse.py:558)。"""
    wf = np.asarray(res["wall_full"], dtype=float)
    w = np.asarray(res["wall"], dtype=float)
    ex = res.get("exit") or {}
    x_F = float(ex.get("x_F", float("nan")))
    x_e = float(w[-1, 0])
    x_end = float(wf[-1, 0])
    sp = float(wf[-1, 0] - wf[-2, 0]) if len(wf) > 1 else float("nan")
    anc = kw.get("axis_anchor")
    # cplus の閉包では x_F は終端 C⁺ 自身の壁点に乗る。丸めで x_F がその点のわずか下になると、切る規則 (wall[:, 0] <= x_F) が
    # その点を落とし、x_e が 1 点間隔だけ短くなる (別 plan の切り分け材料として記録するだけ; 振る舞いは変えない)
    i_n = int(np.argmin(np.abs(wf[:, 0] - x_F))) if np.isfinite(x_F) else 0
    d_n = float(x_F - wf[i_n, 0]) if np.isfinite(x_F) else float("nan")
    return {"x_E": (None if kw.get("x_E") is None else float(kw["x_E"])),
            "x_A": (None if anc is None else float(anc[0])),
            "x_F": x_F, "x_e": x_e, "x_F_minus_x_e_rt": x_F - x_e,
            "x_wall_full_end": x_end, "spacing_end": sp, "margin_to_field_end_rt": x_end - x_F,
            "n_wall_full": int(len(wf)), "n_wall": int(len(w)), "n_beyond_x_F": int(np.sum(wf[:, 0] > x_F)),
            "x_axis_end": (None if kw.get("x_axis_end") is None else float(kw["x_axis_end"])),
            "n_frozen": ex.get("n_frozen"), "pinned": bool(not np.isfinite(x_F) or (x_end - x_F) <= sp),
            "x_F_minus_nearest_node_rt": d_n,
            "cut_dropped_node_at_x_F": bool(np.isfinite(d_n) and -1e-9 * max(1.0, abs(x_F)) <= d_n < 0.0 and x_e < wf[i_n, 0])}


class _MocCapture:
    """design_chain の中の逆 MOC (`runner_axismach.inverse_design`) の各呼び出しを記録する (読むだけ)。design_chain は切る前の壁
    (wall_full) を返さないので、場の終点と x_F の関係はここでしか見えない (設計チェーン本体は変えない)。"""

    def __init__(self):
        self.calls: list = []

    def __enter__(self):
        from ..evaluate import runner_axismach as RA
        self._RA, self._orig = RA, RA.inverse_design

        def wrap(*a, **k):
            r = self._orig(*a, **k)
            self.calls.append(moc_call_record(r, k))
            return r
        RA.inverse_design = wrap
        return self

    def __exit__(self, *exc):
        self._RA.inverse_design = self._orig
        return False


def _chain_field(chain: dict, calls: list) -> dict:
    """design_chain が採った逆 MOC の結果 (x_F が一致する呼び出し) の記録。"""
    xF = float(chain["exit"]["x_F"])
    hit = [c for c in calls if c["x_F"] == xF]
    if not hit:
        raise RuntimeError(f"_chain_field: design_chain の x_F {xF!r} に一致する逆 MOC の呼び出しが記録に無い")
    return dict(hit[-1])


def _evals_log(calls: list, r_t: float, L_target: float | None) -> list:
    """from_length の各評価 (逆 MOC の各呼び出し) の (L_c, x_F, x_e, x_F − x_e [r_t と m]) と場の終点 (card_log 用)。"""
    out = []
    for c in calls:
        L_c = (None if c["x_E"] is None or c["x_A"] is None else c["x_E"] - c["x_A"])
        out.append({"L_c": L_c, "x_F_rt": c["x_F"], "x_e_rt": c["x_e"], "x_F_minus_x_e_rt": c["x_F_minus_x_e_rt"],
                    "x_F_minus_x_e_m": c["x_F_minus_x_e_rt"] * r_t,
                    "x_F_minus_L_total_rt": (None if L_target is None else c["x_F"] - L_target),
                    "x_wall_full_end_rt": c["x_wall_full_end"], "spacing_end_rt": c["spacing_end"],
                    "n_beyond_x_F": c["n_beyond_x_F"], "x_axis_end_rt": c["x_axis_end"], "pinned": c["pinned"],
                    "x_F_minus_nearest_node_rt": c["x_F_minus_nearest_node_rt"],
                    "cut_dropped_node_at_x_F": c["cut_dropped_node_at_x_F"]})
    return out


def _design_chain_of(prob: dict, capture: "_MocCapture | None" = None):
    """問題 dict を一時ファイルに書いて load_problem → design_chain (CFD なし)。"""
    from ..evaluate.runner_axismach import design_chain
    from ..probdef import load_problem
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "problem.yaml"
        path.write_text(yaml.safe_dump(prob, sort_keys=False, allow_unicode=True))
        if capture is None:
            return design_chain(load_problem(path))
        with capture:
            return design_chain(load_problem(path))


class _Retry(Exception):
    """自動調整をしたので最初からやり直す。"""


# ---------------------------------------------------------------------------------------------------- M_knot の自動の選択 (記録)

MKNOT_GATE_HALF_WIDTHS = (0.1, 0.2, 0.3)   # ゲートの窓の半幅の列 (§4.2 窓の端の選択; mknot_select.GATE_HALF_WIDTHS と同じ値を明示して渡す)
MKNOT_GATE_HALF_WIDTH = MKNOT_GATE_HALF_WIDTHS[0]   # 詰めとゲートの最初の窓の半幅 (最終の組の N の文脈 ±(この値 + ΔM_safe) の幅)
_MKNOT_ROW_KEYS = ("M_K", "ok", "reasons", "role", "theta_max_deg", "margin_deg", "min_sin_angle", "mpp_sign_changes", "x_K", "L_c",
                   "moc_gate_ok", "moc_gate_pass", "n_axis")


def _mem_available_gb() -> float | None:
    """/proc/meminfo の MemAvailable [GB] (`free` の available と同じ量)。読めなければ None。"""
    try:
        for ln in Path("/proc/meminfo").read_text().splitlines():
            if ln.startswith("MemAvailable:"):
                return float(ln.split()[1]) / 1024.0 ** 2
    except OSError:
        return None
    return None


def mknot_workers(N: int) -> dict:
    """M_knot の評価の並列の数 (MKNOT_WORKERS: N・2N とも 1、2026-10-14 の取り決め) と、開始前に使えるメモリ (MemAvailable) と
    1 評価のメモリの見積もり (MKNOT_EVAL_MEM_GB_2400 × (N/2400)^MKNOT_EVAL_MEM_EXP; 記録と `_mknot_mem_check` に使う)。"""
    return {"N": int(MKNOT_WORKERS[0]), "2N": int(MKNOT_WORKERS[1]), "available_gb": _mem_available_gb(),
            "per_eval_gb_N": MKNOT_EVAL_MEM_GB_2400 * (N / 2400.0) ** MKNOT_EVAL_MEM_EXP,
            "per_eval_gb_2N": MKNOT_EVAL_MEM_GB_2400 * (2.0 * N / 2400.0) ** MKNOT_EVAL_MEM_EXP}


def _mknot_mem_check(wk: dict, N: int) -> None:
    """2N の 1 評価のメモリの見積もりが使えるメモリを超えるなら、評価を始める前に止める (メモリ不足で落ちる前に理由を出す)。
    n_axis_inv の倍化 (from_length の停止の自動調整) の後は 2N が 4 倍の解像度になる。"""
    if wk["available_gb"] is not None and wk["per_eval_gb_2N"] > wk["available_gb"]:
        raise RuntimeError(f"M_knot: 解像度のゲートの 2N = {2 * N} の 1 評価のメモリの見積もり {wk['per_eval_gb_2N']:.1f} GB が使えるメモリ "
                           f"{wk['available_gb']:.1f} GB を超える — メモリの大きい環境 (AWS) で生成するか、M_knot を数値で与える "
                           "(overrides の geometry.M_knot)")


def _plain(v):
    """記録用: numpy → Python、tuple → list、有限でない浮動小数 → None (YAML に書いて読み直すと同じ dict になるように)。"""
    v = _json_safe(v)
    if isinstance(v, dict):
        return {k: _plain(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_plain(x) for x in v]
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


def _mknot_row(r: dict | None) -> dict | None:
    """評価の行の記録 (掃引の表の 1 行: M_K・成立・理由・θ_max・余裕・最小 sin 角・M″ の符号反転・x_K・MOC のゲート ほか)。"""
    if r is None:
        return None
    out = {k: r.get(k) for k in _MKNOT_ROW_KEYS if k in r}
    out["reasons"] = [str(s)[:200] for s in (r.get("reasons") or [])]
    return out


def _mknot_bounds(bs: list) -> list:
    return [{k: b.get(k) for k in ("kind", "feasible_end", "lo", "hi", "source", "nearest_other", "reasons_infeasible") if k in b}
            for b in bs]


def _mknot_select_record(s: dict) -> dict:
    keep = ("rule", "M_A", "M_d", "settings", "chosen", "distance_to_boundary", "reason", "constrained", "theta_min",
            "theta_min_admissible", "tie_set", "excluded_by_distance", "safe_points", "local_minima", "windows", "warnings", "n_eval",
            "elapsed_s")
    rec = {k: s.get(k) for k in keep}
    rec["boundaries"] = _mknot_bounds(s.get("boundaries") or [])
    return rec


def _mknot_gate_record(g: dict) -> dict:
    # 窓の記録 (§5.1 #3b): window・half_width は最後の窓、n_widen は広げた回数、widen は半幅ごとの選び直しと端だった値、
    # edge_values は最後の窓で端だった値 (θ_max の最小の位置の判定 min_offset_* とは別)
    keep = ("N", "N2", "chosen", "window", "half_width", "n_widen", "widen", "edge_values", "settings", "pass", "checks", "reasons",
            "warnings", "chosen_N", "chosen_N_reason", "chosen_2N", "chosen_2N_reason", "moved_N", "distance_to_conservative_boundary",
            "distance_2N", "distance_input", "min_offset_N", "min_offset_2N", "elapsed_s")
    rec = {k: g.get(k) for k in keep}
    rec["boundaries_all"] = _mknot_bounds(g.get("common_boundaries") or [])
    rec["rows_selected"] = {m: {"at_N": _mknot_row(v["at_N"]), "at_2N": _mknot_row(v["at_2N"])}
                            for m, v in (g.get("rows_selected") or {}).items()}
    return rec


def _mknot_dist(v) -> str:
    """境目からの距離の表示 (境目が無ければ記録は None = 無限大)。"""
    return "∞ (境目なし)" if v is None else f"{v:.4f}"


def _mknot_widen_text(g: dict) -> str:
    """ゲートの窓を広げた記録の表示 (広げていなければ空)。g はゲートの記録 (card_log) か YAML の要約のゲート。"""
    n = g.get("n_widen") or 0
    if not n:
        return ""
    edges = "; ".join(f"半幅 {w['half_width']:g} で " + "・".join(f"{k} {v!r}" for k, v in w["edge_values"].items()) + " が端"
                      for w in g["widen"] if w["edge_values"])
    return f"窓を {n} 回広げた (半幅 {g['widen'][0]['half_width']:g} → {g['half_width']:g}; {edges})"


def _mknot_rule_text(mkr: dict) -> str:
    """M_knot の記録の 1 行の説明 (geometry.M_knot のコメント・対応表)。mkr は card_log の記録 (全部) か YAML の要約。"""
    if mkr.get("mode") != "auto":
        return f"自動 ({mkr.get('rule')}) — {mkr.get('reason')}"
    g = mkr["gate"]
    wt = _mknot_widen_text(g)
    return (f"自動 (規則 {mkr['rule']}): 最終の組 (r_t {mkr['r_t_m']!r} m) の解像度のゲートで選び直した N の値 chosen_N "
            f"(θ_max {mkr['theta_max_deg']:.5f}°、全部の境目から {_mknot_dist(mkr['distance_to_boundary'])} ≥ ΔM_safe {mkr['dM_safe']:g}、"
            f"ゲート {'PASS' if g['pass'] else 'FAIL'} (n_axis_inv {g['N']}/{g['N2']}: {g['chosen_N']!r}/{g['chosen_2N']!r}"
            + (f"; {wt}" if wt else "") + f"); 雛形の M_knot_factor {mkr['recipe_M_knot_factor']:g} は使わない)")


def _mknot_summary_lines(mkr: dict, card_log_name: str | None = None) -> list:
    """M_knot の記録の要約 (YAML のコメント・標準出力)。mkr は card_log の記録 (全部)。card_log_name があれば、掃引の表と全部の
    記録の置き場 (card_log のファイル名と記録の sha256; YAML の _card.M_knot_rule.card_log と同じ値) の行を足す。"""
    if mkr.get("mode") != "auto":
        return [f"M_knot: {_mknot_rule_text(mkr)}"]
    s, g = mkr["selection"], mkr["gate"]
    wt = _mknot_widen_text(g)
    out = [f"M_knot: 自動 (規則 {mkr['rule']}; {mkr['source']}) → geometry.M_knot {mkr['chosen']!r} (θ_max {mkr['theta_max_deg']:.5f}°、"
           f"余裕 {mkr['margin_deg']:.3f}°、網の最小 sin 角 {mkr['min_sin_angle']:.4f})。要約は _card.M_knot_rule、掃引の表と全部の記録は "
           "card_log の M_knot_rule",
           f"  境目 (N・2N を全部、成立側の端): "
           + "、".join(f"{b['kind']} {b['feasible_end']:.6f} ({b['source']})" for b in g["boundaries_all"])
           + f" — 選んだ値から {_mknot_dist(mkr['distance_to_boundary'])} (ΔM_safe {mkr['dM_safe']:g})",
           f"  解像度のゲート (最終の組 r_t {mkr['r_t_m']!r} m、n_axis_inv {g['N']}/{g['N2']}、窓 {g['window'][0]:.6g}〜{g['window'][1]:.6g}): "
           f"{'PASS' if g['pass'] else 'FAIL'} — 選び直し N {g['chosen_N']!r}・2N {g['chosen_2N']!r} (入力 {g['chosen']!r} から "
           f"{g['moved_N']:+.4g})" + (f"; {wt}" if wt else ""),
           f"  探索 (最初の r_t {mkr['first_r_t_m']!r} m、n_axis_inv {mkr['n_axis_inv']}): {s['n_eval']} 点、M_A {s['M_A']:.6g}〜"
           f"M_d {s['M_d']:.6g}、選択 {s['chosen']!r} ({_one_line(s['reason'])[:160]})" + (f"; 警告 {s['warnings']}" if s["warnings"] else "")]
    for e in mkr["events"]:
        out.append(f"  選び直し [{e['stage']}] r_t {e['r_t_m']!r} m・M_K の入力 {e.get('M_K_in')!r} → chosen_N {e.get('chosen_N')!r}: "
                   f"{_one_line(e['action'])[:200]}")
    out.append(f"  M_knot の所要 {mkr['elapsed_s']:.0f} s (並列 N {mkr['workers']['N']}・2N {mkr['workers']['2N']}、開始時の使えるメモリ "
               f"{mkr['workers']['available_gb'] if mkr['workers']['available_gb'] is None else round(mkr['workers']['available_gb'], 1)} GB)")
    if card_log_name is not None:
        out.append(f"  掃引の表と全部の記録: {card_log_name} の M_knot_rule (その記録の sha256 {_mknot_record_sha256(mkr)}; "
                   "YAML の _card.M_knot_rule は要約だけ)")
    return out


def _mknot_record_sha256(rec: dict) -> str:
    """card_log の `M_knot_rule` (記録の全部) の sha256。JSON にしてから測る (キーの順・区切りを固定、ensure_ascii=False; card_log を
    読み直した dict からも同じ値になる)。card_log のファイルは YAML の sha256 (out_sha256) を含むので、ファイルそのものの sha256 は
    YAML に書けない (循環) — 記録の部分を測る。"""
    txt = json.dumps(_json_safe(rec), sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(txt.encode("utf-8")).hexdigest()


def _mknot_yaml_summary(rec: dict, card_log_name: str) -> dict:
    """YAML の `_card.M_knot_rule` に置く要約 (plan tooling-nozzle-mknot-auto §4.4、2026-10-14 諮問 mknot-card-holes の採用): 規則の
    名前・選んだ値 (と出力の点の θ_max・余裕・最小 sin 角・N/2N の成立)・距離・全部の境目・ゲートの判定 (窓を広げた回数・端だった値を
    含む)・r_t の固定点の前後の選び直しの要約・card_log の記録の sha256。掃引の表と行の詳細 (table・ゲートの rows_selected・
    constraints_at_output の行・選び方の記録・選び直しの文) は card_log の `M_knot_rule` だけに置く (掃引の表を YAML に入れると
    1 行 約 52 kB になり、Euler の写しにも入った)。rec は `_mknot_record` の戻り (card_log に書く記録)。
    knot でない軸則の記録 (mode not_applicable) は表を持たないので、そのまま返す。照合は `mknot_card_log_check`。"""
    if rec.get("mode") != "auto":
        return rec

    def gate(x):
        return {"pass": x["pass"], "N": x["N"], "N2": x["N2"], "chosen": x["chosen"], "chosen_N": x["chosen_N"],
                "chosen_2N": x["chosen_2N"], "moved_N": x["moved_N"], "distance_N": x["distance_to_conservative_boundary"],
                "distance_2N": x["distance_2N"], "window": x["window"], "half_width": x["half_width"], "n_widen": x["n_widen"],
                "widen": [{k: w[k] for k in ("half_width", "chosen_N", "chosen_2N", "edge_values")} for w in x["widen"]],
                "checks": x["checks"], "warnings": x["warnings"]}

    s = rec["selection"]
    return _plain({
        "rule": rec["rule"], "mode": rec["mode"], "source": rec["source"], "chosen": rec["chosen"], "chosen_from": rec["chosen_from"],
        "r_t_m": rec["r_t_m"], "first_r_t_m": rec["first_r_t_m"], "n_axis_inv": rec["n_axis_inv"],
        "theta_max_deg": rec["theta_max_deg"], "theta_max_deg_final_chain": rec["theta_max_deg_final_chain"],
        "margin_deg": rec["margin_deg"], "min_sin_angle": rec["min_sin_angle"],
        "ok_at_output": {k: rec["constraints_at_output"][k]["ok"] for k in ("N", "2N")},
        "distance_to_boundary": rec["distance_to_boundary"], "dM_safe": rec["dM_safe"],
        "boundaries": [{k: b.get(k) for k in ("kind", "feasible_end", "lo", "hi", "source")} for b in rec["boundaries"]],
        "gate": gate(rec["gate"]),
        "gate_first": (None if rec["gate_first"] is None else       # 最初の r_t のゲート (最終の組と別のとき) は判定と値だけ
                       {k: rec["gate_first"][k] for k in ("pass", "chosen", "chosen_N", "chosen_2N", "n_widen", "half_width")}),
        "gate_first_reused": rec["gate_first_reused"],
        "selection": {k: s[k] for k in ("M_A", "M_d", "chosen", "distance_to_boundary", "n_eval", "warnings")},
        "events": [{k: e.get(k) for k in ("stage", "r_t_m", "M_K_in", "M_K_used_in_fixed_point", "chosen_N", "chosen_2N", "diff_vs_used",
                                          "gate_pass", "gate_n_widen", "distance")} for e in rec["events"]],
        "redo": rec["redo"], "elapsed_s": rec["elapsed_s"],
        "recipe_M_knot_factor": rec["recipe_M_knot_factor"], "recipe_M_knot_value_unused": rec["recipe_M_knot_value_unused"],
        "card_log": {"file": card_log_name, "key": "M_knot_rule", "sha256": _mknot_record_sha256(rec),
                     "sha256_of": "card_log の M_knot_rule を JSON (キーの順に並べ、区切りに空白なし、ensure_ascii=False) にしたもの",
                     "holds": "掃引の表 (table)・ゲートの全部の記録 (行・理由)・選び方の記録・選び直しの文"},
    })


def mknot_card_log_check(problem_yaml) -> dict:
    """問題 YAML (または Euler の写し) の `_card.M_knot_rule` (要約) を、同じディレクトリの card_log の `M_knot_rule` (全部) と照合する:
    記録の sha256 が要約の値と同じ、かつ要約が card_log の記録から作り直したものと同じ (summary_matches は sha256 の値を除いて比べる:
    card_log の表だけが違えば sha256 だけが合わない)。戻り: {ok (自動でなければ None), card_log, sha256_yaml, sha256_card_log,
    summary_matches}。"""
    p = Path(problem_yaml)
    summ = ((yaml.safe_load(p.read_text()) or {}).get("_card") or {}).get("M_knot_rule")
    if summ is None or summ.get("mode") != "auto":
        return {"ok": None, "card_log": None, "reason": "M_knot が自動の記録ではない (要約・card_log の照合の対象外)"}
    cl = p.with_name(summ["card_log"]["file"])
    full = json.loads(cl.read_text())["M_knot_rule"]
    sha = _mknot_record_sha256(full)

    def wo_sha(d):
        return {**d, "card_log": {**d["card_log"], "sha256": None}}

    same = wo_sha(_mknot_yaml_summary(full, summ["card_log"]["file"])) == wo_sha(summ)
    return {"ok": sha == summ["card_log"]["sha256"] and same, "card_log": str(cl), "sha256_yaml": summ["card_log"]["sha256"],
            "sha256_card_log": sha, "summary_matches": same}


def _mknot_record(card: dict, prob: dict, chain: dict, rt: float, mk: dict) -> dict:
    """`_card.M_knot_rule` (と card_log の M_knot_rule) の記録 (plan tooling-nozzle-mknot-auto §4.4): 規則の名前・掃引の表・選んだ値
    (chosen_N)・境目 (全部)・距離・解像度のゲートの結果・r_t の固定点の前後の選び直し。雛形の M_knot_factor は記録だけ。
    出力する組 (YAML の geometry.M_knot と r_t) の最終の設計チェーンが、ゲートがその組を評価した行と同じ設計か (θ_max が相対 1e-9 で
    一致するか) を確かめる (古い検証を流用していないこと)。違えば RuntimeError。"""
    from ..evaluate.mknot_select import RULE_NAME, _key
    meta = load_recipe(card["recipe"])[0]["_recipe"]
    unused = {"recipe_M_knot_factor": float(meta["M_knot_factor"]),
              "recipe_M_knot_value_unused": float(meta["M_knot_factor"]) * float(card["M_exit"]),
              "note": "雛形の M_knot_factor は使わない (記録だけ; plan tooling-nozzle-mknot-auto §4.4)"}
    g_law = prob["geometry"]
    if "M_knot" not in g_law:
        return _plain({"rule": RULE_NAME, "mode": "not_applicable", "source": card["M_knot_auto_source"],
                       "reason": f"axis_law {g_law.get('axis_law', 'quintic')!r} は knot でない — M_knot を使わない "
                                 "(geometry.M_knot を書かない・探索しない)", **unused})
    g, g1, sel = mk["gate_final"], mk["gate_first"], mk["sel"]
    M_out = float(g_law["M_knot"])
    row_N, row_2N = g["row_chosen_N"], (g["rows_selected"].get(f"{M_out:.12g}") or {}).get("at_2N")
    th_chain = float(chain["qa"]["theta_max_deg"])
    if (M_out != g["chosen_N"] or M_out != mk["M_K"] or row_N is None or row_2N is None or _key(row_N["M_K"]) != _key(M_out)
            or not row_N.get("ok") or not row_2N.get("ok")
            or abs(th_chain - float(row_N["theta_max_deg"])) > 1e-9 * abs(th_chain)):
        raise RuntimeError(f"M_knot: 出力する組 (geometry.M_knot {M_out!r}、r_t {rt!r} m) がゲートの chosen_N {g['chosen_N']!r} の評価と"
                           f"一致しない (最終の設計チェーンの θ_max {th_chain!r}、ゲートの行 "
                           f"{None if row_N is None else row_N.get('theta_max_deg')!r}・成立 N {None if row_N is None else row_N.get('ok')}・"
                           f"2N {None if row_2N is None else row_2N.get('ok')}) — 古い検証を流用しない (§4.3)")
    sel_keys = {_key(r["M_K"]) for r in sel["table"]}
    table = {"select": [_mknot_row(r) for r in sel["table"]],
             "gate_first_N_added": [_mknot_row(r) for r in g1["rows_N"] if _key(r["M_K"]) not in sel_keys],
             "gate_first_2N": [_mknot_row(r) for r in g1["rows_2N"]]}
    if g is not g1:
        table.update(gate_final_N=[_mknot_row(r) for r in g["rows_N"]], gate_final_2N=[_mknot_row(r) for r in g["rows_2N"]])
    rec = {
        "rule": RULE_NAME, "mode": "auto", "source": card["M_knot_auto_source"],
        "chosen": M_out, "chosen_from": "最終の組の解像度のゲートで選び直した N の値 chosen_N (plan §4.2・§4.3)",
        "r_t_m": float(rt), "first_r_t_m": mk["sel_rt"], "n_axis_inv": mk["N"],
        "theta_max_deg": row_N["theta_max_deg"], "margin_deg": row_N.get("margin_deg"), "min_sin_angle": row_N.get("min_sin_angle"),
        "mpp_sign_changes": row_N.get("mpp_sign_changes"), "x_K": row_N.get("x_K"), "moc_gate_ok": row_N.get("moc_gate_ok"),
        "theta_max_deg_final_chain": th_chain,
        "distance_to_boundary": g["distance_to_conservative_boundary"], "dM_safe": sel["settings"]["dM_safe"],
        "boundaries": _mknot_bounds(g["common_boundaries"]),
        "constraints_at_output": {"N": _mknot_row(row_N), "2N": _mknot_row(row_2N)},
        "gate": _mknot_gate_record(g), "gate_first": (None if g is g1 else _mknot_gate_record(g1)),
        "gate_first_reused": g is g1, "selection": _mknot_select_record(sel),
        "events": mk["events"], "redo": mk["redo"], "workers": mk["workers"], "elapsed_s": mk["t"],
        **unused, "table": table}
    return _plain(rec)


def build_problem(card_path, out_yaml, sizing: str = "auto", max_outer: int = 10, rel_tol: float = 1e-6,
                  r_throat_m: float | None = None, ns_prev_run=None, ns_delta_next=None, log=print) -> dict:
    r"""要求カード → 問題 YAML (§4.1 2.)。戻り: card_log.json と同じ dict。

    - 物理スロート径のカード (sizing="auto"): r_t = 径/2 で YAML を書く → `solve_rt_throat(YAML, 径/2)` (CFD 前の δ) → r_t を更新して
      YAML を書き直す、を $|\Delta r_t|/r_t \le$ rel_tol まで繰り返す (max_outer 回で止まらなければ RuntimeError)。内側の
      `solve_rt_throat` の許容差は min(CFD 前の既定 1e-7 m, rel_tol × 径/2) (外側の判定より緩い内側で止めない)。
    - **NS 後の補正** (ns_prev_run、CLI の `--ns-correction PREV_RUN [--delta-next CSV]`; plan tooling-nozzle-deltastar-band-root-zone
      §4.3 3): 上の固定点の内側を `solve_rt_throat(YAML, 径/2, prev_run=PREV_RUN, delta_next=CSV)` (前の run の δ_E の全分布 ×
      (r_t/S_prev)^−0.2、前の設計の x → s → その反復の設計の x) で解く。内側の許容差は min(NS 後の既定 1e-9 m, rel_tol × 径/2)。
      実寸の長さはカードのまま (固定点は同じ)。最終の r_t・最終の設計での補正の表を YAML の横の `<stem>_delta_r.csv` に書き
      (`deltastar_loop.write_ns_correction_table`、寸法を解いたのと同じ関数)、`_card.ns_correction` に出典 (前の run・delta_next・
      壁ファイル・prepare_info の sha256・S_prev・揃えた表の端の差・写し) と表のパス・sha256 を残す。NS は
      `deltastar_loop --delta-r-csv <表>` で立てる (sha256 を照合する)。物理スロート径のカード・sizing="auto"・r_throat_m なしに限る。
    - 設計スロート径のカード: r_t = 径/2 で 1 回。
    - r_throat_m (CLI の --r-throat): 外から与えた r_t で書くだけ (NS 後の解き直しの割り直し)。設計スロート径のカードには使えない。
    - sizing="none" (CLI の --no-solve): 物理スロート径でも解かずに r_t = 径/2 で書く (r_t は未解決とコメントに書く)。
    最後に設計チェーン (`design_chain`) を通し、`estimate_mesh`・`liquefaction_check`・設計壁から測る寸法を YAML のコメントと `_card` に
    書く。寸法が公差外でも書き出しは止めない (戻りの dimensions.pass が False、CLI は終了コード 2)。
    Euler の格子 (ni は明示が無ければ mesh.ni の写し) の節点数も estimate_mesh と同じ上限で判定し、超えたら書き出さずに ValueError。
    NS・Euler の節点数と上限は YAML のコメント・標準出力・card_log の `nodes` に出す (plan tooling-nozzle-mesh-size-guard §4.2)。

    自動調整 (どれも `_card.auto_adjust`・card_log・YAML のコメントに記録し、調整したら r_t の反復から全部やり直す):
    - **from_length が止まる** (STALL_TEXT): まず `n_axis_inv` を 2 倍にする。それでも止まれば、内部許容差 `geometry.L_total_tol` を
      L_TOTAL_TOL_RELAXED (0.05 r_t) に緩め、`n_axis_inv` は雛形の値に戻す (`L_total_tol_relaxed: {from, to, reason}`)。**カードの公差は
      変えない** — 寸法の合否は設計壁から測る検査 (ここ) と成果物の出荷ゲート (STEP・メッシュ) がカードの公差で決める。0.05 r_t でも
      止まれば例外 (形状が成立しない; codex 諮問 2026-10-13)。
    - **x_F が MOC の場の終点に張り付く** (`moc_call_record` の pinned、または終端特性線が壁に届かない例外): `geometry.x_end_margin` を
      0.5 ずつ上げる (上限 5.0、超えても解消しなければ例外)。

    M_knot が自動のカード (plan tooling-nozzle-mknot-auto §4.3; knot の軸則のとき):
    - 最初の r_t (径/2、`--r-throat` ならその値) で `select_mknot` → `resolution_gate`、出力の `chosen_N` を r_t の固定点の M_K に固定する。
      上の自動調整 (n_axis_inv の倍化・L_total_tol の緩和・x_end_margin の変更) が起きたら、その geometry で最初の r_t から選び直す。
    - 収束した最終の組 (r_t, M_K, 実効の geometry) で詰め (±0.1、0.01 刻み) と解像度のゲートを回す (最初の r_t と違えば N の評価は
      ±(0.1 + ΔM_safe) に広げて境目を求める)。`chosen_N` と固定点に使った M_K の差 > 0.01 なら chosen_N で r_t の固定点を前の r_t から
      やり直す (最大 MKNOT_REDO_MAX 回、超えたら RuntimeError)。差 ≤ 0.01 で値が違えば、出力する組 (chosen_N、最後の r_t) で r_t の固定点の
      1 回を回し直して寸法の残差を確かめる (|Δr_t|/r_t ≤ rel_tol でなければやり直しに数える)。出力する組の最終の設計チェーンの θ_max が
      ゲートのその組の評価と一致することも確かめる (`_mknot_record`)。
    - ゲートで選び直した値 (同点の規則の後) が詰めの窓の端なら、同じ組のまま窓を半幅 0.1 → 0.2 → 0.3 に広げて選び直す
      (`MKNOT_GATE_HALF_WIDTHS`; 上限でも端なら MknotUndeterminedError)。
    - 記録は card_log の `M_knot_rule` に全部 (掃引の表を含む)、YAML の `_card.M_knot_rule` (と Euler の写し) には要約だけ
      (`_mknot_yaml_summary`; card_log の記録の sha256 で結び、書いた後に `mknot_card_log_check` で照合する)。
    - 評価の並列は `mknot_workers` (MKNOT_WORKERS: N・2N とも 1)。2N の 1 評価のメモリの見積もりが使えるメモリを超えるなら評価の前に
      RuntimeError (`_mknot_mem_check`)。"""
    from ..feedback.deltastar_loop import SIZING_TOL_M, SIZING_TOL_PRE_CFD_THROAT_M, load_ns_correction, solve_rt_throat
    from ..probdef import load_problem
    t_all = time.time()
    card_path, out = Path(card_path), Path(out_yaml)
    card = load_card(card_path)
    card_sha = _sha256(card_path)
    kind, D = throat_kind(card)
    if sizing not in ("auto", "none"):
        raise ValueError(f"sizing は 'auto' か 'none' (受け取った値: {sizing!r})")
    if r_throat_m is not None and kind == "design":
        raise ValueError("設計スロート径 (diameter_design_m) のカードに --r-throat は使えない (r_t = 径/2 で決まる)")
    if ns_delta_next is not None and ns_prev_run is None:
        raise ValueError("--delta-next は --ns-correction PREV_RUN と一緒に使う")
    if ns_prev_run is not None:
        if kind != "physical" or sizing != "auto" or r_throat_m is not None:
            raise ValueError("--ns-correction は物理スロート径のカードで r_t を解くときだけ (--r-throat・--no-solve・設計スロート径とは併用しない)")
        load_ns_correction(ns_prev_run, ns_delta_next)   # 入力の検査 (使用域の held・表の端・前の設計の壁ファイル) を固定点の前に
    auto: dict = {}                  # 自動調整した geometry のキー → 値 (make が card_to_problem の後に当てる)
    auto_adjust: list = []           # 自動調整の履歴
    relaxed: dict | None = None      # L_total_tol_relaxed
    moc_log: list = []               # design_chain ごとの from_length の各評価 (card_log)
    last_inner_table: bytes | None = None   # NS 後の補正: 最後の内側の solve_rt_throat が書いた補正の表 (最終の表との照合用)
    # M_knot が自動のカード (plan tooling-nozzle-mknot-auto §4.3): 選ぶ前の問題 dict の geometry.M_knot は NaN の仮置き (評価
    # evaluate_mknot は M_knot を差し替える; 設計チェーン・寸法の計算には選んだ値が入るまで渡さない)
    mk_auto = card["M_knot_mode"] == "auto"
    base_geom = card_to_problem(card, D / 2.0, **({"M_knot": float("nan")} if mk_auto else {}))["geometry"]
    mk_knot = mk_auto and str(base_geom.get("axis_law", "quintic")) == "knot"
    mk: dict = {"M_K": None, "key": None, "events": [], "attempt": 0, "aa_seen": 0, "t": 0.0}

    def make(rt, M_K=None):
        if mk_auto:
            v = M_K if M_K is not None else mk["M_K"]
            prob = card_to_problem(card, rt, M_knot=(float("nan") if v is None else v))
        else:
            prob = card_to_problem(card, rt)
        prob["geometry"].update(copy.deepcopy(auto))
        return prob

    def cur(key, rt):
        return make(rt)["geometry"].get(key)

    def on_stall(ex, rt, calls):
        nonlocal relaxed
        L_tgt = float(make(rt)["dv"]["L_total"]["value"])
        best = min((abs(c["x_F"] - L_tgt) for c in calls if np.isfinite(c["x_F"])), default=float("nan"))
        tol_now = float(cur("L_total_tol", rt))
        info = f"{STALL_TEXT}: 最良の |x_F − L_total| {best:.3g} r_t > L_total_tol {tol_now:.3g} r_t (r_t {rt!r} m)"
        if "n_axis_inv" not in auto and relaxed is None:
            n0 = int(base_geom["n_axis_inv"])
            auto["n_axis_inv"] = 2 * n0
            auto_adjust.append({"key": "geometry.n_axis_inv", "from": n0, "to": 2 * n0, "reason": info})
            log(f"[wt_card] {STALL_TEXT} — n_axis_inv を {n0} → {2 * n0} にしてやり直す", flush=True)
            return
        if relaxed is None:
            tol_card = float(card["tolerances"]["length_m"]) / 2.0 / rt
            if tol_card >= L_TOTAL_TOL_RELAXED:
                raise ValueError(f"{info}。L_total_tol はすでに {L_TOTAL_TOL_RELAXED} r_t 以上なので緩めない — 形状が成立しない "
                                 "(形状パラメータを overrides で振る)") from ex
            n_from = auto.pop("n_axis_inv", None)
            auto["L_total_tol"] = L_TOTAL_TOL_RELAXED
            relaxed = {"from": tol_card, "from_text": "カードの長さ公差 / 2 / r_t", "to": L_TOTAL_TOL_RELAXED,
                       "reason": info + (f" (n_axis_inv {int(base_geom['n_axis_inv'])} と {n_from} の両方)" if n_from else ""),
                       "n_axis_inv": int(base_geom["n_axis_inv"]),
                       "note": "カードの公差は変えていない。寸法の合否は設計壁から測る検査と成果物の出荷ゲートがカードの公差で決める"}
            auto_adjust.append({"key": "geometry.L_total_tol", "from": tol_card, "to": L_TOTAL_TOL_RELAXED, "reason": relaxed["reason"]})
            if n_from:
                auto_adjust.append({"key": "geometry.n_axis_inv", "from": n_from, "to": int(base_geom["n_axis_inv"]),
                                    "reason": "内部許容差を緩めたので雛形の値に戻した"})
            log(f"[wt_card] {STALL_TEXT} — 内部許容差 L_total_tol を {tol_card:.3g} → {L_TOTAL_TOL_RELAXED} r_t に緩めてやり直す "
                "(カードの公差は変えない)", flush=True)
            return
        raise ValueError(f"{info}。内部許容差を {L_TOTAL_TOL_RELAXED} r_t に緩めても design_chain が通らない — 形状が成立しない "
                         "(形状パラメータを overrides で振る)") from ex

    def bump_margin(why, rt):
        m0 = float(cur("x_end_margin", rt) if cur("x_end_margin", rt) is not None else 2.3)
        m1 = round(m0 + X_END_MARGIN_STEP, 10)
        if m1 > X_END_MARGIN_MAX + 1e-12:
            raise RuntimeError(f"x_F が逆 MOC の場の終点に張り付く ({why}) — geometry.x_end_margin を上限 {X_END_MARGIN_MAX} まで上げても"
                               "解消しない")
        auto["x_end_margin"] = m1
        auto_adjust.append({"key": "geometry.x_end_margin", "from": m0, "to": m1, "reason": why})
        log(f"[wt_card] x_F が逆 MOC の場の終点に張り付く ({why}) — x_end_margin を {m0} → {m1} にしてやり直す", flush=True)

    def chain_logged(prob, stage, rt):
        """design_chain を逆 MOC の記録つきで回す。from_length の停止・場の不足は自動調整して _Retry。"""
        cap = _MocCapture()
        L_tgt = float(prob["dv"]["L_total"]["value"])
        g = prob["geometry"]
        rec = {"stage": stage, "r_t_m": rt, "n_axis_inv": g.get("n_axis_inv"), "x_end_margin": g.get("x_end_margin"),
               "L_total_tol": g.get("L_total_tol"), "L_total_target_rt": L_tgt}
        try:
            chain = _design_chain_of(prob, capture=cap)
        except (ValueError, RuntimeError) as ex:
            rec.update(evals=_evals_log(cap.calls, rt, L_tgt), outcome=f"{type(ex).__name__}: {str(ex)[:200]}")
            moc_log.append(rec)
            if _is_from_length_stall(ex):
                on_stall(ex, rt, cap.calls)
                raise _Retry() from ex
            if _is_field_short(ex):
                bump_margin(f"終端特性線が壁に届かない ({stage})", rt)
                raise _Retry() from ex
            raise
        field = _chain_field(chain, cap.calls)
        rec.update(evals=_evals_log(cap.calls, rt, L_tgt), outcome="ok", L_c=chain["L_c"], field=field)
        moc_log.append(rec)
        return chain, field

    def inner_tol_of(ns: bool) -> float:
        # 内側の許容差: CFD 前は初期見積もりの既定 (1e-7 m)、NS 後の補正は NS 後の既定 (1e-9 m) と rel_tol × 径/2 の小さい方
        return min(SIZING_TOL_M if ns else SIZING_TOL_PRE_CFD_THROAT_M, rel_tol * D / 2.0)

    def outer_step(k: int, rt: float, ns: bool, inner_tol: float, stage: str | None = None) -> tuple:
        """r_t の固定点の 1 回: r_t で YAML を書く → solve_rt_throat。戻り (r_t の出力, |Δr_t|/r_t, solve_rt_throat の結果, 履歴の行)。"""
        nonlocal last_inner_table
        t0 = time.time()
        prob = make(rt)
        cap = _MocCapture()
        L_tgt = float(prob["dv"]["L_total"]["value"])
        rec = {"stage": stage or f"r_t 反復 {k} (solve_rt_throat{', NS 後の補正' if ns else ''})", "r_t_m": rt,
               "n_axis_inv": prob["geometry"].get("n_axis_inv"),
               "x_end_margin": prob["geometry"].get("x_end_margin"), "L_total_tol": prob["geometry"].get("L_total_tol"),
               "L_total_target_rt": L_tgt}
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / f"{card['name']}_iter{k}.yaml"
            path.write_text(yaml.safe_dump(prob, sort_keys=False, allow_unicode=True))
            try:
                with cap:
                    if ns:
                        tab = Path(td) / "delta_r_iter.csv"
                        res = solve_rt_throat(path, D / 2.0, prev_run=ns_prev_run, delta_next=ns_delta_next, delta_r_out=tab,
                                              tol_R_m=inner_tol)
                        last_inner_table = tab.read_bytes()
                    else:
                        res = solve_rt_throat(path, D / 2.0, tol_R_m=inner_tol)
            except (ValueError, RuntimeError) as ex:
                rec.update(evals=_evals_log(cap.calls, rt, L_tgt), outcome=f"{type(ex).__name__}: {str(ex)[:200]}")
                moc_log.append(rec)
                if _is_from_length_stall(ex):
                    on_stall(ex, rt, cap.calls)
                    raise _Retry() from ex
                if _is_field_short(ex):
                    bump_margin(f"終端特性線が壁に届かない ({stage or f'r_t 反復 {k}'})", rt)
                    raise _Retry() from ex
                raise
        rec.update(evals=_evals_log(cap.calls, rt, L_tgt), outcome="ok")
        moc_log.append(rec)
        rt_new = float(res["r_t_m"])
        rel = abs(rt_new - rt) / rt
        h = {"iter": k, "r_t_in_m": rt, "r_t_out_m": rt_new, "rel_change": rel,
             "inner_residual_m": float(res["residual_m"]), "inner_n_iter": int(res["n_iter"]),
             "physical_throat_r_m": float(res["physical_throat"]["r_m"]),
             "physical_throat_x_m": float(res["physical_throat"]["x_m"]),
             "exit_radius_m": float(res["exit_radius_m"]), "time_s": time.time() - t0}
        if mk_knot:                      # M_knot が自動のときだけ: その反復の M_K と段 (数値の M_knot のカードの記録は今と同じ)
            h.update(M_K=mk["M_K"], stage=stage or "r_t の固定点")
        log(f"[wt_card] {stage or f'r_t 反復 {k}'}: r_t {rt!r} → {rt_new!r} m (|Δ|/r_t {rel:.3e}、内側 {res['n_iter']} 回・残差 "
            f"{res['residual_m']:.3e} m、{time.time() - t0:.1f} s)", flush=True)
        return rt_new, rel, res, h

    def solve_loop(rt_start: float | None = None):
        """r_t を決める。戻り (rt, 履歴, solve_rt_throat の最後の結果, method, method_text, 内側の許容差)。
        rt_start: 固定点の始まりの r_t (既定は径/2; M_knot の選び直しの後のやり直しは前の r_t から)。"""
        if r_throat_m is not None:
            rt = float(r_throat_m)
            return rt, [], None, "external", f"外から与えた r_t (--r-throat {rt!r} m; NS 後の解き直しなど。生成器は解いていない)", None
        if kind == "design":
            return D / 2.0, [], None, "design_throat", "設計スロート径 / 2 (解かない)", None
        if sizing == "none":
            return (D / 2.0, [], None, "unsolved", "**未解決** (--no-solve): 仮に 物理スロート径/2。solve_rt_throat で解いていない", None)
        if mk_knot and not (mk["M_K"] is not None and math.isfinite(mk["M_K"])):
            raise RuntimeError("build_problem: M_knot を選ぶ前に r_t の固定点に入った (内部の誤り)")
        ns = ns_prev_run is not None
        inner_tol = inner_tol_of(ns)
        rt = D / 2.0 if rt_start is None else float(rt_start)
        hist: list = []
        res = None
        for k in range(int(max_outer)):
            rt_new, rel, res, h = outer_step(k, rt, ns, inner_tol)
            hist.append(h)
            rt = rt_new
            if rel <= rel_tol and ns:
                text = (f"NS 後の補正の solve_rt_throat (δ: 前の run {ns_prev_run} の δ_E の全分布 × (r_t/S_prev)^−0.2、前の設計の x → s → "
                        f"今の設計の x) で物理スロート半径 {D / 2.0!r} m に合わせた固定点 (外側 {len(hist)} 回、最後の |Δr_t|/r_t {rel:.2e} ≤ "
                        f"{rel_tol:g}、内側の許容差 {inner_tol:.3g} m)。補正の表は YAML の横の {out.stem}_delta_r.csv "
                        "(出典と sha256 は _card.ns_correction; NS は deltastar_loop --delta-r-csv で立てる)")
                return rt, hist, res, "physical_throat_ns_corrected", text, inner_tol
            if rel <= rel_tol:
                text = (f"CFD 前の solve_rt_throat (δ: 積分法の初期壁、prepare_ns と同じ経路) で物理スロート半径 {D / 2.0!r} m に合わせた "
                        f"固定点 (外側 {len(hist)} 回、最後の |Δr_t|/r_t {rel:.2e} ≤ {rel_tol:g}、内側の許容差 {inner_tol:.3g} m)。"
                        "CFD 前の寸法は初期見積もり — NS 後に deltastar_loop --solve-rt-throat --prev で解き直し、"
                        "その r_t を生成器の --r-throat に渡して長さを割り直す")
                return rt, hist, res, "physical_throat_solved", text, inner_tol
        raise RuntimeError(f"r_t の固定点が {max_outer} 回で |Δr_t|/r_t ≤ {rel_tol:g} に止まらない (履歴 "
                           f"{[(h['r_t_in_m'], h['r_t_out_m'], h['rel_change']) for h in hist]})")

    # ------------------------------------------------------------------ M_knot の自動の選択 (plan tooling-nozzle-mknot-auto §4.3)
    def mlog(msg):
        log(f"[wt_card] M_knot: {msg}", flush=True)

    def mk_key(rt):
        return (float(rt), tuple(sorted((k, repr(v)) for k, v in auto.items())))

    def mk_base(rt, M_K=None) -> dict:
        prob = make(rt, M_K=M_K)
        return {k: v for k, v in prob.items() if k != "_card"}

    def mk_geom(rt) -> dict:
        g = make(rt)["geometry"]
        return {k: g.get(k) for k in ("n_axis_inv", "x_end_margin", "L_total_tol")}

    def mk_event(stage, rt, g, action, t0, M_used=None):
        mk["events"].append(_plain({
            "stage": stage, "attempt": mk["attempt"], "r_t_m": rt, "geometry": mk_geom(rt), "M_K_in": g["chosen"],
            "M_K_used_in_fixed_point": M_used, "chosen_N": g["chosen_N"], "chosen_2N": g["chosen_2N"],
            "diff_vs_used": (None if M_used is None else g["chosen_N"] - M_used), "gate_pass": g["pass"],
            "gate_n_widen": g["n_widen"], "gate_half_width": g["half_width"],
            "distance": g["distance_to_conservative_boundary"], "action": action, "elapsed_s": time.time() - t0}))

    def mk_select_first(rt0):
        """最初の r_t で select_mknot → resolution_gate、出力は chosen_N (§4.2・§4.3)。r_t の固定点の間はこの値で固定する。"""
        from ..evaluate import mknot_select as MS
        t0 = time.time()
        base = mk_base(rt0, M_K=float("nan"))          # M_knot は仮置き (search_range は外し、評価は差し替える)
        N = int(base["geometry"].get("n_axis_inv", 500))
        wk = mknot_workers(N)
        _mknot_mem_check(wk, N)
        rng = MS.search_range(base)
        new_aa = auto_adjust[mk["aa_seen"]:]
        mk["aa_seen"] = len(auto_adjust)
        stage = ("最初の r_t" if mk["key"] is None else
                 "自動調整の後の選び直し (" + "; ".join(f"{a['key']} {a['from']!r} → {a['to']!r}" for a in new_aa) + ")")
        mlog(f"{stage}: r_t {rt0!r} m・n_axis_inv {N} で探索 (M_A {rng['M_A']:.6g}〜M_d {rng['M_d']:.6g}、並列 N {wk['N']}・2N {wk['2N']}、"
             f"使えるメモリ {wk['available_gb'] if wk['available_gb'] is None else round(wk['available_gb'], 1)} GB)")
        sel = MS.select_mknot(MS.make_eval_fn(base, N), rng["M_A"], rng["M_d"], workers=wk["N"], log=mlog)
        bad = [r["M_K"] for r in sel["table"] if r.get("M_A") is not None and r["M_A"] != rng["M_A"]]
        if bad:
            raise RuntimeError(f"M_knot: 評価の行の M_A が探索の原点 {rng['M_A']!r} と違う ({bad[:5]}) — 格子の原点が揃っていない")
        gate = MS.resolution_gate(base, sel["chosen"], N, select_result=sel, half_widths=MKNOT_GATE_HALF_WIDTHS,
                                  workers=(wk["N"], wk["2N"]), log=mlog)
        mk.update(M_K=float(gate["chosen_N"]), key=mk_key(rt0), sel=sel, sel_rt=float(rt0), gate_first=gate, gate_final=None,
                  N=N, workers=wk, redo=0, seedN={MS._key(r["M_K"]): r for r in gate["rows_N"]},
                  seed2={MS._key(r["M_K"]): r for r in gate["rows_2N"]})
        mk_event(stage, rt0, gate, f"select_mknot の選択 {sel['chosen']!r} → ゲートで選び直した chosen_N {gate['chosen_N']!r} を "
                 "r_t の固定点の M_K に固定", t0)
        mlog(f"{stage}: 選択 {sel['chosen']!r} → chosen_N {gate['chosen_N']!r} (ゲート PASS、全部の境目から "
             f"{_mknot_dist(gate['distance_to_conservative_boundary'])}"
             + (f"、{_mknot_widen_text(gate)}" if gate["n_widen"] else "") + f"、{time.time() - t0:.0f} s)")
        mk["t"] += time.time() - t0

    def mk_gate(rt, M_c):
        """(rt, M_c) を中心に詰めと解像度のゲート。同じ r_t・geometry の探索の表があれば再利用し、無ければ N の評価を
        ±(最初の半幅 + ΔM_safe) に広げて境目を求める (窓の候補から ΔM_safe 以内の境目を全部見るため)。ゲートが窓を広げたときの
        窓の外の文脈 ±(広げた半幅 + ΔM_safe) は resolution_gate が N で評価する。"""
        from ..evaluate import mknot_select as MS
        base = mk_base(rt, M_K=M_c)
        N, sel = mk["N"], mk["sel"]
        st = sel["settings"]
        wk = mknot_workers(N)
        _mknot_mem_check(wk, N)
        if mk_key(rt) == mk["key"]:                     # 同じ r_t・geometry: 最初の探索とゲートの評価 (N・2N) を使う
            ctx = {"settings": st, "table": list(mk["seedN"].values())}
            seed2, how = list(mk["seed2"].values()), "最初の r_t と同じ組 (探索とゲートの評価を再利用)"
        else:
            n = int(round((MKNOT_GATE_HALF_WIDTH + st["dM_safe"]) / st["fine"]))
            Ms = [round(M_c + j * st["fine"], 12) for j in range(-n, n + 1)]
            mlog(f"最終の組 r_t {rt!r} m: M_K {M_c!r} の前後 ±{MKNOT_GATE_HALF_WIDTH + st['dM_safe']:g} を n_axis_inv {N} で評価 "
                 f"({len(Ms)} 点、並列 {wk['N']})")
            with MS.Evaluator(MS.make_eval_fn(base, N), wk["N"], mlog) as ev:
                ev.many(Ms, "context")
                ctx = {"settings": st, "table": ev.table()}
            seed2, how = None, "最終の組で詰め直し (N は ±(半幅 + ΔM_safe))"
        g = MS.resolution_gate(base, M_c, N, select_result=ctx, seed_2N=seed2, half_widths=MKNOT_GATE_HALF_WIDTHS,
                               workers=(wk["N"], wk["2N"]), log=mlog)
        if mk_key(rt) == mk["key"]:
            mk["seedN"].update({MS._key(r["M_K"]): r for r in g["rows_N"]})
            mk["seed2"].update({MS._key(r["M_K"]): r for r in g["rows_2N"]})
        return g, how

    def mk_finalize(state):
        """収束した最終の組 (r_t, M_K, 実効の geometry) で詰めと解像度のゲート (§4.3)。chosen_N と r_t の固定点に使った M_K の差が
        0.01 を超えたら chosen_N で r_t の固定点をやり直す (最大 MKNOT_REDO_MAX 回)。差が 0.01 以下で値が違えば、出力する組
        (chosen_N と最後の r_t) で r_t の固定点の 1 回 (寸法の残差) を回し直して確かめる (制約はゲートがその組を N・2N で評価済み)。
        戻り: state (rt, hist, sizing_res, method, method_text, inner_tol)。出力する組は (mk["M_K"], rt)。"""
        rt, hist, sizing_res, method, method_text, inner_tol = state
        solving = method in ("physical_throat_solved", "physical_throat_ns_corrected")
        ns = method == "physical_throat_ns_corrected"
        while True:
            t0 = time.time()
            M_used = mk["M_K"]
            if mk_key(rt) == mk["key"] and M_used == mk["gate_first"]["chosen"]:
                g, how = mk["gate_first"], "最初の r_t のゲートと同じ組 (再利用)"
            else:
                g, how = mk_gate(rt, M_used)
            mk["gate_final"] = g
            cN = float(g["chosen_N"])
            diff = cN - M_used
            stage = f"最終の組 ({how})"
            if abs(diff) > MKNOT_MOVE_TOL + 1e-9:
                mk["redo"] += 1
                if mk["redo"] > MKNOT_REDO_MAX:
                    mk_event(stage, rt, g, f"chosen_N と固定点の M_K の差 {diff:+.4g} > {MKNOT_MOVE_TOL:g} が {MKNOT_REDO_MAX} 回の"
                             "やり直しでも収まらない — 例外", t0, M_used)
                    raise RuntimeError(f"M_knot: 最終の組のゲートの chosen_N {cN!r} が r_t の固定点に使った M_K {M_used!r} から "
                                       f"{diff:+.4g} 動き、r_t の固定点を {MKNOT_REDO_MAX} 回やり直しても {MKNOT_MOVE_TOL:g} 以内に収まらない "
                                       f"(plan tooling-nozzle-mknot-auto §4.3; 記録 {mk['events']})")
                mk["M_K"] = cN
                mk_event(stage, rt, g, f"差 {diff:+.4g} > {MKNOT_MOVE_TOL:g} → chosen_N {cN!r} で r_t の固定点をやり直す "
                         f"({mk['redo']} 回目、始まりは前の r_t)", t0, M_used)
                mlog(f"最終の組: chosen_N {cN!r} が固定点の M_K {M_used!r} から {diff:+.4g} — r_t の固定点をやり直す ({mk['redo']} 回目)")
                mk["t"] += time.time() - t0
                if solving:
                    rt, h2, sizing_res, method, method_text, inner_tol = solve_loop(rt_start=rt)
                    hist = hist + h2
                continue
            if diff == 0.0:
                mk_event(stage, rt, g, "chosen_N = r_t の固定点に使った M_K — その組を出力 (寸法は最終の設計チェーンで測る)", t0, M_used)
                mk["t"] += time.time() - t0
                return rt, hist, sizing_res, method, method_text, inner_tol
            mk["M_K"] = cN
            if not solving:
                mk_event(stage, rt, g, f"差 {diff:+.4g} ≤ {MKNOT_MOVE_TOL:g}: chosen_N {cN!r} を出力 (r_t は固定; 制約はゲートがこの組を "
                         "N・2N で評価済み、寸法は最終の設計チェーンで測る)", t0, M_used)
                mk["t"] += time.time() - t0
                return rt, hist, sizing_res, method, method_text, inner_tol
            rt_v, rel, res_v, h = outer_step(len(hist), rt, ns, inner_tol if inner_tol is not None else inner_tol_of(ns),
                                             stage=f"M_knot の確かめ直し (chosen_N {cN!r}、r_t {rt!r} m)")
            hist = hist + [h]
            if rel <= rel_tol:
                sizing_res = res_v
                mk_event(stage, rt, g, f"差 {diff:+.4g} ≤ {MKNOT_MOVE_TOL:g}: 出力する組 (chosen_N {cN!r}、r_t {rt!r} m) で r_t の固定点の "
                         f"1 回を回し直した (|Δr_t|/r_t {rel:.2e} ≤ {rel_tol:g}; 寸法の残差はこの solve_rt_throat と最終の設計チェーンで測る)",
                         t0, M_used)
                mk["t"] += time.time() - t0
                return rt, hist, sizing_res, method, method_text, inner_tol
            mk["redo"] += 1
            if mk["redo"] > MKNOT_REDO_MAX:
                raise RuntimeError(f"M_knot: 出力する組 (chosen_N {cN!r}、r_t {rt!r} m) の r_t の固定点の 1 回で |Δr_t|/r_t {rel:.2e} > "
                                   f"{rel_tol:g}、やり直しが {MKNOT_REDO_MAX} 回を超えた")
            mk_event(stage, rt, g, f"差 {diff:+.4g} ≤ {MKNOT_MOVE_TOL:g} だが出力する組の r_t の固定点の 1 回で |Δr_t|/r_t {rel:.2e} > "
                     f"{rel_tol:g} — chosen_N で r_t の固定点をやり直す ({mk['redo']} 回目)", t0, M_used)
            mk["t"] += time.time() - t0
            rt, h2, sizing_res, method, method_text, inner_tol = solve_loop(rt_start=rt_v)
            hist = hist + h2

    for attempt in range(_MAX_ATTEMPTS):
        try:
            if mk_knot:
                mk["attempt"] = attempt
                rt0 = float(r_throat_m) if r_throat_m is not None else D / 2.0
                if mk["key"] != mk_key(rt0):           # 最初の r_t、または自動調整で geometry が変わった → 選び直す
                    mk_select_first(rt0)
            rt, hist, sizing_res, method, method_text, inner_tol = solve_loop()
            if mk_knot:
                rt, hist, sizing_res, method, method_text, inner_tol = mk_finalize(
                    (rt, hist, sizing_res, method, method_text, inner_tol))
            prob = make(rt)
            prob["_card"]["r_t_method"] = method
            prob["_card"]["r_t_method_text"] = method_text
            t0 = time.time()
            chain, field = chain_logged(prob, "最終", rt)
            t_chain = time.time() - t0
            if field["pinned"]:
                bump_margin(f"最終の x_F {field['x_F']:.6g} r_t と場の終点 {field['x_wall_full_end']:.6g} r_t の差 "
                            f"{field['margin_to_field_end_rt']:.3g} ≤ 点間隔 {field['spacing_end']:.3g}", rt)
                continue
            break
        except _Retry:
            continue
    else:
        raise RuntimeError(f"自動調整のやり直しが {_MAX_ATTEMPTS} 回を超えた ({auto_adjust})")
    if chain["Lc_solve"] is not None:
        log(f"[wt_card] design_chain {t_chain:.1f} s (L_c {chain['L_c']:.6g}、x_F − L_total {chain['Lc_solve']['xF_residual']:.3e} r_t、"
            f"x_F − x_e {field['x_F_minus_x_e_rt']:.3e} r_t、場の終点まで {field['margin_to_field_end_rt']:.3g} r_t)", flush=True)
    warnings: list = []
    if method == "unsolved":
        warnings.append("r_t: --no-solve のため物理スロート径から r_t を解いていない (仮に 径/2)")
    if relaxed is not None:
        warnings.append(f"from_length: 内部許容差 L_total_tol を {relaxed['from']:.3g} → {relaxed['to']} r_t に緩めた ({STALL_TEXT})。"
                        "カードの公差は変えていない — 寸法の合否を見る")
    est = estimate_mesh(card, prob, chain=chain)
    ov_keys = {o["key"] for o in prob["_card"]["overrides"]}
    for k in _MESH_ESTIMATED + ("wall_first_up_x0", "wall_first_up_x1"):
        if f"mesh.{k}" not in ov_keys:
            prob["mesh"][k] = est[k]
    if "mesh_euler.ni" not in ov_keys:
        prob["mesh_euler"]["ni"] = int(prob["mesh"]["ni"])
    # Euler の格子の節点数 (plan tooling-nozzle-mesh-size-guard §4.2; ni は明示が無ければ mesh.ni の写し)。prepare と同じ Mesh2DParams で数える
    from ..evaluate.runner_axismach import MESH_EULER_DEFAULTS, _mesh_params_from, check_node_budget
    mp_e = _mesh_params_from(prob["mesh_euler"], 1.0, **MESH_EULER_DEFAULTS)
    try:
        nodes_eu = check_node_budget(mp_e, prob["mesh_euler"], "mesh_euler", allow_hint=_card_allow_hint("mesh_euler"))
    except ValueError as ex:
        src = "overrides の mesh_euler.ni" if "mesh_euler.ni" in ov_keys else f"mesh.ni {prob['mesh']['ni']} の写し"
        raise ValueError(f"要求カード → 問題 (Euler の格子 ni {mp_e.ni} は {src}): {ex}") from ex
    nodes = {"mesh": est["nodes"], "mesh_euler": nodes_eu}
    warnings += est["warnings"]
    liq = liquefaction_check(card)
    if liq.get("warning"):
        warnings.append(liq["warning"])
    dims = measure_design_dimensions(card, chain, rt, sizing=sizing_res)
    for r in dims["rows"]:
        if not r["pass"]:
            warnings.append(f"寸法: {r['name']} の差 {r['diff_m']:+.3e} m が公差 {r['tol_m']:g} m を超える ({r['how']})")
    if card["tolerances_default"]:
        warnings.append(f"公差: カードに tolerances が無いので既定 (仮の値) {DEFAULT_TOLERANCES} を使った")
    field_rec = {k: field[k] for k in ("x_F", "x_e", "x_F_minus_x_e_rt", "x_wall_full_end", "spacing_end", "margin_to_field_end_rt",
                                       "n_beyond_x_F", "x_axis_end", "pinned")}
    field_rec["x_F_minus_x_e_m"] = field["x_F_minus_x_e_rt"] * rt
    field_rec["x_end_margin"] = prob["geometry"].get("x_end_margin")
    prob["_card"].update({
        "card_path": str(card_path), "card_sha256": card_sha, "card": _card_record(card),
        "r_t_history": hist, "inner_tol_m": inner_tol, "auto_adjust": auto_adjust,
        "L_total_tol_relaxed": relaxed,
        "moc_field": field_rec,
        "mesh_estimate": {k: est[k] for k in ("note", "wall_first_frac", "wall_first_frac_throat", "ni", "wall_first_frac_raw",
                                               "wall_first_frac_throat_raw", "ni_required_for_ar", "ar_max_limit",
                                               "y1plus_est_max_with_blend", "x_y1plus_est_max_rt", "blend_choice",
                                               "wall_first_up_x0", "wall_first_up_x1", "wall_temperature") if k in est},
        "mesh_estimate_ar_check": est["mesh_check"],
        "liquefaction": {k: liq[k] for k in ("verdict", "T_exit_K", "p_exit_Pa", "p_N2_Pa", "p_sat_N2_Pa", "supersaturated")},
        "dimensions": {"pass": dims["pass"], "rows": [{k: r[k] for k in ("name", "card_m", "measured_m", "diff_m", "tol_m", "pass")}
                                                      for r in dims["rows"]],
                       "x_F_minus_x_e_m": dims["x_F_minus_x_e_m"]},
        "warnings": list(warnings),
        "generated": datetime.date.today().isoformat(),
    })
    mk_rec = _mknot_record(card, prob, chain, rt, mk) if mk_auto else None
    card_log_name = out.name + ".card_log.json"
    if mk_rec is not None:          # YAML には要約だけ (掃引の表と全部の記録は card_log の M_knot_rule; plan §4.4)
        prob["_card"]["M_knot_rule"] = _mknot_yaml_summary(mk_rec, card_log_name)
    out.parent.mkdir(parents=True, exist_ok=True)
    ns_rec = None
    if ns_prev_run is not None:
        # 最終の r_t・最終の設計 (この YAML の設計壁) での補正の表を YAML の横に書く (寸法を解いたのと同じ関数; plan
        # tooling-nozzle-deltastar-band-root-zone §4.3 3)。最後の内側の反復は r_t の入力が最終の r_t と rel_tol 以内で違いうるので、
        # 表は最終の設計の区切りの点で作り直し、内側の表と同じかを記録する
        from ..feedback.deltastar_loop import (design_anchors, ns_correction_record, ns_correction_table,
                                               write_ns_correction_table)
        corr = load_ns_correction(ns_prev_run, ns_delta_next)
        a_new = design_anchors(chain["wall"])
        tb = ns_correction_table(corr, a_new, rt)
        table_path = out.with_name(out.stem + "_delta_r.csv")
        write_ns_correction_table(table_path, tb, corr, rt, a_new)
        ns_rec = _json_safe({
            **ns_correction_record(corr, a_new), "table": str(table_path), "table_sha256": _sha256(table_path),
            "r_t_m": rt, "factor": tb["factor"], "table_x_range_rt": [float(tb["x"][0]), float(tb["x"][-1])],
            "final_table_equals_last_inner": (last_inner_table is not None and table_path.read_bytes() == last_inner_table),
            "use": f"deltastar_loop --problem {out.name} --delta-r-csv {table_path.name} (sha256 を照合する)"})
        prob["_card"]["ns_correction"] = ns_rec
    header = _header(card, prob, card_path, card_sha, method_text, est, liq, dims, warnings, auto_adjust, nodes,
                     mk_full=mk_rec, card_log_name=card_log_name)
    _write_problem(card, prob, out, header)
    load_problem(out)               # 書いたものが読み込みの検査を通ること
    euler_out = write_euler_companion(out)
    mapping = mapping_table(card, prob)
    logd = {"generator": GENERATOR, "generated": datetime.datetime.now().isoformat(timespec="seconds"),
            "card_path": str(card_path), "card_sha256": card_sha, "out": str(out), "out_sha256": _sha256(out),
            "card": _card_record(card), "mapping": mapping, "warnings": warnings,
            "r_t": {"method": method, "text": method_text, "r_t_m": rt, "history": hist, "inner_tol_m": inner_tol,
                    "rel_tol": rel_tol, "max_outer": max_outer,
                    "sizing_result": (None if sizing_res is None else _json_safe({k: v for k, v in sizing_res.items()
                                                                                  if k not in ("iter_history",)}))},
            "auto_adjust": auto_adjust, "L_total_tol_relaxed": relaxed,
            "moc_field": field_rec, "moc_evals": moc_log,
            "design_chain": {"time_s": t_chain, "L_c": chain["L_c"], "Lc_window": chain["Lc_window"],
                             "Lc_solve": chain["Lc_solve"], "x_E": chain["x_E"], "M_knot": chain["M_knot"],
                             "qa": _json_safe(chain["qa"])},
            "mesh_estimate": {k: v for k, v in est.items() if k != "sections"},     # 断面ごとの配列は書かない
            "nodes": nodes,          # NS (mesh)・Euler (mesh_euler) の実効の節点数と上限 (plan tooling-nozzle-mesh-size-guard §4.2)
            "liquefaction": liq, "dimensions": dims, "time_s": time.time() - t_all,
            "euler_companion": (str(euler_out) if euler_out else None)}
    if ns_rec is not None:
        logd["ns_correction"] = ns_rec
    if mk_rec is not None:
        logd["M_knot_rule"] = mk_rec
    logd = _json_safe(logd)
    Path(str(out) + ".card_log.json").write_text(json.dumps(logd, indent=1, ensure_ascii=False))
    if mk_rec is not None and mk_rec.get("mode") == "auto":
        for y in [out] + ([euler_out] if euler_out else []):      # 書いた YAML の要約と card_log の記録が結べること (sha256・作り直した要約)
            chk = mknot_card_log_check(y)
            if not chk["ok"]:
                raise RuntimeError(f"M_knot: {y.name} の _card.M_knot_rule (要約) と card_log の M_knot_rule が一致しない ({chk})")
    return logd


def write_euler_companion(out: Path):
    """最終の雛形 (`physical_wall_repr: single_bspline`・`pw_upstream: poly`) の問題は Euler の prepare が受け付けない (物理壁のキーは
    NS 専用; runner_axismach.prepare)。同じ設計壁で Euler の参照を回すために、その 2 行だけを除いた `<stem>_euler.yaml` を書く
    (2026-10-13 case/70 で判明)。他のキーはそのまま (設計壁・MOC・mesh_euler が同じになる)。該当しなければ None。"""
    from ..probdef import load_problem
    out = Path(out)
    lines = out.read_text().splitlines(keepends=True)
    drop = ("physical_wall_repr:", "pw_upstream:")          # Euler の prepare が拒否する物理壁のキー (runner_axismach.prepare)
    keep = [ln for ln in lines if not ln.lstrip().startswith(drop)]
    if len(keep) == len(lines):
        return None
    eo = out.with_name(out.stem + "_euler" + out.suffix)
    hdr = (f"# Euler の参照用 ({out.name} から geometry.physical_wall_repr と geometry.pw_upstream の行だけを除いた写し; 設計壁・MOC・"
           "mesh_euler は同じ。物理壁のキーは NS 専用で、Euler の prepare は single_bspline と pw_upstream: poly を拒否する)\n")
    eo.write_text(hdr + "".join(keep))
    load_problem(eo)
    return eo


def _card_record(card: dict) -> dict:
    """`_card.card` に残すカードの中身 (既定値を埋めた形)。"""
    return {k: copy.deepcopy(card[k]) for k in ("format", "name", "gas", "conditions", "M_exit", "throat", "inlet",
                                                "divergent_length_m", "wall", "p_ambient", "tolerances", "recipe", "overrides")}


def _json_safe(v):
    if isinstance(v, dict):
        return {str(k): _json_safe(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_json_safe(x) for x in v]
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.ndarray):
        return [_json_safe(x) for x in v.tolist()]
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    return repr(v)


def _header_wall_temperature(card, est) -> list:
    """壁温分布のカードのヘッダの行 (一様・断熱は無し)。"""
    tw = est.get("wall_temperature")
    if tw is None:
        return []
    w = card["wall"]
    h = [f"壁温: {est['inputs']['wall']} → spec.wall_thermal.Tw_profile (表をそのまま写した; {_COORD_TEXT[tw['coord']]})。"
         "NS (wallProfile の CSV) と積分法 (Tw_table) は prepare_ns がその場の寸法で実寸の x に写して読む"]
    if "Tw_profile_csv" in w:
        h.append(f"  壁温の CSV: {w['Tw_profile_csv']['abs_path']} (sha256 {w['Tw_profile_csv']['sha256']}; カードの値 "
                 f"{w['Tw_profile_csv']['path']!r})")
    anc = tw["anchors_m"]
    where = ("区切りの点は設計壁から: 入口端 {x_in_m:.9g}・縮流部の始まり {x_contraction_m:.9g}・スロート {x_throat_m:g}・壁の終点 {x_exit_m:.9g} m"
             .format(**anc) + f" (カードの拡大部長との差 {tw['wall_end_minus_card_exit_m'] * 1e3:+.3f} mm)" if anc is not None else "x_m の表")
    h.append(f"  格子の見積もりは断面ごとの局所の T_w(x) ({where}; 第 1 層を決めた断面の T_w: 下流 {tw['Tw_K_at_wall_first_frac']:.1f} K、"
             f"スロート付近 {tw['Tw_K_at_wall_first_frac_throat']:.1f} K)")
    return h


def _nodes_text(nodes: dict) -> str:
    """NS・Euler の実効の節点数と上限 (YAML のコメントと標準出力に同じ文言で出す; plan tooling-nozzle-mesh-size-guard §4.2)。"""
    def one(label, r):
        src = "既定" if r["max_nodes_source"] == "default" else f"明示: {r['max_nodes_reason']}"
        return f"{label} ni {r['ni']} × nj {r['nj']} = {r['n_nodes']} (上限 {r['max_nodes']}、{src})"
    return ("節点数 (格子を作る前の判定): " + one("NS (mesh)", nodes["mesh"]) + "・" + one("Euler (mesh_euler)", nodes["mesh_euler"])
            + "。上限を超える格子はユーザの了承を取ってから overrides の max_nodes と max_nodes_reason を組で書いて許す")


def _header(card, prob, card_path, card_sha, method_text, est, liq, dims, warnings, auto_adjust, nodes, mk_full=None,
            card_log_name=None) -> list:
    """YAML の見出しのコメント。mk_full: card_log の M_knot_rule (記録の全部; YAML の _card.M_knot_rule は要約なので、選び直しの文・
    所要などはこちらから書く)。card_log_name: card_log のファイル名 (M_knot の記録の置き場の行)。"""
    rt = float(prob["spec"]["r_throat"])
    rec = prob["_card"]["recipe"]
    h = [f"生成: {GENERATOR} ({prob['_card']['generated']})。**手で編集しない** — 要求カードを直して作り直す。",
         f"  カード: {card_path} (sha256 {card_sha[:16]}…)、雛形: {rec['name']} ({rec['path']}, sha256 {rec['sha256'][:16]}…)",
         f"r_t = {rt!r} m: {method_text}",
         "長さの写し: geometry.r_inlet・L_U・L_pipe・dv.L_total はカードの実寸 / r_t。L_total_tol = 長さ公差 / 2 / r_t "
         "(r_t を変えたら生成器で作り直す — 手で r_throat だけ変えると実寸がずれる)",
         f"p_ambient: {prob['_card']['p_ambient_source']}",
         f"格子: mesh.ni {est['ni']}・wall_first_frac {est['wall_first_frac']:g}・wall_first_frac_throat {est['wall_first_frac_throat']} は"
         f" estimate_mesh の**初期値** (y1+ {Y1PLUS_TARGET}、Eckert + 平板乱流 Cf)。合格の根拠ではない — NS 後に check_wall_resolution.py"
         " の VERDICT で判定し、不合格なら細分する (手順書 §2 ⑥)。実際の格子の AR の最大 "
         f"{est['mesh_check']['ar_max']:.0f} (ar_max {est['ar_max_limit']:g})",
         _nodes_text(nodes),
         *_header_wall_temperature(card, est),
         f"液化 (N₂、出口 {liq['T_exit_K']:.2f} K・{liq['p_exit_Pa']:.4g} Pa): {liq['verdict']}",
         "設計壁から測った寸法 (§4.2a の前段; 最終の合格は成果物から測る出荷ゲート): "
         + "、".join(f"{r['name']} {r['measured_m']:.6f} m (差 {r['diff_m']:+.2e}、公差 {r['tol_m']:g}、{'ok' if r['pass'] else 'FAIL'})"
                    for r in dims["rows"]),
         ]
    if dims.get("x_F_minus_x_e_m") is not None:
        h.append(f"  MOC の出口 x_F と設計壁の終点 x_e の差 {dims['x_F_minus_x_e_m'] * 1e3:+.3f} mm (from_length は x_F を合わせ、"
                 "壁とメッシュは x_e で終わる)")
    ns = prob["_card"].get("ns_correction")
    if ns:
        e = ns["table_ends"]
        h.append(f"NS 後の補正: 前の run {ns['prev_run']} の {ns['delta_next']} (sha256 {ns['delta_next_sha256'][:16]}…) の δ_E × "
                 f"(r_t/S_prev)^−0.2 (S_prev {ns['S_prev_m']!r} m、factor {ns['factor']!r})。x は {ns['x_map']}。表の端を前の設計に揃えた差 "
                 f"始点 {e['start']['diff_rt']:+.3e}・終点 {e['end']['diff_rt']:+.3e} r_t (ε_f32 {e['start']['eps_f32_rt']:.2e}・"
                 f"{e['end']['eps_f32_rt']:.2e})")
        h.append(f"  補正の表: {ns['table']} (sha256 {ns['table_sha256'][:16]}…) — NS は deltastar_loop --delta-r-csv で立てる "
                 "(問題 YAML の _card.ns_correction.table_sha256 と照合する)")
    fd = prob["_card"].get("moc_field")
    if fd:
        h.append(f"逆 MOC の場: x_F {fd['x_F']:.6f} r_t、切る前の壁の終点 {fd['x_wall_full_end']:.6f} r_t (x_F の {fd['margin_to_field_end_rt']:.4g} r_t "
                 f"下流、x_F より下流の壁点 {fd['n_beyond_x_F']} 個、最後の点間隔 {fd['spacing_end']:.4g} r_t) — "
                 f"{'**場の終点に張り付き**' if fd['pinned'] else '張り付きなし'} (x_end_margin {fd['x_end_margin']})")
    rl = prob["_card"].get("L_total_tol_relaxed")
    if rl:
        h.append(f"**from_length の内部許容差を緩めた**: L_total_tol {rl['from']:.4g} → {rl['to']} r_t ({_one_line(rl['reason'])[:200]})。"
                 "カードの公差は変えていない — 寸法の合否は上の測定と成果物の出荷ゲートがカードの公差で決める")
    if est.get("blend_choice") not in (None, "recipe"):
        h.append(f"格子の縮流部側のブレンド: {est['blend_choice']} ({est['blend_text']}) に決め直した — y1+ の目安の最大 "
                 f"{est['y1plus_est_max_with_blend']:.2f} @ x = {est['x_y1plus_est_max_rt']:.2f} r_t "
                 f"(雛形の値では {est['blend_candidates'][0]['y1plus_est_max']:.2f} @ x = {est['blend_candidates'][0]['x_y1plus_est_max_rt']:.2f})")
    if mk_full is not None:
        h += _mknot_summary_lines(mk_full, card_log_name)
    for a in auto_adjust:
        h.append(f"自動調整: {a['key']} {a['from']} → {a['to']} ({a['reason'][:120]})")
    for o in prob["_card"]["overrides"]:
        before = "(無し)" if o["before_missing"] else _flow(o["before"])
        h.append(f"上書き (overrides): {o['key']} = {_flow(o['value'])} (雛形/生成値 {before})")
    h.append("雛形に入れていない較正値: initial_line: cfd・Md_moc_offset・deltastar_initializer.cf_scale は手順書 §2 ⑤ の手順で後から足す")
    for w in warnings:
        h.append(f"警告: {w}")
    return h


# ---------------------------------------------------------------------------------------------------- CLI

def _print_report(logd: dict) -> None:
    print("カードの値 → 問題 YAML のキー:")
    for r in logd["mapping"]:
        cv, yv = r["card_value"], r["yaml_value"]
        cv = _flow(cv) if isinstance(cv, (dict, list)) else cv
        yv = _flow(yv) if isinstance(yv, (dict, list)) else yv
        print(f"  {r['card_key']:<28} {str(cv):<26} → {r['yaml_key']:<22} {yv}   ({r['rule']})")
    rt = logd["r_t"]
    print(f"r_t = {rt['r_t_m']!r} m ({rt['method']})")
    for h in rt["history"]:
        print(f"  反復 {h['iter']}: {h['r_t_in_m']!r} → {h['r_t_out_m']!r} m  |Δ|/r_t {h['rel_change']:.3e}  "
              f"物理スロート半径 {h['physical_throat_r_m']!r} m  ({h['time_s']:.1f} s)")
    e = logd["mesh_estimate"]
    print(f"格子 (初期値): ni {e['ni']}  wall_first_frac {e['wall_first_frac']}  wall_first_frac_throat {e['wall_first_frac_throat']}"
          f"  ブレンド {e['blend_choice']} (up_x0/x1 {e['wall_first_up_x0']}/{e['wall_first_up_x1']})  (AR の最大 {e['mesh_check']['ar_max']:.0f}、"
          f"y1+ の目安の最大 {e['y1plus_est_max_with_blend']:.2f} @ x {e['x_y1plus_est_max_rt']:.2f})")
    print(_nodes_text(logd["nodes"]))
    tw = e.get("wall_temperature")
    if tw is not None:
        print(f"壁温 (格子の見積もり): {e['inputs']['wall']}、区切りの点 [m] {tw['anchors_m']}、第 1 層を決めた断面の T_w "
              f"下流 {tw['Tw_K_at_wall_first_frac']:.1f} K・スロート付近 {tw['Tw_K_at_wall_first_frac_throat']:.1f} K")
    f = logd["moc_field"]
    print(f"逆 MOC の場: x_F {f['x_F']:.6f}  x_e {f['x_e']:.6f}  x_F − x_e {f['x_F_minus_x_e_m'] * 1e3:+.3f} mm  場の終点 "
          f"{f['x_wall_full_end']:.4f} r_t ({'張り付き' if f['pinned'] else '張り付きなし'}、x_end_margin {f['x_end_margin']})")
    if logd.get("M_knot_rule") is not None:
        for ln in _mknot_summary_lines(logd["M_knot_rule"], Path(logd["out"]).name + ".card_log.json"):
            print(ln)
    for a in logd["auto_adjust"]:
        print(f"自動調整: {a['key']} {a['from']!r} → {a['to']!r}  ({_one_line(a['reason'])[:160]})")
    if logd.get("L_total_tol_relaxed"):
        rl = logd["L_total_tol_relaxed"]
        print(f"from_length の内部許容差を緩めた: {rl['from']:.4g} → {rl['to']} r_t (カードの公差は不変)")
    print(f"液化: {logd['liquefaction']['verdict']}")
    ns = logd.get("ns_correction")
    if ns:
        e = ns["table_ends"]
        print(f"NS 後の補正: {ns['delta_next']} (前の run {ns['prev_run']}、S_prev {ns['S_prev_m']!r} m)、x の写し {ns['x_map']}、"
              f"表の端の差 始点 {e['start']['diff_rt']:+.3e}・終点 {e['end']['diff_rt']:+.3e} r_t → 表 {ns['table']} "
              f"(sha256 {ns['table_sha256'][:16]}…、x {ns['table_x_range_rt'][0]:.9g}〜{ns['table_x_range_rt'][1]:.9g} r_t、"
              f"最後の内側の表と同一: {ns['final_table_equals_last_inner']})")
    d = logd["dimensions"]
    print(f"寸法 (設計壁から): {'PASS' if d['pass'] else 'FAIL'}")
    for r in d["rows"]:
        print(f"  {r['name']:<26} カード {r['card_m']:.6f} m  測定 {r['measured_m']:.6f} m  差 {r['diff_m']:+.3e} m  "
              f"公差 {r['tol_m']:g}  {'ok' if r['pass'] else 'FAIL'}")
    for w in logd["warnings"]:
        print(f"警告: {w}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="風洞ノズルの要求カード (forge.wt_card/v1) → 問題 YAML (wind_tunnel_axisym_axismach)")
    ap.add_argument("card", help="要求カードの YAML")
    ap.add_argument("--out", required=True, help="書き出す問題 YAML (同じ場所に <out>.card_log.json も書く)")
    ap.add_argument("--r-throat", type=float, default=None, metavar="R",
                    help="外から r_t [m] を与えて解かずに書く (NS 後の解き直しの割り直し)")
    ap.add_argument("--no-solve", action="store_true", help="物理スロート径でも r_t を解かない (仮に 径/2)")
    ap.add_argument("--max-outer", type=int, default=10)
    ap.add_argument("--rel-tol", type=float, default=1e-6)
    ap.add_argument("--ns-correction", default=None, metavar="PREV_RUN",
                    help="NS 後の補正: r_t の固定点の内側を PREV_RUN の δ_E の全分布で解き (solve_rt_throat --prev と同じ)、最終の r_t の補正の表を "
                         "YAML の横の <stem>_delta_r.csv に書く (出典と sha256 は _card.ns_correction)")
    ap.add_argument("--delta-next", default=None, metavar="CSV",
                    help="--ns-correction のとき δ_E を読む delta_r_next.csv (既定 PREV_RUN/delta_r_next.csv)")
    a = ap.parse_args(argv)
    if a.delta_next is not None and a.ns_correction is None:
        ap.error("--delta-next は --ns-correction PREV_RUN と一緒に使う")
    logd = build_problem(a.card, a.out, sizing=("none" if a.no_solve else "auto"), max_outer=a.max_outer, rel_tol=a.rel_tol,
                         r_throat_m=a.r_throat, ns_prev_run=a.ns_correction, ns_delta_next=a.delta_next)
    _print_report(logd)
    print(f"→ {a.out}  ({a.out}.card_log.json)")
    if logd.get("ns_correction"):
        print(f"→ {logd['ns_correction']['table']}  (NS 後の補正の表。deltastar_loop --delta-r-csv に渡す)")
    eo = Path(a.out).with_name(Path(a.out).stem + "_euler" + Path(a.out).suffix)
    if eo.exists():
        print(f"→ {eo}  (Euler の参照用。physical_wall_repr・pw_upstream を除いた写し)")
    return 0 if logd["dimensions"]["pass"] else 2


if __name__ == "__main__":
    sys.exit(main())
