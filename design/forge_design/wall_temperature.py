"""壁温分布 $T_w(x)$ の単一ソース (plan tooling-nozzle-wall-temperature-profile §4.1・§4.2)。

問題 YAML の `spec.wall_thermal.Tw_profile` を検査し、その場の寸法 (scale = 実効の r_t [m]) で実寸の x [m] に写して、
積分法の初期壁 (`thermal_bc.Tw_table`、x は r_t 単位)・ソルバの CSV (`wall_profile_<physID>.csv`、x は m)・成果物の 3 箇所へ渡す。

- 座標は 2 通り (どちらか 1 つ):
  - `s` (区間で正規化): 入口直管 [入口端 → 縮流部の始まり] を [−2, −1]、縮流部 [縮流部の始まり → 設計スロート] を [−1, 0]、
    拡大部 [設計スロート → 壁の終点] を [0, 1] に区間ごとに線形で写す。区切りの点はその設計で実際にできた壁から取る
    (出口は MOC の x_F でなく壁の終点 `d["wall"].x_e`)。表の外側の点 (s < −2・s > 1) は同じ区間の写像を延長する。
  - `x_m` (実寸 [m]、原点は設計スロート): 長さを変えても表は動かない。
- 補間は**実寸の x の 1D 線形** (ソルバと同じ)。`s` の表の点は先に x に写してから x で補間する
  (区切りの点をまたぐ区間では s で線形補間した値とは違う)。端値クランプには頼らず、使う範囲を覆わない表は例外 (§4.1、codex plan M2)。
- 写しは**呼ぶたびに**その場の寸法で作る (表を `spec.r_throat` で換算してキャッシュしない; codex 諮問 2026-10-13)。
- run の中の `wall_temperature_profile.json` と `wall_profile_<physID>.csv` が実行条件の正本 (codex plan M3)。
  `prepare_info.json` の `wall_thermal.files` に両方の sha256 を記録し、起動前・成果物の書き出し前に `load_run_profile` で照合する。

一様の $T_w$ と断熱の問題では、ここを通っても生成物 (bcond・`thermal_bc`) は変更前とバイト同一 (`integral_thermal_bc` は
`Problem.wall_thermal_bc_integral` をそのまま返す)。
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

RUN_JSON = "wall_temperature_profile.json"
FORMAT = "forge.wall_temperature_profile/v1"
COVER_TOL_M = 1e-12          # 範囲の検査の許容差 [m] (区切りの点の写しの丸め分)
S_LO, S_HI = -2.0, 1.0       # s の表が覆うべき範囲 (入口端 〜 壁の終点)


def csv_name(wall_id) -> str:
    """ソルバが読む壁温の CSV の名前 (`wall_profile_<physID>.csv`; `boundaryCond.cpp` の applyWallProfiles)。"""
    if isinstance(wall_id, bool) or int(wall_id) != wall_id:
        raise ValueError(f"wall_id は整数の physID (受け取った値: {wall_id!r})")
    return f"wall_profile_{int(wall_id)}.csv"


def _num_list(v, what: str) -> list:
    """1 次元の有限な数値の列を float の list にする (文字列・真偽値・入れ子・NaN・Inf は理由つきで拒否)。"""
    if isinstance(v, np.ndarray):
        if v.ndim != 1:
            raise ValueError(f"{what} は 1 次元の数値の列 (受け取った配列の次元 {v.ndim})")
        v = v.tolist()
    if not isinstance(v, (list, tuple)):
        raise ValueError(f"{what} は 1 次元の数値の列 (受け取った値: {v!r})")
    out = []
    for i, e in enumerate(v):
        if isinstance(e, (bool, np.bool_)) or not isinstance(e, (int, float, np.integer, np.floating)):
            raise ValueError(f"{what}[{i}] は数値であること (受け取った値: {e!r}; 1 次元の数値の列)")
        f = float(e)
        if not math.isfinite(f):
            raise ValueError(f"{what}[{i}] = {e!r} は有限でない (NaN・Inf は不可)")
        out.append(f)
    return out


def parse_profile(obj) -> dict:
    """`spec.wall_thermal.Tw_profile` を検査して `{"coord": "s" | "x_m", "coords": [...], "Tw_K": [...]}` を返す。

    mapping `{s: [...], Tw_K: [...]}` か `{x_m: [...], Tw_K: [...]}`。座標キーは厳密に 1 つ、未知のキーは不可。
    座標と `Tw_K` は同じ長さの 1 次元の有限な数値列・2 点以上、座標は狭義単調増加、`Tw_K` > 0。
    `s` の表は [−2, 1] (入口端〜壁の終点) を覆うこと (許容 1e-12)。`x_m` の表は設計スロート x = 0 を含むこと
    (代表値 = x = 0 の値を表の中で取るため。ノズル全長を覆う表なら必ず含む)。違反は理由つきの ValueError。"""
    if not isinstance(obj, dict):
        raise ValueError("Tw_profile は mapping ({s: [...], Tw_K: [...]} か {x_m: [...], Tw_K: [...]}) "
                         f"(受け取った値: {obj!r})")
    ck = [k for k in ("s", "x_m") if k in obj]
    if len(ck) != 1:
        raise ValueError(f"Tw_profile の座標キーは s か x_m のどちらか 1 つ (受け取ったキー: {sorted(map(str, obj))})")
    unknown = sorted(str(k) for k in obj if k not in ("s", "x_m", "Tw_K"))
    if unknown:
        raise ValueError(f"Tw_profile の未知のキー {unknown} (使えるのは s か x_m の一方と Tw_K)")
    if "Tw_K" not in obj:
        raise ValueError("Tw_profile に Tw_K (壁温 [K] の列) が無い")
    c = ck[0]
    xs = _num_list(obj[c], f"Tw_profile.{c}")
    Ts = _num_list(obj["Tw_K"], "Tw_profile.Tw_K")
    if len(xs) != len(Ts):
        raise ValueError(f"Tw_profile.{c} ({len(xs)} 点) と Tw_K ({len(Ts)} 点) の長さが違う")
    if len(xs) < 2:
        raise ValueError(f"Tw_profile は 2 点以上 (受け取った点数 {len(xs)})")
    bad = [i for i in range(len(xs) - 1) if not xs[i + 1] > xs[i]]
    if bad:
        i = bad[0]
        raise ValueError(f"Tw_profile.{c} は狭義単調増加であること ({c}[{i}] = {xs[i]!r} → {c}[{i + 1}] = {xs[i + 1]!r})")
    badT = [i for i, t in enumerate(Ts) if not t > 0.0]
    if badT:
        i = badT[0]
        raise ValueError(f"Tw_profile.Tw_K は正の温度 [K] (Tw_K[{i}] = {Ts[i]!r})")
    if c == "s":
        if not (xs[0] <= S_LO + COVER_TOL_M and xs[-1] >= S_HI - COVER_TOL_M):
            raise ValueError(f"Tw_profile.s の表は [−2, 1] (入口端〜壁の終点) を覆うこと (表の範囲 [{xs[0]!r}, {xs[-1]!r}])")
    elif not (xs[0] <= 0.0 <= xs[-1]):
        raise ValueError(f"Tw_profile.x_m の表は設計スロート x = 0 を含むこと (ノズル全長を覆う表; 表の範囲 [{xs[0]!r}, {xs[-1]!r}] m)")
    return {"coord": c, "coords": xs, "Tw_K": Ts}


def representative_Tw(prof: dict) -> float:
    """代表値: 表の座標 0 (s = 0 か x = 0 = 設計スロート) の値 (bcond の `Ts` の既定値・帳簿用; §4.2)。"""
    c = np.asarray(prof["coords"], dtype=float)
    if not (c[0] <= 0.0 <= c[-1]):
        raise ValueError(f"壁温の表が座標 0 (設計スロート) を含まない (表の範囲 [{c[0]!r}, {c[-1]!r}])")
    return float(np.interp(0.0, c, np.asarray(prof["Tw_K"], dtype=float)))


def anchors_m(L_pipe_m: float, L_U_m: float, x_exit_m: float) -> dict:
    """区切りの点 [m] (原点 = 設計スロート): 入口端・縮流部の始まり・設計スロート・壁の終点。"""
    vals = {"L_pipe_m": L_pipe_m, "L_U_m": L_U_m, "x_exit_m": x_exit_m}
    for k, v in vals.items():
        if isinstance(v, bool) or not math.isfinite(float(v)):
            raise ValueError(f"anchors_m: {k} は有限の数値 (受け取った値: {v!r})")
    L_pipe_m, L_U_m, x_exit_m = float(L_pipe_m), float(L_U_m), float(x_exit_m)
    if not (L_pipe_m >= 0.0 and L_U_m > 0.0 and x_exit_m > 0.0):
        raise ValueError(f"anchors_m: L_pipe ≥ 0・L_U > 0・x_exit > 0 であること (L_pipe {L_pipe_m!r}、L_U {L_U_m!r}、x_exit {x_exit_m!r} m)")
    return {"x_in_m": -(L_pipe_m + L_U_m), "x_contraction_m": -L_U_m, "x_throat_m": 0.0, "x_exit_m": x_exit_m}


def anchors_from_design(p, d: dict, scale: float) -> dict:
    """問題 `p` の `geometry.L_pipe`・`L_U` (r_t 単位; 既定は `design_chain` と同じ 0.5・3.5) と設計の結果 `d` の壁の終点
    `d["wall"].x_e` (r_t 単位) に `scale` [m] を掛けて `anchors_m` を呼ぶ。区切りの点はその設計で実際にできた壁から取る
    (出口は MOC の x_F でなく壁の終点)。設計壁の入口端 `x_in` が −(L_U + L_pipe) と合わなければ例外 (別の壁を混ぜない)。"""
    g = p.geometry
    L_pipe, L_U = float(g.get("L_pipe", 0.5)), float(g.get("L_U", 3.5))
    w = d["wall"]
    x_in_w = getattr(w, "x_in", None)
    if x_in_w is not None and abs(float(x_in_w) - (-(L_pipe + L_U))) > 1e-12 * max(1.0, abs(L_pipe + L_U)):
        raise ValueError(f"anchors_from_design: 設計壁の入口端 x_in {float(x_in_w)!r} r_t が −(L_U + L_pipe) = "
                         f"{-(L_pipe + L_U)!r} r_t と合わない (問題と設計の結果が別物)")
    s = float(scale)
    if not (math.isfinite(s) and s > 0.0):
        raise ValueError(f"anchors_from_design: scale は正の有限の長さ [m] (受け取った値: {scale!r})")
    return anchors_m(L_pipe * s, L_U * s, float(w.x_e) * s)


def to_x_m(prof: dict, anchors: dict):
    """表を実寸の x [m] に写す: 戻り (x_m, Tw) の ndarray。`x_m` の表はそのまま、`s` の表は区間ごとの線形写像
    (s ∈ [−2,−1] → [x_in, x_contraction]、[−1,0] → [x_contraction, 0]、[0,1] → [0, x_exit]; 外側は同じ区間の写像を延長)。
    区切りの点 (−2・−1・0・1) は区切りの x にそのまま写る。写した x が狭義単調増加でなければ例外。"""
    c = np.asarray(prof["coords"], dtype=float)
    T = np.asarray(prof["Tw_K"], dtype=float)
    kind = prof["coord"]
    if kind == "x_m":
        x = c.copy()
    elif kind == "s":
        xi, xc = float(anchors["x_in_m"]), float(anchors["x_contraction_m"])
        xt, xe = float(anchors["x_throat_m"]), float(anchors["x_exit_m"])
        if not (xi < xc < xt < xe):
            raise ValueError(f"to_x_m: 区切りの点が狭義単調増加でない (入口端 {xi!r}・縮流部の始まり {xc!r}・スロート {xt!r}・"
                             f"壁の終点 {xe!r} m) — s の表を写せない (入口直管の長さ 0 など)")
        x = np.where(c < -1.0, xi + (c + 2.0) * (xc - xi),
                     np.where(c < 0.0, xc + (c + 1.0) * (xt - xc), xt + c * (xe - xt)))
    else:
        raise ValueError(f"to_x_m: 座標の種類 {kind!r} は未知 (s | x_m)")
    if not np.all(np.diff(x) > 0.0):
        raise ValueError("to_x_m: 写した x [m] が狭義単調増加でない")
    return x, T


def check_covers(x_m, lo_m: float, hi_m: float, what: str) -> None:
    """表の x [m] が [lo_m, hi_m] を覆うか (許容 1e-12 m)。覆わなければ必要な範囲と表の範囲を示して ValueError
    (`np.interp` の端値クランプに頼らない; plan §4.1、codex plan レビュー M2)。"""
    x = np.asarray(x_m, dtype=float)
    if not (x[0] <= float(lo_m) + COVER_TOL_M and x[-1] >= float(hi_m) - COVER_TOL_M):
        raise ValueError(f"壁温の表が{what}を覆わない: 必要な範囲 [{float(lo_m)!r}, {float(hi_m)!r}] m・"
                         f"表の範囲 [{float(x[0])!r}, {float(x[-1])!r}] m (端値で延長しない; 表を広げる)")


def eval_Tw(x_m_tab, Tw_tab, x):
    """実寸の表 (x [m], Tw [K]) を x [m] で線形補間する (ソルバと同じ)。範囲外 (許容 1e-12 m を超える) は例外 —
    呼ぶ側が先に `check_covers` で使う範囲を検査する。"""
    xt = np.asarray(x_m_tab, dtype=float)
    Tt = np.asarray(Tw_tab, dtype=float)
    xq = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(xq)):
        raise ValueError("eval_Tw: 評価点に非有限値がある")
    if np.any(xq < xt[0] - COVER_TOL_M) or np.any(xq > xt[-1] + COVER_TOL_M):
        raise ValueError(f"eval_Tw: 評価点 [{float(np.min(xq))!r}, {float(np.max(xq))!r}] m が表の範囲 "
                         f"[{float(xt[0])!r}, {float(xt[-1])!r}] m の外 (端値で延長しない)")
    return np.interp(xq, xt, Tt)


def _profile_of(p):
    wt = p.wall_thermal
    return wt.get("Tw_profile") if wt["mode"] == "isothermal" else None


def integral_thermal_bc(p, d: dict, scale: float) -> dict:
    """積分法 (`deltastar_integral.integral_bl`) に渡す thermal_bc。**scale は積分法を呼ぶその場の実効の r_t [m]**
    (寸法の逆算の反復中は試しの r_t; plan §4.2 の契約)。`integral_bl` に渡すのと同じ値を渡すこと。

    分布: 区切りの点 → 実寸の x [m] → その scale で積分区間 (入口端〜壁の終点; 積分法の終点 `wall_inv` の最後の x も含める) を
    覆うか検査 → `{"mode": "prescribed_temperature", "Tw_table": [[x/scale, T], ...]}`。
    一様・断熱は `Problem.wall_thermal_bc_integral` そのもの (変更前とバイト同一の dict)。"""
    prof = _profile_of(p)
    if prof is None:
        return p.wall_thermal_bc_integral
    s = float(scale)
    anc = anchors_from_design(p, d, s)
    x_m, Tw = to_x_m(prof, anc)
    hi = anc["x_exit_m"]
    if d.get("wall_inv") is not None:
        hi = max(hi, float(np.asarray(d["wall_inv"])[-1, 0]) * s)     # integral_bl の積分の終点 (= x_F の壁表の最後)
    check_covers(x_m, anc["x_in_m"], hi, f"積分法の積分区間 (入口端〜壁の終点、scale = {s!r} m)")
    return {"mode": "prescribed_temperature", "Tw_table": [[float(x / s), float(t)] for x, t in zip(x_m, Tw)]}


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_run_files(run_dir, p, d: dict, scale: float, wall_id: int = 3) -> dict:
    """run に壁温分布の正本 2 つを書く: `wall_temperature_profile.json` (元の座標の種類・元の表・区切りの点・実寸の表) と
    `wall_profile_<wall_id>.csv` (ソルバの `wallProfile` が読む; ヘッダ `x Ts`・空白区切り、x [m]、値は `repr` の往復できる桁)。
    設計壁 (入口端〜壁の終点) を覆わない表は例外。戻り: `{"json": 名, "json_sha256", "csv": 名, "csv_sha256"}`
    (`prepare_info.json` の `wall_thermal.files` にそのまま記録する)。

    CSV の区切りは**空白** (`boundaryCond.cpp` の applyBoundaryProfiles はヘッダもデータも空白で分割する。`x,Ts` と書くと
    座標列が無いとして起動時に落ちる; case/48 `run_0015_wallprofile` の `wall_profile_4.csv` と同じ書式)。"""
    prof = _profile_of(p)
    if prof is None:
        raise ValueError("write_run_files: 壁温分布 (spec.wall_thermal.Tw_profile) の問題でない")
    rd = Path(run_dir)
    s = float(scale)
    anc = anchors_from_design(p, d, s)
    x_m, Tw = to_x_m(prof, anc)
    check_covers(x_m, anc["x_in_m"], anc["x_exit_m"], f"設計壁 (入口端〜壁の終点、scale = {s!r} m)")
    cname = csv_name(wall_id)
    csv_text = "x Ts\n" + "".join(f"{float(x)!r} {float(t)!r}\n" for x, t in zip(x_m, Tw))
    (rd / cname).write_bytes(csv_text.encode("ascii"))
    rec = {"format": FORMAT, "wall_id": int(wall_id), "scale_m": s,
           "coord": prof["coord"], "coords": [float(v) for v in prof["coords"]], "Tw_K": [float(v) for v in prof["Tw_K"]],
           "Tw_representative_K": representative_Tw(prof),
           "anchors_m": {k: float(v) for k, v in anc.items()},
           "x_m": [float(v) for v in x_m], "Tw_x_K": [float(v) for v in Tw],
           "csv": cname,
           "note": "実行条件の正本 (plan tooling-nozzle-wall-temperature-profile §4.2)。x_m の原点は設計スロート。"
                   "ソルバは csv を壁ノードの x で 1D 線形補間する (範囲外は端値)"}
    (rd / RUN_JSON).write_bytes(json.dumps(rec, indent=1, ensure_ascii=False).encode("utf-8"))
    return {"json": RUN_JSON, "json_sha256": _sha256(rd / RUN_JSON), "csv": cname, "csv_sha256": _sha256(rd / cname)}


def _bcond_wall_profile_flag(bc_path: Path, wall_id: int):
    """bcondConfig.yaml の physID = wall_id の行の (kind, ints.wallProfile)。読めない・無いときは None。"""
    import yaml
    try:
        doc = yaml.safe_load(bc_path.read_text())
    except Exception:          # noqa: BLE001 (読めないことを呼ぶ側で理由にする)
        return None
    if not isinstance(doc, dict):
        return None
    for v in doc.values():
        if isinstance(v, dict) and v.get("physID") == wall_id:
            ints = v.get("ints") if isinstance(v.get("ints"), dict) else {}
            return str(v.get("kind")), ints.get("wallProfile")
    return None


def load_run_profile(run_dir):
    """run の壁温分布の正本を読み、`prepare_info.json` の `wall_thermal.files` に記録した sha256 と JSON・CSV の実ファイルの
    sha256 を照合する。bcondConfig.yaml の壁の行に `wallProfile: 1` があることも確かめる。不一致・欠けは ValueError
    (CSV が書き換えられた / prepare_info と合わない)。戻り: JSON の中身 (`x_m`・`Tw_x_K` を含む)。

    一様・断熱の run (記録に files が無く、JSON も bcond の `wallProfile` も無い) は None。記録が無いのに JSON か
    `wallProfile` があるのは「prepare_info と合わない」として例外 (記録の無い分布でソルバを回さない)。"""
    rd = Path(run_dir)
    pi = rd / "prepare_info.json"
    try:
        info = json.loads(pi.read_text()) if pi.exists() else {}
    except (OSError, ValueError):
        info = None          # 読めない記録: 下で分布の痕跡 (JSON・wallProfile) があれば止め、無ければ一様・断熱として扱う
    wt = info.get("wall_thermal") if isinstance(info, dict) else None
    files = wt.get("files") if isinstance(wt, dict) else None
    bc_path = rd / "bcondConfig.yaml"
    if not files:
        stray = []
        if (rd / RUN_JSON).exists():
            stray.append(RUN_JSON)
        if bc_path.exists() and "wallProfile" in bc_path.read_text():
            stray.append("bcondConfig.yaml の wallProfile")
        if isinstance(wt, dict) and "Tw_profile" in wt:
            stray.append("prepare_info.json の wall_thermal.Tw_profile (files の記録が無い)")
        if stray:
            raise ValueError(f"{rd}: 壁温分布 ({', '.join(stray)}) があるのに prepare_info.json の wall_thermal.files に"
                             " sha256 の記録が無い — prepare_info と合わない")
        return None
    miss = [k for k in ("json", "json_sha256", "csv", "csv_sha256") if k not in files]
    if miss:
        raise ValueError(f"{rd}: prepare_info.json の wall_thermal.files に {miss} が無い — prepare_info と合わない")
    jp, cp = rd / str(files["json"]), rd / str(files["csv"])
    for f in (jp, cp):
        if not f.is_file():
            raise ValueError(f"{rd}: 壁温分布の正本 {f.name} が無い (欠け) — prepare_info と合わない")
    got_j, got_c = _sha256(jp), _sha256(cp)
    if got_j != files["json_sha256"]:
        raise ValueError(f"{jp}: 壁温の表 (JSON) が書き換えられた / prepare_info と合わない "
                         f"(sha256 {got_j[:16]}… ≠ 記録 {str(files['json_sha256'])[:16]}…)")
    if got_c != files["csv_sha256"]:
        raise ValueError(f"{cp}: 壁温の CSV が書き換えられた / prepare_info と合わない "
                         f"(sha256 {got_c[:16]}… ≠ 記録 {str(files['csv_sha256'])[:16]}…)")
    rec = json.loads(jp.read_text())
    if rec.get("format") != FORMAT or rec.get("csv") != files["csv"]:
        raise ValueError(f"{jp}: format {rec.get('format')!r}・csv {rec.get('csv')!r} が記録 ({FORMAT!r}・{files['csv']!r}) と合わない")
    wid = rec.get("wall_id")
    if wid is not None:
        flag = _bcond_wall_profile_flag(bc_path, int(wid)) if bc_path.exists() else None
        if flag is None or flag[0] != "wall_isothermal" or flag[1] != 1:
            raise ValueError(f"{rd}: bcondConfig.yaml の physID {wid} が wall_isothermal + ints: {{wallProfile: 1}} でない "
                             f"({flag}) — 壁温分布の記録と合わない (ソルバが CSV を読まない)")
    return rec
