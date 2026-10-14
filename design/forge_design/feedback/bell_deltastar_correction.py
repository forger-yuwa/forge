r"""③ベル (M) の δ* 補正 — 設計者の選択肢 (plan tooling-nozzle-bell-two-mode §4.4「δ* 補正の提供の形」・§5.1 #7c)。

L100 で検証した一発補正 (§5.1 #7 (7b-4)、`case/68.bell_two_mode/rao_m/run_corrected_wall.py`) を関数にしたもの。縁の M だけ生産の取り方
(MOC の壁 + 凍結源 Euler、codex 2026-10-13 option-design 採用、kδ 据え置きの移植は (7c-0) で A) に替えている。

手順 (`correct_rao`):

1. 補正しない (M) の問題 (ε の非粘性壁 = `geometry.rao_wall_file`) から壁を組む (`runner.build_wall`)。
2. 理論式 `contur_v2` で $\delta_{r,CONTUR}(x)$ (縁 = `metrics.bell_deltastar.moc_edge_table`、断熱壁・NS と同じ μ(T))。
3. $\delta(x)=k_\delta\,\delta_{r,CONTUR}(x)$、出口で $\delta_e$ → $\varepsilon_{inv}=(\sqrt\varepsilon-\delta_e)^2$。
4. 同じ CFD ピンの初期線 (凍結源の Euler、`moc_initial_line`) と同じ L・円弧 kernel で (M) を $\varepsilon_{inv}$ で作り直す (`design_rao`)。
5. 物理壁のずれ `geometry.wall.DisplacementOffset` (x < 1.17 は V-ext = x ∈ [0, 0.41] の 5 次の立ち上げ、L100 の暫定既定)。

出力: `<parts_dir>/<tag>.json`・`<tag>_wall.npz` (`load_rao_wall` で読める形式) と `<tag>_offset.npz`。キャンペーンの要約
(`<parts_dir の親>/<prefix>.json` の `setup_by_tag`) は `write_campaign` が呼び出し側で書く (並列の生成で同じファイルを書き合わないため)。
補正後の壁で理論式を計算し直す反復はしない。kδ は呼び出し側が渡す (較正した条件: 断熱壁・case/68 の排気・燃焼圧 — 条件を変えたら NS で確かめる)。
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

RAMP_VEXT = (0.0, 0.41)            # x < 1.17 の既定 V-ext (codex 2026-10-13 result: L100 の暫定)
X_CAL_START = 1.17                 # 較正の区間の始まり (記録)
N_THEORY = 1501                    # 理論式を出す x の点数 (入口 x_in から出口まで等間隔)


def _sha16(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def wall_record(res: dict, tag: str, L: float, t_s: float, note: str = "") -> dict:
    """`design_rao` の戻りから `load_rao_wall` が読む記録 (`<tag>.json`) を作る (case の `run_rao_t6._record` と同じ鍵)。"""
    rec = {"tag": tag, "L_frac": None, "L_rt": float(L), "n": res["n"], "variant": "deltastar_corrected_inviscid", "dr": res["dr"],
           "m": res["m"], "q": res["q"], "kernel": res["kernel"], "c": res["c"], "accept_rule": res.get("accept_rule"),
           "status": res["status"], "why": res.get("why"), "rejected": res.get("rejected", []), "lost": res.get("lost", []),
           "policy_holes": res.get("policy_holes", []), "time_s": t_s, "note": note}
    if res["status"] != "ok":
        return rec
    E, th, cl, gt = res["E"], res["thrust"], res["closure"], res["g_tilde"]
    rec.update({
        "theta_B_deg": float(np.rad2deg(res["theta_B"])), "s_B": res["s_B"], "eta_D": res["eta_D"], "det_D": res["det_D"],
        "de": res["de"], "branch": res["branch"], "M_D": float(res["D"]["M"]), "theta_D_deg": float(np.rad2deg(res["D"]["th"])),
        "x_D": float(res["D"]["x"]), "r_D": float(res["D"]["r"]), "theta_E_deg": float(np.rad2deg(res["theta_E"])),
        "theta_max_deg": float(np.rad2deg(res["theta_max"])), "x_E": float(E["x"]), "r_E": float(E["r"]), "M_E": float(E["M"]),
        "x_E_minus_L": float(E["x"]) - float(L), "x_B": res["x_B"], "r_B": res["r_B"],
        "F_cs": th["F_cs"], "F_wall": th["F_wall"], "thrust_rel_diff": th["rel_diff"], "mdot_rel_diff": th["mdot_rel_diff"],
        "T_line": th["T_line"], "T_kernel": th["T_kernel"], "T_BE": th["T_BE"], "T_closure": th["T_closure"], "T_DE": th["T_DE"],
        "CF": res["CF"], "isp_ratio": res["isp_ratio"], "closure_norm": float(np.hypot(cl["dx"], cl["dr"])),
        "kappa_B": res["kappa_B"], "kappaW_B": res["kappaW_B"], "dkappa_B": res["dkappa_B"], "g_tilde": gt["g"],
        "goursat_nfail": res["goursat_nfail"], "n_wall_missing": res["n_wall_missing"],
        "streamline_resid_max": res["streamline_resid_max"]})
    return rec


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    if isinstance(v, Path):
        return str(v)
    return v


def correct_rao(problem_yaml, freeze_run, *, k_delta: float, parts_dir, tag: str, freeze_res: str = "res_30000.h5",
                n_line: int = 321, ramp=RAMP_VEXT, repo_root=None) -> dict:
    """§4.4 の手順 1〜5。problem_yaml: 補正しない (M) の問題 (`geometry.bell_mode: moc_rao`・`rao_wall_file` あり・`rao_offset_file`
    なし)。freeze_run: 凍結源の Euler run (CFD ピンの初期線と収縮部の縁)。戻り: 記録 dict (`setup` はキャンペーンの要約へ)。"""
    from ..evaluate.runner import build_wall, rao_bell_spec
    from ..geometry import rao_axisym as ra
    from ..geometry.wall import DisplacementOffset
    from ..metrics.bell_deltastar import contur_bell
    from ..probdef import load_problem
    from .cfd_initial_line import expected_bell_source, moc_initial_line, pinned_factory
    t0 = time.time()
    p = load_problem(problem_yaml)
    if p.geometry.get("bell_mode") != "moc_rao" or p.rao_wall_file is None:
        raise ValueError("correct_rao: geometry.bell_mode moc_rao と rao_wall_file (補正しない ε の壁) が要る")
    if p.rao_offset_file is not None:
        raise ValueError("correct_rao: 既に rao_offset_file がある (補正済みの問題を重ねて補正しない)")
    if float(p.spec.get("p_ambient", 0.0)) != 0.0:
        raise ValueError("correct_rao: 生成器は真空 (p_a = 0) の設計点だけ")
    wall, wrec = build_wall(p)
    spec = rao_bell_spec(p)
    L, eps = float(wall.x_e), float(spec.eps)
    xo = np.linspace(float(wall.x_in) + 1e-6, L, N_THEORY)
    cb = contur_bell(problem_yaml, freeze_run, x_out=xo, edge="moc", moc_wall_npz=p.rao_wall_file)
    x, C = np.asarray(cb["x"], dtype=float), np.asarray(cb["delta_r_ann"], dtype=float)
    if not np.all(np.isfinite(C[x >= 0.0])):
        raise ValueError("correct_rao: 理論式の δ_r に非有限 (x ≥ 0)")
    if abs(x[-1] - L) > 1e-6 * L:
        raise ValueError(f"correct_rao: 理論式の終点 {x[-1]} が壁の終点 {L} と違う")
    kd = float(k_delta)
    dE = kd * float(C[-1])
    eps_inv = (np.sqrt(eps) - dE) ** 2
    gas = ra.as_isentrope(p.exhaust.moc_gas())
    eb = expected_bell_source(p)
    th = pinned_factory(Path(freeze_run), freeze_res, eb["expect"])(eb["R"], float(gas.moc.gamma_throat()))
    line = moc_initial_line(th, gas.moc, int(n_line))
    t1 = time.time()
    res = ra.design_rao(gas, R_t=float(spec.R_t), L=L, eps=eps_inv, kernel="arc", line=line, p_a=0.0)
    parts = Path(parts_dir); parts.mkdir(parents=True, exist_ok=True)
    rec = wall_record(res, tag, L, time.time() - t1,
                      note=f"δ* 補正の非粘性壁 (ε_inv = (√{eps:g} − kδ·δ_r,CONTUR(x_E))²、kδ {kd})。元の壁 {p.rao_wall_file.name}")
    if res["status"] != "ok":
        (parts / f"{tag}.json").write_text(json.dumps(_jsonable(rec), ensure_ascii=False, indent=1))
        raise ValueError(f"correct_rao: ε_inv {eps_inv} の (M) が作れない ({res.get('why')})")
    (parts / f"{tag}.json").write_text(json.dumps(_jsonable(rec), ensure_ascii=False, indent=1))
    np.savez(parts / f"{tag}_wall.npz", BE=res["BE"], kernel_wall=res["kernel_wall"], wall=res["wall"])
    m = x > 0.0
    xs = np.r_[0.0, x[m]]; xs[-1] = L
    ds = kd * np.r_[np.interp(0.0, x, C), C[m]]
    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[3]
    meta = {"k_delta": kd, "eps_phys": eps, "eps_inv": float(eps_inv), "delta_E": dE, "ramp": list(ramp),
            "variant": "vext" if tuple(ramp) == RAMP_VEXT else "custom", "x_cal_start": X_CAL_START,
            "theory": "contur_v2 (metrics.bell_deltastar.contur_bell, edge moc)", "edge": _jsonable(cb.get("edge")),
            "source_wall": str(p.rao_wall_file), "source_wall_sha256_16": _sha16(p.rao_wall_file),
            "inviscid_wall": f"{tag}_wall.npz", "problem": str(problem_yaml),
            "plan": "plans/active/tooling-nozzle-bell-two-mode.md §4.4・§5.1 #7c"}
    off = DisplacementOffset.fit(xs, ds, ramp=tuple(ramp), meta=meta)
    off.save(parts / f"{tag}_offset.npz")
    exh = Path(p.raw["gas"]["exhaust"])
    if not exh.is_absolute():
        exh = Path(p.path).resolve().parent / exh
    try:
        exh_rel = str(exh.resolve().relative_to(root.resolve()))
    except ValueError:
        exh_rel = str(exh.resolve())
    setup = {"gas": "EXH (RocketExhaust.moc_gas)", "exhaust": exh_rel, "R_t": float(spec.R_t), "eps": float(eps_inv),
             "eps_phys": eps, "p_a": 0.0, "kernel": "arc", "accept_rule": "det_D > 0", "L_rt": L, "n_line": int(n_line),
             "initial_line": f"CFD ピン ({Path(freeze_run).name}/{freeze_res}、moc_initial_line n {n_line}・r 等間隔)",
             "deltastar_correction": {"k_delta": kd, "delta_E": dE, "ramp": list(ramp), "offset": f"{tag}_offset.npz"}}
    return {"tag": tag, "setup": setup, "wall_record": {k: rec.get(k) for k in ("theta_B_deg", "theta_E_deg", "theta_max_deg", "x_E",
                                                                              "r_E", "M_E", "F_cs", "thrust_rel_diff", "closure_norm",
                                                                              "det_D")},
            "eps": eps, "eps_inv": float(eps_inv), "delta_E": dE, "k_delta": kd, "offset": off.record(),
            "theory": {"x": x.tolist(), "delta_r": C.tolist()}, "time_s": time.time() - t0,
            "files": {"wall": str(parts / f"{tag}_wall.npz"), "record": str(parts / f"{tag}.json"),
                      "offset": str(parts / f"{tag}_offset.npz")}}


def write_campaign(parts_dir, results: list, extra: dict | None = None) -> Path:
    """`<parts_dir の親>/<prefix>.json` (prefix = parts_dir の名前から `_parts` を除いたもの) に `setup_by_tag` を足す
    (`load_rao_wall` が tag ごとの setup を読む)。既存の要約は読んで合わせる。"""
    parts = Path(parts_dir)
    if not parts.name.endswith("_parts"):
        raise ValueError("write_campaign: parts_dir の名前は <prefix>_parts")
    cpath = parts.parent / (parts.name[:-len("_parts")] + ".json")
    c = json.loads(cpath.read_text()) if cpath.exists() else {}
    sbt = dict(c.get("setup_by_tag") or {})
    for r in results:
        sbt[r["tag"]] = r["setup"]
    c["setup_by_tag"] = sbt
    c.setdefault("setup", next(iter(sbt.values())))
    c["rows"] = {**(c.get("rows") or {}), **{r["tag"]: {k: v for k, v in r.items() if k not in ("setup", "theory")} for r in results}}
    if extra:
        c.update(extra)
    cpath.write_text(json.dumps(_jsonable(c), ensure_ascii=False, indent=1))
    return cpath
