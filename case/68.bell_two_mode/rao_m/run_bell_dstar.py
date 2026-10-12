#!/usr/bin/env python3
r"""③ベル (M) の δ_r 抽出の受理と CONTUR (`contur_v2`) の未較正の比 — plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-1)・(7b-3) 前半。

方針の正本は plan §5.1 #7 の「(7b) 較正の設計 = codex (diagnose) 2026-10-12」と記録
notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md (M1〜M5・m1)。部品は design/forge_design/metrics/bell_deltastar.py。

  extract --tag L100|L70   : (1) 旧 (列のまま) と新 (一定 x) の表、(2) 流量の帳簿 (x_E 欠けの扱い)、(3) 帯位置の感度と抽出不確かさ、
                             (4) 設定・帯位置・感度・受理の判定を JSON/npz に、(6) contur_v2 の未較正の δ_r,CONTUR と比 δ_r,NS/δ_r,CONTUR
                             (感度: 縁の上流 1D・入口 θ0 ×0.25/×4・縁の格子 ×2・rtol 1e-8)。**較正 (kδ の決定) はしない**。
                             L70 は検証用 (kδ は L100 だけで決める — 2026-10-12 主セッション)。
  series <ns_run> [--euler <run>] [--tag ...] : (5) 途中の res_<n>.h5 を順に抽出して δ_r の時系列 CSV
                             (`check_quasisteady.py --series-csv` 用)。B で回す想定。

出力 (tag L100 の例、このディレクトリ): dstar_l100_ns_run0019_constx.{json,npz}、dstar_l100_remap_compare.csv、
dstar_l100_uncert_bands.csv、contur_l100_euler_run0018.{json,npz}、contur_l100_ratio.csv。
旧配列 dstar_l100_ns_run0019.npz (列のまま・不受理) は残す。

python (手元): /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_bell_dstar.py extract --tag L100
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.metrics import bell_deltastar as bd  # noqa: E402

PRESETS = {
    "L100": dict(ns="run_0019_tp_m_ns_l100", euler="run_0018_tp_m_euler_l100_pin", problem="problem_bell_tp_m_ns_l100.yaml",
                 role="較正用 (kδ は L100 だけで決める)"),
    "L70": dict(ns="run_0020_tp_m_ns_l70", euler="run_0023_tp_m_euler_l70_pin", problem="problem_bell_tp_m_ns_l70.yaml",
                role="検証用 (較正には使わない)"),
}
SOURCES = ("design/forge_design/metrics/deltastar.py", "design/forge_design/metrics/bell_deltastar.py",
           "design/forge_design/feedback/deltastar_integral.py", "case/68.bell_two_mode/rao_m/run_bell_dstar.py")


def _sha16(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


def _sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def provenance() -> dict:
    def git(*a):
        try:
            return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True, check=True).stdout.rstrip("\n")
        except Exception as e:                                       # noqa: BLE001
            return f"(git 不可: {e})"
    return {"commit": git("rev-parse", "HEAD"),
            "dirty": [ln for ln in git("status", "--porcelain", "--", "design", "case/68.bell_two_mode/rao_m").splitlines()
                      if not ln.startswith("??") or ln.endswith(".py")],
            "source_sha256_16": {s: _sha16(ROOT / s) for s in SOURCES if (ROOT / s).exists()}}


def _num(rn: str) -> str:
    return rn.split("_")[1]


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items() if not isinstance(v, np.ndarray) or v.ndim == 0}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def wall_check(ns_dir: Path, eu_dir: Path, problem: Path) -> dict:
    """NS・Euler の run の壁 (WALL_DESIGN.json の wall_file_sha256) と手元の壁ファイルが同じか。違えば止める。"""
    import yaml
    wf = (problem.parent / yaml.safe_load(problem.read_text())["geometry"]["rao_wall_file"]).resolve()
    local = _sha(wf)
    shas = {}
    for tag, d in (("ns", ns_dir), ("euler", eu_dir)):
        shas[tag] = json.loads((d / "WALL_DESIGN.json").read_text())["wall_file_sha256"]
    if not (shas["ns"] == shas["euler"] == local):
        raise SystemExit(f"壁ファイルが合わない: 手元 {wf.name} {local[:16]} / NS {shas['ns'][:16]} / Euler {shas['euler'][:16]}")
    return {"wall_file": str(wf.relative_to(ROOT)), "sha256": local, "ns_run_matches": True, "euler_run_matches": True}


def write_csv(path: Path, rows: list, cols: list):
    with open(path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join("" if r.get(c) is None else (repr(float(r[c])) if isinstance(r[c], (float, np.floating)) else str(r[c]))
                              for c in cols) + "\n")


def cmd_extract(a):
    pr = PRESETS[a.tag]
    ns_dir, eu_dir, problem = CASE / pr["ns"], CASE / pr["euler"], CASE / pr["problem"]
    out_dir = Path(a.out_dir) if a.out_dir else HERE
    t = a.tag.lower()
    stem_x = f"dstar_{t}_ns_run{_num(pr['ns'])}_constx"
    stem_c = f"contur_{t}_euler_run{_num(pr['euler'])}"
    wchk = wall_check(ns_dir, eu_dir, problem)
    prov = provenance()

    # (1) 旧 (列のまま) と新 (一定 x)
    cmp_ = bd.remap_comparison(ns_dir, eu_dir)
    write_csv(out_dir / f"dstar_{t}_remap_compare.csv", cmp_["rows"],
              ["station", "x", "raw_column", "raw_constx", "raw_rel_change", "smooth_column", "smooth_constx", "smooth_rel_change"])
    # (2)(3)(4) 一定 x の抽出 + 帯位置の感度 + 受理
    ex = bd.extract_bell(ns_dir, eu_dir)
    u = ex["uncertainty"]
    if not np.array_equal(ex["delta_r_raw"], cmp_["raw_constx"]):
        raise SystemExit("帯位置の感度つきの抽出の主系列が一定 x の比較の値と違う (感度が主系列を変えてはいけない)")
    bcols = (["x_lo", "x_hi", "n", "u_max", "frac_over", "u_outward_max"]
             + [f"rel_f{f:g}_{s}" for f in u["factors"] for s in ("min", "max")])
    write_csv(out_dir / f"dstar_{t}_uncert_bands.csv", u["bands"], bcols)
    reasons = {}
    for r in ex["reason"]:
        for k in (r.split(",") if r else []):
            reasons[k] = reasons.get(k, 0) + 1
    rw = ex["r_wall_ns"]
    yb_st = {k: v for k, v in bd.station_values(ex["x"], ex["band_y_b"] / rw).items()}
    np.savez(out_dir / f"{stem_x}.npz", x=ex["x"], delta_r_raw=ex["delta_r_raw"], delta_r_smooth=ex["delta_r_smooth"],
             delta_r_use=ex["delta_r_use"], delta_r_sens_yb=ex["delta_r_sens_yb"], yb_factors=np.array(u["factors"]),
             band_y_b=ex["band_y_b"], band_retry=ex["band_retry"], band_rms=ex["band_rms"], band_slope=ex["band_slope"],
             band_deficit_share=ex["band_deficit_share"], r_wall_ns=ex["r_wall_ns"], r_wall_euler=ex["r_wall_euler"],
             delta_in=ex["delta_in"], alpha=ex["alpha"], core_rms=ex["core_rms"], mass_deficit=ex["mass_deficit"],
             ok=ex["ok"], hard_ok=ex["hard_ok"], reason=ex["reason"], u=u["u"], u_outward=u["u_outward"], rel_yb=u["rel"],
             delta_r_raw_column=cmp_["raw_column"], delta_r_smooth_column=cmp_["smooth_column"])
    js = {"what": "③ベル (M) の δ_r 抽出 (一定 x の断面・帯位置の感度・受理) — plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-1)",
          "tag": a.tag, "role": pr["role"], "ns_run": f"case/68.bell_two_mode/{pr['ns']}", "euler_run": f"case/68.bell_two_mode/{pr['euler']}",
          "wall": wchk, "settings": _jsonable(ex["settings"]), "massflow": _jsonable(ex["massflow"]),
          "massflow_column_note": "列のまま (列内で x が変わる) の流量の帳簿 — −ρU_r dx の項が落ちるので③では使わない (記録)",
          "massflow_column": _jsonable(ex.get("massflow_column")),
          "n_sections": int(len(ex["x"])), "n_ok": int(ex["ok"].sum()), "n_hard_ok": int(ex["hard_ok"].sum()), "reasons": reasons,
          "uncertainty": _jsonable({k: v for k, v in u.items() if k not in ("u", "u_outward", "rel")}),
          "band_y_b_over_rw_at_stations": yb_st,
          "band_retry_at_stations": bd.station_values(ex["x"], ex["band_retry"].astype(float)),
          "delta_r_raw_at_stations": bd.station_values(ex["x"], ex["delta_r_raw"]),
          "delta_r_smooth_at_stations": bd.station_values(ex["x"], ex["delta_r_smooth"]),
          "remap_compare": _jsonable({k: v for k, v in cmp_.items() if not isinstance(v, np.ndarray)}),
          "files": {"npz": f"{stem_x}.npz", "remap_csv": f"dstar_{t}_remap_compare.csv", "bands_csv": f"dstar_{t}_uncert_bands.csv"},
          "provenance": prov}
    (out_dir / f"{stem_x}.json").write_text(json.dumps(js, indent=1, ensure_ascii=False))
    print(f"[extract] {stem_x}: 断面 {len(ex['x'])}、受理 {u['accepted']} (x ≥ 0 の最大 {u['u_max']:.4f} @ x {u['x_at_u_max']:.3f}、"
          f"5 % 超 {u['n_over']}/{u['n_region']}、外向きだけの最大 {u['u_outward_max_region']:.4f})")
    if a.skip_contur:
        return
    # (6) contur_v2 (未較正)、x は NS の断面
    main = bd.contur_bell(problem, eu_dir, x_out=ex["x"])
    xc = main["x"]
    idx = np.searchsorted(ex["x"], xc - 2e-6)
    idx = np.clip(idx, 0, len(ex["x"]) - 1)
    if np.max(np.abs(ex["x"][idx] - xc)) > 2e-6:
        raise SystemExit("CONTUR の x が NS の断面に対応しない")
    d_ns, d_ns_s = ex["delta_r_raw"][idx], ex["delta_r_smooth"][idx]
    ratio = d_ns / main["delta_r_ann"]
    sens_def = {"edge_1d_upstream": dict(edge="euler_1d_upstream"),
                "theta0_x0.25": dict(theta0_m=0.25 * main["settings"]["theta0_m"]),
                "theta0_x4": dict(theta0_m=4.0 * main["settings"]["theta0_m"]),
                "edge_grid_x2": dict(edge_n_grid=8000), "rtol_1e-8": dict(rtol=1e-8)}
    sens = {}
    for k, kw in sens_def.items():
        try:
            s = bd.contur_bell(problem, eu_dir, x_out=ex["x"], **kw)
            if not np.array_equal(s["x"], xc):
                raise ValueError("x が主の CONTUR と違う")
            sens[k] = s["delta_r_ann"]
        except Exception as e:                                       # noqa: BLE001 — 感度が回らない事実を記録する
            sens[k] = None
            print(f"[contur] 感度 {k} が回らない: {e}")
    st_rows = []
    for key in list(bd.STATIONS) + ["exit"]:
        xs = float(xc[-1]) if key == "exit" else float(key)
        f = lambda v: float(np.interp(xs, xc, v))                    # noqa: E731
        row = {"station": "exit" if key == "exit" else f"{key:g}", "x": xs, "dr_ns_raw": f(d_ns), "dr_ns_smooth": f(d_ns_s),
               "dr_contur_ann": f(main["delta_r_ann"]), "dr_contur_lin": f(main["delta_r_lin"]), "ratio_raw": f(d_ns) / f(main["delta_r_ann"]),
               "ratio_smooth": f(d_ns_s) / f(main["delta_r_ann"]), "M_edge": f(main["M"]), "H": f(main["H"]), "Cf": f(main["Cf"])}
        for k, v in sens.items():
            row[f"sens_{k}_rel"] = (f(v) / f(main["delta_r_ann"]) - 1.0) if v is not None else None
        st_rows.append(row)
    ccols = list(st_rows[0].keys())
    write_csv(out_dir / f"contur_{t}_ratio.csv", st_rows, ccols)
    stats_raw = bd.ratio_stats(xc, d_ns, main["delta_r_ann"])
    stats_smooth = bd.ratio_stats(xc, d_ns_s, main["delta_r_ann"])
    sens_max = {}
    for k, v in sens.items():
        if v is None:
            sens_max[k] = None
            continue
        rel = v / main["delta_r_ann"] - 1.0
        sens_max[k] = {f"x>={lo:g}": float(np.max(np.abs(rel[xc >= lo]))) for lo in (0.0, 0.5, 1.0, 2.0)}
    np.savez(out_dir / f"{stem_c}.npz", x=xc, dstar_n=main["dstar_n"], delta_r_lin=main["delta_r_lin"], delta_r_ann=main["delta_r_ann"],
             M=main["M"], Te=main["Te"], Taw=main["Taw"], H=main["H"], N=main["N"], Cf=main["Cf"], theta=main["theta"],
             delta=main["delta"], cos_phi=main["cos_phi"], r_wall=main["r_wall"], dr_ns_raw=d_ns, dr_ns_smooth=d_ns_s, ratio_raw=ratio,
             **{f"sens_{k}": v for k, v in sens.items() if v is not None})
    jc = {"what": "③ベル (M) の CONTUR contur_v2 の未較正 δ_r と比 δ_r,NS/δ_r,CONTUR — plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-3) 前半。"
                  "較正 (kδ の決定) はしていない。ref_best_const は区間の参考の数値で較正値ではない",
          "tag": a.tag, "role": pr["role"], "problem": f"case/68.bell_two_mode/{pr['problem']}",
          "euler_run": f"case/68.bell_two_mode/{pr['euler']}", "ns_extraction": f"{stem_x}.json", "wall": wchk,
          "settings": _jsonable(main["settings"]), "solve_ivp": _jsonable(main["solve_ivp"]), "edge": _jsonable(main["edge"]),
          "definition": "δ_r,CONTUR = r_w − √(r_w² − 2 r_w δ*_n / cos φ_w) (一定 x の環状の質量等価; M5)。delta_r_lin = δ*_n/cos φ_w (従来の 1 次)",
          "closure_notes": ["回復係数 r = Pr^(1/3)、Pr 0.72 (EdgeConditions の既定; 燃焼室の凍結 Pr ≈ μ c_p/k = 1.0675e-4·3784/0.570 ≈ 0.709)",
                            "速度分布 (z/δ)^(1/N)、N(R_δ) は CONTUR Fig.10、C_f は Eq.75 + F_c (CONTUR)",
                            "入口 θ0 は仮想平板 x_v = L_pipe·r_t (integral_bl の文書の既定; NS の入口は一様流の inlet_Pressure)"],
          "stations": st_rows, "ratio_stats_raw": stats_raw, "ratio_stats_smooth": stats_smooth,
          "sensitivity_max_rel": sens_max, "files": {"npz": f"{stem_c}.npz", "ratio_csv": f"contur_{t}_ratio.csv"},
          "provenance": prov}
    (out_dir / f"{stem_c}.json").write_text(json.dumps(_jsonable(jc), indent=1, ensure_ascii=False))
    print(f"[contur] {stem_c}: 区間 x≥1 の比 中央値 {stats_raw[1]['median']:.4f} 最小 {stats_raw[1]['min']:.4f} 最大 {stats_raw[1]['max']:.4f}")


def cmd_series(a):
    ns_dir = Path(a.ns_run).resolve()
    if a.euler:
        eu_dir = Path(a.euler).resolve()
    elif a.tag:
        eu_dir = CASE / PRESETS[a.tag]["euler"]
    else:
        raise SystemExit("--euler か --tag を与える")
    out_csv = Path(a.out) if a.out else ns_dir / "dstar_series.csv"
    out_npz = out_csv.with_suffix(".npz")
    r = bd.deltastar_series(ns_dir, eu_dir, out_csv, min_step=a.min_step, out_npz=out_npz)
    meta = {"what": "δ_r の時系列 (M3) — plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-1)", "ns_run": str(ns_dir), "euler_run": str(eu_dir),
            "n_snapshots": len(r["rows"]), "steps": [x["step"] for x in r["rows"]], "cols": r["cols"], "command": r["command"],
            "provenance": provenance()}
    out_csv.with_suffix(".json").write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    print(f"[series] {len(r['rows'])} 枚 → {out_csv}")
    print("判定:", r["command"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    e = sp.add_parser("extract")
    e.add_argument("--tag", choices=sorted(PRESETS), required=True)
    e.add_argument("--out-dir", default=None, help="出力先 (既定 このディレクトリ)")
    e.add_argument("--skip-contur", action="store_true")
    s = sp.add_parser("series")
    s.add_argument("ns_run")
    s.add_argument("--euler", default=None, help="Euler 基準の run (既定は --tag の preset)")
    s.add_argument("--tag", choices=sorted(PRESETS), default=None)
    s.add_argument("--out", default=None, help="CSV (既定 <ns_run>/dstar_series.csv; 同名の .npz・.json も書く)")
    s.add_argument("--min-step", type=int, default=None, help="この step 以上のスナップショットだけ")
    a = ap.parse_args(argv)
    {"extract": cmd_extract, "series": cmd_series}[a.cmd](a)


if __name__ == "__main__":
    main()
