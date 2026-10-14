#!/usr/bin/env python3
r"""③ベル (M) モード: CFD ピンの初期線 (run_0015、n 321) で L_k 70・80・90 % の (M) を生成する (L100 は `run_pin_line.py design`)。
記録、合否は親と上位。

計画: plans/active/tooling-nozzle-bell-two-mode.md §5.1 #6 (W3 の続き)・#8 (比較の (M) 壁)。ピン版 V2 は L100 で PASS (run_0018)。

段 (サブコマンド):

- `design` (並列 2 本、各 20〜25 分): L70・L80・L90 (+ 任意の L60 — Hall n 161 では「(M) で作れない L_k」[根なし]) を
  `run_pin_line.py` の L100 と同じ設定 (ガス `RocketExhaust.moc_gas`、R_t 1、ε 50 固定、円弧 kernel、初期線 = CFD ピン
  `moc_initial_line(run_0015/res_30000.h5, n 321)`、Goursat m = q = n = 321・dr 3.2/n の既定) で `design_rao` に通し、
  `rao_pin_parts/L{70,80,90}_pin_n321{.json,_wall.npz}` に同じ書式 (`run_rao_t6._record`) で書く。L60 は記録 (json) だけ。
  **既にある tag は回さず、上書きしない** (止めて理由を出す)。
- `summary` (数十秒): 各本の表 (θ_B・θ_E・θ_max・x_E・det_D・推力 2 会計の差・閉合・Δκ_B・g̃・η_cs・所要時間) と Hall 版との
  θ_B・θ_E の差 (`rao_t6_parts/L70_n321` [n 321]・`L100_n321`、無い L_k は `rao_lk_n161.json` [n 161 — 分解能の差は L70・L100 で
  n 161 と n 321 の Hall の差として併記])。キャンペーンの要約: `rao_pin.json` の既存のキーは変えず、`setup_by_tag` (L_k ごとの
  setup — L100 の `setup` と同じ形で L_rt・L_frac だけが違う) と `lk` (表へのポインタ) を**追記**する。`load_rao_wall` は
  `setup_by_tag` にその壁の tag があればそれを campaign の setup にする (無い壁 [L100・t6] は従来どおり `setup`)。
  各壁が `load_rao_wall` で読め、setup の L・ε・R_t と壁 (x_E・r_E²・1/κ_B) と記録 (L_rt) がその L_k と合うことを確かめる。
  出力 `rao_pin_lk.json`・`rao_pin_lk.csv`。

実行 (手元 PC はユーザの作業機 — 並列 2 本まで、nice):
`nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_pin_lk.py design|summary`
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_pin_line as P  # noqa: E402  (ガス・凍結源・ピン線・記録の書式・保存を共用)

from forge_design.geometry import rao_axisym as ra  # noqa: E402

FRACS = {"L70": 0.7, "L80": 0.8, "L90": 0.9}
FRAC_OPTIONAL = {"L60": 0.6}                                  # 任意: 根なしの確認 (壁は書かない)
LK_JSON = HERE / "rao_pin_lk.json"
LK_CSV = HERE / "rao_pin_lk.csv"
HALL_N321 = {"L70": HERE / "rao_t6_parts" / "L70_n321.json", "L100": HERE / "rao_t6_parts" / "L100_n321.json"}
HALL_N161 = HERE / "rao_lk_n161.json"
TAG_FMT = "{lab}_pin_n321"


def tag_of(lab: str) -> str:
    return TAG_FMT.format(lab=lab)


def L_of(frac: float) -> float:
    return frac * P.L100


def job(lab: str, frac: float) -> dict:
    t0 = time.time()
    tag = tag_of(lab)
    if (P.PARTS / f"{tag}.json").exists() or (P.PARTS / f"{tag}_wall.npz").exists():
        raise FileExistsError(f"{tag} は既にある (上書きしない)")
    ig = P.gas()
    th, match = P.pinned_throat(P.FREEZE, ig)
    line = P.moc_initial_line(th, ig.moc, P.N_LINE)
    L = L_of(frac)
    spec = {"tag": tag, "L_frac": frac, "variant": "pin"}
    res = ra.design_rao(ig, R_t=P.R_T, L=L, eps=P.EPS, kernel="arc", line=line, p_a=P.P_A)
    rec = P.t6._record(res, spec, L, time.time() - t0)
    lr = {"source": "cfd_pinned", "run": P.FREEZE, "res": P.RES, "n": P.N_LINE, "match": match,
          **P._setup(line)["initial_line_detail"]}
    rec.update({"line_source": res.get("line_source"), "initial_line": lr,
                "note": f"CFD ピンの初期線 ({P.FREEZE}/{P.RES}、n {P.N_LINE}・r 等間隔) の (M)、L_k {lab}。Goursat m = q = n (既定)"})
    P._save(P.PARTS, rec, res if lab in FRACS else None)          # L60 (任意) は記録だけ
    print(P.t6._line(rec), flush=True)
    return rec


def cmd_design(workers: int, with_l60: bool) -> None:
    if workers > 2:
        raise SystemExit("--workers は 2 まで (手元 PC はユーザの作業機)")
    jobs = dict(FRACS)
    if with_l60:
        jobs.update(FRAC_OPTIONAL)
    for lab in jobs:
        tag = tag_of(lab)
        if (P.PARTS / f"{tag}.json").exists() or (P.PARTS / f"{tag}_wall.npz").exists():
            raise SystemExit(f"{tag} は既にある (上書きしない) — 止める")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(job, lab, frac): lab for lab, frac in jobs.items()}
        for fu in as_completed(futs):
            try:
                fu.result()
            except Exception as e:  # noqa: BLE001
                print(f"{futs[fu]}: ERROR {type(e).__name__}: {e}\n{traceback.format_exc()}", flush=True)
    print(f"design 完了 ({time.time() - t0:.0f} s)", flush=True)


# ================================================================================================ 表
def setup_for(lab: str, frac: float, base: dict) -> dict:
    """L100 の `setup` (rao_pin.json) と同じ形で、L_rt・L_frac だけをその L_k にしたもの。"""
    s = dict(base)
    s.update({"L_rt": L_of(frac), "L_frac": frac})
    return s


def hall_ref(lab: str, frac: float) -> dict:
    if lab in HALL_N321 and HALL_N321[lab].exists():
        r = json.loads(HALL_N321[lab].read_text())
        return {"source": str(HALL_N321[lab].relative_to(P.ROOT)), "n": r["n"], "theta_B_deg": r["theta_B_deg"],
                "theta_E_deg": r["theta_E_deg"], "theta_max_deg": r["theta_max_deg"]}
    d = json.loads(HALL_N161.read_text())
    r = next(x for x in d["L_k"] if abs(float(x["L_frac"]) - frac) < 1e-9)
    return {"source": str(HALL_N161.relative_to(P.ROOT)) + f" (L_frac {frac})", "n": r["n"], "status": r["status"],
            "theta_B_deg": r.get("theta_B_deg"), "theta_E_deg": r.get("theta_E_deg"), "theta_max_deg": r.get("theta_max_deg")}


def hall_n_effect() -> dict:
    """Hall の n 161 → 321 の θ_B・θ_E の差 (L70・L100 の記録; L80・L90 を n 161 と比べるときの分解能の目安)。"""
    d = json.loads(HALL_N161.read_text())
    out = {}
    for lab, frac in (("L70", 0.7), ("L100", 1.0)):
        r161 = next(x for x in d["L_k"] if abs(float(x["L_frac"]) - frac) < 1e-9)
        r321 = json.loads(HALL_N321[lab].read_text())
        out[lab] = {k: r321[k] - r161[k] for k in ("theta_B_deg", "theta_E_deg", "theta_max_deg")}
    return out


def check_load(lab: str, frac: float, setup: dict) -> dict:
    """壁ファイルが load_rao_wall で読め、campaign の setup・壁・記録の L・ε・R_t がその L_k と合うか。"""
    w, info = ra.load_rao_wall(P.PARTS / f"{tag_of(lab)}_wall.npz")
    L = L_of(frac)
    s = (info.get("campaign") or {}).get("setup") or {}
    chk = {"L_k": L, "x_end": float(w.x_end), "x_end_rel": float(w.x_end) / L - 1.0,
           "eps_wall": float(w.r_end) ** 2, "R_t": info["R_t"], "generator_L_rt": info["generator"]["L_rt"],
           "setup_L_rt": s.get("L_rt"), "setup_L_frac": s.get("L_frac"), "setup_eps": s.get("eps"), "setup_R_t": s.get("R_t"),
           "setup_source": (info.get("campaign") or {}).get("setup_source"),
           "setup_equals_expected": s == setup}
    chk["ok"] = bool(abs(chk["x_end_rel"]) <= 1e-9 and abs(chk["eps_wall"] / P.EPS - 1.0) <= 1e-9 and chk["R_t"] == P.R_T
                     and chk["generator_L_rt"] == L and s.get("L_rt") == L and s.get("L_frac") == frac
                     and s.get("eps") == P.EPS and s.get("R_t") == P.R_T and chk["setup_equals_expected"])
    chk["last_point_minus_E"] = info["last_point_minus_E"]
    chk["fit_info"] = {k: info["fit_info"][k] for k in ("resid_max", "lam", "n_dkappa_sign_changes", "theta_max_BE_deg",
                                                         "kappa_B_plus")}
    return chk


def cmd_summary() -> dict:
    ig = P.gas()
    camp = json.loads(P.OUT.read_text())
    base = camp["setup"]
    setups = {tag_of(lab): setup_for(lab, frac, base) for lab, frac in FRACS.items()
              if (P.PARTS / f"{tag_of(lab)}_wall.npz").exists()}
    # キャンペーンの要約へ追記 (既存のキーは変えない)
    sbt = dict(camp.get("setup_by_tag") or {})
    for k, v in setups.items():
        if k in sbt and sbt[k] != v:
            raise SystemExit(f"rao_pin.json の setup_by_tag[{k}] が既にあり内容が違う (上書きしない)")
        sbt[k] = v
    camp["setup_by_tag"] = sbt
    camp["lk"] = {"script": str(Path(__file__).relative_to(P.ROOT)), "table": str(LK_JSON.relative_to(P.ROOT)),
                  "tags": sorted(setups), "note": "L_k 70・80・90 % のピン版 (M)。setup_by_tag は load_rao_wall が tag ごとに使う"}
    P._dump(P.OUT, camp)
    rows, csv_rows = {}, []
    labs = list(FRACS) + ["L100"] + [lab for lab in FRAC_OPTIONAL if (P.PARTS / f"{tag_of(lab)}.json").exists()]
    fr_all = {**FRACS, "L100": 1.0, **FRAC_OPTIONAL}
    for lab in labs:
        frac = fr_all[lab]
        p = P.PARTS / f"{tag_of(lab)}.json"
        if not p.exists():
            rows[lab] = {"status": "missing"}
            continue
        r = json.loads(p.read_text())
        row = {"tag": r["tag"], "status": r["status"], "why": r.get("why"), "L_frac": frac, "L_rt": r["L_rt"], "time_s": r["time_s"]}
        if r["status"] == "ok":
            dv = P.derived(r, ig)
            row.update({k: r[k] for k in ("theta_B_deg", "theta_E_deg", "theta_max_deg", "x_E", "x_E_minus_L", "det_D", "de",
                                          "theta_E_le_theta_B", "thrust_rel_diff", "mdot_rel_diff", "closure_norm", "dkappa_B",
                                          "kappaW_B", "g_tilde", "M_E", "x_D", "eta_D", "goursat_nfail", "n_wall_missing",
                                          "F_cs", "F_wall", "isp_ratio")})
            row.update({"eta_cs": dv["eta_cs"], "eta_wall": dv["eta_wall"], "mstar_ratio_1d": dv["mstar_ratio_1d"],
                        "n_candidates": len(r.get("candidates") or []), "n_rejected": len(r.get("rejected") or []),
                        "n_lost": len(r.get("lost") or []), "n_policy_holes": len(r.get("policy_holes") or [])})
            h = hall_ref(lab, frac)
            row["hall"] = h
            if h.get("theta_B_deg") is not None:
                row["dtheta_B_deg_vs_hall"] = r["theta_B_deg"] - h["theta_B_deg"]
                row["dtheta_E_deg_vs_hall"] = r["theta_E_deg"] - h["theta_E_deg"]
                row["dtheta_max_deg_vs_hall"] = r["theta_max_deg"] - h["theta_max_deg"]
            if lab in FRACS:
                row["load_check"] = check_load(lab, frac, setups[tag_of(lab)])
        else:
            row["scan"] = {"n_rejected": len(r.get("rejected") or []), "n_lost": len(r.get("lost") or []),
                           "branches": r.get("branches"), "scan_mid_n": len(r.get("scan_mid") or [])}
        rows[lab] = row
        csv_rows.append(row)
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #6 (W3 の続き)・#8", **P.git_info(),
           "setup": {**base, "L_fracs": [fr_all[x] for x in labs], "L_rt": None, "L_frac": None,
                     "note": "L_k ごとの setup は rao_pin.json の setup_by_tag。L100 は run_pin_line.py design の記録"},
           "hall_n161_to_n321": hall_n_effect(), "rows": rows}
    P._dump(LK_JSON, out)
    cols = ("tag", "status", "L_frac", "L_rt", "theta_B_deg", "theta_E_deg", "theta_max_deg", "x_E", "x_E_minus_L", "det_D",
            "thrust_rel_diff", "mdot_rel_diff", "closure_norm", "dkappa_B", "g_tilde", "eta_cs", "eta_wall", "mstar_ratio_1d",
            "dtheta_B_deg_vs_hall", "dtheta_E_deg_vs_hall", "dtheta_max_deg_vs_hall", "hall_source", "time_s", "why")
    with open(LK_CSV, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(cols))
        wr.writeheader()
        for row in csv_rows:
            rr = {k: row.get(k) for k in cols}
            rr["hall_source"] = (row.get("hall") or {}).get("source")
            wr.writerow({k: (f"{v:.10g}" if isinstance(v, float) else ("" if v is None else v)) for k, v in rr.items()})
    for lab, row in rows.items():
        if row.get("status") != "ok":
            print(f"{lab}: {row.get('status')} — {row.get('why')}", flush=True)
            continue
        lc = row.get("load_check") or {}
        print(f"{lab}: θ_B {row['theta_B_deg']:.4f}° θ_E {row['theta_E_deg']:.4f}° θ_max {row['theta_max_deg']:.4f}° "
              f"x_E {row['x_E']:.6f} det_D {row['det_D']:+.4f} 2 会計 {row['thrust_rel_diff']:+.3e} 閉合 {row['closure_norm']:.2e} "
              f"Δκ_B {row['dkappa_B']:+.4f} g̃ {row['g_tilde']:+.4e} η_cs {row['eta_cs']:.6f} ({row['time_s']:.0f} s) | "
              f"Hall (n {row['hall']['n']}) との差 Δθ_B {row.get('dtheta_B_deg_vs_hall', float('nan')):+.4f}° "
              f"Δθ_E {row.get('dtheta_E_deg_vs_hall', float('nan')):+.4f}°"
              + (f" | load {lc.get('ok')} (setup {lc.get('setup_source')})" if lc else ""), flush=True)
    print(f"{LK_JSON.name}・{LK_CSV.name}・{P.OUT.name} (setup_by_tag {sorted(setups)})", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("design", "summary"))
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--with-l60", action="store_true", help="任意の L60 (根なしの確認) も回す")
    a = ap.parse_args()
    if a.step == "design":
        cmd_design(a.workers, a.with_l60)
    else:
        cmd_summary()


if __name__ == "__main__":
    main()
