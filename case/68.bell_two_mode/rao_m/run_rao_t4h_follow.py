#!/usr/bin/env python3
r"""③ベル (M) モード §5.1 #5 (h′) — T4 (停留性) の follow: (h) の A/B を登録の A と読んだ後の k2・k3・L70・負対照と、
k₀ を 1 本下流に変えた k1 (数値を出すだけ、判定は親と上位)。

plan: plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (h′)・§6 V1 の T4。
判断の記録: notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md の「追加諮問 2」。
(h) の台本 `run_rao_t4h.py` (壁の読み込み・前線・当てはめ) をそのまま使う (変えない)。

**初期線**: (h) と同じ — L100 の (M) 壁 (`rao_t6_parts/L100_n321_mq2`) と Hall n 321 で前進 MOC を回した網の前線 k₀ = 64
(壁点 x 0.3006) を r 等間隔の N + 1 = 1282 点へ 3 次補間 (`rao_axisym.resampled_front_line`)。**全ジョブで同じ線**を使い、
各ジョブの壁でも Hall n 321 のマーチを回して前線 k₀ が線の元の前線とビット一致すること (sha) と、k₀ の壁点がその壁の
x_B より上流であることを先に確かめる (一致しなければそのジョブは回さない)。

**ジョブ** (族 sin_c1、ε ∈ {0, ±1, ±2, ±4}×10⁻³ r_t の 7 点、スプラインと台形の両会計、`forward_thrust(..., splits=(x_B,))`):
- `L100_k23`: L100 (`L100_n321_mq2`) の k2・k3。k1 は (h) の `ab_L100_resampled_N1281` (`rao_t4h.json`) を表に載せる。
- `L70_k123`: L70 (`L70_n321_mq2`) の k1・k2・k3。
- `ctrl_L100_n161`・`ctrl_L70_n161` (負対照): 既存の `rao_t4_parts/gen_L100_n161`・`gen_L70_n161` の `ctrl_plus` 壁 (生成器は
  回し直さない。壁の作り方は前任の `run_rao_t4.base_wall` そのもの — σ_D ×1.02 の D⁺ の Rao 線の壁を x = L で切る) の k1。
- `ref_M_L100_n161`・`ref_M_L70_n161` (参照、登録外): 同じ生成器の (M) の壁 (`base_wall(..., "M")`) の k1。負対照の比を
  同じ生成器の壁どうしで見るため (登録の比は mq2 の (M) と n161 の対照)。
- `k0down_L100_k1` (任意): 前線を 1 本下流 (壁点 x ≈ 0.45、L100 の x_B 0.4963 より上流) に変えた N 1281 の線で L100 の k1。
  c₁ε₀ が k₀ 64 の値 (`rao_t4h.json`) から 1e-8 以上動くかを記録。

**判定に使う数値** (判定は親): 各 k の c₁ε₀ (両会計)・2 会計の差・c₂ε₀²・比 |c₁ε₀|/|c₂ε₀²| (線 0.05)・ε*・懸かる推力 / F_1D。
負対照: |c₁ε₀(M)|/|c₁ε₀(対照)| (線 5 %) と、対照の |c₁ε₀| が床 (2 会計差) の何倍か (比の線は 10 倍以上のものだけ — 床は対照
自身の 2 会計差と、(h) の L100 N 1281 の 3.41e-9 の両方で出す)。

出力 (このディレクトリ): `rao_t4h_follow_parts/` (ジョブごと、再開用)、`rao_t4h_follow.csv` (ジョブ × k の 1 行)、
`rao_t4h_follow.json` (全 F(ε)・当てはめ・前線の一致の検査・判定の数値・来歴)。既存の `rao_t4h.*` は上書きしない。
実行: `nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_rao_t4h_follow.py`
(2 並列まで)。`--only` (ジョブ名の部分一致)、`--assemble-only`。
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_rao_t4 as t4  # noqa: E402  (前任の base_wall)
import run_rao_t4h as h  # noqa: E402  ((h) の壁・前線・当てはめ)
from forge_design.geometry import moc_forward as mf  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402

PARTS = HERE / "rao_t4h_follow_parts"
N_LINE = 1281
K0_X = 0.3                            # (h) と同じ前線 (k₀ = 64)
K0_X_DOWN = 0.45                      # 任意: 1 本下流の前線の壁点
LINE_RATIO = 0.05                     # |c₁ε₀| ≤ 0.05 |c₂ε₀²|・|c₁(M)| ≤ 0.05 |c₁(対照)|
CTRL_FLOOR_MULT = 10.0                # 比の線は対照の |c₁ε₀| ≥ 10 × 床 のものだけ
DK0_LINE = 1e-8                       # 任意: k₀ を変えて c₁ε₀ が動いたとみなす幅
LINE_WALL = "L100"                    # 線の元の前線を作る壁


def _front_of(w, L, x_target: float):
    """壁 w で Hall n 321 の前進 MOC を L まで回し、壁点が x_target に最も近い前線 (網・k・sha・情報)。"""
    res = mf.forward_moc(ra.ExtendedWall(w, L), L, h.gas().moc, h.hall(h.N_HALL))
    net, n = res["net"], res["n"]
    xw = net.a["x"][:, n]
    k = int(np.nanargmin(np.abs(xw - x_target)))
    sha = hashlib.sha256(b"".join(np.ascontiguousarray(net.a[f][k]).tobytes() for f in net.FIELDS)
                         + np.float64(net.thr_axis[k]).tobytes()).hexdigest()[:16]
    return net, k, sha, {"k0": k, "x_wall": float(xw[k]), "x_axis": float(net.a["x"][k, 0]),
                         "all_ok_to_axis": bool(np.all(net.status[k] == mf.ST_OK)), "front_sha256_16": sha,
                         "thr_axis": float(net.thr_axis[k])}


_LINES: dict = {}


def line_for(x_target: float):
    """線の元の前線 (L100 の (M) 壁) を N_LINE に再標本化した線と、その前線の情報。"""
    if x_target not in _LINES:
        w, info = h.load_wall(LINE_WALL)
        net, k, sha, fi = _front_of(w, info["L_wall"], x_target)
        ln = ra.resampled_front_line(net, k, N_LINE)
        _LINES[x_target] = (ln, {**fi, "x_B_line_wall": info["x_B"], "line_wall": LINE_WALL, "N": N_LINE})
    return _LINES[x_target]


def wall_for(spec: dict):
    if spec["wall_src"] == "t6":
        return h.load_wall(spec["wall"])
    w, info = t4.base_wall(spec["gen"], spec["which"])
    return w, info


def job(spec: dict) -> dict:
    t0 = time.time()
    ig = h.gas()
    w, info = wall_for(spec)
    L = float(info["L_wall"])
    line, lfi = line_for(spec["k0_x"])
    # 前線の一致: この壁でも Hall n 321 の前線 k₀ が線の元の前線とビット一致し、壁点がこの壁の x_B より上流
    _, k_w, sha_w, fi_w = _front_of(w, L, spec["k0_x"])
    same = (sha_w == lfi["front_sha256_16"]) and (k_w == lfi["k0"])
    upstream = lfi["x_wall"] < float(w.x_B)
    lw = mf._Wall(ra.ExtendedWall(w, L), L)
    codes, reasons, linfo = mf._check_line(line, lw, ig.moc)
    chk = {"front_same_as_line": bool(same), "front_wall": fi_w, "front_line": lfi, "x_B_wall": float(w.x_B),
           "line_wall_upstream_of_x_B": bool(upstream), "check_line_codes": codes, "check_line_info": linfo}
    if not (same and upstream and not codes):
        return {**spec, "status": "error", "why": f"前線の一致・上流・線の検査のどれかが不可 ({same}, {upstream}, {codes} {reasons})",
                "info": info, "front_check": chk, "time_s": time.time() - t0}
    base = ra.forward_thrust(w, L, ig, line, splits=(w.x_B,), gate=True)
    if not base.get("ok"):
        return {**spec, "status": "error", "why": f"基準の壁の推力: {base.get('why')}", "info": info, "front_check": chk}
    fam = {}
    for k in spec["ks"]:
        Fs, Ft, recs = [], [], []
        for e in h.EPSS:
            r = base if e == 0.0 else ra.forward_thrust(ra.PerturbedWall(w, h.FAMILY, k, e), L, ig, line, splits=(w.x_B,), gate=False)
            Fs.append(r.get("F") if r.get("ok") else None)
            Ft.append(r.get("F_trap") if r.get("ok") else None)
            recs.append({"eps": e, "ok": bool(r.get("ok")), "why": r.get("why"), "F": r.get("F"), "F_trap": r.get("F_trap"),
                         "T_wall": r.get("T_wall"), "dod_clean": r.get("dod_clean"),
                         "dod_n_bad": (r.get("dod") or {}).get("n_bad_in_dod"), "end_reason": r.get("end_reason"),
                         "n_nodes_le_L": r.get("n_nodes_le_L")})
        fam[f"{h.FAMILY}_k{k}"] = {"family": h.FAMILY, "k": k, "eps": list(h.EPSS), "F": Fs, "F_trap": Ft, "runs": recs,
                                   "fit": h.fit_cubic(h.EPSS, Fs), "fit_trap": h.fit_cubic(h.EPSS, Ft)}
    return {**spec, "status": "ok", "info": {k: v for k, v in info.items() if k != "fit_info"} | {"fit_info": info.get("fit_info")},
            "front_check": chk, "time_s": time.time() - t0,
            "base": {"F": base["F"], "F_trap": base["F_trap"], "T_line": base["T_line"], "T_wall": base["T_wall"],
                     "T_wall_trap": base["T_wall_trap"], "mdot_line": base["mdot_line"], "gate": base["gate"],
                     "dod": base["dod"], "dod_clean": base["dod_clean"], "n_nodes_le_L": base["n_nodes_le_L"],
                     "L_eval": base["L_eval"]},
            "families": fam}


def jobs() -> list:
    return [
        {"name": "L100_k23", "role": "M", "wall_src": "t6", "wall": "L100", "k0_x": K0_X, "ks": [2, 3]},
        {"name": "L70_k123", "role": "M", "wall_src": "t6", "wall": "L70", "k0_x": K0_X, "ks": [1, 2, 3]},
        {"name": "ctrl_L100_n161", "role": "ctrl", "wall_src": "gen", "gen": "gen_L100_n161", "which": "ctrl_plus", "k0_x": K0_X,
         "ks": [1]},
        {"name": "ctrl_L70_n161", "role": "ctrl", "wall_src": "gen", "gen": "gen_L70_n161", "which": "ctrl_plus", "k0_x": K0_X,
         "ks": [1]},
        {"name": "ref_M_L100_n161", "role": "ref", "wall_src": "gen", "gen": "gen_L100_n161", "which": "M", "k0_x": K0_X, "ks": [1]},
        {"name": "ref_M_L70_n161", "role": "ref", "wall_src": "gen", "gen": "gen_L70_n161", "which": "M", "k0_x": K0_X, "ks": [1]},
        {"name": "k0down_L100_k1", "role": "k0down", "wall_src": "t6", "wall": "L100", "k0_x": K0_X_DOWN, "ks": [1]},
    ]


def run(spec: dict) -> dict:
    try:
        out = job(spec)
    except Exception as e:  # noqa: BLE001
        out = {**spec, "status": "error", "why": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()}
    PARTS.mkdir(exist_ok=True)
    (PARTS / f"{spec['name']}.json").write_text(json.dumps(h._jsonable(out), ensure_ascii=False, indent=1))
    print(summary_line(out), flush=True)
    return out


def summary_line(o: dict) -> str:
    if o.get("status") != "ok":
        return f"{o['name']}: {o.get('status')} — {o.get('why')}"
    fc = o["front_check"]
    s = (f"{o['name']}: F0 {o['base']['F']:.10f} (台形 {o['base']['F_trap']:.10f}) 前線 k₀ {fc['front_line']['k0']} "
         f"(x {fc['front_line']['x_wall']:.4f}) 一致 {fc['front_same_as_line']} 依存域 clean {o['base']['dod_clean']} ({o['time_s']:.0f} s)")
    for key, v in o["families"].items():
        for acc, f in (("spline", v["fit"]), ("trap", v["fit_trap"])):
            s += (f"\n   {key} {acc}: " + (f"c1ε0 {f['c1_eps0']:+.3e} c2ε0² {f['c2_eps0sq']:+.3e} 比 {f['ratio_c1_c2']:.4f} "
                                           f"ε* {f['eps_star']:.2e}" if f.get("ok") else f"不可 {f.get('why')}"))
    return s


# ================================================================================================ 集計
ROW_KEYS = ("c1_eps0", "c2_eps0sq", "c3_eps0cu", "ratio_c1_c2", "eps_star", "hung_over_F1D", "fit_rms")


def _row(name: str, role: str, wall_label: str, k0x: float, o_fam: dict, base: dict, extra: dict) -> dict:
    row = {"job": name, "role": role, "wall": wall_label, "k0_x": k0x, "N": N_LINE, "family": o_fam["family"], "k": o_fam["k"],
           "F0": base["F"], "F0_trap": base["F_trap"], "base_dod_clean": base["dod_clean"],
           "n_pert_dod_clean": sum(1 for r in o_fam["runs"] if r.get("dod_clean")), "n_pert": len(o_fam["runs"]), **extra}
    for acc, f in (("", o_fam["fit"]), ("_trap", o_fam["fit_trap"])):
        if f.get("ok"):
            row.update({f"{k}{acc}": f[k] for k in ROW_KEYS})
            row[f"c2_sign{acc}"] = f["c2_sign"]
            row[f"ratio_le_{LINE_RATIO:g}{acc}"] = bool(f["ratio_c1_c2"] is not None and f["ratio_c1_c2"] <= LINE_RATIO)
    if o_fam["fit"].get("ok") and o_fam["fit_trap"].get("ok"):
        row["acc_diff_c1_eps0"] = abs(o_fam["fit"]["c1_eps0"] - o_fam["fit_trap"]["c1_eps0"])
    return row


def assemble(total_s=None) -> dict:
    parts = {p.stem: json.loads(p.read_text()) for p in sorted(PARTS.glob("*.json"))} if PARTS.exists() else {}
    hj = json.loads((HERE / "rao_t4h.json").read_text())
    hab = hj["parts"]["ab_L100_resampled_N1281"]
    rows = []
    # (h) の L100 k1 (N 1281、k₀ 64) を表の先頭に (再計算しない、rao_t4h.json から)
    rows.append(_row("h:ab_L100_resampled_N1281", "M", "L100_n321_mq2", K0_X, hab["families"]["sin_c1_k1"], hab["base"],
                     {"source": "rao_t4h.json (h)", "wall_info": "L100_n321_mq2"}))
    for name, o in parts.items():
        if o.get("status") != "ok":
            rows.append({"job": name, "role": o.get("role"), "status": o.get("status"), "why": o.get("why")})
            continue
        label = o["wall"] + "_n321_mq2" if o["wall_src"] == "t6" else f"{o['gen']}:{o['which']}"
        fc = o["front_check"]
        for key, v in o["families"].items():
            rows.append(_row(name, o["role"], label, o["k0_x"], v, o["base"],
                             {"source": "rao_t4h_follow", "k0": fc["front_line"]["k0"], "k0_x_wall": fc["front_line"]["x_wall"],
                              "front_same_as_line": fc["front_same_as_line"], "x_B_wall": fc["x_B_wall"],
                              "L_wall": o["info"]["L_wall"], "time_s": o["time_s"]}))
    # 判定に使う数値
    def get(job_, k):
        for r in rows:
            if r.get("job") == job_ and r.get("k") == k and "c1_eps0" in r:
                return r
        return None
    floor_h = hab["families"]["sin_c1_k1"]
    floor_h = abs(floor_h["fit"]["c1_eps0"] - floor_h["fit_trap"]["c1_eps0"])
    per_k = []
    for jb, k in (("h:ab_L100_resampled_N1281", 1), ("L100_k23", 2), ("L100_k23", 3), ("L70_k123", 1), ("L70_k123", 2), ("L70_k123", 3)):
        r = get(jb, k)
        if r is None:
            per_k.append({"job": jb, "k": k, "why": "無い / ok でない"})
            continue
        per_k.append({"job": jb, "wall": r["wall"], "k": k, "c1_eps0": r["c1_eps0"], "c1_eps0_trap": r["c1_eps0_trap"],
                      "acc_diff": r["acc_diff_c1_eps0"], "c2_eps0sq": r["c2_eps0sq"], "c2_eps0sq_trap": r["c2_eps0sq_trap"],
                      "line_0.05_c2eps0sq": LINE_RATIO * abs(r["c2_eps0sq"]),
                      "ratio": r["ratio_c1_c2"], "ratio_trap": r["ratio_c1_c2_trap"],
                      "ratio_le_0.05": r[f"ratio_le_{LINE_RATIO:g}"], "ratio_le_0.05_trap": r[f"ratio_le_{LINE_RATIO:g}_trap"],
                      "acc_diff_le_line": bool(r["acc_diff_c1_eps0"] <= LINE_RATIO * abs(r["c2_eps0sq"])),
                      "eps_star": r["eps_star"], "eps_star_trap": r["eps_star_trap"],
                      "hung_over_F1D": r["hung_over_F1D"], "hung_over_F1D_trap": r["hung_over_F1D_trap"],
                      "c2_neg": bool(r["c2_sign"] < 0 and r["c2_sign_trap"] < 0)})
    ctrl = []
    for cj, mjob, refjob in (("ctrl_L100_n161", "h:ab_L100_resampled_N1281", "ref_M_L100_n161"),
                             ("ctrl_L70_n161", "L70_k123", "ref_M_L70_n161")):
        c, m, rf = get(cj, 1), get(mjob, 1), get(refjob, 1)
        if c is None or m is None:
            ctrl.append({"ctrl": cj, "why": "対照か (M) が無い"})
            continue
        d = {"ctrl": cj, "M": mjob, "ctrl_c1_eps0": c["c1_eps0"], "ctrl_c1_eps0_trap": c["c1_eps0_trap"],
             "ctrl_acc_diff": c["acc_diff_c1_eps0"], "floor_h_L100_N1281": floor_h,
             "ctrl_over_own_floor": abs(c["c1_eps0"]) / c["acc_diff_c1_eps0"] if c["acc_diff_c1_eps0"] > 0 else None,
             "ctrl_over_floor_h": abs(c["c1_eps0"]) / floor_h,
             "ctrl_over_own_floor_trap": abs(c["c1_eps0_trap"]) / c["acc_diff_c1_eps0"] if c["acc_diff_c1_eps0"] > 0 else None,
             "ctrl_over_floor_h_trap": abs(c["c1_eps0_trap"]) / floor_h,
             "M_over_ctrl": abs(m["c1_eps0"]) / abs(c["c1_eps0"]), "M_over_ctrl_trap": abs(m["c1_eps0_trap"]) / abs(c["c1_eps0_trap"])}
        d["ctrl_ge_10x_floor (自身の床・両会計)"] = bool(d["ctrl_over_own_floor"] is not None and d["ctrl_over_own_floor"] >= CTRL_FLOOR_MULT
                                                    and d["ctrl_over_own_floor_trap"] >= CTRL_FLOOR_MULT)
        d["ctrl_ge_10x_floor ((h) の床・両会計)"] = bool(d["ctrl_over_floor_h"] >= CTRL_FLOOR_MULT and d["ctrl_over_floor_h_trap"] >= CTRL_FLOOR_MULT)
        d["M_over_ctrl_le_0.05 (両会計)"] = bool(d["M_over_ctrl"] <= LINE_RATIO and d["M_over_ctrl_trap"] <= LINE_RATIO)
        if rf is not None:
            d.update({"ref_M_n161_c1_eps0": rf["c1_eps0"], "ref_M_n161_c1_eps0_trap": rf["c1_eps0_trap"],
                      "ref_M_n161_acc_diff": rf["acc_diff_c1_eps0"],
                      "refM_over_ctrl (同じ生成器)": abs(rf["c1_eps0"]) / abs(c["c1_eps0"]),
                      "refM_over_ctrl_trap (同じ生成器)": abs(rf["c1_eps0_trap"]) / abs(c["c1_eps0_trap"])})
        ctrl.append(d)
    kd = get("k0down_L100_k1", 1)
    m64 = get("h:ab_L100_resampled_N1281", 1)
    k0down = None
    if kd is not None and m64 is not None:
        k0down = {"k0": kd.get("k0"), "k0_x_wall": kd.get("k0_x_wall"), "c1_eps0": kd["c1_eps0"], "c1_eps0_trap": kd["c1_eps0_trap"],
                  "dc1_eps0_vs_k0_64": kd["c1_eps0"] - m64["c1_eps0"], "dc1_eps0_vs_k0_64_trap": kd["c1_eps0_trap"] - m64["c1_eps0_trap"],
                  "moved_ge_1e-8": bool(abs(kd["c1_eps0"] - m64["c1_eps0"]) >= DK0_LINE or abs(kd["c1_eps0_trap"] - m64["c1_eps0_trap"]) >= DK0_LINE),
                  "F0_vs_k0_64": kd["F0"] - m64["F0"], "acc_diff": kd["acc_diff_c1_eps0"]}
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (h′)・§6 V1 T4",
           "record": "notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md (追加諮問 2)",
           "method": {"line": f"L100_n321_mq2 + Hall n {h.N_HALL} の前進 MOC の前線 (壁点 x ≈ {K0_X}、k₀ 64) を N {N_LINE} に再標本化 "
                              "(全ジョブで同じ線、各壁の前線とビット一致を確認)",
                      "thrust": "rao_axisym.forward_thrust (splits = (x_B,))、会計 = スプライン (F) と台形 (F_trap)",
                      "family": h.FAMILY, "eps": h.EPSS, "eps0": h.EPS0, "F1D_abs": h.F1D_abs(),
                      "lines": {"per_k": "|c₁ε₀| ≤ 0.05 |c₂ε₀²| (k ごと)",
                                "control": "|c₁(M)| ≤ 0.05 |c₁(対照)|、対照の |c₁ε₀| が床 (2 会計差) の 10 倍以上のものだけ",
                                "k0down": "c₁ε₀ が k₀ 64 の値から 1e-8 以上動くか (記録)"},
                      "note": "判定は親と上位。h:ab_L100_resampled_N1281 は (h) の値 (rao_t4h.json) を再掲。ref_M_* は登録外の参照"},
           **h.git_info(), "per_k": per_k, "controls": ctrl, "k0down": k0down, "rows": rows, "parts": parts, "total_s": total_s}
    (HERE / "rao_t4h_follow.json").write_text(json.dumps(h._jsonable(out), ensure_ascii=False, indent=1))
    h.write_csv(HERE / "rao_t4h_follow.csv", rows)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--assemble-only", action="store_true")
    a = ap.parse_args()
    if a.workers > 2:
        raise SystemExit("--workers は 2 まで (手元 PC はユーザの作業機)")
    t0 = time.time()
    PARTS.mkdir(exist_ok=True)
    if not a.assemble_only:
        todo = [j for j in jobs() if (a.only is None or any(s in j["name"] for s in a.only.split(",")))
                and not (PARTS / f"{j['name']}.json").exists()]
        print(f"回すジョブ {[j['name'] for j in todo]} (並列 {a.workers})", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = [ex.submit(run, j) for j in todo]
            for fu in as_completed(futs):
                fu.result()
    out = assemble(time.time() - t0)
    print(json.dumps(h._jsonable({k: out[k] for k in ("per_k", "controls", "k0down")}), ensure_ascii=False, indent=1), flush=True)
    print(f"rao_t4h_follow.csv {len(out['rows'])} 行 ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
