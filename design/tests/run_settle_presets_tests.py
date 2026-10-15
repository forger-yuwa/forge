#!/usr/bin/env python3
r"""収束の判定のプリセット `settle_presets` と節点の量の抽出の一致の試験 (plan tooling-settle-convergence-and-monitor §5.1 #2〜#3、§4.3)。

    design/.venv-opt/bin/python design/tests/run_settle_presets_tests.py

(a) プリセット `wt_nozzle_v1` の解決: case/71 の問題 YAML (`problem_m55_p3_b0.yaml`) と手元の run
    `case/71.m55_combustion_wt/run_0010_ns_b0/` の prepare_info・solverConfig・bcondConfig で、量 6 つ・M_d・格子・間隔・
    入口と出口・残差の列・線の既定 (η 0.03、x (x_E + 2, x_F − 1) r_t)・壁 (等温)・前提の記録・往復 (書いて読める) を確かめる。
    断熱壁の bcond では wall_heat_nodes を外して excluded に理由を書く。足し引き・閾値の上書き・line_eta・x_range_m。
(b) **一致**: 手元の run の全場 `res_86000.h5` 1 枚から `probe_node_ids` の節点だけを抜いた合成の `res_nodes_86000.h5`
    (forge の monitorOutputs と同じ: COORD は全場の /MESH/COORD の写し、VALUE は全場と同じ型) で `settle` が作る
    `exit_core_mach`・`line_mach` (xq は v5c_series と同じ x 4.0〜7.95 m・2001 点; η 0.03 と 0.1)・`exit_flow_angle` が、全場から
    `nozzle_report.load_field` + `eta_line` (+ `stop_rule.recorded_quantities`) で作った値と**ビット同一**。節点の表は forge の書式で
    書いて読み戻せる。η 0.03 は格子の第 1 層 (j = 1) の η と丸めの差で一致・近接するので、列ごとに np.interp が使う対を数えて報告する。
(c) 例外: evaluate.settle の知らないキー・知らないプリセット・閾値の知らない量とキー・quantities_remove の知らない量・
    farfield / 知らない kind / 出口の無い bcond・凝縮・W が間隔の倍数でない・prepare_info の Md の欠け。
"""
from __future__ import annotations

import copy
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
from forge_design.evaluate import settle_presets as SP  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.report.nozzle_report import eta_line, load_field  # noqa: E402
import settle  # noqa: E402

CASE = ROOT / "case/71.m55_combustion_wt"
PROBLEM = CASE / "problem_m55_p3_b0.yaml"
RUN_CANDIDATES = [CASE / "_runs/run_0010_ns_b0", CASE / "run_0010_ns_b0"]
RES = "res_86000.h5"
STEP = 86000
FAIL = []
NCHK = 0
TMP = Path(tempfile.mkdtemp(prefix="settle_presets_"))


def chk(name, cond, detail=""):
    global NCHK
    NCHK += 1
    if not cond:
        FAIL.append(name)
    print(f"  [{'OK ' if cond else 'NG '}] {name}" + (f"  — {detail}" if detail else ""))


def find_run():
    """prepare_info・solverConfig・bcondConfig・nozzle.h5・res_86000.h5 が揃う最寄りの手元の run。"""
    need = ("prepare_info.json", "solverConfig.yaml", "bcondConfig.yaml", "nozzle.h5", RES)
    for r in RUN_CANDIDATES:
        if all((r / f).exists() for f in need):
            return r
    for r in sorted(CASE.glob("**/run_*")):
        if all((r / f).exists() for f in need[:4]) and list(r.glob("res_[0-9]*.h5")):
            return r
    return None


def load_inputs(run):
    info = json.loads((run / "prepare_info.json").read_text())
    cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
    return info, cfg, bc


def raises(fn, word=None):
    try:
        fn()
    except ValueError as e:
        return word is None or word in str(e), str(e)
    return False, "例外が出ない"


# ------------------------------------------------------------------ (a)
def test_resolve(run):
    print("=== (a) プリセット wt_nozzle_v1 の解決 ===")
    prob = yaml.safe_load(PROBLEM.read_text())
    info, cfg, bc = load_inputs(run)
    es = prob.get("evaluate", {}).get("settle") or {"preset": "wt_nozzle_v1"}
    print(f"  問題 YAML {PROBLEM.relative_to(ROOT)} の evaluate.settle: {prob.get('evaluate', {}).get('settle')!r} → {es}")
    crit = SP.resolve(es, info, cfg, bc)
    kinds = [q["kind"] for q in crit["quantities"]]
    q = {x["name"]: x for x in crit["quantities"]}
    chk("a1 量 6 つ (等温壁なので wall_heat_nodes を含む)", kinds == ["residual", "mass_balance", "exit_core_mach", "line_mach",
                                                              "exit_flow_angle", "wall_heat_nodes"], str(kinds))
    chk("a2 M_d・格子・間隔・窓 (prepare_info と solverConfig から)",
        q["line_mach"]["Md"] == 5.5 == float(prob["spec"]["M_design"]) and crit["structured"] == {"ni": info["mesh"]["ni"], "nj": info["mesh"]["nj"],
                                                                                                 "scale_m": info["scale_m"], "mesh_file": "nozzle.h5"}
        and crit["grid"]["out_step_interval"] == cfg["time"]["outStepInterval"] == 2000 and crit["window"] == {"W": 20000, "min_points": 11}
        and crit["hold"] == 40000, f"structured {crit['structured']} grid {crit['grid']}")
    chk("a3 入口・出口は bcond の kind から (inlet_Pressure 1 / outlet_statPress 2)",
        q["mass_balance"]["inlet"] == [1] and q["mass_balance"]["outlet"] == [2], f"{q['mass_balance']['inlet']} {q['mass_balance']['outlet']}")
    chk("a4 残差の列 (SST・化学種 2 つ)", q["residual"]["required"] == ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega",
                                                                  "rms_roY0", "rms_roY1"] and q["residual"]["inactive_ok"] == ["rms_roUz"],
        str(q["residual"]["required"]))
    S = info["scale_m"]
    X, R = SP.mesh_xr(run / "nozzle.h5", info["mesh"]["ni"], info["mesh"]["nj"], S)
    xr = q["line_mach"]["x_range_m"]
    chk("a5 線の既定 η 0.03・x (x_E + 2, x_F − 1) r_t (x_F = physical_wall.domain[1] = 格子の最後の列の軸の x)",
        q["line_mach"]["eta"] == 0.03 and xr == [(info["x_E"] + 2) * S, (info["physical_wall"]["domain"][1] - 1) * S]
        and float(X[-1, 0]) == info["physical_wall"]["domain"][1] and q["line_mach"]["n"] == 2001,
        f"x {xr[0]:.6f}〜{xr[1]:.6f} m, X[-1,0] {float(X[-1, 0])!r}")
    chk("a6 出口コア M の tol = 0.1 × 0.05 % × M_d、閾値はすべて暫定",
        abs(q["exit_core_mach"]["thresholds"]["tol"] - 0.1 * 0.0005 * 5.5) < 1e-18 and all(x["thresholds"]["provisional"] for x in crit["quantities"]),
        f"tol {q['exit_core_mach']['thresholds']['tol']!r}")
    chk("a7 壁: wall (physID 3)、軸対称", q["wall_heat_nodes"]["walls"] == [{"name": "wall", "physID": 3}] and q["wall_heat_nodes"]["axisymmetric"])
    pre = {(p["quantity"], p["requires"]): p["satisfied"] for p in crit["preconditions"]}
    unmet = [k[1] for k, v in pre.items() if not v]
    chk("a8 前提の記録: run_0010 は boundaryMassflux・nodeProbeFile・interfaceDiag が無い (満たされていない)、壁の outputHDFflg 1 は満たす",
        len(unmet) == 3 and any("interfaceDiag" in u for u in unmet) and any("outputHDFflg" in k[1] and v for k, v in pre.items()), str(unmet))
    rd = TMP / "a_run"
    rd.mkdir()
    for fn in ("prepare_info.json", "solverConfig.yaml", "bcondConfig.yaml"):
        shutil.copy(run / fn, rd / fn)
    w = SP.write_settle_criteria(rd, es, info, cfg, bc)
    lc = settle.load_criteria(rd)
    chk("a9 書いた settle_criteria.yaml を settle が読める (status ok、sha256 一致、出典にプリセットの版・sha256・入力の sha256)",
        lc["status"] == "ok" and lc["sha256"] == w["sha256"] and lc["criteria"]["source"]["preset"] == "wt_nozzle_v1"
        and len(lc["criteria"]["source"]["preset_sha256"]) == 64 and set(lc["criteria"]["source"]["inputs_sha256"]) == {
            "prepare_info.json", "solverConfig.yaml", "bcondConfig.yaml"}, lc["status"])
    chk("a9 往復で閾値・座標の範囲・scale_m が変わらない", lc["criteria"]["quantities"] == settle.validate_criteria(crit)["quantities"]
        and lc["criteria"]["structured"]["scale_m"] == info["scale_m"])

    bc2 = copy.deepcopy(bc)
    bc2["wall"]["kind"] = "wall"
    c2 = SP.resolve(es, info, cfg, bc2)
    chk("a10 断熱壁 (kind wall) → wall_heat_nodes を外し excluded に理由", "wall_heat_nodes" not in [x["kind"] for x in c2["quantities"]]
        and c2["excluded"] and c2["excluded"][0]["kind"] == "wall_heat_nodes" and "等温壁が無い" in c2["excluded"][0]["reason"],
        c2["excluded"][0]["reason"][:40] if c2["excluded"] else "")
    bc3 = copy.deepcopy(bc2)
    bc3["wall"]["ints"] = {"wallProfile": 1}
    c3 = SP.resolve(es, info, cfg, bc3)
    chk("a11 ints.wallProfile の壁は等温壁として扱う", "wall_heat_nodes" in [x["kind"] for x in c3["quantities"]])
    es4 = {"preset": "wt_nozzle_v1", "thresholds": {"line_mach": {"tol_pos": 0.02}, "residual": {"provisional": False}},
           "quantities_remove": ["exit_flow_angle"], "line_eta": 0.1, "x_range_m": [4.0, 7.95],
           "quantities_add": [{"kind": "line_mach", "name": "line_axis", "eta": 0.0, "x_range_m": [4.0, 7.95], "n": 2001, "Md": 5.5,
                               "thresholds": {"tol_maxabs": 0.01, "tol_halfamp": 0.005, "tol_pos": 0.01, "provisional": True}}]}
    c4 = SP.resolve(es4, info, cfg, bc)
    q4 = {x["name"]: x for x in c4["quantities"]}
    chk("a12 足し引き・閾値の上書き・line_eta・x_range_m", "exit_flow_angle" not in q4 and "line_axis" in q4
        and q4["line_mach"]["thresholds"]["tol_pos"] == 0.02 and q4["line_mach"]["eta"] == 0.1 and q4["line_mach"]["x_range_m"] == [4.0, 7.95]
        and q4["residual"]["thresholds"]["provisional"] is False
        and any(e["name"] == "exit_flow_angle" and "quantities_remove" in e["reason"] for e in c4["excluded"]))
    return crit


# ------------------------------------------------------------------ (b)
def write_res_nodes(run_out: Path, src_res: Path, ids, step: int, coord_dtype=None):
    """forge の monitorOutputs::writeNodeProbes と同じ形の合成: COORD は全場の /MESH/COORD の写し (coord_dtype を与えると
    その型に丸める = float のビルドの geom_float)、VALUE は同じ型のまま。"""
    names = ("ro", "Ux", "Uy", "Uz", "P", "T", "sonic")
    with h5py.File(src_res, "r") as f:
        coord = f["/MESH/COORD"][:].reshape(-1, 3)
        V = {k: f["/VALUE/" + k][:] for k in names}
    if coord_dtype is not None:
        coord = coord.astype(coord_dtype)
    with h5py.File(run_out / f"res_nodes_{step}.h5", "w") as g:
        g.attrs["step"] = np.int32(step)
        g.attrs["nodeProbeFile"] = "node_probes.txt"
        g["/MESH/COORD"] = coord[ids].ravel()
        g["/MESH/NODE_ID"] = np.asarray(ids, dtype=np.int64)
        for k in names:
            g["/VALUE/" + k] = V[k][ids]
    return {k: str(V[k].dtype) for k in names}, str(coord.dtype)


def test_bit_identity(run, crit):
    print("=== (b) 節点の表から作る量と全場からの抽出のビット一致 ===")
    info, cfg, bc = load_inputs(run)
    ni, nj, S = info["mesh"]["ni"], info["mesh"]["nj"], float(info["scale_m"])
    with h5py.File(run / "nozzle.h5") as f:
        cm = f["/MESH/COORD"][:]
    with h5py.File(run / RES) as f:
        cr = f["/MESH/COORD"][:]
    chk("b0 全場の /MESH/COORD (forge が書く節点の座標) = 格子の h5 の COORD (load_field が読む)", np.array_equal(cm, cr) and cm.dtype == cr.dtype,
        f"dtype {cm.dtype}")
    c = copy.deepcopy(crit)
    lines = [{"kind": "line_mach", "name": f"line_eta{eta}", "eta": eta, "x_range_m": [4.0, 7.95], "n": 2001, "Md": 5.5,
              "thresholds": {"tol_maxabs": 0.01, "tol_halfamp": 0.005, "tol_pos": 0.01}} for eta in (0.03, 0.1)]
    c["quantities"] = [q for q in c["quantities"] if q["kind"] != "line_mach"] + lines
    settle.validate_criteria(c)
    X, R = SP.mesh_xr(run / "nozzle.h5", ni, nj, S)
    t0 = time.time()
    ids = SP.probe_node_ids(ni, nj, X, R, c)
    rd = TMP / "b_run"
    rd.mkdir()
    pf = SP.write_node_probe_file(rd / "node_probes.txt", ids, header="settle_presets.probe_node_ids の試験")
    back = SP.read_node_probe_file(pf)
    chk("b1 節点の表を forge の書式で書いて読み戻せる (1 行 1 番号、# はコメント、昇順・重複なし・範囲内)",
        back == ids.tolist() and ids.min() >= 0 and ids.max() < ni * nj, f"{len(ids)} 節点 (全 {ni * nj} の {100 * len(ids) / (ni * nj):.2f} %)")
    vdt, cdt = write_res_nodes(rd, run / RES, ids, STEP)
    obs = settle._read_nodes(rd / f"res_nodes_{STEP}.h5", STEP, c["structured"])
    print(f"  節点の表 {len(ids)} 節点、合成 res_nodes の COORD {cdt}・VALUE {vdt['Ux']} ({time.time() - t0:.1f} s)")

    F = load_field(run, RES)
    t1 = time.time()
    # 出口コア M: stop_rule.recorded_quantities (座標は load_field と同じ COORD / S、VALUE は float64)
    geo = SR.RunGeometry(run)
    sh = (ni, nj)
    with h5py.File(run / RES) as f:
        Vr = {k: f["/VALUE/" + k][:].reshape(sh) for k in ("ro", "Ux", "Uy", "sonic")}
    ref_core = SR.recorded_quantities(geo.Xf, geo.Rf, Vr, geo.x_E)["exit_core_M"]
    eta_last = F["R"][-1] / F["R"][-1, -1]
    core = (eta_last >= 0.05) & (eta_last <= 0.7)
    ref_core_nr = float(F["V"]["M"][-1][core].mean())
    got_core = settle.exit_core_mach_of(obs, ni, nj, [0.05, 0.7])
    chk("b2 exit_core_mach: settle (節点の表) = stop_rule.recorded_quantities の exit_core_M = nozzle_report の M_mean (ビット同一)",
        got_core == ref_core == ref_core_nr, f"{got_core!r} / {ref_core!r} / {ref_core_nr!r}")
    ref_ang = np.degrees(np.arctan2(F["V"]["Uy"][-1][core], F["V"]["Ux"][-1][core]))
    got_ang = settle.exit_flow_angle_of(obs, ni, nj, [0.05, 0.7])
    chk("b3 exit_flow_angle: 節点ごとにビット同一", np.array_equal(got_ang, ref_ang) and got_ang.dtype == ref_ang.dtype,
        f"{len(got_ang)} 節点、最大 |θ| {np.abs(ref_ang).max():.4f}°")
    xq = np.linspace(4.0 / S, 7.95 / S, 2001)     # v5c_series.py と同じ
    for eta in (0.03, 0.1):
        ref_M = eta_line(F, "M", eta, xq)
        ref = {"maxabs": float(np.max(np.abs(100.0 * (ref_M / 5.5 - 1.0)))), "halfamp": float(100.0 * (ref_M.max() - ref_M.min()) / (2.0 * 5.5))}
        got_M = settle.eta_line_nodes(obs, "M", eta, xq, ni, nj)
        got = settle.line_quantities(got_M, 5.5)
        ref_d = 100.0 * (ref_M / 5.5 - 1.0)
        chk(f"b4 line_mach η {eta}: 線の M・位置ごとの d・maxabs・halfamp がビット同一 (eta_line + v5c_series の式)",
            np.array_equal(got_M, ref_M) and np.array_equal(got["d"], ref_d) and got["maxabs"] == ref["maxabs"] and got["halfamp"] == ref["halfamp"],
            f"maxabs {got['maxabs']:.6f} %、halfamp {got['halfamp']:.6f} %")
        # η が節点の η と丸めの差で一致・近接する列で、np.interp が使う対
        xa = F["X"][:, 0]
        i = np.clip(np.searchsorted(xa, xq) - 1, 0, ni - 2)
        cols = np.unique(np.concatenate([i, i + 1]))
        pairs = {}
        near = 0
        for col in cols:
            xp = F["R"][col] / F["R"][col, -1]
            k = int(np.searchsorted(xp, eta, side="right")) - 1
            pairs[k] = pairs.get(k, 0) + 1
            near += int(np.min(np.abs(xp - eta)) < 1e-9)
        print(f"    η {eta}: 線が使う列 {len(cols)}、np.interp の対 (j, j+1) の j ごとの列の数 {dict(sorted(pairs.items()))}、"
              f"節点の η と |差| < 1e-9 の列 {near}")
    # 全節点の res_nodes でも同じ (部分集合と全節点で計算が変わらない)
    rd2 = TMP / "b_full"
    rd2.mkdir()
    write_res_nodes(rd2, run / RES, np.arange(ni * nj), STEP)
    obs_full = settle._read_nodes(rd2 / f"res_nodes_{STEP}.h5", STEP, c["structured"])
    same = all(np.array_equal(settle.eta_line_nodes(obs_full, "M", eta, xq, ni, nj), settle.eta_line_nodes(obs, "M", eta, xq, ni, nj))
               for eta in (0.03, 0.1))
    chk("b5 全節点の res_nodes と節点の表の res_nodes で線の値が同一", same)
    # 節点の表から η を挟む節点を 1 つ抜くと判定不能 (黙って別の値を返さない)
    col = int(np.clip(np.searchsorted(F["X"][:, 0], xq[1000]) - 1, 0, ni - 2))
    xp = F["R"][col] / F["R"][col, -1]
    k = int(np.searchsorted(xp, 0.1, side="right")) - 1
    drop = col * nj + k + 1
    ids2 = ids[ids != drop]
    rd3 = TMP / "b_drop"
    rd3.mkdir()
    write_res_nodes(rd3, run / RES, ids2, STEP)
    obs3 = settle._read_nodes(rd3 / f"res_nodes_{STEP}.h5", STEP, c["structured"])
    try:
        settle.eta_line_nodes(obs3, "M", 0.1, xq, ni, nj)
        ok = False
        why = "例外が出ない"
    except settle._DataMissing as e:
        ok = True
        why = str(e)
    chk("b6 η を挟む節点が節点の表に無ければデータ不足 (別の対で補間しない)", ok, why[:80])
    # float のビルドを模す: res_nodes の COORD を float32 (geom_float) で書く。structured.mesh_file (格子の h5) があれば
    # load_field と同じ double の COORD を使うので一致し、無ければ (res_nodes の float32 の COORD) 一致しない
    rd4 = TMP / "b_f32"
    rd4.mkdir()
    write_res_nodes(rd4, run / RES, ids, STEP, coord_dtype=np.float32)
    mc = h5py.File(run / "nozzle.h5")["/MESH/COORD"][:].reshape(-1, 3)
    obs_m = settle._read_nodes(rd4 / f"res_nodes_{STEP}.h5", STEP, c["structured"], mc)
    obs_r = settle._read_nodes(rd4 / f"res_nodes_{STEP}.h5", STEP, c["structured"], None)
    ok_m = all(np.array_equal(settle.eta_line_nodes(obs_m, "M", eta, xq, ni, nj), eta_line(F, "M", eta, xq)) for eta in (0.03, 0.1)) \
        and settle.exit_core_mach_of(obs_m, ni, nj, [0.05, 0.7]) == ref_core
    try:
        diff_r = max(float(np.abs(settle.eta_line_nodes(obs_r, "M", eta, xq, ni, nj) - eta_line(F, "M", eta, xq)).max()) for eta in (0.03, 0.1))
    except settle._DataMissing as e:
        diff_r = f"判定不能 ({e})"
    chk("b7 float のビルドを模した res_nodes (COORD float32): mesh_file の COORD を使えばビット同一 (座標の源 mesh_file)",
        ok_m and obs_m["coord_source"] == "mesh_file", f"mesh_file 無しでは線の M の差 {diff_r}")
    bad = mc.copy()
    bad[ids[0], 0] += 1e-3
    try:
        settle._read_nodes(rd4 / f"res_nodes_{STEP}.h5", STEP, c["structured"], bad)
        ok = False
    except settle._DataMissing as e:
        ok = "合わない" in str(e)
    chk("b8 res_nodes の COORD と mesh_file の COORD が合わなければデータ不足 (別の格子)", ok)
    print(f"  (全場の抽出 {time.time() - t1:.1f} s)")


# ------------------------------------------------------------------ (c)
def test_errors(run):
    print("=== (c) 例外 ===")
    info, cfg, bc = load_inputs(run)
    es = {"preset": "wt_nozzle_v1"}
    cases = [
        ("c1 evaluate.settle の知らないキー", lambda: SP.resolve({**es, "windows": 1}, info, cfg, bc), "知らないキー"),
        ("c2 知らないプリセット", lambda: SP.resolve({"preset": "wt_nozzle_v9"}, info, cfg, bc), "知らない"),
        ("c3 プリセットが無い", lambda: SP.resolve({}, info, cfg, bc), "知らない"),
        ("c4 閾値の上書きの知らない量", lambda: SP.resolve({**es, "thresholds": {"exit_mach": {"tol": 1}}}, info, cfg, bc), "量 'exit_mach'"),
        ("c5 閾値の上書きの知らないキー", lambda: SP.resolve({**es, "thresholds": {"line_mach": {"tol": 1}}}, info, cfg, bc), "知らないキー"),
        ("c6 quantities_remove の知らない量", lambda: SP.resolve({**es, "quantities_remove": ["wall_heat"]}, info, cfg, bc), "プリセットに無い"),
        ("c7 quantities_add の知らない種類 (settle の検査)",
         lambda: SP.resolve({**es, "quantities_add": [{"kind": "pitot", "thresholds": {}}]}, info, cfg, bc), "知らない種類"),
    ]
    bcf = copy.deepcopy(bc)
    bcf["inlet"]["kind"] = "farfield"
    cases.append(("c8 farfield (入口か出口か判別できない)", lambda: SP.resolve(es, info, cfg, bcf), "判別できない"))
    bck = copy.deepcopy(bc)
    bck["axis"]["kind"] = "symmetry_plane"
    cases.append(("c9 知らない kind", lambda: SP.resolve(es, info, cfg, bck), "判別できない"))
    bco = copy.deepcopy(bc)
    del bco["outlet"]
    cases.append(("c10 出口が無い", lambda: SP.resolve(es, info, cfg, bco), "出口"))
    cfc = copy.deepcopy(cfg)
    cfc["condensation"] = {"condensation": 1, "nCondSpecies": 1}
    cases.append(("c11 凝縮 (残差の列を決められない)", lambda: SP.resolve(es, info, cfc, bc), "凝縮"))
    cfi = copy.deepcopy(cfg)
    cfi["time"]["outStepInterval"] = 3000
    cases.append(("c12 W 20,000 が間隔 3,000 の倍数でない", lambda: SP.resolve(es, info, cfi, bc), "倍数でない"))
    inf = copy.deepcopy(info)
    del inf["Md"]
    cases.append(("c13 prepare_info に Md が無い", lambda: SP.resolve(es, inf, cfg, bc), "Md"))
    cfu = copy.deepcopy(cfg)
    cfu["time"]["unsteady"] = 1
    cases.append(("c14 非定常の run", lambda: SP.resolve(es, info, cfu, bc), "定常"))
    for name, fn, word in cases:
        ok, msg = raises(fn, word)
        chk(name + " → ValueError", ok, msg[:90])


def main():
    print(f"numpy {np.__version__}, python {sys.version.split()[0]}")
    run = find_run()
    if run is None:
        print("SKIP: 手元に prepare_info・solverConfig・bcondConfig・nozzle.h5・res_*.h5 の揃う case/71 の run が無い")
        print("VERDICT: SKIPPED")
        return 2
    print(f"使う run: {run.relative_to(ROOT)} (候補 {[str(r.relative_to(ROOT)) for r in RUN_CANDIDATES]}; "
          f"_runs/ の写しには h5 が無いので case 直下の run を使う)")
    crit = test_resolve(run)
    test_bit_identity(run, crit)
    test_errors(run)
    shutil.rmtree(TMP, ignore_errors=True)
    print(f"\n{NCHK - len(FAIL)}/{NCHK} OK")
    if FAIL:
        print("FAIL:", *FAIL, sep="\n  ")
    print(f"VERDICT: {'PASS' if not FAIL else 'FAIL'}")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
