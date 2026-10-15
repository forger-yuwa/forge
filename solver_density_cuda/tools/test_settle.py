#!/usr/bin/env python3
r"""`settle.py` / `check_settled.py` の試験 (plan tooling-settle-convergence-and-monitor §6 V2、§4.8 の 1・4)。合成データだけを使う。

    python3 solver_density_cuda/tools/test_settle.py            # 全件 (design/.venv-opt の python でも system python3 でも)

V2 の負例はすべて SETTLED にしないこと、各負例で理由の分類 (閾値の超過 threshold / データ不足 data / 未設定 unset /
非有限 nonfinite) が期待どおりであることを確かめる。正例: 2,000 ごと 0〜40,000 の一定系列は SETTLED、一定系列を 80,001 で止めた
(非格子の終端) 場合も SETTLED。

合成の run: 構造格子 ni 30 × nj 10 (x = i − 2 [r_t]、r_w = 1 + 0.05 i、r = r_w·j/(nj − 1)、scale 0.1 m)。節点の出力は全節点を書く
(判定は部分集合でも全節点でも同じ計算; 部分集合の一致は design/tests/run_settle_presets_tests.py が実データで確かめる)。
**格子点 0 の観測**: forge は本段の開始 (step 0) に残差の行・節点・流量の出力を書かない。V2 の正例「0〜40,000」は全格子点
(0 を含む 21 点) に観測があるという意味なので、合成では格子点 0 の観測 (残差の行 −1・res_nodes_0・流量の行 0・壁の h5) を置く
(`start_obs=True`)。forge と同じ形 (`start_obs=False`) では 40,000 で判定不能・42,000 で SETTLED になることを別に確かめる。
"""
from __future__ import annotations

import io
import json
import math
import shutil
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import settle  # noqa: E402
import check_settled  # noqa: E402

NI, NJ, S = 30, 10, 0.1
MD = 3.0
DT = 2000
RES_COLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"]
RES_BASE = {"rms_ro": 1e-6, "rms_roUx": 1e-4, "rms_roUy": 5e-5, "rms_roUz": 0.0, "rms_roe": 1e-1, "rms_roK": 1e-4, "rms_roOmega": 1e1}
DQ_COLS = ["rms_dq_ro", "rms_dq_roUx", "rms_dq_roUy", "rms_dq_roUz", "rms_dq_roe"]
TMP = Path(tempfile.mkdtemp(prefix="test_settle_"))
FAILS = []
NCHK = 0


def chk(name, cond, detail=""):
    global NCHK
    NCHK += 1
    if not cond:
        FAILS.append(name)
    print(f"  [{'OK ' if cond else 'NG '}] {name}" + (f"  — {detail}" if detail else ""))


# ------------------------------------------------------------------ 合成の run
def grid_xr():
    i = np.arange(NI)[:, None]
    j = np.arange(NJ)[None, :]
    x = (i - 2.0) * np.ones((1, NJ))
    rw = 1.0 + 0.05 * i
    r = rw * j / (NJ - 1)
    return x * S, r * S       # [m] (ni, nj)


def bump(i, c, A):
    """列 c〜c+4 の離散の山 (0, 0.5, 1, 0.5, 0)·A [%]。整数の列だけ動かすので maxabs・halfamp は位置に依らない。"""
    pat = {0: 0.0, 1: 0.5, 2: 1.0, 3: 0.5, 4: 0.0}
    return np.array([pat.get(int(k - c), 0.0) for k in i]) * A


def default_M(s):
    return np.full(NI, MD)


def write_nodes(run, s, M_cols, *, nan=False, ids=None):
    xm, rm = grid_xr()
    a = 300.0
    eta = rm / rm[:, -1:]
    th = np.radians(0.5) * eta
    M = M_cols[:, None] * np.ones((1, NJ))
    ux = (M * a * np.cos(th)).ravel()
    uy = (M * a * np.sin(th)).ravel()
    if nan:
        ux = ux.copy()
        ux[5] = np.nan
    coord = np.stack([xm.ravel(), rm.ravel(), np.zeros(NI * NJ)], axis=1)
    ids = np.arange(NI * NJ) if ids is None else np.asarray(ids)
    with h5py.File(run / f"res_nodes_{s}.h5", "w") as f:
        f.attrs["step"] = np.int32(s)
        f.attrs["nodeProbeFile"] = "node_probes.txt"
        f["/MESH/COORD"] = coord[ids].ravel()
        f["/MESH/NODE_ID"] = ids.astype(np.int64)
        for k, v in (("ro", np.ones(NI * NJ)), ("Ux", ux), ("Uy", uy), ("Uz", np.zeros(NI * NJ)), ("P", np.ones(NI * NJ)),
                     ("T", np.ones(NI * NJ)), ("sonic", np.full(NI * NJ, a))):
            f["/VALUE/" + k] = v[ids]


def write_wall(run, s, Q):
    xm, rm = grid_xr()
    coord = np.stack([xm[:, -1], rm[:, -1], np.zeros(NI)], axis=1)
    conne = []
    for i in range(NI - 1):
        conne += [2, 2, i, i + 1]
    with h5py.File(run / f"res_wall_3_{s}.h5", "w") as f:
        f.attrs["step"] = np.int32(s)
        f["/MESH/COORD"] = coord.ravel()
        f["/MESH/CONNE"] = np.asarray(conne, dtype=np.int32)
        f["/VALUE/iface_Qf_eff"] = np.asarray(Q, dtype=float)
        f["/VALUE/qwall"] = np.zeros(NI)


def write_massflux(run, rows):
    with open(run / "boundary_massflux.csv", "w") as f:
        f.write("# boundary_massflux.csv — 合成\n# 符号: 流出が正\n")
        f.write("step,mdot_inlet_1,raw_inlet_1,mdot_outlet_2,raw_outlet_2,mdot_wall_3,raw_wall_3,mdot_axis_4,raw_axis_4\n")
        for s, mi, mo in rows:
            f.write(f"{s},{-mi!r},{-mi / (2 * math.pi)!r},{mo!r},{mo / (2 * math.pi)!r},0,0,0,0\n")


def write_residual(run, steps, fn, cols=RES_COLS, name="residual_history.csv"):
    with open(run / name, "w") as f:
        f.write("step,inner,phase," + ",".join(cols + DQ_COLS) + "\n")
        for s in steps:
            v = fn(s)
            f.write(f"{s},-1,outer_begin," + ",".join(repr(v.get(c, 0.0)) for c in cols) + ",0,0,0,0,0\n")
            f.write(f"{s},-1,outer_end," + ",".join(repr(v.get(c, 0.0)) for c in cols) + ",0,0,0,0,0\n")


def base_criteria(**over):
    c = {
        "schema": settle.SCHEMA,
        "source": {"made_by": "test_settle.py"},
        "structured": {"ni": NI, "nj": NJ, "scale_m": S},
        "grid": {"out_step_interval": DT, "start_step": 0, "segment": "main"},
        "window": {"W": 20000, "min_points": 11},
        "hold": 40000,
        "quantities": [
            {"kind": "residual", "required": ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega"],
             "inactive_ok": ["rms_roUz"],
             "thresholds": {"slope_dec": math.log10(1.05), "amp_dec": 0.1, "rise_dec": math.log10(1.10), "excursion_dec": 0.5,
                            "provisional": True}},
            {"kind": "mass_balance", "inlet": [1], "outlet": [2],
             "thresholds": {"mean_abs": 5e-4, "swing": 5e-3, "mdot_change_rel": 1e-4, "provisional": True}},
            {"kind": "exit_core_mach", "eta_range": [0.05, 0.7], "thresholds": {"tol": 0.1 * 0.0005 * MD, "provisional": True}},
            {"kind": "line_mach", "eta": 0.03, "x_range_m": [3 * S, 20 * S], "n": 2001, "Md": MD,
             "thresholds": {"tol_maxabs": 0.01, "tol_halfamp": 0.005, "tol_pos": 0.01, "provisional": True}},
            {"kind": "exit_flow_angle", "eta_range": [0.05, 0.7], "thresholds": {"tol": 0.01, "provisional": True}},
            {"kind": "wall_heat_nodes", "walls": [{"name": "wall", "physID": 3}], "axisymmetric": True,
             "thresholds": {"q_tol": 1000.0, "eps_rel": 1e-3, "provisional": True}},
        ],
    }
    c.update(over)
    return c


def build(name, *, end=40000, stop=None, start_obs=True, resid_fn=None, resid_cols=RES_COLS, M_fn=default_M,
          Q_fn=None, mdot_fn=None, skip=(), skip_resid_rows=(), criteria=None, no_massflux=False, nan_nodes_at=None,
          extra_resid_rows=(), stopped_at=True):
    """合成の run を作る。skip: 観測を置かない格子点 (残差の行 s−1・節点・流量・壁)。"""
    run = TMP / name
    if run.exists():
        shutil.rmtree(run)
    run.mkdir()
    grid = [s for s in range(0, end + 1, DT) if (start_obs or s > 0) and s not in skip]
    obs = grid + ([stop] if stop is not None else [])
    last_row = (stop if stop is not None else end) - 1
    resid_fn = resid_fn or (lambda s: dict(RES_BASE))
    rows = set(range(0, last_row + 1, 10)) | {s - 1 for s in obs}
    rows |= set(extra_resid_rows)
    if start_obs:
        rows.add(-1)
    rows = sorted(r for r in rows if r <= last_row and (r + 1) not in skip and r not in skip_resid_rows)
    write_residual(run, rows, resid_fn, cols=resid_cols)
    mdot_fn = mdot_fn or (lambda s: (10.0, 10.0))
    if not no_massflux:
        write_massflux(run, [(s, *mdot_fn(s)) for s in obs])
    Q_fn = Q_fn or (lambda s: 100.0 + np.arange(NI, dtype=float))
    for s in obs:
        write_nodes(run, s, M_fn(s), nan=(nan_nodes_at is not None and s == nan_nodes_at))
        write_wall(run, s, Q_fn(s))
    if stop is not None and stopped_at:
        (run / "STOPPED_AT").write_text(json.dumps({"status": "ok", "step": stop, "residual_step": stop - 1, "files": []}))
    if criteria is not False:
        (run / settle.CRITERIA_FILE).write_text(yaml.safe_dump(criteria or base_criteria(), allow_unicode=True, sort_keys=False))
    return run


def run_cli(run, *args):
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            rc = check_settled.main([str(run), "--no-write", *args])
        except SystemExit as e:
            rc = e.code
    return rc, buf.getvalue()


def qres(res, name):
    for q in res["quantities"]:
        if q["name"] == name:
            return q
    raise KeyError(name)


def cats(q):
    return set(q.get("categories", []))


def failed_checks(q):
    return [c["check"] for c in q.get("checks", []) if not c["ok"]]


# ------------------------------------------------------------------ 試験
def test_positive():
    print("=== 正例 ===")
    run = build("p1_const_0_40000")
    res = settle.judge(run)
    chk("P1 2,000 ごと 0〜40,000 の一定系列は SETTLED", res["verdict"] == "SETTLED",
        f"{res['verdict']} {[(q['name'], q['verdict'], q['reasons'][:1]) for q in res['quantities'] if q['verdict'] != 'SETTLED']}")
    chk("P1 窓 C は端点込み 11 点・直前の窓 P も 11 点", len(res["window_steps"]["C"]) == 11 and len(res["window_steps"]["P"]) == 11
        and res["window_steps"]["P"][0] == 0, str(res["window_steps"]))
    rc, out = run_cli(run)
    chk("P1 CLI の終了コード 0 と VERDICT の行の暫定の注記", rc == 0 and settle.PROVISIONAL_TEXT in out.splitlines()[-1], out.splitlines()[-1])

    run = build("p2_const_stop_80001", end=80000, stop=80001)
    res = settle.judge(run)
    chk("P2 一定系列を 80,001 で止めた (非格子の終端) は SETTLED", res["verdict"] == "SETTLED"
        and res["terminal_step"] == 80000 and res["stop_step"] == 80001 and res["end_source"] == "STOPPED_AT",
        f"{res['verdict']} T={res.get('terminal_step')} stop={res.get('stop_step')} src={res.get('end_source')}")
    chk("P2 停止の観測は窓に入れず別の検査に使う (各量に「停止の観測」の検査がある)",
        all(any("停止の観測" in c["check"] for c in q["checks"]) for q in res["quantities"]),
        str([q["name"] for q in res["quantities"] if not any("停止の観測" in c["check"] for c in q["checks"])]))

    # forge と同じ形: 開始の場 (0) に残差の行・節点・流量が無い
    run = build("p3_forge_like", end=42000, start_obs=False)
    r40 = settle.judge(run, 40000)
    r42 = settle.judge(run, 42000)
    rq = qres(r40, "residual")
    chk("P3 forge と同じ形: 終端 40,000 は直前の窓の格子点 0 の残差が無く判定不能 (データ不足)",
        r40["verdict"] == "INDETERMINATE" and cats(rq) == {settle.CAT_DATA}, f"{r40['verdict']} {rq['reasons'][:1]}")
    chk("P3 forge と同じ形: 終端 42,000 で SETTLED", r42["verdict"] == "SETTLED", r42["verdict"])

    # 壁の熱流: 一定のゼロは正例、許容の内のゼロの交差は正例 (§6 V5 の壁の熱流)
    run = build("p4_wall_zero", Q_fn=lambda s: np.zeros(NI))
    chk("P4 壁の熱量が一定のゼロは SETTLED", qres(settle.judge(run), "wall_heat_nodes")["verdict"] == "SETTLED")
    run = build("p5_wall_zero_cross", Q_fn=lambda s: np.full(NI, 1.0 if (s // DT) % 2 else -1.0))
    q = qres(settle.judge(run), "wall_heat_nodes")
    chk("P5 許容の内のゼロの交差 (±1 W/rad) は SETTLED", q["verdict"] == "SETTLED", f"{q['verdict']} worst {q['values'].get('worst_window')}")

    # --scan: 過渡の後に一定になる系列で最初に SETTLED になる step
    def fn(s):
        v = dict(RES_BASE)
        f = 10 ** (-2.0 * min(s, 30000) / 30000)
        return {k: x * f for k, x in v.items()}
    run = build("p6_scan", end=90000, resid_fn=fn)
    sc = settle.scan(run, DT)
    chk("P6 --scan: 過渡 (0〜30,000 で 2 桁下がる) の後、最初に SETTLED になる step が 30,000 より後",
        sc["first_settled"] is not None and sc["first_settled"] > 30000, f"first {sc['first_settled']}")
    rc, out = run_cli(run, "--scan", "2000")
    chk("P6 CLI --scan の出力に最初の step", f"(--scan 2000): {sc['first_settled']}" in out, out.splitlines()[-2])
    rc, out = run_cli(run, "--scan", "3000")
    chk("P6 --scan が格子の間隔の倍数でない → 終了コード 3", rc == 3, out.splitlines()[-1])


def test_negative_v2():
    print("=== V2 の負例 (SETTLED にしない・理由の分類) ===")
    # 1. 11 点だけの短い履歴
    run = build("n01_short", end=20000)
    res = settle.judge(run)
    q = qres(res, "residual")
    chk("N1 11 点だけの短い履歴 → INDETERMINATE、残差はデータ不足", res["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_DATA},
        f"{res['verdict']} {q['reasons'][:1]}")

    # 2. 窓に大きな欠測 (26,000〜34,000 の観測が無い)
    run = build("n02_gap", skip=tuple(range(26000, 34001, DT)))
    res = settle.judge(run)
    bad = {q["name"]: cats(q) for q in res["quantities"]}
    chk("N2 窓に大きな欠測 → INDETERMINATE、全量がデータ不足", res["verdict"] == "INDETERMINATE"
        and all(c == {settle.CAT_DATA} for c in bad.values()), str(bad))

    # 3. 必要な列が欠ける (rms_ro だけ)
    run = build("n03_ro_only", resid_cols=["rms_ro"])
    res = settle.judge(run)
    q = qres(res, "residual")
    chk("N3 残差の列が rms_ro だけ → INDETERMINATE、残差はデータ不足 (列が無い)", res["verdict"] == "INDETERMINATE"
        and cats(q) == {settle.CAT_DATA} and any("列が無い" in r["text"] for r in q["reasons"]), str(q["reasons"][:1]))

    # 4. 末尾の 1 点だけ 10⁸ 倍
    def spike(s):
        v = dict(RES_BASE)
        if s == 39999:
            v = {k: x * 1e8 for k, x in v.items()}
        return v
    run = build("n04_last_spike", resid_fn=spike)
    res = settle.judge(run)
    q = qres(res, "residual")
    chk("N4 末尾の 1 点だけ 10⁸ 倍 → NOT_SETTLED、残差は閾値の超過 (跳ね)", res["verdict"] == "NOT_SETTLED"
        and cats(q) == {settle.CAT_THRESHOLD} and any("跳ね" in c for c in failed_checks(q)), str(failed_checks(q)[:3]))

    # 5. 短い反復のバースト (格子の行に当たらない 33,100〜33,400 だけ 10³ 倍)
    def burst(s):
        v = dict(RES_BASE)
        if 33100 <= s <= 33400:
            v = {k: x * 1e3 for k, x in v.items()}
        return v
    run = build("n05_burst", resid_fn=burst, extra_resid_rows=range(33100, 33401, 50))
    res = settle.judge(run)
    q = qres(res, "residual")
    fc = failed_checks(q)
    chk("N5 短い反復のバースト → NOT_SETTLED、跳ねの検査だけが落とす (格子の系列の傾き・振幅は通る)",
        res["verdict"] == "NOT_SETTLED" and cats(q) == {settle.CAT_THRESHOLD} and fc and all("跳ね" in c for c in fc), str(fc[:3]))

    # 6. 持続的な上昇 (窓 C の全行 (行 19,999 以降) が直前の窓より 12 % 高い水準で一定 — C の中の傾き・振幅は 0)
    def rise(s):
        v = dict(RES_BASE)
        if s >= 19999:
            v = {k: x * 1.12 for k, x in v.items()}
        return v
    run = build("n06_rise", resid_fn=rise)
    res = settle.judge(run)
    q = qres(res, "residual")
    fc = failed_checks(q)
    chk("N6 持続的な上昇 (+12 %) → NOT_SETTLED、上昇の検査だけが落とす (傾き・振幅・跳ねは通る)", res["verdict"] == "NOT_SETTLED"
        and cats(q) == {settle.CAT_THRESHOLD} and fc and all("上昇" in c for c in fc),
        f"{fc[:2]} rise {q['values']['columns']['rms_ro'].get('rise_dec')}")

    # 7. 活動中の乱流の列が全ゼロ
    def kzero(s):
        v = dict(RES_BASE)
        v["rms_roK"] = 0.0
        return v
    run = build("n07_k_zero", resid_fn=kzero)
    res = settle.judge(run)
    q = qres(res, "residual")
    chk("N7 活動中の乱流の列 (rms_roK) が全ゼロ → INDETERMINATE、データ不足 (恒等ゼロ)", res["verdict"] == "INDETERMINATE"
        and cats(q) == {settle.CAT_DATA} and any("恒等ゼロ" in r["text"] for r in q["reasons"]), str(q["reasons"][:1]))

    # 8. 基準のファイルが無い
    run = build("n08_no_criteria", criteria=False)
    res = settle.judge(run)
    rc, out = run_cli(run)
    chk("N8 基準のファイルが無い → INDETERMINATE (評価対象が未設定)、置き場所と量の種類の一覧", res["verdict"] == "INDETERMINATE"
        and res["reason_category"] == settle.CAT_UNSET and "settle_criteria.yaml" in res["reason"] and "wall_heat_nodes" in res["reason"]
        and rc == 3 and "評価対象が未設定" in out.splitlines()[-1], out.splitlines()[-1])

    # 9. 目的の量が無い (残差と流量だけ)
    c = base_criteria()
    c["quantities"] = c["quantities"][:2]
    run = build("n09_no_target", criteria=c)
    res = settle.judge(run)
    chk("N9 目的の量が無い → INDETERMINATE (評価対象が未設定)", res["verdict"] == "INDETERMINATE"
        and res["reason_category"] == settle.CAT_UNSET and "目的の量" in res["reason"], res.get("reason", "")[:80])
    c = base_criteria()
    c["quantities"] = [q for q in c["quantities"] if q["kind"] != "mass_balance"]
    run = build("n09b_no_mass", criteria=c)
    res = settle.judge(run)
    chk("N9b 流量の収支が無い → INDETERMINATE (評価対象が未設定)", res["verdict"] == "INDETERMINATE"
        and res["reason_category"] == settle.CAT_UNSET, res.get("reason", "")[:80])

    # 10. 知らない種類・キー・閾値の欠け
    for tag, mut, word in (
            ("種類", lambda c: c["quantities"].append({"kind": "exit_pitot", "thresholds": {"tol": 1}}), "知らない種類"),
            ("量のキー", lambda c: c["quantities"][2].update(eta_rng=[0.1, 0.2]), "知らないキー"),
            ("最上位のキー", lambda c: c.update(windw={"W": 1}), "知らないキー"),
            ("閾値の欠け", lambda c: c["quantities"][0]["thresholds"].pop("excursion_dec"), "欠けている"),
            ("閾値のキー", lambda c: c["quantities"][1]["thresholds"].update(mean_abs_rel=1), "知らないキー")):
        c = base_criteria()
        mut(c)
        run = build(f"n10_{len(tag)}_{abs(hash(tag)) % 1000}", criteria=c, end=4000)
        try:
            settle.judge(run)
            raised = None
        except settle.CriteriaError as e:
            raised = str(e)
        rc, out = run_cli(run)
        chk(f"N10 {tag} → CriteriaError ({word})・CLI は終了コード 3 で SETTLED にしない",
            raised is not None and word in raised and rc == 3 and "基準のファイルの誤り" in out.splitlines()[-1], str(raised)[:90])
    # 重複キー
    run = build("n10_dupkey", end=4000)
    txt = (run / settle.CRITERIA_FILE).read_text().replace("hold: 40000", "hold: 40000\nhold: 20000")
    (run / settle.CRITERIA_FILE).write_text(txt)
    try:
        settle.load_criteria(run)
        ok = False
    except settle.CriteriaError as e:
        ok = "重複キー" in str(e)
    chk("N10 YAML の重複キー → CriteriaError (後勝ちにしない)", ok)

    # 11. 入力のファイルが無い (boundary_massflux.csv)
    run = build("n11_no_massflux", no_massflux=True)
    res = settle.judge(run)
    q = qres(res, "mass_balance")
    chk("N11 入力のファイル (boundary_massflux.csv) が無い → INDETERMINATE、流量はデータ不足", res["verdict"] == "INDETERMINATE"
        and cats(q) == {settle.CAT_DATA} and "無い" in q["reasons"][0]["text"], q["reasons"][0]["text"][:80])
    run = build("n11b_no_nodes", end=40000)
    (run / "res_nodes_40000.h5").unlink()
    res = settle.judge(run, 40000)
    chk("N11b 終端の res_nodes が無い → 節点の量 3 つともデータ不足", all(cats(qres(res, n)) == {settle.CAT_DATA}
        for n in ("exit_core_mach", "line_mach", "exit_flow_angle")) and res["verdict"] == "INDETERMINATE")

    # 12. 総熱量は一定で局所の分布が変わる壁
    def swap(s):
        Q = 100.0 + np.arange(NI, dtype=float)
        d = 200.0 * math.sin(s / 7000.0)
        Q[10] += d
        Q[20] -= d
        return Q
    run = build("n12_wall_swap", Q_fn=swap)
    res = settle.judge(run)
    q = qres(res, "wall_heat_nodes")
    sums = [float(np.sum(swap(s))) for s in res["window_steps"]["C"]]
    chk("N12 総熱量は一定 (窓の ΣQ の幅 0) で局所の分布が変わる壁 → NOT_SETTLED、閾値の超過",
        max(sums) - min(sums) < 1e-9 and res["verdict"] == "NOT_SETTLED" and cats(q) == {settle.CAT_THRESHOLD},
        f"ΣQ 幅 {max(sums) - min(sums):.3g}, worst {q['values'].get('worst_window', {}).get('ratio')}")

    # 13. 位置が動くが最大偏差と片振幅が一定の線
    def moving(s):
        c = 8 + (s // 4000) % 6
        return MD * (1.0 + 0.01 * bump(np.arange(NI) + 2, c + 2, 0.5))
    run = build("n13_line_moving", M_fn=moving)
    res = settle.judge(run)
    q = qres(res, "line_mach")
    fc = failed_checks(q)
    v = q["values"]
    chk("N13 位置が動くが maxabs・halfamp が一定の線 → NOT_SETTLED、位置ごとの検査だけが落とす",
        res["verdict"] == "NOT_SETTLED" and cats(q) == {settle.CAT_THRESHOLD} and fc == ["位置ごとの d の幅の最大 [%pt]"],
        f"fail {fc}; maxabs 幅 {max(v['maxabs']) - min(v['maxabs']):.2e}, halfamp 幅 {max(v['halfamp']) - min(v['halfamp']):.2e}, "
        f"pos {v.get('pos_span_max')}")
    chk("N13 出口コア M はこの線の変化に影響されない (SETTLED)", qres(res, "exit_core_mach")["verdict"] == "SETTLED")


def test_negative_extra():
    print("=== そのほかの負例 ===")
    def nanres(s):
        v = dict(RES_BASE)
        if s == 35000:
            v["rms_roe"] = float("nan")
        return v
    run = build("x01_nan_resid", resid_fn=nanres)
    res = settle.judge(run)
    q = qres(res, "residual")
    chk("X1 残差に NaN → DIVERGED (非有限)", res["verdict"] == "DIVERGED" and cats(q) == {settle.CAT_NONFINITE}, str(q["reasons"][:1]))
    run = build("x02_nan_nodes", nan_nodes_at=30000)
    res = settle.judge(run)
    chk("X2 節点の値に NaN → DIVERGED", res["verdict"] == "DIVERGED" and cats(qres(res, "line_mach")) == {settle.CAT_NONFINITE})

    def partzero(s):
        v = dict(RES_BASE)
        if s > 30000:
            v["rms_roOmega"] = 0.0
        return v
    run = build("x03_partly_zero", resid_fn=partzero)
    res = settle.judge(run)
    chk("X3 残差の列の一部がゼロ → INDETERMINATE (データ不足)", res["verdict"] == "INDETERMINATE"
        and cats(qres(res, "residual")) == {settle.CAT_DATA})

    run = build("x04_dup_rows")
    with open(run / "residual_history.csv", "a") as f:
        f.write("30001,-1,outer_end," + ",".join(repr(RES_BASE[c]) for c in RES_COLS) + ",0,0,0,0,0\n")
    res = settle.judge(run, 40000)
    chk("X4 残差の step が逆行・重複 → INDETERMINATE (データ不足)", res["verdict"] == "INDETERMINATE"
        and any("単調でない" in r["text"] or "重複" in r["text"] for r in qres(res, "residual")["reasons"]))

    def extra_col(s):
        v = dict(RES_BASE)
        v["rms_roXi"] = 1e-3
        return v
    run = build("x05_unlisted", resid_fn=extra_col, resid_cols=RES_COLS + ["rms_roXi"])
    res = settle.judge(run)
    chk("X5 基準に無い残差の列が値を持つ → INDETERMINATE (未設定)", res["verdict"] == "INDETERMINATE"
        and settle.CAT_UNSET in cats(qres(res, "residual")))

    def stopspike(s):
        v = dict(RES_BASE)
        if s == 80000:
            v = {k: x * 1e3 for k, x in v.items()}
        return v
    run = build("x06_stop_spike", end=80000, stop=80001, resid_fn=stopspike)
    res = settle.judge(run)
    q = qres(res, "residual")
    chk("X6 停止の観測 (80,001; 残差の行 80,000) だけ 10³ 倍 → NOT_SETTLED (停止の観測の検査)", res["verdict"] == "NOT_SETTLED"
        and all("停止の観測" in c for c in failed_checks(q)) and failed_checks(q), str(failed_checks(q)[:2]))

    def stopline(s):
        M = np.full(NI, MD)
        if s == 80001:
            M[12] *= 1.001
        return M
    run = build("x07_stop_line", end=80000, stop=80001, M_fn=stopline)
    res = settle.judge(run)
    q = qres(res, "line_mach")
    chk("X7 停止の観測で線が動く → NOT_SETTLED (停止の観測の検査だけ)", res["verdict"] == "NOT_SETTLED"
        and failed_checks(q) and all("停止の観測" in c for c in failed_checks(q)), str(failed_checks(q)[:2]))
    run = build("x07b_stop_missing_nodes", end=80000, stop=80001)
    (run / "res_nodes_80001.h5").unlink()
    res = settle.judge(run)
    chk("X7b 停止の step の節点の出力が無い → INDETERMINATE (データ不足)", res["verdict"] == "INDETERMINATE"
        and cats(qres(res, "exit_core_mach")) == {settle.CAT_DATA})

    # 系列の CSV (V3(a) の再生用): 位置ごと・節点ごとは判定不能
    run = build("x08_series")
    rows = ["step,maxabs,halfamp,exit_core_M,r,mdot_in,mdot_out"] + [f"{s},0.2,0.1,{MD},0.0,10.0,10.0" for s in range(0, 40001, DT)]
    (run / "series.csv").write_text("\n".join(rows) + "\n")
    c = base_criteria()
    for q in c["quantities"]:
        if q["kind"] == "line_mach":
            q["series_csv"] = {"path": "series.csv", "step_col": "step", "maxabs_col": "maxabs", "halfamp_col": "halfamp"}
        if q["kind"] == "exit_core_mach":
            q["series_csv"] = {"path": "series.csv", "step_col": "step", "value_col": "exit_core_M"}
        if q["kind"] == "mass_balance":
            q["series_csv"] = {"path": "series.csv", "step_col": "step", "value_col": "r"}
    (run / settle.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    res = settle.judge(run)
    ql = qres(res, "line_mach")
    qm = qres(res, "mass_balance")
    chk("X8 series_csv の line_mach: maxabs・halfamp は判定、位置ごとは判定不能 (データ不足)", ql["input"] == "series_csv"
        and ql["verdict"] == "INDETERMINATE" and cats(ql) == {settle.CAT_DATA} and not failed_checks(ql), str(ql["reasons"][:1]))
    chk("X8 series_csv の mass_balance (r だけ): mdot_change_rel は判定不能 (データ不足)", qm["verdict"] == "INDETERMINATE"
        and cats(qm) == {settle.CAT_DATA})
    chk("X8 series_csv の exit_core_mach は判定できる (SETTLED)", qres(res, "exit_core_mach")["verdict"] == "SETTLED")

    # 節点の表の部分集合が足りない (η を挟む節点の片方が無い) → 線は判定不能
    run = build("x09_subset")
    keep = [i * NJ + j for i in range(NI) for j in range(NJ) if not (j == 1 and i == 12)]
    for s in range(0, 40001, DT):
        write_nodes(run, s, default_M(s), ids=keep)
    res = settle.judge(run)
    q = qres(res, "line_mach")
    chk("X9 節点の表で列 12 の j = 1 が欠ける → 線は判定不能 (隣り合わない)", q["verdict"] == "INDETERMINATE"
        and any("隣り合わない" in r["text"] for r in q["reasons"]), str(q["reasons"][:1]))
    chk("X9 出口の列は揃っているので出口コア M は SETTLED", qres(res, "exit_core_mach")["verdict"] == "SETTLED")

    # 流量の収支の閾値
    run = build("x10_mass_drift", mdot_fn=lambda s: (10.0, 10.0 * (1 - 1e-3 * s / 40000)))
    res = settle.judge(run)
    q = qres(res, "mass_balance")
    chk("X10 出口の流量が 0.1 % ずれていく → NOT_SETTLED (閾値の超過)", res["verdict"] == "NOT_SETTLED"
        and cats(q) == {settle.CAT_THRESHOLD}, str(failed_checks(q)))
    run = build("x11_mass_physid")
    c = base_criteria()
    c["quantities"][1]["outlet"] = [7]
    (run / settle.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    res = settle.judge(run)
    chk("X11 基準の physID の列が CSV に無い → INDETERMINATE (データ不足)", cats(qres(res, "mass_balance")) == {settle.CAT_DATA})

    # 壁の h5 の値に NaN (適用外) → 判定不能
    run = build("x12_wall_nan", Q_fn=lambda s: np.full(NI, np.nan))
    res = settle.judge(run)
    chk("X12 iface_Qf_eff が NaN (適用外) → INDETERMINATE (データ不足)", cats(qres(res, "wall_heat_nodes")) == {settle.CAT_DATA})


def test_verdict_file():
    print("=== SETTLE_VERDICT.json ===")
    run = build("v01_file", end=80000, stop=80001)
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = check_settled.main([str(run)])
    p = run / settle.VERDICT_FILE
    d = json.loads(p.read_text())
    crit_sha = settle._sha256_bytes((run / settle.CRITERIA_FILE).read_bytes())
    files = {Path(f["path"]).name for f in d["files"]}
    chk("V1 SETTLE_VERDICT.json: 基準の中身と sha256・判定区間・終端・停止の step・使ったファイルの一覧",
        rc == 0 and d["criteria"]["sha256"] == crit_sha and d["criteria"]["content"]["schema"] == settle.SCHEMA
        and d["terminal_step"] == 80000 and d["stop_step"] == 80001 and d["grid"]["start_step"] == 0
        and {"residual_history.csv", "boundary_massflux.csv", "res_nodes_80001.h5", "res_wall_3_80000.h5", "STOPPED_AT",
             settle.CRITERIA_FILE} <= files and all(len(f["sha256"]) == 64 for f in d["files"]) and d["provisional"] is True,
        f"files {len(files)}")
    q = {x["name"]: x for x in d["quantities"]}
    chk("V1 各量に verdict・使った値・理由の分類", all({"verdict", "values", "reasons", "categories", "checks"} <= set(x) for x in q.values()))
    chk("V1 壁の面積の取り方を記録", "A_i" in q["wall_heat_nodes"]["values"]["walls"][0]["area_method"])


def test_area():
    print("=== 壁の節点の面積 ===")
    coord = np.array([[0.0, 1.0, 0], [1.0, 1.0, 0], [2.0, 2.0, 0]])
    conne = np.array([2, 2, 0, 1, 2, 2, 1, 2])
    A, _ = settle.wall_node_area(coord, conne, True)
    l2 = math.sqrt(2.0)
    exp = np.array([0.5 * 1.0 * (3 + 1) / 4, 0.5 * 1.0 * (3 + 1) / 4 + 0.5 * l2 * (3 * 1 + 2) / 4, 0.5 * l2 * (3 * 2 + 1) / 4])
    chk("A1 軸対称: A_i = Σ (|e|/2)(3 r_i + r_j)/4 (Σ A = ∫ r dl)", np.allclose(A, exp) and abs(A.sum() - (1.0 + l2 * 1.5)) < 1e-12,
        str(A))
    A2, _ = settle.wall_node_area(coord, conne, False)
    chk("A2 平面: A_i = Σ |e|/2", np.allclose(A2, [0.5, 0.5 + l2 / 2, l2 / 2]))


def main():
    print(f"numpy {np.__version__}, python {sys.version.split()[0]}, 作業 {TMP}")
    test_area()
    test_positive()
    test_negative_v2()
    test_negative_extra()
    test_verdict_file()
    shutil.rmtree(TMP, ignore_errors=True)
    print(f"\n{NCHK - len(FAILS)}/{NCHK} OK")
    if FAILS:
        print("FAIL:", *FAILS, sep="\n  ")
    print(f"VERDICT: {'PASS' if not FAILS else 'FAIL'}")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
