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

**wt_nozzle_v2 の形** (§5.1 #5d; 窓 10,000・保持 0・残差 `series: raw`・壁の熱流 `measure: range`): 毎 step の残差の行 (forge と同じ形)
で V2 の負例 (末尾の跳ね・短いバースト・持続的な上昇・一部ゼロ・欠測) が落ち一定系列が通ること、codex diagnose settle-v2-span10k の
再現 (途中の 1 行の欠け・2,000 ごとの行だけ → raw は判定不能、grid では SETTLED)、区間の端の行、内反復の行の混入、壁の区間の幅
(codex の例 [0.25, 1.75, 1, 1, 1, 1])。**既定 (grid・terminal_diff) は変更前 (`GRID_BASE_COMMIT` の settle.py を git から読む) と
判定の結果がビット同一** (読めなければその検査は NG — 判定不能を合格にしない)。

**道具の選択と本段の起動との照合** (§4.9 の 2・§5.1 #7c; S1〜S11): `select_judge` の mode (main・pre_main・unconfirmed・manual・legacy・invalid)、
照合の不一致 (設定・基準の sha256・記録の fnv の欠け、**起動の記録の欠け・最後の起動の不一致**) は check_settled.py のまま判定不能
(旧の道具へ落とさない)、記録が無い run は stage_manifest.json で前の段と確かめたときだけ旧の道具、`check_settled.py` は既定で照合
(`--no-bind` は VERDICT と結果に書く)、`--select`・`--bind` (result 段レビュー M1)。
**result 段レビュー 2 回目** (§5.1 #9・#10): 基準が無効 (mode invalid; stage を main 以外に差し替え・評価対象が未設定) と照合の不一致
(基準の sha256) は CLI でも判定しない (S11b〜S11f)、v3 の非有限の見張りは区間の毎 step の行 (欠け・重複・逆行なし) と必要な列を要る
(F1〜F10)。負例は修正前 (`PRE_REVIEW2_COMMIT` の settle.py・check_settled.py を git から読む) が SETTLED だったこと (検出力) も確かめる
(読めなければその検査は NG)。
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


def write_residual(run, steps, fn, cols=RES_COLS, name="residual_history.csv", inner_rows=False):
    """inner_rows: 各 step に内反復の行 (inner_begin・inner_end、値は 10⁸ 倍) を足す (outer_end だけを読むことの試験)。"""
    with open(run / name, "w") as f:
        f.write("step,inner,phase," + ",".join(cols + DQ_COLS) + "\n")
        for s in steps:
            v = fn(s)
            txt = ",".join(repr(v.get(c, 0.0)) for c in cols)
            f.write(f"{s},-1,outer_begin," + txt + ",0,0,0,0,0\n")
            if inner_rows:
                big = ",".join(repr(v.get(c, 0.0) * 1e8) for c in cols)
                f.write(f"{s},0,inner_begin," + big + ",0,0,0,0,0\n")
                f.write(f"{s},0,inner_end," + big + ",0,0,0,0,0\n")
            f.write(f"{s},-1,outer_end," + txt + ",0,0,0,0,0\n")


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
            {"kind": "line_mach", "eta": 0.15, "x_range_m": [3 * S, 20 * S], "n": 2001, "Md": MD,
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
          extra_resid_rows=(), stopped_at=True, rows_mode="every10", inner_rows=False, dup_resid_rows=()):
    """合成の run を作る。skip: 観測を置かない格子点 (残差の行 s−1・節点・流量・壁)。
    rows_mode: 残差の outer_end の行 — every10 (10 step ごと + 格子点の行 s−1; 今までの試験)、dense (毎 step; forge と同じ)、
    grid_only (格子点の行 s−1 だけ; codex の「2,000 ごとの行だけ」)。dup_resid_rows: 2 回書く行。"""
    run = TMP / name
    if run.exists():
        shutil.rmtree(run)
    run.mkdir()
    grid = [s for s in range(0, end + 1, DT) if (start_obs or s > 0) and s not in skip]
    obs = grid + ([stop] if stop is not None else [])
    last_row = (stop if stop is not None else end) - 1
    resid_fn = resid_fn or (lambda s: dict(RES_BASE))
    step10 = {"every10": 10, "dense": 1}.get(rows_mode)
    rows = (set(range(0, last_row + 1, step10)) if step10 else set()) | {s - 1 for s in obs}
    rows |= set(extra_resid_rows)
    if start_obs:
        rows.add(-1)
    rows = sorted(r for r in rows if r <= last_row and (r + 1) not in skip and r not in skip_resid_rows)
    rows = sorted(rows + [r for r in dup_resid_rows if r in rows])
    write_residual(run, rows, resid_fn, cols=resid_cols, inner_rows=inner_rows)
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


def test_axis_weight():
    print("=== 線の補間に軸の節点が入る格子 (codex diagnose 2026-10-16 Major 3) ===")
    # 合成の格子の第 1 層は η₁ = 1/9 ≈ 0.111。η 0.03 の線は軸の節点と第 1 層の間で補間され、軸の節点の重みは約 0.73
    crit = base_criteria()
    for q in crit["quantities"]:
        if q["kind"] == "line_mach":
            q["eta"] = 0.03
    run = build("x20_axis_weight", criteria=crit)
    res = settle.judge(run)
    q = qres(res, "line_mach")
    chk("X20 η 0.03 < η₁ の線 → INDETERMINATE (格子と線の定義の不整合)、軸の節点の重みを記録",
        q["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_GRID} and q["values"].get("axis_weight_max", 0) > 0.5
        and res["verdict"] == "INDETERMINATE", f"{q['verdict']} {cats(q)} w {q['values'].get('axis_weight_max')}")
    run = build("x21_axis_weight_ok")
    q = qres(settle.judge(run), "line_mach")
    chk("X21 η 0.15 (第 1・2 層の間) は軸の節点の重み 0 で SETTLED、補間の対を記録",
        q["verdict"] == "SETTLED" and q["values"].get("axis_weight_max") == 0.0 and set(q["values"].get("interp_pairs", {})) == {"1-2"},
        f"{q['verdict']} {q['values'].get('interp_pairs')} {q['values'].get('axis_weight_max')}")


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


def test_select_bind():
    """道具の選択 `select_judge` と `check_settled.py --select` / `--bind` (plan §4.9 の 2、codex diagnose settle-migration-7 M1・M5)。"""
    print("=== 道具の選択と本段の起動との照合 (select_judge・--select・--bind) ===")
    import stage_manifest as SM
    seg = {"stage": "main", "start_record": "MAIN_STAGE_STARTED", "config_fnv_key": "main_config.fnv",
           "criteria_sha256_key": "settle_criteria.sha256"}
    cfg = "time: {last: {nStepOuter: 100000}}\n"

    def mk(name, *, notice=True, launch="match", crit_seg=True):
        run = build(name, end=40000, criteria=base_criteria(grid={"out_step_interval": DT, "start_step": 0,
                                                                  "segment": seg if crit_seg else None}))
        (run / "solverConfig.yaml").write_text(cfg)
        if notice:
            (run / "MAIN_STAGE_STARTED").write_text(json.dumps({
                "event": "main_stage_started",
                "main_config": {"file": "solverConfig.yaml", "sha256": settle._sha256_bytes(cfg.encode()), "nStepOuter": 100000,
                                "fnv": SM.fnv1a64(cfg)},
                "settle_criteria": {"file": settle.CRITERIA_FILE, "sha256": settle._sha256_bytes((run / settle.CRITERIA_FILE).read_bytes())}}))
        if launch is not None:
            fnv = SM.fnv1a64(cfg) if launch == "match" else launch
            (run / "forge_launches.jsonl").write_text(json.dumps({"cfg_fnv": "aaaa"}) + "\n" + json.dumps({"cfg_fnv": fnv}) + "\n")
        return run

    def cli(run, *args):
        buf = io.StringIO()
        with redirect_stdout(buf):
            try:
                rc = check_settled.main([*args, str(run)])
            except SystemExit as e:
                rc = e.code
        return rc, buf.getvalue()

    r = mk("sel_main")
    s = settle.select_judge(r)
    chk("S1 本段 (記録・今の設定・最後の起動・基準の sha256 が一致) → check_settled.py・問題なし",
        s["tool"] == settle.TOOL_NEW and s["mode"] == "main" and s["problems"] == [], str(s["problems"]))
    rc, out = cli(r, "--select")
    chk("S1 --select: 1 行目 check_settled.py・終了コード 0", rc == 0 and out.splitlines()[0] == "check_settled.py", out[:120])
    rc, out = cli(r, "--bind", "--no-write")
    chk("S1 --bind: 照合を表示して判定 (SETTLED・終了コード 0)", rc == 0 and "照合:" in out and "VERDICT: SETTLED" in out, out[-200:])
    rc, out = cli(r, "--no-write")
    chk("S1 引数なし (既定で照合): 照合を表示して判定 (SETTLED・終了コード 0)", rc == 0 and "照合:" in out and "VERDICT: SETTLED" in out
        and "照合を外した" not in out, out[-200:])
    # --- result 段レビュー M1 (§5.1 #7c): 本段の記録がある run は旧の道具へ切り替えない ---
    r = mk("sel_nolaunch", launch=None)
    s = settle.select_judge(r)
    rc, out = cli(r, "--no-write")
    chk("S2 本段の記録があり起動の記録 forge_launches.jsonl が無い → check_settled.py の照合の不一致 (旧の道具へ切り替えない)・"
        "引数なしで INDETERMINATE (終了コード 3)", s["tool"] == settle.TOOL_NEW and s["mode"] == "main"
        and any("起動の記録" in p and "無い" in p for p in s["problems"]) and rc == 3 and "照合の不一致" in out, str(s["problems"]))
    r = mk("sel_pre", notice=False)
    s = settle.select_judge(r)
    rc, out = cli(r, "--no-write")
    chk("S3 本段の記録が無く stage_manifest.json も無い (前の段と確かめられない) → check_settled.py の INDETERMINATE (記録の欠け)",
        s["tool"] == settle.TOOL_NEW and s["mode"] == "unconfirmed" and s["problems"] and rc == 3 and "記録の欠け" in out, out[-200:])
    (r / "stage_manifest.json").write_text(json.dumps({"manifest_version": 3, "stages": [{"order": 0, "tag": "S1_soft"}, {"order": 1, "tag": "S2_mid"}]}))
    s = settle.select_judge(r)
    rc, out = cli(r, "--bind", "--no-write")
    rc2, out2 = cli(r, "--select")
    chk("S3b 本段の記録が無く stage_manifest.json の段 (S1_soft・S2_mid) に本段が無い = 前の段と確かめた → check_convergence.py; "
        "check_settled.py は判定せず INDETERMINATE (本段の起動でない)", s["tool"] == settle.TOOL_OLD and s["mode"] == "pre_main"
        and rc == 3 and "本段の起動でない" in out and out2.splitlines()[0] == "check_convergence.py", out[-160:])
    (r / "stage_manifest.json").write_text(json.dumps({"manifest_version": 3, "stages": [{"order": 0, "tag": "S1_soft"}, {"order": 1, "tag": "S2_mid"},
                                                                                     {"order": 2, "tag": "main"}]}))
    s = settle.select_judge(r)
    chk("S3c 本段の記録が無く stage_manifest.json に本段 (main) がある → 前の段と確かめられない (unconfirmed、旧の道具へ切り替えない)",
        s["tool"] == settle.TOOL_NEW and s["mode"] == "unconfirmed" and any("本段の段" in p for p in s["problems"]), str(s["problems"]))
    (r / "stage_manifest.json").write_text("{not json")
    s = settle.select_judge(r)
    chk("S3d stage_manifest.json が読めない → unconfirmed", s["mode"] == "unconfirmed" and s["tool"] == settle.TOOL_NEW, str(s["problems"]))
    r = mk("sel_not_main", launch="bbbbbbbbbbbbbbbb")
    s = settle.select_judge(r)
    rc, out = cli(r, "--no-write")
    chk("S4 最後の起動の cfg_fnv が本段の記録と違う → check_settled.py の照合の不一致 (以前は not_main で check_convergence.py)・"
        "引数なしで INDETERMINATE", s["tool"] == settle.TOOL_NEW and s["mode"] == "main"
        and any("最後の forge の起動" in p for p in s["problems"]) and rc == 3 and "VERDICT: INDETERMINATE (本段の起動との照合の不一致)" in out,
        str(s["problems"]))
    r = mk("sel_rec_fnv")
    n = json.loads((r / "MAIN_STAGE_STARTED").read_text())
    n["main_config"]["fnv"] = "0123456789abcdef"
    (r / "MAIN_STAGE_STARTED").write_text(json.dumps(n))
    s = settle.select_judge(r)
    rc, out = cli(r, "--no-write")
    chk("S4b レビューの再現: 本段の記録の main_config.fnv だけを変える → 以前は not_main・problems [] → 今は main・照合の不一致・"
        "INDETERMINATE (終了コード 3)", s["tool"] == settle.TOOL_NEW and s["mode"] == "main" and len(s["problems"]) >= 2
        and rc == 3, str(s["problems"])[:300])
    rc, out = cli(r, "--no-bind", "--no-write", "--json")
    jb = json.loads(out[out.index("{"):out.rindex("}") + 1]) if "{" in out else {}
    chk("S4c --no-bind: 照合を外して判定し、VERDICT の行と結果の JSON (binding.skipped) に外したことを書く",
        rc == 0 and "VERDICT: SETTLED" in out and "照合を外した" in out.splitlines()[-1] and (jb.get("binding") or {}).get("skipped") is True,
        out.splitlines()[-1])
    r = mk("sel_cfg_edit", launch=None)
    (r / "solverConfig.yaml").write_text(cfg + "# 書き換え\n")
    s = settle.select_judge(r)
    rc, out = cli(r, "--bind")
    chk("S5 起動の記録なしで今の設定が記録と違う → 本段として照合の不一致 3 件 (起動の記録・fnv・sha256; 旧の道具へ落とさない)、"
        "INDETERMINATE・SETTLE_VERDICT.json を書かない", s["tool"] == settle.TOOL_NEW and s["mode"] == "main" and len(s["problems"]) == 3
        and rc == 3 and "照合の不一致" in out and not (r / settle.VERDICT_FILE).exists(), str(s["problems"])[:200])
    rc, out = cli(r)
    chk("S5b 引数なし (通常のコマンド) でも照合の不一致なら判定せず SETTLE_VERDICT.json を書かない",
        rc == 3 and "照合の不一致" in out and not (r / settle.VERDICT_FILE).exists(), out[-160:])
    r = mk("sel_crit_edit")
    with open(r / settle.CRITERIA_FILE, "a") as f:
        f.write("# 書き換え\n")
    s = settle.select_judge(r)
    chk("S6 基準の sha256 が記録と違う → 照合の不一致", s["mode"] == "main" and any("settle_criteria.sha256" in p for p in s["problems"]),
        str(s["problems"]))
    r = mk("sel_nofnv")
    n = json.loads((r / "MAIN_STAGE_STARTED").read_text())
    del n["main_config"]["fnv"]
    (r / "MAIN_STAGE_STARTED").write_text(json.dumps(n))
    s = settle.select_judge(r)
    chk("S7 記録に main_config.fnv が無い → 照合できない (問題)", s["mode"] == "main" and any("main_config.fnv が無い" in p for p in s["problems"]))
    r = mk("sel_manual", notice=False, launch=None, crit_seg=False)
    s = settle.select_judge(r)
    chk("S8 手で組んだ run (grid.segment が記録を指さず MAIN_STAGE_STARTED も無い) → check_settled.py (照合なし)",
        s["tool"] == settle.TOOL_NEW and s["mode"] == "manual")
    r = mk("sel_manual_rec", crit_seg=False)
    s = settle.select_judge(r)
    chk("S8b 基準が記録を指さなくても run に MAIN_STAGE_STARTED がある → 既定のキーで本段として照合 (一致なら問題なし)",
        s["tool"] == settle.TOOL_NEW and s["mode"] == "main" and s["problems"] == [] and "記録を指さない" in s["why"], str(s["problems"]))
    c = yaml.safe_load((r / settle.CRITERIA_FILE).read_text())
    c["notes"] = "準備の後に grid.segment を外して書き換えた"
    (r / settle.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    s = settle.select_judge(r)
    rc, out = cli(r, "--no-write")
    chk("S8c 基準の segment を外し書き換えても、記録がある run は照合する (基準の sha256 の不一致 → INDETERMINATE)",
        s["mode"] == "main" and any("settle_criteria.sha256" in p for p in s["problems"]) and rc == 3, str(s["problems"])[:200])
    r = mk("sel_legacy")
    (r / settle.CRITERIA_FILE).unlink()
    s = settle.select_judge(r)
    chk("S9 基準のファイルが無い → check_convergence.py (legacy)", s["tool"] == settle.TOOL_OLD and s["mode"] == "legacy")
    r = mk("sel_invalid")
    (r / settle.CRITERIA_FILE).write_text("schema: forge.settle_criteria/v1\nbogus: 1\n")
    s = settle.select_judge(r)
    chk("S10 基準の誤り → check_settled.py (invalid; 判定不能を返す、旧の道具へ落とさない)", s["tool"] == settle.TOOL_NEW and s["mode"] == "invalid"
        and s["problems"])
    r = mk("sel_stage")
    c = yaml.safe_load((r / settle.CRITERIA_FILE).read_text())
    c["grid"]["segment"] = {**seg, "stage": "S1"}
    (r / settle.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    s = settle.select_judge(r)
    chk("S11 grid.segment の stage が main でない → invalid (本段の基準だけを扱う)", s["mode"] == "invalid" and s["tool"] == settle.TOOL_NEW)
    # --- result 段レビュー 2 回目 M1 (§5.1 #9): invalid と照合の問題は mode によらず CLI でも判定しない ---
    rc, out = cli(r)        # 引数なし (通常のコマンド; SETTLE_VERDICT.json を書く経路)
    chk("S11b 本段の記録のある run の基準の 1 項目を差し替え (grid.segment.stage main → S1) → CLI は判定せず "
        "INDETERMINATE (基準が本段の判定に使えない)・終了コード 3・SETTLE_VERDICT.json を書かない",
        rc == 3 and out.splitlines()[-1] == "VERDICT: INDETERMINATE (基準が本段の判定に使えない)" and "判定しない" in out
        and "grid.segment.stage 'S1'" in out and not (r / settle.VERDICT_FILE).exists(), out[-240:])
    rc, out = cli(r, "--no-bind", "--no-write")
    chk("S11c 同じ run で --no-bind を明示 → 今どおり照合を外して判定する (既存の振る舞い; VERDICT の行に外したことを書く)",
        rc == 0 and "VERDICT: SETTLED" in out and "照合を外した" in out.splitlines()[-1], out.splitlines()[-1])
    old_cli, why = _load_at(PRE_REVIEW2_COMMIT, "solver_density_cuda/tools/check_settled.py", f"check_settled_pre_review2_{PRE_REVIEW2_COMMIT}")
    if old_cli is None:
        chk(f"S11d 修正前 ({PRE_REVIEW2_COMMIT}) の check_settled.py を読めた (読めなければ検出力は判定不能 — 合格にしない)", False, why)
    else:
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc_old = old_cli.main([str(r), "--no-write"])
        chk("S11d 検出力: 修正前の check_settled.py は S11b の run で SETTLED・終了コード 0 だった (レビューの再現)",
            rc_old == 0 and "VERDICT: SETTLED" in buf.getvalue(), buf.getvalue().splitlines()[-1] if buf.getvalue() else "")
    r = mk("sel_crit_sha_cli")
    with open(r / settle.CRITERIA_FILE, "a") as f:
        f.write("# 準備の後に書き換えた\n")
    rc, out = cli(r)
    chk("S11e 基準の sha256 が本段の記録と違う → CLI は判定せず INDETERMINATE (照合の不一致)・終了コード 3・SETTLE_VERDICT.json を書かない",
        rc == 3 and out.splitlines()[-1] == "VERDICT: INDETERMINATE (本段の起動との照合の不一致)" and "settle_criteria.sha256" in out
        and not (r / settle.VERDICT_FILE).exists(), out[-240:])
    r = mk("sel_unset_cli")
    c = yaml.safe_load((r / settle.CRITERIA_FILE).read_text())
    c["quantities"] = c["quantities"][:2]         # 目的の量が無い (評価対象が未設定)
    (r / settle.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    s = settle.select_judge(r)
    rc, out = cli(r)
    chk("S11f 評価対象が未設定 (invalid) → 判定せず INDETERMINATE (評価対象が未設定)・終了コード 3、置き場所と量の種類の一覧を出す",
        s["mode"] == "invalid" and rc == 3 and out.splitlines()[-1] == "VERDICT: INDETERMINATE (評価対象が未設定)"
        and "置き場所" in out and "wall_heat_nodes" in out, out[-240:])


def test_mass_denominator():
    """流量の判定の分母 (result 段レビュー M3、§5.1 #7c): 出口の流量 0・符号の反転・逆流で例外にせず理由つきの結果を返す。
    ṁ_in > 0 で ṁ_out ≤ 0 → NOT_SETTLED (r = 1 − ṁ_out/ṁ_in ≥ 1 は測れて収支が取れていない)、ṁ_in ≤ 0 → INDETERMINATE (分母が無い)。"""
    print("=== 流量の判定の分母 (出口 0・符号の反転・逆流; result 段レビュー M3) ===")

    def judged(name, **kw):
        run = build(name, **kw)
        try:
            return run, settle.judge(run, kw.get("stop")), None
        except Exception as e:  # noqa: BLE001
            return run, None, f"{type(e).__name__}: {e}"

    def out_reason(q):
        return any("ṁ_out が正でない" in r["text"] for r in q["reasons"])

    run, res, err = judged("m3a_out_zero", mdot_fn=lambda s: (10.0, 0.0))
    q = qres(res, "mass_balance") if res else {}
    chk("M3a 入口が正・出口がゼロの有限な観測 (codex の再現) → 例外にせず NOT_SETTLED (理由: ṁ_out が正でない)、"
        "ṁ_out の相対の変化は計算しない", err is None and q.get("verdict") == "NOT_SETTLED" and cats(q) == {settle.CAT_THRESHOLD}
        and out_reason(q) and not any("ṁ_out" in c["check"] and "mean" in c["check"] for c in q["checks"]), err or str(failed_checks(q)))
    rc, out = run_cli(run)
    chk("M3a check_settled.py は終了コード 1 (NOT_SETTLED) で VERDICT を出す (Traceback にしない)",
        rc == 1 and "VERDICT: NOT_SETTLED" in out and "Traceback" not in out, out[-160:])
    run, res, err = judged("m3b_out_signflip", mdot_fn=lambda s: (10.0, -10.0 if s == 30000 else 10.0))
    q = qres(res, "mass_balance") if res else {}
    chk("M3b 窓の 1 点だけ出口の流量の符号が反転 → NOT_SETTLED (理由に step 30000)", err is None and q.get("verdict") == "NOT_SETTLED"
        and out_reason(q) and any("30000" in r["text"] for r in q["reasons"]), err or str(q.get("reasons"))[:200])
    def mean_zero(s):   # 窓 20,000〜40,000 の 11 点: +5・−5 を交互に 10 点と最後に 0 → 平均がちょうど 0
        i = (s - 20000) // DT
        return (10.0, 10.0 if i < 0 else (0.0 if i == 10 else (5.0 if i % 2 == 0 else -5.0)))
    run, res, err = judged("m3c_out_mean_zero", mdot_fn=mean_zero)
    q = qres(res, "mass_balance") if res else {}
    chk("M3c 窓の出口の流量の平均がちょうど 0 (+5 と −5 が交互、以前は ZeroDivisionError) → NOT_SETTLED", err is None
        and q.get("verdict") == "NOT_SETTLED" and out_reason(q), err or "")
    run, res, err = judged("m3d_out_reverse", mdot_fn=lambda s: (10.0, -10.0))
    q = qres(res, "mass_balance") if res else {}
    chk("M3d 出口が全点で逆流 (平均が負: 以前は (max − min)/mean が負で ṁ_out の検査を素通り) → NOT_SETTLED", err is None
        and q.get("verdict") == "NOT_SETTLED" and out_reason(q), err or "")
    run, res, err = judged("m3e_stop_out_neg", end=40000, stop=40001, mdot_fn=lambda s: (10.0, -110.0 if s == 40001 else 10.0))
    q = qres(res, "mass_balance") if res else {}
    chk("M3e 停止の観測 (40,001) だけ出口が逆流 (窓 11 点の和 110 と打ち消して平均 0: 以前は停止の観測を加えた除算で "
        "ZeroDivisionError) → 例外にせず NOT_SETTLED (停止の観測の理由)", err is None
        and q.get("verdict") == "NOT_SETTLED" and any("停止の観測の出口の流量" in r["text"] for r in q["reasons"]), err or str(q.get("reasons"))[:200])
    run, res, err = judged("m3f_stop_in_zero", end=40000, stop=40001, mdot_fn=lambda s: (0.0, 10.0) if s == 40001 else (10.0, 10.0))
    q = qres(res, "mass_balance") if res else {}
    chk("M3f 停止の観測だけ入口の流量 0 (r の分母が無い) → INDETERMINATE (データ不足)、例外にしない", err is None
        and q.get("verdict") == "INDETERMINATE" and cats(q) == {settle.CAT_DATA}, err or str(q.get("reasons"))[:200])
    run, res, err = judged("m3g_in_zero", mdot_fn=lambda s: (0.0, 10.0))
    q = qres(res, "mass_balance") if res else {}
    chk("M3g 窓の入口の流量 0 → INDETERMINATE (今どおり; 理由に step)", err is None and q.get("verdict") == "INDETERMINATE"
        and cats(q) == {settle.CAT_DATA}, err or "")
    run, res, err = judged("m3h_positive_stop", end=40000, stop=40001)
    q = qres(res, "mass_balance") if res else {}
    chk("M3h (対照) 正の流量・停止の観測あり → SETTLED、ṁ_in・ṁ_out の相対の変化と停止の観測の検査を今どおり計算",
        err is None and q.get("verdict") == "SETTLED" and sum("(max − min)/mean" in c["check"] for c in q["checks"]) == 4,
        err or str([c["check"] for c in q.get("checks", [])]))
    sc = None
    try:
        sc = settle.scan(TMP / "m3a_out_zero", 2000)
    except Exception as e:  # noqa: BLE001
        sc = {"error": str(e)}
    chk("M3i --scan (格子点ごとの判定; オンラインの見張りと同じ judge) も出口 0 で例外にしない", "error" not in sc
        and sc.get("first_settled") is None, str(sc)[:160])


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


# ------------------------------------------------------------------ wt_nozzle_v2 の形 (§5.1 #5d)
GRID_BASE_COMMIT = "f210850c"     # residual.series・wall_heat_nodes.measure を足す前の settle.py (grid の既定のビット同一の基準)


def v2_criteria(**over):
    """wt_nozzle_v2 と同じ形: 窓 10,000 (min_points 6)・保持 0・残差 series raw・壁の熱流 measure range (閾値は v1 と同じ)。"""
    c = base_criteria()
    c["window"] = {"W": 10000, "min_points": 6}
    c["hold"] = 0
    for q in c["quantities"]:
        if q["kind"] == "residual":
            q["series"] = "raw"
        if q["kind"] == "wall_heat_nodes":
            q["measure"] = "range"
    c.update(over)
    return c


def braw(name, **kw):
    """forge と同じ形 (格子点 0 に観測なし・毎 step の残差の行) の v2 の run。既定の終端 22,000 = 最初に判定できる格子点。"""
    kw.setdefault("start_obs", False)
    kw.setdefault("rows_mode", "dense")
    kw.setdefault("end", 22000)
    kw.setdefault("criteria", v2_criteria())
    return build(name, **kw)


def rewrite_criteria(run, c):
    (run / settle.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))


def grid_of(c):
    """同じ基準で残差を grid (series を消す) にしたもの (対照)。"""
    c = json.loads(json.dumps(c))
    for q in c["quantities"]:
        q.pop("series", None)
    return c


def test_raw_residual():
    print("=== wt_nozzle_v2 の形: 残差の毎 step の系列 (series raw; §5.1 #5d、codex diagnose settle-v2-span10k) ===")
    run = braw("r01_const")
    r22 = settle.judge(run, 22000)
    q = qres(r22, "residual")
    v = q["values"]
    cols = {k: c for k, c in v["columns"].items() if c.get("status") == "judged"}
    chk("R1 一定の系列 (毎 step の行・区間 1 万) は最初に判定できる格子点 22,000 で SETTLED",
        r22["verdict"] == "SETTLED", f"{r22['verdict']} {[(x['name'], x['reasons'][:1]) for x in r22['quantities'] if x['verdict'] != 'SETTLED']}")
    chk("R1 区間は行 12,000〜21,999 (場の (12,000, 22,000]) と直前の区間は行 2,000〜11,999、各 10,000 行を全部使う",
        v.get("series") == "raw" and v["raw_rows"] == {"C": [12000, 21999], "P": [2000, 11999], "n_C": 10000, "n_P": 10000}
        and len(cols) == 6 and all(c["n_C"] == 10000 and c["n_P"] == 10000 and c["series"] == "raw" for c in cols.values()),
        str(v.get("raw_rows")))
    chk("R1 検査の detail に毎 step の系列と行の範囲", all("毎 step の系列" in (c.get("detail") or "") for c in q["checks"] if "跳ね" not in c["check"]))
    chk("R1 窓の格子点 6 点 (端点込み)", len(r22["window_steps"]["C"]) == 6 and r22["window_steps"]["C"][0] == 12000)
    r20 = settle.judge(run, 20000)
    q20 = qres(r20, "residual")
    chk("R1b 終端 20,000 (直前の区間が開始の直後から) は判定不能 (データ不足; 最初に判定できる終端は 22,000)",
        r20["verdict"] == "INDETERMINATE" and cats(q20) == {settle.CAT_DATA} and any("22000" in r["text"] and "開始" in r["text"] for r in q20["reasons"]),
        str(q20["reasons"][:1]))

    def spike(s):
        v = dict(RES_BASE)
        if s == 21999:
            v = {k: x * 1e8 for k, x in v.items()}
        return v
    q = qres(settle.judge(braw("r02_last_spike", resid_fn=spike), 22000), "residual")
    fc = failed_checks(q)
    chk("R2 末尾の 1 行だけ 10⁸ 倍 → NOT_SETTLED、跳ねの検査だけが落とす (毎 step の傾き・振幅・上昇は通る)",
        q["verdict"] == "NOT_SETTLED" and cats(q) == {settle.CAT_THRESHOLD} and fc and all("跳ね" in c for c in fc), str(fc[:3]))

    def burst(s):
        v = dict(RES_BASE)
        if 17100 <= s <= 17400:
            v = {k: x * 1e3 for k, x in v.items()}
        return v
    q = qres(settle.judge(braw("r03_burst", resid_fn=burst), 22000), "residual")
    fc = failed_checks(q)
    chk("R3 短い反復のバースト (行 17,100〜17,400 が 10³ 倍) → NOT_SETTLED、跳ねの検査が落とす", q["verdict"] == "NOT_SETTLED"
        and cats(q) == {settle.CAT_THRESHOLD} and any("跳ね" in c for c in fc), str(fc[:4]))

    def rise(s):
        v = dict(RES_BASE)
        if s >= 12000:
            v = {k: x * 1.12 for k, x in v.items()}
        return v
    q = qres(settle.judge(braw("r04_rise", resid_fn=rise), 22000), "residual")
    fc = failed_checks(q)
    chk("R4 持続的な上昇 (区間の全行が直前の区間より 12 % 高い水準で一定) → NOT_SETTLED、上昇の検査だけが落とす",
        q["verdict"] == "NOT_SETTLED" and cats(q) == {settle.CAT_THRESHOLD} and fc and all("上昇" in c for c in fc),
        f"{fc[:2]} rise {q['values']['columns']['rms_ro'].get('rise_dec')}")

    def partzero(s):
        v = dict(RES_BASE)
        if s > 15000:
            v["rms_roOmega"] = 0.0
        return v
    q = qres(settle.judge(braw("r05_partly_zero", resid_fn=partzero), 22000), "residual")
    chk("R5 残差の列の一部がゼロ → INDETERMINATE (データ不足)", q["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_DATA})

    # codex の再現 (settle-v2-span10k): 途中の 1 行の欠け・2,000 ごとの行だけ。grid では SETTLED、raw では判定不能
    run = braw("r06_one_row_missing", skip_resid_rows=(15555,))
    q = qres(settle.judge(run, 22000), "residual")
    chk("R6 途中の 1 行 (15,555) の欠け → raw は INDETERMINATE (データ不足、欠けの行を書く)", q["verdict"] == "INDETERMINATE"
        and cats(q) == {settle.CAT_DATA} and any("欠け 1 行" in r["text"] and "15555" in r["text"] for r in q["reasons"]), str(q["reasons"][:1]))
    rewrite_criteria(run, grid_of(v2_criteria()))
    chk("R6 (対照) 同じ履歴を grid で判定すると SETTLED (格子点の行しか見ない — codex の指摘の再現)",
        qres(settle.judge(run, 22000), "residual")["verdict"] == "SETTLED")
    run = braw("r07_grid_rows_only", rows_mode="grid_only")
    q = qres(settle.judge(run, 22000), "residual")
    chk("R7 2,000 ごとの行 (格子点の行 s−1) だけの履歴 → raw は INDETERMINATE (区間・直前の区間とも 9,995 行の欠け)",
        q["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_DATA}
        and sum("欠け 9995 行" in r["text"] for r in q["reasons"]) == 2, str([r["text"][:60] for r in q["reasons"]]))
    rewrite_criteria(run, grid_of(v2_criteria()))
    chk("R7 (対照) 同じ履歴を grid で判定すると SETTLED", qres(settle.judge(run, 22000), "residual")["verdict"] == "SETTLED")

    # 出力の境界の step: 区間の端の行が 1 行でも欠ければ判定不能、区間の外の行の欠けは影響しない
    for row, what in ((21999, "区間の最後 (T−1)"), (12000, "区間の最初 (T−W)"), (11999, "直前の区間の最後 (T−W−1)"),
                      (2000, "直前の区間の最初 (T−2W)")):
        q = qres(settle.judge(braw(f"r08_edge_{row}", end=24000, skip_resid_rows=(row,)), 22000), "residual")
        chk(f"R8 出力の境界: 行 {row} ({what}) の欠け → INDETERMINATE", q["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_DATA}
            and any(str(row) in r["text"] for r in q["reasons"]), str(q["reasons"][:1]))
    r = settle.judge(braw("r08_outside", end=24000, skip_resid_rows=(1999, 22000)), 22000)
    chk("R8 出力の境界: 区間の外の行 1,999 (T−2W−1) と 22,000 (T) の欠けは影響しない (SETTLED; 区間は行 2,000〜21,999 ちょうど)",
        r["verdict"] == "SETTLED", f"{r['verdict']} {[x['reasons'][:1] for x in r['quantities'] if x['verdict'] != 'SETTLED']}")

    q = qres(settle.judge(braw("r09_inner_rows", inner_rows=True), 22000), "residual")
    chk("R9 内反復の行 (inner_begin・inner_end、10⁸ 倍) が混じっても outer_end だけを使う (SETTLED・区間 10,000 行)",
        q["verdict"] == "SETTLED" and q["values"]["raw_rows"]["n_C"] == 10000, str(q["reasons"][:1]))
    q = qres(settle.judge(braw("r10_dup_row", dup_resid_rows=(15000,)), 22000), "residual")
    chk("R10 区間の行の重複 (15,000 を 2 回) → INDETERMINATE (重複)", q["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_DATA}
        and any("重複" in r["text"] for r in q["reasons"]), str(q["reasons"][:2]))

    def nanres(s):
        v = dict(RES_BASE)
        if s == 16000:
            v["rms_roe"] = float("nan")
        return v
    q = qres(settle.judge(braw("r11_nan", resid_fn=nanres), 22000), "residual")
    chk("R11 残差に NaN → DIVERGED (非有限)", q["verdict"] == "DIVERGED" and cats(q) == {settle.CAT_NONFINITE})
    r = settle.judge(braw("r12_stop", stop=22001))
    q = qres(r, "residual")
    chk("R12 一定の系列を非格子の 22,001 で止めた → SETTLED (終端 22,000、停止の観測の別検査あり)",
        r["verdict"] == "SETTLED" and r["terminal_step"] == 22000 and r["stop_step"] == 22001
        and any("停止の観測" in c["check"] for c in q["checks"]), f"{r['verdict']} {[x['reasons'][:1] for x in r['quantities'] if x['verdict'] != 'SETTLED']}")
    for key, val, word in (("series", "every", "series"), ("measure", "max", "measure")):
        c = v2_criteria()
        for qq in c["quantities"]:
            if key in qq:
                qq[key] = val
        try:
            settle.validate_criteria(c)
            ok = False
        except settle.CriteriaError as e:
            ok = word in str(e)
        chk(f"R13 {key}: '{val}' (知らない値) → CriteriaError", ok)


def test_wall_range():
    print("=== wt_nozzle_v2 の形: 壁の熱流の区間の幅 (measure range; codex diagnose settle-v2-span10k) ===")
    xm, rm = grid_xr()
    coord = np.stack([xm[:, -1], rm[:, -1], np.zeros(NI)], axis=1)
    conne = []
    for i in range(NI - 1):
        conne += [2, 2, i, i + 1]
    A, _ = settle.wall_node_area(coord, np.asarray(conne), True)
    k = 10
    seq = {12000: 0.25, 14000: 1.75}                 # 区間の格子点 12,000〜22,000 の値 [0.25, 1.75, 1, 1, 1, 1]

    def Qc(s):
        Q = np.ones(NI)
        Q[k] = seq.get(s, 1.0)
        return Q
    c = v2_criteria()
    for q in c["quantities"]:
        if q["kind"] == "wall_heat_nodes":
            q["thresholds"] = {"q_tol": 1.001 / float(A[k]), "eps_rel": 0.0, "provisional": True}    # 節点 10 の許容 1.001
    run = braw("w01_codex_example", Q_fn=Qc, criteria=c)
    r = settle.judge(run, 22000)
    q = qres(r, "wall_heat_nodes")
    w = q["values"].get("worst_window") or {}
    chk("W1 codex の例 [0.25, 1.75, 1, 1, 1, 1]・許容 1.001: range は幅 1.5 で NOT_SETTLED (節点 10、比 1.5/1.001)",
        r["verdict"] == "NOT_SETTLED" and q["verdict"] == "NOT_SETTLED" and cats(q) == {settle.CAT_THRESHOLD} and w.get("node") == k
        and abs(w["span"] - 1.5) < 1e-12 and abs(w["ratio"] - 1.5 / 1.001) < 1e-9 and q["values"].get("measure") == "range",
        f"{r['verdict']} {w}")
    chk("W1 range の記録: 終端からの差の最大 (v1 の量) も残す (比 0.749)", abs((q["values"].get("worst_terminal_diff") or {}).get("ratio", 0) - 0.75 / 1.001) < 1e-9)
    for qq in c["quantities"]:
        qq.pop("measure", None)
    rewrite_criteria(run, c)
    r = settle.judge(run, 22000)
    q = qres(r, "wall_heat_nodes")
    chk("W1 (対照) 同じ系列を既定 (terminal_diff、v1 の判定) で判定すると今どおり SETTLED (比 0.75/1.001 = 0.74925)",
        r["verdict"] == "SETTLED" and abs(q["values"]["worst_window"]["ratio"] - 0.75 / 1.001) < 1e-9 and "measure" not in q["values"],
        f"{r['verdict']} {q['values']['worst_window'].get('ratio')}")
    r = settle.judge(braw("w02_const_stop", stop=22001))
    q = qres(r, "wall_heat_nodes")
    names = [x["check"] for x in q["checks"]]
    chk("W2 range: 一定の壁を 22,001 で止めた → SETTLED、停止の観測の終端からの差と加えた幅を検査",
        r["verdict"] == "SETTLED" and any("停止の観測を加えた節点の熱量の幅" in n for n in names)
        and any("停止の観測 (step 22001)" in n for n in names), str(names))

    def Qs(s):
        Q = 100.0 + np.arange(NI, dtype=float)
        if s == 22001:
            Q[k] += 1.2 * (1000.0 * A[k] + 1e-3 * Q[k])
        return Q
    q = qres(settle.judge(braw("w03_stop_move", stop=22001, Q_fn=Qs)), "wall_heat_nodes")
    fc = failed_checks(q)
    chk("W3 range: 停止の観測だけ節点 10 が許容の 1.2 倍動く → NOT_SETTLED (停止の観測の検査だけ)", q["verdict"] == "NOT_SETTLED"
        and fc and all("停止の観測" in x for x in fc), str(fc))


def _load_base_settle():
    """grid の既定のビット同一の基準: GRID_BASE_COMMIT の settle.py を git から読む (読めなければ None と理由)。"""
    import importlib.util
    import subprocess
    try:
        src = subprocess.run(["git", "-C", str(HERE), "show", f"{GRID_BASE_COMMIT}:solver_density_cuda/tools/settle.py"],
                             capture_output=True, text=True, check=True).stdout
    except Exception as e:  # noqa: BLE001
        return None, f"git show {GRID_BASE_COMMIT}:solver_density_cuda/tools/settle.py が読めない ({type(e).__name__})"
    p = TMP / f"settle_base_{GRID_BASE_COMMIT}.py"
    p.write_text(src)
    spec = importlib.util.spec_from_file_location(f"settle_base_{GRID_BASE_COMMIT}", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m, None


def _strip(res):
    d = json.loads(json.dumps(settle._jsonable(res), ensure_ascii=False))
    for k in ("created", "tool_sha256"):
        d.pop(k, None)
    return d


def test_grid_bit_identity():
    print(f"=== 既定 (series grid・measure terminal_diff) は変更前 ({GRID_BASE_COMMIT}) の settle.py とビット同一 ===")
    base, why = _load_base_settle()
    if base is None:
        chk("G0 変更前の settle.py を読めた (読めなければ判定不能 — 合格にしない)", False, why)
        return

    def spike(s):
        v = dict(RES_BASE)
        if s == 39999:
            v = {k: x * 1e8 for k, x in v.items()}
        return v

    def burst(s):
        v = dict(RES_BASE)
        if 33100 <= s <= 33400:
            v = {k: x * 1e3 for k, x in v.items()}
        return v

    def rise(s):
        v = dict(RES_BASE)
        if s >= 19999:
            v = {k: x * 1.12 for k, x in v.items()}
        return v

    def partzero(s):
        v = dict(RES_BASE)
        if s > 30000:
            v["rms_roOmega"] = 0.0
        return v

    def swap(s):
        Q = 100.0 + np.arange(NI, dtype=float)
        d = 200.0 * math.sin(s / 7000.0)
        Q[10] += d
        Q[20] -= d
        return Q

    def drift(s):
        f = 10 ** (-2.0 * min(s, 30000) / 30000)
        return {k: x * f for k, x in RES_BASE.items()}
    cases = [("g01_const", {}, [None, 30000]), ("g02_stop", {"end": 80000, "stop": 80001}, [None]),
             ("g03_forge_like", {"end": 42000, "start_obs": False}, [40000, 42000]), ("g04_spike", {"resid_fn": spike}, [None]),
             ("g05_burst", {"resid_fn": burst, "extra_resid_rows": range(33100, 33401, 50)}, [None]), ("g06_rise", {"resid_fn": rise}, [None]),
             ("g07_gap", {"skip": tuple(range(26000, 34001, DT))}, [None]), ("g08_partzero", {"resid_fn": partzero}, [None]),
             ("g09_wall_swap", {"Q_fn": swap}, [None]), ("g10_wall_cross", {"Q_fn": lambda s: np.full(NI, 1.0 if (s // DT) % 2 else -1.0)}, [None]),
             ("g11_dense", {"rows_mode": "dense", "end": 42000, "start_obs": False}, [42000]),
             ("g12_grid_only", {"rows_mode": "grid_only", "end": 42000, "start_obs": False}, [42000]),
             ("g13_drift", {"resid_fn": drift, "end": 90000}, [None, 60000])]
    bad = []
    n = 0
    for name, kw, ends in cases:
        run = build(name, **kw)
        for e in ends:
            a = _strip(settle.judge(run, e))
            b = _strip(base.judge(run, e))
            n += 1
            if a != b:
                bad.append(f"{name} end {e}")
    chk(f"G1 既定の判定の結果 (全量の verdict・検査・値・理由・基準の中身) が変更前とビット同一 ({n} 件; 生成の時刻・道具の sha256 を除く)",
        not bad, str(bad))
    run = build("g13_drift_scan", resid_fn=drift, end=90000)
    chk("G2 --scan の結果も変更前と同一", settle.scan(run, DT) == base.scan(run, DT))


# ------------------------------------------------------------------ wt_nozzle_v3 の形 (残差を外し非有限の見張りだけ; §4.4・§5.1 #8)
PRE8_COMMIT = "1de17d21"     # #7c・#8 の前の settle.py (v1・v2 の形の基準の判定のビット同一の基準)


def v3_criteria(**over):
    """wt_nozzle_v3 と同じ形: v2 の形から残差の量を外し、残差の非有限の見張り (residual_finite) を入れて excluded に理由を書く。"""
    c = v2_criteria()
    c["quantities"] = [q for q in c["quantities"] if q["kind"] != "residual"]
    c["quantities"].insert(0, {"kind": "residual_finite", "name": "residual_finite", "thresholds": {}})
    c["excluded"] = [{"kind": "residual", "name": "residual", "reason": "試験: 残差を止める判定から外した (wt_nozzle_v3 の形)"}]
    c.update(over)
    return c


def test_v3():
    """plan §4.4「残差を止める判定から外す・ビルドは float」(2026-10-16 ユーザ決定)・§5.1 #8。"""
    print("=== wt_nozzle_v3 の形: 残差を外し非有限の見張りだけ (§4.4・§5.1 #8) ===")
    osc = lambda s: {k: x * (10.0 ** 0.6 if (s // 7) % 2 else 1.0) for k, x in RES_BASE.items()}   # 振幅 0.6 桁 (float の ρ は 0.585 桁)
    run = braw("v3a_osc", resid_fn=osc, criteria=v3_criteria())
    res = settle.judge(run)
    q = qres(res, "residual_finite")
    chk("V3a 残差の振幅 0.6 桁でも、流量と目的の量が整えば SETTLED (残差は判定しない; 見張りは SETTLED)",
        res["verdict"] == "SETTLED" and q["verdict"] == "SETTLED" and not any(x["kind"] == "residual" for x in res["quantities"])
        and q["values"]["role"].startswith("非有限の見張り") and q["checks"] == [],
        f"{res['verdict']} {[(x['name'], x['verdict']) for x in res['quantities']]}")
    rewrite_criteria(run, v2_criteria())
    res2 = settle.judge(run)
    chk("V3a (対照) 同じ残差を v2 の形の基準で → NOT_SETTLED (残差の振幅)", res2["verdict"] == "NOT_SETTLED"
        and qres(res2, "residual")["verdict"] == "NOT_SETTLED", str(failed_checks(qres(res2, "residual")))[:160])
    nan = lambda s: {k: (float("nan") if s >= 21000 else x) for k, x in RES_BASE.items()}
    run = braw("v3b_nan", resid_fn=nan, criteria=v3_criteria())
    res = settle.judge(run)
    q = qres(res, "residual_finite")
    chk("V3b 残差の非有限 (行 21,000 から) → 見張りが DIVERGED・全体 DIVERGED (分類 nonfinite、最初の行の step)",
        res["verdict"] == "DIVERGED" and q["verdict"] == "DIVERGED" and cats(q) == {settle.CAT_NONFINITE}
        and q["values"]["nonfinite"]["first_step"]["rms_ro"] == 21000, str(q["reasons"])[:200])
    rc, out = run_cli(run)
    chk("V3b check_settled.py は終了コード 2 (DIVERGED)", rc == 2 and "VERDICT: DIVERGED" in out, out[-120:])
    run = braw("v3c_stop_nan", end=22000, stop=22001, criteria=v3_criteria(),
               resid_fn=lambda s: {k: (float("inf") if s == 22000 else x) for k, x in RES_BASE.items()})
    res = settle.judge(run)
    chk("V3c 停止の観測 22,001 の直前の行 (22,000) だけ非有限 → DIVERGED (見張りは停止の step − 1 まで)",
        qres(res, "residual_finite")["verdict"] == "DIVERGED", str(qres(res, "residual_finite")["reasons"])[:160])
    run = braw("v3d_stop_ok", end=22000, stop=22001, criteria=v3_criteria(), resid_fn=osc)
    res = settle.judge(run)
    chk("V3d 停止の観測 22,001 (有限) → SETTLED (終端 22,000)", res["verdict"] == "SETTLED" and res["terminal_step"] == 22000
        and res["stop_step"] == 22001, res["verdict"])
    run = braw("v3e_lastrow", criteria=v3_criteria(), skip_resid_rows=(21999,))
    q = qres(settle.judge(run), "residual_finite")
    chk("V3e 区間の最後の行 (21,999) が無い → 見張りは判定不能 (終端まで届かない)", q["verdict"] == "INDETERMINATE"
        and any("最後の行" in r["text"] for r in q["reasons"]), str(q["reasons"])[:200])
    (run / "residual_history.csv").unlink()
    q = qres(settle.judge(run), "residual_finite")
    chk("V3f 残差の履歴が無い → 見張りは判定不能 (データ不足)", q["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_DATA})
    c = v3_criteria()
    c["quantities"] = [x for x in c["quantities"] if x["kind"] != "residual_finite"]
    run = braw("v3g_nowatch", criteria=c)
    lc = settle.load_criteria(run)
    rc, out = run_cli(run)
    chk("V3g 残差も見張りも無い基準 → INDETERMINATE (評価対象が未設定; 発散を見張れない)・終了コード 3",
        lc["status"] == "unset" and "residual_finite" in lc["reason"] and rc == 3 and "評価対象が未設定" in out, lc.get("reason", "")[:160])
    bad = []
    # 期待の変更 (§5.1 #10、result 段レビュー 2 回目 M2): 以前は「知らないキー required」を負例にしていた。#10 で residual_finite は
    # 必要な列 `required` を受け付けるようにしたので、知らないキーの負例は inactive_ok (residual_finite には無い) に替え、required の形の
    # 誤り (空・rms_dq_*・rms_ で始まらない) を負例に足す。
    for label, mut in (("2 つ", lambda c: c["quantities"].append({"kind": "residual_finite", "name": "rf2", "thresholds": {}})),
                       ("知らないキー inactive_ok", lambda c: c["quantities"][0].update(inactive_ok=["rms_roUz"])),
                       ("required が空", lambda c: c["quantities"][0].update(required=[])),
                       ("required に rms_dq_", lambda c: c["quantities"][0].update(required=["rms_ro", "rms_dq_ro"])),
                       ("required に rms_ でない名前", lambda c: c["quantities"][0].update(required=["ro"])),
                       ("閾値のキー slope_dec", lambda c: c["quantities"][0].update(thresholds={"slope_dec": 0.1})),
                       ("thresholds が無い", lambda c: c["quantities"][0].pop("thresholds"))):
        c = v3_criteria()
        mut(c)
        rewrite_criteria(run, c)
        try:
            settle.load_criteria(run)
            bad.append(label)
        except settle.CriteriaError:
            pass
    chk("V3h residual_finite の誤用 (2 つ・知らないキー・required の形の誤り・閾値のキー・thresholds が無い) は CriteriaError", not bad, str(bad))
    c = v3_criteria()
    c["quantities"][0]["thresholds"] = {"provisional": True}
    rewrite_criteria(run, c)
    chk("V3h thresholds は {} か provisional だけ", settle.load_criteria(run)["status"] == "ok")
    run = braw("v3i_scan", end=30000, criteria=v3_criteria(), resid_fn=osc)
    sc = settle.scan(run, DT)
    import settle_virtual_stop as VS
    fj = VS.first_judgeable(settle.validate_criteria(v3_criteria()))
    chk("V3i 判定の部品 (--scan) は v3 の形で区間 C だけを使うので 12,000 から SETTLED (直前の区間は残差の上昇の検査だけが使う)、"
        "停止規則・再生の最初の判定は量によらず開始 + 2W + 間隔 = 22,000 のまま (v2 と同じ; 早めるかは plan の判断)",
        sc["first_settled"] == 12000 and fj == 22000, f"scan {sc['first_settled']} first_judgeable {fj}")


PRE_REVIEW2_COMMIT = "b64dbfdd"    # result 段レビュー 2 回目 (§5.1 #9・#10) の修正の前 (負例が修正前に落ちないこと = 検出力の確かめ)


def _load_at(commit: str, relpath: str, modname: str):
    """commit の relpath を git から読んでモジュールにする (読めなければ None と理由)。check_settled.py は同じディレクトリの settle を
    import するので、今の settle (sys.modules) と組になる (#9 の修正は check_settled.py の分岐だけなので、修正前の分岐を確かめられる)。"""
    import importlib.util
    import subprocess
    try:
        src = subprocess.run(["git", "-C", str(HERE), "show", f"{commit}:{relpath}"], capture_output=True, text=True, check=True).stdout
    except Exception as e:  # noqa: BLE001
        return None, f"git show {commit}:{relpath} が読めない ({type(e).__name__})"
    p = TMP / f"{modname}.py"
    p.write_text(src)
    spec = importlib.util.spec_from_file_location(modname, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m, None


def _rewrite_residual(run, header: list, rows: list):
    """残差の履歴を書き直す (rows は step,inner,phase,値… の行の文字列のリスト)。"""
    with open(run / "residual_history.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(r + "\n")


def _outer_end_rows(run):
    lines = (run / "residual_history.csv").read_text().splitlines()
    return lines[0].split(","), [ln for ln in lines[1:] if ln.split(",")[2] == "outer_end"]


def test_v3_watch_coverage():
    """result 段レビュー 2 回目 M2 (plan §5.1 #10): v3 の非有限の見張り (residual_finite) は、在る数値だけを見て最終行があれば合格に
    していた。閾値は戻さず、区間の行の被覆 (毎 step の行が重複なく全部・逆行なし) と必要な列を必須にした。負例ごとに修正前
    (PRE_REVIEW2_COMMIT の settle.py) が SETTLED だったこと (検出力) も確かめる。"""
    print("=== v3 の見張りの被覆と必要な列 (result 段レビュー 2 回目 M2、§5.1 #10) ===")
    old, why = _load_at(PRE_REVIEW2_COMMIT, "solver_density_cuda/tools/settle.py", f"settle_pre_review2_{PRE_REVIEW2_COMMIT}")
    if old is None:
        chk(f"F0 修正前 ({PRE_REVIEW2_COMMIT}) の settle.py を読めた (読めなければ検出力は判定不能 — 合格にしない)", False, why)

    def both(run):
        new = settle.judge(run)
        prev = old.judge(run) if old is not None else None
        return new, qres(new, "residual_finite"), (prev["verdict"] if prev is not None else None)

    # 対照 (正例): forge と同じ形 (毎 step の行・全列) は SETTLED
    run = braw("w10_ok", criteria=v3_criteria())
    res, q, _ = both(run)
    chk("F1 対照: 毎 step の行・全列 → SETTLED (見張りの必要な列は既定の 4 列、出所を記録)", res["verdict"] == "SETTLED"
        and q["verdict"] == "SETTLED" and q["values"]["required"] == list(settle.RESIDUAL_FINITE_DEFAULT_REQUIRED)
        and q["values"]["required_source"].startswith("既定"), str(q["reasons"])[:200])
    hdr, rows = _outer_end_rows(run)
    # ① 残差を最終 1 行だけ (codex の再現)
    r1 = braw("w11_lastrow_only", criteria=v3_criteria())
    _rewrite_residual(r1, hdr, rows[-1:])
    res, q, pv = both(r1)
    rc, out = run_cli(r1)
    chk("F2 ① 残差が最終 1 行だけ → 見張りは判定不能 (区間の行の欠け)・全体 INDETERMINATE・CLI 終了コード 3",
        q["verdict"] == "INDETERMINATE" and cats(q) == {settle.CAT_DATA} and any("行欠ける" in r["text"] for r in q["reasons"])
        and res["verdict"] == "INDETERMINATE" and rc == 3, str(q["reasons"])[:200])
    chk("F2 検出力: 修正前は ① で SETTLED", pv == "SETTLED", str(pv))
    # ② rms_ro の列だけ (行は全部)
    r2 = braw("w12_ro_only", criteria=v3_criteria())
    i_ro = hdr.index("rms_ro")
    _rewrite_residual(r2, ["step", "inner", "phase", "rms_ro"], [",".join(ln.split(",")[:3] + [ln.split(",")[i_ro]]) for ln in rows])
    res, q, pv = both(r2)
    chk("F3 ② 残差の列が rms_ro だけ → 見張りは判定不能 (必要な列の欠け: rms_roUx・rms_roUy・rms_roe)",
        q["verdict"] == "INDETERMINATE" and any("必要な残差の列が無い" in r["text"] and "rms_roUx" in r["text"] for r in q["reasons"])
        and res["verdict"] == "INDETERMINATE", str(q["reasons"])[:200])
    chk("F3 検出力: 修正前は ② で SETTLED", pv == "SETTLED", str(pv))
    # ③ 残差の列が 1 つも無い step,inner,phase の 1 行
    r3 = braw("w13_no_cols", criteria=v3_criteria())
    _rewrite_residual(r3, ["step", "inner", "phase"], [",".join(rows[-1].split(",")[:3])])
    res, q, pv = both(r3)
    chk("F4 ③ 残差の列の無い step,inner,phase の 1 行 → 見張りは判定不能 (行の欠けと列の欠けの両方)",
        q["verdict"] == "INDETERMINATE" and any("行欠ける" in r["text"] for r in q["reasons"])
        and any("必要な残差の列が無い" in r["text"] for r in q["reasons"]), str(q["reasons"])[:200])
    chk("F4 検出力: 修正前は ③ で SETTLED", pv == "SETTLED", str(pv))
    # 区間の途中の 1 行の欠け・重複・逆行 (最終行はある)
    r4 = braw("w14_gap", criteria=v3_criteria(), skip_resid_rows=(15000,))
    res, q, pv = both(r4)
    chk("F5 区間の途中の 1 行 (15,000) の欠け → 判定不能 (欠けた step を理由に)", q["verdict"] == "INDETERMINATE"
        and any("1 行欠ける" in r["text"] and "15000" in r["text"] for r in q["reasons"]), str(q["reasons"])[:200])
    chk("F5 検出力: 修正前は途中の欠けで SETTLED", pv == "SETTLED", str(pv))
    r5 = braw("w15_dup", criteria=v3_criteria(), dup_resid_rows=(15000,))
    res, q, pv = both(r5)
    chk("F6 行の重複 (15,000 を 2 回) → 判定不能", q["verdict"] == "INDETERMINATE" and any("重複" in r["text"] for r in q["reasons"]),
        str(q["reasons"])[:200])
    chk("F6 検出力: 修正前は重複で SETTLED", pv == "SETTLED", str(pv))
    r6 = braw("w16_nonmono", criteria=v3_criteria())
    i = next(k for k, ln in enumerate(rows) if ln.startswith("15000,"))
    swapped = rows[:i] + [rows[i + 1], rows[i]] + rows[i + 2:]
    _rewrite_residual(r6, hdr, swapped)
    res, q, pv = both(r6)
    chk("F7 行の逆行 (15,001 の後に 15,000) → 判定不能", q["verdict"] == "INDETERMINATE" and any("単調でない" in r["text"] for r in q["reasons"]),
        str(q["reasons"])[:200])
    chk("F7 検出力: 修正前は逆行で SETTLED", pv == "SETTLED", str(pv))
    # 基準の required (乱流の列を書いたとき): 列が欠ければ判定不能、全部あれば SETTLED
    c = v3_criteria()
    c["quantities"][0]["required"] = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega"]
    r7 = braw("w17_rans_missing", criteria=c, resid_cols=[k for k in RES_COLS if k != "rms_roK"])
    q = qres(settle.judge(r7), "residual_finite")
    chk("F8 基準の required に rms_roK・rms_roOmega (RANS) を書き、履歴に rms_roK が無い → 判定不能 (列の欠け)",
        q["verdict"] == "INDETERMINATE" and any("rms_roK" in r["text"] for r in q["reasons"])
        and q["values"]["required_source"] == "基準の required", str(q["reasons"])[:200])
    r8 = braw("w18_rans_ok", criteria=c)
    q = qres(settle.judge(r8), "residual_finite")
    chk("F9 同じ基準で全列がある → SETTLED", q["verdict"] == "SETTLED" and q["values"]["required"] == c["quantities"][0]["required"],
        str(q["reasons"])[:200])
    # 非有限は今どおり DIVERGED (欠けがあっても DIVERGED を優先)
    r9 = braw("w19_nan_gap", criteria=v3_criteria(), skip_resid_rows=(15000,),
              resid_fn=lambda s: {k: (float("nan") if s == 18000 else x) for k, x in RES_BASE.items()})
    res, q, _ = both(r9)
    chk("F10 非有限 (18,000) と途中の欠け (15,000) → DIVERGED (非有限が優先、今どおり)・CLI 終了コード 2",
        q["verdict"] == "DIVERGED" and res["verdict"] == "DIVERGED" and run_cli(r9)[0] == 2, str(q["reasons"])[:200])


def test_v1_v2_identity_pre8():
    """v1・v2 の形の基準の判定が #7c・#8 の前 (PRE8_COMMIT) の settle.py とビット同一 (正常な流量; #7c の M3 は流量 ≤ 0 だけを変える)。"""
    print(f"=== v1・v2 の形の基準の判定は #7c・#8 の前 ({PRE8_COMMIT}) とビット同一 ===")
    global GRID_BASE_COMMIT
    saved = GRID_BASE_COMMIT
    GRID_BASE_COMMIT = PRE8_COMMIT
    try:
        base, why = _load_base_settle()
    finally:
        GRID_BASE_COMMIT = saved
    if base is None:
        chk("I0 #7c・#8 の前の settle.py を読めた (読めなければ判定不能 — 合格にしない)", False, why)
        return
    spike = lambda s: {k: x * (1e8 if s == 21999 else 1.0) for k, x in RES_BASE.items()}
    cases = [("i01_v1_const", dict(), [None, 30000]), ("i02_v1_stop", dict(end=80000, stop=80001), [None]),
             ("i03_v1_mass_drift", dict(mdot_fn=lambda s: (10.0, 10.0 * (1 - 1e-3 * s / 40000))), [None]),
             ("i04_v2_const", dict(criteria=v2_criteria(), start_obs=False, rows_mode="dense", end=30000), [22000, 30000]),
             ("i05_v2_spike", dict(criteria=v2_criteria(), start_obs=False, rows_mode="dense", end=22000, resid_fn=spike), [None]),
             ("i06_v2_stop", dict(criteria=v2_criteria(), start_obs=False, rows_mode="dense", end=22000, stop=22001), [None]),
             ("i07_v2_wall", dict(criteria=v2_criteria(), start_obs=False, rows_mode="dense", end=22000,
                                  Q_fn=lambda s: 100.0 + np.arange(NI, dtype=float) + (5.0 if s == 20000 else 0.0)), [None])]
    bad, n = [], 0
    for name, kw, ends in cases:
        run = build(name, **kw)
        for e in ends:
            n += 1
            if _strip(settle.judge(run, e)) != _strip(base.judge(run, e)):
                bad.append(f"{name} end {e}")
    chk(f"I1 v1・v2 の形の基準の判定の結果 (全量の verdict・検査・値・理由・基準の中身) が {PRE8_COMMIT} とビット同一 ({n} 件)", not bad, str(bad))


def main():
    print(f"numpy {np.__version__}, python {sys.version.split()[0]}, 作業 {TMP}")
    test_area()
    test_positive()
    test_negative_v2()
    test_negative_extra()
    test_axis_weight()
    test_verdict_file()
    test_select_bind()
    test_mass_denominator()
    test_v3()
    test_v3_watch_coverage()
    test_v1_v2_identity_pre8()
    test_raw_residual()
    test_wall_range()
    test_grid_bit_identity()
    shutil.rmtree(TMP, ignore_errors=True)
    print(f"\n{NCHK - len(FAILS)}/{NCHK} OK")
    if FAILS:
        print("FAIL:", *FAILS, sep="\n  ")
    print(f"VERDICT: {'PASS' if not FAILS else 'FAIL'}")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
