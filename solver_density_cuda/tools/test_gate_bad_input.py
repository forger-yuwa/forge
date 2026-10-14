"""判定ツールが **不正な入力を合格にしない** ことの回帰試験。

2026-09-19 の codex plan レビュー (tooling-convergence-and-wall-resolution-gates, NO-GO/M9) が
現行コードで再現した誤合格をそのまま試験にする。判定ロジックを触るときは必ずこれを通すこと。

    python3 solver_density_cuda/tools/test_gate_bad_input.py

対象と、レビュー時点で**合格してしまっていた**入力:

| ツール | 入力 | 当時の結果 |
| --- | --- | --- |
| `check_convergence.analyze` | `step,phase` のみ (残差列なし) | `ok=True` |
| 〃 | `rms_ro` だけ (他の保存量なし) | `ok=True` |
| 〃 | 19 点が 1 で**最終 1 点だけ 0** | `drop=inf`, `ok=True` |
| `check_quasisteady.classify` | `[1,1,1,1,1,NaN]` | `STEADY` |
| 〃 `classify_series` | step が全点 NaN・値は一定 | `STEADY` |
| `check_mesh_quality` | 体積ゼロの四面体 | AR 1.732 / skew 0.500 = 通常合格 |
| 〃 | 1001 セル中 1 セルが NaN 座標 | `SOFT-PASS`, exit 0 |
| `check_wall_resolution` (2026-10-15, §5.1 #36q) | `--groups a,b` で b の壁ダンプが無い | a だけで判定して PASS |
| 〃 | `--step S` で b の S の壁ダンプが無い (古い step はある) | b は最新の step に戻って PASS |
| 〃 | `--step S` の壁ダンプがどの壁にも無い (古い step だけ) | 最新の step で PASS |

`check_wall_resolution` の合成の run は `make_wall_run` (2 枚の no-slip 壁 + 出口の出力、node 方式・軸対称; forge なし)。
受け渡しの関門 (`deltastar_loop.wall_resolution_verdict`) の受入は同じ合成の run で
`case/45.isobutane_m6_d155/test_h370_launch.py` の O で見る。

`mesh_provenance.py check` (2026-10-15 新設、plan tooling-generator-defaults-and-launch-gate §5.1 #7) も記録の無い・壊れた・
別メッシュの記録を合格にしないことをここで見る。

`restart_field.py` の同一メッシュの検査 (2026-10-15、同 plan §5.1 #20): 座標を ~1e−8 変えた別の入力 (float32 では同じ)・
接続だけ違う入力・接続の無い格子・格子を決められない res が、修正前は**通っていた** (終了 0)。

`launch_gate.py` (投入の関門、2026-10-15、同 plan §5.1 #10): それまで `run_case.sh` は設定・メッシュを照合せずに forge を起動していた
(未知のキーは check_solver_config で WARN・終了 0、AR 6000 の層や bcondConfig の無い run も起動)。検査できない入力・メッシュの FAIL・
禁止キー (`FORGE_LAUNCH_GATE=warn` でも)・`FORGE_LAUNCH_GATE` の綴り違いを「起動してよい」にしない。

入力の完全性 (2026-10-15、plan tooling-convergence-and-wall-resolution-gates §4.6-3、codex plan 段 judge M3 の再現):

| ツール | 入力 | 当時の結果 |
| --- | --- | --- |
| `check_convergence --segment` | 同じ区間の 2 段のうち本段の履歴が無い | 前段だけを連結して PASS |
| 〃 | 区間の段の履歴が空 | その段を飛ばして判定 |
| `check_quasisteady.classify_series` | step = [1,1,1,1,1,1]・値が一定 | `STEADY` |
| 〃 | step が逆行する | 判定した |
| `check_quasisteady` (run) | `--quantity asym` だけ・上下対称でないメッシュ | 評価量ゼロで `STEADY`・終了 0 |
| 〃 `--series-csv` | `--series-cols` が空 | `STEADY`・終了 0 |
"""
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_convergence as cc          # noqa: E402
import check_quasisteady as cq          # noqa: E402

ok = True


def chk(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print("  [%s] %-58s %s" % ("OK " if cond else "NG ", name, detail))


def okof(res):
    """analyze() の戻り (laststep, report, ok, ...) から「合格にしていないか」を判定用に返す。"""
    if res is None:
        return True, "analyze -> None (判定不能)"
    ok_ = res[2]
    return (not ok_), "ok=%s" % ok_


def completeness_cases():
    """入力の完全性 (2026-10-15、plan tooling-convergence-and-wall-resolution-gates §4.6-3): 欠けた段・step の重複と逆行・
    評価量ゼロを合格にしない。"""
    import shutil
    import h5py
    import stage_manifest as sm
    print("=== 入力の完全性: 欠けた段・step が進まない系列・評価量ゼロ (§4.6-3) ===")
    cc_tool = os.path.join(HERE, "check_convergence.py")
    cq_tool = os.path.join(HERE, "check_quasisteady.py")
    cols = "step,phase,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n"
    cfg = 'turbulence: {model: "sst", wallTreatmentSST: 0}\nspace: {convMethod: 1, limiter: 2}\n'

    def hist(path, n, d0, d1):
        with open(path, "w") as f:
            f.write(cols)
            for i in range(n):
                v = 10.0 ** (d0 + (d1 - d0) * i / (n - 1))
                f.write("%d,outer_end,%g,%g,%g,%g,%g\n" % (i, v, v, v, v, v))

    def seg_run(main_hist):
        d = tempfile.mkdtemp(prefix="seg_gate_")
        m = sm.StageManifest(d)
        m.add("S1", cfg, "wall: {}\n", history="residual_history_S1.csv")
        m.add("main", cfg, "wall: {}\n", history="residual_history_main.csv")
        m.write()
        hist(os.path.join(d, "residual_history_S1.csv"), 100, 0.0, -2.0)
        if main_hist == "full":
            hist(os.path.join(d, "residual_history_main.csv"), 100, -2.0, -5.0)
        elif main_hist == "empty":
            with open(os.path.join(d, "residual_history_main.csv"), "w") as f:
                f.write(cols)
        return d

    def cc_seg(d):
        r = subprocess.run([sys.executable, cc_tool, d, "--segment", "--drop", "1.5"], capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr

    d = seg_run("full")
    rc, out = cc_seg(d)
    chk("--segment: 2 段とも履歴あり -> 従来どおり連結して PASS", rc == 0 and "S1 -> main" in out, "rc=%d" % rc)
    shutil.rmtree(d, ignore_errors=True)
    d = seg_run("missing")
    rc, out = cc_seg(d)
    chk("--segment: 本段の履歴が無い -> 判定不能 (前段だけで合格しない)", rc != 0 and "判定不能" in out and "PASS" not in out,
        "rc=%d %s" % (rc, [ln for ln in out.splitlines() if "判定不能" in ln][:1]))
    shutil.rmtree(d, ignore_errors=True)
    d = seg_run("empty")
    rc, out = cc_seg(d)
    chk("--segment: 本段の履歴が空 -> 判定不能", rc != 0 and "判定不能" in out and "PASS" not in out, "rc=%d" % rc)
    shutil.rmtree(d, ignore_errors=True)

    v, det, _ = cq.classify_series([1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1], 0.4, 0.02, 0.05, 4)
    chk("classify_series(step 全点 1・値一定) -> INDETERMINATE", v == "INDETERMINATE", "verdict=%s" % v)
    v, _, _ = cq.classify_series([0, 1, 2, 5, 4, 6], [1, 1, 1, 1, 1, 1], 0.4, 0.02, 0.05, 4)
    chk("classify_series(step が逆行) -> INDETERMINATE", v == "INDETERMINATE", "verdict=%s" % v)
    v, _, _ = cq.classify_series([0, 1, 2, 3, 4, 5], [1, 1, 1, 1, 1, 1], 0.4, 0.02, 0.05, 4)
    chk("classify_series(step が増える・値一定) -> 従来どおり STEADY", v == "STEADY", "verdict=%s" % v)

    # 上下対称でないメッシュで asym だけを要求 (codex plan 段 judge M3: 評価量ゼロ・severity 0 で STEADY だった)
    d = tempfile.mkdtemp(prefix="qs_gate_")
    rng = np.random.default_rng(3)
    cc_ = np.c_[rng.uniform(0, 1, 40), rng.uniform(0.1, 1, 40), np.zeros(40)]
    with h5py.File(os.path.join(d, "mesh.h5"), "w") as f:
        f.create_dataset("CELLS/centCoords", data=cc_.ravel())
    for st in (100, 200, 300, 400, 500):
        with h5py.File(os.path.join(d, "res_%d.h5" % st), "w") as f:
            for k, val in (("P", 1e5), ("Ux", 10.0), ("Uy", 0.0), ("Uz", 0.0), ("sonic", 340.0)):
                f.create_dataset("VALUE/" + k, data=np.full(40, val))
    r = subprocess.run([sys.executable, cq_tool, d, "--quantity", "asym"], capture_output=True, text=True)
    out = r.stdout + r.stderr
    chk("check_quasisteady --quantity asym (非対称メッシュ・評価量ゼロ) -> INDETERMINATE・非ゼロ終了",
        r.returncode != 0 and "INDETERMINATE" in out and "OVERALL: ALL STEADY" not in out, "rc=%d" % r.returncode)
    r = subprocess.run([sys.executable, cq_tool, d, "--quantity", "machmax"], capture_output=True, text=True)
    chk("check_quasisteady --quantity machmax (一定の場) -> 従来どおり ALL STEADY", r.returncode == 0 and "ALL STEADY" in r.stdout,
        "rc=%d" % r.returncode)
    p = write_csv([[i, 1.0] for i in range(10)], ["step", "C_T"])
    r = subprocess.run([sys.executable, cq_tool, "--series-csv", p, "--series-cols", ","], capture_output=True, text=True)
    chk("check_quasisteady --series-cols が空 -> INDETERMINATE・非ゼロ終了", r.returncode != 0 and "INDETERMINATE" in r.stdout,
        "rc=%d" % r.returncode)
    p2 = write_csv([[1, 1.0]] * 10, ["step", "C_T"])
    r = subprocess.run([sys.executable, cq_tool, "--series-csv", p2, "--series-cols", "C_T"], capture_output=True, text=True)
    chk("check_quasisteady --series-csv (step 全行 1・値一定) -> 非ゼロ終了", r.returncode != 0 and "INDETERMINATE" in r.stdout,
        "rc=%d" % r.returncode)
    os.unlink(p); os.unlink(p2)
    shutil.rmtree(d, ignore_errors=True)


def write_csv(rows, cols):
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)
    return path


# --- check_wall_resolution の合成の run (§5.1 #36q) ---------------------------------------------------------------
WALL_VISC, WALL_RO, WALL_DR = 1.8e-5, 1.2, 1.0e-5      # 粘性 (viscMethod 0)・壁の密度・第一内部節点までの距離 [m]


def wall_tau(y1p):
    """y₁⁺ = dr √(ρ τ) / μ になる接線の壁応力 τ。"""
    return (np.asarray(y1p, float) * WALL_VISC / WALL_DR) ** 2 / WALL_RO


def make_wall_run(d, walls, steps, outlet=True, axisym=True):
    """no-slip 壁を持つ node 方式の小さな run を d に作る (forge なし; check_wall_resolution.py がそのまま読める形)。
    walls: {physID: (名前, kind, r, xs, y1p)} — 壁は半径 r 一定の x 方向の線、各壁節点の内側 WALL_DR に第一内部節点。y1p は壁節点ごとの
    y₁⁺ (接線の壁応力 twall_x をそこから決める)。steps: {名前: [step, ...]} — その壁の出力 `res_<名前>_<physID>_<step>.h5` を書く step
    (無い名前は書かない)。outlet: 出口の出力 `res_outlet_2_<step>.h5` (twall_* なし) を全 step に書く。壁の出力の点の順は逆順
    (順序を仮定しない)。要素は Polyline (CONNE code 2)、軸対称なら面積は 2πr を掛ける。"""
    import h5py
    os.makedirs(d, exist_ok=True)
    coords, wall_ids, planes, bc = [], {}, [], {}
    for pid, (_, _, r, xs, _) in walls.items():
        w0 = len(coords)
        coords += [[x, r, 0.0] for x in xs]
        wall_ids[pid] = list(range(w0, w0 + len(xs)))
        i0 = len(coords)
        coords += [[x, r - WALL_DR, 0.0] for x in xs]
        ip = []
        for i, x in enumerate(xs):
            ds = 0.5 * ((x - xs[i - 1]) if i > 0 else 0.0) + 0.5 * ((xs[i + 1] - x) if i < len(xs) - 1 else 0.0)
            ip.append(len(planes))
            planes.append((w0 + i, None, [0.0, ds, 0.0], ds))          # 壁の境界面 (外向き +y)
        for i in range(len(xs)):
            planes.append((w0 + i, i0 + i, [0.0, -1.0, 0.0], 1.0))      # 壁節点と第一内部節点の間の面
        bc[pid] = ip
    st, sv, sa = [], [], []
    for c0, c1, n, a in planes:
        st += ([1, c0, 1, c0] if c1 is None else [2, c0, c1, 2, c0, c1])
        sv.append(n); sa.append(a)
    coords = np.asarray(coords, float)
    with h5py.File(os.path.join(d, "mesh.h5"), "w") as f:
        f.create_dataset("MESH/COORD", data=coords.ravel())
        f.create_dataset("CELLS/centCoords", data=coords.ravel())
        f.create_dataset("PLANES/STRUCT", data=np.asarray(st, np.int32))
        f.create_dataset("PLANES/surfVect", data=np.asarray(sv, float).ravel())
        f.create_dataset("PLANES/surfArea", data=np.asarray(sa, float))
        for pid in walls:
            f.create_dataset("BCONDS/%d/iCells" % pid, data=np.asarray(wall_ids[pid], np.int32))
            f.create_dataset("BCONDS/%d/iPlanes" % pid, data=np.asarray(bc[pid], np.int32))
    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
        f.write('mesh: {discretization: "node", isAxisymmetric: %d, meshFileName: "mesh.h5"}\n'
                "physProp: {viscMethod: 0, visc: %r}\nturbulence: {model: \"sst\", wallTreatmentSST: 0}\n" % (int(axisym), WALL_VISC))
    with open(os.path.join(d, "bcondConfig.yaml"), "w") as f:
        f.write("inlet:  {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 1.0e5, Tt: 300.0}}\n")
        f.write("outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {Ps: 1.0e3}}\n")
        for pid, (nm, kind, *_r) in walls.items():
            f.write("%s: {physID: %d, kind: %s, outputHDFflg: 1, ints: , floats: {Ts: 300.0}}\n" % (nm, pid, kind))
        f.write("axis:   {physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }\n")
    all_steps = sorted(set(s_ for v in steps.values() for s_ in v))
    for pid, (nm, _, r, xs, y1p) in walls.items():
        n = len(xs)
        perm = list(range(n))[::-1]
        loc = {g: k for k, g in enumerate(perm)}                         # 壁の節点の番号 (x の順) -> 出力の中の番号
        conne = []
        for i in range(n - 1):
            conne += [2, 2, loc[i], loc[i + 1]]
        tau = wall_tau(y1p)
        for st_ in steps.get(nm, []):
            with h5py.File(os.path.join(d, "res_%s_%d_%d.h5" % (nm, pid, st_)), "w") as f:
                f.create_dataset("MESH/COORD", data=coords[[wall_ids[pid][i] for i in perm]].ravel())
                f.create_dataset("MESH/CONNE", data=np.asarray(conne, np.int64))
                vals = {"ro": np.full(n, WALL_RO), "Ts": np.full(n, 300.0), "Ps": np.full(n, 1.0e5),
                        "twall_x": tau, "twall_y": np.zeros(n), "twall_z": np.zeros(n)}
                for k, v in vals.items():
                    f.create_dataset("VALUE/" + k, data=np.asarray(v, float)[perm])
    if outlet:
        for st_ in all_steps:
            with h5py.File(os.path.join(d, "res_outlet_2_%d.h5" % st_), "w") as f:
                f.create_dataset("MESH/COORD", data=np.asarray([[0.4, 0.0, 0.0], [0.4, 0.01, 0.0]]).ravel())
                f.create_dataset("MESH/CONNE", data=np.asarray([2, 2, 0, 1], np.int64))
                for k in ("ro", "Ts", "Ps"):
                    f.create_dataset("VALUE/" + k, data=np.ones(2))
    return d


def wall_cases():
    """2 枚の no-slip 壁 (wall_a: wall_isothermal・r 0.01、wall_b: wall・r 0.05) の入力の組 (§5.1 #36q の事前登録の受入)。
    control: 両壁とも y₁⁺ ≤ 1 (0.5・0.6)。inv_point_fail: wall_a の端の 10 点 (x < 1 mm の細かい要素) だけ y₁⁺ 1.5 →
    点数の割合 10/31 で FAIL・面積の割合 約 0.5 % で PASS。inv_area_fail: wall_b の 1 点 (両隣が 5 cm の要素) だけ y₁⁺ 1.5 →
    点数の割合 1/102 で PASS・面積の割合 約 48 % で FAIL。"""
    xa = np.linspace(0.0, 0.1, 21)
    xb = np.linspace(0.2, 0.3, 11)
    xa_c = np.concatenate([np.linspace(0.0, 0.001, 11)[:-1], np.linspace(0.001, 0.1, 21)])
    xb_c = np.concatenate([np.linspace(0.2, 0.205, 100), [0.255, 0.305]])
    ya_c = np.where(np.arange(len(xa_c)) < 10, 1.5, 0.5)
    yb_c = np.where(np.arange(len(xb_c)) == 100, 1.5, 0.6)
    return {
        "control": {3: ("wall_a", "wall_isothermal", 0.01, xa, np.full(len(xa), 0.5)),
                    5: ("wall_b", "wall", 0.05, xb, np.full(len(xb), 0.6))},
        "inv_point_fail": {3: ("wall_a", "wall_isothermal", 0.01, xa_c, ya_c),
                           5: ("wall_b", "wall", 0.05, xb, np.full(len(xb), 0.6))},
        "inv_area_fail": {3: ("wall_a", "wall_isothermal", 0.01, xa, np.full(len(xa), 0.5)),
                          5: ("wall_b", "wall", 0.05, xb_c, yb_c)},
    }


def run_cwr(d, *args):
    r = subprocess.run([sys.executable, os.path.join(HERE, "check_wall_resolution.py"), d] + [str(a) for a in args],
                       capture_output=True, text=True)
    out = r.stdout + r.stderr
    v = [ln.strip() for ln in out.splitlines() if ln.startswith("VERDICT")]
    return r.returncode, (v[-1] if v else out.strip().splitlines()[-1:]), out


def wall_resolution_cases():
    """check_wall_resolution.py の群と時刻 (§5.1 #36q): 指定の群が全部無ければ・指定の step が無ければ判定不能、面積と点数の重み。"""
    import shutil
    print("=== check_wall_resolution: 群の欠け・step の欠け (§5.1 #36q) ===")
    C = wall_cases()
    S, OLD = 2000, 1000
    both = {"wall_a": [OLD, S], "wall_b": [OLD, S]}
    root = tempfile.mkdtemp(prefix="cwr_gate_")
    try:
        d = make_wall_run(os.path.join(root, "control"), C["control"], both)
        rc, v, _ = run_cwr(d, "--step", S, "--groups", "wall_a,wall_b", "--weight", "area")
        chk("対照 (両壁とも停止の step に出力) -> PASS", rc == 0, "rc=%d %s" % (rc, v))
        rc, v, _ = run_cwr(d, "--groups", "wall_a,wall_b", "--weight", "area")
        chk("--step を省けば従来どおり群ごとの最新で PASS", rc == 0, "rc=%d %s" % (rc, v))
        rc, v, _ = run_cwr(d, "--step", S, "--weight", "area")
        chk("--groups なしは出口の出力 (twall なし) を拾って判定不能 (従来どおり)", rc == 2, "rc=%d %s" % (rc, v))

        d = make_wall_run(os.path.join(root, "b_old_only"), C["control"], {"wall_a": [OLD, S], "wall_b": [OLD]})
        rc, v, out = run_cwr(d, "--step", S, "--groups", "wall_a,wall_b", "--weight", "area")
        chk("片方の壁 (wall_b) の停止の step の出力だけ欠く -> 判定不能 (最新の step に戻らない)",
            rc == 2 and "wall_b" in v, "rc=%d %s" % (rc, v))

        d = make_wall_run(os.path.join(root, "b_absent"), C["control"], {"wall_a": [OLD, S]})
        rc, v, _ = run_cwr(d, "--step", S, "--groups", "wall_a,wall_b", "--weight", "area")
        chk("指定の群 wall_b の出力が 1 つも無い -> 判定不能 (見つかった群だけで判定しない)",
            rc == 2 and "wall_b" in v, "rc=%d %s" % (rc, v))

        d = make_wall_run(os.path.join(root, "old_only"), C["control"], {"wall_a": [OLD], "wall_b": [OLD]})
        rc, v, _ = run_cwr(d, "--step", S, "--groups", "wall_a,wall_b", "--weight", "area")
        chk("停止の step がどの壁にも無い (古い step だけ) -> 判定不能", rc == 2, "rc=%d %s" % (rc, v))

        for name, want_pt, want_ar in (("inv_point_fail", 1, 0), ("inv_area_fail", 0, 1)):
            d = make_wall_run(os.path.join(root, name), C[name], both)
            rp, vp, _ = run_cwr(d, "--step", S, "--groups", "wall_a,wall_b", "--weight", "point")
            ra, va, _ = run_cwr(d, "--step", S, "--groups", "wall_a,wall_b", "--weight", "area")
            chk("%s: 点数の割合と面積の割合で合否が逆 (point rc %d・area rc %d)" % (name, want_pt, want_ar),
                rp == want_pt and ra == want_ar, "point rc=%d %s | area rc=%d %s" % (rp, vp, ra, va))
    finally:
        shutil.rmtree(root, ignore_errors=True)


def main():
    print("=== check_convergence: 不正入力 ===")
    # (1) 残差列が 1 つも無い
    p = write_csv([[i, "outer_end"] for i in range(20)], ["step", "phase"])
    chk("残差列なし -> 合格にしない", *okof(cc.analyze(p, 3.0, 0.2)))
    os.unlink(p)

    # (2) rms_ro だけ (運動量・エネルギーが無い)
    ser = [10.0 ** (-i / 4.0) for i in range(40)]
    p = write_csv([[i, "outer_end", v] for i, v in enumerate(ser)],
                  ["step", "phase", "rms_ro"])
    chk("rms_ro だけ -> 合格にしない (必須保存量の欠損)", *okof(cc.analyze(p, 3.0, 0.2)))
    os.unlink(p)

    # (3) 最終 1 点だけ 0 (drop=inf で合格していた)
    ser = [1.0] * 19 + [0.0]
    cols = ["step", "phase"] + ["rms_" + k for k in ("ro", "roUx", "roUy", "roUz", "roe")]
    p = write_csv([[i, "outer_end"] + [v] * 5 for i, v in enumerate(ser)], cols)
    chk("最終 1 点だけ 0 -> 合格にしない", *okof(cc.analyze(p, 3.0, 0.2)))
    os.unlink(p)

    # (4) 正常入力は従来どおり合格する (回帰)
    ser = [10.0 ** (-i / 4.0) for i in range(40)]
    p = write_csv([[i, "outer_end"] + [v] * 5 for i, v in enumerate(ser)], cols)
    res = cc.analyze(p, 3.0, 0.2)
    chk("正常入力 (5 列 10 桁低下) -> 従来どおり合格",
        res is not None and res[2] is True, "ok=%s" % (None if res is None else res[2]))
    os.unlink(p)

    # (5) 遷移モデルの残差列 (rms_roGamma / rms_roReth): NaN は DIVERGED、遷移が有効な run で列が欠けたら判定不能
    tcols = cols + ["rms_roK", "rms_roOmega", "rms_roGamma", "rms_roReth"]
    rows = [[i, "outer_end"] + [v] * 7 + [v if i < 30 else float("nan"), v] for i, v in enumerate(ser)]
    p = write_csv(rows, tcols)
    res = cc.analyze(p, 3.0, 0.2)
    chk("rms_roGamma に NaN -> DIVERGED", res is not None and res[2] is False and res[3] is True, "ok=%s nan=%s" % (res[2], res[3]))
    os.unlink(p)
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
        f.write('turbulence: {model: "sst", wallTreatmentSST: 0, transition: "lm2009"}\n')
    pth = os.path.join(d, "residual_history.csv")
    with open(pth, "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols); w.writerows([[i, "outer_end"] + [v] * 5 for i, v in enumerate(ser)])
    chk("遷移有効の run で rms_roGamma/rms_roReth が無い -> 合格にしない", *okof(cc.analyze(pth, 3.0, 0.2)))
    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
        f.write('turbulence: {model: "sst", wallTreatmentSST: 0, transition: "none"}\n')
    res = cc.analyze(pth, 3.0, 0.2)
    chk("transition: none なら 5 列で従来どおり合格", res is not None and res[2] is True, "ok=%s" % res[2])
    # --from-floor も同じ必須列検査を通す (codex result M2)。参照 = 遷移なしの 5 列、対象 = 遷移有効だが 2 列欠落
    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
        f.write('turbulence: {model: "sst", wallTreatmentSST: 0, transition: "lm2009"}\n')
    fl = {c: 1e-9 for c in cols[2:]}
    r2 = cc.analyze_from_floor(pth, fl, set(), 1e30, 0.2)
    chk("--from-floor: 遷移有効で 2 列欠落 -> 合格にしない", r2 is not None and r2[2] is False, "ok=%s" % (None if r2 is None else r2[2]))
    d2 = tempfile.mkdtemp()
    with open(os.path.join(d2, "solverConfig.yaml"), "w") as f:
        f.write('turbulence: {model: "sst", wallTreatmentSST: 0}\n')
    refp = os.path.join(d2, "residual_history.csv")
    _sh0 = __import__("shutil"); _sh0.copy(pth, refp)
    tfull = os.path.join(d, "residual_history.csv")
    with open(tfull, "w", newline="") as f:
        w = csv.writer(f); w.writerow(tcols); w.writerows([[i, "outer_end"] + [1e-10] * 9 for i in range(40)])
    fl9 = {c: 1e-9 for c in tcols[2:]}
    r3 = cc.analyze_from_floor(tfull, fl9, set(), 1e30, 0.2, refp)
    chk("--from-floor: 参照 (遷移なし) と対象 (遷移あり) の方程式系が違う -> 合格にしない", r3 is not None and r3[2] is False, "ok=%s" % (None if r3 is None else r3[2]))
    _sh0.rmtree(d2)
    import shutil as _sh; _sh.rmtree(d)
    # 段キー: SST と SST+遷移は別の方程式系 -> 別区間
    import stage_manifest as sm
    ka = sm.stage_key('turbulence: {model: "sst", wallTreatmentSST: 0}\nspace: {convMethod: 1, limiter: 2}\n', "")
    kb = sm.stage_key('turbulence: {model: "sst", wallTreatmentSST: 0, transition: "lm2009"}\nspace: {convMethod: 1, limiter: 2}\n', "")
    chk("SST -> SST+遷移 の段切替は別区間", ka != kb)

    print("=== check_quasisteady: 不正入力 ===")
    v, _, _ = cq.classify([0, 1, 2, 3, 4, 5], [1, 1, 1, 1, 1, float("nan")], 0.4, 0.02, 0.05, 4)
    chk("classify([1,1,1,1,1,NaN]) -> STEADY にしない", v != "STEADY", "verdict=%s" % v)

    v, _, _ = cq.classify_series([float("nan")] * 6, [1, 1, 1, 1, 1, 1], 0.4, 0.02, 0.05, 4)
    chk("classify_series(step 全点 NaN) -> STEADY にしない", v != "STEADY", "verdict=%s" % v)

    v, _, _ = cq.classify_series([0, 1, 2, 3, 4, 5], [1, 1, 1, 1, 1, 1], 0.4, 0.02, 0.05, 4)
    chk("正常系列 -> 従来どおり STEADY", v == "STEADY", "verdict=%s" % v)

    completeness_cases()

    print("=== check_mesh_quality: 成立していないメッシュ ===")
    import shutil
    import subprocess
    import h5py
    tool = os.path.join(HERE, "check_mesh_quality.py")

    def mesh_h5(coord, conne, ncell):
        fd, path = tempfile.mkstemp(suffix=".h5")
        os.close(fd)
        with h5py.File(path, "w") as f:
            f.create_dataset("MESH/COORD", data=np.asarray(coord, float).ravel())
            f.create_dataset("MESH/CONNE", data=np.asarray(conne, np.int64))
            f.create_dataset("VALUE/ro", data=np.ones(ncell))
        return path

    def run_tool(path):
        r = subprocess.run([sys.executable, tool, path, "--mode", "3d"],
                           capture_output=True, text=True)
        return r.returncode, (r.stdout + r.stderr)

    # 体積ゼロの四面体 (4 点が同一平面) — CONNE は [code, nodes...]、tetra の code は 6
    zt = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
    p1 = mesh_h5(zt, [6, 0, 1, 2, 3], 1)
    rc, out = run_tool(p1)
    chk("体積ゼロの四面体 -> 非ゼロ終了", rc != 0, "rc=%d %s" % (rc, out.strip().splitlines()[-1:]))
    os.unlink(p1)

    # NaN 座標を含む
    nt = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, float("nan")]]
    p2 = mesh_h5(nt, [6, 0, 1, 2, 3], 1)
    rc, out = run_tool(p2)
    chk("NaN 座標を含む -> 非ゼロ終了", rc != 0, "rc=%d %s" % (rc, out.strip().splitlines()[-1:]))
    os.unlink(p2)

    # 正常な四面体は通る (回帰)
    gt = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    p3 = mesh_h5(gt, [6, 0, 1, 2, 3], 1)
    rc, out = run_tool(p3)
    chk("正常な四面体 -> 従来どおり合格", rc == 0, "rc=%d" % rc)
    os.unlink(p3)

    # 曲面壁の薄い壁層セル。頂点 Jacobian は 8 つとも正だが、
    # 5 四面体分割の和は面の反り (矢高 0.3 mm) が層厚より大きいので負になり、「向きの不整合」と誤判定していた (2026-09-22)。
    # hexahedron の code は 9。隣に平らなセルを置いて符号を比べさせる
    # case/46 接続模型 (第一層 16 µm) の上フィレット壁層から取った実セル (float32 座標)
    shell = [[0, 0.100000001, 0.100000001],
             [0, 0.0988248587, 0.100000001],
             [0, 0.0988088697, 0.0999836549],
             [0, 0.0999836475, 0.0999836475],
             [0.00430586701, 0.100819342, 0.099532254],
             [0.00430583674, 0.0997356549, 0.0997635275],
             [0.0043058279, 0.0997204334, 0.0997476652],
             [0.0043058577, 0.100803949, 0.0995168611]]
    flat = [[1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [1, 0, 1], [2, 0, 1], [2, 1, 1], [1, 1, 1]]
    p4 = mesh_h5(shell + flat, [9] + list(range(8)) + [9] + list(range(8, 16)), 2)
    rc, out = run_tool(p4 + "" )
    chk("曲面壁の薄い六面体 (頂点 Jacobian は全部正) -> FATAL にしない", "FATAL" not in out, out.strip().splitlines()[-1])
    os.unlink(p4)

    # 局所反転した六面体 (1 頂点を反対側へ押し込む)。体積の和は正のままなので、和だけでは見逃す
    inv = [list(v) for v in flat]; inv[0] = [1.9, 0.9, 0.9]
    p5 = mesh_h5(inv + [[v[0] + 2, v[1], v[2]] for v in flat], [9] + list(range(8)) + [9] + list(range(8, 16)), 2)
    rc, out = run_tool(p5)
    chk("局所反転した六面体 -> 非ゼロ終了", rc != 0 and "局所反転" in out, "rc=%d %s" % (rc, out.strip().splitlines()[-2:]))
    os.unlink(p5)

    # --segment で stage_manifest.json が無い run (継続 run で実際に起きた)。
    # 以前は未定義の `worst` を触って UnboundLocalError で落ち、呼び出し側からは
    # 「判定が出ていない」だけに見えて素通りしていた (2026-09-19)。
    cc_tool = os.path.join(HERE, "check_convergence.py")
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "residual_history.csv"), "w") as f:
        f.write("step,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n")
        for i in range(200):
            v = 10.0 ** (-3 - 3 * i / 199.0)
            f.write("%d,%g,%g,%g,%g,%g\n" % (i, v, v, v, v, v))
    r = subprocess.run([sys.executable, cc_tool, d, "--segment"], capture_output=True, text=True)
    out = r.stdout + r.stderr
    chk("--segment で stage_manifest 無し -> 例外でなく非ゼロ終了",
        r.returncode != 0 and "Traceback" not in out,
        "rc=%d %s" % (r.returncode, "Traceback" if "Traceback" in out else out.strip().splitlines()[-1:]))
    shutil.rmtree(d, ignore_errors=True)

    # 段区間のキーが**層流と SST を区別する**か (2026-09-20 codex result M3)。
    # `turbulenceModel:` だけを見ていたため `turbulence: {model: ...}` 書式で同一キーになり、
    # 層流段と SST 段が 1 区間に連結されていた。
    import stage_manifest as sm
    lam = 'turbulence: {model: "none"}\nspace: {convMethod: 1, limiter: 2}\n'
    sst = 'turbulence: {model: "sst", scalarDiffusion: 1}\nspace: {convMethod: 1, limiter: 2}\n'
    chk("層流段と SST 段が別キーになる", sm.stage_key(lam, "") != sm.stage_key(sst, ""),
        "none=%s sst=%s" % (sm.stage_key(lam, "").get("turbulence.model"),
                            sm.stage_key(sst, "").get("turbulence.model")))
    cm = 'turbulence: {model: "sst"}\nspace: {convMethod: 2, limiter: 2}\n'
    chk("同じ乱流モデルで convMethod だけ違えば別キー",
        sm.stage_key(sst, "") != sm.stage_key(cm, ""), "")
    chk("同一 config は同一キー", sm.stage_key(sst, "x") == sm.stage_key(sst, "x"), "")

    # **書式を変えても区別できるか** (2026-09-21 codex result M4)。正規表現は flow 形式 +
    # 二重引用符しか拾えず、block 形式・単引用符では層流と SST が同一キーになっていた。
    blk_l = 'turbulence:\n  model: "none"\nspace:\n  convMethod: 1\n'
    blk_s = 'turbulence:\n  model: "sst"\nspace:\n  convMethod: 1\n'
    chk("block 形式でも層流と SST が別キー", sm.stage_key(blk_l, "") != sm.stage_key(blk_s, ""),
        "none=%s sst=%s" % (sm.stage_key(blk_l, "").get("turbulence.model"),
                            sm.stage_key(blk_s, "").get("turbulence.model")))
    sq_l = "turbulence: {model: 'none'}\n"
    sq_s = "turbulence: {model: 'sst'}\n"
    chk("単引用符でも層流と SST が別キー", sm.stage_key(sq_l, "") != sm.stage_key(sq_s, ""),
        "none=%s sst=%s" % (sm.stage_key(sq_l, "").get("turbulence.model"),
                            sm.stage_key(sq_s, "").get("turbulence.model")))

    # メッシュの来歴 (mesh_provenance.py check、plan tooling-generator-defaults-and-launch-gate §5.1 #7)。
    # 記録が無い・壊れている・別メッシュのものを「照合できた」にしない (詳しい負例は tests/unit/test_mesh_provenance.py)。
    print("=== mesh_provenance: 記録の無い・壊れた・別メッシュの記録 ===")
    import mesh_provenance as mpv
    mp_tool = os.path.join(HERE, "mesh_provenance.py")

    def mp_run(d, shift=0.0):
        with h5py.File(os.path.join(d, "mesh.h5"), "w") as f:
            f.create_dataset("MESH/COORD", data=np.array([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0], float) + shift)
            f.create_dataset("MESH/CONNE", data=np.array([5, 0, 1, 2, 3], np.int64))
            f.create_dataset("VALUE/ro", data=np.ones(1, np.float32))
            f.create_dataset("VALUE/wall_dist", data=np.ones(1))
            f.create_dataset("BCONDS/1/iCells", data=np.array([0], np.int64))
            f["BCONDS/1"].attrs["bcondKind"] = "wall"
        cfg = 'mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5", discretization: "cell"}\n'
        bc = "wall: {physID: 1, kind: wall, outputHDFflg: 0, ints: , floats: }\n"
        for nm, txt in (("solverConfig.yaml", cfg), ("bcondConfig.yaml", bc)):
            with open(os.path.join(d, nm), "w") as fh:
                fh.write(txt)
        return cfg, bc

    def mp_check(d):
        r = subprocess.run([sys.executable, mp_tool, "check", d], capture_output=True, text=True)
        return r.returncode, (r.stdout.strip().splitlines() or [""])[-1]

    da, db = tempfile.mkdtemp(), tempfile.mkdtemp()
    cfg, bc = mp_run(da)
    rc, last = mp_check(da)
    chk("記録の無いメッシュ -> 非ゼロ終了", rc != 0 and last == "VERDICT: FAIL", "rc=%d %s" % (rc, last))
    mpv.write_record(os.path.join(da, "mesh.h5"), mpv.build_record(os.path.join(da, "mesh.h5"), cfg, bc))
    rc, last = mp_check(da)
    chk("正常な記録 -> 従来どおり合格", rc == 0 and last == "VERDICT: PASS", "rc=%d %s" % (rc, last))
    mp_run(db, shift=0.5)
    shutil.copy(os.path.join(da, "mesh.h5.provenance.json"), os.path.join(db, "mesh.h5.provenance.json"))
    rc, last = mp_check(db)
    chk("別メッシュの記録 -> 非ゼロ終了", rc != 0 and last == "VERDICT: FAIL", "rc=%d %s" % (rc, last))
    with open(os.path.join(da, "mesh.h5.provenance.json"), "w") as fh:
        fh.write("{")
    rc, last = mp_check(da)
    chk("壊れた記録 -> 非ゼロ終了", rc != 0 and last == "VERDICT: FAIL", "rc=%d %s" % (rc, last))
    shutil.rmtree(da, ignore_errors=True)
    shutil.rmtree(db, ignore_errors=True)

    # 同一メッシュの restart (restart_field.py、plan tooling-generator-defaults-and-launch-gate §5.1 #20)。
    # 2026-10-15 まで、座標を float32 に丸めて比べ、接続を見ず、座標が無くても続行していた (codex diagnose の再現)。
    # 別の格子を「同一メッシュ」にしない (詳しい負例は tests/unit/test_restart_field_mesh_identity.py)。
    print("=== restart_field: 同一メッシュでない入力 ===")
    rf_tool = os.path.join(HERE, "restart_field.py")
    rroot = tempfile.mkdtemp()
    xs = 0.1234567891 + np.linspace(0.0, 1.3, 5)
    base = np.c_[np.repeat(xs, 3), np.tile(np.linspace(0.0, 0.7, 3), 5), np.zeros(15)]
    quads = [[i * 3 + j, (i + 1) * 3 + j, (i + 1) * 3 + j + 1, i * 3 + j + 1] for i in range(4) for j in range(2)]

    def rf_run(name, coord=base, conne=quads, drop=()):
        d = os.path.join(rroot, name)
        os.makedirs(d)
        with open(os.path.join(d, "solverConfig.yaml"), "w") as fh:
            fh.write('mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5", discretization: "node"}\n'
                     'physProp: {thermalMethod: 0}\n')
        with h5py.File(os.path.join(d, "mesh.h5"), "w") as f:
            f.create_dataset("MESH/COORD", data=np.asarray(coord, np.float64).ravel())
            f.create_dataset("MESH/CONNE", data=np.array([x for n in range(15) for x in (1, n)], np.int32))
            f.create_dataset("VIZMESH/CONNE", data=np.array([x for q in conne for x in [5] + q], np.int32))
            f.create_dataset("BCONDS/1/iCells", data=np.array([0, 3, 6, 9, 12], np.int32))
            f["BCONDS/1"].attrs["bcondKind"] = "wall"
            for k in ("ro", "roUx", "roUy", "roUz", "roe"):
                f.create_dataset("VALUE/" + k, data=np.ones(15, np.float32))
            f.create_dataset("VALUE/wall_dist", data=np.asarray(coord)[:, 1])
            for k in drop:
                del f[k]
        return d

    def rf_res(d):
        """forge の出力と同じ置き方の res (float32 の座標・primal の接続・/VALUE だけ)。"""
        with h5py.File(os.path.join(d, "mesh.h5"), "r") as m, h5py.File(os.path.join(d, "res_10.h5"), "w") as f:
            f.create_dataset("MESH/COORD", data=np.asarray(m["MESH/COORD"]).astype(np.float32))
            f.create_dataset("MESH/CONNE", data=np.asarray(m["VIZMESH/CONNE"]))
            for k in ("ro", "roUx", "roUy", "roUz", "roe"):
                f.create_dataset("VALUE/" + k, data=np.full(15, 2.0, np.float32))
        return os.path.join(d, "res_10.h5")

    def rf_call(src, dst):
        env = dict(os.environ, FORGE_BIN="/nonexistent/forge")
        env.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
        r = subprocess.run([sys.executable, rf_tool, src, dst], capture_output=True, text=True, env=env)
        return r.returncode, (r.stdout + r.stderr).strip().splitlines()

    ra = rf_run("a")
    res_a = rf_res(ra)
    pert = base * (1.0 + 1e-8)
    pert = np.where(pert.astype(np.float32) == base.astype(np.float32), pert, base)   # float32 では元と同じ
    q2 = [list(q) for q in quads]
    q2[0] = q2[0][1:] + q2[0][:1]
    for name, d in (("座標を ~1e-8 変えた別の入力 (float32 では同じ; codex の再現)", rf_run("pert", coord=pert)),
                    ("接続だけ違う入力", rf_run("conne", conne=q2)),
                    ("宛先の格子に /MESH/CONNE が無い", rf_run("noconne", drop=("MESH/CONNE",)))):
        rc, lines = rf_call(res_a, os.path.join(d, "mesh.h5"))
        chk("restart_field: %s -> 非ゼロ終了" % name, rc != 0 and any("REFUSED" in ln for ln in lines),
            "rc=%d %s" % (rc, (lines or [""])[0][:90]))
    lone = os.path.join(rroot, "lone")
    os.makedirs(lone)
    shutil.copy(res_a, lone)
    rc, lines = rf_call(os.path.join(lone, "res_10.h5"), os.path.join(rf_run("copy_a"), "mesh.h5"))
    chk("restart_field: 格子を決められない res (隣に config なし) -> 非ゼロ終了", rc != 0 and any("REFUSED" in ln for ln in lines),
        "rc=%d %s" % (rc, (lines or [""])[0][:90]))
    rc, lines = rf_call(res_a, os.path.join(rroot, "copy_a", "mesh.h5"))
    chk("restart_field: 同じ格子の複製 -> 従来どおり合格", rc == 0 and lines[-1].startswith("VERDICT: OK"),
        "rc=%d %s" % (rc, lines[-1][:60] if lines else ""))
    shutil.rmtree(rroot, ignore_errors=True)

    # 投入の関門 (launch_gate.py、plan tooling-generator-defaults-and-launch-gate §5.1 #10)。2026-10-15 まで run_case.sh は
    # 設定・メッシュを照合せずに forge を起動していた (check_solver_config の未知のキーは WARN・終了 0、品質の検査は呼ばれない)。
    # 検査できない入力・メッシュの FAIL・禁止キーを「起動してよい」にしない (詳しい負例は tests/unit/test_launch_gate.py)。
    print("=== launch_gate: 投入の関門 ===")
    lg_tool = os.path.join(HERE, "launch_gate.py")
    lroot = tempfile.mkdtemp()

    def lg_run(name, cfg_extra="", first=None, bcond=True):
        d = os.path.join(lroot, name)
        os.makedirs(d)
        ni, nj = 11, 4
        ys = np.linspace(0.0, 0.2, nj) if first is None else np.r_[np.linspace(0.0, 0.2 - first, nj - 1), 0.2]
        X, Y = np.meshgrid(np.linspace(0.0, 1.0, ni), ys, indexing="ij")
        nid = lambda i, j: i * nj + j                                                   # noqa: E731
        quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
        wall = [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)]
        with h5py.File(os.path.join(d, "mesh.h5"), "w") as f:
            f.create_dataset("MESH/COORD", data=np.c_[X.ravel(), Y.ravel(), np.zeros(ni * nj)].ravel())
            f.create_dataset("MESH/CONNE", data=np.array([x for n in range(ni * nj) for x in (1, n)], np.int64))
            f.create_dataset("VIZMESH/CONNE", data=np.array([x for q in quads for x in [5] + q], np.int64))
            f["VIZMESH"].attrs["nVizCells"] = len(quads)
            f.create_dataset("VALUE/wall_dist", data=0.2 - Y.ravel())
            f.create_dataset("BCONDS/3/vizBfaceNodes", data=np.array(wall, np.int64).ravel())
            f.create_dataset("BCONDS/3/vizBfaceSizes", data=np.full(len(wall), 2, np.int64))
            f["BCONDS/3"].attrs["bcondKind"] = "wall"
        with open(os.path.join(d, "solverConfig.yaml"), "w") as fh:
            fh.write('mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5", discretization: "node"}\n' + cfg_extra)
        if bcond:
            with open(os.path.join(d, "bcondConfig.yaml"), "w") as fh:
                fh.write("wall: {physID: 3, kind: wall, outputHDFflg: 0, ints: , floats: }\n")
        return d

    def lg_call(d, env_extra=None):
        env = dict(os.environ)
        env.pop("FORGE_LAUNCH_GATE", None)
        env.update(env_extra or {})
        r = subprocess.run([sys.executable, lg_tool, d, "--dry-run", "--forge-bin", os.path.join(lroot, "no_forge")],
                           capture_output=True, text=True, env=env)
        last = [ln for ln in r.stdout.splitlines() if ln.startswith("VERDICT")]
        return r.returncode, (last or [r.stdout.strip()[-120:]])[-1]

    rc, last = lg_call(lg_run("ok"))
    chk("launch_gate: 宣言なし・記録なし・品質 PASS (移行期) -> 起動してよい (WARN)", rc == 0 and last.startswith("VERDICT: WARN"),
        "rc=%d %s" % (rc, last[:60]))
    for name, d, env in (("宣言なし + 強制の未知のキー (check_solver_config では WARN・終了 0)", lg_run("unknown", "space: {notAKey: 1}\n"),
                          {"FORGE_LAUNCH_GATE": "enforce"}),
                         ("宣言なしの起動時に拒否されるキー (移行期でも)", lg_run("rejected", "space: {limiterMatchRecon: 1}\n"), None),
                         ("bcondConfig.yaml が無い (--bcond を渡せない)", lg_run("nobc", bcond=False), None),
                         ("AR 6000 の壁の層 (緩和の上限 5000 を超える)", lg_run("ar6000", first=0.1 / 6000.0), None),
                         ("FORGE_LAUNCH_GATE の綴り違い", lg_run("typo"), {"FORGE_LAUNCH_GATE": "wran"}),
                         ("宣言なし + FORGE_LAUNCH_GATE=enforce", lg_run("enforce"), {"FORGE_LAUNCH_GATE": "enforce"})):
        rc, last = lg_call(d, env)
        chk("launch_gate: %s -> 非ゼロ終了" % name, rc != 0 and last.startswith("VERDICT: FAIL"), "rc=%d %s" % (rc, last[:60]))
    d = lg_run("bnd_warn")
    with open(os.path.join(d, "solverConfig.yaml"), "w") as fh:
        fh.write('mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5", discretization: "node", bndFirstOrder: 1}\n')
    rc, last = lg_call(d, {"FORGE_LAUNCH_GATE": "warn"})
    chk("launch_gate: bndFirstOrder + FORGE_LAUNCH_GATE=warn (禁止は落とさない) -> 非ゼロ終了",
        rc != 0 and last.startswith("VERDICT: FAIL"), "rc=%d %s" % (rc, last[:60]))
    d = lg_run("broken_decl")
    with open(os.path.join(d, "RECIPE.json"), "w") as fh:
        fh.write("{")
    rc, last = lg_call(d)
    chk("launch_gate: 読めない宣言 (RECIPE.json) -> 非ゼロ終了", rc != 0 and last.startswith("VERDICT: FAIL"), "rc=%d %s" % (rc, last[:60]))
    shutil.rmtree(lroot, ignore_errors=True)

    wall_resolution_cases()

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
