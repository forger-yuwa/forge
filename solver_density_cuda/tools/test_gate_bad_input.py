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

    wall_resolution_cases()

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
