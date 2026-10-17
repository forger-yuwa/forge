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
| `check_wall_resolution` (2026-10-17 追加) | `--step` の壁ダンプが無い壁がある | 最新の出力に代えて `PASS`, exit 0 |
| 〃 | `--groups` に存在しない壁の名前 | 飛ばして `PASS`, exit 0 |
| 〃 | メッシュの BCONDS に無い physID の壁ダンプ | Traceback, exit 1 (呼び出し側の FAIL) |
| 〃 (2026-10-17 2 回目) | `--groups` 無しで bcondConfig の壁のダンプが一度も無い | 存在するダンプだけで `PASS`, exit 0 |
| 〃 | `--groups` 無しで bcondConfig を読めない | 警告してダンプだけで `PASS`, exit 0 |
"""
import csv
import os
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

    print("=== check_wall_resolution: 指定の step・壁の出力がそろわない ===")
    # 2026-10-17 codex result M2 (plan tooling-nozzle-card-local-mesh §5.1 #7b): `--step` の出力が無い壁を最新の出力に代え、
    # `--groups` の存在しない名前を飛ばして PASS・終了コード 0 にしていた (case/71 run_0014 で `--step 46170 --groups wall,missing_wall`
    # が step 46171 の wall だけで PASS)。メッシュに無い physID のダンプは KeyError の Traceback で終了コード 1 (= 呼び出し側の FAIL) だった。
    cwr_tool = os.path.join(HERE, "check_wall_resolution.py")
    d = tempfile.mkdtemp()
    # 壁 2 面 (wall_a = physID 1 の DOF 0、wall_b = physID 2 の DOF 2)、それぞれ法線方向 1e-5 m に第一内部点 (y1+ ≈ 0.54 < 1)
    coords = np.array([[0, 0, 0], [0, 0, 1e-5], [1, 0, 0], [1, 0, 1e-5]], float)
    with h5py.File(os.path.join(d, "m.h5"), "w") as f:
        f.create_dataset("MESH/COORD", data=coords.ravel())
        f.create_dataset("CELLS/centCoords", data=coords.ravel())
        f.create_dataset("PLANES/STRUCT", data=np.array([3, 0, 0, 0, 2, 0, 1, 3, 0, 0, 0, 2, 2, 3], np.int32))
        f.create_dataset("PLANES/surfVect", data=np.array([0, 0, -1.0, 0, 0, -1.0]))
        f.create_dataset("PLANES/surfArea", data=np.ones(2))
        for pid, dof, pl in ((1, 0, 0), (2, 2, 1)):
            f.create_dataset("BCONDS/%d/iCells" % pid, data=np.array([dof], np.int32))
            f.create_dataset("BCONDS/%d/iPlanes" % pid, data=np.array([pl], np.int32))
    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
        f.write("mesh: {meshFileName: m.h5}\nphysProp: {viscMethod: 1}\n")

    def dump(name, pid, step, xyz):
        with h5py.File(os.path.join(d, "res_%s_%d_%d.h5" % (name, pid, step)), "w") as f:
            f.create_dataset("MESH/COORD", data=np.asarray(xyz, float).ravel())
            for k, v in (("twall_x", 1.0), ("twall_y", 0.0), ("twall_z", 0.0), ("ro", 1.0), ("Ts", 300.0)):
                f.create_dataset("VALUE/" + k, data=np.full(1, v))

    def cwr(*args):
        r = subprocess.run([sys.executable, cwr_tool, d] + list(args), capture_output=True, text=True)
        out = r.stdout + r.stderr
        last = [ln for ln in out.splitlines() if ln.startswith("VERDICT")] or out.strip().splitlines()[-1:]
        return r.returncode, out, last[-1] if last else ""

    for s in (100, 200):
        dump("wall_a", 1, s, coords[[0]])
    dump("wall_b", 2, 100, coords[[2]])          # wall_b は step 100 だけ
    rc, out, v = cwr("--step", "200")
    chk("--step 200 で wall_b の 200 が無い -> 判定不能 (100 に代えない)", rc == 2 and "INDETERMINATE" in v and "wall_b" in v, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "150", "--groups", "wall_a")
    chk("--step 150 (どの壁にも無い step) -> 判定不能 (最新 200 に代えない)", rc == 2 and "INDETERMINATE" in v, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "200", "--groups", "wall_a,missing_wall")
    chk("--groups に存在しない壁の名前 -> 判定不能", rc == 2 and "INDETERMINATE" in v and "missing_wall" in v, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "200", "--groups", "wall_a")
    chk("--step 200 --groups wall_a (そろっている) -> 従来どおり PASS・評価した壁と step を出す",
        rc == 0 and "VERDICT: PASS" in v and "EVALUATED: wall_a/1@200" in out, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "100")
    chk("--step 100 (全部の壁にある) -> PASS、両方の壁を step 100 で評価",
        rc == 0 and "EVALUATED: wall_a/1@100, wall_b/2@100" in out, "rc=%d %s" % (rc, v))
    rc, out, v = cwr()
    chk("--step 無し -> 従来どおり壁ごとの最新 (wall_a 200・wall_b 100) で PASS",
        rc == 0 and "EVALUATED: wall_a/1@200, wall_b/2@100" in out, "rc=%d %s" % (rc, v))
    dump("wall_c", 7, 200, coords[[0]])          # physID 7 はメッシュの BCONDS に無い
    rc, out, v = cwr("--groups", "wall_c")
    chk("メッシュに無い physID の壁ダンプ -> 判定不能・終了コード 2 (Traceback にしない)",
        rc == 2 and "INDETERMINATE" in v and "Traceback" not in out, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "100")
    chk("bcondConfig.yaml が無い run (--groups 無し) -> 今どおりダンプから数え、そのことを出力に書く",
        rc == 2 and "bcondConfig.yaml が無い -> ある壁ダンプから数える" in out, "rc=%d" % rc)   # wall_c (physID 7) で判定不能のまま

    print("=== check_wall_resolution: --groups 無しで bcondConfig の壁のダンプがそろわない ===")
    # 2026-10-17 codex result 2 回目 M2 (plan tooling-nozzle-card-local-mesh §5.1 #7b の続き): --groups を省くと存在するダンプから
    # 対象を数えたので、bcondConfig の no-slip 壁でダンプが一度も無い壁を見落として PASS・終了コード 0 だった (受け渡しの関門
    # deltastar_loop.wall_resolution_verdict はこの経路)。いまは bcondConfig の kind wall* を全部必須にする
    d0 = d                                       # 同じメッシュ・設定で別の run (ダンプは作り直す)
    d = tempfile.mkdtemp()
    for f_ in ("m.h5", "solverConfig.yaml"):
        shutil.copy2(os.path.join(d0, f_), d)
    shutil.rmtree(d0, ignore_errors=True)
    for s_ in (100, 200):
        dump("wall_a", 1, s_, coords[[0]])
    dump("wall_b", 2, 100, coords[[2]])          # wall_b は step 100 だけ
    bc = ("wall_a: {physID: 1, kind: wall_isothermal, outputHDFflg: 1}\n"
          "wall_b: {physID: 2, kind: wall_adiabatic, outputHDFflg: 1}\n"
          "outlet: {physID: 9, kind: outlet_statPress, outputHDFflg: 1}\n")
    with open(os.path.join(d, "bcondConfig.yaml"), "w") as f:
        f.write(bc)
    rc, out, v = cwr("--step", "200")
    chk("bcond あり・--groups 無し: 存在する壁 wall_b の指定 step 200 が欠ける -> 判定不能 (100 に代えない)",
        rc == 2 and "INDETERMINATE" in v and "wall_b" in v and "EVALUATED" not in out, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "100")
    chk("bcond あり・--groups 無し: 全部の壁に step 100 -> PASS、対象は bcond の壁全部と出力に書く",
        rc == 0 and "EVALUATED: wall_a/1@100, wall_b/2@100" in out and "kind wall*) 全部 = wall_a/1, wall_b/2" in out, "rc=%d %s" % (rc, v))
    rc, out, v = cwr()
    chk("bcond あり・--groups 無し・--step 無し: 壁ごとの最新 (wall_a 200・wall_b 100) で PASS (今どおり)",
        rc == 0 and "EVALUATED: wall_a/1@200, wall_b/2@100" in out, "rc=%d %s" % (rc, v))
    with open(os.path.join(d, "bcondConfig.yaml"), "w") as f:
        f.write(bc + "wall_c: {physID: 3, kind: wall_isothermal, outputHDFflg: 1}\n")   # ダンプが一度も無い壁
    rc, out, v = cwr()
    chk("bcond の壁 wall_c のダンプが一度も無い・--step 無し -> 判定不能・終了コード 2 (修正前は wall_a・wall_b だけで PASS・0)",
        rc == 2 and "INDETERMINATE" in v and "wall_c" in v and "一度も無い" in v and "EVALUATED" not in out, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "100")
    chk("bcond の壁 wall_c のダンプが一度も無い・--step 100 -> 判定不能・終了コード 2 (修正前は PASS・0)",
        rc == 2 and "INDETERMINATE" in v and "wall_c" in v, "rc=%d %s" % (rc, v))
    rc, out, v = cwr("--step", "100", "--groups", "wall_a,wall_b")
    chk("--groups を明示すれば bcond の他の壁は要らない (今どおり PASS; 呼び出し側が壁を選ぶ)",
        rc == 0 and "EVALUATED: wall_a/1@100, wall_b/2@100" in out, "rc=%d %s" % (rc, v))
    with open(os.path.join(d, "bcondConfig.yaml"), "w") as f:
        f.write("wall_a: {physID: [1, 2\n")                 # 読めない YAML
    rc, out, v = cwr("--step", "100")
    chk("bcondConfig.yaml を読めない・--groups 無し -> 必須の壁を作れず判定不能 (前は警告してダンプだけで PASS)",
        rc == 2 and "INDETERMINATE" in v and "必須の壁の集合を作れない" in v, "rc=%d %s" % (rc, v))
    shutil.rmtree(d, ignore_errors=True)

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
