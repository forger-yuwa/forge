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

    print("=== check_quasisteady: 不正入力 ===")
    v, _, _ = cq.classify([0, 1, 2, 3, 4, 5], [1, 1, 1, 1, 1, float("nan")], 0.4, 0.02, 0.05, 4)
    chk("classify([1,1,1,1,1,NaN]) -> STEADY にしない", v != "STEADY", "verdict=%s" % v)

    v, _, _ = cq.classify_series([float("nan")] * 6, [1, 1, 1, 1, 1, 1], 0.4, 0.02, 0.05, 4)
    chk("classify_series(step 全点 NaN) -> STEADY にしない", v != "STEADY", "verdict=%s" % v)

    v, _, _ = cq.classify_series([0, 1, 2, 3, 4, 5], [1, 1, 1, 1, 1, 1], 0.4, 0.02, 0.05, 4)
    chk("正常系列 -> 従来どおり STEADY", v == "STEADY", "verdict=%s" % v)

    print("=== check_mesh_quality: 成立していないメッシュ ===")
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

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
