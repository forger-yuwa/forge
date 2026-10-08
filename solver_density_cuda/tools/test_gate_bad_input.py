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

    floor_events_bad_input()

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def write_floor_csv(d, N, path="steady_implicit", events=None, drop=(), end=True, audit="done", sessions=None,
                    overflow_q=None, uneval_q=None, mismatch_q=None, truncate_last=False):
    """check_floor_events の合成記録。events = {q_index: (nT, nRho, nP)}、drop = 欠けさせる q_index、
    sessions = 先に置く (完走した) session の step 数のリスト (restart / 段階起動の前段)。"""
    import check_floor_events as fe
    cols = fe.COLUMNS
    lines = [",".join(cols)]

    def row(kind, sess, step, inner, q, ev=(0, 0, 0), ovf=0, unev=0, mis=0, note=""):
        r = {c: "0" for c in cols}
        r.update(kind=kind, session=str(sess), step=str(step), inner=str(inner), q_index=str(q), path=path,
                 n_real="100", n_ghost="10", ids_T="", ids_Rho="", ids_P="", ids_near="", ids_truncated="",
                 overflow=str(ovf), note=note)
        r["nT_real"], r["nRho_real"], r["nP_real"] = (str(v) for v in ev)
        r["nT_uneval_real"], r["nT_mismatch_real"] = str(unev), str(mis)
        if ev[0]:
            r["ids_T"] = "7"
        return ",".join(r[c] for c in cols)

    allsess = list(sessions or []) + [N]
    for si, n in enumerate(allsess, start=1):
        last = (si == len(allsess))
        lines.append(row("session_begin", si, 0, -1, -1, note="format=1;path=%s;prior_sessions=%d" % (path, si - 1)))
        lines.append(row("init", si, 0, 0, -1))
        for k in range(1, n + 1):
            q = k - 1
            if last and q in drop:
                continue
            lines.append(row("eos", si, k, 0, q, (events or {}).get(q, (0, 0, 0)) if last else (0, 0, 0),
                             ovf=int(last and q == overflow_q), unev=int(last and q == uneval_q), mis=int(last and q == mismatch_q)))
        if last and not end:
            break
        if audit == "done" and not (last and n in drop):
            lines.append(row("audit", si, n, 0, n, (events or {}).get(n, (0, 0, 0)) if last else (0, 0, 0)))
        lines.append(row("session_end", si, n, -1, -1, note="audit=%s" % audit))
    if truncate_last:
        lines[-1] = lines[-1][:len(lines[-1]) // 2]
    with open(os.path.join(d, "floor_events.csv"), "w") as f:
        f.write("\n".join(lines) + ("" if truncate_last else "\n"))


def floor_events_bad_input():
    """毎更新の床事象の判定 (check_floor_events.judge) が **記録の欠け・区間の欠け・未検証の経路を合格にしない** こと
    (plans/active/tooling-sern-te-wake-grid.md §4「床の判定」、§5.1 #2)。記録が無いことを 0 件として扱わない。"""
    import shutil
    import check_floor_events as fe
    print("=== check_floor_events: 不正入力 ===")
    N = 40
    cases = [
        ("全更新 0 件 → PASS (対照)", {}, "PASS"),
        ("記録ファイルが無い → 判定不能", {"_nofile": True}, "INDETERMINATE"),
        ("途中で 1 件 (q 25)・最終場 (監査) は正常 → FAIL", {"events": {25: (1, 0, 0)}}, "FAIL"),
        ("最後の更新 (監査行) だけ 1 件 → FAIL", {"events": {N: (0, 0, 1)}}, "FAIL"),
        ("区間の step が 1 つ欠けた → 判定不能 (欠落 ≠ 0 件)", {"drop": (30,)}, "INDETERMINATE"),
        ("監査行が欠けた → 判定不能", {"drop": (N,)}, "INDETERMINATE"),
        ("区間の入口 (q = a) が欠けた → 判定不能", {"drop": (20,)}, "INDETERMINATE"),
        ("session_end が無い (異常終了) → 判定不能", {"end": False}, "INDETERMINATE"),
        ("最後の行が書きかけ → 判定不能", {"truncate_last": True}, "INDETERMINATE"),
        ("監査が未対応 (audit=unsupported) → 判定不能", {"audit": "unsupported:condEquilibrium2"}, "INDETERMINATE"),
        ("未検証の経路 (dual_time) → 判定不能", {"path": "dual_time"}, "INDETERMINATE"),
        ("件数の整合が取れない (overflow) → 判定不能", {"overflow_q": 33}, "INDETERMINATE"),
        ("温度の述語を評価していない節点 (二相) → 判定不能", {"uneval_q": 33}, "INDETERMINATE"),
        ("述語と EOS の最終温度が食い違う → 判定不能", {"mismatch_q": 33}, "INDETERMINATE"),
        ("区間の外 (q 5) の事象は合否に使わない → PASS", {"events": {5: (3, 0, 0)}}, "PASS"),
        ("入口 Q_a (q 20) の事象は別枠 → PASS", {"events": {20: (1, 0, 0)}}, "PASS"),
    ]
    for name, kw, want in cases:
        d = tempfile.mkdtemp()
        kw = dict(kw)
        if not kw.pop("_nofile", False):
            write_floor_csv(d, N, **kw)
        r = fe.judge(d, tail=0.5)
        chk(name, r["verdict"] == want and (r["ok"] == (want == "PASS")), "%s %s" % (r["verdict"], (r["reasons"] or [""])[0][:60]))
        shutil.rmtree(d, ignore_errors=True)
    # restart: 前の session (完走 40 step) の後に restart した session が 15 step しかない。判定区間 20 step は
    # restart 後の session の始まりより前に及ぶ → 前の session の記録で埋めずに判定不能
    d = tempfile.mkdtemp()
    write_floor_csv(d, 15, sessions=[40])
    r = fe.judge(d, window_steps=20)
    chk("restart 後の区間の不足 (前の session で埋めない) → 判定不能", r["verdict"] == "INDETERMINATE" and not r["ok"], (r["reasons"] or [""])[0][:60])
    r = fe.judge(d, window_steps=10)
    chk("restart 後の session の中に区間が収まれば判定する → PASS", r["verdict"] == "PASS", r["verdict"])
    shutil.rmtree(d, ignore_errors=True)
    # 記録と res_*.h5 の最終 step が違う (古い記録の取り違え)
    d = tempfile.mkdtemp()
    write_floor_csv(d, N)
    open(os.path.join(d, "res_%d.h5" % (N + 10)), "w").close()
    r = fe.judge(d, tail=0.5)
    chk("記録の最終 step と res_*.h5 が違う → 判定不能", r["verdict"] == "INDETERMINATE", (r["reasons"] or [""])[0][:60])
    shutil.rmtree(d, ignore_errors=True)
    # 見出し違い
    d = tempfile.mkdtemp()
    write_floor_csv(d, N)
    txt = open(os.path.join(d, "floor_events.csv")).read().replace("kind,session", "kind,sess", 1)
    open(os.path.join(d, "floor_events.csv"), "w").write(txt)
    chk("見出しが違う → 判定不能", fe.judge(d)["verdict"] == "INDETERMINATE", "")
    shutil.rmtree(d, ignore_errors=True)
    # CLI は PASS 以外で非ゼロ終了
    import subprocess
    d = tempfile.mkdtemp()
    write_floor_csv(d, N, drop=(30,))
    rr = subprocess.run([sys.executable, os.path.join(HERE, "check_floor_events.py"), d], capture_output=True, text=True)
    chk("CLI: 判定不能は非ゼロ終了で VERDICT: INDETERMINATE", rr.returncode != 0 and "VERDICT: INDETERMINATE" in rr.stdout, "rc=%d" % rr.returncode)
    shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
