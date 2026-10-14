#!/usr/bin/env python3
"""`check_case_gates.py` の判定を**本番関数を直接呼んで**検査する (2026-09-21)。

codex result レビュー (2026-09-20 M2/M4, 2026-09-21 M1/M2/M8) が再現した誤合格の回帰試験。
判定ロジックを触ったら必ず通すこと:

    python3 case/49.plate_annular_cavity_m5/tools/test_case_gates.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_case_gates as g          # noqa: E402
import check_cavity_steady as cs      # noqa: E402

ok = True


def chk(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print("  [%s] %-56s %s" % ("OK " if cond else "FAIL", name, detail))


def main():
    # --- M2: 収支ゲートは**全流束**の残差を見る ---
    fld = {"budget_residual_alldepth": 0.0, "budget_residual_all_alldepth": 0.20}
    chk("対流のみ 0 / 全項 20 % は全項のキーを選ぶ",
        g.budget_key(fld) == "budget_residual_all_alldepth", g.budget_key(fld))
    chk("全項キーが無い JSON は旧版として弾ける (キーが一致しない)",
        g.budget_key({"budget_residual_alldepth": 0.0}) != "budget_residual_all_alldepth")
    chk("収支キーが 1 つも無ければ None", g.budget_key({}) is None)

    # --- M1: 継続 run の「判定不能」は合格にしない ---
    chk("判定不能の出力は不合格",
        g.residual_inconclusive("必須の保存量残差列が無い: rms_roUx  <-- 判定不能", 1))
    chk("継続 run 正常の NOT CONVERGED は不合格にしない",
        not g.residual_inconclusive("=== NOT CONVERGED (stalled/plateau) ===", 1))
    chk("異常終了 (rc>=2) は不合格", g.residual_inconclusive("", 2))

    # --- M8: 準定常の絶対許容が plan §6.4 (2 K / 0.5 mm) と一致する ---
    for k in ("dT_mouth", "dT_mid", "dT_floor"):
        unit, floor_, rel_ = cs.ABS_TOL[k]
        chk("%s の許容は 2 K 固定 (相対緩和なし)" % k,
            unit == "K" and floor_ == 2.0 and rel_ == 0.0, "floor=%g rel=%g" % (floor_, rel_))
    unit, floor_, rel_ = cs.ABS_TOL["zpen_25"]
    chk("zpen_25 の許容は 0.5 mm 固定", unit == "m" and floor_ == 0.5e-3 and rel_ == 0.0)

    # 平均 455 K の量が末尾窓で 6 K 動く反例 (codex M8) が超過になること
    tol = max(cs.ABS_TOL["dT_mid"][1], cs.ABS_TOL["dT_mid"][2] * 455.0)
    chk("平均 455 K で 6 K のうねりは許容超過", 6.0 > tol, "許容 %.3g K" % tol)

    # --- #50: 壁温不連続の run では収支超過を FAIL にしない (が CAVEAT に必ず出す) ---
    import inspect
    src = inspect.getsource(g.main)
    chk("壁温不連続なら収支超過を fails に積まない",
        "if not ok5b and nonuni:" in src and 'caveats.append("エネルギー収支' in src)
    chk("壁温不連続でなければ従来どおり fails に積む",
        src.count('fails.append("エネルギー収支")') == 1)
    chk("is_nonuniform_wall を本番から呼んでいる", "cs.is_nonuniform_wall(rd)" in src)

    # --- [4] 壁解像 (2026-10-15、plan tooling-nozzle-isothermal-wall-chain §5.1 #36s): 受け渡しの関門と同じ規則 ---
    wall_gate_cases()
    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def wall_gate_cases():
    """`wall_gate` を合成の run (solver_density_cuda/tools/test_gate_bad_input.make_wall_run: 2 枚の no-slip 壁 + 出口の出力、
    forge なし) で: 全 no-slip 壁・最終の主スナップショットの step・面積の重みで判定し、`--yplus-frac` を渡す。2026-10-15 までの
    呼び方 (群・step・重みなし) と合否が変わる入力も並べる。"""
    import shutil
    import subprocess
    import tempfile
    tools = Path(__file__).resolve().parents[3] / "solver_density_cuda" / "tools"
    sys.path.insert(0, str(tools))
    import test_gate_bad_input as tgb   # noqa: E402
    C = tgb.wall_cases()
    S, OLD = 2000, 1000
    both = {"wall_a": [OLD, S], "wall_b": [OLD, S]}
    root = tempfile.mkdtemp(prefix="case49_wall_gate_")

    def mk(name, walls, steps, outlet=True, snaps=(OLD, S)):
        d = tgb.make_wall_run(str(Path(root) / name), walls, steps, outlet=outlet)
        import h5py
        for st in snaps:                                        # 主スナップショット (vis_lam の無い空の HDF5; 粘性は config の定数)
            h5py.File(str(Path(d) / ("res_%d.h5" % st)), "w").close()
        return d

    def old_way(d, frac=2.0):
        r = subprocess.run([sys.executable, str(tools / "check_wall_resolution.py"), d, "--over-frac", str(frac)],
                           capture_output=True, text=True)
        return r.returncode
    try:
        d = mk("control", C["control"], both)
        rc, out, w, st = g.wall_gate(d, 2.0)
        chk("[4] 対照: 両壁とも最終の step に出力 -> PASS (全壁・最終の step・面積)",
            rc == 0 and st == S and sorted(w["walls"]) == ["wall_a", "wall_b"] and w["weight"] == "area"
            and "--weight" in w["cmd"] and str(S) in w["cmd"], "rc=%d step=%s %s" % (rc, st, w.get("cmd", [])[2:]))
        d = mk("b_old_only", C["control"], {"wall_a": [OLD, S], "wall_b": [OLD]}, outlet=False)
        rc, out, w, st = g.wall_gate(d, 2.0)
        chk("[4] wall_b の最終の step の出力だけ欠く -> 判定不能 (旧い呼び方は最新に戻って PASS)",
            rc == 2 and "wall_b" in w.get("why", "") and old_way(d) == 0, "rc=%d 旧 rc=%d %s" % (rc, old_way(d), w.get("why", "")[:80]))
        d = mk("inv_area_fail", C["inv_area_fail"], both, outlet=False)
        rc, out, w, st = g.wall_gate(d, 2.0)
        chk("[4] 面積の割合 約 48 % で FAIL (旧い呼び方は点数の割合 1/102 で PASS)", rc == 1 and old_way(d) == 0,
            "rc=%d 旧 rc=%d" % (rc, old_way(d)))
        rc, out, w, st = g.wall_gate(d, 60.0)
        chk("[4] --yplus-frac 60 を渡す -> 同じ入力が PASS (許容の受け渡し)", rc == 0 and "60.0" in w["cmd"], "rc=%d" % rc)
        d = mk("late_snapshot", C["control"], both, snaps=(OLD, S, 3000))
        rc, out, w, st = g.wall_gate(d, 2.0)
        chk("[4] 最終の主スナップショット 3000 に壁の出力が無い -> 判定不能", rc == 2 and st == 3000, "rc=%d step=%s" % (rc, st))
        d = mk("no_snaps", C["control"], both, snaps=())
        rc, out, w, st = g.wall_gate(d, 2.0)
        chk("[4] 主スナップショットが無い -> 判定不能", rc == 2 and st is None, "rc=%d" % rc)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
