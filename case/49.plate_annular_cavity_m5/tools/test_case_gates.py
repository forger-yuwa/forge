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

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
