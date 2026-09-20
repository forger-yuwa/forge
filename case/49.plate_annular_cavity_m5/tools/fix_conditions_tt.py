#!/usr/bin/env python3
"""`conditions.json` の Tt_tp 等を**修正版 setup.py で再計算**して差し替える (2026-09-20)。

`setup.py` の `_solve_T_of_h` がブラケット未検査で、目標エンタルピーが探索上限を超えると
**上限 (3000 K) をそのまま返していた**。M9 の主流全温は 3165.65 K なので、既存 run の
`conditions.json` に 3000.0 が焼き込まれている。`Taw_tp` は範囲内なので無傷。

元の値は `<key>_before_fix` として残す (履歴を消さない)。

usage: python3 tools/fix_conditions_tt.py RUN [RUN ...]   /   --all
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(CASE))
import setup as S  # noqa: E402

KEYS = ("Tt_tp", "cp_tp_at_Tt", "gamma_tp_at_Tt", "cp_rise_pct")


def fix(run):
    p = Path(run) / "conditions.json"
    if not p.exists():
        return "conditions.json なし"
    D = json.loads(p.read_text())
    if D.get("gas_used") != "TP":
        return "TP でない -> skip"
    if "Tt_tp_before_fix" in D:
        return "修正済み"
    g = S.tp_gas(D["dry_air_Y"])
    if g is None:
        return "semi-perfect gas が使えない"
    h0 = S._f(g.h_mass(D["T_inf"])) + 0.5 * D["U_inf"] ** 2
    Tt = S._solve_T_of_h(g, h0)
    old = D.get("Tt_tp")
    if old is not None and abs(old - Tt) < 1e-6:
        return "変更なし (%.2f K)" % Tt
    for k in KEYS:
        if k in D:
            D[k + "_before_fix"] = D[k]
    D["Tt_tp"] = Tt
    cp_tt = S._f(g.cp_mass(Tt))
    D["cp_tp_at_Tt"] = cp_tt
    D["gamma_tp_at_Tt"] = cp_tt / (cp_tt - D["R_tp"])
    D["cp_rise_pct"] = 100.0 * (cp_tt / D["cp_tp_298"] - 1.0)
    D["_Tt_tp_fix"] = ("2026-09-20: _solve_T_of_h のブラケット未検査で上限 %.1f K に"
                       "張り付いていたのを再計算" % (old if old else float("nan")))
    p.write_text(json.dumps(D, indent=2, ensure_ascii=False))
    return "%.2f K -> %.2f K" % (old, Tt)


if __name__ == "__main__":
    args = sys.argv[1:]
    runs = (sorted(str(d) for d in CASE.glob("run_*") if (d / "conditions.json").exists())
            if args == ["--all"] else args)
    for r in runs:
        print("%-46s %s" % (Path(r).name, fix(r)))
