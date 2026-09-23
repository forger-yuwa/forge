#!/usr/bin/env python3
r"""ソルバ内 CHT の**界面収束ゲート (G-if)** を判定する。

plan boundary-conjugate-heat-transfer §6 G-if / V4b(d)。`check_convergence.py` は流体の
保存量残差しか見ないので**別ツールにする** (列も意味も違う)。本ケース系列は流体残差が
`NOT CONVERGED (stalled/plateau)` のまま界面が静定するので、**界面はこちらで判定する**。

入力は run 直下の `conjugate_history.csv` (ソルバが更新ごとに追記する):

    step, physID, n, Tw_mean, Tw_min, Tw_max, dTw_max, res_abs_Wm2, res_max_W, res_rel, q_total[, res_solid]

判定 (**3 つを別々に**、`--n-consec` 回連続で満たすこと):

    ① 未緩和の界面残差 (局所面積で規格化)   res_abs_Wm2 <= eps_abs   [W/m2]
    ② 同 (最大荷重で規格化)                res_rel     <= eps_rel
    ③ **絶対**温度更新量                   dTw_max     <= dT_K      [K]
    (固体内部残差 res_solid があれば ④ res_solid <= tol_solid も要求する)

**許容は事前に登録する**。`--gate-json` でソルバが書いた登録値を読むか、CLI で渡す。
結果を見てから閾値を決めないこと (§6 の合格条件は結果より前に plan に書く)。

**判定不能は合格ではない** (memory: gate tools hardened): 列が無い / 行が無い /
値が非有限 / physID が混ざっているのに `--phys-id` が無い、は `REFUSED` を返す。

使い方:
  python3 solver_density_cuda/tools/check_cht_interface.py <run_dir> \
      --eps-rel 1e-3 --eps-abs 150 --dt-k 1e-2 --n-consec 80
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

REQUIRED = ["step", "physID", "dTw_max", "res_abs_Wm2", "res_rel"]


def fail(run: Path, verdict: str, lines: list[str]) -> int:
    text = "\n".join([f"=== G-if: {run.name} -> {verdict} ==="] + lines + [f"\nVERDICT: {verdict}"])
    print(text)
    (run / "CHT_INTERFACE_VERDICT.txt").write_text(text + "\n")
    return 0 if verdict == "PASS" else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--phys-id", type=int, default=None)
    ap.add_argument("--gate-json", default=None,
                    help="ソルバが起動時に書いた登録値 (eps_rel/eps_abs_Wm2/dT_K/n_consec)")
    ap.add_argument("--eps-rel", type=float, default=None)
    ap.add_argument("--eps-abs", type=float, default=None, help="[W/m2]")
    ap.add_argument("--dt-k", type=float, default=None, help="[K]")
    ap.add_argument("--tol-solid", type=float, default=None, help="固体内部残差 [W/m]")
    ap.add_argument("--n-consec", type=int, default=None)
    a = ap.parse_args()

    run = Path(a.run_dir)
    hist = run / "conjugate_history.csv"
    if not hist.exists():
        return fail(run, "REFUSED", [f"  {hist} が無い (ソルバ内連成を回した run を渡すこと)"])

    gate = {}
    gj = Path(a.gate_json) if a.gate_json else (run / "conjugate_gate.json")
    if gj.exists():
        gate = json.loads(gj.read_text())
    eps_rel = a.eps_rel if a.eps_rel is not None else gate.get("eps_rel")
    eps_abs = a.eps_abs if a.eps_abs is not None else gate.get("eps_abs_Wm2")
    dt_k = a.dt_k if a.dt_k is not None else gate.get("dT_K")
    n_consec = a.n_consec if a.n_consec is not None else gate.get("n_consec")
    tol_solid = a.tol_solid if a.tol_solid is not None else gate.get("tol_solid")
    missing = [k for k, v in (("eps_rel", eps_rel), ("eps_abs", eps_abs),
                              ("dT_K", dt_k), ("n_consec", n_consec)) if v is None]
    if missing:
        return fail(run, "REFUSED",
                    [f"  許容が登録されていない: {', '.join(missing)}",
                     "  **結果を見てから決めないこと**。plan §6 に書いた値を --gate-json か CLI で渡す。"])

    mode = gate.get("mode")
    n_avg = int(gate.get("flux_avg", 1) or 1)

    rows = list(csv.DictReader(open(hist)))
    if not rows:
        return fail(run, "REFUSED", ["  conjugate_history.csv が空 (更新が 1 回も起きていない)"])
    for k in REQUIRED:
        if k not in rows[0]:
            return fail(run, "REFUSED", [f"  列 {k} が無い (古い run か、書式が変わった)"])

    ids = sorted({int(r["physID"]) for r in rows})
    if len(ids) > 1 and a.phys_id is None:
        return fail(run, "REFUSED", [f"  physID が複数ある {ids} — --phys-id で選ぶこと"])
    pid = a.phys_id if a.phys_id is not None else ids[0]
    rows = [r for r in rows if int(r["physID"]) == pid]
    if len(rows) < n_consec:
        return fail(run, "NOT CONVERGED",
                    [f"  physID {pid}: 更新が {len(rows)} 回しかない (n_consec={n_consec} に満たない)"])

    def col(r, k):
        try:
            return float(r[k])
        except (TypeError, ValueError):
            return math.nan

    # **判定は「実際に最後の n_consec 更新」で行う** (codex result 2 巡目 M2)。
    # 条件を満たす行だけを拾って並べ直すと、**末尾の不良データを捨てて過去の良い窓で PASS** にできる
    # (末尾に「バッファ再初期化・残差 1e9」を足しても PASS した)。窓は動かさず、
    # 窓の全行が待機条件を満たしていなければ不合格にする。
    if n_avg > 1:
        if "update" not in rows[0] or "n_filled" not in rows[0]:
            return fail(run, "REFUSED",
                        ["  flux_avg > 1 なのに履歴に update / n_filled 列が無い",
                         "  (古い run。平均窓の待機を検査できないので合格にしない)"])
        # バッファが満ちた (= n_filled が n_avg に達した) 直近の時点を追う。
        # n_filled が減ったら再初期化なので、そこから数え直す。
        fill_at, prev_filled = None, -1
        for r in rows:
            nf = int(float(r["n_filled"]))
            if nf < prev_filled:            # 再初期化
                fill_at = None
            if fill_at is None and nf >= n_avg:
                fill_at = int(float(r["update"]))
            prev_filled = nf
            r["_ready"] = (fill_at is not None
                           and int(float(r["update"])) - fill_at >= 2 * n_avg
                           and nf >= n_avg)
        tail_chk = rows[-n_consec:]
        if not all(r["_ready"] for r in tail_chk):
            nbad = sum(1 for r in tail_chk if not r["_ready"])
            return fail(run, "NOT CONVERGED",
                        [f"  flux_avg={n_avg}: 末尾 {n_consec} 更新のうち {nbad} 行が "
                         f"「充填後 2N={2*n_avg} 更新」の待機を満たしていない",
                         "  (窓は動かさない。条件を満たす行だけ拾うと末尾の不良を捨てて合格にできる)"])

    # **`fem2d` では固体内部残差を必須にする** (codex result 2026-09-23 M2)。
    # 列も許容も無いまま「検査を省略して PASS」を出さない。
    if (mode == "fem2d") or ("res_solid" in rows[0]):
        if "res_solid" not in rows[0]:
            return fail(run, "REFUSED", ["  mode=fem2d なのに履歴に res_solid 列が無い"])
        if tol_solid is None:
            return fail(run, "REFUSED",
                        ["  固体内部残差の許容 (tol_solid) が登録されていない",
                         "  conjugate.gate に tol_solid を書くか --tol-solid で渡すこと"])

    tail = rows[-n_consec:]
    vals = {k: [col(r, k) for r in tail] for k in ("dTw_max", "res_abs_Wm2", "res_rel")}
    if "res_solid" in rows[0]:
        vals["res_solid"] = [col(r, "res_solid") for r in tail]
    if any(not math.isfinite(v) for vs in vals.values() for v in vs):
        return fail(run, "REFUSED", ["  末尾の窓に非有限値がある (NaN/Inf)"])

    checks = [
        ("① res_abs_Wm2 [W/m2]", max(vals["res_abs_Wm2"]), eps_abs),
        ("② res_rel", max(vals["res_rel"]), eps_rel),
        ("③ dTw_max [K]", max(vals["dTw_max"]), dt_k),
    ]
    if tol_solid is not None and "res_solid" in vals:
        checks.append(("④ res_solid [W/m]", max(vals["res_solid"]), tol_solid))

    lines = [f"  physID {pid}: 更新 {len(rows)} 回、末尾 {n_consec} 回で判定 "
             f"(step {tail[0]['step']} .. {tail[-1]['step']})"]
    ok = True
    for name, got, tol in checks:
        good = got <= tol
        ok &= good
        lines.append(f"  {'OK  ' if good else 'NG  '}{name:<22} max {got:.4e}  <= {tol:.4e}")

    # 発散していないか (末尾窓で dTw が単調に増えている)
    d = vals["dTw_max"]
    if len(d) >= 4 and all(d[i] <= d[i+1] for i in range(len(d)-1)) and d[-1] > 2.0 * max(1e-300, d[0]):
        lines.append(f"  NG  dTw_max が末尾窓で単調増加 ({d[0]:.3e} -> {d[-1]:.3e} K)")
        return fail(run, "DIVERGED", lines)

    if "Tw_mean" in rows[0]:
        tw = [col(r, "Tw_mean") for r in tail]
        lines.append(f"  Tw_mean 末尾 {tw[-1]:.4f} K (窓内の振れ {max(tw)-min(tw):.2e} K)")
    if "q_total" in rows[0]:
        q = [col(r, "q_total") for r in tail]
        lines.append(f"  q_total 末尾 {q[-1]:.4f} W/m (窓内の振れ {max(q)-min(q):.2e})")

    return fail(run, "PASS" if ok else "NOT CONVERGED", lines)


if __name__ == "__main__":
    sys.exit(main())
