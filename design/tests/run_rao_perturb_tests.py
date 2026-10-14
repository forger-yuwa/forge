#!/usr/bin/env python3
"""③ベル (M) の摂動 (`geometry.wall.RaoPerturbation`・`geometry.rao_perturb`) と摂動探索の配線の試験 — plan tooling-nozzle-bell-two-mode §5.1 #8c。

合格条件: (1) η_k の端の条件 η(0) = η′(0) = η(1) = 0・η′(1) = kπ(−1)^k、1・2 階が数値微分と一致、(2) runner で a = 0 の壁が補正壁
(`problem_bell_tp_m_ns_l100_corr.yaml`) とビット一致、a ≠ 0 でも出口半径 (ε 50)・スロート・B の r と r′ は不変、出口の角は変わる、
validate が空、(3) driver の設計変数の並びが pert_a1..4 + L、L は固定、(4) `run_soo(with_start=True)` が出発点 (a = 0) を最初に評価する。
実行: design/.venv-opt の Python で `python design/tests/run_rao_perturb_tests.py`。
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.evaluate.runner import _build_wall_rao  # noqa: E402
from forge_design.geometry.wall import RaoPerturbation  # noqa: E402
from forge_design.opt import driver as dr  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
C68 = ROOT / "case" / "68.bell_two_mode"
PERT = C68 / "problem_bell_tp_m_ns_l100_corr_pert.yaml"
CORR = C68 / "problem_bell_tp_m_ns_l100_corr.yaml"
FAIL = 0


def check(name, ok, detail=""):
    global FAIL
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL += 1


def with_coeffs(a, td):
    raw = yaml.safe_load(PERT.read_text())
    for k, v in enumerate(a, start=1):
        raw["dv"][f"pert_a{k}"]["value"] = float(v)
    raw = dr.absolutize_problem_paths(raw, C68)
    p = Path(td) / "p.yaml"
    p.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    return load_problem(p)


def main():
    # (1) 端の条件と微分
    for k in range(1, 5):
        pk = RaoPerturbation([0.0] * (k - 1) + [1.0], 0.0, 1.0)
        v0, d0, v1, d1 = (float(pk(np.array([1e-12]))[0]), float(pk(np.array([1e-12]), 1)[0]),
                          float(pk(np.array([1.0]))[0]), float(pk(np.array([1.0]), 1)[0]))
        check(f"η_{k} の端", abs(v0) < 1e-12 and abs(d0) < 1e-9 and abs(v1) < 1e-12 and abs(d1 - k * np.pi * (-1) ** k) < 1e-9,
              f"(η(0) {v0:.1e}, η′(0) {d0:.1e}, η(1) {v1:.1e}, η′(1) {d1:+.4f} 期待 {k*np.pi*(-1)**k:+.4f})")
    pr = RaoPerturbation([0.01, -0.02, 0.015, -0.005], 0.5, 22.6)
    xs = np.linspace(0.6, 22.5, 300); h = 1e-5
    e1 = np.max(np.abs(pr(xs, 1) - (pr(xs + h) - pr(xs - h)) / (2 * h)))
    e2 = np.max(np.abs(pr(xs, 2) - (pr(xs + h) - 2 * pr(xs) + pr(xs - h)) / h ** 2))
    check("1・2 階が数値微分と一致", e1 < 1e-8 and e2 < 1e-4, f"({e1:.1e}, {e2:.1e})")
    check("区間の外は 0", float(np.max(np.abs(pr(np.array([0.0, 0.3, 0.5]), 0))) ) == 0.0)
    td = tempfile.mkdtemp()
    try:
        # (2) runner
        w0, _ = _build_wall_rao(with_coeffs([0, 0, 0, 0], td))
        wc, _ = _build_wall_rao(load_problem(CORR))
        xg = np.linspace(w0.x_in, w0.x_e, 4001)
        check("a = 0 で補正壁とビット一致", all(np.array_equal(w0.r(xg, d), wc.r(xg, d)) for d in (0, 1, 2)))
        w1, rec = _build_wall_rao(with_coeffs([0.03, -0.03, 0.03, -0.03], td))
        xe = np.array([w1.x_e]); xB = np.array([w1.x_B]); z = np.array([0.0])
        same = (abs(float(w1.r(xe)[0]) - float(wc.r(xe)[0])) < 1e-12 and
                all(abs(float(w1.r(z, d)[0]) - float(wc.r(z, d)[0])) < 1e-15 for d in (0, 1, 2)) and
                all(abs(float(w1.r(xB, d)[0]) - float(wc.r(xB, d)[0])) < 1e-12 for d in (0, 1)))
        th_e = (np.degrees(np.arctan(float(w1.r(xe, 1)[0]))), np.degrees(np.arctan(float(wc.r(xe, 1)[0]))))
        check("a ≠ 0: 出口半径・スロート・B の r と r′ は不変", same, f"(出口 r {float(w1.r(xe)[0]):.9f})")
        check("a ≠ 0: 出口の角が変わる", abs(th_e[0] - th_e[1]) > 1e-3, f"({th_e[1]:.3f}° → {th_e[0]:.3f}°)")
        check("a ≠ 0: validate 空・記録", w1.validate() == [] and rec.get("perturb", {}).get("coeffs") == [0.03, -0.03, 0.03, -0.03])
        # (3) driver の並び
        c = dr.Campaign(PERT, Path(td) / "camp", None, seed=0)
        check("設計変数の並び", c.dv_order == ("pert_a1", "pert_a2", "pert_a3", "pert_a4", "L_over_rt") and
              list(c.fixed_mask) == [False, False, False, False, True], str(c.dv_order))
        corners = [c.feasible(np.r_[np.array(s) * 0.03, c.bounds[4, 0]]) for s in ([1, 1, 1, 1], [-1, -1, -1, -1], [1, -1, 1, -1])]
        check("範囲の角で幾何が成立", all(corners), str(corners))
        # (4) 出発点
        seen = []

        def ev(x, tag):
            x = [float(v) for v in np.asarray(x, dtype=float)]
            seen.append((tag, x))
            row = {"tag": tag, "x": x, "status": "PASS", "eta_cf": 0.98 - 0.01 * sum(v * v for v in x[:4]), "wall_s": 1.0}
            c.rows.append(row)
            with open(c.ledger, "a") as f:
                f.write(json.dumps(row) + "\n")
            return row
        c.evaluate = ev
        c.run_soo(n_doe=6, n_max=10, n_confirm=1, window=3, tol=1e-12, repeat_best=False, with_start=True)
        check("出発点を最初に評価", seen[0][0] == "doe_start" and seen[0][1][:4] == [0.0, 0.0, 0.0, 0.0], str(seen[0]))
    finally:
        shutil.rmtree(td, ignore_errors=True)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
