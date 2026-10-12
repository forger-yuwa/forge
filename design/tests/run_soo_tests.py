#!/usr/bin/env python3
"""③ベル (D) の長さ固定の単目的の再最適化 (`opt.driver.Campaign.run_soo`・`soo_stop`) の試験 — plan tooling-nozzle-bell-two-mode §5.1 #8b。

CFD は回さない: `Campaign.evaluate` を解析関数 (2 変数の滑らかな山 + 一部の失敗) に差し替える。合格条件:
(1) `soo_stop` が停止候補 (DOE の後の PASS ≥ window かつ b_K − b_{K−window} < tol) と予算打切りを、最初に満たした位置で返す、
(2) `run_soo` が DOE → EI → 停止 → 確認 n_confirm 評価 → 最良形の回し直しの順に回り、要約 (stop・G・confirm・repeat) を書く、
固定の L は全評価で値のまま、評価済み点と重複しない、(3) 途中の台帳から再開しても同じ停止位置・同じ要約になる、
(4) MOO の経路 (`run`) は固定の変数を持つ問題を拒否する。
実行: design/.venv-opt の Python で `python design/tests/run_soo_tests.py`。
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.opt import driver as dr  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PROB = ROOT / "case" / "68.bell_two_mode" / "problem_bell_tp_d_fixl80.yaml"
FAIL = 0


def check(name, ok, detail=""):
    global FAIL
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL += 1


def fake_eta(x):
    th, te = x[0], x[1]
    return 0.977 - 2e-5 * (th - 31.0) ** 2 - 4e-5 * (te - 9.0) ** 2


def make(campaign_dir, fail_tags=()):
    c = dr.Campaign(PROB, campaign_dir, None, seed=0)
    c.feasible = lambda x: True
    L = float(c.bounds[2, 0])

    def ev(x, tag):
        x = [float(v) for v in np.asarray(x, dtype=float)]
        assert abs(x[2] - L) == 0.0, "固定の L が動いた"
        ok = tag not in fail_tags
        row = {"tag": tag, "x": x, "status": "PASS" if ok else "FAIL", "eta_cf": fake_eta(x) if ok else None,
               "fail_class": None if ok else "RETRYABLE", "wall_s": 1.0}
        c.rows.append(row)
        with open(c.ledger, "a") as f:
            f.write(json.dumps(row) + "\n")
        return row
    c.evaluate = ev
    return c


def main():
    # (1) soo_stop
    def H(etas, statuses=None):
        out, best, K = [], None, 0
        for i, e in enumerate(etas):
            st = (statuses or ["PASS"] * len(etas))[i]
            if st == "PASS":
                K += 1; best = e if best is None else max(best, e)
            out.append({"n": i + 1, "tag": f"t{i}", "status": st, "eta_cf": e, "K": K, "best": best})
        return out
    # PASS の K 点目までの最良 b_K: 0.95 になるのは K = 5、0.96 は K = 6 から。b_K − b_{K−8} < 1e-4 が初めて成り立つのは K = 14 (b_6 = 0.96)
    h = H([0.90] * 3 + [0.91, 0.95] + [0.96] * 9)
    s = dr.soo_stop(h, 3, 8, 1e-4, 40)
    check("soo_stop: 停止候補の位置", s is not None and s["reason"] == "stop_candidate" and s["n"] == 14, str(s))
    s2 = dr.soo_stop(h[:13], 3, 8, 1e-4, 40)
    check("soo_stop: 条件未成立なら None", s2 is None, str(s2))
    s3 = dr.soo_stop(H([0.9 + 0.001 * i for i in range(12)]), 3, 8, 1e-4, 10)
    check("soo_stop: 予算打切り", s3 is not None and s3["reason"] == "budget_cut" and s3["n"] == 10, str(s3))
    hf = H([0.90] * 3 + [0.96] * 10, ["PASS"] * 3 + ["PASS", "FAIL"] + ["PASS"] * 8)
    s4 = dr.soo_stop(hf, 3, 8, 1e-4, 40)
    check("soo_stop: 失敗は PASS の数に入れない", s4 is not None and s4["n"] == 13 and s4["K"] == 12, str(s4))
    # (2) run_soo
    td = Path(tempfile.mkdtemp())
    try:
        c = make(td / "camp", fail_tags=("ei_002",))
        summ = c.run_soo(n_doe=12, n_max=40, n_confirm=8, window=8, tol=1e-4)
        tags = [r["tag"] for r in c.rows]
        st = summ["stop"]
        n_conf = sum(1 for t in tags if t.startswith("conf_"))
        check("run_soo: DOE 12 点", sum(1 for t in tags if t.startswith("doe_")) == 12)
        check("run_soo: 停止 → 確認 8 評価 → 回し直し", st is not None and n_conf == 8 and tags[-1] == "rep_best",
              f"(stop {st and st['reason']} n {st and st['n']}, conf {n_conf}, 末尾 {tags[-1]})")
        X = np.array([r["x"][:2] for r in c.rows if not r["tag"].startswith("rep_")])
        lo = c.bounds[:2, 0]; span = c.bounds[:2, 1] - c.bounds[:2, 0]
        Xn = (X - lo) / span
        dmin = min(np.linalg.norm(Xn[i] - Xn[j]) for i in range(len(Xn)) for j in range(i))
        check("run_soo: 評価点が重複しない", dmin >= 5e-3, f"(最小距離 {dmin:.2e})")
        check("run_soo: 最良が山の頂上に近い", summ["best"]["eta_cf"] > 0.977 - 2e-4, f"(最良 η {summ['best']['eta_cf']:.6f})")
        check("run_soo: 要約の鍵", all(k in summ for k in ("stop", "G", "confirm", "repeat", "history")) and
              (td / "camp" / "soo_summary.json").exists())
        # (3) 再開: 台帳を確認の途中で切って回し直す
        rows = [json.loads(l) for l in (td / "camp" / "ledger.jsonl").read_text().splitlines()]
        k = next(i for i, r in enumerate(rows) if r["tag"] == "conf_003")
        (td / "camp2").mkdir()
        (td / "camp2" / "ledger.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows[:k]))
        c2 = make(td / "camp2", fail_tags=("ei_002",))
        summ2 = c2.run_soo(n_doe=12, n_max=40, n_confirm=8, window=8, tol=1e-4)
        same = (summ2["stop"] == summ["stop"] and [r["tag"] for r in c2.rows] == tags and
                np.allclose([r["x"] for r in c2.rows], [r["x"] for r in c.rows]))
        check("run_soo: 再開しても同じ停止位置・同じ評価列", same, f"(stop {summ2['stop']})")
        # (4) MOO の経路は固定の変数を拒否
        try:
            make(td / "camp3").run(4, 1, 1)
            check("MOO の経路: 固定の変数を拒否", False)
        except ValueError:
            check("MOO の経路: 固定の変数を拒否", True)
    finally:
        shutil.rmtree(td, ignore_errors=True)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
