#!/usr/bin/env python3
r"""`settle_virtual_stop.py` の試験 (plan tooling-settle-convergence-and-monitor §6「V3(b)・V4 の実施」(1))。合成の run だけを使う。

    python3 solver_density_cuda/tools/test_settle_virtual_stop.py [--quick] [--bench]
    (design/.venv-opt の python でも system python3 でも; 作業ツリーの直下でなくてよい)

合成の run は `test_settle.py` の作り方 (`build`・`base_criteria`; 構造格子 30 × 10、節点・壁・流量・残差) を import して使う。
forge と同じ形 (`start_obs=False`: 格子点 0 に観測が無い) にするので、最初に判定できる格子点は 42,000、保持 40,000 の
停止候補は一定の系列で 82,000。

  (i)   一定の系列 → 結果 A (t_s 82,000、span 82,000〜162,000 の差はすべて 0)。
  (ii)  codex の人工の例 (2 万 step あたり ΔM = 0.0002、許容 0.000275): 42,000〜200,000 の 80 窓すべて SETTLED だが、t_s からの
        累積が許容を超える → 結果 B (出口コア M と線の d; どの量・位置・step かを書く)。
  (iii) 停止候補が無い (出口 M が格子点ごとに振れて SETTLED にならない) → 未成立。
  (iv)  t_s が 120,000 を超える (100,000 まで動き、その後一定 → t_s 160,000) → 未成立。
  (v)   span の途中でデータが切れる (120,000 まで) → 未成立。追加: span の格子点の観測の欠け (判定不能が混じる) → 未成立。
  追加: 超過とデータの欠けが同時 → B (欠けは incomplete に書く; 結果を見る前に決めた順)、(b) だけの超過 (残差の跳ね) → B と
        追加の記録 (目的の量だけなら SETTLED なのに残差で落ちた格子点の数と step)、壁の 1 節点だけのドリフト → B (壁・節点・step)、
        基準のファイルが無い → 未成立、CLI の終了コードと出力のファイル。
  参考の t_s′ (目的の量だけで保持; 判定に使わない): 一定の系列で t_s′ = t_s、(iii) で t_s′ なし、
        (y1) 残差の跳ね (行 50,010) で t_s が 112,000 に遅れ t_s′ は 82,000 (本検証 A・t_s′ からも超過なし)、
        (y2) 同じ跳ねに出口 M の 82,000〜122,000 のゆっくりした変化を重ねる → 本検証は A のまま、t_s′ からは超過あり
        (残差の条件が早すぎる停止を防いだ形)。
  (vi)  保持の数え方が停止規則 v2 と同じ:
        - 写した `hold_step` と `stop_rule._hold_step` が乱数の系列で同じ到達・同じ状態。
        - 本物の v2 の見張り (`watch_main_stage_v2`) を `design/tests/run_stop_rule_v2_tests.py` の偽の forge と run の作り方で回し、
          v2 が置いた要求 (種類・step・条件を満たし始めた step・保持の経過) と、同じ run dir を再生して求めた t_s・要求が一致する:
          保持 38,000 (要求 80,000)、保持の既定 40,000 (82,000)、保持の途中で崩れてやり直す (60,001 で出口 M が跳ぶ → 82,000 から
          数え直して 120,000)、DIVERGED (52,000 で diverged の要求 → 停止候補なし)、未設定 (42,000 で criteria_unset → 停止候補なし)。
        --quick は偽の forge の部分を省く (約 1 分)。
  --bench: 20 万 step・格子 2,000 (格子点 100) の規模の再生の所要を測る (残差の行 60 万 = 3 行/step・14 列、節点 6,000、壁 750 節点;
        case/71 の run_0005 の残差の履歴 600,001 行・120 MB と節点の表 5,335 点に合わせた)。試験の合否には入れない。
"""
from __future__ import annotations

import io
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import redirect_stdout
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import settle  # noqa: E402
import settle_virtual_stop as SVS  # noqa: E402
import test_settle as TS  # noqa: E402  (合成の run の作り方; import で作業の一時ディレクトリ TS.TMP を作る)

QUICK = "--quick" in sys.argv[1:]
BENCH = "--bench" in sys.argv[1:]
FAILS = []
NCHK = 0
NI, NJ, MD, DT = TS.NI, TS.NJ, TS.MD, TS.DT


def chk(name, cond, detail=""):
    global NCHK
    NCHK += 1
    if not cond:
        FAILS.append(name)
    print(f"  [{'OK ' if cond else 'NG '}] {name}" + (f"  — {detail}" if detail else ""), flush=True)


def build(name, **kw):
    """test_settle.build を forge と同じ形 (格子点 0 に観測なし) で呼ぶ。"""
    kw.setdefault("start_obs", False)
    return TS.build(name, **kw)


def run_vs(run, **kw):
    return SVS.analyse(run, **kw)


def a_of(res, name):
    return res["span"]["a"][name]


def grid_rows(run):
    with open(run / SVS.OUT_DIR / "grid_verdicts.csv") as f:
        import csv
        return list(csv.DictReader(f))


# ------------------------------------------------------------------ (i)〜(v) と追加
def test_results():
    print("=== (i) 一定の系列 → 結果 A ===")
    run = build("vs_i_const", end=200000)
    res = run_vs(run)
    chk("(i) 一定の系列 → A、t_s 82,000 (条件を満たし始めた 42,000 + 保持 40,000)、span 82,000〜162,000",
        res["verdict"] == "A" and res["t_s"] == 82000 and res["settled_since"] == 42000 and res["span_end"] == 162000,
        f"{res['verdict']} t_s {res['t_s']} reasons {res['reasons'][:1]} B {res['b_evidence'][:1]} inc {res['incomplete'][:1]}")
    chk("(i) span の 41 格子点 (両端込み) で (a) 4 量すべて差 0・(b) すべて SETTLED",
        res["span"]["n_grid_points"] == 41 and all(a["status"] == "内" and a["n_points"] == 41 and a["max_ratio"] == 0.0
                                                   for a in res["span"]["a"].values()) and res["span"]["b_not_settled"] == [],
        json.dumps({k: (v["status"], v["n_points"], v["max_ratio"]) for k, v in res["span"]["a"].items()}, ensure_ascii=False))
    od = run / SVS.OUT_DIR
    files = ["virtual_stop.json", "VERDICT.txt", "grid_verdicts.csv", "grid_checks.csv", "span_deviation.csv",
             "span_deviation_ts_prime.csv", "deviation_from_ts.png", "grid_verdicts.png"]
    chk("(i) 出力のファイルがそろう (JSON・VERDICT・CSV 4 つ・図 2 枚)", all((od / f).is_file() for f in files),
        str([f for f in files if not (od / f).is_file()]))
    vt = (od / "VERDICT.txt").read_text().splitlines()
    js = json.loads((od / "virtual_stop.json").read_text())
    chk("(i) VERDICT.txt の 1 行目が `VERDICT: A`、JSON の verdict と t_s", vt[0] == "VERDICT: A" and js["verdict"] == "A"
        and js["t_s"] == 82000 and js["files"], vt[0])
    to = res["target_only"]
    chk("(i) 参考の t_s′ (目的の量だけで保持): t_s′ = t_s = 82,000、t_s′ からの検証は超過なし、判定に使わない印",
        to["t_s_prime"] == 82000 and to["relation"] == "t_s′ = t_s" and to["outcome"] == "超過なし" and to["used_for_verdict"] is False
        and js["target_only"]["t_s_prime"] == 82000, f"{to['relation']} {to['outcome']}")
    chk("(i) VERDICT.txt の 2 行目以降に参考の t_s′ の行 (t_s′ = t_s)",
        any(ln.startswith("参考 (判定に使わない)") and "t_s′ = t_s" in ln for ln in vt[1:]), str([ln for ln in vt if "参考" in ln][:1]))
    rows = grid_rows(run)
    chk("(i) grid_verdicts.csv: 42,000〜200,000 の 80 格子点、t_s の行に request:settled", len(rows) == 80
        and [r["step"] for r in rows if "request:settled" in r["v2_event"]] == ["82000"], str(len(rows)))
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = SVS.main([str(run), "--no-plot"])
    chk("(i) CLI: 終了コード 0 と最後の行 VERDICT: A", rc == 0 and buf.getvalue().splitlines()[-1] == "VERDICT: A", buf.getvalue().splitlines()[-1])

    print("=== (ii) codex の人工の例 → 結果 B ===")
    c = TS.base_criteria()
    c["quantities"][2]["thresholds"]["tol"] = 0.000275            # exit_core_mach の許容 (0.1 × 0.05 % × 5.5)
    run = build("vs_ii_drift", end=200000, criteria=c, M_fn=lambda s: np.full(NI, MD + 0.0002 * s / 20000))
    res = run_vs(run)
    rows = grid_rows(run)
    chk("(ii) 窓ごとの判定は 42,000〜200,000 の 80 格子点すべて SETTLED (窓の幅の判定は累積のドリフトを見ない)",
        len(rows) == 80 and all(r["verdict"] == "SETTLED" for r in rows), str(sorted({r["verdict"] for r in rows})))
    a = a_of(res, "exit_core_mach")
    chk("(ii) → B、t_s 82,000、出口コア M の差の最大は span の終わり 162,000 で 0.0008 (許容 0.000275 の 2.91 倍)",
        res["verdict"] == "B" and res["t_s"] == 82000 and a["status"] == "超過" and a["at_step"] == 162000
        and abs(a["diff"] - 0.0008) < 1e-12 and abs(a["max_ratio"] - 0.0008 / 0.000275) < 1e-9,
        f"{res['verdict']} {a['status']} step {a['at_step']} diff {a['diff']} ratio {a['max_ratio']}")
    chk("(ii) 最初の超過は 110,000 (差 0.00028 > 0.000275)、(b) はすべて SETTLED", a["first_exceeded"] == 110000
        and res["span"]["b_not_settled"] == [], f"first {a['first_exceeded']}")
    al = a_of(res, "line_mach")
    chk("(ii) 線の d も超過 (位置ごとの差 100 × 0.0008/3 = 0.0267 %pt > 0.01)、位置を書く", al["status"] == "超過"
        and abs(al["diff"] - 100 * 0.0008 / MD) < 1e-9 and "xq の位置" in al["where"], f"{al['diff']} {al['where']}")
    chk("(ii) 流れ角・壁は内 (差 0)", a_of(res, "exit_flow_angle")["status"] == "内" and a_of(res, "wall_heat_nodes")["status"] == "内")
    chk("(ii) B の根拠に量・step・位置が書かれ、VERDICT.txt の 1 行目が `VERDICT: B`",
        any("exit_core_mach" in b and "162000" in b and "比" in b for b in res["b_evidence"])
        and (run / SVS.OUT_DIR / "VERDICT.txt").read_text().splitlines()[0] == "VERDICT: B", res["b_evidence"][0][:160])

    print("=== (iii) 停止候補が無い → 未成立 ===")
    run = build("vs_iii_no_candidate", end=200000, M_fn=lambda s: np.full(NI, MD + (0.001 if (s // DT) % 2 else 0.0)))
    res = run_vs(run)
    chk("(iii) 出口 M が格子点ごとに 0.001 振れる → どの格子点も SETTLED にならず停止候補なし → 未成立",
        res["verdict"] == "未成立" and res["t_s"] is None and res["v2_emulation"]["request"] is None
        and any("停止候補が無い" in r for r in res["reasons"]), str(res["reasons"][:1]))
    to = res["target_only"]
    vt = (run / SVS.OUT_DIR / "VERDICT.txt").read_text().splitlines()
    chk("(iii) 参考の t_s′ も無い (目的の量が SETTLED にならない) → 「t_s′ なし」を JSON と VERDICT.txt に書く",
        to["t_s_prime"] is None and to["relation"].startswith("t_s′ なし") and to["outcome"] == "t_s′ なし"
        and any(ln.startswith("参考 (判定に使わない)") and "t_s′ なし" in ln for ln in vt), to["relation"])

    print("=== (iv) t_s > 120,000 → 未成立 ===")
    run = build("vs_iv_late", end=200000, M_fn=lambda s: np.full(NI, MD + 1e-6 * min(s, 100000)))
    res = run_vs(run)
    chk("(iv) 100,000 まで動き、その後一定 → 条件を満たし始めた 120,000、t_s 160,000 > 120,000 → 未成立",
        res["verdict"] == "未成立" and res["t_s"] == 160000 and res["settled_since"] == 120000
        and any("160000 > 120000" in r for r in res["reasons"]), str(res["reasons"][:1]))

    print("=== (v) span の途中でデータが切れる → 未成立 ===")
    run = build("vs_v_cut", end=120000)
    res = run_vs(run)
    chk("(v) 一定の系列が 120,000 で切れる (t_s 82,000、span の終わり 162,000) → 未成立 (span の終わりまでデータが無い)",
        res["verdict"] == "未成立" and res["t_s"] == 82000 and res["b_evidence"] == []
        and any("span の終わり" in x for x in res["incomplete"]), str(res["incomplete"][:1]))
    run = build("vs_v_gap", end=200000, skip=(100000,))
    res = run_vs(run)
    chk("(v′) span の格子点 100,000 の観測が無い (判定不能が混じる) → 未成立、(b) INDETERMINATE と (a) の値を作れないを記録",
        res["verdict"] == "未成立" and res["t_s"] == 82000 and any("(b) 格子点 100000 の全体の判定が INDETERMINATE" in x for x in res["incomplete"])
        and any(x.startswith("(a)") and "step 100000" in x for x in res["incomplete"]), str(res["incomplete"][:2]))

    print("=== 追加: 判定の順・(b) だけの超過・壁・基準なし ===")
    c = TS.base_criteria()
    c["quantities"][2]["thresholds"]["tol"] = 0.000275
    run = build("vs_x_drift_cut", end=150000, criteria=c, M_fn=lambda s: np.full(NI, MD + 0.0002 * s / 20000))
    res = run_vs(run)
    chk("(x1) 超過 (110,000〜) とデータの欠け (150,000 で切れる) が同時 → B (超過は反例)、欠けは incomplete に書く",
        res["verdict"] == "B" and any("span の終わり" in x for x in res["incomplete"]) and a_of(res, "exit_core_mach")["first_exceeded"] == 110000,
        f"{res['verdict']} {res['incomplete'][:1]}")

    def spike(s):
        v = dict(TS.RES_BASE)
        if s == 120010:
            v = {k: x * 1e8 for k, x in v.items()}
        return v
    run = build("vs_x_resid_spike", end=200000, resid_fn=spike)
    res = run_vs(run)
    bad = [b["step"] for b in res["span"]["b_not_settled"]]
    ex = res["extra_record"]
    chk("(x2) 残差の行 120,010 だけ 10⁸ 倍 → 終端 122,000〜140,000 の 10 格子点が NOT_SETTLED (跳ね) → B ((b) だけ; (a) はすべて内)",
        res["verdict"] == "B" and bad == list(range(122000, 140001, DT)) and all(v["status"] == "内" for v in res["span"]["a"].values()),
        f"{res['verdict']} {bad}")
    chk("(x2) 追加の記録: 目的の量だけなら SETTLED なのに残差で落ちた格子点 10 (すべて span)、その step",
        ex["failed"]["count"] == 10 and ex["failed"]["by_region"] == {"span": 10}
        and [x["step"] for x in ex["failed"]["steps"]] == list(range(122000, 140001, DT))
        and all(x["residual_mass"]["residual"] == "NOT_SETTLED" for x in ex["failed"]["steps"]), json.dumps(ex["failed"]["by_region"]))

    xm, rm = TS.grid_xr()
    coord = np.stack([xm[:, -1], rm[:, -1], np.zeros(NI)], axis=1)
    conne = []
    for i in range(NI - 1):
        conne += [2, 2, i, i + 1]
    A, _ = settle.wall_node_area(coord, np.asarray(conne), True)
    Q0 = 100.0 + np.arange(NI, dtype=float)
    tol10 = 1000.0 * A[10] + 1e-3 * abs(Q0[10])
    rate = 0.7 * tol10 / 20000.0

    def Qd(s):
        Q = Q0.copy()
        Q[10] += rate * s
        return Q
    run = build("vs_x_wall_drift", end=200000, Q_fn=Qd)
    res = run_vs(run)
    aw = a_of(res, "wall_heat_nodes")
    exp = rate * 80000 / (1000.0 * A[10] + 1e-3 * abs(Q0[10] + rate * 82000))   # 許容は Q_i(t_s) で決まる (約 2.79)
    chk("(x3) 壁の節点 10 だけ 2 万 step あたり許容の 0.7 倍のドリフト: 窓ごとは SETTLED、t_s からの累積 約 2.8 倍 → B (壁・節点 10・162,000)",
        res["verdict"] == "B" and all(r["verdict"] == "SETTLED" for r in grid_rows(run)) and aw["status"] == "超過"
        and aw["at_step"] == 162000 and "節点 10 " in aw["where"] and abs(aw["max_ratio"] - exp) < 1e-9,
        f"{res['verdict']} {aw['where']} {aw['max_ratio']} (期待 {exp})")

    print("=== 追加: 参考の t_s′ (目的の量だけで保持; 判定に使わない) ===")

    def spike50(s):
        v = dict(TS.RES_BASE)
        if s == 50010:
            v = {k: x * 1e8 for k, x in v.items()}
        return v
    run = build("vs_y1_ts_prime_early", end=200000, resid_fn=spike50)
    res = run_vs(run)
    to = res["target_only"]
    ex = res["extra_record"]
    chk("(y1) 残差の行 50,010 の跳ねで 52,000〜70,000 が NOT_SETTLED → 保持がやり直され t_s 112,000、目的の量だけなら t_s′ 82,000",
        res["t_s"] == 112000 and res["settled_since"] == 72000 and to["t_s_prime"] == 82000 and to["settled_since"] == 42000
        and to["relation"] == "t_s′ < t_s (30000 step 早い)", f"t_s {res['t_s']} t_s′ {to['t_s_prime']} {to['relation']}")
    chk("(y1) 一定の系列なので本検証は A、t_s′ からの検証も超過なし (記録だけ)、追加の記録は t_s の前の 10 格子点",
        res["verdict"] == "A" and to["outcome"] == "超過なし" and ex["failed"]["by_region"] == {"before_ts": 10},
        f"{res['verdict']} {to['outcome']} {ex['failed']['by_region']}")
    rows = grid_rows(run)
    r82 = [r for r in rows if r["step"] == "82000"][0]
    chk("(y1) grid_verdicts.csv: 82,000 の行で目的の量だけの事象 request:settled・全体は保持の途中・t_s′ からの区間 span",
        "request:settled" in r82["v2_event_target_only"] and r82["v2_event"] == "" and r82["hold_start"] == "72000"
        and r82["target_verdict"] == "SETTLED" and r82["region_ts_prime"] == "span" and r82["region"] == "before_ts", str(r82)[:200])

    c = TS.base_criteria()
    c["quantities"][2]["thresholds"]["tol"] = 0.000275
    rate = 0.00052 / 40000                       # 2 万 step あたり 0.00026 (許容 0.000275 の内)、82,000〜122,000 だけ動く

    def ramp(s):
        return np.full(NI, MD + rate * min(max(s - 82000, 0), 40000))
    run = build("vs_y2_ts_prime_premature", end=200000, criteria=c, resid_fn=spike50, M_fn=ramp)
    res = run_vs(run)
    to = res["target_only"]
    ae = to["span"]["a"]["exit_core_mach"]
    chk("(y2) 出口 M が 82,000〜122,000 で 0.00052 動く (窓ごとは SETTLED): 本検証は t_s 112,000 からの差 0.00013 で A",
        res["verdict"] == "A" and res["t_s"] == 112000 and abs(a_of(res, "exit_core_mach")["diff"] - 0.00013) < 1e-12,
        f"{res['verdict']} t_s {res['t_s']} diff {a_of(res, 'exit_core_mach')['diff']}")
    chk("(y2) 目的の量だけなら t_s′ 82,000 で止まり、そこからは超過あり (出口コア M の最初の超過 104,000、最大 0.00052 / 0.000275)",
        to["t_s_prime"] == 82000 and to["outcome"] == "超過あり" and ae["status"] == "超過" and ae["first_exceeded"] == 104000
        and abs(ae["max_ratio"] - 0.00052 / 0.000275) < 1e-9 and to["span"]["b_not_settled"] == []
        and any("exit_core_mach" in x for x in to["evidence"]), f"{to['outcome']} first {ae['first_exceeded']} ratio {ae['max_ratio']}")
    tp = (run / SVS.OUT_DIR / "span_deviation_ts_prime.csv").read_text().splitlines()
    vt = (run / SVS.OUT_DIR / "VERDICT.txt").read_text().splitlines()
    chk("(y2) 判定は A のまま (t_s′ の超過は判定に使わない)、span_deviation_ts_prime.csv に t_s′ からの差、VERDICT.txt に参考の超過の行",
        vt[0] == "VERDICT: A" and any(ln.startswith("参考 t_s′ の超過") for ln in vt)
        and any(ln.startswith("104000,span,exit_core_mach,exit_core_mach,ok,") and ",True," in ln for ln in tp), vt[0])

    run = build("vs_x_no_criteria", end=200000, criteria=False)
    res = run_vs(run)
    chk("(x4) 基準のファイルが無い → 未成立 (評価対象が未設定)", res["verdict"] == "未成立" and "評価対象が未設定" in res["reasons"][0],
        res["reasons"][0][:80])
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = SVS.main([str(run), "--no-plot"])
    chk("(x4) CLI: 未成立は終了コード 3、最後の行 `VERDICT: 未成立`", rc == 3 and buf.getvalue().splitlines()[-1] == "VERDICT: 未成立")
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = SVS.main([str(build("vs_ii_drift_cli", end=200000, criteria=c, M_fn=lambda s: np.full(NI, MD + 0.0002 * s / 20000))),
                       "--no-plot"])
    chk("(x5) CLI: B は終了コード 1", rc == 1 and buf.getvalue().splitlines()[-1] == "VERDICT: B")


# ------------------------------------------------------------------ (vi) v2 と同じ保持
def test_hold_unit():
    print("=== (vi-a) 写した hold_step と stop_rule._hold_step ===")
    sys.path.insert(0, str(ROOT / "design"))
    from forge_design.evaluate import stop_rule as SR
    rng = random.Random(20261016)
    same = True
    detail = ""
    for trial in range(400):
        hold = rng.choice([0, 2000, 6000, 38000, 40000])
        seq = [(42000 + 2000 * k, rng.random() < 0.8) for k in range(rng.randint(1, 60))]
        a = {"start": None, "reached": [], "breaks": []}
        b = {"start": None, "reached": [], "breaks": []}
        for s, ok in seq:
            ra = SVS.hold_step(a, ok, s, hold, f"w{s}")
            rb = SR._hold_step(b, ok, s, hold, f"w{s}")
            if ra != rb or a != b:
                same = False
                detail = f"trial {trial} step {s}: {ra} {rb} {a} {b}"
                break
        if not same:
            break
    chk("(vi-a) 乱数の系列 400 本 (保持 0〜40,000) で到達の判定と状態 (start・reached・breaks) が毎格子点で同じ", same, detail)
    chk("(vi-a) 保持の既定が同じ (40,000)", SVS.HOLD_STEPS_DEFAULT == SR.HOLD_STEPS_DEFAULT)


def test_v2_end_to_end():
    print("=== (vi-b) 本物の v2 の見張り (偽の forge) と再生の t_s ===")
    sys.path.insert(0, str(ROOT / "design" / "tests"))
    t0 = time.perf_counter()
    import run_stop_rule_v2_tests as V   # 偽の forge・run の作り方 (import で FORGE_BIN を偽の forge にする)
    from forge_design.evaluate import runner as RN
    from forge_design.evaluate import stop_rule as SR
    if RN._ENV.get("FORGE_BIN") != str(V.FAKE_FORGE):   # runner が先に読み込まれていると本物の forge を起動してしまう
        chk("(vi-b) runner._ENV の FORGE_BIN が偽の forge (本物の forge を起動しない)", False, str(RN._ENV.get("FORGE_BIN")))
        return
    V.FAKE_FORGE.write_text(V.FAKE_SRC.format(py=sys.executable))
    V.FAKE_FORGE.chmod(0o755)
    real = SR.extract_delta_e
    SR.extract_delta_e = V.fake_extract_delta_e
    scen = [
        ("vs_v2_hold38000", {}, {"hold": 38000}, (SVS.REQ_SETTLED, 80000, 42000)),
        ("vs_v2_hold_default", {}, {"hold": None}, (SVS.REQ_SETTLED, 82000, 42000)),
        ("vs_v2_break_restart", {"jump_at": 60001, "jump": 0.01}, {"hold": 38000, "cap": 130000}, (SVS.REQ_SETTLED, 120000, 82000)),
        ("vs_v2_diverged", {"nan_at": 50000}, {"hold": 38000}, (SVS.REQ_DIVERGED, 52000, None)),
        ("vs_v2_unset", {"extra_col": True}, {"hold": 38000}, (SVS.REQ_UNSET, 42000, None)),
    ]
    try:
        for name, sc, kw, (kind, step, since) in scen:
            rd = V.make_run(name, sc, **kw)
            rc, rec, st = V.run_case(rd)
            res = SVS.analyse(rd, plot=False)
            req = res["v2_emulation"]["request"] or {}
            hs = res["v2_emulation"]["hold"]
            v2 = (st or {}).get("request_kind"), (st or {}).get("request_step"), (st or {}).get("settled_since")
            mine = req.get("kind"), req.get("step"), req.get("settled_since")
            chk(f"(vi-b) {name}: v2 の要求 {v2} = 再生の要求 {mine} = 期待 {(kind, step, since)}",
                v2 == mine == (kind, step, since), f"v2 status {(st or {}).get('status')}")
            v2h = (st or {}).get("hold") or {}
            chk(f"(vi-b) {name}: 保持の経過 (start・reached・breaks の start/broken_at/why) が v2 と同じ",
                v2h.get("reached") == hs["reached"] and v2h.get("start") == hs["start"] and v2h.get("hold_steps") == hs["hold_steps"]
                and v2h.get("breaks") == hs["breaks"],
                f"v2 {json.dumps({k: v2h.get(k) for k in ('start', 'reached')})} breaks {len(v2h.get('breaks') or [])}; "
                f"再生 {json.dumps({k: hs[k] for k in ('start', 'reached')})} breaks {len(hs['breaks'])}")
            if kind == SVS.REQ_SETTLED:
                chk(f"(vi-b) {name}: t_s = v2 の要求の格子点 {step}、v2 は止めたので span のデータが無く未成立",
                    res["t_s"] == step and res["verdict"] == "未成立" and any("span の終わり" in x for x in res["incomplete"]),
                    f"t_s {res['t_s']} {res['verdict']}")
            else:
                chk(f"(vi-b) {name}: 停止候補なし (v2 は {kind} で止める) → 未成立", res["t_s"] is None and res["verdict"] == "未成立"
                    and any(kind in r for r in res["reasons"]), str(res["reasons"][:1]))
            if name == "vs_v2_break_restart":
                chk("(vi-b) 保持の途中で崩れる: 42,000 から数えた保持が 62,000 で外れ (breaks)、82,000 から数え直す",
                    [(b["start"], b["broken_at"]) for b in hs["breaks"]] == [(42000, 62000)] and hs["reached"] == [42000, 82000],
                    str(hs["breaks"][:1]))
    finally:
        SR.extract_delta_e = real
    print(f"  (偽の forge 5 本 {time.perf_counter() - t0:.1f} s)")
    if not [f for f in FAILS if f.startswith("(vi-b)")]:
        shutil.rmtree(V.TMP, ignore_errors=True)
    else:
        print(f"  偽の forge の run を残す: {V.TMP}")


# ------------------------------------------------------------------ --bench
def bench():
    print("=== --bench: 20 万 step・格子 2,000 の規模の再生 ===")
    ni, nj, S = 750, 8, 0.1
    end = 200000
    run = TS.TMP / "bench"
    run.mkdir()
    i = np.arange(ni)[:, None]
    j = np.arange(nj)[None, :]
    xm = ((i - 2.0) * 0.02 * np.ones((1, nj))) * S
    rw = 1.0 + 0.002 * i
    rm = (rw * j / (nj - 1)) * S
    coord = np.stack([xm.ravel(), rm.ravel(), np.zeros(ni * nj)], axis=1)
    cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1"]
    dq = ["rms_dq_ro", "rms_dq_roUx", "rms_dq_roUy", "rms_dq_roUz", "rms_dq_roe"]
    t0 = time.perf_counter()
    # 値は 17 桁 (forge の残差の履歴と同じ書式; run_0005 の最初の行の値) — 読み込みの所要を実物の大きさ (約 120 MB) で測る
    vals = ("1.281468768342701e-05,0.01678349617095016,0.006994969463822327,0,14.3520001634708,0.08190914076396126,"
            "87258.16351208636,1.171463426901659e-05,1.099448829037861e-06,0,0,0,0,0")
    with open(run / "residual_history.csv", "w") as f:
        f.write("step,inner,phase," + ",".join(cols + dq) + "\n")
        for s in range(end):
            f.write(f"{s},-1,outer_begin,{vals}\n{s},0,inner_begin,{vals}\n{s},-1,outer_end,{vals}\n")
    a = 300.0
    eta = rm / rm[:, -1:]
    th = np.radians(0.5) * eta
    ux = (5.5 * a * np.cos(th)).ravel()
    uy = (5.5 * a * np.sin(th)).ravel()
    wc = np.stack([xm[:, -1], rm[:, -1], np.zeros(ni)], axis=1)
    conne = []
    for k in range(ni - 1):
        conne += [2, 2, k, k + 1]
    with open(run / "boundary_massflux.csv", "w") as f:
        f.write("step,mdot_inlet_1,raw_inlet_1,mdot_outlet_2,raw_outlet_2\n")
        for s in range(DT, end + 1, DT):
            f.write(f"{s},-10.0,-1.0,10.0,1.0\n")
            with h5py.File(run / f"res_nodes_{s}.h5", "w") as h:
                h.attrs["step"] = np.int32(s)
                h["/MESH/COORD"] = coord.ravel()
                h["/MESH/NODE_ID"] = np.arange(ni * nj, dtype=np.int64)
                for k, v in (("ro", np.ones(ni * nj)), ("Ux", ux), ("Uy", uy), ("Uz", np.zeros(ni * nj)), ("P", np.ones(ni * nj)),
                             ("T", np.ones(ni * nj)), ("sonic", np.full(ni * nj, a))):
                    h["/VALUE/" + k] = v
            with h5py.File(run / f"res_wall_3_{s}.h5", "w") as h:
                h.attrs["step"] = np.int32(s)
                h["/MESH/COORD"] = wc.ravel()
                h["/MESH/CONNE"] = np.asarray(conne, dtype=np.int32)
                h["/VALUE/iface_Qf_eff"] = 100.0 + np.arange(ni, dtype=float)
    crit = {"schema": settle.SCHEMA, "structured": {"ni": ni, "nj": nj, "scale_m": S},
            "grid": {"out_step_interval": DT, "start_step": 0}, "window": {"W": 20000, "min_points": 11}, "hold": 40000,
            "quantities": [
                {"kind": "residual", "required": [c for c in cols if c != "rms_roUz"], "inactive_ok": ["rms_roUz"],
                 "thresholds": {"slope_dec": math.log10(1.05), "amp_dec": 0.1, "rise_dec": math.log10(1.10), "excursion_dec": 0.5}},
                {"kind": "mass_balance", "inlet": [1], "outlet": [2], "thresholds": {"mean_abs": 5e-4, "swing": 5e-3, "mdot_change_rel": 1e-4}},
                {"kind": "exit_core_mach", "eta_range": [0.05, 0.7], "thresholds": {"tol": 0.000275}},
                {"kind": "line_mach", "eta": 0.15, "x_range_m": [0.5 * S, 12.0 * S], "n": 2001, "Md": 5.5,
                 "thresholds": {"tol_maxabs": 0.01, "tol_halfamp": 0.005, "tol_pos": 0.01}},
                {"kind": "exit_flow_angle", "eta_range": [0.05, 0.7], "thresholds": {"tol": 0.01}},
                {"kind": "wall_heat_nodes", "walls": [{"name": "wall", "physID": 3}], "axisymmetric": True,
                 "thresholds": {"q_tol": 1000.0, "eps_rel": 1e-3}}]}
    (run / settle.CRITERIA_FILE).write_text(yaml.safe_dump(crit, allow_unicode=True, sort_keys=False))
    sz = (run / "residual_history.csv").stat().st_size
    print(f"  合成の run を作った: {time.perf_counter() - t0:.1f} s (残差 {sz / 1e6:.0f} MB・{3 * end + 1} 行、節点 {ni * nj}、壁 {ni} 節点、"
          f"格子点 {end // DT})")
    t1 = time.perf_counter()
    res = SVS.analyse(run)
    tt = time.perf_counter() - t1
    print(f"  再生と判定: {tt:.1f} s ({', '.join(f'{k} {v:.1f} s' for k, v in res['timing_s'].items())}); VERDICT {res['verdict']}、"
          f"t_s {res['t_s']}、格子点 {res['grid']['n_points']}")
    return tt


def main():
    print(f"numpy {np.__version__}, python {sys.version.split()[0]} ({sys.executable}), 作業 {TS.TMP}")
    test_results()
    test_hold_unit()
    if not QUICK:
        test_v2_end_to_end()
    if BENCH:
        bench()
    shutil.rmtree(TS.TMP, ignore_errors=True)
    print(f"\n{NCHK - len(FAILS)}/{NCHK} OK")
    if FAILS:
        print("FAIL:", *FAILS, sep="\n  ")
    print(f"VERDICT: {'PASS' if not FAILS else 'FAIL'}")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
