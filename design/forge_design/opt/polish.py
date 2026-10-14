"""前線ポリッシュ: チャンク継続 + η ドリフトゲート (子 plan §5 ステップ 5b)。

残差プラトー (STALLED) は既知のリミットサイクルで「未収束」とは限らないが、
報告量そのものの頭打ちで確認するのが規律 (準定常確認の η 直接版)。各対象 run の
最終 res から同一メッシュ継続 (+nStepOuter) を回し、|Δη| < tol (既定 3e-4 ≈
0.03%) に頭打ちするまで最大 max_chunks 回続ける。継続 run は `<tag>_pK`。

CLI: design/.venv-opt/bin/python -m forge_design.opt.polish <campaign_dir>
     (campaign の pareto.json の全点をポリッシュし polish.jsonl /
      pareto_polished.json を書く)

TP 経路 (plan tooling-nozzle-bell-two-mode §4.1・§5.1 #2b) は `polish_run_tp`: `runner.prepare` 経由の再補間 (CPG 経路の
現行実装。名前に反して同一メッシュの継続ではない) をやめ、restart_field.py による同一メッシュの継続とゲートの再判定にする。
CPG 経路も 2026-10-15 から同一形状の継続 (`runner.continue_cpg`: 前のチャンクの変換済みの格子を写し、段終了ゲート +
restart_field.py。plan tooling-generator-defaults-and-launch-gate §5.1 #19)。η の頭打ちの判定は従来どおり。
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ..evaluate import runner
from ..probdef import load_problem
from .ehvi import nondominated_mask


def polish_run_tp(prob_yaml, run_dir, eta0: float, tol: float = 3e-4, max_chunks: int = 2) -> dict:
    """TP 経路の polish (plan tooling-nozzle-bell-two-mode §4.1、§5.1 #2b): **同一メッシュの継続** (`runner.continue_tp` =
    restart_field.py で親の最終 res を写し、本段をもう 1 段; 段の履歴と manifest を引き継ぐので収束判定は親の本段と連結した
    区間)。チャンクごとに受理ゲート (i)〜(v) (`driver.tp_gates`) を判定し直し、|Δη| < tol **かつ**ゲート PASS で頭打ちとする
    (η の前後差だけでは終了しない)。継続 run は `<tag>_pK`。"""
    from .driver import tp_gates
    run_dir = Path(run_dir)
    prev = run_dir
    eta_prev = float(eta0)
    hist, gates = [eta_prev], []
    settled = False
    for k in range(1, max_chunks + 1):
        nd = run_dir.parent / f"{run_dir.name}_p{k}"
        if not nd.exists():  # 再開可能
            runner.continue_tp(prob_yaml, prev, nd)
        out = runner.collect(prob_yaml, nd)
        g = tp_gates(nd, out)
        eta = float(out["eta_cf"])
        hist.append(eta)
        gates.append({"run": nd.name, "status": g["status"], "fail_class": g["fail_class"], "note": g["note"]})
        drift, eta_prev, prev = abs(eta - eta_prev), eta, nd
        if drift < tol and g["status"] == "PASS":
            settled = True
            break
    return {"eta_final": eta_prev, "eta_hist": hist, "settled": settled, "chunks": len(hist) - 1,
            "gates": gates, "path": "tp"}


def polish_run(prob_yaml, run_dir, eta0: float, tol: float = 3e-4,
               max_chunks: int = 2) -> dict:
    """run_dir の解を η が頭打ちするまでチャンク継続する。TP 経路は `polish_run_tp` (同一メッシュの継続 + ゲート)。"""
    if runner._is_tp(load_problem(prob_yaml)):
        return polish_run_tp(prob_yaml, run_dir, eta0, tol=tol, max_chunks=max_chunks)
    run_dir = Path(run_dir)
    prev = run_dir
    eta_prev = float(eta0)
    hist = [eta_prev]
    settled = False
    for k in range(1, max_chunks + 1):
        nd = run_dir.parent / f"{run_dir.name}_p{k}"
        if not nd.exists():  # 再開可能
            # 同一形状の継続: 前のチャンクの変換済みの格子を写し、最終場を段終了ゲート + restart_field.py で初期場に
            # (`runner.continue_cpg`; 2026-10-15 まで prepare で格子を作り直して interp_field.py。plan
            # tooling-generator-defaults-and-launch-gate §5.1 #19)
            runner.continue_cpg(prob_yaml, prev, nd)
            runner.run_forge(nd)
        out = runner.collect(prob_yaml, nd)
        eta = float(out["eta_cf"])
        hist.append(eta)
        drift, eta_prev, prev = abs(eta - eta_prev), eta, nd
        if drift < tol:
            settled = True
            break
    return {"eta_final": eta_prev, "eta_hist": hist, "settled": settled,
            "chunks": len(hist) - 1}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="パレート前線の η 頭打ちポリッシュ")
    ap.add_argument("campaign_dir")
    ap.add_argument("--tol", type=float, default=3e-4)
    ap.add_argument("--max-chunks", type=int, default=2)
    a = ap.parse_args(argv)
    d = Path(a.campaign_dir)
    rows = [json.loads(l) for l in (d / "ledger.jsonl").read_text().splitlines() if l]
    front = json.load(open(d / "pareto.json"))["pareto"]
    out_rows = []
    for q in front:
        row = next(r for r in rows if r["status"] == "PASS"
                   and np.allclose(r["x"], q["x"], atol=1e-12))
        res = polish_run(d / f"{row['tag']}.yaml", row["run_dir"], row["eta_cf"],
                         tol=a.tol, max_chunks=a.max_chunks)
        rec = {"tag": row["tag"], "x": row["x"], "L_over_rt": q["L_over_rt"],
               "eta_raw": row["eta_cf"], **res}
        out_rows.append(rec)
        print(f"[{row['tag']}] η {row['eta_cf']:.5f} -> {res['eta_final']:.5f} "
              f"(chunks={res['chunks']}, settled={res['settled']})", flush=True)
    (d / "polish.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in out_rows) + "\n")
    F = np.array([[-r["eta_final"], r["L_over_rt"]] for r in out_rows])
    mask = nondominated_mask(F)
    polished = [{"x": out_rows[i]["x"], "eta_cf": out_rows[i]["eta_final"],
                 "L_over_rt": out_rows[i]["L_over_rt"],
                 "settled": out_rows[i]["settled"]}
                for i in np.where(mask)[0]]
    polished.sort(key=lambda r: r["L_over_rt"])
    (d / "pareto_polished.json").write_text(
        json.dumps({"tol": a.tol, "pareto": polished}, indent=1, ensure_ascii=False))
    n_set = sum(r["settled"] for r in out_rows)
    print(f"settled {n_set}/{len(out_rows)}, polished front {len(polished)} 点", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
