"""case/70: 既存の NS の run に停止規則 (design/forge_design/evaluate/stop_rule.py) を当てる — plan tooling-nozzle-cooled-wall-ns-settling
§6 S-stop (c) (forge 0 step、run は読むだけ; 合否は付けない、閾値の妥当性を見る記述)。

対象と Euler 参照 (§6 S-stop (c)): run_0010b × run_0012b、run_0014b × run_0013、run_0014c × run_0013。y_b は各 run の最初の出力 (res_1000) で
`stop_rule.fix_band` (= `extract_and_merge(.., band_select="edge")` と同じ抽出) により固定し、各出力を `stop_rule.snapshot` → `first_stop` の順に
見張りと同じ規則で判定する。REACHED ならその step、届かなければ「未到達 (run の終わりまで)」(見張りの CENSORED に当たる)。

記録 (`_stop_rule_replay/replay.json`): run ごとに REACHED の step か未到達、その時点 (未到達なら最後の出力) の r の窓の平均・振幅、各断面の D_θ・
振幅、出口コア M と run の最後の出口コア M の差。**参考 (判定ではない)**: r の条件を除いた θ_r だけの条件で止まる step とその時点の出口コア M の差
(r が全出力で条件を外れるとき、θ_r が遅いモードを捕まえるかを見るため)。各出力の量と判定の内訳、y_b(x) も残す。
図 `_stop_rule_replay/replay_series.png`: θ_r (最後の出力に対する相対差)・r・出口コア M の系列 (横軸は通算の step: run_0010 / 0014 の 60k から続く)。

usage (case dir、design/.venv-opt):  python stop_rule_replay.py [RUN ...]   (既定は 3 run)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate import stop_rule as SR  # noqa: E402

OUT = C / "_stop_rule_replay"
TARGETS = {"run_0010b_ns_c2b_cont60k": "run_0012b_euler_c2b_cont12k",
           "run_0014b_ns_ns2_cont60k": "run_0013_euler_ns2",
           "run_0014c_ns_ns2_cont120k": "run_0013_euler_ns2"}
# 通算の step の起点 (README の run 一覧: run_0010b = run_0010 の 60k から、run_0014b = run_0014 の 60k から、run_0014c = 通算 120k から)
OFFSET = {"run_0010b_ns_c2b_cont60k": 60000, "run_0014b_ns_ns2_cont60k": 60000, "run_0014c_ns_ns2_cont120k": 120000}


def theta_only_stop(rows, th):
    """参考: r の条件を除き、θ_r の条件 (窓・点数・D_θ・振幅) だけを満たす最初の出力。"""
    for k in range(len(rows)):
        j = SR.judge_rows(rows, k, **th)
        if "theta" in j and all(t["ok"] for t in j["theta"].values()):
            return k
    return None


def at(rows, k, th, last_M):
    j = SR.judge_rows(rows, k, **th)
    out = {"step": rows[k]["step"], "n_window": j["n_window"], "r_mean": j.get("r_mean"), "r_amp": j.get("r_amp"),
           "theta": {key: {kk: t[kk] for kk in ("D", "D_raw", "amp", "mean")} for key, t in j.get("theta", {}).items()},
           "exit_core_M": rows[k]["exit_core_M"], "exit_core_M_last_minus_this": last_M - rows[k]["exit_core_M"],
           "ok": j["ok"], "why": j.get("why")}
    return out


def replay(run_name: str, euler_name: str) -> dict:
    run = C / run_name
    t0 = time.time()
    geo = SR.RunGeometry(run)
    steps = [s for s in SR.res_steps(run) if s > 0]
    od = OUT / run_name
    od.mkdir(parents=True, exist_ok=True)
    th = SR.thresholds()
    try:
        band = SR.fix_band(run, steps[0], C / euler_name, od / "band")
    except SR.Unavailable as e:
        return {"run": run_name, "euler_ref": euler_name, "status": "UNAVAILABLE", "reason": str(e)}
    print(f"[replay] {run_name}: y_b を res_{steps[0]} × {euler_name} で固定 ({band['seconds']} s)", flush=True)
    rows = [SR.snapshot(geo, s, band, th["x_fracs"]) for s in steps]
    for k, r in enumerate(rows):
        r["judge"] = SR.judge_rows(rows, k, **th)
    stop = SR.first_stop(rows, **th)
    last_M = rows[-1]["exit_core_M"]
    k_th = theta_only_stop(rows, th)
    rec = {"run": run_name, "euler_ref": euler_name, "offset_step": OFFSET.get(run_name, 0), "n_outputs": len(rows),
           "first_step": steps[0], "last_step": steps[-1], "thresholds": th,
           "status": (stop["status"] if stop else "未到達 (run の終わりまで; 見張りの CENSORED に当たる)"),
           "reach_step": (stop["step"] if stop and stop["status"] == "REACHED" else None),
           "at_stop_or_last": at(rows, stop["k"] if stop else len(rows) - 1, th, last_M),
           "r_all_outputs": {"min": float(min(r["r"] for r in rows)), "max": float(max(r["r"] for r in rows))},
           "theta_only_reference": ({"step": rows[k_th]["step"], **at(rows, k_th, th, last_M)} if k_th is not None else None),
           "exit_core_M_last": last_M,
           "band": {k: band[k] for k in ("step", "res_sha256", "euler_ref", "extract", "out_dir", "seconds")},
           "band_y_b": {"x": band["x"], "y_b": band["y_b"]},
           "rows": rows, "seconds": round(time.time() - t0, 1)}
    s = rec["at_stop_or_last"]
    print(f"[replay] {run_name}: {rec['status']} {rec['reach_step']}  r 平均 {s['r_mean']:+.4%} 振幅 {s['r_amp']:.4%}  "
          + "  ".join(f"{k} D {v['D']:.4%} 振幅 {v['amp']:.4%}" for k, v in s["theta"].items())
          + f"  出口 M {s['exit_core_M']:.6f} (最後との差 {s['exit_core_M_last_minus_this']:+.2e})"
          + (f"  [θ_r だけなら {rec['theta_only_reference']['step']}]" if rec["theta_only_reference"] else "  [θ_r だけでも届かない]"), flush=True)
    return rec


def figure(recs: list, path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    import os
    font = os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")
    if os.path.exists(font):
        fm.fontManager.addfont(font)
        matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    fig, ax = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    colors = {"run_0010b_ns_c2b_cont60k": "tab:blue", "run_0014b_ns_ns2_cont60k": "tab:orange", "run_0014c_ns_ns2_cont120k": "tab:red"}
    styles = {0.4: "-", 0.7: "--", 0.95: ":"}
    for rec in recs:
        if "rows" not in rec:
            continue
        name = rec["run"]; off = rec["offset_step"]; col = colors.get(name, "k")
        st = np.array([r["step"] for r in rec["rows"]]) + off
        for f in rec["thresholds"]["x_fracs"]:
            v = np.array([r[SR.theta_key(f)] for r in rec["rows"]])
            ax[0].plot(st, 100 * (v / v[-1] - 1), styles.get(f, "-"), color=col, lw=1.2,
                       label=f"{name.split('_ns')[0]} x/x_E {f:g}")
        ax[1].plot(st, 100 * np.array([r["r"] for r in rec["rows"]]), color=col, lw=1.2, label=name.split("_ns")[0])
        ax[2].plot(st, [r["exit_core_M"] for r in rec["rows"]], color=col, lw=1.2, label=name.split("_ns")[0])
        for a in ax:
            if rec["reach_step"] is not None:
                a.axvline(rec["reach_step"] + off, color=col, lw=1.0, ls="-.")
            if rec["theta_only_reference"]:
                a.axvline(rec["theta_only_reference"]["step"] + off, color=col, lw=0.8, ls=":")
    ax[0].set_ylabel("θ_r / θ_r(最後) − 1 [%]")
    ax[0].set_title("停止規則の replay (plan tooling-nozzle-cooled-wall-ns-settling §6 S-stop (c)); 縦の点線 = r を除いた θ_r だけの条件で止まる出力 (参考)")
    ax[1].axhspan(-100 * SR.R_MEAN_MAX, 100 * SR.R_MEAN_MAX, color="0.85", label="|r の平均| ≤ 0.05 %")
    ax[1].set_ylabel("r = (ṁ_in − ṁ_out)/ṁ_in [%]")
    ax[2].set_ylabel("出口コア M")
    ax[2].set_xlabel("通算の step (run_0010 / run_0014 の 0 から)")
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main(argv=None) -> int:
    names = (argv if argv else None) or list(TARGETS)
    OUT.mkdir(exist_ok=True)
    recs = [replay(n, TARGETS[n]) for n in names]
    summary = [{k: r.get(k) for k in ("run", "euler_ref", "status", "reach_step", "at_stop_or_last", "theta_only_reference",
                                       "r_all_outputs", "exit_core_M_last", "first_step", "last_step", "n_outputs", "seconds", "reason")}
               for r in recs]
    p = OUT / "replay.json"
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps({"plan": "plans/active/tooling-nozzle-cooled-wall-ns-settling.md §6 S-stop (c)", "summary": summary, "runs": recs},
                              indent=1, ensure_ascii=False, default=SR._json_default))
    tmp.replace(p)
    figure(recs, OUT / "replay_series.png")
    print(f"[replay] → {p}、{OUT / 'replay_series.png'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
