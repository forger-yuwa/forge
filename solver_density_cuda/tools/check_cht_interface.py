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

**帯を限った判定 (`--band-y YTOP YBOT`)** (plan §5.1 #101 ②): 全域 max しか持たない履歴では
「帯内で n_consec 更新連続」を判定できないので、ソルバの `conjugate.node_log: 1` が書く
`conjugate_iface_log_<physID>.csv` (毎更新・界面全節点の $r_i, Q_{f,i}, \Delta T_i$) と
`conjugate_iface_nodes_<physID>.csv` ($y_i, A_i$) から、**帯内節点だけ**で ① ② ③ を組む:

    ① max_{i∈帯} |r_i|/A_i            ② max_{i∈帯} |r_i| / max_{i∈帯} |Q_{f,i}|
    ③ max_{i∈帯} |ΔT_i|

④ 固体内部残差は**全域のまま** (`conjugate_history.csv` の `res_solid`)。帯は結果を見て動かさないこと
(評価器が別に決めて渡す)。結果は `CHT_INTERFACE_BAND_VERDICT.txt` に書く (全域の VERDICT は上書きしない)。
更新番号が末尾窓で連続していない・帯に節点が無い・履歴と更新数が食い違う、は `REFUSED`。

使い方:
  python3 solver_density_cuda/tools/check_cht_interface.py <run_dir> \
      --eps-rel 1e-3 --eps-abs 150 --dt-k 1e-2 --n-consec 80
  python3 solver_density_cuda/tools/check_cht_interface.py <run_dir> --band-y -0.004166 -0.019779
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

REQUIRED = ["step", "physID", "dTw_max", "res_abs_Wm2", "res_rel"]


VERDICT_FILE = "CHT_INTERFACE_VERDICT.txt"


def fail(run: Path, verdict: str, lines: list[str]) -> int:
    text = "\n".join([f"=== G-if: {run.name} -> {verdict} ==="] + lines + [f"\nVERDICT: {verdict}"])
    print(text)
    (run / VERDICT_FILE).write_text(text + "\n")
    return 0 if verdict == "PASS" else 1


def load_node_log(run: Path, pid: int):
    """節点表と節点ログを読み、**整合を全部検査して**から返す (codex result 8 巡目 M2)。

    検査しないと、座標 NaN の節点が帯判定から黙って落ちる・節点 ID の欠け/重複が配列代入で隠れる・
    更新内で 1 節点だけ別時刻の行が混ざる、がすべて PASS になった (メモリ上の改変で再現済み)。
    戻り値 (nodes[ni,5], log[:,7], updates[:], ni)。不整合は ValueError (理由つき)。
    `eval_v6p.py` の準定常系列もこれを使う。"""
    import numpy as np
    fn_nodes = run / f"conjugate_iface_nodes_{pid}.csv"
    fn_log = run / f"conjugate_iface_log_{pid}.csv"
    if not fn_nodes.exists() or not fn_log.exists():
        raise ValueError(f"{fn_log.name} / {fn_nodes.name} が無い (conjugate.node_log: 1 で回した run を渡すこと)")
    nd = np.loadtxt(fn_nodes, delimiter=",", skiprows=1, ndmin=2)
    if nd.shape[1] < 5 or not np.isfinite(nd).all():
        raise ValueError("節点表に欠損・非有限値がある (座標または面積)")
    ni = nd.shape[0]
    if not np.array_equal(nd[:, 0], np.arange(ni)):
        raise ValueError("節点表の ID が 0..n-1 の連番でない (欠け・重複・非整数)")
    if not (nd[:, 4] > 0).all():
        raise ValueError("節点表に面積 A_i <= 0 の節点がある")
    lg = np.loadtxt(fn_log, delimiter=",", skiprows=1, ndmin=2)
    if lg.shape[1] < 7 or not np.isfinite(lg).all():
        raise ValueError("節点ログに欠損・非有限値がある (NaN/Inf)")
    if not (np.all(lg[:, 0] == np.round(lg[:, 0])) and np.all(lg[:, 2] == np.round(lg[:, 2]))):
        raise ValueError("節点ログの更新番号・節点 ID が整数でない")
    upd = lg[:, 0].astype(np.int64)
    if np.any(np.diff(upd) < 0):
        raise ValueError("節点ログの更新番号が単調でない (再開で番号が戻った / 追記の混線)")
    ups, starts, counts = np.unique(upd, return_index=True, return_counts=True)
    bad = np.where(counts != ni)[0]
    if len(bad):
        raise ValueError(f"更新 {ups[bad[0]]} の行数 {counts[bad[0]]} が界面節点数 {ni} と違う")
    ids = lg[:, 2].astype(np.int64).reshape(len(ups), ni)
    if not (np.sort(ids, axis=1) == np.arange(ni)).all():
        raise ValueError("更新内の節点 ID が 0..n-1 をちょうど 1 回ずつ含まない (欠け・重複)")
    st = lg[:, 1].reshape(len(ups), ni)
    if not (st == st[:, :1]).all():
        k = int(np.where(~(st == st[:, :1]).all(axis=1))[0][0])
        raise ValueError(f"更新 {ups[k]} の中で step が節点ごとに違う (時刻の不一致)")
    return nd, lg, upd, ni


def band_values(run: Path, pid: int, ytop: float, ybot: float, n_consec: int, n_hist: int):
    """帯内節点で ① ② ③ を更新ごとに組む。戻り値 (vals dict, lines) か、(None, 理由 lines)。"""
    import numpy as np
    try:
        nd, lg, upd, ni = load_node_log(run, pid)
    except ValueError as e:
        return None, [f"  {e}"]
    y, A = nd[:, 3], nd[:, 4]
    lo, hi = min(ytop, ybot), max(ytop, ybot)
    inb = (y >= lo - 1e-9) & (y <= hi + 1e-9)       # 帯の端は 1e-9 m に丸めて書かれる
    if not inb.any():
        return None, [f"  帯 y∈[{lo:.6g}, {hi:.6g}] に界面節点が無い"]
    # **帯の節点集合を評価器の帯と照合する**: 評価器が同じ帯を v6p_band.json に書いていれば行数を比べる
    bj = run / "v6p_band.json"
    if bj.exists():
        b = json.loads(bj.read_text())
        if abs(b["ytop"] - hi) < 1e-9 and abs(b["ybot"] - lo) < 1e-9 and int(b["rows"]) != int(inb.sum()):
            return None, [f"  帯内の界面節点 {int(inb.sum())} が評価器の帯の行数 {b['rows']} と違う"]
    ups = np.unique(upd)
    if ups[-1] != n_hist:
        return None, [f"  節点ログの最終更新 {ups[-1]} と履歴の更新数 {n_hist} が食い違う"
                      " (再開で番号が戻った / 書き込み途中)"]
    tail_u = ups[-n_consec:]
    if len(tail_u) < n_consec or np.any(np.diff(tail_u) != 1):
        return None, [f"  節点ログの末尾 {n_consec} 更新が連続していない (欠け・重複)"]
    out = {"res_abs_Wm2": [], "res_rel": [], "dTw_max": [], "step": []}
    for u in tail_u:
        blk = lg[upd == u]
        idx = blk[:, 2].astype(int)
        r = np.empty(ni); Q = np.empty(ni); dT = np.empty(ni)
        r[idx], Q[idx], dT[idx] = blk[:, 3], blk[:, 4], blk[:, 5]
        rb, Qb, Ab, dTb = r[inb], Q[inb], A[inb], dT[inb]
        out["res_abs_Wm2"].append(float(np.max(np.abs(rb) / Ab)))
        qm = float(np.max(np.abs(Qb)))
        out["res_rel"].append(float(np.max(np.abs(rb))) / qm if qm > 0 else math.inf)
        out["dTw_max"].append(float(np.max(np.abs(dTb))))
        out["step"].append(int(blk[0, 1]))
    lines = [f"  帯 y∈[{lo*1e3:.3f}, {hi*1e3:.3f}] mm: 界面 {ni} 節点のうち {int(inb.sum())} 節点",
             f"  節点ログの末尾 {n_consec} 更新 (更新 {tail_u[0]} .. {tail_u[-1]}, "
             f"step {out['step'][0]} .. {out['step'][-1]})"]
    return out, lines


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
    ap.add_argument("--band-y", type=float, nargs=2, metavar=("YTOP", "YBOT"), default=None,
                    help="帯を限った判定 [m] (節点ログを使う。④ は全域のまま)")
    a = ap.parse_args()
    if a.band_y is not None:
        global VERDICT_FILE
        VERDICT_FILE = "CHT_INTERFACE_BAND_VERDICT.txt"

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
    band_lines = []
    if a.band_y is not None:
        if "update" not in rows[0]:
            return fail(run, "REFUSED", ["  履歴に update 列が無い (節点ログと突き合わせられない)"])
        bv, band_lines = band_values(run, pid, a.band_y[0], a.band_y[1], n_consec,
                                     int(float(rows[-1]["update"])))
        if bv is None:
            return fail(run, "REFUSED", band_lines)
        if [int(float(r["step"])) for r in tail] != bv["step"]:
            return fail(run, "REFUSED", ["  履歴の末尾窓と節点ログの末尾窓で step が一致しない"])
        for k in ("dTw_max", "res_abs_Wm2", "res_rel"):
            vals[k] = bv[k]                        # ① ② ③ は帯内。④ res_solid は全域のまま
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
             f"(step {tail[0]['step']} .. {tail[-1]['step']})"] + band_lines
    if a.band_y is not None:
        lines.append("  ① ② ③ は帯内節点、④ res_solid は全域")
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
