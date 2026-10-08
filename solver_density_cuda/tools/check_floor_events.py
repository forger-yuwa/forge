#!/usr/bin/env python3
"""毎更新の EOS 床事象の判定 (記録 `floor_events.csv` を読む; plans/active/tooling-sern-te-wake-grid.md §4「床の判定」・§6 #1)。

forge の `output: {floorEvents: 1}` が run ディレクトリに追記する `floor_events.csv` (書式は
`cuda_forge/floorEvents_d.cu` 冒頭) から、**判定区間 (a, b] の全更新・全実節点で床事象が 0 件か**を判定する。
保存標本 (res_*.h5) の「床近傍 0 節点」は過去の更新事象を復元できないので、その代わりにならない。

時点の契約 (定常陰解法 = path steady_implicit):
  - 外側 step k の EOS は直前の更新の結果 Q_{k-1} を判定する (行 kind=eos, step=k, inner=0, q_index=k-1)。
  - 最後の更新の結果 Q_N は終了時の監査 (行 kind=audit, q_index=N) が判定する。
  - 判定区間 (a, b] は b = N (最後の session の最終 step)、a = N − window_steps (既定は末尾 tail = 0.5)。
    入口 Q_a (行 q_index=a) は**別枠**で記録して合否には使わない。合否は Q_{a+1} … Q_N の行。
  - 初期化の EOS (kind=init)・時間ステップ外の EOS (kind=aux) は更新事象ではない (別枠で要約に載せる)。

判定:
  PASS          区間の全行がそろい、実節点の温度床・密度床・圧力床が全て 0 件
  FAIL          区間内に床事象 (または非有限の入力) がある
  INDETERMINATE 判定不能 = 合格ではない。記録が無い / 読めない行 / 見出し違い / session_end・監査が無い /
                区間の step が欠けた・重複 / 区間が最後の session (restart 後) の始まりより前に及ぶ /
                未検証の経路 (dual_time・explicit) / 温度の述語を評価していない節点 (二相など) /
                述語と EOS の最終温度の食い違い / 件数の整合が取れない (overflow) / 記録と res_*.h5 の最終 step が違う

使い方: check_floor_events.py RUN_DIR [--window-steps N | --tail 0.5] [--allow-path steady_implicit]
終了コード: PASS 0、それ以外 1。
"""
import argparse
import csv
import glob
import os
import re
import sys

FILE_NAME = "floor_events.csv"
COLUMNS = ("kind,session,step,inner,q_index,path,n_real,n_ghost,n_nonfinite_real,n_nonfinite_ghost,"
           "nT_real,nRho_real,nP_real,nT_ghost,nRho_ghost,nP_ghost,nT_uneval_real,nT_mismatch_real,n_near_real,"
           "dRhoE_T_sum,dRhoE_T_max,dRho_sum,dRho_max,dP_sum,dP_max,"
           "ids_T,ids_Rho,ids_P,ids_near,ids_truncated,overflow,note").split(",")
INT_COLS = ("session", "step", "inner", "q_index", "n_real", "n_ghost", "n_nonfinite_real", "n_nonfinite_ghost",
            "nT_real", "nRho_real", "nP_real", "nT_ghost", "nRho_ghost", "nP_ghost", "nT_uneval_real",
            "nT_mismatch_real", "n_near_real", "overflow")
FLOAT_COLS = ("dRhoE_T_sum", "dRhoE_T_max", "dRho_sum", "dRho_max", "dP_sum", "dP_max")
KINDS = ("session_begin", "init", "eos", "aux", "audit", "session_end")
EVENT_KINDS = (("T", "nT_real", "dRhoE_T_sum", "dRhoE_T_max", "ids_T"),
               ("Rho", "nRho_real", "dRho_sum", "dRho_max", "ids_Rho"),
               ("P", "nP_real", "dP_sum", "dP_max", "ids_P"))
DEFAULT_ALLOWED_PATHS = ("steady_implicit",)   # 時点の契約を検証した経路だけ


def _note_dict(s):
    out = {}
    for kv in (s or "").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def last_res_step(run_dir):
    """res_<N>.h5 (体積出力) の最大 N。無ければ None。"""
    steps = []
    for f in glob.glob(os.path.join(run_dir, "res_*.h5")):
        m = re.fullmatch(r"res_(\d+)\.h5", os.path.basename(f))
        if m:
            steps.append(int(m.group(1)))
    return max(steps) if steps else None


def load(path):
    """記録を読む。戻り (rows, errors)。errors が空でなければ判定不能。"""
    errors, rows = [], []
    if not os.path.exists(path):
        return [], [f"記録が無い ({FILE_NAME})"]
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        try:
            header = next(rd)
        except StopIteration:
            return [], ["記録が空"]
        if header != COLUMNS:
            return [], ["見出しが現行の書式と違う"]
        for ln, rec in enumerate(rd, start=2):
            if len(rec) != len(COLUMNS):
                errors.append(f"{ln} 行目: 列数 {len(rec)} (期待 {len(COLUMNS)}) — 書きかけ・壊れた行")
                continue
            r = dict(zip(COLUMNS, rec))
            if r["kind"] not in KINDS:
                errors.append(f"{ln} 行目: 未知の kind {r['kind']!r}")
                continue
            try:
                for c in INT_COLS:
                    r[c] = int(r[c])
                for c in FLOAT_COLS:
                    r[c] = float(r[c])
            except ValueError as e:
                errors.append(f"{ln} 行目: 数値として読めない ({e})")
                continue
            r["_line"] = ln
            rows.append(r)
    return rows, errors


def judge(run_dir, window_steps=None, tail=0.5, allowed_paths=DEFAULT_ALLOWED_PATHS, check_res_step=True):
    """戻り dict: verdict (PASS/FAIL/INDETERMINATE), ok, reasons, window (a, b], 区間内の件数・最初の事象、
    入口・初期化・aux の別枠、床近傍 (補助)、ghost (別集計)、session の情報。"""
    path = os.path.join(run_dir, FILE_NAME)
    out = {"verdict": "INDETERMINATE", "ok": False, "reasons": [], "file": path, "window": None,
           "n_updates_checked": 0, "events": {}, "first_events": [], "entry": None, "init": None, "aux": None,
           "near_window": 0, "ghost_window": {}, "session": None, "path": None, "last_step": None}
    rows, errors = load(path)
    if errors:
        out["reasons"] = ["判定不能: " + e for e in errors[:5]] + ([f"… ほか {len(errors) - 5} 件"] if len(errors) > 5 else [])
        return out
    begins = [r for r in rows if r["kind"] == "session_begin"]
    if not begins:
        out["reasons"] = ["判定不能: session_begin が無い"]
        return out
    sess = begins[-1]["session"]
    if any(r["session"] > sess for r in rows) or [r["session"] for r in begins] != sorted(set(r["session"] for r in begins)):
        out["reasons"] = ["判定不能: session の番号が順に並んでいない"]
        return out
    first_idx = rows.index(begins[-1])
    srows = rows[first_idx:]
    if any(r["session"] != sess for r in srows):
        out["reasons"] = ["判定不能: 最後の session の後に別の session の行がある"]
        return out
    note = _note_dict(begins[-1]["note"])
    out["session"] = {"session": sess, "prior_sessions": int(note.get("prior_sessions", "-1") or -1),
                      "note": note}
    out["path"] = begins[-1]["path"]
    reasons = []
    ends = [r for r in srows if r["kind"] == "session_end"]
    if len(ends) != 1 or srows[-1]["kind"] != "session_end":
        out["reasons"] = ["判定不能: 最後の session に session_end が無い (異常終了または記録の欠け)"]
        return out
    N = ends[0]["step"]
    out["last_step"] = N
    if out["path"] not in allowed_paths:
        out["reasons"] = [f"判定不能: 未検証の経路 {out['path']} (時点の契約を検証した経路: {', '.join(allowed_paths)})"]
        return out
    if _note_dict(ends[0]["note"]).get("audit") != "done":
        out["reasons"] = [f"判定不能: 終了時の監査が無い ({ends[0]['note'] or 'audit 不明'})"]
        return out
    n_real = begins[-1]["n_real"]
    if any(r["n_real"] != n_real for r in srows):
        out["reasons"] = ["判定不能: session の途中で実節点数が変わった"]
        return out
    if check_res_step:
        rs = last_res_step(run_dir)
        if rs is not None and rs != N:
            out["reasons"] = [f"判定不能: 記録の最終 step {N} と res_*.h5 の最終 step {rs} が違う (古い記録・別の起動の記録)"]
            return out
    if window_steps is not None:
        w = int(window_steps)
        if w <= 0:
            out["reasons"] = [f"判定不能: window_steps {w} ≤ 0"]
            return out
        if w > N:
            out["reasons"] = [f"判定不能: 判定区間 {w} step が最後の session (restart 後) の長さ {N} step を超える — "
                              "区間の前半の記録が無い"]
            return out
        a = N - w
    else:
        if not (0.0 < float(tail) <= 1.0):
            out["reasons"] = [f"判定不能: tail {tail} は (0, 1]"]
            return out
        a = N - int(round(float(tail) * N))
    if N <= a:
        out["reasons"] = [f"判定不能: 判定区間 ({a}, {N}] が空"]
        return out
    out["window"] = [a, N]

    eos = [r for r in srows if r["kind"] == "eos"]
    audits = [r for r in srows if r["kind"] == "audit"]
    # 定常陰解法では 1 step に EOS は 1 回 (inner 0)。それ以外の形は時点の契約から外れる
    extra = [r for r in eos if r["inner"] != 0 or r["q_index"] != r["step"] - 1]
    if extra:
        reasons.append(f"判定不能: 1 step に複数の EOS / q_index の不整合 ({len(extra)} 行、最初の step {extra[0]['step']})")
    if len(audits) != 1 or audits[0]["q_index"] != N or audits[0]["step"] != N:
        reasons.append(f"判定不能: 監査行 (q_index = {N}) が 1 行でない ({len(audits)} 行)")
    by_q = {}
    for r in eos:
        by_q.setdefault(r["q_index"], []).append(r)
    for r in audits:
        by_q.setdefault(r["q_index"], []).append(r)
    need = list(range(a, N + 1))   # a = 入口 (別枠)、a+1..N = 判定区間の更新結果
    missing = [q for q in need if q not in by_q]
    dup = [q for q in need if len(by_q.get(q, [])) > 1]
    if missing:
        reasons.append(f"判定不能: 区間の記録が欠けている (q_index {_ranges(missing)}) — 欠けた step を 0 件として扱わない")
    if dup:
        reasons.append(f"判定不能: 同じ更新結果の行が重複 (q_index {_ranges(dup)})")
    if reasons:
        out["reasons"] = reasons
        return out

    win = [by_q[q][0] for q in range(a + 1, N + 1)]
    entry = by_q[a][0]
    out["entry"] = _row_summary(entry)
    out["n_updates_checked"] = len(win)
    inits = [r for r in srows if r["kind"] == "init"]
    auxs = [r for r in srows if r["kind"] == "aux"]
    out["init"] = _sum_rows(inits)
    out["aux"] = _sum_rows(auxs)
    for r in win:
        if r["overflow"] != 0:
            reasons.append(f"判定不能: 件数の整合が取れない行 (q_index {r['q_index']}、overflow)")
        if r["nT_uneval_real"] > 0:
            reasons.append(f"判定不能: 温度の述語を評価していない実節点 {r['nT_uneval_real']} (q_index {r['q_index']}; 二相など)")
        if r["nT_mismatch_real"] > 0:
            reasons.append(f"判定不能: 述語と EOS の最終温度が食い違う実節点 {r['nT_mismatch_real']} (q_index {r['q_index']})")
    if reasons:
        out["reasons"] = reasons[:10] + ([f"… ほか {len(reasons) - 10} 件"] if len(reasons) > 10 else [])
        return out

    ev = {}
    for name, cn, cs, cm, ci in EVENT_KINDS:
        hit = [r for r in win if r[cn] > 0]
        ev[name] = {"n_events": int(sum(r[cn] for r in win)), "n_updates_with_events": len(hit),
                    "sum": float(sum(r[cs] for r in win)), "max": float(max([r[cm] for r in win] or [0.0])),
                    "first_q_index": hit[0]["q_index"] if hit else None,
                    "first_ids": hit[0][ci] if hit else ""}
    out["events"] = ev
    nonfin = [r for r in win if r["n_nonfinite_real"] > 0]
    out["near_window"] = int(sum(r["n_near_real"] for r in win))
    out["ghost_window"] = {k: int(sum(r[c] for r in win)) for k, c in (("T", "nT_ghost"), ("Rho", "nRho_ghost"), ("P", "nP_ghost"),
                                                                        ("nonfinite", "n_nonfinite_ghost"))}
    firsts = []
    for r in win:
        n = r["nT_real"] + r["nRho_real"] + r["nP_real"] + r["n_nonfinite_real"]
        if n > 0:
            firsts.append({"q_index": r["q_index"], "kind": r["kind"], "nT": r["nT_real"], "nRho": r["nRho_real"],
                           "nP": r["nP_real"], "n_nonfinite": r["n_nonfinite_real"],
                           "ids_T": r["ids_T"], "ids_Rho": r["ids_Rho"], "ids_P": r["ids_P"]})
            if len(firsts) >= 5:
                break
    out["first_events"] = firsts
    total = sum(v["n_events"] for v in ev.values())
    if total > 0 or nonfin:
        msg = ", ".join(f"{k} {v['n_events']} 件 ({v['n_updates_with_events']} 更新, 最初 q_index {v['first_q_index']})"
                        for k, v in ev.items() if v["n_events"] > 0)
        if nonfin:
            msg += (", " if msg else "") + f"非有限の入力 {sum(r['n_nonfinite_real'] for r in nonfin)} 件"
        out["verdict"] = "FAIL"
        out["reasons"] = [f"判定区間 ({a}, {N}] に床事象: {msg}"]
        return out
    out["verdict"] = "PASS"; out["ok"] = True
    return out


def _ranges(qs):
    qs = sorted(qs)
    if not qs:
        return ""
    parts, s, p = [], qs[0], qs[0]
    for q in qs[1:]:
        if q == p + 1:
            p = q; continue
        parts.append(f"{s}" if s == p else f"{s}-{p}"); s = p = q
    parts.append(f"{s}" if s == p else f"{s}-{p}")
    txt = ",".join(parts)
    return txt if len(txt) < 120 else txt[:117] + "..."


def _row_summary(r):
    return {"q_index": r["q_index"], "step": r["step"], "nT": r["nT_real"], "nRho": r["nRho_real"], "nP": r["nP_real"],
            "n_nonfinite": r["n_nonfinite_real"], "n_near": r["n_near_real"]}


def _sum_rows(rs):
    return {"n_rows": len(rs), "nT": sum(r["nT_real"] for r in rs), "nRho": sum(r["nRho_real"] for r in rs),
            "nP": sum(r["nP_real"] for r in rs), "n_nonfinite": sum(r["n_nonfinite_real"] for r in rs),
            "nT_ghost": sum(r["nT_ghost"] for r in rs), "nRho_ghost": sum(r["nRho_ghost"] for r in rs),
            "nP_ghost": sum(r["nP_ghost"] for r in rs)}


def main(argv=None):
    ap = argparse.ArgumentParser(description="毎更新の EOS 床事象の判定 (floor_events.csv)")
    ap.add_argument("run_dir")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--window-steps", type=int, default=None, help="判定区間の長さ [step] (末尾から)")
    g.add_argument("--tail", type=float, default=0.5, help="判定区間 = 最後の session の末尾の割合 (既定 0.5)")
    ap.add_argument("--allow-path", action="append", default=None,
                    help=f"受け入れる経路 (既定 {', '.join(DEFAULT_ALLOWED_PATHS)}。未検証の経路を足すときは理由を記録すること)")
    a = ap.parse_args(argv)
    res = judge(a.run_dir, window_steps=a.window_steps, tail=a.tail,
                allowed_paths=tuple(a.allow_path) if a.allow_path else DEFAULT_ALLOWED_PATHS)
    print(f"floor events: {res['file']}")
    if res["session"]:
        print(f"  session {res['session']['session']} (prior {res['session']['prior_sessions']}), path {res['path']}, "
              f"last step {res['last_step']}, window {res['window']}")
    if res["window"] and res["events"]:
        for k, v in res["events"].items():
            print(f"  {k:3s}: {v['n_events']} 件 / {v['n_updates_with_events']} 更新 (Σ {v['sum']:.6g}, max {v['max']:.6g})"
                  + (f", 最初 q_index {v['first_q_index']} ids {v['first_ids']}" if v["n_events"] else ""))
        print(f"  入口 Q_a (別枠): {res['entry']}")
        print(f"  床近傍 T ≤ T_min + 1 K (補助、区間内の延べ): {res['near_window']}、ghost (別集計): {res['ghost_window']}")
        print(f"  初期化の EOS (別枠): {res['init']}、時間ステップ外の EOS (別枠): {res['aux']}")
    for r in res["reasons"]:
        print(f"  - {r}")
    print(f"VERDICT: {res['verdict']}")
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
