"""plan tooling-nozzle-isothermal-wall-chain §5.1 #31・§6「V-g0c-chain の事前登録」(commit 68d56723) の判定。B の case dir で vgc_run.py が最後に呼ぶ。

問い: G0c・FP64・B0 の停止の時点の出力の互換性 (float plan §6.58・§6.61) が、`run_0179/res_20000 → 格子をまたいで移す → S1 → S2 → 本段` という
チェーンの開始の手順へ広がるか。試験 run_0810_vgc_chain、参照 T2a (run_0752_t2_d1)・T2b (run_0754_t2_d2)、対照の出力 run_0804_sp_g0c (記録だけ)。

  1. 前提の検査の記録 (VGC_PRECHECK.json、vgc_run.py の (1)〜(6)) の読み直し。無い・不成立は判定不能。本段の直前の照合 (VGC_MAIN_START.json)、
     DISK_STOP・WATCH_ABORT (基盤の停止) も判定不能。RUN_PROVENANCE.txt のバイナリ・環境変数を登録と照合する。
  2. 段のゲート (記録、収束の合格は求めない): S1・S2 の正常終了 (rc 0、vgc_stages.jsonl)、最終場の保存量の有限性・ρ > 0・P > 0・T > 0
     (vgc_stage_gate.jsonl、段終了ゲートの問題 0 件)、全段の残差の履歴に非有限なし、stage_manifest の最後の区間が本段だけ、
     check_convergence --segment の VERDICT (記録)。
  3. 本段の見張り (m9_watch.json): REACHED は保存した系列から到達を計算し直し (m9_watch.py の done_at・spacing_ok をそのまま使う)、
     到達の出力の sha256 を照合する。CENSORED は系列に到達が無く上限まで回ったことを確かめる。
  4. M1〜M7 (float plan §6.58 と同じ式・許容): sp_judge.py の fields・compare をそのまま使う (ast で関数の定義だけを取り出す。
     取り出した文の sha256 をこの commit の値と照合)。候補 (到達の出力) を T2a・T2b のそれぞれと比べる。run_0804 との比較は記録だけ。
  5. 準定常 (記録): 本段の末尾 2 万 step の 9 出力で、§6.11 の閾値 (θ_r は --drift 0.0001 --osc 0.0001、Q_w は 0.0002) と、
     sp_judge.py と同じ呼び方 (4 量とも 0.0001) の両方。
  6. 時間 (記録): 総時間 = START_EPOCH (準備の前) から残差の履歴が到達の step に達した時刻まで、本段の時間 = MAIN_STAGE_STARTED から同じ時刻まで
     (timing.tsv の 30 秒ごとの記録を本段の行だけで線形補間)。1 step の時間・占有 (forge の数 ≤ 1)。起点の違う対照との比を速さの率にしない。

判定 (登録どおり):
  結果 A = 本段が停止規則に到達し、停止の出力が T2a・T2b のそれぞれと M1〜M7 で許容内
           →「登録した開始の手順で停止規則に到達し、T2 の両方の参照との比較が許容内」(この開始の手順への適用の拡大を支持)。
  結果 B = 本段 20 万 step までに未到達 (CENSORED)・発散 (DIVERGED、または前の段 S1・S2 の発散で本段に入れない)・どれかの許容の超過
           → 拡大を退け、いまの採用の範囲 (登録の Q32 の起点) を保つ。
  判定不能 = 入力の不一致・出力の欠損・見張りの誤作動・基盤の停止 (数値の手法の失敗に数えない)。
  解釈 (台本で決めた細部): 段のゲート (2) は記録だが、A の条件がそろっても段のゲートに不成立があれば A と書かず判定不能にする
  (B の理由があれば B)。前の段の失敗は、非有限・ρ ≤ 0・detectNaN の証拠があれば発散 (B)、無ければ判定不能。
→ _band_ab/cold_pair/vgc_judge.json
言わないこと: 定常解の一致・物理の精度・壁温の感度の保証、壁温を変えた問題・δ* ループの全 pass・別形状・float への一括の適用、時間の比を速さの率にすること。
"""
import argparse
import ast
import csv
import hashlib
import json
import math
import subprocess
import sys
import types
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
HOOKS = {}            # 自己試験用の差し替え (fields・compare)。本番では空
RUN = "run_0810_vgc_chain"
REF = {"T2a": "run_0752_t2_d1", "T2b": "run_0754_t2_d2"}
CTRL = {"S_g0c": "run_0804_sp_g0c"}
CAND = "C"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94")
BUDGET, INTERVAL, CONSEC, MIN_POINTS = 200000, 2500, 2, 9
TAIL = 20000
BIN_SHA = "7e95e5bbd8f50c9e"
MSH_SHA = "285f0fa4ecc0ed1c"
# run_case.sh が RUN_PROVENANCE.txt に書く FORGE_* (今は CUDA_BLOCKSIZE・KERNEL_SYNC・NODE_FX_HALF・AXIS_DIAG_ALPHA・RESID_SNAP・BIN) は、
# 登録 (sp_arm.sh と同じ) の値だけであること。CUDA_BLOCKSIZE=128 は必ずある。BIN は forge_sha256 で照合する
ENV_WANT = {"FORGE_CUDA_BLOCKSIZE": "128", "FORGE_OMEGA_BUDGET": "1", "FORGE_DIAG_COMMIT_LOSS": "500"}
# この commit の sp_judge.py・m9_watch.py の関数の定義の sha256 (先頭 16 桁)。B の写しが違えば、比較・到達の計算し直しは判定不能
SP_FUNCS = {"fields": "198a8d50ea0ed40b", "compare": "230a1a0def551c33"}
M9_FUNCS = {"drift": "72506550c522d561", "window_points": "d16897d68f8af794", "done_at": "2c66bdc29bd942f5", "spacing_ok": "a2e65d791a97ebe3"}


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def jload(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def jlines(p) -> list:
    out = []
    try:
        for ln in Path(p).read_text().splitlines():
            if ln.strip():
                out.append(json.loads(ln))
    except (OSError, ValueError):
        pass
    return out


def extract_funcs(path: Path, want: dict) -> tuple:
    """path の最上位の関数の定義のうち want の名前の文を取り出す。戻り ({名前: 文}, {名前: sha256 の先頭 16 桁}, [食い違い])。"""
    src = path.read_text()
    tree = ast.parse(src)
    seg = {n.name: ast.get_source_segment(src, n) for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in want}
    sha = {k: hashlib.sha256(v.encode()).hexdigest()[:16] for k, v in seg.items()}
    bad = [f"{path.name}:{k} (sha {sha.get(k)}、この commit {v})" for k, v in want.items() if sha.get(k) != v]
    return seg, sha, bad


def residual_nonfinite(p: Path):
    """残差の履歴 CSV の rms_* 列の非有限の数と最初の step。無ければ None。"""
    if not p.is_file():
        return None
    n_bad, first, rows, last = 0, None, 0, None
    with open(p, newline="") as f:
        r = csv.reader(f)
        head = next(r, [])
        cols = [i for i, c in enumerate(head) if c.strip().startswith("rms_")]
        for row in r:
            if not row:
                continue
            rows += 1
            last = row[0]
            for i in cols:
                try:
                    v = float(row[i])
                except (ValueError, IndexError):
                    v = float("nan")
                if not math.isfinite(v):
                    n_bad += 1
                    if first is None:
                        first = row[0]
    return {"file": p.name, "rows": rows, "rms_cols": len(cols), "nonfinite": n_bad, "first_nonfinite_step": first, "last_step": last}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", default=str(HERE), help="run が並ぶ case dir (既定: この台本の場所)")
    ap.add_argument("--code-dir", default=str(HERE), help="sp_judge.py・m9_watch.py・cold_xcheck.py の場所 (既定: この台本の場所)")
    ap.add_argument("--tools", default=None, help="forge の tools (既定: ~/forge-g0cchain の、無ければこの作業ツリーの)")
    ap.add_argument("--run", default=RUN)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    case, code = Path(a.case_dir).resolve(), Path(a.code_dir).resolve()
    tools = Path(a.tools) if a.tools else (Path.home() / "forge-g0cchain/solver_density_cuda/tools")
    if not tools.is_dir():
        tools = HERE.parents[1] / "solver_density_cuda/tools"
    out_p = Path(a.out) if a.out else case / "_band_ab/cold_pair/vgc_judge.json"
    run = case / a.run
    OUT = {"plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #31・§6 V-g0c-chain (事前登録 commit 68d56723)",
           "created": now(), "run": a.run, "case_dir": str(case), "code_dir": str(code), "tools": str(tools),
           "undecidable": [], "b_reasons": [], "gate_failures": [], "precheck": {}, "stages": {}, "watch": {}, "pairs": {},
           "control_record": {}, "quasisteady": {}, "timing": {}, "records": {}}
    und, bres, gfail = OUT["undecidable"], OUT["b_reasons"], OUT["gate_failures"]

    # --- 判定に使う関数 (sp_judge.py の fields・compare、m9_watch.py の到達の計算) ---
    sp_seg, sp_sha, sp_bad = extract_funcs(code / "sp_judge.py", SP_FUNCS)
    m9_seg, m9_sha, m9_bad = extract_funcs(code / "m9_watch.py", M9_FUNCS)
    OUT["records"]["functions"] = {"sp_judge.py": {"sha256": sha256_file(code / "sp_judge.py")[:16], "funcs": sp_sha, "mismatch": sp_bad},
                                   "m9_watch.py": {"sha256": sha256_file(code / "m9_watch.py")[:16], "funcs": m9_sha, "mismatch": m9_bad}}
    m9ns = {"np": np, "a": types.SimpleNamespace(consec=CONSEC, min_points=MIN_POINTS, interval=INTERVAL), "REF": None,
            "KEYS": KEYS + ("Q_w",)}
    if not m9_bad:
        exec(compile("\n\n".join(m9_seg[k] for k in M9_FUNCS), str(code / "m9_watch.py"), "exec"), m9ns)

    # --- 1. 前提の検査の記録・基盤の停止 ---
    pc = jload(run / "VGC_PRECHECK.json")
    if not run.is_dir():
        und.append(f"run {a.run} が無い")
    elif pc is None:
        und.append("前提の検査の記録 VGC_PRECHECK.json が無い")
    else:
        OUT["precheck"] = {"ok": pc.get("ok"), "failures": pc.get("failures"), "items_ok": {k: v.get("ok") for k, v in (pc.get("items") or {}).items()
                                                                                         if isinstance(v, dict) and "ok" in v},
                           "ar_gt5000": ((pc.get("items") or {}).get("3_quality") or {}).get("ar_gt5000"),
                           "sha": (pc.get("items") or {}).get("4_sha"), "code": pc.get("code")}
        if not pc.get("ok"):
            und.append("前提の検査が不成立: " + "; ".join(pc.get("failures") or ["(理由なし)"])[:600])
        elif not (run / "nozzle.msh").is_file() or not sha256_file(run / "nozzle.msh").startswith(MSH_SHA):
            und.append(f"判定の時点の nozzle.msh が G0c ({MSH_SHA}) でない (前提の検査 (1) の読み直し)")
    for f, why in (("DISK_STOP", "ディスクの空きで自分の forge を止めた"), ("WATCH_ABORT", "本段の直前の照合が不成立")):
        if (run / f).is_file():
            und.append(f"{why} ({f}: {(run / f).read_text().strip()[:300]})")
    rc_txt = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").is_file() else None
    run_end = jload(run / "RUN_END.json") or {}
    OUT["records"]["RUN_RC"] = rc_txt
    OUT["records"]["run_end_exception"] = (run_end.get("exception") or "")[-1500:] or None
    mss = jload(run / "MAIN_STAGE_STARTED")
    vms = jload(run / "VGC_MAIN_START.json")
    prov = (run / "RUN_PROVENANCE.txt").read_text() if (run / "RUN_PROVENANCE.txt").is_file() else ""
    if mss is not None:
        if prov:
            sha_l = next((l.split(":", 1)[1].strip() for l in prov.splitlines() if l.startswith("forge_sha256")), "")
            envs = dict(l.split(":", 1)[1].strip().split("=", 1) for l in prov.splitlines() if l.startswith("env ") and "=" in l)
            OUT["records"]["provenance"] = {"forge_sha256": sha_l[:16], "env": envs}
            env_bad = [k for k, v in envs.items() if k != "FORGE_BIN" and ENV_WANT.get(k) != v]
            if not sha_l.startswith(BIN_SHA) or envs.get("FORGE_CUDA_BLOCKSIZE") != "128" or env_bad:
                und.append(f"本段の RUN_PROVENANCE.txt のバイナリ・環境変数が登録と違う ({sha_l[:16]}、{envs})")
        else:
            und.append("本段の RUN_PROVENANCE.txt が無い")

    # --- 2. 段のゲート ---
    stg = jlines(run / "vgc_stages.jsonl")
    gates = jlines(run / "vgc_stage_gate.jsonl")
    ends = {r["tag"]: r for r in stg if r.get("event") == "end"}
    S = OUT["stages"]
    hist = {t: residual_nonfinite(run / f"residual_history_{t}.csv") for t in ("S1_soft", "S2_mid", "main")}
    if hist["main"] is None and mss is not None:
        hist["main"] = residual_nonfinite(run / "residual_history.csv")
    S["residual_history"] = hist
    for t in ("S1_soft", "S2_mid"):
        g = next((x for x in gates if x.get("tag") == t), None)
        e = ends.get(t)
        st = (g or {}).get("stats") or {}
        rec = {"rc": None if e is None else e.get("rc"), "seconds": None if e is None else e["t1"] - e["t0"],
               "gate_problems": None if g is None else g.get("problems"), "stats": st or (g or {}).get("stats_error"),
               "res_sha256": None if g is None else (g.get("res_sha256") or "")[:16]}
        ok = (e is not None and e.get("rc") == 0 and g is not None and not g.get("problems") and bool(st)
              and st.get("nonfinite_total") == 0 and all(st.get(f"{k}_n_nonpositive") == 0 and not st.get(f"{k}_missing") for k in ("ro", "P", "T")))
        rec["ok"] = ok
        S[t] = rec
        h = hist[t]
        if h is None or h["nonfinite"]:
            ok = False
        if not ok and mss is not None:
            gfail.append(f"段 {t}: 正常終了・有限性・ρ・P・T・残差の履歴のどれかが不成立 ({rec}, 残差 {h})")
    if mss is None:
        if run.is_dir() and pc and pc.get("ok") and not und:
            ev = []
            for g in gates:
                st = g.get("stats") or {}
                if g.get("problems") or st.get("nonfinite_total") or any(st.get(f"{k}_n_nonpositive") for k in ("ro", "P", "T")):
                    ev.append(f"{g.get('tag')}: 段終了ゲート {g.get('problems')[:3] if g.get('problems') else ''} 非有限 {st.get('nonfinite_total')}")
            for t in ("S1_soft", "S2_mid"):
                if hist[t] and hist[t]["nonfinite"]:
                    ev.append(f"{t}: 残差の非有限 {hist[t]['nonfinite']} (最初の step {hist[t]['first_nonfinite_step']})")
            nan = sorted(f.name for f in run.glob("res_nan_*.h5"))
            if nan:
                ev.append(f"detectNaN の出力 {nan[:3]}")
            fl = run / "forge_run.log"
            if fl.is_file() and "[detectNaN]" in fl.read_text(errors="replace"):
                ev.append("forge_run.log に [detectNaN]")
            S["prestage_divergence_evidence"] = ev
            if ev:
                bres.append("前の段で発散し本段に入れない: " + "; ".join(ev)[:500])
            else:
                und.append(f"本段が始まらなかった (前の段の失敗に発散の証拠が無い; RUN_RC {rc_txt}; 例外 {(run_end.get('exception') or '')[-300:]!r})")
    else:
        if vms is None or not vms.get("ok"):
            und.append(f"本段の直前の照合の記録 VGC_MAIN_START.json が無いか不成立 ({None if vms is None else vms.get('failures')})")
        if mss.get("res_files_present"):
            und.append(f"本段の通知の時点で前の段の出力が残っていた: {mss.get('res_files_present')}")
        S["main_stage_started"] = {k: mss.get(k) for k in ("time", "time_unix", "main_config", "initial_field", "previous_stages", "prestage_config")}
        if hist["main"] is None:
            gfail.append("本段の残差の履歴が無い")
        elif hist["main"]["nonfinite"]:
            gfail.append(f"本段の残差の履歴に非有限 {hist['main']['nonfinite']} (最初の step {hist['main']['first_nonfinite_step']})")
        # stage_manifest の区間
        try:
            sys.path.insert(0, str(tools))
            import stage_manifest as SM
            man = json.loads((run / "stage_manifest.json").read_text())
            segs = SM.segments(man, SM.load_launches(str(run)))
            tags = [[s["tag"] for s in sg] for sg in segs]
            S["manifest_segments"] = tags
            if not tags or tags[-1] != ["main"]:
                gfail.append(f"stage_manifest の最後の区間が本段だけでない: {tags}")
        except Exception as e:  # noqa: BLE001
            gfail.append(f"stage_manifest の区間を読めない ({e!r})")
        p = subprocess.run([sys.executable, str(tools / "stage_manifest.py"), str(run), "--segments"], capture_output=True, text=True)
        S["stage_manifest_cli"] = {"rc": p.returncode, "out": (p.stdout + p.stderr).strip().splitlines()[-12:]}
        p = subprocess.run([sys.executable, str(tools / "check_convergence.py"), str(run), "--segment"], capture_output=True, text=True)
        lines = (p.stdout + p.stderr).splitlines()
        S["check_convergence_segment"] = {"rc": p.returncode, "overall": next((l for l in lines if "OVERALL" in l), "(行なし)"),
                                          "verdict_lines": [l.strip()[:160] for l in lines if "VERDICT" in l or "<--" in l][:20]}
        # (6) の照合し直し (記録と同じであること)
        if vms and vms.get("files_6_late"):
            S["files_6_late"] = vms["files_6_late"].get("files")

    # --- 3. 本段の見張り ---
    W = OUT["watch"]
    m9 = jload(run / "m9_watch.json") if run.is_dir() else None
    status = None if m9 is None else m9.get("status")
    W.update(status=status, reach_step=None if m9 is None else m9.get("reach_step"), note=None if m9 is None else m9.get("note"),
             rows=None if m9 is None else len(m9.get("rows", [])), RUN_RC=rc_txt,
             termination=("見張りの停止 (REACHED)" if status == "REACHED" else "上限まで回った (CENSORED)" if status == "CENSORED"
                          else f"予期しない終了 ({status}, RUN_RC {rc_txt})" if status else "見張りの記録なし"))
    reached = False
    if mss is not None:
        rows = sorted((m9 or {}).get("rows", []), key=lambda r: r["step"])
        if m9 is None:
            und.append("見張りの記録 m9_watch.json が無い")
        elif m9.get("phase") != "line" or m9.get("budget") != BUDGET or m9.get("consec") != CONSEC or m9.get("interval") != INTERVAL \
                or m9.get("min_points") != MIN_POINTS or m9.get("keep_tail") != 9:
            und.append(f"見張りの引数が §6.11 と違う ({ {k: m9.get(k) for k in ('phase', 'budget', 'consec', 'interval', 'min_points', 'keep_tail')} })")
        elif m9_bad and status in ("REACHED", "CENSORED"):
            und.append("m9_watch.py の到達の計算がこの commit の版と違う: " + "; ".join(m9_bad))
        elif status == "REACHED":
            n = m9.get("reach_step")
            sp_ok = m9ns["spacing_ok"](rows)
            d = m9ns["done_at"](rows)
            W["recheck"] = {"spacing_ok": sp_ok, "done_at": d}
            fpath = run / f"res_{n}.h5"
            if not sp_ok or d != ("REACH", n):
                und.append(f"到達を系列から計算し直すと合わない (見張りの誤作動): 間隔 {sp_ok}、done_at {d}、記録 {n}")
            elif not fpath.is_file():
                und.append(f"到達の出力 {fpath.name} が無い")
            elif sha256_file(fpath) != m9.get("reach_sha256"):
                und.append(f"到達の出力 {fpath.name} の sha256 が見張りの記録と違う")
            else:
                reached = True
        elif status == "CENSORED":
            sp_ok = m9ns["spacing_ok"](rows)
            d = m9ns["done_at"](rows)
            last = rows[-1]["step"] if rows else None
            W["recheck"] = {"spacing_ok": sp_ok, "done_at": d, "last_step": last}
            if not sp_ok or d is not None or last is None or last < BUDGET or rc_txt != "0":
                und.append(f"CENSORED を系列から確かめられない (見張りの誤作動): 間隔 {sp_ok}、done_at {d}、最後の step {last}、RUN_RC {rc_txt}")
            else:
                bres.append(f"本段 {BUDGET} step までに停止規則に到達しない (CENSORED)")
        elif status == "DIVERGED":
            bres.append(f"本段で発散 (DIVERGED、step {m9.get('fail_step')}、{m9.get('note') or ''})")
        else:
            und.append(f"見張りの状態 {status} (予期しない終了・見張りの誤作動は数値の手法の失敗に数えない; RUN_RC {rc_txt}, {m9.get('note') or m9.get('rc')})")

    # --- 4. M1〜M7 (到達したときだけ) ---
    F = {}
    if reached and not und:
        if sp_bad:
            und.append("sp_judge.py の fields・compare がこの commit の版と違う: " + "; ".join(sp_bad))
        else:
            if "fields" in HOOKS and "compare" in HOOKS:          # 自己試験だけ (合成の入力で分岐を確かめる)
                fields = HOOKS["fields"]
                compare = (lambda c, r, _f=HOOKS["compare"]: _f(F, c, r))
            else:
                sys.path.insert(0, str(code))
                import h5py
                import cold_xcheck as XC
                yb_x, yb = XC.common_yb()
                ns = {"np": np, "h5py": h5py, "XC": XC, "HERE": case, "yb_x": yb_x, "yb": yb, "KEYS": KEYS, "F": F}
                exec(compile(sp_seg["fields"] + "\n\n" + sp_seg["compare"], str(code / "sp_judge.py"), "exec"), ns)
                fields, compare = ns["fields"], ns["compare"]
            arms = {CAND: (a.run, int(m9["reach_step"]))}
            for nm, rr in REF.items():
                rm = jload(case / rr / "m9_watch.json")
                if not rm or rm.get("status") != "REACHED":
                    und.append(f"参照 {nm} ({rr}) の見張りの状態 {None if not rm else rm.get('status')}")
                    continue
                arms[nm] = (rr, int(rm["reach_step"]))
            for nm, (rr, n) in arms.items():
                try:
                    F[nm] = fields(rr, n)
                    if F[nm]["nonfinite"]:
                        und.append(f"{nm} ({rr}) の停止の出力に非有限 {F[nm]['nonfinite']}")
                except (OSError, KeyError, SystemExit, ValueError) as e:
                    und.append(f"{nm} ({rr}) の停止の出力を読めない ({e!r})")
            OUT["watch"]["arms"] = {k: {"run": v[0], "step": v[1]} for k, v in arms.items()}
            if not und:
                for r in REF:
                    OUT["pairs"][f"{CAND}/{r}"] = compare(CAND, r)
                over = [k for k, v in OUT["pairs"].items() if not v["許容内"]]
                if over:
                    bres.append(f"M1〜M7 の許容の超過: {over}")
                # 対照 run_0804 との比較 (記録だけ。近さだけで合格にしない)
                for nm, rr in CTRL.items():
                    try:
                        rm = jload(case / rr / "m9_watch.json")
                        F[nm] = fields(rr, int(rm["reach_step"]))
                        OUT["control_record"][f"{CAND}/{nm}"] = compare(CAND, nm)
                    except Exception as e:  # noqa: BLE001
                        OUT["control_record"][f"{CAND}/{nm}"] = f"記録できない ({e!r})"

    # --- 5. 準定常 (記録) ---
    if m9 is not None and status in ("REACHED", "CENSORED") and mss is not None:
        n = int(m9["reach_step"]) if status == "REACHED" else BUDGET
        rows = sorted([x for x in m9.get("rows", []) if n - TAIL <= x["step"] <= n], key=lambda x: x["step"])
        Q = OUT["quasisteady"]
        Q["window"] = [n - TAIL, n]
        Q["points"] = len(rows)
        if len(rows) != 9:
            Q["note"] = f"末尾 2 万 step の出力が 9 でない ({len(rows)})"
        sc = run / "vgc_series_tail.csv"
        try:
            with open(sc, "w") as fh:
                fh.write("step," + ",".join(KEYS + ("Q_w",)) + "\n")
                for x in rows:
                    fh.write(f"{x['step']}," + ",".join(f"{x[k]:.10g}" for k in KEYS + ("Q_w",)) + "\n")
            for tag, cols, thr in (("§6.11 θ_r (0.0001)", KEYS, "0.0001"), ("§6.11 Q_w (0.0002)", ("Q_w",), "0.0002"),
                                   ("sp_judge と同じ (4 量 0.0001)", KEYS + ("Q_w",), "0.0001")):
                q = subprocess.run([sys.executable, str(tools / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", ",".join(cols),
                                    "--tail", "1.0", "--drift", thr, "--osc", thr], capture_output=True, text=True)
                Q[tag] = [l.strip() for l in q.stdout.splitlines() if any(c in l for c in cols) and ":" in l] or (q.stdout + q.stderr).strip().splitlines()[-5:]
        except Exception as e:  # noqa: BLE001
            Q["error"] = repr(e)

    # --- 6. 時間 (記録) ---
    try:
        t0 = float((run / "START_EPOCH").read_text().strip())
        tm = jload(run / "MAIN_STAGE_STARTED") or {}
        t_main = float(tm["time_unix"])
        rows_t = [r for r in csv.reader(open(run / "timing.tsv"), delimiter="\t") if len(r) >= 3]
        allr = [(float(r[0]), int(r[2]), (r[3].strip() if len(r) > 3 else "")) for r in rows_t]
        mr = [(float(r[0]), int(r[1])) for r in rows_t if r[1].strip() and float(r[0]) >= t_main]
        n_end = int(m9["reach_step"]) if status == "REACHED" else (BUDGET if status == "CENSORED" else None)
        T = OUT["timing"]
        T.update(prepare_s=(pc or {}).get("t_prepared", t0) - t0, main_start_after_s=t_main - t0, samples=len(allr), main_samples=len(mr))
        for e in ends.values():
            T[f"stage_{e['tag']}_s"] = e["t1"] - e["t0"]
        if n_end and mr:
            t = np.array([x[0] for x in mr]); s = np.array([x[1] for x in mr])
            if s.max() >= n_end:
                k = int(np.argmax(s >= n_end))
                tr = t[k] if k == 0 or s[k] == s[k - 1] else t[k - 1] + (t[k] - t[k - 1]) * (n_end - s[k - 1]) / (s[k] - s[k - 1])
                j = int(np.argmax(s >= INTERVAL))
                T.update(n_end=n_end, T_total_s=float(tr - t0), T_main_s=float(tr - t_main),
                         per_step_ms=float(1e3 * (tr - t[j]) / max(n_end - s[j], 1)),
                         exclusive=bool(all(nf <= 1 for tt, nf, _ in allr if tt <= tr)),
                         others=sorted({o for tt, _, o in allr if o and tt <= tr}))
            else:
                T["note"] = f"timing.tsv の本段の行が到達の step {n_end} まで無い (最大 {int(s.max())})"
    except Exception as e:  # noqa: BLE001
        OUT["timing"]["error"] = repr(e)

    # --- 判定 ---
    if und:
        verdict = "判定不能: " + "; ".join(und)[:900]
    elif bres:
        verdict = "結果 B (拡大を退け、いまの採用の範囲 = 登録の Q32 の起点を保つ): " + "; ".join(bres)[:900]
    elif gfail:
        verdict = "判定不能 (到達・比較は A の条件だが、段のゲートに不成立があるので A と書かない): " + "; ".join(gfail)[:900]
    elif reached and len(OUT["pairs"]) == len(REF) and all(v["許容内"] for v in OUT["pairs"].values()):
        verdict = "結果 A: 登録した開始の手順で停止規則に到達し、T2 の両方の参照との比較が許容内 (この開始の手順への適用の拡大を支持)"
    else:
        verdict = "判定不能 (想定外の分岐: 到達・比較の記録がそろわない)"
    OUT["verdict"] = verdict
    for k, v in OUT["pairs"].items():
        print(k, "許容内" if v["許容内"] else "超過", {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "許容内"})
    print("時間", json.dumps(OUT["timing"], ensure_ascii=False, default=float))
    print(f"== VERDICT V-g0c-chain: {verdict}")
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    print(f"→ {out_p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
