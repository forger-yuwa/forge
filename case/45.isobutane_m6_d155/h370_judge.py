"""plan tooling-nozzle-isothermal-wall-chain §5.1 #36f・§6「V-h370 の事前登録」(§4.9 codex diagnose 3) の判定。B の case dir で回す
(h370_run.py が呼ぶ)。壁は更新しない (固定壁の pass)。

--euler-qualify: (i) Euler 参照 run_0950_h370_euler_ref の資格 (§6 V-h370 の「Euler」):
  - RUN_RC 0、投入前の検査 (H370_PRECHECK.json) の合格、本段が上限 60,000 step まで回った (最後の出力 res_60000.h5)。
  - NaN/Inf なし: 全段の残差の履歴 (residual_history_{soft,mid,main}.csv) の rms_* と、本段の全出力の VALUE の浮動小数点の量。ρ・P・T > 0。
  - 全残差: check_convergence.py --segment --drop 3 (stage_manifest の最後の区間 = 方程式・BC・離散化が同じ本段の区間、全列) の PASS。
  - 出口 M・流量・参照場の時系列 (本段の全出力) の check_quasisteady.py の VERDICT (CSV 系列モード、閾値は既定のまま): 資格に使う量は
    出口コア M・無次元の流量 (stop_rule.recorded_quantities = nozzle_report.metrics・qs_series.py と同じ定義) と、参照場として軸 (η = 0) と
    η = 0.1 の線の M を x = 40・70・94 r_t で (nozzle_report.eta_line、停止規則の断面と同じ位置)。すべて STEADY であること
    (OSCILLATING も参照の時間の変動なので資格未達)。nozzle_report.metrics の η = 0.1 の線の評価量 (wave・overshoot・slope・range、%) は
    平均が 0 に近く相対の閾値が意味を持たないので記録だけ。閾値を細かくした (1e-4) 判定も記録だけ。
  満たせば本段の最後の出力とその sha256 を _band_ab/half_inlet/h370_euler_ref_pin.json に固定する (既にあって違えば止める; 消さない)。
  参照の時間の変動を後で測るための出力 (check_quasisteady の末尾の窓の始まりと中ほど) と sha も記録する。
  → _band_ab/half_inlet/h370_euler_qualify.json。終了コード 0 資格あり / 1 資格未達 / 2 判定不能。

(既定): (ii) NS run_0951_h370_ns_b0 の §6 V-h370 の判定 → _band_ab/half_inlet/h370_judge.json:
  1. 記録: 投入前の検査 (H370_PRECHECK.json)・本段の直前の照合 (H370_MAIN_START.json)・基盤の停止 (PRECHECK_FAILED・WATCH_ABORT・
     RESOURCE_STOP・DISK_STOP・WRITE_FAILED)・RUN_PROVENANCE.txt のバイナリと環境変数・段の残差の履歴の非有限。
  2. 停止の状態 (STOP_RULE_STATE.json): 状態・停止の出力 s・到達と中断と保持 (stop_rule.json の hold)・要求の断面と実際の列・停止の出力の sha。
  3. 受け渡しの関門 (deltastar_loop.handoff_gate、HANDOFF_GATE.json): 停止の状態 REACHED・check_convergence.py --segment --drop 3 の PASS・
     保存則・局所 y₁⁺ (check_wall_resolution.py --step s)・2 窓 [s − 40000, s − 20000]・[s − 20000, s] の δ_E・δ_target (全使用域で
     最大 − 最小 ≤ 0.01 mm、窓の平均の差 ≤ 0.01 mm)・Euler の参照と ω・抽出の全域の被覆・使用域の held と欠測なし・次の壁の幾何。
     上書き (allow_not_reached・allow_unconverged) は使わない。合格なら関門が抽出の表を run の直下に写す (次の pass は自動でつながない)。
  4. check_quasisteady.py の VERDICT (記録; 既定の相対閾値で寸法の予算の合格を代用しない): 2 窓の r・θ_r (断面ごと)・出口コア M・流量・
     δ_E・δ_target (x = 40・70・94 r_t) の時系列。
  5. 抽出の参照の A/B (同じ NS の保存系列に対して Euler の参照だけを変える): 腕 A = 旧い run_0174_euler_v5d_M_r1 (最後の出力を sha で
     固定)、腕 B = 固定した run_0950。2 窓の全出力で受け渡しと同じ条件 (edge・一定 x の断面への補間・端は欠測・ω 0.5・平滑化の λ・
     昇格の上限; stop_rule.json の extract) の δ_E を比べる。比べるのは腕 A が実際に抽出した列と腕 B が抽出した列の共通の使用域だけ
     (旧い参照の範囲の外 = 出口の区間は抽出しないので比べない; 外挿しない)。
       D_k = 共通の列の max|δ_E,B − δ_E,A| [mm] (出力 k ごと)、合否 = 抽出の例外なし・共通の列に held なし・δ_E・δ_target が有限。
       σ_ref = 腕 B の参照の時間の変動: 停止の出力 s を、固定した出力の代わりに資格の記録の別の出力 (末尾の窓の始まり・中ほど) の Euler で
       抽出したときの max|Δδ_E| [mm]。
       結果 A: 全出力で合否が同じ かつ max_k D_k + σ_ref ≤ 0.01 mm。結果 B: どれかの出力で合否が違う、または min_k D_k − σ_ref > 0.01 mm。
       それ以外 (窓の中で D_k が 0.01 mm をまたぐ = 時間の変動、σ_ref で 0.01 mm をまたぐ = 参照の未整定) と抽出の失敗は判定不能。
     停止の出力 s では両腕の全体の抽出 (run の H370_AB/{A,B}) で帯の位置 (y_b・y*)・感度・held・δ_target の差・次の壁 (腕 B の表、腕 A は
     共通の列だけ腕 A の δ_target に替えた表 = 出口の区間を延ばさない) の真のスロートを記録する。
  6. 設計を終える条件の量 (記録だけ; 最初の pass の判定には使わない): 使用域の未緩和の max|δ_E − δ_in|。
  判定: 「更新に進める条件」= 関門の proceed (と 1 の基盤の不成立が無いこと)、抽出の参照の A/B は結果 A / 結果 B / 判定不能。
usage: python3 h370_judge.py --euler-qualify     python3 h370_judge.py [--no-ab]
       (--case-dir・--design・--tools は試験用; 既定は B の ~/forge-g0cchain、無ければこの作業ツリー)
"""
import argparse
import contextlib
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
HOOKS = {}            # 試験用の差し替え (handoff_gate・delta_snapshot・extract_and_merge・next_wall_check)。本番では空
PLAN = "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36f・§6 V-h370 (§4.9 codex diagnose 3)"
EU_RUN, NS_RUN = "run_0950_h370_euler_ref", "run_0951_h370_ns_b0"
NS_PROBLEM = "problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml"
OLD_EULER = "run_0174_euler_v5d_M_r1"
PIN = "_band_ab/half_inlet/h370_euler_ref_pin.json"
QUAL = "_band_ab/half_inlet/h370_euler_qualify.json"
OUT = "_band_ab/half_inlet/h370_judge.json"
EU_STEPS = 60000
BIN_SHA = "7e95e5bbd8f50c9e"
XS = (40.0, 70.0, 94.0)                       # 参照場・記録の量の位置 [r_t] (停止規則の要求の断面 = x_fracs × x_w)
QUAL_COLS = ["exit_core_M", "mdot_rt"] + [f"M_axis_{x:g}" for x in XS] + [f"M_eta0.1_{x:g}" for x in XS]
REC_COLS = ["wave_eta0.1", "overshoot_eta0.1", "slope_eta0.1", "range_eta0.1"]
QS_TIGHT = ["--drift", "0.0001", "--osc", "0.0001"]       # 記録だけ (既定の tail 0.4)
VERDICT_RE = re.compile(r"^\s+(\S+)\s*:\s(.*?)\s*(STEADY|DRIFTING|OSCILLATING|TRANSIENT-UNSETTLED|NONFINITE)\s*$")
TOL_MM = 0.01                                 # 抽出の参照の A/B (§6 V-h370) と 2 窓の予算
AB_DIR = "H370_AB"
MARKERS = ("PRECHECK_FAILED", "WATCH_ABORT", "RESOURCE_STOP", "DISK_STOP", "WRITE_FAILED")


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


def jwrite(p: Path, d) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(p) + ".tmp")
    tmp.write_text(json.dumps(d, indent=1, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    os.replace(tmp, p)


def res_steps(d: Path) -> list:
    return sorted((int(f.stem.split("_")[1]) for f in Path(d).glob("res_[0-9]*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)))


def residual_nonfinite(p: Path):
    """残差の履歴 CSV の rms_* 列の非有限の数と最初の step (tw370_judge.py と同じ)。無ければ None。"""
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


def run_qs(tools: Path, series_csv: Path, cols: list, extra=()) -> dict:
    """check_quasisteady.py の CSV 系列モード。列ごとの VERDICT を出力の行から読む (行の無い列は missing)。"""
    cmd = [sys.executable, str(tools / "check_quasisteady.py"), "--series-csv", str(series_csv), "--series-cols", ",".join(cols), *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    v = {}
    for ln in r.stdout.splitlines():
        m = VERDICT_RE.match(ln)
        if m and m.group(1) in cols:
            v[m.group(1)] = {"verdict": m.group(3), "detail": m.group(2).strip()}
    return {"cmd": " ".join(cmd[1:]), "rc": r.returncode, "verdicts": v, "missing": [c for c in cols if c not in v],
            "overall": next((ln.strip() for ln in r.stdout.splitlines() if ln.startswith("OVERALL")), None),
            "stderr": r.stderr.strip()[-500:]}


def write_series(path: Path, rows: list, cols: list) -> None:
    with open(path, "w") as fh:
        fh.write("step," + ",".join(cols) + "\n")
        for r in rows:
            fh.write(f"{r['step']}," + ",".join(f"{float(r.get(c, float('nan'))):.12g}" for c in cols) + "\n")


# ===================================================================================================== (i) Euler の資格
def euler_series(run: Path, steps: list) -> tuple:
    """本段の出力ごとの量 (資格の量 QUAL_COLS・記録の量 REC_COLS) と、全 VALUE の非有限・ρ・P・T の最小。戻り (行, 場の問題)。"""
    import h5py
    from forge_design.evaluate import stop_rule as SR
    from forge_design.report import nozzle_report as NR
    info = json.loads((run / "prepare_info.json").read_text())
    ni, nj, S, x_E = int(info["mesh"]["ni"]), int(info["mesh"]["nj"]), float(info["scale_m"]), float(info["x_E"])
    with h5py.File(run / "nozzle.h5", "r") as h:
        nc = h["/MESH/COORD"][:].reshape(-1, 3)
    Xf, Rf = (nc[:, 0] / S).reshape(ni, nj), (nc[:, 1] / S).reshape(ni, nj)       # stop_rule.RunGeometry の Xf・Rf と同じ式
    rows, probs = [], []
    for st in steps:
        name = f"res_{st}.h5"
        row = {"step": st}
        with h5py.File(run / name, "r") as h:
            V = {k: np.asarray(d[:]) for k, d in h["VALUE"].items()}
        nf = {k: int(np.count_nonzero(~np.isfinite(a))) for k, a in V.items() if a.dtype.kind == "f"}
        row["nonfinite_all"] = int(sum(nf.values()))
        for k in ("ro", "P", "T"):
            if k not in V:
                probs.append(f"{name}: VALUE/{k} が無い")
                continue
            row[f"{k}_min"] = float(np.nanmin(V[k])) if V[k].size else float("nan")
        if row["nonfinite_all"]:
            probs.append(f"{name}: 非有限 {row['nonfinite_all']} 個 ({ {k: n for k, n in nf.items() if n} })")
            rows.append(row)
            continue
        if any(row.get(f"{k}_min", 1.0) <= 0 for k in ("ro", "P", "T")):
            probs.append(f"{name}: ρ・P・T の最小 {[row.get(f'{k}_min') for k in ('ro', 'P', 'T')]} (> 0 であること)")
        sh = (ni, nj)
        row.update(SR.recorded_quantities(Xf, Rf, {k: V[k].reshape(sh) for k in ("ro", "Ux", "Uy", "sonic")}, x_E))
        F = NR.load_field(run, name)
        for x in XS:
            row[f"M_axis_{x:g}"] = float(NR.eta_line(F, "M", 0.0, np.array([x]))[0])
            row[f"M_eta0.1_{x:g}"] = float(NR.eta_line(F, "M", 0.1, np.array([x]))[0])
        try:
            m = NR.metrics(run, F)["eta0.1"]
            row.update({"wave_eta0.1": m["wave_pct"], "overshoot_eta0.1": m["overshoot_pct"], "slope_eta0.1": m["slope_pct"],
                        "range_eta0.1": m["range_pct"]})
        except Exception as e:  # noqa: BLE001 — 記録だけの量
            row["metrics_error"] = repr(e)
        rows.append(row)
    return rows, probs


def euler_qualify(case: Path, tools: Path) -> int:
    run = case / EU_RUN
    R = {"plan": PLAN, "created": now(), "run": EU_RUN, "tools": str(tools), "undecidable": [], "failures": [], "qualified": False,
         "definitions": {"qualifying_cols": QUAL_COLS, "record_cols": REC_COLS, "positions_rt": list(XS),
                         "quasisteady": "check_quasisteady.py --series-csv (既定の --tail 0.4・--drift 0.05・--osc 0.10・--min-snaps 4)、全列 STEADY",
                         "exit_core_M・mdot_rt": "forge_design.evaluate.stop_rule.recorded_quantities (nozzle_report.metrics・qs_series.py と同じ定義)",
                         "M_axis・M_eta0.1": "forge_design.report.nozzle_report.eta_line(F, 'M', η, x) (η = 0 は軸、0.1 は metrics の判定の線)",
                         "record_cols": "nozzle_report.metrics の eta0.1 (wave_pct・overshoot_pct・slope_pct・range_pct)、記録だけ",
                         "pin": "本段の最後の出力 (上限 60,000 step) とその sha256"}}
    und, fails = R["undecidable"], R["failures"]
    qual_p, pin_p = case / QUAL, case / PIN

    def finish(rc):
        R["rc"] = rc
        jwrite(qual_p, R)
        print(f"== Euler の資格 ({EU_RUN}): {'資格あり → ' + PIN if rc == 0 else ('資格未達' if rc == 1 else '判定不能')}"
              + (f" — {'; '.join(und + fails)[:1500]}" if und or fails else ""))
        print(f"→ {qual_p}")
        return rc
    if not run.is_dir() or not (run / "prepare_info.json").is_file() or not (run / "nozzle.h5").is_file():
        und.append(f"{EU_RUN} が無いか prepare_info.json・nozzle.h5 が無い")
        return finish(2)
    rc_txt = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").is_file() else None
    R["run_rc"] = rc_txt
    if rc_txt != "0":
        fails.append(f"RUN_RC {rc_txt!r} (0 であること)")
    pc = jload(run / "H370_PRECHECK.json") or {}
    R["precheck_ok"] = pc.get("ok")
    if not pc.get("ok"):
        fails.append("投入前の検査 (H370_PRECHECK.json) が無いか不成立")
    for m in MARKERS:
        if (run / m).is_file():
            fails.append(f"基盤の停止・不成立の印 {m}: {(run / m).read_text()[:200]}")
    steps = [s for s in res_steps(run) if s > 0]
    R["outputs"] = {"n": len(steps), "first": steps[0] if steps else None, "last": steps[-1] if steps else None}
    if not steps:
        und.append("本段の出力 res_<N>.h5 (N > 0) が無い")
        return finish(2)
    if steps[-1] != EU_STEPS:
        fails.append(f"本段の最後の出力 res_{steps[-1]}.h5 が上限 {EU_STEPS} step でない")
    R["residual_nonfinite"] = {t: residual_nonfinite(run / f"residual_history_{t}.csv") for t in ("soft", "mid", "main")}
    for t, v in R["residual_nonfinite"].items():
        if v is None:
            fails.append(f"residual_history_{t}.csv が無い")
        elif v["nonfinite"]:
            fails.append(f"residual_history_{t}.csv の rms_* に非有限 {v['nonfinite']} 個 (最初 step {v['first_nonfinite_step']})")
    # 全残差 (同じ設定の区間)
    cmd = [sys.executable, str(tools / "check_convergence.py"), "--segment", "--drop", "3", str(run)]
    cv = subprocess.run(cmd, capture_output=True, text=True)
    (run / "H370_EULER_CONVERGENCE.txt").write_text(cv.stdout + cv.stderr)
    R["convergence"] = {"cmd": " ".join(cmd[1:]), "rc": cv.returncode,
                        "verdict": [ln.strip() for ln in cv.stdout.splitlines() if ln.startswith("===") or ln.startswith("OVERALL")]}
    if cv.returncode != 0:
        fails.append(f"check_convergence.py --segment --drop 3 が PASS でない (rc {cv.returncode}: {R['convergence']['verdict'][-1:] })")
    # 時系列
    try:
        rows, probs = euler_series(run, steps)
    except Exception as e:  # noqa: BLE001
        und.append(f"出力の量を作れない ({type(e).__name__}: {e})")
        return finish(2)
    fails += probs
    sc = run / "h370_euler_series.csv"
    write_series(sc, [r for r in rows if not r.get("nonfinite_all")], QUAL_COLS + REC_COLS)
    q = run_qs(tools, sc, QUAL_COLS)
    R["quasisteady"] = q
    R["quasisteady_record"] = {"record_cols (既定の閾値)": run_qs(tools, sc, REC_COLS), "tight (1e-4、記録だけ)": run_qs(tools, sc, QUAL_COLS + REC_COLS, QS_TIGHT)}
    if q["missing"]:
        und.append(f"check_quasisteady.py の出力に列 {q['missing']} の VERDICT が無い ({q['stderr'][-200:]})")
    bad = {c: v["verdict"] for c, v in q["verdicts"].items() if v["verdict"] != "STEADY"}
    if bad:
        fails.append(f"check_quasisteady.py の VERDICT が STEADY でない: {bad}")
    last = rows[-1]
    R["last_values"] = {c: last.get(c) for c in QUAL_COLS + REC_COLS}
    # 参照の時間の変動を後で測る出力 (check_quasisteady の末尾の窓 = 既定 tail 0.4 の始まりと中ほど)
    n = len(steps)
    k = max(3, int(math.ceil(0.4 * n)))
    tail = steps[-k:]
    alt = sorted({tail[0], tail[len(tail) // 2]} - {steps[-1]})
    R["sensitivity_snapshots"] = [{"res": f"res_{s}.h5", "step": s, "sha256": sha256_file(run / f"res_{s}.h5")} for s in alt]
    R["quasisteady_tail_window"] = [tail[0], tail[-1]]
    if und:
        return finish(2)
    if fails:
        return finish(1)
    pin = {"plan": PLAN, "created": now(), "run": EU_RUN, "run_abs": str(run.resolve()), "res": f"res_{steps[-1]}.h5", "step": steps[-1],
           "sha256": sha256_file(run / f"res_{steps[-1]}.h5"), "qualify": QUAL,
           "note": "NS の euler_ref は {run: <case dir>/run_0950_h370_euler_ref, res, sha256} (h370_run.py read_pin)"}
    R["pin"] = {k: pin[k] for k in ("run", "res", "step", "sha256")}
    old = jload(pin_p)
    if old is not None:
        if (old.get("run"), old.get("res"), old.get("sha256")) != (pin["run"], pin["res"], pin["sha256"]):
            und.append(f"{PIN} が既にあり、今回の固定 ({pin['res']} {pin['sha256'][:16]}) と違う ({old.get('res')} {str(old.get('sha256'))[:16]}) — 消さない")
            return finish(2)
        R["pin_note"] = "既にある同じ固定をそのまま使う"
    else:
        pin_p.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(pin_p), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, json.dumps(pin, indent=1, ensure_ascii=False).encode())
        os.close(fd)
    R["qualified"] = True
    return finish(0)


# ===================================================================================================== (ii) NS
@contextlib.contextmanager
def snapshot_view(run: Path, step: int):
    """run の prepare_info.json・nozzle.h5・res_<step>.h5 だけを symlink した一時ディレクトリ (抽出器が最新の出力を読むので 1 出力に絞る)。"""
    with tempfile.TemporaryDirectory(prefix="h370_view_") as td:
        v = Path(td) / run.name
        v.mkdir()
        for name in ("prepare_info.json", "nozzle.h5", f"res_{step}.h5"):
            (v / name).symlink_to((run / name).resolve())
        yield v


def _arr(d, k):
    return np.asarray(d.get(k, []), dtype=float)


def ab_pair(step: int, dA: dict, dB: dict, mm: float, src_b: str, x_use) -> dict:
    """出力 1 つの腕 A・B の比較。共通の列 = 両腕とも抽出した列 (delta_snapshot は抽出しなかった列を NaN にする)。"""
    r = {"step": step, "B_source": src_b}
    for arm, d in (("A", dA), ("B", dB)):
        if "error" in d or "delta_E" not in d:
            r[f"error_{arm}"] = d.get("error", "δ_E が無い")
    if any(k.startswith("error") for k in r):
        return r
    eA, eB, tA, tB = _arr(dA, "delta_E"), _arr(dB, "delta_E"), _arr(dA, "delta_target"), _arr(dB, "delta_target")
    hA, hB = np.asarray(dA.get("held", []), dtype=bool), np.asarray(dB.get("held", []), dtype=bool)
    if not (len(eA) == len(eB) == len(tA) == len(tB) == len(hA) == len(hB)):
        r["error_shape"] = "腕 A と B の列の数が違う"
        return r
    com = np.isfinite(eA) & np.isfinite(eB)
    r["n_common"], r["n_excluded_A"], r["n_excluded_B"] = int(com.sum()), int((~np.isfinite(eA)).sum()), int((~np.isfinite(eB)).sum())
    xu = np.asarray(x_use, dtype=float)
    if len(xu) == len(com) and com.any():
        r["x_common"] = [float(xu[com].min()), float(xu[com].max())]
        exA = ~np.isfinite(eA)
        r["x_excluded_A"] = [float(xu[exA].min()), float(xu[exA].max())] if exA.any() else None
    if not com.any():
        r["error_common"] = "共通の列が無い"
        return r
    d = np.abs(eB - eA)[com] * mm
    i = int(np.argmax(d))
    r["D_mm"], r["i_at_D"] = float(d[i]), int(np.flatnonzero(com)[i])
    r["Dt_mm"] = float(np.nanmax(np.abs(tB - tA)[com]) * mm)
    for arm, t, h in (("A", np.asarray(tA), hA), ("B", np.asarray(tB), hB)):
        r[f"pf_{arm}"] = bool(not h[com].any() and np.all(np.isfinite(t[com])))
        r[f"held_common_{arm}"] = int(h[com].sum())
    return r


def classify_ab(per: list, sigma, tol: float = TOL_MM) -> tuple:
    """A/B の結果 (§6 V-h370): 戻り (結果, 理由)。sigma は参照の時間の変動 [mm] (None = 測れない)。"""
    err = [r["step"] for r in per if any(k.startswith("error") for k in r)]
    if not per:
        return "判定不能", "比べる出力が無い"
    if err:
        return "判定不能", f"抽出・比較の失敗 (出力 {err[:6]})"
    if sigma is None:
        return "判定不能", "参照の時間の変動 σ_ref を測れない"
    D = [r["D_mm"] for r in per]
    mis = [r["step"] for r in per if r["pf_A"] != r["pf_B"]]
    dmax, dmin = max(D), min(D)
    why = f"max_k D_k {dmax:.4g} mm・min_k D_k {dmin:.4g} mm・σ_ref {sigma:.4g} mm (許容 {tol} mm)"
    if mis:
        return "結果 B", f"抽出の合否が変わる出力 {mis[:6]} — 旧い参照の流用を退ける ({why})"
    if dmax + sigma <= tol:
        return "結果 A", f"共通の使用域の全断面・判定窓で max|δ_E,B − δ_E,A| ≤ {tol} mm かつ合否が同じ ({why}) — 「旧い参照が共通区間で予算を超えるずれを生む」を退ける (旧い参照の出口の不足は残る)"
    if dmin - sigma > tol:
        return "結果 B", f"全出力で max|δ_E,B − δ_E,A| > {tol} mm ({why}) — 旧い参照の流用を退ける"
    return "判定不能", f"窓の中の時間の変動か参照の未整定で {tol} mm をまたぐ ({why})"


def ab_full(run: Path, s: int, specs: dict, ext: dict, problem: Path, x_w: float, mm: float) -> dict:
    """停止の出力 s の両腕の全体の抽出 (run/H370_AB/<腕>) と、帯の位置・感度・held・δ_target・次の壁の比較 (記録)。"""
    from forge_design.feedback import deltastar_loop as DL
    ext_fn = HOOKS.get("extract_and_merge", DL.extract_and_merge)
    nw_fn = HOOKS.get("next_wall_check", DL.next_wall_check)
    out, tabs = {}, {}
    for arm, spec in specs.items():
        od = run / AB_DIR / arm
        if od.exists():
            shutil.rmtree(od)
        od.mkdir(parents=True)
        try:
            with snapshot_view(run, s) as v:
                summ = ext_fn(v, spec, out_dir=od, **ext)
        except Exception as e:  # noqa: BLE001
            out[arm] = {"error": f"{type(e).__name__}: {e}"[:600]}
            continue
        ne = summ.get("not_extracted") or {}
        out[arm] = {"n_stations": summ.get("n_stations"), "n_held": summ.get("n_held"), "lam": (summ.get("smooth") or {}).get("lam"),
                    "not_extracted": {"n": ne.get("n"), "x_range": [min(ne["x"]), max(ne["x"])] if ne.get("x") else None,
                                      "reasons": sorted(set(ne.get("reason") or []))}}
        tabs[arm] = {"next": np.loadtxt(od / "delta_r_next.csv", delimiter=",", skiprows=1, ndmin=2),
                     "equiv": np.loadtxt(od / "delta_r_equiv.csv", delimiter=",", skiprows=1, ndmin=2),
                     "edge": (np.loadtxt(od / "delta_r_edge.csv", delimiter=",", skiprows=1, ndmin=2) if (od / "delta_r_edge.csv").is_file() else None),
                     "head": (od / "delta_r_next.csv").read_text().splitlines()[0]}
    if set(tabs) != {"A", "B"}:
        out["compare"] = {"error": "どちらかの腕の抽出が失敗"}
        return out
    A, B = tabs["A"], tabs["B"]
    xa, xb = A["next"][:, 0], B["next"][:, 0]
    use = lambda x: (x >= 0.0) & (x <= x_w + 1e-9)  # noqa: E731
    ca = np.isin(xa, xb) & use(xa)
    cb = np.isin(xb, xa) & use(xb)
    cmp = {"n_common_use": int(ca.sum()), "x_common": [float(xa[ca].min()), float(xa[ca].max())] if ca.any() else None,
           "B_use_not_in_A": {"n": int((use(xb) & ~np.isin(xb, xa)).sum()),
                              "x_range": ([float(xb[use(xb) & ~np.isin(xb, xa)].min()), float(xb[use(xb) & ~np.isin(xb, xa)].max())]
                                          if (use(xb) & ~np.isin(xb, xa)).any() else None)}}
    if ca.any():
        ia, ib = np.flatnonzero(ca), np.flatnonzero(cb)
        for name, col in (("delta_E", 3), ("delta_target", 1)):
            d = np.abs(A["next"][ia, col] - B["next"][ib, col]) * mm
            cmp[f"{name}_max_diff_mm"] = float(d.max())
            cmp[f"{name}_x_at_max"] = float(xa[ia][int(np.argmax(d))])
        cmp["held_common"] = {"A": int(A["next"][ia, 4].sum()), "B": int(B["next"][ib, 4].sum())}
        ea, eb = A["equiv"], B["equiv"]
        ja, jb = np.flatnonzero(np.isin(ea[:, 0], xa[ia])), np.flatnonzero(np.isin(eb[:, 0], xb[ib]))
        if len(ja) == len(jb) and len(ja):
            sa = np.max(np.abs(ea[ja][:, 7:9] - ea[ja][:, [1]]), axis=1) / np.maximum(np.abs(ea[ja][:, 1]), 1e-30)
            sb = np.max(np.abs(eb[jb][:, 7:9] - eb[jb][:, [1]]), axis=1) / np.maximum(np.abs(eb[jb][:, 1]), 1e-30)
            cmp["band_sensitivity_rel"] = {"A_max": float(np.nanmax(sa)), "B_max": float(np.nanmax(sb)), "max_abs_diff": float(np.nanmax(np.abs(sa - sb)))}
        if A["edge"] is not None and B["edge"] is not None:
            ga, gb = A["edge"], B["edge"]
            ka, kb = np.flatnonzero(np.isin(ga[:, 0], xa[ia])), np.flatnonzero(np.isin(gb[:, 0], xb[ib]))
            if len(ka) == len(kb) and len(ka):
                for name, col in (("y_b", 5), ("y_star", 2)):
                    d = np.abs(ga[ka, col] - gb[kb, col])
                    fin = np.isfinite(d)
                    cmp[f"band_{name}_max_diff_mm"] = float(np.max(d[fin]) * mm) if fin.any() else None
                    cmp[f"band_{name}_max_rel"] = (float(np.max(d[fin] / np.maximum(np.abs(gb[kb, col][fin]), 1e-30))) if fin.any() else None)
                    cmp[f"band_{name}_n_nonfinite"] = int((~fin).sum())
        # 未緩和の固定点の残差 (設計を終える条件の量、記録だけ): 腕 B の使用域の max|δ_E − δ_in|
        ub = use(xb)
        r = np.abs(B["next"][ub, 3] - B["next"][ub, 2]) * mm
        cmp["B_fixed_point_residual_mm"] = {"max": float(r.max()), "x_at_max": float(xb[ub][int(np.argmax(r))]), "note": "0.05 mm は設計を終える条件 (ここでは記録だけ)"}
    out["compare"] = cmp
    # 次の壁: 腕 B の表 / 腕 A の共通の列の δ_target だけを替えた表 (出口の区間は腕 B のまま = 延ばさない)
    nw = {}
    try:
        nw["B"] = nw_fn(problem, run / AB_DIR / "B" / "delta_r_next.csv")
        hyb = B["next"].copy()
        if ca.any():
            hyb[np.flatnonzero(cb), 1] = A["next"][np.flatnonzero(ca), 1]
        hd = run / AB_DIR / "A_on_common"
        hd.mkdir(parents=True, exist_ok=True)
        np.savetxt(hd / "delta_r_next.csv", hyb, delimiter=",", comments="", header=B["head"] + " [H370: 共通の列の delta_target を腕 A に替えた表]",
                   fmt="%.18e")
        nw["A_on_common"] = nw_fn(problem, hd / "delta_r_next.csv")
        ta, tb = (nw["A_on_common"] or {}).get("throat"), (nw["B"] or {}).get("throat")
        if ta and tb:
            nw["throat_diff"] = {"x_rt": ta["x"] - tb["x"], "r_mm": (ta["r"] - tb["r"]) * mm}
    except Exception as e:  # noqa: BLE001
        nw["error"] = f"{type(e).__name__}: {e}"[:600]
    out["next_wall"] = nw
    return out


def ab_compare(case: Path, run: Path, s: int, rec: dict, spec_b: dict, qual: dict, problem: Path, x_w: float) -> dict:
    """抽出の参照の A/B (上の docstring の 5)。"""
    from forge_design.evaluate import stop_rule as SR
    from forge_design.feedback import deltastar_loop as DL
    snap = HOOKS.get("delta_snapshot", SR.delta_snapshot)
    st, ds = rec.get("settings") or {}, rec.get("delta_series") or {}
    W = int(st.get("window", 0))
    x_use = np.asarray(ds.get("x") or [], dtype=float)
    S = ds.get("scale_m")
    out = {"tol_mm": TOL_MM, "arm_B": spec_b, "windows": [[s - 2 * W, s - W], [s - W, s]], "result": "判定不能", "why": None}
    ext = {k: v for k, v in (ds.get("extract") or DL.handoff_extract_settings()).items() if k != "source"}
    out["extract"] = ext
    old = case / OLD_EULER
    ors = res_steps(old) if old.is_dir() else []
    if not ors or not len(x_use) or S is None or not W:
        out["why"] = f"腕 A の {OLD_EULER} の出力が無いか、stop_rule.json に窓・δ の列が無い"
        return out
    mm = float(S) * 1e3
    spec_a = {"run": str(old.resolve()), "res": f"res_{ors[-1]}.h5", "sha256": sha256_file(old / f"res_{ors[-1]}.h5")}
    out["arm_A"] = spec_a
    rows = {r["step"]: r for r in rec.get("rows") or [] if isinstance(r.get("step"), int)}
    steps = [k for k in sorted(rows) if s - 2 * W <= k <= s]
    per = []
    for k in steps:
        f = run / f"res_{k}.h5"
        if not f.is_file():
            per.append({"step": k, "error_missing": "出力が無い"})
            continue
        sha_k = sha256_file(f)
        dB, src = rows[k].get("delta"), "記録 (stop_rule.json の delta)"
        if not (isinstance(dB, dict) and "delta_E" in dB and len(dB["delta_E"]) == len(x_use) and dB.get("res_sha256") == sha_k):
            src = "抽出し直し (記録が無いか出力の sha が違う)"
            try:
                dB = snap(run, k, spec_b, ext, x_use)
            except Exception as e:  # noqa: BLE001
                dB = {"error": f"{type(e).__name__}: {e}"[:300]}
        try:
            dA = snap(run, k, spec_a, ext, x_use)
        except Exception as e:  # noqa: BLE001
            dA = {"error": f"{type(e).__name__}: {e}"[:300]}
        if k == s:
            out["_dB_s"] = dB
        per.append(ab_pair(k, dA, dB, mm, src, x_use))
    out["per_output"] = per
    # 参照の時間の変動 σ_ref (腕 B の Euler の別の出力で停止の出力 s を抽出し直す)
    sig, sig_rec = None, []
    dbs = out.pop("_dB_s", None)
    snaps = (qual or {}).get("sensitivity_snapshots") or []
    if dbs is not None and "delta_E" in dbs and snaps:
        try:
            e0 = _arr(dbs, "delta_E")
            for sn in snaps:
                d = snap(run, s, {"run": spec_b["run"], "res": sn["res"], "sha256": sn["sha256"]}, ext, x_use)
                e1 = _arr(d, "delta_E")
                m = np.isfinite(e0) & np.isfinite(e1)
                sig_rec.append({"res": sn["res"], "max_mm": float(np.max(np.abs(e1 - e0)[m]) * mm) if m.any() else None})
            vals = [v["max_mm"] for v in sig_rec]
            sig = None if any(v is None for v in vals) else float(max(vals))
        except Exception as e:  # noqa: BLE001
            sig_rec.append({"error": f"{type(e).__name__}: {e}"[:300]})
    out["sigma_ref"] = {"mm": sig, "per_snapshot": sig_rec,
                        "note": "停止の出力 s を、固定した Euler の出力の代わりに資格の記録の別の出力 (check_quasisteady の末尾の窓の始まり・中ほど) で抽出した δ_E の差"}
    ok = [r for r in per if "D_mm" in r]
    if ok:
        out["D_mm"] = {"max": max(r["D_mm"] for r in ok), "min": min(r["D_mm"] for r in ok), "n_outputs": len(per),
                       "n_common": sorted({r["n_common"] for r in ok}), "x_common": ok[-1].get("x_common")}
        i_max = max(ok, key=lambda r: r["D_mm"])
        out["D_mm"]["x_at_max"] = float(x_use[i_max["i_at_D"]])
        out["D_mm"]["step_at_max"] = i_max["step"]
        out["excluded_from_A"] = {"n": sorted({r["n_excluded_A"] for r in ok}), "x_range": ok[-1].get("x_excluded_A"),
                                  "note": "腕 A (旧い参照) が抽出しない列 (範囲の外 = 出口の区間・端の欠測) は比べない (外挿しない)"}
    out["result"], out["why"] = classify_ab(per, sig)
    try:
        out["at_stop"] = ab_full(run, s, {"A": spec_a, "B": spec_b}, ext, problem, x_w, mm)
    except Exception as e:  # noqa: BLE001 — 記録だけ
        out["at_stop"] = {"error": f"{type(e).__name__}: {e}"[:600]}
    return out


def gate_summary(g: dict) -> dict:
    c = g.get("checks") or {}
    pick = lambda d, ks: {k: d.get(k) for k in ks} if isinstance(d, dict) else d  # noqa: E731
    return {"proceed": g.get("proceed"), "blocked_by": g.get("blocked_by"), "stop_rule": g.get("stop_rule"),
            "convergence": pick(g.get("convergence") or {}, ("kind", "rc", "verdict", "cmd")),
            "conservation": pick(c.get("conservation"), ("pass", "r_mean", "r_amp", "n", "why")),
            "wall_resolution": pick(c.get("wall_resolution"), ("pass", "kind", "verdict", "cmd")),
            "two_windows": pick(c.get("two_windows"), ("pass", "windows", "n", "worst", "why")),
            "euler_ref": c.get("euler_ref"), "omega": c.get("omega"),
            "coverage": pick(c.get("coverage"), ("pass", "n_use", "n_held_use", "n_nonfinite_use", "not_extracted", "why")),
            "next_wall": pick(c.get("next_wall"), ("pass", "throat", "upstream_gate_pass", "validate", "why")),
            "made": (g.get("extract") or {}).get("made"), "override": g.get("override")}


def ns_judge(case: Path, tools: Path, do_ab: bool, out_p: Path) -> int:
    from forge_design.evaluate import stop_rule as SR
    from forge_design.feedback import deltastar_loop as DL
    run = case / NS_RUN
    R = {"plan": PLAN, "created": now(), "run": NS_RUN, "case_dir": str(case), "tools": str(tools), "undecidable": [], "records": {},
         "note": "壁は更新しない (固定壁の pass)。0.05 mm の固定点の条件は設計を終える用で、ここでは記録だけ"}
    und = R["undecidable"]

    def finish(rc):
        jwrite(out_p, R)
        print(f"== VERDICT V-h370: {R.get('verdict')}")
        print(f"→ {out_p}")
        return rc
    if not run.is_dir():
        und.append(f"{NS_RUN} が無い")
        R["verdict"] = "判定不能: " + "; ".join(und)
        return finish(2)
    pin, qual = jload(case / PIN), jload(case / QUAL)
    spec = {"run": str(case / pin["run"]), "res": pin["res"], "sha256": pin["sha256"]} if pin else None
    R["euler_ref"] = spec
    if not pin or not (qual or {}).get("qualified"):
        und.append("Euler の参照の固定・資格の記録が無い")
    # --- 1. 記録 ---
    pc, ms = jload(run / "H370_PRECHECK.json") or {}, jload(run / "H370_MAIN_START.json") or {}
    R["records"]["precheck_ok"], R["records"]["main_start_ok"] = pc.get("ok"), ms.get("ok")
    if not pc.get("ok"):
        und.append("投入前の検査 (H370_PRECHECK.json) が無いか不成立")
    if not ms.get("ok"):
        und.append("本段の直前の照合 (H370_MAIN_START.json) が無いか不成立")
    for m in MARKERS:
        if (run / m).is_file():
            und.append(f"基盤の停止・不成立の印 {m}: {(run / m).read_text()[:200]}")
    R["records"]["run_rc"] = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").is_file() else None
    prov = (run / "RUN_PROVENANCE.txt").read_text() if (run / "RUN_PROVENANCE.txt").is_file() else ""
    fsha = next((ln.split(":", 1)[1].strip() for ln in prov.splitlines() if ln.startswith("forge_sha256")), None)
    envs = [ln.split(":", 1)[1].strip() for ln in prov.splitlines() if ln.startswith("env")]
    R["records"]["provenance"] = {"forge_sha256": fsha, "env": envs}
    if not (fsha or "").startswith(BIN_SHA) or "FORGE_CUDA_BLOCKSIZE=128" not in envs:
        und.append(f"RUN_PROVENANCE.txt のバイナリ ({str(fsha)[:16]}) か FORGE_CUDA_BLOCKSIZE が登録と違う")
    R["records"]["residual_nonfinite"] = {t: residual_nonfinite(run / f"residual_history_{t}.csv") for t in ("S1_soft", "S2_mid", "main")}
    # --- 2. 停止の状態 ---
    state = SR.read_state(run) or {}
    rec = jload(run / SR.RECORD) or {}
    status = state.get("status") or "MISSING"
    s = state.get("stop_step") if state.get("stop_step") is not None else rec.get("stop_step")
    rows = [r for r in rec.get("rows") or [] if isinstance(r.get("step"), int)]
    last = rows[-1] if rows else {}
    R["stop"] = {"status": status, "stop_step": s, "last_output": state.get("last_output"), "forge_rc": state.get("forge_rc"),
                 "hold": rec.get("hold"), "note": state.get("note"), "sections_requested": rec.get("sections"),
                 "sections_chosen_x": {k: v for k, v in last.items() if k.startswith("x_theta_r_")}, "recorded_at_stop": rec.get("recorded_at_stop")}
    if status == "REACHED" and s is not None:
        f = run / f"res_{s}.h5"
        R["stop"]["stop_sha256_match"] = f.is_file() and sha256_file(f) == state.get("stop_sha256")
        if not R["stop"]["stop_sha256_match"]:
            und.append(f"停止の出力 res_{s}.h5 が無いか sha が STOP_RULE_STATE.json と違う")
        # 壁の出力の名前は投入前の検査の記録 (wall_output_check) から、無ければ case/45 の bcond の壁 (wall・physID 3)
        wname = (((pc.get("items") or {}).get("N3_wall_output") or {}).get("wall_dump_name") or "res_wall_3_<step>.h5").replace("<step>", str(s))
        R["stop"]["wall_dump_at_stop"] = sorted(p.name for p in run.glob(f"res_*_*_{s}.h5"))
        R["stop"]["wall_dump_required"] = wname
        if not (run / wname).is_file():
            und.append(f"停止の出力 s = {s} の壁の出力 {wname} が無い (局所 y₁⁺ を判定できない)")
    # --- 3. 受け渡しの関門 ---
    gate = None
    if spec is not None:
        try:
            gate = HOOKS.get("handoff_gate", DL.handoff_gate)(run, problem=case / NS_PROBLEM, euler_ref=spec, omega=DL.HANDOFF_OMEGA,
                                                             allow_not_reached=False, allow_unconverged=False, write=True)
            R["gate"] = gate_summary(gate)
        except Exception as e:  # noqa: BLE001
            und.append(f"受け渡しの関門が例外 ({type(e).__name__}: {e})"[:400])
    # --- 4. check_quasisteady (記録) ---
    W = int((rec.get("settings") or {}).get("window", 0) or 0)
    if s is not None and W and rows:
        from forge_design.evaluate.stop_rule import theta_key
        xs_th = [theta_key(f) for f in (rec.get("settings") or {}).get("x_fracs") or []]
        x_use = np.asarray((rec.get("delta_series") or {}).get("x") or [], dtype=float)
        ser = []
        for r in rows:
            if not (s - 2 * W <= r["step"] <= s) or r.get("nonfinite_all"):
                continue
            o = {"step": r["step"], "r": r.get("r", float("nan")), **{k: r.get(k, float("nan")) for k in xs_th},
                 "exit_core_M": r.get("exit_core_M", float("nan")), "mdot_rt": r.get("mdot_rt", float("nan"))}
            d = r.get("delta")
            for x in XS:
                ok = isinstance(d, dict) and "delta_E" in d and len(d["delta_E"]) == len(x_use) and len(x_use)
                o[f"dE_{x:g}"] = float(np.interp(x, x_use, d["delta_E"])) if ok else float("nan")
                o[f"dT_{x:g}"] = float(np.interp(x, x_use, d["delta_target"])) if ok else float("nan")
            ser.append(o)
        cols = ["r"] + xs_th + ["exit_core_M", "mdot_rt"] + [f"d{q}_{x:g}" for q in ("E", "T") for x in XS]
        if ser:
            sc = run / "h370_qs_series.csv"
            write_series(sc, ser, cols)
            R["quasisteady"] = {"window": [s - 2 * W, s], "n": len(ser), "default": run_qs(tools, sc, cols),
                                "tight (1e-4)": run_qs(tools, sc, cols, QS_TIGHT),
                                "note": "記録だけ (§6 V-h370: 既定の相対閾値で寸法の予算の合格を代用しない; 判定は関門の 2 窓)"}
    # --- 5. 抽出の参照の A/B ---
    if do_ab and spec is not None and s is not None and status in ("REACHED", "CENSORED") and rows:
        info = jload(run / "prepare_info.json") or {}
        x_w = float(((info.get("physical_wall") or {}).get("domain") or [0.0, float("inf")])[1])
        try:
            R["ab"] = ab_compare(case, run, int(s), rec, spec, qual, case / NS_PROBLEM, x_w)
            if status != "REACHED":
                R["ab"]["note"] = f"停止の状態 {status} の窓で行った診断の抽出 (更新の許可とは別)"
        except Exception as e:  # noqa: BLE001
            R["ab"] = {"result": "判定不能", "why": f"A/B が例外 ({type(e).__name__}: {e})"[:400]}
    else:
        R["ab"] = {"result": "判定不能" if do_ab else "未実施 (--no-ab)",
                   "why": f"停止の状態 {status}・停止の出力 {s} では A/B の窓を作れない" if do_ab else None}
    # --- 判定 ---
    proceed = bool(gate and gate.get("proceed"))
    upd = ("更新に進める条件をすべて満たす (壁の更新は自動でつながない)" if proceed
           else "更新に進める条件を満たさない: " + "; ".join((gate or {}).get("blocked_by") or ["関門の記録が無い"])[:900])
    R["update_conditions_met"] = proceed and not und
    abr = R["ab"].get("result")
    if und:
        R["verdict"] = "判定不能: " + "; ".join(und)[:900] + f" (記録: 関門 {('進める' if proceed else '止める') if gate else 'なし'}、A/B {abr})"
        return finish(2)
    R["verdict"] = f"{upd} / 抽出の参照の A/B: {abr}" + (f" ({R['ab'].get('why')})" if R["ab"].get("why") else "")
    return finish(0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--euler-qualify", action="store_true", help="(i) Euler の参照の資格の判定と固定")
    ap.add_argument("--no-ab", action="store_true", help="抽出の参照の A/B を回さない (判定の A/B は未実施になる)")
    ap.add_argument("--case-dir", default=str(HERE))
    ap.add_argument("--design", default=None, help="forge_design の design/ (既定: ~/forge-g0cchain の、無ければこの作業ツリーの)")
    ap.add_argument("--tools", default=None, help="forge の tools (既定: ~/forge-g0cchain の、無ければこの作業ツリーの)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    case = Path(a.case_dir).resolve()
    design = Path(a.design) if a.design else Path.home() / "forge-g0cchain/design"
    if not design.is_dir():
        design = HERE.parents[1] / "design"
    tools = Path(a.tools) if a.tools else Path.home() / "forge-g0cchain/solver_density_cuda/tools"
    if not tools.is_dir():
        tools = HERE.parents[1] / "solver_density_cuda/tools"
    if str(design) not in sys.path:
        sys.path.insert(0, str(design))
    if a.euler_qualify:
        return euler_qualify(case, tools)
    return ns_judge(case, tools, not a.no_ab, Path(a.out) if a.out else case / OUT)


if __name__ == "__main__":
    sys.exit(main())
