"""plan tooling-nozzle-isothermal-wall-chain §5.1 #36f・#36g・§6「V-h370 の事前登録」(§4.9 codex diagnose 3、判定の基準の諮問
`notes/reviews/2026-10-14-h370-judge-criteria-diagnose.md`) の判定。B の case dir で回す (h370_run.py が呼ぶ)。壁は更新しない (固定壁の pass)。

--euler-prequalify: (i) Euler 参照 run_0950_h370_euler_ref の**予備の資格** (prequalify; §6 V-h370 の「Euler」、§5.1 #36g (1)):
  - RUN_RC 0、投入前の検査 (H370_PRECHECK.json) の合格、本段が上限 60,000 step まで回った (最後の出力 res_60000.h5)。
  - NaN/Inf なし: 全段の残差の履歴 (residual_history_{soft,mid,main}.csv) の rms_* と、本段の全出力の VALUE の浮動小数点の量。ρ・P・T > 0。
  - 全残差: check_convergence.py --segment --drop 3 (stage_manifest の最後の区間 = 方程式・BC・離散化が同じ本段の区間、全列) の PASS。
  - 出口 M・流量・参照場の時系列 (本段の全出力) の check_quasisteady.py の VERDICT (CSV 系列モード、閾値は既定のまま): 量は
    出口コア M・無次元の流量 (stop_rule.recorded_quantities = nozzle_report.metrics・qs_series.py と同じ定義) と、参照場として軸 (η = 0) と
    η = 0.1 の線の M を x = 40・70・94 r_t で (nozzle_report.eta_line、停止規則の断面と同じ位置)。すべて STEADY であること
    (OSCILLATING も参照の時間の変動なので未達)。nozzle_report.metrics の η = 0.1 の線の評価量 (wave・overshoot・slope・range、%) は
    平均が 0 に近く相対の閾値が意味を持たないので記録だけ。閾値を細かくした (1e-4) 判定も記録だけ。
  これは 6 点の M・出口 M・流量の既定の相対閾値の準定常で、抽出に使う参照場全体の資格には足りない (codex の Major)。**最終の資格**は
  NS の判定の σ_ref (下の 5) で決め、σ_ref を測れない・判定不能・参照の時間の変動が A/B の判定を左右する (参照の未整定) なら
  `update_conditions_met` を偽にする。
  満たせば本段の最後の出力とその sha256 を _band_ab/half_inlet/h370_euler_ref_pin.json に固定する (既にあって違えば止める; 消さない)。
  σ_ref の代わりの参照 (Euler の末尾の窓 = 本段の最後の 20,000 step の出力から 5,000 step ごとの代表 40,000・45,000・50,000・55,000、
  固定の 60,000 と合わせて 5 出力; 選び方は NS の判定の前のここで固定する) とその sha、末尾の窓の全出力の sha と場の差 (記録だけ;
  代表の選び方の裏付け)、スクリーニングの出力 (check_quasisteady の末尾の窓の始まりと中ほど、停止の出力だけに当てる) を記録する。
  → _band_ab/half_inlet/h370_euler_prequalify.json。終了コード 0 予備の資格あり / 1 未達 / 2 判定不能。

(既定): (ii) NS run_0951_h370_ns_b0 の §6 V-h370 の判定 → _band_ab/half_inlet/h370_judge.json:
  1. 記録: 投入前の検査 (H370_PRECHECK.json)・本段の直前の照合 (H370_MAIN_START.json)・基盤の停止 (PRECHECK_FAILED・WATCH_ABORT・
     RESOURCE_STOP・DISK_STOP・WRITE_FAILED)・RUN_PROVENANCE.txt のバイナリと環境変数・段の残差の履歴の非有限。
  2. 停止の状態 (STOP_RULE_STATE.json): 状態・停止の出力 s・到達と中断と保持 (stop_rule.json の hold)・要求の断面と実際の列・停止の出力の sha。
  3. 受け渡しの関門 (deltastar_loop.handoff_gate、HANDOFF_GATE.json): 停止の状態 REACHED・check_convergence.py --segment --drop 3 の PASS・
     保存則・局所 y₁⁺ (check_wall_resolution.py --step s)・2 窓 [s − 40000, s − 20000]・[s − 20000, s] の δ_E・δ_target (全使用域で
     最大 − 最小 ≤ 0.01 mm、窓の平均の差 ≤ 0.01 mm; 入るべき 17 出力を出力の間隔から列挙)・Euler の参照と ω・抽出の全域の被覆・
     使用域の held と欠測なし・次の壁の幾何。上書き (allow_not_reached・allow_unconverged) は使わない。合格なら関門が抽出の表を run の直下に
     写す (次の pass は自動でつながない)。
  4. check_quasisteady.py の VERDICT (記録; 既定の相対閾値で寸法の予算の合格を代用しない): 2 窓の r・θ_r (断面ごと)・出口コア M・流量・
     δ_E・δ_target (x = 40・70・94 r_t) の時系列。
  5. 抽出の参照の A/B と σ_ref (同じ NS の保存系列に対して Euler の参照だけを変える; §5.1 #36g (2)〜(5)):
     - 判定窓の完全性: 停止の時刻 s・出力の間隔 (stop_rule.json の settings.out_step_interval、2,500)・2 窓から入るべき出力
       (stop_rule.required_window_steps、17 出力) を列挙し、stop_rule.json の行と res の有無と照合する。行が無くても res があれば
       抽出し直して補い、res が無ければ判定不能。
     - 比べる列 F を**抽出の前に**座標だけで決める: NS の使用域の列 (stop_rule.json の delta_series.x) のうち、腕 A (旧い
       run_0174_euler_v5d_M_r1、最後の出力を sha で固定) と腕 B (固定した run_0950) の Euler の補間できる範囲 (共通の関数
       deltastar.euler_cover_range、受け渡しの抽出と同じ端の幾何の許容つき) の両方に入る列。F の中の非有限・欠測・held は比較を
       合格にしない (後から有限の列だけを選ばない)。旧い参照の範囲の外 (出口の区間) は F に入らないので比べない (外挿しない)。
     - 出力 k ごと: 各腕の合否 pf (抽出の例外なし・F の中に held・非有限・欠測なし)、両腕とも pf のとき
       D_k = F の max|δ_E,B − δ_E,A| [mm]。抽出の条件は受け渡しと同じ (edge・一定 x の断面への補間・端は欠測・ω 0.5・平滑化の λ・
       昇格の上限; stop_rule.json の extract)。
     - σ_ref (腕 B の参照の時間の変動): 必要な全出力 k について、固定した出力の代わりに予備の資格で固定した代わりの参照 (Euler の末尾の
       窓の代表 4 出力) で抽出し直し、σ_ref = max_k max_e max_F |δ_E(k, e) − δ_E(k, 固定)| [mm] (使用域の全列の値 σ_ref,use も記録)。
       代わりの参照で F・使用域の held・欠測・合否が変われば σ_ref は判定不能、抽出が失敗・代わりの参照が無い・sha が違えば測れない。
       停止の出力 × スクリーニングの 2 出力の差は「スクリーニング」として記録だけ (σ_ref の上限とは扱わない)。
     - 分類 (classify_ab): 必要な出力の欠け・抽出の失敗・σ_ref を測れない・判定不能 → 判定不能。どれかの出力で pf_A ≠ pf_B → 代わりの
       参照の全部でもその違いが保たれれば 結果 B (保たれなければ判定不能)。両腕とも pf でない出力 (F の中の欠測) → 判定不能。
       それ以外: max_k D_k + σ_ref ≤ 0.01 mm → 結果 A、min_k D_k − σ_ref > 0.01 mm → 結果 B、σ_ref を足し引きして初めて 0.01 mm を
       またぐ → 判定不能 (参照の未整定)、D_k だけで窓の中で 0.01 mm をまたぐ → 判定不能 (時間の変動)。
     - 参照の最終の資格 (reference_final): σ_ref を測れる・判定できる・参照の未整定でない (σ_ref が A/B の分類を左右しない) とき資格あり。
       A/B が分類に至らない (出力の欠け・抽出の失敗・F の中の欠測) ときは σ_ref が判定を左右するかを決められないので資格を確かめられない。
     停止の出力 s では両腕の全体の抽出 (run の H370_AB/{A,B}) で帯の位置 (y_b・y*)・感度・held・δ_target の差・次の壁 (腕 B の表、腕 A は
     共通の列だけ腕 A の δ_target に替えた表 = 出口の区間を延ばさない) の真のスロートを記録する。
     抽出の回数の目安 (V-h370、1 抽出 ~116 s = #36e の実寸の格子の実測): 腕 A 17 + σ_ref 17 × 4 + スクリーニング 2 + 停止の出力の全体 2
     ≈ 89 回 ≈ 2.9 h (腕 B は記録を使う)。末尾の窓の全 20 出力を代わりにすると ≈ 361 回 ≈ 11.6 h。
  6. 設計を終える条件の量 (記録だけ; 最初の pass の判定には使わない): 使用域の未緩和の max|δ_E − δ_in|。
  判定: 「更新に進める条件」= 関門の proceed かつ 参照の最終の資格 (と 1 の基盤の不成立が無いこと)、抽出の参照の A/B は
  結果 A / 結果 B / 判定不能。
usage: python3 h370_judge.py --euler-prequalify     python3 h370_judge.py [--no-ab]
       (--no-ab は A/B と σ_ref を回さない = 参照の最終の資格を確かめないので、更新に進める条件は満たさない)
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
PLAN = "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36f・#36g・§6 V-h370 (§4.9 codex diagnose 3、判定の基準の諮問)"
EU_RUN, NS_RUN = "run_0950_h370_euler_ref", "run_0951_h370_ns_b0"
NS_PROBLEM = "problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml"
OLD_EULER = "run_0174_euler_v5d_M_r1"
PIN = "_band_ab/half_inlet/h370_euler_ref_pin.json"
QUAL = "_band_ab/half_inlet/h370_euler_prequalify.json"   # 予備の資格の記録 (§5.1 #36g (1); 旧名 h370_euler_qualify.json)
OUT = "_band_ab/half_inlet/h370_judge.json"
EU_STEPS = 60000
EU_TAIL = 20000                               # σ_ref の Euler の末尾の窓 (本段の最後の 20,000 step、出力 1,000 ごと 21 出力)
SIGMA_REF_EVERY = 5000                        # 末尾の窓の中の代わりの参照の代表の間隔 (40,000・45,000・50,000・55,000 + 固定の 60,000)
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


def sigma_ref_selection(steps: list, last: int = EU_STEPS, tail: int = EU_TAIL, every: int = SIGMA_REF_EVERY) -> dict:
    """σ_ref の代わりの参照の選び方 (§5.1 #36g (2); 予備の資格の時点 = NS の判定の前に固定する)。末尾の窓 [last − tail, last] の出力のうち、
    窓の始まりから every step ごとの出力 (固定の last を除く) を代表にする (V-h370: 40,000・45,000・50,000・55,000)。
    戻り {"tail_window", "tail_steps" (窓の全出力), "representatives", "missing" (代表のうち出力が無いもの), "rule", "why"}。"""
    lo = last - tail
    tail_steps = [s_ for s_ in steps if lo <= s_ <= last]
    want = [lo + k * every for k in range(tail // every + 1) if lo + k * every < last]
    return {"tail_window": [lo, last], "tail_steps": tail_steps, "representatives": [s_ for s_ in want if s_ in tail_steps],
            "missing": [s_ for s_ in want if s_ not in tail_steps],
            "rule": f"末尾の窓 [{lo}, {last}] の {every} step ごとの出力 (固定の {last} を除く) を代わりの参照にする",
            "why": ("必要な NS の 17 出力 × 末尾の窓の全 20 出力は 340 抽出 ≈ 11 h (1 抽出 ~116 s、#36e の実寸の格子の実測) で重いので、"
                    f"{every} step ごとの代表 {len(want)} 出力にした (17 × {len(want)} = {17 * len(want)} 抽出 ≈ {17 * len(want) * 116 / 3600:.1f} h)。"
                    "代表の間の出力だけの変動は見ない (末尾の窓の全出力の場の差 tail_field_dev を記録して裏付ける)")}


def tail_field_dev(run: Path, steps: list, pin_step: int, reps: list) -> dict:
    """記録だけ (代表の選び方の裏付け): 末尾の窓の各出力の ρU_x の固定の出力との差 max|ρU_x − ρU_x,pin| / max|ρU_x,pin|。代表の中の最大と
    代表の外の最大を並べる (代表の外の出力だけが大きく違えば、代表で測る σ_ref はそれを見ない)。"""
    import h5py

    def q(st):
        with h5py.File(run / f"res_{st}.h5", "r") as h:
            return np.asarray(h["/VALUE/ro"][:], dtype=float) * np.asarray(h["/VALUE/Ux"][:], dtype=float)
    q0 = q(pin_step)
    den = float(np.max(np.abs(q0))) or 1.0
    per = [{"step": st, "max_rel": float(np.max(np.abs(q(st) - q0)) / den)} for st in steps if st != pin_step]
    rep = [p_["max_rel"] for p_ in per if p_["step"] in reps]
    oth = [p_ for p_ in per if p_["step"] not in reps]
    o = max(oth, key=lambda p_: p_["max_rel"]) if oth else None
    return {"per_output": per, "max_rel_representatives": max(rep) if rep else None,
            "max_rel_others": o["max_rel"] if o else None, "step_max_others": o["step"] if o else None,
            "note": "記録だけ (合否に使わない)。代表の外の最大が代表の中の最大より大きいときは、σ_ref がその出力の参照の変動を見ていない"}


def euler_prequalify(case: Path, tools: Path) -> int:
    """(i) Euler の参照の予備の資格と固定 (上の docstring)。終了コード 0 予備の資格あり / 1 未達 / 2 判定不能。"""
    run = case / EU_RUN
    R = {"plan": PLAN, "created": now(), "run": EU_RUN, "tools": str(tools), "undecidable": [], "failures": [], "prequalified": False,
         "kind": "予備の資格 (prequalify)。最終の資格は NS の判定の σ_ref (h370_judge.py の既定の経路の reference_final) が決め、"
                 "update_conditions_met に繋ぐ (§5.1 #36g (1))",
         "definitions": {"qualifying_cols": QUAL_COLS, "record_cols": REC_COLS, "positions_rt": list(XS),
                         "quasisteady": "check_quasisteady.py --series-csv (既定の --tail 0.4・--drift 0.05・--osc 0.10・--min-snaps 4)、全列 STEADY",
                         "exit_core_M・mdot_rt": "forge_design.evaluate.stop_rule.recorded_quantities (nozzle_report.metrics・qs_series.py と同じ定義)",
                         "M_axis・M_eta0.1": "forge_design.report.nozzle_report.eta_line(F, 'M', η, x) (η = 0 は軸、0.1 は metrics の判定の線)",
                         "record_cols": "nozzle_report.metrics の eta0.1 (wave_pct・overshoot_pct・slope_pct・range_pct)、記録だけ",
                         "pin": "本段の最後の出力 (上限 60,000 step) とその sha256",
                         "sigma_ref_snapshots": f"σ_ref の代わりの参照: 末尾の窓 (最後の {EU_TAIL} step) の {SIGMA_REF_EVERY} step ごとの代表 (固定を除く)",
                         "screening_snapshots": "スクリーニング (check_quasisteady の末尾の窓の始まりと中ほど; NS の停止の出力だけに当てる記録、σ_ref の上限とは扱わない)"}}
    und, fails = R["undecidable"], R["failures"]
    qual_p, pin_p = case / QUAL, case / PIN

    def finish(rc):
        R["rc"] = rc
        jwrite(qual_p, R)
        print(f"== Euler の予備の資格 ({EU_RUN}): {'予備の資格あり → ' + PIN if rc == 0 else ('未達' if rc == 1 else '判定不能')}"
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
    # σ_ref の代わりの参照 (末尾の窓の代表; NS の判定の前のここで固定する) と末尾の窓の全出力の sha
    sel = sigma_ref_selection(steps)
    R["sigma_ref_selection"] = {k: sel[k] for k in ("tail_window", "representatives", "missing", "rule", "why")}
    R["tail_outputs"] = [{"res": f"res_{s}.h5", "step": s, "sha256": sha256_file(run / f"res_{s}.h5")} for s in sel["tail_steps"]]
    R["sigma_ref_snapshots"] = [o for o in R["tail_outputs"] if o["step"] in sel["representatives"]]
    if sel["missing"]:
        # 本段が上限まで回っていれば出力の欠け (消された) で判定不能、上限より前に終わった run は既に未達 (その帰結)
        (und if steps[-1] == EU_STEPS else fails).append(f"σ_ref の代わりの参照の出力 {sel['missing']} が無い (σ_ref を測れない)")
    elif steps[-1] == EU_STEPS:
        try:
            R["tail_field_dev"] = tail_field_dev(run, sel["tail_steps"], EU_STEPS, sel["representatives"])
        except Exception as e:  # noqa: BLE001 — 記録だけ
            R["tail_field_dev"] = {"error": f"{type(e).__name__}: {e}"[:300]}
    # スクリーニング (check_quasisteady の末尾の窓 = 既定 tail 0.4 の始まりと中ほど; 停止の出力だけに当てる記録)
    n = len(steps)
    k = max(3, int(math.ceil(0.4 * n)))
    tail = steps[-k:]
    alt = sorted({tail[0], tail[len(tail) // 2]} - {steps[-1]})
    R["screening_snapshots"] = [{"res": f"res_{s}.h5", "step": s, "sha256": sha256_file(run / f"res_{s}.h5")} for s in alt]
    R["quasisteady_tail_window"] = [tail[0], tail[-1]]
    if und:
        return finish(2)
    if fails:
        return finish(1)
    pin = {"plan": PLAN, "created": now(), "run": EU_RUN, "run_abs": str(run.resolve()), "res": f"res_{steps[-1]}.h5", "step": steps[-1],
           "sha256": sha256_file(run / f"res_{steps[-1]}.h5"), "prequalify": QUAL,
           "note": "NS の euler_ref は {run: <case dir>/run_0950_h370_euler_ref, res, sha256} (h370_run.py read_pin)。予備の資格で固定し、"
                   "最終の資格は NS の判定の σ_ref"}
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
    R["prequalified"] = True
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


# A/B の分類の理由の符号 (classify_ab の code)。参照の最終の資格 (reference_final) は SIGMA_FAIL なら資格なし、REF_OK なら資格あり、
# それ以外 (A/B が分類に至らない) は σ_ref が判定を左右するかを決められないので確かめられない (更新に進めない)
SIGMA_FAIL = ("U_sigma_unmeasurable", "U_sigma_undecidable", "U_sigma_straddle")
REF_OK = ("A", "B_D", "B_pf", "U_time")


def fixed_columns(x_use, eu_runs: dict, edge_tol: bool) -> dict:
    """比べる列の集合 F を抽出の前に座標だけで決める (§5.1 #36g (3)): NS の使用域の列 x_use のうち、各腕の Euler の補間できる範囲
    (`deltastar.euler_cover_range`; edge_tol は受け渡しの抽出の remap_edge == "missing" と同じ) の全部に入る列。範囲は各 run の格子
    (prepare_info.json・nozzle.h5 の列の x) だけで決まり、場・抽出の結果を見ない。戻り {"mask" (bool の配列), "n", "x_range", "cover",
    "excluded" (腕ごとに範囲の外で外した列の数と x の範囲)}。"""
    from forge_design.metrics.delta_in_check import wall_columns_from_run
    from forge_design.metrics.deltastar import euler_cover_range, euler_covers
    xu = np.asarray(x_use, dtype=float)
    F = np.ones(len(xu), dtype=bool)
    out = {"cover": {}, "excluded": {}, "edge_tol": bool(edge_tol)}
    for arm, rd in eu_runs.items():
        rng = euler_cover_range(wall_columns_from_run(rd)["x"], edge_tol=edge_tol)
        c = euler_covers(xu, rng)
        out["cover"][arm] = rng
        out["excluded"][arm] = {"n": int((~c).sum()), "x_range": [float(xu[~c].min()), float(xu[~c].max())] if (~c).any() else None}
        F &= c
    out.update(mask=F, n=int(F.sum()), n_use=int(len(xu)), x_range=[float(xu[F].min()), float(xu[F].max())] if F.any() else None)
    return out


def arm_check(d, F) -> dict:
    """1 腕 1 出力の抽出 d (`stop_rule.delta_snapshot` の戻り) の、列の集合 F の中の合否 pf: 抽出の例外なし・F の中に held なし・δ_E・δ_target が
    有限 (抽出しなかった列 = 欠測は NaN なのでここで落ちる)。例外・形の違いは {"error"} (合否ではなく判定不能の理由)。"""
    if not isinstance(d, dict) or "error" in d or "delta_E" not in d:
        return {"error": (d or {}).get("error", "δ_E が無い") if isinstance(d, dict) else "記録が無い"}
    e, t = _arr(d, "delta_E"), _arr(d, "delta_target")
    h = np.asarray(d.get("held", []), dtype=bool)
    if not (len(e) == len(t) == len(h) == len(F)):
        return {"error": f"列の数 {len(e)}・{len(t)}・{len(h)} が比べる列の表 {len(F)} と違う"}
    miss = np.zeros(len(F), dtype=bool)
    miss[[j for j in (d.get("missing") or []) if 0 <= int(j) < len(F)]] = True
    nonfin = F & ~(np.isfinite(e) & np.isfinite(t))
    held = F & h
    return {"pf": bool(not nonfin.any() and not held.any()), "n_nonfinite": int(nonfin.sum()), "n_held": int(held.sum()),
            "n_missing": int((F & miss).sum()), "i_bad": [int(j) for j in np.flatnonzero(nonfin | held)[:8]]}


def ab_pair(step: int, dA: dict, dB: dict, mm: float, src_b: str, x_use, F) -> dict:
    """出力 1 つの腕 A・B の比較 (§5.1 #36g (3)): 比べる列は先に決めた F (抽出の後に有限の列だけを選ばない)。各腕の合否 pf_A・pf_B は F の中で
    見る。両腕とも pf のときだけ D_mm = F の max|δ_E,B − δ_E,A| [mm] (Dt_mm は δ_target)。どちらかの抽出の例外は error_<腕>。"""
    r = {"step": step, "B_source": src_b, "n_fixed": int(np.sum(F))}
    xu = np.asarray(x_use, dtype=float)
    for arm, d in (("A", dA), ("B", dB)):
        c = arm_check(d, F)
        if "error" in c:
            r[f"error_{arm}"] = str(c["error"])[:300]
        else:
            r[f"pf_{arm}"] = c["pf"]
            r[f"check_{arm}"] = {k: c[k] for k in ("n_nonfinite", "n_held", "n_missing")}
            if c["i_bad"]:
                r[f"x_bad_{arm}"] = [float(xu[j]) for j in c["i_bad"]]
    if any(k.startswith("error") for k in r) or not np.any(F):
        if not np.any(F):
            r["error_fixed"] = "比べる列 F が空"
        return r
    if r["pf_A"] and r["pf_B"]:
        eA, eB, tA, tB = _arr(dA, "delta_E"), _arr(dB, "delta_E"), _arr(dA, "delta_target"), _arr(dB, "delta_target")
        d = np.abs(eB - eA)[F] * mm
        i = int(np.argmax(d))
        r["D_mm"], r["i_at_D"] = float(d[i]), int(np.flatnonzero(F)[i])
        r["Dt_mm"] = float(np.max(np.abs(tB - tA)[F]) * mm)
    else:
        r["D_mm"] = None
    return r


def _sigma_norm(sigma) -> dict:
    """classify_ab の sigma を {"status", "mm", "per_step", "why"} にそろえる (float = 測れた・per_step なし、None = 測れない)。"""
    if isinstance(sigma, dict):
        return {"status": sigma.get("status"), "mm": sigma.get("mm"), "per_step": sigma.get("per_step"), "why": sigma.get("why")}
    if sigma is None:
        return {"status": "unmeasurable", "mm": None, "per_step": None, "why": "σ_ref が無い"}
    return {"status": "ok", "mm": float(sigma), "per_step": None, "why": None}


def classify_ab(per: list, sigma, required, tol: float = TOL_MM) -> dict:
    """A/B の結果 (§6 V-h370・§5.1 #36g (4)(5)): 戻り {"result" (結果 A / 結果 B / 判定不能), "why", "code"}。
    per: 出力ごとの `ab_pair` の戻り。required: 判定窓に入るべき出力の step (`stop_rule.required_window_steps`; 無ければ判定不能)。
    sigma: `sigma_ref_measure` の戻り (status ok / unmeasurable / undecidable、mm、per_step) か、試験用に float (測れた・per_step なし)・None。
    順に: 必要な出力の欠け → 判定不能、抽出の失敗 → 判定不能、σ_ref を測れない・判定不能 → 判定不能、どれかの出力で pf_A ≠ pf_B →
    代わりの参照の全部 (per_step の alts の pf_F) でも違いが保たれれば 結果 B・確かめられなければ判定不能、両腕とも pf でない出力
    (F の中の欠測・held・非有限) → 判定不能、それ以外は D: max D + σ ≤ tol → 結果 A、min D − σ > tol → 結果 B、σ を足し引きして初めて
    tol をまたぐ → 判定不能 (参照の未整定)、D だけで窓の中で tol をまたぐ → 判定不能 (時間の変動)。"""
    sg = _sigma_norm(sigma)

    def res(result, code, why):
        return {"result": result, "code": code, "why": why}
    if not required:
        return res("判定不能", "U_required", "判定窓に入るべき出力の集合が無い (出力の間隔・窓から列挙できない)")
    if not per:
        return res("判定不能", "U_empty", "比べる出力が無い")
    have = {r["step"] for r in per if "error_missing" not in r}
    lack = [k for k in required if k not in have]
    if lack:
        return res("判定不能", "U_incomplete", f"判定窓に入るべき {len(required)} 出力のうち比べていない出力 (記録も res も無く再抽出で"
                                              f"補えない・渡されていない) {len(lack)}: {lack[:8]}")
    per = [r for r in per if r["step"] in set(required)]
    err = [r["step"] for r in per if any(k.startswith("error") for k in r)]
    if err:
        return res("判定不能", "U_error", f"抽出・比較の失敗 (出力 {err[:8]})")
    if sg["status"] == "unmeasurable" or (sg["status"] == "ok" and sg["mm"] is None):
        return res("判定不能", "U_sigma_unmeasurable", f"参照の時間の変動 σ_ref を測れない ({sg.get('why')})")
    if sg["status"] != "ok":
        return res("判定不能", "U_sigma_undecidable", f"参照の時間の変動 σ_ref が判定不能 — 代わりの参照で held・欠測・合否が変わる ({sg.get('why')})")
    sig = float(sg["mm"])
    pfd = [r for r in per if r["pf_A"] != r["pf_B"]]
    if pfd:
        ps = sg.get("per_step") or {}
        unconf = []
        for r in pfd:
            alts = (ps.get(r["step"]) or ps.get(str(r["step"])) or {}).get("alts") or []
            if not alts or any(a.get("pf_F") is None or a["pf_F"] == r["pf_A"] for a in alts):
                unconf.append(r["step"])
        detail = "; ".join(f"{r['step']}: A {'合' if r['pf_A'] else '否'}・B {'合' if r['pf_B'] else '否'}" for r in pfd[:6])
        if unconf:
            return res("判定不能", "U_pf_unconfirmed", f"抽出の合否が変わる出力 ({detail}) の違いを、代わりの参照の時刻で確かめられない・保たれない "
                                                     f"(出力 {unconf[:6]})")
        return res("結果 B", "B_pf", f"抽出の合否が変わる出力 {len(pfd)} ({detail}); 代わりの参照の全部の時刻でも違いが保たれる "
                                    "— 旧い参照の流用を退ける")
    inc = [r["step"] for r in per if r.get("D_mm") is None]
    if inc:
        return res("判定不能", "U_fixed_cols", f"両腕とも比べる列 F の中に欠測・held・非有限がある出力 {inc[:8]} — 比較を合格にしない")
    D = [r["D_mm"] for r in per]
    dmax, dmin = max(D), min(D)
    why = f"max_k D_k {dmax:.4g} mm・min_k D_k {dmin:.4g} mm・σ_ref {sig:.4g} mm (許容 {tol} mm、出力 {len(per)})"
    if dmax + sig <= tol:
        return res("結果 A", "A", f"比べる列 F の全断面・判定窓の全出力で max|δ_E,B − δ_E,A| + σ_ref ≤ {tol} mm かつ合否が同じ ({why}) — "
                                "「旧い参照が共通区間で予算を超えるずれを生む」を退ける (旧い参照の出口の不足は残る)")
    if dmin - sig > tol:
        return res("結果 B", "B_D", f"全出力で max|δ_E,B − δ_E,A| − σ_ref > {tol} mm ({why}) — 旧い参照の流用を退ける")
    if dmax <= tol or dmin > tol:
        return res("判定不能", "U_sigma_straddle", f"参照の未整定: σ_ref を足し引きして初めて {tol} mm をまたぐ ({why})")
    return res("判定不能", "U_time", f"窓の中の時間の変動で D_k が {tol} mm をまたぐ ({why})")


def reference_final(ab, do_ab: bool) -> dict:
    """参照の最終の資格 (§5.1 #36g (1)): σ_ref を測れる・判定できる・σ_ref が A/B の分類を左右しない (参照の未整定でない) なら資格あり。
    戻り {"qualified": True / False / None (確かめられない), "why", "sigma_ref_mm", "sigma_status", "ab_code"}。None も更新に進めない。"""
    if not do_ab:
        return {"qualified": None, "why": "A/B と σ_ref を回していない (--no-ab) — 参照の最終の資格を確かめていない"}
    if not ab or "classify" not in ab:
        return {"qualified": None, "why": (ab or {}).get("why") or "A/B を回せない — 参照の最終の資格を確かめていない"}
    c, sg = ab["classify"], ab.get("sigma_ref") or {}
    out = {"sigma_ref_mm": sg.get("mm"), "sigma_ref_use_mm": sg.get("mm_use"), "sigma_status": sg.get("status"), "ab_code": c["code"]}
    if c["code"] in SIGMA_FAIL:
        out.update(qualified=False, why=f"資格なし: {c['why']}")
    elif c["code"] in REF_OK:
        out.update(qualified=True, why=f"資格あり: σ_ref {sg.get('mm')} mm は A/B の分類 ({c['result']}) を左右しない")
    else:
        out.update(qualified=None, why=f"確かめられない: A/B が分類に至らない ({c['why']}) ので σ_ref が判定を左右するかを決められない")
    return out


def sigma_ref_measure(run: Path, steps: list, spec_b: dict, alts: list, ext: dict, x_use, F, dB: dict, mm: float, snap) -> dict:
    """σ_ref (§5.1 #36g (2)): 判定窓の全出力 k について、腕 B の固定の Euler の出力の代わりに alts (予備の資格で固定した代わりの参照
    {res, step, sha256}) で抽出し直し、固定の抽出 dB[k] との δ_E の差を比べる列 F (と使用域の全列) で見る。
    status: ok / unmeasurable (代わりの参照が無い・抽出の例外 = sha の不一致を含む・固定の抽出が無い) / undecidable (代わりの参照で
    F・使用域の held・欠測・合否が変わる)。mm = max_k max_e max_F |Δδ_E| [mm]、mm_use は使用域の全列。"""
    out = {"status": "ok", "mm": None, "mm_use": None, "alts": [a.get("res") for a in alts], "per_step": {}, "why": None, "seconds": 0.0}
    if not alts:
        return {**out, "status": "unmeasurable", "why": "代わりの参照 (予備の資格の sigma_ref_snapshots) が無い"}
    FU = np.ones(len(F), dtype=bool)
    worst, worst_u, unmeas, undec = 0.0, 0.0, [], []
    for k in steps:
        d0 = dB.get(k)
        c0F, c0U = arm_check(d0, F), arm_check(d0, FU)
        if "error" in c0F:
            unmeas.append(f"{k}: 固定の参照の抽出が無い ({c0F['error'][:80]})")
            continue
        e0 = _arr(d0, "delta_E")
        rec = {"pf_F_pin": c0F["pf"], "pf_use_pin": c0U["pf"], "alts": []}
        for a in alts:
            try:
                d = snap(run, k, {"run": spec_b["run"], "res": a["res"], "sha256": a["sha256"]}, ext, x_use)
            except Exception as e:  # noqa: BLE001 — sha の不一致・出力が無いを含む
                unmeas.append(f"{k} × {a.get('res')}: {type(e).__name__}: {e}"[:200])
                rec["alts"].append({"res": a.get("res"), "error": f"{type(e).__name__}: {e}"[:200]})
                continue
            out["seconds"] += float(d.get("seconds") or 0.0)
            cF, cU = arm_check(d, F), arm_check(d, FU)
            if "error" in cF:
                unmeas.append(f"{k} × {a.get('res')}: {cF['error'][:80]}")
                rec["alts"].append({"res": a.get("res"), "error": cF["error"][:200]})
                continue
            e1 = _arr(d, "delta_E")
            both = np.isfinite(e0) & np.isfinite(e1)
            dF = np.abs(e1 - e0)[F & both] * mm
            dU = np.abs(e1 - e0)[both] * mm
            ar = {"res": a.get("res"), "pf_F": cF["pf"], "pf_use": cU["pf"], "max_mm_F": float(dF.max()) if dF.size else None,
                  "max_mm_use": float(dU.max()) if dU.size else None,
                  "change": {q: [c0U[q], cU[q]] for q in ("n_held", "n_missing", "n_nonfinite") if c0U[q] != cU[q]}}
            if cF["pf"] != c0F["pf"] or cU["pf"] != c0U["pf"] or ar["change"]:
                undec.append(f"{k} × {a.get('res')}: 合否 F {c0F['pf']}→{cF['pf']}・使用域 {c0U['pf']}→{cU['pf']} {ar['change']}")
            if ar["max_mm_F"] is not None:
                worst = max(worst, ar["max_mm_F"])
            if ar["max_mm_use"] is not None:
                worst_u = max(worst_u, ar["max_mm_use"])
            rec["alts"].append(ar)
        out["per_step"][k] = rec
    out["mm"], out["mm_use"] = (worst, worst_u)
    if unmeas:
        out.update(status="unmeasurable", why="; ".join(unmeas[:4]) + (" …" if len(unmeas) > 4 else ""))
    elif undec:
        out.update(status="undecidable", why="; ".join(undec[:4]) + (" …" if len(undec) > 4 else ""))
    if out["status"] != "ok":
        out["mm_observed"], out["mm"] = out["mm"], None
    return out


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
    """抽出の参照の A/B と σ_ref (上の docstring の 5、§5.1 #36g (2)〜(5))。"""
    from forge_design.evaluate import stop_rule as SR
    from forge_design.feedback import deltastar_loop as DL
    snap = HOOKS.get("delta_snapshot", SR.delta_snapshot)
    st, ds = rec.get("settings") or {}, rec.get("delta_series") or {}
    W = int(st.get("window", 0) or 0)
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
    # --- (4) 判定窓の完全性: 入るべき出力を出力の間隔と 2 窓から列挙する (記録にある行だけから作らない)
    oi = st.get("out_step_interval")
    try:
        required = SR.required_window_steps(int(s), W, oi)
    except (TypeError, ValueError) as e:
        out["why"] = f"判定窓に入るべき出力を列挙できない (settings.out_step_interval {oi!r}: {e})"
        return out
    rows = {r["step"]: r for r in rec.get("rows") or [] if isinstance(r.get("step"), int)}
    comp = {"required": required, "n_required": len(required), "out_step_interval": int(oi),
            "rows_missing": [k for k in required if k not in rows], "res_missing": [k for k in required if not (run / f"res_{k}.h5").is_file()]}
    out["completeness"] = comp
    spec_a = {"run": str(old.resolve()), "res": f"res_{ors[-1]}.h5", "sha256": sha256_file(old / f"res_{ors[-1]}.h5")}
    out["arm_A"] = spec_a
    # --- (3) 比べる列 F を抽出の前に座標と各腕の Euler の補間できる範囲から決める
    try:
        fx = fixed_columns(x_use, {"A": old, "B": Path(spec_b["run"])}, edge_tol=(ext.get("remap_edge") == "missing"))
    except Exception as e:  # noqa: BLE001
        out["why"] = f"比べる列を決められない ({type(e).__name__}: {e})"[:400]
        return out
    F = fx.pop("mask")
    out["columns"] = {**fx, "note": "抽出の前に座標だけで決めた比べる列 (NS の使用域の列 ∩ 両腕の Euler の補間できる範囲)。この中の欠測・held・"
                                    "非有限は比較を合格にしない"}
    if not F.any():
        out["why"] = "比べる列 F が空 (両腕の Euler の範囲が重ならない)"
        return out
    per, dB_by = [], {}
    for k in required:
        f = run / f"res_{k}.h5"
        if not f.is_file():
            per.append({"step": k, "error_missing": f"出力 res_{k}.h5 が無い (再抽出で補えない)"})
            continue
        sha_k = sha256_file(f)
        dB = (rows.get(k) or {}).get("delta")
        if k not in rows:
            src = "抽出し直し (stop_rule.json に行が無い)"
            dB = None
        elif isinstance(dB, dict) and "delta_E" in dB and len(dB["delta_E"]) == len(x_use) and dB.get("res_sha256") == sha_k:
            src = "記録 (stop_rule.json の delta)"
        else:
            src = "抽出し直し (記録が無いか出力の sha が違う)"
            dB = None
        if dB is None:
            try:
                dB = snap(run, k, spec_b, ext, x_use)
            except Exception as e:  # noqa: BLE001
                dB = {"error": f"{type(e).__name__}: {e}"[:300]}
        try:
            dA = snap(run, k, spec_a, ext, x_use)
        except Exception as e:  # noqa: BLE001
            dA = {"error": f"{type(e).__name__}: {e}"[:300]}
        dB_by[k] = dB
        per.append(ab_pair(k, dA, dB, mm, src, x_use, F))
    comp["reextracted_B"] = [r["step"] for r in per if str(r.get("B_source", "")).startswith("抽出し直し")]
    out["per_output"] = per
    # --- (2) σ_ref: 判定窓の全出力 × 予備の資格で固定した代わりの参照
    alts = (qual or {}).get("sigma_ref_snapshots") or []
    sig = sigma_ref_measure(run, [k for k in required if k in dB_by], spec_b, alts, ext, x_use, F, dB_by, mm, snap)
    sig["selection"] = (qual or {}).get("sigma_ref_selection")
    sig["note"] = ("判定窓の全出力を、固定した Euler の出力の代わりに予備の資格の時点で固定した代わりの参照 (末尾の窓の代表) で抽出し直した "
                   "δ_E の差 (mm は比べる列 F、mm_use は使用域の全列)")
    out["sigma_ref"] = sig
    # スクリーニング (停止の出力 × check_quasisteady の末尾の窓の始まり・中ほど): 記録だけ、σ_ref の上限とは扱わない
    scr = []
    d_s = dB_by.get(s)
    for sn in (qual or {}).get("screening_snapshots") or []:
        try:
            d = snap(run, s, {"run": spec_b["run"], "res": sn["res"], "sha256": sn["sha256"]}, ext, x_use)
            e0, e1 = _arr(d_s or {}, "delta_E"), _arr(d, "delta_E")
            m = F & np.isfinite(e0) & np.isfinite(e1) if len(e0) == len(e1) == len(F) else np.zeros(0, bool)
            scr.append({"res": sn["res"], "max_mm_F": float(np.max(np.abs(e1 - e0)[m]) * mm) if m.any() else None})
        except Exception as e:  # noqa: BLE001
            scr.append({"res": sn.get("res"), "error": f"{type(e).__name__}: {e}"[:300]})
    out["screening"] = {"per_snapshot": scr, "label": "スクリーニング (停止の出力 × 2 出力だけ。σ_ref の上限とは扱わない、記録だけ)"}
    ok = [r for r in per if r.get("D_mm") is not None]
    if ok:
        i_max = max(ok, key=lambda r: r["D_mm"])
        out["D_mm"] = {"max": i_max["D_mm"], "min": min(r["D_mm"] for r in ok), "n_outputs": len(per), "n_with_D": len(ok),
                       "x_at_max": float(x_use[i_max["i_at_D"]]), "step_at_max": i_max["step"], "x_fixed": out["columns"]["x_range"]}
    out["classify"] = classify_ab(per, sig, required)
    out["result"], out["why"] = out["classify"]["result"], out["classify"]["why"]
    n_a = sum(1 for r in per if "error_missing" not in r)
    out["extractions"] = {"arm_A": n_a, "arm_B_reextracted": len(comp["reextracted_B"]),
                          "sigma_ref": sum(len(v.get("alts", [])) for v in sig["per_step"].values()), "screening": len(scr),
                          "sigma_ref_seconds": round(sig.get("seconds", 0.0), 1),
                          "note": "停止の出力の全体の抽出 (at_stop、両腕) は別に 2 回"}
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
    if not pin or not (qual or {}).get("prequalified"):
        und.append("Euler の参照の固定・予備の資格の記録が無い")
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
    # --- 参照の最終の資格 (§5.1 #36g (1): σ_ref を測れない・判定不能・参照の未整定なら更新に進めない) ---
    rf = reference_final(R["ab"], do_ab)
    R["reference_final"] = rf
    # --- 判定 ---
    proceed = bool(gate and gate.get("proceed"))
    met = proceed and rf.get("qualified") is True
    blk = ([] if proceed else ((gate or {}).get("blocked_by") or ["関門の記録が無い"])) + \
          ([] if rf.get("qualified") is True else [f"Euler の参照の最終の資格 (σ_ref): {rf.get('why')}"])
    upd = ("更新に進める条件をすべて満たす (関門と参照の最終の資格; 壁の更新は自動でつながない)" if met
           else "更新に進める条件を満たさない: " + "; ".join(blk)[:900])
    R["update_conditions_met"] = bool(met and not und)
    abr = R["ab"].get("result")
    if und:
        R["verdict"] = ("判定不能: " + "; ".join(und)[:900] + f" (記録: 関門 {('進める' if proceed else '止める') if gate else 'なし'}、"
                        f"参照の最終の資格 {rf.get('qualified')}、A/B {abr})")
        return finish(2)
    R["verdict"] = f"{upd} / 抽出の参照の A/B: {abr}" + (f" ({R['ab'].get('why')})" if R["ab"].get("why") else "")
    return finish(0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--euler-prequalify", action="store_true", help="(i) Euler の参照の予備の資格の判定と固定")
    ap.add_argument("--no-ab", action="store_true",
                    help="抽出の参照の A/B と σ_ref を回さない (A/B は未実施、参照の最終の資格を確かめないので更新に進める条件は満たさない)")
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
    if a.euler_prequalify:
        return euler_prequalify(case, tools)
    return ns_judge(case, tools, not a.no_ab, Path(a.out) if a.out else case / OUT)


if __name__ == "__main__":
    sys.exit(main())
