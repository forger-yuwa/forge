"""plan tooling-nozzle-isothermal-wall-chain §5.1 #32・§6「V-tw370 の事前登録」(commit b5df7416) と「登録の訂正」(commit b6aeb263) の判定。
B の case dir で、3 本の腕 (tw370_run.py) が終わった後に回す (最後に終わった腕が自動で回す)。

  腕   : R0 = run_0811_tw370_g0_c1 (G0・上限 1、参照)、C1 = run_0812_tw370_g0c_c1 (G0c・上限 1)、C50 = run_0813_tw370_g0c_c50 (G0c・上限 50)。
  1. 各腕の前提の検査の記録 (TW370_PRECHECK.json)・本段の直前の照合 (TW370_MAIN_START.json)・基盤の停止 (DISK_STOP・WATCH_ABORT・RESOURCE_STOP)・
     RUN_PROVENANCE.txt のバイナリと環境変数。C50 は本段の入力が C1 の本段の入力と同じ (sha) で、前の段が無く、S1・S2 は C1 の記録を引き継いだ
     履歴として扱う (「各 run に S1・S2 が必須」を C50 に当てない; 訂正 4)。
  2. 段のゲート (R0・C1): S1・S2 の正常終了・最終場の有限性・ρ・P・T > 0・残差の履歴の非有限なし・stage_manifest の最後の区間が本段だけ・
     check_convergence --segment (記録)。
  3. 本段の見張り (m9_watch.json): REACHED は保存した系列から到達を計算し直し、到達の出力の sha を照合。CENSORED は系列に到達が無く上限まで回ったこと。
  4. 比較 (3 組 C1/R0・C50/R0・C50/C1、訂正 1): M1〜M7 は float plan §6.58 と同じ量。スカラー (θ_r・δ_loc の x = 40・70・94、代理 Q_w) は
     連続した混合の許容 |a − b| ≤ max(0.001·|b|, a_q)、a_q = 0.001·|q₃₀₀| (q₃₀₀ は run_0754_t2_d2 res_142500 の値、tw370_q300.json に起動の前に
     固定; 訂正 2)。ほか (壁の分布・出口 M・流量) は sp_judge.py の compare のまま (スカラーを除いて呼ぶ)。fields・compare は sp_judge.py の関数の
     定義を ast で取り出す (sha を照合)。
  5. 抽出の検査 (訂正 5): 固定の帯 (cold_xcheck.common_yb) の 1 倍と 1.25 倍の両方で 3 組を判定し、帯で採否が変われば「抽出に依存、判定保留」。
     生成後の r_t (prepare_info の scale_m = 抽出器の scale)、断面の位置 (x = 40・70・94 の列の x が組で一致)、出口の全節点の (x, r) の対応を
     検査し、不一致・抽出の値の非有限は判定不能 (B と分ける)。帯の感度 (θ_r_s125/θ_r − 1・δ_loc_s125/δ_loc − 1) を 300 K の run_0754 と比べ、
     1 ポイント以上ずれたら「抽出帯の要確認」と書く (判定は変えない)。
  6. 記録 (判定には使わない): check_convergence --segment、θ_r・Q_w の check_quasisteady (§6.11 の閾値)、δ_E (x = 40・70・94)・δ_target
     (x = 40・70・94 と全分布)・出口 M の本段の末尾 2 万 step の 9 出力の時系列と check_quasisteady (方法・位置・窓・引数は下の DSTAR・QS_ARGS で
     固定; 訂正 7)。δ_E・δ_target の末尾の変動は r_t と mm で。時間 (共有の GPU の参考値)。
判定:
  結果 A   = 3 本とも停止規則に到達し、3 組すべてで M1〜M7 が許容内、抽出の検査が成立 →「この壁温・固定壁の停止の時点の運用への G0c・B0 の拡大を支持」。
  結果 B   = 有効な腕の発散・20 万 step で未到達・前の段の発散、または有効な組の (帯の両方で同じ) 許容の超過。上限の間 (C50/C1) だけが外れれば
             上限の変更の互換性を否定、G0c の両腕が R0 から外れれば粗い格子と登録の移植の手順を含む運用の不成立 (格子だけが原因とは言わない)。
  判定保留 = 帯 1 倍と 1.25 倍で組の採否が変わる (抽出に依存)。
  判定不能 = 入力の不一致・出力の欠損・見張りの誤作動・資源・基盤の停止・抽出の不成立 (数値の手法の失敗に数えない)。段のゲートの不成立も A と書かない。
  #32 の合格を壁の更新の許可にしない (訂正 7)。
→ _band_ab/cold_pair/tw370_judge.json
usage: python3 tw370_judge.py [--no-dstar-series]          (判定)
       python3 tw370_judge.py --fix-q300                    (起動の前: 300 K の参照の値を固定 → _band_ab/cold_pair/tw370_q300.json)
"""
import argparse
import ast
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import types
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
HOOKS = {}            # 自己試験用の差し替え (fields_ext・dstar_series)。本番では空
PLAN = "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #32・§6 V-tw370 (事前登録 commit b5df7416、登録の訂正 commit b6aeb263)"
ARMS = {"R0": "run_0811_tw370_g0_c1", "C1": "run_0812_tw370_g0c_c1", "C50": "run_0813_tw370_g0c_c50"}
PAIRS = (("C1", "R0"), ("C50", "R0"), ("C50", "C1"))
DONOR, DONOR_RES, DONOR_STEP = "run_0754_t2_d2", "res_142500.h5", 142500
Q300 = "_band_ab/cold_pair/tw370_q300.json"
OUT = "_band_ab/cold_pair/tw370_judge.json"
XS = (40, 70, 94)
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94")
SCALARS = KEYS + ("delta_loc_40", "delta_loc_70", "delta_loc_94", "Q_w")
S125 = {**{f"theta_r_{x}": f"theta_r_s125_{x}" for x in XS}, **{f"delta_loc_{x}": f"delta_loc_s125_{x}" for x in XS}, "Q_w": "Q_w"}
DROP_RES = set(SCALARS[:-1]) | {"Q_w_proxy", "許容内"}       # compare の戻りのうち、スカラーを中和して呼ぶので捨てる項目
TOL_REL, AQ_FAC = 0.001, 0.001
CORR_TOL_RT = 1e-9                                            # 断面の位置・出口の節点の対応の許容 [r_t]
# 300 K の参照 (run_0754、x = 40/70/94) の帯の感度 [%] (登録の値)。帯を 1.25 倍にした差がこれから 1 ポイント以上ずれたら「抽出帯の要確認」
REG_BAND_300K = {"theta_r": (3.4509, 1.6904, -0.0580), "delta_loc": (-12.108, -4.639, 0.134)}
REG_BAND_TOL_PP = 0.001                                       # 固定の時点で作り直した値と登録の値の一致 (登録の値の桁の丸め)
BAND_FLAG_PP = 1.0
BUDGET, INTERVAL, CONSEC, MIN_POINTS = 200000, 2500, 2, 9
TAIL = 20000
BIN_SHA = "7e95e5bbd8f50c9e"
MSH_SHA = {"R0": "9a891968fd9d0e36", "C1": "285f0fa4ecc0ed1c", "C50": "285f0fa4ecc0ed1c"}
CAP = {"R0": 1.0, "C1": 1.0, "C50": 50.0}
ENV_WANT = {"FORGE_CUDA_BLOCKSIZE": "128", "FORGE_OMEGA_BUDGET": "1", "FORGE_DIAG_COMMIT_LOSS": "500"}
SP_FUNCS = {"fields": "198a8d50ea0ed40b", "compare": "230a1a0def551c33"}
M9_FUNCS = {"drift": "72506550c522d561", "window_points": "d16897d68f8af794", "done_at": "2c66bdc29bd942f5", "spacing_ok": "a2e65d791a97ebe3"}
EULER_REF = "run_0174_euler_v5d_M_r1"
# δ_E・δ_target の抽出 (訂正 7 で固定): 生産の抽出の経路 (V-c45 の cold_pair.py extract と同じ) を各出力に当てる
DSTAR = {"function": "forge_design.feedback.deltastar_loop.extract_and_merge(view, run_0174_euler_v5d_M_r1, omega=0.5, band_select='edge')",
         "view": "一時 dir に prepare_info.json・bcondConfig.yaml・solverConfig.yaml の写しと nozzle.h5・res_<N>.h5 のリンク (1 出力ずつ)",
         "delta_E": "delta_r_next.csv の delta_E (未緩和の平滑化抽出、等温壁は符号付き)",
         "delta_target": "delta_r_next.csv の delta_target (= (1 − ω) δ_in + ω δ_E を同じ P-spline で平滑化、ω 0.5 = run_pass の既定; δ_in は固定壁の入力)",
         "positions_rt": list(XS), "delta_target_full": "delta_r_next.csv の x_rt の全点と試験部 [40, 94] の末尾の最大の振れ",
         "remap_constant_x": False, "window": "本段の末尾 2 万 step の 9 出力 (到達の step n の n−20000〜n、未到達は 180000〜200000)",
         "exit_M": "出口の列の M (= |u|/sonic): r/r_e ≤ 0.8 の平均 (軸の行を除く) と軸の行"}
DSTAR_KW = {"omega": 0.5, "band_select": "edge"}
QS_ARGS = ["--tail", "1.0", "--drift", "0.0001", "--osc", "0.0001"]     # δ_E・δ_target・出口 M の check_quasisteady (訂正 7 で固定、記録)


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


def jwrite(p: Path, d) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(p) + ".tmp")
    tmp.write_text(json.dumps(d, indent=1, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    os.replace(tmp, p)


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


# --- 場の量 (sp_judge の fields + 帯 1.25 倍のスカラー・断面の位置・出口の節点) ----------------------------------------------
def make_fields_ext(case: Path, code: Path, sp_seg: dict):
    """sp_judge.py の fields・compare を case dir の cold_xcheck で使えるようにし、(fields_ext, compare, F) を返す。"""
    if "fields_ext" in HOOKS:
        F = {}
        ns = {"np": np, "F": F, "KEYS": KEYS}
        exec(compile(sp_seg["compare"], str(code / "sp_judge.py"), "exec"), ns)
        return HOOKS["fields_ext"], ns["compare"], F
    sys.path.insert(0, str(code))
    import h5py
    import cold_xcheck as XC
    yb_x, yb = XC.common_yb()
    F = {}
    ns = {"np": np, "h5py": h5py, "XC": XC, "HERE": case, "yb_x": yb_x, "yb": yb, "KEYS": KEYS, "F": F}
    exec(compile(sp_seg["fields"] + "\n\n" + sp_seg["compare"], str(code / "sp_judge.py"), "exec"), ns)
    fields = ns["fields"]

    def fields_ext(run, step):
        o = fields(run, step)
        with h5py.File(case / run / "nozzle.h5", "r") as h:
            xy = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)[:, :2]
        ni, nj, S = XC.mesh_info(xy.shape[0])
        X = (xy[:, 0] / S).reshape(ni, nj); R = (xy[:, 1] / S).reshape(ni, nj)
        for xs in XC.XS:
            i = int(np.argmin(np.abs(o["x"] - xs)))
            o[f"theta_r_s125_{int(xs)}"] = float(o["theta_r_s125"][i]); o[f"delta_loc_s125_{int(xs)}"] = float(o["delta_loc_s125"][i])
            o[f"x_station_{int(xs)}"] = float(o["x"][i])
        o["x_exit"] = X[-1]; o["r_exit_nodes"] = R[-1]; o["scale_S"] = float(S); o["ni_nj"] = (int(ni), int(nj))
        pi = jload(case / run / "prepare_info.json") or {}
        o["prepare_scale_m"] = pi.get("scale_m")
        return o
    return fields_ext, ns["compare"], F


def finite_check(o: dict) -> list:
    """抽出の値の非有限 (判定に使う量だけ)。"""
    bad = []
    for k in SCALARS + tuple(S125.values()):
        if not np.isfinite(float(o[k])):
            bad.append(k)
    for k in ("q_w", "tau_w_signed", "x_w", "M_exit", "r_exit", "mdot_in", "mdot_out", "Hdot_in", "Hdot_out"):
        a = np.asarray(o[k], dtype=float)
        if a.size == 0 or not np.all(np.isfinite(a)):
            bad.append(k)
    return bad


def compare_mixed(F, compare, c, r, band: str, aq: dict) -> dict:
    """組 c/r の M1〜M7。スカラーは連続した混合の許容 (band "1" は帯 1 倍、"1.25" は 1.25 倍の θ_r・δ_loc、Q_w は帯に依らない)、
    ほかは sp_judge.compare (スカラーを 1.0 に中和して呼び、スカラーの項目は捨てる)。"""
    A, B = F[c], F[r]
    A1, B1 = dict(A), dict(B)
    for k in SCALARS:
        A1[k] = 1.0; B1[k] = 1.0
    F["__a"], F["__b"] = A1, B1
    try:
        res0 = compare("__a", "__b")
    finally:
        F.pop("__a", None); F.pop("__b", None)
    ns_ok = bool(res0["許容内"])
    sc, ok = {}, ns_ok
    for k in SCALARS:
        src = k if band == "1" else S125[k]
        a, b = float(A[src]), float(B[src])
        aqk = float(aq[band][k])
        tol = max(TOL_REL * abs(b), aqk)
        d = a - b
        okk = bool(abs(d) <= tol)
        sc[k] = {"quantity": src, "a": a, "b": b, "diff": d, "rel_%": (100 * d / abs(b)) if b != 0 else None, "tol": tol, "a_q": aqk,
                 "tol_kind": "相対 0.1 %" if TOL_REL * abs(b) >= aqk else "絶対 a_q", "ok": okk}
        ok &= okk
    return {"band": band, "scalars": sc, "nonscalar": {k: v for k, v in res0.items() if k not in DROP_RES}, "nonscalar_ok": ns_ok, "許容内": bool(ok)}


def band_sens(o: dict) -> dict:
    """帯の感度 [%]: θ_r_s125/θ_r − 1・δ_loc_s125/δ_loc − 1 (x = 40・70・94)。"""
    return {"theta_r": [100 * (float(o[f"theta_r_s125_{x}"]) / float(o[f"theta_r_{x}"]) - 1) for x in XS],
            "delta_loc": [100 * (float(o[f"delta_loc_s125_{x}"]) / float(o[f"delta_loc_{x}"]) - 1) for x in XS]}


def q300_values(o: dict) -> dict:
    vals = {k: float(o[k]) for k in SCALARS}
    vals.update({S125[k]: float(o[S125[k]]) for k in SCALARS})
    aq = {"1": {k: AQ_FAC * abs(vals[k]) for k in SCALARS}, "1.25": {k: AQ_FAC * abs(vals[S125[k]]) for k in SCALARS}}
    bs = band_sens(o)
    match = all(abs(bs[q][i] - REG_BAND_300K[q][i]) <= REG_BAND_TOL_PP for q in bs for i in range(3))
    return {"values": vals, "a_q": aq, "band_300K_%": bs, "registered_band_300K_%": REG_BAND_300K, "registered_band_values_match": bool(match)}


def fix_q300(case: Path, code: Path, out: Path) -> int:
    """起動の前: run_0754_t2_d2 res_142500 の各断面・各量の値を固定 (a_q = 0.001·|q₃₀₀|)。既にあれば作り直して同じことを確かめるだけ。"""
    sp_seg, _, sp_bad = extract_funcs(code / "sp_judge.py", SP_FUNCS)
    if sp_bad:
        print("sp_judge.py の関数がこの commit の版と違う:", sp_bad); return 2
    fields_ext, _, _ = make_fields_ext(case, code, sp_seg)
    res = case / DONOR / DONOR_RES
    m9 = jload(case / DONOR / "m9_watch.json") or {}
    sha = sha256_file(res)
    if m9.get("status") != "REACHED" or m9.get("reach_step") != DONOR_STEP or m9.get("reach_sha256") != sha:
        print(f"ドナーが停止の出力でない: {m9.get('status')} {m9.get('reach_step')} sha 一致 {m9.get('reach_sha256') == sha}"); return 2
    o = fields_ext(DONOR, DONOR_STEP)
    bad = finite_check(o)
    if bad:
        print("ドナーの抽出の値が非有限:", bad); return 2
    rec = {"plan": PLAN, "created": now(), "donor": {"run": DONOR, "res": DONOR_RES, "step": DONOR_STEP, "sha256": sha},
           "rule": "|a − b| ≤ max(0.001·|b|, a_q)、a_q = 0.001·|q₃₀₀| (登録の訂正 2)", **q300_values(o),
           "station_x": {x: o[f"x_station_{x}"] for x in XS}, "scale_S": o["scale_S"],
           "code": {"sp_judge.py": sha256_file(code / "sp_judge.py")[:16], "cold_xcheck.py": sha256_file(code / "cold_xcheck.py")[:16]}}
    old = jload(out)
    if old is not None:
        same = old.get("donor") == rec["donor"] and old.get("values") == rec["values"] and old.get("a_q") == rec["a_q"]
        print(f"{out} は既にある。作り直した値と {'同じ' if same else '違う (書き換えない)'}")
        return 0 if same else 2
    jwrite(out, rec)
    print(json.dumps({k: rec[k] for k in ("donor", "values", "band_300K_%", "registered_band_values_match")}, ensure_ascii=False, indent=1))
    print(f"→ {out}")
    return 0 if rec["registered_band_values_match"] else 3


# --- δ_E・δ_target・出口 M の時系列 (記録) -----------------------------------------------------------------------------------
def dstar_series(case: Path, run: Path, steps: list, euler: Path) -> dict:
    if "dstar_series" in HOOKS:
        return HOOKS["dstar_series"](case, run, steps, euler)
    import h5py
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    info = jload(run / "prepare_info.json")
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    with h5py.File(run / "nozzle.h5", "r") as h:
        r = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)[:, 1].reshape(ni, nj)[-1]
    rr = r / r[-1]
    m80 = rr <= 0.8; m80[0] = False
    out = {"steps": [], "x_rt": None, "delta_E": [], "delta_target": [], "held": [], "M_core": [], "M_axis": [], "errors": []}
    for s in steps:
        try:
            with tempfile.TemporaryDirectory(prefix="tw370_ds_") as td:
                dd = Path(td) / "ns"; dd.mkdir()
                for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"):
                    shutil.copy(run / f, dd / f)
                os.symlink((run / "nozzle.h5").resolve(), dd / "nozzle.h5")
                os.symlink((run / f"res_{s}.h5").resolve(), dd / f"res_{s}.h5")
                extract_and_merge(dd, euler, **DSTAR_KW)
                nx = read_delta_r_next(dd / "delta_r_next.csv")
            with h5py.File(run / f"res_{s}.h5", "r") as h:
                M = (np.hypot(np.asarray(h["VALUE/Ux"][:], float), np.asarray(h["VALUE/Uy"][:], float))
                     / np.asarray(h["VALUE/sonic"][:], float)).reshape(ni, nj)[-1]
        except Exception as e:  # noqa: BLE001
            out["errors"].append(f"res_{s}: {e!r}")
            continue
        if out["x_rt"] is None:
            out["x_rt"] = np.asarray(nx["x_rt"]).tolist()
        elif not np.array_equal(np.asarray(out["x_rt"]), np.asarray(nx["x_rt"])):
            out["errors"].append(f"res_{s}: delta_r_next.csv の x_rt が前の出力と違う")
            continue
        out["steps"].append(int(s)); out["delta_E"].append(np.asarray(nx["delta_E"]).tolist())
        out["delta_target"].append(np.asarray(nx["delta_target"]).tolist()); out["held"].append(int(np.sum(nx["held"])))
        out["M_core"].append(float(np.mean(M[m80]))); out["M_axis"].append(float(M[0]))
    return out


def summarize_dstar(ser: dict, scale_m: float, run: Path, tools: Path) -> dict:
    """末尾の振れ (r_t と mm)・x = 40・70・94 の値・check_quasisteady (QS_ARGS)。"""
    out = {"n_outputs": len(ser["steps"]), "steps": ser["steps"], "errors": ser["errors"], "held_per_output": ser["held"]}
    if len(ser["steps"]) < 2:
        out["note"] = "出力が 2 つ未満 — 振れを作れない"
        return out
    x = np.asarray(ser["x_rt"]); mm = scale_m * 1e3
    cols = {}
    for q in ("delta_E", "delta_target"):
        A = np.asarray(ser[q])                                    # (出力, x)
        ptp = A.max(0) - A.min(0)
        test = (x >= 40.0) & (x <= 94.0)
        k = int(np.argmax(ptp)); kt = int(np.argmax(np.where(test, ptp, -1)))
        out[q] = {"full_max_ptp_rt": float(ptp[k]), "full_max_ptp_mm": float(ptp[k] * mm), "full_at_x_rt": float(x[k]),
                  "test_max_ptp_rt": float(ptp[kt]), "test_max_ptp_mm": float(ptp[kt] * mm), "test_at_x_rt": float(x[kt]),
                  "at": {}}
        for xs in XS:
            v = np.array([np.interp(xs, x, row) for row in A])
            out[q]["at"][xs] = {"last_rt": float(v[-1]), "last_mm": float(v[-1] * mm), "ptp_rt": float(np.ptp(v)), "ptp_mm": float(np.ptp(v) * mm),
                                "series_rt": v.tolist()}
            cols[f"{q}_{xs}"] = v
    for q in ("M_core", "M_axis"):
        v = np.asarray(ser[q]); cols[q] = v
        out[q] = {"last": float(v[-1]), "ptp": float(np.ptp(v)), "series": v.tolist()}
    sc = run / "tw370_dstar_tail.csv"
    with open(sc, "w") as fh:
        fh.write("step," + ",".join(cols) + "\n")
        for i, s in enumerate(ser["steps"]):
            fh.write(f"{s}," + ",".join(f"{cols[c][i]:.12g}" for c in cols) + "\n")
    q = subprocess.run([sys.executable, str(tools / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", ",".join(cols)]
                       + QS_ARGS, capture_output=True, text=True)
    out["check_quasisteady"] = {"args": QS_ARGS, "lines": [l.strip() for l in q.stdout.splitlines() if any(c in l for c in cols) and ":" in l]
                                or (q.stdout + q.stderr).strip().splitlines()[-6:]}
    return out


# --- 腕ごとの状態 ---------------------------------------------------------------------------------------------------------------
def arm_state(name: str, case: Path, tools: Path, m9ns: dict, m9_bad: list, c1_state: dict | None) -> dict:
    """1 本の腕の前提・段・見張りを読む。戻り {und, b, gfail, reached, reach_step, status, ...}。und は腕の判定不能の理由 (その腕の B を無効にする)。"""
    run = case / ARMS[name]
    A = {"run": ARMS[name], "und": [], "b": [], "gfail": [], "reached": False, "reach_step": None, "status": None, "stages": {}, "notes": []}
    und, bres, gfail, S = A["und"], A["b"], A["gfail"], A["stages"]
    pc = jload(run / "TW370_PRECHECK.json")
    if not run.is_dir():
        und.append(f"{name}: run {run.name} が無い"); return A
    if name == "C50" and (run / "C50_NOT_STARTED.json").is_file():
        ns = jload(run / "C50_NOT_STARTED.json") or {}
        A["not_started"] = ns
        if c1_state is not None and c1_state["b"] and not c1_state["und"]:
            A["notes"].append(f"C50: C1 が本段に入らなかったので分岐できない (C1 の失敗の帰結、独立の理由にしない): {ns.get('reason')}")
        else:
            und.append(f"C50: C1 から分岐できない ({ns.get('reason')})")
        return A
    if pc is None:
        und.append(f"{name}: 前提の検査の記録 TW370_PRECHECK.json が無い"); return A
    A["precheck"] = {"ok": pc.get("ok"), "failures": pc.get("failures"),
                     "items_ok": {k: v.get("ok") for k, v in (pc.get("items") or {}).items() if isinstance(v, dict) and "ok" in v},
                     "ar_gt5000": ((pc.get("items") or {}).get("3_quality") or {}).get("ar_gt5000"),
                     "transfer": {k: ((pc.get("items") or {}).get("7_transfer") or {}).get(k) for k in ("tool", "verdict", "correction_after_transfer",
                                                                                                       "nozzle_sha256_after", "nearest", "omega_floor")},
                     "start_identity": ((pc.get("items") or {}).get("8_start_identity") or {}).get("ok"),
                     "sha": (pc.get("items") or {}).get("4_sha"), "resources_start": pc.get("resources_start")}
    if not pc.get("ok"):
        und.append(f"{name}: 前提の検査が不成立: " + "; ".join(pc.get("failures") or ["(理由なし)"])[:500])
        return A
    if not (run / "nozzle.msh").is_file() or not sha256_file(run / "nozzle.msh").startswith(MSH_SHA[name]):
        und.append(f"{name}: 判定の時点の nozzle.msh が {MSH_SHA[name]} でない")
    for f, why in (("DISK_STOP", "ディスクの空きで自分の forge を止めた"), ("WATCH_ABORT", "本段の直前の照合が不成立"),
                   ("RESOURCE_STOP", "資源の不足")):
        if (run / f).is_file():
            und.append(f"{name}: {why} ({f}: {(run / f).read_text().strip()[:300]})")
    rc_txt = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").is_file() else None
    run_end = jload(run / "RUN_END.json") or {}
    A["RUN_RC"] = rc_txt
    A["run_end_exception"] = (run_end.get("exception") or "")[-1200:] or None
    mss = jload(run / "MAIN_STAGE_STARTED")
    vms = jload(run / "TW370_MAIN_START.json")
    A["mss"] = mss
    prov = (run / "RUN_PROVENANCE.txt").read_text() if (run / "RUN_PROVENANCE.txt").is_file() else ""
    if mss is not None:
        if prov:
            sha_l = next((l.split(":", 1)[1].strip() for l in prov.splitlines() if l.startswith("forge_sha256")), "")
            envs = dict(l.split(":", 1)[1].strip().split("=", 1) for l in prov.splitlines() if l.startswith("env ") and "=" in l)
            A["provenance"] = {"forge_sha256": sha_l[:16], "env": envs}
            env_bad = [k for k, v in envs.items() if k != "FORGE_BIN" and ENV_WANT.get(k) != v]
            if not sha_l.startswith(BIN_SHA) or envs.get("FORGE_CUDA_BLOCKSIZE") != "128" or env_bad:
                und.append(f"{name}: 本段の RUN_PROVENANCE.txt のバイナリ・環境変数が登録と違う ({sha_l[:16]}、{envs})")
        else:
            und.append(f"{name}: 本段の RUN_PROVENANCE.txt が無い")
    # --- 段 ---
    stg = jlines(run / "tw370_stages.jsonl")
    gates = jlines(run / "tw370_stage_gate.jsonl")
    ends = {r["tag"]: r for r in stg if r.get("event") == "end"}
    A["stage_ends"] = ends
    hist = {t: residual_nonfinite(run / f"residual_history_{t}.csv") for t in ("S1_soft", "S2_mid", "main")}
    if hist["main"] is None and mss is not None:
        hist["main"] = residual_nonfinite(run / "residual_history.csv")
    S["residual_history"] = hist
    if name in ("R0", "C1"):
        for t in ("S1_soft", "S2_mid"):
            g = next((x for x in gates if x.get("tag") == t), None)
            e = ends.get(t)
            st = (g or {}).get("stats") or {}
            rec = {"rc": None if e is None else e.get("rc"), "seconds": None if e is None else e["t1"] - e["t0"],
                   "gate_problems": None if g is None else g.get("problems"), "stats": st or (g or {}).get("stats_error"),
                   "res_sha256": None if g is None else (g.get("res_sha256") or "")[:16]}
            ok = (e is not None and e.get("rc") == 0 and g is not None and not g.get("problems") and bool(st)
                  and st.get("nonfinite_total") == 0 and all(st.get(f"{k}_n_nonpositive") == 0 and not st.get(f"{k}_missing") for k in ("ro", "P", "T")))
            h = hist[t]
            if h is None or h["nonfinite"]:
                ok = False
            rec["ok"] = ok
            S[t] = rec
            if not ok and mss is not None:
                gfail.append(f"{name} 段 {t}: 正常終了・有限性・ρ・P・T・残差の履歴のどれかが不成立 ({rec}, 残差 {h})")
    else:
        S["inherited_from_C1"] = "S1・S2 は C1 から引き継いだ履歴 (C1_HISTORY/、段のゲートは C1 の腕で見る; 訂正 4)"
        A["c1_history_files"] = sorted(f.name for f in (run / "C1_HISTORY").iterdir()) if (run / "C1_HISTORY").is_dir() else None
    if mss is None:
        if name in ("R0", "C1"):
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
            if ev and not und:
                bres.append(f"{name}: 前の段で発散し本段に入れない: " + "; ".join(ev)[:400])
            elif not und:
                und.append(f"{name}: 本段が始まらなかった (前の段の失敗に発散の証拠が無い; RUN_RC {rc_txt}; 例外 {(run_end.get('exception') or '')[-200:]!r})")
        else:
            und.append(f"C50: 本段が始まらなかった (RUN_RC {rc_txt}; 例外 {(run_end.get('exception') or '')[-200:]!r})")
    else:
        if vms is None or not vms.get("ok"):
            und.append(f"{name}: 本段の直前の照合の記録 TW370_MAIN_START.json が無いか不成立 ({None if vms is None else vms.get('failures')})")
        if mss.get("res_files_present"):
            und.append(f"{name}: 本段の通知の時点で前の段の出力が残っていた: {mss.get('res_files_present')}")
        want_prev = ["S1_soft", "S2_mid"] if name != "C50" else []
        if [s.get("tag") for s in mss.get("previous_stages", [])] != want_prev:
            und.append(f"{name}: 本段の通知の前の段 {[s.get('tag') for s in mss.get('previous_stages', [])]} が {want_prev} でない")
        if name == "C50":
            c1m = (c1_state or {}).get("mss") or {}
            if (mss.get("initial_field") or {}).get("sha256") != (c1m.get("initial_field") or {}).get("sha256") or not c1m:
                und.append("C50: 本段の初期の場の sha256 が C1 の本段の初期の場と違う (受け渡しの不一致)")
            if not A["precheck"].get("start_identity"):
                und.append("C50: 本段の入力と C1 の本段の入力の照合 (8) が記録に無いか不成立")
        S["main_stage_started"] = {k: mss.get(k) for k in ("time", "time_unix", "main_config", "initial_field", "previous_stages", "prestage_config")}
        if hist["main"] is None:
            gfail.append(f"{name}: 本段の残差の履歴が無い")
        elif hist["main"]["nonfinite"]:
            gfail.append(f"{name}: 本段の残差の履歴に非有限 {hist['main']['nonfinite']} (最初の step {hist['main']['first_nonfinite_step']})")
        try:
            if str(tools) not in sys.path:
                sys.path.insert(0, str(tools))
            import stage_manifest as SM
            man = json.loads((run / "stage_manifest.json").read_text())
            segs = SM.segments(man, SM.load_launches(str(run)))
            tags = [[s["tag"] for s in sg] for sg in segs]
            S["manifest_segments"] = tags
            if not tags or tags[-1] != ["main"] or (name == "C50" and tags != [["main"]]):
                gfail.append(f"{name}: stage_manifest の最後の区間が本段だけでない: {tags}")
        except Exception as e:  # noqa: BLE001
            gfail.append(f"{name}: stage_manifest の区間を読めない ({e!r})")
        p = subprocess.run([sys.executable, str(tools / "check_convergence.py"), str(run), "--segment"], capture_output=True, text=True)
        lines = (p.stdout + p.stderr).splitlines()
        S["check_convergence_segment"] = {"rc": p.returncode, "overall": next((l for l in lines if "OVERALL" in l), "(行なし)"),
                                          "verdict_lines": [l.strip()[:160] for l in lines if "VERDICT" in l or "<--" in l][:20]}
    # --- 見張り ---
    m9 = jload(run / "m9_watch.json")
    status = None if m9 is None else m9.get("status")
    A.update(status=status, m9=m9)
    if mss is not None:
        rows = sorted((m9 or {}).get("rows", []), key=lambda r: r["step"])
        if m9 is None:
            und.append(f"{name}: 見張りの記録 m9_watch.json が無い")
        elif m9.get("phase") != "line" or m9.get("budget") != BUDGET or m9.get("consec") != CONSEC or m9.get("interval") != INTERVAL \
                or m9.get("min_points") != MIN_POINTS or m9.get("keep_tail") != 9:
            und.append(f"{name}: 見張りの引数が §6.11 と違う")
        elif m9_bad and status in ("REACHED", "CENSORED"):
            und.append(f"{name}: m9_watch.py の到達の計算がこの commit の版と違う: " + "; ".join(m9_bad))
        elif status == "REACHED":
            n = m9.get("reach_step")
            sp_ok = m9ns["spacing_ok"](rows)
            d = m9ns["done_at"](rows)
            A["recheck"] = {"spacing_ok": sp_ok, "done_at": d}
            fpath = run / f"res_{n}.h5"
            if not sp_ok or d != ("REACH", n):
                und.append(f"{name}: 到達を系列から計算し直すと合わない (見張りの誤作動): 間隔 {sp_ok}、done_at {d}、記録 {n}")
            elif not fpath.is_file():
                und.append(f"{name}: 到達の出力 {fpath.name} が無い")
            elif sha256_file(fpath) != m9.get("reach_sha256"):
                und.append(f"{name}: 到達の出力 {fpath.name} の sha256 が見張りの記録と違う")
            else:
                A["reached"], A["reach_step"] = True, int(n)
        elif status == "CENSORED":
            sp_ok = m9ns["spacing_ok"](rows)
            d = m9ns["done_at"](rows)
            last = rows[-1]["step"] if rows else None
            A["recheck"] = {"spacing_ok": sp_ok, "done_at": d, "last_step": last}
            if not sp_ok or d is not None or last is None or last < BUDGET or rc_txt != "0":
                und.append(f"{name}: CENSORED を系列から確かめられない (見張りの誤作動): 間隔 {sp_ok}、done_at {d}、最後の step {last}、RUN_RC {rc_txt}")
            elif not und:
                bres.append(f"{name}: 本段 {BUDGET} step までに停止規則に到達しない (CENSORED)")
        elif status == "DIVERGED":
            if not und:
                bres.append(f"{name}: 本段で発散 (DIVERGED、step {m9.get('fail_step')}、{m9.get('note') or ''})")
        else:
            und.append(f"{name}: 見張りの状態 {status} (予期しない終了・見張りの誤作動・資源は数値の手法の失敗に数えない; RUN_RC {rc_txt}, "
                       f"{(m9 or {}).get('note') or (m9 or {}).get('rc')})")
    A["termination"] = ("見張りの停止 (REACHED)" if status == "REACHED" else "上限まで回った (CENSORED)" if status == "CENSORED"
                        else f"予期しない終了 ({status}, RUN_RC {rc_txt})" if status else "見張りの記録なし")
    return A


def timing_rec(run: Path, m9: dict, status: str, pc_t_prepared=None) -> dict:
    T = {}
    try:
        t0 = float((run / "START_EPOCH").read_text().strip())
        tm = jload(run / "MAIN_STAGE_STARTED") or {}
        t_main = float(tm["time_unix"])
        rows_t = [r for r in csv.reader(open(run / "timing.tsv"), delimiter="\t") if len(r) >= 3]
        allr = [(float(r[0]), int(r[2]), (r[3].strip() if len(r) > 3 else "")) for r in rows_t]
        mr = [(float(r[0]), int(r[1])) for r in rows_t if r[1].strip() and float(r[0]) >= t_main]
        n_end = int(m9["reach_step"]) if status == "REACHED" else (BUDGET if status == "CENSORED" else None)
        T.update(main_start_after_s=t_main - t0, samples=len(allr), main_samples=len(mr))
        if pc_t_prepared:
            T["prepare_s"] = pc_t_prepared - t0
        if n_end and mr:
            t = np.array([x[0] for x in mr]); s = np.array([x[1] for x in mr])
            if s.max() >= n_end:
                k = int(np.argmax(s >= n_end))
                tr = t[k] if k == 0 or s[k] == s[k - 1] else t[k - 1] + (t[k] - t[k - 1]) * (n_end - s[k - 1]) / (s[k] - s[k - 1])
                j = int(np.argmax(s >= INTERVAL))
                T.update(n_end=n_end, T_total_s=float(tr - t0), T_main_s=float(tr - t_main),
                         per_step_ms=float(1e3 * (tr - t[j]) / max(n_end - s[j], 1)),
                         n_forge_max=int(max((nf for tt, nf, _ in allr if tt <= tr), default=0)),
                         others=sorted({o for tt, _, o in allr if o and tt <= tr}))
            else:
                T["note"] = f"timing.tsv の本段の行が到達の step {n_end} まで無い (最大 {int(s.max())})"
    except Exception as e:  # noqa: BLE001
        T["error"] = repr(e)
    return T


def interpret_b(pairs_fail: set) -> str:
    """結果 B の読み方 (訂正 1)。"""
    if pairs_fail == {"C50/C1"} or ("C50/C1" in pairs_fail and "C1/R0" not in pairs_fail):
        return "上限の間 (C50/C1) が外れた → 上限の変更 (1 → 50) の互換性を否定"
    if {"C1/R0", "C50/R0"} <= pairs_fail:
        return "G0c の両腕が R0 から外れた → 粗い格子と登録の移植の手順を含む運用の不成立 (格子だけが原因とは言わない)"
    if pairs_fail:
        return f"許容の超過の組 {sorted(pairs_fail)}"
    return ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", default=str(HERE))
    ap.add_argument("--code-dir", default=str(HERE), help="sp_judge.py・m9_watch.py・cold_xcheck.py の場所 (既定: この台本の場所)")
    ap.add_argument("--tools", default=None, help="forge の tools (既定: ~/forge-g0cchain の、無ければこの作業ツリーの)")
    ap.add_argument("--design", default=None, help="forge_design の design/ (既定: ~/forge-g0cchain の、無ければこの作業ツリーの)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--fix-q300", action="store_true", help="起動の前: 300 K の参照の値を固定する")
    ap.add_argument("--no-dstar-series", action="store_true", help="δ_E・δ_target の時系列 (記録) を作らない")
    a = ap.parse_args(argv)
    case, code = Path(a.case_dir).resolve(), Path(a.code_dir).resolve()
    tools = Path(a.tools) if a.tools else (Path.home() / "forge-g0cchain/solver_density_cuda/tools")
    if not tools.is_dir():
        tools = HERE.parents[1] / "solver_density_cuda/tools"
    design = Path(a.design) if a.design else (Path.home() / "forge-g0cchain/design")
    if not design.is_dir():
        design = HERE.parents[1] / "design"
    sys.path.insert(0, str(design))
    if a.fix_q300:
        return fix_q300(case, code, case / Q300)
    out_p = Path(a.out) if a.out else case / OUT
    R = {"plan": PLAN, "created": now(), "arms": dict(ARMS), "case_dir": str(case), "code_dir": str(code), "tools": str(tools), "design": str(design),
         "undecidable": [], "b_reasons": [], "gate_failures": [], "notes": [], "per_arm": {}, "extraction": {}, "pairs": {}, "band_dependence": [],
         "band_sensitivity": {}, "quasisteady": {}, "dstar_series": {}, "timing": {}, "records": {}, "fixed": {"DSTAR": DSTAR, "QS_ARGS": QS_ARGS,
                                                                                                           "tol": "|a − b| ≤ max(0.001·|b|, 0.001·|q₃₀₀|)"}}
    und, bres, gfail = R["undecidable"], R["b_reasons"], R["gate_failures"]
    sp_seg, sp_sha, sp_bad = extract_funcs(code / "sp_judge.py", SP_FUNCS)
    m9_seg, m9_sha, m9_bad = extract_funcs(code / "m9_watch.py", M9_FUNCS)
    R["records"]["functions"] = {"sp_judge.py": {"sha256": sha256_file(code / "sp_judge.py")[:16], "funcs": sp_sha, "mismatch": sp_bad},
                                 "m9_watch.py": {"sha256": sha256_file(code / "m9_watch.py")[:16], "funcs": m9_sha, "mismatch": m9_bad}}
    m9ns = {"np": np, "a": types.SimpleNamespace(consec=CONSEC, min_points=MIN_POINTS, interval=INTERVAL), "REF": None, "KEYS": KEYS + ("Q_w",)}
    if not m9_bad:
        exec(compile("\n\n".join(m9_seg[k] for k in M9_FUNCS), str(code / "m9_watch.py"), "exec"), m9ns)

    # --- 1〜3. 腕ごと (C1 を先に読む: C50 の受け渡しの照合と「C1 の失敗の帰結」の判断に使う) ---
    st = {}
    st["C1"] = arm_state("C1", case, tools, m9ns, m9_bad, None)
    st["R0"] = arm_state("R0", case, tools, m9ns, m9_bad, None)
    st["C50"] = arm_state("C50", case, tools, m9ns, m9_bad, st["C1"])
    for k in ("R0", "C1", "C50"):
        A = st[k]
        und += A["und"]; gfail += A["gfail"]; R["notes"] += A["notes"]
        if not A["und"]:
            bres += A["b"]
        R["per_arm"][k] = {kk: v for kk, v in A.items() if kk not in ("m9", "mss")}
        R["per_arm"][k]["m9_rows"] = len((A.get("m9") or {}).get("rows", []))

    # --- 4・5. 比較と抽出の検査 ---
    ok_arms = [k for k in ("R0", "C1", "C50") if st[k]["reached"] and not st[k]["und"]]
    cmp_und = []
    F, fields_ext, compare = {}, None, None
    q300 = jload(case / Q300)
    if len(ok_arms) >= 2:
        if sp_bad:
            cmp_und.append("sp_judge.py の fields・compare がこの commit の版と違う: " + "; ".join(sp_bad))
        donor_sha = sha256_file(case / DONOR / DONOR_RES) if (case / DONOR / DONOR_RES).is_file() else None
        if q300 is None or q300.get("donor", {}).get("sha256") != donor_sha or not q300.get("registered_band_values_match"):
            cmp_und.append(f"{Q300} が無いか、ドナーの sha が違うか、300 K の帯の感度が登録と合わない")
        if not cmp_und:
            fields_ext, compare, F = make_fields_ext(case, code, sp_seg)
            try:
                o300 = fields_ext(DONOR, DONOR_STEP)
                re300 = q300_values(o300)
                R["records"]["q300_recomputed_same"] = (re300["values"] == q300["values"])
                if re300["values"] != q300["values"]:
                    cmp_und.append("q₃₀₀ を作り直すと固定した値と違う (抽出器・入力が起動の前と違う)")
                R["band_sensitivity"]["run_0754 (300 K)"] = re300["band_300K_%"]
            except (OSError, KeyError, SystemExit, ValueError) as e:
                cmp_und.append(f"300 K の参照 {DONOR} の値を作り直せない ({e!r})")
        EX = R["extraction"]
        good = []                                         # 抽出の検査が通った腕 (組はこの腕どうしだけ作る)
        for k in ok_arms if not cmp_und else []:
            try:
                F[k] = fields_ext(ARMS[k], st[k]["reach_step"])
            except (OSError, KeyError, SystemExit, ValueError) as e:
                cmp_und.append(f"{k} の停止の出力を読めない ({e!r})"); continue
            o = F[k]
            ex = {"nonfinite_fields_h5": o.get("nonfinite"), "nonfinite_extracted": finite_check(o), "scale_S": o.get("scale_S"),
                  "prepare_scale_m": o.get("prepare_scale_m"), "station_x": {x: o.get(f"x_station_{x}") for x in XS}, "ni_nj": o.get("ni_nj"),
                  "failures": []}
            EX[k] = ex
            if o.get("nonfinite"):
                ex["failures"].append(f"{k} の停止の出力に非有限 {o['nonfinite']}")
            if ex["nonfinite_extracted"]:
                ex["failures"].append(f"{k} の抽出の値が非有限: {ex['nonfinite_extracted']}")
            if ex["prepare_scale_m"] is None or abs(float(ex["prepare_scale_m"]) - float(ex["scale_S"])) > 1e-12 * abs(float(ex["scale_S"])):
                ex["failures"].append(f"{k} の r_t (prepare_info の scale_m {ex['prepare_scale_m']}) が抽出器の scale {ex['scale_S']} と違う")
            for x in XS:
                if abs(float(o[f"x_station_{x}"]) - x) > 0.5:
                    ex["failures"].append(f"{k} の断面 x = {x} の列が {o[f'x_station_{x}']} (0.5 r_t より離れている)")
            cmp_und += ex["failures"]
            if not ex["failures"]:
                good.append(k)
                R["band_sensitivity"][k] = band_sens(o)
        computable = [(c, r) for c, r in PAIRS if c in good and r in good]
        for c, r in computable:
            tag = f"{c}/{r}"
            a_, b_ = F[c], F[r]
            corr = {"station_dx_rt": {x: float(a_[f"x_station_{x}"]) - float(b_[f"x_station_{x}"]) for x in XS}}
            same_n = np.asarray(a_["x_exit"]).shape == np.asarray(b_["x_exit"]).shape
            corr["exit_n_nodes"] = [int(np.asarray(a_["x_exit"]).size), int(np.asarray(b_["x_exit"]).size)]
            if same_n:
                corr["exit_max_dx_rt"] = float(np.max(np.abs(np.asarray(a_["x_exit"]) - np.asarray(b_["x_exit"]))))
                corr["exit_max_dr_rt"] = float(np.max(np.abs(np.asarray(a_["r_exit_nodes"]) - np.asarray(b_["r_exit_nodes"]))))
            bad = (not same_n or corr["exit_max_dx_rt"] > CORR_TOL_RT or corr["exit_max_dr_rt"] > CORR_TOL_RT
                   or any(abs(v) > CORR_TOL_RT for v in corr["station_dx_rt"].values()))
            EX[f"correspondence {tag}"] = corr
            if bad:
                cmp_und.append(f"{tag}: 断面の位置か出口の節点の (x, r) の対応が {CORR_TOL_RT} r_t を超えて違う ({corr})")
                continue
            try:
                p1 = compare_mixed(F, compare, c, r, "1", q300["a_q"])
                p125 = compare_mixed(F, compare, c, r, "1.25", q300["a_q"])
            except Exception as e:  # noqa: BLE001
                cmp_und.append(f"{tag}: 比較の失敗 ({e!r})"); continue
            R["pairs"][tag] = {"band_1": p1, "band_1.25": p125, "許容内 (帯 1)": p1["許容内"], "許容内 (帯 1.25)": p125["許容内"]}
            if p1["許容内"] != p125["許容内"]:
                R["band_dependence"].append(tag)
        und += cmp_und
        # 帯の感度の 300 K との比較 (記録; R0 は登録の「抽出帯の要確認」)
        for k in ("R0", "C1", "C50"):
            if k in R["band_sensitivity"]:
                bs = R["band_sensitivity"][k]
                dpp = {q: [bs[q][i] - REG_BAND_300K[q][i] for i in range(3)] for q in bs}
                flag = any(abs(v) >= BAND_FLAG_PP for q in dpp for v in dpp[q])
                R["band_sensitivity"][f"{k} − 300 K [pp]"] = dpp
                if flag:
                    R["notes"].append(f"{k}: 抽出帯の要確認 (帯 1.25 倍の感度が 300 K の参照から 1 ポイント以上ずれた、判定は変えない): {dpp}")
    pairs_fail = {t for t, v in R["pairs"].items() if not v["許容内 (帯 1)"] and t not in R["band_dependence"]}
    if pairs_fail:
        bres.append(f"許容の超過 (帯 1 倍と 1.25 倍の両方で): {sorted(pairs_fail)}")

    # --- 6. 記録 ---
    euler = case / EULER_REF
    for k in ("R0", "C1", "C50"):
        A = st[k]; run = case / ARMS[k]; m9 = A.get("m9") or {}
        if A["status"] not in ("REACHED", "CENSORED") or A.get("mss") is None:
            continue
        n = int(m9["reach_step"]) if A["status"] == "REACHED" else BUDGET
        rows = sorted([x for x in m9.get("rows", []) if n - TAIL <= x["step"] <= n], key=lambda x: x["step"])
        Q = R["quasisteady"][k] = {"window": [n - TAIL, n], "points": len(rows)}
        if len(rows) != 9:
            Q["note"] = f"末尾 2 万 step の出力が 9 でない ({len(rows)})"
        try:
            sc = run / "tw370_series_tail.csv"
            with open(sc, "w") as fh:
                fh.write("step," + ",".join(KEYS + ("Q_w",)) + "\n")
                for x in rows:
                    fh.write(f"{x['step']}," + ",".join(f"{x[kk]:.10g}" for kk in KEYS + ("Q_w",)) + "\n")
            for tag, cols, thr in (("§6.11 θ_r (0.0001)", KEYS, "0.0001"), ("§6.11 Q_w (0.0002)", ("Q_w",), "0.0002")):
                q = subprocess.run([sys.executable, str(tools / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", ",".join(cols),
                                    "--tail", "1.0", "--drift", thr, "--osc", thr], capture_output=True, text=True)
                Q[tag] = [l.strip() for l in q.stdout.splitlines() if any(c in l for c in cols) and ":" in l] or (q.stdout + q.stderr).strip().splitlines()[-5:]
        except Exception as e:  # noqa: BLE001
            Q["error"] = repr(e)
        if not a.no_dstar_series:
            try:
                steps = [int(x["step"]) for x in rows if (run / f"res_{int(x['step'])}.h5").is_file()]
                ser = dstar_series(case, run, steps, euler)
                info = jload(run / "prepare_info.json") or {}
                R["dstar_series"][k] = summarize_dstar(ser, float(info.get("scale_m", float("nan"))), run, tools)
                R["dstar_series"][k]["euler_ref"] = {"run": EULER_REF, "files_sha256": {f.name: sha256_file(f)[:16] for f in sorted(euler.glob("*"))
                                                                                       if f.is_file() and (f.suffix == ".h5" or f.name == "prepare_info.json")}}
            except Exception as e:  # noqa: BLE001
                R["dstar_series"][k] = {"error": repr(e)}
        R["timing"][k] = timing_rec(run, m9, A["status"], (jload(run / "TW370_PRECHECK.json") or {}).get("t_prepared"))

    # --- 判定 ---
    reached_all = all(st[k]["reached"] and not st[k]["und"] for k in ("R0", "C1", "C50"))
    if bres:
        verdict = "結果 B (この壁温・固定壁の停止の時点の運用への G0c・B0 の拡大を見送る): " + "; ".join(bres)[:700]
        it = interpret_b(pairs_fail)
        if it:
            verdict += f" — {it}"
        if und:
            verdict += f" (ほかに判定不能の理由 {len(und)} 件: {'; '.join(und)[:300]})"
    elif und:
        verdict = "判定不能: " + "; ".join(und)[:900]
    elif R["band_dependence"]:
        verdict = f"判定保留 (抽出の帯に依存): 帯 1 倍と 1.25 倍で採否が変わる組 {R['band_dependence']}"
    elif gfail:
        verdict = "判定不能 (到達・比較は A の条件だが、段のゲートに不成立があるので A と書かない): " + "; ".join(gfail)[:900]
    elif reached_all and len(R["pairs"]) == len(PAIRS) and all(v["許容内 (帯 1)"] and v["許容内 (帯 1.25)"] for v in R["pairs"].values()):
        verdict = ("結果 A: 3 本とも停止規則に到達し、C1/R0・C50/R0・C50/C1 の 3 組すべてで M1〜M7 が許容内、抽出の検査が成立 "
                   "(この壁温・固定壁の停止の時点の運用への G0c・B0 の拡大を支持。壁の更新の許可ではない)")
    else:
        verdict = "判定不能 (想定外の分岐: 到達・比較の記録がそろわない)"
    R["verdict"] = verdict
    for t, v in R["pairs"].items():
        print(t, "帯 1:", "許容内" if v["許容内 (帯 1)"] else "超過", " 帯 1.25:", "許容内" if v["許容内 (帯 1.25)"] else "超過",
              {k: (round(s["rel_%"], 4) if s["rel_%"] is not None else None, s["tol_kind"]) for k, s in v["band_1"]["scalars"].items()})
    print("時間", json.dumps(R["timing"], ensure_ascii=False, default=float)[:1500])
    print(f"== VERDICT V-tw370: {verdict}")
    jwrite(out_p, R)
    print(f"→ {out_p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
