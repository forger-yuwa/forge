"""plan architecture-float-state-double-geometry §6.62 (試験 P、事前登録・2026-10-13 の資格の訂正の後) の資格 Q1〜Q6 の検査。
腕 P0 = run_0806_tp_p0 (lowMachPrecond 0)、P3 = run_0807_tp_p3 (lowMachPrecond 3)、資格の 1 step の run Q0 = run_0808_tp_q0・Q3 = run_0809_tp_q3。
  Q1 同じ起点・設定・バイナリ: 入力 (_g0c/q32_init.h5 と各 run の nozzle.h5) の sha 6f8fcda6d97482e3、4 本の res_0 の 9 保存量が互いにビット一致
     (入力との差は節点の種類ごとの不一致の数を記録だけ)、P0 と P3 の設定の差が lowMachPrecond の 1 行だけ (資格の run は腕と nStepOuter の 1 行だけ、
     ほかの yaml は 4 本で同一、登録の値 (nStepOuter 20,000・outStepInterval 500・cfl・cfl_pseudo 4・implicitRelax 0.7・precondEps 0.15・
     implicitSolvePrecision 0・extraFields の 8 量) が書いてある)、RUN_PROVENANCE の forge の sha が 75ee4675b0c90f41。
  Q2 実効値 (forge_run.log、forge-25 の cg_pp2_check.py の cfg_val と同じ正規表現): lowMachPrecond 0/3、precondEps 0.15、implicitSolvePrecision 0 (4 本)。
     speciesImplicitCoupling・speciesPrecondDt は記録 (yaml に無ければ forge は出さない (solverConfig.cpp:24-36) → ソースの既定 0・1 と書く)。
  Q3 刻みの効き (cg_pp2_check.py の (3) と同じ式、ε 0.15): step 500 で、島 (壁の x/r_t ∈ [−6, −3.25] の列、壁から 60〜100 節点目) の
     dt_local(P3)/dt_local(P0) を P3 の状態からの期待の倍率 max((|u|+c)/(½(1+β)|u| + c′), 1) で割った値の中央値が 0.85〜1.15、
     かつ島の |c_P3/c_P0 − 1| の中央値 ≤ 0.01。β = (U_r/c)²、U_r = min(c, max(|u|, ε c))、c′ = ½√((1−β)²|u|² + 4U_r²)
     (setDT_d.cu:255-275 の setDTlocal_precond_scale_d と lowMachPrecond_d.cuh:20-37、c′ の法線速度は |u|)。
  Q4 境界の凍結なし (cg_pp2_check.py の (4) と同じ式): 入口・出口・軸・壁 (角を除く) の ρ の step 0 → 2,500 の変化の中央値 / 1 節点内側の変化の中央値
     r_g が、P3 で P0 の 0.1 倍以上。
  Q5 拘束の行 (rlg_an.py の拘束の検査と同じやり方、§6.59 の R2 の成分): Q0・Q3 の stage.h5 で、P0 で d が全節点 0 の成分 (少なくとも 軸の ρv・壁の ρu・ρv・
     強いピンの等温壁の ρE) の P3 の max|d| が 0。前提 (欠ければ判定不能): stage.h5 がそろう、path が blockDPLUR (ライン無し)、commit の自己照合
     q == b + d、2 本の Q0 がビット一致、名前を挙げた 4 成分が P0 で全節点 0。β < 1 は Q3 の res_0 の状態から節点数を記録。
  Q6 組成と物理性: 両腕の全出力 (step 0〜20,000、500 ごと) の tp_series.jsonl で ρ・P・T が有限で正、どの出力でも Y0 の最大 − 最小 ≤ 2.05e−7
     (起点のばらつき。§6.62 の訂正)。detectNaN のダンプ (res_nan_*.h5) があれば不成立。
入力の欠損・読めない量は判定不能 (合格にしない)。不成立 (False) が 1 つでもあれば「不成立」、それ以外で判定不能 (None) があれば「判定不能」。
節点の並び: (ni, nj) = (4467, 121)、j = 0 が軸・j = nj − 1 が壁・i = 0 が入口 (cold_xcheck と同じ)。Q3・Q4 は forge-25 と同じく j を反転して
  jw = 0 を壁にして扱う (島の壁から 60〜100 節点目は j = 20〜60)。
→ _band_ab/cold_pair/tp_check.json"""
import hashlib
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402

OUTD = HERE / "_band_ab" / "cold_pair"
ARMS = {"P0": "run_0806_tp_p0", "P3": "run_0807_tp_p3"}
QUAL = {"P0": "run_0808_tp_q0", "P3": "run_0809_tp_q3"}
LMP = {"P0": 0, "P3": 3}
INPUT = "_g0c/q32_init.h5"
INPUT_SHA = "6f8fcda6d97482e3"
FORGE_SHA = "75ee4675b0c90f41"
EPS = 0.15
Y0_START = 2.05e-7                      # 起点の Y0 のばらつき (B の実測 2.0453e−07、§6.62 の訂正)
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
STEPS = list(range(0, 20001, 500))
NSTEP = {"arm": 20000, "qual": 1}
CFG_REQUIRED = (r"outStepInterval: 500$", r"cfl: 4\.0, cfl_pseudo: 4\.0,", r"implicitRelax: 0\.7,", r"precondEps: 0\.15,", r"implicitSolvePrecision: 0,",
                r"extraFields: \[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local\]\}$")   # 正規表現 (行単位)


def tri(vals):
    """合否の組み合わせ: False が 1 つでもあれば False、そうでなく None があれば None、全部 True なら True。"""
    vals = list(vals)
    if any(v is False for v in vals):
        return False
    if any(v is None for v in vals):
        return None
    return True


def sha16(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:16]


def expected_dt_ratio(vm, c, eps=EPS):
    """setDTlocal_precond_scale_d (setDT_d.cu:255-275) の倍率 max((|u|+c)/(½(1+β)|u| + c′), 1)。c は 1e-8 で下を押さえる (カーネルと同じ)。"""
    c = np.maximum(c, 1e-8)
    ur = np.minimum(c, np.maximum(vm, eps * c))
    beta = (ur / c) ** 2
    cp = 0.5 * np.sqrt((1.0 - beta) ** 2 * vm ** 2 + 4.0 * ur ** 2)
    return np.maximum((vm + c) / np.maximum(0.5 * (1.0 + beta) * vm + cp, 1e-30), 1.0), beta


def vals(path, keys):
    with h5py.File(path, "r") as h:
        miss = [k for k in keys if "VALUE/" + k not in h]
        if miss:
            raise KeyError(f"{path.parent.name}/{path.name} に {miss} が無い")
        return {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}


def grid():
    with h5py.File(HERE / ARMS["P0"] / "nozzle.h5", "r") as h:
        xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)
    ni, nj, S = XC.mesh_info(len(xy))
    return xy, ni, nj, S


def cfg_val(run, key):
    """forge_run.log の実効値 (cg_pp2_check.py の cfg_val と同じ正規表現)。複数あって食い違えば例外。"""
    t = (HERE / run / "forge_run.log").read_text(errors="replace")
    m = sorted(set(re.findall(r"'%s' in '[^']*': (\S+)" % key, t)))
    if len(m) > 1:
        raise ValueError(f"{run}: {key} の実効値が複数 {m}")
    return m[0] if m else None


def prov_sha(run):
    m = re.search(r"forge_sha256\s*:\s*([0-9a-f]+)", (HERE / run / "RUN_PROVENANCE.txt").read_text())
    return m.group(1) if m else None


def node_kinds(ni, nj):
    """節点の種類 (記録用): 壁 j=nj−1、軸 j=0、入口 i=0・出口 i=ni−1 (軸・壁の行を除く)、内部。"""
    k = np.full((ni, nj), "内部", dtype=object)
    k[0, :] = "入口"; k[-1, :] = "出口"; k[:, 0] = "軸"; k[:, -1] = "壁"
    return k.ravel()


def one_line_diff(a, b):
    """2 つの設定の行の差: 行数が同じで違う行の添字の一覧 (行数が違えば None)。"""
    la, lb = a.splitlines(), b.splitlines()
    if len(la) != len(lb):
        return None, la, lb
    return [i for i in range(len(la)) if la[i] != lb[i]], la, lb


def q1(out):
    r = {}
    runs = list(ARMS.values()) + list(QUAL.values())
    # 入力の sha
    try:
        shas = {INPUT: sha16(HERE / INPUT)} | {f"{x}/nozzle.h5": sha16(HERE / x / "nozzle.h5") for x in runs}
        r["input_sha"] = {"values": shas, "ok": all(v == INPUT_SHA for v in shas.values())}
    except OSError as e:
        r["input_sha"] = {"error": str(e), "ok": None}
    # 4 本の res_0 の 9 保存量のビット一致 (P0 の腕を基準)
    try:
        base = vals(HERE / ARMS["P0"] / "res_0.h5", Q9)
        diff = {}
        for x in runs[1:]:
            v = vals(HERE / x / "res_0.h5", Q9)
            bad = {q: (int(np.count_nonzero(v[q] != base[q])) if v[q].shape == base[q].shape else f"形が違う {v[q].shape}")
                   for q in Q9 if v[q].shape != base[q].shape or not np.array_equal(v[q], base[q])}
            diff[x] = bad
        r["res0_bitwise"] = {"not_equal_vs_P0": diff, "ok": all(not b for b in diff.values())}
        # 入力との差 (記録だけ): 節点の種類ごとの不一致の数
        inp = vals(HERE / INPUT, Q9)
        _, ni, nj, _ = grid()
        kinds = node_kinds(ni, nj)
        rec = {}
        for q in Q9:
            ne = inp[q] != base[q]
            rec[q] = {k: int(np.count_nonzero(ne & (kinds == k))) for k in ("壁", "入口", "出口", "軸", "内部")}
        r["res0_vs_input_record"] = rec
    except (OSError, KeyError, SystemExit) as e:
        r["res0_bitwise"] = {"error": str(e), "ok": None}
    # 設定の差
    try:
        txt = {nm: (HERE / x / "solverConfig.yaml").read_text() for nm, x in (("P0", ARMS["P0"]), ("P3", ARMS["P3"]), ("Q0", QUAL["P0"]), ("Q3", QUAL["P3"]))}
        cfg = {}
        idx, la, lb = one_line_diff(txt["P0"], txt["P3"])
        ok_arm = (idx is not None and len(idx) == 1 and "lowMachPrecond: 0, " in la[idx[0]] and "lowMachPrecond: 3, " in lb[idx[0]]
                  and la[idx[0]].replace("lowMachPrecond: 0, ", "X") == lb[idx[0]].replace("lowMachPrecond: 3, ", "X"))
        cfg["P0_vs_P3"] = {"lines": [la[i] for i in idx] + [lb[i] for i in idx] if idx is not None else "行数が違う", "ok": bool(ok_arm)}
        for qn, an in (("Q0", "P0"), ("Q3", "P3")):
            idx, la, lb = one_line_diff(txt[an], txt[qn])
            ok = (idx is not None and len(idx) == 1 and f"nStepOuter: {NSTEP['arm']}}}" in la[idx[0]] and f"nStepOuter: {NSTEP['qual']}}}" in lb[idx[0]]
                  and la[idx[0]].replace(f"nStepOuter: {NSTEP['arm']}}}", "X") == lb[idx[0]].replace(f"nStepOuter: {NSTEP['qual']}}}", "X"))
            cfg[f"{an}_vs_{qn}"] = {"lines": [la[i] for i in idx] + [lb[i] for i in idx] if idx is not None else "行数が違う", "ok": bool(ok)}
        for nm in ("P0", "P3"):
            t = txt[nm]
            miss = [s for s in CFG_REQUIRED if not re.search(s, t, re.M)] + ([] if f"nStepOuter: {NSTEP['arm']}}}" in t else ["nStepOuter"])
            cnt = t.count("lowMachPrecond")
            cfg[f"{nm}_registered_values"] = {"missing": miss, "lowMachPrecond_keys": cnt, "ok": not miss and cnt == 1}
        ys = {x: sorted(p.name for p in (HERE / x).glob("*.yaml") if p.name != "solverConfig.yaml") for x in runs}
        same_set = len({tuple(v) for v in ys.values()}) == 1
        bad = [f"{x}/{f}" for x in runs[1:] for f in ys[runs[0]] if (HERE / x / f).exists() and (HERE / x / f).read_bytes() != (HERE / runs[0] / f).read_bytes()]
        cfg["other_yaml"] = {"files": ys[runs[0]], "same_names": same_set, "differ": bad, "ok": same_set and not bad and bool(ys[runs[0]])}
        r["config"] = cfg | {"ok": all(v["ok"] for v in cfg.values())}
    except OSError as e:
        r["config"] = {"error": str(e), "ok": None}
    # forge の sha
    try:
        ps = {x: prov_sha(x) for x in runs}
        r["forge_sha"] = {"values": ps, "ok": None if any(v is None for v in ps.values()) else all(v.startswith(FORGE_SHA) for v in ps.values())}
    except OSError as e:
        r["forge_sha"] = {"error": str(e), "ok": None}
    r["pass"] = tri(r[k]["ok"] for k in ("input_sha", "res0_bitwise", "config", "forge_sha"))
    out["Q1"] = r


def q2(out):
    r = {"runs": {}}
    oks = []
    for nm, x in (("P0", ARMS["P0"]), ("P3", ARMS["P3"]), ("Q0", QUAL["P0"]), ("Q3", QUAL["P3"])):
        lm = LMP["P0" if nm in ("P0", "Q0") else "P3"]
        try:
            v = {k: cfg_val(x, k) for k in ("lowMachPrecond", "precondEps", "implicitSolvePrecision", "speciesImplicitCoupling", "speciesPrecondDt",
                                           "cfl", "cfl_pseudo", "implicitRelax", "nStepInner", "blockDPLUR")}
        except (OSError, ValueError) as e:
            r["runs"][nm] = {"error": str(e), "ok": None}; oks.append(None); continue
        req = (v["lowMachPrecond"], v["precondEps"], v["implicitSolvePrecision"])
        if None in req:
            ok = None
        else:
            try:
                ok = int(v["lowMachPrecond"]) == lm and abs(float(v["precondEps"]) - 0.15) < 1e-12 and int(v["implicitSolvePrecision"]) == 0
            except ValueError:
                ok = False
        rec = dict(v)
        for k, dflt in (("speciesImplicitCoupling", 0), ("speciesPrecondDt", 1)):
            if v[k] is None:
                rec[k] = f"未記載 (yaml に無い → ソースの既定 {dflt}、solverConfig.cpp:394・396)"
        r["runs"][nm] = rec | {"expected_lowMachPrecond": lm, "ok": ok}
        oks.append(ok)
    r["pass"] = tri(oks)
    out["Q2"] = r


def q3(out):
    r = {}
    try:
        xy, ni, nj, S = grid()
        keys = ["dt_local", "Ux", "Uy", "sonic"]
        with h5py.File(HERE / ARMS["P3"] / "res_500.h5", "r") as h:
            has_uz = "VALUE/Uz" in h
        gp = vals(HERE / ARMS["P3"] / "res_500.h5", keys + (["Uz"] if has_uz else []))
        gc = vals(HERE / ARMS["P0"] / "res_500.h5", ["dt_local", "sonic"])
    except (OSError, KeyError, SystemExit) as e:
        out["Q3"] = {"error": str(e), "pass": None}; return
    rs = lambda a: a.reshape(ni, nj)[:, ::-1]  # noqa: E731  jw = 0 が壁 (forge-25 と同じ向き)
    X = rs(xy[:, 0] / S)
    isl = np.zeros((ni, nj), bool)
    cols = (X[:, 0] >= -6.0) & (X[:, 0] <= -3.25)          # 壁の節点の x で列を選ぶ (cg_pp2_check.py と同じ)
    isl[np.ix_(cols, np.arange(60, 101))] = True            # 壁から 60〜100 節点目 (jw)
    vm = np.sqrt(rs(gp["Ux"]) ** 2 + rs(gp["Uy"]) ** 2 + (rs(gp["Uz"]) ** 2 if has_uz else 0.0))
    cP3 = rs(gp["sonic"]); cP0 = rs(gc["sonic"])
    expect, beta = expected_dt_ratio(vm, cP3)
    meas = rs(gp["dt_local"]) / rs(gc["dt_local"])
    q = meas[isl] / expect[isl]
    close = np.abs(cP3[isl] / cP0[isl] - 1.0)
    fin = np.isfinite(q) & np.isfinite(close)
    r["island"] = {"n_nodes": int(isl.sum()), "n_cols": int(cols.sum()), "nonfinite": int(np.count_nonzero(~fin))}
    if isl.sum() == 0 or not fin.all():
        r["pass"] = None; r["reason"] = "島が空か、島に非有限がある"; out["Q3"] = r; return
    r["median_measured_ratio"] = float(np.median(meas[isl])); r["median_expected_ratio"] = float(np.median(expect[isl]))
    r["median_measured_over_expected"] = float(np.median(q)); r["p10_p90_measured_over_expected"] = [float(np.percentile(q, 10)), float(np.percentile(q, 90))]
    r["median_abs_sonic_rel_diff"] = float(np.median(close))
    r["median_M"] = float(np.median(vm[isl] / np.maximum(cP3[isl], 1e-8))); r["median_beta"] = float(np.median(beta[isl]))
    r["pass"] = bool(0.85 <= r["median_measured_over_expected"] <= 1.15 and r["median_abs_sonic_rel_diff"] <= 0.01)
    out["Q3"] = r


def q4(out):
    r = {}
    try:
        _, ni, nj, _ = grid()

        def rgroups(run):
            a0 = vals(HERE / run / "res_0.h5", ["ro"])["ro"]; a1 = vals(HERE / run / "res_2500.h5", ["ro"])["ro"]
            d = np.abs(a1 - a0).reshape(ni, nj)[:, ::-1]                    # jw = 0 が壁
            g = {"inlet": (d[0, 1:-1], d[1, 1:-1]), "outlet": (d[-1, 1:-1], d[-2, 1:-1]), "axis": (d[1:-1, -1], d[1:-1, -2]), "wall": (d[1:-1, 0], d[1:-1, 1])}
            res = {}
            for k, (b, n) in g.items():
                if not (np.isfinite(b).all() and np.isfinite(n).all()):
                    res[k] = None
                else:
                    res[k] = float(np.median(b) / np.median(n)) if np.median(n) > 0 else None
            return res
        rc, rp = rgroups(ARMS["P0"]), rgroups(ARMS["P3"])
    except (OSError, KeyError, SystemExit) as e:
        out["Q4"] = {"error": str(e), "pass": None}; return
    r["r_g"] = {"P0": rc, "P3": rp}
    r["per_group"] = {k: (None if rc[k] is None or rp[k] is None else bool(rp[k] >= 0.1 * rc[k])) for k in rc}
    r["pass"] = tri(r["per_group"].values())
    out["Q4"] = r


QS = ("ro", "roUx", "roUy", "roUz", "roe")


def load_stage(run):
    p = HERE / run / "stage.h5"
    if not p.exists():
        return None
    with h5py.File(p, "r") as h:
        d = {g: {q: np.asarray(h[f"/{g}/{q}"][:], np.float64) for q in QS} for g in ("Q0", "R", "b", "d", "q")}
        fl = {k: np.asarray(h[f"/cells/{k}_flag"][:]).astype(bool) for k in ("bnode", "axis", "wall", "iso_wall") if f"/cells/{k}_flag" in h}
        at = {k: (h.attrs[k].decode() if isinstance(h.attrs[k], bytes) else (np.asarray(h.attrs[k]).item() if np.asarray(h.attrs[k]).size == 1 else str(h.attrs[k])))
              for k in h.attrs}
    return d, fl, at


def q5(out):
    r = {"pre": {}, "named": {}, "auto": {}}
    und = []
    try:
        S0, S3 = load_stage(QUAL["P0"]), load_stage(QUAL["P3"])
    except (OSError, KeyError) as e:
        out["Q5"] = {"error": str(e), "pass": None}; return
    if S0 is None or S3 is None:
        out["Q5"] = {"error": "stage.h5 がそろわない", "pass": None}; return
    for nm, (d, fl, at) in (("Q0", S0), ("Q3", S3)):
        r["pre"][f"{nm} 属性"] = {k: str(at[k])[:80] for k in ("path", "commit_kernel", "d_source", "step") if k in at}
        path = str(at.get("path", ""))
        if "blockDPLUR" not in path or "line" in path.lower():
            und.append(f"{nm}: path が blockDPLUR (ライン無し) でない ({path})")
        mism = [q for q in QS if not np.array_equal(d["q"][q], (d["b"][q] + d["d"][q]).astype(np.float64))]
        if mism:
            und.append(f"{nm}: commit の自己照合 q != b + d ({mism})")
        if not all(k in fl for k in ("bnode", "axis", "wall", "iso_wall")):
            und.append(f"{nm}: 節点の印がそろわない ({sorted(fl)})")
    if any(not np.array_equal(S3[0]["Q0"][q], S0[0]["Q0"][q]) for q in QS):
        und.append("2 本の Q0 がビット一致しない")
    if und:
        out["Q5"] = r | {"undecidable": und, "pass": None}; return
    d0, fl, _ = S0; d3 = S3[0]
    n = len(d0["d"]["ro"])
    ax, wa, iw, bn = (fl[k] for k in ("axis", "wall", "iso_wall", "bnode"))
    # 名前を挙げた拘束の成分 (§6.59 の R2・port_reg_an.py と同じ)
    cons = {"軸|roUy": (ax, "roUy"), "壁|roUx": (wa & ~ax, "roUx"), "壁|roUy": (wa & ~ax, "roUy"), "等温壁|roe": (iw & ~ax, "roe")}
    oks = []
    for k, (m, q) in cons.items():
        if not m.any():
            r["named"][k] = {"節点数": 0}; und.append(f"{k}: 該当の節点が無い"); oks.append(None); continue
        p0 = float(np.max(np.abs(d0["d"][q][m]))); p3 = float(np.max(np.abs(d3["d"][q][m])))
        r["named"][k] = {"節点数": int(m.sum()), "P0 の max|d|": p0, "P3 の max|d|": p3}
        if p0 != 0.0:
            und.append(f"{k}: P0 で d が 0 でない (拘束の前提が成り立たない)"); oks.append(None); continue
        oks.append(p3 == 0.0)
    # rlg_an.py と同じ: 領域ごとに P0 で d が全節点 0 の成分を拘束とみなし、P3 で 0 か
    REG = {"軸": ax, "壁": wa & ~ax, "その他の境界": bn & ~ax & ~wa, "内部": ~bn & ~ax & ~wa}
    comps = [q for q in QS if not (q == "roUz" and np.all(d0["d"][q] == 0) and np.all(d3["d"][q] == 0))]
    for nm, m in REG.items():
        for q in comps:
            if not m.any():
                continue
            if np.sqrt(np.sum(d0["d"][q][m] ** 2)) == 0.0:
                v = float(np.max(np.abs(d3["d"][q][m])))
                r["auto"][f"{nm}|{q}"] = {"P3 の max|d|": v}
                oks.append(v == 0.0)
    # β < 1 の記録 (Q3 の res_0 の状態、ε 0.15)
    try:
        st = vals(HERE / QUAL["P3"] / "res_0.h5", ["Ux", "Uy", "sonic"])
        if len(st["Ux"]) == n:
            _, beta = expected_dt_ratio(np.hypot(st["Ux"], st["Uy"]), st["sonic"])
            r["beta_record"] = {"β<1 の節点数": int(np.count_nonzero(beta < 1)), "全節点数": n,
                                **{f"{k} の β の最大": float(np.max(beta[m])) for k, (m, _) in cons.items() if m.any()}}
        else:
            r["beta_record"] = f"res_0 と stage.h5 の節点数が違う ({len(st['Ux'])} 対 {n})"
    except (OSError, KeyError) as e:
        r["beta_record"] = f"読めない ({e})"
    if und:
        r["undecidable"] = und
    r["pass"] = tri(oks)
    out["Q5"] = r


def q6(out):
    r = {"Y0_start_limit": Y0_START}
    try:
        with h5py.File(HERE / INPUT, "r") as h:
            r["input_ptp_roY0_over_ro (記録)"] = float(np.ptp(np.asarray(h["VALUE/roY0"][:], float) / np.asarray(h["VALUE/ro"][:], float)))
    except (OSError, KeyError) as e:
        r["input_ptp_roY0_over_ro (記録)"] = f"読めない ({e})"
    oks = []
    for nm, x in ARMS.items():
        a = {}
        nan_dumps = sorted(p.name for p in (HERE / x).glob("res_nan_*.h5"))
        p = HERE / x / "tp_series.jsonl"
        if not p.exists():
            r[nm] = {"error": "tp_series.jsonl が無い", "ok": None}; oks.append(None); continue
        rows = {}
        for line in p.read_text().splitlines():
            if line.strip():
                o = json.loads(line); rows[o["step"]] = o
        missing = [s for s in STEPS if s not in rows]
        bad = []
        for s, o in sorted(rows.items()):
            for k in ("ro", "P", "T"):
                if o.get(f"nonfinite_{k}") != 0 or o.get(f"nonpos_{k}") != 0:
                    bad.append(f"step {s}: {k} 非有限 {o.get(f'nonfinite_{k}')}・0 以下 {o.get(f'nonpos_{k}')}")
            y = o.get("Y0_ptp")
            if y is None or not y <= Y0_START:
                bad.append(f"step {s}: Y0 の最大 − 最小 {y}")
        yv = [(o["Y0_ptp"], s) for s, o in rows.items() if o.get("Y0_ptp") is not None]
        a["outputs"] = len(rows); a["missing_steps"] = missing; a["bad"] = bad[:20]; a["n_bad"] = len(bad); a["nan_dumps"] = nan_dumps
        a["Y0_ptp_max_step"] = list(max(yv)) if yv else None
        a["Y0_ptp_20000"] = rows.get(20000, {}).get("Y0_ptp")
        a["min_ro_P_T"] = [min((o.get(f"min_{k}") for o in rows.values() if o.get(f"min_{k}") is not None), default=None) for k in ("ro", "P", "T")]
        a["ok"] = False if (bad or nan_dumps) else (None if missing else True)
        r[nm] = a; oks.append(a["ok"])
    r["pass"] = tri(oks)
    out["Q6"] = r


def main():
    out = {}
    for f in (q1, q2, q3, q4, q5, q6):
        try:
            f(out)
        except Exception as e:  # noqa: BLE001  想定外の入力は判定不能にする (合格にしない)
            out[f.__name__.upper()] = {"error": f"{type(e).__name__}: {e}", "pass": None}
    rc = {}
    for x in list(ARMS.values()) + list(QUAL.values()):
        p = HERE / x / "RUN_RC"
        rc[x] = p.read_text().strip() if p.exists() else None
    out["RUN_RC (記録)"] = rc
    passes = {k: out[k]["pass"] for k in ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")}
    out["passes"] = passes
    out["failed"] = [k for k, v in passes.items() if v is False]
    out["undecidable"] = [k for k, v in passes.items() if v is None]
    out["status"] = "不成立" if out["failed"] else ("判定不能" if out["undecidable"] else "成立")
    OUTD.mkdir(parents=True, exist_ok=True)
    (OUTD / "tp_check.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    for k, v in passes.items():
        print(f"{k}: {'合格' if v else ('不合格' if v is False else '判定不能')}  "
              + json.dumps({kk: vv for kk, vv in out[k].items() if kk not in ("pass", "res0_vs_input_record")}, ensure_ascii=False, default=float)[:600])
    print(f"== 資格 §6.62: {out['status']} (不合格 {out['failed']}、判定不能 {out['undecidable']}) → _band_ab/cold_pair/tp_check.json")
    return out


if __name__ == "__main__":
    main()
