#!/usr/bin/env python3
"""床事象のカウンタの「計測の閉じ方」の小型 A/B と単体試験 (GPU と小さい node のケースが要る)。

plans/active/tooling-sern-te-wake-grid.md §5.1 #2 の事前登録 (codex diagnose 2026-10-08):
2 更新だけ回し、**最後の commit の 1 内部節点の ρE だけ**を変えて直後に終了する (密度・運動量・組成・種 DB・
前処理は固定)。A 枝は床エネルギーより十分上、B 枝は十分下 (差は c_v(T_min)×1 K 以上かつ格納 ρE の 64 ULP 以上)。
合格 = A 枝 0 件、B 枝は最後の更新 (監査行 q_index 2)・指定節点に温度床 1 件、コピー上の EOS も床到達
(T = T_min)、監査の前後で実配列がビット不変。B 枝の見逃し・更新番号のずれ・EOS との分類違い・A 枝の誤検出は不合格。

同じ枠組みで単体試験も回す: C = 下限近傍 (床の 0.5 K 上: 床事象 0・床近傍 1)、
P = 圧力床だけ (保存量を同じ倍率で縮め ρ ≥ roMin のまま p < pMin)、R = 密度床 (ρ < roMin)。

    python3 test_floor_events_ab.py --case <雛形 run ディレクトリ> --work <作業ディレクトリ> [--forge <forge>]

雛形は solverConfig.yaml / bcondConfig.yaml / 格子 h5 などを含む node のディレクトリ (case/46 の 2D SERN を
runner_sern.prepare で作ったもの、など)。forge は run_case.sh 経由で回す (FORGE_BIN で切り替え)。
"""
import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_floor_events as cfe   # noqa: E402

ok_all = True


def chk(name, cond, detail=""):
    global ok_all
    ok_all = ok_all and bool(cond)
    print("  [%s] %-70s %s" % ("OK " if cond else "NG ", name, detail))
    return bool(cond)


def set_output_key(cfg_text):
    """output 節に floorEvents: 1 を入れる (既存の flow 形式の output 節には足し、無ければ 1 行足す)。"""
    m = re.search(r"(?m)^output:\s*\{([^}]*)\}", cfg_text)
    if m:
        inner = m.group(1).strip()
        if "floorEvents" in inner:
            return re.sub(r"floorEvents:\s*\d+", "floorEvents: 1", cfg_text)
        return cfg_text[:m.start()] + "output: {" + (inner + ", " if inner else "") + "floorEvents: 1}" + cfg_text[m.end():]
    if re.search(r"(?m)^output:", cfg_text):
        raise RuntimeError("block 形式の output 節は未対応 (flow 形式にすること)")
    return cfg_text.rstrip("\n") + "\noutput: {floorEvents: 1}\n"


def make_run(template, rd, nsteps=2):
    if os.path.exists(rd):
        shutil.rmtree(rd)
    os.makedirs(rd)
    for f in os.listdir(template):
        p = os.path.join(template, f)
        if os.path.isdir(p) or re.match(r"res_", f) or f in ("floor_events.csv", "forge_run.log", "residual_history.csv") \
                or f.endswith(".msh") or f.endswith(".xmf"):
            continue
        shutil.copy2(p, rd)
    cfg = open(os.path.join(rd, "solverConfig.yaml")).read()
    cfg = re.sub(r"nStepOuter:\s*\d+", f"nStepOuter: {nsteps}", cfg)
    cfg = re.sub(r"(?m)^(\s*)outStepInterval:\s*\d+", rf"\g<1>outStepInterval: {nsteps}", cfg)
    cfg = set_output_key(cfg)
    open(os.path.join(rd, "solverConfig.yaml"), "w").write(cfg)


def run(rd, forge, env_extra):
    env = dict(os.environ)
    env.update(env_extra)
    if forge:
        env["FORGE_BIN"] = forge
    r = subprocess.run([os.path.join(HERE, "run_case.sh"), rd], env=env, capture_output=True, text=True)
    log = open(os.path.join(rd, "forge_run.log")).read() if os.path.exists(os.path.join(rd, "forge_run.log")) else ""
    m = re.search(r"forge exit=(\d+)", r.stdout)
    return (int(m.group(1)) if m else -1), log


def rows_of(rd):
    rows, err = cfe.load(os.path.join(rd, cfe.FILE_NAME))
    return rows, err


def interior_cold_node(mesh_h5):
    """境界 (BCONDS/*/iCells) に入らない節点のうち、冷たい外部流 (初期場の ρ が小さい側) で壁から最も遠いもの。"""
    with h5py.File(mesh_h5, "r") as f:
        bnd = set()
        for k in f["BCONDS"]:
            bnd.update(int(i) for i in f[f"BCONDS/{k}/iCells"][:])
        ro = f["VALUE/ro"][:].astype(float)
        wd = f["VALUE/wall_dist"][:].astype(float) if "VALUE/wall_dist" in f else np.zeros_like(ro)
    n = len(ro)
    cand = np.array([i for i in range(n) if i not in bnd])
    rmed = np.median(ro[cand])
    cold = cand[ro[cand] <= rmed]          # 2D SERN では外部流 (排気より低密度・低温)
    pick = cold[np.argmax(wd[cold])] if len(cold) else cand[np.argmax(wd[cand])]
    return int(pick), float(ro[pick])


def parse_inject(log):
    m = re.search(r"test inject: step (\d+) node (\d+) mode de roe (\S+) -> (\S+) \(e_floor (\S+) J/kg, cv\(T_floor\) (\S+) J/kg/K, "
                  r"rho (\S+), ek (\S+), de (\S+)\)", log)
    if not m:
        return None
    k = ("step", "node", "roe_old", "roe_new", "e_floor", "cv", "rho", "ek", "de")
    return {kk: (int(v) if kk in ("step", "node") else float(v)) for kk, v in zip(k, m.groups())}


def audit_copy_T(log):
    m = re.search(r"test audit copy: T\[(\d+)\] = (\S+) K \(T_floor (\S+) K\)", log)
    return (float(m.group(2)), float(m.group(3))) if m else (None, None)


def audit_verify(log):
    m = re.search(r"audit verify: (\d+) cell arrays compared byte-wise before/after the audit, (\d+) changed", log)
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)


def common_checks(tag, rd, rc, log, node, nsteps=2):
    chk(f"{tag}: forge rc == 0", rc == 0, f"rc={rc}")
    rows, err = rows_of(rd)
    chk(f"{tag}: 記録が読める", not err, "; ".join(err[:2]))
    narr, nchg = audit_verify(log)
    chk(f"{tag}: 監査の前後で実配列 (全 cell 配列) がビット不変", narr is not None and narr > 20 and nchg == 0, f"{narr} 配列中 {nchg} 変化")
    eos = [r for r in rows if r["kind"] == "eos"]
    aud = [r for r in rows if r["kind"] == "audit"]
    end = [r for r in rows if r["kind"] == "session_end"]
    chk(f"{tag}: eos 行が step 1..{nsteps} (q_index 0..{nsteps - 1}) で 1 行ずつ",
        [(r["step"], r["inner"], r["q_index"]) for r in eos] == [(k, 0, k - 1) for k in range(1, nsteps + 1)],
        str([(r["step"], r["q_index"]) for r in eos]))
    chk(f"{tag}: 監査行が 1 行 (step {nsteps}, q_index {nsteps})", len(aud) == 1 and aud[0]["q_index"] == nsteps and aud[0]["step"] == nsteps)
    chk(f"{tag}: session_end に audit=done", len(end) == 1 and end[0]["note"] == "audit=done", end[0]["note"] if end else "")
    clean = all(r["nT_real"] == 0 and r["nRho_real"] == 0 and r["nP_real"] == 0 and r["n_nonfinite_real"] == 0
                and r["nT_mismatch_real"] == 0 and r["nT_uneval_real"] == 0 and r["overflow"] == 0 for r in eos)
    chk(f"{tag}: 書き換え前の更新 (q_index 0..{nsteps - 1}) は床事象 0・食い違い 0", clean)
    return rows, (aud[0] if aud else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--forge", default=None, help="forge バイナリ (既定は run_case.sh の build/forge)")
    ap.add_argument("--node", type=int, default=None)
    ap.add_argument("--margin-K", type=float, default=5.0, help="A/B の床からの距離 [K 相当] (≥ 1 K を要求)")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    mesh = None
    cfg = open(os.path.join(a.case, "solverConfig.yaml")).read()
    m = re.search(r'meshFileName:\s*"?([^",}\s]+)', cfg)
    mesh = os.path.join(a.case, m.group(1))
    node, ro0 = (a.node, None) if a.node is not None else interior_cold_node(mesh)
    print(f"case {a.case}\nnode {node} (初期場 ρ {ro0})")

    # 雛形の EOS 床の設定 (roMin/pMin) は記録の session_begin から読む (既定値を推測しない)
    # --- 予備: 無書き換えの 2 更新 (c_v(T_min) を知るために de 0 で書き換える = 床ちょうど → 使わない; A を先に回して読む)
    res = {}
    # A: 床より十分上。c_v は記録から読むので、最初は de = margin × 2000 J/kg/K 相当で回し、後で余裕を検査する
    guess_cv = 2000.0
    specs = [("A", f"2:{node}:de:{a.margin_K * guess_cv:.6g}"), ("B", f"2:{node}:de:{-a.margin_K * guess_cv:.6g}")]
    for tag, inj in specs:
        rd = os.path.join(a.work, f"ab_{tag}")
        make_run(a.case, rd)
        rc, log = run(rd, a.forge, {"FORGE_FLOOR_TEST_INJECT": inj, "FORGE_FLOOR_AUDIT_VERIFY": "1"})
        res[tag] = (rd, rc, log)
    # 余裕の検査 (A/B 共通): |de| ≥ c_v(T_min)·1 K かつ ρ|de| ≥ 64 ULP(格納 ρE)
    print("=== 計測の閉じ方の小型 A/B ===")
    for tag in ("A", "B"):
        rd, rc, log = res[tag]
        inj = parse_inject(log)
        if not chk(f"{tag}: 書き換えの記録がある", inj is not None):
            continue
        ulp = float(np.spacing(np.float32(abs(inj["roe_new"]))))
        chk(f"{tag}: |de| ≥ c_v(T_min)·1 K", abs(inj["de"]) >= inj["cv"] * 1.0, f"|de| {abs(inj['de']):.4g} J/kg, c_v {inj['cv']:.4g} J/kg/K")
        chk(f"{tag}: ρ|de| ≥ 64 ULP(ρE)", inj["rho"] * abs(inj["de"]) >= 64 * ulp, f"ρ|de| {inj['rho'] * abs(inj['de']):.4g} J/m³, 64 ULP {64 * ulp:.4g}")
    # A
    rd, rc, log = res["A"]
    rows, aud = common_checks("A", rd, rc, log, node)
    tA, tf = audit_copy_T(log)
    if aud:
        chk("A: 最後の更新 (監査行) に床事象 0 (誤検出なし)", aud["nT_real"] == 0 and aud["nRho_real"] == 0 and aud["nP_real"] == 0
            and aud["nT_mismatch_real"] == 0, f"T {aud['nT_real']} ρ {aud['nRho_real']} p {aud['nP_real']} mismatch {aud['nT_mismatch_real']}")
        chk("A: 床近傍にも入らない (コピー上の EOS の T > T_min + 1 K)", tA is not None and tA > tf + 1.0 and str(node) not in aud["ids_near"].split(";"),
            f"T_copy {tA} K")
    jA = cfe.judge(rd, tail=1.0)
    chk("A: check_floor_events (区間 (0, 2]) は PASS", jA["verdict"] == "PASS", jA["verdict"] + " " + "; ".join(jA["reasons"]))
    # B
    rd, rc, log = res["B"]
    rows, aud = common_checks("B", rd, rc, log, node)
    tB, tf = audit_copy_T(log)
    inj = parse_inject(log)
    if aud:
        chk("B: 最後の更新 (監査行 q_index 2)・指定節点に温度床 1 件", aud["nT_real"] == 1 and aud["ids_T"] == str(node),
            f"nT {aud['nT_real']} ids {aud['ids_T']!r}")
        thermo = cfe._note_dict([r for r in rows if r["kind"] == "session_begin"][0]["note"]).get("thermo")
        if thermo == "tp":
            chk("B: 他の床 (密度・圧力) は 0", aud["nRho_real"] == 0 and aud["nP_real"] == 0, f"ρ {aud['nRho_real']} p {aud['nP_real']}")
        else:
            # CPG の温度下限 tMin (既定 1e-4 K) は P = ρ R tMin ≪ pMin なので、温度床に落ちた節点は圧力床にも必ず入る
            chk("B (CPG): 密度床は 0 (圧力床の同時発生は tMin ≪ pMin/(ρR) から必然)", aud["nRho_real"] == 0,
                f"ρ {aud['nRho_real']} p {aud['nP_real']}")
        chk("B: 述語と EOS の最終温度が一致 (食い違い 0)", aud["nT_mismatch_real"] == 0 and aud["nT_uneval_real"] == 0)
        chk("B: コピー上の EOS も床到達 (T_copy = T_min)", tB is not None and tB == tf, f"T_copy {tB} K, T_min {tf} K")
        if inj:
            exp = inj["rho"] * abs(inj["de"])
            chk("B: 補正量 Δ(ρE) = ρ (e_mix(T_min) − e_in) ≈ ρ|de|", abs(aud["dRhoE_T_sum"] - exp) <= 1e-3 * exp + 64 * float(np.spacing(np.float32(abs(inj["roe_new"])))),
                f"Δ(ρE) {aud['dRhoE_T_sum']:.6g} vs ρ|de| {exp:.6g}")
        chk("B: 床近傍 (補助) にも指定節点が入る", str(node) in aud["ids_near"].split(";"), aud["ids_near"])
    jB = cfe.judge(rd, tail=1.0)
    chk("B: check_floor_events (区間 (0, 2]) は FAIL (最初の事象 q_index 2)", jB["verdict"] == "FAIL" and jB["events"]["T"]["first_q_index"] == 2,
        jB["verdict"] + " " + "; ".join(jB["reasons"]))
    jB1 = cfe.judge(rd, window_steps=1)
    chk("B: 区間 (1, 2] でも FAIL (最後の 1 更新だけでも見逃さない)", jB1["verdict"] == "FAIL", jB1["verdict"])

    # 単体: 下限近傍 C、圧力床だけ P、密度床 R
    print("=== 単体: 下限近傍・圧力床・密度床 ===")
    injA = parse_inject(res["A"][2])
    cv = injA["cv"] if injA else guess_cv
    rd = os.path.join(a.work, "unit_C_near")
    make_run(a.case, rd)
    rc, log = run(rd, a.forge, {"FORGE_FLOOR_TEST_INJECT": f"2:{node}:de:{0.5 * cv:.9g}", "FORGE_FLOOR_AUDIT_VERIFY": "1"})
    rows, aud = common_checks("C", rd, rc, log, node)
    tC, tf = audit_copy_T(log)
    if aud:
        chk("C: 床の 0.5 K 上 → 床事象 0", aud["nT_real"] == 0 and aud["nT_mismatch_real"] == 0, f"nT {aud['nT_real']}")
        chk("C: 床近傍 (T ≤ T_min + 1 K) に指定節点", str(node) in aud["ids_near"].split(";") and aud["n_near_real"] >= 1,
            f"near {aud['n_near_real']} ids {aud['ids_near']!r}, T_copy {tC}")
        chk("C: コピー上の EOS の T は (T_min, T_min + 1 K]", tC is not None and tf < tC <= tf + 1.0, f"T_copy {tC}")
    # P: A の res_2 から節点の ρ と p を読み、ρ s ≥ 1.25 roMin かつ p s ≤ pMin / 1.25 となる倍率 s を選ぶ
    notes = cfe._note_dict([r for r in rows_of(res["A"][0])[0] if r["kind"] == "session_begin"][0]["note"])
    roMin, pMin = float(notes["romin"]), float(notes["pmin"])
    with h5py.File(os.path.join(res["A"][0], "res_2.h5"), "r") as f:
        rho_n, p_n = float(f["VALUE/ro"][node]), float(f["VALUE/P"][node])
    lo, hi = 1.25 * roMin / rho_n, pMin / (1.25 * p_n)
    if not (lo < hi):
        # 保存量を同じ倍率で縮める方法では p/ρ = R T が変わらないので、R T ≥ pMin/roMin の節点では圧力床だけを起こせない
        print(f"  [SKIP] P: ρ s ≥ 1.25 roMin かつ p s ≤ pMin/1.25 の倍率が無い (R T = {p_n / rho_n:.4g} ≥ pMin/roMin/1.56) — この雛形では圧力床単独の試験は不可")
    else:
        s = (lo * hi) ** 0.5
        rd = os.path.join(a.work, "unit_P")
        make_run(a.case, rd)
        rc, log = run(rd, a.forge, {"FORGE_FLOOR_TEST_INJECT": f"2:{node}:rhoscale:{s:.9g}", "FORGE_FLOOR_AUDIT_VERIFY": "1"})
        rows, aud = common_checks("P", rd, rc, log, node)
        if aud:
            chk("P: 圧力床 1 件 (指定節点)、温度床・密度床は 0", aud["nP_real"] == 1 and aud["ids_P"] == str(node)
                and aud["nT_real"] == 0 and aud["nRho_real"] == 0, f"p {aud['nP_real']} T {aud['nT_real']} ρ {aud['nRho_real']}")
            exp = pMin - s * p_n
            chk("P: ΔP = pMin − P_raw ≈ pMin − s p (1 更新分の差を許す 5 %)", abs(aud["dP_sum"] - exp) <= 0.05 * pMin, f"ΔP {aud['dP_sum']:.4g} vs {exp:.4g}")
    rd = os.path.join(a.work, "unit_R")
    s = 0.5 * roMin / rho_n
    make_run(a.case, rd)
    rc, log = run(rd, a.forge, {"FORGE_FLOOR_TEST_INJECT": f"2:{node}:rhoscale:{s:.9g}", "FORGE_FLOOR_AUDIT_VERIFY": "1"})
    rows, aud = common_checks("R", rd, rc, log, node)
    if aud:
        m = re.search(r"mode rhoscale \S+: ro \S+ -> (\S+),", log)
        ro_new = float(m.group(1)) if m else float("nan")
        chk("R: 密度床 1 件 (指定節点)", aud["nRho_real"] == 1 and aud["ids_Rho"] == str(node), f"ρ {aud['nRho_real']} ids {aud['ids_Rho']!r}")
        chk("R: Δρ = roMin − ρ_in (float の丸めまで)", abs(aud["dRho_sum"] - (roMin - ro_new)) <= 1e-6 * roMin, f"Δρ {aud['dRho_sum']:.6g} vs {roMin - ro_new:.6g}")

    print("\nVERDICT: %s" % ("PASS" if ok_all else "FAIL"))
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
