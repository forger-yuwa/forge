#!/usr/bin/env python3
"""R7b ④-0 の再判定の台本 (case/46.sern_design/diag/rejudge_r7b.py) の試験 (plan tooling-nozzle-sern-chain §5.1 R7b「④ の改訂
(2026-10-13)」の ④-0)。全入力のそろった判定区間の run (GPU の接続試験 run_sern_judgment_connection_test.py --gpu --keep が作る
run_conn = 粗い格子の 3D・判定区間の段 80 step) を複写して (元は触らない) 台本の分岐を確かめる。現行の prepare に変換器が要る。

合格 (測る前に固定):
  (a) 全入力がそろえば判定する: status JUDGED。再判定の metrics.json の gates に現行の受理検査の版と判定区間 (#2・#3) が入り、
      旧 metrics.json との A/B (AB_GATES.md)・config の差なし・格子署名と双対の一致・TE_MONITOR の列と最終行の一致・使った YAML の値
      (判定区間のキーは生産 YAML から、te_wake の YAML に無いキーは「(無い)」) が記録され、現行の prepare の格子は消える
  (b) 入力の欠落 (残差・最終場・壁の部品 1 標本・TE_MONITOR.csv・floor_events.csv) を全部列挙して判定不能 (inputs)。metrics.json を
      作らず、現行の prepare にも進まない。prepare の YAML の名前が対応表と違う run も判定不能 (provenance)
  (c) config の許さない差 (cfl)・格子署名の記録の不一致・TE_MONITOR.csv の列の違いで、判定せずに止まり (config / mesh / te_monitor)、
      差を記録する。metrics.json を作らない
  (d) 出力先が既にあれば失敗 (終了コード 2) し、その run のファイル (既存の再判定の出力を含む) を変えない。どの場合も run の元ファイルの
      mtime と sha256 が台本の前後で同じ (台本の自己申告とは別に試験側で測る)。複写元の run も変わらない
  (e) まとめ: 合成の記録 (対応表の 7 本、全部 PASS) で結果 A、Δᵢ・ΔJ が式 (0.5Δ₆ + 0.3Δ₁₀ + 0.2Δ₄、q = C_T_with_shear) どおり、
      時間幅 u = max(a, d, 1e-6)、G₁₀ は未確定。1 本を判定不能・1 本を未実施にすると結果 B でその run を挙げる。
      まとめの出力先が既にあれば失敗 (終了コード 1)

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_rejudge_r7b_tests.py --src-run <run_conn> \\
      [--src-problem <problem_conn.yaml>] [--out DIR] [--keep]
  --src-run を省くと (e) だけ回す。--gpu を付けると接続試験 (GPU・forge が要る) で run_conn を作ってから回す。
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))

from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402

SCRIPT = ROOT / "case" / "46.sern_design" / "diag" / "rejudge_r7b.py"
CONN_TEST = HERE / "run_sern_judgment_connection_test.py"
TAG = "REJUDGE_TEST"
fails = 0


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def load_script():
    spec = importlib.util.spec_from_file_location("rejudge_r7b", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["rejudge_r7b"] = mod          # dataclass が自分のモジュールを引くので先に登録する
    spec.loader.exec_module(mod)
    return mod


def fingerprint(d: Path, exclude: str | None = None) -> dict:
    """{相対パス: (大きさ, mtime_ns, sha256)}。exclude はトップの除外ディレクトリ名。"""
    out = {}
    for dp, dns, fns in os.walk(d):
        if exclude and Path(dp) == d:
            dns[:] = [x for x in dns if x != exclude]
        for fn in fns:
            f = Path(dp) / fn
            st = f.stat()
            h = hashlib.sha256(f.read_bytes()).hexdigest()
            out[str(f.relative_to(d))] = (st.st_size, st.st_mtime_ns, h)
    return out


def run_script(*args):
    r = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def info_of(rd: Path) -> dict:
    p = rd / TAG / "REJUDGE_INFO.json"
    return json.loads(p.read_text()) if p.is_file() else {}


ap = argparse.ArgumentParser()
ap.add_argument("--src-run", default=None, help="全入力のそろった判定区間の run (接続試験の run_conn)。省くと (e) だけ")
ap.add_argument("--src-problem", default=None, help="その run の問題 YAML (既定 <src-run>/../problem_conn.yaml)")
ap.add_argument("--gpu", action="store_true", help="接続試験 (GPU・forge が要る) で run_conn を作ってから回す")
ap.add_argument("--out", default=None)
ap.add_argument("--keep", action="store_true")
a = ap.parse_args()
RJ = load_script()
S = Path(a.out) if a.out else Path(tempfile.mkdtemp(prefix="sern_rejudge_r7b_"))
S.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------------------------------- (e) まとめ (合成の記録)
MEANS = {("m10_on", 0.8): 0.8969471, ("m10_on", 1.0): 0.8969681, ("m6_on", 0.8): 0.8704123, ("m6_on", 1.0): 0.8725505,
         ("m4_off", 0.8): 0.9676994, ("m4_off", 1.0): 0.9676520}
AD = {0.8: (2.0e-6, 3.0e-6), 1.0: (5.0e-7, 4.0e-7)}      # u = 3e-6 (0.8)・1e-6 (1.0、下限)


def fake_rec(spec, verdict="PASS", status="JUDGED", g4_mean=0.8966):
    m = g4_mean if spec.grid == "g4" else MEANS[(spec.op, spec.L_sw)]
    av, dv = AD[spec.L_sw]
    fw = {q: {"mean": m + 0.01 * k, "a": av, "d": dv, "tol": 5e-5, "steady": "STEADY", "ok": True}
          for k, q in enumerate(("C_T_with_shear", "C_T", "C_L", "C_M"))}
    fw["C_T_with_shear"]["mean"] = m
    s = {"status": status, "verdict": verdict if status == "JUDGED" else None, "fail_class": None if verdict == "PASS" else "UNSTEADY",
         "acceptance_version": G.ACCEPTANCE_VERSION, "reasons": [],
         "judgment": {"verdict": "PASS", "fail_class": None, "temperature": "PASS", "force_window": "PASS", "floor_events": "PASS",
                      "floor_window": [10000, 20000]},
         "force_window": fw, "residual_evidence_input_ok": True, "missing_inputs": [],
         "stop": None if status != "INDETERMINATE" else {"kind": "inputs", "reasons": ["入力が欠けている: floor_events.csv"]},
         "run_files_unchanged": True, "error": None}
    return {"status": status, "summary": s, "out_dir": f"/fake/{spec.run}/{TAG}", "end": "2026-10-13", "git": {"head": "x"}}


recs = {s.run: fake_rec(s) for s in RJ.RUNS}
sm = RJ.summarize(RJ.RUNS, recs)
dJ_exp = (0.5 * (MEANS[("m6_on", 1.0)] - MEANS[("m6_on", 0.8)]) + 0.3 * (MEANS[("m10_on", 1.0)] - MEANS[("m10_on", 0.8)])
          + 0.2 * (MEANS[("m4_off", 1.0)] - MEANS[("m4_off", 0.8)]))
dd = sm["design_diff"]
check("(e) 合成の記録 7 本が全部 PASS → 結果 A", sm["result"] == "A" and not sm["failures"], str(sm["failures"])[:300])
check("(e) ΔJ = 0.5Δ₆ + 0.3Δ₁₀ + 0.2Δ₄ (C_T_with_shear) が式どおり (1.0659e-3)",
      dd["dJ"] is not None and abs(dd["dJ"] - dJ_exp) < 1e-12 and abs(dJ_exp - 1.0659e-3) < 1e-7, f"{dd['dJ']} / {dJ_exp}")
check("(e) 時間幅 U = u₀.₈ + u₁.₀、u = max(a, d, 1e-6) (3e-6 + 1e-6)",
      all(abs(v["quantities"]["C_T_with_shear"]["U"] - 4.0e-6) < 1e-15 for v in dd["per_op"].values())
      and abs(dd["U_time_weighted"] - 4.0e-6) < 1e-15, str({k: v["quantities"]["C_T_with_shear"]["U"] for k, v in dd["per_op"].items()}))
g10 = dd["G10"] or {}
ref = g10.get("reference_0.8") or {}
check("(e) G₁₀ は未確定 (g4 は L_sw 0.8 だけ)、参考に L_sw 0.8 の q(g4) − q(g3) を出す",
      g10.get("status") == "未確定" and isinstance(ref, dict) and abs(ref.get("q_g4_minus_g3", 1) - (0.8966 - 0.8969471)) < 1e-12, str(g10))
recs_b = dict(recs)
recs_b["run_1090_tewake_m4_L10"] = fake_rec(next(s for s in RJ.RUNS if s.run == "run_1090_tewake_m4_L10"), status="INDETERMINATE")
recs_b.pop("run_1086_tewake_g4_L10_cont40k")
sm_b = RJ.summarize(RJ.RUNS, recs_b)
fb = {f["run"]: f["reasons"] for f in sm_b["failures"]}
check("(e) 1 本を判定不能・1 本を未実施 → 結果 B でその 2 本を挙げる、ΔJ は未算出",
      sm_b["result"] == "B" and set(fb) == {"run_1090_tewake_m4_L10", "run_1086_tewake_g4_L10_cont40k"}
      and "判定不能" in fb["run_1090_tewake_m4_L10"][0] and "未実施" in fb["run_1086_tewake_g4_L10_cont40k"][0] and sm_b["design_diff"]["dJ"] is None,
      str(fb)[:300])
recs_c = dict(recs)
recs_c["run_1084_tewake_g3_L10_m10"] = fake_rec(next(s for s in RJ.RUNS if s.run == "run_1084_tewake_g3_L10_m10"), verdict="FAIL")
sm_c = RJ.summarize(RJ.RUNS, recs_c)
check("(e) 現行で FAIL の行は結果 B に入り、その作動点の Δ は出さない",
      sm_c["result"] == "B" and "m10_on" not in sm_c["design_diff"]["per_op"] and sm_c["design_diff"]["dJ"] is None)
md, js = S / "sum_e.md", S / "sum_e.json"
RJ.write_summary(sm, md, js)
check("(e) まとめの .md・.json を書く (結果 A・ΔJ の表示)", "結果 A" in md.read_text() and "ΔJ" in md.read_text()
      and json.loads(js.read_text())["result"] == "A")
rc_e, out_e = run_script("--summary-only", "--case-dir", S, "--summary-out", S / "sum_e", "--spec", "x,0.8,m6_on,g3,a.yaml,b.yaml")
check("(e) まとめの出力先が既にあれば失敗 (終了コード 1) し、上書きしない", rc_e == 1 and "既にある" in out_e
      and json.loads(js.read_text())["result"] == "A", out_e[-300:])

# ----------------------------------------------------------------------------------------------- (a)〜(d) 実データの複写
src = Path(a.src_run).resolve() if a.src_run else None
if a.gpu and src is None:
    if not R2.CONVERTER.exists():
        print(f"SKIP (a)〜(d): 変換器 {R2.CONVERTER} が無い")
    else:
        env = {**os.environ, "FORGE_CUDA_BLOCKSIZE": "128", "FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}
        r = subprocess.run([sys.executable, str(CONN_TEST), "--out", str(S / "conn"), "--keep"], env=env, capture_output=True, text=True)
        print(f"接続試験 rc {r.returncode}: {(r.stdout or '').strip().splitlines()[-1:]}")
        src = S / "conn" / "run_conn"
if src is None or not src.is_dir():
    print("SKIP (a)〜(d): --src-run (全入力のそろった判定区間の run、接続試験 run_sern_judgment_connection_test.py --gpu --keep の run_conn) "
          "が無い。--gpu で作れる")
elif not R2.CONVERTER.exists():
    print(f"SKIP (a)〜(d): 変換器 {R2.CONVERTER} が無い (現行の prepare を回せない)")
else:
    src_prob = Path(a.src_problem).resolve() if a.src_problem else src.parent / "problem_conn.yaml"
    info0 = json.loads((src / "prepare_info.json").read_text())
    op = info0["operating_point"]["name"]
    prob_name = Path(info0["problem"]).name             # run の prepare_info.json の problem と同じ名前の te_wake の YAML を置く
    y = yaml.safe_load(src_prob.read_text())
    L_sw = float(y["mesh3d"]["L_sw"])
    N = int(y["evaluate"]["judgment_steps"])
    src_fp0 = fingerprint(src)
    C = S / "case"
    C.mkdir(exist_ok=True)
    (C / "prod.yaml").write_text(src_prob.read_text())          # 生産 YAML の役 (判定区間のキーを持つ)
    yt = yaml.safe_load(src_prob.read_text())
    for k in ("judgment_steps", "judgment_extend_once"):        # te_wake の YAML の役: 判定区間の段のキーが無い (実物と同じ形)
        yt["evaluate"].pop(k, None)
    (C / prob_name).write_text(yaml.safe_dump(yt, allow_unicode=True, sort_keys=False))
    (C / "problem_other_name.yaml").write_text((C / prob_name).read_text())

    def clone(name):
        d = C / name
        shutil.copytree(src, d, symlinks=False)
        return d

    def spec(name, te=prob_name):
        return f"{name},{L_sw},{op},g3,{te},prod.yaml"

    def go(name, te=prob_name, extra=()):
        rd = C / name
        fp0 = fingerprint(rd, exclude=TAG)
        rc, out = run_script("--case-dir", C, "--tag", TAG, "--no-summary", "--spec", spec(name, te), *extra)
        fp1 = fingerprint(rd, exclude=TAG)
        return rc, out, info_of(rd), fp0 == fp1

    # (a) 全入力
    clone("run_ok")
    rc, out, inf, same = go("run_ok")
    o = C / "run_ok" / TAG
    m = json.loads((o / "metrics.json").read_text()) if (o / "metrics.json").is_file() else {}
    g = m.get("gates") or {}
    check("(a) 全入力がそろえば判定する (status JUDGED、終了コード 0)", rc == 0 and inf.get("status") == "JUDGED", out[-600:])
    check("(a) 再判定の metrics.json の gates に現行の受理検査の版と判定区間 (#2・#3) が入る",
          g.get(G.ACCEPTANCE_VERSION_KEY) == G.ACCEPTANCE_VERSION and isinstance(g.get("judgment"), dict)
          and inf.get("summary", {}).get("verdict") == g.get("verdict"), str(g.get("verdict")))
    ab = (o / "AB_GATES.md").read_text() if (o / "AB_GATES.md").is_file() else ""
    check("(a) A/B (旧 metrics.json と現行の gates を項目ごと) を書く", "acceptance_version" in ab and "#3 力の窓" in ab
          and "#1 床事象" in ab and len(inf.get("ab", {}).get("rows", [])) >= 15)
    check("(a) config の差なし (solverConfig.yaml・solverConfig_main.yaml・bcond・化学種)、CONFIG_DIFF.txt を残す",
          inf.get("config", {}).get("disallowed") == [] and "solverConfig.yaml" in inf.get("config", {}).get("compared", {})
          and (o / "CONFIG_DIFF.txt").is_file(), str(inf.get("config"))[:300])
    mi = (inf.get("mesh") or {}).get("items") or {}
    check("(a) 格子署名・双対が run の格子 (いま再計算)・run の記録・現行の prepare で一致",
          (inf.get("mesh") or {}).get("disallowed") == [] and all(mi.get(k, {}).get("current") == mi.get(k, {}).get("run_now")
                                                                  == mi.get(k, {}).get("run_record") is not None
                                                                  for k in ("mesh_signature", "dual_hash")), str(mi)[:300])
    tm = inf.get("te_monitor") or {}
    check("(a) TE_MONITOR.csv の列と最終行が現行の te_monitor.py と一致", (tm.get("columns") or {}).get("same") is True
          and tm.get("final_row_same") is True, str(tm)[:300])
    kp = (inf.get("yaml") or {}).get("keys_from_production") or {}
    check("(a) 使った YAML の値: 判定区間のキーは生産 YAML から (te_wake の YAML に無いキーは「(無い)」)",
          kp.get("judgment_steps", {}).get("used") == N and kp.get("judgment_steps", {}).get("te_wake") == "(無い)"
          and set(kp) == set(RJ.JUDGMENT_KEYS) and (inf.get("yaml") or {}).get("te_wake", {}).get("sha256"), str(kp)[:300])
    prep = o / "prepare_current"
    check("(a) 現行の prepare の格子・中間 (h5・msh・xmf) は消し、設定ファイルは残す",
          prep.is_dir() and not any(f.suffix in (".h5", ".msh", ".xmf") for f in prep.iterdir())
          and (prep / "solverConfig_main.yaml").is_file() and inf.get("current_prepare", {}).get("removed"))
    check("(a)(d) run の元ファイルの mtime・sha256 が台本の前後で同じ (試験側の測定と台本の記録)",
          same and inf.get("run_files_unchanged") is True)

    # (b) 入力の欠落
    rd = clone("run_missing")
    pid = R3._FORCE_PARTS["ramp"][0]
    gone = ["residual_history.csv", f"res_{N}.h5", f"res_ramp_{pid}_{N // 2}.h5", "TE_MONITOR.csv", "floor_events.csv"]
    for f in gone:
        (rd / f).unlink()
    rc, out, inf, same = go("run_missing")
    miss = " / ".join((inf.get("inputs") or {}).get("missing") or [])
    check("(b) 欠けた入力 (残差・最終場・壁の部品 1 標本・TE_MONITOR.csv・floor_events.csv) を全部列挙して判定不能",
          rc == 0 and inf.get("status") == "INDETERMINATE" and (inf.get("stop") or {}).get("kind") == "inputs"
          and all(x in miss for x in ("residual_history.csv", f"res_{N}.h5", "res_ramp_", f"step {N // 2}", "TE_MONITOR.csv", "floor_events.csv")),
          miss[:400])
    check("(b) 判定不能なら metrics.json を作らず、現行の prepare にも進まない",
          not (rd / TAG / "metrics.json").exists() and not (rd / TAG / "prepare_current").exists())
    check("(b)(d) run の元ファイルが台本の前後で同じ", same and inf.get("run_files_unchanged") is True)
    clone("run_prov")
    rc, out, inf, same = go("run_prov", te="problem_other_name.yaml")
    check("(b) prepare の YAML の名前が対応表と違えば判定不能 (provenance)",
          inf.get("status") == "INDETERMINATE" and (inf.get("stop") or {}).get("kind") == "provenance"
          and not (C / "run_prov" / TAG / "metrics.json").exists() and same, str(inf.get("stop"))[:300])

    # (c) 許さない差
    rd = clone("run_cfg")
    t = (rd / "solverConfig.yaml").read_text()
    t2 = t.replace("cfl: 0.25, cfl_pseudo: 0.25", "cfl: 0.3, cfl_pseudo: 0.3")
    assert t2 != t, "試験の前提: solverConfig.yaml に cfl: 0.25, cfl_pseudo: 0.25 がある"
    (rd / "solverConfig.yaml").write_text(t2)
    rc, out, inf, same = go("run_cfg")
    cd = (rd / TAG / "CONFIG_DIFF.txt").read_text() if (rd / TAG / "CONFIG_DIFF.txt").is_file() else ""
    check("(c) config の許さない差 (cfl) で判定せずに止まり、差を記録する",
          inf.get("status") == "INDETERMINATE" and (inf.get("stop") or {}).get("kind") == "config" and "cfl: 0.3" in cd
          and not (rd / TAG / "metrics.json").exists(), str(inf.get("stop"))[:300])
    check("(c) 止めても現行の prepare の格子は消す、run の元ファイルは同じ",
          not any(f.suffix in (".h5", ".msh") for f in (rd / TAG / "prepare_current").iterdir()) and same)
    rd = clone("run_mesh")
    pi = json.loads((rd / "prepare_info.json").read_text())
    pi["mesh_provenance"]["mesh_signature"] = "0" * 64
    (rd / "prepare_info.json").write_text(json.dumps(pi, indent=1))
    rc, out, inf, same = go("run_mesh")
    check("(c) 格子署名の記録が現行の prepare と違えば判定せずに止まる (mesh)",
          inf.get("status") == "INDETERMINATE" and (inf.get("stop") or {}).get("kind") == "mesh"
          and any("記録" in r for r in (inf.get("stop") or {}).get("reasons", [])) and not (rd / TAG / "metrics.json").exists() and same,
          str(inf.get("stop"))[:300])
    rd = clone("run_temon")
    with open(rd / "TE_MONITOR.csv", newline="") as fh:
        rows = list(csv.reader(fh))
    with open(rd / "TE_MONITOR.csv", "w", newline="") as fh:
        csv.writer(fh).writerows([r[:-1] for r in rows])          # 最後の列 (ALL_vol_lt200) を落とす
    rc, out, inf, same = go("run_temon")
    check("(c) TE_MONITOR.csv の列が現行の te_monitor.py と違えば判定せずに止まる (te_monitor)",
          inf.get("status") == "INDETERMINATE" and (inf.get("stop") or {}).get("kind") == "te_monitor"
          and any("列" in r for r in (inf.get("stop") or {}).get("reasons", [])) and not (rd / TAG / "metrics.json").exists() and same,
          str(inf.get("stop"))[:300])

    # (d) 出力先が既にある
    fp_all0 = fingerprint(C / "run_ok")
    rc, out = run_script("--case-dir", C, "--tag", TAG, "--no-summary", "--spec", spec("run_ok"))
    fp_all1 = fingerprint(C / "run_ok")
    check("(d) 出力先が既にあれば失敗 (終了コード 2)、run のファイル (既存の再判定の出力を含む) を変えない",
          rc == 2 and "既にある" in out and fp_all0 == fp_all1, out[-300:])
    # まとめ (記録のある run から)
    rc, out = run_script("--summary-only", "--case-dir", C, "--tag", TAG, "--summary-out", S / "sum_ad",
                         *sum((["--spec", spec(n)] for n in ("run_ok", "run_missing", "run_cfg")), []))
    sj = json.loads((S / "sum_ad.json").read_text()) if (S / "sum_ad.json").is_file() else {}
    st = {r["run"]: r["status"] for r in sj.get("rows", [])}
    check("(a)(b)(c) 既存の記録からまとめを書く (状態を run ごと、判定不能・不合格があれば結果 B)",
          rc == 0 and st.get("run_ok") == "JUDGED" and st.get("run_missing") == "INDETERMINATE" and st.get("run_cfg") == "INDETERMINATE"
          and sj.get("result") == ("A" if not sj.get("failures") else "B") and sj.get("result") == "B", out[-300:])
    check("(d) 複写元の run は変わらない", fingerprint(src) == src_fp0)

print(f"出力: {S}")
if not a.keep and not a.out:
    shutil.rmtree(S, ignore_errors=True)
print(f"VERDICT: {'PASS' if fails == 0 else f'FAIL ({fails})'}")
sys.exit(1 if fails else 0)
