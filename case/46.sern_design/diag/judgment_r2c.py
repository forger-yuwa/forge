#!/usr/bin/env python3
"""plan tooling-sern-te-wake-grid §5.1 R2c・§6 #11 (事前登録 2026-10-12) を実行する台本の Python 部。judgment_r2c.sh から段ごとに呼ぶ。

実 GPU で、標準経路の判定区間の段 (runner_sern.run_judgment) が「未定常だけ → 自動延長 → PASS」を通るかを 1 件確かめる。
起点 = run_1094_tewake_lsw10_m6 の本段の最終場 res_4000.h5 (run_1097 の起点と同じ。L_sw 1.0・m6_on・g3 局所化)。段階起動はやり直さない。

  prepare NEW : (1) 生産 YAML の m6_on で NEW を prepare (格子生成・変換・品質。段階起動は回さない。既存なら止める)
                (2) 格子署名・双対幾何のハッシュを起点 run の prepare_info.json の記録と照合 (違えば止める)
                (3) 本段の config (現行 runner の prepare が生成した solverConfig_main.yaml) と起点 run の solverConfig_main.yaml の差を
                    R2C_CONFIG_DIFF.txt に記録し、nStepOuter・出力系 (outStepInterval・outStepStart・output:) 以外の差があれば止める。
                    bcondConfig.yaml・species_meta.yaml・species_db_external.yaml は差があれば止める。作動点が違えば止める
                (4 の前半) 起点の場の化学種の属性から FORGE_ALLOW_UNVERIFIED_SPECIES の要否を決め (属性なし = 要る、印付き・検証済み = 要らない、
                    記録の破損 = 止める)、起点の場を NEW に複写する (sha256 を照合)。記録は R2C_PREPARE.json
  judge NEW   : (4) runner_sern.run_judgment を te_monitor 付き・強制の分岐なしで呼ぶ (N・W・延長 1 回は問題 YAML から)。
                FORGE_CUDA_BLOCKSIZE=128 が渡されていること、FORGE_ALLOW_UNVERIFIED_SPECIES が prepare の判断と一致することを確かめる。
                記録は R2C_JUDGE.json (環境の実際の値を含む)
  post NEW    : (5) check_convergence (判定区間の段ごとの session、--from-floor 起点 run も)・check_quasisteady (力 4 量・最低温度 4 領域)・
                check_floor_events (窓 W) を回して run に残し、R2C_SUMMARY.txt に段の経過と §6 #11 の分類 (機械的な表示) を書く
  rundir NEW  : NEW の run ディレクトリのパスを出す (シェル用)

終了コード: 0 = 次へ進んでよい、2 = 既存ディレクトリ、3 = 格子署名・双対幾何の不一致、4 = config・作動点の差、5 = 化学種の状態で止める、
6 = 入力・環境の欠け (記録が無い・空き容量・環境変数の不一致)、1 = その他の失敗。止めた理由は NEW/R2C_STOP.txt (NEW がある場合)。
"""
import argparse
import difflib
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "design"
TOOLS = ROOT / "solver_density_cuda" / "tools"
for p_ in (str(DESIGN), str(TOOLS)):
    if p_ not in sys.path:
        sys.path.insert(0, p_)

C46 = ROOT / "case" / "46.sern_design"
DEF_PROBLEM = "problem_3d_prod_3op_wallres_lswx10.yaml"
DEF_SRC_RUN = "run_1094_tewake_lsw10_m6"
DEF_SRC_RES = "res_4000.h5"
BLOCKSIZE = "128"
ALLOW_ENV = "FORGE_ALLOW_UNVERIFIED_SPECIES"
# 本段の config で差を許す行 (nStepOuter と出力系)。それ以外の差は止める
ALLOWED_DIFF = re.compile(r"^\s*(last:\s*\{\s*nStepOuter:\s*\d+\s*\}|nStepOuter:\s*\d+|outStepInterval:\s*\d+|outStepStart:\s*\d+|output:\s*\{.*\})\s*$")
SAME_FILES = ("bcondConfig.yaml", "species_meta.yaml", "species_db_external.yaml")
PREPARE_JSON, JUDGE_JSON, STOP_TXT, SUMMARY_TXT = "R2C_PREPARE.json", "R2C_JUDGE.json", "R2C_STOP.txt", "R2C_SUMMARY.txt"


class Stop(Exception):
    def __init__(self, code: int, msg: str):
        super().__init__(msg)
        self.code = code


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def git_head() -> dict:
    def _g(*a):
        r = subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    st = _g("status", "--porcelain")
    return {"head": _g("rev-parse", "HEAD"), "dirty_entries": None if st is None else len([x for x in st.splitlines() if x.strip()])}


def paths(a):
    case_dir = Path(a.case_dir).resolve() if a.case_dir else C46
    rd = case_dir / a.new
    src = Path(a.src_run) if Path(a.src_run).is_absolute() else case_dir / a.src_run
    prob = Path(a.problem) if Path(a.problem).is_absolute() else case_dir / a.problem
    return case_dir, rd, src, prob


def stop(rd: Path, code: int, msg: str):
    if rd.exists():
        (rd / STOP_TXT).write_text(f"R2C STOP (code {code}) {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n{msg}\n")
    raise Stop(code, msg)


# ----------------------------------------------------------------------------------------------- prepare
def config_diff(new_text: str, src_text: str) -> tuple:
    """本段の config の差 (unified diff の文字列, 許す差の行, 許さない差の行)。"""
    ud = list(difflib.unified_diff(src_text.splitlines(), new_text.splitlines(), "src/solverConfig_main.yaml",
                                   "new/solverConfig_main.yaml", lineterm=""))
    changed = [ln[1:] for ln in ud if ln[:1] in "+-" and not ln.startswith(("+++", "---"))]
    ok = [ln for ln in changed if ALLOWED_DIFF.match(ln)]
    bad = [ln for ln in changed if not ALLOWED_DIFF.match(ln)]
    return "\n".join(ud), ok, bad


def cmd_prepare(a) -> int:
    case_dir, rd, src, prob = paths(a)
    if rd.exists():
        raise Stop(2, f"{rd} は既にある (新しい run 名を渡す。既存の run は上書きしない)")
    for need in (src / "prepare_info.json", src / "solverConfig_main.yaml", src / a.src_res, prob):
        if not need.exists():
            raise Stop(6, f"入力が無い: {need}")
    free = shutil.disk_usage(case_dir).free / 2 ** 30
    if free < a.min_free_gb:
        raise Stop(6, f"{case_dir} の空き {free:.1f} GiB < {a.min_free_gb} GiB (格子・起点の場・判定区間の出力に要る)")
    from forge_design.evaluate import runner_sern as R2
    from forge_design.evaluate import runner_sern3d as R3
    import forge_species as FS
    t0 = time.time()
    info = R3.prepare(prob, rd, op=a.op)                    # (1) 格子生成・変換・品質。段階起動は回さない
    (rd / "sern.msh").unlink(missing_ok=True)                # 中間の gmsh 格子 (再生成できる。容量のため消す)
    rec = {"version": "r2c-v1", "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "prepare_s": round(time.time() - t0, 1),
           "git": git_head(), "problem": str(prob), "op": a.op, "new_run": str(rd), "src_run": str(src), "src_res": a.src_res,
           "forge_bin": None, "ok": False}
    fb = Path(os.environ.get("FORGE_BIN") or R2.FORGE_BUILD / "forge")
    rec["forge_bin"] = {"path": str(fb), "sha256": sha256(fb) if fb.exists() else None}
    rec["converter"] = R2.converter_identity()

    def _write():
        (rd / PREPARE_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False))

    # (2) 格子署名・双対幾何のハッシュ (prepare が最終の格子から再計算した値) を起点 run の記録と照合
    si = json.loads((src / "prepare_info.json").read_text())
    nm, sm = info.get("mesh_provenance") or {}, si.get("mesh_provenance") or {}
    rec["mesh_identity"] = {k: {"new": nm.get(k), "src": sm.get(k)} for k in ("mesh_signature", "dual_hash", "discretization")}
    rec["eval_method_id"] = {"new": (info.get("eval_method") or {}).get(R2.EVAL_METHOD_ID),
                             "src": (si.get("eval_method") or {}).get(R2.EVAL_METHOD_ID)}     # 記録だけ (照合は格子署名・双対で)
    _write()
    if not sm.get("mesh_signature") or not sm.get("dual_hash"):
        stop(rd, 6, f"起点 run の prepare_info.json に格子署名・双対幾何のハッシュの記録が無い ({src})")
    bad = [k for k in ("mesh_signature", "dual_hash") if nm.get(k) != sm.get(k)]
    if bad:
        stop(rd, 3, "格子が起点 run と違う: " + "; ".join(f"{k} new {nm.get(k)} / src {sm.get(k)}" for k in bad))
    so = (si.get("operating_point") or {}).get("name")
    if so != a.op:
        stop(rd, 4, f"起点 run の作動点 {so!r} が {a.op!r} と違う")
    # (3) 本段の config の差 (現行 runner の生成 vs 起点 run)、BC・化学種の設定ファイルの差
    ud, ok_lines, bad_lines = config_diff((rd / "solverConfig_main.yaml").read_text(), (src / "solverConfig_main.yaml").read_text())
    other = {}
    for name in SAME_FILES:
        a_, b_ = rd / name, src / name
        if a_.exists() or b_.exists():
            ta, tb = (a_.read_text() if a_.exists() else None), (b_.read_text() if b_.exists() else None)
            other[name] = "same" if ta == tb else ("\n".join(difflib.unified_diff((tb or "").splitlines(), (ta or "").splitlines(),
                                                                                   f"src/{name}", f"new/{name}", lineterm=""))
                                                   or "片方にだけある")
    with open(rd / "R2C_CONFIG_DIFF.txt", "w") as fh:
        fh.write(f"# 本段の config: 現行 runner の prepare ({rd.name}/solverConfig_main.yaml) と起点 run ({src.name}/solverConfig_main.yaml)\n"
                 f"# 許す差 (nStepOuter・出力系) {len(ok_lines)} 行、許さない差 {len(bad_lines)} 行\n")
        fh.write((ud or "(差なし)") + "\n")
        for name, d in other.items():
            fh.write(f"\n# {name}: {'差なし' if d == 'same' else '差あり'}\n" + ("" if d == "same" else d + "\n"))
    rec["config_diff"] = {"allowed": ok_lines, "disallowed": bad_lines, "files": {k: (v == "same") for k, v in other.items()}}
    _write()
    if bad_lines or any(v != "same" for v in other.values()):
        stop(rd, 4, "config の差が nStepOuter・出力系の外にある (R2C_CONFIG_DIFF.txt): "
             + "; ".join(bad_lines[:6] + [f"{k} が違う" for k, v in other.items() if v != "same"]))
    # (4 の前半) 起点の場の化学種の属性 → FORGE_ALLOW_UNVERIFIED_SPECIES の要否 (属性なしだけ要る。その実行だけ付ける)
    src_res = src / a.src_res
    st = FS.source_species_state(str(src_res), dirs=[str(src)])
    rec["species"] = {"state": st["state"], "attrs": st["attrs"], "why": st["why"]}
    if st["state"] == "broken":
        _write()
        stop(rd, 5, f"起点の場の化学種の記録を検証できない (照合不能): {st['why']}")
    rec["allow_unverified_species"] = st["state"] == "none"
    rec["allow_unverified_species_reason"] = (
        "起点の場に species_hash が無い (属性なし) — restart_by_index とソルバがその実行だけの許可を要する" if st["state"] == "none" else
        f"起点の場は {st['state']} (species_hash あり) — 宛先の解決結果とハッシュが一致すれば許可なしで通る (一致しなければ止まる)")
    # 起点の場を新しい run に置く (判定区間の段の継続元。run_judgment が sern.h5 へ index コピーして消す)
    import h5py
    with h5py.File(src_res, "r") as f, h5py.File(rd / R2.MESH, "r") as g:
        n_src, n_new = int(f["VALUE/ro"].shape[0]), int(g["VALUE/ro"].shape[0])
    if n_src != n_new:
        stop(rd, 3, f"起点の場の節点数 {n_src} が新しい格子の {n_new} と違う")
    if any(rd.glob("res_*.h5")):
        stop(rd, 6, f"{rd} に既に res_*.h5 がある (prepare の直後は無いはず)")
    dst = rd / a.src_res
    shutil.copyfile(src_res, dst)
    h_src, h_dst = sha256(src_res), sha256(dst)
    rec["src_res_sha256"], rec["copied_res_sha256"] = h_src, h_dst
    if h_src != h_dst:
        stop(rd, 6, f"起点の場の複写が一致しない ({h_src} / {h_dst})")
    rec["copied_records"] = []
    for f in sorted(src.glob("resolved_species_*.yaml")):     # 化学種の解決済み記録 (場の隣に置くと継承の照合に使われる)
        if not (rd / f.name).exists():
            shutil.copyfile(f, rd / f.name)
            rec["copied_records"].append(f.name)
    rec["ok"] = True
    _write()
    print(f"R2C PREPARE OK: {rd} (格子署名・双対一致、config の差は {ok_lines or 'なし'}、化学種 {st['state']}、"
          f"{ALLOW_ENV} {'要る' if rec['allow_unverified_species'] else '要らない'})")
    return 0


# ----------------------------------------------------------------------------------------------- judge
def cmd_judge(a) -> int:
    case_dir, rd, src, prob = paths(a)
    pj = rd / PREPARE_JSON
    if not pj.exists() or not json.loads(pj.read_text()).get("ok"):
        raise Stop(6, f"{pj} が無い・ok でない (prepare が通っていない)")
    prep = json.loads(pj.read_text())
    if os.environ.get("FORGE_CUDA_BLOCKSIZE") != BLOCKSIZE:
        raise Stop(6, f"FORGE_CUDA_BLOCKSIZE={os.environ.get('FORGE_CUDA_BLOCKSIZE')!r} (台本が {BLOCKSIZE} を明示して渡す)")
    allow_now = os.environ.get(ALLOW_ENV) == "1"
    if allow_now != bool(prep["allow_unverified_species"]) or (not allow_now and ALLOW_ENV in os.environ):
        raise Stop(6, f"{ALLOW_ENV}={os.environ.get(ALLOW_ENV)!r} が prepare の判断 (要る = {prep['allow_unverified_species']}) と合わない")
    from forge_design.evaluate import runner_sern as R2        # 環境を確かめてから import する (_ENV は import 時の os.environ)
    from forge_design.evaluate import runner_sern3d as R3
    from forge_design.metrics import sern_gates as G
    from forge_design.probdef import load_problem
    snaps = [f.name for _, f in R2._volume_snapshots(rd)]
    if snaps != [prep["src_res"]]:
        raise Stop(6, f"全場のスナップショットが起点の場 1 枚でない: {snaps}")
    o = load_problem(prob).raw.get("opt") or {}
    rec = {"start": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "env": R2.launch_env(), "plan": R2.judgment_plan(load_problem(prob).evaluate),
           "rc_startup": 0, "rc_startup_source": f"{src.name} の段階起動 (本段の最終場 {prep['src_res']} を起点にする。段階起動はやり直さない)",
           "forced_branch": False, "extendable_fn": f"{G.judgment_extendable.__module__}.{G.judgment_extendable.__name__}"}
    rec["plan"]["planned_steps"] = [rec["plan"]["planned_steps"][0], rec["plan"]["planned_steps"][-1], len(rec["plan"]["planned_steps"])]
    (rd / JUDGE_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    if R2.judgment_extendable is not G.judgment_extendable:   # 延長の判断を差し替えない (強制の分岐なし)
        raise Stop(1, "runner_sern.judgment_extendable が sern_gates の関数と違う (差し替えられている)")
    t0 = time.time()
    rc, out = R2.run_judgment(prob, rd, R3.collect, 0, monitor=R3.te_monitor_snapshot, precheck=R3.te_monitor_unavailable,
                              require_residual_pass=bool(o.get("require_residual_pass", False)), poll_s=a.poll_s)
    j = json.loads((rd / R2.JUDGMENT_FILE).read_text())
    rec.update({"end": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "elapsed_s": round(time.time() - t0, 1), "final_rc": rc,
                "judgment": j["verdict"], "reason": j["reason"], "n_attempts": len([x for x in j["attempts"] if "gates" in x])})
    (rd / JUDGE_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print(f"R2C JUDGMENT: {j['verdict']} ({j['reason']}) rc={rc}")
    return 0


# ----------------------------------------------------------------------------------------------- post
def _run_tool(args, out_path: Path, header: str = "") -> int:
    r = subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True)
    out_path.write_text((header + "\n" if header else "") + f"$ {' '.join(map(str, args))}\n" + r.stdout + r.stderr
                        + f"EXIT: {r.returncode}\n")
    return r.returncode


def classify(j: dict) -> str:
    """§6 #11 の分類 (機械的な表示。解釈は呼び出し側)。"""
    att = [x for x in j.get("attempts", []) if "gates" in x]
    if len(att) == 2 and att[0].get("fail_class") in ("JUDGMENT_UNSETTLED", "UNSTEADY") and att[1].get("rc") == 0 \
            and att[1].get("gates") == "PASS" and j.get("verdict") == "PASS":
        return "合格の形: 1 回目が未定常だけで不合格 → 自動延長 → 2 回目 PASS (§6 #11 の合格条件の形。解釈は上位)"
    if len(att) == 1 and j.get("verdict") == "PASS":
        return "(i) 1 回目が PASS (延長経路を通らない) → 事前登録どおり起点を替えて 1 回だけ再試行"
    if len(att) == 2 and (att[1].get("rc") != 0 or att[1].get("fail_class") == "DIVERGED"):
        return "(ii) 2 回目が発散・rc ≠ 0 → 経路の欠陥を疑い、上位に諮るまで先へ進まない"
    if len(att) == 1:
        return f"(ii) 1 回目が未定常以外で不合格 ({att[0].get('fail_class')}) → 上位に諮るまで先へ進まない"
    if len(att) == 2 and j.get("verdict") == "INDETERMINATE" and att[1].get("fail_class") in ("JUDGMENT_UNSETTLED", "UNSTEADY"):
        return "(iii) 2 回目も未定常 (INDETERMINATE) → run_1099 と力の平均・a・d を並べて上位に諮る"
    return f"分類できない (段 {len(att)} 回、JUDGMENT {j.get('verdict')}) → 上位に諮る"


def cmd_post(a) -> int:
    case_dir, rd, src, prob = paths(a)
    from forge_design.evaluate import runner_sern as R2
    from forge_design.probdef import load_problem
    jf = rd / R2.JUDGMENT_FILE
    if not jf.exists():
        raise Stop(6, f"{jf} が無い (judge が終わっていない)")
    j = json.loads(jf.read_text())
    W = int(j["options"]["window_steps"])
    att = [x for x in j["attempts"] if "gates" in x]
    lines = [f"R2C SUMMARY {time.strftime('%Y-%m-%dT%H:%M:%S%z')}  run {rd}", f"JUDGMENT: {j['verdict']} ({j['reason']})",
             f"分類 (§6 #11、機械的な表示): {classify(j)}", ""]
    # 残差 (判定区間の段ごとの session = 同一設定の区間)。--from-floor は起点 run (同じ格子・同じ設定の本段) の末尾床と比べる
    for k in range(1, len(att) + 1):
        csvp = rd / f"residual_judgment{k}.csv"
        if not csvp.exists():
            lines.append(f"残差 {k} 回目: {csvp.name} が無い")
            continue
        with tempfile.TemporaryDirectory() as td:
            shutil.copyfile(csvp, Path(td) / "residual_history.csv")
            hdr = f"# 判定区間 = 判定区間の段 {k} 回目の session (同一設定の {j['options']['steps']} step。{csvp.name})"
            rc1 = _run_tool([TOOLS / "check_convergence.py", td], rd / f"R2C_CONVERGENCE_judgment{k}.txt", hdr)
            rc2 = None
            if (src / "residual_history.csv").exists():
                rc2 = _run_tool([TOOLS / "check_convergence.py", td, "--from-floor", src], rd / f"R2C_CONVERGENCE_judgment{k}_from_floor.txt",
                                hdr + f"\n# --from-floor {src} (起点 run の本段の残差の末尾床)")
        lines.append(f"残差 {k} 回目: check_convergence rc {rc1} / --from-floor rc {rc2} (R2C_CONVERGENCE_judgment{k}*.txt)")
    # 準定常 (力 4 量・最低温度 4 領域、--tail 0.5)。延長したら 1 回目の退避分も
    for tag in ([""] + [f"_judgment{k}" for k in range(1, len(att))]):
        fh, mon = rd / f"force_history{tag}.csv", rd / f"TE_MONITOR{tag}.csv"
        if fh.exists():
            _run_tool([TOOLS / "check_quasisteady.py", "--series-csv", fh, "--series-cols", "C_T,C_T_with_shear,C_L,C_M", "--tail", "0.5"],
                      rd / f"R2C_QUASISTEADY_forces{tag}.txt")
        if mon.exists():
            _run_tool([TOOLS / "check_quasisteady.py", "--series-csv", mon, "--series-cols", "R_TE_Tmin,R_SE_Tmin,RET_Tmin,ALL_Tmin", "--tail", "0.5"],
                      rd / f"R2C_QUASISTEADY_Tmin{tag}.txt")
        for name in (f"R2C_QUASISTEADY_forces{tag}.txt", f"R2C_QUASISTEADY_Tmin{tag}.txt"):
            p_ = rd / name
            if p_.exists():
                ov = [ln.strip() for ln in p_.read_text().splitlines() if ln.strip().startswith("OVERALL")]
                lines.append(f"準定常 {name}: {ov[-1] if ov else '(OVERALL 行なし)'}")
    # 床事象 (最後の session の末尾 W step)
    rcf = _run_tool([TOOLS / "check_floor_events.py", rd, "--window-steps", W], rd / "R2C_FLOOR_EVENTS.txt")
    fv = [ln for ln in (rd / "R2C_FLOOR_EVENTS.txt").read_text().splitlines() if "VERDICT" in ln]
    lines.append(f"床事象 check_floor_events --window-steps {W}: rc {rcf} {fv[-1] if fv else ''}")
    lines.append("")
    for x in att:
        fw, tp, fe = x.get("force_window") or {}, x.get("temperature") or {}, x.get("floor_events") or {}
        lines.append(f"段 {x['attempt']}: 継続元 {x.get('restart_from')} rc {x.get('rc')} gates {x.get('gates')}/{x.get('fail_class')} "
                     f"judgment {x.get('judgment')} 床 {fe.get('verdict')} {fe.get('window')} 更新 {fe.get('n_updates_checked')}")
        for q, v in (fw.get("quantities") or {}).items():
            lines.append(f"   {q:15s} 平均 {v.get('mean'):.7g}  a {v.get('a'):.2e}  d {v.get('d'):.2e}  (許容 {v.get('tol'):.0e})  {v.get('steady')}  "
                         f"{'OK' if v.get('ok') else 'NG'}")
        for R, v in (tp.get("regions") or {}).items():
            lines.append(f"   {R:5s} 最低 {v.get('Tmin_window'):.1f} K (step {v.get('Tmin_step')})  T<150 延べ {v.get('n_lt150_window_sum')}  {v.get('steady')}")
        if x.get("archived"):
            lines.append(f"   退避: {', '.join(x['archived'])}  延長の理由: {x.get('extend_reason')}")
    m = json.loads((rd / "metrics.json").read_text()) if (rd / "metrics.json").exists() else {}
    g = m.get("gates") or {}
    same = bool(att) and g.get("verdict") == att[-1].get("gates") and g.get("fail_class") == att[-1].get("fail_class")
    lines += ["", f"metrics.json (正本 = 最後の判定区間): gates {g.get('verdict')}/{g.get('fail_class')} judgment "
                  f"{(g.get('judgment') or {}).get('verdict')} — JUDGMENT.json の最後の段と{'一致' if same else '不一致'}",
              f"環境 (JUDGMENT.json の env.named): {j.get('env', {}).get('named')}"]
    (rd / SUMMARY_TXT).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


def cmd_rundir(a) -> int:
    print(paths(a)[1])
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="plan tooling-sern-te-wake-grid §5.1 R2c・§6 #11 の台本 (judgment_r2c.sh から呼ぶ)")
    ap.add_argument("stage", choices=("prepare", "judge", "post", "rundir"))
    ap.add_argument("new", help="新しい run の名前 (case/46.sern_design/<new>。既存なら止める)")
    ap.add_argument("--problem", default=DEF_PROBLEM, help=f"問題 YAML (既定 {DEF_PROBLEM}、相対は --case-dir から)")
    ap.add_argument("--op", default="m6_on")
    ap.add_argument("--src-run", default=DEF_SRC_RUN, help=f"起点 run (既定 {DEF_SRC_RUN}、相対は --case-dir から)")
    ap.add_argument("--src-res", default=DEF_SRC_RES, help=f"起点の場 (既定 {DEF_SRC_RES})")
    ap.add_argument("--case-dir", default=None, help="run を置くディレクトリ (既定 case/46.sern_design。台本の試験用)")
    ap.add_argument("--poll-s", type=float, default=30.0, help="判定区間の段でスナップショットを監視する間隔 [s]")
    ap.add_argument("--min-free-gb", type=float, default=6.0, help="prepare の前に要る空き容量 [GiB]")
    a = ap.parse_args(argv)
    try:
        return {"prepare": cmd_prepare, "judge": cmd_judge, "post": cmd_post, "rundir": cmd_rundir}[a.stage](a)
    except Stop as e:
        print(f"R2C STOP ({a.stage}, code {e.code}): {e}", file=sys.stderr)
        return e.code


if __name__ == "__main__":
    raise SystemExit(main())
