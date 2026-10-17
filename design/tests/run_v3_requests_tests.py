#!/usr/bin/env python3
"""§6 V3 の 13a (依頼の課題集の機械の側、GPU なし) の試験 (plan tooling-design-orchestrator §6「V3 の 13a の事前登録」2026-10-18;
codex 諮問 v3)。依頼の原文・期待・固定する条件は `design/tests/fixtures/v3_requests/requests.yaml`、安価な負例と正例の対応表は
`cheap_negatives.json`。**13a の合格は V3 の完了ではない** (13b のモデルを回す側の後)。

  (K) 安価な負例と正例: 対応表の各検査を、その試験を回した出力の行で確かめる (未実行・欠け・skip・FAIL は合格にしない)
  N1・N2: 版 2 のカード → wt_card (承認の前の問題 YAML) → 承認 → 生成し直して approval.check → hand_mesh.py build --approval → check (品質 PASS)
  N3: 版 2 のカード → wt_card → 承認 → 生成し直して approval.check → mesh_stage.check_template (品質 PASS; 格子の部品の成功)
  S1: 旧データの非承認のモードの公差の判定 (試走 A run_0003_ns_p0 の写し): iteration 1 → back_to_correction (V5b・V5c FAIL)、iteration 3 → rejected
  S2: N1 のカードから sources.wall を除く → preview に wall の出所の欠け、write しない (止まった理由が問題 YAML の欠けでない)
  S3: 予算 GPU 0.01 時間・Euler の見積もり 4/60 の campaign、格子の段が済んだ台帳 → driver が Euler の前に予算で止まる (試行を足さない、起動 0 回)

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_v3_requests_tests.py [--only K|N|S]
所要: 40 分程度 (N1〜N3 の生成と格子、(K) の試験の中の run_mesh_limits_tests)。追跡ファイルに変更があると (K) の run_gate_env_tests が止まる (commit してから)。
"""
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DESIGN = REPO / "design"
_FB = REPO / "solver_density_cuda" / ".build-native" / "card4a" / "forge"
if not os.environ.get("FORGE_BIN") and _FB.is_file():
    os.environ["FORGE_BIN"] = str(_FB)          # 準備の経路 (N3) の初期場の種の記録に --resolve-species を持つ forge が要る (起動はしない)
os.environ.setdefault("FORGE_CONVERTER", str(REPO / "solver_density_cuda" / ".build-native" / "card4a" / "convertGmshToForge"))
sys.path.insert(0, str(DESIGN))
from forge_design.campaign import campaign as CP  # noqa: E402
from forge_design.campaign import driver as DR  # noqa: E402
from forge_design.campaign import mesh_stage as MS  # noqa: E402
from forge_design.campaign.ledger import Ledger  # noqa: E402
from forge_design.feedback import wt_qualify as WQ  # noqa: E402
from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

FIX = HERE / "fixtures" / "v3_requests"
REQ = yaml.safe_load((FIX / "requests.yaml").read_text())
CN = json.loads((FIX / "cheap_negatives.json").read_text())
HM = REPO / "case" / "71.m55_combustion_wt" / "_mesh" / "hand_mesh.py"
RATIO = REPO / "case" / "71.m55_combustion_wt" / "_mesh" / "c45" / "h_over_est_ratio_summary.json"
PY = sys.executable
ONLY = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else "KNS"
USER_SECTIONS = ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall", "tolerances", "p_ambient", "budget")
FAIL = 0
RESULTS = {}


def check(name, ok, detail=""):
    global FAIL
    print(("ok  " if ok else "FAIL") + " " + name + (f"  [{str(detail)[:400]}]" if detail else ""), flush=True)
    RESULTS[name] = bool(ok)
    if not ok:
        FAIL += 1


def req(rid):
    return next(r for r in REQ["requests"] if r["id"] == rid)


td = Path(tempfile.mkdtemp(prefix="v3_requests_"))
_clean = {}


def clean_repo() -> Path:
    if "r" not in _clean:
        r = td / "clean_repo"
        r.mkdir()
        for c in (["init", "-q"], ["-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q", "--allow-empty", "-m", "t"]):
            subprocess.run(["git", "-C", str(r), *c], check=True, capture_output=True)
        _clean["r"] = r
    return _clean["r"]


def v2_card(r: dict, dst: Path, drop_source=None, budget=None) -> Path:
    d = yaml.safe_load((REPO / r["card"]["source"]).read_text())
    d["format"] = W.CARD_FORMAT_V2
    d["budget"] = budget or {"gpu_hours": 8.0, "disk_gb": 50.0}
    src = {k: {"kind": "user", "ref": f"依頼 {r['id']} の原文"} for k in USER_SECTIONS}
    if r["card"].get("mesh_limits"):
        d["mesh_limits"] = dict(r["card"]["mesh_limits"])
        src["mesh_limits"] = {"kind": "user", "ref": r["card"]["mesh_limits_ref"]}
    if drop_source:
        src.pop(drop_source)
    d["sources"] = src
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    return dst


def generate(card: Path, out: Path):
    r = subprocess.run([PY, "-m", "forge_design.intake.wt_card", str(card), "--out", str(out)], cwd=DESIGN, capture_output=True, text=True)
    return r.returncode == 0 and out.is_file(), r.stderr[-400:]


def approve(card: Path, prob: Path, out: Path):
    pv = AP.preview(card, prob, repo_dir=clean_repo(), judge_env=False)
    pp = out.with_suffix(".preview.json")
    pp.write_text(json.dumps(pv, ensure_ascii=False))
    try:
        AP.write(pp, "試験の仮の承認 (run_v3_requests_tests)", out)
        return pv, None
    except (AP.ApprovalError, ValueError) as e:
        return pv, e


def quality_lines(rec_path: Path) -> dict:
    rec = json.loads(rec_path.read_text())
    return {t: (rec.get(t) or {}).get("quality_verdict") for t in ("ns", "euler")}


try:
    # ------------------------------------------------------------------------------------------- (K) 安価な負例と正例
    if "K" in ONLY:
        print("== (K) 安価な負例と正例 (対応表)", flush=True)
        outs = {}
        for t in sorted({e["test"] for e in CN["negatives"] + CN["positives"]}):
            r = subprocess.run([PY, str(HERE / t)], capture_output=True, text=True, timeout=3600)
            outs[t] = (r.returncode, r.stdout + r.stderr)
        for grp, items in (("負例", CN["negatives"]), ("正例", CN["positives"])):
            for e in items:
                rc, out = outs[e["test"]]
                lines = [ln for ln in out.splitlines() if (ln.startswith("ok   ") or ln.startswith("FAIL ")) and ln[5:].startswith(e["check"])]
                ok = bool(lines) and all(ln.startswith("ok   ") for ln in lines)
                check(f"(K) {grp} {e['kind']}: {e['test']} の「{e['check']}」が実行されて ok", ok,
                      f"該当 {len(lines)} 行、試験の終了コード {rc}" + ("" if lines else "、未実行 (skip・欠け)"))
        check("(K) 対応表の試験がどれも終了コード 0 (飛ばした・FAIL が無い)", all(rc == 0 for rc, _ in outs.values()),
              {t: rc for t, (rc, _) in outs.items()})

    # ------------------------------------------------------------------------------------------- N1・N2・N3
    gen_n1 = None
    if "N" in ONLY or "S" in ONLY:
        for rid in ("N1", "N2", "N3") if "N" in ONLY else ("N1",):
            r = req(rid)
            print(f"== {rid}: {r['text']}", flush=True)
            w = td / rid
            card = v2_card(r, w / "card_p0.yaml")
            pc = w / "problem_p0_card.yaml"
            ok, err = generate(card, pc)
            check(f"{rid} 承認の前の問題 YAML を生成 (wt_card)", ok, err)
            if not ok:
                continue
            if rid == "N1":
                gen_n1 = (card, pc)
            if "N" not in ONLY:
                continue
            ap = w / "card_approval.json"
            pv, err = approve(card, pc, ap)
            check(f"{rid} 承認の preview に問題が無く承認の記録を書ける", err is None and not pv["problems"], pv["problems"] or err)
            if err is not None:
                continue
            pc2 = w / "regen" / "problem_p0_card.yaml"
            pc2.parent.mkdir()
            ok2, err2 = generate(card, pc2)
            ck = AP.check(ap, card, pc2) if ok2 else {"status": "生成できない", "problems": [err2]}
            check(f"{rid} 承認の後に生成し直した問題 YAML が approval.check OK", ck["status"] == "OK", ck.get("problems", []) + ck.get("diffs", []))
            if rid in ("N1", "N2"):
                prob = w / "problem_p0.yaml"
                rb = subprocess.run([PY, str(HM), "build", "--design", str(DESIGN), "--card", str(card), "--problem", str(pc), "--out-yaml", str(prob),
                                     "--record", str(w / "build_record.json"), "--axis-cap-frac", "0.03", "--ratio-summary", str(RATIO),
                                     "--plan-item", "V3 13a", "--approval", str(ap)], capture_output=True, text=True)
                rc_ = subprocess.run([PY, str(HM), "check", "--design", str(DESIGN), "--card", str(card), "--problem", str(prob),
                                      "--record", str(w / "check_record.json"), "--run-dir", str(w / "chk")], capture_output=True, text=True) \
                    if rb.returncode == 0 else None
                q = quality_lines(w / "check_record.json") if rc_ is not None and rc_.returncode == 0 else {}
                check(f"{rid} hand_mesh.py build --approval → check で NS・Euler の品質 PASS",
                      rb.returncode == 0 and rc_ is not None and rc_.returncode == 0 and all((q.get(t) or "").startswith("VERDICT: PASS") for t in ("ns", "euler")),
                      f"build rc {rb.returncode} {rb.stderr[-300:]} / check {None if rc_ is None else rc_.returncode} {q}")
            else:
                try:
                    rec = MS.check_template(pc, w / "ms", approval=ap, card=card)
                    q = {t: rec[t]["quality"] for t in ("ns", "euler")}
                    err3 = None
                except Exception as e:                                        # noqa: BLE001
                    q, err3 = {}, f"{type(e).__name__}: {e}"
                check(f"{rid} mesh_stage.check_template で NS・Euler の品質 PASS (格子の部品の成功)",
                      err3 is None and all((q.get(t) or "").startswith("VERDICT: PASS") for t in ("ns", "euler")), err3 or q)

    # ------------------------------------------------------------------------------------------- S1・S2・S3
    if "S" in ONLY:
        r = req("S1")
        print(f"== S1: {r['text']}", flush=True)
        ra = REPO / "case" / "90.qs_trial_a" / "run_0003_ns_p0"
        run = td / "S1" / ra.name
        run.mkdir(parents=True)
        for f in ra.iterdir():
            if f.is_dir():
                if f.name in ("stop_rule_v2",):
                    shutil.copytree(f, run / f.name, symlinks=True, copy_function=lambda s, d: os.symlink(Path(s).resolve(), d)
                                    if s.endswith(".h5") else shutil.copy2(s, d))
                continue
            if f.suffix in (".h5", ".xmf"):
                (run / f.name).symlink_to(f.resolve())
            else:
                shutil.copy2(f, run / f.name)
        tol = r["tolerances"]
        for i, e in enumerate(r["expect"]):
            rec = WQ.qualify(run, REPO / "case/90.qs_trial_a/card_p0.yaml", REPO / "case/90.qs_trial_a/problem_p0.yaml",
                             iteration=e["iteration"], max_iterations=e["max_iterations"],
                             cli_tols={k: tol[k] for k in ("exit_core_M_tol_rel", "line_tol_maxabs_pct", "line_tol_halfamp_pct")},
                             wall_txt=td / "S1" / f"wr{i}.txt", series_dir=td / "S1", judge_dir=td / "S1" / f"jr{i}")
            c = rec["checks"]
            want_fail = all(c[k]["verdict"] == "FAIL" for k in e.get("fail", []))
            check(f"S1 iteration {e['iteration']}/{e['max_iterations']} → {e['result']}" + (f" (原因 {e['fail']} の FAIL)" if e.get("fail") else ""),
                  rec["result"] == e["result"] and want_fail and c["run_verdict"]["status"] == "PASS",
                  f"{rec['result']} V5b {c['V5b'].get('verdict')} M {c['V5b'].get('M')} V5c {c['V5c'].get('verdict')} "
                  f"{c['V5c'].get('maxabs_pct')} 解析の成立 {c['run_verdict']['status']} {rec['missing']}")

        r = req("S2")
        print(f"== S2: {r['text']}", flush=True)
        if gen_n1 is None:
            check("S2 前提: N1 の生成した問題 YAML", False, "N1 の生成が無い")
        else:
            card_s2 = v2_card(req("N1"), td / "S2" / "card_p0.yaml", drop_source=r["card"]["drop_source"])
            pv, err = approve(card_s2, gen_n1[1], td / "S2" / "card_approval.json")
            probs = pv["problems"]
            check("S2 sources.wall を除く → preview に wall の出所の欠け、write しない (理由は問題 YAML の欠けでない)",
                  err is not None and any(p.startswith("wall") and "出所" in p for p in probs) and not any("問題 YAML" in p and "無い" in p for p in probs)
                  and not (td / "S2" / "card_approval.json").exists(), probs[:3])

        r = req("S3")
        print(f"== S3: {r['text']}", flush=True)
        if gen_n1 is None:
            check("S3 前提: N1 の生成した問題 YAML", False, "N1 の生成が無い")
        else:
            case = td / "S3"
            card_s3 = v2_card(req("N1"), case / "card_p0.yaml", budget=r["card"]["budget"])
            pc = case / "problem_p0_card.yaml"
            shutil.copy(gen_n1[1], pc)
            pv, err = approve(card_s3, pc, case / "card_approval.json")
            check("S3 前提: 予算 0.01 GPU 時間の承認", err is None, pv["problems"] or err)
            est = {k: {**v, "basis": "V3 S3 の試験の入力 (生産の見積もりではない)"} for k, v in r["campaign_estimates"].items()}
            camp = {"format": "forge.campaign/v1", "name": "v3_s3", "role": "fixer",
                    "approval": {"file": "card_approval.json", "sha256": AP.sha256_file(case / "card_approval.json")},
                    "recipe": {"name": "axismach.wt_card_chain", "version": "v1"}, "case_dir": ".", "ledger": "_ledger", "estimates": est}
            (case / "campaign.yaml").write_text(yaml.safe_dump(camp, sort_keys=False, allow_unicode=True))
            c = CP.load(case / "campaign.yaml")
            with Ledger(c["ledger"]) as L:
                L.init(c, {"gpu_hours": 0.01, "disk_gb": 50.0})
                L.doc["attempts"].append({"id": "0001-mesh_p0", "step": "mesh_p0", "state": "ingested", "run_dir": None,
                                          "interpretation": {"ok": True, "branch": None, "detail": "S3 の前提 (格子の段は済んだ)"}})
                L.doc["ingested"].append("0001-mesh_p0")
                L.save()
            ad = Path(c["ledger"]) / "attempts" / "0001-mesh_p0"
            ad.mkdir(parents=True)
            (ad / "result.json").write_text(json.dumps({"step": "mesh_p0", "ok": True, "branch": None, "detail": "S3 の前提", "outputs": {}}))
            reached = []
            orig = DR.Driver._start_unit
            DR.Driver._start_unit = lambda self, a: reached.append(a["id"]) or (_ for _ in ()).throw(RuntimeError("到達 (起動の代わり)"))
            try:
                try:
                    out = DR.Driver(case / "campaign.yaml", DR.recipe_for(c, REPO), systemd=False, design_dir=DESIGN).run()
                except RuntimeError as e:
                    out = {"status": f"到達 ({e})"}
            finally:
                DR.Driver._start_unit = orig
            led = json.loads((Path(c["ledger"]) / "ledger.json").read_text())
            eu = [a for a in led["attempts"] if a["step"].startswith("euler")]
            check("S3 予算 0.01 GPU 時間・Euler の見積もり 4/60 → Euler の前に予算で止まる (試行を足さない・起動 0 回・格子の試行と止まった記録は残る)",
                  out.get("status") == "stopped" and "予算" in str(out.get("why")) and not eu and not reached
                  and any(a["id"] == "0001-mesh_p0" for a in led["attempts"]) and len(led.get("stops") or []) == 1,
                  f"{out} 試行 {[a['id'] for a in led['attempts']]} stops {len(led.get('stops') or [])}")
            c2 = dict(camp)
            c2.pop("estimates")
            (case / "campaign_noest.yaml").write_text(yaml.safe_dump(c2, sort_keys=False, allow_unicode=True))
            try:
                CP.load(case / "campaign_noest.yaml")
                ok, why = False, "読めた"
            except ValueError as e:
                ok, why = "GPU の段" in str(e), str(e)
            check("S3 GPU の段の見積もりが無い campaign は読み込みで止まる", ok, why[:200])
finally:
    shutil.rmtree(td, ignore_errors=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
