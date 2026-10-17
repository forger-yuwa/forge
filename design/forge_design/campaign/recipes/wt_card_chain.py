"""recipe `axismach.wt_card_chain/v1` (plan tooling-design-orchestrator §4.5.1 の段と分岐、§5.1 #6c)。今の手順書
`procedures/design-chain-quickstart.md` §3〜§9 の流れを、driver が回す段の列にしたもの。

段 (pass k = 0, 1, …):

    p0: mesh_p0 → euler_p0 → ns_p0 → calibrate_p0 (最初の較正; 較正し直しの回数に数えない)
    pk (k ≥ 1): problem_pk → mesh_pk → euler_pk → ns_pk → qualify_pk
        qualify の結果: incomplete・make_deliverables → deliverables_pk (成果物 → 最後の判定) / back_to_correction・recalibrate →
        recalibrate_pk (回数 +1、上限は承認の run_policy) → p(k+1) / accepted → 終わり / ほか → 止まる

- p0 の問題 YAML は承認の記録の問題 YAML (承認の前に `wt_card` で作ったもの)、カードは承認の原本。
- 各段の処理は `python -m forge_design.campaign.recipes.wt_card_chain exec '<段の定義の JSON>'` (試行の service の中で回る)。段の処理は
  試行のディレクトリに `result.json` ({step, ok, branch, detail, verdicts, outputs {パス: sha256}, inputs}) を書く。
- `next_step` は台帳の取り込み済みの試行の `result.json` から位置を決める (driver が落ちても結果のファイルから決まる)。前の段が作った入力は
  使う直前に sha256 を照合し、違えば止まる (変更の要求 `note` に新しい sha256 と理由があれば通す; fixer の役割で直したとき)。
- 格子の段は case/71 の台本 `hand_mesh.py` (`--axis-cap-frac 0.03`、手順書 §3)。節点数の上限は承認の記録 (`--approval`、カードの
  mesh_limits; plan §4.6・§5.1 #8)。一般化は plan §5.1 #10。
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import yaml

NAME, VERSION = "axismach.wt_card_chain", "v1"
MESH_SCRIPT_REL = "case/71.m55_combustion_wt/_mesh/hand_mesh.py"
AXIS_CAP_FRAC = 0.03
RUN_ENV_KEYS = ("FORGE_BIN", "FORGE_CONVERTER", "FORGE_CUDA_BLOCKSIZE", "FREECADCMD", "FORGE_MESH_PY", "FORGE_PPTX_PY")
RUN_ENV_FILES = ("FORGE_BIN", "FORGE_CONVERTER", "FREECADCMD", "FORGE_MESH_PY", "FORGE_PPTX_PY")
RESULT_FILE = "result.json"


def sha256(p) -> str | None:
    p = Path(p)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _approval(campaign) -> dict:
    return json.loads(Path(campaign["approval"]["file"]).read_text())


def _next_run_name(case_dir: Path, step: str, taken: set) -> Path:
    nums = [int(m.group(1)) for d in case_dir.glob("run_*") if (m := re.match(r"run_(\d{4})", d.name))]
    nums += [int(re.match(r"run_(\d{4})", Path(t).name).group(1)) for t in taken if re.match(r"run_(\d{4})", Path(t).name)]
    return case_dir / f"run_{(max(nums) + 1 if nums else 1):04d}_{step}"


class Recipe:
    name, version = NAME, VERSION

    def __init__(self, repo_dir=None, python=None):
        self.repo = Path(repo_dir) if repo_dir else Path(__file__).resolve().parents[4]
        self.design = self.repo / "design"
        self.python = python or sys.executable

    # ---- 走らせる環境 (driver が初回に台帳へ写し、再開で照合する)
    def run_env(self) -> dict:
        out = {}
        for k in RUN_ENV_KEYS:
            v = os.environ.get(k)
            out[k] = {"value": v, "sha256": sha256(v) if (v and k in RUN_ENV_FILES) else None}
        return out

    # ---- 結果のファイルから位置を決める
    def _results(self, doc, campaign) -> list:
        led = Path(campaign["ledger"])
        out = []
        for a in doc["attempts"]:
            if a["id"] not in doc["ingested"]:
                continue
            f = led / "attempts" / a["id"] / RESULT_FILE
            r = json.loads(f.read_text()) if f.is_file() else None
            out.append((a, r))
        return out

    def _accepted_change(self, doc, path: str, sha: str | None) -> bool:
        return any(rq.get("kind") == "note" and rq.get("path") == path and rq.get("input_sha256") == sha and rq.get("reason")
                   for rq in doc.get("requests") or [])

    def _verify_inputs(self, doc, done: dict, paths: list) -> str | None:
        """前の段が作った入力の sha256 を使う直前に照合する。違えば理由 (変更の要求があれば通す)。"""
        made = {}
        for r in done.values():
            made.update((r or {}).get("outputs") or {})
        for p in paths:
            p = str(p)
            if p in made:
                cur = sha256(p)
                if cur != made[p] and not self._accepted_change(doc, p, cur):
                    return f"前の段が作った入力 {p} が変わった (sha256 {str(cur)[:16]}… ≠ 記録 {str(made[p])[:16]}…; 直したなら note で理由と新しい sha256 を出す)"
        return None

    def next_step(self, doc, campaign):
        case = Path(campaign["case_dir"])
        led = Path(campaign["ledger"])
        apr = _approval(campaign)
        rp = (apr.get("protected") or {}).get("run_policy") or {}
        max_it, max_recal = int(rp.get("max_iterations", 3)), int(rp.get("max_recalibrations", 1))
        res = self._results(doc, campaign)
        if res and not (res[-1][1] or {}).get("ok"):
            a, r = res[-1]
            return {"stop": f"段 {a['step']} が合格しなかった ({(r or {}).get('detail') or '結果のファイルが無い'})"}
        done = {}
        for a, r in res:
            if r and r.get("ok"):
                done[a["step"]] = r
        taken = {a.get("run_dir") for a in doc["attempts"] if a.get("run_dir")}
        card = {0: apr["card"]["path"]}
        prob_card = {0: apr["problem"]["path"]}
        recal_done = sum(1 for s in done if s.startswith("recalibrate_p"))
        k = 0
        while True:
            pk = f"p{k}"
            card_k = card.get(k) or str(case / f"card_{pk}.yaml")
            prob_card_k = prob_card.get(k) or str(case / f"problem_{pk}_card.yaml")
            prob_k, prob_eu_k = str(case / f"problem_{pk}.yaml"), str(case / f"problem_{pk}_euler.yaml")
            seq = ([] if k == 0 else [f"problem_{pk}"]) + [f"mesh_{pk}", f"euler_{pk}", f"ns_{pk}"]
            seq += ["calibrate_p0"] if k == 0 else [f"qualify_{pk}"]
            for st in seq:
                if st in done:
                    continue
                return self._spec(st, k, doc, campaign, done, taken, card_k, prob_card_k, prob_k, prob_eu_k, recal_done, max_it)
            if k == 0:
                k = 1
                if k + 1 > max_it:
                    return {"stop": f"反復の上限 {max_it} (承認の run_policy)"}
                continue
            q = done[f"qualify_{pk}"]
            br = q.get("branch")
            if br == "accepted":
                return None
            if br == "make_deliverables":
                if f"deliverables_{pk}" in done:
                    return None if done[f"deliverables_{pk}"].get("branch") == "accepted" else {"stop": "最後の判定が accepted でない"}
                return self._spec(f"deliverables_{pk}", k, doc, campaign, done, taken, card_k, prob_card_k, prob_k, prob_eu_k, recal_done, max_it)
            if br == "recalibrate":
                if f"recalibrate_{pk}" not in done:
                    if recal_done >= max_recal:
                        return {"stop": f"較正し直しの上限 {max_recal} (承認の run_policy)"}
                    return self._spec(f"recalibrate_{pk}", k, doc, campaign, done, taken, card_k, prob_card_k, prob_k, prob_eu_k, recal_done, max_it)
                k += 1
                if k + 1 > max_it:
                    return {"stop": f"反復の上限 {max_it} (承認の run_policy)"}
                continue
            return {"stop": f"判定 {q.get('detail')} (人に返す)"}

    def _spec(self, step, k, doc, campaign, done, taken, card_k, prob_card_k, prob_k, prob_eu_k, recal_done, max_it) -> dict:
        case = Path(campaign["case_dir"])
        led = Path(campaign["ledger"])
        aid = f"{len(doc['attempts']) + 1:04d}-{step}"
        sd = {"step": step, "k": k, "case_dir": str(case), "approval": campaign["approval"]["file"], "card": card_k,
              "problem_card": prob_card_k, "problem": prob_k, "problem_euler": prob_eu_k, "repo": str(self.repo),
              "out": str(led / "attempts" / aid / RESULT_FILE), "iteration": k + 1, "recalibrations_done": recal_done}
        gpu, run_dir, owned, uses = False, None, [], []
        kind = step.split("_p")[0]
        if kind == "problem":
            uses = [card_k]
            owned = [prob_card_k, prob_card_k.replace("_card.yaml", "_card_euler.yaml")]
        elif kind == "mesh":
            uses = [card_k, prob_card_k]
            owned = [prob_k, prob_eu_k, str(case / f"_mesh_p{k}")]
        elif kind == "euler":
            uses = [prob_eu_k]
            run_dir = _next_run_name(case, step, taken)
            gpu = True
        elif kind == "ns":
            uses = [prob_k, done[f"euler_p{k}"]["euler_ref"]["res"]]
            run_dir = _next_run_name(case, step, taken)
            sd["euler_ref"] = done[f"euler_p{k}"]["euler_ref"]
            if k >= 1:
                sd["ic_file"] = done[f"ns_p{k - 1}"]["final_field"]
                uses.append(sd["ic_file"]["path"])
            gpu = True
        elif kind in ("calibrate", "recalibrate"):
            sd["ns_run"] = done[f"ns_p{k}"]["run_dir"]
            sd["new_card"] = str(case / f"card_p{k + 1}.yaml")
            uses = [card_k]
            owned = [sd["new_card"]]
        elif kind in ("qualify", "deliverables"):
            sd["ns_run"] = done[f"ns_p{k}"]["run_dir"]
            sd["prev_run"] = done[f"ns_p{k - 1}"]["run_dir"] if k >= 1 else None
            uses = [card_k, prob_k]
        if run_dir is not None:
            sd["run_dir"] = str(run_dir)
            owned = [str(run_dir)]
        bad = self._verify_inputs(doc, done, [u if isinstance(u, str) else u for u in uses])
        if bad:
            return {"stop": bad}
        env = {k_: v["value"] for k_, v in (doc.get("run_env") or {}).items() if v.get("value")}
        return {"step": step, "command": [self.python, "-m", "forge_design.campaign.recipes.wt_card_chain", "exec", json.dumps(sd, ensure_ascii=False)],
                "gpu": gpu, "run_dir": None if run_dir is None else str(run_dir), "owned_paths": owned, "estimate": {}, "cwd": str(self.design),
                "env": {**env, "PYTHONPATH": str(self.design), "PYTHONDONTWRITEBYTECODE": "1"}}

    def interpret(self, spec, attempt, status) -> dict:
        sd = json.loads(spec["command"][-1])
        f = Path(sd["out"])
        r = json.loads(f.read_text()) if f.is_file() else None
        if r is None:
            return {"ok": False, "detail": f"段 {spec['step']}: 結果のファイルが無い (試行の状態 {status.get('state')})"}
        return {"ok": bool(r.get("ok")) and status.get("state") == "finished", "branch": r.get("branch"), "detail": r.get("detail")}


# ------------------------------------------------------------------------------------------------ 段の処理 (試行の service の中)

def _run(cmd, log_name, cwd, sd) -> subprocess.CompletedProcess:
    r = subprocess.run([str(c) for c in cmd], cwd=str(cwd), capture_output=True, text=True)
    Path(sd["out"]).parent.mkdir(parents=True, exist_ok=True)
    with open(Path(sd["out"]).parent / f"{log_name}.log", "a") as f:
        f.write(f"$ {' '.join(map(str, cmd))}\n{r.stdout}\n{r.stderr}\nrc={r.returncode}\n")
    return r


def _nonfinite(res_h5) -> dict:
    import h5py
    import numpy as np
    with h5py.File(res_h5, "r") as h:
        return {k: int((~np.isfinite(np.asarray(h["VALUE/" + k][:], float))).sum()) for k in h["VALUE"]}


def _last_res(run_dir) -> Path | None:
    rs = [p for p in Path(run_dir).glob("res_*.h5") if re.fullmatch(r"res_\d+\.h5", p.name)]
    return max(rs, key=lambda p: int(p.stem.split("_")[1])) if rs else None


def _approval_check(sd, card, problem) -> dict:
    from ...intake import approval as AP
    return AP.check(sd["approval"], card, problem)


def _card_precheck(sd) -> dict:
    """最初の全格子の生成の前のカードの照合 (承認の原本か試行の入力か・カードだけで決まる保護の条件; plan tooling-design-orchestrator §4.6
    「照合の場所と順」の 1)。問題の段は wt_card (estimate_mesh が全格子を作る) の前、格子の段は hand_mesh.py build の前に呼ぶ。"""
    from ...intake import approval as AP
    try:
        return AP.check_card(sd["approval"], sd["card"])
    except (AP.ApprovalError, OSError, ValueError) as e:
        return {"status": "MISMATCH", "problems": [f"承認の記録を読めない ({e})"]}


def _limits_precheck(sd, problem) -> list:
    """準備の直前の節点数の上限の照合 (§4.6「照合の場所と順」の 3): 実際に使う問題 YAML の mesh・mesh_euler の上限 ≤ 承認の mesh_limits。
    入力の sha256 の変更を note で受け入れても、ここは省かない。問題の列 (空なら OK)。"""
    from ...intake import approval as AP
    try:
        lim = (AP.load_approval(sd["approval"]).get("protected") or {}).get("mesh_limits") or {}
        return AP.mesh_limit_problems(lim, yaml.safe_load(Path(problem).read_text()))
    except (AP.ApprovalError, OSError, ValueError) as e:
        return [f"承認の記録・問題 YAML を読めない ({e})"]


def _card_with_coefficients(src_card, coef_yaml, approval_path, out):
    """係数の断片 (contur_overrides.yaml、`deltastar_initializer.*` の行) をカードの overrides に入れ、承認の記録への参照を足す (手順書 §7.2)。"""
    from ...intake import approval as AP
    d = yaml.safe_load(Path(src_card).read_text())
    coef = yaml.safe_load(Path(coef_yaml).read_text()) or {}
    ov = dict(d.get("overrides") or {})
    for k in [k for k in ov if str(k).startswith("deltastar_initializer.")]:
        del ov[k]
    ov.update(coef)
    d["overrides"] = ov
    d["approval"] = {"file": str(Path(approval_path).resolve()), "sha256": AP.sha256_file(approval_path)}
    if Path(out).exists():
        raise FileExistsError(f"{out} が既にある")
    Path(out).write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))


def card_initializer(card_path) -> dict:
    """カードの overrides の `deltastar_initializer.*` (model を除く較正の係数) → {係数の名前: 値}。"""
    return {k.split(".", 1)[1]: v for k, v in (yaml.safe_load(Path(card_path).read_text()).get("overrides") or {}).items()
            if str(k).startswith("deltastar_initializer.") and k != "deltastar_initializer.model"}


def initializer_mismatch(card_path, prepare_info_path) -> list:
    """カードの較正の係数が準備の記録の initializer に入っているか (型・有限・値; `plan.initializer_problems`)。"""
    from ..plan import initializer_problems
    pi = Path(prepare_info_path)
    if not pi.is_file():
        return ["prepare_info.json が無い"]
    return initializer_problems(card_initializer(card_path), json.loads(pi.read_text()).get("initializer") or {})


def exec_step(sd: dict) -> dict:
    py = sys.executable
    case = Path(sd["case_dir"])
    repo = Path(sd["repo"])
    design = repo / "design"
    kind = sd["step"].split("_p")[0]
    res = {"step": sd["step"], "ok": False, "branch": None, "detail": None, "verdicts": [], "outputs": {}, "inputs": {}}
    gate = [py, "-m", "forge_design.feedback.gate_env"]
    if kind == "problem":
        out = sd["problem_card"]
        pre = _card_precheck(sd)
        if pre["status"] != "OK":
            res.update(detail=f"カードが承認と結び付かない — wt_card を回さない: {'; '.join(pre['problems'])[:400]}")
            return res
        r = _run([py, "-m", "forge_design.intake.wt_card", sd["card"], "--out", out], "wt_card", design, sd)
        ck = _approval_check(sd, sd["card"], out) if r.returncode == 0 else {"status": "未"}
        res.update(ok=r.returncode == 0 and ck["status"] == "OK", detail=f"wt_card rc {r.returncode}・承認の照合 {ck['status']} {ck.get('diffs', [])[:3]}")
        res["outputs"] = {p: sha256(p) for p in (out, out.replace("_card.yaml", "_card_euler.yaml")) if Path(p).is_file()}
    elif kind == "mesh":
        script = repo / MESH_SCRIPT_REL
        md = case / f"_mesh_p{sd['k']}"
        pre = _card_precheck(sd)
        if pre["status"] != "OK":
            res.update(detail=f"カードが承認と結び付かない — hand_mesh を回さない: {'; '.join(pre['problems'])[:400]}")
            return res
        md.mkdir(parents=True, exist_ok=True)
        r1 = _run([py, script, "build", "--design", design, "--card", sd["card"], "--problem", sd["problem_card"], "--out-yaml", sd["problem"],
                   "--record", md / "build_record.json", "--plan-item", "orchestrator recipe v1", "--axis-cap-frac", str(AXIS_CAP_FRAC),
                   "--approval", sd["approval"]],                  # 節点数の上限は承認の mesh_limits (plan §4.6・§5.1 #8)
                  "hand_mesh_build", design, sd)
        r2 = _run([py, script, "check", "--design", design, "--card", sd["card"], "--problem", sd["problem"], "--record", md / "check_record.json",
                   "--run-dir", md / "check"], "hand_mesh_check", design, sd) if r1.returncode == 0 else None
        npass = 0 if r2 is None else len(re.findall(r"^VERDICT: PASS", r2.stdout, flags=re.M))
        ck = _approval_check(sd, sd["card"], sd["problem"]) if r2 is not None and r2.returncode == 0 else {"status": "未"}
        if ck["status"] == "OK" and Path(sd["problem_euler"]).is_file():
            ce = _approval_check(sd, sd["card"], sd["problem_euler"])          # Euler の写しも (上限は mesh_euler)
            ck = ce if ce["status"] != "OK" else ck
        res.update(ok=r1.returncode == 0 and r2 is not None and r2.returncode == 0 and npass >= 2 and ck["status"] == "OK",
                   detail=f"hand_mesh build rc {r1.returncode}・check rc {None if r2 is None else r2.returncode}・VERDICT PASS {npass}・承認の照合 {ck['status']}")
        res["outputs"] = {p: sha256(p) for p in (sd["problem"], sd["problem_euler"]) if Path(p).is_file()}
    elif kind == "euler":
        rd = Path(sd["run_dir"])
        lp = _limits_precheck(sd, sd["problem_euler"])
        if lp:
            res.update(detail=f"使う問題 YAML の節点数の上限が承認と合わない — 準備しない: {'; '.join(lp)[:400]}")
            return res
        r = _run([py, "-m", "forge_design.evaluate.runner_axismach", sd["problem_euler"], rd], "euler", design, sd)
        last = _last_res(rd)
        mq = (rd / "MESH_QUALITY.txt").read_text().splitlines()[-1] if (rd / "MESH_QUALITY.txt").is_file() else ""
        bad = {k: v for k, v in (_nonfinite(last) if last else {}).items() if v}
        res.update(ok=r.returncode == 0 and last is not None and mq.startswith("VERDICT: PASS") and not bad,
                   detail=f"runner rc {r.returncode}・{mq}・最後の場 {None if last is None else last.name}・非有限 {bad or 0}")
        if last is not None:
            res["euler_ref"] = {"run": str(rd), "res": str(last), "sha256": sha256(last)}
            res["outputs"] = {str(last): res["euler_ref"]["sha256"]}
    elif kind == "ns":
        from .. import ic as IC
        from .. import plan as EP
        rd = Path(sd["run_dir"])
        lp = _limits_precheck(sd, sd["problem"])
        if lp:
            res.update(detail=f"使う問題 YAML の節点数の上限が承認と合わない — 準備しない: {'; '.join(lp)[:400]}")
            return res
        ic = None
        if sd.get("ic_file"):
            ic = IC.resolve_ic(sd["ic_file"]["path"], sd["problem"], for_ns=True)
            if ic["snapshot"]["sha256"] != sd["ic_file"]["sha256"]:
                raise ValueError("初期場の sha256 が前の段の記録と違う")
        # 初期壁の係数は計画の確定 (準備の後・起動の前) で照合する (合わなければ forge を起動しない; codex 2026-10-17 v1a-fail)
        pin = EP.fixed_inputs(card=sd["card"], problem=sd["problem"], approval=sd["approval"], ic=ic, euler_ref=sd["euler_ref"]["res"],
                              expected_initializer=card_initializer(sd["card"]))
        pin_f = Path(sd["out"]).parent / "plan_inputs.json"
        pin_f.write_text(json.dumps(pin, ensure_ascii=False, default=str))
        er = sd["euler_ref"]
        cmd = [py, "-m", "forge_design.feedback.deltastar_loop", "--problem", sd["problem"], "--euler-ref", er["run"],
               "--euler-ref-res", Path(er["res"]).name, "--euler-ref-sha256", er["sha256"], "--init-integral", "--run-dir", rd,
               "--plan-inputs", pin_f]
        if ic is not None:
            cmd += ["--ic-file", ic["snapshot"]["path"], "--ic-sha256", ic["snapshot"]["sha256"]]
        env_cmd = cmd
        r = _run(env_cmd, "ns", design, sd)
        sr = json.loads((rd / "STOP_RULE_STATE.json").read_text()) if (rd / "STOP_RULE_STATE.json").is_file() else {}
        hg = json.loads((rd / "HANDOFF_GATE.json").read_text()) if (rd / "HANDOFF_GATE.json").is_file() else {}
        sv = json.loads((rd / "SETTLE_VERDICT.json").read_text()) if (rd / "SETTLE_VERDICT.json").is_file() else {}
        last = _last_res(rd)
        bad = {k: v for k, v in (_nonfinite(last) if last else {}).items() if v}
        vp = EP.verify_after(rd) if (rd / EP.PLAN_FILE).is_file() else {"status": "計画なし"}
        # 初期壁の係数: カードの overrides の deltastar_initializer.* が準備の記録の initializer に入っている (旧の台本 p5_ns.py の検査; §5.1 #6d)
        ini_bad = initializer_mismatch(sd["card"], rd / "prepare_info.json")
        ok = (r.returncode == 0 and sr.get("status") == "REACHED" and hg.get("proceed") is True and sv.get("verdict") == "SETTLED"
              and not bad and vp["status"] == "OK" and not ini_bad)
        res.update(ok=ok, detail=f"deltastar_loop rc {r.returncode}・停止規則 {sr.get('status')}・受け渡しの関門 進む {hg.get('proceed')}・"
                                f"整定 {sv.get('verdict')}・非有限 {bad or 0}・計画の照合 {vp['status']}・初期壁の係数の欠け・違い {ini_bad or 0}")
        res["run_dir"] = str(rd)
        if last is not None:
            res["final_field"] = {"path": str(last), "sha256": sha256(last)}
            res["outputs"] = {str(last): res["final_field"]["sha256"]}
    elif kind in ("calibrate", "recalibrate"):
        ns = Path(sd["ns_run"])
        if kind == "calibrate":
            r = _run([py, "-m", "forge_design.feedback.contur_calibrate", ns], "contur_calibrate", design, sd)
            coef = ns / "contur_calibrate" / "contur_overrides.yaml"
            rc = r.returncode
        else:
            coef = ns / "wt_qualify" / "contur_overrides.yaml"
            rc = 0
        if rc == 0 and coef.is_file():
            _card_with_coefficients(sd["card"], coef, sd["approval"], sd["new_card"])
            ck = _approval_check(sd, sd["new_card"], sd["problem"])
        else:
            ck = {"status": "未"}
        res.update(ok=rc == 0 and coef.is_file() and ck["status"] == "OK",
                   detail=f"較正 rc {rc}・係数 {coef.name if coef.is_file() else '無い'}・新しいカードの承認の照合 {ck['status']} {ck.get('diffs', [])[:3]}")
        res["outputs"] = {sd["new_card"]: sha256(sd["new_card"])} if Path(sd["new_card"]).is_file() else {}
    elif kind in ("qualify", "deliverables"):
        ns = Path(sd["ns_run"])
        if kind == "deliverables":
            rd_ = _run([*gate, "deliverables", "--approval", sd["approval"], "--", ns, "--card", sd["card"], "--problem", sd["problem"]],
                       "deliverables", design, sd)
            if rd_.returncode != 0:
                res.update(ok=False, detail=f"成果物 rc {rd_.returncode}")
                Path(sd["out"]).write_text(json.dumps(res, ensure_ascii=False, indent=1))
                return res
        od = ns / ("wt_qualify" if kind == "qualify" else "wt_qualify_final")
        args = [ns, "--card", sd["card"], "--problem", sd["problem"], "--iteration", str(sd["iteration"]),
                "--recalibrations-done", str(sd["recalibrations_done"]), "--out-dir", od]
        if sd.get("prev_run"):
            args += ["--prev-run", sd["prev_run"]]
        r = _run([*gate, "judge", "--approval", sd["approval"], "--", *args], "judge", design, sd)
        rec = json.loads((od / "WT_QUALIFICATION.json").read_text()) if (od / "WT_QUALIFICATION.json").is_file() else {}
        result, nxt = rec.get("result"), rec.get("next")
        branch = ("accepted" if result == "accepted" else "make_deliverables" if (result, nxt) == ("incomplete", "make_deliverables")
                  else "recalibrate" if (result, nxt) == ("back_to_correction", "recalibrate") else "other")
        res.update(ok=bool(rec) and branch != "other" and (kind == "qualify" or branch == "accepted"), branch=branch,
                   detail=f"判定 rc {r.returncode}・{result}・{nxt}・欠け {rec.get('missing')}", verdicts=[{"tool": "wt_qualify", "result": result, "next": nxt}])
    else:
        raise ValueError(f"知らない段 {sd['step']}")
    Path(sd["out"]).write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))
    return res


RECIPE = Recipe()


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2 or argv[0] != "exec":
        print("使い方: python -m forge_design.campaign.recipes.wt_card_chain exec '<段の定義の JSON>'", file=sys.stderr)
        return 2
    r = exec_step(json.loads(argv[1]))
    print(json.dumps({k: r[k] for k in ("step", "ok", "branch", "detail")}, ensure_ascii=False))
    return 0 if r["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
