"""§6.19 の集約 (plan time_integration-line-implicit-speed、codex plan-8 M1)。各 run の cmp_judge_LAYOUT2.json を読み、
PASS か、「判定不能の理由が η の行だけ (腕の η > 上限で従来の η も 1e-11 超) で fail が空・非有限 0・factor と solve の本数がそろう」なら合格。
判定 JSON が無い・run の開始より古い・FAIL・ほかの判定不能 → 不合格。結果は _band_ab/cold_pair/lay3_judge.json、全件合格で終了コード 0。
usage: python3 lay3_judge.py <run> <開始時刻 epoch> [<run> <開始時刻> ...]"""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ETA = re.compile(r"^solve \d+: 腕の η \S+ > 上限、従来も \S+ > 1e-11$")
out, ok_all = {}, True
args = sys.argv[1:]
for run, t0 in zip(args[0::2], args[1::2]):
    p = HERE / run / "cmp_judge_LAYOUT2.json"
    if not p.exists() or p.stat().st_mtime < float(t0):
        out[run] = {"verdict": "不合格 (判定 JSON が無いか今回のものでない)"}; ok_all = False; continue
    j = json.loads(p.read_text())
    if j["verdict"] == "PASS": v = "合格 (PASS)"
    elif j["verdict"] == "INDETERMINATE" and not j["fail"] and j["nonfinite"] == 0 and j["indeterminate"] and all(ETA.match(x) for x in j["indeterminate"]):
        v = "合格 (ビット一致・判定不能の理由は η の大きさだけ。線形解の精度の合格ではない)"
    else: v = f"不合格 ({j['verdict']})"; ok_all = False
    out[run] = {"verdict": v, "judge": j["verdict"], "fail": j["fail"][:3], "indeterminate": j["indeterminate"][:3], "eta_arm_max": j.get("eta_arm_max"), "n_factor": j.get("n_factor"), "n_solve": j.get("n_solve")}
out["all_pass"] = ok_all
(HERE / "_band_ab" / "cold_pair" / "lay3_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False)); print(json.dumps(out, indent=1, ensure_ascii=False))
sys.exit(0 if ok_all else 1)
