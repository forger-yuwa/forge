#!/usr/bin/env python3
"""recipe `axismach.wt_card_chain/v1` の段と分岐 (`forge_design.campaign.recipes.wt_card_chain.Recipe.next_step`) と README の同期の試験
(plan tooling-design-orchestrator §4.5.1、§5.1 #6c)。GPU・道具の実行なし: 台帳と各段の結果のファイル (result.json) を作って next_step の答えを見る。

  (R1) 空 → mesh_p0、p0 の段の後 → calibrate_p0 (最初の較正)、その後 → problem_p1
  (R2) qualify_p1 が make_deliverables → deliverables_p1、deliverables_p1 が accepted → 終わり (None)
  (R3) qualify_p1 が accepted (成果物が既にある) → 終わり
  (R4) qualify_p1 が recalibrate → recalibrate_p1 (最初の較正は較正し直しの回数に入らない)、その後 → problem_p2
  (R5) 較正し直しの上限 1: qualify_p2 が recalibrate → 止まる
  (R6) 反復の上限 2: recalibrate_p1 の後 → 止まる
  (R7) qualify が escalate など → 止まる、最後の試行が合格しなかった → 止まる
  (R8) 前の段が作った入力 (card_p1) が変わる → 止まる。変更の要求 (note: 新しい sha256 と理由) があれば通る
  (R9) NS の段の定義: k ≥ 1 は前の NS の最後の場を初期場に、Euler の参照を固定、run の名前は次の空き番号
  (I1) 初期壁の係数の照合 (試走 A の NS p1; 係数を 1 % 変えたカードで違いを出す)
  (M1) README の同期: 表の後ろに印つきの行、ほかの行は残る、2 回目は印の間だけ置き換える

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_recipe_tests.py
"""
import copy
import hashlib
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.campaign import readme as RM  # noqa: E402
from forge_design.campaign.recipes import wt_card_chain as WC  # noqa: E402

FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


with tempfile.TemporaryDirectory() as tdn:
    case = Path(tdn) / "case"
    led = case / "_ledger"
    (led / "attempts").mkdir(parents=True)
    (case / "card_p0.yaml").write_text("x: 0\n")
    (case / "problem_p0_card.yaml").write_text("x: 0\n")

    def campaign(max_it=3, max_recal=1):
        apr = case / f"apr_{max_it}_{max_recal}.json"
        apr.write_text(json.dumps({"card": {"path": str(case / "card_p0.yaml")}, "problem": {"path": str(case / "problem_p0_card.yaml")},
                                   "protected": {"run_policy": {"max_iterations": max_it, "max_recalibrations": max_recal}}}))
        return {"case_dir": str(case), "ledger": str(led), "approval": {"file": str(apr), "sha256": "x"}, "name": "t"}

    def doc_with(steps):
        """steps: [(step, result dict)] を順に取り込んだ台帳。"""
        doc = {"attempts": [], "ingested": [], "requests": [], "run_env": {}}
        for i, (st, r) in enumerate(steps, 1):
            aid = f"{i:04d}-{st}"
            ad = led / "attempts" / aid
            ad.mkdir(parents=True, exist_ok=True)
            (ad / "result.json").write_text(json.dumps({"step": st, **r}))
            doc["attempts"].append({"id": aid, "step": st, "run_dir": r.get("run_dir")})
            doc["ingested"].append(aid)
        return doc

    R = WC.Recipe(repo_dir=Path(__file__).resolve().parents[2])
    OK = {"ok": True}
    eu0 = {"ok": True, "euler_ref": {"run": str(case / "run_0001_euler_p0"), "res": str(case / "run_0001_euler_p0" / "res_12000.h5"), "sha256": None}}
    ns0 = {"ok": True, "run_dir": str(case / "run_0002_ns_p0"), "final_field": {"path": str(case / "run_0002_ns_p0" / "res_42060.h5"), "sha256": None}}
    eu1 = {"ok": True, "euler_ref": {"run": str(case / "run_0003_euler_p1"), "res": str(case / "run_0003_euler_p1" / "res_12000.h5"), "sha256": None}}
    ns1 = {"ok": True, "run_dir": str(case / "run_0004_ns_p1"), "final_field": {"path": str(case / "run_0004_ns_p1" / "res_30088.h5"), "sha256": None}}
    p0 = [("mesh_p0", OK), ("euler_p0", eu0), ("ns_p0", ns0)]
    (case / "card_p1.yaml").write_text("x: 1\n")
    cal0 = {"ok": True, "outputs": {str(case / "card_p1.yaml"): hashlib.sha256(b"x: 1\n").hexdigest()}}
    p1 = p0 + [("calibrate_p0", cal0), ("problem_p1", OK), ("mesh_p1", OK), ("euler_p1", eu1), ("ns_p1", ns1)]
    C = campaign()

    def nxt(steps, camp=C):
        r = R.next_step(doc_with(steps), camp)
        return None if r is None else r.get("step", "STOP: " + str(r.get("stop")))

    check("(R1) 空 → mesh_p0", nxt([]) == "mesh_p0")
    check("(R1) p0 の段の後 → calibrate_p0", nxt(p0) == "calibrate_p0", nxt(p0))
    check("(R1) 最初の較正の後 → problem_p1", nxt(p0 + [("calibrate_p0", cal0)]) == "problem_p1")
    q = lambda br: ("qualify_p1", {"ok": br != "other", "branch": br, "detail": br})  # noqa: E731
    check("(R2) make_deliverables → deliverables_p1", nxt(p1 + [q("make_deliverables")]) == "deliverables_p1")
    check("(R2) deliverables_p1 が accepted → 終わり", nxt(p1 + [q("make_deliverables"), ("deliverables_p1", {"ok": True, "branch": "accepted"})]) is None)
    check("(R3) accepted (成果物が既にある) → 終わり", nxt(p1 + [q("accepted")]) is None)
    check("(R4) recalibrate → recalibrate_p1 (最初の較正は回数に入らない)", nxt(p1 + [q("recalibrate")]) == "recalibrate_p1")
    (case / "card_p2.yaml").write_text("x: 2\n")
    rc1 = ("recalibrate_p1", {"ok": True, "outputs": {str(case / "card_p2.yaml"): hashlib.sha256(b"x: 2\n").hexdigest()}})
    check("(R4) 較正し直しの後 → problem_p2", nxt(p1 + [q("recalibrate"), rc1]) == "problem_p2")
    p2 = p1 + [q("recalibrate"), rc1, ("problem_p2", OK), ("mesh_p2", OK), ("euler_p2", eu1), ("ns_p2", ns1),
               ("qualify_p2", {"ok": True, "branch": "recalibrate"})]
    s = nxt(p2)
    check("(R5) 較正し直しの上限 1: qualify_p2 が recalibrate → 止まる", s.startswith("STOP") and "較正し直しの上限" in s, s)
    s = nxt(p1 + [q("recalibrate"), rc1], campaign(max_it=2))
    check("(R6) 反復の上限 2: recalibrate_p1 の後 → 止まる", s.startswith("STOP") and "反復の上限" in s, s)
    s = nxt(p1 + [("qualify_p1", {"ok": True, "branch": "other", "detail": "back_to_correction・escalate"})])
    check("(R7) qualify が escalate → 止まる", s.startswith("STOP"), s)
    s = nxt(p0 + [("calibrate_p0", {"ok": False, "detail": "較正 rc 1"})])
    check("(R7) 最後の試行が合格しなかった → 止まる", s.startswith("STOP") and "合格しなかった" in s, s)
    (case / "card_p1.yaml").write_text("x: 1  # 手で直した\n")
    s = nxt(p0 + [("calibrate_p0", cal0)])
    check("(R8) 前の段が作った入力 (card_p1) が変わる → 止まる", s.startswith("STOP") and "card_p1" in s and "変わった" in s, s)
    d = doc_with(p0 + [("calibrate_p0", cal0)])
    d["requests"].append({"kind": "note", "path": str(case / "card_p1.yaml"), "input_sha256": hashlib.sha256(b"x: 1  # \xe6\x89\x8b\xe3\x81\xa7\xe7\x9b\xb4\xe3\x81\x97\xe3\x81\x9f\n").hexdigest(),
                          "reason": "コメントを足した"})
    r = R.next_step(d, C)
    check("(R8) 変更の要求 (新しい sha256 と理由) があれば通る", r is not None and r.get("step") == "problem_p1", str(r)[:200])
    (case / "card_p1.yaml").write_text("x: 1\n")
    (case / "run_0007_old").mkdir()
    r = R.next_step(doc_with(p1[:-1]), C)
    sd = json.loads(r["command"][-1])
    check("(R9) ns_p1: 前の NS の最後の場を初期場に、Euler の参照を固定、run の名前は次の空き番号",
          r["step"] == "ns_p1" and sd["ic_file"]["path"].endswith("res_42060.h5") and sd["euler_ref"]["res"].endswith("res_12000.h5")
          and Path(r["run_dir"]).name == "run_0008_ns_p1" and r["gpu"] is True, f"{r['run_dir']}")

    print("== 初期壁の係数の照合 (旧の台本 p5_ns.py の検査)")
    QA = Path(__file__).resolve().parents[2] / "case" / "90.qs_trial_a"
    if (QA / "run_0005_ns_p1" / "prepare_info.json").is_file() and (QA / "card_p1.yaml").is_file():
        ok1 = WC.initializer_mismatch(QA / "card_p1.yaml", QA / "run_0005_ns_p1" / "prepare_info.json")
        bad_card = case / "card_bad_coef.yaml"
        import yaml as _y
        dd = _y.safe_load((QA / "card_p1.yaml").read_text())
        dd["overrides"]["deltastar_initializer.delta_scale"] = float(dd["overrides"]["deltastar_initializer.delta_scale"]) * 1.01
        bad_card.write_text(_y.safe_dump(dd, sort_keys=False, allow_unicode=True))
        bad1 = WC.initializer_mismatch(bad_card, QA / "run_0005_ns_p1" / "prepare_info.json")
        check("(I1) 試走 A の NS p1: カードの 3 係数が準備の記録にある、係数を 1 % 変えたカードでは違いを出す", ok1 == [] and len(bad1) == 1 and bad1[0].startswith("delta_scale"), f"{ok1} / {bad1}")
    else:
        print("飛ばした (I1): 試走 A が無い")

    print("== README の同期")
    (case / "README.md").write_text("# case\n\n概要\n\n## 計算 run 一覧\n\n| `run_*` | 目的 | 結果 | 状態 |\n| --- | --- | --- | --- |\n"
                                    "| `run_0001_manual` | 人が書いた行 | x | ref |\n\n## 別の節\n\n本文\n")
    ld = {"campaign_id": "cid-1", "name": "t", "recipe": {"name": WC.NAME, "version": WC.VERSION}, "ingested": ["0001-euler_p0"],
          "attempts": [{"id": "0001-euler_p0", "step": "euler_p0", "run_dir": str(case / "run_0001_euler_p0"), "state": "ingested",
                        "interpretation": {"ok": True, "detail": "runner rc 0"}}]}
    RM.sync_runs(case, ld)
    t1 = (case / "README.md").read_text()
    ld2 = copy.deepcopy(ld)
    ld2["attempts"].append({"id": "0002-ns_p0", "step": "ns_p0", "run_dir": str(case / "run_0002_ns_p0"), "state": "launched"})
    RM.sync_runs(case, ld2)
    t2 = (case / "README.md").read_text()
    tab = t2.split("## 計算 run 一覧")[1].split("## 別の節")[0]
    check("(M1) 表の後ろに印つきの行、人が書いた行と別の節は残る",
          "run_0001_manual" in tab and "run_0001_euler_p0" in tab and "<!-- forge-design:cid-1 begin -->" in tab and "## 別の節\n\n本文" in t2)
    check("(M1) 2 回目は印の間だけ置き換える (行が増え、重複しない)",
          t2.count("run_0001_euler_p0") == 1 and "run_0002_ns_p0" in t2 and t2.count("forge-design:cid-1 begin") == 1
          and t1.split("<!-- forge-design")[0] == t2.split("<!-- forge-design")[0])
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
