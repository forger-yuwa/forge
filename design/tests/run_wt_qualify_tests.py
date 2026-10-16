"""風洞ノズルの設計の判定 `feedback/wt_qualify.py` の試験 (plan tooling-nozzle-card-local-mesh §5.1 #6b・§4.3・§6 V5)。

合成のデータ:
 (1) 公差・窓の出所の順 (CLI → カードの qualification → 問題 YAML の qualification → 問題 YAML の evaluate.settle → run の
     settle_criteria.yaml → 既定 (preset wt_nozzle_v2 の線・plan §6 V5b・V5c))。知らないキーは ValueError。
 (2) V5b・V5c を合成の構造格子の場で: 線形の M で最大偏差・片振幅・出口コア M が式どおり、合否、軸の節点の重み > 1e−6・窓が格子の外は
     判定不能、比較の線 (η 0.1・軸) を併記。
 (3) 計算の成立 (停止規則の状態): v2 の REACHED + SETTLED → PASS、CENSORED・再評価が SETTLED でない・v2 の基準があるのに状態が無い → FAIL、
     状態が無い旧い run → UNVERIFIED、v1 の REACHED → PASS。
 (4) 成果物と寸法の関門: そろって PASS、ファイルの改変・部品の欠け・寸法の FAIL・別のカード → FAIL、無ければ ABSENT、無くてよい部品の指定。
 (5) 結果と次の手 (`decide`): accepted・成果物待ち (make_deliverables)・計算の成立を確かめていない・判定不能・反復の上限 (rejected)・
     壁の形だけの補正 (規則で recalibrate / escalate)・壁解像や寸法を含む補正 (escalate、規則は評価しない)。
 (6) 較正し直しの規則 (合成の run): (a)〜(e) をそろえて recalibrate、1 つずつ欠くと escalate (前の run が無い・前の更新で誤差が増えた・
     さらに減らせない・回数・計算の成立)。
実データ (手元にあるとき):
 (7) case/71 run_0010・run_0011: V5b FAIL (出口コア M 5.502841 / 5.502842)・V5c FAIL (η 0.03 片振幅 0.0573 %、最大偏差 0.0623 % は許容内)・
     壁解像 PASS、較正し直しの規則に当たらない (形の誤差の射影 0.624 → 0.612 mm、1.9 % しか減らない) → back_to_correction・escalate
     (plan §4.3 の 2026-10-15 の判断どおり)。run_0011 は環境変数 FORGE_CASE71_RUN_0011 か case/71/run_0011_ns_v2 にあるとき。
使い方: design/.venv-opt/bin/python design/tests/run_wt_qualify_tests.py  (終了コード 0 全合格 / 1 不合格 / 2 実データを飛ばした)
"""
import hashlib
import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from forge_design.feedback import wt_qualify as WQ  # noqa: E402
import contur_calibrate_fixture as FX  # noqa: E402

FAIL, SKIPPED = [], []
MD = 5.5


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------ (1) 公差・窓の出所
print("== (1) 公差・窓の出所")
info = {"scale_m": 0.1, "x_E": 40.0, "physical_wall": {"domain": [-10.0, 80.0]}, "Md": MD}
with tempfile.TemporaryDirectory() as td:
    run = Path(td)
    t = WQ.resolve_tolerances(run, info, None, None, {})
    # V5b・V5c の片振幅の公差は 2026-10-16 のユーザ決定で 0.05 % → 0.1 % (plan tooling-nozzle-card-local-mesh §6 V5b)
    check("(1) 既定: 線 η 0.03・x = (x_E + 2, x_F − 1) r_t × S = [4.2, 7.9] m・2001 点、V5b 1e−3、V5c 0.1・0.1",
          t["line_eta"]["value"] == 0.03 and np.allclose(t["line_x_range_m"]["value"], [4.2, 7.9]) and t["line_n"]["value"] == 2001
          and t["exit_core_M_tol_rel"]["value"] == 1e-3 and t["line_tol_maxabs_pct"]["value"] == 0.1
          and t["line_tol_halfamp_pct"]["value"] == 0.1 and all("既定" in t[k]["source"] or "plan §6" in t[k]["source"] for k in WQ.QUAL_KEYS),
          f"{t['line_x_range_m']}")
    (run / "settle_criteria.yaml").write_text(yaml.safe_dump({"quantities": [
        {"kind": "line_mach", "eta": 0.04, "x_range_m": [4.0, 7.95], "n": 1001, "Md": MD}]}))
    t = WQ.resolve_tolerances(run, info, None, None, {})
    check("(1) run の settle_criteria.yaml の line_mach の η・窓・点数", t["line_eta"]["value"] == 0.04 and t["line_x_range_m"]["value"] == [4.0, 7.95]
          and t["line_n"]["value"] == 1001 and "settle_criteria" in t["line_x_range_m"]["source"])
    prob = {"evaluate": {"settle": {"preset": "wt_nozzle_v2", "x_range_m": [4.1, 7.8], "line_eta": 0.05}}}
    t = WQ.resolve_tolerances(run, info, None, prob, {})
    check("(1) 問題 YAML の evaluate.settle が run の記録より先", t["line_x_range_m"]["value"] == [4.1, 7.8] and t["line_eta"]["value"] == 0.05
          and "evaluate.settle" in t["line_x_range_m"]["source"])
    prob["qualification"] = {"line_x_range_m": [4.3, 7.7], "exit_core_M_tol_rel": 1e-3}
    t = WQ.resolve_tolerances(run, info, None, prob, {})
    check("(1) 問題 YAML の qualification が evaluate.settle より先", t["line_x_range_m"]["value"] == [4.3, 7.7]
          and t["exit_core_M_tol_rel"]["value"] == 1e-3 and t["exit_core_M_tol_rel"]["source"] == "問題 YAML の qualification")
    card = {"M_exit": MD, "qualification": {"line_x_range_m": [4.4, 7.6], "line_tol_halfamp_pct": 0.06}}
    t = WQ.resolve_tolerances(run, info, card, prob, {})
    check("(1) カードの qualification が問題 YAML より先", t["line_x_range_m"]["value"] == [4.4, 7.6] and t["line_tol_halfamp_pct"]["value"] == 0.06
          and t["line_x_range_m"]["source"] == "カードの qualification")
    t = WQ.resolve_tolerances(run, info, card, prob, {"line_x_range_m": [4.5, 7.5], "line_eta": 0.02})
    check("(1) CLI が最優先", t["line_x_range_m"]["value"] == [4.5, 7.5] and t["line_eta"]["value"] == 0.02 and t["line_eta"]["source"] == "CLI")
    for nm, bad in (("知らないキー", {"qualification": {"halfamp": 0.05}}), ("窓が 2 要素でない", {"qualification": {"line_x_range_m": [4.0]}})):
        try:
            WQ.resolve_tolerances(run, info, bad, None, {})
            check(f"(1) qualification の{nm} → ValueError", False, "例外が出ない")
        except ValueError:
            check(f"(1) qualification の{nm} → ValueError", True)
    (run / "settle_criteria.yaml").unlink()
    t = WQ.resolve_tolerances(run, {"scale_m": 0.1}, None, None, {})
    check("(1) 既定の窓を作れない (x_E が無い) ときは None と理由", t["line_x_range_m"]["value"] is None and "作れない" in t["line_x_range_m"]["source"])


# ------------------------------------------------------------------ (2) V5b・V5c
print("== (2) V5b・V5c (合成の場)")


def field(a_pct=-0.02, b_pct_per_m=0.02, eta1=0.03, S=0.1):
    """x_rt −10〜80 (401 列)、行の η = [0, η₁, 0.06, 0.1, 0.2, …, 1]、M = M_d (1 + (a + b (x_m − 4))/100) (η によらない)。"""
    xs = np.linspace(-10.0, 80.0, 401)
    etas = np.r_[0.0, eta1, 0.06, 0.1, np.linspace(0.2, 1.0, 9)]
    Rw = 1.0 + 0.08 * xs
    X = np.repeat(xs[:, None], len(etas), axis=1)
    R = Rw[:, None] * etas[None, :]
    d = a_pct + b_pct_per_m * (X * S - 4.0)
    M = MD * (1.0 + d / 100.0)
    sonic = np.full_like(M, 300.0)
    V = {"M": M, "Ux": M * sonic, "Uy": np.zeros_like(M), "sonic": sonic, "ro": np.ones_like(M)}
    return {"X": X, "R": R, "V": V, "S": S}


F = field()
v = WQ.v5c_line(F, MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.05)
want_max, want_half = 0.02 * 3.95 - 0.02, 0.02 * 3.95 / 2
check("(2) V5c: 線形の M で最大偏差・片振幅が式どおり (1e−9) で PASS",
      abs(v["maxabs_pct"] - want_max) < 1e-9 and abs(v["halfamp_pct"] - want_half) < 1e-9 and v["verdict"] == "PASS",
      f"{v['maxabs_pct']:.6f} / {v['halfamp_pct']:.6f} %")
check("(2) V5c: 比較の線 η 0.1・軸を併記 (判定に使わない)", set(v["compare"]) == {"eta0.1", "eta0"}
      and abs(v["compare"]["eta0.1"]["halfamp_pct"] - want_half) < 1e-9)
check("(2) V5c: 第 1 層 η 0.03 の格子では軸の重み ≤ 1e−6", v["axis_weight_max"] <= 1e-6, f"{v['axis_weight_max']:.1e}")
v = WQ.v5c_line(field(b_pct_per_m=0.04), MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.05)
check("(2) V5c: 傾きを倍にすると両量とも FAIL", v["verdict"] == "FAIL" and not v["checks"]["maxabs"] and not v["checks"]["halfamp"])
v = WQ.v5c_line(field(a_pct=-0.06, b_pct_per_m=0.03), MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.05)
check("(2) V5c: 片振幅だけ超える (0.059 %) → FAIL、最大偏差 (0.06 %) は許容内", v["verdict"] == "FAIL" and v["checks"]["maxabs"]
      and not v["checks"]["halfamp"], f"{v['maxabs_pct']:.4f} / {v['halfamp_pct']:.4f}")
v = WQ.v5c_line(field(eta1=0.05), MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.05)
check("(2) V5c: 第 1 層の η 0.05 > 線 0.03 → 軸の節点の重み 0.4 で判定不能", v["verdict"] == "INDETERMINATE" and abs(v["axis_weight_max"] - 0.4) < 1e-12,
      v.get("why"))
v = WQ.v5c_line(F, MD, 0.03, [4.0, 9.0], 2001, 0.1, 0.05)
check("(2) V5c: 窓が格子の外 (x_F 8 m を超える) → 判定不能", v["verdict"] == "INDETERMINATE", v.get("why"))
Fn = field()
Fn["V"]["M"][270, 1] = np.nan          # x_rt 50.75 = 5.075 m (窓の中、線 η 0.03 の行)
v = WQ.v5c_line(Fn, MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.05)
check("(2) V5c: 線の M に非有限 → 判定不能", v["verdict"] == "INDETERMINATE")
b = WQ.v5b_exit_core(F, MD, 5e-4, 40.0)
wantM = MD * (1.0 + (-0.02 + 0.02 * (8.0 - 4.0)) / 100.0)
check("(2) V5b: 出口コア M (最後の列の η ∈ [0.05, 0.7] の平均) が式どおり、+0.06 % → FAIL", abs(b["M"] - wantM) < 1e-12 and b["verdict"] == "FAIL",
      f"{b['M']!r}")
b = WQ.v5b_exit_core(field(a_pct=-0.07), MD, 5e-4, 40.0)
check("(2) V5b: +0.01 % → PASS", b["verdict"] == "PASS" and abs(b["dev_pct"] - 0.01) < 1e-9, f"{b['dev_pct']:.6f} %")

# ------------------------------------------------------------------ (3) 計算の成立
print("== (3) 計算の成立")
with tempfile.TemporaryDirectory() as td:
    td = Path(td)

    def flow(name, state=None, v2=False):
        d = td / name
        d.mkdir()
        if state is not None:
            (d / "STOP_RULE_STATE.json").write_text(json.dumps(state))
        if v2:
            (d / "settle_criteria.yaml").write_text("quantities: []\n")
        return WQ.flow_check(d)["status"]

    check("(3) v2 の REACHED + 再評価 SETTLED → PASS", flow("a", {"method": "v2", "status": "REACHED", "reeval": {"verdict": "SETTLED"}}, True) == "PASS")
    check("(3) v2 の CENSORED → FAIL", flow("b", {"method": "v2", "status": "CENSORED", "reeval": {"verdict": "SETTLED"}}, True) == "FAIL")
    check("(3) v2 の再評価が SETTLED でない → FAIL", flow("c", {"method": "v2", "status": "REACHED", "reeval": {"verdict": "NOT_SETTLED"}}, True) == "FAIL")
    check("(3) v2 の基準があるのに状態が無い → FAIL", flow("d", None, True) == "FAIL")
    check("(3) 状態が無い旧い run → UNVERIFIED", flow("e") == "UNVERIFIED")
    check("(3) v1 の REACHED → PASS・CENSORED → FAIL", flow("f", {"version": 1, "status": "REACHED"}) == "PASS"
          and flow("g", {"version": 1, "status": "CENSORED"}) == "FAIL")

# ------------------------------------------------------------------ (4) 成果物と寸法の関門
print("== (4) 成果物と寸法の関門")
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    run = td / "run_x"
    run.mkdir()
    card = td / "card.yaml"
    card.write_text("M_exit: 5.5\n")
    check("(4) 成果物が無い → ABSENT", WQ.deliverables_check(run, None, card)["status"] == "ABSENT")
    from forge_design.export.deliverables import PARTS
    od = run / "deliverables"
    (od / "cad").mkdir(parents=True)
    (od / "cad" / "wall.step").write_text("STEP")
    dim = {"format": "forge.dimension_gate/v1", "verdict": "PASS", "items": {"L_div": {"pass": True}}, "card": "card.yaml", "card_sha256": sha(card)}
    (od / "DIMENSIONS.json").write_text(json.dumps(dim))
    man = {"run": str(run), "verdict": "PASS", "failed_parts": [], "parts": {p: {"ok": True} for p in PARTS}, "skipped": {}, "errors": {},
           "files_sha256": {"DIMENSIONS.json": sha(od / "DIMENSIONS.json"), "cad/wall.step": sha(od / "cad" / "wall.step")}}
    (od / "MANIFEST.json").write_text(json.dumps(man))
    r = WQ.deliverables_check(run, None, card)
    check("(4) そろって PASS (寸法の関門 PASS・カードの sha256 一致)", r["status"] == "PASS" and r["dimensions"]["status"] == "PASS", f"{r.get('problems')}")
    other = td / "other.yaml"
    other.write_text("M_exit: 6.0\n")
    r = WQ.deliverables_check(run, None, other)
    check("(4) 寸法の関門が別のカードで判定された → 寸法 FAIL", r["dimensions"]["status"] == "FAIL", f"{r['dimensions']['problems']}")
    (od / "cad" / "wall.step").write_text("STEP changed")
    r = WQ.deliverables_check(run, None, card)
    check("(4) 列挙したファイルが改変された → FAIL", r["status"] == "FAIL" and any("sha256" in p for p in r["problems"]), f"{r['problems']}")
    (od / "cad" / "wall.step").write_text("STEP")
    man["parts"].pop("fluent")
    man["skipped"] = {"fluent": "--skip で指定"}
    (od / "MANIFEST.json").write_text(json.dumps(man))
    r = WQ.deliverables_check(run, None, card)
    check("(4) 必要な部品が欠けた → FAIL", r["status"] == "FAIL" and any("fluent" in p for p in r["problems"]))
    r = WQ.deliverables_check(run, None, card, optional_parts=("fluent",))
    check("(4) 無くてよい部品に指定すれば PASS", r["status"] == "PASS", f"{r.get('problems')}")
    dim["verdict"] = "FAIL"
    (od / "DIMENSIONS.json").write_text(json.dumps(dim))
    man["files_sha256"]["DIMENSIONS.json"] = sha(od / "DIMENSIONS.json")
    (od / "MANIFEST.json").write_text(json.dumps(man))
    r = WQ.deliverables_check(run, None, card, optional_parts=("fluent",))
    check("(4) 寸法の関門の FAIL → 寸法 FAIL", r["dimensions"]["status"] == "FAIL")

# ------------------------------------------------------------------ (5) 結果と次の手
print("== (5) 結果と次の手 (decide)")
P, Fa, I = {"verdict": "PASS"}, {"verdict": "FAIL"}, {"verdict": "INDETERMINATE", "why": "x"}
fP, fU, fF = {"status": "PASS"}, {"status": "UNVERIFIED", "why": "u"}, {"status": "FAIL", "why": "f"}
dP = {"status": "PASS", "dimensions": {"status": "PASS"}}
dA = {"status": "ABSENT", "dimensions": {"status": "ABSENT"}}
dDF = {"status": "PASS", "dimensions": {"status": "FAIL"}}
calls = []


def rule(apply):
    def f():
        calls.append(apply)
        return {"apply": apply}
    return f


def dec(*a, **k):
    r = WQ.decide(*a, **k)
    return r["result"], r["next"]


check("(5) すべて合格・成果物 PASS・計算の成立 PASS → accepted / done", dec(fP, P, P, P, dP) == ("accepted", "done"))
check("(5) 流れは合格・成果物が無い → incomplete / make_deliverables", dec(fP, P, P, P, dA) == ("incomplete", "make_deliverables"))
check("(5) 成果物が無く計算の成立も UNVERIFIED → incomplete / escalate", dec(fU, P, P, P, dA) == ("incomplete", "escalate"))
check("(5) すべて合格でも計算の成立が UNVERIFIED → incomplete / escalate (accepted にしない)", dec(fU, P, P, P, dP) == ("incomplete", "escalate"))
check("(5) V5c 判定不能 → incomplete / escalate", dec(fP, P, I, P, dP) == ("incomplete", "escalate"))
check("(5) 計算の成立 FAIL → incomplete / escalate (公差が FAIL でも)", dec(fF, Fa, Fa, P, dP) == ("incomplete", "escalate"))
calls.clear()
check("(5) 壁の形 (V5b) だけ FAIL・規則に当たる → back_to_correction / recalibrate", dec(fP, Fa, P, P, dA, rule_fn=rule(True))
      == ("back_to_correction", "recalibrate") and calls == [True])
check("(5) 壁の形 (V5c) だけ FAIL・規則に当たらない → back_to_correction / escalate", dec(fP, P, Fa, P, dA, rule_fn=rule(False))
      == ("back_to_correction", "escalate"))
calls.clear()
check("(5) 反復の上限 (iteration 3 ≥ 3) → rejected / escalate (規則は評価しない)",
      dec(fP, Fa, Fa, P, dA, iteration=3, max_iterations=3, rule_fn=rule(True)) == ("rejected", "escalate") and calls == [])
check("(5) 上限の手前 (iteration 2 < 3) は back_to_correction", dec(fP, Fa, Fa, P, dA, iteration=2, max_iterations=3, rule_fn=rule(True))
      == ("back_to_correction", "recalibrate"))
calls.clear()
check("(5) 壁解像 FAIL を含む補正 → back_to_correction / escalate (規則は評価しない)",
      dec(fP, Fa, P, Fa, dA, rule_fn=rule(True)) == ("back_to_correction", "escalate") and calls == [])
check("(5) 寸法の関門 FAIL → back_to_correction / escalate", dec(fP, P, P, P, dDF, rule_fn=rule(True)) == ("back_to_correction", "escalate"))
r = WQ.decide(fP, P, P, P, dA)
check("(5) 反復の上限を渡さなければ「判定していない」と記録", r["iteration"]["evaluated"] is False and "判定していない" in r["iteration"]["note"])

# ------------------------------------------------------------------ (6) 較正し直しの規則
print("== (6) 較正し直しの規則 (合成の run)")
tgt = (0.68, 0.021, -3.06e-4)
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    cur = FX.make_run(td / "run_cur", coefs=(0.75, 0.01, -3.1e-4), target=tgt, bump=2e-5, v2=True)
    prev = FX.make_run(td / "run_prev", coefs=(0.85, 0.0, -3.2e-4), target=tgt, bump=2e-5)
    near = FX.make_run(td / "run_near", coefs=(0.68, 0.021, -3.06e-4), target=tgt, bump=2e-4, bump_kind="wave")
    v5b, v5c = {"verdict": "FAIL"}, {"verdict": "FAIL"}
    r = WQ.recalibration_rule(cur, v5b, v5c, fP, prev_run=prev)
    c = r["conditions"]
    check("(6) (a)〜(e) がそろえば recalibrate (新しい係数は目標に近い)", r["apply"] and abs(r["calibration"]["coefficients"]["c1_per_m"] - tgt[1]) < 1e-3,
          f"誤差 {c['c_further_reducible']['current_max_abs_mm']:.4f} → {c['c_further_reducible']['projected_max_abs_mm']:.4f} mm、"
          f"前 {c['b_prev_update_reduced_error']['prev_max_abs_mm']:.4f} mm")
    r = WQ.recalibration_rule(cur, v5b, v5c, fP)
    check("(6) 前の run が無い → (b) で escalate", not r["apply"] and not r["conditions"]["b_prev_update_reduced_error"]["ok"])
    r = WQ.recalibration_rule(prev, v5b, v5c, fP, prev_run=cur)
    check("(6) 前の更新で誤差が増えた → (b) で escalate", not r["apply"] and not r["conditions"]["b_prev_update_reduced_error"]["ok"])
    r = WQ.recalibration_rule(near, v5b, v5c, fP, prev_run=prev)
    cc = r["conditions"]["c_further_reducible"]
    check("(6) 係数が既に合っていて残りが形の残差だけ → (c) さらに減らせない → escalate", not r["apply"] and not cc["ok"]
          and r["conditions"]["b_prev_update_reduced_error"]["ok"], f"減る割合 {cc['reduction']:.3f} < {cc['min_reduction']}")
    r = WQ.recalibration_rule(cur, v5b, v5c, fP, prev_run=prev, recalibrations_done=1)
    check("(6) 較正し直しが済んでいる (1 回まで) → (d) で escalate", not r["apply"] and not r["conditions"]["d_count"]["ok"])
    r = WQ.recalibration_rule(cur, v5b, v5c, fU, prev_run=prev)
    check("(6) 計算の成立が UNVERIFIED → (e) で escalate", not r["apply"] and not r["conditions"]["e_flow_verified"]["ok"])
    r = WQ.recalibration_rule(cur, v5b, {"verdict": "INDETERMINATE"}, fP, prev_run=prev)
    check("(6) V5c が判定不能 → (a) で escalate", not r["apply"] and not r["conditions"]["a_both_judged"]["ok"])
    r = WQ.recalibration_rule(cur, v5b, v5c, fP, prev_run=prev, min_reduction=0.999999)
    check("(6) 減る割合の下限を上げると (c) で escalate", not r["apply"] and not r["conditions"]["c_further_reducible"]["ok"])

# ------------------------------------------------------------------ (7) 実データ
print("== (7) case/71 run_0010・run_0011")
C71 = ROOT / "case/71.m55_combustion_wt"
R08 = C71 / "_runs/run_0008_ns_b0"
CARD = C71 / "card_m55_3coef.yaml"
# 2026-10-15 の判定 (V5b ±0.05 %・V5c 片振幅 0.05 % で FAIL → escalate) の再現なので、旧の公差を明示する。2026-10-16 のユーザ決定で
# 既定は 0.1 % になった (plan tooling-nozzle-card-local-mesh §6 V5b)。新しい既定での判定は下の「(7′)」で見る
OLD_TOLS = {"exit_core_M_tol_rel": 5e-4, "line_tol_halfamp_pct": 0.05}
cases = [("run_0010", C71 / "run_0010_ns_b0", {"delta_next": C71 / "_runs/run_0010_ns_b0/judge/extract/delta_r_next.csv.gz",
                                                "cli_tols": {"line_x_range_m": [4.0, 7.95], **OLD_TOLS}}, 5.502841, "UNVERIFIED"),
         ("run_0011", Path(os.environ.get("FORGE_CASE71_RUN_0011", C71 / "run_0011_ns_v2")),
          {"problem": C71 / "problem_m55_p5_settle_v2.yaml", "cli_tols": dict(OLD_TOLS)}, 5.502842, "PASS")]
for nm, run, kw, M_want, flow_want in cases:
    need = [run / "prepare_info.json", run / "nozzle.h5", run / "bcondConfig.yaml", R08 / "prepare_info.json"]
    if not all(p.is_file() for p in need) or not any(run.glob("res_*.h5")):
        SKIPPED.append(f"(7) {nm} ({run} に run が無い)")
        continue
    with tempfile.TemporaryDirectory() as td:
        rec = WQ.qualify(run, CARD, kw.get("problem"), prev_run=R08, prev_delta_next=R08 / "judge/extract/delta_r_next.csv.gz",
                         delta_next=kw.get("delta_next"), cli_tols=kw.get("cli_tols"), wall_txt=Path(td) / "wr.txt")
    c = rec["checks"]
    check(f"(7) {nm}: V5b FAIL (出口コア M {M_want})", c["V5b"]["verdict"] == "FAIL" and round(c["V5b"]["M"], 6) == M_want, f"{c['V5b']['M']!r}")
    check(f"(7) {nm}: V5c FAIL (η 0.03・x 4.0〜7.95 m: 片振幅 0.0573 %、最大偏差 0.0623 % は許容内)",
          c["V5c"]["verdict"] == "FAIL" and round(c["V5c"]["halfamp_pct"], 4) == 0.0573 and round(c["V5c"]["maxabs_pct"], 4) == 0.0623
          and c["V5c"]["checks"]["maxabs"] and not c["V5c"]["checks"]["halfamp"] and c["V5c"]["x_range_m"] == [4.0, 7.95] and c["V5c"]["eta"] == 0.03,
          f"{c['V5c'].get('maxabs_pct')} / {c['V5c'].get('halfamp_pct')}")
    check(f"(7) {nm}: η 0.1 は比較用に併記 (片振幅 0.076 %)", round(c["V5c"]["compare"]["eta0.1"]["halfamp_pct"], 3) == 0.076)
    check(f"(7) {nm}: 壁解像 PASS・計算の成立 {flow_want}", c["wall_resolution"]["verdict"] == "PASS" and c["flow"]["status"] == flow_want,
          f"{c['wall_resolution'].get('verdict_line')} / {c['flow']['status']}")
    rr = rec["recalibration_rule"]
    cnd = rr["conditions"]
    check(f"(7) {nm}: back_to_correction・escalate (較正し直しの規則に当たらない; plan §4.3 の 2026-10-15 の判断)",
          rec["result"] == "back_to_correction" and rec["next"] == "escalate" and not rr["apply"], f"{rec['result']} / {rec['next']}")
    check(f"(7) {nm}: (a) 両方を判定・(b) 前の更新で誤差が縮んだ (2.14 → 0.62 mm) は満たし、(c) さらに減らせない (0.62 → 0.61 mm)",
          cnd["a_both_judged"]["ok"] and cnd["b_prev_update_reduced_error"]["ok"] and not cnd["c_further_reducible"]["ok"]
          and round(cnd["b_prev_update_reduced_error"]["prev_max_abs_mm"], 2) == 2.14
          and round(cnd["c_further_reducible"]["current_max_abs_mm"], 2) == 0.62 and round(cnd["c_further_reducible"]["projected_max_abs_mm"], 2) == 0.61,
          f"{cnd['c_further_reducible'].get('current_max_abs_mm'):.4f} → {cnd['c_further_reducible'].get('projected_max_abs_mm'):.4f} mm")
    if nm == "run_0011":
        check("(7) run_0011: 窓は問題 YAML の evaluate.settle (監視の line_mach と同じ線)、場は停止規則 v2 の stop_step 46202",
              "evaluate.settle" in rec["tolerances"]["line_x_range_m"]["source"] and rec["field"]["step"] == 46202
              and c["V5c"].get("matches_settle_line") is True)
        err = io.StringIO()
        with redirect_stdout(io.StringIO()), redirect_stderr(err):
            rc = WQ.main([str(run), "--card", str(CARD), "--problem", str(C71 / "problem_m55_p5_settle_v2.yaml"), "--prev-run", str(R08),
                          "--prev-delta-next", str(R08 / "judge/extract/delta_r_next.csv.gz"), "--no-write",
                          "--exit-tol-rel", "5e-4", "--line-tol-halfamp", "0.05"])
        check("(7) run_0011: CLI の終了コード 3 (escalate、旧の公差)", rc == 3, err.getvalue()[-200:])
    else:
        check("(7) run_0010: 寸法の関門・成果物 PASS (手元の deliverables)", c["deliverables"]["status"] == "PASS"
              and c["deliverables"]["dimensions"]["status"] == "PASS", f"{c['deliverables'].get('problems')}")
    # (7′) 新しい既定 (V5b ±0.1 %・V5c 片振幅 0.1 %; 2026-10-16 ユーザ決定) では同じ run の V5b・V5c が PASS
    tols_new = {k: v for k, v in (kw.get("cli_tols") or {}).items() if k not in OLD_TOLS}
    with tempfile.TemporaryDirectory() as td:
        rec2 = WQ.qualify(run, CARD, kw.get("problem"), prev_run=R08, prev_delta_next=R08 / "judge/extract/delta_r_next.csv.gz",
                          delta_next=kw.get("delta_next"), cli_tols=tols_new or None, wall_txt=Path(td) / "wr.txt")
    c2 = rec2["checks"]
    check(f"(7′) {nm}: 新しい既定の公差 (±0.1 %・片振幅 0.1 %) では V5b・V5c とも PASS (値は旧の判定と同じ)",
          c2["V5b"]["verdict"] == "PASS" and c2["V5c"]["verdict"] == "PASS" and c2["V5b"]["M"] == c["V5b"]["M"]
          and c2["V5c"]["halfamp_pct"] == c["V5c"]["halfamp_pct"] and rec2["tolerances"]["exit_core_M_tol_rel"]["value"] == 1e-3
          and rec2["tolerances"]["line_tol_halfamp_pct"]["value"] == 0.1, f"判定 {rec2['result']} / {rec2['next']}")

for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if not FAIL else f'{len(FAIL)} FAILURES: {FAIL}'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
