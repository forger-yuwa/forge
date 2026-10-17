"""風洞ノズルの設計の判定 `feedback/wt_qualify.py` の試験 (plan tooling-nozzle-card-local-mesh §5.1 #6b・§4.3・§6 V5)。

合成のデータ:
 (1) 公差・窓の出所の順 (CLI → カードの qualification → 問題 YAML の qualification → 問題 YAML の evaluate.settle → run の
     settle_criteria.yaml → 既定 (preset wt_nozzle_v2 の線・plan §6 V5b・V5c))。知らないキーは ValueError。
 (2) V5b・V5c を合成の構造格子の場で: 線形の M で最大偏差・片振幅・出口コア M が式どおり、合否、軸の節点の重み > 1e−6・窓が格子の外は
     判定不能、比較の線 (η 0.1・軸) を併記。
 (3) 計算の成立 (停止規則の状態): v2 は停止の時の判定と結び付いた (stop_binding を差し替え) ときだけ REACHED + SETTLED → PASS、
     CENSORED・再評価が SETTLED でない → FAIL。結び付かない・状態のファイルだけ (本物の照合)・v2 の基準があるのに状態が無い → INDETERMINATE
     (§5.1 #7a の続き、codex result 2 回目 M1)。状態が無い旧い run → UNVERIFIED、v1 の REACHED → PASS。bound_lines・monitored_line。
 (4) 成果物と寸法の関門: そろって PASS、ファイルの改変・部品の欠け・寸法の FAIL・別のカード → FAIL、無ければ ABSENT、無くてよい部品の指定。
 (4p) 成果物の生成元の照合 (plan §5.1 #6b の続き、2026-10-17 codex result 3 回目 Major): MANIFEST.json の provenance と今の評価の対象を
     内容で照合。記録が無い (古い成果物) → ABSENT、記録の sha256 を 1 つ変えた写し (壁・格子・設定 2 つ・場・壁ダンプ・カード・問題 YAML)・
     書き出しの後に run の壁・全場・設定が変わった・評価の場の step が違う・停止の step が違う → MISMATCH、評価の場・カード・問題 YAML を
     渡さない・生成元の場の記録が無い・壁ダンプの記録が空・書式が違う → UNCHECKABLE。どれも成果物と寸法の関門 INDETERMINATE
     (寸法 FAIL の成果物でも寸法の補正にしない → incomplete・escalate)。別の置き場所への写し (名前も違う・中身は同じ) は MATCH・PASS。
 (5) 結果と次の手 (`decide`): accepted・成果物待ち (make_deliverables)・計算の成立を確かめていない・判定不能・反復の上限 (rejected)・
     壁の形だけの補正 (規則で recalibrate / escalate)・壁解像や寸法を含む補正 (escalate、規則は評価しない)。
 (6) 較正し直しの規則 (合成の run): (a)〜(e) をそろえて recalibrate、1 つずつ欠くと escalate (前の run が無い・前の更新で誤差が増えた・
     さらに減らせない・回数・計算の成立)。
 (8) 線の定常の保証 (plan §5.1 #7a、codex result 2026-10-17 M1): 監視の線との違い (η・x 窓・点数)・line_mach が無い・停止の step と違う・
     停止の時の判定と結び付かない・監視の線が停止の時に SETTLED でない (2 回目 M1) → 保証の欠け。欠けたときは本段の全場の系列 (load_field を差し替えた合成の場) を check_quasisteady で判定し、STEADY なら合否に進み、
     動いている・4 枚未満・本段の開始の記録が無い・評価の場で終わらない → 判定不能。壁解像の照合 (§5.1 #7b、M2): 偽のツールの出力で、
     EVALUATED 行が無い・要求した壁の欠け・要求した step でない・bcond と違う physID・要求していない壁・停止の step でない → 判定不能。
実データ (手元にあるとき):
 (7) case/71 run_0010・run_0011: V5b FAIL (出口コア M 5.502841 / 5.502842)・線 η 0.03 の片振幅 0.0573 %、最大偏差 0.0623 % (許容内) で
     線だけの合否 FAIL・壁解像 PASS、較正し直しの規則に当たらない (形の誤差の射影 0.624 → 0.612 mm、1.9 % しか減らない)
     (plan §4.3 の 2026-10-15 の判断どおり)。run_0011 は監視の線・停止の step なので V5c FAIL → back_to_correction・escalate。run_0010 は
     停止規則 v2 より前の run (線を監視していない・全場 1 枚) なので §5.1 #7a から V5c 判定不能 → incomplete・escalate (規則は直接評価して見る)。
     run_0011 は環境変数 FORGE_CASE71_RUN_0011 か case/71/run_0011_ns_v2 にあるとき。
 (9) case/71 run_0014_chain_v2 の写し (h5 は symlink、元の run には書かない; FORGE_CASE71_RUN_0014 で場所を変えられる): 正しい呼び方は
     今どおり accepted・壁解像 PASS (wall/3@46171)。line_eta 0.06 は incomplete (修正前は accepted)。check_wall_resolution の
     `--step 46170 --groups wall,missing_wall` 等は判定不能・終了コード 2 (修正前は step 46171 の wall だけで PASS・0)。
     2 回目 M1: メモリ上で基準の読み取り値だけ η 0.06・基準のファイルを書き換えた写し (η・流量の閾値)・停止の時の判定の記録の無い写し
     → incomplete (修正前は accepted)。2 回目 M2: bcond に一度もダンプの無い壁を足した写しで --groups 無し → 判定不能・2 (修正前は PASS・0)。
     試走 A (case/90 run_0003_ns_p0・run_0005_ns_p1、読むだけ) の計算の成立は停止の時の判定と結び付いて PASS。
     3 回目 Major (2026-10-17): run_0014 を自分の成果物 (生成元の記録つき; 2026-10-17 に作り直した `deliverables/`) で評価 → accepted、
     その写しを別の置き場所に置いても accepted。run_0010 の成果物 (古い成果物) を渡す → incomplete・終了コード 4 (修正前は accepted・0)、
     run_0010 の成果物に新しい道具の生成元の記録を足した写し → MISMATCH (壁・格子・設定・場の step・停止の step)、作り直す前の
     run_0014 の成果物 (`_deliverables_pre_provenance/`) → ABSENT、記録の格子の sha256 を 1 つ変えた写し → MISMATCH (どれも incomplete)。
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

    def flow(name, state=None, v2=False, bind=None):
        """bind: None = 本物の照合 (deltastar_loop._v2_binding)、dict = stop_binding の戻りに差し替える。"""
        d = td / name
        d.mkdir()
        if state is not None:
            (d / "STOP_RULE_STATE.json").write_text(json.dumps(state))
        if v2:
            (d / "settle_criteria.yaml").write_text("quantities: []\n")
        orig = WQ.stop_binding
        if bind is not None:
            WQ.stop_binding = lambda run_dir: bind
        try:
            return WQ.flow_check(d)
        finally:
            WQ.stop_binding = orig

    # 停止の時の判定に結び付いた (stop_binding の pass) ときだけ状態と再評価で PASS・FAIL を決める (plan §5.1 #7a の続き、
    # 2026-10-17 codex result 2 回目 M1)。以前は状態のファイルの REACHED・SETTLED だけで PASS だった (下の「照合の記録の無い状態」)
    line_c = {"kind": "line_mach", "name": "line_mach", "eta": 0.03, "x_range_m": [4.0, 7.95], "n": 2001, "Md": MD}
    vb = {"verdict": "SETTLED", "criteria": {"content": {"quantities": [line_c]}}, "quantities": [{"name": "line_mach", "verdict": "SETTLED"}]}
    B = {"pass": True, "problems": [], "verdict": vb, "verdict_file": "SETTLE_VERDICT.json", "stop_step": 46000, "criteria_sha256": "x"}
    r = flow("a", {"method": "v2", "status": "REACHED", "stop_step": 46000, "reeval": {"verdict": "SETTLED"}}, True, bind=B)
    check("(3) v2 の REACHED + 結び付いた再評価 SETTLED → PASS、監視の線は停止の時の判定の基準の中身から (SETTLED)",
          r["status"] == "PASS" and r["binding"]["pass"] and r["monitored_lines"] == [{**{k: line_c[k] for k in ("name", "eta", "x_range_m", "n", "Md")},
                                                                                   "verdict": "SETTLED"}], f"{r}")
    check("(3) v2 の CENSORED (結び付いた) → FAIL",
          flow("b", {"method": "v2", "status": "CENSORED", "reeval": {"verdict": "SETTLED"}}, True, bind=B)["status"] == "FAIL")
    check("(3) v2 の結び付いた再評価が SETTLED でない → FAIL",
          flow("c", {"method": "v2", "status": "REACHED", "reeval": {"verdict": "SETTLED"}}, True,
               bind={**B, "verdict": {**vb, "verdict": "NOT_SETTLED"}})["status"] == "FAIL")
    r = flow("c2", {"method": "v2", "status": "REACHED", "reeval": {"verdict": "SETTLED"}}, True,
             bind={**B, "pass": False, "problems": ["基準の sha256 が違う"], "verdict": None})
    check("(3) v2 の REACHED・SETTLED でも停止の時の判定と結び付かない → INDETERMINATE (監視の線なし)",
          r["status"] == "INDETERMINATE" and r["monitored_lines"] == [] and "基準の sha256" in r["why"], r.get("why"))
    r = flow("c3", {"method": "v2", "status": "REACHED", "stop_step": 46000, "reeval": {"verdict": "SETTLED"}}, True)
    check("(3) 状態のファイルだけ (停止の時の判定の記録・起動の記録が無い) → 本物の照合で INDETERMINATE (修正前は PASS)",
          r["status"] == "INDETERMINATE" and not r["binding"]["pass"] and r["binding"]["n_problems"] > 0, r.get("why", "")[:200])
    # 2026-10-17 から、記録 (STOP_RULE_STATE.json) が無い v2 の run は FAIL でなく INDETERMINATE (確かめられない; 結果は同じ incomplete)
    check("(3) v2 の基準があるのに状態が無い → INDETERMINATE (前は FAIL; 結果はどちらも incomplete)", flow("d", None, True)["status"] == "INDETERMINATE")
    check("(3) 状態が無い旧い run → UNVERIFIED", flow("e")["status"] == "UNVERIFIED")
    check("(3) v1 の REACHED → PASS・CENSORED → FAIL", flow("f", {"version": 1, "status": "REACHED"})["status"] == "PASS"
          and flow("g", {"version": 1, "status": "CENSORED"})["status"] == "FAIL")
    lines = WQ.bound_lines({"criteria": {"content": {"quantities": [line_c, {**line_c, "name": "line2", "eta": 0.06}, {"kind": "exit_core_mach", "name": "e"}]}},
                            "quantities": [{"name": "line_mach", "verdict": "SETTLED"}]})
    check("(3) bound_lines: line_mach を全部、量の verdict つき (判定に無い量は None)",
          [ln["name"] for ln in lines] == ["line_mach", "line2"] and lines[0]["verdict"] == "SETTLED" and lines[1]["verdict"] is None)
    check("(3) monitored_line: 評価の線と同じ線を選ぶ (無ければ最初の線、線が無ければ None)",
          WQ.monitored_line({"monitored_lines": lines}, {"eta": 0.06, "x_range_m": [4.0, 7.95], "n": 2001})["name"] == "line2"
          and WQ.monitored_line({"monitored_lines": lines}, {"eta": 0.05, "x_range_m": [4.0, 7.95], "n": 2001})["name"] == "line_mach"
          and WQ.monitored_line({}, {"eta": 0.03}) is None)

# ------------------------------------------------------------------ (4) 成果物と寸法の関門
print("== (4) 成果物と寸法の関門")
with tempfile.TemporaryDirectory() as td:
    import shutil
    td = Path(td)
    run = td / "run_x"
    run.mkdir()
    card = td / "card.yaml"
    card.write_text("M_exit: 5.5\n")
    prob = td / "problem.yaml"
    prob.write_text("spec: {r_throat: 0.1}\n")
    check("(4) 成果物が無い → ABSENT", WQ.deliverables_check(run, None, card)["status"] == "ABSENT")
    from forge_design.export import deliverables as DV
    from forge_design.export.deliverables import PARTS
    # 生成元の照合 (2026-10-17 codex result 3 回目 Major) に使う run の中身 (壁・格子・設定・全場・壁ダンプ・停止規則の状態; 中身は sha256 だけ見る)
    (run / "wall_repr.json").write_text('{"physical_wall_repr": "single_bspline"}')
    (run / "solverConfig.yaml").write_text('mesh: {meshFileName: "nozzle.h5"}\n')
    (run / "bcondConfig.yaml").write_text("wall: {physID: 3, kind: wall_isothermal}\n")
    for nm_ in ("nozzle.h5", "res_100.h5", "res_wall_3_100.h5", "res_200.h5", "res_wall_3_200.h5"):
        (run / nm_).write_bytes(nm_.encode() * 10)
    (run / "STOP_RULE_STATE.json").write_text(json.dumps({"method": "v2", "status": "REACHED", "stop_step": 100}))
    FLD = {"step": 100, "res": "res_100.h5"}
    od = run / "deliverables"
    (od / "cad").mkdir(parents=True)
    (od / "cad" / "wall.step").write_text("STEP")
    dim = {"format": "forge.dimension_gate/v1", "verdict": "PASS", "items": {"L_div": {"pass": True}}, "card": "card.yaml", "card_sha256": sha(card)}
    (od / "DIMENSIONS.json").write_text(json.dumps(dim))
    man = {"run": str(run), "verdict": "PASS", "failed_parts": [], "parts": {p: {"ok": True} for p in PARTS}, "skipped": {}, "errors": {},
           "files_sha256": {"DIMENSIONS.json": sha(od / "DIMENSIONS.json"), "cad/wall.step": sha(od / "cad" / "wall.step")},
           "provenance": DV.provenance_record(run, 100, ["res_wall_3_100.h5"], card, prob)}
    (od / "MANIFEST.json").write_text(json.dumps(man))

    def dchk(card_=card, field=FLD, prob_=prob, **kw):
        return WQ.deliverables_check(run, None, card_, field=field, problem_path=prob_, **kw)

    r = dchk()
    check("(4) そろって PASS (寸法の関門 PASS・カードの sha256 一致・生成元 MATCH)", r["status"] == "PASS" and r["dimensions"]["status"] == "PASS"
          and r["provenance"]["status"] == "MATCH", f"{r.get('problems')}")
    other = td / "other.yaml"
    other.write_text("M_exit: 6.0\n")
    # 寸法の関門のカードの照合 (DIMENSIONS.json の card_sha256) は生成元の照合と別に残る: 生成元の記録のカードが other の成果物
    # (書き出しには other を渡したのに DIMENSIONS.json は card で判定された) を other で評価すると、生成元は一致・寸法は FAIL
    man_o = dict(man, provenance=DV.provenance_record(run, 100, ["res_wall_3_100.h5"], other, prob))
    (od / "MANIFEST.json").write_text(json.dumps(man_o))
    r = dchk(card_=other)
    check("(4) 寸法の関門が別のカードで判定された (生成元は一致) → 寸法 FAIL", r["dimensions"]["status"] == "FAIL"
          and r["provenance"]["status"] == "MATCH", f"{r['dimensions']['problems']}")
    (od / "MANIFEST.json").write_text(json.dumps(man))
    # 2026-10-17 から: 生成元 (card) と違うカード (other) で評価すると生成元の照合で止まり、寸法は INDETERMINATE (前は寸法 FAIL →
    # 寸法の解き直しの補正に回った; 別のカードで作った成果物の寸法はこの評価の判定にしない)
    r = dchk(card_=other)
    check("(4) 生成元のカードと違うカードで評価 → 生成元 MISMATCH・成果物と寸法の関門 INDETERMINATE (寸法の FAIL にしない)",
          r["status"] == "INDETERMINATE" and r["dimensions"]["status"] == "INDETERMINATE" and r["provenance"]["status"] == "MISMATCH"
          and r["dimensions"]["status_before_provenance"] == "FAIL" and any("カード" in p for p in r["provenance"]["problems"]), r.get("why"))
    (od / "cad" / "wall.step").write_text("STEP changed")
    r = dchk()
    check("(4) 列挙したファイルが改変された → FAIL", r["status"] == "FAIL" and any("sha256" in p for p in r["problems"]), f"{r['problems']}")
    (od / "cad" / "wall.step").write_text("STEP")
    man["parts"].pop("fluent")
    man["skipped"] = {"fluent": "--skip で指定"}
    (od / "MANIFEST.json").write_text(json.dumps(man))
    r = dchk()
    check("(4) 必要な部品が欠けた → FAIL", r["status"] == "FAIL" and any("fluent" in p for p in r["problems"]))
    r = dchk(optional_parts=("fluent",))
    check("(4) 無くてよい部品に指定すれば PASS", r["status"] == "PASS", f"{r.get('problems')}")
    dim["verdict"] = "FAIL"
    (od / "DIMENSIONS.json").write_text(json.dumps(dim))
    man["files_sha256"]["DIMENSIONS.json"] = sha(od / "DIMENSIONS.json")
    (od / "MANIFEST.json").write_text(json.dumps(man))
    r = dchk(optional_parts=("fluent",))
    check("(4) 寸法の関門の FAIL → 寸法 FAIL", r["dimensions"]["status"] == "FAIL")
    dim["verdict"] = "PASS"
    (od / "DIMENSIONS.json").write_text(json.dumps(dim))
    man["files_sha256"]["DIMENSIONS.json"] = sha(od / "DIMENSIONS.json")
    man["parts"]["fluent"] = {"ok": True}
    man["skipped"] = {}
    (od / "MANIFEST.json").write_text(json.dumps(man))

    # ---- (4p) 生成元の照合 (plan §5.1 #6b の続き、2026-10-17 codex result 3 回目 Major) ----
    import copy as _copy

    def with_man(m):
        (od / "MANIFEST.json").write_text(json.dumps(m))

    def bad(r, status):
        return r["status"] == "INDETERMINATE" and r["dimensions"]["status"] == "INDETERMINATE" and r["provenance"]["status"] == status

    with_man({k: v for k, v in man.items() if k != "provenance"})
    r = dchk()
    check("(4p) 生成元の記録の無い古い成果物 → 生成元 ABSENT・成果物と寸法 INDETERMINATE・理由に作り直し", bad(r, "ABSENT")
          and "作り直す" in r["why"] and r["status_before_provenance"] == "PASS", r.get("why"))
    for nm_, path_ in (("壁 (wall_repr.json)", ("wall", "sha256")), ("格子 (nozzle.h5)", ("mesh", "sha256")),
                       ("設定 (solverConfig.yaml)", ("config", "solverConfig.yaml")), ("設定 (bcondConfig.yaml)", ("config", "bcondConfig.yaml")),
                       ("場 (res_100.h5)", ("field", "sha256")), ("壁ダンプ", ("field", "wall_dumps", "res_wall_3_100.h5")),
                       ("カード", ("card", "sha256")), ("問題 YAML", ("problem", "sha256"))):
        m2 = _copy.deepcopy(man)
        d_ = m2["provenance"]
        for k in path_[:-1]:
            d_ = d_[k]
        d_[path_[-1]] = "0" * 64
        with_man(m2)
        r = dchk()
        check(f"(4p) 記録の {nm_} の sha256 を 1 つ変えた写し → 生成元 MISMATCH・INDETERMINATE", bad(r, "MISMATCH")
              and len(r["provenance"]["problems"]) == 1 and nm_.split(" ")[0] in r["provenance"]["problems"][0], r["provenance"]["problems"])
    with_man(man)
    for nm_, fn in (("壁 (wall_repr.json)", "wall_repr.json"), ("全場 res_100.h5", "res_100.h5"), ("設定 bcondConfig.yaml", "bcondConfig.yaml")):
        keep = (run / fn).read_bytes()
        (run / fn).write_bytes(keep + b"\n")
        r = dchk()
        (run / fn).write_bytes(keep)
        check(f"(4p) 書き出しの後に run の {nm_} が変わった → MISMATCH", bad(r, "MISMATCH"), r["provenance"]["problems"])
    r = dchk(field={"step": 200, "res": "res_200.h5"})
    check("(4p) 評価の場の step (200) が生成元 (100) と違う → MISMATCH", bad(r, "MISMATCH")
          and any("場の step" in p for p in r["provenance"]["problems"]), r["provenance"]["problems"])
    keep = (run / "STOP_RULE_STATE.json").read_text()
    (run / "STOP_RULE_STATE.json").write_text(json.dumps({"method": "v2", "status": "REACHED", "stop_step": 200}))
    r = dchk()
    (run / "STOP_RULE_STATE.json").write_text(keep)
    check("(4p) 今の run の停止の step (200) が生成元の記録 (100) と違う → MISMATCH", bad(r, "MISMATCH")
          and any("停止の step" in p for p in r["provenance"]["problems"]), r["provenance"]["problems"])
    for nm_, kw in (("評価の場を渡さない", {"field": None}), ("評価にカードを渡さない (生成元はカードあり)", {"card_": None}),
                    ("評価に問題 YAML を渡さない (生成元はあり)", {"prob_": None})):
        r = dchk(**kw)
        check(f"(4p) {nm_} → 照合できない (UNCHECKABLE)・INDETERMINATE", bad(r, "UNCHECKABLE"), r["provenance"]["problems"])
    for nm_, edit in (("生成元の場の記録が無い (CFD 壁面量表を作っていない)", lambda m: m["provenance"].update(field=None, field_why="試験")),
                      ("生成元の壁ダンプの記録が空", lambda m: m["provenance"]["field"].update(wall_dumps={})),
                      ("生成元の記録の書式が違う", lambda m: m["provenance"].update(format="forge.deliverables.provenance/v0"))):
        m2 = _copy.deepcopy(man)
        edit(m2)
        with_man(m2)
        r = dchk()
        check(f"(4p) {nm_} → UNCHECKABLE・INDETERMINATE", bad(r, "UNCHECKABLE"), r["provenance"]["problems"])
    m2 = _copy.deepcopy(man)
    m2["provenance"].update(card=None, problem=None)
    with_man(m2)
    r = dchk(card_=None, prob_=None)
    check("(4p) カード・問題 YAML を書き出しにも評価にも渡していない → その 2 つは比べるものが無い (ほかは一致で MATCH)",
          r["provenance"]["status"] == "MATCH" and r["provenance"]["compared"]["カード"]["recorded"] is None, r["provenance"]["problems"])
    with_man(man)
    # 別の置き場所への写し (run の名前も置き場所も違う・中身は同じ) は内容で一致 → PASS (名前の違いは警告だけ)
    run_c = td / "elsewhere" / "run_copy"
    shutil.copytree(run, run_c)
    od_c = td / "deliv_copy"
    shutil.copytree(od, od_c)
    r = WQ.deliverables_check(run_c, od_c, card, field=FLD, problem_path=prob)
    check("(4p) run と成果物を別の置き場所に写した (名前も違う・中身は同じ) → 生成元 MATCH・PASS、名前の違いは警告だけ",
          r["status"] == "PASS" and r["provenance"]["status"] == "MATCH" and any("名前が違う" in w for w in r["warnings"]), r.get("problems"))
    # 寸法の関門が FAIL の成果物でも、生成元が一致しなければ寸法の補正 (back_to_correction) に回さない
    dim["verdict"] = "FAIL"
    (od / "DIMENSIONS.json").write_text(json.dumps(dim))
    m2 = _copy.deepcopy(man)
    m2["files_sha256"]["DIMENSIONS.json"] = sha(od / "DIMENSIONS.json")
    m2["provenance"]["mesh"]["sha256"] = "0" * 64
    with_man(m2)
    r = dchk()
    d_ = WQ.decide({"status": "PASS"}, {"verdict": "PASS"}, {"verdict": "PASS"}, {"verdict": "PASS"}, r)
    check("(4p) 寸法 FAIL の成果物でも生成元が一致しなければ寸法 INDETERMINATE → decide は incomplete・escalate (寸法の補正にしない)",
          bad(r, "MISMATCH") and r["dimensions"]["status_before_provenance"] == "FAIL" and (d_["result"], d_["next"]) == ("incomplete", "escalate")
          and not d_["corrections"] and any("作り直す" in m_ for m_ in d_["missing"]), f"{d_['result']} {d_['next']} {d_['corrections']}")

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
check("(5) 計算の成立 INDETERMINATE (停止の時の判定と結び付かない) → incomplete / escalate (すべて合格でも・公差が FAIL でも)",
      dec({"status": "INDETERMINATE", "why": "b"}, P, P, P, dP) == ("incomplete", "escalate")
      and dec({"status": "INDETERMINATE", "why": "b"}, Fa, P, P, dA, rule_fn=rule(True)) == ("incomplete", "escalate"))
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
dI = {"status": "INDETERMINATE", "why": "成果物の生成元が今の評価の対象と一致しない — 作り直す", "dimensions": {"status": "INDETERMINATE"}}
check("(5) 成果物の生成元が一致しない (成果物・寸法 INDETERMINATE; 2026-10-17 codex result 3 回目 Major) → incomplete / escalate "
      "(ABSENT 以外は escalate の既存の規則)", dec(fP, P, P, P, dI, rule_fn=rule(True)) == ("incomplete", "escalate"))
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
    # run_0010 は停止規則 v2 より前の run (settle_criteria.yaml・STOP_RULE_STATE.json が無い) で、線 η 0.03 を監視しておらず全場も
    # res_86000 の 1 枚だけ。§5.1 #7a (2026-10-17 codex result M1) から、監視されていない線は全場の系列の準定常 (check_quasisteady) が
    # STEADY のときだけ合否に進むので、線の値と線だけの合否 (FAIL) は今までどおりで V5c は判定不能になる (期待値の変更の理由)
    v5c_want = "FAIL" if nm == "run_0011" else "INDETERMINATE"
    lv = c["V5c"].get("line_verdict", c["V5c"]["verdict"])
    check(f"(7) {nm}: 線 η 0.03・x 4.0〜7.95 m の値 (片振幅 0.0573 %、最大偏差 0.0623 % は許容内) で線だけの合否 FAIL、V5c {v5c_want}",
          c["V5c"]["verdict"] == v5c_want and lv == "FAIL" and round(c["V5c"]["halfamp_pct"], 4) == 0.0573 and round(c["V5c"]["maxabs_pct"], 4) == 0.0623
          and c["V5c"]["checks"]["maxabs"] and not c["V5c"]["checks"]["halfamp"] and c["V5c"]["x_range_m"] == [4.0, 7.95] and c["V5c"]["eta"] == 0.03,
          f"{c['V5c'].get('maxabs_pct')} / {c['V5c'].get('halfamp_pct')} / {c['V5c']['verdict']}")
    if nm == "run_0010":
        ser = c["V5c"]["steadiness"]["series"]
        # 系列は作れない (この run の MAIN_STAGE_STARTED は start_step を持たない旧い形で区間を特定できない; 全場も res_86000 の 1 枚)
        check("(7) run_0010: 判定不能の理由 = 監視の line_mach が無い・停止の step が無い・全場の系列を作れない → incomplete・escalate",
              any("line_mach" in g for g in c["V5c"]["steadiness"]["monitor_gaps"]) and any("停止の step" in g for g in c["V5c"]["steadiness"]["monitor_gaps"])
              and ser["status"] == "INDETERMINATE" and rec["result"] == "incomplete" and rec["next"] == "escalate"
              and rec["recalibration_rule"] is None, f"{c['V5c'].get('why')}")
    check(f"(7) {nm}: η 0.1 は比較用に併記 (片振幅 0.076 %)", round(c["V5c"]["compare"]["eta0.1"]["halfamp_pct"], 3) == 0.076)
    check(f"(7) {nm}: 壁解像 PASS・計算の成立 {flow_want}", c["wall_resolution"]["verdict"] == "PASS" and c["flow"]["status"] == flow_want,
          f"{c['wall_resolution'].get('verdict_line')} / {c['flow']['status']}")
    if nm == "run_0011":
        rr = rec["recalibration_rule"]
        check(f"(7) {nm}: back_to_correction・escalate (較正し直しの規則に当たらない; plan §4.3 の 2026-10-15 の判断)",
              rec["result"] == "back_to_correction" and rec["next"] == "escalate" and not rr["apply"], f"{rec['result']} / {rec['next']}")
    else:
        # run_0010 は判定不能で規則まで進まないので、線だけの合否で規則を直接評価して (a)〜(c) の値を見る (2026-10-15 の判断の再現)
        rr = WQ.recalibration_rule(run, c["V5b"], {"verdict": lv}, c["flow"], delta_next=kw.get("delta_next"), prev_run=R08,
                                   prev_delta_next=R08 / "judge/extract/delta_r_next.csv.gz")
    cnd = rr["conditions"]
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
        # 2026-10-17 から (codex result 3 回目 Major): 手元の run_0010 の成果物は生成元の記録 (MANIFEST.json の provenance) を書く前の
        # 古い成果物なので、生成元 ABSENT → 成果物・寸法の関門は INDETERMINATE (前は PASS)。成果物の中身の照合と寸法の関門は今どおり PASS
        dv_ = c["deliverables"]
        check("(7) run_0010: 手元の deliverables は生成元の記録の無い古い成果物 → 成果物・寸法の関門 INDETERMINATE (中身の照合・寸法は PASS)",
              dv_["status"] == "INDETERMINATE" and dv_["provenance"]["status"] == "ABSENT" and dv_["status_before_provenance"] == "PASS"
              and dv_["dimensions"]["status"] == "INDETERMINATE" and dv_["dimensions"]["status_before_provenance"] == "PASS", f"{dv_.get('problems')}")
    # (7′) 新しい既定 (V5b ±0.1 %・V5c 片振幅 0.1 %; 2026-10-16 ユーザ決定) では同じ run の V5b・V5c が PASS
    tols_new = {k: v for k, v in (kw.get("cli_tols") or {}).items() if k not in OLD_TOLS}
    with tempfile.TemporaryDirectory() as td:
        rec2 = WQ.qualify(run, CARD, kw.get("problem"), prev_run=R08, prev_delta_next=R08 / "judge/extract/delta_r_next.csv.gz",
                          delta_next=kw.get("delta_next"), cli_tols=tols_new or None, wall_txt=Path(td) / "wr.txt")
    c2 = rec2["checks"]
    # run_0010 の V5c は上と同じ理由 (監視されていない線・全場 1 枚) で判定不能、線だけの合否が PASS (§5.1 #7a)
    check(f"(7′) {nm}: 新しい既定の公差 (±0.1 %・片振幅 0.1 %) では V5b・V5c (run_0010 は線だけの合否) とも PASS (値は旧の判定と同じ)",
          c2["V5b"]["verdict"] == "PASS" and c2["V5c"].get("line_verdict", c2["V5c"]["verdict"]) == "PASS"
          and c2["V5c"]["verdict"] == ("PASS" if nm == "run_0011" else "INDETERMINATE") and c2["V5b"]["M"] == c["V5b"]["M"]
          and c2["V5c"]["halfamp_pct"] == c["V5c"]["halfamp_pct"] and rec2["tolerances"]["exit_core_M_tol_rel"]["value"] == 1e-3
          and rec2["tolerances"]["line_tol_halfamp_pct"]["value"] == 0.1, f"判定 {rec2['result']} / {rec2['next']}")

# ------------------------------------------------------------------ (8) 線の定常の保証・壁解像の照合 (合成)
print("== (8) 線の定常の保証 (§5.1 #7a)・壁解像の照合 (§5.1 #7b) (合成)")
from forge_design.report import nozzle_report as NR  # noqa: E402

v0 = WQ.v5c_line(field(), MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.1)
# 監視の線 SL は停止の時の判定に結び付いた基準の line_mach (flow_check の monitored_lines; 量の verdict つき)、flow0 は結び付いた計算の成立
# (2026-10-17 codex result 2 回目 M1 から、結び付き (binding.pass) と量の SETTLED も保証の条件)
SL = {"name": "line_mach", "eta": 0.03, "x_range_m": [4.0, 7.95], "n": 2001, "Md": MD, "verdict": "SETTLED"}
fld0, flow0 = {"step": 46000, "res": "res_46000.h5"}, {"status": "PASS", "stop_step": 46000, "binding": {"pass": True, "problems": []}}
check("(8) 結び付いた監視の線・停止の step と同じ → 保証の欠けなし", WQ.line_monitor_gaps(v0, SL, fld0, flow0) == [])
for nm, sl_, fl_, fw_ in (("η が違う", {**SL, "eta": 0.06}, fld0, flow0), ("x 窓が違う", {**SL, "x_range_m": [4.0, 7.9]}, fld0, flow0),
                          ("点数が違う", {**SL, "n": 1001}, fld0, flow0), ("監視の line_mach が無い", None, fld0, flow0),
                          ("停止の step と違う", SL, fld0, {**flow0, "stop_step": 46171}),
                          ("停止の step の記録が無い", SL, fld0, {"status": "PASS", "binding": {"pass": True}}),
                          ("停止の時の判定と結び付かない (線は同じ)", SL, fld0, {**flow0, "binding": {"pass": False, "problems": ["基準の sha256 が違う"]}}),
                          ("照合の記録が無い (線は同じ)", SL, fld0, {"status": "PASS", "stop_step": 46000}),
                          ("監視の線が停止の時に SETTLED でない", {**SL, "verdict": "NOT_SETTLED"}, fld0, flow0),
                          ("監視の線の量が停止の時の判定に無い", {**SL, "verdict": None}, fld0, flow0)):
    g = WQ.line_monitor_gaps(v0, sl_, fl_, fw_)
    check(f"(8) {nm} → 保証の欠け 1 件", len(g) == 1, f"{g}")


def series_run(td, steps, *, start=0, stale=(), main_started=True):
    """全場の名前だけの run (load_field は差し替える)。"""
    d = Path(td)
    d.mkdir(parents=True, exist_ok=True)
    for s in steps:
        (d / f"res_{s}.h5").write_bytes(b"")
    if main_started:
        (d / "MAIN_STAGE_STARTED").write_text(json.dumps({"start_step": start, "res_files_present": list(stale)}))
    return d


def with_fields(fn_of_step, f):
    """nozzle_report.load_field を step → 合成の場に差し替えて f() を呼ぶ。"""
    orig = NR.load_field
    NR.load_field = lambda run_dir, res: fn_of_step(int(res[len("res_"):-len(".h5")]))
    try:
        return f()
    finally:
        NR.load_field = orig


steps8 = [38000, 40000, 42000, 44000, 46000]
with tempfile.TemporaryDirectory() as td:
    # 定常: 全場が同じ → 系列 STEADY → 線だけの合否のまま (監視の線と違っても)
    run = series_run(Path(td) / "steady", steps8)
    v = WQ.v5c_line(field(), MD, 0.06, [4.0, 7.95], 2001, 0.1, 0.1)
    with_fields(lambda s: field(), lambda: WQ.line_guarantee(v, SL, fld0, flow0, run, MD, series_dir=run))
    ser = v["steadiness"]["series"]
    check("(8) 監視されていない線 (η 0.06)・本段の全場 5 枚が定常 → check_quasisteady STEADY で合否 (PASS) に進む、系列の CSV と VERDICT を記録",
          v["verdict"] == "PASS" and ser["status"] == "STEADY" and ser["verdicts"] == {"maxabs_pct": "STEADY", "halfamp_pct": "STEADY"}
          and ser["steps"] == steps8 and Path(ser["csv"]).is_file() and ser["csv_sha256"] == sha(ser["csv"]) and ser["returncode"] == 0
          and "line_verdict" not in v, f"{ser.get('verdicts')} {ser.get('why')}")
    # 動いている: 傾きが step ごとに 20 % ずつ増える → DRIFTING → 判定不能 (線だけの合否 PASS は line_verdict に残る)
    run = series_run(Path(td) / "drift", steps8)
    v = WQ.v5c_line(field(), MD, 0.06, [4.0, 7.95], 2001, 0.1, 0.1)
    with_fields(lambda s: field(b_pct_per_m=0.02 * (1.0 + 0.2 * (s - 38000) / 2000)),
                lambda: WQ.line_guarantee(v, SL, fld0, flow0, run, MD))
    ser = v["steadiness"]["series"]
    check("(8) 監視されていない線の系列が動いている → check_quasisteady が STEADY でない → V5c 判定不能 (線だけの合否 PASS を残す)",
          v["verdict"] == "INDETERMINATE" and v["line_verdict"] == "PASS" and ser["status"] == "NOT_STEADY" and ser["returncode"] != 0
          and ser["csv"] is None and len(ser["values"]) == 5, f"{ser.get('verdicts')}")
    # 全場が足りない (本段の開始より前の全場・開始の時に残っていた全場・評価の場より後の全場は数えない)
    run = series_run(Path(td) / "few", [1000, 44000, 45000, 46000, 48000], start=2000, stale=["res_45000.h5"])
    v = WQ.v5c_line(field(), MD, 0.06, [4.0, 7.95], 2001, 0.1, 0.1)
    with_fields(lambda s: field(), lambda: WQ.line_guarantee(v, SL, fld0, flow0, run, MD))
    ser = v["steadiness"]["series"]
    check("(8) 本段の全場が 2 枚 (開始前・開始の時に残っていた・評価の場より後を除く) < 4 → 判定不能 (系列を作れない)",
          v["verdict"] == "INDETERMINATE" and ser["status"] == "INDETERMINATE" and ser["steps"] == [44000, 46000], f"{ser.get('why')}")
    run = series_run(Path(td) / "nostart", steps8, main_started=False)
    v = WQ.v5c_line(field(), MD, 0.06, [4.0, 7.95], 2001, 0.1, 0.1)
    with_fields(lambda s: field(), lambda: WQ.line_guarantee(v, SL, fld0, flow0, run, MD))
    check("(8) 本段の開始の記録が無い → 系列の区間を特定できず判定不能", v["verdict"] == "INDETERMINATE"
          and "MAIN_STAGE_STARTED" in v["steadiness"]["series"]["why"])
    run = series_run(Path(td) / "noend", [38000, 40000, 42000, 44000])
    v = WQ.v5c_line(field(), MD, 0.06, [4.0, 7.95], 2001, 0.1, 0.1)
    with_fields(lambda s: field(), lambda: WQ.line_guarantee(v, SL, fld0, flow0, run, MD))
    check("(8) 系列が評価の場の step で終わらない → 判定不能", v["verdict"] == "INDETERMINATE" and "終わらない" in v["steadiness"]["series"]["why"])
    # 監視の線と同じで停止の step → 系列を作らず今どおり
    run = series_run(Path(td) / "mon", [])
    v = WQ.v5c_line(field(), MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.1)
    WQ.line_guarantee(v, SL, fld0, flow0, run, MD)
    check("(8) 監視の線と同じ・停止の step → 今どおり (系列なし、保証の根拠は停止規則の監視)", v["verdict"] == "PASS" and "series" not in v["steadiness"])
    # 停止の step と違う場 (線は同じ) も系列が要る
    v = WQ.v5c_line(field(), MD, 0.03, [4.0, 7.95], 2001, 0.1, 0.1)
    WQ.line_guarantee(v, SL, {"step": 44000}, flow0, run, MD)
    check("(8) 評価の場の step が停止の step と違う (線は同じ) → 系列が作れず判定不能", v["verdict"] == "INDETERMINATE" and v["line_verdict"] == "PASS")

# 壁解像の照合: ツールの出力を差し替えて (偽のツール) 呼び出し側の照合だけを見る
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    run = td / "run_w"
    run.mkdir()
    (run / "bcondConfig.yaml").write_text("inlet: {physID: 1, kind: inlet_Pressure}\nwall: {physID: 3, kind: wall_isothermal}\n"
                                          "wall2: {physID: 5, kind: wall_adiabatic}\naxis: {physID: 4, kind: axis}\n")
    tools = td / "tools"
    tools.mkdir()
    orig_tools = WQ.TOOLS
    WQ.TOOLS = tools

    def fake_wr(out, rc=0, step=46171, stop=46171):
        (tools / "check_wall_resolution.py").write_text(f"import sys\nprint({out!r})\nsys.exit({rc})\n")
        return WQ.wall_resolution(run, step, None, stop_step=stop)

    PASS_TXT = "最大 y1+ = 0.999 ; **目標 1 を超える面積割合 最大 0.0 % (許容 2 %)**\nVERDICT: PASS (壁解像)"
    try:
        r = fake_wr("  wall step 46171\n\nEVALUATED: wall/3@46171, wall2/5@46171\n" + PASS_TXT)
        check("(8) 壁解像: 全部の壁 (bcond の wall*)・要求の step・停止の step で評価した PASS → PASS",
              r["verdict"] == "PASS" and r["evaluated"] == [{"name": "wall", "physID": 3, "step": 46171}, {"name": "wall2", "physID": 5, "step": 46171}]
              and r["groups"] == ["wall", "wall2"], r.get("why"))
        for nm, out, kw in (("EVALUATED 行が無い (照合できない、旧いツールの出力)", PASS_TXT, {}),
                            ("要求した壁 wall2 を評価していない", "EVALUATED: wall/3@46171\n" + PASS_TXT, {}),
                            ("要求した step 46170 でない出力 (46171) を評価した (最新への代え)", "EVALUATED: wall/3@46171, wall2/5@46171\n" + PASS_TXT,
                             {"step": 46170, "stop": None}),
                            ("bcond と違う physID の出力を評価した", "EVALUATED: wall/9@46171, wall2/5@46171\n" + PASS_TXT, {}),
                            ("要求していない壁を評価した", "EVALUATED: wall/3@46171, wall2/5@46171, extra/7@46171\n" + PASS_TXT, {}),
                            ("評価した step 46000 が停止の step 46171 でない", "EVALUATED: wall/3@46000, wall2/5@46000\n" + PASS_TXT,
                             {"step": 46000})):
            r = fake_wr(out, **kw)
            check(f"(8) 壁解像: ツールは PASS だが{nm} → 判定不能", r["verdict"] == "INDETERMINATE" and r["tool_verdict"] == "PASS", r.get("why"))
        r = fake_wr("EVALUATED: wall/3@46171\n" + PASS_TXT.replace("PASS (壁解像)", "FAIL (目標超過)"), rc=1)
        check("(8) 壁解像: ツールの FAIL でも壁がそろわなければ判定不能", r["verdict"] == "INDETERMINATE" and r["tool_verdict"] == "FAIL")
        r = fake_wr("VERDICT: INDETERMINATE (判定不能: 対象の壁・step の出力がそろわない)", rc=2)
        check("(8) 壁解像: ツールの判定不能はそのまま判定不能", r["verdict"] == "INDETERMINATE" and "tool_verdict" not in r)
    finally:
        WQ.TOOLS = orig_tools

# ------------------------------------------------------------------ (9) 実データ run_0014 (codex result M1・M2 の再現)
print("== (9) case/71 run_0014_chain_v2 (写しで; codex result 2026-10-17 の 1・2)")
R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", C71 / "run_0014_chain_v2"))
PROB14 = C71 / "problem_m55_p5_settle_v2.yaml"
need = [R14 / "prepare_info.json", R14 / "nozzle.h5", R14 / "res_46171.h5", R14 / "res_wall_3_46171.h5", R14 / "res_wall_3_46000.h5",
        R14 / "deliverables" / "MANIFEST.json", PROB14, CARD]
if not all(p.is_file() for p in need):
    SKIPPED.append(f"(9) run_0014 ({R14} に run か成果物が無い)")
else:
    import shutil
    import subprocess
    def copy14(dst):
        """run_0014 の写し (h5・xmf は symlink; run_settle_realdata_tests.copy_run と同じ) — 元の run には書かない。"""
        run = Path(dst) / R14.name
        run.mkdir(parents=True)
        for f in R14.iterdir():
            if f.is_dir():
                continue
            if f.suffix in (".h5", ".xmf"):
                (run / f.name).symlink_to(f.resolve())
            else:
                shutil.copy2(f, run / f.name)
        return run

    with tempfile.TemporaryDirectory() as td:
        run = copy14(td)
        common = dict(deliverables=R14 / "deliverables", wall_txt=Path(td) / "wr.txt", series_dir=Path(td))
        rec = WQ.qualify(run, CARD, PROB14, **common)
        c = rec["checks"]
        check("(9) 正しい呼び方 (停止の step 46171・監視の線 η 0.03) は今どおり accepted・V5c PASS (停止規則の監視)・壁解像 PASS (wall/3@46171)",
              rec["result"] == "accepted" and c["V5c"]["verdict"] == "PASS" and "series" not in c["V5c"]["steadiness"]
              and c["wall_resolution"]["verdict"] == "PASS" and c["wall_resolution"]["evaluated"] == [{"name": "wall", "physID": 3, "step": 46171}]
              and round(c["V5c"]["maxabs_pct"], 6) == 0.062305 and round(c["V5c"]["halfamp_pct"], 6) == 0.057285
              and c["deliverables"]["provenance"]["status"] == "MATCH",
              f"{rec['result']} / {rec['missing']}")
        check("(9) 正しい呼び方: 計算の成立は停止の時の判定 (SETTLE_VERDICT.json) と結び付いて PASS、監視の線は停止の時の基準の η 0.03 (SETTLED)",
              c["flow"]["status"] == "PASS" and c["flow"]["binding"]["pass"] and c["flow"]["binding"]["verdict_file"] == "SETTLE_VERDICT.json"
              and c["V5c"]["steadiness"]["monitored_line"]["eta"] == 0.03 and c["V5c"]["steadiness"]["monitored_line"]["verdict"] == "SETTLED",
              f"{c['flow'].get('binding')}")
        # 2026-10-17 codex result 2 回目 M1 の再現: 基準のファイルの読み取り値だけを η 0.06 に変える (メモリ上)。修正前は今の
        # settle_criteria.yaml の読み取り値と評価の線を比べたので、古い停止の記録のまま accepted になった
        orig = WQ._settle_line_from_criteria
        WQ._settle_line_from_criteria = lambda rd: {**orig(rd), "eta": 0.06}
        try:
            rec = WQ.qualify(run, CARD, PROB14, **common)
        finally:
            WQ._settle_line_from_criteria = orig
        c = rec["checks"]
        check("(9) 2 回目 M1 の再現 (メモリ上で基準の読み取り値だけ η 0.06) → 評価の線 η 0.06 は停止の時の監視の線 η 0.03 と違う → V5c 判定不能・"
              "incomplete (修正前は accepted)",
              rec["result"] == "incomplete" and c["V5c"]["eta"] == 0.06 and c["V5c"]["verdict"] == "INDETERMINATE"
              and any("η 0.03" in g for g in c["V5c"]["steadiness"]["monitor_gaps"]) and c["flow"]["status"] == "PASS",
              f"{rec['result']} / {c['V5c'].get('why', '')[:200]}")
    # 基準のファイルを書き換えた写し (η・ほかの値)・停止の時の記録の無い写し: 停止の時の判定と結び付かない → 計算の成立 INDETERMINATE・
    # V5c は監視で保証されず全場 1 枚で系列も作れない → incomplete (修正前はどれも accepted)
    import re as _re
    for nm, edit in (("基準のファイルの line_mach の η を 0.06 に書き換えた写し",
                      lambda r: (r / "settle_criteria.yaml").write_text(_re.sub(r"^  eta: 0\.03$", "  eta: 0.06", (r / "settle_criteria.yaml").read_text(),
                                                                                count=1, flags=_re.M))),
                     ("基準のファイルの流量の収支の閾値 mean_abs を 0.0005 → 0.001 に書き換えた写し (基準としては正しい)",
                      lambda r: (r / "settle_criteria.yaml").write_text((r / "settle_criteria.yaml").read_text().replace("    mean_abs: 0.0005\n",
                                                                                                                      "    mean_abs: 0.001\n", 1))),
                     ("停止の時の判定の記録 (SETTLE_VERDICT.json) の無い写し", lambda r: (r / "SETTLE_VERDICT.json").unlink())):
        with tempfile.TemporaryDirectory() as td:
            run = copy14(td)
            before = sha(run / "settle_criteria.yaml")
            edit(run)
            changed = (not (run / "SETTLE_VERDICT.json").exists()) or sha(run / "settle_criteria.yaml") != before
            rec = WQ.qualify(run, CARD, PROB14, deliverables=R14 / "deliverables", wall_txt=Path(td) / "wr.txt", series_dir=Path(td))
            c = rec["checks"]
            check(f"(9) 2 回目 M1: {nm} → 計算の成立 INDETERMINATE・V5c 判定不能・incomplete・escalate (修正前は accepted)",
                  changed and rec["result"] == "incomplete" and rec["next"] == "escalate" and c["flow"]["status"] == "INDETERMINATE"
                  and not c["flow"]["binding"]["pass"] and c["V5c"]["verdict"] == "INDETERMINATE"
                  and any("結び付かない" in g for g in c["V5c"]["steadiness"]["monitor_gaps"])
                  # 基準を書き換えた写しは基準の sha256 の食い違いで止まる (基準の誤りで止まったのではない)
                  and ("SETTLE_VERDICT" in nm or any("sha256" in p_ and "settle_criteria" in p_ for p_ in c["flow"]["binding"]["problems"])),
                  f"{rec['result']} / {c['flow'].get('why', '')[:160]}")
    # 2 回目 M2 の再現: bcondConfig に壁を足した写し (ダンプが一度も無い)。--groups を省く経路 (受け渡しの関門
    # deltastar_loop.wall_resolution_verdict) は、修正前は存在するダンプ (wall/3) だけで PASS・終了コード 0 だった
    from forge_design.feedback import deltastar_loop as DL  # noqa: E402
    with tempfile.TemporaryDirectory() as td:
        run = copy14(td)
        (run / "bcondConfig.yaml").write_text((run / "bcondConfig.yaml").read_text()
                                              + "missing_wall: {physID: 5, kind: wall_adiabatic, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}\n")
        cwr = str(WQ.TOOLS / "check_wall_resolution.py")
        for args in (["--step", "46171"], []):
            r = subprocess.run([sys.executable, cwr, str(run), *args, "--weight", "area"], capture_output=True, text=True)
            check(f"(9) 2 回目 M2: bcond に一度もダンプの無い壁 missing_wall、--groups 無し {' '.join(args) or '(--step 無し)'} → 判定不能・終了コード 2 "
                  "(修正前は wall/3 だけで PASS・0)",
                  r.returncode == 2 and "VERDICT: INDETERMINATE" in r.stdout and "missing_wall" in r.stdout and "一度も無い" in r.stdout
                  and "EVALUATED" not in r.stdout, r.stdout.strip().splitlines()[-1])
        w = DL.wall_resolution_verdict(run, 46171)
        check("(9) 2 回目 M2: 受け渡しの関門の局所 y₁⁺ (deltastar_loop.wall_resolution_verdict) も UNDETERMINED (修正前は PASS)",
              w["kind"] == "UNDETERMINED" and w["rc"] == 2, f"{w['kind']} {w['verdict']}")
    # 試走 A の写し (case/90; 読むだけ): 停止の時の結び付きは今どおり通る
    for rn in ("run_0003_ns_p0", "run_0005_ns_p1"):
        ra = ROOT / "case/90.qs_trial_a" / rn
        if not (ra / "STOP_RULE_STATE.json").is_file():
            SKIPPED.append(f"(9) 試走 A {rn} ({ra} に run が無い)")
            continue
        f_ = WQ.flow_check(ra)
        check(f"(9) 試走 A {rn}: 計算の成立 PASS (停止の時の判定・基準・入力と結び付く)、監視の線 η 0.03 (SETTLED)",
              f_["status"] == "PASS" and f_["binding"]["pass"] and [(ln["eta"], ln["verdict"]) for ln in f_["monitored_lines"]] == [(0.03, "SETTLED")],
              f"{f_['status']} {f_['binding'].get('problems')}")
    with tempfile.TemporaryDirectory() as td:
        run = copy14(td)
        common = dict(deliverables=R14 / "deliverables", wall_txt=Path(td) / "wr.txt", series_dir=Path(td))
        rec = WQ.qualify(run, CARD, PROB14, cli_tols={"line_eta": 0.06}, **common)
        c = rec["checks"]
        check("(9) M1 の再現: line_eta 0.06 (監視していない線) → V5c 判定不能 (全場 1 枚で系列を作れない)・incomplete・escalate "
              "(修正前は accepted; 線だけの合否 PASS・最大偏差 0.0619 %・片振幅 0.0615 %)",
              rec["result"] == "incomplete" and rec["next"] == "escalate" and c["V5c"]["verdict"] == "INDETERMINATE" and c["V5c"]["line_verdict"] == "PASS"
              and round(c["V5c"]["maxabs_pct"], 4) == 0.0619 and round(c["V5c"]["halfamp_pct"], 4) == 0.0615
              and c["V5c"]["steadiness"]["series"]["steps"] == [46171], f"{rec['result']} / {c['V5c'].get('why')}")
        # 停止の step と違う step の壁解像 (46000 は全場が無く vis_lam を引けないのでツール自身が判定不能になる; 照合は停止の step を
        # 46000 と偽って、ツールの PASS (46171) を呼び出し側が拒むことで見る)
        wr = WQ.wall_resolution(run, 46171, None, stop_step=46000)
        check("(9) 壁解像: ツールは step 46171 で PASS だが停止の step (46000 と仮定) でない → 判定不能",
              wr["verdict"] == "INDETERMINATE" and wr["tool_verdict"] == "PASS" and "停止の step" in wr["why"], wr.get("why"))
        cwr = str(WQ.TOOLS / "check_wall_resolution.py")
        for args, nm in ((["--step", "46170", "--groups", "wall,missing_wall"], "M2 の再現: 無い step・無い壁"),
                         (["--step", "46170", "--groups", "wall"], "無い step (最新の 46171 に代えない)"),
                         (["--step", "46171", "--groups", "wall,missing_wall"], "無い壁の名前")):
            r = subprocess.run([sys.executable, cwr, str(run), *args, "--weight", "area"], capture_output=True, text=True)
            check(f"(9) check_wall_resolution {' '.join(args)} ({nm}) → 判定不能・終了コード 2 (修正前は step 46171 の wall だけで PASS・0)",
                  r.returncode == 2 and "VERDICT: INDETERMINATE" in r.stdout and "EVALUATED" not in r.stdout, r.stdout.strip().splitlines()[-1])
        r = subprocess.run([sys.executable, cwr, str(run), "--step", "46171", "--groups", "wall", "--weight", "area"], capture_output=True, text=True)
        check("(9) check_wall_resolution --step 46171 --groups wall → PASS・0・EVALUATED: wall/3@46171 (今どおり)",
              r.returncode == 0 and "VERDICT: PASS" in r.stdout and "EVALUATED: wall/3@46171" in r.stdout, r.stdout.strip().splitlines()[-1])
    # ---- 3 回目 Major (2026-10-17 codex result): 成果物の生成元と評価の対象の照合 ----
    # 修正前は MANIFEST の中のファイルの sha256 だけを見て、生成元 (run・壁・格子・場) を見なかったので、run_0014 の評価に run_0010 の
    # 成果物 (生成元 run_0010・step 86000・CFD の表は status diagnostic) を渡すと accepted・成果物 PASS だった
    from forge_design.export import deliverables as DV  # noqa: E402
    R10 = C71 / "run_0010_ns_b0"
    OLD14 = R14 / "_deliverables_pre_provenance"
    NEW14 = R14 / "deliverables"

    def link_deliv(src, dst, edit=None):
        """成果物の写し (ファイルは symlink、MANIFEST.json だけ実体で edit で書き換える)。元の成果物には書かない。"""
        for f in sorted(Path(src).rglob("*")):
            if not f.is_file():
                continue
            t = Path(dst) / f.relative_to(src)
            t.parent.mkdir(parents=True, exist_ok=True)
            if f.name == "MANIFEST.json":
                m = json.loads(f.read_text())
                if edit is not None:
                    edit(m)
                t.write_text(json.dumps(m, ensure_ascii=False))
            else:
                t.symlink_to(f.resolve())
        return Path(dst)

    pv14 = json.loads((NEW14 / "MANIFEST.json").read_text()).get("provenance")
    if pv14 is None or not (R10 / "deliverables" / "MANIFEST.json").is_file() or not (OLD14 / "MANIFEST.json").is_file():
        SKIPPED.append(f"(9) 3 回目 Major: 生成元の記録つきの成果物 ({NEW14})・run_0010 の成果物・作り直す前の成果物 ({OLD14}) のどれかが無い")
    else:
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            run = copy14(td / "r")
            common = dict(wall_txt=td / "wr.txt", series_dir=td)
            rec = WQ.qualify(run, CARD, PROB14, deliverables=NEW14, **common)
            dv_ = rec["checks"]["deliverables"]
            check("(9) 3 回目 Major の正例: run_0014 を自分の成果物 (生成元の記録つき) で評価 → 生成元 MATCH・成果物 PASS・accepted",
                  rec["result"] == "accepted" and dv_["status"] == "PASS" and dv_["provenance"]["status"] == "MATCH"
                  and pv14["field"]["step"] == 46171 and pv14["field"]["is_stop_step"] is True, f"{rec['result']} / {dv_.get('why')}")
            rec = WQ.qualify(run, CARD, PROB14, deliverables=link_deliv(NEW14, td / "elsewhere" / "deliv"), **common)
            check("(9) 3 回目 Major: 成果物を別の置き場所に写した (中身は同じ) → 生成元 MATCH・accepted",
                  rec["result"] == "accepted" and rec["checks"]["deliverables"]["provenance"]["status"] == "MATCH", rec["result"])
            err = io.StringIO()
            with redirect_stdout(io.StringIO()), redirect_stderr(err):
                rc = WQ.main([str(run), "--card", str(CARD), "--problem", str(PROB14), "--deliverables", str(R10 / "deliverables"), "--no-write"])
            rec = WQ.qualify(run, CARD, PROB14, deliverables=R10 / "deliverables", **common)
            c = rec["checks"]
            check("(9) 3 回目 Major の再現: run_0014 (停止 46171) の評価に run_0010 の成果物 (生成元の記録が無い古い成果物) → 生成元 ABSENT・"
                  "成果物と寸法 INDETERMINATE・incomplete・escalate・CLI の終了コード 4 (修正前は accepted・0)",
                  rec["result"] == "incomplete" and rec["next"] == "escalate" and rc == 4 and c["deliverables"]["status"] == "INDETERMINATE"
                  and c["deliverables"]["provenance"]["status"] == "ABSENT" and c["deliverables"]["dimensions"]["status"] == "INDETERMINATE"
                  and c["flow"]["status"] == "PASS" and c["V5b"]["verdict"] == "PASS" and c["V5c"]["verdict"] == "PASS"
                  and c["wall_resolution"]["verdict"] == "PASS" and any("作り直す" in m_ for m_ in rec["missing"]),
                  f"{rec['result']} / {rec['next']} / rc {rc}")
            prov10 = DV.provenance_record(R10, 86000, ["res_wall_3_86000.h5"], CARD, PROB14,
                                          field_source="試験: 新しい道具で run_0010 から作ったときと同じ記録")
            rec = WQ.qualify(run, CARD, PROB14, deliverables=link_deliv(R10 / "deliverables", td / "d10",
                                                                        lambda m: m.update(provenance=prov10)), **common)
            pr = rec["checks"]["deliverables"]["provenance"]
            want = ("壁 (wall_repr.json)", "格子 (nozzle.h5)", "設定 (solverConfig.yaml)", "設定 (bcondConfig.yaml)", "場の step", "停止の step")
            check("(9) 3 回目 Major: run_0010 の成果物に新しい道具の生成元の記録 (run_0010・step 86000) を足した写し → 壁・格子・設定 2 つ・"
                  "場の step・停止の step が違う (MISMATCH)・incomplete (カード・問題 YAML は同じ)",
                  rec["result"] == "incomplete" and pr["status"] == "MISMATCH" and all(any(w in p_ for p_ in pr["problems"]) for w in want)
                  and not any(p_.startswith(("カード", "問題 YAML")) for p_ in pr["problems"]), pr["problems"])
            rec = WQ.qualify(run, CARD, PROB14, deliverables=OLD14, **common)
            check("(9) 3 回目 Major: run_0014 自身の作り直す前の成果物 (_deliverables_pre_provenance、生成元の記録が無い) → ABSENT・incomplete "
                  "(修正前は accepted)",
                  rec["result"] == "incomplete" and rec["checks"]["deliverables"]["provenance"]["status"] == "ABSENT"
                  and rec["checks"]["deliverables"]["status_before_provenance"] == "PASS", rec["result"])
            rec = WQ.qualify(run, CARD, PROB14, deliverables=link_deliv(NEW14, td / "d14x",
                                                                        lambda m: m["provenance"]["mesh"].update(sha256="0" * 64)), **common)
            pr = rec["checks"]["deliverables"]["provenance"]
            check("(9) 3 回目 Major: run_0014 の成果物の写しで記録の格子の sha256 を 1 つ変えた → MISMATCH (違いは格子だけ)・incomplete",
                  rec["result"] == "incomplete" and pr["status"] == "MISMATCH" and len(pr["problems"]) == 1 and "格子" in pr["problems"][0],
                  pr["problems"])

for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if not FAIL else f'{len(FAIL)} FAILURES: {FAIL}'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
