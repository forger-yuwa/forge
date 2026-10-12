#!/usr/bin/env python3
"""要求カード (forge.wt_card/v1) の生成器 `forge_design.intake.wt_card` の試験 — plan tooling-nozzle-wt-deliverables §6 V1。

  (a) 正例 2 (物理スロート径 = M5 デモ・設計スロート径)・負例 (未知キー・欠けたキー・R_in ≤ r_t・μ > 20・スロート径の二重指定・
      型違い・負値・組成の和・上書きできないキー ほか)
  (b) 単位換算の往復: 生成した YAML の値 × r_t が実寸と相対 1e-12 (**単位換算の試験であって形状の合格根拠ではない**)
  (c) 固定点: case/45 の生産の実寸 (r_t 0.07666536551630307 m の r_inlet 6.485・L_U 12.0・L_pipe 0.5・設計壁の終点、物理スロート径 =
      生産の prepare_info の物理スロート半径 × 2) のカードで build_problem が 10 回以内に止まる。**--slow のときだけ**
  ほか: M5 デモの生成 YAML が load_problem と design_chain (μ・law の窓・壁 QA) を通り、設計壁から測った寸法がカードと公差内。
  最終の雛形 (axismach_final_v1) の問題も load_problem と design_chain を通る。出口静圧 (semiperfect) を case/45 の記録と照合。
  液化の判定 (未飽和・過飽和・範囲外・N₂ なし)。
  自動調整 (2026-10-13 追加): 縮流部側の第 1 層のブレンドの決め直し (y1+ の目安 ≤ 1・AR ≤ 5000)、x_F が逆 MOC の場の終点に張り付いた
  ときの x_end_margin の引き上げ、from_length が止まったときの内部許容差の明示の緩和 (L_total_tol_relaxed の記録)、
  from_length の各評価の (L_c, x_F, x_e, x_F − x_e) の記録。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_wt_card_tests.py [--slow]
所要: 既定 約 2 分 (M5 デモの r_t の固定点 + design_chain 数回)、--slow で case/45 の固定点が加わる。
"""
import copy
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate.runner_axismach import design_chain  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

FAIL = 0
SLOW = "--slow" in sys.argv[1:]
M5_CARD = DESIGN / "forge_design" / "intake" / "examples" / "m5_air_demo.yaml"
# case/45 の生産 (problem_d155_ns_prod.yaml・生産の run の prepare_info.json・wall_design.csv の値)
C45_RT = 0.07666536551630307
C45_PHYS_THROAT_R = 0.07677180151628452          # prepare_info.json sizing.physical_throat_radius_m
C45_WALL_END_M = 7.300657069598003               # wall_design.csv の最後の x_m (= qa.x_F 95.2276 r_t × r_t)


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def raises(fn, *needles):
    """fn() が ValueError を出し、メッセージに needles を全部含むか。(結果, メッセージ)。"""
    try:
        fn()
    except ValueError as ex:
        msg = str(ex)
        return all(n in msg for n in needles), msg
    return False, "(例外なし)"


BASE = yaml.safe_load(M5_CARD.read_text())


def card(**kw):
    c = copy.deepcopy(BASE)
    for k, v in kw.items():
        if v is None:
            c.pop(k, None)
        else:
            c[k] = v
    return c


# ---------------------------------------------------------------------------------------------- (a) 検査の負例・正例 (軽い)
print("== (a) カードの検査")
ok = True
try:
    v = W.validate_card(BASE)
except ValueError as ex:
    ok, v = False, str(ex)
check("M5 デモのカードが検査を通る (PyYAML が文字列と読む 3.0e6 を数値として受ける)",
      ok and isinstance(v, dict) and v["conditions"]["Pt"] == 3.0e6 and v["conditions"]["Tt"] == 1250.0, str(v)[:200])
check("既定値の補完: p_ambient design・recipe・overrides", ok and v["p_ambient"] == "design" and v["recipe"] == "axismach_explore_v1"
      and v["overrides"] == {})
check("air の preset の質量分率の和が 1", abs(sum(W.GAS_PRESETS["air"].values()) - 1.0) < 1e-12)

neg = [
    ("未知の最上位キー", lambda: W.validate_card(card(diameter_exit_m=0.3)), ("未知のキー", "diameter_exit_m")),
    ("欠けたキー (inlet)", lambda: W.validate_card(card(inlet=None)), ("inlet",)),
    ("欠けた下位キー (inlet.contraction_length_m)",
     lambda: W.validate_card(card(inlet={"radius_m": 0.15, "pipe_length_m": 0.2})), ("contraction_length_m",)),
    ("スロート径の二重指定", lambda: W.validate_card(card(throat={"diameter_physical_m": 0.05, "diameter_design_m": 0.05})),
     ("どちらか 1 つ",)),
    ("型違い (Pt が文字列)", lambda: W.validate_card(card(conditions={"Pt": "3 MPa", "Tt": 1250.0})), ("conditions.Pt", "数値")),
    ("型違い (Tt が真偽値)", lambda: W.validate_card(card(conditions={"Pt": 3.0e6, "Tt": True})), ("conditions.Tt",)),
    ("負値 (pipe_length_m)", lambda: W.validate_card(card(inlet={"radius_m": 0.15, "pipe_length_m": -0.1,
                                                                    "contraction_length_m": 0.3})), ("pipe_length_m", "以上")),
    ("負値 (divergent_length_m)", lambda: W.validate_card(card(divergent_length_m=-1.2)), ("divergent_length_m",)),
    ("M_exit ≤ 1", lambda: W.validate_card(card(M_exit=0.9)), ("M_exit",)),
    ("format の違い", lambda: W.validate_card(card(format="forge.wt_card/v0")), ("format",)),
    ("組成の和が 1 でない", lambda: W.validate_card(card(gas={"species": {"N2": 0.7, "O2": 0.2}})), ("和",)),
    ("未知の preset", lambda: W.validate_card(card(gas={"preset": "argon"})), ("preset",)),
    ("preset と species の併記", lambda: W.validate_card(card(gas={"preset": "air", "species": {"N2": 1.0}})), ("どちらか 1 つ",)),
    ("内蔵 DB に無い種", lambda: W.validate_card(card(gas={"species": {"N2": 0.5, "XE": 0.5}})), ("XE",)),
    ("断熱壁に Tw", lambda: W.validate_card(card(wall={"thermal": "adiabatic", "Tw": 300.0})), ("adiabatic",)),
    ("未知の雛形", lambda: W.validate_card(card(recipe="axismach_v9")), ("recipe",)),
    ("上書きできないキー (geometry.r_inlet)", lambda: W.validate_card(card(overrides={"geometry": {"r_inlet": 7.0}})),
     ("上書きできない",)),
    ("上書きできない節 (spec)", lambda: W.validate_card(card(overrides={"spec.Pt": 1.0e6})), ("上書きできない",)),
    ("R_in ≤ r_t", lambda: W.card_to_problem(card(inlet={"radius_m": 0.02, "pipe_length_m": 0.2, "contraction_length_m": 0.3}),
                                              0.025), ("R_in > r_t",)),
    # plan §4.1 の例の縮流部長 0.40 m は R 2.0・r_t 25 mm で μ = 16²/(2·5) = 25.6 > 20 (仕様書の M5 デモは 0.30 m・μ 14.4)
    ("μ > 20 (縮流部長 0.40 m)", lambda: W.card_to_problem(card(inlet={"radius_m": 0.15, "pipe_length_m": 0.2,
                                                                       "contraction_length_m": 0.40}), 0.025), ("μ", "> 20")),
    ("r_t が負", lambda: W.card_to_problem(card(), -0.025), ("r_t",)),
]
for name, fn, needles in neg:
    r, msg = raises(fn, *needles)
    check(f"負例: {name}", r, msg[:160])
# 設計スロート径のカードに --r-throat (重い処理の前に止まる)
with tempfile.TemporaryDirectory() as td:
    cp = Path(td) / "card.yaml"
    cp.write_text(yaml.safe_dump(card(throat={"diameter_design_m": 0.05})))
    r, msg = raises(lambda: W.build_problem(cp, Path(td) / "p.yaml", r_throat_m=0.0251, log=lambda *a, **k: None), "--r-throat")
    check("負例: 設計スロート径のカードに --r-throat", r, msg[:160])

# ---------------------------------------------------------------------------------------------- 出口静圧・液化
print("== 出口静圧 (semiperfect) と液化の事前警告")
c45 = W.validate_card({"format": W.CARD_FORMAT, "name": "c45_probe",
                       "gas": {"species": {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.0220, "N2": 0.7245}},
                       "conditions": {"Pt": 5.5e6, "Tt": 1600.0}, "M_exit": 6.0,
                       "throat": {"diameter_physical_m": 2.0 * C45_PHYS_THROAT_R},
                       "inlet": {"radius_m": 6.485 * C45_RT, "pipe_length_m": 0.5 * C45_RT, "contraction_length_m": 12.0 * C45_RT},
                       "divergent_length_m": C45_WALL_END_M, "wall": {"thermal": "adiabatic"}})
es = W.exit_state(c45)
# case/45 の生産 YAML のコメント: P/Pt(M6, semi-perfect Tt1600) = 4.067e-4 → p_ambient 2237 Pa、出口 T 233.5 K
check("出口静圧: case/45 の記録 P/Pt 4.067e-4 と 4 桁で一致", abs(es["p_over_Pt"] - 4.067e-4) < 0.5e-7, f"{es['p_over_Pt']:.6e}")
check("出口静温: case/45 の記録 233.5 K と 0.05 K 以内", abs(es["T_K"] - 233.5) < 0.05, f"{es['T_K']:.4f}")
m5 = W.validate_card(BASE)
es5 = W.exit_state(m5)
check("M5 デモの出口 (親の選定 221.8 K・5139 Pa)", abs(es5["T_K"] - 221.8) < 0.05 and abs(es5["p_Pa"] - 5139.0) < 1.0,
      f"{es5['T_K']:.3f} K, {es5['p_Pa']:.1f} Pa")
lq = W.liquefaction_check(m5)
check("液化: M5 デモ (出口 222 K) は範囲の上で液化しない・警告なし", lq["supersaturated"] is False and lq["warning"] is None, lq["verdict"])
cold = {k: v for k, v in BASE.items()}
lq2 = W.liquefaction_check(W.validate_card({**cold, "conditions": {"Pt": 2.0e7, "Tt": 400.0}}))
check("液化: Tt 400 K・Pt 20 MPa (出口 66.6 K) は過飽和で警告", lq2["supersaturated"] is True and lq2["warning"] is not None, lq2["verdict"])
lq3 = W.liquefaction_check(W.validate_card({**cold, "conditions": {"Pt": 1.0e6, "Tt": 400.0}}))
check("液化: Tt 400 K・Pt 1 MPa は未飽和", lq3["supersaturated"] is False and lq3["warning"] is None, lq3["verdict"])
lq4 = W.liquefaction_check(W.validate_card({**cold, "conditions": {"Pt": 1.0e6, "Tt": 370.0}}))
check("液化: 出口 61.6 K は Antoine の範囲外 → 判定できない警告", lq4["supersaturated"] is None and lq4["warning"] is not None, lq4["verdict"])
lq5 = W.liquefaction_check(W.validate_card({**cold, "gas": {"species": {"AR": 1.0}}}))
check("液化: N₂ を含まない組成は判定しない", lq5["x_N2"] is None and "N2 を含まない" in lq5["verdict"], lq5["verdict"])
r, msg = raises(lambda: W.exit_state(W.validate_card({**cold, "conditions": {"Pt": 1.0e6, "Tt": 300.0}})), "表の上限")
check("出口が semiperfect の表の外 (Tt 300 K で M5) は例外", r, msg[:120])

# ---------------------------------------------------------------------------------------------- (b) 単位換算の往復
print("== (b) 単位換算の往復")
for rt in (0.025, 0.024969004795523574, 0.0251234567891):
    prob = W.card_to_problem(m5, rt)
    with tempfile.TemporaryDirectory() as td:
        pth = Path(td) / "p.yaml"
        W._write_problem(m5, prob, pth, ["test"])
        back = yaml.safe_load(pth.read_text())
    g, s = back["geometry"], back["spec"]
    rt_b = float(s["r_throat"])
    pairs = [("r_inlet × r_t = inlet.radius_m", g["r_inlet"] * rt_b, m5["inlet"]["radius_m"]),
             ("L_U × r_t = contraction_length_m", g["L_U"] * rt_b, m5["inlet"]["contraction_length_m"]),
             ("L_pipe × r_t = pipe_length_m", g["L_pipe"] * rt_b, m5["inlet"]["pipe_length_m"]),
             ("L_total × r_t = divergent_length_m", back["dv"]["L_total"]["value"] * rt_b, m5["divergent_length_m"]),
             ("2 L_total_tol × r_t = length_m", 2.0 * g["L_total_tol"] * rt_b, m5["tolerances"]["length_m"]),
             ("2 sizing.target_m = diameter_physical_m", 2.0 * s["sizing"]["target_m"], m5["throat"]["diameter_physical_m"])]
    rel = max(abs(a - b) / abs(b) for _, a, b in pairs)
    check(f"往復 r_t {rt!r}: 長さ 6 量が相対 1e-12", rel <= 1e-12 and rt_b == rt and back == prob, f"最大 {rel:.2e}")
check("往復: Lc_mode from_length・M_knot = 0.45 M_exit・p_ambient = 出口静圧",
      back["geometry"]["Lc_mode"] == "from_length" and abs(back["geometry"]["M_knot"] - 2.25) < 1e-15
      and abs(back["spec"]["p_ambient"] - es5["p_Pa"]) < 1e-9)
# overrides の記録 (雛形の値 → 上書き)、μ は上書き後の R で検査する
po = W.card_to_problem(W.validate_card(card(overrides={"geometry": {"R": 2.5}, "mesh.nj": 81, "geometry.M_knot": 2.4})), 0.025)
ov = {o["key"]: o for o in po["_card"]["overrides"]}
check("overrides: 上書きと記録 (雛形の値つき)", po["geometry"]["R"] == 2.5 and po["mesh"]["nj"] == 81 and po["geometry"]["M_knot"] == 2.4
      and ov["geometry.R"]["before"] == 2.0 and ov["mesh.nj"]["before"] == 97 and abs(po["_card"]["mu_contraction"] - 144 / 12.5) < 1e-12,
      str(ov)[:200])
r, msg = raises(lambda: W.card_to_problem(W.validate_card(card(overrides={"geometry.R": 1.2})), 0.025), "μ")
check("overrides の R で μ > 20 になれば例外 (μ = 144/6 = 24)", r, msg[:120])

# ---------------------------------------------------------------------------------------------- 張り付きの判定・停止の言い方 (軽い)
print("== 逆 MOC の場の終点への張り付きの判定 (合成の配列)")
wf = np.c_[np.arange(0.0, 10.0001, 0.1), np.ones(101), np.zeros(101)]
r_in = W.moc_call_record({"wall_full": wf, "wall": wf[wf[:, 0] <= 5.0], "exit": {"x_F": 5.02}}, {"x_E": 3.0})
r_end = W.moc_call_record({"wall_full": wf, "wall": wf, "exit": {"x_F": 9.95}}, {"x_E": 3.0})
r_at = W.moc_call_record({"wall_full": wf, "wall": wf, "exit": {"x_F": 10.0}}, {"x_E": 3.0})
check("場の内側の x_F (終点まで 4.98) は張り付きでない", r_in["pinned"] is False and abs(r_in["margin_to_field_end_rt"] - 4.98) < 1e-12
      and r_in["n_beyond_x_F"] == 50 and abs(r_in["x_F_minus_x_e_rt"] - 0.02) < 1e-12)
check("終点から 1 点間隔以内 (0.05 ≤ 0.1) と終点ちょうどは張り付き", r_end["pinned"] and r_at["pinned"])
r_drop = W.moc_call_record({"wall_full": wf, "wall": wf[wf[:, 0] <= 5.0 - 1e-14], "exit": {"x_F": 5.0 - 1e-14}}, {"x_E": 3.0})
check("記録: x_F が壁点のわずか下で切る規則がその点を落とした (x_e が 1 点間隔短い) ことを残す",
      r_drop["cut_dropped_node_at_x_F"] and abs(r_drop["x_F_minus_x_e_rt"] - 0.1) < 1e-9 and not r_in["cut_dropped_node_at_x_F"])
check("停止の言い方: 「改善しなくなって止まった」「ノイズ床とは区別できない」", "改善しなくなって止まった" in W.STALL_TEXT
      and "区別できない" in W.STALL_TEXT and "runner_axismach.py:606-616" in W.STALL_TEXT)
src = (DESIGN / "forge_design" / "intake" / "wt_card.py").read_text()
check("生成器の文言に「ノイズ床で止まる」「ノイズ床に到達」の断定が無い", "ノイズ床で" not in src and "ノイズ床に到達" not in src)

# ---------------------------------------------------------------------------------------------- M5 デモ (物理スロート径・正例 1)
print("== (a) 正例 1: M5 デモ (物理スロート径) の build_problem")
with tempfile.TemporaryDirectory() as td:
    out = Path(td) / "problem_m5_air_demo.yaml"
    logd = W.build_problem(M5_CARD, out)
    text = out.read_text()
    hist = logd["r_t"]["history"]
    rt = logd["r_t"]["r_t_m"]
    check("r_t の固定点が 10 回以内・|Δr_t|/r_t ≤ 1e-6", 1 <= len(hist) <= 10 and hist[-1]["rel_change"] <= 1e-6,
          "; ".join(f"{h['r_t_in_m']!r}→{h['r_t_out_m']!r} ({h['rel_change']:.2e})" for h in hist))
    p = load_problem(out)
    check("生成 YAML が load_problem を通る", p.type == "wind_tunnel_axisym_axismach" and float(p.spec["r_throat"]) == rt)
    raw = yaml.safe_load(text)
    check("_card 節 (最上位) が残り load_problem が拒否しない", "_card" in raw and raw["_card"]["format"] == W.CARD_FORMAT
          and raw["_card"]["r_t_method"] == "physical_throat_solved")
    check("spec.sizing {method: throat, target_m: 径/2}", raw["spec"]["sizing"] == {"method": "throat", "target_m": 0.025})
    check("mesh の ni・wall_first_frac・_throat を埋めた (mesh_euler.ni = mesh.ni)",
          isinstance(raw["mesh"]["ni"], int) and raw["mesh"]["wall_first_frac"] > 0 and raw["mesh"]["wall_first_frac_throat"] > 0
          and raw["mesh_euler"]["ni"] == raw["mesh"]["ni"], str({k: raw["mesh"][k] for k in W._MESH_ESTIMATED}))
    check("見積もりは初期値と YAML のコメント・戻り値に書く", "初期値" in text and "初期値" in logd["mesh_estimate"]["note"])
    e = logd["mesh_estimate"]
    check("格子の見積もり: 実際の格子の AR ≤ ar_max 5000", e["mesh_check"]["ar_max"] <= 5000.0, f"{e['mesh_check']['ar_max']:.1f}")
    check("card_log.json (対応表・警告) を書いた", Path(str(out) + ".card_log.json").is_file() and len(logd["mapping"]) >= 12)
    # 生成 YAML を独立に design_chain へ (μ・law の窓・壁 QA は design_chain の中の例外で検査される)
    d = design_chain(p)
    Wd = d["wall"]
    check("生成 YAML が design_chain を通る (law の窓・壁 QA・壁の検査)", not d["qa"]["violations"] and not d["gates"]["violations"]
          and Wd.up.mu <= 20.0, f"μ {Wd.up.mu:.4g}, L_c {d['L_c']:.6g}, 窓 {d['Lc_window']}")
    # 測った寸法 (§4.2a のゲートの前段): 設計壁から
    tl, tdm = BASE["tolerances"]["length_m"], BASE["tolerances"]["diameter_m"]
    x_in, xU, x_e = float(Wd.x_in), -float(Wd.up.L_U), float(Wd.x_e)
    xs = np.linspace(x_in, x_e, 200001)
    rs = Wd.r(xs)
    k_t = int(np.argmin(rs))
    meas = {"入口半径": (float(rs[0]) * rt, BASE["inlet"]["radius_m"], tdm / 2.0),
            "直管長": ((xU - x_in) * rt, BASE["inlet"]["pipe_length_m"], tl),
            "縮流部長 (接合 → r の最小)": ((float(xs[k_t]) - xU) * rt, BASE["inlet"]["contraction_length_m"], tl),
            "スロート→出口長 (設計壁の終点)": (x_e * rt, BASE["divergent_length_m"], tl)}
    for nm, (mv, cv, tol) in meas.items():
        check(f"寸法: {nm} が公差内", abs(mv - cv) <= tol, f"測定 {mv:.6f} m、カード {cv} m、差 {mv - cv:+.3e}、公差 {tol:g}")
    pr = np.abs(Wd.r(np.linspace(x_in, xU, 2001)) - float(rs[0])).max()
    check("寸法: 直管区間で r が一定 (|r − r_in| ≤ 1e-12 r_t)", pr <= 1e-12, f"{pr:.2e}")
    phys = logd["r_t"]["sizing_result"]["physical_throat"]["r_m"]
    check("寸法: 物理スロート径 (CFD 前の δ) が公差内", abs(2.0 * phys - BASE["throat"]["diameter_physical_m"]) <= tdm,
          f"{2.0 * phys!r} m")
    check("寸法: 生成器の判定も PASS", logd["dimensions"]["pass"])
    # 縮流部側のブレンドの決め直し (y1+ の目安 ≤ 1、AR ≤ 5000 から ni)
    bc = e["blend_candidates"]
    check("ブレンド: 雛形の値では y1+ の目安 > 1 なので決め直した", bc[0]["label"] == "recipe" and bc[0]["y1plus_est_max"] > 1.0
          and e["blend_choice"] != "recipe", str([(b["label"], round(b["y1plus_est_max"], 3), round(b["x_y1plus_est_max_rt"], 2)) for b in bc]))
    check("ブレンド: 決め直した配点で y1+ の目安の最大 ≤ 1 (目安 0.8)", e["y1plus_est_max_with_blend"] <= 1.0,
          f"{e['y1plus_est_max_with_blend']:.3f} @ x {e['x_y1plus_est_max_rt']:.2f}")
    check("ブレンド: YAML の mesh の up_x0/x1 が選んだ値", raw["mesh"].get("wall_first_up_x0", "×") == e["wall_first_up_x0"]
          and raw["mesh"].get("wall_first_up_x1", "×") == e["wall_first_up_x1"] and "決め直した" in text)
    check("ブレンド: 選んだ配点で ni を見積もり直し、実際の格子の AR ≤ 5000",
          e["mesh_check"]["ar_max"] <= 5000.0 and e["ni"] >= e["ni_required_for_ar"], f"ni {e['ni']}, AR {e['mesh_check']['ar_max']:.1f}")
    # 逆 MOC の場と各評価の記録
    f = logd["moc_field"]
    check("場: x_F は場の終点に張り付いていない (切る前の壁が x_F より下流に続く)", f["pinned"] is False
          and f["x_wall_full_end"] > f["x_F"] + f["spacing_end"] and raw["_card"]["moc_field"]["pinned"] is False,
          f"x_F {f['x_F']:.6f}, 場の終点 {f['x_wall_full_end']:.4f}, 下流の壁点 {f['n_beyond_x_F']}")
    ev_all = [ev for m in logd["moc_evals"] for ev in m["evals"]]
    need = ("L_c", "x_F_rt", "x_e_rt", "x_F_minus_x_e_rt", "x_F_minus_x_e_m")
    check("card_log: from_length の各評価の (L_c, x_F, x_e, x_F − x_e [r_t, m]) を記録",
          len(logd["moc_evals"]) >= len(hist) + 1 and ev_all and all(all(ev.get(k) is not None for k in need) for ev in ev_all)
          and all(abs(ev["x_F_minus_x_e_m"] - ev["x_F_minus_x_e_rt"] * m["r_t_m"]) < 1e-15 for m in logd["moc_evals"] for ev in m["evals"]),
          f"design_chain {len(logd['moc_evals'])} 回、評価 {len(ev_all)} 回")
    check("最終の x_e は設計壁の終点と一致", abs(logd["moc_field"]["x_e"] - x_e) == 0.0)
    print(f"   r_t の反復: {[(h['r_t_in_m'], h['r_t_out_m'], h['rel_change'], round(h['time_s'], 1)) for h in hist]}")
    print(f"   格子 (初期値): ni {e['ni']} (AR に要る {e['ni_required_for_ar']})、wall_first_frac {e['wall_first_frac']} "
          f"(生 {e['wall_first_frac_raw']:.3e} @ x {e['x_wall_first_frac_rt']:.2f})、_throat {e['wall_first_frac_throat']} "
          f"(生 {e['wall_first_frac_throat_raw']:.3e} @ x {e['x_wall_first_frac_throat_rt']:.2f}、y1 {e['y1_m_throat'] * 1e9:.1f} nm)、"
          f"AR {e['mesh_check']['ar_max']:.0f}、y1+ 目安の最大 {e['y1plus_est_max_with_blend']:.2f}")
    print(f"   液化: {logd['liquefaction']['verdict']}")
    print(f"   場: x_F {f['x_F']!r}、x_e {f['x_e']!r}、x_F − x_e {f['x_F_minus_x_e_m']!r} m、場の終点 {f['x_wall_full_end']!r} r_t "
          f"(x_end_margin {f['x_end_margin']})、拡大部長の差 {[r['diff_m'] for r in logd['dimensions']['rows'] if r['name'] == 'divergent_length'][0]!r} m")
    evs = [(round(ev["L_c"], 6), round(ev["x_F_rt"], 6), round(ev["x_e_rt"], 6), round(ev["x_F_minus_x_e_m"] * 1e3, 4),
            ev["cut_dropped_node_at_x_F"]) for ev in ev_all]
    print(f"   各評価 (L_c, x_F, x_e, x_F − x_e [mm], 切る規則が x_F の壁点を落としたか): {evs}")
    bcs = [(b["label"], b["wall_first_up_x0"], b["wall_first_up_x1"], b["ni"], round(b["y1plus_est_max"], 3), round(b["x_y1plus_est_max_rt"], 2))
           for b in bc]
    print(f"   ブレンド候補 (label, up_x0, up_x1, ni, y1+ の目安の最大, その x): {bcs}")
    print(f"   警告: {logd['warnings']}")
    rt_m5 = rt

# ---------------------------------------------------------------------------------------------- 設計スロート径 (正例 2、CLI 経由)
print("== (a) 正例 2: 設計スロート径 (CLI)")
with tempfile.TemporaryDirectory() as td:
    cp = Path(td) / "m5_design.yaml"
    cd = card(name="m5_air_demo_design", throat={"diameter_design_m": 0.050})
    cp.write_text(yaml.safe_dump(cd, sort_keys=False))
    out = Path(td) / "problem_design.yaml"
    rc = W.main([str(cp), "--out", str(out)])
    logd = json.loads(Path(str(out) + ".card_log.json").read_text())
    raw = yaml.safe_load(out.read_text())
    check("CLI の終了コード 0 (寸法 PASS)", rc == 0)
    check("r_t = 設計スロート径/2 で 1 回 (反復なし、spec.sizing なし)",
          raw["spec"]["r_throat"] == 0.025 and logd["r_t"]["method"] == "design_throat" and not logd["r_t"]["history"]
          and "sizing" not in raw["spec"])
    rows = {r["name"]: r for r in logd["dimensions"]["rows"]}
    check("寸法: 設計スロート径 (設計壁の最小半径 × 2) が公差内", rows["throat_diameter_design"]["pass"],
          f"{rows['throat_diameter_design']['measured_m']!r}")
    check("寸法: 生成器の判定 PASS", logd["dimensions"]["pass"], str([(r["name"], r["diff_m"]) for r in logd["dimensions"]["rows"]]))
    load_problem(out)

# ---------------------------------------------------------------------------------------------- 自動調整: x_end_margin・内部許容差
print("== 自動調整: 逆 MOC の場が足りない (x_end_margin 1.0) → 引き上げ")
with tempfile.TemporaryDirectory() as td:
    cp = Path(td) / "m5_short.yaml"
    cp.write_text(yaml.safe_dump(card(name="m5_short_field", throat={"diameter_design_m": 0.050},
                                      overrides={"geometry.x_end_margin": 1.0}), sort_keys=False))
    out = Path(td) / "problem_short.yaml"
    logd = W.build_problem(cp, out)
    aa = [a for a in logd["auto_adjust"] if a["key"] == "geometry.x_end_margin"]
    raw = yaml.safe_load(out.read_text())
    check("x_end_margin を 0.5 ずつ上げて張り付き (終端特性線が壁に届かない) を解消", aa and aa[0]["from"] == 1.0
          and all(abs(a["to"] - a["from"] - 0.5) < 1e-12 for a in aa) and raw["geometry"]["x_end_margin"] == aa[-1]["to"]
          and logd["moc_field"]["pinned"] is False and raw["_card"]["moc_field"]["x_end_margin"] == aa[-1]["to"],
          str([(a["from"], a["to"]) for a in aa]) + f" 場の終点まで {logd['moc_field']['margin_to_field_end_rt']:.3g} r_t")
    check("引き上げの履歴を YAML のコメントに書いた", "自動調整: geometry.x_end_margin" in out.read_text())
    print(f"   x_end_margin の履歴: {[(a['from'], a['to'], a['reason'][:60]) for a in aa]}")

print("== 自動調整: from_length が止まる (長さ公差 1e-6 m) → n_axis_inv 2 倍 → 内部許容差 0.05 r_t に緩める")
with tempfile.TemporaryDirectory() as td:
    cp = Path(td) / "m5_tight.yaml"
    cp.write_text(yaml.safe_dump(card(name="m5_tight_length", throat={"diameter_design_m": 0.050},
                                      tolerances={"length_m": 1.0e-6, "diameter_m": 2.0e-5}), sort_keys=False))
    out = Path(td) / "problem_tight.yaml"
    rc = W.main([str(cp), "--out", str(out)])
    logd = json.loads(Path(str(out) + ".card_log.json").read_text())
    text = out.read_text()
    raw = yaml.safe_load(text)
    rl = logd["L_total_tol_relaxed"]
    keys = [a["key"] for a in logd["auto_adjust"]]
    check("緩和を card_log の L_total_tol_relaxed {from, to, reason} に記録 (from = 1e-6/2/r_t、to = 0.05)",
          rl is not None and abs(rl["from"] - 1.0e-6 / 2 / 0.025) < 1e-15 and rl["to"] == 0.05 and W.STALL_TEXT in rl["reason"],
          str(rl)[:200])
    check("順序: n_axis_inv を 2 倍 → 内部許容差を緩めて n_axis_inv を雛形に戻す",
          keys[:3] == ["geometry.n_axis_inv", "geometry.L_total_tol", "geometry.n_axis_inv"] and raw["geometry"]["n_axis_inv"] == 2400
          and raw["geometry"]["L_total_tol"] == 0.05, str(keys))
    check("緩和を _card 節と YAML のコメントに記録", raw["_card"]["L_total_tol_relaxed"]["to"] == 0.05
          and "内部許容差を緩めた" in text)
    check("カードの公差は変えない: 寸法はカードの長さ公差 1e-6 m で判定 (FAIL → 終了コード 2、YAML は書く)",
          [r["tol_m"] for r in logd["dimensions"]["rows"] if r["name"] == "divergent_length"] == [1.0e-6]
          and (rc == 2) == (not logd["dimensions"]["pass"]), f"rc {rc}, 寸法 {logd['dimensions']['pass']}")
    load_problem(out)

# ---------------------------------------------------------------------------------------------- 最終の雛形
print("== スロート曲率 (2026-10-13 ユーザ「スロート曲率も指定できるようにして」)")
import copy as _copy
def _th(**kw):
    c = card(); c["throat"] = dict(c["throat"], **kw); return c
p1 = W.card_to_problem(W.validate_card(_th(curvature_radius_m=0.075)), 0.025)
check("曲率半径 0.075 m (r_t 0.025) → geometry.R = 3 (r_t で割る)", abs(p1["geometry"]["R"] - 3.0) < 1e-12 and
      p1["_card"]["throat"]["curvature_source"] == "card throat.curvature_radius_m", p1["geometry"]["R"])
p2 = W.card_to_problem(W.validate_card(_th(curvature_radius_m=0.075)), 0.030)
check("r_t を変えると曲率半径の実寸を保って R を割り直す (0.075/0.030 = 2.5)", abs(p2["geometry"]["R"] - 2.5) < 1e-12, p2["geometry"]["R"])
p3 = W.card_to_problem(W.validate_card(_th(curvature_ratio=3.0)), 0.025)
check("curvature_ratio 3 → geometry.R = 3", p3["geometry"]["R"] == 3.0)
r, msg = raises(lambda: W.validate_card(_th(curvature_radius_m=0.075, curvature_ratio=3.0)), "どちらか 1 つ")
check("曲率を 2 つの書き方で与えると拒否", r, msg)
c4 = _th(curvature_ratio=3.0); c4["overrides"] = {"geometry.R": 2.5}
r, msg = raises(lambda: W.card_to_problem(W.validate_card(c4), 0.025), "両方")
check("throat の曲率と overrides.geometry.R の両方は拒否", r, msg)
r, msg = raises(lambda: W.card_to_problem(W.validate_card(_th(curvature_ratio=0.3)), 0.025), "Hall")
check("R < 0.5 (遷音速解の精度が落ちる範囲) は拒否", r, msg)

print("== 最終の雛形 axismach_final_v1 (M5 デモの r_t で)")
pf = W.card_to_problem(W.validate_card(card(recipe="axismach_final_v1")), rt_m5)
pf["mesh"].update(ni=1000, wall_first_frac=1e-5, wall_first_frac_throat=3e-6)
pf["mesh_euler"]["ni"] = 1000
with tempfile.TemporaryDirectory() as td:
    pth = Path(td) / "pf.yaml"
    W._write_problem(W.validate_card(card(recipe="axismach_final_v1")), pf, pth, ["test"])
    try:
        df = design_chain(load_problem(pth))
        okf, why = (not df["qa"]["violations"] and df["wall_repr"] == "joint"), f"L_c {df['L_c']:.6g}, x_e {df['wall'].x_e:.6g}"
    except ValueError as ex:
        okf, why = False, str(ex)[:200]
    # Euler の参照用の写し (physical_wall_repr・pw_upstream を除く): Euler の prepare が受け付け、設計壁が同じ (2026-10-13 case/70)
    eo = W.write_euler_companion(pth)
    okc, whyc = False, "写しが無い"
    if eo is not None and okf:
        import numpy as _np
        pe = load_problem(eo)
        de = design_chain(pe)
        xs = _np.linspace(df["wall"].x_in, df["wall"].x_e, 20001)
        dr = float(_np.abs(_np.asarray(df["wall"].r(xs)) - _np.asarray(de["wall"].r(xs))).max())
        okc = ("physical_wall_repr" not in pe.geometry and "pw_upstream" not in pe.geometry and dr == 0.0)
        whyc = f"{eo.name}: 設計壁の差 {dr:.3g} r_t"
check("最終の雛形の問題が load_problem と design_chain を通る (joint・単調拘束・analytic/converge)", okf, why)
check("Euler の参照用の写しが物理壁のキーを除き、設計壁がビット同一", okc, whyc)

# ---------------------------------------------------------------------------------------------- (c) case/45 の固定点 (--slow)
if SLOW:
    print("== (c) case/45 の生産の実寸の固定点 (--slow)")
    with tempfile.TemporaryDirectory() as td:
        cp = Path(td) / "c45.yaml"
        c = {"format": W.CARD_FORMAT, "name": "c45_prod_dims",
             "gas": {"species": {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.0220, "N2": 0.7245}},
             "conditions": {"Pt": 5.5e6, "Tt": 1600.0}, "M_exit": 6.0,
             "throat": {"diameter_physical_m": 2.0 * C45_PHYS_THROAT_R},
             "inlet": {"radius_m": 6.485 * C45_RT, "pipe_length_m": 0.5 * C45_RT, "contraction_length_m": 12.0 * C45_RT},
             "divergent_length_m": C45_WALL_END_M, "wall": {"thermal": "adiabatic"}, "recipe": "axismach_final_v1"}
        cp.write_text(yaml.safe_dump(c, sort_keys=False))
        out = Path(td) / "c45.problem.yaml"
        try:
            logd = W.build_problem(cp, out)
            hist = logd["r_t"]["history"]
            okc = 1 <= len(hist) <= 10 and hist[-1]["rel_change"] <= 1e-6
            why = "; ".join(f"{h['r_t_in_m']!r}→{h['r_t_out_m']!r} ({h['rel_change']:.2e}, {h['time_s']:.0f} s)" for h in hist)
        except Exception as ex:          # noqa: BLE001 (失敗の内容を出す)
            okc, why, logd = False, f"{type(ex).__name__}: {str(ex)[:300]}", None
        check("case/45 の実寸のカードで r_t の固定点が 10 回以内に止まる", okc, why)
        if logd is not None:
            rl = logd["L_total_tol_relaxed"]
            ok_yaml = out.is_file()
            text = out.read_text() if ok_yaml else ""
            check("case/45: 内部許容差の緩和が記録され YAML ができる (card_log・_card・コメント)",
                  ok_yaml and rl is not None and rl["to"] == 0.05 and yaml.safe_load(text)["_card"]["L_total_tol_relaxed"] == rl
                  and "内部許容差を緩めた" in text and load_problem(out) is not None, str(rl)[:240])
            print(f"   緩和: {rl}")
            print(f"   自動調整: {[(a['key'], a['from'], a['to']) for a in logd['auto_adjust']]}")
            print(f"   場: {logd['moc_field']}")
            print(f"   r_t {logd['r_t']['r_t_m']!r} m (生産 {C45_RT!r} m; 雛形は Hall 初期線・k_f 1 なので一致はしない)、"
                  f"寸法 {'PASS' if logd['dimensions']['pass'] else 'FAIL'}: "
                  f"{[(r['name'], round(r['diff_m'], 9)) for r in logd['dimensions']['rows']]}")
            print(f"   格子: ni {logd['mesh_estimate']['ni']}、wall_first_frac {logd['mesh_estimate']['wall_first_frac']}、"
                  f"_throat {logd['mesh_estimate']['wall_first_frac_throat']}、警告 {logd['warnings']}")
else:
    print("== (c) case/45 の固定点は --slow のときだけ (飛ばした)")

print(f"\n{'ALL OK' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
