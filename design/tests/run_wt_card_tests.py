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
  M_knot (plan tooling-nozzle-mknot-auto §5.1 #3、§6 V3): カードの既定は自動になったので、上の試験は**数値の M_knot のカード**
  (overrides の geometry.M_knot: 2.25 = 旧既定の 0.45 × M_exit; case/45 は 2.7 = 0.45 × 6) で回し、期待は変えない。
  M5 デモの生成は探索の関数 (search_range・select_mknot・resolution_gate・evaluate_mknot・make_eval_fn) を呼んだら失敗する差し替えの
  中で回す (数値の M_knot では探索しない)。M_knot の書き方 (書かない・auto・数値・不正) の検査と card_to_problem の M_knot=。
  自動の流れ (選ぶ・最終の組のゲート・やり直し・確かめ直し・記録) は探索の関数を合成の結果に差し替えて確かめる (設計スロート径:
  再利用・差 ≤ 0.01 の出力・やり直しが 2 回を超えて例外; 物理スロート径: 差 ≤ 0.01 で出力する組の r_t の固定点の 1 回を回し直す)。
  記録 (§5.1 #3b): card_log の M_knot_rule に全部 (掃引の表を含む)、YAML の _card.M_knot_rule は要約だけ (行 ≤ 4 kB、card_log の記録の
  sha256 で結ぶ; Euler の写しも同じ要約)、カードはゲートに半幅の列 (0.1, 0.2, 0.3) を渡し、窓を広げた記録が要約・コメントに出る。
  実際の探索での統合は §6 V4 (case/70 の c2b の条件、所要 数十分) で別に確かめる。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_wt_card_tests.py [--slow]
所要: 既定 約 10 分 (M5 デモの r_t の固定点 + design_chain 数回 + M_knot の自動の流れの合成の試験 約 7 分; 2026-10-14 の手元、別の計算と
並行)、--slow で case/45 の固定点と物理スロート径の M_knot のやり直しの合成の試験が加わる。
"""
import contextlib
import copy
import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate import mknot_select as MS  # noqa: E402
from forge_design.evaluate.runner_axismach import design_chain  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

FAIL = 0
# M_knot の探索の関数 (数値の M_knot のカードでは呼ばれない; §6 V3)
SEARCH_FUNCS = ("search_range", "select_mknot", "resolution_gate", "evaluate_mknot", "make_eval_fn")
FORBIDDEN_CALLS: list = []


def forbidden(name):
    def f(*a, **k):
        FORBIDDEN_CALLS.append(name)
        raise AssertionError(f"探索の関数 {name} が呼ばれた (数値の M_knot のカードでは呼ばない)")
    return f


@contextlib.contextmanager
def patched_attrs(obj, **kw):
    """obj の属性を差し替え、抜けたら戻す。"""
    old = {k: getattr(obj, k) for k in kw}
    for k, v in kw.items():
        setattr(obj, k, v)
    try:
        yield
    finally:
        for k, v in old.items():
            setattr(obj, k, v)
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
# 数値の M_knot (= 旧既定の 0.45 × M_exit)。M_knot の既定は自動になったので、既存の試験はこのカードで期待を変えずに回す (§6 V3)
MK_NUM = {"geometry.M_knot": 2.25}
BASE_NUM = dict(copy.deepcopy(BASE), overrides=dict(MK_NUM))


def card(**kw):
    """M5 デモのカード (数値の M_knot 2.25) の写し (kw で上書き、None は削除)。"""
    c = copy.deepcopy(BASE_NUM)
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
check("M_knot の既定は自動 (M5 デモのカードは geometry.M_knot を書いていない)", ok and v["M_knot_mode"] == "auto"
      and "書いていない" in v["M_knot_auto_source"], str(v.get("M_knot_mode")))
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
m5 = W.validate_card(BASE_NUM)
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
check("往復: Lc_mode from_length・M_knot = 0.45 M_exit (数値の M_knot 2.25)・p_ambient = 出口静圧",
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
print("== (a) 正例 1: M5 デモ (物理スロート径、数値の M_knot 2.25) の build_problem — 探索の関数を呼んだら失敗する差し替えの中で")
with tempfile.TemporaryDirectory() as td:
    out = Path(td) / "problem_m5_air_demo.yaml"
    m5_num = Path(td) / "m5_air_demo_num.yaml"
    m5_num.write_text(yaml.safe_dump(BASE_NUM, sort_keys=False, allow_unicode=True))
    with patched_attrs(MS, **{k: forbidden(k) for k in SEARCH_FUNCS}):
        logd = W.build_problem(m5_num, out)
    check("数値の M_knot のカード: 探索の関数を 1 度も呼ばない (呼んだら失敗する差し替えで通った)・記録 M_knot_rule を書かない",
          "M_knot_rule" not in logd and "M_knot_rule" not in yaml.safe_load(out.read_text())["_card"] and FORBIDDEN_CALLS == [],
          str(FORBIDDEN_CALLS))
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
                                      overrides={"geometry.x_end_margin": 1.0, **MK_NUM}), sort_keys=False))
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

# ---------------------------------------------------------------------------------------------- M_knot の書き方 (軽い)
print("== M_knot の書き方 (plan tooling-nozzle-mknot-auto §4.4)")
va = W.validate_card(card(overrides={"geometry.M_knot": "auto", "geometry.R": 2.0}))
check("overrides の geometry.M_knot: auto → 自動 (overrides から外して出典を残す; 他の上書きはそのまま)",
      va["M_knot_mode"] == "auto" and va["overrides"] == {"geometry.R": 2.0} and "auto" in va["M_knot_auto_source"], str(va["overrides"]))
va2 = W.validate_card(card(overrides={"geometry": {"M_knot": "auto"}}))
check("入れ子の書き方 overrides: {geometry: {M_knot: auto}} も自動", va2["M_knot_mode"] == "auto" and va2["overrides"] == {})
vn = W.validate_card(card())
check("数値の geometry.M_knot → 固定 (overrides に残す; 探索しない)", vn["M_knot_mode"] == "fixed" and vn["overrides"] == MK_NUM)
for nm, bad in (("文字列", "fast"), ("真偽値", True)):
    r, msg = raises(lambda bad=bad: W.validate_card(card(overrides={"geometry.M_knot": bad})), "geometry.M_knot", "auto")
    check(f"負例: geometry.M_knot に数値でも auto でもない値 ({nm})", r, msg[:160])
mauto = W.validate_card(BASE)
r, msg = raises(lambda: W.card_to_problem(mauto, 0.025), "自動", "M_knot=")
check("自動のカードの card_to_problem は M_knot= が無ければ例外 (黙って 0.45 × M_exit に戻さない)", r, msg[:160])
pa = W.card_to_problem(mauto, 0.025, M_knot=2.069270614)
check("自動のカードの card_to_problem(M_knot=x) → geometry.M_knot = x (overrides の記録には入らない)",
      pa["geometry"]["M_knot"] == 2.069270614 and not any(o["key"] == "geometry.M_knot" for o in pa["_card"]["overrides"]))
r, msg = raises(lambda: W.card_to_problem(vn, 0.025, M_knot=2.1), "M_knot=")
check("数値の M_knot のカードに M_knot= を渡すと例外", r, msg[:160])
pq = W.card_to_problem(W.validate_card(dict(copy.deepcopy(BASE), overrides={"geometry.axis_law": "quintic"})), 0.025)
check("knot でない軸則 (quintic) の自動のカード: geometry.M_knot を書かない (M_knot= なしで通る)", "M_knot" not in pq["geometry"])

# ---------------------------------------------------------------------------------------------- M_knot の自動の流れ (合成)
print("== M_knot の自動の流れ (探索の関数を合成の結果に差し替え; 実際の探索での統合は §6 V4)")
REAL_EVAL = MS.evaluate_mknot
M_A_FAKE = 1.2
SETTINGS = {"coarse": 0.1, "fine": 0.01, "tie_deg": 1e-3, "dM_safe": 0.05, "bisect_tol": 1e-3, "sin_min": 0.02, "margin_min_deg": 1.0,
            "grid_origin": M_A_FAKE + 0.1}


def fake_row(M, role="fake"):
    return {"M_K": float(M), "ok": True, "reasons": [], "role": role, "theta_max_deg": 15.0 + (float(M) - 2.1) ** 2, "margin_deg": 4.0,
            "min_sin_angle": 0.05, "mpp_sign_changes": 1, "x_K": 10.0, "L_c": 40.0, "moc_gate_ok": True, "moc_gate_pass": True,
            "n_axis": 2400, "M_A": M_A_FAKE}


class FakeSearch:
    """select_mknot・resolution_gate・search_range・make_eval_fn を合成の結果に差し替える。moves: ゲートの呼び出しごとの
    chosen_N − 入力 (足りなければ最後の値)。widen: ゲートが窓を広げた回数 (合成; 記録と YAML のコメントの確かめ用)。出力する組の判定
    (最終の設計チェーンの θ_max との照合) のため、chosen_N の N の行だけは本物の evaluate_mknot で評価する。"""

    def __init__(self, c0, moves, widen=0):
        self.c0, self.moves, self.widen, self.gates, self.n_select = c0, list(moves), int(widen), [], 0

    def search_range(self, base, base_dir=None):
        return {"M_A": M_A_FAKE, "M_d": 5.0, "x_A": 0.5, "Mp_A": 0.6, "Mpp_A": 0.03}

    def select_mknot(self, eval_fn, M_A, M_d, **kw):
        self.n_select += 1
        return {"rule": MS.RULE_NAME, "M_A": M_A, "M_d": M_d, "settings": dict(SETTINGS), "chosen": self.c0, "chosen_row": fake_row(self.c0),
                "distance_to_boundary": 0.1, "reason": "合成", "constrained": False, "theta_min": {}, "theta_min_admissible": {},
                "tie_set": [self.c0], "excluded_by_distance": [], "safe_points": [], "local_minima": [], "windows": [], "boundaries": [],
                "warnings": [], "table": [fake_row(self.c0 + 0.01 * j, "fine") for j in range(-3, 4)], "n_eval": 7, "elapsed_s": 0.0}

    def resolution_gate(self, base, chosen, N, select_result=None, seed_2N=None, half_widths=(0.1,), workers=(1, 1), log=None, **kw):
        mv = self.moves[min(len(self.gates), len(self.moves) - 1)]
        cN = round(float(chosen) + mv, 9)
        self.gates.append({"input": float(chosen), "chosen_N": cN, "r_t": float(base["spec"]["r_throat"]), "seed_2N": seed_2N is not None,
                           "n_seed_N": len((select_result or {}).get("table") or []), "half_widths": tuple(half_widths)})
        rowN = dict(REAL_EVAL(base, cN, N), role="window")
        row2 = fake_row(cN, "window")
        hws = list(half_widths)[:self.widen + 1]
        widen = [{"half_width": h, "window": [chosen - h, chosen + h], "chosen_N": cN, "chosen_2N": cN,
                  "edge_values": ({"N": cN} if i < self.widen else {}), "min_offset_N": 0, "min_offset_2N": 0, "action": "合成"}
                 for i, h in enumerate(hws)]
        return {"N": N, "N2": 2 * N, "chosen": float(chosen), "window": [chosen - hws[-1], chosen + hws[-1]], "settings": {},
                "half_width": hws[-1], "n_widen": len(hws) - 1, "widen": widen, "edge_values": {},
                "pass": True, "checks": {"change": True, "chosen_inside_N": True}, "reasons": [], "warnings": [], "chosen_N": cN,
                "chosen_N_reason": "合成",
                "chosen_2N": cN, "chosen_2N_reason": "合成", "moved_N": cN - float(chosen), "distance_to_conservative_boundary": 0.1,
                "distance_2N": 0.1, "distance_input": 0.1, "min_offset_N": 0, "min_offset_2N": 0,
                "common_boundaries": [{"kind": "lower", "feasible_end": 1.97, "lo": 1.969, "hi": 1.97, "source": "N", "nearest_other": None,
                                       "reasons_infeasible": []}], "boundaries_N": [], "boundaries_2N": [], "rule_N": {}, "rule_2N": {},
                "rows_N": [fake_row(chosen + 0.01 * j, "window") for j in range(-10, 11)],
                "rows_2N": [fake_row(chosen + 0.01 * j, "window") for j in range(-10, 11)],
                "row_chosen_N": rowN, "row_chosen_2N": row2, "row_input_N": fake_row(chosen), "row_input_2N": fake_row(chosen),
                "rows_selected": {f"{float(chosen):.12g}": {"at_N": fake_row(chosen), "at_2N": fake_row(chosen)},
                                  f"{cN:.12g}": {"at_N": rowN, "at_2N": row2}}, "elapsed_s": 0.0}

    def make_eval_fn(self, base, n=None, base_dir=None):
        return lambda M: fake_row(M, "context")


def run_fake(fake, card_dict, name, **bp):
    """合成の探索で build_problem を回す。戻り (card_log, YAML の文字列, 例外の文字列)。"""
    with tempfile.TemporaryDirectory() as td:
        cp = Path(td) / f"{name}.yaml"
        cp.write_text(yaml.safe_dump(card_dict, sort_keys=False, allow_unicode=True))
        out = Path(td) / f"problem_{name}.yaml"
        wk1 = {"N": 1, "2N": 1, "available_gb": None, "per_eval_gb_N": 0.0, "per_eval_gb_2N": 0.0}
        with patched_attrs(MS, search_range=fake.search_range, select_mknot=fake.select_mknot, resolution_gate=fake.resolution_gate,
                           make_eval_fn=fake.make_eval_fn), patched_attrs(W, mknot_workers=lambda N: dict(wk1)):
            try:
                logd = W.build_problem(cp, out, log=lambda *a, **k: None, **bp)
            except Exception as ex:          # noqa: BLE001 (例外の試験)
                return None, None, f"{type(ex).__name__}: {ex}"
            return logd, out.read_text(), None


AUTO_DESIGN = card(name="m5_auto_design", throat={"diameter_design_m": 0.050}, overrides={})
t0 = time.time()
fk = FakeSearch(2.1, [0.0], widen=1)
logd, text, err = run_fake(fk, AUTO_DESIGN, "auto_a")
raw = yaml.safe_load(text) if text else {}
mkr = (raw.get("_card") or {}).get("M_knot_rule") or {}
check("自動・設計スロート径・ゲートが動かない: 最初の r_t のゲートを最終の組に再利用し、YAML の geometry.M_knot = 記録の chosen = chosen_N",
      err is None and raw["geometry"]["M_knot"] == 2.1 == mkr["chosen"] and mkr["gate_first_reused"] and fk.n_select == 1
      and len(fk.gates) == 1 and [e["stage"].split(" ")[0] for e in mkr["events"]] == ["最初の", "最終の組"], err or str(fk.gates))
check("記録: 規則の名前・出典 (書いていない)・雛形の係数は記録だけ (0.45 × 5 = 2.25 は使わない)・最終の設計チェーンの θ_max = ゲートの行",
      err is None and mkr["rule"] == "min_theta_max_v1" and "書いていない" in mkr["source"] and mkr["recipe_M_knot_factor"] == 0.45
      and mkr["recipe_M_knot_value_unused"] == 2.25 and mkr["theta_max_deg_final_chain"] == mkr["theta_max_deg"]
      and mkr["ok_at_output"] == {"N": True, "2N": True} and logd["M_knot_rule"]["constraints_at_output"]["N"]["ok"]
      and logd["M_knot_rule"]["constraints_at_output"]["2N"]["ok"], str({k: mkr.get(k) for k in ("rule", "source")}))
CL_A = "problem_auto_a.yaml.card_log.json"
check("記録: YAML の _card.M_knot_rule は card_log の M_knot_rule から作った要約 (card_log の記録の sha256 で結ぶ)・対応表の行・"
      "YAML のコメント (geometry.M_knot と見出し)",
      err is None and W._mknot_yaml_summary(logd["M_knot_rule"], CL_A) == mkr and mkr["card_log"]["file"] == CL_A
      and mkr["card_log"]["sha256"] == W._mknot_record_sha256(logd["M_knot_rule"])
      and any(r["card_key"] == "M_knot (自動)" for r in logd["mapping"])
      and "# M_knot: 自動 (規則 min_theta_max_v1" in text and "  M_knot: 2.1   # 自動 (規則 min_theta_max_v1)" in text
      and f"掃引の表と全部の記録: {CL_A} の M_knot_rule (その記録の sha256 {mkr['card_log']['sha256']}" in text
      and "M_knot は自動" in json.dumps(logd["mapping"], ensure_ascii=False), err or "")
ln_a = [ln for ln in (text or "").splitlines() if ln.startswith("  M_knot_rule:")]
# Euler の写し: この雛形 (explore) では書かれないので、物理壁のキーの行を 1 行足した YAML から write_euler_companion で作って比べる
# (写しは行を落とすだけ; 最終の雛形の実際の写しは生成器が書いた後に mknot_card_log_check で照合する)
eu_a, chk_e = [], {"ok": None}
if text:
    with tempfile.TemporaryDirectory() as td:
        yp = Path(td) / "problem_auto_a.yaml"
        lines = text.splitlines(keepends=True)
        k = next(i for i, ln in enumerate(lines) if ln.startswith("geometry:"))
        yp.write_text("".join(lines[:k + 1] + ["  physical_wall_repr: single_bspline\n"] + lines[k + 1:]))
        (Path(td) / CL_A).write_text(json.dumps(logd, indent=1, ensure_ascii=False))
        eo = W.write_euler_companion(yp)
        eu_a = [ln for ln in eo.read_text().splitlines() if ln.startswith("  M_knot_rule:")]
        chk_e = W.mknot_card_log_check(eo)
check("YAML の大きさ (§4.4 要約だけ): _card.M_knot_rule の行 ≤ 4 kB・掃引の表を持たない、Euler の写しも同じ要約の行で card_log と照合できる",
      err is None and len(ln_a) == 1 and len(ln_a[0].encode()) <= 4096 and "table" not in mkr and eu_a == ln_a and chk_e["ok"] is True,
      f"{len(ln_a[0].encode()) if ln_a else None} byte、YAML の最長の行 {max(len(ln.encode()) for ln in (text or ' ').splitlines())} byte、"
      f"Euler の写しの照合 {chk_e}")
check("ゲート: カードは半幅の列 (0.1, 0.2, 0.3) を渡し、窓を広げた記録 (回数・最後の半幅・端だった値) が YAML の要約・コメントに出る",
      err is None and fk.gates[0]["half_widths"] == (0.1, 0.2, 0.3) and mkr["gate"]["n_widen"] == 1 and mkr["gate"]["half_width"] == 0.2
      and mkr["gate"]["widen"][0]["edge_values"] == {"N": 2.1} and "窓を 1 回広げた" in text
      and [e["gate_n_widen"] for e in mkr["events"]] == [1, 1], err or str(mkr.get("gate")))
check("記録: 掃引の表 (M_K・成立・理由・θ_max・余裕・最小 sin 角・M″ の符号反転・x_K・MOC のゲート) は card_log だけ",
      err is None and all(all(k in r for k in ("M_K", "ok", "reasons", "theta_max_deg", "margin_deg", "min_sin_angle", "mpp_sign_changes",
                                                "x_K", "moc_gate_ok")) for r in logd["M_knot_rule"]["table"]["select"]), "")
fk = FakeSearch(2.1, [0.02, 0.005])
logd, text, err = run_fake(fk, AUTO_DESIGN, "auto_c")
raw = yaml.safe_load(text) if text else {}
mkr = (raw.get("_card") or {}).get("M_knot_rule") or {}
check("自動・設計スロート径・最終の組のゲートで差 0.005 ≤ 0.01: chosen_N 2.125 を出力 (同じ組の評価を再利用して詰め直す)",
      err is None and raw["geometry"]["M_knot"] == 2.125 == mkr["chosen"] and [g["input"] for g in fk.gates] == [2.1, 2.12]
      and fk.gates[1]["seed_2N"] and fk.gates[1]["n_seed_N"] > 0 and mkr["redo"] == 0 and not mkr["gate_first_reused"],
      err or str(fk.gates))
fk = FakeSearch(2.1, [0.02, 0.02, 0.0])
logd, text, err = run_fake(fk, AUTO_DESIGN, "auto_redo")
raw = yaml.safe_load(text) if text else {}
mkr = (raw.get("_card") or {}).get("M_knot_rule") or {}
check("自動・設計スロート径・差 0.02 > 0.01 → やり直して (1 回) 収まった値 2.14 を出力",
      err is None and raw["geometry"]["M_knot"] == 2.14 and mkr["redo"] == 1 and [g["input"] for g in fk.gates] == [2.1, 2.12, 2.14]
      and "やり直す" in logd["M_knot_rule"]["events"][1]["action"], err or str(fk.gates))
fk = FakeSearch(2.1, [0.02])
logd, text, err = run_fake(fk, AUTO_DESIGN, "auto_never")
check("自動・差が毎回 0.02 → やり直しが 2 回を超えたら例外 (YAML を書かない)",
      err is not None and err.startswith("RuntimeError") and "2 回やり直しても" in err and len(fk.gates) == 4, (err or "")[:200])
fk = FakeSearch(2.1, [0.0, 0.004, 0.0])
logd, text, err = run_fake(fk, card(name="m5_auto_phys", overrides={}), "auto_phys")
raw = yaml.safe_load(text) if text else {}
mkr = (raw.get("_card") or {}).get("M_knot_rule") or {}
hist = raw.get("_card", {}).get("r_t_history") or []
ver = [h for h in hist if str(h.get("stage", "")).startswith("M_knot の確かめ直し")]
check("自動・物理スロート径・最終の組の差 0.004 ≤ 0.01: 出力する組 (chosen_N 2.104、最後の r_t) で r_t の固定点の 1 回を回し直し、"
      "寸法ゲート PASS",
      err is None and raw["geometry"]["M_knot"] == 2.104 == mkr["chosen"] and fk.gates[0]["r_t"] == 0.025 and fk.gates[1]["r_t"] != 0.025
      and len(ver) >= 1 and ver[0]["M_K"] == 2.104 and ver[0]["r_t_in_m"] == fk.gates[1]["r_t"] and logd["dimensions"]["pass"]
      and all(h.get("M_K") in (2.1, 2.104) for h in hist),
      err or f"ゲート {fk.gates}、確かめ直し {[(h['r_t_in_m'], h['r_t_out_m'], h['rel_change']) for h in ver]}、events "
             f"{[e['action'][:60] for e in ((logd or {}).get('M_knot_rule') or {}).get('events', [])]}")
if err is None:
    print(f"   物理スロート径の確かめ直し: {[(h['r_t_in_m'], h['r_t_out_m'], h['rel_change']) for h in ver]}、出力 r_t {raw['spec']['r_throat']!r}、"
          f"ゲートの r_t {[g['r_t'] for g in fk.gates]}")
if SLOW:
    fk = FakeSearch(2.1, [0.0, 0.03, 0.0])
    logd, text, err = run_fake(fk, card(name="m5_auto_phys_redo", overrides={}), "auto_phys_redo")
    raw = yaml.safe_load(text) if text else {}
    mkr = (raw.get("_card") or {}).get("M_knot_rule") or {}
    hist = raw.get("_card", {}).get("r_t_history") or []
    check("自動・物理スロート径・最終の組の差 0.03 > 0.01 (--slow): chosen_N 2.13 で r_t の固定点を前の r_t からやり直し、その組でゲートを"
          "回し直して出力",
          err is None and raw["geometry"]["M_knot"] == 2.13 and mkr["redo"] == 1 and [g["input"] for g in fk.gates] == [2.1, 2.1, 2.13]
          and any(h.get("M_K") == 2.13 for h in hist) and fk.gates[2]["r_t"] == raw["spec"]["r_throat"] and logd["dimensions"]["pass"],
          err or str(fk.gates))
print(f"   M_knot の自動の流れ (合成) {time.time() - t0:.0f} s")

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
             "divergent_length_m": C45_WALL_END_M, "wall": {"thermal": "adiabatic"}, "recipe": "axismach_final_v1",
             "overrides": {"geometry.M_knot": 2.7}}          # 数値の M_knot (= 旧既定の 0.45 × 6)
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
