#!/usr/bin/env python3
"""要求カードの壁温分布 (plan tooling-nozzle-wall-temperature-profile §5.1 #3 のカード・格子の見積もり) の試験 — §6 V1 (カード側)・V1c。

  V1  カードの `wall`: `{thermal: isothermal, Tw}`・`{.., Tw_profile: {s|x_m: [...], Tw_K: [...]}}`・`{.., Tw_profile_csv: <path>}` の厳密に 1 つ。
      正例 (s・x_m の表、CSV の s・x_m)。CSV の相対パスは**カードのファイルの所在**が基準 (cwd を別にし、cwd に同名の別の CSV を置いて
      読まれないことを確かめる)、CSV の絶対パスと sha256 が card_log・`_card.card` に残る。生成した YAML の `spec.wall_thermal.Tw_profile`
      が元の表と一致し、`probdef.Problem` を YAML から作って `wall_thermal["Tw_profile"]` が元の表と一致する。
      負例 (`Tw` との併用・CSV と表の併用・どれも無い・未知のキー・非単調・負・NaN・長さ違い・s が [−2, 1] を覆わない・x_m が 0 を含まない・
      CSV の列名違い・CSV が無い・CSV の列数違い・数値でない値・dict のカードで相対パス)。表の検査は `wall_temperature.parse_profile`。
  V1c 格子の見積もり (codex plan レビュー M4): スロートの温度が同じ 500 K で上流・下流が違う 2 つの表 (A: 上流 320/下流 330、
      B: 上流 450/下流 480) を s と x_m の両方で与え、拡大部の中ほど 1 つと縮流部 1 つの断面の局所の T_w と丸め前の第 1 層厚 y1 が、
      試験の中で独立に書き下した式 (区切りの点の写し・線形補間・Eckert・平板 Cf・Sutherland) の期待値と相対 1e-8。s の区切りの点は
      **その設計の壁** (geometry.L_pipe·r_t・L_U·r_t・設計壁の終点 x_e·r_t; plan §4.3 [2026-10-13 訂正]) から試験側で取る。
      誤実装の検出力: (a) 代表温度 (s = 0 の 500 K) で置き換える、(b) 無次元の x [r_t] をそのまま表の x_m に当てる (np.interp の端値)、の
      y1 が期待値から 1e-3 以上ずれることを数値で示す。(c) 区切りの点をカードの長さで取る誤実装は、設計壁の終点がカードの拡大部長と
      ずれると、値のずれ (V1c の相対 1e-8 で落ちる) か、終点が長いときは表が壁を覆わない例外で検出できる。
  ほか: 一様のカードの見積もりは T_w 一定で `wall_temperature` を返さない。設計壁の終点がカードの拡大部長と ±1 mm ずれても (超えても)
      s の表は壁を覆い、区切りの点 = 設計壁の値 (相対 1e-12)・局所の T_w はカードの長さに依らない。x_m の表が壁を覆わなければ例外。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_wt_card_profile_tests.py
所要: 約 30 秒 (design_chain 2 回)。
"""
import copy
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate.ic import _invert_area_ratio_gas, _isentropic_pressure_ratio_gas  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

FAIL = 0
M5_CARD = DESIGN / "forge_design" / "intake" / "examples" / "m5_air_demo.yaml"
BASE = yaml.safe_load(M5_CARD.read_text())
RT = 0.025                                   # 設計スロート径 0.050 m のカード (r_t = 径/2、反復なし)
# デモの表 (case/70 card_m5_air_demo_walltemp.yaml と同じ)
S_DEMO = {"s": [-2.0, -1.0, -0.3, 0.0, 0.2, 0.5, 1.0], "Tw_K": [320.0, 330.0, 420.0, 500.0, 450.0, 380.0, 330.0]}
X_DEMO = {"x_m": [-0.5, -0.3, -0.09, 0.0, 0.24, 0.6, 1.25], "Tw_K": [320.0, 330.0, 420.0, 500.0, 450.0, 380.0, 330.0]}


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


def card(**kw):
    """M5 の例のカードを設計スロート径 0.050 m にした写し (kw で上書き、None は削除)。"""
    c = copy.deepcopy(BASE)
    c["throat"] = {"diameter_design_m": 0.050}
    for k, v in kw.items():
        if v is None:
            c.pop(k, None)
        else:
            c[k] = v
    return c


def iso(**kw):
    return {"thermal": "isothermal", **kw}


def write_yaml_problem(c: dict, prob: dict, path: Path) -> None:
    """card_to_problem の dict に格子の値を仮に入れて YAML に書く (load_problem を通すため)。"""
    prob = copy.deepcopy(prob)
    prob["mesh"].update(ni=1000, wall_first_frac=1e-5, wall_first_frac_throat=3e-6)
    prob["mesh_euler"]["ni"] = 1000
    W._write_problem(c, prob, path, ["test"])


def same_table(prof: dict, tbl: dict) -> bool:
    """probdef が返す {coord, coords, Tw_K} が元の表 {座標: [...], Tw_K: [...]} と一致するか。"""
    coord = next(k for k in ("s", "x_m") if k in tbl)
    return prof["coord"] == coord and list(prof["coords"]) == [float(v) for v in tbl[coord]] \
        and list(prof["Tw_K"]) == [float(v) for v in tbl["Tw_K"]]


def problem_roundtrip(name: str, c: dict, tbl: dict) -> None:
    """card_to_problem → YAML → load_problem で spec.wall_thermal.Tw_profile と Problem.wall_thermal が元の表と一致するか。"""
    prob = W.card_to_problem(c, RT)
    coord = next(k for k in ("s", "x_m") if k in tbl)
    want = {"mode": "isothermal", "Tw_profile": {coord: [float(v) for v in tbl[coord]], "Tw_K": [float(v) for v in tbl["Tw_K"]]}}
    with tempfile.TemporaryDirectory() as td:
        pth = Path(td) / "p.yaml"
        write_yaml_problem(c, prob, pth)
        raw = yaml.safe_load(pth.read_text())
        wt = load_problem(pth).wall_thermal
    rep = float(np.interp(0.0, want["Tw_profile"][coord], want["Tw_profile"]["Tw_K"]))
    check(f"{name}: YAML の spec.wall_thermal = {{mode: isothermal, Tw_profile: 元の表}}", prob["spec"]["wall_thermal"] == want
          and raw["spec"]["wall_thermal"] == want, str(raw["spec"]["wall_thermal"])[:160])
    check(f"{name}: Problem を YAML から作って wall_thermal['Tw_profile'] が元の表と一致 (Tw は代表値 = 座標 0 の値)",
          wt["mode"] == "isothermal" and same_table(wt["Tw_profile"], tbl) and wt["Tw"] == rep, f"{wt['Tw_profile']}, Tw {wt['Tw']}")


# ---------------------------------------------------------------------------------------------- V1 正例
print("== V1 正例: 表 (s・x_m) と CSV (s・x_m)")
v = W.validate_card(card(wall=iso(Tw_profile=S_DEMO)))
check("s の表: 検査を通り {thermal, Tw_profile: 表} に揃う", v["wall"] == {"thermal": "isothermal", "Tw_profile": S_DEMO}, str(v["wall"]))
problem_roundtrip("s の表", v, S_DEMO)
# x_m の表を YAML の文字で書く (整数と、PyYAML が文字列と読む YAML 1.2 の浮動小数 3.2e2 (指数に符号なし) を含む)
xtxt = "x_m: [-0.5, -0.3, -0.09, 0, 0.24, 0.6, 1.25]\nTw_K: [3.2e2, 330, 420, 500, 450, 380, 330]\n"
xraw = yaml.safe_load(xtxt)
v = W.validate_card(card(wall=iso(Tw_profile=xraw)))
check("x_m の表 (整数と 3.2e2 形の数値も受ける): float の表に揃う", isinstance(xraw["Tw_K"][0], str) and v["wall"]["Tw_profile"] == X_DEMO
      and all(isinstance(t, float) for t in v["wall"]["Tw_profile"]["Tw_K"]), str(v["wall"]["Tw_profile"]))
problem_roundtrip("x_m の表", v, X_DEMO)

orig_cwd = os.getcwd()
with tempfile.TemporaryDirectory() as tda, tempfile.TemporaryDirectory() as tdb:
    A, B = Path(tda), Path(tdb)
    (A / "sub").mkdir()
    csv_s = "s,Tw_K\n" + "".join(f"{s!r},{t!r}\n" for s, t in zip(S_DEMO["s"], S_DEMO["Tw_K"]))
    (A / "tw_s.csv").write_text(csv_s)
    # x_m の CSV は BOM・CRLF・空白・空行つき (表計算ソフトの書き出しを想定)
    csv_x = "﻿ x_m , Tw_K \r\n" + "".join(f" {x!r} , {t!r}\r\n" for x, t in zip(X_DEMO["x_m"], X_DEMO["Tw_K"])) + "\r\n"
    (A / "sub" / "tw_x.csv").write_bytes(csv_x.encode("utf-8"))
    # cwd 側に同名の別の CSV (誤ってこちらを読んだら表が違う)
    (B / "tw_s.csv").write_text("s,Tw_K\n-2,999\n1,999\n")
    (B / "sub").mkdir()
    (B / "sub" / "tw_x.csv").write_text("x_m,Tw_K\n-1,999\n2,999\n")
    cards = {"CSV の s (相対パス tw_s.csv)": ("tw_s.csv", A / "tw_s.csv", S_DEMO),
             "CSV の x_m (相対パス sub/tw_x.csv、BOM・CRLF・空白・空行)": ("sub/tw_x.csv", A / "sub" / "tw_x.csv", X_DEMO),
             "CSV の s (絶対パス)": (str(A / "tw_s.csv"), A / "tw_s.csv", S_DEMO)}
    os.chdir(B)
    try:
        for i, (name, (given, real, tbl)) in enumerate(cards.items()):
            cp = A / f"card_{i}.yaml"
            cp.write_text(yaml.safe_dump(card(name=f"csv_{i}", wall=iso(Tw_profile_csv=given)), sort_keys=False))
            v = W.load_card(cp)
            src = v["wall"].get("Tw_profile_csv") or {}
            check(f"{name}: カードの所在が基準で解決し (cwd の同名の CSV は読まない)、表が CSV と一致",
                  v["wall"]["Tw_profile"] == tbl and src.get("abs_path") == str(real.resolve()), f"{src.get('abs_path')}")
            check(f"{name}: 記録に絶対パスと sha256 (ファイルの sha256 と一致)",
                  src == {"path": given, "abs_path": str(real.resolve()), "sha256": hashlib.sha256(real.read_bytes()).hexdigest()})
            problem_roundtrip(name, v, tbl)
        # dict のカードで相対パス (基準が無い) は拒否、絶対パスは通る
        r, msg = raises(lambda: W.validate_card(card(wall=iso(Tw_profile_csv="tw_s.csv"))), "カードのファイルの所在が基準")
        check("負例: dict のカード (所在なし) で CSV の相対パス", r, msg[:160])
        v = W.validate_card(card(wall=iso(Tw_profile_csv=str(A / "tw_s.csv"))))
        check("dict のカードでも CSV の絶対パスは読める", v["wall"]["Tw_profile"] == S_DEMO)

        # ------------------------------------------------------------------------------ build_problem (CLI、cwd は別)
        print("== V1 build_problem (CLI、設計スロート径、CSV の s、cwd はカードと別)")
        cp = A / "card_build.yaml"
        cp.write_text(yaml.safe_dump(card(name="csv_build", wall=iso(Tw_profile_csv="tw_s.csv")), sort_keys=False))
        out = A / "problem_build.yaml"
        rc = W.main([str(cp), "--out", str(out)])
        logd = json.loads(Path(str(out) + ".card_log.json").read_text())
        text = out.read_text()
        raw = yaml.safe_load(text)
        sha = hashlib.sha256((A / "tw_s.csv").read_bytes()).hexdigest()
        src_rec = {"path": "tw_s.csv", "abs_path": str((A / "tw_s.csv").resolve()), "sha256": sha}
        check("CLI の終了コード 0", rc == 0, str(rc))
        check("YAML の spec.wall_thermal は CSV の表を埋め込む (YAML 単体で再現できる)",
              raw["spec"]["wall_thermal"] == {"mode": "isothermal", "Tw_profile": S_DEMO})
        wt = load_problem(out).wall_thermal
        check("Problem.wall_thermal['Tw_profile'] が CSV の表と一致", same_table(wt["Tw_profile"], S_DEMO) and wt["Tw"] == 500.0)
        check("card_log の card.wall と YAML の _card.card.wall に CSV の絶対パスと sha256",
              logd["card"]["wall"]["Tw_profile_csv"] == src_rec and raw["_card"]["card"]["wall"]["Tw_profile_csv"] == src_rec,
              str(logd["card"]["wall"]["Tw_profile_csv"]))
        rows = {r["card_key"]: r for r in logd["mapping"]}
        rr = rows.get("wall.Tw_profile_csv") or {}
        check("対応表: wall.thermal → spec.wall_thermal.mode、wall.Tw_profile_csv → spec.wall_thermal.Tw_profile",
              rows.get("wall.thermal", {}).get("yaml_key") == "spec.wall_thermal.mode" and rr.get("yaml_key") == "spec.wall_thermal.Tw_profile"
              and rr.get("yaml_value") == S_DEMO and rr.get("card_value") == src_rec and "wall" not in rows, str(list(rows))[:200])
        me = logd["mesh_estimate"]
        twr = me.get("wall_temperature") or {}
        gb, rtb, xeb = raw["geometry"], float(raw["spec"]["r_throat"]), float(logd["moc_field"]["x_e"])
        want_anc = {"x_in_m": -(gb["L_pipe"] * rtb + gb["L_U"] * rtb), "x_contraction_m": -(gb["L_U"] * rtb), "x_throat_m": 0.0,
                    "x_exit_m": xeb * rtb}
        got_anc = twr.get("anchors_m") or {}
        check("見積もり: 区切りの点は設計壁から (L_pipe·r_t・L_U·r_t・設計壁の終点 x_e·r_t、相対 1e-12; 出典 design wall)",
              set(got_anc) == set(want_anc) and all(abs(got_anc[k] - v) <= 1e-12 * max(1.0, abs(v)) for k, v in want_anc.items())
              and str(twr.get("anchors_source", "")).startswith("design wall")
              and raw["_card"]["mesh_estimate"]["wall_temperature"]["anchors_m"] == got_anc, str(got_anc))
        check("見積もり: tw_note「等温壁 分布 (s の表 7 点、T_w 320〜500 K)」・card_log に断面ごとの配列を書かない",
              me["inputs"]["wall"] == "等温壁 分布 (s の表 7 点、T_w 320〜500 K)" and "sections" not in me, me["inputs"]["wall"])
        check("YAML のヘッダに CSV の絶対パス・sha256 と局所の T_w(x) の見積もり", src_rec["abs_path"] in text and sha in text
              and "局所の T_w(x)" in text)
        print(f"   card_log mesh_estimate.wall_temperature: {json.dumps({k: twr.get(k) for k in ('anchors_m', 'range_m', 'wall_end_minus_card_exit_m', 'Tw_K_at_wall_first_frac', 'Tw_K_at_wall_first_frac_throat')}, ensure_ascii=False)}")
    finally:
        os.chdir(orig_cwd)

# ---------------------------------------------------------------------------------------------- V1 負例
print("== V1 負例")
with tempfile.TemporaryDirectory() as td:
    T = Path(td)
    (T / "ok.csv").write_text("s,Tw_K\n-2,320\n0,500\n1,330\n")
    (T / "hdr_x.csv").write_text("x,Tw_K\n-0.5,320\n0,500\n1.3,330\n")
    (T / "hdr_semicolon.csv").write_text("s;Tw_K\n-2;320\n0;500\n1;330\n")
    (T / "hdr_tw.csv").write_text("s,Tw\n-2,320\n0,500\n1,330\n")
    (T / "nan.csv").write_text("s,Tw_K\n-2,320\n0,nan\n1,330\n")
    (T / "cols3.csv").write_text("s,Tw_K\n-2,320,1\n0,500\n1,330\n")
    (T / "text.csv").write_text("s,Tw_K\n-2,320\n0,hot\n1,330\n")
    (T / "short.csv").write_text("s,Tw_K\n-1.5,320\n0,500\n1,330\n")
    (T / "empty.csv").write_text("\n\n")

    def ld(wall):
        cp = T / "card.yaml"
        cp.write_text(yaml.safe_dump(card(wall=wall), sort_keys=False))
        return W.load_card(cp)

    neg = [
        ("Tw と Tw_profile の併用", lambda: ld(iso(Tw=350.0, Tw_profile=S_DEMO)), ("どれか 1 つ", "併用", "Tw_profile")),
        ("CSV と表の併用", lambda: ld(iso(Tw_profile=S_DEMO, Tw_profile_csv="ok.csv")), ("併用", "Tw_profile_csv")),
        ("Tw と CSV の併用", lambda: ld(iso(Tw=350.0, Tw_profile_csv="ok.csv")), ("併用",)),
        ("等温壁で Tw・表・CSV のどれも無い", lambda: ld(iso()), ("どれも無い",)),
        ("等温壁の未知のキー", lambda: ld(iso(Tw_profile=S_DEMO, Tw_x=1.0)), ("未知のキー", "Tw_x")),
        ("断熱壁に表", lambda: ld({"thermal": "adiabatic", "Tw_profile": S_DEMO}), ("adiabatic", "未知のキー")),
        ("非単調 (s)", lambda: ld(iso(Tw_profile={"s": [-2.0, -1.0, 0.0, -0.5, 1.0], "Tw_K": [320, 330, 500, 400, 330]})),
         ("狭義単調増加",)),
        ("負の温度", lambda: ld(iso(Tw_profile={"s": [-2.0, 0.0, 1.0], "Tw_K": [320.0, -500.0, 330.0]})), ("正の温度",)),
        ("NaN (YAML の .nan)", lambda: ld(iso(Tw_profile={"s": [-2.0, 0.0, 1.0], "Tw_K": [320.0, float("nan"), 330.0]})), ("有限でない",)),
        ("長さ違い", lambda: ld(iso(Tw_profile={"s": [-2.0, 0.0, 1.0], "Tw_K": [320.0, 500.0]})), ("長さが違う",)),
        ("s が [−2, 1] を覆わない", lambda: ld(iso(Tw_profile={"s": [-1.5, 0.0, 1.0], "Tw_K": [320.0, 500.0, 330.0]})), ("[−2, 1]",)),
        ("x_m の表が設計スロート x = 0 を含まない", lambda: ld(iso(Tw_profile={"x_m": [0.1, 1.3], "Tw_K": [500.0, 330.0]})), ("x = 0",)),
        ("座標キーが 2 つ", lambda: ld(iso(Tw_profile={"s": [-2.0, 1.0], "x_m": [-0.5, 1.3], "Tw_K": [320.0, 330.0]})), ("どちらか 1 つ",)),
        ("表が mapping でない", lambda: ld(iso(Tw_profile=[[-2.0, 320.0], [1.0, 330.0]])), ("mapping",)),
        ("表の値が文字", lambda: ld(iso(Tw_profile={"s": [-2.0, 0.0, 1.0], "Tw_K": [320.0, "hot", 330.0]})), ("数値",)),
        ("CSV の列名違い (x,Tw_K)", lambda: ld(iso(Tw_profile_csv="hdr_x.csv")), ("列名", "s,Tw_K", "x_m,Tw_K")),
        ("CSV の列名違い (セミコロン区切り)", lambda: ld(iso(Tw_profile_csv="hdr_semicolon.csv")), ("列名", "カンマ区切り")),
        ("CSV の列名違い (s,Tw)", lambda: ld(iso(Tw_profile_csv="hdr_tw.csv")), ("列名",)),
        ("CSV が無い", lambda: ld(iso(Tw_profile_csv="missing.csv")), ("が無い", str((T / "missing.csv").resolve()))),
        ("CSV の NaN", lambda: ld(iso(Tw_profile_csv="nan.csv")), ("有限でない",)),
        ("CSV の列数違い", lambda: ld(iso(Tw_profile_csv="cols3.csv")), ("2 列", "2 行目")),
        ("CSV の数値でない値", lambda: ld(iso(Tw_profile_csv="text.csv")), ("数値でない", "3 行目")),
        ("CSV の s が [−2, 1] を覆わない", lambda: ld(iso(Tw_profile_csv="short.csv")), ("[−2, 1]",)),
        ("CSV が空", lambda: ld(iso(Tw_profile_csv="empty.csv")), ("空",)),
        ("CSV のパスが文字列でない", lambda: ld(iso(Tw_profile_csv=3)), ("パス",)),
    ]
    for name, fn, needles in neg:
        r, msg = raises(fn, *needles)
        check(f"負例: {name}", r, msg[:170])

# ---------------------------------------------------------------------------------------------- V1c 格子の見積もり
print("== V1c 格子の見積もり: スロート 500 K が同じで上流・下流が違う 2 つの表")
TABLES = {"A (上流 320/下流 330)": {"s": [-2.0, -1.0, -0.3, 0.0, 0.2, 0.5, 1.0], "Tw_K": [320.0, 330.0, 420.0, 500.0, 450.0, 380.0, 330.0]},
          "B (上流 450/下流 480)": {"s": [-2.0, -1.0, -0.3, 0.0, 0.2, 0.5, 1.0], "Tw_K": [450.0, 455.0, 480.0, 500.0, 495.0, 488.0, 480.0]}}
c0 = W.validate_card(card(wall=iso(Tw_profile=TABLES["A (上流 320/下流 330)"])))
p0 = W.card_to_problem(c0, RT)
chain = W._design_chain_of(p0)
wall = chain["wall"]
L_pipe, L_U, L_div = (float(c0["inlet"]["pipe_length_m"]), float(c0["inlet"]["contraction_length_m"]), float(c0["divergent_length_m"]))
# 区切りの点 [m]: 設計壁から (plan §4.3 [2026-10-13 訂正]; 正) と、カードの長さから (誤実装 (c))
G0 = p0["geometry"]
ANC = {"x_in_m": -(float(G0["L_pipe"]) * RT + float(G0["L_U"]) * RT), "x_contraction_m": -(float(G0["L_U"]) * RT), "x_throat_m": 0.0,
       "x_exit_m": float(wall.x_e) * RT}
ANC_CARD = {"x_in_m": -(L_pipe + L_U), "x_contraction_m": -L_U, "x_throat_m": 0.0, "x_exit_m": L_div}


def s_to_x_m(s: float, anc: dict) -> float:
    """試験側の区切りの写し: s ∈ [−2,−1] → 直管、[−1,0] → 縮流部、[0,1] → スロート→壁の終点 (区間ごとに線形)。"""
    x_in, x_c, x_t, x_e = anc["x_in_m"], anc["x_contraction_m"], anc["x_throat_m"], anc["x_exit_m"]
    if s <= -1.0:
        return x_in + (s + 2.0) * (x_c - x_in)
    if s <= 0.0:
        return x_c + (s + 1.0) * (x_t - x_c)
    return x_t + s * (x_e - x_t)


def same_anchors(a: dict, b: dict, rel: float = 1e-12) -> bool:
    return a is not None and set(a) == set(b) and all(abs(a[k] - b[k]) <= rel * max(1.0, abs(b[k])) for k in b)


def mu_air(T):
    """Sutherland (空気: μ_ref 1.716e-5 Pa s、T_ref 273.15 K、S 110.4 K)。"""
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


# 1 次元の流れ (estimate_mesh と同じ断面の並び: 設計壁の入口端〜終点の 8001 点、面積比 r_w²、亜音速枝 x < 0)
gas = W._gas(c0)
x = np.linspace(float(wall.x_in), float(wall.x_e), 8001)
rw = wall.r(x)
sarc = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(rw)))]) * RT      # 入口端からの壁の弧長 [m]
Mx = _invert_area_ratio_gas(np.maximum(rw, 1.0) ** 2, x >= 0.0, gas)
Te = np.where(Mx >= 1.0, gas.T_of_M(np.maximum(Mx, 1.0)), np.interp(np.minimum(Mx, 1.0), gas._Mu, gas._Tu))
Pe = float(c0["conditions"]["Pt"]) * _isentropic_pressure_ratio_gas(Te, gas)
Rg = float(gas.R)
ue = Mx * np.sqrt(gas.gamma(Te) * Rg * Te)


def y1_of(k: int, Tw: float) -> float:
    """断面 k の y1+ 0.8 の第 1 層 [m] (Eckert T* = T_e(1 + 0.032 M² + 0.58 (T_w/T_e − 1))、C_f* = 0.0592 Re_x*^−0.2)。"""
    Ts = Te[k] * (1.0 + 0.032 * Mx[k] ** 2 + 0.58 * (Tw / Te[k] - 1.0))
    rho_s, rho_w = Pe[k] / (Rg * Ts), Pe[k] / (Rg * Tw)
    Re = rho_s * ue[k] * sarc[k] / mu_air(Ts)
    tau = 0.0592 * Re ** -0.2 * 0.5 * rho_s * ue[k] ** 2
    return 0.8 * mu_air(Tw) / (rho_w * np.sqrt(tau / rho_w))


xm = x * RT
K = {"縮流部の中ほど (x ≈ −0.15 m)": int(np.argmin(np.abs(xm - (-0.5 * L_U)))),
     "拡大部の中ほど (x ≈ 0.6 m)": int(np.argmin(np.abs(xm - 0.5 * L_div)))}
print(f"   設計壁の区切りの点 [m] {ANC} (カードの長さでは {ANC_CARD}; 終点の差 {(ANC['x_exit_m'] - L_div) * 1e3:+.4f} mm)")
v1c = []
v1c_c = []
est_s = {}
for tname, tbl in TABLES.items():
    s_tab = tbl["s"]
    x_tab = [s_to_x_m(s, ANC) for s in s_tab]
    x_tab_card = [s_to_x_m(s, ANC_CARD) for s in s_tab]
    for kind in ("s", "x_m"):
        t_in = tbl if kind == "s" else {"x_m": x_tab, "Tw_K": tbl["Tw_K"]}
        ck = W.validate_card(card(wall=iso(Tw_profile=t_in)))
        pk = W.card_to_problem(ck, RT)
        est = W.estimate_mesh(ck, pk, chain=chain)
        sec = est["sections"]
        check(f"V1c {tname} {kind}: 見積もりの断面が試験と同じ並び (x [r_t] がビット一致)", np.array_equal(sec["x_rt"], x))
        if kind == "s":
            est_s[tname] = est
        else:
            dT = float(np.max(np.abs(sec["Tw_K"] - est_s[tname]["sections"]["Tw_K"])))
            check(f"V1c {tname}: x_m で与えた同じ表と s の表で局所の T_w が一致 (全断面 ≤ 1e-9 K)", dT <= 1e-9, f"{dT:.2e} K")
        for kname, k in K.items():
            Tw_exp = float(np.interp(xm[k], x_tab, tbl["Tw_K"]))         # 試験側: 区切りの写し → 実寸の x で線形補間
            y1_exp = y1_of(k, Tw_exp)
            Tw_imp, y1_imp = float(sec["Tw_K"][k]), float(sec["y1_m"][k])
            Tw_a = float(np.interp(0.0, s_tab, tbl["Tw_K"]))              # 誤実装 (a): 代表温度 (s = 0)
            Tw_b = float(np.interp(x[k], x_tab, tbl["Tw_K"]))             # 誤実装 (b): 無次元の x [r_t] を x_m に当てる (端値)
            y1_a, y1_b = y1_of(k, Tw_a), y1_of(k, Tw_b)
            rT, ry = abs(Tw_imp - Tw_exp) / Tw_exp, abs(y1_imp - y1_exp) / y1_exp
            da, db = abs(y1_a - y1_exp) / y1_exp, abs(y1_b - y1_exp) / y1_exp
            check(f"V1c {tname} {kind} {kname}: 局所の T_w と丸め前の y1 が期待値と相対 1e-8", rT <= 1e-8 and ry <= 1e-8,
                  f"T_w {Tw_imp:.6f} / 期待 {Tw_exp:.6f} K (相対 {rT:.1e})、y1 {y1_imp:.9e} / 期待 {y1_exp:.9e} m (相対 {ry:.1e})")
            check(f"V1c {tname} {kind} {kname}: 誤実装 (a) 代表温度・(b) 無次元の x はどちらも y1 が 1e-3 以上ずれる (検出できる)",
                  da >= 1e-3 and db >= 1e-3, f"(a) T_w {Tw_a:g} K → {da:.2e}、(b) T_w {Tw_b:g} K → {db:.2e}")
            v1c.append((tname, kind, kname, float(xm[k]), Tw_exp, Tw_imp, y1_exp, y1_imp, ry, Tw_a, da, Tw_b, db))
            if kind == "s":
                Tw_c = float(np.interp(xm[k], x_tab_card, tbl["Tw_K"]))      # 誤実装 (c): 区切りの点をカードの長さで取る
                v1c_c.append((tname, kname, Tw_exp, Tw_c, abs(y1_of(k, Tw_c) - y1_exp) / y1_exp))
        twr = est["wall_temperature"]
        check(f"V1c {tname} {kind}: 戻り値に表・見積もりの範囲・第 1 層を決めた断面の T_w (区切りの点は s の表だけ)",
              twr["coord"] == kind and twr["Tw_K_at_wall_first_frac"] == float(sec["Tw_K"][int(np.argmin(np.abs(x - est["x_wall_first_frac_rt"])))])
              and twr["range_m"] == [float(xm[0]), float(xm[-1])]
              and (same_anchors(twr["anchors_m"], ANC) and twr["anchors_source"].startswith("design wall")
                   if kind == "s" else twr["anchors_m"] is None), str(twr["anchors_m"]))
yA = est_s["A (上流 320/下流 330)"]["sections"]["y1_m"]
yB = est_s["B (上流 450/下流 480)"]["sections"]["y1_m"]
check("V1c: スロートの温度が同じ 2 つの表で、選んだ断面の y1 が違う (代表温度では区別できない差)",
      all(abs(yA[k] - yB[k]) / yA[k] > 1e-3 for k in K.values()), str([f"{yA[k]:.4e} vs {yB[k]:.4e}" for k in K.values()]))
# --- codex result レビュー 2026-10-13 m4: s の表は直管長 > 0 が要る (カードで早期に拒否し x_m を案内) ---
_inl = dict(BASE["inlet"]); _inl["pipe_length_m"] = 0.0
r, msg = raises(lambda: W.validate_card(card(inlet=_inl, wall=iso(Tw_profile=S_DEMO))), "直管の長さ > 0", "x_m")
check("m4: 直管長 0 のカードで s の表は拒否 (x_m を案内)", r, msg)
_xm = {"x_m": [-0.3, 0.0, 1.21], "Tw_K": [320.0, 500.0, 330.0]}
try:
    W.validate_card(card(inlet=_inl, wall=iso(Tw_profile=_xm))); ok_ = True; msg = ""
except ValueError as ex:
    ok_ = False; msg = str(ex)[:200]
check("m4: 直管長 0 のカードでも x_m の表は受ける", ok_, msg)

print("   V1c の数値 (表, 座標, 断面, x [m], T_w 期待/実装 [K], y1 期待/実装 [m], 相対差, (a) T_w・ずれ, (b) T_w・ずれ):")
for r in v1c:
    print(f"   {r[0]} {r[1]:>3} {r[2]}: x {r[3]:+.6f} m  T_w {r[4]:.6f}/{r[5]:.6f}  y1 {r[6]:.9e}/{r[7]:.9e}  相対 {r[8]:.1e}  "
          f"(a) {r[9]:g} K → {r[10]:.3e}  (b) {r[11]:g} K → {r[12]:.3e}")
# 誤実装 (c) 区切りの点をカードの長さで取る: このカードは設計壁の終点がカードの拡大部長より短いので、拡大部の断面で値がずれる
dc_div = [r[4] for r in v1c_c if r[1].startswith("拡大部")]
check("V1c 誤実装 (c) 区切りの点をカードの長さ: 拡大部の断面の y1 が V1c の許容 (相対 1e-8) を超えてずれる (値のずれで検出)",
      dc_div and min(dc_div) > 1e-8, str([f"{r[0][:1]} {r[1][:3]}: T_w {r[2]:.6f} → {r[3]:.6f} K、y1 {r[4]:.2e}" for r in v1c_c]))
print("   誤実装 (c) (表, 断面, T_w 期待 → カードの長さの写し [K], y1 のずれ): "
      + "; ".join(f"{r[0][:1]} {r[1][:6]} {r[2]:.6f} → {r[3]:.6f}, {r[4]:.2e}" for r in v1c_c))

# 一様のカード: T_w 一定、wall_temperature を返さない (戻り値・card_log は従来どおり)
cu = W.validate_card(card(wall=iso(Tw=350.0)))
eu = W.estimate_mesh(cu, W.card_to_problem(cu, RT), chain=chain)
check("一様のカード: 局所の T_w は全断面 350 K、wall_temperature を返さない、tw_note は従来の「等温壁 Tw 350.0 K」",
      bool(np.all(eu["sections"]["Tw_K"] == 350.0)) and "wall_temperature" not in eu and eu["inputs"]["wall"] == "等温壁 Tw 350.0 K")
# 設計壁の終点がカードの拡大部長とずれる (カードの拡大部長を ±1 mm; 同じ設計壁で見積もる): s の表は壁を覆い、区切りの点 = 設計壁の値、
# 局所の T_w はカードの長さに依らない。カードの長さで写す誤実装 (c) は、終点が長い側 (−1 mm のカード) では表が壁を覆わず例外になる
from forge_design import wall_temperature as WTP  # noqa: E402
base_Tw = est_s["A (上流 320/下流 330)"]["sections"]["Tw_K"]          # 表 A = S_DEMO
for dL in (-1.0e-3, +1.0e-3):
    cd = W.validate_card(card(wall=iso(Tw_profile=S_DEMO), divergent_length_m=L_div + dL))
    ok_e, why = True, ""
    try:
        ed = W.estimate_mesh(cd, W.card_to_problem(cd, RT), chain=chain)
    except ValueError as ex:
        ok_e, why, ed = False, str(ex)[:200], None
    if ed is not None:
        tw = ed["wall_temperature"]
        dT = float(np.max(np.abs(ed["sections"]["Tw_K"] - base_Tw)))
        ok_e = (same_anchors(tw["anchors_m"], ANC) and abs(tw["wall_end_minus_card_exit_m"] - (ANC["x_exit_m"] - (L_div + dL))) <= 1e-15
                and dT == 0.0)
        why = f"区切りの点 {tw['anchors_m']}、壁の終点 − カードの拡大部長 {tw['wall_end_minus_card_exit_m'] * 1e3:+.4f} mm、T_w の差 {dT:.1e} K"
    check(f"カードの拡大部長 {L_div + dL:.3f} m (設計壁の終点 {ANC['x_exit_m']:.6f} m): s の表は壁を覆い、区切りの点 = 設計壁の値 "
          "(相対 1e-12)・局所の T_w はカードの長さに依らない", ok_e, why)
    xt_card = [s_to_x_m(s, {**ANC_CARD, "x_exit_m": L_div + dL}) for s in S_DEMO["s"]]
    if ANC["x_exit_m"] > L_div + dL:
        r, msg = raises(lambda: WTP.check_covers(xt_card, float(xm[0]), float(xm[-1]), "格子の見積もりの範囲"), "覆わない")
        check(f"誤実装 (c) カードの長さで写すと、設計壁の終点が長い (+{(ANC['x_exit_m'] - L_div - dL) * 1e3:.3f} mm) ので表が壁を覆わず例外",
              r, msg[:200])
    else:
        k = K["拡大部の中ほど (x ≈ 0.6 m)"]
        Tc = float(np.interp(xm[k], xt_card, S_DEMO["Tw_K"]))
        Te_ = float(np.interp(xm[k], [s_to_x_m(s, ANC) for s in S_DEMO["s"]], S_DEMO["Tw_K"]))
        dy = abs(y1_of(k, Tc) - y1_of(k, Te_)) / y1_of(k, Te_)
        check(f"誤実装 (c) カードの長さで写すと、設計壁の終点が短い ({(ANC['x_exit_m'] - L_div - dL) * 1e3:+.3f} mm) ので拡大部の断面の y1 が"
              " 相対 1e-8 を超えてずれる", dy > 1e-8, f"T_w {Te_:.6f} → {Tc:.6f} K、y1 {dy:.2e}")
# x_m の表が設計壁の終点を覆わない
cxs = W.validate_card(card(wall=iso(Tw_profile={"x_m": [-0.5, 0.0, 1.0], "Tw_K": [320.0, 500.0, 330.0]})))
r, msg = raises(lambda: W.estimate_mesh(cxs, W.card_to_problem(cxs, RT), chain=chain), "格子の見積もりの範囲", "覆わない")
check("x_m の表が見積もりの範囲 (入口端〜設計壁の終点) を覆わなければ例外", r, msg[:200])

print(f"\n{'ALL OK' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
