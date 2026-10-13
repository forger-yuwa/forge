#!/usr/bin/env python3
"""帯局所抽出の等価半径の探索を積分域に限る修正 (§4.1)、δ_r の表の端と前の設計 → 今の設計への写し (§4.2)、要求カードの NS 後の補正
(§4.3 3) の試験 — plan tooling-nozzle-deltastar-band-root-zone §6 V1・V2・V4 と §5.1 #3・#4 (codex plan レビュー M1〜M3)。

  V1 (合成断面、`band_local_deficit`):
      (a) 帯の傾きが大きく q_ref がコアで負・全半径の F の末尾 < D → 旧実装は δ_r = r_w、新実装は解析的な真値 (多項式の積分と 3 次方程式の根)
          と相対 1e-6
      (b) 区間の中で F が単調でない (y_b を挟む外側の格子点で q_ref ≪ 0) → `nonmono_zone`・δ_r = NaN
      (c) codex の再現例 r_w 1・q_E 1・q_NS = max(−0.1 + 4y, 0)・y_b 0.1 → 旧は δ_r ≈ −4.98e13 (soft の negative_deficit)、新は `nonpos_ref`・NaN
      ほか: `HARD` に nonpos_ref・nonmono_zone。F が全半径で単調な断面 (既存の合成断面: 円環欠損・べき乗則・コアの波・負の δ_r) は新旧で
      結果がビット同一。
      旧実装はこの試験の中で新実装のソースから復元する (`old_band_local_deficit`: 探索域の参照の検査を外し、全半径の F で np.interp する)。
  M1 (codex plan レビュー M1) の伝播: 使用域 [0, x_e] に held の点がある `delta_r_next.csv` で `solve_rt_throat --prev` が止まる
      (設計チェーンの前に)。使用域の外 (x < 0) だけなら止まらない。
  V4 (表の端と写し): (a)(b) 終点が ε_f32 以内で内側・外側 → 揃う (実データ: run_0010 の終点は内側 1.54e-6 r_t・始点は外側 1.04e-7 r_t、
      run_0009 の終点は外側 2.7e-6 r_t で、どちらも寸法の解き (物理壁の構築) が通る) (c) ε_f32 を超えて内側・外側 → 例外 (d) 負の座標の入口端の
      丸め (内側・外側) (e) 前の設計の記録 (wall_repr.json) が無い・欠けている・S が prepare_info と違う → 例外 (f) 同じ設計への写しが恒等
      (ビット同一) (g) 長さを変えた設計への写しで区切りの点が新しい設計の区切りの点に写る (相対 1e-12)・s → x が `wall_temperature.to_x_m` と
      同じ値 (h) `solve_rt_throat --prev` (同じ設計) の r_t が写しの前後で同じ (case/70 sizing_ab.py の腕 B と同じ入力; 差 ≤ 1e-12 m)。
  --card (要求カードの NS 後の補正、§4.3 3): case/70 の c2b カードを `build_problem(ns_prev_run=...)` で作る (出力は一時ディレクトリ)。
      δ_E は run_0010 の CFD 前の補正 (delta_r_initial.csv の δ_r を同じ 5 次補間で NS の列 [run_0010 × run_0012 の抽出の x] に写したもの;
      表の端は抽出の float32 の丸めのまま) で、r_t が CFD 前のカードの r_t (0.024973561273227432 m; 内側の許容差 2.5e-8 m で止まった残差の
      分を戻した値) を 1e-12 m で再現することも見る。
      (i) 実寸の長さがカードと公差内 (寸法ゲート PASS) (ii) `_card.ns_correction` の sha256 が実ファイルと一致・表の端が最終の設計の入口端・
      出口と一致 (iii) `deltastar_loop --delta-r-csv` の prepare-only が通り、表を 1 値書き換えると run dir を作る前に止まる。
      **実データの δ_E (run_0010 × run_0012 / × run_0003) では、解いた r_t (約 0.02508 m) で設計チェーンの壁の検査 (設計壁区間で半径が非単調)
      が落ちてカードが作れない** — 実装の範囲外の既存の挙動なので、ここでは合成の δ_E で配線を確かめる (plan への報告事項)。
      FORGE_CONVERTER (変換器) が要る。
  --v2 NS_RUN EULER_RUN OUT (V2 後方互換): 旧実装 (計数つき) で band_select="edge" の抽出をして欠陥に当たる呼び出し (全半径の F の末尾 < D・
      探索区間の F が非単調; 計数は case/70 ab_root_zone.py の差し替えと同じ方法) を数え、新実装でも抽出して `delta_r_equiv.csv`・
      `delta_r_next.csv` がバイト同一かを見る (件数 0 の対で同一なら PASS)。探索域の参照が非正の呼び出し (新実装の nonpos_ref) は
      別に数え、旧実装がそこで返した値 (正・負・NaN) の内訳も出す。OUT は書き込み先 (scratch)。
  --case70 OUT: 新実装で case/70 run_0010 × run_0012 を抽出し、ab_root_zone.py の zone の腕 (`_ab_euler_steps/A/zone`) とバイト同一か。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_deltastar_root_zone_tests.py
    [--card] [--v2 NS_RUN EULER_RUN OUT] [--case70 OUT]
所要: 既定 約 1.5 分 (case/70 の設計チェーンと寸法の解き 4 回)。--card で約 3 分 (カードの固定点 + prepare_ns)、--v2・--case70 は
抽出 1 回 約 1 分。case/70 の run が無ければその項目は skip。
"""
from __future__ import annotations

import argparse
import filecmp
import hashlib
import inspect
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
from numpy.polynomial import polynomial as P

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
sys.path.insert(0, str(DESIGN))

import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

C70 = ROOT / "case" / "70.m5_air_wt_demo"
NS10 = C70 / "run_0010_ns_c2b_final"
NS09 = C70 / "run_0009_ns_c2_final"
EU12 = C70 / "run_0012_euler_c2b"
PROB_C2B = C70 / "problem_m5_air_c2b.yaml"
PROB_C2 = C70 / "problem_m5_air_c2.yaml"
CARD_C2B = C70 / "card_m5_air_demo_c2b.yaml"
NEXT_A_ZONE = C70 / "_ab_euler_steps" / "A" / "zone" / "delta_r_next.csv"      # run_0010 × run_0012 (ab_root_zone の zone の腕)
NEXT_B_ZONE = C70 / "_ab_root_zone" / "run_0010" / "zone" / "delta_r_next.csv"  # run_0010 × run_0003 (sizing_ab.py の腕 B の入力)
SIZING_AB = C70 / "_sizing_ab" / "sizing_ab.json"
R_THROAT = 0.025                       # c2b カードの物理スロート半径 [m]
RT_C2B_PRECFD = 0.024973561273227432   # c2b カード (CFD 前の固定点) の r_t [m]

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def skip(name, why):
    print(f"skip {name}  [{why}]", flush=True)


def raises(fn, *needles):
    """fn() が ValueError を出し、メッセージに needles を全部含むか。(結果, メッセージ)。"""
    try:
        fn()
    except ValueError as ex:
        msg = str(ex)
        return all(n in msg for n in needles), msg
    return False, "(例外なし)"


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ 旧実装の復元
# 新実装のソースから、探索域の参照の検査 (nonpos_ref) と区間の探索 (nonmono_zone) を外して、2026-10-13 以前の式 (全半径の F で np.interp)
# に戻す。counters を与えると、ab_root_zone.py と同じ方法で欠陥に当たる呼び出しを数える。
_OLD_SIGN = "    if D < 0:\n"
_OLD_SIGN_COUNTED = ("    _CNT['calls'] += 1\n"
                     "    if not (np.all(np.isfinite(qref[mz])) and np.all(qref[mz] > 0.0)):\n"
                     "        _CNT['nonpos_ref'] += 1\n"
                     "    if D < 0:\n")
_OLD_ROOT = ("            r_eff = float(np.interp(D, F_from_wall, r_desc))\n"
             "            delta_r = rwN - r_eff\n")
_OLD_ROOT_COUNTED = ("            _CNT['root'] += 1\n"
                     "            if float(F_from_wall[-1]) < D:\n"
                     "                _CNT['global_tail'] += 1\n"
                     "            _kz = int(np.searchsorted(rwN - r_desc, y_b, side='right'))\n"
                     "            if np.any(np.diff(F_from_wall[:_kz + 1]) < 0.0):\n"
                     "                _CNT['nonmono_zone'] += 1\n" + _OLD_ROOT)


def _replace_block(src: str, start: str, end: str, new: str) -> str:
    i, j = src.find(start), src.find(end)
    assert src.count(start) == 1 and src.count(end) == 1 and 0 <= i < j, f"復元の目印が見つからない/一意でない: {start!r} … {end!r}"
    return src[:i] + new + src[j + len(end):]


def old_band_local_deficit(counters: dict | None = None):
    src = inspect.getsource(DS.band_local_deficit)
    src = _replace_block(src, "    qz = qref[mz]\n", "    elif D < 0:\n", _OLD_SIGN if counters is None else _OLD_SIGN_COUNTED)
    src = _replace_block(src, "            y_desc = rwN - r_desc\n", "                delta_r = rwN - r_eff\n",
                         _OLD_ROOT if counters is None else _OLD_ROOT_COUNTED)
    assert 'reason.append("nonpos_ref")' not in src and 'reason.append("nonmono_zone")' not in src
    ns = dict(DS.__dict__)
    ns["_CNT"] = counters if counters is not None else {}
    exec(compile(src, "<old band_local_deficit (復元)>", "exec"), ns)
    return ns["band_local_deficit"]


def ones(rr):
    return np.ones_like(np.asarray(rr, dtype=float))


# ================================================================================================ V1
def test_v1():
    print("\n--- V1: 等価半径の探索を積分域に限る (合成断面) ---")
    old = old_band_local_deficit()
    # (a) 帯の傾きが大きく q_ref がコアで負: q_ref = 1 − 5 (y − y_b)、境界層は y < δ99 で q_NS = q_ref·y/δ99 (線形欠損)
    rw, N, yb, d99, a0, b = 1.0, 20000, 0.10003, 0.05, 1.0, -5.0
    r = np.linspace(0.0, rw, N + 1); y = rw - r
    qref = a0 + b * (y - yb)
    qN = np.where(y <= 0.2, qref * np.minimum(y / d99, 1.0), np.maximum(qref, 0.05))
    new_a = DS.band_local_deficit(r, qN, ones, rw, y_b_fixed=yb, refine=1)
    old_a = old(r, qN, ones, rw, y_b_fixed=yb, refine=1)
    # 解析的な真値: D = 2π∫_0^δ99 q_ref (1 − y/δ99)(1 − y) dy、2π∫_0^δ q_ref (1 − y) dy = D の根 (0 < δ < y_b)
    qp = np.array([a0 - b * yb, b])
    G = P.polyint(P.polymul(P.polymul(qp, [1.0, -1.0 / d99]), [1.0, -1.0]))
    D_an = 2 * np.pi * (P.polyval(d99, G) - P.polyval(0.0, G))
    Fi = 2 * np.pi * P.polyint(P.polymul(qp, [1.0, -1.0]))
    F_full_end = float(P.polyval(rw, Fi) - P.polyval(0.0, Fi))
    c = Fi.copy(); c[0] -= D_an
    roots = [float(z.real) for z in P.polyroots(c) if abs(z.imag) < 1e-12 and 0.0 < z.real < yb]
    check("V1(a) 前提: 全半径の F の末尾 < D (コアで q_ref < 0)", len(roots) == 1 and F_full_end < D_an and float(qref[0]) < 0.0,
          f"F 末尾 {F_full_end:.4f}、D {D_an:.6f}、軸の q_ref {float(qref[0]):.3f}")
    check("V1(a) 旧実装は δ_r = r_w (np.interp が軸の半径を返す)", abs(old_a["delta_r"] - rw) < 1e-12, f"旧 δ_r {old_a['delta_r']!r}")
    rel = abs(new_a["delta_r"] / roots[0] - 1.0)
    check("V1(a) 新実装は解析的な真値と相対 1e-6", rel <= 1e-6 and new_a["reason"] == [],
          f"新 {new_a['delta_r']!r}、真値 {roots[0]!r}、相対 {rel:.2e}、D 数値/解析 {new_a['mass_deficit'] / D_an - 1:+.1e}")
    # (b) y_b を挟む外側の格子点で q_ref ≪ 0 → y_b までの F (線形補間の点) が減る → nonmono_zone
    yy = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.099, 0.11, 0.12, 0.13, 0.14, 0.2, 0.4, 0.6, 0.8, 1.0])
    rb = np.sort(1.0 - yy); yv = 1.0 - rb
    qrb = 1.0 - 1000.0 * (yv - 0.1)
    qNb = np.where(yv <= 0.15, qrb, 1.0)
    qNb = np.where(yv < 0.015, 0.5 * qrb, qNb)
    nb = DS.band_local_deficit(rb, qNb, ones, 1.0, y_b_fixed=0.1, refine=1)
    zone = qrb[yv <= 0.1]
    check("V1(b) 前提: 探索域 (y ≤ y_b) の格子点の q_ref は正、y_b の外側の点は負", bool(np.all(zone > 0) and qrb[np.argmin(np.abs(yv - 0.11))] < 0))
    check("V1(b) 区間の中で F が減る → nonmono_zone・NaN", nb is not None and np.isnan(nb["delta_r"]) and nb["reason"] == ["nonmono_zone"],
          f"{None if nb is None else (nb['delta_r'], nb['reason'], nb['mass_deficit'])}")
    # (c) codex の再現例 (notes/reviews/2026-10-13-tooling-nozzle-deltastar-band-root-zone-plan.md M2)
    rc = np.linspace(0.0, 1.0, 401); yc = 1.0 - rc
    qNc = np.maximum(-0.1 + 4.0 * yc, 0.0)
    oc = old(rc, qNc, ones, 1.0, y_b_fixed=0.1)
    nc = DS.band_local_deficit(rc, qNc, ones, 1.0, y_b_fixed=0.1)
    check("V1(c) 旧: δ_r ≈ −4.98e13 (soft の negative_deficit)", abs(oc["delta_r"] / -4.979e13 - 1.0) < 1e-3 and oc["reason"] == ["negative_deficit"],
          f"旧 {oc['delta_r']:.4e} {oc['reason']}、D {oc['mass_deficit']:.6e}")
    check("V1(c) 新: nonpos_ref (hard)・NaN", np.isnan(nc["delta_r"]) and nc["reason"] == ["nonpos_ref"], f"新 {nc['delta_r']} {nc['reason']}")
    src = inspect.getsource(DS.deltastar_from_core_matched_euler)
    check("V1 HARD に nonpos_ref・nonmono_zone", 'HARD = ("no_root", "nan", "no_plateau", "nonpos_ref", "nonmono_zone")' in src)
    # 単調な断面はビット同一 (run_deltastar_tests.py と同じ合成断面: 円環欠損・べき乗則・コアの波・壁半径差・負の δ_r)
    rw = 5.0
    n, ff = 97, 4.5e-5
    lo, hi = 1.0001, 1.5
    for _ in range(100):
        g = 0.5 * (lo + hi)
        lo, hi = (g, hi) if ff * (g ** (n - 1) - 1) / (g - 1) < 1.0 else (lo, g)
    yg = np.concatenate([[0.0], np.cumsum(ff * g ** np.arange(n - 1))]); yg /= yg[-1]
    rg = (rw * (1.0 - yg))[::-1]
    qE = lambda r_: 1.0 + 0.2 * (np.asarray(r_) / rw) ** 2
    cases = []
    for d_true in (0.01, 0.05, 0.2):
        q = 1.3 * qE(rg); q[rg > rw - d_true] = 0.0
        cases += [("円環欠損", rg, q, qE, rw, dict(delta_in=d_true)), ("円環欠損 δ_in×0.3", rg, q, qE, rw, dict(delta_in=0.3 * d_true))]
    yb_ = rw - rg
    bl = qE(rg) * np.where(yb_ < 0.4, (np.maximum(yb_, 0) / 0.4) ** (1 / 7), 1.0)
    cases.append(("べき乗則", rg, bl, qE, rw, dict(delta_in=0.03)))
    wave = qE(rg) * (1.0 + 0.01 * np.sin(2 * np.pi * rg / rw)) * np.where(yb_ < 0.4, (np.maximum(yb_, 0) / 0.4) ** (1 / 7), 1.0)
    cases.append(("コアの波", rg, wave, qE, rw, dict(delta_in=0.03)))
    qx = np.full_like(rg, 2.0); qx[rg >= rw - 0.05] = 3.0
    cases.append(("負の δ_r (冷却壁)", rg, qx, lambda r_: np.full_like(np.asarray(r_, float), 2.0), rw, dict(delta_in=0.0)))
    for yfix in (0.05, 0.11, 0.3, 0.6):
        cases.append((f"帯固定 y_b {yfix}", rg, bl, qE, rw, dict(y_b_fixed=yfix)))
    def eq(u, v):
        if isinstance(u, (list, str)):
            return u == v
        return np.array_equal(np.asarray(u, dtype=float), np.asarray(v, dtype=float), equal_nan=True)
    same = []
    for name, r_, q_, e_, rw_, kw in cases:
        a_, b_ = old(r_, q_, e_, rw_, **kw), DS.band_local_deficit(r_, q_, e_, rw_, **kw)
        ok = (a_ is None and b_ is None) or (a_ is not None and b_ is not None and a_.keys() == b_.keys()
                                             and all(eq(a_[k], b_[k]) for k in a_))
        same.append((name, ok, None if b_ is None else b_["delta_r"]))
    check("V1 F が全半径で単調な断面は新旧でビット同一 (全項目)", all(s[1] for s in same), f"{[(s[0], s[1]) for s in same if not s[1]] or len(same)} 断面")
    neg = [s for s in same if s[0].startswith("負")][0]
    check("V1 正の参照に対する負の δ_r は維持", neg[2] is not None and neg[2] < 0.0, f"{neg[2]!r}")


# ================================================================================================ V4 (単体)
A10 = {"x_in": -20.021173373299316, "x_contraction": -12.01270402397959, "x_throat": 0.0, "x_exit": 48.031377427318134}
S10 = 0.024973561273227432


def test_v4_unit():
    print("\n--- V4: 表の端の照合と揃え・写し (単体) ---")
    eps_e, eps_s = DL.eps_f32(A10["x_exit"], S10), DL.eps_f32(A10["x_in"], S10)
    check("V4 ε_f32 (出口) = plan §3 の値 8.588e-6 r_t", abs(eps_e - 8.58811697323306e-06) < 1e-18, f"{eps_e!r}、入口 {eps_s!r}")
    base = np.linspace(A10["x_in"], A10["x_exit"], 1800)

    def tab(d0=0.0, d1=0.0):
        x = base.copy(); x[0] += d0; x[-1] += d1
        return x
    for nm, d1 in (("(a) 終点が内側 0.5 ε", -0.5 * eps_e), ("(b) 終点が外側 0.5 ε", 0.5 * eps_e), ("(a) 終点が内側 0.999 ε", -0.999 * eps_e)):
        x, rec = DL.snap_table_ends(tab(0.0, d1), A10, S10)
        check(f"V4{nm} → 終点が前の設計の出口に揃う (値の数は同じ)", x[-1] == A10["x_exit"] and len(x) == len(base)
              and np.array_equal(x[:-1], tab(0.0, d1)[:-1]) and rec["end"]["snapped"] and abs(rec["end"]["diff_rt"] - d1) < 1e-13,
              f"差 {rec['end']['diff_rt']:+.3e}")
    for nm, d1 in (("内側", -2.0 * eps_e), ("外側", 2.0 * eps_e)):
        ok, msg = raises(lambda d1=d1: DL.snap_table_ends(tab(0.0, d1), A10, S10), "終点", "ε_f32")
        check(f"V4(c) 終点が ε を超えて{nm} → 例外", ok, msg[:120])
    for nm, d0 in (("外側 (より負)", -0.5 * eps_s), ("内側", 0.5 * eps_s)):
        x, rec = DL.snap_table_ends(tab(d0, 0.0), A10, S10)
        check(f"V4(d) 負の座標の入口端の丸め ({nm}) → 揃う", x[0] == A10["x_in"] and rec["start"]["snapped"], f"差 {rec['start']['diff_rt']:+.3e}")
    ok, msg = raises(lambda: DL.snap_table_ends(tab(-2.0 * eps_s, 0.0), A10, S10), "始点", "ε_f32")
    check("V4(c) 始点が ε を超える → 例外", ok, msg[:120])
    x, rec = DL.snap_table_ends(base, A10, S10)
    check("V4 端が設計と一致する表は揃えても同じ (snapped = False)", np.array_equal(x, base) and not rec["start"]["snapped"] and not rec["end"]["snapped"])
    # (e) 前の設計の記録
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "prepare_info.json").write_text(json.dumps({"scale_m": S10}))
        ok, msg = raises(lambda: DL.prev_design_anchors(td), "wall_repr.json", "無い")
        check("V4(e) 前の設計の壁ファイルが無い → 例外", ok, msg[:100])
        if (NS10 / "wall_repr.json").exists():
            rec = json.loads((NS10 / "wall_repr.json").read_text())
            bad = json.loads(json.dumps(rec)); del bad["design_wall"]["upstream_hermite"]
            (td / "wall_repr.json").write_text(json.dumps(bad))
            ok, msg = raises(lambda: DL.prev_design_anchors(td), "読めない")
            check("V4(e) 壁ファイルに縮流部の記録が無い → 例外", ok, msg[:100])
            (td / "wall_repr.json").write_text(json.dumps(rec))
            (td / "prepare_info.json").write_text(json.dumps({"scale_m": S10 * (1 + 1e-12)}))
            ok, msg = raises(lambda: DL.prev_design_anchors(td), "scale_m")
            check("V4(e) 壁ファイルと prepare_info の S が違う → 例外", ok, msg[:100])
            (td / "prepare_info.json").write_text(json.dumps({"scale_m": S10}))
            a = DL.prev_design_anchors(td)
            check("V4 前の設計の区切りの点 = run_0010 の壁ファイルの値", all(a[k] == A10[k] for k in DL.ANCHOR_KEYS) and a["scale_m"] == S10)
        else:
            skip("V4(e) 壁ファイルの欠け", f"{NS10} が無い")
    # (f) 同じ設計への写し
    xs = np.concatenate([np.linspace(A10["x_in"] - 0.5, A10["x_exit"] + 0.5, 4001), [A10["x_contraction"], -1e-3, 1e-3]])
    xs.sort()
    xm, ident = DL.map_x_between_designs(xs, A10, dict(A10))
    check("V4(f) 同じ設計への写しは恒等 (ビット同一)", ident and np.array_equal(xm, xs))
    xg = DL.s_to_x(DL.x_to_s(xs, A10), A10)
    err = float(np.max(np.abs(xg - xs)))
    check("V4(f) 一般の式 (x → s → x) でも同じ設計なら |誤差| ≤ 1e-14 r_t", err <= 1e-14, f"最大 {err:.2e} r_t")
    # (g) 長さを変えた設計への写し (case/70 の寸法の A/B の B: r_t ×1.004091 で実寸を保つと r_t 単位の長さは 1/1.004091)
    k = 1.0 / 1.0040908933879382
    A_new = {"x_in": A10["x_in"] * k, "x_contraction": A10["x_contraction"] * k, "x_throat": 0.0, "x_exit": A10["x_exit"] * 0.99637}
    pts = np.array([A10[kk] for kk in DL.ANCHOR_KEYS])
    xm, ident = DL.map_x_between_designs(pts, A10, A_new)
    want = np.array([A_new[kk] for kk in DL.ANCHOR_KEYS])
    rel = float(np.max(np.abs(xm - want) / np.maximum(np.abs(want), 1.0)))
    check("V4(g) 長さを変えた設計: 区切りの点が新しい設計の区切りの点に写る (相対 1e-12)", not ident and rel <= 1e-12, f"{rel:.1e}")
    mids = np.array([0.5 * (pts[i] + pts[i + 1]) for i in range(3)])
    xm_mid, _ = DL.map_x_between_designs(mids, A10, A_new)
    want_mid = np.array([0.5 * (want[i] + want[i + 1]) for i in range(3)])
    check("V4(g) 区間の中点は中点に写る (区間ごとの線形)", float(np.max(np.abs(xm_mid - want_mid))) <= 1e-12)
    out_s = np.concatenate([[-2.3], np.linspace(-2.0, 1.0, 301), [1.2]])
    from forge_design.wall_temperature import to_x_m
    xt, _ = to_x_m({"coord": "s", "coords": out_s, "Tw_K": np.ones_like(out_s)},
                   {"x_in_m": A_new["x_in"], "x_contraction_m": A_new["x_contraction"], "x_throat_m": 0.0, "x_exit_m": A_new["x_exit"]})
    check("V4(g) s → x は wall_temperature.to_x_m と同じ値 (区間の外の延長も)", np.array_equal(DL.s_to_x(out_s, A_new), xt))


# ================================================================================================ M1・V4 (実データ)
def _copy_prev(dst: Path, src_run: Path) -> Path:
    dst.mkdir(parents=True, exist_ok=True)
    for f in ("prepare_info.json", "wall_repr.json"):
        shutil.copy2(src_run / f, dst / f)
    return dst


def _write_next_with_held(src_csv: Path, dst_csv: Path, held_idx) -> None:
    lines = src_csv.read_text().splitlines()
    for i in held_idx:
        f = lines[1 + i].split(",")
        f[4] = "1.000000000000000000e+00"
        lines[1 + i] = ",".join(f)
    dst_csv.write_text("\n".join(lines) + "\n")


def test_real(tmp: Path):
    print("\n--- M1 の伝播・V4 (実データ: case/70) ---")
    need = [NS10 / "wall_repr.json", NS10 / "prepare_info.json", NEXT_A_ZONE, NEXT_B_ZONE, PROB_C2B]
    if not all(p.exists() for p in need):
        skip("M1・V4 (実データ)", f"case/70 の入力が無い: {[str(p) for p in need if not p.exists()]}")
        return
    nx = DL.read_delta_r_next(NEXT_A_ZONE)
    x = nx["x_rt"]
    # M1: 使用域に held → 止まる (設計チェーンの前に)
    prev = _copy_prev(tmp / "prev_m1", NS10)
    i2 = int(np.argmin(np.abs(x - 2.0)))
    _write_next_with_held(NEXT_A_ZONE, prev / "delta_r_next.csv", [i2])
    import time
    t0 = time.time()
    ok, msg = raises(lambda: DL.solve_rt_throat(PROB_C2B, R_THROAT, prev_run=prev, delta_r_out=tmp / "m1_stop.csv"), "held", "使用域", "1 点")
    check("M1 使用域 [0, x_e] に held が 1 点 → solve_rt_throat --prev が止まる (設計チェーンの前に)", ok and time.time() - t0 < 5.0,
          f"{msg[:140]} ({time.time() - t0:.1f} s)")
    _write_next_with_held(NEXT_A_ZONE, prev / "delta_r_next.csv", [int(np.argmin(np.abs(x - 48.0)))])
    ok, msg = raises(lambda: DL.load_ns_correction(prev), "held")
    check("M1 出口の近く (x ≈ 48) の held も止まる", ok, msg[:80])
    neg = [int(i) for i in np.where((x > -6.0) & (x < -4.0))[0][:5]]
    _write_next_with_held(NEXT_A_ZONE, prev / "delta_r_next.csv", neg)
    r = DL.solve_rt_throat(PROB_C2B, R_THROAT, prev_run=prev, delta_r_out=tmp / "m1_go.csv")
    nc = r["ns_correction"]
    check("M1 使用域の外 (x < 0) だけの held では止まらない (寸法が解ける; held は記録)",
          np.isfinite(r["r_t_m"]) and nc["held"]["n_held"] == 5 and nc["held"]["n_held_in_use"] == 0,
          f"r_t {r['r_t_m']!r} m、held {nc['held']['n_held']} 点 (x {nc['held']['x_held_outside_use_rt'][:2]}…)")
    e = nc["table_ends"]
    check("V4(a) 実データ (run_0010 × run_0012): 終点は内側 1.54e-6 r_t・始点は外側 1.04e-7 r_t で揃い、物理壁が作れる",
          e["end"]["snapped"] and -1.6e-6 < e["end"]["diff_rt"] < -1.5e-6 and -1.1e-7 < e["start"]["diff_rt"] < -1.0e-7
          and r["physical_wall_repr"] == "SingleBSplinePhysicalWall",
          f"終点 {e['end']['diff_rt']:+.3e}・始点 {e['start']['diff_rt']:+.3e} r_t、壁 {r['physical_wall_repr']}")
    check("V4(f) 同じ設計 (c2b の YAML と run_0010 の壁ファイル) の写しは恒等", nc["x_map"].startswith("identity"), nc["x_map"])
    tbl = np.loadtxt(tmp / "m1_go.csv", delimiter=",", skiprows=1)
    check("V4 書いた表の端 = 設計の入口端・出口 (ビット同一)", tbl[0, 0] == A10["x_in"] and tbl[-1, 0] == A10["x_exit"]
          and np.array_equal(tbl[1:-1, 0], x[1:-1]))
    # V4(b) 外側の丸め: run_0009 の抽出 (終点が出口を 2.7e-6 r_t 超える; codex plan レビュー M3)。δ_E は修正後の抽出 (ab_root_zone の zone の腕;
    # run_0009 の delta_r_next.csv は旧実装の抽出で x ≈ 1〜3 に δ_r = r_w の誤抽出を含み、物理壁の上流のゲートが落ちる)
    next09 = C70 / "_ab_root_zone" / "run_0009" / "zone" / "delta_r_next.csv"
    if (NS09 / "wall_repr.json").exists() and next09.exists() and PROB_C2.exists():
        r9 = DL.solve_rt_throat(PROB_C2, R_THROAT, prev_run=NS09, delta_next=next09, delta_r_out=tmp / "r09.csv")
        e9 = r9["ns_correction"]["table_ends"]["end"]
        check("V4(b) 実データ (run_0009): 終点が外側 2.7e-6 r_t → 揃って寸法が解ける", e9["snapped"] and 2.6e-6 < e9["diff_rt"] < 2.8e-6
              and np.isfinite(r9["r_t_m"]), f"終点 {e9['diff_rt']:+.3e} r_t (ε {e9['eps_f32_rt']:.2e})、r_t {r9['r_t_m']!r}")
    else:
        skip("V4(b) run_0009", "run_0009 か problem_m5_air_c2.yaml が無い")
    # V4(h) 写しの前後で r_t が同じ: 旧の供給 (表の終点だけを x_e に置き換えた写し、sizing_ab.py の腕 B) と新の経路
    from forge_design.evaluate.runner_axismach import delta_r_from_table
    nb = DL.read_delta_r_next(NEXT_B_ZONE)
    xo = np.asarray(nb["x_rt"], dtype=float).copy()
    gap = A10["x_exit"] - xo[-1]
    xo[-1] = A10["x_exit"]
    dEo = np.asarray(nb["delta_E"], dtype=float)
    old_tab = tmp / "h_old.csv"

    def old_supplier(d_, rt):
        fac = (rt / S10) ** -0.2
        np.savetxt(old_tab, np.c_[xo, dEo * fac], delimiter=",", comments="", header="x_rt,delta_r (旧の供給)")
        t = np.loadtxt(old_tab, delimiter=",", skiprows=1)
        return delta_r_from_table(t[:, 0], t[:, 1]), {"kind": "旧の供給 (試験)"}
    r_old = DL.solve_rt_throat(PROB_C2B, R_THROAT, _delta_supplier=old_supplier)
    r_new = DL.solve_rt_throat(PROB_C2B, R_THROAT, prev_run=NS10, delta_next=NEXT_B_ZONE, delta_r_out=tmp / "h_new.csv")
    d_rt = abs(r_new["r_t_m"] - r_old["r_t_m"])
    ref = json.loads(SIZING_AB.read_text())["B"]["r_t_m"] if SIZING_AB.exists() else None
    check("V4(h) solve_rt_throat --prev (同じ設計) の r_t が写しの前後で同じ (差 ≤ 1e-12 m)", d_rt <= 1e-12,
          f"新 {r_new['r_t_m']!r}、旧 {r_old['r_t_m']!r}、差 {d_rt:.1e} m、終点の差 {gap:.3e} r_t、sizing_ab.json の B {ref!r}")
    if ref is not None:
        check("V4(h) sizing_ab.json の腕 B の r_t と同じ", abs(r_new["r_t_m"] - ref) <= 1e-12, f"{r_new['r_t_m'] - ref:.1e} m")


# ================================================================================================ カード (--card)
def test_card(tmp: Path):
    print("\n--- カードの NS 後の補正 (--card) ---")
    need = [CARD_C2B, NS10 / "wall_repr.json", NS10 / "delta_r_initial.csv", NEXT_A_ZONE, EU12]
    if not all(p.exists() for p in need):
        skip("カード", f"入力が無い: {[str(p) for p in need if not p.exists()]}")
        return
    from forge_design.evaluate.runner_axismach import delta_r_from_table, design_chain
    from forge_design.intake import wt_card as W
    from forge_design.probdef import load_problem
    # 合成の δ_E: run_0010 の CFD 前の補正 (delta_r_initial.csv) を同じ 5 次補間で NS の列の x に写す (端は抽出の float32 の丸めのまま)
    ini = np.loadtxt(NS10 / "delta_r_initial.csv", delimiter=",", skiprows=1)
    f_ini = delta_r_from_table(ini[:, 0], ini[:, 1])
    nx = DL.read_delta_r_next(NEXT_A_ZONE)
    xs = np.asarray(nx["x_rt"], dtype=float)
    prev = _copy_prev(tmp / "prev_card", NS10)
    dE = f_ini(xs)
    np.savetxt(prev / "delta_r_next.csv", np.c_[xs, dE, np.zeros_like(xs), dE, np.zeros_like(xs)], delimiter=",", comments="",
               header=",".join(DL.DELTA_R_NEXT_COLUMNS) + " (試験: run_0010 の CFD 前の δ_r を NS の列に写した合成の δ_E)")
    out = tmp / "card" / "problem_c2b_ns.yaml"
    logd = W.build_problem(CARD_C2B, out, ns_prev_run=prev, log=lambda *a, **k: None)
    ns = logd.get("ns_correction") or {}
    rt = logd["r_t"]["r_t_m"]
    check("カード: r_t は NS 後の補正の固定点 (method)", logd["r_t"]["method"] == "physical_throat_ns_corrected", logd["r_t"]["method"])
    # CFD 前のカードの r_t は内側の許容差 2.5e-8 m で止まっていて残差 res (物理スロート半径 − R) を持つ。NS 補正の経路 (許容差 1e-9 m) は
    # 同じ δ なら r_t ≈ r_t,CFD前 · R/(R + res) になるはず
    log_c2b = PROB_C2B.with_name(PROB_C2B.name + ".card_log.json")
    res_pre = (json.loads(log_c2b.read_text())["r_t"]["history"][-1]["inner_residual_m"] if log_c2b.exists() else None)
    if res_pre is not None:
        want = RT_C2B_PRECFD * R_THROAT / (R_THROAT + res_pre)
        check("カード: CFD 前と同じ δ を NS 補正で与えると CFD 前のカードの r_t を再現 (その残差の分を戻して |Δ| ≤ 1e-12 m)",
              abs(rt - want) <= 1e-12, f"r_t {rt!r}、期待 {want!r} (CFD 前 {RT_C2B_PRECFD!r}・残差 {res_pre:.3e} m)、差 {rt - want:+.1e} m、"
              f"外側 {len(logd['r_t']['history'])} 回")
    else:
        check("カード: CFD 前と同じ δ を NS 補正で与えると CFD 前のカードの r_t を再現 (|Δ| ≤ 2.5e-8 m = CFD 前の内側の許容差)",
              abs(rt - RT_C2B_PRECFD) <= 2.5e-8, f"r_t {rt!r}、差 {rt - RT_C2B_PRECFD:+.2e} m")
    rows = {r_["name"]: r_ for r_ in logd["dimensions"]["rows"]}
    lens = ("pipe_length", "contraction_length", "divergent_length")
    check("カード (i) 実寸の長さがカードと公差内・寸法ゲート PASS", logd["dimensions"]["pass"] and all(rows[k]["pass"] for k in lens),
          ", ".join(f"{k} {rows[k]['diff_m'] * 1e3:+.3f} mm" for k in lens) +
          f", 物理スロート径 {rows['throat_diameter_physical']['diff_m'] * 1e3:+.4f} mm")
    tp = Path(ns.get("table", "")) if ns else None
    shas = {"delta_next": (ns.get("delta_next"), ns.get("delta_next_sha256")), "wall_repr": (ns.get("wall_repr"), ns.get("wall_repr_sha256")),
            "prepare_info": (ns.get("prepare_info"), ns.get("prepare_info_sha256")), "table": (ns.get("table"), ns.get("table_sha256"))}
    bad = [k for k, (p_, s_) in shas.items() if not p_ or not Path(p_).exists() or sha256(p_) != s_]
    check("カード (ii) _card.ns_correction の sha256 が実ファイルと一致 (delta_next・壁ファイル・prepare_info・表)", not bad and tp is not None
          and tp == out.with_name(out.stem + "_delta_r.csv"), f"不一致 {bad}" if bad else f"表 {tp.name}")
    import yaml
    card_rec = yaml.safe_load(out.read_text())["_card"]["ns_correction"]
    check("カード (ii) YAML の _card.ns_correction = card_log の記録", card_rec == json.loads(json.dumps(ns)))
    p = load_problem(out); d = design_chain(p)
    a_new = DL.design_anchors(d["wall"])
    t = np.loadtxt(tp, delimiter=",", skiprows=1)
    check("カード (ii) 表の端 = 最終の設計の入口端・出口 (ビット同一)", t[0, 0] == a_new["x_in"] and t[-1, 0] == a_new["x_exit"],
          f"表 [{t[0, 0]!r}, {t[-1, 0]!r}]、設計 [{a_new['x_in']!r}, {a_new['x_exit']!r}]、揃えた差 始点 "
          f"{ns['table_ends']['start']['diff_rt']:+.2e}・終点 {ns['table_ends']['end']['diff_rt']:+.2e} r_t")
    check("カード: 表の大きさ = δ_E × (r_t/S_prev)^−0.2", np.allclose(t[:, 1], dE * (rt / S10) ** -0.2, rtol=0, atol=1e-15)
          and ns["factor"] == (rt / S10) ** -0.2, f"factor {ns['factor']!r}、最後の内側の表と同一 {ns['final_table_equals_last_inner']}")
    # (iii) deltastar_loop --delta-r-csv の prepare-only と、表の書き換えの検出
    if not (os.environ.get("FORGE_CONVERTER") or os.environ.get("FORGE_BIN")):
        skip("カード (iii) prepare-only", "FORGE_CONVERTER / FORGE_BIN が無い")
        return
    run = tmp / "run_prep"
    argv = ["--problem", str(out), "--euler-ref", str(EU12), "--prev", str(NS10), "--run-dir", str(run), "--delta-r-csv", str(tp),
            "--prepare-only"]
    rc = DL.main(argv)
    info = json.loads((run / "prepare_info.json").read_text())
    r_phys = float(info["throat_physical"]["r"]) * float(info["scale_m"])
    check("カード (iii) --delta-r-csv の prepare-only が通る (表の照合済み・物理スロート半径 ≈ カード)",
          rc == 0 and info["delta_r_table"]["checked"] and info["delta_r_table"]["sha256"] == ns["table_sha256"]
          and abs(r_phys - R_THROAT) <= 1e-8, f"物理スロート半径 {r_phys!r} m (差 {r_phys - R_THROAT:+.1e})、{info['dstar_source'][:60]}")
    lines = tp.read_text().splitlines()
    f = lines[900].split(","); f[1] = repr(float(f[1]) * (1.0 + 1e-9)); lines[900] = ",".join(f)
    bad_tab = tmp / "tampered_delta_r.csv"
    bad_tab.write_text("\n".join(lines) + "\n")
    run2 = tmp / "run_prep_tampered"
    argv2 = ["--problem", str(out), "--euler-ref", str(EU12), "--prev", str(NS10), "--run-dir", str(run2), "--delta-r-csv", str(bad_tab),
             "--prepare-only"]
    ok, msg = raises(lambda: DL.main(argv2), "sha256", "_card.ns_correction")
    check("カード (iii) 表を 1 値書き換えると止まる (run dir を作る前)", ok and not run2.exists(), msg[:120])


# ================================================================================================ V2・case/70
def extract_old_new(ns_run, eu_run, out_dir: Path, compare: bool = True) -> dict:
    """旧実装 (計数つき) と新実装で同じ抽出 (band_select="edge") をして、欠陥に当たる呼び出しの件数と出力の同一性を返す。
    欠陥の件数 (`defect_calls`) は ab_root_zone.py と同じ方法 (根の探索の行で、全半径の F の末尾 < D・探索区間の F が非単調)。
    探索域の参照が非正の呼び出し (新実装で nonpos_ref) は別に数え、そのうち旧実装が正の有限値を返した (= 梯子・最終の値に使われうる)
    件数も数える (edge の呼び出しは帯を固定するので再試行の再帰は無い)。"""
    cnt = {"calls": 0, "nonpos_ref": 0, "root": 0, "global_tail": 0, "nonmono_zone": 0}
    new = DS.band_local_deficit
    old = old_band_local_deficit(cnt)
    ns_new = dict(DS.__dict__)
    exec(compile(inspect.getsource(new), "<new band_local_deficit (計数用の写し)>", "exec"), ns_new)
    new_copy = ns_new["band_local_deficit"]
    npos = {"nonpos_ref_new": 0, "nonpos_ref_old_positive": 0, "nonpos_ref_old_negative": 0, "nonpos_ref_old_nan": 0}

    def counted(*a, **k):
        ro = old(*a, **k)
        rn = new_copy(*a, **k)
        if rn is not None and "nonpos_ref" in rn["reason"]:
            npos["nonpos_ref_new"] += 1
            v = None if ro is None else ro["delta_r"]
            key = ("nonpos_ref_old_nan" if v is None or not np.isfinite(v) else
                   "nonpos_ref_old_positive" if v > 0 else "nonpos_ref_old_negative")
            npos[key] += 1
        return ro
    DS.band_local_deficit = counted
    try:
        DL.extract_and_merge(ns_run, eu_run, band_select="edge", out_dir=out_dir / "old")
    finally:
        DS.band_local_deficit = new
    res = {"counts": {**cnt, **npos}, "defect_calls": cnt["global_tail"] + cnt["nonmono_zone"]}
    if compare:
        DL.extract_and_merge(ns_run, eu_run, band_select="edge", out_dir=out_dir / "new")
        res["identical"] = {f: filecmp.cmp(out_dir / "old" / f, out_dir / "new" / f, shallow=False)
                            for f in ("delta_r_equiv.csv", "delta_r_next.csv")}
    return res


def test_v2(ns_run, eu_run, out_dir):
    print(f"\n--- V2 後方互換: {ns_run} × {eu_run} ---")
    res = extract_old_new(Path(ns_run), Path(eu_run), Path(out_dir))
    print(json.dumps(res, ensure_ascii=False), flush=True)
    c = res["counts"]
    if res["defect_calls"] == 0:
        check("V2 欠陥に当たる呼び出し 0 件 (ab_root_zone の計数) の抽出で新旧の delta_r_equiv.csv・delta_r_next.csv がバイト同一",
              all(res["identical"].values()),
              f"{res['identical']}、呼び出し {c['calls']}、探索域の参照が非正 {c['nonpos_ref_new']} 件 (旧が正の値を返したもの "
              f"{c['nonpos_ref_old_positive']}・負 {c['nonpos_ref_old_negative']}・NaN {c['nonpos_ref_old_nan']})")
    else:
        skip("V2 バイト同一", f"欠陥に当たる呼び出しが {res['defect_calls']} 件 (V2 の対象外; 差は記録のみ: {res['identical']})")


def test_case70(out_dir):
    print("\n--- case/70 run_0010 × run_0012 の抽出 vs ab_root_zone の zone の腕 ---")
    ref = NEXT_A_ZONE.parent
    if not (NS10.exists() and EU12.exists() and ref.exists()):
        skip("case/70 の抽出", "run_0010・run_0012・_ab_euler_steps/A/zone のどれかが無い")
        return
    out_dir = Path(out_dir)
    DL.extract_and_merge(NS10, EU12, band_select="edge", out_dir=out_dir)
    same = {f: filecmp.cmp(out_dir / f, ref / f, shallow=False) for f in ("delta_r_next.csv", "delta_r_equiv.csv", "delta_r_edge.csv")}
    check("case/70: 新実装の抽出が ab_root_zone の zone の腕とバイト同一", all(same.values()), f"{same}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--card", action="store_true")
    ap.add_argument("--v2", nargs=3, metavar=("NS_RUN", "EULER_RUN", "OUT"), default=None)
    ap.add_argument("--case70", default=None, metavar="OUT")
    a = ap.parse_args(argv)
    test_v1()
    test_v4_unit()
    with tempfile.TemporaryDirectory() as td:
        test_real(Path(td))
        if a.card:
            test_card(Path(td))
    if a.v2:
        test_v2(*a.v2)
    if a.case70:
        test_case70(a.case70)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
