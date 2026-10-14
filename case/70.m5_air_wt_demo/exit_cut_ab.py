"""出口の終点の切り落としの判別 A/B (plan tooling-nozzle-deltastar-band-root-zone §5.1 #4c、事前登録; forge 0 step)。

同じ逆 MOC の場 (同じ問題・同じ設計変数) で、cplus の閉包の物理出口 x_F での壁の表の切り方だけを変える:
  A = 今の切断 (`wall[wall[:, 0] <= x_F]`; 改修前の `moc_inverse._design_cplus` の式そのまま)
  B = 交点の追加 (`moc_inverse._cut_wall_at_exit`: 切った後、x_F が残った最後の点より 1e-9 r_t 以上下流なら、元の壁で x_F を
      挟む区間で全列を線形補間した点 (x_F, ...) を末尾に足す)
A は design_chain を改修前の切断で回す (`_cut_wall_at_exit` を差し替える)。B は A で解いた L_c を固定して (Lc_mode: explicit)
同じ場を作る (from_length の L_c は x_F だけで決まり切り方に依らないが、比較を同じ場に限るために固定する)。
参考に B を from_length のまま回した L_c も記録する (判定には使わない)。

判定 (事前登録; 結果を見て閾値を変えない):
  - 終点が落ちていた設計 (x_F − x_e,A ≥ 1e-9 r_t): 出口の移動 r_t·(x_e,B − x_e,A) が r_t·(x_F − x_e,A) と 0.05 mm 以内で一致
  - 落ちていない設計 (x_F − x_e,A < 1e-9 r_t): 逆 MOC の壁の表 (design_chain の wall_inv) が A と B でビット同一
対象: case/70 の problem_m5_air_c2b・final・c2・walltemp・ns2 と case/45 の生産 YAML (problem_d155_ns_prod.yaml; 凍結源の run が
この作業ツリーに無いので initial_line_run だけ絶対パスに置き換えた _monotone_gate_ab/_c45_prod_abs_initial_line.yaml を使い、
それが今の生産 YAML と initial_line_run 以外で一致することを確かめる)。

usage (case dir, design/.venv-opt):  python exit_cut_ab.py
出力: _exit_cut_ab/exit_cut_ab.json と標準出力の表
"""
import copy
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.geometry import moc_inverse as MI  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

OUT = C / "_exit_cut_ab"
OUT.mkdir(exist_ok=True)
DESIGNS = [("c70_c2b", C / "problem_m5_air_c2b.yaml"), ("c70_final", C / "problem_m5_air_final.yaml"),
           ("c70_c2", C / "problem_m5_air_c2.yaml"), ("c70_walltemp", C / "problem_m5_air_walltemp.yaml"),
           ("c70_ns2", C / "problem_m5_air_ns2.yaml")]
C45_PROD = C.parents[1] / "case/45.isobutane_m6_d155/problem_d155_ns_prod.yaml"
C45_ABS = C / "_monotone_gate_ab" / "_c45_prod_abs_initial_line.yaml"
GAP_MIN = 1e-9          # [r_t] これ未満は「落ちていない」 (plan §5.1 #4c)
TOL_MM = 0.05           # 出口の移動と r_t·(x_F − x_e,A) の一致の許容 [mm]

_NEW_CUT = MI._cut_wall_at_exit


def _old_cut(wall, x_F):
    """改修前の cplus の閉包の切断 (`wall = wall[wall[:, 0] <= x_F]`)。exit_info には何も足さない。"""
    return wall[wall[:, 0] <= x_F], {}


class _Capture:
    """design_chain の中の逆 MOC (`runner_axismach.inverse_design`) の結果を全部取る (読むだけ)。"""

    def __init__(self, cut):
        self.cut, self.calls = cut, []

    def __enter__(self):
        self._orig = RA.inverse_design
        MI._cut_wall_at_exit = self.cut

        def wrap(*a, **k):
            r = self._orig(*a, **k)
            self.calls.append((r, k))
            return r
        RA.inverse_design = wrap
        return self

    def __exit__(self, *exc):
        RA.inverse_design = self._orig
        MI._cut_wall_at_exit = _NEW_CUT
        return False


def _run(p, cut):
    """design_chain を cut の切り方で回す。戻り: (design_chain の結果 or None, 採用された逆 MOC の結果, 例外の文)。"""
    with _Capture(cut) as cap:
        try:
            d = RA.design_chain(p)
            err = None
        except Exception as ex:     # noqa: BLE001 — 下流 (壁 QA・当てはめ) の例外でも逆 MOC の表は比べる
            d, err = None, f"{type(ex).__name__}: {str(ex)[:300]}"
    if d is not None:
        xF = float(d["exit"]["x_F"])
        hit = [r for r, _ in cap.calls if float(r["exit"]["x_F"]) == xF]
        res = hit[-1] if hit else None
    else:
        res = cap.calls[-1][0] if cap.calls else None
    return d, res, err, cap.calls


def _explicit(p, L_c):
    q = copy.deepcopy(p)
    q.dv = dict(q.dv)
    q.dv["L_c"] = {"value": float(L_c)}
    q.geometry = dict(q.geometry)
    q.geometry["Lc_mode"] = "explicit"
    return q


def _c45_problem():
    """case/45 の生産 YAML。initial_line_run が相対パスで解ければそのまま、解けなければ前の担当の絶対パス版を使う
    (今の生産 YAML と initial_line_run 以外で一致することを確かめる)。"""
    raw = yaml.safe_load(C45_PROD.read_text())
    il = Path(str(raw["geometry"]["initial_line_run"]))
    if (C45_PROD.parent / il).is_dir():
        return load_problem(C45_PROD), {"source": str(C45_PROD.relative_to(C.parents[1])), "initial_line_run": str(il)}
    alt = yaml.safe_load(C45_ABS.read_text())
    a2, r2 = copy.deepcopy(alt), copy.deepcopy(raw)
    a2["geometry"].pop("initial_line_run")
    r2["geometry"].pop("initial_line_run")
    if a2 != r2:
        raise RuntimeError(f"{C45_ABS} が今の生産 YAML と initial_line_run 以外で一致しない")
    if not Path(alt["geometry"]["initial_line_run"]).is_dir():
        raise RuntimeError(f"凍結源 {alt['geometry']['initial_line_run']} が無い")
    return load_problem(C45_ABS), {"source": str(C45_PROD.relative_to(C.parents[1])),
                                   "used": str(C45_ABS.relative_to(C.parents[1])),
                                   "initial_line_run": alt["geometry"]["initial_line_run"],
                                   "equal_except_initial_line_run": True}


def _card_log_x(path: Path):
    """生成器の記録 (card_log.json の moc_field) の x_F・x_e (A の再現の確かめ用)。"""
    lp = Path(str(path) + ".card_log.json")
    if not lp.is_file():
        return None
    f = json.loads(lp.read_text()).get("moc_field") or {}
    return {"x_F": f.get("x_F"), "x_e": f.get("x_e")}


def one(name, p, src) -> dict:
    rt = float(p.spec["r_throat"])
    t0 = time.time()
    dA, resA, errA, callsA = _run(p, _old_cut)
    if resA is None:
        return {"name": name, "source": src, "r_t_m": rt, "error_A": errA}
    L_c = float(dA["L_c"]) if dA is not None else None
    if L_c is None:
        k = callsA[-1][1]
        L_c = float(k["x_E"]) - float(k["axis_anchor"][0])
    dB, resB, errB, _ = _run(_explicit(p, L_c) if p.geometry.get("Lc_mode", "explicit") != "explicit" else p, _NEW_CUT)
    # 参考: B を from_length のまま (L_c を固定しない) 回した L_c
    L_c_B_free = None
    if p.geometry.get("Lc_mode", "explicit") == "from_length":
        dBf, _, errBf, _ = _run(p, _NEW_CUT)
        L_c_B_free = float(dBf["L_c"]) if dBf is not None else errBf
    wA, wB = np.asarray(resA["wall"], dtype=float), np.asarray(resB["wall"], dtype=float)
    wfA, wfB = np.asarray(resA["wall_full"], dtype=float), np.asarray(resB["wall_full"], dtype=float)
    xF_A, xF_B = float(resA["exit"]["x_F"]), float(resB["exit"]["x_F"])
    same_field = bool(xF_A == xF_B and wfA.shape == wfB.shape and wfA.tobytes() == wfB.tobytes())
    x_eA, x_eB = float(wA[-1, 0]), float(wB[-1, 0])
    gap = xF_A - x_eA
    dropped = bool(gap >= GAP_MIN)
    rec_cut = W.moc_call_record(resA, callsA[-1][1])
    moved_mm = rt * (x_eB - x_eA) * 1e3
    expect_mm = rt * gap * 1e3
    bit_same = bool(wA.shape == wB.shape and wA.tobytes() == wB.tobytes())
    if dropped:
        ok = bool(abs(moved_mm - expect_mm) <= TOL_MM)
    else:
        ok = bit_same
    ex = resB["exit"]
    r = {"name": name, "source": src, "r_t_m": rt, "Lc_mode": p.geometry.get("Lc_mode", "explicit"), "L_c_A": L_c,
         "L_c_B_from_length": L_c_B_free, "same_field": same_field, "x_F": xF_A, "x_e_A": x_eA, "x_e_B": x_eB,
         "x_F_minus_x_e_A_rt": gap, "x_F_minus_x_e_A_mm": expect_mm, "moved_mm": moved_mm,
         "match_diff_mm": moved_mm - expect_mm, "dropped": dropped,
         "cut_dropped_node_at_x_F_A": rec_cut["cut_dropped_node_at_x_F"],
         "x_F_minus_nearest_node_rt": rec_cut["x_F_minus_nearest_node_rt"],
         "n_wall_A": int(len(wA)), "n_wall_B": int(len(wB)), "wall_table_bit_identical": bit_same,
         "prefix_identical": bool(len(wB) >= len(wA) and wB[:len(wA)].tobytes() == wA.tobytes()),
         "B_appended_row": (wB[-1].tolist() if len(wB) == len(wA) + 1 else None),
         "B_exit_info": {k: ex.get(k) for k in ("x_F_appended", "x_F_minus_last_node_rt")},
         "r_F": float(ex["r_F"]), "r_end_B_minus_r_F": float(wB[-1, 1] - float(ex["r_F"])),
         "design_x_e_A": (float(dA["wall"].x_e) if dA is not None else None),
         "design_x_e_B": (float(dB["wall"].x_e) if dB is not None else None),
         "error_A": errA, "error_B": errB, "card_log": _card_log_x(Path(src)) if Path(src).is_file() else None,
         "ok": ok, "time_s": time.time() - t0}
    if dA is not None and dB is not None and not dropped:
        xs = np.linspace(float(dA["wall"].x_in), float(dA["wall"].x_e), 20001)
        r["design_wall_bit_identical"] = bool(np.asarray(dA["wall"].r(xs)).tobytes() == np.asarray(dB["wall"].r(xs)).tobytes())
    return r


def _row(r):
    if "x_F" not in r:
        return f"{r['name']:<13} A でエラー: {r.get('error_A')}"
    jud = "ok" if r["ok"] else "NG"
    s = (f"{r['name']:<13} r_t {r['r_t_m']:.9f} | 落ちていた {'yes' if r['dropped'] else 'no ':<3} | x_F − x_e,A {r['x_F_minus_x_e_A_rt']:+.3e} r_t "
         f"= {r['x_F_minus_x_e_A_mm']:+.4f} mm | 移動 {r['moved_mm']:+.4f} mm | 一致の差 {r['match_diff_mm']:+.2e} mm | "
         f"表ビット同一 {r['wall_table_bit_identical']} | 点数 {r['n_wall_A']}→{r['n_wall_B']} | 同じ場 {r['same_field']} | {jud}")
    if r.get("error_A") or r.get("error_B"):
        s += f"\n   例外 A: {r.get('error_A')} / B: {r.get('error_B')}"
    return s


def main():
    out = {"plan": "plans/active/tooling-nozzle-deltastar-band-root-zone.md §5.1 #4c",
           "consult": "plans/accepted/tooling-nozzle-wt-deliverables.md §5.1 #15 (codex 諮問 2026-10-13 2 回目)",
           "A": "wall[wall[:, 0] <= x_F] (改修前の切断)",
           "B": "moc_inverse._cut_wall_at_exit (交点 (x_F, ...) を元の壁の区間で全列線形補間して追加; 間隔 < 1e-9 r_t は追加しない)",
           "rule": {"gap_min_rt": GAP_MIN, "tol_mm": TOL_MM,
                    "dropped": "x_F − x_e,A ≥ 1e-9 r_t → |r_t (x_e,B − x_e,A) − r_t (x_F − x_e,A)| ≤ 0.05 mm",
                    "not_dropped": "x_F − x_e,A < 1e-9 r_t → 逆 MOC の壁の表 (wall_inv) が A と B でビット同一"},
           "designs": []}
    for name, path in DESIGNS:
        r = one(name, load_problem(path), str(path))
        r["source"] = str(path.relative_to(C.parents[1]))
        out["designs"].append(r)
        print(_row(r), flush=True)
    p45, info45 = _c45_problem()
    r = one("c45_prod", p45, str(C45_PROD))
    r["source"] = str(C45_PROD.relative_to(C.parents[1]))
    r["c45_input"] = info45
    out["designs"].append(r)
    print(_row(r), flush=True)
    ok = all(d.get("ok") for d in out["designs"])
    out["judgement"] = {"pass": ok, "n": len(out["designs"]),
                        "n_dropped": sum(1 for d in out["designs"] if d.get("dropped")),
                        "failed": [d["name"] for d in out["designs"] if not d.get("ok")]}
    print("\n判定:", json.dumps(out["judgement"], ensure_ascii=False))
    (OUT / "exit_cut_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
