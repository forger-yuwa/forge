#!/usr/bin/env python3
"""要求カードの M_knot の自動の流れ (規則 coarse_min_theta_v2) の試験 — plan tooling-nozzle-mknot-fast §6 V2 と §4.2 の並列 (§5.1 #3)。

選び方そのもの (`select_coarse_min_theta`) は本物を回し、評価 (`make_eval_fn`) を合成の関数に差し替えてカードの生成 (`build_problem`)
を通す。合成の評価は選択ごとに谷の中心 c を決め、粗い網では [c − 0.3, c + 0.3] で成立・θ_max = 15 + 10 (M_K − c)²、生産の網では c だけ
本物の evaluate_mknot (最終の設計チェーンの θ_max との照合のため; n_axis_inv 600 で 1 回 約 1 s)、他は合成。r_t の固定点の
`solve_rt_throat` は合成 (r_t* = 0.0249 + 1e-4 (M_K − 2)、1 回 1000 s を合成の時計に足す; 停止・場の不足の例外を差し込める)、
最後の設計チェーンは本物。カードは M5 の例 (examples/m5_air_demo.yaml) に n_axis_inv 600・長さ公差 2 mm (L_total_tol 0.04 r_t) を
当てたもの (手元で軽く回すため)。

V2:
  (a) 最初と最終の r_t で選び直す (物理スロート径): 最初の r_t (径/2) と固定点の後の r_t で 1 回ずつ、一致して終了
  (b) M_K が変われば固定点をやり直す: 2 回目のやり直しで一致して成功 / 3 回目が要るなら例外 (YAML を書かない)
  (c) M_K の往復 (2.10 → 2.20 → 2.10): 同じ実効の入力の選択を使い回し、やり直しの上限で例外
  (d) 自動調整のそれぞれで評価を捨てて最初からやり直す: n_axis_inv の倍化・L_total_tol の緩和 (固定点の from_length の停止)、
      x_end_margin の変更 (固定点の場の不足、最終の壁の判定の後の張り付き)
  (e) 最終の 3 点の評価・YAML・card_log の入力 (r_t・n_axis_inv・L_total_tol・x_end_margin・M_K) が対応する (全部の正例で)
  (f) 合成の時計: 選び直しによる固定点のやり直しに既知の時間 (1000 s × 回数) を与え、M_knot_rule.elapsed_s にその全量が入る
      (選択の時間と合わせて一致)、カード全体の時間も並べて記録
  ほか: 設計スロート径 (固定点なし: 最初の選択がそのまま最終)、本物の r_t の固定点で 1 回 (--slow でなくても回す、約 1 分)、
  記録 (粗い掃引の全行・谷ごとの確認の行・選択の理由・N_s・L_s・N_p、YAML の要約の行 ≤ 4 kB、mknot_card_log_check、Euler の写し)、
  並列の数 (mknot_workers: 上限 4・CPU 数・0.7 × MemAvailable / ピーク、FORGE_MKNOT_WORKERS、取れなければ 1、算出 0 は評価の前に例外)。
選択の中の長さ合わせの停止 (plan tooling-nozzle-card-local-mesh §5.1 #6d; plan tooling-nozzle-mknot-fast §4.1 の拡張):
  (6d-a) 谷 1 つで中心が生産の網で未収束 → n_axis_inv を倍にして粗い掃引から選び直して選べる、(6d-b) 谷 2 つで片方だけ未収束 → 発火する
  (今までは黙って他方)、(6d-c) 制約による不成立だけ → 発火しない、(6d-d) 粗い網が全点で未収束 → 発火する、(6d-e) 倍化でも止まり緩和で
  通る (on_stall と同じ順)、(6d-f) 使い切ったら例外。

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_mknot_v2_card_tests.py
所要: 約 3〜5 分 (手元、1 プロセス; FORGE_MKNOT_WORKERS=1 を立てて回す)。
"""
import contextlib
import copy
import json
import os
import sys
import tempfile
import time
from pathlib import Path

import yaml

os.environ["FORGE_MKNOT_WORKERS"] = "1"          # 手元の PC: 評価は 1 プロセスずつ (2026-10-14 の事故の再発防止)

DESIGN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate import mknot_select as MS  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

FAIL = 0
REAL_EVAL = MS.evaluate_mknot
REAL_CHAIN_FIELD = W._chain_field
M5_CARD = DESIGN / "forge_design" / "intake" / "examples" / "m5_air_demo.yaml"
BASE = yaml.safe_load(M5_CARD.read_text())
YAML_LINE_MAX = 4096


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""), flush=True)
    if not cond:
        FAIL += 1


@contextlib.contextmanager
def patched(obj, **kw):
    old = {k: getattr(obj, k) for k in kw}
    for k, v in kw.items():
        setattr(obj, k, v)
    try:
        yield
    finally:
        for k, v in old.items():
            setattr(obj, k, v)


def card(name, **kw):
    """M5 の例のカード (自動の M_knot) に n_axis_inv 600・長さ公差 2 mm を当てた写し (kw で上書き)。"""
    c = copy.deepcopy(BASE)
    c.update(name=name, overrides={"geometry.n_axis_inv": 600}, tolerances={"length_m": 0.002, "diameter_m": 2.0e-5})
    c.update(kw)
    return c


def rstar(M):
    """合成の r_t の固定点の解 (M_K ごとに違う; 同じ M_K なら同じ倍精度)。"""
    return 0.0249 + 1e-4 * (float(M) - 2.0)


class Harness:
    """合成の評価・合成の solve_rt_throat・合成の時計・最終の壁の張り付きの差し込み。
    centers: 選択 (粗い網の評価の関数を作る呼び出し) ごとの谷の中心 (足りなければ最後の値)。stall_on / short_on: solve_rt_throat の
    何回目で from_length の停止 / 場の不足の例外を出すか。pin_on: 最終の設計チェーンの何回目を「場の終点に張り付き」にするか。"""

    def __init__(self, centers, *, stall_on=(), short_on=(), pin_on=(), fake_solve=True, extra=None, unconv=None, coarse_unconv=None,
                 bad=None):
        self.centers, self.stall_on, self.short_on, self.pin_on = list(centers), set(stall_on), set(short_on), set(pin_on)
        self.fake_solve = fake_solve
        # §5.1 #6d の試験: extra = 2 つ目の谷 (中心, θ の底) (合成だけ)。unconv(選択の番号, M_K, N_p, L_p) が真なら生産の網の行を
        # 長さ合わせの未収束 (from_length 未収束) に、coarse_unconv(選択の番号, M_K, N_s) が真なら粗い網の行を同じく未収束に、
        # bad(選択の番号, M_K, N_p) が真なら生産の網の行を制約で不成立 (壁での余裕 0.5°) にする
        self.extra, self.unconv, self.coarse_unconv, self.bad = extra, unconv, coarse_unconv, bad
        self.clock = 0.0
        self.sel: list = []
        self.solves: list = []
        self.n_field = 0

    def now(self):
        return self.clock

    def c_of(self, i):
        return self.centers[min(i, len(self.centers) - 1)]

    def _row(self, M, c, rt, N, L, xem):
        M = float(M)
        ex = self.extra
        if abs(M - c) <= 0.3 + 1e-9:
            r = {"M_K": M, "theta_max_deg": 15.0 + 10.0 * (M - c) ** 2, "margin_deg": 4.0}
        elif ex is not None and abs(M - ex[0]) <= 0.3 + 1e-9:
            r = {"M_K": M, "theta_max_deg": ex[1] + 10.0 * (M - ex[0]) ** 2, "margin_deg": 4.0}
        else:
            r = {"M_K": M, "theta_max_deg": 14.0, "margin_deg": 0.5}
        r.update(min_sin_angle=0.05, mpp_sign_changes=1, mpp_quality_ok=True, n_flip_interior=0, n_crossings_sampled=0,
                 pileup_pairs_interior=0, qa_violations=[], moc_gate_ok=True, moc_gate_reason=None, n_axis=int(N), r_t_m=rt,
                 L_total_tol=float(L), x_end_margin=xem, L_c=21.5, L_residual_rt=0.0, length_stop="合成")
        r["reasons"] = MS.feasibility_reasons(r)
        r["ok"] = not r["reasons"]
        return r

    def make_eval_fn(self, base, n_axis_inv=None, base_dir=None, L_total_tol=None):
        g = base["geometry"]
        rt = float(base["spec"]["r_throat"])
        xem = float(g.get("x_end_margin", 2.3))
        if L_total_tol is not None:                      # 粗い網 = 選択の始まり
            S = {"idx": len(self.sel), "r_t": rt, "N_s": int(n_axis_inv), "L_s": float(L_total_tol), "xem": xem,
                 "c": self.c_of(len(self.sel)), "n_eval": 0, "N_p": None, "L_p": None}
            self.sel.append(S)

            def coarse(M):
                self.clock += 1.0
                S["n_eval"] += 1
                if self.coarse_unconv is not None and self.coarse_unconv(S["idx"], float(M), int(n_axis_inv)):
                    return self._unconv_row(M, rt, n_axis_inv, L_total_tol, xem)
                return self._row(M, S["c"], rt, n_axis_inv, L_total_tol, xem)
            return coarse
        S = self.sel[-1]
        S.update(N_p=int(n_axis_inv), L_p=float(g["L_total_tol"]))

        def prod(M):
            self.clock += 1.0
            S["n_eval"] += 1
            if self.unconv is not None and self.unconv(S["idx"], float(M), int(n_axis_inv), S["L_p"]):
                return self._unconv_row(M, rt, n_axis_inv, S["L_p"], xem)
            if self.bad is not None and self.bad(S["idx"], float(M), int(n_axis_inv)):
                r = self._row(M, S["c"], rt, n_axis_inv, S["L_p"], xem)
                r.update(margin_deg=0.5)
                r["reasons"] = MS.feasibility_reasons(r)
                r["ok"] = not r["reasons"]
                return r
            if abs(float(M) - S["c"]) < 1e-9:
                return REAL_EVAL(base, float(M), n_axis_inv)     # 選ぶ点は本物 (最終の設計チェーンの θ_max と照合される)
            return self._row(M, S["c"], rt, n_axis_inv, S["L_p"], xem)
        return prod

    @staticmethod
    def _unconv_row(M, rt, N, L, xem):
        """長さ合わせが許容に届かずに止まった行 (evaluate_mknot が design_chain の from_length の例外を捕まえた形)。"""
        r = {"M_K": float(M), "error": f"ValueError: Lc_mode: from_length 未収束 (|x_F − L_total| = 0.0172 > tol {float(L):.3g}) — 合成",
             "n_axis": None, "n_axis_input": int(N), "r_t_m": rt, "L_total_tol": float(L), "x_end_margin": xem, "L_residual_rt": -0.0172,
             "length_stop": "from_length 未収束 (改善停止か反復上限; 最良の |x_F − L_total| > L_total_tol)"}
        r["reasons"] = MS.feasibility_reasons(r)
        r["ok"] = not r["reasons"]
        return r

    def solve(self, problem, R_throat_m, prev_run=None, delta_r_out=None, delta_next=None, max_iter=None, tol_R_m=None,
              _delta_supplier=None):
        k = len(self.solves) + 1
        self.clock += 1000.0
        d = yaml.safe_load(Path(problem).read_text())
        g = d["geometry"]
        self.solves.append({"k": k, "r_in": float(d["spec"]["r_throat"]), "M_K": float(g["M_knot"]), "n_axis_inv": g.get("n_axis_inv"),
                            "L_total_tol": g.get("L_total_tol"), "x_end_margin": g.get("x_end_margin")})
        if k in self.stall_on:
            raise ValueError("Lc_mode: from_length 未収束 (|x_F − L_total| = 0.1 > tol 0.04) — 合成")
        if k in self.short_on:
            raise RuntimeError("終端特性線が壁に到達しない (合成)")
        return {"r_t_m": rstar(g["M_knot"]), "residual_m": 0.0, "n_iter": 1,
                "physical_throat": {"r_m": float(R_throat_m), "x_m": 0.0}, "exit_radius_m": 0.1}

    def chain_field(self, chain, calls):
        self.n_field += 1
        f = REAL_CHAIN_FIELD(chain, calls)
        return dict(f, pinned=True) if self.n_field in self.pin_on else f


def run(h, card_dict, name):
    """合成の評価で build_problem を回す。戻り (card_log, YAML の文字列, mknot_card_log_check の結果, 例外の文字列, 所要 s)。"""
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        cp = Path(td) / f"{name}.yaml"
        cp.write_text(yaml.safe_dump(card_dict, sort_keys=False, allow_unicode=True))
        out = Path(td) / f"problem_{name}.yaml"
        solve = h.solve if h.fake_solve else DL.solve_rt_throat
        with patched(MS, make_eval_fn=h.make_eval_fn), patched(DL, solve_rt_throat=solve), \
                patched(W, _clock=h.now, _chain_field=h.chain_field):
            try:
                logd = W.build_problem(cp, out, log=lambda *a, **k: None)
            except Exception as ex:          # noqa: BLE001 (例外の試験)
                return None, None, None, f"{type(ex).__name__}: {ex}", time.time() - t0
            text = out.read_text()
            chk = W.mknot_card_log_check(out)
            # Euler の写し: この雛形 (explore) では書かれないので、物理壁のキーの行を 1 行足した YAML から write_euler_companion で作って比べる
            lines = text.splitlines(keepends=True)
            k = next(i for i, ln in enumerate(lines) if ln.startswith("geometry:"))
            yp = Path(td) / "with_pw.yaml"
            yp.write_text("".join(lines[:k + 1] + ["  physical_wall_repr: single_bspline\n"] + lines[k + 1:]))
            eo = W.write_euler_companion(yp)             # 同じディレクトリなので要約の card_log の名前のまま照合できる
            logd["_test_euler_line"] = [ln for ln in eo.read_text().splitlines() if ln.startswith("  M_knot_rule:")]
            logd["_test_euler_check"] = W.mknot_card_log_check(eo)
    return logd, text, chk, None, time.time() - t0


def summary_line(text):
    ln = [x for x in text.splitlines() if x.startswith("  M_knot_rule:")]
    return ln[0] if len(ln) == 1 else None


def verify_output(tag, h, logd, text, chk):
    """(e) 出力の組と最終の選択の対応・記録の形。"""
    raw = yaml.safe_load(text)
    g, sp = raw["geometry"], raw["spec"]
    full, summ = logd["M_knot_rule"], raw["_card"]["M_knot_rule"]
    M = float(g["M_knot"])
    k = MS.v2_k(M)
    want = {"r_t_m": float(sp["r_throat"]), "n_axis": int(g["n_axis_inv"]), "L_total_tol": float(g["L_total_tol"]),
            "x_end_margin": float(g.get("x_end_margin", 2.3))}
    rows3 = full["final_three"]
    ok3 = ([r["M_K"] for r in rows3] == [MS.v2_M(k - 1), M, MS.v2_M(k + 1)] and all(r["ok"] for r in rows3)
           and all(all(r.get(kk) == vv for kk, vv in want.items()) for r in rows3))
    inp = full["selection"]["inputs"]
    ok_inp = (inp == {"r_t_m": want["r_t_m"], "N_s": want["n_axis"] // 4, "L_s": 0.1, "N_p": want["n_axis"], "L_p": want["L_total_tol"],
                      "x_end_margin": want["x_end_margin"]})
    last = h.sel[-1]
    ok_h = (last["r_t"] == want["r_t_m"] and last["N_p"] == want["n_axis"] and last["L_p"] == want["L_total_tol"]
            and last["xem"] == want["x_end_margin"] and last["c"] == M)
    check(f"{tag} (e): 最終の 3 点 ({MS.v2_M(k - 1)}・{M}・{MS.v2_M(k + 1)}) の評価の入力 = YAML の r_t・n_axis_inv・L_total_tol・x_end_margin、"
          "選択の入力 (N_s = N_p/4・L_s 0.1) と合成の評価の呼び出しも同じ、全部成立",
          ok3 and ok_inp and ok_h and summ["ok_final_three"] == [[r["M_K"], True] for r in rows3],
          f"YAML {want} 入力 {inp} 最後の選択 {({k_: last[k_] for k_ in ('r_t', 'N_p', 'L_p', 'xem', 'c')})} "
          f"行 {[(r['M_K'], r.get('r_t_m'), r.get('n_axis'), r.get('L_total_tol'), r.get('x_end_margin')) for r in rows3]}")
    check(f"{tag} (e): YAML の geometry.M_knot = card_log の記録の chosen = 要約の chosen、最終の設計チェーンの θ_max = 選んだ点の行、"
          "mknot_card_log_check で照合 ok、要約の行 ≤ 4 kB",
          M == full["chosen"] == summ["chosen"] and full["theta_max_deg_final_chain"] == full["theta_max_deg"] and chk["ok"] is True
          and summary_line(text) is not None and len(summary_line(text).encode()) <= YAML_LINE_MAX and summ["rule"] == MS.RULE_NAME_V2,
          f"{M} / {full['chosen']} / {summ['chosen']} 照合 {chk.get('ok')} 行 {len((summary_line(text) or '').encode())} byte")
    return raw, full, summ


# ---------------------------------------------------------------------------------------------------------------------
print("== 並列の数 (§4.2: mknot_workers)")
os.environ.pop("FORGE_MKNOT_WORKERS", None)
with patched(W, _mem_available_gb=lambda: 10.0, _cpu_count=lambda: 12):
    w = W.mknot_workers(600, 2400)
check("使えるメモリ 10 GB・CPU 12: 粗い網 min(12, 4, ⌊7/0.27⌋ = 25) = 4、生産の網 min(12, 4, ⌊7/1.31⌋ = 5) = 4",
      w["coarse"] == 4 and w["prod"] == 4 and w["by_mem_coarse"] == 25 and w["by_mem_prod"] == 5, str(w))
with patched(W, _mem_available_gb=lambda: 10.0, _cpu_count=lambda: 2):
    w = W.mknot_workers(600, 2400)
check("CPU 2 なら 2・2", w["coarse"] == 2 and w["prod"] == 2, str({k: w[k] for k in ("coarse", "prod")}))
with patched(W, _mem_available_gb=lambda: 4.0, _cpu_count=lambda: 12):
    w = W.mknot_workers(600, 2400)
check("使えるメモリ 4 GB: 生産の網は ⌊2.8/1.31⌋ = 2、粗い網は 4", w["coarse"] == 4 and w["prod"] == 2, str({k: w[k] for k in ("coarse", "prod")}))
os.environ["FORGE_MKNOT_WORKERS"] = "1"
with patched(W, _mem_available_gb=lambda: 10.0, _cpu_count=lambda: 12):
    w = W.mknot_workers(600, 2400)
check("FORGE_MKNOT_WORKERS=1 で上限を 1 に下げる", w["coarse"] == 1 and w["prod"] == 1 and w["cap_env"] == 1, str(w["cap"]))
for bad_env in ("0", "x", "-2"):
    os.environ["FORGE_MKNOT_WORKERS"] = bad_env
    try:
        W.mknot_workers(600, 2400)
        check(f"FORGE_MKNOT_WORKERS={bad_env!r} → ValueError", False)
    except ValueError as e:
        check(f"FORGE_MKNOT_WORKERS={bad_env!r} → ValueError", "FORGE_MKNOT_WORKERS" in str(e), str(e)[:100])
os.environ["FORGE_MKNOT_WORKERS"] = "1"
with patched(W, _mem_available_gb=lambda: None, _cpu_count=lambda: 12):
    os.environ.pop("FORGE_MKNOT_WORKERS")
    w = W.mknot_workers(600, 2400)
    os.environ["FORGE_MKNOT_WORKERS"] = "1"
check("MemAvailable が取れなければ 1・1 (上限 4 でも)", w["coarse"] == 1 and w["prod"] == 1 and w["available_gb"] is None, str(w))
with patched(W, _mem_available_gb=lambda: 1.0, _cpu_count=lambda: 12):
    try:
        W.mknot_workers(600, 2400)
        check("使えるメモリ 1 GB: 生産の網の算出 ⌊0.7/1.31⌋ = 0 → 評価の前に RuntimeError", False)
    except RuntimeError as e:
        check("使えるメモリ 1 GB: 生産の網の算出 ⌊0.7/1.31⌋ = 0 → 評価の前に RuntimeError", "算出が 0" in str(e) and "生産の網" in str(e),
              str(e)[:160])
pk = [W.mknot_eval_peak_gb(n) for n in (300, 600, 1200, 2400, 4800)]
check("1 評価のピークの見積もり: 600 以下は 0.27、2400 は 1.31、間は単調、4800 は 1.31 × 2^1.6",
      pk[0] == pk[1] == 0.27 and pk[3] == 1.31 and pk[1] < pk[2] < pk[3] and abs(pk[4] - 1.31 * 2 ** 1.6) < 1e-12,
      str([round(x, 4) for x in pk]))

CARD_PHYS = card("m5_v2_phys")
CARD_DESIGN = card("m5_v2_design", throat={"diameter_design_m": 0.050})
# x_end_margin の試験は長さ公差 4 mm (L_total_tol 0.08 r_t): 2 mm (0.04 r_t) では n_axis_inv 600・x_end_margin 2.8・r_t 0.02491 m の
# from_length が止まる (合成でない本物の停止; 2026-10-15 の手元の確認) — 試験の狙いは x_end_margin の変更なので緩める
CARD_PHYS_LOOSE = card("m5_v2_phys_loose", tolerances={"length_m": 0.004, "diameter_m": 2.0e-5})

print("== 算出 0 なら評価の前に例外 (カード)")
h = Harness([2.10])
with patched(W, _mem_available_gb=lambda: 0.1):
    os.environ.pop("FORGE_MKNOT_WORKERS")
    logd, text, chk, err, dt = run(h, CARD_PHYS, "w0")
    os.environ["FORGE_MKNOT_WORKERS"] = "1"
check("使えるメモリ 0.1 GB: build_problem は評価の前に RuntimeError (評価の関数を 1 度も作らない・固定点も回さない)",
      err is not None and "算出が 0" in err and h.sel == [] and h.solves == [], (err or "")[:200])

print("== (a) 最初と最終の r_t で選び直す (物理スロート径、一致して終了)")
h = Harness([2.10])
logd, text, chk, err, dt = run(h, CARD_PHYS, "a")
if err is None:
    raw, full, summ = verify_output("(a)", h, logd, text, chk)
    check("(a) 選択 2 回: 最初の r_t 0.025 (径/2) と固定点の後の r_t (= 出力の r_t)、一致 → やり直しなし",
          [s["r_t"] for s in h.sel] == [0.025, rstar(2.10)] and raw["spec"]["r_throat"] == rstar(2.10) and full["redo"] == 0
          and full["first_reused"] is False and full["first_r_t_m"] == 0.025 and full["r_t_m"] == rstar(2.10)
          and [e["same"] for e in full["events"]] == [None, True] and len(h.solves) == 2 and raw["geometry"]["M_knot"] == 2.1,
          f"選択の r_t {[s['r_t'] for s in h.sel]} events {[e['same'] for e in full['events']]}")
    n_ev = sum(s["n_eval"] for s in h.sel)
    check("(a) (f) 計時: elapsed_s = 選択の評価の時間 (合成 1 s × 評価の数) だけ (最初の固定点は含めない)、カード全体 = それ + 固定点 2000 s",
          full["elapsed_s"] == float(n_ev) and full["elapsed_breakdown"] == {"select_s": float(n_ev), "refix_s": 0.0}
          and full["card_elapsed_s"] == float(n_ev) + 2000.0 and logd["time_s"] == full["card_elapsed_s"],
          f"elapsed {full['elapsed_s']} 評価 {n_ev} カード {full['card_elapsed_s']} time_s {logd['time_s']}")
    sel_rec = full["selection"]
    grid_n = sel_rec["grid"][2]
    check("(a) 記録: card_log に粗い掃引の全行 (L_c・符号つきの残差・停止の理由つき)・谷ごとの確認の行・選択の理由・N_s・L_s・N_p、"
          "YAML の要約は表を持たない",
          len(sel_rec["coarse_table"]) == grid_n == h.sel[-1]["n_eval"] - 3
          and all(all(kk in r for kk in ("L_c", "L_residual_rt", "length_stop", "n_axis")) for r in sel_rec["coarse_table"])
          and all(r["n_axis"] == 150 for r in sel_rec["coarse_table"]) and len(sel_rec["confirm_table"]) == 3
          and sel_rec["valleys"][0]["steps"][0]["action"].startswith("中心 2.10 と両隣が成立") and "reason" in full
          and (full["N_s"], full["L_s"], full["N_p"]) == (150, 0.1, 600) and (summ["N_s"], summ["L_s"], summ["N_p"]) == (150, 0.1, 600)
          and not ({"selection", "first_selection", "final_three", "timing", "coarse_table"} & set(summ))
          and "length_stop" not in json.dumps(summ, ensure_ascii=False),
          f"粗い掃引 {len(sel_rec['coarse_table'])} 行 / 格子 {grid_n}、確認 {len(sel_rec['confirm_table'])} 行")
    check("(a) 記録: YAML のコメント (geometry.M_knot・見出し) と対応表の行が規則 coarse_min_theta_v2 の文",
          f"# M_knot: 自動 (規則 {MS.RULE_NAME_V2}" in text and f"  M_knot: 2.1   # 自動 (規則 {MS.RULE_NAME_V2})" in text
          and any(r["card_key"] == "M_knot (自動)" and MS.RULE_NAME_V2 in r["rule"] for r in logd["mapping"])
          and W._mknot_rule_text(summ) == W._mknot_rule_text(full), "")
    eu = logd.get("_test_euler_line")
    check("(a) Euler の写しも同じ要約の行で、card_log と照合できる", eu == [summary_line(text)] and logd["_test_euler_check"]["ok"] is True,
          str(logd["_test_euler_check"]))
else:
    check("(a) build_problem が通る", False, err[:300])
print(f"   (a) {dt:.0f} s")

print("== (b) M_K が変われば固定点をやり直す (2 回目で成功・3 回目で例外)")
h = Harness([2.10, 2.20, 2.15])
logd, text, chk, err, dt = run(h, CARD_PHYS, "b_ok")
if err is None:
    raw, full, summ = verify_output("(b)", h, logd, text, chk)
    hist = raw["_card"]["r_t_history"]
    check("(b) 2.10 → 2.20 → 2.15 → 2.15: やり直し 2 回目で一致して 2.15 を出力 (始まりは前の r_t)",
          raw["geometry"]["M_knot"] == 2.15 and full["redo"] == 2 and [s["c"] for s in h.sel] == [2.1, 2.2, 2.15, 2.15]
          and [s["r_t"] for s in h.sel] == [0.025, rstar(2.10), rstar(2.20), rstar(2.15)]
          and [x["M_K"] for x in h.solves] == [2.1, 2.1, 2.2, 2.2, 2.15, 2.15]
          and [x["r_in"] for x in h.solves][2] == rstar(2.10) and [hh["M_K"] for hh in hist] == [2.1, 2.1, 2.2, 2.2, 2.15, 2.15]
          and [e["same"] for e in full["events"]] == [None, False, False, True],
          f"選択 {[s['c'] for s in h.sel]} 固定点 {[x['M_K'] for x in h.solves]} events {[e['same'] for e in full['events']]}")
    n_ev = sum(s["n_eval"] for s in h.sel)
    check("(b) (f) 合成の時計: 固定点のやり直し 4 回 (2 回 × 2) = 4000 s が elapsed_s に全部入る (= 選択の評価の時間 + 4000)、"
          "最初の固定点 2000 s は入らない、カード全体 = 評価 + 6000",
          full["elapsed_s"] == float(n_ev) + 4000.0 and full["elapsed_breakdown"] == {"select_s": float(n_ev), "refix_s": 4000.0}
          and full["card_elapsed_s"] == float(n_ev) + 6000.0 and logd["time_s"] == full["card_elapsed_s"]
          and summ["elapsed_s"] == full["elapsed_s"] and sum(t["elapsed_s"] for t in full["timing"]) == full["elapsed_s"],
          f"elapsed {full['elapsed_s']} 内訳 {full['elapsed_breakdown']} 評価 {n_ev} カード {full['card_elapsed_s']}")
else:
    check("(b) build_problem が通る", False, err[:300])
print(f"   (b ok) {dt:.0f} s")
h = Harness([2.10, 2.20, 2.15, 2.25])
logd, text, chk, err, dt = run(h, CARD_PHYS, "b_ng")
check("(b) 2.10 → 2.20 → 2.15 → 2.25: 3 回目のやり直しが要る → RuntimeError (YAML を書かない)",
      err is not None and err.startswith("RuntimeError") and "2 回やり直しても一致しない" in err and len(h.sel) == 4 and len(h.solves) == 6,
      (err or "")[:240])

print("== (c) M_K の往復")
h = Harness([2.10, 2.20, 2.10])
logd, text, chk, err, dt = run(h, CARD_PHYS, "c")
check("(c) 2.10 → 2.20 → 2.10: 3 回目の固定点の後の r_t は 1 回目と同じ → 同じ実効の入力の選択 (2.20) を使い回し、上限で RuntimeError",
      err is not None and "2 回やり直しても一致しない" in err and len(h.sel) == 3 and [x["M_K"] for x in h.solves] == [2.1, 2.1, 2.2, 2.2, 2.1, 2.1]
      and h.solves[-1]["r_in"] == rstar(2.10), f"選択 {len(h.sel)} 回 {(err or '')[:160]}")

print("== (d) 自動調整で評価を捨てて最初からやり直す")
h = Harness([2.10], stall_on={1})
logd, text, chk, err, dt = run(h, CARD_PHYS, "d_naxis")
if err is None:
    raw, full, summ = verify_output("(d) n_axis_inv", h, logd, text, chk)
    fs = full["first_selection"]
    check("(d) n_axis_inv の倍化: 最初の固定点の from_length の停止 → n_axis_inv 600 → 1200、評価を捨てて最初の r_t から選び直し "
          "(粗い網 300・生産の網 1200 で全部を評価し直す)",
          raw["geometry"]["n_axis_inv"] == 1200 and full["n_restart"] == 1 and len(full["discarded"]) == 1
          and "geometry.n_axis_inv 600 → 1200" in full["discarded"][0]["reason"]
          and [x["discarded"] for x in full["selections"]] == [True, False, False]
          and [(s["N_s"], s["N_p"], s["r_t"]) for s in h.sel] == [(150, 600, 0.025), (300, 1200, 0.025), (300, 1200, rstar(2.10))]
          and h.sel[1]["n_eval"] == h.sel[0]["n_eval"] and fs is not None and all(r["n_axis"] == 300 for r in fs["coarse_table"])
          and full["first_r_t_m"] == 0.025 and h.solves[1]["n_axis_inv"] == 1200,
          f"選択 {[(s['N_s'], s['N_p'], s['r_t']) for s in h.sel]} 捨てた {full['discarded']}")
else:
    check("(d) n_axis_inv の倍化: build_problem が通る", False, err[:300])
print(f"   (d n_axis_inv) {dt:.0f} s")
h = Harness([2.10], stall_on={1, 2})
logd, text, chk, err, dt = run(h, CARD_PHYS, "d_ltol")
if err is None:
    raw, full, summ = verify_output("(d) L_total_tol", h, logd, text, chk)
    check("(d) L_total_tol の緩和: 2 回目の停止で 0.04 → 0.05 r_t (n_axis_inv は 600 に戻す)、そのたびに評価を捨てて選び直す",
          raw["geometry"]["L_total_tol"] == 0.05 and raw["geometry"]["n_axis_inv"] == 600 and full["n_restart"] == 2
          and raw["_card"]["L_total_tol_relaxed"] is not None and "geometry.L_total_tol" in full["discarded"][1]["reason"]
          and [(s["N_p"], s["L_p"]) for s in h.sel] == [(600, 0.04), (1200, 0.04), (600, 0.05), (600, 0.05)]
          and [x["discarded"] for x in full["selections"]] == [True, True, False, False] and full["L_p"] == 0.05,
          f"選択 {[(s['N_p'], s['L_p']) for s in h.sel]} 捨てた {[d['reason'][:60] for d in full['discarded']]}")
else:
    check("(d) L_total_tol の緩和: build_problem が通る", False, err[:300])
print(f"   (d L_total_tol) {dt:.0f} s")
h = Harness([2.10], short_on={1})
logd, text, chk, err, dt = run(h, CARD_PHYS_LOOSE, "d_short")
if err is None:
    raw, full, summ = verify_output("(d) x_end_margin (場の不足)", h, logd, text, chk)
    check("(d) x_end_margin の変更 (固定点で終端特性線が壁に届かない): 2.3 → 2.8、評価を捨てて選び直す",
          raw["geometry"]["x_end_margin"] == 2.8 and full["n_restart"] == 1 and "x_end_margin" in full["discarded"][0]["reason"]
          and [s["xem"] for s in h.sel] == [2.3, 2.8, 2.8] and full["x_end_margin"] == 2.8,
          f"選択の x_end_margin {[s['xem'] for s in h.sel]}")
else:
    check("(d) x_end_margin (場の不足): build_problem が通る", False, err[:300])
print(f"   (d x_end_margin 場の不足) {dt:.0f} s")
h = Harness([2.10], pin_on={1})
logd, text, chk, err, dt = run(h, CARD_PHYS_LOOSE, "d_pin")
if err is None:
    raw, full, summ = verify_output("(d) x_end_margin (最終の壁の判定の後)", h, logd, text, chk)
    check("(d) x_end_margin の変更 (最終の壁の判定の後の張り付き): 一致して出力する直前に 2.3 → 2.8、評価を捨てて最初の r_t から選び直す",
          raw["geometry"]["x_end_margin"] == 2.8 and full["n_restart"] == 1 and h.n_field == 2
          and [(s["xem"], s["r_t"]) for s in h.sel] == [(2.3, 0.025), (2.3, rstar(2.10)), (2.8, 0.025), (2.8, rstar(2.10))]
          and [x["discarded"] for x in full["selections"]] == [True, True, False, False],
          f"選択 {[(s['xem'], s['r_t']) for s in h.sel]}")
else:
    check("(d) x_end_margin (最終の壁の判定の後): build_problem が通る", False, err[:300])
print(f"   (d x_end_margin 張り付き) {dt:.0f} s")

print("== 設計スロート径 (固定点なし: 最初の選択がそのまま最終)")
h = Harness([2.10])
logd, text, chk, err, dt = run(h, CARD_DESIGN, "design")
if err is None:
    raw, full, summ = verify_output("設計スロート径", h, logd, text, chk)
    n_ev = sum(s["n_eval"] for s in h.sel)
    check("設計スロート径: 選択 1 回 (r_t 0.025)、最終は同じ実効の入力の選択を使い回す、固定点なし、elapsed_s = 評価の時間",
          len(h.sel) == 1 and full["first_reused"] is True and full["first_selection"] is None and h.solves == []
          and full["elapsed_s"] == float(n_ev) and "使い回した" in full["events"][-1]["stage"], f"選択 {len(h.sel)}")
else:
    check("設計スロート径: build_problem が通る", False, err[:300])
print(f"   (設計スロート径) {dt:.0f} s")

print("== 本物の r_t の固定点 (solve_rt_throat は本物、評価は合成)")
h = Harness([2.10], fake_solve=False)
logd, text, chk, err, dt = run(h, CARD_PHYS, "real_fp")
if err is None:
    raw, full, summ = verify_output("本物の固定点", h, logd, text, chk)
    hist = raw["_card"]["r_t_history"]
    aa = raw["_card"]["auto_adjust"]
    check("本物の固定点: 最初の r_t 0.025 と固定点の後の r_t で選び、一致して出力 (出力の r_t = 最後の選択の r_t)、寸法 PASS。"
          "本物の固定点で自動調整が起きたら、その回数だけ評価を捨てて最初から選び直している",
          len(h.sel) == 2 + len(aa) and h.sel[-2]["r_t"] == 0.025 and h.sel[-1]["r_t"] == raw["spec"]["r_throat"] != 0.025
          and full["redo"] == 0 and full["n_restart"] == len(aa) and hist[-1]["rel_change"] <= 1e-6 and logd["dimensions"]["pass"],
          f"選択の r_t {[s['r_t'] for s in h.sel]} 自動調整 {[(a['key'], a['from'], a['to']) for a in aa]} 寸法 {logd['dimensions']['pass']}")
else:
    check("本物の固定点: build_problem が通る", False, err[:300])
print(f"   (本物の固定点) {dt:.0f} s")

print("== §5.1 #6d 選択の中の長さ合わせの停止 (plan tooling-nozzle-card-local-mesh; plan tooling-nozzle-mknot-fast §4.1 の拡張)")


def stall_of(full, i):
    x = full["selections"][i]
    return x.get("stall")


h = Harness([2.10], unconv=lambda i, M, N, L: N == 600 and abs(M - 2.10) < 1e-9)
logd, text, chk, err, dt = run(h, CARD_PHYS, "s_a")
if err is None:
    raw, full, summ = verify_output("(6d-a) 谷 1 つの中心が未収束", h, logd, text, chk)
    st = stall_of(full, 0)
    aa = raw["_card"]["auto_adjust"]
    check("(6d-a) 谷 1 つで中心が生産の網 600 で未収束 → 選べない例外を捕まえて n_axis_inv 600 → 1200 (on_stall と同じ)、評価を捨てて "
          "粗い掃引 (300) から選び直して 2.10 を選ぶ。記録: 最初の選択の stall (生産の網・600・M_K 2.1・谷 2.10)・調整・捨てた印",
          raw["geometry"]["M_knot"] == 2.1 and raw["geometry"]["n_axis_inv"] == 1200 and full["n_restart"] == 1
          and [(a["key"], a["from"], a["to"]) for a in aa] == [("geometry.n_axis_inv", 600, 1200)]
          and "M_knot の選択" in aa[0]["reason"] and "生産の網 n_axis_inv 600" in aa[0]["reason"]
          and [(sv["N_s"], sv["N_p"], sv["r_t"]) for sv in h.sel] == [(150, 600, 0.025), (300, 1200, 0.025), (300, 1200, rstar(2.10))]
          and h.sel[1]["n_eval"] == h.sel[0]["n_eval"]
          and st is not None and st["grid"] == "prod" and st["n_axis_inv"] == 600 and [r["M_K"] for r in st["rows"]] == [2.1]
          and st["valleys"] == [{"start": 2.1, "final": None, "discard_reason": st["valleys"][0]["discard_reason"], "unconverged": [2.1]}]
          and st["selected_unused"] is None and st["action"] == [{"key": "geometry.n_axis_inv", "from": 600, "to": 1200}]
          and [x["discarded"] for x in full["selections"]] == [True, False, False] and full["selections"][0]["chosen"] is None
          and "geometry.n_axis_inv 600 → 1200" in full["discarded"][0]["reason"]
          and summ["selections"][0]["stall"]["grid"] == "prod" and summ["selections"][0]["stall"]["M_K"] == [2.1],
          f"選択 {[(sv['N_s'], sv['N_p'], sv['r_t']) for sv in h.sel]} 調整 {[(a['key'], a['from'], a['to']) for a in aa]} stall {st}")
    check("(6d-a) YAML のコメントと標準出力の要約の行に停止と調整が出る",
          "長さ合わせの未収束 (生産の網 n_axis_inv 600" in text and "geometry.n_axis_inv 600 → 1200" in text, "")
else:
    check("(6d-a) build_problem が通る", False, err[:400])
print(f"   (6d-a) {dt:.0f} s")

h = Harness([2.10], extra=(3.00, 30.0), unconv=lambda i, M, N, L: N == 600 and abs(M - 2.10) < 1e-9)
logd, text, chk, err, dt = run(h, CARD_PHYS, "s_b")
if err is None:
    raw, full, summ = verify_output("(6d-b) 谷 2 つで片方だけ未収束", h, logd, text, chk)
    st = stall_of(full, 0)
    check("(6d-b) 谷 2 つ (2.10・3.00) で 2.10 の中心だけ生産の網 600 で未収束: 今までは 2.10 を捨てて 3.00 を黙って選んでいた → 選べた値 3.00 を"
          "使わずに n_axis_inv 1200 で粗い掃引から選び直し、2 つの谷を同じ網で比べて 2.10 を選ぶ",
          raw["geometry"]["M_knot"] == 2.1 and raw["geometry"]["n_axis_inv"] == 1200 and full["n_restart"] == 1
          and st is not None and st["grid"] == "prod" and st["selected_unused"] == 3.0 and [r["M_K"] for r in st["rows"]] == [2.1]
          and [v["unconverged"] for v in st["valleys"]] == [[2.1], []] and full["selections"][0]["chosen"] is None
          and full["n_valleys"] == 2 and full["n_valleys_final"] == 2,
          f"出力 {raw['geometry']['M_knot']} n_axis_inv {raw['geometry']['n_axis_inv']} stall {st and {k: st[k] for k in ('selected_unused', 'valleys')}}")
else:
    check("(6d-b) build_problem が通る", False, err[:400])
print(f"   (6d-b) {dt:.0f} s")

h = Harness([2.10], bad=lambda i, M, N: abs(M - 2.10) < 1e-9)
logd, text, chk, err, dt = run(h, CARD_PHYS, "s_c")
check("(6d-c) 制約による不成立だけ (中心 2.10 が生産の網で壁での余裕 0.5°): 調整しない — 今までどおり MknotSelectionError、選び直しも"
      "固定点も無い",
      err is not None and "MknotSelectionError" in err and "生産の網で最終の代表が残る谷が無い" in err and len(h.sel) == 1 and h.solves == [],
      (err or "")[:200])
h = Harness([2.10], extra=(3.00, 30.0), bad=lambda i, M, N: abs(M - 2.10) < 1e-9)
logd, text, chk, err, dt = run(h, CARD_PHYS, "s_c2")
check("(6d-c) 谷 2 つで片方が制約で不成立: 調整しない (長さの停止ではない) — 他方 3.00 を選ぶ (合成の行なので最後の照合で止まるのは試験の作り)、"
      "選び直しは無い",
      len(h.sel) >= 1 and all(sv["N_p"] == 600 for sv in h.sel) and (err is None or "MknotSelectionError" not in err)
      and (logd is None or not logd["auto_adjust"]),
      f"選択 {[(sv['N_p'], sv['c']) for sv in h.sel]} err {(err or '')[:120]}")

h = Harness([2.10], coarse_unconv=lambda i, M, N: N == 150)
logd, text, chk, err, dt = run(h, CARD_PHYS, "s_d")
if err is None:
    raw, full, summ = verify_output("(6d-d) 粗い網が全点で未収束", h, logd, text, chk)
    st = stall_of(full, 0)
    check("(6d-d) 粗い網 150 が全点で未収束 → n_axis_inv 600 → 1200 (粗い網 300) で選び直す",
          raw["geometry"]["n_axis_inv"] == 1200 and full["n_restart"] == 1 and st is not None and st["grid"] == "coarse"
          and st["n_axis_inv"] == 150 and st["L_total_tol"] == 0.1 and len(st["rows"]) == full["selection"]["grid"][2]
          and [(sv["N_s"], sv["N_p"]) for sv in h.sel] == [(150, 600), (300, 1200), (300, 1200)],
          f"選択 {[(sv['N_s'], sv['N_p']) for sv in h.sel]} stall {st and st['grid']}")
else:
    check("(6d-d) build_problem が通る", False, err[:400])
print(f"   (6d-d) {dt:.0f} s")

h = Harness([2.10], unconv=lambda i, M, N, L: L < 0.05 - 1e-12 and abs(M - 2.10) < 1e-9)
logd, text, chk, err, dt = run(h, CARD_PHYS, "s_e")
if err is None:
    raw, full, summ = verify_output("(6d-e) 倍化でも止まり緩和で通る", h, logd, text, chk)
    aa = raw["_card"]["auto_adjust"]
    check("(6d-e) 600 でも 1200 でも未収束、L_total_tol 0.05 で通る: on_stall と同じ順 (n_axis_inv 600 → 1200 → L_total_tol 0.04 → 0.05 と "
          "n_axis_inv を雛形 600 に戻す)、そのたびに粗い掃引から選び直す",
          raw["geometry"]["L_total_tol"] == 0.05 and raw["geometry"]["n_axis_inv"] == 600 and full["n_restart"] == 2
          and [(a["key"], a["from"], a["to"]) for a in aa] == [("geometry.n_axis_inv", 600, 1200), ("geometry.L_total_tol", 0.04, 0.05),
                                                                ("geometry.n_axis_inv", 1200, 600)]
          and raw["_card"]["L_total_tol_relaxed"] is not None
          and [(sv["N_p"], sv["L_p"]) for sv in h.sel] == [(600, 0.04), (1200, 0.04), (600, 0.05), (600, 0.05)]
          and [x.get("stall") is not None for x in full["selections"]] == [True, True, False, False],
          f"選択 {[(sv['N_p'], sv['L_p']) for sv in h.sel]} 調整 {[(a['key'], a['from'], a['to']) for a in aa]}")
else:
    check("(6d-e) build_problem が通る", False, err[:400])
print(f"   (6d-e) {dt:.0f} s")

h = Harness([2.10], unconv=lambda i, M, N, L: abs(M - 2.10) < 1e-9)
logd, text, chk, err, dt = run(h, CARD_PHYS, "s_f")
check("(6d-f) 倍化・緩和の後も未収束: on_stall と同じく 3 回目で例外 (ValueError、YAML を書かない)、選択は 3 回 (600・1200・600 (0.05))",
      err is not None and err.startswith("ValueError") and "内部許容差を 0.05 r_t に緩めても" in err and "M_knot の選択" in err
      and [(sv["N_p"], sv["L_p"]) for sv in h.sel] == [(600, 0.04), (1200, 0.04), (600, 0.05)] and h.solves == [],
      f"選択 {[(sv['N_p'], sv['L_p']) for sv in h.sel]} {(err or '')[:160]}")

print(f"\n{'ALL OK' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
